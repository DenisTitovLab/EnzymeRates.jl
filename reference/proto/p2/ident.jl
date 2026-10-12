# ABOUTME: P2 identifiable-rank check: BigFloat finite-difference rank of dv/dlog(independents) through the
# ABOUTME: oracle, for flux-free SS mechanisms, hand-built dead-end exchange cases and a calibration sample.
include(joinpath(@__DIR__, "..", "common", "Oracle.jl"))
include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))
using Random, Printf
const T = Thermo
setprecision(BigFloat, 512)

"""Rank of a BigFloat matrix by full-pivot elimination with relative tolerance `tol`; also the
largest dropped and smallest kept pivot (relative)."""
function fp_rank(J; tol = big(10)^-25)
    A = copy(J)          # columns are d v / d log θ, already on one scale; no per-column scaling
    piv = BigFloat[]
    nr, nc = size(A); rows = collect(1:nr); cols = collect(1:nc)
    for k in 1:min(nr, nc)
        best = big(0); bi = 0; bj = 0
        for i in k:nr, j in k:nc
            a = abs(A[i, j]); a > best && (best = a; bi = i; bj = j)
        end
        push!(piv, best)
        best == 0 && break
        A[[k, bi], :] = A[[bi, k], :]; A[:, [k, bj]] = A[:, [bj, k]]
        for i in k+1:nr
            f = A[i, k] / A[k, k]
            f == 0 && continue
            for j in k:nc; A[i, j] -= f * A[k, j]; end
        end
    end
    p1 = piv[1]
    r = count(p -> p > tol * p1, piv)
    kept = r == 0 ? big(0) : piv[r] / p1
    dropped = r < length(piv) ? piv[r + 1] / p1 : big(0)
    r, Float64(kept), Float64(dropped)
end

"""Identifiable rank of `sol`'s parameterization of `m`: FD of v over its independents."""
function ident_rank(m, sol, rng; npts = 0, h = big(10)^-30)
    p = length(sol.indep)
    npts = max(npts, 3p + 6)
    keq = Oracle.logu(rng)
    θ = Dict(i => Oracle.logu(rng) for i in sol.indep)
    concs = [Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites) for _ in 1:npts]
    J = zeros(BigFloat, npts, p)
    for (col, i) in enumerate(sol.indep), sg in (1, -1)
        θp = copy(θ); θp[i] *= exp(sg * h)
        k = T.stored_constants(sol, T.all_values(sol, θp, keq))
        for (r, c) in enumerate(concs)
            J[r, col] += sg * Oracle.steady_state_v(m, k, c) / (2h)
        end
    end
    fp_rank(J)
end

sp(conf, bound...) = Pop.PSpecies(collect(bound), conf, Symbol[], Symbol[],
                                  Symbol(conf, join(string.(getfield.(collect(bound), :name)))))
S = Pop.PMet(:S, :S); P = Pop.PMet(:P, :P); I = Pop.PMet(:I, :Inh); J = Pop.PMet(:J, :Inh)
species = [sp(:E), sp(:E, S), sp(:E, P), sp(:E, I), sp(:E, J)]
core = [Pop.PStep(1, 2, [S], Pop.PMet[], true),       # E + S ⇌ ES
        Pop.PStep(2, 3, Pop.PMet[], Pop.PMet[], false), # ES <--> EP
        Pop.PStep(1, 3, [P], Pop.PMet[], true),       # E + P ⇌ EP
        Pop.PStep(1, 4, [I], Pop.PMet[], true)]       # E + I ⇌ EI
function hand(id, steps)
    used = sort!(unique(vcat([s.from for s in steps], [s.to for s in steps])))
    re = Dict(f => i for (i, f) in enumerate(used))
    st = [Pop.PStep(re[s.from], re[s.to], s.consumed, s.released, s.re) for s in steps]
    mets = unique([:S, :P, [x.name for s in steps for x in [s.consumed; s.released]]...])
    Pop.PMech(; id = id, source = id, substrates = [:S], products = [:P],
              inhibitors = setdiff(mets, [:S, :P]), metabolites = mets, species = species[used],
              steps = st, groups = [[j] for j in eachindex(st)])
end
handbuilt = [
    # ROADMAP: E(I) + J <--> E(J) + I where E(J) binds nothing more (a bridge)
    hand("hand:exchange-bridge", [core; Pop.PStep(4, 5, [J], [I], false)]),
    # the same exchange, but E(J) also releases J to E: a Wegscheider cycle, no bridge
    hand("hand:exchange-in-wegscheider-cycle", [core; Pop.PStep(4, 5, [J], [I], false);
                                                Pop.PStep(1, 5, [J], Pop.PMet[], true)]),
    # plain dead-end SS binding (a bridge)
    hand("hand:dead-end-ss-binding", [core[1:3]; Pop.PStep(1, 4, [I], Pop.PMet[], false)]),
]

rng = Xoshiro(11)
results = []
function report(m; tag = "")
    sol = T.solve(m)
    z = T.zero_flux_ss(m)
    m2, nconv = T.with_flux_free_as_re(m)
    r1, k1, d1 = ident_rank(m, sol, rng)
    line = @sprintf("%-46s count %2d  rank %2d  (kept %.1e, dropped %.1e)", m.id * tag, length(sol.indep), r1, k1, d1)
    if nconv > 0
        sol2 = T.solve(m2)
        r2, k2, d2 = ident_rank(m2, sol2, rng)
        # same rate at matching constants: K = kr/kf (binding) or kf/kr (other) for converted steps
        keq = Oracle.logu(rng); θ = Dict(i => Oracle.logu(rng) for i in sol.indep)
        k = T.stored_constants(sol, T.all_values(sol, θ, keq))
        k2c = [m2.steps[j].re && !m.steps[j].re ?
               [Pop.is_binding(m.steps[j]) ? k[j][2] / k[j][1] : k[j][1] / k[j][2]] : k[j]
               for j in eachindex(k)]
        c = Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites)
        va = Oracle.steady_state_v(m, k, c); vb = Oracle.steady_state_v(m2, k2c, c)
        line *= @sprintf("  | flux-free SS groups made RE: %d -> count %2d  rank %2d  | v(SS) vs v(RE) rel %.1e",
                         nconv, length(sol2.indep), r2, Float64(abs(va - vb) / abs(va)))
    end
    println(line)
    push!(results, (m.id, length(sol.indep), r1))
end

println("— hand-built dead-end cases —")
foreach(report, handbuilt)
println("— whole flux-free SS groups in the populations —")
pops = Dict(p => Pop.load_population(p) for p in Pop.populations())
for (p, i) in (("handwritten", 17), ("specs", 14))
    report(pops[p][i])
end
println("— calibration: every plain spec and handwritten mechanism, 40 enumerated —")
for p in ("specs", "handwritten"), m in pops[p]
    m.allosteric && continue
    m.id in ("handwritten/17", "specs/14") && continue
    report(m)
end
enum = [m for p in ("init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi") for m in pops[p] if !m.allosteric]
for m in enum[sort(randperm(Xoshiro(3), length(enum))[1:40])]
    report(m)
end
gaps = [(id, c - r) for (id, c, r) in results if c != r]
println("mechanisms with count > rank: ", length(gaps), " of ", length(results))
foreach(g -> println("   ", g), gaps)
