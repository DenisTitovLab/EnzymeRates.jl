# ABOUTME: Identifiable rank of handwritten/17 without its dead-end steps (P2 probe).
# ABOUTME: Attributes its remaining phantom; stdlib only.
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Oracle.jl"))
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Thermo.jl"))
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

rng = Xoshiro(17)
m = Pop.load_population("handwritten")[17]
for keep in ([1, 2, 3, 4, 5], [1, 3, 4, 5, 6])
used = sort!(unique(vcat([m.steps[j].from for j in keep], [m.steps[j].to for j in keep])))
re = Dict(f => i for (i, f) in enumerate(used))
st = [Pop.PStep(re[m.steps[j].from], re[m.steps[j].to], m.steps[j].consumed, m.steps[j].released, m.steps[j].re) for j in keep]
core = Pop.PMech(; id = "hw17-core", source = "", substrates = m.substrates, products = m.products,
                 metabolites = m.metabolites, species = m.species[used], steps = st, groups = [[j] for j in eachindex(st)])

sol = T.solve(core)
println("steps ", keep, " count ", length(sol.indep), " rank ", ident_rank(core, sol, rng))
end
