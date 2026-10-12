# ABOUTME: Identifiable rank of every plain mechanism with a flux-free SS step (P2 probe).
# ABOUTME: Shows grouped flux-free members are no phantom; stdlib only.
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

rng = Xoshiro(13)
nm = 0; eq = 0; gaps = String[]
for p in Pop.populations(), m in Pop.load_population(p)
    m.allosteric && continue
    any(T.zero_flux_ss(m)) || continue
    sol = T.solve(m)
    m2, nconv = T.with_flux_free_as_re(m)
    s2 = nconv > 0 ? T.solve(m2) : sol
    r, k, d = ident_rank(m, sol, rng)
    global nm += 1
    length(s2.indep) == r ? (global eq += 1) :
        push!(gaps, @sprintf("%s count %d (after flux-free fix %d) rank %d kept %.1e dropped %.1e", m.id, length(sol.indep), length(s2.indep), r, k, d))
end
println("mechanisms with flux-free SS steps: $nm; count after the whole-group fix == rank: $eq")
foreach(println, gaps)
