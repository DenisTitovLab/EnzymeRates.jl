# ABOUTME: Shared helpers of the P1-engine skeptic: loads Pop/Oracle/Poly/Engine/EngineLean read-only and adds
# ABOUTME: independent checks (oracle error, speed-free tests, exact speed-free predicates, dfree vs populations).
const PROTO = normpath(joinpath(@__DIR__, ".."))
isdefined(Main, :Pop) || include(joinpath(PROTO, "common", "Population.jl"))
isdefined(Main, :Oracle) || include(joinpath(PROTO, "common", "Oracle.jl"))
isdefined(Main, :Poly) || include(joinpath(PROTO, "common", "Poly.jl"))
isdefined(Main, :Engine) || include(joinpath(PROTO, "common", "Engine.jl"))
isdefined(Main, :EngineLean) || include(joinpath(PROTO, "p1_engine", "EngineLean.jl"))

module Sk

using Random, LinearAlgebra, Printf
import ..Pop, ..Oracle, ..Poly, ..Engine, ..EngineLean

const PR = Poly.P{Rational{Int}}

"""Values vector of the bare law `L` (Engine.Law) at per-step constants `k` and concentrations `c`."""
function vals_of(L, m, k, c)
    v = ones(BigFloat, length(L.vars.names))
    for (j, vs) in enumerate(L.stepvars), (i, x) in enumerate(vs)
        v[x] = k[j][i]
    end
    for x in m.metabolites
        v[L.vars.index[x]] = c[x]
    end
    v
end
lawv(L, vals) = Poly.evaluate(L.num, vals) / Poly.evaluate(L.den, vals)
relerr(a, b) = Float64(abs(a - b) / abs(b))

"""Random bare constants: RE steps consistent (from `consistent_constants`), SS steps free."""
function bare_constants(m, rng; Keq = Oracle.logu(rng))
    kc, _, _ = Oracle.consistent_constants(m, rng; Keq)
    [s.re ? kc[j] : [Oracle.logu(rng), Oracle.logu(rng)] for (j, s) in enumerate(m.steps)]
end
random_concs(m, rng) = Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites)

"""SS steps whose kf and kr scaled together (×3) leave the law's value unchanged at `k`, `c`."""
function speed_free_numeric(L, m, k, c; tol = big(10)^-60)
    v0 = lawv(L, vals_of(L, m, k, c))
    [j for (j, s) in enumerate(m.steps) if !s.re &&
     abs(lawv(L, vals_of(L, m, [i == j ? 3 .* q : q for (i, q) in enumerate(k)], c)) - v0) <=
     tol * abs(v0)]
end
"""The same test against the BigFloat mass-action oracle (independent of the engine)."""
function speed_free_oracle(m, k, c; tol = big(10)^-30)
    v0 = Oracle.steady_state_v(m, k, c)
    [j for (j, s) in enumerate(m.steps) if !s.re &&
     abs(Oracle.steady_state_v(m, [i == j ? 3 .* q : q for (i, q) in enumerate(k)], c) - v0) <=
     tol * abs(v0)]
end

"""
Segment multigraph of `m`: `seg[i]` (segment of form i, by union-find over RE steps), and
for every SS step `(j, u, v, ω)` with ω the copies of the first substrate one forward
turnover of the step takes up, RE paths from each segment's lowest-index form included.
"""
function segment_graph(m)
    n = length(m.species)
    p = collect(1:n)
    root(i) = p[i] == i ? i : (p[i] = root(p[i]))
    for s in m.steps
        s.re && (p[root(s.from)] = root(s.to))
    end
    roots = unique(root(i) for i in 1:n)
    seg = [findfirst(==(root(i)), roots) for i in 1:n]
    x = m.substrates[1]
    dx(s) = count(y -> y.name === x, s.consumed) - count(y -> y.name === x, s.released)
    u = fill(typemin(Int), n)                  # x exponent of each form relative to its root
    for r in roots
        u[r] = 0; queue = [r]
        while !isempty(queue)
            a = popfirst!(queue)
            for s in m.steps
                s.re || continue
                if s.from == a && u[s.to] == typemin(Int)
                    u[s.to] = u[a] + dx(s); push!(queue, s.to)
                elseif s.to == a && u[s.from] == typemin(Int)
                    u[s.from] = u[a] - dx(s); push!(queue, s.from)
                end
            end
        end
    end
    edges = [(j, seg[s.from], seg[s.to], dx(s) + u[s.from] - u[s.to])
             for (j, s) in enumerate(m.steps) if !s.re]
    seg, length(roots), edges
end

"""Biconnected blocks of a multigraph (parallel edges distinct; self-loops excluded) as lists
of positions into `es` = [(u, v)]."""
function multiblocks(G, es)
    adj = [Tuple{Int, Int}[] for _ in 1:G]
    for (k, (a, b)) in enumerate(es)
        a == b && continue
        push!(adj[a], (b, k)); push!(adj[b], (a, k))
    end
    disc, low, stack, out, t = zeros(Int, G), zeros(Int, G), Int[], Vector{Int}[], Ref(0)
    function dfs(a, pk)
        t[] += 1; disc[a] = low[a] = t[]
        for (b, k) in adj[a]
            k == pk && continue
            if disc[b] == 0
                push!(stack, k); dfs(b, k)
                low[a] = min(low[a], low[b])
                if low[b] >= disc[a]
                    blk = Int[]
                    while true
                        f = pop!(stack); push!(blk, f); f == k && break
                    end
                    push!(out, blk)
                end
            elseif disc[b] < disc[a]
                push!(stack, k); low[a] = min(low[a], disc[b])
            end
        end
    end
    for a in 1:G
        disc[a] == 0 && dfs(a, 0)
    end
    out
end

"""
Exact speed-free predicates from flux decomposition. An SS step's speed drops out of v exactly
when its net flux is zero for all constants of the regime: `:thermo` (every non-reactive cycle
at detailed balance): a silent self-loop (ω = 0), or a step whose biconnected block of the
segment multigraph holds no cycle with nonzero turnover (ω not a gradient on the block);
`:bare` (independent constants): a silent self-loop or a bridge.
"""
function speed_free_exact(m)
    seg, G, edges = segment_graph(m)
    loops = [e for e in edges if e[2] == e[3]]
    links = [e for e in edges if e[2] != e[3]]
    blocks = multiblocks(G, [(e[2], e[3]) for e in links])
    thermo = Set(e[1] for e in loops if e[4] == 0)
    bare = copy(thermo)
    for blk in blocks
        # ω a gradient on the block ⇔ no reactive cycle in it.
        φ = Dict{Int, Int}()
        bes = [links[k] for k in blk]
        φ[bes[1][2]] = 0
        changed = true
        while changed
            changed = false
            for (_, a, b, w) in bes
                if haskey(φ, a) && !haskey(φ, b)
                    φ[b] = φ[a] + w; changed = true
                elseif haskey(φ, b) && !haskey(φ, a)
                    φ[a] = φ[b] - w; changed = true
                end
            end
        end
        gradient = all(φ[b] - φ[a] == w for (_, a, b, w) in bes)
        gradient && union!(thermo, (e[1] for e in bes))
        length(bes) == 1 && push!(bare, bes[1][1])
    end
    sort!(collect(thermo)), sort!(collect(bare))
end

"""Free form index (`findfirst(Pop.is_free)`) and the lowest-index root of its segment as the
package and the engine choose it."""
function free_and_root(m)
    i = findfirst(Pop.is_free, m.species)
    i === nothing && return nothing, nothing
    seg, _, _ = segment_graph(m)
    S = [f for f in eachindex(m.species) if seg[f] == seg[i]]
    hasres(sp) = !isempty(sp.res_added) || !isempty(sp.res_subtracted)
    r = argmin(f -> (length(m.species[f].bound), hasres(m.species[f]) ? 1 : 0,
                     string(m.species[f].name)), S)
    i, r
end

"""Relative error of dfree/den (engine) against the oracle's steady-state population of the
free form, at bare constants `k` and concentrations `c`."""
function dfree_vs_population(L, m, k, c)
    i, _ = free_and_root(m)
    i === nothing && return NaN
    vals = vals_of(L, m, k, c)
    share = Poly.evaluate(L.dfree, vals) / Poly.evaluate(L.den, vals)
    x = Oracle.steady_state_populations(m, k, c)
    relerr(share, x[i])
end

"""Exact symbolic form `Dict{Vector{Pair{Symbol,Int}},Rational{Int}}` of a Poly polynomial."""
sym(p, V) = Poly.symbolic(p, V)
topoly(pp) = Dict{Vector{Pair{Symbol, Int}}, Rational{Int}}(mono => c for (mono, c) in pp)

"""
Terms of a package polynomial (`Pop.PPoly`) with the dependents of `m` substituted: per term
its coefficient, its exponent of each parameter in `θ` (Rational), of `Keq`, and of each
metabolite. Errors when a symbol is neither a parameter of `θ`, a dependent, `Keq` nor a
metabolite.
"""
function reduced_terms(pp, m, θ::Vector{Symbol})
    dep = Dict(k => Dict(mono) for (k, mono) in m.dependent_monomials)
    pidx = Dict(s => i for (i, s) in enumerate(θ))
    midx = Dict(s => i for (i, s) in enumerate(m.metabolites))
    out = Tuple{Rational{Int}, Vector{Rational{Int}}, Rational{Int}, Vector{Int}}[]
    for (mono, c) in pp
        ep = zeros(Rational{Int}, length(θ)); ek = 0 // 1; ec = zeros(Int, length(m.metabolites))
        add!(s, e) = s === :Keq ? (ek += e) : haskey(pidx, s) ? (ep[pidx[s]] += e) :
                     haskey(midx, s) ? (ec[midx[s]] += Int(e)) :
                     error("symbol $s is not fitted, dependent, Keq or a metabolite")
        for (s, e) in mono
            if haskey(dep, s)
                for (t, f) in dep[s]; add!(t, e * f); end
            else
                add!(s, e)
            end
        end
        push!(out, (c, ep, ek, ec))
    end
    out
end

"""Full-pivoting Gaussian elimination rank of a BigFloat matrix (columns scaled to unit norm; a
column whose norm is below `tol` times the largest is a zero column, so normalizing cannot
amplify its rounding noise), with pivots below `tol` times the largest counted as zero."""
function numeric_rank(A::Matrix{BigFloat}; tol = big(10)^-40)
    A = copy(A)
    norms = [sqrt(sum(abs2, A[:, j])) for j in axes(A, 2)]
    nmax = maximum(norms; init = big(0))
    for j in axes(A, 2)
        norms[j] <= tol * nmax ? (A[:, j] .= 0) : (A[:, j] ./= norms[j])
    end
    r = 0; m, n = size(A); rows = collect(1:m); cols = collect(1:n)
    big0 = maximum(abs, A; init = big(0))
    big0 == 0 && return 0
    while r < min(m, n)
        best, bi, bj = big(0), 0, 0
        for i in (r + 1):m, j in (r + 1):n
            a = abs(A[i, j]); a > best && ((best, bi, bj) = (a, i, j))
        end
        best <= tol * big0 && break
        r += 1
        A[[r, bi], :] = A[[bi, r], :]; A[:, [r, bj]] = A[:, [bj, r]]
        for i in (r + 1):m
            f = A[i, r] / A[r, r]
            f == 0 && continue
            @views A[i, r:end] .-= f .* A[r, r:end]
        end
    end
    r
end

"""
Local identifiable dimension of the package's Reduced law of `m` in its fitted parameters
(Keq and E_total known): the rank of ∂v/∂log θ stacked over random concentration points, in
256-bit BigFloat at a random θ. Returns `(n_fitted, rank)`.
"""
function identifiable_rank(m, rng; npoints = 0)
    θ = m.fitted_params
    nθ = length(θ)
    nθ == 0 && return 0, 0
    setprecision(BigFloat, 256) do
        N = reduced_terms(m.pkg_num, m, θ); D = reduced_terms(m.pkg_den, m, θ)
        θv = [exp(big(log(10.0)) * (4 * big(rand(rng)) - 2)) for _ in 1:nθ]
        Keq = exp(big(log(10.0)) * (4 * big(rand(rng)) - 2))
        tv(t) = BigFloat(numerator(t[1])) / denominator(t[1]) *
                prod((θv[i]^t[2][i] for i in 1:nθ if t[2][i] != 0); init = big(1)) * Keq^t[3]
        Nθ = map(tv, N); Dθ = map(tv, D)
        np = npoints == 0 ? 3 * nθ + 6 : npoints
        J = zeros(BigFloat, np, nθ)
        for p in 1:np
            c = [exp(big(log(10.0)) * (6 * big(rand(rng)) - 3)) for _ in m.metabolites]
            cm(t) = prod((c[i]^t[4][i] for i in eachindex(c) if t[4][i] != 0); init = big(1))
            nv = [Nθ[k] * cm(N[k]) for k in eachindex(N)]
            dv = [Dθ[k] * cm(D[k]) for k in eachindex(D)]
            Nt, Dt = sum(nv), sum(dv)
            for j in 1:nθ
                Nj = sum(nv[k] * N[k][2][j] for k in eachindex(N))
                Dj = sum(dv[k] * D[k][2][j] for k in eachindex(D))
                J[p, j] = (Nj * Dt - Nt * Dj) / Dt^2
            end
        end
        nθ, numeric_rank(J)
    end
end

end # module
