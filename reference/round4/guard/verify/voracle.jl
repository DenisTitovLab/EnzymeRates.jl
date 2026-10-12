# ABOUTME: Independent mass-action steady-state oracle (BigFloat linear solve) over plain step
# ABOUTME: lists, with thermodynamically consistent constants; uses neither slice nor package code.
using LinearAlgebra, Random, Printf

setprecision(BigFloat, 2048)

"""One step `from + cons ⇌ to + rel`; `kf/kr` are the unit-speed constants (an RE step's
ratio; its speed is applied by `rate_at`)."""
struct MStep
    from::String
    to::String
    cons::Vector{Symbol}
    rel::Vector{Symbol}
    re::Bool
    kf::BigFloat
    kr::BigFloat
end

forms_of(steps) = unique!(vcat([s.from for s in steps], [s.to for s in steps]))

"""Terminal strongly connected components of the transition graph (arcs with a positive
rate); more than one means the steady state is not unique."""
function terminal_classes(n, arcs)
    adj = [Int[] for _ in 1:n]; radj = [Int[] for _ in 1:n]
    for (a, b) in arcs; push!(adj[a], b); push!(radj[b], a); end
    seen = falses(n); order = Int[]
    function dfs1(u)
        seen[u] = true
        for v in adj[u]; seen[v] || dfs1(v); end
        push!(order, u)
    end
    for u in 1:n; seen[u] || dfs1(u); end
    comp = zeros(Int, n); c = 0
    function dfs2(u)
        comp[u] = c
        for v in radj[u]; comp[v] == 0 && dfs2(v); end
    end
    for u in reverse(order); comp[u] == 0 && (c += 1; dfs2(u)); end
    out = falses(c)
    for (a, b) in arcs; comp[a] != comp[b] && (out[comp[a]] = true); end
    [[i for i in 1:n if comp[i] == k] for k in 1:c if !out[k]]
end

"""Steady-state net consumption rate of `ref` per enzyme, with RE steps at speed
`Λ·speed[k]` (ratio kept). Returns `(v, populations, forms, nterminal)`; `v` is `NaN` when
the steady state is not unique."""
function rate_at(steps, conc::Dict{Symbol, BigFloat}, Λ, speed, ref::Symbol)
    forms = forms_of(steps)
    idx = Dict(f => i for (i, f) in enumerate(forms))
    n = length(forms)
    M = zeros(BigFloat, n, n)
    fl = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
    arcs = Tuple{Int, Int}[]
    for (k, s) in enumerate(steps)
        a, b = idx[s.from], idx[s.to]
        sc = s.re ? big(Λ) * speed[k] : big(1)
        rf = sc * s.kf * prod((conc[x] for x in s.cons); init = big(1))
        rr = sc * s.kr * prod((conc[x] for x in s.rel); init = big(1))
        M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
        rf > 0 && push!(arcs, (a, b)); rr > 0 && push!(arcs, (b, a))
        push!(fl, (a, b, rf, rr, count(==(ref), s.cons) - count(==(ref), s.rel)))
    end
    term = terminal_classes(n, arcs)
    length(term) > 1 && return (big(NaN), nothing, forms, length(term))
    A = copy(M); A[n, :] .= 1
    rhs = zeros(BigFloat, n); rhs[n] = 1
    x = A \ rhs
    v = sum((w * (rf * x[a] - rr * x[b]) for (a, b, rf, rr, w) in fl if w != 0); init = big(0))
    (v, x, forms, 1)
end

lu10(rng, lo, hi) = big(10)^(lo + (hi - lo) * rand(rng))

"""Constants from random free energies of forms and standard potentials of metabolites, so
every cycle obeys detailed balance (Haldane and Wegscheider). `spec` entries are
`(from, to, cons, rel, re)`. SS steps get a random speed in 10^[-1,1]; RE steps speed 1."""
function consistent(spec, rng; mets = nothing)
    forms = unique!(vcat([s[1] for s in spec], [s[2] for s in spec]))
    names = mets === nothing ? unique!(vcat([vcat(s[3], s[4]) for s in spec]...)) : mets
    G = Dict(f => big(randn(rng)) for f in forms)
    mu = Dict(x => big(randn(rng)) for x in names)
    map(spec) do (f, t, c, r, re)
        lk = G[f] + sum((mu[x] for x in c); init = big(0)) - G[t] -
             sum((mu[x] for x in r); init = big(0))
        a = re ? big(1) : lu10(rng, -1, 1)
        MStep(f, t, collect(Symbol, c), collect(Symbol, r), re, a * exp(lk / 2), a * exp(-lk / 2))
    end, mu
end

"""Overall equilibrium constant implied by `mu` for substrates `subs` → products `prods`."""
keq_of(mu, subs, prods) = exp(sum(mu[x] for x in subs) - sum(mu[x] for x in prods))

fmt(x) = isnan(x) ? "NaN" : @sprintf("%.6g", Float64(x))

"""The law's value on the face where `zs` vanish, as the limit from inside (Λ → ∞ first at
`ε·u`, then ε → 0), along `ndir` random directions `u`; returns the values."""
function interior_limits(steps, base, zs, ref, rng; ndir = 3, Λ = big(10)^80,
                         eps = (big(10)^-20, big(10)^-30))
    sp = ones(BigFloat, length(steps))
    out = Vector{Vector{BigFloat}}()
    for _ in 1:ndir
        u = Dict(z => lu10(rng, -2, 2) for z in zs)
        row = BigFloat[]
        for ε in eps
            c = copy(base); foreach(z -> c[z] = ε * u[z], zs)
            push!(row, rate_at(steps, c, Λ, sp, ref)[1])
        end
        push!(out, row)
    end
    out
end

"""The RE limit on the face itself (`zs` exactly zero) for `npat` random speed patterns
(each RE step's speed in 10^[lo,hi]) at Λ = 1e30 and 1e40."""
function face_limits(steps, base, zs, ref, rng; npat = 20, lo = -3, hi = 3)
    c = copy(base); foreach(z -> c[z] = big(0), zs)
    vals = BigFloat[]; conv = BigFloat[]; nterm = 0
    for _ in 1:npat
        sp = [lu10(rng, lo, hi) for _ in steps]
        v30, _, _, nt = rate_at(steps, c, big(10)^30, sp, ref)
        v40, _, _, _ = rate_at(steps, c, big(10)^40, sp, ref)
        nterm = max(nterm, nt)
        push!(vals, v40); push!(conv, abs(v40 - v30))
    end
    (vals = vals, conv = maximum(conv), nterm = nterm)
end
