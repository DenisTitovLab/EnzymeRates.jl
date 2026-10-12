# ABOUTME: Overshoot helpers: a sup search at zero products, abortive forms, the non-turnover
# ABOUTME: cycle support, the probe-style overshoot class, and a per-law kcat-vs-sup test.

"""Largest v/E_total at products and modifiers zero over log10 substrate concentrations in
[lo, hi]^d: random sampling, then a shrinking compass search from the best five."""
function sup_search(nl, nc, subs, rng; lo = -8.0, hi = 8.0, nsample = 400)
    d = length(subs)
    c = zeros(nc)
    function f(y)
        for (j, i) in enumerate(subs)
            c[i] = 10.0^y[j]
        end
        v = rate(nl, c)
        isfinite(v) ? v : -Inf
    end
    ys = [lo .+ (hi - lo) .* rand(rng, d) for _ in 1:(nsample * d)]
    fs = f.(ys)
    best_y, best_f = copy(ys[argmax(fs)]), maximum(fs)
    for s in sortperm(fs; rev = true)[1:min(5, length(ys))]
        y, fy, step = copy(ys[s]), fs[s], 2.0
        while step > 1e-3
            moved = false
            for j in 1:d, sg in (1.0, -1.0)
                z = copy(y)
                z[j] = clamp(z[j] + sg * step, lo, hi)
                fz = f(z)
                if fz > fy * (1 + 1e-12)
                    y, fy, moved = z, fz, true
                end
            end
            moved || (step /= 2)
        end
        fy > best_f && ((best_y, best_f) = (y, fy))
    end
    best_f, best_y
end

"""Whether form `i` is abortive: it holds a substrate and a product, or a substrate that its
residual already holds."""
function abortive(m, i)
    f = m.forms[i]
    b = [x.name for x in S.mets_of(m, f.bound)]
    r = [x.name for x in S.mets_of(m, f.res)]
    bs = filter(in(m.substrates), b)
    !isempty(bs) && (any(in(m.products), b) || any(in(r), bs))
end

"""Edges on some cycle with no net turnover (own spanning tree and cycle basis)."""
function nonturnover_support(m)
    F, E = length(m.forms), length(m.edges)
    s1 = m.substrates[1]
    adj = [Int[] for _ in 1:F]
    for (k, e) in enumerate(m.edges)
        push!(adj[e.a], k); push!(adj[e.b], k)
    end
    parent_edge = zeros(Int, F); seen = falses(F); seen[1] = true
    queue = [1]; tree = falses(E)
    while !isempty(queue)
        u = popfirst!(queue)
        for k in adj[u]
            e = m.edges[k]
            w = e.a == u ? e.b : e.a
            seen[w] && continue
            seen[w] = true; parent_edge[w] = k; tree[k] = true; push!(queue, w)
        end
    end
    # path from the root to each form as a signed edge vector (direction a -> b positive)
    function path(u)
        z = zeros(Int, E)
        while u != 1
            k = parent_edge[u]
            e = m.edges[k]
            z[k] += e.b == u ? 1 : -1
            u = e.b == u ? e.a : e.b
        end
        z
    end
    take(k) = (e = m.edges[k];
               count(x -> x.name === s1 && !x.modifier, S.mets_of(m, e.cons)) -
               count(x -> x.name === s1 && !x.modifier, S.mets_of(m, e.rel)))
    cyc = Vector{Int}[]
    for k in 1:E
        tree[k] && continue
        e = m.edges[k]
        z = path(e.a) .- path(e.b)
        z[k] += 1
        push!(cyc, z)
    end
    mult = [sum(z[k] * take(k) for k in 1:E) for z in cyc]
    ref = findfirst(!=(0), mult)
    sup = falses(E)
    for (z, c) in zip(cyc, mult)
        v = c == 0 ? z : ref === nothing || z === cyc[ref] ? nothing : mult[ref] .* z .- c .* cyc[ref]
        v === nothing || (sup .|= v .!= 0)
    end
    sup
end

function overshoot_class(m)
    sup = nonturnover_support(m) .& S.catalytic_core(m)
    ab(k) = abortive(m, m.edges[k].a) || abortive(m, m.edges[k].b)
    (any(k -> sup[k] && !m.edges[k].re, eachindex(m.edges)),
     any(k -> sup[k] && ab(k), eachindex(m.edges)),
     any(k -> sup[k] && !m.edges[k].re && ab(k), eachindex(m.edges)))
end

"""Search `npts` parameter points (log-uniform in [10^plo, 10^phi]); per point (kcat, sup,
argmax) and, for plain laws, the oracle at the argmax of the best point."""
function test_law(m, npts, rng; plo = -2.0, phi = 2.0, oracle_check = true)
    ok, why, sol, L = derive_law(m)
    ok || return (ok = false, why = why)
    cls, vars, terms = kcat_class(m, sol)
    cls === :finite || return (ok = false, why = "class $cls")
    names, subs, prods, mods = roles(m)
    best = (ratio = 0.0, kc = NaN, sup = NaN, at = Float64[], pd = nothing)
    for _ in 1:npts
        pd = draw_params(L.vars, rng; lo = plo, hi = phi)
        kc = kcat_value(vars, terms, pd)
        nl = NL(L, pd)
        sup, at = sup_search(nl, length(names), subs, rng)
        sup / kc > best.ratio && (best = (ratio = sup / kc, kc = kc, sup = sup, at = at, pd = pd))
    end
    orc = NaN
    if oracle_check && !S.is_allosteric(m) && best.pd !== nothing
        c = zeros(BigFloat, length(names))
        foreach(((j, i),) -> c[i] = big(10.0)^best.at[j], enumerate(subs))
        vo, _ = oracle(m, edge_constants(m, sol, best.pd), c)
        orc = Float64(abs(vo - best.sup) / abs(vo))
    end
    (ok = true, ratio = best.ratio, kc = best.kc, sup = best.sup, at = best.at, orc = orc)
end

