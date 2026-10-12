# ABOUTME: Lean bare rate law of a plain Pop.PMech (module EngineLean): RE segments, whole-graph spanning
# ABOUTME: trees oriented per root, edge-sum numerator, and exact division of a numerator by (1 − Γ/Keq).
isdefined(Main, :Pop) || include(joinpath(@__DIR__, "..", "common", "Population.jl"))
isdefined(Main, :Poly) || include(joinpath(@__DIR__, "..", "common", "Poly.jl"))

module EngineLean

import ..Pop
using ..Poly: Vars, var!, Mono, mono, mmul, mpow, mexp, P, scale, expand, lowest_shift, term,
              addto!

const PI = P{Int}

uf_root(p, i) = p[i] == i ? i : (p[i] = uf_root(p, p[i]))

"""Undirected spanning trees of the graph on nodes `1:k` with edges `pairs`, as lists of edge
positions (recursive include/exclude with union-find)."""
function spanning_trees(k, pairs)
    trees = Vector{Int}[]
    function rec(i, chosen, p)
        length(chosen) == k - 1 && return push!(trees, copy(chosen))
        length(pairs) - i + 1 < k - 1 - length(chosen) && return
        a, b = uf_root(p, pairs[i][1]), uf_root(p, pairs[i][2])
        if a != b
            q = copy(p); q[a] = b
            push!(chosen, i); rec(i + 1, chosen, q); pop!(chosen)
        end
        rec(i + 1, chosen, p)
    end
    rec(1, Int[], collect(1:k))
    trees
end

"""
    bare_rate_law(m) -> (vars, stepvars, num, den, dfree)

`v = num/den` per unit enzyme (the net consumption of `m.substrates[1]`) with one variable per
step constant (`K#j`, or `kf#j`/`kr#j` in the stored direction) and one per metabolite. The RE
weights follow a breadth-first tree that scans the steps in `order`.
"""
function bare_rate_law(m::Pop.PMech; order = eachindex(m.steps))
    V = Vars()
    conc = Dict(x => var!(V, x) for x in m.metabolites)
    stepvars = [s.re ? [var!(V, Symbol("K#", j))] :
                [var!(V, Symbol("kf#", j)), var!(V, Symbol("kr#", j))]
                for (j, s) in enumerate(m.steps)]
    cmono(ms) = foldl(mmul, (mono(conc[x.name] => 1) for x in ms); init = Mono())
    ratio = [!s.re ? Mono() : Pop.is_binding(s) ?
             mmul(cmono(s.consumed), mono(stepvars[j][1] => -1)) :
             mmul(mono(stepvars[j][1] => 1), mmul(cmono(s.consumed), mpow(cmono(s.released), -1)))
             for (j, s) in enumerate(m.steps)]

    # RE segments, rooted as the package roots them, with each form's weight from its root.
    n = length(m.species)
    p = collect(1:n)
    for s in m.steps
        s.re && (p[uf_root(p, s.from)] = uf_root(p, s.to))
    end
    roots = unique(uf_root(p, i) for i in 1:n)
    segments = [[i for i in 1:n if uf_root(p, i) == r] for r in roots]
    seg = [findfirst(==(uf_root(p, i)), roots) for i in 1:n]
    hasres(sp) = !isempty(sp.res_added) || !isempty(sp.res_subtracted)
    alpha, seen = [Mono() for _ in 1:n], falses(n)
    for S in segments
        r = argmin(i -> (length(m.species[i].bound), hasres(m.species[i]) ? 1 : 0,
                         string(m.species[i].name)), S)
        seen[r] = true; queue = [r]
        while !isempty(queue)
            cur = popfirst!(queue)
            for j in order
                s = m.steps[j]
                s.re || continue
                a, b, w = s.from == cur ? (s.from, s.to, ratio[j]) :
                          s.to == cur ? (s.to, s.from, mpow(ratio[j], -1)) : (0, 0, Mono())
                (a == 0 || seen[b]) && continue
                alpha[b] = mmul(alpha[a], w); seen[b] = true; push!(queue, b)
            end
        end
    end
    isconc(q) = q.first in values(conc)
    for (j, s) in enumerate(m.steps)
        s.re || continue
        filter(isconc, mmul(alpha[s.to], mpow(alpha[s.from], -1))) == filter(isconc, ratio[j]) ||
            error("rate_equation: the rapid-equilibrium steps close a catalytic cycle " *
                  "(through $(m.species[s.from].name) ⇌ $(m.species[s.to].name)), so the " *
                  "mechanism has no finite rate. Make one step of the cycle steady-state.")
    end

    # SS steps between segments: merged weights per segment pair; trees oriented per root.
    G = length(segments)
    pairs, W = Tuple{Int, Int}[], Dict{Tuple{Int, Int}, PI}()   # W[(a, b)]: out of a toward b
    x = conc[m.substrates[1]]
    copies(ms) = count(y -> conc[y.name] == x, ms)
    num_terms = Tuple{Int, Int, Int, Mono, Mono}[]              # (ω, from seg, to seg, fwd, rev)
    for (j, s) in enumerate(m.steps)
        s.re && continue
        a, b = seg[s.from], seg[s.to]
        fwd = mmul(mmul(mono(stepvars[j][1] => 1), cmono(s.consumed)), alpha[s.from])
        rev = mmul(mmul(mono(stepvars[j][2] => 1), cmono(s.released)), alpha[s.to])
        ω = copies(s.consumed) - copies(s.released) + mexp(alpha[s.from], x) - mexp(alpha[s.to], x)
        ω == 0 || push!(num_terms, (ω, a, b, fwd, rev))
        a == b && continue
        minmax(a, b) in pairs || push!(pairs, minmax(a, b))
        W[(a, b)] = get(W, (a, b), PI()) + term(PI, fwd)
        W[(b, a)] = get(W, (b, a), PI()) + term(PI, rev)
    end
    trees = spanning_trees(G, pairs)
    G == 1 || !isempty(trees) || error("the steady-state steps leave the segment graph disconnected")
    D = [sum(expand(PI, oriented(pairs[t], r, W)) for t in trees; init = PI()) for r in 1:G]
    den = sum(sum(term(PI, alpha[i]) for i in segments[g]) * D[g] for g in 1:G)
    num = sum((scale(D[a] * fwd - D[b] * rev, ω) for (ω, a, b, fwd, rev) in num_terms);
              init = PI())
    i_free = findfirst(Pop.is_free, m.species)
    dfree = i_free === nothing ? one(PI) : D[seg[i_free]]
    shift = lowest_shift((num, den), values(conc))
    V, stepvars, num * shift, den * shift, dfree * shift
end

"""Weights of the tree with edges `tp` (segment pairs) oriented toward `r`."""
function oriented(tp, r, W)
    ws, seen, frontier = PI[], Set([r]), [r]
    while !isempty(frontier)
        q = popfirst!(frontier)
        for (a, b) in tp
            c = a == q && b ∉ seen ? b : b == q && a ∉ seen ? a : 0
            c == 0 && continue
            push!(ws, W[(c, q)]); push!(seen, c); push!(frontier, c)
        end
    end
    ws
end

"""Exact quotient of `p` by `1 − x` for a monomial `x` that holds variable `kv` (Keq) to the
power −1, by chain division (terms equal up to powers of `x` form one Laurent polynomial in
`x`); `nothing` when `1 − x` does not divide `p`."""
function divide_one_minus(p::P{T}, x::Mono, kv::Int) where {T}
    chains = Dict{Mono, Vector{Pair{Int, T}}}()
    for (t, c) in p.t
        e = mexp(t, kv)
        push!(get!(chains, mmul(t, mpow(x, e)), Pair{Int, T}[]), -e => c)
    end
    q = P{T}()
    for (key, terms) in chains
        sort!(terms; by = first)
        acc = zero(T)
        for (k, (i, c)) in enumerate(terms)
            acc += c
            k == length(terms) && break
            for j in i:(terms[k + 1].first - 1)
                addto!(q.t, mmul(key, mpow(x, j)), acc)
            end
        end
        iszero(acc) || return nothing
    end
    q
end

end # module
