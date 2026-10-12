# ABOUTME: Bare rate law of a plain Pop.PMech from its reaction graph (module Engine, probe P1).
# ABOUTME: Cha contraction, spanning trees per biconnected block, Hill cycle-flux numerator; stdlib only.
isdefined(Main, :Pop) || include(joinpath(@__DIR__, "Population.jl"))
isdefined(Main, :Poly) || include(joinpath(@__DIR__, "Poly.jl"))

module Engine

import ..Pop
using ..Poly: Poly, Vars, var!, Mono, mono, mmul, mpow, mexp, P, nterms, scale, expand,
              lowest_shift, subst, term

const PI = P{Int}

"""A steady-state step as an edge of the segment multigraph: it joins segment `u` (of its
`from` form) to segment `v`, with monomial weights `fwd` (u → v: kf, the consumed
concentrations, the from form's RE weight) and `rev` (v → u), and `ω`, the copies of the
first substrate one forward turnover takes up, the RE paths at both ends included."""
struct Edge
    step::Int
    u::Int
    v::Int
    fwd::Mono
    rev::Mono
    ω::Int
end

"""Weight of leaving segment `a` along edge `e`."""
wout(e::Edge, a) = a == e.u ? e.fwd : e.rev

"""The biconnected blocks of a segment multigraph: per block its nodes, its edges (indices
into the edge list), its tree sum toward each of its nodes, and for every segment `g` the
block node nearest to `g`."""
struct Blocks
    nodes::Vector{Vector{Int}}
    edges::Vector{Vector{Int}}
    factor::Vector{Dict{Int, PI}}
    attach::Vector{Vector{Int}}
end

"""Union-find root of `i` with path compression."""
uf_root(p, i) = p[i] == i ? i : (p[i] = uf_root(p, p[i]))

"""Components of the nodes `1:G` joined by the node pairs `pairs`, as a root per node."""
function components(G, pairs)
    p = collect(1:G)
    for (a, b) in pairs
        p[uf_root(p, a)] = uf_root(p, b)
    end
    [uf_root(p, i) for i in 1:G]
end

"""Edge-index lists of the biconnected components (Hopcroft–Tarjan) of the simple graph
with adjacency `adj[u] = [(v, key)]`."""
function bicomponents(G, adj)
    disc, low, stack, out = zeros(Int, G), zeros(Int, G), Int[], Vector{Int}[]
    t = Ref(0)
    function dfs(u, pk)
        t[] += 1; disc[u] = low[u] = t[]
        for (v, k) in adj[u]
            k == pk && continue
            if disc[v] == 0
                push!(stack, k); dfs(v, k)
                low[u] = min(low[u], low[v])
                if low[v] >= disc[u]
                    blk = Int[]
                    while true
                        f = pop!(stack); push!(blk, f); f == k && break
                    end
                    push!(out, blk)
                end
            elseif disc[v] < disc[u]
                push!(stack, k); low[u] = min(low[u], disc[v])
            end
        end
    end
    for u in 1:G
        disc[u] == 0 && dfs(u, 0)
    end
    out
end

"""Undirected spanning trees of the graph on `k` local nodes with edges `pairs`, as lists
of edge positions (recursive include/exclude with union-find)."""
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

"""Blocks of the multigraph of `edges` on segments `1:G`. Parallel edges merge into one
simple edge whose directed weights add (Volkenstein–Goldstein rule a); each block's tree
sum toward a node is the sum, over the block's undirected spanning trees, of the product of
the tree's edges oriented toward that node. Errors when the graph is disconnected."""
function blocks(G, edges::Vector{Edge})
    keys = Tuple{Int, Int}[]
    members = Vector{Int}[]
    for (i, e) in enumerate(edges)
        key = minmax(e.u, e.v)
        k = findfirst(==(key), keys)
        k === nothing ? (push!(keys, key); push!(members, [i])) : push!(members[k], i)
    end
    length(unique(components(G, keys))) == 1 ||
        error("the steady-state steps leave the segment graph disconnected")
    adj = [Tuple{Int, Int}[] for _ in 1:G]
    for (k, (a, b)) in enumerate(keys)
        push!(adj[a], (b, k)); push!(adj[b], (a, k))
    end
    W(k, a) = sum(term(PI, wout(edges[i], a)) for i in members[k])   # merged weight out of a
    out = Blocks(Vector{Int}[], Vector{Int}[], Dict{Int, PI}[], Vector{Int}[])
    for blk in bicomponents(G, adj)
        nodes = sort!(unique([x for k in blk for x in keys[k]]))
        loc = Dict(g => i for (i, g) in enumerate(nodes))
        trees = spanning_trees(length(nodes), [(loc[keys[k][1]], loc[keys[k][2]]) for k in blk])
        factor = Dict{Int, PI}()
        for r in nodes
            factor[r] = sum(expand(PI, oriented(t, blk, keys, r, W)) for t in trees)
        end
        rest = components(G, [keys[k] for k in eachindex(keys) if k ∉ blk])
        attach = [nodes[findfirst(a -> rest[a] == rest[g], nodes)] for g in 1:G]
        push!(out.nodes, nodes); push!(out.edges, reduce(vcat, members[blk]))
        push!(out.factor, factor); push!(out.attach, attach)
    end
    out
end

"""The weights of tree `t` (positions into the block's keys `blk`) oriented toward `r`."""
function oriented(t, blk, keys, r, W)
    ws, seen, frontier = PI[], Set([r]), [r]
    while !isempty(frontier)
        p = popfirst!(frontier)
        for i in t
            a, b = keys[blk[i]]
            c = a == p && b ∉ seen ? b : b == p && a ∉ seen ? a : 0
            c == 0 && continue
            push!(ws, W(blk[i], c)); push!(seen, c); push!(frontier, c)
        end
    end
    ws
end

"""Tree sum of the segment graph toward `g`: the product over blocks of each block's tree
sum toward its node nearest to `g` (Volkenstein–Goldstein rule d)."""
tree_sum(B::Blocks, g) =
    expand(PI, [B.factor[b][B.attach[b][g]] for b in eachindex(B.nodes)])

"""Simple cycles of a block (its edges `eidx` into `edges`; parallel edges distinct), each
as the edge positions it uses, from XOR combinations of the fundamental cycles."""
function simple_cycles(nodes, eidx, edges)
    parent = Dict{Int, Int}(nodes[1] => 0)      # node => tree edge position to its parent
    queue, tree = [nodes[1]], falses(length(eidx))
    while !isempty(queue)
        a = popfirst!(queue)
        for (i, j) in enumerate(eidx)
            e = edges[j]
            b = e.u == a ? e.v : e.v == a ? e.u : 0
            (b == 0 || haskey(parent, b)) && continue
            parent[b] = i; tree[i] = true; push!(queue, b)
        end
    end
    up(a) = (e = edges[eidx[parent[a]]]; e.u == a ? e.v : e.u)
    path(a) = (p = Int[]; while parent[a] != 0; push!(p, parent[a]); a = up(a); end; p)
    fundamental = BitVector[]
    for i in findall(!, tree)
        c = falses(length(eidx)); c[i] = true
        for k in [path(edges[eidx[i]].u); path(edges[eidx[i]].v)]
            c[k] = !c[k]
        end
        push!(fundamental, c)
    end
    μ = length(fundamental)
    μ <= 20 || error("block with cycle rank $μ is too large to enumerate")
    cycles = Vector{Int}[]
    for mask in 1:(2^μ - 1)
        c = falses(length(eidx))
        for i in 1:μ
            (mask >> (i - 1)) & 1 == 1 && (c .⊻= fundamental[i])
        end
        used = findall(c)
        deg = Dict{Int, Int}()
        for i in used, x in (edges[eidx[i]].u, edges[eidx[i]].v)
            deg[x] = get(deg, x, 0) + 1
        end
        all(==(2), values(deg)) && length(orient(used, eidx, edges)) == length(used) &&
            push!(cycles, used)
    end
    cycles
end

"""The cycle through edge positions `used` walked from its first edge's `u`: each step as
`(edge index, forward)`. Stops early (shorter result) when the edges form several cycles."""
function orient(used, eidx, edges)
    e1 = edges[eidx[used[1]]]
    walk, left, at = [(eidx[used[1]], true)], used[2:end], e1.v
    while at != e1.u
        k = findfirst(i -> edges[eidx[i]].u == at || edges[eidx[i]].v == at, left)
        k === nothing && break
        e = edges[eidx[left[k]]]
        push!(walk, (eidx[left[k]], e.u == at)); at = e.u == at ? e.v : e.u
        deleteat!(left, k)
    end
    walk
end

"""The bare law `v = num/den` of a plain mechanism over `vars`: one variable per step constant
(`stepvars[j]`: `K#j` for an RE step, `kf#j`/`kr#j` for an SS step in its stored direction) and
one per metabolite (`concvars`); `dfree` is the free enzyme's share of `den`. `stats` holds the
structure and sizes: `segments`, `blocks`, `block_sizes`, `cycles`, `max_n`, term counts
(`den_terms`, `num_terms`, `sigma_terms`, `dsum_terms`, `factor_terms`, and `hill_terms` /
`omega_terms` before cancellation), `re_closing` (RE steps off the BFS tree, whose K the law
lacks), `silent_loops`, `dead_end` (SS steps whose speed the reduced law cannot carry), and the
factors themselves (`sigma`, `factors`, `attach`)."""
struct Law
    vars::Vars
    stepvars::Vector{Vector{Int}}
    concvars::Vector{Int}
    num::PI
    den::PI
    dfree::PI
    stats::Dict{Symbol, Any}
end

"""
    bare_rate_law(m::Pop.PMech) -> Law

`v = num/den` per unit enzyme, the net consumption of `m.substrates[1]`. RE segments
(rooted as the package roots them) contract to nodes weighted by `σ_g = Σ α_i`; the
denominator is `Σ_g σ_g D_g` with `D_g` the product of the block tree sums toward `g`; the
numerator is Hill's `Σ_κ n_κ (Π⁺ − Π⁻) Σ_κ` over the simple cycles κ of the segment
multigraph (self-loops included), checked equal to the edge sum `Σ_e ω_e J_e` when `check`.
Exponents of concentrations are brought to lowest terms.
"""
function bare_rate_law(m::Pop.PMech; check = true)
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

    # (1) RE segments, numbered by their lowest form, and RE weights from the segment root.
    n = length(m.species)
    seg_root = components(n, [(s.from, s.to) for s in m.steps if s.re])
    roots = unique(seg_root)
    segments = [findall(==(r), seg_root) for r in roots]
    seg = [findfirst(==(seg_root[i]), roots) for i in 1:n]
    hasres(sp) = !isempty(sp.res_added) || !isempty(sp.res_subtracted)
    alpha, seen, intree = [Mono() for _ in 1:n], falses(n), falses(length(m.steps))
    for S in segments
        r = argmin(i -> (length(m.species[i].bound), hasres(m.species[i]) ? 1 : 0,
                         string(m.species[i].name)), S)
        seen[r] = true; queue = [r]
        while !isempty(queue)
            cur = popfirst!(queue)
            for (j, s) in enumerate(m.steps)
                s.re || continue
                a, b, w = s.from == cur ? (s.from, s.to, ratio[j]) :
                          s.to == cur ? (s.to, s.from, mpow(ratio[j], -1)) : (0, 0, Mono())
                (a == 0 || seen[b]) && continue
                alpha[b] = mmul(alpha[a], w); seen[b] = intree[j] = true; push!(queue, b)
            end
        end
    end
    isconc(p) = p.first in values(conc)
    for (j, s) in enumerate(m.steps)
        (s.re && !intree[j]) || continue
        filter(isconc, mmul(alpha[s.to], mpow(alpha[s.from], -1))) == filter(isconc, ratio[j]) ||
            error("rate_equation: the rapid-equilibrium steps close a catalytic cycle " *
                  "(through $(m.species[s.from].name) ⇌ $(m.species[s.to].name)), so the " *
                  "mechanism has no finite rate. Make one step of the cycle steady-state.")
    end

    # (2) The segment multigraph of the SS steps; a step within one segment is a self-loop.
    x = conc[m.substrates[1]]
    copies(ms) = count(y -> conc[y.name] == x, ms)
    edges, loops = Edge[], Edge[]
    for (j, s) in enumerate(m.steps)
        s.re && continue
        kf, kr = stepvars[j]
        e = Edge(j, seg[s.from], seg[s.to],
                 mmul(mmul(mono(kf => 1), cmono(s.consumed)), alpha[s.from]),
                 mmul(mmul(mono(kr => 1), cmono(s.released)), alpha[s.to]),
                 copies(s.consumed) - copies(s.released) +
                 mexp(alpha[s.from], x) - mexp(alpha[s.to], x))
        push!(e.u == e.v ? loops : edges, e)
    end

    # (3)–(4) Block tree sums and the denominator Σ_g σ_g D_g.
    G = length(segments)
    B = blocks(G, edges)
    D = [tree_sum(B, g) for g in 1:G]
    sigma = [sum(term(PI, alpha[i]) for i in S) for S in segments]
    den = sum(sigma[g] * D[g] for g in 1:G)
    i_free = findfirst(Pop.is_free, m.species)
    dfree = i_free === nothing ? one(PI) : D[seg[i_free]]

    # (5a) Edge sum: Σ_e ω_e (fwd_e D_u − rev_e D_v), as the package forms it.
    num_omega, omega_terms = zero(PI), 0
    for e in [edges; loops]
        (check && e.ω != 0) || continue
        omega_terms += nterms(D[e.u]) + nterms(D[e.v])
        num_omega += scale(D[e.u] * e.fwd - D[e.v] * e.rev, e.ω)
    end
    # (5b) Hill: Σ_κ n_κ (Π⁺ − Π⁻) Σ_κ with Σ_κ the tree sum of the graph with κ contracted.
    num, hill_terms, ncycles, reactive_blocks, max_n = zero(PI), 0, 0, Set{Int}(), 0
    for b in eachindex(B.nodes), used in simple_cycles(B.nodes[b], B.edges[b], edges)
        ncycles += 1
        walk = orient(used, B.edges[b], edges)
        nκ = sum(f ? edges[i].ω : -edges[i].ω for (i, f) in walk)
        nκ == 0 && continue
        push!(reactive_blocks, b); max_n = max(max_n, abs(nκ))
        Pp = foldl(mmul, (f ? edges[i].fwd : edges[i].rev for (i, f) in walk); init = Mono())
        Pm = foldl(mmul, (f ? edges[i].rev : edges[i].fwd for (i, f) in walk); init = Mono())
        Σκ = contracted_tree_sum(G, edges, Set(x for (i, _) in walk for x in (edges[i].u, edges[i].v)))
        hill_terms += 2 * nterms(Σκ)
        num += scale(term(PI, Pp) - term(PI, Pm), nκ) * Σκ
    end
    for e in loops
        e.ω == 0 && continue
        hill_terms += 2 * nterms(D[e.u]); max_n = max(max_n, abs(e.ω))
        num += scale(term(PI, e.fwd) - term(PI, e.rev), e.ω) * D[e.u]
    end
    check && num != num_omega && error("Hill numerator differs from the edge-sum numerator")

    # (6) Concentrations to lowest terms.
    shift = lowest_shift((num, den), values(conc))
    stats = Dict{Symbol, Any}(
        :segments => G, :forms => n, :ss => length(edges) + length(loops),
        :re => count(s -> s.re, m.steps), :loops => length(loops),
        :silent_loops => [e.step for e in loops if e.ω == 0],
        :re_closing => [j for (j, s) in enumerate(m.steps) if s.re && !intree[j]],
        :blocks => length(B.nodes), :block_sizes => length.(B.nodes),
        :parallel => length(edges) - length(unique(minmax(e.u, e.v) for e in edges)),
        :cycles => ncycles, :den_terms => nterms(den), :num_terms => nterms(num),
        :sigma_terms => sum(nterms, sigma), :dsum_terms => sum(nterms, D), :max_n => max_n,
        :factor_terms => sum(nterms(f) for F in B.factor for f in values(F); init = 0),
        :hill_terms => hill_terms, :omega_terms => omega_terms,
        :dead_end => dead_end_steps(G, edges, loops, B, reactive_blocks),
        :sigma => sigma, :factors => B.factor, :attach => B.attach)
    Law(V, stepvars, collect(values(conc)), num * shift, den * shift, dfree * shift, stats)
end

"""Tree sum toward the contracted node of the segment graph with the segments `κ` merged."""
function contracted_tree_sum(G, edges, κ)
    c = minimum(κ)
    ids = Dict(g => i for (i, g) in enumerate(g for g in 1:G if g ∉ κ || g == c))
    to(g) = ids[g ∈ κ ? c : g]
    E = Edge[Edge(e.step, to(e.u), to(e.v), e.fwd, e.rev, e.ω) for e in edges
             if to(e.u) != to(e.v)]
    tree_sum(blocks(length(ids), E), ids[c])
end

"""SS steps whose speed the law cannot carry: a self-loop with no turnover, and every step
of a branch hanging off an articulation point that holds no turnover cycle (no reactive
cycle in its blocks, no reactive self-loop on its segments)."""
function dead_end_steps(G, edges, loops, B, reactive_blocks)
    out = Set(e.step for e in loops if e.ω == 0)
    for a in 1:G
        count(b -> a in B.nodes[b], eachindex(B.nodes)) >= 2 || continue
        comp = components(G, [(e.u, e.v) for e in edges if e.u != a && e.v != a])
        for r in unique(comp[g] for g in 1:G if g != a)
            branch = Set(g for g in 1:G if g != a && comp[g] == r)
            inside(b) = all(g -> g == a || g in branch, B.nodes[b])
            any(b -> inside(b) && b in reactive_blocks, eachindex(B.nodes)) && continue
            any(e -> e.u in branch && e.ω != 0, loops) && continue
            union!(out, (e.step for e in edges if (e.u in branch || e.u == a) &&
                                                   (e.v in branch || e.v == a)))
            union!(out, (e.step for e in loops if e.u in branch))
        end
    end
    sort!(collect(out))
end

"""The law with each step variable renamed to the package's group-representative name
(`m.step_names`), concentrations brought back to lowest terms: `(num, den, dfree)`."""
function grouped(L::Law, m::Pop.PMech)
    sub = Dict(v => mono(var!(L.vars, m.step_names[j][i]) => 1)
               for (j, vs) in enumerate(L.stepvars) for (i, v) in enumerate(vs))
    num, den, dfree = (subst(p, v -> get(sub, v, nothing)) for p in (L.num, L.den, L.dfree))
    shift = lowest_shift((num, den), L.concvars)
    num * shift, den * shift, dfree * shift
end

end # module
