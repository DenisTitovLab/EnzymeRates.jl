# ABOUTME: The joint log-linear solve over bare per-step constants: cycle and kinetic-group
# ABOUTME: rows of both conformations, a content priority, exact RREF, dependents and folds.

"""A bare constant: slot `slot` of edge `edge`'s constants (`constant_names` order), in the
inactive conformation when `inactive`."""
struct Col
    edge::Int
    slot::Int
    inactive::Bool
end

"""A dependent column as `Keq^keq · Π independent^exponent`, labelled `:kinetic_group` (a
group member tied to its class's representative), `:haldane` (Keq on its right side) or
`:wegscheider`."""
struct Dep
    col::Int
    terms::Vector{Pair{Int, Exp}}
    keq::Exp
    label::Symbol
end

"""
The joint solve of a mechanism. `cols` lists every column in elimination order and
`index` maps every bare constant to its column (an `:EqualAI` edge's inactive constant to
its active one). `indep` holds the independents that reach the law, in parameter order
(active, then inactive; by edge and slot), and `deps` every dependent that reaches it.
`fold` maps each dependent whose right side is one independent to that independent.
"""
struct Solution
    cols::Vector{Col}
    index::Dict{Col, Int}
    indep::Vector{Int}
    deps::Vector{Dep}
    fold::Dict{Int, Int}
end

"""
    inactive_mask(m) -> BitVector

The edges of the inactive conformation's law: the edges outside `:OnlyA` groups that the
flip form, the first free form by content, reaches through such edges (deletion plus
pruning on the single graph). No edge for a plain mechanism.
"""
function inactive_mask(m::Mechanism)
    is_allosteric(m) || return falses(length(m.edges))
    flip = findfirst(is_free, m.forms)
    flip === nothing && error("an allosteric mechanism flips its conformation at a free " *
                              "form (no bound metabolite, no residual); this one has none")
    keep = [m.tags[g] !== :OnlyA for g in m.group]
    seen = falses(length(m.forms))
    seen[flip] = grown = true
    while grown
        grown = false
        for (k, e) in enumerate(m.edges)
            keep[k] && seen[e.a] != seen[e.b] || continue
            seen[e.a] = seen[e.b] = grown = true
        end
    end
    BitVector(keep[k] && seen[e.a] for (k, e) in enumerate(m.edges))
end

"""The kinetic class of column `c`: its group's constant in its slot and state, except that
each member of an `:OnlyA` group keeps its own inactive class."""
class(m::Mechanism, c::Col) = c.inactive && m.tags[m.group[c.edge]] === :OnlyA ?
    (c.edge, c.slot, true, true) : (Int(m.group[c.edge]), c.slot, c.inactive, false)

"""
The columns of `m` in elimination order: inactive columns outside the law (`mask`) first,
then the other inactive columns, then the active ones; within each block the larger key
`(kind, occupied, ss, against)` first, ties broken by the oriented edge's content, larger
first. `kind` is 3 for an isomerization, 1 for a plain binding, 2 for other chemistry;
`occupied` holds when both forms of the step hold a metabolite; `against` marks the rate
constant that drives against the net reaction: the reverse of the physical-forward
direction, except that a binding of a product, whose taking up reverses the reaction,
drives against it when it takes the product up.
"""
function columns(m::Mechanism, fwd, mask)
    cols = Col[]
    for k in eachindex(m.edges), inactive in (false, true)
        inactive && !(is_allosteric(m) && m.tags[m.group[k]] !== :EqualAI) && continue
        foreach(s -> push!(cols, Col(k, s, inactive)), 1:(m.edges[k].re ? 1 : 2))
    end
    takes_product(o) = is_binding(o) &&
        any(i -> count_at(o.cons, i) > 0 && !m.mets[i].modifier &&
                 m.mets[i].name in m.products, eachindex(m.mets))
    keys = map(cols) do c
        o = oriented(m, c.edge, fwd)
        kind = o.cons == 0 && o.rel == 0 ? 3 : is_plain_binding(m, o) ? 1 : 2
        occupied = m.forms[o.a].bound != 0 && m.forms[o.b].bound != 0
        against = !o.re && (c.slot == 2) != takes_product(o)
        block = !c.inactive ? 0 : mask[c.edge] ? 1 : 2
        (block, kind, occupied, !o.re, against, o.a, o.b, o.cons, o.rel)
    end
    cols[sortperm(keys; rev = true)]
end

"""
    cycles(m, fwd) -> (zs, multiples)

The fundamental cycles of a spanning tree of `m`'s graph, each as signed edge counts in
the edges' physical-forward directions, and the multiple of the net reaction each runs
(an error when a cycle's net change is no multiple of it). A modifier copy counts as its
reactant.
"""
function cycles(m::Mechanism, fwd)
    E = length(m.edges)
    path = [zeros(Int, E) for _ in m.forms]                  # signed edges from form 1
    seen, tree, grown = falses(length(m.forms)), falses(E), true
    seen[1] = true
    while grown
        grown = false
        for (k, e) in enumerate(m.edges)
            seen[e.a] == seen[e.b] && continue
            from, to, s = seen[e.a] ? (e.a, e.b, 1) : (e.b, e.a, -1)
            path[to] .= path[from]
            path[to][k] += s
            seen[to] = tree[k] = grown = true
        end
    end
    names = unique!([x.name for x in m.mets])
    at = [searchsortedfirst(names, x.name) for x in m.mets]
    net = [x in m.substrates ? -1 : x in m.products ? 1 : 0 for x in names]
    first_net = findfirst(!=(0), net)
    zs, multiples = Vector{Int}[], Int[]
    for (k, e) in enumerate(m.edges)
        tree[k] && continue
        z = path[e.a] - path[e.b]
        z[k] += 1
        tot = zeros(Int, length(names))
        for (j, f) in enumerate(m.edges), i in eachindex(m.mets)
            tot[at[i]] += z[j] * (count_at(f.rel, i) - count_at(f.cons, i))
        end
        c = tot[first_net] * net[first_net]
        tot == c .* net ||
            error("the cycle closed by $(step_text(m, e)) changes the metabolites by no " *
                  "multiple of the net reaction")
        push!(zs, [fwd[j] ? z[j] : -z[j] for j in 1:E])
        push!(multiples, c)
    end
    zs, multiples
end

"""
    constraint_system(m, mask = inactive_mask(m)) -> (cols, A, b)

The joint system `A·log(constants) = b·log(Keq)` over the columns `cols` (elimination
order). Rows: each fundamental cycle in the active state and, for an allosteric mechanism,
again over the inactive constants (an `:EqualAI` edge reads its active column); each
member of a kinetic class equal to its first member, slot by slot. A cycle row sums, over
its edges in their physical-forward direction, the log of the forward equilibrium ratio:
`kf/kr` at steady state, `1/K` for a binding's dissociation constant, `K` otherwise.
"""
function constraint_system(m::Mechanism, mask = inactive_mask(m))
    fwd = orient(m)
    cols = columns(m, fwd, mask)
    index = Dict(c => i for (i, c) in enumerate(cols))
    rows, b = Vector{Int}[], Int[]
    zs, multiples = cycles(m, fwd)
    states = is_allosteric(m) ? (false, true) : (false,)
    for inactive in states, (z, c) in zip(zs, multiples)
        row = zeros(Int, length(cols))
        col(k, s) = index[Col(k, s, inactive && m.tags[m.group[k]] !== :EqualAI)]
        for (k, zk) in enumerate(z)
            zk == 0 && continue
            o = oriented(m, k, fwd)
            row[col(k, 1)] += o.re && is_binding(o) ? -zk : zk
            o.re || (row[col(k, 2)] -= zk)
        end
        push!(rows, row)
        push!(b, c)
    end
    lead = Dict{Tuple{Int, Int, Bool, Bool}, Int}()
    for (i, c) in enumerate(cols)
        f = get!(lead, class(m, c), i)
        f == i && continue
        row = zeros(Int, length(cols))
        row[i], row[f] = 1, -1
        push!(rows, row)
        push!(b, 0)
    end
    A = isempty(rows) ? zeros(Int, 0, length(cols)) : permutedims(reduce(hcat, rows))
    cols, A, b
end

"""Reduced row echelon form of `M` in place over its first `n` columns; returns the pivot
columns. `Rational{Int}` arithmetic throws `OverflowError` instead of wrapping."""
function rref!(M::Matrix{<:Rational}, n::Int)
    piv = Int[]
    for col in 1:n
        r = length(piv) + 1
        r > size(M, 1) && break
        p = findfirst(i -> !iszero(M[i, col]), r:size(M, 1))
        p === nothing && continue
        p += r - 1
        for k in axes(M, 2)
            M[r, k], M[p, k] = M[p, k], M[r, k]
        end
        f = inv(M[r, col])
        foreach(k -> M[r, k] *= f, axes(M, 2))
        for i in axes(M, 1)
            (i == r || iszero(M[i, col])) && continue
            g = M[i, col]
            foreach(k -> M[i, k] -= g * M[r, k], axes(M, 2))
        end
        push!(piv, col)
    end
    piv
end

"""
    solve(m) -> Solution

The joint solve of `m`'s constraint system (`constraint_system`): the exact reduced row
echelon form in checked `Rational{Int}`, retried in `Rational{BigInt}` on overflow. Pivot
columns are dependents, the others independent; only columns of the law are kept, and
columns outside it, eliminated first, never enter a kept right side.
"""
function solve(m::Mechanism)
    mask = inactive_mask(m)
    cols, A, b = constraint_system(m, mask)
    n = length(cols)
    M = Exp.([A b])
    piv = try
        rref!(M, n)
    catch err
        err isa OverflowError || rethrow()
        M = Rational{BigInt}.([A b])
        rref!(M, n)
    end
    all(iszero, view(M, length(piv) + 1:size(M, 1), n + 1)) ||
        error("the thermodynamic relations contradict each other (0 = log Keq)")
    index = Dict(c => i for (i, c) in enumerate(cols))
    for (i, c) in enumerate(cols)
        is_allosteric(m) && !c.inactive && get!(index, Col(c.edge, c.slot, true), i)
    end
    rep = Dict(class(m, c) => i for (i, c) in enumerate(cols))     # last member kept
    in_law(i) = !cols[i].inactive || mask[cols[i].edge]
    free = setdiff(1:n, piv)
    deps = Dep[]
    for (r, p) in enumerate(piv)
        in_law(p) || continue
        terms = [f => -Exp(M[r, f]) for f in free if !iszero(M[r, f])]
        keq = Exp(M[r, n + 1])
        label = rep[class(m, cols[p])] != p ? :kinetic_group :
                iszero(keq) ? :wegscheider : :haldane
        push!(deps, Dep(p, terms, keq, label))
    end
    indep = sort!(filter(in_law, free);
                  by = i -> (cols[i].inactive, cols[i].edge, cols[i].slot))
    fold = Dict(d.col => d.terms[1].first for d in deps
                if iszero(d.keq) && length(d.terms) == 1 && d.terms[1].second == 1)
    Solution(cols, index, indep, deps, fold)
end
