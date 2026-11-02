# ABOUTME: The segments and cycles of a step mask: rapid-equilibrium segments with the
# ABOUTME: steady-state steps between them, and turnover cycles.

# Sub-graphs: a mask says which steps of a mechanism a sub-graph keeps, such as an
# allosteric mechanism's inactive conformation. Nothing here constructs a mechanism, so no
# check meant for the mechanisms a user writes fires on a sub-graph.

"""
    present_species(m, on) -> BitVector

The nodes of `m` that the steps `on` keeps touch, one flag per node: a sub-graph's law,
segments and network cover these.
"""
function present_species(m::Mechanism, on::BitVector)
    present = falses(length(m.nodes))
    for (k, e) in enumerate(m.steps)
        on[k] && (present[e.from] = present[e.to] = true)
    end
    present
end

"""The distinct metabolite names of `m` (the names of `_table`), sorted: a reactant and its
modifier copy share one, as they share a concentration."""
function _distinct_names(m::Mechanism)
    names = Symbol[]
    for n in m.nodes
        foreach(x -> push!(names, x.name), n.bound)
        append!(names, n.added)
        append!(names, n.subtracted)
    end
    for e in m.steps
        foreach(x -> push!(names, x.name), e.taken_up)
        foreach(x -> push!(names, x.name), e.given_off)
    end
    unique!(sort!(names))
end

"""
    segment_walk(present, out, to, change) -> (; segment_of, segments, offsets, parent,
                                                disagreeing)

The segments of a graph of nodes and transitions, walked breadth-first along the
transitions from each `present` node no earlier walk reached, in node order: `out[u]`
lists the transitions that leave node `u`, in the order the walk takes them; transition
`t` runs to node `to[t]` and takes up `change[n, t]` copies of metabolite name `n` net.
`segment_of[u]` is node `u`'s segment (0 for a node not present), `segments` each one's
nodes in walk order, `offsets[n, u]` the copies of name `n` the walk takes up from the
segment's first node to `u`, `parent[u]` the transition that reached `u` (0 at a first
node), and `disagreeing` the first transition whose offsets disagree, 0 when none does.
Such a transition closes, with its tree path, a cycle that takes up some name net; when
there is none, each offset is a potential and does not depend on the walk.
"""
function segment_walk(present::BitVector, out::Vector{Vector{Int}}, to::Vector{Int},
                      change::Matrix{Int})
    nn = length(present)
    segment_of, segments, parent = zeros(Int, nn), Vector{Int}[], zeros(Int, nn)
    offsets, disagreeing = zeros(Int, size(change, 1), nn), 0
    for root in 1:nn
        present[root] && segment_of[root] == 0 || continue
        members = [root]
        push!(segments, members)
        segment_of[root] = length(segments)
        head = 0
        while head < length(members)
            head += 1
            u = members[head]
            for t in out[u]
                v = to[t]
                if segment_of[v] == 0
                    segment_of[v], parent[v] = length(segments), t
                    push!(members, v)
                    for n in axes(offsets, 1)
                        offsets[n, v] = offsets[n, u] + change[n, t]
                    end
                elseif disagreeing == 0 &&
                       any(n -> offsets[n, v] != offsets[n, u] + change[n, t],
                           axes(offsets, 1))
                    disagreeing = t
                end
            end
        end
    end
    (; segment_of, segments, offsets, parent, disagreeing)
end

"""
    rapid_equilibrium_segments(m, on, re = <m's rapid-equilibrium flags>) -> (; segment_of,
        segments, offsets, bottom, steps)

The rapid-equilibrium segments of the steps of `m` that `on` keeps, joined by the kept
steps that `re` flags (`segment_walk`). `segments[s]` lists segment `s`'s forms
breadth-first from its first form, each form's steps in the order of the groups, and the
segments are numbered by their first forms; `segment_of[f]` is form `f`'s segment, 0 for a
form no kept step touches (`present_species`). `offsets[n, f]` is how many more of
metabolite name `n` form `f` carries than the form of its segment that carries fewest of
it, `n` being the position of the name among the mechanism's distinct metabolite names,
in sorted order. Within a segment, rapid equilibrium fixes each form's weight relative to
any other as a monomial in concentrations, one factor per metabolite taken up or given off
along the path between them; chemistry passes the counts through unchanged.

`steps` lists the kept steps that `re` leaves at steady state and that join two segments,
in step order: the edges of the steady-state graph between the segments.

`bottom` lists the forms of the first segment whose weights all vanish at zero products,
or all vanish at zero substrates; it is empty when no segment does. Each metabolite name
takes its side from `m`'s substrates and products, so a modifier copy, which shares its
reactant's name and concentration, counts on that reactant's side even in a sub-graph
that holds no step taking up or giving off the reactant itself. Once the weights are
scaled to clear every negative exponent, a segment normally has a bottom form whose
weight is free of one side's metabolites — the form the others drain to when that side is
absent. A segment with no such form (e.g. `{E(P), E(Q), E(P, Q)}` with `E + P` and
`E + Q` at steady state: weights `Kq·P : Kp·Q : P·Q`) is empty at zero products, and the
rate there depends on how `E(P, Q)` splits between `E(P)` and `E(Q)` — a ratio of two
fast release rates the RE approximation does not carry. The rate equation is `0/0`
exactly where initial-rate data sit. A segment whose weights vanish only at a corner
mixing a substrate and a product (a mixed abortive complex, or the ping-pong `E` /
`E(; residual)` pair) is harmless: no turnover is possible at that corner.
"""
function rapid_equilibrium_segments(m::Mechanism, on::BitVector,
                                    re::BitVector = BitVector(s.rapid_equilibrium
                                                              for s in m.steps))
    nf = length(m.nodes)
    names = _distinct_names(m)
    nn = length(names)
    row(x::Metabolite) = searchsortedfirst(names, x.name)
    # Each kept rapid-equilibrium step as two transitions, 2j - 1 in its stored direction
    # and 2j in reverse, each listed at the form it leaves.
    kept = [k for k in reduce(vcat, m.groups; init = Int[]) if on[k] && re[k]]
    out, to = [Int[] for _ in 1:nf], zeros(Int, 2 * length(kept))
    change = zeros(Int, nn, 2 * length(kept))
    for (j, k) in enumerate(kept)
        e = m.steps[k]
        t = 2j - 1
        to[t], to[t + 1] = e.to, e.from
        push!(out[e.from], t)
        push!(out[e.to], t + 1)
        for x in e.taken_up
            change[row(x), t] += 1
        end
        for x in e.given_off
            change[row(x), t] -= 1
        end
        for n in 1:nn
            change[n, t + 1] = -change[n, t]
        end
    end
    present = present_species(m, on)
    (; segment_of, segments, offsets) = segment_walk(present, out, to, change)
    for segment in segments, n in axes(offsets, 1)
        low = minimum(f -> offsets[n, f], segment)
        for f in segment
            offsets[n, f] -= low
        end
    end
    # Whether each metabolite name is a product, and whether a substrate.
    products, substrates = falses(nn), falses(nn)
    for (n, x) in enumerate(names)
        products[n] = insorted(x, m.products)
        substrates[n] = insorted(x, m.substrates)
    end
    bottomless(segment, side) = length(segment) >= 2 &&
        all(f -> any(n -> side[n] && offsets[n, f] > 0, 1:nn), segment)
    i = findfirst(seg -> bottomless(seg, products) || bottomless(seg, substrates), segments)
    steps = [k for (k, e) in enumerate(m.steps)
             if on[k] && !re[k] && segment_of[e.from] != segment_of[e.to]]
    (; segment_of, segments, offsets, bottom = i === nothing ? Int[] : segments[i], steps)
end

"""Whether the steps `k` with `on[k]` close a cycle of nonzero weight, step `k` joining
forms `from[k]` and `to[k]` with weight `weight[k]`: giving every form a potential that
rises by each step's weight along it then fails somewhere, which a weighted union-find in
the buffers `parent` and `offset` (one entry per form) detects."""
function has_unbalanced_cycle(from, to, weight, on, parent, offset)
    parent .= eachindex(parent); fill!(offset, 0)
    function root(x)
        p = 0
        while parent[x] != x
            p += offset[x]; x = parent[x]
        end
        x, p
    end
    for k in eachindex(from)
        on[k] || continue
        ra, pa = root(from[k]); rb, pb = root(to[k])
        if ra == rb
            pb == pa + weight[k] || return true
        else
            parent[rb] = ra; offset[rb] = pa + weight[k] - pb
        end
    end
    false
end

"""
    has_turnover_cycle(m::Mechanism, on::BitVector) -> Bool

Whether a cycle of the steps of `m` that `on` keeps turns the reaction over: no potential
of the forms rises along each step by its uptake of the first substrate (a modifier copy
counting as it), which a weighted union-find finds (`has_unbalanced_cycle`). A cycle of a
mechanism whose steps conserve atoms runs the reaction a whole number of times, so its
uptake of the first substrate is zero only when it turns nothing over; a hand-written
cycle that does not is read by that uptake alone. The guard's class 1 rejects a
mechanism whose steps hold no turnover cycle, as its rate is zero everywhere; a turnover
cycle of rapid-equilibrium steps makes the rate infinite (`_has_vmax`).
"""
function has_turnover_cycle(m::Mechanism, on::BitVector)
    # The first substrate, whose net uptake, a modifier copy counting as it, is the rate.
    s1 = m.substrates[1]
    uptake(list::Vector{Metabolite}) = count(x -> x.name === s1, list)
    n = length(m.nodes)
    has_unbalanced_cycle([e.from for e in m.steps], [e.to for e in m.steps],
                         [uptake(e.taken_up) - uptake(e.given_off) for e in m.steps],
                         on, zeros(Int, n), zeros(Int, n))
end
