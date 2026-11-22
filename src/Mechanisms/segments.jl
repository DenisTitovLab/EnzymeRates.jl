# ABOUTME: The segments and cycles of a step mask: rapid-equilibrium segments with the
# ABOUTME: steady-state steps between them, and turnover cycles.

# Sub-graphs: a mask says which steps of a mechanism a sub-graph keeps, such as an
# allosteric mechanism's inactive conformation. Nothing here constructs a mechanism, so no
# check meant for the mechanisms a user writes fires on a sub-graph.

"""
    _present_species(m, on) -> BitVector

The forms of `m` that the steps `on` keeps touch, one flag per form: a sub-graph's law,
segments and network cover these.
"""
function _present_species(m::CatalyticScheme, on::BitVector)
    present = falses(length(m.enzyme_forms))
    for (k, e) in enumerate(m.steps)
        on[k] && (present[e.from] = present[e.to] = true)
    end
    present
end

"""
    rapid_equilibrium_segments(m, on, re = <m's rapid-equilibrium flags>) -> (; segment_of,
        segments, offsets, steps)

The rapid-equilibrium segments of the steps of `m` that `on` keeps, joined by the kept
steps that `re` flags. `segments[s]` lists segment `s`'s forms breadth-first from its first
form, each form's steps in the order of the groups, and the segments are numbered by their
first forms; `segment_of[f]` is form `f`'s segment, 0 for a form no kept step touches
(`_present_species`). `offsets[n, f]` is how many more of metabolite `n` (an index into
`m.metabolites`) form `f` carries than the form of its segment that carries fewest of it.
Within a segment, rapid equilibrium fixes each form's weight relative to any other as a
monomial in concentrations, one factor per metabolite taken up or given off along the path
between them; chemistry passes the counts through unchanged.

`steps` lists the kept steps that `re` leaves at steady state and that join two segments,
in step order: the edges of the steady-state graph between the segments.
"""
function rapid_equilibrium_segments(m::CatalyticScheme, on::BitVector,
                                    re::BitVector = BitVector(s.rapid_equilibrium
                                                              for s in m.steps))
    nf = length(m.enzyme_forms)
    # The kept rapid-equilibrium steps at each form, group by group: `(k, 1)` at the form
    # step `k` leaves in its stored direction, `(k, -1)` at the form it leaves in reverse.
    out = [Tuple{Int, Int}[] for _ in 1:nf]
    for k in reduce(vcat, _groups(m.steps); init = Int[])
        on[k] && re[k] || continue
        push!(out[m.steps[k].from], (k, 1))
        push!(out[m.steps[k].to], (k, -1))
    end
    present = _present_species(m, on)
    segment_of, segments = zeros(Int, nf), Vector{Int}[]
    offsets = zeros(Int, length(m.metabolites), nf)
    for root in 1:nf
        present[root] && segment_of[root] == 0 || continue
        members = [root]
        push!(segments, members)
        segment_of[root] = length(segments)
        head = 0
        while head < length(members)
            head += 1
            u = members[head]
            for (k, direction) in out[u]
                e = m.steps[k]
                v = direction == 1 ? e.to : e.from
                segment_of[v] == 0 || continue
                segment_of[v] = length(segments)
                push!(members, v)
                for n in axes(offsets, 1)
                    offsets[n, v] = offsets[n, u]
                end
                for x in e.taken_up
                    offsets[x, v] += direction
                end
                for x in e.given_off
                    offsets[x, v] -= direction
                end
            end
        end
    end
    for segment in segments, n in axes(offsets, 1)
        low = minimum(f -> offsets[n, f], segment)
        for f in segment
            offsets[n, f] -= low
        end
    end
    steps = [k for (k, e) in enumerate(m.steps)
             if on[k] && !re[k] && segment_of[e.from] != segment_of[e.to]]
    (; segment_of, segments, offsets, steps)
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
    has_turnover_cycle(m::CatalyticScheme, on::BitVector) -> Bool

Whether a cycle of the steps of `m` that `on` keeps turns the reaction over: no potential
of the forms rises along each step by its uptake of the first substrate, which a weighted
union-find finds (`has_unbalanced_cycle`). A cycle of a
mechanism whose steps conserve atoms runs the reaction a whole number of times, so its
uptake of the first substrate is zero only when it turns nothing over; a hand-written
cycle that does not is read by that uptake alone. The guard's class 1 rejects a
mechanism whose steps hold no turnover cycle, as its rate is zero everywhere; a turnover
cycle of rapid-equilibrium steps makes the rate infinite (`_has_vmax`).
"""
function has_turnover_cycle(m::CatalyticScheme, on::BitVector)
    # The first substrate, whose net uptake is the rate.
    s1 = m.substrates[1]
    uptake(list::Vector{Int}) = count(==(s1), list)
    n = length(m.enzyme_forms)
    has_unbalanced_cycle([e.from for e in m.steps], [e.to for e in m.steps],
                         [uptake(e.taken_up) - uptake(e.given_off) for e in m.steps],
                         on, zeros(Int, n), zeros(Int, n))
end
