# ABOUTME: The well-posedness check the mechanism macros, the parser, the enumerator's
# ABOUTME: exits and the derivation entries run; constructors check representation only.

"""
    _check_well_posed(m)

Error unless `m`'s rate law is the unique limit of mass-action kinetics, as its
rapid-equilibrium steps speed up and its tagged constants vanish, on every face where data
sit (spec section 8): the guard (`_guard`) finds no fault, and the error is its message
(`_guard_message`), which names each step at fault by its DSL line. It reads the reaction
graph through masks and builds nothing. The mechanism macros, the parser and every
derivation entry (`_derivation_input`) call it; the enumerator's exits and identify's base
tier keep what `_is_well_posed` accepts; the constructors and the Sig decoder call
neither, so code that builds a mechanism internally never trips a check meant for the
mechanisms it returns. Class 5 reads the derived law and runs where the law is derived
(`_assert_support`).
"""
function _check_well_posed(m::_AnyMech)
    v = _guard(m)
    v.fault === _NO_FAULT || error(_guard_message(m, v))
    nothing
end

"""Whether `_check_well_posed(m)` passes: the enumerator's exits and identify's base tier
keep a child exactly when it does."""
_is_well_posed(m::_AnyMech) = _guard(m).fault === _NO_FAULT

"""
    _bottomless_segment(m, mask, re) -> Vector{Int}

The forms, breadth-first (`_segments`), of the first rapid-equilibrium segment of the
edges `mask` keeps, joined by those `re` flags, whose weights all vanish at zero
products, or all vanish at zero substrates; empty when no segment does. Each metabolite
name takes its side from `m`'s substrates and products, so a modifier copy, which shares
its reactant's name and concentration, counts on that reactant's side even in a
sub-graph that holds no edge taking up or giving off the reactant itself (the inactive
conformation's, `_inactive_mask`).

Within a segment, rapid equilibrium fixes each form's weight relative to any
other as a monomial in concentrations (one factor per RE binding step on the
path between them). Once the weights are scaled to clear every negative
exponent, a segment normally has a bottom form whose weight is free of one
side's metabolites — the form the others drain to when that side is absent. A
segment with no such form (e.g. `{E(P), E(Q), E(P, Q)}` with `E + P` and
`E + Q` at steady state: weights `Kq·P : Kp·Q : P·Q`) is empty at zero
products, and the rate there depends on how `E(P, Q)` splits between `E(P)`
and `E(Q)` — a ratio of two fast release rates the RE approximation does not
carry. The rate equation is `0/0` exactly where initial-rate data sit. A
segment whose weights vanish only at a corner mixing a substrate and a product
(a mixed abortive complex, or the ping-pong `E` / `E(; residual)` pair) is
harmless: no turnover is possible at that corner.
"""
function _bottomless_segment(m::Mechanism, mask::BitVector, re::BitVector)
    s = _segments(m, mask, re)
    name_of = _name_index(m)
    # Whether each metabolite name (`_name_index`) is a product, and whether a substrate.
    products, substrates = falses(size(s.offsets, 1)), falses(size(s.offsets, 1))
    for (i, x) in enumerate(m.table)
        products[name_of[i]] = insorted(x.name, m.products)
        substrates[name_of[i]] = insorted(x.name, m.substrates)
    end
    for segment in s.segments, on_side in (products, substrates)
        length(segment) >= 2 &&
            all(f -> any(n -> on_side[n] && s.offsets[n, f] > 0, eachindex(on_side)),
                segment) && return segment
    end
    Int[]
end

"""`m`, once `_check_well_posed(m)` passes."""
_well_posed(m) = (_check_well_posed(m); m)

"""The concrete mechanism the compiled type `M` encodes (`_concrete`), checked by
`_check_well_posed`: the input of every derivation entry."""
_derivation_input(@nospecialize(M::Type{<:AbstractEnzymeMechanism})) =
    _well_posed(_concrete(M()))

# ─── The guard: verdicts, networks, faces, classes 1 and 2 ─────────────

"""
What the guard finds wrong with a mechanism (spec section 8): nothing (`_NO_FAULT`); no
cycle turns the reaction over (`_NO_CYCLE`, class 1); a cycle of rapid-equilibrium steps
turns it over (`_RE_TURNOVER`, class 2); on a face (class 3), a second closed class
(`_NONUNIQUE`), a rapid-equilibrium segment with no lowest form (`_BOTTOMLESS`), or a law
the limit test rejects where procedure A stopped at a trap fed on the face (`_TRAP`), a
segment whose lowest forms split (`_SPLIT`), a law segment chain with several closed
classes (`_LAW_MULTI`) or an aggregated chain with several (`_MULTILEVEL`); in admission's
interior (class 4, spec section 7), a tagged binding whose two sides the conformation's
other steps join (`_NO_LIMIT`), tagged bindings that close a directed cycle
(`_ARC_CYCLE`), forms free E's component does not reach along the tagged bindings
(`_ISLAND`), or an active conformation that no longer turns the reaction over
(`_NO_ACTIVE_TURNOVER`); on a face, in admission's face condition, a second closed class
of the ε-network (`_EPS_NONUNIQUE`), a stranded basin a step feeds that does not run back
from the law's forms (`_PUMPED`), or tagged steps that shorten the fed forms' cheapest
spanning in-tree where no potential absorbs their turnover (`_BRIDGE`); on a face of the
derived law (class 5), a denominator that loses every term (`_DENOMINATOR_VANISHES`).
"""
@enum _Fault::UInt8 begin
    _NO_FAULT
    _NO_CYCLE
    _RE_TURNOVER
    _NONUNIQUE
    _BOTTOMLESS
    _TRAP
    _SPLIT
    _LAW_MULTI
    _MULTILEVEL
    _NO_LIMIT
    _ARC_CYCLE
    _ISLAND
    _EPS_NONUNIQUE
    _PUMPED
    _BRIDGE
    _NO_ACTIVE_TURNOVER
    _DENOMINATOR_VANISHES
end

"""The faults admission reports (spec section 7); `_admission_message` writes their text."""
const _ADMISSION_FAULTS = (_NO_LIMIT, _ARC_CYCLE, _ISLAND, _NO_ACTIVE_TURNOVER,
                           _EPS_NONUNIQUE, _PUMPED, _BRIDGE)

"""
The limit test's verdict on a face (class 3): not run (`_LIMIT_UNTESTED`); the fast-RE
limit of the steady state equals the law (`_LIMIT_EQUAL`), differs from it
(`_LIMIT_DIFFER`) or moves with the relative speeds of the rapid-equilibrium steps
(`_LIMIT_SPEED`); the law has no value there (`_LIMIT_LAW_NAN`).
"""
@enum _LimitVerdict::UInt8 begin
    _LIMIT_UNTESTED
    _LIMIT_EQUAL
    _LIMIT_DIFFER
    _LIMIT_SPEED
    _LIMIT_LAW_NAN
end

"""
The guard's verdict on a mechanism, an isbits value, so a loop over verdicts keeps its
types concrete: its `fault`, the limit test's verdict where it decided (`limit`), the face
(an index into `_faces`; 0 when the fault belongs to no face), the conformation it lies in
(1 active, 2 inactive, 0 none) and the edge its message leads with (0 none).
"""
struct _Verdict
    fault::_Fault
    limit::_LimitVerdict
    face::Int
    conformation::Int
    edge::Int
end

"""The verdict on a mechanism that passes every check."""
const _WELL_POSED = _Verdict(_NO_FAULT, _LIMIT_UNTESTED, 0, 0, 0)

"""
A network the guard reads: `ncopies` copies of a mechanism's `nforms` forms, node
`(c - 1) * nforms + f` for form `f` in copy `c` (`_node`), `present` marking the nodes it
holds. Transitions come in pairs: `2j - 1` runs an edge (`edge`; 0 for the flip at free
E) in its stored direction in copy `conformation` (0 for the flip), and `2j` runs it in
reverse; `forward` marks the first of each pair. `fast` marks a rapid-equilibrium
transition and `vanishing` one whose rate goes to zero in the ε-network. `change[n, t]`
counts the copies of metabolite name `n` (`_name_index`) that transition `t` takes up less
those it gives off; bit `n - 1` of `takes[t]` is set when it takes up name `n`.
"""
struct _Network
    nforms::Int
    ncopies::Int
    present::BitVector
    from::Vector{Int}
    to::Vector{Int}
    edge::Vector{Int}
    conformation::Vector{Int}
    forward::BitVector
    fast::BitVector
    vanishing::BitVector
    change::Matrix{Int}
    takes::Vector{UInt64}
end

"""The node of form `f` in copy `c` of `net`."""
_node(net::_Network, f::Int, c::Int) = (c - 1) * net.nforms + f

"""Each node's transitions of `net` that `keep` flags and that leave it, in transition
order."""
function _out_transitions(net::_Network, keep::BitVector)
    out = [Int[] for _ in 1:(net.nforms * net.ncopies)]
    for t in eachindex(net.from)
        keep[t] && push!(out[net.from[t]], t)
    end
    out
end

"""Each node's successors in `net` along the transitions `keep` flags, in transition
order; its predecessors when `backward`."""
function _successors(net::_Network, keep::BitVector; backward::Bool = false)
    a, b = backward ? (net.to, net.from) : (net.from, net.to)
    next = [Int[] for _ in 1:(net.nforms * net.ncopies)]
    for t in eachindex(a)
        keep[t] && push!(next[a[t]], b[t])
    end
    next
end

"""The nodes that the nodes `start` reach along the successor lists `next`, `start`
included."""
function _reached(next::Vector{Vector{Int}}, start::Vector{Int})
    seen = falses(length(next))
    seen[start] .= true
    todo = copy(start)
    while !isempty(todo)
        for v in next[pop!(todo)]
            seen[v] || (seen[v] = true; push!(todo, v))
        end
    end
    seen
end

"""
The network over the `present` nodes holding both directions of edge `k` of `m` in copy
`c` wherever `use[k, c]` holds, fast where `fast[k, c]` does, the stored direction
vanishing where `vanish_forward[k, c]` does and the reverse where `vanish_reverse[k, c]`
does; with two copies, the fast flip between the copies of free E, which takes up and
gives off nothing.
"""
function _network(m::Mechanism, present::BitVector, use::BitMatrix, fast::BitMatrix,
                  vanish_forward::BitMatrix, vanish_reverse::BitMatrix)
    nf, ncopies = length(m.forms), size(use, 2)
    name_of = _name_index(m)
    nt = 2 * count(use) + (ncopies == 2 ? 2 : 0)
    from, to, edge, conformation = zeros(Int, nt), zeros(Int, nt), zeros(Int, nt),
                                   zeros(Int, nt)
    forward, isfast, vanishing = falses(nt), falses(nt), falses(nt)
    change, takes = zeros(Int, maximum(name_of; init = 0), nt), zeros(UInt64, nt)
    t = 0
    for c in 1:ncopies, (k, e) in enumerate(m.edges)
        use[k, c] || continue
        a, b = (c - 1) * nf + e.from, (c - 1) * nf + e.to
        for (u, v, gain, loss, gone) in
                ((a, b, e.consumed, e.released, vanish_forward[k, c]),
                 (b, a, e.released, e.consumed, vanish_reverse[k, c]))
            t += 1
            from[t], to[t], edge[t], conformation[t] = u, v, k, c
            forward[t], isfast[t], vanishing[t] = u == a, fast[k, c], gone
            for i in entries(gain)
                change[name_of[i], t] += 1
                takes[t] |= UInt64(1) << (name_of[i] - 1)
            end
            for i in entries(loss)
                change[name_of[i], t] -= 1
            end
        end
    end
    if ncopies == 2
        from[t + 1], to[t + 1], forward[t + 1], isfast[t + 1] = 1, nf + 1, true, true
        from[t + 2], to[t + 2], isfast[t + 2] = nf + 1, 1, true
    end
    _Network(nf, ncopies, present, from, to, edge, conformation, forward, isfast,
             vanishing, change, takes)
end

"""
The law network L of `m`, its rapid-equilibrium flags overridden by `re`: a plain
mechanism's forms and edges as they are; an allosteric mechanism's two copies, copy `c`
holding the edges of its law (`_law_mask`) and the forms they touch, free E always, joined
by the flip at free E.
"""
function _law_network(m::Mechanism, re::BitVector = _flags(m))
    n = length(m.edges)
    fast = falses(n, 1)
    fast[:, 1] .= re
    _network(m, trues(length(m.forms)), trues(n, 1), fast, falses(n, 1), falses(n, 1))
end
function _law_network(am::AllostericMechanism, re::BitVector = _flags(am.catalytic))
    m = am.catalytic
    n, nf = length(m.edges), length(m.forms)
    use, fast, present = falses(n, 2), falses(n, 2), falses(2 * nf)
    for c in 1:2
        mask = _law_mask(am, c)
        present[(c - 1) * nf + 1] = true
        for (k, e) in enumerate(m.edges)
            mask[k] || continue
            use[k, c], fast[k, c] = true, re[k]
            present[(c - 1) * nf + e.from] = present[(c - 1) * nf + e.to] = true
        end
    end
    _network(m, present, use, fast, falses(n, 2), falses(n, 2))
end

"""
The ε-network N of `m`, its rapid-equilibrium flags overridden by `re`: every form and
edge in each copy. In the copy that excludes an edge (`_excludes`), a rapid-equilibrium
plain binding stays fast while its uptake vanishes and its release keeps its rate; any
other excluded edge is slow and vanishes both ways. A plain mechanism's ε-network is its
law network.
"""
_eps_network(m::Mechanism, re::BitVector = _flags(m)) = _law_network(m, re)
function _eps_network(am::AllostericMechanism, re::BitVector = _flags(am.catalytic))
    m = am.catalytic
    n = length(m.edges)
    fast, vanish_forward, vanish_reverse = falses(n, 2), falses(n, 2), falses(n, 2)
    for c in 1:2, (k, e) in enumerate(m.edges)
        excluded, plain = _excludes(am, k, c), _is_plain_binding(m, e)
        fast[k, c] = re[k] && (!excluded || plain)
        vanish_forward[k, c] = excluded
        vanish_reverse[k, c] = excluded && !plain
    end
    _network(m, trues(2 * length(m.forms)), trues(n, 2), fast, vanish_forward,
             vanish_reverse)
end

"""
The faces the guard checks (spec section 8) of `m`, each a mask of the metabolite names it
zeroes (bit `n - 1` for name `n` of `_name_index`), in order: the interior (0), every
product, every substrate, each substrate and then each product alone on a side of two or
more; then each of these with every modifier-only name too. A reactant's `::Mod` copy
shares its reactant's name and so is zeroed with it; a regulatory-site ligand never is.
No face zeroes a substrate together with a product.
"""
function _faces(m::Mechanism)
    name_of = _name_index(m)
    bit(x::Symbol) = UInt64(1) << (name_of[searchsortedfirst(m.table, Met(x, false))] - 1)
    side(xs::Vector{Symbol}) = reduce(|, (bit(x) for x in xs); init = UInt64(0))
    base = UInt64[0, side(m.products), side(m.substrates)]
    length(m.substrates) >= 2 && append!(base, bit.(m.substrates))
    length(m.products) >= 2 && append!(base, bit.(m.products))
    mods = UInt64(0)
    for (i, x) in enumerate(m.table)
        x.modifier && !_is_reactant(m, x.name) && (mods |= UInt64(1) << (name_of[i] - 1))
    end
    mods == 0 ? base : [base; base .| mods]
end

"""Whether face `face` (`_faces`) zeroes metabolite name `n` of `_name_index`."""
_zeroes(face::UInt64, n::Int) = (face >> (n - 1)) & 1 == 1

"""The metabolite names that face `face` of `m` zeroes (`_faces`), in table order."""
function _face_names(m::Mechanism, face::UInt64)
    name_of = _name_index(m)
    Symbol[x.name for (i, x) in enumerate(m.table)
           if (i == 1 || name_of[i] != name_of[i - 1]) && _zeroes(face, name_of[i])]
end

"""The metabolite-name index (`_name_index`) of `m`'s first substrate, whose net uptake
is the rate."""
_turnover_name(m::Mechanism) =
    _name_index(m)[searchsortedfirst(m.table, Met(m.substrates[1], false))]

"""
Whether a cycle of the edges that `mask` keeps turns the reaction over: no potential of the
forms rises along each edge by its uptake of the first substrate (a modifier copy counting
as it), which a weighted union-find finds (`_has_unbalanced_cycle`). Class 1 rejects a
mechanism with none; its rate is zero everywhere.
"""
function _has_turnover_cycle(m::Mechanism, mask::BitVector)
    name_of, x1 = _name_index(m), _turnover_name(m)
    uptake(c::Counts) =
        sum((copies(c, i) for i in eachindex(m.table) if name_of[i] == x1); init = 0)
    n = length(m.forms)
    _has_unbalanced_cycle([e.from for e in m.edges], [e.to for e in m.edges],
                          [uptake(e.consumed) - uptake(e.released) for e in m.edges],
                          mask, zeros(Int, n), zeros(Int, n))
end

"""The tree transitions from node `x` up to its root in the forest `parent` (each node's
tree transition, 0 at a root)."""
function _tree_path(net::_Network, parent::Vector{Int}, x::Int)
    path = Int[]
    while parent[x] != 0
        t = parent[x]
        push!(path, t)
        x = net.from[t] == x ? net.to[t] : net.from[t]
    end
    path
end

"""The transitions of the cycle that the non-tree transition `t` between nodes `u` and
`v` closes in the forest `parent`, in cycle order from `u`."""
function _tree_cycle(net::_Network, parent::Vector{Int}, u::Int, v::Int, t::Int)
    up, down = _tree_path(net, parent, u), _tree_path(net, parent, v)
    while !isempty(up) && !isempty(down) && up[end] == down[end]
        pop!(up)
        pop!(down)
    end
    [up; reverse(down); t]
end

"""
    _fast_segments(net) -> (; segment, segments, offset, cycle)

The rapid-equilibrium segments of `net`, walked breadth-first along its fast transitions
from each present node no earlier walk reached: `segment[u]`, node `u`'s (0 for a node
`net` does not hold); `segments`, each one's nodes in walk order; `offset[n, u]`, the
copies of name `n` (`_name_index`) the walk takes up from the segment's first node to `u`;
and `cycle`, the transitions of the first fast cycle that takes up some name net, which
the first transition whose offsets disagree closes with its tree path, in cycle order.
When `cycle` is empty, each offset is a potential: it does not depend on the walk.
"""
function _fast_segments(net::_Network)
    nn = net.nforms * net.ncopies
    out = _out_transitions(net, net.fast)
    segment, segments, parent = zeros(Int, nn), Vector{Int}[], zeros(Int, nn)
    offset, cycle = zeros(Int, size(net.change, 1), nn), Int[]
    agrees(u, v, t) = all(n -> offset[n, v] == offset[n, u] + net.change[n, t],
                          axes(offset, 1))
    for root in 1:nn
        net.present[root] && segment[root] == 0 || continue
        members = [root]
        push!(segments, members)
        segment[root] = length(segments)
        head = 0
        while head < length(members)
            head += 1
            u = members[head]
            for t in out[u]
                v = net.to[t]
                if segment[v] == 0
                    segment[v], parent[v] = length(segments), t
                    push!(members, v)
                    offset[:, v] .= view(offset, :, u) .+ view(net.change, :, t)
                elseif isempty(cycle) && !agrees(u, v, t)
                    cycle = _tree_cycle(net, parent, u, v, t)
                end
            end
        end
    end
    (; segment, segments, offset, cycle)
end

"""The transitions of the first cycle of fast transitions of `net` that takes up some
metabolite name net (class 2: the rate is infinite), in cycle order; empty when none does
(`_fast_segments`). Class 2 runs once, at the interior, since a fast cycle that runs on a
face runs at the interior too."""
_re_turnover_transitions(net::_Network) = _fast_segments(net).cycle

"""
    _guard(m) -> _Verdict

The well-posedness guard (spec section 8) on `m`, plain or allosteric: class 1 (a cycle
turns the reaction over, `_has_turnover_cycle`); class 2 (no cycle of fast transitions
does, `_re_turnover_transitions`, on the ε-network, where an excluded rapid-equilibrium
binding stays fast while its affinity vanishes); for an allosteric mechanism, admission's
interior (`_admission_interior`); then on every face class 3 (`_check_face`, on the law
network) and, for an allosteric mechanism, admission's face condition (`_admission_face`,
on the ε-network), the first fault winning. `_WELL_POSED` when every check passes. It
reads the reaction graph through masks and builds nothing.
"""
function _guard(m::Mechanism)
    _has_turnover_cycle(m, trues(length(m.edges))) ||
        return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
    net = _law_network(m)
    v = _re_turnover_verdict(net)
    v === _WELL_POSED || return v
    _faces_verdict(m, net)
end
function _guard(am::AllostericMechanism)
    m = am.catalytic
    _has_turnover_cycle(m, trues(length(m.edges))) ||
        return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
    eps = _eps_network(am)
    v = _re_turnover_verdict(eps)
    v === _WELL_POSED || return v
    v = _admission_interior(am)
    v === _WELL_POSED || return v
    law, faces = _law_network(am), _faces(m)
    for i in eachindex(faces)
        v = _check_face(am, law, faces, i)
        v.fault === _NO_FAULT || return v
        v = _admission_face(am, eps, law, faces, i)
        v.fault === _NO_FAULT || return v
    end
    _WELL_POSED
end

"""Class 2's verdict on the ε-network `net`: `_RE_TURNOVER` in the conformation of the
first fast cycle that turns the reaction over, led by its first edge; else
`_WELL_POSED`."""
function _re_turnover_verdict(net::_Network)
    cycle = _re_turnover_transitions(net)
    isempty(cycle) && return _WELL_POSED
    _Verdict(_RE_TURNOVER, _LIMIT_UNTESTED, 0, net.conformation[cycle[1]],
             net.edge[cycle[1]])
end

# ─── Class 3: procedure A on every face ────────────────────────────────

"""Each node's strongly connected component in the digraph with successor lists `next`,
searched from the nodes `nodes` marks, numbered in the order Tarjan's algorithm closes
them; 0 for a node no search reaches."""
function _strong_components(next::Vector{Vector{Int}}, nodes::BitVector)
    n = length(next)
    index, low, comp = zeros(Int, n), zeros(Int, n), zeros(Int, n)
    onstack, stack, work = falses(n), Int[], Tuple{Int, Int}[]
    counter = ncomp = 0
    for s in 1:n
        nodes[s] && index[s] == 0 || continue
        counter += 1
        index[s] = low[s] = counter
        push!(stack, s)
        onstack[s] = true
        push!(work, (s, 1))
        while !isempty(work)
            u, i = work[end]
            if i <= length(next[u])
                work[end] = (u, i + 1)
                v = next[u][i]
                if index[v] == 0
                    counter += 1
                    index[v] = low[v] = counter
                    push!(stack, v)
                    onstack[v] = true
                    push!(work, (v, 1))
                elseif onstack[v]
                    low[u] = min(low[u], index[v])
                end
            else
                pop!(work)
                isempty(work) || (low[work[end][1]] = min(low[work[end][1]], low[u]))
                if low[u] == index[u]
                    ncomp += 1
                    while true
                        v = pop!(stack)
                        onstack[v] = false
                        comp[v] = ncomp
                        v == u && break
                    end
                end
            end
        end
    end
    comp
end

"""The closed classes (bottom strongly connected components, `_strong_components`) of the
digraph with successor lists `next`, over the nodes `nodes` marks, each as its sorted
nodes; a node with no successor is a class of its own."""
function _bottom_sccs(next::Vector{Vector{Int}}, nodes::BitVector)
    comp = _strong_components(next, nodes)
    ncomp = maximum(comp; init = 0)
    leaves = falses(ncomp)
    for u in eachindex(next), v in next[u]
        nodes[u] && comp[u] != comp[v] && (leaves[comp[u]] = true)
    end
    classes = [Int[] for _ in 1:ncomp]
    for u in eachindex(next)
        nodes[u] && push!(classes[comp[u]], u)
    end
    [classes[k] for k in 1:ncomp if !leaves[k]]
end

"""The closed classes of the transitions of `net` that `keep` marks, over its present
nodes."""
_closed_classes(net::_Network, keep::BitVector) =
    _bottom_sccs(_successors(net, keep), net.present)

"""The fast classes (`class_of`) that the alive fast transitions (`fast_next`) reach from
node `w`, its own included."""
function _classes_reached(fast_next::Vector{Vector{Int}}, class_of::Vector{Int}, w::Int)
    seen = falses(length(fast_next))
    seen[w] = true
    todo, found = [w], Int[]
    while !isempty(todo)
        x = pop!(todo)
        class_of[x] > 0 && !(class_of[x] in found) && push!(found, class_of[x])
        for y in fast_next[x]
            seen[y] || (seen[y] = true; push!(todo, y))
        end
    end
    found
end

"""
    _face_structure(net, face) -> NamedTuple

The structure procedure A reads on face `face` of the law network
`net`: `alive`, the transitions that take up no zeroed name; `closed`, the closed classes
of the alive ones; the rapid-equilibrium segments (`segment` per node, `segments`
breadth-first) and each node's uptake of every name from its segment's first node
(`offset`, `_fast_segments`), and `lowest`, the nodes at the segment's componentwise least
uptake of the zeroed names, whose weight stays
positive on the face; `bottomless`, the first segment with no lowest node (0 when none);
the fast classes (closed classes of the alive fast transitions, a node with no alive fast
exit its own), `class_of` each node's (0 for a node that drains elsewhere) and `trap`,
the classes of no lowest node; `aggregated`, the closed classes of the chain over fast
classes that an alive slow transition and then alive fast paths make; and `law_classes`,
the closed classes of the law's segment chain, its alive slow transitions out of lowest
nodes.
"""
function _face_structure(net::_Network, face::UInt64)
    nn, nt = net.nforms * net.ncopies, length(net.from)
    alive = BitVector([net.takes[t] & face == 0 for t in 1:nt])
    closed = _closed_classes(net, alive)
    zeroed = [n for n in axes(net.change, 1) if _zeroes(face, n)]
    (; segment, segments, offset) = _fast_segments(net)
    lowest, bottomless = falses(nn), 0
    for (g, members) in enumerate(segments)
        least = [minimum(offset[n, u] for u in members) for n in zeroed]
        for u in members
            lowest[u] = all(j -> offset[zeroed[j], u] == least[j], eachindex(zeroed))
        end
        bottomless == 0 && !any(u -> lowest[u], members) && (bottomless = g)
    end
    fast_alive = alive .& net.fast
    fast_classes = _closed_classes(net, fast_alive)
    class_of = zeros(Int, nn)
    for (k, C) in enumerate(fast_classes), u in C
        class_of[u] = k
    end
    all(C -> all(u -> lowest[u] == lowest[C[1]], C), fast_classes) ||
        error("_face_structure: a fast class holds lowest and other forms")
    trap = BitVector([!lowest[C[1]] for C in fast_classes])
    fast_next = _successors(net, fast_alive)
    aggregated_next = [Int[] for _ in fast_classes]
    law_next = [Int[] for _ in segments]
    for t in 1:nt
        alive[t] && !net.fast[t] || continue
        u, v = net.from[t], net.to[t]
        if class_of[u] > 0
            for k in _classes_reached(fast_next, class_of, v)
                k == class_of[u] || push!(aggregated_next[class_of[u]], k)
            end
        end
        lowest[u] && segment[u] != segment[v] && push!(law_next[segment[u]], segment[v])
    end
    aggregated = _bottom_sccs(aggregated_next, trues(length(fast_classes)))
    law_classes = _bottom_sccs(law_next, trues(length(segments)))
    (; alive, closed, segment, segments, offset, lowest, bottomless, fast_classes, class_of,
     trap, aggregated, law_classes)
end

"""
Where procedure A settles a face from its structure `s`
(`_face_structure`): a second closed class (`_NONUNIQUE`); a segment with no lowest form
(`_BOTTOMLESS`); `_NO_FAULT` by Case I (every segment holds one fast class, the law's
segment chain one closed class) or A2 (the aggregated chain has one closed class, holding
no trap, and each segment it meets holds one lowest class, the law's chain one closed
class); else the stage at which the limit test decides: the law's chain has several
closed classes (`_LAW_MULTI`), the aggregated chain has (`_MULTILEVEL`), its closed class
holds a trap (`_TRAP`), or a segment it meets splits (`_SPLIT`).
"""
function _face_stage(s)
    length(s.closed) > 1 && return _NONUNIQUE
    s.bottomless > 0 && return _BOTTOMLESS
    classes, lowest = zeros(Int, length(s.segments)), zeros(Int, length(s.segments))
    for (k, C) in enumerate(s.fast_classes)
        g = s.segment[C[1]]
        classes[g] += 1
        s.trap[k] || (lowest[g] += 1)
    end
    one_law = length(s.law_classes) == 1
    one_law && all(==(1), classes) && return _NO_FAULT
    if one_law && length(s.aggregated) == 1
        R = s.aggregated[1]
        !any(k -> s.trap[k], R) &&
            all(k -> lowest[s.segment[s.fast_classes[k][1]]] == 1, R) && return _NO_FAULT
    end
    !one_law ? _LAW_MULTI : length(s.aggregated) > 1 ? _MULTILEVEL :
        any(k -> s.trap[k], s.aggregated[1]) ? _TRAP : _SPLIT
end

"""The copy (1 active, 2 inactive) of node `u` of `net`."""
_copy_of(net::_Network, u::Int) = (u - 1) ÷ net.nforms + 1
"""The form of node `u` of `net`."""
_form_of(net::_Network, u::Int) = (u - 1) % net.nforms + 1

"""The closed class of face structure `s` that holds no active free E (node 1), or its
first when each holds none."""
_away_class(s) = s.closed[something(findfirst(C -> !insorted(1, C), s.closed), 1)]

"""The edges whose transitions leave the nodes `C` and take up a zeroed name on the face
of `s`, in transition order, each once."""
_exits(net::_Network, s, C::Vector{Int}) =
    unique!([net.edge[t] for t in eachindex(net.from)
             if insorted(net.from[t], C) && !insorted(net.to[t], C) && !s.alive[t] &&
                net.edge[t] != 0])

"""
The trap class 3 found on a face (stage `_TRAP`) and the rapid-equilibrium edge that feeds
it or joins it to its segment, as `(trap, edge, joined)`, the trap an index into
`s.fast_classes`. For an alive slow transition out of a class of the aggregated chain's
closed class, the edge is the last fast transition of an alive fast path from its target
into a trap of that class (`joined` false); when the slow transition lands in the trap, or
no transition feeds one (the closed class is a trap alone), it is the rapid-equilibrium
edge that joins the trap to the rest of its segment, whose way out takes up a zeroed name
(`joined` true).
"""
function _trap_feeder(net::_Network, s)
    nn = net.nforms * net.ncopies
    in_R = falses(length(s.fast_classes))
    in_R[s.aggregated[1]] .= true
    trap_of(u) = (q = s.class_of[u]; q > 0 && in_R[q] && s.trap[q] ? q : 0)
    function joining(q)
        for t in eachindex(net.from)
            net.fast[t] && s.class_of[net.from[t]] == q && s.class_of[net.to[t]] != q &&
                return (q, net.edge[t], true)
        end
        (q, 0, true)
    end
    fast_out = _out_transitions(net, s.alive .& net.fast)
    for t in eachindex(net.from)
        k = s.class_of[net.from[t]]
        s.alive[t] && !net.fast[t] && k > 0 && in_R[k] || continue
        v = net.to[t]
        landed = trap_of(v)
        landed > 0 && return joining(landed)
        seen = falses(nn)
        seen[v] = true
        todo = [v]
        while !isempty(todo)
            x = pop!(todo)
            for t2 in fast_out[x]
                y = net.to[t2]
                fed = trap_of(y)
                fed > 0 && s.class_of[x] != fed && return (fed, net.edge[t2], false)
                seen[y] || (seen[y] = true; push!(todo, y))
            end
        end
    end
    joining(s.aggregated[1][findfirst(k -> s.trap[k], s.aggregated[1])])
end

"""
The split class 3 found on a face (stage `_SPLIT`): the first rapid-equilibrium segment
that the aggregated chain's closed class meets and that holds more than one lowest fast
class, those classes (`parts`), and the segment's rapid-equilibrium edges whose
transitions take up a zeroed name (`joins`), the only edges that could join them.
"""
function _split_parts(net::_Network, s)
    lowest_classes = [Int[] for _ in s.segments]
    for (k, C) in enumerate(s.fast_classes)
        s.trap[k] || push!(lowest_classes[s.segment[C[1]]], k)
    end
    segment_of(k) = s.segment[s.fast_classes[k][1]]
    R = s.aggregated[1]
    first_split = findfirst(k -> length(lowest_classes[segment_of(k)]) != 1, R)
    segment = segment_of(R[something(first_split, 1)])
    joins = unique!([net.edge[t] for t in eachindex(net.from)
                     if net.fast[t] && !s.alive[t] && net.edge[t] != 0 &&
                        s.segment[net.from[t]] == segment])
    (; segment, parts = [s.fast_classes[k] for k in lowest_classes[segment]], joins)
end

"""The kinetic-group tags the limit test ties: none for a plain mechanism."""
_tags(::Mechanism) = Symbol[]
_tags(am::AllostericMechanism) = am.cat_allo_states

"""
    _check_face(mech, net, faces, i) -> _Verdict

Class 3 (spec section 8) on face `faces[i]` of `mech`, whose law
network is `net`: procedure A rejects a second closed class or a bottomless segment and
accepts a face whose structure fixes the limit (`_face_stage`); otherwise the limit test
(`_face_limit_check`, seeded by the content hash and the face) decides, and a rejection
carries the stage it stopped at, the edge that feeds the trap or joins it to its segment
(`_trap_feeder`) and the split's first joining edge. A face the limit test accepts reports
`_LIMIT_EQUAL`.
"""
function _check_face(mech::_AnyMech, net::_Network, faces::Vector{UInt64}, i::Int)
    s = _face_structure(net, faces[i])
    stage = _face_stage(s)
    stage === _NO_FAULT && return _Verdict(_NO_FAULT, _LIMIT_UNTESTED, i, 0, 0)
    if stage === _NONUNIQUE
        C = _away_class(s)
        exits = _exits(net, s, C)
        return _Verdict(_NONUNIQUE, _LIMIT_UNTESTED, i, _copy_of(net, C[1]),
                        isempty(exits) ? 0 : exits[1])
    end
    stage === _BOTTOMLESS &&
        return _Verdict(_BOTTOMLESS, _LIMIT_UNTESTED, i,
                        _copy_of(net, s.segments[s.bottomless][1]), 0)
    limit = _face_limit_check(mech, net, s, faces[i], hash(i, hash(mech)))
    limit === _LIMIT_EQUAL && return _Verdict(_NO_FAULT, _LIMIT_EQUAL, i, 0, 0)
    if stage === _TRAP
        trap, edge = _trap_feeder(net, s)
        return _Verdict(_TRAP, limit, i, _copy_of(net, s.fast_classes[trap][1]), edge)
    elseif stage === _SPLIT
        split = _split_parts(net, s)
        return _Verdict(_SPLIT, limit, i, _copy_of(net, s.segments[split.segment][1]),
                        isempty(split.joins) ? 0 : split.joins[1])
    end
    _Verdict(stage, limit, i, 0, 0)
end

"""The first class-3 fault of `mech` on its law network `net`, face by face in the order
of `_faces`; `_WELL_POSED` when every face passes."""
function _faces_verdict(mech::_AnyMech, net::_Network)
    faces = _faces(catalytic(mech))
    for i in eachindex(faces)
        v = _check_face(mech, net, faces, i)
        v.fault === _NO_FAULT || return v
    end
    _WELL_POSED
end

"""
Whether `mech`, its rapid-equilibrium flags overridden by `re`, passes class 2 and then
class 3 on every face: the flip move asks it of a candidate (`_minimal_flips`), which it
reads through the override, building nothing. Class 2 runs first because a candidate's
flags may close a fast cycle its parent's did not.
"""
function _faces_pass(mech::_AnyMech, re::BitVector)
    isempty(_re_turnover_transitions(_eps_network(mech, re))) || return false
    _faces_verdict(mech, _law_network(mech, re)) === _WELL_POSED
end

# ─── Class 3: the limit test mod 2^61 − 1 ───────────────────────────────

"""The prime 2^61 − 1: the limit test computes every rate, weight and limit modulo it."""
const _P61 = (Int64(1) << 61) - 1
"""A limit that grows without bound, outside the residues 0…p − 1."""
const _INFINITE = Int64(-1)
"""A law that has no value on a face, outside the residues 0…p − 1."""
const _UNDEFINED = Int64(-2)

_mulp(a::Int64, b::Int64) = Int64(mod(widemul(a, b), _P61))
_addp(a::Int64, b::Int64) = (s = a + b; s >= _P61 ? s - _P61 : s)
_subp(a::Int64, b::Int64) = (s = a - b; s < 0 ? s + _P61 : s)
_modp(x::Integer) = Int64(mod(x, _P61))

"""`a^k` mod p; a negative `k` raises the inverse."""
function _powp(a::Int64, k::Integer)
    k < 0 && return _powp(_invp(a), -k)
    r = Int64(1)
    while k > 0
        isodd(k) && (r = _mulp(r, a))
        a = _mulp(a, a)
        k >>= 1
    end
    r
end

"""The inverse of `a` mod p."""
_invp(a::Int64) = (a == 0 && throw(DivideError()); _powp(a, _P61 - 2))

"""A deterministic stream of residues in 2…p − 2 (splitmix64), seeded by a mechanism's
content hash and a face, so every process draws the same points."""
mutable struct _ResidueStream
    state::UInt64
end

"""The next residue of `r`."""
function _residue!(r::_ResidueStream)
    r.state += 0x9e3779b97f4a7c15
    z = r.state
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    Int64(2 + (z ⊻ (z >> 31)) % UInt64(_P61 - 3))
end

"""`(A \\ B, det(A))` mod p for square `A`, by Gauss–Jordan elimination: the determinant is
the product of the pivots, its sign flipped by each row swap; `A \\ B` is an empty matrix,
and the determinant 0, when `A` is singular."""
function _solvep(A::Matrix{Int64}, B::Matrix{Int64})
    n = size(A, 1)
    M = hcat(A, B)
    d = Int64(1)
    for j in 1:n
        k = findfirst(i -> M[i, j] != 0, j:n)
        k === nothing && return Matrix{Int64}(undef, 0, 0), Int64(0)
        k += j - 1
        if k != j
            M[[j, k], :] = M[[k, j], :]
            d = _subp(Int64(0), d)
        end
        d = _mulp(d, M[j, j])
        pivot = _invp(M[j, j])
        for c in j:size(M, 2)
            M[j, c] = _mulp(M[j, c], pivot)
        end
        for i in 1:n
            (i == j || M[i, j] == 0) && continue
            f = M[i, j]
            for c in j:size(M, 2)
                M[i, c] = _subp(M[i, c], _mulp(f, M[j, c]))
            end
        end
    end
    M[:, (n + 1):end], d
end

"""The stationary distribution mod p of the chain with off-diagonal rates `R[i, j]`
(`i → j`) and one closed class covering every state, and the determinant of the system that
gives it (the generator's transpose, its last row ones); empty, and 0, when singular."""
function _stationaryp(R::Matrix{Int64})
    n = size(R, 1)
    n == 1 && return Int64[1], Int64(1)
    A = zeros(Int64, n, n)
    for i in 1:n, j in 1:n
        i == j && continue
        A[j, i] = _addp(A[j, i], R[i, j])
        A[i, i] = _subp(A[i, i], R[i, j])
    end
    A[n, :] .= 1
    b = zeros(Int64, n, 1)
    b[n] = 1
    x, d = _solvep(A, b)
    vec(x), d
end

# Polynomials mod p: coefficient vectors, constant term first.
_ptrimp!(a::Vector{Int64}) = (while length(a) > 1 && a[end] == 0; pop!(a); end; a)
_degp(a::Vector{Int64}) =
    (b = _ptrimp!(copy(a)); length(b) == 1 && b[1] == 0 ? -1 : length(b) - 1)
_paddp(a::Vector{Int64}, b::Vector{Int64}) =
    _ptrimp!([_addp(get(a, i, Int64(0)), get(b, i, Int64(0)))
              for i in 1:max(length(a), length(b))])
_psubp(a::Vector{Int64}, b::Vector{Int64}) =
    _ptrimp!([_subp(get(a, i, Int64(0)), get(b, i, Int64(0)))
              for i in 1:max(length(a), length(b))])
function _pmulp(a::Vector{Int64}, b::Vector{Int64})
    c = zeros(Int64, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] = _addp(c[i + j - 1], _mulp(a[i], b[j]))
    end
    _ptrimp!(c)
end
function _pdivremp(a::Vector{Int64}, b::Vector{Int64})
    a, b = _ptrimp!(copy(a)), _ptrimp!(copy(b))
    db = _degp(b)
    db < 0 && throw(DivideError())
    q = zeros(Int64, max(_degp(a) - db + 1, 1))
    lead = _invp(b[end])
    while _degp(a) >= db
        k = _degp(a) - db
        c = _mulp(a[end], lead)
        q[k + 1] = c
        for i in 0:db
            a[i + k + 1] = _subp(a[i + k + 1], _mulp(c, b[i + 1]))
        end
        _ptrimp!(a)
    end
    _ptrimp!(q), a
end
_pevalp(a::Vector{Int64}, x::Int64) =
    foldr((c, acc) -> _addp(c, _mulp(acc, x)), a; init = Int64(0))
"""The lowest-order term of `a` as `(order, coefficient)`; `(typemax(Int), 0)` for 0."""
_plowp(a::Vector{Int64}) =
    (k = findfirst(!=(0), a); k === nothing ? (typemax(Int), Int64(0)) : (k - 1, a[k]))

"""Newton interpolation mod p through the points `(x, y)`."""
function _interpp(x::Vector{Int64}, y::Vector{Int64})
    n = length(x)
    c = copy(y)
    for j in 2:n, i in n:-1:j
        c[i] = _mulp(_subp(c[i], c[i - 1]), _invp(_subp(x[i], x[i - j + 1])))
    end
    p = Int64[c[n]]
    for i in (n - 1):-1:1
        p = _paddp(_pmulp(p, Int64[_subp(Int64(0), x[i]), 1]), Int64[c[i]])
    end
    p
end

"""`(N, D)` with degrees at most `dmax` and `N/D` through the points: the extended
Euclidean algorithm on the interpolant (it needs at least `2 dmax + 1` points)."""
function _ratrecp(x::Vector{Int64}, y::Vector{Int64}, dmax::Int)
    M = Int64[1]
    for xi in x
        M = _pmulp(M, Int64[_subp(Int64(0), xi), 1])
    end
    r0, r1 = M, _interpp(x, y)
    t0, t1 = Int64[0], Int64[1]
    while _degp(r1) > dmax
        q, r = _pdivremp(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, _psubp(t0, _pmulp(q, t1))
    end
    r1, t1
end

"""The rational function of one variable that `f` evaluates mod p, reconstructed from
random points with a degree bound doubled until two more points agree; empty vectors past
degree 1024."""
function _reconstructp(f, rng::_ResidueStream)
    d = 8
    while d <= 1024
        xs = [_residue!(rng) for _ in 1:(2d + 3)]
        ys = [f(x) for x in xs]
        N, D = _ratrecp(xs[1:(end - 2)], ys[1:(end - 2)], d)
        all(k -> (dv = _pevalp(D, xs[k]); dv != 0 &&
                  _mulp(_pevalp(N, xs[k]), _invp(dv)) == ys[k]),
            (length(xs) - 1):length(xs)) && return N, D
        d *= 2
    end
    Int64[], Int64[]
end

"""
A random point mod p of `m`'s model for `ncopies` conformations on
face `face`: log-free energies of every form in each conformation and of every
metabolite name, an integer null-space basis of the kinetic-group ties (members equal in
each conformation that keeps them, an `:EqualAI` member equal across) raised to random
residues; from them each edge's equilibrium constant `K[k, c]` and `L`; a steady-state
speed per kinetic group and conformation (one for both of an `:EqualAI` group), a
rapid-equilibrium speed per edge and conformation, the flip's speed; concentrations
`conc` per name, zero on the face, and `concw`, one on the face, which the weights'
leading coefficients read.
"""
function _model_point(m::Mechanism, tags::Vector{Symbol}, ncopies::Int, face::UInt64,
                      rng::_ResidueStream)
    nf, name_of = length(m.forms), _name_index(m)
    nnames = maximum(name_of; init = 0)
    nv = ncopies * nf + nnames
    excluded(k, c) = !isempty(tags) &&
        ((c == 2 && tags[m.group[k]] === :OnlyA) || (c == 1 && tags[m.group[k]] === :OnlyI))
    function drop!(row, k, c, dir)
        e = m.edges[k]
        row[(c - 1) * nf + e.from] += dir
        row[(c - 1) * nf + e.to] -= dir
        foreach(i -> row[ncopies * nf + name_of[i]] += dir, entries(e.consumed))
        foreach(i -> row[ncopies * nf + name_of[i]] -= dir, entries(e.released))
        row
    end
    ties = Vector{Int}[]
    for members in _group_members(m), c in 1:ncopies
        lead = first(members)
        excluded(lead, c) && continue
        for k in members[2:end]
            push!(ties, drop!(drop!(zeros(Int, nv), k, c, 1), lead, c, -1))
        end
        if c == 2 && tags[m.group[lead]] === :EqualAI
            for k in members
                push!(ties, drop!(drop!(zeros(Int, nv), k, 2, 1), k, 1, -1))
            end
        end
    end
    B = _integer_nullspace(isempty(ties) ? zeros(Int, 0, nv) :
                           permutedims(reduce(hcat, ties)))
    t = [_residue!(rng) for _ in axes(B, 2)]
    x = [reduce(_mulp, (_powp(t[j], B[i, j]) for j in axes(B, 2) if B[i, j] != 0);
                init = Int64(1)) for i in 1:nv]
    K = ones(Int64, length(m.edges), ncopies)
    for (k, e) in enumerate(m.edges), c in 1:ncopies
        num, den = x[(c - 1) * nf + e.from], x[(c - 1) * nf + e.to]
        for i in entries(e.consumed)
            num = _mulp(num, x[ncopies * nf + name_of[i]])
        end
        for i in entries(e.released)
            den = _mulp(den, x[ncopies * nf + name_of[i]])
        end
        K[k, c] = _mulp(num, _invp(den))
    end
    ss = [_residue!(rng) for _ in 1:ngroups(m), _ in 1:ncopies]
    for g in 1:ngroups(m)
        ncopies == 2 && tags[g] === :EqualAI && (ss[g, 2] = ss[g, 1])
    end
    re = [_residue!(rng) for _ in m.edges, _ in 1:ncopies]
    flip = _residue!(rng)
    L = ncopies == 2 ? _mulp(x[1], _invp(x[nf + 1])) : Int64(1)
    conc = [_zeroes(face, n) ? Int64(0) : _residue!(rng) for n in 1:nnames]
    concw = [_zeroes(face, n) ? Int64(1) : conc[n] for n in 1:nnames]
    (; K, ss, re, flip, L, conc, concw)
end

"""`pt` with new rapid-equilibrium and flip speeds: another speed pattern."""
_respeed_point(pt, rng::_ResidueStream) =
    merge(pt, (re = [_residue!(rng) for _ in axes(pt.re, 1), _ in axes(pt.re, 2)],
               flip = _residue!(rng)))

"""The rate mod p of every transition of `net` at the concentrations `conc` (per name):
its speed times `K` in the stored direction, its speed alone in reverse, times each
metabolite it takes up, a fast one times `Λ` unless `Λ` is 0."""
function _rates_p(m::Mechanism, net::_Network, pt, conc::Vector{Int64}, Λ::Int64)
    name_of = _name_index(m)
    r = zeros(Int64, length(net.from))
    for t in eachindex(net.from)
        k, c = net.edge[t], net.conformation[t]
        if k == 0
            r[t] = net.forward[t] ? _mulp(pt.flip, pt.L) : pt.flip
        else
            speed = net.fast[t] ? pt.re[k, c] : pt.ss[m.group[k], c]
            r[t] = net.forward[t] ? _mulp(speed, pt.K[k, c]) : speed
            e = m.edges[k]
            for i in entries(net.forward[t] ? e.consumed : e.released)
                r[t] = _mulp(r[t], conc[name_of[i]])
            end
        end
        net.fast[t] && Λ != 0 && (r[t] = _mulp(r[t], Λ))
    end
    r
end

"""ρ(u) mod p: over the slow transitions out of node `u`, their rate `r` times their
uptake of `x1` plus the fall of its potential `h` along them."""
function _slow_rho(net::_Network, r::Vector{Int64}, h::Vector{Int}, x1::Int)
    rho = zeros(Int64, net.nforms * net.ncopies)
    for t in eachindex(net.from)
        (net.fast[t] || r[t] == 0) && continue
        u = net.from[t]
        rho[u] = _addp(rho[u],
                       _mulp(r[t], _modp(net.change[x1, t] + h[u] - h[net.to[t]])))
    end
    rho
end

"""Equilibrium weights mod p within each part of `parts` over the fast transitions inside
it, from its first node, at rates `r`."""
function _eq_weights_p(net::_Network, parts::Vector{Vector{Int}}, r::Vector{Int64})
    nn = net.nforms * net.ncopies
    part = zeros(Int, nn)
    for (q, P) in enumerate(parts), u in P
        part[u] = q
    end
    next = [Int[] for _ in 1:nn]
    for t in eachindex(net.from)
        u = net.from[t]
        net.fast[t] && part[u] > 0 && part[u] == part[net.to[t]] && push!(next[u], t)
    end
    w = zeros(Int64, nn)
    for P in parts
        w[P[1]] = 1
        done = falses(nn)
        done[P[1]] = true
        todo = [P[1]]
        while !isempty(todo)
            u = pop!(todo)
            for t in next[u]
                v = net.to[t]
                done[v] && continue
                back = isodd(t) ? t + 1 : t - 1
                r[t] != 0 && r[back] != 0 ||
                    error("_face_limit_check: a fast transition inside a part runs one way")
                w[v] = _mulp(w[u], _mulp(r[t], _invp(r[back])))
                done[v] = true
                push!(todo, v)
            end
        end
        all(u -> done[u], P) || error("_face_limit_check: a part is not fast-connected")
    end
    w
end

"""The face chain's limit mod p as every rapid-equilibrium speed grows, by first-order
aggregation over the fast classes of the aggregated chain's one closed class."""
function _limit_by_aggregation(m::Mechanism, net::_Network, s, pt, h::Vector{Int},
                               x1::Int)
    nn = net.nforms * net.ncopies
    r = _rates_p(m, net, pt, pt.conc, Int64(0))
    R = s.aggregated[1]
    place = zeros(Int, length(s.fast_classes))
    place[R] .= eachindex(R)
    nu = _eq_weights_p(net, [s.fast_classes[k] for k in R],
                       _rates_p(m, net, pt, pt.concw, Int64(0)))
    for k in R
        C = s.fast_classes[k]
        inverse = _invp(reduce(_addp, nu[C]; init = Int64(0)))
        nu[C] .= _mulp.(nu[C], inverse)
    end
    transient = [u for u in 1:nn if net.present[u] && s.class_of[u] == 0]
    slot = zeros(Int, nn)
    slot[transient] .= eachindex(transient)
    X = zeros(Int64, length(transient), length(R))
    if !isempty(transient)
        M = zeros(Int64, length(transient), length(transient))
        rhs = zeros(Int64, length(transient), length(R))
        for t in eachindex(net.from)
            u, v = net.from[t], net.to[t]
            net.fast[t] && s.alive[t] && r[t] != 0 && slot[u] > 0 || continue
            M[slot[u], slot[u]] = _addp(M[slot[u], slot[u]], r[t])
            if slot[v] > 0
                M[slot[u], slot[v]] = _subp(M[slot[u], slot[v]], r[t])
            elseif s.class_of[v] > 0 && place[s.class_of[v]] > 0
                q = place[s.class_of[v]]
                rhs[slot[u], q] = _addp(rhs[slot[u], q], r[t])
            end
        end
        X = first(_solvep(M, rhs))
        isempty(X) && error("_face_limit_check: singular absorption system")
    end
    Q = zeros(Int64, length(R), length(R))
    for t in eachindex(net.from)
        u, v = net.from[t], net.to[t]
        !net.fast[t] && s.alive[t] && r[t] != 0 && s.class_of[u] > 0 &&
            place[s.class_of[u]] > 0 || continue
        a, f = place[s.class_of[u]], _mulp(nu[u], r[t])
        if s.class_of[v] > 0
            b = place[s.class_of[v]]
            b == 0 && error("_face_limit_check: the aggregated class is not closed")
            b != a && (Q[a, b] = _addp(Q[a, b], f))
        else
            for b in eachindex(R)
                b == a || (Q[a, b] = _addp(Q[a, b], _mulp(f, X[slot[v], b])))
            end
        end
    end
    mu = first(_stationaryp(Q))
    isempty(mu) && error("_face_limit_check: singular aggregated chain")
    rho = _slow_rho(net, r, h, x1)
    v = Int64(0)
    for (q, k) in enumerate(R), u in s.fast_classes[k]
        v = _addp(v, _mulp(_mulp(mu[q], nu[u]), rho[u]))
    end
    v
end

"""The face chain's limit mod p as every rapid-equilibrium speed grows, by interpolating
the rate's numerator and its tree-sum denominator in Λ (degree below the closed class's
size) on the face chain's closed class; `_INFINITE` when the rate grows with Λ."""
function _limit_by_interpolation(m::Mechanism, net::_Network, s, pt, h::Vector{Int},
                                 x1::Int, rng::_ResidueStream)
    C = s.closed[1]
    n = length(C)
    pos = zeros(Int, net.nforms * net.ncopies)
    pos[C] .= 1:n
    xs, Ds, Ns = Int64[], Int64[], Int64[]
    for _ in 1:(n + 1)
        Λ = _residue!(rng)
        r = _rates_p(m, net, pt, pt.conc, Λ)
        Q = zeros(Int64, n, n)
        for t in eachindex(net.from)
            a, b = pos[net.from[t]], pos[net.to[t]]
            s.alive[t] && r[t] != 0 && a > 0 && b > 0 && a != b || continue
            Q[a, b] = _addp(Q[a, b], r[t])
        end
        x, D = _stationaryp(Q)
        isempty(x) && error("_face_limit_check: singular at a random Λ")
        rho = _slow_rho(net, r, h, x1)
        v = reduce(_addp, (_mulp(x[pos[u]], rho[u]) for u in C); init = Int64(0))
        push!(xs, Λ)
        push!(Ds, D)
        push!(Ns, _mulp(D, v))
    end
    cD, cN = _interpp(xs[1:n], Ds[1:n]), _interpp(xs[1:n], Ns[1:n])
    _pevalp(cD, xs[end]) == Ds[end] && _pevalp(cN, xs[end]) == Ns[end] ||
        error("_face_limit_check: the degree bound failed")
    dD, dN = _degp(cD), _degp(cN)
    dN > dD && return _INFINITE
    dN < dD && return Int64(0)
    _mulp(cN[dD + 1], _invp(cD[dD + 1]))
end

"""The law's value mod p on the face, from the law's segment chain with its one closed
class: each segment's weight on its lowest forms, its alive slow transitions out of
them."""
function _law_on_face(m::Mechanism, net::_Network, s, pt, h::Vector{Int}, x1::Int)
    r = _rates_p(m, net, pt, pt.conc, Int64(0))
    parts = [sort(S; by = u -> !s.lowest[u]) for S in s.segments]
    w = _eq_weights_p(net, parts, _rates_p(m, net, pt, pt.concw, Int64(0)))
    share = zeros(Int64, length(w))
    for S in s.segments
        low = [u for u in S if s.lowest[u]]
        inverse = _invp(reduce(_addp, w[low]; init = Int64(0)))
        share[low] .= _mulp.(w[low], inverse)
    end
    Lc = s.law_classes[1]
    at = zeros(Int, length(s.segments))
    at[Lc] .= eachindex(Lc)
    Q = zeros(Int64, length(Lc), length(Lc))
    for t in eachindex(net.from)
        u = net.from[t]
        sa, sb = s.segment[u], s.segment[net.to[t]]
        !net.fast[t] && s.alive[t] && s.lowest[u] && at[sa] > 0 && sa != sb || continue
        at[sb] == 0 && error("_face_limit_check: the law's closed class is not closed")
        Q[at[sa], at[sb]] = _addp(Q[at[sa], at[sb]], _mulp(share[u], r[t]))
    end
    mu = first(_stationaryp(Q))
    isempty(mu) && error("_face_limit_check: singular law segment chain")
    rho = _slow_rho(net, r, h, x1)
    v = Int64(0)
    for (q, g) in enumerate(Lc), u in s.segments[g]
        s.lowest[u] && (v = _addp(v, _mulp(_mulp(mu[q], share[u]), rho[u])))
    end
    v
end

"""The law mod p at concentrations `conc` with every name present: segment weights and
the slow segment chain over every segment."""
function _law_inside(m::Mechanism, net::_Network, s, pt, conc::Vector{Int64},
                     h::Vector{Int}, x1::Int)
    r = _rates_p(m, net, pt, conc, Int64(0))
    w = _eq_weights_p(net, s.segments, r)
    W = [_invp(reduce(_addp, w[S]; init = Int64(0))) for S in s.segments]
    Q = zeros(Int64, length(s.segments), length(s.segments))
    for t in eachindex(net.from)
        net.fast[t] && continue
        a, b = s.segment[net.from[t]], s.segment[net.to[t]]
        a == b && continue
        Q[a, b] = _addp(Q[a, b], _mulp(_mulp(w[net.from[t]], W[a]), r[t]))
    end
    mu = first(_stationaryp(Q))
    isempty(mu) && error("_face_limit_check: singular segment chain")
    rho = _slow_rho(net, r, h, x1)
    v = Int64(0)
    for (g, S) in enumerate(s.segments), u in S
        v = _addp(v, _mulp(_mulp(mu[g], _mulp(w[u], W[g])), rho[u]))
    end
    v
end

"""The law's value mod p on the face when its segment chain has several closed classes:
its limit as the zeroed concentrations go to zero along two random directions, by
rational reconstruction in their scale; `_UNDEFINED` when a direction meets a pole or the
two differ."""
function _law_by_reconstruction(m::Mechanism, net::_Network, s, pt, h::Vector{Int},
                                x1::Int, face::UInt64, rng::_ResidueStream)
    found = Int64[]
    for _ in 1:2
        direction = [_residue!(rng) for _ in pt.conc]
        law_at(δ) = _law_inside(m, net, s, pt,
                                [_zeroes(face, n) ? _mulp(direction[n], δ) : pt.conc[n]
                                 for n in eachindex(pt.conc)], h, x1)
        num, den = _reconstructp(law_at, rng)
        isempty(den) && error("_face_limit_check: no reconstruction under the degree cap")
        (a, lead_num), (b, lead_den) = _plowp(num), _plowp(den)
        push!(found, a > b ? Int64(0) : a < b ? _INFINITE :
                     _mulp(lead_num, _invp(lead_den)))
    end
    found[1] == _INFINITE || found[1] != found[2] ? _UNDEFINED : found[1]
end

"""
    _face_limit_check(mech, net, s, face, seed) -> _LimitVerdict

The limit test of class 3 on face `face` of `mech`, whose law
network is `net` and face structure `s`: at two random points of the model mod 2^61 − 1,
drawn from `seed`, the fast-RE limit of the face chain (by aggregation when the aggregated
chain has one closed class, else by interpolation in Λ) against the law's face value
(from the law's segment chain when it has one closed class, else by reconstruction as the
zeroed concentrations vanish). On a difference a third limit at new speeds tells
`_LIMIT_DIFFER` from `_LIMIT_SPEED`. A false `_LIMIT_EQUAL` has probability at most twice
the degree over p.
"""
function _face_limit_check(mech::_AnyMech, net::_Network, s, face::UInt64, seed::UInt64)
    m, tags = catalytic(mech), _tags(mech)
    rng = _ResidueStream(seed)
    x1 = _turnover_name(m)
    # Class 2, which runs first, leaves no fast cycle that takes `x1` up: a potential.
    h = s.offset[x1, :]
    limit(pt) = length(s.aggregated) == 1 ? _limit_by_aggregation(m, net, s, pt, h, x1) :
                _limit_by_interpolation(m, net, s, pt, h, x1, rng)
    for _ in 1:2
        pt = _model_point(m, tags, net.ncopies, face, rng)
        law = length(s.law_classes) == 1 ? _law_on_face(m, net, s, pt, h, x1) :
              _law_by_reconstruction(m, net, s, pt, h, x1, face, rng)
        law == _UNDEFINED && return _LIMIT_LAW_NAN
        v = limit(pt)
        # The face chain's fast transitions are the law network's, a subset of the
        # ε-network's, where class 2, which runs first, finds no turnover cycle.
        v == _INFINITE && error("_face_limit_check: the limit grows without bound")
        v == law && continue
        return limit(_respeed_point(pt, rng)) == v ? _LIMIT_DIFFER : _LIMIT_SPEED
    end
    _LIMIT_EQUAL
end

# ─── Admission: the interior (spec section 7) ──────────────────────────

"""
Each form's component under the edges of `am`'s catalytic graph that conformation `c`
(1 active, 2 inactive) keeps (`_excludes`), every form included, labelled by its smallest
form index. Admission reads each plain binding `c` excludes as an arc from its unbound
side's component to its bound side's (`_arcs`).
"""
function _exclusion_components(am::AllostericMechanism, c::Int)
    m = am.catalytic
    parent = collect(eachindex(m.forms))
    top(i) = (while parent[i] != i; i = parent[i]; end; i)
    for (k, e) in enumerate(m.edges)
        _excludes(am, k, c) && continue
        a, b = top(e.from), top(e.to)
        a == b || (parent[max(a, b)] = min(a, b))
    end
    [top(f) for f in eachindex(m.forms)]
end

"""The plain bindings of `am`'s catalytic graph that conformation `c` excludes, in edge
order: admission's arcs, each stored as its uptake, unbound side first."""
_arcs(am::AllostericMechanism, c::Int) =
    [k for (k, e) in enumerate(am.catalytic.edges)
     if _excludes(am, k, c) && _is_plain_binding(am.catalytic, e)]

"""The first edge of `am`'s catalytic graph that conformation `c` excludes; 0 when it
excludes none."""
function _first_excluded(am::AllostericMechanism, c::Int)
    for k in eachindex(am.catalytic.edges)
        _excludes(am, k, c) && return k
    end
    0
end

"""
The arcs of one directed cycle that `arcs` close over the components `comp` (each arc an
edge of `m` read from its unbound side's component to its bound side's), in cycle order;
empty when they close none. Kahn's algorithm removes every component no cycle feeds; a
walk back along the remaining arcs then lands on a cycle.
"""
function _arc_cycle(m::Mechanism, comp::Vector{Int}, arcs::Vector{Int})
    tail(k) = comp[m.edges[k].from]
    head(k) = comp[m.edges[k].to]
    indeg = zeros(Int, length(comp))
    for k in arcs
        indeg[head(k)] += 1
    end
    todo = [r for r in unique(comp) if indeg[r] == 0]
    while !isempty(todo)
        u = pop!(todo)
        for k in arcs
            tail(k) == u && (indeg[head(k)] -= 1) == 0 && push!(todo, head(k))
        end
    end
    start = 0
    for k in arcs
        indeg[tail(k)] > 0 && indeg[head(k)] > 0 && (start = k; break)
    end
    start == 0 && return Int[]
    function back(x)
        for k in arcs
            head(k) == x && indeg[tail(k)] > 0 && return k
        end
        0
    end
    x = head(start)
    for _ in eachindex(comp)
        x = tail(back(x))
    end
    cycle = [back(x)]
    while tail(cycle[end]) != x
        push!(cycle, back(tail(cycle[end])))
    end
    reverse!(cycle)
end

"""The components (`comp` labels) that free E's reaches along `arcs`, edges of `m` read
from their unbound side's component to their bound side's."""
function _reached_components(m::Mechanism, comp::Vector{Int}, arcs::Vector{Int})
    next = [Int[] for _ in comp]
    for k in arcs
        push!(next[comp[m.edges[k].from]], comp[m.edges[k].to])
    end
    _reached(next, [comp[1]])
end

"""
Whether edge `e` of `m` runs, as stored, toward products: from the side that holds more
substrates, bound or free, then fewer products, the order canonical form gives a step that
exchanges several metabolites (`_orientation`). A step that takes up or gives off one
metabolite is stored as its uptake, so a half-reaction that releases a product,
`E(A) <--> E(; residual = A - P) + P`, runs toward products against its stored direction.
A step whose two sides tie, such as an isomerization of free forms, runs as stored.
"""
function _runs_toward_products(m::Mechanism, e::Edge)
    substrate_mask, product_mask = _role_masks(m)
    _side_score(m.forms, e.from, e.consumed, substrate_mask, product_mask) >=
        _side_score(m.forms, e.to, e.released, substrate_mask, product_mask)
end

"""
The edge an island's message names, 0 when free E's component reaches every component
(`reached`): of the edges conformation `c` excludes that join a reached component to an
unreached one, the first chemistry edge (`_is_chemistry`) that, run toward products
(`_runs_toward_products`), enters the unreached one; else the first chemistry edge; else
the first edge.
"""
function _island_edge(am::AllostericMechanism, c::Int, comp::Vector{Int},
                      reached::BitVector)
    m = am.catalytic
    best, rank = 0, 4
    for (k, e) in enumerate(m.edges)
        _excludes(am, k, c) || continue
        a, b = reached[comp[e.from]], reached[comp[e.to]]
        a == b && continue
        r = !_is_chemistry(m, e) ? 3 : a == _runs_toward_products(m, e) ? 1 : 2
        if r < rank
            best, rank = k, r
        end
    end
    best
end

"""
    _admission_interior(am) -> _Verdict

Admission's interior checks (spec section 7), for each conformation `c` whose tags exclude
an edge (`_excludes`), on the components of the forms its other edges join
(`_exclusion_components`) with its tagged plain bindings as arcs (`_arcs`):
- existence: no arc stays inside one component (`_NO_LIMIT`) and the arcs close no
  directed cycle (`_ARC_CYCLE`); equivalent to the exact cone test, since with every other
  constant fixed a direction's equilibrium-ratio exponents are differences of component
  potentials and only the arcs need a sign;
- uniqueness: free E's component reaches every component along the arcs, so every island
  lies behind a tagged plain binding on its bound side and the limit is deletion plus
  pruning (`_ISLAND`);
- turnover: when tags exclude edges of the active conformation, its law network still
  holds a cycle that turns the reaction over (`_NO_ACTIVE_TURNOVER`).
The verdict carries the conformation and the edge its message leads with (face 0).
"""
function _admission_interior(am::AllostericMechanism)
    m = am.catalytic
    for c in 1:2
        _first_excluded(am, c) == 0 && continue
        comp = _exclusion_components(am, c)
        arcs = _arcs(am, c)
        for k in arcs
            comp[m.edges[k].from] == comp[m.edges[k].to] &&
                return _Verdict(_NO_LIMIT, _LIMIT_UNTESTED, 0, c, k)
        end
        cycle = _arc_cycle(m, comp, arcs)
        isempty(cycle) || return _Verdict(_ARC_CYCLE, _LIMIT_UNTESTED, 0, c, cycle[1])
        k = _island_edge(am, c, comp, _reached_components(m, comp, arcs))
        k == 0 || return _Verdict(_ISLAND, _LIMIT_UNTESTED, 0, c, k)
    end
    k = _first_excluded(am, 1)
    k != 0 && !_has_turnover_cycle(m, _active_mask(am)) &&
        return _Verdict(_NO_ACTIVE_TURNOVER, _LIMIT_UNTESTED, 0, 1, k)
    _WELL_POSED
end

"""The edges of a shortest path from form `a` to form `b` over the edges of `am`'s
catalytic graph that conformation `c` keeps, in path order; empty when `a == b`. Callers
pass two forms those edges join."""
function _kept_path(am::AllostericMechanism, c::Int, a::Int, b::Int)
    m = am.catalytic
    via = zeros(Int, length(m.forms))
    seen = falses(length(m.forms))
    seen[a] = true
    queue = [a]
    head = 0
    while !seen[b]
        head += 1
        u = queue[head]
        for (k, e) in enumerate(m.edges)
            _excludes(am, k, c) && continue
            v = e.from == u ? e.to : e.to == u ? e.from : 0
            (v == 0 || seen[v]) && continue
            seen[v] = true
            via[v] = k
            push!(queue, v)
        end
    end
    path = Int[]
    f = b
    while f != a
        k = via[f]
        pushfirst!(path, k)
        f = m.edges[k].from == f ? m.edges[k].to : m.edges[k].from
    end
    path
end

"""
Whether `am` with the kinetic group of each edge in `ks` tagged `tag` passes the guard
(`_guard`). Every fix that retags steps in admission's messages is checked this way.
"""
function _retag_passes(am::AllostericMechanism, ks, tag::Symbol)
    states = copy(cat_allo_states(am))
    states[am.catalytic.group[ks]] .= tag
    _guard(_with(am; states)) === _WELL_POSED
end

"""
The message of admission's verdict `v` on `am`: the cause, naming each step by its DSL
line (`_step_text`), and each fix only when its variant passes the guard
(`_retag_passes`). A fix that tags steps is offered for the inactive conformation only,
since a step tagged `:OnlyI` must bind modifiers alone.
"""
function _admission_message(am::AllostericMechanism, v::_Verdict)
    m = am.catalytic
    c, k = v.conformation, v.edge
    conf = c == 1 ? "active" : "inactive"
    tag = c == 1 ? ":OnlyI" : ":OnlyA"
    line(j) = _step_text(m, j)
    lines(js) = join((_step_text(m, j) for j in js), "; ")
    untag_passes(ks) = all(t -> _retag_passes(am, ks, t), (:EqualAI, :NonequalAI))
    # The fixes whose variants pass: the groups of `untag` (`what`) tagged :EqualAI or
    # :NonequalAI, and the steps `joins` that join `sides` tagged :OnlyA.
    function fixes(untag::Vector{Int}, what::String, joins::Vector{Int}, sides::String)
        found = String[]
        untag_passes(untag) && push!(found, "tag $what :EqualAI or :NonequalAI")
        _retag_passes(am, joins, :OnlyA) &&
            push!(found, "tag :OnlyA the steps that join $sides: $(lines(joins))")
        isempty(found) ? "" : " " * uppercasefirst(join(found, ", or ")) * "."
    end
    if v.fault === _NO_LIMIT
        msg = "AllostericMechanism: in the $conf conformation, steps outside $tag " *
              "groups still join the two sides of the $tag binding $(line(k)), so " *
              "thermodynamics keeps its affinity finite and the binding has no limit " *
              "there."
        c == 2 || return msg
        e = m.edges[k]
        return msg * fixes([k], "its group", _kept_path(am, c, e.from, e.to), "its sides")
    elseif v.fault === _ARC_CYCLE
        cycle = _arc_cycle(m, _exclusion_components(am, c), _arcs(am, c))
        msg = "AllostericMechanism: the $tag bindings $(lines(cycle)) close a cycle in " *
              "the $conf conformation: steps outside $tag groups join each one's bound " *
              "side to the next one's unbound side, so thermodynamics forbids their " *
              "affinities to vanish together."
        c == 2 || return msg
        joins = Int[]
        for (i, j) in enumerate(cycle)
            next = cycle[mod1(i + 1, length(cycle))]
            append!(joins, _kept_path(am, c, m.edges[j].to, m.edges[next].from))
        end
        return msg * fixes(cycle, "them", unique!(joins), "them")
    elseif v.fault === _ISLAND
        comp = _exclusion_components(am, c)
        reached = _reached_components(m, comp, _arcs(am, c))
        e = m.edges[k]
        lost = reached[comp[e.from]] ? comp[e.to] : comp[e.from]
        fixable = c == 2 && untag_passes([k])
        msg = "AllostericMechanism: in the $conf conformation the forms " *
              join([_form_text(m, f) for f in eachindex(m.forms) if comp[f] == lost],
                   ", ") *
              " join free E only through $tag steps that are not bindings toward them, " *
              "such as $(line(k)), so their weight depends on how fast those steps vanish."
        msg = c == 2 && fixable ?
              msg * " Tag the group of $(line(k)) :EqualAI or :NonequalAI." : msg
        # A ping-pong whose inactive conformation runs neither half-reaction strands the
        # covalent forms between the halves; the roadmap lists what would hold them.
        neither = c == 2 && all(j -> !_is_chemistry(m, m.edges[j]) || _excludes(am, j, 2),
                                eachindex(m.edges))
        covalent = any(f -> comp[f] == lost && _has_residual(m.forms[f]),
                       eachindex(m.forms))
        return neither && covalent ?
               msg * " An inactive conformation that runs neither half-reaction of a " *
                     "ping-pong is not supported yet (see the Roadmap page)." : msg
    elseif v.fault in (_EPS_NONUNIQUE, _PUMPED, _BRIDGE)
        eps, law = _eps_network(am), _law_network(am)
        faces = _faces(m)
        face = faces[v.face]
        phrase = _face_phrase(m, face)
        s = _face_sets(eps, law.present, face)
        if v.fault === _EPS_NONUNIQUE
            outside = _class_outside_law(_closed_classes(eps, s.alive), law.present)
            forms = [_form_text(m, _form_of(eps, u)) *
                     (_copy_of(eps, u) == 1 ? " (active)" : " (inactive)") for u in outside]
            return "AllostericMechanism: $(phrase)the forms $(join(forms, ", ")) " *
                   "hold the enzyme they receive and are reached only through tagged " *
                   "steps, so the steady state depends on where the enzyme starts."
        elseif v.fault === _PUMPED
            basin = [u for u in findall(s.stranded) if _copy_of(eps, u) == c]
            exits = unique!([eps.edge[t] for t in eachindex(eps.from)
                             if s.alive[t] && eps.vanishing[t] && eps.from[t] in basin &&
                                !s.stranded[eps.to[t]]])
            out = isempty(exits) ?
                "which no step running there leaves, so the enzyme collects in them" :
                "which leave only through the vanishing steps $(lines(exits)) while " *
                "nothing runs back the way it came, so their weight depends on how fast " *
                "those steps vanish"
            e = m.edges[k]
            inside = s.stranded[_node(eps, e.from, c)] && s.stranded[_node(eps, e.to, c)]
            return "AllostericMechanism: $(phrase)$(line(k)) " *
                   (inside ? "runs one way among" : "feeds") * " the forms " *
                   "$(join([_form_text(m, _form_of(eps, u)) for u in basin], ", ")) of " *
                   "the $conf conformation, $out."
        elseif c == 0
            return "AllostericMechanism: $(phrase)the tagged steps of both " *
                   "conformations together bridge forms that each conformation's other " *
                   "steps keep apart, so the limit depends on how fast those steps vanish."
        end
        nodes = s.fed .& .!s.stranded
        (; van) = _copy_transitions(eps, s.alive, nodes, c)
        ks = _bridging_edges(eps, nodes, s.n0, van)
        bridge = length(ks) == 1 ? "step $(line(k)) bridges" : "steps $(lines(ks)) bridge"
        msg = "AllostericMechanism: $(phrase)the $tag $bridge forms of the $conf " *
              "conformation that its other steps keep apart, so the limit depends on how " *
              "fast the tagged constants vanish."
        # Offer the steady-state fix only where its variant passes the guard: R2funnel
        # keeps its bridge with its tagged S binding steady-state (R2funnel_SS).
        re = _flags(m)
        re[_group_members(m)[m.group[k]]] .= false
        return _guard(_build(_raw(am, re))) === _WELL_POSED ?
               msg * " Make the group of $(line(k)) steady-state." : msg
    elseif v.fault === _NO_ACTIVE_TURNOVER
        ks = [j for j in eachindex(m.edges) if _excludes(am, j, 1)]
        msg = "AllostericMechanism: without its :OnlyI steps $(lines(ks)) the active " *
              "conformation turns no substrate into product, so the enzyme never turns " *
              "the reaction over."
        return untag_passes(ks) ? msg * " Tag their groups :EqualAI or :NonequalAI." : msg
    else
        error("_admission_message: $(v.fault) is no admission fault")
    end
end

# ─── Admission: the face condition ─────────────────────────────────────

"""
    _face_sets(net, law_nodes, face) -> (; alive, fed, n0, stranded)

What admission's face condition reads on face `face` of the ε-network `net`, whose law
network holds the nodes `law_nodes`: the transitions alive there (taking up no zeroed
name), the nodes free E's active copy reaches along them, the alive transitions that do
not vanish, and the stranded basin: fed nodes outside the law network that reach it along
no chain of alive transitions that do not vanish.
"""
function _face_sets(net::_Network, law_nodes::BitVector, face::UInt64)
    alive = BitVector([(net.takes[t] & face) == 0 for t in eachindex(net.takes)])
    n0 = alive .& .!net.vanishing
    fed = _reached(_successors(net, alive), [_node(net, 1, 1)])
    drains = _reached(_successors(net, n0; backward = true),
                      findall(law_nodes .& net.present))
    (; alive, fed, n0, stranded = fed .& .!law_nodes .& .!drains)
end

"""Copy `c`'s alive transitions with both ends in `nodes`: `van`, those that vanish with
their tags, and `held`, those T0+ reads (not vanishing, or fast)."""
function _copy_transitions(net::_Network, alive::BitVector, nodes::BitVector, c::Int)
    van, held = falses(length(net.from)), falses(length(net.from))
    for t in eachindex(net.from)
        alive[t] && net.conformation[t] == c && nodes[net.from[t]] && nodes[net.to[t]] ||
            continue
        van[t] = net.vanishing[t]
        held[t] = !net.vanishing[t] || net.fast[t]
    end
    (; van, held)
end

"""
The cost of a cheapest spanning out-arborescence of the nodes `1:n` rooted at `root` over
`arcs`, each `(from, to, cost)`, by Chu–Liu/Edmonds: take each node's cheapest entering
arc, contract the cycles they close and repeat; -1 when some node has no entering arc.
"""
function _min_arborescence(n::Int, arcs::Vector{NTuple{3, Int}}, root::Int)
    INF = typemax(Int) ÷ 4
    cost = 0
    while true
        inw, pre = fill(INF, n), zeros(Int, n)
        for (u, v, w) in arcs
            u != v && w < inw[v] && (inw[v] = w; pre[v] = u)
        end
        for v in 1:n
            v != root && inw[v] == INF && return -1
        end
        inw[root] = 0
        id, vis = zeros(Int, n), zeros(Int, n)
        cnt = 0
        for v in 1:n
            cost += inw[v]
            x = v
            while vis[x] != v && id[x] == 0 && x != root
                vis[x] = v
                x = pre[x]
            end
            if x != root && id[x] == 0
                cnt += 1
                y = pre[x]
                while y != x
                    id[y] = cnt
                    y = pre[y]
                end
                id[x] = cnt
            end
        end
        cnt == 0 && return cost
        for v in 1:n
            id[v] == 0 && (cnt += 1; id[v] = cnt)
        end
        next = NTuple{3, Int}[]
        for (u, v, w) in arcs
            id[u] != id[v] && push!(next, (id[u], id[v], w - inw[v]))
        end
        arcs, n, root = next, cnt, id[root]
    end
end

"""
    _fewest_slow(net, nodes, keep) -> Int

MINSLOW: the fewest slow transitions in a spanning in-tree of the nodes `nodes` over the
transitions `keep` flags (cost 1 slow, 0 fast), any root; -1 when no in-tree spans them. A
virtual root joined to every node at a cost above any tree's turns the in-tree into
Chu–Liu/Edmonds's out-arborescence (`_min_arborescence`) on reversed transitions; a cost of
two such joins means no single root serves.
"""
function _fewest_slow(net::_Network, nodes::BitVector, keep::BitVector)
    at = zeros(Int, length(nodes))
    n = 0
    for u in eachindex(nodes)
        nodes[u] && (n += 1; at[u] = n)
    end
    arcs = NTuple{3, Int}[]
    for t in eachindex(net.from)
        a, b = at[net.from[t]], at[net.to[t]]
        keep[t] && a != 0 && b != 0 && push!(arcs, (b, a, net.fast[t] ? 0 : 1))
    end
    for u in 1:n
        push!(arcs, (n + 1, u, n + 1))
    end
    c = _min_arborescence(n + 1, arcs, n + 1)
    c < 0 || c >= 2 * (n + 1) ? -1 : c - (n + 1)
end

"""
    _potential_on_sccs(net, keep, x) -> Bool

T0+: whether the uptake of name `x` (a `change` row, `_turnover_name`) is a potential on
every strongly connected part (`_strong_components`) of the transitions `keep` flags, so
no cycle of them turns the reaction over.
"""
function _potential_on_sccs(net::_Network, keep::BitVector, x::Int)
    comp = _strong_components(_successors(net, keep), net.present)
    inside = BitVector([keep[t] && comp[net.from[t]] == comp[net.to[t]]
                        for t in eachindex(net.from)])
    n = length(net.present)
    !_has_unbalanced_cycle(net.from, net.to, view(net.change, x, :), inside, zeros(Int, n),
                           zeros(Int, n))
end

"""Whether the transitions `extra` shorten the cheapest spanning in-tree of `nodes` over
`base` (`_fewest_slow`): R1's failure for one conformation's vanishing transitions, joint
R1's for both conformations'."""
function _shortens_in_tree(net::_Network, nodes::BitVector, base::BitVector,
                           extra::BitVector)
    without = _fewest_slow(net, nodes, base)
    with = _fewest_slow(net, nodes, base .| extra)
    with >= 0 && (without < 0 || with < without)
end

"""
The edges of the vanishing transitions on a cheapest spanning in-tree of `nodes` over
`base` and `van` (`_BRIDGE`), in transition order: an inclusion-minimal set of the
transitions `van` with which `base` still gives that tree's cost (`_fewest_slow`), so
every such tree uses each of them.
"""
function _bridging_edges(net::_Network, nodes::BitVector, base::BitVector, van::BitVector)
    best = _fewest_slow(net, nodes, base .| van)
    kept = copy(van)
    for t in eachindex(van)
        kept[t] || continue
        kept[t] = false
        _fewest_slow(net, nodes, base .| kept) == best || (kept[t] = true)
    end
    unique!([net.edge[t] for t in findall(kept)])
end

"""
Of two or more closed classes of the ε-network, the first that holds no node `law_nodes`
marks, else the second: the class an `_EPS_NONUNIQUE` verdict names. Whenever class 3 has
accepted the face such a class exists, because a closed class of the ε-network that holds
a law node holds a closed class of the law network, and class 3 leaves that network only
one. The fallback arises only from `_admission` on a face that class 3 rejects.
"""
function _class_outside_law(closed::Vector{Vector{Int}}, law_nodes::BitVector)
    for cls in closed
        any(view(law_nodes, cls)) || return cls
    end
    closed[2]
end

"""How a message names face `face` of `m`, as the class-3 messages do: nothing for the
interior, else `at zero ` and its zeroed names (`_face_names`), then a comma."""
_face_phrase(m::Mechanism, face::UInt64) =
    face == 0 ? "" : "at zero " * join(String.(_face_names(m, face)), ", ", " and ") * ", "

"""
    _face_condition(net, law_nodes, face, x) -> _Verdict

Admission's face condition, on face `face` of the ε-network `net` of an allosteric
mechanism, whose law network holds the nodes `law_nodes`; `x` is the `change` row of its
first substrate (`_turnover_name`). In order:
- `_EPS_NONUNIQUE`: the alive transitions leave two closed classes, so at every ε the
  steady state depends on where the enzyme starts; the verdict names the closed class
  outside the law network (`_class_outside_law`);
- `_PUMPED`: a stranded basin (`_face_sets`) is entered by, or holds, an alive transition
  whose reverse is dead, or is entered from outside it and the law network, so its weight
  depends on how fast its exits vanish; a basin entered only reversibly from the law
  network is sealed: its weight vanishes, and it leaves the fed nodes;
- `_BRIDGE`: in some conformation the vanishing transitions shorten the cheapest spanning
  in-tree of the fed nodes (R1 fails, `_shortens_in_tree`) and the first substrate's uptake
  is no potential on the strongly connected parts of its alive transitions that are not
  vanishing or are fast (T0+ fails, `_potential_on_sccs`), led by the first of its
  vanishing steps on a cheapest tree (`_bridging_edges`); or both conformations pass R1 and
  their vanishing transitions together shorten it (joint R1, no conformation or step).
A face whose fed nodes no in-tree spans passes: class 3 decides it. The verdict leaves the
face index 0 (`_admission_face` sets it).
"""
function _face_condition(net::_Network, law_nodes::BitVector, face::UInt64, x::Int)
    s = _face_sets(net, law_nodes, face)
    closed = _closed_classes(net, s.alive)
    length(closed) > 1 &&
        return _Verdict(_EPS_NONUNIQUE, _LIMIT_UNTESTED, 0,
                        _copy_of(net, _class_outside_law(closed, law_nodes)[1]), 0)
    for t in eachindex(net.from)
        u, w = net.from[t], net.to[t]
        s.alive[t] && s.fed[u] && s.stranded[w] || continue
        s.alive[isodd(t) ? t + 1 : t - 1] && (s.stranded[u] || law_nodes[u]) && continue
        return _Verdict(_PUMPED, _LIMIT_UNTESTED, 0, net.conformation[t], net.edge[t])
    end
    nodes = s.fed .& .!s.stranded
    _fewest_slow(net, nodes, s.n0) < 0 && return _WELL_POSED
    passed = 0
    for c in 1:net.ncopies
        (; van, held) = _copy_transitions(net, s.alive, nodes, c)
        any(van) || continue
        if !_shortens_in_tree(net, nodes, s.n0, van)
            passed += 1
        elseif !_potential_on_sccs(net, held, x)
            return _Verdict(_BRIDGE, _LIMIT_UNTESTED, 0, c,
                            first(_bridging_edges(net, nodes, s.n0, van)))
        end
    end
    passed == 2 && _shortens_in_tree(net, nodes, s.n0, s.alive .& net.vanishing) &&
        return _Verdict(_BRIDGE, _LIMIT_UNTESTED, 0, 0, 0)
    _WELL_POSED
end

"""
    _admission_face(am, eps, law, faces, i) -> _Verdict

Admission's face condition (`_face_condition`) on face `faces[i]` of `am`, from its
ε-network `eps` and law network `law`; a rejection carries the face index `i`.
"""
function _admission_face(am::AllostericMechanism, eps::_Network, law::_Network,
                         faces::Vector{UInt64}, i::Int)
    v = _face_condition(eps, law.present, faces[i], _turnover_name(am.catalytic))
    v.fault === _NO_FAULT ? v : _Verdict(v.fault, v.limit, i, v.conformation, v.edge)
end

"""
    _admission(am) -> _Verdict

Admission alone (spec section 7): the interior checks (`_admission_interior`), then the
face condition on every face (`_admission_face`), class 3 not run, as the tagging and
change-tag moves filter their children.
"""
function _admission(am::AllostericMechanism)
    v = _admission_interior(am)
    v.fault === _NO_FAULT || return v
    eps, law = _eps_network(am), _law_network(am)
    faces = _faces(am.catalytic)
    for i in eachindex(faces)
        v = _admission_face(am, eps, law, faces, i)
        v.fault === _NO_FAULT || return v
    end
    _WELL_POSED
end

# ─── Class 5: the support test on the derived law (spec section 8) ─────

"""
    _law_branches(m) -> Vector{Tuple{POLY, POLY}}
    _law_branches(den_A, den_I, d_A, d_I) -> Vector{Tuple{POLY, POLY}}

The derived denominator of `m`'s rate law as branches of positive terms, so that it keeps
a monomial free of a face's names exactly when some branch's two factors both do: a plain
mechanism's denominator in lowest terms (`_raw_symbolic_rate_polys`) times 1; an allosteric
one's two conformations, each denominator times the other conformation's free-enzyme
weight (`_mwc_state_polys`). Powers of a branch and the regulatory sites' factors, whose 1
never vanishes, keep the same faces. The second form pairs the two conformations'
denominators `den_A`, `den_I` and free-enzyme weights `d_A`, `d_I` that the allosteric
derivation already holds (`_num_den_exprs`); the first form derives them.
"""
_law_branches(m::Mechanism) = [(_raw_symbolic_rate_polys(m)[2], poly_one())]
function _law_branches(am::AllostericMechanism)
    _, den_A, _, den_I, d_A, d_I =
        _mwc_state_polys(am, Set{Symbol}(metabolites(am.catalytic)), _kept_members(am))
    _law_branches(den_A, den_I, d_A, d_I)
end
_law_branches(den_A::POLY, den_I::POLY, d_A::POLY, d_I::POLY) =
    [(den_A, d_I), (den_I, d_A)]

"""Whether some monomial of `p` holds none of the names `z` to a positive power."""
function _keeps_monomial(p::POLY, z::Vector{Symbol})
    for mono in keys(p)
        any(((s, e),) -> e > 0 && s in z, mono) || return true
    end
    false
end

"""The first face, each given by its zeroed names, on which no branch of `branches`
(`_law_branches`) keeps a monomial in both factors; 0 when every face keeps one."""
function _vanishing_face(branches::Vector{Tuple{POLY, POLY}},
                         faces::Vector{Vector{Symbol}})
    for (i, z) in enumerate(faces)
        any(b -> _keeps_monomial(b[1], z) && _keeps_monomial(b[2], z), branches) ||
            return i
    end
    0
end

"""
    _support_verdict(m, branches) -> _Verdict

Class 5 (spec section 8): `_DENOMINATOR_VANISHES` on the first face of `m` (`_faces`)
where the denominator `branches` (`_law_branches`) loses every monomial, so the law is 0/0
there; it catches the steady-state dead ends class 3's structure misses.
"""
function _support_verdict(m::Mechanism, branches::Vector{Tuple{POLY, POLY}})
    i = _vanishing_face(branches, [_face_names(m, f) for f in _faces(m)])
    i == 0 ? _WELL_POSED : _Verdict(_DENOMINATOR_VANISHES, _LIMIT_UNTESTED, i, 0, 0)
end
_support_verdict(am::AllostericMechanism, branches::Vector{Tuple{POLY, POLY}}) =
    _support_verdict(am.catalytic, branches)

"""
    _assert_support(m, branches)

Error with class 5's message (`_guard_message`) when the derived denominator of `m`, given
as `branches` (`_law_branches`), loses every monomial on some face (`_support_verdict`):
the law is 0/0 there. Both `_num_den_exprs` methods call it with the polynomials they hold.
"""
function _assert_support(m::_AnyMech, branches::Vector{Tuple{POLY, POLY}})
    v = _support_verdict(m, branches)
    v.fault === _NO_FAULT || error(_guard_message(m, v))
    nothing
end

"""The message of class 5's verdict `v` on a mechanism that `who` names and whose catalytic
graph is `m`: the face where every term of the law's denominator vanishes."""
_support_message(who::String, m::Mechanism, v::_Verdict) =
    "$who: $(_face_phrase(m, _faces(m)[v.face]))every term of the rate law's " *
    "denominator vanishes, so the rate is 0/0 there."

# ─── Guard messages ────────────────────────────────────────────────────

"""
    _guard_message(mech, v) -> String

The error for the guard's verdict `v` on `mech`:
`<Mechanism|AllostericMechanism>: <cause, naming each step by its line>. <fix>.` A fix is
printed only where a test shows that its variant passes the guard.
"""
function _guard_message(mech::_AnyMech, v::_Verdict)
    who = mech isa AllostericMechanism ? "AllostericMechanism" : "Mechanism"
    v.fault === _NO_CYCLE &&
        return "$who: no cycle of steps turns substrates into products, so the rate is " *
               "zero everywhere."
    v.fault === _RE_TURNOVER && return _re_turnover_message(who, mech)
    mech isa AllostericMechanism && v.fault in _ADMISSION_FAULTS &&
        return _admission_message(mech, v)
    v.fault === _DENOMINATOR_VANISHES && return _support_message(who, catalytic(mech), v)
    v.fault === _NO_FAULT && error("_guard_message: verdict $v holds no fault")
    _face_message(who, mech, v)
end

"""Whether making any one step of `cycle`, transitions of `mech`'s ε-network `net`,
steady-state with its kinetic group leaves no fast cycle that turns the reaction over: then
class 2's fix holds whichever step of the cycle is picked."""
function _any_step_cures(mech::_AnyMech, net::_Network, cycle::Vector{Int})
    m = catalytic(mech)
    members = _group_members(m)
    for t in cycle
        re = _flags(m)
        re[members[m.group[net.edge[t]]]] .= false
        isempty(_re_turnover_transitions(_eps_network(mech, re))) || return false
    end
    true
end

"""Class 2's message: the steps of the first fast cycle that turns the reaction over,
each tagged binding among them marked with the conformation where its affinity
vanishes, and the fix where any one of those steps steady-state cures it
(`_any_step_cures`)."""
function _re_turnover_message(who::String, mech::_AnyMech)
    m = catalytic(mech)
    net = _eps_network(mech)
    cycle = _re_turnover_transitions(net)
    c = net.conformation[cycle[1]]
    state = c == 2 ? "inactive" : "active"
    lines = String[]
    for t in cycle
        k = net.edge[t]
        line = _step_text(m, k)
        mech isa AllostericMechanism && _excludes(mech, k, c) &&
            (line *= " (whose affinity vanishes in the $state conformation)")
        push!(lines, line)
    end
    place = mech isa AllostericMechanism && c == 2 ? " in the inactive conformation" : ""
    msg = "$who: the rapid-equilibrium steps $(join(lines, ", ", " and ")) close a " *
          "cycle$place that turns the reaction over, so the mechanism has no finite rate."
    _any_step_cures(mech, net, cycle) ? msg * " Make one step of the cycle steady-state." :
                                        msg
end

"""Where a class-3 fault of `mech` lies: nothing for a plain mechanism; the conformation
`c` of an allosteric one."""
_place(::Mechanism, ::Int) = ""
_place(::AllostericMechanism, c::Int) =
    c == 2 ? " in the inactive conformation" : c == 1 ? " in the active conformation" : ""

"""
The message of a class-3 verdict `v` on `mech`: the face, the
forms and steps at fault, and a fix only for a trap and a bottomless segment, whose
variants the tests show passing.
"""
function _face_message(who::String, mech::_AnyMech, v::_Verdict)
    m = catalytic(mech)
    net = _law_network(mech)
    face = _faces(m)[v.face]
    s = _face_structure(net, face)
    zeroed = "at zero " * join(String.(_face_names(m, face)), ", ", " and ")
    at = zeroed * _place(mech, v.conformation)
    forms(nodes, label = _ -> "") =
        "{" * join([_form_text(m, _form_of(net, u)) * label(u) for u in nodes], ", ") * "}"
    lines(edges) = join([_step_text(m, k) for k in edges], ", ", " and ")
    if v.fault === _NONUNIQUE
        C = _away_class(s)
        return "$who: $at, the forms $(forms(C)) keep enzyme away from free E: every " *
               "step that leaves them takes up a zeroed metabolite " *
               "($(lines(_exits(net, s, C)))), so the steady state is not unique."
    elseif v.fault === _BOTTOMLESS
        return "$who: $at, the rapid-equilibrium segment " *
               "$(forms(s.segments[s.bottomless])) has no form free of the zeroed " *
               "metabolites, so the rate is 0/0 there. Make one of the segment's binding " *
               "steps steady-state."
    elseif v.fault === _TRAP
        trap, edge, joined = _trap_feeder(net, s)
        joined &&
            return "$who: $at, the forms $(forms(s.fast_classes[trap])) hold enzyme, but " *
                   "the rapid-equilibrium step $(_step_text(m, edge)) that joins them to " *
                   "the rest of their segment takes up a zeroed metabolite, so the rate " *
                   "law misweights them. Make that step's kinetic group steady-state."
        return "$who: $at, the rapid-equilibrium step $(_step_text(m, edge)) fills the " *
               "forms $(forms(s.fast_classes[trap])), which no rapid-equilibrium step " *
               "running there empties, so the rate law misweights them. Make that " *
               "step's kinetic group steady-state."
    elseif v.fault === _SPLIT
        split = _split_parts(net, s)
        # Parts that lie in different conformations name each form's conformation.
        if allequal(_copy_of(net, u) for part in split.parts for u in part)
            parts = forms.(split.parts)
        else
            at = zeroed
            state(u) = _copy_of(net, u) == 1 ? " (active)" : " (inactive)"
            parts = [forms(part, state) for part in split.parts]
        end
        return "$who: $at, the forms $(join(parts, " and ")) of one " *
               "rapid-equilibrium segment are joined only by steps that take up a " *
               "zeroed metabolite ($(lines(split.joins))), so the rate law cannot weight " *
               "one against the other."
    elseif v.fault === _LAW_MULTI || v.fault === _MULTILEVEL
        why = v.limit === _LIMIT_LAW_NAN ?
                  "the rate law has no value (0/0), while the steady state does" :
              v.limit === _LIMIT_SPEED ?
                  "the steady state depends on how fast the rapid-equilibrium steps " *
                  "run, which the rate law does not carry" :
                  "the steady state's fast-equilibrium limit differs from the rate law"
        return "$who: $at, $why."
    end
    error("_face_message: no message for $(v.fault)")
end
