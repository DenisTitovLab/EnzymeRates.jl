# ABOUTME: The well-posedness check the mechanism macros, the parser, the enumerator's
# ABOUTME: exits and the derivation entries run; constructors check representation only.

"""
    _check_well_posed(m)

Error when a rapid-equilibrium segment of `m` has no bottom form
(`_assert_bottomed_segments`), in the catalytic graph and, for an allosteric mechanism,
in the inactive conformation's graph (`_inactive_mask`). An
allosteric mechanism must also keep all its free forms in one rapid-equilibrium segment
of the edges outside `:OnlyA` groups (`_assert_one_flip_segment`) and hold no `:OnlyA`
binding that leaves a thermodynamic cycle unsatisfiable (`_onlya_haldane_violation`).
It reads the reaction graph through masks and builds nothing. The mechanism macros and
the parser (`_parse_mechanism`) call it through `_checked_mechanism`; `init_mechanisms`,
`seed_mechanisms`, `expand_mechanisms`, identify's base tier and every derivation entry
(`_derivation_input`) call it too; the constructors and the Sig decoder do not, so
code that builds a mechanism internally never trips a check meant for the mechanisms
it returns.
"""
function _check_well_posed(m::Mechanism)
    _assert_bottomed_segments(m, trues(length(m.edges)))
    nothing
end
function _check_well_posed(am::AllostericMechanism)
    m = am.catalytic
    _assert_bottomed_segments(m, trues(length(m.edges)))
    _assert_one_flip_segment(am)
    offenders = _onlya_haldane_violation(am)
    isempty(offenders) ||
        error("AllostericMechanism: ", _onlya_violation_message(am, offenders))
    _assert_bottomed_segments(m, _inactive_mask(am))
    nothing
end

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

"""
Error when a rapid-equilibrium segment of the edges of `m` that `mask` keeps has no
bottom form (`_bottomless_segment`): its rate is 0/0 where one side's metabolites are
absent. The boundary check (`_check_well_posed`) checks a mechanism's graph and, for an
allosteric mechanism, the inactive conformation's (`_inactive_mask`).
"""
function _assert_bottomed_segments(m::Mechanism, mask::BitVector)
    segment = _bottomless_segment(m, mask, _flags(m))
    isempty(segment) ||
        error("Mechanism: rapid-equilibrium segment {" *
              join([_form_name(m, f) for f in segment], ", ") *
              "} has no form free of one side's " *
              "metabolites, so its rate is undefined when that side is absent " *
              "(0/0). Make one of the segment's binding steps steady-state.")
    nothing
end

"""
Error unless every free form (no bound metabolite, no residual) of the allosteric
mechanism `am`'s catalytic graph lies in one rapid-equilibrium segment of the edges
outside `:OnlyA` groups. Under formulation 1 the conformations interconvert at a free
form. Within one such segment, flipping at any of its free forms gives the same model
with L rescaled; free forms joined any other way flip at different forms in the two
conformations or leave one stranded in the inactive conformation, a different model.
"""
function _assert_one_flip_segment(am::AllostericMechanism)
    m = am.catalytic
    free = [f for f in eachindex(m.forms) if _is_free(m.forms[f])]
    length(free) <= 1 && return nothing
    segment_of = _segments(m, trues(length(m.edges)),
                           _flags(m) .& _outside_onlya(am)).segment_of
    all(f -> segment_of[f] == segment_of[first(free)], free) || error(
        "AllostericMechanism: the free forms " *
        "$(join([_form_name(m, f) for f in free], ", ")) " *
        "are not joined by rapid-equilibrium steps outside :OnlyA groups. The " *
        "conformations interconvert at one free form, so the other free forms must be " *
        "in rapid equilibrium with it: join them by ⇌ steps outside :OnlyA groups, or " *
        "keep one free form.")
    nothing
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
turns it over (`_RE_TURNOVER`, class 2).
"""
@enum _Fault::UInt8 begin
    _NO_FAULT
    _NO_CYCLE
    _RE_TURNOVER
end

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

"""The metabolite names that face `face` of `m` zeroes (`_faces`), in table order."""
function _face_names(m::Mechanism, face::UInt64)
    name_of = _name_index(m)
    Symbol[x.name for (i, x) in enumerate(m.table)
           if (i == 1 || name_of[i] != name_of[i - 1]) &&
              (face >> (name_of[i] - 1)) & 1 == 1]
end

"""The metabolite-name index (`_name_index`) of `m`'s first substrate, whose net uptake
is the rate."""
_turnover_name(m::Mechanism) =
    _name_index(m)[searchsortedfirst(m.table, Met(m.substrates[1], false))]

"""
Whether a cycle of the edges `mask` keeps turns the reaction over: no potential of the
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
The transitions of a cycle of fast transitions of `net` that takes up some metabolite name
net (class 2: the rate is infinite), in cycle order; empty when none does. Every name's
offsets along a breadth-first forest of the fast transitions find the first transition
they disagree on, which closes the cycle with its tree path. Class 2 runs once, at the
interior, since a fast cycle that runs on a face runs at the interior too.
"""
function _re_turnover_transitions(net::_Network)
    nn = net.nforms * net.ncopies
    adjacent = [Int[] for _ in 1:nn]
    for t in eachindex(net.from)
        net.fast[t] && net.forward[t] || continue
        push!(adjacent[net.from[t]], t)
        push!(adjacent[net.to[t]], t)
    end
    offset = zeros(Int, size(net.change, 1), nn)
    parent = zeros(Int, nn)
    seen = falses(nn)
    for root in 1:nn
        net.present[root] && !seen[root] || continue
        seen[root] = true
        queue = [root]
        head = 0
        while head < length(queue)
            head += 1
            u = queue[head]
            for t in adjacent[u]
                t == parent[u] && continue
                v, dir = net.from[t] == u ? (net.to[t], 1) : (net.from[t], -1)
                if !seen[v]
                    seen[v] = true
                    parent[v] = t
                    push!(queue, v)
                    for n in axes(offset, 1)
                        offset[n, v] = offset[n, u] + dir * net.change[n, t]
                    end
                else
                    for n in axes(offset, 1)
                        offset[n, v] == offset[n, u] + dir * net.change[n, t] ||
                            return _tree_cycle(net, parent, u, v, t)
                    end
                end
            end
        end
    end
    Int[]
end

"""
    _guard(m) -> _Verdict

The well-posedness guard (spec section 8) on `m`, plain or allosteric: class 1 (a cycle
turns the reaction over, `_has_turnover_cycle`), then class 2 (no cycle of fast
transitions does, `_re_turnover_transitions`, on the ε-network, where an excluded
rapid-equilibrium binding stays fast while its affinity vanishes). `_WELL_POSED` when
every check passes. It reads the reaction graph through masks and builds nothing.
"""
function _guard(m::Mechanism)
    _has_turnover_cycle(m, trues(length(m.edges))) ||
        return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
    _re_turnover_verdict(_eps_network(m))
end
function _guard(am::AllostericMechanism)
    m = am.catalytic
    _has_turnover_cycle(m, trues(length(m.edges))) ||
        return _Verdict(_NO_CYCLE, _LIMIT_UNTESTED, 0, 0, 0)
    _re_turnover_verdict(_eps_network(am))
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
    error("_guard_message: verdict $v holds no fault")
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
