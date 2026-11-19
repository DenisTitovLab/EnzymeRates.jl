# ABOUTME: The allosteric conformations of a mechanism: the forms the steps a conformation
# ABOUTME: keeps join, and how it switches off each step it excludes (`_switch_offs`).

"""The steps of conformation `c` of `m` (1 active, 2 inactive), one flag per catalytic step:
an allosteric mechanism's `active_steps` or `inactive_steps`; a plain mechanism has only the
active conformation, every step."""
_conformation_steps(am::Mechanisms.AllostericMechanism, c::Int) =
    (am.active_steps, am.inactive_steps)[c]
_conformation_steps(m::Mechanisms.Mechanism, c::Int) =
    (c == 1 ? trues : falses)(length(m.steps))

"""The steps conformation `c` of `m` (1 active, 2 inactive) excludes, one flag per
catalytic step: an allosteric mechanism's `excluded_active` or `excluded_inactive`; a plain
mechanism excludes none."""
_excluded_steps(am::Mechanisms.AllostericMechanism, c::Int) =
    (am.excluded_active, am.excluded_inactive)[c]
_excluded_steps(m::Mechanisms.Mechanism, ::Int) = falses(length(m.steps))

"""
    exclusion_components(am, c::Symbol) -> Vector{Int}

Each form's component under the edges of `am`'s catalytic graph that conformation `c`
(`:active` or `:inactive`) keeps (`am.excluded_active`, `am.excluded_inactive`: where
`Mechanisms.excludes(tag, c)` holds), every form included, labelled by its smallest form
index, so free E's label is 1. A form outside free E's component is cut off in `c`; how
`c` switches off each step it excludes (`_switch_offs`) reads these components.
"""
function exclusion_components(am::Mechanisms.AllostericMechanism, c::Symbol)
    excluded = _excluded_steps(am, c === :active ? 1 : c === :inactive ? 2 :
        throw(ArgumentError("a conformation is :active or :inactive, not :$c")))
    m = am.catalytic
    parent = collect(eachindex(m.nodes))
    for (k, e) in enumerate(m.steps)
        excluded[k] && continue
        a, b = _union_find_root(parent, e.from), _union_find_root(parent, e.to)
        a == b || (parent[max(a, b)] = min(a, b))
    end
    [_union_find_root(parent, f) for f in eachindex(m.nodes)]
end

function _union_find_root(parent::Vector{Int}, i::Int)
    while parent[i] != i
        i = parent[i]
    end
    i
end

"""
How a conformation switches off a catalytic step (`_switch_offs`): it keeps the step
(`_RUNS`); the step loses its affinity (`_LOSES_AFFINITY`): its uptake vanishes and its
release keeps its rate, so a rapid-equilibrium one stays fast; it is stopped, both its
rates vanishing (`_STOPPED`); or it cannot be switched off (`_CANNOT_SWITCH_OFF`).
"""
@enum _SwitchOff::UInt8 _RUNS _LOSES_AFFINITY _STOPPED _CANNOT_SWITCH_OFF

"""
Which excluded steady-state steps that take up one metabolite into a cut-off form lose
their affinity, the rest being stopped (`_switch_offs`): every one
(`_ALL_LOSE_AFFINITY`); one that leads away from free E (`_OUTWARD_LOSE_AFFINITY`); or
none (`_ALL_STOPPED`). The guard tries them in this order (`_guard`).

The steps a conformation keeps join the forms into parts. A part's distance is the fewest
excluded steps that take up one metabolite into a cut-off form, rapid-equilibrium or
steady-state, each followed from its unbound side to its bound side, that lead to it from
free E's part; a part they do not reach is infinitely far. Such a step leads away from
free E when its bound side's part is farther than its unbound side's.
"""
@enum _SlowUptakeReading::UInt8 _ALL_LOSE_AFFINITY _OUTWARD_LOSE_AFFINITY _ALL_STOPPED

"""
    _switch_offs(am, c::Int, reading = _ALL_LOSE_AFFINITY) -> Vector{_SwitchOff}

How conformation `c` (1 active, 2 inactive) of `am` switches off each catalytic step, read
from the reaction graph alone. It keeps a step it does not exclude. An excluded step that
takes up one metabolite and gives off none into a cut-off form, the steps `c` keeps
(`exclusion_components`) joining neither its two sides nor its bound side to free E, loses
its affinity at rapid equilibrium; at steady state it loses its affinity or is stopped as
`reading` says (`_SlowUptakeReading`). Any other excluded step keeps its kind, rapid
equilibrium or steady state, as every step does in both conformations: a steady-state one
is stopped, and a rapid-equilibrium one, whose equilibrium constant thermodynamics pins,
cannot be switched off.
"""
function _switch_offs(am::Mechanisms.AllostericMechanism, c::Int,
                      reading::_SlowUptakeReading = _ALL_LOSE_AFFINITY)
    comp = exclusion_components(am, _CONFORMATIONS[c])
    excluded = _excluded_steps(am, c)
    steps = am.catalytic.steps
    uptake = [excluded[k] && length(e.taken_up) == 1 && isempty(e.given_off) &&
              comp[e.from] != comp[e.to] && comp[e.to] != comp[1]
              for (k, e) in enumerate(steps)]
    # Each part's distance from free E's, breadth-first along the uptakes; a part's label
    # is its smallest form index.
    distance = fill(typemax(Int), length(comp))
    distance[comp[1]] = 0
    queue, head = [comp[1]], 0
    while head < length(queue)
        head += 1
        u = queue[head]
        for (k, e) in enumerate(steps)
            uptake[k] && comp[e.from] == u && distance[comp[e.to]] == typemax(Int) ||
                continue
            distance[comp[e.to]] = distance[u] + 1
            push!(queue, comp[e.to])
        end
    end
    off = fill(_RUNS, length(steps))
    for (k, e) in enumerate(steps)
        excluded[k] || continue
        loses = uptake[k] &&
                (e.rapid_equilibrium || reading === _ALL_LOSE_AFFINITY ||
                 reading === _OUTWARD_LOSE_AFFINITY &&
                 distance[comp[e.to]] > distance[comp[e.from]])
        off[k] = loses ? _LOSES_AFFINITY :
                 e.rapid_equilibrium ? _CANNOT_SWITCH_OFF : _STOPPED
    end
    off
end
