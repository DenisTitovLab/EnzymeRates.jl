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
    _switch_offs(am, c::Int) -> Vector{_SwitchOff}

How conformation `c` (1 active, 2 inactive) of `am` switches off each catalytic step, read
from the reaction graph alone. It keeps a step it does not exclude. An excluded step that
takes up one metabolite and gives off none into a cut-off form, the steps `c` keeps
(`exclusion_components`) joining neither its two sides nor its bound side to free E,
loses its affinity. Any other excluded step keeps its kind, rapid equilibrium or steady
state, as every step does in both conformations: a steady-state one is stopped, and a
rapid-equilibrium one, whose equilibrium constant thermodynamics pins, cannot be switched
off.
"""
function _switch_offs(am::Mechanisms.AllostericMechanism, c::Int)
    comp = exclusion_components(am, _CONFORMATIONS[c])
    excluded = _excluded_steps(am, c)
    off = fill(_RUNS, length(am.catalytic.steps))
    for (k, e) in enumerate(am.catalytic.steps)
        excluded[k] || continue
        into_cut_off = length(e.taken_up) == 1 && isempty(e.given_off) &&
                       comp[e.from] != comp[e.to] && comp[e.to] != comp[1]
        off[k] = into_cut_off ? _LOSES_AFFINITY :
                 e.rapid_equilibrium ? _CANNOT_SWITCH_OFF : _STOPPED
    end
    off
end
