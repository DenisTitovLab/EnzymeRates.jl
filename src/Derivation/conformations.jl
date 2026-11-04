# ABOUTME: The allosteric conformations of a mechanism: the forms the steps a conformation
# ABOUTME: keeps join (`exclusion_components`).

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
index, so free E's label is 1. Admission reads each plain binding `c` excludes as an arc
from its unbound side's component to its bound side's (`_arcs`).
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
