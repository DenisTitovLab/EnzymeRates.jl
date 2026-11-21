# ABOUTME: The allosteric conformations of a mechanism: the forms the steps a conformation
# ABOUTME: keeps join, and the readings of how each switches off the steps it excludes.

"""The steps conformation `c` of `m` (1 active, 2 inactive) excludes by their kinetic
group's allosteric state, one flag per catalytic step: the active conformation those of
`:OnlyI` groups, the inactive one those of `:OnlyA` groups. A plain mechanism excludes
none."""
function _excluded_steps(am::Mechanisms.AllostericMechanism, c::Int)
    state = c == 1 ? :OnlyI : :OnlyA
    (e -> am.allosteric_states[e.group] === state).(am.catalytic.steps)
end
_excluded_steps(m::Mechanisms.Mechanism, ::Int) = falses(length(m.steps))

"""
    exclusion_components(am, c::Symbol) -> Vector{Int}

Each form's component under the edges of `am`'s catalytic graph that conformation `c`
(`:active` or `:inactive`) keeps (all but `_excluded_steps`, the edges of `:OnlyI`
groups in the active conformation, those of `:OnlyA` groups in the inactive), every form
included, labelled by its smallest form index, so free E's label is 1. A form outside free
E's component is cut off in `c`; the guard's built readings of how `c` switches off each
step it excludes (`_built_reading`) read these components.
"""
function exclusion_components(am::Mechanisms.AllostericMechanism, c::Symbol)
    excluded = _excluded_steps(am, c === :active ? 1 : c === :inactive ? 2 :
        throw(ArgumentError("a conformation is :active or :inactive, not :$c")))
    m = am.catalytic
    parent = collect(eachindex(m.enzyme_forms))
    for (k, e) in enumerate(m.steps)
        excluded[k] && continue
        a, b = _union_find_root(parent, e.from), _union_find_root(parent, e.to)
        a == b || (parent[max(a, b)] = min(a, b))
    end
    [_union_find_root(parent, f) for f in eachindex(m.enzyme_forms)]
end

function _union_find_root(parent::Vector{Int}, i::Int)
    while parent[i] != i
        i = parent[i]
    end
    i
end

"""
How a reading switches off a catalytic step in one conformation (`_Reading`): the
conformation keeps the step (`_RUNS`); the step loses the rate of its stored direction
(`_LOSES_FORWARD`) or of its reverse (`_LOSES_REVERSE`), the direction that runs uphill,
into the step's higher side, while the other keeps its rate, so a rapid-equilibrium step
stays fast; it is stopped, both its rates vanishing (`_STOPPED`); or it cannot be switched
off (`_CANNOT_SWITCH_OFF`).
"""
@enum _SwitchOff::UInt8 _RUNS _LOSES_FORWARD _LOSES_REVERSE _STOPPED _CANNOT_SWITCH_OFF

"""How a mechanism's conformations switch off its catalytic steps (`_SwitchOff`):
`reading[k, c]` for step `k` in conformation `c`, 1 active and 2 inactive. A plain
mechanism's reading has its one conformation only."""
const _Reading = Matrix{_SwitchOff}

"""
Which excluded steady-state steps that take up one metabolite into a cut-off form lose
their affinity, the rest being stopped, in the guard's built readings (`_built_reading`):
every one (`_ALL_LOSE_AFFINITY`); those that lead away from free E
(`_OUTWARD_LOSE_AFFINITY`); or none (`_ALL_STOPPED`). The guard tries them in this order,
before it searches the other readings (`_guard_with_reading`).

The steps a conformation keeps join the forms into parts. A part's distance is the fewest
excluded steps that take up one metabolite into a cut-off form, rapid-equilibrium or
steady-state, each followed from its unbound side to its bound side, that lead to it from
free E's part; a part they do not reach is infinitely far. Such a step leads away from
free E when its bound side's part is farther than its unbound side's.
"""
@enum _SlowUptakeReading::UInt8 _ALL_LOSE_AFFINITY _OUTWARD_LOSE_AFFINITY _ALL_STOPPED

"""
    _built_reading(am, r = _ALL_LOSE_AFFINITY) -> _Reading
    _built_reading(m::Mechanisms.Mechanism) -> _Reading

How the conformations of `am` switch off its catalytic steps in the guard's built reading
`r`, read from the reaction graph alone. A conformation keeps a step it does not exclude.
An excluded step that takes up one metabolite and gives off none into a cut-off form, the
steps the conformation keeps (`exclusion_components`) joining neither its two sides nor
its bound side to free E, loses its affinity at rapid equilibrium: its stored direction,
the uptake, vanishes (`_LOSES_FORWARD`). At steady state it loses its affinity or is
stopped as `r` says (`_SlowUptakeReading`). Any other excluded step keeps its kind, rapid
equilibrium or steady state, as every step does in both conformations: a steady-state one
is stopped, and a rapid-equilibrium one is not switched off (`_CANNOT_SWITCH_OFF`).
Thermodynamics pins the equilibrium constant only of such a step inside one group; one
between two groups could lose a direction, as the search tries (`_consistent_columns`),
but the built readings let only an uptake lose. A plain mechanism `m` has one conformation,
which runs every step.
"""
function _built_reading(am::Mechanisms.AllostericMechanism,
                        r::_SlowUptakeReading = _ALL_LOSE_AFFINITY)
    steps = am.catalytic.steps
    reading = fill(_RUNS, length(steps), 2)
    for c in 1:2
        excluded = _excluded_steps(am, c)
        any(excluded) || continue
        comp = exclusion_components(am, _CONFORMATIONS[c])
        uptake = [excluded[k] && length(e.taken_up) == 1 && isempty(e.given_off) &&
                  comp[e.from] != comp[e.to] && comp[e.to] != comp[1]
                  for (k, e) in enumerate(steps)]
        # Each part's distance from free E's, breadth-first along the uptakes; a part's
        # label is its smallest form index.
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
        for (k, e) in enumerate(steps)
            excluded[k] || continue
            loses = uptake[k] &&
                    (e.rapid_equilibrium || r === _ALL_LOSE_AFFINITY ||
                     r === _OUTWARD_LOSE_AFFINITY &&
                     distance[comp[e.to]] > distance[comp[e.from]])
            reading[k, c] = loses ? _LOSES_FORWARD :
                            e.rapid_equilibrium ? _CANNOT_SWITCH_OFF : _STOPPED
        end
    end
    reading
end
_built_reading(m::Mechanisms.Mechanism) = fill(_RUNS, length(m.steps), 1)
