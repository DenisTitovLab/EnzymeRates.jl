# ABOUTME: The facts a mechanism's graph gives: what a step does and its kinetic groups,
# ABOUTME: computed once at construction, and the steps each conformation runs.

"""
What step `e` of `m` does: `:isomerization` when it takes up and gives off nothing,
`:binding` when it is a plain binding (`_is_plain_binding`), `:fused_binding` when it
takes up one metabolite and gives off none but its two forms differ by more than that
metabolite, else `:chemistry`.
"""
_step_kind(m::Mechanism, e::Step) =
    _is_iso(e) ? :isomerization : _is_plain_binding(m, e) ? :binding :
    _is_binding(e) ? :fused_binding : :chemistry

"""Each kinetic group's steps, in step order, each list stored at its length."""
function _groups(steps::Vector{Step})
    G = maximum(e -> e.group, steps; init = 0)
    members = [Vector{Int}(undef, count(e -> e.group == g, steps)) for g in 1:G]
    filled = zeros(Int, G)
    for (k, e) in enumerate(steps)
        members[e.group][filled[e.group] += 1] = k
    end
    members
end

"""
    steps_in_active_conformation(m) -> BitVector

The steps of `m` that run in its active conformation, one flag per catalytic step. A step
runs unless its kinetic group's allosteric state is `:OnlyI` or the step is cut off: the
conformations interconvert at free E alone, so a form that the running steps do not join
to free E holds no enzyme of that conformation, and a step between two such forms does
not run, whatever its state. Every step of a plain mechanism runs.
"""
steps_in_active_conformation(m::Mechanism) = trues(length(m.steps))
steps_in_active_conformation(am::AllostericMechanism) = _steps_in_conformation(am, :OnlyI)

"""
    steps_in_inactive_conformation(m) -> BitVector

The steps of `m` that run in its inactive conformation, one flag per catalytic step: as in
`steps_in_active_conformation`, with the steps of `:OnlyA` kinetic groups excluded in
place of those of `:OnlyI` groups. Every step of a plain mechanism runs.
"""
steps_in_inactive_conformation(m::Mechanism) = trues(length(m.steps))
steps_in_inactive_conformation(am::AllostericMechanism) = _steps_in_conformation(am, :OnlyA)

"""The steps of the catalytic graph of `am` that run in the conformation that excludes the
kinetic groups in allosteric state `state`: the other steps, between forms they join to
free E, form 1."""
function _steps_in_conformation(am::AllostericMechanism, state::Symbol)
    m = am.catalytic
    excluded = (e -> am.allosteric_states[e.group] === state).(m.steps)
    joined = _joined_to_free_E(m, excluded)
    (k -> !excluded[k] && joined[m.steps[k].from]).(eachindex(m.steps))
end

"""Which nodes of `m` the steps not in `excluded` join to free E, node 1."""
function _joined_to_free_E(m::Mechanism, excluded::BitVector)
    reached = falses(length(m.nodes))
    reached[1] = true
    grew = true
    while grew
        grew = false
        for (k, e) in enumerate(m.steps)
            (excluded[k] || reached[e.from] == reached[e.to]) && continue
            reached[e.from] = reached[e.to] = true
            grew = true
        end
    end
    reached
end
