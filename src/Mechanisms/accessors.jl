# ABOUTME: The facts a mechanism's graph gives (what a step does, which forms the steps of a
# ABOUTME: conformation join to free E), computed once at construction, and their reads.

"""What step `k` of `m` does (`_step_kind`; stored as `m.step_kinds`)."""
step_kind(m::Mechanism, k::Int) = m.step_kinds[k]

"""
What step `e` of `m` does: `:isomerization` when it takes up and gives off nothing,
`:binding` when it is a plain binding (`_is_plain_binding`), `:fused_binding` when it
takes up one metabolite and gives off none but its two forms differ by more than that
metabolite, else `:chemistry`.
"""
_step_kind(m::Mechanism, e::Step) =
    _is_iso(e) ? :isomerization : _is_plain_binding(m, e) ? :binding :
    _is_binding(e) ? :fused_binding : :chemistry

"""Each kinetic group's steps, in step order (stored as `m.groups`)."""
groups(m::Mechanism) = m.groups

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

"""
The steps of the conformation of `m` that excludes the steps `excluded`, one flag per step
in `m`'s stored direction: the other steps, between forms they join to free E, form 1. The
conformations interconvert at free E alone, so a form those steps do not join to it holds
no enzyme of that conformation, a second free form included. The law deletes the excluded
steps and prunes the forms this strands; admission decides when that is the limit of the
vanishing constants (spec section 7).
"""
function _conformation_steps(m::Mechanism, excluded::BitVector)
    joined = _joined_to_free_E(m, excluded)
    (k -> !excluded[k] && joined[m.steps[k].from]).(eachindex(m.steps))
end
