# ABOUTME: The graph accessors of a mechanism that compute something: what a step does and
# ABOUTME: which steps share a kinetic group.

"""
What step `k` of `m` does: `:isomerization` when it takes up and gives off nothing,
`:binding` when it is a plain binding (`_is_plain_binding`), `:fused_binding` when it
takes up one metabolite and gives off none but its two forms differ by more than that
metabolite, else `:chemistry`.
"""
function step_kind(m::Mechanism, k::Int)
    e = m.steps[k]
    _is_iso(e) ? :isomerization : _is_plain_binding(m, e) ? :binding :
        _is_binding(e) ? :fused_binding : :chemistry
end

"""Each kinetic group's steps, in step order."""
function groups(m::Mechanism)
    members = [Int[] for _ in 1:maximum(e -> e.group, m.steps; init = 0)]
    for (k, e) in enumerate(m.steps)
        push!(members[e.group], k)
    end
    members
end
