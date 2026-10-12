# ABOUTME: Card C4 follow-up: prints the steps, groups and dependents of chosen plain mechanisms.
# ABOUTME: Usage: julia peek_ids.jl <id...> (ids such as specs/7).
include(joinpath(@__DIR__, "..", "..", "proto", "common", "Population.jl"))
for id in ARGS
    pi = split(id, "/"); m = Pop.load_population(pi[1])[parse(Int, pi[2])]
    println("== ", m.id, "  ", m.source, "  family=", m.family)
    for (j, s) in enumerate(m.steps)
        println("  step $j: ", Pop.step_text(m, s), "   names=", m.step_names[j])
    end
    println("  groups=", m.groups, "  deps=", m.dependents)
end
