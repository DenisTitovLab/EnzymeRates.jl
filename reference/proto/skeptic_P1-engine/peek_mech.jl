# ABOUTME: Prints the steps, groups, fitted parameters and dependents of chosen population mechanisms.
# ABOUTME: Usage: julia peek_mech.jl <population> <index...>
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
pop = Pop.load_population(ARGS[1])
for i in parse.(Int, ARGS[2:end])
    m = pop[i]
    println("== ", m.id, "  ", m.source, "  family=", m.family)
    for (g, grp) in enumerate(m.groups), j in grp
        println("  g$g step $j: ", Pop.step_text(m, m.steps[j]), "   ", m.step_names[j])
    end
    println("  fitted=", m.fitted_params)
    println("  deps=", m.dependents)
end
