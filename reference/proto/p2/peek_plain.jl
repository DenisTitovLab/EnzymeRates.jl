# ABOUTME: Prints the structure (steps, groups, fitted, dependents) of plain mechanisms of a population.
# ABOUTME: Usage: julia peek_plain.jl <population> [indices...]; throwaway P2 probe helper.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
pop = Pop.load_population(ARGS[1])
idx = length(ARGS) > 1 ? parse.(Int, ARGS[2:end]) : eachindex(pop)
for i in idx
    m = pop[i]
    m.allosteric && continue
    println("== ", m.id, "  ", m.source, "  fam=", m.family)
    println("   subs=", m.substrates, " prods=", m.products, " inh=", m.inhibitors)
    g = Pop.group_of_step(m)
    for (j, s) in enumerate(m.steps)
        println("   $j [g$(g[j])] ", Pop.step_text(m, s), "   ", m.step_names[j])
    end
    println("   fitted(", length(m.fitted_params), ")=", m.fitted_params, " ipc=", m.independent_param_count)
    for d in m.dependents; println("   dep ", d); end
    isempty(m.pass2_rename) || println("   pass2=", m.pass2_rename)
end
