# ABOUTME: Prints the structure and package outputs of chosen allosteric mechanisms (probe P5b).
# ABOUTME: Usage: julia peek.jl <population> <index...>
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
pop = Pop.load_population(ARGS[1])
idx = ARGS[2] == "allo" ? findall(m -> m.allosteric, pop) : parse.(Int, ARGS[2:end])
for i in idx
    m = pop[i]
    println("== ", m.id, "  ", m.source, "  n=", m.cat_multiplicity, " family=", m.family)
    println("subs=", m.substrates, " prods=", m.products, " inh=", m.inhibitors, " regs=", m.regulators, " mets=", m.metabolites)
    for (k, sp) in enumerate(m.species)
        println("  form $k: ", sp.name, " conf=", sp.conformation, " bound=", [x.name for x in sp.bound], " res=", sp.res_added, sp.res_subtracted)
    end
    g = Pop.group_of_step(m)
    for (j, s) in enumerate(m.steps)
        println("  step $j g$(g[j]) [$(m.group_tags[g[j]])]: ", Pop.step_text(m, s), "   A=", m.step_names_A[j], " I=", m.step_names_I[j], " Ipkg=", m.step_names_I_pkg[j], " kept=", m.inactive_kept[j], " rev=", m.inactive_reversed[j])
    end
    println("sites=", m.sites)
    println("fitted=", m.fitted_params)
    println("deps=", m.dependents)
    println("p2A=", m.pass2_rename_A, " p2I=", m.pass2_rename_I)
    println("dfreeA=", m.state_polys_A[3], "  dfreeI=", m.state_polys_I[3])
    println(m.rate_equation_string)
end
