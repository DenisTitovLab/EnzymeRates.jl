# ABOUTME: Survey of the allosteric mechanisms in the shared populations (P5a probe).
# ABOUTME: Prints structure, tags, sites and package outputs of selected allosteric records.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
function show_allo(m; polys = false)
    println("== ", m.id, "  ", m.source, " n=", m.cat_multiplicity, " family=", m.family)
    println("   subs=", m.substrates, " prods=", m.products, " inh=", m.inhibitors, " regs=", m.regulators, " mets=", m.metabolites)
    gos = Pop.group_of_step(m)
    for (j, s) in enumerate(m.steps)
        println("   step $j g$(gos[j]) [$(m.group_tags[gos[j]])]: ", Pop.step_text(m, s), "  A=", m.step_names_A[j], " I=", m.step_names_I[j], " Ipkg=", m.step_names_I_pkg[j], " kept=", m.inactive_kept[j])
    end
    println("   sites=", m.sites)
    println("   fitted=", m.fitted_params)
    println("   deps=", m.dependents)
    println("   pass2A=", m.pass2_rename_A, " pass2I=", m.pass2_rename_I)
    if polys
        for (nm, p) in zip(("numA", "denA", "dfreeA"), m.state_polys_A); println("   $nm: ", p); end
        for (nm, p) in zip(("numI", "denI", "dfreeI"), m.state_polys_I); println("   $nm: ", p); end
    end
    for s in m.samples[1:1]
        println("   sample v=", s.v, " kcat=", s.kcat, " params=", s.params, " concs=", s.concs)
    end
end
for popname in ARGS[1:1]
    pop = Pop.load_population(popname)
    idx = length(ARGS) > 1 ? parse.(Int, ARGS[2:end]) : findall(m -> m.allosteric, pop)
    for i in idx
        show_allo(pop[i]; polys = length(ARGS) > 1)
    end
end
