# ABOUTME: Peeks at Pop records: names of residual, inhibitor-copy and allosteric steps.
# ABOUTME: Scratch exploration for the core stage of the vertical slice.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common/Population.jl")
for pop in ("specs", "handwritten", "init_bibi_pp", "allo_uniuni")
    ms = Pop.load_population(pop)
    println("== ", pop, " ", length(ms))
    shown = 0
    for m in ms
        interesting = any(sp -> !isempty(sp.res_added) || any(x -> x.role === :Inh, sp.bound), m.species) || !isempty(m.sites)
        interesting || continue
        shown += 1; shown > 3 && break
        println("-- ", m.id, " ", m.source, " allo=", m.allosteric, " tags=", m.group_tags, " sites=", m.sites)
        println("   subs=", m.substrates, " prods=", m.products, " inh=", m.inhibitors, " regs=", m.regulators, " mets=", m.metabolites)
        g = Pop.group_of_step(m)
        for (j, s) in enumerate(m.steps)
            nm = m.allosteric ? (m.step_names_A[j], m.step_names_I[j]) : m.step_names[j]
            println("   ", j, " g", g[j], " ", Pop.step_text(m, s), "   ", nm)
        end
        println("   fitted: ", m.fitted_params)
    end
end
