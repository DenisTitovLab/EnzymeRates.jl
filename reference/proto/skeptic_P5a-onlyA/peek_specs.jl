# ABOUTME: Skeptic P5a: prints the allosteric specs/handwritten mechanisms with tags and families.
# ABOUTME: Stdlib only; reads the shared populations.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
for name in ("specs", "handwritten", "allo_uniuni")
    for m in Pop.load_population(name)
        m.allosteric || continue
        println(m.id, "  ", m.source, "  n=", m.cat_multiplicity, " tags=", m.group_tags, " fam=", m.family)
        name == "allo_uniuni" && continue
        g = Pop.group_of_step(m)
        for (j, s) in enumerate(m.steps)
            println("     g", g[j], " ", Pop.step_text(m, s))
        end
        for st in m.sites
            println("     site ", st.ligands, " x", st.multiplicity, " ", st.states)
        end
    end
end
