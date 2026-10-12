# ABOUTME: Quick census of the allosteric mechanisms the probe uses: sources, tags, free forms, ping-pong.
# ABOUTME: Stdlib only (probe P5b).
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
for name in ("specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi")
    pop = filter(m -> m.allosteric, Pop.load_population(name))
    nfree2 = [m.id for m in pop if count(Pop.is_free, m.species) >= 2]
    pp = count(m -> :pingpong in m.family, pop)
    println(name, ": ", length(pop), " allosteric; >=2 free forms: ", nfree2, "; pingpong ", pp,
            "; n=", unique(m.cat_multiplicity for m in pop), "; sites ", count(m -> !isempty(m.sites), pop))
    name == "handwritten" && foreach(m -> println("   ", m.id, " ", m.source, " tags=", m.group_tags,
                                                  " re=", [s.re for s in m.steps], " fitted=", m.fitted_params), pop)
end
