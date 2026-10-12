# ABOUTME: Lists every plain mechanism of chosen populations with its id, source and family tags.
# ABOUTME: Throwaway P4 helper; stdlib only.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
for pop in ARGS
    for m in Pop.load_population(pop)
        m.allosteric && continue
        println(m.id, "  ", m.source, "  ", join(m.family, ","), "  nsteps=", length(m.steps), " nforms=", length(m.species))
    end
end
