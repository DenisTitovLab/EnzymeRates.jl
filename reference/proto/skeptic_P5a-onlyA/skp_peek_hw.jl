# ABOUTME: Skeptic P5a: prints the plain handwritten mechanisms with conformations or pins.
# ABOUTME: Stdlib only.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
for m in Pop.load_population("handwritten")
    (:conformations in m.family || occursin("pins", m.source)) || continue
    println(m.id, "  ", m.source, " fam=", m.family)
    g = Pop.group_of_step(m)
    for (j, s) in enumerate(m.steps); println("     g", g[j], " ", Pop.step_text(m, s)); end
    println("   dependents: ", m.dependents)
end
