# ABOUTME: Peek at laws the common-factor diagnostic flags: steps and the normalized law terms.
# ABOUTME: Also tests the candidate factor (K + Q) by substituting Q = -K exactly.
include(joinpath(@__DIR__, "Laws.jl"))
ms = Pop.load_population("init_bibi")
m = only(filter(x -> x.id == "init_bibi/104", ms))
for (g, grp) in enumerate(m.groups), j in grp
    println("   g", g, " ", Pop.step_text(m, m.steps[j]), "   ", m.step_names[j])
end
b = Laws.plain_today(m)
comps = Canon.normalize(b.law)
for (i, c) in enumerate(comps)
    println("comp ", i, " (", length(c), " terms)")
    for (mo, co) in sort!(collect(c); by = string)
        println("   ", co, "  ", join(("$s^$e" for (s, e) in mo), " "))
    end
end
println(m.rate_equation_string)
