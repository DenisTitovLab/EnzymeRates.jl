# ABOUTME: Peek at the two init_bibi mechanism pairs that today's laws merge and the redesign's do not.
# ABOUTME: Prints their steps, groups, fitted names and both dependents (package and joint solve).
include(joinpath(@__DIR__, "Laws.jl"))
ms = Pop.load_population("init_bibi")
for id in ("init_bibi/69", "init_bibi/144", "init_bibi/79", "init_bibi/130")
    m = only(filter(x -> x.id == id, ms))
    println("== ", id, "  ", m.source)
    for (g, grp) in enumerate(m.groups), j in grp
        println("   g", g, " ", Pop.step_text(m, m.steps[j]), "   ", m.step_names[j])
    end
    println("   fitted: ", m.fitted_params)
    println("   pkg deps: ", m.dependents)
    sol = Thermo.solve(m)
    for (n, (q, ts)) in Thermo.dependent_table(sol)
        println("   joint dep: ", n, " = Keq^", q, " * ", ts)
    end
end
