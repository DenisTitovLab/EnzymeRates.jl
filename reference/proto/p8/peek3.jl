# ABOUTME: Prints two same-model laws that the rename-proof key keeps apart (after factor
# ABOUTME: cancellation): their steps, groups, factors removed and normalized terms.
include(joinpath(@__DIR__, "Laws.jl"))
include(joinpath(@__DIR__, "Factor.jl"))
ms = Pop.load_population("init_bibi")
for id in ARGS
    m = only(filter(x -> x.id == id, ms))
    println("== ", id)
    for (g, grp) in enumerate(m.groups), j in grp
        println("   g", g, " ", Pop.step_text(m, m.steps[j]))
    end
    b = Laws.plain_redesign(m)
    law, removed, _ = Factor.cancel(b.law)
    println("   removed: ", [join(("$(c)*" * join(("$s^$e" for (s, e) in mo), "*") for (mo, c) in f), " + ") for (_, f) in removed])
    short = Dict(p => Symbol("p", i) for (i, p) in enumerate(law.params))
    for (p, s) in sort!(collect(short); by = x -> string(x[2])); println("   ", s, " = ", p); end
    for (i, c) in enumerate(Canon.normalize(Canon.rename(law, short)))
        println("   comp ", i, ": ", join(sort!([string(co, "·", join(("$s^$e" for (s, e) in mo), " ")) for (mo, co) in c]), "  +  "))
    end
end
