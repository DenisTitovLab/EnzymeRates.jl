# ABOUTME: Prints two laws whose rate is not positive at zero products and zero modifiers, with
# ABOUTME: their rate at random substrate points with the modifiers at zero and at 1.
include(joinpath(@__DIR__, "klib.jl"))
for (set, idx) in (("specs", 24), ("p2_adv_pkg", 14))
    r = only(filter(r -> r.c.index == idx, records_of(set)))
    m = r.c.mech
    println(r.p.id, " (", set, "/", idx, ") modifiers ", [n for n in S.metabolites(m) if !S.is_reactant(m, n)],
            " sites ", m.sites)
    foreach(l -> println("    ", l), step_lines(m))
    _, _, sol, L = derivable(m)
    println("    kcat class ", kcat_class(m, sol)[1:2])
    rng = Xoshiro(3)
    x, _ = params_at(m, sol, rng)
    nl = NLaw(m, L, x)
    for mv in (0.0, 1.0)
        c = zeros(nl.nc)
        foreach(i -> c[i] = 1.0, nl.subs)
        foreach(i -> c[i] = mv, nl.mods)
        println("    v at substrates 1, products 0, modifiers ", mv, ": ", rate(nl, c))
    end
end
