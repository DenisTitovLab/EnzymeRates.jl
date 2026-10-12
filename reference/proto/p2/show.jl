# ABOUTME: Prints one plain mechanism with my canonical steps, priorities and my solve beside the package's.
# ABOUTME: Usage: julia show.jl <population> <index...>; throwaway P2 probe helper.
include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))
const T = Thermo
pop = Pop.load_population(ARGS[1])
for i in parse.(Int, ARGS[2:end])
    m = pop[i]
    sol = T.solve(m)
    z = T.zero_flux_ss(m)
    println("== ", m.id, " ", m.source, " fam=", m.family, " subs=", m.substrates, " prods=", m.products)
    g = Pop.group_of_step(m)
    for (j, s) in enumerate(m.steps)
        c = sol.cs[j]
        println("  $j [g$(g[j])] ", Pop.step_text(m, s), "  kind=$(c.kind) occ=$(c.occ) fwd=$(c.fwd) flipped=$(c.flipped)", z[j] ? "  FLUX-FREE" : "")
    end
    println("  priority order (eliminated first):")
    for (i, c) in enumerate(sol.cols)
        println("    $i ", c.prio, " ", c.name, i in sol.indep ? "   [indep]" : "")
    end
    for d in sol.deps
        println("  dep ", sol.cols[d.col].name, " = Keq^", d.keq, " * ",
                join([string(sol.cols[f].name, "^", e) for (f, e) in d.terms], " * "), "   ", d.label, d.fold ? " fold" : "")
    end
    println("  pkg fitted: ", m.fitted_params)
    for d in m.dependents; println("  pkg dep ", d); end
end
