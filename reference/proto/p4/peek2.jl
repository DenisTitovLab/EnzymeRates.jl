# ABOUTME: Prints species and steps (with roles) of chosen mechanisms, plus Thermo columns and dependents.
# ABOUTME: Throwaway P4 helper; stdlib only.
include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))
const T = Thermo
pop = Pop.load_population(ARGS[1])
for i in parse.(Int, ARGS[2:end])
    m = pop[i]
    println("== ", m.id, " ", m.source, " fam=", m.family, " subs=", m.substrates, " prods=", m.products, " inh=", m.inhibitors, " mets=", m.metabolites)
    for (k, sp) in enumerate(m.species)
        println("  form $k: ", sp.name, " conf=", sp.conformation, " bound=", [(b.name, b.role) for b in sp.bound], " +res=", sp.res_added, " -res=", sp.res_subtracted)
    end
    for (j, s) in enumerate(m.steps)
        println("  step $j: ", Pop.step_text(m, s), "  cons=", [(x.name, x.role) for x in s.consumed], " rel=", [(x.name, x.role) for x in s.released], " names=", m.step_names[j])
    end
    println("  groups=", m.groups)
    sol = T.solve(m)
    for d in sol.deps
        println("  dep ", sol.cols[d.col].name, " = Keq^", d.keq, " * ",
                join([string(sol.cols[f].name, "^", e) for (f, e) in d.terms], " * "), "   ", d.label, d.fold ? " fold" : "")
    end
    println("  indep: ", [sol.cols[i].name for i in sol.indep])
end
