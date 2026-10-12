# ABOUTME: Smoke test of the solve stage: solves a few converted mechanisms and prints their
# ABOUTME: columns, dependents and parameter lists next to the package's.
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
const DATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/core/data"
include(joinpath(PROTO, "common", "Population.jl"))
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
using Serialization
const S = Slice
conv = deserialize(joinpath(DATA, "specs.jls"))
pops = Pop.load_population("specs")
function show_sol(c, p)
    m = c.mech
    sol = S.solve(m)
    names = S.column_names(m, sol)
    println("== ", p.id, " ", p.source)
    for (k, e) in enumerate(m.edges)
        println("   ", S.step_text(m, S.oriented(m, k, S.orient(m))), "  g", m.group[k],
                S.is_allosteric(m) ? " $(m.tags[m.group[k]])" : "")
    end
    for d in sol.deps
        rhs = join(["$(names[i])^$(e)" for (i, e) in d.terms], " * ")
        println("   dep ", names[d.col], sol.cols[d.col].inactive ? " [I]" : "", " = Keq^", d.keq, " * ", rhs, "  [", d.label, "]")
    end
    println("   fold: ", [names[a] => names[b] for (a, b) in sol.fold])
    println("   fitted: ", S.fitted_params(m))
    println("   pkg fitted: ", p.fitted_params)
    println("   pkg deps: ", p.dependents)
end
for i in [1, 5, 12, 30, 35, 40]
    show_sol(conv[i], pops[conv[i].index])
end
m = conv[5].mech
pars = S.parameters(m)
vals = (; (n => rand() + 0.5 for n in pars)...)
println(S.all_parameters(m, vals))
