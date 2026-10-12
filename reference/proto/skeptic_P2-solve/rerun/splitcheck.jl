# ABOUTME: Shows examples of single-group splits that keep the independent count, with the reason.
# ABOUTME: Throwaway P2 probe script; stdlib only.
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common", "Thermo.jl"))
const T = Thermo
regrouped(m, groups) = Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates,
    products = m.products, inhibitors = m.inhibitors, metabolites = m.metabolites,
    species = m.species, steps = m.steps, groups = groups)
shown = 0; bykind = Dict{String, Int}()
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    sol = T.solve(m); n = length(sol.indep)
    for (gi, g) in enumerate(m.groups)
        length(g) < 2 && continue
        groups = [[gg for (k, gg) in enumerate(m.groups) if k != gi]; [[j] for j in g]]
        s2 = T.solve(regrouped(m, groups))
        length(s2.indep) == n || continue
        k = string(m.steps[g[1]].re ? "RE" : "SS", " size ", length(g))
        bykind[k] = get(bykind, k, 0) + 1
        if shown < 4
            global shown += 1
            println("== ", m.id, " group ", gi, ": ", join([Pop.step_text(m, m.steps[j]) for j in g], " | "),
                    "  count ", n, " -> ", length(s2.indep), "  cycles ", sol.ncycles)
            for d in s2.deps
                d.label === :kinetic_group && continue
                println("     ", T.colname(s2, d.col), " ", d.label, d.fold ? " fold" : "", " = Keq^", d.keq, " ",
                        join([string(T.colname(s2, f), "^", e) for (f, e) in d.terms], " "))
            end
        end
    end
end
println(bykind)
