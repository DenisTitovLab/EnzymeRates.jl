# ABOUTME: Prints the stored steps, groups, Full parameters and names of chosen population mechanisms.
# ABOUTME: Usage: julia peek.jl pop/index ...
const C = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common"
include(joinpath(C, "Population.jl"))
cache = Dict{String, Vector{Pop.PMech}}()
for id in ARGS
    pop, i = split(id, "/")
    m = get!(() -> Pop.load_population(pop), cache, pop)[parse(Int, i)]
    println("== ", m.id, "  ", m.source, "  substrates=", m.substrates, " products=", m.products,
            " inhibitors=", m.inhibitors, " family=", m.family)
    g = Pop.group_of_step(m)
    for (j, s) in enumerate(m.steps)
        println("  [", j, "] g", g[j], "  ", Pop.step_text(m, s), "    ", m.step_names[j])
    end
    println("  full_params = ", m.full_params)
    println("  fitted      = ", m.fitted_params)
    println("  dependents  = ", m.dependents)
    println("  RE cycles (total, open) = ", Pop.re_cycle_closure(m))
end
