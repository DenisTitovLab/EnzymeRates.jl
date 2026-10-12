# ABOUTME: Prints a few Pop records' package outputs (fitted list order, dependents) to inform
# ABOUTME: the solve stage's parameter order; read-only.
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
include(joinpath(PROTO, "common", "Population.jl"))
using Serialization
specs = Pop.load_population("specs")
for p in specs[[1, 5, 12, 30, 35, 40]]
    println(p.id, "  ", p.source, "  allosteric=", p.allosteric)
    for (j, s) in enumerate(p.steps)
        println("   ", Pop.step_text(p, s), "   ", p.allosteric ? (p.step_names_A[j], p.step_names_I[j]) : p.step_names[j])
    end
    println("   fitted: ", p.fitted_params)
    println("   deps: ", p.dependents)
    println("   rename: ", p.pass2_rename, p.pass2_rename_A, p.pass2_rename_I)
    println("   tags: ", p.group_tags, " sites ", p.sites)
end
