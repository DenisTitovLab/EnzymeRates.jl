# ABOUTME: Lists the multi-block plain mechanisms that have no speed-free (dead-end) step.
# ABOUTME: Explains what the extra blocks are in those mechanisms.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
using Serialization
rows = deserialize(joinpath(@__DIR__, "results.jls"))
pops = Dict(p => Dict(m.id => m for m in Pop.load_population(p)) for p in unique(r[:pop] for r in rows))
multi = [r for r in rows if r[:blocks] >= 2]
println("multi-block mechanisms: ", length(multi), "; with a speed-free step: ",
        count(r -> !isempty(r[:speed_free_thermo]), multi))
for r in multi
    isempty(r[:speed_free_thermo]) || continue
    m = pops[r[:pop]][r[:id]]
    println(r[:id], " blocks ", r[:block_sizes], " parallel ", r[:parallel], "  ", r[:source])
    for s in m.steps
        s.re || println("     SS ", Pop.step_text(m, s))
    end
end
