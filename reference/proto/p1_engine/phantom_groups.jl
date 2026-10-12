# ABOUTME: Classifies speed-free SS steps by kinetic group: a group is speed-free only when every step is.
# ABOUTME: Compares that structural prediction with the package-level phantoms found by extras.jl.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
using Serialization
rows = Dict(r[:id] => r for r in deserialize(joinpath(@__DIR__, "results.jls")))
for pop in ["specs", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi", "handwritten"]
    ms = [m for m in Pop.load_population(pop) if !m.allosteric]
    shared = alone = grp_free = 0; examples = String[]
    for m in ms
        free = Set(rows[m.id][:speed_free_thermo])
        isempty(free) && continue
        gof = Pop.group_of_step(m)
        for j in free
            length(m.groups[gof[j]]) > 1 ? (shared += 1) : (alone += 1)
        end
        for g in unique(gof[j] for j in free)
            if all(j -> j in free, m.groups[g])
                grp_free += 1
                length(examples) < 3 && push!(examples, "$(m.id): group $g = steps $(m.groups[g]) " *
                    join([Pop.step_text(m, m.steps[j]) for j in m.groups[g]], "; "))
            end
        end
    end
    println(rpad(pop, 13), " speed-free steps in a multi-step group ", shared, ", alone in their group ", alone,
            "; groups whose every step is speed-free ", grp_free)
    foreach(e -> println("    ", e), examples)
end
m = only(filter(m -> m.id == "specs/14", Pop.load_population("specs")))
println("specs/14 steps: ", [Pop.step_text(m, s) for s in m.steps])
