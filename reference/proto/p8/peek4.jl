# ABOUTME: Prints the structures of same-model pairs that are not renamings (from the LM estimate)
# ABOUTME: so the report can say what differs between them.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
pairs = [("allo_bibi", "allo_bibi/1018", "allo_bibi/1019"), ("allo_bibi", "allo_bibi/1004", "allo_bibi/1438"),
         ("expand_bibi", "expand_bibi/556", "expand_bibi/614"), ("expand_bibi", "expand_bibi/304", "expand_bibi/305")]
for (pop, a, b) in pairs
    ms = Pop.load_population(pop); byid = Dict(m.id => m for m in ms)
    for id in (a, b)
        m = byid[id]
        println("== ", id, " ", m.source, " allosteric=", m.allosteric, " fitted ", length(m.fitted_params))
        for (k, g) in enumerate(m.groups), j in g
            println("   g", k, " ", isempty(m.group_tags) ? "" : m.group_tags[k], " ", Pop.step_text(m, m.steps[j]),
                    m.allosteric ? "  kept in I: $(m.inactive_kept[j])" : "")
        end
    end
end
