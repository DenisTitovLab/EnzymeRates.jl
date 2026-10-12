# ABOUTME: Skeptic P5a: overlap of ping-pong mechanisms with the non-unique / sequestering (iu) limits in
# ABOUTME: the extra populations and in the probe's expand_bibi; prints one example per class.
include(joinpath(@__DIR__, "MWC.jl"))
using Serialization
for (name, rfile) in (("allo_bibi_pp", "res_allo_bibi_pp"), ("allo_terbi", "res_allo_terbi"),
                      ("expand_bibi", "pop_expand_bibi"))
    pop = name == "expand_bibi" ? Pop.load_population(name) : deserialize(joinpath(@__DIR__, "data", "$name.jls"))
    byid = Dict(m.id => m for m in pop)
    recs = filter(r -> r[:conv] === :iu && haskey(r, :cone_v), deserialize(joinpath(@__DIR__, "data", "$rfile.jls")))
    pp = Set(m.id for m in pop if :pingpong in m.family && m.allosteric && :OnlyA in m.group_tags)
    bad = Set(r[:id] for r in recs if r[:cone_spread] >= 1e-20 || r[:cone_vs_del] >= 1e-20 || r[:cone_seq])
    println(rpad(name, 13), " :OnlyA ping-pong ", length(pp), "; not (unique, no seq, = deletion) under iu ", length(bad),
            "; both ", length(intersect(pp, bad)), "; ping-pong that pass ", length(setdiff(pp, bad)),
            "; non-ping-pong that fail ", length(setdiff(bad, pp)))
    for id in sort!(collect(setdiff(pp, bad)))[1:min(end, 1)]
        m = byid[id]; g = Pop.group_of_step(m)
        println("   passing ping-pong example ", id)
        for (j, s) in enumerate(m.steps); println("      g", g[j], " [", m.group_tags[g[j]], "] ", Pop.step_text(m, s)); end
    end
    for id in sort!(collect(setdiff(bad, pp)))[1:min(end, 1)]
        m = byid[id]; g = Pop.group_of_step(m)
        println("   failing non-ping-pong example ", id)
        for (j, s) in enumerate(m.steps); println("      g", g[j], " [", m.group_tags[g[j]], "] ", Pop.step_text(m, s)); end
    end
end
