# ABOUTME: Characterizes the P8 merges: for each merged member, what differs from its class
# ABOUTME: representative (allosteric tags only, grouping, steps), with a few printed examples.
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
using Serialization
steptexts(m) = sort([Pop.step_text(m, s) for s in m.steps])
groupsets(m) = sort([sort([Pop.step_text(m, m.steps[j]) for j in g]) for g in m.groups])
tagmap(m) = isempty(m.group_tags) ? Dict() :
    Dict(sort([Pop.step_text(m, m.steps[j]) for j in g]) => m.group_tags[k] for (k, g) in enumerate(m.groups))
for (tag, pops) in (("init_bibi", ["init_bibi"]), ("init_bibi_pp+expand_bibi", ["init_bibi_pp", "expand_bibi"]),
                    ("allo_bibi", ["allo_bibi"]))
    ms = reduce(vcat, [Pop.load_population(p) for p in pops])
    byid = Dict(m.id => m for m in ms)
    d = deserialize(joinpath(@__DIR__, "data", "$tag.jls"))
    for variant in (:today, :redesign)
        cls = Dict{String, Vector{String}}()
        for r in d.recs
            k = get(r, (variant, :k2c, :key), nothing)
            k === nothing || push!(get!(cls, k, String[]), r[:id])
        end
        counts = Dict{Symbol, Int}(); tagchanges = Dict{Any, Int}(); ex = Dict{Symbol, Vector{String}}()
        for v in values(cls)
            length(v) < 2 && continue
            r = byid[v[1]]
            for id in v[2:end]
                m = byid[id]
                kind = steptexts(m) != steptexts(r) ? :steps_differ :
                       groupsets(m) != groupsets(r) ? :grouping_differs :
                       r.cat_multiplicity != m.cat_multiplicity || r.sites != m.sites ? :sites_differ :
                       :tags_only
                counts[kind] = get(counts, kind, 0) + 1
                if kind === :tags_only
                    tr, tm = tagmap(r), tagmap(m)
                    ch = sort([(tr[g], tm[g]) for g in keys(tr) if tr[g] != tm[g]])
                    tagchanges[ch] = get(tagchanges, ch, 0) + 1
                end
                e = get!(ex, kind, String[])
                length(e) < 2 && push!(e, string(r.id, " ~ ", m.id))
            end
        end
        println(tag, " / ", variant, " (k2c): ", counts)
        for (ch, c) in sort!(collect(tagchanges); by = x -> -x[2])
            println("    tag changes ", ch, " => ", c)
        end
        for (k, e) in ex
            println("    e.g. ", k, ": ", e)
        end
    end
end
# one detailed MWC example
ms = Pop.load_population("allo_bibi"); byid = Dict(m.id => m for m in ms)
for id in ("allo_bibi/682", "allo_bibi/690")
    m = byid[id]
    println("== ", id, " ", m.source, " fitted ", m.fitted_params)
    for (k, g) in enumerate(m.groups), j in g
        println("   g", k, " ", m.group_tags[k], " ", Pop.step_text(m, m.steps[j]), "  kept in I: ",
                m.inactive_kept[j], "  A ", m.step_names_A[j], " I ", m.step_names_I[j])
    end
end
