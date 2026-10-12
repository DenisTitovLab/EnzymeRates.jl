# ABOUTME: For tag-only merges (EqualAI vs OnlyA), checks whether every retagged group is absent
# ABOUTME: from the inactive state in both mechanisms (then the tag cannot matter).
include(joinpath(@__DIR__, "..", "common", "Population.jl"))
isdefined(Main, :Allo) || include(joinpath(@__DIR__, "vendor", "Allo_p5b.jl"))
using Serialization
gkey(m, g) = sort([Pop.step_text(m, m.steps[j]) for j in g])
for (tag, pops) in (("init_bibi_pp+expand_bibi", ["init_bibi_pp", "expand_bibi"]), ("allo_bibi", ["allo_bibi"]))
    ms = reduce(vcat, [Pop.load_population(p) for p in pops]); byid = Dict(m.id => m for m in ms)
    d = deserialize(joinpath(@__DIR__, "data", "$tag.jls"))
    cls = Dict{String, Vector{String}}()
    for r in d.recs
        k = get(r, (:redesign, :k2c, :key), nothing)
        k === nothing || push!(get!(cls, k, String[]), r[:id])
    end
    stats = Dict{Any, Int}()
    for v in values(cls), id in v[2:end]
        r, m = byid[v[1]], byid[id]
        sort([Pop.step_text(r, s) for s in r.steps]) == sort([Pop.step_text(m, s) for s in m.steps]) || continue
        gr = Dict(gkey(r, g) => k for (k, g) in enumerate(r.groups))
        gm = Dict(gkey(m, g) => k for (k, g) in enumerate(m.groups))
        keys(gr) == keys(gm) || continue
        keptr = Allo.inactive_kept(r); keptm = Allo.inactive_kept(m)
        for (g, kr) in gr
            km = gm[g]
            r.group_tags[kr] == m.group_tags[km] && continue
            pruned_r = !any(keptr[j] for j in r.groups[kr]); pruned_m = !any(keptm[j] for j in m.groups[km])
            change = Tuple(sort([r.group_tags[kr], m.group_tags[km]]))
            stats[(change, pruned_r && pruned_m)] = get(stats, (change, pruned_r && pruned_m), 0) + 1
        end
    end
    println(tag, " (redesign k2c), retagged groups: (tag pair, absent from the inactive state in both) => count")
    for (k, c) in sort!(collect(stats); by = string)
        println("   ", k, " => ", c)
    end
end
