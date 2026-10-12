# ABOUTME: Classifies the population mechanisms whose :OnlyA limit is not unique or not deletion (P5a).
# ABOUTME: Usage: julia pop_analyze.jl <population> [conv] [n_examples]
include(joinpath(@__DIR__, "MWC.jl"))
using Serialization, Printf
name = ARGS[1]; conv = Symbol(get(ARGS, 2, "i")); nex = parse(Int, get(ARGS, 3, "3"))
recs = filter(r -> r[:conv] === conv && !haskey(r, :error), deserialize(joinpath(@__DIR__, "data", "pop_$name.jls")))
empty_ids = [r[:id] for r in recs if !haskey(r, :cone_v)]
println("cone empty: ", empty_ids)
pop = Dict(m.id => m for m in Pop.load_population(name))
"""Structural features of the :OnlyA groups of a mechanism."""
function features(m)
    tg = MWC.step_tags(m); ch = MWC.is_chem(m)
    act = [t !== :OnlyA for t in tg]
    keep = MWC.reachable(m, act)
    pruned = setdiff(1:length(m.species), keep)
    # pruned forms that are not the bound side of an :OnlyA plain binding from a kept form
    bound_side = Set(s.to for (j, s) in enumerate(m.steps) if tg[j] === :OnlyA && !ch[j] && s.from in keep)
    island = setdiff(pruned, bound_side)
    (pingpong = :pingpong in m.family,
     onlyA_re_chem = any(j -> tg[j] === :OnlyA && ch[j] && m.steps[j].re, eachindex(tg)),
     onlyA_ss_bind = any(j -> tg[j] === :OnlyA && !ch[j] && !m.steps[j].re, eachindex(tg)),
     mixed_group = any(g -> m.group_tags[g] === :OnlyA && any(ch[j] for j in m.groups[g]) &&
                        any(j -> MWC.Thermo.canonical_steps(m)[j].kind == 1, m.groups[g]), eachindex(m.groups)),
     dead_end = :dead_end_leaf in m.family,
     npruned = length(pruned), island = length(island))
end
classes = Dict{String, Vector{String}}()
for r in recs
    haskey(r, :cone_v) || continue
    cls = r[:cone_spread] >= 1e-20 ? (r[:cone_seq] ? "non-unique+seq" : "non-unique") :
          r[:cone_vs_del] >= 1e-20 ? (r[:cone_seq] ? "unique≠del+seq" : "unique≠del") : "unique=del"
    push!(get!(classes, cls, String[]), r[:id])
end
for (cls, ids) in sort!(collect(classes); by = first)
    fs = [features(pop[id]) for id in ids]
    @printf("%-16s %4d | pingpong %3d | :OnlyA RE chem %3d | :OnlyA SS binding %3d | mixed group %3d | dead-end leaf %3d | pruned forms>0 %3d | island (pruned, not an :OnlyA bound side) %3d\n",
            cls, length(ids), count(f -> f.pingpong, fs), count(f -> f.onlyA_re_chem, fs),
            count(f -> f.onlyA_ss_bind, fs), count(f -> f.mixed_group, fs), count(f -> f.dead_end, fs),
            count(f -> f.npruned > 0, fs), count(f -> f.island > 0, fs))
end
for cls in ("non-unique", "non-unique+seq", "unique≠del", "unique≠del+seq")
    haskey(classes, cls) || continue
    for id in classes[cls][1:min(end, nex)]
        m = pop[id]; gos = Pop.group_of_step(m); ch = MWC.is_chem(m)
        r = only(filter(x -> x[:id] == id, recs))
        println("── $cls example $id  features ", features(m))
        for (j, s) in enumerate(m.steps)
            println("     g", gos[j], " [", m.group_tags[gos[j]], "] ", ch[j] ? "chem " : "bind ", Pop.step_text(m, s))
        end
    end
end
