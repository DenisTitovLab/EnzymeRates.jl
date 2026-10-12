# ABOUTME: For every G1-rejected candidate, whether each offending :OnlyA plain binding shares its kinetic
# ABOUTME: group with a non-plain step (which makes the package's guard drop the group as chemistry).
include(joinpath(@__DIR__, "vlib.jl"))
for name in ["bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    d, mechs = loadrun(name)
    tab = Dict{Tuple{String, Bool}, Int}()
    for r in d.recs
        (r[:status] === :ok && r[:fitted] && r[:admit] !== nothing) || continue
        kind = occursin("still join", r[:admit]) ? "G1_inside" : occursin("close a cycle", r[:admit]) ? "G1_cycle" : continue
        ms = mstruct(mechs[r[:gid]])
        nf = length(ms.forms)
        kept = [ms.tags[s.g] !== :OnlyA for s in ms.steps]
        p = collect(1:nf)
        f(i) = (while p[i] != i; i = p[i]; end; i)
        foreach(k -> kept[k] && (p[f(ms.steps[k].a)] = f(ms.steps[k].b)), eachindex(ms.steps))
        comp = [f(i) for i in 1:nf]
        arcs = [k for (k, s) in enumerate(ms.steps) if !kept[k] && plain_binding(ms, s)]
        # the bindings the verdict rests on: inside a component, or on a directed cycle
        bad = kind == "G1_inside" ? [k for k in arcs if comp[ms.steps[k].a] == comp[ms.steps[k].b]] : arcs
        mixed = all(k -> any(s -> s.g == ms.steps[k].g && !plain_binding(ms, s), ms.steps), bad)
        key = (kind, mixed)
        tab[key] = get(tab, key, 0) + 1
    end
    println("[$name] (kind, every offending binding shares its group with a non-plain step): ",
            join(("$(k[1])/$(k[2]) $v" for (k, v) in sort!(collect(tab); by = first)), ", "))
end
