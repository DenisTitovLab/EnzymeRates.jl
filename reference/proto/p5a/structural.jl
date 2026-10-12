# ABOUTME: A cheap structural test for "deletion equals the true :OnlyA limit" (P5a), checked
# ABOUTME: against the cone verdicts of the population study for conventions (i) and (iu).
include(joinpath(@__DIR__, "MWC.jl"))
using Serialization, Printf

"""Classify the inactive state's pruned forms. Components of pruned forms joined by
non-:OnlyA steps move together; a component is `:safe` when every :OnlyA step on its
boundary is a binding whose bound form lies inside it (its weight vanishes like an
association), `:sequester` when some :OnlyA binding has its unbound form inside and its
bound form outside (dissociation drives mass into it), and `:free` otherwise (only
chemistry on its boundary: its mass follows the vanishing chemistry rates).
`per_group` classifies chemistry per kinetic group (convention i) or per step (iu)."""
function pruned_components(m; per_group)
    tg = MWC.step_tags(m); ch = MWC.is_chem(m; per_group)
    act = [t !== :OnlyA for t in tg]
    keep = MWC.reachable(m, act)
    U = setdiff(1:length(m.species), keep)
    comp = Dict(u => u for u in U)
    root(u) = comp[u] == u ? u : (comp[u] = root(comp[u]))
    for (j, s) in enumerate(m.steps)
        act[j] && s.from in U && s.to in U && (comp[root(s.from)] = root(s.to))
    end
    verdict = Dict{Int, Symbol}()
    for c in unique(root(u) for u in U)
        inside(x) = x in U && root(x) == c
        v = :safe
        for (j, s) in enumerate(m.steps)
            act[j] && continue
            a, b = inside(s.from), inside(s.to)
            a == b && continue                       # not on this component's boundary
            if ch[j]
                v === :sequester || (v = :free)
            elseif a && !b                           # unbound form inside, bound outside
                v = :sequester
            end
        end
        verdict[c] = v
    end
    collect(values(verdict))
end

for name in ("expand_bibi", "allo_bibi")
    pop = Dict(m.id => m for m in Pop.load_population(name))
    recs = deserialize(joinpath(@__DIR__, "data", "pop_$name.jls"))
    for conv in (:i, :iu)
        rs = filter(r -> r[:conv] === conv && haskey(r, :cone_v), recs)
        tab = Dict{Tuple{Symbol, String}, Int}()
        for r in rs
            vs = pruned_components(pop[r[:id]]; per_group = conv === :i)
            pred = any(==(:sequester), vs) ? :sequester : any(==(:free), vs) ? :free : :safe
            truth = r[:cone_spread] >= 1e-20 ? (r[:cone_seq] ? "nonunique+seq" : "nonunique") :
                    r[:cone_vs_del] >= 1e-20 ? (r[:cone_seq] ? "unique≠del+seq" : "unique≠del") : "unique=del"
            tab[(pred, truth)] = get(tab, (pred, truth), 0) + 1
        end
        println("== $name conv ($conv): structural verdict × cone verdict")
        for ((p, t), n) in sort!(collect(tab); by = x -> (string(x[1][1]), x[1][2]))
            @printf("   %-10s %-16s %5d\n", p, t, n)
        end
    end
end
