# ABOUTME: Prints the keys and a few values of P5a's and the skeptic's limit-study records (B1 follow-up).
# ABOUTME: Read-only peek at the serialized formats; stdlib only.
const S = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(S, "proto/common/Population.jl"))
using Serialization
for f in ("proto/p5a/data/pop_allo_bibi.jls", "proto/p5a/data/pop_expand_bibi.jls",
          "proto/p5a/data/results_allo_uniuni.jls", "proto/p5a/data/results_cases.jls",
          "proto/p5a/data/results_allo_bibi_KV.jls",
          "proto/skeptic_P5a-onlyA/data/res_allo_uniuni_inh.jls", "proto/skeptic_P5a-onlyA/data/res_allo_bibi_pp.jls")
    recs = deserialize(joinpath(S, f))
    println("== ", f, ": ", typeof(recs), " n=", length(recs))
    convs = unique(r[:conv] for r in recs)
    println("   convs ", convs, "  ids e.g. ", [r[:id] for r in recs[1:min(3, end)]])
    ks = union((keys(r) for r in recs)...)
    println("   keys ", sort!(collect(ks)))
    r = recs[1]
    for k in (:id, :conv, :cone_spread, :cone_vs_del, :cone_seq, :admitted, :tags, :cone_dim, :cone_n)
        haskey(r, k) && println("   ", k, " = ", r[k])
    end
    if occursin("cases", f)
        println("   ids: ", unique(r[:id] for r in recs))
    end
end
for name in ("allo_uniuni_inh", "allo_unibi", "allo_bibi_pp", "allo_terbi")
    pop = deserialize(joinpath(S, "proto/skeptic_P5a-onlyA/data/$name.jls"))
    println(name, ": ", typeof(pop), " n=", length(pop), " allosteric ", count(m -> m.allosteric, pop),
            " OnlyA ", count(m -> :OnlyA in m.group_tags, pop), " ok ", count(m -> m.ok, pop))
end
for name in ("allo_bibi", "expand_bibi", "allo_uniuni", "specs", "handwritten")
    pop = Pop.load_population(name)
    al = filter(m -> m.allosteric, pop)
    println(name, ": n=", length(pop), " allosteric ", length(al), " OnlyA ", count(m -> :OnlyA in m.group_tags, al),
            " ids(OnlyA) ", name in ("specs", "handwritten") ? [m.id for m in al if :OnlyA in m.group_tags] : "")
end
