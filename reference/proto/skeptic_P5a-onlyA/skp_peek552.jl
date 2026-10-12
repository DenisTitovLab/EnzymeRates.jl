# ABOUTME: Skeptic P5a: prints expand_bibi/552 and its (i)/(iu) light-study records.
# ABOUTME: Stdlib only.
include(joinpath(@__DIR__, "MWC.jl"))
using Serialization
pop = Pop.load_population("expand_bibi")
recs = deserialize(joinpath(@__DIR__, "data", "pop_expand_bibi.jls"))
for id in ("expand_bibi/552", "expand_bibi/559")
    m = only(filter(x -> x.id == id, pop))
    g = Pop.group_of_step(m); kind = [c.kind for c in MWC.Thermo.canonical_steps(m)]
    println("── ", id)
    for (j, s) in enumerate(m.steps)
        println("   g", g[j], " [", m.group_tags[g[j]], "] kind", kind[j], " ", Pop.step_text(m, s))
    end
    for r in filter(r -> r[:id] == id, recs)
        println("   conv ", r[:conv], ": cone_dim ", r[:cone_dim], " n ", r[:cone_n], " spread ", get(r, :cone_spread, "-"),
                " vs del ", get(r, :cone_vs_del, "-"), " seq ", get(r, :cone_seq, "-"))
    end
end
