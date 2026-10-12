# ABOUTME: Skeptic P5a: prints the failing :OnlyA ping-pong mechanisms of allo_bibi_pp with their
# ABOUTME: (iu) verdicts, tag patterns of the chemistry steps, and which forms are pruned.
include(joinpath(@__DIR__, "MWC.jl"))
using Serialization
pop = deserialize(joinpath(@__DIR__, "data", "allo_bibi_pp.jls"))
recs = Dict(r[:id] => r for r in deserialize(joinpath(@__DIR__, "data", "res_allo_bibi_pp.jls")) if r[:conv] === :iu)
pat = Dict{String, Int}()
shown = 0
for m in pop
    :pingpong in m.family || continue
    r = recs[m.id]
    g = Pop.group_of_step(m); kind = [c.kind for c in MWC.Thermo.canonical_steps(m)]
    chem = [string(m.group_tags[g[j]]) for j in eachindex(m.steps) if kind[j] != 1]
    key = join(sort(chem), "+") * (r[:cone_spread] >= 1e-20 ? " nonunique" : "") * (r[:cone_seq] ? " seq" : "")
    pat[key] = get(pat, key, 0) + 1
    if shown < 3
        global shown += 1
        println("── ", m.id, "  spread ", r[:cone_spread], " seq ", r[:cone_seq], " vs del ", r[:cone_vs_del])
        for (j, s) in enumerate(m.steps)
            println("     g", g[j], " [", m.group_tags[g[j]], "] ", kind[j] == 1 ? "bind " : "chem ", Pop.step_text(m, s),
                    m.inactive_kept[j] ? "" : "   (pruned in I)")
        end
    end
end
println("chemistry-step tag patterns of the ping-pong mechanisms (with the iu verdict):")
for (k, v) in sort!(collect(pat); by = first); println("   ", rpad(k, 50), v); end
