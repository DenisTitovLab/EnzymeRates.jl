# ABOUTME: Enumeration and admission: runs the slice's three ported moves on the admitted allosteric
# ABOUTME: mechanisms of the enumerated sets and counts the children that admission rejects, by kind.
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCRATCH, "proto", "common", "Population.jl"))
module Slice
for file in ("types.jl", "names.jl", "poly.jl", "convert.jl", "thermo.jl", "enum.jl",
             "admission.jl")
    include(joinpath(Main.SCRATCH, "slice", "src", file))
end
end
using Serialization, Printf
const S = Slice

kind(msg) = msg === nothing ? "admitted" : startswith(msg, "the conformation flips") ? "flip" :
    occursin("still join", msg) ? "G1 chord" : occursin("close a cycle", msg) ? "G1 cycle" :
    "islands"
const SETS = ("allo_uniuni", "allo_uniuni_inh", "allo_unibi", "allo_bibi", "expand_bibi",
              "allo_bibi_pp", "allo_terbi")
tab, errs, examples = Dict{Any, Int}(), Dict{Any, Int}(), Dict{Any, String}()
t = @elapsed for set in SETS
    conv = deserialize(joinpath(SCRATCH, "slice", "core", "data", set * ".jls"))
    parents = [c.mech for c in conv if S.is_allosteric(c.mech) && S.admit(c.mech) === nothing]
    parents = parents[1:min(end, 300)]
    for m in parents
        declared = unique([:I; m.substrates[1]; m.products[end]])
        for (move, f) in ((:flip, S.flip_children), (:split, S.split_children),
                          (:dead_end, x -> S.dead_end_children(x, declared)))
            children = try
                f(m)
            catch err
                k = (set, move, first(sprint(showerror, err), 60))
                errs[k] = get(errs, k, 0) + 1
                continue
            end
            for c in children
                msg = S.admit(c)
                k = (set, move, kind(msg))
                tab[k] = get(tab, k, 0) + 1
                msg === nothing || haskey(examples, (move, kind(msg))) ||
                    (examples[(move, kind(msg))] = msg)
            end
        end
    end
    @printf("%-16s %d admitted parents (at most 300)\n", set, length(parents))
end
println("children by (set, move, admission) => count:")
foreach(((k, v),) -> println("    ", k, " => ", v), sort!(collect(tab); by = string))
println("move errors: ", isempty(errs) ? "none" : errs)
println("example rejections:")
foreach(((k, v),) -> println("    ", k, ": ", v), examples)
@printf("time %.0f s\n", t)
