# ABOUTME: Feature coverage of the validated plain mechanisms (family tags and graph features).
# ABOUTME: Reads results.jls written by validate_engine.jl.
using Serialization
rows = deserialize(joinpath(@__DIR__, "results.jls"))
tags = Dict{Symbol, Int}()
for r in rows, f in r[:family]; tags[f] = get(tags, f, 0) + 1; end
println("family tags over ", length(rows), " plain mechanisms: ",
        join(["$k=$v" for (k, v) in sort!(collect(tags); by = first)], ", "))
println("self-loops (SS step inside an RE segment): ", count(r -> r[:loops] > 0, rows), " mechanisms; ",
        "parallel SS steps: ", count(r -> r[:parallel] > 0, rows), "; RE cycles: ",
        count(r -> !isempty(r[:re_closing]), rows), "; no free form (dfree = 1 by rule): see synthetic")
