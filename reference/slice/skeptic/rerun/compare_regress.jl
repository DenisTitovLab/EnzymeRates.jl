# ABOUTME: Compares two stage digests (regress.jl outputs) part by part; usage:
# ABOUTME: compare_regress.jl [a.jls b.jls] (default: before and after the integration edits)
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad", "proto", "common", "Population.jl"))
using Serialization
fa, fb = isempty(ARGS) ? ("regress_before.jls", "regress_after.jls") : (ARGS[1], ARGS[2])
a = deserialize(joinpath(@__DIR__, "data", fa))
b = deserialize(joinpath(@__DIR__, "data", fb))
println("records: before ", length(a), ", after ", length(b), ", common ",
        length(intersect(keys(a), keys(b))))
diffs = Dict{Symbol, Vector{Tuple{String, Int}}}()
counts = Dict{Symbol, Int}()
for (k, pa) in a
    pb = Dict(b[k])
    for (part, h) in pa
        counts[part] = get(counts, part, 0) + 1
        get(pb, part, "") == h || push!(get!(diffs, part, Tuple{String, Int}[]), k)
    end
end
for part in sort!(collect(keys(counts)))
    d = get(diffs, part, [])
    println(rpad(part, 14), counts[part], " compared, ", length(d), " differ",
            isempty(d) ? "" : ": " * join(first(sort(d), 8), ", "))
end
