# ABOUTME: Lists the smallest mechanisms whose Hill quotient carries Keq, from census.jls.
# ABOUTME: Throwaway P4 helper; stdlib only.
using Serialization
include(joinpath(@__DIR__, "P4.jl"))
rows = deserialize(joinpath(@__DIR__, "census.jls"))
bad = sort([r for r in rows if !r[:hill_keqfree]]; by = r -> r[:num_terms])
for r in first(bad, 12)
    println(r[:id], "  num_terms=", r[:num_terms], " haldane=", r[:haldane], " fam=", r[:family])
end
println("not Keq-free by population: ", Dict(p => count(r -> startswith(r[:id], p * "/"), bad) for p in unique(first.(split.(getindex.(rows, :id), "/")))))
println("not Keq-free with exactly one Haldane dependent: ", count(r -> length(r[:haldane]) == 1, bad))
