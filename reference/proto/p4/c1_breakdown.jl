# ABOUTME: Breaks down where the literal Cha-contracted criterion (c1) disagrees with the polynomials.
# ABOUTME: Throwaway P4 helper reading census.jls.
using Serialization
include(joinpath(@__DIR__, "P4.jl"))
rows = deserialize(joinpath(@__DIR__, "census.jls"))
noc = [r for r in rows if !r[:const_post]]
wrong = [r for r in noc if r[:crit].c1 == 1]
right = [r for r in noc if r[:crit].c1 != 1]
cls(r) = (laurent = r[:crit].no_minimal > 0, all_ss = :all_ss in r[:family], pingpong = :pingpong in r[:family],
          conformations = :conformations in r[:family], tscc = r[:crit].c3)
tally(v) = sort(collect(Dict(k => count(r -> cls(r) == k, v) for k in unique(cls.(v)))); by = string)
println("no constant term: ", length(noc), "; c1 wrongly predicts one: ", length(wrong), "; c1 right: ", length(right))
println("c1 wrong by class: "); foreach(p -> println("  ", p), tally(wrong))
println("c1 right by class: "); foreach(p -> println("  ", p), tally(right))
println("c1 right ids: ", [r[:id] for r in right])
