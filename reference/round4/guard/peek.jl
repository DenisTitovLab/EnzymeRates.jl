# ABOUTME: Prints a few raw slice-side results of a run (debugging aid).
using Serialization
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
out = deserialize(joinpath(R, "data", "slice_$(ARGS[1]).jls"))
for r in out[1:3]
    println(r.key); haskey(r, :law) && println("  law: ", r.law); haskey(r, :num) && println("  num: ", r.num); haskey(r, :struct_) && println("  struct: ", r.struct_)
end
fat = [r for r in out if haskey(r, :fatal)]
println("fatal: ", length(fat)); foreach(r -> println("  ", r.key, " ", r.fatal), fat[1:min(end, 10)])
nerr = [r for r in out if haskey(r, :num) && r.num isa Dict && haskey(r.num, :err)]
println("num errors: ", length(nerr)); foreach(r -> println("  ", r.key, " ", r.num[:err]), nerr[1:min(end, 10)])
for r in out
    haskey(r, :num) && r.num isa Dict && !haskey(r.num, :err) || continue
    vs = [r.num[f].verdict for f in (:pos, :fwd, :rev)]
    all(==(:agree), vs) && continue
    println(r.key, " ", vs, " red=", get(r.law, :reduced, nothing), " ", r.prov)
end
