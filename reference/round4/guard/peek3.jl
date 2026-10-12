# ABOUTME: Debug: lists the package's non-accept verdicts in a run.
using Serialization
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
pk = deserialize(joinpath(R, "data", "pkg_$(ARGS[1]).jls"))
println("n=", length(pk), " ctor rejects=", count(r -> r.ctor !== nothing, pk), " derive rejects=", count(r -> r.derive !== nothing, pk))
for r in pk
    r.ctor === nothing && r.derive === nothing && continue
    println(r.key, " | ctor=", r.ctor === nothing ? "-" : first(r.ctor, 160), " | derive=", r.derive === nothing ? "-" : first(r.derive, 160))
end
z = [r for r in pk if occursin("* (0)", r.vline) || occursin("(0) /", r.vline)]
println("zero numerators: ", length(z)); foreach(r -> println("  ", r.key, " ", first(r.vline, 150)), z[1:min(end, 5)])
