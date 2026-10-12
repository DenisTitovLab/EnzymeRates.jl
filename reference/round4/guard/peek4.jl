# ABOUTME: Prints every package verdict and v line of a run (small runs only).
using Serialization
const R = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard"
for r in deserialize(joinpath(R, "data", "pkg_$(ARGS[1]).jls"))
    println(r.key, " | ctor=", something(r.ctor, "-"), " | derive=", something(r.derive, "-"), " | ", first(r.vline, 140))
end
