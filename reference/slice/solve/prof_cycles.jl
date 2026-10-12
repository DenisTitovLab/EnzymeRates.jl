# ABOUTME: Times the parts of `Slice.cycles` on one mechanism to find its slow step.
# ABOUTME: Profiling aid for the solve stage only.
const SCRATCH = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCRATCH, "slice", "src", "Slice.jl"))
using Serialization, Printf, Profile
const S = Slice
m = deserialize(joinpath(SCRATCH, "slice", "core", "data", "init_bibi.jls"))[1].mech
fwd = S.orient(m)
t(f) = (f(); n = 2000; t0 = time_ns(); for _ in 1:n; f(); end; (time_ns() - t0) / n / 1e3)
E = length(m.edges)
@printf("isassigned on Int32 index: %.3f µs\n", t(() -> isassigned(Vector{Vector{Int}}(undef, 5), Int32(2))))
@printf("count_at: %.3f µs\n", t(() -> S.count_at(m.edges[1].cons, 3)))
@printf("names: %.3f µs\n", t(() -> unique!([x.name for x in m.mets])))
@printf("cycles: %.3f µs\n", t(() -> S.cycles(m, fwd)))
Profile.clear()
@profile for _ in 1:20000; S.cycles(m, fwd); end
Profile.print(maxdepth = 14, mincount = 50, noisefloor = 2.0)
