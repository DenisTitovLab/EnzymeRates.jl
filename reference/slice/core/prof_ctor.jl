# ABOUTME: Profiles the mechanism constructor on init_bibi raw inputs (allocations and hot
# ABOUTME: frames) to see where its microseconds go; a diagnostic, not part of the slice.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common/Population.jl")
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
using Profile, Printf, Statistics
const S = Slice
raws = [S.pop_raw(p) for p in Pop.load_population("init_bibi")]
f() = for r in raws; S.canonicalize(r); end
f()
@printf("allocations per construction: %.0f, bytes %.0f\n", (@allocations f()) / length(raws),
        (@allocated f()) / length(raws))
t = minimum(@elapsed(f()) for _ in 1:20) / length(raws)
@printf("min over 20 passes: %.0f ns per construction\n", t * 1e9)
Profile.init(n = 10^7, delay = 2e-5); Profile.clear()
@profile for _ in 1:400; f(); end
Profile.print(IOContext(stdout, :displaysize => (200, 160)); format = :flat, sortedby = :count, mincount = 150)
