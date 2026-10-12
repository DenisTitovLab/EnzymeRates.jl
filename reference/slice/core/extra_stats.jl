# ABOUTME: Extra counts for the core report: orientation tiers on package-built records, isbits
# ABOUTME: edges and forms, table sizes, the largest Sigs and kinetic groups that tiers 2-3 tie.
include(joinpath(@__DIR__, "..", "src", "Slice.jl"))
using Serialization, Printf
const S = Slice
const HANDBUILT = ("p1_direct", "p2_hand", "p2_fuzz_dump", "p2_fuzz_dump_flip")
rows = [split(l, '\t') for l in readlines(joinpath(@__DIR__, "data", "MANIFEST.tsv"))[2:end]]
tiers, tiers_hand, tiedgroups, maxmets, maxcopies, sigs = zeros(Int, 4), zeros(Int, 4), 0, 0, 0, Pair{Int, String}[]
for r in rows, c in deserialize(joinpath(@__DIR__, "data", r[1] * ".jls"))
    m = c.mech
    _, tier = S.orientation(m)
    foreach(t -> (r[1] in HANDBUILT ? tiers_hand : tiers)[t + 1] += 1, tier)
    for g in 1:S.ngroups(m)
        count(k -> tier[k] > 1, S.members(m, g)) >= 2 && (global tiedgroups += 1)
    end
    global maxmets = max(maxmets, length(m.mets))
    for f in m.forms, i in eachindex(m.mets)
        global maxcopies = max(maxcopies, S.count_at(f.bound, i))
    end
    push!(sigs, length(S.sig(m)) => "$(r[1])/$(c.index)")
end
@printf("tiers 0/1/2/3 on package-built records %s, on hand-built %s\n", tiers, tiers_hand)
@printf("kinetic groups with two or more members decided by tier 2 or 3: %d\n", tiedgroups)
@printf("largest table %d entries; most copies of one entry in a form %d\n", maxmets, maxcopies)
@printf("isbits: Form %s, Edge %s; sizeof Form %d, Edge %d\n", isbitstype(S.Form),
        isbitstype(S.Edge), sizeof(S.Form), sizeof(S.Edge))
sort!(sigs; by = first, rev = true)
println("longest Sigs: ", first(sigs, 5))
