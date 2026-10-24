# ABOUTME: Verifier run C (package): 200 ter-ter init mechanisms drawn uniformly, in random order,
# ABOUTME: from those outside the unit study's 2,000, written as slice raw data (terter_new2.tsv).
using EnzymeRates, Random, Printf, SHA, TOML
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "phase45/research5/verify")
include(joinpath(VD, "pkgside.jl"))
include(joinpath(SP, "phase3/census/census_common.jl"))
include(joinpath(SP, "phase3/census/neutral_steps.jl"))
T0 = time()
rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end
ms = ER.init_mechanisms(rxn)
theirs_idx = Set(sort(randperm(Xoshiro(20261014), length(ms))[1:2000]))
theirs_id = Set(first(split(l, '\t')) for l in eachline(joinpath(SP, "phase45/research5/unit/data/terter.tsv")))
cand = [i for i in eachindex(ms) if !(i in theirs_idx)]
pick = cand[randperm(Xoshiro(2718), length(cand))[1:400]]      # random order, not sorted
n = 0
idx = Int[]
open(joinpath(VD, "data", "terter_new2.tsv"), "w") do io
    for i in pick
        m = ms[i]
        id = mech_id(census_text(m))
        id in theirs_id && continue
        cxA(m) > 5000 && continue
        println(io, id, '\t', cxA(m), '\t', repr(rawrec(m)))
        push!(idx, i)
        global n += 1
        n == 200 && break
    end
end
@printf("wrote %d ter-ter mechanisms; init indices: median %.0f of %d, quartiles %s; %.0f s, maxrss %.0f MB\n",
        n, sort(idx)[100], length(ms), string(sort(idx)[[50, 150]]), time() - T0, Sys.maxrss() / 2^20)
