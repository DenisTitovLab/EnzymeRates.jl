# ABOUTME: Verifier run B (package + slice): R5.5 round-2 V×τ recounted on the slice's graph, then
# ABOUTME: 200 ter-ter init mechanisms outside the unit study's 2,000 written as slice raw data.
using EnzymeRates, Random, Printf, SHA, TOML
const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "phase45/research5/verify")
include(joinpath(VD, "pkgside.jl"))
include(joinpath(SP, "phase3/census/census_common.jl"))
include(joinpath(SP, "phase3/census/neutral_steps.jl"))
T0 = time()
include(joinpath(SP, "slice/src/Slice.jl"))
include(joinpath(VD, "sliceside.jl"))
const S = Slice
rxnI = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    competitive_inhibitors: A, I
end
function r55_slice()
    seeds = unique!(collect(ER.init_mechanisms(rxnI)))
    base, _ = ER._base_tier(seeds, rxnI)
    seen = Set{UInt64}(hash.(base))
    next1 = Any[]
    for p in [m for m in base if cxA(m) <= 337], c in ER.expand_mechanisms([p], rxnI)
        h = hash(c)
        h in seen && continue
        push!(seen, h)
        cxA(c) <= 337 && push!(next1, c)
    end
    parents = next1[sort(randperm(Xoshiro(20261014), length(next1))[1:1500])]
    n = bad = over_pkg = over_slice = fails = 0
    for p in parents, c in ER.expand_mechanisms([p], rxnI)
        h = hash(c)
        h in seen && continue
        push!(seen, h)
        c isa ER.AllostericMechanism || continue
        a = cxA(c)
        a <= 337 || continue
        i = cxI(c)
        try
            sm = mkmech(S, rawrec(c))
            im = S.inactive_mask(sm)
            sa, si = Int(S.complexity(sm)), any(im) ? Int(S.complexity(sm, im)) : 0
            n += 1
            (sa, si) == (a, i) || (bad += 1; bad <= 5 && println("  mismatch: package ($a, $i) slice ($sa, $si)"))
            over_pkg += a + i > 337
            over_slice += sa + si > 337
        catch e
            fails += 1
            fails <= 3 && println("  slice failed: ", first(sprint(showerror, e), 150))
        end
    end
    @printf("r55 slice recount, round 2 (study seed): %d allosteric within 337 compared, %d V×τ mismatches, %d failures; newly over: package %d, slice %d (%.0f s)\n",
            n, bad, fails, over_pkg, over_slice, time() - T0)
    flush(stdout)
end
r55_slice()
GC.gc()
rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end
ms = ER.init_mechanisms(rxn)
@printf("init_mechanisms(ter_ter): %d in %.0f s, maxrss %.0f MB\n", length(ms), time() - T0, Sys.maxrss() / 2^20)
theirs_idx = Set(sort(randperm(Xoshiro(20261014), length(ms))[1:2000]))
theirs_id = Set(first(split(l, '\t')) for l in eachline(joinpath(SP, "phase45/research5/unit/data/terter.tsv")))
cand = [i for i in eachindex(ms) if !(i in theirs_idx)]
pick = cand[sort(randperm(Xoshiro(99), length(cand))[1:260])]
ids_theirs = Set(mech_id(census_text(ms[i])) for i in sort(collect(theirs_idx)))
@printf("study pick reproduced: %d of %d census ids match unit/data/terter.tsv\n",
        length(intersect(ids_theirs, theirs_id)), length(theirs_id))
n = 0
open(joinpath(VD, "data", "terter_new.tsv"), "w") do io
    for i in pick
        m = ms[i]
        id = mech_id(census_text(m))
        id in theirs_id && (println("overlap ", id); continue)
        cxA(m) > 5000 && continue
        println(io, id, '\t', cxA(m), '\t', repr(rawrec(m)))
        global n += 1
        n == 200 && break
    end
end
@printf("wrote %d ter-ter mechanisms (V×τ ≤ 5000) in %.0f s, maxrss %.0f MB\n", n, time() - T0, Sys.maxrss() / 2^20)
