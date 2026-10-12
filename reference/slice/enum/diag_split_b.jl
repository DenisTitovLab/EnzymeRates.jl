# ABOUTME: Split diagnosis, step 3: rebuilds both split families from per-set gain verdicts (the
# ABOUTME: package's count, the slice's solve, a P5b-style count) and classifies each difference.
include(joinpath(@__DIR__, "load.jl"))
include(joinpath(@__DIR__, "diag_lib.jl"))

meta = deserialize(joinpath(S, "slice/enum/data/diag_split_meta.jls"))
export_ = deserialize(joinpath(S, "slice/enum/data/diag_split_export.jls"))
pkgc = deserialize(joinpath(S, "slice/enum/data/diag_split_pkgcounts.jls"))
issues = deserialize(joinpath(S, "slice/enum/data/diag_split_issues.jls"))
rows = Dict(set => deserialize(joinpath(S, "slice/enum/data/cmp_$set.jls"))
            for set in ("expand_bibi", "allo_bibi"))

at = 0
tally = Dict{String, Int}()
bump(k) = (tally[k] = get(tally, k, 0) + 1)
examples = Dict{String, Vector{String}}()
for mt in meta
    global at
    m = rows[mt.set][mt.row].mech
    units, child = _testhelper_units(m)
    n = length(mt.sets)
    cp = pkgc[at + 1]
    cs = pkgc[at + 2:at + 1 + n]
    fp, fs = issues[at + 1], issues[at + 2:at + 1 + n]
    at += 1 + n
    cands = [child(s) for s in mt.sets]
    mine = [length(Slice.solve(c).indep) for c in cands]
    p5b = [_testhelper_p5b_count(c) for c in cands]
    base_mine, base_p5b = length(Slice.solve(m).indep), _testhelper_p5b_count(m)
    idx = Dict(s => i for (i, s) in enumerate(mt.sets))
    ok(s) = all(u -> isempty(units[u][4]), s) || Slice.bottomless(cands[idx[s]]) === nothing
    fam(cnt, base) = Slice.minimal_sets(_testhelper_partners(m, units)) do s
        haskey(idx, s) || (bump("a set of more than 3 units was tested"); return false)
        ok(s) && cnt[idx[s]] > base
    end
    fpkg, fmine, fp5b = fam(cs, cp), fam(mine, base_mine), fam(p5b, base_p5b)
    Set(idx[s] for s in fpkg) == Set(mt.pkg_sets) || bump("package family not reproduced")
    Set(idx[s] for s in fmine) == Set(mt.mine_sets) || bump("slice family not reproduced")
    # verdict differences on the sets either family tested
    for s in mt.sets
        i = idx[s]
        vp, vm, v5 = cs[i] > cp, mine[i] > base_mine, p5b[i] > base_p5b
        vp == vm && continue
        key = "gain pkg=$vp slice=$vm p5b=$v5" *
              (v5 == vp ? " (slice differs from P5b: thermo's I columns)" :
                          " (package differs from P5b: package count)")
        bump(key)
        bump("  offsets package - slice: parent $(cp - base_mine), child $(cs[i] - mine[i])")
        f = cp - base_mine > cs[i] - mine[i] ? fp : fs[i]
        bump("  excess on the $(cp - base_mine > cs[i] - mine[i] ? "parent" : "child"): " *
             (f.ki1 && f.ki2 ? "known issues 1 and 2" : f.ki1 ? "known issue 1" :
              f.ki2 ? "known issue 2" : "neither flag"))
        ex = get!(examples, key, String[])
        length(ex) < 3 && push!(ex, "$(mt.id) set $(s): Δpkg $(cs[i] - cp), Δslice " *
                                    "$(mine[i] - base_mine), Δp5b $(p5b[i] - base_p5b)")
    end
    bump("parent offset package - slice = $(cp - base_mine), package - P5b = $(cp - base_p5b)")
    fp5b == fpkg ? bump("P5b family = package family") : fp5b == fmine ?
        bump("P5b family = slice family") : bump("P5b family = neither")
end
for (k, v) in sort!(collect(tally); by = first)
    @printf("%6d  %s\n", v, k)
    for e in get(examples, k, String[])
        println("          ", e)
    end
end
