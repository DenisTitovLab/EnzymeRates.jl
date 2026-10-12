# ABOUTME: Verifies classes D and Z on package laws: identity by name and Levenberg–Marquardt family
# ABOUTME: fits both ways at random points, for flagged pairs and unflagged candidates against parents.
include(joinpath(@__DIR__, "vlib.jl"))
const SAMPLES = deserialize(joinpath(VER, "samples.jls"))

rng = Xoshiro(515151)
tally = Dict{String, Dict{String, Int}}()
function bump!(lab, verdict)
    t = get!(tally, lab, Dict{String, Int}())
    t[verdict] = get(t, verdict, 0) + 1
end
verdict(ab, ba) = ab < 1e-8 && ba < 1e-8 ? "same" : ba < 1e-8 ? "first contains second" :
                  ab < 1e-8 ? "second contains first" : "different"

function pair(lab, name, d, mechs, x, y)
    rx, ry = d.recs[x], d.recs[y]
    a, b = plaw(mechs[x]), plaw(mechs[y])
    byname = same_by_name(a, b, rng)
    ab, ba = same_family(a, b, rng)        # ab: x's law inside y's family; ba: y's inside x's
    v = verdict(ab, ba)
    bump!(lab, v)
    @printf("  [%s] %s #%d (np %d rank %d) vs #%d (np %d rank %d): by-name %s | x-in-y %.1e y-in-x %.1e → %s\n",
            lab, name, x, rx[:np], rx[:rank_pkg], y, ry[:np], ry[:rank_pkg],
            isnan(byname) ? "names differ" : @sprintf("%.1e", byname), ab, ba, v)
    v
end

for name in ["uni_uni_allo_2reg", "uni_bi_rxn", "bi_bi_rxn", "bi_bi_pp_rxn",
             "bi_bi_rxn.r2_100000", "bi_bi_pp_rxn.r2_1200"]
    d, mechs = loadrun(name)
    s = SAMPLES[name]
    println("\n== $name")
    for (x, y) in s[:D_pairs]
        pair("D flagged (random)", name, d, mechs, x, y)
    end
    npd = s[:D_pairs_npdiff]
    for (x, y) in npd[sort!(randperm(rng, length(npd))[1:min(8, length(npd))])]
        pair("D flagged, package np/rank differ", name, d, mechs, x, y)
    end
    for (x, p) in s[:Z]
        pair("Z flagged (child, parent)", name, d, mechs, x, p)
    end
    for (x, p) in s[:Z_unflagged][1:min(end, 10)]
        pair("Z unflagged (child, parent)", name, d, mechs, x, p)
    end
    for x in s[:D_unflagged][1:min(end, 10)]
        p = d.recs[x][:parent]
        p == 0 && continue
        pair("D unflagged (candidate, parent)", name, d, mechs, x, p)
    end
    d = mechs = nothing
    GC.gc()
end
println("\nTally:")
for (k, t) in sort!(collect(tally))
    println("  $k: ", join(("$v $n" for (v, n) in sort!(collect(t))), ", "))
end
