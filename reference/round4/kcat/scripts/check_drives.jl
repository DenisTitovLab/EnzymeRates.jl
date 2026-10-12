# ABOUTME: Compares, for the kcat-zero ter-ter laws, the numeric growth verdicts with the exact
# ABOUTME: per-substrate degree gaps, and prints the scans where they disagree.
using Serialization, Printf
const KDATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/data"
rows = [x for l in ("terter_sample", "expand_terter") for x in deserialize(joinpath(KDATA, "rows_$l.jls"))
        if x.ok && x.cls === :zero && !x.allo]
println(length(rows), " plain kcat-zero ter-ter laws")
# per (law, point, substrate): exact gap on the full law vs the final slope of the scan
agree_slope, total = 0, 0
pat = Dict{Tuple{Int, Bool}, Int}()
for x in rows, p in x.pts, j in 1:x.nsubs
    g = x.zp.gap[j]
    s = p.z2.slopes[j]
    k = (g, s < -0.5)
    pat[k] = get(pat, k, 0) + 1
end
println("(exact gap on the full law, final slope log10 v(1e8)/v(1e7) < -0.5) => count over laws x points x substrates: ",
        sort!(collect(pat)))
# the diagonal: the exact total-degree verdict is kcat = 0 for every law here
println("final slope of the all-together scan: ", sort!(collect(Dict(round(p.z2.slopes[end]; digits = 1) => count(y -> round(y.z2.slopes[end]; digits = 1) == round(p.z2.slopes[end]; digits = 1), [q for x in rows for q in x.pts]) for x in rows for p in x.pts))))
bad = [(x, p) for x in rows for p in x.pts if any(j -> (x.zp.gap[j] < 0) != (p.z2.slopes[j] < -0.5), 1:x.nsubs)]
println(length(bad), " (law, point) pairs where gap < 0 and slope < -0.5 disagree for some substrate; first 6:")
for (x, p) in bad[1:min(end, 6)]
    println("  ", x.id, " gaps full ", x.zp.gap, " core ", x.zpc.gap, " slopes ", round.(p.z2.slopes; digits = 2),
            " drives ", p.z2.drives)
end
println("laws where the full-law and core-law gaps differ: ", count(x -> x.zp.gap != x.zpc.gap, rows))
println("laws where some substrate alone has gap < 0 (full law): ", count(x -> any(<(0), x.zp.gap), rows),
        "; core law: ", count(x -> any(<(0), x.zpc.gap), rows))
println("gap patterns (full law): ", sort!(collect(Dict(g => count(y -> y.zp.gap == g, rows) for g in unique(x.zp.gap for x in rows))); by = last, rev = true))
println("gap patterns (core law): ", sort!(collect(Dict(g => count(y -> y.zpc.gap == g, rows) for g in unique(x.zpc.gap for x in rows))); by = last, rev = true))
