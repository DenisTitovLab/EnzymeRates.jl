# ABOUTME: Summarizes zero_rows.jls: kcat classes per population, exact and numeric checks of the
# ABOUTME: kcat-zero laws, per-substrate gaps, traps, leaf dead ends, laws with no forward rate.
using Serialization, Printf, Statistics
const VD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify"
rows = deserialize(joinpath(VD, "data", "zero_rows.jls"))
cnt(f, xs) = (d = Dict{Any, Int}(); foreach(x -> (k = f(x); d[k] = get(d, k, 0) + 1), xs);
              sort!(collect(d); by = string))
ok = filter(x -> x.ok, rows)
println("== classes per population (slice kcat_candidates)")
for pop in ("core", "expand_terter", "terter_sample")
    xs = filter(x -> x.pop == pop, rows)
    oks = filter(x -> x.ok, xs)
    println(pop, ": records ", length(xs), ", derivable ", length(oks), ", classes ",
            cnt(x -> x.cls, oks), ", not derivable ", cnt(x -> first(x.why, 40), filter(x -> !x.ok, xs)))
end

println("\n== independent numeric check of the class (2 points x 2 directions, uniform scaling)")
fin = filter(x -> x.cls === :finite, ok)
relk(p, s) = abs(s[3] / p.kc - 1)
badfin = filter(x -> any(p -> relk(p, p.s1) > 1e-3 || relk(p, p.s2) > 1e-3, x.pts), fin)
println("finite: ", length(fin), " laws; v(1e12 r)/kcat off by > 1e-3 in ", length(badfin))
foreach(x -> println("   ", x.set, "/", x.index, " ", [(p.kc, p.s1[3], p.s2[3]) for p in x.pts]),
        badfin[1:min(5, end)])
zer = filter(x -> x.cls === :zero, ok)
dec(s) = s[3] < 1e-6 * maximum(s[1:2]) || s[3] == 0
println("zero: ", length(zer), " laws (plain ", count(x -> !x.allo, zer), "); numeric v(1e12 r) < 1e-6 max(v(1e4 r), v(1e8 r)) in all scans: ",
        count(x -> all(p -> dec(p.s1) && dec(p.s2), x.pts), zer))
println("   ratio v(1e12)/v(1e8) quantiles: ",
        quantile([p.s1[3] / p.s1[2] for x in zer for p in x.pts], [0, 0.5, 1]))
nou = filter(x -> x.cls === :nounique, ok)
println("nounique: ", length(nou), "; directions give different v(1e12) (rel > 1e-3) at some point: ",
        count(x -> any(p -> abs(p.s1[3] / p.s2[3] - 1) > 1e-3, x.pts), nou))
inf = filter(x -> x.cls === :infinite, ok)
println("infinite: ", length(inf), "; v(1e12)/v(1e8) > 10 in some scan: ",
        count(x -> any(p -> p.s1[3] / p.s1[2] > 10 || p.s2[3] / p.s2[2] > 10, x.pts), inf))
println("other: ", count(x -> x.cls === :other, ok))

println("\n== kcat-zero laws: exact positivity at zero products (and modifiers)")
println("all coefficients positive (every component): ", count(x -> x.zf.pos, zer), " of ", length(zer))
println("numerator keeps a term at P=0: ", count(x -> x.zf.num0, zer), "; denominator keeps a term: ",
        count(x -> x.zf.den0, zer))
println("=> v > 0 at every positive substrate point (exact): ",
        count(x -> x.zf.pos && x.zf.num0 && x.zf.den0, zer))
println("numeric: 20 random substrate points x 2 parameter points all v > 0: ",
        count(x -> all(p -> p.npos == 20, x.pts), zer))
println("uniform-scaling degree gap (full law) < 0: ", count(x -> x.zf.ugap !== nothing && x.zf.ugap < 0, zer),
        "; core law: ", count(x -> x.zc.ugap !== nothing && x.zc.ugap < 0, zer))
println("oracle vs law (plain, P=0, 2 points) worst rel: ",
        maximum(x -> x.orc, filter(x -> !x.allo, zer)))

println("\n== per-substrate gaps deg_X(num) - deg_X(den), X alone growing")
g(v) = v === nothing ? 99 : Int(v)
for (lab, f) in (("full law", x -> x.zf.gap), ("core law", x -> x.zc.gap))
    println(lab, ": any gap > 0 (unbounded growth): ", count(x -> any(v -> v !== nothing && v > 0, f(x)), zer),
            "; any gap undefined: ", count(x -> any(isnothing, f(x)), zer))
    for al in (false, true)
        xs = filter(x -> x.allo == al, zer)
        println("  ", al ? "allosteric" : "plain", " (", length(xs), "): culprit count ",
                cnt(x -> count(v -> v !== nothing && v < 0, f(x)), xs), "; gap values ",
                cnt(x -> sort(g.(f(x))), xs))
    end
end

println("\n== traps (plain kcat-zero laws; law x culprit substrate, one point)")
tr = [(x, t) for x in zer if !x.allo for t in x.traps]
println("traps: ", length(tr), " in ", length(unique(x.set * string(x.index) for (x, _) in tr)), " laws")
println("by (core culprit, full culprit): ", cnt(((x, t),) -> (t.core_culprit, t.full_culprit), tr))
println("share of enzyme held quantiles: ", quantile([t.share for (_, t) in tr], [0, 0.1, 0.5, 1]))
for cc in (true, false)
    ts = [(x, t) for (x, t) in tr if t.core_culprit == cc]
    println(" core culprit = ", cc, ": ", length(ts))
    println("   every trap form holds X (as reactant): ", count(((x, t),) -> t.all_hold_X_reactant, ts),
            "; holds X (any copy): ", count(((x, t),) -> t.all_hold_X, ts),
            "; holds a product or residual: ", count(((x, t),) -> t.all_hold_P, ts))
    println("   exit kinds: ", cnt(((x, t),) -> sort(unique(t.kinds)), ts))
    println("   has an exit that takes up a product: ", count(((x, t),) -> :blocked_takes_up_product in t.kinds, ts))
end
odd = [(x, t) for (x, t) in tr if :other in t.kinds || !t.all_hold_P || !t.all_hold_X || t.share < 0.5]
println("traps with an :other exit, a form without product/residual or X, or share < 0.5: ", length(odd))
for (x, t) in odd[1:min(4, end)]
    println("  ", x.set, "/", x.index, " X=", t.X, " share ", round(t.share; digits = 3), " forms ", t.forms)
    foreach(e -> println("      exit ", e), t.exits)
    foreach(l -> println("      ", l), x.lines)
end
for (x, t) in tr[1:min(2, end)]
    println("example ", x.set, "/", x.index, " X=", t.X, " share ", round(t.share; digits = 4),
            " forms ", t.forms)
    foreach(e -> println("      exit ", e), t.exits)
    foreach(l -> println("      ", l), x.lines)
end

println("\n== allosteric kcat-zero laws")
for x in filter(x -> x.allo, zer)
    println("  ", x.set, "/", x.index, " n=", x.n, " exact v>0: ", x.zf.pos && x.zf.num0 && x.zf.den0,
            " gaps full ", g.(x.zf.gap), " core ", g.(x.zc.gap))
end

println("\n== finite kcat but the full law (products and modifiers zero) falls to 0 under uniform scaling")
leaf = filter(x -> x.zf.ugap !== nothing && x.zf.ugap < 0, fin)
println("finite-kcat laws: ", length(fin), "; full-law uniform gap < 0: ", length(leaf),
        " (no modifiers: ", count(x -> x.nmods == 0, leaf), " of ", count(x -> x.nmods == 0, fin), ")")
println("  by set: ", cnt(x -> x.set, leaf))

println("\n== every derivable law: no forward rate at products and modifiers zero")
nof = filter(x -> !x.zf.num0 || !x.zf.den0, ok)
println("numerator or denominator identically 0 at P=0, M=0: ", length(nof))
foreach(x -> println("  ", x.set, "/", x.index, " num0 ", x.zf.num0, " den0 ", x.zf.den0,
                     " nmods ", x.nmods, " class ", x.cls), nof)
println("coefficients not all positive (any derivable law): ", count(x -> !x.zf.pos, ok))
