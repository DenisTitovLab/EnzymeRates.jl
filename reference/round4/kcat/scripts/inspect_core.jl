# ABOUTME: Prints the edge cases of the core census: laws whose rate vanishes at zero products,
# ABOUTME: full-law 0/0, kcat infinite, full-law kcat zero, and sup far below kcat.
include(joinpath(@__DIR__, "klib.jl"))

rows = filter(x -> x.ok, deserialize(joinpath(KDATA, "rows_core.jls")))
recs = Dict((r.set, r.c.index) => r for set in CORE_SETS for r in records_of(set))

function show_law(x; title = "")
    r = recs[(x.set, x.index)]
    m = r.c.mech
    println("  ", x.id, " (", x.set, "/", x.index, ") ", title)
    foreach(l -> println("      ", l), step_lines(m))
    println("      kcat class ", x.cls, " ", x.msg, " | full law ", x.fcls, " ", x.fmsg)
    if !x.allo
        sol = S.solve(m)
        println("      zero products, full law: ", zero_product_limits(m, S.law(m, sol)),
                "; core law: ", zero_product_limits(m, S.law(m, sol, S.catalytic_core(m))))
    end
    for p in x.pts
        @printf("      kcat %.4g full-law limit %.4g sup %.4g (subs only %.4g) max|v| %.3g (product-free %.3g)\n",
                p.kc, p.kf, p.sup, p.sups, p.maxabs, p.maxabs0)
    end
end

println("== laws with v = 0 at every product-free gauge vector")
foreach(x -> show_law(x), filter(x -> any(p -> !(0 < p.maxabs0 < Inf), x.pts), rows))
println("== full-law class other")
foreach(x -> show_law(x), filter(x -> x.fcls === :other, rows))
println("== kcat infinite (never saturates)")
foreach(x -> show_law(x), filter(x -> x.cls === :infinite, rows))
println("== no unique limit")
foreach(x -> show_law(x), filter(x -> x.cls === :nounique || x.fcls === :nounique, rows))
println("== core finite, full law zero: count per set and three examples")
fz = filter(x -> x.cls === :finite && x.fcls === :zero, rows)
println("  ", sort!(collect(Dict(s => count(x -> x.set == s, fz) for s in unique(x.set for x in fz)))))
foreach(x -> show_law(x), fz[[1, 2, end]])
println("== sup < 0.5 kcat with the full-law limit >= 0.5 kcat")
foreach(x -> show_law(x), filter(x -> x.cls === :finite && x.ratio < 0.5 &&
                                      all(p -> p.kf >= 0.5 * p.kc, x.pts), rows))
println("== three allosteric laws with sup > 1.01 kcat (n = 2), smallest ratios above 1.01")
al = sort(filter(x -> x.allo && x.cls === :finite && x.ratio > 1.01, rows); by = x -> x.ratio)
foreach(x -> show_law(x; title = @sprintf("ratio %.3f", x.ratio)), al[1:3])
println("== plain laws with sup > 1.01 kcat")
pl = sort(filter(x -> !x.allo && x.cls === :finite && x.ratio > 1.01, rows); by = x -> x.ratio)
println("  ", length(pl), " laws; ratios ", round.([x.ratio for x in pl]; digits = 3))
foreach(x -> show_law(x; title = @sprintf("ratio %.3f", x.ratio)), pl[1:3])
