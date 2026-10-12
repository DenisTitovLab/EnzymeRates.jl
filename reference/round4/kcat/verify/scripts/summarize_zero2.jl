# ABOUTME: Breaks down the finite-kcat laws whose full law at saturating substrates (products and
# ABOUTME: modifiers zero) does not return kcat: modifiers, full-law zero, direction dependence, rest.
using Serialization, Printf, Statistics
const VD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify"
rows = deserialize(joinpath(VD, "data", "zero_rows.jls"))
cnt(f, xs) = (d = Dict{Any, Int}(); foreach(x -> (k = f(x); d[k] = get(d, k, 0) + 1), xs);
              sort!(collect(d); by = string))
fin = filter(x -> x.ok && x.cls === :finite, rows)
relk(p, s) = abs(s[3] / p.kc - 1)
bad(x) = any(p -> relk(p, p.s1) > 1e-3 || relk(p, p.s2) > 1e-3, x.pts)
dirdep(x) = any(p -> abs(p.s1[3] / p.s2[3] - 1) > 1e-3, x.pts)
function category(x)
    x.nmods > 0 && return :has_modifiers
    x.zf.ugap !== nothing && x.zf.ugap < 0 && return :full_law_to_zero
    dirdep(x) && return :full_law_direction_dependent
    any(p -> p.s1[3] > p.kc * (1 + 1e-3), x.pts) && return :full_above_kcat
    :full_below_kcat
end
b = filter(bad, fin)
println("finite-kcat laws ", length(fin), "; full law at 1e12 x random direction differs from kcat by > 1e-3: ",
        length(b))
println("categories: ", cnt(category, b))
println("by population x category: ", cnt(x -> (x.pop, category(x)), b))
for c in (:full_below_kcat, :full_above_kcat)
    xs = filter(x -> category(x) === c, b)
    println("\n", c, ": ", length(xs))
    for x in xs[1:min(4, end)]
        println("  ", x.set, "/", x.index, " kcat vs v(1e12 r): ",
                [(round(p.kc; sigdigits = 5), round(p.s1[3]; sigdigits = 5), round(p.s2[3]; sigdigits = 5)) for p in x.pts],
                " ugap full ", x.zf.ugap, " core ", x.zc.ugap)
    end
end
# finite laws without modifiers whose full law saturates in one way: v(1e12 r) vs kcat
clean = filter(x -> x.nmods == 0 && x.zf.ugap == 0 && !dirdep(x), fin)
println("\nfinite, no modifiers, full law uniform gap 0 and direction-independent at the 2 points: ",
        length(clean), "; of these v(1e12 r)/kcat off by > 1e-3: ", count(bad, clean))
r = [p.s1[3] / p.kc for x in clean if bad(x) for p in x.pts]
isempty(r) || println("   ratios v/kcat: ", quantile(r, [0, 0.5, 1]))
