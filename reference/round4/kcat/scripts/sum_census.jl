# ABOUTME: Summarizes the (partial) init_terter census: kcat classes of the core and full law,
# ABOUTME: per-substrate limit signs of kcat-zero and no-unique-limit laws, and the gauge check.
using Serialization, Printf
const KDATA = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/data"
out = deserialize(joinpath(KDATA, "census_terter.jls"))
cnt(f, xs = out) = (d = Dict{Any, Int}(); foreach(x -> (y = f(x); d[y] = get(d, y, 0) + 1), xs); sort!(collect(d); by = last, rev = true))
println("parts ", sort!(unique(x.part for x in out)), ": ", length(out), " derivable laws")
for k in sort!(unique(x.part for x in out))
    xs = filter(x -> x.part == k, out)
    println("  part $k: ", length(xs), " laws; core: ", cnt(x -> x.cls, xs), "; full: ", cnt(x -> x.fcls, xs))
end
println("kcat class (core): ", cnt(x -> x.cls))
println("full-law class:    ", cnt(x -> x.fcls))
println("core => full:      ", cnt(x -> (x.cls, x.fcls)))
println("agreement with the pipeline's verdict (no unique limit or not): ",
        count(x -> (x.cls === :nounique) == occursin("depends on how", x.pipe), out), " of ", length(out))
for c in (:zero, :nounique)
    xs = filter(x -> x.cls === c, out)
    println(c, " (", length(xs), "): sign of deg_X(S*Q) - deg_X(D) per substrate (A, B, D), core law: ", cnt(x -> sign.(x.gap), xs))
    println(c, ": same, full law: ", cnt(x -> sign.(x.gapfull), xs))
    println(c, ": laws where some substrate alone drives v to 0 on the core law: ", count(x -> any(<(0), x.gap), xs),
            "; on the full law: ", count(x -> any(<(0), x.gapfull), xs))
end
println("gauge (3 parameter points x 50 vectors, 10 product-free): laws with a NaN ", count(x -> x.nnan > 0, out),
        ", an Inf ", count(x -> x.ninf > 0, out), ", max|v| not finite > 0 at some point ",
        count(x -> !(0 < x.maxabs < Inf), out), ", product-free max|v| not finite > 0 ",
        count(x -> !(0 < x.maxabs0 < Inf), out))
