# ABOUTME: Prints step lists for the report's examples: no-turnover mutants with a steady-state
# ABOUTME: step, RE-to-SS mutants whose law misses the oracle, and two kcat-zero laws.
using Serialization, Printf
const VD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify"
mr = deserialize(joinpath(VD, "data", "mutant_rows.jls"))
show_lines(x) = foreach(l -> println("      ", l), x.lines)
z = filter(x -> x.law && x.Qzero && any(l -> occursin("<-->", l), x.lines), mr)
println("no-turnover mutants with a steady-state step: ", length(z), " (", count(x -> !x.allo, z), " plain); directions ",
        [x.dir for x in z] |> v -> (d = Dict{Symbol,Int}(); foreach(k -> d[k] = get(d, k, 0) + 1, v); d))
for x in filter(x -> !x.allo, z)[1:min(2, end)]
    o1, o2 = x.oc
    @printf("  %s %s (%s) verdict %s kcat %s; law v = 0; oracle at a positive point %.3g at RE speed 1e30, x%.3g from 1e20\n",
            x.parent, x.how, x.dir, x.verdict, x.cls, o1.o30, o1.growth)
    show_lines(x)
end
for x in filter(x -> x.allo, z)[1:min(1, end)]
    println("  ", x.parent, " ", x.how, " (", x.dir, ") verdict ", x.verdict, " kcat ", x.cls)
    show_lines(x)
end
bad = filter(x -> x.law && !x.allo && !x.Qzero && x.dir === :re_to_ss &&
              any(o -> o.speed > 1e-6 || o.law_vs_oracle > 1e-6, x.oc), mr)
println("\nRE-to-SS plain mutants whose law misses the oracle: ", length(bad))
for x in bad
    o1, o2 = x.oc
    @printf("  %s %s verdict %s: face law %.5g, oracle %.5g (speed pattern changes it by %.2g); positive point law %.5g oracle %.5g\n",
            x.parent, x.how, x.verdict, o2.law, o2.o30, o2.speed, o1.law, o1.o30)
    show_lines(x)
end
zr = deserialize(joinpath(VD, "data", "zero_rows.jls"))
for (s, i) in (("init_terter.part7", 4757), ("allo_bibi", 1023))
    x = zr[findfirst(x -> x.set == s && x.index == i, zr)]
    println("\n", s, "/", i, " class ", x.cls)
    show_lines(x)
end
