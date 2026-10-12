# ABOUTME: Focused queries over a run's joined verdicts (bottomless oracle values, package
# ABOUTME: laws of N≡0 mutants, admission cross-tab), printing examples as step lists.
include(joinpath(@__DIR__, "report.jl"))   # main() runs; reuse its helpers afterwards
println("\n\n======== QUERIES ========")
P(s) = pkg_verdict(get(pk, s.key, nothing))
println("\n## bottomless (face check) records: oracle values")
for s in sl
    haskey(s, :struct_) || continue
    f = something(s.struct_.face_fwd, s.struct_.face_rev, "")
    occursin("bottomless", f) || continue
    example(s; extra = "\n      pkg: $(P(s)); slice: $(slice_verdict(s)); law: $(law_red(s)); face: $f")
    for f in (:pos, :fwd, :rev)
        n = s.num[f]
        @printf("      oracle %s: %s law=%.6g o=%s\n", f, n.verdict, n.vlaw, string(round.(n.o; sigdigits = 7)))
    end
end
println("\n## package law text of N≡0 records (first 6 of each package verdict)")
for v in ("accept", "accept (pkg N≡0)")
    ks = [s for s in sl if law_red(s) == "(b) N≡0" && P(s) == v]
    println("  $v: $(length(ks)); allosteric: $(count(s -> s.allo, ks))")
    for s in ks[1:min(6, end)]
        println("    ", s.key, " allo=", s.allo, " vline: ", first(pk[s.key].vline, 160))
    end
end
table("slice boundary × package verdict", "slice", "package", slice_verdict, P)
table("[allosteric] slice admit × law", "slice", "law", slice_verdict, s -> lawname(law_red(s));
      recs = [s for s in sl if get(s, :allo, false)])
println("\n## slice admit rejections the package accepts (first 6)")
ks = [s for s in sl if slice_verdict(s) == "admit" && startswith(P(s), "accept")]
println("  n = $(length(ks)); by set: ", countmap_([s.set for s in ks]))
for s in ks[1:min(4, end)]
    example(s; extra = "\n      admit: " * first(s.struct_.admit, 200) * "\n      law: $(law_red(s)) oracle: $(oracle_verdict(s))")
end
println("\n## package OnlyA guard rejections the slice admits")
ks = [s for s in sl if P(s) == "ctor: OnlyA guard" && slice_verdict(s) == "pass"]
println("  n = $(length(ks))")
foreach(s -> example(s; extra = "\n      pkg: " * first(pk[s.key].ctor, 200)), ks[1:min(3, end)])
