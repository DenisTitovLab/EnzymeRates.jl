# ABOUTME: Prints flagged enumerator children (parent, flipped steps, child) with the oracle's
# ABOUTME: face values under two RE speed patterns against the law's value.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "approach_lib.jl"))
fl = deserialize(joinpath(R, "data", "enum_flagged.jls"))
from_spec(sp) = begin
    met(x) = S.Met(x[1], x[2])
    form(f) = S.RawForm(f[1], S.Met[met(x) for x in f[2]], collect(Symbol, f[3]), collect(Symbol, f[4]))
    steps = [S.RawStep(form(s[1]), form(s[2]), S.Met[met(x) for x in s[3]], S.Met[met(x) for x in s[4]], s[5]) for s in sp.steps]
    S.Mechanism(sp.substrates, sp.products, steps, sp.group; tags = sp.tags, multiplicity = sp.multiplicity,
                sites = [S.Site(s[1], s[2], s[3]) for s in sp.sites])
end
txt(m) = [S.step_text(m, e) * (S.is_allosteric(m) ? "  ::$(m.tags[m.group[k]])" : "") for (k, e) in enumerate(m.edges)]
rng = Xoshiro(21)
spread = Float64[]
for (i, r) in enumerate(fl)
    p, c = from_spec(r.parent), from_spec(r.child)
    lt = law_test(c)
    num = numeric_checks(c, lt[:law], rng)
    for f in (:fwd, :rev)
        n = num[f]
        n.verdict === :pattern && push!(spread, abs(n.o[2] - n.o[3]) / max(abs(n.o[2]), abs(n.o[3])))
    end
    i in (1, 9, 25) || continue
    println("\n=== $(r.level) $(r.move) child $i: $(r.face)")
    pt, ct = txt(p), txt(c)
    println("  parent:"); foreach(l -> println("    ", l), pt)
    println("  child (lines that differ from the parent marked *):")
    foreach(l -> println("   ", l in pt ? " " : "*", l), ct)
    for f in (:pos, :fwd, :rev)
        n = num[f]
        @printf("  %s: %s law=%.6g oracle Λ=1e10,1e20 (pattern 1), 1e20 (pattern 2) = %s\n", f, n.verdict, n.vlaw, string(round.(n.o; sigdigits = 6)))
    end
end
using Statistics
@printf("\nenumerator children: spread between two speed patterns at the flagged face: n=%d median %.2e max %.2e\n",
        length(spread), median(spread), maximum(spread))
