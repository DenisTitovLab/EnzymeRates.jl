# ABOUTME: Splits the mutant results by toggle direction (RE to SS, the enumerator's flip move;
# ABOUTME: SS to RE; mixed): no-turnover laws, kcat classes, gauge failures, oracle disagreement.
using Serialization, Printf
const VD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify"
rows = deserialize(joinpath(VD, "data", "mutant_rows.jls"))
cnt(f, xs) = (d = Dict{Any, Int}(); foreach(x -> (k = f(x); d[k] = get(d, k, 0) + 1), xs); sort!(collect(d); by = string))
function oc_class(x)
    o1, o2 = x.oc
    any(o -> o.growth > 1e3, (o1, o2)) && return :grows_with_RE_speed
    any(o -> o.speed > 1e-6, (o1, o2)) && return :depends_on_RE_speed_pattern
    o1.law_vs_oracle > 1e-6 && return :wrong_at_positive_point
    o2.law_vs_oracle > 1e-6 && return :wrong_on_forward_face_only
    :matches
end
println("mutants by direction: ", cnt(x -> (x.allo, x.dir), rows))
L = filter(x -> x.law, rows)
for d in (:re_to_ss, :ss_to_re, :mixed)
    xs = filter(x -> x.dir === d, L)
    isempty(xs) && continue
    println("\n== ", d, ": derived ", length(xs), " (plain ", count(x -> !x.allo, xs), ")")
    println("   verdicts: ", cnt(x -> x.verdict, xs))
    println("   v identically 0 (Q = 0): ", count(x -> x.Qzero, xs))
    println("   kcat classes: ", cnt(x -> x.cls, xs))
    fr = filter(x -> !x.Qzero, xs)
    println("   forward-rate laws never saturating: ", count(x -> x.cls === :infinite, fr), " of ", length(fr))
    P = filter(x -> !x.allo && !x.Qzero, xs)
    println("   plain forward-rate laws vs oracle: ", cnt(oc_class, P))
end
