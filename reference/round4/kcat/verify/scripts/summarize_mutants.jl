# ABOUTME: Summarizes mutant_rows.jls: current verdicts, derivability without checks, laws with
# ABOUTME: no turnover or no forward rate, gauge failures, and the law against the oracle.
using Serialization, Printf, Statistics
const VD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify"
rows = deserialize(joinpath(VD, "data", "mutant_rows.jls"))
cnt(f, xs) = (d = Dict{Any, Int}(); foreach(x -> (k = f(x); d[k] = get(d, k, 0) + 1), xs);
              sort!(collect(d); by = string))
println("mutants: ", length(rows), " (allosteric ", count(x -> x.allo, rows), ")")
println("verdict x law derivable without checks: ", cnt(x -> (x.allo, x.verdict, x.law), rows))
println("law errors: ", cnt(x -> first(x.why, 60), filter(x -> !x.law, rows)))
L = filter(x -> x.law, rows)
gfail(x) = any(p -> p.nan > 0 || p.inf > 0 || !(0 < p.maxabs < Inf), x.g)
gfail0(x) = any(p -> !(0 < p.maxabs0 < Inf), x.g)
kind(x) = x.Qzero ? :no_turnover_Q0 : !x.zf.den0 ? :zero_over_zero_at_P0 :
          !x.zf.num0 ? :no_forward_rate_at_P0 : :forward_rate
println("\nderived laws: ", length(L))
println("(allo, verdict, kind): ", cnt(x -> (x.allo, x.verdict, kind(x)), L))
println("gauge fails (NaN/Inf or max|v| not finite > 0 at some parameter point): ",
        cnt(x -> (x.verdict, kind(x)), filter(gfail, L)))
println("product-free gauge vectors give no finite positive |v| at some point: ",
        cnt(x -> (x.verdict, kind(x)), filter(gfail0, L)))
println("kcat class by verdict: ", cnt(x -> (x.verdict, x.cls), L))

P = filter(x -> !x.allo, L)
println("\nplain derived laws: ", length(P))
wrong(o) = o.law_vs_oracle > 1e-6
diverge(o) = o.growth > 1e3
speed(o) = o.speed > 1e-6
function oc_class(x)
    o1, o2 = x.oc
    any(o -> diverge(o), (o1, o2)) && return :oracle_grows_with_RE_speed
    any(o -> speed(o), (o1, o2)) && return :oracle_depends_on_RE_speed_pattern
    wrong(o1) && return :law_wrong_at_positive_point
    wrong(o2) && return :law_wrong_on_forward_face_only
    :law_matches_oracle
end
println("(verdict, kind, oracle class): ")
foreach(p -> println("   ", p), cnt(x -> (x.verdict, kind(x), oc_class(x)), P))
for c in (:oracle_grows_with_RE_speed, :oracle_depends_on_RE_speed_pattern,
          :law_wrong_at_positive_point, :law_wrong_on_forward_face_only)
    xs = filter(x -> oc_class(x) === c && kind(x) === :forward_rate, P)
    isempty(xs) && continue
    println("\n", c, " with a forward rate (", length(xs), "; verdicts ", cnt(x -> x.verdict, xs), "):")
    for x in xs[1:min(3, end)]
        o1, o2 = x.oc
        @printf("  %s %s verdict %s: positive point law %.5g oracle %.5g (growth x%.3g, speed rel %.2g); face law %.5g oracle %.5g (growth x%.3g, speed rel %.2g)\n",
                x.parent, x.how, x.verdict, o1.law, o1.o30, o1.growth, o1.speed, o2.law, o2.o30,
                o2.growth, o2.speed)
        foreach(l -> println("      ", l), x.lines)
    end
end
for k in (:no_turnover_Q0, :no_forward_rate_at_P0, :zero_over_zero_at_P0)
    xs = filter(x -> kind(x) === k, L)
    isempty(xs) && continue
    x = xs[1]
    println("\nexample ", k, " (", length(xs), "; verdicts ", cnt(x -> x.verdict, xs), "): ", x.parent, " ", x.how)
    foreach(l -> println("      ", l), x.lines)
end
