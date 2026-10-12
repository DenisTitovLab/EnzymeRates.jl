# ABOUTME: Checks whether the never-saturating RE/SS mutants' laws agree with the mass-action
# ABOUTME: oracle (genuine never-saturating mechanisms) or come from ill-posed RE limits.
using Serialization, Printf
const VD = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify"
rows = deserialize(joinpath(VD, "data", "mutant_rows.jls"))
cnt(f, xs) = (d = Dict{Any, Int}(); foreach(x -> (k = f(x); d[k] = get(d, k, 0) + 1), xs); sort!(collect(d); by = string))
function oc_class(x)
    o1, o2 = x.oc
    any(o -> o.growth > 1e3, (o1, o2)) && return :oracle_grows_with_RE_speed
    any(o -> o.speed > 1e-6, (o1, o2)) && return :oracle_depends_on_RE_speed_pattern
    o1.law_vs_oracle > 1e-6 && return :law_wrong_at_positive_point
    o2.law_vs_oracle > 1e-6 && return :law_wrong_on_forward_face_only
    :law_matches_oracle
end
P = filter(x -> x.law && !x.allo && !x.Qzero, rows)
println("plain forward-rate mutants by (kcat class, oracle class): ")
foreach(p -> println("   ", p), cnt(x -> (x.cls, oc_class(x)), P))
