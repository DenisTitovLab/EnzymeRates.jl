# ABOUTME: The README usage example, verbatim except the path, run as a check.
# ABOUTME: Stdlib only.
using Random
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto/common/Oracle.jl")          # defines Main.Pop and Main.Oracle
m = only(filter(m -> m.source == "spec:Segel Ordered Bi Bi", Pop.load_population("specs")))
s = m.samples[1]; vals = Pop.all_values(s)    # one package evaluation: fitted, Keq, E_total, dependents
k = Oracle.constants_from_names(m, vals; rename = Dict(m.pass2_rename))
v = Oracle.steady_state_v(m, k, Pop.concs(s)) # BigFloat
@assert abs(v - s.v) <= 1e-12 * abs(s.v)      # package Reduced rate (E_total = 1)
kc, _, _ = Oracle.consistent_constants(m, Random.Xoshiro(1); Keq = 3.0)
c = Oracle.concs_at_ratio(m, Random.Xoshiro(2), 3.0)
@assert abs(Oracle.steady_state_v(m, kc, c)) < 1e-30   # detailed balance at Γ = Keq
println("README example OK: v = ", Float64(v), " vs package ", s.v)
