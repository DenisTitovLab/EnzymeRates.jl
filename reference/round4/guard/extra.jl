# ABOUTME: Extra statistics: the law's values at the random positive point (NaN, Inf, zero),
# ABOUTME: and the cost of each structural check against the support test, on a record sample.
include(joinpath(@__DIR__, "lib.jl"))
include(joinpath(@__DIR__, "approach_lib.jl"))
sl = deserialize(joinpath(R, "data", "slice_full.jls"))
spd = Dict(r.key => r.spec for r in deserialize(joinpath(R, "data", "specs_full.jls")))
cnt = Dict{String, Int}()
bump(k) = cnt[k] = get(cnt, k, 0) + 1
for s in sl
    haskey(s, :num) && s.num isa Dict && !haskey(s.num, :err) || continue
    n = s.num[:pos]
    bump("with numeric checks")
    isnan(n.vlaw) && bump("pos: law NaN")
    isinf(n.vlaw) && bump("pos: law Inf")
    n.vlaw == 0 && bump("pos: law exactly 0")
    (isnan(n.D) || n.D <= 0) && bump("pos: full-law D not > 0")
    (isnan(n.N) || isinf(n.N)) && bump("pos: full-law N not finite")
    for f in (:fwd, :rev)
        isnan(s.num[f].vlaw) && bump("$f: law NaN (D = 0 on the face)")
    end
end
println(sort!(collect(cnt)))

# costs of the structural checks and the support test (warm), on 2000 random records
from_spec(sp) = begin
    met(x) = S.Met(x[1], x[2])
    form(f) = S.RawForm(f[1], S.Met[met(x) for x in f[2]], collect(Symbol, f[3]), collect(Symbol, f[4]))
    steps = [S.RawStep(form(s[1]), form(s[2]), S.Met[met(x) for x in s[3]], S.Met[met(x) for x in s[4]], s[5]) for s in sp.steps]
    S.Mechanism(sp.substrates, sp.products, steps, sp.group; tags = sp.tags, multiplicity = sp.multiplicity,
                sites = [S.Site(s[1], s[2], s[3]) for s in sp.sites])
end
rng = Xoshiro(9)
keys_ = [s.key for s in sl if haskey(s, :law) && haskey(s.law, :reduced)]
pick = keys_[randperm(rng, length(keys_))[1:2000]]
ms = [from_spec(spd[k]) for k in pick]
t(f) = (f(ms[1]); v = Float64[]; for m in ms; push!(v, @elapsed f(m)); end; v)
using Statistics
q(v) = @sprintf("median %.2e s  p95 %.2e s  max %.2e s", median(v), quantile(v, 0.95), maximum(v))
println("n = ", length(ms), " (", count(S.is_allosteric, ms), " allosteric)")
println("bottomless (slice):             ", q(t(m -> S.bottomless(m))))
println("face_check, both faces:         ", q(t(m -> (face_check(m, m.products); face_check(m, m.substrates)))))
println("re_turnover:                    ", q(t(re_turnover)))
println("no turnover (cycles):           ", q(t(m -> S.cycles(m, S.orient(m)))))
println("admit:                          ", q(t(S.admit)))
println("solve:                          ", q(t(S.solve)))
println("solve + reduced law:            ", q(t(m -> S.rate_law(m, S.solve(m)))))
println("support test (law given):       ", q(begin
    laws = [(m, reduced_branches(m, S.solve(m))) for m in ms]
    f((m, rb)) = support_test(rb.N, rb.D, conc_indices(m, rb.vars)...)
    f(laws[1]); [@elapsed f(x) for x in laws] end))
println("full derive (pipeline):         ", q(t(m -> try S.derive(m) catch; end)[1:200]))
