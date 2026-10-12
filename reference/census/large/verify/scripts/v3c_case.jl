# ABOUTME: A closer look at one L-phantom child whose parent did not reproduce it: many starts,
# ABOUTME: more points, a start from the child's values mapped by name with the state tag removed.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "vfit.jl"))
data = load_run(:ter_bi_rxn)
byid = Dict(r[:id] => r for recs in data.rounds for r in recs)
rng = Xoshiro(77)
untag(n) = Symbol(replace(string(n), r"^([Kk])_A_" => s"\1_"))
c, p = byid[(1, 1255)], byid[(0, 1996)]
nc, vc = law_fn(c[:smech]); np_, vp = law_fn(p[:smech])
for trial in 1:12
    θc = Dict(n => 10.0^(2rand(rng) - 1) for n in nc)
    keq = 10.0^(2rand(rng) - 1)
    pts = [Dict(x => 10.0^(2rand(rng) - 1) for x in S.metabolites(c[:smech])) for _ in 1:60]
    ys = [vc(θc, keq, q) for q in pts]
    mapped = Dict(untag(n) => v for (n, v) in θc)
    starts = [Dict(n => get(mapped, n, 1.0) for n in np_);
              [Dict(n => get(mapped, n, 1.0) * 10.0^(rand(rng) - 0.5) for n in np_) for _ in 1:10];
              [Dict(n => 10.0^(4rand(rng) - 2) for n in np_) for _ in 1:40]]
    r1 = lm_fit(np_, vp, keq, pts, ys, starts; iters = 600)
    @printf("  trial %2d: parent fits child: rel RMS %.2e (L = %.3g)\n", trial, r1, θc[:L])
end
