# ABOUTME: Refits the three candidates whose package law did not fit the deletion-model steady state in
# ABOUTME: v10, with more starts and three fresh constant draws, to tell optimizer misses from mismatches.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "moracle.jl"))
rng = Xoshiro(161616)
for (name, g) in (("bi_bi_rxn.r2_100000", 16287), ("bi_bi_rxn.r2_100000", 18971), ("bi_bi_pp_rxn.r2_1200", 13220))
    d, mechs = loadrun(name)
    ms = mstruct(mechs[g])
    pl = plaw(mechs[g])
    println("\n$name #$g np $(length(pl.names)) ", d.recs[g][:text])
    for k in 1:3
        A = aconsts(ms, rng)
        Ip = iconsts(ms, A, rng)
        lnkp = Ip.lnk(Ip.draw(rng))
        L = big(lu10(rng, -1, 1))
        σ = [big(lu10(rng, -1, 1)) for _ in 1:(length(ms.steps) + 1)]
        keq = Float64(exp(BigFloat(sum(A.mu[x] for x in ms.subs) - sum(A.mu[x] for x in ms.prods))))
        npts = 2 * length(pl.names) + 14
        pts = [Dict{Symbol, Float64}(x => lu10(rng, -1, 1) for x in ms.mets) for _ in 1:npts]
        v = [Float64(oracle_rate(ms, A, Ip, lnkp, L, Dict{Symbol, BigFloat}(x => big(y) for (x, y) in c),
                                 big(10)^30, σ)) for c in pts]
        r = fit_values(pl, keq, pts, v, rng; starts = 80, iters = 1500)
        @printf("  draw %d: best residual over 80 starts %.1e\n", k, r)
        flush(stdout)
    end
end
