# ABOUTME: D and Z rechecked with L held fixed under the canonical renaming (every pair by the law,
# ABOUTME: a sample by the oracle), and fingerprint-only matches tested for a shared law by fitting.
include(joinpath(@__DIR__, "vlib.jl"))
include(joinpath(@__DIR__, "vfit.jl"))
setprecision(BigFloat, 256)

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
rel(a, b) = Float64(abs(a - b) / max(abs(a), abs(b), big(1e-60)))
pvals(m, fv, keq) = (; (n => fv[n] for n in S.fitted_params(m))..., Keq = keq,
                     E_total = big(1.0))

"""Law values (and with `orc` pruned-oracle values) of `mx` at random fitted values and of
`my` at the same values moved through the canonical orders, `L` and `Keq` shared, at three
concentration points: the largest relative differences."""
function pair_check(mx, my, rng; orc = false)
    Lx, Ly = S.rate_law(mx), S.rate_law(my)
    cx, cy = S.canonical(Lx), S.canonical(Ly)
    cx.key == cy.key || return (law = NaN, orc = NaN)
    fx = Dict{Symbol, BigFloat}(n => lu10(rng, -1, 1) for n in S.fitted_params(mx))
    fy = Dict{Symbol, BigFloat}(b => fx[a] for (a, b) in zip(cx.order, cy.order))
    haskey(fx, :L) && (fy[:L] = fx[:L])
    for n in S.fitted_params(my)
        haskey(fy, n) || (fy[n] = lu10(rng, -1, 1))
    end
    keq = lu10(rng, -1, 1)
    vx, vy = vals_from(mx, pvals(mx, fx, keq)), vals_from(my, pvals(my, fy, keq))
    wl, wo = 0.0, 0.0
    for _ in 1:3
        conc = random_conc(mx, rng)
        lx, ly = law_value(Lx, lawx(Lx, vx, conc)), law_value(Ly, lawx(Ly, vy, conc))
        wl = max(wl, rel(lx, ly))
        orc || continue
        cxx, cyy = all_columns(mx, vx, rng), all_columns(my, vy, rng)
        ox = first(oracle(mx, cxx[1], cxx[2], vx, conc; Λ = big(10.0)^30,
                          speed = random_speed(mx, rng)))
        oy = first(oracle(my, cyy[1], cyy[2], vy, conc; Λ = big(10.0)^30,
                          speed = random_speed(my, rng)))
        wo = max(wo, rel(ox, oy), rel(ox, lx))
    end
    (law = wl, orc = orc ? wo : NaN)
end

rng = Xoshiro(20261011 + 6)
for name in RX
    data = load_run(name)
    byid = Dict(r[:id] => r for recs in data.rounds for r in recs)
    first_of = Dict{String, Any}()
    for recs in data.rounds, r in recs
        r[:fitted] && get!(first_of, r[:key], r)
    end
    fit = [r for recs in data.rounds for r in recs if r[:fitted]]
    println("="^100, "\n$name")
    D = [r for r in fit if :D in r[:labels]]
    pick = Set(D[sort!(randperm(rng, length(D))[1:min(12, length(D))])])
    worst, worst_o = 0.0, 0.0
    for r in D
        res = pair_check(first_of[r[:key]][:smech], r[:smech], rng; orc = r in pick)
        worst = max(worst, res.law)
        r in pick && (worst_o = max(worst_o, res.orc);
                      @printf("   D %s of %s: law %.1e, oracle %.1e\n", r[:id],
                              first_of[r[:key]][:id], res.law, res.orc))
    end
    @printf("(a) D %d pairs: worst law difference under the canonical renaming %.1e; oracle sample %d: worst %.1e\n",
            length(D), worst, length(pick), worst_o)
    Z = [r for r in fit if :Z in r[:labels]]
    for r in Z
        res = pair_check(byid[r[:parent]][:smech], r[:smech], rng; orc = true)
        @printf("   Z %s vs parent %s: law %.1e, oracle %.1e\n", r[:id], r[:parent], res.law, res.orc)
    end
    # (b) fingerprint-only matches, tested by fitting both ways
    sp = deserialize(joinpath(VER, "data", "v2_splits_$(name).jls"))
    sample = sp[sort!(randperm(rng, length(sp))[1:min(25, length(sp))])]
    tally = Dict{String, Int}()
    for (id, others) in sample
        r, o = byid[id], byid[others[1]]
        a, b = both_ways(r[:smech], o[:smech], rng)
        v = a < 1e-8 && b < 1e-8 ? "same law" : "different"
        tally[v] = get(tally, v, 0) + 1
        @printf("   fp-only %s vs %s (slice %d/%d vs %d/%d): fits %.1e %.1e -> %s\n", id,
                o[:id], r[:rank], r[:n_slice], o[:rank], o[:n_slice], a, b, v)
    end
    println("(b) fingerprint-only matches: $(length(sp)); sampled $(length(sample)): ", tally)
    flush(stdout)
end
