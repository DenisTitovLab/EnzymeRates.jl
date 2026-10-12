# ABOUTME: Every fitted census candidate: the slice law against the mass-action oracle (pruned
# ABOUTME: n=1 network) at one random positive point and two RE speeds; RE turnover diverges.
include(joinpath(@__DIR__, "vlib.jl"))
setprecision(BigFloat, 256)

const RX = (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
results = Dict{Tuple{Symbol, Tuple{Int, Int}}, Any}()
for name in RX
    data = load_run(name)
    rng = Xoshiro(hash((name, :v0)))
    t0 = time()
    tally = Dict{Tuple{Int, Symbol}, Int}()
    worst = Dict{Int, Float64}()
    examples = String[]
    for (ri, recs) in enumerate(data.rounds), r in recs
        r[:fitted] || continue
        m = r[:smech]
        verdict, rel, resid = try
            L = S.rate_law(m)
            vals = random_vals(m, rng)
            conc = random_conc(m, rng)
            cols, x, resid, _ = all_columns(m, vals, rng)
            sp = random_speed(m, rng)
            vl = law_value(L, lawx(L, vals, conc))
            o1, n1 = oracle(m, cols, x, vals, conc; Λ = big(10.0)^14, speed = sp)
            o2, n2 = oracle(m, cols, x, vals, conc; Λ = big(10.0)^22, speed = sp)
            rel = Float64(abs(vl - o2) / max(abs(vl), abs(o2), big(1e-60)))
            v = (n1 != 1 || n2 != 1) ? :nonunique :
                abs(o2) > 1e4 * max(abs(o1), big(1e-60)) ? :diverges :
                rel > 1e-10 ? :mismatch : :agree
            v, rel, Float64(resid)
        catch e
            :error, NaN, NaN
        end
        r[:v0_verdict] = verdict
        k = (ri - 1, verdict)
        tally[k] = get(tally, k, 0) + 1
        verdict === :agree && (worst[ri - 1] = max(get(worst, ri - 1, 0.0), rel))
        if verdict !== :agree && length(examples) < 8
            push!(examples, "  $(r[:id]) $verdict rel $(rel) resid $(resid): " * first(r[:text], 300))
        end
        results[(name, r[:id])] = (verdict = verdict, rel = rel, resid = resid)
    end
    @printf("%s (%.0f s): verdicts %s\n", name, time() - t0, sort!(collect(tally)))
    println("  largest |law - oracle|/|law| among agreeing, per round: ", worst)
    foreach(println, examples)
    flush(stdout)
end
serialize(joinpath(VER, "data", "v0_positive.jls"), results)
