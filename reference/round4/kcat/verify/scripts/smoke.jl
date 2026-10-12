# ABOUTME: Smoke test of the verifier's helpers: loading, the mass-action oracle against the law
# ABOUTME: and at equilibrium, hand-built mechanism parsing, and the pipeline rows' fields.
include(joinpath(@__DIR__, "vlib.jl"))

t0 = time()
recs = [core_records("init_bibi"); core_records("expand_bibi")[1:50]; core_records("specs")]
@printf("loaded %d records in %.1f s\n", length(recs), time() - t0)
rng = Xoshiro(7)
worst, worsteq, n = 0.0, 0.0, 0
for r in recs
    m = r.m
    S.is_allosteric(m) && continue
    ok, why, sol, L = derive_law(m)
    ok || continue
    names, subs, prods, mods = roles(m)
    pd = draw_params(L.vars, rng)
    kv = edge_constants(m, sol, pd)
    c = [10.0^(-1 + 2 * rand(rng)) for _ in names]
    nl = NL(L, pd)
    vl = rate(nl, big.(c))
    vo, _ = oracle(m, kv, c)
    global worst = max(worst, Float64(abs(vl - vo) / abs(vo)))
    # equilibrium: scale the products so that ΠP/ΠS = Keq
    ce = big.(c)
    ratio = prod(ce[prods]) / prod(ce[subs])
    f = (pd[:Keq] / ratio)^(1 / length(prods))
    ce[prods] .*= f
    ve, _ = oracle(m, kv, ce)
    vref, _ = oracle(m, kv, (cc = big.(c); cc[prods] .= 0; cc))
    global worsteq = max(worsteq, Float64(abs(ve) / abs(vref)))
    global n += 1
end
@printf("plain laws %d: oracle vs law worst rel %.2e; |v| at equilibrium over |v| at P=0 worst %.2e\n",
        n, worst, worsteq)

# hand-built parsing
hm = mech([:A, :B], [:P, :Q], ["E + A <=> E(A)", "E(A) + B <--> E(A, B)",
                               "E(A, B) <--> E(P, Q)", "E(P, Q) <--> E(Q) + P",
                               "E(Q) <--> E + Q"])
println("hand-built ordered bi-bi:")
foreach(l -> println("   ", l), steps_text(hm))
ok, why, sol, L = derive_law(hm)
println("  derivable ", ok, " ", why, "; kcat class ", kcat_class(hm, sol)[1])

rows = deserialize(joinpath(IDATA, "pipe_terter11.jls"))
println("pipe row keys: ", keys(rows[1]))
println("pipe row kcat sample: ", [x.kcat for x in rows[1:5] if x.ok])
@printf("done in %.1f s, maxrss %.0f MB\n", time() - t0, Sys.maxrss() / 2^20)
