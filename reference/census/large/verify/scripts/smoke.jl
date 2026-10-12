# ABOUTME: Smoke test of the verification oracle on a few census records: the law against the
# ABOUTME: pruned and vanishing networks at positive points, faces at two speed patterns, timings.
include(joinpath(@__DIR__, "vlib.jl"))

data = load_run(:ter_bi_rxn)
byid = Dict(r[:id] => r for recs in data.rounds for r in recs)
rng = Xoshiro(424242)

function positive_check(r; npts = 3)
    m = r[:smech]
    sol = S.solve(m)
    L = S.rate_law(m, sol)
    println("-- id $(r[:id]) labels $(r[:labels]) admit=$(r[:admit] === nothing ? "ok" : "REJECT")")
    println("   ", first(mech_text(m), 500))
    for ip in 1:npts
        vals = random_vals(m, rng)
        t0 = time()
        cols, x, resid, nfree = all_columns(m, vals, rng)
        tcol = time() - t0
        conc = random_conc(m, rng)
        vl = law_value(L, lawx(L, vals, conc))
        sp = random_speed(m, rng)
        t0 = time()
        vp, nc = oracle(m, cols, x, vals, conc; Λ = big(10.0)^30, speed = sp)
        tor = time() - t0
        line = @sprintf("   [pt %d] resid %.1e free %d | law %s | pruned oracle %s (rel %.1e, %d class, %.3f s, cols %.3f s)",
                        ip, Float64(resid), nfree, fmtb(vl), fmtb(vp),
                        Float64(abs(vp - vl) / abs(vl)), nc, tor, tcol)
        println(line)
        if S.is_allosteric(m)
            st = onlya_structure(m)
            sig = [abs(big(0.0)) + 2 + big(rand(rng)) for _ in m.edges]
            for (lab, off) in (("A", 3), ("B", 0), ("C", -3))
                phi = potentials(st, off, rng)
                vs = BigFloat[]
                for tt in (25, 35)
                    for k in eachindex(sig)
                        o = S.oriented(m, k, S.orient(m))
                        d = phi[st.comp[o.a]] - phi[st.comp[o.b]]
                        sig[k] = abs(d) / 2 + 2 + big(rand(rng))
                    end
                    v, ncl = oracle(m, cols, x, vals, conc; Λ = big(10.0)^12, speed = sp,
                                    mode = :vanish, phi = phi, t = big(tt), sigma = sig)
                    push!(vs, v)
                end
                @printf("      vanish %s (lost %d): t=25 %s t=35 %s rel to law %.2e\n", lab,
                        length(st.lost), fmtb(vs[1]), fmtb(vs[2]),
                        Float64(abs(vs[2] - vl) / abs(vl)))
            end
        end
    end
end

for id in ((0, 5), (1, 8), (1, 9), (1, 1157), (1, 1164))
    positive_check(byid[id])
end

# faces of the plain W_split record (2, 460)
r = byid[(2, 460)]
m = r[:smech]
L = S.rate_law(m)
println("-- faces of $(r[:id]): ", r[:faces])
vals = random_vals(m, rng)
cols, x, resid, _ = all_columns(m, vals, rng)
conc = random_conc(m, rng)
sp1, sp2 = random_speed(m, rng; dec = 3), random_speed(m, rng; dec = 3)
for (fname, zs) in vfaces(m)
    c = copy(conc)
    foreach(z -> c[z] = big(0.0), zs)
    vl = law_value(L, lawx(L, vals, c))
    o1, n1 = oracle(m, cols, x, vals, c; Λ = big(10.0)^14, speed = sp1)
    o2, n2 = oracle(m, cols, x, vals, c; Λ = big(10.0)^20, speed = sp1)
    o3, n3 = oracle(m, cols, x, vals, c; Λ = big(10.0)^20, speed = sp2)
    @printf("   %-6s law %s | oracle Λ1e14 %s, Λ1e20 %s, Λ1e20 pattern2 %s (classes %d %d %d)\n",
            fname, fmtb(vl), fmtb(o1), fmtb(o2), fmtb(o3), n1, n2, n3)
end
