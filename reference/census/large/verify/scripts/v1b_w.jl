# ABOUTME: Follow-ups on W: convergence of each vanishing pattern for a few island-rule cases, and
# ABOUTME: the full-face trap verdicts against many RE speed orderings and parameter points.
include(joinpath(@__DIR__, "vlib.jl"))

rel(a, b) = Float64(abs(a - b) / max(abs(a), abs(b), big(1e-60)))
rng = Xoshiro(20261011 + 11)

println("PART 1: each pattern at t = 30, 45, 60 (Λ = 1e12)")
for (name, ids) in ((:ter_bi_rxn, [(1, 1157), (1, 1275), (2, 577)]),
                    (:ter_ter_rxn, [(1, 2306), (1, 2390)]),
                    (:pyruvate_carboxylase_rxn, [(1, 1799), (1, 1057)]))
    data = load_run(name)
    byid = Dict(r[:id] => r for recs in data.rounds for r in recs)
    for id in ids
        r = byid[id]
        m = r[:smech]
        L = S.rate_law(m)
        st = onlya_structure(m)
        fwd = S.orient(m)
        vals = random_vals(m, rng)
        conc = random_conc(m, rng)
        sp = random_speed(m, rng)
        vl = law_value(L, lawx(L, vals, conc))
        cx = all_columns(m, vals, rng)
        println("-- $name $id law $(fmtb(vl)); ", first(r[:admit], 160))
        println("   ", first(r[:text], 600))
        for (lab, off, sigfix) in (("A", 3, nothing), ("B", 0, nothing), ("B, other speeds", 0, 1),
                                   ("C", -3, nothing))
            phi = potentials(st, off, rng)
            sig = map(eachindex(m.edges)) do k
                o = S.oriented(m, k, fwd)
                abs(phi[st.comp[o.a]] - phi[st.comp[o.b]]) / 2 + 2 +
                    (sigfix === nothing ? big(rand(rng)) : big(3 * rand(rng)))
            end
            vs = [first(oracle(m, cx[1], cx[2], vals, conc; Λ = big(10.0)^12, speed = sp,
                               mode = :vanish, phi = phi, t = big(tt), sigma = sig))
                  for tt in (30, 45, 60)]
            @printf("   %-16s v(t=30,45,60) = %s, %s, %s; |v60 - law|/|law| %.2e\n", lab,
                    fmtb(vs[1]), fmtb(vs[2]), fmtb(vs[3]), rel(vs[3], vl))
        end
    end
end

println("PART 2: every face of the face-flagged records, 20 speed orderings in 10^[-6,6] at 3 parameter points")
for name in (:ter_bi_rxn, :ter_ter_rxn, :pyruvate_carboxylase_rxn)
    data = load_run(name)
    fit = [r for recs in data.rounds for r in recs if r[:fitted]]
    for r in fit
        fs = [f for f in get(r, :faces, []) if any(w -> occursin(w, f[2]), ("split", "trap", "bottomless"))]
        isempty(fs) && continue
        m = r[:smech]
        L = S.rate_law(m)
        println("-- $name $(r[:id]): census ", fs)
        println("   ", first(r[:text], 900))
        for (fname, zs) in vfaces(m)
            gaps, spreads = Float64[], Float64[]
            for _ in 1:3
                vals = random_vals(m, rng)
                conc = random_conc(m, rng)
                foreach(z -> conc[z] = big(0.0), zs)
                cx = all_columns(m, vals, rng)
                vl = law_value(L, lawx(L, vals, conc))
                os = BigFloat[]
                for _ in 1:20
                    v, n = oracle(m, cx[1], cx[2], vals, conc; Λ = big(10.0)^30,
                                  speed = random_speed(m, rng; dec = 6))
                    push!(os, v)
                end
                push!(gaps, maximum(rel(o, vl) for o in os))
                push!(spreads, rel(minimum(os), maximum(os)))
            end
            @printf("   %-14s max |oracle - law|/|law| %s; oracle spread over orderings %s\n",
                    fname, join((@sprintf("%.1e", g) for g in gaps), " "),
                    join((@sprintf("%.1e", s) for s in spreads), " "))
        end
        flush(stdout)
    end
end
