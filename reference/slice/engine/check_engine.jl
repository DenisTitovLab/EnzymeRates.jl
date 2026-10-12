# ABOUTME: Engine-stage checks (a), (b), (d), (e) on every plain record: the law against the
# ABOUTME: BigFloat oracle and the package, rational exponents, units and printed-law values.
include(joinpath(@__DIR__, "lib.jl"))

records = [r for r in load_records() if !S.is_allosteric(r.c.mech)]
rng = Xoshiro(20261010)
rows = NamedTuple[]
fails = String[]
note(s) = (push!(fails, s); length(fails) <= 40 && println(s))
t0 = time()
for r in records
    m, p = r.c.mech, r.p
    sol, L = try
        sol = S.solve(m)
        sol, S.law(m, sol)
    catch err
        note("$(r.set) $(p.id): no law: $(sprint(showerror, err))")
        continue
    end
    # (a) the oracle at 3 random consistent points through all_parameters
    worst_a = big(0)
    for _ in 1:3
        pars = random_params(m, rng)
        k = S.pop_constants(r.c, edge_values(m, S.all_parameters(m, pars)))
        u = (0.3 + 1.7 * rand(rng)) * (rand(rng, Bool) ? 1 : -1)
        conc = Oracle.concs_at_ratio(p, rng, pars.Keq * big(10)^u)
        vo = Oracle.steady_state_v(p, k, conc)
        worst_a = max(worst_a, abs(lawval(L, conc, pars) - vo) / abs(vo))
    end
    worst_a > 1e-30 && note("(a) $(r.set) $(p.id): oracle rel $(Float64(worst_a))")
    # (b) the package's law: its samples, and random points on its exact polynomials
    worst_s, worst_r, nsamp = NaN, NaN, 0
    if !isempty(p.pkg_den)
        worst_s, worst_r = 0.0, big(0)
        for s in p.samples
            pars = my_params(r, sol, pkg_values(p, s))
            v = lawval(L, Dict(s.concs), pars)
            worst_s = max(worst_s, Float64(abs(v - s.v) / abs(v)))
            nsamp += 1
        end
        for _ in 1:2
            pars = random_params(m, rng)
            k = S.pop_constants(r.c, edge_values(m, S.all_parameters(m, pars)))
            vals = Dict{Symbol, BigFloat}(:Keq => pars.Keq)
            for (j, ns) in enumerate(p.step_names), (q, nm) in enumerate(ns)
                vals[nm] = k[j][q]
            end
            conc = Oracle.concs_at_ratio(p, rng, pars.Keq * 10.0^(4 * rand(rng) - 2))
            merge!(vals, conc)
            vp = Oracle.eval_poly(p.pkg_num, vals) / Oracle.eval_poly(p.pkg_den, vals)
            worst_r = max(worst_r, abs(lawval(L, conc, pars) - vp) / abs(vp))
        end
        worst_r > 1e-30 && note("(b) $(r.set) $(p.id): package polynomials rel " *
                                "$(Float64(worst_r))")
        worst_s > 1e-9 && note("(b) $(r.set) $(p.id): package samples rel $worst_s")
    end
    # (d) rational exponents
    halves = any(!isinteger(e) for P in (L.H, L.Q, L.D) for mo in keys(P.terms)
                 for (_, e) in mo)
    # (e) units of every printed term, and the normalizer rule
    u = var_units(m, L, sol)
    uD = unique(mono_units(mo, u) for mo in keys(L.D.terms))
    P = S.printed(m, L)
    rule, lambda = S.normalizer(m, L)
    Dc = S.split_conc(L.D, L.vars)
    punits(p) = p === nothing ? (0 // 1, 0 // 1) : mono_units(first(keys(p.terms)), u)
    u_apart = reduce((a, b) -> a .+ b, (e .* punits(S.content(Dc[mo])[3]) for (mo, e) in lambda);
                     init = (0 // 1, 0 // 1))
    tu(cm, x) = mono_units(S.mmul(cm, x[2]), u) .+ punits(x[3]) .- (x[4] ? u_apart : (0, 0))
    dimless = length(uD) == 1 && all(tu(g...) == (0, 0) for g in P.d)
    S_units = mono_units(S.conc_mono(L.vars, m.substrates), u)
    rate_units = all(tu(g...) .+ S_units == (0, -1) for g in P.q)
    rule === :coefficient || dimless && rate_units ||
        note("(e) $(r.set) $(p.id): $rule terms not dimensionless (D units $(uD))")
    # printed laws evaluated against the law
    pars = random_params(m, rng)
    conc = Oracle.concs_at_ratio(p, rng, pars.Keq * 7)
    cnt = (; (x => conc[x] for x in S.metabolites(m))...)
    vl = lawval(L, conc, pars)
    full, coef = S.rate_equation_string(m, L, sol), S.coefficient_form(m, L)
    e_full = abs(eval_text(full, pars, cnt) / pars.E_total - vl) / abs(vl)
    e_coef = abs(eval_text(coef, pars, cnt; legend_first = true) / pars.E_total - vl) / abs(vl)
    max(e_full, e_coef) > 1e-30 && note("print $(r.set) $(p.id): full $(Float64(e_full)) " *
                                         "coefficient form $(Float64(e_coef))")
    push!(rows, (set = r.set, id = p.id, worst_a = Float64(worst_a), worst_s,
                 worst_r = Float64(worst_r), nsamp, halves, rule, dimless, rate_units,
                 e_full = Float64(e_full), e_coef = Float64(e_coef),
                 qpos = all(>(0), values(L.Q.terms)), nQ = S.nterms(L.Q), nD = S.nterms(L.D),
                 pkg_num = length(p.pkg_num), pkg_den = length(p.pkg_den),
                 sqrt = occursin("sqrt(", full), squared = occursin("²", coef),
                 keqQ = any(mo -> S.mexp(mo, L.vars.index[:Keq]) != 0, keys(L.Q.terms)),
                 big = L isa S.Law{BigInt},
                 maxc = maximum(abs, [collect(values(L.Q.terms)); collect(values(L.D.terms))])))
end
serialize(joinpath(@__DIR__, "data", "engine_rows.jls"), rows)
@printf("%d plain records, %d laws, %.0f s; %d notes\n", length(records), length(rows),
        time() - t0, length(fails))
fin(x) = filter(!isnan, x)
@printf("(a) oracle: %d laws x 3 points, worst relative error %.1e, above 1e-30: %d\n",
        length(rows), maximum(r.worst_a for r in rows), count(r -> r.worst_a > 1e-30, rows))
pk = filter(r -> !isnan(r.worst_r), rows)
@printf("(b) package polynomials: %d laws x 2 random points, worst %.1e; samples: %d, \
        worst %.1e\n", length(pk), maximum(r.worst_r for r in pk), sum(r.nsamp for r in pk),
        maximum(r.worst_s for r in pk))
@printf("    Q positive in %d of %d; Keq in Q in %d (%d of the 806); BigInt laws %d; \
        largest coefficient %d\n", count(r -> r.qpos, rows), length(rows),
        count(r -> r.keqQ, rows), count(r -> r.keqQ && r.set in ("specs", "handwritten",
        "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp", "expand_bibi"), rows),
        count(r -> r.big, rows), maximum(r -> r.maxc, rows))
hv = filter(r -> r.halves, rows)
println("(d) laws with a non-integer exponent: ", length(hv), ": ",
        join(("$(r.set) $(r.id) (sqrt printed: $(r.sqrt))" for r in hv), ", "))
for rule in (:constant, :cleland, :affine, :coefficient)
    rs = filter(r -> r.rule === rule, rows)
    sh = filter(r -> r.set in ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi",
                               "init_bibi_pp", "expand_bibi"), rs)
    @printf("(e) %-12s %4d laws (%3d of the 806): dimensionless denominator terms %d, \
            numerator terms 1/time %d\n", rule, length(rs), length(sh),
            count(r -> r.dimless, rs), count(r -> r.rate_units, rs))
    rule === :coefficient && println("    ", join(("$(r.set) $(r.id)" for r in rs), ", "))
end
@printf("print: full law worst %.1e, coefficient form worst %.1e (BigFloat, %d laws); \
        coefficient names with a superscript: %d\n", maximum(r.e_full for r in rows),
        maximum(r.e_coef for r in rows), length(rows), count(r -> r.squared, rows))
