# ABOUTME: P2 correctness on every plain mechanism: detailed balance of my constants at Γ = Keq, and my
# ABOUTME: dependents vs the package's bare constants and Reduced rate at the same underlying constants.
include(joinpath(@__DIR__, "..", "common", "Oracle.jl"))
include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))
using Random, Printf
const T = Thermo

logu(rng) = Oracle.logu(rng)
rel(a, b) = abs(a - b) / max(abs(b), big(1e-300))

setprecision(BigFloat, 512)
rng = Xoshiro(7)
nm = 0; db_worst = big(0); db_fail = 0
cmp_worst_bare = big(0); cmp_worst_v = big(0); cmp_worst_poly = big(0); cmp_fail = 0; nsamp = 0
worst_id = ""
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.allosteric && continue
    global nm += 1
    sol = T.solve(m)
    # detailed balance: random independents and Keq, v(Γ = Keq) against v(Γ = Keq/10)
    for _ in 1:2
        keq = logu(rng)
        iv = Dict(i => logu(rng) for i in sol.indep)
        k = T.stored_constants(sol, T.all_values(sol, iv, keq))
        v0 = Oracle.steady_state_v(m, k, Oracle.concs_at_ratio(m, rng, keq))
        v1 = Oracle.steady_state_v(m, k, Oracle.concs_at_ratio(m, rng, keq / 10))
        r = abs(v0) / abs(v1)
        r > db_worst && (global db_worst = r)
        r < 1e-30 || (global db_fail += 1; println("DB FAIL ", m.id, " ratio ", Float64(r)))
    end
    # the package's own constants
    for s in m.samples
        global nsamp += 1
        vals = Dict{Symbol, BigFloat}(p => BigFloat(v) for (p, v) in s.params)
        for (d, mono) in m.dependent_monomials
            vals[d] = prod((vals[x]^BigFloat(e) for (x, e) in mono); init = big(1))
        end
        kp = Oracle.constants_from_names(m, vals; rename = Dict(m.pass2_rename))
        bare = T.bare_values(sol, kp)
        keq = vals[:Keq]
        x = T.all_values(sol, Dict(i => bare[i] for i in sol.indep), keq)
        db = maximum(rel(x[i], bare[i]) for i in keys(bare))
        c = Dict{Symbol, BigFloat}(k => BigFloat(v) for (k, v) in s.concs)
        v = Oracle.steady_state_v(m, T.stored_constants(sol, x), c)
        dv = rel(v, BigFloat(s.v))
        allv = merge(vals, c)
        vp = Oracle.eval_poly(m.pkg_num, allv) / Oracle.eval_poly(m.pkg_den, allv)
        dp = rel(v, vp)
        if db > cmp_worst_bare; global cmp_worst_bare = db; global worst_id = m.id; end
        dv > cmp_worst_v && (global cmp_worst_v = dv)
        dp > cmp_worst_poly && (global cmp_worst_poly = dp)
        (db < 1e-100 && dv < 1e-12 && dp < 1e-25) || (global cmp_fail += 1;
            println("CMP FAIL ", m.id, " bare ", Float64(db), " v ", Float64(dv), " poly ", Float64(dp)))
    end
end
@printf("plain mechanisms %d\n", nm)
@printf("detailed balance: worst |v(Keq)|/|v(Keq/10)| = %.2e over %d draws, failures (>= 1e-30) %d\n",
        Float64(db_worst), 2nm, db_fail)
@printf("package samples %d: worst bare-constant rel diff %.2e (%s), worst v vs package Float64 v %.2e, worst v vs exact package polynomials %.2e, failures %d\n",
        nsamp, Float64(cmp_worst_bare), worst_id, Float64(cmp_worst_v), Float64(cmp_worst_poly), cmp_fail)
