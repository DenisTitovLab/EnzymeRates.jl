# ABOUTME: Re-evaluates every best point of the search agent's results_full.jls and its hand-built
# ABOUTME: examples in exact rational arithmetic: R, the saturated limit, thermodynamics, box limits.

include(joinpath(@__DIR__, "vlib.jl"))
using .VL, Serialization, Printf, Random
using .VL: QQ, toQ, qpow10, rate_exact, decode, thermo_check, boxviol, MECHS

setprecision(BigFloat, 512)
const SRC = joinpath(@__DIR__, "..")
const RES = deserialize(joinpath(SRC, "results_full.jls"))

bigc(lc) = toQ(exp10(big(lc)))
function concs(setting, lc)
    A, B = bigc(lc[1]), bigc(lc[2])
    z = zero(QQ)
    setting === :zero && return (A, B, z, z)
    setting === :Q && return (A, B, z, bigc(lc[3]))
    setting === :P && return (A, B, bigc(lc[3]), z)
    (A, B, bigc(lc[3]), bigc(lc[4]))
end

const DIRS48 = ((48, 48), (48, 24), (24, 48))
const DIRS200 = ((200, 200), (200, 100), (100, 200), (300, 300))

vsat_at(m, kfQ, krQ, a, b) =
    rate_exact(m, kfQ, krQ, (qpow10(a), qpow10(b), zero(QQ), zero(QQ)))[1]

f64(x) = Float64(x)

t0 = time()
println("=== agent best points (results_full.jls), exact rational re-evaluation ===")
println("columns: key setting log10Keq | R agent (Float64) | R exact vs max v_sat(1e200 dirs) | ",
        "R exact - 1 | thermo max |log10 dev| | box viol (1e-6..1e6) | v_sat: 1e48 vs 1e200 rel diff, ",
        "dir spread at 1e200")
for (k, r) in sort(collect(RES); by = p -> string(p[1]))
    key, setting, lK = k
    key == "b" && continue   # rapid-equilibrium variant, not rebuilt here
    set = setting === :zero_restricted ? :zero : setting
    m = MECHS[key]
    np = VL.nparams(m)
    θ = r.x[1:np]
    lc = r.x[np+1:end]
    kf, kr = decode(m, θ, lK)
    dev, rows = thermo_check(m, kf, kr, lK)
    bv = boxviol(kf, kr)
    kfQ, krQ = toQ.(kf), toQ.(kr)
    c = concs(set, lc)
    v = rate_exact(m, kfQ, krQ, c)[1]
    vs48 = [vsat_at(m, kfQ, krQ, a, b) for (a, b) in DIRS48]
    vs200 = [vsat_at(m, kfQ, krQ, a, b) for (a, b) in DIRS200]
    vsmax = maximum(vs200)
    R = v / vsmax
    reldiff48 = abs(maximum(vs48) - vsmax) / vsmax
    spread = (maximum(vs200) - minimum(vs200)) / vsmax
    @printf("%-3s %-15s %5.1f | %.10e | %.10e | %+.3e | %.1e | %.1e | %.1e %.1e\n",
            key, setting, lK, r.R, f64(R), f64(BigFloat(R) - 1), f64(dev), bv,
            f64(reldiff48), f64(spread))
    if key == "T"
        an = VL.vsat_T(kf, kr)
        rel = abs(BigFloat(vsmax) - an) / an
        @printf("      T analytic v_sat %.12e vs exact(1e200) rel diff %.1e\n", f64(an), f64(rel))
    end
    if r.R > 1 + 1e-9 && setting !== :zero_restricted
        @printf("      peak: log10 conc %s  v = %.10e  v_sat = %.10e\n",
                string(round.(lc; digits = 4)), f64(v), f64(vsmax))
    end
end
@printf("(%.0f s)\n", time() - t0)

println("\n=== full rate constants at the T, :Q, Keq = 1 optimum ===")
let r = RES[("T", :Q, 0.0)], m = MECHS["T"]
    np = VL.nparams(m)
    kf, kr = decode(m, r.x[1:np], 0.0)
    for (i, e) in enumerate(m.edges)
        lhs = join([String(m.states[e.u]); [String(VL.LIG[l]) for l in e.cons]], " + ")
        rhs = join([String(m.states[e.w]); [String(VL.LIG[l]) for l in e.rel]], " + ")
        @printf("   %-16s kf = %.6e  kr = %.6e\n", "$lhs <-> $rhs", f64(kf[i]), f64(kr[i]))
    end
    k5, km5, k6, k7, k8, k9 = kf[5], kr[5], kf[6], kf[7], kf[8], kf[9]
    T0 = 1 / k5 + (1 + km5 / k5 + k6 / k8 + k7 / k9) / (k6 + k7)
    Tinf = 1 / k5 + (1 + km5 / k5 + k6 / k8) / k6
    @printf("   fraction of turnovers through E(P) at saturation k7/(k6+k7) = %.4f\n", f64(k7 / (k6 + k7)))
    @printf("   T(Q=0)/T(Q->inf) at A,B->inf = %.6e\n", f64(T0 / Tinf))
    @printf("   condition k9 < k6/(1 + k-5/k5 + k6/k8): %.3e < %.3e\n", f64(k9),
            f64(k6 / (1 + km5 / k5 + k6 / k8)))
    kfQ, krQ = toQ.(kf), toQ.(kr)
    vs = VL.vsat_T(kf, kr)
    println("   v / v_sat at P = 0 along Q (A = B = 1e8 and A = B = 1e30):")
    for lq in (0.0, 2.0, 3.0, 3.5, 4.0, 4.035, 4.5, 5.0, 6.0, 8.0)
        v8 = rate_exact(m, kfQ, krQ, (qpow10(8), qpow10(8), zero(QQ), toQ(exp10(big(lq)))))[1]
        v30 = rate_exact(m, kfQ, krQ, (qpow10(30), qpow10(30), zero(QQ), toQ(exp10(big(lq)))))[1]
        @printf("      log10 Q = %5.3f   R(1e8) = %.8e   R(1e30) = %.8e   analytic sat-in-A,B = %.8e\n",
                lq, f64(v8 / vs), f64(v30 / vs), f64(VL.vsat_T(kf, kr, exp10(big(lq))) / vs))
    end
end

println("\n=== hand-built examples ===")
# T example: every constant 1 except E(P) <-> E + P at 1e-3 both ways, Keq = 1
let m = MECHS["T"]
    kfQ = [one(QQ) for _ in 1:9]
    krQ = [one(QQ) for _ in 1:9]
    kfQ[9] = krQ[9] = QQ(1, 1000)
    dev, _ = thermo_check(m, BigFloat.(kfQ), BigFloat.(krQ), 0.0)
    vsx = [vsat_at(m, kfQ, krQ, a, b) for (a, b) in DIRS200]
    an = VL.vsat_T(kfQ, krQ)
    @printf("T example: thermo dev %.1e; exact v_sat (1e200 dirs) = %s; analytic = %s = %.10e\n",
            f64(dev), string(f64.(vsx)), string(an), f64(an))
    println("   exact v_sat (all 1e200 dirs) == analytic to 1e-90? ", all(x -> abs(x - an) / an < QQ(1, big(10)^90), vsx))
    for lq in -4:8
        v = rate_exact(m, kfQ, krQ, (qpow10(8), qpow10(8), zero(QQ), qpow10(lq)))[1]
        @printf("   A = B = 1e8, P = 0, Q = 1e%-3d  v/v_sat = %.10f\n", lq, f64(v / an))
    end
    best, bq = -Inf, NaN
    for lq in 3.0:0.01:5.0
        v = rate_exact(m, kfQ, krQ, (qpow10(8), qpow10(8), zero(QQ), toQ(exp10(big(lq)))))[1]
        R = f64(v / an)
        R > best && ((best, bq) = (R, lq))
    end
    @printf("   max over log10 Q in [3,5] step 0.01 at A = B = 1e8: %.10f at log10 Q = %.2f\n", best, bq)
    Tinf = 1 + (1 + 1 + 1)
    @printf("   analytic T(0)/T(inf) = %s = %.6f\n", string(1 / an / Tinf), f64(1 / an / Tinf))
    for a in (0, 2, 4, 6, 8)
        v = rate_exact(m, kfQ, krQ, (qpow10(a), qpow10(a), zero(QQ), QQ(10)))[1]
        @printf("   A = B = 1e%d, Q = 10:  v/v_sat = %.6f\n", a, f64(v / an))
    end
end

# f example: every constant 1 except E(A,Q) <-> E(A) + Q at 1e-3 both ways, Keq = 1
let m = MECHS["f"]
    kfQ = [one(QQ) for _ in 1:11]
    krQ = [one(QQ) for _ in 1:11]
    kfQ[11] = krQ[11] = QQ(1, 1000)
    dev, _ = thermo_check(m, BigFloat.(kfQ), BigFloat.(krQ), 0.0)
    vsx = [vsat_at(m, kfQ, krQ, a, b) for (a, b) in ((48, 48), (100, 100), (200, 200),
                                                      (200, 100), (100, 200), (300, 300))]
    @printf("f example: thermo dev %.1e; exact v at A,B = (1e48,1e48),(1e100,1e100),(1e200,1e200),(1e200,1e100),(1e100,1e200),(1e300,1e300): %s\n",
            f64(dev), string(round.(f64.(vsx); sigdigits = 12)))
    vs = vsx[end]
    row = [f64(rate_exact(m, kfQ, krQ, (qpow10(a), qpow10(8), zero(QQ), zero(QQ)))[1] / vs)
           for a in -4:8]
    println("   v/v_sat at A = 1e-4..1e8 (x10), B = 1e8: ", string(round.(row; sigdigits = 5)))
    best, ba = -Inf, NaN
    for la in -1.0:0.01:1.5
        v = rate_exact(m, kfQ, krQ, (toQ(exp10(big(la))), qpow10(8), zero(QQ), zero(QQ)))[1]
        R = f64(v / vs)
        R > best && ((best, ba) = (R, la))
    end
    @printf("   max over log10 A in [-1,1.5] step 0.01 (B = 1e8): %.8f at log10 A = %.2f\n", best, ba)
    bestB = -Inf
    for la in -1.0:0.05:1.5, lb in 0.0:0.5:12.0
        v = rate_exact(m, kfQ, krQ, (toQ(exp10(big(la))), toQ(exp10(big(lb))), zero(QQ), zero(QQ)))[1]
        bestB = max(bestB, f64(v / vs))
    end
    @printf("   max over log10 A in [-1,1.5], log10 B in [0,12] (grid): %.8f\n", bestB)
end

# e example and random e sets: v -> 0 as A -> infinity
let m = MECHS["e"]
    kfQ = [one(QQ) for _ in 1:10]
    krQ = [one(QQ) for _ in 1:10]
    dev, _ = thermo_check(m, BigFloat.(kfQ), BigFloat.(krQ), 0.0)
    row = [f64(rate_exact(m, kfQ, krQ, (qpow10(a), qpow10(8), zero(QQ), zero(QQ)))[1])
           for a in -4:8]
    @printf("e example (all constants 1): thermo dev %.1e\n   v at A = 1e-4..1e8 (x10), B = 1e8: %s\n",
            f64(dev), string(round.(row; sigdigits = 4)))
    for (a, b) in ((20, 20), (50, 50), (100, 100), (100, 200), (200, 100))
        v = rate_exact(m, kfQ, krQ, (qpow10(a), qpow10(b), zero(QQ), zero(QQ)))[1]
        @printf("   A = 1e%d, B = 1e%d: v = %.6e, v*A = %.6f\n", a, b, f64(v), f64(v * qpow10(a)))
    end
end
let m = MECHS["e"], seed = 11
    rng = MersenneTwister(seed)
    worst = 0.0
    for trial in 1:20
        θ = vcat(12 .* rand(rng, 10) .- 6, 6 .* rand(rng, 7) .- 3)
        kf, kr = decode(m, θ, 0.0)
        kfQ, krQ = toQ.(kf), toQ.(kr)
        vA = [rate_exact(m, kfQ, krQ, (qpow10(a), qpow10(a), zero(QQ), zero(QQ)))[1] * qpow10(a)
              for a in (100, 200)]
        worst = max(worst, f64(abs(vA[1] - vA[2]) / vA[2]))
    end
    @printf("e, 20 random sets: v*A at A = B = 1e100 vs 1e200 max rel diff %.1e (v ~ c/A, so v_sat = 0)\n",
            worst)
end
@printf("(total %.0f s)\n", time() - t0)
