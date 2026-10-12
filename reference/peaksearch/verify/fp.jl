# ABOUTME: Exact (rational) checks of the first-passage argument for the textbook mechanism and the
# ABOUTME: parallel-route variant d1, the saturated-rate closed forms, and the P = 0, Q > 0 bound.

include(joinpath(@__DIR__, "vlib.jl"))
using .VL, Printf, Random
using .VL: QQ, toQ, qpow10, rate_exact, pseudo, expected_times, MECHS

r10(rng, lo, hi) = toQ(exp10(lo + (hi - lo) * rand(rng)))

"""Random thermodynamically consistent textbook constants, exact: Keq, kf[1:9], K free on the
tree steps 1 2 3 5 6 7; K4 = K1 K3/K2, K8 = Keq/(K1 K3 K5 K6), K9 = Keq/(K1 K3 K5 K7)."""
function random_T(rng; lo = -6, hi = 6)
    Keq = r10(rng, -6, 6)
    kf = [r10(rng, lo, hi) for _ in 1:9]
    K = Vector{QQ}(undef, 9)
    for i in (1, 2, 3, 5, 6, 7)
        K[i] = r10(rng, lo, hi)
    end
    K[4] = K[1] * K[3] / K[2]
    K[8] = Keq / (K[1] * K[3] * K[5] * K[6])
    K[9] = Keq / (K[1] * K[3] * K[5] * K[7])
    kf, kf ./ K, K, Keq
end

const mT = MECHS["T"]
t0 = time()
println("=== first-passage decomposition, textbook, P = Q = 0, exact rationals ===")
let rng = MersenneTwister(20261011), n = 400, ok_ein = 0, ok_ts = 0, ok_sum = 0, ok_bound = 0,
    ok_thermo = 0, minmarg_in = Inf, minmarg_s = Inf, maxR = -Inf
    for trial in 1:n
        kf, kr, K, Keq = random_T(rng)
        dev, _ = VL.thermo_check(mT, BigFloat.(kf), BigFloat.(kr), log10(BigFloat(Keq)))
        dev < 1e-60 && (ok_thermo += 1)
        A, B = r10(rng, -8, 8), r10(rng, -8, 8)
        c = (A, B, zero(QQ), zero(QQ))
        fr, br = pseudo(mT, kf, kr, c)
        # states 1 E, 2 EA, 3 EB, 4 EAB, 5 EPQ, 6 EP, 7 EQ
        Nin = expected_times(mT, fr, br, [1, 2, 3, 4], 1)
        Nin[4] == 1 / kf[5] && (ok_ein += 1)
        τin = sum(Nin)
        tS = sum(expected_times(mT, fr, br, [1, 2, 3, 4, 5], 5))
        πeq = [one(QQ), K[1] * A, K[2] * B, K[1] * K[3] * A * B, K[1] * K[3] * K[5] * A * B]
        krel = kf[6] + kf[7]
        tS_renewal = sum(πeq) / (πeq[5] * krel)
        tS == tS_renewal && (ok_ts += 1)
        τrel = (kf[6] / kf[8] + kf[7] / kf[9]) / krel
        v = rate_exact(mT, kf, kr, c)[1]
        1 / v == τin + tS + τrel && (ok_sum += 1)
        vs = VL.vsat_T(kf, kr)
        v < vs && (ok_bound += 1)
        minmarg_in = min(minmarg_in, Float64((τin - 1 / kf[5]) / τin))
        minmarg_s = min(minmarg_s, Float64((tS - (1 + kr[5] / kf[5]) / krel) / tS))
        maxR = max(maxR, Float64(v / vs))
    end
    @printf("%d random sets (kf, K in 1e-6..1e6, Keq in 1e-6..1e6, A, B in 1e-8..1e8)\n", n)
    @printf("   cycle products consistent (|dev| < 1e-60): %d/%d\n", ok_thermo, n)
    @printf("   mean time in E(A,B) before first catalysis == 1/k5 exactly: %d/%d\n", ok_ein, n)
    @printf("   t_S (fundamental matrix) == renewal formula exactly: %d/%d\n", ok_ts, n)
    @printf("   1/v == tau_in + t_S + tau_rel exactly: %d/%d\n", ok_sum, n)
    @printf("   v < v_sat (closed form) exactly: %d/%d ; max v/v_sat = %.17f\n", ok_bound, n, maxR)
    @printf("   min relative margins: (tau_in - 1/k5)/tau_in = %.3e, (t_S - t_S,sat)/t_S = %.3e\n",
            minmarg_in, minmarg_s)
end
@printf("(%.0f s)\n", time() - t0)

println("\n=== closed forms against the exact limit (A = B = 1e100, and 1e200/1e100 directions) ===")
let rng = MersenneTwister(7), worstT = 0.0, worstQ = 0.0, worstg = 0.0, mg = MECHS["g"]
    for trial in 1:100
        kf, kr, K, Keq = random_T(rng)
        an = VL.vsat_T(kf, kr)
        for (a, b) in ((100, 100), (200, 100), (100, 200))
            v = rate_exact(mT, kf, kr, (qpow10(a), qpow10(b), zero(QQ), zero(QQ)))[1]
            worstT = max(worstT, Float64(abs(v - an) / an))
        end
        Q = r10(rng, -8, 8)
        anQ = VL.vsat_T(kf, kr, Q)
        v = rate_exact(mT, kf, kr, (qpow10(150), qpow10(150), zero(QQ), Q))[1]
        worstQ = max(worstQ, Float64(abs(v - anQ) / anQ))
        # g: 1 E+A 2 E+B 3 EA+B 4 EB+A 5 EAB->EPQ 6 EPQ->EQ+P 7 EQ->E+Q; cycle K1K3 = K2K4,
        # catalytic K1 K3 K5 K6 K7 = Keq
        kfg = [r10(rng, -6, 6) for _ in 1:7]
        Kg = Vector{QQ}(undef, 7)
        for i in (1, 2, 3, 5, 6)
            Kg[i] = r10(rng, -6, 6)
        end
        Kg[4] = Kg[1] * Kg[3] / Kg[2]
        Kg[7] = Keq / (Kg[1] * Kg[3] * Kg[5] * Kg[6])
        krg = kfg ./ Kg
        ang = VL.vsat_g(kfg, krg)
        vg = rate_exact(mg, kfg, krg, (qpow10(150), qpow10(120), zero(QQ), zero(QQ)))[1]
        worstg = max(worstg, Float64(abs(vg - ang) / ang))
    end
    @printf("textbook v_sat closed form vs exact at 1e100..1e200 (3 directions), 100 sets: max rel diff %.1e\n",
            worstT)
    @printf("textbook v_sat(Q) closed form (P = 0) vs exact at A = B = 1e150, 100 sets: max rel diff %.1e\n",
            worstQ)
    @printf("g v_sat closed form vs exact at (1e150, 1e120), 100 sets: max rel diff %.1e\n", worstg)
end

println("\n=== textbook, P = 0, Q > 0: v(A,B,Q) <= v(inf,inf,Q)?  (exact, random sets) ===")
let rng = MersenneTwister(99), n = 400, nok = 0, maxR = -Inf
    for trial in 1:n
        kf, kr, K, Keq = random_T(rng)
        A, B, Q = r10(rng, -8, 8), r10(rng, -8, 8), r10(rng, -8, 8)
        v = rate_exact(mT, kf, kr, (A, B, zero(QQ), Q))[1]
        vsQ = VL.vsat_T(kf, kr, Q)
        v < vsQ && (nok += 1)
        maxR = max(maxR, Float64(v / vsQ))
    end
    @printf("   v < v_sat(Q): %d/%d ; max v/v_sat(Q) = %.17f\n", nok, n, maxR)
end

println("\n=== d1 (parallel ternary complexes, shared E(P,Q)), exact ===")
# states 1 E, 2 EA, 3 EB, 4 EAB(X1), 5 EstarAB(X2), 6 EPQ, 7 EP, 8 EQ
# steps 1 E+A 2 E+B 3 EA+B->X1 4 EB+A->X2 5 X1->EPQ 6 X2->EPQ 7 EPQ->EQ+P 8 EPQ->EP+Q 9 EQ->E+Q 10 EP->E+P
let rng = MersenneTwister(5), m = MECHS["d1"], n = 200, ok1 = 0, ok2 = 0, ok3 = 0, ok4 = 0,
    okth = 0, maxR = -Inf
    for trial in 1:n
        Keq = r10(rng, -6, 6)
        kf = [r10(rng, -6, 6) for _ in 1:10]
        K = Vector{QQ}(undef, 10)
        for i in (1, 2, 3, 4, 5, 7, 8)
            K[i] = r10(rng, -6, 6)
        end
        K[6] = K[1] * K[3] * K[5] / (K[2] * K[4])
        K[9] = Keq / (K[1] * K[3] * K[5] * K[7])
        K[10] = Keq / (K[1] * K[3] * K[5] * K[8])
        kr = kf ./ K
        dev, _ = VL.thermo_check(m, BigFloat.(kf), BigFloat.(kr), log10(BigFloat(Keq)))
        dev < 1e-60 && (okth += 1)
        A, B = r10(rng, -8, 8), r10(rng, -8, 8)
        c = (A, B, zero(QQ), zero(QQ))
        fr, br = pseudo(m, kf, kr, c)
        Nin = expected_times(m, fr, br, [1, 2, 3, 4, 5], 1)
        Nin[4] * kf[5] + Nin[5] * kf[6] == 1 && (ok1 += 1)
        Nin[4] + Nin[5] >= 1 / max(kf[5], kf[6]) && (ok2 += 1)
        τin = sum(Nin)
        πeq = [one(QQ), K[1] * A, K[2] * B, K[1] * K[3] * A * B, K[2] * K[4] * A * B,
               K[1] * K[3] * K[5] * A * B]
        krel = kf[7] + kf[8]
        tS = sum(πeq) / (πeq[6] * krel)
        τrel = (kf[7] / kf[9] + kf[8] / kf[10]) / krel
        v = rate_exact(m, kf, kr, c)[1]
        1 / v == τin + tS + τrel && (ok3 += 1)
        vdirs = [rate_exact(m, kf, kr, (qpow10(a), qpow10(b), zero(QQ), zero(QQ)))[1]
                 for (a, b) in ((200, 100), (100, 200), (150, 150))]
        v < maximum(vdirs) && (ok4 += 1)
        maxR = max(maxR, Float64(v / maximum(vdirs)))
    end
    @printf("%d random sets: thermo ok %d; P(absorb) = N_X1 k5 + N_X2 k6 == 1: %d; time in X1+X2 >= 1/max(kc): %d; 1/v == tau_in + t_S + tau_rel: %d; v < max directional v_sat: %d ; max ratio %.17f\n",
            n, okth, ok1, ok2, ok3, ok4, maxR)
end
@printf("(total %.0f s)\n", time() - t0)
