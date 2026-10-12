# ABOUTME: d1 (two ternary complexes sharing E(P,Q)) against its saturated rate along the same ray,
# ABOUTME: and an explicit textbook slice that falls with A at fixed B with Q = 1, P = 0 (exact).

isdefined(@__MODULE__, :VL) || include(joinpath(@__DIR__, "vlib.jl"))
using .VL, Printf, Random
using .VL: QQ, toQ, qpow10, rate_exact, MECHS
isdefined(@__MODULE__, :de!) || include(joinpath(@__DIR__, "vopt.jl"))

const TD = 60.0
viold(l, lo, hi) = max(0.0, l - hi) + max(0.0, lo - l)

# d1 steps: 1 E+A->EA 2 E+B->EB 3 EA+B->X1 4 EB+A->X2 5 X1->EPQ 6 X2->EPQ 7 EPQ->EQ+P
# 8 EPQ->EP+Q 9 EQ->E+Q 10 EP->E+P.  x = [lkf 1..10; lkr of 1 2 3 4 5 7 8; lKeq; lA; lB]
const D1TREE = (1, 2, 3, 4, 5, 7, 8)
function d1consts(x)
    lkf = ntuple(i -> x[i], 10)
    lkrt = ntuple(j -> x[10+j], 7)
    lK = ntuple(j -> lkf[D1TREE[j]] - lkrt[j], 7)   # K1 K2 K3 K4 K5 K7 K8
    lKeq = x[18]
    lK6 = lK[1] + lK[3] + lK[5] - lK[2] - lK[4]
    lK9 = lKeq - lK[1] - lK[3] - lK[5] - lK[6]
    lK10 = lKeq - lK[1] - lK[3] - lK[5] - lK[7]
    lkr = (lkrt[1], lkrt[2], lkrt[3], lkrt[4], lkrt[5], lkf[6] - lK6, lkrt[6], lkrt[7],
           lkf[9] - lK9, lkf[10] - lK10)
    lkf, lkr
end
"""Saturated d1 turnover along a ray with p = k1 A/(k1 A + k2 B) of enzymes binding A first."""
function vray_d1(kf, kr, p)
    krel = kf[7] + kf[8]
    1 / (p / kf[5] + (1 - p) / kf[6] + (1 + kr[5] / kf[5] + kr[6] / kf[6]) / krel +
         (kf[7] / kf[9] + kf[8] / kf[10]) / krel)
end

let m = MECHS["d1"], ka = VL.build_ka(m)
    println("== d1, P = Q = 0: v(A,B) against the saturated rate along the same ray ==")
    kf = zeros(10); kr = zeros(10); fr = zeros(10); br = zeros(10); w = zeros(8)
    # closed form against the exact limit on random rays
    let rng = MersenneTwister(3), worst = 0.0
        for _ in 1:30
            x = vcat(12 .* rand(rng, 10) .- 6, zeros(7), 6 * rand(rng) - 3, 4 * rand(rng) - 2, 0.0)
            for j in 1:7
                x[10+j] = x[D1TREE[j]] - (6 * rand(rng) - 3)
            end
            lkf, lkr = d1consts(x)
            setprecision(BigFloat, 256) do
                kfq = [toQ(exp10(big(l))) for l in lkf]
                krq = [toQ(exp10(big(l))) for l in lkr]
                A = toQ(exp10(big(x[19]))) * qpow10(150)
                B = qpow10(150)
                vx = rate_exact(m, kfq, krq, (A, B, zero(QQ), zero(QQ)))[1]
                p = kfq[1] * A / (kfq[1] * A + kfq[2] * B)
                an = vray_d1(kfq, krq, p)
                worst = max(worst, Float64(abs(vx - an) / an))
            end
        end
        @printf("   closed-form same-ray v_sat vs exact at 1e150 scale, 30 random rays: max rel diff %.1e\n", worst)
    end
    best = Ref(-Inf); bestx = Ref(Float64[]); bestsup = Ref(-Inf)
    klo, khi = -6.0, 6.0
    lo = vcat(fill(klo, 17), -3.0, -8.0, -8.0)
    hi = vcat(fill(khi, 17), 3.0, 8.0, 8.0)
    F = x -> begin
        lkf, lkr = d1consts(x)
        viol = 0.0
        for i in 1:10
            kf[i] = exp10(lkf[i])
            kr[i] = exp10(lkr[i])
            viol += viold(lkf[i], klo, khi) + viold(lkr[i], klo, khi)
        end
        A, B = exp10(x[19]), exp10(x[20])
        fr[1] = kf[1] * A; br[1] = kr[1]
        fr[2] = kf[2] * B; br[2] = kr[2]
        fr[3] = kf[3] * B; br[3] = kr[3]
        fr[4] = kf[4] * A; br[4] = kr[4]
        for i in 5:10
            fr[i] = kf[i]
        end
        br[5] = kr[5]; br[6] = kr[6]; br[7] = 0.0; br[8] = 0.0; br[9] = 0.0; br[10] = 0.0
        VL.ka_weights!(w, ka, fr, br)
        v = VL.ka_release(m, w, fr, 3)
        p = kf[1] * A / (kf[1] * A + kf[2] * B)
        R = v / vray_d1(kf, kr, p)
        if viol == 0 && isfinite(R)
            R > best[] && (best[] = R; bestx[] = collect(x))
            bestsup[] = max(bestsup[], v / max(vray_d1(kf, kr, 0.0), vray_d1(kf, kr, 1.0)))
        end
        (R > 0 && isfinite(R) ? -log10(R) : 300.0) + 10viol
    end
    de!(F, lo, hi, MersenneTwister(81); NP = 150, tmax = TD)
    x = bestx[]
    @printf("   max v/v_sat(same ray) = %.6e at log10(A,B) = %s; max v / sup over directions of v_sat seen = %.17f\n",
            best[], string(round.(x[19:20]; digits = 3)), bestsup[])
    lkf, lkr = d1consts(x)
    setprecision(BigFloat, 256) do
        kfq = [toQ(exp10(big(l))) for l in lkf]
        krq = [toQ(exp10(big(l))) for l in lkr]
        dev, _ = VL.thermo_check(m, BigFloat.(kfq), BigFloat.(krq), x[18])
        A, B = toQ(exp10(big(x[19]))), toQ(exp10(big(x[20])))
        println("   exact profile along the ray (A, B) = t (A0, B0), v / v_sat(same ray):")
        vlim = rate_exact(m, kfq, krq, (A * qpow10(150), B * qpow10(150), zero(QQ), zero(QQ)))[1]
        for lt in -6:2:12
            v = rate_exact(m, kfq, krq, (A * toQ(exp10(big(lt))), B * toQ(exp10(big(lt))), zero(QQ), zero(QQ)))[1]
            @printf("      t = 1e%-3d  v = %.6e   v/v_sat(ray) = %.6e\n", lt, Float64(v), Float64(v / vlim))
        end
        @printf("   thermo dev %.1e; kc1 = %.3e kc2 = %.3e k1 = %.3e k2 = %.3e; same-ray v_sat = %.6e, sup-direction v_sat = %.6e\n",
                Float64(dev), Float64(kfq[5]), Float64(kfq[6]), Float64(kfq[1]), Float64(kfq[2]),
                Float64(vlim), Float64(max(vray_d1(kfq, krq, 0), vray_d1(kfq, krq, 1))))
    end
end

let m = MECHS["T"], ka = VL.build_ka(m)
    println("\n== textbook, P = 0, Q = 1: a fixed-B slice that falls with A ==")
    kf = zeros(9); kr = zeros(9); fr = zeros(9); br = zeros(9); w = zeros(7)
    TREE = (1, 2, 3, 5, 6, 7)
    function consts!(x)
        lK = ntuple(j -> x[TREE[j]] - x[9+j], 6)
        lkr = (x[10], x[11], x[12], x[4] - (lK[1] + lK[3] - lK[2]), x[13], x[14], x[15],
               x[8] - (x[16] - lK[1] - lK[3] - lK[4] - lK[5]),
               x[9] - (x[16] - lK[1] - lK[3] - lK[4] - lK[6]))
        lkr
    end
    vT(A, B, Q) = begin
        fr[1] = kf[1] * A; br[1] = kr[1]
        fr[2] = kf[2] * B; br[2] = kr[2]
        fr[3] = kf[3] * B; br[3] = kr[3]
        fr[4] = kf[4] * A; br[4] = kr[4]
        fr[5] = kf[5]; br[5] = kr[5]
        fr[6] = kf[6]; br[6] = 0.0
        fr[7] = kf[7]; br[7] = kr[7] * Q
        fr[8] = kf[8]; br[8] = kr[8] * Q
        fr[9] = kf[9]; br[9] = 0.0
        VL.ka_weights!(w, ka, fr, br)
        VL.ka_release(m, w, fr, 3)
    end
    grid = collect(-4.0:0.5:8.0)
    best = Ref(-Inf); bestx = Ref(Float64[])
    lo = vcat(fill(-4.0, 15), -2.0, 0.0)
    hi = vcat(fill(4.0, 15), 2.0, 6.0)
    F = x -> begin
        lkr = consts!(x)
        viol = 0.0
        for i in 1:9
            kf[i] = exp10(x[i])
            kr[i] = exp10(lkr[i])
            viol += viold(x[i], -4.0, 4.0) + viold(lkr[i], -4.0, 4.0)
        end
        B = exp10(x[17])
        vs = [vT(exp10(a), B, 1.0) for a in grid]
        drop = maximum(vs[i] / minimum(vs[i+1:end]) for i in 1:length(vs)-1)
        viol == 0 && isfinite(drop) && drop > best[] && (best[] = drop; bestx[] = collect(x))
        -log10(drop) + 10viol
    end
    de!(F, lo, hi, MersenneTwister(91); NP = 120, tmax = TD / 2)
    x = bestx[]
    lkr = consts!(x)
    setprecision(BigFloat, 256) do
        kfq = [toQ(exp10(big(x[i]))) for i in 1:9]
        krq = [toQ(exp10(big(lkr[i]))) for i in 1:9]
        dev, _ = VL.thermo_check(m, BigFloat.(kfq), BigFloat.(krq), x[16])
        B = toQ(exp10(big(x[17])))
        vsQ = VL.vsat_T(kfq, krq, one(QQ))
        @printf("   constants in 1e-4..1e4, log10 Keq = %.3f, B = %.4g, Q = 1, P = 0 (thermo dev %.1e)\n",
                x[16], Float64(B), Float64(dev))
        @printf("   kf = %s\n   kr = %s\n", string(round.(Float64.(kfq); sigdigits = 4)),
                string(round.(Float64.(krq); sigdigits = 4)))
        for la in -4:1:8
            v = rate_exact(m, kfq, krq, (qpow10(la), B, zero(QQ), one(QQ)))[1]
            @printf("      A = 1e%-3d v = %.6e   v / v_sat(Q = 1) = %.6e\n", la, Float64(v), Float64(v / vsQ))
        end
    end
end
