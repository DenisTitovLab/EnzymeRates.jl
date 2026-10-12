# ABOUTME: Secondary verifier searches: textbook with Q > 0 (same-Q and zero-product saturated rates),
# ABOUTME: monotonicity along rays and fixed-B slices, the ordered-release fixture g, and d2 directions.

include(joinpath(@__DIR__, "vlib.jl"))
using .VL, Printf, Random
using .VL: QQ, toQ, qpow10, rate_exact, MECHS
include(joinpath(@__DIR__, "vopt.jl"))

const TB = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 60.0   # seconds per DE run
const mT = MECHS["T"]
const KAT = VL.build_ka(mT)
const TREE = (1, 2, 3, 5, 6, 7)

function consts_T!(kf, kr, x)
    lK1, lK2, lK3 = x[1] - x[10], x[2] - x[11], x[3] - x[12]
    lK5, lK6, lK7 = x[5] - x[13], x[6] - x[14], x[7] - x[15]
    lKeq = x[16]
    lkr = (x[10], x[11], x[12], x[4] - (lK1 + lK3 - lK2), x[13], x[14], x[15],
           x[8] - (lKeq - lK1 - lK3 - lK5 - lK6), x[9] - (lKeq - lK1 - lK3 - lK5 - lK7))
    for i in 1:9
        kf[i] = exp10(x[i])
        kr[i] = exp10(lkr[i])
    end
    lkr
end
violation(l, klo, khi) = max(0.0, l - khi) + max(0.0, klo - l)

struct WT
    kf::Vector{Float64}
    kr::Vector{Float64}
    fr::Vector{Float64}
    br::Vector{Float64}
    w::Vector{Float64}
end
WT() = WT(zeros(9), zeros(9), zeros(9), zeros(9), zeros(7))

"""Textbook turnover at P = 0 (P release flux, no subtraction) for any A, B, Q."""
function vT!(ws::WT, A, B, Q)
    kf, kr, fr, br = ws.kf, ws.kr, ws.fr, ws.br
    fr[1] = kf[1] * A; br[1] = kr[1]
    fr[2] = kf[2] * B; br[2] = kr[2]
    fr[3] = kf[3] * B; br[3] = kr[3]
    fr[4] = kf[4] * A; br[4] = kr[4]
    fr[5] = kf[5];     br[5] = kr[5]
    fr[6] = kf[6];     br[6] = 0.0
    fr[7] = kf[7];     br[7] = kr[7] * Q
    fr[8] = kf[8];     br[8] = kr[8] * Q
    fr[9] = kf[9];     br[9] = 0.0
    VL.ka_weights!(ws.w, KAT, fr, br)
    VL.ka_release(mT, ws.w, fr, 3)
end

function setT!(ws, x, klo, khi)
    lkr = consts_T!(ws.kf, ws.kr, x)
    viol = 0.0
    for i in 1:9
        viol += violation(x[i], klo, khi) + violation(lkr[i], klo, khi)
    end
    viol
end

mutable struct Best
    val::Float64
    x::Vector{Float64}
    n::Int
end
Best() = Best(-Inf, Float64[], 0)
upd!(b::Best, val, x) = (b.n += 1; val > b.val && ((b.val, b.x) = (val, collect(x))); nothing)
sx(b::Best, r) = isempty(b.x) ? "none" : string(round.(b.x[r]; digits = 3))

function run_de(label, F, lo, hi, seed; NP = 120, tmax = TB)
    gens, restarts = de!(F, lo, hi, MersenneTwister(seed); NP = NP, tmax = tmax)
    gens, restarts
end

function exact_T(x)
    setprecision(BigFloat, 256) do
        q(l) = toQ(exp10(big(l)))
        kf = [q(x[i]) for i in 1:9]
        kr = Vector{QQ}(undef, 9)
        for (j, i) in enumerate(TREE)
            kr[i] = q(x[9+j])
        end
        K(i) = kf[i] / kr[i]
        Keq = q(x[16])
        kr[4] = kf[4] / (K(1) * K(3) / K(2))
        kr[8] = kf[8] / (Keq / (K(1) * K(3) * K(5) * K(6)))
        kr[9] = kf[9] / (Keq / (K(1) * K(3) * K(5) * K(7)))
        kf, kr
    end
end

function main()
    t00 = time()
    ws = WT()

    # ---- S2: textbook, P = 0, Q > 0 -------------------------------------------------------
    println("== S2: textbook, P = 0, Q > 0 ==")
    for (label, klo, khi, cmax, seed) in (("rate constants 1e-9..1e9, A,B,Q 1e-12..1e12", -9.0, 9.0, 12.0, 1),
                                          ("rate constants 1e-6..1e6, A,B,Q 1e-8..1e8", -6.0, 6.0, 8.0, 2))
        bQ = Best()
        b0 = Best()
        lo = vcat(fill(klo, 15), -6.0, -cmax, -cmax, -cmax)
        hi = vcat(fill(khi, 15), 6.0, cmax, cmax, cmax)
        # (i) maximize v(A,B,Q)/v_sat(Q), the same-Q saturated rate
        FQ = x -> begin
            viol = setT!(ws, x, klo, khi)
            A, B, Q = exp10(x[17]), exp10(x[18]), exp10(x[19])
            v = vT!(ws, A, B, Q)
            R = v / VL.vsat_T(ws.kf, ws.kr, Q)
            viol == 0 && isfinite(R) && upd!(bQ, R, x)
            (R > 0 && isfinite(R) ? -log10(R) : 300.0) + 10viol
        end
        run_de(label, FQ, lo, hi, 10seed)
        # (ii) maximize v(A,B,Q)/v_sat(Q = 0)
        F0 = x -> begin
            viol = setT!(ws, x, klo, khi)
            A, B, Q = exp10(x[17]), exp10(x[18]), exp10(x[19])
            v = vT!(ws, A, B, Q)
            R = v / VL.vsat_T(ws.kf, ws.kr, 0.0)
            viol == 0 && isfinite(R) && upd!(b0, R, x)
            (R > 0 && isfinite(R) ? -log10(R) : 300.0) + 10viol
        end
        run_de(label, F0, lo, hi, 10seed + 1)
        @printf("  %s\n    max v(A,B,Q)/v_sat(same Q) = %.17f at log10(A,B,Q) = %s (%d evals)\n",
                label, bQ.val, sx(bQ, 17:19), bQ.n)
        @printf("    max v(A,B,Q)/v_sat(Q = 0)   = %.6e at log10(A,B,Q) = %s; analytic box sup 0.171573*10^(khi-klo) = %.6e\n",
                b0.val, sx(b0, 17:19), 0.17157287525 * exp10(khi - klo))
        isempty(bQ.x) || isempty(b0.x) || let
        kf, kr = exact_T(bQ.x)
        vx = rate_exact(mT, kf, kr, (toQ(exp10(big(bQ.x[17]))), toQ(exp10(big(bQ.x[18]))), zero(QQ),
                                     toQ(exp10(big(bQ.x[19])))))[1]
        @printf("    exact recheck of the same-Q best: 1 - v/v_sat(Q) = %.4e\n",
                Float64(BigFloat(1 - vx / VL.vsat_T(kf, kr, toQ(exp10(big(bQ.x[19])))))))
        kf, kr = exact_T(b0.x)
        vx = rate_exact(mT, kf, kr, (toQ(exp10(big(b0.x[17]))), toQ(exp10(big(b0.x[18]))), zero(QQ),
                                     toQ(exp10(big(b0.x[19])))))[1]
        @printf("    exact recheck of the zero-product best: v/v_sat(0) = %.10e\n",
                Float64(BigFloat(vx / VL.vsat_T(kf, kr))))
        end
        flush(stdout)
    end

    # ---- monotonicity in A at fixed B with Q > 0 (counterexample to "rises monotonically") ----
    println("\n== S2b: fixed-B slices in A at P = 0, Q > 0: largest drop v(A1)/v(A2), A1 < A2 ==")
    let b = Best(), klo = -6.0, khi = 6.0
        grid = collect(-8.0:0.25:8.0)
        vs = zeros(length(grid))
        lo = vcat(fill(klo, 15), -3.0, -8.0, -8.0)
        hi = vcat(fill(khi, 15), 3.0, 8.0, 8.0)
        F = x -> begin
            viol = setT!(ws, x, klo, khi)
            B, Q = exp10(x[17]), exp10(x[18])
            for (i, a) in enumerate(grid)
                vs[i] = vT!(ws, exp10(a), B, Q)
            end
            drop = 1.0
            mr = vs[end]
            for i in length(grid)-1:-1:1
                drop = max(drop, vs[i] / mr)
                mr = min(mr, vs[i])
            end
            viol == 0 && isfinite(drop) && upd!(b, drop, x)
            -log10(drop) + 10viol
        end
        run_de("slice", F, lo, hi, 31; tmax = TB / 2)
        @printf("  largest drop along A (B, Q fixed, A in 1e-8..1e8) = %.6e at log10(B,Q) = %s\n",
                b.val, sx(b, 17:18))
        flush(stdout)
    end

    # ---- S3: rays (A/B fixed) at P = Q = 0 ----------------------------------------------------
    println("\n== S3: textbook, P = Q = 0, monotonicity along rays A = r B and along fixed-B slices ==")
    for (label, klo, khi) in (("k 1e-9..1e9", -9.0, 9.0), ("k 1e-6..1e6", -6.0, 6.0))
        b = Best()
        tgrid = collect(-14.0:0.25:14.0)
        vs = zeros(length(tgrid))
        lo = vcat(fill(klo, 15), -6.0, -12.0)
        hi = vcat(fill(khi, 15), 6.0, 12.0)
        F = x -> begin
            viol = setT!(ws, x, klo, khi)
            ρ = x[17]
            for (i, t) in enumerate(tgrid)
                vs[i] = vT!(ws, exp10(t + ρ / 2), exp10(t - ρ / 2), 0.0)
            end
            drop = 1.0
            mr = vs[end]
            for i in length(tgrid)-1:-1:1
                drop = max(drop, vs[i] / mr)
                mr = min(mr, vs[i])
            end
            viol == 0 && isfinite(drop) && upd!(b, drop, x)
            -log10(drop) + 10viol
        end
        run_de(label, F, lo, hi, 41; tmax = TB)
        @printf("  rays, %s: largest drop v(t1)/v(t2), t1 < t2, along A = 10^(t+ρ/2), B = 10^(t-ρ/2) = %.17f at ρ = %.3f (%d evals)\n",
                label, b.val, isempty(b.x) ? NaN : b.x[17], b.n)
        if b.val > 1 + 1e-9
            kf, kr = exact_T(b.x)
            ρ = b.x[17]
            vv = [rate_exact(mT, kf, kr, (toQ(exp10(big(t + ρ / 2))), toQ(exp10(big(t - ρ / 2))),
                                          zero(QQ), zero(QQ)))[1] for t in tgrid]
            d = maximum(Float64(vv[i] / vv[j]) for i in 1:length(tgrid) for j in i+1:length(tgrid))
            @printf("    exact recheck of the drop: %.17f\n", d)
        end
        flush(stdout)
    end
    let b = Best(), bp = Best(), klo = -6.0, khi = 6.0
        grid = collect(-8.0:0.25:8.0)
        vs = zeros(length(grid))
        lo = vcat(fill(klo, 15), -3.0, -8.0)
        hi = vcat(fill(khi, 15), 3.0, 8.0)
        F = x -> begin
            viol = setT!(ws, x, klo, khi)
            B = exp10(x[17])
            for (i, a) in enumerate(grid)
                vs[i] = vT!(ws, exp10(a), B, 0.0)
            end
            drop = 1.0
            mr = vs[end]
            for i in length(grid)-1:-1:1
                drop = max(drop, vs[i] / mr)
                mr = min(mr, vs[i])
            end
            if viol == 0 && isfinite(drop)
                upd!(b, drop, x)
                upd!(bp, maximum(vs) / VL.vsat_T(ws.kf, ws.kr), x)
            end
            -log10(drop) + 10viol
        end
        run_de("slices", F, lo, hi, 51; tmax = TB / 2)
        @printf("  fixed-B slices (k 1e-6..1e6, A in 1e-8..1e8): largest drop = %.6e at log10 B = %.3f; largest slice max / v_sat seen = %.17f\n",
                b.val, isempty(b.x) ? NaN : b.x[17], bp.val)
        flush(stdout)
    end

    # ---- S4: package fixture g (random binding, ordered release) with one product present ------
    println("\n== S4: g (random binding, ordered release), one product present, R = v / v_sat(P = Q = 0) ==")
    let mg = MECHS["g"], kag = VL.build_ka(mg)
        kf = zeros(7); kr = zeros(7); fr = zeros(7); br = zeros(7); w = zeros(6)
        # x = [lkf 1..7; lkr of 1 2 3 5 6; lKeq; lA; lB; lProduct]
        function setg!(x, klo, khi)
            lK1, lK2, lK3, lK5, lK6 = x[1] - x[8], x[2] - x[9], x[3] - x[10], x[5] - x[11],
                                      x[6] - x[12]
            lkr = (x[8], x[9], x[10], x[4] - (lK1 + lK3 - lK2), x[11], x[12],
                   x[7] - (x[13] - lK1 - lK3 - lK5 - lK6))
            viol = 0.0
            for i in 1:7
                kf[i] = exp10(x[i])
                kr[i] = exp10(lkr[i])
                viol += violation(x[i], klo, khi) + violation(lkr[i], klo, khi)
            end
            viol
        end
        for (prod, label) in ((4, "Q present, P = 0"), (3, "P present, Q = 0"))
            b = Best()
            klo, khi = -9.0, 9.0
            lo = vcat(fill(klo, 12), -6.0, -12.0, -12.0, -12.0)
            hi = vcat(fill(khi, 12), 6.0, 12.0, 12.0, 12.0)
            F = x -> begin
                viol = setg!(x, klo, khi)
                A, B, C = exp10(x[14]), exp10(x[15]), exp10(x[16])
                P = prod == 3 ? C : 0.0
                Q = prod == 4 ? C : 0.0
                fr[1] = kf[1] * A; br[1] = kr[1]
                fr[2] = kf[2] * B; br[2] = kr[2]
                fr[3] = kf[3] * B; br[3] = kr[3]
                fr[4] = kf[4] * A; br[4] = kr[4]
                fr[5] = kf[5];     br[5] = kr[5]
                fr[6] = kf[6];     br[6] = kr[6] * P
                fr[7] = kf[7];     br[7] = kr[7] * Q
                VL.ka_weights!(w, kag, fr, br)
                v = VL.ka_release(mg, w, fr, prod == 4 ? 3 : 4)  # release flux of the absent product
                R = v / VL.vsat_g(kf, kr)
                viol == 0 && isfinite(R) && upd!(b, R, x)
                (R > 0 && isfinite(R) ? -log10(R) : 300.0) + 10viol
            end
            run_de(label, F, lo, hi, 60 + prod; tmax = TB)
            @printf("  %s: max R = %.17f at log10(A,B,product) = %s (%d evals)\n", label, b.val,
                    sx(b, 14:16), b.n)
            flush(stdout)
        end
    end

    # ---- S5: d2 (two ternary complexes, each releasing both products) along fixed directions ---
    println("\n== S5: d2, P = Q = 0: v(A,B) against the saturated rate along the same ray ==")
    let md = MECHS["d2"], kad = VL.build_ka(md)
        kf = zeros(6); kr = zeros(6); fr = zeros(6); br = zeros(6); w = zeros(5)
        b = Best(); bmax = Best()
        klo, khi = -6.0, 6.0
        # x = [lkf 1..6; lkr 1..4; lKeq; lA; lB]; K5 = Keq/(K1 K3), K6 = Keq/(K2 K4)
        lo = vcat(fill(klo, 10), -6.0, -8.0, -8.0)
        hi = vcat(fill(khi, 10), 6.0, 8.0, 8.0)
        F = x -> begin
            lK = (x[1] - x[7], x[2] - x[8], x[3] - x[9], x[4] - x[10])
            lkr = (x[7], x[8], x[9], x[10], x[5] - (x[11] - lK[1] - lK[3]),
                   x[6] - (x[11] - lK[2] - lK[4]))
            viol = 0.0
            for i in 1:6
                kf[i] = exp10(x[i])
                kr[i] = exp10(lkr[i])
                viol += violation(x[i], klo, khi) + violation(lkr[i], klo, khi)
            end
            A, B = exp10(x[12]), exp10(x[13])
            fr[1] = kf[1] * A; br[1] = kr[1]
            fr[2] = kf[2] * B; br[2] = kr[2]
            fr[3] = kf[3] * B; br[3] = kr[3]
            fr[4] = kf[4] * A; br[4] = kr[4]
            fr[5] = kf[5];     br[5] = 0.0
            fr[6] = kf[6];     br[6] = 0.0
            VL.ka_weights!(w, kad, fr, br)
            v = VL.ka_release(md, w, fr, 3)
            p = kf[1] * A / (kf[1] * A + kf[2] * B)
            vray = 1 / (p / kf[5] + (1 - p) / kf[6])
            R = v / vray
            if viol == 0 && isfinite(R)
                upd!(b, R, x)
                upd!(bmax, v / max(kf[5], kf[6]), x)
            end
            (R > 0 && isfinite(R) ? -log10(R) : 300.0) + 10viol
        end
        run_de("d2", F, lo, hi, 71; tmax = TB)
        @printf("  max v(A,B)/v_sat(same ray) = %.6f at log10(A,B) = %s; max v/max(kc1,kc2) seen = %.17f\n",
                b.val, sx(b, 12:13), bmax.val)
        if b.val > 1 + 1e-9
            x = b.x
            setprecision(BigFloat, 256) do
                q(l) = toQ(exp10(big(l)))
                kfq = [q(x[i]) for i in 1:6]
                krq = Vector{QQ}(undef, 6)
                for i in 1:4
                    krq[i] = q(x[6+i])
                end
                K(i) = kfq[i] / krq[i]
                Keq = q(x[11])
                krq[5] = kfq[5] / (Keq / (K(1) * K(3)))
                krq[6] = kfq[6] / (Keq / (K(2) * K(4)))
                A, B = q(x[12]), q(x[13])
                v = rate_exact(md, kfq, krq, (A, B, zero(QQ), zero(QQ)))[1]
                vlim = rate_exact(md, kfq, krq, (A * qpow10(150), B * qpow10(150), zero(QQ), zero(QQ)))[1]
                @printf("    exact: v = %.10e, same-ray limit (x1e150) = %.10e, ratio %.10f; kc1 = %.3e, kc2 = %.3e, k1 = %.3e, k2 = %.3e\n",
                        Float64(v), Float64(vlim), Float64(v / vlim), Float64(kfq[5]), Float64(kfq[6]),
                        Float64(kfq[1]), Float64(kfq[2]))
            end
        end
        flush(stdout)
    end
    @printf("\n(total %.0f s)\n", time() - t00)
end

main()
