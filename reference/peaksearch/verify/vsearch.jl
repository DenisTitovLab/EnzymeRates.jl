# ABOUTME: Independent adversarial search for v(A, B) > v_sat in the textbook random-order bi-bi at
# ABOUTME: P = Q = 0: differential evolution, Nelder–Mead, bilevel sampling; exact recheck of the best.

include(joinpath(@__DIR__, "vlib.jl"))
using .VL, Printf, Random, Serialization
using .VL: QQ, toQ, qpow10, rate_exact, MECHS
include(joinpath(@__DIR__, "vopt.jl"))

const mT = MECHS["T"]
const KAT = VL.build_ka(mT)
const BUDGET = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1500.0   # seconds of search
const OUT = joinpath(@__DIR__, "vsearch_results.jls")

# x = [log10 kf1..kf9; log10 kr of steps 1 2 3 5 6 7; log10 Keq; log10 A; log10 B]
const TREE = (1, 2, 3, 5, 6, 7)

"""Log10 forward and reverse constants from x; reverse constants of steps 4, 8, 9 follow from
K4 = K1 K3/K2, K8 = Keq/(K1 K3 K5 K6), K9 = Keq/(K1 K3 K5 K7)."""
function logconsts(x)
    lkf = ntuple(i -> x[i], 9)
    lkrt = ntuple(j -> x[9+j], 6)
    lK = ntuple(j -> lkf[TREE[j]] - lkrt[j], 6)   # K1 K2 K3 K5 K6 K7
    lKeq = x[16]
    lkr4 = lkf[4] - (lK[1] + lK[3] - lK[2])
    lkr8 = lkf[8] - (lKeq - lK[1] - lK[3] - lK[4] - lK[5])
    lkr9 = lkf[9] - (lKeq - lK[1] - lK[3] - lK[4] - lK[6])
    lkr = (lkrt[1], lkrt[2], lkrt[3], lkr4, lkrt[4], lkrt[5], lkrt[6], lkr8, lkr9)
    lkf, lkr
end

mutable struct WS
    kf::Vector{Float64}
    kr::Vector{Float64}
    fr::Vector{Float64}
    br::Vector{Float64}
    w::Vector{Float64}
end
WS() = WS(zeros(9), zeros(9), zeros(9), zeros(9), zeros(7))

"""R = v(A, B)/v_sat at P = Q = 0 by King–Altman (no subtraction) and the closed-form v_sat;
viol = total excursion of the 18 rate constants outside [klo, khi] (log10 units)."""
function evalR!(ws::WS, x, klo, khi)
    lkf, lkr = logconsts(x)
    viol = 0.0
    for i in 1:9
        viol += max(0.0, lkf[i] - khi) + max(0.0, klo - lkf[i]) +
                max(0.0, lkr[i] - khi) + max(0.0, klo - lkr[i])
        ws.kf[i] = exp10(lkf[i])
        ws.kr[i] = exp10(lkr[i])
    end
    A = exp10(x[17])
    B = exp10(x[18])
    kf, kr, fr, br = ws.kf, ws.kr, ws.fr, ws.br
    fr[1] = kf[1] * A; br[1] = kr[1]
    fr[2] = kf[2] * B; br[2] = kr[2]
    fr[3] = kf[3] * B; br[3] = kr[3]
    fr[4] = kf[4] * A; br[4] = kr[4]
    fr[5] = kf[5];     br[5] = kr[5]
    fr[6] = kf[6];     br[6] = 0.0          # P = 0
    fr[7] = kf[7];     br[7] = 0.0          # Q = 0
    fr[8] = kf[8];     br[8] = 0.0          # Q = 0
    fr[9] = kf[9];     br[9] = 0.0          # P = 0
    VL.ka_weights!(ws.w, KAT, fr, br)
    v = VL.ka_release(mT, ws.w, fr, 3)
    R = v / VL.vsat_T(kf, kr)
    R, viol
end

mutable struct Tracker
    bestR::Float64
    bestx::Vector{Float64}
    nevals::Int
    nflag::Int
    flagged::Vector{Vector{Float64}}
    top::Vector{Tuple{Float64,Vector{Float64}}}
end
Tracker() = Tracker(-Inf, Float64[], 0, 0, Vector{Float64}[], Tuple{Float64,Vector{Float64}}[])

function record!(tr::Tracker, R, x)
    tr.nevals += 1
    isfinite(R) || return
    if R > 1 + 1e-12
        tr.nflag += 1
        length(tr.flagged) < 200 && push!(tr.flagged, collect(x))
    end
    if R > tr.bestR
        tr.bestR = R
        tr.bestx = collect(x)
    end
    if length(tr.top) < 40 || R > tr.top[end][1]
        push!(tr.top, (R, collect(x)))
        sort!(tr.top; by = first, rev = true)
        length(tr.top) > 40 && pop!(tr.top)
    end
end

"""Objective to minimize: -log10 R + 10 * violation; box bounds on x are enforced by the caller
(DE) or penalized (Nelder–Mead via lo/hi)."""
function objective(ws, tr, klo, khi, lo, hi)
    x -> begin
        pen = 0.0
        for j in eachindex(x)
            pen += max(0.0, x[j] - hi[j]) + max(0.0, lo[j] - x[j])
        end
        pen > 0 && return 1e3 + 10pen
        R, viol = evalR!(ws, x, klo, khi)
        viol == 0 && record!(tr, R, x)
        (R > 0 && isfinite(R) ? -log10(R) : 300.0) + 10viol
    end
end

"""Exact R at x in rationals (reverse constants of the chords recomputed exactly)."""
function exactR(x)
    setprecision(BigFloat, 256) do
        q(l) = toQ(exp10(big(l)))
        kf = [q(x[i]) for i in 1:9]
        kr = Vector{QQ}(undef, 9)
        for (j, i) in enumerate(TREE)
            kr[i] = q(x[9+j])
        end
        K(i) = kf[i] / kr[i]
        Keq = q(x[16])
        K4 = K(1) * K(3) / K(2)
        K8 = Keq / (K(1) * K(3) * K(5) * K(6))
        K9 = Keq / (K(1) * K(3) * K(5) * K(7))
        kr[4] = kf[4] / K4
        kr[8] = kf[8] / K8
        kr[9] = kf[9] / K9
        c = (q(x[17]), q(x[18]), zero(QQ), zero(QQ))
        v = rate_exact(mT, kf, kr, c)[1]
        vs = VL.vsat_T(kf, kr)
        vlim = rate_exact(mT, kf, kr, (qpow10(200), qpow10(200), zero(QQ), zero(QQ)))[1]
        v / vs, abs(vlim - vs) / vs
    end
end

# The search agent's targeted guesses (log10 kf of all steps, log10 kr of steps 1 2 3 5 6 7,
# Keq = 1), copied from its guesses.jl as adversarial seeds.
const GUESSES = [
    ([6, -2, -4, 6, 0, 2, 2, 2, 2], [2, 2, 4, -2, 0, 0]),
    ([6, -2, -4, 6, 4, 2, 2, 2, 2], [2, 2, 4, 2, 0, 0]),
    ([6, -2, -4, 6, 0, -2, -2, -2, -2], [2, 2, 4, -2, -4, -4]),
    ([6, 0, -6, 6, 2, 4, 4, 4, 4], [-4, 0, 0, 2, 6, 6]),
    ([4, 0, 0, 4, -2, 2, 2, 2, 2], [0, 0, 6, -4, 0, 0]),
    ([3, 3, 3, 3, 4, -2, -2, 0, 0], [0, 0, 0, 6, 0, 0]),
    ([0, 0, 0, 0, 0, 0, 0, 0, -3], [0, 0, 0, 0, 0, 0]),
]

function main()
    t00 = time()
    ws = WS()
    x = vcat(zeros(16), [1.0, 1.0])
    evalR!(ws, x, -9, 9)
    tt = @elapsed for _ in 1:200_000
        evalR!(ws, x, -9, 9)
    end
    @printf("King–Altman: %d spanning trees; %.2f µs per evaluation\n", VL.ntrees(KAT),
            tt / 200_000 * 1e6)
    # cross-check King–Altman against exact rationals on random points
    let rng = MersenneTwister(1), worst = 0.0, n = 0
        while n < 200
            y = vcat(18 .* rand(rng, 9) .- 9, zeros(6), 12 * rand(rng) - 6,
                     24 .* rand(rng, 2) .- 12)
            for j in 1:6
                y[9+j] = y[TREE[j]] - (12 * rand(rng) - 6)
            end
            R, viol = evalR!(ws, y, -9, 9)
            (viol > 0 || any(i -> abs(y[i]) > 9, 1:15)) && continue
            n += 1
            Rx, _ = exactR(y)
            worst = max(worst, Float64(abs(R - Rx) / Rx))
        end
        @printf("King–Altman Float64 R vs exact rational R, 200 random points: max rel err %.2e\n",
                worst)
    end
    tr = Tracker()
    results = Dict{String,Any}()
    phases = [
        # label, klo, khi, Keq range, log10 A,B range, share of the budget
        ("A1 k 1e-9..1e9, A,B 1e-12..1e12", -9.0, 9.0, (-6.0, 6.0), (-12.0, 12.0), 0.10),
        ("A2 k 1e-9..1e9, A,B 1e-12..1e12 (seed 2)", -9.0, 9.0, (-6.0, 6.0), (-12.0, 12.0), 0.08),
        ("A3 k 1e-9..1e9, A,B <= 1", -9.0, 9.0, (-6.0, 6.0), (-12.0, 0.0), 0.08),
        ("A4 k 1e-3..1e3, A,B 1e-12..1e12", -3.0, 3.0, (-6.0, 6.0), (-12.0, 12.0), 0.08),
        ("A5 k 1e-3..1e3, A,B 1e-6..1", -3.0, 3.0, (-3.0, 3.0), (-6.0, 0.0), 0.08),
        ("A6 k 1e-9..1e9, Keq = 1e6", -9.0, 9.0, (6.0, 6.0), (-12.0, 12.0), 0.05),
        ("A7 k 1e-9..1e9, Keq = 1e-6", -9.0, 9.0, (-6.0, -6.0), (-12.0, 12.0), 0.05),
        ("A8 k 1e-12..1e12, A,B 1e-12..1e12", -12.0, 12.0, (-6.0, 6.0), (-12.0, 12.0), 0.08),
    ]
    for (ip, (label, klo, khi, kq, cab, share)) in enumerate(phases)
        rng = MersenneTwister(1000 + ip)
        lo = vcat(fill(klo, 15), kq[1], cab[1], cab[1])
        hi = vcat(fill(khi, 15), kq[2], cab[2], cab[2])
        tr0 = Tracker()
        f = objective(ws, tr0, klo, khi, lo, hi)
        seeds = Vector{Vector{Float64}}()
        if ip == 1
            for (lkf, lkr) in GUESSES
                push!(seeds, vcat(Float64.(lkf), Float64.(lkr), 0.0, 0.0, 0.0))
            end
        end
        gens, restarts = de!(f, lo, hi, rng; NP = 150, tmax = share * BUDGET, seeds = seeds)
        # Nelder–Mead polish from the phase's 5 best points
        for (R, xb) in tr0.top[1:min(5, length(tr0.top))]
            nelder_mead(f, xb, 0.5; maxevals = 20_000)
        end
        @printf("%-45s evals %.2e  gens %d restarts %d  best R = %.17f  R>1+1e-12: %d  at log10(A,B) = %s\n",
                label, tr0.nevals, gens, restarts, tr0.bestR, tr0.nflag,
                isempty(tr0.bestx) ? "none" : string(round.(tr0.bestx[17:18]; digits = 3)))
        flush(stdout)
        results[label] = (bestR = tr0.bestR, bestx = tr0.bestx, nflag = tr0.nflag,
                          flagged = tr0.flagged, top = tr0.top, nevals = tr0.nevals)
        for (R, xb) in tr0.top
            record!(tr, R, xb)
        end
        append!(tr.flagged, tr0.flagged)
        tr.nflag += tr0.nflag
        tr.nevals += tr0.nevals
        serialize(OUT, results)
    end

    # Bilevel: sampled rate constants, inner sup over (A, B) on a grid plus pattern search
    let rng = MersenneTwister(77), tB = 0.30 * BUDGET, t0 = time(), nθ = 0, nrej = 0,
        bestR = -Inf, bestx = Float64[], ninterior = 0, nflag = 0, klo = -9.0, khi = 9.0,
        grid = collect(-12.0:1.0:12.0)
        y = zeros(18)
        while time() - t0 < tB
            mode = rand(rng)
            if mode < 0.5
                for i in 1:9
                    y[i] = klo + (khi - klo) * rand(rng)
                end
                for j in 1:6
                    y[9+j] = y[TREE[j]] - (12 * rand(rng) - 6)
                end
                y[16] = 12 * rand(rng) - 6
            elseif mode < 0.75
                for i in 1:15
                    y[i] = rand(rng, Bool) ? klo + 3 * rand(rng) : khi - 3 * rand(rng)
                end
                y[16] = 12 * rand(rng) - 6
            else
                src = tr.top[rand(rng, 1:length(tr.top))][2]
                σ = (0.1, 1.0, 3.0)[rand(rng, 1:3)]
                for i in 1:16
                    y[i] = src[i] + σ * randn(rng)
                end
            end
            y[17] = y[18] = 0.0
            _, viol = evalR!(ws, y, klo, khi)
            if viol > 0 || any(i -> y[i] < klo || y[i] > khi, 1:15)
                nrej += 1
                continue
            end
            nθ += 1
            bR, ba, bb = -Inf, 0.0, 0.0
            for a in grid, b in grid
                y[17] = a
                y[18] = b
                R, _ = evalR!(ws, y, klo, khi)
                R > bR && ((bR, ba, bb) = (R, a, b))
            end
            for h in (0.5, 0.25, 0.1, 0.03, 0.01, 0.003, 0.001)
                improved = true
                while improved
                    improved = false
                    for (da, db) in ((h, 0.0), (-h, 0.0), (0.0, h), (0.0, -h), (h, h), (-h, -h))
                        a2 = clamp(ba + da, -12.0, 12.0)
                        b2 = clamp(bb + db, -12.0, 12.0)
                        y[17] = a2
                        y[18] = b2
                        R, _ = evalR!(ws, y, klo, khi)
                        if R > bR
                            bR, ba, bb = R, a2, b2
                            improved = true
                        end
                    end
                end
            end
            (ba < 12 || bb < 12) && (ninterior += 1)
            bR > 1 + 1e-12 && (nflag += 1)
            y[17] = ba
            y[18] = bb
            if bR > bestR
                bestR = bR
                bestx = copy(y)
            end
            record!(tr, bR, y)
        end
        @printf("%-45s θ sets %d (rejected %d)  best inner-sup R = %.17f  R>1+1e-12: %d  sup not at corner (1e12,1e12): %d\n",
                "B  bilevel sampling, k 1e-9..1e9", nθ, nrej, bestR, nflag, ninterior)
        results["B bilevel"] = (bestR = bestR, bestx = bestx, nθ = nθ, nflag = nflag,
                                ninterior = ninterior)
        flush(stdout)
    end
    serialize(OUT, results)

    println("\n-- exact rational recheck of the overall best points --")
    worstgap = Inf
    for (R, xb) in tr.top[1:min(25, length(tr.top))]
        Rx, limerr = exactR(xb)
        gap = BigFloat(1 - Rx)
        worstgap = min(worstgap, Float64(gap))
        @printf("   Float64 R = %.17f   exact 1 - R = %.4e   (closed-form v_sat vs exact 1e200 limit rel diff %.1e)  log10(A,B) = %s\n",
                R, Float64(gap), Float64(limerr), string(round.(xb[17:18]; digits = 3)))
    end
    for xb in tr.flagged[1:min(50, length(tr.flagged))]
        Rx, _ = exactR(xb)
        @printf("   FLAGGED point: exact R - 1 = %+.4e\n", Float64(BigFloat(Rx - 1)))
    end
    @printf("\nTOTAL: evaluations %.3e, best Float64 R = %.17f, points with R > 1 + 1e-12: %d, smallest exact gap 1 - R among rechecked = %.3e, wall %.0f s\n",
            tr.nevals, tr.bestR, tr.nflag, worstgap, time() - t00)
    serialize(OUT, results)
end

main()
