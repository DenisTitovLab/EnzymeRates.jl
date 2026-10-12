# ABOUTME: Targeted parameter guesses, a random-sampling census, a product-activation example and
# ABOUTME: the dead-end abortive-complex analysis for the random-order bi-bi peak question.

include(joinpath(@__DIR__, "search_lib.jl"))

setprecision(BigFloat, 256)
const m = MECHS["T"]
# step indices of the textbook mechanism
# 1 E+A->EA  2 E+B->EB  3 EA+B->EAB  4 EB+A->EAB  5 EAB->EPQ
# 6 EPQ->EQ+P  7 EPQ->EP+Q  8 EQ->E+Q  9 EP->E+P      (tree: 1 2 3 5 6 7; chords 4 8 9)

d(ps...) = Dict{Int, Float64}(ps...)

const GUESSES = [
    ("slow EA+B, fast EB+A, fast EAB dissociation",
     d(1 => 6, 2 => -2, 3 => -4, 4 => 6, 5 => 0, 6 => 2, 7 => 2, 8 => 2, 9 => 2),
     d(1 => 2, 2 => 2, 3 => 4, 5 => -2, 6 => 0, 7 => 0)),
    ("same, fast catalysis",
     d(1 => 6, 2 => -2, 3 => -4, 4 => 6, 5 => 4, 6 => 2, 7 => 2, 8 => 2, 9 => 2),
     d(1 => 2, 2 => 2, 3 => 4, 5 => 2, 6 => 0, 7 => 0)),
    ("same, slow release",
     d(1 => 6, 2 => -2, 3 => -4, 4 => 6, 5 => 0, 6 => -2, 7 => -2, 8 => -2, 9 => -2),
     d(1 => 2, 2 => 2, 3 => 4, 5 => -2, 6 => -4, 7 => -4)),
    ("tight A binding to E, A-first route slow, B-first route fast",
     d(1 => 6, 2 => 0, 3 => -6, 4 => 6, 5 => 2, 6 => 4, 7 => 4, 8 => 4, 9 => 4),
     d(1 => -4, 2 => 0, 3 => 0, 5 => 2, 6 => 6, 7 => 6)),
    ("EAB dissociation >> catalysis, asymmetric routes",
     d(1 => 4, 2 => 0, 3 => 0, 4 => 4, 5 => -2, 6 => 2, 7 => 2, 8 => 2, 9 => 2),
     d(1 => 0, 2 => 0, 3 => 6, 5 => -4, 6 => 0, 7 => 0)),
    ("near-equilibrium catalysis (k-5 >> release)",
     d(1 => 3, 2 => 3, 3 => 3, 4 => 3, 5 => 4, 6 => -2, 7 => -2, 8 => 0, 9 => 0),
     d(1 => 0, 2 => 0, 3 => 0, 5 => 6, 6 => 0, 7 => 0)),
    ("product-activation example: EP -> E + P slow (1e-3), everything else 1",
     d(1 => 0, 2 => 0, 3 => 0, 4 => 0, 5 => 0, 6 => 0, 7 => 0, 8 => 0, 9 => -3),
     d(1 => 0, 2 => 0, 3 => 0, 5 => 0, 6 => 0, 7 => 0)),
]

println("== targeted guesses, textbook mechanism, Keq = 1 ==")
for (label, lkf, lkr) in GUESSES
    θ = theta_from(m, 0.0, lkf, lkr)
    kf, kr, viol = rates(m, θ, 0.0)
    R0, lc0 = inner_sup(m, θ, 0.0, :zero)
    RQ, lcQ = inner_sup(m, θ, 0.0, :Q)
    RP, lcP = inner_sup(m, θ, 0.0, :P)
    @printf("%s\n   box violation %.2g | R(P=Q=0) = %.15f at log10(A,B) = %s\n   R(P=0, Q free) = %.6g at log10(A,B,Q) = %s | R(Q=0, P free) = %.6g at %s\n",
            label, viol, R0, string(round.(lc0; digits = 2)), RQ,
            string(round.(lcQ; digits = 2)), RP, string(round.(lcP; digits = 2)))
end

println("\n== product-activation example: v(A = B = 1e8, P = 0, Q) / v_sat ==")
let lkf = GUESSES[end][2], lkr = GUESSES[end][3]
    θ = theta_from(m, 0.0, lkf, lkr)
    kf, kr, _ = rates(m, big.(θ), 0.0; T = BigFloat)
    vs = maximum(vsat_dirs(m, kf, kr))
    @printf("   v_sat = %.6e\n", Float64(vs))
    for q in -4:1:8
        v = rate(m, kf, kr, (big(1e8), big(1e8), big(0.0), exp10(big(q))))
        @printf("   Q = 1e%-3d  v/v_sat = %.6f\n", q, Float64(v / vs))
    end
    for a in (0, 2, 4, 6, 8)
        v = rate(m, kf, kr, (exp10(big(a)), exp10(big(a)), big(0.0), big(10.0)))
        @printf("   A = B = 1e%d, Q = 10:  v/v_sat = %.6f\n", a, Float64(v / vs))
    end
end

println("\n== random census, textbook mechanism, P = Q = 0 ==")
let rng = MersenneTwister(2026), nacc = 0, ntry = 0, worst = -Inf, nA = 0, nB = 0,
    maxpeakA = 1.0, grid = collect(-8:0.5:8)
    tree = (1, 2, 3, 5, 6, 7)
    while nacc < 5000
        ntry += 1
        lK = rand(rng, (-3.0, 0.0, 3.0))
        lkf = d((i => 12rand(rng) - 6 for i in 1:9)...)
        lkr = d((i => 12rand(rng) - 6 for i in tree)...)
        θ = theta_from(m, lK, lkf, lkr)
        kf, kr, viol = rates(m, θ, lK)
        viol > 0 && continue
        nacc += 1
        vs = maximum(vsat_dirs(m, kf, kr))
        R = maximum(rate(m, kf, kr, (exp10(a), exp10(b), 0.0, 0.0)) / vs
                    for a in grid, b in grid)
        worst = max(worst, R)
        # 1D peaks: in A at fixed B, and in B at fixed A
        for fixed in (-2.0, 0.0, 2.0)
            vA = [rate(m, kf, kr, (exp10(a), exp10(fixed), 0.0, 0.0)) for a in grid]
            vB = [rate(m, kf, kr, (exp10(fixed), exp10(b), 0.0, 0.0)) for b in grid]
            if maximum(vA) > vA[end] * (1 + 1e-6)
                nA += 1
                maxpeakA = max(maxpeakA, maximum(vA) / vA[end])
            end
            maximum(vB) > vB[end] * (1 + 1e-6) && (nB += 1)
        end
    end
    @printf("   %d feasible of %d drawn; max over samples of sup v/v_sat = %.15f\n",
            nacc, ntry, worst)
    @printf("   1D slices with an interior maximum (v_peak > v(1e8) by 1e-6): in A at fixed B %d of %d, in B at fixed A %d of %d; largest A-slice peak/v(A=1e8) = %.4g\n",
            nA, 3nacc, nB, 3nacc, maxpeakA)
end

println("\n== (e) dead-end EQ + A <-> EAQ ==")
let me = MECHS["e"], rng = MersenneTwister(5)
    # θ_e = [kf1..kf9, kf10, g_EA..g_EQ, g_EAQ]; textbook θ_T = θ_e[[1:9; 11:16]]
    for trial in 1:3
        θe = vcat(4 .* (2 .* rand(rng, 10) .- 1), 2 .* (2 .* rand(rng, 7) .- 1))
        kf, kr, _ = rates(me, big.(θe), 0.0; T = BigFloat)
        vs = vsat_dirs(me, kf, kr)
        vals = [Float64(rate(me, kf, kr, (exp10(big(a)), exp10(big(a)), big(0.0), big(0.0))))
                for a in -2:2:12]
        @printf("   trial %d: v(A = B = 1e-2,1e0,...,1e12) = %s ; v at (1e48,1e48) = %.3e\n",
                trial, string(round.(vals; sigdigits = 4)), Float64(vs[1]))
    end
    # sup over everything of v_e(c) / v_sat of the same enzyme without the dead-end step
    mT = MECHS["T"]
    nθ = ntheta(me)
    best = Ref(-Inf); bestx = Ref(Float64[])
    function F(x)
        maximum(abs, x) > 60 && return 1e12
        kf, kr, viol = rates(me, view(x, 1:nθ), 0.0)
        kfT, krT, _ = rates(mT, x[[1:9; 11:16]], 0.0)
        lc = view(x, nθ+1:nθ+2)
        for l in lc
            viol += max(0.0, abs(l) - 8)
        end
        R = rate(me, kf, kr, (exp10(lc[1]), exp10(lc[2]), 0.0, 0.0)) /
            maximum(vsat_dirs(mT, kfT, krT))
        isfinite(R) || return 1e11
        viol == 0 && R > best[] && (best[] = R; bestx[] = collect(x))
        -log10(max(R, 1e-300)) + 100viol
    end
    used = 0
    for r in 1:6
        x0 = vcat(3 .* (2 .* rand(rng, nθ) .- 1), 16 .* rand(rng, 2) .- 8)
        used += cmaes(F, x0, 2.0, rng; λ = (4 + floor(Int, 3log(nθ + 2))) * 2^(r - 1),
                      maxevals = 400_000 - used)
    end
    @printf("   sup v_deadend(c) / v_sat(same enzyme without EAQ) found = %.15f (should be <= 1)\n",
            best[])
end

println("\n== moderate examples of substrate peaks at P = Q = 0 (all K = 1, Keq = 1) ==")
let
    allk(n) = Dict{Int, Float64}(i => 0.0 for i in 1:n)
    tree = Dict{Int, Float64}(i => 0.0 for i in (1, 2, 3, 5, 6, 7, 10))
    for (key, lkf, label) in (
            ("f", merge(allk(11), Dict(11 => -3.0)),
             "(f) all k = 1 except EAQ <-> EA + Q at 1e-3 both ways"),
            ("e", allk(10), "(e) all k = 1 including EQ + A <-> EAQ"))
        mm = MECHS[key]
        θ = theta_from(mm, 0.0, lkf, tree)
        kf, kr, viol = rates(mm, big.(θ), 0.0; T = BigFloat)
        vs = maximum(vsat_dirs(mm, kf, kr))
        row = [Float64(rate(mm, kf, kr, (exp10(big(a)), big(1e8), big(0.0), big(0.0))))
               for a in -4:1:8]
        @printf("%s (box violation %.1g)\n   v_sat = %.4e; v(A = 1e-4..1e8 step 10x, B = 1e8) = %s\n",
                label, Float64(viol), Float64(vs), string(round.(row; sigdigits = 3)))
        if vs > 1e-30
            R, lc = inner_sup(mm, θ, 0.0, :zero)
            @printf("   sup v / v_sat = %.6g at log10(A,B) = %s\n", R,
                    string(round.(lc; digits = 2)))
        end
    end
end
