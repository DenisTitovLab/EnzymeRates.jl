# ABOUTME: Self-check of the solver on the rapid-equilibrium-release variant: net consumption of A
# ABOUTME: and B and the catalytic-step flux agree (RE-step fluxes cancel and are left out).

include(joinpath(@__DIR__, "search_lib.jl"))
setprecision(BigFloat, 256)
let m = MECHS["b"], rng = MersenneTwister(11), worst = big(0.0)
    for _ in 1:200
        θ = [6 * (2rand(rng) - 1) for _ in 1:ntheta(m)]
        lK = 3 * (2rand(rng) - 1)
        kf, kr, _ = rates(m, θ, lK; T = BigFloat)
        c = Tuple(BigFloat[exp10(big(8 * (2rand(rng) - 1))) for _ in 1:4])
        π = steady_state(m, kf, kr, c)
        vA = -netprod(m, kf, kr, c, π, 1)
        vB = -netprod(m, kf, kr, c, π, 2)
        vcut = sum(PeakSolver.stepnet(m, kf, kr, c, π, i) for i in m.cut)
        worst = max(worst, max(abs(vA - vcut), abs(vB - vcut)) / abs(vcut))
    end
    @printf("[b] self-check: max |vA - vcut|, |vB - vcut| over |vcut| = %.2e (BigFloat, P,Q > 0)\n",
            Float64(worst))
end
