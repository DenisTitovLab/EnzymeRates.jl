# ABOUTME: Cross-checks the independent GTH steady-state solver against EnzymeRates.rate_equation
# ABOUTME: for the textbook random-order bi-bi mechanism and every variant the package accepts.

using EnzymeRates, Random, Printf
const ER = EnzymeRates
include(joinpath(@__DIR__, "solver.jl"))
using .PeakSolver
include(joinpath(@__DIR__, "mechs.jl"))

setprecision(BigFloat, 256)

include(joinpath(@__DIR__, "pkgdefs.jl"))

function check(key; npts = 40, seed = 7)
    m = MECHS[key]
    em = try
        PKG[key]()
    catch err
        println("[$key] package rejects the mechanism: ", sprint(showerror, err)[1:min(end, 300)])
        return
    end
    mets = ER.metabolites(em)
    fitted = ER.fitted_params(em)
    @printf("[%s] %s: package fitted params %d, mine (free energies + SS kf) %d\n",
            key, m.name, length(fitted), ntheta(m))
    rng = MersenneTwister(seed)
    worst = Dict(:zero => 0.0, :Q => 0.0, :P => 0.0, :PQ => 0.0)
    worst_full = 0.0
    for ipt in 1:npts
        w = ipt <= npts ÷ 2 ? 3.0 : 6.0   # half moderate, half extreme parameter sets
        θ = [w * (2rand(rng) - 1) for _ in 1:ntheta(m)]
        lK = 3 * (2rand(rng) - 1)
        kf, kr, _ = rates(m, θ, lK; T = BigFloat)
        vals = pkg_values(em, m, kf, kr)
        params = NamedTuple{(fitted..., :Keq, :E_total)}(
            (Float64.(getindex.(Ref(vals), fitted))..., exp10(lK), 1.0))
        for cls in (:zero, :Q, :P, :PQ)
            lc = [3 * (2rand(rng) - 1) for _ in 1:4]
            c = BigFloat[exp10(big(x)) for x in lc]
            cls in (:zero, :Q) && (c[3] = 0)
            cls in (:zero, :P) && (c[4] = 0)
            cnt = NamedTuple{Tuple(mets)}(Tuple(Float64(c[findfirst(==(x), LIGANDS)])
                                                for x in mets))
            v_pkg = ER.rate_equation(em, cnt, params)
            v_me = rate(m, kf, kr, Tuple(c))
            err = abs(v_pkg - Float64(v_me)) / abs(Float64(v_me))
            worst[cls] = max(worst[cls], err)
            if err > 1e-6
                Γ = cls === :PQ ? Float64(c[3] * c[4] / (c[1] * c[2])) / exp10(lK) : 0.0
                @printf("   MISMATCH pt %d %s: pkg %.10e  mine %.10e  rel %.2e  Gamma/Keq %.6f\n",
                        ipt, cls, v_pkg, Float64(v_me), err, Γ)
            end
            # Full mode: every rate constant passed explicitly
            try
                names = ER.parameters(em, ER.Full)
                pf = NamedTuple{Tuple(names)}(Tuple(
                    n === :E_total ? 1.0 : n === :Keq ? exp10(lK) : Float64(vals[n])
                    for n in names))
                v_full = ER.rate_equation(em, cnt, pf, ER.Full)
                worst_full = max(worst_full, abs(v_full - Float64(v_me)) / abs(Float64(v_me)))
            catch err2
                ipt == 1 && cls == :zero &&
                    println("   Full mode unavailable: ", sprint(showerror, err2)[1:min(end, 200)])
            end
        end
    end
    @printf("   worst rel err (Reduced): zero %.2e | P=0,Q>0 %.2e | P>0,Q=0 %.2e | P,Q>0 %.2e\n",
            worst[:zero], worst[:Q], worst[:P], worst[:PQ])
    @printf("   worst rel err (Full):    %.2e\n", worst_full)
end

# Solver self-checks independent of the package: net P production equals net Q production
# and the cut-step flux, and Float64 agrees with BigFloat.
function selfcheck(key; npts = 200, seed = 11)
    m = MECHS[key]
    rng = MersenneTwister(seed)
    worst_cons = big(0.0); worst_f64 = 0.0
    for _ in 1:npts
        θ = [6 * (2rand(rng) - 1) for _ in 1:ntheta(m)]
        lK = 3 * (2rand(rng) - 1)
        kf, kr, _ = rates(m, θ, lK; T = BigFloat)
        c = Tuple(BigFloat[exp10(big(8 * (2rand(rng) - 1))) for _ in 1:4])
        π = steady_state(m, kf, kr, c)
        vP = netprod(m, kf, kr, c, π, 3)
        vQ = netprod(m, kf, kr, c, π, 4)
        vA = -netprod(m, kf, kr, c, π, 1)
        vB = -netprod(m, kf, kr, c, π, 2)
        vcut = sum(PeakSolver.stepnet(m, kf, kr, c, π, i) for i in m.cut)
        sc = max(abs(vP), big(1e-300))
        worst_cons = max(worst_cons, maximum(abs.([vQ, vA, vB, vcut] .- vP)) / sc)
        c0 = (c[1], c[2], big(0.0), big(0.0))
        v0 = rate(m, kf, kr, c0)
        kf64, kr64, _ = rates(m, θ, lK)
        v064 = rate(m, kf64, kr64, Float64.(c0))
        worst_f64 = max(worst_f64, abs(v064 - Float64(v0)) / Float64(v0))
    end
    @printf("[%s] self-check: max |vQ,vA,vB,vcut - vP|/|vP| = %.2e (BigFloat, P,Q>0); Float64 vs BigFloat at P=Q=0: %.2e\n", key, Float64(worst_cons), worst_f64)
end

for key in ("T", "a", "b", "c", "d1", "d2", "e", "f", "g")
    selfcheck(key)
end
for key in ("T", "g", "b", "a", "c", "d1", "d2", "e", "f")
    t = @elapsed check(key)
    @printf("   (%.1f s)\n", t)
end
