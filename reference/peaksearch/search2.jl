# ABOUTME: Multi-seed rerun of the peak search with the least-cancellation turnover formula: many
# ABOUTME: independent IPOP-CMA-ES searches per configuration, best point re-checked in BigFloat.

include(joinpath(@__DIR__, "search_lib.jl"))
using Serialization

setprecision(BigFloat, 256)

# Float64 accuracy of the turnover with both products present, against BigFloat
let rng = MersenneTwister(3), worst_new = 0.0, worst_cut = 0.0
    for key in ("T", "a", "b", "c", "d1", "d2", "e", "f", "g"), _ in 1:300
        m = MECHS[key]
        θ = [12rand(rng) - 6 for _ in 1:ntheta(m)]
        lK = 6rand(rng) - 3
        kfb, krb, _ = rates(m, big.(θ), lK; T = BigFloat)
        kf, kr, _ = rates(m, θ, lK)
        c = Tuple(exp10.(16 .* rand(rng, 4) .- 8))
        πb = steady_state(m, kfb, krb, big.(c))
        vb = turnover(m, kfb, krb, big.(c), πb)
        π = steady_state(m, kf, kr, c)
        vnew = turnover(m, kf, kr, c, π)
        vcut = sum(PeakSolver.stepnet(m, kf, kr, c, π, i) for i in m.cut)
        worst_new = max(worst_new, Float64(abs(vnew - vb) / abs(vb)))
        key == "b" || (worst_cut = max(worst_cut, Float64(abs(vcut - vb) / abs(vb))))
    end
    @printf("Float64 turnover, P,Q > 0, 2700 random points: max rel err %.2e (least-cancellation), %.2e (cut step only, RE variant excluded)\n",
            worst_new, worst_cut)
end

const RESULTS = Dict{Any, Any}()

function multiseed(key, setting, lK, nseeds; clo = -8.0, chi = 8.0, budget = 600_000)
    m = MECHS[key]
    nθ = ntheta(m)
    best = (R = -Inf, x = Float64[])
    nabove = 0
    for s in 1:nseeds
        r = search(m, setting, lK; seed = 1000s + 7, budget = budget, clo = clo, chi = chi)
        isempty(r.x) && continue
        r.R > 1 + 1e-9 && (nabove += 1)
        r.R > best.R && (best = r)
    end
    rb = evalR(m, big.(best.x), setting, lK; T = BigFloat, clo = clo, chi = chi)
    lc = best.x[nθ+1:end]
    @printf("%-3s %-5s log10Keq=%5.1f box[%g,%g] seeds=%d  best R=%.12g  Rbig-1=%+.3e  seeds with R>1+1e-9: %d  at log10 conc %s\n",
            key, setting, lK, clo, chi, nseeds, best.R, Float64(rb.R - 1), nabove,
            string(round.(lc; digits = 3)))
    flush(stdout)
    RESULTS[(key, setting, lK, chi)] = (; best..., Rbig = rb.R)
end

t0 = time()
for key in ("g", "a", "b", "d2"), setting in (:zero, :Q, :P, :PQ), lK in (-3.0, 0.0, 3.0)
    multiseed(key, setting, lK, 20)
end
for key in ("T", "c", "d1"), lK in (-3.0, 0.0, 3.0)
    multiseed(key, :zero, lK, 20)
end
for key in ("T", "c", "d1", "d2", "g"), chi in (2.0, 4.0)
    multiseed(key, :zero, 0.0, 10; chi = chi)
end
for key in ("T", "c", "d1"), setting in (:Q, :P, :PQ)
    multiseed(key, setting, 0.0, 3)
end
@printf("total %.0f s\n", time() - t0)
serialize(joinpath(@__DIR__, "results_multiseed.jls"), RESULTS)
