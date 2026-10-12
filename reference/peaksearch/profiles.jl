# ABOUTME: Post-analysis of the peak search: 1D profiles of v / v_sat through every best point with
# ABOUTME: R > 1, and the analytic product-rescue bound for the textbook mechanism's best points.

include(joinpath(@__DIR__, "search_lib.jl"))
using Serialization

setprecision(BigFloat, 256)
const RES = deserialize(joinpath(@__DIR__, isempty(ARGS) ? "results_full.jls" : ARGS[1]))

"""Mean-cycle-time ratio T(Q=0)/T(Q->inf) at saturating A, B, P = 0 (textbook, steps 5..9)."""
function rescue_bound(kf, kr)
    k5, km5, k6, k7, k8, k9 = kf[5], kr[5], kf[6], kf[7], kf[8], kf[9]
    T0 = 1 / k5 + (km5 / k5 + 1 + k7 / k9 + k6 / k8) / (k6 + k7)
    Tinf = 1 / k5 + (km5 / k5 + 1 + k6 / k8) / k6
    T0 / Tinf
end

for (k, r) in sort(collect(RES); by = p -> string(p[1]))
    key, setting, lK = k
    setting === :zero_restricted && continue
    r.R > 1 + 1e-9 || continue
    m = MECHS[key]
    nθ = ntheta(m)
    θ = r.x[1:nθ]
    lc = r.x[nθ+1:end]
    kf, kr, _ = rates(m, big.(θ), lK; T = BigFloat)
    vs = maximum(vsat_dirs(m, kf, kr))
    @printf("\n%s %s log10Keq=%.0f  R=%.6g  best log10 conc = %s\n", key, setting, lK, r.R,
            string(round.(lc; digits = 3)))
    if key == "T" && setting === :Q
        @printf("   analytic rescue bound T(Q=0)/T(Q->inf) at saturation = %.6g\n",
                Float64(rescue_bound(kf, kr)))
    end
    names = setting === :zero ? ("A", "B") : setting === :Q ? ("A", "B", "Q") :
            setting === :P ? ("A", "B", "P") : ("A", "B", "P", "Q")
    for (j, nm) in enumerate(names)
        row = String[]
        for e in -8:2:8
            y = copy(lc)
            y[j] = e
            v = rate(m, kf, kr, concs(setting, big.(y), BigFloat))
            push!(row, @sprintf("%.3g", Float64(v / vs)))
        end
        @printf("   vary %s (1e-8,1e-6,...,1e8): %s\n", nm, join(row, "  "))
    end
end
