# ABOUTME: Evaluates EnzymeRates.rate_equation at the peak-search optima and the hand-built examples
# ABOUTME: (product activation, abortive-complex peak) and compares the ratio with the GTH solver.

using EnzymeRates, Printf, Serialization, Random
const ER = EnzymeRates
include(joinpath(@__DIR__, "search_lib.jl"))
include(joinpath(@__DIR__, "pkgdefs.jl"))
setprecision(BigFloat, 256)

const RES = deserialize(joinpath(@__DIR__, "results_full.jls"))

function pkg_rate(em, params, c)
    mets = ER.metabolites(em)
    ER.rate_equation(em, NamedTuple{Tuple(mets)}(
        Tuple(Float64(c[findfirst(==(x), LIGANDS)]) for x in mets)), params)
end

function compare(key, θ, lK, cpeak, label)
    m = MECHS[key]
    em = PKG[key]()
    kf, kr, viol = rates(m, big.(θ), lK; T = BigFloat)
    vals = pkg_values(em, m, kf, kr)
    fitted = ER.fitted_params(em)
    params = NamedTuple{(fitted..., :Keq, :E_total)}(
        (Float64.(getindex.(Ref(vals), fitted))..., exp10(lK), 1.0))
    csat = (1e20, 1e20, 0.0, 0.0)
    vp_pkg, vs_pkg = pkg_rate(em, params, cpeak), pkg_rate(em, params, csat)
    vp_me, vs_me = rate(m, kf, kr, big.(cpeak)), rate(m, kf, kr, big.(csat))
    vs_lim = maximum(vsat_dirs(m, kf, kr; dirs = ((80, 80), (80, 40), (40, 80))))
    @printf("%s\n   peak c = %s (box violation %.1g)\n   package: v(peak) = %.10e  v(A=B=1e20) = %.10e  ratio %.10g\n   solver:  v(peak) = %.10e  v(A=B=1e20) = %.10e  ratio %.10g ; limit A,B->inf ratio %.10g\n",
            label, string(round.(Float64.(cpeak); sigdigits = 4)), Float64(viol), vp_pkg, vs_pkg,
            vp_pkg / vs_pkg, Float64(vp_me), Float64(vs_me), Float64(vp_me / vs_me),
            Float64(vp_me / vs_lim))
end

for (key, setting) in (("T", :Q), ("c", :Q), ("d1", :Q), ("f", :zero))
    r = RES[(key, setting, 0.0)]
    nθ = ntheta(MECHS[key])
    compare(key, r.x[1:nθ], 0.0, concs(setting, r.x[nθ+1:end]),
            "$key $setting Keq=1 search optimum")
end

let
    allk(n) = Dict{Int, Float64}(i => 0.0 for i in 1:n)
    tree = Dict{Int, Float64}(i => 0.0 for i in (1, 2, 3, 5, 6, 7))
    θ = theta_from(MECHS["T"], 0.0, merge(allk(9), Dict(9 => -3.0)), tree)
    compare("T", θ, 0.0, (1e8, 1e8, 0.0, 1e4),
            "T example: all constants 1, E(P) <-> E + P at 1e-3 both ways")
    tree10 = merge(tree, Dict(10 => 0.0))
    θf = theta_from(MECHS["f"], 0.0, merge(allk(11), Dict(11 => -3.0)), tree10)
    compare("f", θf, 0.0, (exp10(0.35), 1e8, 0.0, 0.0),
            "f example: all constants 1, E(A,Q) <-> E(A) + Q at 1e-3 both ways")
end
