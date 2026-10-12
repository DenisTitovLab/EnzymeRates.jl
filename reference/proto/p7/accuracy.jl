# ABOUTME: Float64 accuracy of the generated evaluators near equilibrium (Γ = Keq·(1+δ)) and at random
# ABOUTME: points, against BigFloat evaluation of the same law; P7 probe.
include(joinpath(@__DIR__, "Encode.jl"))
include(joinpath(@__DIR__, "Select.jl"))
using Printf, Random, Statistics
setprecision(BigFloat, 256)
const VARS = (("flat", false, (false, false)), ("flatr", false, (true, false)), ("flath", false, (true, true)),
              ("fact", true, (true, false)), ("facth", true, (true, true)))
errs = Dict(v[1] => Dict(δ => Float64[] for δ in (0.0, 1e-3, 1e-6, 1e-9)) for v in VARS)
set = selected()
for (k, (pop, m)) in enumerate(set)
    rng = Xoshiro(k)
    fl, fx = P7.flat_law(m), P7.factored_law(m)
    Ms = Dict(v => Gen.define!(Symbol("acc_", v, "_", m.id), encode(fac ? fx : fl, (o..., k)), m.metabolites,
                               [m.fitted_params; :Keq; :E_total]) for (v, fac, o) in VARS)
    s = m.samples[1]
    p = NamedTuple{Tuple(first.(s.params))}(Tuple(Float64.(last.(s.params))))
    keq = Dict(s.params)[:Keq]
    for δ in (0.0, 1e-3, 1e-6, 1e-9), _ in 1:20
        c = Dict(x => 10.0^(4 * rand(rng) - 2) for x in m.metabolites)
        if δ > 0       # set the first product so that Γ = Keq·(1+δ)
            Γ = prod(c[x] for x in m.products) / prod(c[x] for x in m.substrates)
            c[m.products[1]] *= keq * (1 + δ) / Γ
        end
        cn = NamedTuple{Tuple(m.metabolites)}(Tuple(c[x] for x in m.metabolites))
        x = P7.sample_vals(fl, Pop.PSample(collect(pairs(cn)), s.params, s.dep_values, s.v, s.full_params, s.v_full, s.kcat))
        ref = P7.value(fl, x)
        for (v, _, _) in VARS
            got = Base.invokelatest(Gen.rate_equation, Ms[v], cn, p)
            push!(errs[v][δ], Float64(abs(got - ref) / abs(ref)))
        end
    end
end
println("relative error vs BigFloat over $(length(set)) mechanisms × 20 points (median / max)")
for (v, _, _) in VARS
    @printf("%-6s", v)
    for δ in (0.0, 1e-3, 1e-6, 1e-9)
        e = errs[v][δ]
        @printf("  δ=%-6s %.1e / %.1e", δ == 0 ? "rand" : string(δ), median(e), maximum(e))
    end
    println()
end
