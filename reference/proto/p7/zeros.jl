# ABOUTME: Zero concentrations: no IR term carries a negative concentration exponent (all populations),
# ABOUTME: and the generated kernels stay finite and exact at zero products / one zero metabolite.
include(joinpath(@__DIR__, "Encode.jl"))
include(joinpath(@__DIR__, "Select.jl"))
setprecision(BigFloat, 256)
negc = 0; nlaws = 0
for pop in Pop.populations(), m in Pop.load_population(pop)
    m.ok || continue
    for law in (P7.flat_law(m), P7.factored_law(m), P7.factored_law(m; expand_den = true))
        global nlaws += 1
        c0 = law.np + length(law.deps)
        any(e < 0 && v > c0 for nd in law.nodes for t in nd for (v, e) in t.mono) && (global negc += 1)
    end
end
println("laws with a negative concentration exponent: $negc of $nlaws")
bad = 0; checked = 0; maxerr = 0.0
for (k, (pop, m)) in enumerate(selected())
    s = m.samples[1]
    p = NamedTuple{Tuple(first.(s.params))}(Tuple(Float64.(last.(s.params))))
    fl = P7.flat_law(m)
    laws = (("flat", fl, (false, false, k)), ("facth", P7.factored_law(m), (true, true, k)),
            ("halh", P7.factored_law(m; expand_den = true), (true, true, k)))
    Ms = [(n, Gen.define!(Symbol("z_", n, "_", m.id), encode(l, o), m.metabolites, [m.fitted_params; :Keq; :E_total])) for (n, l, o) in laws]
    base = Dict(s.concs)
    cases = [[x => 0.0 for x in m.products]; [[x => 0.0] for x in m.metabolites]...]
    cases = [Dict(m.products .=> 0.0); [Dict(x => 0.0) for x in m.metabolites]]
    for z in cases
        c = merge(base, z)
        cn = NamedTuple{Tuple(m.metabolites)}(Tuple(c[x] for x in m.metabolites))
        ref = P7.value(fl, P7.sample_vals(fl, Pop.PSample(collect(pairs(cn)), s.params, s.dep_values, s.v, s.full_params, s.v_full, s.kcat)))
        for (n, M) in Ms
            global checked += 1
            v = Base.invokelatest(Gen.rate_equation, M, cn, p)
            ok = isfinite(v) && (iszero(ref) ? abs(v) < 1e-300 : abs(v - ref) <= 1e-12 * abs(ref))
            ok || (global bad += 1; println("  $(m.id) $n zeros=$(collect(keys(z))) v=$v ref=$(Float64(ref))"))
        end
    end
end
println("zero-concentration evaluations: $checked, wrong or non-finite: $bad")
