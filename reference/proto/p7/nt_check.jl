# ABOUTME: The NamedTuple wrapper with permuted and extra fields: same value, zero allocations.
# ABOUTME: Throwaway check for the P7 report.
include(joinpath(@__DIR__, "Encode.jl"))
m = only(filter(x -> x.id == "expand_bibi/626", Pop.load_population("expand_bibi")))
M = Gen.define!(:ntc, encode(P7.factored_law(m; blocks = false), (true, true, 0)), m.metabolites, [m.fitted_params; :Keq; :E_total])
s = m.samples[1]
c = NamedTuple{Tuple(first.(s.concs))}(Tuple(last.(s.concs)))
p = NamedTuple{Tuple(first.(s.params))}(Tuple(last.(s.params)))
cr = NamedTuple{Tuple(reverse(first.(s.concs)))}(Tuple(reverse(last.(s.concs))))
pr = merge((extra = 1.0,), NamedTuple{Tuple(reverse(first.(s.params)))}(Tuple(reverse(last.(s.params)))))
f(M, c, p) = (Gen.rate_equation(M, c, p); @allocated Gen.rate_equation(M, c, p))
v1 = Base.invokelatest(Gen.rate_equation, M, c, p); v2 = Base.invokelatest(Gen.rate_equation, M, cr, pr)
println("same value: ", v1 == v2, " (", v1, " vs package ", s.v, "); allocs ", Base.invokelatest(f, M, c, p), " / ", Base.invokelatest(f, M, cr, pr))
