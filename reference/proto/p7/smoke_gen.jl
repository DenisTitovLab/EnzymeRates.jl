# ABOUTME: Smoke test of the P7 codegen: kernels for a few mechanisms, values, allocations, inlining.
# ABOUTME: Throwaway.
include(joinpath(@__DIR__, "Encode.jl"))
using InteractiveUtils
ms = [Pop.load_population("specs"); Pop.load_population("expand_bibi")]
for id in ("specs/5", "specs/27", "specs/35", "specs/40", "expand_bibi/626")
    m = only(filter(x -> x.id == id, ms))
    for (vname, law, opt) in (("flat", P7.flat_law(m), (false, false, 0)), ("flatr", P7.flat_law(m), (true, false, 0)),
                               ("fact", P7.factored_law(m), (true, false, 0)), ("facth", P7.factored_law(m), (true, true, 0)))
        L = encode(law, opt)
        M = Gen.define!(Symbol(id, "_", vname), L, m.metabolites, [m.fitted_params; :Keq; :E_total])
        errs = map(m.samples) do s
            c = NamedTuple{Tuple(first.(s.concs))}(Tuple(last.(s.concs)))
            p = NamedTuple{Tuple(first.(s.params))}(Tuple(last.(s.params)))
            v = Base.invokelatest(Gen.rate_equation, M, c, p)
            abs(v - s.v) / abs(s.v)
        end
        s = m.samples[1]
        c = NamedTuple{Tuple(first.(s.concs))}(Tuple(last.(s.concs)))
        p = NamedTuple{Tuple(first.(s.params))}(Tuple(last.(s.params)))
        f(M, c, p) = (Gen.rate_equation(M, c, p); @allocated Gen.rate_equation(M, c, p))
        a = Base.invokelatest(f, M, c, p)
        b = Gen.body(typeof(L).parameters...)
        nops = Ref(0); walk(x) = (x isa Expr && x.head === :call && (nops[] += 1); x isa Expr && foreach(walk, x.args)); walk(b)
        println(rpad(id, 18), rpad(vname, 6), " maxrelerr=", maximum(errs), " allocs=", a, " ops=", nops[], " ints=", encoded_ints(law))
    end
end
# inlining check: wrapper calls the kernel through `invoke`
m = only(filter(x -> x.id == "specs/5", ms))
L = encode(P7.flat_law(m), (false, false, 0))
M = Gen.define!(:inl, L, m.metabolites, [m.fitted_params; :Keq; :E_total])
s = m.samples[1]
c = NamedTuple{Tuple(first.(s.concs))}(Tuple(last.(s.concs)))
p = NamedTuple{Tuple(first.(s.params))}(Tuple(last.(s.params)))
ct = Base.invokelatest(code_typed, Gen.rate_equation, Tuple{typeof(M), typeof(c), typeof(p)})
println(ct[1].first.code)
println(Gen.body(typeof(L).parameters...))
