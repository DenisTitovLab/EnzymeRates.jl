include(joinpath(@__DIR__, "Encode.jl"))
ms = [Pop.load_population("specs"); Pop.load_population("expand_bibi")]
for id in ARGS
    m = only(filter(x -> x.id == id, ms))
    for (vname, mk, opt) in (("flat", P7.flat_law, (false, false, 0)), ("flatr", P7.flat_law, (true, false, 0)),
                             ("fact", P7.factored_law, (true, false, 0)), ("facth", P7.factored_law, (true, true, 0)))
        t0 = @elapsed law = mk(m)
        t1 = @elapsed L = encode(law, opt)
        t2 = @elapsed b = Gen.body(Gen.decode(typeof(L).parameters[1]))
        t3 = @elapsed M = Gen.define!(Symbol(id, vname), L, m.metabolites, [m.fitted_params; :Keq; :E_total])
        errs = map(m.samples) do s
            c = NamedTuple{Tuple(first.(s.concs))}(Tuple(last.(s.concs)))
            p = NamedTuple{Tuple(first.(s.params))}(Tuple(last.(s.params)))
            abs(Base.invokelatest(Gen.rate_equation, M, c, p) - s.v) / abs(s.v)
        end
        s = m.samples[1]
        c = NamedTuple{Tuple(first.(s.concs))}(Tuple(last.(s.concs)))
        p = NamedTuple{Tuple(first.(s.params))}(Tuple(last.(s.params)))
        f(M, c, p) = (Gen.rate_equation(M, c, p); @allocated Gen.rate_equation(M, c, p))
        a = Base.invokelatest(f, M, c, p)
        nops = Ref(0); walk(x) = (x isa Expr && x.head === :call && (nops[] += 1); x isa Expr && foreach(walk, x.args)); walk(b)
        println(rpad(id, 18), rpad(vname, 6), " law ", round(t0, digits=3), " enc ", round(t1, digits=3), " body ", round(t2, digits=3),
                " def ", round(t3, digits=3), " maxrelerr=", maximum(errs), " allocs=", a, " ops=", nops[], " ints=", length(code(law, opt)))
        flush(stdout)
    end
end
using InteractiveUtils
m = only(filter(x -> x.id == "specs/5", ms))
L = encode(P7.flat_law(m), (false, false, 0))
M = Gen.define!(:inl, L, m.metabolites, [m.fitted_params; :Keq; :E_total])
s = m.samples[1]
c = NamedTuple{Tuple(first.(s.concs))}(Tuple(last.(s.concs)))
p = NamedTuple{Tuple(first.(s.params))}(Tuple(last.(s.params)))
ct = Base.invokelatest(code_typed, Gen.rate_equation, Tuple{typeof(M), typeof(c), typeof(p)})
println(ct[1].first)
