# ABOUTME: Evaluates the generated code and the law of two records with each metabolite at
# ABOUTME: zero, to explain NaN relative errors in the zero check; usage: peek_zero.jl
include(joinpath(@__DIR__, "lib.jl"))
for (set, id) in (("specs", "specs/24"), ("p2_fuzz_dump_flip", "fuzz/262"))
    r = only(filter(r -> r.p.id == id, records_of(set)))
    m = r.c.mech
    L = S.rate_law(m)
    em = S.compile(m)
    rng = Xoshiro(3)
    c = (; (x => 10.0^(4rand(rng) - 2) for x in S.metabolites(m))...)
    p = (; (x => (x === :E_total ? 1.0 : 10.0^(4rand(rng) - 2)) for x in S.parameters(m))...)
    println(id, " ", r.p.source, " metabolites ", S.metabolites(m))
    for x in S.metabolites(m)
        c0 = merge(c, NamedTuple{(x,)}((0.0,)))
        v = S.rate_equation(em, c0, p)
        exact = law_value(L, Dict(pairs(c0)), Dict(pairs(p)))
        println("  $x = 0: code $v, law $(Float64(exact))")
    end
end
