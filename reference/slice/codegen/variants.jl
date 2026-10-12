# ABOUTME: Experiment: codegen variants on the largest plain laws. A: factored, dependents
# ABOUTME: substituted (src); B: factored, dependents as locals; C: expanded D; each timed.
include(joinpath(@__DIR__, "variants_lib.jl"))

@noinline barrier(f::F, c, p) where {F} = f(c, p)
function t1n(f::F, c, p) where {F}
    barrier(f, c, p)
    best, acc = Inf, 0.0
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += barrier(f, c, p)
        end
        best = min(best, t / 10_000)
    end
    best
end

recs = [r for r in load_records(("init_bibi", "expand_bibi", "init_bibi_pp", "specs"))
        if !S.is_allosteric(r.c.mech)]
big = sort(recs; by = r -> -S.nterms(S.law(r.c.mech).D))[1:parse(Int, ARGS[1])]
rng = Xoshiro(7)
for r in big
    m = r.c.mech
    sol = S.solve(m)
    L = S.law(m, sol)
    c, p = rand_concs(m, rng), rand_params(m, rng)
    vb = Float64(law_big(L, c, p))
    out = String[]
    law7 = P7.factored_law(r.p; blocks = false)
    b7 = Gen.body(Gen.decode(Tuple(code(law7, (true, true, 0)))))
    s7 = r.p.samples[1]
    c7 = Tuple(Float64(Dict(s7.concs)[x]) for x in r.p.metabolites)
    p7 = Tuple(Float64(Dict(s7.params)[x]) for x in law7.vars.names[1:law7.np])
    f7 = @eval function $(gensym(:p7))(c, p)
        $b7
    end
    Base.invokelatest(f7, c7, p7)
    t7 = Base.invokelatest(t1n, f7, c7, p7)
    out7 = @sprintf("P7 segh %5.1f ns ops %4d", 1e9 * t7, countops(b7))
    for (lab, P) in (("A", program_a(L)), ("B", program_b(m, sol, L)), ("S", S.program(m)),
                     ("C", program_c(L)))
        b = S.body(P)
        f = @eval function $(gensym(:v))(concs, params)
            $b
        end
        s0 = snapshot()
        v = Base.invokelatest(f, c, p)
        s1 = snapshot()
        t = Base.invokelatest(t1n, f, c, p)
        abs(v - vb) / abs(vb) < 1e-12 || println("  $lab value off: $v vs $vb")
        push!(out, @sprintf("%s %5.1f ns ops %4d comp %3.0f ms jit %5d", lab, 1e9 * t,
                            countops(b), (s1.compile - s0.compile) / 1e6, s1.jit - s0.jit))
    end
    @printf("%-18s D %3d | %s | %s\n", r.p.id, S.nterms(L.D), join(out, " | "), out7)
    flush(stdout)
end
