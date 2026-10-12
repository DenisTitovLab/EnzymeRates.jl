# ABOUTME: Calibration: one law timed through the generated rate_equation and through an
# ABOUTME: eval'd function with the same body, interleaved, to compare the two harnesses.
include(joinpath(@__DIR__, "lib.jl"))

@noinline barrier(em, c, p) = S.rate_equation(em, c, p)
@noinline fbarrier(f, c, p) = f(c, p)
function tloop(g, a, c, p)
    g(a, c, p)
    best = Inf
    for _ in 1:5
        acc = 0.0
        t = @elapsed for _ in 1:10_000
            acc += g(a, c, p)
        end
        best = min(best, t / 10_000)
        isfinite(acc) || error()
    end
    1e9 * best
end
for (set, id) in Iterators.partition(ARGS, 2)
    r = only(filter(x -> x.p.id == id, load_records((set,))))
    m = r.c.mech
    em = S.compile(m)
    rng = Xoshiro(3)
    c, p = rand_concs(m, rng), rand_params(m, rng)
    b = S.body(S.program(m))
    f = @eval function $(gensym(:v))(concs, params)
        $b
    end
    Base.invokelatest(f, c, p)
    for rep in 1:3
        t_gen = tloop(barrier, em, c, p)
        t_eval = Base.invokelatest(tloop, fbarrier, f, c, p)
        t_gen2 = Base.invokelatest(tloop, barrier, em, c, p)
        @printf("%-18s rep %d: generated %.1f ns, eval'd %.1f ns, generated via invokelatest %.1f ns\n",
                id, rep, t_gen, t_eval, t_gen2)
    end
end
