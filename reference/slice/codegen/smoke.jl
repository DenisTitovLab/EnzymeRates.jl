# ABOUTME: Smoke test of the generated rate_equation on the plain specs: value against the
# ABOUTME: BigFloat law, allocations and one timing; prints one generated body.
include(joinpath(@__DIR__, "lib.jl"))
rng = Xoshiro(1)
recs = [r for r in load_records(("specs",)) if !S.is_allosteric(r.c.mech)]
for (i, r) in enumerate(recs)
    m = r.c.mech
    em = S.compile(m)
    p, c = rand_params(m, rng), rand_concs(m, rng)
    L = S.law(m)
    t = @elapsed v = S.rate_equation(em, c, p)
    vb = law_big(L, c, p)
    a = allocs_of(S.rate_equation, em, c, p)
    t1 = t1_of(S.rate_equation, em, c, p)
    @printf("%-10s first %.3fs v %.6e rel %.1e allocs %d t1 %.1f ns\n", r.p.id, t, v,
            Float64(abs(v - vb) / abs(vb)), a, t1 * 1e9)
    i == 1 && println(S.body(S.program(m)))
end
