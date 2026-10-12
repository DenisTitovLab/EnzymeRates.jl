# ABOUTME: Re-checks the laws with the largest oracle errors at a faster rapid-equilibrium scale
# ABOUTME: (Λ = 1e60 against 1e40): an error that shrinks with Λ is the oracle's, not the law's.
include(joinpath(@__DIR__, "lib.jl"))

rows = [x for l in readdir(IDATA) if startswith(l, "pipe_")
        for x in deserialize(joinpath(IDATA, l)) if x.ok]
worst = first(sort(rows; by = x -> -x.oracle), 12)
rng = Xoshiro(5)
for x in worst
    r = only(filter(r -> r.c.index == x.index, records_of(x.set)))
    m = r.c.mech
    sol = S.solve(m)
    L = S.rate_law(m, sol)
    e40, e60 = 0.0, 0.0
    for _ in 1:3
        v, cv, conc = random_point(r, sol, rng)
        law = law_value(L, conc, v)
        o40 = oracle_value(r, sol, v, cv, conc)
        o60 = S.is_allosteric(m) ? o40 :
              Oracle.steady_state_v(r.p, step_constants(r, sol, cv, false), conc;
                                    Λ = big(10)^60)
        e40, e60 = max(e40, relerr(law, o40)), max(e60, relerr(law, o60))
    end
    @printf("%-20s %-22s pipeline %.1e; now at Λ = 1e40 %.1e, at Λ = 1e60 %.1e\n", x.set,
            x.id, x.oracle, e40, e60)
end
