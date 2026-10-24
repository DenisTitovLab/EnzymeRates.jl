# ABOUTME: Hand-built unit-rule fixtures (ratio form, fractional λ, RE Theorell–Chance) through
# ABOUTME: every rule, printed in full and in coefficients, with value and unit checks.
include(joinpath(@__DIR__, "ulib.jl"))

met(x) = S.Met(x, false)
form(c, b...) = S.RawForm(c, [met(x) for x in b], Symbol[], Symbol[])
step(a, b, cons, rel, re) = S.RawStep(a, b, met.(cons), met.(rel), re)
E, F = form(:E), form(:F)

fixtures = [
    # Two half-reactions between two free forms, each A → P: every term holds A or P.
    "F1 two-cycle uni-uni (ratio form)" => S.Mechanism([:A], [:P],
        [step(E, F, [:A], [:P], false), step(F, E, [:A], [:P], false)], [1, 2]),
    # F1 plus a dead end on F that takes up two A at once: terms A, P, A³, A²P.
    "F2 two-cycle uni-uni + F(A,A) dead end (fractional λ)" => S.Mechanism([:A], [:P],
        [step(E, F, [:A], [:P], false), step(F, E, [:A], [:P], false),
         step(F, form(:F, :A, :A), [:A, :A], Symbol[], true)], [1, 2, 3]),
    # Rapid-equilibrium Theorell–Chance step between steady-state A and Q bindings.
    "F3 RE Theorell–Chance (ratio form)" => S.Mechanism([:A, :B], [:P, :Q],
        [step(E, form(:E, :A), [:A], Symbol[], false),
         step(form(:E, :A), form(:E, :Q), [:B], [:P], true),
         step(form(:E, :Q), E, Symbol[], [:Q], false)], [1, 2, 3]),
    # F2 with the dead end as two single bindings: terms A, P, A², A³, AP, A²P.
    "F4 two-cycle uni-uni + F(A) + F(A,A) dead ends" => S.Mechanism([:A], [:P],
        [step(E, F, [:A], [:P], false), step(F, E, [:A], [:P], false),
         step(F, form(:F, :A), [:A], Symbol[], true),
         step(form(:F, :A), form(:F, :A, :A), [:A], Symbol[], true)], [1, 2, 3, 4]),
]

for (name, m) in fixtures
    println("=== ", name)
    try
        sol = S.solve(m)
        law = S.rate_law(m, sol)
        D = S.split_conc(law.D, law.vars)
        println("D monomials: ", join(sort([isempty(mo) ? "1" : join(("$(law.vars.names[v])^$e" for (v, e) in mo), "*") for mo in keys(D)]), ", "))
        r = analyze(m, Xoshiro(3))
        for (t, x) in r.rules
            println("rules: spec=", x.spec_kind, " λ=[", x.spec_lambda, "] slice=", x.slice_kind,
                    " cleland=", x.cleland, " sameC=", x.same_cleland, " orders=", x.same_order)
        end
        @printf("errors: full %.1e coef %.1e full2 %.1e slice %.1e/%.1e ref %.1e; units full %s coef %s (%s) slice %s/%s (%s); artifacts %d\n",
                r.e_full, r.e_coef, r.e_full2, r.e_sfull, r.e_scoef, r.e_ref, r.units_full, r.units_coef,
                r.units_note, r.units_sfull, r.units_scoef, r.units_snote, r.art)
        println(r.full, "\n--- coefficient form\n", r.coef)
    catch err
        println("ERROR: ", sprint(showerror, err))
        Base.show_backtrace(stdout, catch_backtrace())
        println()
    end
end
