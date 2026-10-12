# ABOUTME: Prints the fitted names and keys the adversarial run saved for the naming cases, and
# ABOUTME: the cube-root law's printed Haldane line.
include(joinpath(@__DIR__, "..", "advlib.jl"))
include(joinpath(@__DIR__, "..", "adv_cases.jl"))
rows = deserialize(joinpath(SCR, "slice", "skeptic", "data", "adv.jls"))
for x in rows
    x.outcome === :derived || continue
    x.name in ("substrate_inhibition_mod", "inhibitor_two_sites", "two_copies_inhibitor",
               "allo_n3_sites", "allo_equal_all_regulator") && println(rpad(x.name, 26), x.fitted)
end
for name in ("cube_root", "coef_name_collision")
    cs = only(c for c in CASES if c.name == name)
    m = S.Converted(cs.thunk(), 0).mech
    t = S.rate_equation_string(m)
    println("== ", name)
    foreach(l -> (occursin("^(", l) || occursin("Haldane", l) || startswith(l, "N =")) &&
                 println("  ", first(l, 160)), split(t, "\n"))
    name == "coef_name_collision" && foreach(l -> startswith(l, "d_A_B") && println("  ", l),
                                             split(S.coefficient_form(m), "\n"))
end
p = build("dsl", [:S], [:P], [st(E(), E(:S); c = [:S], g = 1), st(E(:S), E(); r = [:P], g = 2),
                               st(E(:S), E(:S, Mod(:S)); c = [Mod(:S)], re = true, g = 3)])
m = S.Converted(p, 0).mech
println("DSL text of the copy step: ", S.step_text(m, m.edges[end]))
