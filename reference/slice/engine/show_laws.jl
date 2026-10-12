# ABOUTME: Prints the full law and the coefficient form of the report's example mechanisms,
# ABOUTME: with each law's normalizer rule, steps and kcat candidates.
include(joinpath(@__DIR__, "lib.jl"))
ids = isempty(ARGS) ? ["specs/1", "specs/5", "specs/20", "specs/7", "adv2/1", "specs/11",
                       "handwritten/21", "specs/10"] : ARGS
recs = load_records()
for id in ids
    r = first(filter(r -> r.p.id == id, recs))
    m = r.c.mech
    L = S.law(m)
    rule, _ = S.normalizer(m, L)
    println("=== ", r.set, " ", id, " ", r.p.source, " (normalizer: ", rule, ")")
    fwd = S.orient(m)
    for k in eachindex(m.edges)
        println("    ", S.step_text(m, S.oriented(m, k, fwd)))
    end
    println("--- full law")
    println(S.rate_equation_string(m, L))
    println("--- coefficient form")
    println(S.coefficient_form(m, L))
    println()
end
