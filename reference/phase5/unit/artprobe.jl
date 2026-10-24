# ABOUTME: Prints where the printing-artifact patterns match in the slice-layout texts of given
# ABOUTME: mechanisms; args: <pop> <id>... (ids from summarize.jl's artifact list).
include(joinpath(@__DIR__, "ulib.jl"))
include(joinpath(@__DIR__, "pops.jl"))

pop, ids = ARGS[1], Set(ARGS[2:end])
for (id, text, subs, prods) in population(pop)
    id in ids || continue
    m = parse_neutral(text, subs, prods)
    sol = S.solve(m)
    law = S.rate_law(m, sol)
    allo = law isa S.AlloLaw
    confs = allo ? filter(x -> x !== nothing, [law.active, law.inactive]) : [law]
    norms = [spec_norm(m, L) for L in confs]
    full2, coef2 = allo ? allo_texts(m, law, sol, norms; layout = :slice, refstyle = :product) :
                          plain_texts(m, law, sol, norms[1]; layout = :slice, refstyle = :product)
    println("=== ", id)
    for (label, t) in (("full", full2), ("coefficients", coef2)), r in ARTIFACTS,
        x in eachmatch(r, t)
        i = x.offset
        println("  ", label, " ", r, ": …", t[max(1, prevind(t, i, 40)):min(end, nextind(t, i, 40))], "…")
    end
end
