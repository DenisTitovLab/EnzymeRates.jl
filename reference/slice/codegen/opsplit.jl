# ABOUTME: Experiment: operation counts per statement of variant B against P7's segh body on
# ABOUTME: one law, to locate where the generated bodies differ in size.
include(joinpath(@__DIR__, "variants_lib.jl"))
set, id = ARGS
r = only(filter(x -> x.p.id == id, load_records((set,))))
m = r.c.mech
sol = S.solve(m)
L = S.law(m, sol)
law7 = P7.factored_law(r.p; blocks = false)
b7 = Gen.body(Gen.decode(Tuple(code(law7, (true, true, 0)))))
bB = S.body(program_b(m, sol, L))
for (lab, b) in (("P7", b7), ("B", bB))
    println("== $lab: total ops $(countops(b))")
    for ex in b.args
        ex isa Expr || continue
        n = countops(ex)
        n > 3 && println("  ", n, "  ", first(string(ex), 60))
    end
end
println("P7 deps: ", length(law7.deps), "; mine: ", length(sol.deps), " (folds ",
        length(sol.fold), ")")
println("P7 Q terms: ", [length(nd) for nd in law7.nodes])
for ex in bB.args
    ex isa Expr && ex.head === :(=) && countops(ex) <= 12 && println("  ", ex)
end
for d in sol.deps
    println("  dep col $(d.col) $(S.column_names(m, sol)[d.col]) label $(d.label) keq $(d.keq) terms ",
            [(S.column_names(m, sol)[f], e) for (f, e) in d.terms], haskey(sol.fold, d.col) ? " fold" : "")
end
println("P7 dependents: ", r.p.dependent_monomials)
