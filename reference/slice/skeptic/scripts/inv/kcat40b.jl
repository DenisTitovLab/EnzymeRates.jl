# ABOUTME: specs/40: the slice's own law at saturating substrates against its kcat, and the
# ABOUTME: leading faces kcat_terms builds, to locate the difference.
include(joinpath(@__DIR__, "..", "spotlib.jl"))
r = only(x for x in skeptic_records("specs") if x.p.id == "specs/40")
m = r.c.mech
sol = S.solve(m)
rng = Xoshiro(7)
params = (; (x => (x === :E_total ? big(1) : Oracle.logu(rng)) for x in S.parameters(m))...)
em = S.compile(m)
println("substrates ", m.substrates, " products ", m.products)
for X in (big(10)^10, big(10)^20, big(10)^40)
    for u in (big(1), big(7))
        c = (; (x => (x == m.substrates[1] ? X * u : x in m.substrates ? X : big(0))
                for x in S.metabolites(m))...)
        @printf("X=1e%d ratio %d: law %.12e\n", round(Int, log10(X)), Int(u),
                S.rate_equation(em, c, params))
    end
end
println("kcat ", Float64(S.kcat(m, params)))
A = S.allo_law(m, sol, S.kcat_mask(m))
nc = length(S.metabolites(m))
for (nm, p) in (("H", A.H), ("fA", A.fA), ("fI", A.fI), ("QA", A.active.Q), ("DA", A.active.D),
                ("QI", S.inactive_terms(A)[1]), ("DI", S.inactive_terms(A)[2]))
    cs = S.split_conc(p, A.vars)
    degs = sort(unique([sum((e for (v, e) in mo if A.vars.names[v] in m.substrates); init = 0//1)
                        for mo in keys(cs) if all(t -> A.vars.names[t.first] in m.substrates, mo)]))
    println(rpad(nm, 3), " terms ", S.nterms(p), "; substrate-only degrees ", degs)
end
println("n = ", A.n)
