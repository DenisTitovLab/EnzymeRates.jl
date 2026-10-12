# ABOUTME: Diagnoses kcat failures: prints the concentration supports of the core law's
# ABOUTME: numerator and denominator at zero products, and the top-degree face's ratios.
include(joinpath(@__DIR__, "..", "solve", "lib.jl"))
want = ARGS
recs = filter(r -> r.p.id in want, load_records(["specs", "handwritten", "expand_bibi"]))
rng = Xoshiro(3)
for r in recs
    m = r.c.mech
    sol = S.solve(m)
    core = S.catalytic_core(m)
    println("=== ", r.p.id, " core edges ", count(core), " of ", length(core))
    L = S.law(m, sol, core)
    nc = length(S.metabolites(m))
    names = L.vars.names
    prods = [i for i in 1:nc if names[i] in m.products]
    keep(p) = S.Poly(filter(t -> all(v -> !(v.first in prods), t.first), p.terms))
    N = keep(L.Q * S.conc_mono(L.vars, m.substrates)); D = keep(L.D)
    dn, nn = S.split_conc(D, L.vars), S.split_conc(N, L.vars)
    deg(mo) = sum(last, mo; init = 0//1)
    td = maximum(deg, keys(dn)); tn = maximum(deg, keys(nn))
    println("top degree D ", td, " N ", tn)
    pars = random_params(m, rng)
    vals = [L.vars.conc[i] ? big(1) : BigFloat(getfield(pars, n)) for (i, n) in enumerate(names)]
    for (lab, g) in (("D", dn), ("N", nn))
        for mo in sort(collect(keys(g)); by = deg)
            println("  ", lab, " ", S.term_text(1, mo, names), "  terms ", S.nterms(g[mo]),
                    deg(mo) == td ? "  value " * string(Float64(S.evaluate(g[mo], vals))) : "")
        end
    end
    for mo in keys(dn)
        deg(mo) == td || continue
        haskey(nn, mo) && println("  ratio at ", S.term_text(1, mo, names), " = ",
                                  Float64(S.evaluate(nn[mo], vals) / S.evaluate(dn[mo], vals)))
    end
end
