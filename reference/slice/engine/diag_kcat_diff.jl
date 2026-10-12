# ABOUTME: Explains kcat differences from the package: per product-free pattern, the ratio of
# ABOUTME: numerator to denominator coefficients in the full law and in the core law.
include(joinpath(@__DIR__, "lib.jl"))
ids = ARGS
recs = filter(r -> r.p.id in ids && !S.is_allosteric(r.c.mech),
              load_records(["expand_bibi", "p1_adv_pkg", "specs", "handwritten"]))
for r in recs
    m, p = r.c.mech, r.p
    sol = S.solve(m)
    pars = my_params(r, sol, pkg_values(p, p.samples[1]))
    println("=== ", r.set, " ", p.id, " ", p.source, "  package kcat ", p.samples[1].kcat,
            "  mine ", Float64(S.kcat(m, pars)))
    fwd = S.orient(m)
    core = S.catalytic_core(m)
    for k in eachindex(m.edges)
        println(core[k] ? "   " : " X ", S.step_text(m, S.oriented(m, k, fwd)))
    end
    for (lab, mask) in (("full", trues(length(m.edges))), ("core", core))
        L = S.law(m, sol, mask)
        vals = [L.vars.conc[i] ? big(1) : BigFloat(getfield(pars, n))
                for (i, n) in enumerate(L.vars.names)]
        prods = Set(L.vars.index[x] for x in m.products)
        n = S.split_conc(S.Poly(filter(t -> t.second > 0, (L.H * L.Q).terms)), L.vars)
        d = S.split_conc(L.D, L.vars)
        ks = sort!([k for k in keys(d) if !any(v -> v.first in prods, k)];
                   by = k -> sum(last, k; init = 0 // 1))
        println("  ", lab, " law: product-free denominator patterns (n/d where n exists)")
        for k in ks
            println("    ", rpad(S.term_text(1, k, L.vars.names), 16), " d ",
                    Float64(S.evaluate(d[k], vals)), haskey(n, k) ?
                    "  n/d " * string(Float64(S.evaluate(n[k], vals) / S.evaluate(d[k], vals))) : "")
        end
    end
end
