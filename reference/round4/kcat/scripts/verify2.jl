# ABOUTME: Lists the steady-state random-order laws (no SS cycle through product-holding forms)
# ABOUTME: whose rate exceeds kcat in the harder search of verify.jl, with their steps.
include(joinpath(@__DIR__, "klib.jl"))

conc_dict(m, c) = Dict{Symbol, BigFloat}(n => big(c[i]) for (i, n) in enumerate(S.metabolites(m)))

rows = [x for l in ("core", "expand_terter", "terter_sample")
        for x in deserialize(joinpath(KDATA, "rows_$l.jls")) if x.ok]
over = filter(x -> x.cls === :finite && x.ratio > 1.01, rows)
ssub = filter(x -> x.cls === :finite && !x.allo && x.feat.ss_nt_sub > 0 && x.feat.ss_nt_prod == 0, rows)
need = Set((x.set, x.index) for x in ssub)
recs = Dict{Tuple{String, Int}, Any}()
for set in unique(first.(need))
    if set in CORE_SETS || set == "expand_terter"
        foreach(r -> (r.set, r.c.index) in need && (recs[(r.set, r.c.index)] = r), records_of(set))
    else
        ps = deserialize(joinpath(IDATA, "pop_$set.jls"))
        for (s, i) in need
            s == set && (recs[(s, i)] = (set = s, p = ps[i], c = S.Converted(ps[i], i)))
        end
        ps = nothing
        GC.gc()
    end
end
println(length(over), " sup > 1.01 kcat laws, ", length(ssub), " plain laws whose SS non-turnover cycles avoid product-holding forms")

# 2. harder search in random-order SS mechanisms without SS cycles through product-holding forms
function wide_params(m, sol, rng)
    names = S.column_names(m, sol)
    x = Dict{Symbol, BigFloat}(n => big(10.0^(-3 + 6 * rand(rng))) for n in S.parameters(m))
    cv = Dict{Int, BigFloat}(i => x[names[i]] for i in sol.indep)
    for d in sol.deps
        cv[d.col] = prod((S.pow(cv[f], e) for (f, e) in d.terms); init = S.pow(x[:Keq], d.keq))
    end
    x, cv
end
function hard_search(m; npts = 30, rng)
    _, _, sol, L = derivable(m)
    _, _, vars, terms = kcat_class(m, sol)
    best = 0.0
    for _ in 1:npts
        xx, _ = wide_params(m, sol, rng)
        nl = NLaw(m, L, xx)
        sup, _ = sup_search(nl, rng; mods = false, lo = -6.0, hi = 6.0, nrand = 300 * length(nl.subs))
        best = max(best, sup / kcat_value(vars, terms, xx))
    end
    best
end
ratios = [hard_search(recs[(x.set, x.index)].c.mech; rng = Xoshiro(hash((x.set, x.index, :hard)))) for x in ssub]
@printf("random-order SS without product-holding SS cycles: %d laws, 30 points each (constants 1e-3..1e3, box 1e-6..1e6): sup > 1.01 kcat in %d; max ratio %.4f\n",
        length(ssub), count(>(1.01), ratios), maximum(ratios))
for (x, q) in zip(ssub, ratios)
    q > 1.01 || continue
    println("  overshoot: ", x.id, " (", x.set, "/", x.index, ") max sup/kcat ", round(q; digits = 3), " feat ", x.feat)
    foreach(l -> println("      ", l), step_lines(recs[(x.set, x.index)].c.mech))
end
