# ABOUTME: The engine stage's check (f), copied to run on the integrated sources: kcat against the
# ABOUTME: package's _kcat_forward on every plain record that has it; writes only integrate/data.
include(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad", "slice", "engine", "lib.jl"))

"""The package's rule on law `L` at `vals`: the largest ratio of the positive numerator terms
to the denominator terms that share a product-free concentration monomial."""
function pattern_kcat(m, L, vals)
    prods = Set(L.vars.index[x] for x in m.products)
    num = L.H * L.Q
    n = S.split_conc(S.Poly(filter(t -> t.second > 0, num.terms)), L.vars)
    d = S.split_conc(L.D, L.vars)
    ks = [k for k in keys(n) if haskey(d, k) && !any(v -> v.first in prods, k)]
    isempty(ks) ? big(NaN) : maximum(S.evaluate(n[k], vals) / S.evaluate(d[k], vals) for k in ks)
end

records = [r for r in load_records() if !S.is_allosteric(r.c.mech) && !isempty(r.p.samples)]
exported = deserialize(joinpath("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad", "slice", "engine", "data", "pkg_kcat_times.jls"))
rows, t0 = NamedTuple[], time()
for r in records
    m, p = r.c.mech, r.p
    s = p.samples[1]
    x = get(exported, p.id, nothing)
    r.set in ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp",
              "expand_bibi") || (x = nothing)
    pkg_err = x === nothing ? any(e -> e.first === :kcat, p.stage_errors) : !isempty(x.msg)
    pkgk = x === nothing ? s.kcat : x.kcat
    isnan(pkgk) && !pkg_err && continue
    sol = S.solve(m)
    pars = my_params(r, sol, pkg_values(p, s))
    mine, msg, corners = try
        S.kcat(m, pars), "", length(S.kcat_terms(m, sol)[2])
    catch err
        big(NaN), sprint(showerror, err), 0
    end
    L = S.law(m, sol)
    vals = [L.vars.conc[i] ? big(1) : BigFloat(getfield(pars, n))
            for (i, n) in enumerate(L.vars.names)]
    pat = pattern_kcat(m, L, vals)
    core = S.catalytic_core(m)
    whole = try
        vars, ts = S.kcat_terms(m, sol, trues(length(m.edges)))
        xs = [vars.conc[i] ? big(1) : BigFloat(getfield(pars, n))
              for (i, n) in enumerate(vars.names)]
        Float64(maximum(S.evaluate(a, xs) / S.evaluate(b, xs) for (a, b) in ts))
    catch err
        NaN
    end
    x === nothing || isnan(s.kcat) || x.kcat == s.kcat ||
        println("exported kcat differs from the Pop sample at ", p.id)
    push!(rows, (set = r.set, id = p.id, pkg = pkgk, mine = Float64(mine), msg,
                 pat = Float64(pat), pkg_err, deadend = count(!, core), nedges = length(core),
                 corners, whole))
end
serialize(joinpath(@__DIR__, "data", "stage_engine_kcat_rows.jls"), rows)
rel(a, b) = abs(a - b) / max(abs(a), abs(b))
ok = filter(r -> !r.pkg_err, rows)
@printf("%d plain records with a package kcat or a package kcat error (%.0f s)\n",
        length(rows), time() - t0)
@printf("package pattern rule rebuilt on my law reproduces the package: %d of %d (1e-10)\n",
        count(r -> rel(r.pat, r.pkg) < 1e-10, ok), length(ok))
same = count(r -> isempty(r.msg) && rel(r.mine, r.pkg) < 1e-10, ok)
@printf("kcat equal to the package (1e-10): %d of %d\n", same, length(ok))
@printf("    of the equal ones, with dead-end branches removed: %d\n",
        count(r -> isempty(r.msg) && rel(r.mine, r.pkg) < 1e-10 && r.deadend > 0, ok))
nw = count(r -> !isnan(r.whole) && rel(r.whole, r.pkg) >= 1e-10 && isempty(r.msg) &&
           rel(r.mine, r.pkg) < 1e-10, ok)
@printf("    of the equal ones, the limit with the dead-end branches kept differs: %d \
        (0 in %d of them)\n", nw, count(r -> r.whole == 0 && isempty(r.msg) &&
                                         rel(r.mine, r.pkg) < 1e-10, ok))
@printf("    corners (modifier sets) per mechanism: max %d; mechanisms with more than one: %d\n",
        maximum(r.corners for r in rows), count(r -> r.corners > 1, rows))
println("differences:")
for r in ok
    isempty(r.msg) && rel(r.mine, r.pkg) < 1e-10 && continue
    @printf("  %-12s %-16s pkg %.6g mine %.6g ratio %.4f dead-end edges %d/%d  %s\n", r.set,
            r.id, r.pkg, r.mine, r.mine / r.pkg, r.deadend, r.nedges, r.msg)
end
println("package errors (\"no kcat components\"):")
for r in filter(r -> r.pkg_err, rows)
    println("  ", r.set, " ", r.id, ": ", isempty(r.msg) ? "mine $(r.mine)" : r.msg)
end
