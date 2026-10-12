# ABOUTME: Engine-stage check (g): derivation time and term counts on the 806 shared plain
# ABOUTME: mechanisms against the package's (export_pkg.jl's timings, this machine).
include(joinpath(@__DIR__, "lib.jl"))

const SHARED = ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi", "init_bibi_pp",
                "expand_bibi")
pkgt = deserialize(joinpath(@__DIR__, "data", "pkg_kcat_times.jls"))
records = [r for r in load_records(SHARED) if !S.is_allosteric(r.c.mech)]
"""The number of simple spanning trees of `m`'s segment graph (what the engine expands)."""
function tree_count(m)
    node, G = S.segment_nodes(m, trues(length(m.edges)))
    length(S.spanning_trees(G, unique(minmax(node[e.a], node[e.b]) for e in m.edges
                                      if !e.re && node[e.a] != node[e.b])))
end
best(f; reps = 5) = (f(); minimum(@elapsed(f()) for _ in 1:reps) * 1e6)
rows = NamedTuple[]
for r in records
    m = r.c.mech
    sol = S.solve(m)
    L = S.law(m, sol)
    b = S.bare_law(m)
    t_solve = best(() -> S.solve(m))
    t_law = best(() -> S.law(m, sol))
    t_print = best(() -> (S.rate_equation_string(m, L, sol); S.coefficient_form(m, L)); reps = 2)
    pt = get(pkgt, r.p.id, nothing)
    push!(rows, (id = r.p.id, t_solve, t_law, t_print,
                 pkg_raw = pt === nothing ? NaN : pt.t_polys / 1e3,
                 pkg_dep = pt === nothing ? NaN : pt.t_dep / 1e3,
                 bare_num = S.nterms(b.num), bare_den = S.nterms(b.den), nQ = S.nterms(L.Q),
                 nD = S.nterms(L.D), pkg_num = length(r.p.pkg_num),
                 pkg_den = length(r.p.pkg_den), cx = Int(S.complexity(m)),
                 trees = tree_count(m)))
end
serialize(joinpath(@__DIR__, "data", "cost_rows.jls"), rows)
q(v, p) = sort(v)[clamp(round(Int, p * length(v)), 1, length(v))]
for (lab, f) in (("solve", r -> r.t_solve), ("law (engine)", r -> r.t_law),
                 ("solve + law", r -> r.t_solve + r.t_law), ("printing (both forms)", r -> r.t_print),
                 ("package raw polys", r -> r.pkg_raw),
                 ("package dependents", r -> r.pkg_dep),
                 ("package raw polys + dependents", r -> r.pkg_raw + r.pkg_dep))
    v = filter(!isnan, f.(rows))
    @printf("%-32s n %3d  median %8.1f µs  p90 %8.1f  max %9.1f\n", lab, length(v), q(v, 0.5),
            q(v, 0.9), maximum(v))
end
ratio = [(r.t_solve + r.t_law) / (r.pkg_raw + r.pkg_dep) for r in rows if !isnan(r.pkg_raw)]
@printf("(solve + law) / (package raw polys + dependents): median %.2f, p90 %.2f, max %.2f\n",
        q(ratio, 0.5), q(ratio, 0.9), maximum(ratio))
@printf("terms: package numerator %d, denominator %d; mine Q %d, D %d, Q + D %d (package %d)\n",
        sum(r.pkg_num for r in rows), sum(r.pkg_den for r in rows), sum(r.nQ for r in rows),
        sum(r.nD for r in rows), sum(r.nQ + r.nD for r in rows),
        sum(r.pkg_num + r.pkg_den for r in rows))
@printf("D terms equal to the package's denominator terms in %d of %d; fewer in %d, more in %d\n",
        count(r -> r.nD == r.pkg_den, rows), length(rows), count(r -> r.nD < r.pkg_den, rows),
        count(r -> r.nD > r.pkg_den, rows))
@printf("Q terms vs half the package's numerator terms: equal %d, fewer %d, more %d\n",
        count(r -> 2r.nQ == r.pkg_num, rows), count(r -> 2r.nQ < r.pkg_num, rows),
        count(r -> 2r.nQ > r.pkg_num, rows))
big3 = sort(rows; by = r -> -r.nD)[1:5]
for r in big3
    @printf("  %-16s law %7.1f µs solve %6.1f µs | Q %4d D %4d | package num %4d den %4d | \
            bare num %5d den %5d | complexity %d, simple trees %d\n", r.id, r.t_law, r.t_solve,
            r.nQ, r.nD, r.pkg_num, r.pkg_den, r.bare_num, r.bare_den, r.cx, r.trees)
end
@printf("largest complexity %d (cap %d on trees x segments); largest simple tree count %d\n",
        maximum(r.cx for r in rows), S.MAX_TERMS, maximum(r.trees for r in rows))
