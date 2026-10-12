# ABOUTME: The integration stage's result tables from the saved runs: per population pass
# ABOUTME: criteria, normalizer rules and kcat outcomes, and the ter-ter scaling; usage: final_summary.jl
include(joinpath(@__DIR__, "lib.jl"))

load(prefix) = [x for l in readdir(IDATA) if startswith(l, prefix)
                for x in deserialize(joinpath(IDATA, l))]
pipe, bench, keyrows = load("pipe_"), load("bench_"), load("keys_")
classes = load("classes_")

const GROUPS = [("specs", ["specs"]), ("handwritten", ["handwritten"]),
                ("init_uniuni", ["init_uniuni"]), ("init_unibi", ["init_unibi"]),
                ("init_bibi", ["init_bibi"]), ("init_bibi_pp", ["init_bibi_pp"]),
                ("expand_bibi", ["expand_bibi"]), ("allo_uniuni", ["allo_uniuni"]),
                ("allo_bibi", ["allo_bibi"]), ("allo_uniuni_inh", ["allo_uniuni_inh"]),
                ("allo_unibi", ["allo_unibi"]), ("allo_bibi_pp", ["allo_bibi_pp"]),
                ("allo_terbi", ["allo_terbi"]),
                ("P5a hand-built", ["p5a_adv_cases", "p5a_adv_cases2", "p5a_adv_cases3",
                                    "p5a_adv_cases4", "p5a_cases"]),
                ("P1 adversarial", ["p1_adv_pkg", "p1_adv_pkg2", "p1_direct"]),
                ("p1_enum_inh", ["p1_enum_inh"]),
                ("P2 adversarial", ["p2_adv_pkg", "p2_adv_pkg2", "p2_fuzz_dump",
                                    "p2_fuzz_dump_flip", "p2_hand"]),
                ("init_terter", ["init_terter.part$k" for k in 1:11]),
                ("expand_terter", ["expand_terter"])]
e1(x) = isnan(x) ? "-" : @sprintf("%.1e", x)
"""The largest of `f` over `rows`, NaN values left out (NaN when none is left)."""
mx(f, rows) = (v = filter(!isnan, Float64[f(x) for x in rows]); isempty(v) ? NaN : maximum(v))

println("== Pass criteria per population")
@printf("%-16s %7s %7s %5s %8s %4s | %6s %3s %6s %4s %7s %4s %7s | %11s %7s %14s %10s\n",
        "population", "records", "laws", "rej", "oracle", ">30", "comp", "all",
        "maxns", ">120", "valerr", "zero", "bigerr", "2nd", "var", "same/key/nm",
        "ren=/rev")
for (name, sets) in GROUPS
    P = filter(x -> x.set in sets, pipe)
    isempty(P) && continue
    ok = filter(x -> x.ok, P)
    B = filter(x -> x.set in sets, bench)
    K = filter(x -> x.set in sets, keyrows)
    n2 = length(K) + count(x -> haskey(x, :second), B)
    orc = maximum(x -> x.oracle, ok; init = 0.0)
    s(f) = sum(getfield(x, f) for x in K; init = 0)
    @printf("%-16s %7d %7d %5d %8s %4d | %6d %3d %6.1f %4d %7s %4d %7s | %11s %7s %14s %10s\n",
            name, length(P), length(ok), length(P) - length(ok), e1(orc),
            count(x -> x.oracle > 1e-30, ok), length(B), count(x -> x.allocs > 0, B),
            mx(x -> x.t1n_ns, B), count(x -> x.t1n_ns >= 120, B),
            e1(mx(x -> x.cond_err, B)), count(x -> !x.zeros_ok, B), e1(mx(x -> x.bigerr, B)),
            "$(count(x -> !x.second, K) + count(x -> haskey(x, :second) && !x.second, B))/$n2",
            "$(s(:badvar))/$(s(:nvar))", "$(s(:nsame))/$(s(:keyfail))/$(s(:namechange))",
            "$(count(x -> !x.same_order, K))/$(count(x -> !x.rev_key, K))/$(length(K))")
end
println("columns: laws derived; rej = rejected (admission, cap or solve); oracle = worst relative ",
        "error against the BigFloat oracles (>30: above 1e-30); comp = laws compiled; all = ",
        "with allocations; maxns = slowest T1n; valerr = worst Float64 error over the Haldane ",
        "condition number; zero = zero-concentration failures; bigerr = BigFloat inputs ",
        "through the generated code; 2nd = keys that differ in a second process / keys ",
        "recomputed (key checks and compile runs); var = rewritten inputs whose mechanism or ",
        "key differs / tried; same/key/nm = same-model splits, key changes, fitted-name ",
        "changes; ren=/rev = order-preserving / order-reversing renamings that change the ",
        "key, of the records checked")

println("\n== Rejections")
for (name, sets) in GROUPS
    P = filter(x -> x.set in sets && !x.ok, pipe)
    isempty(P) && continue
    kinds = Dict{String, Int}()
    foreach(x -> (k = first(replace(x.msg, r"\{[^}]*\}" => "{…}",
                                    r"E\([^)]*\)" => "E(…)"), 64);
                  kinds[k] = get(kinds, k, 0) + 1), P)
    println("  ", name, ": ", join(("$n × $k…" for (k, n) in kinds), "; "))
end

println("\n== Normalizer rules and kcat outcomes per population")
for (name, sets) in GROUPS
    ok = filter(x -> x.set in sets && x.ok, pipe)
    isempty(ok) && continue
    rules = Dict{Symbol, Int}()
    foreach(x -> foreach(r -> rules[r] = get(rules, r, 0) + 1, x.rules), ok)
    kc = Dict{String, Int}()
    foreach(x -> (k = x.kcat isa String ? replace(first(x.kcat, 48), r"with .*" => "with …") :
                      "a kcat";
                  kc[k] = get(kc, k, 0) + 1), ok)
    println("  ", rpad(name, 16), " rules ", join(("$k $v" for (k, v) in rules), ", "),
            " | kcat: ", join(("$v × $k" for (k, v) in kc), "; "))
end

println("\n== K1 classes (every merge verified through the oracles)")
for c in classes
    @printf("  %-44s laws %6d classes %6d merges %5d worst %8s failures %d spanning %d\n",
            first(c.pool, 44), c.laws, c.classes, c.merges, e1(c.worst), c.failures, c.cross)
    for (v, p) in c.p8
        @printf("      P8 %-8s common %5d: classes mine %5d, P8 %5d; pairs only mine %d, \
                only P8 %d\n", v, p.common, p.mine, p.p8, p.only_mine, p.only_p8)
    end
end

println("\n== Ter-ter scaling (init_terter and expand_terter laws)")
T = filter(x -> startswith(x.set, "init_terter") || x.set == "expand_terter", pipe)
ok = [merge(x, (t_stages = x.t_solve + x.t_law + x.t_key + x.t_print + x.t_code + x.t_kcat,))
      for x in T if x.ok]
timed = filter(x -> !isnan(x.t_stages), ok)
TB = Dict((x.set, x.index) => x for x in bench
          if startswith(x.set, "init_terter") || x.set == "expand_terter")
edges = [0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
@printf("%-11s %6s %6s %6s %8s %7s %7s %7s %7s %7s %7s | %5s %8s %7s %7s %5s\n", "D terms",
        "laws", "trees", "rref", "solve", "law", "key", "print", "code", "kcat", "sum ms", "comp",
        "comp ms", "KB", "ns", ">120")
for i in 1:length(edges) - 1
    v = filter(x -> edges[i] < x.nD <= edges[i + 1], ok)
    isempty(v) && continue
    b = [TB[(x.set, x.index)] for x in v if haskey(TB, (x.set, x.index))]
    med(f) = (w = [1e3 * getfield(x, f) for x in v if !isnan(getfield(x, f))];
              isempty(w) ? NaN : median(w))
    @printf("%5d-%-5d %6d %6d %6d %8.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.2f | %5d %8.1f %7.1f %7.1f %5d\n",
            edges[i] + 1, edges[i + 1], length(v), median(sum(x.trees) for x in v),
            median(prod(x.rref) for x in v), med(:t_solve), med(:t_law), med(:t_key),
            med(:t_print), med(:t_code), med(:t_kcat),
            med(:t_stages), length(b),
            isempty(b) ? NaN : median(1e3 * x.compile for x in b),
            isempty(b) ? NaN : median(x.jit / 1024 for x in b),
            isempty(b) ? NaN : median(x.t1n_ns for x in b), count(x -> x.t1n_ns >= 120, b))
end
println("medians per bin; trees = spanning trees over both masks; rref = rows x columns of ",
        "the joint system; comp = laws compiled in the bin; KB = native code of the first call")
"""Least-squares slope of log(y) on log(x) over positive pairs."""
function slope(xs, ys)
    p = [(log(a), log(b)) for (a, b) in zip(xs, ys) if a > 0 && b > 0]
    mx, my = mean(first.(p)), mean(last.(p))
    sum((a - mx) * (b - my) for (a, b) in p) / sum((a - mx)^2 for (a, _) in p)
end
println("log-log slope against D terms (1 = linear):")
for f in (:t_solve, :t_law, :t_key, :t_print, :t_code, :t_kcat, :t_stages)
    @printf("  %-8s %.2f\n", f, slope([x.nD for x in timed], [getfield(x, f) for x in timed]))
end
@printf("  derive   %.2f (every law, one call of the entry point)\n",
        slope([x.nD for x in ok], [x.t_derive for x in ok]))
bb = [TB[(x.set, x.index)] for x in ok if haskey(TB, (x.set, x.index))]
@printf("  compile  %.2f (against ops), ns per call %.2f (against ops)\n",
        slope([x.ops for x in bb], [x.compile for x in bb]),
        slope([x.ops for x in bb], [x.t1n_ns for x in bb]))
@printf("  largest law: %d D terms, %d ops; slowest derive %.0f ms; slowest timed key %.0f ms\n",
        maximum(x -> x.nD, ok), maximum(x -> x.ops, ok), 1e3 * maximum(x -> x.t_derive, ok),
        1e3 * maximum(x -> x.t_key, timed))
cap = filter(x -> !x.ok && occursin("spanning-tree", x.msg), T)
println("  rejected by the spanning-tree cap: ", length(cap))
# The package's own size guard (V×τ with parallel steps, `_eq_complexity`) on the same laws.
need = Set((x.set, x.index) for x in cap)
cx = [S.complexity(r.c.mech) for set in unique(x.set for x in cap)
      for r in records_of(set) if (set, r.c.index) in need]
@printf("  their V×τ (the package's _eq_complexity): min %d, median %d; above the package's \
        limit of 5000: %d of %d\n", isempty(cx) ? 0 : minimum(cx), isempty(cx) ? 0 : median(cx),
        count(>(5000), cx), length(cx))
