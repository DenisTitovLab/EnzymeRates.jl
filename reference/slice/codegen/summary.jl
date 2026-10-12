# ABOUTME: Summarizes the codegen-stage benchmarks (mine and the package's, matched by record
# ABOUTME: id): allocations, ns per call, compile time and native code, medians and maxima.
using Serialization, Printf, Statistics
const DATA = joinpath(@__DIR__, "data")
loadrows(prefix) = reduce(vcat, [deserialize(joinpath(DATA, f)) for f in readdir(DATA)
                                 if startswith(f, prefix)]; init = NamedTuple[])
mine = loadrows("bench_") |> rs -> filter(r -> !haskey(r, :exported_ns), rs)
pkg = loadrows("bench_pkg_")
q(v, p) = isempty(v) ? NaN : sort(v)[clamp(ceil(Int, p * length(v)), 1, length(v))]
stats(v) = @sprintf("median %7.1f  p90 %7.1f  max %7.1f", q(v, 0.5), q(v, 0.9), maximum(v))
const SHARED = ("specs", "handwritten", "init_uniuni", "init_unibi", "init_bibi",
                "init_bibi_pp", "expand_bibi", "allo_uniuni", "allo_bibi")

println("== mine: $(length(mine)) laws ($(count(r -> !r.allo, mine)) plain, ",
        "$(count(r -> r.allo, mine)) MWC); allocations > 0: $(count(r -> r.allocs > 0, mine))")
for (lab, sel) in (("plain", r -> !r.allo), ("MWC", r -> r.allo))
    rs = filter(sel, mine)
    fresh = filter(r -> r.jit > 0, rs)         # an equal mechanism seen earlier compiles nothing
    println("  $lab ($(length(rs)) laws, $(length(fresh)) compiled fresh):")
    for (name, f) in (("T1 ns", r -> r.t1_ns), ("T1n ns", r -> r.t1n_ns), ("T2 ns", r -> r.t2_ns))
        println("    ", rpad(name, 10), stats(f.(rs)), "   over 120 ns: ", count(r -> f(r) > 120, rs))
    end
    println("    ", rpad("compile ms", 10), stats(1e3 .* getfield.(fresh, :compile_s)))
    println("    ", rpad("native KB", 10), stats(getfield.(fresh, :jit) ./ 1024))
    println("    ", rpad("gen ms", 10), stats(1e3 .* getfield.(rs, :t_gen)))
    println("    ", rpad("ops", 10), stats(Float64.(getfield.(rs, :ops))))
end
for id in ("expand_bibi/626", "expand_bibi/642", "expand_bibi/658")
    r = only(filter(x -> x.id == id, mine))
    p = only(filter(x -> x.id == id, pkg))
    @printf("  %s: mine T1 %.1f T1n %.1f T2 %.1f ns, compile %.0f ms, native %.1f KB, ops %d | \
            package T1 %.1f T1n %.1f T2 %.1f ns, compile %.0f ms, native %.1f KB, ops %d\n",
            id, r.t1_ns, r.t1n_ns, r.t2_ns, 1e3 * r.compile_s, r.jit / 1024, r.ops, p.t1_ns,
            p.t1n_ns, p.t2_ns, 1e3 * p.compile_s, p.jit / 1024, p.ops)
end
worst = sort(mine; by = r -> -r.t1n_ns)[1:5]
println("  slowest T1n: ", join(("$(r.id) $(round(r.t1n_ns; digits = 1))" for r in worst), ", "))

println("== package: $(length(pkg)) laws; allocations > 0: $(count(r -> r.allocs > 0, pkg))")
bykey = Dict(r.id => r for r in mine if r.set in SHARED)
pairs_ = [(bykey[p.id], p) for p in pkg if haskey(bykey, p.id)]
println("  matched with mine: $(length(pairs_))")
for (lab, sel) in (("plain", x -> !x[2].allo), ("MWC", x -> x[2].allo))
    ps = filter(sel, pairs_)
    isempty(ps) && continue
    fresh = filter(x -> x[1].jit > 0 && x[2].jit > 0, ps)
    println("  $lab ($(length(ps)) matched):")
    for (name, f) in (("T1 ns", r -> r.t1_ns), ("T1n ns", r -> r.t1n_ns), ("T2 ns", r -> r.t2_ns))
        m, k = [f(a) for (a, _) in ps], [f(b) for (_, b) in ps]
        rat = [f(a) / f(b) for (a, b) in ps if f(b) > 0.5]
        println("    ", rpad(name, 8), "mine ", stats(m), " | package ", stats(k),
                " | over 120 ns mine $(count(>(120), m)) package $(count(>(120), k))")
        println("    ", rpad("", 8), "ratio mine/package ", stats(rat))
    end
    for (name, f) in (("compile ms", r -> 1e3 * r.compile_s), ("native KB", r -> r.jit / 1024))
        m, k = [f(a) for (a, _) in fresh], [f(b) for (_, b) in fresh]
        rat = [f(a) / f(b) for (a, b) in fresh]
        println("    ", rpad(name, 10), "mine ", stats(m), " | package ", stats(k))
        @printf("    %s ratio mine/package %s; totals mine %.1f package %.1f\n", rpad("", 10),
                stats(rat), sum(m), sum(k))
    end
end
