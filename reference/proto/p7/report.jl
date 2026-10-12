# ABOUTME: Tables of the P7 benchmark: per-variant summaries over the measurement set, the largest
# ABOUTME: mechanisms side by side, ratios to the package, and the kernel sharing of the shared runs.
using Serialization, Printf, Statistics
const D = joinpath(@__DIR__, "data")
load(v, m) = (f = joinpath(D, "bench_$(v)_$(m).jls"); isfile(f) ? deserialize(f) : nothing)
med(x) = median(x); q(x, p) = quantile(x, p)
variants = [("pkg", "unique"), ("flat", "unique"), ("flatr", "unique"), ("flath", "unique"),
            ("fact", "unique"), ("facth", "unique"), ("hal", "unique"), ("halh", "unique"), ("segh", "unique")]
R = Dict(v => load(v, m) for (v, m) in variants if load(v, m) !== nothing)
println("== Summary over the measurement set (n per variant), unique kernels (no sharing)")
@printf("%-6s %4s %6s | %7s %7s %7s %4s | %7s %7s | %8s %8s %8s | %7s %7s %8s | %6s %6s\n",
        "var", "n", "allocs", "T1 med", "T1 p90", "T1 max", ">120", "T2 med", "T2 max",
        "comp med", "comp max", "comp sum", "jit med", "jit max", "jit sum", "stmts", "ops")
for (v, _) in variants
    haskey(R, v) || continue
    r = R[v]
    t1 = [x.t1_ns for x in r]; t2 = [x.t2_ns for x in r]; cs = [x.compile_s for x in r]; j = [x.jit for x in r]
    @printf("%-6s %4d %6d | %7.1f %7.1f %7.1f %4d | %7.1f %7.1f | %8.4f %8.3f %8.2f | %7.0f %7d %8d | %6.0f %6.0f\n",
            v, length(r), maximum(x.allocs for x in r), med(t1), q(t1, 0.9), maximum(t1), count(>(120), t1),
            med(t2), maximum(t2), med(cs), maximum(cs), sum(cs), med(j), maximum(j), sum(j),
            med([x.nstmts for x in r]), med([x.nops for x in r]))
end
if haskey(R, "pkg")
    println("\n== Ratios to the package (median and range over mechanisms; <1 is better)")
    P = Dict(x.id => x for x in R["pkg"])
    for (v, _) in variants
        (v == "pkg" || !haskey(R, v)) && continue
        for (lab, f) in (("T1", x -> x.t1_ns), ("T2", x -> x.t2_ns), ("compile", x -> x.compile_s), ("jit", x -> x.jit), ("stmts", x -> x.nstmts))
            rat = [f(x) / f(P[x.id]) for x in R[v] if haskey(P, x.id)]
            big = [f(x) / f(P[x.id]) for x in R[v] if haskey(P, x.id) && P[x.id].nstmts > 900]
            @printf("  %-6s %-8s median %.2f  [%.2f, %.2f]   large (pkg stmts>900, n=%d) median %.2f [%.2f, %.2f]\n",
                    v, lab, med(rat), minimum(rat), maximum(rat), length(big), med(big), minimum(big), maximum(big))
        end
    end
end
println("\n== Distinct large laws side by side (T1 ns / compile ms / jit KB / typed stmts)")
ids = unique([x.id for x in R["flat"] if x.nstmts > 900])
seen = Set{UInt}()
for id in ids
    f = only(filter(x -> x.id == id, R["flat"]))
    hash(f.code[11:end]) in seen && continue
    push!(seen, hash(f.code[11:end]))
    print(rpad(id, 16))
    for (v, _) in variants
        haskey(R, v) || continue
        x = only(filter(y -> y.id == id, R[v]))
        @printf(" | %-5s %5.1f %5.0f %4.0f %5d", v, x.t1_ns, 1000 * x.compile_s, x.jit / 1024, x.nstmts)
    end
    println()
end
println("\n== Law construction time per mechanism (derivation side, warm; s): median / max")
for (v, _) in variants
    (v == "pkg" || !haskey(R, v)) && continue
    t = [x.t_law for x in R[v]]
    @printf("  %-6s %.5f / %.4f\n", v, med(t), maximum(t))
end
if haskey(R, "pkg")
    t = [x.t_gen for x in R["pkg"]]
    @printf("  pkg generator (derivation + Expr, warm): %.5f / %.4f\n", med(t), maximum(t))
end
for v in ("facth", "flat")
    S = load(v, "shared")
    S === nothing && continue
    n = length(unique(x.code[[1:9; 11:end]] for x in S))
    println("\n== Shared run of $v (salt 0): $(length(S)) mechanisms, $n distinct kernels")
end
