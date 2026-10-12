# ABOUTME: Extra statistics for the P7 report: Sig size vs law encoding size, op-count ratios to flat,
# ABOUTME: and compile/JIT totals for the largest 60 laws.
using Serialization, Statistics, Printf
D = joinpath(@__DIR__, "data")
P = deserialize(joinpath(D, "bench_pkg_unique.jls"))
sig = [x.sig_bytes for x in P]
@printf("package Sig (summarysize of the mechanism type): median %d B, max %d B\n", median(sig), maximum(sig))
for v in ("flat", "facth", "halh")
    R = deserialize(joinpath(D, "bench_$(v)_unique.jls"))
    ints = [x.ints for x in R]
    @printf("%-6s law code: median %d Ints (%.1f KB as Int64, %.1f KB as Int16), max %d Ints (%.1f KB / %.1f KB)\n",
            v, median(ints), 8median(ints) / 1024, 2median(ints) / 1024, maximum(ints), 8maximum(ints) / 1024, 2maximum(ints) / 1024)
end
F = Dict(x.id => x for x in deserialize(joinpath(D, "bench_flat_unique.jls")))
big = Set(x.id for x in P if x.nstmts > 900)
for v in ("flatr", "flath", "fact", "facth", "hal", "halh", "segh")
    R = deserialize(joinpath(D, "bench_$(v)_unique.jls"))
    rs = [(x.t1_ns / F[x.id].t1_ns, x.compile_s / F[x.id].compile_s, x.jit / F[x.id].jit, x.nstmts / F[x.id].nstmts) for x in R if x.id in big]
    @printf("%-6s vs flat on the 60 large laws: T1 %.2f, compile %.2f, JIT %.2f, stmts %.2f (medians)\n", v,
            median(first.(rs)), median(getindex.(rs, 2)), median(getindex.(rs, 3)), median(last.(rs)))
end
for v in ("pkg", "flat", "flath", "facth", "halh", "segh")
    R = deserialize(joinpath(D, "bench_$(v)_unique.jls"))
    L = [x for x in R if x.id in big]
    @printf("%-6s large 60: compile sum %.2f s, JIT sum %.0f KB, T1 max %.1f ns; specs 42: compile sum %.2f s, JIT sum %.0f KB, T1 max %.1f ns\n",
            v, sum(x.compile_s for x in L), sum(x.jit for x in L) / 1024, maximum(x.t1_ns for x in L),
            sum(x.compile_s for x in R if startswith(x.id, "specs")), sum(x.jit for x in R if startswith(x.id, "specs")) / 1024,
            maximum(x.t1_ns for x in R if startswith(x.id, "specs")))
end
