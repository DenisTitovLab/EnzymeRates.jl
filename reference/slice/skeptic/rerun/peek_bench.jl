# ABOUTME: Prints the bench rows of a saved run whose zero-concentration check failed or whose
# ABOUTME: time is the largest; usage: peek_bench.jl label
include(joinpath(@__DIR__, "lib.jl"))
rows = deserialize(joinpath(IDATA, "bench_$(ARGS[1]).jls"))
for x in filter(x -> !x.zeros_ok, rows)
    println("zeros: ", x.set, " ", x.id, " zero_err ", x.zero_err, " err ", x.err)
end
for x in first(sort(rows; by = x -> -x.t1n_ns), 5)
    println("slow: ", x.set, " ", x.id, " ", round(x.t1n_ns; digits = 1), " ns, ops ", x.ops)
end
