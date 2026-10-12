# ABOUTME: Lists the matched laws where my generated code is slower than the package's (T1n,
# ABOUTME: T2), with both timings and op counts, to tell real regressions from timer noise.
using Serialization, Printf
const DATA = joinpath(@__DIR__, "data")
loadrows(prefix) = reduce(vcat, [deserialize(joinpath(DATA, f)) for f in readdir(DATA)
                                 if startswith(f, prefix)]; init = NamedTuple[])
mine = filter(r -> !haskey(r, :exported_ns), loadrows("bench_"))
pkg = loadrows("bench_pkg_")
by = Dict(r.id => r for r in mine)
slow = [(by[p.id], p) for p in pkg if haskey(by, p.id) &&
        (by[p.id].t1n_ns > 1.2 * p.t1n_ns || by[p.id].t2_ns > 1.2 * p.t2_ns)]
println(length(slow), " of ", length(pkg), " matched laws slower by > 20% in T1n or T2")
for (a, b) in sort(slow; by = x -> -x[1].t1n_ns / x[2].t1n_ns)[1:min(end, 25)]
    @printf("%-20s mine T1 %5.1f T1n %5.1f T2 %5.1f ops %4d | package T1 %5.1f T1n %5.1f T2 %5.1f ops %4d\n",
            a.id, a.t1_ns, a.t1n_ns, a.t2_ns, a.ops, b.t1_ns, b.t1n_ns, b.t2_ns, b.ops)
end
big = [(by[p.id], p) for p in pkg if haskey(by, p.id) && p.t1n_ns > 10]
@printf("matched laws with package T1n > 10 ns: %d; mine slower in %d\n", length(big),
        count(x -> x[1].t1n_ns > x[2].t1n_ns, big))
