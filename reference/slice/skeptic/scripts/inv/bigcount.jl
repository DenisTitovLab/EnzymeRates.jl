# ABOUTME: Counts the laws that needed BigInt (Int overflow retried) per pipeline batch, and shows
# ABOUTME: whether the laws timed before the GC crashes are among them.
include(joinpath(@__DIR__, "..", "indep.jl"))
for f in sort(readdir(IDATA_IN))
    startswith(f, "pipe_") || continue
    rows = filter(x -> x.ok, deserialize(joinpath(IDATA_IN, f)))
    nb = count(x -> x.big, rows)
    @printf("%-24s laws %6d, BigInt laws %5d\n", f, length(rows), nb)
    for x in rows
        x.id in ("expand_bibi/626", "expand_bibi/642", "expand_bibi/658", "init_terter/62921",
                 "init_terter/62932", "init_terter/63069") && println("    ", x.id, " big = ", x.big)
    end
end
