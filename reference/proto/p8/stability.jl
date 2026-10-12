# ABOUTME: Cross-process stability of the P8 keys: compares the keys two separate Julia runs saved.
# ABOUTME: Reads data/prev/<pop>.jls (earlier process) and data/<pop>.jls (later process).
using Serialization
tot, same = 0, 0
for pop in ("init_bibi_pp", "expand_bibi", "allo_bibi", "allo_uniuni")
    a = deserialize(joinpath(@__DIR__, "data", "prev", "$pop.jls")).recs
    b = deserialize(joinpath(@__DIR__, "data", "$pop.jls")).recs
    @assert [r[:id] for r in a] == [r[:id] for r in b]
    n, s = 0, 0
    for (ra, rb) in zip(a, b), v in (:today, :redesign), kk in (:k1, :k2)
        haskey(ra, (v, kk, :key)) || continue
        n += 1; s += ra[(v, kk, :key)] == rb[(v, kk, :key)]
    end
    println(pop, ": ", s, " of ", n, " keys identical across the two processes")
    global tot += n; global same += s
end
println("all: ", same, " of ", tot)
