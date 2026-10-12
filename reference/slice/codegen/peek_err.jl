# ABOUTME: Lists the records with the largest Float64 errors in the saved MWC value rows.
using Serialization
for f in filter(x -> startswith(x, "allo_values_"), readdir(joinpath(@__DIR__, "data")))
    rows = deserialize(joinpath(@__DIR__, "data", f))
    for r in sort(rows; by = r -> -r.e64)[1:3]
        println(f, ": ", r.set, " ", r.id, " e64 ", r.e64, " oracle ", r.eor, " n ", r.n)
    end
end
