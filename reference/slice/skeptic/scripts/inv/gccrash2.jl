# ABOUTME: Second crash probe: deserialize the saved pipeline rows (NamedTuples holding BigInts),
# ABOUTME: drop them and force full collections, with and without the slice loaded.
const SCR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
include(joinpath(SCR, "proto", "common", "Population.jl"))
include(joinpath(SCR, "slice", "src", "Slice.jl"))
using Serialization
const NFILES = parse(Int, get(ENV, "NFILES", "11"))
function main()
    d = joinpath(SCR, "slice", "integrate", "data")
    rows = [x for k in 1:NFILES for x in deserialize(joinpath(d, "pipe_terter$k.jls"))]
    println("rows ", length(rows), "; complexity field type ", typeof(rows[1].complexity))
    s = sum(x -> x.ok ? x.complexity : big(0), rows)
    println("sum of complexities ", s)
    rows = nothing
    for i in 1:4
        GC.gc()
        println("full GC $i ok")
    end
end
main()
