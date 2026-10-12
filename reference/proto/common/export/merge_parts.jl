# ABOUTME: Joins the chunk files data/<name>.part<k>.jls of one population into data/<name>.jls.
# ABOUTME: Usage: julia merge_parts.jl <name> <nchunks>; also joins the chunk summaries.
include(joinpath(@__DIR__, "..", "Population.jl"))
name, n = ARGS[1], parse(Int, ARGS[2])
all = Pop.PMech[]
for k in 1:n
    append!(all, Pop.load_population("$name.part$k"))
end
length(unique(m.id for m in all)) == length(all) || error("duplicate ids across chunks")
Pop.save_population(name, all)
open(joinpath(Pop.DATA_DIR, "$name.summary.txt"), "w") do io
    println(io, "population $name: $(length(all)) mechanisms joined from $n chunks")
    for k in 1:n
        println(io, "---- chunk $k ----")
        write(io, read(joinpath(Pop.DATA_DIR, "$name.part$k.summary.txt")))
    end
end
for k in 1:n
    rm(joinpath(Pop.DATA_DIR, "$name.part$k.jls")); rm(joinpath(Pop.DATA_DIR, "$name.part$k.summary.txt"))
end
println("joined $name: $(length(all)) mechanisms")
