# ABOUTME: Chooses the ter-ter laws to compile: the largest tenth of init_terter by denominator
# ABOUTME: terms and 2000 others at random, split into batches; usage: select_terter.jl
include(joinpath(@__DIR__, "lib.jl"))
rows = [x for k in 1:11 for x in deserialize(joinpath(IDATA, "pipe_terter$k.jls")) if x.ok]
order = sortperm([x.nD for x in rows]; rev = true)
top = order[1:cld(length(rows), 10)]
rest = setdiff(eachindex(rows), top)
pick = [top; rest[randperm(Xoshiro(20261013), length(rest))[1:2000]]]
@printf("%d derived init_terter laws; largest tenth %d (D terms %d to %d), plus 2000 others\n",
        length(rows), length(top), rows[top[end]].nD, rows[top[1]].nD)
nb = 9
for b in 1:nb
    sel = Dict{String, Vector{Int}}()
    for i in pick[b:nb:end]
        push!(get!(sel, rows[i].set, Int[]), rows[i].index)
    end
    serialize(joinpath(IDATA, "terter_bench_batch$b.jls"), sel)
end
println("batches: ", nb, " of about ", cld(length(pick), nb), " laws")
