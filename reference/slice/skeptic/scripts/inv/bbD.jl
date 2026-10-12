# ABOUTME: Crash probe without the slice: load the pipeline rows as the timing script does, keep a
# ABOUTME: subset referenced, then allocate heavily (generic Julia work) and force a full GC.
include(joinpath(@__DIR__, "..", "indep.jl"))
function work(n)
    acc = 0
    for i in 1:n
        d = Dict{Vector{Pair{Int, Rational{Int}}}, Int}()
        for j in 1:50
            d[[k => (j // (k + 1)) for k in 1:(j % 7 + 1)]] = j
        end
        acc += length(d)
    end
    acc
end
function main()
    rows = [x for k in 1:11 for x in deserialize(joinpath(IDATA_IN, "pipe_terter$k.jls")) if x.ok]
    top = sort(rows; by = x -> -x.nD)[1:1]
    rng = Xoshiro(3)
    big = filter(x -> 1001 <= x.nD <= 5000, rows)
    pick = [top; big[randperm(rng, length(big))[1:1]]]
    rows = nothing
    GC.gc()
    println("work ", work(get(ENV, "WORK", "200000") |> x -> parse(Int, x)))
    GC.gc()
    println("full GC ok; picked ", length(pick))
end
main()
