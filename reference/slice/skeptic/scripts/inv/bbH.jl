# ABOUTME: Crash probe: deserialized pipeline rows (holding BigInts) kept alive while BigInt
# ABOUTME: determinants run (as complexity() does), then a full GC; KEEP=all|some|none.
include(joinpath(@__DIR__, "..", "indep.jl"))
using LinearAlgebra
function dets(n)
    rng = Xoshiro(1)
    acc = big(0)
    for _ in 1:n
        k = rand(rng, 3:9)
        L = zeros(BigInt, k, k)
        for _ in 1:2k
            a, b = rand(rng, 1:k), rand(rng, 1:k)
            a == b && continue
            L[a, a] += 1; L[b, b] += 1; L[a, b] -= 1; L[b, a] -= 1
        end
        acc += k * det(L[2:end, 2:end])
    end
    acc
end
function main()
    rows = [x for k in 1:11 for x in deserialize(joinpath(IDATA_IN, "pipe_terter$k.jls")) if x.ok]
    keep = get(ENV, "KEEP", "some")
    sub = keep == "all" ? rows : keep == "some" ? rows[1:1000:end] : []
    rows = nothing
    GC.gc()
    println("dets ", dets(parse(Int, get(ENV, "N", "20000"))))
    GC.gc()
    println("full GC ok; kept ", length(sub))
end
main()
