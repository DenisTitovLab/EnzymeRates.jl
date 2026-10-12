# ABOUTME: Tests whether det of a BigInt matrix built like complexity()'s (zeros, then +=) mutates
# ABOUTME: shared BigInt objects in place (aliasing), a candidate cause of the GC crash.
using LinearAlgebra
A = zeros(BigInt, 3, 3)
println("zeros(BigInt) cells alias one object: ", A[1, 1] === A[2, 2], "; zero(BigInt) cached: ",
        zero(BigInt) === zero(BigInt))
println("det method: ", which(det, (Matrix{BigInt},)))
# the Laplacian of a 4-cycle plus a chord, built as complexity() builds it
G = 5
lap = zeros(BigInt, G, G)
for (a, b) in ((1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (1, 3))
    lap[a, a] += 1; lap[b, b] += 1; lap[a, b] -= 1; lap[b, a] -= 1
end
before = deepcopy(lap)
M = lap[2:end, 2:end]
z0 = A[3, 3]
d = det(M)
println("det = ", d, " (expected spanning trees of that graph)")
println("lap unchanged after det: ", lap == before, "; untouched zero cell still 0: ", z0, " ", A[3, 3])
println("zero(BigInt) after: ", zero(BigInt), "; big(0) = ", big(0))
# repeated: many Laplacians, compare det against Rational elimination
using Random
rng = Xoshiro(1)
bad = 0
for t in 1:2000
    n = rand(rng, 3:9)
    L = zeros(BigInt, n, n)
    for _ in 1:rand(rng, n:3n)
        a, b = rand(rng, 1:n), rand(rng, 1:n)
        a == b && continue
        L[a, a] += 1; L[b, b] += 1; L[a, b] -= 1; L[b, a] -= 1
    end
    d1 = det(L[2:end, 2:end])
    d2 = round(BigInt, det(Rational{BigInt}.(L[2:end, 2:end])))
    global bad += d1 != d2
end
println("2000 random Laplacians: det(BigInt) differs from the rational determinant in ", bad)
GC.gc(); GC.gc()
println("GC ok")
