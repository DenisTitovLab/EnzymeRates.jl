# ABOUTME: Standalone reproducer (no Slice): main_b's shape — a non-unrolled loop whose body
# ABOUTME: union-splits a field of a Union{L{Int}, L{BigInt}} result, then runs a full GC.
struct P{T}
    terms::Dict{Int, T}
end
struct L{T}
    names::Vector{Symbol}
    H::P{T}
    D::P{T}
    shift::Vector{Int}
end
nterms(p::P) = length(p.terms)
mkp(::Type{T}, n) where {T} = P{T}(Dict{Int, T}(i => T(i) for i in 1:n))
@noinline make(flag::Bool, n) =
    flag ? L{BigInt}([:a], mkp(BigInt, 2), mkp(BigInt, n), [1]) :
           L{Int}([:a], mkp(Int, 2), mkp(Int, n), [1])
@noinline churn(n) = sum(length, [Vector{Int}(undef, 8) for _ in 1:n])   # young garbage

function main(rounds)
    for r in 1:rounds
        println(r, " ", nterms(make(false, 1000).D))
        churn(200_000)
        println("full gc"); flush(stdout)
        GC.gc(true)
        println("full gc ok"); flush(stdout)
    end
end
const ROUNDS = parse(Int, get(ENV, "ROUNDS", "3"))
get(ENV, "IR", "0") == "1" || main(ROUNDS)
