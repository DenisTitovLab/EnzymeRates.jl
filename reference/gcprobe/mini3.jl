# ABOUTME: Standalone reproducer (no Slice): main_b's shape — a loop whose body passes a
# ABOUTME: Union{L{Int}, L{BigInt}} result's field to a function, then runs a full GC.
struct P{T}
    terms::Dict{Int, T}
end
struct L{T}
    D::P{T}
end
nterms(p::P) = length(p.terms)
@noinline make(flag::Bool, n) =
    flag ? L{BigInt}(P{BigInt}(Dict{Int, BigInt}(i => big(i) for i in 1:n))) :
           L{Int}(P{Int}(Dict{Int, Int}(i => i for i in 1:n)))
@noinline churn(n) = sum(length, [Vector{Int}(undef, 8) for _ in 1:n])   # young garbage

function main()
    for r in 1:2
        println(r, " ", nterms(make(false, 1000).D))
        churn(200_000)
        println("full gc"); flush(stdout)
        GC.gc(true)
        println("full gc ok"); flush(stdout)
    end
end
get(ENV, "IR", "0") == "1" || main()
