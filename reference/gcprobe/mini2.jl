# ABOUTME: Standalone reproducer (no Slice): a loop union-splits a field of a Union{L{Int},
# ABOUTME: L{BigInt}} result; JuliaLICM hoists the phi's box, whose field store has no barrier.
struct P{T}
    terms::Dict{Int, T}
end
struct L{T}
    D::P{T}
end
@noinline make(flag::Bool, n) =
    flag ? L{BigInt}(P{BigInt}(Dict{Int, BigInt}(i => big(i) for i in 1:n))) :
           L{Int}(P{Int}(Dict{Int, Int}(i => i for i in 1:n)))
@noinline churn(n) = sum(length, [Vector{Int}(undef, 8) for _ in 1:n])   # young garbage

function main()
    for r in 1:4
        x = make(false, 1000)
        n = length(x.D.terms)          # union-split getfield: phi of Union{P{Int}, P{BigInt}}
        churn(200_000)                 # quick collections while the Dict hangs off the box only
        println(r, " ", n); flush(stdout)
        GC.gc(true)
        println("full gc ", r, " ok"); flush(stdout)
    end
end
get(ENV, "IR", "0") == "1" || main()
