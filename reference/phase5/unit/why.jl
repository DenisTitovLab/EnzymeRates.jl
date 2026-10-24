# ABOUTME: Classifies why the spec rule's N differs from Cleland's in the census rows: the basis
# ABOUTME: holds a lower-degree substrate term, a square, or a product term; args: pops.
using Serialization, Printf
include(joinpath(@__DIR__, "ulib.jl"))

const DATA = joinpath(@__DIR__, "data")

"""The basis monomials of a `lamtext` string (`A*B:1//1 C:-1//1`) as name lists."""
basis(t) = [split(first(split(p, ":")), "*") for p in split(t)]

for pop in ARGS
    rows = reduce(vcat, [deserialize(joinpath(DATA, f)) for f in readdir(DATA)
                         if startswith(f, pop * "_") && endswith(f, ".jls")]; init = Any[])
    C = [(r.id, x) for r in rows if r.ok for (_, x) in r.rules
         if !x.constant && x.cleland && x.same_cleland !== true]
    n = Dict{String, Int}()
    for (id, x) in C
        b = basis(x.spec_lambda)
        clel = basis(x.cleland_lambda)
        nsub = maximum(length, clel)              # n, the degree of ΠS
        low = any(m -> length(m) < nsub - 1 && !any(occursin("^", s) for s in m), b)
        sq = any(m -> any(occursin("^", s) for s in m), b)
        subs = Set(s for m in clel for s in m)
        pr = any(m -> any(s -> !(first(split(s, "^")) in subs), m), b)
        key = join(filter(!isempty, [low ? "lower-degree substrate term" : "",
                                     sq ? "a power" : "", pr ? "a non-substrate" : ""]), " + ")
        key = isempty(key) ? "other" : key
        n[key] = get(n, key, 0) + 1
    end
    println(pop, ": ", length(C), " differ from Cleland; ",
            join(("$k $v" for (k, v) in sort!(collect(n))), ", "))
    lens = Dict{String, Int}()
    for (id, x) in C
        k = string(length(basis(x.spec_lambda)), " terms in λ")
        lens[k] = get(lens, k, 0) + 1
    end
    println("   ", join(("$k: $v" for (k, v) in sort!(collect(lens))), ", "))
end

# Under Cleland-first only laws without every Cleland monomial reach the graded basis: how
# often each order variant changes N there.
for pop in ARGS
    rows = reduce(vcat, [deserialize(joinpath(DATA, f)) for f in readdir(DATA)
                         if startswith(f, pop * "_") && endswith(f, ".jls")]; init = Any[])
    G = [x for r in rows if r.ok for (_, x) in r.rules if !x.constant && !x.cleland]
    println(pop, ": ", length(G), " no-constant conformations without Cleland's monomials; ",
            join(("$o changes N in $(count(x -> x.same_order[o] !== true, G))" for o in ORDERS), ", "))
end
