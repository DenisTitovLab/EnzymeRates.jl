# ABOUTME: Standalone attempt at the F6 shape without Slice: a non-inlined function returns
# ABOUTME: Union{L{Int}, L{BigInt}} of a 9-pointer struct built by an allocation-heavy callee.
struct P{T}; terms::Dict{Vector{Int}, T}; end
struct V; names::Vector{Symbol}; conc::BitVector; index::Dict{Symbol, Int}; end
struct L{T}
    vars::V; H::P{T}; Q::P{T}; D::P{T}; sigma::Vector{P{T}}; trees::Vector{P{T}}; shift::Vector{Int}
end
mkp(::Type{T}, n, s) where {T} = P{T}(Dict{Vector{Int}, T}([i, s, i % 7] => T(i) for i in 1:n))
@noinline function build(::Type{T}, n) where {T}
    vars = V([Symbol("x", i) for i in 1:8], trues(8), Dict(Symbol("x", i) => i for i in 1:8))
    H = mkp(T, 2, 0); Q = mkp(T, n, 1); D = mkp(T, n, 2)
    sigma = [mkp(T, n ÷ 4, k) for k in 1:4]
    trees = [mkp(T, n ÷ 4, k) for k in 5:8]
    L{T}(vars, H, Q, D, sigma, trees, collect(1:3))
end
@noinline lawT(::Type{T}, n) where {T} = build(T, n)
@noinline w(n, flag::Bool) = flag ? lawT(BigInt, n) : lawT(Int, n)
check(x::L) = length(x.Q.terms) + length(x.D.terms) + sum(p -> length(p.terms), x.sigma) +
              sum(p -> length(p.terms), x.trees) + length(x.vars.names)
function main(N, n)
    keep = Any[]
    acc = 0
    for i in 1:N
        x = w(n, false)
        acc += check(x)
        push!(keep, x)
        length(keep) > 20 && popfirst!(keep)
        i % 50 == 0 && (GC.gc(true); println(i, " ", acc); flush(stdout))
    end
    println("finished ", acc)
end
main(parse(Int, get(ENV, "N", "2000")), parse(Int, get(ENV, "NT", "3000")))
