# ABOUTME: Shared definitions for the skeptic's standalone variants of repro_upstream.jl:
# ABOUTME: a union-returning maker whose `.D` field is boxed by the caller.
struct P{T}
    terms::Dict{Int, T}
end
struct L{T}
    names::Vector{Symbol}
    D::P{T}
end
@noinline make(flag::Bool, n) =
    flag ? L{Float64}([:a], P(Dict{Int, Float64}(i => i for i in 1:n))) :
           L{Int}([:a], P(Dict{Int, Int}(i => i for i in 1:n)))
nterms(p::P) = length(p.terms)
