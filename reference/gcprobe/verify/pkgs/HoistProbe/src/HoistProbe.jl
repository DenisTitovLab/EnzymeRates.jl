# ABOUTME: Positive control for the remark scan: the repro_upstream.jl pattern inside a package
# ABOUTME: whose `main` is precompiled, to test whether --pkgimages=no exposes its remarks.
module HoistProbe
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
function main(rounds, flag::Bool)
    for r in 1:rounds
        println(r, " ", nterms(make(flag, 1000).D)); flush(stdout)
        GC.gc(true)
        println("full gc ", r, " ok"); flush(stdout)
    end
    println("finished")
end
precompile(main, (Int, Bool))
end
