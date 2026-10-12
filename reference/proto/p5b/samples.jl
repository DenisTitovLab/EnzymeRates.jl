# ABOUTME: Prints sample printed MWC laws from the P5b derivation next to the package's v line.
# ABOUTME: Usage: julia samples.jl <population> <id...>
include(joinpath(@__DIR__, "Allo.jl"))
P = Dict(m.id => m for m in Pop.load_population(ARGS[1]))
for id in ARGS[2:end]
    m = P[id]; md = Allo.derive(m)
    println("== $id $(m.source)  case $(Allo.norm_case(md))")
    println(Allo.law_text(md))
    println("-- package v line:\n", split(m.rate_equation_string, "\n")[end])
end
