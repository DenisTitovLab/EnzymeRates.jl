using Test
using EnzymeRates
using LinearAlgebra
using Random

include("test_helpers.jl")

@testset "EnzymeRates.jl" begin
    include("test_accessors.jl")
    include("test_types.jl")
    include("test_dsl.jl")
    include("test_haldane_wegscheider.jl")
    include("test_mechanisms.jl")
    include("test_ode_steadystate.jl")
    include("test_identifiability.jl")
    include("test_identifiability_validation.jl")
    include("test_aqua_jet.jl")
end
