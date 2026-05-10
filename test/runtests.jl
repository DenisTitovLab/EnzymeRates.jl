using Test
using EnzymeRates
using LinearAlgebra
using Random

include("mechanism_definitions_for_test_enzyme_derivation.jl")
@testset "EnzymeRates.jl" begin
    include("test_accessors.jl")
    include("test_types.jl")
    include("test_dsl.jl")
    include("test_sym_poly.jl")
    include("test_rate_eq_derivation.jl")
    include("test_fitting.jl")
    include("test_mechanism_enumeration.jl")
    include("test_identify_rate_equation.jl")
    include("test_eq_hash_dedup.jl")
    include("test_readme_runs.jl")
    include("test_aqua_jet.jl")
end
