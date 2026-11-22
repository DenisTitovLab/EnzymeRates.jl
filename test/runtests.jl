# ABOUTME: Top-level test entry point that wires up shared fixtures and
# ABOUTME: includes every EnzymeRates test file.
using Test
using EnzymeRates
using LinearAlgebra
using Random

include("testhelpers.jl")
@testset "EnzymeRates.jl" begin
    include("test_compile_budget.jl")
    include("Mechanisms/test_enzyme_reaction_type.jl")
    include("Mechanisms/test_mechanism_type.jl")
    include("DSL/test_dsl.jl")
    include("Derivation/test_rate_equation_derivation.jl")
    include("Derivation/test_rate_equation_dedup_key.jl")
    include("test_oracles.jl")
    include("Derivation/test_rate_equation_correctness_guard.jl")
    include("Derivation/test_allosteric_rate_equation_correctness_guard.jl")
    include("Derivation/test_thermodynamic_constraints.jl")
    include("Derivation/test_memo.jl")
    include("Data/test_rate_data.jl")
    include("Fitting/test_fitting.jl")
    include("Fitting/test_fitting_scale.jl")
    include("Fitting/test_fit_search.jl")
    include("Enumeration/test_mechanism_enumeration.jl")
    include("test_identify_rate_equation.jl")
    include("Plotting/test_plot_rate_equation_on_data.jl")
    include("Plotting/test_plot_rate_equation_on_data_makie.jl")
    if v"1.13.0-" <= VERSION < v"1.14.0-"
        include("test_type_stability.jl")
    else
        @test_skip "zero JET reports is required on Julia 1.13 only"
    end
    include("test_aqua_jet.jl")
    include("test_boundaries.jl")
end
