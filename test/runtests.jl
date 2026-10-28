# ABOUTME: Top-level test entry point that wires up shared fixtures and
# ABOUTME: includes every EnzymeRates test file.
using Test
using EnzymeRates
using LinearAlgebra
using Random

include("mechanism_definitions_for_test_enzyme_derivation.jl")
@testset "EnzymeRates.jl" begin
    include("test_compile_budget.jl")
    include("Mechanisms/test_accessors.jl")
    include("Mechanisms/test_graph_accessors.jl")
    include("Mechanisms/test_write.jl")
    include("Mechanisms/test_types.jl")
    include("Mechanisms/test_reaction_graph.jl")
    include("Mechanisms/test_dsl.jl")
    include("Derivation/test_rate_eq_derivation.jl")
    include("Derivation/test_poly.jl")
    include("Derivation/test_law.jl")
    include("Derivation/test_key.jl")
    include("Derivation/test_print.jl")
    include("Derivation/test_oracles.jl")
    include("Derivation/test_admission.jl")
    include("Derivation/test_joint_solve.jl")
    include("Derivation/test_pivot_priority_regression.jl")
    include("Derivation/test_allosteric_collapse.jl")
    include("Derivation/test_allosteric_golden.jl")
    include("Derivation/allosteric_ground_truth.jl")
    include("Data/test_rate_data.jl")
    include("Fitting/test_fitting.jl")
    include("Fitting/test_fitting_scale.jl")
    include("Enumeration/test_mechanism_enumeration.jl")
    include("test_identify_rate_equation.jl")
    include("test_aqua_jet.jl")
    include("test_boundaries.jl")
end
