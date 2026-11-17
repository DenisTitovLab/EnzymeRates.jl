# ABOUTME: Top-level test entry point that wires up shared fixtures and
# ABOUTME: includes every EnzymeRates test file.
using Test
using EnzymeRates
using LinearAlgebra
using Random

include("mechanism_definitions_for_test_enzyme_derivation.jl")

# Throwaway diagnosis: time, garbage collection and memory of each test file.
_testhelper_gib(bytes) = round(bytes / 2^30; digits = 2)
function _testhelper_timed_include(file)
    run = @timed include(file)
    @info "TIMING $file: $(round(run.time; digits = 1)) s, " *
          "$(round(run.gctime; digits = 1)) s in GC, " *
          "max RSS $(_testhelper_gib(Sys.maxrss())) GiB, " *
          "live heap $(_testhelper_gib(Base.gc_live_bytes())) GiB, " *
          "free $(_testhelper_gib(Sys.free_memory())) of " *
          "$(_testhelper_gib(Sys.total_memory())) GiB"
end
@info "TIMING start: total memory $(_testhelper_gib(Sys.total_memory())) GiB, " *
      "free $(_testhelper_gib(Sys.free_memory())) GiB, threads $(Threads.nthreads()), " *
      "CPU threads $(Sys.CPU_THREADS)"
@testset verbose = true "EnzymeRates.jl" begin
    _testhelper_timed_include("test_compile_budget.jl")
    _testhelper_timed_include("Mechanisms/test_accessors.jl")
    _testhelper_timed_include("Mechanisms/test_graph_accessors.jl")
    _testhelper_timed_include("Mechanisms/test_types.jl")
    _testhelper_timed_include("Mechanisms/test_reaction_graph.jl")
    _testhelper_timed_include("Mechanisms/test_segments.jl")
    _testhelper_timed_include("DSL/test_dsl.jl")
    _testhelper_timed_include("Derivation/test_rate_eq_derivation.jl")
    _testhelper_timed_include("Derivation/test_poly.jl")
    _testhelper_timed_include("Derivation/test_law.jl")
    _testhelper_timed_include("Derivation/test_key.jl")
    _testhelper_timed_include("Derivation/test_print.jl")
    _testhelper_timed_include("Derivation/test_oracles.jl")
    _testhelper_timed_include("Derivation/test_admission.jl")
    _testhelper_timed_include("Derivation/test_joint_solve.jl")
    _testhelper_timed_include("Derivation/test_names.jl")
    _testhelper_timed_include("Derivation/test_conformations.jl")
    _testhelper_timed_include("Derivation/test_pivot_priority_regression.jl")
    _testhelper_timed_include("Derivation/test_allosteric_collapse.jl")
    _testhelper_timed_include("Derivation/test_allosteric_golden.jl")
    _testhelper_timed_include("Derivation/allosteric_ground_truth.jl")
    _testhelper_timed_include("Data/test_rate_data.jl")
    _testhelper_timed_include("Fitting/test_fitting.jl")
    _testhelper_timed_include("Fitting/test_fitting_scale.jl")
    _testhelper_timed_include("Fitting/test_fit_search.jl")
    _testhelper_timed_include("Enumeration/test_mechanism_enumeration.jl")
    _testhelper_timed_include("test_identify_rate_equation.jl")
    if v"1.13.0-" <= VERSION < v"1.14.0-"
        _testhelper_timed_include("test_type_stability.jl")
    else
        @test_skip "zero JET reports is required on Julia 1.13 only"
    end
    _testhelper_timed_include("test_aqua_jet.jl")
    _testhelper_timed_include("test_boundaries.jl")
end
