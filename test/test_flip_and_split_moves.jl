# ABOUTME: Tests for the RE→SS group-set flip move and the context-bipartition split move:
# ABOUTME: eligibility, gain proofs, minimal sets, reachability, and the once-per-parent counter.
using Test
using EnzymeRates
using LinearAlgebra
using Random

const _bibi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end
const _uni_uni_rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

@testset "_independent_param_count matches fitted_params on the test specs" begin
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism ?
            EnzymeRates.AllostericMechanism(spec.mechanism) :
            EnzymeRates.Mechanism(spec.mechanism)
        @test EnzymeRates._independent_param_count(m) ==
              length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(m)))
    end
end
