# ABOUTME: Tests for beam-based mechanism enumeration runtime functions.
# ABOUTME: Validates parameter counting and fingerprinting against compiled mechanisms.

@testset "Beam Enumeration" begin
    @testset "Runtime functions" begin
        @testset "_runtime_param_count matches @generated" begin
            for (name, rxn) in [("uni-uni", uni_uni),
                                ("uni-bi", uni_bi),
                                ("bi-bi", bi_bi),
                                ("bi-bi pp", bi_bi_ping_pong)]
                topos = EnzymeRates._catalytic_topologies(rxn)
                for spec in topos
                    runtime_pc = EnzymeRates._runtime_param_count(spec)
                    m = compile_mechanism(spec)
                    generated_pc = length(parameters(m))
                    @test runtime_pc == generated_pc
                end
            end
        end

        @testset "_runtime_denominator_monomials" begin
            topos = EnzymeRates._catalytic_topologies(bi_bi)
            for spec in topos
                monos = EnzymeRates._runtime_denominator_monomials(spec)
                @test !isempty(monos)
            end

            spec1 = topos[1]
            fp1 = EnzymeRates._runtime_denominator_monomials(spec1)
            fp2 = EnzymeRates._runtime_denominator_monomials(spec1)
            @test fp1 == fp2

            # Verify runtime fingerprint matches old _concentration_fingerprint
            for spec in topos
                runtime_fp = EnzymeRates._runtime_denominator_monomials(spec)
                partition = EnzymeRates._compute_re_partition_from_steps(spec.steps)
                old_fp = EnzymeRates._concentration_fingerprint(spec.steps, partition)
                @test runtime_fp == old_fp
            end
        end
    end
end
