# ABOUTME: Tests for beam-based mechanism enumeration runtime functions.
# ABOUTME: Validates parameter counting and fingerprinting against compiled mechanisms.

@testset "Beam Enumeration" begin
    @testset "expand_mechanisms_by_one_param" begin

        @testset "RE→SS move" begin
            topos = EnzymeRates._catalytic_topologies(uni_uni)
            spec = topos[1]
            n_re = count(s -> s.is_equilibrium, spec.steps)
            @test n_re == 2

            results = EnzymeRates.expand_mechanisms_by_one_param(
                [spec], uni_uni)
            re_to_ss = filter(
                r -> r.param_count == spec.param_count + 1, results)
            @test length(re_to_ss) == n_re

            for r in re_to_ss
                r_n_re = count(s -> s.is_equilibrium, r.steps)
                @test r_n_re == n_re - 1
                @test r.param_count == spec.param_count + 1
            end

            topos_bb = EnzymeRates._catalytic_topologies(bi_bi)
            spec_bb = topos_bb[1]
            n_re_bb = count(s -> s.is_equilibrium, spec_bb.steps)
            @test n_re_bb == 4

            results_bb = EnzymeRates.expand_mechanisms_by_one_param(
                [spec_bb], bi_bi)
            re_to_ss_bb = filter(
                r -> r.param_count == spec_bb.param_count + 1, results_bb)
            @test length(re_to_ss_bb) == n_re_bb

            for r in re_to_ss_bb
                m = compile_mechanism(r)
                @test length(parameters(m)) == r.param_count
            end
        end

    end  # expand_mechanisms_by_one_param testset

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
