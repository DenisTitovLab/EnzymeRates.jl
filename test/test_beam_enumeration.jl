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
            re_to_ss_bb = filter(results_bb) do r
                r.param_count == spec_bb.param_count + 1 &&
                    length(r.steps) == length(spec_bb.steps)
            end
            @test length(re_to_ss_bb) == n_re_bb

            for r in re_to_ss_bb
                m = compile_mechanism(r)
                @test length(parameters(m)) == r.param_count
            end
        end

        @testset "Remove equivalence constraint" begin
            topos = EnzymeRates._catalytic_topologies(uni_uni)
            base = topos[1]

            de_steps = copy(base.steps)
            push!(de_steps, EnzymeRates.StepSpec([:E, :I], [:E_I], true))
            push!(de_steps, EnzymeRates.StepSpec([:E_S, :I], [:E_I_S], true))
            push!(de_steps, EnzymeRates.StepSpec([:E_I, :S], [:E_I_S], true))

            constraint = (Symbol("K5"), 1, [(Symbol("K4"), 1)])
            constraints = [constraint]

            spec_with_constraint = EnzymeRates.MechanismSpec(
                uni_uni_dead_end_I, de_steps, constraints, 0)
            pc = EnzymeRates._runtime_param_count(spec_with_constraint)
            spec_with_constraint = EnzymeRates.MechanismSpec(
                uni_uni_dead_end_I, de_steps, constraints, pc)

            results = EnzymeRates.expand_mechanisms_by_one_param(
                [spec_with_constraint], uni_uni_dead_end_I)

            remove_constraint = filter(results) do r
                length(r.param_constraints) <
                    length(spec_with_constraint.param_constraints)
            end

            @test length(remove_constraint) == 1
            @test isempty(remove_constraint[1].param_constraints)
            @test remove_constraint[1].param_count ==
                spec_with_constraint.param_count + 1
        end

        @testset "Add dead-end binding (+1)" begin
            # Uni-uni + dead-end I: capacity=1
            # I can bind to E and E_S (both have 0 and 1
            # bound entity respectively — but E_S has 1
            # and cap=1, so only E is eligible for uni-uni)
            topos = EnzymeRates._catalytic_topologies(uni_uni)
            spec = topos[1]

            results = EnzymeRates.expand_mechanisms_by_one_param(
                [spec], uni_uni_dead_end_I)
            de_results = filter(results) do r
                length(r.steps) > length(spec.steps)
            end

            # At least one dead-end candidate
            @test !isempty(de_results)

            # Each candidate has param_count = spec + 1
            for r in de_results
                @test r.param_count == spec.param_count + 1
            end

            # Verify param_count matches compiled mechanism
            for r in de_results
                m = compile_mechanism(r)
                @test length(parameters(m)) == r.param_count
            end

            # Bi-bi + dead-end I: more forms, more
            # opportunities
            topos_bb = EnzymeRates._catalytic_topologies(bi_bi)
            spec_bb = topos_bb[1]

            results_bb = EnzymeRates.expand_mechanisms_by_one_param(
                [spec_bb], bi_bi_dead_end_I)
            de_results_bb = filter(results_bb) do r
                length(r.steps) > length(spec_bb.steps)
            end

            @test !isempty(de_results_bb)

            for r in de_results_bb
                @test r.param_count == spec_bb.param_count + 1
                m = compile_mechanism(r)
                @test length(parameters(m)) == r.param_count
            end
        end

    end  # expand_mechanisms_by_one_param testset

    @testset "Dead-end binding helpers" begin
        # _bound_entities: track what's bound at each form
        topos = EnzymeRates._catalytic_topologies(bi_bi)
        spec = topos[1]  # ordered bi-bi: E→EA→EAB→EPQ→EQ→E

        bound = EnzymeRates._bound_entities(spec)
        @test bound[:E] == Set{Symbol}()
        @test bound[:E_A] == Set([:A])
        @test bound[:E_A_B] == Set([:A, :B])
        @test bound[:E_P_Q] == Set([:P, :Q])
        @test bound[:E_Q] == Set([:Q])

        # _binding_capacity: max(n_subs, n_prods)
        @test EnzymeRates._binding_capacity(bi_bi) == 2
        @test EnzymeRates._binding_capacity(ter_ter) == 3
        @test EnzymeRates._binding_capacity(uni_uni) == 1
        @test EnzymeRates._binding_capacity(uni_bi) == 2
    end

    @testset "Binding capacity limit" begin
        # Uni-uni: capacity = 1
        # E_S already has S bound (capacity 1), so I cannot
        # bind to E_S for uni-uni (would make 2 bound entities)
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        spec = topos[1]
        results = EnzymeRates.expand_mechanisms_by_one_param(
            [spec], uni_uni_dead_end_I)
        de_results = filter(
            r -> length(r.steps) > length(spec.steps),
            results)

        for r in de_results
            bound = EnzymeRates._bound_entities(r)
            cap = EnzymeRates._binding_capacity(
                uni_uni_dead_end_I)
            for (form, entities) in bound
                @test length(entities) <= cap
            end
        end
    end

    @testset "Multi-level dead-end binding" begin
        # Bi-bi: capacity = 2
        # After adding I to E (creating E_I), further
        # metabolites can bind to E_I since capacity=2
        # allows 2 entities.
        topos = EnzymeRates._catalytic_topologies(bi_bi)
        spec = topos[1]

        # First expansion: add I to get mechanisms with E_I
        results1 = EnzymeRates.expand_mechanisms_by_one_param(
            [spec], bi_bi_dead_end_I)
        de_with_I = filter(results1) do r
            length(r.steps) > length(spec.steps)
        end
        @test !isempty(de_with_I)

        # Second expansion: expand the I-containing
        # mechanisms. Some should have further binding to
        # dead-end forms.
        results2 = EnzymeRates.expand_mechanisms_by_one_param(
            de_with_I, bi_bi_dead_end_I)
        multi_level = filter(results2) do r
            length(r.steps) > maximum(
                length(d.steps) for d in de_with_I)
        end
        # Multi-level dead-end forms should exist for bi-bi
        @test !isempty(multi_level)

        # All forms respect capacity
        cap = EnzymeRates._binding_capacity(bi_bi_dead_end_I)
        for r in multi_level
            bound = EnzymeRates._bound_entities(r)
            for (_, entities) in bound
                @test length(entities) <= cap
            end
        end
    end

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

    @testset "Deduplication within levels" begin
        # Two different catalytic topologies for bi-bi that
        # produce equivalent rate equations should deduplicate
        topos = EnzymeRates._catalytic_topologies(bi_bi)
        @test length(topos) >= 2

        # Expand all by one param (RE→SS)
        all_expanded = EnzymeRates.expand_mechanisms_by_one_param(
            topos, bi_bi)

        # Deduplication should reduce count (strict for bi-bi
        # since different topologies produce equivalent expansions)
        deduped = EnzymeRates._deduplicate_specs(
            all_expanded, bi_bi)
        @test length(deduped) < length(all_expanded)

        # All deduped specs should still have valid param_counts
        for spec in deduped
            @test spec.param_count > 0
        end

        # Deduplication is idempotent
        deduped2 = EnzymeRates._deduplicate_specs(
            deduped, bi_bi)
        @test length(deduped2) == length(deduped)
    end

    @testset "expand_mechanisms_same_param_count" begin
        # Ordered bi-bi: E→EA→EAB→EPQ→EQ→E
        # E_Q + A ⇌ E_A_Q where K_A = K_A_catalytic → +0
        topos = EnzymeRates._catalytic_topologies(bi_bi)
        spec = topos[1]  # ordered bi-bi
        original_pc = spec.param_count

        results = EnzymeRates.expand_mechanisms_same_param_count(
            [spec], bi_bi)

        # All results should have same param_count
        for r in results
            @test r.param_count == original_pc
        end

        # Should find at least one +0 dead-end
        @test !isempty(results)

        # Verify param_count matches compiled mechanism
        for r in results
            m = compile_mechanism(r)
            @test length(parameters(m)) == r.param_count
        end

        # Fixed-point: calling again on the union should
        # produce no additional mechanisms
        all_input = vcat([spec], results)
        results2 = EnzymeRates.expand_mechanisms_same_param_count(
            all_input, bi_bi)
        all_specs = vcat(all_input, results2)
        deduped = EnzymeRates._deduplicate_specs(
            all_specs, bi_bi)
        @test length(deduped) == length(
            EnzymeRates._deduplicate_specs(all_input, bi_bi))
    end
end
