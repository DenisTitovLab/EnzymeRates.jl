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

    @testset "expand_mechanisms_by_two_params" begin
        # Uni-uni + allosteric R
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        spec = topos[1]

        results = EnzymeRates.expand_mechanisms_by_two_params(
            [spec], uni_uni_allosteric_R)

        @test !isempty(results)

        # All results should be AllostericMechanismSpec
        for r in results
            @test r isa EnzymeRates.AllostericMechanismSpec
        end

        # Each result has exactly one non-TR-equivalent
        # metabolite (by construction)
        for r in results
            t_mets = EnzymeRates._collect_t_state_metabolites(r)
            n_non_equiv = length(t_mets) -
                length(r.tr_equiv_metabolites)
            @test n_non_equiv == 1
        end

        # Each result should compile and have correct param count
        for r in results
            m = compile_mechanism(r)
            compiled_pc = length(parameters(m))
            runtime_pc = EnzymeRates._runtime_param_count(r)
            @test runtime_pc == compiled_pc
        end

        # All allosteric param_counts should exceed base
        for r in results
            @test EnzymeRates._runtime_param_count(r) >
                spec.param_count
        end

        # No results for reaction without allosteric regulators
        results_no_allo = EnzymeRates.expand_mechanisms_by_two_params(
            [spec], uni_uni)
        @test isempty(results_no_allo)

        # Bi-bi + allosteric R1, R2: should produce results
        # with different site partitions
        topos_bb = EnzymeRates._catalytic_topologies(bi_bi)
        results_bb = EnzymeRates.expand_mechanisms_by_two_params(
            topos_bb, bi_bi_allosteric_R1_R2;
            max_catalytic_n=2)
        @test !isempty(results_bb)

        for r in first(results_bb, 5)
            m = compile_mechanism(r)
            @test EnzymeRates._runtime_param_count(r) ==
                length(parameters(m))
        end

        # Multiple site partitions should be represented
        # for 2 regulators (Bell number B(2)=2: {R1,R2}
        # and {R1},{R2})
        site_configs = Set(
            length(r.allosteric_reg_sites)
            for r in results_bb)
        @test length(site_configs) >= 2
    end

    @testset "Remove TR equivalence (+1, allosteric)" begin
        # Create an allosteric mechanism with TR-equivalent mets
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        spec = topos[1]

        allo_results = EnzymeRates.expand_mechanisms_by_two_params(
            [spec], uni_uni_allosteric_R)
        @test !isempty(allo_results)

        # Pick one with TR-equivalent metabolites
        allo_spec = first(
            r for r in allo_results
            if !isempty(r.tr_equiv_metabolites))

        n_equiv = length(allo_spec.tr_equiv_metabolites)
        @test n_equiv >= 1

        # expand_mechanisms_by_one_param should handle
        # AllostericMechanismSpec and produce TR-removal
        # candidates
        results = EnzymeRates.expand_mechanisms_by_one_param(
            [allo_spec], uni_uni_allosteric_R)

        tr_removed = filter(results) do r
            r isa EnzymeRates.AllostericMechanismSpec &&
            length(r.tr_equiv_metabolites) ==
                n_equiv - 1
        end

        # One candidate per TR-equivalent metabolite
        @test length(tr_removed) == n_equiv

        # Each TR-removed candidate should have one fewer
        # TR-equiv metabolite and higher param count
        base_pc = EnzymeRates._runtime_param_count(allo_spec)
        for r in tr_removed
            pc = EnzymeRates._runtime_param_count(r)
            @test pc == base_pc + 1
        end

        # Each removed metabolite should be from the original
        # TR-equiv set
        for r in tr_removed
            removed = setdiff(
                Set(allo_spec.tr_equiv_metabolites),
                Set(r.tr_equiv_metabolites))
            @test length(removed) == 1
            @test first(removed) in allo_spec.tr_equiv_metabolites
        end
    end

    @testset "Allosteric expand_by_one_param: base moves" begin
        # RE→SS on the base should produce an allosteric spec
        # with a modified base
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        spec = topos[1]

        allo_results = EnzymeRates.expand_mechanisms_by_two_params(
            [spec], uni_uni_allosteric_R)
        allo_spec = first(allo_results)

        results = EnzymeRates.expand_mechanisms_by_one_param(
            [allo_spec], uni_uni_allosteric_R)

        # Should have some results from base RE→SS moves
        base_re_count = count(
            s -> s.is_equilibrium, allo_spec.base.steps)

        # Results include TR-removal AND base moves
        @test length(results) >= 1

        # Verify param counts match compiled
        for r in first(results, 3)
            m = compile_mechanism(r)
            @test EnzymeRates._runtime_param_count(r) ==
                length(parameters(m))
        end
    end

    @testset "TR-equivalence covers SS rate constants" begin
        # When all RE binding metabolites are TR-equivalent,
        # SS rate constants should also be TR-equivalent
        # (k3f_T = k3f, so k3f_T is dependent, not independent)
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        spec = topos[1]

        # All metabolites TR-equivalent
        allo_all_equiv = EnzymeRates.AllostericMechanismSpec(
            spec, 1, [[:R]], [1], [:S, :P, :R])
        m_all = compile_mechanism(allo_all_equiv)
        params_all = parameters(m_all)

        # k3f_T should NOT appear (it should be dependent)
        @test :k3f_T ∉ params_all

        # Base has 3 indep (K1, K2, k3f) + Keq + E_total = 5
        # Allosteric adds: L + K_R_reg1 = +2
        # Total should be 7
        @test length(params_all) == 7

        # Runtime param count should match compiled
        runtime_pc = EnzymeRates._runtime_param_count(allo_all_equiv)
        @test runtime_pc == length(params_all)

        # With one metabolite NOT TR-equivalent (S):
        # S's binding K1_T becomes independent, AND
        # SS k3f_T becomes independent (since not all RE mets are equiv)
        # → delta from base = +4 (L, K_R_reg1, K1_T, k3f_T)
        allo_partial = EnzymeRates.AllostericMechanismSpec(
            spec, 1, [[:R]], [1], [:P, :R])
        m_partial = compile_mechanism(allo_partial)
        params_partial = parameters(m_partial)

        @test :k3f_T ∈ params_partial
        runtime_pc_partial = EnzymeRates._runtime_param_count(
            allo_partial)
        @test runtime_pc_partial == length(params_partial)
    end

    @testset "_runtime_param_count for AllostericMechanismSpec" begin
        # Verify runtime param count matches compiled for
        # allosteric specs from old pipeline
        all_specs = collect(
            EnzymeRates.old_enumerate_mechanisms(
                uni_uni_allosteric_R))
        allo_specs = filter(
            s -> s isa EnzymeRates.AllostericMechanismSpec,
            all_specs)
        @test !isempty(allo_specs)
        for s in first(allo_specs, 5)
            m = compile_mechanism(s)
            @test EnzymeRates._runtime_param_count(s) ==
                length(parameters(m))
        end
    end
end
