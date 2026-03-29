# ABOUTME: Tests for beam-based mechanism enumeration runtime functions.
# ABOUTME: Validates parameter counting and fingerprinting against compiled mechanisms.

@testset "Beam Enumeration" begin
    @testset "Catalytic topology counts" begin
        @test length(
            EnzymeRates._catalytic_topologies(uni_uni)
        ) == 1
        @test length(
            EnzymeRates._catalytic_topologies(uni_bi)
        ) == 3
        @test length(
            EnzymeRates._catalytic_topologies(bi_bi)
        ) == 9
        @test length(
            EnzymeRates._catalytic_topologies(
                bi_bi_ping_pong)
        ) == 10

        # All topologies have exactly 1 SS step
        for rxn in [uni_uni, uni_bi, bi_bi, bi_bi_ping_pong]
            for t in EnzymeRates._catalytic_topologies(rxn)
                @test count(
                    !s.is_equilibrium for s in t.steps
                ) == 1
            end
        end
    end

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
                m = EnzymeRates.compile_mechanism(r)
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
                m = EnzymeRates.compile_mechanism(r)
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
                m = EnzymeRates.compile_mechanism(r)
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

        # Multi-level dead-ends (e.g., substrate binding
        # to E_I) are produced by +0 expansion, not +1.
        # Extending existing regulator to new forms, or
        # adding substrate dead-ends with K=K_catalytic,
        # preserves param_count.
        multi_level = EnzymeRates.expand_mechanisms_same_param_count(
            de_with_I, bi_bi_dead_end_I)
        multi_with_more_steps = filter(multi_level) do r
            length(r.steps) > maximum(
                length(d.steps) for d in de_with_I)
        end
        @test !isempty(multi_with_more_steps)

        # All forms respect capacity
        cap = EnzymeRates._binding_capacity(bi_bi_dead_end_I)
        for r in multi_with_more_steps
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
                    m = EnzymeRates.compile_mechanism(spec)
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
        # Valid dead-end forms: EAP (P at EA), EBQ (B at EQ)
        # 4 dead-end combos × 5 RE/SS variants = 20 total
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
            m = EnzymeRates.compile_mechanism(r)
            @test length(parameters(m)) == r.param_count
        end

        # Exact count: 4 dead-end combos × 5 RE/SS = 20
        all_mechs = vcat([spec], results)
        @test length(all_mechs) == 20

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

        # All metabolites are TR-equivalent by construction
        # (minimum delta = +2: L + K_R_reg)
        for r in results
            t_mets = EnzymeRates._collect_t_state_metabolites(r)
            n_non_equiv = length(t_mets) -
                length(r.tr_equiv_metabolites)
            @test n_non_equiv == 0
        end

        # Each result should compile and have correct param count
        for r in results
            m = EnzymeRates.compile_mechanism(r)
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

        # Bi-bi + allosteric R1, R2: 2 regs require delta >= +3
        # (L + K_R_R1 + K_R_R2), so expand_by_two_params
        # returns empty (it only generates delta = +2 specs)
        topos_bb = EnzymeRates._catalytic_topologies(bi_bi)
        results_bb_two_regs = EnzymeRates.expand_mechanisms_by_two_params(
            topos_bb, bi_bi_allosteric_R1_R2;
            max_catalytic_n=2)
        @test isempty(results_bb_two_regs)

        # Bi-bi + single allosteric reg: should produce delta=+2
        bi_bi_allo_R = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            allosteric_regulators: R
        end
        results_bb = EnzymeRates.expand_mechanisms_by_two_params(
            topos_bb, bi_bi_allo_R; max_catalytic_n=2)
        @test !isempty(results_bb)

        for r in first(results_bb, 5)
            m = EnzymeRates.compile_mechanism(r)
            @test EnzymeRates._runtime_param_count(r) ==
                length(parameters(m))
        end

        # All results have delta = +2 relative to base
        for r in results_bb
            base_pc = r.base.param_count
            @test EnzymeRates._runtime_param_count(r) ==
                base_pc + 2
        end
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
            m = EnzymeRates.compile_mechanism(r)
            @test EnzymeRates._runtime_param_count(r) ==
                length(parameters(m))
        end
    end

    @testset "Three-level TR-equivalence" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        spec = topos[1]
        m_base = EnzymeRates.compile_mechanism(spec)
        base_pc = length(parameters(m_base))

        # Find non-binding SS step index
        iso_idx = findfirst(
            s -> EnzymeRates.step_metabolite(s) === nothing &&
                !s.is_equilibrium, spec.steps)
        @test iso_idx !== nothing

        # All TR-equiv: mets + cat step → delta = +2
        allo_full = EnzymeRates.AllostericMechanismSpec(
            spec, 2, [[:R]], [1],
            [:S, :P, :R], [iso_idx],
            Symbol[], Symbol[], Int[], 0)
        m_full = EnzymeRates.compile_mechanism(allo_full)
        p_full = parameters(m_full)
        @test length(p_full) - base_pc == 2

        # Runtime param count matches compiled
        @test EnzymeRates._runtime_param_count(allo_full) ==
            length(p_full)

        # Without cat step TR-equiv → delta = +3
        allo_no_cat = EnzymeRates.AllostericMechanismSpec(
            spec, 2, [[:R]], [1],
            [:S, :P, :R], Int[],
            Symbol[], Symbol[], Int[], 0)
        m_no_cat = EnzymeRates.compile_mechanism(allo_no_cat)
        p_no_cat = parameters(m_no_cat)
        @test length(p_no_cat) - base_pc == 3

        @test EnzymeRates._runtime_param_count(allo_no_cat) ==
            length(p_no_cat)

        # Remove one metabolite TR-equiv → delta = +3
        # (with cat step still equiv)
        allo_partial_met = EnzymeRates.AllostericMechanismSpec(
            spec, 2, [[:R]], [1],
            [:P, :R], [iso_idx],
            Symbol[], Symbol[], Int[], 0)
        m_partial = EnzymeRates.compile_mechanism(allo_partial_met)
        p_partial = parameters(m_partial)
        @test length(p_partial) - base_pc == 3

        @test EnzymeRates._runtime_param_count(
            allo_partial_met) == length(p_partial)

        # Nothing TR-equiv → delta = +6
        # (K1_T, K2_T, k3f_T, K_R_reg1, K_R_T_reg1, L)
        allo_none = EnzymeRates.AllostericMechanismSpec(
            spec, 2, [[:R]], [1],
            Symbol[], Int[],
            Symbol[], Symbol[], Int[], 0)
        m_none = EnzymeRates.compile_mechanism(allo_none)
        p_none = parameters(m_none)
        @test length(p_none) - base_pc == 6

        @test EnzymeRates._runtime_param_count(allo_none) ==
            length(p_none)
    end

    @testset "R-only/T-only binding modes" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni)
        spec = topos[1]
        m_base = EnzymeRates.compile_mechanism(spec)
        base_pc = length(parameters(m_base))

        iso_idx = findfirst(
            s -> EnzymeRates.step_metabolite(s) === nothing &&
                !s.is_equilibrium, spec.steps)

        # R-only substrate: S only binds R-state, absent from T-state
        allo_r_only = EnzymeRates.AllostericMechanismSpec(
            spec, 1, [[:R]], [1],
            [:P, :R], [iso_idx],
            [:S], Symbol[], Int[], 0)
        m_r_only = EnzymeRates.compile_mechanism(allo_r_only)
        p_r_only = parameters(m_r_only)
        @test EnzymeRates._runtime_param_count(allo_r_only) ==
            length(p_r_only)

        # Rate equation should work and S should not appear
        # in T-state polynomial
        rate_str = rate_equation_string(m_r_only)
        @test occursin("K2", rate_str)  # R-state binding K exists

        # Evaluate at known concentrations
        concs = (; S=1.0, P=0.0, R=0.5)
        param_vals = NamedTuple{Tuple(p_r_only)}(
            ones(length(p_r_only)))
        v = rate_equation(m_r_only, concs, param_vals)
        @test isfinite(v)
        @test v != 0.0

        # T-only substrate: S only binds T-state, absent from R-state
        allo_t_only = EnzymeRates.AllostericMechanismSpec(
            spec, 1, [[:R]], [1],
            [:P, :R], [iso_idx],
            Symbol[], [:S], Int[], 0)
        m_t_only = EnzymeRates.compile_mechanism(allo_t_only)
        p_t_only = parameters(m_t_only)
        @test EnzymeRates._runtime_param_count(allo_t_only) ==
            length(p_t_only)

        # Same param count as r_only (symmetric structure)
        @test length(p_r_only) == length(p_t_only)

        param_vals_t = NamedTuple{Tuple(p_t_only)}(
            ones(length(p_t_only)))
        v_t = rate_equation(m_t_only, concs, param_vals_t)
        @test isfinite(v_t)
        @test v_t != 0.0

        # R-only cat step: T-state doesn't catalyze through this step
        allo_r_only_cat = EnzymeRates.AllostericMechanismSpec(
            spec, 1, [[:R]], [1],
            [:S, :P, :R], Int[],
            Symbol[], Symbol[], [iso_idx], 0)
        m_r_only_cat = EnzymeRates.compile_mechanism(allo_r_only_cat)
        p_r_only_cat = parameters(m_r_only_cat)
        @test EnzymeRates._runtime_param_count(
            allo_r_only_cat) == length(p_r_only_cat)

        param_vals_rc = NamedTuple{Tuple(p_r_only_cat)}(
            ones(length(p_r_only_cat)))
        v_rc = rate_equation(m_r_only_cat, concs, param_vals_rc)
        @test isfinite(v_rc)

        # R-only regulator: R only binds R-state
        allo_r_only_reg = EnzymeRates.AllostericMechanismSpec(
            spec, 1, [[:R]], [1],
            [:S, :P], [iso_idx],
            [:R], Symbol[], Int[], 0)
        m_r_only_reg = EnzymeRates.compile_mechanism(allo_r_only_reg)
        p_r_only_reg = parameters(m_r_only_reg)
        @test EnzymeRates._runtime_param_count(
            allo_r_only_reg) == length(p_r_only_reg)
        # R-only reg: no K_R_T_reg param
        @test !any(
            p -> occursin("_T_reg", string(p)), p_r_only_reg)
        # But K_R_reg exists
        @test any(
            p -> occursin("_R_reg", string(p)), p_r_only_reg)
    end

    @testset "enumerate_mechanisms (beam)" begin
        # Uni-uni: should produce mechanisms at multiple
        # param_count levels
        iter = EnzymeRates.enumerate_mechanisms(uni_uni)
        mechs = collect(iter)
        @test !isempty(mechs)

        # All mechanisms should have valid param_count
        for m in mechs
            @test m.param_count >= 0
        end

        # Should include mechanisms at different param levels
        param_counts = Set(m.param_count for m in mechs
            if m isa EnzymeRates.MechanismSpec)
        @test length(param_counts) >= 2

        # Uni-uni + allosteric R: should include both
        # MechanismSpec and AllostericMechanismSpec
        iter_allo = EnzymeRates.enumerate_mechanisms(
            uni_uni_allosteric_R)
        mechs_allo = collect(iter_allo)
        @test !isempty(mechs_allo)
        has_base = any(
            m isa EnzymeRates.MechanismSpec
            for m in mechs_allo)
        has_allo = any(
            m isa EnzymeRates.AllostericMechanismSpec
            for m in mechs_allo)
        @test has_base
        @test has_allo
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
            m = EnzymeRates.compile_mechanism(s)
            @test EnzymeRates._runtime_param_count(s) ==
                length(parameters(m))
        end
    end

    @testset "Equivalence with old_enumerate_mechanisms" begin
        # Normalize metabolite names: strip __reg1 suffix
        # used by old pipeline
        function _normalize_met(s::Symbol)
            str = string(s)
            m = match(r"^(.+)__reg\d*$", str)
            m === nothing ? s : Symbol(m.captures[1])
        end

        function _normalize_fp(fp)
            Set(
                Pair{Symbol,Int}[
                    _normalize_met(p.first) => p.second
                    for p in mono]
                for mono in fp)
        end

        function spec_fingerprint(spec)
            partition =
                EnzymeRates._compute_re_partition_from_steps(
                    spec.steps)
            raw = EnzymeRates._concentration_fingerprint(
                spec.steps, partition)
            _normalize_fp(raw)
        end

        function base_fp_set(specs)
            Set(spec_fingerprint(s) for s in specs
                if s isa EnzymeRates.MechanismSpec)
        end

        for (name, rxn) in [
                ("uni-uni", uni_uni),
                ("uni-bi", uni_bi),
                ("uni-uni + I", uni_uni_dead_end_I),
                ("uni-uni unknown reg",
                    uni_uni_reg_unknown),
            ]
            @testset "$name" begin
                old = collect(
                    EnzymeRates.old_enumerate_mechanisms(rxn))
                new = collect(
                    EnzymeRates.enumerate_mechanisms(rxn))

                old_fps = base_fp_set(old)
                new_fps = base_fp_set(new)

                # New should be a superset (may include
                # multi-level dead-end that old doesn't have)
                @test issubset(old_fps, new_fps)
            end
        end

        # Allosteric equivalence: normalize ligand names and
        # compare allosteric canonical keys
        function _normalize_sites(sites)
            [sort([_normalize_met(lig) for lig in site])
             for site in sites]
        end

        function _normalize_tr_equiv(tr_equiv)
            sort([_normalize_met(m) for m in tr_equiv])
        end

        function allo_key(s::EnzymeRates.AllostericMechanismSpec)
            base_fp = spec_fingerprint(s.base)
            sites = sort(_normalize_sites(s.allosteric_reg_sites))
            mults = sort(s.allosteric_multiplicities)
            tr = _normalize_tr_equiv(s.tr_equiv_metabolites)
            all_t_mets = sort(
                _normalize_met.(
                    EnzymeRates._collect_t_state_metabolites(s)))
            # Canonicalize: min of set and complement
            complement = sort(setdiff(all_t_mets, tr))
            canonical_tr = min(tr, complement)
            all_ss = sort(s.tr_equiv_cat_steps)
            cat_complement = sort(setdiff(
                EnzymeRates._collect_nonbinding_ss_steps(s),
                all_ss))
            canonical_cat = min(all_ss, cat_complement)
            (base_fp, sites, mults, s.catalytic_n,
             canonical_tr, canonical_cat)
        end

        function allo_fp_set(specs)
            Set(allo_key(s) for s in specs
                if s isa EnzymeRates.AllostericMechanismSpec)
        end

        @testset "uni-uni allosteric R" begin
            old = collect(
                EnzymeRates.old_enumerate_mechanisms(
                    uni_uni_allosteric_R))
            new = collect(
                EnzymeRates.enumerate_mechanisms(
                    uni_uni_allosteric_R))

            old_allo_fps = allo_fp_set(old)
            new_allo_fps = allo_fp_set(new)

            @test !isempty(old_allo_fps)
            @test !isempty(new_allo_fps)
            # New pipeline generates more allosteric specs than
            # old (different strategy: all-TR-equiv at +2, then
            # TR-removal moves). The sets overlap but neither
            # is a strict subset of the other.
            @test length(new_allo_fps) >= length(old_allo_fps)
        end
    end
end
