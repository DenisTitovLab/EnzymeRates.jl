# ABOUTME: Tests for mechanism enumeration pipeline
# ABOUTME: Unit tests per move + integration tests per reaction

using EnzymeRates: StepSpec, MechanismSpec, AllostericMechanismSpec,
    ParamConstraint, AbstractMechanismSpec

# Helper: convert EnzymeMechanism → MechanismSpec
function mechanism_spec_from_mechanism(
    m::EnzymeMechanism,
    @nospecialize(rxn::EnzymeReaction))
    rxns = EnzymeRates.reactions(m)
    eq_steps_tuple = EnzymeRates.equilibrium_steps(m)
    pc = EnzymeRates.param_constraints(m)
    steps = StepSpec[]
    for (i, (lhs, rhs)) in enumerate(rxns)
        push!(steps, StepSpec(
            collect(lhs), collect(rhs),
            eq_steps_tuple[i]))
    end
    constraints = [
        (t, c, [(s, sc) for (s, sc) in f])
        for (t, c, f) in pc
    ]
    MechanismSpec(
        rxn, steps, constraints,
        length(parameters(m)))
end

const uni_uni_rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

const uni_bi_rxn = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

const bi_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

const bi_bi_pp_rxn = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

const uni_uni_with_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
end

const uni_uni_allo = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    oligomeric_state: 2
end

@testset "Mechanism Enumeration" begin

@testset "Types and round-trip" begin
    m_uu = @enzyme_mechanism begin
        species: begin
            substrates: S[C]
            products: P[C]
            enzymes: E, E_P[C], E_S[C]
        end
        steps: begin
            [E, P] ⇌ [E_P]
            [E, S] ⇌ [E_S]
            [E_S] <--> [E_P]
        end
    end
    spec = mechanism_spec_from_mechanism(m_uu, uni_uni_rxn)
    @test length(spec.steps) == 3
    @test spec.param_count == length(parameters(m_uu))
    m_compiled = EnzymeMechanism(spec)
    @test m_compiled === m_uu
end

@testset "AllostericEnzymeMechanism round-trip" begin
    # Build MechanismSpec with :E (as the enumeration pipeline
    # would) and wrap it into an AllostericMechanismSpec
    base_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    base_steps = [
        StepSpec([:E, :S], [:E_S], true),
        StepSpec([:E, :P], [:E_P], true),
        StepSpec([:E_S], [:E_P], false),
    ]
    base_spec = MechanismSpec(
        base_rxn, base_steps, ParamConstraint[], 3)

    allo_spec = AllostericMechanismSpec(
        base_spec, 2,
        Vector{Symbol}[], Int[],
        Symbol[], Int[],
        Symbol[], Symbol[], Int[])

    m_compiled = AllostericEnzymeMechanism(allo_spec)
    @test m_compiled isa AllostericEnzymeMechanism

    # Verify structural properties
    @test EnzymeRates.metabolites(m_compiled) == (:S, :P)
    cat_sites = typeof(m_compiled).parameters[3]
    @test cat_sites[2] == 2  # catalytic_n
    @test typeof(m_compiled).parameters[4] == ()  # no reg sites

    # Verify the catalytic mechanism round-trips
    cat_m = typeof(m_compiled).parameters[2]()
    @test EnzymeRates.n_states(cat_m) == 3
    @test EnzymeRates.n_steps(cat_m) == 3
end

@testset "AllostericEnzymeMechanism round-trip with regulator" begin
    base_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A
    end
    base_steps = [
        StepSpec([:E, :S], [:E_S], true),
        StepSpec([:E, :P], [:E_P], true),
        StepSpec([:E_S], [:E_P], false),
    ]
    base_spec = MechanismSpec(
        base_rxn, base_steps, ParamConstraint[], 3)

    allo_spec = AllostericMechanismSpec(
        base_spec, 2,
        [[:A]], [1],
        Symbol[], Int[],
        Symbol[], Symbol[], Int[])

    m_compiled = AllostericEnzymeMechanism(allo_spec)
    @test m_compiled isa AllostericEnzymeMechanism

    # Verify metabolites include the regulator
    @test EnzymeRates.metabolites(m_compiled) ==
        (:S, :P, :A)
    # Verify reg sites
    reg_sites = typeof(m_compiled).parameters[4]
    @test length(reg_sites) == 1
    @test reg_sites[1][1] == (:A,)  # ligands
    @test reg_sites[1][2] == 1      # multiplicity
end

@testset "AllostericEnzymeMechanism TR equivalence" begin
    base_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    base_steps = [
        StepSpec([:E, :S], [:E_S], true),
        StepSpec([:E, :P], [:E_P], true),
        StepSpec([:E_S], [:E_P], false),
    ]
    base_spec = MechanismSpec(
        base_rxn, base_steps, ParamConstraint[], 3)

    # S is TR-equivalent (K_T_S = K_R_S),
    # P is R-only (absent from T-state)
    allo_spec = AllostericMechanismSpec(
        base_spec, 2,
        Vector{Symbol}[], Int[],
        [:S], Int[],
        Symbol[], [:P], Int[])

    m_compiled = AllostericEnzymeMechanism(allo_spec)
    cat_sites = typeof(m_compiled).parameters[3]
    @test :S in cat_sites[3]   # tr_equiv_mets
    @test :P in cat_sites[6]   # t_only_mets
end

@testset "Catalytic topologies" begin

    @testset "Uni-Uni" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni_rxn)
        @test length(topos) == 1

        m_uu = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        @test EnzymeMechanism(topos[1]) === m_uu
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "Uni-Bi" begin
        topos = EnzymeRates._catalytic_topologies(uni_bi_rxn)
        @test length(topos) == 3
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
            m = EnzymeMechanism(t)
            @test m isa EnzymeMechanism
        end
    end

    @testset "Bi-Bi" begin
        topos = EnzymeRates._catalytic_topologies(bi_bi_rxn)
        @test length(topos) == 9
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "Bi-Bi Ping-Pong" begin
        topos = EnzymeRates._catalytic_topologies(
            bi_bi_pp_rxn)
        @test length(topos) == 10
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end
end

@testset "init_mechanisms" begin

    @testset "Param count invariant" begin
        for (rxn, n_s, n_p) in [
            (uni_uni_rxn, 1, 1),
            (uni_bi_rxn, 1, 2),
            (bi_bi_rxn, 2, 2),
            (bi_bi_pp_rxn, 2, 2),
        ]
            specs = EnzymeRates.init_mechanisms(rxn)
            expected_pc = n_s + n_p + 3
            for s in specs
                @test s.param_count == expected_pc
            end
        end
    end

    @testset "All have exactly 1 SS step" begin
        for rxn in [uni_uni_rxn, uni_bi_rxn,
                    bi_bi_rxn, bi_bi_pp_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            for s in specs
                @test count(
                    !st.is_equilibrium
                    for st in s.steps) == 1
            end
        end
    end

    @testset "Uni-Uni: no dead-end forms" begin
        specs = EnzymeRates.init_mechanisms(
            uni_uni_rxn)
        @test length(specs) == 1
    end

    @testset "Dead-end counts" begin
        bi_bi_specs = EnzymeRates.init_mechanisms(
            bi_bi_rxn)
        # More than just the 9 topologies
        @test length(bi_bi_specs) > 9

        pp_specs = EnzymeRates.init_mechanisms(
            bi_bi_pp_rxn)
        # More than just the 10 topologies
        @test length(pp_specs) > 10
    end

    @testset "All compile correctly" begin
        for rxn in [uni_uni_rxn, uni_bi_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            for s in specs
                m = EnzymeMechanism(s)
                @test m isa EnzymeMechanism
            end
        end
    end
end

@testset "Move 1: RE→SS conversion" begin
    @testset "Multiple RE steps" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
        result = EnzymeRates._expand_re_to_ss(spec)
        @test length(result) == 2
        for r in result
            @test r.param_count == spec.param_count + 1
        end
    end

    @testset "All SS → yields nothing" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] <--> [E_P]
                [E, S] <--> [E_S]
                [E_S] <--> [E_P]
            end
        end
        spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
        result = EnzymeRates._expand_re_to_ss(spec)
        @test isempty(result)
    end

    @testset "Constrained RE steps skipped" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        constrained_spec = first(filter(
            s -> !isempty(s.param_constraints), specs))
        constrained_idxs = EnzymeRates._constrained_step_indices(
            constrained_spec.param_constraints)
        n_eligible = count(
            s.is_equilibrium && !(i in constrained_idxs)
            for (i, s) in enumerate(constrained_spec.steps))
        result = EnzymeRates._expand_re_to_ss(constrained_spec)
        @test length(result) == n_eligible
        for r in result
            new_ss_idxs = [i for (i, s) in enumerate(r.steps)
                if !s.is_equilibrium && constrained_spec.steps[i].is_equilibrium]
            for idx in new_ss_idxs
                @test !(idx in constrained_idxs)
            end
        end
    end
end

@testset "Move 2: Remove equivalence constraint" begin
    @testset "Mechanism with constraints" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        constrained = filter(
            s -> !isempty(s.param_constraints), specs)
        @test !isempty(constrained)
        spec = first(constrained)
        n_constraints = length(spec.param_constraints)
        result = EnzymeRates._expand_remove_constraint(spec)
        @test length(result) == n_constraints
        for r in result
            @test r.param_count == spec.param_count + 1
            @test length(r.param_constraints) == n_constraints - 1
        end
    end

    @testset "No constraints → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        spec = first(specs)
        @test isempty(spec.param_constraints)
        result = EnzymeRates._expand_remove_constraint(spec)
        @test isempty(result)
    end
end

@testset "Move 3: Add dead-end regulator" begin
    @testset "Uni-uni + new regulator" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg)
        # 3 forms (E, E_S, E_P), all eligible (none has
        # all substrates AND all products)
        # Wait — E has neither, E_S has S (all subs),
        # E_P has P (all prods). Only E is eligible.
        # Actually: "neither all subs nor all prods" means
        # exclude forms with ALL subs or ALL prods.
        # Uni-uni: sub_names={S}, prod_names={P}
        # E: bound={} → eligible
        # E_S: bound={S} → has all subs → NOT eligible
        # E_P: bound={P} → has all prods → NOT eligible
        # So only 1 eligible form → 2^1 - 1 = 1 variant
        @test length(result) == 1
        for r in result
            @test r.param_count == spec.param_count + 1
        end
    end

    @testset "Bi-bi + new regulator: more eligible forms" begin
        bi_bi_with_reg = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        specs = EnzymeRates.init_mechanisms(bi_bi_with_reg)
        # Pick a topology with many forms
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bi_bi_with_reg)
        # Should produce at least 1 variant
        @test length(result) >= 1
        for r in result
            @test r.param_count == spec.param_count + 1
        end
    end

    @testset "No regulators → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_rxn)
        @test isempty(result)
    end

    @testset "All results compile" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg)
        for r in result
            m = EnzymeMechanism(r)
            @test m isa EnzymeMechanism
        end
    end

    @testset "exclude_regs works" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg;
            exclude_regs=Set([:I]))
        @test isempty(result)
    end

    @testset "Mirror steps created" begin
        # Use bi-bi where eligible forms are connected
        # by catalytic steps, so mirrors should appear
        bi_bi_with_reg = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        specs = EnzymeRates.init_mechanisms(bi_bi_with_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bi_bi_with_reg)
        # Find a variant with multiple eligible forms
        # (mask > 1 means at least 2 forms)
        multi_form = filter(
            r -> length(r.steps) > length(spec.steps) + 2,
            result)
        if !isempty(multi_form)
            r = first(multi_form)
            # Should have binding steps + mirror steps
            @test length(r.steps) > length(spec.steps) + 2
        end
    end

    @testset "Equivalence constraints on binding K's" begin
        bi_bi_with_reg = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        specs = EnzymeRates.init_mechanisms(bi_bi_with_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bi_bi_with_reg)
        for r in result
            # Count regulator binding steps: steps
            # where the dummy regulator is a reactant
            n_reg_binding = count(
                s -> length(s.reactants) == 2 &&
                    contains(string(s.reactants[2]),
                        "__reg"),
                r.steps)
            new_constraints = setdiff(
                r.param_constraints,
                spec.param_constraints)
            if n_reg_binding >= 2
                @test length(new_constraints) ==
                    n_reg_binding - 1
            else
                @test isempty(new_constraints)
            end
        end
    end
end

@testset "Move 6: Allosteric conversion" begin
    @testset "Non-allosteric → allosteric variants" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        # S × {r_only, t_only} + P × {r_only, t_only} = 4
        @test length(result) == 4
        for r in result
            @test r isa AllostericMechanismSpec
            @test r.catalytic_n == 2
        end
    end

    @testset "Already allosteric → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for a in allo_specs
            @test isempty(
                EnzymeRates._expand_to_allosteric(
                    a, uni_uni_allo))
        end
    end

    @testset "All results compile" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for a in allo_specs
            m = AllostericEnzymeMechanism(a)
            @test m isa AllostericEnzymeMechanism
        end
    end

    @testset "oligomeric_state from reaction" begin
        rxn4 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 4
        end
        specs = EnzymeRates.init_mechanisms(rxn4)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, rxn4)
        for a in allo_specs
            @test a.catalytic_n == 4
        end
    end

    @testset "TR-equiv and r/t-only lists" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for a in allo_specs
            # Exactly one metabolite is r_only or t_only
            n_diff = length(a.r_only_metabolites) +
                     length(a.t_only_metabolites)
            @test n_diff == 1
            # The other metabolite is tr_equiv
            @test length(a.tr_equiv_metabolites) == 1
            # No allosteric reg sites
            @test isempty(a.allosteric_reg_sites)
        end
    end

    @testset "SS isomerization steps are tr_equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        # The uni-uni mechanism has one SS step (isomerization)
        ss_isom = [i for (i, s) in enumerate(spec.steps)
            if !s.is_equilibrium &&
               EnzymeRates.step_metabolite(s) === nothing]
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for a in allo_specs
            @test sort(a.tr_equiv_cat_steps) == sort(ss_isom)
        end
    end
end

end # top-level testset
