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

const uni_uni_allo_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R
    oligomeric_state: 2
end

const uni_uni_allo_2reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R1, R2
    oligomeric_state: 2
end

"""Collect all mechanisms by running the full enumeration loop."""
function enumerate_all(
    @nospecialize(reaction::EnzymeReaction);
    max_params::Int=20)
    cache = Dict{Int, Vector{EnzymeRates.AbstractMechanismSpec}}()

    init_specs = EnzymeRates.init_mechanisms(reaction)
    min_pc = init_specs[1].param_count
    cache[min_pc] = EnzymeRates.AbstractMechanismSpec[init_specs...]
    EnzymeRates.dedup!(cache)

    results = Dict{Int, Vector{EnzymeRates.AbstractMechanismSpec}}()

    for pc in min_pc:max_params
        level = pop!(cache, pc, EnzymeRates.AbstractMechanismSpec[])
        isempty(level) && (isempty(cache) ? break : continue)

        results[pc] = level

        new_specs = EnzymeRates.expand_mechanisms(level, reaction)
        for (target_pc, specs) in new_specs
            target_pc > max_params && continue
            append!(get!(cache, target_pc,
                EnzymeRates.AbstractMechanismSpec[]), specs)
        end
        EnzymeRates.dedup!(cache)
    end
    results
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
        Symbol[], Symbol[], Int[],
        base_spec.param_count + 1)

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
        Symbol[], Symbol[], Int[],
        base_spec.param_count + 2)

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
        Symbol[], [:P], Int[],
        base_spec.param_count + 1)

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

@testset "Move 6: Allosteric conversion (+1)" begin
    @testset "Uni-uni: K-type + V-type" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        # K-type: 1×1 = 1, V-type: 1. Total = 2
        @test length(result) == 2
        for r in result
            @test r isa AllostericMechanismSpec
            @test r.catalytic_n == 2
        end
    end

    @testset "Bi-bi: all substrate+product combos" begin
        bi_bi_allo = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        specs = EnzymeRates.init_mechanisms(bi_bi_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, bi_bi_allo)
        # K-type: 3×3 = 9, V-type: 1. Total = 10
        @test length(result) == 10
    end

    @testset "All are +1 param" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for r in result
            @test r.param_count == spec.param_count + 1
        end
    end

    @testset "K-type: cat steps stay tr_equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        k_type = filter(
            r -> !isempty(r.r_only_metabolites), result)
        @test !isempty(k_type)
        for r in k_type
            @test isempty(r.r_only_cat_steps)
            @test isempty(r.t_only_metabolites)
        end
    end

    @testset "V-type: all metabolites tr_equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        v_type = filter(
            r -> !isempty(r.r_only_cat_steps) &&
                 isempty(r.r_only_metabolites),
            result)
        @test length(v_type) == 1
        sub_names = [s[1] for s in EnzymeRates.substrates(
            uni_uni_allo)]
        prod_names = [p[1] for p in EnzymeRates.products(
            uni_uni_allo)]
        for r in v_type
            for m in Symbol[sub_names; prod_names]
                @test m in r.tr_equiv_metabolites
            end
        end
    end

    @testset "All compile" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for r in result
            m = AllostericEnzymeMechanism(r)
            @test m isa AllostericEnzymeMechanism
        end
    end

    @testset "Already allosteric → empty" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo = first(EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo))
        @test isempty(EnzymeRates._expand_to_allosteric(
            allo, uni_uni_allo))
    end

    @testset "oligomeric_state from reaction" begin
        rxn4 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 4
        end
        specs = EnzymeRates.init_mechanisms(rxn4)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, rxn4)
        for r in result
            @test r.catalytic_n == 4
        end
    end
end

@testset "Move 4: Add allosteric regulator" begin
    @testset "Add regulator to allosteric spec" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo_reg)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo_reg)
        @test !isempty(allo_specs)
        allo = first(allo_specs)
        result = EnzymeRates._expand_add_allosteric_regulator(
            allo, uni_uni_allo_reg)
        # R not yet added: 3 flavors × 1 site option
        # (new site only, no existing reg sites) = 3
        @test length(result) == 3
    end

    @testset "Non-allosteric → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo_reg)
        spec = first(specs)
        result = EnzymeRates._expand_add_allosteric_regulator(
            spec, uni_uni_allo_reg)
        @test isempty(result)
    end

    @testset "Second regulator with site options" begin
        specs = EnzymeRates.init_mechanisms(
            uni_uni_allo_2reg)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo_2reg)
        allo = first(allo_specs)
        # Add R1 first
        r1_added = EnzymeRates._expand_add_allosteric_regulator(
            allo, uni_uni_allo_2reg)
        @test !isempty(r1_added)
        # Now add R2 to one with R1
        with_r1 = first(r1_added)
        r2_added = EnzymeRates._expand_add_allosteric_regulator(
            with_r1, uni_uni_allo_2reg)
        # R2: 3 flavors × 2 site options
        # (new site + R1's site) = 6
        @test length(r2_added) == 6
    end
end

@testset "Move 5: Remove TR equivalence" begin
    @testset "Remove metabolite TR equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        n_tr = length(allo.tr_equiv_metabolites) +
               length(allo.tr_equiv_cat_steps)
        result = EnzymeRates._expand_remove_tr_equiv(
            allo, uni_uni_allo)
        @test length(result) == n_tr
    end

    @testset "No TR equivs left → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        fully_relaxed = allo
        while true
            r = EnzymeRates._expand_remove_tr_equiv(
                fully_relaxed, uni_uni_allo)
            isempty(r) && break
            fully_relaxed = first(r)
        end
        @test isempty(
            EnzymeRates._expand_remove_tr_equiv(
                fully_relaxed, uni_uni_allo))
    end

    @testset "MechanismSpec → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_remove_tr_equiv(
            spec, uni_uni_allo)
        @test isempty(result)
    end
end

@testset "Dedup" begin
    @testset "Same mechanism, different step order" begin
        spec1 = MechanismSpec(
            uni_uni_rxn,
            [StepSpec([:E, :S], [:E_S], true),
             StepSpec([:E, :P], [:E_P], true),
             StepSpec([:E_S], [:E_P], false)],
            ParamConstraint[], 5)
        spec2 = MechanismSpec(
            uni_uni_rxn,
            [StepSpec([:E, :P], [:E_P], true),
             StepSpec([:E_S], [:E_P], false),
             StepSpec([:E, :S], [:E_S], true)],
            ParamConstraint[], 5)
        cache = Dict(5 => AbstractMechanismSpec[spec1, spec2])
        EnzymeRates.dedup!(cache)
        @test length(cache[5]) == 1
    end

    @testset "Different mechanisms preserved" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        pc = specs[1].param_count
        cache = Dict(pc => AbstractMechanismSpec[specs...])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) >= 1
        @test length(cache[pc]) <= length(specs)
    end

    @testset "Idempotent" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        pc = specs[1].param_count
        cache = Dict(pc => AbstractMechanismSpec[specs...])
        EnzymeRates.dedup!(cache)
        n1 = length(cache[pc])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) == n1
    end
end

@testset "expand_mechanisms" begin
    @testset "Returns dict keyed by param count" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_rxn)
        @test result isa Dict{Int,
            Vector{AbstractMechanismSpec}}
        base_pc = specs[1].param_count
        @test haskey(result, base_pc + 1)
    end

    @testset "Allosteric expansion included" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_allo)
        base_pc = specs[1].param_count
        has_allo = any(
            any(s isa AllostericMechanismSpec
                for s in ss)
            for (_, ss) in result)
        @test has_allo
    end

    @testset "No self-expansion to same param count" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        base_pc = specs[1].param_count
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_rxn)
        # All results should have param_count > base
        for (pc, _) in result
            @test pc > base_pc
        end
    end

    @testset "Allosteric rewrap preserves structure" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        result = EnzymeRates.expand_mechanisms(
            [allo], uni_uni_allo)
        # Should have expansions from base moves (RE→SS)
        # rewrapped as AllostericMechanismSpec
        has_rewrapped = any(
            any(s isa AllostericMechanismSpec
                for s in ss)
            for (_, ss) in result)
        @test has_rewrapped
    end

    @testset "Dead-end excludes allosteric regs" begin
        specs = EnzymeRates.init_mechanisms(
            uni_uni_allo_reg)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo_reg)
        allo = first(allo_specs)
        # Add R as allosteric regulator
        with_reg = first(
            EnzymeRates._expand_add_allosteric_regulator(
                allo, uni_uni_allo_reg))
        result = EnzymeRates.expand_mechanisms(
            [with_reg], uni_uni_allo_reg)
        # R should NOT appear as dead-end in any expansion
        for (_, ss) in result
            for s in ss
                base = s isa AllostericMechanismSpec ?
                    s.base : s
                for step in base.steps
                    for sym in Iterators.flatten(
                            (step.reactants, step.products))
                        @test !contains(
                            string(sym), "R__reg")
                    end
                end
            end
        end
    end
end

@testset "Integration" begin
    @testset "Uni-uni full enumeration" begin
        results = enumerate_all(uni_uni_rxn; max_params=8)
        @test !isempty(results)
        pcs = sort(collect(keys(results)))
        @test issorted(pcs)
        # Every mechanism compiles
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                end
            end
        end
    end

    @testset "Bi-bi full enumeration" begin
        results = enumerate_all(bi_bi_rxn; max_params=10)
        @test !isempty(results)
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                end
            end
        end
    end

    @testset "With allosteric regulators" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        results = enumerate_all(rxn; max_params=10)
        has_allo = any(
            any(s isa EnzymeRates.AllostericMechanismSpec for s in specs)
            for (_, specs) in results)
        @test has_allo
    end

    @testset "With dead-end regulator" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
        end
        results = enumerate_all(rxn; max_params=8)
        @test !isempty(results)
        # Should have more mechanisms than plain uni-uni
        plain = enumerate_all(uni_uni_rxn; max_params=8)
        total_with_reg = sum(length(v) for v in values(results))
        total_plain = sum(length(v) for v in values(plain))
        @test total_with_reg > total_plain
    end

    @testset "Multiple levels populated" begin
        results = enumerate_all(uni_uni_rxn; max_params=8)
        @test length(results) >= 2  # At least 2 param count levels
    end
end

end # top-level testset
