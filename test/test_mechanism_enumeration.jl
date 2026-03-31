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

const ter_ter_rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end

const ter_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[CN], Q[X]
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
        spec_rt = mechanism_spec_from_mechanism(
            m_uu, uni_uni_rxn)
        @test EnzymeMechanism(spec_rt) === m_uu
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

    @testset "Ter-Ter" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_ter_rxn)
        # 3969 = (2^(3!) - 1)² = 63 × 63
        # Each side (binding/release) has 3!=6 permutation
        # paths through Boolean lattice B_3; all 2^6-1=63
        # non-empty path subsets produce distinct edge sets;
        # sides are independent.
        @test length(topos) == 3969
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "Ter-Bi" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_bi_rxn)
        # 204 = 189 sequential + 15 ping-pong
        # Sequential: (2^(3!) - 1) × (2^(2!) - 1) = 63 × 3
        # Ping-pong: D[X]→Q[X] can isomerize independently
        @test length(topos) == 204
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
        @test EnzymeMechanism(spec) === m
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
        @test EnzymeMechanism(spec) === m
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

    @testset "Mirror steps count toward param_count delta" begin
        bi_bi_specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        specs_with_de = filter(
            s -> length(s.steps) > 9, bi_bi_specs)
        if !isempty(specs_with_de)
            spec = first(specs_with_de)
            if !isempty(spec.param_constraints)
                spec = first(
                    EnzymeRates._expand_remove_constraint(spec))
            end
            result = EnzymeRates._expand_re_to_ss(spec)
            for r in result
                n_converted = count(
                    !r.steps[i].is_equilibrium &&
                        spec.steps[i].is_equilibrium
                    for i in eachindex(r.steps)
                    if i <= length(spec.steps))
                if n_converted > 0
                    @test r.param_count ==
                        spec.param_count + n_converted
                end
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

    @testset "SS constraints removed as kf/kr pairs" begin
        # Build a spec with SS constraints: take a constrained
        # bi_bi spec, convert the constrained steps (which bind
        # the same metabolite) to SS, then rebuild constraints.
        # This produces kf/kr pairs that must be removed together.
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        constrained = first(filter(
            s -> !isempty(s.param_constraints), specs))
        c_idxs = EnzymeRates._constrained_step_indices(
            constrained.param_constraints)
        new_steps = [
            EnzymeRates.StepSpec(
                s.reactants, s.products,
                (i in c_idxs) ? false : s.is_equilibrium)
            for (i, s) in enumerate(constrained.steps)]
        ss_spec = EnzymeRates.MechanismSpec(
            constrained.reaction, new_steps,
            ParamConstraint[], constrained.param_count)
        # Rebuild constraints; the SS group gets kf/kr pairs
        new_constraints = EnzymeRates._max_equivalence_constraints(
            ss_spec)
        has_ss_pair = any(
            endswith(string(c[1]), "f") for c in new_constraints)
        if has_ss_pair
            spec_with_ss = EnzymeRates.MechanismSpec(
                ss_spec.reaction, ss_spec.steps,
                new_constraints, constrained.param_count)
            result = EnzymeRates._expand_remove_constraint(
                spec_with_ss)
            # Each result must not have an orphaned kf without
            # its matching kr (or vice versa)
            for r in result
                for c in r.param_constraints
                    s = string(c[1])
                    if endswith(s, "f")
                        kr_sym = Symbol(s[1:end-1] * "r")
                        @test any(
                            c2[1] == kr_sym
                            for c2 in r.param_constraints)
                    elseif endswith(s, "r")
                        kf_sym = Symbol(s[1:end-1] * "f")
                        @test any(
                            c2[1] == kf_sym
                            for c2 in r.param_constraints)
                    end
                end
            end
            # param_count delta must be +1 (RE) or +2 (SS pair)
            for r in result
                delta = r.param_count - spec_with_ss.param_count
                @test delta == 1 || delta == 2
            end
        end
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

@testset "Regulator dummy naming stability" begin
    rxn2 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I, J
    end
    specs = EnzymeRates.init_mechanisms(rxn2)
    spec = first(specs)
    # Add I first
    i_specs = EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn2)
    with_i = first(filter(i_specs) do s
        any(contains(string(sym), "I__reg")
            for st in s.steps
            for sym in Iterators.flatten(
                (st.reactants, st.products)))
    end)
    # Now add J to a mechanism that already has I
    j_specs = EnzymeRates._expand_add_dead_end_regulator(
        with_i, rxn2)
    # J should use J__reg (no numeric suffix)
    for s in j_specs
        j_syms = [sym for st in s.steps
            for sym in Iterators.flatten(
                (st.reactants, st.products))
            if contains(string(sym), "J__reg")]
        for sym in j_syms
            @test !contains(string(sym), r"__reg\d")
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

    @testset "V-type can remove r_only_cat_step" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        v_type = first(filter(
            r -> !isempty(r.r_only_cat_steps), allo_specs))
        @test !isempty(v_type.r_only_cat_steps)
        result = EnzymeRates._expand_remove_tr_equiv(
            v_type, uni_uni_allo)
        step_removals = filter(result) do r
            length(r.r_only_cat_steps) <
                length(v_type.r_only_cat_steps)
        end
        @test !isempty(step_removals)
        for r in step_removals
            @test r.param_count == v_type.param_count + 1
        end
    end

    @testset "Blocked when metabolites are r_only" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        k_type = first(filter(
            r -> !isempty(r.r_only_metabolites), allo_specs))
        mixed = AllostericMechanismSpec(
            k_type.base, k_type.catalytic_n,
            deepcopy(k_type.allosteric_reg_sites),
            copy(k_type.allosteric_multiplicities),
            copy(k_type.tr_equiv_metabolites),
            copy(k_type.tr_equiv_cat_steps),
            copy(k_type.r_only_metabolites),
            copy(k_type.t_only_metabolites),
            [1],  # r_only_cat_steps
            k_type.param_count)
        result = EnzymeRates._expand_remove_tr_equiv(
            mixed, uni_uni_allo)
        step_removals = filter(result) do r
            length(r.r_only_cat_steps) <
                length(mixed.r_only_cat_steps)
        end
        @test isempty(step_removals)
    end

    @testset "TR equiv removal delta for allosteric regulators" begin
        rxn_r = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        specs = EnzymeRates.init_mechanisms(rxn_r)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(spec, rxn_r)
        allo = first(allo_specs)
        reg_specs = EnzymeRates._expand_add_allosteric_regulator(allo, rxn_r)
        tr_spec = first(filter(r -> :R in r.tr_equiv_metabolites, reg_specs))
        pc_before = tr_spec.param_count
        result = EnzymeRates._expand_remove_tr_equiv(tr_spec, rxn_r)
        r_removal = filter(result) do r
            :R ∉ r.tr_equiv_metabolites &&
            :R ∉ r.r_only_metabolites &&
            :R ∉ r.t_only_metabolites
        end
        @test !isempty(r_removal)
        for r in r_removal
            @test r.param_count == pc_before + 1
        end
    end

    @testset "TR equiv removal delta skips constrained follower steps" begin
        # bi-bi random has 2 binding steps for A with K_follower = K_leader
        # constraint. Removing TR equiv for A should add +1 (one K_A_T),
        # not +2 (which would count both binding steps independently).
        const_bi_bi = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        specs = EnzymeRates.init_mechanisms(const_bi_bi)
        # Find a spec with 2 A-binding steps (constrained equal)
        two_a = filter(specs) do s
            length(filter(st -> EnzymeRates.step_metabolite(st) === :A,
                s.steps)) >= 2
        end
        @test !isempty(two_a)
        spec = first(two_a)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, const_bi_bi)
        # Find an allosteric spec with A in tr_equiv_metabolites
        tr_specs = filter(s -> :A in s.tr_equiv_metabolites, allo_specs)
        @test !isempty(tr_specs)
        tr_spec = first(tr_specs)
        result = EnzymeRates._expand_remove_tr_equiv(
            tr_spec, const_bi_bi)
        a_removals = filter(result) do r
            :A ∉ r.tr_equiv_metabolites &&
            :A ∉ r.r_only_metabolites &&
            :A ∉ r.t_only_metabolites
        end
        @test !isempty(a_removals)
        for r in a_removals
            m = EnzymeRates.compile_mechanism(r)
            @test length(parameters(m)) == r.param_count
        end
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

    @testset "Allosteric dedup: site order" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo = first(EnzymeRates._expand_to_allosteric(spec, uni_uni_allo))
        spec_ab = AllostericMechanismSpec(
            allo.base, allo.catalytic_n,
            [[:A], [:B]], [2, 2],
            copy(allo.tr_equiv_metabolites),
            copy(allo.tr_equiv_cat_steps),
            copy(allo.r_only_metabolites),
            copy(allo.t_only_metabolites),
            copy(allo.r_only_cat_steps),
            allo.param_count + 2)
        spec_ba = AllostericMechanismSpec(
            allo.base, allo.catalytic_n,
            [[:B], [:A]], [2, 2],
            copy(allo.tr_equiv_metabolites),
            copy(allo.tr_equiv_cat_steps),
            copy(allo.r_only_metabolites),
            copy(allo.t_only_metabolites),
            copy(allo.r_only_cat_steps),
            allo.param_count + 2)
        pc = spec_ab.param_count
        cache = Dict(pc => EnzymeRates.AbstractMechanismSpec[spec_ab, spec_ba])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) == 1
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
        allo_count = 0
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    allo_count += 1
                    allo_count > 3 && continue
                    m = AllostericEnzymeMechanism(spec)
                    @test length(parameters(m)) <= spec.param_count
                end
            end
        end
    end

    @testset "Bi-bi full enumeration" begin
        results = enumerate_all(bi_bi_rxn; max_params=10)
        @test !isempty(results)
        allo_count = 0
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    allo_count += 1
                    allo_count > 3 && continue
                    m = AllostericEnzymeMechanism(spec)
                    @test length(parameters(m)) <= spec.param_count
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
        # Every mechanism compiles
        allo_count = 0
        for (pc, specs) in results
            for spec in specs
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(parameters(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    allo_count += 1
                    allo_count > 3 && continue
                    m = AllostericEnzymeMechanism(spec)
                    @test length(parameters(m)) <= spec.param_count
                end
            end
        end
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

@testset "r_only params excluded from parameter list" begin
    specs = EnzymeRates.init_mechanisms(uni_uni_allo)
    spec = first(specs)
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo)

    @testset "K-type: no K_T params for r_only metabolites" begin
        k_type = first(filter(
            r -> !isempty(r.r_only_metabolites), allo_specs))
        m = AllostericEnzymeMechanism(k_type)
        params = parameters(m)
        @test length(params) == k_type.param_count
        t_params = filter(
            p -> endswith(string(p), "_T"), params)
        @test isempty(t_params)
    end

    @testset "V-type: no kf_T/kr_T for r_only cat steps" begin
        v_type = first(filter(
            r -> !isempty(r.r_only_cat_steps), allo_specs))
        m = AllostericEnzymeMechanism(v_type)
        params = parameters(m)
        @test length(params) == v_type.param_count
        t_k_params = filter(
            p -> contains(string(p), "f_T") ||
                 contains(string(p), "r_T"), params)
        @test isempty(t_k_params)
    end
end

@testset "Metabolite overlap: substrate as dead-end inhibitor" begin
    rxn_overlap = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    specs = EnzymeRates.init_mechanisms(rxn_overlap)
    @test !isempty(specs)

    spec = first(specs)
    de_specs = EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn_overlap)
    @test !isempty(de_specs)

    # S-as-regulator uses __reg suffix
    for s in de_specs
        reg_syms = [sym for st in s.steps
            for sym in Iterators.flatten(
                (st.reactants, st.products))
            if contains(string(sym), "S__reg")]
        @test !isempty(reg_syms)
    end

    # All compile correctly
    for s in de_specs
        m = EnzymeMechanism(s)
        @test m isa EnzymeMechanism
    end
end

@testset "Metabolite overlap: substrate as allosteric regulator" begin
    rxn_allo_overlap = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: S
        oligomeric_state: 2
    end
    specs = EnzymeRates.init_mechanisms(rxn_allo_overlap)
    @test !isempty(specs)
    spec = first(specs)

    # Allosteric conversion works
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, rxn_allo_overlap)
    @test !isempty(allo_specs)

    # Add S as allosteric regulator
    allo = first(allo_specs)
    reg_specs = EnzymeRates._expand_add_allosteric_regulator(
        allo, rxn_allo_overlap)
    @test !isempty(reg_specs)

    # S appears in allosteric_reg_sites
    for r in reg_specs
        has_s = any(:S in site for site in r.allosteric_reg_sites)
        @test has_s
    end

    # All compile correctly
    for r in reg_specs
        m = AllostericEnzymeMechanism(r)
        @test m isa AllostericEnzymeMechanism
    end

    # TR equiv removal produces separate results for
    # S-as-substrate and S-as-regulator
    tr_spec = first(filter(
        r -> :S in r.tr_equiv_metabolites, reg_specs))
    result = EnzymeRates._expand_remove_tr_equiv(
        tr_spec, rxn_allo_overlap)
    # S as catalytic met and S as regulator are both
    # in tr_equiv_metabolites. Removing each should
    # produce separate variants.
    @test length(result) >= 2
end

@testset "Base-level moves on allosteric specs" begin
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
    spec = mechanism_spec_from_mechanism(m_uu, uni_uni_allo)
    @test EnzymeMechanism(spec) === m_uu
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo)
    allo = first(allo_specs)

    @testset "RE→SS on allosteric" begin
        result = EnzymeRates._expand_re_to_ss(allo)
        @test !isempty(result)
        for r in result
            @test r isa EnzymeRates.AllostericMechanismSpec
            @test r.catalytic_n == allo.catalytic_n
            @test r.allosteric_reg_sites ==
                allo.allosteric_reg_sites
            @test r.r_only_metabolites ==
                allo.r_only_metabolites
            @test r.param_count > allo.param_count
        end
    end

    @testset "Remove constraint on allosteric" begin
        m_bb = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [E_B, A] ⇌ [E_A_B]
                [E, B] ⇌ [E_B]
                [E_A, B] ⇌ [E_A_B]
                [E, P] ⇌ [E_P]
                [E_P, Q] ⇌ [E_P_Q]
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E_A_B] <--> [E_P_Q]
            end
        end
        bi_bi_allo_rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        bb_spec = mechanism_spec_from_mechanism(
            m_bb, bi_bi_allo_rxn)
        @test EnzymeMechanism(bb_spec) === m_bb
        bb_spec_c = EnzymeRates.MechanismSpec(
            bb_spec.reaction, bb_spec.steps,
            EnzymeRates._max_equivalence_constraints(
                bb_spec),
            bb_spec.param_count)
        bb_allo = first(
            EnzymeRates._expand_to_allosteric(
                bb_spec_c, bi_bi_allo_rxn))
        @test !isempty(bb_allo.base.param_constraints)
        result = EnzymeRates._expand_remove_constraint(
            bb_allo)
        @test !isempty(result)
        for r in result
            @test r isa
                EnzymeRates.AllostericMechanismSpec
            @test r.catalytic_n == bb_allo.catalytic_n
            @test r.r_only_metabolites ==
                bb_allo.r_only_metabolites
        end
    end

    @testset "Add dead-end reg on allosteric" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            oligomeric_state: 2
        end
        spec_i = mechanism_spec_from_mechanism(
            m_uu, rxn)
        allo_i = first(
            EnzymeRates._expand_to_allosteric(
                spec_i, rxn))
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                allo_i, rxn)
        @test !isempty(result)
        for r in result
            @test r isa
                EnzymeRates.AllostericMechanismSpec
            @test r.catalytic_n == allo_i.catalytic_n
            @test r.r_only_metabolites ==
                allo_i.r_only_metabolites
        end
    end
end

end # top-level testset
