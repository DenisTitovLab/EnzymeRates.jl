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

end # top-level testset
