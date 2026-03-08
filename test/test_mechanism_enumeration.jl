# Tests for the staged mechanism enumeration pipeline
# Tests each stage independently (StageExpansionTestSpec),
# end-to-end (EnumerationTestSpec), and property-based tests.

using Random

# ── Test spec structs ──────────────────────────────────────────

"""
Tests a single hand-built MechanismSpec through each pipeline
stage, comparing output counts at each stage.
"""
Base.@kwdef struct StageExpansionTestSpec
    name::String
    reaction::Any
    base_mechanism::EnzymeRates.MechanismSpec
    dead_end_regs::Vector{Symbol} = Symbol[]
    allosteric_regs::Vector{Symbol} = Symbol[]
    catalytic_n::Int = 0

    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int = 0
    expected_n_tr_equiv::Int = 0
    expected_n_oem_dedup::Int = 0
end

"""
Tests end-to-end from EnzymeReaction through full pipeline,
comparing output count at each stage.
"""
Base.@kwdef struct EnumerationTestSpec
    name::String
    reaction::Any
    catalytic_n::Int = 0

    expected_n_forms::Int
    expected_n_catalytic::Int
    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int = 0
    expected_n_tr_equiv::Int = 0
    expected_n_oem_dedup::Int = 0
    expected_n_total::Int
end

# ── Helper: run pipeline stage by stage (full, all partitions) ─

function _run_full_pipeline_stages(rxn; catalytic_n::Int=0,
                                   max_re_groups::Int=7)
    catalytic = EnzymeRates._catalytic_topologies(rxn)

    roles = EnzymeRates.regulator_roles(rxn)
    fixed_de = Symbol[r[1] for r in roles if r[2] == :dead_end]
    fixed_allo = Symbol[r[1] for r in roles
                        if r[2] == :allosteric]
    unknown = Symbol[r[1] for r in roles if r[2] == :unknown]
    n_unknown = length(unknown)

    n_ress = 0; n_gm = 0; n_ea = 0
    n_de = 0; n_eq = 0; n_dd = 0
    n_allo = 0; n_tr = 0; n_oem_dd = 0

    for reg_mask in 0:(1 << n_unknown) - 1
        de = Symbol[fixed_de;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 0]]
        al = Symbol[fixed_allo;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 1]]

        ress = EnzymeRates._expand_ress_variants(
            catalytic, rxn; max_re_groups)
        n_ress += length(ress)

        gm = EnzymeRates._expand_general_modifiers(
            ress, rxn; allosteric_regs=al)
        n_gm += length(gm)

        ea = EnzymeRates._expand_essential_activators(
            gm, rxn; allosteric_regs=al)
        n_ea += length(ea)

        with_de = EnzymeRates._expand_dead_end_inhibitors(
            ea, rxn; dead_end_regs=de)
        n_de += length(with_de)

        with_eq = EnzymeRates._expand_equivalence_constraints(
            with_de, rxn)
        n_eq += length(with_eq)

        dd = EnzymeRates._deduplicate(with_eq, rxn)
        n_dd += length(dd)

        if catalytic_n > 0
            oem = EnzymeRates._expand_allosteric(
                dd, rxn; catalytic_n, allosteric_regs=al)
            n_allo += length(oem)
            oem = EnzymeRates._expand_tr_equivalence(oem, rxn)
            n_tr += length(oem)
            oem = EnzymeRates._deduplicate_oem(oem, rxn)
            n_oem_dd += length(oem)
        end
    end

    (catalytic=length(catalytic), ress=n_ress,
     general_modifier=n_gm, essential_activator=n_ea,
     dead_end=n_de, equivalence=n_eq, dedup=n_dd,
     allosteric=n_allo, tr_equiv=n_tr, oem_dedup=n_oem_dd)
end

# ── Test reactions ─────────────────────────────────────────────

const _rxn_uu = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

const _rxn_uu_reg = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    regulators: I
end

const _rxn_ub = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

const _rxn_ub_reg = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    regulators: I
end

const _rxn_bu = @enzyme_reaction begin
    substrates: A[X], B[Y]
    products: P[XY]
end

const _rxn_ub_allo = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: I
end

const _rxn_bb = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

# ── Build hand-built MechanismSpecs ──────────────────────────

# Uni-Uni forms:
#   1: E_0_0   (free enzyme)
#   2: E_S_0   (substrate bound)
#   3: E_0_P   (product bound)
# Adjacency: (1,2)=S, (1,3)=P, (2,3)=isomerization
# Triangle: 3 edges, 3 forms, 1 cycle, 1 thermo constraint
# All RE except isomerization (2,3) which is SS:
#   n_RE=2, n_SS=1, param_count = 2 + 2*1 - 1 + 2 = 5
const _hand_uu = EnzymeRates.MechanismSpec(
    _rxn_uu,
    [(1, 2), (1, 3), (2, 3)],  # S-binding, P-binding, isom
    3,                          # n_catalytic_edges
    [true, true, false],        # RE, RE, SS
    EnzymeRates.ParamConstraint[],
    5,                          # param_count
)

# Bi-Uni forms:
#   1: E_0_0_0  (free enzyme)
#   2: E_A_0_0  (A bound)
#   3: E_0_B_0  (B bound)
#   4: E_A_B_0  (A+B bound)
#   5: E_0_0_P  (product bound)
# Adjacency: (1,2)=A, (1,3)=B, (1,5)=P, (2,4)=B, (3,4)=A,
#            (4,5)=isomerization
# One ordered path: E→EA→EAB→EP→E (4 edges, 4 forms, 1 cycle)
# SS = isomerization (4,5), rest RE
#   n_RE=3, n_SS=1, param_count = 3 + 2*1 - 1 + 2 = 6
const _hand_bu = EnzymeRates.MechanismSpec(
    _rxn_bu,
    [(2, 4), (1, 2), (4, 5), (1, 5)],  # B-bind, A-bind, isom, P-rel
    4,                                   # n_catalytic_edges
    [true, true, false, true],           # RE, RE, SS, RE
    EnzymeRates.ParamConstraint[],
    6,                                   # param_count
)

function _first_catalytic(rxn)
    EnzymeRates._catalytic_topologies(rxn)[1]
end

# ── StageExpansionTestSpec instances ───────────────────────────

const STAGE_EXPANSION_SPECS = [
    # Uni-Uni: hand-built 3-edge triangle (E, ES, EP), no regs
    StageExpansionTestSpec(
        name="Uni-Uni hand-built, no regs",
        reaction=_rxn_uu,
        base_mechanism=_hand_uu,
        expected_n_ress=3,
        expected_n_general_modifier=3,
        expected_n_essential_activator=3,
        expected_n_dead_end=3,
        expected_n_equivalence=3,
        expected_n_dedup=1,
    ),
    # Bi-Uni: hand-built ordered path, no regs
    StageExpansionTestSpec(
        name="Bi-Uni hand-built, no regs",
        reaction=_rxn_bu,
        base_mechanism=_hand_bu,
        expected_n_ress=7,
        expected_n_general_modifier=7,
        expected_n_essential_activator=7,
        expected_n_dead_end=7,
        expected_n_equivalence=7,
        expected_n_dedup=2,
    ),
    # Uni-Bi first topology, dead-end reg I
    StageExpansionTestSpec(
        name="Uni-Bi first topology, dead-end I",
        reaction=_rxn_ub_reg,
        base_mechanism=_first_catalytic(_rxn_ub_reg),
        dead_end_regs=[:I],
        expected_n_ress=7,
        expected_n_general_modifier=7,
        expected_n_essential_activator=7,
        expected_n_dead_end=112,
        expected_n_equivalence=140,
        expected_n_dedup=74,
    ),
    # Uni-Bi first topology, allosteric reg I
    StageExpansionTestSpec(
        name="Uni-Bi first topology, allosteric I",
        reaction=_rxn_ub_reg,
        base_mechanism=_first_catalytic(_rxn_ub_reg),
        allosteric_regs=[:I],
        expected_n_ress=7,
        expected_n_general_modifier=14,
        expected_n_essential_activator=21,
        expected_n_dead_end=21,
        expected_n_equivalence=42,
        expected_n_dedup=28,
    ),
]

# ── EnumerationTestSpec instances ──────────────────────────────

const ENUMERATION_SPECS = [
    EnumerationTestSpec(
        name="Uni-Uni, no regulators",
        reaction=_rxn_uu,
        expected_n_forms=3,
        expected_n_catalytic=1,
        expected_n_ress=3,
        expected_n_general_modifier=3,
        expected_n_essential_activator=3,
        expected_n_dead_end=3,
        expected_n_equivalence=3,
        expected_n_dedup=1,
        expected_n_total=1,
    ),
    EnumerationTestSpec(
        name="Uni-Uni + 1 unconstrained reg",
        reaction=_rxn_uu_reg,
        expected_n_forms=6,
        expected_n_catalytic=1,
        # 2 partitions (de/allo), each 3 RE/SS = 6 total
        expected_n_ress=6,
        expected_n_general_modifier=9,
        expected_n_essential_activator=12,
        expected_n_dead_end=33,
        expected_n_equivalence=48,
        expected_n_dedup=28,
        expected_n_total=28,
    ),
    EnumerationTestSpec(
        name="Uni-Bi, no regulators",
        reaction=_rxn_ub,
        expected_n_forms=11,
        expected_n_catalytic=3,
        expected_n_ress=41,
        expected_n_general_modifier=41,
        expected_n_essential_activator=41,
        expected_n_dead_end=41,
        expected_n_equivalence=41,
        expected_n_dedup=9,
        expected_n_total=9,
    ),
    EnumerationTestSpec(
        name="Uni-Bi + 1 reg",
        reaction=_rxn_ub_reg,
        expected_n_forms=22,
        expected_n_catalytic=3,
        # 2 partitions: de gets 41, allo gets 41 = 82
        expected_n_ress=82,
        expected_n_general_modifier=123,
        expected_n_essential_activator=164,
        expected_n_dead_end=1211,
        expected_n_equivalence=1606,
        expected_n_dedup=697,
        expected_n_total=697,
    ),
    EnumerationTestSpec(
        name="Bi-Uni, no regulators",
        reaction=_rxn_bu,
        expected_n_forms=5,
        expected_n_catalytic=3,
        expected_n_ress=41,
        expected_n_general_modifier=41,
        expected_n_essential_activator=41,
        expected_n_dead_end=41,
        expected_n_equivalence=74,
        expected_n_dedup=23,
        expected_n_total=23,
    ),
    EnumerationTestSpec(
        name="Uni-Bi + allosteric, catalytic_n=2",
        reaction=_rxn_ub_allo,
        catalytic_n=2,
        expected_n_forms=22,
        expected_n_catalytic=3,
        expected_n_ress=41,
        expected_n_general_modifier=82,
        expected_n_essential_activator=123,
        expected_n_dead_end=123,
        expected_n_equivalence=246,
        expected_n_dedup=126,
        expected_n_allosteric=252,
        expected_n_tr_equiv=252,
        expected_n_oem_dedup=252,
        expected_n_total=378,
    ),
    EnumerationTestSpec(
        name="Bi-Bi, no regulators",
        reaction=_rxn_bb,
        expected_n_forms=11,
        expected_n_catalytic=9,
        expected_n_ress=527,
        expected_n_general_modifier=527,
        expected_n_essential_activator=527,
        expected_n_dead_end=527,
        expected_n_equivalence=958,
        expected_n_dedup=207,
        expected_n_total=207,
    ),
]

# ── Tests ──────────────────────────────────────────────────────

@testset "Mechanism Enumeration Pipeline" begin

    # ── StageExpansionTestSpec tests ──────────────────────────

    @testset "Stage expansion: $(s.name)" for s in STAGE_EXPANSION_SPECS
        spec = s.base_mechanism
        rxn = s.reaction

        ress = EnzymeRates._expand_ress_variants([spec], rxn)
        @test length(ress) == s.expected_n_ress

        gm = EnzymeRates._expand_general_modifiers(
            ress, rxn; allosteric_regs=s.allosteric_regs)
        @test length(gm) == s.expected_n_general_modifier

        ea = EnzymeRates._expand_essential_activators(
            gm, rxn; allosteric_regs=s.allosteric_regs)
        @test length(ea) == s.expected_n_essential_activator

        de = EnzymeRates._expand_dead_end_inhibitors(
            ea, rxn; dead_end_regs=s.dead_end_regs)
        @test length(de) == s.expected_n_dead_end

        eq = EnzymeRates._expand_equivalence_constraints(
            de, rxn)
        @test length(eq) == s.expected_n_equivalence

        dd = EnzymeRates._deduplicate(eq, rxn)
        @test length(dd) == s.expected_n_dedup

        if s.catalytic_n > 0
            oem = EnzymeRates._expand_allosteric(
                dd, rxn; catalytic_n=s.catalytic_n,
                allosteric_regs=s.allosteric_regs)
            @test length(oem) == s.expected_n_allosteric

            oem = EnzymeRates._expand_tr_equivalence(
                oem, rxn)
            @test length(oem) == s.expected_n_tr_equiv

            oem = EnzymeRates._deduplicate_oem(oem, rxn)
            @test length(oem) == s.expected_n_oem_dedup
        end
    end

    # ── EnumerationTestSpec tests ────────────────────────────

    @testset "End-to-end: $(s.name)" for s in ENUMERATION_SPECS
        rxn = s.reaction

        _, forms = EnzymeRates.enumerate_enzyme_forms(rxn)
        @test length(forms) == s.expected_n_forms

        counts = _run_full_pipeline_stages(
            rxn; catalytic_n=s.catalytic_n)
        @test counts.catalytic == s.expected_n_catalytic
        @test counts.ress == s.expected_n_ress
        @test counts.general_modifier ==
            s.expected_n_general_modifier
        @test counts.essential_activator ==
            s.expected_n_essential_activator
        @test counts.dead_end == s.expected_n_dead_end
        @test counts.equivalence == s.expected_n_equivalence
        @test counts.dedup == s.expected_n_dedup

        if s.catalytic_n > 0
            @test counts.allosteric == s.expected_n_allosteric
            @test counts.tr_equiv == s.expected_n_tr_equiv
            @test counts.oem_dedup == s.expected_n_oem_dedup
        end

        total = collect(EnzymeRates.enumerate_mechanisms(
            rxn; catalytic_n=s.catalytic_n))
        @test length(total) == s.expected_n_total
    end

    # ── Property-based tests ─────────────────────────────────

    @testset "Catalytic topology properties" begin
        rxn = _rxn_ub
        catalytic = EnzymeRates._catalytic_topologies(rxn)
        @test length(catalytic) == 3
        for spec in catalytic
            @test spec.n_catalytic_edges == length(spec.edges)
            # Exactly one SS step (first isomerization)
            @test count(.!spec.equilibrium_steps) == 1
        end
    end

    @testset "RE/SS expansion properties" begin
        rxn = _rxn_ub
        catalytic = EnzymeRates._catalytic_topologies(rxn)
        spec = catalytic[1]
        ress = EnzymeRates._expand_ress_variants(
            [spec], rxn)
        @test length(ress) > 0
        for s in ress
            @test any(.!s.equilibrium_steps)
            partition = EnzymeRates._compute_re_partition(
                s.edges, s.equilibrium_steps)
            @test 2 <= length(partition) <= 7
        end
    end

    @testset "Dead-end expansion properties" begin
        rxn = _rxn_ub_reg
        catalytic = EnzymeRates._catalytic_topologies(rxn)
        de_specs = EnzymeRates._expand_dead_end_inhibitors(
            catalytic, rxn; dead_end_regs=[:I])
        @test length(de_specs) > length(catalytic)
        for s in de_specs
            @test length(s.edges) >= s.n_catalytic_edges
        end

        # No dead-end regs: passthrough
        no_de = EnzymeRates._expand_dead_end_inhibitors(
            catalytic, rxn; dead_end_regs=Symbol[])
        @test length(no_de) == length(catalytic)
    end

    @testset "Deduplication reduces count" begin
        rxn = _rxn_ub_reg
        catalytic = EnzymeRates._catalytic_topologies(rxn)
        de_specs = EnzymeRates._expand_dead_end_inhibitors(
            catalytic, rxn; dead_end_regs=[:I])
        ress = EnzymeRates._expand_ress_variants(
            [de_specs[1]], rxn)
        with_eq = EnzymeRates._expand_equivalence_constraints(
            ress, rxn)
        deduped = EnzymeRates._deduplicate(with_eq, rxn)
        @test length(deduped) <= length(with_eq)
    end

    @testset "Stage monotonicity" begin
        rxn = _rxn_ub_reg
        counts = _run_full_pipeline_stages(rxn)
        @test counts.dead_end >= counts.ress
        @test counts.equivalence >= counts.dead_end
        @test counts.dedup <= counts.equivalence
    end

    @testset "Equivalence constraints add variants" begin
        rxn = _rxn_ub_reg
        catalytic = EnzymeRates._catalytic_topologies(rxn)
        de_specs = EnzymeRates._expand_dead_end_inhibitors(
            catalytic, rxn; dead_end_regs=[:I])
        spec_idx = findfirst(
            s -> length(s.edges) > s.n_catalytic_edges,
            de_specs)
        if spec_idx !== nothing
            s = de_specs[spec_idx]
            ress = EnzymeRates._expand_ress_variants(
                [s], rxn)
            with_eq =
                EnzymeRates._expand_equivalence_constraints(
                    ress, rxn)
            @test length(with_eq) >= length(ress)
        end
    end

    @testset "Regulator roles affect partitioning" begin
        rxn_unknown = @enzyme_reaction begin
            substrates: A[X], B[Y]
            products: P[XY]
            regulators: I
        end
        rxn_dead_end = @enzyme_reaction begin
            substrates: A[X], B[Y]
            products: P[XY]
            dead_end_inhibitors: I
        end
        rxn_allosteric = @enzyme_reaction begin
            substrates: A[X], B[Y]
            products: P[XY]
            allosteric_regulators: I
        end

        c_unk = _run_full_pipeline_stages(rxn_unknown)
        c_de = _run_full_pipeline_stages(rxn_dead_end)
        c_al = _run_full_pipeline_stages(rxn_allosteric)

        @test c_unk.catalytic == c_de.catalytic
        @test c_unk.catalytic == c_al.catalytic
        @test c_de.dedup >= c_al.dedup
    end

    # ── param_count accuracy ─────────────────────────────────

    @testset "param_count accuracy (EM)" begin
        rxn = _rxn_ub
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(rxn))
        @test length(all_specs) > 0
        for s in all_specs
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end
    end

    @testset "param_count accuracy (EM with reg, sampled)" begin
        rxn = _rxn_ub_reg
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(rxn))
        rng = Random.MersenneTwister(42)
        n = min(20, length(all_specs))
        sample = all_specs[randperm(rng, length(all_specs))[
            1:n]]
        for s in sample
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end
    end

    @testset "param_count accuracy per stage" begin
        rxn = _rxn_ub  # Uni-Bi, no regulators
        rng = Random.MersenneTwister(99)

        # Stage 1: catalytic topologies
        catalytic = EnzymeRates._catalytic_topologies(rxn)
        for s in catalytic
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end

        # Stage 2: RE/SS expansion
        ress = EnzymeRates._expand_ress_variants(
            catalytic, rxn)
        for s in ress[1:min(10, length(ress))]
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end

        # Stage 3: equivalence constraints (no regs => no
        # dead-end stage changes)
        eq = EnzymeRates._expand_equivalence_constraints(
            ress, rxn)
        for s in eq[1:min(10, length(eq))]
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end

        # With regulators: dead-end + equivalence stages
        rxn_r = _rxn_ub_reg
        cat_r = EnzymeRates._catalytic_topologies(rxn_r)
        ress_r = EnzymeRates._expand_ress_variants(
            cat_r, rxn_r)
        de_r = EnzymeRates._expand_dead_end_inhibitors(
            ress_r, rxn_r; dead_end_regs=[:I])
        sample_de = de_r[randperm(rng, length(de_r))[
            1:min(10, length(de_r))]]
        for s in sample_de
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end

        eq_r = EnzymeRates._expand_equivalence_constraints(
            de_r, rxn_r)
        sample_eq = eq_r[randperm(rng, length(eq_r))[
            1:min(10, length(eq_r))]]
        for s in sample_eq
            m = compile_mechanism(s)
            @test s.param_count == length(parameters(m))
        end
    end

    # ── OEM properties ───────────────────────────────────────

    @testset "OEM expansion properties" begin
        rxn = _rxn_ub
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
                rxn; catalytic_n=2))
        em = filter(
            s -> s isa EnzymeRates.MechanismSpec, all_specs)
        oem = filter(
            s -> s isa EnzymeRates.OligomericMechanismSpec,
            all_specs)
        @test length(em) > 0
        @test length(oem) > 0
        for s in oem
            @test s.catalytic_n == 2
        end
    end

    @testset "OEM with allosteric regulators" begin
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
                _rxn_ub_allo; catalytic_n=2))
        oem = filter(
            s -> s isa EnzymeRates.OligomericMechanismSpec,
            all_specs)
        @test length(oem) > 0
        for s in oem
            @test !isempty(s.allosteric_reg_sites)
            @test !isempty(s.allosteric_multiplicities)
        end
    end

    # ── compile_mechanism round-trip ─────────────────────────

    @testset "compile_mechanism produces valid mechanisms" begin
        rxn = _rxn_ub
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(rxn))
        for s in all_specs
            m = compile_mechanism(s)
            @test m isa EnzymeMechanism
            @test length(metabolites(m)) > 0
            @test length(parameters(m)) > 0
        end
    end

    @testset "compile_mechanism OEM" begin
        rxn = _rxn_ub
        all_specs = collect(
            EnzymeRates.enumerate_mechanisms(
                rxn; catalytic_n=2))
        oem = filter(
            s -> s isa EnzymeRates.OligomericMechanismSpec,
            all_specs)
        for s in oem[1:min(3, length(oem))]
            m = compile_mechanism(s)
            @test m isa OligomericEnzymeMechanism
            @test length(parameters(m)) > 0
        end
    end
end
