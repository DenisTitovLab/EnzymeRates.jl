# ABOUTME: Tests for mechanism enumeration pipeline
# ABOUTME: Unit tests per move + integration tests per reaction

using EnzymeRates: StepSpec, MechanismSpec, AllostericMechanismSpec,
    AbstractMechanismSpec

# Helper: convert EnzymeMechanism → MechanismSpec.
# `m` and `rxn` are dual inputs: `m` carries the catalytic structure,
# `rxn` carries the reaction-level metadata (atoms, declared regulators).
# Validates internally that they're consistent: substrates and products
# names must match exactly; m's regulators must be a subset of rxn's
# declared regulators (allowing the hybrid case where the mechanism
# doesn't yet bind every declared regulator).
function mechanism_spec_from_mechanism_and_rxn(
    m::EnzymeMechanism,
    @nospecialize(rxn::EnzymeReaction))
    m_subs = Set(EnzymeRates.substrates(m))
    rxn_subs = Set(s[1] for s in EnzymeRates.substrates(rxn))
    m_subs == rxn_subs ||
        error("mechanism_spec_from_mechanism_and_rxn: substrate names " *
              "disagree — m=$m_subs, rxn=$rxn_subs")
    m_prods = Set(EnzymeRates.products(m))
    rxn_prods = Set(p[1] for p in EnzymeRates.products(rxn))
    m_prods == rxn_prods ||
        error("mechanism_spec_from_mechanism_and_rxn: product names " *
              "disagree — m=$m_prods, rxn=$rxn_prods")
    m_regs = Set(EnzymeRates.regulators(m))
    rxn_regs = Set(EnzymeRates.regulators(rxn))
    m_regs ⊆ rxn_regs ||
        error("mechanism_spec_from_mechanism_and_rxn: m has regulators " *
              "$(setdiff(m_regs, rxn_regs)) not declared in rxn")

    rxns = EnzymeRates.reactions(m)
    steps = StepSpec[]
    for (lhs, rhs, is_eq, gnum) in rxns
        push!(steps, StepSpec(
            collect(lhs), collect(rhs), is_eq, gnum))
    end
    MechanismSpec(
        rxn, steps,
        length(EnzymeRates.fitted_params(m)))
end

# Helper: build an AllostericMechanismSpec from a compiled
# AllostericEnzymeMechanism and a reaction. Symmetric to
# mechanism_spec_from_mechanism_and_rxn — `m` carries the catalytic
# structure and tags; `rxn` carries reaction-level metadata. The helper
# validates internally that they're consistent: substrate/product names
# match exactly; m's regulators (catalytic + allosteric) are a subset of
# rxn's declared regulators; oligomeric_state(rxn) ==
# catalytic_multiplicity(m). AllostericMechanismSpec uses dense Dict
# storage — every kinetic group and every regulator ligand has an
# explicit entry, so the spec build itself is a pure pass-through.
function allosteric_spec_from_mechanism_and_rxn(
    m::AllostericEnzymeMechanism,
    @nospecialize(rxn::EnzymeReaction))
    cm = EnzymeRates.catalytic_mechanism(m)
    m_subs = Set(EnzymeRates.substrates(m))
    rxn_subs = Set(s[1] for s in EnzymeRates.substrates(rxn))
    m_subs == rxn_subs ||
        error("allosteric_spec_from_mechanism_and_rxn: substrate names " *
              "disagree — m=$m_subs, rxn=$rxn_subs")
    m_prods = Set(EnzymeRates.products(m))
    rxn_prods = Set(p[1] for p in EnzymeRates.products(rxn))
    m_prods == rxn_prods ||
        error("allosteric_spec_from_mechanism_and_rxn: product names " *
              "disagree — m=$m_prods, rxn=$rxn_prods")
    # Allosteric m: regulators(m) returns allosteric-only;
    # regulators(catalytic_mechanism(m)) returns catalytic-side dead-ends.
    m_regs = Set(EnzymeRates.regulators(cm)) ∪ Set(EnzymeRates.regulators(m))
    rxn_regs = Set(EnzymeRates.regulators(rxn))
    m_regs ⊆ rxn_regs ||
        error("allosteric_spec_from_mechanism_and_rxn: m has regulators " *
              "$(setdiff(m_regs, rxn_regs)) not declared in rxn")
    EnzymeRates.oligomeric_state(rxn) == EnzymeRates.catalytic_multiplicity(m) ||
        error("allosteric_spec_from_mechanism_and_rxn: oligomeric_state " *
              "disagrees — m=$(EnzymeRates.catalytic_multiplicity(m)), " *
              "rxn=$(EnzymeRates.oligomeric_state(rxn))")

    base_spec = mechanism_spec_from_mechanism_and_rxn(cm, rxn)

    cat_n = EnzymeRates.catalytic_multiplicity(m)

    group_tags = Dict{Int, Symbol}()
    for g in EnzymeRates.kinetic_groups(EnzymeRates.catalytic_mechanism(m))
        group_tags[g] = EnzymeRates.cat_allo_state(m, g)
    end

    n_reg_sites = length(EnzymeRates.regulatory_sites(m))
    reg_sites = Vector{Symbol}[]
    multiplicities = Int[]
    reg_ligand_tags = Dict{Symbol, Symbol}()
    for i in 1:n_reg_sites
        ligs = collect(EnzymeRates.regulatory_site_ligands(m, i))
        push!(reg_sites, ligs)
        push!(multiplicities,
            EnzymeRates.regulatory_site_multiplicity(m, i))
        for lig in ligs
            reg_ligand_tags[lig] = EnzymeRates.reg_allo_state(m, i, lig)
        end
    end

    AllostericMechanismSpec(
        base_spec, cat_n, reg_sites, multiplicities,
        group_tags, reg_ligand_tags,
        length(EnzymeRates.fitted_params(m)))
end

@testset "AllostericMechanismSpec constructor density validation" begin
    # The constructor rejects sparse Dicts: every kinetic group used in
    # base.steps must have a group_tags entry; every ligand listed in
    # allosteric_reg_sites must have a reg_ligand_tags entry. This guards
    # the dense-storage invariant that both spec and compiled mechanism
    # use throughout the pipeline.
    rxn_uu = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    base = MechanismSpec(rxn_uu,
        [StepSpec([:E, :P], [:E_P], true, 1),
         StepSpec([:E, :S], [:E_S], true, 2),
         StepSpec([:E_S], [:E_P], false, 3)],
        3)

    # Missing group_tags entry (group 3 omitted).
    @test_throws ErrorException AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[], Int[],
        Dict(1 => :EqualRT, 2 => :EqualRT),  # group 3 missing
        Dict{Symbol, Symbol}(),
        4)

    # Missing reg_ligand_tags entry (ligand R declared in reg_sites but
    # not in reg_ligand_tags).
    @test_throws ErrorException AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[[:R]], Int[2],
        Dict(1 => :EqualRT, 2 => :EqualRT, 3 => :EqualRT),
        Dict{Symbol, Symbol}(),  # :R missing
        4)

    # Both Dicts complete → constructor succeeds.
    valid = AllostericMechanismSpec(
        base, 2,
        Vector{Symbol}[[:R]], Int[2],
        Dict(1 => :EqualRT, 2 => :EqualRT, 3 => :EqualRT),
        Dict(:R => :EqualRT),
        4)
    @test valid isa AllostericMechanismSpec
end

@testset "spec-from-mechanism helpers reject inconsistent inputs" begin
    # Both helpers validate that mechanism and reaction agree on
    # substrates/products/regulators (and oligomeric_state for the
    # allosteric variant). Mismatches throw ErrorException at the helper
    # level, before any spec is constructed.

    m_uu = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
        end
    end

    # Substrate name mismatch
    rxn_wrong_sub = @enzyme_reaction begin
        substrates: X[C]   # not :S
        products: P[C]
    end
    @test_throws ErrorException mechanism_spec_from_mechanism_and_rxn(m_uu, rxn_wrong_sub)

    # Product name mismatch
    rxn_wrong_prod = @enzyme_reaction begin
        substrates: S[C]
        products: Y[C]   # not :P
    end
    @test_throws ErrorException mechanism_spec_from_mechanism_and_rxn(m_uu, rxn_wrong_prod)

    # Regulator subset rule: m_with_I has :I bound; rxn lacks :I → reject.
    m_with_I = @enzyme_mechanism begin
        substrates: S
        products: P
        regulators: I
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
            E + I ⇌ E_I
        end
    end
    rxn_no_I = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    @test_throws ErrorException mechanism_spec_from_mechanism_and_rxn(m_with_I, rxn_no_I)

    # Allosteric helper: oligomeric_state mismatch
    m_allo_2 = @allosteric_mechanism begin
        substrates: S
        products: P
        site(:catalytic, 2): begin
            steps: begin
                E + P ⇌ E_P    :: EqualRT
                E + S ⇌ E_S    :: EqualRT
                E_S <--> E_P   :: EqualRT
            end
        end
    end
    rxn_oligo_4 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 4   # mismatch — m_allo_2 has catalytic_n=2
    end
    @test_throws ErrorException allosteric_spec_from_mechanism_and_rxn(m_allo_2, rxn_oligo_4)
end

@testset "allosteric_spec_from_mechanism_and_rxn round-trip" begin
    # K-type uni-uni: catalytic 2-mer, all bindings :EqualRT, iso :EqualRT,
    # no regulators. Round-trip must be lossless: spec → AllostericEnzymeMechanism
    # rebuilds to the same singleton type as the macro produced.
    m1 = @allosteric_mechanism begin
        substrates: S
        products: P
        site(:catalytic, 2): begin
            steps: begin
                E + S ⇌ E_S       :: EqualRT
                E + P ⇌ E_P       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    rxn1 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    spec1 = allosteric_spec_from_mechanism_and_rxn(m1, rxn1)
    @test AllostericEnzymeMechanism(spec1) === m1

    # Mixed group tags: one :OnlyR, one :EqualRT, one :NonequalRT.
    # Dense storage: every group has an explicit entry in group_tags.
    m2 = @allosteric_mechanism begin
        substrates: S
        products: P
        site(:catalytic, 2): begin
            steps: begin
                E + S ⇌ E_S       :: OnlyR
                E + P ⇌ E_P       :: EqualRT
                E_S <--> E_P      :: NonequalRT
            end
        end
    end
    spec2 = allosteric_spec_from_mechanism_and_rxn(m2, rxn1)
    @test AllostericEnzymeMechanism(spec2) === m2
    @test spec2.group_tags == Dict(1 => :OnlyR, 2 => :EqualRT, 3 => :NonequalRT)

    # With one allosteric regulator at its own site, tag :OnlyT.
    rxn3 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    m3 = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R::OnlyT
        site(:catalytic, 2): begin
            steps: begin
                E + S ⇌ E_S       :: EqualRT
                E + P ⇌ E_P       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
    end
    spec3 = allosteric_spec_from_mechanism_and_rxn(m3, rxn3)
    @test AllostericEnzymeMechanism(spec3) === m3
    @test spec3.reg_ligand_tags == Dict(:R => :OnlyT)

    # Two regulators at the same site, one :OnlyR, one :NonequalRT.
    # Dense storage: both ligands appear in reg_ligand_tags.
    rxn4 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R1, R2
        oligomeric_state: 2
    end
    m4 = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R1::OnlyR, R2::NonequalRT
        site(:catalytic, 2): begin
            steps: begin
                E + S ⇌ E_S       :: EqualRT
                E + P ⇌ E_P       :: EqualRT
                E_S <--> E_P      :: EqualRT
            end
        end
        site(:regulatory, 2): begin
            ligands: R1, R2
        end
    end
    spec4 = allosteric_spec_from_mechanism_and_rxn(m4, rxn4)
    @test AllostericEnzymeMechanism(spec4) === m4
    @test spec4.reg_ligand_tags == Dict(:R1 => :OnlyR, :R2 => :NonequalRT)
    @test spec4.allosteric_multiplicities == [2]
    @test spec4.allosteric_reg_sites == [[:R1, :R2]]
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

# Pyruvate carboxylase: Pyr + HCO3 + ATP = OAA + ADP + Pi
# Mechanism: ATP+HCO3 → ADP+Pi+CO2_residual,
#            then Pyr+CO2 → OAA
const pyruvate_carboxylase_rxn = @enzyme_reaction begin
    substrates: Pyr[C3H3O3], HCO3[HCO3], ATP[C10H16N5O13P3]
    products: OAA[C4H3O5], ADP[C10H15N5O10P2], Pi[H2PO4]
end

# Pyruvate dehydrogenase: Pyr + NAD + CoA = AcCoA + NADH + CO2
# Mechanism: Pyr → CO2+residual, CoA+residual → AcCoA+residual,
#            NAD+residual → NADH
const pyruvate_dehydrogenase_rxn = @enzyme_reaction begin
    substrates: Pyr[C3H3O3], NAD[C21H28N7O14P2], CoA[C21H36N7O16P3S]
    products: AcCoA[C23H38N7O17P3S], NADH[C21H29N7O14P2], CO2[CO2]
end

"""Post-condition every spec must satisfy after any expansion move.
Catches the broad class of "expansion move forgot to update a
parallel data structure" bugs.

For AllostericMechanismSpec, asserts the dense-storage invariants:
- every kinetic group used in `base.steps` has an explicit entry in
  `group_tags`,
- every ligand listed in `allosteric_reg_sites` has an explicit
  entry in `reg_ligand_tags`.

Universal invariant — holds after every `_expand_*` move's output,
including `_expand_change_allo_state` (which now writes
`:NonequalRT` explicitly rather than `delete!`ing).
"""
function _assert_spec_invariants(spec::MechanismSpec)
    @test spec.n_fit_params_estimate >= 0
end

function _assert_spec_invariants(spec::AllostericMechanismSpec)
    @test spec.n_fit_params_estimate >= 0
    used_groups = Set(s.kinetic_group for s in spec.base.steps)
    for g in used_groups
        @test haskey(spec.group_tags, g)
    end
    for site in spec.allosteric_reg_sites
        for lig in site
            @test haskey(spec.reg_ligand_tags, lig)
        end
    end
end

"""Collect all mechanisms by running the full enumeration loop."""
function enumerate_all(
    @nospecialize(reaction::EnzymeReaction);
    max_params::Int=20)
    cache = Dict{Int, Vector{EnzymeRates.AbstractMechanismSpec}}()

    init_specs = EnzymeRates.init_mechanisms(reaction)
    for spec in init_specs
        push!(get!(cache, spec.n_fit_params_estimate,
            EnzymeRates.AbstractMechanismSpec[]), spec)
    end
    EnzymeRates.dedup!(cache)

    results = Dict{Int, Vector{EnzymeRates.AbstractMechanismSpec}}()

    for pc in minimum(keys(cache)):max_params
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

# ═══════════════════════════════════════════════════════════════════════
# 1. Support functions (no spec input)
# ═══════════════════════════════════════════════════════════════════════

# ─── _catalytic_topologies ──────────────────────────────────────────────
@testset "_catalytic_topologies" begin

    @testset "Uni-Uni" begin
        topos = EnzymeRates._catalytic_topologies(uni_uni_rxn)
        @test length(topos) == 1

        m_uu = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec_rt = mechanism_spec_from_mechanism_and_rxn(
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
        # 11 = 9 sequential + 2 empty-residual ping-pong
        @test length(topos) == 11
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
        @test length(topos) == 283
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "Ter-Bi" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_bi_rxn)
        # 51 = 39 sequential + 6 nonempty-residual +
        # 6 empty-residual ping-pong
        @test length(topos) == 51
        for t in topos
            @test count(
                !s.is_equilibrium for s in t.steps) == 1
        end
    end

    @testset "empty-residual ping-pong" begin
        ter_ter = @enzyme_reaction begin
            substrates: A[C], B[N], D[X]
            products: P[C], Q[N], R[X]
        end
        topos = EnzymeRates._catalytic_topologies(ter_ter)
        has_estar = any(topos) do spec
            any(spec.steps) do s
                any(
                    sym -> startswith(
                        string(sym), "Estar"),
                    Iterators.flatten(
                        (s.reactants, s.products)),
                )
            end
        end
        @test has_estar
    end

    @testset "weak-ordering combining" begin
        # For bi-bi: 9 sequential + 2 empty-residual
        # ping-pong = 11
        bi_bi_rxn_test = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
        end
        topos = EnzymeRates._catalytic_topologies(
            bi_bi_rxn_test)
        @test length(topos) == 11

        topos_tt = EnzymeRates._catalytic_topologies(
            ter_ter_rxn)
        @test length(topos_tt) == 283
    end

    @testset "isomerization constraints" begin
        topos = EnzymeRates._catalytic_topologies(
            ter_ter_rxn)

        met_names = Set([:A, :B, :D, :P, :Q, :R])
        sub_names_set = Set([:A, :B, :D])

        # C5: max bound metabolites = max(3,3) = 3
        for spec in topos
            for s in spec.steps
                for sym_list in (s.reactants, s.products)
                    for sym in sym_list
                        str = replace(
                            string(sym),
                            "Estar" => "E")
                        parts = split(str, "_")
                        n_mets = count(
                            p -> Symbol(p) ∈ met_names,
                            parts)
                        @test n_mets <= 3
                    end
                end
            end
        end

        # C7: every iso source must contain at
        # least one substrate
        for spec in topos
            for s in spec.steps
                if length(s.reactants) == 1 &&
                        length(s.products) == 1
                    src = string(s.reactants[1])
                    src_parts = Symbol.(split(
                        replace(src, "Estar" => "E"),
                        "_"))
                    has_sub = any(
                        p -> p ∈ sub_names_set,
                        src_parts)
                    @test has_sub
                end
            end
        end

        # C8: iso product forms should not contain
        # substrate names
        for spec in topos
            for s in spec.steps
                if length(s.reactants) == 1 &&
                        length(s.products) == 1
                    dst = string(s.products[1])
                    if startswith(dst, "Estar_")
                        dst_parts = Symbol.(
                            split(dst[7:end], "_"))
                        for p in dst_parts
                            @test p ∉ sub_names_set
                        end
                    end
                end
            end
        end
    end

    @testset "pyruvate carboxylase mechanism" begin
        topos = EnzymeRates._catalytic_topologies(
            pyruvate_carboxylase_rxn)

        # Known mechanism: ATP+HCO3 → ADP+Pi (CO2 residual),
        # then Pyr+CO2 → OAA
        found = false
        for spec in topos
            iso_steps = [
                (sort(s.reactants), sort(s.products))
                for s in spec.steps
                if length(s.reactants) == 1 &&
                    length(s.products) == 1
            ]
            has_atp_hco3_iso = any(iso_steps) do (r, p)
                r == [Symbol("E_ATP_HCO3")] &&
                    p == [Symbol("Estar_ADP_Pi")]
            end
            has_pyr_iso = any(iso_steps) do (r, p)
                r == [Symbol("Estar_Pyr")] &&
                    p == [Symbol("E_OAA")]
            end
            if has_atp_hco3_iso && has_pyr_iso
                found = true
                break
            end
        end
        @test found

        # 312 = 169 seq + 143 pp
        @test length(topos) == 312
        seq_count = count(topos) do spec
            !any(spec.steps) do s
                any(sym -> startswith(string(sym), "Estar"),
                    Iterators.flatten(
                        (s.reactants, s.products)))
            end
        end
        pp_count = length(topos) - seq_count
        @test seq_count == 169
        @test pp_count == 143
    end

    @testset "pyruvate dehydrogenase mechanism" begin
        topos = EnzymeRates._catalytic_topologies(
            pyruvate_dehydrogenase_rxn)

        # Known mechanism: Pyr→CO2 (residual C2H3O),
        # CoA+residual→AcCoA (residual H),
        # NAD+residual→NADH (no residual)
        found = false
        for spec in topos
            iso_steps = [
                (sort(s.reactants), sort(s.products))
                for s in spec.steps
                if length(s.reactants) == 1 &&
                    length(s.products) == 1
            ]
            has_pyr = any(iso_steps) do (r, p)
                r == [:E_Pyr] && p == [:Estar_CO2]
            end
            has_coa = any(iso_steps) do (r, p)
                r == [:Estar_CoA] && p == [:Estar_AcCoA]
            end
            has_nad = any(iso_steps) do (r, p)
                r == [:Estar_NAD] && p == [:E_NADH]
            end
            if has_pyr && has_coa && has_nad
                found = true
                break
            end
        end
        @test found

        # 334 = 169 seq + 165 pp
        @test length(topos) == 334
        seq_count = count(topos) do spec
            !any(spec.steps) do s
                any(sym -> startswith(string(sym), "Estar"),
                    Iterators.flatten(
                        (s.reactants, s.products)))
            end
        end
        pp_count = length(topos) - seq_count
        @test seq_count == 169
        @test pp_count == 165
    end

    @testset "quad-quad: C6 forces ping-pong" begin
        # Quad-quad reaction: 4 subs, 4 prods
        # With C6 (iso ≤ 3×3), 4→4 sequential iso is blocked
        # All topologies must use ping-pong (at least 2 iso steps)
        quad_rxn = @enzyme_reaction begin
            substrates: A[C], B[N], D[X], F[Y]
            products: P[C], Q[N], R[X], S[Y]
        end
        topos = EnzymeRates._catalytic_topologies(quad_rxn)
        @test length(topos) > 0
        # Every topology must have ≥ 2 iso steps (no 4→4)
        for spec in topos
            n_iso = count(spec.steps) do s
                length(s.reactants) == 1 &&
                    length(s.products) == 1
            end
            @test n_iso >= 2
        end
    end

end

# ─── _competition_patterns ──────────────────────────────────────────────
@testset "_competition_patterns" begin
    # Uni-uni: 1×1, only 1 pattern (single edge)
    pats_11 = EnzymeRates._competition_patterns(
        Set([:S]), Set([:P]))
    @test length(pats_11) == 1
    @test pats_11[1] == Set([(:S, :P)])

    # Uni-bi: 1×2, S competes with both P and Q
    pats_12 = EnzymeRates._competition_patterns(
        Set([:S]), Set([:P, :Q]))
    @test length(pats_12) == 1
    @test pats_12[1] == Set([(:S, :P), (:S, :Q)])

    # Bi-uni: symmetric
    pats_21 = EnzymeRates._competition_patterns(
        Set([:A, :B]), Set([:P]))
    @test length(pats_21) == 1
    @test pats_21[1] == Set([(:A, :P), (:B, :P)])

    # Bi-bi: 7 patterns
    pats_22 = EnzymeRates._competition_patterns(
        Set([:A, :B]), Set([:P, :Q]))
    @test length(pats_22) == 7
    # Every pattern covers all vertices
    for pat in pats_22
        for s in [:A, :B]
            @test any(
                p -> (s, p) in pat, [:P, :Q])
        end
        for p in [:P, :Q]
            @test any(
                s -> (s, p) in pat, [:A, :B])
        end
    end
    # Invalid: {A↔P, B↔P} leaves Q uncovered
    @test Set([(:A, :P), (:B, :P)]) ∉ pats_22

    # Ter-ter: 265 patterns
    pats_33 = EnzymeRates._competition_patterns(
        Set([:A, :B, :C]),
        Set([:P, :Q, :R]))
    @test length(pats_33) == 265
    for pat in pats_33
        for s in [:A, :B, :C]
            @test any(
                p -> (s, p) in pat,
                [:P, :Q, :R])
        end
        for p in [:P, :Q, :R]
            @test any(
                s -> (s, p) in pat,
                [:A, :B, :C])
        end
    end
end

# ─── _inhibitor_competition_patterns ────────────────────────────────────
@testset "_inhibitor_competition_patterns" begin
    # Uni-uni, no existing inhibitors
    pats = EnzymeRates._inhibitor_competition_patterns(
        Set([:S]), Set([:P]), Symbol[])
    @test length(pats) == 1
    @test pats[1] == (Set([:S]), Set([:P]), Set{Symbol}())

    # Bi-bi, no existing inhibitors: 3×3 = 9
    pats_bb = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]), Symbol[])
    @test length(pats_bb) == 9

    # Ter-ter, no existing inhibitors: 7×7 = 49
    pats_tt = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B, :C]), Set([:P, :Q, :R]), Symbol[])
    @test length(pats_tt) == 49

    # Bi-bi, 1 existing inhibitor: 9 × 2 = 18
    pats_1i = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]), [:I1__reg])
    @test length(pats_1i) == 18

    # Bi-bi, 2 existing inhibitors: 9 × 4 = 36
    pats_2i = EnzymeRates._inhibitor_competition_patterns(
        Set([:A, :B]), Set([:P, :Q]),
        [:I1__reg, :I2__reg])
    @test length(pats_2i) == 36
end

# ─── _forms_with_binding_step ───────────────────────────────────────────
@testset "_forms_with_binding_step" begin
    # Uni-uni: S binds to E, P binds to E
    m_uu = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
        end
    end
    spec_uu = mechanism_spec_from_mechanism_and_rxn(
        m_uu, uni_uni_rxn)
    @test EnzymeRates._forms_with_binding_step(
        spec_uu.steps, :S) == Set([:E])
    @test EnzymeRates._forms_with_binding_step(
        spec_uu.steps, :P) == Set([:E])

    # Bi-bi random: B binds to E and E_A
    m_bb = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E_A
            E_B + A ⇌ E_A_B
            E + B ⇌ E_B
            E_A + B ⇌ E_A_B
            E + P ⇌ E_P
            E_P + Q ⇌ E_P_Q
            E + Q ⇌ E_Q
            E_Q + P ⇌ E_P_Q
            E_A_B <--> E_P_Q
        end
    end
    spec_bb = mechanism_spec_from_mechanism_and_rxn(
        m_bb, bi_bi_rxn)
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :B) == Set([:E, :E_A])
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :A) == Set([:E, :E_B])
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :P) == Set([:E, :E_Q])
    @test EnzymeRates._forms_with_binding_step(
        spec_bb.steps, :Q) == Set([:E, :E_P])
end

# ─── _substrate_product_dead_end_opportunities ──────────────────────────
@testset "_substrate_product_dead_end_opportunities" begin
    # Random ter-ter topology has 27 possible
    # dead-end forms. With diagonal competition
    # {A↔P, B↔Q, D↔R}:
    #   1S+1P: 6 allowed, 3 forbidden (the 3
    #     diagonal pairs A-P, B-Q, D-R)
    #   2S+1P: 3 allowed (E_A_B_R, E_A_D_Q,
    #     E_B_D_P), 6 forbidden
    #   1S+2P: 3 allowed (E_A_Q_R, E_B_P_R,
    #     E_D_P_Q), 6 forbidden
    #   Total: 12 allowed out of 27
    topos = EnzymeRates._catalytic_topologies(
        ter_ter_rxn)
    # Pick a sequential topo with most forms
    _, idx = findmax(
        length(EnzymeRates.all_form_names(t))
        for t in topos)
    random_topo = topos[idx]
    bound =
        EnzymeRates._bound_metabolites_at_forms(
            random_topo, ter_ter_rxn)
    sub_names = Set([:A, :B, :D])
    prod_names = Set([:P, :Q, :R])
    cat_forms =
        EnzymeRates.all_form_names(random_topo)
    _sp_de_opps =
        EnzymeRates._substrate_product_dead_end_opportunities
    de_opps = _sp_de_opps(
        bound, cat_forms,
        sub_names, prod_names)
    # Group dead-end forms
    de_forms = Dict{Symbol,
        Vector{Tuple{Symbol, Symbol}}}()
    for (f, m) in de_opps
        de_name =
            EnzymeRates._dead_end_form_name(
                f, bound[f], m)
        push!(get!(de_forms, de_name,
            Tuple{Symbol, Symbol}[]), (f, m))
    end
    de_form_names =
        sort(collect(keys(de_forms)))
    @test length(de_form_names) == 27

    # Build de_bound mapping
    de_bound = Dict{Symbol, Set{Symbol}}()
    for de_name in de_form_names
        f, m = first(de_forms[de_name])
        de_bound[de_name] = union(
            bound[f], Set([m]))
    end

    # Apply diagonal competition filter
    diagonal =
        Set([(:A, :P), (:B, :Q), (:D, :R)])
    allowed = Symbol[]
    for de_name in de_form_names
        mets = de_bound[de_name]
        de_subs = intersect(mets, sub_names)
        de_prods =
            intersect(mets, prod_names)
        has_conflict = any(
            (s, p) in diagonal
            for s in de_subs
            for p in de_prods)
        has_conflict ||
            push!(allowed, de_name)
    end
    @test length(allowed) == 12

    # Verify specific allowed forms
    @test :E_A_Q in allowed   # 1S+1P
    @test :E_A_R in allowed
    @test :E_B_P in allowed
    @test :E_B_R in allowed
    @test :E_D_P in allowed
    @test :E_D_Q in allowed
    @test :E_A_B_R in allowed # 2S+1P
    @test :E_A_D_Q in allowed
    @test :E_B_D_P in allowed
    @test :E_A_Q_R in allowed # 1S+2P
    @test :E_B_P_R in allowed
    @test :E_D_P_Q in allowed

    # Verify specific forbidden forms
    @test :E_A_P ∉ allowed    # A↔P diagonal
    @test :E_B_Q ∉ allowed    # B↔Q diagonal
    @test :E_D_R ∉ allowed    # D↔R diagonal
end

# ─── _expand_substrate_product_dead_ends ────────────────────────────────
@testset "_expand_substrate_product_dead_ends" begin

    @testset "Uni-Uni: no dead-end forms" begin
        # 3 forms: E, E_S[C], E_P[C]. E_S has all subs,
        # E_P has all prods. No mixed dead-end possible.
        # → 0 dead-end forms, 1 variant (bare topology)
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, uni_uni_rxn)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], uni_uni_rxn)
        @test length(result) == 1
    end

    @testset "Bi-Bi random: 4 dead-end forms" begin
        # 7 forms: E, E_A, E_B, E_A_B, E_P, E_Q, E_P_Q
        # Eligible dead-end forms (mixed sub+prod binding):
        #   E_A: +P→E_A_P(mixed✓), +Q→E_A_Q(mixed✓)
        #   E_B: +P→E_B_P(mixed✓), +Q→E_B_Q(mixed✓)
        #   E_P: +A→E_A_P(same), +B→E_B_P(same)
        #   E_Q: +A→E_A_Q(same), +B→E_B_Q(same)
        # 4 unique mixed-substrate-product forms across competition patterns.
        # Competition patterns for bi-bi (2 subs × 2 prods): 7 patterns
        # (the count from _competition_patterns(2, 2)). Each pattern produces
        # a distinct dead-end-form set:
        #   {A↔P, B↔Q}: forbids E_A_P, E_B_Q → emits {E_A_Q, E_B_P}
        #   {A↔Q, B↔P}: forbids E_A_Q, E_B_P → emits {E_A_P, E_B_Q}
        #   ... (one set per pattern, all distinct)
        #   {A↔P, A↔Q, B↔P, B↔Q}: forbids all → emits {} (bare topology)
        # All 7 sets are distinct → 7 variants after dedup.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_B + A ⇌ E_A_B
                E + B ⇌ E_B
                E_A + B ⇌ E_A_B
                E + P ⇌ E_P
                E_P + Q ⇌ E_P_Q
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, bi_bi_rxn)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], bi_bi_rxn)
        # 4 unique dead-end forms, 7 competition patterns,
        # all 7 produce distinct dead-end sets → 7 variants
        @test length(result) == 7
    end

    @testset "Uni-Bi ordered: no dead-end forms" begin
        # 4 forms: E, E_S, E_P_Q, E_Q
        # E+P→E_P: single-product → rejected (need mixed)
        # E_Q+S→E_S_Q: has all subs → rejected
        # → 0 dead-end forms, 1 variant
        m = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E + S ⇌ E_S
                E_S <--> E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, uni_bi_rxn)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], uni_bi_rxn)
        @test length(result) == 1
    end

    @testset "Bi-Bi Ping-Pong: 5 dead-end forms → 7 variants" begin
        # Forms: E, E_A, Estar, Estar_A_P, Estar_B, E_Q
        # 5 dead-end forms total (E-side: E_A_P, E_A_Q, E_B_Q; Estar-side:
        # Estar_B_P, Estar_B_Q). 7 competition patterns; each yields a
        # distinct dead-end-form set after dedup → 7 variants.
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                Estar + B ⇌ Estar_B
                E + Q ⇌ E_Q
                Estar + P ⇌ Estar_A_P
                E_A <--> Estar_A_P
                Estar_B ⇌ E_Q
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(
            m, bi_bi_pp_rxn)
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], bi_bi_pp_rxn)
        # 5 dead-end forms (E_A_P, E_A_Q, E_B_Q from
        # E-side + Estar_B_P, Estar_B_Q from
        # Estar-side), competition-filtered
        @test length(result) == 7
    end

    @testset "Dead-end filtering by competition" begin

        # Shared bi-bi random mechanism for multiple tests
        m_bb = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_B + A ⇌ E_A_B
                E + B ⇌ E_B
                E_A + B ⇌ E_A_B
                E + P ⇌ E_P
                E_P + Q ⇌ E_P_Q
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        spec_bb = mechanism_spec_from_mechanism_and_rxn(
            m_bb, bi_bi_rxn)

        # Note: the variant count assertion (length == 7) is covered by
        # the "Bi-Bi random: 4 dead-end forms" sub-testset higher in this
        # same parent testset. The sub-testsets below probe the SHAPE of
        # those 7 variants — which forms appear in each, and which
        # patterns produce empty/full dead-end sets.

        @testset "Bi-bi random: complete competition → bare topology" begin
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec_bb], bi_bi_rxn)
            # Complete pattern {A↔P,A↔Q,B↔P,B↔Q} forbids
            # all dead-end forms → 1 variant has no dead-end
            # steps (same step count as original)
            bare = filter(
                r -> length(r.steps) == length(spec_bb.steps),
                result)
            @test length(bare) == 1
        end

        @testset "Bi-bi random: diagonal has exactly 2 dead-end forms" begin
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec_bb], bi_bi_rxn)
            # Diagonal patterns {A↔P,B↔Q} and {A↔Q,B↔P}
            # each allow exactly 2 dead-end forms.
            two_de = filter(result) do r
                de_forms = setdiff(
                    EnzymeRates.all_form_names(r),
                    EnzymeRates.all_form_names(spec_bb))
                length(de_forms) == 2
            end
            @test length(two_de) == 2  # diagonal + anti-diagonal
        end

        @testset "Ter-ter per-topology (OOM on full init)" begin
            # Test that competition filtering works
            # on representative ter-ter topologies.
            topos = EnzymeRates._catalytic_topologies(
                ter_ter_rxn)
            @test length(topos) == 283
            # Test first (random, most forms) and last topology
            for topo in [topos[1], topos[end]]
                result =
                    EnzymeRates._expand_substrate_product_dead_ends(
                        [topo], ter_ter_rxn)
                # Competition patterns reduce 2^27 to
                # ≤265 variants per topology
                @test length(result) > 0
                @test length(result) <= 265
                for spec in result
                    # Verify n_fit_params_estimate set correctly
                    @test spec.n_fit_params_estimate >=
                        length(EnzymeRates.substrates(
                            ter_ter_rxn)) +
                        length(EnzymeRates.products(
                            ter_ter_rxn)) + 1
                end
            end
        end
    end

end

# ═══════════════════════════════════════════════════════════════════════
# Testsets covering non-enumeration features (atom balance from
# @enzyme_reaction; AllostericEnzymeMechanism accessor identity)
# ═══════════════════════════════════════════════════════════════════════

@testset "test reaction atom balance" begin
    for rxn in [pyruvate_carboxylase_rxn,
                pyruvate_dehydrogenase_rxn]
        sub_atoms = Dict{Symbol,Int}()
        for (_, atoms) in EnzymeRates.substrates(rxn)
            for (a, c) in atoms
                sub_atoms[a] = get(sub_atoms, a, 0) + c
            end
        end
        prod_atoms = Dict{Symbol,Int}()
        for (_, atoms) in EnzymeRates.products(rxn)
            for (a, c) in atoms
                prod_atoms[a] = get(prod_atoms, a, 0) + c
            end
        end
        @test sub_atoms == prod_atoms
    end
end

@testset "AllostericEnzymeMechanism TR equivalence" begin
    base_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    # Steps with kinetic_group: 1 = S binding (RE),
    # 2 = P binding (RE), 3 = iso (SS).
    base_steps = [
        StepSpec([:E, :S], [:E_S], true, 1),
        StepSpec([:E, :P], [:E_P], true, 2),
        StepSpec([:E_S], [:E_P], false, 3),
    ]
    base_spec = MechanismSpec(base_rxn, base_steps, 3)

    # S binding (group 1) is `:EqualRT` (K_T_S = K_R_S),
    # P binding (group 2) is `:OnlyR` (absent from T-state).
    # Iso (group 3) is `:NonequalRT` (default).
    allo_spec = AllostericMechanismSpec(
        base_spec, 2,
        Vector{Symbol}[], Int[],
        Dict(1 => :EqualRT, 2 => :OnlyR, 3 => :NonequalRT),
        Dict{Symbol, Symbol}(),
        base_spec.n_fit_params_estimate + 1)

    m_compiled = AllostericEnzymeMechanism(allo_spec)
    @test EnzymeRates.cat_allo_state(m_compiled, 1) == :EqualRT
    @test EnzymeRates.cat_allo_state(m_compiled, 2) == :OnlyR
end

# ═══════════════════════════════════════════════════════════════════════
# 2. Initialization (compile_mechanism + init_mechanisms)
# ═══════════════════════════════════════════════════════════════════════

# ─── compile_mechanism / EnzymeMechanism round-trip ────────────────────
@testset "compile_mechanism round-trip" begin
    # Round-trip lossless invariant: for any mechanism built via the DSL,
    # mechanism_spec_from_mechanism_and_rxn ∘ EnzymeMechanism (== compile_mechanism)
    # returns the same singleton type. Validates the helper AND the
    # constructor's bidirectional consistency.

    # uni-uni
    m_uu = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
        end
    end
    spec_uu = mechanism_spec_from_mechanism_and_rxn(m_uu, uni_uni_rxn)
    @test EnzymeMechanism(spec_uu) === m_uu

    # bi-bi sequential
    m_seq = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E_A
            E_A + B ⇌ E_A_B
            E + Q ⇌ E_Q
            E_Q + P ⇌ E_P_Q
            E_A_B <--> E_P_Q
        end
    end
    spec_seq = mechanism_spec_from_mechanism_and_rxn(m_seq, bi_bi_rxn)
    @test EnzymeMechanism(spec_seq) === m_seq

    # bi-bi ping-pong
    m_pp = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E_A
            Estar + B ⇌ Estar_B
            E + Q ⇌ E_Q
            Estar + P ⇌ Estar_A_P
            E_A <--> Estar_A_P
            Estar_B ⇌ E_Q
        end
    end
    spec_pp = mechanism_spec_from_mechanism_and_rxn(m_pp, bi_bi_pp_rxn)
    @test EnzymeMechanism(spec_pp) === m_pp

    # uni-uni with dead-end inhibitor (regulator strip in the round-trip)
    m_uu_i = @enzyme_mechanism begin
        substrates: S
        products: P
        regulators: I
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
            E + I ⇌ E_I
        end
    end
    spec_uu_i = mechanism_spec_from_mechanism_and_rxn(m_uu_i, uni_uni_with_reg)
    @test EnzymeMechanism(spec_uu_i) === m_uu_i
end

# ─── init_mechanisms ───────────────────────────────────────────────────
@testset "init_mechanisms" begin

    @testset "min param count floor: subs + prods + 1" begin
        # Floor invariant: n_fit_params_estimate ≥ n_subs + n_prods + 1.
        # Derivation: every init mechanism has 1 SS step (the iso) plus
        # n_subs RE binding groups for substrates and n_prods for products.
        # The kinetic-group count = n_subs + n_prods (RE) + 1 (SS), and
        # n_thermo subtracts off based on cycle counts. The minimum after
        # subtractions equals exactly n_subs + n_prods + 1 for the simplest
        # topology with no dead-ends.
        for (rxn, n_s, n_p) in [
            (uni_uni_rxn, 1, 1),
            (uni_bi_rxn, 1, 2),
            (bi_bi_rxn, 2, 2),
            (bi_bi_pp_rxn, 2, 2),
        ]
            specs = EnzymeRates.init_mechanisms(rxn)
            min_pc = n_s + n_p + 1
            for s in specs
                @test s.n_fit_params_estimate >= min_pc
            end
        end
    end

    @testset "n_fit_params_estimate matches fitted_params for uni-uni init" begin
        # Uni-uni: 3 forms (E, E_S, E_P), 3 steps. n_thermo = 3 - 3 + 1 = 1
        # (one independent thermodynamic constraint = Keq).
        # Formula: n_re_groups + 2*n_ss_groups - n_thermo = 2 + 2 - 1 = 3.
        # length(fitted_params(m)) for uni-uni init = 3 (K1, K2, k3f).
        # Estimate must equal actual on the simplest case.
        init_specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        @test !isempty(init_specs)
        spec = first(init_specs)
        m = EnzymeRates.compile_mechanism(spec)
        n_actual = length(EnzymeRates.fitted_params(m))
        @test spec.n_fit_params_estimate == n_actual
    end

    @testset "n_fit_params_estimate upper-bound for dead-end mirrors" begin
        # When dead-end mirror cycles exist, the formula can underestimate
        # the true thermodynamic-constraint count, so the floor in
        # _apply_equivalence_grouping ensures pc ≥ n_subs + n_prods + 1.
        # This guards the upper-bound invariant: estimate ≥ actual.
        # Cap compiled specs to keep @generated cost bounded.
        cap = 30
        init_specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        for spec in init_specs[1:min(cap, end)]
            m = EnzymeRates.compile_mechanism(spec)
            @test spec.n_fit_params_estimate >=
                length(EnzymeRates.fitted_params(m))
        end
        expanded = EnzymeRates.expand_mechanisms(
            init_specs, uni_uni_with_reg)
        expanded_specs = EnzymeRates.AbstractMechanismSpec[]
        for (_, specs) in expanded
            append!(expanded_specs, specs)
        end
        for spec in expanded_specs[1:min(cap, end)]
            m = EnzymeRates.compile_mechanism(spec)
            @test spec.n_fit_params_estimate >=
                length(EnzymeRates.fitted_params(m))
        end
    end

    @testset "exactly 1 SS step per init spec" begin
        # init_mechanisms produces minimum-parameter mechanisms — exactly
        # one isomerization step, which is SS by construction. Subsequent
        # RE→SS expansions add more SS steps; init never does.
        for rxn in [uni_uni_rxn, uni_bi_rxn,
                    bi_bi_rxn, bi_bi_pp_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            for s in specs
                @test count(st -> !st.is_equilibrium, s.steps) == 1
            end
        end
    end

    @testset "Same-metabolite RE bindings share kinetic_group" begin
        # _apply_equivalence_grouping collapses all RE binding steps for
        # the same metabolite into one kinetic group (one shared K).
        # For uni-uni + dead-end inhibitor, the inhibitor's mirror cycles
        # mean :I binds at multiple forms — these mirror bindings must
        # share a single kinetic_group (one K_I).
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            dead_end_inhibitors: I
        end
        specs = EnzymeRates.init_mechanisms(rxn)
        @test !isempty(specs)
        for spec in specs
            by_metabolite = Dict{Symbol, Vector{EnzymeRates.StepSpec}}()
            for step in spec.steps
                step.is_equilibrium || continue
                length(step.reactants) == 2 || continue
                met = step.reactants[2]
                push!(get!(by_metabolite, met,
                           EnzymeRates.StepSpec[]), step)
            end
            for (_met, steps) in by_metabolite
                length(steps) >= 2 || continue
                groups = Set(s.kinetic_group for s in steps)
                @test length(groups) == 1
            end
        end
    end

    @testset "Uni-uni: exactly 1 init mechanism" begin
        # Uni-uni topology: 1 catalytic topology × 1 dead-end variant
        # (none possible — see test_expand_substrate_product_dead_ends
        # uni-uni case). Hence init produces exactly 1 spec.
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        @test length(specs) == 1
    end

    @testset "Init compiles for all small reactions" begin
        # Every init spec must compile to a valid EnzymeMechanism, and
        # the actual fitted-param count must respect the upper-bound
        # invariant. Tests first 5 specs per reaction to cap @generated cost.
        for rxn in [uni_uni_rxn, bi_bi_rxn, bi_bi_pp_rxn]
            specs = EnzymeRates.init_mechanisms(rxn)
            for spec in first(specs, 5)
                m = EnzymeMechanism(spec)
                @test m isa EnzymeMechanism
                @test length(EnzymeRates.fitted_params(m)) <=
                    spec.n_fit_params_estimate
            end
        end
    end

    @testset "Drops unbound regulators from spec→type" begin
        # init_mechanisms produces specs without dead-end regulators bound.
        # When compiled to EnzymeMechanism, the regulator must NOT appear
        # in the regulators tuple — only the catalytic mechanism is built.
        # After expand_mechanisms adds the dead-end regulator, it should
        # appear.
        init_specs = EnzymeRates.init_mechanisms(uni_uni_with_reg)
        @test !isempty(init_specs)
        for spec in init_specs
            m = EnzymeRates.EnzymeMechanism(spec)
            @test :I ∉ EnzymeRates.regulators(m)
        end

        expanded = EnzymeRates.expand_mechanisms(init_specs, uni_uni_with_reg)
        found_with_reg = false
        for (_, specs) in expanded
            for spec in specs
                m = EnzymeRates.EnzymeMechanism(spec)
                if :I in EnzymeRates.regulators(m)
                    found_with_reg = true
                    break
                end
            end
            found_with_reg && break
        end
        @test found_with_reg
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 3. Base-spec expansion moves (polymorphic over Mechanism/AllostericMechanismSpec)
# ═══════════════════════════════════════════════════════════════════════

# ─── _expand_re_to_ss ──────────────────────────────────────────────────
@testset "_expand_re_to_ss" begin

    @testset "MechanismSpec — uni-uni: 2 RE binding groups → 2 variants" begin
        # SEED: uni-uni with 3 singleton kinetic groups.
        # Group 1 = E+P binding (RE), group 2 = E+S binding (RE),
        # group 3 = iso E_S↔E_P (SS).
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_rxn)

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: RE→SS fires per all-RE kinetic group atomically.
        # The seed has 2 all-RE groups (P-binding, S-binding). The iso
        # group is already SS so it's excluded. → 2 variants.
        @test length(result) == 2

        # 2. Δ params: each conversion replaces 1 RE param (K) with 2 SS
        # params (kf, kr). For a plain MechanismSpec, _re_to_ss_delta = +1
        # (ratchet of 1 K → kf + kr is a +1 net under the kinetic-group
        # accounting in _n_fit_params_estimate_from_steps).
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability — implicit in item 4's equivalence-style call.

        # 4. structural change — equivalence-style (N=2 ≤ 6).
        # Variant A: P-binding flipped to SS (group 1 RE→SS).
        # Variant B: S-binding flipped to SS (group 2 RE→SS).
        # No third variant exists because the iso group was already SS.
        v_p_flipped = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P <--> E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        v_s_flipped = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S <--> E_S
                E_S <--> E_P
            end
        end
        expected = Set([v_p_flipped, v_s_flipped])
        @test Set(EnzymeRates.compile_mechanism(r) for r in result) == expected

        # 5. preservation: reaction unchanged; non-flipped steps remain
        # in their original RE/SS state with the same kinetic_group.
        for r in result
            @test r.reaction === spec.reaction
            # The seed has only singleton kinetic groups, so flipping a
            # group means flipping exactly its one step. (For multi-step
            # groups the per-step count would equal the group's size; see
            # the bi-bi multi-step testset for atomic-flip property check.)
            n_newly_ss = count(zip(spec.steps, r.steps)) do (s_old, s_new)
                s_old.is_equilibrium && !s_new.is_equilibrium &&
                    s_old.kinetic_group == s_new.kinetic_group
            end
            @test n_newly_ss == 1
        end
    end

    @testset "MechanismSpec — bi-bi sequential: 4 RE binding groups → 4 variants" begin
        # SEED: bi-bi sequential. 2 binding groups (one for A, one for B
        # via parens to share kinetic group; same for P, Q). But here
        # we use the simplest sequential bi-bi where each metabolite has
        # its own singleton group — that gives 4 RE groups + 1 SS iso.
        # Sequential ordered: E + A → E_A + B → E_A_B ↔ E_P_Q → E + P/Q
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, bi_bi_rxn)

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: 4 all-RE singleton groups (A, B, Q, P bindings).
        # Iso group is already SS. → 4 variants, one per RE group.
        @test length(result) == 4

        # 2. Δ params: +1 per variant (plain MechanismSpec, no allosteric).
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability — implicit in item 4's equivalence-style call.

        # 4. structural change — N=4 ≤ 6 → equivalence-style.
        # Each variant flips exactly one of the four RE binding groups.
        v_a = @enzyme_mechanism begin
            substrates: A, B; products: P, Q
            steps: begin
                E + A <--> E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        v_b = @enzyme_mechanism begin
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B <--> E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        v_q = @enzyme_mechanism begin
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q <--> E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        v_p = @enzyme_mechanism begin
            substrates: A, B; products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P <--> E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        expected = Set([v_a, v_b, v_q, v_p])
        @test Set(EnzymeRates.compile_mechanism(r) for r in result) == expected

        # 5. preservation: each result has exactly one step's is_equilibrium
        # flipped from true to false, with kinetic_group preserved.
        for r in result
            @test r.reaction === spec.reaction
            n_newly_ss = count(zip(spec.steps, r.steps)) do (s_old, s_new)
                s_old.is_equilibrium && !s_new.is_equilibrium &&
                    s_old.kinetic_group == s_new.kinetic_group
            end
            @test n_newly_ss == 1
        end
    end

    @testset "MechanismSpec — bi-bi multi-step kinetic group: atomic conversion" begin
        # SEED: bi-bi random where A binds at two forms (E and E_B) and
        # those two RE binding steps share kinetic_group 1 (parenthesized).
        # B-binding shares group 2 (E and E_A). P shares group 3.
        # Q shares group 4. Iso = group 5 (SS).
        # When _expand_re_to_ss fires on group 1, BOTH A-binding steps
        # flip atomically (same group → same kinetic params).
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E_A, E_B + A ⇌ E_A_B)
                (E + B ⇌ E_B, E_A + B ⇌ E_A_B)
                (E + P ⇌ E_P, E_Q + P ⇌ E_P_Q)
                (E + Q ⇌ E_Q, E_P + Q ⇌ E_P_Q)
                E_A_B <--> E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, bi_bi_rxn)

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: 4 multi-step RE groups (A, B, P, Q each with 2 steps).
        # Iso group SS. → 4 variants, each flipping 2 steps atomically.
        @test length(result) == 4

        # 2. Δ params: +1 per variant.
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style check: in each variant, exactly one kinetic
        # group has ALL its steps now SS (atomic conversion). All other
        # groups retain their original state.
        for r in result
            groups = Dict{Int, Vector{Bool}}()
            for st in r.steps
                push!(get!(groups, st.kinetic_group, Bool[]), st.is_equilibrium)
            end
            # Exactly one group: all-false (newly SS, was multi-step RE).
            n_all_ss_multi = count(((_, vs),) ->
                length(vs) >= 2 && all(==(false), vs), groups)
            @test n_all_ss_multi == 1
        end

        # 5. preservation: reaction unchanged.
        for r in result
            @test r.reaction === spec.reaction
        end
    end

    @testset "MechanismSpec — bi-bi ping-pong: 5 RE groups → 5 variants" begin
        # SEED: bi-bi ping-pong topology with Estar (residual) form.
        # 3 singleton RE groups (A binding, Q release on E-side; B binding
        # via Estar; A→Estar iso; Estar→E iso). Iso steps are SS.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                Estar + B ⇌ Estar_B
                E + Q ⇌ E_Q
                Estar + P ⇌ Estar_A_P
                E_A <--> Estar_A_P
                Estar_B ⇌ E_Q
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, bi_bi_pp_rxn)

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: RE groups = E+A, Estar+B, E+Q, Estar+P, Estar_B → E_Q
        # (5 groups). The iso E_A↔Estar_A_P is SS so excluded. Among the
        # remaining 5, count how many are RE in the seed: E+A, Estar+B,
        # E+Q, Estar+P are RE; Estar_B → E_Q is RE; iso is SS. So 5 RE
        # groups → 5 variants. (Note: derivation depends on exact group
        # numbering in the macro — verify via the seed's compiled reactions
        # tuple if it differs.)
        @test length(result) == 5

        # 2. Δ params: +1 each (plain MechanismSpec).
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 5. preservation
        for r in result
            @test r.reaction === spec.reaction
        end
    end

    @testset "MechanismSpec — all-SS catalytic seed: empty (negative)" begin
        # When every catalytic step is already SS, _expand_re_to_ss has no
        # all-RE group to fire on → empty result.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P <--> E_P
                E + S <--> E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_rxn)
        @test isempty(EnzymeRates._expand_re_to_ss(spec))
    end

    @testset "AllostericMechanismSpec — :EqualRT group: Δ=+1" begin
        # SEED: uni-uni with all groups :EqualRT. Each catalytic group's
        # R/T tag is :EqualRT (one shared K_R = K_T). When RE→SS converts
        # an :EqualRT group, the new (kf, kr) pair is also state-shared,
        # so Δ = +1 (the EqualRT/OnlyR/OnlyT cheap-tag delta).
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: EqualRT
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: 2 all-RE groups (P-binding, S-binding); iso is SS.
        # _expand_re_to_ss fires per group; same as plain. → 2 variants.
        @test length(result) == 2

        # 2. Δ params: :EqualRT is a cheap tag → +1 per variant.
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability — must produce AllostericEnzymeMechanism.
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: exactly one group now SS; all group_tags
        # preserved including the converted group's :EqualRT tag (move
        # MUST NOT change R/T-state semantics).
        for r in result
            n_newly_ss = count(zip(spec.base.steps, r.base.steps)) do (s_old, s_new)
                s_old.is_equilibrium && !s_new.is_equilibrium
            end
            @test n_newly_ss == 1
            @test r.group_tags == spec.group_tags  # tags untouched
        end

        # 5. preservation: catalytic_n, reg sites, reg ligand tags untouched.
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
            @test r.allosteric_multiplicities == spec.allosteric_multiplicities
            @test r.reg_ligand_tags == spec.reg_ligand_tags
            @test r.base.reaction === spec.base.reaction
        end
    end

    @testset "AllostericMechanismSpec — :OnlyR group: Δ=+1" begin
        # SEED: uni-uni allosteric with one group :OnlyR (the S-binding
        # group, group 2), others :EqualRT. The S-binding step is
        # non-functional in the T-state by construction. After RE→SS
        # converts the :OnlyR group, the new (kf, kr) pair lives in the
        # R-state only; T-state contributes no kf_T/kr_T because the
        # group is :OnlyR.
        # Δ derivation: _re_to_ss_delta returns 1 when the group's tag is
        # NOT :NonequalRT. :OnlyR is a "cheap tag" → Δ = +1, same as
        # :EqualRT. (CLAUDE.md: catalytic groups cannot be :OnlyT.)
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: OnlyR
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: 2 all-RE groups (P-binding :EqualRT, S-binding :OnlyR);
        # iso group is SS so excluded. → 2 variants.
        @test length(result) == 2

        # 2. Δ params: both groups carry "cheap" tags (:EqualRT, :OnlyR).
        # Per `_re_to_ss_delta`, +1 for any non-:NonequalRT tag.
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: exactly one group flipped to SS; ALL group_tags
        # preserved including the converted group's :OnlyR (move MUST NOT
        # change R/T-state semantics).
        for r in result
            n_newly_ss = count(zip(spec.base.steps, r.base.steps)) do (s_old, s_new)
                s_old.is_equilibrium && !s_new.is_equilibrium
            end
            @test n_newly_ss == 1
            @test r.group_tags == spec.group_tags
        end

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
            @test r.allosteric_multiplicities == spec.allosteric_multiplicities
            @test r.reg_ligand_tags == spec.reg_ligand_tags
            @test r.base.reaction === spec.base.reaction
        end
    end

    @testset "AllostericMechanismSpec — :NonequalRT group: Δ=+2" begin
        # SEED: uni-uni with one :NonequalRT group, others :EqualRT.
        # When RE→SS converts the :NonequalRT group, BOTH the R-state K
        # and the T-state K_T must split into (kf, kr) and (kf_T, kr_T).
        # Δ for :NonequalRT = 2 × base = +2.
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: NonequalRT
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)

        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: 2 RE groups (group 1 = P-binding :EqualRT,
        # group 2 = S-binding :NonequalRT). → 2 variants.
        @test length(result) == 2

        # 2. Δ params: depends on tag of the converted group.
        # P-binding :EqualRT → +1; S-binding :NonequalRT → +2.
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 2]

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: exactly one RE group flipped per variant; tag
        # of the flipped group preserved.
        for r in result
            @test r.group_tags == spec.group_tags
        end

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
        end
    end

    @testset "Substrate-as-dead-end-inhibitor overlap (S used as both)" begin
        # SEED: uni-uni where S is BOTH a substrate AND a dead-end inhibitor.
        # The reaction declares dead_end_inhibitors: S, and the seed has
        # been pre-expanded to bind :S as inhibitor (giving rise to S__reg
        # binding steps). The base RE→SS move shouldn't be confused by
        # the metabolite-overlap — it operates on kinetic groups, not names.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: S
        end
        # Build via init + dead-end expansion to get the S/__reg overlap form.
        init_specs = EnzymeRates.init_mechanisms(rxn)
        @test length(init_specs) == 1   # uni-uni: 1 catalytic topology
        seed_spec = first(init_specs)
        de_specs = EnzymeRates._expand_add_dead_end_regulator(seed_spec, rxn)
        @test !isempty(de_specs)
        spec = first(de_specs)

        # Move
        result = EnzymeRates._expand_re_to_ss(spec)

        # 1. count: each all-RE kinetic group can flip. The seed after
        # add-dead-end has groups: substrate-binding (RE), product-binding
        # (RE), iso (SS), and dead-end-S__reg-binding (RE). → 3 RE groups → 3 variants.
        @test length(result) == 3

        # 2. Δ params: +1 each (plain MechanismSpec).
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style: each variant flips exactly one RE step from
        # is_equilibrium=true to false. Across the 3 variants, the
        # flipped step's kinetic_group covers 3 distinct values (one per
        # RE group: substrate-S binding, product-P binding, dead-end
        # S__reg binding). This proves the move treats the substrate-
        # name and inhibitor-dummy-name kinetic groups as independent.
        flipped_groups = Int[]
        for r in result
            flipped = [s_new.kinetic_group
                for (s_old, s_new) in zip(spec.steps, r.steps)
                if s_old.is_equilibrium && !s_new.is_equilibrium]
            @test length(flipped) == 1
            push!(flipped_groups, only(flipped))
        end
        @test length(unique(flipped_groups)) == 3

        # 5. preservation: reaction === spec.reaction.
        for r in result
            @test r.reaction === spec.reaction
        end
    end

end

# ─── _expand_split_kinetic_group ───────────────────────────────────────
@testset "_expand_split_kinetic_group" begin

    @testset "MechanismSpec — bi-bi: 4 multi-step groups → 8 splits" begin
        # SEED: bi-bi random with 4 multi-step kinetic groups (A, B, P, Q
        # each with 2 binding steps shared via parens). Iso is singleton.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E_A, E_B + A ⇌ E_A_B)
                (E + B ⇌ E_B, E_A + B ⇌ E_A_B)
                (E + P ⇌ E_P, E_Q + P ⇌ E_P_Q)
                (E + Q ⇌ E_Q, E_P + Q ⇌ E_P_Q)
                E_A_B <--> E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, bi_bi_rxn)

        result = EnzymeRates._expand_split_kinetic_group(spec)

        # 1. count: split picks one step out of a multi-step group into
        # a fresh group. 4 groups × 2 steps each → 4 × 2 = 8 split variants.
        @test length(result) == 8

        # 2. Δ params: each split adds +1 (RE plain, no allosteric).
        # The new group inherits is_equilibrium=true; one extra K param.
        for r in result
            @test r.n_fit_params_estimate ==
                spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style (N=8 > 6): the result has exactly one new
        # kinetic_group integer (max+1) and exactly one step now belongs
        # to it. The remaining group has size n_old - 1.
        for r in result
            old_groups = Set(s.kinetic_group for s in spec.steps)
            new_groups = Set(s.kinetic_group for s in r.steps)
            extra = setdiff(new_groups, old_groups)
            @test length(extra) == 1
            new_g = first(extra)
            n_in_new = count(s.kinetic_group == new_g for s in r.steps)
            @test n_in_new == 1
        end

        # 5. preservation: reaction unchanged; total step count unchanged.
        for r in result
            @test r.reaction === spec.reaction
            @test length(r.steps) == length(spec.steps)
        end
    end

    @testset "MechanismSpec — all singleton groups: empty (negative)" begin
        # SEED: uni-uni init, every kinetic group is a singleton (size 1).
        # Splitting requires a multi-step group → no eligible group → empty.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_rxn)
        @test isempty(EnzymeRates._expand_split_kinetic_group(spec))
    end

    @testset "MechanismSpec — mixed RE/SS group sizes: deltas differ" begin
        # SEED: bi-bi random where the A-binding kinetic group has been
        # converted to SS (size-2 group, both steps SS) and the B-binding
        # kinetic group is RE (size-2 group, both steps RE). The remaining
        # P-binding (size-2, RE), Q-binding (size-2, RE), and iso (singleton,
        # SS) groups are unchanged. Total: 9 steps, 5 kinetic groups.
        # Splitting an RE multi-step group adds +1 (one new K). Splitting
        # an SS multi-step group adds +2 (one new kf AND one new kr).
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A <--> E_A, E_B + A <--> E_A_B)
                (E + B ⇌ E_B, E_A + B ⇌ E_A_B)
                (E + P ⇌ E_P, E_Q + P ⇌ E_P_Q)
                (E + Q ⇌ E_Q, E_P + Q ⇌ E_P_Q)
                E_A_B <--> E_P_Q
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, bi_bi_rxn)

        result = EnzymeRates._expand_split_kinetic_group(spec)

        # 1. count: 4 multi-step groups (A SS×2, B RE×2, P RE×2, Q RE×2).
        # Each can split per member → 4 × 2 = 8 variants.
        @test length(result) == 8

        # 2. Δ params: split on an RE group → +1 (one new K). Split on an
        # SS group → +2 (one new kf + one new kr). The seed has 1 SS multi-
        # step group (A) and 3 RE multi-step groups (B, P, Q). Per-member
        # splits give 2 SS splits at +2 each and 6 RE splits at +1 each.
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 1, 1, 1, 1, 2, 2]

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style: each result introduces exactly one new
        # kinetic group (max-id + 1) with exactly one step in it.
        old_max = maximum(s.kinetic_group for s in spec.steps)
        for r in result
            new_max = maximum(s.kinetic_group for s in r.steps)
            @test new_max == old_max + 1
            @test count(s -> s.kinetic_group == new_max, r.steps) == 1
        end

        # 5. preservation
        for r in result
            @test r.reaction === spec.reaction
            @test length(r.steps) == length(spec.steps)
        end
    end

    @testset "AllostericMechanismSpec — split inherits parent's tag" begin
        # SEED: bi-bi allosteric where one multi-step group is :NonequalRT
        # and another is :EqualRT. Split must produce results where the
        # NEW group inherits the parent's tag — splitting is a parameter-
        # relaxation move and MUST NOT change R/T-state semantics.
        m_seed = @allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            site(:catalytic, 2): begin
                steps: begin
                    (E + A ⇌ E_A, E_B + A ⇌ E_A_B)        :: NonequalRT
                    (E + B ⇌ E_B, E_A + B ⇌ E_A_B)        :: EqualRT
                    E + P ⇌ E_P             :: EqualRT
                    E_P + Q ⇌ E_P_Q         :: EqualRT
                    E + Q ⇌ E_Q             :: EqualRT
                    E_Q + P ⇌ E_P_Q         :: EqualRT
                    E_A_B <--> E_P_Q        :: EqualRT
                end
            end
        end
        bi_bi_allo_rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, bi_bi_allo_rxn)

        result = EnzymeRates._expand_split_kinetic_group(spec)

        # 1. count: 2 multi-step groups × 2 members each = 4 variants.
        @test length(result) == 4

        # 2. Δ params: depends on tag and is_equilibrium of parent group.
        # :NonequalRT RE split: +2 (base 1 × 2 for NonequalRT factor).
        # :EqualRT RE split: +1.
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 2, 2]

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: the new group's tag equals its parent group's tag.
        pre_groups = Set(s.kinetic_group for s in spec.base.steps)
        for r in result
            post_groups = Set(s.kinetic_group for s in r.base.steps)
            new_groups = setdiff(post_groups, pre_groups)
            @test length(new_groups) == 1
            new_g = first(new_groups)
            # Identify parent: the only pre-group whose count dropped.
            pre_counts = Dict(g => count(s -> s.kinetic_group == g,
                                         spec.base.steps)
                              for g in pre_groups)
            post_counts = Dict(g => count(s -> s.kinetic_group == g,
                                          r.base.steps)
                               for g in pre_groups)
            old_g = only(g for g in pre_groups
                         if post_counts[g] < pre_counts[g])
            @test r.group_tags[new_g] == spec.group_tags[old_g]
            for g in pre_groups
                @test r.group_tags[g] == spec.group_tags[g]
            end
        end

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
            @test r.reg_ligand_tags == spec.reg_ligand_tags
        end
    end

    @testset "AllostericMechanismSpec — SS multi-step :NonequalRT split: Δ=+4" begin
        # SEED: bi-bi allosteric where one multi-step group is BOTH SS AND
        # :NonequalRT. _split_group_delta returns 4 for this case (factor 2
        # for SS × factor 2 for :NonequalRT R/T-state pair).
        m_seed = @allosteric_mechanism begin
            substrates: A, B
            products: P, Q
            site(:catalytic, 2): begin
                steps: begin
                    (E + A <--> E_A, E_B + A <--> E_A_B)        :: NonequalRT
                    (E + B ⇌ E_B, E_A + B ⇌ E_A_B)             :: EqualRT
                    (E + P ⇌ E_P, E_Q + P ⇌ E_P_Q)             :: EqualRT
                    (E + Q ⇌ E_Q, E_P + Q ⇌ E_P_Q)             :: EqualRT
                    E_A_B <--> E_P_Q                            :: EqualRT
                end
            end
        end
        bi_bi_allo_rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, bi_bi_allo_rxn)

        result = EnzymeRates._expand_split_kinetic_group(spec)

        # 1. count: 4 multi-step groups (A-binding SS×2 :NonequalRT,
        # B-binding RE×2 :EqualRT, P-binding RE×2 :EqualRT, Q-binding RE×2 :EqualRT).
        # 4 × 2 members = 8 variants.
        @test length(result) == 8

        # 2. Δ params:
        # - SS × :NonequalRT split: factor 2 (SS) × factor 2 (NonequalRT) = +4
        #   → 2 variants × +4
        # - RE × :EqualRT split: factor 1 × factor 1 = +1
        #   → 6 variants × +1
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 1, 1, 1, 1, 4, 4]

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. tag inheritance: split's new group inherits parent's tag.
        pre_groups = Set(s.kinetic_group for s in spec.base.steps)
        for r in result
            post_groups = Set(s.kinetic_group for s in r.base.steps)
            new_g = only(setdiff(post_groups, pre_groups))
            # find the parent group
            pre_counts = Dict(g => count(s -> s.kinetic_group == g, spec.base.steps)
                              for g in pre_groups)
            post_counts = Dict(g => count(s -> s.kinetic_group == g, r.base.steps)
                               for g in pre_groups)
            old_g = only(g for g in pre_groups
                         if post_counts[g] < pre_counts[g])
            @test r.group_tags[new_g] == spec.group_tags[old_g]
        end
    end
end

# ─── _expand_add_dead_end_regulator ────────────────────────────────────
@testset "_expand_add_dead_end_regulator" begin

    @testset "Uni-uni + I: 1 variant (equivalence-style)" begin
        # SEED: uni-uni init, no regulators bound yet.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
        end
        # Hybrid: m_seed has no regulators yet, but rxn declares :I as a
        # dead-end inhibitor. The helper allows m's regulators to be a
        # subset of rxn's, so this constructs a spec attached to the
        # full rxn directly.
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn)

        result = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)

        # 1. count: eligible forms for I-binding = forms NOT bound by all subs
        # AND NOT bound by all prods. E is the only eligible form (E_S has
        # all subs → ineligible; E_P has all prods → ineligible). Inhibitor
        # competition patterns for uni-uni: 1 (S × P × no-existing-inh = 1).
        # → 1 form set → 1 variant.
        @test length(result) == 1

        # 2. Δ params: +1 (one new K_I parameter for the dead-end binding group).
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability — implicit in item 4's equivalence-style call.

        # 4. equivalence-style (N=1).
        # Expected: same uni-uni catalytic + a new I-binding RE step.
        # (The spec's step uses the :I__reg dummy form name; compile_mechanism
        # strips the __reg suffix, so the compiled mechanism has bare :I and
        # :E_I — which is what the @enzyme_mechanism literal below produces.)
        expected = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: I
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
                E + I ⇌ E_I
            end
        end
        @test EnzymeRates.compile_mechanism(first(result)) === expected

        # 5. preservation
        @test first(result).reaction === spec.reaction
    end

    @testset "Uni-uni no regulators → empty (negative)" begin
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_rxn)
        @test isempty(EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_rxn))
    end

    @testset "exclude_regs kwarg suppresses regulator addition (negative)" begin
        # SEED: uni-uni with two dead-end inhibitors I and J available.
        # Compares baseline (no kwarg) to filtered (exclude_regs=Set([:I]))
        # and verifies :I appears in baseline results but is absent when
        # excluded. Passing exclude_regs=Set([:I,:J]) yields empty.
        rxn_ij = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I, J
        end
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn_ij)

        # Without exclude_regs: 2 eligible regs × 1 form pattern → 2 variants.
        # Property: both I and J appear across the result set.
        baseline = EnzymeRates._expand_add_dead_end_regulator(spec, rxn_ij)
        @test length(baseline) == 2
        has_i_baseline = any(any(contains(string(sym), "I__reg")
                                 for st in s.steps
                                 for sym in Iterators.flatten((st.reactants, st.products)))
                             for s in baseline)
        has_j_baseline = any(any(contains(string(sym), "J__reg")
                                 for st in s.steps
                                 for sym in Iterators.flatten((st.reactants, st.products)))
                             for s in baseline)
        @test has_i_baseline && has_j_baseline

        # With exclude_regs=Set([:I]): only J is eligible → 1 variant, only J in result.
        excluded = EnzymeRates._expand_add_dead_end_regulator(spec, rxn_ij; exclude_regs=Set([:I]))
        @test length(excluded) == 1
        has_i_excluded = any(any(contains(string(sym), "I__reg")
                                 for st in s.steps
                                 for sym in Iterators.flatten((st.reactants, st.products)))
                             for s in excluded)
        @test !has_i_excluded
        has_j_excluded = any(any(contains(string(sym), "J__reg")
                                 for st in s.steps
                                 for sym in Iterators.flatten((st.reactants, st.products)))
                             for s in excluded)
        @test has_j_excluded

        # With exclude_regs=Set([:I, :J]): no eligible regs → empty.
        @test isempty(EnzymeRates._expand_add_dead_end_regulator(
            spec, rxn_ij; exclude_regs=Set([:I, :J])))
    end

    @testset "Sequential bi-bi + I: 4 distinct form sets" begin
        # SEED: bi-bi sequential. Eligible forms: E, E_A, E_Q.
        # Inhibitor competition patterns enumerate which subs/prods I
        # competes with; combined with which forms have a binding step
        # for those mets, produces exactly 4 distinct form sets:
        #   {E, E_Q}: ({A},{P}), ({A,B},{P})
        #   {E, E_A}: ({B},{Q}), ({B},{P,Q})
        #   {E_A, E_Q}: ({B},{P})
        #   {E}:      ({A},{Q}), ({A},{P,Q}), ({A,B},{Q}), ({A,B},{P,Q})
        # → 4 unique form sets after dedup.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn)
        result = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)

        # 1. count
        @test length(result) == 4

        # 2. Δ params
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style — each variant adds at least one I__reg binding
        # step; all I__reg binding steps in a single variant share the same
        # new kinetic_group.
        for r in result
            i_binding_groups = unique(
                s.kinetic_group for s in r.steps
                if length(s.reactants) == 2 &&
                    contains(string(s.reactants[2]), "I__reg"))
            @test length(i_binding_groups) == 1
        end

        # 5. preservation
        for r in result
            @test r.reaction === spec.reaction
        end
    end

    @testset "Bi-bi random + I: 9 variants (property-style)" begin
        # SEED: bi-bi random. Eligible forms = E, E_A, E_B, E_P, E_Q (E_A_B
        # has all subs → ineligible; E_P_Q has all prods → ineligible).
        # Inhibitor competition patterns: 9 (3 sub subsets × 3 prod subsets,
        # no existing inhibitors). Each pattern produces a distinct active-
        # form set after dedup → exactly 9 variants.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E_A, E_B + A ⇌ E_A_B)
                (E + B ⇌ E_B, E_A + B ⇌ E_A_B)
                (E + P ⇌ E_P, E_Q + P ⇌ E_P_Q)
                (E + Q ⇌ E_Q, E_P + Q ⇌ E_P_Q)
                E_A_B <--> E_P_Q
            end
        end
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn)
        result = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)

        # 1. count
        @test length(result) == 9

        # 2. Δ params
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style: each variant has at least one I__reg binding
        # step; all I__reg binding steps in a single variant share the same
        # kinetic_group (one K_I, not multiple).
        for r in result
            i_binding_groups = unique(
                s.kinetic_group for s in r.steps
                if length(s.reactants) == 2 &&
                    contains(string(s.reactants[2]), "I__reg"))
            @test length(i_binding_groups) == 1
        end

        # 5. preservation
        for r in result
            @test r.reaction === spec.reaction
        end
    end

    @testset "Bi-bi PP + I: 3 variants" begin
        # SEED: bi-bi ping-pong. Forms: E, E_A, Estar, Estar_A_P, Estar_B,
        # E_Q. Estar_A_P has all prods → ineligible. Estar_B has all subs
        # → ineligible. Eligible: E, E_A, Estar, E_Q (4 forms). Inhibitor
        # competition patterns × dedup → 3 unique form sets.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                Estar + B ⇌ Estar_B
                E + Q ⇌ E_Q
                Estar + P ⇌ Estar_A_P
                E_A <--> Estar_A_P
                Estar_B ⇌ E_Q
            end
        end
        rxn = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn)
        result = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)

        # 1. count
        @test length(result) == 3

        # 2. Δ params
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style (N=3, but mechanisms differ in atom-bracket
        # placement — equivalence at the EnzymeMechanism level is
        # cumbersome). Each variant adds one new kinetic_group used only
        # by the I__reg binding steps.
        for r in result
            i_binding_groups = unique(
                s.kinetic_group for s in r.steps
                if length(s.reactants) == 2 &&
                    contains(string(s.reactants[2]), "I__reg"))
            @test length(i_binding_groups) == 1
        end

        # 5. preservation
        for r in result
            @test r.reaction === spec.reaction
        end
    end

    @testset "Two regulators chain: J__reg dummy naming has no numeric suffix" begin
        # SEED: uni-uni with rxn declaring two dead-end inhibitors I, J.
        # Step A: add I via the move.
        # Step B: pick a variant with I bound, call the move again to add J.
        # Property: J's dummy symbol is exactly :J__reg (no numeric suffix
        # like :J__reg2). Catches the regulator-naming-stability regression.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I, J
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn)

        # Step A: 2 eligible regs (I, J), 1 form each → 2 variants total.
        i_or_j_specs = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)
        with_i = first(filter(i_or_j_specs) do s
            any(contains(string(sym), "I__reg")
                for st in s.steps
                for sym in Iterators.flatten(
                    (st.reactants, st.products)))
        end)

        # Step B
        j_specs = EnzymeRates._expand_add_dead_end_regulator(with_i, rxn)
        @test !isempty(j_specs)

        # Property: J__reg with no numeric suffix.
        for s in j_specs
            j_syms = [sym for st in s.steps
                      for sym in Iterators.flatten(
                          (st.reactants, st.products))
                      if contains(string(sym), "J__reg")]
            @test !isempty(j_syms)
            for sym in j_syms
                @test !occursin(r"__reg\d", string(sym))
            end

            # Adding J must NOT remove I-binding steps from the spec.
            i_syms = [sym for st in s.steps
                      for sym in Iterators.flatten(
                          (st.reactants, st.products))
                      if contains(string(sym), "I__reg")]
            @test !isempty(i_syms)
        end
    end

    @testset "Two regulators competition: 17 variants (property-style)" begin
        # SEED: bi-bi random with rxn declaring I1, I2 as dead-end inhibitors.
        # Step A: add both regs in one call → 9 + 9 = 18 variants.
        # Step B: pick variant where I1 binds at multiple forms; call the
        # move again to add I2. Inhibitor competition patterns with one
        # existing inhibitor: 9 × 2 = 18 raw patterns, but I2's active-form
        # set deduplication after intersection with eligible forms yields
        # 17 unique form sets → 17 variants.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E_A, E_B + A ⇌ E_A_B)
                (E + B ⇌ E_B, E_A + B ⇌ E_A_B)
                (E + P ⇌ E_P, E_Q + P ⇌ E_P_Q)
                (E + Q ⇌ E_Q, E_P + Q ⇌ E_P_Q)
                E_A_B <--> E_P_Q
            end
        end
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I1, I2
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn)

        # Step A
        result1 = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)
        @test length(result1) == 18

        # Pick variant where I1 binds at multiple forms.
        multi = filter(result1) do r
            inh_forms = filter(
                f -> contains(string(f), "I1__reg"),
                collect(EnzymeRates.all_form_names(r)))
            length(inh_forms) >= 2
        end
        @test !isempty(multi)
        spec_i1 = first(multi)

        # Step B: add I2 to a spec that already has I1.
        result2 = EnzymeRates._expand_add_dead_end_regulator(spec_i1, rxn)
        @test length(result2) == 17

        # Property: at least one variant has I1 + I2 coexisting on the same
        # form (non-competing pattern), and at least one has I2 forms that
        # never coexist with I1 (fully-competing pattern).
        has_coexist = any(result2) do r
            any(
                f -> contains(string(f), "I1__reg") &&
                     contains(string(f), "I2__reg"),
                collect(EnzymeRates.all_form_names(r)))
        end
        @test has_coexist
        has_compete = any(result2) do r
            forms = collect(EnzymeRates.all_form_names(r))
            has_i2 = any(
                f -> contains(string(f), "I2__reg"),
                forms)
            no_coexist = !any(
                f -> contains(string(f), "I1__reg") &&
                     contains(string(f), "I2__reg"),
                forms)
            has_i2 && no_coexist
        end
        @test has_compete
    end

    @testset "Substrate-as-dead-end-inhibitor overlap (S used as both)" begin
        # SEED: uni-uni where :S is both a substrate AND declared as a
        # dead-end inhibitor. The move treats the substrate-:S binding step
        # and the inhibitor-:S__reg binding step as independent — the dummy
        # :S__reg is what enters the new dead-end binding step.
        rxn_overlap = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: S
        end
        # Init from rxn_overlap (1 catalytic topology for uni-uni).
        init_specs = EnzymeRates.init_mechanisms(rxn_overlap)
        @test length(init_specs) == 1
        spec = first(init_specs)

        result = EnzymeRates._expand_add_dead_end_regulator(spec, rxn_overlap)

        # 1. count: eligible forms = E (the form not bound by all subs and
        # not bound by all prods). Inhibitor competition patterns for uni-
        # uni with no existing inhibitors: 1. → 1 variant.
        @test length(result) == 1

        # 2. Δ params
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa EnzymeMechanism
        end

        # 4. property-style: variant introduces an :S__reg-named binding
        # step (proves the substrate-:S vs inhibitor-:S__reg distinction
        # is preserved by the move's dummy-symbol rewriting).
        for r in result
            reg_syms = [sym for st in r.steps
                        for sym in Iterators.flatten(
                            (st.reactants, st.products))
                        if contains(string(sym), "S__reg")]
            @test !isempty(reg_syms)
            # Exactly one new kinetic group introduced (the S__reg binding).
            new_groups = setdiff(
                Set(s.kinetic_group for s in r.steps),
                Set(s.kinetic_group for s in spec.steps))
            @test length(new_groups) == 1
        end

        # 5. preservation
        for r in result
            @test r.reaction === spec.reaction
        end
    end

    @testset "AllostericMechanismSpec input: dead-end binding tagged :EqualRT" begin
        # SEED: uni-uni allosteric (catalytic_n=2) with mixed regs in rxn:
        # :I dead-end, :R allosteric. The move must (a) exclude :R from
        # eligible_regs (allosteric ligand), (b) tag the new dead-end
        # binding kinetic group :EqualRT.
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            allosteric_regulators: R
            oligomeric_state: 2
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)

        result = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)

        # 1. count: same as plain uni-uni + I (1 eligible form E),
        # because :R is excluded as allosteric.
        @test length(result) == 1

        # 2. Δ params: +1 (EqualRT gives a single shared K parameter).
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability — assert via item 4 plus type check below.
        for r in result
            _assert_spec_invariants(r)
            @test r isa AllostericMechanismSpec
            @test EnzymeRates.compile_mechanism(r) isa
                AllostericEnzymeMechanism
        end

        # 4. property-style: exactly one new kinetic group, tagged :EqualRT.
        for r in result
            new_groups = setdiff(
                Set(s.kinetic_group for s in r.base.steps),
                Set(s.kinetic_group for s in spec.base.steps))
            @test length(new_groups) == 1
            new_g = first(new_groups)
            @test r.group_tags[new_g] == :EqualRT
            # Pre-existing group tags must be unchanged.
            for (g, t) in spec.group_tags
                @test r.group_tags[g] == t
            end
        end

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
            @test r.reg_ligand_tags == spec.reg_ligand_tags
        end
    end

    @testset "Allosteric-only regulator → empty (negative)" begin
        # SEED: uni-uni allosteric with rxn declaring only :R as an
        # allosteric regulator (no dead-end inhibitors). All declared
        # regulators are allosteric ligands → eligible_regs is empty
        # → result is empty.
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo_reg)
        @test isempty(EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_allo_reg))
    end

    @testset "Plain MechanismSpec: I__reg bindings share one kinetic group" begin
        # SEED: bi-bi random + I. Variants where I binds at multiple forms
        # MUST keep all I__reg binding steps in a single shared kinetic
        # group (one K_I parameter, not one per form). This is invariant
        # across all variants.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ E_A, E_B + A ⇌ E_A_B)
                (E + B ⇌ E_B, E_A + B ⇌ E_A_B)
                (E + P ⇌ E_P, E_Q + P ⇌ E_P_Q)
                (E + Q ⇌ E_Q, E_P + Q ⇌ E_P_Q)
                E_A_B <--> E_P_Q
            end
        end
        rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn)
        result = EnzymeRates._expand_add_dead_end_regulator(spec, rxn)

        # Pick variants with at least 2 I__reg binding steps (multi-form
        # inhibitor binding).
        multi = filter(result) do r
            count(r.steps) do s
                length(s.reactants) == 2 &&
                    contains(string(s.reactants[2]), "I__reg")
            end >= 2
        end
        @test !isempty(multi)
        for r in multi
            i_binding_groups = unique(
                s.kinetic_group for s in r.steps
                if length(s.reactants) == 2 &&
                    contains(string(s.reactants[2]), "I__reg"))
            @test length(i_binding_groups) == 1
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 4. Allosteric expansion moves (AllostericMechanismSpec only)
# ═══════════════════════════════════════════════════════════════════════

# ─── _expand_to_allosteric ─────────────────────────────────────────────
@testset "_expand_to_allosteric" begin

    @testset "MechanismSpec — uni-uni: 4 variants (equivalence-style)" begin
        # SEED: uni-uni init, 3 singleton kinetic groups.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)

        result = EnzymeRates._expand_to_allosteric(spec, uni_uni_allo)

        # 1. count: _expand_to_allosteric emits the all-:EqualRT baseline
        # once plus one :OnlyR variant per kinetic group. 3 groups → 1 + 3 = 4.
        @test length(result) == 4

        # 2. Δ params: +1 per variant (just L, the conformation equilibrium).
        # All other tag deltas are zero relative to the all-:EqualRT baseline.
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability — implicit in item 4's equivalence-style call.

        # 4. equivalence-style (N=4 ≤ 6). 4 expected mechanisms:
        # all-:EqualRT baseline + one :OnlyR variant per group.
        v_baseline = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_g1_OnlyR = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: OnlyR
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_g2_OnlyR = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: OnlyR
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_g3_OnlyR = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: OnlyR
                end
            end
        end
        @test Set(EnzymeRates.compile_mechanism(r) for r in result) ==
            Set([v_baseline, v_g1_OnlyR, v_g2_OnlyR, v_g3_OnlyR])

        # 5. preservation: base spec's reaction and step content unchanged
        # (the move only attaches allosteric tags; the catalytic mechanism is
        # carried through unmodified).
        for r in result
            @test r.base.reaction === spec.reaction
            @test r.base.steps == spec.steps
            @test r.base.n_fit_params_estimate == spec.n_fit_params_estimate
            @test r.catalytic_n == 2
        end
    end

    @testset "AllostericMechanismSpec → empty (negative)" begin
        # Already-allosteric specs cannot be re-converted. The
        # _expand_to_allosteric specialization on AllostericMechanismSpec
        # returns an empty vector.
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)
        result_via_dispatch = EnzymeRates._expand_to_allosteric(spec, uni_uni_allo)
        @test isempty(result_via_dispatch)
        @test result_via_dispatch isa Vector{AllostericMechanismSpec}
    end

    @testset "oligomeric_state from reaction" begin
        # The catalytic_n of the result is taken from the reaction's
        # oligomeric_state, not hardcoded to 2.
        rxn4 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 4
        end
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, rxn4)
        result = EnzymeRates._expand_to_allosteric(spec, rxn4)
        @test !isempty(result)
        for r in result
            @test r.catalytic_n == 4
        end
    end

    @testset "Bi-bi sequential: 5 groups → 6 variants" begin
        # _expand_to_allosteric emits 1 baseline + 1 :OnlyR per group.
        # 5 groups → 1 + 5 = 6 variants.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                E_A + B ⇌ E_A_B
                E + Q ⇌ E_Q
                E_Q + P ⇌ E_P_Q
                E_A_B <--> E_P_Q
            end
        end
        bi_bi_allo_rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, bi_bi_allo_rxn)
        result = EnzymeRates._expand_to_allosteric(spec, bi_bi_allo_rxn)
        @test length(result) == 6
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end
    end

    @testset "Bi-bi ping-pong: 6 groups → 7 variants" begin
        # SEED: bi-bi ping-pong topology mapped to allosteric.
        # 6 kinetic groups (5 RE binding + 1 iso). Plus the second iso step
        # makes 7 groups total. Move emits 1 baseline + 6 :OnlyR variants.
        m_seed = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E_A
                Estar + B ⇌ Estar_B
                E + Q ⇌ E_Q
                Estar + P ⇌ Estar_A_P
                E_A <--> Estar_A_P
                Estar_B ⇌ E_Q
            end
        end
        bi_bi_pp_allo_rxn = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            oligomeric_state: 2
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, bi_bi_pp_allo_rxn)
        result = EnzymeRates._expand_to_allosteric(spec, bi_bi_pp_allo_rxn)

        # 1. count: 6 kinetic groups (5 binding + 1 iso) → 7 variants
        # (1 baseline + 6 :OnlyR per group). Verify the actual group count
        # from the seed before accepting this number.
        n_groups = length(unique(s.kinetic_group for s in spec.steps))
        @test length(result) == n_groups + 1

        # 2. Δ params: +1 (just L) per variant.
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end
    end
end

# ─── _expand_add_allosteric_regulator ──────────────────────────────────
@testset "_expand_add_allosteric_regulator" begin

    @testset "Allosteric uni-uni + first allo regulator R: 3 variants" begin
        # SEED: uni-uni allosteric, all groups :EqualRT, no allosteric
        # regulator added yet.
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo_reg)

        result = EnzymeRates._expand_add_allosteric_regulator(
            spec, uni_uni_allo_reg)

        # 1. count: R is the only un-added allosteric regulator. There
        # are 0 existing reg sites. The move enumerates non-:EqualRT
        # tags {:OnlyR, :OnlyT, :NonequalRT} × site options (new site
        # only, since 0 existing): 3 × 1 = 3 variants. The :EqualRT
        # branch is gated to "existing site with at least one non-:EqualRT
        # ligand" → not applicable here (no existing sites). → 3.
        @test length(result) == 3

        # 2. Δ params: cost of new R-binding K (+1) plus per-tag delta vs
        # :EqualRT base. :OnlyR/:OnlyT cheap → +1 total. :NonequalRT → +2
        # (K_R + K_T). So deltas across the 3 variants: [1, 1, 2].
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 2]

        # 3. compilability — implicit in item 4's equivalence-style call.

        # 4. equivalence-style (N=3 ≤ 6).
        v_onlyR = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::OnlyR
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_onlyT = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::OnlyT
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        v_neq = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::NonequalRT
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        @test Set(EnzymeRates.compile_mechanism(r) for r in result) ==
            Set([v_onlyR, v_onlyT, v_neq])

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.group_tags == spec.group_tags
            @test r.base.reaction === spec.base.reaction
        end
    end

    @testset "Non-allosteric MechanismSpec → empty (negative)" begin
        # MechanismSpec is not an allosteric spec. The move specializes on
        # AllostericMechanismSpec; a plain MechanismSpec dispatches to the
        # fallback that returns an empty vector.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo_reg)

        result = EnzymeRates._expand_add_allosteric_regulator(
            spec, uni_uni_allo_reg)

        # 1. count: 0 — non-allosteric input → empty fallback.
        @test isempty(result)
        @test result isa Vector{AllostericMechanismSpec}
    end

    @testset "Two regulators with site options: count = 7" begin
        # SEED: allosteric uni-uni with R1 already added as :OnlyR.
        # R2 is un-added; existing site has one non-:EqualRT ligand (R1).
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R1::OnlyR
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo_2reg)

        result = EnzymeRates._expand_add_allosteric_regulator(
            spec, uni_uni_allo_2reg)

        # 1. count: 3 non-:EqualRT tag flavors {:OnlyR, :OnlyT, :NonequalRT}
        # × 2 site options (new site OR R1's existing site) = 6.
        # Plus 1 variant for :EqualRT at R1's (non-:EqualRT) existing
        # site — gated on existing site having a non-:EqualRT ligand (R1
        # qualifies). → 6 + 1 = 7.
        @test length(result) == 7

        # 2. Δ params: derivation from _expand_add_allosteric_regulator
        # source. delta_cost = _allo_lig_state_delta(:EqualRT, tag) + 1.
        # Non-:EqualRT branch (3 tags × 2 sites = 6 variants):
        #   :OnlyR (cost 1): (1-1) + 1 = +1 → 2 variants × +1
        #   :OnlyT (cost 1): (1-1) + 1 = +1 → 2 variants × +1
        #   :NonequalRT (cost 2): (2-1) + 1 = +2 → 2 variants × +2
        # :EqualRT-at-existing branch (1 variant; gated on R1 ≠ :EqualRT):
        #   :EqualRT (cost 1): (1-1) + 1 = +1 → 1 variant × +1
        # Sorted multiset: 5 ones + 2 twos = [1, 1, 1, 1, 1, 2, 2].
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 1, 1, 1, 2, 2]

        # 4. structural: every result has R2 in allosteric_reg_sites.
        for r in result
            has_r2 = any(:R2 in site for site in r.allosteric_reg_sites)
            @test has_r2
        end

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.group_tags == spec.group_tags
            @test r.base.reaction === spec.base.reaction
        end
    end

    @testset "EqualRT ligand reachable at existing reg site" begin
        # SEED: same as seed 3 — R1 is :OnlyR at site 1. Adding R2 must
        # produce at least one variant where R2 is :EqualRT at site 1
        # (same site as R1). The :EqualRT branch fires only when the
        # existing site has a non-:EqualRT ligand.
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R1::OnlyR
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo_2reg)

        result = EnzymeRates._expand_add_allosteric_regulator(
            spec, uni_uni_allo_2reg)

        # 1. count: 7 (same derivation as seed 3).
        @test length(result) == 7

        # 4. structural: at least one result has R2::EqualRT at site 1
        # (the same site where R1 lives).
        target = findfirst(result) do s
            get(s.reg_ligand_tags, :R2, nothing) == :EqualRT &&
                :R2 in s.allosteric_reg_sites[1]
        end
        @test target !== nothing

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.group_tags == spec.group_tags
            @test r.base.reaction === spec.base.reaction
        end
    end

    @testset "Adding :EqualRT R2 at site with :OnlyT R1" begin
        # SEED: allosteric uni-uni with R1::OnlyT already present at site 1.
        # Adding R2 should enumerate non-:EqualRT tags (×2 sites: new + existing) +
        # :EqualRT at existing site (because R1 is non-:EqualRT, the
        # :EqualRT-at-existing branch fires). Total: 3×2 + 1 = 7 variants.
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R1::OnlyT
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: EqualRT
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo_2reg)

        result = EnzymeRates._expand_add_allosteric_regulator(spec, uni_uni_allo_2reg)

        # 1. count: 3 non-:EqualRT tags × 2 site options + 1 :EqualRT-at-existing
        # = 7 variants.
        @test length(result) == 7

        # 2. Δ params: same multiset as the analogous :OnlyR seed
        # (5 ones + 2 twos = [1,1,1,1,1,2,2]).
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 1, 1, 1, 2, 2]

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: at least one variant has R2 :EqualRT at site 1
        # (the EqualRT-at-existing branch, gated on R1 being non-:EqualRT).
        has_eq_at_site1 = any(result) do r
            :R2 in r.allosteric_reg_sites[1] && r.reg_ligand_tags[:R2] == :EqualRT
        end
        @test has_eq_at_site1
    end

    @testset "Substrate-as-allosteric-regulator overlap" begin
        # SEED: allosteric uni-uni where S is both substrate and allosteric
        # regulator. Reaction declares S in allosteric_regulators.
        rxn_allo_overlap = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: S
            oligomeric_state: 2
        end
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn_allo_overlap)

        result = EnzymeRates._expand_add_allosteric_regulator(
            spec, rxn_allo_overlap)

        # 1. count: S is the only un-added allosteric regulator. 0 existing
        # reg sites → 3 non-:EqualRT tags × 1 site option = 3. No :EqualRT
        # branch (no existing sites). → 3.
        @test length(result) == 3

        # 2. Δ params: same as seed 1 — [1, 1, 2].
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 2]

        # 4. structural: S appears in allosteric_reg_sites of every result.
        # S retains its catalytic role in the base spec (substrates of rxn).
        for r in result
            has_s = any(:S in site for site in r.allosteric_reg_sites)
            @test has_s
        end

        # 4b. compilability (explicit since dual-role is unusual).
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 5. preservation
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.group_tags == spec.group_tags
            @test r.base.reaction === spec.base.reaction
        end
    end

    @testset "Two regulators at different sites" begin
        # SEED: allosteric uni-uni with R1 already present at site 1; we add R2.
        # R2 can go to: a NEW site (site_idx=0) OR R1's existing site (site_idx=1).
        # The cross-site placement (site_idx ≥ 2 in the source) requires 2+
        # existing sites; this seed has only 1 existing, so we focus on the
        # site_idx=0 vs site_idx=1 distinction.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R1, R2
            oligomeric_state: 2
        end
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R1::OnlyR
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: EqualRT
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)

        result = EnzymeRates._expand_add_allosteric_regulator(spec, rxn)

        # 1. count: 3 non-:EqualRT tags × 2 site options + 1 :EqualRT-at-existing
        # = 7 variants (same as the existing "Two regulators with site options"
        # but verifying via two-different-sites placement specifically).
        @test length(result) == 7

        # 4. property-style: separate the new-site and existing-site placements.
        new_site_variants = filter(r -> length(r.allosteric_reg_sites) == 2, result)
        existing_site_variants = filter(r -> length(r.allosteric_reg_sites) == 1, result)
        @test length(new_site_variants) == 3   # 3 tags × new site
        @test length(existing_site_variants) == 4  # 3 tags + 1 :EqualRT
    end

    @testset "Product-as-allosteric-regulator overlap" begin
        # SEED: uni-uni allosteric where product P is ALSO declared as an
        # allosteric regulator. Adding :P as allo regulator should produce
        # 3 tag variants × 1 site option = 3 variants. Verifies the move
        # treats name-overlapping ligand and product as independent.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: P
            oligomeric_state: 2
        end
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: EqualRT
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)

        result = EnzymeRates._expand_add_allosteric_regulator(spec, rxn)

        # 1. count: 3 non-:EqualRT tags × 1 new site = 3 variants.
        @test length(result) == 3

        # 2. Δ params: [1, 1, 2] (sorted: 2 cheap + 1 :NonequalRT).
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 2]

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property: :P appears in allosteric_reg_sites
        for r in result
            @test any(:P in site for site in r.allosteric_reg_sites)
        end
    end

    @testset "All declared regs already present → empty (negative)" begin
        # SEED: allosteric uni-uni with R already added. The reaction declares
        # only R as a regulator. eligible_regs computation excludes R because
        # it's in existing_allo → result is empty.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R::OnlyR
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: EqualRT
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, rxn)
        @test isempty(EnzymeRates._expand_add_allosteric_regulator(spec, rxn))
    end

end

# ─── _expand_change_allo_state ─────────────────────────────────────────
@testset "_expand_change_allo_state" begin

    @testset "Allosteric uni-uni all-:EqualRT: 3 group-tag relaxations" begin
        # SEED: uni-uni allosteric with all 3 groups tagged :EqualRT.
        # Each non-:NonequalRT entry contributes ONE relaxation variant
        # (flip its value to :NonequalRT in the dense Dict).
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)

        result = EnzymeRates._expand_change_allo_state(spec, uni_uni_allo)

        # 1. count: 3 group_tags entries non-:NonequalRT (one per group,
        # all :EqualRT) + 0 reg_ligand_tags entries (no regulators) →
        # 3 relaxation variants.
        @test length(result) == 3

        # 2. Δ params: each removal converts :EqualRT → :NonequalRT.
        # _allo_state_delta(:EqualRT, :NonequalRT, is_re):
        #   For RE (is_re=true): factor 1 × (cost(NonequalRT) - cost(EqualRT))
        #   = 1 × (2 - 1) = +1.
        # For SS iso group (is_re=false): factor 2 × (2 - 1) = +2.
        # Two RE binding groups → +1 each (2 variants). One SS iso group →
        # +2 (1 variant). Deltas: [1, 1, 2].
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 2]

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: in each result, exactly one group's tag flipped
        # from non-:NonequalRT to :NonequalRT (relaxation move).
        for r in result
            relaxed = [g for g in keys(spec.group_tags)
                       if spec.group_tags[g] != :NonequalRT &&
                          r.group_tags[g] == :NonequalRT]
            @test length(relaxed) == 1
        end

        # 5. preservation: catalytic_n, reg sites, base.reaction unchanged.
        for r in result
            @test r.catalytic_n == spec.catalytic_n
            @test r.allosteric_reg_sites == spec.allosteric_reg_sites
            @test r.base.reaction === spec.base.reaction
        end
    end

    @testset "Fully relaxed → empty (negative)" begin
        # SEED: every group_tag and reg_ligand_tag is :NonequalRT.
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: NonequalRT
                    E + S ⇌ E_S    :: NonequalRT
                    E_S <--> E_P   :: NonequalRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)
        # All group_tags == :NonequalRT, no reg_ligand_tags → no eligible
        # entries → empty result.
        @test isempty(EnzymeRates._expand_change_allo_state(spec, uni_uni_allo))
    end

    @testset "MechanismSpec → empty" begin
        # Plain MechanismSpec (no allosteric conversion) dispatches to the
        # MechanismSpec specialization which returns empty.
        spec = first(EnzymeRates.init_mechanisms(uni_uni_allo))
        result = EnzymeRates._expand_change_allo_state(spec, uni_uni_allo)
        @test isempty(result)
        @test result isa Vector{AllostericMechanismSpec}
    end

    @testset "Allosteric regulator tag removal delta" begin
        # SEED: uni-uni allosteric with one regulator R tagged :OnlyR.
        # Move yields 3 group-tag relaxations + 1 reg-ligand-tag relaxation
        # = 4 variants total.
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::OnlyR
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo_reg)

        result = EnzymeRates._expand_change_allo_state(spec, uni_uni_allo_reg)

        # 1. count: 3 group-tag relaxations + 1 reg-ligand-tag relaxation = 4.
        @test length(result) == 4

        # 2. filter for the reg-ligand-removal variant (R flips to :NonequalRT).
        r_removal = filter(
            r -> get(r.reg_ligand_tags, :R, :NonequalRT) == :NonequalRT &&
                 r.group_tags == spec.group_tags,
            result)
        @test length(r_removal) == 1

        # 3. delta for :OnlyR → :NonequalRT: cost(:NonequalRT) - cost(:OnlyR)
        # = 2 - 1 = +1.
        @test only(r_removal).n_fit_params_estimate ==
              spec.n_fit_params_estimate + 1
    end

    @testset ":OnlyT regulator-ligand relaxation" begin
        # SEED: uni-uni allosteric with one regulator R tagged :OnlyT.
        # _expand_change_allo_state should produce variants for each
        # non-:NonequalRT entry, including the :OnlyT ligand. Δ for the
        # ligand-relaxation variant: cost(:NonequalRT) - cost(:OnlyT) = +1.
        m_seed = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: R::OnlyT
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P       :: EqualRT
                    E + S ⇌ E_S       :: EqualRT
                    E_S <--> E_P      :: EqualRT
                end
            end
        end
        spec = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo_reg)

        result = EnzymeRates._expand_change_allo_state(spec, uni_uni_allo_reg)

        # 1. count: 3 group_tags entries (all :EqualRT) + 1 reg_ligand_tags
        # entry (:OnlyT) → 4 variants.
        @test length(result) == 4

        # 2. Δ params: 2 RE-binding-group EqualRT relaxations (+1 each),
        # 1 SS-iso EqualRT relaxation (+2), 1 reg-ligand :OnlyT relaxation
        # (cost(:NonequalRT) - cost(:OnlyT) = 2 - 1 = +1). Sorted: [1, 1, 1, 2].
        deltas = sort([r.n_fit_params_estimate -
                       spec.n_fit_params_estimate for r in result])
        @test deltas == [1, 1, 1, 2]

        # 3. compilability
        for r in result
            @test EnzymeRates.compile_mechanism(r) isa AllostericEnzymeMechanism
        end

        # 4. property-style: exactly one ligand-relaxation variant has the
        # R tag flipped to :NonequalRT.
        n_r_relaxed = count(r -> r.reg_ligand_tags[:R] == :NonequalRT, result)
        @test n_r_relaxed == 1
    end

end

# ═══════════════════════════════════════════════════════════════════════
# 5. Composition (dedup!, expand_mechanisms)
# ═══════════════════════════════════════════════════════════════════════

# ─── dedup! ────────────────────────────────────────────────────────────
@testset "Dedup" begin
    @testset "Same mechanism, different step order" begin
        spec1 = MechanismSpec(
            uni_uni_rxn,
            [StepSpec([:E, :S], [:E_S], true, 1),
             StepSpec([:E, :P], [:E_P], true, 2),
             StepSpec([:E_S], [:E_P], false, 3)],
            5)
        spec2 = MechanismSpec(
            uni_uni_rxn,
            [StepSpec([:E, :P], [:E_P], true, 2),
             StepSpec([:E_S], [:E_P], false, 3),
             StepSpec([:E, :S], [:E_S], true, 1)],
            5)
        cache = Dict(5 => AbstractMechanismSpec[spec1, spec2])
        EnzymeRates.dedup!(cache)
        @test length(cache[5]) == 1
    end

    @testset "Different mechanisms preserved" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        pc = first(specs).n_fit_params_estimate
        cache = Dict(pc => AbstractMechanismSpec[specs...])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) >= 1
        @test length(cache[pc]) <= length(specs)
    end

    @testset "Idempotent" begin
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        pc = first(specs).n_fit_params_estimate
        cache = Dict(pc => AbstractMechanismSpec[specs...])
        EnzymeRates.dedup!(cache)
        n1 = length(cache[pc])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) == n1
    end

    @testset "Allosteric dedup: site order" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        base = first(specs)
        used_groups = sort!(collect(
            Set(s.kinetic_group for s in base.steps)))
        group_tags = Dict{Int, Symbol}(
            g => :NonequalRT for g in used_groups)
        lig_tags = Dict{Symbol, Symbol}(
            :A => :NonequalRT, :B => :NonequalRT)
        spec_ab = AllostericMechanismSpec(
            base, 2,
            [[:A], [:B]], [2, 2],
            copy(group_tags), copy(lig_tags),
            base.n_fit_params_estimate + 2)
        spec_ba = AllostericMechanismSpec(
            base, 2,
            [[:B], [:A]], [2, 2],
            copy(group_tags), copy(lig_tags),
            base.n_fit_params_estimate + 2)
        pc = spec_ab.n_fit_params_estimate
        cache = Dict(pc => EnzymeRates.AbstractMechanismSpec[spec_ab, spec_ba])
        EnzymeRates.dedup!(cache)
        @test length(cache[pc]) == 1
    end
end

# ─── expand_mechanisms ─────────────────────────────────────────────────
@testset "expand_mechanisms" begin
    @testset "Returns dict keyed by param count" begin
        # SEED: uni-uni RE-only, 3 singleton kinetic groups → base_pc = 3.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_rxn)
        base_pc = spec.n_fit_params_estimate
        result = EnzymeRates.expand_mechanisms([spec], uni_uni_rxn)
        @test result isa Dict{Int, Vector{AbstractMechanismSpec}}
        @test haskey(result, base_pc + 1)
    end

    @testset "Allosteric expansion included" begin
        # SEED: uni-uni RE-only attached to an oligomeric reaction.
        # expand_mechanisms must include allosteric variants in its output.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)
        result = EnzymeRates.expand_mechanisms([spec], uni_uni_allo)
        allo_count = sum(count(s -> s isa AllostericMechanismSpec, ss)
                         for (_, ss) in result)
        # _expand_to_allosteric on a uni-uni seed with 3 kinetic groups
        # produces n_groups+1=4 allosteric variants (one per group +
        # one for L-only). Other moves do not produce allosteric output
        # from a plain MechanismSpec, so at least 4 exist.
        @test allo_count >= 4
    end

    @testset "No self-expansion to same param count" begin
        # SEED: uni-uni RE-only, base_pc = 3.
        # All keys in the result dict must be strictly greater than base_pc.
        m_seed = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m_seed, uni_uni_rxn)
        base_pc = spec.n_fit_params_estimate
        result = EnzymeRates.expand_mechanisms([spec], uni_uni_rxn)
        for (pc, _) in result
            @test pc > base_pc
        end
    end

    @testset "Allosteric rewrap preserves structure" begin
        # SEED: all-:EqualRT allosteric uni-uni. Passing this to
        # expand_mechanisms must produce AllostericMechanismSpec expansions
        # (RE→SS rewrapped as allosteric, etc.).
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        allo = allosteric_spec_from_mechanism_and_rxn(m_seed, uni_uni_allo)
        result = EnzymeRates.expand_mechanisms([allo], uni_uni_allo)
        allo_results = filter(s -> s isa AllostericMechanismSpec,
                              vcat([ss for (_, ss) in result]...))
        @test !isempty(allo_results)
        # Every rewrapped allosteric result must preserve the input's
        # catalytic_n and base.reaction — base.steps may differ
        # (a base move may have changed them) but the allosteric-side
        # metadata is preserved.
        for r in allo_results
            @test r.catalytic_n == allo.catalytic_n
            @test r.base.reaction === allo.base.reaction
        end
    end

    @testset "Dead-end excludes allosteric regs" begin
        # SEED: allosteric uni-uni with R already added as an allosteric
        # regulator (:OnlyR). expand_mechanisms must never add R as a
        # dead-end inhibitor (R__reg must not appear in any expansion).
        m_seed = @allosteric_mechanism begin
            substrates: S; products: P
            allosteric_regulators: R::OnlyR
            site(:catalytic, 2): begin
                steps: begin
                    E + P ⇌ E_P    :: EqualRT
                    E + S ⇌ E_S    :: EqualRT
                    E_S <--> E_P   :: EqualRT
                end
            end
        end
        with_reg = allosteric_spec_from_mechanism_and_rxn(
            m_seed, uni_uni_allo_reg)
        result = EnzymeRates.expand_mechanisms([with_reg], uni_uni_allo_reg)
        for (_, ss) in result
            for s in ss
                base = s isa AllostericMechanismSpec ? s.base : s
                for step in base.steps
                    for sym in Iterators.flatten(
                            (step.reactants, step.products))
                        @test !contains(string(sym), "R__reg")
                    end
                end
            end
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 6. Integration (enumerate_all)
# ═══════════════════════════════════════════════════════════════════════

# ─── enumerate_all ─────────────────────────────────────────────────────
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
                    @test length(EnzymeRates.fitted_params(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    allo_count += 1
                    allo_count > 3 && continue
                    m = AllostericEnzymeMechanism(spec)
                    @test length(EnzymeRates.fitted_params(m)) <=
                        spec.n_fit_params_estimate
                end
            end
        end
    end

    @testset "Bi-bi full enumeration" begin
        # Sample-based: with full per-group tag enumeration,
        # bi-bi at pc=10 has ~190k specs. Test the upper-bound
        # invariant on a sample of each n_fit_params_estimate bucket.
        results = enumerate_all(bi_bi_rxn; max_params=8)
        @test !isempty(results)
        allo_count = 0
        for (pc, specs) in results
            sample = first(specs, 5)
            for spec in sample
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(EnzymeRates.fitted_params(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    allo_count += 1
                    allo_count > 3 && continue
                    m = AllostericEnzymeMechanism(spec)
                    @test length(EnzymeRates.fitted_params(m)) <=
                        spec.n_fit_params_estimate
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
        results = enumerate_all(rxn; max_params=8)
        has_allo = any(
            any(s isa EnzymeRates.AllostericMechanismSpec for s in specs)
            for (_, specs) in results)
        @test has_allo
        # Sample-based; full enumeration is large.
        allo_count = 0
        for (pc, specs) in results
            sample = first(specs, 5)
            for spec in sample
                if spec isa EnzymeRates.MechanismSpec
                    m = EnzymeMechanism(spec)
                    @test length(EnzymeRates.fitted_params(m)) <= pc
                elseif spec isa EnzymeRates.AllostericMechanismSpec
                    allo_count += 1
                    allo_count > 3 && continue
                    m = AllostericEnzymeMechanism(spec)
                    @test length(EnzymeRates.fitted_params(m)) <=
                        spec.n_fit_params_estimate
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

        # Param-count buckets must form a near-contiguous range.
        # The maximum single-move delta is 4 (SS :NonequalRT split),
        # so consecutive bucket keys can be separated by at most 4.
        pcs = sort(collect(keys(results)))
        @test all(pcs[i+1] - pcs[i] <= 4 for i in 1:length(pcs)-1)
    end
end

# ═══════════════════════════════════════════════════════════════════════
# Testsets covering downstream concerns (canonicalization parameter-naming;
# move-on-allosteric polymorphism). Adjacent to enumeration but tests
# rate-equation-derivation and AllostericEnzymeMechanism integration.
# ═══════════════════════════════════════════════════════════════════════

@testset "Tagged groups exclude T-state params" begin
    specs = EnzymeRates.init_mechanisms(uni_uni_allo)
    spec = first(specs)
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo)

    @testset ":OnlyR binding group: no K_T param" begin
        only_r = first(filter(allo_specs) do r
            any(((g, t),) -> t == :OnlyR &&
                begin
                    is_iso = all(
                        EnzymeRates.step_metabolite(s) === nothing
                        for s in r.base.steps
                        if s.kinetic_group == g)
                    !is_iso
                end,
                r.group_tags)
        end)
        m = AllostericEnzymeMechanism(only_r)
        params = parameters(m)
        @test length(EnzymeRates.fitted_params(m)) ==
            only_r.n_fit_params_estimate
        t_params = filter(
            p -> endswith(string(p), "_T"), params)
        @test isempty(t_params)
    end

    @testset ":OnlyR iso group: no kf_T/kr_T param" begin
        only_r_iso = first(filter(allo_specs) do r
            any(((g, t),) -> t == :OnlyR &&
                all(
                    !s.is_equilibrium &&
                    EnzymeRates.step_metabolite(s) === nothing
                    for s in r.base.steps
                    if s.kinetic_group == g),
                r.group_tags)
        end)
        m = AllostericEnzymeMechanism(only_r_iso)
        params = parameters(m)
        @test length(EnzymeRates.fitted_params(m)) ==
            only_r_iso.n_fit_params_estimate
        t_k_params = filter(
            p -> contains(string(p), "f_T") ||
                 contains(string(p), "r_T"), params)
        @test isempty(t_k_params)
    end

    @testset "t_state_dead with :NonequalRT: K_T in body must be in parameters(Full)" begin
        # K-type allosteric uni-uni: catalytic step is :OnlyR (so
        # `_t_state_dead == true`), but binding steps are :NonequalRT.
        # Bug B.4/A.1: `_all_t_state_names` returns empty when
        # `_t_state_dead == true`, so K1_T/K2_T leak unrenamed into
        # the canonical hash string → cache miss for structurally-
        # equivalent specs whose rep-step indices differ.
        m = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    E_c + S ⇌ E_S    :: NonequalRT
                    E_c + P ⇌ E_P    :: NonequalRT
                    E_S <--> E_P     :: OnlyR
                end
            end
        end
        @test EnzymeRates._t_state_dead(m)
        params_full = parameters(m, Full)
        # K1_T and K2_T are referenced in `den_T` of the body
        # (the binding partition function for :NonequalRT groups
        # is built regardless of `t_state_dead` since `den_T`
        # always appears in the denominator).
        @test :K1_T in params_full
        @test :K2_T in params_full

        # Canonicalizer invariant: every parameter token in the
        # body must be renamed away. After canonicalization, no
        # raw `_T` suffixed names should survive.
        canon, _ = EnzymeRates._canonicalize_rate_eq_with_map(m)
        @test !occursin(r"\bK\d+_T\b", canon)
        @test !occursin(r"\bk\d+[fr]_T\b", canon)
    end
end

@testset "Base-level moves on allosteric specs" begin
    m_uu = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
        end
    end
    spec = mechanism_spec_from_mechanism_and_rxn(m_uu, uni_uni_allo)
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo)
    allo = first(allo_specs)

    @testset "RE→SS on allosteric" begin
        result = EnzymeRates._expand_re_to_ss(allo)
        @test !isempty(result)
        for r in result
            _assert_spec_invariants(r)
            @test r isa EnzymeRates.AllostericMechanismSpec
            @test r.catalytic_n == allo.catalytic_n
            @test r.allosteric_reg_sites ==
                allo.allosteric_reg_sites
            @test r.group_tags == allo.group_tags
            @test r.n_fit_params_estimate > allo.n_fit_params_estimate
        end
    end

end

end # top-level testset
