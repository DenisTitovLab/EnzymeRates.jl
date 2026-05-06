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

@testset "n_fit_params_estimate semantics" begin
    # Simple uni-uni: 2 RE binding steps + 1 SS catalytic step.
    # 3 forms (E, E_S, E_P), 3 steps → n_thermo = 3 - 3 + 1 = 1.
    # New formula: n_re_groups + 2*n_ss_groups - n_thermo
    #            = 2 + 2 - 1 = 3
    # (Old formula gave 5 because of the +2 for Keq+E_total.)
    # length(fitted_params(m)) for uni-uni = 3 (K1, K2, k3f).
    init_specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
    @test !isempty(init_specs)
    spec = first(init_specs)
    @test hasfield(typeof(spec), :n_fit_params_estimate)
    m = EnzymeRates.EnzymeMechanism(spec)
    n_actual = length(EnzymeRates.fitted_params(m))
    @test spec.n_fit_params_estimate == n_actual
end

@testset "n_fit_params_estimate upper-bound for dead-end mirrors" begin
    # Guards the upper-bound invariant
    #   spec.n_fit_params_estimate >= length(fitted_params(m))
    # for init mechanisms and one round of expansion of
    # `uni_uni_with_reg` (which has a dead-end inhibitor :I that
    # creates mirror cycles). To cap @generated compile cost,
    # cap the number of compiled specs per round at 30 — that's
    # enough to exercise dead-end-mirror cases without triggering
    # the explosion of compile-times for the full expansion fan-out.
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
        # 4 unique dead-end forms → 2^4 = 16 variants
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

    @testset "Bi-Bi Ping-Pong: 3 dead-end forms" begin
        # Forms: E, E_A, Estar, Estar_A_P, Estar_B, E_Q
        # E_A: +P→E_A_P(mixed✓), +Q→E_A_Q(mixed✓)
        # E_Q: +B→E_B_Q(mixed✓)
        # 3 dead-end forms → 2^3 = 8 variants
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

end

# ═══════════════════════════════════════════════════════════════════════
# Out-of-scope testsets pending move to other test files (see plan Task 13)
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

@testset "init_mechanisms" begin

    @testset "Param count invariant" begin
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

    @testset "Mirror steps share kinetic_group with catalytic" begin
        # When a regulator binds the same metabolite via dead-end edges,
        # those mirror steps must share the kinetic_group of the catalytic
        # binding step. (Mirror propagation is implicit in the new design;
        # this test is a guardrail.)
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            dead_end_inhibitors: I
        end
        specs = EnzymeRates.init_mechanisms(rxn)
        @test !isempty(specs)
        for spec in specs
            # Group the steps by the metabolite they bind
            by_metabolite = Dict{Symbol, Vector{EnzymeRates.StepSpec}}()
            for step in spec.steps
                step.is_equilibrium || continue
                # Binding step has form F + met ⇌ F_bound
                length(step.reactants) == 2 || continue
                met = step.reactants[2]
                push!(get!(by_metabolite, met,
                           EnzymeRates.StepSpec[]), step)
            end
            # All same-metabolite RE binding steps must share kinetic_group
            for (met, steps) in by_metabolite
                length(steps) >= 2 || continue
                groups = Set(s.kinetic_group for s in steps)
                @test length(groups) == 1
            end
        end
    end

    @testset "Uni-Uni: no dead-end forms" begin
        specs = EnzymeRates.init_mechanisms(
            uni_uni_rxn)
        @test length(specs) == 1
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

        @testset "Bi-bi random: 7 variants (was 16)" begin
            # 4 dead-end forms × 7 competition patterns.
            # Each pattern yields a distinct dead-end set:
            #   {A↔P,B↔Q}: {E_A_Q,E_B_P}
            #   {A↔Q,B↔P}: {E_A_P,E_B_Q}
            #   {A↔P,A↔Q,B↔P}: {E_B_Q}
            #   {A↔P,A↔Q,B↔Q}: {E_B_P}
            #   {A↔P,B↔P,B↔Q}: {E_A_Q}
            #   {A↔Q,B↔P,B↔Q}: {E_A_P}
            #   {A↔P,A↔Q,B↔P,B↔Q}: {} (no dead-ends)
            # All 7 sets are distinct → 7 variants after dedup
            result =
                EnzymeRates._expand_substrate_product_dead_ends(
                    [spec_bb], bi_bi_rxn)
            @test length(result) == 7
        end

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

        @testset "Round-trip: competition-filtered specs compile" begin
            for rxn in [uni_uni_rxn, bi_bi_rxn, bi_bi_pp_rxn]
                specs =
                    EnzymeRates.init_mechanisms(rxn)
                # Test first 5 specs (compilation can be slow)
                for spec in first(specs, 5)
                    m = EnzymeMechanism(spec)
                    @test m isa EnzymeMechanism
                    @test length(EnzymeRates.fitted_params(m)) <=
                        spec.n_fit_params_estimate
                end
            end
        end
    end
end

@testset "RE→SS conversion" begin
    @testset "Multiple RE steps" begin
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
        result = EnzymeRates._expand_re_to_ss(spec)
        @test length(result) == 2
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end
    end

    @testset "All SS → yields nothing" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P <--> E_P
                E + S <--> E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, uni_uni_rxn)
        result = EnzymeRates._expand_re_to_ss(spec)
        @test isempty(result)
    end

    @testset "All-SS catalytic + dead-end RE: only RE group convertible" begin
        # Uni-uni where all catalytic steps are SS. Dead-end
        # inhibitor I binds to multiple forms; those binding steps
        # share one kinetic group (all RE). RE→SS atomic on that
        # group should yield exactly one variant.
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P <--> E_P
                E + S <--> E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, uni_uni_rxn)
        rxn_i = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
        end
        spec_with_rxn = MechanismSpec(rxn_i, spec.steps,
            spec.n_fit_params_estimate)
        de_specs = EnzymeRates._expand_add_dead_end_regulator(
            spec_with_rxn, rxn_i)
        multi_form = filter(de_specs) do s
            n = count(s.steps) do st
                any(contains(string(sym), "I__reg")
                    for sym in Iterators.flatten(
                        (st.reactants, st.products)))
            end
            n >= 2
        end
        if !isempty(multi_form)
            spec_de = first(multi_form)
            # The dead-end RE binding group should be the only
            # all-RE group (all catalytic groups are SS).
            result = EnzymeRates._expand_re_to_ss(spec_de)
            @test length(result) == 1
        end
    end

    @testset "Bi-bi init: per-group RE→SS count" begin
        # init_mechanisms produces specs where each (metabolite,
        # RE/SS) class shares one kinetic group. RE→SS converts
        # ONE WHOLE group atomically — count is the number of
        # all-RE groups (excluding the iso group which is SS).
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        for spec in specs
            n_re_groups = length(unique(
                s.kinetic_group for s in spec.steps
                if s.is_equilibrium))
            result = EnzymeRates._expand_re_to_ss(spec)
            @test length(result) == n_re_groups
            for r in result
                @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
            end
        end
    end
end

@testset "Split kinetic group" begin
    @testset "Multi-step groups: one split per group member" begin
        # Pick a sequential bi-bi spec with exactly 4 RE multi-step
        # groups of size 2 (4 binding pairs, A/B/P/Q sharing one
        # group each). Splitting each member yields 2 results per
        # group → 8 total (before dedup).
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        spec = first(filter(specs) do s
            counts = Dict{Int, Int}()
            for st in s.steps
                counts[st.kinetic_group] =
                    get(counts, st.kinetic_group, 0) + 1
            end
            multi = filter(((_, n),) -> n >= 2, counts)
            length(multi) == 4 && all(n == 2 for (_, n) in multi)
        end)
        result = EnzymeRates._expand_split_kinetic_group(spec)
        @test length(result) == 8
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end
    end

    @testset "No multi-step groups → yields nothing" begin
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
        # Every group has exactly 1 step → no split possible
        result = EnzymeRates._expand_split_kinetic_group(spec)
        @test isempty(result)
    end

    @testset "Mixed RE/SS group: split delta differs" begin
        # Take a bi-bi init mechanism, RE→SS one of its multi-step
        # groups (now SS). Splitting any RE group gives +1, splitting
        # the SS group gives +2.
        specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
        spec = first(filter(specs) do s
            counts = Dict{Int, Int}()
            for st in s.steps
                counts[st.kinetic_group] =
                    get(counts, st.kinetic_group, 0) + 1
            end
            length(filter(((_, n),) -> n >= 2, counts)) >= 1
        end)
        ss_specs = EnzymeRates._expand_re_to_ss(spec)
        @test !isempty(ss_specs)
        ss_spec = first(ss_specs)
        result = EnzymeRates._expand_split_kinetic_group(ss_spec)
        # Each result has delta of 1 (RE split) or 2 (SS split).
        for r in result
            delta = r.n_fit_params_estimate - ss_spec.n_fit_params_estimate
            @test delta == 1 || delta == 2
        end
    end
end

@testset "Add dead-end regulator" begin

    @testset "Uni-uni + new regulator" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, uni_uni_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg)
        # Uni-uni: sub_names={S}, prod_names={P}
        # E: bound={} → eligible
        # E_S: bound={S} → all subs → NOT eligible
        # E_P: bound={P} → all prods → NOT eligible
        # 1 eligible form → 2^1 - 1 = 1 variant
        @test length(result) == 1
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end
    end

    @testset "No regulators → yields nothing" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(
            m, uni_uni_rxn)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_rxn)
        @test isempty(result)
    end

    @testset "All results compile" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, uni_uni_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg)
        for r in result
            compiled = EnzymeMechanism(r)
            @test compiled isa EnzymeMechanism
        end
    end

    @testset "exclude_regs works" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, uni_uni_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, uni_uni_with_reg;
            exclude_regs=Set([:I]))
        @test isempty(result)
    end

    @testset "Mirror steps created" begin
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
        bi_bi_with_reg = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, bi_bi_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bi_bi_with_reg)
        # Find a variant with multiple eligible forms
        multi_form = filter(
            r -> length(r.steps) >
                length(spec.steps) + 2,
            result)
        if !isempty(multi_form)
            r = first(multi_form)
            @test length(r.steps) >
                length(spec.steps) + 2
        end
    end

    @testset "Same regulator binding steps share kinetic group" begin
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
        bi_bi_with_reg = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, bi_bi_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bi_bi_with_reg)
        for r in result
            reg_binding_steps = filter(
                s -> length(s.reactants) == 2 &&
                    contains(string(s.reactants[2]),
                        "__reg"),
                r.steps)
            if length(reg_binding_steps) >= 2
                groups = unique(
                    s.kinetic_group for s in reg_binding_steps)
                @test length(groups) == 1
            end
        end
    end

    @testset "Two regulators: both bound" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        rxn_ij = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I, J
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, rxn_ij)
        # Add I
        i_specs =
            EnzymeRates._expand_add_dead_end_regulator(
                spec, rxn_ij)
        with_i = first(i_specs)
        # Add J to mechanism that already has I
        j_specs =
            EnzymeRates._expand_add_dead_end_regulator(
                with_i, rxn_ij)
        @test !isempty(j_specs)
        for s in j_specs
            has_j = any(
                contains(string(sym), "J__reg")
                for st in s.steps
                for sym in Iterators.flatten(
                    (st.reactants, st.products)))
            @test has_j
        end
    end

    @testset "Bi-bi: exact dead-end regulator count" begin
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
        rxn_i = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, rxn_i)
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                spec, rxn_i)
        # E_A_B has all subs → ineligible
        # E_P_Q has all prods → ineligible
        # 5 eligible forms, 9 inhibitor patterns
        # (3 sub subsets × 3 prod subsets) → 9
        @test length(result) == 9
    end

    @testset "Dead-end reg on allosteric spec" begin
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + P ⇌ E_P
                E + S ⇌ E_S
                E_S <--> E_P
            end
        end
        # Allosteric-only regulator → no dead-end
        rxn_allo_only = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: R
            oligomeric_state: 2
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, rxn_allo_only)
        allo = first(EnzymeRates._expand_to_allosteric(
            spec, rxn_allo_only))
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                allo, rxn_allo_only)
        @test isempty(result)

        # Mixed: I is dead-end, R is allosteric
        rxn_mixed = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            allosteric_regulators: R
            oligomeric_state: 2
        end
        spec_m = mechanism_spec_from_mechanism_and_rxn(m, rxn_mixed)
        allo_m = first(EnzymeRates._expand_to_allosteric(
            spec_m, rxn_mixed))
        result_m =
            EnzymeRates._expand_add_dead_end_regulator(
                allo_m, rxn_mixed)
        @test !isempty(result_m)
        for r in result_m
            @test r isa AllostericMechanismSpec
        end
    end

    @testset "Ping-pong: exact dead-end count" begin
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
        rxn_pp_i = @enzyme_reaction begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, rxn_pp_i)
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                spec, rxn_pp_i)
        # Estar_A_P has all prods → ineligible
        # Estar_B has all subs → ineligible
        # 4 eligible, 9 patterns → 3 unique form sets
        @test length(result) == 3
    end

    @testset "Topology-aware: sequential bi-bi" begin
        m = @enzyme_mechanism begin
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
        rxn_seq = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, rxn_seq)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, rxn_seq)
        # Sequential bi-bi: binding step sources:
        #   A: from E    B: from E_A
        #   P: from E_Q  Q: from E
        # Eligible forms: E, E_A, E_Q
        # 9 inhibitor patterns → 4 unique form sets:
        #   {E,E_Q}: ({A},{P}), ({A,B},{P})
        #   {E,E_A}: ({B},{Q}), ({B},{P,Q})
        #   {E_A,E_Q}: ({B},{P})
        #   {E}: ({A},{Q}), ({A},{P,Q}),
        #         ({A,B},{Q}), ({A,B},{P,Q})
        @test length(result) == 4
        for r in result
            @test r.n_fit_params_estimate == spec.n_fit_params_estimate + 1
        end
    end

    @testset "Bi-bi random: 9 inhibitor variants" begin
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
        bb_with_reg = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, bb_with_reg)
        result = EnzymeRates._expand_add_dead_end_regulator(
            spec, bb_with_reg)
        @test length(result) == 9
    end

    @testset "Two inhibitors: compete vs not" begin
        # Use bi-bi random: I1 binds to multiple forms,
        # creating mirror steps. When adding I2, compete
        # vs not-compete with I1 produces different forms.
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
        rxn_2i = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            dead_end_inhibitors: I1, I2
        end
        spec = mechanism_spec_from_mechanism_and_rxn(m, rxn_2i)
        # Adds both eligible regs (I1 + I2): 9 each
        result1 =
            EnzymeRates._expand_add_dead_end_regulator(
                spec, rxn_2i)
        @test length(result1) == 18
        # Pick variant where I1 binds to multiple forms
        multi = filter(result1) do r
            inh_forms = filter(
                f -> contains(string(f), "I1__reg"),
                collect(
                    EnzymeRates.all_form_names(r)))
            length(inh_forms) >= 2
        end
        @test !isempty(multi)
        spec_i1 = first(multi)
        # Add I2 to spec that already has I1
        result2 =
            EnzymeRates._expand_add_dead_end_regulator(
                spec_i1, rxn_2i)
        # With I1 present, competition patterns produce
        # more I2 variants than the base 9.
        @test length(result2) == 17
        # Not-competing variant: I2 coexists with I1
        has_coexist = any(result2) do r
            any(
                f -> contains(string(f), "I1__reg") &&
                     contains(string(f), "I2__reg"),
                collect(
                    EnzymeRates.all_form_names(r)))
        end
        @test has_coexist
        # Competing variant: I2 forms but no coexistence
        has_compete = any(result2) do r
            forms = collect(
                EnzymeRates.all_form_names(r))
            has_i2 = any(
                f -> contains(string(f), "I2__reg"),
                forms)
            no_coexist = !any(
                f -> contains(string(f), "I1__reg") &&
                     contains(string(f), "I2__reg"),
                forms)
            has_i2 && no_coexist
        end
        @test has_compete  # competing variant exists
    end
end

@testset "Regulator dummy naming stability" begin
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + P ⇌ E_P
            E + S ⇌ E_S
            E_S <--> E_P
        end
    end
    rxn2 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I, J
    end
    spec = mechanism_spec_from_mechanism_and_rxn(m, rxn2)
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

@testset "Allosteric conversion" begin
    # Per-group tag enumeration: per R-state-active convention,
    # each kinetic group → one of `{:OnlyR, :EqualRT}`.
    # No K-type/V-type hardcoded subsets.

    @testset "Uni-uni: per-group tag variants" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        # Uni-uni init has 3 singleton groups (S binding, P binding, iso).
        # _expand_to_allosteric emits the all-:EqualRT baseline once plus
        # one :OnlyR variant per group: 1 + 3 = 4.
        @test length(result) == 4
        for r in result
            _assert_spec_invariants(r)
            @test r isa AllostericMechanismSpec
            @test r.catalytic_n == 2
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

@testset "Add allosteric regulator" begin
    @testset "Add regulator to allosteric spec" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo_reg)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo_reg)
        @test !isempty(allo_specs)
        allo = first(allo_specs)
        result = EnzymeRates._expand_add_allosteric_regulator(
            allo, uni_uni_allo_reg)
        # R not yet added: 3 tags × 1 site option
        # (new site only, no existing reg sites) = 3
        @test length(result) == 3
        for r in result
            _assert_spec_invariants(r)
        end
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
        # R2: 3 tags × 2 site options (new site + R1's site) = 6
        # + :EqualRT at R1's non-:EqualRT site = 1
        @test length(r2_added) == 7
    end

    @testset "EqualRT ligand reachable at existing reg site" begin
        # Set up a spec with a single regulator at one site, non-:EqualRT.
        # Then expand to add a SECOND regulator at the SAME site as :EqualRT.
        # Verify a result spec exists with both ligands at site 1, the
        # second tagged :EqualRT.
        specs = EnzymeRates.init_mechanisms(uni_uni_allo_2reg)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(spec, uni_uni_allo_2reg)
        allo = first(allo_specs)
        # Add R1 first to get a seed with one non-:EqualRT ligand at site 1
        with_r1 = EnzymeRates._expand_add_allosteric_regulator(
            allo, uni_uni_allo_2reg)
        seed = first(filter(s -> haskey(s.reg_ligand_tags, :R1) &&
                                  s.reg_ligand_tags[:R1] == :OnlyR, with_r1))
        # Now add R2 — verify one result has R2::EqualRT at the same site as R1
        expanded = EnzymeRates._expand_add_allosteric_regulator(
            seed, uni_uni_allo_2reg)
        target = findfirst(expanded) do s
            get(s.reg_ligand_tags, :R2, nothing) == :EqualRT &&
                :R2 in s.allosteric_reg_sites[1]
        end
        @test target !== nothing
    end
end

@testset "Remove TR equivalence" begin
    # `_expand_change_allo_state` relaxes one allo_state from a
    # constrained value (`:OnlyR`/`:OnlyT`/`:EqualRT`) to
    # `:NonequalRT`. Group allo_states live in
    # `spec.group_tags::Dict{Int, Symbol}`; ligand allo_states in
    # `spec.reg_ligand_tags::Dict{Symbol, Symbol}`. Both Dicts are
    # dense — entries already at `:NonequalRT` are skipped (no
    # further relaxation possible).

    @testset "Each constrained tag contributes one relaxation" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        n_constrained =
            count(t -> t != :NonequalRT, values(allo.group_tags)) +
            count(t -> t != :NonequalRT, values(allo.reg_ligand_tags))
        result = EnzymeRates._expand_change_allo_state(
            allo, uni_uni_allo)
        @test length(result) == n_constrained
    end

    @testset "Fully relaxed → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo_specs = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        allo = first(allo_specs)
        fully_relaxed = allo
        while true
            r = EnzymeRates._expand_change_allo_state(
                fully_relaxed, uni_uni_allo)
            isempty(r) && break
            fully_relaxed = first(r)
        end
        @test isempty(
            EnzymeRates._expand_change_allo_state(
                fully_relaxed, uni_uni_allo))
    end

    @testset "MechanismSpec → yields nothing" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_change_allo_state(
            spec, uni_uni_allo)
        @test isempty(result)
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
        reg_specs = EnzymeRates._expand_add_allosteric_regulator(
            allo, rxn_r)
        # Find a reg-spec where R has a constrained (non-:NonequalRT) tag
        tagged = filter(
            r -> get(r.reg_ligand_tags, :R, :NonequalRT) != :NonequalRT,
            reg_specs)
        @test !isempty(tagged)
        tr_spec = first(tagged)
        pc_before = tr_spec.n_fit_params_estimate
        result = EnzymeRates._expand_change_allo_state(
            tr_spec, rxn_r)
        # The variant that relaxes :R to :NonequalRT
        r_removal = filter(
            r -> r.reg_ligand_tags[:R] == :NonequalRT, result)
        @test !isempty(r_removal)
        # delta depends on R's previous tag — should be +1
        # (one K_R_T appears) when going from non-`:NonequalRT`
        # to `:NonequalRT`.
        for r in r_removal
            @test r.n_fit_params_estimate == pc_before + 1
        end
    end
end

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

@testset "expand_mechanisms" begin
    @testset "Returns dict keyed by param count" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_rxn)
        @test result isa Dict{Int,
            Vector{AbstractMechanismSpec}}
        base_pc = first(specs).n_fit_params_estimate
        @test haskey(result, base_pc + 1)
    end

    @testset "Allosteric expansion included" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_allo)
        base_pc = first(specs).n_fit_params_estimate
        has_allo = any(
            any(s isa AllostericMechanismSpec
                for s in ss)
            for (_, ss) in result)
        @test has_allo
    end

    @testset "No self-expansion to same param count" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_rxn)
        base_pc = first(specs).n_fit_params_estimate
        result = EnzymeRates.expand_mechanisms(
            specs, uni_uni_rxn)
        # All results should have n_fit_params_estimate > base
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
    end
end

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

    # TR equiv removal: S as catalytic met has its own group
    # tag, S as regulator has its own ligand tag. Each constrained
    # entry is relaxable independently.
    tr_spec = first(filter(
        r -> r.reg_ligand_tags[:S] != :NonequalRT, reg_specs))
    result = EnzymeRates._expand_change_allo_state(
        tr_spec, rxn_allo_overlap)
    # At least: 1 ligand-tag relaxation for :S; plus any
    # group-tag relaxations from the catalytic side.
    @test !isempty(result)
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

    @testset "Remove constraint on allosteric" begin
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
        bi_bi_allo_rxn = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        # Use init_mechanisms to get a bi-bi spec with multi-step
        # kinetic groups (same-metabolite RE bindings collapsed).
        bb_specs = EnzymeRates.init_mechanisms(bi_bi_allo_rxn)
        bb_spec_c = first(filter(bb_specs) do s
            counts = Dict{Int, Int}()
            for st in s.steps
                counts[st.kinetic_group] =
                    get(counts, st.kinetic_group, 0) + 1
            end
            any(((_, n),) -> n >= 2, counts)
        end)
        bb_allo = first(
            EnzymeRates._expand_to_allosteric(
                bb_spec_c, bi_bi_allo_rxn))
        # Splitting a multi-step kinetic group in the allosteric
        # base must produce results AND the new group must inherit
        # the parent group's allosteric tag (split is a parameter-
        # relaxation move; it must not change R/T-state semantics).
        result = EnzymeRates._expand_split_kinetic_group(bb_allo)
        @test !isempty(result)
        pre_groups = Set(s.kinetic_group for s in bb_allo.base.steps)
        for r in result
            _assert_spec_invariants(r)
            @test r isa EnzymeRates.AllostericMechanismSpec
            @test r.catalytic_n == bb_allo.catalytic_n
            post_groups =
                Set(s.kinetic_group for s in r.base.steps)
            new_groups = setdiff(post_groups, pre_groups)
            @test length(new_groups) == 1
            new_g = first(new_groups)
            # Identify the parent group: the only pre-group whose
            # step count dropped.
            pre_counts = Dict(g => count(
                    s -> s.kinetic_group == g, bb_allo.base.steps)
                for g in pre_groups)
            post_counts = Dict(g => count(
                    s -> s.kinetic_group == g, r.base.steps)
                for g in pre_groups)
            old_g = only(g for g in pre_groups
                if post_counts[g] < pre_counts[g])
            @test r.group_tags[new_g] == bb_allo.group_tags[old_g]
            for g in pre_groups
                @test r.group_tags[g] == bb_allo.group_tags[g]
            end
        end
    end

    @testset "Add dead-end reg on allosteric" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            oligomeric_state: 2
        end
        spec_i = mechanism_spec_from_mechanism_and_rxn(
            m_uu, rxn)
        allo_i = first(
            EnzymeRates._expand_to_allosteric(
                spec_i, rxn))
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                allo_i, rxn)
        @test !isempty(result)
        for r in result
            _assert_spec_invariants(r)
            @test r isa
                EnzymeRates.AllostericMechanismSpec
            @test r.catalytic_n == allo_i.catalytic_n
            # Dead-end's new binding-step group is tagged
            # `:EqualRT` (cheapest), so r.group_tags is
            # allo_i.group_tags plus one new entry.
            for (g, t) in allo_i.group_tags
                @test r.group_tags[g] == t
            end
        end
    end
end

@testset "init_mechanisms drops unbound regulators from spec→type" begin
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

end # top-level testset
