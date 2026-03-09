# ABOUTME: Test specifications for mechanism enumeration pipeline
# ABOUTME: Defines spec types, reactions, helper functions, and builder functions

using Random

# ── Helper: convert EnzymeMechanism → MechanismSpec ──────────

"""
Convert a compiled EnzymeMechanism back to a MechanismSpec.
Matches mechanism species to enumerated form indices by walking
the reaction graph using the adjacency from enumerate_enzyme_forms.
"""
function mechanism_spec_from_mechanism(
    m::EnzymeMechanism, @nospecialize(rxn::EnzymeReaction);
    n_catalytic_edges::Int=0)
    site_defs, forms = EnzymeRates.enumerate_enzyme_forms(rxn)
    adj = EnzymeRates._build_adjacency(site_defs, forms)

    mech_forms = EnzymeRates.enzyme_forms(m)
    mf_names = Set(n for (n, _) in mech_forms)

    # Free enzyme (first listed, no atoms) → form index 1
    name_to_idx = Dict{Symbol,Int}()
    first_name, first_atoms = mech_forms[1]
    isempty(first_atoms) || error(
        "First enzyme form $first_name has atoms $first_atoms")
    name_to_idx[first_name] = 1

    # Iteratively resolve unknown forms via known ones + reactions
    rxns = EnzymeRates.reactions(m)
    for _ in 1:length(rxns)
        for (lhs, rhs) in rxns
            el = [s for s in lhs if s ∈ mf_names]
            er = [s for s in rhs if s ∈ mf_names]
            length(el) == 1 && length(er) == 1 || continue
            # Determine known/unknown sides
            known, unknown = if haskey(name_to_idx, el[1]) &&
                    !haskey(name_to_idx, er[1])
                el[1], er[1]
            elseif haskey(name_to_idx, er[1]) &&
                    !haskey(name_to_idx, el[1])
                er[1], el[1]
            else
                continue
            end
            kidx = name_to_idx[known]
            met = [s for s in Iterators.flatten((lhs, rhs))
                   if s ∉ mf_names]
            met_sym = isempty(met) ? nothing : met[1]
            for j in 1:length(forms)
                j == kidx && continue
                key = minmax(kidx, j)
                haskey(adj, key) || continue
                adj[key] == met_sym || continue
                name_to_idx[unknown] = j
                break
            end
        end
        length(name_to_idx) == length(mech_forms) && break
    end

    length(name_to_idx) == length(mech_forms) || error(
        "Could not resolve all enzyme forms: " *
        "resolved=$(keys(name_to_idx))")

    eq_steps_tuple = EnzymeRates.equilibrium_steps(m)
    pc = EnzymeRates.param_constraints(m)

    edges = Tuple{Int,Int}[]
    for (lhs, rhs) in rxns
        enz_lhs = [s for s in lhs if s ∈ mf_names]
        enz_rhs = [s for s in rhs if s ∈ mf_names]
        length(enz_lhs) == 1 && length(enz_rhs) == 1 ||
            error("Expected exactly 1 enzyme on each side")
        push!(edges, (name_to_idx[enz_lhs[1]],
                      name_to_idx[enz_rhs[1]]))
    end

    eq_steps = collect(Bool, eq_steps_tuple)
    constraints = [
        (t, c, [(s, sc) for (s, sc) in f])
        for (t, c, f) in pc
    ]

    n_cat = n_catalytic_edges > 0 ?
        n_catalytic_edges : length(edges)
    EnzymeRates.MechanismSpec(
        rxn, edges, n_cat, eq_steps,
        constraints, length(parameters(m)))
end

# ── Test spec types ──────────────────────────────────────────

"""
Tests a single base MechanismSpec through each pipeline stage
independently (not chained). Each stage runs on [base_mechanism]
and the count is compared to the expected value.
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
comparing output count at each stage across all regulator
partitions.
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

# ── Reaction definitions ─────────────────────────────────────
# 8 logical reactions. Regulated reactions have :unknown version
# (for EnumerationTestSpec) and explicit-role versions (for
# StageExpansionTestSpec).

# 1. Uni-Uni, no regulators
const uni_uni = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

# 2. Uni-Uni + 1 regulator
const uni_uni_reg_unknown = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    regulators: I
end
const uni_uni_dead_end_I = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
end
const uni_uni_allosteric_I = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: I
end

# 3a. Uni-Bi, no regulators
const uni_bi = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

# 3b. Uni-Bi + 1 regulator
const uni_bi_reg_unknown = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    regulators: I
end
const uni_bi_dead_end_I = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    dead_end_inhibitors: I
end
const uni_bi_allosteric_I = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: I
end

# 4. Uni-Bi + allosteric regulator (OEM, catalytic_n=2)
const uni_bi_allosteric_I_oem = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: I
end

# 5. Bi-Bi, no regulators
const bi_bi = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

# 6. Bi-Bi Ping-Pong, no regulators
const bi_bi_ping_pong = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

# 7. Bi-Bi Ping-Pong + 1 regulator
const bi_bi_ping_pong_reg_unknown = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    regulators: I
end
const bi_bi_ping_pong_dead_end_I = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    dead_end_inhibitors: I
end
const bi_bi_ping_pong_allosteric_I = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
    allosteric_regulators: I
end

# 8. Bi-Bi + 2 regulators
const bi_bi_two_regs_unknown = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    regulators: I, J
end
const bi_bi_dead_end_I_allosteric_J = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
    allosteric_regulators: J
end

# ── StageExpansionTestSpec builders ──────────────────────────

"""
Build StageExpansionTestSpecs for reactions without regulators.
Uses the first catalytic topology as the base mechanism for each.
"""
function build_no_reg_stage_expansion_specs()
    specs = StageExpansionTestSpec[]

    # --- Uni-Uni (3 forms, 3 edges: 2 RE binding + 1 SS isomerization) ---
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, ES[C], EP[C]
            end
            steps: begin
                [E, S] ⇌ [ES]
                [E, P] ⇌ [EP]
                [ES] <--> [EP]
            end
        end
        base = mechanism_spec_from_mechanism(m, uni_uni)
        push!(specs, StageExpansionTestSpec(;
            name="Uni-Uni (no reg)",
            reaction=uni_uni,
            base_mechanism=base,
            # 2 RE binding edges, each can be toggled RE→SS;
            # at least 1 binding edge must be SS → 2^2 - 1 = 3
            expected_n_ress=3,
            # no regulators → passthrough
            expected_n_general_modifier=1,
            # no regulators → passthrough
            expected_n_essential_activator=1,
            # no regulators → passthrough
            expected_n_dead_end=1,
            # S and P bind different atoms → no equiv groups → passthrough
            expected_n_equivalence=1,
            # single mechanism → no duplicates to remove
            expected_n_dedup=1,
        ))
    end

    # --- Bi-Bi (5 forms, 5 edges: 4 RE binding + 1 SS isomerization) ---
    # First topology: E+B→EB, EB+A→EAB, EAB→EPQ(SS),
    #                 E+P→EP, EP+Q→EPQ
    let
        base = EnzymeRates._catalytic_topologies(bi_bi)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Bi-Bi (no reg)",
            reaction=bi_bi,
            base_mechanism=base,
            # 4 RE binding edges, each can be toggled RE→SS;
            # at least 1 binding edge must be SS → 2^4 - 1 = 15
            expected_n_ress=15,
            # no regulators → passthrough
            expected_n_general_modifier=1,
            # no regulators → passthrough
            expected_n_essential_activator=1,
            # no regulators → passthrough
            expected_n_dead_end=1,
            # each metabolite (A,B,P,Q) binds once → no equiv groups
            expected_n_equivalence=1,
            # single mechanism → no duplicates to remove
            expected_n_dedup=1,
        ))
    end

    # --- Bi-Bi Ping-Pong (5 forms, 5 edges: 4 RE binding + 1 SS) ---
    # First topology: E+B→EB, EB+A→EAB, EAB→EPQ(SS),
    #                 E+Q→EQ, EQ+P→EPQ
    # Transferred atom group X creates intermediate enzyme forms
    let
        base = EnzymeRates._catalytic_topologies(bi_bi_ping_pong)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Bi-Bi Ping-Pong (no reg)",
            reaction=bi_bi_ping_pong,
            base_mechanism=base,
            # 4 RE binding edges, each can be toggled RE→SS;
            # at least 1 binding edge must be SS → 2^4 - 1 = 15
            expected_n_ress=15,
            # no regulators → passthrough
            expected_n_general_modifier=1,
            # no regulators → passthrough
            expected_n_essential_activator=1,
            # no regulators → passthrough
            expected_n_dead_end=1,
            # each metabolite (A,B,P,Q) binds once → no equiv groups
            expected_n_equivalence=1,
            # single mechanism → no duplicates to remove
            expected_n_dedup=1,
        ))
    end

    return specs
end

"""
Build StageExpansionTestSpecs for single-regulator reactions.
Uses the first catalytic topology as the base mechanism for each.
Two variants per reaction: dead-end and allosteric.
"""
function build_single_reg_stage_expansion_specs()
    specs = StageExpansionTestSpec[]

    # --- Uni-Uni + dead-end I ---
    # Base: 3 forms (E, ES, EP), 3 edges (2 RE binding + 1 SS)
    let
        base = EnzymeRates._catalytic_topologies(uni_uni_dead_end_I)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Uni-Uni (dead-end I)",
            reaction=uni_uni_dead_end_I,
            base_mechanism=base,
            dead_end_regs=[:I],
            # 2 RE binding edges toggled; 2^2 - 1 = 3 (at least 1 SS)
            expected_n_ress=3,
            # no allosteric regs → passthrough
            expected_n_general_modifier=1,
            # no allosteric regs → passthrough
            expected_n_essential_activator=1,
            # 3 catalytic forms × I can bind each → 2^3 = 8 subsets
            expected_n_dead_end=8,
            # S and P bind different atoms → no equiv groups
            expected_n_equivalence=1,
            # single mechanism → no duplicates
            expected_n_dedup=1,
            # no allosteric regs → passthrough
            expected_n_allosteric=1,
            expected_n_tr_equiv=1,
            expected_n_oem_dedup=1,
        ))
    end

    # --- Uni-Uni + allosteric I ---
    # Base: 3 forms (E, ES, EP), 3 edges (2 RE binding + 1 SS)
    let
        base = EnzymeRates._catalytic_topologies(uni_uni_allosteric_I)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Uni-Uni (allosteric I)",
            reaction=uni_uni_allosteric_I,
            base_mechanism=base,
            allosteric_regs=[:I],
            # 2 RE binding edges toggled; 2^2 - 1 = 3 (at least 1 SS)
            expected_n_ress=3,
            # 1 allosteric reg → original + 1 general modifier = 2
            expected_n_general_modifier=2,
            # 1 allosteric reg → original + 1 essential activator = 2
            expected_n_essential_activator=2,
            # no dead-end regs → passthrough
            expected_n_dead_end=1,
            # no equiv groups (S, P bind different atoms)
            expected_n_equivalence=1,
            # single mechanism → no duplicates
            expected_n_dedup=1,
            # 1 reg → 1 partition {I} × 2 multiplicities (m=1,2) = 2
            expected_n_allosteric=2,
            # passthrough (no TR equivalence yet)
            expected_n_tr_equiv=2,
            # no duplicates among 2 distinct OEM specs
            expected_n_oem_dedup=2,
        ))
    end

    # --- Uni-Bi + dead-end I ---
    # Base: 4 forms, 4 edges (3 RE binding + 1 SS isomerization)
    let
        base = EnzymeRates._catalytic_topologies(uni_bi_dead_end_I)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Uni-Bi (dead-end I)",
            reaction=uni_bi_dead_end_I,
            base_mechanism=base,
            dead_end_regs=[:I],
            # 3 RE binding edges toggled; 2^3 - 1 = 7 (at least 1 SS)
            expected_n_ress=7,
            # no allosteric regs → passthrough
            expected_n_general_modifier=1,
            # no allosteric regs → passthrough
            expected_n_essential_activator=1,
            # 4 catalytic forms × I can bind each → 2^4 = 16 subsets
            expected_n_dead_end=16,
            # S, P, Q bind different sites → no equiv groups
            expected_n_equivalence=1,
            # single mechanism → no duplicates
            expected_n_dedup=1,
            # no allosteric regs → passthrough
            expected_n_allosteric=1,
            expected_n_tr_equiv=1,
            expected_n_oem_dedup=1,
        ))
    end

    # --- Uni-Bi + allosteric I ---
    # Base: 4 forms, 4 edges (3 RE binding + 1 SS isomerization)
    let
        base = EnzymeRates._catalytic_topologies(uni_bi_allosteric_I)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Uni-Bi (allosteric I)",
            reaction=uni_bi_allosteric_I,
            base_mechanism=base,
            allosteric_regs=[:I],
            # 3 RE binding edges toggled; 2^3 - 1 = 7 (at least 1 SS)
            expected_n_ress=7,
            # 1 allosteric reg → original + 1 general modifier = 2
            expected_n_general_modifier=2,
            # 1 allosteric reg → original + 1 essential activator = 2
            expected_n_essential_activator=2,
            # no dead-end regs → passthrough
            expected_n_dead_end=1,
            # no equiv groups (S, P, Q bind different sites)
            expected_n_equivalence=1,
            # single mechanism → no duplicates
            expected_n_dedup=1,
            # 1 reg → 1 partition {I} × 2 multiplicities (m=1,2) = 2
            expected_n_allosteric=2,
            # passthrough
            expected_n_tr_equiv=2,
            # no duplicates among 2 distinct OEM specs
            expected_n_oem_dedup=2,
        ))
    end

    # --- Bi-Bi Ping-Pong + dead-end I ---
    # Base: 5 forms, 5 edges (4 RE binding + 1 SS isomerization)
    let
        base = EnzymeRates._catalytic_topologies(
            bi_bi_ping_pong_dead_end_I)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Bi-Bi Ping-Pong (dead-end I)",
            reaction=bi_bi_ping_pong_dead_end_I,
            base_mechanism=base,
            dead_end_regs=[:I],
            # 4 RE binding edges toggled; 2^4 - 1 = 15 (at least 1 SS)
            expected_n_ress=15,
            # no allosteric regs → passthrough
            expected_n_general_modifier=1,
            # no allosteric regs → passthrough
            expected_n_essential_activator=1,
            # 5 catalytic forms × I can bind each → 2^5 = 32 subsets
            expected_n_dead_end=32,
            # A, B, P, Q each bind once → no equiv groups
            expected_n_equivalence=1,
            # single mechanism → no duplicates
            expected_n_dedup=1,
            # no allosteric regs → passthrough
            expected_n_allosteric=1,
            expected_n_tr_equiv=1,
            expected_n_oem_dedup=1,
        ))
    end

    # --- Bi-Bi Ping-Pong + allosteric I ---
    # Base: 5 forms, 5 edges (4 RE binding + 1 SS isomerization)
    let
        base = EnzymeRates._catalytic_topologies(
            bi_bi_ping_pong_allosteric_I)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Bi-Bi Ping-Pong (allosteric I)",
            reaction=bi_bi_ping_pong_allosteric_I,
            base_mechanism=base,
            allosteric_regs=[:I],
            # 4 RE binding edges toggled; 2^4 - 1 = 15 (at least 1 SS)
            expected_n_ress=15,
            # 1 allosteric reg → original + 1 general modifier = 2
            expected_n_general_modifier=2,
            # 1 allosteric reg → original + 1 essential activator = 2
            expected_n_essential_activator=2,
            # no dead-end regs → passthrough
            expected_n_dead_end=1,
            # A, B, P, Q each bind once → no equiv groups
            expected_n_equivalence=1,
            # single mechanism → no duplicates
            expected_n_dedup=1,
            # 1 reg → 1 partition {I} × 2 multiplicities (m=1,2) = 2
            expected_n_allosteric=2,
            # passthrough
            expected_n_tr_equiv=2,
            # no duplicates among 2 distinct OEM specs
            expected_n_oem_dedup=2,
        ))
    end

    return specs
end

"""
Build StageExpansionTestSpecs for multi-regulator reactions and
OEM (oligomeric enzyme mechanism) variants.
Uses the first catalytic topology as the base mechanism for each.
"""
function build_multi_reg_stage_expansion_specs()
    specs = StageExpansionTestSpec[]

    # --- Uni-Bi + allosteric I (OEM, catalytic_n=2) ---
    # Base: 4 forms, 4 edges (3 RE binding + 1 SS isomerization)
    # Same catalytic base as Uni-Bi allosteric, but with OEM stages
    let
        base = EnzymeRates._catalytic_topologies(
            uni_bi_allosteric_I_oem)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Uni-Bi (allosteric I, OEM n=2)",
            reaction=uni_bi_allosteric_I_oem,
            base_mechanism=base,
            allosteric_regs=[:I],
            catalytic_n=2,
            # 3 RE binding edges toggled; 2^3 - 1 = 7
            expected_n_ress=7,
            # 1 allosteric reg → original + 1 general modifier = 2
            expected_n_general_modifier=2,
            # 1 allosteric reg → original + 1 essential activator = 2
            expected_n_essential_activator=2,
            # no dead-end regs → passthrough
            expected_n_dead_end=1,
            # S, P, Q bind different sites → no equiv groups
            expected_n_equivalence=1,
            # single mechanism → no duplicates
            expected_n_dedup=1,
            # 1 reg → 1 partition {I} × 2 multiplicities (m=1,2) = 2
            expected_n_allosteric=2,
            # passthrough (no TR equivalence yet)
            expected_n_tr_equiv=2,
            # no duplicates among 2 distinct OEM specs
            expected_n_oem_dedup=2,
        ))
    end

    # --- Bi-Bi + dead-end I + allosteric J ---
    # Base: 5 forms, 5 edges (4 RE binding + 1 SS isomerization)
    # catalytic_n=0 (not OEM), so OEM stages are skipped
    let
        base = EnzymeRates._catalytic_topologies(
            bi_bi_dead_end_I_allosteric_J)[1]
        push!(specs, StageExpansionTestSpec(;
            name="Bi-Bi (dead-end I, allosteric J)",
            reaction=bi_bi_dead_end_I_allosteric_J,
            base_mechanism=base,
            dead_end_regs=[:I],
            allosteric_regs=[:J],
            # 4 RE binding edges toggled; 2^4 - 1 = 15
            expected_n_ress=15,
            # 1 allosteric reg J → original + 1 general modifier = 2
            expected_n_general_modifier=2,
            # 1 allosteric reg J → original + 1 essential activator = 2
            expected_n_essential_activator=2,
            # 5 catalytic forms × I can bind each → 2^5 = 32 subsets
            expected_n_dead_end=32,
            # A, B, P, Q each bind once → no equiv groups
            expected_n_equivalence=1,
            # single mechanism → no duplicates
            expected_n_dedup=1,
        ))
    end

    return specs
end

function build_stage_expansion_specs()
    return vcat(
        build_no_reg_stage_expansion_specs(),
        build_single_reg_stage_expansion_specs(),
        build_multi_reg_stage_expansion_specs(),
    )
end

# ── Helper: run pipeline stage by stage (all partitions) ─────

"""
Run the full enumeration pipeline stage by stage, summing
counts across all regulator partitions. Returns a NamedTuple
with per-stage totals.
"""
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

# ── EnumerationTestSpec builder ──────────────────────────────

"""
Build EnumerationTestSpecs for end-to-end pipeline testing.
Each spec tests a reaction through the full enumeration pipeline
including regulator partitioning for :unknown regulators.
"""
function build_enumeration_specs()
    specs = EnumerationTestSpec[]

    push!(specs, EnumerationTestSpec(;
        name="Uni-Uni, no regs",
        reaction=uni_uni,
        expected_n_forms=3,
        expected_n_catalytic=1,
        expected_n_ress=3,
        expected_n_general_modifier=3,
        expected_n_essential_activator=3,
        expected_n_dead_end=3,
        expected_n_equivalence=3,
        expected_n_dedup=1,
        expected_n_total=1,
    ))

    push!(specs, EnumerationTestSpec(;
        name="Uni-Uni + 1 unknown reg",
        reaction=uni_uni_reg_unknown,
        expected_n_forms=6,
        expected_n_catalytic=1,
        # 2 partitions (de/allo), each 3 RE/SS = 6
        expected_n_ress=6,
        expected_n_general_modifier=9,
        expected_n_essential_activator=12,
        expected_n_dead_end=33,
        expected_n_equivalence=48,
        expected_n_dedup=28,
        expected_n_total=28,
    ))

    push!(specs, EnumerationTestSpec(;
        name="Uni-Bi, no regs",
        reaction=uni_bi,
        expected_n_forms=11,
        expected_n_catalytic=3,
        expected_n_ress=41,
        expected_n_general_modifier=41,
        expected_n_essential_activator=41,
        expected_n_dead_end=41,
        expected_n_equivalence=41,
        expected_n_dedup=9,
        expected_n_total=9,
    ))

    push!(specs, EnumerationTestSpec(;
        name="Uni-Bi + 1 unknown reg",
        reaction=uni_bi_reg_unknown,
        expected_n_forms=22,
        expected_n_catalytic=3,
        expected_n_ress=82,
        expected_n_general_modifier=123,
        expected_n_essential_activator=164,
        expected_n_dead_end=1211,
        expected_n_equivalence=1606,
        expected_n_dedup=697,
        expected_n_total=697,
    ))

    push!(specs, EnumerationTestSpec(;
        name="Uni-Bi + allosteric, OEM n=2",
        reaction=uni_bi_allosteric_I_oem,
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
    ))

    push!(specs, EnumerationTestSpec(;
        name="Bi-Bi, no regs",
        reaction=bi_bi,
        expected_n_forms=11,
        expected_n_catalytic=9,
        expected_n_ress=527,
        expected_n_general_modifier=527,
        expected_n_essential_activator=527,
        expected_n_dead_end=527,
        expected_n_equivalence=958,
        expected_n_dedup=207,
        expected_n_total=207,
    ))

    push!(specs, EnumerationTestSpec(;
        name="Bi-Bi Ping-Pong, no regs",
        reaction=bi_bi_ping_pong,
        expected_n_forms=17,
        expected_n_catalytic=10,
        expected_n_ress=558,
        expected_n_general_modifier=558,
        expected_n_essential_activator=558,
        expected_n_dead_end=558,
        expected_n_equivalence=989,
        expected_n_dedup=210,
        expected_n_total=210,
    ))

    push!(specs, EnumerationTestSpec(;
        name="Bi-Bi Ping-Pong + 1 unknown reg",
        reaction=bi_bi_ping_pong_reg_unknown,
        expected_n_forms=34,
        expected_n_catalytic=10,
        expected_n_ress=1116,
        expected_n_general_modifier=1674,
        expected_n_essential_activator=2232,
        expected_n_dead_end=50250,
        expected_n_equivalence=106646,
        expected_n_dedup=35113,
        expected_n_total=35113,
    ))

    # Bi-Bi + 2 unknown regs: too large to enumerate in
    # reasonable time/memory (4 partitions × 9 topologies ×
    # combinatorial dead-end expansion on 11+ forms).
    # Skipped from end-to-end testing.

    specs
end
