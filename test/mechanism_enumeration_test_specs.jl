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

    # Free enzyme (no atoms) → form index 1
    name_to_idx = Dict{Symbol,Int}()
    for (mf_name, mf_atoms) in mech_forms
        if isempty(mf_atoms)
            name_to_idx[mf_name] = 1
        end
    end

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

# 3. Uni-Bi + 1 regulator
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
