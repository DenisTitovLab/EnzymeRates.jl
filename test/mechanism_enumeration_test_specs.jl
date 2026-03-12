# ABOUTME: Reaction definitions and helpers for mechanism enumeration tests.
# ABOUTME: Defines EnzymeReaction constants and pipeline stage runner.

using Random

# ── Helper: convert EnzymeMechanism → MechanismSpec ──────────

"""
Convert a compiled EnzymeMechanism back to a MechanismSpec.
Reconstructs StepSpec entries from the mechanism's reactions.
"""
function mechanism_spec_from_mechanism(
    m::EnzymeMechanism,
    @nospecialize(rxn::EnzymeReaction))
    rxns = EnzymeRates.reactions(m)
    eq_steps_tuple = EnzymeRates.equilibrium_steps(m)
    pc = EnzymeRates.param_constraints(m)

    steps = EnzymeRates.StepSpec[]
    for (i, (lhs, rhs)) in enumerate(rxns)
        push!(steps, EnzymeRates.StepSpec(
            collect(lhs), collect(rhs),
            eq_steps_tuple[i]))
    end

    constraints = [
        (t, c, [(s, sc) for (s, sc) in f])
        for (t, c, f) in pc
    ]

    EnzymeRates.MechanismSpec(
        rxn, steps, constraints,
        length(parameters(m)))
end

# ── Reaction definitions ─────────────────────────────────────

# --- No regulators ---

const uni_uni = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

const uni_bi = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
end

const bi_bi = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
end

const bi_bi_ping_pong = @enzyme_reaction begin
    substrates: A[CX], B[N]
    products: P[C], Q[NX]
end

const ter_ter = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end

# --- Dead-end inhibitor variants ---

const uni_uni_dead_end_I = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
end

const uni_bi_dead_end_I = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    dead_end_inhibitors: I
end

const bi_bi_dead_end_I = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
end

const uni_uni_dead_end_I_J = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I, J
end

# --- Allosteric regulator variants ---

const uni_uni_allosteric_R = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: R
end

const uni_bi_allosteric_R = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    allosteric_regulators: R
end

# --- Mixed ---

const bi_bi_dead_end_I_allosteric_R = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
    allosteric_regulators: R
end

const bi_bi_allosteric_R1_R2 = @enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    allosteric_regulators: R1, R2
end

# --- Unknown role (for end-to-end partitioning) ---

const uni_uni_reg_unknown = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    regulators: I
end

const uni_bi_reg_unknown = @enzyme_reaction begin
    substrates: S[AB]
    products: P[A], Q[B]
    regulators: I
end

# ── Helper: run pipeline stage by stage (all partitions) ─────

function _run_full_pipeline_stages(rxn; catalytic_n::Int=0,
                                   max_re_groups::Int=7)
    catalytic = EnzymeRates._catalytic_topologies(rxn)

    roles = EnzymeRates.regulator_roles(rxn)
    fixed_de = Symbol[r[1] for r in roles
                       if r[2] == :dead_end]
    fixed_allo = Symbol[r[1] for r in roles
                         if r[2] == :allosteric]
    unknown = Symbol[r[1] for r in roles
                      if r[2] == :unknown]
    n_unknown = length(unknown)

    n_ress = 0; n_de = 0; n_eq = 0; n_dd = 0
    n_allo = 0; n_tr = 0; n_allo_dd = 0

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

        with_de = EnzymeRates._expand_dead_end_inhibitors(
            ress, rxn; dead_end_regs=de)
        n_de += length(with_de)

        with_eq =
            EnzymeRates._expand_equivalence_constraints(
                with_de, rxn)
        n_eq += length(with_eq)

        dd = EnzymeRates._deduplicate(with_eq, rxn)
        n_dd += length(dd)

        if !isempty(al)
            cn = catalytic_n > 0 ? catalytic_n : 1
            allo = EnzymeRates._expand_allosteric(
                dd, rxn; catalytic_n=cn,
                allosteric_regs=al)
            n_allo += length(allo)
            allo = EnzymeRates._expand_tr_equivalence(
                allo, rxn)
            n_tr += length(allo)
            allo = EnzymeRates._deduplicate_allosteric(
                allo, rxn)
            n_allo_dd += length(allo)
        end
    end

    (catalytic=length(catalytic), ress=n_ress,
     dead_end=n_de, equivalence=n_eq, dedup=n_dd,
     allosteric=n_allo, tr_equiv=n_tr,
     allosteric_dedup=n_allo_dd)
end
