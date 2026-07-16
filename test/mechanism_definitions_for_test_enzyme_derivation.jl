# ABOUTME: Shared mechanism test specifications used across the derivation tests.
# ABOUTME: Each mechanism is defined inline with its expected properties.

using LinearAlgebra
using Random

# ── MechanismTestSpec struct ────────────────────────────────────────────────

"""
Data-driven test specification for enzyme mechanisms.
Contains all expected properties for comprehensive testing.
"""
Base.@kwdef struct MechanismTestSpec
    # Core data
    name::String                          # Human-readable name for test labels
    mechanism::Any                        # EnzymeMechanism or AllostericEnzymeMechanism
    metabolite_names::Vector{Symbol}      # For param/conc generation

    # Structural expectations
    expected_n_states::Int
    expected_n_steps::Int
    expected_n_metabolites::Int

    # Constraint expectations
    expected_n_haldane_constraints::Int       # RHS references Keq (catalytic-cycle closure)
    expected_n_mirror_constraints::Int        # RHS is a single Symbol (allosteric :EqualAI rename)
    # NOTE: largely vestigial under structural naming — mostly 0, only nonzero for :EqualAI
    # reg ligands; candidate for repurposing-or-removal in the structural-naming cleanup.
    expected_n_wegscheider_constraints::Int   # RHS Expr without Keq (multi-cycle futile-cycle closure)
    expected_n_independent_params::Int        # 2*n_steps - n_constraints

    # Test configuration (optional)
    run_ode_test::Bool = true
    reference_rtol::Float64 = 1e-8
    ode_rtol::Float64 = 1e-6

    # Optional analytical rate function for extra validation (unicyclic mechanisms)
    # Signature: (all_params::NamedTuple, concs::NamedTuple) -> Float64
    analytical_rate_fn::Union{Function,Nothing} = nothing

    # Optional analytical kcat function for kcat/rescaling validation
    # Signature: (params::NamedTuple) -> Float64
    # params contains fitted_params + Keq + E_total
    analytical_kcat_fn::Union{Function,Nothing} = nothing

    # Factored form validation (optional — for mechanisms with known factored forms)
    # Expected numerator/denominator from rate_equation_string output
    expected_factored_num::Union{String,Nothing} = nothing
    expected_factored_denom::Union{String,Nothing} = nothing

    # As-written (source) catalytic-step / regulatory-site order, captured
    # before the constructor canonicalizes. Positional textbook oracles number
    # k1,k2,… (and reg1,reg2,…) by source order; the bridge in
    # `analytical_oracle_params` uses these to remap onto canonical stored
    # order. Populated only for oracle-bearing fixtures via `@..._src` macros.
    source_steps::Union{Vector{Vector{EnzymeRates.Step}},Nothing} = nothing
    source_reg_sites::Union{Vector{EnzymeRates.RegulatorySite},Nothing} = nothing
end

# Companion macros: expand exactly like `@enzyme_mechanism` /
# `@allosteric_mechanism` but ALSO return the as-written source-order step
# groups (and, for allosteric, regulatory sites), captured before the
# constructor canonicalizes. Used to populate `MechanismTestSpec.source_steps`
# / `source_reg_sites` so positional oracles bridge to canonical stored order.
macro enzyme_mechanism_src(block)
    EnzymeRates._reject_allosteric_syntax!(block)
    mech_expr, groups_expr = EnzymeRates._parse_plain_mechanism_body(block)
    esc(:(($mech_expr, $groups_expr)))
end

macro allosteric_mechanism_src(block)
    mech_expr, groups_expr, reg_sites_expr =
        EnzymeRates._parse_allosteric_mechanism_body(block)
    esc(:(($mech_expr, $groups_expr, $reg_sites_expr)))
end

# Build an AllostericEnzymeMechanism binding catalytic allosteric states to the
# catalytic steps AS WRITTEN (source order). `cm_src` is an
# `@enzyme_mechanism_src` result `(cm, source_groups)`. Routing through
# AllostericMechanism canonicalizes catalytic steps and their allosteric tags
# together, so each tag stays on its intended step (the 3-arg
# AllostericEnzymeMechanism ctor instead indexes tags by canonical group
# order). `cat_sites`/`reg_sites` use the same tuple shapes that ctor accepts:
# `cat_sites = (multiplicity, cat_allo_states)`,
# `reg_sites = ((ligands, multiplicity, ligand_states), …)`.
function allo_from_source(cm_src, cat_sites, reg_sites)
    cm, src = cm_src
    mult, cat_states = cat_sites
    sites = EnzymeRates.RegulatorySite[
        EnzymeRates.RegulatorySite(
            EnzymeRates.AllostericRegulator[
                EnzymeRates.AllostericRegulator(l) for l in ligs],
            m, collect(Symbol, states))
        for (ligs, m, states) in reg_sites]
    EnzymeRates.AllostericEnzymeMechanism(
        EnzymeRates.AllostericMechanism(
            EnzymeRates.reaction(EnzymeRates.Mechanism(cm)),
            src, collect(Symbol, cat_states), mult, sites))
end

# ── Mechanism test specifications ───────────────────────────────────────────

function build_mechanism_test_specs()
    specs = MechanismTestSpec[]

    # 1. Uni-Uni (simplest): E + S ⇌ ES ⇌ E + P
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        push!(specs, MechanismTestSpec(
            name="Uni-Uni",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P],
            expected_n_states=2,
            expected_n_steps=2,
            expected_n_metabolites=2,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=3,
            analytical_rate_fn=(p, c) -> begin
                num = p.k1f * p.k2f * c.S - p.k1r * p.k2r * c.P
                denom = p.k1r + p.k2f + p.k1f * c.S + p.k2r * c.P
                p.Et * num / denom
            end
        ))
    end

    # 2. Segel Uni Uni (replaces Three-Step Iso): E + A ⇌ EA ⇌ EP ⇌ E + P
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-8
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A
            products: P
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P)
                E(P) <--> E + P
            end
        end

        # Segel Eq. IX-8: Uni Uni steady-state rate
        function rate_uni_uni(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, Etotal) = params
            (; A, P) = concs
            num = k1f * k2f * k3f * A - k1r * k2r * k3r * P
            denom = (k1r * k2r + k1r * k3f + k2f * k3f) +
                    k1f * (k2r + k2f + k3f) * A +
                    k3r * (k1r + k2r + k2f) * P
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Uni Uni",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :P],
            expected_n_states=3,
            expected_n_steps=3,
            expected_n_metabolites=2,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=5,
            analytical_rate_fn=(p, c) -> rate_uni_uni(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 3. Segel Iso Uni Uni (new): E + A ⇌ EA ⇌ EP ⇌ F + P, F ⇌ E
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-45
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A
            products: P
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(P)
                E(P) <--> F + P
                F <--> E
            end
        end

        # Segel Eq. IX-45: Iso Uni Uni steady-state rate
        function rate_iso_uni_uni(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal) = params
            (; A, P) = concs
            num = k1f * k2f * k3f * k4f * A - k1r * k2r * k3r * k4r * P
            denom = (k4f + k4r) * (k1r * k3f + k1r * k2r + k2f * k3f) +
                    k1f * (k2f * k3f + k2f * k4f + k2r * k4f + k3f * k4f) * A +
                    k3r * (k1r * k2r + k1r * k4r + k2f * k4r + k2r * k4r) * P +
                    k1f * k3r * (k2f + k2r) * A * P
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Iso Uni Uni",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :P],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=2,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=7,
            analytical_rate_fn=(p, c) -> rate_iso_uni_uni(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 4. Segel Ordered Uni Bi (replaces Seq Uni-Bi): E + A ⇌ (EA≡EPQ) ⇌ EQ + P ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-60
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end

        # Segel Eq. IX-60: Ordered Uni Bi steady-state rate
        function rate_ordered_uni_bi(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, Etotal) = params
            (; A, P, Q) = concs
            num = k1f * k2f * k3f * A - k1r * k2r * k3r * P * Q
            denom = k3f * (k2f + k1r) +
                    k1f * (k2f + k3f) * A +
                    k1r * k2r * P +
                    k3r * (k2f + k1r) * Q +
                    k1f * k2r * A * P +
                    k2r * k3r * P * Q
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Ordered Uni Bi",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :P, :Q],
            expected_n_states=3,
            expected_n_steps=3,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=5,
            analytical_rate_fn=(p, c) -> rate_ordered_uni_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 5. Segel Ordered Bi Bi (replaces Seq Bi-Bi):
    #    E + A ⇌ EA + B ⇌ (EAB≡EPQ) ⇌ EQ + P ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-87
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end

        # Segel Eq. IX-87: Ordered Bi Bi steady-state rate
        function rate_ordered_bi_bi(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal) = params
            (; A, B, P, Q) = concs
            num = k1f * k2f * k3f * k4f * A * B - k1r * k2r * k3r * k4r * P * Q
            denom = k1r * k4f * (k2r + k3f) +
                    k1f * k4f * (k2r + k3f) * A +
                    k2f * k3f * k4f * B +
                    k1r * k2r * k3r * P +
                    k1r * k4r * (k2r + k3f) * Q +
                    k1f * k2f * (k3f + k4f) * A * B +
                    k1f * k2r * k3r * A * P +
                    k2f * k3f * k4r * B * Q +
                    k3r * k4r * (k1r + k2r) * P * Q +
                    k1f * k2f * k3r * A * B * P +
                    k2f * k3r * k4r * B * P * Q
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Ordered Bi Bi",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=7,
            analytical_rate_fn=(p, c) -> rate_ordered_bi_bi(merge(p, (Etotal=p.Et,)), c),
            analytical_kcat_fn=p -> p.k3f * p.k4f / (p.k3f + p.k4f),
        ))
    end

    # 7. Segel Ping Pong Bi Bi (replaces Ping-Pong Bi-Bi):
    #    E + A ⇌ (EA≡FP) ⇌ F + P, F + B ⇌ (FB≡EQ) ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-140
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) <--> F + P
                F + B <--> F(B)
                F(B) <--> E + Q
            end
        end

        # Segel Eq. IX-140: Ping Pong Bi Bi steady-state rate
        function rate_ping_pong_bi_bi(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, Etotal) = params
            (; A, B, P, Q) = concs
            num = k1f * k2f * k3f * k4f * A * B - k1r * k2r * k3r * k4r * P * Q
            denom = k1f * k2f * (k3r + k4f) * A +
                    k3f * k4f * (k1r + k2f) * B +
                    k1r * k2r * (k3r + k4f) * P +
                    k3r * k4r * (k1r + k2f) * Q +
                    k1f * k3f * (k2f + k4f) * A * B +
                    k1f * k2r * (k3r + k4f) * A * P +
                    k3f * k4r * (k1r + k2f) * B * Q +
                    k2r * k4r * (k1r + k3r) * P * Q
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Ping Pong Bi Bi",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=7,
            analytical_rate_fn=(p, c) ->
                rate_ping_pong_bi_bi(
                    merge(p, (Etotal=p.Et,)), c),
            analytical_kcat_fn=p -> p.k2f * p.k4f / (p.k2f + p.k4f),
        ))
    end

    # 8. Segel Ordered Ter Bi (replaces Seq Ter-Bi):
    #     E + A ⇌ EA + B ⇌ EAB + C ⇌ (EABC≡EPQ) ⇌ EQ + P ⇌ E + Q
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-195
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B, C
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) + C <--> E(A, B, C)
                E(A, B, C) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end

        # Segel Eq. IX-195: Ordered Ter Bi steady-state rate
        function rate_ordered_ter_bi(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, Etotal) = params
            (; A, B, C, P, Q) = concs
            num = k1f * k2f * k3f * k4f * k5f * A * B * C -
                  k1r * k2r * k3r * k4r * k5r * P * Q
            denom = # constant
                k1r * k2r * k5f * (k3r + k4f) +
                # single substrate/product
                k1f * k2r * k5f * (k3r + k4f) * A +
                k1r * k3f * k4f * k5f * C +
                k1r * k2r * k3r * k4r * P +
                k1r * k2r * k5r * (k3r + k4f) * Q +
                # two substrates/products
                k1f * k2f * k5f * (k3r + k4f) * A * B +
                k1f * k3f * k4f * k5f * A * C +
                k2f * k3f * k4f * k5f * B * C +
                k1f * k2r * k3r * k4r * A * P +
                k1r * k3f * k4f * k5r * C * Q +
                k4r * k5r * (k1r * k2r + k1r * k3r + k2r * k3r) * P * Q +
                # three substrates/products
                k1f * k2f * k3f * (k4f + k5f) * A * B * C +
                k1f * k2f * k3r * k4r * A * B * P +
                k2f * k3f * k4f * k5r * B * C * Q +
                k1r * k3f * k4r * k5r * C * P * Q +
                k2f * k3r * k4r * k5r * B * P * Q +
                # four substrates/products
                k1f * k2f * k3f * k4r * A * B * C * P +
                k2f * k3f * k4r * k5r * B * C * P * Q
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Ordered Ter Bi",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :C, :P, :Q],
            expected_n_states=5,
            expected_n_steps=5,
            expected_n_metabolites=5,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=9,
            analytical_rate_fn=(p, c) -> rate_ordered_ter_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 9. Segel Ordered Ter Ter (replaces Seq Ter-Ter):
    #     E + A ⇌ EA + B ⇌ EAB + C ⇌ (EABC≡EPQR) ⇌ EQR + P ⇌ ER + Q ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-261
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) + C <--> E(A, B, C)
                E(A, B, C) <--> E(Q, R) + P
                E(Q, R) <--> E(R) + Q
                E(R) <--> E + R
            end
        end

        # Segel Eq. IX-261: Ordered Ter Ter steady-state rate
        function rate_ordered_ter_ter(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
            (; A, B, C, P, Q, R) = concs
            num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
                  k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R
            denom = # constant
                k1r * k2r * k5f * k6f * (k3r + k4f) +
                # single substrate/product
                k1f * k2r * k5f * k6f * (k3r + k4f) * A +
                k1r * k3f * k4f * k5f * k6f * C +
                k1r * k2r * k3r * k4r * k6f * P +
                k1r * k2r * k5f * k6r * (k3r + k4f) * R +
                # two substrates/products
                k1f * k2f * k5f * k6f * (k3r + k4f) * A * B +
                k1f * k3f * k4f * k5f * k6f * A * C +
                k1f * k2r * k3r * k4r * k6f * A * P +
                k2f * k3f * k4f * k5f * k6f * B * C +
                k1r * k3f * k4f * k5f * k6r * C * R +
                k1r * k2r * k3r * k4r * k5r * P * Q +
                k1r * k2r * k3r * k4r * k6r * P * R +
                k1r * k2r * k5r * k6r * (k3r + k4f) * Q * R +
                # three substrates/products
                k1f * k2f * k3f * (k4f * k5f + k4f * k6f + k5f * k6f) * A * B * C +
                k1f * k2f * k3r * k4r * k6f * A * B * P +
                k1f * k2r * k3r * k4r * k5r * A * P * Q +
                k2f * k3f * k4f * k5f * k6r * B * C * R +
                k1r * k3f * k4f * k5r * k6r * C * Q * R +
                k4r * k5r * k6r * (k1r * k2r + k1r * k3r + k2r * k3r) * P * Q * R +
                # four substrates/products
                k1f * k2f * k3f * k4r * k6f * A * B * C * P +
                k1f * k2f * k3f * k4f * k5r * A * B * C * Q +
                k1f * k2f * k3r * k4r * k5r * A * B * P * Q +
                k2f * k3f * k4f * k5r * k6r * B * C * Q * R +
                k2f * k3r * k4r * k5r * k6r * B * P * Q * R +
                k1r * k3f * k4r * k5r * k6r * C * P * Q * R +
                # five substrates/products
                k1f * k2f * k3f * k4r * k5r * A * B * C * P * Q +
                k2f * k3f * k4r * k5r * k6r * B * C * P * Q * R
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Ordered Ter Ter",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :C, :P, :Q, :R],
            expected_n_states=6,
            expected_n_steps=6,
            expected_n_metabolites=6,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=11,
            analytical_rate_fn=(p, c) ->
                rate_ordered_ter_ter(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 10. Segel Bi Uni Uni Uni Ping Pong Ter Bi (new):
    #     E + A ⇌ EA + B ⇌ (EAB≡FP) ⇌ F + P, F + C ⇌ (FC≡EQ) ⇌ E + Q
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-228
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B, C
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> F + P
                F + C <--> F(C)
                F(C) <--> E + Q
            end
        end

        # Segel Eq. IX-228: Bi Uni Uni Uni Ping Pong Ter Bi steady-state rate
        function rate_bi_uni_uni_uni_ping_pong_ter_bi(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, Etotal) = params
            (; A, B, C, P, Q) = concs
            num = k1f * k2f * k3f * k4f * k5f * A * B * C -
                  k1r * k2r * k3r * k4r * k5r * P * Q
            denom = # single substrate/product
                k1r * k4f * k5f * (k2r + k3f) * C +
                k1r * k2r * k3r * (k4r + k5f) * P +
                k1r * k4r * k5r * (k2r + k3f) * Q +
                # two substrates/products
                k1f * k2f * k3f * (k4r + k5f) * A * B +
                k1f * k4f * k5f * (k2r + k3f) * A * C +
                k1f * k2r * k3r * (k4r + k5f) * A * P +
                k2f * k3f * k4f * k5f * B * C +
                k2f * k3f * k4r * k5r * B * Q +
                k1r * k4f * k5r * (k2r + k3f) * C * Q +
                k3r * k5r * (k1r * k2r + k1r * k4r + k2r * k4r) * P * Q +
                # three substrates/products
                k1f * k2f * k4f * (k3f + k5f) * A * B * C +
                k1f * k2f * k3r * (k4r + k5f) * A * B * P +
                k2f * k3f * k4f * k5r * B * C * Q +
                k2f * k3r * k4r * k5r * B * P * Q
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Bi Uni Uni Uni PP Ter Bi",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :C, :P, :Q],
            expected_n_states=5,
            expected_n_steps=5,
            expected_n_metabolites=5,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=9,
            analytical_rate_fn=(p, c) ->
                rate_bi_uni_uni_uni_ping_pong_ter_bi(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 11. Random-order Bi-Bi (branched): Two substrate binding orders converge
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B <--> E(A, B)
                E(B) + A <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        push!(specs, MechanismTestSpec(
            name="Random-order Bi-Bi",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=6,
            expected_n_steps=7,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=1,
            expected_n_independent_params=12,
            analytical_rate_fn=nothing
        ))
    end

    # 11b. Random-order, chemistry-SS Bi-Bi: two substrate binding orders,
    #      a steady-state chemistry step, single product. ODE is ground truth.
    let
        m = @enzyme_mechanism begin
            substrates: S1, S2
            products: P
            steps: begin
                E + S1 ⇌ E(S1)
                E + S2 ⇌ E(S2)
                E(S1) + S2 <--> E(S1, S2)
                E(S2) + S1 <--> E(S1, S2)
                E(S1, S2) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        push!(specs, MechanismTestSpec(
            name="Numerator: random chem-SS (RE/SS)",
            mechanism=m,
            metabolite_names=[:S1, :S2, :P],
            expected_n_states=5,
            expected_n_steps=6,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=1,
            expected_n_independent_params=7,
            analytical_rate_fn=nothing
        ))
    end

    # 11c. Catalytic isomerization RE, binding+release SS: random-order binding
    #      and product release around a rapid-equilibrium chemistry step.
    let
        m = @enzyme_mechanism begin
            substrates: S1, S2
            products: P1, P2
            steps: begin
                E + S1 ⇌ E(S1)
                E + S2 ⇌ E(S2)
                E(S1) + S2 <--> E(S1, S2)
                E(S2) + S1 <--> E(S1, S2)
                E(S1, S2) ⇌ E(P1, P2)
                E(P1, P2) <--> E(P1) + P2
                E(P1, P2) <--> E(P2) + P1
                E(P1) ⇌ E + P1
                E(P2) ⇌ E + P2
            end
        end
        push!(specs, MechanismTestSpec(
            name="Numerator: RE-chemistry",
            mechanism=m,
            metabolite_names=[:S1, :S2, :P1, :P2],
            expected_n_states=7,
            expected_n_steps=9,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=2,
            expected_n_independent_params=10,
            analytical_rate_fn=nothing
        ))
    end

    # 11d. Single-RE-segment with a redundant SS binding branch (one branch
    #      RE, one SS). ODE is ground truth.
    let
        m = @enzyme_mechanism begin
            substrates: S1, S2
            products: P
            steps: begin
                E + S1 <--> E(S1)
                E + S2 ⇌ E(S2)
                E(S1) + S2 ⇌ E(S1, S2)
                E(S2) + S1 ⇌ E(S1, S2)
                E(S1, S2) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        push!(specs, MechanismTestSpec(
            name="Numerator: redundant SS-bind",
            mechanism=m,
            metabolite_names=[:S1, :S2, :P],
            expected_n_states=5,
            expected_n_steps=6,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=1,
            expected_n_independent_params=6,
            analytical_rate_fn=nothing
        ))
    end

    # 12. Segel Bi Uni Uni Bi Ping Pong Ter Ter (new):
    #     E + A ⇌ EA + B ⇌ (EAB≡FP) ⇌ F + P, F + C ⇌ (FC≡EQR) ⇌ ER + Q ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-278
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> F + P
                F + C <--> F(C)
                F(C) <--> E(R) + Q
                E(R) <--> E + R
            end
        end

        # Segel Eq. IX-278: Bi Uni Uni Bi Ping Pong Ter Ter steady-state rate
        function rate_bi_uni_uni_bi_ping_pong_ter_ter(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
            (; A, B, C, P, Q, R) = concs
            num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
                  k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R
            denom = # single substrate/product
                k1r * k4f * k5f * k6f * (k2r + k3f) * C +
                k1r * k2r * k3r * k6f * (k4r + k5f) * P +
                # two substrates/products
                k1f * k2f * k3f * k6f * (k4r + k5f) * A * B +
                k1f * k4f * k5f * k6f * (k2r + k3f) * A * C +
                k1f * k2r * k3r * k6f * (k4r + k5f) * A * P +
                k2f * k3f * k4f * k5f * k6f * B * C +
                k1r * k4f * k5f * k6r * (k2r + k3f) * C * R +
                k1r * k2r * k3r * k4r * k5r * P * Q +
                k1r * k2r * k3r * k6r * (k4r + k5f) * P * R +
                k1r * k4r * k5r * k6r * (k2r + k3f) * Q * R +
                # three substrates/products
                k1f * k2f * k4f * (k3f * k5f + k3f * k6f + k5f * k6f) * A * B * C +
                k1f * k2f * k3r * k6f * (k4r + k5f) * A * B * P +
                k1f * k2f * k3f * k4r * k5r * A * B * Q +
                k1f * k2r * k3r * k4r * k5r * A * P * Q +
                k2f * k3f * k4f * k5f * k6r * B * C * R +
                k2f * k3f * k4r * k5r * k6r * B * Q * R +
                k1r * k4f * k5r * k6r * (k2r + k3f) * C * Q * R +
                k3r * k5r * k6r * (k1r * k2r + k1r * k4r + k2r * k4r) * P * Q * R +
                # four substrates/products
                k1f * k2f * k3f * k4f * k5r * A * B * C * Q +
                k1f * k2f * k3r * k4r * k5r * A * B * P * Q +
                k2f * k3f * k4f * k5r * k6r * B * C * Q * R +
                k2f * k3r * k4r * k5r * k6r * B * P * Q * R
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Bi Uni Uni Bi PP Ter Ter",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :C, :P, :Q, :R],
            expected_n_states=6,
            expected_n_steps=6,
            expected_n_metabolites=6,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=11,
            analytical_rate_fn=(p, c) ->
                rate_bi_uni_uni_bi_ping_pong_ter_ter(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 13. Segel Bi Bi Uni Uni Ping Pong Ter Ter (new):
    #     E + A ⇌ EA + B ⇌ (EAB≡FPQ) ⇌ FQ + P, FQ ⇌ F + Q, F + C ⇌ (FC≡ER) ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-288
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> F(Q) + P
                F(Q) <--> F + Q
                F + C <--> F(C)
                F(C) <--> E + R
            end
        end

        # Segel Eq. IX-288: Bi Bi Uni Uni Ping Pong Ter Ter steady-state rate
        function rate_bi_bi_uni_uni_ping_pong_ter_ter(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
            (; A, B, C, P, Q, R) = concs
            num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
                  k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R
            denom = # single substrate/product
                k1r * k4f * k5f * k6f * (k2r + k3f) * C +
                k1r * k4f * k5r * k6r * (k2r + k3f) * R +
                # two substrates/products
                k1f * k2f * k3f * k4f * (k5r + k6f) * A * B +
                k1f * k4f * k5f * k6f * (k2r + k3f) * A * C +
                k2f * k3f * k4f * k5f * k6f * B * C +
                k2f * k3f * k4f * k5r * k6r * B * R +
                k1r * k2r * k3r * k5f * k6f * C * P +
                k1r * k4f * k5f * k6r * (k2r + k3f) * C * R +
                k1r * k2r * k3r * k4r * (k5r + k6f) * P * Q +
                k1r * k2r * k3r * k5r * k6r * P * R +
                k1r * k4r * k5r * k6r * (k2r + k3f) * Q * R +
                # three substrates/products
                k1f * k2f * k5f * (k3f * k4f + k3f * k6f + k4f * k6f) * A * B * C +
                k1f * k2f * k3f * k4r * (k5r + k6f) * A * B * Q +
                k1f * k2r * k3r * k5f * k6f * A * C * P +
                k1f * k2r * k3r * k4r * (k5r + k6f) * A * P * Q +
                k2f * k3f * k4f * k5f * k6r * B * C * R +
                k2f * k3f * k4r * k5r * k6r * B * Q * R +
                k1r * k2r * k3r * k5f * k6r * C * P * R +
                k3r * k4r * k6r * (k1r * k2r + k1r * k5r + k2r * k5r) * P * Q * R +
                # four substrates/products
                k1f * k2f * k3r * k5f * k6f * A * B * C * P +
                k1f * k2f * k3r * k4r * (k5r + k6f) * A * B * P * Q +
                k2f * k3r * k4r * k5r * k6r * B * P * Q * R
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Bi Bi Uni Uni PP Ter Ter",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :C, :P, :Q, :R],
            expected_n_states=6,
            expected_n_steps=6,
            expected_n_metabolites=6,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=11,
            analytical_rate_fn=(p, c) ->
                rate_bi_bi_uni_uni_ping_pong_ter_ter(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 14. Segel Hexa Uni Ping Pong (new):
    #     E + A ⇌ (EA≡FP) ⇌ F + P, F + B ⇌ (FB≡GQ) ⇌ G + Q, G + C ⇌ (GC≡ER) ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-308
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                E + A <--> E(A)
                E(A) <--> F + P
                F + B <--> F(B)
                F(B) <--> G + Q
                G + C <--> G(C)
                G(C) <--> E + R
            end
        end

        # Segel Eq. IX-308: Hexa Uni Ping Pong steady-state rate
        function rate_hexa_uni_ping_pong(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, k4f, k4r, k5f, k5r, k6f, k6r, Etotal) = params
            (; A, B, C, P, Q, R) = concs
            num = k1f * k2f * k3f * k4f * k5f * k6f * A * B * C -
                  k1r * k2r * k3r * k4r * k5r * k6r * P * Q * R
            denom = # two substrates/products (9 classes, each with 2 terms → factored)
                k1f * k2f * k3f * k4f * (k5r + k6f) * A * B +
                k1f * k2f * k5f * k6f * (k3r + k4f) * A * C +
                k1f * k2f * k3r * k4r * (k5r + k6f) * A * Q +
                k3f * k4f * k5f * k6f * (k1r + k2f) * B * C +
                k3f * k4f * k5r * k6r * (k1r + k2f) * B * R +
                k1r * k2r * k5f * k6f * (k3r + k4f) * C * P +
                k1r * k2r * k3r * k4r * (k5r + k6f) * P * Q +
                k1r * k2r * k5r * k6r * (k3r + k4f) * P * R +
                k3r * k4r * k5r * k6r * (k1r + k2f) * Q * R +
                # three substrates/products (8 classes)
                k1f * k3f * k5f * (k2f * k4f + k2f * k6f + k4f * k6f) * A * B * C +
                k1f * k2f * k3f * k4r * (k5r + k6f) * A * B * Q +
                k1f * k2r * k5f * k6f * (k3r + k4f) * A * C * P +
                k1f * k2r * k3r * k4r * (k5r + k6f) * A * P * Q +
                k3f * k4f * k5f * k6r * (k1r + k2f) * B * C * R +
                k3f * k4r * k5r * k6r * (k1r + k2f) * B * Q * R +
                k1r * k2r * k5f * k6r * (k3r + k4f) * C * P * R +
                k2r * k4r * k6r * (k1r * k3r + k1r * k5r + k3r * k5r) * P * Q * R
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Hexa Uni Ping Pong",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :C, :P, :Q, :R],
            expected_n_states=6,
            expected_n_steps=6,
            expected_n_metabolites=6,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=11,
            analytical_rate_fn=(p, c) ->
                rate_hexa_uni_ping_pong(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # ── Rapid-Equilibrium (RE) Mechanisms ──────────────────────────────────────

    # 15. RE Uni-Uni: E + A ⇌_RE EA <-->_SS E + P
    #     Rapid-equilibrium substrate binding, steady-state catalysis
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A
            products: P
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E + P
            end
        end

        # Rapid-equilibrium Michaelis-Menten (K1 = Kd = [E][A]/[EA]):
        # rate = E_t * (k2f * A/K1 - k2r * P) / (1 + A/K1)
        function rate_re_uni_uni(params, concs)
            (; K1, k2f, k2r, Et) = params
            (; A, P) = concs
            num = k2f * A / K1 - k2r * P
            denom = 1.0 + A / K1
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="RE Uni-Uni",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :P],
            expected_n_states=2,
            expected_n_steps=2,
            expected_n_metabolites=2,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=2,
            analytical_rate_fn=(p, c) -> rate_re_uni_uni(merge(p, (Et=p.Et,)), c),
            analytical_kcat_fn=p -> p.k2f,
        ))
    end

    # 16. RE Ordered Bi-Bi: substrate binding is RE, catalysis and product release are SS
    #     E + A ⇌_RE EA, EA + B ⇌_RE EAB, (EAB≡EPQ) <-->_SS EQ + P, EQ <-->_SS E + Q
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end

        # Groups: {E, EA, EAB} (RE group) and {EQ} (singleton)
        # K1, K2 are Kd (dissociation constants):
        #   σ = 1 + A/K1 + A*B/(K1*K2)
        # num = k3f*k4f*A*B/(K1*K2) - k3r*k4r*P*Q
        # denom = (1+A/K1+A*B/(K1*K2))*(k3r*P+k4f) + k3f*A*B/(K1*K2) + k4r*Q
        function rate_re_ordered_bi_bi(params, concs)
            (; K1, K2, k3f, k3r, k4f, k4r, Et) = params
            (; A, B, P, Q) = concs
            num = k3f * k4f * A * B / (K1 * K2) - k3r * k4r * P * Q
            sigma1 = 1.0 + A / K1 + A * B / (K1 * K2)
            R12 = k3f * A * B / (K1 * K2) + k4r * Q
            R21 = k3r * P + k4f
            denom = sigma1 * R21 + R12
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="RE Ordered Bi-Bi",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=5,
            analytical_rate_fn=(p, c) -> rate_re_ordered_bi_bi(merge(p, (Et=p.Et,)), c)
        ))
    end

    # 17. RE Random Bi-Bi: binding to free enzyme is RE, other steps are SS
    #     E + A ⇌_RE EA, E + B ⇌_RE EB, EA + B <-->_SS EAB, EB + A <-->_SS EAB,
    #     EAB <-->_SS EPQ, EPQ <-->_SS EQ + P, EQ <-->_SS E + Q
    let
        m, src = @enzyme_mechanism_src begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E + B ⇌ E(B)
                E(A) + B <--> E(A, B)
                E(B) + A <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        push!(specs, MechanismTestSpec(
            name="RE Random Bi-Bi",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=6,
            expected_n_steps=7,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=1,
            expected_n_independent_params=10,
            analytical_rate_fn=nothing
        ))
    end

    # ── Classical Inhibitor/Activator Mechanisms (factored form tests) ────────

    # 18. Competitive inhibitor: E + S ⇌ ES, ES ⇌ EP (SS), EP ⇌ E + P, E + R ⇌ ER
    #     Dead-end inhibitor R binds free enzyme only.
    #     No Cartesian product structure → flat sum denominator.
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                E + S ⇌ E(S)      # K1
                E(S) <--> E(P)    # k2f, k2r (SS)
                E(P) ⇌ E + P      # K3
                E + R ⇌ E(R)      # K4
            end
        end

        # v = Et * (k2f*S/K1 - k2r*P/K3) / (1 + S/K1 + P/K3 + R/K4)
        function rate_competitive_inh(params, concs)
            (; K1, k2f, k2r, K3, K4, Et) = params
            (; S, P, R) = concs
            num = k2f * S / K1 - k2r * P / K3
            denom = 1.0 + S / K1 + P / K3 + R / K4
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Competitive Inhibitor",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :R],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=4,
            analytical_rate_fn=(p, c) -> rate_competitive_inh(
                merge(p, (Et=p.Et,)), c),
            analytical_kcat_fn=p -> p.k2f,
            # Textbook: flat sum denominator (no Cartesian product structure)
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E",
            expected_factored_denom=
            "1 + P / K_P_E + R / K_Rinh_E + S / K_S_E",
        ))
    end

    # 19. Non-competitive inhibitor: R binds both free E and ES with same K
    #     Forms: E, E_S, E_P, E_R, E_S_R; SS: E_S↔E_P; K5=K4
    #     Denom factors as (1+R/K4)*(1+S/K1) + P/K3
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                E + S ⇌ E(S)        # K1
                E(S) <--> E(P)      # k2f, k2r (SS)
                E + P ⇌ E(P)        # K3
                (E + R ⇌ E(R), E(S) + R ⇌ E(S, R))  # K4 = K5 (R binding shared)
                E(R) + S ⇌ E(S, R)  # K6
            end
        end

        # v = Et * (k2f*S/K1 - k2r*P/K3) / ((1+R/K4)*(1+S/K1) + P/K3)
        function rate_noncompetitive_inh(params, concs)
            (; K1, k2f, k2r, K3, K4, Et) = params
            (; S, P, R) = concs
            num = k2f * S / K1 - k2r * P / K3
            denom = (1.0 + R / K4) * (1.0 + S / K1) + P / K3
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Non-competitive Inhibitor",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :R],
            expected_n_states=5,
            expected_n_steps=6,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=4,
            analytical_rate_fn=(p, c) -> rate_noncompetitive_inh(
                merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E",
            expected_factored_denom=
            "1 + P / K_P_E + R / K_Rinh_E + S / K_S_E + R * S / (K_Rinh_E * K_S_E)",
        ))
    end

    # 20. Uncompetitive inhibitor: R binds ES only (not free E)
    #     Forms: E, E_S, E_P, E_S_R; SS: E_S↔E_P; No extra constraints
    #     Denom: 1 + P/K3 + S/K1*(1+R/K4)
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                E + S ⇌ E(S)        # K1
                E(S) <--> E(P)      # k2f, k2r (SS)
                E + P ⇌ E(P)        # K3
                E(S) + R ⇌ E(S, R)  # K4
            end
        end

        # v = Et * (k2f*S/K1 - k2r*P/K3) / (1 + P/K3 + S/K1*(1+R/K4))
        function rate_uncompetitive_inh(params, concs)
            (; K1, k2f, k2r, K3, K4, Et) = params
            (; S, P, R) = concs
            num = k2f * S / K1 - k2r * P / K3
            denom = 1.0 + P / K3 + (S / K1) * (1.0 + R / K4)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Uncompetitive Inhibitor",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :R],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=4,
            analytical_rate_fn=(p, c) -> rate_uncompetitive_inh(
                merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E",
            expected_factored_denom=
            "1 + P / K_P_E + S / K_S_E + R * S / (K_Rinh_ES * K_S_E)",
        ))
    end

    # 21. Essential activator: R must bind before S can bind
    #     Forms: E, E_R, E_S_R, E_P_R; SS: E_S_R↔E_P_R
    #     Num: R/K4 * (k2f*S/K1 - k2r*P/K3)
    #     Denom: 1 + R/K4*(1+S/K1+P/K3)
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                E(R) + S ⇌ E(S, R)    # K1
                E(S, R) <--> E(P, R)  # k2f, k2r (SS)
                E(R) + P ⇌ E(P, R)    # K3
                E + R ⇌ E(R)          # K4
            end
        end

        # v = Et * R/K4 * (k2f*S/K1 - k2r*P/K3) / (1 + R/K4*(1+S/K1+P/K3))
        function rate_essential_activator(params, concs)
            (; K1, k2f, k2r, K3, K4, Et) = params
            (; S, P, R) = concs
            num = (R / K4) * (k2f * S / K1 - k2r * P / K3)
            denom = 1.0 + (R / K4) * (1.0 + S / K1 + P / K3)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Essential Activator",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :R],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=4,
            analytical_rate_fn=(p, c) -> rate_essential_activator(
                merge(p, (Et=p.Et,)), c),
            analytical_kcat_fn=p -> p.k2f,
            expected_factored_num=
            "k_ERinhS_to_EPRinh * R * S / (K_Rinh_E * K_S_ERinh) - k_EPRinh_to_ERinhS * P * R / (K_P_ERinh * K_Rinh_E)",
            expected_factored_denom=
            "1 + R / K_Rinh_E + P * R / (K_P_ERinh * K_Rinh_E) + R * S / (K_Rinh_E * K_S_ERinh)",
        ))
    end

    # 22. Non-essential activator (general modifier):
    #     R modifies catalysis but isn't required. Two parallel SS cycles.
    #     Forms: E, E_S, E_P, E_R, E_S_R, E_P_R; SS: E_S↔E_P, E_S_R↔E_P_R
    #     K8=K7, K9=K7 (R binding independent of S/P)
    #     K4=K1 and K6=K3 are implied by Wegscheider relations
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                E + S ⇌ E(S)            # K1
                E(S) <--> E(P)          # k2f, k2r (SS)
                E + P ⇌ E(P)            # K3
                E(R) + S ⇌ E(S, R)      # K4
                E(S, R) <--> E(P, R)    # k5f, k5r (SS)
                E(R) + P ⇌ E(P, R)      # K6
                (E + R ⇌ E(R),          # K7 = K8 = K9 (R binding shared)
                 E(S) + R ⇌ E(S, R),
                 E(P) + R ⇌ E(P, R))
            end
        end

        # v = Et * ((k2f*S/K1-k2r*P/K3) + R/K7*(k5f*S/K1-k5r*P/K3))
        #         / ((1+S/K1+P/K3)*(1+R/K7))
        function rate_nonessential_activator(params, concs)
            (; K1, k2f, k2r, K3, k5f, k5r, K7, Et) = params
            (; S, P, R) = concs
            num = (k2f * S / K1 - k2r * P / K3) +
                  (R / K7) * (k5f * S / K1 - k5r * P / K3)
            denom = (1.0 + S / K1 + P / K3) * (1.0 + R / K7)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Non-essential Activator",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :R],
            expected_n_states=6,
            expected_n_steps=9,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=2,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=5,
            analytical_rate_fn=(p, c) -> rate_nonessential_activator(
                merge(p, (Et=p.Et,)), c),
            analytical_kcat_fn=p -> max(p.k2f, p.k5f),
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E + k_ERinhS_to_EPRinh * R * S / (K_Rinh_E * K_S_E) - (k_EP_to_ES * P / K_P_E + k_EPRinh_to_ERinhS * P * R / (K_P_E * K_Rinh_E))",
            expected_factored_denom=
            "1 + P / K_P_E + R / K_Rinh_E + S / K_S_E + P * R / (K_P_E * K_Rinh_E) + R * S / (K_Rinh_E * K_S_E)",
        ))
    end

    # 23. Non-essential activator + competitive inhibitor:
    #     Combines non-essential activation with competitive inhibition.
    #     A modifies catalysis but isn't required (binds E, E_S, E_P with same K).
    #     I binds only free E (competitive dead-end).
    #     Forms: E, E_S, E_P, E_A, E_S_A, E_P_A, E_I
    #     SS: E_S↔E_P, E_S_A↔E_P_A
    #     Constraints: K8=K7, K9=K7 (A binding K independent of S/P)
    #     Wegscheider gives K4=K1, K6=K3
    #     Denom: (1+S/K1+P/K3)*(1+A/K7) + I/K10
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: A, I
            steps: begin
                E + S ⇌ E(S)            # K1
                E(S) <--> E(P)          # k2f, k2r (SS)
                E + P ⇌ E(P)            # K3
                E(A) + S ⇌ E(A, S)      # K4
                E(A, S) <--> E(A, P)    # k5f, k5r (SS)
                E(A) + P ⇌ E(A, P)      # K6
                (E + A ⇌ E(A),          # K7 = K8 = K9 (A binding shared)
                 E(S) + A ⇌ E(A, S),
                 E(P) + A ⇌ E(A, P))
                E + I ⇌ E(I)            # K10
            end
        end

        # v = Et * [(k2f*S/K1 - k2r*P/K3) + (A/K7)*(k5f*S/K1 - k5r*P/K3)]
        #         / [(1+S/K1+P/K3)*(1+A/K7) + I/K10]
        function rate_activator_inhibitor(params, concs)
            (; K1, k2f, k2r, K3, k5f, k5r, K7, K10, Et) = params
            (; S, P, A, I) = concs
            num = (k2f * S / K1 - k2r * P / K3) +
                  (A / K7) * (k5f * S / K1 - k5r * P / K3)
            denom = (1.0 + S / K1 + P / K3) * (1.0 + A / K7) +
                    I / K10
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Activator + Competitive Inhibitor",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :A, :I],
            expected_n_states=7,
            expected_n_steps=10,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=2,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=6,
            analytical_rate_fn=(p, c) ->
                rate_activator_inhibitor(
                    merge(p, (Et=p.Et,)), c),
            # Denom has both multiplicative (activator) and additive
            # (inhibitor) structure
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E + k_EAinhS_to_EAinhP * A * S / (K_Ainh_E * K_S_E) - (k_EP_to_ES * P / K_P_E + k_EAinhP_to_EAinhS * A * P / (K_Ainh_E * K_P_E))",
            expected_factored_denom=
            "1 + A / K_Ainh_E + I / K_Iinh_E + P / K_P_E + S / K_S_E + A * P / (K_Ainh_E * K_P_E) + A * S / (K_Ainh_E * K_S_E)",
        ))
    end

    # 24. MWC (Monod-Wyman-Changeux) dimer with R/T conformational states
    #     Two conformational states: R (relaxed, active) and T (tense, inactive).
    #     Two identical substrate binding sites per state. Only R state is
    #     catalytically active (SS isomerization S→P at each site).
    #     R ⇌ T isomerization at base form only (classical MWC).
    #     R forms: R_00, R_S0, R_0S, R_SS, R_P0, R_0P, R_SP, R_PS, R_PP (9)
    #     T forms: T_00, T_S0, T_0S, T_SS, T_P0, T_0P, T_SP, T_PS, T_PP (9)
    #     Total: 18 forms, 23 steps
    #     Denominator: (1 + S/Ks_R + P/Kp_R)^2 + L*(1 + S/Ks_T + P/Kp_T)^2
    #     Numerator: (
    #                   (k_cat_R*S/Ks_R - k_rev_R*P/Kp_R)*(1 + S/Ks_R + P/Kp_R) +
    #                   L*(k_cat_T*S/Ks_T - k_rev_T*P/Kp_T)*(1 + S/Ks_T + P/Kp_T)
    #                 )
    # The explicit flat homodimer fixture is intentionally omitted; the
    # AllostericEnzymeMechanism form below is the canonical representation.

    # 24B. AllostericEnzymeMechanism MWC Dimer
    #      2 catalytic sites × 2 conformations (R/T). No explicit equality constraints
    #      needed — symmetric subunits are captured by the site multiplicity.
    #      Conformational equilibrium: L (= K37 in the EnzymeMechanism above).
    let
        m, src, src_reg = @allosteric_mechanism_src begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)    :: NonequalAI
                E + P ⇌ E(P)    :: NonequalAI
                E(S) <--> E(P)  :: NonequalAI
            end
        end

        function rate_mwc_dimer_oligo(params, concs)
            (; K1, K2, k3f, k3r, K1_T, K2_T, k3f_T, k3r_T, L, Et) = params
            (; S, P) = concs
            r_flux   = k3f * S / K1 - k3r * P / K2
            t_flux   = k3f_T * S / K1_T - k3r_T * P / K2_T
            r_factor = 1.0 + S / K1 + P / K2
            t_factor = 1.0 + S / K1_T + P / K2_T
            return Et * (r_flux * r_factor + L * t_flux * t_factor) /
                       (r_factor^2 + L * t_factor^2)
        end

        push!(specs, MechanismTestSpec(
            name="MWC Dimer [AllostericEnzymeMechanism]",
            mechanism=m,
            source_steps=src,
            source_reg_sites=src_reg,
            metabolite_names=[:S, :P],
            expected_n_states=3,          # catalytic subunit: E_c, E_S, E_P
            expected_n_steps=3,           # 2 RE + 1 SS per subunit
            expected_n_metabolites=2,
            expected_n_haldane_constraints=2,         # k3r per conformation × 2
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=7,
            run_ode_test=false,
            analytical_rate_fn=rate_mwc_dimer_oligo,
            expected_factored_num=
            "(k_A_ES_to_EP * S / K_A_S_E - k_A_EP_to_ES * P / K_A_P_E) * (1 + P / K_A_P_E + S / K_A_S_E)" *
            " + L * (S * k_I_ES_to_EP / K_I_S_E - P * k_I_EP_to_ES / K_I_P_E) * (1 + P / K_I_P_E + S / K_I_S_E)",
            expected_factored_denom=
            "(1 + P / K_A_P_E + S / K_A_S_E) ^ 2 + L * (1 + P / K_I_P_E + S / K_I_S_E) ^ 2",
        ))
    end

    # ── Edge-case factoring tests ─────────────────────────────────────────────
    # These mechanisms test factoring patterns not covered by the classical
    # inhibitor/activator mechanisms above.

    # The explicit flat homodimer inhibitor fixture is intentionally omitted;
    # the AllostericEnzymeMechanism form below is canonical.

    # 25B. AllostericEnzymeMechanism Homodimer + Non-competitive Inhibitor
    #      I binds all enzyme forms independently with the same Ki (enzyme-level).
    #      sigma = Q_cat^2 * (1 + I/K_I_reg1)  (multiplicative factor).
    let
        m, src, src_reg = @allosteric_mechanism_src begin
            substrates: S
            products: P
            allosteric_regulators: I::NonequalAI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)    :: NonequalAI
                E + P ⇌ E(P)    :: NonequalAI
                E(S) <--> E(P)  :: NonequalAI
            end
            regulatory_site(multiplicity = 1): begin
                ligands: I
            end
        end

        function rate_homodimer_noncomp_inh_oligo(params, concs)
            (; K1, K2, k3f, k3r, K1_T, K2_T, k3f_T, k3r_T,
               L, K_I_reg1, K_I_T_reg1, Et) = params
            (; S, P, I) = concs
            r_flux   = k3f * S / K1 - k3r * P / K2
            t_flux   = k3f_T * S / K1_T - k3r_T * P / K2_T
            r_factor = 1.0 + S / K1 + P / K2
            t_factor = 1.0 + S / K1_T + P / K2_T
            num   = r_flux * r_factor * (1.0 + I / K_I_reg1) +
                    L * t_flux * t_factor * (1.0 + I / K_I_T_reg1)
            denom = r_factor^2 * (1.0 + I / K_I_reg1) +
                    L * t_factor^2 * (1.0 + I / K_I_T_reg1)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Homodimer + Non-competitive Inhibitor [AllostericEnzymeMechanism]",
            mechanism=m,
            source_steps=src,
            source_reg_sites=src_reg,
            metabolite_names=[:S, :P, :I],
            expected_n_states=3,
            expected_n_steps=3,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=2,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=9,
            run_ode_test=false,
            analytical_rate_fn=rate_homodimer_noncomp_inh_oligo,
            expected_factored_num=
            "(k_A_ES_to_EP * S / K_A_S_E - k_A_EP_to_ES * P / K_A_P_E) * (1 + P / K_A_P_E + S / K_A_S_E) * (1 + I / K_A_Ireg)" *
            " + L * (S * k_I_ES_to_EP / K_I_S_E - P * k_I_EP_to_ES / K_I_P_E) * (1 + P / K_I_P_E + S / K_I_S_E) * (1 + I / K_I_Ireg)",
            expected_factored_denom=
            "(1 + P / K_A_P_E + S / K_A_S_E) ^ 2 * (1 + I / K_A_Ireg)" *
            " + L * (1 + P / K_I_P_E + S / K_I_S_E) ^ 2 * (1 + I / K_I_Ireg)",
        ))
    end

    # The explicit flat MWC inhibitor fixture is intentionally omitted; the
    # AllostericEnzymeMechanism form below is canonical.

    # 26B. AllostericEnzymeMechanism MWC Dimer + Independent Inhibitor
    #      The Wegscheider constraint K80 = K47*K37/K38 (R_00I ⇌ T_00I equilibrium)
    #      is automatically satisfied by the conformational assembly formula — no
    #      explicit constraint needed in the AllostericEnzymeMechanism DSL.
    let
        m, src, src_reg = @allosteric_mechanism_src begin
            substrates: S
            products: P
            allosteric_regulators: I::NonequalAI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)    :: NonequalAI
                E + P ⇌ E(P)    :: NonequalAI
                E(S) <--> E(P)  :: NonequalAI
            end
            regulatory_site(multiplicity = 1): begin
                ligands: I
            end
        end

        function rate_mwc_dimer_inh_oligo(params, concs)
            (; K1, K2, k3f, k3r, K1_T, K2_T, k3f_T, k3r_T,
               L, K_I_reg1, K_I_T_reg1, Et) = params
            (; S, P, I) = concs
            r_flux   = k3f * S / K1 - k3r * P / K2
            t_flux   = k3f_T * S / K1_T - k3r_T * P / K2_T
            r_factor = 1.0 + S / K1 + P / K2
            t_factor = 1.0 + S / K1_T + P / K2_T
            num   = r_flux * r_factor * (1.0 + I / K_I_reg1) +
                    L * t_flux * t_factor * (1.0 + I / K_I_T_reg1)
            denom = r_factor^2 * (1.0 + I / K_I_reg1) +
                    L * t_factor^2 * (1.0 + I / K_I_T_reg1)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="MWC Dimer + Independent Inhibitor [AllostericEnzymeMechanism]",
            mechanism=m,
            source_steps=src,
            source_reg_sites=src_reg,
            metabolite_names=[:S, :P, :I],
            expected_n_states=3,
            expected_n_steps=3,
            expected_n_metabolites=3,
            expected_n_haldane_constraints=2,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,  # thermodynamic consistency is automatic
            expected_n_independent_params=9,
            run_ode_test=false,
            analytical_rate_fn=rate_mwc_dimer_inh_oligo,
            expected_factored_num=
            "(k_A_ES_to_EP * S / K_A_S_E - k_A_EP_to_ES * P / K_A_P_E) * (1 + P / K_A_P_E + S / K_A_S_E) * (1 + I / K_A_Ireg)" *
            " + L * (S * k_I_ES_to_EP / K_I_S_E - P * k_I_EP_to_ES / K_I_P_E) * (1 + P / K_I_P_E + S / K_I_S_E) * (1 + I / K_I_Ireg)",
            expected_factored_denom=
            "(1 + P / K_A_P_E + S / K_A_S_E) ^ 2 * (1 + I / K_A_Ireg)" *
            " + L * (1 + P / K_I_P_E + S / K_I_S_E) ^ 2 * (1 + I / K_I_Ireg)",
        ))
    end

    # 27. Two Competitive Inhibitors (monomer)
    #     Tests: multiple additive dead-end terms (flat sum denominator).
    #     I1 and I2 both bind only free E (competitive with S/P and each other).
    #     Denom: 1 + S/K1 + P/K3 + I1/K4 + I2/K5
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
                E + I1 ⇌ E(I1)
                E + I2 ⇌ E(I2)
            end
        end

        function rate_two_comp_inh(p, c)
            (; K1, k2f, k2r, K3, K4, K5, Et) = p
            (; S, P, I1, I2) = c
            num = k2f * S / K1 - k2r * P / K3
            denom = 1.0 + S / K1 + P / K3 + I1 / K4 + I2 / K5
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Two Competitive Inhibitors",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=5,
            expected_n_steps=5,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=5,
            analytical_rate_fn=(p, c) ->
                rate_two_comp_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E",
            expected_factored_denom=
            "1 + I1 / K_I1inh_E + I2 / K_I2inh_E + P / K_P_E + S / K_S_E",
        ))
    end

    # 28. Two Non-competitive Inhibitors (monomer)
    #     Tests: triple multiplicative product factoring.
    #     I1 and I2 bind independently to all forms (E, ES, EP) at separate sites.
    #     Denom: (1+S/K1+P/K3) * (1+I1/K4) * (1+I2/K7)
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                # S binding (shared K1)
                (E + S ⇌ E(S),
                 E(I1) + S ⇌ E(S, I1),
                 E(I2) + S ⇌ E(S, I2),
                 E(I1, I2) + S ⇌ E(S, I1, I2))
                E(S) <--> E(P)
                # P binding (shared K3)
                (E + P ⇌ E(P),
                 E(I1) + P ⇌ E(P, I1),
                 E(I2) + P ⇌ E(P, I2),
                 E(I1, I2) + P ⇌ E(P, I1, I2))
                # I1 binding (shared K4)
                (E + I1 ⇌ E(I1),
                 E(S) + I1 ⇌ E(S, I1),
                 E(P) + I1 ⇌ E(P, I1),
                 E(I2) + I1 ⇌ E(I1, I2),
                 E(S, I2) + I1 ⇌ E(S, I1, I2),
                 E(P, I2) + I1 ⇌ E(P, I1, I2))
                # I2 binding (shared K7)
                (E + I2 ⇌ E(I2),
                 E(S) + I2 ⇌ E(S, I2),
                 E(P) + I2 ⇌ E(P, I2),
                 E(I1) + I2 ⇌ E(I1, I2),
                 E(S, I1) + I2 ⇌ E(S, I1, I2),
                 E(P, I1) + I2 ⇌ E(P, I1, I2))
            end
        end

        # Param names use kinetic-group representative-step indices:
        # K1=S-binding, k5f=iso, K6=P-binding, K10=I1-binding, K16=I2-binding.
        function rate_two_noncomp_inh(p, c)
            (; K1, k5f, k5r, K6, K10, K16, Et) = p
            (; S, P, I1, I2) = c
            num = k5f * S / K1 - k5r * P / K6
            denom = (1.0 + S / K1 + P / K6) *
                    (1.0 + I1 / K10) * (1.0 + I2 / K16)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Two Non-competitive Inhibitors",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=12,
            expected_n_steps=21,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=5,
            analytical_rate_fn=(p, c) ->
                rate_two_noncomp_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E",
            expected_factored_denom=
            "1 + I1 / K_I1inh_E + I2 / K_I2inh_E + P / K_P_E + S / K_S_E + I1 * I2 / (K_I1inh_E * K_I2inh_E) + I1 * P / (K_I1inh_E * K_P_E) + I1 * S / (K_I1inh_E * K_S_E) + I2 * P / (K_I2inh_E * K_P_E) + I2 * S / (K_I2inh_E * K_S_E) + I1 * I2 * P / (K_I1inh_E * K_I2inh_E * K_P_E) + I1 * I2 * S / (K_I1inh_E * K_I2inh_E * K_S_E)",
        ))
    end

    # 29. Non-competitive + Competitive Inhibitor (monomer)
    #     Tests: multiplicative product + additive dead-end term.
    #     I1 binds all forms independently (non-competitive).
    #     I2 binds only free E (competitive).
    #     Denom: (1+S/K1+P/K3)*(1+I1/K4) + I2/K9
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                # S binding (shared K1)
                (E + S ⇌ E(S), E(I1) + S ⇌ E(S, I1))
                E(S) <--> E(P)
                # P binding (shared K3)
                (E + P ⇌ E(P), E(I1) + P ⇌ E(P, I1))
                # I1 binding (shared K4)
                (E + I1 ⇌ E(I1),
                 E(S) + I1 ⇌ E(S, I1),
                 E(P) + I1 ⇌ E(P, I1))
                E + I2 ⇌ E(I2)
            end
        end

        # Param names use kinetic-group representative-step indices:
        # K1=S, k3f=iso, K4=P, K6=I1-binding, K9=I2-dead-end.
        function rate_noncomp_comp_inh(p, c)
            (; K1, k3f, k3r, K4, K6, K9, Et) = p
            (; S, P, I1, I2) = c
            num = k3f * S / K1 - k3r * P / K4
            denom = (1.0 + S / K1 + P / K4) *
                    (1.0 + I1 / K6) + I2 / K9
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Non-competitive + Competitive Inhibitor",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=7,
            expected_n_steps=9,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=5,
            analytical_rate_fn=(p, c) ->
                rate_noncomp_comp_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E",
            expected_factored_denom=
            "1 + I1 / K_I1inh_E + I2 / K_I2inh_E + P / K_P_E + S / K_S_E + I1 * P / (K_I1inh_E * K_P_E) + I1 * S / (K_I1inh_E * K_S_E)",
        ))
    end

    # 30. Uncompetitive + Competitive Inhibitor (monomer)
    #     Tests: mixed additive structure with nested multiplicative term.
    #     I1 binds only ES (uncompetitive). I2 binds only free E (competitive).
    #     Denom: 1 + I2/K5 + P/K3 + (S/K1)*(1+I1/K4)
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
                E(S) + I1 ⇌ E(S, I1)
                E + I2 ⇌ E(I2)
            end
        end

        function rate_uncomp_comp_inh(p, c)
            (; K1, k2f, k2r, K3, K4, K5, Et) = p
            (; S, P, I1, I2) = c
            num = k2f * S / K1 - k2r * P / K3
            denom = 1.0 + I2 / K5 + P / K3 +
                    (S / K1) * (1.0 + I1 / K4)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Uncompetitive + Competitive Inhibitor",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=5,
            expected_n_steps=5,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=5,
            analytical_rate_fn=(p, c) ->
                rate_uncomp_comp_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E",
            expected_factored_denom=
            "1 + I2 / K_I2inh_E + P / K_P_E + S / K_S_E + I1 * S / (K_I1inh_ES * K_S_E)",
        ))
    end

    # 31. Two Same-site (Competing) Non-competitive Inhibitors (monomer)
    #     Tests: multiplicative product with flat-sum allosteric factor.
    #     I1 and I2 compete for the same allosteric site X (mutually exclusive).
    #     Both bind independently of S/P (non-competitive).
    #     Denom: (1+S/K1+P/K3) * (1+I1/K4+I2/K9)
    let
        m, src = @enzyme_mechanism_src begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                # S binding (shared K1)
                (E + S ⇌ E(S),
                 E(I1) + S ⇌ E(S, I1),
                 E(I2) + S ⇌ E(S, I2))
                E(S) <--> E(P)
                # P binding (shared K3)
                (E + P ⇌ E(P),
                 E(I1) + P ⇌ E(P, I1),
                 E(I2) + P ⇌ E(P, I2))
                # I1 binding (shared K4)
                (E + I1 ⇌ E(I1),
                 E(S) + I1 ⇌ E(S, I1),
                 E(P) + I1 ⇌ E(P, I1))
                # I2 binding (shared K9)
                (E + I2 ⇌ E(I2),
                 E(S) + I2 ⇌ E(S, I2),
                 E(P) + I2 ⇌ E(P, I2))
            end
        end

        # Param names use kinetic-group representative-step indices:
        # K1=S, k4f=iso, K5=P, K8=I1, K11=I2.
        function rate_two_samesite_inh(p, c)
            (; K1, k4f, k4r, K5, K8, K11, Et) = p
            (; S, P, I1, I2) = c
            num = k4f * S / K1 - k4r * P / K5
            denom = (1.0 + S / K1 + P / K5) *
                    (1.0 + I1 / K8 + I2 / K11)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Two Same-site Inhibitors",
            mechanism=m,
            source_steps=src,
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=9,
            expected_n_steps=13,
            expected_n_metabolites=4,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=5,
            analytical_rate_fn=(p, c) ->
                rate_two_samesite_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E",
            expected_factored_denom=
            "1 + I1 / K_I1inh_E + I2 / K_I2inh_E + P / K_P_E + S / K_S_E + I1 * P / (K_I1inh_E * K_P_E) + I1 * S / (K_I1inh_E * K_S_E) + I2 * P / (K_I2inh_E * K_P_E) + I2 * S / (K_I2inh_E * K_S_E)",
        ))
    end

    # 32. MWC Tetramer — Random-Order Bi-Bi RE + Two Allosteric Sites
    #
    # Reaction: S1 + S2 ⇌ P1 + P2  (Bi-Bi)
    #
    # Catalytic site (×4 subunits, 2 conformations):
    #   Site A: S1 and P1 bind competitively (Kd = K1, K2 in R state)
    #   Site B: S2 and P2 bind competitively (Kd = K3, K4 in R state)
    #   Sites A and B are independent — cross-occupancy (S1+P2, P1+S2) is allowed.
    #   All binding steps are RE; only the ternary complex interconversion is SS:
    #     E_S1S2 <--> E_P1P2  (k13f, k13r; k13r derived from shared Keq via Haldane)
    #   Wegscheider constraints enforce site independence:
    #     K5=K3, K6=K1  (two routes to E_S1S2)
    #     K7=K4, K8=K1  (two routes to E_S1P2)
    #     K9=K3, K10=K2 (two routes to E_P1S2)
    #     K11=K4, K12=K2 (two routes to E_P1P2)
    #
    # Reg site type 1 (×4 per tetramer): R1 and R2 bind competitively
    # Reg site type 2 (×4 per tetramer): R3 binds exclusively
    #
    # Analytical rate equation (all Kd = dissociation constants):
    #
    #   R-state catalytic site:
    #     Q_A_R   = 1 + S1/K1 + P1/K2
    #     Q_B_R   = 1 + S2/K3 + P2/K4
    #     Q_cat_R = Q_A_R * Q_B_R          (site A ⊗ site B independence)
    #     N_cat_R = k13f*S1*S2/(K1*K3) - k13r*P1*P2/(K2*K4)
    #
    #   T-state: same with _T suffix on every K and k parameter.
    #
    #   Regulatory partition functions (star topology → direct sum):
    #     Q_reg1_R = 1 + R1/K_R1_reg1 + R2/K_R2_reg1
    #     Q_reg2_R = 1 + R3/K_R3_reg2
    #     (T-state: K_R1_T_reg1, K_R2_T_reg1, K_R3_T_reg2)
    #
    #   MWC assembly (n_cat=4, n_reg1=4, n_reg2=4):
    #     Q_R = Q_cat_R^4 * Q_reg1_R^4 * Q_reg2_R^4
    #     Q_T = Q_cat_T^4 * Q_reg1_T^4 * Q_reg2_T^4
    #     Z   = Q_R + L * Q_T
    #     Rate = Et * (N_cat_R * Q_cat_R^3 * Q_reg1_R^4 * Q_reg2_R^4
    #                + L * N_cat_T * Q_cat_T^3 * Q_reg1_T^4 * Q_reg2_T^4) / Z
    #
    # Independent parameters (17, excluding Keq and Et):
    #   R-state catalytic: K1, K2, K3, K4, k13f
    #   T-state catalytic: K1_T, K2_T, K3_T, K4_T, k13f_T
    #   Reg site 1 R: K_R1_reg1, K_R2_reg1
    #   Reg site 1 T: K_R1_T_reg1, K_R2_T_reg1
    #   Reg site 2 R: K_R3_reg2
    #   Reg site 2 T: K_R3_T_reg2
    #   Conformational equilibrium: L
    let
        m, src, src_reg = @allosteric_mechanism_src begin
            substrates: S1, S2
            products: P1, P2
            allosteric_regulators: R1::NonequalAI, R2::NonequalAI, R3::NonequalAI
            catalytic_multiplicity: 4
            catalytic_steps: begin
                # S1 binding (shared K)
                (E + S1 ⇌ E(S1),
                 E(S2) + S1 ⇌ E(S1, S2),
                 E(P2) + S1 ⇌ E(S1, P2)) :: NonequalAI
                # P1 binding (shared K)
                (E + P1 ⇌ E(P1),
                 E(S2) + P1 ⇌ E(P1, S2),
                 E(P2) + P1 ⇌ E(P1, P2)) :: NonequalAI
                # S2 binding (shared K)
                (E + S2 ⇌ E(S2),
                 E(S1) + S2 ⇌ E(S1, S2),
                 E(P1) + S2 ⇌ E(P1, S2)) :: NonequalAI
                # P2 binding (shared K)
                (E + P2 ⇌ E(P2),
                 E(S1) + P2 ⇌ E(S1, P2),
                 E(P1) + P2 ⇌ E(P1, P2)) :: NonequalAI
                E(S1, S2) <--> E(P1, P2) :: NonequalAI
            end
            regulatory_site(multiplicity = 4): begin
                ligands: R1, R2
            end
            regulatory_site(multiplicity = 4): begin
                ligands: R3
            end
        end

        # Parameter naming convention:
        #   R-state catalytic Kd: K1=Kd(S1), K2=Kd(P1), K3=Kd(S2), K4=Kd(P2)
        #   R-state SS rate:      k13f (k13r is Haldane-derived)
        #   T-state catalytic:    K1_T, K2_T, K3_T, K4_T, k13f_T, k13r_T
        #   Reg site 1, R state:  K_R1_reg1, K_R2_reg1
        #   Reg site 1, T state:  K_R1_T_reg1, K_R2_T_reg1
        #   Reg site 2, R state:  K_R3_reg2
        #   Reg site 2, T state:  K_R3_T_reg2
        #   Shared: L (= [E_T]/[E_R] for bare enzyme), Keq, Et
        #
        # Param naming uses kinetic-group representative-step indices:
        # K1=S1-binding (group 1, rep step 1), K4=P1-binding (rep step 4),
        # K7=S2-binding (rep step 7), K10=P2-binding (rep step 10),
        # k13f/k13r=catalysis SS (rep step 13).
        function rate_mwc_tetramer_bi_bi(params, concs)
            (; K1, K4, K7, K10, k13f, k13r,
               K1_T, K4_T, K7_T, K10_T, k13f_T, k13r_T,
               K_R1_reg1, K_R2_reg1, K_R1_T_reg1, K_R2_T_reg1,
               K_R3_reg2, K_R3_T_reg2,
               L, Et) = params
            (; S1, S2, P1, P2, R1, R2, R3) = concs

            # R-state: catalytic site factors by site-A ⊗ site-B independence
            Q_A_R   = 1.0 + S1 / K1 + P1 / K4
            Q_B_R   = 1.0 + S2 / K7 + P2 / K10
            Q_cat_R = Q_A_R * Q_B_R
            N_cat_R = k13f * S1 * S2 / (K1 * K7) - k13r * P1 * P2 / (K4 * K10)

            # R-state: regulatory site partition functions (star topology → direct sum)
            Q_reg1_R = 1.0 + R1 / K_R1_reg1 + R2 / K_R2_reg1
            Q_reg2_R = 1.0 + R3 / K_R3_reg2

            # T-state: same structure with _T parameters
            Q_A_T   = 1.0 + S1 / K1_T + P1 / K4_T
            Q_B_T   = 1.0 + S2 / K7_T + P2 / K10_T
            Q_cat_T = Q_A_T * Q_B_T
            N_cat_T = k13f_T * S1 * S2 / (K1_T * K7_T) - k13r_T * P1 * P2 / (K4_T * K10_T)

            Q_reg1_T = 1.0 + R1 / K_R1_T_reg1 + R2 / K_R2_T_reg1
            Q_reg2_T = 1.0 + R3 / K_R3_T_reg2

            # MWC assembly (n_cat=4, n_reg1=4, n_reg2=4)
            Q_R = Q_cat_R^4 * Q_reg1_R^4 * Q_reg2_R^4
            Q_T = Q_cat_T^4 * Q_reg1_T^4 * Q_reg2_T^4
            Z   = Q_R + L * Q_T

            num = N_cat_R * Q_cat_R^3 * Q_reg1_R^4 * Q_reg2_R^4 +
                  L * N_cat_T * Q_cat_T^3 * Q_reg1_T^4 * Q_reg2_T^4

            return Et * num / Z
        end

        push!(specs, MechanismTestSpec(
            name="MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites",
            mechanism=m,
            source_steps=src,
            source_reg_sites=src_reg,
            metabolite_names=[:S1, :S2, :P1, :P2, :R1, :R2, :R3],
            expected_n_states=9,           # catalytic subunit states
            expected_n_steps=13,           # catalytic subunit steps
            expected_n_metabolites=7,
            expected_n_haldane_constraints=2,          # one k13r per conformation (R and T)
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,      # site-independence constraints are in param_constraints of CM
            expected_n_independent_params=17,
            run_ode_test=false,
            analytical_rate_fn=rate_mwc_tetramer_bi_bi,
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end

    # ── PFK-1 hand-verified mechanism ───────────────────────────────────────
    # Reaction: F6P + ATP ⇌ F16BP + ADP, 4 catalytic subunits, 2 conformations.
    # F6P binding is :OnlyA — T-state can't bind F6P, so the T-state cycle is
    # broken in both directions and N_cat_T = 0. ATP appears as both
    # substrate and allosteric regulator (different tags per context).
    let
        m, src, src_reg = @allosteric_mechanism_src begin
            substrates: F6P, ATP
            products:   F16BP, ADP
            allosteric_regulators: Pi::EqualAI, ATP::OnlyI, ADP::OnlyA, Citrate::OnlyI, F26BP::NonequalAI

            catalytic_multiplicity: 4
            catalytic_steps: begin
                (E + F6P ⇌ E(F6P), E(ATP) + F6P ⇌ E(F6P, ATP))           :: OnlyA
                (E + ATP ⇌ E(ATP), E(F6P) + ATP ⇌ E(F6P, ATP))           :: EqualAI
                E(F6P, ATP) <--> E(F16BP, ADP)                            :: OnlyA
                (E(F16BP, ADP) ⇌ E(ADP) + F16BP, E(F16BP) ⇌ E + F16BP)   :: EqualAI
                (E(F16BP, ADP) ⇌ E(F16BP) + ADP, E(ADP) ⇌ E + ADP)       :: EqualAI
            end

            regulatory_site(multiplicity = 4): begin
                ligands: Pi, ATP
            end
        end

        function pfk_rate_analytical(params, concs)
            (; K1, K3, k5f, K6, K8,
               K_Pi_reg1, K_ATP_T_reg1,
               K_ADP_reg2, K_Citrate_T_reg3,
               K_F26BP_reg4, K_F26BP_T_reg4,
               L, Keq, Et) = params
            (; F6P, ATP, F16BP, ADP, Pi, Citrate, F26BP) = concs
            k5r = k5f * K6 * K8 / (Keq * K1 * K3)

            Q_cat_R = 1 + F6P/K1 + ATP/K3 + F6P*ATP/(K1*K3) +
                      F16BP/K6 + ADP/K8 + F16BP*ADP/(K6*K8)
            Q_cat_T = 1 + ATP/K3 + F16BP/K6 + ADP/K8 + F16BP*ADP/(K6*K8)

            # N_cat_T = 0: T-state cycle is broken (F6P binding :OnlyA), so
            # the Cha-nominal reverse term -k5r*F16BP*ADP/(K6*K8) is
            # non-physical at steady state.
            N_cat_R = k5f * F6P * ATP / (K1 * K3) - k5r * F16BP * ADP / (K6 * K8)
            N_cat_T = 0.0

            Q_reg1_R = 1 + Pi / K_Pi_reg1
            Q_reg1_T = 1 + Pi / K_Pi_reg1 + ATP / K_ATP_T_reg1
            Q_reg2_R = 1 + ADP / K_ADP_reg2
            Q_reg2_T = 1
            Q_reg3_R = 1
            Q_reg3_T = 1 + Citrate / K_Citrate_T_reg3
            Q_reg4_R = 1 + F26BP / K_F26BP_reg4
            Q_reg4_T = 1 + F26BP / K_F26BP_T_reg4

            Q_R = Q_cat_R^4 * Q_reg1_R^4 * Q_reg2_R^4 * Q_reg3_R^4 * Q_reg4_R^4
            Q_T = Q_cat_T^4 * Q_reg1_T^4 * Q_reg2_T^4 * Q_reg3_T^4 * Q_reg4_T^4
            num = N_cat_R * Q_cat_R^3 * Q_reg1_R^4 * Q_reg2_R^4 * Q_reg3_R^4 * Q_reg4_R^4 +
                  L * N_cat_T * Q_cat_T^3 * Q_reg1_T^4 * Q_reg2_T^4 * Q_reg3_T^4 * Q_reg4_T^4
            return Et * num / (Q_R + L * Q_T)
        end

        push!(specs, MechanismTestSpec(
            name="PFK-1",
            mechanism=m,
            source_steps=src,
            source_reg_sites=src_reg,
            metabolite_names=[:F6P, :ATP, :F16BP, :ADP, :Pi, :Citrate, :F26BP],
            expected_n_states=7,
            expected_n_steps=9,
            expected_n_metabolites=7,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=1,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=12,
            run_ode_test=false,
            analytical_rate_fn=pfk_rate_analytical,
            # F6P binding (group 3) is :OnlyA → the Glucose·ATP saturating
            # pattern is unreachable in T-state, so kcat = k5f
            # for every regulator corner. Regression test for the
            # `t_pattern_dead` branch in `_kcat_forward`.
            analytical_kcat_fn = p -> p.k5f,
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end

    # ── Hexokinase hand-verified mechanism ──────────────────────────────────
    #
    # G6P has THREE distinct binding sites on HK:
    #
    #   1. Catalytic site (glucose pocket): G6P competes with Glucose, but
    #      binds together with ATP at the nucleotide pocket. This is the
    #      product-release equilibrium — group 4 (K7).
    #
    #   2. Inhibitory site (= nucleotide pocket): G6P competes with ATP/ADP
    #      at the same pocket, but can co-bind with Glucose or with G6P at
    #      the catalytic site — group 6 (K12). G6P here EXCLUDES ATP/ADP at
    #      the same site, so there are no E_ATP_G6Pi or E_ADP_G6Pi forms.
    #
    #   3. Allosteric site: G6P competes with Pi at a separate regulatory
    #      site — regulatory_site(multiplicity = 2) with G6P::OnlyI and Pi::EqualAI.
    #
    # ATP binding (group 2) is :OnlyA — T-state can't bind ATP, so the
    # catalytic cycle is broken and N_cat_T = 0.
    let
        m, src, src_reg = @allosteric_mechanism_src begin
            substrates: Glucose, ATP
            products:   G6P, ADP
            allosteric_regulators: G6P::OnlyI, Pi::EqualAI
            catalytic_inhibitors:  G6P

            catalytic_multiplicity: 2
            catalytic_steps: begin
                # Group 1 (Glucose binding at catalytic site, EqualAI)
                (E + Glucose ⇌ E(Glucose),
                 E(ATP) + Glucose ⇌ E(Glucose, ATP),
                 E(G6P::Inh) + Glucose ⇌ E(Glucose, G6P::Inh))    :: EqualAI
                # Group 2 (ATP binding at nucleotide pocket, OnlyA —
                # T-state can't bind ATP)
                (E + ATP ⇌ E(ATP),
                 E(Glucose) + ATP ⇌ E(Glucose, ATP))              :: OnlyA
                # Group 3 (catalysis SS, OnlyA)
                E(Glucose, ATP) <--> E(G6P, ADP)                  :: OnlyA
                # Group 4 (G6P binding/release at catalytic site, EqualAI)
                (E(G6P, ADP) ⇌ E(ADP) + G6P,
                 E(G6P) ⇌ E + G6P,
                 E(G6P, G6P::Inh) ⇌ E(G6P::Inh) + G6P)            :: EqualAI
                # Group 5 (ADP release, EqualAI)
                (E(G6P, ADP) ⇌ E(G6P) + ADP,
                 E(ADP) ⇌ E + ADP)                                :: EqualAI
                # Group 6 (G6P binding at INHIBITORY site, EqualAI) —
                # G6P at site 2 competes with ATP/ADP, can co-bind with
                # Glucose at site 1 or G6P at site 1.
                (E + G6P::Inh ⇌ E(G6P::Inh),
                 E(Glucose) + G6P::Inh ⇌ E(Glucose, G6P::Inh),
                 E(G6P) + G6P::Inh ⇌ E(G6P, G6P::Inh))            :: EqualAI
            end

            regulatory_site(multiplicity = 2): begin
                ligands: G6P, Pi
            end
        end

        # Param naming follows kinetic-group representative-step indices:
        #   K1  (Glucose binding, group 1 rep step 1)
        #   K4  (ATP binding, group 2 rep step 4)
        #   k6f (catalysis, group 3 rep step 6)
        #   K7  (G6P at catalytic site, group 4 rep step 7)
        #   K10 (ADP release, group 5 rep step 10)
        #   K12 (G6P at inhibitory site, group 6 rep step 12) — single K
        #        for all three E_G6Pi-form bindings.
        function hk_rate_analytical(params, concs)
            (; K1, K4, k6f, K7, K10, K12,
               K_Pi_reg1, K_G6P_T_reg1,
               L, Keq, Et) = params
            (; Glucose, ATP, G6P, ADP, Pi) = concs
            k6r = k6f * K7 * K10 / (Keq * K1 * K4)

            # R-state catalytic partition function (10 enzyme forms).
            Q_cat_R = 1 +
                      Glucose / K1 +
                      ATP / K4 +
                      Glucose * ATP / (K1 * K4) +
                      G6P * ADP / (K7 * K10) +
                      G6P / K7 +
                      ADP / K10 +
                      G6P / K12 +
                      Glucose * G6P / (K1 * K12) +
                      G6P^2 / (K7 * K12)
            # T-state: ATP group :OnlyA → zero ATP terms.
            Q_cat_T = 1 +
                      Glucose / K1 +
                      G6P * ADP / (K7 * K10) +
                      G6P / K7 +
                      ADP / K10 +
                      G6P / K12 +
                      Glucose * G6P / (K1 * K12) +
                      G6P^2 / (K7 * K12)

            # N_cat_T = 0: T-state cycle is broken (ATP binding :OnlyA).
            N_cat_R = k6f * Glucose * ATP / (K1 * K4) -
                      k6r * G6P * ADP / (K7 * K10)
            N_cat_T = 0.0

            Q_reg1_R = 1 + Pi / K_Pi_reg1
            Q_reg1_T = 1 + Pi / K_Pi_reg1 + G6P / K_G6P_T_reg1

            Q_R = Q_cat_R^2 * Q_reg1_R^2
            Q_T = Q_cat_T^2 * Q_reg1_T^2
            num = N_cat_R * Q_cat_R * Q_reg1_R^2 +
                  L * N_cat_T * Q_cat_T * Q_reg1_T^2
            return Et * num / (Q_R + L * Q_T)
        end

        push!(specs, MechanismTestSpec(
            name="HK",
            mechanism=m,
            source_steps=src,
            source_reg_sites=src_reg,
            metabolite_names=[:Glucose, :ATP, :G6P, :ADP, :Pi],
            expected_n_states=10,    # 7 catalytic-cycle + 3 G6Pi dead-end
            expected_n_steps=14,     # 3+2+1+3+2+3
            expected_n_metabolites=5,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=1,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=9,
            run_ode_test=false,
            analytical_rate_fn=hk_rate_analytical,
            # ATP binding (group 2) is :OnlyA → the Glucose·ATP saturating
            # pattern is unreachable in T-state, so kcat = k6f
            # for every regulator corner. Regression test for the
            # `t_pattern_dead` branch in `_kcat_forward`.
            analytical_kcat_fn = p -> p.k6f,
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end

    # ── Pyruvate kinase (PK) hand-verified mechanism ──────────────────────────
    # Reaction: PEP + ADP ⇌ Pyruvate + ATP, 4 catalytic subunits.
    # PEP binding is :NonequalAI (independent K_R and K_T) so the T-state
    # cycle is alive. Catalysis (groups 2-5) are :EqualAI; k5r and k5r_T both
    # derive from the shared k5f via per-state Haldanes (R-state uses K1, T-state
    # uses K1_T). Reg sites have MISMATCHED multiplicities:
    #   ATP::OnlyI at mult 2
    #   F16BP::OnlyA at mult 4 (matches catalytic mult)
    # This exercises the symmetric all-reg-sites contribution to both numerator
    # and denominator.
    # Independent parameters (9): K1, K1_T, K3, k5f, K6, K8,
    # K_ATP_T_reg1, K_F16BP_reg2, L
    let
        m, src, src_reg = @allosteric_mechanism_src begin
            substrates: PEP, ADP
            products:   Pyruvate, ATP
            allosteric_regulators: ATP::OnlyI, F16BP::OnlyA

            catalytic_multiplicity: 4
            catalytic_steps: begin
                (E + PEP ⇌ E(PEP),
                 E(ADP) + PEP ⇌ E(PEP, ADP))                          :: OnlyA
                (E + ADP ⇌ E(ADP),
                 E(PEP) + ADP ⇌ E(PEP, ADP))                          :: EqualAI
                E(PEP, ADP) <--> E(Pyruvate, ATP)                     :: OnlyA
                (E(Pyruvate, ATP) ⇌ E(ATP) + Pyruvate,
                 E(Pyruvate) ⇌ E + Pyruvate)                          :: EqualAI
                (E(Pyruvate, ATP) ⇌ E(Pyruvate) + ATP,
                 E(ATP) ⇌ E + ATP)                                    :: EqualAI
            end

            regulatory_site(multiplicity = 2): begin
                ligands: ATP
            end
            regulatory_site(multiplicity = 4): begin
                ligands: F16BP
            end
        end

        # Param mapping:
        #   K1  : PEP binding (group 1, :OnlyA — R-state only)
        #   K3  : ADP binding (group 2, :EqualAI)
        #   k5f : catalysis SS forward rate (group 3, :OnlyA)
        #   K6  : Pyruvate release (group 4, :EqualAI)
        #   K8  : ATP release (group 5, :EqualAI)
        #
        # PEP binding AND catalysis are :OnlyA (the exclusive-binding K-system):
        # the T-state cannot bind PEP, so it cannot catalyze — the T-catalytic
        # cycle is dead (N_T = 0) and the T-state is a clean binding partition with
        # no PEP term. k5r derives via the R-state Haldane k5r = k5f·K6·K8/(Keq·K1·K3).
        # At saturation the R-state dominates, so forward kcat = k5f.
        function pk_rate_analytical(params, concs)
            (; K1, K3, k5f, K6, K8,
               K_ATP_T_reg1, K_F16BP_reg2,
               L, Keq, Et) = params
            (; PEP, ADP, Pyruvate, ATP, F16BP) = concs

            k5r = k5f * K6 * K8 / (Keq * K1 * K3)

            Q_cat_R = 1 + PEP/K1 + ADP/K3 + PEP*ADP/(K1 * K3) +
                      Pyruvate/K6 + ATP/K8 + Pyruvate*ATP/(K6 * K8)
            Q_cat_T = 1 + ADP/K3 + Pyruvate/K6 + ATP/K8 +
                      Pyruvate*ATP/(K6 * K8)                    # no PEP (:OnlyA)

            N_R = k5f * PEP * ADP / (K1 * K3) - k5r * Pyruvate * ATP / (K6 * K8)
            # N_T = 0: catalysis :OnlyA — the dead T-state contributes no flux.

            Q_reg1_R = 1                                     # ATP::OnlyI, no R term
            Q_reg1_T = 1 + ATP / K_ATP_T_reg1
            Q_reg2_R = 1 + F16BP / K_F16BP_reg2               # F16BP::OnlyA, no T term
            Q_reg2_T = 1

            num = N_R * Q_cat_R^3 * Q_reg1_R^2 * Q_reg2_R^4
            den = Q_cat_R^4 * Q_reg1_R^2 * Q_reg2_R^4 +
                  L * Q_cat_T^4 * Q_reg1_T^2 * Q_reg2_T^4
            return Et * num / den
        end

        push!(specs, MechanismTestSpec(
            name="PK",
            mechanism=m,
            source_steps=src,
            source_reg_sites=src_reg,
            metabolite_names=[:PEP, :ADP, :Pyruvate, :ATP, :F16BP],
            expected_n_states=7,           # E, E_PEP, E_ADP, E_PEP_ADP, E_Pyr_ATP, E_Pyr, E_ATP
            expected_n_steps=9,
            expected_n_metabolites=5,
            expected_n_haldane_constraints=1,
            # PEP binding and catalysis are :OnlyA (the T-state is pruned), so
            # there is no mirror/collapse — the T-catalytic cycle is simply dead
            # (N_T = 0). One R-state Haldane derives k5r; no Wegscheider tie.
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=8,
            run_ode_test=false,
            analytical_rate_fn=pk_rate_analytical,
            analytical_kcat_fn = p -> p.k5f,
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end

    # ── m_all: NonequalAI coverage on substrate + catalysis + product ───────
    # Two-substrate two-product reaction with explicit :NonequalAI on
    # S1 binding (group 1), catalysis (group 3), and P1 release (group 4),
    # and a 2-ligand mixed-state reg site (R1::NonequalAI + R2::EqualAI).
    # Catalysis :NonequalAI combined with :NonequalAI substrate yields
    # independent R-state and T-state Haldanes (one per state).
    # Independent parameters (12): K1, K1_T, K3, k5f, k5f_T, K6, K6_T, K8,
    # K_R1_reg1, K_R1_T_reg1, K_R2_reg1, L
    let
        m, src, src_reg = @allosteric_mechanism_src begin
            substrates: S1, S2
            products:   P1, P2
            allosteric_regulators: R1::NonequalAI, R2::EqualAI

            catalytic_multiplicity: 2
            catalytic_steps: begin
                (E + S1 ⇌ E(S1),
                 E(S2) + S1 ⇌ E(S1, S2))    :: NonequalAI
                (E + S2 ⇌ E(S2),
                 E(S1) + S2 ⇌ E(S1, S2))    :: EqualAI
                E(S1, S2) <--> E(P1, P2)    :: NonequalAI
                (E(P1, P2) ⇌ E(P2) + P1,
                 E(P1) ⇌ E + P1)            :: NonequalAI
                (E(P1, P2) ⇌ E(P1) + P2,
                 E(P2) ⇌ E + P2)            :: EqualAI
            end

            regulatory_site(multiplicity = 2): begin
                ligands: R1, R2
            end
        end

        # Param mapping (kinetic-group representative-step convention):
        #   K1, K1_T : S1 binding (group 1, NonequalAI)
        #   K3       : S2 binding (group 2, EqualAI)
        #   k5f, k5f_T : catalysis SS (group 3, NonequalAI)
        #   K6, K6_T : P1 release (group 4, NonequalAI)
        #   K8       : P2 release (group 5, EqualAI)
        function m_all_rate_analytical(params, concs)
            (; K1, K1_T, K3, k5f, k5f_T, K6, K6_T, K8,
               K_R1_reg1, K_R1_T_reg1, K_R2_reg1,
               L, Keq, Et) = params
            (; S1, S2, P1, P2, R1, R2) = concs

            k5r   = k5f   * K6   * K8 / (Keq * K1   * K3)
            k5r_T = k5f_T * K6_T * K8 / (Keq * K1_T * K3)

            Q_cat_R = 1 + S1/K1   + S2/K3 + S1*S2/(K1   * K3) +
                      P1/K6   + P2/K8 + P1*P2/(K6   * K8)
            Q_cat_T = 1 + S1/K1_T + S2/K3 + S1*S2/(K1_T * K3) +
                      P1/K6_T + P2/K8 + P1*P2/(K6_T * K8)

            N_R = k5f   * S1 * S2 / (K1   * K3) - k5r   * P1 * P2 / (K6   * K8)
            N_T = k5f_T * S1 * S2 / (K1_T * K3) - k5r_T * P1 * P2 / (K6_T * K8)

            Q_reg1_R = 1 + R1/K_R1_reg1   + R2/K_R2_reg1
            Q_reg1_T = 1 + R1/K_R1_T_reg1 + R2/K_R2_reg1

            num = N_R * Q_cat_R   * Q_reg1_R^2 + L * N_T * Q_cat_T   * Q_reg1_T^2
            den = Q_cat_R^2 * Q_reg1_R^2       + L *      Q_cat_T^2 * Q_reg1_T^2

            return Et * num / den
        end

        push!(specs, MechanismTestSpec(
            name="m_all",
            mechanism=m,
            source_steps=src,
            source_reg_sites=src_reg,
            metabolite_names=[:S1, :S2, :P1, :P2, :R1, :R2],
            expected_n_states=7,
            expected_n_steps=9,
            expected_n_metabolites=6,
            expected_n_haldane_constraints=2,
            # structural naming: :EqualAI catalytic groups share one symbol (no rename); only :EqualAI reg ligands emit a mirror
            expected_n_mirror_constraints=1,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=12,
            run_ode_test=false,
            analytical_rate_fn=m_all_rate_analytical,
            analytical_kcat_fn=nothing,      # cat is :NonequalAI → kcat L-dependent
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end

    # ── m_OnlyA_prod: catalysis :OnlyA → T-state truly inactive ─────────────
    # Exercises `_i_state_num_zero` detection when the catalytic conversion is
    # :OnlyA: the T-state cannot run the reaction (forward or reverse), so its
    # cycle is dead — N_T is 0 and the L*num_T branch is dropped. Analytical
    # kcat = 2·k2f/(1+L) — L-dependent because the saturating R-state
    # pattern (S only) IS reachable in T-state (Q_cat_T at sat S = S/K1,
    # same as Q_cat_R), so B_T ≠ 0 and the 1/(1+L) factor appears.
    let
        m, src, src_reg = @allosteric_mechanism_src begin
            substrates: S
            products:   P

            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)    :: EqualAI      # group 1, K1
                E(S) <--> E(P)  :: OnlyA        # group 2, k2f catalysis :OnlyA
                E + P ⇌ E(P)    :: OnlyA        # group 3, K3 (P binding)
            end
        end

        # Param mapping:
        #   K1   : S binding (group 1, EqualAI)
        #   k2f  : catalysis SS (group 2, EqualAI, k2r derived via Haldane)
        #   K3   : P release (group 3, OnlyA)
        function m_OnlyA_prod_rate_analytical(params, concs)
            (; K1, k2f, K3, L, Keq, Et) = params
            (; S, P) = concs

            k2r = k2f * K3 / (Keq * K1)

            Q_cat_R = 1 + S/K1 + P/K3
            Q_cat_T = 1 + S/K1                    # E(P) unreachable in T-state (catalysis + P-binding :OnlyA)

            N_R = k2f * S/K1 - k2r * P/K3
            # N_T = 0 forced (t_state_dead via catalysis group 2 :OnlyA)

            num = N_R * Q_cat_R                   # L*N_T*Q_cat_T term elided
            den = Q_cat_R^2 + L * Q_cat_T^2

            return Et * num / den
        end

        push!(specs, MechanismTestSpec(
            name="m_OnlyA_prod",
            mechanism=m,
            source_steps=src,
            source_reg_sites=src_reg,
            metabolite_names=[:S, :P],
            expected_n_states=3,                  # E, E_S, E_P
            expected_n_steps=3,
            expected_n_metabolites=2,
            expected_n_haldane_constraints=1,
            expected_n_mirror_constraints=0,
            expected_n_wegscheider_constraints=0,
            expected_n_independent_params=4,
            run_ode_test=false,
            analytical_rate_fn=m_OnlyA_prod_rate_analytical,
            # kcat at saturating S, zero P:
            #   A_R = k2f/K1², B_R = 1/K1², B_T = 1/K1² (T-state pattern same as R)
            #   kcat = A_R / (B_R + L · B_T) = k2f / (1 + L)
            analytical_kcat_fn = p -> p.k2f / (1 + p.L),
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end

    # The three LDH i-state mechanisms that exposed the Bug-2 fitted_params
    # leak (a Haldane-dependent reverse rate that landed in the independent
    # set). Written as @allosteric_mechanism for readability; RE-containing
    # with a products=0 boundary, so the ODE cross-check is skipped. Their
    # golden PARAMS_REDUCED confirm the dependent reverse rate is excluded.
    push!(specs, MechanismTestSpec(
        name="LDH i-state NonequalAI 6-group",
        mechanism=(@allosteric_mechanism begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            catalytic_multiplicity: 4
            catalytic_steps: begin
                (E + Lactate <--> E(Lactate),
                 E(NADH) + Lactate <--> E(Lactate, NADH)) :: NonequalAI
                (E + NAD ⇌ E(NAD),
                 E(Lactate) + NAD ⇌ E(Lactate, NAD),
                 E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)) :: EqualAI
                (E + NADH <--> E(NADH),
                 E(Lactate) + NADH <--> E(Lactate, NADH),
                 E(Pyruvate) + NADH <--> E(NADH, Pyruvate)) :: OnlyA
                (E + Pyruvate ⇌ E(Pyruvate),
                 E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
                 E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)) :: EqualAI
                E(NAD) + Lactate ⇌ E(Lactate, NAD) :: EqualAI
                E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
            end
        end),
        metabolite_names=[:NADH, :Pyruvate, :Lactate, :NAD],
        expected_n_states=9, expected_n_steps=13, expected_n_metabolites=4,
        expected_n_haldane_constraints=1, expected_n_mirror_constraints=0,
        expected_n_wegscheider_constraints=2, expected_n_independent_params=9,
        run_ode_test=false))

    push!(specs, MechanismTestSpec(
        name="LDH i-state NonequalAI 5-group",
        mechanism=(@allosteric_mechanism begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            catalytic_multiplicity: 4
            catalytic_steps: begin
                (E + Lactate ⇌ E(Lactate),
                 E(NAD) + Lactate ⇌ E(Lactate, NAD),
                 E(NADH) + Lactate ⇌ E(Lactate, NADH)) :: NonequalAI
                (E + NAD ⇌ E(NAD),
                 E(Lactate) + NAD ⇌ E(Lactate, NAD)) :: EqualAI
                (E + NADH <--> E(NADH),
                 E(Lactate) + NADH <--> E(Lactate, NADH)) :: OnlyA
                E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate) :: EqualAI
                E(NADH, Pyruvate) <--> E(Lactate, NAD) :: OnlyA
            end
        end),
        metabolite_names=[:NADH, :Pyruvate, :Lactate, :NAD],
        expected_n_states=7, expected_n_steps=9, expected_n_metabolites=4,
        expected_n_haldane_constraints=1, expected_n_mirror_constraints=0,
        expected_n_wegscheider_constraints=0, expected_n_independent_params=8,
        run_ode_test=false))

    push!(specs, MechanismTestSpec(
        name="LDH i-state EqualAI-NonequalAI 6-group",
        mechanism=(@allosteric_mechanism begin
            substrates: NADH, Pyruvate
            products: Lactate, NAD
            catalytic_multiplicity: 4
            catalytic_steps: begin
                E + NAD ⇌ E(NAD) :: EqualAI
                (E + NADH ⇌ E(NADH),
                 E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)) :: EqualAI
                (E + Pyruvate ⇌ E(Pyruvate),
                 E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
                 E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)) :: EqualAI
                (E(NAD) + Lactate ⇌ E(Lactate, NAD),
                 E(NADH) + Lactate ⇌ E(Lactate, NADH)) :: EqualAI
                E(NADH, Pyruvate) <--> E(Lactate, NAD) :: EqualAI
                E(Pyruvate) + NAD <--> E(NAD, Pyruvate) :: NonequalAI
            end
        end),
        metabolite_names=[:NADH, :Pyruvate, :Lactate, :NAD],
        expected_n_states=8, expected_n_steps=10, expected_n_metabolites=4,
        expected_n_haldane_constraints=1, expected_n_mirror_constraints=0,
        expected_n_wegscheider_constraints=2, expected_n_independent_params=8,
        run_ode_test=false))

    # NOTE: a multi-:OnlyA derivation/perf spec was intentionally NOT added here.
    # The representative multi-:OnlyA mechanism triggers the pre-existing allosteric
    # MWC L-term leak (its inactive graph fragments), so its derivation is
    # known-incorrect until that bug is fixed — see
    # docs/superpowers/specs/2026-07-13-allosteric-mwc-derivation-known-issues.md.
    # The enumeration move that makes multi-:OnlyA reachable is validated by its own
    # tests in test_mechanism_enumeration.jl; the n=1 mass-action ground truth for the
    # multi-:OnlyA derivation lives (as an @test_broken gate) in allosteric_ground_truth.jl.

    return specs
end

const MECHANISM_TEST_SPECS = build_mechanism_test_specs()

"""Look up a mechanism test spec by name."""
function _spec_by_name(name)
    for s in MECHANISM_TEST_SPECS
        s.name == name && return s
    end
    error("MechanismTestSpec not found: $name")
end

const pfk_mechanism = _spec_by_name("PFK-1").mechanism
const pfk_rate_analytical = _spec_by_name("PFK-1").analytical_rate_fn
const hk_mechanism = _spec_by_name("HK").mechanism
const hk_rate_analytical = _spec_by_name("HK").analytical_rate_fn
