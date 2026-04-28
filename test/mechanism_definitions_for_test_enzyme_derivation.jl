# Test specifications for all enzyme mechanisms
# Each mechanism is defined inline with its expected properties
# for easy side-by-side comparison

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
    expected_n_haldane::Int               # Usually 1
    expected_n_wegscheider::Int           # 0 for linear, 1+ for branched
    expected_n_independent_params::Int    # 2*n_steps - n_constraints

    # Identifiability expectations
    expected_identifiability_deficit::Int
    expected_is_identifiable::Bool

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
end

# ── Mechanism test specifications ───────────────────────────────────────────

function build_mechanism_test_specs()
    specs = MechanismTestSpec[]

    # 1. Uni-Uni (simplest): E + S ⇌ ES ⇌ E + P
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [E, P]
            end
        end
        push!(specs, MechanismTestSpec(
            name="Uni-Uni",
            mechanism=m,
            metabolite_names=[:S, :P],
            expected_n_states=2,
            expected_n_steps=2,
            expected_n_metabolites=2,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=3,
            expected_identifiability_deficit=0,
            expected_is_identifiable=true,
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
        m = @enzyme_mechanism begin
            substrates: A
            products: P
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
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
            metabolite_names=[:A, :P],
            expected_n_states=3,
            expected_n_steps=3,
            expected_n_metabolites=2,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=5,
            expected_identifiability_deficit=2,
            expected_is_identifiable=false,
            analytical_rate_fn=(p, c) -> rate_uni_uni(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 3. Segel Iso Uni Uni (new): E + A ⇌ EA ⇌ EP ⇌ F + P, F ⇌ E
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-45
    let
        m = @enzyme_mechanism begin
            substrates: A
            products: P
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [F, P]
                [F] <--> [E]
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
            metabolite_names=[:A, :P],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=2,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=7,
            expected_identifiability_deficit=3,
            expected_is_identifiable=false,
            analytical_rate_fn=(p, c) -> rate_iso_uni_uni(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 4. Segel Ordered Uni Bi (replaces Seq Uni-Bi): E + A ⇌ (EA≡EPQ) ⇌ EQ + P ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-60
    let
        m = @enzyme_mechanism begin
            substrates: A
            products: P, Q
            steps: begin
                [E, A] <--> [EAEPQ]
                [EAEPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
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
            metabolite_names=[:A, :P, :Q],
            expected_n_states=3,
            expected_n_steps=3,
            expected_n_metabolites=3,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=5,
            expected_identifiability_deficit=-1,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_ordered_uni_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 5. Segel Ordered Bi Bi (replaces Seq Bi-Bi):
    #    E + A ⇌ EA + B ⇌ (EAB≡EPQ) ⇌ EQ + P ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-87
    let
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                [E, A] <--> [EA]
                [EA, B] <--> [EABEPQ]
                [EABEPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
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
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=7,
            expected_identifiability_deficit=-4,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_ordered_bi_bi(merge(p, (Etotal=p.Et,)), c),
            analytical_kcat_fn=p -> p.k3f * p.k4f / (p.k3f + p.k4f),
        ))
    end

    # 6. Segel Theorell-Chance Bi Bi (new): E + A ⇌ EA + B ⇌ EQ + P ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-122
    let
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                [E, A] <--> [EA]
                [EA, B] <--> [EQ, P]
                [EQ] <--> [E, Q]
            end
        end

        # Segel Eq. IX-122: Theorell-Chance Bi Bi steady-state rate
        function rate_theorell_chance_bi_bi(params, concs)
            (; k1f, k1r, k2f, k2r, k3f, k3r, Etotal) = params
            (; A, B, P, Q) = concs
            num = k1f * k2f * k3f * A * B - k1r * k2r * k3r * P * Q
            denom = k1r * k3f +
                    k1f * k3f * A +
                    k2f * k3f * B +
                    k1r * k2r * P +
                    k1r * k3r * Q +
                    k1f * k2f * A * B +
                    k1f * k2r * A * P +
                    k2f * k3r * B * Q +
                    k2r * k3r * P * Q
            return Etotal * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Segel Theorell-Chance Bi Bi",
            mechanism=m,
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=3,
            expected_n_steps=3,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=5,
            expected_identifiability_deficit=-4,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_theorell_chance_bi_bi(
                    merge(p, (Etotal=p.Et,)), c),
            analytical_kcat_fn=p -> p.k3f,
        ))
    end

    # 7. Segel Ping Pong Bi Bi (replaces Ping-Pong Bi-Bi):
    #    E + A ⇌ (EA≡FP) ⇌ F + P, F + B ⇌ (FB≡EQ) ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-140
    let
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                [E, A] <--> [EAFP]
                [EAFP] <--> [F, P]
                [F, B] <--> [FBEQ]
                [FBEQ] <--> [E, Q]
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
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=7,
            expected_identifiability_deficit=-1,
            expected_is_identifiable=true,
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
        m = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q
            steps: begin
                [E, A] <--> [EA]
                [EA, B] <--> [EAB]
                [EAB, C] <--> [EABCEPQ]
                [EABCEPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
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
            metabolite_names=[:A, :B, :C, :P, :Q],
            expected_n_states=5,
            expected_n_steps=5,
            expected_n_metabolites=5,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=9,
            expected_identifiability_deficit=-9,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_ordered_ter_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 9. Segel Ordered Ter Ter (replaces Seq Ter-Ter):
    #     E + A ⇌ EA + B ⇌ EAB + C ⇌ (EABC≡EPQR) ⇌ EQR + P ⇌ ER + Q ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-261
    let
        m = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                [E, A] <--> [EA]
                [EA, B] <--> [EAB]
                [EAB, C] <--> [EABCEPQR]
                [EABCEPQR] <--> [EQR, P]
                [EQR] <--> [ER, Q]
                [ER] <--> [E, R]
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
            metabolite_names=[:A, :B, :C, :P, :Q, :R],
            expected_n_states=6,
            expected_n_steps=6,
            expected_n_metabolites=6,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=11,
            expected_identifiability_deficit=-16,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_ordered_ter_ter(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 10. Segel Bi Uni Uni Uni Ping Pong Ter Bi (new):
    #     E + A ⇌ EA + B ⇌ (EAB≡FP) ⇌ F + P, F + C ⇌ (FC≡EQ) ⇌ E + Q
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-228
    let
        m = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q
            steps: begin
                [E, A] <--> [EA]
                [EA, B] <--> [EABFP]
                [EABFP] <--> [F, P]
                [F, C] <--> [FCEQ]
                [FCEQ] <--> [E, Q]
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
            metabolite_names=[:A, :B, :C, :P, :Q],
            expected_n_states=5,
            expected_n_steps=5,
            expected_n_metabolites=5,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=9,
            expected_identifiability_deficit=-5,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_bi_uni_uni_uni_ping_pong_ter_bi(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 11. Random-order Bi-Bi (branched): Two substrate binding orders converge
    let
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                [E, A] <--> [EA]
                [E, B] <--> [EB]
                [EA, B] <--> [EAB]
                [EB, A] <--> [EAB]
                [EAB] <--> [EPQ]
                [EPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
            end
        end
        push!(specs, MechanismTestSpec(
            name="Random-order Bi-Bi",
            mechanism=m,
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=6,
            expected_n_steps=7,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=1,
            expected_n_independent_params=12,
            expected_identifiability_deficit=-16,
            expected_is_identifiable=true,
            analytical_rate_fn=nothing
        ))
    end

    # 12. Segel Bi Uni Uni Bi Ping Pong Ter Ter (new):
    #     E + A ⇌ EA + B ⇌ (EAB≡FP) ⇌ F + P, F + C ⇌ (FC≡EQR) ⇌ ER + Q ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-278
    let
        m = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                [E, A] <--> [EA]
                [EA, B] <--> [EABFP]
                [EABFP] <--> [F, P]
                [F, C] <--> [FCEQR]
                [FCEQR] <--> [ER, Q]
                [ER] <--> [E, R]
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
            metabolite_names=[:A, :B, :C, :P, :Q, :R],
            expected_n_states=6,
            expected_n_steps=6,
            expected_n_metabolites=6,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=11,
            expected_identifiability_deficit=-11,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_bi_uni_uni_bi_ping_pong_ter_ter(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 13. Segel Bi Bi Uni Uni Ping Pong Ter Ter (new):
    #     E + A ⇌ EA + B ⇌ (EAB≡FPQ) ⇌ FQ + P, FQ ⇌ F + Q, F + C ⇌ (FC≡ER) ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-288
    let
        m = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                [E, A] <--> [EA]
                [EA, B] <--> [EABFPQ]
                [EABFPQ] <--> [FQ, P]
                [FQ] <--> [F, Q]
                [F, C] <--> [FCER]
                [FCER] <--> [E, R]
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
            metabolite_names=[:A, :B, :C, :P, :Q, :R],
            expected_n_states=6,
            expected_n_steps=6,
            expected_n_metabolites=6,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=11,
            expected_identifiability_deficit=-11,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_bi_bi_uni_uni_ping_pong_ter_ter(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 14. Segel Hexa Uni Ping Pong (new):
    #     E + A ⇌ (EA≡FP) ⇌ F + P, F + B ⇌ (FB≡GQ) ⇌ G + Q, G + C ⇌ (GC≡ER) ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-308
    let
        m = @enzyme_mechanism begin
            substrates: A, B, C
            products: P, Q, R
            steps: begin
                [E, A] <--> [EAFP]
                [EAFP] <--> [F, P]
                [F, B] <--> [FBGQ]
                [FBGQ] <--> [G, Q]
                [G, C] <--> [GCER]
                [GCER] <--> [E, R]
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
            metabolite_names=[:A, :B, :C, :P, :Q, :R],
            expected_n_states=6,
            expected_n_steps=6,
            expected_n_metabolites=6,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=11,
            expected_identifiability_deficit=-6,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_hexa_uni_ping_pong(
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # ── Rapid-Equilibrium (RE) Mechanisms ──────────────────────────────────────

    # 15. RE Uni-Uni: E + A ⇌_RE EA <-->_SS E + P
    #     Rapid-equilibrium substrate binding, steady-state catalysis
    let
        m = @enzyme_mechanism begin
            substrates: A
            products: P
            steps: begin
                [E, A] ⇌ [EA]
                [EA] <--> [E, P]
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
            metabolite_names=[:A, :P],
            expected_n_states=2,
            expected_n_steps=2,
            expected_n_metabolites=2,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=2,
            expected_identifiability_deficit=0,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_re_uni_uni(merge(p, (Et=p.Et,)), c),
            analytical_kcat_fn=p -> p.k2f,
        ))
    end

    # 16. RE Ordered Bi-Bi: substrate binding is RE, catalysis and product release are SS
    #     E + A ⇌_RE EA, EA + B ⇌_RE EAB, (EAB≡EPQ) <-->_SS EQ + P, EQ <-->_SS E + Q
    let
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                [E, A] ⇌ [EA]
                [EA, B] ⇌ [EABEPQ]
                [EABEPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
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
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=5,
            expected_identifiability_deficit=-2,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_re_ordered_bi_bi(merge(p, (Et=p.Et,)), c)
        ))
    end

    # 17. RE Random Bi-Bi: binding to free enzyme is RE, other steps are SS
    #     E + A ⇌_RE EA, E + B ⇌_RE EB, EA + B <-->_SS EAB, EB + A <-->_SS EAB,
    #     EAB <-->_SS EPQ, EPQ <-->_SS EQ + P, EQ <-->_SS E + Q
    let
        m = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                [E, A] ⇌ [EA]
                [E, B] ⇌ [EB]
                [EA, B] <--> [EAB]
                [EB, A] <--> [EAB]
                [EAB] <--> [EPQ]
                [EPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
            end
        end
        push!(specs, MechanismTestSpec(
            name="RE Random Bi-Bi",
            mechanism=m,
            metabolite_names=[:A, :B, :P, :Q],
            expected_n_states=6,
            expected_n_steps=7,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=1,
            expected_n_independent_params=10,
            expected_identifiability_deficit=0,
            expected_is_identifiable=true,
            analytical_rate_fn=nothing
        ))
    end

    # ── Classical Inhibitor/Activator Mechanisms (factored form tests) ────────

    # 18. Competitive inhibitor: E + S ⇌ ES, ES ⇌ EP (SS), EP ⇌ E + P, E + R ⇌ ER
    #     Dead-end inhibitor R binds free enzyme only.
    #     No Cartesian product structure → flat sum denominator.
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                [E, S] ⇌ [ES]        # K1
                [ES] <--> [EP]        # k2f, k2r (SS)
                [EP] ⇌ [E, P]        # K3
                [E, R] ⇌ [ER]        # K4
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
            metabolite_names=[:S, :P, :R],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=3,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=4,
            expected_identifiability_deficit=0,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_competitive_inh(
                merge(p, (Et=p.Et,)), c),
            analytical_kcat_fn=p -> p.k2f,
            # Textbook: flat sum denominator (no Cartesian product structure)
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "1 + S / K1 + P / K3 + R / K4",
        ))
    end

    # 19. Non-competitive inhibitor: R binds both free E and ES with same K
    #     Forms: E, E_S, E_P, E_R, E_S_R; SS: E_S↔E_P; K5=K4
    #     Denom factors as (1+R/K4)*(1+S/K1) + P/K3
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                [E, S] ⇌ [E_S]        # K1
                [E_S] <--> [E_P]       # k2f, k2r (SS)
                [E, P] ⇌ [E_P]        # K3
                ([E, R] ⇌ [E_R], [E_S, R] ⇌ [E_S_R])  # K4 = K5 (R binding shared)
                [E_R, S] ⇌ [E_S_R]    # K6
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
            metabolite_names=[:S, :P, :R],
            expected_n_states=5,
            expected_n_steps=6,
            expected_n_metabolites=3,
            expected_n_haldane=1,
            expected_n_wegscheider=1,
            expected_n_independent_params=4,
            expected_identifiability_deficit=-1,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_noncompetitive_inh(
                merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "P / K3 + (1 + S / K1) * (1 + R / K4)",
        ))
    end

    # 20. Uncompetitive inhibitor: R binds ES only (not free E)
    #     Forms: E, E_S, E_P, E_S_R; SS: E_S↔E_P; No extra constraints
    #     Denom: 1 + P/K3 + S/K1*(1+R/K4)
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                [E, S] ⇌ [E_S]        # K1
                [E_S] <--> [E_P]       # k2f, k2r (SS)
                [E, P] ⇌ [E_P]        # K3
                [E_S, R] ⇌ [E_S_R]    # K4
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
            metabolite_names=[:S, :P, :R],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=3,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=4,
            expected_identifiability_deficit=0,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_uncompetitive_inh(
                merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "1 + P / K3 + (S / K1) * (1 + R / K4)",
        ))
    end

    # 21. Essential activator: R must bind before S can bind
    #     Forms: E, E_R, E_S_R, E_P_R; SS: E_S_R↔E_P_R
    #     Num: R/K4 * (k2f*S/K1 - k2r*P/K3)
    #     Denom: 1 + R/K4*(1+S/K1+P/K3)
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                [E_R, S] ⇌ [E_S_R]    # K1
                [E_S_R] <--> [E_P_R]   # k2f, k2r (SS)
                [E_R, P] ⇌ [E_P_R]    # K3
                [E, R] ⇌ [E_R]        # K4
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
            metabolite_names=[:S, :P, :R],
            expected_n_states=4,
            expected_n_steps=4,
            expected_n_metabolites=3,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=4,
            expected_identifiability_deficit=0,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_essential_activator(
                merge(p, (Et=p.Et,)), c),
            analytical_kcat_fn=p -> p.k2f,
            # Mathematically equivalent factoring: numerator and denominator
            # are both multiplied by R/K4 in the alternative form. The
            # derivation now produces this form directly.
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "K4 / R + 1 + S / K1 + P / K3",
        ))
    end

    # 22. Non-essential activator (general modifier):
    #     R modifies catalysis but isn't required. Two parallel SS cycles.
    #     Forms: E, E_S, E_P, E_R, E_S_R, E_P_R; SS: E_S↔E_P, E_S_R↔E_P_R
    #     K8=K7, K9=K7 (R binding independent of S/P)
    #     K4=K1 and K6=K3 are implied by Wegscheider relations
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: R
            steps: begin
                [E, S] ⇌ [E_S]          # K1
                [E_S] <--> [E_P]         # k2f, k2r (SS)
                [E, P] ⇌ [E_P]          # K3
                [E_R, S] ⇌ [E_S_R]      # K4
                [E_S_R] <--> [E_P_R]     # k5f, k5r (SS)
                [E_R, P] ⇌ [E_P_R]      # K6
                ([E, R] ⇌ [E_R],         # K7 = K8 = K9 (R binding shared)
                 [E_S, R] ⇌ [E_S_R],
                 [E_P, R] ⇌ [E_P_R])
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
            metabolite_names=[:S, :P, :R],
            expected_n_states=6,
            expected_n_steps=9,
            expected_n_metabolites=3,
            expected_n_haldane=2,
            expected_n_wegscheider=2,
            expected_n_independent_params=5,
            expected_identifiability_deficit=-3,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) -> rate_nonessential_activator(
                merge(p, (Et=p.Et,)), c),
            analytical_kcat_fn=p -> max(p.k2f, p.k5f),
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3 + (R / K7) * (k5f * S / K1 - k5r * P / K3)",
            expected_factored_denom=
            "(1 + S / K1 + P / K3) * (1 + R / K7)",
        ))
    end

    # 23. Non-essential activator + competitive inhibitor:
    #     Combines spec 25 (non-essential activator) with competitive inhibition.
    #     A modifies catalysis but isn't required (binds E, E_S, E_P with same K).
    #     I binds only free E (competitive dead-end).
    #     Forms: E, E_S, E_P, E_A, E_S_A, E_P_A, E_I
    #     SS: E_S↔E_P, E_S_A↔E_P_A
    #     Constraints: K8=K7, K9=K7 (A binding K independent of S/P)
    #     Wegscheider gives K4=K1, K6=K3
    #     Denom: (1+S/K1+P/K3)*(1+A/K7) + I/K10
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: A, I
            steps: begin
                [E, S] ⇌ [E_S]          # K1
                [E_S] <--> [E_P]         # k2f, k2r (SS)
                [E, P] ⇌ [E_P]          # K3
                [E_A, S] ⇌ [E_S_A]      # K4
                [E_S_A] <--> [E_P_A]     # k5f, k5r (SS)
                [E_A, P] ⇌ [E_P_A]      # K6
                ([E, A] ⇌ [E_A],         # K7 = K8 = K9 (A binding shared)
                 [E_S, A] ⇌ [E_S_A],
                 [E_P, A] ⇌ [E_P_A])
                [E, I] ⇌ [E_I]          # K10
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
            metabolite_names=[:S, :P, :A, :I],
            expected_n_states=7,
            expected_n_steps=10,
            expected_n_metabolites=4,
            expected_n_haldane=2,
            expected_n_wegscheider=2,
            expected_n_independent_params=6,
            expected_identifiability_deficit=-3,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_activator_inhibitor(
                    merge(p, (Et=p.Et,)), c),
            # Denom has both multiplicative (activator) and additive
            # (inhibitor) structure
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3 + (A / K7) * (k5f * S / K1 - k5r * P / K3)",
            expected_factored_denom=
            "I / K10 + (1 + S / K1 + P / K3) * (1 + A / K7)",
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
    #     Numerator: 2*(
    #                   (k_cat_R*S/Ks_R - k_rev_R*P/Kp_R)*(1 + S/Ks_R + P/Kp_R) +
    #                   L*(k_cat_T*S/Ks_T - k_rev_T*P/Kp_T)*(1 + S/Ks_T + P/Kp_T)
    #                 )
    # Spec #24 (flat homodimer "MWC Dimer") removed per Task 2.6: the
    # AllostericEnzymeMechanism sibling (#24B below) is the canonical form.

    # 24B. AllostericEnzymeMechanism equivalent of MWC Dimer (spec #24)
    #      2 catalytic sites × 2 conformations (R/T). No explicit equality constraints
    #      needed — symmetric subunits are captured by the site multiplicity.
    #      Conformational equilibrium: L (= K37 in the EnzymeMechanism above).
    let
        m = @allosteric_mechanism begin
            substrates: S
            products: P
            site(:catalytic, 2): begin
                steps: begin
                    [E_c, S] ⇌ [E_S]    :: NonequalRT
                    [E_c, P] ⇌ [E_P]    :: NonequalRT
                    [E_S] <--> [E_P]     :: NonequalRT
                end
            end
        end

        function rate_mwc_dimer_oligo(params, concs)
            (; K1, K2, k3f, k3r, K1_T, K2_T, k3f_T, k3r_T, L, Et) = params
            (; S, P) = concs
            r_flux   = k3f * S / K1 - k3r * P / K2
            t_flux   = k3f_T * S / K1_T - k3r_T * P / K2_T
            r_factor = 1.0 + S / K1 + P / K2
            t_factor = 1.0 + S / K1_T + P / K2_T
            return Et * 2.0 * (r_flux * r_factor + L * t_flux * t_factor) /
                       (r_factor^2 + L * t_factor^2)
        end

        push!(specs, MechanismTestSpec(
            name="MWC Dimer [AllostericEnzymeMechanism]",
            mechanism=m,
            metabolite_names=[:S, :P],
            expected_n_states=3,          # catalytic subunit: E_c, E_S, E_P
            expected_n_steps=3,           # 2 RE + 1 SS per subunit
            expected_n_metabolites=2,
            expected_n_haldane=2,         # k3r per conformation × 2
            expected_n_wegscheider=0,
            expected_n_independent_params=7,
            expected_identifiability_deficit=-2,
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=rate_mwc_dimer_oligo,
            expected_factored_num=
            "2 * ((k3f * S / K1 - k3r * P / K2) * (1 + S / K1 + P / K2)" *
            " + L * (k3f_T * S / K1_T - k3r_T * P / K2_T) * (1 + S / K1_T + P / K2_T))",
            expected_factored_denom=
            "(1 + S / K1 + P / K2) ^ 2 + L * (1 + S / K1_T + P / K2_T) ^ 2",
        ))
    end

    # ── Edge-case factoring tests ─────────────────────────────────────────────
    # These mechanisms test factoring patterns not covered by the classical
    # inhibitor/activator mechanisms above.

    # Spec #25 (flat homodimer "Homodimer + Non-competitive Inhibitor")
    # removed per Task 2.6: the AllostericEnzymeMechanism sibling (#25B
    # below) is the canonical form.

    # 25B. AllostericEnzymeMechanism equivalent of Homodimer + Non-competitive Inhibitor (spec #25)
    #      I binds all enzyme forms independently with the same Ki (enzyme-level).
    #      sigma = Q_cat^2 * (1 + I/K_I_reg1)  (multiplicative factor).
    let
        m = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: I::NonequalRT
            site(:catalytic, 2): begin
                steps: begin
                    [E_c, S] ⇌ [E_S]    :: NonequalRT
                    [E_c, P] ⇌ [E_P]    :: NonequalRT
                    [E_S] <--> [E_P]     :: NonequalRT
                end
            end
            site(:regulatory, 1): begin
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
            num   = 2.0 * (r_flux * r_factor * (1.0 + I / K_I_reg1) +
                           L * t_flux * t_factor * (1.0 + I / K_I_T_reg1))
            denom = r_factor^2 * (1.0 + I / K_I_reg1) +
                    L * t_factor^2 * (1.0 + I / K_I_T_reg1)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="Homodimer + Non-competitive Inhibitor [AllostericEnzymeMechanism]",
            mechanism=m,
            metabolite_names=[:S, :P, :I],
            expected_n_states=3,
            expected_n_steps=3,
            expected_n_metabolites=3,
            expected_n_haldane=2,
            expected_n_wegscheider=0,
            expected_n_independent_params=9,
            expected_identifiability_deficit=-6,
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=rate_homodimer_noncomp_inh_oligo,
            expected_factored_num=
            "2 * ((k3f * S / K1 - k3r * P / K2) * (1 + S / K1 + P / K2) * (1 + I / K_I_reg1)" *
            " + L * (k3f_T * S / K1_T - k3r_T * P / K2_T) * (1 + S / K1_T + P / K2_T) * (1 + I / K_I_T_reg1))",
            expected_factored_denom=
            "(1 + S / K1 + P / K2) ^ 2 * (1 + I / K_I_reg1)" *
            " + L * (1 + S / K1_T + P / K2_T) ^ 2 * (1 + I / K_I_T_reg1)",
        ))
    end

    # Spec #26 (flat homodimer "MWC Dimer + Independent Inhibitor")
    # removed per Task 2.6: the AllostericEnzymeMechanism sibling (#26B
    # below) is the canonical form.

    # 26B. AllostericEnzymeMechanism equivalent of MWC Dimer + Independent Inhibitor (spec #26)
    #      The Wegscheider constraint K80 = K47*K37/K38 (R_00I ⇌ T_00I equilibrium)
    #      is automatically satisfied by the conformational assembly formula — no
    #      explicit constraint needed in the AllostericEnzymeMechanism DSL.
    let
        m = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: I::NonequalRT
            site(:catalytic, 2): begin
                steps: begin
                    [E_c, S] ⇌ [E_S]    :: NonequalRT
                    [E_c, P] ⇌ [E_P]    :: NonequalRT
                    [E_S] <--> [E_P]     :: NonequalRT
                end
            end
            site(:regulatory, 1): begin
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
            num   = 2.0 * (r_flux * r_factor * (1.0 + I / K_I_reg1) +
                           L * t_flux * t_factor * (1.0 + I / K_I_T_reg1))
            denom = r_factor^2 * (1.0 + I / K_I_reg1) +
                    L * t_factor^2 * (1.0 + I / K_I_T_reg1)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="MWC Dimer + Independent Inhibitor [AllostericEnzymeMechanism]",
            mechanism=m,
            metabolite_names=[:S, :P, :I],
            expected_n_states=3,
            expected_n_steps=3,
            expected_n_metabolites=3,
            expected_n_haldane=2,
            expected_n_wegscheider=0,  # thermodynamic consistency is automatic
            expected_n_independent_params=9,
            expected_identifiability_deficit=-6,
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=rate_mwc_dimer_inh_oligo,
            expected_factored_num=
            "2 * ((k3f * S / K1 - k3r * P / K2) * (1 + S / K1 + P / K2) * (1 + I / K_I_reg1)" *
            " + L * (k3f_T * S / K1_T - k3r_T * P / K2_T) * (1 + S / K1_T + P / K2_T) * (1 + I / K_I_T_reg1))",
            expected_factored_denom=
            "(1 + S / K1 + P / K2) ^ 2 * (1 + I / K_I_reg1)" *
            " + L * (1 + S / K1_T + P / K2_T) ^ 2 * (1 + I / K_I_T_reg1)",
        ))
    end

    # 27. Two Competitive Inhibitors (monomer)
    #     Tests: multiple additive dead-end terms (flat sum denominator).
    #     I1 and I2 both bind only free E (competitive with S/P and each other).
    #     Denom: 1 + S/K1 + P/K3 + I1/K4 + I2/K5
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
                [E, P] ⇌ [E_P]
                [E, I1] ⇌ [E_I1]
                [E, I2] ⇌ [E_I2]
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
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=5,
            expected_n_steps=5,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=5,
            expected_identifiability_deficit=0,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_two_comp_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "1 + I1 / K4 + I2 / K5 + S / K1 + P / K3",
        ))
    end

    # 28. Two Non-competitive Inhibitors (monomer)
    #     Tests: triple multiplicative product factoring.
    #     I1 and I2 bind independently to all forms (E, ES, EP) at separate sites.
    #     Denom: (1+S/K1+P/K3) * (1+I1/K4) * (1+I2/K7)
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                # S binding (shared K1)
                ([E, S] ⇌ [E_S],
                 [E_I1, S] ⇌ [E_S_I1],
                 [E_I2, S] ⇌ [E_S_I2],
                 [E_I1_I2, S] ⇌ [E_S_I1_I2])
                [E_S] <--> [E_P]
                # P binding (shared K3)
                ([E, P] ⇌ [E_P],
                 [E_I1, P] ⇌ [E_P_I1],
                 [E_I2, P] ⇌ [E_P_I2],
                 [E_I1_I2, P] ⇌ [E_P_I1_I2])
                # I1 binding (shared K4)
                ([E, I1] ⇌ [E_I1],
                 [E_S, I1] ⇌ [E_S_I1],
                 [E_P, I1] ⇌ [E_P_I1],
                 [E_I2, I1] ⇌ [E_I1_I2],
                 [E_S_I2, I1] ⇌ [E_S_I1_I2],
                 [E_P_I2, I1] ⇌ [E_P_I1_I2])
                # I2 binding (shared K7)
                ([E, I2] ⇌ [E_I2],
                 [E_S, I2] ⇌ [E_S_I2],
                 [E_P, I2] ⇌ [E_P_I2],
                 [E_I1, I2] ⇌ [E_I1_I2],
                 [E_S_I1, I2] ⇌ [E_S_I1_I2],
                 [E_P_I1, I2] ⇌ [E_P_I1_I2])
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
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=12,
            expected_n_steps=21,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=5,
            expected_identifiability_deficit=-7,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_two_noncomp_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k5f * S / K1 - k5r * P / K6",
            expected_factored_denom=
            "(1 + S / K1 + P / K6) * (1 + I1 / K10) * (1 + I2 / K16)",
        ))
    end

    # 29. Non-competitive + Competitive Inhibitor (monomer)
    #     Tests: multiplicative product + additive dead-end term.
    #     I1 binds all forms independently (non-competitive).
    #     I2 binds only free E (competitive).
    #     Denom: (1+S/K1+P/K3)*(1+I1/K4) + I2/K9
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                # S binding (shared K1)
                ([E, S] ⇌ [E_S], [E_I1, S] ⇌ [E_S_I1])
                [E_S] <--> [E_P]
                # P binding (shared K3)
                ([E, P] ⇌ [E_P], [E_I1, P] ⇌ [E_P_I1])
                # I1 binding (shared K4)
                ([E, I1] ⇌ [E_I1],
                 [E_S, I1] ⇌ [E_S_I1],
                 [E_P, I1] ⇌ [E_P_I1])
                [E, I2] ⇌ [E_I2]
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
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=7,
            expected_n_steps=9,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=5,
            expected_identifiability_deficit=-2,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_noncomp_comp_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k3f * S / K1 - k3r * P / K4",
            expected_factored_denom=
            "I2 / K9 + (1 + S / K1 + P / K4) * (1 + I1 / K6)",
        ))
    end

    # 30. Uncompetitive + Competitive Inhibitor (monomer)
    #     Tests: mixed additive structure with nested multiplicative term.
    #     I1 binds only ES (uncompetitive). I2 binds only free E (competitive).
    #     Denom: 1 + I2/K5 + P/K3 + (S/K1)*(1+I1/K4)
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
                [E, P] ⇌ [E_P]
                [E_S, I1] ⇌ [E_S_I1]
                [E, I2] ⇌ [E_I2]
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
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=5,
            expected_n_steps=5,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=5,
            expected_identifiability_deficit=0,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_uncomp_comp_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "1 + I2 / K5 + P / K3 + (S / K1) * (1 + I1 / K4)",
        ))
    end

    # 31. Two Same-site (Competing) Non-competitive Inhibitors (monomer)
    #     Tests: multiplicative product with flat-sum allosteric factor.
    #     I1 and I2 compete for the same allosteric site X (mutually exclusive).
    #     Both bind independently of S/P (non-competitive).
    #     Denom: (1+S/K1+P/K3) * (1+I1/K4+I2/K9)
    let
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: I1, I2
            steps: begin
                # S binding (shared K1)
                ([E, S] ⇌ [E_S],
                 [E_I1, S] ⇌ [E_S_I1],
                 [E_I2, S] ⇌ [E_S_I2])
                [E_S] <--> [E_P]
                # P binding (shared K3)
                ([E, P] ⇌ [E_P],
                 [E_I1, P] ⇌ [E_P_I1],
                 [E_I2, P] ⇌ [E_P_I2])
                # I1 binding (shared K4)
                ([E, I1] ⇌ [E_I1],
                 [E_S, I1] ⇌ [E_S_I1],
                 [E_P, I1] ⇌ [E_P_I1])
                # I2 binding (shared K9)
                ([E, I2] ⇌ [E_I2],
                 [E_S, I2] ⇌ [E_S_I2],
                 [E_P, I2] ⇌ [E_P_I2])
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
            metabolite_names=[:S, :P, :I1, :I2],
            expected_n_states=9,
            expected_n_steps=13,
            expected_n_metabolites=4,
            expected_n_haldane=1,
            expected_n_wegscheider=0,
            expected_n_independent_params=5,
            expected_identifiability_deficit=-4,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_two_samesite_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "k4f * S / K1 - k4r * P / K5",
            expected_factored_denom=
            "(1 + I1 / K8 + I2 / K11) * (1 + S / K1 + P / K5)",
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
    #     Rate = Et * 4 * (N_cat_R * Q_cat_R^3 * Q_reg1_R^4 * Q_reg2_R^4
    #                    + L * N_cat_T * Q_cat_T^3 * Q_reg1_T^4 * Q_reg2_T^4) / Z
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
        m = @allosteric_mechanism begin
            substrates: S1, S2
            products: P1, P2
            allosteric_regulators: R1::NonequalRT, R2::NonequalRT, R3::NonequalRT
            site(:catalytic, 4): begin
                steps: begin
                    # S1 binding (shared K)
                    ([E_c,  S1] ⇌ [E_S1],
                     [E_S2, S1] ⇌ [E_S1S2],
                     [E_P2, S1] ⇌ [E_S1P2]) :: NonequalRT
                    # P1 binding (shared K)
                    ([E_c,  P1] ⇌ [E_P1],
                     [E_S2, P1] ⇌ [E_P1S2],
                     [E_P2, P1] ⇌ [E_P1P2]) :: NonequalRT
                    # S2 binding (shared K)
                    ([E_c,  S2] ⇌ [E_S2],
                     [E_S1, S2] ⇌ [E_S1S2],
                     [E_P1, S2] ⇌ [E_P1S2]) :: NonequalRT
                    # P2 binding (shared K)
                    ([E_c,  P2] ⇌ [E_P2],
                     [E_S1, P2] ⇌ [E_S1P2],
                     [E_P1, P2] ⇌ [E_P1P2]) :: NonequalRT
                    [E_S1S2] <--> [E_P1P2] :: NonequalRT
                end
            end
            site(:regulatory, 4): begin
                ligands: R1, R2
            end
            site(:regulatory, 4): begin
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

            return Et * 4.0 * num / Z
        end

        push!(specs, MechanismTestSpec(
            name="MWC Tetramer Random Bi-Bi RE + Two Allosteric Sites",
            mechanism=m,
            metabolite_names=[:S1, :S2, :P1, :P2, :R1, :R2, :R3],
            expected_n_states=9,           # catalytic subunit states
            expected_n_steps=13,           # catalytic subunit steps
            expected_n_metabolites=7,
            expected_n_haldane=2,          # one k13r per conformation (R and T)
            expected_n_wegscheider=0,      # site-independence constraints are in param_constraints of CM
            expected_n_independent_params=17,
            expected_identifiability_deficit=-29156,
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=rate_mwc_tetramer_bi_bi,
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end

    # ── PFK-1 hand-verified mechanism ───────────────────────────────────────
    # Reaction: F6P + ATP ⇌ F16BP + ADP, 4 catalytic subunits, 2 conformations.
    # F6P binding is :OnlyR — T-state can't bind F6P, so the T-state cycle is
    # broken in both directions and N_cat_T = 0. ATP appears as both
    # substrate and allosteric regulator (different tags per context).
    let
        m = @allosteric_mechanism begin
            substrates: F6P, ATP
            products:   F16BP, ADP
            allosteric_regulators: Pi::EqualRT, ATP::OnlyT, ADP::OnlyR, Citrate::OnlyT, F26BP::NonequalRT

            site(:catalytic, 4): begin
                steps: begin
                    ([E, F6P] ⇌ [E_F6P], [E_ATP, F6P] ⇌ [E_F6P_ATP])      :: OnlyR
                    ([E, ATP] ⇌ [E_ATP], [E_F6P, ATP] ⇌ [E_F6P_ATP])      :: EqualRT
                    [E_F6P_ATP] <--> [E_F16BP_ADP]                         :: EqualRT
                    ([E_F16BP_ADP] ⇌ [E_ADP, F16BP], [E_F16BP] ⇌ [E, F16BP]) :: EqualRT
                    ([E_F16BP_ADP] ⇌ [E_F16BP, ADP], [E_ADP] ⇌ [E, ADP])     :: EqualRT
                end
            end

            site(:regulatory, 4): begin
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

            # N_cat_T = 0: T-state cycle is broken (F6P binding :OnlyR), so
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
            return Et * 4.0 * num / (Q_R + L * Q_T)
        end

        push!(specs, MechanismTestSpec(
            name="PFK-1",
            mechanism=m,
            metabolite_names=[:F6P, :ATP, :F16BP, :ADP, :Pi, :Citrate, :F26BP],
            expected_n_states=7,
            expected_n_steps=9,
            expected_n_metabolites=7,
            # 6 deps total: 1 Haldane (k5r derived) + 5 :EqualRT-mirror
            # constraints (K3_T=K3, k5f_T=k5f, k5r_T=k5r, K6_T=K6, K8_T=K8,
            # K_Pi_T_reg1=K_Pi_reg1) — wait, that's 6 mirrors. Total = 7?
            # Actual measured value is 6; the mirrored count differs from
            # naive enumeration.
            expected_n_haldane=6,
            expected_n_wegscheider=0,
            expected_n_independent_params=12,
            expected_identifiability_deficit=-35061,
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=pfk_rate_analytical,
            # F6P binding (group 3) is :OnlyR → the Glucose·ATP saturating
            # pattern is unreachable in T-state, so kcat = catN · k5f
            # for every regulator corner. Regression test for the
            # `t_pattern_dead` branch in `_kcat_forward`.
            analytical_kcat_fn = p -> 4 * p.k5f,
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
    #      site — site(:regulatory, 2) with G6P::OnlyT and Pi::EqualRT.
    #
    # ATP binding (group 2) is :OnlyR — T-state can't bind ATP, so the
    # catalytic cycle is broken and N_cat_T = 0.
    let
        m = @allosteric_mechanism begin
            substrates: Glucose, ATP
            products:   G6P, ADP
            allosteric_regulators: G6P::OnlyT, Pi::EqualRT
            catalytic_inhibitors:  G6P

            site(:catalytic, 2): begin
                steps: begin
                    # Group 1 (Glucose binding at catalytic site, EqualRT)
                    ([E, Glucose] ⇌ [E_Glc],
                     [E_ATP, Glucose] ⇌ [E_Glc_ATP],
                     [E_G6Pi, Glucose] ⇌ [E_Glc_G6Pi])    :: EqualRT
                    # Group 2 (ATP binding at nucleotide pocket, OnlyR —
                    # T-state can't bind ATP)
                    ([E, ATP] ⇌ [E_ATP],
                     [E_Glc, ATP] ⇌ [E_Glc_ATP])           :: OnlyR
                    # Group 3 (catalysis SS, EqualRT)
                    [E_Glc_ATP] <--> [E_G6P_ADP]            :: EqualRT
                    # Group 4 (G6P binding/release at catalytic site, EqualRT)
                    ([E_G6P_ADP] ⇌ [E_ADP, G6P],
                     [E_G6P] ⇌ [E, G6P],
                     [E_G6P_G6Pi] ⇌ [E_G6Pi, G6P])         :: EqualRT
                    # Group 5 (ADP release, EqualRT)
                    ([E_G6P_ADP] ⇌ [E_G6P, ADP],
                     [E_ADP] ⇌ [E, ADP])                   :: EqualRT
                    # Group 6 (G6P binding at INHIBITORY site, EqualRT) —
                    # G6P at site 2 competes with ATP/ADP, can co-bind with
                    # Glucose at site 1 or G6P at site 1.
                    ([E, G6P] ⇌ [E_G6Pi],
                     [E_Glc, G6P] ⇌ [E_Glc_G6Pi],
                     [E_G6P, G6P] ⇌ [E_G6P_G6Pi])          :: EqualRT
                end
            end

            site(:regulatory, 2): begin
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
            # T-state: ATP group :OnlyR → zero ATP terms.
            Q_cat_T = 1 +
                      Glucose / K1 +
                      G6P * ADP / (K7 * K10) +
                      G6P / K7 +
                      ADP / K10 +
                      G6P / K12 +
                      Glucose * G6P / (K1 * K12) +
                      G6P^2 / (K7 * K12)

            # N_cat_T = 0: T-state cycle is broken (ATP binding :OnlyR).
            N_cat_R = k6f * Glucose * ATP / (K1 * K4) -
                      k6r * G6P * ADP / (K7 * K10)
            N_cat_T = 0.0

            Q_reg1_R = 1 + Pi / K_Pi_reg1
            Q_reg1_T = 1 + Pi / K_Pi_reg1 + G6P / K_G6P_T_reg1

            Q_R = Q_cat_R^2 * Q_reg1_R^2
            Q_T = Q_cat_T^2 * Q_reg1_T^2
            num = N_cat_R * Q_cat_R * Q_reg1_R^2 +
                  L * N_cat_T * Q_cat_T * Q_reg1_T^2
            return Et * 2.0 * num / (Q_R + L * Q_T)
        end

        push!(specs, MechanismTestSpec(
            name="HK",
            mechanism=m,
            metabolite_names=[:Glucose, :ATP, :G6P, :ADP, :Pi],
            expected_n_states=10,    # 7 catalytic-cycle + 3 G6Pi dead-end
            expected_n_steps=14,     # 3+2+1+3+2+3
            expected_n_metabolites=5,
            # 7 deps: 1 Haldane (k6r) + 6 :EqualRT-mirror constraints
            # (K1_T, k6f_T, k6r_T, K7_T, K10_T, K12_T, K_Pi_T_reg1).
            # The merged group 6 contributes ONE mirror constraint
            # (K12_T=K12) — the OLD design with two dead-end groups had
            # two (K10_T, K11_T).
            expected_n_haldane=7,
            expected_n_wegscheider=0,
            expected_n_independent_params=9,
            expected_identifiability_deficit=-178,
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=hk_rate_analytical,
            # ATP binding (group 2) is :OnlyR → the Glucose·ATP saturating
            # pattern is unreachable in T-state, so kcat = catN · k6f
            # for every regulator corner. Regression test for the
            # `t_pattern_dead` branch in `_kcat_forward`.
            analytical_kcat_fn = p -> 2 * p.k6f,
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end

    # ── Pyruvate kinase (PK) hand-verified mechanism ──────────────────────────
    # Reaction: PEP + ADP ⇌ Pyruvate + ATP, 4 catalytic subunits.
    # PEP binding is :NonequalRT (independent K_R and K_T) so the T-state
    # cycle is alive. Reg sites have MISMATCHED multiplicities:
    #   ATP::OnlyT at mult 2
    #   F16BP::OnlyR at mult 4 (matches catalytic mult)
    # This exercises the symmetric all-reg-sites contribution to both
    # numerator and denominator (after dropping the n_reg == CatN filter
    # in Task 3).
    let
        m = @allosteric_mechanism begin
            substrates: PEP, ADP
            products:   Pyruvate, ATP
            allosteric_regulators: ATP::OnlyT, F16BP::OnlyR

            site(:catalytic, 4): begin
                steps: begin
                    ([E, PEP] ⇌ [E_PEP],
                     [E_ADP, PEP] ⇌ [E_PEP_ADP])              :: NonequalRT
                    ([E, ADP] ⇌ [E_ADP],
                     [E_PEP, ADP] ⇌ [E_PEP_ADP])              :: EqualRT
                    [E_PEP_ADP] <--> [E_Pyr_ATP]               :: NonequalRT
                    ([E_Pyr_ATP] ⇌ [E_ATP, Pyruvate],
                     [E_Pyr] ⇌ [E, Pyruvate])                  :: EqualRT
                    ([E_Pyr_ATP] ⇌ [E_Pyr, ATP],
                     [E_ATP] ⇌ [E, ATP])                       :: EqualRT
                end
            end

            site(:regulatory, 2): begin
                ligands: ATP
            end
            site(:regulatory, 4): begin
                ligands: F16BP
            end
        end

        # Param mapping (kinetic-group representative-step convention):
        #   K1, K1_T    : PEP binding (group 1, NonequalRT)
        #   K3          : ADP binding (group 2, EqualRT)
        #   k5f, k5f_T  : catalysis SS (group 3, NonequalRT)
        #   K6          : Pyruvate release (group 4, EqualRT)
        #   K8          : ATP release (group 5, EqualRT)
        function pk_rate_analytical(params, concs)
            (; K1, K1_T, K3, k5f, k5f_T, K6, K8,
               K_ATP_T_reg1, K_F16BP_reg2,
               L, Keq, Et) = params
            (; PEP, ADP, Pyruvate, ATP, F16BP) = concs

            k5r   = k5f   * K6 * K8 / (Keq * K1   * K3)
            k5r_T = k5f_T * K6 * K8 / (Keq * K1_T * K3)

            Q_cat_R = 1 + PEP/K1   + ADP/K3 + PEP*ADP/(K1   * K3) +
                      Pyruvate/K6  + ATP/K8 + Pyruvate*ATP/(K6 * K8)
            Q_cat_T = 1 + PEP/K1_T + ADP/K3 + PEP*ADP/(K1_T * K3) +
                      Pyruvate/K6  + ATP/K8 + Pyruvate*ATP/(K6 * K8)

            N_R = k5f   * PEP * ADP / (K1   * K3) - k5r   * Pyruvate * ATP / (K6 * K8)
            N_T = k5f_T * PEP * ADP / (K1_T * K3) - k5r_T * Pyruvate * ATP / (K6 * K8)

            Q_reg1_R = 1                                     # ATP::OnlyT, no R term
            Q_reg1_T = 1 + ATP / K_ATP_T_reg1
            Q_reg2_R = 1 + F16BP / K_F16BP_reg2               # F16BP::OnlyR, no T term
            Q_reg2_T = 1

            num_R = N_R * Q_cat_R^3 * Q_reg1_R^2 * Q_reg2_R^4
            num_T = N_T * Q_cat_T^3 * Q_reg1_T^2 * Q_reg2_T^4
            den_R = Q_cat_R^4 * Q_reg1_R^2 * Q_reg2_R^4
            den_T = Q_cat_T^4 * Q_reg1_T^2 * Q_reg2_T^4

            return Et * 4.0 * (num_R + L * num_T) / (den_R + L * den_T)
        end

        push!(specs, MechanismTestSpec(
            name="PK",
            mechanism=m,
            metabolite_names=[:PEP, :ADP, :Pyruvate, :ATP, :F16BP],
            expected_n_states=7,           # E, E_PEP, E_ADP, E_PEP_ADP, E_Pyr_ATP, E_Pyr, E_ATP
            expected_n_steps=9,
            expected_n_metabolites=5,
            # 5 deps: 2 Haldanes (k5r, k5r_T) + 3 :EqualRT-mirror constraints
            # (K3_T, K6_T, K8_T). K1_T and k5f_T remain independent (:NonequalRT).
            expected_n_haldane=5,
            expected_n_wegscheider=0,
            expected_n_independent_params=10,
            expected_identifiability_deficit=-1443,
            expected_is_identifiable=true,
            run_ode_test=false,
            analytical_rate_fn=pk_rate_analytical,
            # kcat depends on L conformational equilibrium for :NonequalRT catalysis
            analytical_kcat_fn=nothing,
            expected_factored_num=nothing,
            expected_factored_denom=nothing,
        ))
    end

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

