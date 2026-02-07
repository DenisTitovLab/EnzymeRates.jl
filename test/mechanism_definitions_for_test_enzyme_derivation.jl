# Test specifications for all enzyme mechanisms
# Each mechanism is defined inline with its expected properties for easy side-by-side comparison

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
    mechanism::EnzymeMechanism            # The mechanism instance
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
    analytical_rate_fn::Union{Function, Nothing} = nothing
end

# ── Mechanism test specifications ───────────────────────────────────────────

function build_mechanism_test_specs()
    specs = MechanismTestSpec[]

    # 1. Uni-Uni (simplest): E + S ⇌ ES ⇌ E + P
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products:   P[C]
                enzymes:    E, ES[C]
            end
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [E, P]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Uni-Uni",
            mechanism = m,
            metabolite_names = [:S, :P],
            expected_n_states = 2,
            expected_n_steps = 2,
            expected_n_metabolites = 2,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 3,
            expected_identifiability_deficit = 0,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> begin
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
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
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
            name = "Segel Uni Uni",
            mechanism = m,
            metabolite_names = [:A, :P],
            expected_n_states = 3,
            expected_n_steps = 3,
            expected_n_metabolites = 2,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 5,
            expected_identifiability_deficit = 2,
            expected_is_identifiable = false,
            analytical_rate_fn = (p, c) -> rate_uni_uni(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 3. Segel Iso Uni Uni (new): E + A ⇌ EA ⇌ EP ⇌ F + P, F ⇌ E
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-45
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C], F
            end
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
            name = "Segel Iso Uni Uni",
            mechanism = m,
            metabolite_names = [:A, :P],
            expected_n_states = 4,
            expected_n_steps = 4,
            expected_n_metabolites = 2,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 7,
            expected_identifiability_deficit = 3,
            expected_is_identifiable = false,
            analytical_rate_fn = (p, c) -> rate_iso_uni_uni(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 4. Segel Ordered Uni Bi (replaces Seq Uni-Bi): E + A ⇌ (EA≡EPQ) ⇌ EQ + P ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-60
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[CN]
                products:   P[C], Q[N]
                enzymes:    E, EAEPQ[CN], EQ[N]
            end
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
            name = "Segel Ordered Uni Bi",
            mechanism = m,
            metabolite_names = [:A, :P, :Q],
            expected_n_states = 3,
            expected_n_steps = 3,
            expected_n_metabolites = 3,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 5,
            expected_identifiability_deficit = -1,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_ordered_uni_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 5. Segel Ordered Bi Bi (replaces Seq Bi-Bi): E + A ⇌ EA + B ⇌ (EAB≡EPQ) ⇌ EQ + P ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-87
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products:   P[C], Q[N]
                enzymes:    E, EA[C], EABEPQ[CN], EQ[N]
            end
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
            name = "Segel Ordered Bi Bi",
            mechanism = m,
            metabolite_names = [:A, :B, :P, :Q],
            expected_n_states = 4,
            expected_n_steps = 4,
            expected_n_metabolites = 4,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 7,
            expected_identifiability_deficit = -4,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_ordered_bi_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 6. Segel Theorell-Chance Bi Bi (new): E + A ⇌ EA + B ⇌ EQ + P ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-122
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products:   P[C], Q[N]
                enzymes:    E, EA[C], EQ[N]
            end
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
            name = "Segel Theorell-Chance Bi Bi",
            mechanism = m,
            metabolite_names = [:A, :B, :P, :Q],
            expected_n_states = 3,
            expected_n_steps = 3,
            expected_n_metabolites = 4,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 5,
            expected_identifiability_deficit = -4,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_theorell_chance_bi_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 7. Segel Ping Pong Bi Bi (replaces Ping-Pong Bi-Bi):
    #    E + A ⇌ (EA≡FP) ⇌ F + P, F + B ⇌ (FB≡EQ) ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-140
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[CX], B[N]
                products:   P[C], Q[NX]
                enzymes:    E, EAFP[CX], F[X], FBEQ[NX]
            end
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
            name = "Segel Ping Pong Bi Bi",
            mechanism = m,
            metabolite_names = [:A, :B, :P, :Q],
            expected_n_states = 4,
            expected_n_steps = 4,
            expected_n_metabolites = 4,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 7,
            expected_identifiability_deficit = -1,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_ping_pong_bi_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 8. Segel Ordered Ter Bi (replaces Seq Ter-Bi):
    #     E + A ⇌ EA + B ⇌ EAB + C ⇌ (EABC≡EPQ) ⇌ EQ + P ⇌ E + Q
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-195
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[X], B[Y], C[Z]
                products:   P[XZ], Q[Y]
                enzymes:    E, EA[X], EAB[XY], EABCEPQ[XYZ], EQ[Y]
            end
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
            name = "Segel Ordered Ter Bi",
            mechanism = m,
            metabolite_names = [:A, :B, :C, :P, :Q],
            expected_n_states = 5,
            expected_n_steps = 5,
            expected_n_metabolites = 5,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 9,
            expected_identifiability_deficit = -9,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_ordered_ter_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 9. Segel Ordered Ter Ter (replaces Seq Ter-Ter):
    #     E + A ⇌ EA + B ⇌ EAB + C ⇌ (EABC≡EPQR) ⇌ EQR + P ⇌ ER + Q ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-261
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[X], B[Y], C[Z]
                products:   P[X], Q[Y], R[Z]
                enzymes:    E, EA[X], EAB[XY], EABCEPQR[XYZ], EQR[YZ], ER[Z]
            end
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
            name = "Segel Ordered Ter Ter",
            mechanism = m,
            metabolite_names = [:A, :B, :C, :P, :Q, :R],
            expected_n_states = 6,
            expected_n_steps = 6,
            expected_n_metabolites = 6,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 11,
            expected_identifiability_deficit = -16,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_ordered_ter_ter(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 10. Segel Bi Uni Uni Uni Ping Pong Ter Bi (new):
    #     E + A ⇌ EA + B ⇌ (EAB≡FP) ⇌ F + P, F + C ⇌ (FC≡EQ) ⇌ E + Q
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-228
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[XW], B[Y], C[Z]
                products:   P[X], Q[WYZ]
                enzymes:    E, EA[XW], EABFP[XWY], F[WY], FCEQ[WYZ]
            end
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
            name = "Segel Bi Uni Uni Uni PP Ter Bi",
            mechanism = m,
            metabolite_names = [:A, :B, :C, :P, :Q],
            expected_n_states = 5,
            expected_n_steps = 5,
            expected_n_metabolites = 5,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 9,
            expected_identifiability_deficit = -5,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_bi_uni_uni_uni_ping_pong_ter_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 11. Random-order Bi-Bi (branched): Two substrate binding orders converge
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products:   P[C], Q[N]
                enzymes:    E, EA[C], EB[N], EAB[CN], EPQ[CN], EQ[N]
            end
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
            name = "Random-order Bi-Bi",
            mechanism = m,
            metabolite_names = [:A, :B, :P, :Q],
            expected_n_states = 6,
            expected_n_steps = 7,
            expected_n_metabolites = 4,
            expected_n_haldane = 1,
            expected_n_wegscheider = 1,
            expected_n_independent_params = 12,
            expected_identifiability_deficit = -16,
            expected_is_identifiable = true,
            analytical_rate_fn = nothing
        ))
    end

    # 12. Doubly Branched: E + S ⇌ EA, EA ⇌ EB, EA ⇌ EC, EB ⇌ E + P, EC ⇌ E + P
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products:   P[C]
                enzymes:    E, EA[C], EB[C], EC[C]
            end
            steps: begin
                [E, S] <--> [EA]
                [EA] <--> [EB]
                [EA] <--> [EC]
                [EB] <--> [E, P]
                [EC] <--> [E, P]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Doubly Branched",
            mechanism = m,
            metabolite_names = [:S, :P],
            expected_n_states = 4,
            expected_n_steps = 5,
            expected_n_metabolites = 2,
            expected_n_haldane = 1,
            expected_n_wegscheider = 1,
            expected_n_independent_params = 8,
            expected_identifiability_deficit = 5,
            expected_is_identifiable = false,
            analytical_rate_fn = nothing
        ))
    end

    # 13. Segel Bi Uni Uni Bi Ping Pong Ter Ter (new):
    #     E + A ⇌ EA + B ⇌ (EAB≡FP) ⇌ F + P, F + C ⇌ (FC≡EQR) ⇌ ER + Q ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-278
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[XW], B[Y], C[Z]
                products:   P[X], Q[Y], R[WZ]
                enzymes:    E, EA[XW], EABFP[XWY], F[WY], FCEQR[WYZ], ER[WZ]
            end
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
            name = "Segel Bi Uni Uni Bi PP Ter Ter",
            mechanism = m,
            metabolite_names = [:A, :B, :C, :P, :Q, :R],
            expected_n_states = 6,
            expected_n_steps = 6,
            expected_n_metabolites = 6,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 11,
            expected_identifiability_deficit = -11,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_bi_uni_uni_bi_ping_pong_ter_ter(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 14. Segel Bi Bi Uni Uni Ping Pong Ter Ter (new):
    #     E + A ⇌ EA + B ⇌ (EAB≡FPQ) ⇌ FQ + P, FQ ⇌ F + Q, F + C ⇌ (FC≡ER) ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-288
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[XW], B[YZ], C[V]
                products:   P[X], Q[Z], R[VWY]
                enzymes:    E, EA[XW], EABFPQ[XWYZ], FQ[WYZ], F[WY], FCER[VWY]
            end
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
            name = "Segel Bi Bi Uni Uni PP Ter Ter",
            mechanism = m,
            metabolite_names = [:A, :B, :C, :P, :Q, :R],
            expected_n_states = 6,
            expected_n_steps = 6,
            expected_n_metabolites = 6,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 11,
            expected_identifiability_deficit = -11,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_bi_bi_uni_uni_ping_pong_ter_ter(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 15. Segel Hexa Uni Ping Pong (new):
    #     E + A ⇌ (EA≡FP) ⇌ F + P, F + B ⇌ (FB≡GQ) ⇌ G + Q, G + C ⇌ (GC≡ER) ⇌ E + R
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-308
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[XW], B[YV], C[Z]
                products:   P[X], Q[V], R[WYZ]
                enzymes:    E, EAFP[XW], F[W], FBGQ[WYV], G[WY], GCER[WYZ]
            end
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
            name = "Segel Hexa Uni Ping Pong",
            mechanism = m,
            metabolite_names = [:A, :B, :C, :P, :Q, :R],
            expected_n_states = 6,
            expected_n_steps = 6,
            expected_n_metabolites = 6,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 11,
            expected_identifiability_deficit = -6,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> rate_hexa_uni_ping_pong(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    return specs
end

const MECHANISM_TEST_SPECS = build_mechanism_test_specs()
