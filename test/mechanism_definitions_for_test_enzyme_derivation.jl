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

    # 1. Uni-Uni: E + S ⇌ ES ⇌ E + P
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
                f = [p.k1f * c.S, p.k2f]
                r = [p.k1r, p.k2r * c.P]
                _unicyclic_flux(f, r, p.Et)
            end
        ))
    end

    # 2. Seq Uni-Bi: E + S1 ⇌ ES1 ⇌ EP1P2 ⇌ EP2 + P1 ⇌ E + P2
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S1[CH]
                products:   P1[C], P2[H]
                enzymes:    E, ES1[CH], EP1P2[CH], EP2[H]
            end
            steps: begin
                [E, S1] <--> [ES1]
                [ES1] <--> [EP1P2]
                [EP1P2] <--> [EP2, P1]
                [EP2] <--> [E, P2]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Seq Uni-Bi",
            mechanism = m,
            metabolite_names = [:S1, :P1, :P2],
            expected_n_states = 4,
            expected_n_steps = 4,
            expected_n_metabolites = 3,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 7,
            expected_identifiability_deficit = 1,
            expected_is_identifiable = false,
            analytical_rate_fn = (p, c) -> begin
                f = [p.k1f * c.S1, p.k2f, p.k3f, p.k4f]
                r = [p.k1r, p.k2r, p.k3r * c.P1, p.k4r * c.P2]
                _unicyclic_flux(f, r, p.Et)
            end
        ))
    end

    # 3. Seq Bi-Uni: E + S1 ⇌ ES1 + S2 ⇌ ES1S2 ⇌ EP1 ⇌ E + P1
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S1[C], S2[H]
                products:   P1[CH]
                enzymes:    E, ES1[C], ES1S2[CH], EP1[CH]
            end
            steps: begin
                [E, S1] <--> [ES1]
                [ES1, S2] <--> [ES1S2]
                [ES1S2] <--> [EP1]
                [EP1] <--> [E, P1]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Seq Bi-Uni",
            mechanism = m,
            metabolite_names = [:S1, :S2, :P1],
            expected_n_states = 4,
            expected_n_steps = 4,
            expected_n_metabolites = 3,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 7,
            expected_identifiability_deficit = 1,
            expected_is_identifiable = false,
            analytical_rate_fn = (p, c) -> begin
                f = [p.k1f * c.S1, p.k2f * c.S2, p.k3f, p.k4f]
                r = [p.k1r, p.k2r, p.k3r, p.k4r * c.P1]
                _unicyclic_flux(f, r, p.Et)
            end
        ))
    end

    # 4. Seq Bi-Bi: E + S1 ⇌ ES1 + S2 ⇌ ES1S2 ⇌ EP1P2 ⇌ EP2 + P1 ⇌ E + P2
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S1[C], S2[H]
                products:   P1[C], P2[H]
                enzymes:    E, ES1[C], ES1S2[CH], EP1P2[CH], EP2[H]
            end
            steps: begin
                [E, S1] <--> [ES1]
                [ES1, S2] <--> [ES1S2]
                [ES1S2] <--> [EP1P2]
                [EP1P2] <--> [EP2, P1]
                [EP2] <--> [E, P2]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Seq Bi-Bi",
            mechanism = m,
            metabolite_names = [:S1, :S2, :P1, :P2],
            expected_n_states = 5,
            expected_n_steps = 5,
            expected_n_metabolites = 4,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 9,
            expected_identifiability_deficit = -2,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> begin
                f = [p.k1f * c.S1, p.k2f * c.S2, p.k3f, p.k4f, p.k5f]
                r = [p.k1r, p.k2r, p.k3r, p.k4r * c.P1, p.k5r * c.P2]
                _unicyclic_flux(f, r, p.Et)
            end
        ))
    end

    # 5. Seq Bi-Ter
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S1[C], S2[HN]
                products:   P1[C], P2[H], P3[N]
                enzymes:    E, ES1[C], ES1S2[CHN], EP1P2P3[CHN], EP2P3[HN], EP3[N]
            end
            steps: begin
                [E, S1] <--> [ES1]
                [ES1, S2] <--> [ES1S2]
                [ES1S2] <--> [EP1P2P3]
                [EP1P2P3] <--> [EP2P3, P1]
                [EP2P3] <--> [EP3, P2]
                [EP3] <--> [E, P3]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Seq Bi-Ter",
            mechanism = m,
            metabolite_names = [:S1, :S2, :P1, :P2, :P3],
            expected_n_states = 6,
            expected_n_steps = 6,
            expected_n_metabolites = 5,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 11,
            expected_identifiability_deficit = -7,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> begin
                f = [p.k1f * c.S1, p.k2f * c.S2, p.k3f, p.k4f, p.k5f, p.k6f]
                r = [p.k1r, p.k2r, p.k3r, p.k4r * c.P1, p.k5r * c.P2, p.k6r * c.P3]
                _unicyclic_flux(f, r, p.Et)
            end
        ))
    end

    # 6. Seq Ter-Bi
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S1[C], S2[H], S3[N]
                products:   P1[CH], P2[N]
                enzymes:    E, ES1[C], ES1S2[CH], ES1S2S3[CHN], EP1P2[CHN], EP2[N]
            end
            steps: begin
                [E, S1] <--> [ES1]
                [ES1, S2] <--> [ES1S2]
                [ES1S2, S3] <--> [ES1S2S3]
                [ES1S2S3] <--> [EP1P2]
                [EP1P2] <--> [EP2, P1]
                [EP2] <--> [E, P2]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Seq Ter-Bi",
            mechanism = m,
            metabolite_names = [:S1, :S2, :S3, :P1, :P2],
            expected_n_states = 6,
            expected_n_steps = 6,
            expected_n_metabolites = 5,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 11,
            expected_identifiability_deficit = -7,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> begin
                f = [p.k1f * c.S1, p.k2f * c.S2, p.k3f * c.S3, p.k4f, p.k5f, p.k6f]
                r = [p.k1r, p.k2r, p.k3r, p.k4r, p.k5r * c.P1, p.k6r * c.P2]
                _unicyclic_flux(f, r, p.Et)
            end
        ))
    end

    # 7. Seq Ter-Ter
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S1[C], S2[H], S3[N]
                products:   P1[C], P2[H], P3[N]
                enzymes:    E, ES1[C], ES1S2[CH], ES1S2S3[CHN], EP1P2P3[CHN], EP2P3[HN], EP3[N]
            end
            steps: begin
                [E, S1] <--> [ES1]
                [ES1, S2] <--> [ES1S2]
                [ES1S2, S3] <--> [ES1S2S3]
                [ES1S2S3] <--> [EP1P2P3]
                [EP1P2P3] <--> [EP2P3, P1]
                [EP2P3] <--> [EP3, P2]
                [EP3] <--> [E, P3]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Seq Ter-Ter",
            mechanism = m,
            metabolite_names = [:S1, :S2, :S3, :P1, :P2, :P3],
            expected_n_states = 7,
            expected_n_steps = 7,
            expected_n_metabolites = 6,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 13,
            expected_identifiability_deficit = -14,
            expected_is_identifiable = true,
            analytical_rate_fn = (p, c) -> begin
                f = [p.k1f * c.S1, p.k2f * c.S2, p.k3f * c.S3, p.k4f, p.k5f, p.k6f, p.k7f]
                r = [p.k1r, p.k2r, p.k3r, p.k4r, p.k5r * c.P1, p.k6r * c.P2, p.k7r * c.P3]
                _unicyclic_flux(f, r, p.Et)
            end
        ))
    end

    # 8. Ping-Pong Bi-Bi: E + A ⇌ EA ⇌ FP ⇌ F + P; F + B ⇌ FB ⇌ EQ ⇌ E + Q
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C2N], B[C3]
                products:   P[C2], Q[C3N]
                enzymes:    E, EA[C2N], FP[C2N], F[N], FB[C3N], EQ[C3N]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [FP]
                [FP] <--> [F, P]
                [F, B] <--> [FB]
                [FB] <--> [EQ]
                [EQ] <--> [E, Q]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Ping-Pong Bi-Bi",
            mechanism = m,
            metabolite_names = [:A, :P, :B, :Q],
            expected_n_states = 6,
            expected_n_steps = 6,
            expected_n_metabolites = 4,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 11,
            expected_identifiability_deficit = 3,
            expected_is_identifiable = false,
            analytical_rate_fn = nothing
        ))
    end

    # 9. Random-order Bi-Bi (branched): Two substrate binding orders converge
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

    # 10. Doubly Branched: E + S ⇌ EA, EA ⇌ EB, EA ⇌ EC, EB ⇌ E + P, EC ⇌ E + P
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

    # 11. Three-Step Isomerization: E + S ⇌ ES ⇌ ES' ⇌ E + P
    # 3 steps, 6 k's, 1 Haldane → 5 independent k's
    # Rate equation has same form as uni-uni: (a*S - b*P)/(1 + c*S + d*P)
    # So: 3 identifiable coefficients, deficit = 5 - 3 = 2
    let
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products:   P[C]
                enzymes:    E, ES[C], ES2[C]
            end
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [ES2]
                [ES2] <--> [E, P]
            end
        end
        push!(specs, MechanismTestSpec(
            name = "Three-Step Isomerization",
            mechanism = m,
            metabolite_names = [:S, :P],
            expected_n_states = 3,
            expected_n_steps = 3,
            expected_n_metabolites = 2,
            expected_n_haldane = 1,
            expected_n_wegscheider = 0,
            expected_n_independent_params = 5,
            expected_identifiability_deficit = 2,
            expected_is_identifiable = false,
            analytical_rate_fn = (p, c) -> begin
                f = [p.k1f * c.S, p.k2f, p.k3f]
                r = [p.k1r, p.k2r, p.k3r * c.P]
                _unicyclic_flux(f, r, p.Et)
            end
        ))
    end

    return specs
end

const MECHANISM_TEST_SPECS = build_mechanism_test_specs()
