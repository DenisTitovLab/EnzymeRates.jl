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
    mechanism::Any                        # EnzymeMechanism or OligomericEnzymeMechanism
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

    # Factored form validation (optional — for mechanisms with known factored forms)
    # Expected numerator/denominator from rate_equation_string output
    # When broken=true, uses @test_broken (known bugs, will alert when fixed)
    expected_factored_num::Union{String,Nothing} = nothing
    expected_factored_denom::Union{String,Nothing} = nothing
    factored_num_broken::Bool = false
    factored_denom_broken::Bool = false
end

# ── Mechanism test specifications ───────────────────────────────────────────

function build_mechanism_test_specs()
    specs = MechanismTestSpec[]

    # 1. Uni-Uni (simplest): E + S ⇌ ES ⇌ E + P
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                enzymes:E, ES[C]
            end
            steps:begin
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
            species:begin
                substrates:A[C]
                products:P[C]
                enzymes:E, EA[C], EP[C]
            end
            steps:begin
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
            species:begin
                substrates:A[C]
                products:P[C]
                enzymes:E, EA[C], EP[C], F
            end
            steps:begin
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
            species:begin
                substrates:A[CN]
                products:P[C], Q[N]
                enzymes:E, EAEPQ[CN], EQ[N]
            end
            steps:begin
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
            species:begin
                substrates:A[C], B[N]
                products:P[C], Q[N]
                enzymes:E, EA[C], EABEPQ[CN], EQ[N]
            end
            steps:begin
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
            analytical_rate_fn=(p, c) -> rate_ordered_bi_bi(merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 6. Segel Theorell-Chance Bi Bi (new): E + A ⇌ EA + B ⇌ EQ + P ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-122
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:A[C], B[N]
                products:P[C], Q[N]
                enzymes:E, EA[C], EQ[N]
            end
            steps:begin
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
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 7. Segel Ping Pong Bi Bi (replaces Ping-Pong Bi-Bi):
    #    E + A ⇌ (EA≡FP) ⇌ F + P, F + B ⇌ (FB≡EQ) ⇌ E + Q
    #    Reference: Segel, Enzyme Kinetics, Eq. IX-140
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:A[CX], B[N]
                products:P[C], Q[NX]
                enzymes:E, EAFP[CX], F[X], FBEQ[NX]
            end
            steps:begin
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
                    merge(p, (Etotal=p.Et,)), c)
        ))
    end

    # 8. Segel Ordered Ter Bi (replaces Seq Ter-Bi):
    #     E + A ⇌ EA + B ⇌ EAB + C ⇌ (EABC≡EPQ) ⇌ EQ + P ⇌ E + Q
    #     Reference: Segel, Enzyme Kinetics, Eq. IX-195
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:A[X], B[Y], C[Z]
                products:P[XZ], Q[Y]
                enzymes:E, EA[X], EAB[XY], EABCEPQ[XYZ], EQ[Y]
            end
            steps:begin
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
            species:begin
                substrates:A[X], B[Y], C[Z]
                products:P[X], Q[Y], R[Z]
                enzymes:E, EA[X], EAB[XY], EABCEPQR[XYZ], EQR[YZ], ER[Z]
            end
            steps:begin
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
            species:begin
                substrates:A[XW], B[Y], C[Z]
                products:P[X], Q[WYZ]
                enzymes:E, EA[XW], EABFP[XWY], F[WY], FCEQ[WYZ]
            end
            steps:begin
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
            species:begin
                substrates:A[C], B[N]
                products:P[C], Q[N]
                enzymes:E, EA[C], EB[N], EAB[CN], EPQ[CN], EQ[N]
            end
            steps:begin
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
            species:begin
                substrates:A[XW], B[Y], C[Z]
                products:P[X], Q[Y], R[WZ]
                enzymes:E, EA[XW], EABFP[XWY], F[WY], FCEQR[WYZ], ER[WZ]
            end
            steps:begin
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
            species:begin
                substrates:A[XW], B[YZ], C[V]
                products:P[X], Q[Z], R[VWY]
                enzymes:E, EA[XW], EABFPQ[XWYZ], FQ[WYZ], F[WY], FCER[VWY]
            end
            steps:begin
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
            species:begin
                substrates:A[XW], B[YV], C[Z]
                products:P[X], Q[V], R[WYZ]
                enzymes:E, EAFP[XW], F[W], FBGQ[WYV], G[WY], GCER[WYZ]
            end
            steps:begin
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
            species:begin
                substrates:A[C]
                products:P[C]
                enzymes:E, EA[C]
            end
            steps:begin
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
            analytical_rate_fn=(p, c) -> rate_re_uni_uni(merge(p, (Et=p.Et,)), c)
        ))
    end

    # 16. RE Ordered Bi-Bi: substrate binding is RE, catalysis and product release are SS
    #     E + A ⇌_RE EA, EA + B ⇌_RE EAB, (EAB≡EPQ) <-->_SS EQ + P, EQ <-->_SS E + Q
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:A[C], B[N]
                products:P[C], Q[N]
                enzymes:E, EA[C], EABEPQ[CN], EQ[N]
            end
            steps:begin
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
            species:begin
                substrates:A[C], B[N]
                products:P[C], Q[N]
                enzymes:E, EA[C], EB[N], EAB[CN], EPQ[CN], EQ[N]
            end
            steps:begin
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
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:R
                enzymes:E, ES[C], EP[C], ER[X]
            end
            steps:begin
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
            # Textbook: flat sum denominator (no Cartesian product structure)
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "1 + S / K1 + P / K3 + R / K4",
            factored_num_broken=false,
            factored_denom_broken=false,
        ))
    end

    # 19. Non-competitive inhibitor: R binds both free E and ES with same K
    #     Forms: E, E_S, E_P, E_R, E_S_R; SS: E_S↔E_P; K5=K4
    #     Denom factors as (1+R/K4)*(1+S/K1) + P/K3
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:R
                enzymes:E, E_S[C], E_P[C], E_R[X], E_S_R[CX]
            end
            steps:begin
                [E, S] ⇌ [E_S]        # K1
                [E_S] <--> [E_P]       # k2f, k2r (SS)
                [E, P] ⇌ [E_P]        # K3
                [E, R] ⇌ [E_R]        # K4
                [E_S, R] ⇌ [E_S_R]    # K5
                [E_R, S] ⇌ [E_S_R]    # K6
            end
            constraints:begin
                K5 = K4
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
            factored_num_broken=false,
            factored_denom_broken=false,
        ))
    end

    # 20. Uncompetitive inhibitor: R binds ES only (not free E)
    #     Forms: E, E_S, E_P, E_S_R; SS: E_S↔E_P; No extra constraints
    #     Denom: 1 + P/K3 + S/K1*(1+R/K4)
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:R
                enzymes:E, E_S[C], E_P[C], E_S_R[CX]
            end
            steps:begin
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
            factored_num_broken=false,
            factored_denom_broken=false,
        ))
    end

    # 21. Essential activator: R must bind before S can bind
    #     Forms: E, E_R, E_S_R, E_P_R; SS: E_S_R↔E_P_R
    #     Num: R/K4 * (k2f*S/K1 - k2r*P/K3)
    #     Denom: 1 + R/K4*(1+S/K1+P/K3)
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:R
                enzymes:E, E_R[X], E_S_R[CX], E_P_R[CX]
            end
            steps:begin
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
            expected_factored_num=
            "(R / K4) * (k2f * S / K1 - k2r * P / K3)",
            factored_num_broken=false,
            expected_factored_denom=
            "1 + (R / K4) * (1 + S / K1 + P / K3)",
            factored_denom_broken=false,
        ))
    end

    # 22. Non-essential activator (general modifier):
    #     R modifies catalysis but isn't required. Two parallel SS cycles.
    #     Forms: E, E_S, E_P, E_R, E_S_R, E_P_R; SS: E_S↔E_P, E_S_R↔E_P_R
    #     K8=K7, K9=K7 (R binding independent of S/P)
    #     K4=K1 and K6=K3 are implied by Wegscheider relations
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:R
                enzymes:E, E_S[C], E_P[C], E_R[X], E_S_R[CX], E_P_R[CX]
            end
            steps:begin
                [E, S] ⇌ [E_S]          # K1
                [E_S] <--> [E_P]         # k2f, k2r (SS)
                [E, P] ⇌ [E_P]          # K3
                [E_R, S] ⇌ [E_S_R]      # K4
                [E_S_R] <--> [E_P_R]     # k5f, k5r (SS)
                [E_R, P] ⇌ [E_P_R]      # K6
                [E, R] ⇌ [E_R]          # K7
                [E_S, R] ⇌ [E_S_R]      # K8
                [E_P, R] ⇌ [E_P_R]      # K9
            end
            constraints:begin
                K8 = K7
                K9 = K7
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
            expected_factored_num=
            "k2f * S / K1 - k2r * P / K3 + (R / K7) * (k5f * S / K1 - k5r * P / K3)",
            factored_num_broken=false,
            expected_factored_denom=
            "(1 + S / K1 + P / K3) * (1 + R / K7)",
            factored_denom_broken=false,
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
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:A, I
                enzymes:E, E_S[C], E_P[C], E_A[X],
                E_S_A[CX], E_P_A[CX], E_I[Y]
            end
            steps:begin
                [E, S] ⇌ [E_S]          # K1
                [E_S] <--> [E_P]         # k2f, k2r (SS)
                [E, P] ⇌ [E_P]          # K3
                [E_A, S] ⇌ [E_S_A]      # K4
                [E_S_A] <--> [E_P_A]     # k5f, k5r (SS)
                [E_A, P] ⇌ [E_P_A]      # K6
                [E, A] ⇌ [E_A]          # K7
                [E_S, A] ⇌ [E_S_A]      # K8
                [E_P, A] ⇌ [E_P_A]      # K9
                [E, I] ⇌ [E_I]          # K10
            end
            constraints:begin
                K8 = K7
                K9 = K7
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
            factored_num_broken=false,
            expected_factored_denom=
            "I / K10 + (1 + S / K1 + P / K3) * (1 + A / K7)",
            factored_denom_broken=false,
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
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                enzymes:R_00, T_00,
                R_S0[C], R_0S[C], R_SS[C2], R_SP[C2], R_PS[C2], R_P0[C], R_0P[C], R_PP[C2],
                T_S0[C], T_0S[C], T_SS[C2], T_SP[C2], T_PS[C2], T_P0[C], T_0P[C], T_PP[C2]
            end
            steps:begin
                # S binding to R (RE)
                [R_00, S] ⇌ [R_S0]         # K1 (K_R)
                [R_00, S] ⇌ [R_0S]         # K2
                [R_S0, S] ⇌ [R_SS]         # K3
                [R_0S, S] ⇌ [R_SS]         # K4
                [R_0P, S] ⇌ [R_SP]         # K5
                [R_P0, S] ⇌ [R_PS]         # K6
                # Catalysis on R (SS)
                [R_S0] <--> [R_P0]          # k7f, k7r
                [R_0S] <--> [R_0P]          # k8f, k8r
                [R_SS] <--> [R_SP]          # k9f, k9r
                [R_SS] <--> [R_PS]          # k10f, k10r
                [R_SP] <--> [R_PP]          # k11f, k11r
                [R_PS] <--> [R_PP]          # k12f, k12r
                # P binding to R (RE)
                [R_00, P] ⇌ [R_P0]         # K13 (K_P)
                [R_00, P] ⇌ [R_0P]         # K14
                [R_P0, P] ⇌ [R_PP]         # K15
                [R_0P, P] ⇌ [R_PP]         # K16
                [R_S0, P] ⇌ [R_SP]         # K17
                [R_0S, P] ⇌ [R_PS]         # K18
                # S binding to T (TE)
                [T_00, S] ⇌ [T_S0]         # K19 (K_T)
                [T_00, S] ⇌ [T_0S]         # K20
                [T_S0, S] ⇌ [T_SS]         # K21
                [T_0S, S] ⇌ [T_SS]         # K22
                [T_0P, S] ⇌ [T_SP]         # K23
                [T_P0, S] ⇌ [T_PS]         # K24
                # Catalysis on T (SS)
                [T_S0] <--> [T_P0]          # k25f, k25r
                [T_0S] <--> [T_0P]          # k26f, k26r
                [T_SS] <--> [T_SP]          # k27f, k27r
                [T_SS] <--> [T_PS]          # k28f, k28r
                [T_SP] <--> [T_PP]          # k29f, k29r
                [T_PS] <--> [T_PP]          # k30f, k30r
                # P binding to T (TE)
                [T_00, P] ⇌ [T_P0]         # K31 (K_TP)
                [T_00, P] ⇌ [T_0P]         # K32
                [T_P0, P] ⇌ [T_PP]         # K33
                [T_0P, P] ⇌ [T_PP]         # K34
                [T_S0, P] ⇌ [T_SP]         # K35
                [T_0S, P] ⇌ [T_PS]         # K36
                # R ↔ T isomerization (RE, K37 = L = [T]/[R])
                [R_00] ⇌ [T_00]            # K37
            end
            constraints:begin
                # S binding to R: all equal
                K2 = K1
                K3 = K1
                K4 = K1
                K5 = K1
                K6 = K1
                # Catalysis forward: all equal
                k8f = k7f
                k9f = k7f
                k10f = k7f
                k11f = k7f
                k12f = k7f
                # P binding to R: all equal
                K14 = K13
                K15 = K13
                K16 = K13
                K17 = K13
                K18 = K13
                # S binding to T: all equal
                K20 = K19
                K21 = K19
                K22 = K19
                K23 = K19
                K24 = K19
                # Catalysis forward: all equal
                k26f = k25f
                k27f = k25f
                k28f = k25f
                k29f = k25f
                k30f = k25f
                # P binding to T: all equal
                K32 = K31
                K33 = K31
                K34 = K31
                K35 = K31
                K36 = K31
            end
        end

        # MWC dimer rate equation:
        # sigma = (1 + S/K_R + P/K_P)^2 + L*(1 + S/K_T + P/K_T)^2
        # flux = 2*((k_cat_R*S/K_R - k_rev_R*P/K_P) + L*(k_cat_T*S/K_T - k_rev_T*P/K_T))
        # v = Et * flux / sigma
        # K37 = L (Ka convention, non-binding isomerization)
        function rate_mwc_dimer(params, concs)
            (; K1, k7f, k7r, k25f, k25r, K13, K19, K37, K31, Et) = params
            (; S, P) = concs
            r_flux = k7f * S / K1 - k7r * P / K13
            t_flux = k25f * S / K19 - k25r * P / K31
            r_factor = 1.0 + S / K1 + P / K13
            t_factor = 1.0 + S / K19 + P / K31
            num = 2.0 * (r_flux * r_factor + K37 * t_flux * t_factor)
            denom = r_factor^2 + K37 * t_factor^2
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="MWC Dimer",
            mechanism=m,
            metabolite_names=[:S, :P],
            expected_n_states=18,
            expected_n_steps=37,
            expected_n_metabolites=2,
            expected_n_haldane=12,
            expected_n_wegscheider=0,
            expected_n_independent_params=7,
            expected_identifiability_deficit=-2,
            expected_is_identifiable=true,
            analytical_rate_fn=rate_mwc_dimer,
            expected_factored_num=
            "2 * (1 + S / K1 + P / K13) * (k7f * S / K1 - k7r * P / K13) + 2 * K37 * (1 + S / K19 + P / K31) * (k25f * S / K19 - k25r * P / K31)",
            factored_num_broken=false,
            expected_factored_denom=
            "(1 + S / K1 + P / K13) ^ 2 + K37 * (1 + S / K19 + P / K31) ^ 2",
            factored_denom_broken=false,
        ))
    end

    # 24B. OligomericEnzymeMechanism equivalent of MWC Dimer (spec #24)
    #      2 catalytic sites × 2 conformations (R/T). No explicit equality constraints
    #      needed — symmetric subunits are captured by the site multiplicity.
    #      Conformational equilibrium: L (= K37 in the EnzymeMechanism above).
    if isdefined(EnzymeRates, :OligomericEnzymeMechanism)
        let
            m = @enzyme_mechanism begin
                metabolites: S[C], P[C]
                conformations: 2    # R (active) and T (tense)
                site(:catalytic, 2):begin
                    species: begin
                        substrates: S[C]
                        products: P[C]
                        enzymes: E_c, E_S[C], E_P[C]
                    end
                    steps: begin
                        [E_c, S] ⇌ [E_S]    # K1 (R), K1_T (T)
                        [E_c, P] ⇌ [E_P]    # K2 (R), K2_T (T)
                        [E_S] <--> [E_P]     # k3f, k3r Haldane (R); k3f_T, k3r_T Haldane (T)
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
                name="MWC Dimer [OligomericEnzymeMechanism]",
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
    end

    # ── Edge-case factoring tests ─────────────────────────────────────────────
    # These mechanisms test factoring patterns not covered by the classical
    # inhibitor/activator mechanisms above.

    # 25. Homodimer + Non-competitive Inhibitor
    #     Tests: power^2 × multiplicative linear factor.
    #     I binds all 9 homodimer forms independently (same Ki).
    #     Denom: (1+S/K1+P/K13)^2 * (1+I/K19)
    #     Num: 2*(1+S/K1+P/K13)*(k7f*S/K1 - k7r*P/K13)
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:I
                enzymes:E_00,
                E_0S[C], E_S0[C], E_0P[C], E_P0[C],
                E_SS[C2], E_SP[C2], E_PS[C2], E_PP[C2],
                E_00I[X],
                E_0SI[CX], E_S0I[CX], E_0PI[CX], E_P0I[CX],
                E_SSI[C2X], E_SPI[C2X], E_PSI[C2X], E_PPI[C2X]
            end
            steps:begin
                [E_00, S] ⇌ [E_0S]
                [E_00, S] ⇌ [E_S0]
                [E_0S, S] ⇌ [E_SS]
                [E_S0, S] ⇌ [E_SS]
                [E_0P, S] ⇌ [E_SP]
                [E_P0, S] ⇌ [E_PS]
                [E_S0] <--> [E_P0]
                [E_0S] <--> [E_0P]
                [E_SS] <--> [E_SP]
                [E_SS] <--> [E_PS]
                [E_SP] <--> [E_PP]
                [E_PS] <--> [E_PP]
                [E_00, P] ⇌ [E_P0]
                [E_00, P] ⇌ [E_0P]
                [E_P0, P] ⇌ [E_PP]
                [E_0P, P] ⇌ [E_PP]
                [E_S0, P] ⇌ [E_SP]
                [E_0S, P] ⇌ [E_PS]
                [E_00, I] ⇌ [E_00I]
                [E_S0, I] ⇌ [E_S0I]
                [E_0S, I] ⇌ [E_0SI]
                [E_SS, I] ⇌ [E_SSI]
                [E_P0, I] ⇌ [E_P0I]
                [E_0P, I] ⇌ [E_0PI]
                [E_SP, I] ⇌ [E_SPI]
                [E_PS, I] ⇌ [E_PSI]
                [E_PP, I] ⇌ [E_PPI]
                [E_00I, S] ⇌ [E_0SI]
                [E_00I, S] ⇌ [E_S0I]
                [E_0SI, S] ⇌ [E_SSI]
                [E_S0I, S] ⇌ [E_SSI]
                [E_0PI, S] ⇌ [E_SPI]
                [E_P0I, S] ⇌ [E_PSI]
                [E_00I, P] ⇌ [E_P0I]
                [E_00I, P] ⇌ [E_0PI]
                [E_P0I, P] ⇌ [E_PPI]
                [E_0PI, P] ⇌ [E_PPI]
                [E_S0I, P] ⇌ [E_SPI]
                [E_0SI, P] ⇌ [E_PSI]
            end
            constraints:begin
                K2 = K1
                K3 = K1
                K4 = K1
                K5 = K1
                K6 = K1
                k8f = k7f
                k9f = k7f
                k10f = k7f
                k11f = k7f
                k12f = k7f
                K14 = K13
                K15 = K13
                K16 = K13
                K17 = K13
                K18 = K13
                K20 = K19
                K21 = K19
                K22 = K19
                K23 = K19
                K24 = K19
                K25 = K19
                K26 = K19
                K27 = K19
                K28 = K1
                K29 = K1
                K30 = K1
                K31 = K1
                K32 = K1
                K33 = K1
                K34 = K13
                K35 = K13
                K36 = K13
                K37 = K13
                K38 = K13
                K39 = K13
            end
        end

        function rate_homodimer_noncomp_inh(p, c)
            (; K1, k7f, k7r, K13, K19, Et) = p
            (; S, P, I) = c
            flux = k7f * S / K1 - k7r * P / K13
            factor = 1.0 + S / K1 + P / K13
            return Et * 2.0 * flux * factor /
                   (factor^2 * (1.0 + I / K19))
        end

        push!(specs, MechanismTestSpec(
            name="Homodimer + Non-competitive Inhibitor",
            mechanism=m,
            metabolite_names=[:S, :P, :I],
            expected_n_states=18,
            expected_n_steps=39,
            expected_n_metabolites=3,
            expected_n_haldane=6,
            expected_n_wegscheider=0,
            expected_n_independent_params=4,
            expected_identifiability_deficit=-11,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_homodimer_noncomp_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "2 * (1 + S / K1 + P / K13) * (k7f * S / K1 - k7r * P / K13)",
            expected_factored_denom=
            "(1 + S / K1 + P / K13) ^ 2 * (1 + I / K19)",
        ))
    end

    # 25B. OligomericEnzymeMechanism equivalent of Homodimer + Non-competitive Inhibitor (spec #25)
    #      I binds all enzyme forms independently with the same Ki (enzyme-level).
    #      sigma = Q_cat^2 * (1 + I/K_I_reg1)  (multiplicative factor).
    if isdefined(EnzymeRates, :OligomericEnzymeMechanism)
        let
            m = @enzyme_mechanism begin
                metabolites: S[C], P[C], I[X]
                site(:catalytic, 2):begin
                    species: begin
                        substrates: S[C]
                        products: P[C]
                        enzymes: E_c, E_S[C], E_P[C]
                    end
                    steps: begin
                        [E_c, S] ⇌ [E_S]    # K1
                        [E_c, P] ⇌ [E_P]    # K2
                        [E_S] <--> [E_P]     # k3f, k3r via Haldane
                    end
                end
                site(:regulatory, 1):begin  # enzyme-level: sigma *= (1 + I/K_I_reg1)
                    ligands: I
                end
            end

            function rate_homodimer_noncomp_inh_oligo(params, concs)
                (; K1, K2, k3f, k3r, K_I_reg1, Et) = params
                (; S, P, I) = concs
                flux   = k3f * S / K1 - k3r * P / K2
                factor = 1.0 + S / K1 + P / K2
                return Et * 2.0 * flux * factor /
                           (factor^2 * (1.0 + I / K_I_reg1))
            end

            push!(specs, MechanismTestSpec(
                name="Homodimer + Non-competitive Inhibitor [OligomericEnzymeMechanism]",
                mechanism=m,
                metabolite_names=[:S, :P, :I],
                expected_n_states=3,
                expected_n_steps=3,
                expected_n_metabolites=3,
                expected_n_haldane=1,
                expected_n_wegscheider=0,
                expected_n_independent_params=4,
                expected_identifiability_deficit=-11,
                expected_is_identifiable=true,
                run_ode_test=false,
                analytical_rate_fn=rate_homodimer_noncomp_inh_oligo,
                expected_factored_num=
                "2 * (k3f * S / K1 - k3r * P / K2) * (1 + S / K1 + P / K2)",
                expected_factored_denom=
                "(1 + S / K1 + P / K2) ^ 2 * (1 + I / K_I_reg1)",
            ))
        end
    end

    # 26. MWC Dimer + Independent Inhibitor
    #     Tests: multi-group power^2 × per-group linear factor.
    #     I binds independently to all R and T forms (different Ki per conformation).
    #     Denom: (1+S/K1+P/K13)^2*(1+I/K38) + K37*(1+S/K19+P/K31)^2*(1+I/K47)
    #     Num: 2*[(1+S/K1+P/K13)*(k7f*S/K1-k7r*P/K13)
    #            + K37*(1+S/K19+P/K31)*(k25f*S/K19-k25r*P/K31)]
    #     K37 is Wegscheider-derived: K37 = K47*K80/K38
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:I
                enzymes:R_00, R_S0[C], R_0S[C], R_SS[C2], R_P0[C], R_0P[C], R_SP[C2], R_PS[C2], R_PP[C2], T_00, T_S0[C], T_0S[C], T_SS[C2], T_P0[C], T_0P[C], T_SP[C2], T_PS[C2], T_PP[C2], R_00I[Y], R_S0I[CY], R_0SI[CY], R_SSI[C2Y], R_P0I[CY], R_0PI[CY], R_SPI[C2Y], R_PSI[C2Y], R_PPI[C2Y], T_00I[Y], T_S0I[CY], T_0SI[CY], T_SSI[C2Y], T_P0I[CY], T_0PI[CY], T_SPI[C2Y], T_PSI[C2Y], T_PPI[C2Y]
            end
            steps:begin
                [R_00, S] ⇌ [R_S0]
                [R_00, S] ⇌ [R_0S]
                [R_0S, S] ⇌ [R_SS]
                [R_S0, S] ⇌ [R_SS]
                [R_0P, S] ⇌ [R_SP]
                [R_P0, S] ⇌ [R_PS]
                [R_S0] <--> [R_P0]
                [R_0S] <--> [R_0P]
                [R_SS] <--> [R_SP]
                [R_SS] <--> [R_PS]
                [R_SP] <--> [R_PP]
                [R_PS] <--> [R_PP]
                [R_00, P] ⇌ [R_P0]
                [R_00, P] ⇌ [R_0P]
                [R_P0, P] ⇌ [R_PP]
                [R_0P, P] ⇌ [R_PP]
                [R_S0, P] ⇌ [R_SP]
                [R_0S, P] ⇌ [R_PS]
                [T_00, S] ⇌ [T_S0]
                [T_00, S] ⇌ [T_0S]
                [T_S0, S] ⇌ [T_SS]
                [T_0S, S] ⇌ [T_SS]
                [T_0P, S] ⇌ [T_SP]
                [T_P0, S] ⇌ [T_PS]
                [T_S0] <--> [T_P0]
                [T_0S] <--> [T_0P]
                [T_SS] <--> [T_SP]
                [T_SS] <--> [T_PS]
                [T_SP] <--> [T_PP]
                [T_PS] <--> [T_PP]
                [T_00, P] ⇌ [T_P0]
                [T_00, P] ⇌ [T_0P]
                [T_P0, P] ⇌ [T_PP]
                [T_0P, P] ⇌ [T_PP]
                [T_S0, P] ⇌ [T_SP]
                [T_0S, P] ⇌ [T_PS]
                [R_00] ⇌ [T_00]
                [R_00, I] ⇌ [R_00I]
                [R_S0, I] ⇌ [R_S0I]
                [R_0S, I] ⇌ [R_0SI]
                [R_SS, I] ⇌ [R_SSI]
                [R_P0, I] ⇌ [R_P0I]
                [R_0P, I] ⇌ [R_0PI]
                [R_SP, I] ⇌ [R_SPI]
                [R_PS, I] ⇌ [R_PSI]
                [R_PP, I] ⇌ [R_PPI]
                [T_00, I] ⇌ [T_00I]
                [T_S0, I] ⇌ [T_S0I]
                [T_0S, I] ⇌ [T_0SI]
                [T_SS, I] ⇌ [T_SSI]
                [T_P0, I] ⇌ [T_P0I]
                [T_0P, I] ⇌ [T_0PI]
                [T_SP, I] ⇌ [T_SPI]
                [T_PS, I] ⇌ [T_PSI]
                [T_PP, I] ⇌ [T_PPI]
                [R_00I, S] ⇌ [R_S0I]
                [R_00I, S] ⇌ [R_0SI]
                [R_0SI, S] ⇌ [R_SSI]
                [R_S0I, S] ⇌ [R_SSI]
                [R_0PI, S] ⇌ [R_SPI]
                [R_P0I, S] ⇌ [R_PSI]
                [R_00I, P] ⇌ [R_P0I]
                [R_00I, P] ⇌ [R_0PI]
                [R_P0I, P] ⇌ [R_PPI]
                [R_0PI, P] ⇌ [R_PPI]
                [R_S0I, P] ⇌ [R_SPI]
                [R_0SI, P] ⇌ [R_PSI]
                [T_00I, S] ⇌ [T_S0I]
                [T_00I, S] ⇌ [T_0SI]
                [T_0SI, S] ⇌ [T_SSI]
                [T_S0I, S] ⇌ [T_SSI]
                [T_0PI, S] ⇌ [T_SPI]
                [T_P0I, S] ⇌ [T_PSI]
                [T_00I, P] ⇌ [T_P0I]
                [T_00I, P] ⇌ [T_0PI]
                [T_P0I, P] ⇌ [T_PPI]
                [T_0PI, P] ⇌ [T_PPI]
                [T_S0I, P] ⇌ [T_SPI]
                [T_0SI, P] ⇌ [T_PSI]
                [R_00I] ⇌ [T_00I]
            end
            constraints:begin
                K2 = K1
                K3 = K1
                K4 = K1
                K5 = K1
                K6 = K1
                k8f = k7f
                k9f = k7f
                k10f = k7f
                k11f = k7f
                k12f = k7f
                K14 = K13
                K15 = K13
                K16 = K13
                K17 = K13
                K18 = K13
                K20 = K19
                K21 = K19
                K22 = K19
                K23 = K19
                K24 = K19
                k26f = k25f
                k27f = k25f
                k28f = k25f
                k29f = k25f
                k30f = k25f
                K32 = K31
                K33 = K31
                K34 = K31
                K35 = K31
                K36 = K31
                K39 = K38
                K40 = K38
                K41 = K38
                K42 = K38
                K43 = K38
                K44 = K38
                K45 = K38
                K46 = K38
                K48 = K47
                K49 = K47
                K50 = K47
                K51 = K47
                K52 = K47
                K53 = K47
                K54 = K47
                K55 = K47
                K56 = K1
                K57 = K1
                K58 = K1
                K59 = K1
                K60 = K1
                K61 = K1
                K62 = K13
                K63 = K13
                K64 = K13
                K65 = K13
                K66 = K13
                K67 = K13
                K68 = K19
                K69 = K19
                K70 = K19
                K71 = K19
                K72 = K19
                K73 = K19
                K74 = K31
                K75 = K31
                K76 = K31
                K77 = K31
                K78 = K31
                K79 = K31
            end
        end

        # K37 is Wegscheider-dependent: K37 = K47*K80/K38
        # k7r, k25r are Haldane-dependent
        function rate_mwc_dimer_inh(p, c)
            (; K1, k7f, k7r, K13, K19, k25f, k25r, K31, K37,
               K38, K47, Et) = p
            (; S, P, I) = c
            r_flux = k7f * S / K1 - k7r * P / K13
            t_flux = k25f * S / K19 - k25r * P / K31
            r_factor = 1.0 + S / K1 + P / K13
            t_factor = 1.0 + S / K19 + P / K31
            num = 2.0 * (r_flux * r_factor +
                         K37 * t_flux * t_factor)
            denom = r_factor^2 * (1.0 + I / K38) +
                    K37 * t_factor^2 * (1.0 + I / K47)
            return Et * num / denom
        end

        push!(specs, MechanismTestSpec(
            name="MWC Dimer + Independent Inhibitor",
            mechanism=m,
            metabolite_names=[:S, :P, :I],
            expected_n_states=36,
            expected_n_steps=80,
            expected_n_metabolites=3,
            expected_n_haldane=12,
            expected_n_wegscheider=1,
            expected_n_independent_params=9,
            expected_identifiability_deficit=-6,
            expected_is_identifiable=true,
            analytical_rate_fn=(p, c) ->
                rate_mwc_dimer_inh(merge(p, (Et=p.Et,)), c),
            expected_factored_num=
            "2 * (1 + S / K1 + P / K13) * (k7f * S / K1 - k7r * P / K13) + 2 * K37 * (1 + S / K19 + P / K31) * (k25f * S / K19 - k25r * P / K31)",
            expected_factored_denom=
            "(1 + S / K1 + P / K13) ^ 2 * (1 + I / K38) + K37 * (1 + S / K19 + P / K31) ^ 2 * (1 + I / K47)",
        ))
    end

    # 26B. OligomericEnzymeMechanism equivalent of MWC Dimer + Independent Inhibitor (spec #26)
    #      The Wegscheider constraint K80 = K47*K37/K38 (R_00I ⇌ T_00I equilibrium)
    #      is automatically satisfied by the conformational assembly formula — no
    #      explicit constraint needed in the OligomericEnzymeMechanism DSL.
    if isdefined(EnzymeRates, :OligomericEnzymeMechanism)
        let
            m = @enzyme_mechanism begin
                metabolites: S[C], P[C], I[Y]
                conformations: 2    # R and T
                site(:catalytic, 2):begin
                    species: begin
                        substrates: S[C]
                        products: P[C]
                        enzymes: E_c, E_S[C], E_P[C]
                    end
                    steps: begin
                        [E_c, S] ⇌ [E_S]    # K1 (R), K1_T (T)
                        [E_c, P] ⇌ [E_P]    # K2 (R), K2_T (T)
                        [E_S] <--> [E_P]     # k3f, k3r Haldane (R); k3f_T, k3r_T Haldane (T)
                    end
                end
                site(:regulatory, 1):begin  # enzyme-level; conformation-specific Ki
                    ligands: I    # K_I_reg1 (R state), K_I_T_reg1 (T state)
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
                num   = 2.0 * (r_flux * r_factor + L * t_flux * t_factor)
                denom = r_factor^2 * (1.0 + I / K_I_reg1) +
                        L * t_factor^2 * (1.0 + I / K_I_T_reg1)
                return Et * num / denom
            end

            push!(specs, MechanismTestSpec(
                name="MWC Dimer + Independent Inhibitor [OligomericEnzymeMechanism]",
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
                "2 * ((k3f * S / K1 - k3r * P / K2) * (1 + S / K1 + P / K2)" *
                " + L * (k3f_T * S / K1_T - k3r_T * P / K2_T) * (1 + S / K1_T + P / K2_T))",
                expected_factored_denom=
                "(1 + S / K1 + P / K2) ^ 2 * (1 + I / K_I_reg1)" *
                " + L * (1 + S / K1_T + P / K2_T) ^ 2 * (1 + I / K_I_T_reg1)",
            ))
        end
    end

    # 27. Two Competitive Inhibitors (monomer)
    #     Tests: multiple additive dead-end terms (flat sum denominator).
    #     I1 and I2 both bind only free E (competitive with S/P and each other).
    #     Denom: 1 + S/K1 + P/K3 + I1/K4 + I2/K5
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:I1, I2
                enzymes:E, E_S[C], E_P[C], E_I1[X], E_I2[Y]
            end
            steps:begin
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
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:I1, I2
                enzymes:E, E_S[C], E_P[C],
                E_I1[X], E_S_I1[CX], E_P_I1[CX],
                E_I2[Y], E_S_I2[CY], E_P_I2[CY],
                E_I1_I2[XY], E_S_I1_I2[CXY], E_P_I1_I2[CXY]
            end
            steps:begin
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
                [E, P] ⇌ [E_P]
                [E, I1] ⇌ [E_I1]
                [E_S, I1] ⇌ [E_S_I1]
                [E_P, I1] ⇌ [E_P_I1]
                [E, I2] ⇌ [E_I2]
                [E_S, I2] ⇌ [E_S_I2]
                [E_P, I2] ⇌ [E_P_I2]
                [E_I1, I2] ⇌ [E_I1_I2]
                [E_S_I1, I2] ⇌ [E_S_I1_I2]
                [E_P_I1, I2] ⇌ [E_P_I1_I2]
                [E_I1, S] ⇌ [E_S_I1]
                [E_I1, P] ⇌ [E_P_I1]
                [E_I2, S] ⇌ [E_S_I2]
                [E_I2, P] ⇌ [E_P_I2]
                [E_I1_I2, S] ⇌ [E_S_I1_I2]
                [E_I1_I2, P] ⇌ [E_P_I1_I2]
                [E_I2, I1] ⇌ [E_I1_I2]
                [E_S_I2, I1] ⇌ [E_S_I1_I2]
                [E_P_I2, I1] ⇌ [E_P_I1_I2]
            end
            constraints:begin
                K5 = K4
                K6 = K4
                K19 = K4
                K20 = K4
                K21 = K4
                K8 = K7
                K9 = K7
                K10 = K7
                K11 = K7
                K12 = K7
                K13 = K1
                K15 = K1
                K17 = K1
                K14 = K3
                K16 = K3
                K18 = K3
            end
        end

        function rate_two_noncomp_inh(p, c)
            (; K1, k2f, k2r, K3, K4, K7, Et) = p
            (; S, P, I1, I2) = c
            num = k2f * S / K1 - k2r * P / K3
            denom = (1.0 + S / K1 + P / K3) *
                    (1.0 + I1 / K4) * (1.0 + I2 / K7)
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
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "(1 + S / K1 + P / K3) * (1 + I1 / K4) * (1 + I2 / K7)",
        ))
    end

    # 29. Non-competitive + Competitive Inhibitor (monomer)
    #     Tests: multiplicative product + additive dead-end term.
    #     I1 binds all forms independently (non-competitive).
    #     I2 binds only free E (competitive).
    #     Denom: (1+S/K1+P/K3)*(1+I1/K4) + I2/K9
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:I1, I2
                enzymes:E, E_S[C], E_P[C],
                E_I1[X], E_S_I1[CX], E_P_I1[CX],
                E_I2[Y]
            end
            steps:begin
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
                [E, P] ⇌ [E_P]
                [E, I1] ⇌ [E_I1]
                [E_S, I1] ⇌ [E_S_I1]
                [E_P, I1] ⇌ [E_P_I1]
                [E_I1, S] ⇌ [E_S_I1]
                [E_I1, P] ⇌ [E_P_I1]
                [E, I2] ⇌ [E_I2]
            end
            constraints:begin
                K5 = K4
                K6 = K4
                K7 = K1
                K8 = K3
            end
        end

        function rate_noncomp_comp_inh(p, c)
            (; K1, k2f, k2r, K3, K4, K9, Et) = p
            (; S, P, I1, I2) = c
            num = k2f * S / K1 - k2r * P / K3
            denom = (1.0 + S / K1 + P / K3) *
                    (1.0 + I1 / K4) + I2 / K9
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
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "I2 / K9 + (1 + S / K1 + P / K3) * (1 + I1 / K4)",
        ))
    end

    # 30. Uncompetitive + Competitive Inhibitor (monomer)
    #     Tests: mixed additive structure with nested multiplicative term.
    #     I1 binds only ES (uncompetitive). I2 binds only free E (competitive).
    #     Denom: 1 + I2/K5 + P/K3 + (S/K1)*(1+I1/K4)
    let
        m = @enzyme_mechanism begin
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:I1, I2
                enzymes:E, E_S[C], E_P[C],
                E_S_I1[CX],
                E_I2[Y]
            end
            steps:begin
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
            species:begin
                substrates:S[C]
                products:P[C]
                regulators:I1, I2
                enzymes:E, E_S[C], E_P[C],
                E_I1[X], E_S_I1[CX], E_P_I1[CX],
                E_I2[X], E_S_I2[CX], E_P_I2[CX]
            end
            steps:begin
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
                [E, P] ⇌ [E_P]
                [E, I1] ⇌ [E_I1]
                [E_S, I1] ⇌ [E_S_I1]
                [E_P, I1] ⇌ [E_P_I1]
                [E_I1, S] ⇌ [E_S_I1]
                [E_I1, P] ⇌ [E_P_I1]
                [E, I2] ⇌ [E_I2]
                [E_S, I2] ⇌ [E_S_I2]
                [E_P, I2] ⇌ [E_P_I2]
                [E_I2, S] ⇌ [E_S_I2]
                [E_I2, P] ⇌ [E_P_I2]
            end
            constraints:begin
                K5 = K4
                K6 = K4
                K7 = K1
                K8 = K3
                K10 = K9
                K11 = K9
                K12 = K1
                K13 = K3
            end
        end

        function rate_two_samesite_inh(p, c)
            (; K1, k2f, k2r, K3, K4, K9, Et) = p
            (; S, P, I1, I2) = c
            num = k2f * S / K1 - k2r * P / K3
            denom = (1.0 + S / K1 + P / K3) *
                    (1.0 + I1 / K4 + I2 / K9)
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
            "k2f * S / K1 - k2r * P / K3",
            expected_factored_denom=
            "(1 + I1 / K4 + I2 / K9) * (1 + S / K1 + P / K3)",
        ))
    end

    # 32. MWC Tetramer — Random-Order Bi-Bi RE + Two Allosteric Sites
    #
    # Reaction: S1 + S2 ⇌ P1 + P2  (Bi-Bi)
    #
    # Catalytic site (×4 subunits, NConf=2):
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
    if isdefined(EnzymeRates, :OligomericEnzymeMechanism)
        let
            m = @enzyme_mechanism begin
                metabolites: S1[C], S2[N], P1[C], P2[N], R1[X], R2[Y], R3[Z]
                conformations: 2    # NConf=2: R (active) and T (tense) states
                site(:catalytic, 4):begin
                    species: begin
                        substrates: S1[C], S2[N]
                        products: P1[C], P2[N]
                        enzymes: E_c,
                                 E_S1[C], E_P1[C],
                                 E_S2[N], E_P2[N],
                                 E_S1S2[CN], E_P1P2[CN],
                                 E_S1P2[CN], E_P1S2[CN]
                    end
                    steps: begin
                        [E_c,  S1] ⇌ [E_S1]         # K1
                        [E_c,  P1] ⇌ [E_P1]         # K2
                        [E_c,  S2] ⇌ [E_S2]         # K3
                        [E_c,  P2] ⇌ [E_P2]         # K4
                        [E_S1, S2] ⇌ [E_S1S2]       # K5  = K3
                        [E_S2, S1] ⇌ [E_S1S2]       # K6  = K1
                        [E_S1, P2] ⇌ [E_S1P2]       # K7  = K4
                        [E_P2, S1] ⇌ [E_S1P2]       # K8  = K1
                        [E_P1, S2] ⇌ [E_P1S2]       # K9  = K3
                        [E_S2, P1] ⇌ [E_P1S2]       # K10 = K2
                        [E_P1, P2] ⇌ [E_P1P2]       # K11 = K4
                        [E_P2, P1] ⇌ [E_P1P2]       # K12 = K2
                        [E_S1S2] <--> [E_P1P2]       # k13f, k13r (k13r via Haldane)
                    end
                    constraints: begin
                        K5=K3; K6=K1; K7=K4; K8=K1; K9=K3; K10=K2; K11=K4; K12=K2
                    end
                end
                site(:regulatory, 4):begin
                    ligands: R1, R2
                end
                site(:regulatory, 4):begin
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
            function rate_mwc_tetramer_bi_bi(params, concs)
                (; K1, K2, K3, K4, k13f, k13r,
                   K1_T, K2_T, K3_T, K4_T, k13f_T, k13r_T,
                   K_R1_reg1, K_R2_reg1, K_R1_T_reg1, K_R2_T_reg1,
                   K_R3_reg2, K_R3_T_reg2,
                   L, Et) = params
                (; S1, S2, P1, P2, R1, R2, R3) = concs

                # R-state: catalytic site factors by site-A ⊗ site-B independence
                Q_A_R   = 1.0 + S1 / K1 + P1 / K2
                Q_B_R   = 1.0 + S2 / K3 + P2 / K4
                Q_cat_R = Q_A_R * Q_B_R
                N_cat_R = k13f * S1 * S2 / (K1 * K3) - k13r * P1 * P2 / (K2 * K4)

                # R-state: regulatory site partition functions (star topology → direct sum)
                Q_reg1_R = 1.0 + R1 / K_R1_reg1 + R2 / K_R2_reg1
                Q_reg2_R = 1.0 + R3 / K_R3_reg2

                # T-state: same structure with _T parameters
                Q_A_T   = 1.0 + S1 / K1_T + P1 / K2_T
                Q_B_T   = 1.0 + S2 / K3_T + P2 / K4_T
                Q_cat_T = Q_A_T * Q_B_T
                N_cat_T = k13f_T * S1 * S2 / (K1_T * K3_T) - k13r_T * P1 * P2 / (K2_T * K4_T)

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
                expected_factored_num=
                "4 * ((k13f * S1 * S2 / (K1 * K3) - k13r * P1 * P2 / (K2 * K4))" *
                " * ((1 + S1 / K1 + P1 / K2) * (1 + S2 / K3 + P2 / K4)) ^ 3" *
                " * (1 + R1 / K_R1_reg1 + R2 / K_R2_reg1) ^ 4" *
                " * (1 + R3 / K_R3_reg2) ^ 4" *
                " + L * (k13f_T * S1 * S2 / (K1_T * K3_T) - k13r_T * P1 * P2 / (K2_T * K4_T))" *
                " * ((1 + S1 / K1_T + P1 / K2_T) * (1 + S2 / K3_T + P2 / K4_T)) ^ 3" *
                " * (1 + R1 / K_R1_T_reg1 + R2 / K_R2_T_reg1) ^ 4" *
                " * (1 + R3 / K_R3_T_reg2) ^ 4)",
                expected_factored_denom=
                "((1 + S1 / K1 + P1 / K2) * (1 + S2 / K3 + P2 / K4)) ^ 4" *
                " * (1 + R1 / K_R1_reg1 + R2 / K_R2_reg1) ^ 4" *
                " * (1 + R3 / K_R3_reg2) ^ 4" *
                " + L * ((1 + S1 / K1_T + P1 / K2_T) * (1 + S2 / K3_T + P2 / K4_T)) ^ 4" *
                " * (1 + R1 / K_R1_T_reg1 + R2 / K_R2_T_reg1) ^ 4" *
                " * (1 + R3 / K_R3_T_reg2) ^ 4",
            ))
        end
    end

    return specs
end

const MECHANISM_TEST_SPECS = build_mechanism_test_specs()
