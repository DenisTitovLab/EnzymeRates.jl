@testset "Identifiability Analysis" begin

    @testset "Mode types" begin
        # Test that mode types are exported and usable
        @test Raw isa RawMode
        @test HaldaneWegscheider isa HaldaneWegscheiderMode
        @test IdentifiableHaldaneWegscheider isa IdentifiableHaldaneWegscheiderMode
    end

    @testset "Monomial extraction" begin
        # Test internal monomial extraction on Uni-Uni
        m, _ = make_uni_uni()
        M = typeof(m)
        A, k_params, monomials = EnzymeRates._monomial_exponent_matrix(M)

        @test length(k_params) > 0
        @test length(monomials) > 0
        @test size(A, 1) == length(monomials)
        @test size(A, 2) == length(k_params)

        # All exponents should be non-negative
        @test all(A .>= 0)
    end

    @testset "Uni-Uni identifiability" begin
        m, _ = make_uni_uni()

        # Uni-Uni should be fully identifiable
        @test is_identifiable(m)

        # No non-identifiable directions
        NS = non_identifiable_directions(m)
        @test size(NS, 2) == 0

        # Check combinations list
        combos = identifiable_combinations(m)
        indep = independent_parameters(m)
        @test length(combos) == length(indep)
    end

    @testset "Seq Bi-Bi identifiability" begin
        m, _ = make_seq_bibi()

        # Sequential mechanisms should be fully identifiable
        @test is_identifiable(m)

        NS = non_identifiable_directions(m)
        @test size(NS, 2) == 0
    end

    @testset "Ping-Pong Bi-Bi identifiability" begin
        m, _ = make_pingpong_bibi()

        # Check that analysis completes without error
        combos = identifiable_combinations(m)
        NS = non_identifiable_directions(m)

        # Number of identifiable combinations + null space dimension
        # should equal number of independent parameters
        indep = independent_parameters(m)
        @test length(combos) + size(NS, 2) == length(indep)
    end

    @testset "Random Bi-Bi identifiability" begin
        m, _ = make_random_bibi()

        # Check that analysis completes
        combos = identifiable_combinations(m)
        NS = non_identifiable_directions(m)

        indep = independent_parameters(m)
        @test length(combos) + size(NS, 2) == length(indep)

        # Verify that rate_equation_string generation works
        s = rate_equation_string(m)
        @test occursin("v = E_total", s)
    end

    @testset "Mode-dispatched parameters" begin
        m, _ = make_uni_uni()

        # Raw mode: all k's + E_total
        raw_params = parameters(m, Raw)
        @test :k1f in raw_params
        @test :k1r in raw_params
        @test :k2f in raw_params
        @test :k2r in raw_params
        @test :E_total in raw_params
        @test :Keq ∉ raw_params

        # HaldaneWegscheider mode: independent k's + Keq + E_total
        hw_params = parameters(m, HaldaneWegscheider)
        @test :Keq in hw_params
        @test :E_total in hw_params
        # HW has Keq instead of some k's - total count may be same but different parameters
        @test Set(hw_params) != Set(raw_params)  # Different parameter sets
        # At least one raw k-parameter should be missing from HW (dependent)
        raw_ks = filter(p -> startswith(string(p), "k"), raw_params)
        hw_ks = filter(p -> startswith(string(p), "k"), hw_params)
        @test length(hw_ks) < length(raw_ks)  # Some k's are dependent

        # IdentifiableHaldaneWegscheider mode (default): identifiable + Keq + E_total
        ident_params = parameters(m)
        @test :Keq in ident_params
        @test :E_total in ident_params

        # For fully identifiable mechanism, HW and Identifiable should have same k's
        @test Set(hw_params) == Set(ident_params)
    end

    @testset "Mode-dispatched rate_equation" begin
        m, met_names = make_uni_uni()
        rng = Random.MersenneTwister(12345)

        # Generate parameters for each mode
        # Raw mode needs all k's
        raw_params = (k1f=1.0, k1r=0.5, k2f=2.0, k2r=0.3, E_total=0.1)
        concs = (S=1.0, P=0.5)

        # HaldaneWegscheider needs independent k's + Keq
        Keq = (raw_params.k1f * raw_params.k2f) / (raw_params.k1r * raw_params.k2r)
        hw_params = (k1f=1.0, k2f=2.0, k2r=0.3, Keq=Keq, E_total=0.1)

        # Identifiable (same as HW for fully identifiable mechanism)
        ident_params = hw_params

        # All modes should give the same rate
        v_raw = rate_equation(m, raw_params, concs, Raw)
        v_hw = rate_equation(m, hw_params, concs, HaldaneWegscheider)
        v_ident = rate_equation(m, ident_params, concs, IdentifiableHaldaneWegscheider)
        v_default = rate_equation(m, ident_params, concs)  # default is Identifiable

        @test v_raw ≈ v_hw rtol=1e-10
        @test v_hw ≈ v_ident rtol=1e-10
        @test v_ident ≈ v_default rtol=1e-10
    end

    @testset "Mode-dispatched rate_equation_string" begin
        m, _ = make_uni_uni()

        # Raw mode: no constraints shown
        raw_str = rate_equation_string(m, Raw)
        @test occursin("v = E_total", raw_str)
        @test !occursin("Keq", raw_str)  # Raw doesn't use Keq

        # HaldaneWegscheider mode: shows Haldane constraints
        hw_str = rate_equation_string(m, HaldaneWegscheider)
        @test occursin("v = E_total", hw_str)
        @test occursin("k1r", hw_str)  # Shows dependent param expression

        # IdentifiableHaldaneWegscheider mode (default)
        ident_str = rate_equation_string(m)
        @test occursin("v = E_total", ident_str)
    end

    @testset "Scaling invariance verification" begin
        # For each non-identifiable direction, verify that scaling
        # parameters in that direction doesn't change the rate

        for (make_fn, name) in [
            (make_uni_uni, "Uni-Uni"),
            (make_seq_bibi, "Seq Bi-Bi"),
            (make_random_bibi, "Random Bi-Bi"),
        ]
            @testset "$name scaling invariance" begin
                m, met_names = make_fn()
                NS = non_identifiable_directions(m)
                indep = collect(independent_parameters(m))

                rng = Random.MersenneTwister(12345)

                # Generate random parameters and concentrations
                new_params, concs, _ = random_independent_params_concs(m, met_names; rng=rng)

                # Compute baseline rate using HaldaneWegscheider mode
                v_base = rate_equation(m, new_params, concs, HaldaneWegscheider)

                # For each null space direction, scale and check invariance
                for j in 1:size(NS, 2)
                    # Use a random scaling factor
                    c = 0.5 + rand(rng)  # Random c between 0.5 and 1.5

                    # Build scaled parameters
                    scaled_vals = Float64[]
                    for (i, k) in enumerate(indep)
                        base_val = new_params[k]
                        exp = NS[i, j]
                        push!(scaled_vals, base_val * c^exp)
                    end
                    # Add Keq and E_total unchanged
                    scaled_keys = (indep..., :Keq, :E_total)
                    scaled_vals_full = (scaled_vals..., new_params.Keq, new_params.E_total)
                    scaled_params = NamedTuple{scaled_keys}(scaled_vals_full)

                    v_scaled = rate_equation(m, scaled_params, concs, HaldaneWegscheider)

                    # Rate should be invariant (within numerical tolerance)
                    @test v_base ≈ v_scaled rtol=1e-10
                end
            end
        end
    end

    @testset "Integer null space consistency" begin
        # Verify that _integer_nullspace used in identifiability
        # produces valid null vectors

        for (make_fn, _) in [
            (make_uni_uni, "Uni-Uni"),
            (make_seq_bibi, "Seq Bi-Bi"),
            (make_random_bibi, "Random Bi-Bi"),
        ]
            m, _ = make_fn()
            M = typeof(m)
            A, _, _ = EnzymeRates._monomial_exponent_matrix(M)

            n_monomials, n_params = size(A)
            if n_monomials > 1
                # Build difference matrix
                D = zeros(Int, n_monomials - 1, n_params)
                for i in 1:(n_monomials - 1)
                    D[i, :] = A[i + 1, :] - A[1, :]
                end

                NS = EnzymeRates._integer_nullspace(D)

                # Verify D * NS = 0
                for j in 1:size(NS, 2)
                    @test all(D * NS[:, j] .== 0)
                end
            end
        end
    end

    @testset "Sequential mechanisms are identifiable" begin
        # All simple sequential mechanisms should be fully identifiable
        for (make_fn, name) in [
            (make_uni_uni, "Uni-Uni"),
            (make_seq_unibi, "Seq Uni-Bi"),
            (make_seq_biuni, "Seq Bi-Uni"),
            (make_seq_bibi, "Seq Bi-Bi"),
            (make_seq_biter, "Seq Bi-Ter"),
            (make_seq_terbi, "Seq Ter-Bi"),
            (make_seq_terter, "Seq Ter-Ter"),
        ]
            @testset "$name" begin
                m, _ = make_fn()
                @test is_identifiable(m)
                @test size(non_identifiable_directions(m), 2) == 0
            end
        end
    end

    @testset "Default mode is IdentifiableHaldaneWegscheider" begin
        m, met_names = make_uni_uni()
        rng = Random.MersenneTwister(54321)
        params, concs, _ = random_independent_params_concs(m, met_names; rng=rng)

        # Default should equal explicit IdentifiableHaldaneWegscheider
        v_default = rate_equation(m, params, concs)
        v_explicit = rate_equation(m, params, concs, IdentifiableHaldaneWegscheider)
        @test v_default == v_explicit

        # Same for parameters
        @test parameters(m) == parameters(m, IdentifiableHaldaneWegscheider)

        # Same for rate_equation_string
        @test rate_equation_string(m) == rate_equation_string(m, IdentifiableHaldaneWegscheider)
    end

end
