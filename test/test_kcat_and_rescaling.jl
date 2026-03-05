# Tests for kcat computation and parameter rescaling
# Validates _is_ss_rate_constant, _kcat_forward, and rescale_parameter_values

@testset "kcat and rescaling" begin

# ── 1a. _is_ss_rate_constant classification ──────────────────
@testset "_is_ss_rate_constant" begin
    # SS rate constants: lowercase k followed by digit
    for sym in (:k1f, :k2r, :k3f_T, :k10f)
        @test EnzymeRates._is_ss_rate_constant(sym)
    end
    # Non-SS parameters
    for sym in (:K1, :K2, :K_I_reg1, :Keq, :L, :E_total)
        @test !EnzymeRates._is_ss_rate_constant(sym)
    end
end

# ── Helper: generate random reduced params for a mechanism ───
function random_reduced_params(m; rng=Random.default_rng())
    fp = fitted_params(m)
    vals = Tuple(0.1 + 9.9 * rand(rng) for _ in fp)
    Keq_val = 0.1 + 9.9 * rand(rng)
    E_total_val = 0.1 + 9.9 * rand(rng)
    keys_out = (fp..., :Keq, :E_total)
    vals_out = (vals..., Keq_val, E_total_val)
    NamedTuple{keys_out}(vals_out)
end

# ── 1b. Analytical kcat tests per-mechanism ──────────────────
@testset "analytical kcat" begin
    specs = build_mechanism_test_specs()

    # Find mechanisms by name
    find_spec(name) = first(s for s in specs if s.name == name)

    # RE Uni-Uni (1 RE + 1 SS): kcat = k2f
    let spec = find_spec("RE Uni-Uni")
        m = spec.mechanism
        rng = Random.MersenneTwister(42)
        params = random_reduced_params(m; rng)
        kcat = EnzymeRates._kcat_forward(m, params)
        @test kcat ≈ params.k2f
    end

    # Segel Ordered Bi Bi (4 SS): kcat = k3f*k4f/(k3f+k4f)
    let spec = find_spec("Segel Ordered Bi Bi")
        m = spec.mechanism
        rng = Random.MersenneTwister(43)
        params = random_reduced_params(m; rng)
        kcat = EnzymeRates._kcat_forward(m, params)
        @test kcat ≈ params.k3f * params.k4f / (params.k3f + params.k4f)
    end

    # Segel Theorell-Chance Bi Bi (3 SS): kcat = k3f
    let spec = find_spec("Segel Theorell-Chance Bi Bi")
        m = spec.mechanism
        rng = Random.MersenneTwister(44)
        params = random_reduced_params(m; rng)
        kcat = EnzymeRates._kcat_forward(m, params)
        @test kcat ≈ params.k3f rtol=1e-10
    end

    # Segel Ping Pong Bi Bi (4 SS): kcat = k2f*k4f/(k2f+k4f)
    let spec = find_spec("Segel Ping Pong Bi Bi")
        m = spec.mechanism
        rng = Random.MersenneTwister(45)
        params = random_reduced_params(m; rng)
        kcat = EnzymeRates._kcat_forward(m, params)
        @test kcat ≈ params.k2f * params.k4f / (params.k2f + params.k4f)
    end

    # Competitive Inhibitor: kcat = k2f
    let spec = find_spec("Competitive Inhibitor")
        m = spec.mechanism
        rng = Random.MersenneTwister(46)
        params = random_reduced_params(m; rng)
        kcat = EnzymeRates._kcat_forward(m, params)
        @test kcat ≈ params.k2f
    end

    # Essential Activator: kcat = k2f
    let spec = find_spec("Essential Activator")
        m = spec.mechanism
        rng = Random.MersenneTwister(47)
        params = random_reduced_params(m; rng)
        kcat = EnzymeRates._kcat_forward(m, params)
        @test kcat ≈ params.k2f
    end

    # Non-essential Activator: kcat = max(k2f, k5f)
    let spec = find_spec("Non-essential Activator")
        m = spec.mechanism
        # Case 1: k5f > k2f → kcat = k5f
        rng = Random.MersenneTwister(48)
        params = random_reduced_params(m; rng)
        # Force k5f > k2f
        params_hi = merge(params, (k5f=10.0, k2f=1.0))
        kcat = EnzymeRates._kcat_forward(m, params_hi)
        @test kcat ≈ 10.0

        # Case 2: k5f < k2f → kcat = k2f
        params_lo = merge(params, (k5f=1.0, k2f=10.0))
        kcat = EnzymeRates._kcat_forward(m, params_lo)
        @test kcat ≈ 10.0
    end
end

# ── 1c. Comprehensive V=1 test over all mechanism specs ──────
@testset "V≈1 after rescaling: $(spec.name)" for spec in
        build_mechanism_test_specs()
    m = spec.mechanism
    rng = Random.MersenneTwister(100)
    params = random_reduced_params(m; rng)

    # kcat should be positive
    kcat_orig = EnzymeRates._kcat_forward(m, params)
    @test kcat_orig > 0

    # Rescale so kcat = 1
    norm = rescale_parameter_values(m, params)
    kcat_norm = EnzymeRates._kcat_forward(m, norm)
    @test kcat_norm ≈ 1.0 rtol=1e-10

    # K values (non-SS params) unchanged
    for k in keys(params)
        if !EnzymeRates._is_ss_rate_constant(k)
            @test norm[k] == params[k]
        end
    end

    # V ≈ 1 at saturating substrates, products=0
    # Build concentration tuple: substrates=BIG, products=0,
    # regulators: try all 2^R corners
    met_names = metabolites(m)
    sub_names = Symbol[s[1] for s in EnzymeRates.substrates(m)]
    prod_names = Symbol[p[1] for p in EnzymeRates.products(m)]
    reg_names = Symbol[]
    for r in EnzymeRates.regulators(m)
        push!(reg_names, r isa Tuple ? r[1] : r)
    end
    n_reg = length(reg_names)

    # Set E_total = 1 for the V=1 check
    norm_e1 = merge(norm, (E_total=1.0,))

    BIG = 1e6
    max_rate = 0.0
    for mask in 0:(2^n_reg - 1)
        conc_dict = Dict{Symbol,Float64}()
        for s in sub_names
            conc_dict[s] = BIG
        end
        for p in prod_names
            conc_dict[p] = 0.0
        end
        for (i, r) in enumerate(reg_names)
            conc_dict[r] = ((mask >> (i - 1)) & 1) == 1 ? BIG : 0.0
        end
        concs = NamedTuple{Tuple(met_names)}(
            Tuple(conc_dict[n] for n in met_names))
        v = rate_equation(m, concs, norm_e1)
        max_rate = max(max_rate, v)
    end
    @test max_rate ≈ 1.0 rtol=1e-3
end

# ── 1d. Scale invariance ─────────────────────────────────────
@testset "scale invariance" begin
    specs = build_mechanism_test_specs()
    m = specs[1].mechanism  # Uni-Uni
    rng = Random.MersenneTwister(200)
    params = random_reduced_params(m; rng)

    α = 3.7
    scaled_params = NamedTuple{keys(params)}(Tuple(
        EnzymeRates._is_ss_rate_constant(k) ? v * α : v
        for (k, v) in zip(keys(params), values(params))
    ))
    kcat_orig = EnzymeRates._kcat_forward(m, params)
    kcat_scaled = EnzymeRates._kcat_forward(m, scaled_params)
    @test kcat_scaled ≈ α * kcat_orig rtol=1e-10
end

# ── 1e. Rate proportionality ─────────────────────────────────
@testset "rate proportionality" begin
    specs = build_mechanism_test_specs()
    # Test with first 3 small mechanisms
    for spec in specs[1:min(3, end)]
        m = spec.mechanism
        rng = Random.MersenneTwister(300)
        params = random_reduced_params(m; rng)
        norm = rescale_parameter_values(m, params)
        kcat_orig = EnzymeRates._kcat_forward(m, params)

        met_names = metabolites(m)
        conc_vals = Tuple(0.5 + rand(rng) for _ in met_names)
        concs = NamedTuple{Tuple(met_names)}(conc_vals)

        v_orig = rate_equation(m, concs, params)
        v_norm = rate_equation(m, concs, norm)
        @test v_norm / v_orig ≈ 1.0 / kcat_orig rtol=1e-8
    end
end

# ── 1f. Custom kcat target ───────────────────────────────────
@testset "custom kcat target" begin
    specs = build_mechanism_test_specs()
    m = specs[1].mechanism  # Uni-Uni
    rng = Random.MersenneTwister(400)
    params = random_reduced_params(m; rng)

    norm42 = rescale_parameter_values(m, params; kcat=42.0)
    @test EnzymeRates._kcat_forward(m, norm42) ≈ 42.0 rtol=1e-10
end

end  # @testset "kcat and rescaling"
