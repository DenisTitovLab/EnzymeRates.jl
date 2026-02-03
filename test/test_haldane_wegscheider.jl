@testset "Haldane & Wegscheider conditions" begin
    m_uu, mets_uu = make_uni_uni()
    m_bb, mets_bb = make_seq_bibi()
    m_ro, mets_ro = make_random_bibi()
    m_db, mets_db = make_doubly_branched()

    @testset "Constraint counting" begin
        # Uni-Uni: 2 steps, 2 states → 1 constraint (1 Haldane)
        @test length(dependent_parameters(m_uu)) == 1
        @test length(independent_parameters(m_uu)) == 3  # 2*2 - 1 = 3

        # Seq Bi-Bi: 5 steps, 5 states → 1 constraint (1 Haldane)
        @test length(dependent_parameters(m_bb)) == 1
        @test length(independent_parameters(m_bb)) == 9  # 2*5 - 1 = 9

        # Random-order Bi-Bi: 7 steps, 6 states → 2 constraints
        @test length(dependent_parameters(m_ro)) == 2
        @test length(independent_parameters(m_ro)) == 12  # 2*7 - 2 = 12
    end

    @testset "Free-enzyme-binding steps are independent" begin
        # Uni-Uni: k1f (free enzyme binding) and k2r (free enzyme releasing) should be independent
        indep = Set(independent_parameters(m_uu))
        dep = Set(d[1] for d in dependent_parameters(m_uu))
        @test :k1f in indep
        @test :k2r in indep

        # Random-order Bi-Bi: steps 1,2,7 involve free enzyme
        indep_ro = Set(independent_parameters(m_ro))
        @test :k1f in indep_ro
        @test :k1r in indep_ro
        @test :k2f in indep_ro
        @test :k2r in indep_ro
        @test :k7f in indep_ro
        @test :k7r in indep_ro
    end

    @testset "Numerical equivalence with all-params formula" begin
        rng = Random.MersenneTwister(7001)

        # Uni-Uni
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m_uu, mets_uu; rng=rng)
            @test rate_equation(m_uu, new_params, concs) ≈ reference_qssa(m_uu, all_params, concs) rtol=1e-10
        end

        # Random-order Bi-Bi
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m_ro, mets_ro; rng=rng)
            @test rate_equation(m_ro, new_params, concs) ≈ reference_qssa(m_ro, all_params, concs) rtol=1e-10
        end
    end

    @testset "Equilibrium (rate = 0 at Keq)" begin
        Keq = 5.0
        p = (k1f=3.2, k2f=2.5, k2r=1.1, Keq=Keq, E_total=1.0)
        S_eq = 1.0
        P_eq = Keq * S_eq
        v = rate_equation(m_uu, p, (S=S_eq, P=P_eq))
        @test abs(v) < 1e-12
    end

    @testset "Performance" begin
        p = (k1f=3.2, k2f=2.5, k2r=1.1, Keq=5.0, E_total=1.0)
        c = (S=1.0, P=0.5)
        allocs, t = test_rate_equation_performance(m_uu, p, c)
        @test allocs == 0
        @test t < 100e-9
    end

    @testset "Multiple Haldane verification (random-order Bi-Bi)" begin
        rng = Random.MersenneTwister(7002)
        for _ in 1:10
            _, _, all_params = random_independent_params_concs(m_ro, mets_ro; rng=rng)
            # Path 1 (A first): E→EA→EAB→EPQ→EQ→E = steps 1,4,5,6,7
            Keq_path1 = (all_params.k1f * all_params.k4f * all_params.k5f * all_params.k6f * all_params.k7f) /
                        (all_params.k1r * all_params.k4r * all_params.k5r * all_params.k6r * all_params.k7r)
            # Path 2 (B first): E→EB→EAB→EPQ→EQ→E = steps 2,3,5,6,7
            Keq_path2 = (all_params.k2f * all_params.k3f * all_params.k5f * all_params.k6f * all_params.k7f) /
                        (all_params.k2r * all_params.k3r * all_params.k5r * all_params.k6r * all_params.k7r)
            @test Keq_path1 ≈ Keq_path2 rtol=1e-10
        end
    end

    @testset "all_parameters returns all k's" begin
        @test all_parameters(m_uu) == (:k1f, :k1r, :k2f, :k2r)
        @test Set(independent_parameters(m_uu)) ∪ Set(d[1] for d in dependent_parameters(m_uu)) == Set(all_parameters(m_uu))
    end

    @testset "Roundtrip consistency" begin
        mechanisms = [
            (m_uu, mets_uu, false),
            (m_bb, mets_bb, false),
            (m_ro, mets_ro, true),
            (m_db, mets_db, true),
        ]

        for (m, met_names, is_branched) in mechanisms
            rng = Random.MersenneTwister(8001)
            for _ in 1:10
                new_params, _, all_params = random_independent_params_concs(m, met_names; rng=rng)
                if !is_branched
                    ns = n_steps(m)
                    kf_prod = prod(all_params[Symbol("k$(i)f")] for i in 1:ns)
                    kr_prod = prod(all_params[Symbol("k$(i)r")] for i in 1:ns)
                    @test kf_prod / kr_prod ≈ new_params.Keq rtol=1e-10
                else
                    v_new = rate_equation(m, new_params, NamedTuple{Tuple(met_names)}(Tuple(0.1 + 9.9 * rand(rng) for _ in met_names)))
                    @test isfinite(v_new)
                end
            end
        end
    end

    uu_and_ro = [("Uni-Uni", m_uu, mets_uu), ("Random Bi-Bi", m_ro, mets_ro)]

    @testset "Sensitivity / perturbation: $label" for (label, m, met_names) in uu_and_ro
        rng = Random.MersenneTwister(8002)
        new_params, concs, _ = random_independent_params_concs(m, met_names; rng=rng)
        v_base = rate_equation(m, new_params, concs)

        indep = independent_parameters(m)
        for k in indep
            perturbed = merge(new_params, NamedTuple{(k,)}((new_params[k] * 1.01,)))
            v_pert = rate_equation(m, perturbed, concs)
            @test v_pert != v_base
        end

        perturbed_keq = merge(new_params, (Keq=new_params.Keq * 1.01,))
        v_keq = rate_equation(m, perturbed_keq, concs)
        @test v_keq > v_base

        perturbed_et = merge(new_params, (E_total=new_params.E_total * 1.01,))
        v_et = rate_equation(m, perturbed_et, concs)
        @test v_et / v_base ≈ 1.01 rtol=1e-10

        dep = dependent_parameters(m)
        dep_keys = Tuple(d[1] for d in dep)
        dep_vals = Tuple(999.0 for _ in dep)
        params_with_dep = merge(new_params, NamedTuple{dep_keys}(dep_vals))
        @test rate_equation(m, params_with_dep, concs) == v_base
        for (dk, _) in dep
            perturbed_dep = merge(new_params, NamedTuple{(dk,)}((0.001,)))
            @test rate_equation(m, perturbed_dep, concs) == v_base
            perturbed_dep2 = merge(new_params, NamedTuple{(dk,)}((1e8,)))
            @test rate_equation(m, perturbed_dep2, concs) == v_base
        end
    end

    @testset "Extreme parameter values: $label" for (label, m, met_names) in uu_and_ro
        rng = Random.MersenneTwister(8003)

        indep = independent_parameters(m)
        param_keys = (indep..., :Keq, :E_total)
        param_vals = Tuple(exp(rand(rng) * 16 * log(10) - 8 * log(10)) for _ in param_keys)
        new_params = NamedTuple{param_keys}(param_vals)
        conc_vals = Tuple(exp(rand(rng) * 16 * log(10) - 8 * log(10)) for _ in met_names)
        concs = NamedTuple{Tuple(met_names)}(conc_vals)

        all_params = compute_all_params(m, new_params)

        v = rate_equation(m, new_params, concs)
        @test isfinite(v)

        v_ref = reference_qssa(m, all_params, concs)
        @test v ≈ v_ref rtol=1e-4

        Keq = new_params.Keq
        if length(met_names) == 2
            S_eq = 1.0
            P_eq = Keq * S_eq
            eq_concs = NamedTuple{Tuple(met_names)}((S_eq, P_eq))
        else
            A_eq = 1.0; B_eq = 1.0
            P_eq = sqrt(Keq); Q_eq = sqrt(Keq)
            eq_concs = NamedTuple{Tuple(met_names)}((A_eq, B_eq, P_eq, Q_eq))
        end
        v_eq = rate_equation(m, new_params, eq_concs)
        @test abs(v_eq) / (abs(v) + eps()) < 1e-4
    end

    @testset "Doubly-branched mechanism (Wegscheider cycle)" begin
        # 1. Constraint counting
        @test length(dependent_parameters(m_db)) == 2
        @test length(independent_parameters(m_db)) == 8  # 2*5 - 2 = 8

        # 2. Numerical equivalence
        rng = Random.MersenneTwister(8004)
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m_db, mets_db; rng=rng)
            v_new = rate_equation(m_db, new_params, concs)
            v_ref = reference_qssa(m_db, all_params, concs)
            @test v_new ≈ v_ref rtol=1e-10
        end

        # 3. Roundtrip: both catalytic paths give same Keq
        rng = Random.MersenneTwister(8005)
        for _ in 1:10
            _, _, all_params = random_independent_params_concs(m_db, mets_db; rng=rng)
            Keq_path1 = (all_params.k1f * all_params.k2f * all_params.k4f) /
                        (all_params.k1r * all_params.k2r * all_params.k4r)
            Keq_path2 = (all_params.k1f * all_params.k3f * all_params.k5f) /
                        (all_params.k1r * all_params.k3r * all_params.k5r)
            @test Keq_path1 ≈ Keq_path2 rtol=1e-10
        end

        # 4. Equilibrium: rate ≈ 0 at Keq-balanced concentrations
        new_params, _, _ = random_independent_params_concs(m_db, mets_db; rng=rng)
        Keq = new_params.Keq
        S_eq = 1.0; P_eq = Keq * S_eq
        v_eq = rate_equation(m_db, new_params, (S=S_eq, P=P_eq))
        @test abs(v_eq) < 1e-12

        # 5. Performance: zero allocations
        p = (;(k => 1.0 + i * 0.1 for (i, k) in enumerate(independent_parameters(m_db)))..., Keq=5.0, E_total=1.0)
        c = (S=1.0, P=0.5)
        allocs, t = test_rate_equation_performance(m_db, p, c)
        @test allocs == 0
    end
end
