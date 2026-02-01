@testset "Uni-Uni" begin
    function rate_seq_uniuni(p, c)
        S = c.S; P = c.P
        f = [p.k1f * S, p.k2f]
        r = [p.k1r, p.k2r * P]
        return _unicyclic_flux(f, r, p.Et)
    end

    species = (
        ((:S, ((:C, 1),)),),          # substrates
        ((:P, ((:C, 1),)),),          # products
        (),                            # regulators
        ((:E, ()), (:ES, ((:C, 1),))),  # enzymes
    )
    rxns = (
        ((:E, :S), (:ES,)),
        ((:ES,), (:E, :P)),
    )
    m = EnzymeMechanism(species, rxns)

    @testset "Structure" begin
        @test n_states(m) == 2
        @test length(enzyme_forms(m)) == 2
        @test length(metabolites(m)) == 2
        @test Set(e[1] for e in enzyme_forms(m)) == Set([:E, :ES])
        @test Set(mt[1] for mt in metabolites(m)) == Set([:S, :P])
    end

    @testset "Validation" begin
        species_bad = (
            ((:S, ((:C, 2),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 2),))),
        )
        rxns_bad = (
            ((:E, :S), (:ES,)),
            ((:ES,), (:E, :P)),
        )
        @test_throws ErrorException EnzymeMechanism(species_bad, rxns_bad)
    end

    @testset "Independent params" begin
        @test n_independent_params(m) == 3
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(1001)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S, :P]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            p_pkg = merge(params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_uniuni(p, concs) rtol=1e-12
        end
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(42)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S, :P]; rng=rng)
            v_ka = rate_equation(m, params, concs)
            v_ref = reference_qssa(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-10
        end
    end

    @testset "Haldane relation" begin
        params = (k1f=3.2, k1r=0.8, k2f=2.5, k2r=1.1, E_total=1.0)
        Keq = params.k1f * params.k2f / (params.k1r * params.k2r)
        S_eq = 1.0
        P_eq = Keq * S_eq
        v = rate_equation(m, params, (S=S_eq, P=P_eq))
        @test abs(v) < 1e-12
    end

    @testset "Performance" begin
        allocs, t = test_rate_equation_performance(m,
            (k1f=3.2, k1r=0.8, k2f=2.5, k2r=1.1, E_total=1.0), (S=1.0, P=0.5))
        @test allocs == 0
        @test t < 100e-9
    end
end

@testset "Seq Uni-Bi" begin
    function rate_seq_unibi(p, c)
        S1 = c.S1; P1 = c.P1; P2 = c.P2
        f = [p.k1f * S1, p.k2f, p.k3f, p.k4f]
        r = [p.k1r, p.k2r, p.k3r * P1, p.k4r * P2]
        return _unicyclic_flux(f, r, p.Et)
    end

    species = (
        ((:S1, ((:C, 1), (:H, 1))),),
        ((:P1, ((:C, 1),)), (:P2, ((:H, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1), (:H, 1))), (:EP1P2, ((:C, 1), (:H, 1))), (:EP2, ((:H, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)),
        ((:ES1,), (:EP1P2,)),
        ((:EP1P2,), (:EP2, :P1)),
        ((:EP2,), (:E, :P2)),
    )
    m = EnzymeMechanism(species, rxns)

    @testset "Structure" begin
        @test n_states(m) == 4
        @test length(metabolites(m)) == 3
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(1002)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :P1, :P2]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            p_pkg = merge(params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_unibi(p, concs) rtol=1e-10
        end
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(123)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :P1, :P2]; rng=rng)
            v_ka = rate_equation(m, params, concs)
            v_ref = reference_qssa(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-8
        end
    end

    @testset "Haldane relation" begin
        params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6, E_total=1.0)
        Keq = prod(params[Symbol("k$(i)f")] for i in 1:4) /
              prod(params[Symbol("k$(i)r")] for i in 1:4)
        S1_eq = 1.0
        P_prod = Keq * S1_eq
        P1_eq = sqrt(P_prod)
        P2_eq = sqrt(P_prod)
        v = rate_equation(m, params, (S1=S1_eq, P1=P1_eq, P2=P2_eq))
        @test abs(v) < 1e-10
    end

    @testset "Performance" begin
        allocs, t = test_rate_equation_performance(m,
            (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6, E_total=1.0),
            (S1=1.0, P1=0.5, P2=0.3))
        @test allocs == 0
        @test t < 100e-9
    end
end

@testset "Ping-Pong Bi-Bi" begin
    species = (
        ((:A, ((:C, 2), (:N, 1))), (:B, ((:C, 3),))),
        ((:P, ((:C, 2),)), (:Q, ((:C, 3), (:N, 1)))),
        (),
        ((:E, ()), (:EA, ((:C, 2), (:N, 1))), (:FP, ((:C, 2), (:N, 1))), (:F, ((:N, 1),)), (:FB, ((:C, 3), (:N, 1))), (:EQ, ((:C, 3), (:N, 1)))),
    )
    rxns = (
        ((:E, :A), (:EA,)), ((:EA,), (:FP,)), ((:FP,), (:F, :P)),
        ((:F, :B), (:FB,)), ((:FB,), (:EQ,)), ((:EQ,), (:E, :Q)),
    )
    m = EnzymeMechanism(species, rxns)

    @testset "Structure" begin
        @test n_states(m) == 6
        @test Set(e[1] for e in enzyme_forms(m)) == Set([:E, :EA, :FP, :F, :FB, :EQ])
    end

    @testset "Validation" begin
    end

    @testset "Numerator structure" begin
        s = rate_equation_string(m)
        @test occursin("A", s)
        @test occursin("B", s)
    end

    @testset "Denominator structure" begin
        base_params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3,
            k4f=2.5, k4r=0.6, k5f=1.8, k5r=0.2, k6f=3.5, k6r=0.7, E_total=1.0)
        eps_val = 1e-10
        concs = (A=eps_val, P=eps_val, B=eps_val, Q=eps_val)
        v = rate_equation(m, base_params, concs)
        pf = prod(base_params[Symbol("k$(i)f")] for i in 1:6)
        pr = prod(base_params[Symbol("k$(i)r")] for i in 1:6)
        num = pf * eps_val * eps_val - pr * eps_val * eps_val
        denom = num / v
        @test abs(denom) < 1e-5

        s = rate_equation_string(m)
        @test occursin("A", s)
        @test occursin("B", s)
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(456)
        for _ in 1:10
            params, concs = random_params_concs(m, [:A, :P, :B, :Q]; rng=rng)
            v_ka = rate_equation(m, params, concs)
            v_ref = reference_qssa(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-8
        end
    end

    @testset "Performance" begin
        allocs, t = test_rate_equation_performance(m,
            (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3,
             k4f=2.5, k4r=0.6, k5f=1.8, k5r=0.2, k6f=3.5, k6r=0.7, E_total=1.0),
            (A=1.0, P=0.5, B=0.8, Q=0.3))
        @test allocs == 0
        @test t < 100e-9
    end
end

@testset "Seq Bi-Uni" begin
    function rate_seq_biuni(p, c)
        S1 = c.S1; S2 = c.S2; P1 = c.P1
        f = [p.k1f * S1, p.k2f * S2, p.k3f, p.k4f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r * P1]
        return _unicyclic_flux(f, r, p.Et)
    end

    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1),))),
        ((:P1, ((:C, 1), (:H, 1))),),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1))), (:EP1, ((:C, 1), (:H, 1)))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)),
        ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2,), (:EP1,)),
        ((:EP1,), (:E, :P1)),
    )
    m = EnzymeMechanism(species, rxns)

    @testset "Structure" begin
        @test n_states(m) == 4
        @test length(metabolites(m)) == 3
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2001)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :P1]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            p_pkg = merge(params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_biuni(p, concs) rtol=1e-10
        end
    end

    @testset "Performance" begin
        allocs, t = test_rate_equation_performance(m,
            (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6, E_total=1.0),
            (S1=1.0, S2=0.8, P1=0.5))
        @test allocs == 0
        @test t < 100e-9
    end
end

@testset "Seq Bi-Bi" begin
    function rate_seq_bibi(p, c)
        S1 = c.S1; S2 = c.S2
        P1 = c.P1; P2 = c.P2
        f = [p.k1f * S1, p.k2f * S2, p.k3f, p.k4f, p.k5f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r * P1, p.k5r * P2]
        return _unicyclic_flux(f, r, p.Et)
    end

    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1),))),
        ((:P1, ((:C, 1),)), (:P2, ((:H, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1))), (:EP1P2, ((:C, 1), (:H, 1))), (:EP2, ((:H, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)),
        ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2,), (:EP1P2,)),
        ((:EP1P2,), (:EP2, :P1)),
        ((:EP2,), (:E, :P2)),
    )
    m = EnzymeMechanism(species, rxns)

    @testset "Structure" begin
        @test n_states(m) == 5
        @test length(metabolites(m)) == 4
    end

    @testset "Denominator has S1 and S2" begin
        s = rate_equation_string(m)
        @test occursin(r"\bS1\b", s)
        @test occursin(r"\bS2\b", s)
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2002)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            p_pkg = merge(params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_bibi(p, concs) rtol=1e-10
        end
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(789)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2]; rng=rng)
            v_ka = rate_equation(m, params, concs)
            v_ref = reference_qssa(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-8
        end
    end

    @testset "Haldane relation" begin
        params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6, k5f=2.5, k5r=0.7, E_total=1.0)
        Keq = prod(params[Symbol("k$(i)f")] for i in 1:5) /
              prod(params[Symbol("k$(i)r")] for i in 1:5)
        S1_eq = 1.0; S2_eq = 1.0
        P_prod = Keq * S1_eq * S2_eq
        P1_eq = sqrt(P_prod); P2_eq = sqrt(P_prod)
        v = rate_equation(m, params, (S1=S1_eq, S2=S2_eq, P1=P1_eq, P2=P2_eq))
        @test abs(v) < 1e-10
    end

    @testset "Performance" begin
        allocs, t = test_rate_equation_performance(m,
            (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6, k5f=2.5, k5r=0.7, E_total=1.0),
            (S1=1.0, S2=0.8, P1=0.5, P2=0.3))
        @test allocs == 0
        @test t < 100e-9
    end
end

@testset "Seq Bi-Ter" begin
    function rate_seq_biter(p, c)
        S1 = c.S1; S2 = c.S2
        P1 = c.P1; P2 = c.P2; P3 = c.P3
        f = [p.k1f * S1, p.k2f * S2, p.k3f, p.k4f, p.k5f, p.k6f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r * P1, p.k5r * P2, p.k6r * P3]
        return _unicyclic_flux(f, r, p.Et)
    end

    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1), (:N, 1)))),
        ((:P1, ((:C, 1),)), (:P2, ((:H, 1),)), (:P3, ((:N, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1), (:N, 1))),
         (:EP1P2P3, ((:C, 1), (:H, 1), (:N, 1))), (:EP2P3, ((:H, 1), (:N, 1))), (:EP3, ((:N, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)),
        ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2,), (:EP1P2P3,)),
        ((:EP1P2P3,), (:EP2P3, :P1)),
        ((:EP2P3,), (:EP3, :P2)),
        ((:EP3,), (:E, :P3)),
    )
    m = EnzymeMechanism(species, rxns)

    @testset "Structure" begin
        @test n_states(m) == 6
        @test length(metabolites(m)) == 5
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2005)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2, :P3]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            p_pkg = merge(params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_biter(p, concs) rtol=1e-10
        end
    end

    @testset "Performance" begin
        allocs, t = test_rate_equation_performance(m,
            (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3,
             k4f=4.0, k4r=0.6, k5f=2.5, k5r=0.7, k6f=1.8, k6r=0.2, E_total=1.0),
            (S1=1.0, S2=0.8, P1=0.5, P2=0.3, P3=0.2))
        @test allocs == 0
        @test t < 100e-9
    end
end

@testset "Seq Ter-Bi" begin
    function rate_seq_terbi(p, c)
        S1 = c.S1; S2 = c.S2; S3 = c.S3
        P1 = c.P1; P2 = c.P2
        f = [p.k1f * S1, p.k2f * S2, p.k3f * S3, p.k4f, p.k5f, p.k6f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r, p.k5r * P1, p.k6r * P2]
        return _unicyclic_flux(f, r, p.Et)
    end

    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1),)), (:S3, ((:N, 1),))),
        ((:P1, ((:C, 1), (:H, 1))), (:P2, ((:N, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1))),
         (:ES1S2S3, ((:C, 1), (:H, 1), (:N, 1))), (:EP1P2, ((:C, 1), (:H, 1), (:N, 1))), (:EP2, ((:N, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)),
        ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2, :S3), (:ES1S2S3,)),
        ((:ES1S2S3,), (:EP1P2,)),
        ((:EP1P2,), (:EP2, :P1)),
        ((:EP2,), (:E, :P2)),
    )
    m = EnzymeMechanism(species, rxns)

    @testset "Structure" begin
        @test n_states(m) == 6
        @test length(metabolites(m)) == 5
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2003)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            p_pkg = merge(params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_terbi(p, concs) rtol=1e-10
        end
    end

    @testset "Performance" begin
        allocs, t = test_rate_equation_performance(m,
            (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3,
             k4f=4.0, k4r=0.6, k5f=2.5, k5r=0.7, k6f=1.8, k6r=0.2, E_total=1.0),
            (S1=1.0, S2=0.8, S3=0.6, P1=0.5, P2=0.3))
        @test allocs == 0
        @test t < 100e-9
    end
end

@testset "Seq Ter-Ter" begin
    function rate_seq_terter(p, c)
        S1 = c.S1; S2 = c.S2; S3 = c.S3
        P1 = c.P1; P2 = c.P2; P3 = c.P3
        f = [p.k1f * S1, p.k2f * S2, p.k3f * S3, p.k4f, p.k5f, p.k6f, p.k7f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r, p.k5r * P1, p.k6r * P2, p.k7r * P3]
        return _unicyclic_flux(f, r, p.Et)
    end

    species = (
        ((:S1, ((:C, 1),)), (:S2, ((:H, 1),)), (:S3, ((:N, 1),))),
        ((:P1, ((:C, 1),)), (:P2, ((:H, 1),)), (:P3, ((:N, 1),))),
        (),
        ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1))),
         (:ES1S2S3, ((:C, 1), (:H, 1), (:N, 1))), (:EP1P2P3, ((:C, 1), (:H, 1), (:N, 1))),
         (:EP2P3, ((:H, 1), (:N, 1))), (:EP3, ((:N, 1),))),
    )
    rxns = (
        ((:E, :S1), (:ES1,)),
        ((:ES1, :S2), (:ES1S2,)),
        ((:ES1S2, :S3), (:ES1S2S3,)),
        ((:ES1S2S3,), (:EP1P2P3,)),
        ((:EP1P2P3,), (:EP2P3, :P1)),
        ((:EP2P3,), (:EP3, :P2)),
        ((:EP3,), (:E, :P3)),
    )
    m = EnzymeMechanism(species, rxns)

    @testset "Structure" begin
        @test n_states(m) == 7
        @test length(metabolites(m)) == 6
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2004)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2, :P3]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            p_pkg = merge(params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_terter(p, concs) rtol=1e-10
        end
    end

    @testset "Performance" begin
        allocs, t = test_rate_equation_performance(m,
            (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3,
             k4f=4.0, k4r=0.6, k5f=2.5, k5r=0.7, k6f=1.8, k6r=0.2, k7f=3.0, k7r=0.9, E_total=1.0),
            (S1=1.0, S2=0.8, S3=0.6, P1=0.5, P2=0.3, P3=0.2))
        @test allocs == 0
        @test t < 100e-9
    end
end

@testset "Random-order Bi-Bi (branched)" begin
    species = (
        ((:A, ((:C, 1),)), (:B, ((:N, 1),))),
        ((:P, ((:C, 1),)), (:Q, ((:N, 1),))),
        (),
        ((:E, ()), (:EA, ((:C, 1),)), (:EB, ((:N, 1),)),
         (:EAB, ((:C, 1), (:N, 1))), (:EPQ, ((:C, 1), (:N, 1))), (:EQ, ((:N, 1),))),
    )
    rxns = (
        ((:E, :A), (:EA,)),
        ((:E, :B), (:EB,)),
        ((:EA, :B), (:EAB,)),
        ((:EB, :A), (:EAB,)),
        ((:EAB,), (:EPQ,)),
        ((:EPQ,), (:EQ, :P)),
        ((:EQ,), (:E, :Q)),
    )
    m = EnzymeMechanism(species, rxns)

    @testset "Structure" begin
        @test n_states(m) == 6
        @test length(metabolites(m)) == 4
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(3001)
        for _ in 1:20
            params, concs = random_params_concs(m, [:A, :B, :P, :Q]; rng=rng)
            v_ka = rate_equation(m, params, concs)
            v_ref = reference_qssa(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-8
        end
    end

    @testset "Equilibrium (zero flux)" begin
        params = (k1f=2.0, k1r=0.5, k2f=1.5, k2r=0.3, k3f=3.0, k3r=0.4,
                  k4f=2.5, k4r=0.6, k5f=1.8, k5r=0.2, k6f=3.5, k6r=0.7, k7f=2.0, k7r=0.8, E_total=1.0)
        rng = Random.MersenneTwister(3002)
        for _ in 1:10
            _, concs = random_params_concs(m, [:A, :B, :P, :Q]; rng=rng)
            v_ka = rate_equation(m, params, concs)
            v_ref = reference_qssa(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-10
        end
    end

    @testset "Performance" begin
        allocs, t = test_rate_equation_performance(m,
            (k1f=2.0, k1r=0.5, k2f=1.5, k2r=0.3, k3f=3.0, k3r=0.4,
             k4f=2.5, k4r=0.6, k5f=1.8, k5r=0.2, k6f=3.5, k6r=0.7, k7f=2.0, k7r=0.8, E_total=1.0),
            (A=1.0, B=0.8, P=0.5, Q=0.3))
        @test allocs == 0
        @test t < 100e-9
    end
end
