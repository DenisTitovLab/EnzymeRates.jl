@testset "Uni-Uni reversible" begin
    E  = Species(:E,  enzyme)
    ES = Species(:ES, enzyme)
    S  = Species(:S,  metabolite, Dict(:C => 1))
    P  = Species(:P,  metabolite, Dict(:C => 1))

    m = EnzymeMechanism([[E, S] => [ES], [ES] => [E, P]])

    @testset "Structure" begin
        @test n_states(m) == 2
        @test length(enzyme_forms(m)) == 2
        @test length(metabolites(m)) == 2
        @test Set(s.name for s in enzyme_forms(m)) == Set([:E, :ES])
        @test Set(s.name for s in metabolites(m)) == Set([:S, :P])
    end

    @testset "Validation" begin
        @test validate(m) == true
        # Mutate atoms to break conservation
        S_bad = Species(:S, metabolite, Dict(:C => 2))
        m_bad = EnzymeMechanism([[E, S_bad] => [ES], [ES] => [E, P]])
        @test validate(m_bad) == false
    end

    @testset "Independent params" begin
        # 2 steps, 4 raw params, 1 Haldane constraint
        @test n_independent_params(m) == 3
    end

    @testset "Spot check v ≈ 0.9091" begin
        params = (k1f=3.2, k1r=0.8, k2f=2.5, k2r=1.1)
        concs = (S=0.7, P=0.3)
        fn = rate_function(m)
        v = fn(params, concs)
        @test v ≈ 0.9091 atol=0.001
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(42)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S, :P]; rng=rng)
            fn = rate_function(m)
            v_ka = fn(params, concs)
            v_ref = reference_king_altman(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-10
        end
    end

    @testset "Rate equation string" begin
        s = rate_equation_string(m)
        @test occursin("k1f", s)
        @test occursin("k2f", s)
        @test occursin("S", s)
        @test occursin("P", s)
    end

    @testset "Haldane relation" begin
        params = (k1f=3.2, k1r=0.8, k2f=2.5, k2r=1.1)
        Keq = params.k1f * params.k2f / (params.k1r * params.k2r)
        S_eq = 1.0
        P_eq = Keq * S_eq
        fn = rate_function(m)
        v = fn(params, (S=S_eq, P=P_eq))
        @test abs(v) < 1e-12
    end
end

@testset "Uni-Bi ordered" begin
    E   = Species(:E,   enzyme)
    ES  = Species(:ES,  enzyme)
    EPQ = Species(:EPQ, enzyme)
    EQ  = Species(:EQ,  enzyme)
    S   = Species(:S,   metabolite, Dict(:C => 1, :H => 1))
    P   = Species(:P,   metabolite, Dict(:C => 1))
    Q   = Species(:Q,   metabolite, Dict(:H => 1))

    m = EnzymeMechanism([
        [E, S] => [ES], [ES] => [EPQ], [EPQ] => [EQ, P], [EQ] => [E, Q]
    ])

    @testset "Structure" begin
        @test n_states(m) == 4
        @test Set(s.name for s in enzyme_forms(m)) == Set([:E, :ES, :EPQ, :EQ])
    end

    @testset "Validation" begin
        @test validate(m) == true
    end

    @testset "Spot check v ≈ 0.4625" begin
        params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6)
        concs = (S=0.8, P=0.5, Q=0.3)
        fn = rate_function(m)
        v = fn(params, concs)
        @test v ≈ 0.4625 atol=0.001
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(123)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S, :P, :Q]; rng=rng)
            fn = rate_function(m)
            v_ka = fn(params, concs)
            v_ref = reference_king_altman(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-8
        end
    end

    @testset "Haldane relation" begin
        params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6)
        Keq = (params.k1f * params.k2f * params.k3f * params.k4f) /
              (params.k1r * params.k2r * params.k3r * params.k4r)
        S_eq = 1.0
        PQ = Keq * S_eq
        P_eq = sqrt(PQ)
        Q_eq = sqrt(PQ)
        fn = rate_function(m)
        v = fn(params, (S=S_eq, P=P_eq, Q=Q_eq))
        @test abs(v) < 1e-10
    end
end

@testset "Ping-Pong Bi-Bi" begin
    E  = Species(:E,  enzyme)
    EA = Species(:EA, enzyme)
    FP = Species(:FP, enzyme)
    F  = Species(:F,  enzyme)
    FB = Species(:FB, enzyme)
    EQ = Species(:EQ, enzyme)

    A = Species(:A, metabolite, Dict(:C => 2, :N => 1))
    P_met = Species(:P, metabolite, Dict(:C => 2))
    B = Species(:B, metabolite, Dict(:C => 3))
    Q = Species(:Q, metabolite, Dict(:C => 3, :N => 1))

    m = EnzymeMechanism([
        [E, A] => [EA], [EA] => [FP], [FP] => [F, P_met],
        [F, B] => [FB], [FB] => [EQ], [EQ] => [E, Q]
    ])

    @testset "Structure" begin
        @test n_states(m) == 6
        @test Set(s.name for s in enzyme_forms(m)) == Set([:E, :EA, :FP, :F, :FB, :EQ])
    end

    @testset "Validation" begin
        @test validate(m) == true
    end

    @testset "Numerator structure" begin
        s = rate_equation_string(m)
        @test occursin("A", s)
        @test occursin("B", s)
    end

    @testset "Denominator structure" begin
        fn = rate_function(m)
        base_params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3,
                       k4f=2.5, k4r=0.6, k5f=1.8, k5r=0.2, k6f=3.5, k6r=0.7)
        eps_val = 1e-10
        concs = (A=eps_val, P=eps_val, B=eps_val, Q=eps_val)
        v = fn(base_params, concs)
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
            fn = rate_function(m)
            v_ka = fn(params, concs)
            v_ref = reference_king_altman(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-8
        end
    end
end

@testset "Bi-Bi Sequential ordered" begin
    E   = Species(:E,   enzyme)
    EA  = Species(:EA,  enzyme)
    EAB = Species(:EAB, enzyme)
    EPQ = Species(:EPQ, enzyme)
    EQ  = Species(:EQ,  enzyme)

    A = Species(:A, metabolite, Dict(:C => 2))
    B = Species(:B, metabolite, Dict(:C => 3))
    P_met = Species(:P, metabolite, Dict(:C => 2))
    Q = Species(:Q, metabolite, Dict(:C => 3))

    m = EnzymeMechanism([
        [E, A] => [EA], [EA, B] => [EAB], [EAB] => [EPQ],
        [EPQ] => [EQ, P_met], [EQ] => [E, Q]
    ])

    @testset "Structure" begin
        @test n_states(m) == 5
        @test Set(s.name for s in enzyme_forms(m)) == Set([:E, :EA, :EAB, :EPQ, :EQ])
    end

    @testset "Validation" begin
        @test validate(m) == true
    end

    @testset "Denominator has A*B term" begin
        s = rate_equation_string(m)
        denom_start = findfirst('/', s)
        denom_str = s[denom_start+1:end]
        terms = split(denom_str, " + ")
        has_AB = any(term -> occursin(r"\bA\b", term) && occursin(r"\bB\b", term), terms)
        @test has_AB
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(789)
        for _ in 1:10
            params, concs = random_params_concs(m, [:A, :B, :P, :Q]; rng=rng)
            fn = rate_function(m)
            v_ka = fn(params, concs)
            v_ref = reference_king_altman(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-8
        end
    end

    @testset "Haldane relation" begin
        params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6, k5f=2.5, k5r=0.7)
        Keq = prod(params[Symbol("k$(i)f")] for i in 1:5) /
              prod(params[Symbol("k$(i)r")] for i in 1:5)
        A_eq = 1.0; B_eq = 1.0
        PQ = Keq * A_eq * B_eq
        P_eq = sqrt(PQ); Q_eq = sqrt(PQ)
        fn = rate_function(m)
        v = fn(params, (A=A_eq, B=B_eq, P=P_eq, Q=Q_eq))
        @test abs(v) < 1e-10
    end
end
