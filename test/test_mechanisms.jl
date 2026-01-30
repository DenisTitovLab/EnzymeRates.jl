@testset "Uni-Uni" begin
    function rate_seq_uniuni(p, c)
        S = _getc(c, :S); P = _getc(c, :P)
        f = [p.k1f * S, p.k2f]
        r = [p.k1r, p.k2r * P]
        return _unicyclic_flux(f, r, p.Et)
    end

    E = Species(:E, enzyme)
    ES = Species(:ES, enzyme)
    S = Species(:S, metabolite, Dict(:C => 1))
    P = Species(:P, metabolite, Dict(:C => 1))

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
        S_bad = Species(:S, metabolite, Dict(:C => 2))
        m_bad = EnzymeMechanism([[E, S_bad] => [ES], [ES] => [E, P]])
        @test validate(m_bad) == false
    end

    @testset "Independent params" begin
        @test n_independent_params(m) == 3
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(1001)
        fn = rate_function(m)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S, :P]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            c_pkg = merge(concs, (E_total=Et,))
            @test fn(params, c_pkg) ≈ rate_seq_uniuni(p, concs) rtol=1e-12
        end
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

@testset "Seq Uni-Bi" begin
    function rate_seq_unibi(p, c)
        S1 = _getc(c, :S1); P1 = _getc(c, :P1); P2 = _getc(c, :P2)
        f = [p.k1f * S1, p.k2f, p.k3f, p.k4f]
        r = [p.k1r, p.k2r, p.k3r * P1, p.k4r * P2]
        return _unicyclic_flux(f, r, p.Et)
    end

    E     = Species(:E,     enzyme)
    ES1   = Species(:ES1,   enzyme)
    EP1P2 = Species(:EP1P2, enzyme)
    EP2   = Species(:EP2,   enzyme)
    S1    = Species(:S1, metabolite, Dict(:C => 1, :H => 1))
    P1    = Species(:P1, metabolite, Dict(:C => 1))
    P2    = Species(:P2, metabolite, Dict(:H => 1))

    m = EnzymeMechanism([
        [E, S1] => [ES1],
        [ES1] => [EP1P2],
        [EP1P2] => [EP2, P1],
        [EP2] => [E, P2]
    ])

    @testset "Structure" begin
        @test n_states(m) == 4
        @test length(metabolites(m)) == 3
        @test validate(m) == true
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(1002)
        fn = rate_function(m)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :P1, :P2]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            c_pkg = merge(concs, (E_total=Et,))
            @test fn(params, c_pkg) ≈ rate_seq_unibi(p, concs) rtol=1e-10
        end
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(123)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :P1, :P2]; rng=rng)
            fn = rate_function(m)
            v_ka = fn(params, concs)
            v_ref = reference_king_altman(m, params, concs)
            @test v_ka ≈ v_ref rtol=1e-8
        end
    end

    @testset "Haldane relation" begin
        params = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6)
        Keq = prod(params[Symbol("k$(i)f")] for i in 1:4) /
              prod(params[Symbol("k$(i)r")] for i in 1:4)
        S1_eq = 1.0
        P_prod = Keq * S1_eq
        P1_eq = sqrt(P_prod)
        P2_eq = sqrt(P_prod)
        fn = rate_function(m)
        v = fn(params, (S1=S1_eq, P1=P1_eq, P2=P2_eq))
        @test abs(v) < 1e-10
    end
end

@testset "Ping-Pong Bi-Bi" begin
    E = Species(:E, enzyme)
    EA = Species(:EA, enzyme)
    FP = Species(:FP, enzyme)
    F = Species(:F, enzyme)
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

@testset "Seq Bi-Uni" begin
    function rate_seq_biuni(p, c)
        S1 = _getc(c, :S1); S2 = _getc(c, :S2); P1 = _getc(c, :P1)
        f = [p.k1f * S1, p.k2f * S2, p.k3f, p.k4f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r * P1]
        return _unicyclic_flux(f, r, p.Et)
    end

    E      = Species(:E,      enzyme)
    ES1    = Species(:ES1,    enzyme)
    ES1S2  = Species(:ES1S2,  enzyme)
    EP1    = Species(:EP1,    enzyme)
    S1     = Species(:S1, metabolite, Dict(:C => 1))
    S2     = Species(:S2, metabolite, Dict(:H => 1))
    P1     = Species(:P1, metabolite, Dict(:C => 1, :H => 1))

    m = EnzymeMechanism([
        [E, S1] => [ES1],
        [ES1, S2] => [ES1S2],
        [ES1S2] => [EP1],
        [EP1] => [E, P1]
    ])

    @testset "Structure" begin
        @test n_states(m) == 4
        @test length(metabolites(m)) == 3
        @test validate(m) == true
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2001)
        fn = rate_function(m)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :P1]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            c_pkg = merge(concs, (E_total=Et,))
            @test fn(params, c_pkg) ≈ rate_seq_biuni(p, concs) rtol=1e-10
        end
    end
end

@testset "Seq Bi-Bi" begin
    function rate_seq_bibi(p, c)
        S1 = _getc(c, :S1); S2 = _getc(c, :S2)
        P1 = _getc(c, :P1); P2 = _getc(c, :P2)
        f = [p.k1f * S1, p.k2f * S2, p.k3f, p.k4f, p.k5f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r * P1, p.k5r * P2]
        return _unicyclic_flux(f, r, p.Et)
    end

    E      = Species(:E,      enzyme)
    ES1    = Species(:ES1,    enzyme)
    ES1S2  = Species(:ES1S2,  enzyme)
    EP1P2  = Species(:EP1P2,  enzyme)
    EP2    = Species(:EP2,    enzyme)
    S1     = Species(:S1, metabolite, Dict(:C => 1))
    S2     = Species(:S2, metabolite, Dict(:H => 1))
    P1     = Species(:P1, metabolite, Dict(:C => 1))
    P2     = Species(:P2, metabolite, Dict(:H => 1))

    m = EnzymeMechanism([
        [E, S1] => [ES1],
        [ES1, S2] => [ES1S2],
        [ES1S2] => [EP1P2],
        [EP1P2] => [EP2, P1],
        [EP2] => [E, P2]
    ])

    @testset "Structure" begin
        @test n_states(m) == 5
        @test length(metabolites(m)) == 4
        @test validate(m) == true
    end

    @testset "Denominator has S1*S2 term" begin
        s = rate_equation_string(m)
        denom_start = findfirst('/', s)
        denom_str = s[denom_start+1:end]
        terms = split(denom_str, " + ")
        has_S1S2 = any(term -> occursin(r"\bS1\b", term) && occursin(r"\bS2\b", term), terms)
        @test has_S1S2
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2002)
        fn = rate_function(m)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            c_pkg = merge(concs, (E_total=Et,))
            @test fn(params, c_pkg) ≈ rate_seq_bibi(p, concs) rtol=1e-10
        end
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(789)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2]; rng=rng)
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
        S1_eq = 1.0; S2_eq = 1.0
        P_prod = Keq * S1_eq * S2_eq
        P1_eq = sqrt(P_prod); P2_eq = sqrt(P_prod)
        fn = rate_function(m)
        v = fn(params, (S1=S1_eq, S2=S2_eq, P1=P1_eq, P2=P2_eq))
        @test abs(v) < 1e-10
    end
end

@testset "Seq Bi-Ter" begin
    function rate_seq_biter(p, c)
        S1 = _getc(c, :S1); S2 = _getc(c, :S2)
        P1 = _getc(c, :P1); P2 = _getc(c, :P2); P3 = _getc(c, :P3)
        f = [p.k1f * S1, p.k2f * S2, p.k3f, p.k4f, p.k5f, p.k6f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r * P1, p.k5r * P2, p.k6r * P3]
        return _unicyclic_flux(f, r, p.Et)
    end

    E        = Species(:E,        enzyme)
    ES1      = Species(:ES1,      enzyme)
    ES1S2    = Species(:ES1S2,    enzyme)
    EP1P2P3  = Species(:EP1P2P3,  enzyme)
    EP2P3    = Species(:EP2P3,    enzyme)
    EP3      = Species(:EP3,      enzyme)
    S1       = Species(:S1, metabolite, Dict(:C => 1))
    S2       = Species(:S2, metabolite, Dict(:H => 1, :N => 1))
    P1       = Species(:P1, metabolite, Dict(:C => 1))
    P2       = Species(:P2, metabolite, Dict(:H => 1))
    P3       = Species(:P3, metabolite, Dict(:N => 1))

    m = EnzymeMechanism([
        [E, S1] => [ES1],
        [ES1, S2] => [ES1S2],
        [ES1S2] => [EP1P2P3],
        [EP1P2P3] => [EP2P3, P1],
        [EP2P3] => [EP3, P2],
        [EP3] => [E, P3]
    ])

    @testset "Structure" begin
        @test n_states(m) == 6
        @test length(metabolites(m)) == 5
        @test validate(m) == true
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2005)
        fn = rate_function(m)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2, :P3]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            c_pkg = merge(concs, (E_total=Et,))
            @test fn(params, c_pkg) ≈ rate_seq_biter(p, concs) rtol=1e-10
        end
    end
end

@testset "Seq Ter-Bi" begin
    function rate_seq_terbi(p, c)
        S1 = _getc(c, :S1); S2 = _getc(c, :S2); S3 = _getc(c, :S3)
        P1 = _getc(c, :P1); P2 = _getc(c, :P2)
        f = [p.k1f * S1, p.k2f * S2, p.k3f * S3, p.k4f, p.k5f, p.k6f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r, p.k5r * P1, p.k6r * P2]
        return _unicyclic_flux(f, r, p.Et)
    end

    E        = Species(:E,        enzyme)
    ES1      = Species(:ES1,      enzyme)
    ES1S2    = Species(:ES1S2,    enzyme)
    ES1S2S3  = Species(:ES1S2S3,  enzyme)
    EP1P2    = Species(:EP1P2,    enzyme)
    EP2      = Species(:EP2,      enzyme)
    S1       = Species(:S1, metabolite, Dict(:C => 1))
    S2       = Species(:S2, metabolite, Dict(:H => 1))
    S3       = Species(:S3, metabolite, Dict(:N => 1))
    P1       = Species(:P1, metabolite, Dict(:C => 1, :H => 1))
    P2       = Species(:P2, metabolite, Dict(:N => 1))

    m = EnzymeMechanism([
        [E, S1] => [ES1],
        [ES1, S2] => [ES1S2],
        [ES1S2, S3] => [ES1S2S3],
        [ES1S2S3] => [EP1P2],
        [EP1P2] => [EP2, P1],
        [EP2] => [E, P2]
    ])

    @testset "Structure" begin
        @test n_states(m) == 6
        @test length(metabolites(m)) == 5
        @test validate(m) == true
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2003)
        fn = rate_function(m)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            c_pkg = merge(concs, (E_total=Et,))
            @test fn(params, c_pkg) ≈ rate_seq_terbi(p, concs) rtol=1e-10
        end
    end
end

@testset "Seq Ter-Ter" begin
    function rate_seq_terter(p, c)
        S1 = _getc(c, :S1); S2 = _getc(c, :S2); S3 = _getc(c, :S3)
        P1 = _getc(c, :P1); P2 = _getc(c, :P2); P3 = _getc(c, :P3)
        f = [p.k1f * S1, p.k2f * S2, p.k3f * S3, p.k4f, p.k5f, p.k6f, p.k7f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r, p.k5r * P1, p.k6r * P2, p.k7r * P3]
        return _unicyclic_flux(f, r, p.Et)
    end

    E          = Species(:E,          enzyme)
    ES1        = Species(:ES1,        enzyme)
    ES1S2      = Species(:ES1S2,      enzyme)
    ES1S2S3    = Species(:ES1S2S3,    enzyme)
    EP1P2P3    = Species(:EP1P2P3,    enzyme)
    EP2P3      = Species(:EP2P3,      enzyme)
    EP3        = Species(:EP3,        enzyme)
    S1         = Species(:S1, metabolite, Dict(:C => 1))
    S2         = Species(:S2, metabolite, Dict(:H => 1))
    S3         = Species(:S3, metabolite, Dict(:N => 1))
    P1         = Species(:P1, metabolite, Dict(:C => 1))
    P2         = Species(:P2, metabolite, Dict(:H => 1))
    P3         = Species(:P3, metabolite, Dict(:N => 1))

    m = EnzymeMechanism([
        [E, S1] => [ES1],
        [ES1, S2] => [ES1S2],
        [ES1S2, S3] => [ES1S2S3],
        [ES1S2S3] => [EP1P2P3],
        [EP1P2P3] => [EP2P3, P1],
        [EP2P3] => [EP3, P2],
        [EP3] => [E, P3]
    ])

    @testset "Structure" begin
        @test n_states(m) == 7
        @test length(metabolites(m)) == 6
        @test validate(m) == true
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2004)
        fn = rate_function(m)
        for _ in 1:20
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2, :P3]; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(params, (Et=Et,))
            c_pkg = merge(concs, (E_total=Et,))
            @test fn(params, c_pkg) ≈ rate_seq_terter(p, concs) rtol=1e-10
        end
    end
end
