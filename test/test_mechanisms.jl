@testset "Uni-Uni" begin
    function rate_seq_uniuni(p, c)
        S = c.S; P = c.P
        f = [p.k1f * S, p.k2f]
        r = [p.k1r, p.k2r * P]
        return _unicyclic_flux(f, r, p.Et)
    end

    m, met_names = make_uni_uni()

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

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(1001)
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(all_params, (Et=Et,))
            p_pkg = merge(new_params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_uniuni(p, concs) rtol=1e-12
        end
    end

    @testset "Reference comparison" begin
        test_reference_comparison(m, met_names; seed=42, rtol=1e-10)
    end

    @testset "Haldane relation" begin
        Keq = 3.2 * 2.5 / (0.8 * 1.1)
        params = (k1f=3.2, k2f=2.5, k2r=1.1, Keq=Keq, E_total=1.0)
        S_eq = 1.0
        P_eq = Keq * S_eq
        v = rate_equation(m, params, (S=S_eq, P=P_eq))
        @test abs(v) < 1e-12
    end

    @testset "Performance" begin
        test_performance(m, met_names)
    end
end

@testset "Seq Uni-Bi" begin
    function rate_seq_unibi(p, c)
        S1 = c.S1; P1 = c.P1; P2 = c.P2
        f = [p.k1f * S1, p.k2f, p.k3f, p.k4f]
        r = [p.k1r, p.k2r, p.k3r * P1, p.k4r * P2]
        return _unicyclic_flux(f, r, p.Et)
    end

    m, met_names = make_seq_unibi()

    @testset "Structure" begin
        @test n_states(m) == 4
        @test length(metabolites(m)) == 3
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(1002)
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(all_params, (Et=Et,))
            p_pkg = merge(new_params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_unibi(p, concs) rtol=1e-10
        end
    end

    @testset "Reference comparison" begin
        test_reference_comparison(m, met_names; seed=123)
    end

    @testset "Haldane relation" begin
        all_p = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6, E_total=1.0)
        Keq = prod(all_p[Symbol("k$(i)f")] for i in 1:4) /
              prod(all_p[Symbol("k$(i)r")] for i in 1:4)
        params = make_independent_params(m, all_p, Keq)
        S1_eq = 1.0
        P_prod = Keq * S1_eq
        P1_eq = sqrt(P_prod)
        P2_eq = sqrt(P_prod)
        v = rate_equation(m, params, (S1=S1_eq, P1=P1_eq, P2=P2_eq))
        @test abs(v) < 1e-10
    end

    @testset "Performance" begin
        test_performance(m, met_names; seed=123)
    end
end

@testset "Ping-Pong Bi-Bi" begin
    m, met_names = make_pingpong_bibi()

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
        all_p = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3,
            k4f=2.5, k4r=0.6, k5f=1.8, k5r=0.2, k6f=3.5, k6r=0.7, E_total=1.0)
        Keq = prod(all_p[Symbol("k$(i)f")] for i in 1:6) /
              prod(all_p[Symbol("k$(i)r")] for i in 1:6)
        base_params = make_independent_params(m, all_p, Keq)
        eps_val = 1e-10
        concs = (A=eps_val, P=eps_val, B=eps_val, Q=eps_val)
        v = rate_equation(m, base_params, concs)
        pf = prod(all_p[Symbol("k$(i)f")] for i in 1:6)
        pr = prod(all_p[Symbol("k$(i)r")] for i in 1:6)
        num = pf * eps_val * eps_val - pr * eps_val * eps_val
        denom = num / v
        @test abs(denom) < 1e-5

        s = rate_equation_string(m)
        @test occursin("A", s)
        @test occursin("B", s)
    end

    @testset "Reference comparison" begin
        test_reference_comparison(m, met_names; seed=456)
    end

    @testset "Performance" begin
        test_performance(m, met_names; seed=456)
    end
end

@testset "Seq Bi-Uni" begin
    function rate_seq_biuni(p, c)
        S1 = c.S1; S2 = c.S2; P1 = c.P1
        f = [p.k1f * S1, p.k2f * S2, p.k3f, p.k4f]
        r = [p.k1r, p.k2r, p.k3r, p.k4r * P1]
        return _unicyclic_flux(f, r, p.Et)
    end

    m, met_names = make_seq_biuni()

    @testset "Structure" begin
        @test n_states(m) == 4
        @test length(metabolites(m)) == 3
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2001)
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(all_params, (Et=Et,))
            p_pkg = merge(new_params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_biuni(p, concs) rtol=1e-10
        end
    end

    @testset "Performance" begin
        test_performance(m, met_names; seed=2001)
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

    m, met_names = make_seq_bibi()

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
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(all_params, (Et=Et,))
            p_pkg = merge(new_params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_bibi(p, concs) rtol=1e-10
        end
    end

    @testset "Reference comparison" begin
        test_reference_comparison(m, met_names; seed=789)
    end

    @testset "Haldane relation" begin
        all_p = (k1f=2.0, k1r=0.5, k2f=3.0, k2r=0.4, k3f=1.5, k3r=0.3, k4f=4.0, k4r=0.6, k5f=2.5, k5r=0.7, E_total=1.0)
        Keq = prod(all_p[Symbol("k$(i)f")] for i in 1:5) /
              prod(all_p[Symbol("k$(i)r")] for i in 1:5)
        params = make_independent_params(m, all_p, Keq)
        S1_eq = 1.0; S2_eq = 1.0
        P_prod = Keq * S1_eq * S2_eq
        P1_eq = sqrt(P_prod); P2_eq = sqrt(P_prod)
        v = rate_equation(m, params, (S1=S1_eq, S2=S2_eq, P1=P1_eq, P2=P2_eq))
        @test abs(v) < 1e-10
    end

    @testset "Performance" begin
        test_performance(m, met_names; seed=789)
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

    m, met_names = make_seq_biter()

    @testset "Structure" begin
        @test n_states(m) == 6
        @test length(metabolites(m)) == 5
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2005)
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(all_params, (Et=Et,))
            p_pkg = merge(new_params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_biter(p, concs) rtol=1e-10
        end
    end

    @testset "Performance" begin
        test_performance(m, met_names; seed=2005)
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

    m, met_names = make_seq_terbi()

    @testset "Structure" begin
        @test n_states(m) == 6
        @test length(metabolites(m)) == 5
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2003)
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(all_params, (Et=Et,))
            p_pkg = merge(new_params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_terbi(p, concs) rtol=1e-10
        end
    end

    @testset "Performance" begin
        test_performance(m, met_names; seed=2003)
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

    m, met_names = make_seq_terter()

    @testset "Structure" begin
        @test n_states(m) == 7
        @test length(metabolites(m)) == 6
    end

    @testset "Expected rate equation" begin
        rng = Random.MersenneTwister(2004)
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            Et = 0.1 + 9.9 * rand(rng)
            p = merge(all_params, (Et=Et,))
            p_pkg = merge(new_params, (E_total=Et,))
            @test rate_equation(m, p_pkg, concs) ≈ rate_seq_terter(p, concs) rtol=1e-10
        end
    end

    @testset "Performance" begin
        test_performance(m, met_names; seed=2004)
    end
end

@testset "Random-order Bi-Bi (branched)" begin
    m, met_names = make_random_bibi()

    @testset "Structure" begin
        @test n_states(m) == 6
        @test length(metabolites(m)) == 4
    end

    @testset "Reference comparison" begin
        rng = Random.MersenneTwister(3001)
        for _ in 1:20
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            v_ka = rate_equation(m, new_params, concs)
            v_ref = reference_qssa(m, all_params, concs)
            @test v_ka ≈ v_ref rtol=1e-8
        end
    end

    @testset "Equilibrium (zero flux)" begin
        rng = Random.MersenneTwister(3002)
        params, _, all_p = random_independent_params_concs(m, met_names; rng=rng)
        for _ in 1:10
            _, concs, _ = random_independent_params_concs(m, met_names; rng=rng)
            v_ka = rate_equation(m, params, concs)
            v_ref = reference_qssa(m, all_p, concs)
            @test v_ka ≈ v_ref rtol=1e-10
        end
    end

    @testset "Performance" begin
        test_performance(m, met_names; seed=5555)
    end
end

@testset "rate_equation_string formatting" begin
    # Helper: parse rate_equation_string and evaluate numerically
    function _eval_rate_string(s, params, concs)
        # Extract just the equation line (after "v = ")
        eq_line = last(split(s, "v = "; limit=2))
        bindings = vcat(
            ["$k = $(params[k])" for k in keys(params)],
            ["$k = $(concs[k])" for k in keys(concs)],
        )
        code = "let $(join(bindings, ", "))\n  $eq_line\nend"
        eval(Meta.parse(code))
    end

    @testset "General formatting properties (Uni-Uni)" begin
        m, _ = make_uni_uni()
        s = rate_equation_string(m)

        @test !occursin("params.", s)
        @test !occursin("concs.", s)
        @test !occursin("+ -", s)
        @test !occursin("- -", s)
        @test occursin("v = E_total * (", s)
        @test occursin(") / (", s)
    end

    @testset "Uni-Uni numerical equivalence" begin
        m, met_names = make_uni_uni()
        s = rate_equation_string(m)
        rng = Random.MersenneTwister(9001)
        for _ in 1:10
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            @test rate_equation(m, new_params, concs) ≈ _eval_rate_string(s, all_params, concs) rtol=1e-10
        end
    end

    @testset "Seq Uni-Bi" begin
        m, met_names = make_seq_unibi()
        s = rate_equation_string(m)
        @test occursin("v = E_total * (", s)
        @test !occursin("params.", s)
        @test !occursin("concs.", s)
        @test !occursin("+ -", s)
        rng = Random.MersenneTwister(9002)
        for _ in 1:10
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            @test rate_equation(m, new_params, concs) ≈ _eval_rate_string(s, all_params, concs) rtol=1e-10
        end
    end

    @testset "Seq Bi-Uni" begin
        m, met_names = make_seq_biuni()
        s = rate_equation_string(m)
        @test occursin("v = E_total * (", s)
        @test !occursin("+ -", s)
        rng = Random.MersenneTwister(9003)
        for _ in 1:10
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            @test rate_equation(m, new_params, concs) ≈ _eval_rate_string(s, all_params, concs) rtol=1e-10
        end
    end

    @testset "Seq Bi-Bi" begin
        m, met_names = make_seq_bibi()
        s = rate_equation_string(m)
        @test occursin("v = E_total * (", s)
        @test occursin(") / (", s)
        @test !occursin("+ -", s)
        @test occursin("S1", s)
        @test occursin("S2", s)
        @test occursin("P1", s)
        @test occursin("P2", s)
        rng = Random.MersenneTwister(9004)
        for _ in 1:10
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            @test rate_equation(m, new_params, concs) ≈ _eval_rate_string(s, all_params, concs) rtol=1e-10
        end
    end

    @testset "Ping-Pong Bi-Bi" begin
        m, met_names = make_pingpong_bibi()
        s = rate_equation_string(m)
        @test occursin("v = E_total * (", s)
        @test !occursin("+ -", s)
        @test occursin("A", s)
        @test occursin("B", s)
        rng = Random.MersenneTwister(9005)
        for _ in 1:10
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            @test rate_equation(m, new_params, concs) ≈ _eval_rate_string(s, all_params, concs) rtol=1e-10
        end
    end

    @testset "Random-order Bi-Bi (branched)" begin
        m, met_names = make_random_bibi()
        s = rate_equation_string(m)
        @test occursin("v = E_total * (", s)
        @test !occursin("+ -", s)
        @test !occursin("- -", s)
        @test occursin("A", s)
        @test occursin("B", s)
        rng = Random.MersenneTwister(9006)
        for _ in 1:10
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            @test rate_equation(m, new_params, concs) ≈ _eval_rate_string(s, all_params, concs) rtol=1e-10
        end
    end

    @testset "Seq Ter-Ter" begin
        m, met_names = make_seq_terter()
        s = rate_equation_string(m)
        @test occursin("v = E_total * (", s)
        @test !occursin("+ -", s)
        rng = Random.MersenneTwister(9007)
        for _ in 1:10
            new_params, concs, all_params = random_independent_params_concs(m, met_names; rng=rng)
            @test rate_equation(m, new_params, concs) ≈ _eval_rate_string(s, all_params, concs) rtol=1e-10
        end
    end

    @testset "k-constants sorted before metabolites" begin
        m, _ = make_uni_uni()
        s = rate_equation_string(m)
        @test occursin("k1f * k2f * S", s)
        @test occursin("k1r * k2r * P", s)
    end

    @testset "Denominator has no subtraction (Uni-Uni)" begin
        m, _ = make_uni_uni()
        s = rate_equation_string(m)
        denom_start = findlast(") / (", s)
        denom = s[denom_start.stop+1:end-1]
        @test !occursin(" - ", denom)
    end
end
