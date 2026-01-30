using OrdinaryDiffEqFIRK

function build_ode_rhs(m::EnzymeMechanism, params, concs)
    forms = enzyme_forms(m)
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))

    # Precompute pseudo-first-order rates for each step
    step_data = []
    for (step_idx, (lhs, rhs)) in enumerate(m.steps)
        e_lhs = [s for s in lhs if s.role == enzyme][1]
        e_rhs = [s for s in rhs if s.role == enzyme][1]
        i = name_to_idx[e_lhs.name]
        j = name_to_idx[e_rhs.name]

        m_lhs = [s for s in lhs if s.role == metabolite]
        m_rhs = [s for s in rhs if s.role == metabolite]

        kf = Float64(params[Symbol("k$(step_idx)f")])
        kr = Float64(params[Symbol("k$(step_idx)r")])

        rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1].name]
        rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1].name]

        push!(step_data, (i, j, rf, rr))
    end

    n = length(forms)
    function rhs!(du, u, p, t)
        fill!(du, 0.0)
        for (i, j, rf, rr) in step_data
            flux = rf * u[i] - rr * u[j]
            du[i] -= flux
            du[j] += flux
        end
    end
    return rhs!
end

function ode_steady_state_flux(m::EnzymeMechanism, params, concs; E_total=1.0)
    forms = enzyme_forms(m)
    n = length(forms)

    u0 = zeros(n)
    u0[1] = E_total

    rhs! = build_ode_rhs(m, params, concs)
    prob = ODEProblem(rhs!, u0, (0.0, 1e6))
    sol = solve(prob, RadauIIA9(); abstol=1e-12, reltol=1e-12)
    u_ss = sol.u[end]

    # Flux through step 1
    lhs, rhs = m.steps[1]
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))
    e_lhs = [s for s in lhs if s.role == enzyme][1]
    e_rhs = [s for s in rhs if s.role == enzyme][1]
    i = name_to_idx[e_lhs.name]
    j = name_to_idx[e_rhs.name]

    m_lhs = [s for s in lhs if s.role == metabolite]
    m_rhs = [s for s in rhs if s.role == metabolite]

    step_idx = 1
    kf = Float64(params[Symbol("k$(step_idx)f")])
    kr = Float64(params[Symbol("k$(step_idx)r")])
    rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1].name]
    rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1].name]

    return rf * u_ss[i] - rr * u_ss[j]
end

@testset "ODE steady-state validation" begin
    @testset "Uni-Uni" begin
        E = Species(:E, enzyme)
        ES = Species(:ES, enzyme)
        S = Species(:S, metabolite, Dict(:C => 1))
        P = Species(:P, metabolite, Dict(:C => 1))

        m = EnzymeMechanism([[E, S] => [ES], [ES] => [E, P]])
        fn = rate_function(m)

        rng = Random.MersenneTwister(2001)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S, :P]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = fn(params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Uni-Bi" begin
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
        fn = rate_function(m)

        rng = Random.MersenneTwister(3001)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :P1, :P2]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = fn(params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Ping-Pong Bi-Bi" begin
        E  = Species(:E,  enzyme)
        EA = Species(:EA, enzyme)
        FP = Species(:FP, enzyme)
        F  = Species(:F,  enzyme)
        FB = Species(:FB, enzyme)
        EQ = Species(:EQ, enzyme)
        A     = Species(:A, metabolite, Dict(:C => 2, :N => 1))
        P_met = Species(:P, metabolite, Dict(:C => 2))
        B     = Species(:B, metabolite, Dict(:C => 3))
        Q     = Species(:Q, metabolite, Dict(:C => 3, :N => 1))

        m = EnzymeMechanism([
            [E, A] => [EA], [EA] => [FP], [FP] => [F, P_met],
            [F, B] => [FB], [FB] => [EQ], [EQ] => [E, Q]
        ])
        fn = rate_function(m)

        rng = Random.MersenneTwister(3002)
        for _ in 1:10
            params, concs = random_params_concs(m, [:A, :P, :B, :Q]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = fn(params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Bi-Uni" begin
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
        fn = rate_function(m)

        rng = Random.MersenneTwister(3003)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = fn(params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Bi-Bi" begin
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
        fn = rate_function(m)

        rng = Random.MersenneTwister(3004)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = fn(params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Bi-Ter" begin
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
        fn = rate_function(m)

        rng = Random.MersenneTwister(3005)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2, :P3]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = fn(params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Ter-Bi" begin
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
        fn = rate_function(m)

        rng = Random.MersenneTwister(3006)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = fn(params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Ter-Ter" begin
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
        fn = rate_function(m)

        rng = Random.MersenneTwister(3007)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2, :P3]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = fn(params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Random-order Bi-Bi (branched)" begin
        E   = Species(:E,   enzyme)
        EA  = Species(:EA,  enzyme)
        EB  = Species(:EB,  enzyme)
        EAB = Species(:EAB, enzyme)
        EPQ = Species(:EPQ, enzyme)
        EQ  = Species(:EQ,  enzyme)
        A   = Species(:A, metabolite, Dict(:C => 1))
        B   = Species(:B, metabolite, Dict(:N => 1))
        P   = Species(:P, metabolite, Dict(:C => 1))
        Q   = Species(:Q, metabolite, Dict(:N => 1))

        m = EnzymeMechanism([
            [E, A]   => [EA],
            [E, B]   => [EB],
            [EA, B]  => [EAB],
            [EB, A]  => [EAB],
            [EAB]    => [EPQ],
            [EPQ]    => [EQ, P],
            [EQ]     => [E, Q]
        ])
        fn = rate_function(m)

        rng = Random.MersenneTwister(3008)
        for _ in 1:10
            params, concs = random_params_concs(m, [:A, :B, :P, :Q]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = fn(params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end
end
