using OrdinaryDiffEqFIRK

function build_ode_rhs(m, params, concs)
    forms = enzyme_forms(m)
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))

    # Precompute pseudo-first-order rates for each step
    step_data = []
    for (step_idx, (lhs, rhs)) in enumerate(steps(m))
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

function ode_steady_state_flux(m, params, concs; E_total=1.0)
    forms = enzyme_forms(m)
    n = length(forms)
    ref_name, nu_ref = _reference_metabolite(m)

    u0 = zeros(n)
    u0[1] = E_total

    rhs! = build_ode_rhs(m, params, concs)
    prob = ODEProblem(rhs!, u0, (0.0, 1e6))
    sol = solve(prob, RadauIIA9(); abstol=1e-12, reltol=1e-12)
    u_ss = sol.u[end]

    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))
    v = 0.0
    for (step_idx, (lhs, rhs)) in enumerate(steps(m))
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

        flux = rf * u_ss[i] - rr * u_ss[j]
        if !isempty(m_lhs) && m_lhs[1].name == ref_name
            v += flux
        elseif !isempty(m_rhs) && m_rhs[1].name == ref_name
            v -= flux
        end
    end

    return v / abs(nu_ref)
end

@testset "ODE steady-state validation" begin
    @testset "Uni-Uni" begin
        E = Species(:E, enzyme, Dict{Symbol,Int}())
        ES = Species(:ES, enzyme, Dict(:C => 1))
        S = Species(:S, metabolite, Dict(:C => 1))
        P = Species(:P, metabolite, Dict(:C => 1))

        steps = [[E, S] => [ES], [ES] => [E, P]]
        m = mechanism_from_species([S], [P], Species[], [E, ES], steps)
        rng = Random.MersenneTwister(2001)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S, :P]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Uni-Bi" begin
        E     = Species(:E,     enzyme, Dict{Symbol,Int}())
        ES1   = Species(:ES1,   enzyme, Dict(:C => 1, :H => 1))
        EP1P2 = Species(:EP1P2, enzyme, Dict(:C => 1, :H => 1))
        EP2   = Species(:EP2,   enzyme, Dict(:H => 1))
        S1    = Species(:S1, metabolite, Dict(:C => 1, :H => 1))
        P1    = Species(:P1, metabolite, Dict(:C => 1))
        P2    = Species(:P2, metabolite, Dict(:H => 1))

        steps = [
            [E, S1] => [ES1],
            [ES1] => [EP1P2],
            [EP1P2] => [EP2, P1],
            [EP2] => [E, P2]
        ]
        m = mechanism_from_species([S1], [P1, P2], Species[], [E, ES1, EP1P2, EP2], steps)
        rng = Random.MersenneTwister(3001)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :P1, :P2]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Ping-Pong Bi-Bi" begin
        E  = Species(:E,  enzyme, Dict{Symbol,Int}())
        EA = Species(:EA, enzyme, Dict(:C => 2, :N => 1))
        FP = Species(:FP, enzyme, Dict(:C => 2, :N => 1))
        F  = Species(:F,  enzyme, Dict(:N => 1))
        FB = Species(:FB, enzyme, Dict(:C => 3, :N => 1))
        EQ = Species(:EQ, enzyme, Dict(:C => 3, :N => 1))
        A     = Species(:A, metabolite, Dict(:C => 2, :N => 1))
        P_met = Species(:P, metabolite, Dict(:C => 2))
        B     = Species(:B, metabolite, Dict(:C => 3))
        Q     = Species(:Q, metabolite, Dict(:C => 3, :N => 1))

        steps = [
            [E, A] => [EA], [EA] => [FP], [FP] => [F, P_met],
            [F, B] => [FB], [FB] => [EQ], [EQ] => [E, Q]
        ]
        m = mechanism_from_species([A, B], [P_met, Q], Species[], [E, EA, FP, F, FB, EQ], steps)
        rng = Random.MersenneTwister(3002)
        for _ in 1:10
            params, concs = random_params_concs(m, [:A, :P, :B, :Q]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Bi-Uni" begin
        E      = Species(:E,      enzyme, Dict{Symbol,Int}())
        ES1    = Species(:ES1,    enzyme, Dict(:C => 1))
        ES1S2  = Species(:ES1S2,  enzyme, Dict(:C => 1, :H => 1))
        EP1    = Species(:EP1,    enzyme, Dict(:C => 1, :H => 1))
        S1     = Species(:S1, metabolite, Dict(:C => 1))
        S2     = Species(:S2, metabolite, Dict(:H => 1))
        P1     = Species(:P1, metabolite, Dict(:C => 1, :H => 1))

        steps = [
            [E, S1] => [ES1],
            [ES1, S2] => [ES1S2],
            [ES1S2] => [EP1],
            [EP1] => [E, P1]
        ]
        m = mechanism_from_species([S1, S2], [P1], Species[], [E, ES1, ES1S2, EP1], steps)
        rng = Random.MersenneTwister(3003)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Bi-Bi" begin
        E      = Species(:E,      enzyme, Dict{Symbol,Int}())
        ES1    = Species(:ES1,    enzyme, Dict(:C => 1))
        ES1S2  = Species(:ES1S2,  enzyme, Dict(:C => 1, :H => 1))
        EP1P2  = Species(:EP1P2,  enzyme, Dict(:C => 1, :H => 1))
        EP2    = Species(:EP2,    enzyme, Dict(:H => 1))
        S1     = Species(:S1, metabolite, Dict(:C => 1))
        S2     = Species(:S2, metabolite, Dict(:H => 1))
        P1     = Species(:P1, metabolite, Dict(:C => 1))
        P2     = Species(:P2, metabolite, Dict(:H => 1))

        steps = [
            [E, S1] => [ES1],
            [ES1, S2] => [ES1S2],
            [ES1S2] => [EP1P2],
            [EP1P2] => [EP2, P1],
            [EP2] => [E, P2]
        ]
        m = mechanism_from_species([S1, S2], [P1, P2], Species[], [E, ES1, ES1S2, EP1P2, EP2], steps)
        rng = Random.MersenneTwister(3004)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Bi-Ter" begin
        E        = Species(:E,        enzyme, Dict{Symbol,Int}())
        ES1      = Species(:ES1,      enzyme, Dict(:C => 1))
        ES1S2    = Species(:ES1S2,    enzyme, Dict(:C => 1, :H => 1, :N => 1))
        EP1P2P3  = Species(:EP1P2P3,  enzyme, Dict(:C => 1, :H => 1, :N => 1))
        EP2P3    = Species(:EP2P3,    enzyme, Dict(:H => 1, :N => 1))
        EP3      = Species(:EP3,      enzyme, Dict(:N => 1))
        S1       = Species(:S1, metabolite, Dict(:C => 1))
        S2       = Species(:S2, metabolite, Dict(:H => 1, :N => 1))
        P1       = Species(:P1, metabolite, Dict(:C => 1))
        P2       = Species(:P2, metabolite, Dict(:H => 1))
        P3       = Species(:P3, metabolite, Dict(:N => 1))

        steps = [
            [E, S1] => [ES1],
            [ES1, S2] => [ES1S2],
            [ES1S2] => [EP1P2P3],
            [EP1P2P3] => [EP2P3, P1],
            [EP2P3] => [EP3, P2],
            [EP3] => [E, P3]
        ]
        m = mechanism_from_species([S1, S2], [P1, P2, P3], Species[], [E, ES1, ES1S2, EP1P2P3, EP2P3, EP3], steps)
        rng = Random.MersenneTwister(3005)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2, :P3]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Ter-Bi" begin
        E        = Species(:E,        enzyme, Dict{Symbol,Int}())
        ES1      = Species(:ES1,      enzyme, Dict(:C => 1))
        ES1S2    = Species(:ES1S2,    enzyme, Dict(:C => 1, :H => 1))
        ES1S2S3  = Species(:ES1S2S3,  enzyme, Dict(:C => 1, :H => 1, :N => 1))
        EP1P2    = Species(:EP1P2,    enzyme, Dict(:C => 1, :H => 1, :N => 1))
        EP2      = Species(:EP2,      enzyme, Dict(:N => 1))
        S1       = Species(:S1, metabolite, Dict(:C => 1))
        S2       = Species(:S2, metabolite, Dict(:H => 1))
        S3       = Species(:S3, metabolite, Dict(:N => 1))
        P1       = Species(:P1, metabolite, Dict(:C => 1, :H => 1))
        P2       = Species(:P2, metabolite, Dict(:N => 1))

        steps = [
            [E, S1] => [ES1],
            [ES1, S2] => [ES1S2],
            [ES1S2, S3] => [ES1S2S3],
            [ES1S2S3] => [EP1P2],
            [EP1P2] => [EP2, P1],
            [EP2] => [E, P2]
        ]
        m = mechanism_from_species([S1, S2, S3], [P1, P2], Species[], [E, ES1, ES1S2, ES1S2S3, EP1P2, EP2], steps)
        rng = Random.MersenneTwister(3006)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Ter-Ter" begin
        E          = Species(:E,          enzyme, Dict{Symbol,Int}())
        ES1        = Species(:ES1,        enzyme, Dict(:C => 1))
        ES1S2      = Species(:ES1S2,      enzyme, Dict(:C => 1, :H => 1))
        ES1S2S3    = Species(:ES1S2S3,    enzyme, Dict(:C => 1, :H => 1, :N => 1))
        EP1P2P3    = Species(:EP1P2P3,    enzyme, Dict(:C => 1, :H => 1, :N => 1))
        EP2P3      = Species(:EP2P3,      enzyme, Dict(:H => 1, :N => 1))
        EP3        = Species(:EP3,        enzyme, Dict(:N => 1))
        S1         = Species(:S1, metabolite, Dict(:C => 1))
        S2         = Species(:S2, metabolite, Dict(:H => 1))
        S3         = Species(:S3, metabolite, Dict(:N => 1))
        P1         = Species(:P1, metabolite, Dict(:C => 1))
        P2         = Species(:P2, metabolite, Dict(:H => 1))
        P3         = Species(:P3, metabolite, Dict(:N => 1))

        steps = [
            [E, S1] => [ES1],
            [ES1, S2] => [ES1S2],
            [ES1S2, S3] => [ES1S2S3],
            [ES1S2S3] => [EP1P2P3],
            [EP1P2P3] => [EP2P3, P1],
            [EP2P3] => [EP3, P2],
            [EP3] => [E, P3]
        ]
        m = mechanism_from_species([S1, S2, S3], [P1, P2, P3], Species[], [E, ES1, ES1S2, ES1S2S3, EP1P2P3, EP2P3, EP3], steps)
        rng = Random.MersenneTwister(3007)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2, :P3]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Random-order Bi-Bi (branched)" begin
        E   = Species(:E,   enzyme, Dict{Symbol,Int}())
        EA  = Species(:EA,  enzyme, Dict(:C => 1))
        EB  = Species(:EB,  enzyme, Dict(:N => 1))
        EAB = Species(:EAB, enzyme, Dict(:C => 1, :N => 1))
        EPQ = Species(:EPQ, enzyme, Dict(:C => 1, :N => 1))
        EQ  = Species(:EQ,  enzyme, Dict(:N => 1))
        A   = Species(:A, metabolite, Dict(:C => 1))
        B   = Species(:B, metabolite, Dict(:N => 1))
        P   = Species(:P, metabolite, Dict(:C => 1))
        Q   = Species(:Q, metabolite, Dict(:N => 1))

        steps = [
            [E, A]   => [EA],
            [E, B]   => [EB],
            [EA, B]  => [EAB],
            [EB, A]  => [EAB],
            [EAB]    => [EPQ],
            [EPQ]    => [EQ, P],
            [EQ]     => [E, Q]
        ]
        m = mechanism_from_species([A, B], [P, Q], Species[], [E, EA, EB, EAB, EPQ, EQ], steps)
        rng = Random.MersenneTwister(3008)
        for _ in 1:10
            params, concs = random_params_concs(m, [:A, :B, :P, :Q]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end
end
