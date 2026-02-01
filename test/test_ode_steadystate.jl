using OrdinaryDiffEqFIRK

function build_ode_rhs(m::EnzymeMechanism{Species, Reactions}, params, concs) where {Species, Reactions}
    enzs = enzyme_forms(m)
    enz_names = Tuple(e[1] for e in enzs)
    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))
    enz_set = Set(enz_names)

    step_data = []
    for (step_idx, (lhs, rhs)) in enumerate(Reactions)
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        i = name_to_idx[e_lhs]
        j = name_to_idx[e_rhs]

        m_lhs = [s for s in lhs if s ∉ enz_set]
        m_rhs = [s for s in rhs if s ∉ enz_set]

        kf = Float64(params[Symbol("k$(step_idx)f")])
        kr = Float64(params[Symbol("k$(step_idx)r")])

        rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1]]
        rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1]]

        push!(step_data, (i, j, rf, rr))
    end

    n = length(enzs)
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

function ode_steady_state_flux(m::EnzymeMechanism{Species, Reactions}, params, concs) where {Species, Reactions}
    E_total = params.E_total
    enzs = enzyme_forms(m)
    n = length(enzs)
    enz_names = Tuple(e[1] for e in enzs)
    enz_set = Set(enz_names)
    ref_name, nu_ref = _reference_metabolite(m)

    u0 = zeros(n)
    u0[1] = E_total

    rhs! = build_ode_rhs(m, params, concs)
    prob = ODEProblem(rhs!, u0, (0.0, 1e6))
    sol = solve(prob, RadauIIA9(); abstol=1e-12, reltol=1e-12)
    u_ss = sol.u[end]

    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))
    v = 0.0
    for (step_idx, (lhs, rhs)) in enumerate(Reactions)
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        i = name_to_idx[e_lhs]
        j = name_to_idx[e_rhs]

        m_lhs = [s for s in lhs if s ∉ enz_set]
        m_rhs = [s for s in rhs if s ∉ enz_set]

        kf = Float64(params[Symbol("k$(step_idx)f")])
        kr = Float64(params[Symbol("k$(step_idx)r")])
        rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1]]
        rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1]]

        flux = rf * u_ss[i] - rr * u_ss[j]
        if !isempty(m_lhs) && m_lhs[1] == ref_name
            v += flux
        elseif !isempty(m_rhs) && m_rhs[1] == ref_name
            v -= flux
        end
    end

    return v / abs(nu_ref)
end

@testset "ODE steady-state validation" begin
    @testset "Uni-Uni" begin
        species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
        m = EnzymeMechanism(species, rxns)
        rng = Random.MersenneTwister(2001)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S, :P]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Uni-Bi" begin
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
        rng = Random.MersenneTwister(3001)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :P1, :P2]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Ping-Pong Bi-Bi" begin
        species = (
            ((:A, ((:C, 2), (:N, 1))), (:B, ((:C, 3),))),
            ((:P, ((:C, 2),)), (:Q, ((:C, 3), (:N, 1)))),
            (),
            ((:E, ()), (:EA, ((:C, 2), (:N, 1))), (:FP, ((:C, 2), (:N, 1))),
             (:F, ((:N, 1),)), (:FB, ((:C, 3), (:N, 1))), (:EQ, ((:C, 3), (:N, 1)))),
        )
        rxns = (
            ((:E, :A), (:EA,)), ((:EA,), (:FP,)), ((:FP,), (:F, :P)),
            ((:F, :B), (:FB,)), ((:FB,), (:EQ,)), ((:EQ,), (:E, :Q)),
        )
        m = EnzymeMechanism(species, rxns)
        rng = Random.MersenneTwister(3002)
        for _ in 1:10
            params, concs = random_params_concs(m, [:A, :P, :B, :Q]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Bi-Uni" begin
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
        rng = Random.MersenneTwister(3003)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Bi-Bi" begin
        species = (
            ((:S1, ((:C, 1),)), (:S2, ((:H, 1),))),
            ((:P1, ((:C, 1),)), (:P2, ((:H, 1),))),
            (),
            ((:E, ()), (:ES1, ((:C, 1),)), (:ES1S2, ((:C, 1), (:H, 1))),
             (:EP1P2, ((:C, 1), (:H, 1))), (:EP2, ((:H, 1),))),
        )
        rxns = (
            ((:E, :S1), (:ES1,)),
            ((:ES1, :S2), (:ES1S2,)),
            ((:ES1S2,), (:EP1P2,)),
            ((:EP1P2,), (:EP2, :P1)),
            ((:EP2,), (:E, :P2)),
        )
        m = EnzymeMechanism(species, rxns)
        rng = Random.MersenneTwister(3004)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Bi-Ter" begin
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
        rng = Random.MersenneTwister(3005)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :P1, :P2, :P3]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Ter-Bi" begin
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
        rng = Random.MersenneTwister(3006)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end

    @testset "Seq Ter-Ter" begin
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
        rng = Random.MersenneTwister(3007)
        for _ in 1:10
            params, concs = random_params_concs(m, [:S1, :S2, :S3, :P1, :P2, :P3]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
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
        rng = Random.MersenneTwister(3008)
        for _ in 1:10
            params, concs = random_params_concs(m, [:A, :B, :P, :Q]; rng=rng)
            v_ode = ode_steady_state_flux(m, params, concs)
            v_ka = rate_equation(m, params, concs)
            @test v_ode ≈ v_ka rtol=1e-6
        end
    end
end
