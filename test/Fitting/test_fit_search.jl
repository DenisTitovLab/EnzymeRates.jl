# ABOUTME: Tests for the fit's search around the data: the concentration and rate scales,
# ABOUTME: the parameter centres they give, the starts and box built on them, and the fits.
using Optimization
using Optimization.SciMLBase: build_solution, ReturnCode, DefaultOptimizationCache
using OptimizationCMAEvolutionStrategy
using Random

# An optimizer stub that records each problem's start and box, and returns the start.
struct _testhelper_RecordingStubOpt
    seen::Vector{NamedTuple}
end
_testhelper_RecordingStubOpt() = _testhelper_RecordingStubOpt(NamedTuple[])
Optimization.allowsbounds(::_testhelper_RecordingStubOpt) = true
function Optimization.SciMLBase.__solve(prob::Optimization.OptimizationProblem,
        opt::_testhelper_RecordingStubOpt; kwargs...)
    push!(opt.seen, (u0 = copy(prob.u0), lb = copy(prob.lb), ub = copy(prob.ub)))
    build_solution(DefaultOptimizationCache(prob.f, prob.p), opt, prob.u0,
                   prob.f(prob.u0, prob.p); retcode = ReturnCode.Success)
end

# Tests that the centres of `mechanism` move with the law: scaling each metabolite by its
# own factor, `Keq` by those factors to their stoichiometry and each fitted parameter by the
# move of its centre multiplies the rate by exactly the factor on the rate reference.
function _testhelper_check_centres_match_law(mechanism, rng)
    m = _testhelper_lift(mechanism)
    compiled = EnzymeRates.compile_mechanism(m)
    mets = Tuple(metabolites(m))
    catalytic = EnzymeRates.Mechanisms.catalytic(m)
    s = Dict(x => exp(randn(rng)) for x in mets)
    f = Dict(x => exp(2 * randn(rng)) for x in mets)
    sf = Dict(x => s[x] * f[x] for x in mets)
    r, g = exp(randn(rng)), exp(2 * randn(rng))
    keq_factor = prod(f[x] for x in catalytic.products) /
                 prod(f[x] for x in catalytic.substrates)
    centre = EnzymeRates.Fitting._log_centre
    params = _testhelper_random_independent_params(m; rng)
    scaled = merge(params, (Keq = params.Keq * keq_factor,), NamedTuple(
        EnzymeRates.Derivation.name(p, m) => params[EnzymeRates.Derivation.name(p, m)] *
            exp(centre(p, m, sf, g * r) - centre(p, m, s, r))
        for p in EnzymeRates.Derivation.parameter_info(m)))
    concs = NamedTuple{mets}(Tuple(exp(randn(rng)) for _ in mets))
    concs_scaled = NamedTuple{mets}(Tuple(concs[x] * f[x] for x in mets))
    @test rate_equation(compiled, concs_scaled, scaled) ≈
          g * rate_equation(compiled, concs, params) rtol = 1e-9
end

# A noiseless table of `mechanism` at the fitted constants `params` and `Keq`: one row, in
# one group, per point of the grid `axes`, a NamedTuple of each metabolite's concentrations.
function _testhelper_noiseless_table(mechanism, params, Keq, axes)
    mets = keys(axes)
    rows = vec(collect(Iterators.product(values(axes)...)))
    p = merge(params, (Keq = Keq, E_total = 1.0))
    (; group = fill("G1", length(rows)),
       Rate = [rate_equation(mechanism, NamedTuple{mets}(row), p) for row in rows],
       NamedTuple{mets}(Tuple(getindex.(rows, i) for i in eachindex(mets)))...)
end

@testset "Fit search around the data" begin
    log_centre = EnzymeRates.Fitting._log_centre
    log_centres = EnzymeRates.Fitting._log_centres
    concentration_scales = EnzymeRates.Fitting._concentration_scales
    rate_scale = EnzymeRates.Fitting._rate_scale
    optimizer = CMAEvolutionStrategyOpt()

    # Fitted: k_E_S_to_ES, k_ES_to_E_S, k_ES_to_E_P, all of which the data determine.
    steady_state = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end
    steady_state_params = (k_E_S_to_ES = 1e7, k_ES_to_E_S = 50.0, k_ES_to_E_P = 200.0)
    steady_state_axes = (S = [1e-6, 1e-5, 1e-4, 1e-3], P = [1e-5, 3e-4, 1e-3])

    # Fitted: K_ES_to_E_S, k_ES_to_E_P.
    rapid_binding = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end

    @testset "Facts match the law" begin
        rng = Random.MersenneTwister(97)
        @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
            _testhelper_check_centres_match_law(spec.mechanism, rng)
        end
        inline_laws = (
            two_in_steady_state = @enzyme_mechanism(begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A + B <--> E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(P, Q) <--> E(Q) + P
                    E(Q) <--> E + Q
                end
            end),
            two_in_rapid = @enzyme_mechanism(begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A + B ⇌ E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(P, Q) <--> E(Q) + P
                    E(Q) <--> E + Q
                end
            end),
            two_same_in_rapid = @enzyme_mechanism(begin
                substrates: A
                products: P
                steps: begin
                    E + A <--> F + P
                    F + A <--> E + P
                    F + A + A ⇌ F(A, A)
                end
            end),
            two_out_steady_state = @enzyme_mechanism(begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E(A) + B <--> E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(P, Q) <--> E + P + Q
                end
            end),
            modifier_copy = @enzyme_mechanism(begin
                substrates: S
                products: P
                modifiers: S
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                    E(P) + S::Mod ⇌ E(P, S::Mod)
                end
            end),
            modifier_copy_steady_state = @enzyme_mechanism(begin
                substrates: S
                products: P
                modifiers: S
                steps: begin
                    E + S <--> E(S)
                    E(S) <--> E(P)
                    E(P) <--> E + P
                    E(P) + S::Mod <--> E(P, S::Mod)
                end
            end),
        )
        @testset "$label" for (label, law) in pairs(inline_laws)
            _testhelper_check_centres_match_law(law, rng)
        end
    end

    @testset "Each parameter's centre" begin
        Derivation = EnzymeRates.Derivation
        s = Dict(:S => 10.0, :P => 100.0, :R => 1000.0)
        r = 4.0
        allosteric = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: R
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)     :: EqualAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
            regulatory_site(1): R::OnlyI
        end)
        expected = (K_ES_to_E_S = log(10.0), K_EP_to_E_P = log(100.0),
                    k_A_ES_to_EP = log(4.0), K_I_Rreg = log(1000.0), L = 0.0)
        for p in Derivation.parameter_info(allosteric)
            @test log_centre(p, allosteric, s, r) ≈ expected[Derivation.name(p, allosteric)]
        end

        # A Theorell–Chance step's rate constants each read what they take up, and its
        # equilibrium constant reads what it gives off over what it takes up.
        s = Dict(:A => 2.0, :B => 10.0, :P => 100.0, :Q => 3.0)
        tc_steady = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) + B <--> E(Q) + P
                E(Q) ⇌ E + Q
            end
        end)
        step = only(k for (k, e) in enumerate(tc_steady.steps) if !e.rapid_equilibrium)
        forward = Derivation.StepParameter(step, 1, false)
        backward = Derivation.StepParameter(step, 2, false)
        @test log_centre(forward, tc_steady, s, r) ≈ log(4.0) - log(10.0)
        @test log_centre(backward, tc_steady, s, r) ≈ log(4.0) - log(100.0)
        tc_rapid = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B ⇌ E(Q) + P
                E(Q) <--> E + Q
            end
        end)
        step = only(k for (k, e) in enumerate(tc_rapid.steps) if e.rapid_equilibrium)
        @test log_centre(Derivation.StepParameter(step, 1, false), tc_rapid, s, r) ≈
              log(100.0) - log(10.0)
    end

    @testset "Issue #97 regression" begin
        K_true = 0.2e-6
        table = _testhelper_noiseless_table(rapid_binding,
            (K_ES_to_E_S = K_true, k_ES_to_E_P = 5.0), 1000.0,
            (S = [2e-8, 1e-7, 2e-7, 5e-7, 2e-6, 1e-5], P = [1e-7, 1e-5]))
        for absolute in (false, true)
            Random.seed!(97)
            fp = FittingProblem(rapid_binding, table; Keq = 1000.0, absolute)
            fit = fit_rate_equation(fp, optimizer; n_restarts = 2, maxtime = 5.0)
            @test fit.params.K_ES_to_E_S ≈ K_true rtol = 0.01
            @test abs(log(fit.params.K_ES_to_E_S) - (-15.0)) > 0.01
        end
    end

    @testset "Unit independence" begin
        # Each law, its true constants in M and the power of the concentration unit that
        # each fitted constant carries.
        laws = (
            (steady_state, steady_state_params,
             (k_E_S_to_ES = -1, k_ES_to_E_S = 0, k_ES_to_E_P = 0)),
            (rapid_binding, (K_ES_to_E_S = 2e-5, k_ES_to_E_P = 200.0),
             (K_ES_to_E_S = 1, k_ES_to_E_P = 0)),
        )
        for (law, params, powers) in laws, absolute in (true, false)
            table = _testhelper_noiseless_table(law, params, 2.0, steady_state_axes)
            fits = map((1.0, 1e3, 1e6)) do unit
                scaled = merge(table, (S = table.S .* unit, P = table.P .* unit))
                Random.seed!(11)
                fp = FittingProblem(law, scaled; Keq = 2.0, absolute)
                fit_rate_equation(fp, optimizer; n_restarts = 2, maxtime = 10.0).params
            end
            for (fit, unit) in zip(fits, (1.0, 1e3, 1e6)), k in keys(powers)
                @test log(fit[k]) - powers[k] * log(unit) ≈ log(fits[1][k]) atol = 1e-3
            end
        end
    end

    @testset "Both modes" begin
        table = _testhelper_noiseless_table(steady_state, steady_state_params, 2.0,
                                            steady_state_axes)
        problems = (
            FittingProblem(steady_state, table; Keq = 2.0, absolute = true),
            FittingProblem(steady_state, table; Keq = 2.0,
                           reference = _testhelper_row_reference(table)),
        )
        for fp in problems
            Random.seed!(23)
            fit = fit_rate_equation(fp, optimizer; n_restarts = 3, maxtime = 10.0)
            @test fit.loss < 1e-6
            for k in keys(steady_state_params)
                @test fit.params[k] ≈ steady_state_params[k] rtol = 0.01
            end
        end
    end

    @testset "Two-metabolite law" begin
        two_in = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A + B <--> E(A, B)
                E(A, B) <--> E + P + Q
            end
        end
        params = (k_E_A_B_to_EAB = 1e10, k_EAB_to_E_A_B = 100.0, k_EAB_to_E_P_Q = 300.0)
        @test issetequal(keys(params), parameters(two_in, FittedParameters()))
        table = _testhelper_noiseless_table(two_in, params, 20.0,
            (A = [1e-5, 3e-5, 1e-4, 3e-4], B = [1e-4, 1e-3], P = [1e-4, 1e-3],
             Q = [1e-4, 1e-3]))
        Random.seed!(29)
        fit = fit_rate_equation(FittingProblem(two_in, table; Keq = 20.0, absolute = true),
                                optimizer; n_restarts = 3, maxtime = 10.0)
        @test fit.loss < 1e-6
        for k in keys(params)
            @test fit.params[k] ≈ params[k] rtol = 0.01
        end
    end

    @testset "Always-zero metabolite" begin
        @testset "the fit runs through the positive concentrations' geometric mean" begin
            table = _testhelper_noiseless_table(rapid_binding,
                (K_ES_to_E_S = 1e-4, k_ES_to_E_P = 5.0), 2.0,
                (S = [1e-5, 1e-4, 1e-3], P = [0.0]))
            fp = FittingProblem(rapid_binding, table; Keq = 2.0)
            @test concentration_scales(fp)[:P] ≈ 1e-4 rtol = 1e-12
            Random.seed!(31)
            fit = fit_rate_equation(fp, optimizer; n_restarts = 2, maxtime = 5.0)
            @test isfinite(fit.loss)
            @test fit.params.K_ES_to_E_S ≈ 1e-4 rtol = 0.01
        end

        @testset "the fallback pools the mechanism's metabolites' concentrations" begin
            two_in = @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A + B <--> E(A, B)
                    E(A, B) <--> E + P + Q
                end
            end
            # The other column, `Dose`, is no metabolite, so the fallback leaves it out.
            table = (group = fill("G1", 4), Rate = [1.0, 2.0, 3.0, 4.0],
                     A = [1.0, 0.0, 0.0, 0.0], B = [10.0, 10.0, 10.0, 10.0],
                     P = zeros(4), Q = zeros(4), Dose = fill(1e6, 4))
            s = concentration_scales(FittingProblem(two_in, table; Keq = 2.0))
            @test s[:A] == 1.0
            @test s[:B] ≈ 10.0 rtol = 1e-12
            @test s[:P] ≈ 10.0^0.8 rtol = 1e-12
            @test s[:Q] == s[:P]
        end

        @testset "a table with no positive concentration centres on 1" begin
            table = (group = fill("G1", 3), Rate = [0.5, 0.8, 1.0],
                     S = zeros(3), P = zeros(3))
            fp = FittingProblem(rapid_binding, table; Keq = 2.0)
            @test concentration_scales(fp) == Dict(:S => 1.0, :P => 1.0)
        end

        # A cross-validation fold builds its problem from its training rows, which can
        # hold a metabolite only at zero while the whole table does not.
        @testset "a training subset with a metabolite zero only in its rows" begin
            params = (K_ES_to_E_S = 1e-4, k_ES_to_E_P = 5.0)
            table = _testhelper_noiseless_table(rapid_binding, params, 2.0,
                (S = [1e-5, 1e-4, 1e-3], P = [0.0, 1e-2]))
            train = table.P .== 0.0
            @test any(>(0.0), table.P) && any(train)
            train_table = (group = table.group[train], Rate = table.Rate[train],
                           S = table.S[train], P = table.P[train])
            fp = FittingProblem(rapid_binding, train_table; Keq = 2.0)
            @test concentration_scales(fp)[:P] ≈ 1e-4 rtol = 1e-12
            whole = FittingProblem(rapid_binding, table; Keq = 2.0)
            @test concentration_scales(whole)[:P] ≈ 1e-2 rtol = 1e-12
            Random.seed!(37)
            fit = fit_rate_equation(fp, optimizer; n_restarts = 2, maxtime = 5.0)
            @test isfinite(fit.loss)
            @test fit.params.K_ES_to_E_S ≈ 1e-4 rtol = 0.01
        end
    end

    @testset "Integer columns" begin
        integers = (group = fill("G1", 4), Rate = [3, 5, -7, 9],
                    S = [1, 2, 5, 10], P = [1, 1, 0, 2])
        floats = (group = integers.group, Rate = Float64.(integers.Rate),
                  S = Float64.(integers.S), P = Float64.(integers.P))
        fp_int = FittingProblem(rapid_binding, integers; Keq = 2.0)
        fp_float = FittingProblem(rapid_binding, floats; Keq = 2.0)
        @test log_centres(fp_int) == log_centres(fp_float)
        Random.seed!(41)
        fit = fit_rate_equation(fp_int, optimizer; n_restarts = 1, maxiters = 20)
        @test isfinite(fit.loss)
    end

    @testset "Search box" begin
        table = (group = ["G1", "G1", "G2", "G2"], Rate = [2.0, 8.0, -2.0, 8.0],
                 S = [1e-9, 1e-7, 1e-9, 1e-7], P = [1e-5, 1e-5, 1e-5, 1e-5])
        fp = FittingProblem(steady_state, table; Keq = 2.0)
        @test rate_scale(fp) ≈ 4.0 rtol = 1e-12
        s = concentration_scales(fp)
        @test s[:S] ≈ 1e-8 rtol = 1e-12
        @test s[:P] ≈ 1e-5 rtol = 1e-12
        names = parameters(steady_state, FittedParameters())
        centre = log_centres(fp)
        expected = (k_E_S_to_ES = log(4.0) - log(1e-8), k_ES_to_E_S = log(4.0),
                    k_ES_to_E_P = log(4.0))
        @test centre ≈ [expected[k] for k in names]

        @testset "defaults sit 15 log units either side of the centre" begin
            stub = _testhelper_RecordingStubOpt()
            Random.seed!(43)
            fit_rate_equation(fp, stub; n_restarts = 300)
            @test length(stub.seen) == 300
            @test all(seen -> seen.lb == centre .- 15.0 && seen.ub == centre .+ 15.0,
                      stub.seen)
            @test all(seen -> all(seen.lb .<= seen.u0 .<= seen.ub), stub.seen)
            # Each restart draws anew, with a standard deviation of 2 around the centre.
            offsets = reduce(vcat, [seen.u0 .- centre for seen in stub.seen])
            @test allunique(offsets)
            @test abs(sum(offsets) / length(offsets)) < 0.3
            @test 1.8 < sqrt(sum(abs2, offsets) / length(offsets)) < 2.2
        end

        @testset "starts clamp into an explicit box, which holds absolute log values" begin
            stub = _testhelper_RecordingStubOpt()
            lb, ub = fill(-1.0, 3), fill(1.0, 3)
            Random.seed!(47)
            fit_rate_equation(fp, stub; n_restarts = 1, lb, ub)
            seen = only(stub.seen)
            @test seen.lb == lb && seen.ub == ub
            # The first constant's centre lies near 20, far above this box.
            @test seen.u0[findfirst(==(:k_E_S_to_ES), names)] == 1.0
            @test all(lb .<= seen.u0 .<= ub)
        end
    end
end
