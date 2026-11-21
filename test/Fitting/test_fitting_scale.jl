# ABOUTME: Tests for the scales of a fit: the reference anchors of a relative fit, the
# ABOUTME: constants its scale multiplies, the starting values, and test_rate_equation.
using EnzymeRates: test_rate_equation
using Random
using OptimizationCMAEvolutionStrategy

@testset "Reference turnover scale" begin

    # Uni-uni with a steady-state binding; at Keq = 2 these constants give kcat = 5 and
    # the Michaelis constant Km. Fitted: k_ES_to_E_P, k_ES_to_E_S, k_E_S_to_ES.
    mm = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end
    _testhelper_mm_params(Km) = (k_ES_to_E_P = 5.0, k_ES_to_E_S = 5.0,
                                 k_E_S_to_ES = 10.0 / Km)
    mm_concs = [(S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
                (S = 5.0, P = 0.1), (S = 10.0, P = 0.1), (S = 1.0, P = 0.5)]
    mm_data = (group = fill("G1", length(mm_concs)),
               Rate = [rate_equation(mm, c, _testhelper_mm_params(1.0), 2.0)
                       for c in mm_concs],
               S = [c.S for c in mm_concs], P = [c.P for c in mm_concs])

    @testset "relative mode accepts no reference" begin
        fp_none = FittingProblem(mm, mm_data, 2.0)
        @test !fp_none.rate_data.absolute && isempty(fp_none.rate_data.reference)
        # Absolute data set the scale themselves: no reference is needed, and one passed
        # is an error.
        fp_abs = FittingProblem(mm, mm_data, 2.0; absolute = true)
        @test fp_abs.rate_data.absolute && isempty(fp_abs.rate_data.reference)
        @test_throws(ErrorException("reference applies only to relative mode; absolute " *
                                    "turnover data set the rate scale themselves"),
                     FittingProblem(mm, mm_data, 2.0; absolute = true,
                                    reference = [(S = 10.0,) => 3.0]))
        # The problem stores the anchors with Float64 turnovers; integers convert.
        fp = FittingProblem(mm, mm_data, 2.0; reference = [(S = 10,) => 3])
        @test fp.rate_data.reference == [(S = 10,) => 3.0]
        @test fp.rate_data.reference isa Vector{Pair{NamedTuple, Float64}}
        # The concrete-mechanism constructor passes the reference on.
        fp_conc = FittingProblem(EnzymeRates.Mechanism(mm), mm_data, 2.0;
                                 reference = [(S = 10.0,) => 3.0])
        @test fp_conc.rate_data.reference == [(S = 10.0,) => 3.0]
    end

    @testset "reference anchors are validated at construction" begin
        _testhelper_problem(reference) =
            FittingProblem(mm, mm_data, 2.0; reference)
        @test_throws(ErrorException("reference must be a vector of (concentrations => " *
                                    "turnover) pairs; got (S = 10.0,) => 3.0"),
                     _testhelper_problem((S = 10.0,) => 3.0))
        @test_throws(ErrorException("reference is empty; pass at least one " *
                                    "(concentrations => turnover) anchor"),
                     _testhelper_problem(Pair{NamedTuple, Float64}[]))
        @test_throws(ErrorException("reference: concentrations must be a NamedTuple; " *
                                    "got (10.0,)"),
                     _testhelper_problem([(10.0,) => 3.0]))
        @test_throws(ErrorException("reference: X is not a metabolite of the mechanism " *
                                    "(metabolites: S, P)"),
                     _testhelper_problem([(X = 1.0,) => 3.0]))
        @test_throws(ErrorException("reference: concentration of S must be a finite " *
                                    "number ≥ 0; got -1.0"),
                     _testhelper_problem([(S = -1.0,) => 3.0]))
        for (bad, shown) in ((0.0, "0.0"), (NaN, "NaN"), (Inf, "Inf"))
            @test_throws(ErrorException("reference: turnover must be a finite, nonzero " *
                                        "number; got $shown"),
                         _testhelper_problem([(S = 10.0,) => bad]))
        end
    end

    @testset "relative fit matches the reference turnover" begin
        # The scale holds at any fit quality, so a short fit is enough.
        fp = FittingProblem(mm, mm_data, 2.0; reference = [(S = 10.0,) => 3.0])
        fit = fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts = 1,
                                maxiters = 50)
        @test keys(fit) == (:params, :loss, :retcode, :at_bound)
        @test isfinite(fit.loss)
        @test rate_equation(mm, (S = 10.0, P = 0.0), fit.params, 2.0) ≈ 3.0 rtol = 1e-8
        # Several anchors combine by the geometric mean of measured over predicted
        # turnover: after the scale, the log ratios average to zero.
        anchors = [(S = 10.0,) => 3.0, (S = 1.0,) => 2.0]
        fit2 = fit_rate_equation(
            FittingProblem(mm, mm_data, 2.0; reference = anchors),
            CMAEvolutionStrategyOpt(); n_restarts = 1, maxiters = 50)
        log_ratios = [log(t / rate_equation(mm, merge((S = 0.0, P = 0.0), c), fit2.params,
                                            2.0))
                      for (c, t) in anchors]
        @test sum(log_ratios) / length(log_ratios) ≈ 0 atol = 1e-12
    end

    @testset "an unlisted metabolite sits at zero in an anchor" begin
        # Random-order bi-bi: with only A listed, B = P = Q = 0 and the turnover is 0,
        # which no scale reproduces, so the fit fails without throwing.
        bibi = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B <--> E(A, B)
                E(B) + A <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(P) + Q
                E(P, Q) <--> E(Q) + P
                E(P) <--> E + P
                E(Q) <--> E + Q
            end
        end
        data = (group = ["G1", "G1", "G1"], Rate = [0.5, 0.8, 1.0], A = [1.0, 2.0, 4.0],
                B = [1.0, 1.0, 1.0], P = [0.1, 0.1, 0.1], Q = [0.1, 0.1, 0.1])
        fp = FittingProblem(bibi, data, 2.0; reference = [(A = 1.0,) => 1.0])
        fit = fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts = 1,
                                maxiters = 20)
        @test fit.loss == Inf
        @test fit.retcode === :ReferenceFailed
        @test keys(fit.params) == parameters(bibi)
        @test all(v -> isfinite(v) && v > 0, values(fit.params))
    end

    @testset "reference scaling leaves L and binding constants" begin
        # MWC uni-uni, catalysis in the active conformation only; X binds the inactive
        # conformation, Y both with different constants.
        two_site = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: X, Y
            oligomeric_state: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: OnlyA
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::OnlyI
            regulatory_site(1): Y::NonequalAI
        end
        p, Keq = _testhelper_random_independent_params(two_site;
                                                       rng = Random.MersenneTwister(3))
        concs = (S = 2.0, P = 0.5, X = 0.3, Y = 1.5)
        predicted = rate_equation(two_site, concs, p, Keq)
        q = EnzymeRates.Fitting._scale_rate_constants(two_site, p, 7.0 / predicted)
        @test keys(q) == keys(p)
        rate_constants = [k for k in keys(p) if startswith(String(k), "k_")]
        @test !isempty(rate_constants)
        @test :L in keys(p)
        # Every k_… constant changes by one factor, measured over predicted turnover.
        for k in rate_constants
            @test q[k] / p[k] ≈ 7.0 / predicted rtol = 1e-12
        end
        # Every other entry — binding constants, L, regulator constants — stays exactly as
        # it was.
        for k in keys(p)
            k in rate_constants || @test q[k] === p[k]
        end
    end

    @testset "default scale sets the fastest data row to 1" begin
        fp = FittingProblem(mm, mm_data, 2.0)
        fit = fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts = 1,
                                maxiters = 50)
        @test isfinite(fit.loss) && fit.retcode !== :ReferenceFailed
        peak = maximum(abs(rate_equation(mm, c, fit.params, 2.0)) for c in mm_concs)
        @test peak ≈ 1 rtol = 1e-10
    end

    @testset "default scale fails softly" begin
        zero_data = (group = ["G1", "G1", "G1"], Rate = [0.5, 0.8, 1.0],
                     S = [0.0, 0.0, 0.0], P = [0.0, 0.0, 0.0])
        fit = fit_rate_equation(FittingProblem(mm, zero_data, 2.0),
                                CMAEvolutionStrategyOpt(); n_restarts = 1, maxiters = 20)
        @test fit.retcode === :ReferenceFailed
        @test fit.loss == Inf
        @test keys(fit.params) == parameters(mm)
        @test all(v -> isfinite(v) && v > 0, values(fit.params))
    end

    @testset "test_rate_equation is the fitting loss" begin
        fit = fit_rate_equation(FittingProblem(mm, mm_data, 2.0),
                                CMAEvolutionStrategyOpt(); n_restarts = 1, maxiters = 50)
        @test keys(fit.params) == parameters(mm)
        @test test_rate_equation(mm, fit.params, mm_data, 2.0) ≈ fit.loss
        # Names beyond the mechanism's fitted parameters are ignored.
        @test test_rate_equation(mm, merge(fit.params, (Keq = 99.0, E_total = 1.0)),
                                 mm_data, 2.0) ≈ fit.loss
        # The concrete mechanism gives the same loss.
        @test test_rate_equation(EnzymeRates.Mechanism(mm), fit.params, mm_data, 2.0) ≈
              fit.loss
        # Absolute mode scores the magnitude: the true constants fit the data exactly,
        # constants ten times as fast do not.
        truth = _testhelper_mm_params(1.0)
        exact = test_rate_equation(mm, truth, mm_data, 2.0; absolute = true)
        @test exact ≈ 0 atol = 1e-20
        fast = NamedTuple{keys(truth)}(10 .* values(truth))
        @test test_rate_equation(mm, fast, mm_data, 2.0; absolute = true) > 1
        without = Base.structdiff(fit.params, (k_ES_to_E_S = 0.0,))
        @test_throws ErrorException("params is missing k_ES_to_E_S") test_rate_equation(
            mm, without, mm_data, 2.0)
    end

end

@testset "starting values are unchanged" begin
    # Each spec's log centres at the typical rate and concentrations of
    # `_testhelper_typical_table`, against the values pinned before Fitting read its
    # constants' units from `parameter_units`.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        fp = FittingProblem(m, _testhelper_typical_table(m), 2.0)
        pinned = _testhelper_typical_log_centres[spec.name]
        @test first.(pinned) == collect(parameters(m))
        @test EnzymeRates.Fitting._log_centres(fp) ≈ last.(pinned) atol = 1e-12
    end
end

@testset "a relative fit's scale multiplies exactly the rate constants" begin
    # The steady-state rate constants, named k_…, take the factor; every other fitted
    # constant (a K, a regulator's constant, L) keeps its value.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        ks = parameters(m)
        params = NamedTuple{ks}(Tuple(Float64.(eachindex(ks))))
        scaled = EnzymeRates.Fitting._scale_rate_constants(m, params, 10.0)
        @test keys(scaled) == ks
        @test [scaled[k] / params[k] for k in ks] ==
              [startswith(String(k), "k_") ? 10.0 : 1.0 for k in ks]
    end
end
