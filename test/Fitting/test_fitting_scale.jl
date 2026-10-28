# ABOUTME: Tests for the reference scale of relative fits: the reference anchors a problem
# ABOUTME: takes, the fit they scale, and rescale_parameter_values.
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
                                 k_E_S_to_ES = 10.0 / Km, Keq = 2.0)
    mm_concs = [(S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
                (S = 5.0, P = 0.1), (S = 10.0, P = 0.1), (S = 1.0, P = 0.5)]
    mm_data = (group = fill("G1", length(mm_concs)),
               Rate = [rate_equation(mm, c, merge(_testhelper_mm_params(1.0),
                                                  (E_total = 1.0,))) for c in mm_concs],
               S = [c.S for c in mm_concs], P = [c.P for c in mm_concs])

    @testset "reference is required in relative mode" begin
        required = "relative mode needs reference = [(concentrations...) => turnover, " *
                   "...] to set the rate scale; pass absolute = true for absolute " *
                   "per-enzyme turnover data"
        @test_throws ErrorException(required) FittingProblem(mm, mm_data; Keq = 2.0)
        # Absolute data set the scale themselves: no reference is needed, and one passed
        # is an error.
        fp_abs = FittingProblem(mm, mm_data; Keq = 2.0, absolute = true)
        @test fp_abs.absolute && isempty(fp_abs.reference)
        @test_throws(ErrorException("reference applies only to relative mode; absolute " *
                                    "turnover data set the rate scale themselves"),
                     FittingProblem(mm, mm_data; Keq = 2.0, absolute = true,
                                    reference = [(S = 10.0,) => 3.0]))
        # The problem stores the anchors with Float64 turnovers; integers convert.
        fp = FittingProblem(mm, mm_data; Keq = 2.0, reference = [(S = 10,) => 3])
        @test fp.reference == [(S = 10,) => 3.0]
        @test fp.reference isa Vector{Pair{NamedTuple, Float64}}
        # The concrete-mechanism constructor passes the reference on.
        fp_conc = FittingProblem(EnzymeRates.Mechanism(mm), mm_data; Keq = 2.0,
                                 reference = [(S = 10.0,) => 3.0])
        @test fp_conc.reference == [(S = 10.0,) => 3.0]
    end

    @testset "reference anchors are validated at construction" begin
        _testhelper_problem(reference) =
            FittingProblem(mm, mm_data; Keq = 2.0, reference)
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
        fp = FittingProblem(mm, mm_data; Keq = 2.0, reference = [(S = 10.0,) => 3.0])
        fit = fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts = 1,
                                maxiters = 50)
        @test keys(fit) == (:params, :loss, :retcode)
        @test isfinite(fit.loss)
        @test rate_equation(mm, (S = 10.0, P = 0.0),
                            merge(fit.params, (Keq = 2.0, E_total = 1.0))) ≈ 3.0 rtol = 1e-8
        # Several anchors combine by the geometric mean of measured over predicted
        # turnover, as rescale_parameter_values combines them.
        anchors = [(S = 10.0,) => 3.0, (S = 1.0,) => 2.0]
        fit2 = fit_rate_equation(
            FittingProblem(mm, mm_data; Keq = 2.0, reference = anchors),
            CMAEvolutionStrategyOpt(); n_restarts = 1, maxiters = 50)
        q, _ = rescale_parameter_values(mm, merge(fit2.params, (Keq = 2.0,)), anchors)
        @test all(k -> isapprox(fit2.params[k], q[k]; rtol = 1e-12), keys(fit2.params))
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
        fp = FittingProblem(bibi, data; Keq = 2.0, reference = [(A = 1.0,) => 1.0])
        fit = fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts = 1,
                                maxiters = 20)
        @test fit.loss == Inf
        @test fit.retcode === :ReferenceFailed
        @test keys(fit.params) == EnzymeRates.fitted_params(bibi)
        @test all(v -> isfinite(v) && v > 0, values(fit.params))
    end

    @testset "reference scaling leaves L and binding constants" begin
        # MWC uni-uni, catalysis in the active conformation only; X binds the inactive
        # conformation, Y both with different constants.
        two_site = @allosteric_mechanism begin
            substrates: S
            products: P
            allosteric_regulators: X, Y
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: OnlyA
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): X::OnlyI
            regulatory_site(1): Y::NonequalAI
        end
        p = _testhelper_random_independent_params(two_site;
                                                  rng = Random.MersenneTwister(3))
        concs = (S = 2.0, P = 0.5, X = 0.3, Y = 1.5)
        predicted = rate_equation(two_site, concs, merge(p, (E_total = 1.0,)))
        q, spread = rescale_parameter_values(two_site, p, [concs => 7.0])
        @test spread == 0.0
        @test keys(q) == keys(p)
        rate_constants = [k for k in keys(p) if startswith(String(k), "k_")]
        @test !isempty(rate_constants)
        @test :L in keys(p)
        # Every k_… constant changes by one factor, measured over predicted turnover.
        for k in rate_constants
            @test q[k] / p[k] ≈ 7.0 / predicted rtol = 1e-12
        end
        # Every other entry — binding constants, L, regulator constants, Keq, E_total —
        # stays exactly as it was.
        for k in keys(p)
            k in rate_constants || @test q[k] === p[k]
        end
    end

    @testset "rescale_parameter_values: measured turnovers" begin
        p = _testhelper_mm_params(1.0)
        _testhelper_mm_turnover(q, concs) = rate_equation(
            mm, merge((S = 0.0, P = 0.0), concs), merge(q, (E_total = 1.0,)))
        # One anchor: the rescaled law predicts the measured turnover there; P, unlisted,
        # is zero.
        q, spread = rescale_parameter_values(mm, p, [(S = 10.0,) => 42.0])
        @test _testhelper_mm_turnover(q, (S = 10.0,)) ≈ 42.0 rtol = 1e-12
        @test spread == 0.0
        @test q.Keq == 2.0
        # Two anchors combine by the geometric mean of measured over predicted turnover
        # (predicted 50/11 at S = 10 and 5/2 at S = 1); the spread is the standard
        # deviation of the two log ratios.
        r1, r2 = 42.0 / (50 / 11), 10.0 / (5 / 2)
        q2, spread2 = rescale_parameter_values(mm, p,
                                               [(S = 10.0,) => 42.0, (S = 1.0,) => 10.0])
        @test q2.k_ES_to_E_P ≈ 5.0 * sqrt(r1 * r2) rtol = 1e-12
        @test spread2 ≈ abs(log(r1) - log(r2)) / sqrt(2) rtol = 1e-12
        # Integer concentrations and turnovers work.
        @test rescale_parameter_values(mm, p, [(S = 10,) => 42])[1] == q
        # A reverse anchor: the product alone runs the reaction backwards.
        qr, _ = rescale_parameter_values(mm, p, [(P = 1.0,) => -2.0])
        @test _testhelper_mm_turnover(qr, (P = 1.0,)) ≈ -2.0 rtol = 1e-12
    end

    @testset "rescale_parameter_values: errors" begin
        p = _testhelper_mm_params(1.0)
        nokeq = (k_ES_to_E_P = 5.0, k_ES_to_E_S = 5.0, k_E_S_to_ES = 10.0)
        keq_message = "params must include Keq, which the dependent rate constants " *
                      "need; merge it into a fit's params"
        @test_throws ErrorException(keq_message) rescale_parameter_values(
            mm, nokeq, [(S = 1.0,) => 1.0])
        # Without anchors there is no scale to set.
        @test_throws MethodError rescale_parameter_values(mm, p)
        @test_throws(ErrorException("kcats is empty; pass at least one " *
                                    "(concentrations => turnover) anchor"),
                     rescale_parameter_values(mm, p, Pair{NamedTuple, Float64}[]))
        @test_throws(ErrorException("kcats: concentrations must be a NamedTuple; got " *
                                    "(1.0,)"),
                     rescale_parameter_values(mm, p, [(1.0,) => 1.0]))
        @test_throws(ErrorException("kcats: X is not a metabolite of the mechanism " *
                                    "(metabolites: S, P)"),
                     rescale_parameter_values(mm, p, [(X = 1.0,) => 1.0]))
        @test_throws(ErrorException("kcats: concentration of S must be a finite number " *
                                    "≥ 0; got -1.0"),
                     rescale_parameter_values(mm, p, [(S = -1.0,) => 1.0]))
        for (bad, shown) in ((0.0, "0.0"), (NaN, "NaN"), (Inf, "Inf"))
            @test_throws(ErrorException("kcats: turnover must be a finite, nonzero " *
                                        "number; got $shown"),
                         rescale_parameter_values(mm, p, [(S = 1.0,) => bad]))
        end
        anchor_message = "an anchor needs a finite, nonzero predicted turnover of the " *
                         "measured turnover's sign"
        # No substrate and no product: the predicted turnover is zero. The message names
        # the failing anchor and its prediction.
        @test_throws(ErrorException("kcats: at (S = 0.0,) the mechanism predicts " *
                                    "turnover 0.0; " * anchor_message),
                     rescale_parameter_values(mm, p,
                                              [(S = 1.0,) => 1.0, (S = 0.0,) => 1.0]))
        @test_throws(anchor_message, rescale_parameter_values(mm, p, [(S = 0.0,) => 1.0]))
        # The product alone runs the reaction backwards, against a forward measurement.
        @test_throws(anchor_message, rescale_parameter_values(mm, p, [(P = 1.0,) => 1.0]))
        # These constants overflow the predicted turnover at S = 1e300 to Inf, which no
        # finite scale reproduces.
        overflow = (k_ES_to_E_P = 1e308, k_ES_to_E_S = 1.0, k_E_S_to_ES = 1.0, Keq = 2.0)
        @test_throws(anchor_message,
                     rescale_parameter_values(mm, overflow, [(S = 1e300,) => 1.0]))
        # Every product of two constants underflows to zero, and so does the predicted
        # turnover.
        tiny = (k_ES_to_E_P = 1e-200, k_ES_to_E_S = 1e-200, k_E_S_to_ES = 1e-200,
                Keq = 2.0)
        @test_throws(anchor_message,
                     rescale_parameter_values(mm, tiny, [(S = 1.0,) => 1.0]))
    end

end
