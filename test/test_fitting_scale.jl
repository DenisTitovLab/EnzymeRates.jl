# ABOUTME: Tests for the peak-turnover scale of relative fits: the search for a law's
# ABOUTME: peak turnover over concentrations, and rescale_parameter_values.
using Random

const _testhelper_TurnoverFunction = EnzymeRates.FunctionWrapper{Float64,
                                                                 Tuple{Vector{Float64}}}

# The largest turnover of `m` at `p` on the grid of whole decades from 10^lo_exp to
# 10^hi_exp for every metabolite; undefined points count as -Inf.
function _testhelper_grid_peak(m, p, lo_exp, hi_exp)
    mets = metabolites(m)
    turnover_params = merge(p, (E_total = 1.0,))
    maximum(Iterators.product(ntuple(_ -> lo_exp:hi_exp, length(mets))...)) do e
        v = rate_equation(m, NamedTuple{mets}(Tuple(10.0 .^ e)), turnover_params)
        isfinite(v) ? v : -Inf
    end
end

@testset "Peak turnover scale" begin

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

    # S also binds E(S) as a dead-end inhibitor. At zero product,
    # v = k·x/(1 + x + x²·K_S/K_I) with x = S/K_S, maximal at S = √(K_S·K_I), where
    # v = k/(1 + 2√(K_S/K_I)).
    si = @enzyme_mechanism begin
        substrates: S
        products:   P
        regulators: S
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
            E(S) + S::Inh ⇌ E(S, S::Inh)
        end
    end
    si_params = (K_EP_to_E_P = 1.0, K_ESSinh_to_ES_Sinh = 4.0, K_ES_to_E_S = 1.0,
                 k_ES_to_EP = 3.0, Keq = 2.0)

    # Random-order steady-state bi-bi.
    rob = @enzyme_mechanism begin
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

    # Allosteric uni-uni whose two regulators bind sites of multiplicity 2 and 1.
    two_site = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: X::OnlyI, Y::NonequalAI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(multiplicity = 2): begin
            ligands: X
        end
        regulatory_site(multiplicity = 1): begin
            ligands: Y
        end
    end

    @testset "a saturating law peaks at the range's upper substrate corner" begin
        peak = EnzymeRates._peak_turnover(mm, _testhelper_mm_params(1.0), (1e-12, 1e3))
        # v = kcat·S/(Km + S) at S = 1e3, with the product at the range's lower bound.
        @test peak.value ≈ 5 * 1e3 / (1 + 1e3) rtol = 1e-12
        @test peak.concs.S ≈ 1e3 rtol = 1e-12
        @test peak.concs.P ≈ 1e-12 rtol = 1e-12
        @test peak.at_edge === false
        # The range bounds the search: capped at S = 10 the peak is the turnover there.
        capped = EnzymeRates._peak_turnover(mm, _testhelper_mm_params(1.0), (1e-12, 10.0))
        @test capped.value ≈ 5 * 10 / (1 + 10) rtol = 1e-12
        @test capped.concs.S ≈ 10.0 rtol = 1e-12
    end

    @testset "a law still rising at the range's upper edge is flagged" begin
        # Hyperbolic in S: v(1e3)/v(1e2) = 10·(Km + 1e2)/(Km + 1e3) exceeds 1.1 exactly
        # when Km > 1e3/89 ≈ 11.24.
        flagged, unflagged = _testhelper_mm_params(11.5), _testhelper_mm_params(10.5)
        @test EnzymeRates._peak_turnover(mm, flagged, (1e-12, 1e3)).at_edge
        @test !EnzymeRates._peak_turnover(mm, unflagged, (1e-12, 1e3)).at_edge
        # Capped at S = 10 with Km = 1: v(10)/v(1) = (10/11)/(1/2) ≈ 1.82.
        p = _testhelper_mm_params(1.0)
        @test EnzymeRates._peak_turnover(mm, p, (1e-12, 10.0)).at_edge
    end

    @testset "substrate inhibition peaks inside the range" begin
        peak = EnzymeRates._peak_turnover(si, si_params, (1e-12, 1e3))
        @test peak.value ≈ 3 / (1 + 2 * sqrt(1 / 4)) rtol = 1e-9
        @test peak.concs.S ≈ 2.0 rtol = 1e-6
        @test peak.at_edge === false
    end

    @testset "a product can raise the peak above the product-free turnover" begin
        # Every constant 1 except E(P) ⇌ E + P at 1e-3 both ways (Keq = 1): Q rebinding
        # E(P) reroutes the slow release of P through E(Q) + P, so the turnover peaks at
        # an intermediate Q, 125.62 times the turnover at saturating A and B without
        # product.
        pn = EnzymeRates.fitted_params(rob)
        p = merge(NamedTuple{pn}(Tuple(k === :k_E_P_to_EP ? 1e-3 : 1.0 for k in pn)),
                  (Keq = 1.0,))
        peak = EnzymeRates._peak_turnover(rob, p, (1e-12, 1e8))
        corner = rate_equation(rob, (A = 1e8, B = 1e8, P = 1e-12, Q = 1e-12),
                               merge(p, (E_total = 1.0,)))
        @test peak.value ≈ 0.249991171308 rtol = 1e-9
        @test peak.value / corner ≈ 125.62056353 rtol = 1e-8
        @test peak.concs.Q ≈ 14120.94 rtol = 1e-4
        @test peak.concs.A ≈ 1e8 rtol = 1e-12
        @test peak.concs.P ≈ 1e-12 rtol = 1e-12
    end

    @testset "an allosteric law can peak below saturating substrate" begin
        # The active conformation turns over fast but binds S weakly; the inactive one
        # binds S tightly and turns over slowly. At catalytic multiplicity 2, saturating
        # S drives the enzyme into the inactive conformation, so the turnover peaks near
        # S = 1 at about 5 and falls toward 0.2 at the range's upper edge.
        mwc = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        p = (K_A_ES_to_E_S = 100.0, K_EP_to_E_P = 1.0, K_I_ES_to_E_S = 0.01,
             k_A_ES_to_EP = 1000.0, k_I_ES_to_EP = 0.1, L = 1e-4, Keq = 2.0)
        peak = EnzymeRates._peak_turnover(mwc, p, (1e-12, 1e3))
        @test peak.value ≈ 5.00100010001 rtol = 1e-9
        @test 0.5 < peak.concs.S < 2
        @test peak.value > 20 * rate_equation(mwc, (S = 1e3, P = 1e-12),
                                              merge(p, (E_total = 1.0,)))
        @test peak.at_edge === false
    end

    @testset "a law that needs bound product peaks at a positive product level" begin
        # Catalysis runs only on the product-bound form E(P), so the forward rate needs P
        # and vanishes as P goes to zero.
        ep = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + P ⇌ E(P)
                E(P) + S ⇌ E(P, S)
                E(P, S) <--> E(P, P)
                E(P, P) ⇌ E(P) + P
            end
        end
        pn = EnzymeRates.fitted_params(ep)
        p = merge(NamedTuple{pn}(ntuple(_ -> 1.5, length(pn))), (Keq = 2.0,))
        peak = EnzymeRates._peak_turnover(ep, p, (1e-12, 1e3))
        @test peak.value ≈ 1.49227389554 rtol = 1e-9
        @test peak.concs.S ≈ 1e3 rtol = 1e-12
        @test 1e-3 < peak.concs.P < 1e3
    end

    @testset "the peak reaches every whole-decade grid point's turnover" begin
        for (m, n_draws) in ((rob, 5), (two_site, 3))
            pn = EnzymeRates.fitted_params(m)
            for seed in 1:n_draws
                rng = Random.MersenneTwister(seed)
                p = merge(NamedTuple{pn}(Tuple(0.2 + 5 * rand(rng) for _ in pn)),
                          (Keq = 2.0,))
                peak = EnzymeRates._peak_turnover(m, p, (1e-12, 1e3))
                @test peak.value >= _testhelper_grid_peak(m, p, -12, 3) * (1 - 1e-12)
                # The reported point reaches the reported value.
                @test rate_equation(m, peak.concs, merge(p, (E_total = 1.0,))) ==
                      peak.value
            end
        end
    end

    @testset "the peak is per enzyme: E_total in params is ignored" begin
        p = _testhelper_mm_params(1.0)
        @test EnzymeRates._peak_turnover(mm, merge(p, (E_total = 7.0,)),
                                         (1e-12, 1e3)).value ==
              EnzymeRates._peak_turnover(mm, p, (1e-12, 1e3)).value
    end

    @testset "the search is deterministic" begin
        rng = Random.MersenneTwister(7)
        pn = EnzymeRates.fitted_params(rob)
        p = merge(NamedTuple{pn}(Tuple(0.2 + 5 * rand(rng) for _ in pn)), (Keq = 2.0,))
        @test EnzymeRates._peak_turnover(rob, p, (1e-12, 1e3)) ===
              EnzymeRates._peak_turnover(rob, p, (1e-12, 1e3))
    end

    @testset "the search skips undefined points and flags a rising upper edge" begin
        # A concave bowl with its top at (-1, 2), undefined wherever u[1] > 0.
        bowl = _testhelper_TurnoverFunction(
            u -> u[1] > 0 ? NaN : -((u[1] + 1)^2 + (u[2] - 2)^2))
        u, v, at_edge = EnzymeRates._maximize_on_box(bowl, 2, -12.0, 3.0)
        @test v ≈ 0 atol = 1e-12
        @test u ≈ [-1.0, 2.0] atol = 1e-5
        @test at_edge === false
        # Rising up to the upper edge: f(3) = 3 exceeds 1.1·f(2).
        ramp = _testhelper_TurnoverFunction(u -> u[1])
        u, v, at_edge = EnzymeRates._maximize_on_box(ramp, 1, -12.0, 3.0)
        @test u == [3.0] && v == 3.0 && at_edge === true
        # A best value that is not positive never counts as rising at the edge: this ramp
        # is -1 at the upper edge and rises only 0.1% over the last decade.
        negative = _testhelper_TurnoverFunction(u -> -1.0 - 1e-3 * (3 - u[1]))
        @test EnzymeRates._maximize_on_box(negative, 1, -12.0, 3.0)[3] === false
        # A function undefined everywhere has no peak.
        nowhere = _testhelper_TurnoverFunction(u -> NaN)
        @test EnzymeRates._maximize_on_box(nowhere, 2, -12.0, 3.0)[2] == -Inf
    end

    @testset "peak_range must be a tuple (lo, hi) with 0 < lo < hi < Inf" begin
        @test EnzymeRates.DEFAULT_PEAK_RANGE === (1e-12, 1e3)
        @test EnzymeRates._peak_range((1, 1000)) === (1.0, 1000.0)
        @test EnzymeRates._peak_range((1e-12, 1e3)) === (1e-12, 1e3)
        for bad in ((0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (1e-12, Inf), (NaN, 1.0),
                    (-1.0, 1.0), 1.0, (1e-12, 1e3, 1e6), [1e-12, 1e3])
            @test_throws(ErrorException("peak_range must be a tuple (lo, hi) with " *
                                        "0 < lo < hi < Inf; got $(repr(bad))"),
                         EnzymeRates._peak_range(bad))
        end
    end

    @testset "rescale_parameter_values: peak scale" begin
        p = _testhelper_mm_params(1.0)
        q, spread = rescale_parameter_values(mm, p)
        @test spread === nothing
        @test keys(q) == keys(p)
        # Every fitted constant here is a rate constant; each is divided by the peak
        # turnover 5·1e3/(1 + 1e3).
        @test q.k_ES_to_E_P ≈ 1.001 rtol = 1e-12
        @test q.k_ES_to_E_S ≈ 1.001 rtol = 1e-12
        @test q.k_E_S_to_ES ≈ 2.002 rtol = 1e-12
        @test q.Keq == 2.0
        @test EnzymeRates._peak_turnover(mm, q, (1e-12, 1e3)).value ≈ 1 rtol = 1e-12
        # E_total passes through and does not enter the scale.
        qE, _ = rescale_parameter_values(mm, merge(p, (E_total = 3.0,)))
        @test qE.E_total == 3.0 && qE.k_ES_to_E_P == q.k_ES_to_E_P
        # peak_range sets the range: capped at S = 10 the peak turnover is 5·10/11.
        q10, _ = rescale_parameter_values(mm, p; peak_range = (1e-12, 10.0))
        @test q10.k_ES_to_E_P ≈ 1.1 rtol = 1e-12
        # Binding constants stay; only the rate constant of the substrate-inhibition law
        # scales, by 1/1.5.
        q_si, _ = rescale_parameter_values(si, si_params)
        @test q_si.k_ES_to_EP ≈ 2.0 rtol = 1e-9
        @test (q_si.K_EP_to_E_P, q_si.K_ESSinh_to_ES_Sinh, q_si.K_ES_to_E_S) ==
              (1.0, 4.0, 1.0)
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
        @test_throws ErrorException(keq_message) rescale_parameter_values(mm, nokeq)
        @test_throws ErrorException(keq_message) rescale_parameter_values(
            mm, nokeq, [(S = 1.0,) => 1.0])
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
        # No substrate and no product: the predicted turnover is zero.
        @test_throws("an anchor needs a nonzero predicted turnover of the measured " *
                     "turnover's sign",
                     rescale_parameter_values(mm, p, [(S = 0.0,) => 1.0]))
        # The product alone runs the reaction backwards, against a forward measurement.
        @test_throws("an anchor needs a nonzero predicted turnover of the measured " *
                     "turnover's sign",
                     rescale_parameter_values(mm, p, [(P = 1.0,) => 1.0]))
        # Every product of two constants underflows to zero: no peak scale exists.
        tiny = (k_ES_to_E_P = 1e-200, k_ES_to_E_S = 1e-200, k_E_S_to_ES = 1e-200,
                Keq = 2.0)
        @test_throws(ErrorException("the mechanism has no positive forward turnover " *
                                    "within peak_range (1.0e-12, 1000.0); the peak scale " *
                                    "is undefined"),
                     rescale_parameter_values(mm, tiny))
    end

end
