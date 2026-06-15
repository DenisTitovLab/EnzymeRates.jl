# ABOUTME: Tests for fitting rate equations to data via FittingProblem and
# ABOUTME: fit_rate_equation, including loss evaluation and parameter recovery.
using Tables

@testset "Fitting" begin

    # ── Helper: build a Uni-Uni mechanism ─────────────────────────────────────
    uni_uni = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end

    # ── Synthetic data generator ──────────────────────────────────────────────
    function make_synthetic_data(
            mechanism, true_params, concs_list;
            groups=fill("G1", length(concs_list)),
            scale=1.0,
    )
        rates = Float64[]
        for (i, concs) in enumerate(concs_list)
            r = rate_equation(mechanism, concs, true_params) * scale
            push!(rates, r)
        end
        met_names = metabolites(mechanism)
        cols = Dict{Symbol, Vector}()
        cols[:group] = groups
        cols[:Rate] = rates
        for mn in met_names
            cols[mn] = [concs[mn] for concs in concs_list]
        end
        return (; (k => cols[k] for k in (:group, :Rate, met_names...))...)
    end

    # ── Test 1: Mechanism-level accessors ─────────────────────────────────────
    @testset "Mechanism-level accessors" begin
        all_param_syms = parameters(uni_uni)
        expected_fitted = Tuple(p for p in all_param_syms if p !== :E_total && p !== :Keq)

        @test EnzymeRates.fitted_params(uni_uni) == expected_fitted
        @test metabolites(uni_uni) == (:S, :P)
    end

    # ── Test 2: FittingProblem construction ───────────────────────────────────
    @testset "Construction" begin
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [
            (S = 1.0, P = 0.1),
            (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1),
            (S = 1.0, P = 0.5),
            (S = 2.0, P = 0.5),
        ]
        data = make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        @test length(fp.log_abs_rates) == 5
        @test length(fp.group_point_indexes) == 1
        @test length(fp.group_point_indexes[1]) == 5
    end

    # ── Test 3: Loss function correctness ─────────────────────────────────────
    @testset "Loss at true params is zero" begin
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [
            (S = 1.0, P = 0.1),
            (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1),
            (S = 1.0, P = 0.5),
            (S = 2.0, P = 0.5),
        ]
        data = make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        pn = EnzymeRates.fitted_params(uni_uni)
        x_true = [log(true_params[p]) for p in pn]
        l = EnzymeRates.loss!(x_true, fp)
        @test l ≈ 0.0 atol=1e-20
    end

    # ── Test 4: Per-group centering invariance ───────────────────────────────
    @testset "Centering invariance" begin
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [
            (S = 1.0, P = 0.1),
            (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1),
            (S = 1.0, P = 0.5),
            (S = 2.0, P = 0.5),
        ]

        # Data with scale=1
        data1 = make_synthetic_data(uni_uni, true_params, concs_list;
            groups=fill("G1", 5), scale=1.0)
        fp1 = FittingProblem(uni_uni, data1; Keq=Keq_val)

        # Data with scale=10 (simulates different E_total)
        data2 = make_synthetic_data(uni_uni, true_params, concs_list;
            groups=fill("G1", 5), scale=10.0)
        fp2 = FittingProblem(uni_uni, data2; Keq=Keq_val)

        # For any x, loss should be the same (centering removes the uniform scale)
        np = length(EnzymeRates.fitted_params(uni_uni))
        @test all(1:10) do _
            x = randn(np) .* 2.0
            isapprox(EnzymeRates.loss!(x, fp1), EnzymeRates.loss!(x, fp2); rtol=1e-12)
        end
    end

    # ── Test 5: Multi-group centering invariance ─────────────────────────────
    @testset "Multi-group centering invariance" begin
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [
            (S = 1.0, P = 0.1),
            (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1),
            (S = 1.0, P = 0.5),
            (S = 2.0, P = 0.5),
        ]

        # Two groups, each independently scaled
        data1 = make_synthetic_data(uni_uni, true_params, concs_list;
            groups=["G1","G1","G1","G2","G2"],
            scale=1.0)
        fp1 = FittingProblem(uni_uni, data1; Keq=Keq_val)

        # Scale group1 by 5x and group2 by 100x
        rates2 = copy(data1.Rate)
        rates2[1:3] .*= 5.0
        rates2[4:5] .*= 100.0
        data2 = merge(data1, (Rate = rates2,))
        fp2 = FittingProblem(uni_uni, data2; Keq=Keq_val)

        np = length(EnzymeRates.fitted_params(uni_uni))
        @test all(1:10) do _
            x = randn(np) .* 2.0
            isapprox(EnzymeRates.loss!(x, fp1), EnzymeRates.loss!(x, fp2); rtol=1e-12)
        end
    end

    # ── Absolute mode: uncentered loss (scale_k_to_kcat=nothing) ──────────────
    @testset "Absolute mode uncentered loss" begin
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0, Keq = Keq_val, E_total = 1.0)
        concs_list = [
            (S = 1.0, P = 0.1), (S = 2.0, P = 0.1), (S = 5.0, P = 0.1),
            (S = 1.0, P = 0.5), (S = 2.0, P = 0.5),
        ]
        data = make_synthetic_data(uni_uni, true_params, concs_list)
        pn = EnzymeRates.fitted_params(uni_uni)
        x_true = [log(true_params[p]) for p in pn]

        fp_rel = FittingProblem(uni_uni, data; Keq=Keq_val)                        # default 1.0
        fp_abs = FittingProblem(uni_uni, data; Keq=Keq_val, scale_k_to_kcat=nothing)

        # At true params both modes are ~0 (predictions match data exactly).
        @test EnzymeRates.loss!(x_true, fp_rel) ≈ 0.0 atol=1e-20
        @test EnzymeRates.loss!(x_true, fp_abs) ≈ 0.0 atol=1e-20

        # Scale every rate by 3 (a pure per-group offset). Relative loss is
        # invariant (centering removes it); absolute loss sees it: every
        # residual becomes log(3), so absolute loss = log(3)^2.
        data3 = merge(data, (Rate = data.Rate .* 3.0,))
        fp_rel3 = FittingProblem(uni_uni, data3; Keq=Keq_val)
        fp_abs3 = FittingProblem(uni_uni, data3; Keq=Keq_val, scale_k_to_kcat=nothing)
        @test EnzymeRates.loss!(x_true, fp_rel3) ≈ 0.0 atol=1e-20
        @test EnzymeRates.loss!(x_true, fp_abs3) ≈ log(3.0)^2 rtol=1e-8
    end

    # ── scale_k_to_kcat validation ────────────────────────────────────────────
    @testset "scale_k_to_kcat validation" begin
        ok_data = (group = ["G1"], Rate = [1.0], S = [1.0], P = [0.1])
        @test_throws ErrorException FittingProblem(uni_uni, ok_data; Keq=1.0, scale_k_to_kcat=0.0)
        @test_throws ErrorException FittingProblem(uni_uni, ok_data; Keq=1.0, scale_k_to_kcat=-5.0)
        @test FittingProblem(uni_uni, ok_data; Keq=1.0, scale_k_to_kcat=nothing) isa FittingProblem
        @test FittingProblem(uni_uni, ok_data; Keq=1.0) isa FittingProblem  # default 1.0
    end

    # ── Test 6: Sign-mismatch penalty ─────────────────────────────────────────
    # Regression test for all-mismatch groups: when every prediction in a
    # group is a sign mismatch, centering must not zero every deviation.
    # The loss must be nonzero to distinguish a bad mechanism from a perfect one.
    # Three sub-cases exercise each path that sets buf[i] = 10.0:
    #   (i)  pred == 0.0            (S=0, P=0)
    #   (ii) pred < 0, Rate > 0     (S=0, P>0: only reverse term survives → pred always negative)
    #   (iii) pred > 0, Rate < 0    (S>0, P=0: only forward term survives → pred always positive)
    # In all cases every point in the group is a mismatch so the expected
    # loss is: (0 from centering + 100.0 × n_mismatch) / n_data = 100.0
    @testset "All-mismatch group not cancelled by centering" begin
        Keq_val = 2.0
        pn = EnzymeRates.fitted_params(uni_uni)
        x = randn(length(pn))

        @testset "zero prediction (S=0, P=0)" begin
            data = (
                group = fill("G1", 5),
                Rate  = [1.0, 2.0, 3.0, 4.0, 5.0],
                S     = zeros(5),
                P     = zeros(5),
            )
            fp = FittingProblem(uni_uni, data; Keq=Keq_val)
            l = EnzymeRates.loss!(x, fp)
            @test l > 0.0
            @test l ≈ 100.0
        end

        @testset "sign mismatch: pred<0, Rate>0 (S=0, P>0)" begin
            # With S=0 the forward numerator term vanishes; only the reverse
            # term (∝ -P) remains, so pred < 0 for any positive parameters.
            data = (
                group = fill("G1", 5),
                Rate  = [1.0, 2.0, 3.0, 4.0, 5.0],
                S     = zeros(5),
                P     = [0.1, 0.2, 0.3, 0.4, 0.5],
            )
            fp = FittingProblem(uni_uni, data; Keq=Keq_val)
            l = EnzymeRates.loss!(x, fp)
            @test l > 0.0
            @test l ≈ 100.0
        end

        @testset "sign mismatch: pred>0, Rate<0 (S>0, P=0)" begin
            # With P=0 the reverse numerator term vanishes; only the forward
            # term (∝ S) remains, so pred > 0 for any positive parameters.
            data = (
                group = fill("G1", 5),
                Rate  = [-1.0, -2.0, -3.0, -4.0, -5.0],
                S     = [0.1, 0.2, 0.3, 0.4, 0.5],
                P     = zeros(5),
            )
            fp = FittingProblem(uni_uni, data; Keq=Keq_val)
            l = EnzymeRates.loss!(x, fp)
            @test l > 0.0
            @test l ≈ 100.0
        end
    end

    # ── Test 7: Zero allocations ──────────────────────────────────────────────
    @testset "Zero allocations" begin
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [(S = Float64(i), P = 0.1) for i in 1:20]
        data = make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        x = randn(length(EnzymeRates.fitted_params(uni_uni)))
        EnzymeRates.loss!(x, fp)  # warmup
        allocs = @allocated EnzymeRates.loss!(x, fp)
        @test allocs == 0

        # Absolute mode (uncentered branch) is equally allocation-free.
        fp_abs = FittingProblem(uni_uni, data; Keq=Keq_val, scale_k_to_kcat=nothing)
        EnzymeRates.loss!(x, fp_abs)  # warmup
        allocs_abs = @allocated EnzymeRates.loss!(x, fp_abs)
        @test allocs_abs == 0
    end

    # ── Test 8: Speed benchmark ───────────────────────────────────────────────
    @testset "Speed" begin
        # Build a larger mechanism: Ordered Bi-Bi
        # Decomposed ordered bi-bi: the central complex EAB↔EPQ becomes an
        # explicit iso step E(A, B) <--> E(P, Q). 5 steps total (was 4 in
        # the lumped-central-complex form); fitter recovers k1-k5 params.
        ordered_bi_bi = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end

        Keq_val = 1.5
        bb_params_syms = parameters(ordered_bi_bi)
        # Generate random params
        param_vals = ntuple(i -> 1.0 + 9.0 * rand(), length(bb_params_syms))
        bb_true_params = NamedTuple{bb_params_syms}(param_vals)
        bb_true_params = merge(bb_true_params, (Keq = Keq_val, E_total = 1.0))

        # 500 synthetic datapoints
        n_points = 500
        concs_list = [(A = 0.1 + 9.9*rand(), B = 0.1 + 9.9*rand(),
                       P = 0.1 + 9.9*rand(), Q = 0.1 + 9.9*rand()) for _ in 1:n_points]
        groups = [string("G", div(i-1, 50)+1) for i in 1:n_points]

        rates = [rate_equation(ordered_bi_bi, c, bb_true_params) for c in concs_list]
        met_names_bb = metabolites(ordered_bi_bi)
        data = (
            group = groups,
            Rate = rates,
            (mn => [c[mn] for c in concs_list] for mn in met_names_bb)...
        )
        fp = FittingProblem(ordered_bi_bi, data; Keq=Keq_val)

        x = randn(length(EnzymeRates.fitted_params(ordered_bi_bi)))
        EnzymeRates.loss!(x, fp)  # warmup/compile

        # Minimum over several batches defeats the GC/scheduling inflation a
        # single mean suffers (matches test_rate_equation_performance); the
        # accumulator's finiteness check keeps the calls from being elided.
        best_us = Inf
        acc = 0.0
        for _ in 1:5
            acc = 0.0
            t = @elapsed for _ in 1:2000; acc += EnzymeRates.loss!(x, fp); end
            best_us = min(best_us, t / 2000 * 1e6)
        end
        isfinite(acc) || error("loss! produced a non-finite result")
        @test best_us < 50  # 500 datapoints; ~5 μs typical local, min-of-batches strips CI noise
    end

    # ── Test 9: scale_k_to_kcat normalization + retcode ────────────────
    @testset "scale_k_to_kcat normalization" begin
        using OptimizationBBO
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
            (S = 0.5, P = 0.5), (S = 1.0, P = 0.5), (S = 2.0, P = 0.5),
        ]
        data = make_synthetic_data(uni_uni, true_params, concs_list)

        # Default scale_k_to_kcat=1.0: returned params have kcat ≈ 1.
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)
        result = fit_rate_equation(fp, BBO_adaptive_de_rand_1_bin_radiuslimited();
            n_restarts=3, maxtime=5.0)
        full = merge(result.params, (Keq = Keq_val, E_total = 1.0))
        @test EnzymeRates._kcat_forward(uni_uni, full) ≈ 1.0 rtol=0.01
        @test result.retcode isa Symbol

        # Custom target set on the FittingProblem.
        fp42 = FittingProblem(uni_uni, data; Keq=Keq_val, scale_k_to_kcat=42.0)
        result2 = fit_rate_equation(fp42, BBO_adaptive_de_rand_1_bin_radiuslimited();
            n_restarts=3, maxtime=5.0)
        full2 = merge(result2.params, (Keq = Keq_val, E_total = 1.0))
        @test EnzymeRates._kcat_forward(uni_uni, full2) ≈ 42.0 rtol=0.01
        @test result2.retcode isa Symbol

        # scale_k_to_kcat=nothing: raw params (no rescale), retcode still present.
        fpN = FittingProblem(uni_uni, data; Keq=Keq_val, scale_k_to_kcat=nothing)
        result3 = fit_rate_equation(fpN, BBO_adaptive_de_rand_1_bin_radiuslimited();
            n_restarts=3, maxtime=5.0)
        @test haskey(result3, :params)
        @test result3.retcode isa Symbol
    end

    # ── Test: solver-option forwarding (named commons + solver_kwargs) ──
    @testset "solver kwarg forwarding" begin
        using OptimizationCMAEvolutionStrategy
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0,
            Keq = Keq_val, E_total = 1.0)
        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
        ]
        data = make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        # Default (empty) solver_kwargs runs on a solver that rejects unknown
        # options — no solver-specific option is force-injected.
        res = fit_rate_equation(fp, CMAEvolutionStrategyOpt();
            n_restarts=1, maxtime=1.0)
        @test isfinite(res.loss)        # a finite loss means a fit actually ran

        # solver_kwargs is forwarded verbatim: an option no optimizer
        # recognizes surfaces as an error (proves the bag reaches `solve`);
        # a bogus name keeps this independent of any real option's support.
        @test_throws Exception fit_rate_equation(
            fp, CMAEvolutionStrategyOpt();
            n_restarts=1, maxtime=1.0,
            solver_kwargs=(; not_a_real_solver_option=1))

        # Merge semantics: the same key in both a named common option and
        # solver_kwargs does NOT raise a duplicate-keyword error (a naive
        # double-splat would). The override value taking effect is guaranteed
        # by `merge(common, solver_kwargs)` in fit_rate_equation; here we
        # assert the call succeeds and produces a finite-loss fit.
        res2 = fit_rate_equation(
            fp, CMAEvolutionStrategyOpt();
            n_restarts=1, maxtime=60.0, solver_kwargs=(; maxtime=1.0))
        @test isfinite(res2.loss)

        # Clean break: popsize/verbose are no longer accepted named kwargs.
        @test_throws Exception fit_rate_equation(
            fp, CMAEvolutionStrategyOpt(); n_restarts=1, maxtime=1.0, popsize=200)
        @test_throws Exception fit_rate_equation(
            fp, CMAEvolutionStrategyOpt(); n_restarts=1, maxtime=1.0, verbose=-9)
    end

    # ── Test 10: Validation errors ─────────────────────────────────────────────
    @testset "Validation errors" begin
        # Missing Rate column
        data_no_rate = (group = ["G1"], S = [1.0], P = [0.1])
        @test_throws ErrorException FittingProblem(uni_uni, data_no_rate; Keq=1.0)

        # Missing metabolite column
        data_no_met = (group = ["G1"], Rate = [1.0], S = [1.0])
        @test_throws ErrorException FittingProblem(uni_uni, data_no_met; Keq=1.0)

        # Zero rate
        data_zero = (group = ["G1"], Rate = [0.0], S = [1.0], P = [0.1])
        @test_throws ErrorException FittingProblem(uni_uni, data_zero; Keq=1.0)

        # Missing group column
        data_no_grp = (Rate = [1.0], S = [1.0], P = [0.1])
        @test_throws ErrorException FittingProblem(uni_uni, data_no_grp; Keq=1.0)
    end

end
