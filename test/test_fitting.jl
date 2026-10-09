# ABOUTME: Tests for fitting rate equations to data via FittingProblem and
# ABOUTME: fit_rate_equation, including loss evaluation and parameter recovery.
using Tables

@testset "Fitting" begin

    # ── Helper: build a Uni-Uni mechanism ─────────────────────────────────────
    # Fitted: k_E_S_to_ES, k_ES_to_E_S, k_ES_to_E_P. The true points below put
    # k_ES_to_E_P = 5 at Keq = 2, so the Haldane-derived product rebinding rate is
    # k_E_P_to_ES = k_E_S_to_ES·k_ES_to_E_P/(Keq·k_ES_to_E_S) = 1.
    uni_uni = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end

    # True point and 5-point concentration grid shared by the uni-uni testsets below.
    Keq_val = 2.0
    true_params = (k_E_S_to_ES = 10.0, k_ES_to_E_S = 25.0, k_ES_to_E_P = 5.0,
                   Keq = Keq_val, E_total = 1.0)
    concs5 = [
        (S = 1.0, P = 0.1),
        (S = 2.0, P = 0.1),
        (S = 5.0, P = 0.1),
        (S = 1.0, P = 0.5),
        (S = 2.0, P = 0.5),
    ]

    # ── Helper: an allosteric uni-uni with an inactive-state regulator R ──────
    allo = @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R::OnlyI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: EqualAI
            E(S) <--> E(P)   :: OnlyA
            E(P) ⇌ E + P     :: EqualAI
        end
    end
    allo_fitted = EnzymeRates.fitted_params(allo)
    allo_params = merge(
        NamedTuple{allo_fitted}(ntuple(i -> 1.0 + 0.1 * i, length(allo_fitted))),
        (Keq = Keq_val, E_total = 1.0))
    allo_concs = [(S = Float64(i), P = 0.1, R = 0.5) for i in 1:20]

    # ── Synthetic data generator ──────────────────────────────────────────────
    function _testhelper_make_synthetic_data(
            mechanism, true_params, concs_list;
            groups=fill("G1", length(concs_list)),
    )
        rates = Float64[]
        for (i, concs) in enumerate(concs_list)
            r = rate_equation(mechanism, concs, true_params)
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

    # ── BlackBoxOptim fit: the only warning is its convergence ────────────────
    # BlackBoxOptim can converge on these fits before maxtime: an absolute loss has a
    # single optimum, and a centered population can collapse along its scale direction.
    # It then stops because the search has converged, a reason Optimization does not
    # recognize, so Optimization warns. Runs `fit`, checks that every warning it logs
    # is that one, and returns its result.
    function _testhelper_fit_capturing_convergence(fit)
        logs, result = Test.collect_test_logs(fit)
        @test all(l -> l.level < Base.CoreLogging.Warn ||
                       occursin("probably search has converged", string(l.message)), logs)
        result
    end

    # ── Mechanism-level accessors ─────────────────────────────────────────────
    @testset "Mechanism-level accessors" begin
        all_param_syms = parameters(uni_uni)
        expected_fitted = Tuple(p for p in all_param_syms if p !== :E_total && p !== :Keq)

        @test EnzymeRates.fitted_params(uni_uni) == expected_fitted
        @test metabolites(uni_uni) == (:S, :P)
    end

    # ── FittingProblem construction ───────────────────────────────────────────
    @testset "Construction" begin
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs5)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        @test length(fp.log_abs_rates) == 5
        @test length(fp.group_point_indexes) == 1
        @test length(fp.group_point_indexes[1]) == 5
    end

    # ── Multi-group centering invariance ─────────────────────────────────────
    @testset "Multi-group centering invariance" begin
        # Two groups, each independently scaled
        data1 = _testhelper_make_synthetic_data(uni_uni, true_params, concs5;
            groups=["G1","G1","G1","G2","G2"])
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
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs5)
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
        # Integer Keq and scale_k_to_kcat convert to the Float64 fields.
        fp = FittingProblem(uni_uni, ok_data; Keq=2, scale_k_to_kcat=3)
        @test fp.Keq === 2.0 && fp.scale_k_to_kcat === 3.0
    end

    # ── Sign-mismatch penalty ─────────────────────────────────────────────────
    # Regression test for all-mismatch groups: when every prediction in a
    # group is a sign mismatch, centering must not zero every deviation.
    # The loss must be nonzero to distinguish a bad mechanism from a perfect one.
    # Three sub-cases exercise each kind of mismatch that sets buf[i] = 10.0:
    #   (i)  pred == 0.0            (S=0, P=0)
    #   (ii) pred < 0, Rate > 0     (S=0, P>0: only reverse term survives → pred always negative)
    #   (iii) pred > 0, Rate < 0    (S>0, P=0: only forward term survives → pred always positive)
    # In all cases every point in the group is a mismatch so the expected
    # loss is: (0 from centering + 100.0 × n_mismatch) / n_data = 100.0
    @testset "All-mismatch group not cancelled by centering" begin
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

    # ── Zero allocations ──────────────────────────────────────────────────────
    # Measured inside a function so the count excludes boxing the Float64 return
    # of a dynamically dispatched call (Julia < 1.12 `@allocated` counts it).
    loss_allocs(x, fp) = @allocated EnzymeRates.loss!(x, fp)

    @testset "Zero allocations" begin
        concs_list = [(S = Float64(i), P = 0.1) for i in 1:20]
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        x = randn(length(EnzymeRates.fitted_params(uni_uni)))
        loss_allocs(x, fp)  # warmup
        allocs = loss_allocs(x, fp)
        @test allocs == 0

        # Absolute mode (uncentered branch) is equally allocation-free.
        fp_abs = FittingProblem(uni_uni, data; Keq=Keq_val, scale_k_to_kcat=nothing)
        loss_allocs(x, fp_abs)  # warmup
        allocs_abs = loss_allocs(x, fp_abs)
        @test allocs_abs == 0
    end

    # loss! builds each data point's concentration NamedTuple from metabolites(m), so an
    # allosteric mechanism's metabolite names must be a compile-time constant as well.
    @testset "Zero allocations: allosteric mechanism" begin
        data = _testhelper_make_synthetic_data(allo, allo_params, allo_concs)
        fp = FittingProblem(allo, data; Keq=Keq_val)

        x = randn(length(allo_fitted))
        loss_allocs(x, fp)  # warmup
        allocs = loss_allocs(x, fp)
        @test allocs == 0
    end

    # ── Optimizer boundary ────────────────────────────────────────────────────
    # fit_rate_equation hands Optimization.jl one problem type whatever the mechanism and
    # its data, so solve and the solver compile once per optimizer, not once per
    # mechanism. The objective must still return loss! and allocate nothing, for a Vector
    # and for the column view of its population that CMA-ES passes.
    function _testhelper_optimization_problem(fp)
        n = length(EnzymeRates.fitted_params(fp.mechanism))
        EnzymeRates._optimization_problem(fp, zeros(n), fill(-15.0, n), fill(15.0, n))
    end
    objective_allocs(prob, x) = @allocated prob.f(x, prob.p)

    @testset "Optimizer boundary" begin
        fp_uni = FittingProblem(uni_uni,
            _testhelper_make_synthetic_data(uni_uni, true_params, concs5); Keq=Keq_val)
        fp_allo = FittingProblem(allo,
            _testhelper_make_synthetic_data(allo, allo_params, allo_concs); Keq=Keq_val)

        @testset "one problem type for every mechanism" begin
            @test typeof(_testhelper_optimization_problem(fp_uni)) ===
                  typeof(_testhelper_optimization_problem(fp_allo))
        end

        @testset "objective returns loss! and allocates nothing" begin
            for fp in (fp_uni, fp_allo)
                prob = _testhelper_optimization_problem(fp)
                population = randn(length(prob.u0), 3)
                x = population[:, 2]
                for x_in in (x, view(population, :, 2))
                    @test prob.f(x_in, prob.p) == EnzymeRates.loss!(x, fp)
                    objective_allocs(prob, x_in)  # warmup
                    @test objective_allocs(prob, x_in) == 0
                end
            end
        end
    end

    # ── Speed benchmark ───────────────────────────────────────────────────────
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

        bb_Keq = 1.5
        bb_params_syms = parameters(ordered_bi_bi)
        # Generate random params
        param_vals = ntuple(i -> 1.0 + 9.0 * rand(), length(bb_params_syms))
        bb_true_params = NamedTuple{bb_params_syms}(param_vals)
        bb_true_params = merge(bb_true_params, (Keq = bb_Keq, E_total = 1.0))

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
        fp = FittingProblem(ordered_bi_bi, data; Keq=bb_Keq)

        x = randn(length(EnzymeRates.fitted_params(ordered_bi_bi)))
        EnzymeRates.loss!(x, fp)  # warmup/compile

        # Minimum over several batches defeats the GC/scheduling inflation a
        # single mean suffers (matches _testhelper_test_rate_equation_performance); the
        # accumulator's finiteness check keeps the calls from being elided.
        best_us = Inf
        acc = 0.0
        for _ in 1:5
            acc = 0.0
            t = @elapsed for _ in 1:2000; acc += EnzymeRates.loss!(x, fp); end
            best_us = min(best_us, t / 2000 * 1e6)
        end
        isfinite(acc) || error("loss! produced a non-finite result")
        @test best_us < 70  # 500 datapoints; ~5 μs typical local, min-of-batches strips CI noise
    end

    # ── fit_rate_equation: scale_k_to_kcat anchors kcat; nothing leaves it raw ──
    # Anchoring rescales every fitted point to the target kcat, so the assertions hold at
    # any fit quality and a half-second budget is enough.
    @testset "scale_k_to_kcat normalization + retcode" begin
        using OptimizationBBO
        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
            (S = 0.5, P = 0.5), (S = 1.0, P = 0.5), (S = 2.0, P = 0.5),
        ]
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs_list)
        opt = BBO_adaptive_de_rand_1_bin_radiuslimited()

        # Default target 1.0 (best of 2 restarts) and a custom target 7.0: the returned
        # params are anchored so kcat ≈ target.
        for (fp_kwargs, n_restarts, target) in
                (((;), 2, 1.0), ((; scale_k_to_kcat=7.0), 1, 7.0))
            fp = FittingProblem(uni_uni, data; Keq=Keq_val, fp_kwargs...)
            res = _testhelper_fit_capturing_convergence() do
                fit_rate_equation(fp, opt; n_restarts, maxtime=0.5)
            end
            @test keys(res.params) == EnzymeRates.fitted_params(uni_uni)
            @test isfinite(res.loss)
            @test res.retcode isa Symbol
            full = merge(res.params, (Keq = Keq_val, E_total = 1.0))
            @test EnzymeRates._kcat_forward(uni_uni, full) ≈ target rtol=0.01
        end

        # scale_k_to_kcat=nothing: params returned verbatim (data fixes the scale).
        fpN = FittingProblem(uni_uni, data; Keq=Keq_val, scale_k_to_kcat=nothing)
        resN = _testhelper_fit_capturing_convergence() do
            fit_rate_equation(fpN, opt; n_restarts=1, maxtime=0.5)
        end
        @test keys(resN.params) == EnzymeRates.fitted_params(uni_uni)
        @test resN.retcode isa Symbol
    end

    # ── solver-option forwarding (named commons + solver_kwargs) ──
    @testset "solver kwarg forwarding" begin
        using OptimizationCMAEvolutionStrategy
        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
        ]
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs_list)
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

    # ── maxtime forwarded from fit_rate_equation to Optimization.solve (§6) ─
    @testset "maxtime forwarded to Optimization.solve" begin
        using Optimization
        using Optimization.SciMLBase: build_solution, ReturnCode, DefaultOptimizationCache

        # A stub optimizer that records the `maxtime` kwarg it receives, so the
        # forwarding path (fit_rate_equation -> the `common = (; maxtime,
        # maxiters)` merge -> Optimization.solve) can be asserted end-to-end
        # without depending on a real solver's behavior.
        mutable struct _testhelper_MaxtimeStubOpt
            maxtime_seen::Union{Nothing, Real}
        end
        _testhelper_MaxtimeStubOpt() = _testhelper_MaxtimeStubOpt(nothing)
        Optimization.allowsbounds(::_testhelper_MaxtimeStubOpt) = true
        function Optimization.SciMLBase.__solve(
                prob::Optimization.OptimizationProblem,
                opt::_testhelper_MaxtimeStubOpt; kwargs...)
            opt.maxtime_seen = kwargs[:maxtime]
            u = zeros(length(prob.u0))
            cache = DefaultOptimizationCache(prob.f, prob.p)
            build_solution(cache, opt, u, prob.f(u, prob.p); retcode = ReturnCode.Success)
        end

        concs_list = [(S = 1.0, P = 0.1), (S = 2.0, P = 0.1)]
        data = _testhelper_make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        stub = _testhelper_MaxtimeStubOpt()
        fit_rate_equation(fp, stub; n_restarts=1, maxtime=1.23)
        @test stub.maxtime_seen == 1.23
    end

    # ── Validation errors ──────────────────────────────────────────────────────
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

        # A non-finite, missing or non-numeric rate: the error names the Rate column and
        # the row, and shows the value as Julia prints it, so a String rate is quoted.
        for (bad, shown) in ((NaN, "NaN"), (Inf, "Inf"), (-Inf, "-Inf"),
                             (missing, "missing"), ("1.0", "\"1.0\""))
            data_bad = (group = ["G1", "G1"], Rate = [1.0, bad], S = [1.0, 2.0],
                        P = [0.1, 0.1])
            @test_throws(
                ErrorException("Rate at row 2 must be a finite number; got $shown"),
                FittingProblem(uni_uni, data_bad; Keq=1.0))
        end

        # A concentration that is not a finite number ≥ 0: the error names the column
        # and the row. Zero is a valid concentration.
        for (bad, shown) in ((NaN, "NaN"), (Inf, "Inf"), (-1.0, "-1.0"),
                             (missing, "missing"), ("1.0", "\"1.0\""))
            data_bad = (group = ["G1", "G1"], Rate = [1.0, 2.0], S = [1.0, bad],
                        P = [0.1, 0.1])
            @test_throws(
                ErrorException(
                    "Concentration S at row 2 must be a finite number ≥ 0; got $shown"),
                FittingProblem(uni_uni, data_bad; Keq=1.0))
        end
        data_zero_S = (group = ["G1", "G1"], Rate = [1.0, 2.0], S = [0.0, 1.0],
                       P = [0.1, 0.1])
        @test FittingProblem(uni_uni, data_zero_S; Keq=1.0) isa FittingProblem

        # A Keq or a scale_k_to_kcat that is not positive and finite
        data_ok = (group = ["G1"], Rate = [1.0], S = [1.0], P = [0.1])
        for (bad, shown) in ((0, "0"), (-1, "-1"), (Inf, "Inf"), (NaN, "NaN"))
            @test_throws(ErrorException("Keq must be positive and finite; got $shown"),
                FittingProblem(uni_uni, data_ok; Keq=bad))
            @test_throws(
                ErrorException("scale_k_to_kcat must be positive and finite (or " *
                               "nothing); got $shown"),
                FittingProblem(uni_uni, data_ok; Keq=1.0, scale_k_to_kcat=bad))
        end
    end

end
