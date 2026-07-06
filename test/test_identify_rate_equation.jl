# ABOUTME: Tests for identify_rate_equation pipeline.
# ABOUTME: Covers construction, helpers, and full pipeline
# ABOUTME: with mechanism recovery including allosteric path.

using DataFrames
using CSV
using Random
using Statistics
using OptimizationCMAEvolutionStrategy
using Optimization
using Optimization.SciMLBase: build_solution, ReturnCode, DefaultOptimizationCache

@testset "identify_rate_equation" begin

    # ── Shared test setup ────────────────────────────
    # Allosteric K-type uni-uni with regulator R:
    # S, P bind only R-state; R binds only T-state.
    # T-state cannot catalyze (K-type allosteric).
    test_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: R
    end

    # Build the constrained allosteric mechanism: K-type
    # allosteric with S, P only in R-state (`:OnlyA` group
    # tags) and R only in T-state (`:OnlyI` ligand tag).
    _base = first(EnzymeRates.init_mechanisms(test_rxn))
    _cat_allo_states = Symbol[]
    for g in EnzymeRates.kinetic_groups(_base)
        rep = EnzymeRates.rep_step(_base, g)
        met = EnzymeRates.bound_metabolite(rep)
        tag = (met isa EnzymeRates.Reactant) ? :OnlyA : :NonequalAI
        push!(_cat_allo_states, tag)
    end
    _site = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:R)], 1, [:OnlyI])
    _am = EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(_base), copy(EnzymeRates.steps(_base)),
        _cat_allo_states, 1, [_site])
    test_mechanism = EnzymeRates.AllostericEnzymeMechanism(_am)

    Keq_val = 2.0
    # 5 fitted params: K_A_P_E, K_A_S_E, k_A_ES_to_EP, K_I_Rreg, L
    true_params = (
        K_A_P_E = 1.0, K_A_S_E = 0.5, k_A_ES_to_EP = 5.0,
        K_I_Rreg = 2.0, L = 0.1,
        Keq = Keq_val, E_total = 1.0)

    function make_test_data(
        mechanism, params;
        n_per_group=10, n_groups=5
    )
        groups = String[]
        rates = Float64[]
        S_vals = Float64[]
        P_vals = Float64[]
        R_vals = Float64[]
        for g in 1:n_groups
            for _ in 1:n_per_group
                s = 0.1 + 9.9 * rand()
                p = 0.1 + 9.9 * rand()
                r = 0.1 + 9.9 * rand()
                concs = (S = s, P = p, R = r)
                v = rate_equation(
                    mechanism, concs, params)
                push!(groups, "G$g")
                push!(rates, v)
                push!(S_vals, s)
                push!(P_vals, p)
                push!(R_vals, r)
            end
        end
        return (group = groups, Rate = rates,
                S = S_vals, P = P_vals,
                R = R_vals)
    end

    Random.seed!(42)
    test_data = make_test_data(
        test_mechanism, true_params)

    cmaes_opt = CMAEvolutionStrategyOpt()

    # ── Construction validation ──────────────────────
    @testset "construction" begin
        prob = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val)
        @test prob.reaction === test_rxn
        @test prob.Keq == Keq_val
        @test length(
            unique(prob.data.group)) == 5

        # scale_k_to_kcat field: default 1.0, settable to nothing, validated.
        @test prob.scale_k_to_kcat == 1.0
        prob_abs = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, scale_k_to_kcat=nothing)
        @test prob_abs.scale_k_to_kcat === nothing
        @test_throws ErrorException IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, scale_k_to_kcat=0.0)

        # Missing metabolite column
        @test_throws(
            ErrorException,
            IdentifyRateEquationProblem(
                test_rxn,
                (group = ["G1"], Rate = [1.0]);
                Keq=1.0))

        # Missing regulator column
        @test_throws(
            ErrorException,
            IdentifyRateEquationProblem(
                test_rxn,
                (group = ["G1", "G2"],
                 Rate = [1.0, 2.0],
                 S = [1.0, 2.0],
                 P = [0.1, 0.1]);
                Keq=1.0))

        # Zero rate
        @test_throws(
            ErrorException,
            IdentifyRateEquationProblem(
                test_rxn,
                (group = ["G1", "G2"],
                 Rate = [0.0, 1.0],
                 S = [1.0, 1.0],
                 P = [0.1, 0.1],
                 R = [1.0, 1.0]);
                Keq=1.0))

        # Need >= 2 groups
        @test_throws(
            ErrorException,
            IdentifyRateEquationProblem(
                test_rxn,
                (group = ["G1", "G1"],
                 Rate = [1.0, 2.0],
                 S = [1.0, 2.0],
                 P = [0.1, 0.1],
                 R = [1.0, 1.0]);
                Keq=1.0))
    end

    # ── Helper unit tests (cheap, no fitting) ────────
    @testset "_rows_to_dataframe" begin
        rows = [(
            n_params = 3,
            loss = 0.5,
            mechanism_type = "test",
            rate_equation = "v = ...",
            retcode = "Success",
            error = missing,
            fitted_param_names = (:a, :b),
            fitted_param_values = (1.0, 2.0),
            eq_hash = "0123456789abcdef",
            fit_inherited = false,
        )]
        df = EnzymeRates._rows_to_dataframe(
            rows)
        @test nrow(df) == 1
        @test "a" in names(df)
        @test "b" in names(df)
        @test "eq_hash" in names(df)
        @test "retcode" in names(df)
        @test "error" in names(df)
        @test "fit_inherited" in names(df)
        @test df.fit_inherited == [false]
        @test !("fit_inherited_from_estimate" in names(df))

        # Empty rows
        df2 = EnzymeRates._rows_to_dataframe(
            NamedTuple[])
        @test nrow(df2) == 0
    end

    @testset "_rows_to_dataframe with failure row" begin
        rows = [
            (n_params = 3, loss = 0.5, mechanism_type = "M",
             rate_equation = "v = ...", retcode = "Success", error = missing,
             fitted_param_names = (:a,), fitted_param_values = (1.0,),
             eq_hash = "0123456789abcdef", fit_inherited = false),
            (n_params = missing, loss = missing, mechanism_type = "M",
             rate_equation = missing, retcode = missing,
             error = "StackOverflowError: ", fitted_param_names = (),
             fitted_param_values = (), eq_hash = missing,
             fit_inherited = missing),
        ]
        df = EnzymeRates._rows_to_dataframe(rows)
        @test nrow(df) == 2
        @test ismissing(df.loss[2])
        @test df.error[2] == "StackOverflowError: "
        @test ismissing(df.retcode[2])
        @test ismissing(df.eq_hash[2])
        @test df.fit_inherited[1] == false
        @test ismissing(df.fit_inherited[2])
        @test "a" in names(df)              # param column still built from row 1
        @test ismissing(df.a[2])            # failure row contributes no param value
    end

    @testset "failure row preserves round-trippable mechanism" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        m = first(EnzymeRates.init_mechanisms(rxn))
        f = EnzymeRates.FitFailure(m, "boom")
        row = EnzymeRates._failure_row(f)
        # Round-trippable parametric Sig, not the bare concrete type name.
        @test row.mechanism_type == string(typeof(EnzymeRates.compile_mechanism(m)))
        @test row.mechanism_type != "EnzymeRates.Mechanism"
        T = Core.eval(EnzymeRates, Meta.parse(row.mechanism_type))
        @test EnzymeRates.Mechanism(T()) == m
        @test row.error == "boom"
    end

    @testset "_rate_eq_dedup_key" begin
        base = "(; K_a, k_b) = params\n(; A) = concs\n" *
               "# Haldane constraints:\nk_r = (1/Keq)*K_a\nv = k_b*A/K_a"
        # differs only in a comment header + a substituted-into-v provenance line:
        a = "# Wegscheider constraints:\nK_x = K_a  (substituted into v)\n" * base
        b = "# Wegscheider constraints:\nK_y = K_a  (substituted into v)\n" * base
        @test EnzymeRates._rate_eq_dedup_key(a) ==
              EnzymeRates._rate_eq_dedup_key(b)
        # differs in a Haldane definition -> different key:
        c = replace(base, "k_r = (1/Keq)*K_a" => "k_r = (2/Keq)*K_a")
        @test EnzymeRates._rate_eq_dedup_key(base) !=
              EnzymeRates._rate_eq_dedup_key(c)
        # differs in the v= line -> different key:
        d = replace(base, "v = k_b*A/K_a" => "v = k_b*A/(K_a + A)")
        @test EnzymeRates._rate_eq_dedup_key(base) !=
              EnzymeRates._rate_eq_dedup_key(d)
    end

    # ── Run pipeline ONCE, test everything ───────────
    prob = IdentifyRateEquationProblem(
        test_rxn, test_data; Keq=Keq_val)
    save_dir = mktempdir()
    # Smoke-test settings: greedy beam (min_beam_width=1 +
    # tightest thresholds) so only the strictly-best mechanism
    # passes per level. Tests verify the pipeline runs and
    # produces correct shape — they don't require an exhaustive
    # search. Light n_restarts/maxtime keep each fit under ~1s.
    results = identify_rate_equation(prob;
        min_beam_width=1,
        loss_rel_threshold=1.0,
        loss_abs_threshold=0.0,
        max_param_count=8,
        n_cv_candidates=1,
        save_dir=save_dir,
        optimizer=cmaes_opt,
        n_restarts=1, maxtime=1.0)

    @testset "mechanism recovery" begin
        # The best mechanism should fit the noiseless data with near-zero loss.
        # (Whether the search recovers the *most parsimonious* mechanism is a
        # search-quality property that needs heavy, seeded fits to test
        # reliably — too slow and too stochastic for CI, so it is not asserted
        # here.)
        fp_best = FittingProblem(
            results.best, test_data;
            Keq=Keq_val)
        fit_best = fit_rate_equation(
            fp_best, cmaes_opt;
            n_restarts=3, maxtime=10.0)
        @test fit_best.loss < 0.01
    end

    @testset "results structure" begin
        @test results isa
            IdentifyRateEquationResults
        @test results.best isa
            EnzymeRates.AbstractEnzymeMechanism
        @test nrow(results.cv_results) > 0
        @test "cv_score" in names(
            results.cv_results)
        @test all(
            isfinite,
            results.cv_results.cv_score)
        @test "n_params" in names(
            results.cv_results)
        @test "loss" in names(
            results.cv_results)
        @test "eq_hash" in names(
            results.cv_results)
        # LOOCV candidate dedup invariant: within each n_params
        # bucket, each eq_hash should appear at most once (the
        # `seen_hashes in continue` filter in
        # `_cv_model_selection`). This catches a regression where
        # duplicates would enter LOOCV and waste compute / bias
        # the per-bucket "best".
        for gdf in groupby(
                results.cv_results, :n_params)
            @test allunique(gdf.eq_hash)
        end

        # Diagnostic columns from paired 1-SE + permutation rule.
        @test "mean_log_loss_diff" in
            names(results.cv_results)
        @test "se_paired" in names(results.cv_results)
        @test "permutation_p" in names(results.cv_results)

        # n_min bucket = the bucket with lowest cv_score after rep
        # selection (which equals lowest mean log-fold-loss). Its
        # rows have all three diagnostics = 0.0.
        n_min_val = results.cv_results.n_params[
            argmin(results.cv_results.cv_score)]
        n_min_rows = filter(row -> row.n_params == n_min_val,
                            results.cv_results)
        @test all(==(0.0),
                  n_min_rows.mean_log_loss_diff)
        @test all(==(0.0), n_min_rows.se_paired)
        @test all(==(0.0), n_min_rows.permutation_p)

        # Per-fold columns named by held-out group label.
        groups = unique(prob.data.group)
        for g in groups
            @test Symbol("cv_fold_$g") in
                propertynames(results.cv_results)
        end

        # CSV-roundtrip: cv_results must be CSV-serializable. Verify by
        # writing and re-reading; check column-name preservation and row count.
        buf = IOBuffer()
        CSV.write(buf, results.cv_results)
        seekstart(buf)
        roundtrip = CSV.read(buf, DataFrame)
        for col in (:n_params, :loss, :cv_score,
                    :mean_log_loss_diff, :se_paired,
                    :permutation_p)
            @test col in propertynames(roundtrip)
        end
        @test nrow(roundtrip) == nrow(results.cv_results)
    end

    @testset "best mechanism computes rates" begin
        req = rate_equation_string(results.best)
        @test length(req) > 0
    end

    @testset "CSV output (new schema)" begin
        files = sort(filter(f -> endswith(f, ".csv"), readdir(save_dir)))
        @test "initial_mechanisms.csv" in files
        @test isfile(joinpath(save_dir, "progress.log"))
        @test filesize(joinpath(save_dir, "progress.log")) > 0
        log_text = read(joinpath(save_dir, "progress.log"), String)
        @test occursin("new fits", log_text)
        @test occursin("skipped (>", log_text)
        @test occursin("best loss by n_params:", log_text)
        @test !any(startswith(f, "params_estimate_") for f in files)
        iters = filter(f -> startswith(f, "equation_search_iteration_"), files)
        @test !isempty(iters)
        nums = sort(parse.(Int, replace.(iters,
            "equation_search_iteration_" => "", ".csv" => "")))
        @test nums == collect(1:length(nums))      # sequential, no gaps
        init_df = CSV.read(joinpath(save_dir, "initial_mechanisms.csv"), DataFrame)
        @test nrow(init_df) == length(unique!(
            collect(EnzymeRates.init_mechanisms(prob.reaction))))
        @test "eq_hash" in names(init_df)
        @test !("fit_inherited_from_estimate" in names(init_df))
        for f in files
            df_file = CSV.read(joinpath(save_dir, f), DataFrame)
            @test "eq_hash" in names(df_file)
            @test all(length.(string.(skipmissing(df_file.eq_hash))) .== 16)
            @test all(<=(8), skipmissing(df_file.n_params))     # max_param_count=8
        end
    end

    @testset "loocv_results.csv and best_equation.csv saved" begin
        files = filter(f -> endswith(f, ".csv"), readdir(save_dir))
        @test "loocv_results.csv" in files
        @test "best_equation.csv" in files

        cvf = CSV.read(joinpath(save_dir, "loocv_results.csv"), DataFrame)
        @test nrow(cvf) == nrow(results.cv_results)
        @test "cv_score" in names(cvf)

        bestf = CSV.read(joinpath(save_dir, "best_equation.csv"), DataFrame)
        @test nrow(bestf) == 1
        best_hash = string(
            EnzymeRates._rate_eq_dedup_key(rate_equation_string(results.best)),
            base=16, pad=16)
        @test string(bestf.eq_hash[1]) == best_hash
    end

    @testset "save_dir non-empty check" begin
        # Should error before any fitting starts (the save_dir
        # validation runs up-front), so the heavy settings would
        # never matter — but use the same lean settings as
        # above for consistency.
        @test_throws(
            ErrorException,
            identify_rate_equation(prob;
                min_beam_width=1,
                loss_rel_threshold=1.0,
                loss_abs_threshold=0.0,
                max_param_count=8,
                n_cv_candidates=1,
                save_dir=save_dir,
                optimizer=cmaes_opt,
                n_restarts=1, maxtime=1.0))
    end

    @testset "_cv_fold_loss over all groups: per-fold scores, floored at eps" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
        end
        m = EnzymeRates.EnzymeMechanism(
            first(EnzymeRates.init_mechanisms(rxn)))

        # 3 groups × 2 rows each so per-fold fits aren't degenerate
        data = DataFrame(
            S    = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            P    = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            Rate = [0.5, 0.8, 1.0, 1.1, 1.2, 1.3],
            group = [1, 1, 2, 2, 3, 3],
        )
        prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)

        groups = unique(prob.data.group)
        scores = [EnzymeRates._cv_fold_loss(m, prob, g;
            optimizer=CMAEvolutionStrategyOpt(),
            n_restarts=2, maxtime=2.0, maxiters=500) for g in groups]

        @test scores isa Vector{Float64}
        # Require the success path: fitting MUST converge on
        # this trivial uni-uni fixture in 2s × 2 restarts. A
        # length-0 (full failure) result would let the per-fold
        # eps-floor + isfinite assertions below pass vacuously.
        @test length(scores) == 3
        @test all(s -> s >= eps(Float64), scores)
        @test all(isfinite, scores)
    end

    @testset "_cv_fold_loss is loud on fit failure" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
        end
        m = EnzymeRates.EnzymeMechanism(
            first(EnzymeRates.init_mechanisms(rxn)))
        data = DataFrame(
            S    = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            P    = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            Rate = [0.5, 0.8, 1.0, 1.1, 1.2, 1.3],
            group = [1, 1, 2, 2, 3, 3],
        )
        prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
        # An unrecognized kwarg (`beam_fraction`) makes the fold's
        # `fit_rate_equation` call throw; `_cv_fold_loss` must propagate that
        # error, not swallow it (a corrupted CV must abort model selection).
        @test_throws Exception EnzymeRates._cv_fold_loss(
            m, prob, first(unique(prob.data.group));
            optimizer=CMAEvolutionStrategyOpt(),
            n_restarts=1, maxtime=1.0, beam_fraction=0.5)
    end

    @testset "_cv_fold_loss: one fold, floored + finite" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
        end
        m = EnzymeRates.EnzymeMechanism(
            first(EnzymeRates.init_mechanisms(rxn)))
        data = DataFrame(
            S = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            P = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            Rate = [0.5, 0.8, 1.0, 1.1, 1.2, 1.3],
            group = [1, 1, 2, 2, 3, 3])
        prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
        kw = (; optimizer=CMAEvolutionStrategyOpt(),
                n_restarts=2, maxtime=2.0, maxiters=500)

        one = EnzymeRates._cv_fold_loss(m, prob, 2; kw...)
        @test one isa Float64
        @test one >= eps(Float64) && isfinite(one)
    end

end

@testset "csv writers" begin
    rows = [(
        n_params = 5, loss = 1.0, mechanism_type = "M",
        rate_equation = "v = 1", retcode = "Success", error = missing,
        fitted_param_names = (:K_a,),
        fitted_param_values = (2.0,), eq_hash = "abc",
        fit_inherited = false,
    )]
    mktempdir() do tmp
        EnzymeRates._save_initial_csv(tmp, rows)
        @test isfile(joinpath(tmp, "initial_mechanisms.csv"))
        EnzymeRates._save_iteration_csv(tmp, rows, 3)
        @test isfile(joinpath(tmp, "equation_search_iteration_3.csv"))
        df = CSV.read(joinpath(tmp, "equation_search_iteration_3.csv"), DataFrame)
        @test df.n_params == [5]
        @test "eq_hash" in names(df)
        # dir-creation branch: save_dir does not exist yet
        subdir = joinpath(tmp, "made")
        EnzymeRates._save_initial_csv(subdir, rows)
        @test isfile(joinpath(subdir, "initial_mechanisms.csv"))
    end
end

@testset "beam selection: loss thresholds + min_beam_width floor" begin
    losses = [1.0, 1.5, 2.5, 5.0, 10.0]
    sel = EnzymeRates._select_beam(
        losses;
        loss_rel_threshold=2.0,
        loss_abs_threshold=0.0,
        min_beam_width=1)
    @test sort(sel) == [1, 2]

    sel = EnzymeRates._select_beam(
        losses;
        loss_rel_threshold=2.0,
        loss_abs_threshold=0.0,
        min_beam_width=4)
    @test sort(sel) == [1, 2, 3, 4]

    losses_small = [1e-6, 0.005, 0.05]
    sel = EnzymeRates._select_beam(
        losses_small;
        loss_rel_threshold=2.0,
        loss_abs_threshold=0.01,
        min_beam_width=1)
    @test sort(sel) == [1, 2]

    sel = EnzymeRates._select_beam(
        [Inf, Inf, Inf];
        loss_rel_threshold=2.0,
        loss_abs_threshold=0.01,
        min_beam_width=5)
    @test isempty(sel)

    sel = EnzymeRates._select_beam(
        [1.0, NaN, 2.0];
        loss_rel_threshold=2.5,
        loss_abs_threshold=0.0,
        min_beam_width=1)
    @test sort(sel) == [1, 3]

    # `_select_beam` returns indices in INPUT order, not loss
    # order. Verify with a deliberately-shuffled input.
    sel = EnzymeRates._select_beam(
        [5.0, 1.0, 10.0, 2.0];
        loss_rel_threshold=2.5,
        loss_abs_threshold=0.0,
        min_beam_width=1)
    @test sel == [2, 4]   # input-order, not [2, 4] sorted by loss
end

@testset "all base fits fail: failure CSV written, then raises" begin
    # A `solver_kwargs` option no optimizer recognizes is forwarded
    # verbatim to `Optimization.solve` and makes every fit throw (a bogus
    # name keeps this independent of whether any real option is supported).
    # Per-mechanism fit failures are isolated in `_process_batch`, so the
    # base tier is then empty and the pipeline raises. The contract under
    # test is that the all-base-fail path persists the failure rows to
    # `initial_mechanisms.csv` before raising (for cluster debugging).
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (group = ["G1", "G1", "G2", "G2"],
            Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0],
            P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    tmp = mktempdir()
    @test_throws ErrorException identify_rate_equation(
        prob; solver_kwargs=(; not_a_real_solver_option=1),
        optimizer=CMAEvolutionStrategyOpt(),
        n_restarts=1, maxtime=1.0, save_dir=tmp)
    # Failure rows were written before the re-raise: a CSV exists whose rows
    # are all failures (non-missing `error`, missing `eq_hash`).
    @test isfile(joinpath(tmp, "initial_mechanisms.csv"))
    fail_df = CSV.read(joinpath(tmp, "initial_mechanisms.csv"), DataFrame)
    @test nrow(fail_df) >= 1
    @test all(.!ismissing.(fail_df.error))
    @test all(ismissing.(fail_df.eq_hash))
end

@testset "_onesided_permutation_p" begin
    # All-zero diffs: every sign flip yields perm_mean = observed = 0,
    # so count_ge = 2^n. p = 1.0.
    @test EnzymeRates._onesided_permutation_p(
        [0.0, 0.0, 0.0]) == 1.0

    # All-positive equal diffs: only the identity permutation matches
    # observed; every flipped variant gives a smaller mean.
    # p = 1/2^4 = 0.0625.
    @test EnzymeRates._onesided_permutation_p(
        [1.0, 1.0, 1.0, 1.0]) ≈ 1/16

    # All-negative diffs: observed = -1, all flips ≥ -1 → count_ge = 2^n.
    # p = 1.0.
    @test EnzymeRates._onesided_permutation_p(
        [-1.0, -1.0, -1.0]) == 1.0

    # Mixed-sign 8-fold fixture: 256 exact perms, p strictly in (0, 1).
    diffs = [0.10, -0.05, 0.08, -0.02,
             0.06, -0.04, 0.03, -0.01]
    p_exact = EnzymeRates._onesided_permutation_p(diffs)
    @test 0 < p_exact < 1

    # Force Monte Carlo path (exact_threshold=0) on the same diffs;
    # results must agree within sampling SE. With 10^6 samples and
    # p ≈ 0.5, SE on count_ge/N is √(0.25/10^6) ≈ 5e-4.
    p_mc = EnzymeRates._onesided_permutation_p(
        diffs;
        exact_threshold = 0,
        mc_samples = 10^6,
        rng = MersenneTwister(42),
    )
    @test abs(p_exact - p_mc) < 0.01

    # Determinism: a seeded RNG must produce bit-identical output
    # across runs. This proves the `rng` kwarg threads through both
    # the exact (no-op) and MC paths.
    p1 = EnzymeRates._onesided_permutation_p(
        diffs; exact_threshold = 0, mc_samples = 10^4,
        rng = MersenneTwister(7))
    p2 = EnzymeRates._onesided_permutation_p(
        diffs; exact_threshold = 0, mc_samples = 10^4,
        rng = MersenneTwister(7))
    @test p1 == p2

    # Regression: 20-fold all-positive equal-spaced diffs. Only the
    # identity permutation reproduces observed; all sign-flipped
    # variants give smaller s. Correct p = 1/2^20. The bug was that
    # `mean(diffs)` (pairwise sum) produced an observed value 1 ULP
    # larger than the loop's sequential sum, dropping the identity
    # from the count and returning 0.0 instead.
    long_diffs = collect(0.1:0.1:2.0)   # n=20, all positive
    @test EnzymeRates._onesided_permutation_p(long_diffs) ==
          1 / 2^20

    # Also verify n=16: the smallest n that triggers Julia's pairwise
    # path. With these inputs the pairwise/sequential sums happen to
    # agree, so this case passes pre-fix too — included to lock in the
    # boundary.
    @test EnzymeRates._onesided_permutation_p(
        collect(0.1:0.1:1.6)) == 1 / 2^16
end

@testset "_select_best_n_params: paired SE math" begin
    # Simple two-bucket case: n=7 best (lowest cv_score). n=5 paired
    # diffs = [+0.5, +0.5, +0.5, +0.5, +0.5, +0.5] (uniform offset).
    # mean_diff = 0.5, std_diff = 0, se_paired = 0. mean_diff > 0
    # → 1-SE rejects. Permutation: all-positive diffs → only the
    # identity perm reproduces observed; p = 1/2^6 = 0.015625 < 0.16
    # → perm rejects. Both fail → best_n = n_min = 7.
    cv_df = DataFrame(
        n_params       = [5, 7],
        cv_score       = [0.6, 0.1],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.6, 0.7, 0.5, 0.6, 0.5, 0.7]),
            exp.([0.1, 0.2, 0.0, 0.1, 0.0, 0.2]),
        ],
    )
    res = EnzymeRates._select_best_n_params(cv_df)
    @test res.n_min == 7
    @test res.best_n == 7
    # n_min self-comparison: hardcoded literal zeros, so === holds
    # regardless of FP noise in the input fold scores.
    @test res.diagnostics[7] ===
          (mean_log_loss_diff=0.0, se_paired=0.0,
           permutation_p=0.0)
    d5 = res.diagnostics[5]
    @test d5.mean_log_loss_diff ≈ 0.5
    # FP roundoff through exp/log makes std(diffs) ≈ 1e-17, not 0.0
    # exactly. `≈ 0.0` at default tolerance fails — use atol.
    @test isapprox(d5.se_paired, 0.0; atol = 1e-10)

    # Mixed-sign small diffs: n=7 best, n=5 paired diffs
    # = [0.0, 0.0, 0.03, -0.01]. mean = 0.005, std ≈ 0.01732,
    # se_paired = 0.01732/sqrt(4) = 0.00866. 0.005 ≤ 0.00866 → 1-SE
    # passes. Mixed-sign → permutation_p > 0.16 in 16 flips. Accept.
    cv_df2 = DataFrame(
        n_params       = [5, 7],
        cv_score       = [0.115, 0.110],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.10, 0.12, 0.13, 0.11]),
            exp.([0.10, 0.12, 0.10, 0.12]),
        ],
    )
    res2 = EnzymeRates._select_best_n_params(cv_df2)
    @test res2.n_min == 7
    @test res2.best_n == 5
    @test res2.diagnostics[5].mean_log_loss_diff ≈ 0.005
    @test res2.diagnostics[5].se_paired ≈
          std([0.0, 0.0, 0.03, -0.01]) / sqrt(4)

    # Multi-row bucket: rep is the row with lowest cv_score.
    # Bucket-3: row-A has fold scores giving mean_log = 0.135;
    # row-B has fold scores giving mean_log = 0.115 → row-B is rep.
    cv_df3 = DataFrame(
        n_params       = [3, 3, 7],
        cv_score       = [0.135, 0.115, 0.110],
        loss           = [0.0, 0.0, 0.0],
        cv_fold_scores = [
            exp.([0.13, 0.135, 0.14, 0.135]),
            exp.([0.11, 0.12, 0.115, 0.115]),
            exp.([0.09, 0.13, 0.10, 0.12]),
        ],
    )
    res3 = EnzymeRates._select_best_n_params(cv_df3)
    @test res3.diagnostics[3].mean_log_loss_diff ≈
          mean(log.(exp.([0.11, 0.12, 0.115, 0.115])) .-
               log.(exp.([0.09, 0.13, 0.10, 0.12])))

    # Single-fold case: n_folds_min = 1 → return n_min, no comparisons.
    cv_df4 = DataFrame(
        n_params       = [3, 5],
        cv_score       = [0.5, 0.1],
        loss           = [0.0, 0.0],
        cv_fold_scores = [exp.([0.5]), exp.([0.1])],
    )
    @test EnzymeRates._select_best_n_params(cv_df4).best_n == 5

    # Single-bucket cv_df: only one n_params value → return it as both
    # n_min and best_n. Diagnostics has exactly one entry, all zeros.
    cv_df_single = DataFrame(
        n_params       = [4],
        cv_score       = [0.15],
        loss           = [0.0],
        cv_fold_scores = [exp.([0.1, 0.2, 0.15, 0.12])],
    )
    res_single = EnzymeRates._select_best_n_params(cv_df_single)
    @test res_single.n_min == 4
    @test res_single.best_n == 4
    @test length(res_single.diagnostics) == 1
    @test res_single.diagnostics[4].mean_log_loss_diff == 0.0

    # Tie in mean log-fold-loss across buckets: n_min resolves to
    # smallest n_params (parsimony tiebreak). Both buckets have
    # identical fold scores → identical log-means → tie.
    cv_df_tie = DataFrame(
        n_params       = [3, 5, 7],
        cv_score       = [0.115, 0.115, 0.115],
        loss           = [0.0, 0.0, 0.0],
        cv_fold_scores = [
            exp.([0.10, 0.12, 0.11, 0.13]),
            exp.([0.10, 0.12, 0.11, 0.13]),
            exp.([0.10, 0.12, 0.11, 0.13]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df_tie).n_min == 3

    # Larger-than-best bucket: n_min is the middle bucket. Both the
    # smaller (n=3) and larger (n=7) buckets get diagnostics, but only
    # the smaller can be selected (loop is over `smaller_ns`).
    # Construct: n=5 has lowest log-mean. n=3 has uniform +0.5 offset
    # (FF — both gates fail). n=7 has uniform +0.3 offset (would also
    # fail every gate, but the loop never visits it).
    cv_df_three = DataFrame(
        n_params       = [3, 5, 7],
        cv_score       = [0.6, 0.1, 0.4],
        loss           = [0.0, 0.0, 0.0],
        cv_fold_scores = [
            exp.([0.6, 0.6, 0.6, 0.6, 0.6, 0.6]),
            exp.([0.1, 0.1, 0.1, 0.1, 0.1, 0.1]),
            exp.([0.4, 0.4, 0.4, 0.4, 0.4, 0.4]),
        ],
    )
    res_three = EnzymeRates._select_best_n_params(cv_df_three)
    @test res_three.n_min == 5
    @test res_three.best_n == 5
    # Larger bucket has populated diagnostics with positive mean_diff.
    @test haskey(res_three.diagnostics, 7)
    @test res_three.diagnostics[7].mean_log_loss_diff ≈ 0.3
    # Smaller bucket also has populated diagnostics (rejected by gates).
    @test haskey(res_three.diagnostics, 3)
    @test res_three.diagnostics[3].mean_log_loss_diff ≈ 0.5
end

@testset "_select_best_n_params: AND-combiner truth table" begin
    # Cell 1: pass-pass — mixed-sign small diffs, simpler accepted.
    # n=7 has strictly lower log-mean (0.11333) than n=3 (0.115) so n=7
    # is n_min; diffs ≈ [0, 0, 0, 0, +0.04, -0.03] → mean ≈ 0.00167,
    # se_paired ≈ 0.00910 (1-SE pass at default). perm_p = 0.5 > 0.16
    # (perm pass at default) → both gates pass → simpler bucket selected.
    cv_df_pp = DataFrame(
        n_params       = [3, 7],
        cv_score       = [0.115, 0.110],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.10, 0.13, 0.11, 0.12, 0.13, 0.10]),
            exp.([0.10, 0.13, 0.11, 0.12, 0.09, 0.13]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df_pp).best_n == 3

    # Cell 4: fail-fail — uniform large positive diffs, simpler rejected.
    cv_df_ff = DataFrame(
        n_params       = [3, 7],
        cv_score       = [1.5, 0.1],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([1.5, 1.5, 1.5, 1.5, 1.5, 1.5]),
            exp.([0.1, 0.1, 0.1, 0.1, 0.1, 0.1]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df_ff).best_n == 7

    # Cell 2: 1-SE pass + perm fail — crank perm_p_threshold to 0.99
    # so any non-degenerate p fails the perm gate while 1-SE still
    # passes on the mean=0 fixture.
    @test EnzymeRates._select_best_n_params(
        cv_df_pp; perm_p_threshold = 0.99).best_n == 7

    # Cell 3: 1-SE fail + perm pass — needs strictly positive mean
    # with mixed signs and enough variance that perm p stays above
    # threshold. Hand-computed fixture:
    #   n=3 log-folds = [0.105, 0.098, 0.103, 0.099]
    #   n=7 log-folds = [0.100, 0.100, 0.100, 0.100]
    #   diffs = [0.005, -0.002, 0.003, -0.001]
    #   mean = 0.00125; std ≈ 0.003304; se_paired ≈ 0.001652.
    # 1-SE default: 0.00125 ≤ 1.0*0.001652 ✓ pass.
    # Force fail: se_threshold=0.5 → require 0.00125 ≤ 0.000826 ✗.
    # Permutation (16 exact perms): 5 perms have perm_mean ≥ 0.00125
    # → p = 5/16 = 0.3125 > 0.16 ✓ pass.
    cv_df_marginal = DataFrame(
        n_params       = [3, 7],
        cv_score       = [0.10125, 0.100],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.105, 0.098, 0.103, 0.099]),
            exp.([0.100, 0.100, 0.100, 0.100]),
        ],
    )
    @test EnzymeRates._select_best_n_params(
        cv_df_marginal).best_n == 3   # both pass default
    @test EnzymeRates._select_best_n_params(
        cv_df_marginal; se_threshold = 0.5).best_n == 7
end

@testset "_select_best_n_params: edge cases" begin
    # Length mismatch between buckets → error. Fold-count must be
    # uniform across buckets because pairs (same held-out group →
    # same fold index) are required for the paired diffs.
    cv_df_mismatch = DataFrame(
        n_params       = [3, 7],
        cv_score       = [0.116, 0.113],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.105, 0.108, 0.115, 0.135]),                 # 4
            exp.([0.10, 0.105, 0.11, 0.115, 0.12, 0.125]),      # 6
        ],
    )
    @test_throws ErrorException EnzymeRates._select_best_n_params(
        cv_df_mismatch)
end

@testset "cv_results: exotic group labels survive CSV roundtrip" begin
    # The column-flattening step in _cv_model_selection does:
    #   col = Symbol("cv_fold_$g")
    #   cv_df[!, col] = [v[i] for v in cv_df.cv_fold_scores]
    # Exotic group labels (containing =, ,, spaces) must produce
    # valid Symbol column names that survive CSV.write/CSV.read.
    df = DataFrame(n_params = [3, 5])
    exotic_groups = ["a=b", "c,d", "x y"]
    fold_scores = [[0.1, 0.2, 0.3], [0.05, 0.1, 0.15]]
    for (i, g) in enumerate(exotic_groups)
        col = Symbol("cv_fold_$g")
        df[!, col] = [v[i] for v in fold_scores]
    end
    @test Symbol("cv_fold_a=b") in propertynames(df)
    @test Symbol("cv_fold_c,d") in propertynames(df)
    @test Symbol("cv_fold_x y") in propertynames(df)

    buf = IOBuffer()
    CSV.write(buf, df)
    seekstart(buf)
    roundtrip = CSV.read(buf, DataFrame)
    @test "cv_fold_a=b" in names(roundtrip)
    @test "cv_fold_c,d" in names(roundtrip)
    @test "cv_fold_x y" in names(roundtrip)
end

@testset "_default_save_dir" begin
    mktempdir() do tmp
        cd(tmp) do
            d1 = EnzymeRates._default_save_dir()
            @test occursin(r"^\d{4}_\d{2}_\d{2}_results$", d1)
            mkpath(d1)
            d2 = EnzymeRates._default_save_dir()
            @test d2 == d1 * "_2"
            mkpath(d2)
            @test EnzymeRates._default_save_dir() == d1 * "_3"
        end
    end
end

@testset "rate-eq dedup-key partition stability" begin
    # Representative reactions exercising the dedup key's edge cases:
    # - uni_uni: trivial structural equivalence
    # - bi_bi:   substituted-into-v ties across multiple kinetic groups
    # ter-ter intentionally omitted — `rate_equation_string` derivation is
    # extremely slow for mechanisms with >~30 enzyme forms (CLAUDE.md
    # "Known Issues"), and the dedup key renders that string per
    # candidate. The bi-bi enumeration already covers every structural
    # symmetry the dedup key collapses.
    test_reactions = [
        ("uni_uni", @enzyme_reaction(begin
            substrates: S[C]
            products:   P[C]
        end)),
        ("bi_bi", @enzyme_reaction(begin
            substrates: A[C], B[N]
            products:   P[C], Q[N]
        end)),
    ]

    # Expected partition sizes per reaction = the number of DISTINCT rate
    # equations the init-level enumeration produces. The 55 bi_bi init
    # mechanisms are all structurally distinct AND each yields a distinct
    # `rate_equation_string`, so the comment-stripped string key produces
    # exactly 55 classes (zero over- and zero under-collapse): clean
    # topologies have distinct enzyme-form sets, hence distinct rate
    # equations.
    # If these counts change in a future commit, the dedup key's
    # equivalence classes (or the enumeration) have shifted — investigate.
    expected_n_classes = Dict(
        "uni_uni" => 1,
        "bi_bi"   => 55,
    )

    for (label, reaction) in test_reactions
        # init_mechanisms only — skip expand_mechanisms. The init level
        # already produces multiple structurally-equivalent variants
        # (mirror-step orderings, kinetic-group renumberings) that
        # exercise the dedup key's collapse rules. expand_mechanisms
        # adds variants at higher param counts whose dedup-key
        # behavior is the same modulo size, at exponential compile cost.
        all_mechs = EnzymeRates.init_mechanisms(reaction)

        new_buckets = Dict{UInt64, Vector{Int}}()
        for (i, m) in enumerate(all_mechs)
            em = EnzymeRates.compile_mechanism(m)
            h = EnzymeRates._rate_eq_dedup_key(rate_equation_string(em))
            push!(get!(new_buckets, h, Int[]), i)
            # Determinism: same input, same key across invocations.
            @test EnzymeRates._rate_eq_dedup_key(rate_equation_string(em)) === h
        end

        @test length(new_buckets) == expected_n_classes[label]
    end
end

@testset "_select_beam best_override" begin
    losses = [1.0, 1.5, 3.0]
    kw = (loss_rel_threshold=1.2, loss_abs_threshold=0.0, min_beam_width=1)
    # without override: best = min = 1.0, cutoff = 1.2 -> only index 1
    @test EnzymeRates._select_beam(losses; kw...) == [1]
    # override best = 2.0 -> cutoff 2.4 -> indices 1 and 2
    @test EnzymeRates._select_beam(losses; kw..., best_override=2.0) == [1, 2]
    # min_beam_width still honored
    @test EnzymeRates._select_beam(losses;
        loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        min_beam_width=2, best_override=0.0) == [1, 2]
end

@testset "_select_beam parsimony_cutoff" begin
    # Floor guarantee: a parsimony_cutoff below every loss admits nothing via
    # the loss filter, yet min_beam_width still keeps the top-k by loss.
    losses = [1.0, 1.5, 2.5, 5.0, 10.0]
    @test EnzymeRates._select_beam(losses;
        loss_rel_threshold=2.0, loss_abs_threshold=0.0,
        min_beam_width=2, parsimony_cutoff=0.5) == [1, 2]

    # Tightening: a parsimony_cutoff stricter than the rel/abs cutoff lowers the
    # combined cutoff to 2.0, so indices 1 and 2 (losses 1.0, 1.5) pass and
    # index 3 (2.5) is dropped. Without it, rel=10 would admit all four.
    losses = [1.0, 1.5, 2.5, 5.0]
    @test EnzymeRates._select_beam(losses;
        loss_rel_threshold=10.0, loss_abs_threshold=0.0,
        min_beam_width=1, parsimony_cutoff=2.0) == [1, 2]

    # No-op: parsimony_cutoff=nothing reproduces the parsimony-free selection.
    kw = (loss_rel_threshold=2.0, loss_abs_threshold=0.0, min_beam_width=1)
    @test EnzymeRates._select_beam(losses; kw..., parsimony_cutoff=nothing) ==
          EnzymeRates._select_beam(losses; kw...)

    # Interaction: min() picks the smaller cutoff. With best_override=2.0 the
    # rel cutoff is 2.4 (admits 1,2); a tighter parsimony_cutoff=1.0 overrides
    # it down to just the single best.
    losses = [1.0, 1.5, 3.0]
    ov = (loss_rel_threshold=1.2, loss_abs_threshold=0.0,
          min_beam_width=1, best_override=2.0)
    @test EnzymeRates._select_beam(losses; ov...) == [1, 2]
    @test EnzymeRates._select_beam(losses; ov..., parsimony_cutoff=1.0) == [1]
end

@testset "_select_count! cumulative per-count floor" begin
    expanded = Dict{Int,Int}()
    # Sweep 1 at count 5: rel cutoff admits only the best (loss 1.0); the
    # floor budget (3) tops it up to the top 3 by loss. expanded[5] -> 3.
    sel1 = EnzymeRates._select_count!(expanded, 5, [1.0, 2.0, 3.0, 4.0, 5.0];
        loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        min_beam_width=3, best_override=1.0)
    @test sort(sel1) == [1, 2, 3]
    @test expanded[5] == 3

    # Sweep 2 at count 5: budget spent (3 of 3). New mechanisms all above the
    # cutoff -> the floor admits NONE (unlike the old per-sweep floor, which
    # would grant a fresh 3). expanded[5] stays 3.
    sel2 = EnzymeRates._select_count!(expanded, 5, [10.0, 11.0, 12.0];
        loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        min_beam_width=3, best_override=1.0)
    @test isempty(sel2)
    @test expanded[5] == 3

    # A cutoff-passer is still admitted after the floor is spent.
    sel3 = EnzymeRates._select_count!(expanded, 5, [1.0, 20.0];
        loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        min_beam_width=3, best_override=1.0)
    @test sel3 == [1]
    @test expanded[5] == 4
end

@testset "§1 _parsimony_cutoff = threshold * min over all counts < c" begin
    f = EnzymeRates._parsimony_cutoff
    @test f(Dict(5=>0.02), 5, 1.01) === nothing            # no count < c
    @test f(Dict(5=>0.02,6=>0.05,7=>0.03), 8, 1.01) ≈ 1.01*0.02   # min over <c, not c-1
    @test f(Dict(5=>0.02), 7, 1.01) ≈ 1.01*0.02            # count gap: c-1=6 absent
    @test f(Dict(5=>0.01,6=>0.04), 7, 1.01) ≈ 1.01*0.01    # non-monotone → true min
end

@testset "_progress" begin
    mktempdir() do tmp
        # show_progress=true: writes to progress.log AND to stdout.
        out_file = joinpath(tmp, "stdout.txt")
        open(out_file, "w") do io
            redirect_stdout(io) do
                EnzymeRates._progress(tmp, true, "stage one")
            end
        end
        @test occursin("stage one", read(out_file, String))
        @test isfile(joinpath(tmp, "progress.log"))
        @test occursin("stage one", read(joinpath(tmp, "progress.log"), String))

        # show_progress=false: writes neither.
        out_file2 = joinpath(tmp, "stdout2.txt")
        open(out_file2, "w") do io
            redirect_stdout(io) do
                EnzymeRates._progress(tmp, false, "silent line")
            end
        end
        @test !occursin("silent line", read(out_file2, String))
        @test !occursin("silent line", read(joinpath(tmp, "progress.log"), String))
    end

    # show_progress=false has no side effect: a non-existent save_dir is NOT
    # created (the early return precedes the mkpath).
    mktempdir() do tmp2
        ghost = joinpath(tmp2, "ghost_dir")
        EnzymeRates._progress(ghost, false, "no side effects")
        @test !isdir(ghost)
    end

    # _batch_summary reports the four reconciling buckets with the right
    # success/non-Success denominator.
    mech = first(EnzymeRates.init_mechanisms(@enzyme_reaction begin
        substrates: S[C]; products: P[C] end))
    row = (n_params=3, loss=0.5, mechanism_type="M", rate_equation="v",
           retcode="Success", error=missing, fitted_param_names=(:K,),
           fitted_param_values=(1.0,), eq_hash="abc", fit_inherited=false)
    e_succ = EnzymeRates.BatchEntry(mech, 3, 0.5, :Success, hash(:a), row)
    e_mt   = EnzymeRates.BatchEntry(mech, 3, 0.9, :MaxTime, hash(:b), row)
    f      = EnzymeRates.FitFailure(mech, "StackOverflowError: ")
    s = EnzymeRates._batch_summary([e_succ, e_mt], [f]; n_skipped=4, max_param_count=8)
    @test occursin("2 new fits + 0 inherited + 4 skipped (>8 params) + 1 errored", s)
    @test occursin("Success 50.0%", s)                 # 1 of 2 fitted
    @test occursin("non-Success retcode 50.0%", s)     # e_mt is :MaxTime
    @test !occursin("best loss", s)                    # best loss moved to its own line
end

@testset "_best_loss_line" begin
    line = EnzymeRates._best_loss_line(
        Dict(5 => 0.01751, 6 => 0.009316), Set([6]))
    @test occursin("best loss by n_params:", line)
    @test occursin("5:0.01751 ", line)          # unimproved: no star
    @test occursin("6:0.009316*", line)         # improved: starred
    @test occursin("(* improved)", line)

    quiet = EnzymeRates._best_loss_line(Dict(5 => 0.01751), Set{Int}())
    @test occursin("(no improvement)", quiet)
    @test !occursin("*", quiet)
end

@testset "_process_batch" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = DataFrame(
        S = [1.0, 2.0, 3.0, 4.0],
        P = [0.1, 0.2, 0.3, 0.4],
        Rate = [0.5, 0.8, 1.0, 1.1],
        group = [1, 1, 2, 2],
    )
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    ms = unique!(collect(EnzymeRates.init_mechanisms(rxn)))

    entries, failures = EnzymeRates._process_batch(ms, prob;
        optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=20, n_restarts=1, maxtime=1.0)
    @test !isempty(entries)
    @test isempty(failures)
    @test all(e -> e isa EnzymeRates.BatchEntry, entries)
    @test all(e -> e.retcode isa Symbol, entries)
    @test all(e -> e.n_params == length(e.row.fitted_param_names), entries)
    @test all(e -> occursin(r"^[0-9a-f]{16}$", e.row.eq_hash), entries)

    # cap filter: nothing over the cap is fit (and it is not a failure).
    capped_entries, capped_failures = EnzymeRates._process_batch(ms, prob;
        optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=0, n_restarts=1, maxtime=1.0)
    @test isempty(capped_entries)
    @test isempty(capped_failures)

    # config error (solver rejects an option) → every fit throws → all
    # failures, no entries; each failure carries a non-empty error string.
    fail_entries, fail_failures = EnzymeRates._process_batch(ms, prob;
        optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=20, n_restarts=1, maxtime=1.0,
        solver_kwargs=(; not_a_real_solver_option=1))
    @test isempty(fail_entries)
    @test !isempty(fail_failures)
    @test all(f -> f isa EnzymeRates.FitFailure, fail_failures)
    @test all(f -> !isempty(f.error), fail_failures)
end

@testset "_ingest! and cv pool" begin
    mk(n, loss, h) = EnzymeRates.BatchEntry(
        first(EnzymeRates.init_mechanisms(@enzyme_reaction begin
            substrates:S[C]; products:P[C] end)),
        n, loss, :Success, hash(h),
        (n_params=n, loss=loss, mechanism_type="M",
         rate_equation="v", retcode="Success", error=missing,
         fitted_param_names=(:K,), fitted_param_values=(1.0,),
         eq_hash=string(hash(h),base=16,pad=16), fit_inherited=false))
    frontier = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    cv_pool  = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    best     = Dict{Int,Float64}()
    # two distinct equations + one duplicate-eq with worse loss, n_cv=2
    EnzymeRates._ingest!(frontier, cv_pool, best,
        [mk(5,2.0,:a), mk(5,1.0,:b), mk(5,3.0,:a)]; n_cv_candidates=2)
    @test length(frontier[5]) == 3            # frontier keeps ALL
    @test best[5] == 1.0                       # running min
    @test length(cv_pool[5]) == 2              # bounded, distinct eq_hash
    # the kept :a entry is the lower-loss one (2.0, not 3.0);
    # BatchEntry.eq_hash is a UInt64, so compare against hash(:a), not hex:
    a = only(filter(e -> e.eq_hash == hash(:a), cv_pool[5]))
    @test a.loss == 2.0
    # n=0 must not panic on the empty pool (n_cv_candidates is public)
    @test EnzymeRates._offer_cv!(EnzymeRates.BatchEntry[], mk(5,1.0,:a), 0) ==
          EnzymeRates.BatchEntry[]
end

@testset "identify runs on a solver that rejects popsize" begin
    # identify_rate_equation must run end-to-end with only default
    # solver_kwargs on a solver that does not accept solver-specific
    # options — it injects no solver-specific option of its own.
    # (CMAEvolutionStrategy rejects unknown options such as popsize.)
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (group = ["G1", "G1", "G2", "G2"],
            Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0],
            P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    tmp = mktempdir()
    results = identify_rate_equation(prob;
        optimizer=CMAEvolutionStrategyOpt(),
        min_beam_width=1, loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        max_param_count=6, n_cv_candidates=1, n_restarts=1, maxtime=1.0,
        save_dir=tmp, show_progress=false)
    @test results isa IdentifyRateEquationResults
end

@testset "all-cap-skipped expansion batch is reported (M2)" begin
    # uni-uni base mechanism has 3 params; every child has 4. With
    # max_param_count=3 the base fits but the whole expansion batch is
    # cap-skipped — no rows, no CSV — so it must still emit a progress line.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0], P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    tmp = mktempdir()
    identify_rate_equation(prob;
        optimizer=CMAEvolutionStrategyOpt(),
        min_beam_width=1, loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        max_param_count=3, n_cv_candidates=1, n_restarts=1, maxtime=1.0,
        save_dir=tmp)
    log_text = read(joinpath(tmp, "progress.log"), String)
    @test occursin(r"all \d+ skipped \(>3 params\)", log_text)
    # The all-skip batch produced no rows, so no iteration CSV was written.
    @test !any(startswith(f, "equation_search_iteration_") for f in readdir(tmp))
end

@testset "loss_parsimony_threshold threads through identify_rate_equation" begin
    # An unknown keyword throws at the call boundary (see the removed-kwargs
    # test), so a clean end-to-end run with an explicit non-default value
    # proves the keyword is accepted and forwarded to the beam.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (group = ["G1", "G1", "G2", "G2"],
            Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0],
            P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    results = identify_rate_equation(prob;
        optimizer=CMAEvolutionStrategyOpt(),
        min_beam_width=1, loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        loss_parsimony_threshold=2.0,
        max_param_count=6, n_cv_candidates=1, n_restarts=1, maxtime=1.0,
        save_dir=mktempdir(), show_progress=false)
    @test results isa IdentifyRateEquationResults
end

@testset "removed kwargs error at the identify boundary" begin
    # popsize/verbose are no longer named kwargs and there is no catch-all,
    # so they are rejected immediately at the call boundary (before any
    # fitting or CSV write) — distinct from a solver-rejected solver_kwargs
    # option, which fails inside fitting.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (group = ["G1", "G1", "G2", "G2"],
            Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0],
            P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    @test_throws Exception identify_rate_equation(
        prob; popsize=200, optimizer=CMAEvolutionStrategyOpt(),
        n_restarts=1, maxtime=1.0, save_dir=mktempdir())
    @test_throws Exception identify_rate_equation(
        prob; verbose=-9, optimizer=CMAEvolutionStrategyOpt(),
        n_restarts=1, maxtime=1.0, save_dir=mktempdir())
end

# ── §2 fit-dedup by eq_hash ──────────────────────────────────────────────────
# A stub optimizer that counts `solve` invocations and returns a canned
# log-space optimum (`uval` for every coordinate), so a batch's fits can be
# counted exactly and the raw→rescale path exercised deterministically.
mutable struct _CountingStubOpt
    count::Int
    uval::Float64
    throwit::Bool
end
_CountingStubOpt(; uval=0.0, throwit=false) = _CountingStubOpt(0, uval, throwit)
Optimization.allowsbounds(::_CountingStubOpt) = true
function Optimization.SciMLBase.__solve(
        prob::Optimization.OptimizationProblem, opt::_CountingStubOpt; kwargs...)
    opt.count += 1
    opt.throwit && error("stub solver forced failure")
    u = fill(opt.uval, length(prob.u0))
    cache = DefaultOptimizationCache(prob.f, prob.p)
    build_solution(cache, opt, u, prob.f(u, prob.p); retcode = ReturnCode.Success)
end

# Two structurally-distinct bi-bi mechanisms that render the SAME reduced rate
# equation (eq_hash 78546b5f56b20e15): identical 7-edge topology, differing only
# in kinetic-group partitioning. Captured from the bi-bi enumeration (init ∪
# expand) and reconstructed from their EnzymeMechanism{Sig} type, so the test is
# self-contained and cheap (no full child enumeration, which would compile ~800
# distinct rate equations). The test re-verifies the collision at run time.
const _DEDUP_SIG1 =
    "EnzymeMechanism{(((((:Substrate, :A), ((:C, 1),)), ((:Substrate, :B), ((:N" *
    ", 1),)), ((:Product, :P), ((:C, 1),)), ((:Product, :Q), ((:N, 1),))), (), " *
    "(1,)), (((((), :E, ((), ())), (((:Substrate, :A),), :E, ((), ())), (:Subst" *
    "rate, :A), true), ((((:Product, :Q),), :E, ((), ())), (((:Substrate, :A), " *
    "(:Product, :Q)), :E, ((), ())), (:Substrate, :A), true)), ((((), :E, ((), " *
    "())), (((:Product, :Q),), :E, ((), ())), (:Product, :Q), true), ((((:Subst" *
    "rate, :A),), :E, ((), ())), (((:Substrate, :A), (:Product, :Q)), :E, ((), " *
    "())), (:Product, :Q), true)), (((((:Substrate, :A),), :E, ((), ())), (((:S" *
    "ubstrate, :A), (:Substrate, :B)), :E, ((), ())), (:Substrate, :B), true),)" *
    ", (((((:Substrate, :A), (:Substrate, :B)), :E, ((), ())), (((:Product, :P)" *
    ", (:Product, :Q)), :E, ((), ())), nothing, false),), (((((:Product, :Q),)," *
    " :E, ((), ())), (((:Product, :P), (:Product, :Q)), :E, ((), ())), (:Produc" *
    "t, :P), true),)))}"

const _DEDUP_SIG2 =
    "EnzymeMechanism{(((((:Substrate, :A), ((:C, 1),)), ((:Substrate, :B), ((:N" *
    ", 1),)), ((:Product, :P), ((:C, 1),)), ((:Product, :Q), ((:N, 1),))), (), " *
    "(1,)), (((((), :E, ((), ())), (((:Substrate, :A),), :E, ((), ())), (:Subst" *
    "rate, :A), true),), ((((), :E, ((), ())), (((:Product, :Q),), :E, ((), ())" *
    "), (:Product, :Q), true), ((((:Substrate, :A),), :E, ((), ())), (((:Substr" *
    "ate, :A), (:Product, :Q)), :E, ((), ())), (:Product, :Q), true)), (((((:Su" *
    "bstrate, :A),), :E, ((), ())), (((:Substrate, :A), (:Substrate, :B)), :E, " *
    "((), ())), (:Substrate, :B), true),), (((((:Substrate, :A), (:Substrate, :" *
    "B)), :E, ((), ())), (((:Product, :P), (:Product, :Q)), :E, ((), ())), noth" *
    "ing, false),), (((((:Product, :Q),), :E, ((), ())), (((:Substrate, :A), (:" *
    "Product, :Q)), :E, ((), ())), (:Substrate, :A), true),), (((((:Product, :Q" *
    "),), :E, ((), ())), (((:Product, :P), (:Product, :Q)), :E, ((), ())), (:Pr" *
    "oduct, :P), true),)))}"

@testset "fit-dedup by eq_hash in _process_batch" begin
    recon(sig) = EnzymeRates.Mechanism(Core.eval(EnzymeRates, Meta.parse(sig))())
    m1 = recon(_DEDUP_SIG1)
    m2 = recon(_DEDUP_SIG2)
    em1 = EnzymeRates.compile_mechanism(m1)
    em2 = EnzymeRates.compile_mechanism(m2)
    key = EnzymeRates._rate_eq_dedup_key(rate_equation_string(em1))
    # Preconditions: distinct structure, identical eq_hash + fitted-param set
    # (same names AND order — the reuse maps params by name directly).
    @test m1 != m2
    @test key == EnzymeRates._rate_eq_dedup_key(rate_equation_string(em2))
    @test EnzymeRates.fitted_params(em1) == EnzymeRates.fitted_params(em2)

    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            A = [1.0, 2.0, 1.0, 2.0], B = [0.5, 0.5, 1.0, 1.0],
            P = [0.1, 0.2, 0.1, 0.2], Q = [0.3, 0.3, 0.4, 0.4])
    prob = IdentifyRateEquationProblem(EnzymeRates.reaction(m1), data; Keq=2.0)
    pair = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1, m2]

    memo = Dict{UInt64, NamedTuple}()
    opt = _CountingStubOpt(; uval = log(5.0))
    entries, failures = EnzymeRates._process_batch(pair, prob;
        optimizer=opt, max_param_count=20, n_restarts=1, maxtime=1.0, memo)

    # The shared equation is fit exactly ONCE (n_restarts=1 → one solve).
    @test opt.count == 1
    @test length(entries) == 2
    @test isempty(failures)
    # loss + retcode are equation properties → eq_hash-invariant → identical.
    @test entries[1].loss == entries[2].loss
    @test entries[1].retcode === entries[2].retcode === :Success
    # Representative fit first (false); the duplicate is inherited (true).
    @test [e.row.fit_inherited for e in entries] == [false, true]

    # The memo stores the fit; every row sharing this eq_hash copies it verbatim,
    # so both rows carry identical params — same equation ⟹ same fit + rescaling.
    fit = memo[key]
    fkeys = EnzymeRates.fitted_params(em1)
    for e in entries
        @test e.row.fitted_param_values == Tuple(fit.params[k] for k in fkeys)
    end
    @test entries[1].row.fitted_param_values == entries[2].row.fitted_param_values
    # scale_k_to_kcat=1.0 anchored kcat: the copied params are the rescaled fit,
    # not the raw 5.0 the stub optimizer returned.
    @test !all(v -> v ≈ 5.0, entries[1].row.fitted_param_values)

    # Cross-batch memo hit: a later batch with the same eq_hash refits NOTHING.
    single = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1]
    reused, _ = EnzymeRates._process_batch(single, prob;
        optimizer=opt, max_param_count=20, n_restarts=1, maxtime=1.0, memo)
    @test opt.count == 1                                # no new solve
    @test [e.row.fit_inherited for e in reused] == [true]

    # A representative whose fit throws fails ALL its duplicates
    # (all-or-nothing per equation).
    opt_bad = _CountingStubOpt(; throwit=true)
    bad_entries, bad_failures = EnzymeRates._process_batch(pair, prob;
        optimizer=opt_bad, max_param_count=20, n_restarts=1, maxtime=1.0,
        memo = Dict{UInt64, NamedTuple}())
    @test isempty(bad_entries)
    @test length(bad_failures) == 2
    @test all(f -> f isa EnzymeRates.FitFailure, bad_failures)
end


# Confirmed LDH renaming-dup pair (same graph, tied kinetic-group split): merged
# form (8 groups) vs split form (9 groups). Currently different eq_hash; the
# pre-fit canonical-partition merge must collapse them.
const _CANON_SIG_MERGED =
    "EnzymeMechanism{(((((:Product, :Lactate), ((:C, 3), (:H, 6), (:O, " *
    "3))), ((:Product, :NAD), ((:C, 21), (:H, 27), (:N, 7), (:O, 14), (" *
    ":P, 2))), ((:Substrate, :NADH), ((:C, 21), (:H, 29), (:N, 7), (:O," *
    " 14), (:P, 2))), ((:Substrate, :Pyruvate), ((:C, 3), (:H, 4), (:O," *
    " 3)))), (), (4,)), (((((), :E, ((), ())), (((:Product, :Lactate),)" *
    ", :E, ((), ())), (:Product, :Lactate), true), ((((:Product, :NAD)," *
    "), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, (" *
    "(), ())), (:Product, :Lactate), true), ((((:Substrate, :NADH),), :" *
    "E, ((), ())), (((:Product, :Lactate), (:Substrate, :NADH)), :E, ((" *
    "), ())), (:Product, :Lactate), true)), ((((), :E, ((), ())), (((:P" *
    "roduct, :NAD),), :E, ((), ())), (:Product, :NAD), true), ((((:Prod" *
    "uct, :Lactate),), :E, ((), ())), (((:Product, :Lactate), (:Product" *
    ", :NAD)), :E, ((), ())), (:Product, :NAD), true)), ((((), :E, (()," *
    " ())), (((:Substrate, :NADH),), :E, ((), ())), (:Substrate, :NADH)" *
    ", true), ((((:Product, :Lactate),), :E, ((), ())), (((:Product, :L" *
    "actate), (:Substrate, :NADH)), :E, ((), ())), (:Substrate, :NADH)," *
    " true)), ((((), :E, ((), ())), (((:Substrate, :Pyruvate),), :E, ((" *
    "), ())), (:Substrate, :Pyruvate), true),), (((((:Product, :NAD),)," *
    " :E, ((), ())), (((:Product, :NAD), (:Substrate, :Pyruvate)), :E, " *
    "((), ())), (:Substrate, :Pyruvate), true), ((((:Substrate, :NADH)," *
    "), :E, ((), ())), (((:Substrate, :NADH), (:Substrate, :Pyruvate))," *
    " :E, ((), ())), (:Substrate, :Pyruvate), true)), (((((:Substrate, " *
    ":NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (((:Product, :Lac" *
    "tate), (:Product, :NAD)), :E, ((), ())), nothing, false),), (((((:" *
    "Substrate, :Pyruvate),), :E, ((), ())), (((:Substrate, :NADH), (:S" *
    "ubstrate, :Pyruvate)), :E, ((), ())), (:Substrate, :NADH), true),)" *
    ", (((((:Substrate, :Pyruvate),), :E, ((), ())), (((:Product, :NAD)" *
    ", (:Substrate, :Pyruvate)), :E, ((), ())), (:Product, :NAD), true)" *
    ",)))}"
const _CANON_SIG_SPLIT =
    "EnzymeMechanism{(((((:Product, :Lactate), ((:C, 3), (:H, 6), (:O, " *
    "3))), ((:Product, :NAD), ((:C, 21), (:H, 27), (:N, 7), (:O, 14), (" *
    ":P, 2))), ((:Substrate, :NADH), ((:C, 21), (:H, 29), (:N, 7), (:O," *
    " 14), (:P, 2))), ((:Substrate, :Pyruvate), ((:C, 3), (:H, 4), (:O," *
    " 3)))), (), (4,)), (((((), :E, ((), ())), (((:Product, :Lactate),)" *
    ", :E, ((), ())), (:Product, :Lactate), true), ((((:Product, :NAD)," *
    "), :E, ((), ())), (((:Product, :Lactate), (:Product, :NAD)), :E, (" *
    "(), ())), (:Product, :Lactate), true)), ((((), :E, ((), ())), (((:" *
    "Product, :NAD),), :E, ((), ())), (:Product, :NAD), true), ((((:Pro" *
    "duct, :Lactate),), :E, ((), ())), (((:Product, :Lactate), (:Produc" *
    "t, :NAD)), :E, ((), ())), (:Product, :NAD), true)), ((((), :E, (()" *
    ", ())), (((:Substrate, :NADH),), :E, ((), ())), (:Substrate, :NADH" *
    "), true), ((((:Product, :Lactate),), :E, ((), ())), (((:Product, :" *
    "Lactate), (:Substrate, :NADH)), :E, ((), ())), (:Substrate, :NADH)" *
    ", true)), ((((), :E, ((), ())), (((:Substrate, :Pyruvate),), :E, (" *
    "(), ())), (:Substrate, :Pyruvate), true),), (((((:Product, :NAD),)" *
    ", :E, ((), ())), (((:Product, :NAD), (:Substrate, :Pyruvate)), :E," *
    " ((), ())), (:Substrate, :Pyruvate), true), ((((:Substrate, :NADH)" *
    ",), :E, ((), ())), (((:Substrate, :NADH), (:Substrate, :Pyruvate))" *
    ", :E, ((), ())), (:Substrate, :Pyruvate), true)), (((((:Substrate," *
    " :NADH),), :E, ((), ())), (((:Product, :Lactate), (:Substrate, :NA" *
    "DH)), :E, ((), ())), (:Product, :Lactate), true),), (((((:Substrat" *
    "e, :NADH), (:Substrate, :Pyruvate)), :E, ((), ())), (((:Product, :" *
    "Lactate), (:Product, :NAD)), :E, ((), ())), nothing, false),), (((" *
    "((:Substrate, :Pyruvate),), :E, ((), ())), (((:Substrate, :NADH), " *
    "(:Substrate, :Pyruvate)), :E, ((), ())), (:Substrate, :NADH), true" *
    "),), (((((:Substrate, :Pyruvate),), :E, ((), ())), (((:Product, :N" *
    "AD), (:Substrate, :Pyruvate)), :E, ((), ())), (:Product, :NAD), tr" *
    "ue),)))}"

@testset "canonical-partition dedup collapses renaming-dups in _process_batch" begin
    recon(sig) = EnzymeRates.Mechanism(Core.eval(EnzymeRates, Meta.parse(sig))())
    m1 = recon(_CANON_SIG_MERGED)   # merged, 8 kinetic groups
    m2 = recon(_CANON_SIG_SPLIT)    # split, 9 groups — same graph, Wegscheider-tied
    em1 = EnzymeRates.compile_mechanism(m1)
    em2 = EnzymeRates.compile_mechanism(m2)
    # Precondition: same rate function, but the RAW dedup key currently DIFFERS —
    # the renaming-dup the pre-fit canonicalization must collapse.
    @test m1 != m2
    @test EnzymeRates._rate_eq_dedup_key(rate_equation_string(em1)) !=
          EnzymeRates._rate_eq_dedup_key(rate_equation_string(em2))

    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            NADH = [1.0, 2.0, 1.0, 2.0], Pyruvate = [0.5, 0.5, 1.0, 1.0],
            Lactate = [0.1, 0.2, 0.1, 0.2], NAD = [0.3, 0.3, 0.4, 0.4])
    prob = IdentifyRateEquationProblem(EnzymeRates.reaction(m1), data; Keq=2.0)
    pair = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1, m2]

    memo = Dict{UInt64, NamedTuple}()
    opt = _CountingStubOpt(; uval = log(5.0))
    entries, failures = EnzymeRates._process_batch(pair, prob;
        optimizer=opt, max_param_count=20, n_restarts=1, maxtime=1.0, memo)

    # Canonicalized to one form → fit exactly ONCE; both rows share it.
    @test opt.count == 1
    @test length(entries) == 2
    @test isempty(failures)
    @test entries[1].row.eq_hash == entries[2].row.eq_hash
    @test [e.row.fit_inherited for e in entries] == [false, true]
end

@testset "LOOCV eq_hash-uniqueness guard (§4)" begin
    # _cv_model_selection dedups candidates by eq_hash per n_params bucket before
    # LOOCV: same-equation twins collapse to ONE candidate (lowest loss kept), so
    # folds are never wasted on, or biased by, textually identical equations.
    recon(sig) = EnzymeRates.Mechanism(Core.eval(EnzymeRates, Meta.parse(sig))())
    m1 = recon(_DEDUP_SIG1); m2 = recon(_DEDUP_SIG2)   # distinct structure, same eq_hash
    em1 = EnzymeRates.compile_mechanism(m1)
    fkeys = EnzymeRates.fitted_params(em1)
    h = string(EnzymeRates._rate_eq_dedup_key(rate_equation_string(em1)), base=16, pad=16)
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            A = [1.0, 2.0, 1.0, 2.0], B = [0.5, 0.5, 1.0, 1.0],
            P = [0.1, 0.2, 0.1, 0.2], Q = [0.3, 0.3, 0.4, 0.4])
    prob = IdentifyRateEquationProblem(EnzymeRates.reaction(m1), data; Keq=2.0)
    mechs = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1, m2]
    mkrow(loss) = (n_params=length(fkeys), loss=loss, mechanism_type="M",
        rate_equation="v", retcode="Success", error=missing,
        fitted_param_names=Tuple(fkeys), fitted_param_values=Tuple(fill(1.0, length(fkeys))),
        eq_hash=h, fit_inherited=false)
    df = EnzymeRates._rows_to_dataframe([mkrow(0.5), mkrow(0.2)])  # m1 loss .5, m2 loss .2
    res = EnzymeRates._cv_model_selection(mechs, df, prob;
        n_cv_candidates=5, optimizer=_CountingStubOpt(; uval=log(5.0)),
        se_threshold=1.0, perm_p_threshold=1.0, save_dir=mktempdir(),
        show_progress=false, n_restarts=1, maxtime=1.0)
    # The two same-eq_hash twins collapsed to a single LOOCV candidate…
    @test nrow(res.cv_results) == 1
    @test res.cv_results.eq_hash[1] == h
    # …and the lower-loss twin (0.2) was the one kept.
    @test res.cv_results.loss[1] == 0.2

    # `_offer_cv!` likewise keeps at most one entry per eq_hash: a repeat hash
    # updates its own slot to the lower loss, never consuming a second.
    mech = first(EnzymeRates.init_mechanisms(@enzyme_reaction begin
        substrates: S[C]; products: P[C] end))
    mkentry(loss, h) = EnzymeRates.BatchEntry(mech, 5, loss, :Success, h,
        (n_params=5, loss=loss, mechanism_type="M", rate_equation="v",
         retcode="Success", error=missing, fitted_param_names=(:K,),
         fitted_param_values=(1.0,), eq_hash=string(h, base=16, pad=16),
         fit_inherited=false))
    pool = EnzymeRates.BatchEntry[]
    for (loss, h) in [(1.0, UInt64(1)), (0.5, UInt64(1)), (2.0, UInt64(2))]
        EnzymeRates._offer_cv!(pool, mkentry(loss, h), 5)
    end
    @test allunique([e.eq_hash for e in pool])
    @test length(pool) == 2
    @test only(filter(e -> e.eq_hash == UInt64(1), pool)).loss == 0.5
end

@testset "_cv_model_selection flatten reproduces serial LOOCV" begin
    # Deterministic stub optimizer → identical fits whether folds run serially
    # or across the flattened (candidate, fold) grid.
    recon(sig) = EnzymeRates.Mechanism(Core.eval(EnzymeRates, Meta.parse(sig))())
    m1 = recon(_DEDUP_SIG1)
    em1 = EnzymeRates.compile_mechanism(m1)
    fkeys = EnzymeRates.fitted_params(em1)
    h = string(EnzymeRates._rate_eq_dedup_key(rate_equation_string(em1)),
               base=16, pad=16)
    data = (group = ["G1", "G1", "G2", "G2", "G3", "G3"],
            Rate = [0.5, 0.8, 1.0, 1.1, 0.9, 1.2],
            A = [1.0, 2.0, 1.0, 2.0, 1.5, 2.5], B = [0.5, 0.5, 1.0, 1.0, 0.7, 0.7],
            P = [0.1, 0.2, 0.1, 0.2, 0.15, 0.25], Q = [0.3, 0.3, 0.4, 0.4, 0.35, 0.35])
    prob = IdentifyRateEquationProblem(EnzymeRates.reaction(m1), data; Keq=2.0)
    mechs = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1]
    mkrow(loss) = (n_params=length(fkeys), loss=loss, mechanism_type="M",
        rate_equation="v", retcode="Success", error=missing,
        fitted_param_names=Tuple(fkeys),
        fitted_param_values=Tuple(fill(1.0, length(fkeys))),
        eq_hash=h, fit_inherited=false)
    df = EnzymeRates._rows_to_dataframe([mkrow(0.5)])

    res = EnzymeRates._cv_model_selection(mechs, df, prob;
        n_cv_candidates=5, optimizer=_CountingStubOpt(; uval=log(5.0)),
        se_threshold=1.0, perm_p_threshold=1.0, save_dir=mktempdir(),
        show_progress=false, n_restarts=1, maxtime=1.0)

    groups = unique(prob.data.group)
    flat = [res.cv_results[1, Symbol("cv_fold_$g")] for g in groups]
    m1c = EnzymeRates.compile_mechanism(m1)
    serial = [EnzymeRates._cv_fold_loss(m1c, prob, g;
        optimizer=_CountingStubOpt(; uval=log(5.0)), n_restarts=1, maxtime=1.0)
        for g in groups]
    @test flat == serial
end

@testset "_scatter_fold_scores places each fold at (candidate, group)" begin
    groups = ["G1", "G2", "G3"]
    # 2 candidates, distinct per-candidate scores, deliberately shuffled
    flat = [(1, "G2", 0.12), (2, "G1", 0.20), (1, "G1", 0.10),
            (2, "G3", 0.23), (1, "G3", 0.13), (2, "G2", 0.21)]
    out = EnzymeRates._scatter_fold_scores(flat, 2, groups)
    @test out[1] == [0.10, 0.12, 0.13]
    @test out[2] == [0.20, 0.21, 0.23]
end
