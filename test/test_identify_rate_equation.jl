# ABOUTME: Tests for identify_rate_equation pipeline.
# ABOUTME: Covers construction, helpers, and full pipeline
# ABOUTME: with mechanism recovery including allosteric path.

using DataFrames
using CSV
using Random
using Statistics
using OptimizationPyCMA

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

    pycma_opt = PyCMAOpt()

    # ── Construction validation ──────────────────────
    @testset "construction" begin
        prob = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val)
        @test prob.reaction === test_rxn
        @test prob.Keq == Keq_val
        @test length(
            unique(prob.data.group)) == 5

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
            fitted_param_names = (:a, :b),
            fitted_param_values = (1.0, 2.0),
            eq_hash = "0123456789abcdef",
            fit_inherited_from_estimate = missing,
        )]
        df = EnzymeRates._rows_to_dataframe(
            rows)
        @test nrow(df) == 1
        @test "a" in names(df)
        @test "b" in names(df)
        @test "eq_hash" in names(df)
        @test "fit_inherited_from_estimate" in names(df)

        # Empty rows
        df2 = EnzymeRates._rows_to_dataframe(
            NamedTuple[])
        @test nrow(df2) == 0
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
        pmap_function=map,
        optimizer=pycma_opt,
        n_restarts=1, maxtime=1.0)

    @testset "mechanism recovery" begin
        # The best mechanism should have the same
        # number of independent parameters as the
        # generating mechanism (or fewer)
        best_np = length(
            parameters(results.best, Reduced))
        gen_np = length(
            parameters(test_mechanism, Reduced))
        @test best_np <= gen_np

        # The best mechanism should fit the
        # noiseless data with near-zero loss
        fp_best = FittingProblem(
            results.best, test_data;
            Keq=Keq_val)
        fit_best = fit_rate_equation(
            fp_best, pycma_opt;
            n_restarts=3, maxtime=10.0,
            popsize=200)
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

    @testset "CSV output" begin
        csv_files = filter(
            f -> endswith(f, ".csv"),
            readdir(save_dir))
        @test length(csv_files) > 0
        # Per-level filename uses the estimate-level naming.
        for fname in csv_files
            @test startswith(fname, "params_estimate_")
        end

        first_csv = CSV.read(
            joinpath(
                save_dir, csv_files[1]),
            DataFrame)
        @test "n_params" in names(first_csv)
        @test "loss" in names(first_csv)
        @test "mechanism_type" in names(
            first_csv)
        @test nrow(first_csv) > 0

        # eq_hash column: 16-char hex, no missing values.
        for fname in csv_files
            df_file = CSV.read(
                joinpath(save_dir, fname), DataFrame)
            @test "eq_hash" in names(df_file)
            @test all(.!ismissing.(df_file.eq_hash))
            @test all(length.(df_file.eq_hash) .== 16)
        end

        # Cross-level fit-inheritance chain: rows with non-
        # missing `fit_inherited_from_estimate` must point to
        # an existing level whose CSV contains a row with the
        # same `eq_hash`.
        all_rows_by_level = Dict{Int, DataFrame}()
        for fname in csv_files
            est = parse(Int, replace(fname,
                "params_estimate_" => "", ".csv" => ""))
            all_rows_by_level[est] = CSV.read(
                joinpath(save_dir, fname), DataFrame)
        end
        for (_, df_lvl) in all_rows_by_level
            for row in eachrow(df_lvl)
                ismissing(row.fit_inherited_from_estimate) &&
                    continue
                src = row.fit_inherited_from_estimate
                @test haskey(all_rows_by_level, src)
                @test row.eq_hash in
                    all_rows_by_level[src].eq_hash
            end
        end
        # NOTE: a true positive-path test of the cross-level fit
        # cache (asserting at least one inherited row is produced)
        # would require either a richer fixture (e.g., a bi-bi
        # mechanism prone to Haldane-collapse hash hits) or a
        # `_beam_search` unit test with a recording fit-wrapper.
        # The greedy-beam smoke settings here (min_beam_width=1)
        # don't reliably produce inherited rows. A recording fit-wrapper
        # would make that positive path deterministic; the loop above
        # catches invalid `fit_inherited_from_estimate` references when present.
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
                pmap_function=map,
                optimizer=pycma_opt,
                n_restarts=1, maxtime=1.0))
    end

    @testset "_loocv returns per-fold scores, floored at eps" begin
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

        scores = EnzymeRates._loocv(
            m, prob;
            optimizer=PyCMAOpt(),
            n_restarts=2, maxtime=2.0,
            maxiters=500, popsize=40, verbose=-9)

        @test scores isa Vector{Float64}
        # Require the success path: fitting MUST converge on
        # this trivial uni-uni fixture in 2s × 2 restarts. A
        # length-0 (full failure) result would let the per-fold
        # eps-floor + isfinite assertions below pass vacuously.
        @test length(scores) == 3
        @test all(s -> s >= eps(Float64), scores)
        @test all(isfinite, scores)
    end

end

@testset "save_level_csv uses estimate-level filename" begin
    mktempdir() do tmp
        rows = [(n_params=3, loss=1.0,
                 mechanism_type="m1", rate_equation="eq1",
                 fitted_param_names=(:K1, :K2, :k3f),
                 fitted_param_values=(1.0, 2.0, 3.0),
                 eq_hash="0123456789abcdef",
                 fit_inherited_from_estimate=missing)]
        # Caller passes the estimate-level pc (e.g., 5) — could
        # diverge from the row's actual n_params=3 due to Haldane
        # reduction. Filename must reflect the estimate.
        EnzymeRates._save_level_csv(tmp, rows, 5)
        @test isfile(joinpath(tmp, "params_estimate_5.csv"))
        @test !isfile(joinpath(tmp, "params_5.csv"))
        df = CSV.read(joinpath(tmp, "params_estimate_5.csv"),
                      DataFrame)
        @test df.n_params == [3]
    end
end

@testset "csv writers" begin
    rows = [(
        n_params = 5, loss = 1.0, mechanism_type = "M",
        rate_equation = "v = 1", fitted_param_names = (:K_a,),
        fitted_param_values = (2.0,), eq_hash = "abc",
        fit_inherited_from_estimate = missing,
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

@testset "beam_fraction kwarg removed: passing it errors" begin
    # The removed `beam_fraction` kwarg was replaced by
    # `loss_rel_threshold` + `loss_abs_threshold` + `min_beam_width`.
    # No alias / deprecation shim — passing it must error.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (group = ["G1", "G1", "G2", "G2"],
            Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0],
            P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    @test_throws MethodError identify_rate_equation(
        prob; beam_fraction=0.5,
        optimizer=PyCMAOpt(),
        n_restarts=1, maxtime=1.0)
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
    # All rows have empty fold scores → error.
    cv_df_empty = DataFrame(
        n_params       = [3, 5],
        cv_score       = [Inf, Inf],
        loss           = [0.0, 0.0],
        cv_fold_scores = [Float64[], Float64[]],
    )
    @test_throws ErrorException EnzymeRates._select_best_n_params(
        cv_df_empty)

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

    # Partial bucket failure: one row in bucket-3 failed (empty fold
    # scores), another row valid → bucket retained via rep selection.
    cv_df_partial = DataFrame(
        n_params       = [3, 3, 5],
        cv_score       = [Inf, 0.115, 0.115],
        loss           = [0.0, 0.0, 0.0],
        cv_fold_scores = [
            Float64[],
            exp.([0.10, 0.12, 0.11, 0.13]),
            exp.([0.10, 0.12, 0.11, 0.13]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df_partial).best_n == 3
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

@testset "canonical-hash partition stability" begin
    # Representative reactions exercising the canonicalizer's edge cases:
    # - uni_uni: trivial structural equivalence
    # - bi_bi:   substituted-into-v ties across multiple kinetic groups
    # ter-ter intentionally omitted — `rate_equation` compilation is
    # extremely slow for mechanisms with >~30 enzyme forms (CLAUDE.md
    # "Known Issues"), and the canonical hasher invokes that path per
    # candidate. The bi-bi enumeration already covers every structural
    # symmetry the canonicalizer collapses.
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
    # equations the init-level enumeration produces. After the catalytic-
    # topology connectivity fix, the 69 bi_bi init mechanisms are all
    # structurally distinct AND each yields a distinct `rate_equation_string`,
    # and the canonical hasher produces exactly 69 classes (verified directly:
    # distinct hashes == distinct rate-equation strings == 69, i.e. zero over-
    # and zero under-collapse). The fix removed the dangling-form / binding-
    # order rapid-equilibrium twins that previously collapsed the (buggy) 77
    # mechanisms to 21 classes: clean topologies have distinct enzyme-form
    # sets, hence distinct rate equations.
    # If these counts change in a future commit, the canonical hasher's
    # equivalence classes (or the enumeration) have shifted — investigate.
    expected_n_classes = Dict(
        "uni_uni" => 1,
        "bi_bi"   => 69,
    )

    for (label, reaction) in test_reactions
        # init_mechanisms only — skip expand_mechanisms. The init level
        # already produces multiple structurally-equivalent variants
        # (mirror-step orderings, kinetic-group renumberings) that
        # exercise the canonicalizer's collapse rules. expand_mechanisms
        # adds variants at higher param counts whose canonicalizer
        # behavior is the same modulo size, at exponential compile cost.
        all_mechs = EnzymeRates.init_mechanisms(reaction)

        new_buckets = Dict{UInt64, Vector{Int}}()
        for (i, m) in enumerate(all_mechs)
            em = EnzymeRates.compile_mechanism(m)
            h = EnzymeRates._canonical_rate_eq_hash(em)
            push!(get!(new_buckets, h, Int[]), i)
            # Determinism: same input, same hash across invocations.
            @test EnzymeRates._canonical_rate_eq_hash(em) === h
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
