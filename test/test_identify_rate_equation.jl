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

# The uni-uni reaction S[C] → P[C], and a 4-row, two-group problem on it (Keq = 10) whose
# data is a NamedTuple with String groups or a DataFrame with Int groups.
const _testhelper_uni_rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
end

# The seed of that reaction: rapid-equilibrium bindings around a steady-state
# isomerization.
const _testhelper_uni_seed = EnzymeRates.Mechanism(@enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end)

# A BatchEntry of the uni seed with rate-law key `key`; selection reads only the key, the
# loss and the count.
_testhelper_entry(key, loss, n = 3) =
    EnzymeRates.BatchEntry(_testhelper_uni_seed, n, loss, :Success, key, key, (;))

function _testhelper_uni_prob(data_form)
    columns = (S = [1.0, 2.0, 3.0, 4.0], P = [0.1, 0.2, 0.3, 0.4],
               Rate = [0.5, 0.8, 1.0, 1.1])
    data = data_form === DataFrame ?
           DataFrame(merge(columns, (group = [1, 1, 2, 2],))) :
           merge((group = ["G1", "G1", "G2", "G2"],), columns)
    return IdentifyRateEquationProblem(_testhelper_uni_rxn, data; Keq=10.0)
end

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

    # The constrained allosteric mechanism: K-type allosteric
    # with S, P only in R-state (`:OnlyA` binding groups) and R
    # only in T-state (`:OnlyI` ligand tag).
    test_mechanism = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        steps: begin
            E + S ⇌ E(S)     :: OnlyA
            E(S) <--> E(P)   :: NonequalAI
            E(P) ⇌ E + P     :: OnlyA
        end
        regulatory_site(1): R::OnlyI
    end

    Keq_val = 2.0
    # 5 fitted params: K_A_EP_to_E_P, K_A_ES_to_E_S, k_A_ES_to_EP, K_I_Rreg, L
    true_params = (
        K_A_EP_to_E_P = 1.0, K_A_ES_to_E_S = 0.5, k_A_ES_to_EP = 5.0,
        K_I_Rreg = 2.0, L = 0.1,
        Keq = Keq_val, E_total = 1.0)

    function _testhelper_make_test_data(
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
    test_data = _testhelper_make_test_data(
        test_mechanism, true_params)

    cmaes_opt = CMAEvolutionStrategyOpt()

    # ── Construction validation ──────────────────────
    @testset "construction" begin
        prob = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val)
        @test prob.reaction === test_rxn
        @test typeof(prob) === IdentifyRateEquationProblem{typeof(prob.data)}
        @test prob.Keq == Keq_val
        @test length(
            unique(prob.data.group)) == 5

        # absolute and peak_range: relative mode and the default range unless set.
        @test prob.absolute === false
        @test prob.peak_range === (1e-12, 1e3)
        prob_abs = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, absolute=true, peak_range=(1e-9, 1e2))
        @test prob_abs.absolute === true && prob_abs.peak_range === (1e-9, 1e2)
        @test_throws MethodError IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, scale_k_to_kcat=1.0)
        # An integer Keq converts to the Float64 field.
        @test IdentifyRateEquationProblem(test_rxn, test_data; Keq=10).Keq === 10.0

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

        # A non-finite or missing rate, a concentration that is not a finite number
        # ≥ 0, a Keq that is not positive and finite, or a peak_range out of order: each
        # error names the column (with the row), Keq or peak_range.
        two_groups = (group = ["G1", "G2"], Rate = [1.0, 2.0],
                      S = [1.0, 1.0], P = [0.1, 0.1], R = [1.0, 1.0])
        for (bad, shown) in ((NaN, "NaN"), (Inf, "Inf"), (-Inf, "-Inf"),
                             (missing, "missing"))
            @test_throws(
                ErrorException("Rate at row 2 must be a finite number; got $shown"),
                IdentifyRateEquationProblem(
                    test_rxn, merge(two_groups, (Rate = [1.0, bad],)); Keq=1.0))
        end
        for (bad, shown) in ((NaN, "NaN"), (Inf, "Inf"), (-1.0, "-1.0"),
                             (missing, "missing"), ("1.0", "\"1.0\""))
            @test_throws(
                ErrorException(
                    "Concentration S at row 2 must be a finite number ≥ 0; got $shown"),
                IdentifyRateEquationProblem(
                    test_rxn, merge(two_groups, (S = [1.0, bad],)); Keq=1.0))
        end
        @test IdentifyRateEquationProblem(
            test_rxn, merge(two_groups, (S = [0.0, 1.0],)); Keq=1.0) isa
              IdentifyRateEquationProblem
        for (bad, shown) in ((0, "0"), (-1, "-1"), (Inf, "Inf"), (NaN, "NaN"))
            @test_throws(ErrorException("Keq must be positive and finite; got $shown"),
                IdentifyRateEquationProblem(test_rxn, two_groups; Keq=bad))
        end
        @test_throws(
            ErrorException("peak_range must be a tuple (lo, hi) with 0 < lo < hi < Inf; " *
                           "got (1.0, 0.5)"),
            IdentifyRateEquationProblem(test_rxn, two_groups; Keq=1.0,
                                        peak_range=(1.0, 0.5)))

        # A metabolite named after a required column would read that column, here the
        # group numbers, as its concentrations.
        group_named_rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: group
        end
        @test_throws(
            ErrorException("Metabolite group has the name of the required group " *
                           "column; rename the metabolite"),
            IdentifyRateEquationProblem(group_named_rxn,
                (group = [1, 2], Rate = [1.0, 2.0], S = [1.0, 1.0], P = [0.1, 0.1]);
                Keq=1.0))

        # Every search fits allosteric mechanisms, whose rate equations take the
        # conformational constant L, so a metabolite named L would fail every allosteric
        # candidate; the reaction refuses it up front.
        @test_throws(
            ErrorException("EnzymeReaction: L is a reserved name; E_total, Keq, L, N, " *
                           "params, concs and v are kept for a rate equation's constants " *
                           "and local names"),
            @enzyme_reaction begin
                substrates: S[C]
                products: L[C]
            end)

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
            mechanism = "test",
            mechanism_id = "0123456789abcdef",
            parent_id = "fedcba9876543210",
            rate_equation = "v = ...",
            retcode = "Success",
            error = missing,
            params = (a = 1.0, b = 2.0),
            key = repeat("0123456789abcdef", 4),
            fit_source_id = "0123456789abcdef",
            L_phantom = false,
            fit_inherited = false,
            peak_at_range_edge = true,
        )]
        df = EnzymeRates._rows_to_dataframe(
            rows)
        @test nrow(df) == 1
        # The schema: the fit's columns, then one column per fitted parameter.
        @test names(df) == ["n_params", "parent_n_params", "loss", "mechanism",
                            "mechanism_id", "parent_id", "rate_equation", "retcode",
                            "error", "key", "fit_inherited", "fit_source_id",
                            "L_phantom", "peak_at_range_edge", "a", "b"]
        @test df.key == [repeat("0123456789abcdef", 4)]
        @test df.fit_source_id == ["0123456789abcdef"]
        @test df.L_phantom == [false]
        @test (df.mechanism, df.mechanism_id, df.parent_id) ==
              (["test"], ["0123456789abcdef"], ["fedcba9876543210"])
        @test "a" in names(df)
        @test "b" in names(df)
        @test "key" in names(df)
        @test "retcode" in names(df)
        @test "error" in names(df)
        @test "fit_inherited" in names(df)
        @test df.fit_inherited == [false]
        @test df.peak_at_range_edge == [true]
        @test !("fit_inherited_from_estimate" in names(df))

        # Empty rows
        df2 = EnzymeRates._rows_to_dataframe(
            NamedTuple[])
        @test nrow(df2) == 0

        # A failure row (every fit field missing) beside a fitted row.
        rows_with_failure = [
            (n_params = 3, loss = 0.5, mechanism = "M", mechanism_id = "0000000000000001",
             rate_equation = "v = ...", retcode = "Success", error = missing,
             params = (a = 1.0,), key = "0123456789abcdef",
             fit_source_id = "0000000000000001", fit_inherited = false, L_phantom = false,
             peak_at_range_edge = false),
            (n_params = missing, loss = missing, mechanism = "M",
             mechanism_id = "0000000000000001", rate_equation = missing, retcode = missing,
             error = "StackOverflowError: ", params = (;), key = missing,
             fit_source_id = missing, fit_inherited = missing, L_phantom = missing,
             peak_at_range_edge = missing),
        ]
        df_with_failure = EnzymeRates._rows_to_dataframe(rows_with_failure)
        @test nrow(df_with_failure) == 2
        @test ismissing(df_with_failure.loss[2])
        @test df_with_failure.error[2] == "StackOverflowError: "
        @test ismissing(df_with_failure.retcode[2])
        @test ismissing(df_with_failure.key[2])
        @test ismissing(df_with_failure.fit_source_id[2])
        @test ismissing(df_with_failure.L_phantom[2])
        @test df_with_failure.fit_inherited[1] == false
        @test ismissing(df_with_failure.fit_inherited[2])
        @test df_with_failure.peak_at_range_edge[1] == false
        @test ismissing(df_with_failure.peak_at_range_edge[2])
        @test all(ismissing, df_with_failure.parent_id)
        # A row without the field records missing.
        @test ismissing(EnzymeRates._rows_to_dataframe([(n_params = 1, loss = 0.1,
            mechanism = "M", mechanism_id = "0", rate_equation = "v", retcode = "Success",
            error = missing, params = (;), key = "0", fit_source_id = "0",
            L_phantom = false, fit_inherited = false)
            ]).peak_at_range_edge[1])
        @test "a" in names(df_with_failure)    # param column still built from row 1
        @test ismissing(df_with_failure.a[2])  # failure row contributes no param value
    end

    @testset "a failure row holds its mechanism's text and id" begin
        m = _testhelper_uni_seed
        f = EnzymeRates.FitFailure(m, "boom")
        row = EnzymeRates._failure_row(f)
        # The text of the mechanism, which reads back as it, and the text's id; a failure
        # row records no parent.
        @test row.mechanism == sprint(show, m)
        @test EnzymeRates._parse_mechanism(row.mechanism) == m
        @test row.mechanism_id == EnzymeRates._mechanism_id(row.mechanism)
        @test ismissing(row.parent_id)
        @test row.error == "boom"
        @test ismissing(row.peak_at_range_edge)
        @test ismissing(row.key) && ismissing(row.fit_inherited) &&
              ismissing(row.fit_source_id) && ismissing(row.L_phantom)
    end

    @testset "the key of a law is the SHA-256 of its canonical text" begin
        # identify fits one law per key (spec section 10). The key reads the law, not its
        # printed text: it is the 64-hex SHA-256 of the canonical text, the same on every
        # call, and another law gets another key.
        m = _testhelper_uni_seed
        id = EnzymeRates._mechanism_id(sprint(show, m))
        c = EnzymeRates._canonical(EnzymeRates._law(m), id)
        again = EnzymeRates._canonical(EnzymeRates._law(m), id)
        @test (again.key, again.text, again.order) == (c.key, c.text, c.order)
        @test c.key == bytes2hex(EnzymeRates.sha256(c.text))
        @test occursin(r"^[0-9a-f]{64}$", c.key)
        ss_release = EnzymeRates.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) <--> E + P
            end
        end)
        other = EnzymeRates._mechanism_id(sprint(show, ss_release))
        @test EnzymeRates._canonical(EnzymeRates._law(ss_release), other).key != c.key
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
    results = redirect_stdout(devnull) do
        identify_rate_equation(prob;
            min_beam_width=1,
            loss_rel_threshold=1.0,
            loss_abs_threshold=0.0,
            max_param_count=8,
            n_cv_candidates=1,
            save_dir=save_dir,
            optimizer=cmaes_opt,
            n_restarts=1, maxtime=1.0)
    end

    # The selected mechanism is the 1-SE-rule row of cv_results.
    best_row = results.cv_results[
        EnzymeRates._select_best_row(results.cv_results), :]
    # The rate-law key of a mechanism's text: the key PASS1 gives it.
    _testhelper_key_of(text) =
        EnzymeRates._canonical(EnzymeRates._law(EnzymeRates._parse_mechanism(text)),
                               EnzymeRates._mechanism_id(text)).key

    @testset "mechanism recovery" begin
        # The best mechanism should fit the noiseless data with near-zero loss.
        # (Whether the search recovers the *most parsimonious* mechanism is a
        # search-quality property that needs heavy, seeded fits to test
        # reliably — too slow and too stochastic for CI, so it is not asserted
        # here.)
        @test best_row.loss < 0.01
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
        @test "key" in names(results.cv_results)
        @test "fit_source_id" in names(results.cv_results)
        # LOOCV candidate dedup invariant: within each n_params bucket, each key appears
        # at most once (the cv pool `_offer_cv!` keeps one slot per key). This catches a
        # regression where duplicates would enter LOOCV and waste compute / bias the
        # per-bucket "best". Each row's key is its mechanism's.
        for gdf in groupby(results.cv_results, :n_params)
            @test allunique(gdf.key)
        end
        @test all(r -> r.key == _testhelper_key_of(r.mechanism),
                  eachrow(results.cv_results))

        # cv_score_se = standard error of the row's per-fold scores.
        @test "cv_score_se" in names(results.cv_results)
        fold_cols = [Symbol("cv_fold_$g") for g in unique(prob.data.group)]
        for row in eachrow(results.cv_results)
            folds = [row[c] for c in fold_cols]
            @test row.cv_score ≈ mean(folds)
            @test row.cv_score_se ≈ std(folds) / sqrt(length(folds))
        end

        # The best row's text reads back as the selected mechanism.
        @test EnzymeRates.compile_mechanism(EnzymeRates._parse_mechanism(
            best_row.mechanism)) === results.best

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
        for col in (:n_params, :loss, :cv_score, :cv_score_se)
            @test col in propertynames(roundtrip)
        end
        @test nrow(roundtrip) == nrow(results.cv_results)
    end

    @testset "best mechanism computes rates" begin
        req = rate_equation_string(results.best)
        @test length(req) > 0
    end

    @testset "CSV output: every file's columns, texts and ids" begin
        files = sort(filter(f -> endswith(f, ".csv"), readdir(save_dir)))
        @test "initial_mechanisms.csv" in files
        @test isfile(joinpath(save_dir, "progress.log"))
        @test filesize(joinpath(save_dir, "progress.log")) > 0
        log_text = read(joinpath(save_dir, "progress.log"), String)
        @test startswith(log_text, "EnzymeRates v$(pkgversion(EnzymeRates))\n")
        @test occursin("new fits", log_text)
        @test occursin("skipped (>", log_text)
        @test occursin("best loss by n_params:", log_text)
        # The base tier and every iteration log the same four-line block: a header, the
        # pre-fit summary, the post-fit summary and the best-loss line.
        _testhelper_block(header) = Regex(
            "^" * header * "\\n  \\d+ new fits \\+ [^\\n]*\\n" *
            "  \\d+ errored \\| Success [^\\n]*\\n" *
            "  best loss by n_params: ", "m")
        @test occursin(_testhelper_block("Fitting \\d+ initial mechanisms…"), log_text)
        @test occursin(
            _testhelper_block("Iteration 1: \\d+ parents → \\d+ children"), log_text)
        @test !any(startswith(f, "params_estimate_") for f in files)
        iters = filter(f -> startswith(f, "equation_search_iteration_"), files)
        @test !isempty(iters)
        nums = sort(parse.(Int, replace.(iters,
            "equation_search_iteration_" => "", ".csv" => "")))
        @test nums == collect(1:length(nums))      # sequential, no gaps
        # An id is 16 hex digits, which CSV.read would take for a number when they look
        # like one (decimal digits, or digits around one `e` or `f`), so the ids are read
        # as text.
        _testhelper_read_csv(f) = CSV.read(joinpath(save_dir, f), DataFrame;
            types = Dict(:mechanism_id => String, :parent_id => String,
                         :fit_source_id => String))
        # Iteration 1 holds an L phantom, so the L_phantom branch below runs.
        @test any(f -> any(isequal(true), _testhelper_read_csv(f).L_phantom), files)
        init_df = _testhelper_read_csv("initial_mechanisms.csv")
        @test nrow(init_df) == length(unique!(
            collect(EnzymeRates.init_mechanisms(prob.reaction))))
        @test "key" in names(init_df)
        @test !("fit_inherited_from_estimate" in names(init_df))
        texts = Dict{String, Set{String}}()
        for f in files
            df_file = _testhelper_read_csv(f)
            @test "key" in names(df_file) && "fit_source_id" in names(df_file) &&
                  "L_phantom" in names(df_file)
            @test all(k -> occursin(r"^[0-9a-f]{64}$", k), skipmissing(df_file.key))
            # Each fitted row's key is its mechanism's, its fit's source is a 16-hex id,
            # and the row inherits its fit exactly when the source is another mechanism.
            for r in eachrow(df_file)
                # A failure row has no flag; a flagged row is an L phantom that reads
                # another mechanism's fit with L = 0.
                if ismissing(r.key)
                    @test ismissing(r.L_phantom)
                    continue
                end
                if r.L_phantom
                    m = EnzymeRates._parse_mechanism(r.mechanism)
                    @test EnzymeRates._is_phantom(EnzymeRates._law(m))
                    @test r.L == 0.0 && r.fit_inherited
                end
                @test r.key == _testhelper_key_of(r.mechanism)
                @test occursin(r"^[0-9a-f]{16}$", r.fit_source_id)
                @test r.fit_inherited == (r.fit_source_id != r.mechanism_id)
            end
            @test all(<=(8), skipmissing(df_file.n_params))     # max_param_count=8
            # Diagnostic columns: every saved CSV carries them.
            @test "parent_n_params" in names(df_file)
            @test "parent_id" in names(df_file)
            @test "peak_at_range_edge" in names(df_file)
            # Each row holds its mechanism's text, which reads back as the mechanism that
            # prints it, and the text's id.
            @test all(EnzymeRates._mechanism_id.(df_file.mechanism) .==
                      df_file.mechanism_id)
            @test all(t -> sprint(show, EnzymeRates._parse_mechanism(t)) == t,
                      df_file.mechanism)
            for (id, t) in zip(df_file.mechanism_id, df_file.mechanism)
                push!(get!(texts, id, Set{String}()), t)
            end
        end
        # Across the run's files, each id names one text.
        @test all(==(1) ∘ length, values(texts))
        # Base tier is parentless.
        @test all(ismissing, init_df.parent_n_params)
        @test all(ismissing, init_df.parent_id)
        # Iteration 1's parents are exactly the base-tier mechanisms; each child
        # records the id of the parent it expanded from, with the parent's own param
        # count.
        iter1 = _testhelper_read_csv("equation_search_iteration_1.csv")
        parent_np = Dict(zip(init_df.mechanism_id, init_df.n_params))
        @test any(!ismissing, iter1.parent_id)   # children have parents
        for r in eachrow(iter1)
            ismissing(r.parent_id) && continue
            @test r.parent_id in keys(parent_np)
            @test r.parent_n_params == parent_np[r.parent_id]
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
        # The CSV's law is the selected mechanism's printed law in named coefficients:
        # PASS1 prints it from the concrete mechanism, rate_equation_string from the
        # compiled one.
        @test bestf.rate_equation[1] ==
              rate_equation_string(results.best; coefficients = true)
        # The key reads back as text, and it is the selected mechanism's key.
        @test bestf.key[1] isa AbstractString
        @test bestf.key[1] == _testhelper_key_of(best_row.mechanism)
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

    @testset "_cv_fold_loss: finite per-fold scores, loud on fit failure" begin
        m = EnzymeRates.EnzymeMechanism(_testhelper_uni_seed)

        # 3 groups × 2 rows each so per-fold fits aren't degenerate
        data = DataFrame(
            S    = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            P    = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            Rate = [0.5, 0.8, 1.0, 1.1, 1.2, 1.3],
            group = [1, 1, 2, 2, 3, 3],
        )
        prob = IdentifyRateEquationProblem(_testhelper_uni_rxn, data; Keq=10.0)

        groups = unique(prob.data.group)
        scores = [EnzymeRates._cv_fold_loss(m, prob, g;
            optimizer=CMAEvolutionStrategyOpt(),
            n_restarts=2, maxtime=2.0, maxiters=500) for g in groups]

        @test scores isa Vector{Float64}
        # Require the success path: fitting MUST converge on
        # this trivial uni-uni fixture in 2s × 2 restarts. A
        # length-0 (full failure) result would let the per-fold
        # non-negativity + isfinite assertions below pass vacuously.
        @test length(scores) == 3
        @test all(s -> s >= 0.0, scores)
        @test all(isfinite, scores)

        # An unrecognized kwarg (`beam_fraction`) makes the fold's
        # `fit_rate_equation` call throw; `_cv_fold_loss` must propagate that
        # error, not swallow it (a corrupted CV must abort model selection).
        @test_throws Exception EnzymeRates._cv_fold_loss(
            m, prob, first(groups);
            optimizer=CMAEvolutionStrategyOpt(),
            n_restarts=1, maxtime=1.0, beam_fraction=0.5)
    end

end

# A fold's problem holds the fold's rows in the full fit's problem type, so the fold
# fits reuse the code compiled for the full fit.
@testset "_fold_problem: the fold's rows in the full fit's problem type" begin
    em = EnzymeRates.compile_mechanism(_testhelper_uni_seed)
    for data_form in (NamedTuple, DataFrame)
        prob = _testhelper_uni_prob(data_form)
        full = FittingProblem(em, prob.data; Keq=prob.Keq, absolute=prob.absolute,
                              peak_range=prob.peak_range)
        fold = EnzymeRates._fold_problem(em, prob, prob.data.group .== prob.data.group[3])
        @test typeof(fold) === typeof(full)
        @test fold.data.S == [3.0, 4.0]
        @test fold.data.P == [0.3, 0.4]
        @test fold.data.Rate == [1.0, 1.1]
        @test fold.Keq === 10.0 && fold.absolute === false &&
              fold.peak_range === (1e-12, 1e3)
        # A fold copies the problem's mode and range.
        prob_abs = IdentifyRateEquationProblem(prob.reaction, prob.data; Keq=prob.Keq,
                                               absolute=true, peak_range=(1e-9, 1e2))
        fold_abs = EnzymeRates._fold_problem(em, prob_abs,
                                             prob.data.group .== prob.data.group[3])
        @test fold_abs.absolute === true && fold_abs.peak_range === (1e-9, 1e2)
    end
end

# The number of method instances, inferred or compiled, of the methods EnzymeRates
# defines, keyword-call methods included.
function _testhelper_enzymerates_method_instances()
    ms = [m for m in methods(Core.kwcall) if m.module === EnzymeRates]
    for n in names(EnzymeRates; all = true)
        f = isdefined(EnzymeRates, n) ? getfield(EnzymeRates, n) : nothing
        f isa Function && append!(ms, m for m in methods(f) if m.module === EnzymeRates)
    end
    sum(m -> count(_ -> true, Base.specializations(m)), ms)
end

# LOOCV fits a candidate's folds after the candidate's full fit, and a fold's problem has
# the full fit's type. `_cv_fold_loss` and `_fold_problem` take the mechanism
# `@nospecialize`, so the fold then infers and compiles no EnzymeRates method for it.
@testset "_cv_fold_loss: a fold compiles nothing for its mechanism" begin
    rxn = @enzyme_reaction begin
        substrates: Sfold[C]
        products: Pfold[C]
    end
    # Metabolite names no other test uses, so no earlier test has fit these mechanisms.
    re_binding = @enzyme_mechanism begin
        substrates: Sfold
        products: Pfold
        steps: begin
            E + Sfold ⇌ E(Sfold)
            E(Sfold) <--> E(Pfold)
            E(Pfold) ⇌ E + Pfold
        end
    end
    ss_binding = @enzyme_mechanism begin
        substrates: Sfold
        products: Pfold
        steps: begin
            E + Sfold <--> E(Sfold)
            E(Sfold) <--> E(Pfold)
            E(Pfold) ⇌ E + Pfold
        end
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            Sfold = [1.0, 2.0, 3.0, 4.0], Pfold = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    optimizer = CMAEvolutionStrategyOpt()
    added = Int[]
    # Typed as the abstract mechanism, as in the pipeline, so each call dispatches on the
    # mechanism's run-time type.
    for em in EnzymeRates.AbstractEnzymeMechanism[re_binding, ss_binding]
        full = FittingProblem(em, prob.data; Keq=prob.Keq, absolute=prob.absolute,
                              peak_range=prob.peak_range)
        fit_rate_equation(full, optimizer; n_restarts=1, maxiters=10)
        n = _testhelper_enzymerates_method_instances()
        EnzymeRates._cv_fold_loss(em, prob, "G1"; optimizer, n_restarts=1, maxiters=10)
        push!(added, _testhelper_enzymerates_method_instances() - n)
    end
    # The first mechanism's fold compiles the fold path for every mechanism.
    @test added[2] == 0
end

# The CSV mechanism_id column: the first 16 hex digits of the SHA-256 of a mechanism's
# text, so every process and Julia version gives one id. The ids pinned here come from
# another SHA-256 implementation (Python's hashlib), applied to the texts that
# test_types.jl pins for these three mechanisms.
@testset "_mechanism_id cuts the SHA-256 of a mechanism's text to 16 hex digits" begin
    plain = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end
    ping_pong = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E + Q ⇌ E(Q)
        end
    end
    allosteric = @allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: I, J
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)     :: EqualAI
            E(S) <--> E(P)   :: OnlyA
            E(P) ⇌ E + P     :: EqualAI
        end
        regulatory_site(3): I::OnlyI, J::EqualAI
    end
    @test [EnzymeRates._mechanism_id(sprint(show, em))
           for em in (plain, ping_pong, allosteric)] ==
          ["f4ce3f6878d7229f", "dbcebb00ac48b73f", "3d6c468d7b38224b"]
end

@testset "csv writer" begin
    rows = [(
        n_params = 5, loss = 1.0, mechanism = "M", mechanism_id = "00000000000000aa",
        rate_equation = "v = 1", retcode = "Success", error = missing,
        params = (K_a = 2.0,), key = "abc", fit_source_id = "00000000000000aa",
        L_phantom = false,
        fit_inherited = false,
        peak_at_range_edge = true,
    )]
    mktempdir() do tmp
        # `texts` holds the run's ids, each with its mechanism's text.
        texts = Dict{String, String}()
        EnzymeRates._write_rows_csv!(tmp, "equation_search_iteration_3.csv", rows, texts)
        @test isfile(joinpath(tmp, "equation_search_iteration_3.csv"))
        df = CSV.read(joinpath(tmp, "equation_search_iteration_3.csv"), DataFrame)
        @test df.n_params == [5]
        @test df.K_a == [2.0]
        @test df.peak_at_range_edge == [true]
        @test "key" in names(df) && df.key == ["abc"]
        @test texts == Dict("00000000000000aa" => "M")
        # dir-creation branch: save_dir does not exist yet. A row may repeat an id with
        # its text (an expansion failure's row repeats its parent's).
        subdir = joinpath(tmp, "made")
        EnzymeRates._write_rows_csv!(subdir, "initial_mechanisms.csv", rows, texts)
        @test isfile(joinpath(subdir, "initial_mechanisms.csv"))
        # An id that names a second text is an error, raised before the file is written.
        clash = [merge(rows[1], (mechanism = "N",))]
        @test_throws(ErrorException("two mechanisms share the mechanism_id " *
                                    "00000000000000aa, the first 16 hex digits of the " *
                                    "SHA-256 of their texts, so the run's CSVs cannot " *
                                    "tell them apart:\nM\nN"),
                     EnzymeRates._write_rows_csv!(tmp, "equation_search_iteration_4.csv",
                                                  clash, texts))
        @test !isfile(joinpath(tmp, "equation_search_iteration_4.csv"))
    end
end

@testset "_select_count!: thresholds, floor, best loss, parsimony" begin
    # One call at count 5 with a fresh floor budget: `best` is the count's best loss
    # over the whole search, and `others` holds the best losses of other counts.
    # Each mechanism its own key, so the key ranks are the mechanism ranks.
    _testhelper_selected(losses, best; rel, add = 0.0, width = 1,
                         others = Dict{Int,Float64}(), parsimony = 1.0) =
        EnzymeRates._select_count!(Dict{Int,Int}(), merge(Dict(5 => best), others), 5,
            losses, string.(eachindex(losses)); loss_rel_threshold = rel,
            loss_abs_threshold = add, loss_parsimony_threshold = parsimony,
            min_beam_width = width, selected = Set{String}())

    losses = [1.0, 1.5, 2.5, 5.0, 10.0]
    @test _testhelper_selected(losses, 1.0; rel = 2.0) == [1, 2]
    @test _testhelper_selected(losses, 1.0; rel = 2.0, width = 4) == [1, 2, 3, 4]

    # The additive term keeps a near-zero best loss from collapsing the cutoff.
    @test _testhelper_selected([1e-6, 0.005, 0.05], 1e-6; rel = 2.0, add = 0.01) == [1, 2]

    # Indices come back in INPUT order, not loss order.
    @test _testhelper_selected([5.0, 1.0, 10.0, 2.0], 1.0; rel = 2.5) == [2, 4]

    # The relative cutoff uses the count's best loss, which can differ from this
    # sweep's minimum.
    losses = [1.0, 1.5, 3.0]
    @test _testhelper_selected(losses, 1.0; rel = 1.2) == [1]         # cutoff 1.2
    @test _testhelper_selected(losses, 2.0; rel = 1.2) == [1, 2]      # cutoff 2.4
    # floor still honored
    @test _testhelper_selected(losses, 0.0; rel = 1.0, width = 2) == [1, 2]

    # Floor guarantee: a parsimony cutoff below every loss admits nothing via
    # the loss filter, yet min_beam_width still keeps the top-k by loss.
    losses = [1.0, 1.5, 2.5, 5.0, 10.0]
    @test _testhelper_selected(losses, 1.0; rel = 2.0, width = 2,
                               others = Dict(4 => 0.5)) == [1, 2]

    # Tightening: a parsimony cutoff stricter than the rel/abs cutoff lowers the
    # combined cutoff to 2.0, so indices 1 and 2 (losses 1.0, 1.5) pass and
    # index 3 (2.5) is dropped. Without it, rel=10 would admit all four.
    losses = [1.0, 1.5, 2.5, 5.0]
    @test _testhelper_selected(losses, 1.0; rel = 10.0, others = Dict(4 => 2.0)) == [1, 2]

    # No-op: with no smaller count fit yet the parsimony term is dropped, whatever its
    # threshold, and a larger count is no parsimony reference.
    @test _testhelper_selected(losses, 1.0; rel = 2.0, parsimony = 0.0) ==
          _testhelper_selected(losses, 1.0; rel = 2.0, parsimony = Inf) == [1, 2]
    @test _testhelper_selected(losses, 1.0; rel = 2.0, others = Dict(6 => 0.1)) == [1, 2]

    # Interaction: min() picks the smaller cutoff. With best loss 2.0 the
    # rel cutoff is 2.4 (admits 1,2); a tighter parsimony cutoff of 1.0 lowers
    # it to just the single best.
    losses = [1.0, 1.5, 3.0]
    @test _testhelper_selected(losses, 2.0; rel = 1.2) == [1, 2]
    @test _testhelper_selected(losses, 2.0; rel = 1.2, others = Dict(4 => 1.0)) == [1]
end

@testset "all base fits fail: failure CSV written, then raises" begin
    # A `solver_kwargs` option no optimizer recognizes is forwarded
    # verbatim to `Optimization.solve` and makes every fit throw (a bogus
    # name keeps this independent of whether any real option is supported).
    # Per-mechanism fit failures are isolated in `_process_batch`, so the
    # base tier is then empty and the pipeline raises. The contract under
    # test is that the all-base-fail path persists the failure rows to
    # `initial_mechanisms.csv` before raising (for cluster debugging).
    prob = _testhelper_uni_prob(NamedTuple)
    tmp = mktempdir()
    @test_throws ErrorException redirect_stdout(devnull) do
        identify_rate_equation(
            prob; solver_kwargs=(; not_a_real_solver_option=1),
            optimizer=CMAEvolutionStrategyOpt(),
            n_restarts=1, maxtime=1.0, save_dir=tmp)
    end
    # Failure rows were written before the re-raise: a CSV exists whose rows
    # are all failures (non-missing `error`, missing `key`).
    @test isfile(joinpath(tmp, "initial_mechanisms.csv"))
    fail_df = CSV.read(joinpath(tmp, "initial_mechanisms.csv"), DataFrame)
    @test nrow(fail_df) >= 1
    @test all(.!ismissing.(fail_df.error))
    @test all(ismissing.(fail_df.key))
    # The post-fit summary is logged before the raise.
    @test occursin(
        "\n  $(nrow(fail_df)) errored | Success 0.0% | non-Success retcode 0.0%\n",
        read(joinpath(tmp, "progress.log"), String))
end

@testset "_select_best_row: 1-SE rule on the best equation's fold scores" begin
    # Best row = lowest cv_score (n=7, 0.10, SE 0.02) → cutoff 0.12. n=5 (0.115)
    # passes, n=3 (0.13) does not → smallest passing n_params is 5.
    cv_df = DataFrame(
        n_params    = [3, 5, 7],
        cv_score    = [0.13, 0.115, 0.10],
        cv_score_se = [0.001, 0.001, 0.02],
    )
    @test EnzymeRates._select_best_row(cv_df) == 2

    # The cutoff uses the best row's SE only: a simpler row's own large SE
    # does not rescue it.
    cv_df_own_se = DataFrame(
        n_params    = [3, 7],
        cv_score    = [0.5, 0.10],
        cv_score_se = [10.0, 0.02],
    )
    @test EnzymeRates._select_best_row(cv_df_own_se) == 2

    # A score exactly at the cutoff passes (≤).
    cv_df_edge = DataFrame(
        n_params    = [3, 7],
        cv_score    = [0.75, 0.5],
        cv_score_se = [0.0, 0.25],
    )
    @test EnzymeRates._select_best_row(cv_df_edge) == 1

    # se_threshold scales the cutoff: 0.10 + 0.5 * 0.02 = 0.11 < 0.115.
    @test EnzymeRates._select_best_row(cv_df; se_threshold = 0.5) == 3
    @test EnzymeRates._select_best_row(cv_df; se_threshold = 2.0) == 1

    # Within the selected n_params the lowest-cv_score row wins, not the first
    # passing row and not the lowest training loss.
    cv_df_multi = DataFrame(
        n_params    = [5, 5, 5, 7],
        cv_score    = [0.118, 0.112, 0.30, 0.10],
        cv_score_se = [0.001, 0.001, 0.001, 0.02],
        loss        = [0.01, 0.05, 0.02, 0.03],
    )
    @test EnzymeRates._select_best_row(cv_df_multi) == 2

    # Rows with more parameters than the best row are never selected, even
    # inside the cutoff.
    cv_df_larger = DataFrame(
        n_params    = [5, 9],
        cv_score    = [0.10, 0.101],
        cv_score_se = [0.02, 0.02],
    )
    @test EnzymeRates._select_best_row(cv_df_larger) == 1

    # Ties on cv_score resolve to the smallest n_params (parsimony), and the
    # tied simpler row is the one returned.
    cv_df_tie = DataFrame(
        n_params    = [7, 3, 5],
        cv_score    = [0.115, 0.115, 0.115],
        cv_score_se = [0.0, 0.0, 0.0],
    )
    @test EnzymeRates._select_best_row(cv_df_tie) == 2

    # Single row.
    cv_df_single = DataFrame(
        n_params = [4], cv_score = [0.15], cv_score_se = [0.01])
    @test EnzymeRates._select_best_row(cv_df_single) == 1
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

@testset "the key partitions the bi-bi init laws into 239 classes" begin
    # The 239 bi-bi init mechanisms (55 seeds and their 184 merged and Theorell–Chance
    # variants) have 239 laws no renaming of fitted constants makes equal, so the key
    # gives 239 classes. On a fixed sample the same law gives the same key twice.
    reaction = @enzyme_reaction(begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end)
    buckets = Dict{String, Vector{Int}}()
    for (i, m) in enumerate(EnzymeRates.init_mechanisms(reaction))
        id = EnzymeRates._mechanism_id(sprint(show, m))
        c = EnzymeRates._canonical(EnzymeRates._law(m), id)
        @test !c.capped
        push!(get!(buckets, c.key, Int[]), i)
        i % 25 == 1 && @test EnzymeRates._canonical(EnzymeRates._law(m), id).key == c.key
    end
    @test length(buckets) == 239
end

@testset "_select_count! cumulative per-count floor" begin
    expanded = Dict{Int,Int}()
    best = Dict(5 => 1.0)
    # Each mechanism its own key; the selected keys persist across sweeps, as in the beam.
    kw = (loss_rel_threshold=1.0, loss_abs_threshold=0.0,
          loss_parsimony_threshold=1.0, min_beam_width=3,
          selected=Set{String}())
    # Sweep 1 at count 5: rel cutoff admits only the best (loss 1.0); the
    # floor budget (3) tops it up to the top 3 by loss. expanded[5] -> 3.
    sel1 = EnzymeRates._select_count!(expanded, best, 5, [1.0, 2.0, 3.0, 4.0, 5.0],
                                      ["a", "b", "c", "d", "e"]; kw...)
    @test sort(sel1) == [1, 2, 3]
    @test expanded[5] == 3

    # Sweep 2 at count 5: budget spent (3 of 3). New mechanisms all above the
    # cutoff -> the floor admits NONE (unlike the old per-sweep floor, which
    # would grant a fresh 3). expanded[5] stays 3.
    sel2 = EnzymeRates._select_count!(expanded, best, 5, [10.0, 11.0, 12.0],
                                      ["f", "g", "h"]; kw...)
    @test isempty(sel2)
    @test expanded[5] == 3

    # A cutoff-passer is still admitted after the floor is spent.
    sel3 = EnzymeRates._select_count!(expanded, best, 5, [1.0, 20.0], ["i", "j"]; kw...)
    @test sel3 == [1]
    @test expanded[5] == 4
end

@testset "§1 parsimony cutoff = threshold * min over all counts < c" begin
    # No floor and a loose relative cutoff (10 × the count's best), so the parsimony
    # cutoff alone decides: 1.01 × the best loss over the counts below c.
    _testhelper_selected(best_loss_by_count, c, losses) = EnzymeRates._select_count!(
        Dict{Int,Int}(), best_loss_by_count, c, losses, string.(eachindex(losses));
        loss_rel_threshold=10.0, loss_abs_threshold=0.0, loss_parsimony_threshold=1.01,
        min_beam_width=0, selected=Set{String}())
    # No count < c: no parsimony term (else 0.15 > 1.01*0.02 would be dropped).
    @test _testhelper_selected(Dict(5=>0.02), 5, [0.02, 0.15]) == [1, 2]
    # min over <c, not c-1: the cutoff is 1.01*0.02, not 1.01*0.03.
    @test _testhelper_selected(Dict(5=>0.02, 6=>0.05, 7=>0.03, 8=>0.02), 8,
                               [0.0201, 0.0203, 0.03]) == [1]
    # count gap: c-1=6 absent, the cutoff is 1.01*0.02.
    @test _testhelper_selected(Dict(5=>0.02, 7=>0.02), 7, [0.0201, 0.0203]) == [1]
    # non-monotone → true min: the cutoff is 1.01*0.01, not 1.01*0.04.
    @test _testhelper_selected(Dict(5=>0.01, 6=>0.04, 7=>0.01), 7, [0.0100, 0.0102]) == [1]
end

@testset "_select_count! charges the budget once per key" begin
    # Spec section 10: the beam counts each key once. Five mechanisms on two keys, the
    # three of key a sharing their fit's loss. With a budget of two, ranking mechanisms
    # would spend it on two of a's three and leave b out; ranking keys spends one unit on
    # each key and expands all five mechanisms.
    expanded = Dict{Int,Int}()
    selected = Set{String}()
    sel = EnzymeRates._select_count!(expanded, Dict(5 => 1.0), 5,
        [1.0, 1.0, 1.0, 2.0, 2.0], ["a", "a", "a", "b", "b"];
        loss_rel_threshold = 1.0, loss_abs_threshold = 0.0, loss_parsimony_threshold = 1.0,
        min_beam_width = 2, selected)
    @test sel == [1, 2, 3, 4, 5]
    @test expanded[5] == 2
    @test selected == Set(["a", "b"])
end

@testset "a late member of a selected key expands without budget" begin
    # Key a is selected in the first sweep at loss 2.0, then the count's best. A later
    # sweep brings a structurally new member of a after another key lowered the best to
    # 1.0: the cutoff (1.0) and the spent budget would refuse a mechanism at 2.0, but a's
    # verdict stands, so the member expands and the budget stays spent on a alone.
    expanded = Dict{Int,Int}()
    selected = Set{String}()
    kw = (loss_rel_threshold = 1.0, loss_abs_threshold = 0.0,
          loss_parsimony_threshold = 1.0, min_beam_width = 1, selected)
    @test EnzymeRates._select_count!(expanded, Dict(5 => 2.0), 5, [2.0], ["a"];
                                     kw...) == [1]
    @test expanded[5] == 1
    @test EnzymeRates._select_count!(expanded, Dict(5 => 1.0), 5, [2.0, 3.0],
                                     ["a", "c"]; kw...) == [1]
    @test expanded[5] == 1
    @test selected == Set(["a"])
end

@testset "a late member of a rejected key stays out" begin
    # Key b loses on rank in the first sweep and is rejected. A later member of b is
    # judged again and refused, while a later member of the selected key a expands.
    expanded = Dict{Int,Int}()
    selected = Set{String}()
    kw = (loss_rel_threshold = 1.0, loss_abs_threshold = 0.0,
          loss_parsimony_threshold = 1.0, min_beam_width = 1, selected)
    @test EnzymeRates._select_count!(expanded, Dict(5 => 1.0), 5, [1.0, 1.5],
                                     ["a", "b"]; kw...) == [1]
    @test EnzymeRates._select_count!(expanded, Dict(5 => 1.0), 5, [1.5, 1.0],
                                     ["b", "a"]; kw...) == [2]
    @test expanded[5] == 1
    @test selected == Set(["a"])
end

@testset "a later sweep expands a selected key's member, judges a refused key again and \
          keeps the spent budget" begin
    # Two sweeps of the beam at target 3 that share one key set. Sweep 1 selects key a
    # (loss 1.0, under the cutoff 1.3) and refuses b (5.0), since a spends count 3's
    # budget of one. Sweep 2 pops two counts. At count 3 (cutoff 0.65, budget spent), a2
    # expands by a's verdict though its loss is over the cutoff, the new key c passes the
    # cutoff, and b2 is judged again and refused; at count 2, d is within the budget. The
    # parents are compared as a set, so the frontier's bucket order cannot matter.
    expanded, best, selected = Dict{Int,Int}(), Dict(3 => 1.0), Set{String}()
    kw = (loss_rel_threshold = 1.3, loss_abs_threshold = 0.0,
          loss_parsimony_threshold = Inf, min_beam_width = 1)
    a, b = _testhelper_entry("a", 1.0, 3), _testhelper_entry("b", 5.0, 3)
    frontier = Dict(3 => [a, b])
    @test EnzymeRates._sweep_parents!(frontier, 3, expanded, best, selected; kw...) == [a]
    @test isempty(frontier)
    @test expanded == Dict(3 => 1)
    @test selected == Set(["a"])

    a2, b2 = _testhelper_entry("a", 1.0, 3), _testhelper_entry("b", 5.0, 3)
    c, d = _testhelper_entry("c", 0.5, 3), _testhelper_entry("d", 0.2, 2)
    frontier = Dict(3 => [a2, b2, c], 2 => [d])
    best[3], best[2] = 0.5, 0.2
    parents = EnzymeRates._sweep_parents!(frontier, 3, expanded, best, selected; kw...)
    @test length(parents) == 3
    @test Set(objectid.(parents)) == Set(objectid.([a2, c, d]))
    @test isempty(frontier)
    @test expanded == Dict(3 => 2, 2 => 1)
    @test selected == Set(["a", "c", "d"])
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

    # _prefit_summary: the five pre-fit buckets, no errored, no percentages.
    pre = EnzymeRates._prefit_summary(2, 0, 4, 2, 3;
        max_param_count=8, eq_complexity_filter=337)
    @test occursin("2 new fits + 0 inherited + 3 skipped (already fit) + " *
                   "4 skipped (>8 params) + 2 skipped (>337 complexity)", pre)
    @test !occursin("errored", pre)
    @test !occursin("Success", pre)

    # _postfit_summary: errored + success/non-Success over the fitted set.
    mech = _testhelper_uni_seed
    row = (n_params=3, loss=0.5, mechanism="M", rate_equation="v",
           retcode="Success", error=missing, params=(K = 1.0,), key="abc",
           fit_inherited=false)
    e_succ = EnzymeRates.BatchEntry(mech, 3, 0.5, :Success, "a", "a", row)
    e_mt   = EnzymeRates.BatchEntry(mech, 3, 0.9, :MaxTime, "b", "b", row)
    f      = EnzymeRates.FitFailure(mech, "StackOverflowError: ")
    post = EnzymeRates._postfit_summary([e_succ, e_mt], [f])
    @test occursin("1 errored", post)
    @test occursin("Success 50.0%", post)              # 1 of 2 fitted
    @test occursin("non-Success retcode 50.0%", post)  # e_mt is :MaxTime
    @test !occursin("best loss", post)
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
    ms = unique!(collect(EnzymeRates.init_mechanisms(_testhelper_uni_rxn)))
    prob = _testhelper_uni_prob(DataFrame)

    entries, failures = EnzymeRates._process_batch(ms, prob;
        optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=20, n_restarts=1, maxtime=1.0)
    @test !isempty(entries)
    @test isempty(failures)
    @test all(e -> e isa EnzymeRates.BatchEntry, entries)
    @test all(e -> e.retcode isa Symbol, entries)
    @test all(e -> e.n_params == length(e.row.params), entries)
    @test all(e -> occursin(r"^[0-9a-f]{64}$", e.key) && e.key == e.row.key == e.fit_key,
              entries)

    # cap filter: nothing over the cap is fit (and it is not a failure), and the
    # pre-fit summary logs every mechanism as a param-count skip.
    capped_log = String[]
    capped_entries, capped_failures = EnzymeRates._process_batch(ms, prob;
        optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=0, eq_complexity_filter=337, n_restarts=1, maxtime=1.0,
        log = msg -> push!(capped_log, msg))
    @test isempty(capped_entries)
    @test isempty(capped_failures)
    @test capped_log == ["0 new fits + 0 inherited + 0 skipped (already fit) + " *
                         "$(length(ms)) skipped (>0 params) + 0 skipped (>337 complexity)"]

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

    # seen set: a structure already produced is skipped, not reprocessed.
    m = first(ms)
    seen = Set{UInt64}()
    e1, f1, ps1, cs1, ss1 = EnzymeRates._process_batch([m], prob;
        optimizer=CMAEvolutionStrategyOpt(), max_param_count=20,
        n_restarts=1, maxtime=1.0, memo=Dict{String,EnzymeRates._Fit}(), seen)
    @test ss1 == 0 && length(e1) == 1

    # Same structure again in a later batch → seen-skipped, no new entry.
    seen_log = String[]
    e2, f2, ps2, cs2, ss2 = EnzymeRates._process_batch([m], prob;
        optimizer=CMAEvolutionStrategyOpt(), max_param_count=20, eq_complexity_filter=337,
        n_restarts=1, maxtime=1.0, memo=Dict{String,EnzymeRates._Fit}(), seen,
        log = msg -> push!(seen_log, msg))
    @test ss2 == 1 && isempty(e2) && isempty(f2)
    @test seen_log == ["0 new fits + 0 inherited + 1 skipped (already fit) + " *
                       "0 skipped (>20 params) + 0 skipped (>337 complexity)"]

    # A structure repeated within one batch is fit once; the repeat is seen-skipped.
    e3, f3, ps3, cs3, ss3 = EnzymeRates._process_batch([m, m], prob;
        optimizer=CMAEvolutionStrategyOpt(), max_param_count=20,
        n_restarts=1, maxtime=1.0, memo=Dict{String,EnzymeRates._Fit}(), seen=Set{UInt64}())
    @test ss3 == 1 && length(e3) == 1 && isempty(f3)
end

# A random-order ter-ter (all binding/release orders, all SS): V×τ ≈ 5.9M, far
# above any complexity threshold. Shared by the filter and derivation-guard tests.
_testhelper_random_terter() = @enzyme_mechanism begin
    substrates: S1, S2, S3
    products:   P1, P2, P3
    steps: begin
        E + S1 <--> E(S1); E + S2 <--> E(S2); E + S3 <--> E(S3)
        E(S1) + S2 <--> E(S1, S2); E(S1) + S3 <--> E(S1, S3); E(S2) + S1 <--> E(S1, S2)
        E(S2) + S3 <--> E(S2, S3); E(S3) + S1 <--> E(S1, S3); E(S3) + S2 <--> E(S2, S3)
        E(S1, S2) + S3 <--> E(S1, S2, S3); E(S1, S3) + S2 <--> E(S1, S2, S3)
        E(S2, S3) + S1 <--> E(S1, S2, S3); E(S1, S2, S3) <--> E(P1, P2, P3)
        E(P1, P2, P3) <--> E(P1, P2) + P3; E(P1, P2, P3) <--> E(P1, P3) + P2
        E(P1, P2, P3) <--> E(P2, P3) + P1; E(P1, P2) <--> E(P1) + P2
        E(P1, P2) <--> E(P2) + P1; E(P1, P3) <--> E(P1) + P3
        E(P1, P3) <--> E(P3) + P1; E(P2, P3) <--> E(P2) + P3
        E(P2, P3) <--> E(P3) + P2; E(P1) <--> E + P1
        E(P2) <--> E + P2; E(P3) <--> E + P3
    end
end

@testset "eq_complexity_filter: random bi-bi passes, ter-ter skipped (default 337)" begin
    rxn = @enzyme_reaction begin
        substrates: S1[C2H4], S2[C2H2]
        products:   P1[C2H2], P2[C2H4]
    end
    data = DataFrame(
        S1 = [1.0, 2.0, 1.0, 2.0], S2 = [1.0, 1.0, 2.0, 2.0],
        P1 = [0.1, 0.2, 0.1, 0.2], P2 = [0.1, 0.1, 0.2, 0.2],
        Rate = [0.5, 0.8, 0.9, 1.1], group = [1, 1, 2, 2])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=2.0)
    random_bibi = @enzyme_mechanism begin        # V×τ = 336 — the ceiling; must pass
        substrates: S1, S2
        products:   P1, P2
        steps: begin
            E + S1 <--> E(S1)
            E + S2 <--> E(S2)
            E(S1) + S2 <--> E(S1, S2)
            E(S2) + S1 <--> E(S1, S2)
            E(S1, S2) <--> E(P1, P2)
            E(P1, P2) <--> E(P1) + P2
            E(P1, P2) <--> E(P2) + P1
            E(P1) <--> E + P1
            E(P2) <--> E + P2
        end
    end
    # At the default 337, a random-order bi-bi (V×τ = 336) passes, while a
    # random-order ter-ter (V×τ ≈ 5.9M) is complexity-skipped in PASS-1 — before
    # fitting, so its metabolite mismatch with the bi-bi problem is never reached.
    batch = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[
        EnzymeRates.Mechanism(random_bibi),
        EnzymeRates.Mechanism(_testhelper_random_terter())]
    entries, failures, n_param_skip, n_cx_skip = EnzymeRates._process_batch(
        batch, prob; optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=20, eq_complexity_filter=337, n_restarts=1, maxtime=1.0)
    @test n_cx_skip == 1                             # ter-ter skipped
    @test n_param_skip == 0
    @test length(entries) + length(failures) == 1   # random bi-bi (336) reached fitting
    @test all(e -> EnzymeRates._eq_complexity(e.mech) <= 337, entries)

    # A tighter cap catches the bi-bi too — the filter is tunable and both skips
    # are counted.
    _, _, _, cx_tight = EnzymeRates._process_batch(
        batch, prob; optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=20, eq_complexity_filter=100, n_restarts=1, maxtime=1.0)
    @test cx_tight == 2
end

@testset "derivation guard (filter > MAX_RATE_EQUATION_TERMS) → FitFailure, not crash" begin
    # A random-order ter-ter is far above MAX_RATE_EQUATION_TERMS (V×τ ≈ 5.9M).
    # With eq_complexity_filter raised above it, it is NOT complexity-skipped, so
    # it reaches derivation — where the MAX_RATE_EQUATION_TERMS guard aborts it.
    # That error must be caught and recorded as a FitFailure (→ CSV), never fatal.
    rxn = @enzyme_reaction begin
        substrates: S1[C2H4], S2[C2H2], S3[C2]
        products:   P1[C2H2], P2[C2H4], P3[C2]
    end
    data = DataFrame(
        S1 = [1.0, 2.0], S2 = [1.0, 2.0], S3 = [1.0, 2.0],
        P1 = [0.1, 0.2], P2 = [0.1, 0.2], P3 = [0.1, 0.2],
        Rate = [0.5, 0.8], group = [1, 2])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=2.0)
    m = EnzymeRates.Mechanism(_testhelper_random_terter())
    @test EnzymeRates._eq_complexity(m) > EnzymeRates.MAX_RATE_EQUATION_TERMS
    batch = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m]
    entries, failures, n_param_skip, n_cx_skip = EnzymeRates._process_batch(
        batch, prob; optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=100, eq_complexity_filter=typemax(Int),
        n_restarts=1, maxtime=1.0)
    @test n_cx_skip == 0                     # NOT complexity-skipped (filter above V×τ)
    @test isempty(entries)                   # never fit
    @test length(failures) == 1              # recorded, not a crash
    @test occursin("polynomial terms", failures[1].error)
end

@testset "eq_complexity_filter counts both conformations" begin
    # The 336-term random-order bi-bi at one catalytic subunit with its chemistry and its
    # product bindings at free E :OnlyA, a K-type the tagging move emits: the inactive
    # conformation keeps the substrate square E, E(S1), E(S1, S2), E(S2), four segments and
    # four spanning trees, so the law counts 336 + 16 = 352. The default filter skips it,
    # where the active count alone would let it through. A filter of 700 lets it through
    # to the parameter cap, which skips its 16 fitted constants without a fit.
    rxn = @enzyme_reaction begin
        substrates: S1[C2H4], S2[C2H2]
        products:   P1[C2H2], P2[C2H4]
    end
    data = DataFrame(
        S1 = [1.0, 2.0, 1.0, 2.0], S2 = [1.0, 1.0, 2.0, 2.0],
        P1 = [0.1, 0.2, 0.1, 0.2], P2 = [0.1, 0.1, 0.2, 0.2],
        Rate = [0.5, 0.8, 0.9, 1.1], group = [1, 1, 2, 2])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=2.0)
    ktype = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S1, S2
        products: P1, P2
        catalytic_multiplicity: 1
        steps: begin
            E + P1 <--> E(P1) :: OnlyA
            E + P2 <--> E(P2) :: OnlyA
            E + S1 <--> E(S1) :: EqualAI
            E + S2 <--> E(S2) :: EqualAI
            E(P1) + P2 <--> E(P1, P2) :: EqualAI
            E(P2) + P1 <--> E(P1, P2) :: EqualAI
            E(S1) + S2 <--> E(S1, S2) :: EqualAI
            E(S1, S2) <--> E(P1, P2) :: OnlyA
            E(S2) + S1 <--> E(S1, S2) :: EqualAI
        end
    end)
    cm = ktype.catalytic
    @test _testhelper_laplacian_tree_count(cm, EnzymeRates._active_mask(ktype)) == 336
    @test _testhelper_laplacian_tree_count(cm, EnzymeRates._inactive_mask(ktype)) == 16
    @test EnzymeRates._eq_complexity(ktype) == 352
    @test EnzymeRates._independent_param_count(ktype) == 16
    batch = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[ktype]
    _, _, n_param_skip, n_cx_skip = EnzymeRates._process_batch(
        batch, prob; optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=15, eq_complexity_filter=337, n_restarts=1, maxtime=1.0)
    @test (n_cx_skip, n_param_skip) == (1, 0)
    _, _, n_param_wide, n_cx_wide = EnzymeRates._process_batch(
        batch, prob; optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=15, eq_complexity_filter=700, n_restarts=1, maxtime=1.0)
    @test (n_cx_wide, n_param_wide) == (0, 1)
end

@testset "_ingest! and cv pool" begin
    frontier = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    cv_pool  = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    best     = Dict{Int,Float64}()
    # two distinct equations + one duplicate-eq with worse loss, n_cv=2
    improved = EnzymeRates._ingest!(frontier, cv_pool, best,
        [_testhelper_entry("a", 2.0, 5), _testhelper_entry("b", 1.0, 5),
         _testhelper_entry("a", 3.0, 5)];
        n_cv_candidates=2)
    @test improved == Set([5])                 # count 5 first appeared
    @test length(frontier[5]) == 3            # frontier keeps ALL
    @test best[5] == 1.0                       # running min
    @test length(cv_pool[5]) == 2              # bounded, distinct keys
    # the kept "a" entry is the lower-loss one (2.0, not 3.0)
    a = only(filter(e -> e.key == "a", cv_pool[5]))
    @test a.loss == 2.0
    # A later batch reports only the counts whose best strictly dropped or first
    # appeared: count 5 ties its best (1.0), count 6 first appears, and count 7's
    # entry is worse than its best (0.5).
    best[7] = 0.5
    @test EnzymeRates._ingest!(frontier, cv_pool, best,
        [_testhelper_entry("c", 1.0, 5), _testhelper_entry("d", 4.0, 6),
         _testhelper_entry("e", 0.6, 7)];
        n_cv_candidates=2) == Set([6])
    @test EnzymeRates._ingest!(frontier, cv_pool, best,
        [_testhelper_entry("f", 0.9, 5)]; n_cv_candidates=2) == Set([5])
    @test isempty(EnzymeRates._ingest!(frontier, cv_pool, best,
        EnzymeRates.BatchEntry[]; n_cv_candidates=2))

    # A fit whose loss is not finite (fit_rate_equation returns Inf when no restart is
    # finite) is skipped: it joins neither the frontier nor the cv pool and sets no best
    # loss, even at a count it is the first to reach. A non-finite fit beside a finite one
    # at count 5 leaves a free slot rather than enter LOOCV, and counts 6 and 7, which
    # only non-finite fits reach, get no LOOCV candidate at all.
    frontier_nf = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    cv_pool_nf  = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    best_nf     = Dict{Int,Float64}()
    @test EnzymeRates._ingest!(frontier_nf, cv_pool_nf, best_nf,
        [_testhelper_entry("g", Inf, 5), _testhelper_entry("h", 1.0, 5),
         _testhelper_entry("i", NaN, 6), _testhelper_entry("j", Inf, 7)];
        n_cv_candidates=2) == Set([5])
    @test Set(keys(frontier_nf)) == Set(keys(cv_pool_nf)) == Set([5])
    @test [e.key for e in frontier_nf[5]] == ["h"]
    @test [e.key for e in cv_pool_nf[5]] == ["h"]
    @test best_nf == Dict(5 => 1.0)

    # `_offer_cv!` keeps at most one entry per key: a repeat key updates its own slot to
    # the lower loss, never consuming a second.
    pool = EnzymeRates.BatchEntry[]
    for (loss, h) in [(1.0, "a"), (0.5, "a"), (2.0, "b")]
        EnzymeRates._offer_cv!(pool, _testhelper_entry(h, loss, 5), 5)
    end
    @test allunique([e.key for e in pool])
    @test length(pool) == 2
    @test only(filter(e -> e.key == "a", pool)).loss == 0.5
end

@testset "all-cap-skipped expansion batch is reported (M2)" begin
    # uni-uni base mechanism has 3 params; every child has 4. With
    # max_param_count=3 the base fits but the whole expansion batch is
    # cap-skipped — no rows, no CSV — so it must still log its header and
    # pre-fit summary.
    prob = _testhelper_uni_prob(NamedTuple)
    tmp = mktempdir()
    # An explicit non-default loss_parsimony_threshold proves the keyword is accepted
    # (an unknown keyword throws at the call boundary).
    redirect_stdout(devnull) do
        identify_rate_equation(prob;
            optimizer=CMAEvolutionStrategyOpt(),
            min_beam_width=1, loss_rel_threshold=1.0, loss_abs_threshold=0.0,
            loss_parsimony_threshold=2.0,
            max_param_count=3, n_cv_candidates=1, n_restarts=1, maxtime=1.0,
            save_dir=tmp)
    end
    log_text = read(joinpath(tmp, "progress.log"), String)
    @test occursin(Regex(
        "^Iteration 1: \\d+ parents → (\\d+) children\\n  0 new fits \\+ 0 inherited \\+ " *
        "0 skipped \\(already fit\\) \\+ \\1 skipped \\(>3 params\\) \\+ " *
        "0 skipped \\(>337 complexity\\)\\nCross-validating", "m"), log_text)
    # The all-skip batch produced no rows, so no iteration CSV was written.
    @test !any(startswith(f, "equation_search_iteration_") for f in readdir(tmp))
end

@testset "_required_regulators" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: A, B
        competitive_inhibitors: I, J
        oligomeric_state: 2
    end
    # Empty optional lists: every declared regulator is required, split by kind.
    ra, rc = EnzymeRates._required_regulators(rxn, Symbol[], Symbol[])
    @test ra == Set([:A, :B])
    @test rc == Set([:I, :J])
    # All regulators optional: both required sets empty (so the beam keeps its
    # unregulated init seed).
    ra_opt, rc_opt = EnzymeRates._required_regulators(rxn, [:A, :B], [:I, :J])
    @test isempty(ra_opt)
    @test isempty(rc_opt)
    # Mixed: one of each kind optional, the other required.
    ra_mix, rc_mix = EnzymeRates._required_regulators(rxn, [:A], [:J])
    @test ra_mix == Set([:B])
    @test rc_mix == Set([:I])

    # A dual-role name (ATP declared as BOTH allosteric regulator and
    # competitive inhibitor) is required/optional per role, independently:
    # marking ATP optional as a competitive inhibitor leaves it required as an
    # allosteric regulator.
    dual_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: ATP
        competitive_inhibitors: ATP
        oligomeric_state: 2
    end
    ra_dual, rc_dual = EnzymeRates._required_regulators(dual_rxn, Symbol[], [:ATP])
    @test :ATP in ra_dual        # still required as an allosteric regulator
    @test !(:ATP in rc_dual)     # optional as a competitive inhibitor
end

@testset "removed kwargs error at the identify boundary" begin
    # popsize/verbose are no longer named kwargs and there is no catch-all,
    # so they are rejected immediately at the call boundary (before any
    # fitting or CSV write) — distinct from a solver-rejected solver_kwargs
    # option, which fails inside fitting.
    prob = _testhelper_uni_prob(NamedTuple)
    @test_throws Exception identify_rate_equation(
        prob; popsize=200, optimizer=CMAEvolutionStrategyOpt(),
        n_restarts=1, maxtime=1.0, save_dir=mktempdir())
    @test_throws Exception identify_rate_equation(
        prob; verbose=-9, optimizer=CMAEvolutionStrategyOpt(),
        n_restarts=1, maxtime=1.0, save_dir=mktempdir())
end

@testset "save_dir holding only a progress.log is refused" begin
    # A run that crashed before its first CSV leaves only progress.log behind; a new
    # run there would append its log to the old one. The lean beam settings keep the
    # run short should the guard ever let it through.
    prob = _testhelper_uni_prob(NamedTuple)
    crashed_log = "EnzymeRates v$(pkgversion(EnzymeRates))\n" *
                  "Enumerating initial mechanisms…\n"
    mktempdir() do tmp
        write(joinpath(tmp, "progress.log"), crashed_log)
        @test_throws(
            ErrorException("save_dir already contains results (CSV files or " *
                "progress.log). Use an empty directory to avoid mixing results."),
            identify_rate_equation(prob; optimizer=CMAEvolutionStrategyOpt(),
                min_beam_width=1, loss_rel_threshold=1.0, loss_abs_threshold=0.0,
                max_param_count=3, n_cv_candidates=1, n_restarts=1, maxtime=1.0,
                save_dir=tmp))
        @test readdir(tmp) == ["progress.log"]
        @test read(joinpath(tmp, "progress.log"), String) == crashed_log
    end
end

@testset "n_cv_candidates below 1 is refused before any fit" begin
    # The lean beam settings keep the run short should the check ever let it through.
    prob = _testhelper_uni_prob(NamedTuple)
    mktempdir() do tmp
        @test_throws(ErrorException("n_cv_candidates must be ≥ 1; got 0"),
            identify_rate_equation(prob; optimizer=CMAEvolutionStrategyOpt(),
                min_beam_width=1, loss_rel_threshold=1.0, loss_abs_threshold=0.0,
                max_param_count=3, n_cv_candidates=0, n_restarts=1, maxtime=1.0,
                save_dir=tmp))
        @test isempty(readdir(tmp))       # raised before writing anything
    end
end

# ── §2 fit-dedup by key ──────────────────────────────────────────────────────
# A stub optimizer that counts `solve` invocations and returns a canned
# log-space optimum (`uval` for every coordinate), so a batch's fits can be
# counted exactly and the raw→rescale path exercised deterministically. With
# `infinite = true` every restart reports an infinite objective, as a fit with no
# finite loss does. `throwit = true` throws on every call, `throw_first = true` on the
# first call only; the error message names the call by its number.
mutable struct _testhelper_CountingStubOpt
    count::Int
    uval::Float64
    throwit::Bool
    infinite::Bool
    throw_first::Bool
end
_testhelper_CountingStubOpt(; uval=0.0, throwit=false, infinite=false, throw_first=false) =
    _testhelper_CountingStubOpt(0, uval, throwit, infinite, throw_first)
Optimization.allowsbounds(::_testhelper_CountingStubOpt) = true
function Optimization.SciMLBase.__solve(
        prob::Optimization.OptimizationProblem, opt::_testhelper_CountingStubOpt; kwargs...)
    opt.count += 1
    (opt.throwit || opt.throw_first && opt.count == 1) &&
        error("stub solver forced failure $(opt.count)")
    u = fill(opt.uval, length(prob.u0))
    cache = DefaultOptimizationCache(prob.f, prob.p)
    opt.infinite && return build_solution(cache, opt, u, Inf; retcode = ReturnCode.Failure)
    build_solution(cache, opt, u, prob.f(u, prob.p); retcode = ReturnCode.Success)
end

# Two allosteric uni-uni mechanisms that compile to distinct types yet render one rate
# equation: S and P bind only the active conformation (:OnlyA), so the inactive one keeps
# no form but E and the tag of the chemistry step (:NonequalAI or :OnlyA) never reaches
# the law. Written out directly so the tests are self-contained and cheap. The tests
# re-verify the collision at run time.
function _testhelper_dedup_twins()
    m1 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)       :: OnlyA
            E(S) <--> E(P)     :: NonequalAI
            E(P) ⇌ E + P       :: OnlyA
        end
    end)
    m2 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)       :: OnlyA
            E(S) <--> E(P)     :: OnlyA
            E(P) ⇌ E + P       :: OnlyA
        end
    end)
    return m1, m2
end

@testset "fit-dedup by key in _process_batch" begin
    m1, m2 = _testhelper_dedup_twins()
    em1 = EnzymeRates.compile_mechanism(m1)
    em2 = EnzymeRates.compile_mechanism(m2)
    c1, c2 = _testhelper_pass1_canonical.((m1, m2))
    # Preconditions: distinct structure and compiled types, one key and one law under
    # the canonical map, and the same fitted names.
    @test m1 != m2
    @test typeof(em1) != typeof(em2)
    @test c1.key == c2.key
    @test _testhelper_same_law_under_map(EnzymeRates._law(m1), EnzymeRates._law(m2),
                                         c1.order, c2.order;
                                         rng = Random.MersenneTwister(1))
    @test EnzymeRates.fitted_params(em1) == EnzymeRates.fitted_params(em2)

    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 1.0, 2.0], P = [0.1, 0.2, 0.1, 0.2])
    # The range stops below the stub fit's binding constants (5), so the fit still rises at
    # the range's upper edge and carries the edge flag.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    prob = IdentifyRateEquationProblem(rxn, data; Keq=2.0, peak_range=(0.1, 1.0))
    pair = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1, m2]

    memo = Dict{String, EnzymeRates._Fit}()
    opt = _testhelper_CountingStubOpt(; uval = log(5.0))
    # Each logged line is stored with the solve count at the moment it was logged.
    batch_log = Tuple{String,Int}[]
    entries, failures = EnzymeRates._process_batch(pair, prob;
        optimizer=opt, max_param_count=20, eq_complexity_filter=337, n_restarts=1,
        maxtime=1.0, memo, log = msg -> push!(batch_log, (msg, opt.count)))

    # The shared law is fit exactly ONCE (n_restarts=1 → one solve).
    @test opt.count == 1
    # The pre-fit summary counts one new fit and one inherited row, and is logged
    # before the fit runs.
    @test batch_log == [("1 new fits + 1 inherited + 0 skipped (already fit) + " *
                         "0 skipped (>20 params) + 0 skipped (>337 complexity)", 0)]
    @test length(entries) == 2
    @test isempty(failures)
    # loss + retcode belong to the law's fit → one key → identical.
    @test entries[1].loss == entries[2].loss
    @test entries[1].retcode === entries[2].retcode === :Success
    # Representative fit first (false); the duplicate is inherited (true), and both
    # name the representative as the fit's source.
    @test [e.row.fit_inherited for e in entries] == [false, true]
    @test [e.row.fit_source_id for e in entries] == fill(entries[1].row.mechanism_id, 2)

    # The memo stores the fit at canonical positions; every row reads it through its own
    # canonical order, so both rows carry the fit's values under their own names.
    fit = memo[c1.key]
    @test fit.source == m1
    fkeys = EnzymeRates.fitted_params(em1)
    for (e, c) in zip(entries, (c1, c2))
        @test keys(e.row.params) == fkeys
        @test Float64[e.row.params[x] for x in [c.order; :L]] == fit.values
        # The edge flag belongs to the law's fit, so each row copies it.
        @test e.row.peak_at_range_edge === fit.peak_at_range_edge
    end
    @test fit.peak_at_range_edge === true
    @test entries[1].row.params == entries[2].row.params
    # The peak scale: the copied params are the rescaled fit, not the raw
    # 5.0 the stub optimizer returned.
    @test !all(v -> v ≈ 5.0, entries[1].row.params)

    # Cross-batch memo hit: a later batch with the same key refits NOTHING.
    single = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1]
    reuse_log = String[]
    reused, _ = EnzymeRates._process_batch(single, prob;
        optimizer=opt, max_param_count=20, eq_complexity_filter=337, n_restarts=1,
        maxtime=1.0, memo, log = msg -> push!(reuse_log, msg))
    @test opt.count == 1                                # no new solve
    # The mechanism that was fitted reads its own fit in a later batch: not inherited.
    @test [e.row.fit_inherited for e in reused] == [false]
    @test reuse_log == ["0 new fits + 1 inherited + 0 skipped (already fit) + " *
                        "0 skipped (>20 params) + 0 skipped (>337 complexity)"]

    # A representative whose fit throws fails ALL the mechanisms with its key
    # (all-or-nothing per law).
    opt_bad = _testhelper_CountingStubOpt(; throwit=true)
    bad_entries, bad_failures = EnzymeRates._process_batch(pair, prob;
        optimizer=opt_bad, max_param_count=20, n_restarts=1, maxtime=1.0,
        memo = Dict{String, EnzymeRates._Fit}())
    @test isempty(bad_entries)
    @test length(bad_failures) == 2
    @test all(f -> f isa EnzymeRates.FitFailure, bad_failures)
    @test [f.mech for f in bad_failures] == [m1, m2]   # the mechanisms as handed in
end

# A stub optimizer that counts `solve` calls and returns the log-space point
# u_i = log(i + 1), so each fitted constant gets a value of its own and a row read through
# the wrong names gives another rate.
mutable struct _testhelper_SpreadStubOpt
    count::Int
end
_testhelper_SpreadStubOpt() = _testhelper_SpreadStubOpt(0)
Optimization.allowsbounds(::_testhelper_SpreadStubOpt) = true
function Optimization.SciMLBase.__solve(
        prob::Optimization.OptimizationProblem, opt::_testhelper_SpreadStubOpt; kwargs...)
    opt.count += 1
    u = [log(i + 1.0) for i in eachindex(prob.u0)]
    cache = DefaultOptimizationCache(prob.f, prob.p)
    build_solution(cache, opt, u, prob.f(u, prob.p); retcode = ReturnCode.Success)
end

# Two allosteric uni-bi mechanisms with one law and different fitted names: A and
# Q bind only the active conformation (:OnlyA), so the inactive one keeps free E alone, and
# P's binding at E(Q), which that conformation never reaches, is :EqualAI in the first and
# :OnlyA in the second. The first names its constant without a conformation
# (`K_EPQ_to_EQ_P`), the second as the active conformation's (`K_A_EPQ_to_EQ_P`), so the
# two sort it to different fitted positions. Returns the two and a problem on their
# reaction.
function _testhelper_unibi_pair()
    equal = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A
        products: P, Q
        catalytic_multiplicity: 2
        steps: begin
            E + A <--> E(A)          :: OnlyA
            E(A) <--> E(P, Q)        :: OnlyA
            E(Q) + P ⇌ E(P, Q)       :: EqualAI
            E + Q ⇌ E(Q)             :: OnlyA
        end
    end)
    only_active = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A
        products: P, Q
        catalytic_multiplicity: 2
        steps: begin
            E + A <--> E(A)          :: OnlyA
            E(A) <--> E(P, Q)        :: OnlyA
            E(Q) + P ⇌ E(P, Q)       :: OnlyA
            E + Q ⇌ E(Q)             :: OnlyA
        end
    end)
    rxn = @enzyme_reaction begin
        substrates: A[CN]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            A = [1.0, 2.0, 1.0, 2.0], P = [0.1, 0.2, 0.1, 0.2], Q = [0.3, 0.3, 0.4, 0.4])
    equal, only_active, IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
end

@testset "two same-key mechanisms read one fit under their own names" begin
    # Spec section 10: identify fits one representative per key and stores the fit at
    # canonical positions; the other mechanism reads it through its own labelling, so its
    # row carries its own names and both rows compute one rate.
    a, b, prob = _testhelper_unibi_pair()
    ca, cb = _testhelper_pass1_canonical.((a, b))
    fa, fb = EnzymeRates.fitted_params(a), EnzymeRates.fitted_params(b)
    @test ca.key == cb.key
    @test Set(fa) != Set(fb)
    # The canonical map sends a constant to another fitted position, so a fit copied by
    # position would hand the second row wrong values.
    @test [findfirst(==(x), fa) for x in ca.order] !=
          [findfirst(==(x), fb) for x in cb.order]
    opt = _testhelper_SpreadStubOpt()
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[a, b], prob;
        optimizer = opt, max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test opt.count == 1
    @test isempty(failures) && length(entries) == 2
    @test count(e -> e.row.fit_inherited, entries) == 1
    @test [e.key for e in entries] == [ca.key, ca.key]
    @test keys(entries[1].row.params) == fa && keys(entries[2].row.params) == fb
    rng = Random.MersenneTwister(10)
    for _ in 1:3
        concs = (A = exp(randn(rng)), P = exp(randn(rng)), Q = exp(randn(rng)))
        rates = [rate_equation(EnzymeRates.compile_mechanism(e.mech), concs,
                               merge(e.row.params, (Keq = 2.0, E_total = 1.0)))
                 for e in entries]
        @test isapprox(rates[1], rates[2]; rtol = 1e-12)
    end
end

@testset "a fit is stored by canonical position and read through each row's order" begin
    # The memo's values sit at the canonical positions of the representative's law (its
    # `_Canonical.order`, then L); each row's own canonical order reads them back.
    a, b, prob = _testhelper_unibi_pair()
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[a, b], prob;
        optimizer = _testhelper_SpreadStubOpt(), max_param_count = 20, memo,
        n_restarts = 1, maxtime = 1.0)
    fit = only(values(memo))
    @test fit.source == a && fit.source_id == entries[1].row.mechanism_id
    @test fit.loss == entries[1].loss == entries[2].loss
    for e in entries
        c = EnzymeRates._canonical(EnzymeRates._law(e.mech), e.row.mechanism_id)
        @test Float64[e.row.params[x] for x in [c.order; :L]] == fit.values
    end
end

@testset "_offer_cv! keeps one slot per key" begin
    pool = EnzymeRates.BatchEntry[]
    EnzymeRates._offer_cv!(pool, _testhelper_entry("a", 1.0), 2)
    EnzymeRates._offer_cv!(pool, _testhelper_entry("b", 2.0), 2)
    # A full pool: a repeat key with a lower loss takes its own slot, not the worst one.
    EnzymeRates._offer_cv!(pool, _testhelper_entry("a", 0.5), 2)
    @test [(e.key, e.loss) for e in pool] == [("a", 0.5), ("b", 2.0)]
    # A repeat key with a higher loss changes nothing; a new key evicts the worst.
    EnzymeRates._offer_cv!(pool, _testhelper_entry("a", 3.0), 2)
    EnzymeRates._offer_cv!(pool, _testhelper_entry("c", 1.5), 2)
    @test [(e.key, e.loss) for e in pool] == [("a", 0.5), ("c", 1.5)]
end

@testset "LOOCV fits each fit key's source once per fold" begin
    # Spec section 10: rows that read one fit share their fold losses, so LOOCV fits
    # that fit's source once per fold. Handed to model selection beside a second law, the
    # dedup twins (one key) cost two fits per fold, not three, and share their fold columns.
    m1, m2 = _testhelper_dedup_twins()
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        oligomeric_state: 2
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 1.0, 2.0], P = [0.1, 0.2, 0.1, 0.2])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[
            m1, m2, _testhelper_uni_seed], prob;
        optimizer = _testhelper_CountingStubOpt(; uval = log(5.0)), max_param_count = 20,
        memo, n_restarts = 1, maxtime = 1.0)
    @test length(unique(e.fit_key for e in entries)) == 2
    opt = _testhelper_CountingStubOpt(; uval = log(5.0))
    res = EnzymeRates._cv_model_selection(entries, memo, prob; optimizer = opt,
        se_threshold = 1.0, save_dir = mktempdir(), show_progress = false,
        n_restarts = 1, maxtime = 1.0)
    groups = unique(prob.data.group)
    @test opt.count == 2 * length(groups)
    folds(r) = [r[Symbol("cv_fold_$g")] for g in groups]
    twins = [r for r in eachrow(res.cv_results) if r.key == entries[1].key]
    @test length(twins) == 2 && folds(twins[1]) == folds(twins[2])
end

@testset "best_equation.csv reads back with its key and ids" begin
    # In a one-row CSV, CSV.jl guesses a column's type from one value: a
    # 16-hex id of digits and one `e` would read as a number, so the documented read names
    # the three id columns `String`; the 64-hex key reads back as text without help.
    prob = _testhelper_uni_prob(NamedTuple)
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[_testhelper_uni_seed],
        prob; optimizer = _testhelper_CountingStubOpt(; uval = 0.0), max_param_count = 20,
        memo, n_restarts = 1, maxtime = 1.0)
    e = only(entries)
    ids = (mechanism_id = "1234567e89012345", parent_id = "0123456789012345",
           fit_source_id = "1234567e89012345")
    pool = [EnzymeRates.BatchEntry(e.mech, e.n_params, e.loss, e.retcode, e.key, e.fit_key,
                                   merge(e.row, ids))]
    mktempdir() do dir
        EnzymeRates._cv_model_selection(pool, memo, prob;
            optimizer = _testhelper_CountingStubOpt(; uval = 0.0), se_threshold = 1.0,
            save_dir = dir, show_progress = false, n_restarts = 1, maxtime = 1.0)
        path = joinpath(dir, "best_equation.csv")
        plain = CSV.read(path, DataFrame)
        @test plain.key[1] isa AbstractString
        @test plain.key[1] == e.key && occursin(r"^[0-9a-f]{64}$", plain.key[1])
        @test plain.L_phantom[1] === false
        typed = CSV.read(path, DataFrame; types = Dict(:mechanism_id => String,
                                                       :parent_id => String,
                                                       :fit_source_id => String))
        @test (typed.mechanism_id[1], typed.parent_id[1], typed.fit_source_id[1]) ==
              Tuple(ids)
    end
end

@testset "a law at the leaf cap is fitted under its mechanism's key" begin
    # Five two-step cycles E → X → E, each taking up A and giving off P twice, alike up to
    # the label of X: every permutation of the cycles maps the law onto itself, so a
    # canonical labelling would compare 5! = 120 labellings and the search stops at the
    # leaf cap. The mechanism is fitted under a key of its own, the SHA-256 of
    # "mechanism|" and its mechanism_id, and the pre-fit summary names it. (The key search
    # takes about 0.4 s here once compiled, so it runs once, inside the batch.)
    rxn = @enzyme_reaction begin
        substrates: A[C]
        products: P[C]
    end
    pinwheel = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> F + P
            F + A <--> E + P
            E + A <--> G + P
            G + A <--> E + P
            E + A <--> H + P
            H + A <--> E + P
            E + A <--> J + P
            J + A <--> E + P
            E + A <--> M + P
            M + A <--> E + P
        end
    end)
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            A = [1.0, 2.0, 3.0, 4.0], P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 10.0)
    id = EnzymeRates._mechanism_id(sprint(show, pinwheel))
    batch_log = String[]
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[pinwheel], prob;
        optimizer = _testhelper_CountingStubOpt(), max_param_count = 20,
        eq_complexity_filter = 337, n_restarts = 1, maxtime = 1.0,
        log = msg -> push!(batch_log, msg))
    @test isempty(failures)
    @test only(entries).key == only(entries).row.key ==
          bytes2hex(EnzymeRates.sha256("mechanism|" * id))
    @test batch_log == ["1 new fits + 0 inherited + 0 skipped (already fit) + " *
                        "0 skipped (>20 params) + 0 skipped (>337 complexity)\n" *
                        "  leaf cap reached for $id; fitted under its own key"]
end

# The merged uni-uni and its K-type child at one catalytic subunit
# (test_mechanism_enumeration.jl, "expand_mechanisms on a merged uni-uni"): the child's
# inactive conformation is free E alone and its certificate maps it onto the parent's
# law, so it is an L phantom whose law at L = 0 has the parent's key.
function _testhelper_phantom_pair()
    parent = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    phantom = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(P)    :: OnlyA
            E + P ⇌ E(P)       :: OnlyA
        end
    end)
    parent, phantom
end

@testset "a phantom inherits the fit of whoever holds its L→0 key" begin
    # Spec section 13: the phantom is not fitted. It reads the fit of its law at
    # L = 0 from the parent, which holds that key and is fitted in the same batch, with
    # L = 0, the parent's loss and its own parameter count; its row names the parent as
    # the fit's source, and the two rates agree. It then joins the frontier and the cv pool,
    # and a sweep whose width floor admits one key expands it.
    parent, phantom = _testhelper_phantom_pair()
    prob = _testhelper_uni_prob(NamedTuple)
    @test EnzymeRates._is_phantom(EnzymeRates._law(phantom))
    opt = _testhelper_SpreadStubOpt()
    batch_log = String[]
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[phantom, parent],
        prob; optimizer = opt, max_param_count = 20, eq_complexity_filter = 337,
        n_restarts = 1, maxtime = 1.0, log = msg -> push!(batch_log, msg))
    @test opt.count == 1 && isempty(failures)
    @test batch_log == ["1 new fits + 1 inherited + 0 skipped (already fit) + " *
                        "0 skipped (>20 params) + 0 skipped (>337 complexity)"]
    p = only(filter(e -> e.mech == phantom, entries))
    f = only(filter(e -> e.mech == parent, entries))
    @test p.row.L_phantom && !f.row.L_phantom
    @test p.row.params.L == 0.0
    @test p.loss == f.loss && p.retcode === f.retcode
    @test (p.n_params, f.n_params) == (3, 2)
    @test p.fit_key == f.key != p.key
    @test p.row.fit_inherited && !f.row.fit_inherited
    @test p.row.fit_source_id == f.row.mechanism_id == f.row.fit_source_id
    @test p.row.peak_at_range_edge == f.row.peak_at_range_edge
    rng = Random.MersenneTwister(13)
    for _ in 1:3
        concs = (S = exp(randn(rng)), P = exp(randn(rng)))
        _testhelper_rate(e) =
            rate_equation(EnzymeRates.compile_mechanism(e.mech), concs,
                          merge(e.row.params, (Keq = 10.0, E_total = 1.0)))
        @test isapprox(_testhelper_rate(p), _testhelper_rate(f); rtol = 1e-12)
    end
    frontier = Dict{Int, Vector{EnzymeRates.BatchEntry}}()
    cv_pool = Dict{Int, Vector{EnzymeRates.BatchEntry}}()
    best = Dict{Int, Float64}()
    EnzymeRates._ingest!(frontier, cv_pool, best, entries; n_cv_candidates = 5)
    @test p in frontier[3] && p in cv_pool[3]
    parents = EnzymeRates._sweep_parents!(frontier, 3, Dict{Int, Int}(), best,
        Set{String}(); loss_rel_threshold = 1.3, loss_abs_threshold = 0.001,
        loss_parsimony_threshold = 0.99, min_beam_width = 1)
    @test p in parents
end

@testset "a phantom whose L→0 key has no fit is fitted unflagged" begin
    # Census #3268, the RE→SS flip of B's binding in the K-type of an ordered
    # bi-bi, is an L phantom whose law at L = 0 no mechanism of the batch or the memo
    # holds (the plain flip move never emits it: `_chain_flank_groups` keeps B's group
    # at rapid equilibrium). It is fitted as an ordinary candidate under its own key,
    # unflagged, its L fitted.
    m3268 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)              :: OnlyA
            E + Q ⇌ E(Q)              :: OnlyA
            E(A) + B <--> E(A, B)     :: EqualAI
            E(A, B) <--> E(P, Q)      :: OnlyA
            E(Q) + P ⇌ E(P, Q)        :: EqualAI
        end
    end)
    @test EnzymeRates._is_phantom(EnzymeRates._law(m3268))
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            A = [1.0, 2.0, 1.0, 2.0], B = [0.5, 0.5, 1.0, 1.0],
            P = [0.1, 0.2, 0.1, 0.2], Q = [0.3, 0.3, 0.4, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
    opt = _testhelper_SpreadStubOpt()
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m3268], prob;
        optimizer = opt, max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    e = only(entries)
    @test opt.count == 1 && isempty(failures)
    @test !e.row.L_phantom && !e.row.fit_inherited
    @test e.fit_key == e.key
    @test e.n_params == 7 && e.row.params.L != 0.0
end

@testset "a phantom whose source lacks its constants is fitted itself" begin
    # The uni-uni at one catalytic subunit with S's and P's bindings :EqualAI, the
    # chemistry :OnlyA and a regulatory site whose one ligand R is :EqualAI. The site's
    # factor is the same in both conformations, so the certificate holds, and the law at
    # L = 0 is the law of the plain uni-uni beside it in the batch, which lacks R's
    # constant. Reading that fit would leave R's constant without a value, so the variant
    # is fitted itself under its own key, unflagged, its L fitted. The enumerator never
    # writes an all-:EqualAI site, so the variant is written by hand.
    plain = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
        end
    end)
    regulated = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        catalytic_multiplicity: 1
        steps: begin
            E + S ⇌ E(S)       :: EqualAI
            E(S) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)       :: EqualAI
        end
        regulatory_site(1): R::EqualAI
    end)
    A = EnzymeRates._law(regulated)
    @test EnzymeRates._is_phantom(A)
    @test EnzymeRates._canonical(A.active, "regulated").key ==
          EnzymeRates._canonical(EnzymeRates._law(plain), "plain").key
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0], P = [0.1, 0.2, 0.3, 0.4], R = [0.5, 1.0, 0.5, 1.0])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 10.0)
    opt = _testhelper_SpreadStubOpt()
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[regulated, plain],
        prob; optimizer = opt, max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test opt.count == 2 && isempty(failures)
    e = only(filter(e -> e.mech == regulated, entries))
    @test !e.row.L_phantom && !e.row.fit_inherited
    @test e.fit_key == e.key
    @test e.n_params == 5 && e.row.params.L != 0.0
end

@testset "a phantom reads its source's fold losses" begin
    # Equal families give equal fold optima, so LOOCV fits the source once per
    # fold and the phantom takes those fold losses. Beside its source the phantom, one
    # parameter larger, loses the 1-SE tie.
    parent, phantom = _testhelper_phantom_pair()
    prob = _testhelper_uni_prob(NamedTuple)
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[parent, phantom],
        prob; optimizer = _testhelper_SpreadStubOpt(), max_param_count = 20, memo,
        n_restarts = 1, maxtime = 1.0)
    opt = _testhelper_CountingStubOpt(; uval = 0.0)
    res = EnzymeRates._cv_model_selection(entries, memo, prob; optimizer = opt,
        se_threshold = 1.0, save_dir = mktempdir(), show_progress = false,
        n_restarts = 1, maxtime = 1.0)
    groups = unique(prob.data.group)
    @test opt.count == length(groups)
    folds(r) = [r[Symbol("cv_fold_$g")] for g in groups]
    @test res.cv_results.n_params == [2, 3]
    @test folds(res.cv_results[1, :]) == folds(res.cv_results[2, :])
    @test res.best === EnzymeRates.compile_mechanism(parent)
end

@testset "a phantom chosen by the 1-SE rule reports L = 0 and evaluates" begin
    # The pool holds the phantom without its source, so the 1-SE rule
    # selects it. best_equation.csv reports L = 0, the flag and the source's id, and its
    # parameters, read back from the CSV, give the source's rates through the phantom's
    # own rate equation; rescale_parameter_values accepts them.
    parent, phantom = _testhelper_phantom_pair()
    prob = _testhelper_uni_prob(NamedTuple)
    memo = Dict{String, EnzymeRates._Fit}()
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[parent, phantom],
        prob; optimizer = _testhelper_SpreadStubOpt(), max_param_count = 20, memo,
        n_restarts = 1, maxtime = 1.0)
    p = only(filter(e -> e.mech == phantom, entries))
    f = only(filter(e -> e.mech == parent, entries))
    mktempdir() do dir
        res = EnzymeRates._cv_model_selection([p], memo, prob;
            optimizer = _testhelper_CountingStubOpt(; uval = 0.0), se_threshold = 1.0,
            save_dir = dir, show_progress = false, n_restarts = 1, maxtime = 1.0)
        @test res.best === EnzymeRates.compile_mechanism(phantom)
        best = CSV.read(joinpath(dir, "best_equation.csv"), DataFrame;
            types = Dict(:mechanism_id => String, :parent_id => String,
                         :fit_source_id => String))
        @test best.L[1] == 0.0 && best.L_phantom[1] === true
        @test best.fit_source_id[1] == f.row.mechanism_id
        fp = EnzymeRates.fitted_params(res.best)
        params = merge(NamedTuple{fp}(Tuple(Float64(best[1, x]) for x in fp)),
                       (Keq = prob.Keq, E_total = 1.0))
        source = merge(f.row.params, (Keq = prob.Keq, E_total = 1.0))
        rng = Random.MersenneTwister(14)
        for _ in 1:3
            concs = (S = exp(randn(rng)), P = exp(randn(rng)))
            v = rate_equation(res.best, concs, params)
            v0 = rate_equation(EnzymeRates.compile_mechanism(parent), concs, source)
            @test isfinite(v) && isapprox(v, v0; rtol = 1e-12)
        end
        scaled, spread = rescale_parameter_values(res.best, params)
        @test spread === nothing && scaled.L == 0.0
    end
end

@testset "a phantom whose source's fit fails" begin
    # When the source's fit throws, in the phantom's batch or an earlier one, the phantom is
    # fitted itself, unflagged; when that fit throws too, it fails with its own error. When
    # the source's fit reaches no finite loss, the phantom's row carries the source's loss
    # and retcode, and `_ingest!` keeps both out of the frontier and the cv pool.
    parent, phantom = _testhelper_phantom_pair()
    prob = _testhelper_uni_prob(NamedTuple)
    batch = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[parent, phantom]
    memo = Dict{String, EnzymeRates._Fit}()
    lines = String[]
    opt = _testhelper_CountingStubOpt(; throw_first = true)
    entries, failures = EnzymeRates._process_batch(batch, prob; optimizer = opt,
        max_param_count = 20, memo, n_restarts = 1, maxtime = 1.0,
        log = msg -> push!(lines, msg))
    @test [f.mech for f in failures] == [parent]
    @test opt.count == 2
    e = only(entries)
    @test e.mech == phantom
    @test e.row.L_phantom === false
    @test e.row.fit_inherited === false
    @test e.row.fit_source_id == e.row.mechanism_id
    @test e.fit_key == e.key
    @test isfinite(e.loss)
    @test haskey(memo, e.key)
    @test lines[2:end] == ["source fit threw for L phantom $(e.row.mechanism_id); " *
                           "it reads its own key's fit"]
    opt = _testhelper_CountingStubOpt(; throwit = true)
    entries, failures = EnzymeRates._process_batch(batch, prob; optimizer = opt,
        max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test isempty(entries)
    @test [f.mech for f in failures] == [parent, phantom]
    @test opt.count == 2
    @test occursin("stub solver forced failure 1", failures[1].error)
    @test occursin("stub solver forced failure 2", failures[2].error)
    @test ismissing(EnzymeRates._failure_row(failures[2]).L_phantom)
    memo = Dict{String, EnzymeRates._Fit}()
    EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[parent], prob;
        optimizer = _testhelper_CountingStubOpt(; throwit = true), max_param_count = 20,
        memo, n_restarts = 1, maxtime = 1.0)
    @test isempty(memo)
    entries, _ = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[phantom], prob;
        optimizer = _testhelper_CountingStubOpt(), max_param_count = 20, memo,
        n_restarts = 1, maxtime = 1.0)
    e = only(entries)
    @test e.row.L_phantom === false
    @test e.row.fit_inherited === false
    @test e.row.fit_source_id == e.row.mechanism_id
    @test e.fit_key == e.key
    opt = _testhelper_CountingStubOpt(; infinite = true)
    entries, failures = EnzymeRates._process_batch(batch, prob; optimizer = opt,
        max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test opt.count == 1 && isempty(failures)
    @test all(e -> e.loss == Inf && e.retcode === :NoFiniteLoss, entries)
    @test only(filter(e -> e.mech == phantom, entries)).row.L_phantom
    frontier = Dict{Int, Vector{EnzymeRates.BatchEntry}}()
    pool = Dict{Int, Vector{EnzymeRates.BatchEntry}}()
    best = Dict{Int, Float64}()
    @test isempty(EnzymeRates._ingest!(frontier, pool, best, entries; n_cv_candidates = 5))
    @test isempty(frontier) && isempty(pool) && isempty(best)
end

@testset "phantoms demoted after their source's fit throws share one fit per key" begin
    # Two L phantoms of the uni seed, tag-set siblings with one law and one key. When the
    # seed's fit throws, both are demoted and one fit serves their key; a key that `memo`
    # already holds is read, not fitted again, so a phantom with a usable fit never fails.
    m1 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)       :: OnlyA
            E(S) <--> E(P)     :: NonequalAI
            E(P) ⇌ E + P       :: OnlyA
        end
    end)
    m2 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)       :: OnlyA
            E(S) <--> E(P)     :: OnlyA
            E(P) ⇌ E + P       :: OnlyA
        end
    end)
    seed = _testhelper_uni_seed
    prob = _testhelper_uni_prob(NamedTuple)
    batch(ms...) = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[ms...]
    # One batch: the seed's fit throws, and m1 alone is fitted for both phantoms.
    memo = Dict{String, EnzymeRates._Fit}()
    lines = String[]
    opt = _testhelper_CountingStubOpt(; throw_first = true)
    entries, failures = EnzymeRates._process_batch(batch(seed, m1, m2), prob;
        optimizer = opt, max_param_count = 20, memo, n_restarts = 1, maxtime = 1.0,
        log = msg -> push!(lines, msg))
    @test opt.count == 2
    @test [f.mech for f in failures] == [seed]
    @test [e.mech for e in entries] == [m1, m2]
    e1, e2 = entries
    @test e1.key == e2.key
    @test collect(keys(memo)) == [e1.key]
    @test e1.row.fit_inherited === false
    @test e2.row.fit_inherited === true
    @test e2.row.fit_source_id == e1.row.mechanism_id
    @test lines[2:end] == [join(("source fit threw for L phantom $(e.row.mechanism_id); " *
                                 "it reads its own key's fit" for e in (e1, e2)), "\n  ")]
    # m1's key is fitted in an earlier batch; m2, demoted beside the seed, reads that fit,
    # which stays as it was.
    memo = Dict{String, EnzymeRates._Fit}()
    EnzymeRates._process_batch(batch(m1), prob; optimizer = _testhelper_CountingStubOpt(),
        max_param_count = 20, memo, n_restarts = 1, maxtime = 1.0)
    before = memo[e1.key]
    opt = _testhelper_CountingStubOpt(; throw_first = true)
    entries, failures = EnzymeRates._process_batch(batch(seed, m2), prob;
        optimizer = opt, max_param_count = 20, memo, n_restarts = 1, maxtime = 1.0)
    @test opt.count == 1
    @test [f.mech for f in failures] == [seed]
    e = only(entries)
    @test e.mech == m2
    @test e.row.fit_inherited === true
    @test e.row.fit_source_id == e1.row.mechanism_id
    @test memo[e1.key] === before
    # As above with every fit throwing: m2 still reads m1's fit and does not fail.
    opt = _testhelper_CountingStubOpt(; throwit = true)
    entries, failures = EnzymeRates._process_batch(batch(seed, m2), prob;
        optimizer = opt, max_param_count = 20, memo, n_restarts = 1, maxtime = 1.0)
    @test opt.count == 1
    @test [f.mech for f in failures] == [seed]
    @test only(entries).mech == m2
end


# LDH renaming-dup pair (same graph, tied kinetic-group split): merged form (8 groups) vs
# split form (9 groups). The split form's extra binding K is tied straight back by a
# Wegscheider cycle, so the two forms have one law: one key, one fit.
@testset "renaming-dup pair: one key, one fit" begin
    # Merged, 8 kinetic groups.
    m1 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: NADH, Pyruvate
        products:   Lactate, NAD
        steps: begin
            (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD),
             E(NADH) + Lactate ⇌ E(Lactate, NADH))
            (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD))
            (E + NADH ⇌ E(NADH), E(Lactate) + NADH ⇌ E(Lactate, NADH))
            E + Pyruvate ⇌ E(Pyruvate)
            (E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
             E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
            E(NADH, Pyruvate) <--> E(Lactate, NAD)
            E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)
            E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)
        end
    end)
    # Split, 9 groups — same graph, Wegscheider-tied.
    m2 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: NADH, Pyruvate
        products:   Lactate, NAD
        steps: begin
            (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD))
            (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD))
            (E + NADH ⇌ E(NADH), E(Lactate) + NADH ⇌ E(Lactate, NADH))
            E + Pyruvate ⇌ E(Pyruvate)
            (E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate),
             E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate))
            E(NADH) + Lactate ⇌ E(Lactate, NADH)
            E(NADH, Pyruvate) <--> E(Lactate, NAD)
            E(Pyruvate) + NADH ⇌ E(NADH, Pyruvate)
            E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate)
        end
    end)
    @test m1 != m2
    @test EnzymeRates._independent_param_count(m1) ==
          EnzymeRates._independent_param_count(m2)
    c1, c2 = _testhelper_pass1_canonical.((m1, m2))
    @test c1.key == c2.key
    @test _testhelper_same_law_under_map(EnzymeRates._law(m1), EnzymeRates._law(m2),
                                         c1.order, c2.order;
                                         rng = Random.MersenneTwister(2))
    # identify fits the pair once; the second form inherits the first's fit.
    rxn = @enzyme_reaction begin
        substrates: NADH[N1H1], Pyruvate[C1]
        products:   Lactate[C1H1], NAD[N1]
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            NADH = [1.0, 2.0, 1.0, 2.0], Pyruvate = [0.5, 0.5, 1.0, 1.0],
            Lactate = [0.1, 0.2, 0.1, 0.2], NAD = [0.3, 0.3, 0.4, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
    opt = _testhelper_CountingStubOpt(; uval = 0.0)
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1, m2], prob;
        optimizer = opt, max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test opt.count == 1 && isempty(failures)
    @test [e.row.fit_inherited for e in entries] == [false, true]
    @test [e.key for e in entries] == [c1.key, c1.key]
end


# LDH ALLOSTERIC split/merge pair: split form (8 cat groups) vs merged (7). The cycles tie
# the split form's extra constant back, so both are one law, with one key.
@testset "allosteric split/merge pair: the split form is the merged form's model" begin
    # Split form: 8 catalytic groups.
    am1 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products:   Lactate, NAD
        catalytic_multiplicity: 4
        steps: begin
            (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD))   :: EqualAI
            (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD),
             E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate))                           :: EqualAI
            (E + NADH <--> E(NADH), E(Lactate) + NADH <--> E(Lactate, NADH),
             E(Pyruvate) + NADH <--> E(NADH, Pyruvate))                      :: OnlyA
            E + Pyruvate ⇌ E(Pyruvate)                                       :: NonequalAI
            E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate)                             :: EqualAI
            E(NADH) + Lactate ⇌ E(Lactate, NADH)                             :: EqualAI
            E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)                           :: EqualAI
            E(NADH, Pyruvate) <--> E(Lactate, NAD)                           :: OnlyA
        end
    end)
    # Merged form: 7 groups.
    am2 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products:   Lactate, NAD
        catalytic_multiplicity: 4
        steps: begin
            (E + Lactate ⇌ E(Lactate), E(NAD) + Lactate ⇌ E(Lactate, NAD),
             E(NADH) + Lactate ⇌ E(Lactate, NADH))                           :: EqualAI
            (E + NAD ⇌ E(NAD), E(Lactate) + NAD ⇌ E(Lactate, NAD),
             E(Pyruvate) + NAD ⇌ E(NAD, Pyruvate))                           :: EqualAI
            (E + NADH <--> E(NADH), E(Lactate) + NADH <--> E(Lactate, NADH),
             E(Pyruvate) + NADH <--> E(NADH, Pyruvate))                      :: OnlyA
            E + Pyruvate ⇌ E(Pyruvate)                                       :: NonequalAI
            E(NAD) + Pyruvate ⇌ E(NAD, Pyruvate)                             :: EqualAI
            E(NADH) + Pyruvate ⇌ E(NADH, Pyruvate)                           :: EqualAI
            E(NADH, Pyruvate) <--> E(Lactate, NAD)                           :: OnlyA
        end
    end)
    @test am1 != am2
    # one law under the canonical map, so one key
    c1, c2 = _testhelper_pass1_canonical.((am1, am2))
    @test c1.key == c2.key
    @test _testhelper_same_law_under_map(EnzymeRates._law(am1), EnzymeRates._law(am2),
                                         c1.order, c2.order;
                                         rng = Random.MersenneTwister(3))
    # The split constant is tied back: the square E, E(Lactate), E(Lactate, NADH),
    # E(NADH), closed by the shared NADH group, makes K_ELactateNADH_to_ENADH_Lactate equal
    # the free-enzyme Lactate constant in the active conformation, and the inactive one
    # binds no NADH. Both forms fit the same constants, as many as the test's own count
    # of the dimension.
    @test EnzymeRates._independent_param_count(am1) ==
          EnzymeRates._independent_param_count(am2) == 7
    @test _testhelper_declared_dimension(am1) == _testhelper_declared_dimension(am2) == 7
    @test EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(am1)) ==
          EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(am2))
end


@testset "_process_batch failures report the ORIGINAL mechanism" begin
    # A mechanism whose derivation throws: T binds free E only in a dead end, so the
    # catalytic cycle turns S alone into P, and the derivation's thermodynamic-cycle
    # check (`_solve`, which `_process_batch` calls after the guard) raises
    # "Cycle 1 produces metabolite change not proportional to net reaction". The
    # FitFailure must carry the ORIGINAL `m0`.
    rxn_bad = @enzyme_reaction begin
        substrates: S[C], T[N]
        products:   P[CN]
    end
    m_bad = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S, T
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
            E + T ⇌ E(T)
        end
    end)
    data_bad = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
                S = [1.0, 2.0, 1.0, 2.0], T = [0.5, 0.5, 1.0, 1.0],
                P = [0.1, 0.2, 0.1, 0.2])
    prob_bad = IdentifyRateEquationProblem(rxn_bad, data_bad; Keq=2.0)
    e1, f1 = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m_bad],
        prob_bad; optimizer=_testhelper_CountingStubOpt(), max_param_count=20,
        n_restarts=1, maxtime=1.0, memo=Dict{String,EnzymeRates._Fit}())
    @test isempty(e1)
    @test length(f1) == 1 && f1[1] isa EnzymeRates.FitFailure
    @test f1[1].mech == m_bad                     # the mechanism as handed in
end

@testset "PASS1 runs the guard before it derives a law" begin
    # PASS1 derives the law from the concrete mechanism, so it runs the guard first, as
    # every derivation entry does. Rapid-equilibrium steps that turn S into P beside
    # a steady-state route fail class 2: PASS1 records the guard's message, not the
    # engine's internal error, and fits nothing.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    beside = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) ⇌ E(P)
            E + P ⇌ E(P)
            E(S) <--> F(P)
            F(P) <--> E + P
        end
    end
    data = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 1.0, 2.0], P = [0.1, 0.2, 0.1, 0.2])
    prob = IdentifyRateEquationProblem(rxn, data; Keq = 2.0)
    opt = _testhelper_CountingStubOpt()
    entries, failures = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[beside], prob;
        optimizer = opt, max_param_count = 20, n_restarts = 1, maxtime = 1.0)
    @test isempty(entries) && opt.count == 0
    @test only(failures).mech == beside
    @test occursin("E + S ⇌ E(S) and E(S) ⇌ E(P) close a cycle", only(failures).error)
    @test !occursin("internal error", only(failures).error)
end

@testset "_expand_parents records an expansion error instead of aborting" begin
    # expand_mechanisms asserts its input conserves atoms; this mechanism's
    # chemistry step does not (T binds free E only in a dead end, so the step makes
    # P's N from nothing), and the assertion raises. _expand_parents must catch that
    # and return the parent as a FitFailure (so the beam records it in CSV and
    # continues), not propagate and abort the search.
    rxn_bad = @enzyme_reaction begin
        substrates: S[C], T[N]
        products:   P[CN]
    end
    m_bad = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S, T
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E + P ⇌ E(P)
            E + T ⇌ E(T)
        end
    end)
    # Sanity: expansion of this mechanism genuinely raises.
    @test_throws ErrorException EnzymeRates.expand_mechanisms(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m_bad],
        rxn_bad)
    entry(m) = EnzymeRates.BatchEntry(m, 0, 0.0, :Success, "", "", (mechanism="M",))
    kids, parent_of, failures = EnzymeRates._expand_parents([entry(m_bad)], rxn_bad)
    @test isempty(kids) && isempty(parent_of)
    failure = only(failures)
    @test failure isa EnzymeRates.FitFailure
    @test failure.mech == m_bad                    # the ORIGINAL parent
    # The recorded error is the assertion's: the chemistry step ES → EP makes an N.
    @test occursin("atom-non-conserving step ES → EP", failure.error)
    @test occursin("= Dict(:N => 1)", failure.error)
    # A well-formed parent, a seed of the reaction, expands with no failure.
    good = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S, T
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) + T ⇌ E(S, T)
            E(S, T) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    gkids, _, gfail = EnzymeRates._expand_parents([entry(good)], rxn_bad)
    @test isempty(gfail) && !isempty(gkids)
end

@testset "_expand_parents matches each parent's expand_mechanisms" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        competitive_inhibitors: I
    end
    M = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}
    # The seed's children as parents: their own children overlap, so the
    # first-parent dedup and attribution are exercised.
    parents = EnzymeRates.expand_mechanisms(M[EnzymeRates.init_mechanisms(rxn)...], rxn)
    to_expand = [EnzymeRates.BatchEntry(m, 3, 0.0, :Success, "", "", (mechanism="M",))
                 for m in parents]

    children, parent_of, failures = EnzymeRates._expand_parents(to_expand, rxn)
    per_parent = [EnzymeRates.expand_mechanisms(M[pe.mech], rxn) for pe in to_expand]
    @test isempty(failures)
    # Children in first-parent order, each once; several parents share a child.
    @test children == unique(reduce(vcat, per_parent))
    @test length(children) < sum(length, per_parent)
    # Each child maps to the first parent that produced it.
    @test keys(parent_of) == Set(children)
    for child in children
        @test parent_of[child] === to_expand[findfirst(kids -> child in kids, per_parent)]
    end
end

@testset "_base_tier expands degenerate seeds instead of fitting them" begin
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    seeds = unique!(collect(EnzymeRates.init_mechanisms(rxn)))
    degenerate = filter(EnzymeRates._degenerate, seeds)
    @test length(degenerate) == 7                  # the seven ping-pong seeds
    base, failures = EnzymeRates._base_tier(seeds, rxn)
    @test isempty(failures)
    @test !any(EnzymeRates._degenerate, base)
    cures = unique!([c for m in degenerate for c in EnzymeRates._expand_re_to_ss(m)
                     if !EnzymeRates._degenerate(c)])
    @test length(cures) == 21                      # three per ping-pong seed
    # Only a flip changes which steps are at rapid equilibrium, so every other move leaves
    # a degenerate seed's child degenerate.
    kids = EnzymeRates.expand_mechanisms(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[degenerate...], rxn)
    @test Set(filter(!EnzymeRates._degenerate, kids)) == Set(cures)
    @test Set(base) == Set(vcat(filter(!EnzymeRates._degenerate, seeds), cures))
    @test length(base) == 278
end

@testset "_base_tier replaces a degenerate seed by the flips that cure it" begin
    # The undecorated ping-pong seed. Its steady-state chemistry E(A) → E(P; res) is a
    # qualifying chain, so the flip rule keeps its flanks, the A and P bindings, at rapid
    # equilibrium. The rapid-equilibrium chemistry E(B; res) ⇌ E(Q) is left by the B and
    # the Q binding, both rapid equilibrium, so the seed is degenerate. Each of the other
    # three groups alone raises the segment count, and each flip cures the seed: it turns
    # that chemistry steady state, or one of its two exits.
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: A[CX], B[N]
        products: P[C], Q[NX]
    end
    seed = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E + Q ⇌ E(Q)
        end
    end)
    chemistry_flipped = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) <--> E(Q)
            E + Q ⇌ E(Q)
        end
    end)
    b_flipped = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B <--> E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E + Q ⇌ E(Q)
        end
    end)
    q_flipped = ER.Mechanism(@enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E + Q <--> E(Q)
        end
    end)
    @test ER._degenerate(seed)
    base, failures = ER._base_tier([seed], rxn)
    @test isempty(failures)
    @test length(base) == 3
    @test Set(base) == Set([chemistry_flipped, b_flipped, q_flipped])
end

@testset "_base_tier records a degenerate seed's expansion error" begin
    # The chemistry step E(S) ⇌ E(P) makes P's N from nothing (T binds free E only in a
    # dead end), so expand_mechanisms' atom-conservation assertion raises. With the P
    # binding steady state the seed has no maximal rate in products (turning that
    # binding rapid equilibrium closes an all-RE turnover cycle), so it is degenerate and
    # expanded.
    rxn_bad = @enzyme_reaction begin
        substrates: S[C], T[N]
        products:   P[CN]
    end
    m_bad = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S, T
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) ⇌ E(P)
            E + P <--> E(P)
            E + T ⇌ E(T)
        end
    end)
    @test EnzymeRates._degenerate(m_bad)
    base, failures = EnzymeRates._base_tier([m_bad], rxn_bad)
    @test isempty(base)
    @test length(failures) == 1 && failures[1] isa EnzymeRates.FitFailure
    @test failures[1].mech == m_bad
    # The recorded error is the assertion's: the chemistry step ES → EP makes an N.
    @test occursin("atom-non-conserving step ES → EP", failures[1].error)
    @test occursin("= Dict(:N => 1)", failures[1].error)
end

@testset "_base_tier drops a flip child whose rate law the guard rejects" begin
    # The guard admits the seed. Its chemistry E(S) ⇌ E(P) sits in rapid equilibrium with
    # both sides, so the seed is degenerate and expanded, into six flips. Three of them
    # (the S binding, the chemistry, or the pair E(P) ⇌ F(P) and F + P ⇌ E(P) made
    # steady state) leave the :OnlyA step E <--> F bridging forms of the inactive
    # conformation at zero P. That is admission's face condition, which the flip move's
    # own check (`_faces_pass`, classes 2 and 3) does not run, so `_base_tier`'s guard
    # filter drops them. Two of the other flips keep the chemistry in rapid equilibrium
    # with both sides, which makes them degenerate. The base tier keeps the last flip.
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    seed = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) ⇌ E(P)         :: EqualAI
            F + P ⇌ E(P)        :: EqualAI
            F <--> E            :: OnlyA
            F + S ⇌ F(S)        :: EqualAI
            F(S) <--> F(P)      :: EqualAI
            F + P ⇌ F(P)        :: EqualAI
            F(P) ⇌ E(P)         :: EqualAI
        end
    end)
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    dropped = [
        @_testhelper_unguarded(@allosteric_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)     :: OnlyA
                E(S) ⇌ E(P)         :: EqualAI
                F + P ⇌ E(P)        :: EqualAI
                F <--> E            :: OnlyA
                F + S ⇌ F(S)        :: EqualAI
                F(S) <--> F(P)      :: EqualAI
                F + P ⇌ F(P)        :: EqualAI
                F(P) ⇌ E(P)         :: EqualAI
            end
        end),
        @_testhelper_unguarded(@allosteric_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)        :: OnlyA
                E(S) <--> E(P)      :: EqualAI
                F + P ⇌ E(P)        :: EqualAI
                F <--> E            :: OnlyA
                F + S ⇌ F(S)        :: EqualAI
                F(S) <--> F(P)      :: EqualAI
                F + P ⇌ F(P)        :: EqualAI
                F(P) ⇌ E(P)         :: EqualAI
            end
        end),
        @_testhelper_unguarded(@allosteric_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)        :: OnlyA
                E(S) ⇌ E(P)         :: EqualAI
                F + P <--> E(P)     :: EqualAI
                F <--> E            :: OnlyA
                F + S ⇌ F(S)        :: EqualAI
                F(S) <--> F(P)      :: EqualAI
                F + P ⇌ F(P)        :: EqualAI
                F(P) <--> E(P)      :: EqualAI
            end
        end)]
    kept = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)        :: OnlyA
            E(S) ⇌ E(P)         :: EqualAI
            F + P <--> E(P)     :: EqualAI
            F <--> E            :: OnlyA
            F + S ⇌ F(S)        :: EqualAI
            F(S) <--> F(P)      :: EqualAI
            F + P <--> F(P)     :: EqualAI
            F(P) ⇌ E(P)         :: EqualAI
        end
    end)
    @test ER._is_well_posed(seed)
    @test ER._degenerate(seed)
    kids = ER._expand_re_to_ss(seed)
    @test length(kids) == 6
    @test Set(k for k in kids if !ER._is_well_posed(k)) == Set(dropped)
    @test all(k -> ER._guard(k).fault === ER._BRIDGE, dropped)
    # The guard filter drops them, not the degeneracy test.
    @test !any(ER._degenerate, dropped)
    @test ER._guard_message(dropped[1], ER._guard(dropped[1])) ==
          "AllostericMechanism: at zero P, the :OnlyA step E <--> F bridges forms of " *
          "the inactive conformation that its other steps keep apart, so admission " *
          "cannot show that the limit does not depend on how fast that step vanishes."
    base, failures = ER._base_tier([seed], rxn)
    @test isempty(failures)
    @test length(base) == 1
    @test Set(base) == Set([kept])
end

@testset "_base_tier keeps every flip child the guard admits" begin
    # Every step is rapid equilibrium, so the catalytic cycle E, E(S), E(P), F turns over
    # at rapid equilibrium: the seed is degenerate and expanded. The six flips of two steps
    # of the ring cut it, and each passes every face of the guard: the conformations flip
    # at E, so F needs no rapid-equilibrium link to it. No failure is recorded.
    ER = EnzymeRates
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    seed = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) ⇌ E(P)       :: EqualAI
            E(P) ⇌ F + P      :: EqualAI
            F ⇌ E             :: EqualAI
        end
    end
    # Untyped: inside `T[...]` the parser reads `begin` as an index keyword.
    kept = [
        ER.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)   :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ F + P      :: EqualAI
                F ⇌ E             :: EqualAI
            end
        end),
        ER.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)   :: EqualAI
                E(S) ⇌ E(P)       :: EqualAI
                E(P) <--> F + P   :: EqualAI
                F ⇌ E             :: EqualAI
            end
        end),
        ER.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) <--> F + P   :: EqualAI
                F ⇌ E             :: EqualAI
            end
        end),
        ER.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)   :: EqualAI
                E(S) ⇌ E(P)       :: EqualAI
                E(P) ⇌ F + P      :: EqualAI
                F <--> E          :: EqualAI
            end
        end),
        ER.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ F + P      :: EqualAI
                F <--> E          :: EqualAI
            end
        end),
        ER.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) ⇌ E(P)       :: EqualAI
                E(P) <--> F + P   :: EqualAI
                F <--> E          :: EqualAI
            end
        end)]
    @test ER._degenerate(seed)
    @test length(ER._expand_re_to_ss(seed)) == 6
    @test all(ER._is_well_posed, ER._expand_re_to_ss(seed))
    base, failures = ER._base_tier([seed], rxn)
    @test isempty(failures)
    @test length(base) == 6
    @test Set(base) == Set(kept)
end

@testset "_degenerate reads the active conformation's law" begin
    # An activator X that binds only the inactive conformation opens a route whose
    # chemistry only the inactive law holds. The active law, three steady-state steps,
    # has a maximal rate both ways; read over every step, the route's product release
    # turns rapid equilibrium and closes a turnover cycle with the X complex's substrate
    # binding, so the whole graph has no reverse maximal rate. The guard admits the
    # mechanism.
    route = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> E(S)         :: EqualAI
            E(S) <--> E(P)          :: EqualAI
            E(P) <--> E + P         :: EqualAI
            E + X ⇌ E(X)            :: OnlyI
            E(X) + S ⇌ E(S, X)      :: EqualAI
            E(S, X) <--> E(X) + P   :: EqualAI
        end
    end)
    @test !EnzymeRates._degenerate(route)
    @test EnzymeRates._degenerate(EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> E(S)
            E(S) <--> E(P)
            E(P) <--> E + P
            E + X ⇌ E(X)
            E(X) + S ⇌ E(S, X)
            E(S, X) <--> E(X) + P
        end
    end))
    # The maximal-rate test turns rapid equilibrium every step that takes up or gives
    # off a substrate, so masking only the flags would bring the excluded steps back:
    # here the X route's two Theorell–Chance steps, each taking up S and giving off P,
    # would close a turnover cycle. Read through the edge mask, the active law keeps
    # only its three steady-state steps, which have a maximal rate both ways.
    tc_route = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> E(S)           :: EqualAI
            E(S) <--> E(P)            :: EqualAI
            E(P) <--> E + P           :: EqualAI
            E + X ⇌ E(X)              :: OnlyI
            E(X) + S <--> F(X) + P    :: EqualAI
            F(X) + S <--> E(X) + P    :: EqualAI
        end
    end
    cm, on = tc_route.catalytic, EnzymeRates._active_mask(tc_route)
    re = EnzymeRates._flags(cm)
    @test !EnzymeRates._has_vmax(cm, re .& on, cm.substrates)
    @test EnzymeRates._has_vmax(cm, re, cm.substrates, on)
    @test !EnzymeRates._degenerate(tc_route)
    # Read over every step, the X route's merged complex E(P, X) sits in rapid
    # equilibrium with both sides (condition C), while the maximal rate exists both
    # ways: the steady-state F(X) <--> E(X) keeps the route from closing a
    # rapid-equilibrium turnover cycle. Read through the edge mask, the route lies
    # outside the active law.
    merged_route = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> E(S)           :: EqualAI
            E(S) <--> E(P)            :: EqualAI
            E(P) <--> E + P           :: EqualAI
            E + X ⇌ E(X)              :: OnlyI
            E(X) + S ⇌ E(P, X)        :: EqualAI
            E(P, X) ⇌ F(X) + P        :: EqualAI
            F(X) <--> E(X)            :: EqualAI
        end
    end
    cm, on = merged_route.catalytic, EnzymeRates._active_mask(merged_route)
    re = EnzymeRates._flags(cm)
    subs = EnzymeRates._names_mask(cm, cm.substrates)
    @test EnzymeRates._form_text.(Ref(cm), EnzymeRates._merged_complexes(cm, subs)) ==
          ["E(P, X)"]
    @test isempty(EnzymeRates._merged_complexes(cm, subs, on))
    @test EnzymeRates._chemistry_equilibrates_both_sides(cm, re)
    @test !EnzymeRates._chemistry_equilibrates_both_sides(cm, re, on)
    @test EnzymeRates._has_vmax(cm, re, cm.substrates) &&
          EnzymeRates._has_vmax(cm, re, cm.products)
    @test !EnzymeRates._degenerate(merged_route)
    # The same route through a rapid-equilibrium Theorell–Chance step. Read over every
    # step, that step takes up S and gives off P in rapid equilibrium (condition C), and
    # the maximal rate again exists both ways; read through the edge mask, the step lies
    # outside the active law.
    re_tc_route = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: X
        steps: begin
            E + S <--> E(S)           :: EqualAI
            E(S) <--> E(P)            :: EqualAI
            E(P) <--> E + P           :: EqualAI
            E + X ⇌ E(X)              :: OnlyI
            E(X) + S ⇌ F(X) + P       :: EqualAI
            F(X) <--> E(X)            :: EqualAI
        end
    end
    cm, on = re_tc_route.catalytic, EnzymeRates._active_mask(re_tc_route)
    re = EnzymeRates._flags(cm)
    @test EnzymeRates._chemistry_equilibrates_both_sides(cm, re)
    @test !EnzymeRates._chemistry_equilibrates_both_sides(cm, re, on)
    @test EnzymeRates._has_vmax(cm, re, cm.substrates) &&
          EnzymeRates._has_vmax(cm, re, cm.products)
    @test !EnzymeRates._degenerate(re_tc_route)
end

@testset "_cv_model_selection: flatten reproduces serial LOOCV per candidate" begin
    # Three candidates with distinct keys, each the source of its own fit. A deterministic
    # stub optimizer gives identical fits whether folds run serially or across the flattened
    # (candidate, fold) grid, so each candidate's fold columns must hold its own serial fold
    # losses; a transposed grid would swap them. Group labels containing `=`, `,` and spaces
    # must come out as valid fold-column names that survive the CSV write.
    m1 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products:   P, Q
        steps: begin
            (E + A ⇌ E(A), E(Q) + A ⇌ E(A, Q))
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
        end
    end)
    m3 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products:   P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q ⇌ E(Q)
        end
    end)
    # m3 with a steady-state Q release: one more fitted parameter.
    m6 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products:   P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + P ⇌ E(P, Q)
            E + Q <--> E(Q)
        end
    end)
    data = (group = ["a=b", "a=b", "c,d", "c,d", "x y", "x y"],
            Rate = [0.5, 0.8, 1.0, 1.1, 0.9, 1.2],
            A = [1.0, 2.0, 1.0, 2.0, 1.5, 2.5], B = [0.5, 0.5, 1.0, 1.0, 0.7, 0.7],
            P = [0.1, 0.2, 0.1, 0.2, 0.15, 0.25], Q = [0.3, 0.3, 0.4, 0.4, 0.35, 0.35])
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    prob = IdentifyRateEquationProblem(rxn, data; Keq=2.0)
    # The pool as `_beam_search` hands it over: each candidate with its row and key, and
    # the memo holding each key's fit with the candidate as its source.
    memo = Dict{String, EnzymeRates._Fit}()
    function _testhelper_pool_entry(m, loss)
        em = EnzymeRates.compile_mechanism(m)
        fkeys = EnzymeRates.fitted_params(em)
        text = sprint(show, m)
        id = EnzymeRates._mechanism_id(text)
        key = EnzymeRates._canonical(EnzymeRates._law(m), id).key
        memo[key] = EnzymeRates._Fit(ones(length(fkeys)), loss, :Success, false, m, id)
        row = (n_params=length(fkeys), loss=loss, mechanism=text, mechanism_id=id,
               rate_equation="v", retcode="Success", error=missing,
               params=NamedTuple{fkeys}(ntuple(_ -> 1.0, length(fkeys))), key,
               fit_inherited=false, fit_source_id=id, L_phantom=false)
        EnzymeRates.BatchEntry(m, length(fkeys), loss, :Success, key, key, row)
    end
    # The pool goes in against (n_params, loss) order and must come back sorted by it:
    # the 6-parameter candidate has the lowest loss but sorts last.
    pool = [_testhelper_pool_entry(m6, 0.1), _testhelper_pool_entry(m1, 0.5),
            _testhelper_pool_entry(m3, 0.2)]
    save_dir = mktempdir()
    stub() = _testhelper_CountingStubOpt(; uval=log(5.0))
    res = EnzymeRates._cv_model_selection(pool, memo, prob;
        optimizer=stub(), se_threshold=1.0, save_dir, show_progress=false,
        n_restarts=1, maxtime=1.0)
    @test nrow(res.cv_results) == 3
    @test [e.n_params for e in pool] == [6, 5, 5]
    @test res.cv_results.n_params == [5, 5, 6]
    @test res.cv_results.loss == [0.2, 0.5, 0.1]

    groups = unique(prob.data.group)
    folds(r) = [r[Symbol("cv_fold_$g")] for g in groups]
    for r in eachrow(res.cv_results)
        m = only(e.mech for e in pool if sprint(show, e.mech) == r.mechanism)
        serial = [EnzymeRates._cv_fold_loss(EnzymeRates.compile_mechanism(m), prob, g;
            optimizer=stub(), n_restarts=1, maxtime=1.0) for g in groups]
        @test folds(r) == serial
        @test r.cv_score == mean(serial)
        @test r.cv_score_se == std(serial) / sqrt(length(groups))
    end
    @test folds(res.cv_results[1, :]) != folds(res.cv_results[2, :])

    # Every group label is a fold column in the saved table, holding the same scores.
    saved = CSV.read(joinpath(save_dir, "loocv_results.csv"), DataFrame)
    for g in groups
        @test saved[!, "cv_fold_$g"] == res.cv_results[!, Symbol("cv_fold_$g")]
    end
end
