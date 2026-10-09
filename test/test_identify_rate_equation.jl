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

    # Build the constrained allosteric mechanism: K-type
    # allosteric with S, P only in R-state (`:OnlyA` group
    # tags) and R only in T-state (`:OnlyI` ligand tag).
    _base = first(EnzymeRates.init_mechanisms(test_rxn))
    _cat_allo_states = Symbol[]
    for g in eachindex(EnzymeRates.steps(_base))
        rep = first(EnzymeRates.steps(_base)[g])
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

        # scale_k_to_kcat field: default 1.0, settable to nothing, validated.
        @test prob.scale_k_to_kcat == 1.0
        prob_abs = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, scale_k_to_kcat=nothing)
        @test prob_abs.scale_k_to_kcat === nothing
        @test_throws ErrorException IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, scale_k_to_kcat=0.0)
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
        # ≥ 0, or a Keq or scale_k_to_kcat that is not positive and finite: each error
        # names the column (with the row), Keq or scale_k_to_kcat.
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
            @test_throws(
                ErrorException("scale_k_to_kcat must be positive and finite (or " *
                               "nothing); got $shown"),
                IdentifyRateEquationProblem(test_rxn, two_groups; Keq=1.0,
                                            scale_k_to_kcat=bad))
        end

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
        # candidate; the problem refuses it up front.
        L_named_rxn = @enzyme_reaction begin
            substrates: S[C]
            products: L[C]
        end
        @test_throws(
            ErrorException("Metabolite L has the name of the conformational constant L " *
                           "of the allosteric mechanisms the search fits; rename the " *
                           "metabolite"),
            IdentifyRateEquationProblem(L_named_rxn,
                (group = ["G1", "G2"], Rate = [1.0, 2.0], S = [1.0, 1.0], L = [0.1, 0.1]);
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
            params = (a = 1.0, b = 2.0),
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

        # A failure row (every fit field missing) beside a fitted row.
        rows_with_failure = [
            (n_params = 3, loss = 0.5, mechanism_type = "M",
             rate_equation = "v = ...", retcode = "Success", error = missing,
             params = (a = 1.0,), eq_hash = "0123456789abcdef", fit_inherited = false),
            (n_params = missing, loss = missing, mechanism_type = "M",
             rate_equation = missing, retcode = missing,
             error = "StackOverflowError: ", params = (;), eq_hash = missing,
             fit_inherited = missing),
        ]
        df_with_failure = EnzymeRates._rows_to_dataframe(rows_with_failure)
        @test nrow(df_with_failure) == 2
        @test ismissing(df_with_failure.loss[2])
        @test df_with_failure.error[2] == "StackOverflowError: "
        @test ismissing(df_with_failure.retcode[2])
        @test ismissing(df_with_failure.eq_hash[2])
        @test df_with_failure.fit_inherited[1] == false
        @test ismissing(df_with_failure.fit_inherited[2])
        @test "a" in names(df_with_failure)    # param column still built from row 1
        @test ismissing(df_with_failure.a[2])  # failure row contributes no param value
    end

    @testset "failure row preserves round-trippable mechanism" begin
        m = first(EnzymeRates.init_mechanisms(_testhelper_uni_rxn))
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
        @test "eq_hash" in names(
            results.cv_results)
        # LOOCV candidate dedup invariant: within each n_params
        # bucket, each eq_hash should appear at most once (the cv
        # pool `_offer_cv!` keeps one slot per eq_hash). This
        # catches a regression where duplicates would enter LOOCV
        # and waste compute / bias the per-bucket "best".
        for gdf in groupby(
                results.cv_results, :n_params)
            @test allunique(gdf.eq_hash)
        end

        # cv_score_se = standard error of the row's per-fold scores.
        @test "cv_score_se" in names(results.cv_results)
        fold_cols = [Symbol("cv_fold_$g") for g in unique(prob.data.group)]
        for row in eachrow(results.cv_results)
            folds = [row[c] for c in fold_cols]
            @test row.cv_score ≈ mean(folds)
            @test row.cv_score_se ≈ std(folds) / sqrt(length(folds))
        end

        @test best_row.mechanism_type == string(typeof(results.best))

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

    @testset "CSV output (new schema)" begin
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
            # Diagnostic parent columns: every saved CSV carries them.
            @test "parent_n_params" in names(df_file)
            @test "parent_mechanism_type" in names(df_file)
        end
        # Base tier is parentless.
        @test all(ismissing, init_df.parent_n_params)
        @test all(ismissing, init_df.parent_mechanism_type)
        # Iteration 1's parents are exactly the base-tier mechanisms; each child
        # records the round-trippable type of the parent it expanded from, with
        # the parent's own param count.
        iter1 = CSV.read(joinpath(save_dir,
                                  "equation_search_iteration_1.csv"), DataFrame)
        parent_np = Dict(string(t) => n
                         for (t, n) in zip(init_df.mechanism_type, init_df.n_params))
        @test any(!ismissing, iter1.parent_mechanism_type)   # children have parents
        for r in eachrow(iter1)
            ismissing(r.parent_mechanism_type) && continue
            @test string(r.parent_mechanism_type) in keys(parent_np)
            @test r.parent_n_params == parent_np[string(r.parent_mechanism_type)]
        end
        # The stored parent type round-trips back to a mechanism.
        pt = String(first(skipmissing(iter1.parent_mechanism_type)))
        @test EnzymeRates.Mechanism(
            Core.eval(EnzymeRates, Meta.parse(pt))()) isa EnzymeRates.Mechanism
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

    @testset "_cv_fold_loss: finite per-fold scores, loud on fit failure" begin
        m = EnzymeRates.EnzymeMechanism(
            first(EnzymeRates.init_mechanisms(_testhelper_uni_rxn)))

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
    em = EnzymeRates.compile_mechanism(
        first(EnzymeRates.init_mechanisms(_testhelper_uni_rxn)))
    for data_form in (NamedTuple, DataFrame)
        prob = _testhelper_uni_prob(data_form)
        full = FittingProblem(em, prob.data; Keq=prob.Keq,
                              scale_k_to_kcat=prob.scale_k_to_kcat)
        fold = EnzymeRates._fold_problem(em, prob, prob.data.group .== prob.data.group[3])
        @test typeof(fold) === typeof(full)
        @test fold.data.S == [3.0, 4.0]
        @test fold.data.P == [0.3, 0.4]
        @test fold.data.Rate == [1.0, 1.1]
        @test fold.Keq === 10.0 && fold.scale_k_to_kcat === 1.0
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
        full = FittingProblem(em, prob.data; Keq=prob.Keq,
                              scale_k_to_kcat=prob.scale_k_to_kcat)
        fit_rate_equation(full, optimizer; n_restarts=1, maxiters=10)
        n = _testhelper_enzymerates_method_instances()
        EnzymeRates._cv_fold_loss(em, prob, "G1"; optimizer, n_restarts=1, maxiters=10)
        push!(added, _testhelper_enzymerates_method_instances() - n)
    end
    # The first mechanism's fold compiles the fold path for every mechanism.
    @test added[2] == 0
end

# The CSV mechanism_type column holds a mechanism's type as string(typeof(em)) prints it
# where EnzymeRates is loaded with `using`, and the string parses back to that type.
@testset "_mechanism_type_string prints the mechanism type" begin
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
        allosteric_regulators: A::OnlyA, I::OnlyI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: EqualAI
            E(S) <--> E(P)   :: OnlyA
            E(P) ⇌ E + P     :: EqualAI
        end
        regulatory_site(multiplicity = 4): begin
            ligands: A, I
        end
    end
    for em in (plain, ping_pong, allosteric)
        s = EnzymeRates._mechanism_type_string(em)
        @test s == string(typeof(em))
        @test Core.eval(EnzymeRates, Meta.parse(s)) === typeof(em)
    end
end

@testset "csv writer" begin
    rows = [(
        n_params = 5, loss = 1.0, mechanism_type = "M",
        rate_equation = "v = 1", retcode = "Success", error = missing,
        params = (K_a = 2.0,), eq_hash = "abc",
        fit_inherited = false,
    )]
    mktempdir() do tmp
        EnzymeRates._write_rows_csv(tmp, "equation_search_iteration_3.csv", rows)
        @test isfile(joinpath(tmp, "equation_search_iteration_3.csv"))
        df = CSV.read(joinpath(tmp, "equation_search_iteration_3.csv"), DataFrame)
        @test df.n_params == [5]
        @test df.K_a == [2.0]
        @test "eq_hash" in names(df)
        # dir-creation branch: save_dir does not exist yet
        subdir = joinpath(tmp, "made")
        EnzymeRates._write_rows_csv(subdir, "initial_mechanisms.csv", rows)
        @test isfile(joinpath(subdir, "initial_mechanisms.csv"))
    end
end

@testset "_select_count!: thresholds, floor, best loss, parsimony" begin
    # One call at count 5 with a fresh floor budget: `best` is the count's best loss
    # over the whole search, and `others` holds the best losses of other counts.
    _testhelper_selected(losses, best; rel, add = 0.0, width = 1,
                         others = Dict{Int,Float64}(), parsimony = 1.0) =
        EnzymeRates._select_count!(Dict{Int,Int}(), merge(Dict(5 => best), others), 5,
            losses; loss_rel_threshold = rel, loss_abs_threshold = add,
            loss_parsimony_threshold = parsimony, min_beam_width = width)

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
    # are all failures (non-missing `error`, missing `eq_hash`).
    @test isfile(joinpath(tmp, "initial_mechanisms.csv"))
    fail_df = CSV.read(joinpath(tmp, "initial_mechanisms.csv"), DataFrame)
    @test nrow(fail_df) >= 1
    @test all(.!ismissing.(fail_df.error))
    @test all(ismissing.(fail_df.eq_hash))
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

@testset "rate-eq dedup-key partition stability" begin
    # bi_bi exercises the dedup key's edge cases: substituted-into-v ties across
    # multiple kinetic groups. uni_uni has one init mechanism, hence one class, so
    # it adds nothing to a partition test.
    # ter-ter intentionally omitted — `rate_equation_string` derivation is
    # extremely slow for mechanisms with >~30 enzyme forms (CLAUDE.md
    # "Known Issues"), and the dedup key renders that string per
    # candidate. The bi-bi enumeration already covers every structural
    # symmetry the dedup key collapses.
    reaction = @enzyme_reaction(begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end)

    # Expected partition size = the number of DISTINCT rate equations the
    # init-level enumeration produces. The 239 bi_bi init mechanisms (55 seeds
    # and their 184 merged and Theorell–Chance variants) are all structurally
    # distinct AND each yields a distinct `rate_equation_string`, so the
    # comment-stripped string key produces exactly 239 classes (zero over- and
    # zero under-collapse): clean topologies have distinct enzyme-form sets,
    # hence distinct rate equations. The 8 ordered/random pairs and the 2
    # Theorell–Chance pairs among the variants share a family, but each pair's
    # members are structurally distinct, so they render distinct equations.
    # If this count changes in a future commit, the dedup key's
    # equivalence classes (or the enumeration) have shifted — investigate.
    expected_n_classes = 239

    # init_mechanisms only — skip expand_mechanisms. The init level
    # already produces multiple structurally-equivalent variants
    # (mirror-step orderings, kinetic-group renumberings) that
    # exercise the dedup key's collapse rules. expand_mechanisms
    # adds variants at higher param counts whose dedup-key
    # behavior is the same modulo size, at exponential compile cost.
    buckets = Dict{UInt64, Vector{Int}}()
    for (i, m) in enumerate(EnzymeRates.init_mechanisms(reaction))
        em = EnzymeRates.compile_mechanism(m)
        h = EnzymeRates._rate_eq_dedup_key(rate_equation_string(em))
        push!(get!(buckets, h, Int[]), i)
        # Determinism on a fixed sample: same input, same key across invocations.
        i % 25 == 1 && @test EnzymeRates._rate_eq_dedup_key(rate_equation_string(em)) === h
    end

    @test length(buckets) == expected_n_classes
end

@testset "_select_count! cumulative per-count floor" begin
    expanded = Dict{Int,Int}()
    best = Dict(5 => 1.0)
    kw = (loss_rel_threshold=1.0, loss_abs_threshold=0.0,
          loss_parsimony_threshold=1.0, min_beam_width=3)
    # Sweep 1 at count 5: rel cutoff admits only the best (loss 1.0); the
    # floor budget (3) tops it up to the top 3 by loss. expanded[5] -> 3.
    sel1 = EnzymeRates._select_count!(expanded, best, 5, [1.0, 2.0, 3.0, 4.0, 5.0]; kw...)
    @test sort(sel1) == [1, 2, 3]
    @test expanded[5] == 3

    # Sweep 2 at count 5: budget spent (3 of 3). New mechanisms all above the
    # cutoff -> the floor admits NONE (unlike the old per-sweep floor, which
    # would grant a fresh 3). expanded[5] stays 3.
    sel2 = EnzymeRates._select_count!(expanded, best, 5, [10.0, 11.0, 12.0]; kw...)
    @test isempty(sel2)
    @test expanded[5] == 3

    # A cutoff-passer is still admitted after the floor is spent.
    sel3 = EnzymeRates._select_count!(expanded, best, 5, [1.0, 20.0]; kw...)
    @test sel3 == [1]
    @test expanded[5] == 4
end

@testset "§1 parsimony cutoff = threshold * min over all counts < c" begin
    # No floor and a loose relative cutoff (10 × the count's best), so the parsimony
    # cutoff alone decides: 1.01 × the best loss over the counts below c.
    _testhelper_selected(best_loss_by_count, c, losses) = EnzymeRates._select_count!(
        Dict{Int,Int}(), best_loss_by_count, c, losses; loss_rel_threshold=10.0,
        loss_abs_threshold=0.0, loss_parsimony_threshold=1.01, min_beam_width=0)
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
    mech = first(EnzymeRates.init_mechanisms(_testhelper_uni_rxn))
    row = (n_params=3, loss=0.5, mechanism_type="M", rate_equation="v",
           retcode="Success", error=missing, params=(K = 1.0,), eq_hash="abc",
           fit_inherited=false)
    e_succ = EnzymeRates.BatchEntry(mech, 3, 0.5, :Success, hash(:a), row)
    e_mt   = EnzymeRates.BatchEntry(mech, 3, 0.9, :MaxTime, hash(:b), row)
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
    @test all(e -> occursin(r"^[0-9a-f]{16}$", e.row.eq_hash), entries)

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
        n_restarts=1, maxtime=1.0, memo=Dict{UInt64,NamedTuple}(), seen)
    @test ss1 == 0 && length(e1) == 1

    # Same structure again in a later batch → seen-skipped, no new entry.
    seen_log = String[]
    e2, f2, ps2, cs2, ss2 = EnzymeRates._process_batch([m], prob;
        optimizer=CMAEvolutionStrategyOpt(), max_param_count=20, eq_complexity_filter=337,
        n_restarts=1, maxtime=1.0, memo=Dict{UInt64,NamedTuple}(), seen,
        log = msg -> push!(seen_log, msg))
    @test ss2 == 1 && isempty(e2) && isempty(f2)
    @test seen_log == ["0 new fits + 0 inherited + 1 skipped (already fit) + " *
                       "0 skipped (>20 params) + 0 skipped (>337 complexity)"]

    # A structure repeated within one batch is fit once; the repeat is seen-skipped.
    e3, f3, ps3, cs3, ss3 = EnzymeRates._process_batch([m, m], prob;
        optimizer=CMAEvolutionStrategyOpt(), max_param_count=20,
        n_restarts=1, maxtime=1.0, memo=Dict{UInt64,NamedTuple}(), seen=Set{UInt64}())
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

@testset "_ingest! and cv pool" begin
    _testhelper_mk(n, loss, h) = EnzymeRates.BatchEntry(
        first(EnzymeRates.init_mechanisms(_testhelper_uni_rxn)),
        n, loss, :Success, hash(h),
        (n_params=n, loss=loss, mechanism_type="M",
         rate_equation="v", retcode="Success", error=missing,
         params=(K = 1.0,), eq_hash=string(hash(h),base=16,pad=16),
         fit_inherited=false))
    frontier = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    cv_pool  = Dict{Int,Vector{EnzymeRates.BatchEntry}}()
    best     = Dict{Int,Float64}()
    # two distinct equations + one duplicate-eq with worse loss, n_cv=2
    improved = EnzymeRates._ingest!(frontier, cv_pool, best,
        [_testhelper_mk(5,2.0,:a), _testhelper_mk(5,1.0,:b), _testhelper_mk(5,3.0,:a)];
        n_cv_candidates=2)
    @test improved == Set([5])                 # count 5 first appeared
    @test length(frontier[5]) == 3            # frontier keeps ALL
    @test best[5] == 1.0                       # running min
    @test length(cv_pool[5]) == 2              # bounded, distinct eq_hash
    # the kept :a entry is the lower-loss one (2.0, not 3.0);
    # BatchEntry.eq_hash is a UInt64, so compare against hash(:a), not hex:
    a = only(filter(e -> e.eq_hash == hash(:a), cv_pool[5]))
    @test a.loss == 2.0
    # A later batch reports only the counts whose best strictly dropped or first
    # appeared: count 5 ties its best (1.0), count 6 first appears, and count 7's
    # entry is worse than its best (0.5).
    best[7] = 0.5
    @test EnzymeRates._ingest!(frontier, cv_pool, best,
        [_testhelper_mk(5,1.0,:c), _testhelper_mk(6,4.0,:d), _testhelper_mk(7,0.6,:e)];
        n_cv_candidates=2) == Set([6])
    @test EnzymeRates._ingest!(frontier, cv_pool, best,
        [_testhelper_mk(5,0.9,:f)]; n_cv_candidates=2) == Set([5])
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
        [_testhelper_mk(5,Inf,:g), _testhelper_mk(5,1.0,:h), _testhelper_mk(6,NaN,:i),
         _testhelper_mk(7,Inf,:j)];
        n_cv_candidates=2) == Set([5])
    @test Set(keys(frontier_nf)) == Set(keys(cv_pool_nf)) == Set([5])
    @test [e.eq_hash for e in frontier_nf[5]] == [hash(:h)]
    @test [e.eq_hash for e in cv_pool_nf[5]] == [hash(:h)]
    @test best_nf == Dict(5 => 1.0)

    # `_offer_cv!` keeps at most one entry per eq_hash: a repeat hash updates its own
    # slot to the lower loss, never consuming a second.
    pool = EnzymeRates.BatchEntry[]
    for (loss, h) in [(1.0, :a), (0.5, :a), (2.0, :b)]
        EnzymeRates._offer_cv!(pool, _testhelper_mk(5, loss, h), 5)
    end
    @test allunique([e.eq_hash for e in pool])
    @test length(pool) == 2
    @test only(filter(e -> e.eq_hash == hash(:a), pool)).loss == 0.5
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
        allosteric_regulators: ATP(1)
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

# ── §2 fit-dedup by eq_hash ──────────────────────────────────────────────────
# A stub optimizer that counts `solve` invocations and returns a canned
# log-space optimum (`uval` for every coordinate), so a batch's fits can be
# counted exactly and the raw→rescale path exercised deterministically.
mutable struct _testhelper_CountingStubOpt
    count::Int
    uval::Float64
    throwit::Bool
end
_testhelper_CountingStubOpt(; uval=0.0, throwit=false) =
    _testhelper_CountingStubOpt(0, uval, throwit)
Optimization.allowsbounds(::_testhelper_CountingStubOpt) = true
function Optimization.SciMLBase.__solve(
        prob::Optimization.OptimizationProblem, opt::_testhelper_CountingStubOpt; kwargs...)
    opt.count += 1
    opt.throwit && error("stub solver forced failure")
    u = fill(opt.uval, length(prob.u0))
    cache = DefaultOptimizationCache(prob.f, prob.p)
    build_solution(cache, opt, u, prob.f(u, prob.p); retcode = ReturnCode.Success)
end

# Two structurally-distinct bi-bi mechanisms that render the SAME reduced rate
# equation: identical 7-edge topology, differing only in kinetic-group partitioning.
# Written out directly so the tests are self-contained and cheap (no full child
# enumeration, which would compile ~800 distinct rate equations). The tests re-verify
# the collision at run time.
function _testhelper_dedup_twins()
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
    m2 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products:   P, Q
        steps: begin
            E + A ⇌ E(A)
            (E + Q ⇌ E(Q), E(A) + Q ⇌ E(A, Q))
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(Q) + A ⇌ E(A, Q)
            E(Q) + P ⇌ E(P, Q)
        end
    end)
    return m1, m2
end

@testset "fit-dedup by eq_hash in _process_batch" begin
    m1, m2 = _testhelper_dedup_twins()
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
    opt = _testhelper_CountingStubOpt(; uval = log(5.0))
    # Each logged line is stored with the solve count at the moment it was logged.
    batch_log = Tuple{String,Int}[]
    entries, failures = EnzymeRates._process_batch(pair, prob;
        optimizer=opt, max_param_count=20, eq_complexity_filter=337, n_restarts=1,
        maxtime=1.0, memo, log = msg -> push!(batch_log, (msg, opt.count)))

    # The shared equation is fit exactly ONCE (n_restarts=1 → one solve).
    @test opt.count == 1
    # The pre-fit summary counts one new fit and one inherited row, and is logged
    # before the fit runs.
    @test batch_log == [("1 new fits + 1 inherited + 0 skipped (already fit) + " *
                         "0 skipped (>20 params) + 0 skipped (>337 complexity)", 0)]
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
        @test keys(e.row.params) == fkeys
        @test e.row.params == fit.params
    end
    @test entries[1].row.params == entries[2].row.params
    # scale_k_to_kcat=1.0 anchored kcat: the copied params are the rescaled fit,
    # not the raw 5.0 the stub optimizer returned.
    @test !all(v -> v ≈ 5.0, entries[1].row.params)

    # Cross-batch memo hit: a later batch with the same eq_hash refits NOTHING.
    single = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m1]
    reuse_log = String[]
    reused, _ = EnzymeRates._process_batch(single, prob;
        optimizer=opt, max_param_count=20, eq_complexity_filter=337, n_restarts=1,
        maxtime=1.0, memo, log = msg -> push!(reuse_log, msg))
    @test opt.count == 1                                # no new solve
    @test [e.row.fit_inherited for e in reused] == [true]
    @test reuse_log == ["0 new fits + 1 inherited + 0 skipped (already fit) + " *
                        "0 skipped (>20 params) + 0 skipped (>337 complexity)"]

    # A representative whose fit throws fails ALL its duplicates
    # (all-or-nothing per equation).
    opt_bad = _testhelper_CountingStubOpt(; throwit=true)
    bad_entries, bad_failures = EnzymeRates._process_batch(pair, prob;
        optimizer=opt_bad, max_param_count=20, n_restarts=1, maxtime=1.0,
        memo = Dict{UInt64, NamedTuple}())
    @test isempty(bad_entries)
    @test length(bad_failures) == 2
    @test all(f -> f isa EnzymeRates.FitFailure, bad_failures)
    @test [f.mech for f in bad_failures] == [m1, m2]   # the mechanisms as handed in
end


# LDH renaming-dup pair (same graph, tied kinetic-group split): merged form
# (8 groups) vs split form (9 groups). Different eq_hash, one model: the split
# form's extra binding K is tied straight back by a Wegscheider cycle.
@testset "renaming-dup pair: same independent count, different eq_hash" begin
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
    em1 = EnzymeRates.compile_mechanism(m1)
    em2 = EnzymeRates.compile_mechanism(m2)
    # Precondition: the two forms are distinct mechanisms and render different
    # dedup keys, so `eq_hash` alone never collapses this pair.
    @test m1 != m2
    @test EnzymeRates._rate_eq_dedup_key(rate_equation_string(em1)) !=
          EnzymeRates._rate_eq_dedup_key(rate_equation_string(em2))

    # The pair is one model under the constraint solve: the split form's extra
    # binding K is tied straight back, so both carry the same independent
    # parameters. `eq_hash` does not see through which tied name survives, which
    # is why the split move rejects a candidate by parameter count rather than by
    # equation string.
    @test EnzymeRates._independent_param_count(m1) ==
          EnzymeRates._independent_param_count(m2)
end


# LDH ALLOSTERIC split/merge pair: split form (8 cat groups) vs merged (7).
# Different eq_hash, and the split form carries one parameter the merged form
# cannot express.
@testset "allosteric split/merge pair: the split form carries one more param" begin
    # Split form: 8 catalytic groups.
    am1 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: NADH, Pyruvate
        products:   Lactate, NAD
        catalytic_multiplicity: 4
        catalytic_steps: begin
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
        catalytic_steps: begin
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
    _testhelper_key(m) = EnzymeRates._rate_eq_dedup_key(
        rate_equation_string(EnzymeRates.compile_mechanism(m)))
    @test am1 != am2
    # the two render different equations
    @test _testhelper_key(am1) != _testhelper_key(am2)
    # The split form is not a reparameterization of the merged one: it carries
    # K_ELactateNADH_to_ENADH_Lactate on top of the merged form's parameters, in the
    # independent count and in the fitted set alike. A finite-difference rank of ∂v/∂θ,
    # measured outside this file, agrees (9 against 8).
    @test EnzymeRates._independent_param_count(am1) ==
          EnzymeRates._independent_param_count(am2) + 1
    @test length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(am1))) ==
          length(EnzymeRates.fitted_params(EnzymeRates.compile_mechanism(am2))) + 1
end


@testset "_process_batch failures report the ORIGINAL mechanism" begin
    # A mechanism whose derivation throws — its chemistry step consumes an atom
    # of the never-binding substrate T, so the thermodynamic-cycle check inside
    # compile_mechanism raises "Cycle 1 produces metabolite change not
    # proportional to net reaction". The FitFailure must carry the ORIGINAL `m0`.
    rxn_bad = @enzyme_reaction begin
        substrates: S[C], T[N]
        products:   P[CN]
    end
    e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
    m_bad = EnzymeRates.Mechanism(rxn_bad, [
        [EnzymeRates.Step(e, e_s, [EnzymeRates.Substrate(:S)],
                          EnzymeRates.Metabolite[], true)],
        [EnzymeRates.Step(e_s, e_p, EnzymeRates.Metabolite[],
                          EnzymeRates.Metabolite[], false)],
        [EnzymeRates.Step(e, e_p, [EnzymeRates.Product(:P)],
                          EnzymeRates.Metabolite[], true)],
    ])
    data_bad = (group = ["G1", "G1", "G2", "G2"], Rate = [0.5, 0.8, 1.0, 1.1],
                S = [1.0, 2.0, 1.0, 2.0], T = [0.5, 0.5, 1.0, 1.0],
                P = [0.1, 0.2, 0.1, 0.2])
    prob_bad = IdentifyRateEquationProblem(rxn_bad, data_bad; Keq=2.0)
    e1, f1 = EnzymeRates._process_batch(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m_bad],
        prob_bad; optimizer=_testhelper_CountingStubOpt(), max_param_count=20,
        n_restarts=1, maxtime=1.0, memo=Dict{UInt64, NamedTuple}())
    @test isempty(e1)
    @test length(f1) == 1 && f1[1] isa EnzymeRates.FitFailure
    @test f1[1].mech == m_bad                     # the mechanism as handed in
end

@testset "_expand_parents records an expansion error instead of aborting" begin
    # expand_mechanisms asserts its input conserves atoms; this mechanism's
    # chemistry step does not (T is a declared substrate that never binds, so the
    # step loses an N), and the assertion raises. _expand_parents must catch that
    # and return the parent as a FitFailure (so the beam records it in CSV and
    # continues), not propagate and abort the search.
    rxn_bad = @enzyme_reaction begin
        substrates: S[C], T[N]
        products:   P[CN]
    end
    e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
    m_bad = EnzymeRates.Mechanism(rxn_bad, [
        [EnzymeRates.Step(e, e_s, [EnzymeRates.Substrate(:S)],
                          EnzymeRates.Metabolite[], true)],
        [EnzymeRates.Step(e_s, e_p, EnzymeRates.Metabolite[],
                          EnzymeRates.Metabolite[], false)],
        [EnzymeRates.Step(e, e_p, [EnzymeRates.Product(:P)],
                          EnzymeRates.Metabolite[], true)],
    ])
    # Sanity: expansion of this mechanism genuinely raises.
    @test_throws ErrorException EnzymeRates.expand_mechanisms(
        Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m_bad],
        rxn_bad)
    entry(m) = EnzymeRates.BatchEntry(m, 0, 0.0, :Success, UInt64(0), (mechanism_type="M",))
    kids, parent_of, failures = EnzymeRates._expand_parents([entry(m_bad)], rxn_bad)
    @test isempty(kids) && isempty(parent_of)
    failure = only(failures)
    @test failure isa EnzymeRates.FitFailure
    @test failure.mech == m_bad                    # the ORIGINAL parent
    # The recorded error is the assertion's: the chemistry step ES → EP loses T's N.
    @test occursin("atom-non-conserving step ES → EP", failure.error)
    @test occursin("= Dict(:N => 1)", failure.error)
    # A well-formed parent expands with no failure.
    good = first(EnzymeRates.init_mechanisms(rxn_bad))
    gkids, _, gfail = EnzymeRates._expand_parents([entry(good)], rxn_bad)
    @test isempty(gfail) && !isempty(gkids)
end

@testset "_expand_parents matches each parent's expand_mechanisms" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I
    end
    M = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}
    # The seed's children as parents: their own children overlap, so the
    # first-parent dedup and attribution are exercised.
    parents = EnzymeRates.expand_mechanisms(M[EnzymeRates.init_mechanisms(rxn)...], rxn)
    to_expand = [EnzymeRates.BatchEntry(m, 3, 0.0, :Success, hash(m), (mechanism_type="M",))
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
    attach(em) = ER.Mechanism(rxn, ER.steps(ER.Mechanism(em)))
    seed = attach(@enzyme_mechanism begin
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
    chemistry_flipped = attach(@enzyme_mechanism begin
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
    b_flipped = attach(@enzyme_mechanism begin
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
    q_flipped = attach(@enzyme_mechanism begin
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
    # The chemistry step E(S) ⇌ E(P) loses the N of T, a declared substrate that never
    # binds, so expand_mechanisms' atom-conservation assertion raises. With the P binding
    # steady state the seed has no maximal rate in products (turning that binding rapid
    # equilibrium closes an all-RE turnover cycle), so it is degenerate and expanded.
    rxn_bad = @enzyme_reaction begin
        substrates: S[C], T[N]
        products:   P[CN]
    end
    e   = EnzymeRates.Species(EnzymeRates.Metabolite[], :E)
    e_s = EnzymeRates.Species([EnzymeRates.Substrate(:S)], :E)
    e_p = EnzymeRates.Species([EnzymeRates.Product(:P)], :E)
    m_bad = EnzymeRates.Mechanism(rxn_bad, [
        [EnzymeRates.Step(e, e_s, [EnzymeRates.Substrate(:S)],
                          EnzymeRates.Metabolite[], true)],
        [EnzymeRates.Step(e_s, e_p, EnzymeRates.Metabolite[],
                          EnzymeRates.Metabolite[], true)],
        [EnzymeRates.Step(e, e_p, [EnzymeRates.Product(:P)],
                          EnzymeRates.Metabolite[], false)],
    ])
    @test EnzymeRates._degenerate(m_bad)
    base, failures = EnzymeRates._base_tier([m_bad], rxn_bad)
    @test isempty(base)
    @test length(failures) == 1 && failures[1] isa EnzymeRates.FitFailure
    @test failures[1].mech == m_bad
    # The recorded error is the assertion's: the chemistry step ES → EP loses T's N.
    @test occursin("atom-non-conserving step ES → EP", failures[1].error)
    @test occursin("= Dict(:N => 1)", failures[1].error)
end

@testset "_cv_model_selection: flatten reproduces serial LOOCV per candidate" begin
    # Two candidates with distinct equations. A deterministic stub optimizer gives
    # identical fits whether folds run serially or across the flattened (candidate,
    # fold) grid, so each candidate's fold columns must hold its own serial fold
    # losses; a transposed grid would swap them. Group labels containing `=`, `,`
    # and spaces must come out as valid fold-column names that survive the CSV write.
    m1, _ = _testhelper_dedup_twins()
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
    cands = Union{EnzymeRates.Mechanism, EnzymeRates.AllostericMechanism}[m6, m1, m3]
    data = (group = ["a=b", "a=b", "c,d", "c,d", "x y", "x y"],
            Rate = [0.5, 0.8, 1.0, 1.1, 0.9, 1.2],
            A = [1.0, 2.0, 1.0, 2.0, 1.5, 2.5], B = [0.5, 0.5, 1.0, 1.0, 0.7, 0.7],
            P = [0.1, 0.2, 0.1, 0.2, 0.15, 0.25], Q = [0.3, 0.3, 0.4, 0.4, 0.35, 0.35])
    prob = IdentifyRateEquationProblem(EnzymeRates.reaction(m1), data; Keq=2.0)
    function _testhelper_mkrow(m, loss)
        em = EnzymeRates.compile_mechanism(m)
        fkeys = EnzymeRates.fitted_params(em)
        (n_params=length(fkeys), loss=loss, mechanism_type=string(typeof(em)),
         rate_equation="v", retcode="Success", error=missing,
         params=NamedTuple{fkeys}(ntuple(_ -> 1.0, length(fkeys))),
         eq_hash=string(EnzymeRates._rate_eq_dedup_key(rate_equation_string(em)),
                        base=16, pad=16),
         fit_inherited=false)
    end
    # `cands` and the rows are parallel and go in against (n_params, loss) order; they
    # must come back sorted by it: the 6-parameter row has the lowest loss but sorts last.
    df = EnzymeRates._rows_to_dataframe([_testhelper_mkrow(m6, 0.1),
        _testhelper_mkrow(m1, 0.5), _testhelper_mkrow(m3, 0.2)])
    save_dir = mktempdir()
    stub() = _testhelper_CountingStubOpt(; uval=log(5.0))
    res = EnzymeRates._cv_model_selection(cands, df, prob;
        optimizer=stub(), se_threshold=1.0, save_dir, show_progress=false,
        n_restarts=1, maxtime=1.0)
    @test nrow(res.cv_results) == 3
    @test df.n_params == [6, 5, 5]
    @test res.cv_results.n_params == [5, 5, 6]
    @test res.cv_results.loss == [0.2, 0.5, 0.1]

    groups = unique(prob.data.group)
    folds(r) = [r[Symbol("cv_fold_$g")] for g in groups]
    for r in eachrow(res.cv_results)
        m = only(c for c in cands
                 if string(typeof(EnzymeRates.compile_mechanism(c))) == r.mechanism_type)
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
