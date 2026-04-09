# ABOUTME: Tests for identify_rate_equation pipeline.
# ABOUTME: Covers IdentifyRateEquationProblem construction and _beam_search.

using DataFrames
using CSV
using Statistics
using OptimizationBBO

@testset "identify_rate_equation" begin

    # ── Shared test setup ─────────────────────────────────
    test_mechanism = @enzyme_mechanism begin
        species: begin
            substrates: S[C]
            products:   P[C]
            enzymes:    E, ES[C]
        end
        steps: begin
            [E, S] <--> [ES]
            [ES] <--> [E, P]
        end
    end

    test_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end

    Keq_val = 2.0
    true_params = (
        k1f = 10.0, k2f = 5.0, k2r = 1.0,
        Keq = Keq_val, E_total = 1.0)

    function make_test_data(
        mechanism, params;
        n_per_group=5, n_groups=3
    )
        groups = String[]
        rates = Float64[]
        S_vals = Float64[]
        P_vals = Float64[]
        for g in 1:n_groups
            for _ in 1:n_per_group
                s = 0.1 + 9.9 * rand()
                p = 0.1 + 9.9 * rand()
                concs = (S = s, P = p)
                r = rate_equation(
                    mechanism, concs, params) *
                    (0.5 + rand())
                push!(groups, "G$g")
                push!(rates, r)
                push!(S_vals, s)
                push!(P_vals, p)
            end
        end
        return (group = groups, Rate = rates,
                S = S_vals, P = P_vals)
    end

    Random.seed!(42)
    test_data = make_test_data(
        test_mechanism, true_params)

    bbo_opt = BBO_adaptive_de_rand_1_bin_radiuslimited()

    # ── IdentifyRateEquationProblem construction ──────────
    @testset "IdentifyRateEquationProblem construction" begin
        prob = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val)
        @test prob.reaction === test_rxn
        @test prob.Keq == Keq_val
        @test length(unique(prob.data.group)) == 3

        # Missing metabolite column
        bad_data = (group = ["G1"], Rate = [1.0])
        @test_throws ErrorException IdentifyRateEquationProblem(
            test_rxn, bad_data; Keq=1.0)

        # Zero rate
        bad_data2 = (
            group = ["G1", "G2"],
            Rate = [0.0, 1.0],
            S = [1.0, 1.0], P = [0.1, 0.1])
        @test_throws ErrorException IdentifyRateEquationProblem(
            test_rxn, bad_data2; Keq=1.0)

        # Need >= 2 groups
        bad_data3 = (
            group = ["G1", "G1"],
            Rate = [1.0, 2.0],
            S = [1.0, 2.0], P = [0.1, 0.1])
        @test_throws ErrorException IdentifyRateEquationProblem(
            test_rxn, bad_data3; Keq=1.0)
    end

    # ── Helper functions ──────────────────────────────────
    @testset "_build_result_row" begin
        fp = FittingProblem(
            test_mechanism, test_data; Keq=Keq_val)
        fit = fit_rate_equation(
            fp, bbo_opt;
            n_restarts=1, maxtime=3.0)
        row = EnzymeRates._build_result_row(
            test_mechanism, fit)
        @test row.n_params == length(
            EnzymeRates.fitted_params(test_mechanism))
        @test row.loss == fit.loss
        @test row.mechanism_type isa String
        @test row.rate_equation isa String
        @test length(row.fitted_param_names) ==
            length(row.fitted_param_values)
    end

    @testset "_rows_to_dataframe" begin
        rows = [(
            n_params = 3,
            loss = 0.5,
            mechanism_type = "test",
            rate_equation = "v = ...",
            fitted_param_names = (:a, :b),
            fitted_param_values = (1.0, 2.0),
        )]
        df = EnzymeRates._rows_to_dataframe(rows)
        @test nrow(df) == 1
        @test "a" in names(df) || :a in names(df)
        @test "b" in names(df) || :b in names(df)

        # Empty rows
        df2 = EnzymeRates._rows_to_dataframe(
            NamedTuple[])
        @test nrow(df2) == 0
    end

    # ── _beam_search ──────────────────────────────────────
    @testset "_beam_search" begin
        prob = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val)

        specs, df = EnzymeRates._beam_search(prob;
            min_beam_width=200,
            beam_fraction=0.1,
            max_param_count=8,
            save_dir=nothing,
            pmap_function=map,
            optimizer=bbo_opt,
            n_restarts=2, maxtime=5.0)

        @test length(specs) > 0
        @test nrow(df) > 0
        @test nrow(df) == length(specs)

        @test "n_params" in names(df)
        @test "loss" in names(df)
        @test "mechanism_type" in names(df)
        @test "rate_equation" in names(df)

        @test all(isfinite, df.loss)
        @test all(>=(0), df.loss)

        @test length(unique(df.n_params)) >= 1
    end

    # ── _loocv ────────────────────────────────────────────
    @testset "_loocv" begin
        prob = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val)

        cv_score = EnzymeRates._loocv(
            test_mechanism, prob;
            optimizer=bbo_opt,
            n_restarts=2, maxtime=5.0)

        @test isfinite(cv_score)
        @test cv_score >= 0.0
    end

    # ── _cv_model_selection ───────────────────────────────
    @testset "_cv_model_selection" begin
        prob = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val)

        specs, df = EnzymeRates._beam_search(prob;
            min_beam_width=200,
            beam_fraction=0.1,
            max_param_count=8,
            save_dir=nothing,
            pmap_function=map,
            optimizer=bbo_opt,
            n_restarts=2, maxtime=5.0)

        results = EnzymeRates._cv_model_selection(
            specs, df, prob;
            n_cv_candidates=3, pmap_function=map,
            optimizer=bbo_opt,
            n_restarts=2, maxtime=5.0)

        @test results isa IdentifyRateEquationResults
        @test results.best isa
            EnzymeRates.AbstractEnzymeMechanism
        @test nrow(results.cv_results) > 0
        @test "cv_score" in names(results.cv_results)
        @test all(
            isfinite, results.cv_results.cv_score)
    end

end
