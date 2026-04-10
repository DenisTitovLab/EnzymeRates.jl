# ABOUTME: Tests for identify_rate_equation pipeline.
# ABOUTME: Covers construction, helpers, and full pipeline
# ABOUTME: with exact mechanism recovery.

using DataFrames
using CSV
using Statistics
using OptimizationPyCMA

@testset "identify_rate_equation" begin

    # ── Shared test setup ────────────────────────────
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

    # Plain reaction (no regulators) — matches the
    # generating mechanism for exact recovery testing
    test_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end

    # Reaction with regulator — for construction
    # validation only
    test_rxn_with_reg = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I
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
                    mechanism, concs,
                    params) *
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

    pycma_opt = PyCMAOpt()

    # ── Construction validation ──────────────────────
    @testset "construction" begin
        prob = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val)
        @test prob.reaction === test_rxn
        @test prob.Keq == Keq_val
        @test length(
            unique(prob.data.group)) == 3

        # Missing metabolite column
        @test_throws(
            ErrorException,
            IdentifyRateEquationProblem(
                test_rxn,
                (group = ["G1"], Rate = [1.0]);
                Keq=1.0))

        # Zero rate
        @test_throws(
            ErrorException,
            IdentifyRateEquationProblem(
                test_rxn,
                (group = ["G1", "G2"],
                 Rate = [0.0, 1.0],
                 S = [1.0, 1.0],
                 P = [0.1, 0.1]);
                Keq=1.0))

        # Need >= 2 groups
        @test_throws(
            ErrorException,
            IdentifyRateEquationProblem(
                test_rxn,
                (group = ["G1", "G1"],
                 Rate = [1.0, 2.0],
                 S = [1.0, 2.0],
                 P = [0.1, 0.1]);
                Keq=1.0))

        # Missing regulator column
        @test_throws(
            ErrorException,
            IdentifyRateEquationProblem(
                test_rxn_with_reg,
                (group = ["G1", "G2"],
                 Rate = [1.0, 2.0],
                 S = [1.0, 2.0],
                 P = [0.1, 0.1]);
                Keq=1.0))

        # Regulator column present — should work
        prob_reg = IdentifyRateEquationProblem(
            test_rxn_with_reg,
            (group = ["G1", "G2"],
             Rate = [1.0, 2.0],
             S = [1.0, 2.0],
             P = [0.1, 0.1],
             I = [0.5, 0.5]);
            Keq=1.0)
        @test prob_reg.reaction ===
            test_rxn_with_reg
    end

    # ── Helper unit tests (cheap, no fitting) ────────
    @testset "_build_result_row" begin
        fp = FittingProblem(
            test_mechanism, test_data;
            Keq=Keq_val)
        fit = fit_rate_equation(
            fp, pycma_opt;
            n_restarts=1, maxtime=3.0,
            popsize=200)
        row = EnzymeRates._build_result_row(
            test_mechanism, fit)
        @test row.n_params == length(
            EnzymeRates.fitted_params(
                test_mechanism))
        @test row.loss == fit.loss
        @test row.mechanism_type isa String
        @test row.rate_equation isa String
        @test length(
            row.fitted_param_names) ==
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
        df = EnzymeRates._rows_to_dataframe(
            rows)
        @test nrow(df) == 1
        @test "a" in names(df)
        @test "b" in names(df)

        # Empty rows
        df2 = EnzymeRates._rows_to_dataframe(
            NamedTuple[])
        @test nrow(df2) == 0
    end

    # ── Run pipeline ONCE, test everything ───────────
    prob = IdentifyRateEquationProblem(
        test_rxn, test_data; Keq=Keq_val)
    save_dir = mktempdir()
    results = identify_rate_equation(prob;
        min_beam_width=200,
        beam_fraction=0.1,
        max_param_count=8,
        n_cv_candidates=3,
        save_dir=save_dir,
        pmap_function=map,
        optimizer=pycma_opt,
        popsize=200,
        n_restarts=3, maxtime=10.0)

    @testset "exact mechanism recovery" begin
        @test typeof(results.best) ==
            typeof(test_mechanism)
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

        first_csv = CSV.read(
            joinpath(
                save_dir, csv_files[1]),
            DataFrame)
        @test "n_params" in names(first_csv)
        @test "loss" in names(first_csv)
        @test "mechanism_type" in names(
            first_csv)
        @test nrow(first_csv) > 0
    end

    @testset "save_dir non-empty check" begin
        @test_throws(
            ErrorException,
            identify_rate_equation(prob;
                min_beam_width=200,
                beam_fraction=0.1,
                max_param_count=8,
                n_cv_candidates=3,
                save_dir=save_dir,
                pmap_function=map,
                optimizer=pycma_opt,
                n_restarts=1, maxtime=1.0))
    end

end
