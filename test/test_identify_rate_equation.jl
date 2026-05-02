# ABOUTME: Tests for identify_rate_equation pipeline.
# ABOUTME: Covers construction, helpers, and full pipeline
# ABOUTME: with mechanism recovery including allosteric path.

using DataFrames
using CSV
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
        regulators: R
    end

    # Build the constrained allosteric mechanism: K-type
    # allosteric with S, P only in R-state (`:OnlyR` group
    # tags) and R only in T-state (`:OnlyT` ligand tag).
    _init = EnzymeRates.init_mechanisms(test_rxn)
    _base_spec = _init[1]
    # Find S-binding and P-binding kinetic groups by their
    # metabolite (single-step groups in init).
    _g_s = first(s.kinetic_group for s in _base_spec.steps
                 if EnzymeRates.step_metabolite(s) === :S)
    _g_p = first(s.kinetic_group for s in _base_spec.steps
                 if EnzymeRates.step_metabolite(s) === :P)
    _allo_spec =
        EnzymeRates.AllostericMechanismSpec(
            _base_spec,
            1,                # catalytic_n
            [[:R]],           # reg sites
            [1],              # multiplicities
            Dict(_g_s => :OnlyR, _g_p => :OnlyR),
            Dict(:R => :OnlyT),
            8)                # param_count
    test_mechanism =
        EnzymeRates.AllostericEnzymeMechanism(
            _allo_spec)

    Keq_val = 2.0
    # 5 identifiable params + 1 ghost (k3f_T)
    true_params = (
        K1 = 1.0, K2 = 0.5, k3f = 5.0,
        k3f_T = 1.0, K_R_T_reg1 = 2.0, L = 0.1,
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
