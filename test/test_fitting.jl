using Tables

@testset "Fitting" begin

    # ── Helper: build a Uni-Uni mechanism ─────────────────────────────────────
    uni_uni = @enzyme_mechanism begin
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

    # ── Synthetic data generator ──────────────────────────────────────────────
    function make_synthetic_data(
            mechanism, true_params, concs_list;
            articles=fill("A1", length(concs_list)),
            figs=fill("F1", length(concs_list)),
            scale=1.0,
    )
        rates = Float64[]
        for (i, concs) in enumerate(concs_list)
            r = rate_equation(mechanism, concs, true_params) * scale
            push!(rates, r)
        end
        met_names = metabolites(mechanism)
        cols = Dict{Symbol, Vector}()
        cols[:Article] = articles
        cols[:Fig] = figs
        cols[:Rate] = rates
        for mn in met_names
            cols[mn] = [concs[mn] for concs in concs_list]
        end
        return (; (k => cols[k] for k in (:Article, :Fig, :Rate, met_names...))...)
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
        true_params = (k1f = 10.0, k2f = 5.0, k2r = 1.0, Keq = Keq_val, E_total = 1.0)

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
        @test length(fp.fig_point_indexes) == 1  # all same (Article, Fig)
        @test length(fp.fig_point_indexes[1]) == 5
    end

    # ── Test 3: Loss function correctness ─────────────────────────────────────
    @testset "Loss at true params is zero" begin
        Keq_val = 2.0
        true_params = (k1f = 10.0, k2f = 5.0, k2r = 1.0, Keq = Keq_val, E_total = 1.0)

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

    # ── Test 4: Per-figure centering invariance ───────────────────────────────
    @testset "Centering invariance" begin
        Keq_val = 2.0
        true_params = (k1f = 10.0, k2f = 5.0, k2r = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [
            (S = 1.0, P = 0.1),
            (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1),
            (S = 1.0, P = 0.5),
            (S = 2.0, P = 0.5),
        ]

        # Data with scale=1
        data1 = make_synthetic_data(uni_uni, true_params, concs_list;
            articles=fill("A1", 5), figs=fill("F1", 5), scale=1.0)
        fp1 = FittingProblem(uni_uni, data1; Keq=Keq_val)

        # Data with scale=10 (simulates different E_total)
        data2 = make_synthetic_data(uni_uni, true_params, concs_list;
            articles=fill("A1", 5), figs=fill("F1", 5), scale=10.0)
        fp2 = FittingProblem(uni_uni, data2; Keq=Keq_val)

        # For any x, loss should be the same (centering removes the uniform scale)
        np = length(EnzymeRates.fitted_params(uni_uni))
        @test all(1:10) do _
            x = randn(np) .* 2.0
            isapprox(EnzymeRates.loss!(x, fp1), EnzymeRates.loss!(x, fp2); rtol=1e-12)
        end
    end

    # ── Test 5: Multi-figure centering invariance ─────────────────────────────
    @testset "Multi-figure centering invariance" begin
        Keq_val = 2.0
        true_params = (k1f = 10.0, k2f = 5.0, k2r = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [
            (S = 1.0, P = 0.1),
            (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1),
            (S = 1.0, P = 0.5),
            (S = 2.0, P = 0.5),
        ]

        # Two figures, each independently scaled
        data1 = make_synthetic_data(uni_uni, true_params, concs_list;
            articles=["A1","A1","A1","A2","A2"],
            figs=["F1","F1","F1","F1","F1"],
            scale=1.0)
        fp1 = FittingProblem(uni_uni, data1; Keq=Keq_val)

        # Scale fig1 by 5x and fig2 by 100x
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

    # ── Test 6: Sign-mismatch penalty ─────────────────────────────────────────
    @testset "Sign mismatch penalty" begin
        Keq_val = 2.0
        # Create data with positive rates
        data = (
            Article = ["A1", "A1", "A1"],
            Fig = ["F1", "F1", "F1"],
            Rate = [1.0, 2.0, 3.0],
            S = [1.0, 2.0, 3.0],
            P = [0.1, 0.1, 0.1],
        )
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        # Use params that produce negative predictions (very large k2r relative to k2f)
        pn = EnzymeRates.fitted_params(uni_uni)
        np = length(pn)
        x_bad = zeros(np)
        for (i, p) in enumerate(pn)
            if p == :k2r
                x_bad[i] = 15.0  # exp(15) ≈ 3.3M
            else
                x_bad[i] = -15.0  # exp(-15) ≈ 3e-7
            end
        end
        l = EnzymeRates.loss!(x_bad, fp)
        @test isfinite(l)
        @test l > 0.0
    end

    # ── Test 7: Zero allocations ──────────────────────────────────────────────
    @testset "Zero allocations" begin
        Keq_val = 2.0
        true_params = (k1f = 10.0, k2f = 5.0, k2r = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [(S = Float64(i), P = 0.1) for i in 1:20]
        data = make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        x = randn(length(EnzymeRates.fitted_params(uni_uni)))
        EnzymeRates.loss!(x, fp)  # warmup
        allocs = @allocated EnzymeRates.loss!(x, fp)
        @test allocs == 0
    end

    # ── Test 8: Speed benchmark ───────────────────────────────────────────────
    @testset "Speed" begin
        # Build a larger mechanism: Ordered Bi-Bi
        ordered_bi_bi = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products:   P[C], Q[N]
                enzymes:    E, EA[C], EABEPQ[CN], EQ[N]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA, B] <--> [EABEPQ]
                [EABEPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
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
        articles = [string("A", div(i-1, 50)+1) for i in 1:n_points]
        figs = [string("F", mod(i-1, 5)+1) for i in 1:n_points]

        rates = [rate_equation(ordered_bi_bi, c, bb_true_params) for c in concs_list]
        met_names_bb = metabolites(ordered_bi_bi)
        data = (
            Article = articles,
            Fig = figs,
            Rate = rates,
            (mn => [c[mn] for c in concs_list] for mn in met_names_bb)...
        )
        fp = FittingProblem(ordered_bi_bi, data; Keq=Keq_val)

        x = randn(length(EnzymeRates.fitted_params(ordered_bi_bi)))
        EnzymeRates.loss!(x, fp)  # warmup

        # Time 1000 calls
        t = @elapsed for _ in 1:1000; EnzymeRates.loss!(x, fp); end
        avg_us = t / 1000 * 1e6
        @test avg_us < 50  # < 50 μs per call for 500 datapoints (~5 μs typical)
    end

    # ── Test 9: Validation errors ─────────────────────────────────────────────
    @testset "Validation errors" begin
        # Missing Rate column
        data_no_rate = (Article = ["A1"], Fig = ["F1"], S = [1.0], P = [0.1])
        @test_throws ErrorException FittingProblem(uni_uni, data_no_rate; Keq=1.0)

        # Missing metabolite column
        data_no_met = (Article = ["A1"], Fig = ["F1"], Rate = [1.0], S = [1.0])
        @test_throws ErrorException FittingProblem(uni_uni, data_no_met; Keq=1.0)

        # Zero rate
        data_zero = (Article = ["A1"], Fig = ["F1"], Rate = [0.0], S = [1.0], P = [0.1])
        @test_throws ErrorException FittingProblem(uni_uni, data_zero; Keq=1.0)

        # Missing Article column
        data_no_art = (Fig = ["F1"], Rate = [1.0], S = [1.0], P = [0.1])
        @test_throws ErrorException FittingProblem(uni_uni, data_no_art; Keq=1.0)
    end

end
