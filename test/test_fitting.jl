using Tables

@testset "Fitting" begin

    # ── Helper: build a Uni-Uni mechanism ─────────────────────────────────────
    uni_uni = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> ES
            ES <--> E + P
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
        @test length(fp.group_point_indexes) == 1
        @test length(fp.group_point_indexes[1]) == 5
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

    # ── Test 4: Per-group centering invariance ───────────────────────────────
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
        true_params = (k1f = 10.0, k2f = 5.0, k2r = 1.0, Keq = Keq_val, E_total = 1.0)

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

    # ── Test 6: Sign-mismatch penalty ─────────────────────────────────────────
    @testset "Sign mismatch penalty" begin
        Keq_val = 2.0
        # Create data with positive rates
        data = (
            group = ["G1", "G1", "G1"],
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

    # ── Test 6b: All-mismatch group has positive loss ─────────────────────────
    # Regression test: previously, when all predictions in a group were sign-
    # mismatches, the centering step zeroed every deviation (10-mean(10)=0).
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
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> EA
                EA + B <--> EABEPQ
                EABEPQ <--> EQ + P
                EQ <--> E + Q
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
        EnzymeRates.loss!(x, fp)  # warmup

        # Time 1000 calls
        t = @elapsed for _ in 1:1000; EnzymeRates.loss!(x, fp); end
        avg_us = t / 1000 * 1e6
        @test avg_us < 50  # < 50 μs per call for 500 datapoints (~5 μs typical)
    end

    # ── Test 9: kcat normalization in fit_rate_equation ───────────────
    @testset "kcat normalization" begin
        using OptimizationBBO
        Keq_val = 2.0
        true_params = (k1f = 10.0, k2f = 5.0, k2r = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
            (S = 0.5, P = 0.5), (S = 1.0, P = 0.5), (S = 2.0, P = 0.5),
        ]
        data = make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        # Default kcat=1.0: returned params should have kcat ≈ 1
        result = fit_rate_equation(fp, BBO_adaptive_de_rand_1_bin_radiuslimited();
            n_restarts=3, maxtime=5.0)
        full = merge(result.params, (Keq = Keq_val, E_total = 1.0))
        @test EnzymeRates._kcat_forward(uni_uni, full) ≈ 1.0 rtol=0.01

        # Custom kcat target
        result2 = fit_rate_equation(fp, BBO_adaptive_de_rand_1_bin_radiuslimited();
            n_restarts=3, maxtime=5.0, kcat=42.0)
        full2 = merge(result2.params, (Keq = Keq_val, E_total = 1.0))
        @test EnzymeRates._kcat_forward(uni_uni, full2) ≈ 42.0 rtol=0.01

        # kcat=nothing: raw params (no normalization guarantee)
        result3 = fit_rate_equation(fp, BBO_adaptive_de_rand_1_bin_radiuslimited();
            n_restarts=3, maxtime=5.0, kcat=nothing)
        @test haskey(result3, :params)
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
