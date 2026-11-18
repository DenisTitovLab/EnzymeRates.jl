# ABOUTME: Tests for the Makie-free panel preparation behind plot_rate_equation_on_data.
# ABOUTME: Expected values come from rate_equation called directly on written-out tables.

using Test
using EnzymeRates
using EnzymeRates.Plotting: rate_equation_panels

function _testhelper_plot_panels(table; rate = _testhelper_plot_m,
                                 params = _testhelper_plot_params, absolute = false,
                                 unit = "M", xscale = identity)
    rate_equation_panels(rate, params, table; absolute, concentration_unit = unit, xscale)
end

# Group "exp_1": A at four values for each B in (1e-3, 2e-3), Q = 5e-5, P = 0.
const _testhelper_plot_rows = [(a, b, 0.0, 5e-5) for b in (1e-3, 2e-3)
                               for a in (1e-4, 2e-4, 5e-4, 1e-3)]
_testhelper_plot_exp1(; factor = 3.0) =
    _testhelper_plot_table(fill("exp_1", 8), _testhelper_plot_rows; factor)

@testset "rate_equation_panels" begin
    @testset "x metabolite and series" begin
        panels = _testhelper_plot_panels(_testhelper_plot_exp1())
        @test length(panels) == 1
        p = only(panels)
        @test p.title == "exp 1"
        @test p.x_label == "A (mM)"
        @test p.legend_title == "B"
        @test p.fixed_concentrations == ["Q 50 µM"]
        @test [s.label for s in p.series] == ["1 mM", "2 mM"]
        @test p.series[1].x ≈ [0.1, 0.2, 0.5, 1.0]
        @test p.rate_power == 0
    end

    @testset "tie goes to metabolites order" begin
        rows = [(a, b, 0.0, 5e-5) for a in (1e-4, 2e-4) for b in (1e-3, 2e-3)]
        p = only(_testhelper_plot_panels(
            _testhelper_plot_table(fill("g", 4), rows)))
        @test startswith(p.x_label, string(first(metabolites(_testhelper_plot_m))))
        @test p.legend_title == "B"
    end

    @testset "group scale" begin
        for (absolute, factor) in ((false, 3.0), (true, 1.0))
            table = _testhelper_plot_exp1(; factor)
            p = only(_testhelper_plot_panels(table; absolute))
            for (s, b) in zip(p.series, (1e-3, 2e-3))
                @test length(s.curve_x) == 200
                @test s.curve_x[1] == 0
                @test s.curve_x[end] ≈ 1.0
                expected = [factor * _testhelper_plot_rate(
                    (A = x * 1e-3, B = b, P = 0.0, Q = 5e-5)) for x in s.curve_x]
                @test s.curve_rate ≈ expected rtol = 1e-10
            end
        end
    end

    @testset "group scale is the geometric mean of absolute ratios" begin
        rows = [(1e-4, 1e-3, 0.0, 5e-5), (1e-3, 1e-3, 0.0, 5e-5)]
        predicted = [_testhelper_plot_rate((A = r[1], B = r[2], P = r[3], Q = r[4]))
                     for r in rows]
        # Ratios 2 and 8 give 4 (the arithmetic mean and the median are 5); ratios 2 and
        # -8 give 4 as well, because the ratios enter as absolute values.
        for ratios in ([2.0, 8.0], [2.0, -8.0])
            table = _testhelper_plot_table(fill("g", 2), rows; rates = ratios .* predicted)
            s = only(only(_testhelper_plot_panels(table)).series)
            expected = [4 * _testhelper_plot_rate(
                (A = x * 1e-3, B = 1e-3, P = 0.0, Q = 5e-5)) for x in s.curve_x]
            @test s.curve_rate ≈ expected rtol = 1e-10
        end
    end

    @testset "panels follow first appearance" begin
        rows = [(1e-4, 1e-3, 0.0, 0.0), (2e-4, 1e-3, 0.0, 0.0),
            (3e-4, 1e-3, 0.0, 0.0)]
        panels = _testhelper_plot_panels(_testhelper_plot_table(["b", "a", "b"], rows))
        @test [p.title for p in panels] == ["b", "a"]
    end

    @testset "small rates" begin
        table = _testhelper_plot_exp1()
        scaled = merge(table, (Rate = table.Rate .* (3e-5 / maximum(table.Rate)),))
        p = only(_testhelper_plot_panels(scaled))
        @test p.rate_power == -5
        @test maximum(maximum(s.rate) for s in p.series) ≈ 3.0
    end

    @testset "rate limits" begin
        hi_of(p) = maximum(maximum(s.rate) for s in p.series)
        lo_of(p) = minimum(minimum(s.rate) for s in p.series)
        curve_hi_of(p) = maximum(maximum(s.curve_rate) for s in p.series)
        curve_lo_of(p) = minimum(minimum(s.curve_rate) for s in p.series)
        # In absolute mode each curve is 1 / factor times the rates it passes through: at
        # factor 3 the data set the limit, at 0.1 the curves stop at 1.6 times the data,
        # and at 0.8 the limit is 5% past the curves.
        forward(factor) = only(_testhelper_plot_panels(_testhelper_plot_exp1(; factor);
            absolute = true))
        p = forward(3.0)
        @test p.rate_limits == (-0.05 * 1.12hi_of(p), 1.12hi_of(p))
        p = forward(0.1)
        @test p.rate_limits == (-0.05 * 1.6hi_of(p), 1.6hi_of(p))
        p = forward(0.8)
        @test p.rate_limits == (-0.05 * 1.05curve_hi_of(p), 1.05curve_hi_of(p))

        # Only reverse rates: A = B = 0.
        reverse_rows = [(0.0, 0.0, p_, 5e-4) for p_ in (1e-4, 2e-4, 5e-4)]
        backward(factor) = only(_testhelper_plot_panels(
            _testhelper_plot_table(fill("rev", 3), reverse_rows; factor); absolute = true))
        p = backward(3.0)
        @test hi_of(p) < 0
        @test p.rate_limits == (1.12lo_of(p), -0.05 * 1.12lo_of(p))
        p = backward(0.1)
        @test p.rate_limits == (1.6lo_of(p), -0.05 * 1.6lo_of(p))
        p = backward(0.8)
        @test p.rate_limits == (1.05curve_lo_of(p), -0.05 * 1.05curve_lo_of(p))

        # P Q / (A B) is 1, 2, 10 and 20 against Keq = 4: two rows run forward, two back,
        # and each side follows the rule.
        mixed_rows = [(1e-4, 1e-4, p_, 1e-3) for p_ in (1e-5, 2e-5, 1e-4, 2e-4)]
        for (factor, margin) in ((3.0, 1.12), (0.1, 1.6))
            p = only(_testhelper_plot_panels(_testhelper_plot_table(fill("mixed", 4),
                mixed_rows; factor); absolute = true))
            @test lo_of(p) < 0 < hi_of(p)
            @test p.rate_limits == (margin * lo_of(p), margin * hi_of(p))
        end
    end

    @testset "two varying metabolites" begin
        # Series (B, P) in table order (2 µM, 0), (1 µM, 10 µM), (1 µM, 0); rows alternate.
        b_and_p = [(2e-6, 0.0), (1e-6, 1e-5), (1e-6, 0.0)]
        rows = [(a, b, p_, 0.0) for a in (1e-4, 2e-4, 5e-4, 1e-3) for (b, p_) in b_and_p]
        table = _testhelper_plot_table(fill("g", 12), rows; rates = collect(1.0:12.0))
        p = only(_testhelper_plot_panels(table))
        @test p.x_label == "A (mM)"
        @test p.legend_title == "B, P"
        @test isempty(p.fixed_concentrations)
        @test [s.label for s in p.series] == ["1 µM, 0", "1 µM, 10 µM", "2 µM, 0"]
        @test [s.rate for s in p.series] ==
              [[3.0, 6.0, 9.0, 12.0], [2.0, 5.0, 8.0, 11.0], [1.0, 4.0, 7.0, 10.0]]
        @test all(s.x ≈ [0.1, 0.2, 0.5, 1.0] for s in p.series)
        # Each curve is drawn at its own series' B and P, all with one group scale.
        scales = [s.curve_rate[end] / _testhelper_plot_rate((A = 1e-3, B = b, P = p_,
                      Q = 0.0)) for (s, (b, p_)) in zip(p.series, b_and_p[[3, 2, 1]])]
        @test all(≈(first(scales)), scales)
    end

    @testset "a one-row group draws a curve through its point" begin
        p = only(_testhelper_plot_panels(
            _testhelper_plot_table(["g"], [(1e-4, 1e-3, 0.0, 0.0)]; rates = [2.0])))
        s = only(p.series)
        @test p.x_label == "A (µM)"
        @test length(s.curve_x) == 200
        @test s.curve_x[end] ≈ 100
        @test s.curve_rate[end] ≈ 2.0
    end

    @testset "concentration units" begin
        rows = [(1e-4, 2.0, 0.0, 0.0), (2e-4, 2.0, 0.0, 0.0)]
        table = _testhelper_plot_table(fill("g", 2), rows)
        p = only(_testhelper_plot_panels(table; unit = "mM"))
        @test p.x_label == "A (mM)"
        @test p.series[1].x == [1e-4, 2e-4]
        @test p.fixed_concentrations == ["B 2 mM"]

        # Rounding to four digits comes first: 999.96 µM reads 1 mM, 0.99996 M reads 1 M.
        for (value, text) in ((1.25e-3, "1.25 mM"), (1.08e-3, "1.08 mM"),
            (5e-9, "5 nM"), (2.5, "2.5 M"), (1.0, "1 M"), (9.9996e-4, "1 mM"),
            (0.99996, "1 M"))
            rows = [(1e-4, value, 0.0, 0.0), (2e-4, value, 0.0, 0.0)]
            p = only(_testhelper_plot_panels(
                _testhelper_plot_table(fill("g", 2), rows)))
            @test p.fixed_concentrations == ["B $text"]
        end

        for (largest, label) in ((2.0, "A (M)"), (0.99996, "A (M)"),
            (9.9996e-4, "A (mM)"))
            rows = [(largest / 2, 1e-3, 0.0, 0.0), (largest, 1e-3, 0.0, 0.0)]
            p = only(_testhelper_plot_panels(
                _testhelper_plot_table(fill("g", 2), rows)))
            @test p.x_label == label
        end
        rows = [(1.0, 1e-3, 0.0, 0.0), (2.0, 1e-3, 0.0, 0.0)]
        p = only(_testhelper_plot_panels(_testhelper_plot_table(fill("g", 2), rows)))
        @test p.series[1].x == [1.0, 2.0]
    end

    @testset "no curve" begin
        p = only(_testhelper_plot_panels(
            _testhelper_plot_table(["g"], [(0.0, 1e-3, 0.0, 0.0)]; rates = [1.0])))
        @test length(p.series) == 1
        @test length(p.series[1].x) == 1
        @test isempty(p.series[1].curve_x) && isempty(p.series[1].curve_rate)

        p = only(_testhelper_plot_panels(
            _testhelper_plot_table(["g", "g"], [(0.0, 1e-3, 0.0, 0.0),
                (0.0, 1e-3, 0.0, 0.0)]; rates = [1.0, 2.0])))
        @test p.x_label == "A (M)"
        @test all(isempty(s.curve_x) for s in p.series)

        # At Keq-equilibrium every prediction is zero.
        eq = [(1.0, 1.0, 2.0, 2.0), (1.0, 1.0, 2.0, 2.0)]
        table = _testhelper_plot_table(["g", "g"], eq; rates = [1.0, 2.0])
        p = only(_testhelper_plot_panels(table))
        @test all(isempty(s.curve_x) && isempty(s.curve_rate) for s in p.series)
    end

    @testset "non-string groups" begin
        rows = [(1e-4, 1e-3, 0.0, 0.0), (2e-4, 1e-3, 0.0, 0.0),
            (1e-4, 1e-3, 0.0, 0.0), (2e-4, 1e-3, 0.0, 0.0)]
        panels = _testhelper_plot_panels(_testhelper_plot_table([1, 1, 2, 2], rows))
        @test [p.title for p in panels] == ["1", "2"]
    end

    @testset "missing and NaN groups" begin
        rows = [(a, 1e-3, 0.0, 0.0) for a in (1e-4, 3e-4, 2e-4, 4e-4)]
        for (groups, titles) in (([missing, "a", missing, "a"], ["missing", "a"]),
                                 ([NaN, 1.0, NaN, 1.0], ["NaN", "1.0"]))
            panels = _testhelper_plot_panels(_testhelper_plot_table(groups, rows))
            @test [p.title for p in panels] == titles
            @test [only(p.series).x for p in panels] ≈ [[100, 200], [300, 400]]
        end
    end

    @testset "log x" begin
        p = only(_testhelper_plot_panels(_testhelper_plot_exp1(); xscale = log10))
        s = p.series[1]
        @test length(s.curve_x) == 200
        @test s.curve_x[1] ≈ 0.1
        @test s.curve_x[end] ≈ 1.0
        @test s.curve_x[2] / s.curve_x[1] ≈ s.curve_x[3] / s.curve_x[2]

        rows = [(0.0, 1e-3, 0.0, 0.0), (1e-4, 1e-3, 0.0, 0.0)]
        table = _testhelper_plot_table(fill("zero_group", 2), rows; rates = [1.0, 2.0])
        err = try
            _testhelper_plot_panels(table; xscale = log10)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("zero_group", err.msg) && occursin("cannot show zero", err.msg)
    end

    @testset "metabolites not in the mechanism" begin
        # A at three values with B fixed, once without and once with 2 mM of I, a
        # metabolite the mechanism does not bind and which halves the measured rate.
        rows = [(a, 1e-3, 0.0, 0.0) for _ in 1:2 for a in (1e-4, 2e-4, 5e-4)]
        inhibitor = [fill(0.0, 3); fill(2e-3, 3)]
        predicted = [_testhelper_plot_rate((A = r[1], B = r[2], P = r[3], Q = r[4]))
                     for r in rows]
        rates = predicted .* [fill(2.0, 3); fill(0.5, 3)]
        table = merge(_testhelper_plot_table(fill("g", 6), rows; rates), (I = inhibitor,))

        p = only(_testhelper_plot_panels(table))
        @test p.x_label == "A (µM)"
        @test p.legend_title == "I"
        @test p.fixed_concentrations == ["B 1 mM"]
        @test [s.label for s in p.series] == ["0", "2 mM"]
        @test p.series[1].rate ≈ rates[1:3]
        @test p.series[2].rate ≈ rates[4:6]
        # The rate equation does not respond to I, so both series share one curve, scaled
        # by the geometric mean of the ratios 2 and 0.5.
        expected = [_testhelper_plot_rate((A = x * 1e-6, B = 1e-3, P = 0.0, Q = 0.0))
                    for x in p.series[1].curve_x]
        @test p.series[1].curve_rate ≈ expected rtol = 1e-10
        @test p.series[2].curve_rate == p.series[1].curve_rate

        # I held at one nonzero concentration stands with the fixed metabolites.
        held = merge(table, (I = fill(2e-3, 6),))
        p = only(_testhelper_plot_panels(held))
        @test p.fixed_concentrations == ["B 1 mM", "I 2 mM"]
        @test length(p.series) == 1

        # I alone varies: it is the x axis and the curve is flat.
        rows = fill((1e-4, 1e-3, 0.0, 0.0), 4)
        table = merge(_testhelper_plot_table(fill("g", 4), rows;
            rates = [4.0, 3.0, 2.0, 1.0]), (I = [0.0, 1e-3, 2e-3, 4e-3],))
        p = only(_testhelper_plot_panels(table))
        @test p.x_label == "I (mM)"
        @test p.fixed_concentrations == ["A 100 µM", "B 1 mM"]
        s = only(p.series)
        @test s.x ≈ [0.0, 1.0, 2.0, 4.0]
        @test s.curve_x[end] ≈ 4.0
        @test allequal(s.curve_rate)
        @test first(s.curve_rate) ≈ exp(sum(log, [4.0, 3.0, 2.0, 1.0]) / 4)

        # In a tie for the x axis the mechanism's metabolites come first, wherever the
        # other column stands in the table.
        tied = merge((I = [1e-3, 2e-3, 3e-3],), _testhelper_plot_table(fill("g", 3),
            [(a, 1e-3, 0.0, 0.0) for a in (1e-4, 2e-4, 5e-4)]))
        @test only(_testhelper_plot_panels(tied)).x_label == "A (µM)"

        # Every column but group and Rate is checked as a concentration.
        negative = merge(table, (I = [0.0, 1e-3, 2e-3, -4e-3],))
        message = "Concentration I at row 4 must be a finite number ≥ 0; got -0.004"
        @test_throws ErrorException(message) _testhelper_plot_panels(negative)
        @test_throws ErrorException _testhelper_plot_panels(
            merge(table, (Article = fill("Smith 1990", 4),)))
        without_B = Base.structdiff(table, NamedTuple{(:B,)})
        @test_throws ErrorException("Missing metabolite column: B") _testhelper_plot_panels(
            without_B)
    end

    @testset "a rate function in place of a mechanism" begin
        # Michaelis–Menten with a competitive inhibitor, its constants in a vector.
        rate(concs, p) = p[1] * concs.S / (p[2] * (1 + concs.I / p[3]) + concs.S)
        params = [2.0, 3e-4, 1e-3]
        rows = [(S = s, I = i) for i in (0.0, 1e-3) for s in (1e-4, 2e-4, 5e-4)]
        table = (group = fill("g", 6), S = [r.S for r in rows], I = [r.I for r in rows],
                 Rate = [3 * rate(r, params) for r in rows])
        panels(; table = table, keywords...) =
            _testhelper_plot_panels(table; rate, params, keywords...)

        for (absolute, factor) in ((false, 3.0), (true, 1.0))
            p = only(panels(; absolute))
            @test p.x_label == "S (µM)"
            @test p.legend_title == "I"
            @test [s.label for s in p.series] == ["0", "1 mM"]
            for (s, i) in zip(p.series, (0.0, 1e-3))
                expected = [factor * rate((S = x * 1e-6, I = i), params) for x in s.curve_x]
                @test s.curve_rate ≈ expected rtol = 1e-10
            end
        end
        # The function receives every column but group and Rate, in the table's order.
        received = Set{Any}()
        recording(concs, p) = (push!(received, keys(concs)); rate(concs, p))
        swapped = (group = table.group, Rate = table.Rate, I = table.S, S = table.I)
        p = only(_testhelper_plot_panels(swapped; rate = recording, params))
        @test received == Set([(:I, :S)])
        @test p.x_label == "I (µM)"

        @test_throws ErrorException panels(table = merge(table, (Note = fill("a", 6),)))
        @test_throws ErrorException("Missing required column: group") panels(
            table = Base.structdiff(table, NamedTuple{(:group,)}))
    end

    @testset "input checks" begin
        table = _testhelper_plot_exp1()
        no_group = Base.structdiff(table, NamedTuple{(:group,)})
        @test_throws ErrorException _testhelper_plot_panels(no_group)
        for absent in (first(keys(_testhelper_plot_params)), :Keq)
            params = Base.structdiff(_testhelper_plot_params, NamedTuple{(absent,)})
            @test_throws ErrorException("params is missing $absent") begin
                _testhelper_plot_panels(table; params)
            end
        end
        @test_throws ErrorException("Keq must be positive and finite; got -4.0") begin
            _testhelper_plot_panels(table;
                params = merge(_testhelper_plot_params, (Keq = -4.0,)))
        end

        no_rows = map(column -> column[1:0], table)
        @test_throws ErrorException("Table has no rows") _testhelper_plot_panels(no_rows)
    end
end

# runtests.jl includes this file before the drawing tests load CairoMakie.
@testset "needs a Makie backend" begin
    @test Base.get_extension(EnzymeRates, :EnzymeRatesMakieExt) === nothing
    message = "plot_rate_equation_on_data needs a Makie backend: run `using CairoMakie` " *
              "(or another Makie backend) first"
    @test_throws ErrorException(message) plot_rate_equation_on_data(_testhelper_plot_m,
        _testhelper_plot_params, _testhelper_plot_exp1())
    # A table the preparation would reject gets the same error, so the check comes first.
    @test_throws ErrorException(message) plot_rate_equation_on_data(_testhelper_plot_m,
        _testhelper_plot_params, (; Rate = [1.0]))
end
