# ABOUTME: Tests for the Makie-free panel preparation behind plot_rate_equation_on_data.
# ABOUTME: Expected values come from rate_equation called directly on written-out tables.

using Test
using EnzymeRates
using EnzymeRates: rate_equation_panels

function _testhelper_plot_panels(table; absolute = false, unit = "M", xscale = identity)
    rate_equation_panels(_testhelper_plot_m, _testhelper_plot_params, table;
        Keq = _testhelper_plot_keq, absolute, concentration_unit = unit, xscale)
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
        p = only(_testhelper_plot_panels(_testhelper_plot_exp1(; factor = 3.0);
            absolute = true))
        @test p.rate_limits == (-0.05 * 1.12hi_of(p), 1.12hi_of(p))

        p = only(_testhelper_plot_panels(_testhelper_plot_exp1(; factor = 0.1);
            absolute = true))
        @test p.rate_limits == (-0.05 * 1.6hi_of(p), 1.6hi_of(p))

        rows = [(0.0, 0.0, p_, 5e-4) for p_ in (1e-4, 2e-4, 5e-4)]
        p = only(_testhelper_plot_panels(
            _testhelper_plot_table(fill("rev", 3), rows); absolute = true))
        lo = minimum(minimum(s.rate) for s in p.series)
        @test lo < 0
        @test p.rate_limits[1] ≤ 1.12lo
        @test p.rate_limits[2] == -0.05 * p.rate_limits[1]
    end

    @testset "concentration units" begin
        rows = [(1e-4, 2.0, 0.0, 0.0), (2e-4, 2.0, 0.0, 0.0)]
        table = _testhelper_plot_table(fill("g", 2), rows)
        p = only(_testhelper_plot_panels(table; unit = "mM"))
        @test p.x_label == "A (mM)"
        @test p.series[1].x == [1e-4, 2e-4]
        @test p.fixed_concentrations == ["B 2 mM"]

        for (value, text) in ((1.25e-3, "1.25 mM"), (1.08e-3, "1.08 mM"),
            (5e-9, "5 nM"))
            rows = [(1e-4, value, 0.0, 0.0), (2e-4, value, 0.0, 0.0)]
            p = only(_testhelper_plot_panels(
                _testhelper_plot_table(fill("g", 2), rows)))
            @test p.fixed_concentrations == ["B $text"]
        end
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

    @testset "input checks" begin
        table = _testhelper_plot_exp1()
        no_group = Base.structdiff(table, NamedTuple{(:group,)})
        @test_throws ErrorException _testhelper_plot_panels(no_group)
        short = Base.structdiff(_testhelper_plot_params,
            NamedTuple{(first(keys(_testhelper_plot_params)),)})
        err = try
            rate_equation_panels(_testhelper_plot_m, short, table;
                Keq = _testhelper_plot_keq, absolute = false,
                concentration_unit = "M", xscale = identity)
            nothing
        catch e
            e
        end
        @test err isa ErrorException && occursin("params is missing", err.msg)
    end
end

# runtests.jl includes this file before the drawing tests load CairoMakie.
@testset "needs a Makie backend" begin
    @test Base.get_extension(EnzymeRates, :EnzymeRatesMakieExt) === nothing
    message = "plot_rate_equation_on_data needs a Makie backend: run `using CairoMakie` " *
              "(or another Makie backend) first"
    @test_throws ErrorException(message) plot_rate_equation_on_data(_testhelper_plot_m,
        _testhelper_plot_params, _testhelper_plot_exp1(); Keq = _testhelper_plot_keq)
    # A table the preparation would reject gets the same error, so the check comes first.
    @test_throws ErrorException(message) plot_rate_equation_on_data(_testhelper_plot_m,
        _testhelper_plot_params, (; Rate = [1.0]); Keq = _testhelper_plot_keq)
end
