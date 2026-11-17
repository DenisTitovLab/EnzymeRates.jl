# ABOUTME: Tests the figure plot_rate_equation_on_data draws through its Makie extension.
# ABOUTME: Draws written-out tables with CairoMakie and reads the figure's blocks and plots.

using Test
using EnzymeRates
using EnzymeRates: rate_equation_panels
using CairoMakie

# "exp_1": A at four values for each B, Q fixed. "exp_2": A varies, B fixed. "exp_3": B
# varies, A fixed.
const _testhelper_makie_table = _testhelper_plot_table(
    [fill("exp_1", 8); fill("exp_2", 3); fill("exp_3", 3)],
    [[(a, b, 0.0, 5e-5) for b in (1e-3, 2e-3) for a in (1e-4, 2e-4, 5e-4, 1e-3)];
     [(a, 1e-3, 0.0, 0.0) for a in (1e-4, 2e-4, 5e-4)];
     [(1e-3, b, 0.0, 0.0) for b in (1e-4, 2e-4, 5e-4)]])

_testhelper_makie_plot(table = _testhelper_makie_table; keywords...) =
    plot_rate_equation_on_data(_testhelper_plot_m, _testhelper_plot_params, table;
        Keq = _testhelper_plot_keq, keywords...)

_testhelper_blocks(T, fig) = filter(b -> b isa T, fig.content)
_testhelper_figure_size(fig) = Tuple(widths(fig.scene.viewport[]))
_testhelper_plots(T, ax) = filter(p -> p isa T, ax.scene.plots)

# One group of the three-group table, drawn alone.
function _testhelper_makie_group(name)
    rows = _testhelper_makie_table.group .== name
    _testhelper_makie_plot(map(column -> column[rows], _testhelper_makie_table))
end

@testset "plot_rate_equation_on_data" begin
    @testset "returns a figure with one axis per group" begin
        fig = _testhelper_makie_plot()
        @test fig isa Makie.Figure
        panel_axes = _testhelper_blocks(Makie.Axis, fig)
        @test length(panel_axes) == 3
        @test [ax.title[] for ax in panel_axes] == ["exp 1", "exp 2", "exp 3"]
        @test [ax.xlabel[] for ax in panel_axes] == ["A (mM)", "A (µM)", "B (µM)"]
        panels = rate_equation_panels(_testhelper_plot_m, _testhelper_plot_params,
            _testhelper_makie_table; Keq = _testhelper_plot_keq, absolute = false,
            concentration_unit = "M", xscale = identity)
        @test [ax.limits[][2] for ax in panel_axes] == [p.rate_limits for p in panels]
    end

    @testset "grid" begin
        fig = _testhelper_makie_plot(columns = 2)
        @test _testhelper_figure_size(fig) == (640, 400)
        @test [ax.ylabel[] for ax in _testhelper_blocks(Makie.Axis, fig)] ==
              ["Rate", "", "Rate"]
        fig = _testhelper_makie_plot()
        @test _testhelper_figure_size(fig) == (960, 200)
        @test [ax.ylabel[] for ax in _testhelper_blocks(Makie.Axis, fig)] ==
              ["Rate", "", ""]
    end

    @testset "figure and axis keywords" begin
        fig = _testhelper_makie_plot(figure = (; size = (500, 300), fontsize = 8))
        @test _testhelper_figure_size(fig) == (500, 300)
        ax = first(_testhelper_blocks(Makie.Axis, fig))
        @test ax.titlesize[] == 8
        @test ax.xticklabelsize[] == 7

        panel_axes = _testhelper_blocks(Makie.Axis, _testhelper_makie_plot())
        @test !any(ax.xgridvisible[] for ax in panel_axes)
        panel_axes = _testhelper_blocks(Makie.Axis,
            _testhelper_makie_plot(axis = (; xgridvisible = true)))
        @test all(ax.xgridvisible[] for ax in panel_axes)
    end

    @testset "tick labels carry no trailing zeros" begin
        ax = first(_testhelper_blocks(Makie.Axis, _testhelper_makie_plot()))
        values = [-5.0, -0.0, 0.0, 0.25, 0.5, 1.0, 1500.0, 0.1 + 0.2]
        labels = ["−5", "0", "0", "0.25", "0.5", "1", "1500", "0.3"]
        @test ax.xtickformat[](values) == labels
        @test ax.ytickformat[](values) == labels
        ax = first(_testhelper_blocks(Makie.Axis,
            _testhelper_makie_plot(axis = (; xscale = log10))))
        @test ax.xtickformat[] === Makie.automatic
    end

    @testset "log x" begin
        panel_axes = _testhelper_blocks(Makie.Axis,
            _testhelper_makie_plot(axis = (; xscale = log10)))
        @test all(ax.xscale[] === log10 for ax in panel_axes)

        zero_x = _testhelper_plot_table(fill("zero_x", 2),
            [(0.0, 1e-3, 0.0, 0.0), (1e-4, 1e-3, 0.0, 0.0)]; rates = [1.0, 2.0])
        @test_throws ArgumentError _testhelper_makie_plot(zero_x; axis = (; xscale = log10))
    end

    @testset "rejected arguments" begin
        no_rows = map(column -> column[1:0], _testhelper_makie_table)
        @test_throws ErrorException("Table has no rows") _testhelper_makie_plot(no_rows)
        for columns in (0, -1)
            @test_throws(ArgumentError("columns must be at least 1; got $columns"),
                _testhelper_makie_plot(; columns))
        end
    end

    @testset "global theme unchanged" begin
        before = Makie.current_default_theme().fontsize[]
        _testhelper_makie_plot(figure = (; fontsize = 7))
        @test Makie.current_default_theme().fontsize[] == before
        @test Axis(Figure()[1, 1]).xgridvisible[]
    end

    @testset "saves" begin
        fig = _testhelper_makie_plot()
        mktempdir() do dir
            for file in ("panels.png", "panels.pdf")
                path = joinpath(dir, file)
                save(path, fig)
                @test filesize(path) > 0
            end
        end
    end

    @testset "legend and labels" begin
        fig = _testhelper_makie_group("exp_1")
        legend = only(_testhelper_blocks(Makie.Legend, fig))
        title, entries = only(legend.entrygroups[])
        @test title == "B"
        @test [e.label[] for e in entries] == ["1 mM", "2 mM"]
        @test [l.text[] for l in _testhelper_blocks(Makie.Label, fig)] == ["Q 50 µM"]

        fig = _testhelper_makie_group("exp_2")
        @test isempty(_testhelper_blocks(Makie.Legend, fig))
        @test [l.text[] for l in _testhelper_blocks(Makie.Label, fig)] == ["B 1 mM"]
    end

    @testset "series styles cycle by index" begin
        # Nine series: past seven colours and eight markers, both start again.
        rows = [(a, b, 0.0, 0.0) for b in (1:9) .* 1e-4 for a in (1:10) .* 1e-4]
        table = _testhelper_plot_table(fill("nine", 90), rows)
        ax = only(_testhelper_blocks(Makie.Axis, _testhelper_makie_plot(table)))
        @test length(_testhelper_plots(Lines, ax)) == 9
        scatters = _testhelper_plots(Scatter, ax)
        markers = [:circle, :rect, :diamond, :utriangle, :dtriangle, :pentagon, :hexagon,
                   :star5, :circle]
        @test [s.marker[] for s in scatters] == Makie.to_spritemarker.(markers)
        @test [s.color[] for s in scatters] == Makie.wong_colors()[[1:7; 1:2]]
    end

    @testset "power of ten note" begin
        ax = only(_testhelper_blocks(Makie.Axis, _testhelper_makie_group("exp_1")))
        @test isempty(_testhelper_plots(Makie.Text, ax))

        forward = [(a, 1e-3, 0.0, 0.0) for a in (1e-4, 2e-4, 5e-4)]
        backward = [(0.0, 0.0, p, 5e-4) for p in (1e-4, 2e-4, 5e-4)]
        for (rows, rates, align) in ((forward, [3e-5, 5e-5, 8e-5], (:left, :top)),
                                     (backward, [-3e-5, -5e-5, -8e-5], (:left, :bottom)))
            table = _testhelper_plot_table(fill("small", 3), rows; rates)
            ax = only(_testhelper_blocks(Makie.Axis, _testhelper_makie_plot(table)))
            note = only(_testhelper_plots(Makie.Text, ax))
            @test note.align[] == align
        end
    end
end
