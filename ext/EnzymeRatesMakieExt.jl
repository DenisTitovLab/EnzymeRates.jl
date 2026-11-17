# ABOUTME: Draws the panels plot_rate_equation_on_data prepares, one Makie axis per group.
# ABOUTME: Julia loads it with Makie; every number arrives in display units, ready to draw.
module EnzymeRatesMakieExt

using EnzymeRates: RateEquationPanel
using Makie

const _MARKERS = [:circle, :rect, :diamond, :utriangle, :dtriangle, :pentagon, :hexagon,
                  :star5]
const _LINE_WIDTH = 1.6
const _MARKER_STYLE = (; markersize = 7, strokecolor = :white, strokewidth = 0.5)

# plot_rate_equation_on_data calls this by name through `Base.get_extension`.
function draw_rate_equation_panels(panels::Vector{RateEquationPanel}; columns::Integer,
                                   figure::NamedTuple, axis::NamedTuple)
    fs = get(figure, :fontsize, 10)
    used = min(columns, length(panels))
    rows = cld(length(panels), used)
    linear_x = get(axis, :xscale, identity) === identity
    with_theme(_panel_theme(fs, linear_x)) do
        fig = Figure(; size = (320used, 200rows), figure_padding = 10, figure...)
        for (i, panel) in enumerate(panels)
            row, column = fldmod1(i, used)
            _draw_panel!(GridLayout(fig[row, column]), panel, fs;
                ylabel = column == 1 ? "Rate" : "", axis)
        end
        colgap!(fig.layout, 12)
        rowgap!(fig.layout, 10)
        fig
    end
end

function _panel_theme(fs, linear_x)
    ticks = WilkinsonTicks(3; k_min = 3, k_max = 4)
    Theme(; fontsize = fs,
        Axis = (; titlefont = :regular, titlesize = fs, titlegap = 2,
            xgridvisible = false, ygridvisible = false, topspinevisible = false,
            rightspinevisible = false, xticklabelsize = fs - 1, yticklabelsize = fs - 1,
            xticksize = 3, yticksize = 3, xlabelpadding = 1, ylabelpadding = 2,
            xticks = linear_x ? ticks : Makie.automatic, yticks = ticks,
            xtickformat = linear_x ? _tick_labels : Makie.automatic,
            ytickformat = _tick_labels))
end

# "0.5" and "1" where Makie writes "0.50" and "1.0". Adding 0.0 turns -0.0 into 0.0.
_tick_labels(values) = [replace(chopsuffix(string(round(v; sigdigits = 4) + 0.0), ".0"),
                                "-" => "−") for v in values]

function _draw_panel!(layout, panel::RateEquationPanel, fs; ylabel, axis)
    ax = Axis(layout[1, 1]; title = panel.title, xlabel = panel.x_label, ylabel,
        limits = (nothing, panel.rate_limits), axis...)
    colors = Makie.wong_colors()
    entries = map(enumerate(panel.series)) do (i, s)
        color = colors[mod1(i, length(colors))]
        marker = _MARKERS[mod1(i, length(_MARKERS))]
        isempty(s.curve_x) || lines!(ax, s.curve_x, s.curve_rate; color,
            linewidth = _LINE_WIDTH)
        scatter!(ax, s.x, s.rate; color, marker, _MARKER_STYLE...)
        point = MarkerElement(; color, marker, _MARKER_STYLE...)
        isempty(s.curve_x) ? [point] :
            [LineElement(; color, linewidth = _LINE_WIDTH), point]
    end
    panel.rate_power == 0 || _power_of_ten_note!(ax, panel, fs)

    side = GridLayout(layout[1, 2]; tellheight = false, valign = :top, halign = :left)
    fixed = !isempty(panel.fixed_concentrations)
    fixed && Label(side[1, 1], join(panel.fixed_concentrations, "\n"); fontsize = fs - 1,
        justification = :left, halign = :left, color = :gray25, lineheight = 1.05)
    isempty(panel.legend_title) || Legend(side[fixed ? 2 : 1, 1], entries,
        [s.label for s in panel.series], panel.legend_title; framevisible = false,
        labelsize = fs - 1, titlesize = fs - 1, titlefont = :regular, titlehalign = :left,
        halign = :left, patchsize = (1.8(fs - 1), 0.9(fs - 1)), rowgap = 0, titlegap = 1,
        padding = (0, 0, 0, 0), patchlabelgap = 3)
    rowgap!(side, 3)
    colgap!(layout, 2)
end

# The note sits on the data's side of zero, where the axis has room for it.
function _power_of_ten_note!(ax, panel, fs)
    low, high = panel.rate_limits
    top = high > -low
    text!(ax, 0.02, top ? 0.98 : 0.02; space = :relative,
        align = (:left, top ? :top : :bottom),
        text = rich("×10", superscript(string(panel.rate_power))), fontsize = fs - 1,
        color = :gray35)
end

end
