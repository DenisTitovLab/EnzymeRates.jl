# ABOUTME: plot_rate_equation_on_data and the plain description of each measurement group's
# ABOUTME: panel it prepares; ext/EnzymeRatesMakieExt.jl draws the panels with Makie.

struct RateEquationSeries
    label::String
    x::Vector{Float64}
    rate::Vector{Float64}
    curve_x::Vector{Float64}
    curve_rate::Vector{Float64}
end

struct RateEquationPanel
    title::String
    x_label::String
    rate_power::Int
    fixed_concentrations::Vector{String}
    legend_title::String
    series::Vector{RateEquationSeries}
    rate_limits::Tuple{Float64, Float64}
end

# The divisor and unit label that show `largest` (a concentration in `unit`): under "M" the
# prefix that fits it, otherwise the unit as given. The prefix fits the value rounded to the
# four digits it prints with, so 999.96 µM reads 1 mM and not 1000 µM.
function _display_unit(largest::Real, unit::AbstractString)
    unit == "M" || return 1.0, String(unit)
    shown = round(largest; sigdigits = 4)
    shown >= 1 ? (1.0, "M") :
    shown >= 1e-3 ? (1e-3, "mM") :
    shown >= 1e-6 ? (1e-6, "µM") :
    shown > 0 ? (1e-9, "nM") : (1.0, "M")
end

function _concentration_text(value::Real, unit::AbstractString)
    value == 0 && return "0"
    divisor, label = _display_unit(value, unit)
    text = string(round(value / divisor; sigdigits = 4))
    "$(chopsuffix(text, ".0")) $label"
end

# The y limits that show a group's measured rates and, up to 1.6 times them, its curves.
function _rate_limits(rate, curve_rate)
    hi, lo = maximum(rate), minimum(rate)
    finite = filter(isfinite, curve_rate)
    chi, clo = isempty(finite) ? (hi, lo) : (maximum(finite), minimum(finite))
    top = hi > 0 ? max(1.12hi, min(1.6hi, 1.05chi)) : 0.0
    bottom = lo < 0 ? min(1.12lo, max(1.6lo, 1.05clo)) : 0.0
    bottom == 0 ? (-0.05top, top) : top == 0 ? (bottom, -0.05bottom) : (bottom, top)
end

# `exp(mean(log(|measured / predicted|)))` over the rows with a finite nonzero prediction,
# or `nothing` when there is none.
function _group_scale(measured, predicted)
    usable = findall(p -> isfinite(p) && p != 0, predicted)
    isempty(usable) ? nothing :
    exp(mean(log(abs(measured[i] / predicted[i])) for i in usable))
end

const _AnyMechanism = Union{Mechanisms.AbstractEnzymeMechanism, Mechanisms.Mechanism,
                            Mechanisms.AllostericMechanism}

"""
    rate_equation_panels(predict, listed_first, Keq, table; absolute, concentration_unit,
        xscale)

One `RateEquationPanel` per distinct `group` of `table`, in order of first appearance, with
the measured points and the curves of `predict`, the rate at a `NamedTuple` of one row's
concentrations, as `plot_rate_equation_on_data` draws them. The concentration columns
`listed_first` come ahead of the table's others, and `Keq` is checked with the table.
`xscale` other than `identity` starts each curve at the group's smallest x and refuses a
group with an x of zero.
"""
function rate_equation_panels(predict, listed_first, Keq::Real, table; absolute::Bool,
                              concentration_unit::AbstractString, xscale)
    columns = Tables.columntable(table)
    shown = (listed_first...,
             filter(!in((:group, :Rate, listed_first...)), keys(columns))...)
    data = Data.RateData(columns, shown, Keq; absolute).data
    [_rate_equation_panel(predict, data, shown, findall(isequal(g), data.group),
         absolute, concentration_unit, xscale)
     for g in unique(data.group)]
end

function _rate_equation_panel(predict, data, mets, rows, absolute, unit, xscale)
    columns = NamedTuple{mets}(map(k -> data[k][rows], mets))
    row_concentrations = [map(column -> column[i], columns) for i in eachindex(rows)]

    x_metabolite = argmax(k -> length(unique(columns[k])), mets)
    others = filter(!=(x_metabolite), mets)
    varying = filter(k -> !allequal(columns[k]), others)
    fixed = filter(k -> k ∉ varying && first(columns[k]) != 0, others)

    group_name = string(data.group[first(rows)])
    x = columns[x_metabolite]
    xscale === identity || any(iszero, x) && throw(ArgumentError(
        "Group $group_name has $x_metabolite = 0, and a log x axis cannot show zero"))
    largest_x = maximum(x)
    divisor, x_unit = _display_unit(largest_x, unit)

    measured = data.Rate[rows]
    scale = absolute ? 1.0 : _group_scale(measured, predict.(row_concentrations))
    drawn = scale !== nothing && largest_x > 0
    largest_rate = maximum(abs, measured)
    rate_power = largest_rate < 0.01 ? floor(Int, log10(largest_rate)) : 0
    unit_of_rate = 10.0^rate_power

    curve_x = !drawn ? Float64[] :
        xscale === identity ? collect(range(0, largest_x, 200)) :
        exp.(range(log(minimum(x)), log(largest_x), 200))
    series_of_row = [Tuple(concs[k] for k in varying) for concs in row_concentrations]
    series = map(sort!(unique(series_of_row))) do key
        in_series = findall(isequal(key), series_of_row)
        concs = row_concentrations[first(in_series)]
        curve_rate = [scale / unit_of_rate * predict(merge(concs, (; x_metabolite => cx)))
                      for cx in curve_x]
        RateEquationSeries(join((_concentration_text(v, unit) for v in key), ", "),
            x[in_series] / divisor, measured[in_series] / unit_of_rate,
            curve_x / divisor, curve_rate)
    end
    all_curves = reduce(vcat, (s.curve_rate for s in series))
    RateEquationPanel(replace(group_name, "_" => " "), "$x_metabolite ($x_unit)",
        rate_power, ["$k $(_concentration_text(first(columns[k]), unit))" for k in fixed],
        join(varying, ", "), series, _rate_limits(measured / unit_of_rate, all_curves))
end

"""
    plot_rate_equation_on_data(m, params, table, Keq; absolute = false,
        concentration_unit = "M", columns = 5, figure = (;), axis = (;)) -> Makie.Figure
    plot_rate_equation_on_data(rate::Function, params, table; absolute = false,
        concentration_unit = "M", columns = 5, figure = (;), axis = (;)) -> Makie.Figure

Draws the rate equation of the mechanism `m`, with the equilibrium constant `Keq`, or a
rate function `rate` of your own, at `params` over the measured rates of `table`, one
panel per `group` in order of first appearance, `columns` panels per row, and returns the
`Makie.Figure`. It fits nothing. It needs a Makie backend: run `using CairoMakie` (or
another backend) first.

```julia
using EnzymeRates, CairoMakie
fig = plot_rate_equation_on_data(m, fit.params, data, 2.0)
save("fit.pdf", fig)

michaelis_menten(concentrations, params) =
    params.Vmax * concentrations.S / (params.Km + concentrations.S)
fig = plot_rate_equation_on_data(michaelis_menten, (Vmax = 2.0, Km = 1e-4), data)
```

# Arguments
- `m`: a mechanism, in any form `FittingProblem` accepts, which draws its rate law,
  `compile_rate_equation(m)(concentrations, params, Keq)`; or `rate`, a function, which
  draws `rate(concentrations, params)`. `concentrations` is a `NamedTuple` of one row's
  concentrations, one per concentration column of `table`.
- `params`: a `NamedTuple`. With a mechanism, it holds every name of `parameters(m)`, as a
  fit's `params` does; other names are ignored. A function receives it as given.
- `table`: a Tables.jl table, such as a `DataFrame`, a `NamedTuple` of columns or the
  result of `CSV.File`, with a `group` column, a `Rate` column and concentration
  columns. Every other column is a concentration, so remove columns that hold anything
  else, such as a pH or a source name. A mechanism needs a column for each of its
  metabolites and ignores the concentrations of metabolites it does not bind; the panels
  still show them, so that measurements differing only in one of them are not drawn as
  one series. The table must pass the checks of `FittingProblem`.
- `Keq`: the reaction's equilibrium constant, with a mechanism; a rate function holds its
  own.
- `absolute`: `false` for rates relative within each group, `true` for absolute
  per-enzyme turnover rates.
- `concentration_unit`: the unit of the table's concentrations. With `"M"`, each x axis
  shows M, mM, µM or nM to fit its largest value, and each legend value picks its own
  prefix. Any other unit is shown as given, and the numbers are not converted.
- `columns`: panels per row.
- `figure`: keywords for `Makie.Figure`, such as `(; size = (510, 330), fontsize = 7.5)`.
  Every text size follows `fontsize` (10 by default); the default size is 320 × 200 per
  panel.
- `axis`: keywords for every panel's `Makie.Axis`, such as `(; xgridvisible = true)`.

Values in `figure` and `axis` override the built-in defaults.

# A panel
The x axis is the metabolite with the most distinct concentrations in the group; a tie goes
to the earlier one in `metabolites(m)`, then to the earlier column of `table`. Each
distinct combination of the other metabolites' concentrations is a series with its own
colour and marker: the measured `Rate` as points and the rate equation as a curve, a
mechanism's per unit enzyme. Where the rate equation does not respond to a metabolite,
the curve along it is flat, and series that differ only in it share a curve. With
`absolute = false`, each group's curves are multiplied by one factor, the geometric mean
of `|Rate / predicted|` over the group's rows; no parameter changes. Beside the panel
stand the metabolites held at one nonzero concentration, then a legend headed by those
that vary between series.

The y axis, in the table's units of `Rate`, includes zero and shows curves up to 1.6 times
the largest rate. When a group's largest absolute rate is below 0.01, the axis counts in a
power of ten, noted in the panel's corner. The x axis is linear from zero unless `axis`
sets `xscale`, such as `log10`; then each curve runs from the group's smallest x to its
largest.

# Errors
- `ErrorException` when no Makie backend is loaded, or when `params` lacks a fitted
  parameter of the mechanism.
- `ArgumentError` when `columns` is less than 1.
- `ArgumentError` when `axis` sets an `xscale` other than `identity` and a group has an x
  of zero, which a log axis cannot show.
- The errors `FittingProblem` raises for a table it rejects, with every column but `group`
  and `Rate` checked as a concentration.
"""
function plot_rate_equation_on_data(@nospecialize(m::_AnyMechanism),
                                    @nospecialize(params::NamedTuple), table, Keq::Real;
                                    absolute::Bool = false,
                                    concentration_unit::AbstractString = "M",
                                    columns::Integer = 5, figure::NamedTuple = (;),
                                    axis::NamedTuple = (;))
    makie_extension = _makie_extension(columns)
    rate = Derivation.compile_rate_equation(m)
    for name in Derivation.parameters(rate)
        haskey(params, name) || error("params is missing $name")
    end
    panels = rate_equation_panels(concs -> rate(concs, params, Keq),
        Mechanisms.metabolites(rate), Keq, table; absolute, concentration_unit,
        xscale = get(axis, :xscale, identity))
    makie_extension.draw_rate_equation_panels(panels; columns, figure, axis)
end

function plot_rate_equation_on_data(@nospecialize(rate::Function),
                                    @nospecialize(params::NamedTuple), table;
                                    absolute::Bool = false,
                                    concentration_unit::AbstractString = "M",
                                    columns::Integer = 5, figure::NamedTuple = (;),
                                    axis::NamedTuple = (;))
    makie_extension = _makie_extension(columns)
    # A rate function holds its own Keq, so the table is checked with a placeholder.
    panels = rate_equation_panels(concs -> rate(concs, params), (), 1.0, table; absolute,
        concentration_unit, xscale = get(axis, :xscale, identity))
    makie_extension.draw_rate_equation_panels(panels; columns, figure, axis)
end

# The extension that draws the panels, its `draw_rate_equation_panels`, once a Makie
# backend is loaded and `columns` is at least 1: both are checked before the table. It is
# not a function declared here with its method in the extension: without Makie that call
# has no method, which JET's package test reports, and with Makie JET would follow it into
# the drawing, so test_type_stability.jl would depend on what is loaded.
function _makie_extension(columns::Integer)
    makie_extension = Base.get_extension(parentmodule(Plotting), :EnzymeRatesMakieExt)
    makie_extension === nothing && error(
        "plot_rate_equation_on_data needs a Makie backend: run `using CairoMakie` " *
        "(or another Makie backend) first")
    columns >= 1 || throw(ArgumentError("columns must be at least 1; got $columns"))
    makie_extension
end
