# ABOUTME: Turns a mechanism, parameters and a rate table into one plain description per
# ABOUTME: measurement group, for plotting; it holds no plotting code.

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
# prefix that fits it, otherwise the unit as given.
function _display_unit(largest::Real, unit::AbstractString)
    unit == "M" || return 1.0, String(unit)
    largest >= 1e-3 ? (1e-3, "mM") :
    largest >= 1e-6 ? (1e-6, "µM") :
    largest > 0 ? (1e-9, "nM") : (1.0, "M")
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

"""
    rate_equation_panels(m, params, table; Keq, absolute, concentration_unit, xscale)

One `RateEquationPanel` per distinct `group` of `table`, in order of first appearance, with
the measured points and the rate equation of `m` at `params` for plotting. `xscale`
other than `identity` starts each curve at the group's smallest x and refuses a group with
an x of zero.
"""
function rate_equation_panels(@nospecialize(m), @nospecialize(params::NamedTuple), table;
                              Keq::Real, absolute::Bool,
                              concentration_unit::AbstractString, xscale)
    fp = Fitting.FittingProblem(m, table; Keq, absolute)
    mech, data = fp.mechanism, fp.rate_data.data
    for name in Derivation.parameters(mech, Derivation.FittedParameters())
        haskey(params, name) || error("params is missing the fitted parameter $name")
    end
    turnover_params = merge(params, (Keq = fp.rate_data.Keq, E_total = 1.0))
    mets = Mechanisms.metabolites(mech)
    [_rate_equation_panel(mech, turnover_params, mets, data, findall(==(g), data.group),
         absolute, concentration_unit, xscale) for g in unique(data.group)]
end

function _rate_equation_panel(mech, params, mets, data, rows, absolute, unit, xscale)
    group_name = string(data.group[first(rows)])
    concentration(k) = data[k][rows]
    predict(concs) = Derivation.rate_equation(mech, NamedTuple{mets}(concs), params)

    distinct = [length(unique(concentration(k))) for k in mets]
    x_metabolite = mets[argmax(distinct)]
    others = filter(!=(x_metabolite), mets)
    varying = filter(k -> length(unique(concentration(k))) > 1, others)
    fixed = filter(k -> k ∉ varying && first(concentration(k)) != 0, others)

    x = concentration(x_metabolite)
    xscale === identity || any(iszero, x) && throw(ArgumentError(
        "Group $group_name has $x_metabolite = 0, and a log x axis cannot show zero"))
    largest_x = maximum(x)
    divisor, x_unit = _display_unit(largest_x, unit)

    measured = data.Rate[rows]
    scale = absolute ? 1.0 :
        _group_scale(measured, [predict(Tuple(Float64(data[k][i]) for k in mets))
                                for i in rows])
    drawn = scale !== nothing && largest_x > 0
    largest_rate = maximum(abs, measured)
    rate_power = largest_rate < 0.01 ? floor(Int, log10(largest_rate)) : 0
    unit_of_rate = 10.0^rate_power

    curve_x = !drawn ? Float64[] :
        xscale === identity ? collect(range(0, largest_x, 200)) :
        exp.(range(log(minimum(x)), log(largest_x), 200))
    series_keys = sort!(unique([Tuple(Float64(concentration(k)[i]) for k in varying)
                                for i in eachindex(rows)]))
    series = map(series_keys) do key
        in_series = [i for i in eachindex(rows)
                     if Tuple(Float64(concentration(k)[i]) for k in varying) == key]
        base = Dict(k => Float64(concentration(k)[first(in_series)]) for k in others)
        curve_rate = [scale / unit_of_rate *
                      predict(Tuple(k == x_metabolite ? cx : base[k] for k in mets))
                      for cx in curve_x]
        RateEquationSeries(join((_concentration_text(v, unit) for v in key), ", "),
            x[in_series] / divisor, measured[in_series] / unit_of_rate,
            curve_x / divisor, curve_rate)
    end
    all_curves = reduce(vcat, (s.curve_rate for s in series))
    RateEquationPanel(replace(group_name, "_" => " "), "$x_metabolite ($x_unit)",
        rate_power,
        ["$k $(_concentration_text(first(concentration(k)), unit))" for k in fixed],
        join(varying, ", "), series,
        _rate_limits(measured / unit_of_rate, all_curves))
end
