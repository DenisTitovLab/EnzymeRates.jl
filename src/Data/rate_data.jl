# ABOUTME: RateData: a rate table validated against a metabolite list, with its Keq and
# ABOUTME: the reference anchors that set a relative fit's scale.

"""
    RateData(table, metabolite_names; Keq, absolute = false, reference = nothing,
        owner = "mechanism")

A rate table that has passed its checks, with the constants that go with it.

# Fields
- `data`: `Tables.columntable(table)`
- `Keq`: fixed equilibrium constant
- `absolute`: `false` selects relative data, whose scale the `reference` sets; `true`
  selects absolute per-enzyme turnover
- `reference`: the anchors `concentrations => turnover` that set a relative fit's scale,
  with Float64 turnovers; empty in absolute mode and when a relative table gives none

`Keq` must be positive and finite, the table needs a `group` column, a `Rate` column and
one column per name in `metabolite_names`, no metabolite may be named `group` or `Rate`,
every concentration must be a finite number ≥ 0 (zero is valid), and every rate must be a
finite, nonzero number (the loss takes its log). In relative mode `reference` is optional
and, when given, checked against `metabolite_names` of the `owner` (`"mechanism"` or
`"reaction"`); in absolute mode it must be `nothing`.
"""
struct RateData{D<:NamedTuple}
    data::D
    Keq::Float64
    absolute::Bool
    reference::Vector{Pair{NamedTuple, Float64}}
end

function RateData(table, metabolite_names; Keq::Real, absolute::Bool = false,
        reference = nothing, owner::String = "mechanism")
    data = _rate_table(table, metabolite_names, Keq)
    RateData{typeof(data)}(data, Keq, absolute,
        _reference(reference, absolute, metabolite_names, owner))
end

# The validated `Tables.columntable(table)`; see `RateData`.
function _rate_table(table, mnames, Keq)
    0 < Keq < Inf || error("Keq must be positive and finite; got $Keq")
    data = Tables.columntable(table)
    for req in (:group, :Rate)
        req in mnames && error("Metabolite $req has the name of the required $req " *
                               "column; rename the metabolite")
        req in keys(data) || error("Missing required column: $req")
    end
    for m in mnames
        m in keys(data) || error("Missing metabolite column: $m")
        i = findfirst(c -> !(c isa Real && isfinite(c) && c >= 0), data[m])
        i === nothing || error("Concentration $m at row $i must be a finite number ≥ 0; " *
                               "got $(repr(data[m][i]))")
    end
    i = findfirst(r -> !(r isa Real && isfinite(r)), data.Rate)
    i === nothing ||
        error("Rate at row $i must be a finite number; got $(repr(data.Rate[i]))")
    i = findfirst(iszero, data.Rate)
    i === nothing || error("Zero rate at row $i: log(0) is undefined")
    data
end

# The reference anchors a problem stores; see `RateData`.
function _reference(reference, absolute::Bool, mets, owner::String)
    if absolute
        reference === nothing || error(
            "reference applies only to relative mode; absolute turnover data set the " *
            "rate scale themselves")
        return Pair{NamedTuple, Float64}[]
    end
    reference === nothing && return Pair{NamedTuple, Float64}[]
    _check_anchors(reference, mets, owner)
    Pair{NamedTuple, Float64}[concs => Float64(t) for (concs, t) in reference]
end

"""
    _check_anchors(anchors, mets, owner)

Error unless `anchors`, the `reference`, is a nonempty
vector of `concentrations => turnover` pairs whose concentrations are a `NamedTuple` of
finite numbers ≥ 0 keyed by names in `mets`, the metabolites of the `owner`, and whose
turnover is a finite, nonzero number.
"""
function _check_anchors(anchors, mets, owner::String)
    anchors isa AbstractVector{<:Pair} || error("reference must be a vector of " *
        "(concentrations => turnover) pairs; got $(repr(anchors))")
    isempty(anchors) &&
        error("reference is empty; pass at least one (concentrations => turnover) anchor")
    for (concs, measured) in anchors
        concs isa NamedTuple ||
            error("reference: concentrations must be a NamedTuple; got $(repr(concs))")
        for (k, c) in pairs(concs)
            k in mets || error("reference: $k is not a metabolite of the $owner " *
                               "(metabolites: $(join(mets, ", ")))")
            c isa Real && isfinite(c) && c >= 0 || error("reference: concentration of " *
                "$k must be a finite number ≥ 0; got $(repr(c))")
        end
        measured isa Real && isfinite(measured) && measured != 0 || error(
            "reference: turnover must be a finite, nonzero number; got $(repr(measured))")
    end
end

