# ABOUTME: Beam-search pipeline to identify the best rate equation for an enzyme.
# ABOUTME: Enumerates mechanisms, fits each to data, selects optimal complexity via CV.

using DataFrames
using CSV
using Statistics

"""
    IdentifyRateEquationProblem{R, D}

Holds the reaction, experimental data, and equilibrium constant for rate equation
identification.

# Fields
- `reaction`: the `EnzymeReaction` instance
- `data`: `NamedTuple` of column vectors with `:group`, `:Rate`, and metabolite columns
- `Keq`: fixed equilibrium constant
"""
struct IdentifyRateEquationProblem{R<:EnzymeReaction, D<:NamedTuple}
    reaction::R
    data::D
    Keq::Float64
end

function IdentifyRateEquationProblem(
    reaction::EnzymeReaction, table; Keq::Real
)
    data = Tables.columntable(table)
    col_names = keys(data)

    for req in (:group, :Rate)
        req in col_names ||
            error("Missing required column: $req")
    end
    # Extract metabolite names from reaction
    # (substrates/products are (name, atoms) pairs)
    mnames = tuple(
        [s[1] for s in substrates(reaction)]...,
        [p[1] for p in products(reaction)]...,
    )
    for m in mnames
        m in col_names ||
            error("Missing metabolite column: $m")
    end

    # Validate non-zero rates
    for i in eachindex(data.Rate)
        data.Rate[i] == 0 &&
            error("Zero rate at row $i: log(0) is undefined")
    end

    # Validate at least 2 groups for CV
    n_groups = length(unique(data.group))
    n_groups >= 2 || error(
        "Need at least 2 unique groups for " *
        "cross-validation, got $n_groups")

    IdentifyRateEquationProblem{typeof(reaction),typeof(data)}(
        reaction, data, Float64(Keq))
end

"""
    IdentifyRateEquationResults

Results from `identify_rate_equation`.

# Fields
- `best`: the best `AbstractEnzymeMechanism`
  (lowest loss at optimal param count)
- `cv_results`: `DataFrame` with LOOCV results for top
  candidates per param count
"""
struct IdentifyRateEquationResults
    best::AbstractEnzymeMechanism
    cv_results::DataFrame
end

"""
    identify_rate_equation(prob; kwargs...)

Find the best rate equation for the given reaction and data
using beam search.

# Keyword Arguments
- `min_beam_width::Int = 200`: minimum mechanisms to keep
- `beam_fraction::Float64 = 0.1`: fraction to keep
- `max_param_count::Int = 20`: stop expanding beyond this
- `optimizer`: Optimization.jl optimizer
- `n_restarts::Int = 10`: multi-start restarts per fit
- `maxtime::Real = 60.0`: max time per fit (seconds)
- `n_cv_candidates::Int = 5`: LOOCV top N per param count
- `save_dir`: directory for per-level CSV files
- `pmap_function::Function = map`: parallelism function
"""
function identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    min_beam_width::Int = 200,
    beam_fraction::Float64 = 0.1,
    max_param_count::Int = 20,
    optimizer = nothing,
    n_cv_candidates::Int = 5,
    save_dir::Union{Nothing,String} = nothing,
    pmap_function::Function = map,
    kwargs...
)
    specs, df = _beam_search(prob;
        min_beam_width, beam_fraction, max_param_count,
        save_dir, pmap_function, optimizer, kwargs...)

    return _cv_model_selection(specs, df, prob;
        n_cv_candidates, pmap_function, optimizer,
        kwargs...)
end

# Compile a mechanism spec to the appropriate type.
_compile(spec::MechanismSpec) = EnzymeMechanism(spec)
_compile(spec::AllostericMechanismSpec) =
    AllostericEnzymeMechanism(spec)

function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width=200, beam_fraction=0.1,
    max_param_count=20, save_dir=nothing,
    pmap_function=map, optimizer=nothing, kwargs...
)
    error("Not implemented yet")
end

function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates=5, pmap_function=map,
    optimizer=nothing, kwargs...
)
    error("Not implemented yet")
end
