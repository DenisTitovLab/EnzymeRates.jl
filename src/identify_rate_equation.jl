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

"""
Build a result row NamedTuple from a fitted mechanism.
"""
function _build_result_row(mechanism, fit_result)
    pnames = fitted_params(mechanism)
    return (
        n_params = length(pnames),
        loss = fit_result.loss,
        mechanism_type = _mechanism_type_string(
            mechanism),
        rate_equation = rate_equation_string(mechanism),
        fitted_param_names = pnames,
        fitted_param_values = Tuple(
            fit_result.params[p] for p in pnames),
    )
end

"""
Convert a mechanism to an eval-able type string.
"""
function _mechanism_type_string(
    m::AbstractEnzymeMechanism
)
    return string(typeof(m))
end

"""
Save results for one beam level to a CSV file.
"""
function _save_level_csv(
    save_dir::String, rows, param_count::Int
)
    isdir(save_dir) || mkpath(save_dir)
    path = joinpath(
        save_dir, "params_$(param_count).csv")

    all_pnames = Set{Symbol}()
    for row in rows
        for p in row.fitted_param_names
            push!(all_pnames, p)
        end
    end
    sorted_pnames = sort(collect(all_pnames))

    df = DataFrame(
        n_params = [r.n_params for r in rows],
        loss = [r.loss for r in rows],
        mechanism_type = [
            r.mechanism_type for r in rows],
        rate_equation = [
            r.rate_equation for r in rows],
    )
    for pn in sorted_pnames
        df[!, pn] = [
            pn in r.fitted_param_names ?
                r.fitted_param_values[
                    findfirst(
                        ==(pn),
                        r.fitted_param_names)
                ] : missing
            for r in rows
        ]
    end

    CSV.write(path, df; append=isfile(path))
end

"""
Convert result row NamedTuples to a DataFrame.
Row order is preserved (no sorting) to maintain
alignment with the specs vector.
"""
function _rows_to_dataframe(rows)
    isempty(rows) && return DataFrame()

    all_pnames = Set{Symbol}()
    for row in rows
        for p in row.fitted_param_names
            push!(all_pnames, p)
        end
    end
    sorted_pnames = sort(collect(all_pnames))

    df = DataFrame(
        n_params = [r.n_params for r in rows],
        loss = [r.loss for r in rows],
        mechanism_type = [
            r.mechanism_type for r in rows],
        rate_equation = [
            r.rate_equation for r in rows],
    )
    for pn in sorted_pnames
        df[!, pn] = [
            pn in r.fitted_param_names ?
                r.fitted_param_values[
                    findfirst(
                        ==(pn),
                        r.fitted_param_names)
                ] : missing
            for r in rows
        ]
    end

    return df
end

function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width=200, beam_fraction=0.1,
    max_param_count=20, save_dir=nothing,
    pmap_function=map, optimizer=nothing, kwargs...
)
    specs = init_mechanisms(prob.reaction)
    all_specs = AbstractMechanismSpec[]
    all_rows = NamedTuple[]

    while !isempty(specs)
        filter!(
            s -> s.param_count <= max_param_count,
            specs)
        isempty(specs) && break

        # Fit all specs in parallel; catch failures
        results = pmap_function(specs) do spec
            try
                m = _compile(spec)
                fp = FittingProblem(
                    m, prob.data; Keq=prob.Keq)
                fit = fit_rate_equation(
                    fp, optimizer; kwargs...)
                (spec=spec,
                 row=_build_result_row(m, fit),
                 ok=true)
            catch e
                @warn(
                    "Mechanism compilation/fitting" *
                    " failed",
                    exception=(e, catch_backtrace()))
                (spec=spec, row=nothing, ok=false)
            end
        end
        filter!(r -> r.ok, results)
        isempty(results) && break

        new_specs = [r.spec for r in results]
        new_rows = [r.row for r in results]

        append!(all_specs, new_specs)
        append!(all_rows, new_rows)

        # Save per-param-count CSV if requested
        if save_dir !== nothing
            by_pc = Dict{Int,Vector{eltype(
                new_rows)}}()
            for row in new_rows
                push!(
                    get!(by_pc, row.n_params,
                        eltype(new_rows)[]),
                    row)
            end
            for (pc, rows) in by_pc
                _save_level_csv(save_dir, rows, pc)
            end
        end

        # Beam select: keep top mechanisms by loss
        perm = sortperm(
            [r.row.loss for r in results])
        beam_size = max(
            ceil(Int,
                beam_fraction * length(results)),
            min_beam_width)
        beam_size = min(beam_size, length(results))
        beam_specs = [results[perm[i]].spec
                      for i in 1:beam_size]

        # Expand to next level
        cache = expand_mechanisms(
            beam_specs, prob.reaction)
        dedup!(cache)

        specs = reduce(vcat,
            [v for (k, v) in cache
             if k <= max_param_count];
            init=AbstractMechanismSpec[])
    end

    df = _rows_to_dataframe(all_rows)
    return all_specs, df
end

function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates=5, pmap_function=map,
    optimizer=nothing, kwargs...
)
    error("Not implemented yet")
end
