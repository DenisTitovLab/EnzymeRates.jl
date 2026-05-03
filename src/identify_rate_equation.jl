# ABOUTME: Beam-search pipeline to identify the best rate equation.
# ABOUTME: Enumerates mechanisms, fits each, selects via CV.

using DataFrames
using CSV
using Statistics

"""
    IdentifyRateEquationProblem{R, D}

Holds the reaction, experimental data, and equilibrium
constant for rate equation identification.

# Fields
- `reaction`: the `EnzymeReaction` instance
- `data`: `NamedTuple` of column vectors with `:group`,
  `:Rate`, and metabolite columns
- `Keq`: fixed equilibrium constant
"""
struct IdentifyRateEquationProblem{
    R<:EnzymeReaction, D<:NamedTuple
}
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
        regulators(reaction)...,
    )
    for m in mnames
        m in col_names ||
            error(
                "Missing metabolite column: $m")
    end

    # Validate non-zero rates
    for i in eachindex(data.Rate)
        data.Rate[i] == 0 &&
            error(
                "Zero rate at row $i: " *
                "log(0) is undefined")
    end

    # Validate at least 2 groups for CV
    n_groups = length(unique(data.group))
    n_groups >= 2 || error(
        "Need at least 2 unique groups for " *
        "cross-validation, got $n_groups")

    IdentifyRateEquationProblem{
        typeof(reaction),typeof(data)
    }(reaction, data, Float64(Keq))
end

"""
    IdentifyRateEquationResults

Results from `identify_rate_equation`.

# Fields
- `best`: the best `AbstractEnzymeMechanism`
  (lowest loss at optimal param count)
- `cv_results`: `DataFrame` with LOOCV results for
  top candidates per param count
"""
struct IdentifyRateEquationResults
    best::AbstractEnzymeMechanism
    cv_results::DataFrame
end

"""
    identify_rate_equation(prob; kwargs...)

Find the best rate equation for the given reaction
and data using beam search.

# Keyword Arguments
- `min_beam_width::Int = 200`: minimum mechanisms
  to keep per level
- `beam_fraction::Float64 = 0.1`: fraction to keep
- `max_param_count::Int = 20`: stop expanding beyond
- `optimizer`: Optimization.jl optimizer (required).
  Recommended: `PyCMAOpt()` from OptimizationPyCMAES.
- `n_restarts::Int = 10`: multi-start restarts per fit
- `maxtime::Real = 60.0`: max time per fit (seconds)
- `maxiters::Int = 10_000_000`: max iterations per
  optimizer run (forwarded to `Optimization.solve`)
- `popsize::Int = 200`: population size for optimizer
  (forwarded to `Optimization.solve`)
- `verbose::Int = -9`: optimizer verbosity
  (forwarded to `Optimization.solve`)
- `n_cv_candidates::Int = 5`: LOOCV top N per
  param count
- `save_dir`: directory for per-level CSV files
- `pmap_function::Function = pmap`: parallelism
  function (Distributed.pmap by default)
- Extra kwargs are forwarded to `fit_rate_equation`
  and then to `Optimization.solve`.
"""
function identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    # Beam search
    min_beam_width::Int = 200,
    beam_fraction::Float64 = 0.1,
    max_param_count::Int = 20,
    # Fitting
    optimizer,
    n_restarts::Int = 10,
    maxtime::Real = 60.0,
    maxiters::Int = 10_000_000,
    popsize::Int = 200,
    verbose::Int = -9,
    # Model selection
    n_cv_candidates::Int = 5,
    # Output & parallelism
    save_dir::Union{Nothing,String} = nothing,
    pmap_function::Function = pmap,
    # Extra fitting/optimizer kwargs
    optim_kwargs...
)
    fitting_kwargs = (;
        n_restarts, maxtime,
        maxiters, popsize, verbose,
        optim_kwargs...)

    if save_dir !== nothing && isdir(save_dir)
        existing = filter(
            f -> endswith(f, ".csv"),
            readdir(save_dir))
        isempty(existing) || error(
            "save_dir already contains CSV " *
            "files. Use an empty directory " *
            "to avoid mixing results.")
    end

    specs, df = _beam_search(prob;
        min_beam_width, beam_fraction,
        max_param_count, save_dir,
        pmap_function, optimizer,
        fitting_kwargs...)

    return _cv_model_selection(
        specs, df, prob;
        n_cv_candidates, pmap_function,
        optimizer, fitting_kwargs...)
end

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
        rate_equation = rate_equation_string(
            mechanism),
        fitted_param_names = pnames,
        fitted_param_values = Tuple(
            fit_result.params[p]
            for p in pnames),
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

"""
Save results for one beam level to a CSV file.
"""
function _save_level_csv(
    save_dir::String, rows, n_fit_params::Int
)
    isdir(save_dir) || mkpath(save_dir)
    path = joinpath(
        save_dir, "params_$(n_fit_params).csv")
    df = _rows_to_dataframe(rows)
    CSV.write(path, df; append=isfile(path))
end

function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, beam_fraction,
    max_param_count, save_dir, pmap_function,
    optimizer, kwargs...
)
    # Initialize cache by param count
    cache = Dict{
        Int,Vector{AbstractMechanismSpec}
    }()
    for spec in init_mechanisms(prob.reaction)
        push!(
            get!(cache, spec.n_fit_params_estimate,
                AbstractMechanismSpec[]),
            spec)
    end
    dedup!(cache)

    all_specs = AbstractMechanismSpec[]
    all_rows = NamedTuple[]

    isempty(cache) && return (
        all_specs,
        _rows_to_dataframe(all_rows))

    min_pc = minimum(keys(cache))
    for pc in min_pc:max_param_count
        level = pop!(
            cache, pc,
            AbstractMechanismSpec[])
        isempty(level) &&
            (isempty(cache) ? break : continue)

        # Fit all specs at this level in parallel
        results = pmap_function(level) do spec
            try
                m = compile_mechanism(spec)
                fp = FittingProblem(
                    m, prob.data;
                    Keq=prob.Keq)
                fit = fit_rate_equation(
                    fp, optimizer; kwargs...)
                (spec=spec,
                 row=_build_result_row(m, fit),
                 ok=true)
            catch e
                @debug(
                    "Mechanism compilation/" *
                    "fitting failed",
                    exception=(
                        e, catch_backtrace()))
                (spec=spec,
                 row=nothing, ok=false)
            end
        end
        filter!(r -> r.ok, results)
        isempty(results) && continue

        append!(
            all_specs,
            [r.spec for r in results])
        append!(
            all_rows,
            [r.row for r in results])

        # Save CSV for this param count
        if save_dir !== nothing
            _save_level_csv(
                save_dir,
                [r.row for r in results], pc)
        end

        # Beam select within this level
        perm = sortperm(
            [r.row.loss for r in results])
        beam_size = max(
            ceil(Int,
                beam_fraction *
                length(results)),
            min_beam_width)
        beam_size = min(
            beam_size, length(results))
        beam_specs = [results[perm[i]].spec
                      for i in 1:beam_size]

        # Expand beam to next levels
        new_cache = expand_mechanisms(
            beam_specs, prob.reaction)
        for (target_pc, specs) in new_cache
            target_pc > max_param_count &&
                continue
            append!(
                get!(cache, target_pc,
                    AbstractMechanismSpec[]),
                specs)
        end
        dedup!(cache)
    end

    df = _rows_to_dataframe(all_rows)
    return all_specs, df
end

"""
Subset a columnar NamedTuple by a boolean mask.
"""
function _subset_data(data::NamedTuple, mask)
    return map(col -> col[mask], data)
end

"""
Evaluate loss of a mechanism on data with given
params.
"""
function _evaluate_loss(
    mechanism, data, params, Keq
)
    pnames = fitted_params(mechanism)
    x = [log(params[p]) for p in pnames]
    fp = FittingProblem(mechanism, data; Keq=Keq)
    return loss!(x, fp)
end

"""
    _loocv(mechanism, prob; optimizer, kwargs...)

Leave-one-group-out cross-validation. `optimizer`
is extracted explicitly because `fit_rate_equation`
takes it as a positional argument.
"""
function _loocv(
    mechanism::AbstractEnzymeMechanism,
    prob::IdentifyRateEquationProblem;
    optimizer, kwargs...
)
    groups = unique(prob.data.group)
    scores = Float64[]

    try
        for held_out in groups
            train_mask =
                prob.data.group .!= held_out
            test_mask =
                prob.data.group .== held_out

            train_data = _subset_data(
                prob.data, train_mask)
            test_data = _subset_data(
                prob.data, test_mask)

            fp_train = FittingProblem(
                mechanism, train_data;
                Keq=prob.Keq)
            fit = fit_rate_equation(
                fp_train, optimizer;
                kwargs...)

            test_loss = _evaluate_loss(
                mechanism, test_data,
                fit.params, prob.Keq)
            push!(scores, test_loss)
        end
    catch e
        @debug("LOOCV failed",
            exception=(e, catch_backtrace()))
        return Inf
    end

    result = mean(scores)
    return isfinite(result) ? result : Inf
end

function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function,
    optimizer, kwargs...
)
    # specs[i] corresponds to df[i, :] (same append
    # order in _beam_search, df is NOT sorted)
    isempty(specs) && error(
        "No mechanisms were successfully " *
        "fitted during beam search")
    df_indexed = copy(df)
    df_indexed.spec_idx = 1:nrow(df_indexed)

    # Group by n_params, take top n_cv_candidates
    # by loss
    candidate_indices = Int[]
    for gdf in groupby(df_indexed, :n_params)
        sorted = sort(gdf, :loss)
        n_take = min(
            n_cv_candidates, nrow(sorted))
        for i in 1:n_take
            push!(candidate_indices,
                sorted[i, :spec_idx])
        end
    end

    candidate_specs = specs[candidate_indices]
    candidate_rows = df[candidate_indices, :]

    # LOOCV each candidate in parallel
    cv_scores = pmap_function(
        candidate_specs
    ) do spec
        m = compile_mechanism(spec)
        _loocv(m, prob; optimizer, kwargs...)
    end

    # Build CV results DataFrame
    cv_df = copy(candidate_rows)
    cv_df.cv_score = collect(cv_scores)
    cv_df.spec_idx = candidate_indices

    # Best param count by CV score
    best_cv_per_pc = combine(
        groupby(cv_df, :n_params),
        :cv_score => minimum => :best_cv)
    best_pc_row = best_cv_per_pc[
        argmin(best_cv_per_pc.best_cv), :]
    best_param_count = best_pc_row.n_params

    # Best mechanism = lowest loss at that
    # param count
    at_best_pc = filter(
        row -> row.n_params == best_param_count,
        cv_df)
    sort!(at_best_pc, :loss)
    best_idx = at_best_pc[1, :spec_idx]
    best_mechanism = compile_mechanism(
        specs[best_idx])

    select!(cv_df, Not(:spec_idx))
    return IdentifyRateEquationResults(
        best_mechanism, cv_df)
end
