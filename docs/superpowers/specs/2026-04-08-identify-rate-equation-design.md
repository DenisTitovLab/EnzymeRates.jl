# Design: `identify_rate_equation` — Beam Search Pipeline

## Goal

Implement the central `identify_rate_equation` function that takes an `EnzymeReaction` + kinetic data and finds the best-fitting rate equation with optimal complexity. Uses beam search over the mechanism enumeration space with cross-validation for model selection.

## Types

### `IdentifyRateEquationProblem`

```julia
struct IdentifyRateEquationProblem{R<:EnzymeReaction, D<:NamedTuple}
    reaction::R
    data::D          # columnar: :group, :Rate, + metabolite columns
    Keq::Float64
end
```

Constructor validates:
- `data` has `:group`, `:Rate` columns + one column per metabolite in `reaction`
- At least 2 unique `group` values (required for LOOCV)
- All `Rate` values are non-zero

### `IdentifyRateEquationResults`

```julia
struct IdentifyRateEquationResults
    best::AbstractEnzymeMechanism
    cv_results::DataFrame   # top N per param count with cv_score column
end
```

The `best` mechanism is selected as: the mechanism with lowest training loss at the param count that has the best (lowest) CV score.

`cv_results` DataFrame columns: `n_params`, `loss`, `cv_score`, `mechanism_type`, `rate_equation`, plus one column per fitted parameter.

## Public API

### `identify_rate_equation`

```julia
function identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    # Beam search
    max_beam_width::Int = 200,
    beam_fraction::Float64 = 0.1,
    max_param_count::Int = 20,
    # Fitting
    optimizer = PyCMAOpt(),
    n_restarts::Int = 10,
    maxtime::Real = 60.0,
    maxiters::Int = 10_000_000,
    popsize::Int = 200,
    # Model selection
    n_cv_candidates::Int = 5,
    # Output
    save_dir::Union{Nothing,String} = nothing,
    # Parallelism
    pmap_function::Function = map,
) → IdentifyRateEquationResults
```

## Internal Structure

Two internal functions. `identify_rate_equation` collects all its keyword arguments and passes them through as `kwargs...` to both internal functions. Each internal function destructures the kwargs it needs and passes the rest along (e.g., to `fit_rate_equation`). This avoids fragile per-kwarg copying.

### `_beam_search`

```julia
function _beam_search(
    prob::IdentifyRateEquationProblem;
    max_beam_width, beam_fraction, max_param_count,
    save_dir, pmap_function, optimizer, kwargs...
) → Vector{AbstractMechanismSpec}, DataFrame
```

Returns: the specs of all fitted mechanisms (needed for LOOCV recompilation) and a DataFrame of all fitting results.

Algorithm:

```
specs = init_mechanisms(prob.reaction)
all_results = []   # accumulates (spec, result_row) pairs

while !isempty(specs)
    filter!(s -> s.param_count <= max_param_count, specs)
    isempty(specs) && break

    # Fit all specs in parallel
    fitted = pmap_function(specs) do spec
        m = _compile(spec)
        fp = FittingProblem(m, prob.data; Keq=prob.Keq)
        fit = fit_rate_equation(fp, optimizer; kwargs...)  # n_restarts, maxtime, etc. flow through
        (spec=spec, row=_build_result_row(spec, m, fit))
    end

    append!(all_results, fitted)

    # Save per-param-count CSV files if save_dir set
    if save_dir !== nothing
        _save_level_results(save_dir, fitted)
    end

    # Beam select: keep top max(beam_fraction * n, max_beam_width) by loss
    sort!(fitted, by = r -> r.row.loss)
    beam_size = max(ceil(Int, beam_fraction * length(fitted)), max_beam_width)
    beam = fitted[1:min(beam_size, end)]
    beam_specs = [r.spec for r in beam]

    # Expand to next level
    cache = expand_mechanisms(beam_specs, prob.reaction)
    dedup!(cache)
    specs = reduce(vcat, [v for (k,v) in cache if k <= max_param_count]; init=eltype(values(cache))[])
end

return [r.spec for r in all_results], DataFrame([r.row for r in all_results])
```

### `_cv_model_selection`

```julia
function _cv_model_selection(
    specs::Vector{<:AbstractMechanismSpec},
    df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function, kwargs...
) → IdentifyRateEquationResults
```

Algorithm:

```
# Group by n_params, take top n_cv_candidates per group by loss
candidates = []
for group in groupby(df, :n_params)
    sorted = sort(group, :loss)
    top = sorted[1:min(n_cv_candidates, nrow(sorted)), :]
    append!(candidates, zip(matching_specs, eachrow(top)))
end

# LOOCV each candidate in parallel
groups = unique(prob.data.group)
cv_results = pmap_function(candidates) do (spec, row)
    m = _compile(spec)
    cv_score = _loocv(m, prob, groups; kwargs...)
    merge(row, (cv_score=cv_score,))
end

cv_df = DataFrame(cv_results)

# Select best: param count with best CV, then lowest loss at that count
best_param_count = cv_df[argmin(cv_df.cv_score), :n_params]
best_row = filter(r -> r.n_params == best_param_count, cv_df) |> 
           r -> sort(r, :loss) |> first
best_spec = ... # matching spec
best_mechanism = _compile(best_spec)

return IdentifyRateEquationResults(best_mechanism, cv_df)
```

### `_loocv`

```julia
function _loocv(
    mechanism::AbstractEnzymeMechanism,
    prob::IdentifyRateEquationProblem,
    groups; kwargs...
) → Float64
```

Leave-one-group-out cross-validation:

```
scores = Float64[]
for held_out in groups
    train_mask = prob.data.group .!= held_out
    test_mask = prob.data.group .== held_out
    train_data = _subset_data(prob.data, train_mask)
    test_data = _subset_data(prob.data, test_mask)

    fp_train = FittingProblem(mechanism, train_data; Keq=prob.Keq)
    fit = fit_rate_equation(fp_train, optimizer; kwargs...)

    # Evaluate on held-out group
    test_loss = _evaluate_loss(mechanism, test_data, fit.params, prob.Keq)
    push!(scores, test_loss)
end
return mean(scores)
```

### `_compile`

Dispatches to the right constructor:

```julia
_compile(spec::MechanismSpec) = EnzymeMechanism(spec)
_compile(spec::AllostericMechanismSpec) = AllostericEnzymeMechanism(spec)
```

## CSV Output

When `save_dir` is provided, one CSV file per param count level:

```
save_dir/
  params_7.csv
  params_8.csv
  params_9.csv
  ...
```

Each CSV has columns: `n_params, loss, mechanism_type, rate_equation, K_S, K_P, k1f, ...`

Written after each beam level completes. Parameter columns vary per mechanism — missing params get `missing` values.

## FittingProblem Migration

Migrate `FittingProblem` to use a `group` column instead of `Article` + `Fig`:

- Replace `fig_point_indexes` with `group_point_indexes`
- Grouping becomes `unique(data.group)` instead of `unique(zip(data.Article, data.Fig))`
- Update constructor validation to require `:group` instead of `:Article` and `:Fig`
- Update all existing tests to use `:group` column

## Parallelism

The `pmap_function` parameter accepts any function with `map` semantics:

- Default: `map` (serial execution)
- Single machine: `ThreadsX.map` or similar
- HPC cluster: `Distributed.pmap`

The parallelism unit is "fit one mechanism" — each call to the map function receives one spec, compiles it, fits it, and returns a result row. No shared mutable state between map calls.

For LOOCV in Phase 2, the parallelism unit is "LOOCV one candidate" — each candidate's full leave-one-group-out loop runs on one worker.

## Dependencies

- `DataFrames.jl` — for result DataFrame
- `CSV.jl` — for incremental saves
- `Statistics.jl` — for `mean` in LOOCV

These should be added to Project.toml.

## Testing Strategy

### Unit tests for `_beam_search`
- Use a simple uni-uni reaction with synthetic data from a known mechanism
- Verify that the known mechanism appears in results
- Verify beam selection keeps best mechanisms

### Unit tests for `_loocv`
- Verify that LOOCV score is computed correctly for a known mechanism
- Verify with 2-3 groups

### Integration test for `identify_rate_equation`
- Uni-uni reaction with synthetic data generated from a specific mechanism
- Verify that `results.best` recovers the correct mechanism type
- Verify that `results.cv_results` has cv_score for top candidates
- Verify CSV files are written when `save_dir` is provided

### FittingProblem migration tests
- Existing test_fitting.jl tests updated to use `:group` column
- Verify identical behavior

## Exports

Add to `src/EnzymeRates.jl`:
```julia
export IdentifyRateEquationProblem, IdentifyRateEquationResults
export identify_rate_equation
```

## File Layout

- `src/identify_rate_equation.jl` — `IdentifyRateEquationProblem`, `IdentifyRateEquationResults`, `identify_rate_equation`, `_beam_search`, `_cv_model_selection`, `_loocv`, `_compile`, helpers
- `test/test_identify_rate_equation.jl` — all tests for the new functionality
