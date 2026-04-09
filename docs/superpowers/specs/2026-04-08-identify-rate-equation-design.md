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
- `data` has `:group`, `:Rate` columns + one column per substrate/product name in `reaction`
  (metabolite names extracted via `substrates(reaction)` and `products(reaction)`, which return
  `(name, atoms)` pairs — use `s[1]` to get the name `Symbol`)
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
    min_beam_width::Int = 200,
    beam_fraction::Float64 = 0.1,
    max_param_count::Int = 20,
    # Fitting — optimizer is passed positionally to fit_rate_equation
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

Two internal functions. `identify_rate_equation` collects all its keyword arguments and passes them through as `kwargs...` to both internal functions. Each internal function destructures the kwargs it needs and passes the rest along (e.g., to `fit_rate_equation`). Note: `optimizer` is extracted explicitly because `fit_rate_equation` takes it as a positional argument — it cannot flow through `kwargs...`.

### `_beam_search`

```julia
function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, beam_fraction, max_param_count,
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

    # Fit all specs in parallel; wrap in try/catch to handle compilation failures
    fitted = pmap_function(specs) do spec
        try
            m = _compile(spec)
            fp = FittingProblem(m, prob.data; Keq=prob.Keq)
            fit = fit_rate_equation(fp, optimizer; kwargs...)
            (spec=spec, row=_build_result_row(spec, m, fit), ok=true)
        catch e
            @warn "Mechanism compilation/fitting failed" exception=e
            (spec=spec, row=nothing, ok=false)
        end
    end
    filter!(r -> r.ok, fitted)

    append!(all_results, fitted)

    # Save per-param-count CSV files if save_dir set
    if save_dir !== nothing
        _save_level_results(save_dir, fitted)
    end

    # Beam select: keep at least min_beam_width, or beam_fraction of total
    sort!(fitted, by = r -> r.row.loss)
    beam_size = max(ceil(Int, beam_fraction * length(fitted)), min_beam_width)
    beam_size = min(beam_size, length(fitted))
    beam = fitted[1:beam_size]
    beam_specs = [r.spec for r in beam]

    # Expand to next level
    cache = expand_mechanisms(beam_specs, prob.reaction)
    dedup!(cache)
    specs = reduce(vcat, [v for (k,v) in cache if k <= max_param_count]; init=AbstractMechanismSpec[])
end

return [r.spec for r in all_results], DataFrame([r.row for r in all_results])
```

### `_cv_model_selection`

```julia
function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function, optimizer, kwargs...
) → IdentifyRateEquationResults
```

Algorithm:

Note: `specs` and `df` rows are positionally aligned (both built in the same
order by `_beam_search`). We track the original index via a `spec_idx` column
to map back from sorted/grouped DataFrames to the correct spec.

```
# Track original row→spec mapping (df may be unsorted)
df_indexed = copy(df)
df_indexed.spec_idx = 1:nrow(df_indexed)

# Group by n_params, take top n_cv_candidates per group by loss
candidate_indices = Int[]
for gdf in groupby(df_indexed, :n_params)
    sorted = sort(gdf, :loss)
    for i in 1:min(n_cv_candidates, nrow(sorted))
        push!(candidate_indices, sorted[i, :spec_idx])
    end
end

# LOOCV each candidate in parallel
candidate_specs = specs[candidate_indices]
cv_scores = pmap_function(candidate_specs) do spec
    m = _compile(spec)
    _loocv(m, prob; optimizer, kwargs...)
end

# Build CV results, select best param count, then best loss at that count
cv_df = df[candidate_indices, :]
cv_df.cv_score = collect(cv_scores)
...
best_mechanism = _compile(specs[best_idx])
return IdentifyRateEquationResults(best_mechanism, cv_df)
```

### `_loocv`

```julia
function _loocv(
    mechanism::AbstractEnzymeMechanism,
    prob::IdentifyRateEquationProblem;
    optimizer, kwargs...
) → Float64
```

Leave-one-group-out cross-validation. `optimizer` is extracted explicitly
from kwargs because `fit_rate_equation` takes it as a positional argument.

```
groups = unique(prob.data.group)
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

If the same param count appears across multiple beam iterations, results are **appended** to the existing file (not overwritten).

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
