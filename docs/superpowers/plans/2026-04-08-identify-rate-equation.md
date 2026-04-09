# identify_rate_equation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `identify_rate_equation()` — a beam-search pipeline that enumerates enzyme mechanisms, fits each to data, and selects the optimal rate equation via cross-validation.

**Architecture:** Two-phase pipeline. Phase 1 (`_beam_search`) iteratively expands and fits mechanisms, keeping the top N by training loss. Phase 2 (`_cv_model_selection`) runs leave-one-group-out CV on the top candidates per param count to select optimal complexity. Parallelism is user-pluggable via a `pmap_function` argument.

**Tech Stack:** Julia, DataFrames.jl, CSV.jl, Statistics.jl, Optimization.jl (PyCMAOpt default)

**Spec:** `docs/superpowers/specs/2026-04-08-identify-rate-equation-design.md`

---

### Task 1: Migrate FittingProblem from Article/Fig to group column

**Files:**
- Modify: `src/fitting.jl`
- Modify: `test/test_fitting.jl`

This task changes `FittingProblem` to use a single `group` column instead of `Article` + `Fig` for grouping data points. Also widens the type parameter from `EnzymeMechanism` to `AbstractEnzymeMechanism` so it works with allosteric mechanisms.

- [ ] **Step 1: Update test_fitting.jl — change make_synthetic_data helper**

In `test/test_fitting.jl`, replace the `make_synthetic_data` helper function (lines 19-39) to use `group` instead of `articles`/`figs`:

```julia
    function make_synthetic_data(
            mechanism, true_params, concs_list;
            groups=fill("G1", length(concs_list)),
            scale=1.0,
    )
        rates = Float64[]
        for (i, concs) in enumerate(concs_list)
            r = rate_equation(mechanism, concs, true_params) * scale
            push!(rates, r)
        end
        met_names = metabolites(mechanism)
        cols = Dict{Symbol, Vector}()
        cols[:group] = groups
        cols[:Rate] = rates
        for mn in met_names
            cols[mn] = [concs[mn] for concs in concs_list]
        end
        return (; (k => cols[k] for k in (:group, :Rate, met_names...))...)
    end
```

- [ ] **Step 2: Update all test data to use group instead of Article/Fig**

In `test/test_fitting.jl`, make the following replacements throughout the file:

Test "Construction" (lines 51-68): no change needed — uses default `make_synthetic_data` args.

Test "Centering invariance" (lines 92-120): replace `articles=fill("A1", 5), figs=fill("F1", 5)` with `groups=fill("G1", 5)` in both calls to `make_synthetic_data`.

Test "Multi-figure centering invariance" (lines 122-154): replace the two groups test. Change:
```julia
        data1 = make_synthetic_data(uni_uni, true_params, concs_list;
            articles=["A1","A1","A1","A2","A2"],
            figs=["F1","F1","F1","F1","F1"],
            scale=1.0)
```
to:
```julia
        data1 = make_synthetic_data(uni_uni, true_params, concs_list;
            groups=["G1","G1","G1","G2","G2"],
            scale=1.0)
```

Test "Sign mismatch penalty" (lines 156-183): change the inline data:
```julia
        data = (
            group = ["G1", "G1", "G1"],
            Rate = [1.0, 2.0, 3.0],
            S = [1.0, 2.0, 3.0],
            P = [0.1, 0.1, 0.1],
        )
```

Test "All-mismatch figure" (lines 195-245): change all three sub-test inline data blocks to use `group = fill("G1", 5)` instead of `Article = fill("A1", 5), Fig = fill("F1", 5)`.

Test "Speed" (lines 262-310): replace:
```julia
        articles = [string("A", div(i-1, 50)+1) for i in 1:n_points]
        figs = [string("F", mod(i-1, 5)+1) for i in 1:n_points]
```
with:
```julia
        groups = [string("G", div(i-1, 50)+1) for i in 1:n_points]
```
And update the data tuple to use `group = groups` instead of `Article = articles, Fig = figs`.

Test "Validation errors" (lines 344-361): update all inline data tuples to use `group` instead of `Article`/`Fig`. For example:
```julia
        data_no_rate = (group = ["G1"], S = [1.0], P = [0.1])
```
Remove the "Missing Article column" test and replace with "Missing group column":
```julia
        data_no_group = (Rate = [1.0], S = [1.0], P = [0.1])
        @test_throws ErrorException FittingProblem(uni_uni, data_no_group; Keq=1.0)
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: Tests fail because `FittingProblem` still expects `Article`/`Fig` columns.

- [ ] **Step 4: Update FittingProblem struct and constructor in fitting.jl**

In `src/fitting.jl`, make these changes:

Change the struct (line 16):
```julia
struct FittingProblem{M<:AbstractEnzymeMechanism, D<:NamedTuple}
```

Update the docstring and field description (lines 3-15): replace `fig_point_indexes` with `group_point_indexes`, mention `group` column instead of `Article`/`Fig`.

Change the constructor signature (line 35):
```julia
function FittingProblem(mechanism::AbstractEnzymeMechanism, table; Keq::Real)
```

In the constructor body, replace the column validation (lines 42-44):
```julia
    for req in (:group, :Rate)
        req in col_names || error("Missing required column: $req")
    end
```

Replace the grouping logic (lines 59-71):
```julia
    # Build group_point_indexes by grouping on :group column
    group_col = data.group
    group_map = Dict{eltype(group_col), Vector{Int}}()
    for i in 1:n
        key = group_col[i]
        if haskey(group_map, key)
            push!(group_map[key], i)
        else
            group_map[key] = [i]
        end
    end
    group_point_indexes = collect(values(group_map))
```

Update the struct field name in the constructor call (line 76):
```julia
    FittingProblem{typeof(mechanism), typeof(data)}(
        mechanism, data, group_point_indexes, Float64(Keq),
        log_abs_rates, log_ratios_buffer
    )
```

Rename the field in the struct definition (line 19):
```julia
    group_point_indexes::Vector{Vector{Int}}
```

- [ ] **Step 5: Update loss! to use group_point_indexes**

In `src/fitting.jl`, line 124, change:
```julia
    @inbounds for fig_idx in fp.fig_point_indexes
```
to:
```julia
    @inbounds for fig_idx in fp.group_point_indexes
```

- [ ] **Step 6: Update test assertions for new field name**

In `test/test_fitting.jl`, Test "Construction" (line 66-67), change:
```julia
        @test length(fp.group_point_indexes) == 1
        @test length(fp.group_point_indexes[1]) == 5
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass.

- [ ] **Step 8: Commit**

```bash
git add src/fitting.jl test/test_fitting.jl
git commit -m "Migrate FittingProblem from Article/Fig to group column

Also widen type parameter from EnzymeMechanism to AbstractEnzymeMechanism."
```

---

### Task 2: Add DataFrames.jl, CSV.jl, Statistics.jl dependencies

**Files:**
- Modify: `Project.toml`

- [ ] **Step 1: Add dependencies**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using Pkg; Pkg.add(["DataFrames", "CSV", "Statistics"])'
```

- [ ] **Step 2: Verify the packages resolve**

```bash
julia --project -e 'using DataFrames, CSV, Statistics; println("OK")'
```
Expected: prints "OK"

- [ ] **Step 3: Commit**

```bash
git add Project.toml Manifest.toml
git commit -m "Add DataFrames, CSV, Statistics dependencies"
```

---

### Task 3: Create identify_rate_equation.jl with types and stubs

**Files:**
- Create: `src/identify_rate_equation.jl`
- Modify: `src/EnzymeRates.jl`

- [ ] **Step 1: Create the file with ABOUTME, types, and function stubs**

Create `src/identify_rate_equation.jl`:

```julia
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

function IdentifyRateEquationProblem(reaction::EnzymeReaction, table; Keq::Real)
    data = Tables.columntable(table)
    col_names = keys(data)

    for req in (:group, :Rate)
        req in col_names || error("Missing required column: $req")
    end
    # Extract metabolite names from reaction (substrates/products are (name, atoms) pairs)
    mnames = tuple([s[1] for s in substrates(reaction)]...,
                   [p[1] for p in products(reaction)]...)
    for m in mnames
        m in col_names || error("Missing metabolite column: $m")
    end

    # Validate non-zero rates
    for i in eachindex(data.Rate)
        data.Rate[i] == 0 && error("Zero rate at row $i: log(0) is undefined")
    end

    # Validate at least 2 groups for CV
    n_groups = length(unique(data.group))
    n_groups >= 2 || error("Need at least 2 unique groups for cross-validation, got $n_groups")

    IdentifyRateEquationProblem{typeof(reaction), typeof(data)}(
        reaction, data, Float64(Keq))
end

"""
    IdentifyRateEquationResults

Results from `identify_rate_equation`.

# Fields
- `best`: the best `AbstractEnzymeMechanism` (lowest loss at optimal param count)
- `cv_results`: `DataFrame` with LOOCV results for top candidates per param count
"""
struct IdentifyRateEquationResults
    best::AbstractEnzymeMechanism
    cv_results::DataFrame
end

"""
    identify_rate_equation(prob::IdentifyRateEquationProblem; kwargs...)

Find the best rate equation for the given reaction and data using beam search.

# Keyword Arguments
- `min_beam_width::Int = 200`: minimum number of mechanisms to keep per beam level
- `beam_fraction::Float64 = 0.1`: fraction of mechanisms to keep (beam = max of this and min_beam_width)
- `max_param_count::Int = 20`: stop expanding beyond this parameter count
- `optimizer`: Optimization.jl optimizer (default: PyCMAOpt())
- `n_restarts::Int = 10`: multi-start restarts per mechanism fit
- `maxtime::Real = 60.0`: max time per mechanism fit (seconds)
- `maxiters::Int = 10_000_000`: max iterations per optimizer run
- `popsize::Int = 200`: population size for population-based optimizers
- `n_cv_candidates::Int = 5`: number of top mechanisms per param count to LOOCV
- `save_dir::Union{Nothing,String} = nothing`: directory for per-level CSV files
- `pmap_function::Function = map`: parallelism function (e.g., `Distributed.pmap`)
"""
function identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    min_beam_width::Int = 200,
    beam_fraction::Float64 = 0.1,
    max_param_count::Int = 20,
    optimizer = BBO_adaptive_de_rand_1_bin_radiuslimited(),
    n_cv_candidates::Int = 5,
    save_dir::Union{Nothing,String} = nothing,
    pmap_function::Function = map,
    kwargs...
)
    specs, df = _beam_search(prob;
        min_beam_width, beam_fraction, max_param_count,
        save_dir, pmap_function, optimizer, kwargs...)

    return _cv_model_selection(specs, df, prob;
        n_cv_candidates, pmap_function, optimizer, kwargs...)
end
```

- [ ] **Step 2: Add include and exports to EnzymeRates.jl**

In `src/EnzymeRates.jl`, uncomment and update the export lines:
```julia
export IdentifyRateEquationProblem, IdentifyRateEquationResults
export identify_rate_equation
```

Add the include after `mechanism_enumeration.jl`:
```julia
include("identify_rate_equation.jl")
```

- [ ] **Step 3: Verify the module loads**

```bash
julia --project -e 'using EnzymeRates; println("OK")'
```
Expected: Fails because `_beam_search` and `_cv_model_selection` are not defined yet. That's fine — we'll add stubs next.

- [ ] **Step 4: Add function stubs so the module loads**

Append to `src/identify_rate_equation.jl`:

```julia
# Compile a mechanism spec to the appropriate mechanism type.
_compile(spec::MechanismSpec) = EnzymeMechanism(spec)
_compile(spec::AllostericMechanismSpec) = AllostericEnzymeMechanism(spec)

function _beam_search(prob::IdentifyRateEquationProblem;
    min_beam_width=200, beam_fraction=0.1, max_param_count=20,
    save_dir=nothing, pmap_function=map, optimizer=nothing, kwargs...)
    error("Not implemented yet")
end

function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates=5, pmap_function=map, optimizer=nothing, kwargs...)
    error("Not implemented yet")
end
```

- [ ] **Step 5: Verify the module loads now**

```bash
julia --project -e 'using EnzymeRates; println("OK")'
```
Expected: prints "OK"

- [ ] **Step 6: Commit**

```bash
git add src/identify_rate_equation.jl src/EnzymeRates.jl
git commit -m "Add IdentifyRateEquationProblem, IdentifyRateEquationResults types and stubs"
```

---

### Task 4: Implement _beam_search

**Files:**
- Modify: `src/identify_rate_equation.jl`
- Create: `test/test_identify_rate_equation.jl`
- Modify: `test/runtests.jl`

- [ ] **Step 1: Write test for _beam_search**

Create `test/test_identify_rate_equation.jl`:

```julia
using DataFrames
using CSV
using Statistics

@testset "identify_rate_equation" begin

    # ── Shared test setup ────────────────────────────────────────────────────
    # Simple Uni-Uni mechanism for testing
    test_mechanism = @enzyme_mechanism begin
        species: begin
            substrates: S[C]
            products:   P[C]
            enzymes:    E, ES[C]
        end
        steps: begin
            [E, S] <--> [ES]
            [ES] <--> [E, P]
        end
    end

    test_rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end

    Keq_val = 2.0
    true_params = (k1f = 10.0, k2f = 5.0, k2r = 1.0, Keq = Keq_val, E_total = 1.0)

    # Generate synthetic data with multiple groups
    function make_test_data(mechanism, params, rxn; n_per_group=5, n_groups=3)
        groups = String[]
        rates = Float64[]
        S_vals = Float64[]
        P_vals = Float64[]
        for g in 1:n_groups
            for _ in 1:n_per_group
                s = 0.1 + 9.9 * rand()
                p = 0.1 + 9.9 * rand()
                concs = (S = s, P = p)
                r = rate_equation(mechanism, concs, params) * (0.5 + rand())
                push!(groups, "G$g")
                push!(rates, r)
                push!(S_vals, s)
                push!(P_vals, p)
            end
        end
        return (group = groups, Rate = rates, S = S_vals, P = P_vals)
    end

    Random.seed!(42)
    test_data = make_test_data(test_mechanism, true_params, test_rxn)

    # ── Test: IdentifyRateEquationProblem construction ────────────────────────
    @testset "IdentifyRateEquationProblem construction" begin
        prob = IdentifyRateEquationProblem(test_rxn, test_data; Keq=Keq_val)
        @test prob.reaction === test_rxn
        @test prob.Keq == Keq_val
        @test length(unique(prob.data.group)) == 3

        # Validation: missing column
        bad_data = (group = ["G1"], Rate = [1.0])
        @test_throws ErrorException IdentifyRateEquationProblem(test_rxn, bad_data; Keq=1.0)

        # Validation: zero rate
        bad_data2 = (group = ["G1", "G2"], Rate = [0.0, 1.0], S = [1.0, 1.0], P = [0.1, 0.1])
        @test_throws ErrorException IdentifyRateEquationProblem(test_rxn, bad_data2; Keq=1.0)

        # Validation: need >= 2 groups
        bad_data3 = (group = ["G1", "G1"], Rate = [1.0, 2.0], S = [1.0, 2.0], P = [0.1, 0.1])
        @test_throws ErrorException IdentifyRateEquationProblem(test_rxn, bad_data3; Keq=1.0)
    end

    # ── Test: _beam_search returns results ────────────────────────────────────
    @testset "_beam_search" begin
        using OptimizationBBO
        prob = IdentifyRateEquationProblem(test_rxn, test_data; Keq=Keq_val)

        specs, df = EnzymeRates._beam_search(prob;
            min_beam_width=200, beam_fraction=0.1, max_param_count=8,
            save_dir=nothing, pmap_function=map,
            optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
            n_restarts=2, maxtime=5.0)

        # Should return non-empty results
        @test length(specs) > 0
        @test nrow(df) > 0
        @test nrow(df) == length(specs)

        # DataFrame should have required columns
        @test "n_params" in names(df)
        @test "loss" in names(df)
        @test "mechanism_type" in names(df)
        @test "rate_equation" in names(df)

        # All losses should be finite and non-negative
        @test all(isfinite, df.loss)
        @test all(>=(0), df.loss)

        # Should have mechanisms at multiple param counts
        @test length(unique(df.n_params)) >= 1
    end

end
```

- [ ] **Step 2: Add test file to runtests.jl**

In `test/runtests.jl`, add before the closing `end`:
```julia
    include("test_identify_rate_equation.jl")
```

- [ ] **Step 3: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: Fails with "Not implemented yet" from `_beam_search` stub.

- [ ] **Step 4: Implement _beam_search**

In `src/identify_rate_equation.jl`, replace the `_beam_search` stub with:

```julia
"""
Build a result row NamedTuple from a fitted mechanism.
"""
function _build_result_row(spec, mechanism, fit_result)
    pnames = fitted_params(mechanism)
    param_dict = Dict{Symbol,Float64}()
    for p in pnames
        param_dict[p] = fit_result.params[p]
    end
    return (
        n_params = length(pnames),
        loss = fit_result.loss,
        mechanism_type = _mechanism_type_string(mechanism),
        rate_equation = rate_equation_string(mechanism),
        fitted_param_names = pnames,
        fitted_param_values = Tuple(fit_result.params[p] for p in pnames),
    )
end

"""
Convert a mechanism to an eval-able type string for CSV storage.
"""
function _mechanism_type_string(m::AbstractEnzymeMechanism)
    return string(typeof(m))
end

"""
Save results for one beam level to a CSV file in save_dir.
"""
function _save_level_csv(save_dir::String, rows, param_count::Int)
    isdir(save_dir) || mkpath(save_dir)
    path = joinpath(save_dir, "params_$(param_count).csv")

    # Collect all parameter names across rows
    all_pnames = Set{Symbol}()
    for row in rows
        for p in row.fitted_param_names
            push!(all_pnames, p)
        end
    end
    sorted_pnames = sort(collect(all_pnames))

    # Build DataFrame
    df = DataFrame(
        n_params = [r.n_params for r in rows],
        loss = [r.loss for r in rows],
        mechanism_type = [r.mechanism_type for r in rows],
        rate_equation = [r.rate_equation for r in rows],
    )
    for pn in sorted_pnames
        df[!, pn] = [
            pn in r.fitted_param_names ?
                r.fitted_param_values[findfirst(==(pn), r.fitted_param_names)] :
                missing
            for r in rows
        ]
    end

    CSV.write(path, df; append=isfile(path))
end

function _beam_search(prob::IdentifyRateEquationProblem;
    min_beam_width, beam_fraction, max_param_count,
    save_dir, pmap_function, optimizer, kwargs...
)
    specs = init_mechanisms(prob.reaction)
    all_specs = AbstractMechanismSpec[]
    all_rows = NamedTuple[]

    while !isempty(specs)
        filter!(s -> s.param_count <= max_param_count, specs)
        isempty(specs) && break

        # Fit all specs in parallel; catch compilation/fitting failures
        results = pmap_function(specs) do spec
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
        filter!(r -> r.ok, results)
        isempty(results) && break

        new_specs = [r.spec for r in results]
        new_rows = [r.row for r in results]

        append!(all_specs, new_specs)
        append!(all_rows, new_rows)

        # Save per-param-count CSV if requested
        if save_dir !== nothing
            by_pc = Dict{Int, Vector{eltype(new_rows)}}()
            for row in new_rows
                push!(get!(by_pc, row.n_params, eltype(new_rows)[]), row)
            end
            for (pc, rows) in by_pc
                _save_level_csv(save_dir, rows, pc)
            end
        end

        # Beam select: keep at least min_beam_width, or beam_fraction of total
        perm = sortperm([r.row.loss for r in results])
        beam_size = max(ceil(Int, beam_fraction * length(results)), min_beam_width)
        beam_size = min(beam_size, length(results))
        beam_specs = [results[perm[i]].spec for i in 1:beam_size]

        # Expand to next level
        cache = expand_mechanisms(beam_specs, prob.reaction)
        dedup!(cache)

        specs = reduce(vcat,
            [v for (k, v) in cache if k <= max_param_count];
            init=AbstractMechanismSpec[])
    end

    # Build results DataFrame
    df = _rows_to_dataframe(all_rows)
    return all_specs, df
end

"""
Convert a vector of result row NamedTuples to a DataFrame.
"""
function _rows_to_dataframe(rows)
    isempty(rows) && return DataFrame()

    # Collect all parameter names
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
        mechanism_type = [r.mechanism_type for r in rows],
        rate_equation = [r.rate_equation for r in rows],
    )
    for pn in sorted_pnames
        df[!, pn] = [
            pn in r.fitted_param_names ?
                r.fitted_param_values[findfirst(==(pn), r.fitted_param_names)] :
                missing
            for r in rows
        ]
    end

    # Do NOT sort here — row order must match all_specs order for positional
    # indexing in _cv_model_selection.
    return df
end
```

- [ ] **Step 5: Run test to verify it passes**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: `_beam_search` test passes. This may take a while due to fitting.

- [ ] **Step 6: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl test/runtests.jl
git commit -m "Implement _beam_search for mechanism enumeration and fitting pipeline"
```

---

### Task 5: Implement _cv_model_selection and _loocv

**Files:**
- Modify: `src/identify_rate_equation.jl`
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write test for _loocv**

In `test/test_identify_rate_equation.jl`, add inside the main `@testset` block:

```julia
    # ── Test: _loocv ─────────────────────────────────────────────────────────
    @testset "_loocv" begin
        prob = IdentifyRateEquationProblem(test_rxn, test_data; Keq=Keq_val)

        # LOOCV with the true mechanism should give a finite score
        using OptimizationBBO
        cv_score = EnzymeRates._loocv(test_mechanism, prob;
            optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
            n_restarts=2, maxtime=5.0)

        @test isfinite(cv_score)
        @test cv_score >= 0.0
    end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: Fails because `_loocv` is not defined.

- [ ] **Step 3: Implement _loocv**

In `src/identify_rate_equation.jl`, add:

```julia
"""
Subset a columnar NamedTuple by a boolean mask.
"""
function _subset_data(data::NamedTuple, mask)
    return map(col -> col[mask], data)
end

"""
Evaluate loss of a mechanism on test data with given parameters.
"""
function _evaluate_loss(mechanism, data, params, Keq)
    pnames = fitted_params(mechanism)
    x = [log(params[p]) for p in pnames]
    fp = FittingProblem(mechanism, data; Keq=Keq)
    return loss!(x, fp)
end

"""
    _loocv(mechanism, prob; optimizer, kwargs...) → Float64

Leave-one-group-out cross-validation. Fits the mechanism on all groups except
the held-out group, evaluates loss on the held-out group, and returns the mean
CV score across all held-out groups.

`optimizer` is extracted explicitly because `fit_rate_equation` takes it as a
positional argument.
"""
function _loocv(mechanism::AbstractEnzymeMechanism,
    prob::IdentifyRateEquationProblem; optimizer, kwargs...
)
    groups = unique(prob.data.group)
    scores = Float64[]

    for held_out in groups
        train_mask = prob.data.group .!= held_out
        test_mask = prob.data.group .== held_out

        train_data = _subset_data(prob.data, train_mask)
        test_data = _subset_data(prob.data, test_mask)

        fp_train = FittingProblem(mechanism, train_data; Keq=prob.Keq)
        fit = fit_rate_equation(fp_train, optimizer; kwargs...)

        test_loss = _evaluate_loss(mechanism, test_data, fit.params, prob.Keq)
        push!(scores, test_loss)
    end

    return mean(scores)
end
```

- [ ] **Step 4: Run test to verify it passes**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: `_loocv` test passes.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Implement _loocv for leave-one-group-out cross-validation"
```

- [ ] **Step 6: Write test for _cv_model_selection**

In `test/test_identify_rate_equation.jl`, add:

```julia
    # ── Test: _cv_model_selection ─────────────────────────────────────────────
    @testset "_cv_model_selection" begin
        using OptimizationBBO
        prob = IdentifyRateEquationProblem(test_rxn, test_data; Keq=Keq_val)

        specs, df = EnzymeRates._beam_search(prob;
            min_beam_width=200, beam_fraction=0.1, max_param_count=8,
            save_dir=nothing, pmap_function=map,
            optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
            n_restarts=2, maxtime=5.0)

        results = EnzymeRates._cv_model_selection(specs, df, prob;
            n_cv_candidates=3, pmap_function=map,
            optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
            n_restarts=2, maxtime=5.0)

        @test results isa IdentifyRateEquationResults
        @test results.best isa EnzymeRates.AbstractEnzymeMechanism
        @test nrow(results.cv_results) > 0
        @test "cv_score" in names(results.cv_results)
        @test all(isfinite, results.cv_results.cv_score)
    end
```

- [ ] **Step 7: Run test to verify it fails**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: Fails with "Not implemented yet" from `_cv_model_selection` stub.

- [ ] **Step 8: Implement _cv_model_selection**

In `src/identify_rate_equation.jl`, replace the `_cv_model_selection` stub with:

```julia
function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function, optimizer, kwargs...
)
    # specs[i] corresponds to df[i, :] — they're appended in the same order
    # in _beam_search (df is NOT sorted). Use positional indexing.
    df_indexed = copy(df)
    df_indexed.spec_idx = 1:nrow(df_indexed)

    # Group by n_params, take top n_cv_candidates per group by loss
    candidate_indices = Int[]
    for gdf in groupby(df_indexed, :n_params)
        sorted = sort(gdf, :loss)
        n_take = min(n_cv_candidates, nrow(sorted))
        for i in 1:n_take
            push!(candidate_indices, sorted[i, :spec_idx])
        end
    end

    candidate_specs = specs[candidate_indices]
    candidate_rows = df[candidate_indices, :]

    # LOOCV each candidate in parallel
    cv_scores = pmap_function(candidate_specs) do spec
        m = _compile(spec)
        _loocv(m, prob; optimizer, kwargs...)
    end

    # Build CV results DataFrame
    cv_df = copy(candidate_rows)
    cv_df.cv_score = collect(cv_scores)
    cv_df.spec_idx = candidate_indices

    # Select best: find param count with best (lowest) CV score
    best_cv_per_pc = combine(groupby(cv_df, :n_params),
        :cv_score => minimum => :best_cv)
    best_pc_row = best_cv_per_pc[argmin(best_cv_per_pc.best_cv), :]
    best_param_count = best_pc_row.n_params

    # Best mechanism = lowest loss at that param count
    at_best_pc = filter(row -> row.n_params == best_param_count, cv_df)
    sort!(at_best_pc, :loss)
    best_idx = at_best_pc[1, :spec_idx]
    best_mechanism = _compile(specs[best_idx])

    select!(cv_df, Not(:spec_idx))
    return IdentifyRateEquationResults(best_mechanism, cv_df)
end
```

- [ ] **Step 9: Run test to verify it passes**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: `_cv_model_selection` test passes.

- [ ] **Step 10: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Implement _cv_model_selection with LOOCV-based model selection"
```

---

### Task 6: Implement identify_rate_equation integration and CSV save

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write integration test**

In `test/test_identify_rate_equation.jl`, add:

```julia
    # ── Test: identify_rate_equation end-to-end ──────────────────────────────
    @testset "identify_rate_equation end-to-end" begin
        using OptimizationBBO
        prob = IdentifyRateEquationProblem(test_rxn, test_data; Keq=Keq_val)

        results = identify_rate_equation(prob;
            min_beam_width=200, beam_fraction=0.1, max_param_count=8,
            n_cv_candidates=3, save_dir=nothing, pmap_function=map,
            optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
            n_restarts=2, maxtime=5.0)

        @test results isa IdentifyRateEquationResults
        @test results.best isa EnzymeRates.AbstractEnzymeMechanism
        @test nrow(results.cv_results) > 0
    end

    # ── Test: CSV save ───────────────────────────────────────────────────────
    @testset "CSV save" begin
        using OptimizationBBO
        prob = IdentifyRateEquationProblem(test_rxn, test_data; Keq=Keq_val)

        save_dir = mktempdir()
        results = identify_rate_equation(prob;
            min_beam_width=200, beam_fraction=0.1, max_param_count=8,
            n_cv_candidates=3, save_dir=save_dir, pmap_function=map,
            optimizer=BBO_adaptive_de_rand_1_bin_radiuslimited(),
            n_restarts=2, maxtime=5.0)

        # Check CSV files were created
        csv_files = filter(f -> endswith(f, ".csv"), readdir(save_dir))
        @test length(csv_files) > 0

        # Check a CSV file is readable and has expected columns
        first_csv = CSV.read(joinpath(save_dir, csv_files[1]), DataFrame)
        @test "n_params" in names(first_csv)
        @test "loss" in names(first_csv)
        @test "mechanism_type" in names(first_csv)
        @test nrow(first_csv) > 0
    end
```

- [ ] **Step 2: Run test to verify it passes**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: All tests pass including integration and CSV tests.

- [ ] **Step 3: Commit**

```bash
git add test/test_identify_rate_equation.jl
git commit -m "Add integration tests for identify_rate_equation and CSV save"
```

---

### Task 7: Run full test suite and clean up

**Files:**
- Modify: `src/identify_rate_equation.jl` (if cleanup needed)
- Modify: `test/test_identify_rate_equation.jl` (if cleanup needed)

- [ ] **Step 1: Run the full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: All tests pass, including Aqua and JET checks.

- [ ] **Step 2: Fix any Aqua/JET issues**

If Aqua reports stale dependencies or JET reports type instabilities, fix them. Common issues:
- DataFrames/CSV/Statistics may need to be added to `[deps]` (not just extras) if they're used in src/
- Any unresolved type annotations

- [ ] **Step 3: Re-read identify_rate_equation.jl for dead code and simplification**

Read the full file, looking for:
- Unused functions
- Duplicated logic between `_save_level_csv` and `_rows_to_dataframe`
- Any code that can be simplified

- [ ] **Step 4: Commit any cleanup**

```bash
git add -u
git commit -m "Clean up identify_rate_equation implementation"
```

- [ ] **Step 5: Run full test suite one final time**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: All tests pass.
