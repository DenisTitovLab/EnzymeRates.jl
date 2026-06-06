# Rate-Scaling, Fit Robustness, and Progress Output Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a single `scale_k_to_kcat` knob (relative-centered vs absolute-turnover loss + rescale), record solver retcodes, make fit/CV failures loud instead of silent, and emit cluster-visible progress output — all in the `fit → identify_rate_equation` path.

**Architecture:** `scale_k_to_kcat::Union{Real,Nothing}` becomes a field on `FittingProblem` and `IdentifyRateEquationProblem` (set at construction, like `Keq`): a `Real` (default `1.0`) means relative data (per-group-centered loss + rescale SS k's to that kcat); `nothing` means absolute per-enzyme turnover (uncentered loss, raw k's). `fit_rate_equation` reads the field and returns `(params, loss, retcode)`. `_process_batch` captures per-mechanism exceptions into a `(entries, failures)` split; an all-failed base tier re-raises the real exception; `_loocv` raises instead of silently returning empty. A `_progress` helper writes flushed stdout + an appended `<save_dir>/progress.log`.

**Tech Stack:** Julia, Optimization.jl (+ SciMLBase `ReturnCode`, transitively), DataFrames.jl, CSV.jl. Tests use `OptimizationBBO` and `OptimizationPyCMA` (test-only deps via `[extras]`/`[targets]`).

---

## Running tests

This repo has **no per-file test runner** — `[extras]`/`[targets]` build a temp test env. Every "Run" step uses the full suite:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Cold runs pay precompile + JIT. Look for the named `@testset` in the output. TDD here is **per task** (write the task's tests → run → implement → run → commit), not per individual assertion, because each run is the whole suite.

The spec this implements: `docs/superpowers/specs/2026-06-05-rate-scaling-fit-robustness-progress-design.md`.

---

## File Structure

- **Modify `src/fitting.jl`** — `FittingProblem` struct + constructors gain `scale_k_to_kcat`; `loss!` branches centering on it; `fit_rate_equation` drops its `kcat` kwarg, reads the field, returns `retcode`. (Tasks 1, 2)
- **Modify `src/rate_eq_derivation.jl`** — rename `rescale_parameter_values` kwarg `kcat` → `scale_k_to_kcat`. (Task 2)
- **Modify `src/identify_rate_equation.jl`** — `IdentifyRateEquationProblem` gains `scale_k_to_kcat`; `BatchEntry` gains `retcode`; new `FitFailure` struct, `_failure_row`, `_exc_string`, `_progress`, `_batch_summary` helpers; `_rows_to_dataframe` gains `retcode`/`error` columns; `_process_batch` returns `(entries, failures)`; `_beam_search` re-raises an all-failed base tier and writes failure rows; `_loocv` is loud; `_cv_model_selection`/`identify_rate_equation` thread `show_progress`/`save_dir` and emit stages. (Tasks 3–8)
- **Modify `test/test_fitting.jl`** — absolute-mode loss tests, validation test, migrated kcat-normalization test. (Tasks 1, 2)
- **Modify `test/test_rate_eq_derivation.jl`** — rename the one explicit `rescale_parameter_values(...; kcat=...)` call. (Task 2)
- **Modify `test/test_identify_rate_equation.jl`** — construction field test, migrated `_process_batch`/`_ingest!`/`_rows_to_dataframe`/csv-writer fixtures, loud-`_loocv` test, `_progress` test, progress.log assertion. (Tasks 3–8)
- **Modify `.claude/CLAUDE.md`** — doc sync. (Task 9)

---

## Task 1: `scale_k_to_kcat` field on `FittingProblem` + `loss!` centering branch

**Files:**
- Modify: `src/fitting.jl:18-25` (struct), `src/fitting.jl:38-83` (constructor), `src/fitting.jl:88-90` (forwarding constructor), `src/fitting.jl:134-147` (`loss!` pass 2)
- Test: `test/test_fitting.jl`

- [ ] **Step 1: Write the failing tests**

Add these two `@testset`s to `test/test_fitting.jl`, immediately after the existing "Multi-group centering invariance" testset (after line 150, before "Sign mismatch penalty"):

```julia
    # ── Absolute mode: uncentered loss (scale_k_to_kcat=nothing) ──────────────
    @testset "Absolute mode uncentered loss" begin
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0, Keq = Keq_val, E_total = 1.0)
        concs_list = [
            (S = 1.0, P = 0.1), (S = 2.0, P = 0.1), (S = 5.0, P = 0.1),
            (S = 1.0, P = 0.5), (S = 2.0, P = 0.5),
        ]
        data = make_synthetic_data(uni_uni, true_params, concs_list)
        pn = EnzymeRates.fitted_params(uni_uni)
        x_true = [log(true_params[p]) for p in pn]

        fp_rel = FittingProblem(uni_uni, data; Keq=Keq_val)                        # default 1.0
        fp_abs = FittingProblem(uni_uni, data; Keq=Keq_val, scale_k_to_kcat=nothing)

        # At true params both modes are ~0 (predictions match data exactly).
        @test EnzymeRates.loss!(x_true, fp_rel) ≈ 0.0 atol=1e-20
        @test EnzymeRates.loss!(x_true, fp_abs) ≈ 0.0 atol=1e-20

        # Scale every rate by 3 (a pure per-group offset). Relative loss is
        # invariant (centering removes it); absolute loss sees it: every
        # residual becomes log(3), so absolute loss = log(3)^2.
        data3 = merge(data, (Rate = data.Rate .* 3.0,))
        fp_rel3 = FittingProblem(uni_uni, data3; Keq=Keq_val)
        fp_abs3 = FittingProblem(uni_uni, data3; Keq=Keq_val, scale_k_to_kcat=nothing)
        @test EnzymeRates.loss!(x_true, fp_rel3) ≈ 0.0 atol=1e-20
        @test EnzymeRates.loss!(x_true, fp_abs3) ≈ log(3.0)^2 rtol=1e-8
    end

    # ── scale_k_to_kcat validation ────────────────────────────────────────────
    @testset "scale_k_to_kcat validation" begin
        ok_data = (group = ["G1"], Rate = [1.0], S = [1.0], P = [0.1])
        @test_throws ErrorException FittingProblem(uni_uni, ok_data; Keq=1.0, scale_k_to_kcat=0.0)
        @test_throws ErrorException FittingProblem(uni_uni, ok_data; Keq=1.0, scale_k_to_kcat=-5.0)
        @test FittingProblem(uni_uni, ok_data; Keq=1.0, scale_k_to_kcat=nothing) isa FittingProblem
        @test FittingProblem(uni_uni, ok_data; Keq=1.0) isa FittingProblem  # default 1.0
    end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL in the "Fitting" testset — `FittingProblem` has no `scale_k_to_kcat` keyword (`MethodError`/unknown kwarg).

- [ ] **Step 3: Implement**

In `src/fitting.jl`, add the field to the struct (replace lines 18-25):

```julia
struct FittingProblem{M<:AbstractEnzymeMechanism, D<:NamedTuple}
    mechanism::M
    data::D
    group_point_indexes::Vector{Vector{Int}}
    Keq::Float64
    scale_k_to_kcat::Union{Float64,Nothing}
    log_abs_rates::Vector{Float64}
    log_ratios_buffer::Vector{Float64}
end
```

Add a field line to the struct docstring (after the `Keq` bullet, currently `src/fitting.jl:14`):

```julia
- `scale_k_to_kcat`: a positive `Real` selects relative mode (per-group-centered
  loss); `nothing` selects absolute per-enzyme-turnover mode (uncentered loss)
```

Update the main constructor signature + body (replace lines 38-39 and the final `FittingProblem{...}(...)` call at 79-82):

```julia
function FittingProblem(mechanism::AbstractEnzymeMechanism, table;
        Keq::Real, scale_k_to_kcat::Union{Real,Nothing}=1.0)
    scale_k_to_kcat !== nothing && scale_k_to_kcat <= 0 && error(
        "scale_k_to_kcat must be positive (or nothing); got $scale_k_to_kcat")
```

```julia
    sk = scale_k_to_kcat === nothing ? nothing : Float64(scale_k_to_kcat)
    FittingProblem{typeof(mechanism), typeof(data)}(
        mechanism, data, group_point_indexes, Float64(Keq), sk,
        log_abs_rates, log_ratios_buffer
    )
end
```

Update the forwarding constructor (replace lines 88-90):

```julia
FittingProblem(mechanism::Union{Mechanism, AllostericMechanism}, table;
        Keq::Real, scale_k_to_kcat::Union{Real,Nothing}=1.0) =
    FittingProblem(compile_mechanism(mechanism), table;
        Keq=Keq, scale_k_to_kcat=scale_k_to_kcat)
```

Replace `loss!` pass 2 (lines 134-147) with the branch:

```julia
    # Pass 2: loss. Relative (scale_k_to_kcat isa Real) → per-group mean-
    # centered, removing each group's arbitrary scale. Absolute
    # (scale_k_to_kcat === nothing) → uncentered: the y-axis is absolute
    # per-enzyme turnover, so the absolute magnitude is meaningful.
    total_loss = 0.0
    if fp.scale_k_to_kcat === nothing
        @inbounds for i in 1:n_data
            total_loss += buf[i] * buf[i]
        end
    else
        @inbounds for grp_idx in fp.group_point_indexes
            n_grp = length(grp_idx)
            mean_lr = 0.0
            for j in grp_idx
                mean_lr += buf[j]
            end
            mean_lr /= n_grp
            for j in grp_idx
                d = buf[j] - mean_lr
                total_loss += d * d
            end
        end
    end
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — new testsets green; existing "Centering invariance", "Multi-group centering invariance", and "Zero allocations" (default `scale_k_to_kcat=1.0` → `else` branch) stay green.

- [ ] **Step 5: Commit**

```bash
git add src/fitting.jl test/test_fitting.jl
git commit -m "Add scale_k_to_kcat field to FittingProblem; absolute (uncentered) loss mode" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `fit_rate_equation` reads `fp.scale_k_to_kcat`, returns `retcode`; rename `rescale_parameter_values` kwarg

**Files:**
- Modify: `src/rate_eq_derivation.jl:941-957` (`rescale_parameter_values`)
- Modify: `src/fitting.jl:161-218` (`fit_rate_equation`)
- Test/migrate: `test/test_fitting.jl:303-333` (kcat-normalization testset), `test/test_rate_eq_derivation.jl:903`

- [ ] **Step 1: Write/migrate the tests**

Replace the entire "kcat normalization" testset in `test/test_fitting.jl` (lines 303-333) with:

```julia
    # ── Test 9: scale_k_to_kcat normalization + retcode ────────────────
    @testset "scale_k_to_kcat normalization" begin
        using OptimizationBBO
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0, Keq = Keq_val, E_total = 1.0)

        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
            (S = 0.5, P = 0.5), (S = 1.0, P = 0.5), (S = 2.0, P = 0.5),
        ]
        data = make_synthetic_data(uni_uni, true_params, concs_list)

        # Default scale_k_to_kcat=1.0: returned params have kcat ≈ 1.
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)
        result = fit_rate_equation(fp, BBO_adaptive_de_rand_1_bin_radiuslimited();
            n_restarts=3, maxtime=5.0)
        full = merge(result.params, (Keq = Keq_val, E_total = 1.0))
        @test EnzymeRates._kcat_forward(uni_uni, full) ≈ 1.0 rtol=0.01
        @test result.retcode isa Symbol

        # Custom target set on the FittingProblem.
        fp42 = FittingProblem(uni_uni, data; Keq=Keq_val, scale_k_to_kcat=42.0)
        result2 = fit_rate_equation(fp42, BBO_adaptive_de_rand_1_bin_radiuslimited();
            n_restarts=3, maxtime=5.0)
        full2 = merge(result2.params, (Keq = Keq_val, E_total = 1.0))
        @test EnzymeRates._kcat_forward(uni_uni, full2) ≈ 42.0 rtol=0.01

        # scale_k_to_kcat=nothing: raw params (no rescale), retcode still present.
        fpN = FittingProblem(uni_uni, data; Keq=Keq_val, scale_k_to_kcat=nothing)
        result3 = fit_rate_equation(fpN, BBO_adaptive_de_rand_1_bin_radiuslimited();
            n_restarts=3, maxtime=5.0)
        @test haskey(result3, :params)
        @test result3.retcode isa Symbol
    end
```

In `test/test_rate_eq_derivation.jl`, change line 903 from:

```julia
        norm_custom = rescale_parameter_values(m, params; kcat=kcat_target)
```

to:

```julia
        norm_custom = rescale_parameter_values(m, params; scale_k_to_kcat=kcat_target)
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `fit_rate_equation` still has a `kcat` kwarg (so `fp42`/`fpN` don't influence rescale), `result.retcode` doesn't exist, and `rescale_parameter_values(...; scale_k_to_kcat=...)` is an unknown kwarg.

- [ ] **Step 3: Implement**

In `src/rate_eq_derivation.jl`, replace the `rescale_parameter_values` docstring + function (lines 941-957):

```julia
"""
    rescale_parameter_values(m, params::NamedTuple; scale_k_to_kcat=1.0)

Rescale SS rate constants so that `_kcat_forward(m, result) ≈ scale_k_to_kcat`.
Non-SS parameters (K's, Keq, E_total, L, regulatory K's) are unchanged.
"""
function rescale_parameter_values(
    m::_AnyMechanism, params::NamedTuple; scale_k_to_kcat=1.0,
)
    kcat_current = _kcat_forward(m, params)
    scale = scale_k_to_kcat / kcat_current
    ss_names = _ss_rate_constant_names(m)
    NamedTuple{keys(params)}(Tuple(
        k in ss_names ? v * scale : v
        for (k, v) in zip(keys(params), values(params))
    ))
end
```

In `src/fitting.jl`, replace the `fit_rate_equation` docstring + function (lines 161-218). New docstring + body:

```julia
"""
    fit_rate_equation(fp::FittingProblem, optimizer;
        n_restarts=20, maxtime=60.0,
        lb=fill(-15.0, length(fitted_params(fp.mechanism))),
        ub=fill(15.0, length(fitted_params(fp.mechanism))),
        kwargs...)

Fit rate constants by minimizing `loss!` using Optimization.jl.

Runs `n_restarts` independent optimizations from random initial points and returns
the best result.

Rescaling is driven by `fp.scale_k_to_kcat`: when it is a `Real`, the returned
parameters are rescaled so that `_kcat_forward(mechanism, params) ≈
fp.scale_k_to_kcat`; when it is `nothing`, the raw (unrescaled) parameters are
returned (the data fixes the absolute scale).

Returns a NamedTuple `(params, loss, retcode)` where:
- `params`: fitted rate constants as a NamedTuple
- `loss`: the best loss value achieved
- `retcode`: the `Symbol` form of the best restart's `sol.retcode`
  (e.g. `:Success`, `:MaxTime`). `:MaxTime` means the optimizer hit its time
  budget before converging — treat the fit as un-converged.
"""
function fit_rate_equation(fp::FittingProblem, optimizer;
    n_restarts::Int=20,
    maxtime::Real=60.0,
    lb=fill(-15.0, length(fitted_params(fp.mechanism))),
    ub=fill(15.0, length(fitted_params(fp.mechanism))),
    kwargs...
)
    obj = Optimization.OptimizationFunction((x, p) -> loss!(x, p))
    np = length(fitted_params(fp.mechanism))

    best_x = zeros(np)
    best_loss = Inf
    best_retcode = :Default

    for _ in 1:n_restarts
        x0 = clamp.(randn(np) .* 2.0, lb, ub)
        prob = Optimization.OptimizationProblem(obj, x0, fp; lb=lb, ub=ub)
        sol = Optimization.solve(prob, optimizer; maxtime=maxtime, kwargs...)
        if sol.objective < best_loss
            best_loss = sol.objective
            best_x .= sol.u
            best_retcode = Symbol(sol.retcode)
        end
    end

    pnames = fitted_params(fp.mechanism)
    result_params = NamedTuple{pnames}(ntuple(i -> exp(best_x[i]), Val(length(pnames))))
    if fp.scale_k_to_kcat !== nothing
        full = merge(result_params, (Keq = fp.Keq, E_total = 1.0))
        rp = rescale_parameter_values(
            fp.mechanism, full; scale_k_to_kcat=fp.scale_k_to_kcat,
        )
        result_params = NamedTuple{pnames}(
            ntuple(i -> rp[pnames[i]], Val(length(pnames))),
        )
    end
    return (params = result_params, loss = best_loss, retcode = best_retcode)
end
```

(`Symbol(sol.retcode)` returns the bare enum name like `:Success`; SciMLBase is loaded transitively via Optimization, so no import is needed.)

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — `scale_k_to_kcat normalization` and the renamed `rescale_parameter_values` call green; `test_kcat_rescaling` still green.

- [ ] **Step 5: Commit**

```bash
git add src/fitting.jl src/rate_eq_derivation.jl test/test_fitting.jl test/test_rate_eq_derivation.jl
git commit -m "fit_rate_equation reads fp.scale_k_to_kcat and returns retcode; rename rescale kwarg" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `scale_k_to_kcat` field on `IdentifyRateEquationProblem`

**Files:**
- Modify: `src/identify_rate_equation.jl:21-70` (struct + constructor)
- Test: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write the failing test**

Add to the "construction" testset in `test/test_identify_rate_equation.jl` (after line 90, inside the testset):

```julia
        # scale_k_to_kcat field: default 1.0, settable to nothing, validated.
        @test prob.scale_k_to_kcat == 1.0
        prob_abs = IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, scale_k_to_kcat=nothing)
        @test prob_abs.scale_k_to_kcat === nothing
        @test_throws ErrorException IdentifyRateEquationProblem(
            test_rxn, test_data; Keq=Keq_val, scale_k_to_kcat=0.0)
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `prob.scale_k_to_kcat` undefined / unknown kwarg.

- [ ] **Step 3: Implement**

In `src/identify_rate_equation.jl`, add the field to the struct (replace lines 21-27):

```julia
struct IdentifyRateEquationProblem{
    R<:EnzymeReaction, D<:NamedTuple
}
    reaction::R
    data::D
    Keq::Float64
    scale_k_to_kcat::Union{Float64,Nothing}
end
```

Update the constructor signature (line 29-31) and final construction (lines 67-69):

```julia
function IdentifyRateEquationProblem(
    reaction::EnzymeReaction, table; Keq::Real,
    scale_k_to_kcat::Union{Real,Nothing}=1.0
)
    scale_k_to_kcat !== nothing && scale_k_to_kcat <= 0 && error(
        "scale_k_to_kcat must be positive (or nothing); got $scale_k_to_kcat")
```

```julia
    sk = scale_k_to_kcat === nothing ? nothing : Float64(scale_k_to_kcat)
    IdentifyRateEquationProblem{
        typeof(reaction),typeof(data)
    }(reaction, data, Float64(Keq), sk)
end
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — construction testset green; existing identify tests unaffected (default `1.0`).

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add scale_k_to_kcat field to IdentifyRateEquationProblem" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: `BatchEntry.retcode` + `retcode`/`error` row columns in `_rows_to_dataframe`

**Files:**
- Modify: `src/identify_rate_equation.jl:317-324` (`BatchEntry`), `src/identify_rate_equation.jl:238-265` (`_rows_to_dataframe`)
- Migrate fixtures: `test/test_identify_rate_equation.jl:137-159` (`_rows_to_dataframe`), `:370-389` (csv writers), `:898-922` (`_ingest!`)

- [ ] **Step 1: Write/migrate the tests**

In `test/test_identify_rate_equation.jl`, update the `_rows_to_dataframe` testset fixture (lines 138-153) to include the new fields and assert the new columns:

```julia
        rows = [(
            n_params = 3,
            loss = 0.5,
            mechanism_type = "test",
            rate_equation = "v = ...",
            retcode = "Success",
            error = missing,
            fitted_param_names = (:a, :b),
            fitted_param_values = (1.0, 2.0),
            eq_hash = "0123456789abcdef",
        )]
        df = EnzymeRates._rows_to_dataframe(
            rows)
        @test nrow(df) == 1
        @test "a" in names(df)
        @test "b" in names(df)
        @test "eq_hash" in names(df)
        @test "retcode" in names(df)
        @test "error" in names(df)
        @test !("fit_inherited_from_estimate" in names(df))
```

Update the "csv writers" fixture (lines 371-375) to include the new fields:

```julia
    rows = [(
        n_params = 5, loss = 1.0, mechanism_type = "M",
        rate_equation = "v = 1", retcode = "Success", error = missing,
        fitted_param_names = (:K_a,),
        fitted_param_values = (2.0,), eq_hash = "abc",
    )]
```

Update the `_ingest!`-test `mk` helper (lines 899-905) so its `BatchEntry` call takes the new `retcode` field and its row carries `retcode`/`error`:

```julia
    mk(n, loss, h) = EnzymeRates.BatchEntry(
        first(EnzymeRates.init_mechanisms(@enzyme_reaction begin
            substrates:S[C]; products:P[C] end)),
        n, loss, :Success, hash(h),
        (n_params=n, loss=loss, mechanism_type="M",
         rate_equation="v", retcode="Success", error=missing,
         fitted_param_names=(:K,),
         fitted_param_values=(1.0,), eq_hash=string(hash(h),base=16,pad=16)))
```

Add a new testset (a mixed fitted + failure-row DataFrame) right after the `_rows_to_dataframe` testset (after line 159):

```julia
    @testset "_rows_to_dataframe with failure row" begin
        rows = [
            (n_params = 3, loss = 0.5, mechanism_type = "M",
             rate_equation = "v = ...", retcode = "Success", error = missing,
             fitted_param_names = (:a,), fitted_param_values = (1.0,),
             eq_hash = "0123456789abcdef"),
            (n_params = missing, loss = missing, mechanism_type = "M",
             rate_equation = missing, retcode = missing,
             error = "StackOverflowError: ", fitted_param_names = (),
             fitted_param_values = (), eq_hash = missing),
        ]
        df = EnzymeRates._rows_to_dataframe(rows)
        @test nrow(df) == 2
        @test ismissing(df.loss[2])
        @test df.error[2] == "StackOverflowError: "
        @test ismissing(df.eq_hash[2])
        @test "a" in names(df)              # param column still built from row 1
        @test ismissing(df.a[2])            # failure row contributes no param value
    end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `_rows_to_dataframe` has no `retcode`/`error` columns; `BatchEntry` has no 6-arg constructor.

- [ ] **Step 3: Implement**

In `src/identify_rate_equation.jl`, replace the `BatchEntry` struct (lines 317-324):

```julia
"""One fitted mechanism: its own params + retcode + eq_hash + the CSV row."""
struct BatchEntry
    mech::Union{Mechanism, AllostericMechanism}
    n_params::Int
    loss::Float64
    retcode::Symbol
    eq_hash::UInt64
    row::NamedTuple
end
```

In `_rows_to_dataframe`, replace the `DataFrame(...)` construction (lines 248-254) to add the two columns:

```julia
    df = DataFrame(
        n_params = [r.n_params for r in rows],
        loss = [r.loss for r in rows],
        mechanism_type = [r.mechanism_type for r in rows],
        rate_equation = [r.rate_equation for r in rows],
        retcode = [r.retcode for r in rows],
        error = [r.error for r in rows],
        eq_hash = [r.eq_hash for r in rows],
    )
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — `_rows_to_dataframe`, the new failure-row testset, csv writers, and `_ingest!` testsets green.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add retcode to BatchEntry; retcode/error columns in _rows_to_dataframe" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: `_process_batch` captures failures → `(entries, failures)`

**Files:**
- Modify: `src/identify_rate_equation.jl:317-366` (add `FitFailure`, `_exc_string`, `_failure_row`; rewrite `_process_batch`)
- Migrate: `test/test_identify_rate_equation.jl:871-896` (`_process_batch` testset)

- [ ] **Step 1: Write/migrate the tests**

Replace the `_process_batch` testset (lines 871-896) with one that destructures `(entries, failures)` and exercises both the success and all-fail paths:

```julia
@testset "_process_batch" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = DataFrame(
        S = [1.0, 2.0, 3.0, 4.0],
        P = [0.1, 0.2, 0.3, 0.4],
        Rate = [0.5, 0.8, 1.0, 1.1],
        group = [1, 1, 2, 2],
    )
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    ms = EnzymeRates._dedup_flat!(collect(EnzymeRates.init_mechanisms(rxn)))

    entries, failures = EnzymeRates._process_batch(ms, prob;
        pmap_function=map, optimizer=PyCMAOpt(),
        max_param_count=20, n_restarts=1, maxtime=1.0)
    @test !isempty(entries)
    @test isempty(failures)
    @test all(e -> e isa EnzymeRates.BatchEntry, entries)
    @test all(e -> e.retcode isa Symbol, entries)
    @test all(e -> e.n_params == length(e.row.fitted_param_names), entries)
    @test all(e -> occursin(r"^[0-9a-f]{16}$", e.row.eq_hash), entries)

    # cap filter: nothing over the cap is fit (and it is not a failure).
    capped_entries, capped_failures = EnzymeRates._process_batch(ms, prob;
        pmap_function=map, optimizer=PyCMAOpt(),
        max_param_count=0, n_restarts=1, maxtime=1.0)
    @test isempty(capped_entries)
    @test isempty(capped_failures)

    # config error (unknown optimizer kwarg) → every fit throws → all failures,
    # no entries; each failure carries a non-empty error string.
    fail_entries, fail_failures = EnzymeRates._process_batch(ms, prob;
        pmap_function=map, optimizer=PyCMAOpt(),
        max_param_count=20, n_restarts=1, maxtime=1.0, beam_fraction=0.5)
    @test isempty(fail_entries)
    @test !isempty(fail_failures)
    @test all(f -> f isa EnzymeRates.FitFailure, fail_failures)
    @test all(f -> !isempty(f.error), fail_failures)
end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `_process_batch` returns a single vector, not `(entries, failures)`; `FitFailure` is undefined.

- [ ] **Step 3: Implement**

In `src/identify_rate_equation.jl`, add the `FitFailure` struct and helpers just after the `BatchEntry` struct (after line 324):

```julia
"""A mechanism that threw during compile/fit; kept for the CSV + summary."""
struct FitFailure
    mech::Union{Mechanism, AllostericMechanism}
    error::String
end

# Compact, CSV-safe rendering of a thrown exception: type + truncated message.
_exc_string(e) = first(sprint(showerror, e), 200)

# CSV row for a mechanism that threw. Same NamedTuple schema as a fitted row,
# with `missing` wherever the value is unavailable (compile/fit never produced it).
function _failure_row(f::FitFailure)
    (n_params = missing,
     loss = missing,
     mechanism_type = string(typeof(f.mech)),
     rate_equation = missing,
     retcode = missing,
     error = f.error,
     fitted_param_names = (),
     fitted_param_values = (),
     eq_hash = missing)
end
```

Replace the `_process_batch` docstring + function (lines 326-366):

```julia
"""
Compile + cap-check + fit every mechanism in `mechs`, one `pmap` pass with
compile and fit fused on the same worker. Returns `(entries, failures)`:
`entries::Vector{BatchEntry}` are the fitted mechanisms (each keeping its OWN
fitted params, `retcode`, and `eq_hash`); `failures::Vector{FitFailure}` are
mechanisms that threw (StackOverflow/OOM at compile, an unsupported-kwarg error,
etc.) — captured WITH the exception text, never silently swallowed. A mechanism
whose actual fitted-param count exceeds `max_param_count` is dropped BEFORE
fitting (a cap skip, not a failure). No dedup here — `mechs` is already
structurally deduped by the caller (`_dedup_flat!`).
"""
function _process_batch(
    mechs, prob::IdentifyRateEquationProblem;
    pmap_function, optimizer, max_param_count, kwargs...
)
    results = pmap_function(mechs) do m
        try
            em = compile_mechanism(m)
            fkeys = fitted_params(em)
            n = length(fkeys)
            n > max_param_count && return nothing
            eq_text = rate_equation_string(em)
            key = _rate_eq_dedup_key(eq_text)
            fp = FittingProblem(em, prob.data;
                Keq=prob.Keq, scale_k_to_kcat=prob.scale_k_to_kcat)
            fit = fit_rate_equation(fp, optimizer; kwargs...)
            row = (
                n_params = n,
                loss = fit.loss,
                mechanism_type = string(typeof(em)),
                rate_equation = eq_text,
                retcode = string(fit.retcode),
                error = missing,
                fitted_param_names = fkeys,
                fitted_param_values =
                    Tuple(fit.params[k] for k in fkeys),
                eq_hash = string(key, base=16, pad=16),
            )
            BatchEntry(m, n, fit.loss, fit.retcode, key, row)
        catch e
            FitFailure(m, _exc_string(e))
        end
    end
    entries  = BatchEntry[r for r in results if r isa BatchEntry]
    failures = FitFailure[r for r in results if r isa FitFailure]
    return entries, failures
end
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL still — `_beam_search` (Task 6) and `_cv_model_selection` still call `_process_batch` expecting a single vector. That is fixed in Task 6; the `_process_batch` *unit* testset itself should pass. (If you want the suite green before Task 6, do Tasks 5 and 6 back-to-back; they are coupled by the return-type change.)

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "_process_batch captures per-mechanism failures into (entries, failures)" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: `_beam_search` consumes `(entries, failures)`; re-raise an all-failed base tier; write failure rows

**Files:**
- Modify: `src/identify_rate_equation.jl:421-429` (base tier), `:450-468` (expansion block)
- Existing test `beam_fraction kwarg removed` (`test/test_identify_rate_equation.jl:439-459`) already asserts the all-fail re-raise — it must stay green.

- [ ] **Step 1: (no new test — covered by existing `beam_fraction` + `_process_batch` tests)**

The contract "config error → all base fits fail → `identify_rate_equation` raises `ErrorException`" is already pinned by the `beam_fraction kwarg removed` testset (lines 439-459). This task makes that path consume the new return shape and embed the captured exception in the message.

- [ ] **Step 2: Run to verify the current breakage**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `_beam_search` destructures `_process_batch`'s result as a single vector (`base_entries = _process_batch(...)`), now a `Tuple`, so `isempty(base_entries)` and `[e.row for e in base_entries]` break.

- [ ] **Step 3: Implement**

In `src/identify_rate_equation.jl`, replace the base-tier block (lines 421-429):

```julia
    # ── Base tier: fit ALL init mechanisms (no bucketing — siblings) ──
    base = _dedup_flat!(collect(init_mechanisms(prob.reaction)))
    base_entries, base_failures = _process_batch(base, prob;
        pmap_function, optimizer, max_param_count, kwargs...)
    if isempty(base_entries)
        isempty(base_failures) && return (
            Union{Mechanism, AllostericMechanism}[],
            _rows_to_dataframe(NamedTuple[]))
        error("Every base-tier fit failed ($(length(base_failures)) " *
              "mechanisms). This usually indicates an optimizer/solver " *
              "configuration problem (e.g. an unsupported kwarg). First " *
              "failure: $(base_failures[1].error)")
    end
    _save_initial_csv(save_dir,
        vcat([e.row for e in base_entries],
             [_failure_row(f) for f in base_failures]))
    _ingest!(frontier, cv_pool, best_loss_by_count,
             base_entries; n_cv_candidates)
```

Replace the expansion block (lines 450-468):

```julia
        if !isempty(to_expand)
            # Typed for dispatch: expand_mechanisms needs a concrete
            # Vector{<:Union{Mechanism, AllostericMechanism}} eltype.
            parents = Union{Mechanism, AllostericMechanism}[
                e.mech for e in to_expand]
            children = _dedup_flat!(
                expand_mechanisms(parents, prob.reaction))
            child_entries, child_failures = _process_batch(children, prob;
                pmap_function, optimizer, max_param_count, kwargs...)
            if !isempty(child_entries) || !isempty(child_failures)
                # Count only iterations that produced rows, so the
                # equation_search_iteration_N CSVs are gap-free.
                iteration += 1
                _save_iteration_csv(save_dir,
                    vcat([e.row for e in child_entries],
                         [_failure_row(f) for f in child_failures]),
                    iteration)
            end
            !isempty(child_entries) && _ingest!(
                frontier, cv_pool, best_loss_by_count,
                child_entries; n_cv_candidates)
        end
```

In `test/test_identify_rate_equation.jl`, update the "CSV output (new schema)" assertions that read columns now possibly carrying `missing` from failure rows (lines 310-312) to skip missing:

```julia
            @test "eq_hash" in names(df_file)
            @test all(length.(string.(skipmissing(df_file.eq_hash))) .== 16)
            @test all(<=(8), skipmissing(df_file.n_params))     # max_param_count=8
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — `_process_batch` testset, `beam_fraction` re-raise testset, the full identify integration, and CSV-output assertions all green.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "_beam_search: re-raise all-failed base tier; write failure rows to CSVs" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: `_loocv` loud + thread `scale_k_to_kcat`; `_cv_model_selection` cleanup

**Files:**
- Modify: `src/identify_rate_equation.jl:495-502` (`_evaluate_loss`), `:504-563` (`_loocv`), `:776-784` (`_cv_model_selection` cv_score)
- Test: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write the failing test**

Add a new testset to `test/test_identify_rate_equation.jl` right after the existing "_loocv returns per-fold scores, floored at eps" testset (after line 366, still inside the outer `@testset "identify_rate_equation"`):

```julia
    @testset "_loocv is loud on fit failure" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
        end
        m = EnzymeRates.EnzymeMechanism(
            first(EnzymeRates.init_mechanisms(rxn)))
        data = DataFrame(
            S    = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            P    = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            Rate = [0.5, 0.8, 1.0, 1.1, 1.2, 1.3],
            group = [1, 1, 2, 2, 3, 3],
        )
        prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
        # An unsupported optimizer kwarg makes every fold fit throw; _loocv
        # must NOT swallow it (the old behavior returned Float64[]).
        @test_throws Exception EnzymeRates._loocv(
            m, prob; optimizer=PyCMAOpt(),
            n_restarts=1, maxtime=1.0, beam_fraction=0.5)
    end
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — current `_loocv` catches the exception and returns `Float64[]`, so nothing is thrown and `@test_throws` fails.

- [ ] **Step 3: Implement**

In `src/identify_rate_equation.jl`, replace `_evaluate_loss` (lines 491-502):

```julia
"""
Evaluate loss of a mechanism on data with given
params and scaling mode.
"""
function _evaluate_loss(
    mechanism, data, params, Keq, scale_k_to_kcat
)
    pnames = fitted_params(mechanism)
    x = [log(params[p]) for p in pnames]
    fp = FittingProblem(mechanism, data; Keq=Keq, scale_k_to_kcat=scale_k_to_kcat)
    return loss!(x, fp)
end
```

Replace the `_loocv` docstring + function (lines 504-563):

```julia
"""
    _loocv(mechanism, prob; optimizer, kwargs...)

Leave-one-group-out cross-validation. Returns `Vector{Float64}` of per-fold
test losses (one per held-out group). Each score is floored at `eps(Float64)`
so `log(score)` is finite — the selection rules operate in log space.

Loud by design: a fold fit that throws propagates, and a non-finite fold test
loss raises (naming the held-out group). A corrupted CV invalidates model
selection, and CV is cheap to recompute from the saved search CSVs — so the
pipeline aborts rather than silently dropping the candidate.
"""
function _loocv(
    mechanism::AbstractEnzymeMechanism,
    prob::IdentifyRateEquationProblem;
    optimizer, kwargs...
)
    groups = unique(prob.data.group)
    scores = Float64[]

    for held_out in groups
        train_mask = prob.data.group .!= held_out
        test_mask  = prob.data.group .== held_out

        train_data = _subset_data(prob.data, train_mask)
        test_data  = _subset_data(prob.data, test_mask)

        fp_train = FittingProblem(mechanism, train_data;
            Keq=prob.Keq, scale_k_to_kcat=prob.scale_k_to_kcat)
        fit = fit_rate_equation(fp_train, optimizer; kwargs...)

        test_loss = _evaluate_loss(mechanism, test_data,
            fit.params, prob.Keq, prob.scale_k_to_kcat)
        # A non-finite fold loss means the fit is unusable; aborting model
        # selection is correct (re-run CV from the saved CSVs after fixing
        # the fit). max(NaN, eps) === NaN, so the floor below would not catch it.
        isfinite(test_loss) || error(
            "LOOCV produced a non-finite test loss for held-out group " *
            "$held_out — the fit is unusable; aborting model selection.")
        # Floor at eps so log(score) is finite. The centered-residuals loss
        # can be exactly 0 (e.g. a single-row held-out group).
        push!(scores, max(test_loss, eps(Float64)))
    end

    scores
end
```

In `_cv_model_selection`, replace the cv_score computation + the all-non-finite backstop (lines 776-784) with the simplified form (every `_loocv` now returns full finite scores or raises):

```julia
    cv_df = copy(candidate_rows)
    cv_df.cv_fold_scores = collect(fold_scores_per_candidate)
    cv_df.cv_score = [mean(log.(v)) for v in cv_df.cv_fold_scores]
```

(Delete the `all(!isfinite, cv_df.cv_score) && error(...)` block entirely. Leave `_select_best_n_params`'s own empty-fold handling untouched — its direct unit tests at lines 718-757 still feed it hand-built empty/partial fold-score vectors.)

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — new loud testset green; the existing "_loocv returns per-fold scores" success test still green; `_select_best_n_params` edge-case tests still green; full identify integration green.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "_loocv raises instead of swallowing; thread scale_k_to_kcat; simplify cv_score" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: `_progress` helper + `show_progress` wiring + stage output

**Files:**
- Modify: `src/identify_rate_equation.jl` — add `_progress` + `_batch_summary`; `identify_rate_equation` signature + call sites (`:162-212`); `_beam_search` signature + stages (`:410-479`); `_cv_model_selection` signature + stages (`:733-825`)
- Test: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Write the failing tests**

Add a standalone testset at the very end of `test/test_identify_rate_equation.jl` (after the final top-level testset, e.g. after line 869's `_select_beam best_override`):

```julia
@testset "_progress" begin
    mktempdir() do tmp
        # show_progress=true: writes to progress.log AND to stdout.
        out_file = joinpath(tmp, "stdout.txt")
        open(out_file, "w") do io
            redirect_stdout(io) do
                EnzymeRates._progress(tmp, true, "stage one")
            end
        end
        @test occursin("stage one", read(out_file, String))
        @test isfile(joinpath(tmp, "progress.log"))
        @test occursin("stage one", read(joinpath(tmp, "progress.log"), String))

        # show_progress=false: writes neither.
        out_file2 = joinpath(tmp, "stdout2.txt")
        open(out_file2, "w") do io
            redirect_stdout(io) do
                EnzymeRates._progress(tmp, false, "silent line")
            end
        end
        @test !occursin("silent line", read(out_file2, String))
        @test !occursin("silent line", read(joinpath(tmp, "progress.log"), String))
    end

    # _batch_summary reports the three buckets with the right denominator.
    mech = first(EnzymeRates.init_mechanisms(@enzyme_reaction begin
        substrates: S[C]; products: P[C] end))
    row = (n_params=3, loss=0.5, mechanism_type="M", rate_equation="v",
           retcode="Success", error=missing, fitted_param_names=(:K,),
           fitted_param_values=(1.0,), eq_hash="abc")
    e_succ = EnzymeRates.BatchEntry(mech, 3, 0.5, :Success, hash(:a), row)
    e_mt   = EnzymeRates.BatchEntry(mech, 3, 0.9, :MaxTime, hash(:b), row)
    f      = EnzymeRates.FitFailure(mech, "StackOverflowError: ")
    s = EnzymeRates._batch_summary([e_succ, e_mt], [f])
    @test occursin("2 fitted", s)
    @test occursin("Success 33.3%", s)            # 1 of 3 total
    @test occursin("non-Success retcode 33.3%", s)
    @test occursin("errored 33.3%", s)
end
```

Also add a progress.log assertion to the existing "CSV output (new schema)" testset (after line 295's `@test "initial_mechanisms.csv" in files`):

```julia
        @test isfile(joinpath(save_dir, "progress.log"))
        @test filesize(joinpath(save_dir, "progress.log")) > 0
```

- [ ] **Step 2: Run to verify failure**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `_progress`/`_batch_summary` undefined; no `progress.log` written.

- [ ] **Step 3: Implement**

In `src/identify_rate_equation.jl`, add these helpers near the top (after the `IdentifyRateEquationResults` struct, before `identify_rate_equation`, around line 87):

```julia
"""
Emit one progress line to flushed stdout AND append it to
`<save_dir>/progress.log`. Gated by `show_progress`. The explicit flush makes
lines appear in a redirected cluster job log (otherwise stdout is block-buffered
when redirected and withholds output until the process exits).
"""
function _progress(save_dir::AbstractString, show_progress::Bool, msg::AbstractString)
    show_progress || return nothing
    println(msg)
    flush(stdout)
    isdir(save_dir) || mkpath(save_dir)
    open(joinpath(save_dir, "progress.log"), "a") do io
        println(io, msg)
    end
    nothing
end

"""
One-line summary of a fitted batch: count + the three retcode/error buckets
(% Success / % non-Success retcode / % errored) and the best loss. Denominator
is fitted + errored mechanisms for the batch.
"""
function _batch_summary(entries::Vector{BatchEntry}, failures::Vector{FitFailure})
    n_fit   = length(entries)
    n_err   = length(failures)
    total   = n_fit + n_err
    n_succ  = count(e -> e.retcode === :Success, entries)
    n_other = n_fit - n_succ
    pct(x)  = total == 0 ? 0.0 : round(100 * x / total; digits=1)
    best    = isempty(entries) ? NaN : minimum(e -> e.loss, entries)
    string(n_fit, " fitted | Success ", pct(n_succ),
           "% | non-Success retcode ", pct(n_other),
           "% | errored ", pct(n_err), "% | best loss ",
           isnan(best) ? "n/a" : round(best; sigdigits=4))
end
```

Add `show_progress::Bool = true` to `identify_rate_equation`'s kwargs (in the "Output & parallelism" group, after `save_dir`, around line 181):

```julia
    save_dir::String = _default_save_dir(),
    show_progress::Bool = true,
    pmap_function::Function = pmap,
```

Pass `show_progress` to `_beam_search` (replace lines 201-206) and capture/announce the `_cv_model_selection` result (replace lines 208-212):

```julia
    mechanisms, df = _beam_search(prob;
        min_beam_width, loss_rel_threshold,
        loss_abs_threshold,
        max_param_count, save_dir, show_progress,
        pmap_function, optimizer, n_cv_candidates,
        fitting_kwargs...)

    result = _cv_model_selection(
        mechanisms, df, prob;
        n_cv_candidates, se_threshold, perm_p_threshold,
        pmap_function, optimizer, save_dir, show_progress,
        fitting_kwargs...)
    _progress(save_dir, show_progress, "Done. Results saved to $save_dir")
    return result
```

Add `show_progress` to `_beam_search`'s signature (replace line 411-415's signature head):

```julia
function _beam_search(
    prob::IdentifyRateEquationProblem;
    min_beam_width, loss_rel_threshold, loss_abs_threshold,
    max_param_count, save_dir, show_progress, pmap_function,
    optimizer, n_cv_candidates, kwargs...
)
```

Add the base-tier progress lines. In the base-tier block from Task 6, insert the enumerate line before `base = ...` and the summary line after `_save_initial_csv(...)`:

```julia
    # ── Base tier: fit ALL init mechanisms (no bucketing — siblings) ──
    _progress(save_dir, show_progress, "Enumerating initial mechanisms…")
    base = _dedup_flat!(collect(init_mechanisms(prob.reaction)))
    _progress(save_dir, show_progress,
        "Fitting $(length(base)) initial mechanisms…")
    base_entries, base_failures = _process_batch(base, prob;
        pmap_function, optimizer, max_param_count, kwargs...)
```

…and after the `_save_initial_csv(...)` call in that block:

```julia
    _progress(save_dir, show_progress,
        "Base tier: " * _batch_summary(base_entries, base_failures))
```

Add the per-iteration progress line in the expansion block from Task 6 (inside the `if !isempty(child_entries) || !isempty(child_failures)` branch, after `_save_iteration_csv(...)`):

```julia
                _progress(save_dir, show_progress,
                    "Iteration $iteration (target n_params=$target): " *
                    "$(length(parents)) parents → $(length(children)) children | " *
                    _batch_summary(child_entries, child_failures))
```

Add `save_dir`/`show_progress` to `_cv_model_selection`'s signature (replace lines 733-740):

```julia
function _cv_model_selection(
    mechs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function, optimizer,
    se_threshold::Float64,
    perm_p_threshold::Float64,
    save_dir, show_progress,
    kwargs...
)
```

Emit the CV stage line just before the `fold_scores_per_candidate = pmap_function(...)` call (around line 767):

```julia
    _progress(save_dir, show_progress,
        "Cross-validating $(length(candidate_mechs)) candidate equations (LOOCV)…")
    fold_scores_per_candidate = pmap_function(
```

Emit the selection line just before `return IdentifyRateEquationResults(...)` (around line 824):

```julia
    _progress(save_dir, show_progress,
        "Selected: $(string(typeof(best_mechanism))), n_params=$(sel.best_n)")
    select!(cv_df, Not(:cv_fold_scores))
    return IdentifyRateEquationResults(best_mechanism, cv_df)
```

- [ ] **Step 4: Run to verify pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — `_progress` testset green; progress.log assertion in the integration test green; full suite green.

- [ ] **Step 5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add cluster-visible progress output (flushed stdout + progress.log)" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: Doc sync + final full-suite verification

**Files:**
- Modify: `.claude/CLAUDE.md`

- [ ] **Step 1: Update CLAUDE.md**

In the `src/fitting.jl` Source-Layout bullet, append:

```markdown
`FittingProblem` carries `scale_k_to_kcat::Union{Real,Nothing}`: a `Real`
(default `1.0`) gives relative data (per-group-centered `loss!` + rescale SS k's
to that kcat); `nothing` gives absolute per-enzyme turnover (uncentered loss, raw
k's). `fit_rate_equation` reads it from `fp` (no kwarg) and returns `(params,
loss, retcode)`.
```

In the `src/rate_eq_derivation.jl` Source-Layout bullet (and the Vmax section's `rescale_parameter_values(m, params; kcat=1.0)` reference), change the kwarg name to `scale_k_to_kcat` (keep `_kcat_forward`/`_kcat_components` as-is — those name the computed turnover, not a scale target).

In the `src/identify_rate_equation.jl` Source-Layout bullet, append:

```markdown
`IdentifyRateEquationProblem` carries `scale_k_to_kcat` (threaded into every
`FittingProblem`). `_process_batch` returns `(entries, failures)` —
`FitFailure`s carry the captured exception text; an all-failed base tier
re-raises it. `_loocv` is loud (raises on a fold exception or non-finite fold
loss). `identify_rate_equation` takes `show_progress::Bool=true`; `_progress`
writes flushed stdout + `<save_dir>/progress.log`. CSVs + `cv_results` carry
`retcode` and `error` columns; CSVs also contain failure rows.
```

- [ ] **Step 2: Run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — entire suite green, including the `rate_equation` perf gate (alloc-free, <100 ns), compile-budget gate, Aqua, and JET (no stale deps, no new exports, derivation untouched).

- [ ] **Step 3: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "Doc sync: scale_k_to_kcat, retcode/error columns, loud failures, progress.log" \
  -m "Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review notes (for the implementer)

- **Return-type coupling (Tasks 5 ↔ 6):** changing `_process_batch` to `(entries, failures)` breaks `_beam_search` until Task 6 lands. Run them back-to-back; the suite is only green again after Task 6.
- **`fit_rate_equation` is accessed by field everywhere** (`.params`, `.loss`) — adding `.retcode` to the returned NamedTuple is non-breaking except for the kcat-normalization test (migrated in Task 2).
- **Failure rows carry `missing`** in `n_params`/`eq_hash`/etc.; the migrated CSV-output assertions use `skipmissing`. Failure rows never enter the frontier/`cv_pool` (only `entries` are ingested), so beam/CV logic never sees `missing`.
- **`_select_best_n_params` keeps its empty-fold handling** — it is unit-tested directly with empty/partial fold-score vectors (lines 718-757). Only `_cv_model_selection`'s `isempty(v) ? Inf` branch + all-non-finite backstop are removed (made dead by the loud `_loocv`).
