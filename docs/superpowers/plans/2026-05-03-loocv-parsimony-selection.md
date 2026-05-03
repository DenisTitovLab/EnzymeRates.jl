# LOOCV Parsimony Selection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the strict-argmin model selection in `_cv_model_selection` with a parsimony-aware rule that combines the 1-SE rule and a Wilcoxon signed-rank test, both in log-loss space, taking the more parsimonious of the two — matching the proven approach from the predecessor package `DataDrivenEnzymeRateEqs.jl`.

**Architecture:**

1. `_loocv` returns `Vector{Float64}` of per-fold scores (was: scalar mean). Each score is floored at `eps(Float64)` so `log(score)` is always well-defined.
2. `_cv_model_selection` builds an internal DataFrame with both `cv_score::Float64` (mean) and `cv_fold_scores::Vector{Float64}` (raw per-fold). Per-bucket representative = row with lowest `cv_score` (= lowest mean fold-loss).
3. Selection picks one representative per `n_params` bucket. Both 1-SE and Wilcoxon operate on representatives' per-fold scores so Wilcoxon's pairing assumption is preserved (same fold = same data, paired correctly).
4. `_cv_fold_scores` is dropped before returning `IdentifyRateEquationResults` so `cv_results` stays CSV-safe.
5. `p_value_threshold::Float64 = 0.4` is the new kwarg on `identify_rate_equation` (matches DDEEqs default — parsimony-permissive).

**Tech Stack:** Julia 1.9+; new dep on `HypothesisTests.jl@0.11` (for `ExactSignedRankTest`, `pvalue`).

---

## File Structure

**Modify:**
- `Project.toml` — add `HypothesisTests` to `[deps]` and `[compat]`
- `src/EnzymeRates.jl` — `using HypothesisTests: ExactSignedRankTest, pvalue`
- `src/identify_rate_equation.jl` — `_loocv` return-shape change + eps floor; new `_find_best_n_params_1se`, `_find_best_n_params_wilcoxon`, `_select_best_n_params`; `_cv_model_selection` rework; `p_value_threshold` kwarg threading; docstring update
- `test/test_identify_rate_equation.jl` — unit tests for the three selection helpers + edge cases + `_loocv` shape test

No new files. Each helper goes inside `identify_rate_equation.jl` next to `_cv_model_selection`.

---

## Verification convention

Each task's red/green steps run a single REPL invocation that includes the test file directly, capturing the relevant testset output:

```bash
julia --project -e 'using Test, EnzymeRates, DataFrames, OptimizationPyCMA
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
include("test/test_identify_rate_equation.jl")' 2>&1 | tail -50
```

The `OptimizationPyCMA` import is required because `test/test_identify_rate_equation.jl:8` does `using OptimizationPyCMA`. If the package isn't in your global env, run `Pkg.add("OptimizationPyCMA")` once first. (`Pkg.test` would handle this automatically but takes ~2× longer due to recompilation.)

After every task's commit, also run `Pkg.test()` once to ensure nothing else regressed.

---

## Task 1: Per-fold LOOCV scores + eps floor

**Goal:** Change `_loocv` to return `Vector{Float64}` (one score per fold). Floor each score at `eps(Float64)` so downstream `log(score)` is finite. On failure, return an empty vector.

**Files:**
- Modify: `src/identify_rate_equation.jl` — function `_loocv` (currently lines 650-690)
- Modify: `src/identify_rate_equation.jl` — caller in `_cv_model_selection` (currently around line 720-732)
- Modify: `test/test_identify_rate_equation.jl` — add a test inside the existing `@testset "identify_rate_equation"` block

- [ ] **Step 1.1: Write the failing test**

Append inside the existing `@testset "identify_rate_equation"` block in `test/test_identify_rate_equation.jl`:

```julia
@testset "_loocv returns per-fold scores, floored at eps" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    init = EnzymeRates.init_mechanisms(rxn)
    m = EnzymeRates.EnzymeMechanism(first(init))

    # 3 groups × 2 rows each so per-fold fits aren't degenerate
    data = DataFrame(
        S    = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        P    = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        Rate = [0.5, 0.8, 1.0, 1.1, 1.2, 1.3],
        group = [1, 1, 2, 2, 3, 3],
    )
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)

    scores = EnzymeRates._loocv(
        m, prob;
        optimizer=PyCMAOpt(),
        n_restarts=2, maxtime=2.0,
        maxiters=500, popsize=40, verbose=-9)

    @test scores isa Vector{Float64}
    # When fitting succeeds, length == n_groups; when it fails
    # `_loocv` returns Float64[]. Either is acceptable shape-wise.
    @test length(scores) ∈ (0, 3)
    # Every reported score must be ≥ eps (the floor) and finite.
    for s in scores
        @test s >= eps(Float64)
        @test isfinite(s)
    end
end
```

- [ ] **Step 1.2: Run test to verify it fails**

```bash
julia --project -e 'using Test, EnzymeRates, DataFrames, OptimizationPyCMA
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
include("test/test_identify_rate_equation.jl")' 2>&1 | tail -30
```
Expected: FAIL because `_loocv` currently returns `Float64`, not `Vector{Float64}`.

- [ ] **Step 1.3: Modify `_loocv` to return per-fold scores with eps floor**

Replace the body of `_loocv` (currently the `function _loocv(...)` ... `end` block defined just before `_cv_model_selection`):

```julia
"""
    _loocv(mechanism, prob; optimizer, kwargs...)

Leave-one-group-out cross-validation. Returns `Vector{Float64}`
of per-fold test losses (one per held-out group). Each score is
floored at `eps(Float64)` so `log(score)` is finite — the
selection rules in `_select_best_n_params` operate in log space.
On any internal failure (compile error, fit failure, etc.),
returns an empty `Float64[]`.
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
            # Floor at eps so log() is finite. test_loss can be
            # exactly 0 when the held-out group has only 1 row
            # (centering kills the residual variance), and when a
            # perfect fit happens to align.
            push!(scores,
                  max(test_loss, eps(Float64)))
        end
    catch e
        @debug("LOOCV failed",
            exception=(e, catch_backtrace()))
        return Float64[]
    end

    scores
end
```

The `mean(scores)` line at the end is removed; `scores` (the vector) is returned directly. Failure returns `Float64[]` instead of `Inf`.

- [ ] **Step 1.4: Update the caller in `_cv_model_selection`**

Find the block that begins:

```julia
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
```

Replace with:

```julia
    # LOOCV each candidate in parallel — each result is the
    # candidate's Vector{Float64} of per-fold scores (or empty
    # on failure).
    fold_scores_per_candidate = pmap_function(
        candidate_specs
    ) do spec
        m = compile_mechanism(spec)
        _loocv(m, prob; optimizer, kwargs...)
    end

    # Build CV results DataFrame. `cv_score` holds the mean
    # (used for sorting/display). `cv_fold_scores` holds the
    # raw per-fold vector — INTERNAL only; dropped before
    # returning to the user (see end of this function) so
    # `IdentifyRateEquationResults.cv_results` stays
    # CSV-serialisable.
    cv_df = copy(candidate_rows)
    cv_df.cv_fold_scores =
        collect(fold_scores_per_candidate)
    cv_df.cv_score = [isempty(v) ? Inf : mean(v)
                      for v in cv_df.cv_fold_scores]
    cv_df.spec_idx = candidate_indices
```

The `cv_score` formula `mean(v)` reproduces the previous scalar return.

- [ ] **Step 1.5: Run test to verify it passes**

Same command as Step 1.2. Expected: PASS.

Run `Pkg.test()` to ensure no regression elsewhere:
```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30
```

- [ ] **Step 1.6: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Return per-fold LOOCV scores; floor at eps for log safety"
```

---

## Task 2: HypothesisTests.jl dependency

**Goal:** Add `HypothesisTests@0.11` (verified version 0.11.7 installs cleanly) and import the two needed symbols.

**Files:**
- Modify: `Project.toml`
- Modify: `src/EnzymeRates.jl`

- [ ] **Step 2.1: Add to Project.toml**

In `[deps]`, add:
```toml
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
```

In `[compat]`, add:
```toml
HypothesisTests = "0.11"
```

- [ ] **Step 2.2: Import the symbols**

In `src/EnzymeRates.jl`, add near the existing `using` lines (group with other `using` declarations):

```julia
using HypothesisTests: ExactSignedRankTest, pvalue
```

- [ ] **Step 2.3: Verify package precompiles**

```bash
julia --project -e 'using Pkg; Pkg.instantiate(); using EnzymeRates; println("OK")'
```
Expected: prints `OK` (after any necessary install).

- [ ] **Step 2.4: Commit**

```bash
git add Project.toml src/EnzymeRates.jl
git commit -m "Add HypothesisTests dep for Wilcoxon signed-rank test"
```

---

## Task 3: 1-SE rule (log-loss space, per-bucket representative)

**Goal:** Implement `_find_best_n_params_1se(cv_df)`. Returns the smallest `n_params` whose representative bucket-mean log-loss is within 1 SE of the lowest representative bucket-mean.

**Per-bucket representative:** within each `n_params` bucket, pick the row with the lowest `cv_score` (= lowest mean fold-loss). Use only that row's `cv_fold_scores` for the bucket. This preserves Wilcoxon pairing in Task 4 (each fold = same held-out group across buckets).

**Edge cases handled:**
- Empty `cv_df` (all rows had failed LOOCV → empty fold-scores) → throw `ErrorException("no finite LOOCV scores in cv_df")`.
- Bucket with `n_folds == 1` → `std` returns `NaN`. Special-case: return `n_min` (no SE gating possible).
- Buckets with empty `cv_fold_scores` → drop those rows before bucket-grouping.

**Files:**
- Modify: `src/identify_rate_equation.jl` (add helper before `_cv_model_selection`)
- Modify: `test/test_identify_rate_equation.jl` (add unit tests)

- [ ] **Step 3.1: Write the failing tests**

Append a new top-level testset (sibling of `@testset "identify_rate_equation"`) to `test/test_identify_rate_equation.jl`:

```julia
@testset "_find_best_n_params_1se" begin
    # Case 1: 3 buckets, bucket-7 best, bucket-5 within 1 SE,
    # bucket-3 way outside.
    # bucket-7: mean=0.115, std≈0.01291, SE≈0.00645
    # threshold ≈ 0.12145
    # bucket-5 mean = 0.1165 (within); bucket-3 mean = 0.5375 (out)
    cv_df = DataFrame(
        n_params       = [3, 5, 7],
        cv_score       = [0.5375, 0.1165, 0.115],
        cv_fold_scores = [
            exp.([0.50, 0.55, 0.52, 0.58]),
            exp.([0.115, 0.117, 0.118, 0.116]),
            exp.([0.10, 0.12, 0.11, 0.13]),
        ],
    )
    @test EnzymeRates._find_best_n_params_1se(cv_df) == 5

    # Case 2: simpler bucket within 1 SE of best — pick simpler.
    cv_df2 = DataFrame(
        n_params       = [3, 5],
        cv_score       = [0.115, 0.115],
        cv_fold_scores = [
            exp.([0.10, 0.12, 0.11, 0.13]),
            exp.([0.10, 0.12, 0.11, 0.13]),
        ],
    )
    @test EnzymeRates._find_best_n_params_1se(cv_df2) == 3

    # Case 3: single bucket → returns it.
    cv_df3 = DataFrame(
        n_params       = [4],
        cv_score       = [0.15],
        cv_fold_scores = [exp.([0.1, 0.2, 0.15])],
    )
    @test EnzymeRates._find_best_n_params_1se(cv_df3) == 4

    # Case 4: multiple rows per bucket — representative is the
    # row with lowest cv_score. Bucket-3 row-A has mean 0.135,
    # row-B has mean 0.115 (the rep). Bucket-7 rep has mean
    # 0.110. Rep-bucket-3 mean (0.115) is within 1 SE
    # (≈0.00645) of rep-bucket-7 mean (0.110)? gap=0.005 ≤
    # 0.00645 ✓ → returns 3.
    cv_df4 = DataFrame(
        n_params       = [3, 3, 7],
        cv_score       = [0.135, 0.115, 0.110],
        cv_fold_scores = [
            exp.([0.13, 0.135, 0.14, 0.135]),  # row-A worse
            exp.([0.11, 0.12, 0.115, 0.115]),  # row-B rep
            exp.([0.105, 0.115, 0.110, 0.110]),  # rep
        ],
    )
    @test EnzymeRates._find_best_n_params_1se(cv_df4) == 3

    # Case 5: single-fold bucket → SE undefined → returns n_min.
    cv_df5 = DataFrame(
        n_params       = [3, 5],
        cv_score       = [0.5, 0.1],
        cv_fold_scores = [exp.([0.5]), exp.([0.1])],
    )
    @test EnzymeRates._find_best_n_params_1se(cv_df5) == 5
end
```

Note: every `cv_score` in fixtures equals `mean(log.(cv_fold_scores))`, matching what `_cv_model_selection` would produce after the new `cv_score` formula in Task 1 — but the helper itself only consumes `cv_fold_scores` and `n_params`, not `cv_score`. (We use `cv_score` in the case-4 fixture only to make the row-selection logic explicit for readers.)

Actually, the helper picks the rep per bucket by `cv_score` — so `cv_score` IS load-bearing in case 4. Verified by hand: row-A `cv_score=0.135` > row-B `cv_score=0.115`, so row-B is the rep.

- [ ] **Step 3.2: Run test to verify it fails**

```bash
julia --project -e 'using Test, EnzymeRates, DataFrames, OptimizationPyCMA
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
include("test/test_identify_rate_equation.jl")' 2>&1 | tail -30
```
Expected: FAIL with `UndefVarError: _find_best_n_params_1se`.

- [ ] **Step 3.3: Implement `_find_best_n_params_1se`**

Add to `src/identify_rate_equation.jl` immediately before the `function _cv_model_selection(` line:

```julia
"""
    _find_best_n_params_1se(cv_df) → Int

1-SE rule on log-transformed per-fold LOOCV scores. For each
`n_params` bucket, picks a single representative row (lowest
`cv_score` = lowest mean fold-loss) and uses ONLY that row's
`cv_fold_scores`. This preserves the fold-pairing relied on by
`_find_best_n_params_wilcoxon` and avoids deflating SE by
mixing fits from independent mechanisms.

Returns the smallest `n_params` whose representative mean
log-loss is within one standard error of the bucket with the
lowest representative mean log-loss. Standard error is
`std(log_losses_at_min) / sqrt(n_folds)`.

If `n_folds == 1` for the best bucket (SE undefined), returns
`n_min` (no widening possible).

Drops rows whose `cv_fold_scores` is empty (LOOCV-failure rows)
before grouping. Errors if no valid bucket remains.
"""
function _find_best_n_params_1se(cv_df::DataFrame)
    # Drop failed rows
    valid = filter(row -> !isempty(row.cv_fold_scores), cv_df)
    isempty(valid) && error(
        "no finite LOOCV scores in cv_df")
    # Per-bucket representative = row with lowest cv_score
    sorted = sort(valid, [:n_params, :cv_score])
    reps = combine(groupby(sorted, :n_params), first)
    # Compute log-mean per representative
    log_means = Dict{Int, Float64}()
    log_scores = Dict{Int, Vector{Float64}}()
    for row in eachrow(reps)
        ls = log.(row.cv_fold_scores)
        log_means[row.n_params] = mean(ls)
        log_scores[row.n_params] = ls
    end
    n_min = argmin(n -> log_means[n], keys(log_means))
    losses_at_min = log_scores[n_min]
    n_folds = length(losses_at_min)
    n_folds == 1 && return n_min
    se = std(losses_at_min) / sqrt(n_folds)
    threshold = log_means[n_min] + se
    candidates = [n for n in keys(log_means)
                  if n <= n_min && log_means[n] <= threshold]
    minimum(candidates)
end
```

- [ ] **Step 3.4: Run test to verify it passes**

Same command as Step 3.2. Expected: PASS.

`Pkg.test()` to verify no regression.

- [ ] **Step 3.5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add _find_best_n_params_1se (1-SE rule, log-space, per-bucket rep)"
```

---

## Task 4: Wilcoxon signed-rank rule (log-loss space, paired)

**Goal:** Implement `_find_best_n_params_wilcoxon(cv_df, p_threshold)`. Pairs same-fold scores between the representative of each smaller bucket and the n_min representative.

**Pairing requirement:** `ExactSignedRankTest` is paired — `losses_smaller[i]` is paired with `losses_at_min[i]`. Since both come from the same fold (same held-out group, by construction), pairing is meaningful — provided we use single representatives (not concatenations across mechanisms).

**Files:**
- Modify: `src/identify_rate_equation.jl`
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 4.1: Write the failing tests**

Append to `test/test_identify_rate_equation.jl`:

```julia
@testset "_find_best_n_params_wilcoxon" begin
    # Case 1: bucket-3 has mixed-sign small log-loss diffs vs
    # bucket-7 → high p-value → Wilcoxon picks 3.
    # Both samples have 6 paired observations; diffs ≈
    # [+0.005, -0.002, -0.005, +0.005, +0.002, -0.005] →
    # verified p_value = 1.0 (well above 0.4).
    cv_df = DataFrame(
        n_params       = [3, 5, 7],
        cv_score       = [0.124, 1.125, 0.125],
        cv_fold_scores = [
            exp.([0.105, 0.108, 0.115, 0.135, 0.142, 0.145]),
            exp.([1.10, 1.11, 1.12, 1.13, 1.14, 1.15]),
            exp.([0.10, 0.11, 0.12, 0.13, 0.14, 0.15]),
        ],
    )
    @test EnzymeRates._find_best_n_params_wilcoxon(
        cv_df, 0.4) == 3

    # Case 2: bucket-3 is significantly worse (large uniform
    # offset). Diffs are all +1.0 (six positive ranks) →
    # p_value = 0.03125 < 0.4 → Wilcoxon falls through, returns
    # n_min = 5.
    cv_df2 = DataFrame(
        n_params       = [3, 5],
        cv_score       = [1.125, 0.125],
        cv_fold_scores = [
            exp.([1.10, 1.11, 1.12, 1.13, 1.14, 1.15]),
            exp.([0.10, 0.11, 0.12, 0.13, 0.14, 0.15]),
        ],
    )
    @test EnzymeRates._find_best_n_params_wilcoxon(
        cv_df2, 0.4) == 5

    # Case 3: single bucket → returns it.
    cv_df3 = DataFrame(
        n_params       = [4],
        cv_score       = [0.15],
        cv_fold_scores = [exp.([0.1, 0.2, 0.15])],
    )
    @test EnzymeRates._find_best_n_params_wilcoxon(
        cv_df3, 0.4) == 4

    # Case 4: multiple rows per bucket — representative used.
    # Bucket-3 row-A is worse; row-B is the rep with diffs ≈
    # mixed-sign small → p>0.4 → returns 3.
    cv_df4 = DataFrame(
        n_params       = [3, 3, 7],
        cv_score       = [1.125, 0.124, 0.125],
        cv_fold_scores = [
            exp.([1.10, 1.11, 1.12, 1.13, 1.14, 1.15]),
            exp.([0.105, 0.108, 0.115, 0.135, 0.142, 0.145]),
            exp.([0.10, 0.11, 0.12, 0.13, 0.14, 0.15]),
        ],
    )
    @test EnzymeRates._find_best_n_params_wilcoxon(
        cv_df4, 0.4) == 3
end
```

The numerics in Case 1 and Case 4 were verified empirically — `pvalue(ExactSignedRankTest(losses_close, losses_min))` with the given vectors returns `1.0`, and `pvalue(...)` for Case 2's uniform-offset returns `0.03125`. Don't change the numbers without re-verifying.

- [ ] **Step 4.2: Run test to verify it fails**

Same command as Step 3.2. Expected: FAIL.

- [ ] **Step 4.3: Implement `_find_best_n_params_wilcoxon`**

Add to `src/identify_rate_equation.jl` immediately after `_find_best_n_params_1se`:

```julia
"""
    _find_best_n_params_wilcoxon(cv_df, p_threshold) → Int

Wilcoxon signed-rank rule on log-transformed per-fold LOOCV
scores. For each `n_params` bucket strictly below the argmin
representative bucket, runs a paired signed-rank test comparing
that bucket's representative per-fold log-losses to the best
bucket's representative log-losses. Returns the smallest
`n_params` whose `pvalue > p_threshold` (NOT significantly
worse). Returns `n_min` if no smaller bucket qualifies.

`p_threshold = 0.4` is the parsimony-permissive default
matching `DataDrivenEnzymeRateEqs.jl`. Lower thresholds (e.g.
0.05) require stronger evidence to accept simpler models.

Pairing semantics: the i-th element of each bucket's per-fold
score vector corresponds to the same held-out group, so
pairing `losses_smaller[i]` with `losses_at_min[i]` is
meaningful. Per-bucket representatives (lowest cv_score row)
ensure both vectors come from a single mechanism.

Skips comparisons where fold-counts differ between the two
representatives.
"""
function _find_best_n_params_wilcoxon(
    cv_df::DataFrame, p_threshold::Float64,
)
    valid = filter(row -> !isempty(row.cv_fold_scores), cv_df)
    isempty(valid) && error(
        "no finite LOOCV scores in cv_df")
    sorted = sort(valid, [:n_params, :cv_score])
    reps = combine(groupby(sorted, :n_params), first)

    log_means = Dict{Int, Float64}()
    log_scores = Dict{Int, Vector{Float64}}()
    for row in eachrow(reps)
        ls = log.(row.cv_fold_scores)
        log_means[row.n_params] = mean(ls)
        log_scores[row.n_params] = ls
    end
    n_min = argmin(n -> log_means[n], keys(log_means))
    losses_at_min = log_scores[n_min]
    n_folds = length(losses_at_min)

    smaller_ns = sort([n for n in keys(log_means) if n < n_min])
    for n in smaller_ns
        losses = log_scores[n]
        length(losses) == n_folds || continue
        try
            p = pvalue(ExactSignedRankTest(
                losses, losses_at_min))
            p > p_threshold && return n
        catch e
            @debug("Wilcoxon test failed for n_params=$n",
                   exception=(e, catch_backtrace()))
            continue
        end
    end
    n_min
end
```

- [ ] **Step 4.4: Run test to verify it passes**

Same command as Step 3.2. Expected: PASS. Then `Pkg.test()`.

- [ ] **Step 4.5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add _find_best_n_params_wilcoxon (paired signed-rank, log-space)"
```

---

## Task 5: Combined `_select_best_n_params`

**Goal:** Combine both methods — return `min(n_se, n_wilcoxon)`.

**Files:**
- Modify: `src/identify_rate_equation.jl`
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 5.1: Write the failing test**

Append to `test/test_identify_rate_equation.jl`:

```julia
@testset "_select_best_n_params" begin
    # Wilcoxon picks 3 (mixed-sign diffs, p>0.4); 1-SE picks
    # 7 (bucket-5 mean is far above bucket-7 mean + SE because
    # bucket-5 has high uniform offset). min = 3.
    cv_df = DataFrame(
        n_params       = [3, 5, 7],
        cv_score       = [0.124, 1.125, 0.125],
        cv_fold_scores = [
            exp.([0.105, 0.108, 0.115, 0.135, 0.142, 0.145]),
            exp.([1.10, 1.11, 1.12, 1.13, 1.14, 1.15]),
            exp.([0.10, 0.11, 0.12, 0.13, 0.14, 0.15]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df, 0.4) == 3

    # Both methods agree: bucket-3 == bucket-5 → simpler.
    cv_df2 = DataFrame(
        n_params       = [3, 5],
        cv_score       = [0.115, 0.115],
        cv_fold_scores = [
            exp.([0.10, 0.12, 0.11, 0.13]),
            exp.([0.10, 0.12, 0.11, 0.13]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df2, 0.4) == 3
end
```

- [ ] **Step 5.2: Implement `_select_best_n_params`**

Add to `src/identify_rate_equation.jl` immediately after the two finder helpers:

```julia
"""
    _select_best_n_params(cv_df, p_threshold) → Int

Parsimony-aware model selection. Returns the more parsimonious
of the 1-SE rule and the Wilcoxon signed-rank rule. Both
operate on log-transformed per-fold LOOCV scores via per-bucket
representative selection (lowest cv_score row in each bucket).

Replaces strict-`argmin` over CV means: tiny CV improvements
from added parameters no longer justify the complexity bump.
"""
function _select_best_n_params(
    cv_df::DataFrame, p_threshold::Float64,
)
    n_se = _find_best_n_params_1se(cv_df)
    n_wx = _find_best_n_params_wilcoxon(
        cv_df, p_threshold)
    min(n_se, n_wx)
end
```

- [ ] **Step 5.3: Run test to verify it passes**

Same command as Step 3.2. Expected: PASS. Then `Pkg.test()`.

- [ ] **Step 5.4: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add _select_best_n_params combiner (min of 1-SE and Wilcoxon)"
```

---

## Task 6: Wire into pipeline + drop internal column + kwarg

**Goal:** Replace strict-`argmin` block in `_cv_model_selection` with `_select_best_n_params`. Drop `cv_fold_scores` from the returned `cv_df` so `cv_results` is CSV-safe. Add `p_value_threshold::Float64 = 0.4` kwarg to `identify_rate_equation` and thread through. Update docstring.

**Files:**
- Modify: `src/identify_rate_equation.jl`

- [ ] **Step 6.1: Replace the strict-argmin selection block**

In `_cv_model_selection`, find the block that begins:

```julia
    # Best param count by CV score
    best_cv_per_pc = combine(
        groupby(cv_df, :n_params),
        :cv_score => minimum => :best_cv)
    best_pc_row = best_cv_per_pc[
        argmin(best_cv_per_pc.best_cv), :]
    best_param_count = best_pc_row.n_params
```

Replace with:

```julia
    # Parsimony-aware selection: 1-SE rule + Wilcoxon
    # signed-rank test, both in log-loss space; take the more
    # parsimonious answer.
    best_param_count = _select_best_n_params(
        cv_df, p_value_threshold)
```

- [ ] **Step 6.2: Drop `cv_fold_scores` before returning**

In `_cv_model_selection`, find the line:

```julia
    select!(cv_df, Not(:spec_idx))
```

Replace with:

```julia
    # Drop internal columns so user-facing cv_results stays
    # CSV-serialisable. cv_fold_scores is a Vector column —
    # CSV.write would stringify each cell.
    select!(cv_df, Not([:spec_idx, :cv_fold_scores]))
```

- [ ] **Step 6.3: Add `p_value_threshold` kwarg to `_cv_model_selection`**

Find the signature:

```julia
function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function,
    optimizer, kwargs...
)
```

Replace with:

```julia
function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function,
    optimizer, p_value_threshold, kwargs...
)
```

- [ ] **Step 6.4: Add `p_value_threshold` kwarg to `identify_rate_equation`**

Find the signature of `identify_rate_equation` (a few lines below the docstring near the top of the same file). It looks like:

```julia
function identify_rate_equation(
    prob::IdentifyRateEquationProblem;
    min_beam_width::Int = 50,
    loss_rel_threshold::Float64 = 2.0,
    loss_abs_threshold::Float64 = 0.01,
    max_param_count::Int = 20,
    optimizer,
    n_restarts::Int = 10,
    maxtime::Real = 60.0,
    n_cv_candidates::Int = 5,
    save_dir = nothing,
    pmap_function::Function = pmap,
    kwargs...
)
```

(Exact ordering and names are as in the current file. The "kwargs..." slurps the remaining kwargs forwarded to `fit_rate_equation` / `Optimization.solve` — preserve it.)

Add `p_value_threshold::Float64 = 0.4` between `n_cv_candidates` and `save_dir`:

```julia
    n_cv_candidates::Int = 5,
    p_value_threshold::Float64 = 0.4,
    save_dir = nothing,
```

Find the call site that invokes `_cv_model_selection`. It looks like:

```julia
    return _cv_model_selection(
        all_specs, all_rows_df, prob;
        n_cv_candidates, pmap_function,
        optimizer, kwargs...)
```

(Exact form may vary — check around the end of the function body.)

Add `p_value_threshold` to the kwargs:

```julia
    return _cv_model_selection(
        all_specs, all_rows_df, prob;
        n_cv_candidates, p_value_threshold,
        pmap_function, optimizer, kwargs...)
```

- [ ] **Step 6.5: Update the `identify_rate_equation` docstring**

Find the `# Keyword Arguments` section in the docstring. Insert (between `n_cv_candidates` and `save_dir` entries):

```
- `p_value_threshold::Float64 = 0.4`: parsimony threshold for
  the Wilcoxon signed-rank test in model selection. Smaller
  values demand stronger evidence to accept simpler models.
  Default 0.4 matches DataDrivenEnzymeRateEqs.jl convention
  (parsimony-permissive).
```

Add a new section after the existing `# Beam selection` section:

```
# Model selection (LOOCV)

The best `n_params` is chosen by the more parsimonious of two
methods, both operating on log of per-fold LOOCV scores:

1. **1-SE rule**: smallest `n_params` whose representative
   bucket-mean log-loss is within one standard error of the
   lowest representative bucket-mean log-loss.
2. **Wilcoxon signed-rank test**: smallest `n_params` whose
   representative per-fold log-losses are NOT statistically
   significantly worse than the best bucket
   (`pvalue > p_value_threshold`).

Per-bucket representative = the row with the lowest mean
fold-loss in that `n_params` bucket. Final pick is
`min(n_1se, n_wilcoxon)`. Within the chosen `n_params`, the
mechanism with lowest training loss wins.
```

- [ ] **Step 6.6: Run the full suite**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -50
```
Expected: all tests pass.

If the existing e2e LDH test asserts a specific best mechanism that the new selection no longer picks, update its assertion to assert the new (more parsimonious) outcome — explicitly note that strict-`argmin` was the old behavior.

- [ ] **Step 6.7: Commit**

```bash
git add src/identify_rate_equation.jl
git commit -m "Replace strict-argmin with parsimony-aware LOOCV selection"
```

---

## Task 7: Edge case tests

**Goal:** Cover edge cases not yet exercised: empty cv_df, mixed (some rows fail, some succeed), bucket with only failed rows.

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 7.1: Add edge case tests**

Append to `test/test_identify_rate_equation.jl`:

```julia
@testset "selection helpers: edge cases" begin
    # All rows have empty fold-scores → error.
    cv_df_empty = DataFrame(
        n_params       = [3, 5],
        cv_score       = [Inf, Inf],
        cv_fold_scores = [Float64[], Float64[]],
    )
    @test_throws ErrorException EnzymeRates._find_best_n_params_1se(
        cv_df_empty)
    @test_throws ErrorException EnzymeRates._find_best_n_params_wilcoxon(
        cv_df_empty, 0.4)
    @test_throws ErrorException EnzymeRates._select_best_n_params(
        cv_df_empty, 0.4)

    # Mixed: bucket-3 has one failed row, bucket-5 has a valid
    # row. Bucket-3 dropped → only bucket-5 remains → returns 5.
    cv_df_mixed = DataFrame(
        n_params       = [3, 5],
        cv_score       = [Inf, 0.115],
        cv_fold_scores = [Float64[],
                          exp.([0.10, 0.12, 0.11, 0.13])],
    )
    @test EnzymeRates._find_best_n_params_1se(
        cv_df_mixed) == 5
    @test EnzymeRates._find_best_n_params_wilcoxon(
        cv_df_mixed, 0.4) == 5
    @test EnzymeRates._select_best_n_params(
        cv_df_mixed, 0.4) == 5

    # Bucket-3 has TWO rows: one failed, one valid. The valid
    # row is the rep. Bucket should be retained → returns 3
    # via 1-SE (rep means equal).
    cv_df_partial = DataFrame(
        n_params       = [3, 3, 5],
        cv_score       = [Inf, 0.115, 0.115],
        cv_fold_scores = [Float64[],
                          exp.([0.10, 0.12, 0.11, 0.13]),
                          exp.([0.10, 0.12, 0.11, 0.13])],
    )
    @test EnzymeRates._find_best_n_params_1se(
        cv_df_partial) == 3

    # Wilcoxon length-mismatch: smaller bucket has 4 folds,
    # n_min has 6 → length check skips smaller bucket → returns
    # n_min.
    cv_df_mismatch = DataFrame(
        n_params       = [3, 7],
        cv_score       = [0.124, 0.125],
        cv_fold_scores = [
            exp.([0.105, 0.108, 0.115, 0.135]),       # 4 folds
            exp.([0.10, 0.11, 0.12, 0.13, 0.14, 0.15]),  # 6
        ],
    )
    @test EnzymeRates._find_best_n_params_wilcoxon(
        cv_df_mismatch, 0.4) == 7
end
```

- [ ] **Step 7.2: Run tests**

Same command as Step 3.2. Expected: PASS. Then `Pkg.test()`.

- [ ] **Step 7.3: Commit**

```bash
git add test/test_identify_rate_equation.jl
git commit -m "Edge case tests for LOOCV selection helpers"
```

---

## Task 8: Final cleanup and verification

- [ ] **Step 8.1: Re-read the diff**

```bash
git log --oneline main..HEAD
git diff main..HEAD --stat
git diff main..HEAD -- src/identify_rate_equation.jl | head -200
```

Look for:
- Dead code (e.g. an unreferenced helper).
- Stale references to `cv_score => minimum => :best_cv` or `argmin(best_cv_per_pc.best_cv)` that should have been removed.
- Mentions of `cv_fold_scores` outside the internal selection path (the column shouldn't appear in user-facing returns).
- Docstring claims that don't match the code (e.g. mentioning a default that drifted).

- [ ] **Step 8.2: Run the full suite**

```bash
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -30
```

All tests must pass.

- [ ] **Step 8.3: Smoke-test the new kwarg via the public API**

```bash
julia --project -e '
using EnzymeRates
m = first(methods(identify_rate_equation))
println("identify_rate_equation kwarg names:")
println(Base.kwarg_decl(m))
@assert :p_value_threshold in Base.kwarg_decl(m)
println("OK")
'
```
Expected: prints kwarg list and `OK`.

- [ ] **Step 8.4: Final commit if anything cleaned up**

```bash
git status --short
# If any changes:
git add -p
git commit -m "Cleanup: ..."
```

---

## Notes for the implementer

- **`_loocv` API change is breaking.** The only internal caller is `_cv_model_selection` (rewritten in Task 1.4). There are no external callers — the function is `_`-prefixed.
- **`cv_results` schema is preserved.** Existing columns are unchanged; the internal `cv_fold_scores` column is added during selection then dropped. Users see exactly the same columns as before.
- **Per-bucket representative semantics.** When `n_cv_candidates > 1` and a bucket has multiple rows, only the lowest-mean-fold-loss row is used for selection. This is by design — DDEEqs does the same — to preserve Wilcoxon's pairing assumption.
- **Floor at `eps(Float64)`.** This prevents `log(0)` issues when held-out groups have ≤1 row (the centered-residual loss collapses to 0). The floor is applied in `_loocv`, so all downstream callers see strictly-positive scores.
- **`HypothesisTests@0.11`.** Verified against 0.11.7. `ExactSignedRankTest` and `pvalue` are stable in this major version.
- **`argmin(f, iter)` requires Julia 1.7+.** Project compat is 1.9. Safe.
- **Determinism.** `Dict` iteration order is non-deterministic in Julia, but `argmin(f, keys(log_means))` is deterministic (it returns the unique argmin; ties only arise with exact float equality, which is vanishingly rare). If exact ties prove problematic, switch to `argmin(n -> (log_means[n], n), keys(log_means))` for a deterministic-by-`n` tiebreak.
- **Per-task verification.** Use the focused REPL command for fast unit test iteration; run `Pkg.test()` at the end of each task to catch regressions.
