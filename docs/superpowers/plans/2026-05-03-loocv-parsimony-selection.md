# LOOCV Parsimony Selection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the strict-argmin model selection in `_cv_model_selection` with a parsimony-aware rule that combines the 1-SE rule (in log-loss space) and a Wilcoxon signed-rank test, taking the more parsimonious of the two — matching the proven approach from the predecessor package `DataDrivenEnzymeRateEqs.jl`.

**Architecture:** Modify `_loocv` to return per-fold scores instead of just the mean. Store per-fold scores alongside the mean in `cv_df`. Add a new `_select_best_n_params(cv_df, p_value_threshold)` helper that implements both methods on log-transformed per-fold scores and returns `min(n_se, n_wilcoxon)`. Replace the strict-`argmin` block in `_cv_model_selection` with a call to this helper. Add `p_value_threshold::Float64 = 0.4` kwarg to `identify_rate_equation`.

**Tech Stack:** Julia 1.11+; new dep on `HypothesisTests.jl` (for `ExactSignedRankTest`, `pvalue`).

---

## File Structure

**Modify:**
- `Project.toml` — add `HypothesisTests` to `[deps]` and `[compat]`
- `src/EnzymeRates.jl` — `using HypothesisTests: ExactSignedRankTest, pvalue`
- `src/identify_rate_equation.jl` — `_loocv` return shape; new `_find_best_n_params_1se`, `_find_best_n_params_wilcoxon`, `_select_best_n_params` helpers; `_cv_model_selection` integration; `p_value_threshold` kwarg threading; docstring update
- `test/test_identify_rate_equation.jl` — unit tests for the three selection helpers + edge cases

No new files needed. Each helper goes inside `identify_rate_equation.jl` next to `_cv_model_selection`.

---

## Task 1: Per-fold LOOCV scores

**Goal:** Change `_loocv` to return `Vector{Float64}` (one score per fold) instead of `Float64` (mean). Failure case returns an empty vector.

**Files:**
- Modify: `src/identify_rate_equation.jl:650-690` (function `_loocv`)
- Modify: `src/identify_rate_equation.jl:720-732` (caller in `_cv_model_selection`)

- [ ] **Step 1.1: Write the failing test**

Add to `test/test_identify_rate_equation.jl` inside the existing `@testset "identify_rate_equation"`:

```julia
@testset "_loocv returns per-fold scores" begin
    # Trivial mechanism + tiny dataset with 2 groups:
    # _loocv should return 2 scores (one per held-out group).
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    init = EnzymeRates.init_mechanisms(rxn)
    m = EnzymeRates.compile_mechanism(first(init))

    data = DataFrame(
        S = [1.0, 2.0, 3.0, 4.0],
        P = [0.1, 0.2, 0.3, 0.4],
        rate = [0.5, 0.8, 1.0, 1.1],
        group = [1, 1, 2, 2],
    )
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)

    scores = EnzymeRates._loocv(
        m, prob;
        optimizer=PyCMAOpt(),
        n_restarts=2, maxtime=2.0,
        maxiters=100, popsize=20, verbose=-9)
    @test scores isa Vector{Float64}
    @test length(scores) == 2  # 2 unique groups
    @test all(isfinite, scores)
end
```

- [ ] **Step 1.2: Run test to verify it fails**

```
julia --project -e 'using Pkg; Pkg.test(test_args=["--testset=_loocv returns per-fold scores"])'
```
Expected: FAIL — test asserts `scores isa Vector{Float64}` but current `_loocv` returns `Float64`.

- [ ] **Step 1.3: Modify `_loocv` to return per-fold scores**

Replace the body of `_loocv` (lines 650-690) with:

```julia
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
        return Float64[]
    end

    scores
end
```

The change: return `scores` (the vector) instead of `mean(scores)`; return `Float64[]` instead of `Inf` on failure. The mean is computed downstream by callers.

- [ ] **Step 1.4: Update the caller in `_cv_model_selection`**

Replace lines 720-732 (the `cv_scores = pmap_function(...)` block plus the cv_df construction) with:

```julia
    # LOOCV each candidate in parallel — returns Vector{Float64}
    # of per-fold scores per candidate.
    cv_fold_scores = pmap_function(
        candidate_specs
    ) do spec
        m = compile_mechanism(spec)
        _loocv(m, prob; optimizer, kwargs...)
    end

    # Build CV results DataFrame. `cv_score` holds the mean for
    # display/sort; `cv_fold_scores` holds the raw per-fold vector
    # used by `_select_best_n_params`.
    cv_df = copy(candidate_rows)
    cv_df.cv_fold_scores = collect(cv_fold_scores)
    cv_df.cv_score = [isempty(v) ? Inf : mean(v)
                      for v in cv_df.cv_fold_scores]
    cv_df.spec_idx = candidate_indices
```

- [ ] **Step 1.5: Run test to verify it passes**

Same command as Step 1.2. Expected: PASS.

- [ ] **Step 1.6: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Return per-fold LOOCV scores; expose cv_fold_scores in cv_df"
```

---

## Task 2: HypothesisTests.jl dependency

**Goal:** Add `HypothesisTests` to project deps and import the two needed symbols.

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

(Verify the latest compat-compatible version via `julia --project -e 'using Pkg; Pkg.add("HypothesisTests"); Pkg.status()'` if 0.11 is wrong.)

- [ ] **Step 2.2: Import the symbols**

In `src/EnzymeRates.jl`, add near the existing `using` lines:

```julia
using HypothesisTests: ExactSignedRankTest, pvalue
```

- [ ] **Step 2.3: Verify package precompiles**

```
julia --project -e 'using EnzymeRates; println("OK")'
```
Expected: prints `OK`.

- [ ] **Step 2.4: Commit**

```bash
git add Project.toml src/EnzymeRates.jl
git commit -m "Add HypothesisTests dep for Wilcoxon signed-rank test"
```

---

## Task 3: 1-SE rule (log-loss space)

**Goal:** Implement `_find_best_n_params_1se(cv_df)`. Returns the smallest `n_params` whose mean log-loss is within 1 SE of the lowest mean log-loss.

**Files:**
- Modify: `src/identify_rate_equation.jl` (add helper next to `_cv_model_selection`, around line 690)
- Modify: `test/test_identify_rate_equation.jl` (add unit tests)

- [ ] **Step 3.1: Write the failing tests**

Add a new testset to `test/test_identify_rate_equation.jl`:

```julia
@testset "_find_best_n_params_1se" begin
    # Synthetic cv_df: 3 buckets, n_params = [3, 5, 7].
    # Per-fold scores chosen so:
    #   bucket-7 has lowest mean log-loss
    #   bucket-5's mean is within 1 SE of bucket-7's
    #   bucket-3's mean is far above 1 SE
    # Expected: 1-SE rule picks bucket-5.
    cv_df = DataFrame(
        n_params = [3, 5, 7],
        cv_fold_scores = [
            exp.([0.5, 0.6, 0.55, 0.62]),  # mean log-loss ~0.57
            exp.([0.12, 0.15, 0.13, 0.14]),  # mean log-loss ~0.135
            exp.([0.10, 0.12, 0.11, 0.13]),  # mean log-loss ~0.115
        ],
    )
    @test EnzymeRates._find_best_n_params_1se(cv_df) == 5

    # Simpler is within 1 SE: pick simpler.
    cv_df2 = DataFrame(
        n_params = [3, 5],
        cv_fold_scores = [
            exp.([0.10, 0.11, 0.12, 0.13]),
            exp.([0.10, 0.11, 0.12, 0.13]),
        ],
    )
    @test EnzymeRates._find_best_n_params_1se(cv_df2) == 3

    # Single bucket: returns it.
    cv_df3 = DataFrame(
        n_params = [4],
        cv_fold_scores = [exp.([0.1, 0.2, 0.15])],
    )
    @test EnzymeRates._find_best_n_params_1se(cv_df3) == 4
end
```

- [ ] **Step 3.2: Run test to verify it fails**

```
julia --project -e 'using Pkg; Pkg.test(test_args=["--testset=_find_best_n_params_1se"])'
```
Expected: FAIL with `UndefVarError: _find_best_n_params_1se`.

- [ ] **Step 3.3: Implement `_find_best_n_params_1se`**

Add to `src/identify_rate_equation.jl` immediately before `_cv_model_selection`:

```julia
"""
    _find_best_n_params_1se(cv_df) → Int

1-SE rule (log-loss space). Returns the smallest `n_params` whose
mean log-loss is within one standard error of the bucket with
the lowest mean log-loss. Operates in log space to dampen
outlier folds.

Standard error is `std(log_losses_at_min) / sqrt(n_folds)`.
Only buckets with `n_params ≤ n_min` are considered (more
complex models can't be more parsimonious).

Multiple `cv_df` rows may share the same `n_params`; their
fold-scores are concatenated for that bucket's stats. Empty
fold-score vectors (LOOCV-failure rows) are dropped.
"""
function _find_best_n_params_1se(cv_df::DataFrame)
    by_n = Dict{Int, Vector{Float64}}()
    for row in eachrow(cv_df)
        isempty(row.cv_fold_scores) && continue
        append!(get!(by_n, row.n_params, Float64[]),
                log.(row.cv_fold_scores))
    end
    isempty(by_n) && error(
        "no finite LOOCV scores in cv_df")

    means = Dict(n => mean(v) for (n, v) in by_n)
    n_min = argmin(n -> means[n], keys(means))
    losses_at_min = by_n[n_min]
    se_at_min = std(losses_at_min) /
                sqrt(length(losses_at_min))

    threshold = means[n_min] + se_at_min
    candidates = [n for (n, m) in means
                  if n <= n_min && m <= threshold]
    minimum(candidates)
end
```

Notes:
- `argmin(f, iter)` requires Julia 1.7+; project compat is 1.9+.
- `std` returns `NaN` if `length == 1`; in that case the threshold equals the mean, and only the n_min bucket itself qualifies → `minimum([n_min]) == n_min`. Safe.

- [ ] **Step 3.4: Run test to verify it passes**

Same command as Step 3.2. Expected: PASS.

- [ ] **Step 3.5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add _find_best_n_params_1se helper (1-SE rule, log-space)"
```

---

## Task 4: Wilcoxon signed-rank rule (log-loss space)

**Goal:** Implement `_find_best_n_params_wilcoxon(cv_df, p_threshold)`. Returns the smallest `n_params` whose per-fold log-losses are NOT statistically significantly worse than the best bucket's log-losses.

**Files:**
- Modify: `src/identify_rate_equation.jl` (add helper)
- Modify: `test/test_identify_rate_equation.jl` (add unit tests)

- [ ] **Step 4.1: Write the failing tests**

Append to `test/test_identify_rate_equation.jl`:

```julia
@testset "_find_best_n_params_wilcoxon" begin
    # Synthetic cv_df where bucket-3's per-fold losses are NOT
    # significantly different from bucket-7 (the best), so the
    # Wilcoxon rule picks bucket-3.
    losses = [0.10, 0.11, 0.12, 0.13, 0.14, 0.15]
    cv_df = DataFrame(
        n_params = [3, 5, 7],
        cv_fold_scores = [
            exp.(losses .+ 0.001),  # nearly identical to best
            exp.(losses .+ 1.0),    # much worse
            exp.(losses),           # best
        ],
    )
    @test EnzymeRates._find_best_n_params_wilcoxon(
        cv_df, 0.4) == 3

    # When the simpler buckets ARE significantly worse, return
    # n_min.
    cv_df2 = DataFrame(
        n_params = [3, 5],
        cv_fold_scores = [
            exp.([2.0, 2.1, 2.2, 2.3, 2.4, 2.5]),  # much worse
            exp.([0.1, 0.11, 0.12, 0.13, 0.14, 0.15]),  # best
        ],
    )
    @test EnzymeRates._find_best_n_params_wilcoxon(
        cv_df2, 0.4) == 5

    # Single bucket: returns it.
    cv_df3 = DataFrame(
        n_params = [4],
        cv_fold_scores = [exp.([0.1, 0.2, 0.15])],
    )
    @test EnzymeRates._find_best_n_params_wilcoxon(
        cv_df3, 0.4) == 4
end
```

- [ ] **Step 4.2: Run test to verify it fails**

```
julia --project -e 'using Pkg; Pkg.test(test_args=["--testset=_find_best_n_params_wilcoxon"])'
```
Expected: FAIL with `UndefVarError`.

- [ ] **Step 4.3: Implement `_find_best_n_params_wilcoxon`**

Add to `src/identify_rate_equation.jl` next to `_find_best_n_params_1se`:

```julia
"""
    _find_best_n_params_wilcoxon(cv_df, p_threshold) → Int

Wilcoxon signed-rank rule (log-loss space). For each bucket
strictly below the argmin-mean-log-loss bucket, perform a paired
Wilcoxon signed-rank test comparing per-fold log-losses to the
best bucket's per-fold log-losses. Return the smallest `n_params`
whose p-value exceeds `p_threshold` (i.e. NOT significantly worse
than the best). Defaults to the best bucket if no smaller bucket
qualifies.

`p_threshold = 0.4` is parsimony-permissive: simpler models are
accepted unless their p-value is below 0.4 (strong evidence
they're worse). Smaller thresholds (e.g. 0.05) require stronger
evidence to accept simpler models.

Multiple `cv_df` rows may share the same `n_params`; their
fold-scores are concatenated. The signed-rank test requires
matched-length samples — if rows in one bucket have different
fold-counts, fall back to skipping that bucket.
"""
function _find_best_n_params_wilcoxon(
    cv_df::DataFrame, p_threshold::Float64,
)
    by_n = Dict{Int, Vector{Float64}}()
    for row in eachrow(cv_df)
        isempty(row.cv_fold_scores) && continue
        append!(get!(by_n, row.n_params, Float64[]),
                log.(row.cv_fold_scores))
    end
    isempty(by_n) && error(
        "no finite LOOCV scores in cv_df")

    means = Dict(n => mean(v) for (n, v) in by_n)
    n_min = argmin(n -> means[n], keys(means))
    losses_at_min = by_n[n_min]
    n_folds = length(losses_at_min)

    smaller_ns = sort([n for n in keys(by_n) if n < n_min])
    for n in smaller_ns
        losses = by_n[n]
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

Same command as Step 4.2. Expected: PASS.

- [ ] **Step 4.5: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add _find_best_n_params_wilcoxon helper"
```

---

## Task 5: Combined `_select_best_n_params`

**Goal:** Combine both methods, returning `min(n_se, n_wilcoxon)`.

**Files:**
- Modify: `src/identify_rate_equation.jl`
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 5.1: Write the failing test**

```julia
@testset "_select_best_n_params" begin
    # Wilcoxon picks 3, 1-SE picks 5 → min returns 3.
    losses = [0.10, 0.11, 0.12, 0.13, 0.14, 0.15]
    cv_df = DataFrame(
        n_params = [3, 5, 7],
        cv_fold_scores = [
            exp.(losses .+ 0.001),
            exp.([0.105, 0.115, 0.125, 0.135,
                  0.145, 0.155]),
            exp.(losses),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df, 0.4) == 3

    # Both methods agree on the same simpler bucket.
    cv_df2 = DataFrame(
        n_params = [3, 5],
        cv_fold_scores = [
            exp.([0.10, 0.11, 0.12, 0.13, 0.14, 0.15]),
            exp.([0.10, 0.11, 0.12, 0.13, 0.14, 0.15]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df2, 0.4) == 3
end
```

- [ ] **Step 5.2: Implement `_select_best_n_params`**

Add to `src/identify_rate_equation.jl` after the two finders:

```julia
"""
    _select_best_n_params(cv_df, p_threshold) → Int

Parsimony-aware model selection. Returns the more parsimonious
of the 1-SE rule and the Wilcoxon signed-rank test result. Both
methods operate on log-transformed per-fold LOOCV scores.

This replaces strict-`argmin` over CV means: tiny CV improvements
from added parameters are no longer enough to win.
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

```
julia --project -e 'using Pkg; Pkg.test(test_args=["--testset=_select_best_n_params"])'
```
Expected: PASS.

- [ ] **Step 5.4: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Add _select_best_n_params combiner"
```

---

## Task 6: Wire into `_cv_model_selection`

**Goal:** Replace strict `argmin` block in `_cv_model_selection` with `_select_best_n_params` call. Add `p_value_threshold` kwarg to `identify_rate_equation` and thread through.

**Files:**
- Modify: `src/identify_rate_equation.jl`

- [ ] **Step 6.1: Write a failing integration test**

Add to `test/test_identify_rate_equation.jl` (in the existing e2e testset, or as a new focused test):

```julia
@testset "identify_rate_equation prefers parsimony" begin
    # Synthetic data fit-able by a simple uni-uni mechanism.
    # The pipeline should select a small n_params answer even
    # if more-complex mechanisms achieve marginally better CV.
    # Concretely: assert the chosen mechanism has fewer params
    # than the absolute-CV-best mechanism in cv_results.
    #
    # (Full set-up follows the existing e2e test pattern in
    # `test_identify_rate_equation.jl`. Use a tiny dataset with
    # 3-4 groups and a uni-uni reaction so the run is fast.)
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    # ... seed data so a 3-param fit is good enough; build prob;
    # call identify_rate_equation; assert
    # length(fitted_params(result.best)) <=
    # min n_params bucket whose CV is within tolerance.
    # (Exact dataset construction mirrors the existing e2e pattern.)
end
```

(The exact data construction should mirror the existing e2e test in `test_identify_rate_equation.jl`; the key invariant to assert is that `result.best` has `length(fitted_params(result.best)) <= the n_params bucket selected by `_select_best_n_params`.)

- [ ] **Step 6.2: Run integration test to verify it fails**

It will fail because the current code uses strict-`argmin`.

- [ ] **Step 6.3: Replace strict `argmin` block in `_cv_model_selection`**

In `_cv_model_selection`, replace lines 744-750 (the `best_cv_per_pc` / `argmin` block):

```julia
    # Best param count by CV score
    best_cv_per_pc = combine(
        groupby(cv_df, :n_params),
        :cv_score => minimum => :best_cv)
    best_pc_row = best_cv_per_pc[
        argmin(best_cv_per_pc.best_cv), :]
    best_param_count = best_pc_row.n_params
```

With:

```julia
    # Parsimony-aware model selection: 1-SE rule and Wilcoxon
    # signed-rank test, both in log-loss space, take the more
    # parsimonious answer.
    best_param_count = _select_best_n_params(
        cv_df, p_value_threshold)
```

Update `_cv_model_selection`'s signature to accept `p_value_threshold`:

```julia
function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function,
    optimizer, p_value_threshold, kwargs...
)
```

- [ ] **Step 6.4: Add the kwarg to `identify_rate_equation`**

In `identify_rate_equation` (around line 326 — the kwarg list):

Add `p_value_threshold::Float64 = 0.4` to the kwarg list:

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
    p_value_threshold::Float64 = 0.4,
    save_dir = nothing,
    pmap_function::Function = pmap,
    kwargs...
)
```

(Replace the existing kwarg signature in-place. The exact ordering and surrounding kwargs may differ; preserve the rest.)

Then thread it through to `_cv_model_selection`:

```julia
    return _cv_model_selection(
        all_specs, all_rows_df, prob;
        n_cv_candidates, pmap_function,
        optimizer, p_value_threshold, kwargs...)
```

- [ ] **Step 6.5: Update the docstring**

In the `identify_rate_equation` docstring (around line 270), in the `# Keyword Arguments` section, ADD between `n_cv_candidates` and `save_dir`:

```julia
- `p_value_threshold::Float64 = 0.4`: parsimony threshold for
  the Wilcoxon signed-rank test in model selection. Smaller
  values require stronger evidence to accept simpler models.
  Default 0.4 matches DataDrivenEnzymeRateEqs.jl convention
  (parsimony-permissive).
```

In the `# Beam selection` section (or below), ADD a new section:

```julia
# Model selection (LOOCV)

The best `n_params` is chosen by the more parsimonious of two
methods, both operating on log of per-fold LOOCV scores:

1. **1-SE rule**: smallest `n_params` whose mean log-loss is
   within one standard error of the lowest mean log-loss.
2. **Wilcoxon signed-rank test**: smallest `n_params` whose
   per-fold log-losses are NOT statistically significantly
   worse than the best bucket (`pvalue > p_value_threshold`).

Final pick is `min(n_1se, n_wilcoxon)`. Within the chosen
`n_params`, the mechanism with lowest training loss wins.
```

- [ ] **Step 6.6: Run integration test to verify it passes**

```
julia --project -e 'using Pkg; Pkg.test(test_args=["--testset=identify_rate_equation prefers parsimony"])'
```
Expected: PASS.

- [ ] **Step 6.7: Run the full test suite**

```
julia --project -e 'using Pkg; Pkg.test()'
```
Expected: all tests pass.

If the existing e2e LDH test asserts a specific "best" mechanism that the new selection no longer picks, update its assertion to match the new (more parsimonious) answer.

- [ ] **Step 6.8: Commit**

```bash
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Replace strict-argmin with parsimony-aware LOOCV selection"
```

---

## Task 7: Edge case tests

**Goal:** Cover the documented edge cases in `_find_best_n_params_1se` and `_find_best_n_params_wilcoxon`.

**Files:**
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 7.1: Add edge case tests**

Append to `test/test_identify_rate_equation.jl`:

```julia
@testset "selection helpers: edge cases" begin
    # All fold scores empty → error.
    cv_df_empty = DataFrame(
        n_params = [3, 5],
        cv_fold_scores = [Float64[], Float64[]],
    )
    @test_throws ErrorException EnzymeRates._find_best_n_params_1se(
        cv_df_empty)
    @test_throws ErrorException EnzymeRates._find_best_n_params_wilcoxon(
        cv_df_empty, 0.4)

    # Mixed: one bucket has empty fold-scores → that bucket
    # is dropped, others considered.
    cv_df_mixed = DataFrame(
        n_params = [3, 5],
        cv_fold_scores = [Float64[],
                          exp.([0.1, 0.11, 0.12])],
    )
    @test EnzymeRates._find_best_n_params_1se(
        cv_df_mixed) == 5

    # Single fold per bucket: SE = NaN → 1-SE rule falls
    # back to n_min.
    cv_df_single = DataFrame(
        n_params = [3, 5],
        cv_fold_scores = [
            exp.([0.5]),
            exp.([0.1]),
        ],
    )
    @test EnzymeRates._find_best_n_params_1se(
        cv_df_single) == 5

    # Multiple rows per bucket (same n_params): per-fold scores
    # are concatenated.
    cv_df_multi = DataFrame(
        n_params = [3, 3, 5],
        cv_fold_scores = [
            exp.([0.10, 0.11, 0.12]),
            exp.([0.13, 0.14, 0.15]),
            exp.([0.10, 0.11, 0.12]),
        ],
    )
    # Combined bucket-3 has [0.10..0.15] in log-space —
    # mean ~0.125, vs bucket-5 mean = 0.11 → 1-SE rule
    # depends on bucket-5's SE, but should not error.
    result = EnzymeRates._find_best_n_params_1se(cv_df_multi)
    @test result in (3, 5)
end
```

- [ ] **Step 7.2: Run tests**

```
julia --project -e 'using Pkg; Pkg.test(test_args=["--testset=selection helpers: edge cases"])'
```
Expected: PASS.

- [ ] **Step 7.3: Commit**

```bash
git add test/test_identify_rate_equation.jl
git commit -m "Edge case tests for LOOCV selection helpers"
```

---

## Task 8: Final cleanup and verification

**Goal:** Ensure the diff is clean, the docstring is complete, and all tests pass.

- [ ] **Step 8.1: Re-read the diff with fresh eyes**

```bash
git log --oneline main..HEAD
git diff main..HEAD --stat
```

Look for: dead code, leftover debug prints, missing kwarg propagation, stale docstrings.

- [ ] **Step 8.2: Run the full test suite**

```
julia --project -e 'using Pkg; Pkg.test()'
```

All tests must pass (24,100+ existing tests + new ones).

- [ ] **Step 8.3: Verify package exports**

```
julia --project -e 'using EnzymeRates; @assert :p_value_threshold in
    Base.kwarg_decl.(methods(identify_rate_equation))[1]'
```

Confirms the new kwarg is on the exported function.

- [ ] **Step 8.4: Final commit (if any cleanup applied)**

```bash
git status --short
# If any changes:
git add -p  # review hunks
git commit -m "Cleanup: …"
```

---

## Notes for the implementer

- The new `_loocv` API change (returns `Vector{Float64}`) is breaking for any direct caller. The only internal caller is `_cv_model_selection`. There are no external callers of `_loocv` (it's prefixed with `_`).
- `cv_df.cv_score` is preserved (mean of fold scores) for backward compat with existing CSV output and downstream consumers. The `cv_fold_scores` column is added — extending the schema, not breaking it.
- `HypothesisTests.jl`'s `ExactSignedRankTest` requires non-tied input pairs; ties are resolved by the package automatically (mid-rank handling). Our log-fold-losses can have ties when fits converge to identical log-likelihoods, but this is rare with float scores.
- `argmin(f, iter)` was added in Julia 1.7. Project compat is already 1.9+. Safe.
- The `1-SE rule` historically uses `1*SE`; some literature uses `0.5*SE` or other multipliers. We're matching DDEEqs's `1*SE` choice exactly. Don't make this configurable unless asked.
