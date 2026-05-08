# Paired 1-SE + Permutation CV Selection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the unpaired 1-SE rule + Wilcoxon signed-rank combiner in `_cv_model_selection` with a single paired 1-SE check AND-combined with a one-sided sign-flip permutation test. Compress the selection block from ~270 lines to ~180. Drop `HypothesisTests.jl`.

**Architecture:** Two new internal helpers (`_onesided_permutation_p`, rewritten `_select_best_n_params`) replace four old ones (`_per_bucket_log_stats`, `_find_best_n_params_1se`, `_find_best_n_params_wilcoxon`, old `_select_best_n_params`). The new `_select_best_n_params` does rep-per-bucket selection, paired-diagnostics computation, and AND-combiner policy in a single pass over `cv_df`. `_cv_model_selection` consumes its NamedTuple result, populates new `mean_log_loss_diff` / `se_paired` / `permutation_p` columns plus per-fold `cv_fold_<group>` columns, and writes a CSV-serializable `cv_results` DataFrame.

**Tech Stack:** Julia 1.9+; DataFrames.jl; Statistics (`mean`, `std`); Random (`default_rng`, `MersenneTwister` for tests). No external statistics packages.

---

## File Structure

| File | Purpose | Touched By |
|---|---|---|
| `src/identify_rate_equation.jl` | Selection helpers + `_cv_model_selection` + `identify_rate_equation` API | Tasks 1, 2, 3 |
| `src/EnzymeRates.jl` | Module imports | Task 1 (add `Random`), Task 4 (drop HypothesisTests) |
| `Project.toml` | Dependency manifest | Task 4 |
| `test/test_identify_rate_equation.jl` | Unit + integration tests for the selection block | Tasks 1, 2, 3 |

Layout decisions:
- `_onesided_permutation_p` lives in `identify_rate_equation.jl` next to `_select_best_n_params` — small, related, no reuse elsewhere.
- The new `_select_best_n_params` returns a NamedTuple bundling four pieces (`best_n`, `n_min`, `bucket_fold_scores`, `diagnostics`), avoiding a struct definition for an internal type.
- `_cv_model_selection` consumes that NamedTuple to populate `cv_df`. No new helper function — inline the column flattening (~12 lines) where it's used once.

---

## Task 1: Add `_onesided_permutation_p` (TDD, pure addition)

**Files:**
- Modify: `src/EnzymeRates.jl` (add `using Random`)
- Modify: `src/identify_rate_equation.jl` (add new function below `_loocv`, around line 702)
- Modify: `test/test_identify_rate_equation.jl` (add new testset)

- [ ] **Step 1.1: Add `Random` import**

In `src/EnzymeRates.jl`, change line 26 from:

```julia
using Distributed
```

to:

```julia
using Distributed
using Random
```

(`HypothesisTests` line stays for now — Task 4 removes it.)

- [ ] **Step 1.2: Write the failing test**

Append to `test/test_identify_rate_equation.jl` (just below the existing `_select_best_n_params` testset; line ordering within the file isn't load-bearing):

```julia
@testset "_onesided_permutation_p" begin
    # All-zero diffs: every sign flip yields perm_mean = observed = 0,
    # so count_ge = 2^n. p = 1.0.
    @test EnzymeRates._onesided_permutation_p(
        [0.0, 0.0, 0.0]) == 1.0

    # All-positive equal diffs: only the identity permutation matches
    # observed; every flipped variant gives a smaller mean.
    # p = 1/2^4 = 0.0625.
    @test EnzymeRates._onesided_permutation_p(
        [1.0, 1.0, 1.0, 1.0]) ≈ 1/16

    # All-negative diffs: observed = -1, all flips ≥ -1 → count_ge = 2^n.
    # p = 1.0.
    @test EnzymeRates._onesided_permutation_p(
        [-1.0, -1.0, -1.0]) == 1.0

    # Mixed-sign 8-fold fixture: 256 exact perms, p strictly in (0, 1).
    diffs = [0.10, -0.05, 0.08, -0.02,
             0.06, -0.04, 0.03, -0.01]
    p_exact = EnzymeRates._onesided_permutation_p(diffs)
    @test 0 < p_exact < 1

    # Force Monte Carlo path (exact_threshold=0) on the same diffs;
    # results must agree within sampling SE. With 10^6 samples and
    # p ≈ 0.5, SE on count_ge/N is √(0.25/10^6) ≈ 5e-4.
    p_mc = EnzymeRates._onesided_permutation_p(
        diffs;
        exact_threshold = 0,
        mc_samples = 10^6,
        rng = MersenneTwister(42),
    )
    @test abs(p_exact - p_mc) < 0.01
end
```

Note: this requires `using Random` at the top of `test/test_identify_rate_equation.jl` so `MersenneTwister` resolves. Check whether it's already imported:

```bash
grep -n "^using" /home/denis.linux/.julia/dev/EnzymeRates/test/test_identify_rate_equation.jl
```

If `Random` isn't already imported, add `using Random` to that file's `using` block.

- [ ] **Step 1.3: Run the test to verify it fails**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
julia --project -e 'using Pkg; Pkg.test(test_args=["_onesided_permutation_p"])'
```

Expected: tests fail with `UndefVarError: _onesided_permutation_p not defined` (or similar). The function does not yet exist.

If the test runner doesn't filter by name, run the full test file:

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
julia --project=test -e 'include("test/test_identify_rate_equation.jl")'
```

and look for the `_onesided_permutation_p` testset failures.

- [ ] **Step 1.4: Implement the function**

In `src/identify_rate_equation.jl`, immediately after the closing `end` of `_loocv` (currently line 702) and before the `"""\n Per-bucket setup ...` docstring, insert:

```julia
"""
    _onesided_permutation_p(diffs; exact_threshold=20,
                             mc_samples=10^6,
                             rng=Random.default_rng()) → Float64

One-sided p-value `Pr(perm_mean ≥ observed)` for paired-difference vector
`diffs` under the sign-flip null. Exact enumeration of `2^n` sign patterns
when `length(diffs) ≤ exact_threshold` (default 20 → up to ~10^6 perms);
Monte Carlo with `mc_samples` random sign-flips otherwise.

Both default branches do ~10^6 inner iterations. The `exact_threshold` and
`mc_samples` kwargs are exposed primarily for tests (forcing the MC branch
on small fixtures and using seeded RNGs).

Errors on empty `diffs` (caller's invariant: a bucket comparison always has
at least one fold).
"""
function _onesided_permutation_p(
    diffs::Vector{Float64};
    exact_threshold::Int = 20,
    mc_samples::Int = 10^6,
    rng = Random.default_rng(),
)
    n = length(diffs)
    n == 0 && error("_onesided_permutation_p: empty diffs vector")
    observed = mean(diffs)

    if n <= exact_threshold
        total = 1 << n   # 2^n
        count_ge = 0
        for mask in 0:(total - 1)
            s = 0.0
            @inbounds for i in 1:n
                bit = (mask >> (i - 1)) & 1
                s += bit == 1 ? -diffs[i] : diffs[i]
            end
            count_ge += (s / n >= observed)
        end
        return count_ge / total
    else
        count_ge = 0
        @inbounds for _ in 1:mc_samples
            s = 0.0
            for i in 1:n
                s += rand(rng, Bool) ? -diffs[i] : diffs[i]
            end
            count_ge += (s / n >= observed)
        end
        return count_ge / mc_samples
    end
end
```

- [ ] **Step 1.5: Run tests to verify they pass**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -50
```

Expected: the new `_onesided_permutation_p` testset shows 5 passed (or whatever the actual `@test` count is). All other previously-passing tests still pass. (The Wilcoxon and 1-SE tests still pass because we haven't touched the old code yet.)

- [ ] **Step 1.6: Commit**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
git add src/EnzymeRates.jl src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Add _onesided_permutation_p helper

Sign-flip-null one-sided p-value with exact enumeration up to n=20 and
Monte Carlo (10^6 samples) above. Pure addition; not yet wired into the
selection pipeline.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Replace selection helpers with new `_select_best_n_params`

**Files:**
- Modify: `src/identify_rate_equation.jl` (delete lines 704–824; insert new function; update `_cv_model_selection` lines 826–912; update `identify_rate_equation` kwargs ~line 332 and docstring ~line 314–331)
- Modify: `test/test_identify_rate_equation.jl` (delete the four old testsets at lines ~521–723; add three new testsets)

This task is the heart of the rewrite. Steps are individually small (each is ~3 minutes) but there are several of them.

- [ ] **Step 2.1: Write failing tests for new `_select_best_n_params` (paired SE math)**

Append to `test/test_identify_rate_equation.jl`. The old testsets at lines 521–723 will be deleted in Step 2.4 — leave them alone for now so the test file still compiles and the new tests fail because the new signature doesn't exist yet.

```julia
@testset "_select_best_n_params: paired SE math" begin
    using Statistics  # for std/mean if not already in scope

    # Simple two-bucket case: n=7 best (lowest cv_score). n=5 paired
    # diffs = [+0.5, +0.5, +0.5, +0.5, +0.5, +0.5] (uniform offset).
    # mean_diff = 0.5, std_diff = 0, se_paired = 0. mean_diff > 0
    # → 1-SE rejects. Permutation: all-positive diffs → only the
    # identity perm reproduces observed; p = 1/2^6 = 0.015625 < 0.16
    # → perm rejects. Both fail → best_n = n_min = 7.
    cv_df = DataFrame(
        n_params       = [5, 7],
        cv_score       = [0.6, 0.1],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.6, 0.7, 0.5, 0.6, 0.5, 0.7]),
            exp.([0.1, 0.2, 0.0, 0.1, 0.0, 0.2]),
        ],
    )
    res = EnzymeRates._select_best_n_params(cv_df)
    @test res.n_min == 7
    @test res.best_n == 7
    @test res.diagnostics[7] ===
          (mean_log_loss_diff=0.0, se_paired=0.0,
           permutation_p=0.0)
    d5 = res.diagnostics[5]
    @test d5.mean_log_loss_diff ≈ 0.5
    @test d5.se_paired ≈ 0.0

    # Mixed-sign small diffs: n=7 best, n=5 paired diffs
    # = [0.0, 0.0, 0.03, -0.01]. mean = 0.005, std ≈ 0.01732,
    # se_paired = 0.01732/sqrt(4) = 0.00866. 0.005 ≤ 0.00866 → 1-SE
    # passes. Mixed-sign → permutation_p > 0.16 in 16 flips. Accept.
    cv_df2 = DataFrame(
        n_params       = [5, 7],
        cv_score       = [0.115, 0.110],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.10, 0.12, 0.13, 0.11]),
            exp.([0.10, 0.12, 0.10, 0.12]),
        ],
    )
    res2 = EnzymeRates._select_best_n_params(cv_df2)
    @test res2.n_min == 7
    @test res2.best_n == 5
    @test res2.diagnostics[5].mean_log_loss_diff ≈ 0.005
    @test res2.diagnostics[5].se_paired ≈
          std([0.0, 0.0, 0.03, -0.01]) / sqrt(4)

    # Multi-row bucket: rep is the row with lowest cv_score.
    # Bucket-3: row-A has fold scores giving mean_log = 0.135;
    # row-B has fold scores giving mean_log = 0.115 → row-B is rep.
    cv_df3 = DataFrame(
        n_params       = [3, 3, 7],
        cv_score       = [0.135, 0.115, 0.110],
        loss           = [0.0, 0.0, 0.0],
        cv_fold_scores = [
            exp.([0.13, 0.135, 0.14, 0.135]),
            exp.([0.11, 0.12, 0.115, 0.115]),
            exp.([0.09, 0.13, 0.10, 0.12]),
        ],
    )
    res3 = EnzymeRates._select_best_n_params(cv_df3)
    @test res3.diagnostics[3].mean_log_loss_diff ≈
          mean(log.(exp.([0.11, 0.12, 0.115, 0.115])) .-
               log.(exp.([0.09, 0.13, 0.10, 0.12])))

    # Single-fold case: n_folds_min = 1 → return n_min, no comparisons.
    cv_df4 = DataFrame(
        n_params       = [3, 5],
        cv_score       = [0.5, 0.1],
        loss           = [0.0, 0.0],
        cv_fold_scores = [exp.([0.5]), exp.([0.1])],
    )
    @test EnzymeRates._select_best_n_params(cv_df4).best_n == 5
end

@testset "_select_best_n_params: AND-combiner truth table" begin
    # Cell 1: pass-pass — mixed-sign small diffs, simpler accepted.
    cv_df_pp = DataFrame(
        n_params       = [3, 7],
        cv_score       = [0.115, 0.110],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.10, 0.13, 0.11, 0.12, 0.13, 0.10]),
            exp.([0.10, 0.13, 0.11, 0.12, 0.10, 0.13]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df_pp).best_n == 3

    # Cell 4: fail-fail — uniform large positive diffs, simpler rejected.
    cv_df_ff = DataFrame(
        n_params       = [3, 7],
        cv_score       = [1.5, 0.1],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([1.5, 1.5, 1.5, 1.5, 1.5, 1.5]),
            exp.([0.1, 0.1, 0.1, 0.1, 0.1, 0.1]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df_ff).best_n == 7

    # Cell 2: 1-SE pass + perm fail — crank perm_p_threshold to 0.99
    # so any non-degenerate p fails the perm gate while 1-SE still
    # passes on the mean=0 fixture.
    @test EnzymeRates._select_best_n_params(
        cv_df_pp; perm_p_threshold = 0.99).best_n == 7

    # Cell 3: 1-SE fail + perm pass — needs strictly positive mean
    # with mixed signs and enough variance that perm p stays above
    # threshold. Hand-computed fixture:
    #   n=3 log-folds = [0.105, 0.098, 0.103, 0.099]
    #   n=7 log-folds = [0.100, 0.100, 0.100, 0.100]
    #   diffs = [0.005, -0.002, 0.003, -0.001]
    #   mean = 0.00125; std ≈ 0.003304; se_paired ≈ 0.001652.
    # 1-SE default: 0.00125 ≤ 1.0*0.001652 ✓ pass.
    # Force fail: se_threshold=0.5 → require 0.00125 ≤ 0.000826 ✗.
    # Permutation (16 exact perms): 5 perms have perm_mean ≥ 0.00125
    # → p = 5/16 = 0.3125 > 0.16 ✓ pass.
    cv_df_marginal = DataFrame(
        n_params       = [3, 7],
        cv_score       = [0.10125, 0.100],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.105, 0.098, 0.103, 0.099]),
            exp.([0.100, 0.100, 0.100, 0.100]),
        ],
    )
    @test EnzymeRates._select_best_n_params(
        cv_df_marginal).best_n == 3   # both pass default
    @test EnzymeRates._select_best_n_params(
        cv_df_marginal; se_threshold = 0.5).best_n == 7
end

@testset "_select_best_n_params: edge cases" begin
    # All rows have empty fold scores → error.
    cv_df_empty = DataFrame(
        n_params       = [3, 5],
        cv_score       = [Inf, Inf],
        loss           = [0.0, 0.0],
        cv_fold_scores = [Float64[], Float64[]],
    )
    @test_throws ErrorException EnzymeRates._select_best_n_params(
        cv_df_empty)

    # Length mismatch between buckets → error (was silent skip in
    # the old Wilcoxon path; we now surface it).
    cv_df_mismatch = DataFrame(
        n_params       = [3, 7],
        cv_score       = [0.116, 0.113],
        loss           = [0.0, 0.0],
        cv_fold_scores = [
            exp.([0.105, 0.108, 0.115, 0.135]),                 # 4
            exp.([0.10, 0.105, 0.11, 0.115, 0.12, 0.125]),      # 6
        ],
    )
    @test_throws ErrorException EnzymeRates._select_best_n_params(
        cv_df_mismatch)

    # Partial bucket failure: one row in bucket-3 failed (empty fold
    # scores), another row valid → bucket retained via rep selection.
    cv_df_partial = DataFrame(
        n_params       = [3, 3, 5],
        cv_score       = [Inf, 0.115, 0.115],
        loss           = [0.0, 0.0, 0.0],
        cv_fold_scores = [
            Float64[],
            exp.([0.10, 0.12, 0.11, 0.13]),
            exp.([0.10, 0.12, 0.11, 0.13]),
        ],
    )
    @test EnzymeRates._select_best_n_params(cv_df_partial).best_n == 3
end
```

- [ ] **Step 2.2: Run tests to verify they fail**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -80
```

Expected: the three new testsets all fail with `MethodError` for the new signature `_select_best_n_params(cv_df)` (no `p_threshold` positional argument). Old testsets still pass.

- [ ] **Step 2.3: Delete the old selection helpers**

In `src/identify_rate_equation.jl`, delete lines 704–824 (inclusive). That's the entire block from `"""\nPer-bucket setup shared by ..."""` through the closing `end` of the old `_select_best_n_params(cv_df, p_threshold)`.

To verify the cut is correct, after deleting, the file should go directly from `_loocv`'s closing `end` (and the `_onesided_permutation_p` block added in Task 1) into `_cv_model_selection` at what was line 826.

Reading guide for the cut:
- Line 704–735: `_per_bucket_log_stats` — DELETE
- Line 737–759: `_find_best_n_params_1se` — DELETE
- Line 761–804: `_find_best_n_params_wilcoxon` — DELETE
- Line 806–824: old `_select_best_n_params` — DELETE

- [ ] **Step 2.4: Delete the old testsets that referenced the deleted helpers**

In `test/test_identify_rate_equation.jl`, delete the four old testsets at lines ~521–723 (inclusive). The exact testsets to remove:
- `@testset "_find_best_n_params_1se" begin … end`
- `@testset "_find_best_n_params_wilcoxon" begin … end`
- `@testset "_select_best_n_params" begin … end` (the OLD one — keep the three new testsets you added in Step 2.1)
- `@testset "selection helpers: edge cases" begin … end`

Verify by searching:

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
grep -c "_find_best_n_params_1se\|_find_best_n_params_wilcoxon\|_per_bucket_log_stats" \
    test/test_identify_rate_equation.jl src/identify_rate_equation.jl
```

Expected: `0` references after the cut.

- [ ] **Step 2.5: Implement the new `_select_best_n_params`**

In `src/identify_rate_equation.jl`, in the gap left by Step 2.3 (between `_loocv`/`_onesided_permutation_p` and `_cv_model_selection`), insert:

```julia
"""
    _select_best_n_params(cv_df; se_threshold=1.0,
                          perm_p_threshold=0.16) → NamedTuple

Paired 1-SE rule AND-combined with a one-sided sign-flip permutation test
on log-transformed per-fold LOOCV scores. Returns:

  best_n::Int                — selected `n_params`
  n_min::Int                 — bucket with lowest mean log-fold-loss
  bucket_fold_scores::Dict   — rep's raw eps-floored fold scores per bucket
  diagnostics::Dict{Int, NamedTuple{(:mean_log_loss_diff, :se_paired,
                                     :permutation_p)}}
                              — `n_min` bucket has all three = 0.0

Per-bucket representative = the row with the lowest `cv_score` in that
`n_params` bucket. For each smaller bucket, computes paired diffs vs the
rep of `n_min`'s log-folds, then accepts iff BOTH:

  mean(diffs) ≤ se_threshold * std(diffs)/sqrt(n_folds)   (paired 1-SE)
  permutation_p > perm_p_threshold                         (perm test)

Iterates smaller buckets in ascending `n_params`; returns the first that
passes. Falls through to `n_min` if none pass.

Errors if:
  * no bucket has any non-empty `cv_fold_scores` row;
  * any bucket's fold-score length differs from `n_min`'s.

When `n_folds_min == 1`, returns `n_min` immediately (SE undefined).
"""
function _select_best_n_params(
    cv_df::DataFrame;
    se_threshold::Float64 = 1.0,
    perm_p_threshold::Float64 = 0.16,
)
    valid = filter(row -> !isempty(row.cv_fold_scores), cv_df)
    isempty(valid) && error(
        "no finite LOOCV scores in cv_df")
    sorted = sort(valid, [:n_params, :cv_score])
    reps = combine(groupby(sorted, :n_params), first)

    bucket_fold_scores = Dict(
        row.n_params => collect(row.cv_fold_scores)
        for row in eachrow(reps))
    log_scores = Dict(n => log.(v)
                      for (n, v) in bucket_fold_scores)
    log_means = Dict(n => mean(ls) for (n, ls) in log_scores)
    n_min = argmin(n -> log_means[n], keys(log_means))
    n_folds_min = length(log_scores[n_min])

    DiagT = NamedTuple{
        (:mean_log_loss_diff, :se_paired, :permutation_p),
        Tuple{Float64, Float64, Float64}}
    diagnostics = Dict{Int, DiagT}()
    diagnostics[n_min] = (
        mean_log_loss_diff = 0.0,
        se_paired = 0.0,
        permutation_p = 0.0,
    )

    smaller_ns = sort([n for n in keys(log_means) if n < n_min])
    larger_ns  = sort([n for n in keys(log_means) if n > n_min])

    for n in vcat(smaller_ns, larger_ns)
        ls = log_scores[n]
        length(ls) == n_folds_min || error(
            "fold-count mismatch for n_params=$n: " *
            "got $(length(ls)), expected $n_folds_min " *
            "(n_min=$n_min)")
        diffs = ls .- log_scores[n_min]
        md  = mean(diffs)
        sep = n_folds_min == 1 ? 0.0 :
              std(diffs) / sqrt(n_folds_min)
        p   = _onesided_permutation_p(diffs)
        diagnostics[n] = (
            mean_log_loss_diff = md,
            se_paired = sep,
            permutation_p = p,
        )
    end

    best_n = n_min
    if n_folds_min > 1
        for n in smaller_ns
            d = diagnostics[n]
            if d.mean_log_loss_diff <= se_threshold * d.se_paired &&
               d.permutation_p > perm_p_threshold
                best_n = n
                break
            end
        end
    end

    return (
        best_n = best_n,
        n_min = n_min,
        bucket_fold_scores = bucket_fold_scores,
        diagnostics = diagnostics,
    )
end
```

- [ ] **Step 2.6: Run new selection tests in isolation**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -100
```

Expected: the three new `_select_best_n_params` testsets pass. The full test run will FAIL elsewhere because `_cv_model_selection` still calls the old signature `_select_best_n_params(cv_df, p_value_threshold)` — that's fixed in Step 2.7.

- [ ] **Step 2.7: Update `_cv_model_selection` and the public API**

Replace the entire body of `_cv_model_selection` (currently lines 826–912 — note the offset shifts after Steps 2.3 and 2.5). Replace with:

```julia
function _cv_model_selection(
    specs::Vector, df::DataFrame,
    prob::IdentifyRateEquationProblem;
    n_cv_candidates, pmap_function, optimizer,
    se_threshold::Float64,
    perm_p_threshold::Float64,
    kwargs...
)
    isempty(specs) && error(
        "No mechanisms were successfully " *
        "fitted during beam search")

    # Pick top n_cv_candidates per (n_params, eq_hash) bucket.
    candidate_indices = Int[]
    df_idx = DataFrame(
        spec_idx = 1:nrow(df),
        n_params = df.n_params,
        loss = df.loss,
        eq_hash = df.eq_hash,
    )
    for gdf in groupby(df_idx, :n_params)
        seen_hashes = Set{String}()
        sorted = sort(gdf, :loss)
        for row in eachrow(sorted)
            row.eq_hash in seen_hashes && continue
            push!(seen_hashes, row.eq_hash)
            push!(candidate_indices, row.spec_idx)
            length(seen_hashes) >= n_cv_candidates && break
        end
    end

    candidate_specs = specs[candidate_indices]
    candidate_rows = df[candidate_indices, :]

    fold_scores_per_candidate = pmap_function(
        candidate_specs
    ) do spec
        m = compile_mechanism(spec)
        _loocv(m, prob; optimizer, kwargs...)
    end

    cv_df = copy(candidate_rows)
    cv_df.cv_fold_scores = collect(fold_scores_per_candidate)
    cv_df.cv_score = [isempty(v) ? Inf : mean(log.(v))
                      for v in cv_df.cv_fold_scores]

    all(!isfinite, cv_df.cv_score) && error(
        "All LOOCV scores are non-finite — every fold's fit " *
        "failed. The pipeline cannot select a best mechanism. " *
        "Inspect optimizer settings (n_restarts, maxtime), " *
        "data quality, or compile failures (run with " *
        "ENV[\"JULIA_DEBUG\"] = \"EnzymeRates\").")

    sel = _select_best_n_params(
        cv_df;
        se_threshold = se_threshold,
        perm_p_threshold = perm_p_threshold,
    )

    # Best mechanism = lowest training `loss` within sel.best_n.
    at_best_pc_idx = findall(==(sel.best_n), cv_df.n_params)
    isempty(at_best_pc_idx) && error(
        "internal: best_n=$(sel.best_n) has no rows in cv_df")
    sort_perm = sortperm(cv_df.loss[at_best_pc_idx])
    best_row_idx = at_best_pc_idx[sort_perm[1]]
    best_spec = candidate_specs[best_row_idx]
    best_mechanism = compile_mechanism(best_spec)

    # Populate diagnostic columns. A bucket may be absent from
    # diagnostics only if every row in it had empty fold scores.
    cv_df.mean_log_loss_diff = [
        haskey(sel.diagnostics, n) ?
            sel.diagnostics[n].mean_log_loss_diff : missing
        for n in cv_df.n_params
    ]
    cv_df.se_paired = [
        haskey(sel.diagnostics, n) ?
            sel.diagnostics[n].se_paired : missing
        for n in cv_df.n_params
    ]
    cv_df.permutation_p = [
        haskey(sel.diagnostics, n) ?
            sel.diagnostics[n].permutation_p : missing
        for n in cv_df.n_params
    ]

    # Flatten per-fold scores into one column per held-out group.
    # Group order matches `_loocv`'s iteration over
    # `unique(prob.data.group)`.
    groups = unique(prob.data.group)
    for (i, g) in enumerate(groups)
        col = Symbol("cv_fold_$g")
        cv_df[!, col] = [
            isempty(v) ? missing : v[i]
            for v in cv_df.cv_fold_scores
        ]
    end

    select!(cv_df, Not(:cv_fold_scores))
    return IdentifyRateEquationResults(best_mechanism, cv_df)
end
```

Then update `identify_rate_equation`'s kwargs (currently around lines 332–354). Find:

```julia
    # Model selection
    n_cv_candidates::Int = 5,
    p_value_threshold::Float64 = 0.4,
```

Replace with:

```julia
    # Model selection
    n_cv_candidates::Int = 5,
    se_threshold::Float64 = 1.0,
    perm_p_threshold::Float64 = 0.16,
```

And find the call to `_cv_model_selection` at the bottom of `identify_rate_equation` (currently line 377–380):

```julia
    return _cv_model_selection(
        specs, df, prob;
        n_cv_candidates, p_value_threshold,
        pmap_function, optimizer, fitting_kwargs...)
```

Replace with:

```julia
    return _cv_model_selection(
        specs, df, prob;
        n_cv_candidates, se_threshold, perm_p_threshold,
        pmap_function, optimizer, fitting_kwargs...)
```

Now update the docstring (currently lines 270–331, the `identify_rate_equation` block doc). Find:

```julia
- `p_value_threshold::Float64 = 0.4`: parsimony threshold for
  the Wilcoxon signed-rank test in model selection. Smaller
  values demand stronger evidence to accept simpler models.
  Default 0.4 matches DataDrivenEnzymeRateEqs.jl convention
  (parsimony-permissive).
```

Replace with:

```julia
- `se_threshold::Float64 = 1.0`: paired 1-SE multiplier for
  model selection. Simpler-model bucket accepted iff its mean
  paired log-loss difference vs the best bucket is `≤
  se_threshold * std(diffs)/sqrt(n_folds)`. Default 1.0 is the
  textbook "1-SE rule".
- `perm_p_threshold::Float64 = 0.16`: minimum one-sided
  permutation p-value for model selection. Simpler-model
  bucket accepted iff `p > perm_p_threshold` under the
  sign-flip null. Default 0.16 matches paired 1-SE empirically.
```

And find the model-selection docstring section (currently lines 314–331):

```julia
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

Replace with:

```julia
# Model selection (LOOCV)

For each `n_params` bucket below `n_min` (lowest mean
log-fold-loss) the rule computes paired log-loss differences
between the bucket's representative and `n_min`'s, then
accepts the simpler bucket iff BOTH:

1. **Paired 1-SE rule**: `mean(diffs) ≤ se_threshold *
   std(diffs)/sqrt(n_folds)`.
2. **One-sided permutation test**:
   `permutation_p > perm_p_threshold`, where p is computed by
   exact enumeration when `n_folds ≤ 20` and Monte Carlo
   (10⁶ samples) otherwise.

Returns the smallest passing `n_params`; falls through to
`n_min` if none pass. Within the chosen bucket the mechanism
with lowest training loss wins. Per-bucket representative =
the row with the lowest `cv_score` in that bucket. Diagnostic
columns `mean_log_loss_diff`, `se_paired`, `permutation_p`
are surfaced in `cv_results`; the `n_min` bucket has 0.0 in
all three.
```

- [ ] **Step 2.8: Update integration test for new cv_results columns**

In `test/test_identify_rate_equation.jl`, find the `results structure` testset (currently around line 215). Append after the existing `cv_score` checks:

```julia
        # Diagnostic columns from paired 1-SE + permutation rule.
        @test "mean_log_loss_diff" in
            names(results.cv_results)
        @test "se_paired" in names(results.cv_results)
        @test "permutation_p" in names(results.cv_results)

        # n_min bucket = the bucket with lowest cv_score after rep
        # selection (which equals lowest mean log-fold-loss). Its
        # rows have all three diagnostics = 0.0.
        n_min_val = results.cv_results.n_params[
            argmin(results.cv_results.cv_score)]
        n_min_rows = filter(row -> row.n_params == n_min_val,
                            results.cv_results)
        @test all(==(0.0),
                  n_min_rows.mean_log_loss_diff)
        @test all(==(0.0), n_min_rows.se_paired)
        @test all(==(0.0), n_min_rows.permutation_p)

        # Per-fold columns named by held-out group label.
        groups = unique(prob.data.group)
        for g in groups
            @test Symbol("cv_fold_$g") in
                propertynames(results.cv_results)
        end
```

- [ ] **Step 2.9: Run all tests**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -120
```

Expected: every testset passes. If the integration test (`test_identify_rate_equation.jl`'s top-level `identify_rate_equation` block) now picks a different `best` mechanism than before — that's by design (the new selection rule is more discriminating than the old min-of-tests). The existing assertion `@test best_np <= gen_np` should still hold since the new rule is more conservative toward simpler models, not more aggressive. If it fails, the test fixture's noiseless data may be exposing a real selection-rule difference; investigate before patching.

- [ ] **Step 2.10: Commit**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
git add src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "$(cat <<'EOF'
Replace 1-SE/Wilcoxon combiner with paired 1-SE + permutation

New _select_best_n_params merges rep selection, paired diagnostics, and
the AND-combiner policy in one pass. Replaces _per_bucket_log_stats,
_find_best_n_params_1se, _find_best_n_params_wilcoxon, and the old
_select_best_n_params combiner.

Public API: p_value_threshold replaced by se_threshold (default 1.0) and
perm_p_threshold (default 0.16). cv_results gains mean_log_loss_diff,
se_paired, permutation_p, and per-fold cv_fold_<group> columns.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Drop `HypothesisTests.jl` dependency

**Files:**
- Modify: `src/EnzymeRates.jl` (remove `using HypothesisTests` line)
- Modify: `Project.toml` (remove from `[deps]` and `[compat]`)

- [ ] **Step 3.1: Verify `HypothesisTests` is no longer referenced in source**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
grep -rn "HypothesisTests\|ExactSignedRankTest\|pvalue" src/ test/
```

Expected: zero matches. (If any remain, Task 2 left a reference unswept; fix that first before continuing.)

- [ ] **Step 3.2: Remove the `using` line**

In `src/EnzymeRates.jl`, delete this line (currently line 27):

```julia
using HypothesisTests: ExactSignedRankTest, pvalue
```

- [ ] **Step 3.3: Remove from `Project.toml`**

In `Project.toml`, delete the line in `[deps]`:

```
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
```

And the line in `[compat]`:

```
HypothesisTests = "0.11"
```

- [ ] **Step 3.4: Run full test suite (Aqua's stale-deps check is the gate)**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
julia --project -e 'using Pkg; Pkg.test()' 2>&1 | tail -60
```

Expected: every testset including Aqua passes. If Aqua flags an undeclared dep instead, double-check Step 3.1's grep — something still imports `HypothesisTests`.

- [ ] **Step 3.5: Commit**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates && \
git add src/EnzymeRates.jl Project.toml
git commit -m "$(cat <<'EOF'
Drop HypothesisTests.jl dependency

Wilcoxon signed-rank test was removed from the selection pipeline; no
other call site remained.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Acceptance Criteria

After all three tasks:

1. `src/identify_rate_equation.jl`'s selection block (`_loocv` through end of `_cv_model_selection`) is approximately 180 lines (down from 270).
2. `_find_best_n_params_1se`, `_find_best_n_params_wilcoxon`, `_per_bucket_log_stats` no longer exist in source or tests.
3. `_select_best_n_params(cv_df; se_threshold, perm_p_threshold)` returns a NamedTuple with `best_n`, `n_min`, `bucket_fold_scores`, `diagnostics`.
4. `_onesided_permutation_p(diffs; ...)` returns a Float64 in [0, 1]; switches between exact and Monte Carlo at `length(diffs) == 20`.
5. `identify_rate_equation` accepts `se_threshold::Float64 = 1.0` and `perm_p_threshold::Float64 = 0.16`; rejects `p_value_threshold` (no shim).
6. `IdentifyRateEquationResults.cv_results` includes columns `mean_log_loss_diff`, `se_paired`, `permutation_p` and per-fold columns `cv_fold_<group_label>` for each unique held-out group; the `n_min` bucket has 0.0 in all three diagnostic columns.
7. `cv_results` is CSV-serializable via `CSV.write` (no Vector columns).
8. `HypothesisTests.jl` removed from `Project.toml` and `src/EnzymeRates.jl`.
9. `julia --project -e 'using Pkg; Pkg.test()'` passes (including Aqua and JET).
