# Paired 1-SE + Permutation CV Model Selection

## Problem

The current model-selection code in `src/identify_rate_equation.jl` (lines 643–913)
combines two flawed tests with a `min()` combiner:

1. **Unpaired 1-SE rule** (`_find_best_n_params_1se`) computes
   `std(log_folds_best) / sqrt(n_folds)` — the SE of the best model's marginal
   fold variance. Group-difficulty variance (some held-out groups are inherently
   harder) dominates this SE, inflating it past any plausible model difference.
   On the LDH dataset: marginal SE = 0.31, total spread across param counts =
   0.25 — every smaller model passes.

2. **Wilcoxon signed-rank** (`_find_best_n_params_wilcoxon`) discards magnitudes,
   testing rank-based median equality. Inconsistent with `cv_score` (a mean) and
   under-powered when one model is better on most folds but by small amounts vs
   another that wins fewer folds by larger amounts.

3. **`min()` combiner** picks the more parsimonious of the two test outputs, so
   the most aggressive test always wins. On LDH: 1-SE picks n=5, Wilcoxon picks
   n=6, final = n=5, even though n=7 is unambiguously the best by mean log-loss.

The selection block also weighs in at ~270 lines, with three near-parallel
helpers (`_per_bucket_log_stats`, `_find_best_n_params_1se`,
`_find_best_n_params_wilcoxon`, `_select_best_n_params`) and a `_BucketDiagnostic`
scaffolding pattern that exists only to feed two-rule duplication.

## Solution Overview

Replace the two tests with a single rule combining a **paired 1-SE check** and a
**one-sided sign-flip permutation test**, both operating on paired log-loss
differences `d_i = log(fold_loss_n[i]) − log(fold_loss_n_min[i])` between a
candidate's representative LOOCV fold scores and the best bucket's. Accept the
simpler bucket only if **both** tests agree the difference is plausibly
explained by noise.

Pairing matters: fold-to-fold variance is dominated by group difficulty, which
cancels in the per-fold differences. Paired SE on LDH: 0.15 (n=6 vs n=7), 0.23
(n=5 vs n=7) — 2× tighter than the marginal SE, enough to actually discriminate.

The unified rule lets us compress the selection pipeline from ~270 lines to
~180, mostly by eliminating duplicate machinery.

## Selection Rule (per smaller-`n` bucket)

For each `n_params` bucket below `n_min` (the bucket with lowest mean log-fold-
loss), compute on its representative row's per-fold log-losses:

```
diffs = log_scores[n] .- log_scores[n_min]
mean_diff = mean(diffs)
se_paired = std(diffs) / sqrt(n_folds)
permutation_p = Pr(perm_mean ≥ mean_diff)   # one-sided, sign-flip null
```

Accept the simpler bucket iff **both**:
- `mean_diff ≤ se_threshold * se_paired`        (default `se_threshold = 1.0`)
- `permutation_p > perm_p_threshold`            (default `perm_p_threshold = 0.16`)

Iterate smaller buckets in ascending `n_params`; return the first that passes.
Fall through to `n_min` if none pass.

**Edge cases**:
- `n_folds == 1` → SE undefined, return `n_min` immediately.
- A smaller bucket's fold-score vector length ≠ `n_min`'s → error (per-bucket
  fold sets must align for pairing to be meaningful).
- `se_paired == 0` (paired diffs all equal): `mean_diff ≤ 0` accepts; `> 0`
  rejects. Permutation p degenerates to 1.0 if all diffs are 0; otherwise to
  1/2ⁿ if all are equal-sign.
- Single-bucket `cv_df` (only one `n_params` value): no smaller buckets to
  consider; return that bucket as both `n_min` and `best_n` with empty
  comparisons.
- Tie in `mean(log_scores)` across buckets: break ties by smallest `n_params`
  (parsimony). Without this rule, `argmin` over a `Dict` is order-dependent.

## Permutation Test

`_onesided_permutation_p(diffs::Vector{Float64};
                         exact_threshold::Int = 20,
                         mc_samples::Int = 10^6,
                         rng = Random.default_rng())`:

- **Exact enumeration when `length(diffs) ≤ exact_threshold`** (default 20,
  ≤ 2²⁰ ≈ 10⁶ sign patterns). Iterate bitmasks 0…2ⁿ−1; for each mask, the sign
  of `diffs[i]` is flipped iff bit `i−1` is set. Count perms with
  `perm_mean ≥ observed`. Return `count / 2^n_folds`.
- **Monte Carlo with `mc_samples` samples otherwise** (default 10⁶). Random
  `Bool` per fold; same count-and-divide.

Both default branches do approximately equal work (~10⁶ inner iterations).
The `exact_threshold` and `mc_samples` kwargs exist primarily for tests — to
force the MC branch on small fixtures and to use seeded RNGs.

## API Change

`identify_rate_equation` keyword arguments:

| Change | Old | New |
|---|---|---|
| Remove | `p_value_threshold::Float64 = 0.4` | — |
| Add | — | `se_threshold::Float64 = 1.0` |
| Add | — | `perm_p_threshold::Float64 = 0.16` |

Breaking rename, no deprecation shim (per CLAUDE.md). The "Model selection
(LOOCV)" docstring section is rewritten to describe the paired 1-SE +
permutation AND-combiner.

`HypothesisTests.jl` is dropped from `Project.toml` (deps + compat) and from
`src/EnzymeRates.jl`'s `using` line.

## Code Structure

Replace the four current selection helpers with two:

```julia
# Pr(perm_mean ≥ observed_mean) under sign-flip null.
# Exact when length(diffs) ≤ exact_threshold; Monte Carlo otherwise.
# kwargs default to (20, 10⁶, default_rng()); exposed for tests.
_onesided_permutation_p(diffs::Vector{Float64};
                        exact_threshold::Int = 20,
                        mc_samples::Int = 10^6,
                        rng = Random.default_rng()) → Float64

# Single-pass: rep-per-bucket selection, paired diagnostics, AND-combiner policy.
# Returns a NamedTuple with:
#   best_n::Int                        — selected n_params
#   n_min::Int                         — bucket with lowest mean log-fold-loss
#   diagnostics::Dict{Int,@NamedTuple{
#       mean_log_loss_diff::Float64,
#       se_paired::Float64,
#       permutation_p::Float64}}       — n_min bucket has all three = 0.0
_select_best_n_params(cv_df::DataFrame;
                      se_threshold::Float64,
                      perm_p_threshold::Float64) → NamedTuple
```

`_per_bucket_log_stats`, `_find_best_n_params_1se`, `_find_best_n_params_wilcoxon`,
and the existing `_select_best_n_params` are removed. The `_BucketDiagnostic`
struct is not introduced; the diagnostics Dict's value type is a NamedTuple.

`_cv_model_selection` consumes the NamedTuple and:
1. Picks the actual best mechanism (lowest training `loss`) within `best_n`.
2. Stamps the diagnostic columns on `cv_df` from the diagnostics Dict.
   *Diagnostics describe the bucket's representative*, so every row in the
   same `n_params` bucket gets the same values regardless of its own fold
   scores. Non-rep rows therefore see `mean_log_loss_diff` that does not
   match `cv_score[row] − cv_score[n_min_rep]`. This is intentional — the
   selection rule operates on the rep, and surfacing the rep's diagnostic
   uniformly across the bucket is the simplest faithful representation.
3. Flattens raw per-fold scores into `cv_fold_<group_label>` columns from the
   working `cv_fold_scores` Vector field on `cv_df` (using the same group
   ordering as `unique(prob.data.group)`).
4. Drops the working `cv_fold_scores` Vector field and `spec_idx` from the
   user-facing DataFrame.

## `cv_results` DataFrame Shape

**Kept columns**: `n_params`, `loss`, `mechanism_type`, `rate_equation`,
`eq_hash`, `fit_inherited_from_estimate`, `cv_score` (mean log-fold-loss), plus
the existing per-mechanism fitted-param columns.

**Removed**: internal `cv_fold_scores` Vector column, `spec_idx`.

**Added — per-fold columns**: one per held-out group, named `cv_fold_<group_label>`
in the order returned by `unique(prob.data.group)`. Cell value is the raw
eps-floored fold loss (same units `_loocv` returns). If a row's LOOCV failed
entirely (empty fold-scores), all per-fold cells = `missing`.

**Added — diagnostic columns**: `mean_log_loss_diff`, `se_paired`,
`permutation_p`. Computed once per bucket from its representative; every row in
the same `n_params` bucket gets the same diagnostic values. The `n_min` bucket
gets `0.0` in all three (sentinel — `permutation_p` is undefined when comparing
the bucket against itself, but `0.0` is unambiguous since a real test could
never produce p ≤ 0).

Group labels become `Symbol`s in column names (`Symbol("cv_fold_$g")`).
Anything that round-trips through `Symbol` works (Int, String, Symbol). Exotic
labels (with `=`, `,`, etc.) get bracketed-access column names; `DataFrames`
and `CSV.write` handle these natively.

## Test Plan

In `test/test_identify_rate_equation.jl`:

**Remove** the existing testsets (currently lines ~521–723):
- `_find_best_n_params_1se`
- `_find_best_n_params_wilcoxon`
- The current `_select_best_n_params` testset
- Cases inside `selection helpers: edge cases` that exercise the removed helpers

**Add** three new testsets:

### 1. `_select_best_n_params: paired SE math`

Fixtures with hand-computed `mean(diffs)` and `std(diffs)/sqrt(n_folds)`.
Assert returned `mean_log_loss_diff` and `se_paired` match within ~1e-12.
Coverage:
- Mixed-sign diffs (typical case).
- Paired diffs all equal (se_paired = 0): `mean_diff ≤ 0` accepts, `> 0` rejects.
- Single-fold (n_folds = 1): returns `n_min` immediately, no comparisons.
- Multi-row bucket: rep is the row with lowest `cv_score`.
- Single-bucket cv_df: returns the bucket as both `n_min` and `best_n`.
- Tie in mean log-fold-loss: `n_min` resolves to the smallest `n_params`.
- Larger-than-best bucket: diagnostics populated, never selected.

### 2. `_onesided_permutation_p: exact vs Monte Carlo`

Small fixture (n_folds = 8 → 256 exact perms). Compare exact mode against MC
mode by passing `exact_threshold = 0` to force MC, with a seeded RNG and 10⁶
samples; assert agreement within `5/√(10⁶) ≈ 0.005`. Trivial cases:
- All positive diffs → p = 1/2ⁿ exactly.
- All zero diffs → p = 1.0.
- MC determinism: two calls with `MersenneTwister(seed)` produce bit-identical
  output (proves the `rng` kwarg is threaded all the way through).

### 3. `_select_best_n_params: AND-combiner truth table`

Four fixtures hitting each cell of {1-SE pass/fail} × {perm pass/fail}. Assert
the simpler bucket is selected only in the (pass, pass) cell; `n_min` in the
other three.

### Edge cases preserved (in renamed/expanded `selection helpers: edge cases`)

- All LOOCV-failed rows → error.
- Partial bucket failure (one row empty, another valid) → bucket retained via
  rep selection.
- **Length mismatch between fold-score vectors → error** (changed from
  Wilcoxon's silent skip).

### Integration check (extends existing `results structure` testset)

- Assert `mean_log_loss_diff`, `se_paired`, `permutation_p` columns exist.
- Assert the n_min bucket's rows have all three == 0.0.
- Assert per-fold columns named `cv_fold_<group>` exist and align with
  `unique(prob.data.group)`.
- CSV roundtrip: `CSV.write(buf, results.cv_results)` → `CSV.read(buf,
  DataFrame)` preserves all columns including `missing` cells. Acceptance
  criterion #7 ("CSV-serializable, no Vector columns") needs an actual test.
- Exotic group labels: separate small fixture with `group = ["a=b", "c,d"]`
  passed through `_cv_model_selection`'s column-flattening step. Verify the
  resulting columns survive `CSV.write` → `CSV.read` round-trip with
  bracketed access.

`_loocv returns per-fold scores` testset is unchanged.

## Line-Count Estimate

Selection block (current 643–913, ~270 lines) → ~180 lines after rewrite.
Reduction comes from:

- `_per_bucket_log_stats` (~30 lines) merged into `_select_best_n_params`.
- `_find_best_n_params_1se` (~20) and `_find_best_n_params_wilcoxon` (~45) removed.
- Combiner `_select_best_n_params` (~20) replaced by ~12-line policy.
- `_cv_model_selection` shaves ~10 lines via `enumerate` instead of `spec_idx`
  aliasing and removal of the `select! Not(:cv_fold_scores)` step.
- 5–10 lines of explanatory prose in `_loocv` and remaining docstrings trimmed.

Inherent complexity preserved:
- `_loocv` per-fold try/catch with eps-floor and NaN-as-failure.
- Candidate-dedup loop in `_cv_model_selection` (unique-eq-hash per bucket).
- All-Inf check on `cv_score` (protects against silent meaningless results).

## Out of Scope

- Changes to `_beam_search` (lines 471–617) — different concern.
- Changes to fold-score units (still raw eps-floored loss; selection internally
  applies `log()`).
- A `verbosity` knob for printing per-bucket diagnostics during selection — the
  diagnostic columns supersede this need.
