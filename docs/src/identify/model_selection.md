# Model selection

After beam search fits all candidate mechanisms, `identify_rate_equation`
selects the simplest one that generalizes. It does this by leave-one-group-out
cross-validation followed by a two-test rule that guards against choosing an
over-simple model whose better CV score is within noise.

## Leave-one-group-out cross-validation

Each unique value in the `:group` column is one fold. A `group` collects
measurements that share the same total enzyme concentration `E_total`, so
leaving one group out estimates how well a mechanism predicts a new
experimental condition — a new enzyme batch, a different dilution, or a
distinct assay plate.

For each fold, `_loocv` (`src/identify_rate_equation.jl`) fits the mechanism
on all other groups and evaluates the test loss on the held-out group. It then
floors each fold score at `eps(Float64)` so the log is finite even when a
single-row held-out group fits exactly.

### The CV score

The CV score for a mechanism is the **mean of the log per-fold losses**:

```julia
cv_score = mean(log.(fold_scores))
```

This is not the mean of raw fold losses. Working in log space puts all fold
scores on a comparable scale regardless of the absolute rate magnitudes and
penalises extreme fold misses more gracefully than a linear mean would.

### Loud CV

A fold fit that throws propagates immediately. A non-finite fold test loss
raises with the name of the held-out group. The pipeline aborts rather than
silently dropping the candidate, because a corrupted CV invalidates model
selection. The saved CSVs contain every fitted mechanism; re-running CV from
them is cheap once the root cause is resolved.

## Model-selection rule

Let `n_min` be the parameter count with the lowest mean log CV score (with a
parsimony tiebreak to the smaller count). `_select_best_n_params`
(`src/identify_rate_equation.jl`) then checks each simpler bucket (ascending
`n_params < n_min`). A simpler bucket is accepted only if it passes **both**:

1. **Paired 1-SE rule** [Hastie2009](@cite): the mean of the paired
   log-fold-loss differences (`simpler − n_min`) must not exceed
   `se_threshold × std(diffs) / √n_folds`. The default `se_threshold=1.0`
   is the textbook one-standard-error rule.

2. **One-sided sign-flip permutation test**: the p-value `Pr(perm_mean ≥
   observed)` under the sign-flip null must exceed `perm_p_threshold`. The
   default `perm_p_threshold=0.16` matches the paired 1-SE criterion
   empirically.

The function returns the **smallest** `n_params` that passes both tests, or
`n_min` if none pass. When there is only one fold (one group), the SE is
undefined and the loop is skipped, so `n_min` is always returned.

Both `se_threshold` and `perm_p_threshold` are tunable kwargs of
`identify_rate_equation`.

### Within-bucket selection

Within the chosen bucket, the mechanism with the lowest **training loss** wins
(`src/identify_rate_equation.jl`). CV scores rank buckets; training loss
resolves ties inside a bucket.

## The permutation test, exactly

For a small paired-difference vector the permutation p-value is computed by
exact enumeration of all `2^n` sign-flips — pure and reproducible:

```jldoctest
julia> using EnzymeRates

julia> EnzymeRates._onesided_permutation_p([0.1, 0.2, 0.3, 0.4])
0.0625
```

With four elements, there are `2^4 = 16` sign patterns; exactly one of them
(the identity) has a mean at or above the observed mean of 0.25, so the
p-value is `1/16 = 0.0625`.

For `n > 20` folds the function switches to a Monte Carlo estimate
(`mc_samples = 10^6` random sign-flips).

`_onesided_permutation_p` and `_select_best_n_params` are internal helpers
shown here only to make the rule concrete. Users call `identify_rate_equation`
and read `results.cv_results`.

## The `cv_results` DataFrame

`IdentifyRateEquationResults` exposes exactly two fields:

- `best::AbstractEnzymeMechanism` — the selected mechanism.
- `cv_results::DataFrame` — LOOCV results for every candidate that entered
  cross-validation.

`cv_results` columns:

| Column | Description |
|--------|-------------|
| `n_params` | Actual fitted-parameter count. |
| `loss` | Training loss (used for within-bucket ranking). |
| `mechanism_type` | Julia type name of the compiled mechanism. |
| `rate_equation` | Full symbolic rate-equation string. |
| `retcode` | Optimizer return code (`"Success"` = converged). |
| `error` | Exception text if the fit errored; otherwise `missing`. |
| `eq_hash` | Hex hash of the comment-stripped rate equation. Two mechanisms with the same `eq_hash` compute the same rate function. |
| one per fitted parameter | Fitted parameter value, or `missing` if the mechanism lacks that parameter. |
| `cv_score` | Mean of log per-fold losses (lower is better). |
| `mean_log_loss_diff` | Mean paired log-fold-loss difference vs the `n_min` bucket's representative. `0.0` for `n_min`. |
| `se_paired` | Paired standard error: `std(diffs) / √n_folds`. `0.0` for `n_min`. |
| `permutation_p` | One-sided sign-flip p-value. `0.0` for `n_min`. |
| `cv_fold_<group>` | Per-fold test loss for held-out group `<group>`, one column per group. |

Only the top `n_cv_candidates` distinct equations per parameter count enter
LOOCV (default `n_cv_candidates=5`, a kwarg of `identify_rate_equation`).
Distinctness is by `eq_hash`.

---

*This page supersedes the README's model-selection description, which
incorrectly described the CV score as a plain mean of losses and selection
as a plain argmin.*
