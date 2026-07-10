# Best mechanism selection

After the rate-equation search fits all candidate mechanisms, `identify_rate_equation`
selects the simplest one that generalizes: the rate equation with the fewest
parameters that best predicts kinetic data it was never fit on. This is the
[bias–variance tradeoff](https://en.wikipedia.org/wiki/Bias%E2%80%93variance_tradeoff)
— adding parameters always improves the fit to the data in hand, but past a
point the extra freedom tracks noise rather than signal and predicts new
measurements worse, while too few parameters miss real structure.
Cross-validation finds the balance by scoring each candidate on data held out
from fitting.

Selection happens in two phases. **While** the search runs, a *beam* keeps only
the most promising mechanisms at each parameter count and expands those into the
next count. **After** the search, leave-one-group-out cross-validation and a
two-test rule pick the single best among the fitted candidates, guarding against
an over-simple model whose better CV score is within noise.

## Running the search

[The enumeration engine](@ref) generates every candidate the beam fits.
`init_mechanisms` builds the minimum-parameter starting mechanisms, and
`expand_mechanisms` grows each survivor through single moves —
splitting a shared rate constant, flipping a rapid-equilibrium step to steady
state, adding a regulator, or making the enzyme allosteric. Every move yields
slightly more complex children, so the beam meets candidates in order of
increasing complexity.

The search walks parameter counts in ascending order. At each count it fits
every candidate, keeps a *beam* of the most promising, and expands only those
survivors into the next count. A mechanism at parameter count `n` stays in the
beam if **either** condition holds:

- **Loss cutoff** — its loss is within

  ```
  min(loss_rel_threshold * best(n) + loss_abs_threshold,
      loss_parsimony_threshold * best(<n))
  ```

  where `best(n)` is the lowest loss at parameter count `n` and `best(<n)` is
  the lowest loss over *all* smaller counts — an added parameter must beat the
  best simpler model of any size. The `best(<n)` term is dropped at the base
  count (no smaller level exists yet).
- **Width floor** — a cumulative per-count budget. The floor keeps expanding
  mechanisms at a parameter count until `min_beam_width` of them have been
  expanded *over the whole search*, then stops; after that only the loss cutoff
  admits at that count. The budget is spent once, not re-granted each time the
  count is revisited.

The four knobs, all tunable kwargs of `identify_rate_equation`:

| Keyword | Default | Meaning |
|---------|---------|---------|
| `loss_rel_threshold` | `2.0` | Relative tolerance: keep losses within this factor of `best(n)`, the best at the same parameter count. |
| `loss_abs_threshold` | `0.01` | Additive tolerance: guards against `best(n)` approaching zero (simulated / very-low-loss data), where a purely multiplicative cutoff would collapse the beam to the single best mechanism. |
| `loss_parsimony_threshold` | `1.01` | An added parameter must earn its keep: keep expanding only if the loss is within this factor of `best(<n)`, the best model of any smaller parameter count. Set to `Inf` to disable it. |
| `min_beam_width` | `50` | Cumulative floor: expand at least this many mechanisms per parameter count over the whole search. Spent once, then only the loss cutoff admits at that count. |

The search ends when the frontier empties — when no surviving mechanism yields a
new candidate to fit. Three limits drain it: `max_param_count` drops any child
with more fitted parameters than the cap, `eq_complexity_filter` drops any child
whose rate equation exceeds the complexity bound (both before fitting), and the
enumeration itself runs dry — every move on every survivor reaches only
mechanisms already fit or already dropped. The complexity bound roughly counts
the terms in a rate equation's denominator; its default (`337`) passes a fully
steady-state random-order bi-bi and skips anything denser. The beam then hands
its pool of fitted candidates to the cross-validation rule below.

## Leave-one-group-out cross-validation

Each unique value in the `:group` column is one fold. A `group` is one
independent experiment — rates measured at a single enzyme amount across varying
metabolite concentrations. For each fold the mechanism is fit on all other
groups and scored on the held-out group, with the fold score floored at
`eps(Float64)` so the log stays finite even when a single-row group fits
exactly. Leaving one group out this way estimates how well the mechanism
predicts a new experiment — a new enzyme batch, a different dilution, or a
distinct assay plate.

The CV score for a mechanism is the **mean of the log per-fold losses**:

```julia
cv_score = mean(log.(fold_scores))
```

This is not the mean of raw fold losses. Working in log space puts all fold
scores on a comparable scale regardless of the absolute rate magnitudes and
penalises extreme fold misses more gracefully than a linear mean would.

## The cross-validation selection rule

Let `n_min` be the parameter count with the lowest mean log CV score (with a
parsimony tiebreak to the smaller count). The selection rule then checks each
simpler bucket (ascending `n_params < n_min`). A simpler bucket is accepted
only if it passes **both**:

1. **Paired 1-SE rule** [Hastie2009](@cite): the mean of the paired
   log-fold-loss differences (`simpler − n_min`) must not exceed
   `se_threshold × std(diffs) / √n_folds`. The default `se_threshold=1.0`
   is the textbook one-standard-error rule.

2. **One-sided sign-flip permutation test**: the p-value `Pr(perm_mean ≥
   observed)` under the sign-flip null must exceed `perm_p_threshold`. The
   default `perm_p_threshold=0.16` matches the paired 1-SE criterion: under a
   normal approximation a ±1-SE band covers 68.3% of the distribution, leaving
   1 − 0.683 = 0.317 in the two tails and ≈ 0.16 in the single tail the
   one-sided test uses. The test asks how often a random sign-flip of the
   per-fold loss differences would look at least as favorable to the simpler
   model; a high p-value means the simpler model is not meaningfully worse.

The function returns the **smallest** `n_params` that passes both tests, or
`n_min` if none pass. When there is only one fold (one group), the SE is
undefined and the loop is skipped, so `n_min` is always returned. Within the
chosen bucket, the mechanism with the lowest **training loss** wins — CV scores
rank buckets, training loss resolves ties inside a bucket.

Both `se_threshold` and `perm_p_threshold` are tunable kwargs of
`identify_rate_equation`.

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
