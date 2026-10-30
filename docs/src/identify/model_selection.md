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
`init_mechanisms` builds the starting mechanisms (or `seed_mechanisms` starts
from fully-regulated ones when a reaction's regulators are required), and
`expand_mechanisms` grows each survivor through single moves — splitting a shared
rate constant, flipping a rapid-equilibrium group to steady state, adding a
regulator, or making the enzyme allosteric. Every move yields slightly more
complex children, so the beam meets candidates in order of increasing complexity.

The beam first fits every starting mechanism, each at its own parameter count.
A degenerate starting mechanism is expanded instead of fitted: its rate at zero
products ignores some substrate or never saturates, or the same holds for
products, so it is not worth a fit, but it is the only parent of mechanisms that
are not degenerate. Every seed with two chemistry isomerizations, such as a
ping-pong seed, is one: a seed holds one steady-state step, so its other
isomerization sits at rapid equilibrium with the steps on both its sides. The
beam replaces such a mechanism at once by its children that are not degenerate
and fits them beside the other starting mechanisms (`Enumeration.starting_mechanisms`). Only a
flip of a rapid-equilibrium group to steady state can cure it, by turning that
isomerization steady state or every step on one of its sides; the other moves
keep which steps are at rapid equilibrium, so their children stay degenerate and
are not fitted. The children's rows, like the starting mechanisms' rows, name no
parent, and the degenerate mechanism gets a row of its own only when its
expansion fails. A degenerate mechanism that a move builds is fitted like any
other.

The beam ranks laws, not mechanisms. Mechanisms whose laws agree up to renaming the fitted
constants share a key, count once at the key's lowest loss, and share one fit, which each
reads under its own names. An L phantom's key is the exception: a member met before the fit
of its law at `L = 0` is fitted itself, and one met after reads that fit, so the key can
carry two losses. A sweep expands every member it meets of a selected key; a member the beam
passed over in an earlier sweep has left the frontier. A key the beam passes over is judged
again when a later sweep meets another member, and at the same loss it stays out.

The search walks parameter counts in ascending order. At each count it fits
every candidate, keeps a *beam* of the most promising, and expands only those
survivors into the next count. A key at parameter count `n` enters the beam if
**either** condition holds:

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
  mechanisms at a parameter count until `min_beam_width` distinct keys have been
  expanded *over the whole search*, then stops; after that only the loss cutoff
  admits at that count. The budget is spent once, not re-granted each time the
  count is revisited.

The four knobs, all tunable kwargs of `identify_rate_equation`:

| Keyword | Default | Meaning |
|---------|---------|---------|
| `loss_rel_threshold` | `1.3` | Relative tolerance: keep losses within this factor of `best(n)`, the best at the same parameter count. |
| `loss_abs_threshold` | `0.001` | Additive tolerance: guards against `best(n)` approaching zero (simulated / very-low-loss data), where a purely multiplicative cutoff would collapse the beam to the single best key. The default `0.001` corresponds to about a 1.03× error — well below experimental error in enzyme kinetics. |
| `loss_parsimony_threshold` | `0.99` | An added parameter must earn its keep: the loss cutoff admits a key only if its loss is at most this factor times `best(<n)`, the best model of any smaller parameter count (default `0.99`: an added parameter must improve on every simpler model by at least 1%). Set to `Inf` to disable it. |
| `min_beam_width` | `50` | Cumulative floor: expand at least this many distinct laws (keys) per parameter count over the whole search. Spent once, then only the loss cutoff admits at that count. |

The search ends when the frontier empties — when no surviving mechanism yields a
new candidate to fit. Three limits drain it: `max_param_count` drops any child
with more fitted parameters than the cap, `eq_complexity_filter` drops any child
whose rate equation exceeds the complexity bound (both before fitting), and the
enumeration itself runs dry — every move on every survivor reaches only
mechanisms already fit or already dropped. The complexity bound roughly counts
the terms in a rate equation's denominator, both conformations together for an
allosteric mechanism; its default (`337`) passes a fully steady-state random-order bi-bi
and skips anything denser, and a one-subunit allosteric variant of a law near it can
exceed it. The beam then hands
its pool of fitted candidates to the cross-validation rule below.

## Leave-one-group-out cross-validation

Each unique value in the `:group` column is one fold. A `group` is one
independent experiment — rates measured at a single enzyme amount across varying
metabolite concentrations. For each fold the mechanism is fit on all other
groups and scored on the held-out group. Leaving one group out this way
estimates how well the mechanism predicts a new experiment — a new enzyme
batch, a different dilution, or a distinct assay plate.

The CV score for a mechanism is the **mean of the per-fold losses**:

```julia
cv_score = mean(fold_scores)
```

## The cross-validation selection rule

Selection follows the one-standard-error rule [Hastie2009](@cite). The
**best equation** is the candidate with the lowest `cv_score` (with a parsimony
tiebreak to the smaller `n_params`). Its fold scores set the yardstick:

```julia
cv_score_se = std(fold_scores) / sqrt(n_folds)
cutoff      = cv_score_best + se_threshold * cv_score_se_best
```

The selected equation is the lowest-`cv_score` candidate at the **smallest**
`n_params` that has a candidate at or below the cutoff. A simpler equation wins
whenever it predicts held-out groups as well as the best equation to within the
best equation's fold-to-fold standard error; an added parameter must lower the
CV score by more than that standard error to be kept. The best equation always
passes its own cutoff, so an equation with more parameters than the best one is
never selected.

The default `se_threshold=1.0` is the textbook rule; it is a tunable kwarg of
`identify_rate_equation`. Larger values favor simpler equations.

## The `cv_results` DataFrame

`IdentifyRateEquationResults` exposes exactly two fields:

- `best::AbstractEnzymeMechanism` — the selected mechanism.
- `cv_results::DataFrame` — LOOCV results for every candidate that entered
  cross-validation.

`cv_results` columns:

| Column | Description |
|--------|-------------|
| `n_params` | Actual fitted-parameter count. |
| `parent_n_params` | Fitted-parameter count of the parent the mechanism was expanded from; `missing` for a starting mechanism. |
| `loss` | Training loss (ranks which equations enter LOOCV). |
| `mechanism` | The mechanism as the `@enzyme_mechanism` or `@allosteric_mechanism` block that builds it, the text `show` prints. `EnzymeRates.DSL.parse_mechanism(text)` reads it back without evaluating it. |
| `mechanism_id` | The first 16 hex digits of the SHA-256 of `mechanism`, the same in every process and Julia version. Read the ids as text, `CSV.read(file, DataFrame; types = Dict(:mechanism_id => String, :parent_id => String, :fit_source_id => String))`: CSV.jl reads an id that looks like a number (decimal digits, or digits around one `e` or `f`) as a number. |
| `parent_id` | `mechanism_id` of the parent the mechanism was expanded from; `missing` for a starting mechanism. |
| `rate_equation` | The printed law in named coefficients (`rate_equation_string(m; coefficients = true)`). |
| `retcode` | Optimizer return code (`"Success"` = converged; `"ReferenceFailed"` = the law predicts zero, non-finite or opposite-sign turnover at a `reference` anchor, so the fit has infinite loss and leaves the search). |
| `error` | Exception text if the mechanism's derivation, compilation, fit or expansion threw; otherwise `missing`. |
| `key` | The SHA-256 (64 hex digits) of the law's canonical form: its terms with the fitted constants labelled canonically and the concentrations by name. Mechanisms share a key exactly when their laws agree up to renaming the fitted constants. A law whose canonical search reaches 100 labellings takes its mechanism's own key, the SHA-256 of `mechanism|` and its `mechanism_id`. `missing` in a failure row. |
| `fit_inherited` | `true` when the row reads another mechanism's fit: one with the same key, under the row's own names, or, for an L phantom, the fit of the law it reduces to; `false` for the mechanism actually fit. |
| `fit_source_id` | `mechanism_id` of the mechanism whose fit the row reads, the row's own when it was fit. |
| `L_phantom` | `true` for a one-subunit allosteric variant whose inactive conformation turns nothing over, when an exact rescaling of its other constants absorbs its `L`, its active conformation's law holds every fitted constant but `L`, and a mechanism with that law was fitted earlier in the search or is fitted in the variant's own iteration: the data cannot determine `L`, so the row reads that fit, with `L = 0`, and keeps its own `n_params`. Such a variant is fitted like any other candidate, `L` included, and its flag is `false`, when no mechanism with that law was fitted earlier in the search or is fitted in its own iteration, or when that mechanism's fit throws. |
| one per fitted parameter | Fitted parameter value, scaled to the `reference` turnovers in relative mode, or `missing` if the mechanism lacks that parameter. |
| `cv_score` | Mean of per-fold losses (lower is better). |
| `cv_score_se` | Standard error of the per-fold losses: `std(fold_scores) / √n_folds`. |
| `cv_fold_<group>` | Per-fold test loss for held-out group `<group>`, one column per group. |

Only the top `n_cv_candidates` distinct equations per parameter count enter
LOOCV (default `n_cv_candidates=5`, at least 1, a kwarg of
`identify_rate_equation`). Distinctness is by `key`: LOOCV fits each
key's source once per fold, and every candidate with that key, an L phantom of it
included, reads those fold losses. A fit whose training loss
is not finite never enters LOOCV, nor is it expanded.

---
