# Log the fit plan before fitting, results after

Date: 2026-07-18
Status: Approved design, ready for planning
Branch: `parallelize-seed-and-expansion` (rides the same PR as ideas 1 & 6)

## Problem

Each search batch prints its whole summary only after fitting finishes:

```
Iteration 4 (child n_params 8-11): 1003 parents → 49683 children
  49405 new fits + 252 inherited + 0 skipped (already fit) + 0 skipped (>13 params) + 26 skipped (>337 complexity) + 0 errored | Success 100.0% | non-Success retcode 0.0%
  best loss by n_params: …
```

The parent count, child count, n_params range, and the new/inherited/skip
counts are all known *before* fitting — the slowest step. Printing them only
afterward means an HPC run shows nothing for the whole fitting window, then
dumps everything at once. Only the errored count and the success/retcode
percentages actually need the fit results.

## Goal

Print what is known before fitting, before fitting; print the rest after.

```
Iteration 4 (child n_params 8-11): 1003 parents → 49683 children          ← before fitting
  49405 new fits + 252 inherited + 0 skipped (already fit) + 0 skipped (>13 params) + 26 skipped (>337 complexity)   ← before fitting
  0 errored | Success 100.0% | non-Success retcode 0.0%                    ← after fitting
  best loss by n_params: 5:0.01751 6:0.0139 …                              ← after fitting
```

The base tier gets the same split. This changes only what prints and when; it
does not change what is enumerated or fit.

### "new fits" semantics (decided)

`new fits` is printed before fitting, so it counts the distinct new equations
being fit this batch (the fit representatives), not the ones that succeed. The
post-fit `errored` line reports how many of those threw; `new fits = successful +
errored`. In the common zero-error case the numbers read exactly as today.

## Design

`_process_batch` already runs two passes: PASS-1 (`pmap`) compiles, caps, and
deduplicates to pick one fit representative per equation; PASS-2 (`pmap`) fits
the representatives and builds the rows. The pre-fit counts are all PASS-1
outputs; the errored/success counts are PASS-2 outputs. Split the function at
that seam.

### Functions

- **`_compile_batch(mechs, prob; max_param_count, eq_complexity_filter, memo,
  fitted)` → `(compiled, reps, rep_idx, n_fitted_skip)`** — PASS-1: the
  already-fit filter (mutating `fitted`), the compile/cap `pmap`, and
  representative selection. Everything `_process_batch` does before "Pick one
  representative per `eq_hash`".
- **`_fit_batch(compiled, reps, rep_idx, prob, memo; optimizer, parent_of,
  kwargs...)` → `(entries, failures)`** — PASS-2: fit the representatives,
  populate `memo`, build one row per compiled mechanism.
- **`_process_batch`** stays, as a thin wrapper: call `_compile_batch`, then
  `_fit_batch`, and return the current 5-tuple `(entries, failures,
  n_param_skip, n_cx_skip, n_fitted_skip)`. Its existing tests
  (`test/test_identify_rate_equation.jl:1123-1243`) pass unchanged.
- **`_prefit_summary(child_count, n_new, n_inherited, n_param_skip, n_cx_skip,
  n_fitted_skip; max_param_count, eq_complexity_filter)` → String** — the first
  five buckets: `"N new fits + M inherited + a skipped (already fit) + b skipped
  (>P params) + c skipped (>X complexity)"`.
- **`_postfit_summary(entries, failures)` → String** — `"d errored | Success
  X% | non-Success retcode Y%"`, where `d = length(failures)` and the
  percentages are over the fitted set (`entries`), matching today's
  computation.
- **`_batch_summary` is removed** — the two new summaries replace it. Its test
  (`test/test_identify_rate_equation.jl:1089-1123`) is rewritten to cover
  `_prefit_summary` and `_postfit_summary`.

### Pre-fit counts (from PASS-1 output)

- `child_count` = number of mechanisms in the batch.
- n_params range = min/max `n_params` over the compiled `NamedTuple` records.
- `n_new` = `length(reps)` (distinct new equations to fit).
- `n_inherited` = (compiled `NamedTuple` count) − `n_new` (equations reusing a
  memoized or in-batch fit).
- `n_param_skip` = `count(===(nothing), compiled)`;
  `n_cx_skip` = `count(===(:complexity_skip), compiled)`;
  `n_fitted_skip` from `_compile_batch`.

`child_count = n_new + n_inherited + n_param_skip + n_cx_skip + n_fitted_skip +
(compile failures)`; compile failures surface post-fit in `errored`.

### Call sites (both in `_beam_search`)

Base tier (`identify_rate_equation.jl:791-796`) and the iteration block
(`:842-891`) each change from one `_process_batch` call to:

```
compiled, reps, rep_idx, n_fitted_skip = _compile_batch(batch, prob; …)
n_param_skip = count(===(nothing), compiled)
n_cx_skip    = count(===(:complexity_skip), compiled)
# derive n_new, n_inherited, n_params range; print the pre-fit line(s)
entries, failures = _fit_batch(compiled, reps, rep_idx, prob, memo; …)
# ingest; print the post-fit line + best-loss line
```

The iteration's "print only if rows were produced" guard
(`:851`, `!isempty(child_entries) || !isempty(child_failures)`) and its
`iteration += 1` move earlier: whether rows will be produced is known after
PASS-1 (any compiled `NamedTuple`, any compile failure, or any expansion
failure). The pre-fit line prints in the row-producing case; the existing
"Expanded P parents → C children | all skipped" line stays for the no-row case.

## Scope

Out of scope: the parallelization already on this branch (done), and any change
to what is enumerated or fit. `rate_equation` is untouched.

## Testing

1. `_prefit_summary` renders the five buckets correctly for a batch with new,
   inherited, and each skip kind; `_postfit_summary` renders errored + the two
   percentages, including a non-zero-error case. (Rewrite of the current
   `_batch_summary` test.)
2. The existing `_process_batch` tests stay green (the wrapper preserves its
   contract), confirming `_compile_batch` + `_fit_batch` compose to the same
   result.
3. Full suite green.
4. A search run's `progress.log` shows the two-line-before / two-line-after
   ordering (a small `identify_rate_equation` run, or inspection of the emitted
   lines).
