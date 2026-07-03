# LDH Search Improvements — Design

## Context

The 2026-07-02 LDH run (`docs/ldh_hpc_results/2026_07_02_results`, driver
`docs/ldh_hpc_results/identify_ldh.jl`) exposed five problems in
`identify_rate_equation`'s beam search, progress reporting, and
cross-validation. Reconstructing the run from its saved CSVs pinned the root
causes:

- **The beam floor is re-granted every sweep.** `min_beam_width=50` keeps the
  top 50 mechanisms per parameter count *per sweep*, not per count over the
  whole run. The advancing-`target` sweep revisits each count many times, so
  ~50 mechanisms per count expand every iteration regardless of loss. The
  parsimony cutoff never prunes.
- **The run wasted 13 iterations.** No parameter count improved its best loss
  after iteration 9, yet iterations 10–22 ran ~21,600 genuine fits that changed
  nothing. The search stopped only when children hit the `max_param_count` cap.
- **"fitted" in the log counts reused fits.** At iteration 9, 34,742 of 48,958
  "fitted" mechanisms inherited a memoized fit and 14,216 were genuine. At
  iteration 22, all 4 "fitted" were inherited.
- **A third of late-iteration children vanish silently.** A child whose fitted
  parameter count exceeds `max_param_count` is dropped before fitting and left
  out of the summary. Iteration 13 dropped 1,413 of its 3,896 children this way.
- **Cross-validation used 45 of ~1,000 cores.** LOOCV parallelizes over 45
  candidates, each running its 18 folds serially, so utilization stalled at
  39–63% for 2.5 hours.
- **RAM grew from ~10 GB to over 50 GB**, risking OOM. The cause is unconfirmed;
  a local profiling run will diagnose it before any fix.

## Scope

Five changes. Changes 1, 2, 3, and 5 are fully specified below. Change 4
(memory) begins with a diagnostic run; its fix follows the findings and will be
specified then. All code changes live in `src/identify_rate_equation.jl` except
the documentation in change 5.

`min_beam_width`'s default stays 50. This design fixes what the floor *means*;
retuning the default is a separate decision, taken after the cumulative version
runs.

---

## 1. Cumulative per-count beam floor

**Problem.** `_select_beam` ranks each sweep's fresh batch on its own
(`rank ≤ min_beam_width`). A count swept N times receives N independent floors
of 50. Because expansion feeds each count new children across many sweeps, the
floor keeps ~50 per count expanding every iteration, and the loss cutoff
(`min(loss_rel_threshold·best(c)+loss_abs_threshold, loss_parsimony_threshold·best(<c))`)
never gets to prune. The search runs to the parameter cap instead of stopping
when it stops improving.

**Change.** Track one integer per parameter count, `expanded[c]` — the number
of count-`c` mechanisms expanded so far — beside the existing per-count best
loss. In the per-count selection loop of `_beam_search`:

1. Pass `max(0, min_beam_width − expanded[c])` to `_select_beam` as the floor.
2. Add the number selected to `expanded[c]`.

A mechanism is admitted iff `loss ≤ cutoff` **or** its rank within this sweep's
count-`c` batch `≤ min_beam_width − expanded[c]`. Once `expanded[c]` reaches
`min_beam_width`, only the loss cutoff admits at that count.

`expanded[c]` counts every admitted mechanism, cutoff-passers included. So
`min_beam_width` means "expand at least this many per count over the whole
search" — the standard beam-width reading. `_select_beam` itself does not
change; the caller manages the budget by passing the remaining floor.

**The advancing-`target` sweep stays verbatim.** Expansion is non-monotone in
parameter count: `Δparams ≥ 0` always, but `Δ = 0` on ~16% of LDH moves (a move
adds structure and a Haldane/Wegscheider constraint absorbs the added
parameter). A count-`c` mechanism can therefore spawn a count-`c` child, and the
re-sweep exists to catch these stragglers. The cumulative floor depends on it:
stragglers draw from the same count's shared budget rather than a fresh 50.
Collapsing the sweep to one pass per count would drop those ~16% of children — a
correctness regression, not a simplification.

**Effect.** Each count's floor budget is spent once. When every count's budget
is spent and no new child passes the loss cutoff, no mechanism is selected, the
frontier drains, and the search ends on convergence. On the LDH data this ends
near iteration 9 rather than 22, cutting ~21,600 pointless fits.

**Tests.**
- Cumulative budget: across two sweeps of the same count with the same
  `min_beam_width`, floor-admitted mechanisms total at most `min_beam_width`;
  the second sweep's floor budget equals `min_beam_width` minus the first
  sweep's admissions.
- Cutoff still admits after the budget is spent: with `expanded[c] ≥
  min_beam_width`, a mechanism under the loss cutoff is still selected.
- Termination: a search over a fixture converges before `max_param_count` when
  later counts stop improving.

---

## 2. Progress-log reporting

**Problem.** The summary line reports one ambiguous "fitted" number (genuine
fits plus inherited copies), hides cap-skipped children entirely, and prints a
single per-batch best loss with no parameter-count context.

**Change — reconcile every child.** Report the four buckets so they sum to the
child count:

```
children = new fits + inherited + skipped (>max_param_count) + errored
```

- `new fits` — entries with `fit_inherited == false` (a genuine optimizer run).
- `inherited` — entries with `fit_inherited == true` (a memoized fit reused).
- `skipped` — children dropped before fitting for exceeding `max_param_count`;
  `length(children) − length(entries) − length(failures)`, computed at the call
  site, which knows the child count.
- `errored` — `failures` (compile, render, or representative fit threw).

`_batch_summary` derives `new`/`inherited` from `entries` (each `BatchEntry`
carries `row.fit_inherited`) and takes `skipped` from the caller. The
`Success` / `non-Success retcode` percentages stay, computed over the fitted
set.

**Change — best loss by parameter count (Option C).** On every iteration line
print the running best loss for each parameter count seen so far, and mark the
counts that improved this iteration:

```
Iteration 8  (child n_params 7-13): 3127 parents → 39501 children
  14595 new fits + 24552 inherited + 342 skipped (>13 params) + 12 errored
  Success 99.8% | non-Success retcode 0.1%
  best loss by n_params: 5:0.01751 ... 12:0.00970 13:0.00932*   (* improved)

Iteration 13 (child n_params 8-13): 300 parents → 3896 children
  1439 new fits + 1043 inherited + 1413 skipped (>13 params) + 1 errored
  Success 99.6% | non-Success retcode 0.4%
  best loss by n_params: 5:0.01751 ... 12:0.00970 13:0.00932   (no improvement)
```

The map reads from `best_loss_by_count`, already maintained. To mark
improvements, snapshot `best_loss_by_count` before `_ingest!` and compare after.

**Reorder.** Emit the progress line *after* `_ingest!`, so the map reflects the
current iteration's children. The line moves below the ingest call; the CSV
save and iteration counter stay where they are.

The base tier gets the same four-bucket line and best-loss map.

**Tests.**
- `_batch_summary` reports `new`/`inherited` split matching the entries'
  `fit_inherited` flags, and the four buckets sum to the child count.
- The best-loss line lists every seen count and stars exactly the counts whose
  best improved this iteration; a no-improvement iteration prints the
  no-improvement marker.

---

## 3. Cross-validation parallelism

**Problem.** `_cv_model_selection` runs `pmap` over the 45 candidate equations,
and each candidate's `_loocv` runs its 18 folds serially on one worker. The grid
holds 45 × 18 = 810 independent fits but exposes only 45 to the scheduler, and
load imbalance drops utilization further.

**Change.** Flatten the parallelism to `(candidate, fold)`.

1. Extract the per-fold body of `_loocv` into
   `_cv_fold_loss(mech, prob, held_out; optimizer, kwargs...)`: subset train and
   test data, fit on train, evaluate the test loss, keep the non-finite check
   (naming the held-out group), and floor at `eps`.
2. In `_cv_model_selection`, build the full `(candidate_index, group)` grid and
   `pmap` over all 810 tasks at once. Each task compiles its candidate and calls
   `_cv_fold_loss`.
3. Reduce the results into one per-fold score vector per candidate, ordered by
   `unique(prob.data.group)`, then feed `_select_best_n_params` unchanged.

`_loocv` remains as a thin serial loop over `_cv_fold_loss`, so its existing
tests keep passing. Loud failures stay loud: a throwing fold or non-finite loss
propagates out of `pmap` and aborts model selection, as today. Per-task
compilation adds no real cost — the `@generated` rate equation is cached per
worker by mechanism signature.

**Tests.**
- Flattened CV reproduces the serial `_loocv` fold scores for a fixture
  (equivalence).
- A fold that throws, and a fold that returns a non-finite loss, each abort
  model selection with the group named.

---

## 4. Memory — diagnose first

**Problem.** Master or worker RAM climbs from ~10 GB to over 50 GB across the
run. `GC.gc()` may not help if the growth is compiled code rather than data.

**Diagnostic run (local).** Run a reduced LDH search locally — lower
`max_param_count`, a handful of workers, the real optimizer — and sample memory
each iteration:

- Master: `Base.gc_live_bytes()` and `Base.summarysize` of `frontier`, `memo`,
  and `cv_pool`.
- Each worker: `Base.gc_live_bytes()` via `remotecall`, to separate data
  retention from per-signature code-cache growth.

**Leading hypotheses (to confirm or reject, not to fix blind).**
- Code-cache growth: `rate_equation` is `@generated` per `EnzymeMechanism{Sig}`,
  so every distinct mechanism compiles specialized code each worker caches and
  never frees. `GC.gc()` would not reclaim it.
- The per-iteration DataFrame: `_rows_to_dataframe` builds one column per
  distinct parameter name, `missing`-filled — iteration 9's CSV is 962 MB, so
  its in-memory form is multiple GB on the master.

The fix is chosen after the run isolates the dominant term, and specified in a
follow-up section of this doc.

---

## 5. Documentation

Align every description of the search filter with the cumulative behavior.
`README.md` only links to the docs, so no filter prose lives there; the targets
are:

- `docs/src/identify/model_selection.md` — the primary description (the "beam"
  section). Correct `best(n-1)` to `best(<c)` (the best over *all* smaller
  counts), and replace "hard floor per parameter count" / per-rank wording with
  the cumulative floor: `min_beam_width` expands at least this many per count
  over the whole search, spent once, not re-granted per sweep.
- The `identify_rate_equation` docstring and the `_select_beam` /
  `_parsimony_cutoff` comments in `src/identify_rate_equation.jl` — same two
  corrections.
- `docs/src/getting_started.md` and `docs/src/identify/tutorial.md` — the
  default-settings sentences, if they describe the floor's behavior.

The generated `docs/build/*` files regenerate from `docs/src`; leave them.

---

## Implementation order

Land the correctness fixes and the docs first; profile last.

1. Change 1 (cumulative floor) — the core correctness fix.
2. Change 2 (reporting) — depends on change 1's `expanded[c]` and the ingest
   reorder.
3. Change 3 (CV parallelism) — independent of 1–2.
4. Change 5 (docs) — describes the shipped cumulative behavior, so it follows
   changes 1–3.
5. Change 4's diagnostic run, then its fix once the dominant term is identified.

Each change follows TDD: a failing test first, then the change.

## Recorded decisions

- Cumulative floor tracked as two numbers per count, `(best_loss, expanded)`;
  no per-count top-K set (a hard budget bounds work at `min_beam_width`, which
  top-K does not).
- The advancing-`target` sweep is load-bearing for non-monotone expansion and
  stays unchanged.
- Best-loss reporting uses the running map with per-iteration improvement marks
  (Option C).
- `min_beam_width` default stays 50; retuning is deferred.
- Memory fix waits on the local diagnostic; no blind `GC.gc()`.
