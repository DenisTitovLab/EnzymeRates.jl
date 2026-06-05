# Rate-Scaling Modes, Fit Robustness, and Progress Output — Design

**Date:** 2026-06-05
**Status:** design approved, pending spec review → implementation plan

---

## Goal

Four coupled improvements to the `fit → identify_rate_equation` pipeline. All changes live in
`src/fitting.jl` and `src/identify_rate_equation.jl`; `rate_equation` and the `@generated`
derivation are **not** touched.

1. **`scale_k_to_kcat` (one knob).** A `Union{Real,Nothing}` data property on the problem
   structs: a `Real` (default `1.0`) = relative y-axis (`loss!` centers per group,
   `fit_rate_equation` rescales SS k's so kcat = the value); `nothing` = absolute per-enzyme
   turnover (uncentered loss, raw fitted k's). One field replaces a `:relative`/`:absolute` enum
   plus a separate rescale value, and makes the two nonsensical center/rescale combinations
   unrepresentable. Enables absolute-turnover datasets and removes the arbitrary-scale rate
   constants the old `kcat=nothing` path reported.
2. **Solver outcome recording.** `fit_rate_equation` returns the best restart's `retcode`;
   it flows into the CSVs and `cv_results`. No solver exceptions are swallowed here.
3. **Loud failures, no silent drops.** Replace the `@debug`-swallowing `try/catch` in
   `_process_batch` with capture+count+report; a whole-base-tier failure re-raises the real
   exception (catches config bugs like an unsupported optimizer kwarg). `_loocv` failures
   become loud (a corrupted CV invalidates model selection; CV is cheap to recompute from the
   saved CSVs).
4. **Cluster-visible progress.** Master-level stage logging to flushed stdout **and** an
   appended `<save_dir>/progress.log`, gated by `show_progress`.

## Why

- **Scaling.** `loss!` (`fitting.jl:108`) subtracts each group's mean log-ratio, making the
  loss invariant to per-group `E_total` — so only equation *shape* (K's and ratios of k's) is
  identifiable and the absolute k-scale is free. The `scale_k_to_kcat` rescale then pins the
  scale to a *convention*; with no convention the reported SS rate constants are arbitrary.
  Datasets whose y-axis is genuine per-enzyme turnover (1/time) carry real absolute scale that
  centering throws away. Absolute mode preserves it (and incidentally makes single-group
  datasets fittable — centering degenerates a lone group to zero loss).
- **Silent failures.** `_process_batch` (`identify_rate_equation.jl:338`) catches every
  compile/fit exception → `@debug` → `nothing`, and `_loocv` (`:522`) catches → `Float64[]`.
  A deterministic config error (e.g. an optimizer that rejects `popsize`) therefore drops
  *every* mechanism silently and surfaces only as the downstream "All LOOCV scores non-finite"
  error, with no cause. `fit_rate_equation` never inspects `sol.retcode`, so a budget-truncated
  (non-converged) fit is reported with the same confidence as a converged one.
- **Progress.** There is no progress output anywhere in `src/`. A long run (especially a
  cluster batch job) is opaque. The classic gotcha Denis hit before: stdout is line-buffered on
  a TTY but **block**-buffered when redirected to a job log, so progress lines never appear
  until the buffer fills or the job ends — an explicit `flush` is the fix.

---

## Architecture

### 1. `scale_k_to_kcat` — one knob for relative vs absolute

**A single `scale_k_to_kcat::Union{Real,Nothing}` field on the problem structs decides both the
loss shape and the output rescale.** There is no separate `rate_scaling` enum — the two
sensible modes are exactly the two values this field can take:

| `scale_k_to_kcat` | loss (`loss!`)        | output (`fit_rate_equation`)        | meaning                          |
|-------------------|-----------------------|-------------------------------------|----------------------------------|
| a `Real` (def 1.0)| per-group **centered**| **rescale** SS k's so kcat = value  | relative y-axis (current default)|
| `nothing`         | **uncentered**        | **no** rescale; raw fitted k's      | absolute per-enzyme turnover     |

The other two (center+no-rescale, uncenter+rescale) are nonsense — center+no-rescale is the old
`kcat=nothing` "arbitrary k" footgun, and uncenter+rescale discards the absolute scale you have.
Collapsing to one `nothing`-or-value field makes them unrepresentable.

- `FittingProblem{M,D}` and `IdentifyRateEquationProblem{R,D}` each gain a field
  `scale_k_to_kcat::Union{Real,Nothing}`. Both constructors gain
  `scale_k_to_kcat::Union{Real,Nothing} = 1.0` (preserves today's default: centered +
  rescale-to-1.0). A `Real` value is validated `> 0` (kcat is a positive turnover).
- It sits exactly where `Keq` sits — a problem property set at construction — so
  `fit_rate_equation` and `identify_rate_equation` **drop the kwarg** and read `fp`/`prob`.
  `IdentifyRateEquationProblem.scale_k_to_kcat` threads into **every** `FittingProblem` the
  pipeline builds — `_process_batch` (`:346`), `_loocv` (`:534`), `_evaluate_loss` (`:500`) — so
  training, beam-search, and CV share one loss definition. **Caveat to document loudly:** this
  field controls the *loss* (centering), not just output scaling — `nothing` changes how the fit
  is computed, not merely how k's are reported.

**`loss!` — center iff `fp.scale_k_to_kcat !== nothing`, in pass 2 only.** Pass 1 (fill `buf`
with per-point log-residuals or the sign-mismatch sentinel `10.0`) is unchanged. Pass 2:

```julia
total_loss = 0.0
if fp.scale_k_to_kcat === nothing            # absolute: y is per-enzyme turnover
    @inbounds for i in 1:n_data
        total_loss += buf[i] * buf[i]        # uncentered
    end
else                                          # relative
    @inbounds for grp_idx in fp.group_point_indexes   # current per-group mean-centering
        ...
    end
end
```

`fp.scale_k_to_kcat === nothing` is read once into a local at the top of `loss!` (small-union
field, cheap). The sign-mismatch sentinel and the post-hoc `+ 100.0 * n_mismatch` penalty are
**unchanged** in both modes (the penalty's "centering would cancel a uniform sentinel" rationale
applies only to the centered branch, but a flat extra penalty on sign mismatches is harmless and
beneficial uncentered too). **`rate_equation` is not called differently and is not modified** —
the sub-100 ns / 0-alloc perf gate is unaffected (`loss!` is not itself perf-gated, and the
added work is one branch).

**`fit_rate_equation` — rescale driven by `fp.scale_k_to_kcat`; return retcode.**

- **Drops** its `kcat` kwarg entirely (reads `fp.scale_k_to_kcat`); takes only optimizer kwargs.
- `fp.scale_k_to_kcat isa Real` → `rescale_parameter_values(fp.mechanism, full;
  scale_k_to_kcat = fp.scale_k_to_kcat)`. `=== nothing` → **no** rescale, return raw fitted params.
- **Return `(params, loss, retcode)`** (was `(params, loss)`). See §2.
- Binding constants (K's) are unaffected by rescaling — kcat is homogeneous degree-1 in SS k's
  and independent of RE K's — so K's are always in real concentration units; only SS rate
  constants are normalized, and only when `scale_k_to_kcat isa Real`.

**`rescale_parameter_values` (public API) — kwarg renamed `kcat` → `scale_k_to_kcat::Real`** for
one-concept-one-name consistency (comment 2). It genuinely scales k's to a target kcat, so the
name fits. **Not renamed:** `_kcat_forward`, `_kcat_components`, `analytical_kcat_fn`, and any
`kcat` naming the *computed turnover quantity* — those are the value, not a scale target.

### 2. Solver outcome recording

- In `fit_rate_equation`'s multi-start loop, track the **best-objective restart's** `retcode`
  alongside `best_loss`/`best_x`. Store it as a short `Symbol` (e.g. `:Success`, `:MaxTime`)
  derived from `sol.retcode`. (Exact extraction from the SciMLBase `ReturnCode` enum — `Symbol`
  vs a name helper — is confirmed in the implementation plan; `Optimization` re-exports
  SciMLBase so no new dep.)
- `MaxTime` is treated as a genuine "did not converge" outcome, **not** success: with
  `maxiters = 10_000_000` the only binding budget is `maxtime`, so a `MaxTime` stop means
  CMA-ES's own convergence criteria did not fire. No special handling beyond recording — there
  is **no** warning on a non-`Success` best fit (it is recorded in the CSV/`cv_results` and the
  user filters).
- `BatchEntry` gains `retcode::Symbol`. The row `NamedTuple` gains a `retcode` (String for CSV)
  and an `error` (see §3) field. `_rows_to_dataframe` adds `retcode` and `error` columns; these
  propagate to `cv_results` (built from candidate rows).

### 3. Loud failures — no silent drops

**`_process_batch` — capture, don't swallow.** The per-mechanism `pmap` closure keeps a
`try/catch`, but the `catch` **captures** the exception (type + truncated message) instead of
discarding it. The function returns **two** lists:

```
_process_batch(mechs, prob; …) → (entries::Vector{BatchEntry}, failures::Vector{FitFailure})
```

- A mechanism that compiles, caps OK, and fits → a `BatchEntry` (with `retcode`, `error =
  missing` in its row).
- A mechanism that throws (StackOverflow/OOM at compile, an unsupported-kwarg `MethodError`,
  etc.) → a `FitFailure` record `(mech, error::String, n_params::Union{Int,Missing})` where
  `error` is `"<ExceptionType>: <truncated msg>"` and `n_params` is filled if compile got far
  enough, else `missing`.
- Rows written to CSV are **fitted rows + failure rows**. A failure row has `retcode = missing`,
  `loss = missing`, params `missing`, and `error` set. A fitted row has `error = missing`. This
  lets the user tell a *compile blowup* (`StackOverflowError` — expected for huge mechanisms)
  from a *config bug* (`MethodError` — fix the kwargs) directly in the CSV.
- **Row schema:** a failure row is the **same** `NamedTuple` schema as a fitted row, with
  `missing` in every field unavailable at the failure point: `loss`, `retcode`, `mechanism_type`,
  `rate_equation`, `eq_hash` all `missing`; `n_params` filled iff compile succeeded; `error`
  set; `fitted_param_names = ()` and `fitted_param_values = ()` (so it contributes no param
  columns). `_rows_to_dataframe` must tolerate `missing` in `n_params`/`loss`/`eq_hash`/
  `retcode` (columns widen to `Union{…,Missing}`); its `all_pnames` gather over
  `fitted_param_names` is unaffected by the empty tuples.

**Whole-base-tier failure re-raises.** A deterministic config error fails the *first* fit and
therefore the entire **base** tier uniformly. So:

- Base tier: if `isempty(base_entries) && !isempty(base_failures)` → **raise** an informative
  error embedding the first captured exception
  (`"Every base-tier fit failed (N mechanisms) — usually an optimizer/config problem. First
  failure: <ExceptionType>: <msg>"`). This is where the popsize-style bug surfaces loudly.
- Expansion tiers: a tier whose children all fail is **not** fatal (they are likely just
  too-large mechanisms; a config error would already have aborted the base tier). Record the
  failures in the iteration CSV + summary and continue with the rest of the frontier.

`_ingest!` consumes only `entries` (failures never enter the frontier/`cv_pool`/`best_loss`).

**`_loocv` — loud.** Remove the `try/catch` **and** the silent empty-return-on-non-finite:

- Any exception during a fold propagates.
- A non-finite (`Inf`/`NaN`) fold test loss **raises**, naming the candidate and the held-out
  group.
- The legitimate zero-loss floor at `eps(Float64)` stays (a single-row held-out group can give
  exactly 0); only `Inf`/`NaN` raises.

Since `_loocv` runs inside `_cv_model_selection`'s `pmap`, a raise aborts model selection — but
the search CSVs are already on disk, so the user re-runs only the cheap CV step. This **deletes
code** in `_cv_model_selection` / `_select_best_n_params`: the `isempty(v) ? Inf` `cv_score`
path, the `valid = filter(!isempty, …)` + "no finite LOOCV scores" error, and the
`all(!isfinite, cv_df.cv_score)` backstop all become dead (every `_loocv` now returns full
finite scores or raises) and are removed.

### 4. Progress output

A small helper, gated by a new `show_progress::Bool = true` kwarg on `identify_rate_equation`
(threaded into `_beam_search` and `_cv_model_selection`):

```julia
# Stage line → flushed stdout (visible in REPL AND redirected cluster logs) and appended to
# <save_dir>/progress.log (a durable artifact you can `tail -f`). Gated by show_progress.
function _progress(save_dir::AbstractString, show_progress::Bool, msg::AbstractString)
    show_progress || return nothing
    println(msg); flush(stdout)
    isdir(save_dir) || mkpath(save_dir)
    open(joinpath(save_dir, "progress.log"), "a") do io
        println(io, msg)
    end
    nothing
end
```

- `progress.log` is not a `.csv`, so it does not trip the "save_dir already contains CSV files"
  guard. `show_progress` gates both stdout and the file (file-only is a non-goal/YAGNI).
- **Master-level stages only** (emitted between `pmap` batches; no per-mechanism live bar — that
  needs worker→master plumbing, deferred):
  1. start — reaction summary, `scale_k_to_kcat` (with "relative"/"absolute turnover" gloss),
     `save_dir`;
  2. `"Enumerating initial mechanisms…"` → `"N initial mechanisms"`;
  3. base fit — `"Fitting N initial mechanisms…"` → batch summary line (below);
  4. each expansion iteration `i` — `"Iteration i (target n_params=c): expanding B → C
     children; fitting…"` → batch summary line;
  5. CV — `"Cross-validating K candidates (LOOCV)…"`;
  6. selection — `"Selected: <mechanism type>, n_params=c"`;
  7. end — overall summary (cumulative).
- **Batch summary line** reports the three buckets Denis asked for, no `MaxTime` singling-out:
  `"<F> fitted | Success <a>% | non-Success retcode <b>% | errored <c>% | best loss <x>"`,
  where the denominator is `F_fitted + n_failures` for that batch (and cumulative for the end
  summary). `Success%` compares the recorded `retcode` to `:Success`; `non-Success retcode%` is
  fits that ran with any other retcode; `errored%` is the captured `Failure`s.

### API surface (summary of changes)

- `FittingProblem(mechanism, table; Keq, scale_k_to_kcat=1.0)` — new field + kwarg
  (`Union{Real,Nothing}`).
- `IdentifyRateEquationProblem(reaction, table; Keq, scale_k_to_kcat=1.0)` — new field + kwarg.
- `fit_rate_equation(fp, optimizer; …)` — **drops** the `kcat` kwarg (reads `fp.scale_k_to_kcat`),
  **returns `(params, loss, retcode)`**.
- `rescale_parameter_values(m, params; scale_k_to_kcat=1.0)` — kwarg renamed from `kcat`.
- `identify_rate_equation(prob; show_progress=true, …)` — `show_progress` is the one new kwarg;
  `scale_k_to_kcat` is **not** a kwarg here (it lives on `prob`, like `Keq`).
- CSVs (`initial_mechanisms.csv`, `equation_search_iteration_N.csv`) and `cv_results` gain
  `retcode` and `error` columns; CSVs additionally contain failure rows. `<save_dir>/progress.log`
  is written when `show_progress`.

---

## Tests — TDD, one failing test (set) per new/changed unit, written first

- **`scale_k_to_kcat` validation** — a `Real` value passes; `0`/negative errors; `nothing`
  passes (both constructors).
- **`loss!` (`nothing` vs `Real`)** —
  - single-group crafted dataset: a `Real` `scale_k_to_kcat` gives loss == 0 (centering
    degenerates) while `nothing` gives loss > 0 equal to the hand-computed Σ(log|pred|−log|meas|)² / n;
  - two-group dataset with one group's rates scaled by a constant: a `Real` gives loss
    **invariant** to the scaling, `nothing` **changes** by the expected amount;
  - sign-mismatch sentinel/penalty behaves identically (a mismatch point contributes in both
    branches).
- **`fit_rate_equation` / `rescale_parameter_values`** — on a fixture mechanism:
  - `scale_k_to_kcat isa Real` → `_kcat_forward(result) ≈ scale_k_to_kcat` (e.g. 1.0, and a
    non-default value); `rescale_parameter_values(…; scale_k_to_kcat=…)` honors the renamed kwarg;
  - `scale_k_to_kcat === nothing` → result params equal the raw fitted params (no rescale);
  - `fit_rate_equation` returns a 3-tuple whose `retcode` is a `Symbol`.
- **`scale_k_to_kcat` threading** — the value set on `IdentifyRateEquationProblem` reaches the
  `FittingProblem`s built in `_process_batch`, `_loocv`, and `_evaluate_loss` (same loss in
  training, beam search, and CV).
- **`_process_batch`** — on a tiny dataset:
  - all-fit case returns `entries` with `retcode` set and `error = missing`, empty `failures`;
  - a mechanism forced to throw (spy optimizer / injected error) appears in `failures` with
    `error = "<Type>: …"` and is **absent** from `entries`; its CSV row has `error` set,
    `retcode`/`loss`/params `missing`;
  - over-`max_param_count` mechanism is never fit (existing behavior preserved);
  - whole-batch throw → `entries` empty, `failures` non-empty.
- **base-tier all-fail re-raise** — a `_beam_search` whose optimizer throws on every fit raises
  an error whose message contains the captured exception type/message (simulate the popsize
  case via an optimizer stub that rejects a kwarg).
- **`_loocv` loud** — a fold whose fit throws propagates the error; a fold yielding a non-finite
  test loss raises naming the candidate + held-out group; a legitimate zero test loss is floored
  to `eps` and does **not** raise.
- **`_progress`** — `show_progress=true` writes the expected line to `progress.log` and to
  stdout (capture via `redirect_stdout`); `show_progress=false` writes neither; the batch
  summary line reports the three percentage buckets with the correct denominator.
- **CSV/`cv_results` columns** — `retcode` and `error` columns are present; a failure row is
  representable (missing loss/params) alongside fitted rows.
- **end-to-end** — small dataset, both a `Real` and a `nothing` `scale_k_to_kcat` run to
  completion; `equation_search_iteration_N.csv` + `initial_mechanisms.csv` + `progress.log` exist.

**Unchanged / must stay green:** the `rate_equation` perf gate (alloc-free, <100 ns) and the
flat-string / Expr-shape regression tests in `test/test_rate_eq_derivation.jl` (derivation
untouched); the compile-budget gate (~750 bi-bi init); Aqua/JET; export count (18 — no export
changes). **Migration:** update callers that destructure `fit_rate_equation`'s old 2-tuple
return; rename every `kcat=…` kwarg meaning *scale-target* to `scale_k_to_kcat=…` — at
`fit_rate_equation` call sites (now via the `FittingProblem` field, not the dropped kwarg) and
`rescale_parameter_values` call sites. Audit existing `kcat=nothing` usages: under the new design
`nothing` *also* uncenters the loss (not just skips rescale), so any test asserting the old
centered-but-unscaled behavior must be updated. Do **not** rename `kcat` where it names the
computed turnover (`_kcat_forward`, `analytical_kcat_fn`, kcat-value docstrings).

---

## Doc sync

- `.claude/CLAUDE.md`:
  - `src/fitting.jl` description — `FittingProblem` carries `scale_k_to_kcat::Union{Real,Nothing}`;
    `loss!` centers only when it is a `Real`; `fit_rate_equation` reads it from `fp` (no kwarg)
    and returns `(params, loss, retcode)`. `rescale_parameter_values` kwarg `kcat` →
    `scale_k_to_kcat`.
  - `src/identify_rate_equation.jl` description / `identify_rate_equation` docstring —
    `scale_k_to_kcat` on `IdentifyRateEquationProblem`; `show_progress` kwarg;
    `retcode`/`error` CSV+`cv_results` columns and failure rows; `_loocv` is loud;
    `progress.log` artifact.

## Deferred / non-goals

- **Absolute mode with a known per-group `E_total` column** (the "bulk rate ÷ measured enzyme"
  variant). Chose per-enzyme turnover (`E_total = 1`, drop centering). Revisit if a dataset
  needs per-group enzyme amounts.
- **Per-mechanism live progress** (a fit counter / progress bar) — needs worker→master plumbing
  (`RemoteChannel` / `progress_pmap`). Master-level stages only for now.
- **File-only progress** (decoupling the `progress.log` write from the stdout gate). YAGNI.
- **"Warn loudly if >X% of fits hit a non-Success retcode" threshold** — recorded per fit +
  summarized; no abort, no threshold knob. YAGNI.
- **Preferring a converged restart over a lower-objective `MaxTime` restart** in
  `fit_rate_equation` — record the best-objective restart's retcode; do not change the selection
  criterion.
