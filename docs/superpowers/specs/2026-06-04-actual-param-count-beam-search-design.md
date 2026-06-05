# Actual-Param-Count Beam Search — Design

**Date:** 2026-06-04
**Status:** design approved, pending spec review → implementation plan

---

## Goal

Replace the parameter-*estimate*-driven beam search with one keyed on the **actual**
fitted-parameter count (`length(fitted_params(compile_mechanism(m)))`), learned from the
compile pass. Three coupled changes:

1. **Param-count axis.** Delete `_n_fit_params_estimate` entirely. Fit all `init_mechanisms`
   up front; run a single ascending-by-actual-count expansion loop seeded by that base set.
2. **Bounded memory.** Make `save_dir` mandatory; stop accumulating the full result set in
   memory; keep only what's needed downstream (a bounded top-N-per-param-count CV pool, the
   live frontier, a tiny seen-equation set). Purge per-iteration rows after each CSV write.
3. **Dedup simplification.** Delete the ~200-line canonical-rate-equation-hash machinery and
   the cross-mechanism fit-projection layer. Dedup on a **comment-stripped rate-equation
   string key** (cheaper *and* more correct — see probe results).

## Why (probe evidence, LDH `NADH+Pyruvate ⇌ Lactate+NAD`, `oligomeric_state:4`)

- `init_mechanisms` returns 69 mechanisms (post-#41), all distinct equations; `fitted_params`
  count ∈ {5,6}. They are **sibling base topologies**, not a complexity ladder — bucketing
  them by parameter count is the wrong axis (compulsory-ordered vs rapid-equilibrium-random
  differ in count by a thermodynamic identification, not by elaboration).
- `_n_fit_params_estimate` collapsed after #40's equivalence-grouping (a 5-param mechanism
  estimates as 1): the estimate counts kinetic *groups* while the thermodynamic term still
  uses the full *step* count, so the most-merged mechanisms estimate far below their real
  fitted count. That mislabels the per-level CSVs and applies `max_param_count` in
  estimate-space (admitting far more complex mechanisms than intended) — the observed
  `identify_ldh.jl` regression.
- **Dedup probe (632 mechanisms across LDH non-allosteric + allosteric expansion):**
  - Structural dedup (`dedup!`) removes the bulk of duplicates (76 / 258 / 12 per round) but
    cannot catch structurally-distinct-but-equation-identical mechanisms.
  - In LDH expand round 1: 512 structural-unique children → canonical hash finds 504 distinct
    (catches 8); full rate-eq string finds 508 (catches 4); **comment-stripped string finds
    491 (catches 21).** The 8 canonical-dups are cosmetic: 4 byte-identical strings, 4 differ
    only in the Wegscheider provenance comment (same `v=`, same Haldane, same free params).
  - The canonical hash is **over-conservative**: it keeps 13 pairs distinct that have
    byte-identical `v=` + Haldane + fitted-param set (it hashes eliminated-param `dep_exprs`
    that don't affect the rate). The comment-stripped key merges those correctly.
  - **Comment-stripped key: zero false merges across all 632 mechanisms**, strictly ⊇ the
    canonical hash's catches. Allosteric sample (120) had no duplicates; all three methods
    agreed (safety confirmed; allosteric *catch* behavior untested but can only under-merge,
    never wrong-merge).

---

## Architecture

All changes in `src/identify_rate_equation.jl` and `src/mechanism_enumeration.jl`.

### 1. `expand_mechanisms` → flat vector (`mechanism_enumeration.jl`)

- `expand_mechanisms(mechs, rxn)` returns `Vector{Union{Mechanism, AllostericMechanism}}`.
- `_add_expansions_mech!` pushes into a flat vector; `_push_mech!` becomes a plain `push!`
  (delete the `pc = _n_fit_params_estimate(m)` line and the `Dict` bucketing).
- Bucketing by parameter count is no longer enumeration's job.

### 2. Delete `_n_fit_params_estimate`

- Remove both overloads (`mechanism_enumeration.jl:1084`, `:1102`) and all call sites.
- Remove the `n_subs + n_prods + 1` floor logic wherever it floored the estimate.
- Update tests that assert estimate deltas (see Tests).

### 3. Comment-stripped equation dedup key (`identify_rate_equation.jl`)

Replace the canonical-hash machinery with:

```julia
# Equation identity key: the rendered string with provenance comments removed
# (the choice of which dependent K was eliminated is cosmetic — it is already
# substituted into v). Two mechanisms with the same key compute the identical
# rate function. Used as a CSV column tag and the CV-candidate selection key.
function _rate_eq_dedup_key(eq_text::AbstractString)
    kept = Iterators.filter(split(eq_text, '\n')) do ln
        l = strip(ln)
        !startswith(l, "#") && !occursin("(substituted into v)", l)
    end
    hash(join(kept, '\n'))
end
```

`eq_hash` is **not** a collapse key. Every structurally-distinct mechanism is fit and its
row kept (with its own fitted params — no projection/copying). `eq_hash` serves two purposes
only: a CSV column for downstream grouping in analysis, and the LOOCV candidate-selection key
(pick N *distinct* equations per param count, not N best mechanisms). The search loop never
dedups by `eq_hash`; the only mechanism-collapsing dedup is the structural `_dedup_flat`.

**Delete entirely:** `_canonical_rate_eq_hash_data`, `_canonical_rate_eq_hash`,
`_canonicalize_for_hash`, `_build_name_map`, `_dep_exprs_canonical`, `_synth_dep_a_names`,
`_expr_canonical_via_name_map` (and any helper used only by them), `_project_cached_params`,
the whole `_CachedFitResult` struct (replaced by `BatchEntry`), and the inherited-row path.
(Grep `_canonical_rate_eq_hash`, `_project_cached_params`, `name_map`, `canon_to_rep`,
`_CachedFitResult` to confirm no remaining src callers.)

The `eq_hash` CSV column stays, now holding `_rate_eq_dedup_key(eq_text)` rendered as a hex
string. The `fit_inherited_from_estimate` column is **removed** (no inherited rows, no
cross-mechanism projection).

### 4. `_process_batch` — one pmap pass, compile+fit fused per worker

**Cluster rationale (comment 2):** compile and fit happen on the **same worker** for each
mechanism, in a single `pmap`. Splitting compile into its own pass would re-pay the
`@generated` JIT on whichever worker later draws the fit (Julia can't ship a compiled
specialization between processes) — every fitted mechanism would compile twice. Fusing them
means each surviving mechanism is compiled exactly once, on the worker that fits it.

```
_process_batch(mechs, prob; pmap_function, optimizer, max_param_count, kwargs...)
    → Vector{BatchEntry}     # one per fitted mechanism
```

`mechs` is already structurally pre-deduped on the master (`_dedup_flat`). There is **no**
`seen_keys` and **no** equation collapse — every structurally-distinct mechanism that fits
within the cap gets its own `BatchEntry`.

Each worker, for one mechanism, runs the full chain locally:
1. `em = compile_mechanism(m)`; `n = length(fitted_params(em))`. (compile failure → skip)
2. **Cap filter, before fitting:** `n > max_param_count` → skip.
3. **Fit** `FittingProblem(em, prob.data; Keq)`; compute `eq_hash =
   _rate_eq_dedup_key(rate_equation_string(em))`.
4. Return `(mech=m, n_params=n, eq_hash, loss, params, eq_text)`.

On the **master**, after the `pmap` returns: drop skips, build one row + `BatchEntry` per
fitted mechanism (each carrying its *own* fitted params and `eq_hash`). No collapse, no
projection.

This is why fusing compile+fit is clean here: with no cross-worker dedup to coordinate, each
worker is a self-contained compile→cap→fit unit, and the only shared state is the immutable
`prob`. The race condition from earlier designs is gone because there is no shared write at
all. Equation-level redundancy (the ~1.6% same-equation-different-structure mechanisms) is
intentionally kept — both rows are written (tagged by `eq_hash`) and both are eligible to
expand, since their structural expansions differ.

`BatchEntry` is a small new struct
`(mech, n_params::Int, loss::Float64, eq_hash::UInt64, row::NamedTuple)`.

`_dedup_flat(ms)::Vector` is a new helper applying the same structural canonicalization +
`unique!` as `dedup!` to a flat vector (factor the per-bucket body out of `dedup!(::Dict)`
and call it from both). It is the **only** mechanism-collapsing dedup — cheap, pre-compile,
dropping exact structural twins before the expensive compile pass.

### 5. `_beam_search` rewrite — ascending actual-count loop

State:
- `frontier::Dict{Int, Vector{BatchEntry}}` — unexpanded fit entries, keyed by **actual**
  `n_params`. The live work queue. Holds *all* structurally-distinct mechanisms (no eq-dedup).
- `cv_pool::Dict{Int, Vector{BatchEntry}}` — top `n_cv_candidates` **distinct equations**
  (by `eq_hash`, lowest loss each) per param count, accumulated across iterations. The *only*
  thing retained for final CV.
- `best_loss_by_count::Dict{Int, Float64}` — per-count best-ever loss, the beam-cutoff
  reference (a count's best may be an already-expanded entry, so it's tracked separately).
- `target::Int` — a monotonically-increasing sweep pointer (the highest count tier admitted
  so far).

No `seen_keys` / persistent fit cache. **Expansion moves are not param-count-monotonic** — a
move can produce a child with the *same* (Δ=0, ≈16% of LDH expansions) or even fewer params
than its parent, because a Haldane/Wegscheider constraint absorbs the added structural
parameter. Termination instead rests on **irreversibility**: every move (`_expand_re_to_ss`,
`_expand_split_kinetic_group`, `_expand_add_dead_end_regulator`, `_expand_to_allosteric`,
`_expand_add_allosteric_regulator`, `_expand_change_allo_state`) is add-only with no inverse,
so structure strictly elaborates along every path, the reachable structure space under
`max_param_count` is finite, and the search cannot cycle.

Flow:
```
base = _dedup_flat(collect(init_mechanisms(prob.reaction)))   # structural pre-dedup
base_entries = _process_batch(base, prob; max_param_count, …)
_save_initial_csv(save_dir, [e.row for e in base_entries])    # mandatory
_ingest!(frontier, cv_pool, best_loss_by_count, base_entries; n_cv_candidates)

iteration = 0
target = minimum(keys(frontier))                  # lowest count present
while !isempty(frontier):
    # Sweep this tier PLUS every same-or-lower-count straggler into one batch.
    group = BatchEntry[]
    for c in collect(keys(frontier))
        c <= target && append!(group, pop!(frontier, c))
    end
    to_expand = BatchEntry[]
    for c in unique(e.n_params for e in group)    # beam-select PER count
        ec  = [e for e in group if e.n_params == c]
        sel = _select_beam([e.loss for e in ec];
                 best_override = best_loss_by_count[c], thresholds…, min_beam_width)
        append!(to_expand, ec[sel])               # non-selected are pruned (dropped)
    end

    if !isempty(to_expand)
        children = _dedup_flat(
            expand_mechanisms([e.mech for e in to_expand], prob.reaction))
        child_entries = _process_batch(children, prob; max_param_count, …)
        iteration += 1
        _save_iteration_csv(save_dir, [e.row for e in child_entries], iteration)
        _ingest!(frontier, cv_pool, best_loss_by_count, child_entries; n_cv_candidates)
    end

    isempty(frontier) && break
    target = max(target + 1, minimum(keys(frontier)))   # advance; jump over gaps
```

- **Advancing-sweep / straggler handling:** `target` only increases, so a Δ≤0 child (same or
  lower count than the tier that produced it) lands at a count `≤ target` and is swept into
  the **next** iteration's batch — pooled with the next higher tier, never expanded in a lonely
  small pass. This is robust to Δ<0 for free. Δ=0 is common (~16% on LDH), so this matters in
  practice, not just defensively. A child's anomalous low `n_params` is visible directly in
  `equation_search_iteration_N.csv` (no separate log).
- **`_select_beam` change:** add a `best_override` kwarg so the cutoff is computed against the
  per-count best-ever loss (`best_loss_by_count[c]`), not just the min of the current group —
  so a straggler swept in later is judged against its tier's real best.
- **Termination:** irreversible structural elaboration + the `max_param_count` cap (over-cap
  children are dropped at the fit filter) ⇒ a finite reachable structure space ⇒ the loop
  ends. `target` strictly increases each iteration, guaranteeing forward progress. (No
  seen-set: YAGNI, justified by the irreversibility of every current move.)

### 6. `_ingest!` — bounded memory (comment 2)

For each new `BatchEntry`:
- push into `frontier[n_params]` (unexpanded work — **all** structurally-distinct mechanisms,
  no eq-dedup, so both members of a same-equation pair stay eligible to expand);
- update `best_loss_by_count[n_params]` to the running min;
- offer it to `cv_pool[n_params]`, which keeps the top `n_cv_candidates` **distinct
  equations** by loss: if the entry's `eq_hash` is already present, keep the lower-loss one;
  else if fewer than `n_cv_candidates` distinct equations are held, add it; else if it beats
  the worst-loss held entry, replace that entry. **Mechanisms that fall out are dropped.**

The `cv_pool`'s distinct-`eq_hash` rule is exactly the LOOCV requirement ("N different
equations, not N best mechanisms") computed incrementally, so `_cv_model_selection` consumes
the pool directly. Steady-state memory ≈ `(#distinct param counts) × n_cv_candidates`
mechanisms + the live frontier (bounded by beam width × expansion fan-out) — independent of
iteration count. (The frontier carries the ~1.6% same-equation redundancy; it does not
compound unboundedly because each generation is structurally pre-deduped and the beam is
size-capped.)

### 7. CSV outputs + `_cv_model_selection`

- `save_dir::String` always points at a real directory — it gets a **computed default** so a
  run always writes results. Default = `_default_save_dir()`:

  ```julia
  # First non-existent "<date>_results[_N]" dir in the cwd, e.g.
  # 2026_06_04_results, then 2026_06_04_results_2, _3, …
  function _default_save_dir()
      base = string(Dates.format(Dates.today(), "yyyy_mm_dd"), "_results")
      isdir(base) || return base
      n = 2
      while isdir(string(base, "_", n)); n += 1; end
      string(base, "_", n)
  end
  ```

  Evaluated per-call as the default argument, so each run picks the next free suffix. An
  explicitly-passed `save_dir` is used verbatim; the "directory already contains CSVs" guard
  (`identify_rate_equation.jl:190-198`) still fires for an explicit dir that holds prior CSVs
  (the auto-suffixed default never collides). Requires `using Dates` (already a stdlib dep).
- `_save_level_csv(…, pc)` → split into `_save_initial_csv(dir, rows)` writing
  `initial_mechanisms.csv` and `_save_iteration_csv(dir, rows, iteration)` writing
  `equation_search_iteration_$(iteration).csv` (`iteration` is a sequential counter, **not** a
  param count — real count is the `n_params` column).
- `_cv_model_selection` receives the already-bounded `cv_pool`. Its distinct-`eq_hash`
  candidate selection (`:1061-1069`) is **kept** — LOOCV must evaluate N *different equations*
  per param count, not N best mechanisms — but it now operates on the small pool (which was
  built by the same distinct-equation rule during ingest) instead of the full result set, and
  reads the comment-stripped `eq_hash`.
- `max_param_count` is documented as a cap on **actual fitted** parameters; base topologies
  are never excluded by it (they are the seed), only expansion depth is capped. Default base
  set is **uncapped** (YAGNI; add an explicit `max_base_topologies` that `log()`s drops only
  if a real bottleneck appears).

---

## Tests — TDD, one unit test (set) per new/changed function

Every new or modified function gets its own failing test written **first**, then just enough
code to pass (CLAUDE.md TDD). The implementation plan sequences each function as
test→implement. Per-function unit tests:

- **`_rate_eq_dedup_key`** — (a) two rate-equation strings differing only in `# …` headers or
  `(substituted into v)` lines hash equal; (b) strings differing in a Haldane definition or a
  `v=` term hash *unequal*; (c) on the LDH expand-round-1 fixture the 8 known equation-dups
  (4 byte-identical + 4 comment-only) share a key and there are **no false merges** (every
  key collision has byte-identical params-decl + Haldane + `v=`).
- **`_dedup_flat`** — a vector with a structural twin collapses to one; an empty vector and an
  all-distinct vector are unchanged; result matches `dedup!` on the single-bucket dict form.
- **`_default_save_dir`** — returns `<today>_results` when absent; returns `…_results_2`,
  `…_3` when prior dirs exist (use a temp cwd with pre-made dirs).
- **`_save_initial_csv` / `_save_iteration_csv`** — write the expected filenames; the
  iteration writer name encodes the sequential counter, not a param count; round-trips the
  rows DataFrame.
- **`_select_beam` (`best_override`)** — with an override lower than the group min, fewer
  entries pass the cutoff than without it; `min_beam_width` still honored; no override
  reproduces current behavior (regression-pins existing callers).
- **`_process_batch`** — on a tiny dataset: (a) compiles+fits and returns one `BatchEntry` per
  fitted mechanism, each with its own `eq_hash` and fitted params; (b) an
  over-`max_param_count` mechanism is **never fit** (assert via a spy `pmap`/optimizer or a
  fit-count probe); (c) two structurally-distinct same-equation mechanisms in one batch yield
  **two** `BatchEntry`s sharing one `eq_hash` (no collapse).
- **`_ingest!`** — pushes every entry to `frontier[n]` (both members of a same-equation pair
  stay); updates `best_loss_by_count` to the min; `cv_pool[n]` keeps the top
  `n_cv_candidates` **distinct `eq_hash`** by loss (a same-equation entry replaces only its
  own hash's slot, never consuming a second slot) and drops overflow mechanisms.
- **`_cv_model_selection`** — given a pool where one param count holds several mechanisms with
  a repeated `eq_hash`, the selected CV candidates are N *distinct* equations (no equation
  CV-evaluated twice); fed by the bounded pool, not the full result set.
- **`expand_mechanisms`** — returns a flat `Vector{Union{Mechanism,AllostericMechanism}}`
  (contract change from `Dict`).
- **`_beam_search` integration** — end-to-end on a small LDH-like dataset with an explicit
  `save_dir` and a low `max_param_count`. Asserts: the run **terminates** (this is the real
  test of the irreversibility/advancing-sweep termination argument — LDH expansion has ~16%
  Δ=0 children, so a real run exercises same-count stragglers); output dir contains exactly
  `initial_mechanisms.csv` + sequential `equation_search_iteration_N.csv` (no gaps);
  `initial_mechanisms.csv` row count == `length(_dedup_flat(init_mechanisms(rxn)))`; every CSV
  row carries an `eq_hash` column; every `n_params` ≤ `max_param_count`; the returned CV pool
  size ≤ `(#param counts) × n_cv_candidates` (memory bound — confirms rows are not all
  accumulated). Because `target` only advances, no `equation_search_iteration_N.csv` is a
  lonely single-mechanism file purely from a Δ=0 child.

**Delete:**
- `_n_fit_params_estimate` delta tests (`test/test_mechanism_enumeration.jl` ~1466/1521/1692/1835)
  and the floor-applying sites (~968/1064/1176/1252).
- Tests consuming `expand_mechanisms` as a `Dict` / asserting bucket keys.
- Canonical-hash tests and `_project_cached_params` tests.
- Tests asserting `params_estimate_*` filenames.

**Unchanged / must stay green:** `rate_equation` perf gate (alloc-free, <100 ns) and the
flat-string / Expr-shape regression tests in `test/test_rate_eq_derivation.jl` — the
derivation is **not** touched. The compile-budget gate (~750 bi-bi init) and export count
(18) hold (no new `Sig` types; no export changes).

---

## Doc sync

- `.claude/CLAUDE.md`: `expand_mechanisms → Vector{…}` (not `Dict`); note bucketing by actual
  `n_params` lives in `identify_rate_equation.jl`; remove the "parameter naming chokepoint /
  canonical hash" references that describe the deleted machinery; update the `identify_rate_equation`
  docstring (CSV names, `max_param_count` = actual fitted cap, `save_dir` required).

## Deferred / non-goals

- Positive allosteric equation-dedup catch-behavior (probe found no allosteric dups to
  exercise; risk is under-merge only). Revisit if a real allosteric run shows duplicate rows.
- `max_base_topologies` guard for very large reactions (uncapped until a bottleneck appears).
