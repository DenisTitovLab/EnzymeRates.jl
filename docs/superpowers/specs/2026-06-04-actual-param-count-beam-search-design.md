# Actual-Param-Count Beam Search — Design

**Date:** 2026-06-04
**Status:** design approved, pending spec review → implementation plan

Supersedes the recommendation in `.claude/2026-06-03-base-tier-fit-all-and-iteration-naming.md`
(that plan kept `_n_fit_params_estimate` as a pre-filter and kept the canonical-hash
machinery; this design removes both, based on the empirical probe results recorded below).

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
  estimates as 1), which is the root cause of the `identify_ldh.jl` regression
  (`.claude/2026-06-03-identify-ldh-regression-findings.md`).
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
# Dedup key for a rate equation: the rendered string with provenance comments
# removed (the choice of which dependent K was eliminated is cosmetic — it is
# already substituted into v). Two mechanisms with the same key compute the
# identical rate function. Hashed for compact storage / CSV column.
function _rate_eq_dedup_key(eq_text::AbstractString)
    kept = Iterators.filter(split(eq_text, '\n')) do ln
        l = strip(ln)
        !startswith(l, "#") && !occursin("(substituted into v)", l)
    end
    hash(join(kept, '\n'))
end
```

**Delete entirely:** `_canonical_rate_eq_hash_data`, `_canonical_rate_eq_hash`,
`_canonicalize_for_hash`, `_build_name_map`, `_dep_exprs_canonical`, `_synth_dep_a_names`,
`_expr_canonical_via_name_map` (and any helper used only by them), `_project_cached_params`,
the whole `_CachedFitResult` struct (replaced by `BatchEntry` + `seen_keys`), and the
inherited-row path. (Grep `_canonical_rate_eq_hash`, `_project_cached_params`, `name_map`,
`canon_to_rep`, `_CachedFitResult` to confirm no remaining src callers.)

The `eq_hash` CSV column stays, now holding `_rate_eq_dedup_key(eq_text)` rendered as a hex
string. The `fit_inherited_from_estimate` column is **removed** (no inherited rows).

### 4. `_process_batch` (extracted, reusable)

Extract today's inlined compile→dedup→fit→row body (`identify_rate_equation.jl:668-762`) into:

```
_process_batch(mechs, prob, seen_keys; pmap_function, optimizer, max_param_count, kwargs...)
    → Vector{BatchEntry}     # one per NEW distinct equation
```

Steps:
1. **Parallel compile** each mechanism → `(mech, n_params, eq_text, eq_key)` (drop compile
   failures). `n_params = length(fitted_params(m))`; `eq_key = _rate_eq_dedup_key(eq_text)`.
2. **Cap filter (right before fit):** drop entries with `n_params > max_param_count`.
3. **Dedup:** drop entries whose `eq_key ∈ seen_keys` (already fit in an earlier batch) and
   collapse same-`eq_key` entries within this batch to one representative. Add survivors'
   keys to `seen_keys`.
4. **Parallel fit** one representative per surviving `eq_key`.
5. Build one row + `BatchEntry(mech, n_params, loss, eq_key, row)` per representative.

`seen_keys::Set{UInt64}` is the only persistent dedup state — tiny (one hash per distinct
equation), replaces the old params-carrying `fit_cache`. No cross-mechanism param projection.

`BatchEntry` is a small new struct `(mech, n_params::Int, loss::Float64, eq_key::UInt64, row::NamedTuple)`.

`_dedup_flat(ms)::Vector` is a new helper applying the same structural canonicalization +
`unique!` as `dedup!` to a flat vector (factor the per-bucket body out of `dedup!(::Dict)`
and call it from both). It is the cheap pre-compile dedup that drops structural twins before
the expensive compile pass; `_rate_eq_dedup_key` is the post-compile equation dedup.

### 5. `_beam_search` rewrite — ascending actual-count loop

State:
- `seen_keys::Set{UInt64}` — equations already fit/recorded.
- `frontier::Dict{Int, Vector{BatchEntry}}` — unexpanded fit entries, keyed by **actual**
  `n_params`. The live work queue.
- `cv_pool::Dict{Int, Vector{BatchEntry}}` — top `n_cv_candidates` by loss per param count,
  accumulated across iterations. The *only* thing retained for final CV.
- `best_loss_by_count::Dict{Int, Float64}` — beam-cutoff reference per count.
- `high_water::Int` — highest count whose tier has been processed.

Flow:
```
base = _dedup_flat(collect(init_mechanisms(prob.reaction)))   # structural pre-dedup
base_entries = _process_batch(base, prob, seen_keys; max_param_count, …)
_save_initial_csv(save_dir, [e.row for e in base_entries])    # mandatory
_ingest!(frontier, cv_pool, best_loss_by_count, base_entries; n_cv_candidates)

iteration = 0
while any nonempty vector in frontier:
    # sweep the next tier PLUS any stragglers below high_water into one batch
    counts   = sort(collect(keys(frontier) with nonempty unexpanded))
    k        = first(counts)                       # min unexpanded count
    tiers    = [c for c in counts if c == k || c < high_water]   # k + stragglers
    to_expand = BatchEntry[]
    for c in tiers
        grp = pop!(frontier, c)
        sel = _select_beam([e.loss for e in grp];
                 best_override = best_loss_by_count[c], thresholds…, min_beam_width)
        append!(to_expand, grp[sel])              # non-selected are pruned (dropped)
    end
    high_water = max(high_water, maximum(tiers))

    children = _dedup_flat(expand_mechanisms([e.mech for e in to_expand], prob.reaction))
    child_entries = _process_batch(children, prob, seen_keys; max_param_count, …)
    iteration += 1
    _save_iteration_csv(save_dir, [e.row for e in child_entries], iteration)
    _ingest!(frontier, cv_pool, best_loss_by_count, child_entries; n_cv_candidates)
    # per-iteration rows are now only in the CSV + cv_pool; nothing else retained

return cv_pool flattened → (candidate_mechs, candidate_df)
```

- **Straggler handling (comment 1):** a child landing at a count `< high_water` is swept into
  the *next* iteration's combined batch via the `c < high_water` clause — it never triggers a
  lonely 1-mechanism fit pass. With today's strictly-additive expansion moves this branch
  never fires; a straggler, if it ever appears, is visible in the CSV as a row whose
  `n_params` is below the iteration's stage (no separate log needed).
- **`_select_beam` change:** add a `best_override` kwarg so the cutoff is computed against the
  best-ever loss at that count (`best_loss_by_count[c]`), not just the min of the current
  group — so a late straggler is judged against its tier's real best.
- **Termination:** structural + equation dedup (`seen_keys`) guarantee no equation is
  re-expanded; finitely many distinct mechanisms ≤ cap ⇒ the loop terminates.

### 6. `_ingest!` — bounded memory (comment 2)

For each new `BatchEntry`:
- push into `frontier[n_params]` (unexpanded work);
- update `best_loss_by_count[n_params]`;
- insert into `cv_pool[n_params]`, then truncate that vector to the `n_cv_candidates` lowest
  losses — **dropping the `Mechanism` objects that fall out.** (Equation-dedup already
  happened via `seen_keys`, so the pool holds distinct equations.)

Steady-state memory ≈ `(#distinct param counts) × n_cv_candidates` mechanisms + the live
frontier (bounded by expansion fan-out × beam width) + `seen_keys` (hashes only) —
independent of iteration count.

### 7. CSV outputs + `_cv_model_selection`

- `save_dir::String` is **required** (drop the `Union{Nothing,String}`; keep the
  "directory already contains CSVs" guard at `identify_rate_equation.jl:190-198`).
- `_save_level_csv(…, pc)` → split into `_save_initial_csv(dir, rows)` writing
  `initial_mechanisms.csv` and `_save_iteration_csv(dir, rows, iteration)` writing
  `equation_search_iteration_$(iteration).csv` (`iteration` is a sequential counter, **not** a
  param count — real count is the `n_params` column).
- `_cv_model_selection` now receives the already-bounded `cv_pool` (its internal
  "top-N per `(n_params, eq_hash)`" selection at `:1061-1069` becomes a no-op / is removed —
  the pool is already that set, keyed on the comment-stripped `eq_hash`).
- `max_param_count` is documented as a cap on **actual fitted** parameters; base topologies
  are never excluded by it (they are the seed), only expansion depth is capped. Default base
  set is **uncapped** (YAGNI; add an explicit `max_base_topologies` that `log()`s drops only
  if a real bottleneck appears).

---

## Tests

**Delete:**
- `_n_fit_params_estimate` delta tests (`test/test_mechanism_enumeration.jl` ~1466/1521/1692/1835)
  and the floor-applying sites (~968/1064/1176/1252).
- Tests consuming `expand_mechanisms` as a `Dict` / asserting bucket keys.
- Canonical-hash tests and `_project_cached_params` tests.
- Tests asserting `params_estimate_*` filenames.

**Add:**
- `expand_mechanisms` returns a flat `Vector`.
- `_rate_eq_dedup_key`: on the LDH expand-round-1 set, collapses the 8 known equation-dups
  (4 byte-identical + 4 comment-only) and produces **no false merges** (every collision has a
  byte-identical params-decl + Haldane + `v=`). Small, fixture-based.
- End-to-end on a tiny dataset with `save_dir`: output dir contains exactly
  `initial_mechanisms.csv` + sequential `equation_search_iteration_N.csv` (no gaps);
  `initial_mechanisms.csv` row count == `length(_dedup_flat(init_mechanisms(rxn)))`; every
  `n_params` ≤ `max_param_count`; no over-cap equation is ever fit.
- Memory bound: after a multi-iteration run, the in-memory CV pool size ≤
  `(#param counts) × n_cv_candidates` (assert we did not accumulate all rows).

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
