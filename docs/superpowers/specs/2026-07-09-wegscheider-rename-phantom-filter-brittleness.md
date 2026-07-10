# Wegscheider rename + phantom filter — derivation-machine brittleness (tech debt)

Date: 2026-07-09
Status: NOTE — not urgent. The current code is correct and the full suite is green
(PR #62). This records why the construct is brittle and the cleanest way to remove it.

## What the construct is

Two coupled pieces in the thermodynamic-constraint derivation:

1. **Wegscheider rename** (`_build_wegscheider_rename_map` for `Mechanism`,
   `_state_wegscheider_rename_map` for `AllostericMechanism`). It detects a
   *single-symbol* binding-K Wegscheider tie (`K_b = K_a`, both binding K's) and
   folds `K_b → K_a` by renaming. Its real job is **canonical form**: it feeds
   `_merge_tied_kinetic_groups`, so split and merged encodings of a rate-equivalent
   graph collapse to one partition (the load-bearing Canonical Step Form used for
   structural dedup).

2. **Phantom filter** in `_dependent_param_exprs`. Because the folded symbol `K_b`
   is still a kinetic-group representative, `_raw_param_symbols` / `_state_all_params`
   still emit it as a matrix **column**. Its cycle-incidence contributions are mapped
   onto `K_a`'s column (via `step_name = get(rename, …)`), so `K_b` becomes a
   **zero-column**: never pivoted, it lands in `indep` as a fittable dummy that
   appears nowhere in the rate equation. The filter drops it
   (`indep = filter(p -> get(rename,p,p) == p, …)`).

## Why it's brittle

- **The column set and the rename disagree.** The rename asserts `K_b` *is* `K_a`,
  yet `K_b` is still emitted as a column. The phantom is the residue of that
  disagreement, patched after the fact by the filter.
- **The phantom is a silent over-parameterization.** Nothing in the natural test
  surface catches it — a dummy fitted param breaks neither callability nor detailed
  balance; it only inflates the parameter count (hurting the parsimony-based model
  selection that is the package's whole point) and destabilizes finite-restart fits.
- **It bit us in PR #62.** The allosteric path has *per-state* renames and shared
  `:EqualAI` symbols. The naive non-allosteric one-liner over-dropped a symbol folded
  in one state's rename but still referenced by the other state's polynomial
  (`UndefVarError: K_Lactate_ENADH`); the shipped filter needs a reference guard
  (`get(rename,p,p) == p || p in refs`) that only exists because of the underlying
  fold-then-filter design.
- **The rename serves two conflated roles** — construction-time canonical dedup and
  solve-time folding — through one mechanism, so a change for one purpose perturbs the
  other.

## Cleanest fix (recommended)

**Apply the rename when building the column set, not after solving.** If
`_state_all_params` / `_raw_param_symbols` exclude renamed-away symbols
(`get(rename, s, s) != s`), then `K_b` is never a column, never a zero-column, never a
phantom — and the phantom filter (with its allosteric reference guard) can be deleted
outright. This is **output-equivalent**: the rate equation already uses `K_a` (the
rename is applied to the polynomials), so `fitted_params` and the rendered equation are
unchanged — no golden rebaseline. It must be applied **per state** for the allosteric
path (exclude from `cols_A` those the A-rename folds, from `cols_I` those the I-rename
folds) so a shared `:EqualAI` symbol folded in only one state survives via the other
state's column — exactly the case that broke the naive filter.

Open verification before implementing: confirm nothing else relies on the folded
symbol being present in the column set (the constraint matrix, the priority vector),
and that `_merge_tied_kinetic_groups` still sees what it needs (it consumes the rename
map directly, not the column set, so it should be unaffected).

## Alternative (more invasive)

PR #62's sentinel-free pivot means the solver can now pivot a single-symbol binding
tie **directly** (emitting an explicit `K_b = K_a` dependent), so the rename's
*solve-time* role is redundant — the fold was originally a workaround for the `-1`
sentinel dropping that row. One could stop passing the rename to the solve and keep it
only for construction-time canonical dedup. This removes the phantom too, but it
**changes output** (folded single name → explicit Wegscheider line, `K_b` used in the
equation), so it rebaselines the non-allosteric textbook goldens and is worth doing
only alongside a deliberate naming pass. Prefer the column-set fix above.
