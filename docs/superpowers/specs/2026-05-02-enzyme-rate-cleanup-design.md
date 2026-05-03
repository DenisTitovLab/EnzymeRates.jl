# EnzymeRates Cleanup — Validation, Display, Pipeline & Dedup Fixes

**Date:** 2026-05-02
**Author:** Denis + Claude (brainstorming)

---

## Process Requirement: TDD for Every Fix

**Every fix in this spec MUST follow strict TDD per `CLAUDE.md`:**

1. Write a failing regression test that captures the bug or the new
   requirement.
2. Run the test suite to confirm the test fails for the expected
   reason (not a typo, missing import, etc.).
3. Write the **smallest** code change that turns the test green.
4. Run the test suite again to confirm pass + no regression in any
   other test.
5. Refactor only after green.

**No implementation lands without a regression test asserting the bug
cannot recur.** When a fix is split across multiple files, each file's
behavior change has its own targeted test before the implementation
changes that file.

Where existing tests cover the affected surface (e.g.,
`test_types.jl`, `test_mechanism_enumeration.jl`,
`test_identify_rate_equation.jl`, `test_dsl.jl`), the regression tests
go into the same file. New test files are created only when a fix
introduces a wholly new surface.

The plan that follows MUST list, for each step, the failing test that
gates it. Implementation steps without a corresponding test step are
forbidden.

---

## Goal

Fix seven issues observed during real LDH identification runs and ad-
hoc mechanism construction:

1. Mechanisms emerge from enumeration with regulators listed in the
   type but never bound by any step.
2. `AllostericEnzymeMechanism` display of `cat_allo_states` is
   ambiguous about which group each tag refers to.
3. `identify_rate_equation` save files are named by an upper-bound
   estimate that disagrees with the `n_params` column in the rows.
4. A single saved file contains rows with several different actual
   `n_params` values.
5. `beam_fraction` selects mechanisms by rank rather than by closeness
   to the best loss; loss-threshold selection is wanted instead.
6. LDH identification produces ~3-4× redundant mechanisms — different
   structural specs that compile to the same rate equation.
7. `@enzyme_reaction` accepts atom-imbalanced reactions silently.

---

## Scope

In: src changes to `types.jl`, `dsl.jl`, `mechanism_enumeration.jl`,
`identify_rate_equation.jl`; matching test additions; updates to any
existing test/README example that becomes invalid under the new atom
requirement.

Out: optimizer changes (the kcat redundant direction noted during
brainstorming is left for a future spec); spec-level rate-equation
dedup (decided against — see §5).

---

## 1. Strict regulator validation (Issue #1)

### Behavior change

The 2-arg `EnzymeMechanism(mets, rxns)` constructor in
`types.jl:70-153` currently allows regulators in `mets[3]` to be
absent from every step (see lines 127-138 and the comment "Regulators
are optional bindings"). This produces types like

    EnzymeMechanism{((:S,), (:P,), (:A,)),
        (((:E, :P), (:E_P,), true, 1),
         ((:E, :S), (:E_S,), true, 2),
         ((:E_S,), (:E_P,), false, 3))}

where `:A` is in the regulator slot but never appears in any
reaction step.

**New rule:** every name in `mets[3]` must appear in some step. If
not, error with a message naming the offending regulator(s).

### Spec → type path

`EnzymeMechanism(spec::MechanismSpec; ...)` in
`mechanism_enumeration.jl:1286-…` currently builds the regulator
tuple via

    regs = Tuple(r for r in regulators(rxn) if r ∉ auto_exclude)

regardless of whether `r` appears in any of `spec.steps`. **Change
to:** intersect with the set of names that actually appear (with or
without the `__reg` suffix) on any step's reactants/products. The
strict constructor then accepts the result without erroring.

The allosteric `AllostericEnzymeMechanism(spec::AllostericMechanismSpec)`
path goes through `EnzymeMechanism(spec.base)` for the catalytic
mechanism — the same fix covers it. Allosteric regulators bound at
regulatory sites are unaffected (they live in `RegSites`, not in
`mets[3]`).

### Tests (TDD — write first, then implement)

In `test_types.jl`:

- 2-arg constructor errors when a regulator is unbound in steps.
- 2-arg constructor accepts a mechanism where all regulators bind.
- 2-arg constructor accepts an empty regulator tuple.

In `test_mechanism_enumeration.jl`:

- Initial mechanisms emitted by `init_mechanisms(reaction)` for a
  reaction with regulators do NOT include those regulators in the
  resulting `EnzymeMechanism` type when no expansion has added a
  binding step for them.
- After `_expand_add_dead_end_regulator` adds a binding step, the
  regulator is present in the resulting type.

---

## 2. Atom balance & mandatory atoms in `@enzyme_reaction` (Issue #7)

### Behavior change

In `EnzymeReaction(subs, prods, regs; oligomeric_state)`
(`types.jl:17`):

1. Every substrate and every product MUST carry a non-empty atom
   tuple. Reject (with a clear error naming the offending species)
   any whose `atoms` field is empty.
2. Sum atom counts element-by-element across substrates and across
   products.
3. Error if any element appears with different totals on the two
   sides, or appears on only one side.

Regulators are unaffected — atoms are dropped at parse time in
`dsl.jl:130`.

### Macro and parsing

`@enzyme_reaction` already supports atom brackets like `S[C6H12O6]`.
The atom mandatory-ness is enforced inside the runtime constructor —
no changes needed in `dsl.jl` parsing.

### Migration

Existing tests, README examples, and any in-repo demo that uses bare
substrate or product symbols in `@enzyme_reaction` must be updated
to declare atoms. This is part of the work for this spec.

Audit checklist:

- `test/mechanism_definitions_for_test_enzyme_derivation.jl`
- `test/test_dsl.jl`
- `test/test_types.jl`
- `test/test_accessors.jl`
- `test/test_fitting.jl`
- `test/test_identify_rate_equation.jl`
- `test/test_readme_runs.jl`
- `README.md`

For each example, choose plausible atom annotations (do not invent
biochemistry — copy from real reactions where the example mirrors
one, otherwise pick a minimal balanced placeholder like `[C]` on
both sides if the example is purely structural).

### Tests (TDD)

In `test_dsl.jl` and/or `test_types.jl`:

- `EnzymeReaction(((:S, ()),), ((:P, ((:C, 1),)),))` errors —
  substrate has no atoms.
- `EnzymeReaction(…)` errors when LHS has `(C, 6)` and RHS has
  `(C, 5)`.
- `EnzymeReaction(…)` errors when LHS declares carbon but RHS
  doesn't list it.
- `EnzymeReaction(…)` accepts a balanced 6C-6C reaction.
- A multi-substrate, multi-product example balances correctly when
  each side's element totals match.

---

## 3. `AllostericEnzymeMechanism` display refactor (Issue #2)

### Behavior change

Replace the current `Base.show(io, m::AllostericEnzymeMechanism)`
implementation in `types.jl:426-444`:

- Drop the standalone `cat_allo_states: [...]` line entirely.
- For the catalytic mechanism, force multi-line display (do not use
  the chain shortcut even if the topology is linear).
- Group steps that share a `kinetic_group` and emit them with a
  single tag using the `@allosteric_mechanism` macro syntax:
  - Single-step group:
    `  [E, S] ⇌ [E_S] :: EqualRT`
  - Multi-step group:
    `  ([E, S] ⇌ [E_S], [E_P, S] ⇌ [E_PS]) :: EqualRT`
  - SS step in own group:
    `  [E_S] <--> [E_P] :: OnlyR`
- Regulator-site rendering (existing `name::Tag` per ligand) is
  unchanged.

### Helper

Add `_format_allo_step_groups(cm::EnzymeMechanism, m::AllostericEnzymeMechanism)`:

1. Collect kinetic groups in their first-appearance order.
2. For each group, gather all steps with that group, render each
   step body (`[lhs] ⇌ [rhs]` / `[lhs] <--> [rhs]`).
3. If the group has one step → emit `body :: Tag`.
4. If multiple → emit `(body1, body2, …) :: Tag`.

The tag is `cat_allo_state(m, g)`.

### Tests (TDD)

In `test_types.jl` (or `test_dsl.jl`'s pretty-print section):

- `repr(m)` for an allosteric mechanism with one step per group
  contains `[E, S] ⇌ [E_S] :: EqualRT`.
- `repr(m)` for an allosteric mechanism with two steps in one group
  contains `(…, …) :: EqualRT`.
- `repr(m)` does NOT contain the substring `cat_allo_states:`.
- `repr(m)` for an allosteric mechanism whose catalytic mechanism
  has a linear topology still uses one step per line (no chain).

---

## 4. `identify_rate_equation` file naming and beam selection (Issues #3, #4, #5)

### File naming and `n_params` (Issues #3 + #4)

Currently `_save_level_csv` in `identify_rate_equation.jl:233-241`
names files by the spec's upper-bound `param_count` estimate, while
the row's `n_params` column is `length(fitted_params(m))` — the
actual independent rate-constant count. The two diverge for two
reasons:

- The estimate adds `+2` for `Keq` and `E_total`, neither of which
  is a fitted parameter.
- The estimate over-counts when Haldane reduction eliminates more
  parameters than the formula assumes.

**Fix has two parts:**

**Part A — rename and re-define `param_count` so its semantics
match the row's `n_params` column.**

The field `param_count` on `MechanismSpec` and
`AllostericMechanismSpec` is renamed to `n_fit_params_estimate`. Its
value drops the `+2` for `Keq` and `E_total`, so the estimate is the
predicted count of *independent rate constants only*, the same
quantity reported in the row's `n_params` column post-fit. The
formula in `init_mechanisms` (`mechanism_enumeration.jl:841`) becomes

    n_fit_params_estimate = n_re + 2 * n_ss - n_thermo

instead of the current

    param_count = n_re + 2 * n_ss - n_thermo + 2

Every reference to `.param_count` in `mechanism_enumeration.jl`
(~70 occurrences across src + tests) and the helper
`_param_count(spec)` / `_param_count_from_steps(steps)` are
renamed accordingly: `_n_fit_params_estimate(spec)`,
`_n_fit_params_estimate_from_steps(steps)`. Internal cache
keys and dictionary types
(`Dict{Int, Vector{AbstractMechanismSpec}}`) keep their `Int`
type — they're still bucketing by the renamed field's value.

The `expand_mechanisms` function and every `_expand_*` move that
adds to the count keeps its delta arithmetic — the deltas are
relative to the renamed field and don't change. A move that
previously bumped `param_count` by `+1` still bumps
`n_fit_params_estimate` by `+1`.

CLAUDE.md note about `param_count` being upper-bound is updated to
reference the new name and to clarify that it estimates fit
parameters only (matching the final `n_params` column convention).

**Part B — bucket save files by *actual*
`length(fitted_params(m))`.**

The estimate-vs-actual gap remains (Haldane reduction can still
collapse declared groups), so file naming continues to be by actual
count, not by the estimate. This naturally fixes #4:

- After the per-spec compilation+fit, group results by
  `r.row.n_params` (which comes from `fitted_params(m)`).
- For each group, call `_save_level_csv` with that group's rows and
  that param count.

The result: filename `params_3.csv` contains rows where every
`n_params == 3`, where `3` matches both the saved column AND the
estimate semantics declared at enumeration time (modulo Haldane
reduction collapsing it lower).

The kcat-induced redundant fit direction (one degree of optimizer
freedom is degenerate due to per-group centering) is **not**
subtracted from `n_params`; the count remains raw. A future spec can
revisit constraining kcat during fitting.

### Beam selection (Issue #5)

Replace the current rank-based logic in
`identify_rate_equation.jl:316-326`:

    perm = sortperm([r.row.loss for r in results])
    beam_size = max(ceil(Int, beam_fraction * length(results)),
                    min_beam_width)
    beam_size = min(beam_size, length(results))
    beam_specs = [results[perm[i]].spec for i in 1:beam_size]

with a loss-threshold rule:

- Compute `best_loss = minimum(loss_i)` over results that
  successfully fit.
- Mechanisms qualify for the beam if
  `loss ≤ loss_rel_threshold * best_loss + loss_abs_threshold`.
- Floor: always keep at least `min_beam_width` mechanisms (sorted by
  loss ascending; fill from the bottom of the sorted list).

Drop the `beam_fraction` keyword. New keyword arguments on
`identify_rate_equation`:

- `loss_rel_threshold::Float64 = 2.0`
- `loss_abs_threshold::Float64 = 0.01`

`min_beam_width::Int = 50` (default lowered from 200 — with the
loss-threshold rule now doing most of the trimming, a smaller floor
is appropriate; the previous 200 floor existed because rank was the
only mechanism for capping beam size).

The naming follows the convention used by ODE solvers and
optimization-error tolerances (`reltol`/`abstol`,
`rel_tol`/`abs_tol`) so the meaning is recognizable without
re-reading the docs.

The `identify_rate_equation` docstring must explicitly include the
selection formula:

    A mechanism qualifies for the next-level beam if either:
      • its loss ≤ loss_rel_threshold * best_loss + loss_abs_threshold,
      • OR its rank by loss is ≤ min_beam_width (the floor).

with a one-line note explaining why the additive term exists (so
that simulated / very-low-loss data, where best_loss can approach 0,
still admits structurally similar mechanisms within
`loss_abs_threshold` rather than collapsing to the single best).

### Tests (TDD)

In `test_mechanism_enumeration.jl`:

- For a small reaction, compute `init_mechanisms(reaction)` and
  assert every spec's `n_fit_params_estimate` equals
  `length(fitted_params(EnzymeMechanism(spec)))` for the simplest
  cases (where Haldane reduction does not collapse the declared
  groups further).
- Construct a spec where Haldane is known to collapse one group;
  assert `n_fit_params_estimate ≥ length(fitted_params(m))` (upper
  bound holds) and that the gap matches the expected number of
  collapsed groups.
- For an expanded spec (`expand_mechanisms` output), assert the
  expansion's delta on `n_fit_params_estimate` matches the move's
  documented `+1` / `+2` / etc. — same delta arithmetic as before
  the rename.

In `test_identify_rate_equation.jl`:

- Mock fit results with deliberate `n_params` values matching what
  `fitted_params(m)` would return; assert the saved CSV filename
  matches the actual count, not the estimate. (Use a tempdir.)
- Save two specs whose actual `n_params` differs but whose estimate
  is the same; assert they land in different files.
- Beam selection with `loss_rel_threshold=2.0`,
  `loss_abs_threshold=0.0`: a result with `loss = 2.5*best` is
  excluded; one with `loss = 1.9*best` is included.
- Beam selection with `loss_abs_threshold=0.01`: when `best_loss` is
  near zero (e.g. `1e-6`), a result with `loss = 0.005` is included
  even though it's 5000× the best (additive term saves it).
- `min_beam_width=5` floor: with only 3 results, all 3 are kept
  even if some are far above the threshold.
- `beam_fraction` removed: passing the kwarg raises
  `MethodError` (Julia's standard behavior for unknown kwargs).
  No alias / deprecation shim; this is an internal-version
  cleanup, not a public-API deprecation.

---

## 5. Post-compile rate-equation dedup (Issue #6)

### Background — what we found

LDH identification analysis showed ~3-4× redundancy in saved
mechanisms. The dominant pattern: structurally distinct specs that
all compile to identical rate equations. Two sources:

- **Pattern A (~80%):** alternative-binding-order RE edges. Multiple
  specs have the same form set, same kinetic-group structure, same
  SS steps, and differ only in which RE edges to doubly-bound forms
  are included. Under the rapid-equilibrium approximation, the RE
  concentration of each form is path-independent —
  `[F] = [E] · ∏(M_i / K_i)` over bound metabolites — so any
  spanning subset of the RE edges yields the same equation.
- **Pattern B:** specs declare distinct kinetic groups that
  Wegscheider closure forces to be equal. The reduced equation has
  fewer independent K's than the spec declared. Different "which
  declared groups end up forced equal" patterns collapse to the
  same equation.

### Decision: dedup at fit, NOT at enumeration

The 9 alternative-edge variants in Pattern A are equivalent under
their current all-RE-plus-one-SS form, but DIFFERENT once an
expansion converts an RE step to SS — the resulting rate equation
depends on which path was kept (since `[E_Lactate] ≠ [E_NAD]` and
the SS rate uses one of them as the source-form concentration).
Removing duplicates at enumeration time would lose downstream
coverage of these RE→SS variants.

So: dedup *fitting*, keep all specs for *expansion*.

### Implementation — global rate-equation hash cache

The cache is **persistent across all beam levels** so that a spec at
a higher param-count level whose Haldane-reduced rate equation
matches one fit at an earlier level reuses the prior fit instead of
duplicating it. This catches Pattern B (Haldane collapse to a
previously seen equation) for free, in addition to Pattern A
within-level duplicates.

The cache lives in `_beam_search` outside the level loop:

    fit_cache = Dict{
        UInt256,                    # full SHA-256 of canonical text
        @NamedTuple{
            loss::Float64,
            params::NamedTuple,
            first_seen_n_params::Int,
            first_seen_eq_hash::String,  # 8-char display hash
        }
    }()

Per-level processing in four stages — keep parallelism by batching
compile and fit separately:

1. **Compile + hash all specs at this level (parallel via
   `pmap_function`).** Each task returns
   `(spec, mechanism, eq_hash_full::UInt256, eq_hash_short::String,
   n_actual::Int)`.
2. **Bucket specs by `eq_hash_full` within the level.** Serial,
   cheap.
3. **Identify hashes that are NEW relative to `fit_cache`.** Fit
   one representative per new hash in parallel via `pmap_function`.
   Insert each result into `fit_cache`.
4. **Build CSV rows — one row per (within-level) hash group.**
   Columns:
   - existing: `n_params, loss, mechanism_type, rate_equation,
     fitted_param_names, fitted_param_values`.
   - new: `eq_hash::String` (8-char), `n_equivalent::Int`
     (within-level group size),
     `fit_inherited_from_n_params::Union{Int, Missing}` (`missing`
     if first-fit at this level; the originating `n_params`
     otherwise — diagnostic for Pattern B / Haldane collapse).

5. **Pass all specs (every member of every hash group, both newly
   fit and inherited) to `expand_mechanisms`.** Structure matters
   for downstream RE→SS expansion regardless of fit-cache hits.

6. In `_cv_model_selection`, dedup candidates by `eq_hash_full`
   globally — within each `n_params` bucket, run LOOCV once per
   unique hash, since identical rate equations give identical CV
   scores.

### Why not pure serial iteration?

Denis's original phrasing was "iteratively when a particular
mechanism is about to be fit, check if the same mechanism was
already fit." A pure serial loop would lose `pmap` parallelism for
fits — the dominant cost. The four-stage batch above preserves the
"check cache before fitting" semantics while keeping fits
parallelizable: hashes are deduplicated before the parallel fit
launches, so no two workers ever fit the same equation.

### Canonical rate-equation hash

`rate_equation_string(m)` returns a multi-line string with the
parameter destructure, intermediate dependent-parameter lines (e.g.,
`k10r = (1 / Keq) * …`), and the `v = …` formula. The canonicalizer
must rename consistently across all of these:

1. Drop the `(; … ) = params` and `(; … ) = concs` destructure
   lines (they are derived from the body and can vary in ordering).
2. Rename parameter symbols by *first-appearance order* in the
   canonical text:
   - SS forward rate constants: `k{N}f → kf_1, kf_2, …` in the
     order they first appear; same for `k{N}r → kr_1, kr_2, …`.
     Always indexed (no special-case "drop the index when only one
     exists"), so the rule is uniform regardless of mechanism size.
   - Equilibrium constants: `K{N} → K_1, K_2, …` in first-appearance
     order, independent of the `k`-renumbering.
   - `Keq`, `E_total`, and metabolite names are NOT renamed.
3. Whitespace-normalize (collapse runs of whitespace; trim).
4. Hash via SHA-256, take first 8 hex chars for `eq_hash` (display)
   plus the full 256-bit value for the dedup key.

The renaming is order-sensitive; the canonical text must be produced
deterministically (a single pass through the body, replacing each
unseen parameter with the next available canonical name).

### Tests (TDD)

In `test_identify_rate_equation.jl`:

- Two hand-constructed specs with the same form set, kinetic-group
  structure, and SS step but different RE edge subsets compile to
  rate-equation strings whose canonical hashes are equal. (This is
  the regression test for Pattern A.)
- A third spec whose kinetic-group structure has a Wegscheider
  redundancy hashes to the same value as a spec without the
  redundancy. (Pattern B.)
- Two specs with different form sets hash to different values.
- Two specs differing only in `k{N}f` vs `k{N+1}f` parameter names
  (representative-step indexing changing) hash to the same value.
- Mock `_beam_search` over a small reaction with known duplicates:
  - The number of fits performed equals the number of distinct
    hashes, not the number of specs.
  - Every spec from every hash group is forwarded to
    `expand_mechanisms` (use a recording mock for
    `expand_mechanisms`).
  - The saved CSV has one row per hash with `eq_hash`,
    `n_equivalent`, and `fit_inherited_from_n_params` populated.
- **Cross-level cache regression (the new test for Denis's
  refinement):** construct two specs at different beam levels (i.e.,
  different `n_fit_params_estimate`) that compile to the same
  canonical rate-equation hash. Run `_beam_search` and assert:
  - Only one fit is performed across the two levels (use a
    recording wrapper around `fit_rate_equation`).
  - The level-2 row has
    `fit_inherited_from_n_params == <level-1 n_params>`; the
    level-1 row has `missing` in that column.
  - The global cache survives the level transition (no per-level
    reset).
- **Within-level no-inheritance test:** two specs at the same level
  with the same hash both have `fit_inherited_from_n_params ==
  missing` (first-fit at this level, neither inherits).
- LOOCV regression: candidate count entering LOOCV equals the count
  of distinct `eq_hash` among top-N-by-loss-per-`n_params`.

---

## Cross-cutting test additions

Beyond per-fix tests, add an **end-to-end LDH regression test** in
`test_identify_rate_equation.jl` that runs `identify_rate_equation`
on a small LDH-like reaction (small enough to finish in the test
time budget) and asserts:

- The saved files' filenames match the `n_params` column inside
  every row.
- The canonical hash is unique within each saved file (no duplicate
  rate equations leak through).
- The canonical hash is unique **across all saved files combined**
  in the sense that any hash appearing in level-N also appears in
  exactly one earlier level (when present at all) — i.e., the
  cross-level cache is enforcing global dedup.
- `n_equivalent` summed across all rows in a file ≥ the row count
  (hash-group multiplicities are recorded).
- Beam threshold honoured: in a contrived run with deliberately
  wide loss spread, no row above threshold and rank `> min_beam_width`
  is forwarded to the next level.

---

## File-by-file change inventory

- `src/types.jl` — atom-balance check + mandatory atoms in
  `EnzymeReaction`; strict regulator check in `EnzymeMechanism`;
  `Base.show` refactor for `AllostericEnzymeMechanism`.
- `src/mechanism_enumeration.jl` — derive regulators from steps in
  `EnzymeMechanism(spec)` and `AllostericEnzymeMechanism(spec)`;
  rename `param_count` field on `MechanismSpec` and
  `AllostericMechanismSpec` to `n_fit_params_estimate` and drop the
  `+2` for `Keq` + `E_total` from the formula in `init_mechanisms`;
  rename helpers (`_param_count` → `_n_fit_params_estimate`, etc.)
  and update all expansion-move delta arithmetic to use the new
  field name (deltas themselves unchanged).
- `src/identify_rate_equation.jl` — new beam selection rule, drop
  `beam_fraction`, add `loss_rel_threshold` /
  `loss_abs_threshold`; bucket-by-actual-`n_params` save logic;
  post-compile hash dedup; LOOCV per unique hash.
- `src/dsl.jl` — no behavior change (atom mandatoriness lives in the
  runtime constructor); update inline doctests/examples that became
  invalid.
- `test/*.jl` — regression tests per fix as detailed above; update
  examples that used bare-symbol reactions.
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` —
  audit and update reaction definitions to declare atoms.
- `README.md` — update any reaction example.

---

## Out of scope (noted for future work)

- **Constraining `kcat` during fitting.** The optimizer currently
  has a redundant degenerate direction (uniform SS-k scaling)
  absorbed by per-group centering and post-hoc rescaling. A future
  spec could replace the post-hoc rescale with a true constraint,
  reducing the search dimension and making `n_params − 1` the
  effective complexity measure.
- **Spec-level enumeration dedup.** The K-multiset / Haldane-reduced
  signature was discussed and decided against (Q7 in brainstorming):
  it would lose downstream RE→SS expansion coverage.
- **Compilation cost optimization.** If post-compile dedup
  diagnostics (the `n_equivalent` column) reveal compile time is a
  real bottleneck, a future spec can revisit a partial spec-level
  pre-filter. Until then, observed compile time of ~1 s/mechanism
  is acceptable.
