# Codebase cleanup after the concrete-types + structural-parameter-names refactor

Date: 2026-05-31
Branch: refactor-to-concrete-types-instead-of-symbols

## Goal

The branch transitioned the package to concrete structs (`Mechanism` /
`AllostericMechanism` / `Step` / `Species`) and to structural parameter names
(`:K_ATP_E` instead of `:K1`). Documentation, code comments, memory, and the
docs tree still describe the *in-progress* state: stage/phase labels, deleted
helper functions, a "future" parameter-naming refactor that is now done, and a
large `docs/superpowers/` tree of execution plans for work that is complete.

This cleanup makes the repo's prose and memory describe the code **as it is
now**, with no references to prior code status, deleted files, or plan-execution
stages. It deletes refactor-era artifacts and the tooling that only existed to
police the refactor.

This is documentation/comment/memory hygiene plus targeted dead-tooling
removal. It changes **no production code paths** and must leave the full test
suite green.

## Scope decisions (settled with Denis)

1. **docs/superpowers/**: delete the whole tree; rescue only the 2 specs that
   describe *pending future work* into a new kept `docs/design/` folder.
2. **Operational scripts** (`check_test_integrity.sh`, `test_timing_report.sh`)
   and their data files existed only to enforce the refactor's test-deletion
   discipline. The refactor is done → **retire the scripts and their data files
   entirely**.
3. **Memory**: prune memories that only captured now-finished migration state;
   keep durable invariants and strip stage/date language from them. Denis
   reviewed the per-memory table below.
4. **CLAUDE.md parameter-naming section**: rewrite to the structural reality the
   code emits (verified against `name()` output), not the old `:K1` index names.

## Out of scope

- README.md — already clean (one benign "the one we used to generate it").
- Any production-code behavior change.
- `test_types.jl:654` R/T-rejection loop — **investigated, verified correct**.
  It is a `@test_throws` guard asserting the old `:OnlyR/:OnlyT/:EqualRT/
  :NonequalRT` symbols are *rejected* by `RegulatorySite`. Not stale. No change.

## Work streams

### Stream 1 — docs/superpowers/ deletion + rescue

- `git mv docs/superpowers/specs/2026-05-29-direction-symmetry-constraint-resolution.md docs/design/`
- `git mv docs/superpowers/specs/2026-05-29-nonequalai-rank-validity.md docs/design/`
- `git rm -r` the remainder of `docs/superpowers/` (all `plans/`, all other
  `specs/`, `handoff-2026-05-24-finish-refactor.md`,
  `2026-05-30-refactor-audit-findings.md`, `refactor-deleted-tests.md`,
  `refactor-test-timings-main-baseline.txt`).
- Remove the now-pointless `.gitignore` entry for
  `docs/superpowers/scratch-refactor-audit-notes.md`.
- End state of `docs/`: `docs/design/` holds the 2 rescued specs (+ this
  cleanup spec). Nothing else.

### Stream 2 — retire refactor-policing tooling

- `git rm scripts/check_test_integrity.sh scripts/test_timing_report.sh`.
- Their data files (`refactor-deleted-tests.md`, `refactor-test-timings-main-
  baseline.txt`) are removed in Stream 1.
- Verified: no CI workflow and no tracked file outside `docs/superpowers/`
  (being deleted) references either script.

### Stream 3 — CLAUDE.md

- **Rewrite the "Parameter naming convention" section.** It currently describes
  index names (`:K1`, `:K6`, `:k{rep}f`, `:K{rep}_T`). The value-context
  `name(p::Parameter, m)` now renders **structural** names via
  `_render_binding`/`_render_iso` (species-based, e.g. `:K_ATP_E`).
  - Before writing: compile a representative mechanism and read the actual
    `name()` output for binding, iso, SS-forward/reverse, and inactive-state
    parameters; also confirm what the type/index companion
    `name(::Type{P}, idx::Int)` emits today. Rewrite the section to match.
  - Update the dependent sections that repeat the `:K1` convention
    (`_dependent_param_exprs`, `rate_equation_string`, `_kcat_forward` prose).
- **Chokepoint section (L236):** delete "The future parameter-naming refactor
  (`:K1` → `:K_ATP`) changes a single function body." — that work is done.
  Keep the chokepoint description itself (still accurate).
- **L222:** delete the parenthetical "(The earlier singleton-typed
  `EnzymeReactionLegacy{S,P,R,N}` was retired in Stage 7d.1.)".
- **L265:** repoint the two follow-up spec paths from
  `docs/superpowers/specs/` → `docs/design/`.
- **L315:** delete "from Stage 6.1".
- Sweep for any remaining `Stage N` / `Phase N` / `PR #` anchors and remove the
  temporal framing while preserving the technical statement.

### Stream 4 — source + test comments

- `src/mechanism_enumeration.jl:1650-1652` — "the legacy parametric form … this
  file forwards from legacy": `EnzymeReactionLegacy` has 0 src references.
  Verify the boundary adapter is gone and rewrite the docstring to describe the
  current `EnzymeReaction`-only dispatch, with no "legacy" framing.
- `src/mechanism_enumeration.jl:1909, 1936` — "`_T` is a legacy cache-token
  literal from the old R/T allosteric naming": `_T` is the *current*
  inactive-state suffix (CLAUDE.md: `:K{rep}_T`). Rewrite to state what the
  suffix is and why, with no historical reference.
- `src/rate_eq_derivation.jl:992` — trim "that each previously rebuilt this set"
  to drop the history clause.
- Re-grep `src/` for any other `legacy` / `previously` / `used to` / `formerly`
  history comments surfaced during the pass and fix them in the same stream.
- `test/test_rate_eq_derivation.jl:980` — the cited spec
  (`2026-05-11-rate-eq-emission-perf-fix-design.md`) is deleted; drop the path,
  keep the n-ary-varargs boxing rationale.
- `test/test_rate_eq_derivation.jl:1413` — repoint the cited spec to
  `docs/design/2026-05-29-nonequalai-rank-validity.md` (rescued).

### Stream 5 — memory (Denis-reviewed dispositions)

| Memory | Action | Notes |
|---|---|---|
| `project-new-sig-lossy-for-opaque-species` | **Delete** | Built on `_legacy_step_tuple_from_sig` (deleted); opaque path removed. |
| `project-ss-dissociation-reconstruction-rule` | **Delete** | Built on `_step_tuple_from_sig` (deleted). |
| `project-lumped-complex-migrates-by-rename` | **Delete** | Migration done; cites deleted `_legacy_step_tuple` + stale `types.jl:1338`. Durable "rename preserves rate-eq" insight no longer actionable. |
| `project-derivation-backend-opaque-tuples-deferred` | **Rewrite** | Reduce to the 2 durable facts: (a) the `Sig`/`@generated` bridge stays (Denis-confirmed, infeasible to remove without losing the 0-alloc/<100ns gate); (b) the unresolved canonical-hash group-order vs `fitted_params` rep-naming tension. Drop the migration narrative and "deferred" framing. Re-title without "deferred". |
| `project_n_fit_params_estimate_undercounts` | **Rewrite** | Verify the SS under-count still holds in current code; update `:NonequalRT` → `:NonequalAI`; drop "deferred to Task 7d.3" and the broken `[[project-test-helpers-from-6beta-1]]` link. Keep the "assert against ground-truth `fitted_params(compile_mechanism(m))`" guidance. |
| `reference-test-invocation-and-integrity-baseline` | **Rewrite** | Remove the entire integrity-check section (scripts retired). Fix the broken per-file snippet (delete the `using EnzymeRates: EnzymeReactionLegacy` line — the symbol is gone). Keep: `Pkg.test()` works, `--project=test` does not, ~26.9k test count with CMA-ES ±1 fluctuation. |
| `feedback-test-translation-assertion-strength` | **Keep, de-temporalize** | Durable principle (don't weaken assertions; use ground-truth `fitted_params`). Drop the 6β task numbers / spec-§ refs; state the principle evergreen. |
| `project-iso-steps-zero-metabolite-stoich` | **Keep, de-temporalize** | Invariant still true (`_thermodynamic_constraints`/`_step_sides`). Drop the "deferred derivation-backend refactor" and 2026-05-30-plan references. |
| `project-metabolites-fitted-params-stay-generated` | **Keep, de-temporalize** | Hot-path invariant still true. Drop the "2026-05-30 refactor plan Cluster A / Q-005" framing; keep the invariant + why + how-to-apply. |

For every "Rewrite"/"Keep" memory: verify each `file:line` / function-name
citation against current code before re-asserting it (per the memory-staleness
rule).

### Stream 6 — MEMORY.md index

Rewrite the index to match the post-cleanup set: drop the 3 deleted entries,
update the one-line hooks for the 3 rewritten entries (fix `NonequalRT`
wording), leave the 3 kept entries' hooks accurate.

## Verification

- `git grep -n "Stage [0-9]\|Phase [0-9]\|\blegacy\b\|EnzymeReactionLegacy\|
  NonequalRT\|docs/superpowers"` over tracked files returns only legitimate
  hits (e.g. `identify_rate_equation.jl`'s genuine Stage-1/2/3 pipeline labels,
  the `test_types.jl` R/T-rejection guard).
- `git grep -n "docs/superpowers"` returns nothing.
- The 2 rescued spec paths resolve under `docs/design/`.
- Full suite green: `julia --project -e 'using Pkg; Pkg.test()'` (the doc/comment
  edits touch no code paths; this confirms the comment edits didn't disturb any
  string-shape/Expr-shape regression test, and that retiring the scripts left
  no dangling test reference).
- Memory files and `MEMORY.md` index are internally consistent (no dangling
  `[[links]]`, every index line points to an existing file).

## Commit shape

Small, reviewable commits along stream boundaries:
1. `docs: delete refactor-era docs/superpowers tree; rescue 2 live specs to docs/design`
2. `chore: retire refactor test-integrity + timing scripts`
3. `docs: rewrite CLAUDE.md parameter-naming section to structural names; drop stage anchors`
4. `docs: remove historical framing from src/test comments`
5. memory edits (committed per the memory workflow)
