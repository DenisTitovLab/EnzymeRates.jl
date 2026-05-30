# Concrete-Types Refactor Audit — Workflow Design

**Date:** 2026-05-30
**Branch:** `refactor-to-concrete-types-instead-of-symbols`
  (unpushed; main = 7,136 src LOC; this branch = 9,402 src LOC across 9 files)
**Status:** design pending Denis review

## Goal

Design a workflow that systematically investigates the entire codebase
(src/ plus blocking tests) line-by-line and produces a single findings
document capturing every simplification opportunity the concrete-types and
structural-parameter-names refactors enabled but did not complete. The
finished refactor branch should be substantially smaller than current —
Denis's hypothesis is that up to half of src may be removable. The audit's
job is to test that hypothesis with evidence.

## What counts as a finding (all four categories in scope)

- **Dead code / leftover legacy** — helpers with no production callers,
  validators that no longer fire, stale comments referencing removed paths.
- **Collapse duplicated logic** — the same algorithm or dispatch expressed
  in two or more places (A/I rename across 3 files, rep-step → param-name
  in multiple spots, symbol-tuple round-trips that regenerate Step/Species
  structure).
- **Architectural simplifications** — moves the refactor enabled but did
  not finish: demote `EnzymeMechanism{Sig}` / `AllostericEnzymeMechanism`
  to an internal compile artifact, replace symbol-tuple plumbing with a
  struct-native analysis context, fold compile-time accessors that aren't
  hot-path back into plain Mechanism field access.
- **Test-surface cleanup** — tests exercising private helpers (e.g.
  `_onlyA_parameters`, `_I_rename_parameters`), `name_map` string-key tests,
  accessor perf gates — but only where they block a src simplification.

## Hard constraints (from CLAUDE.md)

These are non-negotiable. Any finding that would violate one is rejected
or held until the constraint is renegotiated with Denis.

1. **`rate_equation` performance** — must remain 0 allocations and <100 ns
   per call across every spec in `MECHANISM_TEST_SPECS`. Enforced by
   `test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`.
2. **No net test-coverage reduction** — tests may be moved or replaced
   (behavior coverage replacing private-helper coverage), but not net
   deleted. Each deletion logged per the existing
   `docs/superpowers/refactor-deleted-tests.md` convention.
3. **No backward-compatibility shims** without explicit Denis approval
   (per CLAUDE.md "Writing code").
4. **Opaque bound-form rejection is load-bearing** — the
   `_assert_no_opaque_terms` / `_is_conformation_shape` guard in `dsl.jl`
   stays (per `2026-05-26-finish-refactor-legacy-sig-removal-design.md`,
   "Dropped: Bucket B′").
5. **Single struct family in the front end is already achieved** — the
   audit must not propose reintroducing parallel representations on the
   front-end (DSL + enumeration + public Mechanism surface).
6. **The derivation back-end's opaque-tuple round-trip is deferred** —
   `_species_name_from_sig`, `_step_tuple_from_sig`, and consumers in
   King-Altman / Wegscheider are the genuine next-phase refactor (per
   `project-derivation-backend-opaque-tuples-deferred` memory). The audit
   may surface architectural findings here, but they go in a separate
   cluster flagged as "next-phase" — not part of this branch's
   simplification scope unless Denis says otherwise.

## Workflow — three passes

### Pass 1 — Catalog (serial, by Claude)

Read every src file top-to-bottom, in `EnzymeRates.jl` inclusion order
(this is the natural data-flow order: types → DSL → derivation →
enumeration → identify). For each function / method / struct / macro,
flag suspects matching any of the table below.

| Pattern | Suspect criterion |
|---|---|
| Dead code | Non-exported symbol; visible call sites only in other suspect code, in removed-style branches, or in test-only "internal helper" assertions |
| Duplication | Same algorithm, regex, dispatch, or formula appears in 2+ functions or files |
| Symbol-tuple plumbing | Code converts Sig → opaque Symbol tuple → back, OR regenerates structure (substrate list, product list, kinetic-group rep, bound metabolite) that Step / Species / Mechanism already carries directly |
| Stale comment | Mentions "legacy", "old Sig", "previous path", "deprecated", "OLD:", "moved", "renamed from" |
| Compile-time accessor | `@generated` walking Sig where plain Mechanism field access would work. "Not on the hot path" means: not called inside `rate_equation`'s body or any function it inlines into. Verified by reading the `rate_equation` `@generated` body (the only function bound by `test_rate_equation_performance`); the accessor perf test `test/test_accessors.jl` is **explicitly negotiable** per the prior-agent audit |
| Test-private helper | A `_`-prefixed helper that has a direct test assertion against its return value (so the test constrains the helper's signature, blocking refactor) |
| String-keyed projection | `name_map` and related — projections built from rendered Symbol strings that structural parameter keys obsoleted (`src/identify_rate_equation.jl:300-390, 421-496`; `src/mechanism_enumeration.jl:1840-2030` per previous-agent audit) |

Suspects go in a scratch working file
`docs/superpowers/scratch-refactor-audit-notes.md` (NOT committed; added
to .gitignore for the audit). Each suspect row:
`<file>:<lines> | <category guess> | <one-line summary> | <on-encounter confidence H/M/L>`.

**On-encounter confidence criteria:**
- **H** — match is unambiguous from reading the local code (e.g., a helper
  whose only usage site is inside a `if false`-style dead branch, or a
  comment block literally saying "TODO: delete after Stage 7d").
- **M** — match is likely but verification is needed (e.g., a helper that
  *looks* uncalled from its file but I haven't grepped the codebase yet).
- **L** — match is a hunch (e.g., "this `@generated` function might be
  collapsible to field access but I'd need to compare it to
  `rate_equation`'s call graph").

No verification in Pass 1. The goal is breadth: capture every suspect; let
Pass 2 cull.

**Meta-completeness check at end of Pass 1.** Before handing off to Pass 2,
verify every line range of every src file is accounted for (read or
explicitly skipped with reason). The check is: for each file, walk the
scratch entries and confirm the union of `<lines>` ranges plus
"intentionally not flagged" regions covers `[1, file_LOC]`. This catches
silent omissions (a whole function block skipped because I lost my place).

### Pass 2 — Verify (parallel subagents)

Suspects are batched into verification queries dispatched to parallel
subagents. Typical query shapes:

- *"Find all callers of `<symbol>` in `src/` and `test/`; report
  `<file>:<line>` for each, separating definitions, call sites, and
  comment / string references."*
- *"Confirm `<func A>` is called only from `<file>`; report any external
  callers."*
- *"Compare bodies of `<func A>` and `<func B>` — list every line they
  share, every line that differs, and whether either has hidden
  side-effects."*
- *"Walk every `@generated` body in `src/rate_eq_derivation.jl` and
  classify: (a) builds the rate-equation Expr (hot path); (b) computes a
  derived structure that's reused at runtime; (c) exposes a Mechanism field
  through a Sig-walk that could be a plain `getfield`."*

Each verified suspect is promoted to a finding with the format below.

**Post-verification confidence (for the finding):**
- **High** — every claim in the recommendation is backed by grep / read
  evidence; no caller is unaccounted for; the LOC saving is from counted
  lines (not estimated).
- **Medium** — the recommendation rests on one judgment call (e.g.
  "these two functions are equivalent up to symbol renaming, which I'm
  asserting without running both"). Must be re-verified during impl.
- **Low** — the recommendation requires significant design choice (e.g.
  "demote this type — but the exact API needs design"). Treated as a
  pointer rather than a ready-to-implement change.

Suspects with partial evidence get demoted to Low-confidence findings;
suspects that fail verification are dropped (and the drop logged in the
scratch file with one line explaining why).

### Pass 3 — Synthesize (serial, by Claude)

- Cluster findings by shared architectural move (e.g., demoting
  `EnzymeMechanism{Sig}` touches 9-12 findings across 3 files — one cluster).
- Build the dependency graph: `F-B` depends on `F-A` iff `F-B`'s
  recommendation assumes `F-A` has landed (e.g., collapsing an accessor
  depends on demoting the type that exposes it).
- Sum LOC savings per cluster.
- For each finding, scan tests for blocking coverage (`test/...:LL-MM`).
- Integrate the previous agent's 5 findings: for each, mark
  - **Validated** — my audit produces ≥1 finding that covers the same
    recommendation, with file:line evidence the previous agent did not
    cite. Disposition: list the covering F-IDs.
  - **Refined** — my audit covers the recommendation but with a different
    scope, ordering, or implementation strategy. Disposition: list the
    covering F-IDs and one-sentence delta.
  - **Contradicted** — my audit finds the previous recommendation
    incorrect or unsafe (e.g., a guard the prev agent missed). Disposition:
    cite evidence; document the contradiction in §6.
  - **Superseded** — my audit replaces the previous recommendation with a
    materially different one (e.g., chose a different alternative).
    Disposition: cite the new F-ID and explain why.
  Carry each into the new numbering system rather than referring back to
  the prompt.
- Decide the suggested sequencing (see §"Findings Doc Structure" below).
- Write the findings doc.

## Finding format

```
#### F-NNN  <one-line title in imperative voice>

**Location:** src/<file>:<L1>-<L2> [, src/<file>:<L>-<L>, ...]
**Category:** Dead code | Duplication | Architectural | Test-surface
**Confidence:** High | Medium | Low
**LOC saving:** ~N (or "unknown; verify in impl")
**Depends on:** F-NNN, F-NNN (or "none")
**Blocking tests:** test/<file>:<L>-<L> (or "none"). A test is *blocking*
   if it would fail or no longer compile after the recommendation is
   applied — e.g., it asserts on the exact return value of a helper marked
   for deletion, or names a Symbol the rename would change. Tests that
   merely call public API are not blocking.
**Recommendation:** <one paragraph: what to remove / collapse / refactor,
   in plain Julia terms, and what behavior must be preserved>
```

## Findings doc structure

Output location: `docs/superpowers/2026-05-30-refactor-audit-findings.md`
(top-level `docs/superpowers/`, **not** `specs/` — this is a deliverable
output, not a design). Naming follows the existing pattern set by
`docs/superpowers/refactor-deleted-tests.md`.

```
# Concrete-Types Refactor Audit — Findings

**Date:** 2026-05-30
**Branch:** refactor-to-concrete-types-instead-of-symbols
**Baseline:** 9,402 src LOC across 9 files
**Workflow:** docs/superpowers/specs/2026-05-30-refactor-audit-workflow-design.md

## §1 Executive summary

- N total findings (D dead-code, Du duplication, A architectural, T test-surface)
- Estimated LOC savings: ~Z lines (P% of baseline)
- Hypothesis test: half-of-src claim is supported / partially supported / not supported
- Top 3 architectural moves (one-line each)

## §2 Findings (ordered by src file, then by line range)

### src/types.jl
#### F-001 ...
#### F-002 ...

### src/dsl.jl
...

[and so on for all 9 src files]

## §3 Dependency clusters

- **Cluster A — <name>**: F-NNN, F-NNN, F-NNN. Total ~N LOC.
- **Cluster B — <name>**: ...

## §4 Suggested sequencing

1. **First wave (no-deps, high-confidence dead code)** — F-..., F-..., ...
2. **Second wave (duplication collapse)** — Cluster B, Cluster C, F-...
3. **Third wave (architectural)** — Cluster A, ...
4. **Next-phase (out of scope for this branch)** — Cluster X (derivation
   back-end opaque-tuple removal)

## §5 Hard constraints tracked

- rate_equation perf budget: every finding touching the derivation marks
  whether it crosses the budget; impl plan must include perf evidence
- Test coverage: every finding deleting a test names the behavior-test
  replacement (or marks the helper as truly internal with no behavior loss)
- Opaque-form rejection: load-bearing; not in scope

## §6 Previous agent's 5 findings — disposition

- (Prev #1) Demote singleton Sig types → Cluster A (F-NNN..F-NNN). Validated.
- (Prev #2) Replace symbol-tuple plumbing with struct-native context → F-NNN. ...
- (Prev #3) Collapse parameter/state machinery → Cluster B (F-NNN..F-NNN). ...
- (Prev #4) Delete `name_map` projection → F-NNN. Refined: chose alternative ...
- (Prev #5) Simplify enumeration moves → Cluster F (F-NNN..F-NNN). ...
```

## Working notes file (throwaway)

`docs/superpowers/scratch-refactor-audit-notes.md` is the running Pass 1
catalog. It is added to `.gitignore` for the duration of the audit and
deleted after the findings doc lands. The findings doc is the committed
deliverable; the scratch is throwaway by design.

## Hypothesis testing

The half-of-src target (Denis's hypothesis) is tested at Pass 3. The
thresholds below are **heuristic** — they exist to trigger a conversation,
not to gate the audit:

- If estimated savings ≥ 40% of baseline: claim **supported**; proceed to
  impl plans.
- If 20%-40%: claim **partially supported**; impl plans proceed but
  Denis should consider whether further refactor passes are warranted.
- If < 20%: claim **not supported**; pause and discuss with Denis whether
  the audit missed something or whether the hypothesis was too optimistic.

The audit reports the honest number, not a contorted one (matches the
"LOC gate: renegotiate to measured" stance in
`2026-05-26-finish-refactor-legacy-sig-removal-design.md`).

## Out of scope

- Refactoring tests beyond those that block a src simplification
- Modifying `rate_equation` performance characteristics
- Replacing the `@generated`-driven King-Altman / Cha derivation
  (different refactor; deferred per `project-derivation-backend-opaque-
  tuples-deferred` memory)
- Pruning the test suite for general redundant coverage (only
  private-helper / structural-key tests are in scope)
- Implementing the simplifications themselves — that happens in
  follow-on impl plans, one per dependency cluster

## Sequencing after this spec

1. Denis reviews this spec
2. `writing-plans` skill produces the audit execution plan (the step-by-
   step "run Pass 1, verify, synthesize" plan)
3. Audit executes; findings doc lands
4. Subsequent impl plans (one per dependency cluster) drive the actual
   simplifications

## Success criteria for this audit

The audit succeeds if:

- Every src file has been read line-by-line and every suspect captured
- Every finding cites file:lines and a verifiable evidence chain
- Findings are organized for execution: clusters identified, dependencies
  named, blocking tests flagged, sequencing proposed
- The previous agent's 5 findings are each accounted for with disposition
- The half-of-src hypothesis is tested with the honest measured number

The audit does **not** need to find exactly half of src removable. The
honest measured number is the deliverable, whatever it is.
