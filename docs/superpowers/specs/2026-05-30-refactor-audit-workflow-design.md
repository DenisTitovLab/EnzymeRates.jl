# Concrete-Types Refactor Audit — Workflow Design

**Date:** 2026-05-30
**Branch:** `refactor-to-concrete-types-instead-of-symbols`
  (unpushed; main = 7,136 total src LOC; this branch = 8,456 total src
  LOC across 9 files, 5,706 non-comment non-doc — the audit's target
  metric — measured 2026-05-30 by `wc -l` and an awk filter dropping
  blanks, `#` comments, and `"""…"""` docstring blocks.)
**Status:** design pending Denis review

## Goal

Design a workflow that systematically investigates the entire codebase
(src/ plus blocking tests) line-by-line and produces a single findings
document capturing every simplification opportunity the concrete-types and
structural-parameter-names refactors enabled but did not complete.

**Dual goal:**

1. **Reduce non-comment non-doc LOC in src.** Denis's hypothesis is that up
   to half of src may be removable. The audit measures and reports the
   honest number against this hypothesis. The metric is non-comment
   non-doc LOC — comment and docstring LOC do not count toward the
   reduction target.
2. **Simplify the code.** Where two pieces of logic do the same thing,
   collapse them. Where a `@generated` accessor exists to walk a Sig that
   could be a plain `getfield`, collapse it. Where an indirection serves
   no purpose, remove it. Simplifications that do not reduce LOC still
   count (e.g., flattening a 50-line `@generated` body into a 50-line
   plain function with no Sig walk is a simplification win).

These two goals reinforce each other but are tracked separately: each
finding records both its LOC saving (goal 1) and a one-line
"simplification gain" (goal 2).

## What counts as a finding

Five categories, all in scope. The first four reduce non-comment non-doc
LOC; the fifth doesn't but is in scope as a simplification / hygiene pass.

- **Dead code / leftover legacy** — helpers with no production callers,
  validators that no longer fire, branches for paths that no longer exist.
- **Collapse duplicated logic** — the same algorithm or dispatch expressed
  in two or more places (A/I rename across 3 files, rep-step → param-name
  in multiple spots, symbol-tuple round-trips that regenerate Step/Species
  structure that the value-context already carries).
- **Architectural simplifications** — moves the refactor enabled but did
  not finish. The big ones include:
  - Demote `EnzymeMechanism{Sig}` / `AllostericEnzymeMechanism{...}` to
    an internal compile artifact (front-end already on one struct family
    per `2026-05-26-finish-refactor-legacy-sig-removal-design.md` — still
    investigate whether the singleton-type bridge can collapse further).
  - **Replace the derivation back-end's symbol-tuple round-trip with a
    struct-native analysis context.** `_species_name_from_sig`,
    `_step_tuple_from_sig`, and the King-Altman / Wegscheider consumers
    in `rate_eq_derivation.jl` and
    `thermodynamic_constr_for_rate_eq_derivation.jl` regenerate opaque
    Symbol tuples that the `Mechanism` value already carries. This is
    **in scope for this audit** (previously deferred per the
    `project-derivation-backend-opaque-tuples-deferred` memory — Denis
    reopened it).
  - Fold compile-time accessors that aren't on the `rate_equation`
    hot path back into plain Mechanism field access.
  - Strict-parser-replaces-guard: `_assert_no_opaque_terms` is a
    post-hoc patch over a permissive parser at `dsl.jl:259-262`. Refusing
    non-conformation-shaped bare Symbols at the original parse site
    makes the guard redundant — investigate.
- **Test-surface cleanup** — tests exercising private helpers (e.g.
  `_onlyA_parameters`, `_I_rename_parameters`), `name_map` string-key
  tests, accessor perf gates — but only where they block a src
  simplification.
- **Comment / doc hygiene** *(does not reduce non-comment non-doc src LOC;
  in scope as a simplification / readability pass)*:
  - Comments referencing previous specs, stages, or phases (e.g. "from
    Stage 7d", "see 2026-05-26-finish-refactor", "Phase 2 enumerator")
    that no longer point to live design context — flag for removal.
  - Block comments used where a `"""docstring"""` should be (e.g. a
    function explained by a `# ABOUTME` line above it instead of a
    docstring on the function itself) — flag for conversion to a
    proper docstring.
  - Stale "OLD:" / "renamed from" / "previously called X" / "moved from
    Y" markers — flag for removal per CLAUDE.md's "evergreen comments"
    rule.

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
4. **The front end remains on one struct family** — the
   `2026-05-26-finish-refactor-legacy-sig-removal-design.md` work
   unified DSL + enumeration + the public `Mechanism` surface onto one
   struct family. Audit findings may simplify that family further but
   must not reintroduce parallel front-end representations.

## Current state — facts, not constraints

The following are observations about the post-refactor state, **not**
guards on the audit. They orient the auditor; findings that simplify or
reverse them are legitimate.

- **Opaque-form rejection guard at `dsl.jl:567,1130`**
  (`_assert_no_opaque_terms`) — currently rejects bare-Symbol bound forms
  like `:ES`. The previous spec called this "load-bearing"; the audit
  should test that claim. The guard is a post-hoc walk that duplicates a
  check the parser at `dsl.jl:259-262` could enforce directly. If the
  parser is tightened, the guard can be deleted.
- **Derivation back-end still uses opaque Symbol tuples** —
  `_species_name_from_sig` (`types.jl:~1408`) and `_step_tuple_from_sig`
  rebuild opaque `(lhs, rhs, is_eq, g)` Symbol tuples at `@generated`
  time, and the King-Altman / Wegscheider consumers
  (`rate_eq_derivation.jl:142-143,445`,
  `thermodynamic_constr_for_rate_eq_derivation.jl:88-90,207-209`) match
  enzyme forms by Symbol identity. **This is now in scope for the
  audit** (Denis reopened it).

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
| Compile-time accessor | `@generated` walking Sig where plain Mechanism field access would work. "Not on the hot path" means: not called inside `rate_equation`'s body or any function it inlines into. Verified by reading the `rate_equation` `@generated` body (the only function bound by `test_rate_equation_performance`); the accessor perf test `test/test_accessors.jl` is **explicitly negotiable** per the prior-agent audit |
| Test-private helper | A `_`-prefixed helper that has a direct test assertion against its return value (so the test constrains the helper's signature, blocking refactor) |
| String-keyed projection | `name_map` and related — projections built from rendered Symbol strings that structural parameter keys obsoleted (`src/identify_rate_equation.jl:300-390, 421-496`; `src/mechanism_enumeration.jl:1840-2030`) |
| Permissive parser + post-hoc guard | A parser branch that accepts a wider input language than intended, paired with a later validator that rejects the slack. Collapsing the parser to be strict makes the guard redundant. Example: `_assert_no_opaque_terms` over the bare-Symbol branch at `dsl.jl:259-262` |
| **Stale spec/stage comments** *(doc category)* | Comment mentions "legacy", "old Sig", "previous path", "deprecated", "OLD:", "moved", "renamed from", "from Stage Nx", "see YYYY-MM-DD-...", "Phase N", "per <past spec doc>" — flag for **removal**, per CLAUDE.md evergreen-comments rule. Does not reduce non-comment non-doc LOC |
| **Comment used as docstring** *(doc category)* | A function or struct whose explanation lives in a block of `#`-comments immediately above it (often an `ABOUTME` line + several `#`-explanation lines) rather than a `"""docstring"""` attached to the definition. Flag for **conversion** to docstring. Does not reduce non-comment non-doc LOC |

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
- Decide the suggested sequencing (see §"Findings Doc Structure" below).
- Write the findings doc.

The previous agent's 5 findings (cited in Denis's prompt) are **reference
material only**. The audit produces its own findings from scratch with
full file:line evidence; no integration / disposition step is required.

## Finding format

```
#### F-NNN  <one-line title in imperative voice>

**Location:** src/<file>:<L1>-<L2> [, src/<file>:<L>-<L>, ...]
**Category:** Dead code | Duplication | Architectural | Test-surface | Doc hygiene
**Confidence:** High | Medium | Low
**LOC saving (non-comment non-doc):** ~N (or "0; doc category" or
   "unknown; verify in impl")
**Simplification gain:** <one line: what the code does more clearly /
   directly after the change, even if LOC is unchanged>
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
**Baseline (non-comment non-doc src LOC):** measured at Pass-3 start
**Workflow:** docs/superpowers/specs/2026-05-30-refactor-audit-workflow-design.md

## §1 Executive summary

- N total findings (D dead-code, Du duplication, A architectural,
  T test-surface, H doc-hygiene)
- Estimated savings: ~Z non-comment non-doc src LOC (P% of baseline)
- Simplification themes (one-line each): <top 3-5 architectural moves
  and what they replace>
- Hypothesis test: half-of-src claim is supported / partially supported /
  not supported

## §2 Findings (ordered by src file, then by line range)

### src/types.jl
#### F-001 ...
#### F-002 ...

### src/dsl.jl
...

[and so on for all 9 src files]

## §3 Dependency clusters

- **Cluster A — <name>**: F-NNN, F-NNN, F-NNN. Total ~N LOC + simplification.
- **Cluster B — <name>**: ...

## §4 Suggested sequencing

1. **First wave (no-deps, high-confidence dead code)** — F-..., F-..., ...
2. **Second wave (duplication collapse)** — Cluster B, Cluster C, F-...
3. **Third wave (architectural)** — Cluster A, derivation-back-end cluster, ...
4. **Doc-hygiene sweep** — all doc-category findings batched into one
   commit (no behavior change)

## §5 Hard constraints tracked

- rate_equation perf budget: every finding touching the derivation marks
  whether it crosses the budget; impl plan must include perf evidence
- Test coverage: every finding deleting a test names the behavior-test
  replacement (or marks the helper as truly internal with no behavior loss)
- Front-end struct-family unification: must not be reintroduced as parallel
```

## Working notes file (throwaway)

`docs/superpowers/scratch-refactor-audit-notes.md` is the running Pass 1
catalog. It is added to `.gitignore` for the duration of the audit and
deleted after the findings doc lands. The findings doc is the committed
deliverable; the scratch is throwaway by design.

## Hypothesis testing

The half-of-src target (Denis's hypothesis) is measured against
**non-comment non-doc src LOC**, baseline computed at Pass 3 start (a
plain count, e.g. `grep -vE '^\s*(#|""")' src/*.jl | wc -l`-style with
a Julia-aware filter; the exact command goes in the impl plan).

The thresholds below are **heuristic** — they exist to trigger a
conversation, not to gate the audit:

- If estimated savings ≥ 40% of baseline: claim **supported**; proceed to
  impl plans.
- If 20%-40%: claim **partially supported**; impl plans proceed but
  Denis should consider whether further refactor passes are warranted.
- If < 20%: claim **not supported**; pause and discuss with Denis whether
  the audit missed something or whether the hypothesis was too optimistic.

The audit reports the honest number, not a contorted one (matches the
"LOC gate: renegotiate to measured" stance in
`2026-05-26-finish-refactor-legacy-sig-removal-design.md`).

Doc-hygiene findings (stale spec/stage comments, comment-as-docstring) do
not count toward the LOC-savings number, but their existence and rough
count is reported in the executive summary as a separate line.

## Out of scope

- Refactoring tests beyond those that block a src simplification
- Modifying `rate_equation` performance characteristics
- Reintroducing parallel front-end representations (the
  `2026-05-26-finish-refactor-legacy-sig-removal-design.md` unification
  is a constraint to preserve, though further simplification within the
  unified family is in scope)
- Pruning the test suite for general redundant coverage (only
  private-helper / structural-key tests are in scope)
- Implementing the simplifications themselves — that happens in
  follow-on impl plans, one per dependency cluster

(Note: replacing the symbol-tuple round-trip in the derivation back-end
IS in scope per Denis's direction, even though prior planning had
deferred it.)

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
  (Pass-1 meta-completeness check passes)
- Every finding cites file:lines and a verifiable evidence chain
- Findings are organized for execution: clusters identified, dependencies
  named, blocking tests flagged, sequencing proposed
- The dual goal is reported honestly: (a) measured non-comment non-doc
  src LOC saving as a number and a percent; (b) a one-paragraph
  description of how the code is simpler after the audit's
  recommendations land
- The half-of-src hypothesis is tested against the honest measured number

The audit does **not** need to find exactly half of src removable. The
honest measured number is the deliverable, whatever it is.
