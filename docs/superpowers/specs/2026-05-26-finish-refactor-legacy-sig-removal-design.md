# Finish the Concrete-Types Refactor — Deletion, Accessor Discipline, Layout

> **Bucket order (locked):** B (dsl dead chain) → B′ (opaque guard) →
> A (legacy Sig) → C (stale memory) → D (accessor discipline) →
> E (layout reorg). Deletions first, accessors next, pure-movement last.

**Date:** 2026-05-26
**Branch:** `refactor-to-concrete-types-instead-of-symbols` (unpushed; main = `7,136` src LOC)
**Status:** design approved by Denis; ready for implementation plan.

## Goal

Close the final success criterion of the Symbol→concrete-struct refactor
(continuation spec §3.1 #1, *"one struct family, no parallel
representations"*) by deleting the now-dead legacy opaque-Symbol code paths.
The refactor is ~90% done: the DSL emits `EnzymeMechanism(Mechanism(...))`,
the enumerator builds decomposed `Mechanism`/`AllostericMechanism` directly,
the form-name string helpers and `_RawSpec` working-rep are gone, and the
five accessors are `@generated` and zero-alloc on the decomposed path. What
remains is genuinely dead code kept alive only by an unused legacy
constructor and an unused legacy parsing branch.

The deletion buckets (A/B/B′) are **pure deletion / collapse** — no new
behavior. Two further maintainability passes follow the deletions: routing
all data access through accessor functions (D) and reorganizing file layout
(E). All three are behavior-preserving: success = same observable behavior,
cleaner structure.

## Precondition (proven 2026-05-26)

The roadmap required proving the bi-bi end-to-end derivation green before
touching legacy code (the legacy Sig path was kept as insurance for opaque
central complexes). **Proven GREEN:** a bi-bi reaction (`A[C], B[N] →
P[C], Q[N]`) produces 77 `init_mechanisms`, and all 77 derive a rate
equation via `compile_mechanism` + `rate_equation_string` with zero "Cycle N
not proportional" failures. This confirms the enumerator routes entirely
through the new decomposed Sig and does not depend on the legacy path.

This probe becomes a committed regression test (see Bucket A, Step 0).

## Scope decisions (locked with Denis)

- **Sequence: Bucket B first, then Bucket A.** B (dsl.jl dead chain) is
  low-risk and independent; land it green as a warm-up before the
  higher-stakes A (legacy Sig).
- **Opaque-rejection guard: REMOVE** (`_assert_no_opaque_terms` /
  `_is_conformation_shape` + tests) — with the caveat below.
- **Out of scope:** the `_n_fit_params_estimate` SS `:NonequalRT` undercount
  (deferred, orthogonal counting bug); contorting code to hit a specific LOC
  number.
- **LOC gate: renegotiate to measured.** Delete everything genuinely dead,
  then record the honest achieved total as the new baseline. LOC is
  non-gating (continuation spec §8). The deletions (A/B/B′) remove ~410 LOC
  (8,781 → ~8,370); the accessor pass (D) is roughly LOC-neutral (a few added
  accessor definitions); the layout pass (E) is net-zero. Record the final
  measured number; do not contort to hit a target.

## Bucket B — dsl.jl dual-grammar dead chain (~80 LOC)

The DSL once emitted a legacy `rxns` tuple-Expr (opaque Symbol shape) *in
parallel* with the decomposed `side_terms_per_step` records. Both callers of
`_parse_steps_block_with_groups` (`src/dsl.jl:618`, `:1204`) discard the
`rxns` slot via `_,`. The entire `rxns` subtree is therefore dead:

| Function | Lines | Role |
|---|---|---|
| `_parse_single_step` | 949–959 | builds one legacy step tuple-Expr; output only ever pushed to discarded `rxns` |
| `_parse_step_side_symbols` | 246–258 | only called by `_parse_single_step` |
| `_step_side_term_to_symbol` | 264–274 | only called by `_parse_step_side_symbols` |
| `_synthesize_species_name` | 389–391 | only called by `_step_side_term_to_symbol` |

Plus: the `rxns = Expr(:tuple)` construction and its `push!(rxns.args, ...)`
sites inside `_parse_steps_block_with_groups` (L828, 847–848, 861–862, 882),
and the unused `.sym` synthesis block inside `_call_form_term_info`
(L369–381, only consumed by the dead `_step_side_term_to_symbol`).

**Change:** `_parse_steps_block_with_groups` returns only
`side_terms_per_step` (no-tag) / `(tags, side_terms_per_step)` (tag); the two
call sites drop the leading `_,`. Delete the four dead functions and the
`.sym` synthesis. Update the docstrings on the touched functions to describe
the decomposed-only behavior (no "legacy"/"mirrors" temporal references).

**Verification:** full suite green; `test_dsl.jl` green; per-commit perf /
compile-budget / chokepoint gates green.

## Bucket B′ — remove opaque-rejection guard (~20 LOC)

Delete `_assert_no_opaque_terms` (L435–451) and `_is_conformation_shape`
(L428–429) and their call sites (L620, L1208), plus the dedicated rejection
tests in `test_dsl.jl` (~L321–342).

**Caveat / hard gate before deleting:** first confirm that opaque forms
(`:ES`, `E_S`) still raise *some* parse error once the explicit guard is
gone (the decomposed grammar only knows call-forms / bare conformations).
Write a throwaway probe: `@enzyme_mechanism` with an opaque step must still
error. If it would instead be **silently accepted** (parsed into a wrong
mechanism), removing the guard is a correctness regression — STOP and keep
the guard instead, returning to Denis. Only proceed with removal if opaque
forms still fail loudly.

**Test integrity:** removing the rejection tests is a coverage reduction —
requires a `docs/superpowers/refactor-deleted-tests.md` §2.1 entry **in the
same commit**, stating that opaque rejection is now enforced by grammar parse
failure rather than an explicit guard (cite the probe result).

## Bucket A — legacy Sig path (~310 LOC, `src/types.jl`)

Reachable from production only via the 1-arg `EnzymeMechanism(Mechanism(...))`
lift, which routes through the new Sig — **not** the legacy path. But the
2-arg `EnzymeMechanism(metabolites, reactions)` constructor **is** called
~16× in `test/test_types.jl` (opaque tuples), where it carries a body of
**validation logic absent from the decomposed `Mechanism` constructor**
(which does only `source_idx` bookkeeping). So Bucket A is not pure
dead-code removal — it must decide the fate of those validators + tests.

**Validator audit (decided: port the meaningful, drop the moot/superseded).**
The decomposed-world invariant home is `_assert_mechanism_invariants`
(`src/mechanism_enumeration.jl:2122`), already called ~80× in tests — *not*
the `Mechanism` constructor, which is deliberately permissive (it is called
constantly by the enumerator; the design intentionally moved structural
invariants out of it — see the `test_types.jl:322–331` note that connectivity
was dropped as a downstream concern).

| Legacy validator | Verdict | Reason |
|---|---|---|
| step is a 4-tuple / `is_eq::Bool` / `gnum::Int` | **drop (moot)** | `Step` is a typed struct — unrepresentable invalid shape |
| exactly one enzyme form per side | **drop (moot)** | `Step.from_species`/`to_species` are single `Species` |
| ≤1 metabolite per side | **drop (moot)** | `Step.bound_metabolite` is a single field |
| empty reactions | **already covered** | `_assert_mechanism_invariants`: `isempty(flat)` |
| every metabolite appears in some step | **port** to `_assert_mechanism_invariants` | representable, not yet checked |
| kinetic-group: RE/SS mixing, same-metabolite, iso-singleton | **port** the not-yet-covered parts to `_assert_mechanism_invariants` | bound/iso consistency already checked there; RE/SS-mix + same-met are not |
| enzyme-form connectivity / orphans | **drop (superseded by design)** | `test_types.jl:322–331` documents the decomposed design intentionally dropped this; Wegscheider handles it downstream |
| stoichiometry rank test | **drop (superseded by design)** | a rank test in the hot constructor is a perf/design regression; not enforced anywhere on the decomposed path today |

Each dropped test gets a `refactor-deleted-tests.md` §2.1 entry citing the
row's reason. Positive construction tests (`n_steps`, `substrates`,
kinetic-group sharing) are **migrated** to decomposed grammar, not deleted.

**Step 0 (first):** commit the bi-bi exit-gate probe as a permanent test
(e.g. in `test/test_mechanism_enumeration.jl` or `test_identify_rate_equation.jl`):
bi-bi `init_mechanisms` → every mechanism `compile_mechanism` +
`rate_equation_string` derives without error. This locks the precondition.

**Deletions:**
1. `_is_new_sig` (L1292–1308) + its doc-comment block (L1280–1287).
2. Collapse the `@generated` accessors that branch on `_is_new_sig(Sig)`:
   drop the `if _is_new_sig(Sig) … end` guard and the trailing legacy body,
   keeping the new-shape body. Sites: `substrates`, `products`, `regulators`,
   `metabolites`, and the others at types.jl ~L1325–1620 (~12 branches).
3. `_mechanism_from_legacy_sig` (L727–~795) + simplify `Mechanism(em)`
   (L711–714) to `_mechanism_from_sig(Sig)`.
4. 2-arg `EnzymeMechanism(metabolites, reactions)` constructor + any
   validators that become dead with it.
5. Dead `_legacy_step_tuple` / `_species_sym` (confirm zero refs first).
6. Rename `_legacy_step_tuple_from_sig` → `_step_tuple_from_sig` (it
   reconstructs the *new* Sig — the "legacy" name is a misnomer; see memory
   `project-ss-dissociation-reconstruction-rule`). Update its callers.
7. Stale doc comments referencing the 2-arg shorthand: types.jl L641, 649,
   796, 1275, 1291; dsl.jl L629; rate_eq_derivation.jl L189.

**Test handling:** the ~16 `test_types.jl` sites using the 2-arg
constructor are resolved per the audit table above — positive construction
tests migrated to decomposed grammar; `@test_throws` validation tests either
re-pointed at `_assert_mechanism_invariants` (ported validators) or deleted
with a §2.1 entry (moot/superseded validators), in the same commit.

**Verification:** full suite green; `check_test_integrity.sh main` EXIT=0;
per-commit perf / compile-budget / chokepoint gates green; `@testset` count
≥ main (minus the documented §2.1 removals).

## Bucket C — housekeeping

- Delete the stale memory `project-accessor-allocates-on-decomposed-sig`
  (resolved by `bc1c592`; `test_accessors.jl` is on decomposed grammar and
  asserts `== 0` allocs).

## Bucket D — accessor discipline (after deletions)

Route all data access through accessor functions instead of direct
`.field` access, for maintainability. Accessors already exist for most
fields (`bound`/`conformation`/`residual` on `Species`;
`from_species`/`to_species`/`bound_metabolite` on `Step`;
`reaction`/`steps`/`regulatory_sites` on `Mechanism`/`AllostericMechanism`).
Direct field access is currently pervasive (~900 sites: `.steps` 216×,
`.reaction` 103×, `.bound` 85×, `.atoms` 83×, `.ligands` 49×, …).

**Scope (locked with Denis): everywhere except definitions.** Convert all
`.field` reads to accessor calls across `src/`, *including* inside
`types.jl` — except the accessor definitions themselves, the struct
constructors, and `@generated`-internal Sig-tuple indexing (those legitimately
touch fields / tuple slots directly). No enforcement test (rely on review).

**No perf risk:** trivial accessors inline to zero cost; the hot
`rate_equation` path uses the separate `@generated` `EnzymeMechanism{Sig}`
accessors (compile-time), not these concrete-struct fields. The
`test_rate_equation_performance` gate is re-verified at the end regardless.

**Step 1 — accessor inventory + add missing.** Enumerate every field of
every concrete struct (`Species`, `Step`, `Mechanism`, `AllostericMechanism`,
`EnzymeReaction`, `ReactantAtoms`, `RegulatorMults`, `RegulatorySite`,
`Residual`, `Metabolite` family). For each field, confirm an accessor exists;
add the missing ones (candidates: `is_equilibrium`/`source_idx` on `Step`;
`cat_allo_states`/`catalytic_multiplicity` on `AllostericMechanism`; `atoms`
on `ReactantAtoms`; `ligands` on `RegulatorySite`/`RegulatorMults`; `name` on
the `Metabolite` family — note `name` is already overloaded for `Parameter`,
so the metabolite method is an added signature). Resolve naming collisions
deliberately: `AllostericMechanism.cat_steps` is exposed as `steps(m)` (no
separate `cat_steps` accessor) — map `.cat_steps` → `steps(m)`.

**Step 2 — mechanical conversion, file by file**, full suite green after
each file. Smallest reasonable diff per commit; commit per file or per type.

## Bucket E — file layout reorganization (last)

Reorder each `src/` file so the reader meets the important things first:
major public types and functions near the top, functions that do similar
things grouped adjacently, less-important helpers at the end.

**Sequenced last** because it's pure movement and is easiest to review
against a stable, already-deleted-and-accessored function set.

**Julia ordering constraints (must be preserved):**
- A `struct` must be defined before any method signature or other struct
  field that references its type. Type/struct blocks stay first in
  `types.jl`.
- A macro must be defined before its use site is parsed. `@enzyme_reaction`
  etc. and their helpers keep their relative order in `dsl.jl`.
- `const`s used at top level must precede their use.
- Plain function *call* resolution is order-independent within the module, so
  helper functions may move freely subject to the above.

**Verification:** this is behavior-preserving movement — full suite green
plus all gates after each file. A failure means something was genuinely
order-dependent (surface it, don't paper over).

## Final verification (continuation spec §3)

- §3.1 #1 one struct family — `grep -rnE "_is_new_sig|_mechanism_from_legacy_sig|EnzymeMechanism\(metabolites" src/` → empty.
- §3.1 #2 no Symbol-string dispatch — `test_chokepoint.jl` green; no
  form-name parse-back helpers (already true).
- §3.1 #3 test integrity — `bash scripts/check_test_integrity.sh main; echo "EXIT=$?"` → 0 (no pipe).
- §3.1 #4 perf — `test_rate_equation_performance` 0-alloc/<100ns green; 3
  compile-budget gates green.
- §3.2 LOC — record measured `wc -l src/*.jl`; renegotiate the target to the
  honest number.
- Full suite green — read the `Test Summary` line, not the notification exit
  code.
- Docs — README + CLAUDE.md re-read for any opaque-form / legacy-Sig
  references invalidated by the deletions; update continuation spec §11
  status to "done".

## Ship

Push the unpushed branch; open the PR with honest LOC numbers (celebrate the
single-struct-family closure; link the deferred parameter-naming refactor as
the next step). PR body ends with the Claude Code attribution footer.

## Conventions (every commit)

- TDD where applicable (this is mostly deletion — the discipline is: delete,
  run the suite, confirm still green). No `--amend`; stay on the branch.
- No deletion/weakening of a test without a `refactor-deleted-tests.md` §2.1
  entry in the same commit.
- No temporal-context comments in code ("legacy", "previously", "Stage N").
- Commit footer: `src delta: -X / +Y net Z, cumulative: ±W` (vs main 7,136);
  end with `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`.

## Self-review notes

- **Riskiest is Bucket A's accessor collapse** — 12 `@generated` branches;
  a wrong collapse silently changes an accessor's output. Mitigated: the
  new-shape body is already exercised by the full suite (DSL emits new Sig),
  so a regression surfaces immediately. Collapse one accessor, run
  `test_accessors.jl`, repeat.
- **Bucket B′ has a real correctness gate** — do not remove the opaque guard
  if opaque forms would become silently accepted.
- **`_assert_mechanism_invariants` is NOT dead** (used ~80× in tests) — do
  not delete it despite its zero `src/` callers.
- **`mechanism_enumeration.jl` is large but not dead** — no Symbol-cleanup
  deletions there; its bulk is the topology backtracker + canonical-hash
  infrastructure. Leave it alone.
- **Don't trust `Pkg.test` notification exit codes** — read the summary line.
