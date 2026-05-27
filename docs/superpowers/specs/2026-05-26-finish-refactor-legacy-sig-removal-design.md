# Finish the Concrete-Types Refactor — Deletion, Accessor Discipline, Layout

> **Bucket order (locked):** B (dsl dead chain) → A (legacy Sig) →
> C (stale memory) → D (accessor discipline) → E (layout reorg).
> Deletions first, accessors next, pure-movement last.
> **Bucket B′ was dropped** — the opaque-rejection guard is load-bearing
> (see "Dropped: Bucket B′" below).

**Date:** 2026-05-26
**Branch:** `refactor-to-concrete-types-instead-of-symbols` (unpushed; main = `7,136` src LOC)
**Status:** design approved by Denis; revised post-review 2026-05-26.

## Goal (honest scope)

Delete the now-dead **legacy** opaque-Symbol code paths (the `(mets, rxns)`
Sig shape + its constructor and accessor branches), route field access
through accessors, and reorganize file layout. All behavior-preserving.

**This does NOT achieve "one struct family, no parallel representations"
(continuation spec §3.1 #1).** A 2026-05-26 adversarial review established
that the opaque Symbol representation survives as the **derivation
back-end's working format**: `_species_name_from_sig` (types.jl:1408) still
synthesizes form-name Symbols by string-join (`Symbol(join(parts, "_"))`),
`_step_tuple_from_sig` rebuilds opaque `(lhs, rhs, is_eq, g)` tuples at
`@generated` time, and the King-Altman/Wegscheider derivation
(`rate_eq_derivation.jl:142–143,445`,
`thermodynamic_constr_for_rate_eq_derivation.jl:88–90,207–209`) consumes
`reactions(m)`/`enzyme_forms(m)` and matches enzyme forms by Symbol
identity. The Sig↔opaque-tuple round-trip is the single largest remaining
parallel representation, and it is **out of scope here** — it is the
genuine next phase (see "Deferred — derivation back-end unification").

So this spec's accurate claim is narrower: it unifies the **front end**
(DSL + enumeration + the public `Mechanism`/`AllostericMechanism` surface)
onto one struct family and removes the *legacy* Sig. The PR must say this
plainly rather than claiming full closure of §3.1 #1.

The deletion buckets (A/B) are **pure deletion / collapse**. Two further
maintainability passes follow: routing data access through accessors (D)
and reorganizing file layout (E). All are behavior-preserving.

## Dropped: Bucket B′ (opaque-rejection guard)

The earlier plan removed `_assert_no_opaque_terms` / `_is_conformation_shape`.
**Review proved this is a correctness regression and the bucket is dropped.**
Without the guard, an opaque step like `ES <--> E + P` is **silently
accepted** — `ES` parses as a bare conformation via `_term_bare_enzyme`
(`dsl.jl:419`), the product is dropped, and a structurally wrong mechanism is
built (the new Sig is lossy for opaque Species). `_is_conformation_shape`'s
regex rejects `ES` (uppercase `S`); the guard is the *only* thing turning
that into an error. **Keep the guard and its `test_dsl.jl` tests.**

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
  non-gating (continuation spec §8). The deletions (A/B) remove ~390 LOC
  (8,781 → ~8,390); the accessor pass (D) is roughly LOC-neutral (a few added
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
and the dead `:call`-branch `.sym` *synthesis* inside `_call_form_term_info`
(L369–378, the `parts`/`join` Symbol-building consumed only by the dead
`_synthesize_species_name`).

**`_StepSideTerm.sym` field stays — do NOT drop it.** Review verified the
field is read on the live emission path (`_build_step_expr` →
`bound_met_term.sym` at `dsl.jl:705`; `_split_side` at 720/726/734) and by
the opaque guard at 445. For `:metabolite`/`:bare_enzyme` terms `sym` is the
real name. Only the `:call`-branch synthesis is dead; give the `:call`
`_StepSideTerm` a non-synthesized first arg (the bare conformation name).

**Change:** `_parse_steps_block_with_groups` returns only
`side_terms_per_step` (no-tag) / `(tags, side_terms_per_step)` (tag); the two
call sites drop the leading `_,`. Delete the four dead functions and the
`:call` `.sym` synthesis. Update the docstrings on the touched functions to
describe the decomposed-only behavior (no "legacy"/"mirrors" references).

**Verification:** full suite green; `test_dsl.jl` green; per-commit perf /
compile-budget / chokepoint gates green.

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
| every **substrate/product** appears in some step | **port** to `_assert_mechanism_invariants` | representable, not yet checked. **Regulators are EXCLUDED** — `_drop_unbound_regulators` (types.jl:671) intentionally lets `init_mechanisms` declare an inhibitor no step binds (bound later by `expand_mechanisms`); checking regulators would error on every dead-end-inhibitor init mechanism (review-confirmed). Mirror `_drop_unbound_regulators`: substrates/products only. |
| kinetic-group: RE/SS mixing, same-metabolite, iso-singleton | **port** the not-yet-covered parts to `_assert_mechanism_invariants` | bound/iso consistency already checked there; RE/SS-mix + same-met are not |
| enzyme-form connectivity / orphans | **drop (superseded by design)** | `test_types.jl:322–331` documents the decomposed design intentionally dropped this; Wegscheider handles it downstream |
| stoichiometry rank test | **drop (superseded by design)** | a rank test in the hot constructor is a perf/design regression; not enforced anywhere on the decomposed path today |

Each dropped test gets a `refactor-deleted-tests.md` §2.1 entry citing the
row's reason. Positive construction tests (`n_steps`, `substrates`,
kinetic-group sharing) are **migrated** to decomposed grammar, not deleted.

**The ported invariant must not reject any currently-valid enumerated
mechanism.** After porting, run the full `test_mechanism_enumeration.jl` —
all ~80 existing `_assert_mechanism_invariants` calls must stay green. If one
errors, the ported check is wrong for the decomposed world — STOP and
re-examine (do not loosen the test).

**Ordering within Bucket A (review-corrected — the original order produced a
red suite mid-bucket).** The accessor collapse must come *after* all
legacy-Sig construction is gone, because `test_types.jl` calls accessors on
legacy-Sig mechanisms. Order:
1. **Exit-gate test** (subset — see below).
2. **Port validators** into `_assert_mechanism_invariants`.
3. **Migrate ALL legacy-Sig test construction** in `test_types.jl` to
   decomposed grammar + re-point/drop the validator `@test_throws`.
4. **Delete the 2-arg constructor + its now-orphaned validators.**
5. **Collapse the ~12 accessor branches** (safe now: no legacy Sig exists).
6. **Delete `_is_new_sig` + `_mechanism_from_legacy_sig`** (now unused) +
   fix stale doc comments.

**Exit-gate test (cost-bounded).** A test deriving `rate_equation_string`
for *all 77* bi-bi `init_mechanisms` costs ~86 s cold (review-measured) and
risks the compile-budget gate. So: assert `length(init_mechanisms(bi_bi)) ==
77` (cheap) AND derive only a small fixed subset (e.g. the first 5 by form
count). The precondition is "enumerator routes through the new Sig" — a
subset proves it; full-77 derivation is not worth ~80 s on every cold run.

**Deletions (line numbers will drift — re-grep before editing):**
- `_is_new_sig` (~L1292–1308) + doc block (~L1280–1287).
- The ~12 `@generated` accessor `if _is_new_sig(Sig) … end` guards +
  trailing legacy bodies (types.jl ~L1325–1620): collapse to the new body.
- `_mechanism_from_legacy_sig` (~L727–795) + simplify `Mechanism(em)` to
  `_mechanism_from_sig(Sig)`.
- 2-arg `EnzymeMechanism(metabolites, reactions)` (~L803–888) + the
  validators that become dead with it (`_validate_kinetic_groups`,
  `_validate_enzyme_connectivity`, `_validate_stoichiometry`,
  `_pretty_reaction`).
- Stale doc comments referencing the 2-arg shorthand (types.jl ~L641, 649,
  796, 1275, 1291; rate_eq_derivation.jl ~L189).
- *Already done:* `_legacy_step_tuple_from_sig` → `_step_tuple_from_sig`
  rename (commit 388944c); no `_legacy_*` names remain. Just confirm.

**Test handling — full scope (review-corrected).** `test_types.jl`
constructs legacy-Sig mechanisms BOTH via the 2-arg constructor AND directly
via `EnzymeMechanism{(mets, rxns)}()` (lines 9, 33, 67, 413). The migration
must cover **all** of these (grep both `EnzymeMechanism(((` and
`EnzymeMechanism{(`), including the "EnzymeMechanism struct + accessors"
testset (~2–37), the "stoich_matrix" testset (~60–74), and the
`AllostericEnzymeMechanism` validator testset whose `cm_bad`/`cm` fixtures
are legacy-Sig (~372–417, 449). Per the audit table: migrate positive tests
to decomposed grammar; re-point ported-validator `@test_throws` at
`_assert_mechanism_invariants`; delete moot/superseded `@test_throws` with
§2.1 entries.

**Integrity caveat (review-flagged).** `check_test_integrity.sh` counts
`@testset` headings, not `@test_throws` *inside* a surviving testset — so
deleting individual assertions inside a kept testset passes the gate
silently. EXIT=0 is necessary but NOT sufficient here: **manually confirm a
§2.1 entry exists for each dropped `@test_throws`**, and that each is either
moot or re-pointed.

**Verification:** full suite green; `check_test_integrity.sh main` EXIT=0
(+ the manual assertion-level check above); per-commit perf / compile-budget
/ chokepoint gates green.

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
every concrete struct. Review verified **most accessors already exist** —
`name` (Metabolite family, types.jl:23–26), `is_equilibrium` (184), `ligands`
(120), `atoms` (258), `catalytic_multiplicity` (469), plus
`bound`/`conformation`/`residual`/`from_species`/`to_species`/`bound_metabolite`/`reaction`/`steps`/`regulatory_sites`.
**Genuinely missing (add these):** `source_idx` (Step), `cat_allo_states`
(AllostericMechanism), and `allo_states` on `RegulatorySite` (read 6× as
`site.allo_states` in `rate_eq_derivation.jl`). Resolve naming collisions
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

## Final verification

- **Legacy Sig gone** — `grep -rnE "_is_new_sig|_mechanism_from_legacy_sig|EnzymeMechanism\(metabolites" src/` → empty.
- **§3.1 #2 (partial — front end only)** — `test_chokepoint.jl` green. **Note
  honestly:** `_species_name_from_sig`/`_step_tuple_from_sig` still synthesize
  opaque form-names + tuples for the derivation back-end; this criterion is
  NOT fully closed (see Deferred). Do not claim otherwise.
- **§3.1 #3 test integrity** — `bash scripts/check_test_integrity.sh main; echo "EXIT=$?"` → 0 (no pipe) + the manual assertion-level check (Bucket A).
- **§3.1 #4 perf** — `test_rate_equation_performance` 0-alloc/<100ns green; 3
  compile-budget gates green (re-measure after the exit-gate test lands).
- **LOC** — record measured `wc -l src/*.jl`; renegotiate the target to the
  honest number (note: branch is +1,645 over main *before* deletions; even
  after, it lands well above main — the cost of structured types).
- **Full suite green** — read the `Test Summary` line, not the exit code.
- **Docs** — README + CLAUDE.md re-read for opaque-form / legacy-Sig
  references invalidated by the deletions.

## Deferred — derivation back-end unification (the genuine next phase)

The headline "one struct family, no parallel representations" is **not**
achieved by this spec. The opaque Symbol representation survives as the
derivation back-end's working format: `_species_name_from_sig` synthesizes
form-name Symbols by string-join, `_step_tuple_from_sig` rebuilds opaque
`(lhs, rhs, is_eq, g)` tuples at `@generated` time, and King-Altman/Wegscheider
(`rate_eq_derivation.jl`, `thermodynamic_constr_for_rate_eq_derivation.jl`)
run on those tuples via `reactions(m)`/`enzyme_forms(m)`, matching forms by
Symbol identity. Making the `@generated` derivation consume `Step`/`Species`
structurally is the real remaining unification — **high-risk** (the
non-negotiable `rate_equation` 0-alloc/<100ns gate constrains any change to
the `@generated` path) and large. **It gets its own brainstorm → spec →
plan**, not this PR.

## Ship

Push the unpushed branch; open the PR. **Honest framing:** this removes the
*legacy* Sig and unifies the front end onto one struct family; it does NOT
close §3.1 #1 — the Sig↔opaque-tuple derivation back-end remains a deliberate
second representation, scoped as the next phase. Honest LOC numbers. Link the
deferred derivation-back-end unification and the parameter-naming refactor as
next steps. PR body ends with the Claude Code attribution footer.

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
