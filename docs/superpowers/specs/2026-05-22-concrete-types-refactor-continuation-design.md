# Concrete-types refactor — continuation design (Stages 6β through 7e)

**Branch:** `refactor-to-concrete-types-instead-of-symbols`
**Status:** spec; supersedes the remaining-stages portion of
[`2026-05-20-concrete-types-refactor-design.md`](2026-05-20-concrete-types-refactor-design.md)
**Date:** 2026-05-22
**Prerequisite reading:** the original design doc above, plus
`.claude/CLAUDE.md` (project rules).

## 1. Why a continuation spec

The original 2026-05-20 spec has been partly executed (Stages 0–6α) and
partly deviated from. The branch is at 9,743 src LOC (+2,607 vs main).
The spec's headline success metric (§3: ≤3,600 src LOC) is unreachable
without removing the parameter-indexing scheme — a change that's been
deliberately deferred to a future refactor.

This document re-anchors the remaining work with:

- Revised success criteria that reflect the actual architectural cost.
- An explicit staging plan for Stages 6β through 7e, including the
  three new stages (7a–7d) that absorb the work the original spec
  deferred under "approach C" and "option D".
- The list of legacy infrastructure the continuation will delete.
- Explicit declaration of what is **out of scope**, so the
  follow-up refactor's surface is clear.

Everything not contradicted here — test-integrity rule, perf gates,
the chokepoint architecture, the no-parallel-representations principle
— carries forward unchanged from the original spec.

## 2. Status snapshot (2026-05-22)

- Branch: `refactor-to-concrete-types-instead-of-symbols`, unpushed.
- Tip: `895447c` (Stage 6α.3 complete).
- Tags: `stage-1-complete`, `stage-2-complete`, `stage-3-complete`,
  `stage-4-complete-checkpoint`.
- Tests: 26,904 pass; `test_rate_equation_performance` 0-alloc/<100ns
  green; all 3 compile-budget gates green.
- A previous Stage 6β.1 attempt deleted ~2,940 test lines without
  preserving coverage. Denis rejected it; the branch is now at the
  pre-deletion state. The mistake informs §4 below.

## 3. Revised success criteria

The original spec ranked LOC reduction as the primary success metric.
Reality has shifted the weighting:

### 3.1 Primary criteria (gating)

1. **One concrete struct family.** No parallel representations of the
   same concept remain. After Stage 7d, the codebase has exactly one
   way to express a Step, a Species, a Mechanism, a Reaction.
2. **No Symbol-string dispatch.** Parameter names are rendered through
   the chokepoint `name(p::Parameter, m::Mechanism)`. Form names are
   not parsed back into structure. The seven remaining direct
   `Symbol("K$idx")` constructions are routed through the chokepoint
   (Stage 7a).
3. **Test integrity preserved.** Per original §2 — no test deleted,
   commented, weakened, or `@test_skip`/`@test_broken` tagged.
   `bash scripts/check_test_integrity.sh main` PASSES at the final
   commit.
4. **Performance gates green.** `rate_equation` 0-alloc/<100ns per
   call invariant holds. Compile-budget gates remain at 2× main.

### 3.2 Secondary criterion (tracked, not gating)

5. **LOC reduction.** Target: **≤ 8,200 src LOC** at the final commit
   (≈ +15% vs main, ≈ −16% vs current). Tracked per-commit with the
   `src delta: -X / +Y net Z, cumulative: ±W` line.

   The original ≤3,600 target was based on an architecture vision
   that didn't survive contact with Julia type-parameter constraints
   (the `{Sig}` collapse deviation) and with the legacy-lift
   compromise (opaque-form Species). Reaching ≤3,600 requires
   removing parameter-name indexing, which is deferred (see §7).

### 3.3 Foundation for the deferred refactor

A future refactor will replace positional parameter names (`:K1`,
`:k10f`) with semantic names (`:K_ATP`). To keep that refactor
contained, this continuation establishes the foundation:

- `Step.source_idx` is the only field carrying the index; it is
  presentation metadata (Step `==`/`hash` ignore it) — already in
  place.
- All Symbol rendering of parameter names flows through the
  chokepoint `name(p::Parameter, m::Mechanism)` — Stage 7a closes
  the remaining gaps.
- A regression test asserts the chokepoint exclusivity.

After this continuation lands, the parameter-naming refactor becomes
"change the body of `name(p, m)` plus a single sweep of test fixtures
that hardcode `:K1`-style symbols."

## 4. Test-migration rule (clarifies original §2 in light of Stage 6β.1)

The original §2 prohibits test deletion. The Stage 6β.1 mistake was
deleting spec-overload testsets *before* verifying that equivalent
Mechanism-overload testsets covered the same surface.

**Rule for the continuation:** when retiring a tested code path, every
deleted testset must be replaced by a Mechanism-form testset *in the
same commit*, with assertions that cover the same surface (same input
shapes, same expected counts, same property checks). The deletion
commit must include — in the commit body or the test diff —
testset-by-testset evidence that the new tests subsume the old ones.

This is a tightening of §2.1's "narrow exception" log requirement: the
log is still required for any deletion, AND the replacement must be in
the same commit (not "we'll add it later"). If you can't fit both in
one commit, split the work: migration commit first (adds Mechanism
tests), deletion commit second (removes spec tests).

## 5. Continuation stages

Each stage's commit messages end with the LOC line per original spec
§3:
```
src delta: -X / +Y net Z, cumulative: ±W
```

### Stage 6β — Migrate spec testsets to Mechanism (opaque-form)

**Goal:** retire `MechanismSpec`/`StepSpec`/`AllostericMechanismSpec`
from `test_mechanism_enumeration.jl` and from `src/` without losing
assertion coverage.

**Approach:** translate each `MechanismSpec(rxn, [StepSpec(...)], n)`
construction to `Mechanism(rxn, [[Step(Species([], :E_S),
Species([], :E), Substrate(:S), true)], ...])`. Use **opaque-form
Species** (`Species([], :ES)` — empty bound list, opaque name) so
today's `_expand_*` moves work unchanged. Decomposed Species lands in
Stage 7b.

**Commits (suggested, all on `refactor-to-concrete-types-instead-of-symbols`):**

1. Migrate test helpers (`enumerate_all`, `_assert_*_invariants`,
   `mechanism_spec_from_mechanism_and_rxn`).
2. Migrate `_expand_re_to_ss` testsets.
3. Migrate `_expand_split_kinetic_group` testsets.
4. Migrate `_expand_add_dead_end_regulator` testsets.
5. Migrate `_expand_to_allosteric` testsets.
6. Migrate `_expand_add_allosteric_regulator` testsets.
7. Migrate `_expand_change_allo_state` testsets.
8. Migrate `_canonicalize!`/`_dedup_key`/dedup-integration testsets.
9. Migrate `expand_mechanisms` and `Integration` testsets.
10. Delete spec types (`MechanismSpec`, `StepSpec`,
    `AllostericMechanismSpec`), their `_expand_*` overloads,
    `dedup!(::Dict{Int,Vector{<:AbstractMechanismSpec}})`,
    `_canonicalize!`, `_dedup_key`, `_init_mechanism_specs`,
    `_spec_from_mechanism`, `_mechanism_from_spec`,
    `EnzymeMechanism(::MechanismSpec)` adapter,
    `expand_mechanisms(::Vector{<:AbstractMechanismSpec}, ...)` overload,
    `_n_fit_params_estimate(::AbstractMechanismSpec)` if separate.

**LOC delta target:** Stage 6β cumulative ≈ −500 to −800 in
`src/mechanism_enumeration.jl`.

**Symbol form-name helpers stay** (`_form_name`, `_parse_bound`,
`_bound_mets_from_form_name`, `_dead_end_form_name`, `_atoms_dict`,
`_is_estar_form`, `_can_pingpong`, `_subtract_atoms`). They're called
by the native Mechanism `_expand_*` moves to interpret opaque-form
Species. They get deleted in Stage 7b.

### Stage 6 — `identify_rate_equation` canonicalizer rewrite

**Goal:** as originally planned (plan §6 / lines 3027-3450 of
`docs/superpowers/plans/2026-05-20-concrete-types-refactor.md`).
Replace the regex-pipeline canonical hash with a structural hash over
the `Mechanism` form.

**LOC delta target:** −300 to −500 in `src/identify_rate_equation.jl`.

**No revision to the original plan section.** Execute as written.

### Stage 7a — Chokepoint hygiene

**Goal:** route every parameter-name rendering through
`name(p::Parameter, m::Mechanism)`. After this stage, the chokepoint
is the sole `Parameter → Symbol` mapping in `src/`.

**Inventory of sites bypassing the chokepoint today:**
- `src/rate_eq_derivation.jl` lines 132, 148, 627, 1632 — four
  `Symbol("K$idx")` constructions.
- `src/types.jl` lines 1840, 1841 — `Symbol("K_$(lig_name)_T_reg$site_idx")`
  and `Symbol("K_$(lig_name)_reg$site_idx")`.

**Approach:** at each site, identify the Parameter struct value that
should produce the Symbol (a `K`, `Kr`, `Kreg`, etc., with the right
state and site fields), call `name(p, m)` instead. If the Parameter
isn't readily available at the call site, build it from the Step and
Mechanism in scope.

**Regression test:** add a test that scans `src/*.jl` for direct
`Symbol("K…")`/`Symbol("k…")`/`Symbol("L")`/`Symbol("V…")`
construction outside of `name(p, m)` implementations and fails if any
slip back in. Likely a grep-based test in `test/test_aqua_jet.jl` or a
new `test/test_chokepoint.jl`.

**LOC delta target:** ≈ 0 net (substitution, not deletion). Value is
architectural: unblocks the future parameter-naming refactor.

### Stage 7b — DSL native emission with decomposed Species

**Goal:** macros `@enzyme_mechanism` and `@allosteric_mechanism` emit
`EnzymeMechanism(Mechanism(...))` with **decomposed Species** —
i.e., `Species([Substrate(:S)], :E)` instead of opaque
`Species([], :E_S)`. After this stage, the legacy 2-tuple form
`(metabolites_tuple, reactions_tuple)` and the opaque-form helpers
are deletable.

**DSL grammar change (per original spec §2.2 "approach B"):**

Today's grammar uses bare metabolite-bound form names:
```julia
steps: begin
    E + S => ES
    ES => EP
    EP => E + P
end
```

New grammar uses function-call decomposed-Species notation:
```julia
steps: begin
    E + S => E(S)
    E(S) => E(P)
    E(P) => E + P
end
```

`E(S)` parses to `Species([Substrate(:S)], :E)`. `E(S, ATP)` parses
to `Species([Substrate(:S), Substrate(:ATP)], :E)`. Iso conformations
keep their bare name (`Estar`).

**Approach (dual-grammar transition):**

The DSL parser accepts both grammars during migration, distinguished
by syntax (`ES` is a `Symbol`; `E(S)` is a `Call` expression). Old
fixtures keep working while new fixtures roll in file-by-file. After
all fixtures have migrated, the old-grammar branch is deleted in a
single cleanup commit.

1. Extend the DSL parser to recognise `E(S)` / `E(S, ATP)` /
   `Estar` (bare-name iso) as decomposed-Species notation, alongside
   the existing bare-form path. Emit `Mechanism(...)` with
   `Species([Substrate(:S)], :E)` for decomposed entries.
2. Switch the macros to emit `EnzymeMechanism(Mechanism(...))`
   directly when all step entries are decomposed; fall back to the
   legacy `EnzymeMechanism(metabolites, reactions)` shape when any
   step is bare-form.
3. Migrate `test/mechanism_definitions_for_test_enzyme_derivation.jl`
   (37 DSL invocations, the bulk of the analytical fixture
   formulas). Analytical formulas referencing `:K1`-style param names
   stay unchanged — they belong to the deferred refactor (§7).
4. Migrate the remaining DSL-using test files: `test/test_dsl.jl`
   (36), `test/test_types.jl` (15), `test/test_rate_eq_derivation.jl`
   (14), `test/test_accessors.jl` (3), `test/test_identify_rate_equation.jl`
   (3), `test/test_fitting.jl` (2), `test/test_compile_budget.jl` (1).
5. Rewrite the native `_expand_*` moves in
   `src/mechanism_enumeration.jl` to read decomposed Species fields
   directly. Safe at this point because every Species reaching the
   moves comes through new-grammar DSL or through Stage 6β-migrated
   construction sites (which can themselves migrate to decomposed
   form once the moves accept it).
6. Delete the old-grammar parser branch and the legacy `(metabolites,
   reactions)` macro emission path.
7. Delete opaque-form helpers (`_form_name`, `_parse_bound`,
   `_bound_mets_from_form_name`, `_dead_end_form_name`, `_atoms_dict`,
   `_is_estar_form`, `_can_pingpong`, `_subtract_atoms`) and
   `_mechanism_from_legacy_sig` plus the `_is_new_sig` branch in
   `Mechanism(em::EnzymeMechanism{Sig})`.

**LOC delta target:** net ≈ −200 to −250 (DSL rewrite ≈ neutral;
opaque-form helper deletion ≈ −200; `_mechanism_from_legacy_sig`
deletion ≈ −60; fixture rewrites are tests, not src).

### Stage 7c — Enumeration dispatches on `EnzymeReaction`

**Goal:** the ~20 sites in `src/mechanism_enumeration.jl` that
currently dispatch on `EnzymeReactionLegacy` move to dispatch on
`EnzymeReaction` (the concrete struct). After this stage, no internal
call path needs `_to_legacy_reaction`.

**Affected dispatch points (inventory):**
- `src/mechanism_enumeration.jl` ~20 method signatures
- `src/identify_rate_equation.jl`:
  - `IdentifyRateEquationProblem{R<:EnzymeReactionLegacy, D<:NamedTuple}`
    → `{R<:EnzymeReaction, D<:NamedTuple}`
  - `identify_rate_equation(reaction::EnzymeReactionLegacy, ...)` →
    `(reaction::EnzymeReaction, ...)`

**Compile-budget consideration:** `test_compile_budget.jl` (lines 77,
147, etc.) currently constructs `EnzymeReactionLegacy` directly to
exercise per-arity specialization paths. After 7c, those tests
construct `EnzymeReaction` instead. The compile-budget gate remains at
2× main. Since main's measurements used `EnzymeReactionLegacy`, the
budget values may need re-baselining at this stage. If the new
non-singleton `EnzymeReaction` produces *better* compile metrics
(fewer specialized methods), reduce the budget accordingly so the gate
stays meaningful. **Do not raise the budget** — spec §3 carry-over.

**LOC delta target:** small net (dispatch swap; deletions land in 7d).

### Stage 7d — Delete legacy infrastructure

**Goal:** delete all the now-unreachable legacy code paths. After
this stage, `_to_legacy_reaction` has no callers and goes; the
`EnzymeReactionLegacy` struct has no users and goes; the dual-Sig
accessor branches in `src/types.jl` collapse to single bodies.

**Deletion inventory:**

| Item | File | Approx LOC |
|------|------|-----------|
| `EnzymeReactionLegacy{S,P,R,N}` struct + outer constructor | `types.jl` 628–693 | ~70 |
| `EnzymeReactionLegacy` accessors (`substrates`, `products`, `regulators`, `regulator_roles`, `oligomeric_state`, `Base.show`) | `types.jl` 1451–1488, 1221 | ~80 |
| `_to_legacy_reaction(r::EnzymeReaction)` | `types.jl` 697–730 | ~35 |
| 2-arg `EnzymeMechanism(metabolites, reactions)` constructor + atom-balance check | `types.jl` 749–870 | ~120 |
| Dual-Sig branches in 12 `EnzymeMechanism` accessors (`_is_new_sig` guards) | `types.jl` (scattered) | ~150 |
| `_n_fit_params_estimate(::AbstractMechanismSpec)` (if not already gone in 6β) | `mechanism_enumeration.jl` | ~10 |
| **Total** | | **~455** |

**LOC delta target:** −450 to −550. Stage 7d is the second-largest
deletion stage of the continuation (after Stage 6β's spec-type
removal).

### Stage 7e — Cleanup + docs + PR

**Goal:** as in original spec §7 (Tasks 7.1–7.5). Dead-code sweep,
README updates, CLAUDE.md updates, final PR.

**Adjustments from the original:**
- Update README's DSL examples to the new decomposed-Species grammar.
- Update CLAUDE.md "Source Layout" and "Key Architecture Decisions" to
  match the final architecture: single `Mechanism`/`Step`/`Species`
  family, single `EnzymeReaction`, chokepoint `name(p, m)`, no
  opaque-form helpers, no spec types.
- The PR description ships honest LOC numbers — celebrate the
  architectural simplification, acknowledge the +15% LOC cost as the
  price of structured types, and link the future parameter-naming
  refactor as the next step.

**Final verification (per original §7.5):**
- All tests green
- `test_rate_equation_performance` 0-alloc/<100ns gate green
- All compile-budget gates green (possibly re-baselined in 7c)
- `bash scripts/check_test_integrity.sh main` PASSES
- @testset count grew vs main (new tests added in 1.A, 6β, 7a)
- `wc -l src/*.jl` ≤ 8,200

## 6. Stage-by-stage LOC budget summary

| Stage | Target Δ | Cumulative target |
|-------|----------|-------------------|
| Start | — | +2,607 |
| 6β | −500 to −800 | ≈ +1,900 |
| 6 | −300 to −500 | ≈ +1,500 |
| 7a | ≈ 0 | ≈ +1,500 |
| 7b | −200 to −250 | ≈ +1,300 |
| 7c | small | ≈ +1,300 |
| 7d | −450 to −550 | ≈ +800 |
| 7e | −100 | ≈ +700 |
| **End** | | **≈ 7,800 (range 7,500–8,200)** |

If any stage exceeds budget, the next stage absorbs the slack. The
≤8,200 final target is the gate; intermediate per-stage targets are
guides.

## 7. Out of scope (deferred to a follow-up refactor)

**Parameter naming change** (`:K1` → `:K_ATP` style).

This refactor is deliberately left out for the reasons documented in
§3.3. Its scope is roughly:

- Change `name(p::K, m)`, `name(p::Kreg, m)`, and the other Parameter
  `name` methods to emit semantic names.
- Remove `Step.source_idx` field and `_rep_idx_for_step` helper.
- Rewrite analytical formulas in
  `test/mechanism_definitions_for_test_enzyme_derivation.jl` (~2,306
  LOC, much of it formula text) to use new param names.
- Update ~57 hand-written `:K1`-style assertions across test files.
- Update external-API docstrings explaining parameter names.

After this continuation lands, that refactor's surface is contained:
all parameter-name rendering already flows through `name(p, m)`
(Stage 7a), the `source_idx` field is already isolated, and the
mechanical fixture rewrites are the main cost.

## 8. Non-negotiables (carry-over from original spec)

The following clauses of the 2026-05-20 spec are unchanged and apply
to every commit of the continuation:

- **§2 test integrity** — no test deletion, weakening,
  `@test_skip`/`@test_broken`, or commenting out. Mechanical syntax
  adaptation only. §4 of this document tightens the migration rule
  further.
- **§2.1 deleted-tests log** — if any deletion fits the narrow
  exception (helper deletion with no surviving referent), append
  to `docs/superpowers/refactor-deleted-tests.md` (create the file
  on first use). Plus the §4 same-commit-replacement requirement.
- **Performance gates** — `rate_equation` 0-alloc/<100ns invariant;
  compile-budget gate at 2× main; `loss!` runtime no regression vs
  main baseline.
- **Test-integrity script** — `bash scripts/check_test_integrity.sh
  main` PASSES at every commit.
- **No `--amend`** — always create new commits.
- **No temporal-context comments** — no "Stage N", "previously",
  "legacy", "will be" in code comments. Documentation in this spec
  and the plan is fine; code comments must be evergreen.

## 9. Existing-deviation reaffirmation

The following deviations from the original spec are now formalized as
the continuation's working assumptions (do not roll them back):

1. **`Step.source_idx::Int` field exists.** Step `==`/`hash` ignore it.
   It's presentation metadata for parameter naming.
2. **`AllostericEnzymeMechanism{Sig}` collapse not attempted.** Stays
   as `{CM, CS, RS}` (3-param). Julia rejects DataType inside
   value-tuple type parameters; the converter
   `AllostericMechanism(::AllostericEnzymeMechanism)` carries the
   load.
3. **Iso-step canonicalization in `Step` constructor** —
   deterministic lex-order direction. Causes the Segel Iso Uni Uni
   fixture's `k4f`/`k4r` confusion under opaque-form lift; resolved
   when fixtures move to decomposed Species in Stage 7b.
4. **`regulators:` DSL alias dropped** — 8 test sites updated to
   `competitive_inhibitors:`.
5. **Memory note correction** —
   `project_dedup_pass2_dead_code` was wrong; pass-2 Wegscheider
   absorption IS load-bearing for inhibitor/activator mechanisms.
   `_build_kinetic_rename_map` is kept.
