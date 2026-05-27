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

**Goal:** retire the *test-side* construction of
`MechanismSpec`/`StepSpec`/`AllostericMechanismSpec` from
`test_mechanism_enumeration.jl` without losing assertion coverage,
plus delete the spec types' public surface and accessor helpers in
`src/mechanism_enumeration.jl` that have no remaining `src/` callers
after the migration. The **heavy internal enumeration pipeline** —
`init_mechanisms(::EnzymeReactionLegacy)`, `_apply_equivalence_grouping`,
`_expand_substrate_product_dead_ends`, and helpers — still
**constructs `MechanismSpec` values internally** because the algorithm
operates on the spec form; that rewrite is Stage 7d.

After 6β: the spec types remain alive in `src/` but only as internal
implementation detail of the enumeration; no test or external caller
constructs them; the `Mechanism`-returning entry points
(`init_mechanisms(::EnzymeReaction) → Vector{Mechanism}`) are the sole
public API.

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
10. Delete only the spec types' **test-facing surface** in
    `src/mechanism_enumeration.jl`: the six `_expand_*`
    `AbstractMechanismSpec` overloads;
    `dedup!(::Dict{Int,Vector{<:AbstractMechanismSpec}})`;
    `_canonicalize!(::AbstractMechanismSpec)`;
    `_dedup_key(::AbstractMechanismSpec)`;
    `_n_fit_params_estimate(::AbstractMechanismSpec)` if a separate
    method from the Mechanism version.
    For each, verify zero src/ callers via `grep -rn` BEFORE deleting;
    if any have callers, defer to Stage 7d.

    **Keep alive (deleted in Stage 7d):** `MechanismSpec`, `StepSpec`,
    `AllostericMechanismSpec`, `AbstractMechanismSpec` structs;
    `_init_mechanism_specs(::EnzymeReaction)` AND
    `_init_mechanism_specs(::EnzymeReactionLegacy)` (both still feed
    the live `init_mechanisms(::EnzymeReaction)` adapter at
    `mechanism_enumeration.jl:3390`); `_spec_from_mechanism`,
    `_mechanism_from_spec`; `EnzymeMechanism(spec::MechanismSpec)`
    adapter; `AllostericEnzymeMechanism(spec::AllostericMechanismSpec)`;
    `_apply_equivalence_grouping`;
    `_expand_substrate_product_dead_ends` if it still returns specs;
    `_catalytic_topologies` and `_atoms_dict` (called by the heavy
    pipeline on the Legacy path).

**LOC delta target:** Stage 6β cumulative ≈ −300 to −500 in
`src/mechanism_enumeration.jl` (test-surface deletions only;
heavy-pipeline deletions land in 7d).

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

**Goal:** route every parameter-name rendering through the
`name(...)` chokepoint, in two complementary forms:
- **Value-context** (existing): `name(p::Parameter, m::Mechanism)`
  takes a fully constructed Parameter (which carries a `step::Step`)
  and a Mechanism.
- **Type/index-context** (new — added by this stage):
  `name(::Type{P}, idx::Int) where P<:Parameter` takes a Parameter
  type and a kinetic-group rep-index, returning the same Symbol the
  value-context form would produce for that group.

The two methods delegate to a shared private formatter. The reason
for splitting: rate-equation derivation in `src/rate_eq_derivation.jl`
runs inside `@generated` functions that walk the `EnzymeMechanism{Sig}`
type at compile time — at those sites only an `Int` group index is in
scope; constructing a `Step` value to satisfy the value-context
signature would require walking the Sig manually, defeating the
chokepoint. The type/index variant fits the @generated call shape.

After this stage, the chokepoint (these two methods) is the sole
`Parameter → Symbol` mapping in `src/`. The future parameter-naming
refactor changes both bodies in lockstep — still a contained edit.

**Inventory of sites bypassing the chokepoint today:**

`src/rate_eq_derivation.jl` — 10 sites (in @generated context, no Step
in scope; use the type/index variant):
- Lines ~132, 148, 627, 1632 — `Symbol("K$idx")` (Kd rendering, R-state)
- Lines ~134, 135, 631, 632, 1636, 1637 — `Symbol("k$(idx)f")` /
  `Symbol("k$(idx)r")` (Kfor / Krev rendering)

`src/types.jl` — 2 sites (value context with Mechanism in scope; use
the value-context variant):
- Line 1840 — `Symbol("K_$(lig_name)_T_reg$site_idx")` (Kreg T-state)
- Line 1841 — `Symbol("K_$(lig_name)_reg$site_idx")` (Kreg R-state)

**Approach:**
1. Implement the type/index-context variant in `src/types.jl`. The
   variant applies only to step-indexed Parameter types that
   the @generated callers actually use: `Kd`, `Kiso`, `Kfor`,
   `Krev`, `Kon`, `Koff` — each gets a `name(::Type{P}, idx::Int)`
   (and a `name(::Type{P}, idx::Int, state::Symbol)` for T-state)
   returning the same positional Symbol the value-context form
   computes for a step at that rep-index. **`Kreg` does NOT get a
   type/index companion** — Kreg names need a ligand symbol and a
   site, not just an integer index, so it remains value-context only
   (taking a `Kreg` value with `site::RegulatorySite,
   ligand::AllostericRegulator, state::Symbol` fields).

   Both call paths route through a private `_param_symbol(...)`
   formatter family — one method per Parameter shape (3-arg for the
   step-indexed types, 4-arg for Kreg). "Shared formatter" means the
   Symbol-building logic lives in one place per Parameter type, not
   that all parameters share a single `_param_symbol` arity.
2. At each `rate_eq_derivation.jl` site, replace `Symbol("K$idx")`
   with `name(Kd, idx)` (and similarly `Symbol("k$(idx)f")` →
   `name(Kfor, idx)`, `Symbol("k$(idx)r")` → `name(Krev, idx)`).
3. At the two `types.jl` Kreg sites (which sit *inside* the
   `name(p::Kreg, m::AllostericMechanism)` chokepoint body at
   `src/types.jl:1837`), the goal is to factor the Symbol-building
   logic out of the function body into the shared `_param_symbol`
   formatter so the chokepoint regression test sees no top-level
   Symbol literals. Add a `_param_symbol(::Type{Kreg}, site_idx::Int,
   lig_name::Symbol, state::Symbol)` method that returns the same
   Symbol the current inline code does, then rewrite the chokepoint
   body to look up `site_idx = findfirst(==(p.site),
   m.regulatory_sites)` and `lig_name = name(p.ligand)`, and
   delegate to `_param_symbol`. **Do not construct a new `Kreg` —**
   `Kreg(site::RegulatorySite, ligand::AllostericRegulator,
   state::Symbol)` takes struct values, not primitives.

**Regression test:** new `test/test_chokepoint.jl`. Scans `src/*.jl`
line-by-line; any direct `Symbol("K…"` / `Symbol("k…"` /
`Symbol("V…"` / `Symbol("L…"` literal *outside* the body of a `name`
method is a failure. The regex must match both digit-literal forms
(`Symbol("K1")`) and interpolation forms (`Symbol("K$idx")`,
`Symbol("k$(idx)f")`) — the matcher pattern is `Symbol\("[KkVL][_a-zA-Z0-9$]`.

**LOC delta target:** ≈ +20 (new type/index name methods; substitution
at call sites). Value is architectural.

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
6. **Narrow the bare-Symbol arm** (don't delete it outright — a few
   conformation entries like `:E`, `:Estar` are still bare Symbols
   in the new grammar). Allowed bare entries: any Symbol whose
   string form contains no `_` character followed by an uppercase
   metabolite-name character — i.e., conformations like `:E`,
   `:Estar`, `:Estar2`. Any Symbol matching `_[A-Z]` (e.g., `:E_S`,
   `:E_AB`) raises a clear error: "looks like an opaque bound-form
   name; migrate to E(S) or E(A, B) call notation." This catches
   any fixture that slipped migration; without it, an unmigrated
   `:E_S` would silently re-parse as a pure conformation, changing
   semantics.
7. Delete the legacy `(metabolites, reactions)` macro emission path
   (the always-decomposed branch is now the only emission path).
8. Delete opaque-form helpers (`_form_name`, `_parse_bound`,
   `_bound_mets_from_form_name`, `_dead_end_form_name`, `_atoms_dict`,
   `_is_estar_form`, `_can_pingpong`, `_subtract_atoms`) and
   `_mechanism_from_legacy_sig` plus the `_is_new_sig` branch in
   `Mechanism(em::EnzymeMechanism{Sig})`.

**LOC delta target:** net ≈ −200 to −250 (DSL rewrite ≈ neutral;
opaque-form helper deletion ≈ −200; `_mechanism_from_legacy_sig`
deletion ≈ −60; fixture rewrites are tests, not src).

### Stage 7c — Enumeration dispatches on `EnzymeReaction`

**Goal:** the **public-facing** dispatch sites in
`src/mechanism_enumeration.jl` and `src/identify_rate_equation.jl`
that currently take `EnzymeReactionLegacy` move to dispatch on
`EnzymeReaction` (the concrete struct). The **heavy enumeration
pipeline** (`init_mechanisms(::EnzymeReactionLegacy)`,
`_catalytic_topologies`, `_atoms_dict`,
`_apply_equivalence_grouping`, `_expand_substrate_product_dead_ends`)
still dispatches on `EnzymeReactionLegacy` because its body operates
on the singleton type parameters; that conversion is Task 7d.0's
job. So `_to_legacy_reaction` remains called from inside the public
`init_mechanisms(::EnzymeReaction)` adapter until 7d.0 retires it.

**Affected dispatch points (inventory):**
- `src/mechanism_enumeration.jl` ~20 method signatures
- `src/identify_rate_equation.jl`:
  - `IdentifyRateEquationProblem{R<:EnzymeReactionLegacy, D<:NamedTuple}`
    → `{R<:EnzymeReaction, D<:NamedTuple}`
  - `identify_rate_equation(reaction::EnzymeReactionLegacy, ...)` →
    `(reaction::EnzymeReaction, ...)`

**Body rewrites needed for non-heavy-pipeline methods.** For methods
whose dispatch DOES swap to `::EnzymeReaction` in 7c, the function
*bodies* often need updating too: the `EnzymeReactionLegacy`
accessors return tuple shapes
(`substrates(::EnzymeReactionLegacy)` → `Tuple{Tuple{Symbol,
Tuple{...}}}`), while `EnzymeReaction` accessors return reactant
structs (`substrates(::EnzymeReaction)` reads from `r.reactants`).
Internals that destructure `(name, atoms)` tuples must be rewritten
to use the struct accessors.

**Do NOT** rewrite heavy-pipeline methods that stay on Legacy
(`_atoms_dict`, `_catalytic_topologies`, `_apply_equivalence_grouping`,
`_expand_substrate_product_dead_ends`) — their bodies are correct as
written for the Legacy shape and will get rewritten when Task 7d.0
ports the entire heavy path to `EnzymeReaction`. Audit each method
at signature-swap time: if the swap is to `::EnzymeReaction`, body
rewrites needed; if the method stays on `::EnzymeReactionLegacy`,
leave it alone.

**Compile-budget framing:** `test_compile_budget.jl` lines 77, 147,
etc. construct `EnzymeReactionLegacy` directly precisely to detect
per-arity specialization regressions in the `{S,P,R,N}` parametric
form. The non-parametric concrete `EnzymeReaction` cannot exhibit
per-arity specialization at all — that's the *point* of switching:
compilation of uni-uni should not specialize again for ter-ter
because they share the same dispatch. The warmup-reuse test (lines
154–194) is therefore expected to pass trivially after 7c (`t_warm ≈
t_cold`) for the right reasons. This is a success, not a regression.

Concretely for Stage 7c:
- Re-baseline the trace-compile budgets to the new (lower) counts;
  set 2× the new baseline as the gate.
- The warmup-reuse test's ratio bound may need to be relaxed *up*
  (because `t_warm` and `t_cold` are both dominated by tiny constant
  overheads that don't compose into the original "warmup wins"
  pattern) or the test may need a comment update describing that
  it now verifies "no per-arity specialization happens" by direct
  measurement rather than by relative speedup.
- **Do not raise the budgets in a way that would mask a real
  regression** — if you raise a bound, document why the previous bound
  was inappropriate under the new architecture.

**LOC delta target:** small net (dispatch swap + body rewrites; struct
deletions land in 7d).

### Stage 7d — Delete legacy infrastructure

**Goal:** delete all the now-unreachable legacy code paths. After
this stage, `_to_legacy_reaction` has no callers and goes; the
`EnzymeReactionLegacy` struct has no users and goes; the dual-Sig
accessor branches in `src/types.jl` collapse to single bodies.

**Deletion inventory:**

| Item | File | Approx LOC |
|------|------|-----------|
| `EnzymeReactionLegacy{S,P,R,N}` struct (line 628) + outer constructor (lines 654–697) | `types.jl` 628, 654–697 | ~75 |
| `EnzymeReactionLegacy` accessors (`substrates`, `products`, `regulators`, `regulator_roles`, `oligomeric_state`, `Base.show`) | `types.jl` 1221, 1451–1488 | ~80 |
| `_to_legacy_reaction(r::EnzymeReaction)` | `types.jl` 700–723 | ~25 |
| 2-arg `EnzymeMechanism(metabolites, reactions)` constructor | `types.jl` 877–~970 (verify) | ~100 |
| Dual-Sig branches in 12 `EnzymeMechanism` accessors (`_is_new_sig` guards) + `_is_new_sig` helper | `types.jl` (scattered) | ~150 |
| `MechanismSpec` / `StepSpec` / `AllostericMechanismSpec` / `AbstractMechanismSpec` structs + helpers preserved by 6β | `mechanism_enumeration.jl` | ~150 |
| `_init_mechanism_specs(::EnzymeReactionLegacy)` (the heavy-pipeline entry point) | `mechanism_enumeration.jl` 3414 | ~5 |
| `_apply_equivalence_grouping`, `_expand_substrate_product_dead_ends` spec-form helpers (after Task 7d.0 rewrites `init_mechanisms(::EnzymeReactionLegacy)` to produce Mechanism directly) | `mechanism_enumeration.jl` | ~80 |
| **Total** | | **~665** |

**Task 7d.0 (prerequisite):** rewrite the heavy
`init_mechanisms(::EnzymeReactionLegacy)` at
`mechanism_enumeration.jl:1410` to build `Vector{Mechanism}`
directly. Today the body builds `MechanismSpec[]` via
`_apply_equivalence_grouping`; the rewrite either (a) keeps an
internal scratch representation and converts to `Mechanism` at the
end, or (b) reworks the algorithm to manipulate `Mechanism` /
`Vector{Vector{Step}}` throughout. This is the deletion-gating change
— without it, `MechanismSpec` cannot be retired.

**LOC delta target:** −650 to −750. Stage 7d is the single largest
deletion stage of the continuation.

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

Cumulative is relative to main (7,136). Start of continuation:
+2,607 (current 9,743 src LOC).

| Stage | Target Δ | Cumulative |
|-------|----------|------------|
| Start | — | +2,607 |
| 6β | −300 to −500 | ≈ +2,200 |
| 6 | −300 to −500 | ≈ +1,800 |
| 7a | +10 to +30 | ≈ +1,820 |
| 7b | −200 to −250 | ≈ +1,600 |
| 7c | ≈ 0 | ≈ +1,600 |
| 7d | −650 to −750 | ≈ +900 |
| 7e | −100 | ≈ +800 |
| **End** | | **≈ +800 (≈ 7,940; gate ≤ 8,200)** |

If a stage misses its target, the next stage absorbs the slack. The
≤8,200 final src LOC is the only gate; per-stage targets are guides.
Spec types' deletion shifted from 6β to 7d (sequencing fix — the heavy
enumeration pipeline still constructs `MechanismSpec` internally
until rewritten in Task 7d.0), which is why 6β's delta range is now
smaller and 7d's larger than prior drafts.

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
4. **`@enzyme_reaction`'s `regulators:` alias dropped** — the
   `@enzyme_reaction` macro rejects `regulators:` (only
   `dead_end_inhibitors:`, `competitive_inhibitors:`,
   `allosteric_regulators:` accepted, per
   `src/dsl.jl` `_VALID_REACTION_LABELS`). `@enzyme_mechanism` and
   `@allosteric_mechanism` still accept `regulators:` (~70 fixture
   sites use it) — no change planned for those macros in this
   continuation.
5. **Memory note correction** —
   `project_dedup_pass2_dead_code` was wrong; pass-2 Wegscheider
   absorption IS load-bearing for inhibitor/activator mechanisms.
   `_build_kinetic_rename_map` is kept.

## 10. Finishing-phase addendum (2026-05-25)

The previous session landed Stages 6β–7e but deferred three coupled
cleanups: the opaque-form enumeration helpers, the 5 opaque-form Segel
fixtures, and the legacy DSL emission + dual-Sig accessor branches. This
addendum records the decisions that close them out.

### 10.1 Lumped central-complex fixtures migrate by rename, not re-derivation

The 5 fixtures (`:EABEPQ`, `:EAFP`/`:FBEQ`, `:EABCEPQ`, `:EABCEPQR`) use
lumped central-complex node names. The deferred-work framing assumed
moving them off opaque grammar required either re-deriving the textbook
analytical formula for an explicit-iso split form, or replacing the
analytical-comparison test with a synthetic-data fit.

**Both are unnecessary.** A spike (Segel Ordered Bi Bi built both ways)
showed the derived `rate_equation_string` is byte-for-byte identical
between the opaque `EABEPQ` form and the decomposed single-node `E(A,B)`
form: same `n_states`, `n_steps`, `fitted_params`, Haldane constraint,
and denominator. Migration is therefore a **mechanical rename** of the
lumped node to a single decomposed node (`EABEPQ → E(A,B)`,
`EAFP → E(A)`, `FBEQ → F(B)`, …); the textbook formula is preserved.

Why it works: `stoich_matrix` (`src/types.jl`) treats enzyme forms as
opaque species *by name* (no atom-tracking of bound contents), so a
fused catalytic-release step `E(A,B) <--> E(Q) + P` — where `P` was never
bound in `E(A,B)` — is a valid stoichiometric column. Release direction
is inferred by `_legacy_step_tuple`'s bound-list-size fallback
(`length(from_bound) > length(to_bound)`). This keeps the lumped 4-step
form 4 steps in decomposed grammar instead of forcing the 5-step
explicit-iso split.

### 10.2 No opaque-form escape hatch (Stage 7b.6 confirmed, reasoning corrected)

The opaque bound-form grammar (`:EABEPQ`, `:E_S`) is removed from the
macros. Users must write decomposed forms (`E(A,B)`, `E(P,Q)`)
explicitly. §5 Stage 7b.6's "reject opaque bound-form Symbols at parse
time" decision stands; only its premise changes — rejecting opaque forms
does NOT force fixture re-derivation (per §10.1), so there is no cost to
deleting the opaque path outright.

### 10.3 Revised finishing sequence

Stages collapse from the handoff's A→B→C into two phases:

**Phase 1 — fixture rename + legacy-path removal** (low risk):
1. Rename the 5 lumped fixtures to single decomposed nodes; their
   analytical-rate tests stay green.
2. DSL rejects opaque bound-form Symbols with a clear migration error.
3. Delete the legacy DSL emission path, the 2-arg
   `EnzymeMechanism(metabolites, reactions)` constructor, and
   `_mechanism_from_legacy_sig`.
4. Collapse the 12 `_is_new_sig` dual-Sig accessor branches to the
   new-shape body; delete `_is_new_sig` and `_legacy_step_tuple` if dead.

**Phase 2 — enumeration Symbol→struct rewrite** (high risk, last):
Rewrite the topology backtracker + dead-end enumeration to operate on
`Step`/decomposed `Species` instead of opaque `Symbol` working-rep;
delete the parse-back helpers (`_parse_bound`,
`_bound_mets_from_form_name`, `_dead_end_form_name`, `_is_estar_form`)
and the topology-backtracker atom helpers as their callers disappear.

All §8 non-negotiables (test integrity, perf gates, compile-budget,
no `--amend`, no temporal-context comments, chokepoint exclusivity)
apply unchanged. Verify the full suite + integrity check + compile-budget
+ `rate_equation` 0-alloc/<100ns gate after every commit.

## 11. Execution findings (2026-05-25) — Phase 1 is harder than the spike implied

The §10 plan assumed the only opaque fixtures were the 5 lumped
single-metabolite ones, migratable by rename. Executing it surfaced
three architectural constraints and a much larger fixture set. These are
established facts (spike + full-suite diagnostic), not hypotheses; the
revised plan must be built on them.

### 11.1 The new Step-based Sig is lossy for opaque Species

The new Sig encodes each `Step` as `(from_species, to_species,
bound_metabolite, is_equilibrium)` and reconstructs the reaction
`(lhs, rhs)` via `_legacy_step_tuple_from_sig`, which infers
binding-vs-release direction from bound-list membership:

- **Fully-decomposed** Species (explicit bound lists) → direction is
  recoverable losslessly (metabolite in `to`-bound but not `from` ⇒
  binds; in `from` but not `to` ⇒ releases).
- **Opaque** Species (empty bound lists, e.g. `Species([], :EABEPQ)`) →
  the metabolite is in *neither* bound list, so direction is *guessed*
  (default: binding). Product-**release** steps of opaque central
  complexes (`EABEPQ <--> EQ + P`) get mis-reconstructed as binding
  (`EABEPQ + P → EQ`), corrupting stoichiometry. The derivation then
  fails with `Cycle N produces metabolite change not proportional to
  net reaction`.

**Consequence:** disabling opaque rejection to route opaque fixtures
through the new Sig as opaque Species does **not** yield green — the
opaque central complexes crash. The legacy Sig avoided this by storing
the explicit reaction tuples; the new Sig cannot. Therefore every opaque
fixture **must be decomposed** (or deleted); there is no opaque-Species
shortcut, and legacy-path removal is coupled to completing fixture
decomposition.

### 11.2 Three mechanism shapes and how they resolve

Full-suite diagnostic (rejection off) → 54 errors across ~7 fixtures:

| Shape | Fixtures | Resolution |
|---|---|---|
| Bind+release in one step (`EA + B <--> EQ + P`) | Theorell-Chance | **Un-representable** — `Step` has one `bound_metabolite`, can't encode two. **Delete** (Denis-approved §2.1). |
| Opaque central complex, single-metabolite steps | Ordered Uni Bi (`EAEPQ`), RE Ordered Bi-Bi (`EABEPQ`), Segel Bi-Uni-Uni-Uni PP Ter Bi (`EABFP`/`FCEQ`), Bi-Uni-Uni-Bi PP Ter Ter (`FCEQR`), Bi-Bi-Uni-Uni PP Ter Ter (`EABFPQ`/`FCER`), Hexa Uni PP (`EAFP`/`FBGQ`/`GCER`) | **Decompose** by renaming each lumped node to a single decomposed node (the spike-validated pattern). ~6 fixtures. Tractable. |
| Multi-site metabolite | HK — G6P binds as **Product** (catalytic site), **CompetitiveInhibitor** (nucleotide pocket), and **AllostericRegulator** (regulatory site) simultaneously; forms `E_G6P` vs `E_G6Pi` vs `E_G6P_G6Pi` | **Build a DSL feature** to express "metabolite bound at site/role X vs Y" in decomposed grammar (the model already distinguishes via `Metabolite` subtype in the bound list — `Product(:G6P)` vs `CompetitiveInhibitor(:G6P)` — but the grammar + one-role-per-name `role_of` can't). Denis: keep HK, modify code. |

### 11.3 Revised Phase 1 sequence (supersedes §10.3 Phase 1)

1. **Design + build the multi-site-binding DSL feature** (HK's blocker;
   warrants its own brainstorm/spec — syntax for per-binding role +
   `role_of`/`_build_step_expr` changes). Hardest; do first.
2. **Decompose all opaque fixtures** to decomposed grammar: the ~6
   single-metabolite lumped fixtures by rename; HK via the new feature;
   delete Theorell-Chance (§2.1). Suite green.
3. **Delete the legacy Sig path** (2-arg `EnzymeMechanism(metabolites,
   reactions)`, `_mechanism_from_legacy_sig`, the 12 `_is_new_sig`
   accessor branches) — now safe, nothing produces legacy-shape Sig.
4. **Re-enable opaque rejection** in both macros as the final gate.

### 11.4 Process note

Delegating exploratory fixture migration to autonomous subagents failed
(one reverted the architecture to escape a blocker instead of
escalating). This work is exploratory, not mechanical: drive it directly
or with tightly-scoped, closely-monitored single-fixture subagents, and
escalate architectural blockers rather than working around them.

### 11.5 Status (branch `refactor-to-concrete-types-instead-of-symbols`)

Phase 2 enumeration internals are complete: `init_mechanisms` and
`expand_mechanisms` emit decomposed `Mechanism` / `AllostericMechanism`
(built from `Step` / `Species`) end-to-end. The entire opaque-form
`_Raw*` working-representation family and all form-name string helpers
(`_dead_end_form_name`, `_is_estar_form`, `_bound_mets_from_form_name`,
`_parse_bound`, `_mechanism_spec_from_steps`, `_stepspec_from_step`,
`_mechanism_from_raw`, `_raw_from_mechanism`, `_form_name`, etc.) are
deleted. §3.1 #1 (one struct family in the enumerator) and #2 (no
form-name parse-back) hold for the enumeration path. The 2-arg
`EnzymeMechanism(metabolites, reactions)` legacy Sig path still exists in
`src/types.jl` (unreachable from the DSL); its removal is the separate
§11.3 step 3, not yet done.
