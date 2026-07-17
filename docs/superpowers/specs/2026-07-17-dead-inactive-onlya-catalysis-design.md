# Dead-inactive `:OnlyA` catalysis — an `:OnlyA` binding kills all inactive catalysis

**Date:** 2026-07-17
**Status:** Design, approved (design decisions settled with Denis; see "Decisions").
Resolves the err2 open question left by
`docs/superpowers/specs/2026-07-16-mwc-derivation-targeted-fixes-design.md`
(deliverable 2, reported BLOCKED) and the ping-pong `:OnlyA` `_kcat_forward` crash in
`docs/superpowers/findings/2026-07-16-pingpong-onlya-kcat-bug.md` — including its
V-system class.

## Summary

When the inactive conformation cannot bind a catalytic substrate/product (an `:OnlyA`
binding), it cannot complete the catalytic cycle, so the correct MWC reading is that it
is **catalytically dead**: every chemical (isomerisation) step is `:OnlyA`. The current
enumeration instead promotes the *fewest* chemical steps that satisfy the Haldane
relation — one — which for a multi-chemical-step (ping-pong) mechanism produces a
**partial-catalysis** form that has a kinetic sink and crashes `_kcat_forward`. The
thermodynamically-complete dead-inactive form is never generated, because it is
non-minimal.

This design makes `_expand_to_allosteric` emit, for every mechanism, all valid
dead-inactive `:OnlyA`-binding combinations directly — each with **all** chemical steps
`:OnlyA` — and never a partial-catalysis form. It is enforced **in the enumeration logic
only** (no constructor guard).

## Decisions

Settled with Denis before writing this spec:

1. **Dead-inactive is the only allosteric-catalysis form under an `:OnlyA` binding.** Any
   `:OnlyA` catalytic binding ⟹ every chemical (iso) step `:OnlyA` — a fully
   catalytically-dead inactive conformation (the standard MWC tense state). No
   partial-catalysis form is generated. Scope is **all chemical steps in the mechanism**,
   not a per-cycle subset.
2. **`_expand_to_allosteric` emits all valid `:OnlyA`-binding combos directly.** For each
   mechanism it emits every non-empty thermodynamically-valid subset of binding groups as
   `:OnlyA` (bare — the `:OnlyA` binding's metabolite reveals `L`), each with all chemical
   steps `:OnlyA`; plus the empty subset paired with a declared regulator (the V-type —
   now also all-chem-`:OnlyA`, not one). This subsumes the completion search and the
   V-system case. Bounded at `2^B` (B = number of binding groups): measured **15 for
   ping-pong bi-bi** (all valid), **3 for uni-uni** — small because pinning all chemical
   steps `:OnlyA` removes the chemical-tag dimension.
3. **`_expand_promote_catalytic_to_onlya` is dropped — contingent on an empirical coverage
   check** (see Architecture §3). With to_allosteric emitting every combo, promote's
   incremental `:OnlyA`-adding is redundant; but dropping it is verified not to lose
   coverage before it lands, because a silent loss is exactly the failure this work
   exists to prevent.
4. **Enforcement = enumeration logic only.** No `AllostericMechanism` constructor
   backstop. A partial-catalysis form is never *enumerated*, but a hand-written
   `@allosteric_mechanism` (or a `Sig` round-trip of one) that encodes a partial still
   constructs — and still crashes `_kcat_forward`. Accepted: the HPC failure is an
   enumeration problem, and production never hand-builds these.

## The problem, verified

`_valid_onlya_completions` (`src/mechanism_enumeration.jl:1776`) searches subsets of the
`:EqualAI` groups by **increasing size** and returns every valid vector at the *first*
size that satisfies `_onlya_haldane_violation` (`:1790`). For the PFK-P ping-pong,
promoting the F6P binding (g6) to `:OnlyA` returns four **size-1** completions and stops:

```
1. {bind(ADP),   bind(F6P)}      opposing-binding repair (no chemical step)
2. {chem1,       bind(F6P)}      PARTIAL — 1 of 2 chemical steps  (chem1-only)
3. {chem2,       bind(F6P)}      PARTIAL — 1 of 2 chemical steps  ← err2, the sink
4. {bind(F16BP), bind(F6P)}      opposing-binding repair (no chemical step)
```

The dead-inactive form — F6P binding + **both** chemical steps `:OnlyA` — is size-2, so
it is **never generated** (`scratchpad/completions.jl`).

Why the partials are wrong (measured, `scratchpad/err2_thermo.jl`, shipped derivation,
Keq = 3):

| variant | guard | `d_free_I` | v(norm) | v(F16BP→0) | v(equil) | kcat |
|---|---|---|---|---|---|---|
| chem2-only (err2, enumerated) | admit | `F16BP·k/K` | 0.0498 | **0.0** | 0.0 | **CRASH** |
| **both chem `:OnlyA` (dead-inactive)** | admit | **1** | 0.1045 | **0.145** | 0.0 | **1.3** |
| chem1-only | admit | 1 | 0.0966 | 0.133 | 0.0 | 1.3 |

- Every variant is thermodynamically consistent (v = 0 at equilibrium) — so this is
  **not** a detailed-balance inconsistency, and PR #70's guard correctly admits them all.
  It is a weaker, consistent-but-degenerate failure the Haldane guard cannot see.
- Only chem2-only **sinks**: v → 0 as F16BP → 0 while v(norm) is finite. The covalent
  intermediate `E~` in the inactive conformation is fed by the live chemical step and can
  only drain by rebinding a product, unavailable at zero product, so all enzyme
  accumulates there and the forward turnover at products = 0 is genuinely zero — which is
  why `_kcat_forward` finds no saturating pattern and throws. The throw is a *true report*
  about the partial form, not a derivation bug.
- The opposing-binding repairs (1, 4) keep the inactive cycle formally balanced but the
  inactive state still cannot bind F6P, so it still cannot turn over — the finding doc
  counted 14 of these "balanced K-system" forms among the 112 crashes. The V-system class
  (14 more) is the same disease from a lone `:OnlyA` chemical step.

The dead-inactive form is the only sound one: `d_free_I = 1`, kcat finite (1.3), **no
sink**, v = 0 at equilibrium — the canonical MWC tense state (binds ligands, does no
chemistry).

## Why dead-inactive derives cleanly (dependency on the residual-island fix)

The dead-inactive form's `d_free_I = 1` is produced by the residual-island stranding
already merged on `mwc-targeted-fixes` (the err1 fix, commit `29eee7e`,
`_reachable_from_free` seeds only forms with empty `bound` **and** empty `residual`).
With all chemical steps `:OnlyA`, the inactive graph drops them, the covalent branch is
reachable only through those (now-dropped) steps, so it strands and the inactive
conformation collapses to `{E, E·ATP, E·ADP}` with `d_free_I = 1` (verified,
`scratchpad/prereq.jl`). **This design depends on the err1 stranding fix being present**
(branch off `mwc-targeted-fixes`, or off `main` once it merges).

## Architecture

All changes are in `src/mechanism_enumeration.jl`, plus tests. No `@generated` codegen
changes, so the `rate_equation` 0-allocation / <120 ns contract is untouched.

### 1. `_expand_to_allosteric` — emit all dead-inactive `:OnlyA`-binding combos

Rewrite the per-group emission (`:1834-1859`) into a per-subset emission. Let `iso` be
the iso (chemical) groups and `bind` the binding groups. For each `catalytic_multiplicity`:

- **Non-empty binding subset `S ⊆ bind`:** build the tag vector with every group in `S`
  and every group in `iso` set to `:OnlyA`, the rest `:EqualAI`; keep it iff
  `_onlya_haldane_violation` passes (drops the binding-Wegscheider-inconsistent subsets on
  random-order squares); emit it bare (no regulator — the `:OnlyA` binding's metabolite
  reveals `L`). This is the K-type family, now spanning every valid subset.
- **Empty binding subset (V-type):** build the tag vector with every `iso` group `:OnlyA`
  and all bindings `:EqualAI`; with no `:OnlyA` binding, `L` folds into `kcat` and is
  unobservable, so emit it only paired with a declared allosteric regulator — one variant
  per `(regulator, tag ∈ {:OnlyA, :OnlyI})`, exactly as the current V-type path does, but
  with **all** chemical steps `:OnlyA` instead of the single promoted one.

The observability rules are preserved (all-`:EqualAI` never emitted; a bare form needs an
`:OnlyA` binding to reveal `L`; a V-type needs a regulator). The four-way minimal-subset
completion is gone.

`_valid_onlya_completions` and its `_each_subset` helper are removed if no caller remains
(see §3); the dead-inactive tag vector is constructed directly and validated with the
existing `_onlya_haldane_violation`.

### 2. `_expand_change_allo_state` — filter partial-creating relaxations

This move relaxes a constrained tag (`:EqualAI`/`:OnlyA`/`:OnlyI`) to `:NonequalAI` and
already drops relaxations that violate `_onlya_haldane_violation` (`:2020`). It must
**also** drop any relaxation that yields a partial-catalysis form — an `:OnlyA` binding
present while some iso group is not `:OnlyA` — because relaxing one chemical step of a
dead-inactive mechanism recreates exactly the err2 sink, and with no constructor guard
nothing downstream catches it. Add a `_partial_onlya_catalysis(cat_steps, cat_allo_states)
→ Bool` predicate (an `:OnlyA` binding group coexisting with a non-`:OnlyA` iso group) and
skip any relaxation whose result satisfies it. Relaxing the `:OnlyA` **binding** itself to
`:NonequalAI` is still allowed (it leaves all chemical steps `:OnlyA`, no `:OnlyA` binding
— a valid "dead inactive that binds the ligand with different affinity" form).

### 3. `_expand_promote_catalytic_to_onlya` — drop, contingent on a coverage diff

With §1 emitting every dead-inactive combo, promote's incremental `:OnlyA`-adding is
redundant: promoting a binding to `:OnlyA` on any mechanism now yields a dead-inactive
form that to_allosteric already emits (on the modify-then-allosterize path). Drop the move
**only after** confirming coverage is preserved:

- Enumerate a spread of reactions (uni-uni, ordered bi-bi, ping-pong bi-bi, a
  ter-substrate) to a fixed depth **both ways** — current move set, and new §1 with
  promote removed — and diff the set of enumerated `AllostericMechanism`s carrying any
  `:OnlyA` tag.
- If the diff is empty, drop promote. If any mechanism is lost, **keep** promote (it is
  harmlessly redundant once its completion uses the dead-inactive construction — promoting
  a binding forces all-chem-`:OnlyA`, a form to_allosteric also emits, so dedup absorbs it)
  and record which mechanisms needed it.

This is the one place a silent coverage loss could hide, so it is gated on measurement,
not reasoning.

### Relationship to the existing guards

`_onlya_haldane_violation` (the Stiemke ε-feasibility check completed on
`mwc-targeted-fixes`) stays, unchanged, as the per-subset validity check in §1. It is
**complementary**, not subsumed: with all chemical steps `:OnlyA` the inactive *catalytic*
cycle is gone, but an `:OnlyA` binding on a random-order **binding** square can still leave
a binding-only Wegscheider cycle unsatisfiable, which is exactly what it catches (and what
prunes some subsets in §1).

## Scope and migration

- **Single-chemical-step mechanisms are effectively unchanged.** For uni-uni, ordered
  bi-bi, and the LDH i-state targets, "all chemical steps `:OnlyA`" equals "the one
  chemical step `:OnlyA`", so the emitted `:OnlyA`-binding subsets match the current
  completion's — with one cosmetic difference: the current opposing-binding form (all
  bindings `:OnlyA`, chemical step `:EqualAI`) becomes all-`:OnlyA`, but that form's
  inactive conformation binds nothing and so is inert regardless of the chemical tag, so
  the two are **rate-equivalent** and dedup to one. Confirm this empirically (the golden
  byte-identical and the enumeration counts unchanged for these reactions); the rule only
  changes the equations of **multi-chemical-step (ping-pong)** mechanisms — the crash
  class.
- **Enumeration counts change** for ping-pong reactions: the partial/opposing/V-system
  forms are replaced by the dead-inactive combo family (≤ `2^B` per topology). Re-measure
  and update any affected `expected_n_*` counts; a change must be explained by this rule.
- **The allosteric golden reference:** verify byte-identical. The static
  `MECHANISM_TEST_SPEC`s contain no allosteric ping-pong-with-`:OnlyA` mechanism (the
  Segel ping-pongs are non-allosteric), so the golden is expected unchanged; if a block
  moves, stop and explain.
- **err1, err2, and the V-system crashes cease to be enumerated** as partial forms; the
  dead-inactive combos replace them. err1's stranding fix stays (load-bearing for the
  dead-inactive form and for any hand-written residual-island mechanism).

## Testing

- **All-combos gate.** For a ping-pong bi-bi, `_expand_to_allosteric` emits exactly the
  set of valid dead-inactive `:OnlyA`-binding subsets (15 measured), each with both
  chemical steps `:OnlyA` and no partial. Assert the count and that every emitted `:OnlyA`
  mechanism has all iso groups `:OnlyA`.
- **Derivation / ground-truth gate.** A dead-inactive ping-pong `:OnlyA` mechanism derives
  with `d_free_I = 1`, finite `rate_equation` and finite `_kcat_forward`, **no sink** (v
  finite as a product → 0), and v = 0 at equilibrium — matched against a two-conformation
  mass-action ground truth in `test/allosteric_ground_truth.jl` (the inactive conformation
  is catalytically dead, contributing only `L` to the denominator).
- **`_expand_change_allo_state` gate.** Relaxing a chemical step of a dead-inactive
  `:OnlyA`-binding mechanism is dropped; relaxing the `:OnlyA` binding to `:NonequalAI` is
  retained.
- **Promote-drop coverage diff** (§3): the both-ways enumeration diff is empty for the
  tested reactions (or promote is kept with the mechanisms that needed it recorded).
- **Regression.** Single-chemical-step allosteric mechanisms unchanged; the 12+
  `allosteric_ground_truth.jl` gates green; the allosteric golden byte-identical; the
  `rate_equation` perf gate (0 alloc, <120 ns) green; enumeration-count tests updated to
  measured values with the rule as the explanation.

## Risks

- **Enumeration coverage (promote drop).** The central risk, gated on the §3 diff — drop
  only if measured lossless.
- **Combo growth on large reactions.** `2^B` is 15 for ping-pong bi-bi but grows with
  binding groups; a ter-substrate random-order mechanism has more, partly offset by
  Wegscheider pruning. The plan must measure the combo count for the largest reaction in
  play and confirm the search stays tractable; if not, revisit whether all subsets are
  wanted or only a physically-motivated subset.
- **No constructor backstop (decision 4).** `_expand_change_allo_state` is the only place a
  relaxation could recreate a partial, so its filter is load-bearing — a missed case there
  re-introduces an enumerated sink with no error. The `_expand_change_allo_state` gate
  guards against that.
- **Multi-cycle mechanisms.** "All chemical steps in the mechanism" over-restricts a
  hypothetical mechanism with two independent catalytic cycles where one should stay
  active. The plan confirms the enumeration produces no such mechanism (ping-pong is a
  single linear cycle); if one exists, revisit decision 1's scope.

## Evidence

Scratchpad scripts: `completions.jl` (the four minimal completions; dead-inactive absent),
`combos.jl` (15 valid dead-inactive subsets for ping-pong, 3 for uni-uni),
`err2_thermo.jl` (the four-variant sink/consistency table), `err2_scope.jl` (the
inactive-graph trap trace), `prereq.jl` (dead-inactive `d_free_I = 1` via stranding),
`a_form.jl` (the active-form cycle).
