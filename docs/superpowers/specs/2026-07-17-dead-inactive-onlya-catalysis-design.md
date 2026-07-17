# Dead-inactive `:OnlyA` catalysis — an `:OnlyA` binding kills all inactive catalysis

**Date:** 2026-07-17
**Status:** Design, approved (3 design decisions settled with Denis; see "Decisions").
Resolves the err2 open question left by
`docs/superpowers/specs/2026-07-16-mwc-derivation-targeted-fixes-design.md`
(deliverable 2, reported BLOCKED) and the ping-pong `:OnlyA` `_kcat_forward` crash in
`docs/superpowers/findings/2026-07-16-pingpong-onlya-kcat-bug.md`.

## Summary

When the enumerator promotes a catalytic **substrate/product binding** to `:OnlyA`
(the inactive conformation cannot bind that metabolite), the inactive conformation can
no longer complete the catalytic cycle. The correct MWC reading is that the inactive
conformation is then **catalytically dead**: every chemical (isomerisation) step
becomes `:OnlyA`. The current enumeration instead promotes the *fewest* chemical steps
that satisfy the Haldane relation — one — which for a multi-chemical-step (ping-pong)
mechanism produces a **partial-catalysis** form that has a kinetic sink and crashes
`_kcat_forward`. The thermodynamically-complete dead-inactive form is never generated,
because it is non-minimal.

This design changes the enumeration completion so that an `:OnlyA` binding yields the
single dead-inactive form (all chemical steps `:OnlyA`), and it never emits a
partial-catalysis form. It is enforced **in the enumeration logic only** — no
constructor guard.

## Decisions

Settled with Denis before writing this spec:

1. **Completion rule = dead-inactive only.** Promoting a catalytic binding to `:OnlyA`
   yields exactly one completion: that binding plus **all** chemical (iso) steps
   `:OnlyA`. The current partial-chemical repairs *and* the opposing-binding repairs are
   dropped. Further `:OnlyA` bindings are reached by iterating the promote move, not by
   emitting binding subsets at once.
2. **Scope = all chemical steps in the mechanism** (a fully catalytically-dead inactive
   conformation — the standard MWC tense state), not a per-cycle subset.
3. **Enforcement = enumeration logic only.** No `AllostericMechanism` constructor
   backstop. A partial-catalysis form is never *enumerated*, but a hand-written
   `@allosteric_mechanism` (or a `Sig` round-trip of one) that encodes a partial still
   constructs — and still crashes `_kcat_forward`. That is accepted: the HPC failure is
   an enumeration problem, and production never hand-builds these.

## The problem, verified

`_valid_onlya_completions` (`src/mechanism_enumeration.jl:1776`) searches subsets of the
`:EqualAI` groups by **increasing size** and returns every valid vector at the *first*
size that satisfies `_onlya_haldane_violation` (`isempty(found) || return found`,
`:1790`). For the PFK-P ping-pong, promoting the F6P binding (g6) to `:OnlyA` returns
four **size-1** completions and stops:

```
1. {bind(ADP),   bind(F6P)}      opposing-binding repair (no chemical step)
2. {chem1,       bind(F6P)}      PARTIAL — 1 of 2 chemical steps  (chem1-only)
3. {chem2,       bind(F6P)}      PARTIAL — 1 of 2 chemical steps  ← err2, the sink
4. {bind(F16BP), bind(F6P)}      opposing-binding repair (no chemical step)
```

The dead-inactive form — F6P binding + **both** chemical steps `:OnlyA` — is size-2, so
it is **never generated**. Measured (`scratchpad/completions.jl`).

Why the partials are wrong (measured, `scratchpad/err2_thermo.jl`, shipped derivation,
Keq = 3):

| variant | guard | `d_free_I` | v(norm) | v(F16BP→0) | v(equil) | kcat |
|---|---|---|---|---|---|---|
| chem2-only (err2, enumerated) | admit | `F16BP·k/K` | 0.0498 | **0.0** | 0.0 | **CRASH** |
| **both chem `:OnlyA` (dead-inactive)** | admit | **1** | 0.1045 | **0.145** | 0.0 | **1.3** |
| chem1-only | admit | 1 | 0.0966 | 0.133 | 0.0 | 1.3 |

- Every variant is thermodynamically consistent (v = 0 at the equilibrium metabolite
  ratio) — so this is **not** a detailed-balance inconsistency, and PR #70's guard
  correctly admits them all. It is a weaker, consistent-but-degenerate failure the
  Haldane guard cannot see.
- Only chem2-only **sinks**: v → 0 as F16BP → 0 while v(norm) is finite. The covalent
  intermediate `E~` in the inactive conformation is fed by the live chemical step and
  can only drain by rebinding a product, which is unavailable at zero product, so all
  enzyme accumulates there and the forward turnover at products = 0 is genuinely zero —
  which is why `_kcat_forward` finds no saturating pattern and throws. The throw is a
  *true report* about the partial form, not a derivation bug.
- The opposing-binding repairs (1, 4) keep the inactive cycle formally balanced but the
  inactive state still cannot bind F6P, so it still cannot turn over — the finding doc
  counted 14 of these "balanced K-system" forms among the 112 crashes.

The dead-inactive form is the only sound one: `d_free_I = 1`, kcat finite (1.3), **no
sink** (v stays finite as F16BP → 0), v = 0 at equilibrium. It is the canonical MWC
tense state — an inactive conformation that binds ligands (possibly with different
affinity) but does no chemistry.

## Why dead-inactive derives cleanly (dependency on the residual-island fix)

The dead-inactive form's `d_free_I = 1` is produced by the residual-island stranding
already merged on `mwc-targeted-fixes` (the err1 fix, commit `29eee7e`,
`_reachable_from_free` seeds only forms with empty `bound` **and** empty `residual`).
With all chemical steps `:OnlyA`, the inactive graph drops them, the covalent branch is
reachable only through those (now-dropped) chemical steps, so it strands and the
inactive conformation collapses to `{E, E·ATP, E·ADP}` with `d_free_I = 1` (verified,
`scratchpad/prereq.jl`). **This design therefore depends on the err1 stranding fix being
present** (branch off `mwc-targeted-fixes`, or off `main` once it merges).

## Architecture

All changes are in `src/mechanism_enumeration.jl`, plus tests. No `@generated` codegen
changes, so the `rate_equation` 0-allocation / <120 ns contract is untouched.

### 1. A partial-catalysis predicate (helper, used by the moves — not the constructor)

Add `_partial_onlya_catalysis(cat_steps, cat_allo_states) → Bool`: true when some
catalytic group is an `:OnlyA` **binding** (its representative step binds a metabolite)
while some catalytic **iso** group is not `:OnlyA`. This is the structural signature of
a live-catalysis inactive conformation under an `:OnlyA` binding — the sink class. It is
a cheap tag/step scan, no linear algebra. Per decision 3 it is **not** wired into the
`AllostericMechanism` constructor; it is used only inside the enumeration moves below.

### 2. `_valid_onlya_completions` — dead-inactive completion

Rewrite the completion so that when `tags` contains an `:OnlyA` **binding**, the
completion promotes **every** iso group to `:OnlyA` (preserving existing `:OnlyA` /
`:NonequalAI` tags on bindings), and returns that single vector. Concretely:

- If `_onlya_haldane_violation(tags) === nothing` and `tags` is not a partial (no
  `:OnlyA` binding, or already all-chem-`:OnlyA`), return `[tags]` unchanged (today's
  fast path for the non-`:OnlyA`-binding cases).
- If `tags` has an `:OnlyA` binding, set every iso group to `:OnlyA` and return that one
  vector. With all chemical steps `:OnlyA` the inactive cycle is broken, so the Haldane
  relation is satisfied by construction — no subset search. (A residual binding-only
  Wegscheider inconsistency is still possible on a random-order binding square; keep the
  final `_onlya_haldane_violation` check and return `[]` if the dead-inactive vector
  still violates it — that is a genuinely impossible mechanism.)

The four-way minimal-subset search is removed. This is where the partials and the
opposing-binding repairs disappear.

### 3. `_expand_change_allo_state` — filter partial-creating relaxations

This move relaxes a constrained tag (`:EqualAI`/`:OnlyA`/`:OnlyI`) to `:NonequalAI` and
already drops relaxations that violate `_onlya_haldane_violation` (`:2020`). It must
**also** drop any relaxation that yields a partial (`_partial_onlya_catalysis`), because
relaxing one chemical step of a dead-inactive mechanism recreates exactly the err2 sink
and — with no constructor guard — nothing downstream would catch it. Relaxing the
`:OnlyA` **binding** itself to `:NonequalAI` is still allowed (it leaves all chemical
steps `:OnlyA`, no `:OnlyA` binding, no partial — a valid "dead inactive that binds the
ligand with different affinity" form).

### Relationship to the existing guards

`_onlya_haldane_violation` (the Stiemke ε-feasibility check completed on
`mwc-targeted-fixes`) stays, unchanged. It is **complementary**, not subsumed: with all
chemical steps `:OnlyA` the inactive *catalytic* cycle is gone, but an `:OnlyA` binding
on a random-order **binding** square can still leave a binding-only Wegscheider cycle
unsatisfiable, which is exactly what that guard catches. The new completion keeps calling
it as the final validity check.

## Scope and migration

- **Single-chemical-step mechanisms are unchanged.** For uni-uni, ordered bi-bi, and the
  LDH i-state targets, "all chemical steps `:OnlyA`" equals "the one chemical step
  `:OnlyA`" equals today's minimal completion. The rule only changes **multi-chemical-step
  (ping-pong)** mechanisms — precisely the crash class.
- **Enumeration counts change for ping-pong reactions:** the four partial/opposing
  completions per promoted binding collapse toward one dead-inactive form (more `:OnlyA`
  bindings still reachable by iterating). Re-measure and update any affected
  `expected_n_*` counts; a change must be explained by this rule.
- **The allosteric golden reference:** verify byte-identical. The static
  `MECHANISM_TEST_SPEC`s contain no allosteric ping-pong-with-`:OnlyA` mechanism (the
  Segel ping-pongs are non-allosteric), so the golden is expected unchanged; if a block
  moves, stop and explain.
- **err1 and err2 both cease to be enumerated** as partial forms; the dead-inactive form
  replaces them. err1's stranding fix stays (it is load-bearing for the dead-inactive
  form and for any hand-written residual-island mechanism).

## Open sub-decision (flagged, not decided here): V-system partials

Decision 1's rule triggers on an `:OnlyA` **binding**. There is a symmetric case it does
not cover: a **V-system** — `_expand_to_allosteric` (`:1848-1853`) promotes a single
chemical step to `:OnlyA` paired with a regulator, with the other chemical step(s) left
`:EqualAI`. For a ping-pong that is a partial-catalysis V-system, and the finding doc
counted 14 "V-system" forms among the 112 crashes. The clean generalisation is "inactive
catalysis is all-or-nothing" — promoting **any** catalytic step (binding or chemical) to
`:OnlyA` forces all chemical steps `:OnlyA`. Whether to extend the rule to the V-system
trigger is left for a follow-up: the sink physics for a lone-`:OnlyA`-chemical V-system
has not been reproduced firsthand here (an earlier hand-trace suggested a chem1-only
V-system may *not* sink, unlike the binding-triggered case), so it needs its own
measurement before the rule is broadened. This spec deliberately scopes to the
binding trigger Denis specified.

## Testing

- **Enumeration gate.** `_valid_onlya_completions` on a ping-pong with a promoted
  substrate binding returns exactly the dead-inactive vector (all iso groups `:OnlyA`),
  and no partial. Assert against the concrete PFK-P case in `scratchpad/completions.jl`.
- **Derivation / ground-truth gate.** The dead-inactive ping-pong `:OnlyA` mechanism
  derives with `d_free_I = 1`, finite `rate_equation` and finite `_kcat_forward`, **no
  sink** (v finite as a product → 0), and v = 0 at the equilibrium metabolite ratio —
  matched against a two-conformation mass-action ground truth in
  `test/allosteric_ground_truth.jl` (build the oracle as in the ping-pong gate added on
  `mwc-targeted-fixes`; the inactive conformation is catalytically dead, so it contributes
  only `L` to the denominator).
- **`_expand_change_allo_state` gate.** Relaxing a chemical step of a dead-inactive
  `:OnlyA`-binding mechanism is dropped (produces no partial); relaxing the `:OnlyA`
  binding to `:NonequalAI` is retained.
- **Regression.** Single-chemical-step allosteric mechanisms unchanged; the 12+
  `allosteric_ground_truth.jl` gates green; the allosteric golden byte-identical; the
  `rate_equation` perf gate (0 alloc, <120 ns) green; enumeration-count tests updated to
  measured values with the rule as the explanation.

## Risks

- **Enumeration coverage.** Collapsing four completions to one removes the opposing-binding
  hypotheses (balanced K-system). Those were physically broken (inactive cannot bind the
  essential substrate), so this is intended — but confirm no *legitimate* distinct
  hypothesis is lost. Multiple `:OnlyA` bindings remain reachable by iterating the promote
  move; verify the iterated path still reaches the multi-`:OnlyA` mechanisms the search
  needs.
- **No constructor backstop (decision 3).** `_expand_change_allo_state` is now the only
  place a relaxation could recreate a partial, so its filter is load-bearing — a missed
  case there re-introduces an enumerated sink with no error. The `_expand_change_allo_state`
  gate above is the guard against that.
- **Multi-cycle mechanisms.** "All chemical steps in the mechanism" would over-restrict a
  mechanism with two independent catalytic cycles where one should stay active. The plan
  must confirm the enumeration produces no such mechanism (ping-pong is a single linear
  cycle); if one exists, revisit decision 2.
- **V-system partials remain** (flagged above) — the ping-pong V-system crash class is not
  addressed by this spec and needs its own measurement and decision.

## Evidence

Scratchpad scripts: `completions.jl` (the four minimal completions; dead-inactive absent),
`err2_thermo.jl` (the four-variant sink/consistency table), `err2_scope.jl` (the
inactive-graph trap trace), `prereq.jl` (dead-inactive `d_free_I = 1` via stranding),
`a_form.jl` (the active-form cycle).
