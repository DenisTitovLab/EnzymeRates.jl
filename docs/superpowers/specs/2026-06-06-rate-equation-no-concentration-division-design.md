# Representation-Correct, Division-Free Rate-Equation Derivation — Design

**Date:** 2026-06-06 (rev. 2026-06-07 after review)
**Status:** design proposed, pending spec review → implementation plan

---

## Goal

A mechanism's rate equation must depend only on *what the mechanism is*, never on how its
steps are written, and must be finite wherever a real measurement is defined. Several
coupled bugs break this; this PR fixes them together (they must ship together because the
mutating dedup makes *every* equation fail, not just the ping-pong ones).

Four parts:

1. **Enumeration represents ping-pong with residuals, not conformations**, and only emits
   **atom-conserving** steps.
2. **Derivation references each rapid-equilibrium segment to its free enzyme**, and reduces
   the result to lowest terms over concentrations (**minimal concentration-GCD**) so no
   concentration ever sits in a denominator.
3. **Mechanisms are canonical by construction**; dedup stops mutating; derivation no longer
   depends on step index. Textbook oracles are bridged so they need no rewrite.
4. **Regression tests** lock in atom conservation, division-freeness at zero concentration,
   and representation-independence.

*Simplification:* the PR removes overengineering debt — the `is_estar` conformation plumbing
(Part 1) and the mutating `_canonicalize_mechanism!` (Part 3; canonicalization moves into the
constructors, dedup becomes `unique!`). `normalize` was investigated and is **kept** — it is
load-bearing for readability (Part 2).

## Evidence (LDH: `substrates NADH, Pyruvate; products Lactate, NAD; oligomeric_state 4`)

- Deduped `init` set (what `identify_rate_equation` fits): **69/69** equations divide by a
  concentration. Reverse initial velocity (`NADH=Pyruvate=0`) makes **65/69** non-finite.
- The division is *not* intrinsic: it comes from the derivation referencing each RE segment
  to `group[1]` (after `_dedup_flat!` sorts steps, a bound complex). **Free-enzyme reference
  cleans 55/69** and yields the textbook form, e.g.
  `den = 1 + NADH/K_NADH_E + NAD/K_NAD_E + Lactate·NADH/(K_Lactate_ENAD·K_NADH_E) + …`.
- The residual 14 are a single RE group with **two zero-bound forms** (`E` and `Estar`),
  both with **empty residual** — i.e. enumeration tagged a second *conformation* with no
  covalent residue. They are also **atom-non-conserving**: e.g. the step
  `ENADH → EstarLactate` turns bound NADH into bound Lactate with nothing to balance it
  (would need a chemically-absurd `C18H23N7O11P2` residual). 0/69 forms carry any residual.
- A correct ping-pong (Segel Ping-Pong Bi-Bi fixture) derives **division-free**; with
  residuals it is a valid minimum-parameter `init` mechanism (RE steps), still a single RE
  group, so it still needs the concentration-GCD step.

---

## Part 1 — Enumeration: residuals (not conformations) + atom-conserving steps

**Today (bug):** `_make_species` (`src/mechanism_enumeration.jl:13-20`) builds
`Species(mets, is_estar ? :Estar : :E)` with the 2-arg constructor — an **empty residual** —
and the ping-pong-isomerize branch (`backtrack!` Option 2, ~lines 333-379) passes
`is_estar=true`. So enumeration encodes ping-pong as a bare `:Estar` *conformation* with no
residual, and emits atom-non-conserving steps. (The residual atoms are tracked in the
`residual_atoms` bookkeeping dict but never stored on the form.) This is the documented
"empty-residual ping-pong" (CLAUDE.md constraint C4), now being reversed.

**Change:**

- **Residuals, not conformations.** `_make_species` always uses conformation `:E` and
  stores the covalent residue as a `Residual`. The residual at any intermediate is, at the
  metabolite level, `Residual(added = consumed_subs, subtracted = released_prods)` — both
  already threaded through `backtrack!` / `_release_products!`. The `:Estar` tag and the
  `is_estar` plumbing are removed.
- **Conformations stay a `Species` field** and derivation keeps supporting them, but
  enumeration generates **no** conformation-changing steps in this PR. (Adding conformation
  steps — with rules for what must occur between conformations — is a separate future PR.)
- **Atom-conservation invariant in the `Step` constructor** (`src/types.jl`). Define a
  species' *units* as the signed metabolite multiset `bound ∪ residual.added −
  residual.subtracted`. After the constructor's canonical-direction normalization, a valid
  Step must satisfy:
  - binding step: `to.units − from.units == {bound_metabolite}`
  - iso step (`bound_metabolite === nothing`): `to.units − from.units == ∅`
  This is cheap (metabolite counting, no reaction-atom lookup) and errors on violation, so
  the broken `ENADH → EstarLactate` (`to−from = {Lactate}−{NADH} ≠ ∅`) can never be built.
- **Admissible residuals (ping-pong well-formedness).** A residual is generated only if it is
  both **produced** by a step (a substrate consumed → ≥1 product released, residual formed)
  **and consumed** by a step (residual + a substrate → ≥1 product released, optionally a new
  residual). This keeps every residual a genuine, fully-cycled covalent intermediate — never
  dangling — and is what bounds which ping-pong topologies enumeration emits.

**Effect:** the previously-broken ping-pong forms become valid, atom-conserving,
residual-bearing forms (conformation `:E`), still RE / `init` / single-RE-group. The
exact `init` mechanism counts may change; **we recompute them empirically during
implementation** (they may stay near 69 if the 14 become valid residual ping-pong, or
change — to be measured, not assumed).

## Part 2 — Derivation: free-enzyme reference + minimal concentration-GCD

**Reference choice.** `_compute_alpha` (`src/rate_eq_derivation.jl:250-310`) seeds its
per-segment BFS at `group[1]` (the first form in step order — order-dependent, and a bound
complex after dedup). Change the seed to the **free enzyme of the segment**: the form with
the fewest bound metabolites, tie-broken toward no residual (the true apoenzyme), then a
deterministic key. This is chosen by *content*, so the derivation becomes
**order-independent on its own** (Part 3 reinforces it but is no longer load-bearing for
correctness), and it yields interpretable dissociation constants (`1 + [S]/K + …`). It
cleans every single-valley mechanism (the 55) with no further work.

**Concentration-GCD.** A `G=1` ping-pong has two zero-bound forms (apo `E` and the
residual form) in one segment; referenced to apo `E`, the residual form's population is
genuinely concentration-coupled (`[E·resid]/[E] ∝ [NADH]/(K·[NAD])`), so a `1/[conc]` term
survives. As the final step of `_raw_symbolic_rate_polys`, reduce the numerator/denominator
to **lowest terms over concentrations**: for each metabolite, shift its exponent so the
minimum across all terms of `N ∪ D` is 0. No term can go negative → no concentration
denominator. New helper `_reduce_conc_lowest_terms(num, den, conc_set)` in
`src/sym_poly_for_rate_eq_derivation.jl`; reduce over **concentrations only** (never
parameters — that could drop a fitted parameter). For sequential mechanisms this is the
identity (they already have a constant term); for ping-pong it clears the coupling and
yields the standard no-constant-term ping-pong form.

**`normalize` is kept — verified load-bearing for readability.** It re-expresses the
King-Altman denominator as a sum of fractional enzyme-form populations (free enzyme = 1) —
the readable `1 + [S]/K` form. Removing it is value-correct (0 mismatches across the LDH set)
but yields high-rate-constant-power equations (e.g. `K_NADH_E^3·K_NAD_E^3·… + …`),
unreadable. With a free-enzyme reference it gives the textbook form for sequential
mechanisms; the concentration-GCD post-pass then clears the residual concentration division
that `normalize` introduces for ping-pong.

- **kcat is unaffected.** `_kcat_forward` (`:699`) reads the same polys and takes ratios at
  matching concentration patterns; a common monomial factor cancels, so kcat is invariant.
  Existing kcat/rescaling/scale-invariance tests are the gate.
- **Allosteric path:** `_allosteric_num_den_exprs` (`:1518`) assembles from the per-state
  catalytic polys (now clean) by multiply/add only, so no concentration denominator appears.
- **Note for a later PR:** once we can see the real reduced equations, we may revisit the
  GCD form for biochemical readability. Out of scope here — minimal concentration-GCD only.

## Part 3 — Canonical by construction + non-mutating dedup + oracle bridge

- **Canonicalize in the constructors.** Move the step/group sorting in
  `_canonicalize_mechanism!` (`src/mechanism_enumeration.jl:1670-1693`) into the `Mechanism`
  / `AllostericMechanism` constructors (`src/types.jl`), beside the existing iso-direction
  canonicalization. Every mechanism is then canonical at construction; two step orderings
  produce the identical struct. `_canonicalize_mechanism!` is deleted.
- **Non-mutating dedup.** `_dedup_flat!` (`:1712`) collapses to `unique!(mechs)` — pure
  comparison, no mutation.
- **Derivation independent of step index.** With Part 2's content-based reference and
  canonical construction, nothing in derivation depends on the position of a step.
- **Oracle bridge (no oracle rewrites).** Textbook oracles key parameters by flat step
  position (`positional_params`, `test/test_rate_eq_derivation.jl:113`). Keep oracles in
  as-written order; record each oracle's as-written→canonical **permutation** (recoverable
  by matching as-written steps to the canonical `steps(m)`, since `Step` compares
  structurally) and have `positional_params` apply it. The `@enzyme_mechanism` block stays
  the human-readable record of "step 1, step 2, …".

## Part 4 — Regression tests

- **Division-freeness (the headline test).** Inside `run_all_tests(spec)` over
  `MECHANISM_TEST_SPECS` (already compiles `rate_equation`, true piggyback): random positive
  params/concs; for each metabolite set just that one to 0 and assert `rate_equation` is
  **finite and nonzero** (nonzero guards the spurious-`0.0`-from-`1/Inf` regression). Covers
  allosteric fixtures too.
- **Enumeration coverage** on the existing `bi_bi_pp_rxn` fixture
  (`test/test_mechanism_enumeration.jl:102`; `A[CX], B[N] → P[C], Q[NX]` — ping-pong-capable,
  the minimal analogue of the LDH reaction that surfaced these bugs): its full `init` set,
  each metabolite zeroed → finite/nonzero, folded into the existing enumeration loop.
- **Atom conservation**: every Step in every `init`/`expand` mechanism satisfies the
  signed-metabolite-unit invariant (Part 1) — and the `Step` constructor errors otherwise.
- **Canonicalization / representation-independence**: the same mechanism written in two step
  orders constructs to the identical struct; `_dedup_flat!` does not mutate its inputs and
  still collapses duplicates.
- **Must stay green unchanged**: `test_rate_equation_performance` (allocation-free,
  sub-100 ns); kcat / rescaling / scale-invariance; Aqua/JET. The 20
  `expected_factored_num`/`expected_factored_denom` snapshots and Expr-shape/flat-string
  tests will change to the reduced forms (regenerate + eyeball); oracle tests pass via the
  Part 3 bridge.

## Non-goals

- **Conformation-changing steps** in enumeration (a future PR with its own rules); the
  `Species.conformation` field and derivation support for it stay.
- **Biochemical refinement of the GCD form** (a future PR after we see the equations).
- **Steady-state ping-pong at `init`** — SS steps are added in `expand_mechanisms`; `init`
  ping-pong is rapid-equilibrium and legitimately needs the concentration-GCD.
- No change to the parameter API, thermodynamic-constraint machinery, or fitting loop.

## Open questions / risks

- Exact `init`/`expand` mechanism counts after Part 1 — measure, then update the count
  assertions and the CLAUDE.md "verified topology counts" note.
- Threading the metabolite-level residual through `backtrack!`/`_release_products!` and
  confirming derivation handles residual-bearing forms end-to-end (the `name` chokepoint
  already renders residuals).
- `normalize` stays (verified load-bearing). Mild related debt, possible follow-up cleanup
  (out of scope unless we decide otherwise): it's a `G==1`-only special-case, so
  `sigma_num`/`sigma_den` are each computed-but-sometimes-unused and `G>1` mechanisms get a
  different (un-normalized) form; unifying the fractional-population form across all `G` would
  be cleaner. In this PR, verify the GCD post-pass clears ping-pong without disturbing the 55.

## TDD order

1. Add the division-freeness + atom-conservation regression tests (fail today).
2. Part 1: residual-bearing enumeration + `Step` atom-conservation invariant; recompute and
   update mechanism counts.
3. Part 2: free-enzyme reference + `_reduce_conc_lowest_terms`; make division tests pass;
   confirm kcat/perf/Aqua green; regenerate snapshots.
4. Part 3: canonicalize in constructors; `_dedup_flat!`→`unique!`; oracle permutation bridge;
   add representation-independence tests.
5. Confirm on the LDH `identify_rate_equation` run that reverse-direction mechanisms now fit.
