# NonequalAI Validity — Rank/Nullspace Algorithm — Design Knowledge

**Date:** 2026-05-29
**Status:** Design knowledge for a future follow-up PR. Not implemented here.
**Companions:**
- `2026-05-29-direction-symmetry-constraint-resolution.md` — *how* a
  constraint is resolved (derivation mechanics; complementary to this doc).

## What this algorithm is for

Determine **which `:NonequalAI` tag configurations are non-degenerate** — i.e.
configurations in which every `:NonequalAI` group's active/inactive split is a
*genuinely free* parameter, not one the thermodynamics forces to zero. This is
a **validity / enumeration** concern, separate from the **symmetry principle**
(which decides *how* a constraint is resolved). Without this validation, a
degenerate config can still emit a phantom `K_T` fit parameter forced to equal
`K_A`.

It serves two consumers:
1. **Constructor validator** — reject (or normalize) a degenerate
   `AllostericMechanism` so the fitter never explores a phantom dimension.
2. **Enumeration filter** — propose only non-degenerate `:NonequalAI`
   extensions.

## The test

In log space, subtracting the active and inactive forms of each thermodynamic
cycle (live in both states) gives the **split constraint**
`Σ_g c_g · d_g = 0`, with `d_g = log(symbol_A,g) − log(symbol_I,g)` and `c_g`
the cycle's incidence of group `g` (rows assembled into `C_live`; cycles
broken in the inactive state by an `:OnlyA` member are dropped).

**Free-absorber column set** (the crucial correction):

```
F = {NonequalAI RE groups}  ∪  {all SS groups present in the live cycles}
```

SS groups are included **regardless of tag** because an SS step's *reverse
rate* is a dependent (derived) parameter whose split is free — it absorbs
constraints without contradicting any tag. `:EqualAI` RE groups are pinned
(`d_g = 0`) and excluded.

A `:NonequalAI` group `g` is **non-degenerate** iff its split can be nonzero in
some solution of the homogeneous system restricted to `F`:

```
valid(g)  ⟺  rank(C_live[:, F]) == rank(C_live[:, F ∖ {g}])
```

(i.e. column `g` is linearly dependent on the other free-absorber columns). A
whole configuration is valid iff every `:NonequalAI` group passes. Equivalently,
the valid `:NonequalAI` sets are the **unions of circuits** of the column
matroid of `C_live[:, F]`.

### Why `F` must include SS groups

The naive version restricted `F` to `:NonequalAI` columns only and therefore
declared PK invalid (lone NonequalAI PEP). That was wrong: the catalytic
reverse `k5r`/`k5r_T` is an SS-group absorber that makes PEP's split free.
With SS groups in `F`, PK is correctly **valid**.

## What it catches

1. **Lone NonequalAI binding K trapped in a pure-RE Wegscheider loop** — the
   canonical fixture is the existing **"Random-order Bi-Bi"** test mechanism
   (`test/mechanism_definitions_for_test_enzyme_derivation.jl`, `n_wegscheider=1`,
   `n_mirror=0`: its binding steps are ungrouped, so the square
   `E→EA→EAB→EB→E` is a genuine Wegscheider cycle). Make it allosteric with
   **one** binding group `:NonequalAI`, the rest `:EqualAI`. That group `g1`
   sits in two cycles: the catalytic Haldane (SS reverse as absorber) **and**
   the pure-RE Wegscheider square (no SS step). Over `F = {g1, SS-reverse}` the
   Wegscheider row is `[c·g1, 0]` — the SS reverse is not in the binding loop,
   so it cannot help — forcing `d_g1 = 0` → **degenerate** → reject. The SS
   reverse rescues a lone NonequalAI in a *Haldane* cycle (why PK is valid) but
   never in a pure-RE *Wegscheider* loop. A permissive implementation can
   silently promote the Wegscheider-dependent `:EqualAI` binding K to absorb
   the split (computes zero at equilibrium, over-parametrized); this rejection
   replaces that behavior. The symmetry rewrite cannot rescue it
   either (no speed DOF in a pure-RE loop). Only this rank test flags it.
2. **Full-rank multi-cycle (including Haldane-containing) mechanisms** —
   interlocking Haldane + Wegscheider cycles where the free-absorber columns are
   full-rank, so some `:NonequalAI` split is forced to zero *even with SS
   absorbers present*. Example shape: three groups with cycle basis
   `[[1,1,0],[1,0,1],[0,1,1]]` (rank 3 over 3 splits) → every split pinned →
   the "all three NonequalAI" config is degenerate though a naive per-cycle
   count rule would accept it.

The rank test is **basis-invariant** (depends only on the cycle space); a
per-cycle "≥2 NonequalAI" count rule is basis-dependent and strictly more
permissive (accepts degenerate configs). Use the rank test, not the count.

## Concrete cycle machinery

`_thermodynamic_cycles(m::Mechanism) → (C, b)`: integer nullspace of the
enzyme-form incidence matrix, computed **directly on the `Mechanism` struct**
(no `Sig` lift / no `@generated` derivation), following the
`_n_fit_params_estimate(m::Mechanism)` precedent. `b` classifies Haldane vs
Wegscheider. This is shared with the symmetry-resolution implementation — the
expensive piece is building the cycle basis; the rank test on top is cheap.
The cycle basis depends only on the catalytic skeleton (not the tags), so it is
cacheable per skeleton across the tag-variants the enumerator spins out.

Regulators are **out of `C_live` by construction**: the basis is built from the
catalytic step graph only; regulator-site thermodynamics live in the MWC
rate-equation assembly, and dead-end inhibitors are pendant branches that add no
cycle.

## Enumeration

Replace per-group `:NonequalAI` flips with **minimal +1-parameter extensions**:
from the current `:NonequalAI` set `N`, emit each minimal superset that raises
the free-split count by exactly 1 and stays rank-valid.
- Single absorbing cycle ⇒ the minimal step is a pair (or a lone NonequalAI
  when an SS reverse absorbs it — e.g. PK-style, +1).
- Coupled multi-loop ⇒ the smallest relative circuit.
- The ladder (+2, +3, …) is climbed over successive enumeration rounds.

**Parameter count:** the `:NonequalAI` contribution is the **split-space
dimension** `|N| − rank(C_live[:, N])` (per live-cycle accounting with SS
absorbers), replacing the incorrect "+1 per NonequalAI group."

## Placement & performance

- Shared predicate `_nonequalai_splits_free(am) → Bool`; constructor throws on
  `false`, enumeration pre-filters with it (so the throw is a backstop, not a
  hot path).
- No caching initially; benchmark a full enumeration and add skeleton-keyed
  memoization of the cycle basis only if measurably needed (YAGNI).
- The validator error message names the degenerate group(s) and the remedy
  (tag `:EqualAI`, or give it a co-cyclic absorber / NonequalAI partner).

## Tests

Anchor on real mechanisms with multi-cycle structure (random-order bi-bi, PK,
random-order ter, pyruvate dehydrogenase/carboxylase); derive expected
`C_live`, ranks, and valid/degenerate `:NonequalAI` sets with **independent**
linear algebra (generic `nullspace` / fraction-free elimination), plus by-hand
cross-checks for the small mechanisms. Where a config is degenerate, also assert
the physical symptom (a phantom parameter that does not change the fitted rate
curve, or — if mis-resolved — a nonzero rate at chemical equilibrium).

## Implementation Target

Implement as a future feature alongside or after the symmetry-principle
rewrite. Both features share `_thermodynamic_cycles`.
