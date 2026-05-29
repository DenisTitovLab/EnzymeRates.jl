# NonequalAI Validity — Rank/Nullspace Algorithm — Design Knowledge

**Date:** 2026-05-29
**Status:** Design knowledge for a FUTURE follow-up PR (after the parent
structural-parameter-names refactor lands). Not implemented here.
**Companions:**
- `2026-05-29-equalai-nonequalai-coupling-design.md` — the *contained* fix
  that ships first (makes current configs *compute*; does **not** address
  degeneracy).
- `2026-05-29-direction-symmetry-constraint-resolution.md` — *how* a
  constraint is resolved (derivation mechanics; complementary to this doc).

## What this algorithm is for

Determine **which `:NonequalAI` tag configurations are non-degenerate** — i.e.
configurations in which every `:NonequalAI` group's active/inactive split is a
*genuinely free* parameter, not one the thermodynamics forces to zero. This is
a **validity / enumeration** concern, separate from:
- the **contained fix** (which makes configs compute but would happily emit a
  phantom `K_T` fit parameter forced to equal `K_A`), and
- the **symmetry principle** (which decides *how* a constraint is resolved).

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

### Why `F` must include SS groups (the bug in the first draft)

An earlier version restricted `F` to `:NonequalAI` columns only and therefore
declared PK invalid (lone NonequalAI PEP). That was wrong: the catalytic
reverse `k5r`/`k5r_T` is an SS-group absorber that makes PEP's split free.
With SS groups in `F`, PK is correctly **valid**.

## What it catches (and the contained fix does not)

1. **Pure-RE Wegscheider loop, lone NonequalAI** — a random-order binding loop
   with no SS step. `F` = just that one NonequalAI binding group; its column is
   independent of nothing → split forced to zero → **degenerate**. The
   contained fix would emit a phantom `K_T`. The symmetry rewrite cannot rescue
   it either (no speed DOF in a pure-RE loop). Only this test flags it.
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

## Tests (true TDD — pre-derive expected values independently)

Anchor on real mechanisms with multi-cycle structure (random-order bi-bi, PK,
random-order ter, pyruvate dehydrogenase/carboxylase); derive expected
`C_live`, ranks, and valid/degenerate `:NonequalAI` sets with **independent**
linear algebra (generic `nullspace` / fraction-free elimination), plus by-hand
cross-checks for the small mechanisms. Where a config is degenerate, also assert
the physical symptom (a phantom parameter that does not change the fitted rate
curve, or — if mis-resolved — a nonzero rate at chemical equilibrium).

## Sequencing (decided)

Follow-up PR, after the parent structural-parameter-names refactor — naturally
alongside or after the symmetry-principle rewrite (they share
`_thermodynamic_cycles`). **Not** part of the contained fix.
