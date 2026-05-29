# EqualAI × NonequalAI Coupling — Design

**Date:** 2026-05-29
**Branch:** `refactor-to-concrete-types-instead-of-symbols`
**Status:** Approved design, pending implementation plan

## Problem

MWC allosteric mechanisms tag each catalytic kinetic group with an
allosteric state: `:OnlyA`, `:EqualAI`, or `:NonequalAI` (`:OnlyI` is
rejected for catalytic groups). The enumeration's `_expand_change_allo_state`
flips one group to `:NonequalAI` at a time, which manufactures *mixed*
configurations such as `[:EqualAI, :NonequalAI, :EqualAI]`.

These mixed configurations are **physically inconsistent**. The active (A)
and inactive (I) states of an MWC enzyme share the same enzyme-form graph
and the same steps — they differ only in parameter values. Each
thermodynamic cycle (Haldane or Wegscheider) therefore imposes the *same*
relation in both states. When a single group inside a live cycle is
`:NonequalAI` while the rest are `:EqualAI`, the two-state relation forces
the shared (`:EqualAI`) symbols to take different values in the A and I
branches — a contradiction. Symptom on the current branch: `m_mixed`
evaluates to rate ≈ 2.615 at chemical equilibrium (should be 0), and PK
reports the wrong Haldane/mirror constraint counts.

## The invariant (correctness condition)

Work in log space. For a thermodynamic cycle with incidence row `c` over
the kinetic groups:

```
Σ_g c_g · log(symbol_A,g) = b · log(Keq)      (active state)
Σ_g c_g · log(symbol_I,g) = b · log(Keq)      (inactive state, if live)
```

- `g` indexes a kinetic group; `symbol_*,g` is its per-state equilibrium
  symbol (an RE group's Kd, or an SS group's `kf/kr`).
- `c_g` is the cycle's traversal incidence of group `g`'s step (±, from the
  integer nullspace of the enzyme-form incidence matrix). Non-members: 0.
- `b` is how many net overall reactions one trip performs: `b ≠ 0` ⇒ Haldane,
  `b = 0` ⇒ Wegscheider.

Subtracting the two equations (RHS cancels) and defining the **split**
`d_g = log(symbol_A,g) − log(symbol_I,g)`:

```
Σ_g c_g · d_g = 0      for every cycle live in both states.
```

By tag: `:EqualAI` ⇒ `d_g = 0` (one shared symbol); `:NonequalAI` ⇒ `d_g`
is meant to be a free, usable parameter; `:OnlyA` member ⇒ the group has no
I-symbol, so any cycle containing it has **no I-constraint** — that cycle
is *dead in I* and its row is dropped.

Stacking the live cycles into `C_live` and restricting to the `:NonequalAI`
columns `N` (EqualAI columns contribute `d = 0`):

```
C_live[:, N] · d_N = 0
```

**A configuration is valid iff every NonequalAI group's split is
exercisable** — i.e. no `d_g` is forced to zero. Group `g` is forced to
zero iff its column is linearly independent of the other NonequalAI
columns, so the test is:

```
valid  ⟺  ∀ g ∈ N:  rank(C_live[:, N]) == rank(C_live[:, N ∖ {g}])
```

Equivalently: the valid NonequalAI sets are exactly the **unions of
circuits** of the column matroid of `C_live`.

### Why the rank test, not a per-cycle count

"Each live cycle has 0 or ≥2 NonequalAI members" is a *necessary but not
sufficient* shadow of the rank test:

- **rank-pass ⟹ count-pass**, but **count-pass ⇏ rank-pass**. The count
  rule wrongly *accepts* some degenerate configs (it never wrongly rejects).
- It is **basis-dependent** — the answer depends on which basis of the
  cycle space the nullspace algorithm returns; the rank test depends only
  on the (physical) cycle space.
- **Minimal divergence:** three NonequalAI groups `a,b,c` with cycle basis
  `[[1,1,0],[1,0,1],[0,1,1]]` (full rank 3 over 3 splits) — every cycle has
  2 NonequalAI members (count says valid) but all splits are forced to zero
  (rank says degenerate, invalid).

The rank test and the count rule run on the *same* `C_live`; the cost is
building `C_live`, not the rank-vs-count step on top. So the rank test is
both strictly more correct and effectively free — chosen on principle.

### Consequences

- **No structural shortcut is safe.** Even "flip the whole coupling
  cluster to NonequalAI" can be invalid (the full-rank triple above), so
  the rank predicate is the single source of truth; enumeration proposes,
  the predicate disposes.
- **OnlyA rescues a lone NonequalAI** automatically — its cycle row is
  dropped from `C_live`, so a single NonequalAI in an otherwise-dead cycle
  is unconstrained and valid.

## Design

### Component 1 — concrete cycle basis
`_thermodynamic_cycles(m::Mechanism) → (C, b)`: integer nullspace of the
enzyme-form incidence matrix, computed directly on `Mechanism.steps` (no
`Sig` lift / no `@generated` derivation), following the
`_n_fit_params_estimate(m::Mechanism)` precedent of computing cycle info on
the concrete struct. `b` carries the Haldane/Wegscheider classification via
the existing `classify_cycle` logic ported off the singleton accessors.
This is the one shared computation.

### Component 2 — `C_live` + rank predicate
- `C_live`: rows of `C` with no `:OnlyA` group in their support (OnlyA ⇒
  cycle dead in I ⇒ drop the row). After dropping, OnlyA columns vanish.
- `_nonequalai_splits_free(am::AllostericMechanism) → Bool`: build `C_live`,
  restrict to NonequalAI columns, apply the rank test.
- **Regulator orthogonality is asserted, not assumed.** Regulator Kreg's
  do not appear in catalytic cycles, and dead-end binding is a pendant
  branch that adds no cycle. A guard errors if a regulator symbol ever
  appears in `C`, so the assumption cannot silently rot. Regulator-ligand
  allo-states stay per-ligand.

### Component 3 — validator placement
The `AllostericMechanism` constructor calls `_nonequalai_splits_free` and
throws a clear message on `false` (DSL / direct-construction safety).
Enumeration **pre-filters** with the same predicate, so it never constructs
an invalid config — the constructor throw is a backstop, not a hot path.
**No caching initially.** Benchmark a full enumeration; add skeleton-keyed
memoization of the cycle basis only if the redundant per-construction
compute is measurably too slow (YAGNI).

Error message names the offending groups and the two remedies (make the
lone group EqualAI, or give it a co-cyclic NonequalAI partner).

### Component 4 — enumeration move
Replace `_expand_change_allo_state`'s per-group NonequalAI flip with
**minimal +1-parameter extensions**: from the current NonequalAI set `N`,
emit each minimal superset that raises the free-split count
(`|N| − rank(C_live[:,N])`) by exactly 1 and stays rank-valid.

- Single-loop mechanisms ⇒ the minimal first step is any **pair** (e.g. a
  4-group loop yields the 6 pairs, each +1 param).
- Coupled multi-loop mechanisms ⇒ the minimal step is the smallest relative
  circuit (e.g. a shared-group triple that adds +1).
- The ladder (+2, +3, …) is climbed over successive expansion rounds, as
  the beam search re-feeds frontier mechanisms — same pattern as RE→SS and
  split-group.

All-EqualAI base and per-group `:OnlyA` flips in `_expand_to_allosteric`
are unchanged. Regulator-ligand flips unchanged.

### Component 5 — parameter count (Piece 1 only)
The allosteric NonequalAI contribution in `_n_fit_params_estimate` becomes
the **split-space dimension** `|N| − rank(C_live[:,N])`, replacing the
incorrect "+1 per NonequalAI group" (which over-counts: a 4-group loop
fully NonequalAI adds +3, not +4; a shared-group triple adds +1, not +3).

**Deferred (Piece 2, separate follow-up):** making the *base* catalytic
`_n_fit_params_estimate(m::Mechanism)` exact via `rank(A_merged)` (so it
equals `length(fitted_params(compile_mechanism(m)))`). It is cheap in
compute but changes the function's contract from upper-bound to exact,
making the call-site `n_subs+n_prods+1` floor dead and flipping many
`length(fitted_params) <= estimate` test assertions to `==`. Kept out of
this correctness fix to keep the diff surgical and the hand-back clean.

### Component 6 — tests
- **New:** constructor rejects a mixed config (e.g. `[:EqualAI, :NonequalAI,
  :EqualAI]` on a single loop) with the documented message.
- **New:** rank-divergence coverage — a mechanism where the count rule and
  rank test disagree (full-rank shared-group triple), asserting the
  validator rejects it.
- **`m_mixed`:** rewrite as a Wegscheider-consistent config (all-NonequalAI
  on its loop) so the rate-eq math is exercised and equilibrium → 0, **and**
  add a separate `@test_throws` for the rejected mixed config. Keeps both
  paths covered.
- **PK Constraints golden:** recompute `n_haldane_constraints` /
  `n_mirror_constraints` from the dep machinery on the corrected PK config
  and match the implementation to truth — no blind golden edits.
- **Enumeration variant-count tests:** update to the rank-filtered move
  counts.
- **No deletions**, no `@testset` removal, no `MECHANISM_TEST_SPECS` removal.

### Component 7 — docs
CLAUDE.md "Allosteric state taxonomy" gains a paragraph on the coupling
rule; the parent structural-parameter-names spec gets a postscript noting
the constraint.

## Out of scope / tripwires
- Do **not** touch the structural-naming chokepoint (`_state_tag`,
  `_render_binding`, `_render_iso`) — parent session's work.
- Do **not** change `_flip_to_inactive`'s semantics.
- Leave the deferred `canonical-hash partition stability` test (21 vs 23)
  failing — parent session's territory.
- Do not modify the `positional_params` shim.

## Done when
- Full suite green except the deferred hash-partition failure.
- Validator rejects mixed configs with the documented message.
- Enumeration never produces an invalid mechanism.
- CLAUDE.md + parent spec updated.
