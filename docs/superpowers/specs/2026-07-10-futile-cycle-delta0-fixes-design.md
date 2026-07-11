# Futile-cycle delta-0 edges: canonicalization fix, dead-end-SS reduction, change_allo filter

## Context

The LDH HPC run with four competitive inhibitors (`docs/ldh_hpc_results/2026_07_09_results_2`)
finishes the useful search by iteration 7 but then spins in a **period-6 limit cycle** to
wall-clock: from iteration 13 on it produces 0 new fits, re-expanding a fixed set of ~755
equations. The complexity filter (#63) stopped the earlier segfault; this is the remaining
non-termination.

The cycle is sustained by **delta-0 expansion edges** — expansion moves that produce a child with
the same fitted-parameter count as its parent, so the beam circulates at fixed counts 10–13 instead
of climbing past the cap and draining. A systematic sweep of every delta-0 edge across the 26
low-parameter reproducer parents (33 distinct edges) shows they are **not one phenomenon**. They
fall into three classes with three distinct root causes:

| Class | # edges | Move | What the child is | Root cause |
|-------|---------|------|-------------------|------------|
| **A**  | 5  | `split` | same eq_hash as parent | non-idempotent canonicalization |
| **B2** | 8  | `split` | **same rate function**, different eq_hash | parent is over-parameterized |
| **C**  | 20 | `change_allo_state` | genuinely different function, no new fitted param | pointless (Wegscheider-reverted) relaxation |

Only **6 of the 26 parents** are over-parameterized (all `np=10`, identifiable rank 7); the other 20
are fully identifiable. So the fixes must be targeted, not uniform.

### Evidence (all reproduced under c201d50; scripts in the session scratchpad)

- **A**: `_canonical_mechanism` is not idempotent. For a split no-op the raw split canonicalizes to 9
  groups (`C(raw) != C(parent)`, so the split guard keeps it), but a **second** pass merges back to
  the parent's 8 groups. Cause: the merge key for a `:NonequalAI` group is
  `(:NonequalAI, A-state-fold, I-state-fold)`; on pass 1 `rename_I` is empty, so two groups match on
  the A-fold but differ on the I-fold and don't merge — merging the A-state tie then changes the
  I-graph so pass 2's I-state solve finds the tie and merges. The A/I Wegscheider maps are coupled.
- **B2**: the split child and its parent compute the **identical** rate function
  (`max |v_parent − v_child| / |v| = 8.4e-16` over 5000 wide-range points, identical fitted-param
  names) yet carry different eq_hashes — `_rate_eq_dedup_key` hashes rendered *text*, so it cannot
  see the equivalence. This is only possible because the parent is over-parameterized: its 10 fitted
  params have identifiable rank **7** (a clean 7-order singular-value cliff, robust to sampling). The
  three redundant directions (from the Jacobian null space, always the same six parents) are:
  1. `kon_A_Pyruvateinh_E` + `koff_A_Pyruvateinh_E` move together with `v` unchanged → only the ratio
     `K` is identifiable — the SS binding's *speed* is redundant;
  2. `kon_A_Pyruvate_EPyruvateinh` + `koff_A_Pyruvate_EPyruvateinh` → same;
  3. **`L`** — the allosteric constant is non-identifiable because the mechanism's `:NonequalAI`
     Pyruvate affinities are all pinned equal to their active-state values, so the allostery is
     **degenerate**: relaxing those groups added no identifiable freedom. This is the same "pointless
     relaxation" as class C, baked into the mechanism (see the trace note below).
- **C**: the `change_allo_state` delta-0 child is a **genuinely different** function, not a
  parent-duplicate — cross-fitting confirms the child cannot reproduce the parent's data
  (residual 1.3, not 0). The relaxation changed the mechanism but added no identifiable parameter
  (the freed parameter is Wegscheider-reverted). These children are correctly eq_hash-deduped in the
  run (they return "inherited"), but they are *structurally distinct*, so with no seen-set they are
  re-expanded forever.

### Trace note — the inactive-state Wegscheider box is NOT missed (2026-07-10)

An early hypothesis was that the derivation misses the inactive-state box that pins `K_I_Pyruvate_E`.
Tracing `_combined_state_dependent_exprs` on the B2 reproducer disproves it: the box cycle is
enumerated in both states (A raw `cycle3`, I raw `cycle2`), and the inactive state already handles it
— `_state_wegscheider_rename_map(am, :I)` returns `Dict(:K_I_Pyruvate_E => :K_I_Pyruvate_ENAD)`,
folding the two inactive Pyruvate affinities onto one representative (the box row is then a satisfied
`0 = 0`), and the phantom filter in `_dependent_param_exprs` drops the folded `K_I_Pyruvate_E` from
`fitted_params` (the parent's fitted set contains `L`, not `K_I_Pyruvate_E`). So the box is enforced.
The residual `L` non-identifiability is degenerate allostery, addressed by Fix 3, not a missed box.

## Goals

Stop the futile cycle by removing all three delta-0-edge classes at their source, and along the way
fix a correctness bug (over-parameterized mechanisms report inflated parameter counts, which can
distort parsimony-based model selection). No search-level stop rule or seen-set is introduced; each
fix removes edges at the move/derivation that creates them.

## Non-goals

- Text-to-function dedup (hashing evaluated rate values). Removing the over-parameterization makes
  equivalent mechanisms parameterize identically, so the existing textual eq_hash suffices.
- A global structural seen-set / blanket delta≥1 monotonicity guard. The three targeted fixes remove
  the edges at the move, so the beam's monotonic-drain termination is restored without a guard.

## Fix 1 — Idempotent canonicalization (class A)

`_canonical_mechanism` must reach a fixed point. Iterate `_merge_tied_kinetic_groups` until the
partition stops changing (with a small max-iteration safety bound; convergence is 2 passes in every
case observed). This makes the split self-loop guard `child == _canonical_mechanism(parent)` catch
the delta-0 no-op splits it currently misses, and it collapses renaming variants more broadly (the
same non-idempotency inflates the structural-multiplicity count across the whole search).

**Design decisions**
- Iterate inside `_canonical_mechanism` (both the `Mechanism` and `AllostericMechanism` methods)
  rather than inside `_merge_tied_kinetic_groups`, so every caller benefits and the merge stays
  single-pass and testable.
- Re-baseline: some mechanisms' canonical form changes, so allosteric golden fixtures and eq_hashes
  need a reviewed re-baseline.

**Test**: `_canonical_mechanism(_canonical_mechanism(m)) == _canonical_mechanism(m)` for all
`MECHANISM_TEST_SPECS`; the reproducer A split no-op is dropped by `_expand_split_kinetic_group`.

## Fix 2 — A dead-end steady-state binding carries only its equilibrium constant `K`

A binding to a complex with no onward catalytic edge (a competitive-inhibitor / dead-end leaf) has,
at steady state, `[EI] = [E][I]/K` — only `K = koff/kon` is identifiable, not `kon` and `koff`
separately. The derivation currently fits both, which is exactly the two SS-speed redundant
directions in every over-parameterized parent. Such a binding must contribute one parameter (`K`,
the rapid-equilibrium form), not two.

**Design decision (resolve in planning): where to enforce it.**
- *Enumeration filter* (recommended) — do not generate the steady-state form for a dead-end binding;
  it is always identifiability-equivalent to, and parameter-dominated by, the rapid-equilibrium form.
  Likely a guard in `_expand_re_to_ss` (skip dead-end bindings) and wherever dead-end regulator
  bindings are emitted steady-state. Simpler, and it also cuts the wasted enumeration.
- *Derivation reduction* — recognize in the constraint solve that a dead-end binding's `kon`/`koff`
  appear only as their ratio and collapse them. More general but harder.

To confirm before implementing: that both redundant SS bindings are structural dead-ends (a form
whose only edges are the binding and its reverse).

**Test**: identifiable rank rises by 2 for the six reproducer parents; equilibrium-flux oracle
(`v = 0` at `Q = Keq`) still holds.

## Fix 3 — change_allo_state delta-0 filter (class C, and the `L`-degenerate part of B2)

Analogous to the split self-loop guard: `_expand_change_allo_state` drops a relaxation that does not
increase the fitted-parameter count — i.e., the freed inactive-state parameter is reverted (made
dependent) by Wegscheider. A relaxation that frees no identifiable parameter is not a meaningful
refinement; it produces a same-complexity variant that can never win parsimony.

This removes class C (20 edges) directly. It also removes the **third** over-parameterization
direction in B2: a mechanism only acquires a degenerate `:NonequalAI` tag (which makes `L`
non-identifiable) through a delta-0 relaxation, so filtering those relaxations prevents the
over-parameterized B2 parents from being generated at all. Fix 2 handles the remaining two
(dead-end-SS) directions.

**Design decision (resolve in planning): the delta check.**
- *Compile + `fitted_params`* per candidate child, compare count to the parent. Simple and exactly
  correct; cost is one derivation per candidate (~7 candidates per parent).
- *Structural Wegscheider check* — decide from the per-state rename map whether the freed parameter
  is dependent, without a full compile. Cheaper, matches the split guard's spirit, but must be shown
  equivalent to the param-count test.

Recommend starting with the param-count test for correctness, then optimizing to the structural check
if expansion throughput requires it.

**Policy note**: these delta-0 relaxations are genuinely-different, often fully-identifiable models,
so dropping them forgoes some same-count allosteric variants. This is intended: they add no
identifiable degree of freedom over the parent and cannot be selected over it — the same rationale as
the existing split no-op guard.

**Test**: the reproducer change_allo_state delta-0 child is dropped; a relaxation that *does* free an
identifiable parameter (delta ≥ 1) is kept; no legitimate delta≥1 allosteric mechanism disappears
from `MECHANISM_TEST_SPECS`.

## Verification (whole change)

- **Confirm the B2 attribution**: after Fix 2 + Fix 3, the six reproducer parents are fully
  identifiable (rank == fitted count) — i.e., `L`'s non-identifiability really was the degenerate
  allostery that Fix 3 prevents, not a third independent cause. Run before locking the plan.
- The 33 reproducer delta-0 edges are all removed: re-run the systematic sweep (`sweep.jl`) and
  confirm 0 delta-0 edges remain across the 26 parents.
- Full `Pkg.test()` green, including the allosteric golden re-baseline, the `rate_equation`
  0-allocation / sub-120 ns perf gate, and the parameter-naming chokepoint guard.
- Re-run the LDH four-inhibitor search (or a bounded local proxy) and confirm it terminates by
  draining the frontier rather than hitting wall-clock (the beam-loss table was already final by
  iteration 7).

## Sequencing

Fix 1, Fix 2, Fix 3 are independent and can land as separate commits. Recommended order: Fix 1
(smallest, also improves dedup broadly), then Fix 3 (removes the largest edge class and the
`L`-degenerate B2 direction), then Fix 2 (the dead-end-SS reduction with its golden re-baseline).
Re-run the sweep after each to measure the residual edge count.

## Open questions carried into planning

1. **B2 attribution check** — verify Fix 2 + Fix 3 make the six parents fully identifiable (above).
2. **Fix 2 placement** — confirm the two redundant SS bindings are structural dead-ends and choose
   enumeration-filter vs derivation-reduction.
3. **Fix 3 check** — param-count (compile) vs structural-Wegscheider implementation.
