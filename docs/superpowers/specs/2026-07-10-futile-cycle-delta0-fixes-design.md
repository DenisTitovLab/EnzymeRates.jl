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
fall into three classes, and each fix below owns exactly one class:

| Class | # edges | Move | What the child is | Root cause | Fix |
|-------|---------|------|-------------------|------------|-----|
| **A**  | 5  | `split` | same eq_hash as parent | non-idempotent canonicalization | Fix 1 |
| **B2** | 8  | `split` | **same rate function**, different eq_hash | parent is over-parameterized | Fix 2 |
| **C**  | 20 | `change_allo_state` | genuinely different function, no new fitted param | pointless (Wegscheider-reverted) relaxation | Fix 3 |

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
  three redundant directions (from the Jacobian null space, identical across all six parents) are:
  1. `kon_A_Pyruvateinh_E` + `koff_A_Pyruvateinh_E` move together with `v` unchanged → only the ratio
     `K` is identifiable — the SS binding's *speed* is redundant;
  2. `kon_A_Pyruvate_EPyruvateinh` + `koff_A_Pyruvate_EPyruvateinh` → same;
  3. **`L`** — the allosteric constant, but *entangled with* the two dead-end-SS speeds above, not a
     separate cause. Converting the dead-end bindings to rapid-equilibrium (Fix 2) resolves all three
     directions at once (verified: `np 10→8`, `rank 7→8`, fully identifiable).

  The two SS-speed directions trace to an **enumeration-intent bug**. `_expand_add_dead_end_regulator`
  creates competitive-inhibitor bindings as RE and creates each catalytic mirror step (the same step in
  an inhibitor-bound species) preserving its type, in the same kinetic group. But `_expand_re_to_ss`
  flips *any* all-RE kinetic group to SS with no guard, so it (i) makes an inhibitor binding
  steady-state — `g5` (`E→E·Pyruvateinh`, `INHIBITOR`, SS), while the sibling `g3` inhibitor binding
  stays RE, showing it's `re_to_ss`, not systematic — and (ii) once `split` has separated a catalytic
  step's inhibitor-bound mirror into its own group, flips that mirror independently of its
  inhibitor-free base: Pyruvate binds RE in the inhibitor-free branch (`g4`/`g7`) but SS in the
  inhibitor-bound branch (`g9`, `E·Pyruvateinh→E·Pyruvate·Pyruvateinh`). Both are dead-end SS bindings
  whose speed is never identifiable — the two null directions. So the over-parameterization is a
  *symptom* of the enumeration emitting mechanisms it never intended.
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
`fitted_params` (the parent's fitted set contains `L`, not `K_I_Pyruvate_E`). So the box is enforced,
and the residual `L` non-identifiability is entangled with the dead-end-SS speeds (resolved by Fix 2),
not a missed box.

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

## Fix 2 — Enumeration: inhibitor bindings are RE-only, and catalytic-step mirrors don't diverge (class B2)

The fix is at the enumeration move that broke the intent, `_expand_re_to_ss`, not in the derivation.
Two invariants:

1. **Competitive-inhibitor bindings are RE-only.** A binding onto a dead-end inhibitor complex has
   only its dissociation constant identifiable, and it is always a dead-end, so `_expand_re_to_ss`
   must never flip an inhibitor binding to steady-state. (This alone fixes `g5`.)
2. **A catalytic step and its inhibitor-bound mirror share RE/SS type — they never flip
   independently.** The same-group case is already clean (`re_to_ss` flips the whole group at once).
   When a `split` has separated the mirror into a different kinetic group, and the mirror's
   inhibitor-free counterpart is a genuine catalytic-cycle step (**not** a dead-end species carrying
   both a substrate and a product), the two flip **together**. The mirror lookup borrows
   `_expand_add_dead_end_regulator`'s machinery: map a species to its inhibitor-added / inhibitor-removed
   counterpart and find the step between them (the inverse of the `de_species_map` / mirror-step
   construction that move already does).

The type-lock is on the **step type**, not identifiability. Biochemically the binding mechanism is the
same whether or not an inhibitor is bound elsewhere, so a productive base legitimately flipped to SS
carries its dead-end mirror to SS too. Forcing the mirror to stay RE would fuse inhibitor-bound and
inhibitor-free forms into one RE group — biochemically wrong. The resulting "both-SS" mechanism may
not be fully identifiable, but it is parameter-dominated (cannot win parsimony) and is **not** a
delta-0 cycle edge, so it is accepted.

This closes the B2 delta-0 edges: for the reproducers the base binding is RE, so the mirror follows to
RE, giving the both-RE mechanism — verified fully identifiable (`np 10→8`, `rank=8`), with no label
degeneracy left for a split to duplicate. It is also a **correctness fix**: the enumeration was
emitting type-inconsistent, never-identifiable mechanisms regardless of the futile cycle.

**Test**: `_expand_re_to_ss` never yields a steady-state inhibitor binding; a catalytic step and its
inhibitor-bound mirror always share type after any `re_to_ss` (no divergence); the six reproducer
parents flip to the fully-identifiable both-RE form (`np 10→8`, `rank=8`); equilibrium-flux oracle
(`v = 0` at `Q = Keq`) still holds.

## Fix 3 — change_allo_state delta-0 filter (class C)

Analogous to the split self-loop guard: `_expand_change_allo_state` drops a relaxation that does not
increase the fitted-parameter count — i.e., the freed inactive-state parameter is reverted (made
dependent) by Wegscheider. A relaxation that frees no identifiable parameter is not a meaningful
refinement; it produces a same-complexity variant that can never win parsimony.

This is independent of B2: the class-C parents are fully identifiable, and their change_allo delta-0
children are genuinely-different functions (not over-parameterization). Fix 3 removes those 20 edges;
Fix 2 owns B2; there is no cross-dependency.

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
the existing split no-op guard. (Note the relaxation must genuinely free no parameter: flipping a
legitimate `:OnlyA` group all the way to `:EqualAI` *removes* an identifiable direction, so the
filter keys on "no fitted-param increase from this specific relaxation", not on de-allostering.)

**Test**: the reproducer change_allo_state delta-0 child is dropped; a relaxation that *does* free an
identifiable parameter (delta ≥ 1) is kept; no legitimate delta≥1 allosteric mechanism disappears
from `MECHANISM_TEST_SPECS`.

## Verification (whole change)

- **B2 attribution — confirmed.** Converting the two dead-end SS bindings of the B2 reproducer to
  rapid-equilibrium: `np=10, rank=7` → `np=8, rank=8` (fully identifiable). Fix 2 alone resolves all
  three redundant directions; Fix 3 alone leaves `rank=7`. So Fix 2 owns B2 entirely.
- The 33 reproducer delta-0 edges are all removed: re-run the systematic sweep (`sweep.jl`) and
  confirm 0 delta-0 edges remain across the 26 parents.
- Full `Pkg.test()` green, including the allosteric golden re-baseline, the `rate_equation`
  0-allocation / sub-120 ns perf gate, and the parameter-naming chokepoint guard.
- Re-run the LDH four-inhibitor search (or a bounded local proxy) and confirm it terminates by
  draining the frontier rather than hitting wall-clock (the beam-loss table was already final by
  iteration 7).

## Sequencing

Fix 1, Fix 2, Fix 3 are independent and can land as separate commits. Recommended order: Fix 1
(smallest, also improves dedup broadly), then Fix 3 (removes the largest edge class), then Fix 2 (the
`_expand_re_to_ss` invariants, with a golden re-baseline since it changes which mechanisms the move
emits). Re-run the sweep after each to measure the residual edge count.

## Open questions carried into planning

1. **Fix 2 mirror predicate** — pin down "inhibitor-free counterpart is a catalytic-cycle step, not a
   dead-end species with both a substrate and a product bound" as a concrete structural test, and how
   `_expand_re_to_ss` walks the mirror class to flip linked groups together.
2. **Fix 3 check** — param-count (compile) vs structural-Wegscheider implementation.
