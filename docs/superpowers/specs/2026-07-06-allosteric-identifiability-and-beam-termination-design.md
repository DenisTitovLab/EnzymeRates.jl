# Allosteric identifiability and beam-search termination — design

Date: 2026-07-06
Status: approved — implementation in progress

## Summary

An LDH HPC run (bi-bi, `oligomeric_state:4`, `min_beam_width=50`,
`max_param_count=13`) surfaced three defects. All three trace to a single theme:
**the allosteric machinery produces parameters and mechanisms that the data cannot
distinguish, and the beam search then fails to terminate over them.**

| # | Defect | Layer | Root cause |
|---|--------|-------|-----------|
| 1 | Beam search never terminates; rewrites identical CSVs until the disk fills | enumeration | `_expand_split_kinetic_group` emits a Wegscheider-tied no-op split that `_canonical_mechanism` (PR #59) merges back to the parent — a self-loop |
| 2 | A Haldane-dependent parameter appears in `fitted_params` (a phantom fit dimension) | derivation | `indep` is a hand-assembled union, not the complement of `dep` |
| 3 | The enumerator emits allosteric mechanisms that are empirically identical to simpler ones | enumeration | `_expand_to_allosteric` emits degenerate tag patterns |

Defect 1 has its own root cause — a self-loop in the split/canonicalize interaction
introduced by PR #59 — but the wide beam is what exposes it. A narrow or quickly-draining
beam (v0.1.4, or any small search) hits the `max_param_count` cap before the self-loop can
sustain a fixpoint; the wide v0.1.5 beam, swollen by the defect-2 and defect-3 duplicates,
never drains. So fixing 1 is necessary on its own, and fixing 2 and 3 removes the bloat
that surfaced it.

The scope of this spec is these three fixes plus the tests that would have caught
each. Two LDH-run questions are already resolved and out of scope: the drop in
fitted+inherited equation counts (the intended effect of the PR #59 dedup and PR #60
EqualAI collapse) and the 14 six-param initial mechanisms (legitimate iso/segment
topologies, fully identifiable).

## Bug 1 — beam search never terminates

### Root cause

Kinetic grouping and Wegscheider both force parameters to share a value — grouping by
modeling choice, Wegscheider by thermodynamics. They overlap whenever a group's members
are Wegscheider-tied. `_expand_split_kinetic_group` (`src/mechanism_enumeration.jl:1249`,
and the `AllostericMechanism` method at `:1260`) knows only about grouping: it splits
every group of two or more steps, including splits that separate constants a Wegscheider
cycle forces equal. Such a split is a no-op — same reduced equation, same `n_params` —
and PR #59's `_canonical_mechanism` merges it straight back.

The beam expands canonical mechanisms (`_process_batch` PASS 1 stores
`mech = _canonical_mechanism(m0)`), so it expands a canonical parent `P`, the split adds
a group, and canonicalization returns `P`: **`canonical(split(P)) == P`**. The canonical
expansion graph has a self-loop at every mechanism with a splittable Wegscheider-tied
group. `_beam_search` (`src/identify_rate_equation.jl:679`) loops
`while !isempty(frontier)` and re-ingests the self-loop child, which carries the parent's
`eq_hash` and so reads as inherited ("0 new fits") while being re-queued and re-expanded
forever.

Measured on the LDH mechanisms: expanding the canonical pool yields **2736 splits that
canonicalize back to their parent** — 100% of the self-loops — against 1632 genuine
splits. This is new in v0.1.5: `_canonical_mechanism` is PR #59. Before it, a no-op split
stayed a distinct mechanism and the split chain terminated by exhausting partitions;
v0.1.4's beam also drained through the `max_param_count` cap before any fixpoint. The
wider v0.1.5 beam does not drain, so the self-loop runs unbounded.

### Fix

The split move canonicalizes its output and drops the no-ops: for each candidate
`child`, compute `Q = _canonical_mechanism(child)`, and keep `Q` only when
`Q != _canonical_mechanism(m)`. This is memory-free, at the source, and the direct
statement of the rule above — do not un-share what Wegscheider re-shares.

This is also a *simplification*, not an added filter, because the split move is the only
source of non-canonical mechanisms. Measured across the expansion moves: init mechanisms
are canonical (0/69), and every non-split move produces canonical children
(`re_to_ss` 0/3192, `to_allosteric` 0/2438, `change_allo_state` 0/2238, …); only
`split_kinetic_group` does not (2736/4368, all of them the no-op self-loops — genuine
splits are already canonical). So once the split move returns canonical output, **all of
`expand_mechanisms` is canonical**, and the `_canonical_mechanism(m0)` call in
`_process_batch` PASS 1 (`src/identify_rate_equation.jl`) is redundant — use `m0`
directly. With PASS 1 no longer canonicalizing, the split move becomes the sole caller of
`_canonical_mechanism` and `_merge_tied_kinetic_groups`, so both move from
`src/rate_eq_derivation.jl` to `src/mechanism_enumeration.jl`, beside
`_expand_split_kinetic_group` (they keep referencing `_build_wegscheider_rename_map`,
which stays in place). The `eq_hash` dedup in PASS 1 stays: it handles the separate case
of distinct graphs that render to the same equation, and same-canonical mechanisms are now
structurally `==`, so `unique!` collapses them before fitting. Net: less code, and
canonicalization runs only on split children instead of every child.

It guarantees termination. Every non-split move and every genuine split that adds a
parameter strictly raises `n_params`; the only parameter-preserving edges are genuine
"delta-0" splits — `canonical(child)` differs from the parent yet carries the same
`n_params` — and those form a **directed acyclic graph** (verified: 855 nodes, 48 edges,
no cycle over `n_params` 5–9). With no cycle on any edge, the finite expansion graph is a
DAG and the frontier drains.

It costs no diversity. A dropped split is byte-identical to its parent. Distinct
mechanisms that merely render to the same equation are untouched — they keep their own
frontier entries and only share a *fit*.

Rejected alternative: keying the frontier on `eq_hash` ("ingest only new fits"). It masks
the loop instead of removing the pointless splits, and it drops distinct mechanisms that
share an equation — a real loss of search diversity.

### Test gap and regression

The only end-to-end `identify_rate_equation` test
(`test/test_identify_rate_equation.jl:245`) runs the narrowest possible beam —
`min_beam_width=1`, `loss_rel_threshold=1.0`, `loss_abs_threshold=0.0`, "only the
strictly-best mechanism passes per level" — on a uni-uni reaction. That beam drains
through the cap and terminates trivially. No test runs a wide beam.

Add two regressions:

- A structural invariant: every mechanism `expand_mechanisms` returns is canonical
  (`_canonical_mechanism(child) == child`). Cheap, it pins the split-move fix, and it
  guards the PASS-1 canonicalization removal — any future move that emitted a
  non-canonical mechanism would fail here. It also implies no self-loop survives (a
  self-loop child is non-canonical).
- A wide beam (`min_beam_width` well above 1, loose thresholds) on a reaction rich enough
  to generate same-count expansion, asserting the search terminates within a bounded
  iteration count. Under the current code this hangs; under the fix it terminates.

## Bug 2 — a dependent parameter leaks into `fitted_params`

### Root cause

The derivation partitions parameters into a dependent set (`dep`, each expressed
through a Haldane or Wegscheider constraint) and an independent set (`indep`, returned
as `fitted_params`). A dependent parameter must never appear in `indep`.

The generic path enforces this by construction. `_dependent_param_exprs_kernel`
(`src/thermodynamic_constr_for_rate_eq_derivation.jl:432`) returns
`indep = all_params ∖ keys(dep)` — a true complement, so overlap is impossible.

The allosteric override does not. `_dependent_param_exprs(::AllostericEnzymeMechanism)`
(`src/rate_eq_derivation.jl:1555`) hand-assembles `indep` as a union of five
independently-computed lists (`:1603`):

```julia
merged_indep = (indep_A..., indep_I_list..., reg_params_a..., reg_params_i_indep..., :L)
```

The segments carry inconsistent filters. `indep_A` already excludes `dep_A`. But
`indep_I_list` (`:1583`) filters on `p ∉ a_set && p in S_I && p ∉ collapse_targets` and
omits `p ∉ keys(dep)`. The reverse catalytic rate is dependent in the A-state (in
`dep_A`, hence outside `indep_A` and `a_set`), yet the dead I-state references it as an
unpinned weight (hence in `indep_I` and `S_I`). The filter re-admits it, and it lands
in both `dep` and `indep`. At runtime the rate body destructures it from `params` and
then overwrites it with its Haldane value — a fit dimension with zero effect on the
loss.

### Fix

Restore the invariant `indep ∩ keys(dep) == ∅` by construction, without the
regressions a naive full-set complement would introduce.

The obvious move — `indep = full_param_set ∖ keys(dep)` over
`_enumerate_parameters_full_allosteric` — fails two ways. That enumerator over-emits
unreferenced dead-I mirrors (harmless in the `Full` name list, but new phantom fit
dimensions here), and it orders reg-I before reg-A where `merged_indep` orders reg-A
first. Both would corrupt `fitted_params` order and `eq_hash`, which the Canonical Step
Form treats as load-bearing.

Instead, keep the ordered, `S_I`-filtered candidate list `merged_indep` already builds
and apply one uniform `p ∉ keys(dep)` filter across all of its segments, replacing the
per-segment filters. This is an ordered complement over the body-parameter set. It
preserves order and the `S_I` filter, and it closes the symmetric leak (a parameter
dependent in the I-state but independent in the A-state) as well as the one observed.

### Tests

- Assert `isempty(intersect(keys(dep), indep))` for every entry in
  `MECHANISM_TEST_SPECS`. Existing tests check `allunique(fitted_params)` (the symbol
  appears once, so it passes) and count `dep` and `indep` separately; none intersect
  them.
- Promote the LDH mechanisms that trigger the leak (currently in
  `LDH_ISTATE_FAILURE_MECHS`, checked only for `isfinite`/no-undefined-symbols/perf)
  into `MECHANISM_TEST_SPECS`, so they get golden parameter lists and the invariant.
  Regenerate the golden after the fix; the phantom parameter drops from their
  `fitted_params`.

## Bug 3 — the enumerator emits indistinguishable allosteric mechanisms

### Root cause

`_expand_to_allosteric` (`src/mechanism_enumeration.jl:1554`) emits, per catalytic
multiplicity, an all-`:EqualAI` baseline (`:1560`) plus one variant per group with that
group flipped to `:OnlyA` (`:1565`). Two of these are empirically indistinguishable
from a simpler mechanism — the conformational constant `L` does not affect the rate:

- **All `:EqualAI`.** The two conformations are identical, so `(1 + L)` factors out of
  numerator and denominator and cancels. The rate equals the non-allosteric parent's.
- **`:OnlyA` on the catalytic step, all bindings `:EqualAI`.** The inactive state binds
  identically but cannot catalyze, so `L` enters only as a global `1/(1 + L)` scale that
  folds into `kcat`: `v = kcat/(1 + L) · shape`.

Both still list `L` in `fitted_params`, so they pass every count-based test.

### The distinguishability principle

An allosteric split is worth enumerating only when the data can resolve the
conformational shift. What resolves it depends on where the states differ:

- **Binding differs (`:OnlyA`/`:NonequalAI` on a binding step) — K-type.** The substrate
  itself binds the states differently, so varying substrate concentration reveals `L`.
  Identifiable on its own.
- **Catalysis differs (`:OnlyA`/`:NonequalAI` on the catalytic step) — V-type.** The
  substrate binds both states alike and cannot reveal the shift; only a regulator's
  concentration can. Identifiable only alongside an `:OnlyA` or `:OnlyI` regulator.
- **Nothing differs (all `:EqualAI`).** No observable resolves the states, with or
  without a regulator. Never worth enumerating.

The all-`:EqualAI` baseline is not a required bridge to the good mechanisms: the
distinguishable per-group-`:OnlyA` variants are emitted directly as siblings, and
`_expand_change_allo_state` relaxes an `:OnlyA` group to `:NonequalAI` as readily as it
would relax the baseline. So dropping the degenerate emissions costs no reachability.

### Fix

Emit from `_expand_to_allosteric` only distinguishable mechanisms:

- Drop the all-`:EqualAI` baseline.
- Keep the binding-step `:OnlyA` variants (K-type).
- Emit the catalytic-step `:OnlyA` variant only paired with an `:OnlyA` or `:OnlyI`
  regulator, as a combined move (V-type). For a reaction that declares no regulators,
  the V-type does not exist and the variant is simply not emitted.

`_expand_change_allo_state` already collapses a related case — a lone binding group
relaxed to `:NonequalAI` while its partners stay `:EqualAI`, forbidden by a
thermodynamic cycle (`K_I = K_A`, Δ = 0), covered by `test/test_allosteric_collapse.jl`.
Confirm the same principle governs both moves.

### Tests

- Update the Δ-count assertions that baked in the phantom `L`
  (`test/test_mechanism_enumeration.jl:2837` and `:2892`, currently `[1,1,2,2,2,2]`),
  and the assertions that require the all-`:EqualAI` mechanism to be produced (`:2929`,
  `:2968`). The comments at `:3565` and `:3622` already flag the degeneracy as a
  deferred follow-up; this is that follow-up.
- Assert the distinguishability rule structurally on the expansion moves: no all-`:EqualAI`
  mechanism is emitted, a catalytic-step `:OnlyA` mechanism appears only paired with an
  `:OnlyA`/`:OnlyI` regulator, and binding-step `:OnlyA` mechanisms are emitted bare.

A blanket rank-based identifiability invariant (`rank(∂v/∂θ) == length(fitted_params)`
over emitted mechanisms) is *out of scope*. It belongs to a later identifiability pass,
along with a known residual it would flag: a rare iso/segment mechanism that confounds an
isomerization constant with a binding constant, outside the catalysis/`:EqualAI` rules
here. This spec relies on the structural distinguishability rule, not a rank check.

## Sequencing

1. **Bug 1 first.** It is independent and results-preserving, and it is the defect that
   makes the tool unusable. Land it with its wide-beam termination regression.
2. **Bug 2 next.** Small, high-certainty, and it corrects the parameter counts that bias
   model selection against allosteric mechanisms.
3. **Bug 3 last.** The largest change (the V-type combined move and the distinguishability
   rule), and the one that most reduces the beam bloat behind bug 1.

## Open items

- The V-type reduction reinterprets `kcat` as an effective turnover `kcat/(1 + L)` when
  no regulator is present; under the fix such mechanisms are not emitted, so no
  reinterpretation reaches the user. Confirm this holds for every regulator-free
  reaction.
- Deferred to a later identifiability pass (out of scope here): a blanket rank-based
  identifiability invariant over emitted mechanisms, and the iso/segment
  isomerization↔binding confounding residual it would surface.
