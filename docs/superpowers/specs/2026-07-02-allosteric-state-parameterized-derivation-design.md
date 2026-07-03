# Allosteric D1 — State-Parameterized I-State Re-Derivation — Design

**Date:** 2026-07-02
**Status:** Approved design; implementation pending.
**Scope:** D1, the keystone of a five-part allosteric refactor (Denis's points 1+2+3, plus the point-4b per-monomer normalization).
**Companions:**
- `2026-05-29-direction-symmetry-constraint-resolution.md` — D2 (point 5), the follow-up that replaces single-pivot elimination with symmetric speed/ratio resolution. Builds on D1's clean kernel.
- `2026-05-29-nonequalai-rank-validity.md` — D3 (point 6), the rank/split-freedom validity filter. Builds on D1's factored cycle basis.
- D4a (point 4a, `1^n` printing cleanup) and the `rate_equation_string` polish are explicitly **out of D1 scope**.

## 1. Context and root cause

The King–Altman/Cha derivation engine is already shared: an allosteric mechanism's A-state numerator/denominator come from the plain, non-allosteric `_raw_symbolic_rate_polys(CM)`, and the emitted rate equation already has the MWC shape `E_total · CatN · (num_A + L·num_I) / (den_A + L·den_I)`.

The excess complexity comes from one design choice: **the I-state is produced by *renaming* A-state symbols** (`:NonequalAI`→`:I`, and dropping `:OnlyA` monomials) **rather than *re-deriving* on an I-tagged copy of the catalytic mechanism.** Because the dependent-parameter reduction (`_dependent_param_exprs_kernel`) runs only in A-state names, every I-state Haldane/Wegscheider relation must be *reconstructed* in I-names after the fact. That reconstruction is the bulk of the allosteric code: the Case-A/Case-B synthesized-dependent apparatus, three overlapping rename maps, five near-identical tag-walking loops, a Symbol-vs-Parameter dual accounting, and a re-inlined copy of I-state polynomial construction inside `_kcat_forward`. It is also the source of the §5a-class desync bugs (four functions that must stay in lockstep).

## 2. Approach — tag the input, do not rename the output

`_raw_symbolic_rate_polys(mech, step_params, rename_map, subs, prods)` and `_dependent_param_exprs_kernel` determine every polynomial and constraint symbol from the `step_params` they receive, through the `name(p, mech)` chokepoint. So we build **state-tagged `step_params` and a per-state catalytic graph up front**, then run the existing engine once per conformational state:

- **A-run:** `:NonequalAI` and `:OnlyA` groups tagged `:A`; `:EqualAI` groups tagged `:EqualAI`; `:OnlyA` groups present. Yields `num_A`, `den_A`, and the A-state dependent-parameter assignments — all in native A-state names.
- **I-run:** `:NonequalAI` groups tagged `:I`; `:EqualAI` groups tagged `:EqualAI`; **`:OnlyA` steps pruned from the graph.** Yields `num_I`, `den_I`, and the I-state dependent assignments natively in I-state names.

The MWC assembly then combines them into the standard form, `E_total · (num_A + L·num_I) / (den_A + L·den_I)` after the §5 per-monomer normalization (or `E_total · CatN · (…)` before it, in commit 1), with the reg-site partition factors and `L` unchanged.

Two alternatives were considered and rejected:
- **Factor-the-rename:** keep derive-once, collapse the rename bookkeeping into a single helper. Smaller payoff; retains the Symbol-vs-Parameter dual accounting and the synthesized-dependent concept.
- **Single combined A+I enzyme-form graph:** run King–Altman over one graph containing both conformations linked by `L`. Breaks the `num_A + L·num_I` factorization and risks the 0-allocation `rate_equation` contract.

## 3. Per-state handling

- **EqualAI sharing is automatic.** `_state_tag(:EqualAI)` renders the empty string, so an `:EqualAI` parameter produces the same bare Symbol in both runs (e.g. `k_EADPPEP_to_EATPPyruvate`, `K_ADP_E`). No bridge assignments; printed names are unchanged from today.
- **NonequalAI splits naturally.** A `:NonequalAI` group renders `K_A_…` in the A-run and `K_I_…` in the I-run, so the two states get genuinely distinct symbols and genuinely distinct per-state constraints.
- **OnlyA is graph construction, not monomial deletion.** An `:OnlyA` step is absent from the I-run's graph, so King–Altman re-derives the correct broken-cycle I-state law natively. This replaces the current `N_I = 0` patch (and its duplicate inside `_kcat_forward`) with a derived result. `:OnlyA` catalytic mechanisms must be added to `MECHANISM_TEST_SPECS` to pin this (Denis's point 3).
- **Shared Keq needs no special handling.** `Keq` is user-provided and never eliminated; each state's Haldane merely *references* the one `:Keq` literal to express that state's dependent reverse rate. PK already carries two Haldane constraints today (one per state), so per-state re-derivation reproduces the same count. There is no double-elimination to guard against.

## 4. The Case-B residue — one naming rule

One irreducible case survives: an **`:EqualAI` catalytic group whose Haldane-dependent reverse rate references a `:NonequalAI` symbol** (PK's central step: `k5r` references PEP's `K1`/`K1_T`). Its *value* differs between states, so it needs a distinct I-name even though its group tag is `:EqualAI`; the group tag alone cannot express this.

In the re-derivation model this shrinks from a reconstruction subsystem to a single local rule applied at derivation time:

> After the I-run, any bare/`:EqualAI`-named dependent parameter whose derived expression references a state-differing symbol takes its I-form name (e.g. `k_EATPPyruvate_to_EADPPEP` → `k_I_EATPPyruvate_to_EADPPEP`).

This is one function with no cross-function lockstep, and it reproduces today's naming exactly (today's `_synthesized_dep_i_names` output). It will be revisited and likely absorbed by the D2 symmetry rewrite and D3 rank work; D1 keeps it as a small, self-contained rule.

## 4a. A/I constraint consistency and the degeneracy boundary (D3 interface)

D1 solves the constraint kernel once per state, which raises the question of whether the two solves can be mutually contradictory or can silently collapse an intended `:NonequalAI` split (`K_A = K_I`). They cannot, for three structural reasons:

1. **Identical partition.** The A- and I-state graphs share topology (identical except for `:OnlyA`-pruned steps), and the kernel's pivot selection (`_step_priority`) depends only on graph structure — not on parameter values or state tags. So the two runs choose the **same dependent/independent partition**: the A-run's dependent at a position is `K_A_…`, the I-run's is `K_I_…` (or the shared symbol for `:EqualAI`). The split is structurally identical across states, only A/I-relabeled. `K_A` and `K_I` are always distinct symbols — D1 never collapses a `:NonequalAI` split.

2. **Correct equilibrium by construction.** Each state's Haldane is re-derived natively, so at chemical equilibrium `N_A = 0` and `N_I = 0` independently, giving total rate `E_total · (0 + L·0) / den = 0` **regardless of the tag configuration**. Per-state re-derivation therefore never produces a thermodynamically contradictory equation or a nonzero equilibrium rate. The kernel's own "contradictory mechanism" error can only fire on a graph the A-run would also reject; `:OnlyA` pruning only removes constraints.

3. **The Case-B rule keeps each state individually correct.** When a shared `:EqualAI` dependent is fixed by a cycle that also contains a `:NonequalAI` symbol, the two runs assign it different values; the §4 rule gives the I-form a distinct name so each state's equation uses the right value. Gaussian elimination expresses every dependent purely in terms of *independent* parameters, so the check is complete and non-transitive (a dependent never references another dependent).

What D1 does **not** do is detect **degeneracy**: a `:NonequalAI` split that the combined A/I constraints render non-identifiable. The canonical case is a lone `:NonequalAI` binding K in a pure-RE Wegscheider loop with no SS step to absorb the split. There, `K_A` and `K_I` remain distinct free symbols but the observable rate is invariant to their difference — the split is silently absorbed by the shared `:EqualAI` dependent, producing a **phantom** (over-parametrized) DOF. This is not a collapse and not a wrong equation (equilibrium is still zero, per point 2); it is over-parametrization. Detecting and rejecting it is exactly the job of **D3** (`2026-05-29-nonequalai-rank-validity.md`): its `_nonequalai_splits_free` rank test rejects or normalizes the degenerate config so the phantom never reaches the fitter.

**Decision: detection/rejection of degeneracy is deferred to D3.** D1 preserves today's behavior (it emits the phantom, as today's `_add_case_b_renames!` does) and is thermodynamically safe in the interim. No current allosteric mechanism exercises this case (there is no allosteric-Wegscheider fixture yet; existing `:OnlyA`/`:OnlyI` usage is on regulatory sites). To make the seam explicit rather than silent, D1 adds a **characterization test** (§8): the existing "Random-order Bi-Bi" mechanism made allosteric with one `:NonequalAI` binding group, asserting D1's current phantom output with a comment marking it as the config D3 will reject — which also hands D3 its canonical fixture.

## 5. Per-monomer normalization (point 4b) — decided

`E_total` is redefined as the **active-site (protomer) concentration**, and the leading `CatN` coefficient is **removed**.

Rationale: the current `v = E_total · CatN · (…)/(…)` is the standard MWC rate *per oligomer*, i.e. `E_total` is implicitly the oligomer concentration. Experimental enzyme concentration is reported per active site. With `E_total` = active-site concentration, `E_oligomer = E_total / CatN`, and `CatN` cancels exactly:

```
v = E_total · (num_A + L·num_I) / (den_A + L·den_I)
```

The `Q^(CatN−1)` / `Q^CatN` binding-statistics powers are unchanged — only the leading coefficient is a normalization artifact. This also makes allosteric `kcat` identical in meaning to non-allosteric (per-active-site turnover): PK's saturating `kcat` becomes `k5f`, not `4·k5f`.

This is a deliberate, documented convention change: fitted `kcat` shifts by a factor of `CatN` relative to today. It removes the "leading 2/4 coefficient" half of point 4 by correct normalization; the `1^n` half is a separate string-only concern deferred to D4a.

## 6. What deletes vs stays

**Deletes** (rename/reconstruction bookkeeping, ~250–400 LOC): `_A_rename_parameters`, `_a_to_i_rename`, `_I_rename_parameters`, `_raw_symbolic_rate_polys_allosteric`, `_dependent_param_exprs_allosteric`, `_i_state_num_den_polys`, `_synthesized_dep_i_names`, `_add_case_b_renames!`, `_i_state_referenced_syms`, most of the `_dependent_param_exprs(::AllostericEnzymeMechanism)` closure, and the re-inlined I-state polynomial construction inside `_kcat_forward` (which instead consumes the shared per-state derivation).

**Stays** (irreducible MWC): `_reg_site_expr` and `Kreg` naming; `Lallo`; OnlyA handling (reframed from monomial-deletion to graph construction); the one Case-B naming rule (§4); and the `num_A + L·num_I` MWC combine.

**Refactored, not deleted:** `_kcat_forward(::AllostericEnzymeMechanism)` and `parameters(m, Full)` for allosteric both consume the shared per-state derivation instead of re-deriving.

## 7. Two-commit TDD sequencing

1. **Structure-preserving re-derivation** (keeps `CatN`). Tag-the-input per-state derivation replaces the rename bookkeeping. Success criterion: **byte-identical** `rate_equation_string` and `parameters(Full)` / `parameters(Reduced)` for every allosteric mechanism that exists in `MECHANISM_TEST_SPECS` **today** (none of which has an `:OnlyA` *catalytic* group — existing OnlyA/OnlyI usage is on regulatory sites, whose assembly is unchanged), measured against today's output captured as a golden reference. The risky refactor is validated against unchanged references. New `:OnlyA`-catalytic fixtures (§8) have no prior reference and are validated by hand-derived analytical formulas instead.
2. **Per-monomer normalization** (drops `CatN`). Isolated, reviewable diff: remove the prefactor, update the analytical formulas (`4*`/`2*` removed; `analytical_kcat_fn` → per-active-site, e.g. PK `p -> p.k5f`), and document `E_total` = active-site concentration.

## 8. Test strategy

- **Golden reference (commit 1).** Before touching derivation, capture current `rate_equation_string`/`parameters` for the allosteric `MECHANISM_TEST_SPECS` entries that exist today. The refactor must reproduce them byte-for-byte. This is the failing-test-first anchor.
- **New `:OnlyA` catalytic fixtures.** Add mechanisms with an `:OnlyA` catalytic group (substrate-side and product-side) to `MECHANISM_TEST_SPECS`, with hand-derived analytical rate functions, so the natively re-derived broken-cycle I-state law is pinned. These are new coverage for point 3.
- **Degeneracy characterization fixture (§4a).** Add the "Random-order Bi-Bi" mechanism made allosteric with one `:NonequalAI` binding group, asserting D1's current phantom output and confirming equilibrium flux is zero. Marked as the config D3 will reject; hands D3 its canonical fixture.
- **Per-monomer assertions (commit 2).** Analytical formulas and `analytical_kcat_fn` updated to per-active-site; the derivation must match.
- **0-allocation / <100 ns `rate_equation` contract** is a hard gate on both commits (`test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`). The `@generated` body folds all state-tag bookkeeping and any residual `1`-factors at compile time, so the contract is expected to hold; it must be verified, not assumed.
- **Parameter-naming chokepoint guard** (`test/test_types.jl` AST walker) must stay green — all new symbols flow through `name(p, m)`.

## 9. Scope boundaries and open verification points

**Out of scope for D1:** the D2 direction-symmetry rewrite (point 5), the D3 NonequalAI rank-validity filter (point 6), and the D4a `1^n` string cleanup (point 4a). D1 is model-preserving except for the deliberate per-monomer normalization.

**To verify during implementation (not blocking the design):**
- **Wegscheider in allosteric.** No allosteric-Wegscheider fixture exists today; confirm per-state re-derivation reproduces current Wegscheider handling for a mechanism that has one. The degenerate pure-RE Wegscheider case is analyzed in §4a and its detection is deferred to D3; D1 only adds the characterization test.
- **`parameters(m, Full)` order.** The Full-mode enumeration must reproduce today's name list (post the already-merged `filter!` dedup) so the accessor contract is preserved.
- **`_kcat_forward` reg-corner assembly.** Confirm the shared derivation reproduces the analytic kcat across state × reg-corner combinations for the existing `analytical_kcat_fn` fixtures.
