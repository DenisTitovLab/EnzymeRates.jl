# Strict `:EqualAI` Allosteric Derivation — Affinity/Speed Collapse (No Case-B) — Design

**Date:** 2026-07-06
**Status:** Implemented (branch `strict-equalai-allosteric`); full suite green (21864/21864).
**Supersedes:** `2026-07-05-allosteric-nonequalai-itwin-orphan-fix-design.md`. That fix repaired a crash that is a *symptom* of the mechanism this spec removes (Case-B); under the strict rule the offending configs collapse instead of deriving. Only the all-`:EqualAI` regulator-guard removal from that branch carries forward.
**Reuses (not the filter):** the split-constraint linear algebra sketched in the deprecated `2026-05-29-nonequalai-rank-validity.md`. That document's *filter* (reject degenerate configs) was dropped; its *math* (cross-state split constraints over the cycle space) is the basis of the collapse here — but the resolution is **collapse to a consistent degenerate equation**, never rejection.

> **Design update (as built — "Option 3").** The mechanic below was refined during implementation. The original framing collapsed one `δ_g` per `:NonequalAI` *group* and treated every steady-state group as always-free. That is wrong for a steady-state **binding** whose affinity is forbidden by a Wegscheider box (`m_ro`): the box does not close and equilibrium flux is nonzero. The mechanic now decomposes each reversible step into an **affinity** (`Kd` / `kon·koff⁻¹`, cycle-constrained) and — for steady-state steps — a **speed** (`kon·koff`, always free); it collapses only *forbidden affinities*, keyed on the base **free/derived (`indep_A`) partition** rather than on binding-vs-catalytic step type. §2 describes the mechanic as built. Two consequences: the deferred **direction-symmetry (D2)** work is **dropped** — its "shared speed, different affinity" resolution is not physically motivated and strict + affinity-collapse supersedes it; and the collapse is **rational (no `√`)**, because it collapses ratios rather than sharing speeds.

## 1. Context and root cause

The allosteric derivation runs the shared King–Altman/Haldane–Wegscheider engine once per conformational state (A/I). An `:EqualAI` group is meant to render one shared symbol in both states; a `:NonequalAI` group renders distinct `K_A_…`/`K_I_…`.

**The defect ("Case-B").** When an `:EqualAI`-tagged *dependent* parameter's derived expression references a `:NonequalAI` symbol, its value genuinely differs between states, so the current code silently **promotes** it to a distinct I-name (`_case_b_rename_map`, `src/rate_eq_derivation.jl:1277`, applied in `_state_dependent_exprs(am,:I)` ~1301 and `_state_rate_polys(am,:I)` ~1169). Two shapes hit it: a catalytic Haldane's SS **reverse rate** (PK: `k5r → k5r_T`, referencing the `:NonequalAI` PEP binding) and a binding-Wegscheider box **pivot** (`K_B_EA → K_I_B_EA`).

This is not transparent. A user tags a step `:EqualAI` expecting it shared, and a parameter of that step silently un-shares. Which parameter absorbs the difference is decided by the (structural, tag-blind) elimination pivot, not the user. The physics is unavoidable — in a thermodynamic cycle you cannot make one affinity differ between conformations without a partner differing — so Case-B was *making the tags mean something the user did not write*: "catalysis is identical" is false the moment the substrate binding it is tied to differs.

**Decision (Denis).** `:EqualAI` means genuinely shared, without exception. A `:NonequalAI` split is honored only to the extent that it forces no `:EqualAI` parameter to differ. The forbidden component **collapses** to a valid, thermodynamically consistent, degenerate equation (`K_A = K_I` along that direction) — never an error. Real allosteric models therefore tag the coupled steps explicitly: the catalytic step `:NonequalAI` (so its reverse differs natively), the coupled binding edges `:NonequalAI`, or `:OnlyA` (substrate binds only the active state — the classic K-system, handled by graph pruning, not Case-B). Detecting and pre-empting collapsing configs in the enumerator is **out of scope** (a later PR).

## 2. The rule and the collapse mechanic

**Rule.** Across the two states, every `:EqualAI` parameter — independent or derived — is shared. Each reversible step carries two independent quantities: an **affinity** (`Kd`, or `kon/koff` for a steady-state binding) that thermodynamic cycles constrain, and — for a steady-state step — a **speed** (`kon·koff`) that no cycle constrains and is therefore always free. The `:NonequalAI` *affinity* splits are constrained to the subspace in which all `:EqualAI` parameters stay shared; the orthogonal complement collapses. A steady-state binding's **speed split is always free** — only its affinity can be forbidden.

**Collapsibility is keyed on the free/derived partition, not the step type.** A group contributes a collapsible affinity iff it has an *independent* affinity in the base thermodynamic derivation (`indep_A`): a rapid-equilibrium `Kd`, or a steady-state binding whose forward AND reverse rate constants are both independent (a non-pivot binding). A steady-state group with only one independent rate constant — its reverse is a derived cycle pivot (every catalytic step, and any binding chosen as a Wegscheider pivot) — has its affinity already absorbed by that derived reverse; its split is always free. Keying on `indep_A` rather than "binding vs catalytic" is load-bearing: a `:NonequalAI` binding that happens to be a box pivot is structurally identical to a catalytic step (one free rate constant, one derived), and the free/derived test classifies it correctly with no special case and no edge.

**Mechanic (`_split_resolution` + `_collapse_mirror_exprs`).**

1. Compute the base per-state derivation, yielding `indep_A` (the free rate constants) and each `:EqualAI` dependent's expression. Classify each `:NonequalAI` group by `nfree(g)` = its rate constants in `indep_A`: RE `Kd` or 2-free SS binding ⟹ collapsible affinity; 1-free SS ⟹ absorbed (always free).
2. Build the affinity-constraint matrix `M` over the collapsible-affinity columns from the thermodynamic cycle incidence `C` (an RE `Kd` enters as `−C`; a steady-state affinity `kon/koff` as `+C` — the same affinity coordinate `−C·δ_Kd`, opposite raw sign). Order the absorbed columns first so they eliminate before the collapsible columns partition — a derived relation then references only free collapsible affinities.
3. The RREF partition of `M` gives the honorable affinity-split space: free columns keep genuine DOF; each pivot column's affinity is **derived** from the free ones.
4. Emit, in terms of each group's effective dissociation constant `effK` (an RE `Kd`; `koff/kon` for a 2-free SS binding), tying `effK_I_g/effK_A_g = ∏(effK_I_f/effK_A_f)^a`:
   - **RE group** — derive the I `Kd`: `K_I_g = K_A_g·∏(effK_I_f/effK_A_f)^a`.
   - **SS group** — derive the I *reverse rate*, keep the forward (speed) free: `koff_I_g = koff_A_g·(kon_I_g/kon_A_g)·∏(effK_I_f/effK_A_f)^a`.
   A trivial nullspace is full collapse (`K_I=K_A`, or `koff_I=koff_A·kon_I/kon_A`); all rational (no `√`), built by the existing `build_power_expr`. Derived I-symbols drop from `fitted_params`; `:EqualAI` dependents stay a single shared symbol (no promotion).

This single **split resolution** feeds both the polynomials and the dependent-assignment builder, replacing Case-B's two separate rename sites and their coupling.

**Behavior by case.**
- Single `:NonequalAI` binding + `:EqualAI` catalysis → full collapse: RE `K_I = K_A`; SS affinity shared (`koff_I = koff_A·kon_I/kon_A`). Degenerate — use `:OnlyA` or tag catalysis for a real model.
- Two coupled bindings + `:EqualAI` catalysis → **one honorable DOF survives** (`K_I_P = K_A_P·K_I_S/K_A_S`), catalysis genuinely shared. A valid K-system, preserved. The coupled bindings may be `:OnlyA` **or** several `:NonequalAI` steps moving together so `K_A_P/K_A_S = K_I_P/K_I_S`.
- Steady-state binding whose affinity is forbidden (`m_ro`: a `:NonequalAI` SS binding in a Wegscheider box) → its affinity collapses but its **speed stays free**: the two conformations bind with the same `Kd` but different kinetics — an identifiable steady-state allosteric DOF, preserved (where full collapse would flatten it away).
- Catalysis (or any step) tagged `:NonequalAI` → its reverse differs natively; no `:EqualAI` dependent forced, no collapse.

**Verified consistent.** Equilibrium net flux is machine-zero on every collapsed mechanism (the hard invariant — `m_ro` moved from `−0.016` to `4e-18`), and the fully-collapsed single-split case reproduces the all-`:EqualAI` equation. The 0-flux-at-equilibrium invariant is a hard test on every collapsed mechanism.

## 3. What gets removed

**Delete (Case-B cluster, all `src/rate_eq_derivation.jl`, all Case-B-exclusive):** `_case_b_rename_map` (1277), `_i_nonequalai_syms` (1256), `_state_i_case_b_renames` (1290), `_dep_inactive_name` (1386). Remove both application sites (`_state_dependent_exprs(am,:I)` ~1301-1308 and `_state_rate_polys(am,:I)` ~1169-1171) together — removing one desyncs the polynomials from the dependent names.

**Orphan-fix closure — N/A (built clean off main).** The superseded orphan-fix branch added a transitive-closure loop to `_i_state_referenced_syms` to keep an inner-edge I-twin alive. This work was built fresh off `main`, where that closure never landed — `_i_state_referenced_syms` is already the poly-only form. Nothing to delete; the orphan it guarded no longer arises, because an inner-edge I-twin is now either a genuine free split (coupled tagging / free SS speed) or a collapsed derived symbol.

**Simplifies:** `_state_rate_polys` and `_state_dependent_exprs` lose their `:I` special-case branches and become single-path. The base `S_I` poly-referenced-symbol filter **stays** (still needed for genuine `:NonequalAI` distinct I-names).

**Keep (shared, only their Case-B call removed):** `_expr_references_any`, `_flip_to_inactive`/`_force_inactive`/`_param_for_symbol`.

**Pattern to follow for the mirror emission:** the `:EqualAI` regulator mirror already emits `K_I_reg = K_A_reg` in `_dependent_param_exprs` (~1467-1476) and `_build_dep_assignments` (~1559-1566), dropping `K_I_reg` from `fitted_params`. The catalytic collapse mirror follows the same convention, extended to a derived *expression* (not just `= K_A_sym`) for the coupled case.

**Guard the merge.** The dep merge `for (k,v) in dep_I; k in S_I && (dep[k]=v)` (`_dependent_param_exprs` ~1447-1453) silently overwrites the A-state dependent with the I-state RHS when a shared `:EqualAI` dependent has the same bare name in both states. The split resolution must run so the shared dependent carries one consistent value and the collapse mirrors are emitted before any dependent that reads them.

## 4. Fixture retagging

Deriving all nine allosteric `MECHANISM_TEST_SPECS` shows **only PK relies on Case-B** (its `k5r_T`). `m_all` is already strict-compliant (catalysis already `:NonequalAI`; its lone `:EqualAI`-forced difference is a regulator that correctly mirrors). PFK-1, HK, `m_OnlyA_prod` use `:OnlyA`; the four all-`:NonequalAI` fixtures are native.

**PK → Option B (PEP `:OnlyA`).** Tag PK's PEP binding `:OnlyA` (exclusive-binding K-system, graph pruning, not Case-B): the T-state cannot bind PEP, the T-catalytic cycle is dead (`N_T = 0`), catalysis stays `:EqualAI`, and `kcat = k5f` is preserved. `fitted_params` 9 → 8. **Retag the SOURCE PEP step in the `@allosteric_mechanism_src` block, not a stored index** (canonicalization reorders groups). Rewrite `pk_rate_analytical` (`test/…:2222` — zero `N_T`, strip all T-state PEP / `K1_T` terms from `Q_cat_T` and the `L·num_T` flux), keep `analytical_kcat_fn = p -> p.k5f`, set `expected_n_independent_params` 9 → 8, and regenerate PK's golden lines (`test/reference/allosteric_golden_reference.txt` ~25-28).

## 5. Test strategy

- **Collapse consistency (new, the anchor).** For each collapse shape — single `:NonequalAI` binding + `:EqualAI` catalysis (full collapse), and two coupled bindings + `:EqualAI` catalysis (one honorable DOF) — assert: derivation succeeds (no error), `rate_equation` evaluates finite, equilibrium flux ≈ 0, the full-collapse case emits `K_I_g = K_A_g` and its `fitted_params` matches the all-`:EqualAI` baseline, and the honorable-DOF case emits the coupled mirror `K_I_P = K_A_P·K_I_S/K_A_S`, keeps one split free, and its rate responds to that surviving split (identifiable).
- **Branch suite flips (`test/test_allosteric_wegscheider.jl`).** The two testsets that assumed Case-B — "I-twin retained for inner-independent edges" and "retained split is identifiable" — become: single-inner-edge → asserts collapse (`K_I_A_EB = K_A_A_EB`, no new free split); two-coupled-edge (steps 3 **and** 4 `:NonequalAI`) → asserts the identifiable coupled split. The robustness testsets (derives / evaluates / equilibrium-flux ≈ 0 across placements, `:OnlyA`, regulators, multi-tag) stay green *because* the collapse is a graceful mirror, not a rejection.
- **PK golden + oracles.** `test_allosteric_golden.jl` (byte-exact) and the `analytical_rate`/`analytical_kcat` oracles must pass with the retagged PK; regenerate and review the golden.
- **Delete the Case-B unit tests** that assert the old promotion (`test_rate_eq_derivation.jl` ~1660-1663 for `_state_i_case_b_renames`, and the `_dep_inactive_name` testset ~1727-1762).
- **Hard gates:** `rate_equation` 0-allocation / < 100 ns (`test_rate_equation_performance`) — the collapse is compile-time only; mirrors are as cheap as the existing reg mirrors. `fitted_params`/`metabolites` stay `@generated`. The parameter-naming chokepoint AST walker stays green (all symbols via `name(p,m)`).
- **Golden churn is intended.** Every mechanism where Case-B fired changes output; this is a deliberate semantic change to what `:EqualAI` means, not a regression.

## 6. Non-goals

- **No enumerator change.** Pre-empting collapsing configs in `_expand_change_allo_state` is a separate PR.
- **No rejection.** Collapsing configs derive a valid degenerate equation for teaching; they never error.
- **No D3 validity filter.** The split-constraint solve is used to *resolve* the derivation, not to reject mechanisms.
- **Direction-symmetry (D2) is dropped, not deferred.** D2's model-changing part resolves a forced state-difference by sharing the *speed* and letting the *affinity* differ (`k_for_I·k_rev_I = k_for_A·k_rev_A`, which needs `√`). That resolution is not physically motivated (there is no reason to conserve `kon·koff` while `Kd` differs), and strict `:EqualAI` + affinity-collapse is the opposite, physically-grounded policy (collapse the forbidden affinity, keep the free speed) — so it supersedes D2 rather than composing with it. D2's cosmetic part (single-state direction-invariant base parametrization) does not change any model and is not pursued.

## 7. Sequencing and open verification points

- **Build fresh; do not stack on the orphan closure** — it is being deleted. This work supersedes the orphan-fix branch; carry only the all-`:EqualAI` regulator-guard removal.
- **Free-split basis choice.** The RREF partition (which coupled affinity is free vs derived) is a structural pivot choice; it is deterministic and order-stable (absorbed columns first, then collapsible by group index), so `fitted_params` and the golden are reproducible.
- **Coverage of transitive references.** Verify the split resolution mirrors *every* `:NonequalAI` I-symbol a shared dependent's RHS can reach (including a Wegscheider dependent that chains through an inner binding constant), so no undefined symbol survives into the generated body — the failure mode the deleted orphan closure was guarding.
- **`parameters(m, Full)` unaffected.** Full mode over-emits raw names without Haldane reduction; the collapse is a Reduced-mode concept. Confirm the Full accessor still lists the expected names.
