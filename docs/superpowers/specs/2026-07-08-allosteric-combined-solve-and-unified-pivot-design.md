# Allosteric combined constraint solve + unified pivot priority — design

Status: DESIGN (spike-validated core; approved, moving to implementation plan)
Date: 2026-07-08
Branch: `allosteric-combined-constraint-solve`
Supersedes: `docs/superpowers/specs/2026-07-07-allosteric-combined-constraint-solve-design.md`
(the 2026-07-07 attempt's box/hybrid routes are abandoned — see §2 for why).

Delivery: **one change** — the shared pivot-sentinel fix (§3.2) and the allosteric
combined-solve rewrite (§3.1) ship together, with a single reviewed rebaseline pass
across non-allosteric, enumeration, and allosteric baselines. The two fixes share a
root cause (§2) and both rebaseline the allosteric goldens; bundling avoids
rebaselining them twice. Sole ownership (the parallel Issue 2-4 investigation that
produced `2026-07-08-wegscheider-pivot-nonallosteric.md` is complete and stands
down). After this lands, an LDH HPC rerun checks whether Issues 2-4 remain.

## 1. Problem

The allosteric dependent-parameter derivation
(`_dependent_param_exprs(::AllostericEnzymeMechanism)`) runs two independent
per-state solves (`_state_dependent_exprs` for `:A` and `:I`), then reconciles
cross-state affinity "splits" with a separate pass (`_split_resolution` →
`_collapse_mirror_exprs`) plus a hand-guarded merge and an `S_I` reference gate.
Every Issue-1 defect lives in those seams: the LDH v0.1.6 run produced 10,180
`UndefVarError` rows, and the shipped surgical fix (`c59f1d7`) leaves reproducers
2 and 3 as `@test_broken` — callable but detailed-balance-wrong (`|v| ≈ 0.04,
0.16` at `Q = Keq`).

The goal: replace the fragmented per-state-solve-plus-reconcile with **one
combined constraint solve**, and remove the class of bug at its source.

## 2. Root cause (found by the 2026-07-08 spike)

Two findings, both measured.

**2a. `_split_resolution` classifies splits from the A-cycle space alone.** It
imposes the affinity-split constraint `δα = 0` from every A-state cycle,
implicitly assuming each A-cycle has a live I-mirror. That is wrong in both
directions on partially-dead-I mechanisms:

- *Over-collapse:* a ligand that binds the inactive state where I cannot turn it
  over (dead-end `E_I·X`) sits in no cycle, so `K_I` is thermodynamically free —
  yet `_split_resolution` collapses it to `K_A`. Confirmed on HEAD: a uni-uni
  mechanism with `:OnlyA` catalysis and a `:NonequalAI` substrate binding
  collapses `K_I_S_E = K_A_S_E`, silently dropping an identifiable parameter
  (detailed balance can't see it; it's an identifiability loss).
- *Under-collapse:* a split pinned by a constraint that only appears once the
  I-Haldane rows are included stays free — reproducers 2/3, the `@test_broken`.

The correct constraint is the **joint** A+I cycle space: a split is pinned iff a
cycle through it is live in *both* states — Haldane *and* Wegscheider. Stacking
`[A-rows; I-rows]` and solving once encodes exactly that.

**2b. The `-1` pivot sentinel silently drops real constraints.**
`_solve_dependent_set` initializes `best_col, best_pri = 0, -1` and selects with a
strict `>` (`thermodynamic_constr_for_rate_eq_derivation.jl:396-403`).
`_step_priority` returns `-1` for free-enzyme RE binding (`:85-91`). So a binding-K
column can *never* be chosen as a pivot. Any constraint row whose surviving
columns are all binding-K's — `rhs = 0` — therefore finds no eligible pivot and is
dropped as "redundant," even though it is a real constraint. This strikes twice:

1. A cross-state affinity split reduces (`I-row − A-row`) to `[K_A:-1, K_I:+1] = 0`
   — both binding K's at `-1` → dropped → the split is never pinned.
2. An A-state binding-K Wegscheider tie (random-order binding, e.g.
   `K_A_NAD_E = K_A_NAD_EPyruvate`) is also all-binding-K at `-1` → dropped →
   unenforced. It even defeats `_state_wegscheider_rename_map`, whose internal
   kernel call drops the same row.

**Why the prior rewrite stalled.** It read the symptom ("`K_I` stays free") as "the
constraint is missing because the I-cycle is dead," and built box rows (to *add* a
constraint) and the mirror (to *re-apply* one). But the constraint was in the
matrix — the solver was throwing it away. Box over-collapsed (it re-derived the
A-only classification) and the hybrid mis-ordered live-cycle splits. Both fought a
symptom of the pivot sentinel.

**2c. The sentinel is not allosteric-only — it over-counts plain non-allosteric
mechanisms and drives the Issue 2-4 enumeration defect.** A separate investigation
(`docs/superpowers/specs/2026-07-08-wegscheider-pivot-nonallosteric.md`) hit the
same `best_pri = -1` in `_dependent_param_exprs_kernel` on a non-allosteric LDH
`Mechanism`: `fitted_params` reports 7, should be 6 — it drops the Wegscheider tie
`K_Lactate_E = K_Lactate_ENADH` forced by a closed all-RE cycle, because the tie's
only pivot column is a free-enzyme binding K at `-1`. `_build_wegscheider_rename_map`
does **not** recover it for that grouping (so this is a genuine bug, not a cosmetic
fold). The over-count is the root of the LDH beam's non-monotone `split −1` edges
(Issue 2-4 "futile enumeration"): a split re-partitions a group so its canonical
child recovers a constraint the parent's grouping lost, giving the child one fewer
param. `best_pri = typemin(Int)` fixes it: the split delta histogram goes
`{−1:64, 0:54, +1:4680, +2:496}` → `{+1:3184, +2:496}` (strictly monotone) and the
canonical pool shrinks 1513→1389. This means the pivot fix is one shared root cause
across Issue-1 (allosteric) and Issues 2-4 (non-allosteric enumeration), and it
will change many non-allosteric `fitted_params` counts and enumeration baselines —
not a no-op.

## 3. Design

### 3.1 One combined constraint solve

Replace the per-state-plus-reconcile pipeline with
`_combined_state_dependent_exprs(am)`:

1. Assemble each state's system with the already-shipped Phase-1 factoring:
   `_assemble_constraints(_state_mechanism(am, state),
   _state_wegscheider_rename_map(am, state); step_params, all_params)` for
   `state ∈ (:A, :I)`.
2. Union the columns. A shared `:EqualAI` group carries the same bare `Symbol` in
   both states and coincides; a `:NonequalAI` group contributes distinct
   `K_A_…` / `K_I_…` columns.
3. Stack the two constraint blocks `[A_A; A_I]` over the combined column space,
   `vcat` the right-hand sides, and call `_solve_dependent_set` **once**.

Cross-state ties emerge as `I-row − A-row` (the `log Keq` cancels), so honorable
splits fall in the nullspace and stay free while forbidden splits are pinned —
for live *and* dead I-cycles, with no separate classification pass.

Both consumers share this one solve: `_dependent_param_exprs` (the fitted-param
list) and `_build_dep_assignments` (the body/`rate_equation_string` constraint
lines). Regulator-site affinities complete no catalytic cycle and stay
independent (an `:EqualAI` regulator's I-name mirrors its shared A-name); `L` is
always independent.

### 3.2 Unified, sentinel-free pivot priority

One shared priority scheme for allosteric and non-allosteric mechanisms, with **no
forbidden / never-pivot value**:

- **Keep the preference ordering** (internal-iso > metabolite-step >
  free-enzyme-binding; reverse > forward within a step). It is load-bearing: it
  decides which parameters are dependent (reverse rates, Haldane-derived) vs
  fitted, and it is also the `argmin` naming representative (`_group_rep`). This
  ordering must be preserved exactly — shift scores by a constant, never re-rank.
- **Remove the sentinel.** Every parameter gets a real, finite priority in that
  order. `_solve_dependent_set` initializes `best_pri` strictly below every real
  priority (`typemin`), so any column can be chosen as a last-resort pivot; a row
  is dropped **only** when it is genuinely `0 = 0` (all coefficients zero after
  elimination). Concretely, either (a) `best_pri = typemin(Int)`, or (b) shift
  `_step_priority` scores to be non-negative and keep `best_pri = -1` as a pure
  "none-found" marker below all parameters. Implementation will do both so no
  parameter carries `-1` *and* the marker can never collide with a real score.
- **Audit for other magic values.** `_step_priority` is the only known source of
  `-1`; the implementation confirms there is no other never-pick path (per Denis's
  rule: no ability anywhere to mark a parameter un-pivotable).

**Collapse direction (allosteric-only input, not a second function).** In the
combined solve, among two columns of equal type-priority the I-side must be
preferred as the dependent, so a pinned split collapses onto the free A-side
(`K_I = K_A`) and a steady-state binding's reverse `koff_I` outranks its forward
`kon_I` (affinity collapses, speed stays free). This is a deterministic
state-based offset applied to the shared priorities when the two state systems are
stacked — the I-state ranks above the A-state. Non-allosteric mechanisms have no
I-state and never supply this component, so the same scheme reduces to the
type-only ordering for them.

The shared solver drops only true `0 = 0` rows; the shared type ordering is
identical for both paths; the only allosteric-specific piece is the I-above-A
state input.

### 3.3 Real allosteric parameter names

Allosteric parameters render with real, structural names (`K_A_ATP_E`,
`kon_I_S_E`, `K_A_NAD_EPyruvate`) — Denis derives these mechanisms, so the
positional textbook encoding is unnecessary. Non-allosteric mechanisms keep
positional names (`K1`, `K2`) verbatim from Segel and the source articles. Naming
is decided upstream of the solver through the `name(p, m)` chokepoint, so it is
orthogonal to §3.2 (the solver permutes `Symbol`s and is name-agnostic). The
combined solve already emits real allosteric names; the work is to regenerate the
allosteric goldens under them and confirm no positional encoding remains in the
allosteric path.

### 3.4 What gets deleted

Subsumed by the combined solve — remove after confirming no remaining usages:

- `_split_resolution` and `SplitResolution`
- `_collapse_mirror_exprs`
- `_i_state_referenced_syms` (the `S_I` gate)
- native per-state `_state_dependent_exprs`
- the hand-guarded merge in `_dependent_param_exprs` and #61's line-1610 filter
- likely `_flat_expr_syms` and `_partition_constraint_lines!` (verify)

Kept: `_thermodynamic_constraints`, `_assemble_constraints`,
`_solve_dependent_set`, `_state_mechanism` / `_state_step_params` /
`_state_all_params` / `_state_rate_polys` / `_state_wegscheider_rename_map`, and
the `_state_*` polynomial family. The `_state_wegscheider_rename_map` /
`_build_wegscheider_rename_map` folds become cosmetic-only (correctness is now
guaranteed by the sentinel-free solver); keep them for rendering, especially the
non-allosteric textbook names.

## 4. Invariants and ground-truth gates (never re-baselined)

- **Detailed balance:** `v = 0` at `Q = Keq` for arbitrary parameters
  (`test_haldane_equilibrium`, tol `1e-10`), for every allosteric spec and the
  three reproducers.
- **Dependency graph:** acyclic; every RHS symbol defined (`_dep_graph_is_sound`);
  `rate_equation` callable without `UndefVarError`.
- **Reference oracles:** `test_reference_qssa`, `test_ode_steadystate`,
  `analytical_rate_fn`.
- **Performance contract (hard gate):** `rate_equation` allocation-free and
  sub-120 ns for every mechanism in `MECHANISM_TEST_SPECS`
  (`test_rate_equation_performance`). The combined solve changes the compiled body;
  this must stay green.

## 5. Test plan

1. **New TDD test — dead-I identifiability.** A uni-uni mechanism with `:OnlyA`
   catalysis and a `:NonequalAI` substrate binding: assert both `K_A_S_E` and
   `K_I_S_E` are identifiable (both present/derivable, `K_I_S_E` moves the rate),
   and detailed balance holds. RED on HEAD (over-collapses), GREEN after the fix.
2. **Flip `@test_broken` → `@test`** for reproducers 2/3 detailed balance
   (`test/test_rate_eq_derivation.jl:1301-1305`).
3. **`test_allosteric_collapse.jl`:** confirm the mirror-string and free/collapsed
   assertions still hold (the spike shows they do for the uni cases); update only
   labeling that legitimately changed.
4. **Delete `test_split_resolution.jl`** (tests deleted internals).
5. **Golden re-baseline with review.** Regenerate
   `test/reference/allosteric_golden_reference.txt` under real names; review the
   diff line-by-line (it encodes the identifiable dimension and the naming) — do
   not regenerate blindly.
6. **Non-allosteric regression + rebaseline.** The sentinel fix (§3.2) is a real
   non-allo correctness fix (§2c), so expect non-allo `fitted_params` counts,
   derivation goldens, and enumeration baselines to change. Add the non-allo
   Wegscheider regression (the LDH `K_Lactate_E = K_Lactate_ENADH` tie: assert
   `fitted_params` = 6, not 7) and the enumeration-monotonicity regression (no
   `split` move reduces the fitted-param count — the histogram has no `≤ 0` delta;
   repro in `2026-07-08-wegscheider-pivot-nonallosteric.md`). Rebaseline every
   affected non-allo derivation/enumeration test to the corrected values, each
   change understood, not accepted blindly.
7. **Performance gate** (§4).

## 6. Evidence (2026-07-08 spike)

Reconstructed the plain combined solve (from the abandoned `hybrid-incomplete.diff`,
minus the mirror), added the §3.2 priority remap, measured. Patch saved at
`docs/superpowers/plans/2026-07-08-combined-solve-spike.patch` (net −34 lines
before any §3.4 deletion).

| Case | Result | `\|v\|` at `Q=Keq` |
|---|---|---|
| uni live forbidden (`NonequalAI,EqualAI,EqualAI`) | collapses `K_I_S_E=K_A_S_E` | 1.7e-17 |
| uni two `NonequalAI` bindings | exactly 1 honorable DOF | 5.6e-18 |
| uni `NonequalAI` catalysis | `K_I_S` free + identifiable | 1.1e-17 |
| **dead-I** (`NonequalAI,OnlyA,EqualAI`) | **`K_I_S` free + identifiable** | 1.1e-17 |
| SS binding + `EqualAI` catalysis | `koff_I` collapses, `kon_I` free | 5.9e-17 |
| SS `NonequalAI` catalysis | both free | 5.0e-17 |
| reproducer 1 | detbal OK | 2.5e-18 |
| reproducer 2 (`@test_broken`) | **detbal OK** (was 0.007) | 1.9e-18 |
| reproducer 3 (`@test_broken`) | **detbal OK** (was 0.14) | 4.3e-17 |

Param counts dropped to the true identifiable dimension (repro 2: 12→10, repro 3:
13→12). Not yet run: the full suite, goldens, non-allo regression, and the perf
gate — those are the implementation's job (§5).

## 7. Risks and scope

- **Shared-kernel change (§3.2) touches the load-bearing non-allosteric path and
  changes real output** (over-count fix, §2c) — not a no-op. It will break many
  hardcoded non-allo `fitted_params`/derivation/enumeration tests. Mitigation:
  preserve the exact relative order; treat the full non-allo suite as a
  bug-fix rebaseline (§5.6), each change understood. Delivered as one change under
  sole ownership (see Delivery note), so the fix lands once.
- **Golden re-baseline** is expected (labeling/pivot direction, real names).
  Mitigation: line-by-line review (§5.5).
- **Performance contract** could regress if the combined-solve body allocates or
  slows. Mitigation: the §4 perf gate; if it regresses, STOP and discuss.

## 8. Resolved decisions

- **Fractional split coefficients — ACCEPT.** A cycle traversed with multiplicity
  can yield a rational exponent (e.g. `K_I = K_A · (…)^(1//2)`, a geometric-mean
  affinity). The rational solver renders it natively; the old `_split_resolution`
  errored. Decision (Denis): allow `1//2` rational exponents — no known mechanism
  forces one, and no reason to forbid them.
- **I-above-A offset encoding — lexicographic tuple.** Encode the pivot priority as
  `(is_I_state, type_priority)` compared lexicographically, so there is no magic
  offset constant. (Denis: no strong preference; tuple is fine.)
