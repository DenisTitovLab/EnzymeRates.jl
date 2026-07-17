# Solve-then-limit MWC derivation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the graph-deletion allosteric MWC derivation with solve-then-limit — derive both conformations in full, solve the thermodynamic constraints on the coupled system, normalize each per state, combine, then apply `:OnlyA` as a limit — so consistency holds by construction.

**Architecture:** All work is in the `@generated` allosteric path of `src/rate_eq_derivation.jl`. The single-conformation King–Altman (`_raw_symbolic_rate_polys`) and PR #70's constructor guard are reused. The n=1 mass-action harness `test/allosteric_ground_truth.jl` is the acceptance gate and is *built/strengthened first* so every derivation change is validated against it.

**Tech Stack:** Julia 1.12, EnzymeRates.jl. Tests via `TestEnv` / per-file drivers.

**Spec:** `docs/superpowers/specs/2026-07-16-mwc-solve-then-limit-derivation-design.md`. Read it before Task 1. This plan assumes `main` contains PR #70 (branch is rebased on it).

## Global Constraints

- 92-character lines, 4-space indentation.
- `rate_equation` MUST stay allocation-free and < 120 ns per call (`test_rate_equation_performance` in `test/test_rate_eq_derivation.jl`). This is non-negotiable; if a change forces an allocation, STOP and reconsider — do not weaken the gate.
- All `Parameter → Symbol` rendering flows through `name(p, m)`. No `Symbol("K…")`/`Symbol("k…")` literals (AST-walker test at `test/test_types.jl`).
- The **non-allosteric** derivation path must stay byte-identical — verify with the golden.
- **The oracle is the definition of correctness.** No derived equation ships unless it matches the strengthened `mwc_ground_truth_flux` to rtol 1e-4 and gives `v = 0` at the equilibrium metabolite ratio.
- Run Julia FOREGROUND and WAIT; never background-and-yield (a prior session orphaned a run for 40 min). The full `Pkg.test()` OOMs in a memory-limited sandbox — verify per-file; run the full suite on an adequate-memory machine.
- Commit after each task. Never skip a pre-commit hook.

---

## Background the implementer needs

**Current shape.** `_allosteric_num_den_exprs(M_type)` (`rate_eq_derivation.jl:1662`) calls `_state_rate_polys(am, :A)` and `_state_rate_polys(am, :I)`, each returning `(num, den, d_free)`. For `:I`, `_state_rate_polys` derives on the **graph-deleted** subgraph (`_state_mechanism` → `_state_allo_mechanism`, `:1214`, which drops `:OnlyA` groups). It then combines A/I with a three-way branch on `d_free_A` vs `d_free_I` (raw / divide / cross-weight). The constraint solve is `_combined_state_dependent_exprs` (`:1439`). `_kcat_forward` (`:941`) repeats the cross-weight for the saturating limit.

**Target shape.** The `:I` polys derive on the **full, undeleted** graph with I-scoped parameters; the constraint solve runs on that full two-state system; each state is normalized to its free-enzyme weight independently; the plain combine `(N_A·D_A^{n-1} + L·N_I·D_I^{n-1})/(D_A^n + L·D_I^n)` is taken; then the `:OnlyA` limit (`K_I → ∞` RE / `k_I → 0` SS) drops the vanishing monomials. The guard is the same limit applied to the solved constraints (reject on a contradiction). Per-state normalization is **algebraically the shipped cross-weight** — the recipe is known-correct; the work is re-expressing it decoupled and dropping the graph-deletion.

**Known-correct fixed points** (do not regress): the shipped equations for LDH i-state, multi-`:OnlyA`, metabolite-bearing-D, and `:NonequalAI` ping-pong all match the independent limit oracle to ≤ 3e-16 today. The rewrite must reproduce those, not change them.

---

### Task 1: Strengthen the ground-truth oracle (the acceptance gate)

Build the gate before touching the derivation, so every later task is checkable.

**Files:**
- Modify: `test/allosteric_ground_truth.jl`

**Interfaces:**
- Produces: an oracle that encodes `:OnlyA` as the **limit** (`K_I` large, `k_I → 0`), plus a ping-pong value gate and an `N_I ≠ 0` multi-protomer gate.

- [ ] **Step 1: Add a limit-based encoding helper.** Add a comment-documented convention at the top of the file: `:OnlyA` RE binding → set its inactive `K_I = 1e10` (× the active `K`); `:OnlyA` SS step → set its inactive rate `k_I = 0` (dead inactive cycle, `N_I = 0` — the only thermodynamically legal SS limit, verified in the spec). Every `*_onlyA_flux` network in this file must build the **full** two-conformation network (both S-binding and catalysis present in I) and realize `:OnlyA` by these values, NOT by omitting edges.

- [ ] **Step 2: Rewrite one existing gate (`uni_onlyA_flux`) to the limit encoding and confirm it still matches.** Keep the derived-vs-oracle assertion; it must still pass (for a valid `S:OnlyA` mechanism the limit equals deletion, so the number is unchanged). Run:
```
cd /home/denis.linux/.julia/dev/EnzymeRates
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/allosteric_ground_truth.jl")' 2>&1 | tail -20
```
Expected: the `OnlyA MWC derivation matches...` testset still green (≤1e-4).

- [ ] **Step 3: Do the same for `multi_onlyA_flux` and `metab_dfree_onlyA_flux`.** Same invariant: limit encoding, numbers unchanged.

- [ ] **Step 4: Add a ping-pong value gate.** Port the prototype's `oracle_Eflip` from `scratchpad/redesign_pingpong.jl` (an 8-form E-flip two-conformation ping-pong, `:NonequalAI` catalysis) as a testset that asserts the shipped `rate_equation` for the ping-pong `@allosteric_mechanism` (the one in `scratchpad/compare_pingpong.jl`, verified 3e-16) matches it. Map params by meaning as documented there. This closes the harness's long-deferred ping-pong gap and pins the case the rewrite must not break.

- [ ] **Step 5: Add an `N_I ≠ 0` multi-protomer gate.** Build an explicit 2-protomer concerted-MWC network with an **active** inactive conformation (`:NonequalAI` catalysis, `N_I ≠ 0`) and assert `rate_equation` at `catalytic_multiplicity: 2` matches it — this is the `^n` cross term the prototype left untested.

- [ ] **Step 6: Run and commit.**
```
julia --project -e 'using TestEnv; TestEnv.activate(); include("test/allosteric_ground_truth.jl")' 2>&1 | tail -20
git add test/allosteric_ground_truth.jl
git commit -m "Strengthen MWC oracle: limit-based :OnlyA, ping-pong value gate, N_I!=0 n=2 gate"
```
Expected: all gates green against the *current shipped* derivation (it is already correct). These gates now define the target for Tasks 2-7.

---

### Task 2: Derive the inactive state on the full (undeleted) graph

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_state_rate_polys` (`:1290`), `_state_mechanism` / `_state_allo_mechanism` usage.
- Test: `test/test_rate_eq_derivation.jl` (new testset).

**Interfaces:**
- Produces: `_state_rate_polys(am, :I)` derives on the full I-graph with I-scoped parameter names, returning `(num_I, den_I, d_free_I)` on the **undeleted** topology.

- [ ] **Step 1: Write the failing test.** For an `S:OnlyA` uni-uni, assert the I-state poly derivation includes the S-bound form (before the limit) — i.e. the undeleted graph has the same form-count as the A-state. Concretely, assert `EnzymeRates.n_states` of the I-derivation graph equals the A-graph's (the deletion currently makes it smaller).
```julia
@testset "I-state derives on the full graph (undeleted)" begin
    m = @allosteric_mechanism begin
        substrates: S ; products: P ; catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S) :: OnlyA ; E(S) <--> E(P) :: OnlyA ; E + P ⇌ E(P) :: EqualAI
        end
    end
    am = EnzymeRates.AllostericMechanism(m)
    a_forms = # count forms in _state_mechanism(am, :A)
    i_forms = # count forms in the NEW full-graph I derivation
    @test i_forms == a_forms
end
```
Fill the form-count via the accessors used in `_compute_re_groups` (read `_raw_symbolic_rate_polys`).

- [ ] **Step 2: Run it and see it fail** (the deleted I-graph is smaller).

- [ ] **Step 3: Implement.** Change the `:I` branch of `_state_rate_polys` to derive on the full mechanism with I-scoped step params (reuse `_state_step_params(am, :I)` for the naming, but the graph is the full `Mechanism`, not `_state_allo_mechanism(am, :I)`). Do NOT yet change the combine or apply the limit — the extra I-terms will be dropped by the limit in Task 5. Keep `d_free_I` returned.

- [ ] **Step 4: Run the test; confirm form-counts match.** Other allosteric tests WILL break (the combine now sees un-limited I-polys) — that is expected; Tasks 3-5 fix the combine. Record which break; do not fix them here.

- [ ] **Step 5: Commit.**
```
git commit -am "Derive the inactive MWC state on the full undeleted graph"
```

---

### Task 3: Solve the thermodynamic constraints on the full coupled system

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_combined_state_dependent_exprs` (`:1439`), `_dependent_param_exprs`.
- Test: `test/test_rate_eq_derivation.jl`.

**Interfaces:**
- Produces: `_combined_state_dependent_exprs(am)` solves the stacked A+I constraints on the **full** graph; `:OnlyA`/`:NonequalAI` I-params are genuine unknowns it resolves.

- [ ] **Step 1: Write the failing test.** For the `S:OnlyA`+`cat:OnlyA` family-A mechanism, assert the solve returns dependent expressions for the I-state reverse rates that reference the I-forward rate and `Keq` (a real Haldane row), rather than the deletion producing no I-constraint. Assert `_dep_graph_is_sound(typeof(m))` (reuse the helper at `test_rate_eq_derivation.jl:1232`).

- [ ] **Step 2: Run it; see how the current (deleted) solve differs.**

- [ ] **Step 3: Implement.** `_combined_state_dependent_exprs` already stacks A and I rows (`:1439`); ensure the I rows now come from the full I-graph (a consequence of Task 2) and that the solve treats the full I-parameter set. Confirm the priority/pivot machinery still resolves a unique dependent set.

- [ ] **Step 4: Run; confirm the dependent-param graph is sound.**

- [ ] **Step 5: Commit.**
```
git commit -am "Solve MWC thermodynamic constraints on the full two-state system"
```

---

### Task 4: The limit-the-constraints guard

**Files:**
- Create/modify: `src/rate_eq_derivation.jl` — a helper applied to the solved constraints; wire the check where the derivation begins.
- Test: `test/test_rate_eq_derivation.jl`.

**Interfaces:**
- Produces: a function that takes the solved constraint relations, applies the `:OnlyA` limits, and errors if any reduces to a contradiction (a lone `∞`, or `K = ∞` forced onto an `:EqualAI`/shared step).

- [ ] **Step 1: Write the failing test.** Assert the guard rejects the same set PR #70's `_onlya_haldane_violation` rejects, PLUS a ter-uni multi-cycle case #70's sign heuristic admits. Use the ter-uni witness from `scratchpad/` (the `EB→EAB` lone-`:OnlyA`-edge cube). Assert the limit-the-constraints check flags it while `_onlya_haldane_violation` returns `nothing`.

- [ ] **Step 2: Run; see the sign heuristic miss the ter-uni case.**

- [ ] **Step 3: Implement.** After the constraint solve, substitute the `:OnlyA` limits (`log K_I → +∞` for RE, the SS `k_I → 0`) into each solved relation and test that the infinite exponents cancel (a finite/`0=0` equality). Error with a message naming the offending cycle on any contradiction. This is exact (Stiemke feasibility); it supersedes the per-row heuristic for the derivation.

- [ ] **Step 4: Run; confirm it catches both the #70 set and the ter-uni case, and accepts all valid mechanisms (family A, TR, ping-pong).**

- [ ] **Step 5: Commit.**
```
git commit -am "Complete MWC consistency guard: limit the solved constraints"
```

---

### Task 5: Per-state normalization, plain combine, and the :OnlyA limit (the crux)

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_allosteric_num_den_exprs` (`:1662-1740`).
- Test: `test/allosteric_ground_truth.jl` (Task 1 gates are the acceptance test).

**Interfaces:**
- Produces: `_allosteric_num_den_exprs` returns the combined `(num, den)` Exprs via decoupled per-state normalization + plain combine + `:OnlyA` limit, with NO coupled cross-weight branch.

This is the load-bearing task. The recipe is **known-correct** (per-state normalization is algebraically the shipped cross-weight `d_free_I·Q_A + L·d_free_A·Q_I`); the work is expressing it decoupled and taking the limit. The Task 1 oracle is the definition of done — implement until every gate passes.

- [ ] **Step 1: Confirm the failing state.** Run the Task 1 gates; several fail now (Task 2 changed the I-polys, the old combine is inconsistent with them).

- [ ] **Step 2: Implement per-state normalization.** Normalize each state's `(num, den)` to its own free-enzyme weight `d_free`. Metabolite-free `d_free` (monomial): divide (reuse `_is_metabolite_free_monomial`/`_invert_monomial` as the fast path). Metabolite-bearing `d_free`: keep polynomial form via the cross-weight identity (multiply each state's contribution by the *other* state's weight) — this is what makes it a polynomial, not a rational. Do NOT branch on `d_free_A == d_free_I`; normalize each state unconditionally.

- [ ] **Step 3: Plain combine + `:OnlyA` limit.** Combine the normalized states as `(N_A·D_A^{n-1} + L·N_I·D_I^{n-1})/(D_A^n + L·D_I^n)`. Then apply the `:OnlyA` limit: drop every monomial carrying an `:OnlyA` `K_I`-inverse (RE) or `:OnlyA` `k_I` (SS). Verify the limit commutes with the `^n` power (it does — `N_I = 0` for a dead inactive kills the cross term at any `n`).

- [ ] **Step 4: Run the full Task 1 oracle.** Every gate — uni/multi-`:OnlyA`, metabolite-bearing-D, `:NonequalAI`, TR, ping-pong value, `N_I≠0` n=2 — must match to ≤1e-4 with `v=0` at equilibrium. Iterate against the oracle until green. If a case cannot be made to match, STOP and report — that is a real design problem, not a tuning issue.

- [ ] **Step 5: Confirm the non-allosteric path is byte-identical** (the golden; regenerate only if allosteric blocks changed).

- [ ] **Step 6: Commit.**
```
git commit -am "MWC combine via decoupled per-state normalization + :OnlyA limit"
```

---

### Task 6: Recompute `_kcat_forward` without the coupled cross-weight

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_kcat_forward` (`:941`).
- Test: `test/allosteric_ground_truth.jl` (the kcat assertions) + `test/test_rate_eq_derivation.jl` (kcat-consistency + perf).

- [ ] **Step 1: Confirm the failing state** — `_kcat_forward` still cross-weights; make it use the Task-5 normalized combine's saturating limit.

- [ ] **Step 2: Implement.** Recompute the saturating (substrate → ∞, products = 0) limit of the new normalized combine. Reuse the group-key machinery, but on the per-state-normalized polys, not the cross-weighted ones.

- [ ] **Step 3: Run.** The `isfinite(_kcat_forward(...))` gate assertions from Task 1 and the `kcat consistent with rate_equation` testset must pass. **The perf gate (`test_rate_equation_performance`, 0 alloc, <120 ns) must stay green** — `_kcat_forward` runs post-optimization but the `rate_equation` body must not regress.
```
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_rate_eq_derivation.jl")' 2>&1 | tail -20
```

- [ ] **Step 4: Commit.**
```
git commit -am "Recompute kcat from the normalized combine"
```

---

### Task 7: Delete the graph-deletion machinery

**Files:**
- Modify: `src/rate_eq_derivation.jl`.
- Test: per-file suite.

Only now that the new path is green — delete the dead code, so nothing is removed while still referenced.

- [ ] **Step 1: Remove `_state_allo_mechanism`'s pruning** (`:1214`) and the `d_free_A == d_free_I` three-way branch selection in `_allosteric_num_den_exprs`. Grep for every caller first; the `:I` state now derives on the full graph, so the pruning function may become entirely unused.

- [ ] **Step 2: Remove PR #70's per-row sign guard from the constructor** IF the Task-4 limit-the-constraints check is wired to run at construction; otherwise keep `_onlya_haldane_violation` as the cheap enumeration pre-filter and note that in a comment. (Design sub-decision — default to KEEPING it as the pre-filter unless the derivation-time check is confirmed cheap enough for enumeration.)

- [ ] **Step 3: Run per-file suites** — `test_rate_eq_derivation.jl`, `test_mechanism_enumeration.jl`, `test_identify_rate_equation.jl`, `allosteric_ground_truth.jl` — each green in isolation.

- [ ] **Step 4: Report the LOC delta** (`git diff --stat main..HEAD -- src/rate_eq_derivation.jl`). Expect net-negative; the removed lines should be exactly the delete-then-solve surface.

- [ ] **Step 5: Commit.**
```
git commit -am "Delete graph-deletion machinery; solve-then-limit is the only path"
```

---

### Task 8: Migration and full-suite verification

**Files:**
- Modify: `test/reference/allosteric_golden_reference.txt`, any spec whose expected counts moved.

- [ ] **Step 1: Regenerate the allosteric golden.** Run `_allosteric_golden_lines()` and rewrite `test/reference/allosteric_golden_reference.txt`. Diff it: allosteric REDUCED_STRING/PARAMS blocks may change (to the correct-but-possibly-different-form equations); the non-allosteric blocks must be byte-identical. Any changed block must be justified by an oracle match — do not accept a change you cannot explain.

- [ ] **Step 2: Re-measure any `expected_n_*` counts** in `test/mechanism_definitions_for_test_enzyme_derivation.jl` for allosteric specs whose derivation changed; update to measured values (fitted-param counts should be unchanged — the limit drops the same params — so a change is a finding to explain).

- [ ] **Step 3: Full suite (adequate-memory machine).** `julia --project -e 'using Pkg; Pkg.test()'` — 0 failures. In a memory-limited sandbox, verify per-file instead and note that the monolithic run needs the real machine.

- [ ] **Step 4: Perf gate green** (0 alloc, <120 ns).

- [ ] **Step 5: Update the docs** — `docs/src/deriving/mwc_allostery.md` num/den section, if the rendered form changed; `docs/src/identify/enumeration_engine.md` if the guard moved.

- [ ] **Step 6: Commit.**
```
git commit -am "Regenerate goldens and finalize solve-then-limit derivation"
```

---

## Self-review

- **Spec coverage:** motivation/robustness (Tasks 2-3, 7 remove the constraint-losing surface); solve-then-limit (Tasks 2-5); per-state normalization (Task 5, with the metabolite-bearing recipe); guard = limit-the-constraints (Task 4); oracle = limit-based (Task 1); `N_I≠0` n≥2 (Tasks 1, 5); code deletion moderate (Task 7 + LOC report); perf contract (every task). Covered.
- **The crux (Task 5)** is deliberately "iterate against the oracle" rather than exact code, because the decoupled symbolic re-expression is genuine implementation work; the recipe (cross-weight identity) and the acceptance test (Task 1 gates) are both pinned, so it is not a placeholder — it is a well-gated open task.
- **Ordering:** gate-first (Task 1), then derivation inside-out (full graph → solve → guard → combine → kcat), delete dead code last (Task 7), migrate last (Task 8). Nothing is deleted while referenced.
