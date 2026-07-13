# Allosteric MWC Free-Enzyme Normalization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Fix the pre-existing allosteric MWC denominator bug (a spurious catalytic rate constant leaks into the `L`-term when the inactive graph fragments into more segments than the active), validated against a first-principles `n=1` mass-action ground truth.

**Architecture:** Cross-weight each conformation's King–Altman terms by the *other* conformation's free-enzyme spanning-tree weight (`den = D_I^n·Q_A^n + L·D_A^n·Q_I^n`) in `_allosteric_num_den_exprs` and `_kcat_forward`. Gate the fix on an `n=1` two-conformation mass-action steady-state solver built specifically for this validation. Ship on the current branch together with the multi-`:OnlyA` enumeration work.

**Tech Stack:** Julia, EnzymeRates.jl. `julia --project`.

## Global Constraints

- **The fix is accepted ONLY when the fixed `rate_equation` equals the `n=1` mass-action ground truth** for the gate mechanisms. Dimensional homogeneity and kcat-rescaling are necessary, not sufficient.
- **`rate_equation` performance is non-negotiable** (`allocs == 0`, `t < 120e-9` for every `MECHANISM_TEST_SPECS` entry). The fix multiplies polynomial factors into num/den; if it makes `rate_equation` allocate or exceed 120 ns, **STOP and flag Denis**.
- **STOP conditions** (report, do not force): the fix cannot be made to match the ground truth; the ground-truth harness cannot pass its own self-validation (L=0 / all-`:EqualAI`); the blast radius is far larger than "allosteric mechanisms with a form-disconnecting `:OnlyA`/`:OnlyI` ligand" (e.g. a non-allosteric or NonequalAI golden changes — that would signal the fix touched the wrong path).
- **The fix must live only in the MWC assembly** (`_allosteric_num_den_exprs`, `_kcat_forward`). The non-allosteric `N/Q` path must be byte-identical after the change (there the `[1/time]^(G−1)` factor cancels; renormalizing it would be a bug).
- **Reference prototype:** `scratchpad/gt.jl` (this session) is a working `n=1` ground truth for the uni-uni exemplar — 6-species linear steady-state, confirms cross-weighting == ground truth, inactive flux == 0. Task 1 generalizes its solver.
- Style: 92-char lines, 4-space indent; `# ABOUTME:` header on new files. Commit trailers on every commit:
  ```
  Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01VFCCyUsAM2fdaXLq5Pdtw1
  ```
- Branch: `catalytic-onlya-promote-move` (already checked out; enumeration + spec commits already present).

---

### Task 1: `n=1` two-conformation mass-action ground-truth harness (the gate)

**Files:**
- Create: `test/allosteric_ground_truth.jl` — the solver + explicit gate-mechanism networks + self-validation.

**Interfaces:**
- Produces: `mwc_ground_truth_flux(network, params)::Float64` — pseudo-first-order steady-state net catalytic flux for an explicit two-conformation network. A `network` is a list of species names and a list of directed edges `(from, to, rate_expr)` where `rate_expr` is a closure of `params`. Also produces the explicit networks `gt_uni_onlyA()`, `gt_biuni_multi_onlyA()`, `gt_uni_onlyI()`.
- Consumes: `LinearAlgebra`.

- [ ] **Step 1: Write the solver + the uni-uni `:OnlyA` network + self-validation tests.**

Generalize `scratchpad/gt.jl`. The solver builds the rate matrix `M` (edge `a→b` at rate `r`: `M[b,a]+=r; M[a,a]-=r`), replaces one row with the conservation `Σc=E_total`, solves `M c = e`, and returns the net active catalytic flux `Σ (k_fwd·c[reactant] − k_rev·c[product])` over active catalytic edges. Fast RE bindings and flips use a large separation constant `FAST=1e7`; catalysis is `O(1)`. Detailed-balance flip ratio `[X_I]/[X_A] = L·∏(K_A_i/K_I_i)` (free enzyme `L`; a state bearing an `:OnlyA` ligand → `0`, i.e. no flip edge; dead-end forms reached only by catalysis have no binding/flip edge).

```julia
# ABOUTME: n=1 two-conformation mass-action ground truth for allosteric MWC rate equations.
# ABOUTME: Small explicit networks solved as a linear steady state; the acceptance gate for the normalization fix.
using LinearAlgebra

"Net active catalytic flux at steady state for an explicit two-conformation network.
`species` is a Vector{Symbol}; `edges` is Vector of (from::Symbol,to::Symbol,rate::Float64);
`cat_edges` is Vector of (reactant::Symbol,product::Symbol,kfwd::Float64,krev::Float64) whose
NET flux (over ACTIVE-conformation edges only) is the reaction velocity."
function mwc_ground_truth_flux(species, edges, cat_edges_active, Etot)
    n = length(species); idx = Dict(s=>i for (i,s) in enumerate(species))
    M = zeros(n,n)
    for (a,b,r) in edges; M[idx[b],idx[a]] += r; M[idx[a],idx[a]] -= r; end
    A = copy(M); A[1,:] .= 1.0; rhs = zeros(n); rhs[1] = Etot
    c = A \ rhs
    sum(kf*c[idx[r]] - kr*c[idx[p]] for (r,p,kf,kr) in cat_edges_active)
end
```

Then the uni-uni `:OnlyA` network builder (mirrors `gt.jl`: species `E_A,ES_A,EP_A,E_I,EP_I,ES_I`; S `:OnlyA`, catalysis `:EqualAI`, P `:EqualAI`; `ES_I` dead-end via catalysis only), and self-validation tests:

```julia
@testset "ground-truth harness self-validation" begin
    KA,KP,k,L,Keq,S,P = 1.3,0.9,2.1,0.7,3.0,1.1,0.6
    # (a) L=0 : all-active, equals the single-conformation (non-allosteric) uni-uni rate
    f0 = uni_onlyA_flux(KA,KP,k,L=0.0,Keq=Keq,S=S,P=P)
    kr = k*KP/(Keq*KA); nonallo = (k*S/KA - kr*P/KP)/(1 + S/KA + P/KP)
    @test isapprox(f0, nonallo; rtol=1e-4)
    # (b) all-:EqualAI (K_I=K_A everywhere) : L-independent, equals the base rate
    fa = uni_equalAI_flux(KA,KP,k,L,Keq,S,P)
    @test isapprox(fa, nonallo; rtol=1e-4)
    @test isapprox(fa, uni_equalAI_flux(KA,KP,k,5.0,Keq,S,P); rtol=1e-4)   # L-independent
end
```

- [ ] **Step 2: Run the self-validation; confirm it passes.**

Run: `julia --project test/allosteric_ground_truth.jl`
Expected: PASS. If either sanity fails, the harness model is wrong — **STOP and diagnose; do not proceed to gate the fix with an unvalidated harness.**

- [ ] **Step 3: Add the `:OnlyA` gate assertion against the CURRENT (buggy) code — RED.**

Build the compiled uni-uni `:OnlyA` mechanism (`@allosteric_mechanism`, `catalytic_multiplicity: 1`), evaluate `EnzymeRates.rate_equation` at random params, and assert it equals `uni_onlyA_flux`. With the current code this must **FAIL** (the bug). Capture it — this is the RED the fix will turn green.

```julia
@testset "current OnlyA derivation vs ground truth (RED until fixed)" begin
    # build am, compile, map fitted_params -> (KA,KP,k,L,Keq); random trials
    # @test isapprox(rate_equation(cem, concs, prm), uni_onlyA_flux(...); rtol=1e-4)
end
```

- [ ] **Step 4: Commit** (`test/allosteric_ground_truth.jl` only). Message: "Add n=1 MWC mass-action ground-truth harness (self-validated)".

---

### Task 2: Surface `D[g_free]` (free-enzyme spanning-tree weight) per state

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_raw_symbolic_rate_polys` (`:331`, the `Mechanism` method that computes the `D` array at `:365`) and `_state_rate_polys` (`:1248`).

**Interfaces:**
- Produces: `_state_rate_polys(am, state)` returns `(num, den, d_free)` where `d_free` is the polynomial `D[g_free]` — the spanning-tree weight of the segment holding the free resting enzyme (empty `bound`, empty `residual`). For a single-segment state `d_free = poly_one()`.
- Consumes: the existing `D` array (`:365`), `groups`, `form_to_group`, and `_free_enz_set` (identifies the free resting enzyme forms).

- [ ] **Step 1: Write a failing test for `D[g_free]` on the exemplar.**

In `test/test_rate_eq_derivation.jl` (near the allosteric derivation tests), assert the free-enzyme weight for the uni-uni exemplar: active state `d_free == 1` (single segment), inactive state `d_free == k_ES_to_EP` (the catalytic forward rate — the segment holding `E_I` is spanned through the catalytic edge). Use whatever accessor Task 2 exposes; render the poly to compare (e.g. via `_poly_to_expr` or by evaluating at unit params).

- [ ] **Step 2: Run it; confirm it fails** (accessor/return value not present yet).

- [ ] **Step 3: Implement.**

In `_raw_symbolic_rate_polys(mech::Mechanism, ...)`: after the `D` array is built (`:365`), identify `g_free` = the group `form_to_group[f]` for the free resting-enzyme form `f` (from `_free_enz_set(mech)` / the form with empty bound + residual), and return `D[g_free]` alongside `num, den`. Thread the extra return value through `_state_rate_polys` (`:1248`). Update the other `_raw_symbolic_rate_polys` caller(s) if the tuple arity changes (grep for callers; the non-allosteric path may call it and must keep working — return the extra value but ignore it there).

- [ ] **Step 4: Run; confirm pass** (`d_free` values correct: active 1, inactive `k_ES_to_EP`).

- [ ] **Step 5: Commit** (`src/rate_eq_derivation.jl`, `test/test_rate_eq_derivation.jl`). Message: "Surface free-enzyme segment weight D[g_free] per allosteric state".

---

### Task 3: Cross-weight the MWC combination (the fix)

**Files:**
- Modify: `src/rate_eq_derivation.jl` — `_allosteric_num_den_exprs` (`:1595`, uses `_mwc_combine` at `:1647`) and `_kcat_forward` (`:881`, uses `_mwc_power_pair`/`_mwc_combine` at `:1002`–`:1051`).
- Test: `test/allosteric_ground_truth.jl` (flip the Task-1 RED assertions to GREEN).

**Interfaces:**
- Consumes: `d_free` from Task 2; `_power_expr` (`:1521`), `_mwc_combine` (`:1530`).

- [ ] **Step 1: Turn the Task-1 `:OnlyA` gate assertion into the RED/GREEN target** and add the multi-`:OnlyA` bi-uni and `:OnlyI` gate assertions (all currently failing).

- [ ] **Step 2: Confirm RED** (current derivation ≠ ground truth for `:OnlyA`, multi-`:OnlyA`, `:OnlyI`).

- [ ] **Step 3: Implement the cross-weighting.**

In `_allosteric_num_den_exprs`, with `CatN = n`, `D_A`/`D_I` the free-enzyme weights (Task 2) rendered to Exprs:
- Denominator: `den_A_term = D_I^n · (Q_A^n · reg_A)`, `den_I_term = D_A^n · (Q_I^n · reg_I)`; combine `den_A_term + L·den_I_term`.
- Numerator: `num_A_term = D_I^n · (N_A · Q_A^(n-1) · reg_A)`; `num_I_term = D_A^n · (N_I · Q_I^(n-1) · reg_I)` — **dropped entirely when `N_I = 0`** (the existing `isempty(num_i_poly)` branch: keep only the active term).
- Use `_power_expr(D_expr, CatN)` for `D_A^n`/`D_I^n`. The `D_A^n·D_I^n` common factor cancels in `num/den`; do NOT attempt to cancel it symbolically — leave both as polynomial factors (the reduction/simplification the engine already runs handles it, and the ground truth checks numeric equality regardless).

Apply the **same** cross-weighting in `_kcat_forward` so kcat stays consistent with `rate_equation`.

- [ ] **Step 4: Confirm GREEN** — `rate_equation` == ground truth for uni `:OnlyA`, multi-`:OnlyA` bi-uni, and `:OnlyI`, across random trials (rtol 1e-4). If any cannot be made to pass, **STOP and report** with the mismatch.

- [ ] **Step 5: Perf check.** Run `test_rate_equation_performance` on the multi-`:OnlyA` spec: `allocs == 0`, `t < 120e-9`. If it regresses, **STOP and flag Denis** (non-negotiable contract).

- [ ] **Step 6: Commit** (`src/rate_eq_derivation.jl`, `test/allosteric_ground_truth.jl`). Message: "Cross-weight MWC free-enzyme normalization (fixes L-term dimensional leak)".

---

### Task 4: Dimensional guard, ground-truth gate specs, golden re-baseline, full suite

**Files:**
- Modify: `test/test_rate_eq_derivation.jl` — dimensional-homogeneity guard over allosteric specs; wire the ground-truth harness as gate assertions.
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` — re-baselined expected counts; the multi-`:OnlyA` spec (from the enumeration branch) now derives correctly.
- Regenerate: `test/reference/allosteric_golden_reference.txt`.

- [ ] **Step 1: Dimensional-homogeneity guard.** For each allosteric spec, assert: rate constants ×τ → `v ×τ`; concs and dissociation-K's ×λ together → `v` unchanged; `E_total ×μ` → `v ×μ`. This must now pass for every allosteric spec (it currently fails for the buggy multi-`:OnlyA` shape).

- [ ] **Step 2: Re-derive the affected goldens and JUSTIFY each change.** Run the full derivation suite; for every allosteric golden that changed (`allosteric_golden_reference.txt`, LDH i-state `expected_n_independent_params`, the enumeration branch's multi-`:OnlyA` spec counts, D1), confirm the change is a correctness re-baseline (the new equation matches the ground truth or is dimensionally clean where the old was not). Regenerate `allosteric_golden_reference.txt`. **If a NON-allosteric or a NonequalAI-only golden changed, STOP** — the fix touched the wrong path.

- [ ] **Step 3: Update the multi-`:OnlyA` `MECHANISM_TEST_SPECS` golden counts** (Task-3 of the enumeration work) to the corrected derivation, and confirm the kcat-rescaling sub-test now passes for it (it was the original failure).

- [ ] **Step 4: Full suite.** `julia --project -e 'using Pkg; Pkg.test()'` — all green. Confirm especially: perf contract, allosteric golden/collapse, the ground-truth gate assertions, kcat-rescaling, and the enumeration reachability tests.

- [ ] **Step 5: Commit** the re-baselined goldens + guards. Message: "Re-baseline allosteric goldens + dimensional guard under corrected normalization".

---

### Final verification

- [ ] Full `Pkg.test()` green (foreground, monitored — never background+yield).
- [ ] Ground-truth gate passes for uni `:OnlyA`, multi-`:OnlyA`, `:OnlyI` (and ping-pong / `n=2` if the harness was extended to them; if not, note in the report that those remain validated only by dimensional homogeneity + the cross-weighting's polynomial-preservation argument).
- [ ] `git diff` on the non-allosteric derivation path is empty (fix is MWC-assembly-only).
- [ ] Branch ready for Denis's review.
