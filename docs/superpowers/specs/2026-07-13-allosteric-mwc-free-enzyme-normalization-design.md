# Allosteric MWC free-enzyme normalization: fix the cross-state denominator weighting

**Date:** 2026-07-13
**Status:** Fix and acceptance gate both validated by a first-principles `n=1` mass-action
ground truth (this session). Pending plan + implementation.

## Context

The allosteric rate-equation derivation combines the two MWC conformations into

```
den = Q_A^n + L · Q_I^n            (schematically; n = catalytic_multiplicity)
```

where `Q_A`, `Q_I` are the King–Altman denominators of the active and inactive
conformational sub-mechanisms and `L` is the conformational constant. For a mechanism
whose inactive graph fragments into more rapid-equilibrium segments than its active
graph, this combination is **wrong**: it leaves a spurious catalytic rate constant in the
`L`-term.

Discovered while adding the multi-`:OnlyA` enumeration move (branch
`catalytic-onlya-promote-move`). The minimal exemplar is an `n=1` uni-uni — S binding
`:OnlyA`, catalysis `:EqualAI`, P binding `:EqualAI` — which derives to

```
den = 1 + S/K_A + P/K_P + L·(k + k_r·P/K_P + k·P/K_P)     ← bare rate constants k, k_r
num = k·S/K_A − k_r·P/K_P
```

The active partition `1 + S/K_A + P/K_P` is dimensionless; the inactive `L`-term carries
bare catalytic rate constants `k ≡ k_ES_to_EP`, `k_r ≡ k_EP_to_ES`. That makes the rate
non-homogeneous in the steady-state rate constants and breaks `rate_equation`'s
kcat-rescaling contract (`v_norm/v_orig ≈ 1/kcat` fails). A single-`:OnlyA` bi-uni of the
same shape — reachable before the enumeration move — fails identically, so the defect is
**pre-existing**, not introduced by the move; the move only makes it far more reachable.

## Root cause

The King–Altman denominator of a sub-mechanism is `Q = Σ_g σ_g · D[g]`, summing over
rapid-equilibrium segments `g`: `σ_g` is the within-segment binding partition
(dimensionless), and `D[g]` is the Matrix–Tree spanning-tree weight of segment `g` (a
product of `G−1` steady-state edge weights, `G` = segment count). Every steady-state edge
weight is `[1/time]`, so **`Q` is dimensionally homogeneous of degree `[1/time]^(G−1)`**.

Each state's `Q` is therefore expressed **relative to its own free-enzyme segment's
spanning-tree weight `D[g_free]`**, not to free enzyme = 1:

- Active state of the exemplar: one RE segment (`G_A = 1`), so `D_A = 1` and `Q_A` is
  dimensionless.
- Inactive state: dropping the `:OnlyA` S binding disconnects the S-bound catalytic form
  `E(S)_I` — reachable only through the steady-state catalytic edge — into its own
  segment (`G_I = 2`), so `D_I = k` and `Q_I = k·(dimensionless partition)`.

`Q_A^n + L·Q_I^n` then adds a `[1/time]^0` term to a `[1/time]^n` term. With a
dimensionless `L` this cannot be homogeneous, and the derivation "pays for" the mismatch
by leaving the catalytic rate constant in the `L`-term. **The defect is purely that `Q_A`
and `Q_I` are combined on inconsistent normalization bases.**

### What the dead-end is, and what it is *not*

The inactive S-bound form `E(S)_I` is a **populated dead-end**: with S `:OnlyA` there is
no inactive S-binding step, so `E(S)_I` is fed only by reverse catalysis from `E(P)_I` and
its single graph edge is that catalytic one. A pendant node carries **zero net steady-state
flux**, so:

- It contributes **nothing to the numerator** — the inactive state runs no productive
  cycle. The code's `:OnlyA` numerator `k·S/K_A − k_r·P/K_P` (no `L` term) is *correct*;
  `N_I = 0`.
- It **does** contribute to the denominator (it holds real enzyme), pinned by the
  catalytic balance `[E(S)_I] = (k_r/k)·[E(P)_I]`, so the physical inactive partition is
  `Q_I = 1 + P/K_P + (k_r/k)·P/K_P`.

**`:OnlyA` is not the `K_I→∞` limit of `:NonequalAI`.** That limit is singular: it keeps
the inactive S-binding *step* present (just weak), so `E(S)_I` can still release S and
completes a flux-carrying inactive cycle — different physics from `:OnlyA`, which removes
the step entirely. The `:NonequalAI`-limit form (with a spurious `−L·k_r·P/K_P` numerator
term and a pure-binding `Q_I = 1 + P/K_P`) is therefore the **wrong** target and must not
be used as the gate. (Classical pure-binding `Q_I` is recovered only when *both* the S
binding and the catalysis are `:OnlyA`, which makes `E(S)_I` wholly unreachable.)

This was confirmed by the `n=1` mass-action ground truth below.

The standalone non-allosteric rate is unaffected: there `v = E_total·N/Q` and the
`[1/time]^(G−1)` factor cancels between numerator and denominator. Only the MWC
combination — two partitions added against a dimensionless `L` — needs a common basis.

## The fix

Normalize each conformation's partition to its own free enzyme,
`P_state = Q_state / D[g_free_state]` (the physical `[total enzyme]/[free E]`), and combine
as `P_A^n + L·P_I^n`. Dividing by `D[g_free]` would produce rational functions when
`D[g_free]` carries a metabolite (ping-pong), so clear the common `D_A^n·D_I^n` and keep
everything polynomial by **cross-weighting** each conformation's terms with the *other*
conformation's free-enzyme weight:

```
den = D_I^n · Q_A^n            +  L · D_A^n · Q_I^n
num = D_I^n · N_A · Q_A^(n-1)  +  L · D_A^n · N_I · Q_I^(n-1)
```

- `Q_A`, `Q_I` — per-state King–Altman denominators (as today).
- `N_A`, `N_I` — per-state King–Altman numerators (`N_I = 0` for a dead inactive state,
  already detected by `_i_state_num_zero`; then the `L·num_I` term drops).
- `D_A ≡ D[g_free_A]`, `D_I ≡ D[g_free_I]` — spanning-tree weights of the segment holding
  the free resting enzyme `E` (empty `bound`, empty `residual`), read from the `D` array
  already computed in `_raw_symbolic_rate_polys` (`rate_eq_derivation.jl:365`); identify
  `g_free` via the existing `_free_enz_set`.
- Regulatory-site factors multiply each state's terms exactly as today.

The `D_A^n`/`D_I^n` common factor cancels in `num/den`, so the result is the physical,
dimensionally homogeneous `P_A^n + L·P_I^n` with a bare conformational `L`.

### Validation — this fix is confirmed, not hypothesized

An independent `n=1` two-conformation mass-action steady-state solve (6 species; S `:OnlyA`
dead-end inactive graph, `:EqualAI` catalysis, `:EqualAI` product; fast RE bindings, fast
detailed-balance conformational flips, slow catalysis) reproduces the **cross-weighted**
rate to 4+ digits across random parameters, with inactive net flux `= 0` (dead-end
confirmed). It rejects both the current buggy derivation and the `:NonequalAI`-limit form.
So the cross-weighting fix is the correct answer, `N_I = 0` for the dead inactive state,
and `E(S)_I` belongs in `Q_I`.

## Why `n=1` is a sufficient gate

`_allosteric_num_den_exprs` computes the per-state polynomials `N_A, Q_A, N_I, Q_I`
**once**; `catalytic_multiplicity` (`n`) enters only as the exponent afterward
(`Q^n`, `Q^(n-1)`). The cross-state normalization bug — and its fix — live entirely in
those per-state polys and their combination, which are `n`-independent. So the fix is fully
exercised at `n=1`, where the two-conformation mass-action system is small (single-digit
species) and exactly solvable. The oligomer `Q^n` structure is the standard concerted MWC
form and is checked separately at `n=2` (a dimer, still small). **The full tetramer is
never modeled** — that is the combinatorial explosion the ODE approach rightly avoids, and
it is not where the bug lives.

## Acceptance gate — `n=1` two-conformation mass-action ground truth

The existing `ode_steady_state_flux` cannot serve as the gate: for an allosteric mechanism
it builds the ODE from the **base single-conformation** reaction set (`enzyme_forms`
returns the base forms, no inactive conformation, no `L`), so it models the wrong system —
which is why every allosteric spec sets `run_ode_test=false`. A new ground truth is
required.

Build an **`n=1` two-conformation mass-action ground truth** (an "unfold"):

- Duplicate every catalytic form into an active (`_A`) and inactive (`_I`) copy.
- Per-conformation binding steps set by the tag: `:EqualAI` binds both (K_A = K_I);
  `:NonequalAI` binds both (K_A ≠ K_I); `:OnlyA` binds active only (no inactive binding
  step → inactive forms holding that ligand exist only if reachable by catalysis, and are
  dead-ends); `:OnlyI` binds inactive only.
- Conformational flips between corresponding `_A`/`_I` forms, fast, with detailed-balance
  ratio `[X_I]/[X_A] = L·∏(K_A_i/K_I_i)` over the ligands in state `X` (so free enzyme is
  `L`, an `:OnlyA`-ligand-bearing state is `0` — no flip).
- Catalytic steps per conformation set by the catalytic tag; reverse rate from the
  existing Haldane relation.
- Solve the pseudo-first-order steady state (linear solve, fast/slow rate separation;
  exact, non-stiff), and read the net catalytic flux.

**Self-validation of the harness (before it may gate anything):**
- `L = 0` → flux equals the **non-allosteric active-mechanism** `rate_equation` (an
  ODE-validated path).
- all-`:EqualAI` → flux equals the non-allosteric base `rate_equation`, independent of `L`.

Only once those hold does the harness gate the fix.

**The fix is accepted only when the fixed `rate_equation` equals this ground truth** for:
1. an allosteric multi-`:OnlyA` mechanism (dead inactive state, `G_I > G_A`),
2. an allosteric `:OnlyI` mechanism (the mirror case),
3. an allosteric ping-pong mechanism (concentration-bearing `D[g_free]` — the case that
   forbids dividing by `D[g_free]` and requires cross-weighting), and
4. an `n=2` version of (1) (exercises the concerted `Q^n` power).

Dimensional homogeneity and the kcat-rescaling contract are **necessary but not
sufficient**; ground-truth equality is the gate.

## Secondary guard — dimensional homogeneity

Add a cheap, permanent, per-allosteric-spec check (runs on every spec, unlike the `n=1`
ground truth): scaling all steady-state rate constants by `τ` must scale `v` by `τ`;
scaling all concentrations and all dissociation constants by `λ` together must leave `v`
unchanged; scaling `E_total` by `μ` must scale `v` by `μ`. This catches the whole bug class
on every mechanism, including ones too large to unfold. It is a regression guard, not the
acceptance gate.

## Scope / blast radius

The change alters every allosteric mechanism whose inactive (or active) graph fragments
into a different segment count than its counterpart — any mechanism where an
`:OnlyA`/`:OnlyI` binding disconnects a downstream catalytic form. Verify empirically which
`MECHANISM_TEST_SPECS` are affected (a mechanism whose `:OnlyA` ligand does *not*
disconnect a form — e.g. one already kept connected by other bindings, or by SS bindings
like LDH's `kon`/`koff` — is unchanged and still passes kcat-rescaling today). Expect
possible updates to: `test/reference/allosteric_golden_reference.txt` (regenerated), LDH
i-state `MECHANISM_TEST_SPECS` expected param counts if the fix removes a spuriously
independent constant, the multi-`:OnlyA` spec and D1 golden on the enumeration branch, and
any selection golden that depended on the old (incorrect) equations. Each re-baseline must
be justified by the ground-truth gate, not accepted merely because it changed.

## Files

- `src/rate_eq_derivation.jl` — expose `D[g_free]` per state; cross-weight in
  `_allosteric_num_den_exprs` and the `_kcat_forward` combine.
- `test/` — a new `n=1` two-conformation mass-action ground-truth harness (self-validated),
  ODE-gate assertions for the four gate mechanisms, and the dimensional-homogeneity guard.
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` — gate specs (multi-`:OnlyA`,
  `:OnlyI`, ping-pong, `n=2`); re-baselined expected counts.
- `test/reference/allosteric_golden_reference.txt` — regenerated.

## Non-goals

- Changing the King–Altman engine or the non-allosteric derivation.
- Pruning the inactive dead-end complex — it is correct enzyme mass; the bug is
  normalization, not membership.
- Modeling the full oligomer ODE — `n=1` (plus an `n=2` power check) is sufficient because
  the bug is `n`-independent.
- Estimating `Keq` from data (always user-supplied).

## Sequencing

Denis's decision: land this derivation fix **together with** the multi-`:OnlyA` enumeration
work on the same branch (`catalytic-onlya-promote-move`), since they are coupled — the
enumeration branch's failing multi-`:OnlyA` golden/perf spec and D1 golden re-baseline
correctly only once this fix is in. Implement the fix (gate first), then re-baseline all
affected goldens, then the full suite is green on the combined branch.
