# Allosteric MWC free-enzyme normalization: fix the cross-state denominator weighting

**Date:** 2026-07-13
**Status:** Design approved (brainstorm with Denis), pending plan + implementation.

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
`catalytic-onlya-promote-move`). An ordered bi-uni with two `:OnlyA` substrate bindings
derives to

```
den = (1 + A/K_A + P/K_P + A·B/(K_A·K_B))^2
    + L · (k_EAB_to_EP + k_EAB_to_EP·P/K_P + k_EP_to_EAB·P/K_P)^2
```

The active partition is dimensionless; the inactive partition carries a bare
`k_EAB_to_EP` (a catalytic rate constant). That extra factor makes the rate
non-homogeneous in the steady-state rate constants and breaks `rate_equation`'s
kcat-rescaling contract (`v_norm/v_orig ≈ 1/kcat` fails). A single-`:OnlyA` bi-uni of the
same topology — reachable before the enumeration move — fails identically, so the defect
is **pre-existing**, not introduced by the move; the move only makes it far more
reachable.

## Root cause

The King–Altman denominator of a sub-mechanism is `Q = Σ_g σ_g · D[g]`, summing over
rapid-equilibrium segments `g`: `σ_g` is the within-segment binding partition
(dimensionless), and `D[g]` is the Matrix–Tree spanning-tree weight of segment `g` (a
product of `G−1` steady-state edge weights, where `G` is the segment count). Every
steady-state edge weight is `[1/time]`, so **`Q` is dimensionally homogeneous of degree
`[1/time]^(G−1)`** (confirmed empirically by a concentration-unit-scaling test on the
Segel Ping Pong Bi Bi derivation: the rate scales exactly with the concentration unit
once bimolecular constants are identified by whether their direction consumes a
metabolite, not by the `kon`/`koff` name).

Each state's `Q` is therefore expressed **relative to its own free-enzyme segment's
spanning-tree weight `D[g_free]`**, not to free enzyme = 1:

- Active state of the exemplar: one RE segment (`G_A = 1`), so `D_A = 1` and `Q_A` is
  dimensionless.
- Inactive state: dropping the `:OnlyA` substrate bindings disconnects the
  substrate-bound catalytic form into its own segment, reached only through the
  steady-state catalytic edge (`G_I = 2`), so `D_I = k_EAB_to_EP` and
  `Q_I = k_EAB_to_EP · (dimensionless partition)`.

`Q_A^n + L·Q_I^n` then adds a `[1/time]^0` term to a `[1/time]^n` term. With a
dimensionless `L` this cannot be homogeneous, and the derivation "pays for" the mismatch
by leaving the `k_EAB_to_EP` factor in the `L`-term. **The inactive complex is a
legitimate dead-end species and belongs in the denominator; nothing needs pruning. The
defect is purely that `Q_A` and `Q_I` are combined on inconsistent normalization bases.**

The standalone non-allosteric rate is unaffected: there `v = E_total·N/Q` and the
`[1/time]^(G−1)` factor cancels between numerator and denominator, so a dimensionful `Q`
is fine. Only the MWC combination — which adds two partitions against a dimensionless
`L` — requires them on a common basis.

## The fix

Normalize each conformation's partition to its own free enzyme,
`P_state = Q_state / D[g_free_state]` (the physical partition `[total enzyme]/[free E]`),
and combine as `P_A^n + L·P_I^n`, giving the **physical** `L`. Dividing by `D[g_free]`
would produce rational functions when `D[g_free]` carries a metabolite (ping-pong), so
clear the common `D_A^n·D_I^n` and keep everything polynomial by **cross-weighting** each
conformation's terms with the *other* conformation's free-enzyme weight:

```
den = D_I^n · Q_A^n            +  L · D_A^n · Q_I^n
num = D_I^n · N_A · Q_A^(n-1)  +  L · D_A^n · N_I · Q_I^(n-1)
```

- `Q_A`, `Q_I` — per-state King–Altman denominators (as today).
- `N_A`, `N_I` — per-state King–Altman numerators (`N_I = 0` for a dead inactive state,
  already detected by `_i_state_num_zero`; then the `L·num_I` term drops).
- `D_A ≡ D[g_free_A]`, `D_I ≡ D[g_free_I]` — spanning-tree weights of the segment holding
  the free resting enzyme `E` (no bound metabolite, no residual), read from the `D` array
  already computed in `_raw_symbolic_rate_polys` (`rate_eq_derivation.jl:365`).
- Regulatory-site factors multiply each state's terms exactly as today.

### Properties

- **Polynomial.** `Q`, `N`, `D[g]` are all polynomials the engine already produces — no
  rational-function machinery. This is what the concentration-dependent-`D` (ping-pong)
  case needs and could not get by dividing.
- **Dimensionally homogeneous, `L` dimensionless.** Both denominator terms carry
  `[1/time]^(n(G_A+G_I−2))`; `v = E_total·num/den ~ conc/time`.
- **Physical `L`.** The cross-weighted form is the cleared `P_A^n + L·P_I^n`, so `L` is
  the bare conformational constant for every mechanism, ping-pong included; a
  concentration inside a ping-pong `D[g_free]` rides along as a polynomial factor and is
  never divided out.
- **Reduces to the reported-bug fix.** Multi-`:OnlyA`: `D_A = 1`, `D_I = k_EAB_to_EP`, so
  `den = k_EAB^n·Q_A^n + L·Q_I^n`, identical to `Q_A^n + L·(Q_I/k_EAB)^n` up to the
  `k_EAB^n` that cancels in `num/den`.
- **Non-allosteric path untouched.** The change lives only in
  `_allosteric_num_den_exprs` (the MWC assembly). Standalone `N/Q` is never renormalized.

## Scope / blast radius

The change alters every allosteric mechanism whose inactive (or active) graph fragments
into a different segment count than its counterpart — i.e. any mechanism where an
`:OnlyA`/`:OnlyI` binding disconnects a downstream catalytic form. That includes LDH and
PK, whose catalytic reactant holds an `:OnlyA` ligand, so their current allosteric
goldens encode the un-normalized denominator and **will change**. This is a correctness
re-baseline, not a regression. Expect updates to: the allosteric golden reference
(`test/reference/allosteric_golden_reference.txt`), the LDH i-state
`MECHANISM_TEST_SPECS` expected param counts if the fix removes a spuriously-independent
constant, and any selection golden that depended on the old (incorrect) equations.

## Acceptance gate — ODE ground truth

The fix is accepted only when the derived MWC `rate_equation` equals the mass-action
steady-state flux (`ode_steady_state_flux`, the validator that caught the
kcat-singularity defect) for:

1. an allosteric multi-`:OnlyA` mechanism (dead inactive state, `G_I > G_A`), and
2. an allosteric ping-pong mechanism (which also establishes whether ping-pong even
   presents a `G_I ≠ G_A` mismatch, and exercises a concentration-bearing `D[g_free]`).

Dimensional homogeneity and the kcat-rescaling contract are necessary checks but not
sufficient; ODE equality is the gate.

## Implementation notes

- **Surface `D[g_free_state]`.** `_state_rate_polys` / `_raw_symbolic_rate_polys` compute
  `D[g]` internally; expose the free-enzyme segment's weight. Identify `g_free` via the
  existing `_free_enz_set` (the free resting enzyme: empty `bound`, empty `residual`).
- **Cross-weight at assembly.** In `_allosteric_num_den_exprs`, multiply the active
  numerator/denominator terms by `D_I^n` and the inactive terms by `D_A^n` before the
  `_mwc_combine`. Reuse `_power_expr`.
- **Dead inactive state.** When `N_I = 0`, the numerator keeps only
  `D_I^n · N_A · Q_A^(n-1)`; the denominator keeps both terms (`Q_I` is still enzyme
  mass). This mirrors today's `isempty(num_i_poly)` branch.
- **Kcat path.** `_kcat_forward` combines the same per-state polys via
  `_mwc_power_pair`/`_mwc_combine`; apply the same cross-weighting so kcat stays
  consistent with `rate_equation` (the kcat-rescaling test compares them).

## Non-goals

- Changing the King–Altman engine or the non-allosteric derivation.
- Pruning the inactive dead-end complex — it is a correct species; the bug is
  normalization, not membership.
- Estimating `Keq` from data (always user-supplied).
- The `catalytic-onlya-promote-move` enumeration branch itself, which is correct; it is
  merely the change that surfaced this defect. That branch's failing multi-`:OnlyA`
  golden/perf spec depends on this fix and should land after it (or together, once the
  goldens are re-baselined).

## Files

- `src/rate_eq_derivation.jl` — expose `D[g_free]` per state; cross-weight in
  `_allosteric_num_den_exprs` and the `_kcat_forward` combine.
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` — ODE-gated multi-`:OnlyA`
  and ping-pong allosteric specs; re-baselined expected counts.
- `test/reference/allosteric_golden_reference.txt` — regenerated.
- `test/test_rate_eq_derivation.jl` — ODE-equality assertion for the two gate mechanisms.
