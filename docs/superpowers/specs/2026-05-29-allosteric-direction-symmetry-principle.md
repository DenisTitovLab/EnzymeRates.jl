# Allosteric Direction-Symmetry Principle — Design Knowledge

**Date:** 2026-05-29
**Status:** Design knowledge for a FUTURE follow-up PR (after the parent
structural-parameter-names refactor lands). Not implemented here. This
document exists so the reasoning is not lost.
**Companion:** `2026-05-29-equalai-nonequalai-coupling-design.md` (the
*contained* fix that ships first).

## Why this document exists

While investigating the EqualAI × NonequalAI coupling bug we worked through
the thermodynamics of MWC allosteric mechanisms carefully and uncovered a
deeper, more principled way to resolve thermodynamic constraints in
two-state (active/inactive) enzymes. It is too large and too entangled with
the parent refactor's derivation core to implement now, but it is the
correct long-term foundation. Denis's decision: ship the contained fix now,
finish the parent refactor, then do this as its own (large) PR.

## 1. The physical model (corrected)

Both MWC states (A = active, I = inactive) share the **same enzyme-form
graph and the same steps**. They differ only in parameter *values*, with
the per-group allosteric tag (`:OnlyA`, `:EqualAI`, `:NonequalAI`)
controlling which symbols are shared, independent, or state-restricted.

Each thermodynamic cycle (Haldane, net `b≠0`; or Wegscheider, `b=0`)
imposes, in log space, the same relation in both states:

```
Σ_g c_g · log(symbol_state,g) = b · log(Keq)
```

where `c_g` is the cycle's traversal incidence of group `g`'s step.
Subtracting the A and I equations for any cycle live in both states gives
the **split constraint**:

```
Σ_g c_g · d_g = 0,   d_g = log(symbol_A,g) − log(symbol_I,g)
```

## 2. Key realization — dependent parameters absorb splits

Thermodynamic constraints **pin ratio-type quantities** (equilibrium
constants, `k_f/k_r`). Each independent cycle eliminates exactly one
*dependent* parameter (a pivot), computed per state from the others.

A **dependent parameter's split is free** — it is derived, not a shared
user value. In particular, an **SS step's reverse rate** is (almost
always) a free absorber: it has no "shared affinity" meaning, so it can
differ between states without contradicting any tag.

This is why **a single NonequalAI binding group in a Haldane cycle is
perfectly valid** (hand-verified mechanism PK: PEP binding `:NonequalAI`,
catalysis `:EqualAI`; the catalytic reverse `k5r`/`k5r_T` absorbs the
PEP-affinity split). An earlier "≥2 NonequalAI per cycle" rule was wrong
because it omitted dependent (SS-reverse) absorbers from the accounting.

### Degeneracy / rank test (still useful)

Free-absorbing column set `F = {NonequalAI RE groups} ∪ {SS groups in the
cycle}`. A NonequalAI group `g` buys a *genuine* free parameter iff
`rank(C_live[:,F]) == rank(C_live[:,F∖{g}])`. The **only** genuinely
degenerate case is a **lone NonequalAI binding K in a pure-RE Wegscheider
loop with no SS step to absorb the split** — there is no free absorber, so
the split is forced to zero and the NonequalAI tag is vacuous. Such configs
should be rejected (or normalized) as over-parametrized.

## 3. The direction-invariance defect

For an SS catalytic step, write `f = log k_f`, `r = log k_r`. The
thermodynamic constraint pins the **ratio** `v = f − r`; the **speed**
`u = f + r` (= `2·log √(k_f·k_r)`) is a free, constraint-free kinetic DOF.

When a NonequalAI binding K forces `v` to differ between states, the
difference must be absorbed somehow. The current behavior ("drop `k_rev_T`")
keeps `k_f` shared and dumps the entire difference on `k_r`. That is the
classic **MWC K-system** hypothesis (forward kcat state-independent) — but
it **privileges the user's arbitrary forward/reverse labeling**:

> Defining the reaction `S → P` and sharing `k_f` is *not the same model* as
> defining the **same enzyme** `P → S` and sharing `k_f` (now the other
> direction). The labels swap `f ↔ r`, turning "share forward" into "share
> reverse." **Same enzyme, two different rate equations, decided only by how
> the user typed the reaction.** A real reproducibility bug.

## 4. The symmetric resolution (option c)

Share the **speed** `u = f + r` (direction-symmetric) between states, and
let the constraint-forced **ratio** `v` differ. Because `u` is symmetric
under `f ↔ r`, the model is **invariant to the forward/reverse labeling**.

### Concrete: one NonequalAI binding K

Mechanism `E+S⇌ES` (`K_S`, NonequalAI), `ES⇌EP` (`k_f,k_r`, EqualAI),
`EP⇌E+P` (`K_P`, EqualAI). Haldane: `k_f/k_r = Keq·K_S/K_P`. Let
`ρ = K_S^I/K_S^A`, and share speed `s = √(k_f·k_r)`. All four rate
constants become derived:

```
A:  k_f^A = s·√(Keq·K_S^A/K_P),   k_r^A = s/√(Keq·K_S^A/K_P)
I:  k_f^I = k_f^A·√ρ,             k_r^I = k_r^A/√ρ
```

The state-difference splits **symmetrically**: forward up by `√ρ`, reverse
down by `√ρ`, product `k_f·k_r` shared. (Contrast option (a):
`k_f^I = k_f^A`, `k_r^I = k_r^A/ρ`.) Same parameter count: one speed `s`
plus the binding K's.

The `√` is **intrinsic** to symmetry — half the log-ratio-difference goes
each way. Denis has accepted `√` in the rate equation (positive parameter
combinations; smooth and fast).

## 5. The general principle

Split every reversible step's two log-rate-constants into:
- **ratio** `v = log k_f − log k_r` — *antisymmetric* under forward↔reverse;
  this is what thermodynamic constraints pin.
- **speed** `u = log k_f + log k_r` — *symmetric*; a free kinetic DOF no
  constraint touches.

**Principle:** constraints act only on the antisymmetric (ratio)
coordinates; resolve each constraint's forced state-difference by the
**unique minimum-norm (symmetric) distribution** across the participating
steps' ratios, and keep the symmetric speeds shared (EqualAI).

Consequences:
- **Removes all pivot arbitrariness.** The current single-pivot elimination
  is a basis/label choice; the min-norm symmetric distribution is
  basis- and direction-independent — one canonical model per tag config.
- **Subsumes the "absorber choice" question.** When multiple SS steps share
  a cycle, instead of enumerating "which reverse absorbs," the symmetric
  distribution spreads the difference equally — one canonical model. The
  asymmetric alternatives remain reachable by explicitly tagging a step
  `:NonequalAI` (spending a parameter to assert that step carries the
  state-dependence).
- **Consistent with the pure-RE rejection.** RE steps at equilibrium have
  only a ratio (`Kd`), no speed DOF, so a lone-NonequalAI pure-RE
  Wegscheider loop has nothing symmetric to absorb the split → genuinely
  degenerate → rejected.
- **Unifies the taxonomy:** speeds are shared kinetic DOF; ratios are
  thermodynamically pinned and symmetrically split; `:NonequalAI` buys
  explicit asymmetry; absorber-less loops are rejected.

## 6. Implementation notes & costs (for the follow-up PR)

1. **`√` / rational powers in the rate equation.** Derived rate constants
   gain `√(Keq·K/K)` terms; the rate equation is no longer purely rational.
   - Verify the `rate_equation` perf invariant survives `sqrt`
     (`allocs == 0`, `t < 100 ns` in `test/test_rate_eq_derivation.jl`).
     `sqrt` is a few ns; expected fine.
   - Confirm `build_power_expr` emits rational (`1//2`) exponents (the
     elimination already uses `Rational{BigInt}`, so likely yes).
   - Every hand-written analytical formula in `MECHANISM_TEST_SPECS`
     (including `pk_rate_analytical`) must be updated to the `√` form.
2. **Re-founds `_dependent_param_exprs_kernel`.** Single-pivot Gaussian
   elimination → symmetric (pseudoinverse / projection) resolution that
   splits ratio-differences in the (speed, ratio) basis. This is the *core*
   of the derivation and the **parent refactor's active code** — must be
   coordinated, not smuggled in.
3. **kcat / Vmax normalization interacts.** kcat is defined from the SS
   forward rates, which option (c) makes state-dependent (`k_f^A ≠ k_f^I`).
   Re-check the normalization invariants (kcat homogeneity, scale
   invariance) under speed-sharing.
4. **Parameter naming chokepoint.** The free SS symbol becomes a *speed*
   with derived per-state `k_f/k_r`; new symbols/derivations flow through
   `name(...)` / `_param_symbol`. This is the parent refactor's chokepoint
   territory.

## 7. Sequencing (decided)

1. **Now:** contained fix (companion doc) — repair the synth-dep overwrite
   so PK / `m_mixed` compute consistently under the current option-(a)
   convention; unblock the parent session. No representation change.
2. **Parent refactor finishes** (structural parameter names).
3. **Then:** this principle as its own (large) PR — re-found the constraint
   resolution on the symmetric (speed, ratio) basis, accept `√` in rate
   equations, update analytical formulas, re-verify kcat normalization, and
   fold the rank-based degeneracy rejection into the new foundation.
