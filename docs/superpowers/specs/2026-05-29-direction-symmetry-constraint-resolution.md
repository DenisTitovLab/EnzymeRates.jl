# Direction-Symmetry Principle for Thermodynamic Constraint Resolution — Design Knowledge

**Date:** 2026-05-29
**Status:** Design knowledge for a future follow-up PR. Not implemented here.
This document exists so the reasoning is not lost.
**Companion:** `2026-05-29-nonequalai-rank-validity.md` covers the related
rank/nullspace validity algorithm.

## Why this document exists

While investigating the EqualAI × NonequalAI coupling bug we worked through
the thermodynamics of MWC allosteric mechanisms carefully and uncovered a
deeper, more principled way to resolve thermodynamic (Haldane / Wegscheider)
constraints — one that applies to **every reversible mechanism, not just
allosteric ones**.

**The core idea (general):** rewrite dependent-parameter removal so the
constraint is resolved in a **direction-symmetry-invariant** manner — it must
not matter whether the user wrote the reaction forward or reverse. The
constraint pins a *ratio* (`k_for/k_rev`, or an equilibrium constant); keep
the direction-*symmetric* combination (the *speed* `√(k_for·k_rev)`) as the
free parameter and derive both rate constants from speed + the pinned ratio,
distributing any forced difference by the unique minimum-norm symmetric
solution. The current code instead drops one rate constant (e.g. `k_rev`),
which silently privileges the user's arbitrary forward/reverse labeling.

This is a major change (it touches the constraint-resolution core *and*
every hand-written analytical formula, which gain `√` terms) but also a
major **unification and correctness fix**. It should be implemented as its
own large feature.

The allosteric (A/I) case below is the most *visible* manifestation (it
changes the model, not just the parametrization), but the principle is
general — see "The principle is general" below.

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

### Degeneracy / rank test (companion document)

A NonequalAI group buys a *genuine* free parameter only when its split can be
nonzero given the free absorbers (`F = {NonequalAI RE groups} ∪ {SS groups in
the cycle}`). The clearest degenerate case is a **lone NonequalAI binding K in
a pure-RE Wegscheider loop with no SS step to absorb the split**; full-rank
multi-cycle mechanisms can also force a split to zero even with SS absorbers
present. This **validity/degeneracy algorithm is captured in its own document**,
`2026-05-29-nonequalai-rank-validity.md` — a distinct, durable concern (config
validity / enumeration) that remains necessary even after this
symmetric-resolution rewrite, since a pure-RE loop has no speed DOF to absorb a
split regardless of how constraints are resolved.

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

## 5a. The principle is general — non-allosteric mechanisms too

The forward/reverse asymmetry is **not** created by the allosteric tags. It
is present in *every* reversible mechanism the moment a Haldane/Wegscheider
constraint eliminates a dependent parameter. "Drop `k_rev`, keep `k_for`
free" privileges the forward direction even in a single-state,
non-allosteric mechanism: write the same enzyme backwards and the code keeps
the *other* rate constant free, yielding a different parametrization and a
different-looking analytical formula for the identical enzyme.

The symmetric resolution applies uniformly — keep the speed
`√(k_for·k_rev)` free, derive both rate constants from speed + the pinned
ratio. Two regimes of the *same* principle:

- **Non-allosteric (single state):** the symmetric choice does **not** change
  the model (same achievable rate curves — it is a reparametrization), but it
  makes the canonical parametrization and the rate-equation **form**
  direction-invariant, and parameter reporting reproducible. Cost: `√` enters
  even simple rate equations, and every analytical formula in the suite must
  be rewritten.
- **Allosteric (A/I):** the symmetric distribution of the forced state-split
  additionally changes the **model** itself (different I-state rate
  equation), and is the genuinely-more-correct, direction-invariant choice
  (`k_for_I · k_rev_I = k_for_A · k_rev_A`).

So the follow-up is not "an allosteric tweak" — it is a **rewrite of
dependent-parameter removal** for all mechanisms, with the allosteric case as
the most visible payoff. `_dependent_param_exprs_kernel` is the single
chokepoint for both.

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
   of the derivation and must be coordinated with the current constraint
   machinery.
3. **kcat / Vmax normalization interacts.** kcat is defined from the SS
   forward rates, which option (c) makes state-dependent (`k_f^A ≠ k_f^I`).
   Re-check the normalization invariants (kcat homogeneity, scale
   invariance) under speed-sharing.
4. **Parameter naming chokepoint.** The free SS symbol becomes a *speed*
   with derived per-state `k_f/k_r`; new parameter symbols must flow through
   the `Parameter` family and `name(p, m)` naming chokepoint.

## 7. Implementation Target

Implement this principle as its own large feature: re-found constraint
resolution on the symmetric (speed, ratio) basis, accept `√` in rate
equations, update analytical formulas, re-verify kcat normalization, and fold
the rank-based degeneracy rejection into the same foundation.
