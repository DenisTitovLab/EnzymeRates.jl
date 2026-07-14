# Allosteric MWC derivation — known correctness issues and the path to a correct fix

**Date:** 2026-07-13
**Status:** KNOWN ISSUES, tracked. The cross-weighting fix attempt was **reverted** (it
introduced regressions worse than the original bug). The `n=1` mass-action ground-truth
harness built during the investigation is **kept** as the acceptance gate for the eventual
correct fix. This document supersedes the fix design in
`2026-07-13-allosteric-mwc-free-enzyme-normalization-design.md` (whose proposed
cross-weighting is now known to be wrong).

## The core bug (pre-existing, unfixed)

The allosteric MWC derivation combines the two conformations as

```
den = Q_A^n + L·Q_I^n,   num = N_A·Q_A^(n-1) + L·N_I·Q_I^(n-1)     (n = catalytic_multiplicity)
```

where `Q_A/Q_I`, `N_A/N_I` are the per-conformation King–Altman denominator/numerator
polynomials. This is **wrong when the inactive graph fragments** — when an `:OnlyA`/`:OnlyI`
binding drops an edge so the inactive free-enzyme segment is reached differently than the
active one. Then `Q_A` and `Q_I` are expressed relative to *structurally different*
free-enzyme references (`D[g_free]`), so combining them on a shared `L`-weighted basis is
inconsistent and leaks a bare catalytic rate constant into the `L`-term (a dimensional
inhomogeneity that breaks the kcat-rescaling contract).

**Reachability:** the bug affects any allosteric mechanism whose inactive graph fragments —
single-`:OnlyA` mechanisms already reach it (pre-existing), and the multi-`:OnlyA`
enumeration move (shipped separately, Tasks 1–3) makes it *more* reachable. The real LDH
i-state 5-/6-group specs are in this class.

## Root-cause understanding (validated)

Each state's King–Altman `Q = Σ_g σ_g·D[g]` is expressed relative to its own free-enzyme
segment spanning-tree weight `D[g_free]`. The physical partition is `P = Q/D[g_free]`
(`[total enzyme]/[free E]`), and the correct MWC combination is on the `P` basis, not the
raw `Q` basis. When the inactive graph fragments, `D_A ≠ D_I` *structurally* (different
monomials — an edge dropped), so `Q_A` and `Q_I` need re-basing. When they differ only by a
conformational `k_A↔k_I` (`:NonequalAI`) substitution with identical topology, they are
*already* on a common basis and must NOT be re-based.

**The discriminator** (validated on all cases below): compare `D_A` and `D_I`'s monomial
support after collapsing every rate constant to one placeholder (`_free_enz_fragments`,
implemented and validated during the investigation). Different support ⇒ topological
fragmentation ⇒ re-base needed. Same support ⇒ rate-constant-only difference ⇒ no re-base.
Dimension alone is NOT the discriminator: fragmenting an *SS* binding (real LDH) keeps the
same `[1/time]` degree yet still needs re-basing; only fragmenting an *RE* binding changes
the degree.

## Why the cross-weighting fix failed (3 confirmed failure modes)

The attempted fix multiplied the active terms by `D_I^n` and inactive by `D_A^n`
(`den = D_I^n·Q_A^n + L·D_A^n·Q_I^n`, numerator likewise). Validated on uni-`:OnlyA`,
multi-`:OnlyA`, and a metabolite-bearing-`D` candidate — but three failure modes surfaced as
the ground truths widened:

1. **`:NonequalAI` catalysis over-correction.** Both conformations productive, identical
   topology, `D_A≠D_I` only via `k_A≠k_I`. Cross-weighting mistakes the rate-constant
   difference for a normalization difference and corrupts the correct `(k_A+L·k_I)` coupling
   that `Q_A+L·Q_I` already produces. Pre-fix 6/6 → post-fix 1/6. *Partially* addressed by
   the `_free_enz_fragments` guard (don't cross-weight same-support `D`).

2. **Metabolite-bearing `D` numerator regression (real LDH).** The numerator clearing is
   mathematically wrong: clearing `P_A^n + L·P_I^n` correctly requires `D_A·D_I^n` on the
   A-numerator term, but the code used `D_I^n`. For a metabolite-bearing `D_I` (real LDH,
   where `D_I ∝ Pyruvate`), this injects a spurious `Pyruvate` factor into the
   product-forming numerator, so the reverse flux dies at `Pyruvate=0`. Confirmed:
   `v(Pyruvate=0)` = `−0.0278`/`−0.0041` (correct, pre-fix) → `−0.0` (wrong, post-fix) for
   the LDH 5-/6-group. Caught by `test_zero_metabolite_finite`. The metabolite-free ground
   truths (`D_I` a bare rate constant) never exercised this.

3. **Ping-pong breaks enumeration/identify.** A fail-loud guard was added for mechanisms
   with multiple free-enzyme forms (ping-pong) because their cross-weighting basis is
   unresolved. But the promote move makes multi-free-form (`:OnlyA`, ping-pong-shaped)
   mechanisms reachable *dynamically* during enumeration/`identify_rate_equation`, so the
   guard throws uncaught and breaks the search. "No static spec triggers it" was true but
   irrelevant — the search reaches them at runtime.

**Conclusion:** cross-weighting is an approximation that keeps needing case-by-case patches
and still fails. The pattern (each widened ground truth reveals a new failure) indicates it
is not the right foundation.

## The correct answer (direction for the eventual fix)

The truly-correct object is the **full coupled two-conformation King–Altman** at `n=1`: build
the explicit two-conformation graph (both conformations, fast conformational flips with
detailed-balance ratios, `:OnlyA`/`:OnlyI` bindings absent from the excluded conformation,
`:NonequalAI` as distinct `k_A/k_I`), and solve its steady state. This is exactly what the
`n=1` mass-action ground-truth harness computes, and it is correct in *every* case tested. It
is feasible because the per-state polynomials and their combination are `n`-independent — `n`
enters only as the `Q^n` exponent (the concerted structure), so validating at `n=1` (and
`n=2` for the power) suffices and the tetramer is never modeled.

Whether the derivation can be refactored to compute the coupled `n=1` King–Altman and then
apply the concerted power, or whether a corrected per-state normalization exists that passes
all the gates below, is the open design question. The `_free_enz_fragments` discriminator and
the corrected numerator clearing (`D_A·D_I^n`, not `D_I^n`) are leads, not a validated
solution — do NOT re-adopt cross-weighting without passing the full gate.

## Acceptance gate for the eventual fix — KEEP AND EXTEND

`test/allosteric_ground_truth.jl` (kept) is an `n=1` two-conformation mass-action steady-state
solver with self-validated ground truths. Any correct fix MUST match it for:

1. uni-`:OnlyA` (RE binding, dead inactive) — dimensional fragmentation.
2. multi-`:OnlyA` bi-uni — dimensional fragmentation, 2 groups.
3. metabolite-bearing `D` (SS binding, dead inactive) — the **LDH regime**, same-dimension
   structural fragmentation.
4. `:NonequalAI` catalysis (both productive, identical topology) — must NOT re-base;
   `Q_A+L·Q_I` is correct.
5. **[add]** a reversible bi-bi where a substrate is absent from one reaction direction and
   `v(that substrate = 0) ≠ 0` — the LDH `test_zero_metabolite_finite` case (failure mode 2).
6. **[add]** an allosteric ping-pong (multiple free-enzyme forms) — failure mode 3; the fix
   must produce a correct equation, or the pipeline must handle these mechanisms without a
   hard error.

Each ground truth self-validates first (`L=0` → non-allosteric active rate; all-`:EqualAI` →
base rate, `L`-independent) before it may gate the fix. Dimensional homogiety and
kcat-rescaling are necessary but NOT sufficient — they missed failure modes 2 and 3.

## What was reverted / kept

- **Reverted:** the cross-weighting fix + `_free_enz_fragments` guard + `D[g_free]` surfacing
  in `src/rate_eq_derivation.jl` — the derivation returns to the pre-fix state (original
  `:OnlyA` L-term leak present, but LDH correct and enumeration/identify working).
- **Kept:** the multi-`:OnlyA` enumeration move (Tasks 1–3 — correct, independent); the `n=1`
  mass-action ground-truth harness + the `:NonequalAI`-catalysis (candidate-1) guard (the
  reusable acceptance gate).
- **Known-issue markers:** the ground-truth gates that fail against the pre-fix derivation
  (uni-`:OnlyA`, multi-`:OnlyA`, metabolite-bearing-`D`) are `@test_broken` — they document
  the bug and will flip to passing when the correct fix lands.

## Meta-lesson

The `n=1` mass-action ground-truth gate is what exposed all three failure modes; hand-analysis
(across two sessions and Denis) repeatedly reached plausible-but-wrong conclusions. Any future
work on this derivation must be gated by the ground truth, not by dimensional/kcat arguments
alone.
