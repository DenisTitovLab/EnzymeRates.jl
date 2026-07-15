# Allosteric MWC free-enzyme normalization — formulation-1 design

**Date:** 2026-07-14
**Status:** Design, validated. Supersedes the reverted cross-weighting attempt and the
`2026-07-13` known-issues framing. The fix below matches every `n=1` ground truth in
`test/allosteric_ground_truth.jl` with a single uniform rule and no discriminator.

## The bug

The MWC derivation combines the two conformations as `den = Q_A^n + L·Q_I^n`,
`num = N_A·Q_A^(n-1) + L·N_I·Q_I^(n-1)`, where `Q_state` and `N_state` are the
King–Altman denominator and numerator polynomials of each conformation's catalytic
graph. Each `Q_state` is the raw Matrix–Tree polynomial `Σ_g σ_g·D[g]`, expressed
relative to the free-enzyme segment's spanning-tree weight `D_state[g_free]`. That weight
is `1` only when the graph has a single rapid-equilibrium segment. When an `:OnlyA` /
`:OnlyI` binding drops an edge, the inactive graph fragments into more segments, so
`D_I[g_free] ≠ 1` and `Q_I` carries that weight as a factor. Combining the two raw
polynomials therefore adds partition functions written on different bases, and the
mismatch leaks a bare rate constant (or a metabolite factor) into the `L`-term.

Reproduced on the uni-uni `:OnlyA` mechanism: the current L-term is
`k_ES_to_EP + k_EP_to_ES·P/K_P_E + k_ES_to_EP·P/K_P_E`, whose bare `k_ES_to_EP` sits
beside the dimensionless `1` of `Q_A`. Against the `n=1` mass-action ground truth the
current derivation is off by up to 30%; dividing `Q_I` by `D_I[g_free] = k_ES_to_EP`
matches it to 3e-7.

## The model: formulation 1

The MWC model this package derives is **formulation 1**: only the free enzyme flips
between the active and inactive conformations, with free-enzyme ratio `L`. Equivalently,
every form may flip with a Wegscheider-fixed ratio `L'` (which equals `L` only when its
bindings are `:EqualAI`), but the enzyme commits to one conformation for the length of a
catalytic cycle. Physically, each enzyme molecule catalyzes entirely in one conformation.

Formulation 2 — every catalytic intermediate flips mid-cycle — routes turnover through
the faster conformation and gives a slightly higher rate. The two models agree for every
dead-inactive case and for `:EqualAI` catalysis; they diverge only for live
`:NonequalAI` catalysis, by 0.1–3%. We derive formulation 1.

## The fix: uniform per-state normalization

Normalize each conformation to its own free-enzyme weight before combining:

```
den = (Q_A/D_A)^n + L·(Q_I/D_I)^n
num = (N_A/D_A)·(Q_A/D_A)^(n-1) + L·(N_I/D_I)·(Q_I/D_I)^(n-1)
```

where `D_A ≡ D_A[g_free]`, `D_I ≡ D_I[g_free]`. The value is the same however we render
it; the rendering rule below keeps the common cases in clean standard-MWC form and falls
back to a polynomial cross-weight only when division would put a concentration in a
denominator.

**Rendering rule**, applied per conformation:

1. **`D_A == D_I` exactly → raw combine.** Identical conformations (all-`:EqualAI`, or
   `:EqualAI` catalysis) produce the same `D` polynomial in both states, so the factor
   cancels and `Q_A^n + L·Q_I^n` already equals formulation 1. Skipping the normalization
   here avoids a redundant `D^n` factor that would bloat multi-segment equations like
   ping-pong. This is a simplification, not a correctness switch — the three rates are
   identical when `D_A == D_I`. Test it by exact polynomial equality, not the reverted
   `_free_enz_fragments` heuristic.

2. **`D` a metabolite-free monomial → divide.** Replace `Q → Q/D`, `N → N/D`. Dividing a
   `POLY` by a single rate-constant term just subtracts exponents, so the result stays a
   Laurent `POLY` and renders in standard-MWC form with a leading `1`. Uni- and
   multi-`:OnlyA` land here: `Q_I/D_I = 1 + P/K_P_E + (k_EP_to_ES/k_ES_to_EP)·P/K_P_E`,
   and `den = Q_A + L·(Q_I/D_I)`. The only denominators introduced are ratios of rate
   constants, computed once per call — allocation-free.

3. **`D` carries a metabolite, or is a multi-term sum → cross-weight.** Division would
   either place a concentration in a denominator (forbidden) or produce a rational the
   `POLY` type cannot hold. Clear both sides by `(D_A·D_I)^n` instead — exact, and it
   keeps the equation a polynomial:

   ```
   den = D_I^n·Q_A^n           + L·D_A^n·Q_I^n
   num = D_I^n·N_A·Q_A^(n-1)    + L·D_A^n·N_I·Q_I^(n-1)     (drop the L-term when N_I = 0)
   ```

   The LDH metabolite-D and `:NonequalAI` cases land here. The A-term numerator factor is
   `D_I^n`, **not** the known-issues doc's `D_A·D_I^n` — that came from clearing numerator
   and denominator by different factors; the same `(D_A·D_I)^n` on both gives `D_I^n·N_A`.

The rule carries no fuzzy discriminator: rules 2 and 3 pick a *rendering* of the same
value, distinguished only by whether `D` divides cleanly.

## Validation completed this session

All checks use the `n=1` two-conformation mass-action harness in
`test/allosteric_ground_truth.jl`; scripts live in the session scratchpad.

- **Uniform rule matches all four ground truths** to ~1e-7 — uni-`:OnlyA`,
  multi-`:OnlyA`, LDH metabolite-D, and a formulation-1 (free-flip) `:NonequalAI`
  reference. One rule, no discriminator.
- **Formulations 1 and 2 are distinct**, and only for `:NonequalAI` catalysis. In
  formulation 2 the catalytic-intermediate flip rungs carry nonzero net enzyme current
  (0.035 on `E(A,B)`/`E(P)`) when `k_A ≠ k_I`, and zero current when `k_A = k_I`. The
  normalized combine equals a free-flip-only network to 5e-8.
- **A metabolite inside `D[g_free]` behaves correctly at zero concentration.** When the
  inactive is dead and `D_I ∝ B`, the enzyme traps and `v → 0`; the cleared polynomial
  returns exactly `0` where the rational form returns `NaN`. When the inactive is
  productive and `D` stays finite, reverse flux survives (`v ≠ 0`). No spurious zero.
- **Ping-pong has one free form.** `E(; residual = A - P)` carries a residual, so
  `_segment_root` picks `E` uniquely; `D[g_free]` is unambiguous. The reverted "multiple
  free forms" guard read an outdated macro form and guarded a non-problem. The allosteric
  ping-pong derives without a crash, and its all-`:EqualAI` form equals the
  non-allosteric ping-pong to 1e-17, `L`-independent.

## Work to do

**Derivation** (`src/rate_eq_derivation.jl`)
- Surface `D_state[g_free]` per conformation. The reverted commit `0e4c556` built this;
  resurrect that part.
- In `_allosteric_num_den_exprs`, apply the rendering rule: raw when `D_A == D_I`, divide
  `Q/D` when `D` is a metabolite-free monomial, cross-weight otherwise. Apply the same rule
  in `_kcat_forward`, which the kcat-rescaling test compares against.
- Drop the reverted `_free_enz_fragments` guard and the ping-pong fail-loud guard; neither
  is needed.

**Gates** (`test/allosteric_ground_truth.jl`)
- Flip `@test_broken` → `@test` for uni-`:OnlyA`, multi-`:OnlyA`, and LDH metabolite-D.
- Point the `:NonequalAI` gate at the free-flip-only reference network; its expected value
  shifts by 0.1–3%.
- Add a metabolite-at-zero gate covering both the trap (`v → 0`) and the surviving reverse
  flux (`v ≠ 0`).
- Add an allosteric ping-pong gate. Build the eight-form free-flip network and self-validate
  it (`L = 0` → active-only rate; `k_I = k_A` → base rate, `L`-independent) before it gates
  the derivation.

**Performance** (`test/test_rate_eq_derivation.jl`)
- Confirm `rate_equation` stays allocation-free and under 120 ns. The reverted cross-weight
  held this bound; re-verify, since the uniform rule cross-weights more mechanisms.

**Goldens**
- Regenerate `allosteric_golden_reference.txt` and any other allosteric golden strings. Only
  mechanisms whose inactive fragments (or whose catalysis is `:NonequalAI`) change.

**Docs** (`docs/src/deriving/mwc_allostery.md`)
- Rewrite the model description to state formulation 1 and contrast the mid-cycle-crossing
  alternative. Replace the raw `num`/`den` block (currently the buggy combine) with the
  cross-weighted form. The `@example` blocks execute the code, so this must land in the same
  change as the derivation fix. Write it with the elements-of-style skill.

## Deferred

- **`n > 1`.** Every ground truth is `n=1`. The concerted power `P^n` follows from the
  commit-when-free model: each subunit runs a committed-conformation cycle with partition
  `P_state`, and the oligomer weights `P_A^n` against `L·P_I^n`. Not directly gated.
- **The reverted numerator/ping-pong machinery.** Reuse the `D[g_free]` surfacing; rewrite
  the combine logic (uniform, `D_A == D_I` skip, corrected numerator, no discriminator).

## Acceptance gate

A correct fix matches `test/allosteric_ground_truth.jl` for: uni-`:OnlyA`, multi-`:OnlyA`,
LDH metabolite-D, formulation-1 `:NonequalAI`, metabolite-at-zero (trap and survive), and
allosteric ping-pong. Each ground truth self-validates (`L = 0` → active-only rate;
identical conformations → base rate, `L`-independent) before it gates the derivation.
`rate_equation` stays allocation-free and under 120 ns.
