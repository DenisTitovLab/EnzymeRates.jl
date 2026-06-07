# Rate Equations Without Division by Concentration — Design

**Date:** 2026-06-06
**Status:** design proposed, pending spec review → implementation plan

---

## Goal

Make every rate equation the package derives a clean ratio of polynomials in the
metabolite concentrations — **no concentration ever in a denominator**. Today many
derived equations contain `1/[conc]` terms (e.g. `… / (NADH * Pyruvate)`). When a data
point sets that metabolite to `0` (which real datasets do constantly, since each
measurement uses only a subset of metabolites), `rate_equation` evaluates to `NaN` /
`Inf` (and, for some topologies, a spurious `0.0`), so the mechanism cannot be fit and
is silently dropped from identification.

The rate is mathematically finite at those points — the `1/[conc]` form is a
representation artifact, not a property of the rate. We fix the **derivation** so it
emits the equation in a form that is finite wherever the rate is defined.

## Evidence (LDH, `substrates: NADH, Pyruvate; products: Lactate, NAD; oligomeric_state: 4`)

Measured on the **deduped** `init_mechanisms` set — i.e. exactly what
`identify_rate_equation` fits (it runs `_dedup_flat!`):

- **69 / 69** mechanisms contain a division by some concentration
  (by `NADH` 17, `Pyruvate` 38, `Lactate`&`NAD` 4, `NADH`&`Pyruvate` 10).
- Current `rate_equation` (random positive params), count returning non-finite:
  - single `Pyruvate = 0`: 48/69 · single `NADH = 0`: 27/69 · single `Lactate`/`NAD = 0`: 4/69 each
  - **reverse** initial velocity (`NADH = Pyruvate = 0`): **65/69**
  - **forward** initial velocity (`Lactate = NAD = 0`): 4/69
- Baseline (all concentrations > 0): all 69 finite — the equations are correct for
  nonzero concentrations; only the zero-concentration evaluation breaks.

(Before `_dedup_flat!`, only 14/69 show the problem — which is why an un-canonicalized
spot check undercounts. The pipeline canonicalizes step order first, which moves the
derivation's reference form onto a bound complex and exposes the issue in every
equation.)

## Root cause

`rate_equation` and `rate_equation_string` both flow through
`_raw_symbolic_rate_polys` (`src/rate_eq_derivation.jl`), which derives the
numerator/denominator polynomials by the Cha/King-Altman method. Within a
rapid-equilibrium group, each enzyme form's relative population is expressed
**relative to a chosen reference form** — the first form in canonical step order. After
`_dedup_flat!` sorts the steps, that reference is typically a fully-loaded complex
(e.g. `E·NADH·Pyruvate`). Expressing the free enzyme relative to a loaded complex
requires dividing out the bound substrates:

```
[E] / [E·NADH·Pyruvate]  =  (K_NADH · K_Pyruvate) / (NADH · Pyruvate)
```

The `normalize` step (`src/rate_eq_derivation.jl:373-389`) then bakes this into the
final equation: for the single-RE-group case (`G == 1`, which every unicyclic init
mechanism is) it divides the numerator by `sigma_den[1]` and each denominator term by
`alpha_den[i]`. Those divisors carry concentrations whenever the reference is a bound
complex, producing `1/[conc]` in both numerator and denominator. At `[conc] = 0` both
blow up and `∞/∞ = NaN` — a **removable** singularity (the limit is finite).

Choosing a different reference cannot fix this in general: multi-"valley" RE groups
(free `E` and a covalent `Estar` reachable from each other only *through* a binding
step) force an unbinding hop no matter where the traversal starts. (Verified
experimentally — reseeding the BFS at the free enzyme left 14/69 still divided.)

## The fix (derivation-level)

A fraction's value is unchanged by multiplying numerator and denominator by the same
thing: `N/D = (N·X)/(D·X)`. We use that freedom to emit the equation with no
concentration in a denominator, as the **reduced (lowest-terms) polynomial ratio**.

Two changes to `_raw_symbolic_rate_polys`:

1. **Stop dividing by the bound-reference population.** Remove the `normalize` branch.
   Build the denominator from `sigma_num[g]` for all groups (the path already used when
   `G > 1`), and leave the numerator undivided. So built, `N` and `D` are plain
   polynomials with **non-negative** concentration powers — but they may carry a common
   concentration factor (from the reference choice), which would give `0/0`.

2. **Reduce the pair to lowest terms over concentrations.** Divide both `N` and `D` by
   their concentration **monomial GCD**: for each metabolite symbol, the smallest power
   it has across all terms of `N ∪ D`; divide it out of every term. Because we only ever
   divide by the *minimum* power present, no term can go negative — no concentration
   denominator can appear. The least-loaded form's term becomes concentration-free,
   giving the constant denominator term that keeps the rate finite at `[conc] = 0`.

Reduction is over **concentrations only**, never parameters. (Reducing over parameters
could divide out a parameter common to every term, silently dropping a fitted parameter
from the equation. Concentration-only reduction cannot do this — verified: 0 parameters
dropped across all 69 LDH mechanisms.)

### Validation (deduped LDH set, 69 mechanisms)

With "skip `normalize` + concentration-only lowest-terms reduction":

- **0** correctness mismatches vs. the current derivation for random concentrations > 0
  (69 mechanisms × 3 random draws) — identical rate values.
- **0** equations with division by concentration.
- **0** fitted parameters dropped.
- **0** non-finite results for every realistic zero-pattern: each metabolite zeroed
  one-by-one, forward (`products = 0`), reverse (`substrates = 0`).

### Code changes

- `src/sym_poly_for_rate_eq_derivation.jl`: add helper
  `_reduce_conc_lowest_terms(num::POLY, den::POLY, conc::Set{Symbol}) -> (POLY, POLY)`
  (next to the existing POLY helpers). For each symbol in `conc`, compute the minimum
  exponent over all monomials of both polys; build the GCD monomial from the nonzero
  minimums; divide both polys by it via the existing `_poly_div_mono`. Returns the inputs
  unchanged when the GCD is empty.
- `src/rate_eq_derivation.jl` `_raw_symbolic_rate_polys` (the 5-arg method, lines
  333-399): delete the `normalize` flag and the `_poly_div_mono(num, sigma_den[1])`
  division; the denominator loop always uses `sigma_num[g]`; as the final step (after the
  `rename_map` pass) call `_reduce_conc_lowest_terms(num, den, conc_set)`, where
  `conc_set` is the mechanism's metabolite symbols.
- `src/rate_eq_derivation.jl` `_compute_alpha` (lines 250-310): `sigma_den` is now dead
  (its only uses were the deleted `normalize` lines). Stop computing and returning it;
  return `(alpha_num, alpha_den, sigma_num)`.

### kcat is unaffected (no extra work, but must stay green)

`_kcat_forward` (`src/rate_eq_derivation.jl:699`) reads the same `_raw_symbolic_rate_polys`
output, groups monomials by **concentration pattern**, and computes kcat as the **ratio**
`num_k / den_k` at matching patterns. Skipping `normalize` and reducing to lowest terms
both multiply `N` and `D` by a common monomial; the concentration patterns of `N` and `D`
shift by the same amount (matching preserved, patterns stay distinct) and the parameter
parts cancel in the ratio. kcat is therefore invariant. The existing kcat / rescaling /
scale-invariance tests are the gate.

### Allosteric path

`_allosteric_num_den_exprs` (`src/rate_eq_derivation.jl:1518`) gets its per-state
catalytic polys from `_raw_symbolic_rate_polys_allosteric` → `_raw_symbolic_rate_polys`,
so each state is already reduced (no concentration denominators). It then assembles the
MWC numerator/denominator by **multiplication and addition only** (`make_num_term` /
`make_den_term`, `_nest_binary`, `_power_expr`), which cannot introduce a concentration
denominator. The allosteric regression tests (below) are the gate; if a residual
*common* concentration factor between the combined numerator and denominator ever leaves
a `0/0` at a realistic pattern, a follow-up can apply the same reduction to the combined
form. Not expected to be needed.

## Residual: all-metabolites-zero (out of scope)

14 / 69 mechanisms still give `0/0` when **every** metabolite is `0` simultaneously.
These are multi-valley topologies where free `E` and `Estar` interconvert only through
metabolite binding, so no enzyme form has a concentration-independent population and the
reduced denominator has no constant term. At all-zero there is genuinely no steady-state
enzyme distribution — `0/0` is the honest answer, and it is never a real data row (no
substrate and no product present means nothing was measured). We leave it as-is and
document it. If such a row ever reached fitting it would surface as a loud non-finite
loss (the pipeline already raises on non-finite folds), not silent corruption. The
regression tests below deliberately do **not** assert finiteness at all-metabolites-zero.

## Testing

Per the agreed scope: every shared fixture exercised by fitting, plus enumerated
mechanisms, each metabolite zeroed one-by-one, asserting a finite rate.

1. **Shared fixtures** (`test/test_rate_eq_derivation.jl`, inside the
   `run_all_tests(spec)` loop over `MECHANISM_TEST_SPECS`). For each fixture: draw random
   positive fitted params (`random_reduced_params`) and random positive concentrations;
   for each metabolite, set just that one to `0.0` (others positive) and assert
   `rate_equation(m, concs, params)` is **finite and nonzero**. The nonzero check guards
   the spurious-`0.0`-from-`1/Inf` regression Denis flagged; with all other metabolites
   positive a genuine zero is measure-zero. This covers the non-allosteric and allosteric
   fixtures (allosteric path coverage).

2. **Enumerated mechanisms** (new test, POLY-level so it is cheap and needs no
   `@generated` compile). For the LDH reaction above, run `init_mechanisms` +
   `_dedup_flat!`; for each mechanism assert `_raw_symbolic_rate_polys` yields polynomials
   with **no negative concentration exponent**, and that numeric evaluation of `num/den`
   (random positive params) is finite for each metabolite zeroed one-by-one, for
   `products = 0`, and for `substrates = 0`. Include at least one additional reaction
   shape (e.g. a uni-uni or uni-bi) so the invariant is checked beyond bi-bi. Optionally
   include one round of `expand_mechanisms` on the LDH base set.

3. **Snapshot updates.** The 20 `expected_factored_num` / `expected_factored_denom`
   snapshots in `mechanism_definitions_for_test_enzyme_derivation.jl` and the
   Expr-shape / flat-string regression tests in `test_rate_eq_derivation.jl` will change
   to the cleared forms; regenerate and eyeball them (the leading numerator term is no
   longer a bare rate constant — it now carries the substrate monomial).

4. **Must stay green unchanged.** `test_rate_equation_performance` (allocation-free,
   sub-100 ns) — the cleared equations have similar or fewer operations (fewer concentration
   divisions); kcat / rescaling / scale-invariance tests; Aqua/JET.

## Non-goals

- **Not** abandoning the Cha rapid-equilibrium lumping for full steady-state King-Altman.
  That would explode parameter counts and compile cost (the whole point of the RE
  treatment is small, few-parameter equations). The `1/[conc]` artifact is inherent to
  the Cha *intermediate* algebra for these topologies; we cancel it in the final reduced
  form rather than avoid forming it.
- **Not** guaranteeing finiteness at the chemically-degenerate all-metabolites-zero point.
- **Not** changing the parameter API (`parameters` / `fitted_params`), the
  thermodynamic-constraint machinery, or the fitting loop.

## TDD order (for the plan)

1. Add the enumeration POLY-level test (fails today: many negative concentration
   exponents) and the fixture zero-concentration test.
2. Add `_reduce_conc_lowest_terms` + rewire `_raw_symbolic_rate_polys`; drop `sigma_den`
   from `_compute_alpha`.
3. Make the new tests pass; regenerate snapshots; confirm kcat/perf/Aqua/JET green.
4. Confirm on the LDH `identify_rate_equation` run that reverse-direction mechanisms now
   fit instead of dropping out.
