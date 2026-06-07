# Representation-Independent Rate-Equation Derivation — Design

**Date:** 2026-06-06 (revised 2026-06-07 after review)
**Status:** design proposed, pending spec review → implementation plan

---

## Goal

A mechanism's derived rate equation should depend only on *what the mechanism is*,
never on *how its steps happen to be written or ordered*. Two coupled problems break
that today:

1. **Division by concentration.** Many derived equations contain `1/[conc]` terms (e.g.
   `… / (NADH * Pyruvate)`). When a data point sets that metabolite to `0` — which real
   datasets do constantly, since each measurement uses only a subset of metabolites —
   `rate_equation` returns `NaN`/`Inf` (and, for some topologies, a spurious `0.0`), so
   the mechanism cannot be fit and is silently dropped from identification.
2. **Mutating canonicalization.** `_dedup_flat!` rewrites mechanisms in place to compare
   them; that mutation is dangerous and is also what turns problem 1 from a 14/69 issue
   into a 69/69 issue (it relocates the derivation's reference form onto a bound complex).

The fix has two complementary parts: **(A)** emit each equation as a reduced polynomial
ratio with no concentration in any denominator, making the *derivation* independent of
step order; and **(B)** make mechanisms canonical *by construction*, making the *struct*
unique so dedup compares without mutating and nothing downstream depends on step order.

## Evidence (LDH, `substrates: NADH, Pyruvate; products: Lactate, NAD; oligomeric_state: 4`)

- **As-written `init_mechanisms` (no dedup mutation):** 14/69 equations divide by a
  concentration.
- **After `_dedup_flat!` (what `identify_rate_equation` fits):** **69/69** divide by a
  concentration (`NADH` 17, `Pyruvate` 38, `Lactate`&`NAD` 4, `NADH`&`Pyruvate` 10). The
  dedup canonicalization sorts steps so a bound complex becomes the reference form.
- **Ideal reference (free enzyme, BFS reseeded):** still 14/69 — multi-valley topologies
  (free `E` and covalent `Estar` reachable only *through* a binding step) cannot avoid the
  division by any reference choice.
- Current `rate_equation` on the deduped set, count non-finite: reverse initial velocity
  (`NADH = Pyruvate = 0`) **65/69**; forward (`Lactate = NAD = 0`) 4/69; single
  `Pyruvate = 0` 48/69; single `NADH = 0` 27/69. Baseline (all conc > 0): all finite — the
  equations are correct for nonzero concentrations; only the zero-concentration evaluation
  breaks.

These show why both parts are needed: choosing/fixing the reference (Part B's
canonicalization) caps the symptom at the 14 irreducible cases; the reduction (Part A)
fixes all of them and removes the order-dependence entirely.

---

## Part A — Reduced polynomial form (eliminates division by concentration)

### Root cause

`rate_equation` and `rate_equation_string` both flow through `_raw_symbolic_rate_polys`
(`src/rate_eq_derivation.jl`). Within a rapid-equilibrium group, the Cha method expresses
each enzyme form's relative population **relative to a chosen reference form** — the first
form in step order. When that reference is a loaded complex, expressing the free enzyme
requires dividing out the bound substrates:

```
[E] / [E·NADH·Pyruvate]  =  (K_NADH · K_Pyruvate) / (NADH · Pyruvate)
```

The `normalize` step (`src/rate_eq_derivation.jl:373-389`, the `G == 1` case every
unicyclic init mechanism hits) bakes this in: it divides the numerator by `sigma_den[1]`
and each denominator term by `alpha_den[i]`, and those divisors carry concentrations
whenever the reference is loaded. At `[conc] = 0` numerator and denominator each blow up
and `∞/∞ = NaN` — a **removable** singularity (the limit is finite).

### The fix

A fraction's value is unchanged by multiplying numerator and denominator by the same
thing. The current code already divides (that is what `normalize` does) — it just divides
by the **wrong** factor (`sigma_den`, one form's population, which contains concentrations
the other terms don't all share, manufacturing the `1/[conc]` terms). Divide by the
**right** factor instead — the one common to every term:

1. **Drop `normalize`.** Build the denominator from `sigma_num[g]` for all groups (the
   path already used when `G > 1`) and leave the numerator undivided. So built, `N` and `D`
   are plain polynomials with **non-negative** concentration powers — but they share a
   common concentration factor (from the reference choice), which would give `0/0`.
2. **Reduce the pair to lowest terms over concentrations.** Divide both `N` and `D` by
   their concentration monomial GCD: for each metabolite symbol, the smallest power it has
   across all terms of `N ∪ D`; divide it out of every term. Because we only ever divide by
   the minimum power present, no term can go negative — no concentration denominator
   appears — and the least-loaded form's term becomes the constant denominator term that
   keeps the rate finite at `[conc] = 0`.

The GCD step is load-bearing, not cosmetic: skipping `normalize` *without* it leaves the
common factor and reproduces the original failures exactly (reverse 65/69). Reduction is
over **concentrations only**, never parameters — reducing a parameter common to every
term would silently drop a fitted parameter from the equation; concentration-only
reduction cannot (verified: 0 parameters dropped).

### Validation (deduped LDH set, 69 mechanisms)

- **0** rate-value mismatches vs. the current derivation for random conc > 0 (69 × 3 draws).
- **0** equations with division by concentration.
- **0** fitted parameters dropped.
- **0** non-finite for every realistic zero-pattern (each metabolite zeroed, forward,
  reverse).

### Code changes

- `src/sym_poly_for_rate_eq_derivation.jl`: add
  `_reduce_conc_lowest_terms(num::POLY, den::POLY, conc::Set{Symbol}) -> (POLY, POLY)`.
- `src/rate_eq_derivation.jl` `_raw_symbolic_rate_polys` (5-arg, lines 333-399): delete the
  `normalize` flag and `_poly_div_mono(num, sigma_den[1])`; the denominator loop always uses
  `sigma_num[g]`; final step (after the `rename_map` pass) call
  `_reduce_conc_lowest_terms(num, den, conc_set)` with the mechanism's metabolite symbols.
- `src/rate_eq_derivation.jl` `_compute_alpha` (250-310): drop the now-dead `sigma_den`
  (its only uses were the deleted `normalize` lines); return `(alpha_num, alpha_den, sigma_num)`.

### kcat is unaffected

`_kcat_forward` (`src/rate_eq_derivation.jl:699`) groups by concentration pattern and takes
the **ratio** `num_k/den_k` at matching patterns. Skipping `normalize` and reducing both
multiply `N`/`D` by a common monomial; concentration patterns shift uniformly (matching
preserved, patterns stay distinct) and the parameter parts cancel in the ratio. kcat is
invariant. The existing kcat / rescaling / scale-invariance tests are the gate.

### Allosteric path

`_allosteric_num_den_exprs` (`src/rate_eq_derivation.jl:1518`) takes its per-state
catalytic polys from `_raw_symbolic_rate_polys` (now reduced) and assembles the MWC
numerator/denominator by multiplication and addition only — which cannot introduce a
concentration denominator. The allosteric regression tests are the gate; a residual common
concentration factor in the combined form is not expected but would be caught there.

### Residual: all-metabolites-zero (out of scope)

14/69 still give `0/0` when *every* metabolite is `0` (multi-valley topologies with no
concentration-independent enzyme population, so no constant denominator term). There is
genuinely no steady state there and it is never a real data row; we document it and the
regression tests do not assert finiteness at all-zero.

---

## Part B — Canonical by construction + non-mutating dedup

### Change

Move the canonicalization currently performed by `_canonicalize_mechanism!`
(`src/mechanism_enumeration.jl:1670-1693` — sort steps within each group, sort groups by
their representative step, and for `AllostericMechanism` permute the parallel
`cat_allo_states` and sort `regulatory_sites`) **into the `Mechanism` and
`AllostericMechanism` constructors** (`src/types.jl`), alongside the iso-direction
canonicalization they already do (`_canonicalize_iso_groups`). Then:

- Every `Mechanism`/`AllostericMechanism` is canonical the moment it exists; writing the
  same steps in any order yields the **identical struct**.
- `_dedup_flat!` (`src/mechanism_enumeration.jl:1712`) collapses to `unique!(mechs)` —
  pure comparison, **no mutation**.
- `_canonicalize_mechanism!` is deleted.

### Why both Part A and Part B

Part A already makes the *rate* identical regardless of step order (the reduced form is
the canonical rational function), so "rewriting steps doesn't change the derivation" is
secured by Part A alone. Part B adds "rewriting steps doesn't change the *struct*"
(uniqueness), removes the dangerous in-place mutation from the dedup pass, and forecloses
a class of future order-dependence bugs. With Part A in place, the canonicalization need
only be *deterministic* (for dedup), not *derivation-aware* — the two concerns are fully
decoupled.

---

## Part C — Oracle bridge (canonizer info)

Hand-derived textbook oracles key parameters by **flat step position** (`K1`, `k2f`, …;
`positional_params`, `test/test_rate_eq_derivation.jl:113`). Once the constructor
canonicalizes step order, an oracle written in textbook (as-written) order would index the
wrong steps. Rather than rewrite oracles, **bridge** them:

- Oracles stay in as-written textbook order — the `@enzyme_mechanism` block remains the
  readable record of which step is "1", "2", ….
- Each oracle carries its as-written→canonical step **permutation** ("canonizer info"),
  recoverable by matching the as-written flat steps to the canonical `steps(m)` (`Step`
  compares structurally).
- `positional_params` applies the permutation so textbook `K1` resolves to the parameter
  of whatever canonical slot that step landed in.

No oracle formula changes; `positional_params` gains a permutation argument (or computes it
from a recorded as-written step list on the test spec). Exact plumbing is deferred to the
plan.

---

## Testing

1. **Shared fixtures** (`run_all_tests(spec)` over `MECHANISM_TEST_SPECS`, which already
   compiles `rate_equation` per fixture — true piggyback). For each fixture: random positive
   params and concentrations; for each metabolite, set just that one to `0.0` and assert
   `rate_equation` is **finite and nonzero** (the nonzero check guards the
   spurious-`0.0`-from-`1/Inf` regression). Covers allosteric fixtures too.
2. **Enumerated mechanisms** (folded into the existing enumeration loop). Run full
   `init_mechanisms` on a **small** reaction and `rate_equation`-test each (each metabolite
   zeroed → finite/nonzero), capping `rate_equation` compile cost by the established
   simplest-N-by-form-count pattern. **Open:** confirm what `bibi_ping_pong` refers to — a
   specific small fixture/reaction, or "a bi-bi reaction" (whose full init is 69, too heavy
   to compile in full → cap by form count). Will use the named fixture if one exists.
3. **Canonicalization** (`test/test_mechanism_enumeration.jl` / `test_types.jl`): the same
   mechanism written in two different step orders constructs to the **identical struct**;
   `_dedup_flat!` does not mutate its inputs (inputs equal their pre-call selves) and still
   collapses duplicates.
4. **Snapshot updates.** The 20 `expected_factored_num`/`expected_factored_denom` snapshots
   and the Expr-shape / flat-string regression tests change to the reduced forms; regenerate
   and eyeball. The textbook oracle tests must pass via the Part C bridge.
5. **Must stay green unchanged.** `test_rate_equation_performance` (allocation-free,
   sub-100 ns); kcat / rescaling / scale-invariance; Aqua/JET.

## Non-goals

- **Not** abandoning Cha rapid-equilibrium lumping for full steady-state King-Altman
  (parameter-count and compile-cost explosion). The `1/[conc]` artifact is inherent to the
  Cha intermediate algebra here; we cancel it in the reduced form rather than avoid forming it.
- **Not** guaranteeing finiteness at the all-metabolites-zero point.
- **Not** changing the parameter API, the thermodynamic-constraint machinery, or the
  fitting loop.

## TDD order (for the plan)

1. Add the zero-concentration tests (fixtures + small enumeration) — fail today.
2. Part A: add `_reduce_conc_lowest_terms`; rewire `_raw_symbolic_rate_polys`; drop
   `sigma_den` from `_compute_alpha`. Make the new tests pass; regenerate snapshots; confirm
   kcat/perf/Aqua/JET green.
3. Part B: add the canonicalization "same mechanism, two orders → identical struct" and
   "dedup does not mutate" tests — fail today. Move canonicalization into the constructors;
   reduce `_dedup_flat!` to `unique!`; delete `_canonicalize_mechanism!`.
4. Part C: make the textbook oracle tests pass through the permutation bridge.
5. Confirm on the LDH `identify_rate_equation` run that reverse-direction mechanisms now fit
   instead of dropping out.
