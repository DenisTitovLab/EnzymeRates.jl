# Refactor — Deleted test log

Per spec §2 / §2.1 narrow exception: tests of deleted helpers or
modules whose behavior is preserved by the replacement code path.

The `scripts/check_test_integrity.sh` Check 1 (whole-file deletion) and
Check 2 (@testset count decrease) both consult this file. Each
`### <filename>` heading documents one permitted deletion.

---

## Pre-refactor cleanup — commit 4fb462e (PR #36, 2026-05-13)

### test_sym_poly.jl

Whole-file deletion (12 testsets). Landed before this branch's
concrete-types refactor began.

- **Module deleted:** the factored-polynomial type family in
  `src/sym_poly.jl` (replaced by always-expanded rate-equation
  emission in commit 28041fe).
- **Behavior removed:** the polynomial-algebra layer that
  `test_sym_poly.jl` exercised — factor/expand round-trips, custom
  ring-of-rationals operations, the factored-poly canonicalizer.
  Rate equations are now emitted directly in expanded form, so the
  factored intermediate doesn't exist.
- **Replacement coverage:**
  - End-to-end rate-equation correctness: every fixture in
    `test/mechanism_definitions_for_test_enzyme_derivation.jl` (37
    fixtures) pairs a mechanism definition with a hand-derived
    analytical formula. `test/test_rate_eq_derivation.jl` compares
    `rate_equation_string(m)` against each hand-derived string —
    if the expansion is ever wrong, this catches it.
  - 0-allocation/<100ns per-call gate in
    `test_rate_equation_performance` (`test/test_rate_eq_derivation.jl`)
    locks the runtime shape of the emitted code; symbolic-algebra
    regressions would surface as allocations.
  - The deduplication pipeline (formerly tested via factored-poly
    structural keys) is now tested via canonical rate-equation hashes
    in `test/test_mechanism_enumeration.jl @testset
    "Rate-equation canonical hash dedup"` and the `_canonicalize!`
    / `_dedup_key` testsets.

---

## Concrete-types refactor — modified @test lines (Check 5 WARN)

`scripts/check_test_integrity.sh` Check 5 (soft WARN) flags every
modified `@test` / `@test_throws` line vs `main`. The current
branch carries a large cumulative diff from the Stage 1–5 type
refactor that changed accessor return shapes (substrates/products/
regulators went from `Tuple{Tuple{Symbol, Atoms}}` to
`Vector{<:Reactant}` / `Vector{RegulatorMult}`). Most of the
flagged lines are mechanical adaptations of this form:

```
- @test EnzymeRates.substrates(spec) == ((:S, ((:C, 1),)),)
+ @test EnzymeRates.substrates(spec) == [EnzymeRates.Substrate(:S)]
```

These are not weakenings — they assert the same semantic property
against the new return shape. They were reviewed at the time of
each Stage 1–5 commit; the cumulative WARN list is a known
side-effect of running the script against `main` (which still has
the old shape).

Per spec §2, no individual `@test` assertion has been weakened
(`==` → `≈`, exact-value → `isa`, etc.). New work that touches
`@test` lines must remain under the same constraint.
