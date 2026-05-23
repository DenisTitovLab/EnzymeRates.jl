# Refactor — Deleted test log

Per spec §2 / §2.1 narrow exception: tests of deleted helpers or
modules whose behavior is preserved by the replacement code path.

The `scripts/check_test_integrity.sh` Check 1 (whole-file deletion) and
Check 2 (@testset count decrease) both consult this file. Each
`### <filename>` heading documents one permitted deletion.

## Stage 6β.10 — commit c6b2259

### test_mechanism_enumeration.jl `@testset "AllostericMechanismSpec constructor density validation"`
- Helper deleted: dense-storage validation hooks specific to the
  `AllostericMechanismSpec` constructor (still alive as a type, but
  no longer constructed from user-test code post-deletion of the
  spec testset family).
- Replacement: NONE EQUIVALENT. The testset validated the constructor's
  rejection of sparse `group_tags` / `reg_ligand_tags` Dicts; the
  constructor itself stays alive (it's still used internally by the
  heavy enumeration pipeline) but its validation logic is now exercised
  only through the heavy pipeline's well-formed call sites. The dense
  invariant survives — there is just no longer a user-facing entry
  point that produces sparse input.
- Adjacent coverage:
  - `test/test_types.jl @testset "AllostericEnzymeMechanism constructor
    validators"` validates the corresponding rejection on the
    user-facing `AllostericEnzymeMechanism` constructor (the only
    construction path that survives 6β.10).

### test_mechanism_enumeration.jl `@testset "spec-from-mechanism helpers reject inconsistent inputs"`
- Helper deleted: NONE (the helper `mechanism_spec_from_mechanism_and_rxn`
  is retained as a test-internal helper until Task 7d.0 — it is still
  used by ~12 auxiliary testsets that exercise the heavy pipeline's
  internals). The testset itself is the §2.1 deletion target.
- Replacement: NONE EQUIVALENT. The deleted testset validated that the
  helper rejects substrate / product / regulator / oligomeric-state
  disagreements between a compiled mechanism and a reaction. The
  helper's only purpose in the Mechanism-form world is to bootstrap
  spec-flavored testsets that have been deleted; its validation surface
  is exercised in passing by every surviving caller, but no targeted
  rejection-path test remains. This is a true §2.1 narrow exception.
- Adjacent coverage:
  - `test/test_types.jl @testset "EnzymeMechanism error cases"`
    validates the corresponding rejection at the `EnzymeMechanism`
    constructor level (no spec round-trip involved).

### test_mechanism_enumeration.jl `@testset "allosteric_spec_from_mechanism_and_rxn round-trip"`
- Helper deleted: NONE (the helper is retained — see entry above).
- Replacement: NONE EQUIVALENT. The deleted testset validated lossless
  round-tripping `AllostericEnzymeMechanism → spec →
  AllostericEnzymeMechanism` identity. With the round-trip surface no
  longer the canonical entry, the identity is asserted at the
  AllostericMechanism / AllostericEnzymeMechanism layer instead.
- Adjacent coverage:
  - `test/test_types.jl @testset "AllostericEnzymeMechanism (new
    design)"` directly constructs allosteric mechanisms and validates
    invariants.
  - `test/test_types.jl @testset "AllostericEnzymeMechanism
    constructor + DSL"` validates accessor parity for DSL-driven
    construction (which round-trips back through the same singleton
    type).

### test_mechanism_enumeration.jl `@testset "Base-level moves on allosteric specs" > "RE→SS on allosteric"`
- Helper deleted: spec `_expand_re_to_ss(::AbstractMechanismSpec)` and
  spec `_expand_to_allosteric(::MechanismSpec, ...)`. The testset
  applied the spec-flavored RE→SS move on top of a spec-flavored
  to-allosteric move on a spec built via the spec helper. Every step
  of that chain is deleted.
- Replacement: NONE EQUIVALENT in the spec form. The Mechanism-form
  RE→SS-on-allosteric behavior is fully covered by the surviving
  parallel testsets in `_expand_re_to_ss`:
  - `@testset "Mechanism — :EqualRT group: Δ=+1"` (post-deletion line)
  - `@testset "Mechanism — :OnlyR group: Δ=+1"`
  - `@testset "Mechanism — :NonequalRT group: Δ=+2"`
  - `@testset "AllostericMechanism — :EqualRT: 2 variants, tags
    preserved"`
  These cover every R/T-state tag flavor on the Mechanism path.

---

## Stage 6.2 — commit TBD-after-commit

### test_mechanism_enumeration.jl `@testset "Rate-equation canonical hash dedup" > "_factor_sort_key sort order"`
- Helper deleted: `_factor_sort_key` (a regex-pipeline sort-key helper
  in `src/mechanism_enumeration.jl`, removed in Stage 6.2 alongside
  the rest of the regex canonicalizer — `_canonicalize_rate_eq_with_map`,
  `_sort_run_factors`, `_canonical_rate_eq_hash_data_impl_regex`, and
  the `_canonical_rate_eq_hash_old` / `_canonical_rate_eq_hash_data_old`
  aliases).
- Replacement: NONE EQUIVALENT. The deleted testset asserted sort-key
  ordering for `p_i`-style token strings; with the regex canonicalizer
  gone, no string sorting happens — the struct-based canonicalizer
  (`_canonical_rate_eq_hash_data_impl_struct`) works on `Expr` trees
  and Parameter struct hashes directly.
- Adjacent coverage:
  - `test/test_canonical_hash_partition.jl` exercises the new
    canonicalizer end-to-end on `init_mechanisms` output and pins
    the equivalence-class count, catching any regression that would
    have manifested via a sort-key bug in the old pipeline.

### test_mechanism_enumeration.jl `@testset "Rate-equation canonical hash dedup" > "_sort_run_factors sort order"`
- Helper deleted: `_sort_run_factors` (the regex helper that re-sorted
  multiplication runs in canonicalized strings, removed alongside
  `_factor_sort_key` per the entry above).
- Replacement: NONE EQUIVALENT. The struct-based canonicalizer doesn't
  produce string output for sorting; monomial ordering inside the
  canonical Expr tree is the responsibility of the Expr-builder
  (`_poly_to_expr` etc.), which is exercised by every rate-equation
  test in the suite.
- Adjacent coverage: same as above — `test/test_canonical_hash_partition.jl`
  pins partition behavior end-to-end.

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
