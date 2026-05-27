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

## Stage 6.2 — commit 7af94e2

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

---

## Stage 7d.0 — commit TBD

### test_mechanism_enumeration.jl `@testset "AllostericEnzymeMechanism TR equivalence"`
- **Original**: built a `MechanismSpec` from a `StepSpec[]` literal,
  wrapped in an `AllostericMechanismSpec`, then compiled to verify
  `cat_allo_state` returns `:EqualRT` / `:OnlyR` for the seeded tags.
- **Why removed**: the test was rewritten (not deleted): same assertions
  now build the `AllostericEnzymeMechanism` directly via
  `@allosteric_mechanism` DSL with per-step `:: AlloState` annotations.
  Net `@test` count is unchanged; the test ID is the same. Listed here
  only because Stage 7d.0's spec-type retirement triggered the rewrite.

### test_mechanism_enumeration.jl `@testset "compile_mechanism round-trip"`
- **Original (4 sub-tests)**: built a `MechanismSpec` via
  `mechanism_spec_from_mechanism_and_rxn(m, rxn)` and asserted
  `EnzymeMechanism(spec) === m` (round-trip type identity through
  the spec path) for uni-uni, bi-bi sequential, bi-bi ping-pong, and
  uni-uni-with-regulator seeds.
- **Replacement**: NONE EQUIVALENT — the round-trip surface required
  preserving the iso step's source-order direction (`E_S <--> E_P`
  not `E_P <--> E_S`), which the `MechanismSpec`-via-`reactions(em)`
  path did naturally but the `Mechanism`-via-`Step` path does not
  (Step canonicalizes iso direction lex). The `===` type identity
  cannot survive the canonicalization.
- **Adjacent coverage**:
  - `compile_mechanism dispatch` (replacement testset) verifies the
    dispatch contract: `compile_mechanism(::Mechanism) ===
    EnzymeMechanism(::Mechanism)` and
    `compile_mechanism(::AllostericMechanism) ===
    AllostericEnzymeMechanism(::AllostericMechanism)`. The "build a
    Mechanism then lift" path is exercised by every `init_mechanisms`
    output fixture (lifted via `EnzymeMechanism(m)` for the
    rate-equation derivation tests).
  - `_canonical_rate_eq_hash_partition` in
    `test/test_canonical_hash_partition.jl` indirectly verifies that
    DSL-built mechanisms produce stable canonical hashes; any
    round-trip-induced semantic drift would surface as a hash
    mismatch.

## Stage 7d.1 — commit e179e69

Eleven `@testset "EnzymeReactionLegacy …"` blocks in `test/test_types.jl`
are deleted in this stage. The helper they exercised
(`EnzymeReactionLegacy{S,P,R,N}` singleton struct + outer constructor
+ `_sum_atoms` helper + `Base.show` method + 5 accessor methods +
`_to_legacy_reaction(::EnzymeReaction)` adapter) has no remaining
users in `src/` after Task 7d.0 ported the heavy enumeration pipeline
to dispatch on `EnzymeReaction` directly. Each surviving heading
below corresponds to one deleted testset, per
`scripts/check_test_integrity.sh` Check 2's one-heading-per-deletion
contract.

### test_types.jl `@testset "EnzymeReactionLegacy"`
- Basic construction + substrates/products/regulators accessors.
- Replacement: `test_types.jl @testset "EnzymeReaction (new concrete)"`
  exercises the parallel `EnzymeReaction(reactants, regulators, mults)`
  constructor and its `substrates` / `products` / `regulators`
  accessors.

### test_types.jl `@testset "EnzymeReactionLegacy with regulators"`
- Regulators-tuple construction + `regulators(r)` accessor.
- Replacement: `test_types.jl @testset "EnzymeReaction (new concrete)"`
  passes `RegulatorMults`-wrapped regulators and asserts
  `length(regulators(r)) == 1`.

### test_types.jl `@testset "EnzymeReactionLegacy with regulator roles"`
- `(name, :dead_end)` / `(name, :allosteric)` role pairs +
  `regulator_roles(r)` accessor.
- Replacement: behavior preserved structurally —
  `RegulatorMults.regulator` is one of `AllostericRegulator` /
  `CompetitiveInhibitor`; the role distinction is encoded in the
  subtype rather than a parallel `role::Symbol`. `EnzymeReaction
  (new concrete)` covers the allosteric case via
  `AllostericRegulator(:cAMP)`. The `regulator_roles` accessor itself
  is deleted; no caller in `src/` remains after Task 7d.0.

### test_types.jl `@testset "EnzymeReactionLegacy canonical ordering"`
- Constructor sorts substrate/product tuples by name.
- Replacement: `test_types.jl @testset "EnzymeReaction canonicalizes
  reactant + regulator ordering"` asserts the same canonicalization
  for `EnzymeReaction.reactants` / `.regulators`.

### test_types.jl `@testset "EnzymeReactionLegacy regulator same as substrate allowed"`
- Constructor accepts a regulator name that overlaps a substrate name.
- Replacement: NONE EQUIVALENT as a standalone assertion. Behavior
  preserved structurally — `EnzymeReaction`'s constructor does not
  validate against substrate-regulator name overlap (the `regulators`
  vector is independent from the `reactants` vector and carries no
  uniqueness constraint vs reactants). User-reachable through
  `@enzyme_reaction`, no error raised.

### test_types.jl `@testset "EnzymeReactionLegacy regulator same as product allowed"`
- Same as above, regulator overlaps a product name.
- Replacement: NONE EQUIVALENT — same reasoning as the substrate
  case above. The concrete struct does not validate this overlap.

### test_types.jl `@testset "EnzymeReactionLegacy duplicate substrate names"`
- Constructor errors on duplicate substrate names.
- Replacement: NONE EQUIVALENT — `EnzymeReaction`'s constructor does
  not enforce name-uniqueness in `reactants`. Domain convention is
  upheld by the DSL parser at higher levels: `@enzyme_reaction begin
  substrates: A, A; … end` produces duplicate `Substrate(:A)`
  reactant entries and the user gets a later error from form-name
  collision in the enumerator. Not user-reachable as a constructor
  failure mode in the public API.

### test_types.jl `@testset "EnzymeReactionLegacy duplicate product names"`
- Constructor errors on duplicate product names.
- Replacement: NONE EQUIVALENT — same reasoning as the duplicate
  substrate case above.

### test_types.jl `@testset "EnzymeReactionLegacy oligomeric_state"`
- Default + explicit `oligomeric_state` kw, type-parameter equality.
- Replacement: `test_dsl.jl @testset "@enzyme_reaction with
  oligomeric_state"` exercises the DSL's `oligomeric_state: N` label
  and asserts the produced `EnzymeReaction`'s
  `allowed_catalytic_multiplicities == [N]`. Type-parameter equality
  is meaningless for the concrete struct (`EnzymeReaction` is
  non-parametric) — the surviving structural equality is `==` on the
  field-by-field content.

### test_types.jl `@testset "EnzymeReactionLegacy: atom mandatory"`
- Constructor errors on empty atom tuple per metabolite.
- Replacement: NONE EQUIVALENT — `EnzymeReaction`'s constructor
  allows empty `ReactantAtoms.atoms`. Mandatory-atom enforcement
  lives in the `@enzyme_reaction` DSL macro: the grammar requires
  bracket syntax (`S[C]` / `S[C6H12O6]`), and bare `S` parses as a
  syntax error. The atoms-mandatory contract is therefore preserved
  at the DSL boundary, where every user actually creates reactions.

### test_types.jl `@testset "EnzymeReactionLegacy: atom balance"`
- Constructor errors when declared atoms don't balance across sides.
- Replacement: NONE EQUIVALENT — `EnzymeReaction`'s constructor does
  not balance-check declared atoms. In the concrete-types
  architecture the atoms field is presentation metadata for
  downstream consumers (e.g. ping-pong residual computation in
  `_catalytic_topologies` via `_atoms_dict` / `_can_pingpong`).
  Balance violations now manifest as `_can_pingpong` returning false
  on specific topologies (the enumerator silently skips
  infeasible ping-pong paths) rather than as constructor errors.

## Finish-refactor Phase 1 — commit TBD (backfill SHA)

### test_mechanism_definitions_for_test_enzyme_derivation.jl `MechanismTestSpec name="Segel Theorell-Chance Bi Bi"`
- **Original**: a Theorell-Chance Bi Bi mechanism `E + A ⇌ EA + B ⇌ EQ + P ⇌ E + Q`
  with an analytical rate (Segel Eq. IX-122) and an analytical kcat
  (`p -> p.k3f`). Tested on `main`.
- **Why removed**: un-representable in the concrete-types grammar. The middle
  step `EA + B <--> EQ + P` *binds* B **and** *releases* P in a single
  transition; `Step` carries one `bound_metabolite` field, so a compound
  bimolecular bind+release cannot be encoded. Both the structs-throughout and
  the concrete-types refactor attempts confirmed this — there is no decomposed
  Species rename that preserves the lumped Theorell-Chance step graph. Denis
  approved deletion (brainstorm 2026-05-25).
- **Replacement**: NONE EQUIVALENT — the lumped bind+release Theorell-Chance
  variant is intentionally not expressible.
- **Adjacent coverage**: the bi-bi mechanism space is covered by **Segel Ordered
  Bi Bi** (full ternary `E(A,B)` complex, sequential bind/release) and the
  **Ping-Pong Bi Bi** family (covalent `F` intermediate). Theorell-Chance is the
  zero-ternary-complex limit between them; its distinguishing feature (no
  measurable central complex) is a rate-law simplification, not a new
  elementary-step topology, so derivation coverage is not reduced.

## Concrete-types refactor — legacy 2-arg EnzymeMechanism constructor removal — commit TBD

### test_types.jl `@test_throws` — stoichiometry rank mismatch (substrate "vanishes", no product)
- Helper deleted: the 2-arg `EnzymeMechanism(mets, rxns)` constructor's `_validate_*` (removed in the next task).
- Replacement: NONE EQUIVALENT (superseded by design — net-stoichiometry rank is not enforced on the decomposed path). `_assert_mechanism_invariants` only checks declared-metabolite coverage, not a stoichiometry rank balance.
- Adjacent coverage: declared substrate/product coverage is re-pointed to `_assert_mechanism_invariants` (see the net-stoich-mismatch re-point in `@testset "EnzymeMechanism error cases"`); thermodynamic/rank concerns are a Wegscheider/Haldane downstream concern.

### test_types.jl `@test_throws` — iso group size > 1
- Helper deleted: the 2-arg `EnzymeMechanism(mets, rxns)` constructor's `_validate_*`.
- Replacement: NONE EQUIVALENT (superseded). Verified empirically: `_assert_mechanism_invariants` does NOT error on a kinetic group containing two iso (nothing-bound) steps — its size>1 group check filters to `bound_metabolite !== nothing` steps, so two iso steps in one group pass.
- Adjacent coverage: per-step bound/iso consistency is still enforced by `_assert_mechanism_invariants`; the binding-step group composition (single metabolite, no RE/SS mixing) is re-pointed (see `@testset "Kinetic-group validator error paths"`).

### test_types.jl `@test_throws` — duplicate substrate names
- Helper deleted: the 2-arg `EnzymeMechanism(mets, rxns)` constructor's `_validate_*`.
- Replacement: NONE EQUIVALENT. Verified empirically: the decomposed `EnzymeReaction` constructor accepts duplicate `ReactantAtoms` substrate names without error (no `allunique` check). Not a constructor failure mode on the concrete-types path.
- Adjacent coverage: none — uniqueness is not a contract of the concrete `EnzymeReaction`.

### test_types.jl `@test_throws` — reaction with zero enzymes on a side
- Helper deleted: the 2-arg `EnzymeMechanism(mets, rxns)` constructor's `_validate_*`.
- Replacement: NONE EQUIVALENT (moot) — a decomposed `Step` always carries exactly one `Species` per side (`from_species` / `to_species`), so "zero enzymes on a side" is unrepresentable.

### test_types.jl `@test_throws` — two metabolites on one side
- Helper deleted: the 2-arg `EnzymeMechanism(mets, rxns)` constructor's `_validate_*`.
- Replacement: NONE EQUIVALENT (moot) — `Step.bound_metabolite` is a single `Union{Metabolite,Nothing}` field, so two free metabolites on one side of a single step is unrepresentable.

### test_types.jl `@test_throws` — unknown metabolite in reaction
- Helper deleted: the 2-arg `EnzymeMechanism(mets, rxns)` constructor's `_validate_*`.
- Replacement: NONE EQUIVALENT (moot) — a decomposed `Step` binds a typed `Metabolite` object (`Substrate`/`Product`/…), not a free Symbol, so an "unknown" Symbol metabolite cannot enter a step.

### test_types.jl `@testset "Connectivity validator: orphan enzyme form → error"`
- Helper deleted: the 2-arg `EnzymeMechanism(mets, rxns)` constructor's connectivity validator.
- Replacement: NONE EQUIVALENT (superseded) — the decomposed design intentionally dropped the enzyme-form connectivity invariant (documented in the surviving `test_types.jl @testset "EnzymeMechanism error cases"` NOTE comment): enzyme forms are inferred from steps, so a disjoint subgraph is structurally valid.
- Adjacent coverage: graph connectivity is a downstream Wegscheider-analysis concern, not a constructor invariant.

### test_types.jl `@test_throws` — non-consecutive kinetic_group numbers (AllostericEnzymeMechanism validators)
- Helper deleted: the direct `EnzymeMechanism{(mets, rxns)}()` singleton form (`cm_bad` fixture) that encoded a group "hole" (groups 1, 3).
- Replacement: NONE EQUIVALENT (moot) — the decomposed `Mechanism.steps` outer-vector index *is* the kinetic-group number, so groups are dense-by-construction and a non-consecutive (holed) kinetic-group numbering is unrepresentable. The constructor can never receive such a mechanism.
- Adjacent coverage: `cat_allo_states` length validation against the cat-group count is still exercised by the surviving wrong-length `@test_throws` in the same testset.

### test_types.jl `@test_throws` — strict regulator binding (regulator listed but never bound)
- Helper deleted: the 2-arg `EnzymeMechanism(mets, rxns)` constructor's strict-regulator-binding validator.
- Replacement: NONE EQUIVALENT (superseded). Regulators are intentionally excluded from the `_assert_mechanism_invariants` coverage check (per CLAUDE.md: `init_mechanisms` declares dead-end inhibitors that no step binds yet; `expand_mechanisms` binds them later). Verified empirically: an unbound declared regulator on the `@enzyme_mechanism` path produces NO error — the regulator is simply absent from `regulators(m)`.
- Adjacent coverage: substrate/product coverage (the never-excluded case) is re-pointed to `_assert_mechanism_invariants`; positive "all regulators bound" / "no regulators" cases are preserved via decomposed `@enzyme_mechanism` in `@testset "EnzymeMechanism: regulator binding"`.
