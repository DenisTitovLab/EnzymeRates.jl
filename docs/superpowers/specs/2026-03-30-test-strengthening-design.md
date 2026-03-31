# Mechanism Enumeration Test Strengthening & Self-Contained File

## Goal

Make `mechanism_enumeration.jl` self-contained (no dependency on old files), delete old enumeration pipelines and tests, and strengthen the new test suite with exact counts, explicit mechanism definitions, and comprehensive edge case coverage.

## Source Changes

### 1. Self-contained `mechanism_enumeration.jl`

Copy from `old_mechanism_enumeration.jl` into `mechanism_enumeration.jl`:

- **Types** (~100 lines): `ParamConstraint`, `AbstractMechanismSpec`, `StepSpec` (with `==`, `hash`), `MechanismSpec`, `AllostericMechanismSpec`
- **Step helpers** (~25 lines): `step_metabolite`, `step_forms`, `all_form_names`
- **Topology generation** (~450 lines): `_form_name`, `_atoms_dict`, `_can_pingpong`, `_subtract_atoms`, `_add_atoms`, `_catalytic_topologies`
- **Dead-end helpers** (~200 lines): `_bound_metabolites_at_forms`, `_dead_end_form_name`, `_substrate_product_dead_end_opportunities`, `_expand_substrate_product_dead_ends`
- **Compilation** (~150 lines): `_compile_enzyme_mechanism` (with `_clean_met`, form_name_map collision handling), wrapped as `EnzymeMechanism(spec::MechanismSpec)`

### 2. AllostericMechanismSpec dispatch for expansion moves

Add `AllostericMechanismSpec` methods to three expansion functions that currently only accept `MechanismSpec`. Each method delegates to `.base` and rewraps:

```julia
function _expand_re_to_ss(spec::AllostericMechanismSpec)
    [_rewrap_allosteric(spec, new_base)
     for new_base in _expand_re_to_ss(spec.base)]
end

function _expand_remove_constraint(spec::AllostericMechanismSpec)
    [_rewrap_allosteric(spec, new_base)
     for new_base in _expand_remove_constraint(spec.base)]
end

function _expand_add_dead_end_regulator(
    spec::AllostericMechanismSpec, reaction;
    exclude_regs=Set{Symbol}())
    # Merge allosteric reg names into exclude_regs
    allo_regs = Set{Symbol}()
    for site in spec.allosteric_reg_sites
        for lig in site; push!(allo_regs, lig); end
    end
    all_excluded = union(exclude_regs, allo_regs)
    [_rewrap_allosteric(spec, new_base)
     for new_base in _expand_add_dead_end_regulator(
         spec.base, reaction; exclude_regs=all_excluded)]
end
```

### 3. Simplify `_add_expansions!` for AllostericMechanismSpec

The allosteric `_add_expansions!` method currently has explicit `_rewrap_allosteric` calls for base-level moves. With the new dispatch methods, it simplifies to the same structure as the `MechanismSpec` version. This is a pure refactor — no behavior change.

`_expand_to_allosteric` already has an `AllostericMechanismSpec` method that returns empty (already-allosteric specs don't convert again). `_expand_add_allosteric_regulator` and `_expand_remove_tr_equiv` also already have `AllostericMechanismSpec` methods. With the new dispatch methods for `_expand_re_to_ss`, `_expand_remove_constraint`, and `_expand_add_dead_end_regulator`, ALL expansion functions accept both spec types via multiple dispatch.

This means we can collapse the two `_add_expansions!` methods into one that works for both types:

```julia
function _add_expansions!(result, spec::AbstractMechanismSpec, reaction)
    for s in _expand_re_to_ss(spec)
        _push_to_dict!(result, s)
    end
    for s in _expand_remove_constraint(spec)
        _push_to_dict!(result, s)
    end
    for s in _expand_add_dead_end_regulator(spec, reaction)
        _push_to_dict!(result, s)
    end
    for s in _expand_to_allosteric(spec, reaction)
        _push_to_dict!(result, s)
    end
    for s in _expand_add_allosteric_regulator(spec, reaction)
        _push_to_dict!(result, s)
    end
    for s in _expand_remove_tr_equiv(spec, reaction)
        _push_to_dict!(result, s)
    end
end
```

Each function returns empty for inapplicable types (e.g., `_expand_to_allosteric` returns `[]` for `AllostericMechanismSpec`, `_expand_add_allosteric_regulator` returns `[]` for `MechanismSpec`).

### 4. Delete old files

- Delete `src/old_mechanism_enumeration.jl`
- Delete `src/old_beam_enumeration.jl`
- Delete `test/old_test_mechanism_enumeration.jl`
- Delete `test/old_test_beam_enumeration.jl`
- Remove their includes from `src/EnzymeRates.jl` and `test/runtests.jl`

## Test Changes

### T1: Universal round-trip after every `@enzyme_mechanism`

Every mechanism defined with `@enzyme_mechanism` in the test file immediately gets a round-trip check:

```julia
m = @enzyme_mechanism begin ... end
spec = mechanism_spec_from_mechanism(m, rxn)
@test EnzymeMechanism(spec) === m
```

Remove the separate "Types and round-trip", "AllostericEnzymeMechanism round-trip", and "AllostericEnzymeMechanism round-trip with regulator" testsets.

### T2: Catalytic topologies — add ter-ter and ter-bi

Add reaction constants at top of file:
```julia
const ter_ter_rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[C], Q[N], R[X]
end

const ter_bi_rxn = @enzyme_reaction begin
    substrates: A[C], B[N], D[X]
    products: P[CN], Q[X]
end
```

Add testsets:
- Ter-Ter: `@test length(topos) == 3969`, all with 1 SS step
- Ter-Bi: `@test length(topos) == 204`, all with 1 SS step

With comments explaining derivation (counts verified by old test suite):
- Ter-ter: `(2^(3!) - 1)^2 = 63 × 63 = 3969`
- Ter-bi: `63 × 3 + 15 = 204` (189 sequential + 15 ping-pong)

### T3: init_mechanisms — add ter-ter/ter-bi, exact dead-end counts

Add ter-ter and ter-bi to "Param count invariant" (`n_s + n_p + 3`) and "All have exactly 1 SS step".

Replace "Dead-end counts" with exact per-topology counts using `@enzyme_mechanism`-defined mechanisms. Taken from old Stage 3a tests:

**Uni-uni** (passthrough):
```julia
# 3 forms: E, E_S[C], E_P[C]. E_S has all subs, E_P has
# all prods. No mixed sub+prod dead-end forms possible.
# → 0 dead-end forms, 1 variant (bare topology)
```

**Bi-bi random** (16 variants):
```julia
# 7 forms: E, E_A, E_B, E_A_B, E_P, E_Q, E_P_Q
# Eligible dead-end forms (mixed sub+prod binding):
#   E_A+P→E_A_P✓, E_A+Q→E_A_Q✓
#   E_B+P→E_B_P✓, E_B+Q→E_B_Q✓
#   E_P+A→E_A_P(same), E_P+B→E_B_P(same)
#   E_Q+A→E_A_Q(same), E_Q+B→E_B_Q(same)
# 4 unique dead-end forms → 2^4 = 16 variants
```

**Uni-bi ordered** (passthrough):
```julia
# 4 forms: E, E_S, E_P_Q, E_Q
# E+P→E_P would be single-product → rejected
# E_Q+S→E_S_Q has all subs → rejected
# → 0 dead-end forms, 1 variant
```

**Bi-bi ping-pong** (8 variants):
```julia
# Forms: E, E_A, Estar, Estar_A_P, Estar_B, E_Q
# E_A: +P→E_A_P(mixed✓), +Q→E_A_Q(mixed✓)
# E_Q: +B→E_B_Q(mixed✓)
# 3 dead-end forms → 2^3 = 8 variants
```

Remove "All compile correctly" testset.

### T4: Rename testsets — drop "Move N:" prefixes

- "Move 1: RE→SS conversion" → "RE→SS conversion"
- "Move 2: Remove equivalence constraint" → "Remove equivalence constraint"
- "Move 3: Add dead-end regulator" → "Add dead-end regulator"
- "Move 4: Add allosteric regulator" → "Add allosteric regulator"
- "Move 5: Remove TR equivalence" → "Remove TR equivalence"
- "Move 6: Allosteric conversion (+1)" → "Allosteric conversion"

### T5: RE→SS conversion — strengthen

All test mechanisms defined with `@enzyme_mechanism`, round-tripped.

- **All-SS with dead-end RE binding**: Define a mechanism where all catalytic steps are SS and a dead-end inhibitor binds to 2+ forms (so binding steps have K equivalence constraints). Verify `_expand_re_to_ss` yields nothing — the only RE steps are the constrained dead-end binding steps, which are skipped.
- **Allosteric mechanism**: Define an allosteric spec, call `_expand_re_to_ss` on it directly (new dispatch), verify results are `AllostericMechanismSpec` with correct allosteric fields preserved.
- **Bi-bi with exact count**: Define a bi-bi mechanism, verify exact number of results matches number of unconstrained RE steps. Comment explains which steps are eligible.
- **Mirror step param_count**: Existing test preserved (verifies +1 per converted mirror).

### T6: Remove equivalence constraint — strengthen

All test mechanisms defined with `@enzyme_mechanism`, round-tripped.

- **Mechanism with constraints**: Define a bi-bi random mechanism (has K_A and K_B constrained across forms). Verify exact count = number of removable constraint groups.
- **Same metabolite as substrate and regulator**: Define mechanism where S is both substrate and dead-end inhibitor. Add S as regulator via `_expand_add_dead_end_regulator`. Verify that substrate/product K constraints and regulator K constraints exist independently and removing one doesn't affect the other.
- **Multiple equivalence constraints**: Mechanism with 3+ constraint groups. Verify exact count and that each result has exactly one fewer constraint group.
- **SS constraint pairing**: Existing test preserved (kf/kr removed together).

### T7: Add dead-end regulator — strengthen

All test mechanisms defined with `@enzyme_mechanism`, round-tripped.

- **Two regulators, both bound**: Define reaction with two regulators. Add first, then add second. Verify forms exist where both are bound (mirror steps connecting E_I_J forms).
- **Bi-bi exact count**: Define bi-bi mechanism with one regulator. Verify exact count = `2^n_eligible - 1` where n_eligible is the number of eligible catalytic forms. Comment explains eligibility.
- **Allosteric mechanism**: Test on `AllostericMechanismSpec` via new dispatch. Allosteric-only regulator yields nothing. Dead-end-eligible regulator produces results.
- **Ping-pong mechanism**: Define ping-pong with regulator. Verify exact dead-end count.

### T8: All tests use `@enzyme_mechanism`, not init_mechanisms indexing

Every test (except catalytic topology counts and init_mechanisms param count invariant) defines its input mechanism explicitly. The pattern:

```julia
m = @enzyme_mechanism begin
    species: begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        enzymes: E, E_A[C], E_A_B[CN], E_B[N],
                 E_P[C], E_P_Q[CN], E_Q[N]
    end
    steps: begin
        [E, A] ⇌ [E_A]
        [E_B, A] ⇌ [E_A_B]
        [E, B] ⇌ [E_B]
        [E_A, B] ⇌ [E_A_B]
        [E, P] ⇌ [E_P]
        [E_P, Q] ⇌ [E_P_Q]
        [E, Q] ⇌ [E_Q]
        [E_Q, P] ⇌ [E_P_Q]
        [E_A_B] <--> [E_P_Q]
    end
end
spec = mechanism_spec_from_mechanism(m, bi_bi_rxn)
@test EnzymeMechanism(spec) === m  # round-trip
# ... use spec for tests
```

### T9: Tests for allosteric dispatch of base-level moves

New testset: "Base-level moves on allosteric mechanisms"

- **RE→SS on allosteric**: Define an allosteric spec. Call `_expand_re_to_ss(allo_spec)`. Verify results are `AllostericMechanismSpec`. Verify allosteric fields (reg_sites, tr_equiv, r_only) preserved. Verify param_count delta.
- **Remove constraint on allosteric**: Same pattern.
- **Add dead-end regulator on allosteric**: Same pattern. Verify allosteric regulators excluded from dead-end candidates.

## TDD Order

All implementation follows TDD:
1. Write failing test for allosteric dispatch methods
2. Implement the dispatch methods
3. Copy types and functions from old files
4. Verify all existing tests still pass
5. Strengthen tests per T1-T8
6. Delete old files
7. Final full test suite run
