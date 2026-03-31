# Test Strengthening & Self-Contained File Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `mechanism_enumeration.jl` self-contained, delete old enumeration files, and comprehensively strengthen tests.

**Architecture:** Copy ~900 lines from old files into new file, add allosteric dispatch methods for base-level moves, unify `_add_expansions!`, rewrite tests with explicit `@enzyme_mechanism` definitions and exact counts.

**Tech Stack:** Julia, EnzymeRates.jl

**Spec:** `docs/superpowers/specs/2026-03-30-test-strengthening-design.md`

---

## File Map

| File | Action | What |
|------|--------|------|
| `src/mechanism_enumeration.jl` | Modify | Copy types/functions from old files, add allosteric dispatch, unify `_add_expansions!` |
| `test/test_mechanism_enumeration.jl` | Rewrite | Strengthen all tests per spec |
| `src/old_mechanism_enumeration.jl` | Delete | No longer needed |
| `src/old_beam_enumeration.jl` | Delete | No longer needed |
| `test/old_test_mechanism_enumeration.jl` | Delete | Replaced by new tests |
| `test/old_test_beam_enumeration.jl` | Delete | Replaced by new tests |
| `src/EnzymeRates.jl` | Modify | Remove old includes |
| `test/runtests.jl` | Modify | Remove old test includes |

---

### Task 1: Add AllostericMechanismSpec dispatch methods (TDD)

**Files:**
- Modify: `src/mechanism_enumeration.jl`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing tests for allosteric dispatch**

Add to `test/test_mechanism_enumeration.jl`:

```julia
@testset "Base-level moves on allosteric specs" begin
    # Create an allosteric spec from uni-uni
    m_uu = @enzyme_mechanism begin
        species: begin
            substrates: S[C]
            products: P[C]
            enzymes: E, E_P[C], E_S[C]
        end
        steps: begin
            [E, P] ⇌ [E_P]
            [E, S] ⇌ [E_S]
            [E_S] <--> [E_P]
        end
    end
    spec = mechanism_spec_from_mechanism(m_uu, uni_uni_allo)
    @test EnzymeMechanism(spec) === m_uu
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo)
    allo = first(allo_specs)

    @testset "RE→SS on allosteric" begin
        result = EnzymeRates._expand_re_to_ss(allo)
        @test !isempty(result)
        for r in result
            @test r isa AllostericMechanismSpec
            # Allosteric fields preserved
            @test r.catalytic_n == allo.catalytic_n
            @test r.allosteric_reg_sites ==
                allo.allosteric_reg_sites
            @test r.r_only_metabolites ==
                allo.r_only_metabolites
            @test r.param_count == allo.param_count + 1
        end
    end

    @testset "Remove constraint on allosteric" begin
        # Need a spec with constraints — use bi-bi
        m_bb = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        bi_bi_allo = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        bb_spec = mechanism_spec_from_mechanism(
            m_bb, bi_bi_allo)
        @test EnzymeMechanism(bb_spec) === m_bb
        bb_spec_c = EnzymeRates.MechanismSpec(
            bb_spec.reaction, bb_spec.steps,
            EnzymeRates._max_equivalence_constraints(
                bb_spec),
            bb_spec.param_count)
        bb_allo = first(
            EnzymeRates._expand_to_allosteric(
                bb_spec_c, bi_bi_allo))
        result = EnzymeRates._expand_remove_constraint(
            bb_allo)
        @test !isempty(result)
        for r in result
            @test r isa AllostericMechanismSpec
            @test r.catalytic_n == bb_allo.catalytic_n
        end
    end

    @testset "Add dead-end reg on allosteric" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            oligomeric_state: 2
        end
        spec_i = mechanism_spec_from_mechanism(
            m_uu, rxn)
        allo_i = first(
            EnzymeRates._expand_to_allosteric(
                spec_i, rxn))
        result =
            EnzymeRates._expand_add_dead_end_regulator(
                allo_i, rxn)
        @test !isempty(result)
        for r in result
            @test r isa AllostericMechanismSpec
            @test r.catalytic_n == allo_i.catalytic_n
        end
    end
end
```

- [ ] **Step 2: Run to verify failure**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — no method `_expand_re_to_ss(::AllostericMechanismSpec)`.

- [ ] **Step 3: Implement dispatch methods**

In `src/mechanism_enumeration.jl`, add after each `MechanismSpec` method:

```julia
function _expand_re_to_ss(spec::AllostericMechanismSpec)
    [_rewrap_allosteric(spec, new_base)
     for new_base in _expand_re_to_ss(spec.base)]
end

function _expand_remove_constraint(
    spec::AllostericMechanismSpec)
    [_rewrap_allosteric(spec, new_base)
     for new_base in
         _expand_remove_constraint(spec.base)]
end

function _expand_add_dead_end_regulator(
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction);
    exclude_regs::Set{Symbol}=Set{Symbol}())
    allo_regs = Set{Symbol}()
    for site in spec.allosteric_reg_sites
        for lig in site
            push!(allo_regs, lig)
        end
    end
    all_excluded = union(exclude_regs, allo_regs)
    [_rewrap_allosteric(spec, new_base)
     for new_base in _expand_add_dead_end_regulator(
         spec.base, reaction;
         exclude_regs=all_excluded)]
end
```

- [ ] **Step 4: Unify `_add_expansions!`**

Replace both `_add_expansions!` methods with a single one dispatching on `AbstractMechanismSpec`:

```julia
function _add_expansions!(result, spec::AbstractMechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
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

Delete the old two-method `_add_expansions!` and the explicit `_rewrap_allosteric` calls that were in the allosteric version.

- [ ] **Step 5: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add allosteric dispatch for base-level moves, unify _add_expansions!"
```

---

### Task 2: Make `mechanism_enumeration.jl` self-contained

**Files:**
- Modify: `src/mechanism_enumeration.jl`

- [ ] **Step 1: Copy types from `old_mechanism_enumeration.jl`**

Copy to the TOP of `mechanism_enumeration.jl` (after ABOUTME):
- `ParamConstraint` type alias (line 11)
- `AbstractMechanismSpec` abstract type (line 15)
- `StepSpec` struct with `==` and `hash` (lines 18-31)
- `MechanismSpec` struct (lines 40-45)
- `AllostericMechanismSpec` struct with 10 fields (lines 60-71)
- `step_metabolite`, `step_forms` helpers (lines 75-79)
- `all_form_names` both methods (lines 82-98)

- [ ] **Step 2: Copy topology generation**

Copy after the types:
- `_form_name` (line 356)
- `_atoms_dict` (line 370)
- `_can_pingpong` (line 393)
- `_subtract_atoms` (line 404)
- `_add_atoms` (line 417)
- `_catalytic_topologies` (lines 435-862)

- [ ] **Step 3: Copy dead-end helpers**

- `_bound_metabolites_at_forms` (line 925)
- `_dead_end_form_name` (line 979)
- `_substrate_product_dead_end_opportunities` (line 1122)
- `_expand_substrate_product_dead_ends` (lines 1001-1335)

- [ ] **Step 4: Copy compilation**

- `_compile_enzyme_mechanism` (lines 2067-2215) including `_clean_met` and form_name_map
- Ensure `EnzymeMechanism(spec::MechanismSpec) = _compile_enzyme_mechanism(spec)` is present

- [ ] **Step 5: Run tests to verify nothing broke**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

All existing tests must pass with the duplicated definitions (Julia allows redefining methods in the same module — last definition wins).

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "Copy types and functions into mechanism_enumeration.jl for self-containment"
```

---

### Task 3: Delete old files

**Files:**
- Delete: `src/old_mechanism_enumeration.jl`, `src/old_beam_enumeration.jl`
- Delete: `test/old_test_mechanism_enumeration.jl`, `test/old_test_beam_enumeration.jl`
- Modify: `src/EnzymeRates.jl`, `test/runtests.jl`

- [ ] **Step 1: Remove includes**

In `src/EnzymeRates.jl`, remove:
```julia
include("old_mechanism_enumeration.jl")
include("old_beam_enumeration.jl")
```

In `test/runtests.jl`, remove:
```julia
include("old_test_mechanism_enumeration.jl")
include("old_test_beam_enumeration.jl")
```

- [ ] **Step 2: Delete old files**

```bash
git rm src/old_mechanism_enumeration.jl
git rm src/old_beam_enumeration.jl
git rm test/old_test_mechanism_enumeration.jl
git rm test/old_test_beam_enumeration.jl
```

- [ ] **Step 3: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

If any test fails, it means a dependency was missed in Task 2. Fix and re-run.

- [ ] **Step 4: Commit**

```bash
git add src/EnzymeRates.jl test/runtests.jl
git commit -m "Delete old enumeration pipelines and tests"
```

---

### Task 4: Rewrite test file — reactions, round-trips, topology tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Add reaction constants**

At top of file, ensure all needed reactions are defined:

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

Plus ensure `uni_uni_rxn`, `uni_bi_rxn`, `bi_bi_rxn`, `bi_bi_pp_rxn`, `uni_uni_allo` are present (most already exist).

- [ ] **Step 2: Remove old round-trip testsets**

Delete "Types and round-trip", "AllostericEnzymeMechanism round-trip", and "AllostericEnzymeMechanism round-trip with regulator" testsets. Their coverage is replaced by universal round-trips after every `@enzyme_mechanism`.

- [ ] **Step 3: Add universal round-trip pattern**

For every existing `@enzyme_mechanism` definition in the file, add immediately after:
```julia
spec = mechanism_spec_from_mechanism(m, rxn)
@test EnzymeMechanism(spec) === m
```

- [ ] **Step 4: Add Ter-Ter and Ter-Bi topology tests**

Add to "Catalytic topologies" testset:

```julia
@testset "Ter-Ter" begin
    topos = EnzymeRates._catalytic_topologies(ter_ter_rxn)
    # 3969 = 63 × 63 = (2^(3!) - 1)²
    # Each side (binding/release) has 3!=6 permutation
    # paths through Boolean lattice B_3; all 2^6-1=63
    # non-empty path subsets produce distinct edge sets;
    # sides are independent.
    @test length(topos) == 3969
    for t in topos
        @test count(!s.is_equilibrium for s in t.steps) == 1
    end
end

@testset "Ter-Bi" begin
    topos = EnzymeRates._catalytic_topologies(ter_bi_rxn)
    # 204 = 189 sequential + 15 ping-pong
    # Sequential: (2^(3!) - 1) × (2^(2!) - 1) = 63 × 3
    # Ping-pong: D[X]→Q[X] can isomerize independently
    @test length(topos) == 204
    for t in topos
        @test count(!s.is_equilibrium for s in t.steps) == 1
    end
end
```

- [ ] **Step 5: Add ter-ter/ter-bi to init_mechanisms tests**

Add to "Param count invariant" and "All have exactly 1 SS step":
```julia
(ter_ter_rxn, 3, 3),
(ter_bi_rxn, 3, 2),
```

- [ ] **Step 6: Run tests, commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "Add ter-ter/ter-bi tests, universal round-trips, remove redundant testsets"
```

---

### Task 5: Rewrite dead-end count tests with exact values

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Replace "Dead-end counts" testset**

Replace the existing naive `> N` tests with exact per-topology tests using `@enzyme_mechanism` definitions. Copy the mechanism definitions and expected counts from old Stage 3a tests:

```julia
@testset "Dead-end substrate/product expansion" begin

    @testset "Uni-Uni: no dead-end forms" begin
        # 3 forms: E, E_S[C], E_P[C]. E_S has all subs,
        # E_P has all prods. No mixed dead-end possible.
        # → 0 dead-end forms, 1 variant (bare topology)
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, E_P[C], E_S[C]
            end
            steps: begin
                [E, P] ⇌ [E_P]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P]
            end
        end
        spec = mechanism_spec_from_mechanism(m, uni_uni_rxn)
        @test EnzymeMechanism(spec) === m
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], uni_uni_rxn)
        @test length(result) == 1
    end

    @testset "Bi-Bi random: 4 dead-end forms" begin
        # 7 forms: E, E_A, E_B, E_A_B, E_P, E_Q, E_P_Q
        # E_A: +P→E_A_P(mixed✓), +Q→E_A_Q(mixed✓)
        # E_B: +P→E_B_P(mixed✓), +Q→E_B_Q(mixed✓)
        # E_P: +A→E_A_P(same), +B→E_B_P(same)
        # E_Q: +A→E_A_Q(same), +B→E_B_Q(same)
        # 4 unique dead-end forms → 2^4 = 16 variants
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products: P[C], Q[N]
                enzymes: E, E_A[C], E_A_B[CN],
                    E_B[N], E_P[C], E_P_Q[CN], E_Q[N]
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
        @test EnzymeMechanism(spec) === m
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], bi_bi_rxn)
        @test length(result) == 16
    end

    @testset "Uni-Bi ordered: no dead-end forms" begin
        # 4 forms: E, E_S, E_P_Q, E_Q
        # E+P→E_P: single-product → rejected
        # E_Q+S→E_S_Q: has all subs → rejected
        # → 0 dead-end forms, 1 variant
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[AB]
                products: P[A], Q[B]
                enzymes: E, E_P_Q[AB], E_Q[B], E_S[AB]
            end
            steps: begin
                [E, Q] ⇌ [E_Q]
                [E_Q, P] ⇌ [E_P_Q]
                [E, S] ⇌ [E_S]
                [E_S] <--> [E_P_Q]
            end
        end
        spec = mechanism_spec_from_mechanism(m, uni_bi_rxn)
        @test EnzymeMechanism(spec) === m
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], uni_bi_rxn)
        @test length(result) == 1
    end

    @testset "Bi-Bi Ping-Pong: 3 dead-end forms" begin
        # Forms: E, E_A, Estar, Estar_A_P, Estar_B, E_Q
        # E_A: +P→E_A_P(mixed✓), +Q→E_A_Q(mixed✓)
        # E_Q: +B→E_B_Q(mixed✓)
        # 3 dead-end forms → 2^3 = 8 variants
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[CX], B[N]
                products: P[C], Q[NX]
                enzymes: E, E_A[CX], E_Q[NX],
                    Estar[X], Estar_A_P[CX], Estar_B[NX]
            end
            steps: begin
                [E, A] ⇌ [E_A]
                [Estar, B] ⇌ [Estar_B]
                [E, Q] ⇌ [E_Q]
                [Estar, P] ⇌ [Estar_A_P]
                [E_A] <--> [Estar_A_P]
                [Estar_B] ⇌ [E_Q]
            end
        end
        spec = mechanism_spec_from_mechanism(
            m, bi_bi_pp_rxn)
        @test EnzymeMechanism(spec) === m
        result =
            EnzymeRates._expand_substrate_product_dead_ends(
                [spec], bi_bi_pp_rxn)
        @test length(result) == 8
    end
end
```

- [ ] **Step 2: Remove "All compile correctly" testset**

- [ ] **Step 3: Run tests, commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "Exact dead-end counts with @enzyme_mechanism definitions and derivation comments"
```

---

### Task 6: Rename testsets, strengthen RE→SS and constraint removal tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Rename all "Move N:" testsets**

Find and replace:
- "Move 1: RE→SS conversion" → "RE→SS conversion"
- "Move 2: Remove equivalence constraint" → "Remove equivalence constraint"
- "Move 3: Add dead-end regulator" → "Add dead-end regulator"
- "Move 4: Add allosteric regulator" → "Add allosteric regulator"
- "Move 5: Remove TR equivalence" → "Remove TR equivalence"
- "Move 6: Allosteric conversion (+1)" → "Allosteric conversion"

- [ ] **Step 2: Strengthen RE→SS tests**

Replace init_mechanisms indexing with `@enzyme_mechanism` definitions. Add:

- **All-SS with constrained dead-end RE binding**: Define a uni-uni mechanism where all catalytic steps are SS. Add a dead-end inhibitor binding to 2+ forms (so binding K's are constrained). Verify `_expand_re_to_ss` yields nothing.
- **Bi-bi with exact count**: Define bi-bi random with `@enzyme_mechanism`. Count eligible RE steps (not constrained). Verify exact result count.

- [ ] **Step 3: Strengthen constraint removal tests**

Replace init_mechanisms indexing with `@enzyme_mechanism` definitions. Add:

- **Multiple constraints with exact count**: Bi-bi random with K_A and K_B each constrained. Verify exact count of results and that each has one fewer constraint.
- **Substrate=regulator constraint independence**: Define mechanism where S is both substrate and dead-end inhibitor. Add S as regulator. Verify substrate K constraints and regulator K constraints are separate. Removing one type doesn't affect the other.

- [ ] **Step 4: Run tests, commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "Rename testsets, strengthen RE→SS and constraint removal tests"
```

---

### Task 7: Strengthen dead-end regulator tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Replace init_mechanisms indexing with @enzyme_mechanism**

All dead-end regulator tests use explicit mechanism definitions with round-trips.

- [ ] **Step 2: Add two-regulator test**

Define reaction with I and J. Add I first, then J. Verify that forms with both regulators bound exist (mirror steps where I and J both extend the same catalytic step endpoints).

- [ ] **Step 3: Fix bi-bi exact count**

Define bi-bi mechanism explicitly. Count eligible forms for regulator binding (forms where neither all subs nor all prods are bound). For bi-bi random with 7 forms:
```julia
# Eligible forms: E(0 bound), E_A(1 sub), E_B(1 sub),
#   E_P(1 prod), E_Q(1 prod) = 5 eligible
# (E_A_B has all subs → ineligible)
# (E_P_Q has all prods → ineligible)
# 2^5 - 1 = 31 non-empty subsets
```
Verify exact count = 31.

- [ ] **Step 4: Add allosteric test**

- Allosteric spec with allosteric-only regulator → `_expand_add_dead_end_regulator` yields nothing
- Allosteric spec where one regulator is dead-end-eligible → yields results

- [ ] **Step 5: Add ping-pong test with exact count**

Define ping-pong mechanism, count eligible forms, verify exact dead-end regulator count.

- [ ] **Step 6: Run tests, commit**

```bash
julia --project -e 'using Pkg; Pkg.test()'
git add test/test_mechanism_enumeration.jl
git commit -m "Strengthen dead-end regulator tests with exact counts and edge cases"
```

---

### Task 8: Final cleanup and verification

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Verify no init_mechanisms indexing remains**

Search for `first(filter(` and `specs[` patterns in the test file. Replace any remaining with `@enzyme_mechanism` definitions.

- [ ] **Step 2: Verify all @enzyme_mechanism have round-trips**

Search for `@enzyme_mechanism` — each must be followed by `mechanism_spec_from_mechanism` + `@test EnzymeMechanism(spec) ===`.

- [ ] **Step 3: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 4: Commit any remaining fixes**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Final test cleanup: all mechanisms explicitly defined with round-trips"
```

---

## Self-Review

**Spec coverage:**
- ✅ Self-contained file (Task 2-3)
- ✅ Allosteric dispatch + unified _add_expansions! (Task 1)
- ✅ Universal round-trips (Task 4)
- ✅ Ter-ter/ter-bi topologies (Task 4)
- ✅ Exact dead-end counts with derivations (Task 5)
- ✅ Renamed testsets (Task 6)
- ✅ Strengthened RE→SS, constraint, dead-end reg tests (Tasks 6-7)
- ✅ No init_mechanisms indexing (Task 8)
- ✅ Delete old files (Task 3)

**TDD followed:**
- Task 1: failing test → implement → verify ✅
- Tasks 2-3: copy then delete with test verification ✅
- Tasks 4-8: test strengthening with continuous verification ✅
