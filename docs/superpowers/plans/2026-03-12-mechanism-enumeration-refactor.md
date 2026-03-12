# Mechanism Enumeration Refactor Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the enumerate-all-forms search-based mechanism enumeration with a constructive builder using a simplified step-based representation.

**Architecture:** The `MechanismSpec` type changes from edge-index pairs with external form lookups to inline `StepSpec` structs (reactants, products, is_equilibrium). Each pipeline stage operates directly on step data — no adjacency recomputation. The constructive catalytic topology builder replaces DFS cycle-finding over pre-enumerated forms.

**Tech Stack:** Julia, EnzymeRates.jl internal types (`EnzymeReaction`, `EnzymeMechanism`, `AllostericEnzymeMechanism`)

**Spec:** `docs/superpowers/specs/2026-03-11-mechanism-enumeration-refactor-design.md`

**Test command:** `julia --project -e 'using Pkg; Pkg.test()'`

**Targeted test command:** `julia --project -e 'using Test, EnzymeRates; include("test/mechanism_enumeration_test_specs.jl"); include("test/test_mechanism_enumeration.jl")'`

---

## Chunk 1: Foundation + Stage 1 + compile_mechanism

### Task 1: Define new types and helpers

**Files:**
- Modify: `src/mechanism_enumeration.jl` (replace type definitions at top of file)

- [ ] **Step 1: Write StepSpec and new MechanismSpec types**

Replace `SiteDefinition`, `EnzymeFormSpec`, old `MechanismSpec` with:

```julia
"""Elementary step in canonical binding direction (metabolite on LHS)."""
struct StepSpec
    reactants::Vector{Symbol}   # [:E, :S] or [:EAB]
    products::Vector{Symbol}    # [:ES] or [:EPQ]
    is_equilibrium::Bool
end

struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{StepSpec}
    param_constraints::Vector{ParamConstraint}
    param_count::Int
end
```

Keep `ParamConstraint`, `AbstractMechanismSpec`, `AllostericMechanismSpec` (update `base` field type).

- [ ] **Step 2: Write helper functions for StepSpec**

```julia
"""Return the metabolite for a step, or nothing for isomerization."""
step_metabolite(s::StepSpec) =
    length(s.reactants) == 2 ? s.reactants[2] : nothing

"""Return (from_form, to_form) for a step."""
step_forms(s::StepSpec) = (s.reactants[1], s.products[1])

"""Collect all unique form names from a MechanismSpec."""
function all_form_names(spec::MechanismSpec)
    forms = Set{Symbol}()
    for s in spec.steps
        push!(forms, s.reactants[1])
        push!(forms, s.products[1])
    end
    forms
end

```

Also make `all_form_names` accept a `Vector{StepSpec}` (used by `_expand_ress_variants` which operates on steps before constructing a full `MechanismSpec`):
```julia
function all_form_names(steps::Vector{StepSpec})
    forms = Set{Symbol}()
    for s in steps
        push!(forms, s.reactants[1])
        push!(forms, s.products[1])
    end
    forms
end
```

Do NOT add `is_bound_at` — it's not needed by any task (YAGNI). Dead-end expansion will track bound metabolites via form names directly.

- [ ] **Step 3: Replace `MechanismSpec` in-place with the new definition**

Replace the old `MechanismSpec` struct with the new step-based version. This WILL break all existing code that references old fields (`edges`, `n_catalytic_edges`, `equilibrium_steps`). That's intentional — each subsequent task rewrites one consumer. The old types `SiteDefinition`, `EnzymeFormSpec` stay temporarily (they are only used by functions being deleted in Task 12). `EnumerationStage` subtypes and `MechanismIterator` can be deleted now (only consumed by `enumerate_mechanisms` which is rewritten in Task 9).

**Transition strategy:** Since old and new `MechanismSpec` can't coexist under the same name, we replace in-place and accept that tests/functions referencing old fields will be broken until their task rewrites them. The package will still load (broken functions aren't called at load time), but tests won't pass until all tasks complete.

- [ ] **Step 4: Verify package loads**

Run: `julia --project -e 'using EnzymeRates'`
Expected: loads without error (existing tests should still pass since old code is preserved)

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "Replace MechanismSpec with step-based representation"
```

---

### Task 2: Constructive catalytic topology builder

**Files:**
- Modify: `src/mechanism_enumeration.jl` (replace `_catalytic_topologies`, remove `enumerate_enzyme_forms`, `_classify_edge`, `_build_adjacency`, `_is_binding_direction`)

**Reference:** Spec §Stage 1 for algorithm table and ping-pong rules.

- [ ] **Step 1: Write a test for Uni-Uni catalytic topology**

The existing test at `test/test_mechanism_enumeration.jl:8-32` already covers this. Verify it will work with the new types by checking what it asserts:
- `length(topos) == 1`
- `compile_mechanism(topos[1]) === m_uu` (needs compile_mechanism — Task 3)
- `t.n_catalytic_edges == length(t.edges)` → update to `length(t.steps)`
- `count(.!t.equilibrium_steps) >= 1` → update to `count(!s.is_equilibrium for s in t.steps) >= 1`

Don't modify tests yet — Stage 1 tests will be updated in Task 3 after compile_mechanism is ready.

- [ ] **Step 2: Implement `_catalytic_topologies` with constructive builder**

The function signature stays the same: `_catalytic_topologies(reaction::EnzymeReaction) -> Vector{MechanismSpec}`

Implementation approach:
1. Extract substrate/product names and atom compositions from `reaction`
2. Define a recursive `backtrack!` function that tracks:
   - Current form name (Symbol)
   - Accumulated atoms on enzyme (Dict{Symbol,Int})
   - Set of bound substrates
   - Set of released products
   - Whether enzyme has a residual (atoms from prior ping-pong)
   - Path of steps taken so far
3. At each state, try all valid next steps per the spec's algorithm table
4. When a complete cycle (back to E, all substrates consumed, all products released) is found, record the path
5. After backtracking, take power-set combinations of paths
6. For each combination, collect unique steps, compute param_count, create `MechanismSpec`

Form naming convention: construct from bound metabolites, e.g., `:E`, `:E_S`, `:E_S_Q`, `:E_A_B`. Use a consistent naming function.

Key helpers needed:
- `_form_name(bound_subs, bound_prods, residual_atoms)` → Symbol
- `_atoms_dict(reaction, metabolite)` → Dict{Symbol,Int} for substrate/product atoms
- `_can_pingpong(accumulated_atoms, product_atoms)` → Bool (product atoms ⊆ accumulated)

- [ ] **Step 3: Verify `_catalytic_topologies` returns correct count for Uni-Uni**

Run: `julia --project -e '
using EnzymeRates
rxn = @enzyme_reaction substrates: S[C] products: P[C]
topos = EnzymeRates._catalytic_topologies(rxn)
println("Uni-Uni topologies: ", length(topos))
println("Steps per topo: ", [length(t.steps) for t in topos])
'`
Expected: 1 topology with 3 steps

- [ ] **Step 4: Verify Bi-Bi topology count**

Run same approach for Bi-Bi (A[C] + B[N] → P[C] + Q[N]):
Expected: 9 topologies (4 sequential + 4 mixed + 1 fully random)

- [ ] **Step 5: Verify Bi-Bi Ping-Pong topology count**

Run for Bi-Bi-PP (A[CX] + B[N] → P[C] + Q[NX]):
Expected: 10 topologies

- [ ] **Step 6: Verify Ter-Ter does NOT OOM**

Run for Ter-Ter (A[C] + B[N] + D[X] → P[C] + Q[N] + R[X]):
Expected: completes without OOM, returns topologies (count TBD — the constructive builder should handle 3+3)

- [ ] **Step 7: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "Implement constructive catalytic topology builder"
```

---

### Task 3: compile_mechanism + test infrastructure update

**Files:**
- Modify: `src/mechanism_enumeration.jl` (rewrite `compile_mechanism` for new `MechanismSpec`)
- Modify: `test/mechanism_enumeration_test_specs.jl` (update `mechanism_spec_from_mechanism`, update field references)
- Modify: `test/test_mechanism_enumeration.jl` (update field references for Stage 1 tests)

- [ ] **Step 1: Implement `compile_mechanism(spec::MechanismSpec)`**

Algorithm:
1. Collect form names from steps
2. Compute form atoms via BFS from free enzyme E:
   - Use `spec.reaction` to get substrate/product atom definitions
   - At each binding step, if metabolite is substrate → add atoms; if product → subtract atoms
   - Isomerization steps don't change atoms
   - Regulator binding (dummy names) don't change atoms
3. Build species tuple: `(substrates, products, regulators, enzyme_forms_with_atoms)`
4. Build reactions tuple from steps (convert Vector{Symbol} → Tuple{Symbol,...})
5. Strip `__regN` suffixes from metabolite names
6. Build equilibrium_steps and param_constraints tuples
7. Return `EnzymeMechanism(species, reactions, eq_steps, constraints)`

- [ ] **Step 2: Implement `compile_mechanism(spec::AllostericMechanismSpec)`**

This stays mostly the same — compile the base `MechanismSpec`, then wrap with allosteric info. Update to use `spec.base.steps` instead of `spec.base.edges`.

- [ ] **Step 3: Update `mechanism_spec_from_mechanism` in test specs**

This helper converts `EnzymeMechanism` → `MechanismSpec`. Rewrite to construct `StepSpec` entries from the mechanism's reactions:

```julia
function mechanism_spec_from_mechanism(m::EnzymeMechanism, rxn; kwargs...)
    steps = StepSpec[]
    for (i, rxn_step) in enumerate(reactions(m))
        push!(steps, StepSpec(
            collect(rxn_step[1]),
            collect(rxn_step[2]),
            equilibrium_steps(m)[i]))
    end
    # ... build MechanismSpec
end
```

- [ ] **Step 4: Update Stage 1 test assertions for new field names**

In `test/test_mechanism_enumeration.jl`, update references:
- `t.n_catalytic_edges == length(t.edges)` → `true` (at Stage 1, all steps are catalytic — this invariant is implicit)
- `count(.!t.equilibrium_steps) >= 1` → `count(!s.is_equilibrium for s in t.steps) >= 1`
- `t.edges` → `t.steps` where used for comparison
- `s.n_catalytic_edges` → remove or replace as appropriate

- [ ] **Step 5: Flip Ter-Ter @test_broken (line 268)**

Change `@test_broken false  # ter_ter OOMs` to an actual test:
```julia
@testset "Ter-Ter" begin
    topos = EnzymeRates._catalytic_topologies(ter_ter)
    @test length(topos) >= 1  # verify it completes without OOM; exact count TBD after implementation — harden to == once known
    for t in topos
        @test count(!s.is_equilibrium for s in t.steps) >= 1
    end
end
```

- [ ] **Step 6: Run Stage 1 tests**

Run: `julia --project -e '
using Test, EnzymeRates, Random
include("test/mechanism_enumeration_test_specs.jl")
include("test/test_mechanism_enumeration.jl")
'`
Expected: Stage 1 tests pass (other stages will fail — expected)

- [ ] **Step 7: Commit**

```bash
git add src/mechanism_enumeration.jl test/mechanism_enumeration_test_specs.jl test/test_mechanism_enumeration.jl
git commit -m "Add compile_mechanism for new types, update Stage 1 tests"
```

---

## Chunk 2: Stages 2-3

### Task 4: Stage 2 — RE/SS assignment

**Files:**
- Modify: `src/mechanism_enumeration.jl` (rewrite `_expand_ress_variants`)
- Modify: `test/test_mechanism_enumeration.jl` (flip @test_broken lines 281, 302, 316; update field references)

- [ ] **Step 1: Flip @test_broken for RE/SS tests**

Lines 281, 302, 316: change `@test_broken length(result) == 2^n_steps - 1` to `@test length(result) == 2^n_steps - 1`. Remove the `@test length(result) == 2^n_binding - 1` lines (those tested the buggy behavior).

Update field references in Stage 2 tests:
- `topo.equilibrium_steps` → `[s.is_equilibrium for s in topo.steps]`
- `s.equilibrium_steps` → `[s2.is_equilibrium for s2 in s.steps]`
- `s.edges` → `s.steps`
- `s.n_catalytic_edges` → remove

- [ ] **Step 2: Run tests to verify they fail**

Run targeted Stage 2 tests. Expected: FAIL (old implementation doesn't exist anymore or produces wrong counts)

- [ ] **Step 3: Implement `_compute_re_partition_from_steps` helper**

This helper is needed by `_expand_ress_variants` below. Adapt existing `_compute_re_partition` (union-find on form indices connected by RE edges) to work with `StepSpec` form names:

```julia
function _compute_re_partition_from_steps(steps::Vector{StepSpec})
    # Collect forms from RE steps only
    # Union-find: connect forms linked by RE steps
    # Return Vector{Vector{Symbol}} — groups of connected form names
end
```

- [ ] **Step 4: Implement `_expand_ress_variants`**

```julia
function _expand_ress_variants(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    max_re_groups::Int=7,
)
    result = MechanismSpec[]
    for spec in specs
        n = length(spec.steps)
        n == 0 && continue
        # Enumerate all 2^n - 1 RE/SS assignments (exclude all-RE)
        for mask in 1:(1 << n) - 1  # mask=0 is all-RE, skip
            steps = [StepSpec(s.reactants, s.products,
                        (mask >> (i-1)) & 1 == 0)
                     for (i, s) in enumerate(spec.steps)]
            any(!s.is_equilibrium for s in steps) || continue
            # Check RE-connected groups ≤ max_re_groups
            partition = _compute_re_partition_from_steps(steps)
            length(partition) > max_re_groups && continue
            # Compute param_count using the same formula as current code:
            # n_thermo = n_total_edges - n_forms + 1 (full-graph cycle rank)
            # This counts ALL independent cycles, not just RE cycles.
            # Thermodynamic constraints apply to every cycle regardless
            # of RE/SS status (Haldane/Wegscheider relations).
            n_re = count(s.is_equilibrium for s in steps)
            n_ss = n - n_re
            n_forms = length(all_form_names(steps))
            n_thermo = n - n_forms + 1
            pc = n_re + 2 * n_ss - n_thermo + 2
            push!(result, MechanismSpec(spec.reaction, steps,
                spec.param_constraints, pc))
        end
    end
    result
end
```

Needs helper: `_compute_re_partition_from_steps(steps)` — union-find on form names connected by RE steps. Adapt existing `_compute_re_partition` to work with `StepSpec` form names instead of integer indices.

- [ ] **Step 4: Run Stage 2 tests**

Expected: Stage 2 tests pass with correct counts (2^n_steps - 1)

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Implement RE/SS assignment with all-edge bitmask"
```

---

### Task 5: Stage 3a — Substrate/product dead-end expansion

**Files:**
- Modify: `src/mechanism_enumeration.jl` (add `_expand_substrate_product_dead_ends`)
- Modify: `test/test_mechanism_enumeration.jl` (flip @test_broken lines 349, 358, 367; write actual tests)

- [ ] **Step 1: Write tests for substrate/product dead-ends**

Replace `@test_broken false` at lines 349, 358, 367 with actual tests:

```julia
@testset "Uni-Uni: passthrough (no off-cycle forms)" begin
    topo = EnzymeRates._catalytic_topologies(uni_uni)[1]
    ress = EnzymeRates._expand_ress_variants([topo], uni_uni)
    result = EnzymeRates._expand_substrate_product_dead_ends(
        ress, uni_uni)
    @test length(result) == length(ress)  # no off-cycle forms
end

@testset "Bi-Bi: 4 off-cycle forms" begin
    topo = EnzymeRates._catalytic_topologies(bi_bi)
    # Find the fully random topology (most forms)
    random_topo = topo[findfirst(t ->
        length(all_form_names(t)) == 7, topo)]
    ress = EnzymeRates._expand_ress_variants(
        [random_topo], bi_bi)
    result = EnzymeRates._expand_substrate_product_dead_ends(
        [ress[1]], bi_bi)
    @test length(result) == 2^4  # 4 off-cycle forms
end
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL (`_expand_substrate_product_dead_ends` not defined)

- [ ] **Step 3: Implement `_expand_substrate_product_dead_ends`**

For each spec, for each catalytic form F where neither all substrates nor all products are bound:
1. Find available substrate/product metabolites not already bound at F
2. Filter: F + metabolite must not be a catalytic form, must satisfy validity constraint
3. For each power-set of available dead-end metabolites at F, create dead-end form + steps
4. Add binding steps (always RE) and mirror steps (inherit RE/SS)

Key helpers:
- `_bound_metabolites(spec, form)` → Set{Symbol} of metabolites bound at form
- `_is_eligible_for_dead_end(spec, form, reaction)` → Bool (neither all subs nor all prods)
- `_create_dead_end_form_name(base_form, added_metabolites)` → Symbol

- [ ] **Step 4: Run tests**

Expected: Stage 2.5 tests pass

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add substrate/product dead-end expansion"
```

---

### Task 6: Stage 3b — Regulator dead-end expansion

**Files:**
- Modify: `src/mechanism_enumeration.jl` (add `_expand_regulator_dead_ends`, `_expand_dead_end`)
- Modify: `test/test_mechanism_enumeration.jl` (flip @test_broken lines 391, 409, 427, 439; update field refs)

- [ ] **Step 1: Flip @test_broken for regulator dead-end tests**

Lines 391, 409, 427, 439: change `@test_broken length(result) == N` to `@test length(result) == N`. Remove the lines testing buggy counts.

Update field references:
- `length(s.edges) - s.n_catalytic_edges` → compute dead-end step count differently (count steps with dead-end forms)
- `s.n_catalytic_edges` → remove

Rename all test references from `_expand_dead_end_inhibitors` to `_expand_dead_end` (the function is being renamed to reflect it now handles substrate/product dead-ends too, not just inhibitors). **Important:** this changes test semantics — `_expand_dead_end` now also does substrate/product dead-end expansion (Task 5). Verify that existing test expected values still hold under the combined function (for reactions without off-cycle substrate/product forms, like Uni-Uni, the result is unchanged).

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL (eligibility not enforced, or function not implemented yet)

- [ ] **Step 3: Implement `_expand_regulator_dead_ends`**

Similar to substrate/product dead-ends but:
- Uses dummy metabolite names (`:R__reg1`)
- Only considers dead-end regulator metabolites from `dead_end_regs`
- Same eligibility check (neither all subs nor all prods bound)

- [ ] **Step 4: Implement `_expand_dead_end` combining 3a and 3b**

`_expand_dead_end` calls both `_expand_substrate_product_dead_ends` and `_expand_regulator_dead_ends` internally to determine *available* dead-end metabolites per form, then does a single combined power-set expansion over all of them (per the spec §3c). The two helper functions are used for computing available metabolites and creating individual dead-end steps — not applied sequentially as separate stages.

```julia
function _expand_dead_end(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    dead_end_regs::Vector{Symbol}=Symbol[],
)
    # For each spec, for each eligible catalytic form:
    #   1. Collect available substrate/product dead-end mets (from 3a logic)
    #   2. Collect available regulator dead-end mets (from 3b logic)
    #   3. Single power-set over ALL available dead-end mets
    #   4. For each combination, create dead-end forms + binding + mirror steps
end
```

- [ ] **Step 5: Run Stage 3 tests**

Expected: Regulator dead-end tests pass with correct counts (eligibility enforced)

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add regulator dead-end expansion with eligibility check"
```

---

## Chunk 3: Stages 4-5 + Pipeline

### Task 7: Stage 4 — Equivalence constraints

**Files:**
- Modify: `src/mechanism_enumeration.jl` (rewrite `_expand_equivalence_constraints`)
- Modify: `test/test_mechanism_enumeration.jl` (flip @test_broken lines 533, 540; update field refs)

- [ ] **Step 1: Flip @test_broken for equivalence tests**

Lines 533, 540: these are placeholder `@test_broken false` with no test logic. Replace them with actual test assertions:

- **Line 533** (substrate/regulator same metabolite): Create a reaction where a metabolite is both a substrate and a dead-end regulator. Verify equivalence grouping keeps them separate (dummy name `:X__reg1` prevents grouping with catalytic `:X`).
- **Line 540** (substrate/product dead-end equivalence): Create a Bi-Bi reaction, expand dead-ends, then expand equivalence. Verify that substrate dead-end binding edges are included in equivalence groups with their catalytic counterparts.

Both tests need the reaction definition + full pipeline through Stage 4 to be meaningful.

- [ ] **Step 2: Run tests to verify they fail**

- [ ] **Step 3: Implement `_expand_equivalence_constraints`**

With the new representation, equivalence grouping is simple:
1. For each step with a metabolite, record `(metabolite_name, is_equilibrium)` → step index
2. Group by `(metabolite, re_ss_status)` where metabolite is non-product
3. Filter groups to ≥2 steps with same RE/SS status
4. For each power-set of groups (each bit = constrain/unconstrain), generate variant
5. Build param_constraints and adjust param_count

Dummy regulatory names (`__reg1`) naturally separate catalytic from regulatory equivalence groups.

- [ ] **Step 4: Run tests**

Expected: Equivalence constraint tests pass

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Implement equivalence constraints with dummy-name grouping"
```

---

### Task 8: Stage 5 — Deduplication adaptation

**Files:**
- Modify: `src/mechanism_enumeration.jl` (adapt `_deduplicate`, `_concentration_fingerprint`, `_compute_re_partition` to new types)
- Modify: `test/test_mechanism_enumeration.jl` (update field refs in Stage 5 tests)

- [ ] **Step 1: Update Stage 5 test field references**

Update references to `edges`, `equilibrium_steps`, `n_catalytic_edges` in Stage 5 tests.

Also flip @test_broken at line 626 (Bi-Bi dedup undercounts) — this cascading fix from the Stage 2 isom toggle fix should now produce 5 distinct fingerprints:
```julia
@test length(deduped) == 5  # was @test_broken
```
Remove the buggy-baseline assertion at line 627 (`@test length(deduped) == 3`).

Also update the `MechanismSpec` constructor call at line 632-634 (`MechanismSpec(topo.reaction, topo.edges, topo.n_catalytic_edges, ...)`) to use new field names.

- [ ] **Step 2: Adapt `_concentration_fingerprint` for step-based input**

`_compute_re_partition_from_steps` was already defined in Task 4. Reuse it here.


The core algorithm (spanning arborescences) stays the same. Adapt the interface:
- Build form-to-index mapping from step form names
- Build edge metabolite info from `step_metabolite(s)`
- Binding direction is always LHS (canonical form) — no `_is_binding_direction` needed
- Feed indexed data into existing arborescence algorithm

- [ ] **Step 4: Adapt `_deduplicate`**

Same algorithm: fingerprint + constraint descriptor as key. Update to use step-based accessors.

- [ ] **Step 5: Run Stage 5 tests**

Expected: all pass, including the previously-broken dedup undercount

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Adapt deduplication to step-based representation"
```

---

### Task 9: Pipeline orchestration + allosteric adaptation

**Files:**
- Modify: `src/mechanism_enumeration.jl` (rewrite `enumerate_mechanisms`, adapt allosteric stages)
- Modify: `test/test_mechanism_enumeration.jl` (update cross-stage and allosteric test references)
- Modify: `test/test_enzyme_derivation.jl` (update `WithDeadEnd()` references)

- [ ] **Step 1: Adapt allosteric stages (6-8) for new MechanismSpec**

Stages 6-8 access `spec.base.edges`, `spec.base.equilibrium_steps`, `spec.base.n_catalytic_edges`. Update to use `spec.base.steps`:
- `_expand_allosteric`: references `spec.base` — update to new fields
- `_expand_tr_equivalence`: accesses edges and adjacency — update to use step metabolites
- `_deduplicate_allosteric`: uses `_allosteric_canonical_key` — update field references

- [ ] **Step 2: Rewrite `enumerate_mechanisms`**

Replace `EnumerationStage` dispatch with simple pipeline:
```julia
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    max_re_groups::Int=7,
    catalytic_n::Int=0,
)
    catalytic = _catalytic_topologies(reaction)
    # ... regulator partitioning loop
    # Stage 2: RE/SS
    # Stage 3: dead-end
    # Stage 4: equivalence
    # Stage 5: dedup
    # Stages 6-8: allosteric (if applicable)
    # Return Vector{AbstractMechanismSpec}
end
```

Remove the `stage=` kwarg entirely (delete `EnumerationStage` types). Only `test/test_enzyme_derivation.jl` uses it (lines 924, 984) — replace with direct pipeline calls.

- [ ] **Step 3: Update `test/test_enzyme_derivation.jl`**

Replace `enumerate_mechanisms(rxn; stage=EnzymeRates.WithDeadEnd(), ...)` at lines 924 and 984 with direct pipeline calls:
```julia
specs = EnzymeRates._expand_dead_end(
    EnzymeRates._expand_ress_variants(
        EnzymeRates._catalytic_topologies(rxn), rxn),
    rxn; dead_end_regs=...)
```

- [ ] **Step 4: Update `_run_full_pipeline_stages` in `test/mechanism_enumeration_test_specs.jl`**

This helper (line ~192-251) calls `_expand_dead_end_inhibitors` and references old fields (`s.edges`, `s.equilibrium_steps`, `s.n_catalytic_edges`). Update to use `_expand_dead_end` and new step-based fields. Also update the `_compute_re_partition` call at test line ~926 to use the new step-based API.

- [ ] **Step 5: Update compile_mechanism round-trip tests**

Test lines ~1003-1018 reference old fields (`s.n_catalytic_edges`, `s.edges`, `s.equilibrium_steps`). Update to use new step-based fields. Update the `mechanism_spec_from_mechanism` call to not pass `n_catalytic_edges`.

- [ ] **Step 6: Run allosteric and cross-stage tests**

Expected: Stages 6-8 tests pass, cross-stage property tests pass

- [ ] **Step 7: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl test/test_enzyme_derivation.jl test/mechanism_enumeration_test_specs.jl
git commit -m "Adapt allosteric stages and pipeline orchestration"
```

---

## Chunk 4: Bug Fixes + End-to-End + Cleanup

### Task 10: T/R mirror dedup fix

**Files:**
- Modify: `src/mechanism_enumeration.jl` (fix `_allosteric_canonical_key`)
- Modify: `test/test_mechanism_enumeration.jl` (flip @test_broken line 851)

- [ ] **Step 1: Flip @test_broken at line 851**

Change to: `@test length(deduped) < length(tr)`. Remove the buggy-baseline assertion at line 852 (`@test length(deduped) == length(tr)`).

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL (complementary TR sets not recognized as mirrors)

- [ ] **Step 3: Fix `_allosteric_canonical_key`**

The canonical key should map complementary TR-equiv sets to the same key. For metabolites {A, B, P, Q, R1, R2}, the TR-equiv set {A, B} and its complement {P, Q, R1, R2} should produce the same key.

Fix: use `min(sorted_tr_equiv, sorted_complement)` as the canonical form:

```julia
function _allosteric_canonical_key(spec::AllostericMechanismSpec)
    base_key = (spec.base.steps,
                spec.base.param_constraints)
    pairs = collect(zip(spec.allosteric_reg_sites,
                        spec.allosteric_multiplicities))
    sort!(pairs)
    sorted_sites = [p[1] for p in pairs]
    sorted_mults = [p[2] for p in pairs]

    # Canonical TR-equiv: min of set and complement
    tr = sort(spec.tr_equiv_metabolites)
    # Collect all metabolites with T-state params
    all_t_mets = _collect_t_state_metabolites(spec)
    complement = sort(setdiff(all_t_mets, tr))
    canonical_tr = min(tr, complement)

    (base_key, spec.catalytic_n, sorted_sites, sorted_mults,
     canonical_tr)
end
```

Need helper `_collect_t_state_metabolites(spec)` that returns all metabolites with T-state parameters (same logic as `_expand_tr_equivalence` uses to determine `t_mets`).

- [ ] **Step 4: Run test**

Expected: PASS — mirrors now detected and removed

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Fix T/R mirror dedup to recognize complementary sets"
```

---

### Task 11: param_count formula fix

**Files:**
- Modify: `src/mechanism_enumeration.jl` (fix param_count computation in `_expand_equivalence_constraints`)
- Modify: `test/test_mechanism_enumeration.jl` (flip @test_broken line 1157)

- [ ] **Step 1: Flip @test_broken at line 1157**

Change `@test_broken n_match == 126` to `@test n_match == 126`. Remove the buggy-baseline assertion at line 1156 (`@test n_match == 90`).

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL (36 of 126 constrained specs have wrong param_count)

- [ ] **Step 3: Diagnose the root cause**

The issue: when equivalence constraints make two binding constants equal (e.g., K_S1 = K_S2), this can make a Wegscheider thermodynamic constraint redundant (because the constraint was: K_S1 * K_X = K_S2 * K_Y, which with K_S1=K_S2 simplifies to K_X = K_Y, reducing n_thermo by 1).

The current formula `param_count = n_re + 2*n_ss - n_thermo + 2` computes `n_thermo` from the cycle structure alone, not accounting for constraint-induced reductions.

- [ ] **Step 4: Implement fix**

After computing base param_count, check if any equivalence constraint makes a thermodynamic constraint redundant:
- For each independent cycle, the Wegscheider constraint is a product relationship among the K/k values on the cycle
- If an equivalence constraint forces two params in the cycle to be equal, check if this makes the Wegscheider constraint redundant
- If so, decrement n_thermo by 1 for each redundant constraint

This requires cycle analysis. Use the existing cycle detection from the dedup infrastructure.

- [ ] **Step 5: Run test**

Expected: PASS — all 126 constrained specs match

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Fix param_count for equivalence-constraint/Wegscheider interaction"
```

---

### Task 12: End-to-end verification + dead code cleanup

**Files:**
- Modify: `src/mechanism_enumeration.jl` (delete dead code)
- Modify: `test/test_mechanism_enumeration.jl` (update end-to-end expected counts)
- Modify: `test/mechanism_enumeration_test_specs.jl` (cleanup)

- [ ] **Step 1: Run full end-to-end tests to see what changed**

Run: `julia --project -e 'using Pkg; Pkg.test()'`

Capture the actual mechanism counts for each end-to-end test. The counts will have changed due to:
- RE/SS toggle fix (more variants)
- Dead-end eligibility fix (fewer variants)
- Substrate/product dead-end expansion (more variants)
- param_count fix (changes dedup winners)

- [ ] **Step 2: Hand-verify new counts are correct**

For each end-to-end test (Uni-Uni, Uni-Bi, Bi-Bi, Bi-Bi-PP, Uni-Uni+reg, Uni-Bi+reg):
- Trace through the pipeline mentally or with debug prints
- Verify the new counts make biochemical sense

Also harden the Ter-Ter topology count: replace the weak `>= 1` assertion from Task 3 with an exact `==` count now that the constructive builder's output is known.
- Update test assertions with new correct counts

- [ ] **Step 3: Update end-to-end test assertions**

Replace old expected counts with verified new counts in the end-to-end testset (lines 1038-1088).

- [ ] **Step 4: Delete dead code**

Remove from `src/mechanism_enumeration.jl` (deferred from Task 1 — old types kept alive until all consumers were rewritten):
- `SiteDefinition`, `EnzymeFormSpec`, `EnumerationStage` subtypes, `MechanismIterator` (old types)
- `enumerate_enzyme_forms` (replaced by constructive builder)
- `_classify_edge` (no longer needed)
- `_build_adjacency` (no longer needed)
- `_is_binding_direction` (canonical form handles this)
- `_find_dead_end` (dead-end forms constructed directly)
- `_dead_end_catalytic_map` (counterparts known at construction time)
- `_propagate_de_eq_steps!` (RE/SS inherited directly)
- `_find_first_isomerization` (no longer pinning isom to SS)
- `_edge_site_role` (replaced by dummy names)
- Old `EnumerationStage` types and `MechanismIterator` if not already removed

KEPT (adapted, not deleted): `_compute_re_partition` (→ `_compute_re_partition_from_steps`), `_concentration_fingerprint`, `_spanning_arborescence_monomials`, `_constraint_descriptor`, `_set_partitions`, `_partition_mult_count`

- [ ] **Step 5: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL tests pass (including Aqua, JET, test_enzyme_derivation)

- [ ] **Step 6: Verify no @test_broken remains in mechanism enumeration tests**

Run: `grep -n '@test_broken' test/test_mechanism_enumeration.jl`
Expected: no matches (all bugs fixed, all features implemented)

- [ ] **Step 7: Final commit**

```bash
git add -A  # after git status review
git commit -m "Complete mechanism enumeration refactor: cleanup dead code, update counts"
```
