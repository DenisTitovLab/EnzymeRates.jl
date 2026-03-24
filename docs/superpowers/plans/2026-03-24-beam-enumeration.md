# Beam-Based Mechanism Enumeration Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace eager mechanism enumeration with a level-by-level beam search that expands mechanisms by param_count, enabling enumeration of large reactions (ter-ter + 3 regulators) without OOM.

**Architecture:** Three pure expansion functions (`expand_mechanisms_same_param_count`, `expand_mechanisms_by_one_param`, `expand_mechanisms_by_two_params`) produce mechanism candidates at +0/+1/+2 param_count deltas. Parameter counting and deduplication use runtime versions of existing `@generated` algorithms (`_runtime_param_count`, `_runtime_denominator_monomials`) — no reimplementation, no JIT per mechanism. A new `enumerate_mechanisms` orchestrates these in a level-by-level loop, materializing one level at a time. The old pipeline is preserved as `old_enumerate_mechanisms` for correctness verification.

**Tech Stack:** Julia, EnzymeRates.jl internal types (`MechanismSpec`, `StepSpec`, `AllostericMechanismSpec`, `ParamConstraint`)

**Spec:** `docs/superpowers/specs/2026-03-24-beam-enumeration-design.md`

---

## File Structure

| File | Responsibility |
|------|---------------|
| `src/beam_enumeration.jl` | All beam search logic: three expansion functions, `_deduplicate_specs`, dead-end binding helpers, and `enumerate_mechanisms` orchestrator. Calls runtime entry points for param counting and dedup fingerprints. |
| `src/thermodynamic_constr_for_rate_eq_derivation.jl` | Refactor: extract runtime entry points from `_enzyme_incidence_matrix`, `_thermodynamic_constraints`, `_dependent_param_exprs`. Existing type-dispatch code becomes thin wrappers. |
| `src/rate_eq_derivation.jl` | Refactor: extract runtime entry point from `_raw_symbolic_rate_polys`. Existing type-dispatch code becomes thin wrapper. Dedup fingerprints use this + kinetic symbol stripping. |
| `src/mechanism_enumeration.jl` | Existing file — rename `enumerate_mechanisms` → `old_enumerate_mechanisms`. No other changes. |
| `src/EnzymeRates.jl` | Add `include("beam_enumeration.jl")` after mechanism_enumeration.jl |
| `test/test_beam_enumeration.jl` | All beam enumeration tests |
| `test/old_test_mechanism_enumeration.jl` | Renamed from `test/test_mechanism_enumeration.jl` (function reference updated) |
| `test/runtests.jl` | Update includes for renamed + new test files |

---

## Chunk 1: Rename Old Pipeline and Scaffold New Files

### Task 1: Rename old pipeline function and test file

**Files:**
- Modify: `src/mechanism_enumeration.jl` — rename `enumerate_mechanisms` → `old_enumerate_mechanisms`
- Rename: `test/test_mechanism_enumeration.jl` → `test/old_test_mechanism_enumeration.jl`
- Modify: `test/old_test_mechanism_enumeration.jl` — update any references to `enumerate_mechanisms` → `old_enumerate_mechanisms`
- Modify: `test/runtests.jl` — update include path

- [ ] **Step 1: Rename `enumerate_mechanisms` to `old_enumerate_mechanisms`**

In `src/mechanism_enumeration.jl`, find the function definition at line 2166:
```julia
function enumerate_mechanisms(
```
Rename to:
```julia
function old_enumerate_mechanisms(
```

- [ ] **Step 2: Rename the test file**

```bash
git mv test/test_mechanism_enumeration.jl test/old_test_mechanism_enumeration.jl
```

- [ ] **Step 3: Update test file references**

In `test/old_test_mechanism_enumeration.jl`, find any calls to `enumerate_mechanisms` and rename to `old_enumerate_mechanisms`. Search for:
```julia
enumerate_mechanisms(
```
The test file uses `_run_full_pipeline_stages` (a local helper) and the end-to-end test calls `EnzymeRates.enumerate_mechanisms`. Update these to `EnzymeRates.old_enumerate_mechanisms`.

- [ ] **Step 4: Update runtests.jl**

In `test/runtests.jl`, change:
```julia
include("test_mechanism_enumeration.jl")
```
to:
```julia
include("old_test_mechanism_enumeration.jl")
```

- [ ] **Step 5: Run tests to verify rename didn't break anything**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass (the old tests still exercise the same code, just renamed)

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/old_test_mechanism_enumeration.jl test/runtests.jl
git commit -m "Rename enumerate_mechanisms → old_enumerate_mechanisms

Preserve the current eager pipeline and its tests for correctness
verification against the new beam-based implementation."
```

### Task 2: Scaffold beam_enumeration.jl and test file

**Files:**
- Create: `src/beam_enumeration.jl`
- Create: `test/test_beam_enumeration.jl`
- Modify: `src/EnzymeRates.jl` — add include
- Modify: `test/runtests.jl` — add include

- [ ] **Step 1: Refactor existing functions to have runtime entry points**

The key functions in `src/thermodynamic_constr_for_rate_eq_derivation.jl`
currently take `Type{<:EnzymeMechanism}` and call `m = M()` to extract data.
Refactor each to have a runtime core that takes data directly, then have the
type-dispatch version call the runtime core.

Functions to refactor:

**In `src/thermodynamic_constr_for_rate_eq_derivation.jl`** (for param counting):

1. `_enzyme_incidence_matrix(M)` → runtime core takes `(enz_names, steps_lhs, steps_rhs)`
2. `_thermodynamic_constraints(M)` → runtime core takes incidence matrix + stoich data
3. `_dependent_param_exprs(M)` → runtime core takes `(enz_names, steps_lhs, steps_rhs, eq_steps, constraints, sub_names, prod_names, met_names)` plus data for `_binding_K_symbols` and pivot priority

**In `src/rate_eq_derivation.jl`** (for dedup fingerprints):

4. `_raw_symbolic_rate_polys(M)` → runtime core takes `(enz_names, steps_lhs, steps_rhs, eq_steps, constraints, sub_names, prod_names, reg_names)`. Returns `(num::POLY, denom_terms)`. The `@generated` `rate_equation` body calls this at compile time; dedup calls it at runtime and strips kinetic symbols to get concentration fingerprints.

Helper functions already called by `_raw_symbolic_rate_polys` that take data
(not types): `_compute_re_groups`, `_compute_alpha`, `_split_reaction_side`,
`_ss_contrib`. These don't need refactoring — they already take runtime data.
The refactoring is at the top-level entry point that currently calls `m = M()`.

Pattern for each:
```julia
# Runtime core: takes data, no types
function _enzyme_incidence_matrix(
    enz_names::Vector{Symbol},
    steps_lhs::Vector,
    steps_rhs::Vector,
)
    # ... existing algorithm, unchanged ...
end

# Type-dispatch wrapper: extracts data, calls runtime core
function _enzyme_incidence_matrix(M::Type{<:EnzymeMechanism})
    m = M()
    enz_names = [e[1] for e in enzyme_forms(m)]
    rxns = reactions(m)
    _enzyme_incidence_matrix(
        enz_names,
        [r[1] for r in rxns],
        [r[2] for r in rxns])
end
```

The implementer should:
1. Read each existing function carefully
2. Identify every accessor call (`reactions(m)`, `equilibrium_steps(m)`, etc.)
3. Replace with data parameters in the runtime core
4. Keep the type-dispatch version as a thin wrapper
5. Run existing tests to verify no regression

- [ ] **Step 2: Create src/beam_enumeration.jl with ABOUTME and runtime wrappers**

```julia
# ABOUTME: Level-by-level beam search for mechanism enumeration.
# ABOUTME: Expands mechanisms by param_count using three expansion functions.

# ─── Runtime Parameter Counting ───────────────────────────────

"""
    _runtime_param_count(spec) → Int

Compute parameter count at runtime (<1ms, no JIT) by calling
the runtime version of _dependent_param_exprs with data
extracted from the MechanismSpec.
"""
function _runtime_param_count(spec::MechanismSpec)
    enz_names = collect(all_form_names(spec.steps))
    steps_lhs = [Tuple(s.reactants) for s in spec.steps]
    steps_rhs = [Tuple(s.products) for s in spec.steps]
    eq_steps = Tuple(s.is_equilibrium for s in spec.steps)

    rxn = spec.reaction
    sub_names = Symbol[s[1] for s in substrates(rxn)]
    prod_names = Symbol[p[1] for p in products(rxn)]
    reg_names = collect(regulators(rxn))
    met_names = unique(vcat(
        sub_names, prod_names, reg_names))

    _, indep = _dependent_param_exprs(
        enz_names, steps_lhs, steps_rhs, eq_steps,
        spec.param_constraints,
        sub_names, prod_names, met_names)
    # indep params + Keq + E_total
    length(indep) + 2
end

# ─── Runtime Deduplication Fingerprint ────────────────────────

"""
    _runtime_denominator_monomials(spec) → Set{MONO}

Compute concentration monomials at runtime by calling the
runtime version of _raw_symbolic_rate_polys (King-Altman/Cha)
and stripping kinetic symbols from denominator monomials.
"""
function _runtime_denominator_monomials(
    spec::MechanismSpec,
)
    enz_names = collect(all_form_names(spec.steps))
    steps_lhs = [Tuple(s.reactants) for s in spec.steps]
    steps_rhs = [Tuple(s.products) for s in spec.steps]
    eq_steps = Tuple(s.is_equilibrium for s in spec.steps)
    rxn = spec.reaction
    sub_names = Symbol[s[1] for s in substrates(rxn)]
    prod_names = Symbol[p[1] for p in products(rxn)]
    reg_names = collect(regulators(rxn))

    _, denom_terms = _raw_symbolic_rate_polys(
        enz_names, steps_lhs, steps_rhs, eq_steps,
        spec.param_constraints,
        sub_names, prod_names, reg_names)

    # Expand denom POLY and strip kinetic symbols
    _strip_to_concentration_fingerprint(denom_terms)
end

"""
    _is_kinetic_symbol(s::Symbol) → Bool

True for rate constant symbols: K1, k2f, k3r, etc.
"""
function _is_kinetic_symbol(s::Symbol)
    occursin(r"^[kK]\d+[fr]?$", string(s))
end

"""
    _strip_to_concentration_fingerprint(denom_terms) → Set{MONO}

Extract concentration monomials from denominator polynomial
terms by stripping kinetic symbols (K's and k's).
"""
function _strip_to_concentration_fingerprint(
    denom_terms,
)
    fingerprint = Set{MONO}()
    for dt in denom_terms
        # dt has sigma (POLY) and cofactor (POLY)
        # Expand sigma * cofactor
        expanded = poly_mul(
            _expand_factored(dt.sigma), dt.cofactor)
        for mono in keys(expanded)
            stripped = MONO(
                filter(p -> !_is_kinetic_symbol(p.first),
                    mono))
            push!(fingerprint, stripped)
        end
    end
    fingerprint
end
```

- [ ] **Step 2: Create test/test_beam_enumeration.jl with ABOUTME and reaction definitions**

```julia
# ABOUTME: Tests for beam-based mechanism enumeration.
# ABOUTME: Tests individual expansion functions and level-by-level orchestration.

@testset "Beam Mechanism Enumeration" begin

end  # top-level testset
```

Copy the reaction definitions from `test/old_test_mechanism_enumeration.jl` lines 40-156 (the `const uni_uni = ...` through `const bi_bi_ping_pong_reg_unknown = ...` block) into the test file before the `@testset` block — these are the shared test reactions.

- [ ] **Step 3: Add include to EnzymeRates.jl**

In `src/EnzymeRates.jl`, after line 35 (`include("mechanism_enumeration.jl")`), add:
```julia
include("beam_enumeration.jl")
```

- [ ] **Step 4: Add include to runtests.jl**

In `test/runtests.jl`, after the `old_test_mechanism_enumeration.jl` include, add:
```julia
include("test_beam_enumeration.jl")
```

- [ ] **Step 5: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass (new test file has an empty testset)

- [ ] **Step 6: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl src/EnzymeRates.jl test/runtests.jl
git commit -m "Scaffold beam_enumeration.jl with runtime functions"
```

### Task 3: Test runtime functions

**Files:**
- Modify: `test/test_beam_enumeration.jl`

Validates that runtime functions produce correct results by comparing against `@generated` versions on compiled mechanisms.

- [ ] **Step 1: Write runtime function tests**

Inside the `@testset "Beam Mechanism Enumeration"` block:

```julia
@testset "Runtime functions" begin
    @testset "_runtime_param_count matches @generated" begin
        # Verify runtime param count matches compile_mechanism
        # + parameters() for all catalytic topologies
        for (name, rxn) in [("uni-uni", uni_uni),
                            ("uni-bi", uni_bi),
                            ("bi-bi", bi_bi),
                            ("bi-bi pp", bi_bi_ping_pong)]
            topos = EnzymeRates._catalytic_topologies(rxn)
            for spec in topos
                runtime_pc = EnzymeRates._runtime_param_count(
                    spec)
                m = compile_mechanism(spec)
                generated_pc = length(parameters(m))
                @test runtime_pc == generated_pc
            end
        end
    end

    @testset "_runtime_denominator_monomials" begin
        # Verify fingerprints are non-empty and consistent
        topos = EnzymeRates._catalytic_topologies(bi_bi)
        for spec in topos
            monos = EnzymeRates._runtime_denominator_monomials(
                spec)
            @test !isempty(monos)
        end

        # Two identical mechanisms should have same fingerprint
        spec1 = topos[1]
        fp1 = EnzymeRates._runtime_denominator_monomials(spec1)
        fp2 = EnzymeRates._runtime_denominator_monomials(spec1)
        @test fp1 == fp2

        # Verify runtime fingerprint matches old
        # _concentration_fingerprint for all topologies
        for spec in topos
            runtime_fp =
                EnzymeRates._runtime_denominator_monomials(
                    spec)
            partition =
                EnzymeRates._compute_re_partition_from_steps(
                    spec.steps)
            old_fp =
                EnzymeRates._concentration_fingerprint(
                    spec.steps, partition)
            @test runtime_fp == old_fp
        end
    end
end
```

- [ ] **Step 2: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All runtime function tests pass

- [ ] **Step 3: Commit**

```bash
git add test/test_beam_enumeration.jl
git commit -m "Add runtime function tests: param_count and denominator monomials"
```

---

## Chunk 2: expand_mechanisms_by_one_param — RE→SS Move

### Task 4: RE→SS expansion — test and implement

**Files:**
- Modify: `src/beam_enumeration.jl`
- Modify: `test/test_beam_enumeration.jl`

The simplest +1 move: flip one RE step to SS.

- [ ] **Step 1: Write failing test for RE→SS**

```julia
@testset "expand_mechanisms_by_one_param" begin

@testset "RE→SS move" begin
    # Uni-uni has 1 topology: 3 steps, 2 RE + 1 SS
    topos = EnzymeRates._catalytic_topologies(uni_uni)
    spec = topos[1]
    n_re = count(s -> s.is_equilibrium, spec.steps)
    @test n_re == 2  # sanity: 2 RE steps to flip

    results = EnzymeRates.expand_mechanisms_by_one_param(
        [spec], uni_uni)
    re_to_ss = filter(
        r -> r.param_count == spec.param_count + 1,
        results)

    # Should produce exactly 2 candidates (one per RE step)
    @test length(re_to_ss) == n_re

    # Each candidate should have one fewer RE step
    for r in re_to_ss
        r_n_re = count(s -> s.is_equilibrium, r.steps)
        @test r_n_re == n_re - 1
        @test r.param_count == spec.param_count + 1
    end

    # Bi-bi topo 1: 5 steps, 4 RE + 1 SS → 4 candidates
    topos_bb = EnzymeRates._catalytic_topologies(bi_bi)
    spec_bb = topos_bb[1]
    n_re_bb = count(s -> s.is_equilibrium, spec_bb.steps)
    @test n_re_bb == 4

    results_bb = EnzymeRates.expand_mechanisms_by_one_param(
        [spec_bb], bi_bi)
    re_to_ss_bb = filter(
        r -> r.param_count == spec_bb.param_count + 1,
        results_bb)
    @test length(re_to_ss_bb) == n_re_bb

    # Verify param_count matches compiled mechanism
    for r in re_to_ss_bb
        m = compile_mechanism(r)
        @test length(parameters(m)) == r.param_count
    end
end

end  # expand_mechanisms_by_one_param testset
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `expand_mechanisms_by_one_param` not defined

- [ ] **Step 3: Implement expand_mechanisms_by_one_param with RE→SS only**

In `src/beam_enumeration.jl`:

```julia
# ─── Expansion: +1 Parameter ─────────────────────────────────

"""
    expand_mechanisms_by_one_param(specs, reaction) → Vector{MechanismSpec}

Generate mechanism candidates with param_count + 1 by applying:
1. RE→SS conversion (one step at a time)
2. Remove one equivalence constraint
3. Add dead-end binding (+1 configurations)
4. Remove TR equivalence (allosteric only)
"""
function expand_mechanisms_by_one_param(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = MechanismSpec[]
    for spec in specs
        _expand_re_to_ss!(result, spec, reaction)
    end
    result
end

"""Convert each RE step to SS, producing one candidate per RE step."""
function _expand_re_to_ss!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    for (i, step) in enumerate(spec.steps)
        step.is_equilibrium || continue
        new_steps = copy(spec.steps)
        new_steps[i] = StepSpec(
            step.reactants, step.products, false)
        candidate = MechanismSpec(
            spec.reaction, new_steps,
            copy(spec.param_constraints), 0)
        candidate = MechanismSpec(
            spec.reaction, new_steps,
            copy(spec.param_constraints),
            _runtime_param_count(candidate))
        push!(result, candidate)
    end
end
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: RE→SS tests pass

- [ ] **Step 5: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl
git commit -m "Add RE→SS move to expand_mechanisms_by_one_param"
```

---

## Chunk 3: expand_mechanisms_by_one_param — Remove Equivalence Constraint Move

### Task 5: Remove equivalence constraint — test and implement

**Files:**
- Modify: `src/beam_enumeration.jl`
- Modify: `test/test_beam_enumeration.jl`

To test this, we need a mechanism WITH constraints. We'll manually construct one by running the old pipeline's dead-end expansion on a uni-uni + I reaction to get a mechanism with equivalence constraints, or construct one by hand.

- [ ] **Step 1: Write failing test**

```julia
@testset "Remove equivalence constraint" begin
    # Build a mechanism with a known constraint.
    # Uni-uni + dead-end I: E + S ⇌ ES → EP ⇌ E + P
    # Adding I to both E and ES with K_I shared creates
    # a constraint. We construct this manually.
    topos = EnzymeRates._catalytic_topologies(uni_uni)
    base = topos[1]

    # Add dead-end I binding to E and ES with shared K
    # Steps: base steps + [E, I] ⇌ [EI] + [ES, I] ⇌ [ESI]
    # Mirror step: [EI, S] ⇌ [ESI] (inherits RE from E+S)
    de_steps = copy(base.steps)
    push!(de_steps, StepSpec([:E, :I], [:E_I], true))
    push!(de_steps, StepSpec([:E_S, :I], [:E_I_S], true))
    push!(de_steps, StepSpec(
        [:E_I, :S], [:E_I_S], true))

    # Constraint: K for step 6 (ES+I) = K for step 5 (E+I)
    # i.e., K_I is the same regardless of S binding
    constraint = (Symbol("K6"), 1,
        [(Symbol("K5"), 1)])
    constraints = [constraint]

    spec_with_constraint = EnzymeRates.MechanismSpec(
        uni_uni_dead_end_I, de_steps, constraints, 0)
    pc = EnzymeRates._runtime_param_count(
        spec_with_constraint)
    spec_with_constraint = EnzymeRates.MechanismSpec(
        uni_uni_dead_end_I, de_steps, constraints, pc)

    results = EnzymeRates.expand_mechanisms_by_one_param(
        [spec_with_constraint], uni_uni_dead_end_I)

    # Filter to only "remove constraint" candidates
    # (exclude RE→SS candidates)
    remove_constraint = filter(results) do r
        length(r.param_constraints) <
            length(spec_with_constraint.param_constraints)
    end

    # Should produce exactly 1 candidate (drop the one
    # constraint)
    @test length(remove_constraint) == 1
    @test isempty(
        remove_constraint[1].param_constraints)
    @test remove_constraint[1].param_count ==
        spec_with_constraint.param_count + 1
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — no "remove constraint" candidates produced

- [ ] **Step 3: Implement _expand_remove_constraint!**

Add to `src/beam_enumeration.jl`:

```julia
"""Remove one equivalence constraint, producing one candidate per constraint."""
function _expand_remove_constraint!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    for (i, _) in enumerate(spec.param_constraints)
        new_constraints = [
            spec.param_constraints[j]
            for j in eachindex(spec.param_constraints)
            if j != i]
        candidate = MechanismSpec(
            spec.reaction, copy(spec.steps),
            new_constraints, 0)
        push!(result, MechanismSpec(
            spec.reaction, copy(spec.steps),
            new_constraints,
            _runtime_param_count(candidate)))
    end
end
```

And add the call in `expand_mechanisms_by_one_param`:
```julia
function expand_mechanisms_by_one_param(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = MechanismSpec[]
    for spec in specs
        _expand_re_to_ss!(result, spec, reaction)
        _expand_remove_constraint!(result, spec, reaction)
    end
    result
end
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl
git commit -m "Add remove-equivalence-constraint move"
```

---

## Chunk 4: Dead-End Binding Helpers and +1 Dead-End Move

### Task 6: Dead-end binding helpers — bound metabolite tracking and form naming

**Files:**
- Modify: `src/beam_enumeration.jl`
- Modify: `test/test_beam_enumeration.jl`

These helpers determine what metabolites are bound at each enzyme form and compute dead-end form names. Written from scratch (not reusing old pipeline code).

- [ ] **Step 1: Write failing tests for helpers**

```julia
@testset "Dead-end binding helpers" begin
    # _bound_entities_at_form: track what's bound at each form
    topos = EnzymeRates._catalytic_topologies(bi_bi)
    spec = topos[1]  # ordered bi-bi: E→EA→EAB→EPQ→EQ→E

    bound = EnzymeRates._bound_entities(spec)
    @test bound[:E] == Set{Symbol}()
    @test bound[:E_A] == Set([:A])
    @test bound[:E_A_B] == Set([:A, :B])
    @test bound[:E_P_Q] == Set([:P, :Q])
    @test bound[:E_Q] == Set([:Q])

    # _binding_capacity: max(n_subs, n_prods)
    @test EnzymeRates._binding_capacity(bi_bi) == 2
    @test EnzymeRates._binding_capacity(ter_ter) == 3
    @test EnzymeRates._binding_capacity(uni_uni) == 1
    @test EnzymeRates._binding_capacity(uni_bi) == 2
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — functions not defined

- [ ] **Step 3: Implement helpers**

In `src/beam_enumeration.jl`:

```julia
# ─── Dead-End Binding Helpers ─────────────────────────────────

"""
    _bound_entities(spec) → Dict{Symbol, Set{Symbol}}

For each enzyme form in the mechanism, return the set of metabolites/regulators
bound to it. Derived from step reactants/products: if a step binds metabolite M
to form F producing form FM, then FM has M plus everything F had.
"""
function _bound_entities(spec::MechanismSpec)
    bound = Dict{Symbol, Set{Symbol}}()
    # Initialize all forms with empty sets
    for s in spec.steps
        for f in (s.reactants[1], s.products[1])
            haskey(bound, f) || (bound[f] = Set{Symbol}())
        end
    end

    # Iterative propagation: repeat until stable
    changed = true
    while changed
        changed = false
        for s in spec.steps
            from_form = s.reactants[1]
            to_form = s.products[1]
            met = step_metabolite(s)

            # Forward: binding step adds metabolite
            expected_to = if met !== nothing
                union(bound[from_form], Set([met]))
            else
                copy(bound[from_form])
            end
            if !issubset(expected_to, bound[to_form])
                union!(bound[to_form], expected_to)
                changed = true
            end

            # Reverse: unbinding step
            expected_from = if met !== nothing
                setdiff(bound[to_form], Set([met]))
            else
                copy(bound[to_form])
            end
            if !issubset(expected_from, bound[from_form])
                union!(bound[from_form], expected_from)
                changed = true
            end
        end
    end
    bound
end

"""
    _binding_capacity(reaction) → Int

Maximum number of entities that can bind at the catalytic site.
Equal to max(n_substrates, n_products).
"""
function _binding_capacity(
    @nospecialize(reaction::EnzymeReaction),
)
    max(length(substrates(reaction)),
        length(products(reaction)))
end
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl
git commit -m "Add dead-end binding helpers: _bound_entities, _binding_capacity"
```

### Task 7: Add dead-end binding (+1) move — test and implement

**Files:**
- Modify: `src/beam_enumeration.jl`
- Modify: `test/test_beam_enumeration.jl`

This is the most complex +1 move. For each unbound metabolite/regulator at each eligible form, generate dead-end binding configurations with maximal equivalence constraints, keeping only those at +1 param.

- [ ] **Step 1: Write failing test**

```julia
@testset "Add dead-end binding (+1)" begin
    # Uni-uni + dead-end I: E ⇌ ES → EP ⇌ E + P
    # Regulator I can bind to E and/or ES (but not both
    # of S,P are bound, since capacity=1 for uni-uni)
    topos = EnzymeRates._catalytic_topologies(uni_uni)
    spec = topos[1]

    results = EnzymeRates.expand_mechanisms_by_one_param(
        [spec], uni_uni_dead_end_I)
    de_results = filter(results) do r
        # Dead-end candidates have more steps than source
        length(r.steps) > length(spec.steps)
    end

    # Each candidate should have param_count = spec + 1
    for r in de_results
        @test r.param_count == spec.param_count + 1
    end

    # At least one candidate should exist (I binding to
    # some form)
    @test !isempty(de_results)

    # Verify param_count matches compiled mechanism
    for r in de_results
        m = compile_mechanism(r)
        @test length(parameters(m)) == r.param_count
    end

    # Bi-bi + dead-end I: more forms, more opportunities
    topos_bb = EnzymeRates._catalytic_topologies(bi_bi)
    spec_bb = topos_bb[1]  # ordered bi-bi

    results_bb = EnzymeRates.expand_mechanisms_by_one_param(
        [spec_bb], bi_bi_dead_end_I)
    de_results_bb = filter(results_bb) do r
        length(r.steps) > length(spec_bb.steps)
    end

    for r in de_results_bb
        @test r.param_count == spec_bb.param_count + 1
        m = compile_mechanism(r)
        @test length(parameters(m)) == r.param_count
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — no dead-end candidates produced

- [ ] **Step 3: Implement _expand_add_dead_end!**

This is the most substantial implementation. Add to `src/beam_enumeration.jl`:

```julia
"""
    _dead_end_form_name(bound_set, new_met) → Symbol

Generate a canonical form name for a dead-end complex.
Sorts all bound metabolites alphabetically.
"""
function _beam_dead_end_form_name(
    bound::Set{Symbol}, new_met::Symbol,
)
    all_bound = sort!(collect(union(bound, Set([new_met]))))
    Symbol("E_" * join(all_bound, "_"))
end

"""
    _dead_end_opportunities(spec, reaction) →
        Vector{Tuple{Symbol, Symbol}}

Find (form, metabolite) pairs where dead-end binding is possible.
A binding is valid when:
- The metabolite is not already bound at the form
- The result doesn't exceed binding capacity
- The result is not a catalytic form (would create a shortcut)
- For substrate/product dead-ends: must bind at least one sub
  AND one prod (mixed binding)
"""
function _dead_end_opportunities(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    bound = _bound_entities(spec)
    cap = _binding_capacity(reaction)
    cat_forms = all_form_names(spec)
    all_forms = collect(keys(bound))

    sub_names = Set(s[1] for s in substrates(reaction))
    prod_names = Set(p[1] for p in products(reaction))
    all_mets = union(sub_names, prod_names)

    # Dead-end regulators from reaction (includes unknown-role)
    roles = regulator_roles(reaction)
    de_regs = Symbol[
        r[1] for r in roles
        if r[2] == :dead_end || r[2] == :unknown]

    opportunities = Tuple{Symbol, Symbol}[]

    for form in sort!(collect(keys(bound)))
        fb = bound[form]
        # Count entities toward capacity (subs + prods +
        # dead-end regs)
        n_bound = length(fb)
        n_bound >= cap && continue

        # Substrate/product dead-end opportunities
        for met in sort!(collect(all_mets))
            met in fb && continue
            # Result must not be a catalytic form
            de_name = _beam_dead_end_form_name(fb, met)
            de_name in cat_forms && continue
            push!(opportunities, (form, met))
        end

        # Dead-end regulator opportunities
        for reg in de_regs
            reg in fb && continue
            push!(opportunities, (form, reg))
        end
    end
    opportunities
end

"""
    _expand_add_dead_end!(result, spec, reaction)

Add dead-end binding configurations at exactly +1 param.
Uses shared helper `_expand_dead_end_at_delta!`.
"""
function _expand_add_dead_end!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    _expand_dead_end_at_delta!(result, spec, reaction, 1)
end

"""
    _expand_dead_end_at_delta!(result, spec, reaction, delta)

Shared dead-end expansion logic. For each metabolite/regulator,
tries binding to 1, 2, ... forms with maximal equivalence
constraints. Keeps only configurations where computed param_count
equals spec.param_count + delta.
"""
function _expand_dead_end_at_delta!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
    delta::Int,
)
    opps = _dead_end_opportunities(spec, reaction)
    isempty(opps) && return

    bound = _bound_entities(spec)
    cat_forms = all_form_names(spec)

    by_met = Dict{Symbol, Vector{Symbol}}()
    for (form, met) in opps
        push!(get!(by_met, met, Symbol[]), form)
    end

    target_pc = spec.param_count + delta

    for (met, forms) in by_met
        n = length(forms)
        for mask in 1:(1 << n) - 1
            selected = Symbol[
                forms[i] for i in 1:n
                if (mask >> (i - 1)) & 1 == 1]

            new_steps = copy(spec.steps)
            new_constraints = copy(spec.param_constraints)

            first_step_idx = nothing
            for (j, form) in enumerate(selected)
                de_name = _beam_dead_end_form_name(
                    bound[form], met)
                push!(new_steps, StepSpec(
                    [form, met], [de_name], true))
                step_idx = length(new_steps)

                if j == 1
                    first_step_idx = step_idx
                else
                    push!(new_constraints, (
                        Symbol("K$step_idx"), 1,
                        [(Symbol("K$first_step_idx"), 1)]))
                end
            end

            # Add catalytic equivalence for +0 delta
            if delta == 0
                _add_catalytic_equivalence!(
                    new_steps, new_constraints,
                    spec.steps, selected, bound, met)
            end

            _add_mirror_steps!(
                new_steps, new_constraints, spec.steps,
                selected, bound, met, cat_forms)

            candidate = MechanismSpec(
                spec.reaction, new_steps,
                new_constraints, 0)
            pc = _runtime_param_count(candidate)
            pc == target_pc || continue

            push!(result, MechanismSpec(
                spec.reaction, new_steps,
                new_constraints, pc))
        end
    end
end

"""Add mirror steps for dead-end binding.

For each existing step [F1, M] ⇌ [F2], if both F1 and F2 have
dead-end extensions with the same metabolite, add a mirror step
[F1_de, M] ⇌ [F2_de] inheriting RE/SS from the original.
"""
function _add_mirror_steps!(
    new_steps, new_constraints, orig_steps,
    selected_forms, bound, de_met, cat_forms,
)
    for s in orig_steps
        from = s.reactants[1]
        to = s.products[1]
        met = step_metabolite(s)

        from in selected_forms || continue
        to in selected_forms || continue
        de_met in bound[from] && continue
        de_met in bound[to] && continue

        from_de = _beam_dead_end_form_name(
            bound[from], de_met)
        to_de = _beam_dead_end_form_name(
            bound[to], de_met)

        if met !== nothing
            push!(new_steps, StepSpec(
                [from_de, met], [to_de],
                s.is_equilibrium))
        else
            push!(new_steps, StepSpec(
                [from_de], [to_de],
                s.is_equilibrium))
        end

        # Add equivalence constraint: mirror step params
        # match catalytic step params
        mirror_idx = length(new_steps)
        orig_idx = findfirst(
            x -> x == s, orig_steps)
        if orig_idx !== nothing
            if s.is_equilibrium
                push!(new_constraints, (
                    Symbol("K$mirror_idx"), 1,
                    [(Symbol("K$orig_idx"), 1)]))
            else
                # SS: constrain both kf and kr
                push!(new_constraints, (
                    Symbol("k$(mirror_idx)f"), 1,
                    [(Symbol("k$(orig_idx)f"), 1)]))
                push!(new_constraints, (
                    Symbol("k$(mirror_idx)r"), 1,
                    [(Symbol("k$(orig_idx)r"), 1)]))
            end
        end
    end
end
```

Update `expand_mechanisms_by_one_param` to call the new function:
```julia
function expand_mechanisms_by_one_param(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = MechanismSpec[]
    for spec in specs
        _expand_re_to_ss!(result, spec, reaction)
        _expand_remove_constraint!(result, spec, reaction)
        _expand_add_dead_end!(result, spec, reaction)
    end
    result
end
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl
git commit -m "Add dead-end binding (+1) move with maximal constraints"
```

### Task 8: Test binding capacity limit and multi-level dead-end

**Files:**
- Modify: `test/test_beam_enumeration.jl`

- [ ] **Step 1: Write binding capacity and multi-level tests**

```julia
@testset "Binding capacity limit" begin
    # Uni-uni: capacity = 1
    # Dead-end forms can have at most 1 bound entity
    # E_S already has S bound (capacity 1), so I cannot
    # bind to E_S for uni-uni (would make 2 bound entities)
    topos = EnzymeRates._catalytic_topologies(uni_uni)
    spec = topos[1]
    results = EnzymeRates.expand_mechanisms_by_one_param(
        [spec], uni_uni_dead_end_I)
    de_results = filter(r -> length(r.steps) > length(spec.steps), results)

    for r in de_results
        bound = EnzymeRates._bound_entities(r)
        cap = EnzymeRates._binding_capacity(uni_uni_dead_end_I)
        for (form, entities) in bound
            @test length(entities) <= cap
        end
    end
end

@testset "Multi-level dead-end binding" begin
    # Bi-bi: capacity = 2
    # After adding I to E (creating EI), S can bind to EI
    # (creating EIS) since capacity=2 allows 2 entities.
    # This tests that dead-end forms are valid binding
    # targets.
    topos = EnzymeRates._catalytic_topologies(bi_bi)
    spec = topos[1]

    # First expansion: add I to get mechanisms with EI
    results1 = EnzymeRates.expand_mechanisms_by_one_param(
        [spec], bi_bi_dead_end_I)
    de_with_I = filter(results1) do r
        length(r.steps) > length(spec.steps)
    end
    @test !isempty(de_with_I)

    # Second expansion: expand the I-containing mechanisms
    # Some should have substrate binding to dead-end forms
    results2 = EnzymeRates.expand_mechanisms_by_one_param(
        de_with_I, bi_bi_dead_end_I)
    multi_level = filter(results2) do r
        # Has more steps than any single-level dead-end
        length(r.steps) > maximum(
            length(d.steps) for d in de_with_I)
    end
    # Multi-level dead-end forms should exist for bi-bi
    # (capacity allows 2 bound entities)
    @test !isempty(multi_level)

    # All forms respect capacity
    cap = EnzymeRates._binding_capacity(bi_bi_dead_end_I)
    for r in multi_level
        bound = EnzymeRates._bound_entities(r)
        for (_, entities) in bound
            @test length(entities) <= cap
        end
    end
end
```

- [ ] **Step 2: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass (if dead-end implementation considers dead-end forms as targets). If tests fail, fix the implementation to include dead-end forms in `_dead_end_opportunities`.

- [ ] **Step 3: Commit**

```bash
git add test/test_beam_enumeration.jl
git commit -m "Add tests for binding capacity limit and multi-level dead-end"
```

---

## Chunk 5: Deduplication and +0 Expansion

### Task 9: Deduplication within levels — test and implement

**Files:**
- Modify: `src/beam_enumeration.jl`
- Modify: `test/test_beam_enumeration.jl`

Reimplement deduplication from scratch using concentration fingerprints. The fingerprint computation uses `_concentration_fingerprint` and `_compute_re_partition_from_steps` which are already in `mechanism_enumeration.jl` — these are shared graph-theoretic utilities, not pipeline logic.

- [ ] **Step 1: Write failing test**

```julia
@testset "Deduplication within levels" begin
    # Two different catalytic topologies for bi-bi that
    # produce equivalent rate equations should deduplicate
    topos = EnzymeRates._catalytic_topologies(bi_bi)
    @test length(topos) >= 2

    # Expand all by one param (RE→SS)
    all_expanded = EnzymeRates.expand_mechanisms_by_one_param(
        topos, bi_bi)

    # Deduplication should reduce count
    deduped = EnzymeRates._deduplicate_specs(
        all_expanded, bi_bi)
    @test length(deduped) <= length(all_expanded)
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `_deduplicate_specs` not defined

- [ ] **Step 3: Implement _deduplicate_specs**

In `src/beam_enumeration.jl`:

```julia
# ─── Deduplication ────────────────────────────────────────────

"""
    _deduplicate_specs(specs, reaction) → Vector{MechanismSpec}

Remove duplicate mechanisms within a level. Two mechanisms are
duplicates if they have the same concentration fingerprint and
constraint descriptor. Keeps the one with lowest param_count.
"""
function _deduplicate_specs(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    isempty(specs) && return specs

    best = Dict{_DedupKey, MechanismSpec}()
    for spec in specs
        steps = spec.steps
        fp = _runtime_denominator_monomials(spec)

        groups = Dict{Tuple{Symbol,Bool}, Vector{Int}}()
        for (i, s) in enumerate(steps)
            met = step_metabolite(s)
            met === nothing && continue
            key = (met, s.is_equilibrium)
            push!(get!(groups, key, Int[]), i)
        end
        valid_groups = sort!(
            [sort!(g) for (_, g) in groups
             if length(g) >= 2];
            by=first)

        constraint_mask = _constraints_to_mask(
            spec.param_constraints, valid_groups, steps)
        desc = _constraint_descriptor(
            steps, valid_groups, constraint_mask)

        dedup_key = (fp, desc)
        if !haskey(best, dedup_key) ||
                spec.param_count < best[dedup_key].param_count
            best[dedup_key] = spec
        end
    end
    collect(values(best))
end
```

Note: This reuses `_concentration_fingerprint`, `_compute_re_partition_from_steps`, `_constraint_descriptor`, `_constraints_to_mask`, and `_DedupKey` from `mechanism_enumeration.jl` — these are shared utilities.

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl
git commit -m "Add _deduplicate_specs for within-level deduplication"
```

### Task 10: expand_mechanisms_same_param_count — test and implement

**Files:**
- Modify: `src/beam_enumeration.jl`
- Modify: `test/test_beam_enumeration.jl`

+0 dead-end additions where equivalence constraints exactly cancel the new parameters.

- [ ] **Step 1: Write failing test**

```julia
@testset "expand_mechanisms_same_param_count" begin
    # Ordered bi-bi: E→EA→EAB→EPQ→EQ→E
    # E_Q + A ⇌ E_QA where K_A = K_A_catalytic → +0 params
    topos = EnzymeRates._catalytic_topologies(bi_bi)
    spec = topos[1]  # ordered bi-bi
    original_pc = spec.param_count

    results = EnzymeRates.expand_mechanisms_same_param_count(
        [spec], bi_bi)

    # All results should have same param_count
    for r in results
        @test r.param_count == original_pc
    end

    # Should find at least one +0 dead-end
    # (E_Q + A ⇌ E_A_Q is a valid +0 addition)
    @test !isempty(results)

    # Verify param_count matches compiled mechanism
    for r in results
        m = compile_mechanism(r)
        @test length(parameters(m)) == r.param_count
    end

    # Fixed-point guarantee: calling again on the union
    # should produce no additional mechanisms
    all_input = vcat([spec], results)
    results2 = EnzymeRates.expand_mechanisms_same_param_count(
        all_input, bi_bi)
    all_specs = vcat(all_input, results2)
    deduped = EnzymeRates._deduplicate_specs(all_specs, bi_bi)
    @test length(deduped) == length(
        EnzymeRates._deduplicate_specs(all_input, bi_bi))
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — function not defined

- [ ] **Step 3: Implement expand_mechanisms_same_param_count**

```julia
# ─── Expansion: Same Param Count (+0) ────────────────────────

"""
    expand_mechanisms_same_param_count(specs, reaction)
        → Vector{MechanismSpec}

Add dead-end configurations that result in +0 net parameter change.
Iterates to a fixed point: each pass may create new forms that
enable further +0 additions. Returns the union of input specs
and all discovered +0 variants (deduplicated).
"""
function expand_mechanisms_same_param_count(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    all_specs = copy(specs)
    prev_count = 0
    while length(all_specs) != prev_count
        prev_count = length(all_specs)
        new_zero = MechanismSpec[]
        for spec in all_specs
            _expand_dead_end_zero!(
                new_zero, spec, reaction)
        end
        append!(all_specs, new_zero)
        all_specs = _deduplicate_specs(
            all_specs, reaction)
    end
    # Return only the new specs (not input)
    filter(s -> s ∉ specs, all_specs)
end

"""Generate dead-end configurations with +0 param change.
Uses shared helper `_expand_dead_end_at_delta!`."""
function _expand_dead_end_zero!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    _expand_dead_end_at_delta!(result, spec, reaction, 0)
end

"""
Add equivalence constraints between new dead-end binding step
and existing catalytic step binding the same metabolite.
"""
function _add_catalytic_equivalence!(
    new_steps, new_constraints, orig_steps,
    selected_forms, bound, de_met,
)
    # Find catalytic steps that bind the same metabolite
    for (i, s) in enumerate(orig_steps)
        s.is_equilibrium || continue
        cat_met = step_metabolite(s)
        cat_met == de_met || continue

        # Find new steps binding de_met to selected forms
        for (j, ns) in enumerate(new_steps)
            j <= length(orig_steps) && continue
            ns_met = step_metabolite(ns)
            ns_met == de_met || continue

            # Add constraint: new step K = catalytic K
            push!(new_constraints, (
                Symbol("K$j"), 1,
                [(Symbol("K$i"), 1)]))
            break  # One constraint per catalytic match
        end
    end
end
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl
git commit -m "Add expand_mechanisms_same_param_count (+0 dead-end)"
```

---

## Chunk 6: Allosteric Expansion (+2) and TR Equivalence

### Task 11: expand_mechanisms_by_two_params — test and implement

**Files:**
- Modify: `src/beam_enumeration.jl`
- Modify: `test/test_beam_enumeration.jl`

Converts a `MechanismSpec` to `AllostericMechanismSpec` with L + one K_T≠K_R.

- [ ] **Step 1: Write failing test**

```julia
@testset "expand_mechanisms_by_two_params" begin
    # Uni-uni + allosteric R
    topos = EnzymeRates._catalytic_topologies(uni_uni)
    spec = topos[1]

    results = EnzymeRates.expand_mechanisms_by_two_params(
        [spec], uni_uni_allosteric_R)

    @test !isempty(results)

    # All results should be AllostericMechanismSpec
    for r in results
        @test r isa EnzymeRates.AllostericMechanismSpec
    end

    # Each result should have exactly one non-TR-equivalent
    # metabolite
    for r in results
        t_mets = EnzymeRates._collect_t_state_metabolites(r)
        n_non_equiv = length(t_mets) -
            length(r.tr_equiv_metabolites)
        @test n_non_equiv == 1
    end

    # Verify with compilation and param count
    for r in results
        m = compile_mechanism(r)
        # AllostericMechanismSpec param_count = base +
        # allosteric additions
        base_pc = r.base.param_count
        n_t_mets = length(
            EnzymeRates._collect_t_state_metabolites(r))
        n_equiv = length(r.tr_equiv_metabolites)
        # L + (n_t_mets - n_equiv) K_T params
        allo_pc = 1 + (n_t_mets - n_equiv)
        @test length(parameters(m)) == base_pc + allo_pc
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — function not defined

- [ ] **Step 3: Implement expand_mechanisms_by_two_params**

```julia
# ─── Expansion: +2 Parameters (Allosteric) ───────────────────

"""
    expand_mechanisms_by_two_params(specs, reaction;
        max_catalytic_n=4) → Vector{AllostericMechanismSpec}

Convert base mechanisms to allosteric with L + one K_T≠K_R.
Generates all variants of which metabolite has K_T≠K_R,
all catalytic_n values, and all regulator site partitions.
Keeps only those at exactly +2 params.
"""
function expand_mechanisms_by_two_params(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    max_catalytic_n::Int=4,
)
    roles = regulator_roles(reaction)
    allo_regs = Symbol[
        r[1] for r in roles
        if r[2] == :allosteric || r[2] == :unknown]
    isempty(allo_regs) && return AllostericMechanismSpec[]

    result = AllostericMechanismSpec[]
    partitions = _set_partitions(allo_regs)

    for spec in specs
        # Collect metabolites that would have T-state params
        t_mets = Symbol[]
        for s in spec.steps
            s.is_equilibrium || continue
            met = step_metabolite(s)
            met !== nothing && met ∉ t_mets &&
                push!(t_mets, met)
        end
        for reg in allo_regs
            reg ∉ t_mets && push!(t_mets, reg)
        end

        for partition in partitions
            n_groups = length(partition)
            for cn in 1:max_catalytic_n
                for combo in Iterators.product(
                        ntuple(_ -> 1:cn, n_groups)...)
                    # Try each single metabolite as
                    # non-TR-equivalent
                    for (mi, non_equiv_met) in
                            enumerate(t_mets)
                        tr_equiv = Symbol[
                            m for m in t_mets
                            if m != non_equiv_met]
                        allo = AllostericMechanismSpec(
                            spec, cn, partition,
                            collect(combo), tr_equiv)
                        # Verify param_count = base + 2
                        n_non_equiv = length(t_mets) -
                            length(tr_equiv)
                        allo_addition = 1 + n_non_equiv
                        allo_addition == 2 || continue
                        push!(result, allo)
                    end
                end
            end
        end
    end
    result
end
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl
git commit -m "Add expand_mechanisms_by_two_params (allosteric +2)"
```

### Task 12: Remove TR equivalence move (allosteric +1) — test and implement

**Files:**
- Modify: `src/beam_enumeration.jl`
- Modify: `test/test_beam_enumeration.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "Remove TR equivalence (+1, allosteric)" begin
    # Create an allosteric mechanism with 2 TR-equivalent mets
    topos = EnzymeRates._catalytic_topologies(uni_uni)
    spec = topos[1]

    allo_results = EnzymeRates.expand_mechanisms_by_two_params(
        [spec], uni_uni_allosteric_R)
    @test !isempty(allo_results)

    # Pick one with TR-equivalent metabolites
    allo_spec = first(
        r for r in allo_results
        if !isempty(r.tr_equiv_metabolites))

    n_equiv = length(allo_spec.tr_equiv_metabolites)

    # expand_mechanisms_by_one_param should handle
    # AllostericMechanismSpec and produce TR-removal candidates
    results = EnzymeRates.expand_mechanisms_by_one_param(
        [allo_spec], uni_uni_allosteric_R)

    tr_removed = filter(results) do r
        r isa EnzymeRates.AllostericMechanismSpec &&
        length(r.tr_equiv_metabolites) ==
            n_equiv - 1
    end

    # One candidate per TR-equivalent metabolite
    @test length(tr_removed) == n_equiv
end
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL — `expand_mechanisms_by_one_param` doesn't accept `AllostericMechanismSpec`

- [ ] **Step 3: Implement AllostericMechanismSpec support in expand_mechanisms_by_one_param**

Add an overload and the TR removal function:

```julia
"""expand_mechanisms_by_one_param for AllostericMechanismSpec."""
function expand_mechanisms_by_one_param(
    specs::Vector{AllostericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = AllostericMechanismSpec[]
    for spec in specs
        _expand_remove_tr_equiv!(result, spec, reaction)
    end
    result
end

"""Remove one TR equivalence: make one metabolite's K_T ≠ K_R."""
function _expand_remove_tr_equiv!(
    result::Vector{AllostericMechanismSpec},
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    for (i, met) in enumerate(spec.tr_equiv_metabolites)
        new_equiv = [
            spec.tr_equiv_metabolites[j]
            for j in eachindex(spec.tr_equiv_metabolites)
            if j != i]
        push!(result, AllostericMechanismSpec(
            spec.base, spec.catalytic_n,
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities,
            new_equiv))
    end
end
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl
git commit -m "Add remove-TR-equivalence move for allosteric mechanisms"
```

---

## Chunk 7: enumerate_mechanisms Orchestrator and Equivalence Tests

### Task 13: enumerate_mechanisms orchestrator — test and implement

**Files:**
- Modify: `src/beam_enumeration.jl`
- Modify: `test/test_beam_enumeration.jl`

The main loop that ties everything together.

- [ ] **Step 1: Write failing test**

```julia
@testset "enumerate_mechanisms (beam)" begin
    # Uni-uni: should produce mechanisms at multiple
    # param_count levels
    iter = EnzymeRates.enumerate_mechanisms(uni_uni)
    mechs = collect(iter)
    @test !isempty(mechs)

    # All mechanisms should have valid param_count
    for m in mechs
        @test m.param_count >= 0
    end

    # Should include mechanisms at different param levels
    param_counts = Set(m.param_count for m in mechs
        if m isa EnzymeRates.MechanismSpec)
    @test length(param_counts) >= 2
end
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL — `enumerate_mechanisms` not defined (old one was renamed)

- [ ] **Step 3: Implement enumerate_mechanisms**

```julia
# ─── Orchestrator ─────────────────────────────────────────────

"""
    enumerate_mechanisms(reaction; max_param_count=nothing,
        max_catalytic_n=4) → MechanismIterator

Enumerate all valid mechanisms for a reaction, expanding
level-by-level by param_count. At each level:
1. Merge catalytic seeds + cached +2 specs + expanded +1 specs
2. Apply expand_mechanisms_same_param_count to fixed point
3. Deduplicate
4. Yield this level's mechanisms
5. Expand by +1 and +2 for next levels
"""
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    max_param_count::Union{Nothing,Int}=nothing,
    max_catalytic_n::Int=4,
)
    catalytic = _catalytic_topologies(reaction)

    # Group catalytic topologies by param_count
    seeds_by_pc = Dict{Int, Vector{MechanismSpec}}()
    for spec in catalytic
        push!(get!(seeds_by_pc, spec.param_count,
            MechanismSpec[]), spec)
    end

    min_pc = minimum(keys(seeds_by_pc))
    max_pc = if max_param_count !== nothing
        max_param_count
    else
        # Heuristic upper bound: enough levels to cover
        # all RE→SS + dead-end + allosteric expansions
        min_pc + 50
    end

    cache = Dict{Int, Vector{AbstractMechanismSpec}}()
    all_results = AbstractMechanismSpec[]
    current_plus_one = MechanismSpec[]

    for pc in min_pc:max_pc
        # Assemble this level
        level = MechanismSpec[]
        append!(level, get(seeds_by_pc, pc, MechanismSpec[]))
        append!(level, current_plus_one)

        # Add cached specs (extract both types before deletion)
        cached_at_pc = get(cache, pc, AbstractMechanismSpec[])
        for spec in cached_at_pc
            if spec isa MechanismSpec
                push!(level, spec)
            end
        end

        # Check for allosteric specs at this level
        allo_at_level = AllostericMechanismSpec[]
        for spec in cached_at_pc
            if spec isa AllostericMechanismSpec
                push!(allo_at_level, spec)
            end
        end
        delete!(cache, pc)

        # Process allosteric specs even if level is empty
        if !isempty(allo_at_level)
            allo_deduped = _deduplicate_allosteric(
                allo_at_level, reaction)
            append!(all_results, allo_deduped)

            allo_plus_one = expand_mechanisms_by_one_param(
                allo_deduped, reaction)
            for spec in allo_plus_one
                push!(get!(cache, pc + 1,
                    AbstractMechanismSpec[]), spec)
            end
        end

        isempty(level) && continue

        # +0 expansion (fixed point handled internally)
        new_zero = expand_mechanisms_same_param_count(
            level, reaction)
        append!(level, new_zero)
        level = _deduplicate_specs(level, reaction)

        # Yield this level
        append!(all_results, level)

        # Expand +1
        current_plus_one = expand_mechanisms_by_one_param(
            level, reaction)
        # Filter to actual +1
        filter!(
            s -> s.param_count == pc + 1,
            current_plus_one)

        # Expand +2 (allosteric)
        plus_two = expand_mechanisms_by_two_params(
            level, reaction; max_catalytic_n)
        for spec in plus_two
            push!(get!(cache, pc + 2,
                AbstractMechanismSpec[]), spec)
        end

        # (allosteric handling already done above)

        # Termination: if no new mechanisms from any source
        has_future_seeds = any(
            haskey(seeds_by_pc, k) for k in (pc+1):max_pc)
        has_future_cache = any(
            haskey(cache, k) for k in (pc+1):max_pc)
        if isempty(current_plus_one) &&
                !has_future_seeds && !has_future_cache
            break
        end
    end

    MechanismIterator(all_results, length(all_results))
end
```

- [ ] **Step 4: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add src/beam_enumeration.jl test/test_beam_enumeration.jl
git commit -m "Add enumerate_mechanisms beam search orchestrator"
```

### Task 14: Equivalence test against old_enumerate_mechanisms

**Files:**
- Modify: `test/test_beam_enumeration.jl`

For small reactions where `old_enumerate_mechanisms` works, verify the beam search produces the same mechanism set (by concentration fingerprint).

- [ ] **Step 1: Write equivalence test**

```julia
@testset "Equivalence with old_enumerate_mechanisms" begin
    # Helper: extract fingerprint set from mechanisms
    function fingerprint_set(specs)
        fps = Set()
        for spec in specs
            if spec isa EnzymeRates.MechanismSpec
                steps = spec.steps
                partition =
                    EnzymeRates._compute_re_partition_from_steps(
                        steps)
                fp = EnzymeRates._concentration_fingerprint(
                    steps, partition)
                push!(fps, fp)
            end
        end
        fps
    end

    for (name, rxn) in [
            ("uni-uni", uni_uni),
            ("uni-bi", uni_bi),
            ("uni-uni + I", uni_uni_dead_end_I),
            ("uni-uni unknown reg", uni_uni_reg_unknown),
        ]
        @testset "$name" begin
            old = collect(
                EnzymeRates.old_enumerate_mechanisms(rxn))
            new = collect(
                EnzymeRates.enumerate_mechanisms(rxn))

            old_base = filter(
                s -> s isa EnzymeRates.MechanismSpec, old)
            new_base = filter(
                s -> s isa EnzymeRates.MechanismSpec, new)

            old_fps = fingerprint_set(old_base)
            new_fps = fingerprint_set(new_base)

            # New should be a superset (may include
            # multi-level dead-end that old doesn't have)
            @test issubset(old_fps, new_fps)
        end
    end
end
```

- [ ] **Step 2: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All equivalence tests pass. If failures occur, debug by comparing the mechanisms produced by each pipeline.

- [ ] **Step 3: Commit**

```bash
git add test/test_beam_enumeration.jl
git commit -m "Add equivalence test: beam vs old_enumerate_mechanisms"
```

### Task 15: Final integration test and cleanup

**Files:**
- Modify: `test/test_beam_enumeration.jl`
- Modify: `src/beam_enumeration.jl` (ABOUTME, dead code removal)

- [ ] **Step 1: Add catalytic topology tests (copied from old tests)**

Copy the Stage 1 catalytic topology tests from `test/old_test_mechanism_enumeration.jl` into `test/test_beam_enumeration.jl`. These test `_catalytic_topologies` which is shared code, but we want to ensure the beam tests independently validate the seed generation. Copy the uni-uni, uni-bi, bi-bi, and bi-bi ping-pong topology count tests.

- [ ] **Step 2: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass (old and new)

- [ ] **Step 3: Review src/beam_enumeration.jl for dead code and cleanup**

Read through the file. Remove any unused functions. Ensure ABOUTME comments are present. Verify line length ≤ 92 chars, 4-space indentation.

- [ ] **Step 4: Commit**

```bash
git add test/test_beam_enumeration.jl src/beam_enumeration.jl
git commit -m "Add catalytic topology tests and cleanup beam_enumeration.jl"
```

- [ ] **Step 5: Run full test suite one final time**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass

- [ ] **Step 6: Final commit if any cleanup needed**

```bash
git add -A
git commit -m "Final cleanup for beam enumeration implementation"
```
