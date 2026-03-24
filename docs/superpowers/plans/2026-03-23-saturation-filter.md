# Saturation Filter Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Filter out non-saturating RE/SS mechanism variants at the earliest possible point in the enumeration pipeline (stage 2), so that only biochemically reasonable rate equations are generated.

**Architecture:** Add a `_saturates_all_metabolites` predicate function that checks whether every substrate/product has at least one binding step that is either RE, or SS connecting different RE segments. Call this inside the existing mask loop in `_expand_ress_variants`, after the RE partition is computed but before allocating the `MechanismSpec`. Update all downstream test expected counts.

**Tech Stack:** Julia, EnzymeRates.jl internal API

---

## Background

The Cha method (rapid-equilibrium + steady-state hybrid) produces the rate equation denominator from two sources:
1. **Sigma (intra-segment)**: RE equilibria within each RE group — metabolites from RE binding steps appear here
2. **Cofactors (inter-segment)**: King-Altman determinants of the reduced graph — metabolites from SS steps connecting different RE segments appear here

A metabolite is absent from the denominator (and thus non-saturating) when ALL its binding steps are SS steps connecting forms within the SAME RE segment. Such steps become self-loops in the reduced Laplacian, and their metabolite contribution is lost from the cofactor determinant.

## Files

- **Modify:** `src/mechanism_enumeration.jl` — add `_saturates_all_metabolites` helper, integrate into `_expand_ress_variants`
- **Modify:** `test/test_mechanism_enumeration.jl` — add new tests for the helper, update expected counts throughout

---

### Task 1: Add `_saturates_all_metabolites` helper with tests

**Files:**
- Modify: `src/mechanism_enumeration.jl` (near line 100, after `_compute_re_partition_from_steps`)
- Modify: `test/test_mechanism_enumeration.jl` (new testset after "Stage 2" section, around line 640)

- [ ] **Step 1: Write failing tests for `_saturates_all_metabolites`**

Add a new testset in `test/test_mechanism_enumeration.jl` after the "Stage 2" testset (after line 640):

```julia
@testset "Saturation check" begin
    # ── Uni-Uni 3-step cycle ──
    # E+P⇌EP (RE), E+S⇌ES (RE), ES⇌EP (SS) → saturates
    @test EnzymeRates._saturates_all_metabolites(
        [EnzymeRates.StepSpec([:E, :P], [:E_P], true),
         EnzymeRates.StepSpec([:E, :S], [:E_S], true),
         EnzymeRates.StepSpec([:E_S], [:E_P], false)],
        [[:E, :E_P, :E_S]],  # one RE group
        Set([:S, :P]))

    # E+P⇌EP (SS), E+S⇌ES (RE), ES⇌EP (RE)
    # → P binding is SS within same RE group → P doesn't saturate
    @test !EnzymeRates._saturates_all_metabolites(
        [EnzymeRates.StepSpec([:E, :P], [:E_P], false),
         EnzymeRates.StepSpec([:E, :S], [:E_S], true),
         EnzymeRates.StepSpec([:E_S], [:E_P], true)],
        [[:E, :E_S, :E_P]],  # one RE group
        Set([:S, :P]))

    # All-SS Uni-Uni: each form is its own group → saturates
    @test EnzymeRates._saturates_all_metabolites(
        [EnzymeRates.StepSpec([:E, :P], [:E_P], false),
         EnzymeRates.StepSpec([:E, :S], [:E_S], false),
         EnzymeRates.StepSpec([:E_S], [:E_P], false)],
        [[:E], [:E_P], [:E_S]],  # three groups
        Set([:S, :P]))

    # Ordered Bi-Bi, SS step 1 only (E+A→EA):
    # all forms in one RE group, A binding is SS within it
    @test !EnzymeRates._saturates_all_metabolites(
        [EnzymeRates.StepSpec([:E, :A], [:E_A], false),
         EnzymeRates.StepSpec([:E_A, :B], [:E_A_B], true),
         EnzymeRates.StepSpec([:E, :Q], [:E_Q], true),
         EnzymeRates.StepSpec([:E_Q, :P], [:E_P_Q], true),
         EnzymeRates.StepSpec([:E_A_B], [:E_P_Q], true)],
        [[:E, :E_A, :E_A_B, :E_P_Q, :E_Q]],
        Set([:A, :B, :P, :Q]))

    # Ordered Bi-Bi, SS step 5 (isomerization):
    # all binding steps are RE → all saturate
    @test EnzymeRates._saturates_all_metabolites(
        [EnzymeRates.StepSpec([:E, :A], [:E_A], true),
         EnzymeRates.StepSpec([:E_A, :B], [:E_A_B], true),
         EnzymeRates.StepSpec([:E, :Q], [:E_Q], true),
         EnzymeRates.StepSpec([:E_Q, :P], [:E_P_Q], true),
         EnzymeRates.StepSpec([:E_A_B], [:E_P_Q], false)],
        [[:E, :E_A, :E_A_B, :E_P_Q, :E_Q]],
        Set([:A, :B, :P, :Q]))

    # Bi-Bi with SS steps 1+2 (E+A, EA+B):
    # Two RE segments: {EA} and {E,EQ,EPQ,EAB}
    # A: SS step1 connects E(seg2)→EA(seg1) = different → ok
    # B: SS step2 connects EA(seg1)→EAB(seg2) = different → ok
    @test EnzymeRates._saturates_all_metabolites(
        [EnzymeRates.StepSpec([:E, :A], [:E_A], false),
         EnzymeRates.StepSpec([:E_A, :B], [:E_A_B], false),
         EnzymeRates.StepSpec([:E, :Q], [:E_Q], true),
         EnzymeRates.StepSpec([:E_Q, :P], [:E_P_Q], true),
         EnzymeRates.StepSpec([:E_A_B], [:E_P_Q], true)],
        [[:E_A], [:E, :E_Q, :E_P_Q, :E_A_B]],
        Set([:A, :B, :P, :Q]))

    # Empty metabolites set → trivially saturates
    @test EnzymeRates._saturates_all_metabolites(
        [EnzymeRates.StepSpec([:E, :S], [:E_S], false)],
        [[:E], [:E_S]],
        Set{Symbol}())
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — `_saturates_all_metabolites` is not defined

- [ ] **Step 3: Implement `_saturates_all_metabolites`**

Add in `src/mechanism_enumeration.jl` after `_compute_re_partition_from_steps` (around line 134), before the "Concentration Fingerprint" section:

```julia
"""
    _saturates_all_metabolites(steps, partition, metabolites)
        -> Bool

Check that every metabolite in `metabolites` has at least one
binding step that is either RE, or SS connecting forms in
different RE segments. A metabolite whose only binding steps
are SS within the same RE segment will not appear in the rate
equation denominator and thus will not exhibit saturation.
"""
function _saturates_all_metabolites(
    steps::Vector{StepSpec},
    partition::Vector{Vector{Symbol}},
    metabolites::Set{Symbol},
)
    isempty(metabolites) && return true
    form_to_group = Dict{Symbol, Int}()
    for (gi, group) in enumerate(partition)
        for f in group
            form_to_group[f] = gi
        end
    end
    for met in metabolites
        has_saturating_step = false
        for s in steps
            step_metabolite(s) == met || continue
            if s.is_equilibrium ||
                    form_to_group[s.reactants[1]] !=
                    form_to_group[s.products[1]]
                has_saturating_step = true
                break
            end
        end
        has_saturating_step || return false
    end
    true
end
```

- [ ] **Step 4: Run tests to verify the new testset passes**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: The "Saturation check" tests PASS. Other tests still pass at their original counts (filter not integrated yet).

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add _saturates_all_metabolites predicate for saturation filtering"
```

---

### Task 2: Integrate saturation filter into `_expand_ress_variants`

**Files:**
- Modify: `src/mechanism_enumeration.jl:864-903` (the `_expand_ress_variants` function)

- [ ] **Step 1: Write a failing test for filtered stage 2 output**

Add a test inside the existing "Stage 2: RE/SS expansion" testset (around line 536, after the Uni-Uni count test) that verifies saturation filtering:

```julia
    @testset "Saturation filtering" begin
        # All single-SS variants of a simple cycle where
        # the SS binding step's metabolite is in the same
        # RE group should be filtered out. Only the
        # isomerization-SS variant survives.
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
        topo = mechanism_spec_from_mechanism(
            m_uu, uni_uni)
        result = EnzymeRates._expand_ress_variants(
            [topo], uni_uni)
        # Every surviving variant must saturate
        for s in result
            partition =
                EnzymeRates._compute_re_partition_from_steps(
                    s.steps)
            sub_names = Set(
                n[1] for n in EnzymeRates.substrates(
                    uni_uni))
            prod_names = Set(
                n[1] for n in EnzymeRates.products(
                    uni_uni))
            @test EnzymeRates._saturates_all_metabolites(
                s.steps, partition,
                union(sub_names, prod_names))
        end
    end
```

- [ ] **Step 2: Run test to verify it fails (filter not yet active)**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL — 2 of 7 current Uni-Uni variants don't saturate, so the assertion fails. This is our red test.

- [ ] **Step 3: Integrate the filter into `_expand_ress_variants`**

Modify `_expand_ress_variants` in `src/mechanism_enumeration.jl`. Add a precomputation of metabolite names before the mask loop, and the saturation check after the RE partition check.

The function should look like this (changes marked with comments):

```julia
function _expand_ress_variants(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    max_re_groups::Int=7,
)
    # Precompute substrate/product metabolite names
    met_names = Set{Symbol}(
        s[1] for s in substrates(reaction))
    union!(met_names,
        Set{Symbol}(p[1] for p in products(reaction)))

    result = MechanismSpec[]
    for spec in specs
        n = length(spec.steps)
        n == 0 && continue

        for mask in 1:(1 << n) - 1
            steps = [
                StepSpec(
                    s.reactants, s.products,
                    (mask >> (i - 1)) & 1 == 0,
                )
                for (i, s) in enumerate(spec.steps)
            ]
            # At least one step must be SS
            any(!s.is_equilibrium for s in steps) ||
                continue
            # Check RE-connected groups ≤ max_re_groups
            partition =
                _compute_re_partition_from_steps(steps)
            length(partition) > max_re_groups && continue
            # Every substrate/product must saturate
            _saturates_all_metabolites(
                steps, partition, met_names) ||
                continue

            n_re = count(s.is_equilibrium for s in steps)
            n_ss = n - n_re
            n_forms = length(all_form_names(steps))
            n_thermo = n - n_forms + 1
            pc = n_re + 2 * n_ss - n_thermo + 2

            push!(result, MechanismSpec(
                spec.reaction, steps,
                spec.param_constraints, pc,
            ))
        end
    end
    result
end
```

- [ ] **Step 4: Run tests to verify the saturation filtering test passes**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: "Saturation filtering" test PASSES. Many other tests FAIL due to changed counts — this is expected and will be fixed in Task 3.

- [ ] **Step 5: Commit (WIP — tests have known failures due to count changes)**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "WIP: integrate saturation filter into _expand_ress_variants"
```

---

### Task 3: Update stage 2 test expected counts

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

The saturation filter reduces stage 2 output. Update these assertions:

- [ ] **Step 1: Update Uni-Uni stage 2 count**

Line 529: `@test length(result) == 7` → `@test length(result) == 5`

Two variants filtered: single-SS on P-binding (P unsaturated) and single-SS on S-binding (S unsaturated).

- [ ] **Step 2: Update Uni-Bi stage 2 count**

Line 559: `@test length(result) == 2^n_steps - 1` → `@test length(result) == 12`

4-step cycle: 15 total, 3 single-SS binding variants filtered (S, P, Q each unsaturated).

- [ ] **Step 3: Update Bi-Bi sequential stage 2 count**

Line 591: `@test length(result) == 2^n_steps - 1` → `@test length(result) == 27`

5-step cycle: 31 total, 4 single-SS binding variants filtered (A, B, P, Q each unsaturated).

- [ ] **Step 4: Update Bi-Bi random max_re_groups counts**

Line 626: `511` → `483`
Line 629: `241` → `213`
Line 632: `379` → `351`
Line 638: `89` → `69`

- [ ] **Step 5: Run stage 2 tests to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: Stage 2 tests pass. Other failures remain.

- [ ] **Step 6: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Update stage 2 test counts for saturation filter"
```

---

### Task 4: Update stage 5 (dedup) test expected counts

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Update Uni-Uni dedup count**

Line 1213: `@test length(deduped) == 3` → `@test length(deduped) == 1`

Only the all-binding-RE variant (fingerprint `{1, [S], [P]}`) survives. The `{1, [S]}` and `{1, [P]}` variants had one binding as SS within one RE group → filtered.

Update the comment above it (lines 1205-1212) to describe the surviving fingerprint:
```julia
        # 1 distinct concentration fingerprint survives:
        #   {1,[S],[P]}: both binding steps RE — classic
        #     reversible Michaelis-Menten denominator
```

- [ ] **Step 2: Update Bi-Bi sequential single-SS counts**

Line 1270: `@test length(single_ss) == 5` → `@test length(single_ss) == 1`
Line 1274: `@test length(deduped) == 5` → `@test length(deduped) == 1`

Only the isomerization-SS variant survives. Update comment at line 1269 accordingly:
```julia
        # Only isomerization-SS survives saturation filter
        @test length(single_ss) == 1
```

- [ ] **Step 3: Update Bi-Bi random dedup counts**

Line 1334: single-SS count stays `@test length(single_ss) == 9` — all 9 random-order single-SS variants saturate because every metabolite binds at 2+ steps (at least one always RE).

Line 1398: `@test length(ress) == 511` → `@test length(ress) == 483`
Line 1401: `@test length(deduped) == 146` → `@test length(deduped) == 126`

- [ ] **Step 4: Update Uni-Uni + I regulator dedup counts**

Line 1369: `@test length(eq) == 14` → `@test length(eq) == 10`
Line 1370: `@test length(deduped) == 6` → `@test length(deduped) == 2`

The Uni-Uni base topology shrinks from 7→5 at stage 2, which propagates through dead-end expansion and equivalence constraints.

- [ ] **Step 5: Run dedup tests to verify**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: Dedup tests pass.

- [ ] **Step 6: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Update dedup test counts for saturation filter"
```

---

### Task 5: Update end-to-end test expected counts

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Update end-to-end pipeline counts**

Line 2067: `@test length(result) == 3` → `@test length(result) == 1` (Uni-Uni)
Line 2073: `@test length(result) == 56` → `@test length(result) == 39` (Uni-Bi)
Line 2079: `@test length(stats.value) == 63762` → `@test length(stats.value) == 62034` (Bi-Bi)
Line 2088: `@test length(stats.value) == 64276` → `@test length(stats.value) == 62548` (Bi-Bi PP)
Line 2098: `@test length(result) == 17` → `@test length(result) == 7` (Uni-Uni + unknown reg)
Line 2105: `@test length(result) == 1012` → `@test length(result) == 683` (Uni-Bi + unknown reg)

- [ ] **Step 2: Update param_count accuracy test counts**

Line 2113: `@test length(all_specs) == 56` → `@test length(all_specs) == 39` (Uni-Bi)
Line 2118: `@test n_match == 56` → `@test n_match == 39`

- [ ] **Step 3: Update Sampled Bi-Bi param_count count**

Search for the "Sampled Bi-Bi specs" testset (around line 2139). It collects all Bi-Bi mechanisms and checks the total count:
Line ~2143: `@test length(all_specs) == 63762` → `@test length(all_specs) == 62034`

- [ ] **Step 4: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL tests pass. If any assertion still fails due to a count mismatch, compute the correct count by running the pipeline with the saturation filter and update accordingly.

- [ ] **Step 5: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Update end-to-end test counts for saturation filter"
```

---

### Task 6: Update docstring and CLAUDE.md

**Files:**
- Modify: `src/mechanism_enumeration.jl` (docstring for `_expand_ress_variants`)
- Modify: `.claude/CLAUDE.md` (pipeline documentation)

- [ ] **Step 1: Update `_expand_ress_variants` docstring**

Replace the docstring (lines 855-862) to mention saturation filtering:

```julia
"""
    _expand_ress_variants(specs, reaction; max_re_groups=7)
        -> Vector{MechanismSpec}

Enumerate all RE/SS assignment combinations for mechanism
steps. Filters out: (1) all-RE assignments (at least one
step must be SS), (2) assignments exceeding max_re_groups,
(3) assignments where any substrate/product fails to
saturate (all its binding steps are SS within the same
RE segment).
"""
```

- [ ] **Step 2: Update CLAUDE.md pipeline documentation**

In the "Mechanism enumeration staged pipeline" section, update stage 2 description to mention saturation filtering. Add after the pipeline order bullet:

```
- Stage 2 also filters non-saturating RE/SS assignments: every substrate/product must have at least one binding step that is RE or SS connecting different RE segments
```

- [ ] **Step 3: Run tests to verify nothing broke**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
git add src/mechanism_enumeration.jl .claude/CLAUDE.md
git commit -m "Document saturation filter in docstrings and CLAUDE.md"
```
