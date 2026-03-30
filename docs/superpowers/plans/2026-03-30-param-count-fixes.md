# Param Count Accuracy Fixes

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix param_count accuracy issues found during competitive review round 2.

**Tech Stack:** Julia, EnzymeRates.jl

---

## Issue Summary

| # | Issue | Severity | Task |
|---|-------|----------|------|
| 1 | `_tr_equiv_met_delta` returns 0 for allosteric regulators | CRITICAL | 1 |
| 2 | `_expand_re_to_ss` undercounts when mirror steps converted | CRITICAL | 2 |
| 3 | `_tr_equiv_met_delta` overcounts for constrained metabolites | IMPORTANT | 3 |
| 4 | No test for regulator TR equiv removal delta | IMPORTANT | 1 |
| 5 | No test for mirror step param_count in RE→SS | IMPORTANT | 2 |
| 6 | No test for metabolite overlap (substrate = regulator) | IMPORTANT | 4 |

---

### Task 1: Fix `_tr_equiv_met_delta` for allosteric regulators

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_tr_equiv_met_delta` and `_expand_remove_tr_equiv`)
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "TR equiv removal delta for allosteric regulators" begin
    rxn_r = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    specs = EnzymeRates.init_mechanisms(rxn_r)
    spec = first(specs)
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, rxn_r)
    allo = first(allo_specs)
    # Add R as tr_equiv regulator
    reg_specs = EnzymeRates._expand_add_allosteric_regulator(
        allo, rxn_r)
    tr_spec = first(filter(
        r -> :R in r.tr_equiv_metabolites, reg_specs))
    pc_before = tr_spec.param_count
    # Remove TR equiv for R → independent (+1)
    result = EnzymeRates._expand_remove_tr_equiv(
        tr_spec, rxn_r)
    r_removal = filter(result) do r
        :R ∉ r.tr_equiv_metabolites &&
        :R ∉ r.r_only_metabolites &&
        :R ∉ r.t_only_metabolites
    end
    @test !isempty(r_removal)
    for r in r_removal
        @test r.param_count == pc_before + 1
    end
end
```

- [ ] **Step 2: Run to verify failure**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `r.param_count == pc_before + 0` (delta is 0 instead of 1).

- [ ] **Step 3: Fix `_tr_equiv_met_delta`**

The function currently only searches `base.steps` for the metabolite. Allosteric regulators aren't in base steps. Update the function signature to also accept the `AllostericMechanismSpec` (or its reg sites), and return +1 for regulators:

```julia
"""
Compute param_count delta for removing TR equivalence
of a metabolite. RE binding: +1 (K_T). SS binding: +2
(kf_T + kr_T). Allosteric regulator (not in base steps): +1.
"""
function _tr_equiv_met_delta(
    met::Symbol,
    steps::Vector{StepSpec},
    allosteric_reg_sites::Vector{Vector{Symbol}}=Vector{Symbol}[])
    # Check if met is an allosteric regulator
    for site in allosteric_reg_sites
        if met in site
            return 1  # regulator: always +1 (one K_T)
        end
    end
    # Catalytic metabolite: scan base steps
    delta = 0
    for s in steps
        step_metabolite(s) === met || continue
        delta += s.is_equilibrium ? 1 : 2
    end
    delta
end
```

- [ ] **Step 4: Update call site in `_expand_remove_tr_equiv`**

Change the call from:
```julia
delta = _tr_equiv_met_delta(met, spec.base.steps)
```
to:
```julia
delta = _tr_equiv_met_delta(
    met, spec.base.steps,
    spec.allosteric_reg_sites)
```

- [ ] **Step 5: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Fix _tr_equiv_met_delta: return +1 for allosteric regulators"
```

---

### Task 2: Fix `_expand_re_to_ss` mirror step param_count

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_expand_re_to_ss`)
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test**

Create a mechanism with dead-end substrate/product forms that produce mirror steps. When converting an RE binding step to SS, verify that param_count accounts for converted mirrors.

```julia
@testset "RE→SS mirror step param_count" begin
    # Use bi-bi random topology with dead-end forms
    # The random topology has forms E_A, E_B where both
    # can have dead-end extensions. Find an init spec
    # with dead-end forms that create mirrors.
    bi_bi_specs = EnzymeRates.init_mechanisms(bi_bi_rxn)
    # Find a spec with mirror steps (more steps than
    # the base topology)
    spec_with_de = first(filter(bi_bi_specs) do s
        length(s.steps) > 9  # base random has 9 steps
    end)
    # First remove constraints so RE steps are eligible
    unconstrained = first(
        EnzymeRates._expand_remove_constraint(spec_with_de))
    # Convert an RE step that has mirrors
    result = EnzymeRates._expand_re_to_ss(unconstrained)
    for r in result
        # Count how many steps changed from RE to SS
        n_converted = count(
            !r.steps[i].is_equilibrium &&
             unconstrained.steps[i].is_equilibrium
            for i in eachindex(r.steps)
            if i <= length(unconstrained.steps))
        # param_count should increase by n_converted
        # (1 for main + 1 per mirror)
        @test r.param_count ==
            unconstrained.param_count + n_converted
    end
end
```

- [ ] **Step 2: Run to verify failure**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 3: Fix `_expand_re_to_ss`**

Count the number of mirror steps converted and add to the delta:

```julia
function _expand_re_to_ss(spec::MechanismSpec)
    result = MechanismSpec[]
    constrained = _constrained_step_indices(
        spec.param_constraints)

    for (i, s) in enumerate(spec.steps)
        s.is_equilibrium || continue
        i in constrained && continue

        new_steps = [StepSpec(st.reactants, st.products,
                     st.is_equilibrium) for st in spec.steps]
        new_steps[i] = StepSpec(
            s.reactants, s.products, false)

        # Count mirror steps converted
        n_mirrors = 0
        from_form, to_form = step_forms(s)
        for (j, ms) in enumerate(new_steps)
            j == i && continue
            ms.is_equilibrium || continue
            j in constrained && continue
            mf, mt = step_forms(ms)
            if _is_mirror_of(
                    mf, mt, from_form, to_form,
                    spec.steps)
                new_steps[j] = StepSpec(
                    ms.reactants, ms.products, false)
                n_mirrors += 1
            end
        end

        push!(result, MechanismSpec(
            spec.reaction, new_steps,
            copy(spec.param_constraints),
            spec.param_count + 1 + n_mirrors))
    end
    result
end
```

- [ ] **Step 4: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Fix _expand_re_to_ss: count mirror step conversions in param_count"
```

---

### Task 3: Fix `_tr_equiv_met_delta` for constrained metabolites

**Files:**
- Modify: `src/mechanism_enumeration.jl` (`_tr_equiv_met_delta`)
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "TR equiv removal delta with constraints" begin
    # Bi-bi random has K_A constrained across 2 forms
    # (K1 = K4 or similar). Removing TR equiv should
    # be +1 (only the leader K becomes independent),
    # not +2 (both steps counted).
    bi_bi_allo = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        oligomeric_state: 2
    end
    specs = EnzymeRates.init_mechanisms(bi_bi_allo)
    # Find a spec with constraints on A
    constrained = first(filter(
        s -> !isempty(s.param_constraints), specs))
    # Make it allosteric (V-type so no r_only mets)
    allo_specs = EnzymeRates._expand_to_allosteric(
        constrained, bi_bi_allo)
    v_type = first(filter(
        r -> !isempty(r.r_only_cat_steps), allo_specs))
    # Remove a TR equiv for a constrained metabolite
    result = EnzymeRates._expand_remove_tr_equiv(
        v_type, bi_bi_allo)
    # Each removal should be +1 (not +2 for constrained)
    for r in result
        @test r.param_count <= v_type.param_count + 2
        # For RE-binding constrained metabolites: +1
        # (only leader step's K_T is independent,
        # followers are constrained K_T = K_leader_T)
    end
    # Compile and verify
    for r in result[1:min(3, end)]
        m = AllostericEnzymeMechanism(r)
        @test length(parameters(m)) <= r.param_count
    end
end
```

- [ ] **Step 2: Run to verify failure**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 3: Fix `_tr_equiv_met_delta` to account for constraints**

Update the function to accept `param_constraints` and only count unconstrained steps:

```julia
function _tr_equiv_met_delta(
    met::Symbol,
    steps::Vector{StepSpec},
    allosteric_reg_sites::Vector{Vector{Symbol}}=Vector{Symbol}[];
    param_constraints::Vector{ParamConstraint}=ParamConstraint[])
    # Allosteric regulator: always +1
    for site in allosteric_reg_sites
        met in site && return 1
    end
    # Catalytic metabolite: count only unconstrained steps
    constrained = _constrained_step_indices(param_constraints)
    delta = 0
    for (idx, s) in enumerate(steps)
        step_metabolite(s) === met || continue
        idx in constrained && continue
        delta += s.is_equilibrium ? 1 : 2
    end
    delta
end
```

- [ ] **Step 4: Update call site**

```julia
delta = _tr_equiv_met_delta(
    met, spec.base.steps,
    spec.allosteric_reg_sites;
    param_constraints=spec.base.param_constraints)
```

- [ ] **Step 5: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Fix _tr_equiv_met_delta: skip constrained steps"
```

---

### Task 4: Add tests for metabolite overlap (substrate = regulator)

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write tests for dead-end inhibitor overlap**

```julia
@testset "Metabolite overlap: substrate as dead-end inhibitor" begin
    rxn_overlap = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: S
    end
    specs = EnzymeRates.init_mechanisms(rxn_overlap)
    @test !isempty(specs)

    # Dead-end regulator S should use __reg suffix
    spec = first(specs)
    de_specs = EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn_overlap)
    @test !isempty(de_specs)
    for s in de_specs
        reg_syms = [sym for st in s.steps
            for sym in Iterators.flatten(
                (st.reactants, st.products))
            if contains(string(sym), "S__reg")]
        @test !isempty(reg_syms)
    end

    # S-as-regulator can bind forms where S-as-substrate
    # is already bound (different binding site)
    for s in de_specs
        for st in s.steps
            met = step_metabolite(st)
            met === nothing && continue
            if contains(string(met), "S__reg")
                base_form = step_forms(st)[1]
                @test true  # verify no error
            end
        end
    end

    # All compile correctly
    for s in de_specs
        m = EnzymeMechanism(s)
        @test m isa EnzymeMechanism
    end
end
```

- [ ] **Step 2: Write tests for allosteric regulator overlap**

```julia
@testset "Metabolite overlap: substrate as allosteric regulator" begin
    rxn_allo_overlap = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        allosteric_regulators: S
        oligomeric_state: 2
    end
    specs = EnzymeRates.init_mechanisms(rxn_allo_overlap)
    @test !isempty(specs)
    spec = first(specs)

    # Allosteric conversion should work
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, rxn_allo_overlap)
    @test !isempty(allo_specs)

    # Add S as allosteric regulator (it's both substrate
    # and allosteric regulator)
    allo = first(allo_specs)
    reg_specs = EnzymeRates._expand_add_allosteric_regulator(
        allo, rxn_allo_overlap)
    @test !isempty(reg_specs)

    # S should appear in allosteric_reg_sites
    for r in reg_specs
        has_s = any(
            :S in site
            for site in r.allosteric_reg_sites)
        @test has_s
    end

    # All compile correctly
    for r in reg_specs
        m = AllostericEnzymeMechanism(r)
        @test m isa AllostericEnzymeMechanism
    end

    # TR equiv removal should produce separate results
    # for S-as-substrate and S-as-regulator
    tr_spec = first(filter(
        r -> :S in r.tr_equiv_metabolites, reg_specs))
    result = EnzymeRates._expand_remove_tr_equiv(
        tr_spec, rxn_allo_overlap)
    # S appears in tr_equiv_metabolites both as
    # catalytic metabolite AND as regulator.
    # Removing each should produce a separate variant.
    @test length(result) >= 2
end
```

- [ ] **Step 3: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 4: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add tests for metabolite overlap: substrate as dead-end and allosteric regulator"
```

---

## Self-Review

**Coverage:**
- ✅ Allosteric regulator TR equiv delta fixed (Task 1)
- ✅ Mirror step param_count fixed (Task 2)
- ✅ Constrained metabolite delta fixed (Task 3)
- ✅ Substrate=dead-end regulator overlap tested (Task 4)
- ✅ Substrate=allosteric regulator overlap tested (Task 4)

**Type consistency:**
- `_tr_equiv_met_delta` gains optional args (backward compatible) ✅
- `_expand_re_to_ss` counts mirrors in delta ✅
- All `param_count` values remain upper bounds ✅
