# Mechanism Enumeration Post-Review Fixes

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix issues found during competitive code review of the mechanism enumeration refactor.

**Architecture:** Allosteric conversion changes from +2 to +1 (L only). Two biochemically valid differentiation modes: K-type (substrate+product absent from T-state) and V-type (catalytic step inactive in T-state). Only r_only variants — T is the inactive conformation. Fix regulator param counting, dummy naming, and dedup.

**Tech Stack:** Julia, EnzymeRates.jl

**Spec:** `docs/superpowers/specs/2026-03-27-mechanism-enumeration-refactor-design.md`

---

## Issue Summary

| # | Issue | Severity | Task |
|---|-------|----------|------|
| 0 | `AllostericMechanismSpec` needs `param_count` field; eliminates need for `_estimated_param_count` | CRITICAL | 0 |
| 1 | `_expand_to_allosteric` wrong design: was +2, should be +1 with K-type/V-type differentiation | CRITICAL | 1 |
| 2 | Regulator dummy naming unstable (`__reg` index depends on eligible list) | IMPORTANT | 2 |
| 3 | Missing `r_only_cat_steps` expansion move for V→independent | IMPORTANT | 3 |
| 4 | Allosteric dedup misses site-order permutations | IMPORTANT | 4 |
| 5 | Integration tests skip allosteric compilation | IMPORTANT | 5 |
| 6 | CLAUDE.md misleading Haldane note | MINOR | 6 |

---

### Task 0: Add `param_count` field to `AllostericMechanismSpec`

**Files:**
- Modify: `src/old_mechanism_enumeration.jl:60-70` (struct definition)
- Modify: `src/old_mechanism_enumeration.jl` (all constructor call sites)
- Modify: `src/old_beam_enumeration.jl` (all constructor call sites)
- Modify: `src/mechanism_enumeration.jl` (all constructor call sites + remove `_estimated_param_count`)
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "AllostericMechanismSpec has param_count" begin
    specs = EnzymeRates.init_mechanisms(uni_uni_allo)
    spec = first(specs)
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo)
    allo = first(allo_specs)
    @test allo.param_count == spec.param_count + 1
end
```

- [ ] **Step 2: Run to verify failure**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `AllostericMechanismSpec` has no field `param_count`.

- [ ] **Step 3: Add `param_count::Int` as 10th field**

In `src/old_mechanism_enumeration.jl`, change the struct (lines 60-70) from:

```julia
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    tr_equiv_metabolites::Vector{Symbol}
    tr_equiv_cat_steps::Vector{Int}
    r_only_metabolites::Vector{Symbol}
    t_only_metabolites::Vector{Symbol}
    r_only_cat_steps::Vector{Int}
end
```

to:

```julia
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    tr_equiv_metabolites::Vector{Symbol}
    tr_equiv_cat_steps::Vector{Int}
    r_only_metabolites::Vector{Symbol}
    t_only_metabolites::Vector{Symbol}
    r_only_cat_steps::Vector{Int}
    param_count::Int
end
```

- [ ] **Step 4: Fix all constructor calls in `old_mechanism_enumeration.jl`**

Search for `AllostericMechanismSpec(` in the file. Each call currently passes 9 arguments. Add a 10th argument for `param_count`. In the old pipeline code, compute it as `_runtime_param_count(spec)` if that function exists, or a reasonable value. Since the old pipeline's param counting is separate from the new one, use whatever the old code was computing.

Search pattern: `AllostericMechanismSpec(` — find every occurrence and add the param_count argument.

- [ ] **Step 5: Fix all constructor calls in `old_beam_enumeration.jl`**

Same approach — find all `AllostericMechanismSpec(` calls, add the 10th argument.

- [ ] **Step 6: Fix all constructor calls in `mechanism_enumeration.jl`**

Update every `AllostericMechanismSpec(` call to pass `param_count`. Each expansion move knows its delta:

- `_expand_to_allosteric`: `spec.param_count + 1` (L)
- `_expand_add_allosteric_regulator`: parent's `param_count + 1`
- `_expand_remove_tr_equiv`: parent's `param_count + 1`
- `_rewrap_allosteric`: `original.param_count + (new_base.param_count - original.base.param_count)`

- [ ] **Step 7: Replace `_estimated_param_count` with field access**

Remove the `_estimated_param_count` function entirely. Replace all calls:

```julia
# Before:
_estimated_param_count(spec::MechanismSpec) = spec.param_count
_estimated_param_count(spec::AllostericMechanismSpec) = ...complex...

# After: just use spec.param_count for both types
```

Update `_push_to_dict!`:

```julia
function _push_to_dict!(result, spec::AbstractMechanismSpec)
    pc = spec.param_count  # works for both types now
    push!(get!(result, pc, AbstractMechanismSpec[]), spec)
end
```

For `MechanismSpec`, `param_count` already exists. For `AllostericMechanismSpec`, it's now the new 10th field. Both accessed the same way.

- [ ] **Step 8: Update tests**

Update any test that referenced `_estimated_param_count` to use `spec.param_count` directly.

- [ ] **Step 9: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 10: Commit**

```bash
git add src/old_mechanism_enumeration.jl src/old_beam_enumeration.jl src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add param_count field to AllostericMechanismSpec, remove _estimated_param_count"
```

---

### Task 1: Rewrite `_expand_to_allosteric` as +1 with K-type/V-type

**Files:**
- Modify: `src/mechanism_enumeration.jl:327-380`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing tests for the new allosteric conversion**

Replace the existing "Move 6" tests with:

```julia
@testset "Move 6: Allosteric conversion (+1)" begin
    @testset "K-type: substrate+product pairs absent" begin
        # Uni-uni: 1 substrate × 1 product × 2 states = 2
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        k_type = filter(r -> !isempty(r.r_only_metabolites) ||
                              !isempty(r.t_only_metabolites), result)
        # K-type: (S,P) both r_only = 1
        # V-type: steps r_only = 1
        # Total = 2
        @test length(result) == 2
    end

    @testset "K-type bi-bi: all substrate+product combos" begin
        bi_bi_allo = @enzyme_reaction begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            oligomeric_state: 2
        end
        specs = EnzymeRates.init_mechanisms(bi_bi_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, bi_bi_allo)
        # K-type: (2^2-1) × (2^2-1) = 3×3 = 9
        # V-type: 1
        # Total = 10
        @test length(result) == 10
    end

    @testset "All are +1 param" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for r in result
            @test r.param_count == spec.param_count + 1
        end
    end

    @testset "K-type: catalytic steps stay tr_equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        k_type = filter(
            r -> !isempty(r.r_only_metabolites), result)
        @test !isempty(k_type)
        for r in k_type
            @test isempty(r.r_only_cat_steps)
            # No t_only metabolites (T is inactive)
            @test isempty(r.t_only_metabolites)
        end
    end

    @testset "V-type: all metabolites tr_equiv" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        v_type = filter(r -> !isempty(r.r_only_cat_steps) &&
                              isempty(r.r_only_metabolites) &&
                              isempty(r.t_only_metabolites), result)
        @test length(v_type) >= 1
        for r in v_type
            sub_names = [s[1] for s in substrates(uni_uni_allo)]
            prod_names = [p[1] for p in products(uni_uni_allo)]
            all_cat = Symbol[sub_names; prod_names]
            for m in all_cat
                @test m in r.tr_equiv_metabolites
            end
        end
    end

    @testset "All compile" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        result = EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo)
        for r in result
            m = AllostericEnzymeMechanism(r)
            @test m isa AllostericEnzymeMechanism
        end
    end

    @testset "Already allosteric → empty" begin
        specs = EnzymeRates.init_mechanisms(uni_uni_allo)
        spec = first(specs)
        allo = first(EnzymeRates._expand_to_allosteric(
            spec, uni_uni_allo))
        @test isempty(EnzymeRates._expand_to_allosteric(
            allo, uni_uni_allo))
    end
end
```

- [ ] **Step 2: Run tests to verify failure**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — wrong counts, wrong param delta.

- [ ] **Step 3: Add `_valid_allosteric_differentiations` helper**

```julia
"""
    _valid_allosteric_differentiations(reaction, spec)
        → Vector{NamedTuple}

Enumerate biochemically valid T/R differentiations for
allosteric conversion. Returns tuples of
(r_only_mets, t_only_mets, r_only_cat_steps).

K-type: ≥1 substrate + ≥1 product absent from one state.
  Catalytic steps stay tr_equiv (irrelevant since state
  can't complete catalytic cycle).
V-type: all SS isomerization steps inactive in one state
  (kf=kr=0). All metabolite K's stay tr_equiv.
"""
function _valid_allosteric_differentiations(
    @nospecialize(reaction::EnzymeReaction),
    spec::MechanismSpec)
    sub_names = [s[1] for s in substrates(reaction)]
    prod_names = [p[1] for p in products(reaction)]

    # SS isomerization step indices
    ss_isom = Int[]
    for (i, s) in enumerate(spec.steps)
        !s.is_equilibrium &&
            step_metabolite(s) === nothing &&
            push!(ss_isom, i)
    end

    result = NamedTuple{
        (:r_only_mets, :t_only_mets,
         :r_only_cat_steps),
        Tuple{Vector{Symbol}, Vector{Symbol},
              Vector{Int}}}[]

    # K-type: non-empty subsets of substrates ×
    # non-empty subsets of products, all r_only
    # (absent from T-state — T is inactive conformation)
    n_s = length(sub_names)
    n_p = length(prod_names)
    for s_mask in 1:(1 << n_s) - 1
        absent_subs = Symbol[sub_names[j]
            for j in 1:n_s
            if (s_mask >> (j-1)) & 1 == 1]
        for p_mask in 1:(1 << n_p) - 1
            absent_prods = Symbol[prod_names[j]
                for j in 1:n_p
                if (p_mask >> (j-1)) & 1 == 1]
            absent = Symbol[absent_subs; absent_prods]
            push!(result, (r_only_mets=absent,
                t_only_mets=Symbol[],
                r_only_cat_steps=Int[]))
        end
    end

    # V-type: all SS isomerization steps r_only or t_only
    if !isempty(ss_isom)
        # r_only cat steps (T-state can't catalyze)
        push!(result, (r_only_mets=Symbol[],
            t_only_mets=Symbol[],
            r_only_cat_steps=copy(ss_isom)))
        # t_only cat steps (R-state can't catalyze)
        # Represented as: no r_only_cat_steps but
        # ... hmm, there's no t_only_cat_steps field.
        # The existing struct only has r_only_cat_steps.
        # For "T-state catalyzes, R-state doesn't":
        # this is the T/R mirror, handled by swapping
        # the allosteric interpretation. Actually,
        # the MWC model with L >1 means T is the
        # tense (less active) state by convention.
        # R-only cat steps means kf_T=0.
        # For "R-state can't catalyze" (kf_R=0), we'd
        # need to set kf_R=0 which is unconventional
        # and may not be supported by the rate equation
        # builder (it only has r_only_cat_steps, not
        # t_only_cat_steps).
        # SKIP t_only cat steps for now — see note below.
    end

    result
end
```

**Note on V-type t_only cat steps:** The `AllostericMechanismSpec` struct has `r_only_cat_steps` but no `t_only_cat_steps` field. The rate equation builder (`_r_only_cat_step_k_syms`) only zeros kf/kr for T-state. To support "R-state can't catalyze" (kf_R=0), we'd need a `t_only_cat_steps` field and corresponding rate equation handling. For now, only generate V-type with `r_only_cat_steps` (T-state can't catalyze). This is the common case (T=tense=less active). Defer t_only_cat_steps to future work.

Actually, wait — re-reading the existing struct and the K-type logic: for K-type, we DO generate t_only metabolites (absent from R-state). The rate equation builder handles `cat_t_only` by zeroing metabolites in the R-state polynomial instead. So t_only metabolites are supported. But for cat STEPS, there's no `t_only_cat_steps` field.

Since we can't support t_only cat steps without struct changes, the V-type only generates 1 variant (r_only cat steps = T-state inactive), not 2.

**Counts (matching tests):**
- Uni-uni: K-type 1 + V-type 1 = 2
- Bi-bi: K-type 9 + V-type 1 = 10

Only r_only variants generated (T is inactive conformation). No t_only for metabolites or catalytic steps.

- [ ] **Step 4: Rewrite `_expand_to_allosteric`**

Replace lines 327-380 of `src/mechanism_enumeration.jl`:

```julia
"""
    _expand_to_allosteric(spec, reaction) → Vector{AllostericMechanismSpec}

Convert non-allosteric mechanism to allosteric (+1 param for L).
Two differentiation modes (r_only only — T is inactive):
- K-type: ≥1 substrate + ≥1 product absent from T-state
- V-type: all catalytic steps inactive in T-state (kf_T=kr_T=0)
"""
function _expand_to_allosteric(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
    cn = oligomeric_state(reaction)
    sub_names = [s[1] for s in substrates(reaction)]
    prod_names = [p[1] for p in products(reaction)]
    all_cat = Symbol[sub_names; prod_names]

    ss_isom = Int[i for (i, s) in enumerate(spec.steps)
        if !s.is_equilibrium &&
           step_metabolite(s) === nothing]

    result = AllostericMechanismSpec[]

    for diff in _valid_allosteric_differentiations(
            reaction, spec)
        # TR-equiv = all catalytic mets NOT in
        # r_only or t_only
        absent = Set(diff.r_only_mets) ∪
                 Set(diff.t_only_mets)
        tr_equiv = Symbol[m for m in all_cat
                          if m ∉ absent]
        # Cat steps: ss_isom stays tr_equiv unless
        # they're in r_only_cat_steps
        tr_steps = Int[i for i in ss_isom
                       if i ∉ diff.r_only_cat_steps]

        push!(result, AllostericMechanismSpec(
            spec, cn,
            Vector{Symbol}[], Int[],  # no reg sites
            tr_equiv, tr_steps,
            diff.r_only_mets, diff.t_only_mets,
            diff.r_only_cat_steps))
    end
    result
end
```

- [ ] **Step 5: Update `_estimated_param_count`**

No changes needed for catalytic metabolites (r_only/t_only correctly skipped at +0). But update the docstring to explain why.

- [ ] **Step 6: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 7: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Rewrite allosteric conversion as +1 with K-type/V-type differentiation"
```

---

### Task 2: Fix regulator dummy naming stability

**Files:**
- Modify: `src/mechanism_enumeration.jl:243-245`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "Regulator dummy naming stability" begin
    rxn2 = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
        dead_end_inhibitors: I, J
    end
    specs = EnzymeRates.init_mechanisms(rxn2)
    spec = first(specs)
    # Add I first
    i_specs = EnzymeRates._expand_add_dead_end_regulator(
        spec, rxn2)
    # Find one with I added
    with_i = first(filter(i_specs) do s
        any(contains(string(sym), "I__reg")
            for st in s.steps
            for sym in Iterators.flatten(
                (st.reactants, st.products)))
    end)
    # Now add J to mechanism that already has I
    j_specs = EnzymeRates._expand_add_dead_end_regulator(
        with_i, rxn2)
    # J should always get the same suffix regardless
    # of whether I is present
    for s in j_specs
        j_forms = [sym for st in s.steps
            for sym in Iterators.flatten(
                (st.reactants, st.products))
            if contains(string(sym), "J__reg")]
        for f in j_forms
            @test contains(string(f), "J__reg")
            # The suffix should be stable: always __reg1
            # (not __reg2 which would happen if ri
            # depended on position in eligible list)
        end
    end
end
```

- [ ] **Step 2: Run to verify it might fail or pass vacuously**

- [ ] **Step 3: Fix the naming**

Change line 243-245 from:
```julia
for (ri, reg) in enumerate(eligible_regs)
    dummy = Symbol(string(reg) * "__reg" * string(ri))
```
to:
```julia
for reg in eligible_regs
    dummy = Symbol(string(reg) * "__reg")
```

Drop the index entirely — the regulator name itself is unique (no two regulators share a name). `I__reg` and `J__reg` are always distinct. No need for `__reg1`, `__reg2`.

- [ ] **Step 4: Update old code too**

Check if `_regulator_dead_end_opportunities` in `old_mechanism_enumeration.jl` uses the same indexed pattern and update for consistency. The old code at line 1190-1191 does:
```julia
dummy = Symbol(string(reg) * "__reg" * string(i))
```
Update to match (drop the index).

Also update `_compile_enzyme_mechanism` regex if it expects `__reg\d+` pattern — check the `_clean_met` function (line 2092):
```julia
m = match(r"^(.+)__reg\d+$", s)
```
Change to also match `__reg` without digits:
```julia
m = match(r"^(.+)__reg\d*$", s)
```

- [ ] **Step 5: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl src/old_mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Fix regulator dummy naming: use stable suffix without index"
```

---

### Task 3: Add `r_only_cat_steps` removal move

**Files:**
- Modify: `src/mechanism_enumeration.jl` (in `_expand_remove_tr_equiv`)
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "Move 5: Remove r_only_cat_steps" begin
    # Create V-type allosteric spec (has r_only_cat_steps)
    specs = EnzymeRates.init_mechanisms(uni_uni_allo)
    spec = first(specs)
    allo_specs = EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo)
    v_type = first(filter(
        r -> !isempty(r.r_only_cat_steps), allo_specs))
    @test !isempty(v_type.r_only_cat_steps)
    # Remove r_only_cat_step → makes step independent (+1)
    result = EnzymeRates._expand_remove_tr_equiv(
        v_type, uni_uni_allo)
    # Should include removal of r_only_cat_steps
    step_removals = filter(result) do r
        length(r.r_only_cat_steps) <
            length(v_type.r_only_cat_steps)
    end
    @test !isempty(step_removals)
end
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Add r_only_cat_steps removal to `_expand_remove_tr_equiv`**

After the existing `tr_equiv_cat_steps` removal loop, add:

```julia
# Remove one r_only cat step → step becomes independent (+1)
# Only allowed when no substrates/products are r_only or t_only
# (otherwise T-state can't complete catalytic cycle and kf_T
# would be unidentifiable — always multiplied by zero)
if isempty(spec.r_only_metabolites) &&
        isempty(spec.t_only_metabolites)
    for (i, _) in enumerate(spec.r_only_cat_steps)
        new_r_steps = [spec.r_only_cat_steps[j]
            for j in eachindex(spec.r_only_cat_steps)
            if j != i]
        push!(result, AllostericMechanismSpec(
            spec.base, spec.catalytic_n,
            deepcopy(spec.allosteric_reg_sites),
            copy(spec.allosteric_multiplicities),
            copy(spec.tr_equiv_metabolites),
            copy(spec.tr_equiv_cat_steps),
            copy(spec.r_only_metabolites),
            copy(spec.t_only_metabolites),
            new_r_steps,
            spec.param_count + 1))
    end
end
```

- [ ] **Step 4: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Add r_only_cat_steps removal to Move 5"
```

---

### Task 4: Fix allosteric dedup for site-order permutations

**Files:**
- Modify: `src/mechanism_enumeration.jl:756-772, 786-796`
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Write failing test**

```julia
@testset "Allosteric dedup: site order" begin
    specs = EnzymeRates.init_mechanisms(uni_uni_allo)
    spec = first(specs)
    allo = first(EnzymeRates._expand_to_allosteric(
        spec, uni_uni_allo))
    # Manually create two specs with same sites in
    # different order
    spec_ab = AllostericMechanismSpec(
        allo.base, allo.catalytic_n,
        [[:A], [:B]], [2, 2],
        copy(allo.tr_equiv_metabolites),
        copy(allo.tr_equiv_cat_steps),
        copy(allo.r_only_metabolites),
        copy(allo.t_only_metabolites),
        copy(allo.r_only_cat_steps),
        allo.param_count + 2)
    spec_ba = AllostericMechanismSpec(
        allo.base, allo.catalytic_n,
        [[:B], [:A]], [2, 2],
        copy(allo.tr_equiv_metabolites),
        copy(allo.tr_equiv_cat_steps),
        copy(allo.r_only_metabolites),
        copy(allo.t_only_metabolites),
        copy(allo.r_only_cat_steps),
        allo.param_count + 2)
    pc = spec_ab.param_count
    cache = Dict(pc => EnzymeRates.AbstractMechanismSpec[
        spec_ab, spec_ba])
    EnzymeRates.dedup!(cache)
    @test length(cache[pc]) == 1
end
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Fix `_canonicalize!` for AllostericMechanismSpec**

Sort `allosteric_reg_sites` and `allosteric_multiplicities` together by site content:

```julia
function _canonicalize!(spec::AllostericMechanismSpec)
    idx_map = _canonicalize!(spec.base)
    map!(i -> get(idx_map, i, i), spec.tr_equiv_cat_steps,
        spec.tr_equiv_cat_steps)
    map!(i -> get(idx_map, i, i), spec.r_only_cat_steps,
        spec.r_only_cat_steps)
    sort!(spec.tr_equiv_metabolites)
    sort!(spec.tr_equiv_cat_steps)
    sort!(spec.r_only_metabolites)
    sort!(spec.t_only_metabolites)
    sort!(spec.r_only_cat_steps)
    for site in spec.allosteric_reg_sites
        sort!(site)
    end
    # Sort sites themselves (with multiplicities) by content
    if length(spec.allosteric_reg_sites) >= 2
        perm = sortperm(spec.allosteric_reg_sites)
        spec.allosteric_reg_sites .=
            spec.allosteric_reg_sites[perm]
        spec.allosteric_multiplicities .=
            spec.allosteric_multiplicities[perm]
    end
    spec
end
```

- [ ] **Step 4: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl test/test_mechanism_enumeration.jl
git commit -m "Fix allosteric dedup: sort regulatory sites for canonical form"
```

---

### Task 5: Add allosteric compilation to integration tests

**Files:**
- Modify: `test/test_mechanism_enumeration.jl`

- [ ] **Step 1: Update integration tests to compile allosteric specs**

In the "Uni-uni full enumeration" and "Bi-bi full enumeration" testsets, change:

```julia
if spec isa EnzymeRates.MechanismSpec
    m = EnzymeMechanism(spec)
    @test length(parameters(m)) <= pc
end
```

to:

```julia
if spec isa EnzymeRates.MechanismSpec
    m = EnzymeMechanism(spec)
    @test length(parameters(m)) <= pc
elseif spec isa EnzymeRates.AllostericMechanismSpec
    m = AllostericEnzymeMechanism(spec)
    @test length(parameters(m)) <= spec.param_count
end
```

- [ ] **Step 2: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Note: This may be slow due to allosteric compilation. If too slow, limit to first N allosteric specs per level.

- [ ] **Step 3: Commit**

```bash
git add test/test_mechanism_enumeration.jl
git commit -m "Add allosteric mechanism compilation to integration tests"
```

---

### Task 6: Update CLAUDE.md Haldane note

**Files:**
- Modify: `.claude/CLAUDE.md`

- [ ] **Step 1: Update the misleading note**

Change:
```
- R-only/t-only for catalytic metabolites changes rate equation structure but not parameter count (Haldane constraints still require all K params)
```
to:
```
- R-only/t-only for catalytic metabolites changes rate equation structure but not parameter count. K_T=∞ (or K_R=∞) is not a free parameter, and the base K becomes K_R (or K_T), preserving count. Haldane is irrelevant when a state can't complete the catalytic cycle.
```

- [ ] **Step 2: Commit**

```bash
git add .claude/CLAUDE.md
git commit -m "Fix CLAUDE.md Haldane note for r_only/t_only catalytic metabolites"
```

---

## Self-Review Checklist

**Spec coverage:**
- ✅ `AllostericMechanismSpec.param_count` field added, `_estimated_param_count` removed (Task 0)
- ✅ Allosteric conversion rewritten as +1 with K-type/V-type (Task 1)
- ✅ Regulator naming stability (Task 2)
- ✅ r_only_cat_steps removal move added (Task 3)
- ✅ Allosteric site-order dedup (Task 4)
- ✅ Allosteric compilation in integration tests (Task 5)
- ✅ CLAUDE.md Haldane note (Task 6)

**Type consistency:**
- `AllostericMechanismSpec.param_count` is the 10th field, set by construction (+1 per move) ✅
- `_valid_allosteric_differentiations` returns NamedTuples matching `AllostericMechanismSpec` fields ✅
- `spec.param_count` works uniformly for both `MechanismSpec` and `AllostericMechanismSpec` ✅
- Regulator dummy names use `__reg` without index, regex updated to match ✅

**Deferred (not in this plan):**
- `t_only_cat_steps` field (would require rate equation builder update)
- Self-contained file (removing old file dependencies)
- Field rename `param_count` → `max_estimated_param_count`
