# Allosteric Mechanism Redesign Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace `OligomericEnzymeMechanism` with `AllostericEnzymeMechanism`, unify GM/EA/MWC into a single MWC model, implement T/R equivalence expansion, and restructure the enumeration pipeline from 10 stages to 8.

**Architecture:** Add `AbstractEnzymeMechanism` supertype, rename `OligomericEnzymeMechanism` → `AllostericEnzymeMechanism` (removing `NConf` — always 2 conformations), delete general modifier and essential activator stages from the enumeration pipeline, restructure the orchestrator into 2 independent phases, and implement T/R equivalence expansion (currently a no-op stub).

**Tech Stack:** Julia, `@generated` functions for compile-time rate equation derivation, King-Altman/Cha method for symbolic algebra.

**Key reference files:**
- Design doc: `docs/plans/2026-03-09-allosteric-mechanism-redesign.md`
- Previous design: `docs/plans/2026-03-07-mechanism-enumeration-redesign.md`

---

## Task 1: Add AbstractEnzymeMechanism and rename type in types.jl

**Files:**
- Modify: `src/types.jl:290-303` (type definition), `src/types.jl:538-565` (accessors)

**Step 1: Add AbstractEnzymeMechanism and rename struct**

In `src/types.jl`, before the `EnzymeMechanism` constructor (around line 220), add:

```julia
abstract type AbstractEnzymeMechanism end
```

Make `EnzymeMechanism` a subtype — find the struct definition and add `<: AbstractEnzymeMechanism`.

Replace the `OligomericEnzymeMechanism` definition (lines 290-303) with:

```julia
"""
    AllostericEnzymeMechanism{Metabolites, CatalyticMech, CatalyticN, RegSites}

Singleton type for allosteric enzymes (MWC model, always 2 conformations).

- `Metabolites`: tuple of `Symbol` names from `metabolites:` block
- `CatalyticMech`: `EnzymeMechanism` type for one catalytic subunit
- `CatalyticN`: number of catalytic sites per enzyme molecule
- `RegSites`: tuple of `((ligand_syms...,), multiplicity)` pairs
"""
struct AllostericEnzymeMechanism{
    Metabolites, CatalyticMech, CatalyticN, RegSites,
} <: AbstractEnzymeMechanism end
```

Note: `NConf` type parameter is removed entirely (always 2 conformations).

**Step 2: Update accessor methods**

Replace all accessor methods (lines 538-565). Change every
`OligomericEnzymeMechanism{M,CM,N,RS,NC}` to
`AllostericEnzymeMechanism{M,CM,N,RS}` and remove the `NC` where
clause. Section header comment: change "OligomericEnzymeMechanism"
to "AllostericEnzymeMechanism".

```julia
# ─── AllostericEnzymeMechanism Accessors ────────────────────────

"""Delegate structural accessors to the CatalyticMech singleton."""
n_states(::AllostericEnzymeMechanism{M,CM,N,RS}) where {M,CM,N,RS} =
    n_states(CM())
n_steps(::AllostericEnzymeMechanism{M,CM,N,RS}) where {M,CM,N,RS} =
    n_steps(CM())
equilibrium_steps(::AllostericEnzymeMechanism{M,CM,N,RS}) where {M,CM,N,RS} =
    equilibrium_steps(CM())
substrates(::AllostericEnzymeMechanism{M,CM,N,RS}) where {M,CM,N,RS} =
    substrates(CM())
products(::AllostericEnzymeMechanism{M,CM,N,RS}) where {M,CM,N,RS} =
    products(CM())
@generated function regulators(
    ::AllostericEnzymeMechanism{M,CM,N,RS},
) where {M,CM,N,RS}
    ligs = Symbol[]
    for (ligands, _) in RS
        for lig in ligands
            lig in ligs || push!(ligs, lig)
        end
    end
    Tuple(ligs)
end
param_constraints(::AllostericEnzymeMechanism) = ()

"""Return all metabolite names (catalytic + regulatory) from the Metabolites type param."""
metabolites(::AllostericEnzymeMechanism{Mets,CM,N,RS}) where {Mets,CM,N,RS} = Mets
```

**Step 3: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: Many failures (downstream code still references old name). This is expected — we fix them in subsequent tasks.

**Step 4: Commit**

```bash
git add src/types.jl
git commit -m "Add AbstractEnzymeMechanism, rename OligomericEnzymeMechanism → AllostericEnzymeMechanism"
```

---

## Task 2: Update exports and module file

**Files:**
- Modify: `src/EnzymeRates.jl:4`

**Step 1: Update export**

Change line 4 from:
```julia
export EnzymeReaction, EnzymeMechanism, OligomericEnzymeMechanism
```
to:
```julia
export EnzymeReaction, EnzymeMechanism, AllostericEnzymeMechanism
```

**Step 2: Commit**

```bash
git add src/EnzymeRates.jl
git commit -m "Update export: OligomericEnzymeMechanism → AllostericEnzymeMechanism"
```

---

## Task 3: Update DSL (@enzyme_mechanism)

**Files:**
- Modify: `src/dsl.jl:331-419`

**Step 1: Update `_parse_oligomeric_mechanism`**

Rename function to `_parse_allosteric_mechanism`. Remove `nconf`
variable entirely. Update the return expression to use
`AllostericEnzymeMechanism` without `NConf`:

Changes to make:
1. Line 331: docstring — replace "OligomericEnzymeMechanism" with "AllostericEnzymeMechanism"
2. Line 335: rename function `_parse_oligomeric_mechanism` → `_parse_allosteric_mechanism`
3. Line 337: delete `nconf = 1`
4. Lines 363-364: delete the `elseif label == :conformations` branch entirely (remove 2 lines)
5. Lines 414-419: update the return expression:

```julia
    :(let _cm = $cm_expr
        AllostericEnzymeMechanism{$mets_tuple, typeof(_cm), $catalytic_n, $reg_sites_expr}()
    end)
```

**Step 2: Update call site**

Find where `_parse_oligomeric_mechanism` is called (around line 218)
and rename to `_parse_allosteric_mechanism`.

**Step 3: Decide on `conformations:` DSL keyword**

The `conformations:` keyword in the DSL is no longer needed (always 2).
Two options:
- Remove it entirely (breaking change for existing mechanism definitions)
- Keep it but ignore the value with a deprecation warning

Since CLAUDE.md says "get Denis's explicit approval before implementing
ANY backward compatibility", remove it entirely. Any existing
`conformations: 2` in test mechanism definitions will need updating
(Task 8).

**Step 4: Commit**

```bash
git add src/dsl.jl
git commit -m "Update DSL: rename to _parse_allosteric_mechanism, remove NConf/conformations"
```

---

## Task 4: Update rate equation derivation — rename and remove NConf

**Files:**
- Modify: `src/rate_eq_derivation.jl` (many locations)
- Modify: `src/sym_poly_for_rate_eq_derivation.jl:555-619`

This is the largest single task. It involves:
1. Replacing `_AnyMechanism` with `AbstractEnzymeMechanism`
2. Renaming all `OligomericEnzymeMechanism` → `AllostericEnzymeMechanism`
3. Removing `NConf` from all type parameter lists and function signatures
4. Removing all `NConf == 1` code paths (always 2 conformations)
5. Renaming `_build_oligomeric_rate_body` → `_build_allosteric_rate_body` (and similar helpers)

**Step 1: Replace `_AnyMechanism`**

Line 3: Replace
```julia
const _AnyMechanism = Union{EnzymeMechanism, OligomericEnzymeMechanism}
```
with:
```julia
const _AnyMechanism = AbstractEnzymeMechanism
```

This preserves all existing dispatch signatures that use `_AnyMechanism`
while changing the underlying type to the new abstract supertype.

**Step 2: Rename all `OligomericEnzymeMechanism` references**

Global search-replace in `src/rate_eq_derivation.jl`:
- `OligomericEnzymeMechanism` → `AllostericEnzymeMechanism`
- Remove `,NConf` from all `where {Mets,CM,CatN,RS,NConf}` clauses
- Remove `NConf` from all function argument lists

Key locations (approximate line numbers — verify after earlier edits):
- Line 824+: `_kcat_forward` dispatch
- Lines 1025-1039: `_binding_K_symbols`
- Lines 1052-1082: `_dependent_param_exprs`
- Lines 1111-1134: `_oligomeric_dep_assignments`
- Lines 1141-1190: `_oligomeric_num_den_exprs`
- Lines 1193-1211: `_build_oligomeric_rate_body`
- Lines 1215-1220: `rate_equation` dispatch
- Lines 1224-1247: `rate_equation_string`
- Lines 1251-1259: `structural_identifiability_deficit`

**Step 3: Remove NConf == 1 code paths**

In each function that branches on `NConf`:

`_kcat_forward` (lines 839-969):
- Remove the `if NConf == 1` early return branch
- Keep only the `NConf == 2` logic (which becomes the only path)

`_binding_K_symbols` (lines 1025-1039):
- Remove `NConf == 2 ?` ternary — always include T-state K's
- Remove `NConf == 2 ?` ternary for reg K's — always include T-state reg K's

`_dependent_param_exprs` for AllostericEnzymeMechanism (lines 1052-1082):
- Remove `if NConf == 1` early return
- Keep only the NConf==2 path (which becomes unconditional)

`_oligomeric_dep_assignments` (lines 1111-1134):
- Remove `NConf` parameter from function signature
- Remove `if NConf == 2` conditional — always generate T assignments

`_oligomeric_num_den_exprs` (lines 1141-1190):
- Remove `NConf` parameter from function signature
- Remove `NConf == 1 && return ...` early return
- Always generate T-state expressions

`_build_oligomeric_rate_body` (lines 1193-1211):
- Rename to `_build_allosteric_rate_body`
- Remove `NConf` parameter
- Update call to `_oligomeric_dep_assignments` (no NConf arg)
- Update call to `_oligomeric_num_den_exprs` (no NConf arg)
- Update `M_type` construction (no NConf)

`rate_equation` dispatch (lines 1215-1220):
- Update type signature (no NConf)
- Call `_build_allosteric_rate_body` instead of `_build_oligomeric_rate_body`

`rate_equation_string` (lines 1224-1247):
- Update type signature (no NConf)
- Update `_oligomeric_dep_assignments` call (no NConf)
- Update `_oligomeric_num_den_exprs` call (no NConf)

`structural_identifiability_deficit` (lines 1251-1259):
- Update type signature (no NConf)
- Update `_count_oligomeric_rate_monomials` call (no NConf)

**Step 4: Update sym_poly helper**

In `src/sym_poly_for_rate_eq_derivation.jl`, update
`_count_oligomeric_rate_monomials` (lines 570-619):
- Remove `NConf` parameter
- Remove `if NConf == 2` conditional — always include T-state
- Rename to `_count_allosteric_rate_monomials`

**Step 5: Run tests**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: Failures in tests that still reference old names or use
`conformations: 2` in DSL. Mechanism enumeration tests may also fail
due to `OligomericMechanismSpec` references.

**Step 6: Commit**

```bash
git add src/rate_eq_derivation.jl src/sym_poly_for_rate_eq_derivation.jl
git commit -m "Update rate equation derivation: rename to AllostericEnzymeMechanism, remove NConf"
```

---

## Task 5: Update mechanism enumeration — rename spec types

**Files:**
- Modify: `src/mechanism_enumeration.jl:64-75` (spec type), `:1175-1254` (allosteric stages), `:1260-1320` (compile_mechanism), `:1327-1420` (orchestrator)

**Step 1: Rename OligomericMechanismSpec → AllostericMechanismSpec**

Replace struct definition (lines 64-75):
```julia
"""
    AllostericMechanismSpec <: AbstractMechanismSpec

Allosteric enzyme mechanism specification (MWC model, always 2
conformations). Wraps a base MechanismSpec with allosteric annotation.
"""
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    tr_equivalence::Vector{Bool}
end
```

**Step 2: Update `_expand_allosteric`**

Replace function (lines 1175-1204). Change return type and constructor
calls from `OligomericMechanismSpec` to `AllostericMechanismSpec`:

```julia
function _expand_allosteric(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    catalytic_n::Int=2,
    allosteric_regs::Vector{Symbol}=Symbol[],
)
    result = AllostericMechanismSpec[]
    if isempty(allosteric_regs)
        for spec in specs
            push!(result, AllostericMechanismSpec(
                spec, catalytic_n,
                Vector{Symbol}[], Int[], Bool[]))
        end
        return result
    end
    partitions = _set_partitions(allosteric_regs)
    for spec in specs
        n_re_binding = count(spec.equilibrium_steps)
        for partition in partitions
            n_groups = length(partition)
            for combo in Iterators.product(
                    ntuple(_ -> 1:catalytic_n, n_groups)...)
                push!(result, AllostericMechanismSpec(
                    spec, catalytic_n, partition,
                    collect(combo),
                    fill(false, n_re_binding)))
            end
        end
    end
    result
end
```

Note: `tr_equivalence` is initialized to `fill(false, n_re_binding)` —
all T/R parameters independent by default. `_expand_tr_equivalence`
(Task 6) will enumerate the 2^n variants.

**Step 3: Update `_expand_tr_equivalence` signature**

Change parameter and return types from `OligomericMechanismSpec` to
`AllostericMechanismSpec`. Keep the passthrough behavior for now
(Task 6 implements it):

```julia
function _expand_tr_equivalence(
    specs::Vector{AllostericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    specs  # passthrough until Task 6
end
```

**Step 4: Rename `_deduplicate_oem` → `_deduplicate_allosteric`**

Update function (lines 1231-1254). Change name, parameter types,
and key function:

```julia
function _deduplicate_allosteric(
    specs::Vector{AllostericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    seen = Dict{Any, AllostericMechanismSpec}()
    for spec in specs
        key = _allosteric_canonical_key(spec)
        if !haskey(seen, key)
            seen[key] = spec
        end
    end
    collect(values(seen))
end

function _allosteric_canonical_key(spec::AllostericMechanismSpec)
    base_key = (spec.base.edges, spec.base.equilibrium_steps,
                spec.base.param_constraints)
    tr = spec.tr_equivalence
    (base_key, sort(spec.allosteric_reg_sites),
     spec.allosteric_multiplicities, tr)
end
```

**Step 5: Update `compile_mechanism` for AllostericMechanismSpec**

Replace the `OligomericMechanismSpec` dispatch (lines 1301-1317):

```julia
function compile_mechanism(spec::AllostericMechanismSpec)
    cm = compile_mechanism(spec.base)

    # Build metabolites tuple: catalytic metabolites + regulator symbols
    cat_mets = metabolites(cm)
    reg_syms = Symbol[]
    for site in spec.allosteric_reg_sites
        for s in site
            s in reg_syms || s in cat_mets || push!(reg_syms, s)
        end
    end
    mets = (cat_mets..., reg_syms...)

    # Build RegSites type parameter
    reg_sites = Tuple(
        (Tuple(site), spec.allosteric_multiplicities[i])
        for (i, site) in enumerate(spec.allosteric_reg_sites))

    AllostericEnzymeMechanism{mets, typeof(cm), spec.catalytic_n, reg_sites}()
end
```

Note: No `NConf` parameter — removed from type. The `2` that was
hardcoded as the last type parameter is no longer needed.

**Step 6: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "Rename OligomericMechanismSpec → AllostericMechanismSpec, update pipeline functions"
```

---

## Task 6: Delete GM/EA stages and restructure orchestrator

**Files:**
- Modify: `src/mechanism_enumeration.jl:788-967` (delete GM/EA), `:1327-1420` (restructure orchestrator)

**Step 1: Delete `_expand_general_modifiers`**

Delete the entire function (lines 788-878, approximately 90 lines).

**Step 2: Delete `_expand_essential_activators`**

Delete the entire function (lines 889-967, approximately 80 lines).

**Step 3: Restructure `enumerate_mechanisms` orchestrator**

Replace the orchestrator (lines 1327-1420) with the new 2-phase
structure. Key changes:
- Remove stages 3 and 4 (GM/EA) from the inner loop
- Restructure into Phase 1 (base pipeline) and Phase 2 (allosteric)
- Phase 2 runs independently with its own dedup
- Always emit allosteric variants when allosteric regs present
  (regardless of catalytic_n — use cn=1 when catalytic_n=0)

```julia
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    stage::EnumerationStage=FullEnumeration(),
    max_re_groups::Int=7,
    catalytic_n::Int=0,
)
    # Stage 1: Catalytic topologies
    catalytic = _catalytic_topologies(reaction)
    stage isa Catalytic && return catalytic

    # Regulator partitioning
    roles = regulator_roles(reaction)
    fixed_dead_end = Symbol[r[1] for r in roles if r[2] == :dead_end]
    fixed_allosteric = Symbol[r[1] for r in roles
                              if r[2] == :allosteric]
    unknown = Symbol[r[1] for r in roles if r[2] == :unknown]
    n_unknown = length(unknown)

    if stage isa WithDeadEnd
        all_de = MechanismSpec[]
        for reg_mask in 0:(1 << n_unknown) - 1
            de_regs = Symbol[fixed_dead_end;
                [unknown[i] for i in 1:n_unknown
                 if (reg_mask >> (i - 1)) & 1 == 0]]
            append!(all_de, _expand_dead_end_inhibitors(
                catalytic, reaction; dead_end_regs=de_regs))
        end
        return all_de
    end

    all_base = MechanismSpec[]
    all_allosteric = AllostericMechanismSpec[]

    for reg_mask in 0:(1 << n_unknown) - 1
        de_regs = Symbol[fixed_dead_end;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 0]]
        allo_regs = Symbol[fixed_allosteric;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 1]]

        # Phase 1: base mechanism pipeline (chained)
        base = _expand_ress_variants(
            catalytic, reaction; max_re_groups)
        base = _expand_dead_end_inhibitors(
            base, reaction; dead_end_regs=de_regs)
        base = _expand_equivalence_constraints(base, reaction)
        base = _deduplicate(base, reaction)
        append!(all_base, base)

        # Phase 2: allosteric expansion (independent)
        if !isempty(allo_regs)
            cn = catalytic_n > 0 ? catalytic_n : 1
            allo = _expand_allosteric(base, reaction;
                catalytic_n=cn, allosteric_regs=allo_regs)
            allo = _expand_tr_equivalence(allo, reaction)
            allo = _deduplicate_allosteric(allo, reaction)
            append!(all_allosteric, allo)
        end
    end

    total = length(all_base) + length(all_allosteric)
    inner = Iterators.flatten((all_base, all_allosteric))
    MechanismIterator(inner, total)
end
```

**Step 4: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "Delete GM/EA stages, restructure orchestrator into 2-phase pipeline"
```

---

## Task 7: Implement `_expand_tr_equivalence`

**Files:**
- Modify: `src/mechanism_enumeration.jl` (the `_expand_tr_equivalence` function)

**Step 1: Write the failing test**

In `test/test_mechanism_enumeration.jl`, add a focused test for TR
equivalence (exact location will depend on test restructuring in Task
8, but the logic is):

```julia
@testset "TR equivalence expansion" begin
    # A base mechanism with 2 RE binding edges should produce
    # 2^2 = 4 TR equivalence variants
    base = EnzymeRates._catalytic_topologies(uni_uni)[1]
    allo = EnzymeRates._expand_allosteric(
        [base], uni_uni_allosteric_I;
        catalytic_n=1, allosteric_regs=[:I])
    @test length(allo) > 0
    n_before = length(allo)
    expanded = EnzymeRates._expand_tr_equivalence(
        allo, uni_uni_allosteric_I)
    n_re = count(base.equilibrium_steps)
    @test length(expanded) == n_before * 2^n_re
    # Verify TR equivalence vectors are all distinct combos
    tr_combos = Set(s.tr_equivalence for s in expanded)
    @test length(tr_combos) == 2^n_re
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: FAIL because `_expand_tr_equivalence` is still a passthrough.

**Step 3: Implement `_expand_tr_equivalence`**

Replace the stub with:

```julia
"""
    _expand_tr_equivalence(specs, reaction)
        -> Vector{AllostericMechanismSpec}

Enumerate T/R parameter equivalence masks. For each RE binding
parameter, K_S_T can equal K_S_R (equivalent, fewer params) or be
independent (more params). Produces 2^n variants per input spec.
"""
function _expand_tr_equivalence(
    specs::Vector{AllostericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = AllostericMechanismSpec[]
    for spec in specs
        n = length(spec.tr_equivalence)
        for mask in 0:(1 << n) - 1
            tr_eq = [((mask >> (i-1)) & 1) == 1 for i in 1:n]
            push!(result, AllostericMechanismSpec(
                spec.base, spec.catalytic_n,
                spec.allosteric_reg_sites,
                spec.allosteric_multiplicities, tr_eq))
        end
    end
    result
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: The TR equivalence test passes. Other tests may still fail
(test spec restructuring in Task 8).

**Step 5: Commit**

```bash
git add src/mechanism_enumeration.jl
git commit -m "Implement _expand_tr_equivalence: enumerate T/R parameter equivalence masks"
```

---

## Task 8: Update test infrastructure

**Files:**
- Modify: `test/mechanism_enumeration_test_specs.jl` (spec types, builders, helper)
- Modify: `test/test_mechanism_enumeration.jl` (test assertions)
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl` (OEM type refs)
- Modify: `test/test_enzyme_derivation.jl` (OEM dispatch)

This task updates ALL tests to work with the renamed types and
restructured pipeline. Expected counts will change because GM/EA
stages are removed and TR equivalence is now implemented.

**Step 1: Update StageExpansionTestSpec**

In `test/mechanism_enumeration_test_specs.jl`, update the struct
(lines 99-116):

```julia
Base.@kwdef struct StageExpansionTestSpec
    name::String
    reaction::Any
    base_mechanism::EnzymeRates.MechanismSpec
    dead_end_regs::Vector{Symbol} = Symbol[]
    allosteric_regs::Vector{Symbol} = Symbol[]
    catalytic_n::Int = 0

    expected_n_ress::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int = 0
    expected_n_tr_equiv::Int = 0
    expected_n_allosteric_dedup::Int = 0
end
```

Remove: `expected_n_general_modifier`, `expected_n_essential_activator`,
`expected_n_oem_dedup`.
Add: `expected_n_allosteric_dedup` (replacing `expected_n_oem_dedup`).

**Step 2: Update EnumerationTestSpec**

Similar changes (lines 123-140):

```julia
Base.@kwdef struct EnumerationTestSpec
    name::String
    reaction::Any
    catalytic_n::Int = 0

    expected_n_forms::Int
    expected_n_catalytic::Int
    expected_n_ress::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int = 0
    expected_n_tr_equiv::Int = 0
    expected_n_allosteric_dedup::Int = 0
    expected_n_total::Int
end
```

Remove: `expected_n_general_modifier`, `expected_n_essential_activator`,
`expected_n_oem_dedup`.

**Step 3: Update spec builders**

In `build_no_reg_stage_expansion_specs()`,
`build_single_reg_stage_expansion_specs()`, and
`build_multi_reg_stage_expansion_specs()`:

- Remove all `expected_n_general_modifier=...` and
  `expected_n_essential_activator=...` lines from every spec
  instantiation
- Rename `expected_n_oem_dedup` → `expected_n_allosteric_dedup`
- **Recalculate expected values**: Since GM/EA stages are removed,
  the pipeline produces fewer mechanisms. The expected counts for
  stages that remain (ress, dead_end, equivalence, dedup) are
  unchanged. Allosteric counts change because TR equivalence now
  expands.

For allosteric specs, the new expected counts need to account for
TR equivalence. For a base mechanism with `n_re` RE binding edges
and 1 allosteric regulator:
- `expected_n_allosteric` = (set partitions) × (multiplicities) =
  same as before (1 partition × catalytic_n multiplicities)
- `expected_n_tr_equiv` = `expected_n_allosteric` × 2^n_re
  (each allosteric spec gets 2^n_re TR equivalence variants)
- `expected_n_allosteric_dedup` ≤ `expected_n_tr_equiv`
  (dedup may reduce)

**Important**: The exact expected counts must be determined
empirically by running the pipeline in the Julia REPL. Do NOT
guess — compute them by running:

```julia
using EnzymeRates
rxn = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    allosteric_regulators: I
end
base = EnzymeRates._catalytic_topologies(rxn)
allo = EnzymeRates._expand_allosteric(base, rxn;
    catalytic_n=1, allosteric_regs=[:I])
tr = EnzymeRates._expand_tr_equivalence(allo, rxn)
dd = EnzymeRates._deduplicate_allosteric(tr, rxn)
println("allosteric=$(length(allo)), tr=$(length(tr)), dedup=$(length(dd))")
```

**Step 4: Update `_run_full_pipeline_stages` helper**

Remove calls to `_expand_general_modifiers` and
`_expand_essential_activators`. Remove corresponding accumulator
variables (`n_gm`, `n_ea`). Update return NamedTuple to remove
`general_modifier` and `essential_activator` keys. Rename
`oem_dedup` → `allosteric_dedup`.

Update the allosteric section to use the new 2-phase structure:
- Always call `_expand_allosteric` when allosteric regs present
  (not just when `catalytic_n > 0`)
- Use `cn = catalytic_n > 0 ? catalytic_n : 1`

**Step 5: Update `build_enumeration_specs`**

Remove `expected_n_general_modifier` and
`expected_n_essential_activator` from all specs. Add
`expected_n_allosteric`, `expected_n_tr_equiv`,
`expected_n_allosteric_dedup` where applicable. Recalculate
`expected_n_total` (will differ from current values since GM/EA
variants are no longer produced, but allosteric+TR variants are
added for reactions with allosteric regs).

**Again**: determine exact values empirically in the REPL.

**Step 6: Update test_mechanism_enumeration.jl**

Update the stage expansion test loop (lines 31-69):
- Remove `_expand_general_modifiers` and `_expand_essential_activators`
  assertions
- Keep ress, dead_end, equivalence, dedup assertions
- Update allosteric section to not require `catalytic_n > 0`
  (allosteric expansion now runs for all allosteric regs)

Update end-to-end tests (lines 72-107):
- Remove `general_modifier` and `essential_activator` count assertions
- Add allosteric count assertions for reactions with allosteric regs

Update property tests (lines 110-225):
- Remove stage monotonicity assertions for GM/EA
- Update remaining monotonicity checks (dead_end >= ress, etc.)

Update OEM property tests (lines 320-377):
- Replace `OligomericMechanismSpec` with
  `EnzymeRates.AllostericMechanismSpec`
- Replace `OligomericEnzymeMechanism` with `AllostericEnzymeMechanism`

Update combinatorial cross-checks (lines 380-444):
- Remove assertions for `expected_n_general_modifier` and
  `expected_n_essential_activator`
- Add assertions for TR equivalence counts
  (e.g., `expected_n_tr_equiv == expected_n_allosteric * 2^n_re`)

**Step 7: Update mechanism_definitions_for_test_enzyme_derivation.jl**

Replace all `isdefined(EnzymeRates, :OligomericEnzymeMechanism)` guards
(lines 1472, 1656, 1917) with
`isdefined(EnzymeRates, :AllostericEnzymeMechanism)`.

Remove `conformations: 2` lines from all `@enzyme_mechanism` blocks
within these guards (the DSL no longer accepts this keyword).

**Step 8: Update test_enzyme_derivation.jl**

Replace `OligomericEnzymeMechanism` with `AllostericEnzymeMechanism`
in the helper function (line 223-226):

```julia
function compute_all_params(m::EnzymeRates.AllostericEnzymeMechanism)
    ...
end
```

**Step 9: Run full test suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL tests pass.

**Step 10: Commit**

```bash
git add test/
git commit -m "Update tests for allosteric mechanism redesign"
```

---

## Task 9: Clean up dead code and verify

**Files:**
- Modify: `src/mechanism_enumeration.jl` (remove any remaining dead code)
- Check: `src/rate_eq_derivation.jl` (verify no orphaned helpers)

**Step 1: Search for remaining old references**

Run in the REPL or via grep:
```bash
grep -rn "OligomericEnzymeMechanism\|OligomericMechanismSpec\|_expand_general_modifiers\|_expand_essential_activators\|_deduplicate_oem\|_oem_canonical_key\|NConf\|oligomeric" src/ test/
```

Fix any remaining references.

**Step 2: Check for dead helper functions**

After deleting GM/EA, check if any helper functions in
`mechanism_enumeration.jl` are now unused:
- `_find_dead_end` — still used by `_expand_dead_end_inhibitors`
- `reg_bound_map` construction helpers — if they were only used by
  GM/EA, they can be deleted

**Step 3: Update CLAUDE.md**

Update the architecture documentation in `.claude/CLAUDE.md`:
- Replace `OligomericEnzymeMechanism` references with
  `AllostericEnzymeMechanism`
- Remove NConf references
- Update pipeline stage list (10 → 8 stages)
- Remove GM/EA stage descriptions
- Add note about T/R equivalence implementation
- Update `_AnyMechanism` description

**Step 4: Run full test suite one more time**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: ALL tests pass with no warnings.

**Step 5: Commit**

```bash
git add -A  # after checking git status
git commit -m "Clean up dead code and update documentation for allosteric redesign"
```

---

## Task Summary

| Task | Description | Est. Complexity |
|------|-------------|-----------------|
| 1 | Add AbstractEnzymeMechanism, rename type in types.jl | Small |
| 2 | Update exports | Trivial |
| 3 | Update DSL | Small |
| 4 | Update rate equation derivation (rename + remove NConf) | Large |
| 5 | Rename enumeration spec types | Medium |
| 6 | Delete GM/EA, restructure orchestrator | Medium |
| 7 | Implement TR equivalence | Small |
| 8 | Update all tests | Large |
| 9 | Clean up and verify | Small |

## Dependency Graph

```
Task 1 (types) → Task 2 (exports) → Task 3 (DSL)
                                         ↓
                                    Task 4 (rate eq)
                                         ↓
                                    Task 5 (enum specs)
                                         ↓
                                    Task 6 (pipeline restructure)
                                         ↓
                                    Task 7 (TR equiv)
                                         ↓
                                    Task 8 (tests)
                                         ↓
                                    Task 9 (cleanup)
```

Tasks must be executed sequentially — each depends on the previous.

## Important Notes

- **Expected count determination**: Task 8 requires computing new
  expected counts empirically in the Julia REPL. Do NOT guess these
  values. Run the actual pipeline functions and record the outputs.
- **NConf removal**: The `conformations:` DSL keyword is removed, not
  deprecated. Any existing mechanism definitions using it must be
  updated.
- **TR equivalence**: This is genuinely new functionality (was a no-op
  stub). It needs careful testing — verify that the 2^n expansion
  produces the right number of variants and that dedup correctly
  identifies duplicates.
- **Backward compatibility**: `OligomericEnzymeMechanism` is an
  exported name. Removing it is a breaking change. Since this is a
  pre-1.0 package, this is acceptable.
