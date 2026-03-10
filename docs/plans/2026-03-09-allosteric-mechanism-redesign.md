# Allosteric Mechanism Redesign

## Motivation

The current mechanism enumeration pipeline has 10 stages, three of which
(general modifier, essential activator, allosteric/MWC expansion) perform
related graph manipulations to handle regulators that affect catalysis.
These stages are complex, produce large intermediate mechanism counts,
and share significant duplicated code.

This redesign unifies all three into a single allosteric model based on
MWC (Monod-Wyman-Changeux) theory. General modifier and essential
activator are limiting cases of MWC that emerge during fitting (L→0 or
L→∞), not separate enumeration variants.

### Key insight: GM and EA are MWC limits

For MWC with CatN=1, one regulator R binding only the T-state
(K_R_R = ∞):

```
v = (N_R + L·N_T·(1 + R/K_R_T)) / (Q_R + L·Q_T·(1 + R/K_R_T))
```

When L << 1 (R-state dominates), the constant "1" inside (1+R/K_R_T)
contributes L·Q_T which is negligible compared to Q_R. So:

```
L·Q_T·(1 + R/K_R_T) ≈ L·Q_T·R/K_R_T = Q_T·R/(K_R_T/L) = Q_T·R/K'_R
```

Setting K'_R = K_R_T/L absorbs the allosteric constant into an apparent
regulator binding constant. The rate equation becomes
indistinguishable from a general modifier:

```
v ≈ (N_R + N_T·R/K'_R) / (Q_R + Q_T·R/K'_R)
```

Similarly, L→∞ gives essential activator behavior (T-state dominates,
only R-bound state catalyzes effectively). Model selection (cross-
validation) with different parameter counts handles the choice between
these limiting cases.

## Goals

1. Replace `OligomericEnzymeMechanism` with `AllostericEnzymeMechanism`
   (always MWC, always 2 conformations)
2. Eliminate general modifier and essential activator as separate
   enumeration stages — MWC subsumes both
3. Implement T/R equivalence expansion (currently a no-op stub) to
   enumerate parameter-equivalence variants that reduce param count
4. Simplify the enumeration pipeline from 10 stages to 8 (5 base + 3
   allosteric)
5. Add `AbstractEnzymeMechanism` supertype for extensibility

## Type System

### Runtime types (src/types.jl)

```julia
abstract type AbstractEnzymeMechanism end

# Unchanged — base catalytic mechanism + dead-end inhibitors
struct EnzymeMechanism{Species, Reactions, EqSteps, PC} <: AbstractEnzymeMechanism
end

# Replaces OligomericEnzymeMechanism
# Always MWC with 2 conformations (T + R)
struct AllostericEnzymeMechanism{Mets, BaseMech<:EnzymeMechanism,
                                 CatN, RegSites} <: AbstractEnzymeMechanism
end
```

Type parameters:
- `Mets`: metabolite symbols tuple (for dispatch)
- `BaseMech`: the base `EnzymeMechanism` type (catalytic cycle +
  dead-end inhibitors)
- `CatN`: number of catalytic subunits (1 for monomeric, ≥2 for
  oligomeric)
- `RegSites`: tuple of `(ligands::Tuple{Symbol...}, multiplicity::Int)`
  pairs encoding which regulators bind at each allosteric site and
  their Hill-like multiplicity

`NConf` is removed as a type parameter — always 2 conformations.

`_AnyMechanism` Union is replaced by `AbstractEnzymeMechanism` dispatch.

### Enumeration spec types (src/mechanism_enumeration.jl)

```julia
abstract type AbstractMechanismSpec end

# Unchanged — base catalytic + dead-end graph
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    edges::Vector{Tuple{Int,Int}}
    n_catalytic_edges::Int
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{ParamConstraint}
    param_count::Int
end

# Replaces OligomericMechanismSpec
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    tr_equivalence::Vector{Bool}   # per RE binding param: true = K_T == K_R
end
```

The `tr_equivalence` vector has one entry per RE binding parameter in
the base mechanism. `true` means K_S_T = K_S_R (same binding affinity
in both conformations — fewer parameters). `false` means they are
independent (more parameters). `_expand_tr_equivalence` enumerates all
2^n combinations.

## Enumeration Pipeline

### Overview

```
Phase 1 — base mechanism enumeration (chained, 5 stages):
  1. _catalytic_topologies(rxn)
  2. _expand_ress_variants(specs, rxn)
  3. _expand_dead_end_inhibitors(specs, rxn; dead_end_regs)
  4. _expand_equivalence_constraints(specs, rxn)
  5. _deduplicate(specs, rxn)
  → base_specs::Vector{MechanismSpec}  (deduplicated)

Phase 2 — allosteric expansion (independent, 3 stages):
  6. _expand_allosteric(base_specs; catalytic_n, allosteric_regs)
  7. _expand_tr_equivalence(allosteric_specs, rxn)
  8. _deduplicate_allosteric(allosteric_specs, rxn)
  → allosteric_specs::Vector{AllostericMechanismSpec}  (deduplicated)

Output: AbstractMechanismSpec[base_specs; allosteric_specs]
```

Phase 1 is unchanged from the current pipeline (stages 1, 2, 5, 6, 7
renumbered). Phase 2 replaces current stages 3, 4, 8, 9, 10.

Each phase produces independently valid, deduplicated mechanism specs.
The final output is their concatenation.

### Phase 1 stages (unchanged logic)

**Stage 1: `_catalytic_topologies(rxn)`**
Find valid catalytic cycles via DFS. Each topology gets one SS
isomerization edge, all others RE. Output: `Vector{MechanismSpec}`.

**Stage 2: `_expand_ress_variants(specs, rxn)`**
Enumerate RE/SS assignments for catalytic edges. Dead-end edges
inherit RE/SS from catalytic counterparts. Filters by RE group count
(2 ≤ G ≤ max_re_groups).

**Stage 3: `_expand_dead_end_inhibitors(specs, rxn; dead_end_regs)`**
Add dead-end inhibitor complexes. For each catalytic form, regulator
can bind creating dead-end state. Enumerates all subsets of form ×
regulator binding. Dead-end edges inherit RE/SS from catalytic
counterparts.

**Stage 4: `_expand_equivalence_constraints(specs, rxn)`**
Enumerate parameter constraint masks for equivalent binding edges
(e.g., K_S for E+S⇌ES equals K_S for ER+S⇌ESR in dead-end complex).

**Stage 5: `_deduplicate(specs, rxn)`**
Deduplicate by (concentration fingerprint, constraint descriptor).
Keeps mechanism with fewest parameters.

### Phase 2 stages

**Stage 6: `_expand_allosteric(base_specs; catalytic_n, allosteric_regs)`**

For each base spec, enumerate allosteric variants:
- Partition allosteric regulators into sites via `_set_partitions`
  (Bell number enumeration — regulators at the same site share a
  denominator factor `(1 + R1/K_R1 + R2/K_R2)^m`)
- For each partition, enumerate multiplicity combos
  (`1:catalytic_n` per site)
- All `tr_equivalence` entries initialized to `false` (all T/R
  parameters independent — maximum parameter count)

```julia
function _expand_allosteric(base_specs, rxn;
                            catalytic_n, allosteric_regs)
    result = AllostericMechanismSpec[]
    isempty(allosteric_regs) && return result
    partitions = _set_partitions(allosteric_regs)
    for spec in base_specs
        n_re_binding = count(spec.equilibrium_steps)  # RE binding params
        for partition in partitions
            n_groups = length(partition)
            for combo in Iterators.product(
                    ntuple(_ -> 1:catalytic_n, n_groups)...)
                push!(result, AllostericMechanismSpec(
                    spec, catalytic_n, partition, collect(combo),
                    fill(false, n_re_binding)))
            end
        end
    end
    result
end
```

**Stage 7: `_expand_tr_equivalence(allosteric_specs, rxn)`**

For each `AllostericMechanismSpec`, enumerate all 2^n combinations of
T/R parameter equivalence. Each RE binding parameter can be either
independent (K_S_T ≠ K_S_R, more params) or equivalent (K_S_T = K_S_R,
fewer params).

```julia
function _expand_tr_equivalence(specs, rxn)
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

This replaces the current no-op stub. The 2^n expansion is manageable
because it operates on allosteric specs (already filtered by
partition/multiplicity), not on full mechanism graphs.

**Stage 8: `_deduplicate_allosteric(allosteric_specs, rxn)`**

Deduplicate on `(base_fingerprint, sorted_reg_sites, multiplicities,
tr_equivalence)`. The base fingerprint comes from the already-deduped
phase 1 output, so this mainly removes duplicates from T/R equivalence
expansion (e.g., different masks that produce the same effective
parameter set).

### Orchestrator

```julia
function enumerate_mechanisms(rxn; catalytic_n=0)
    roles = regulator_roles(rxn)
    fixed_de = [r[1] for r in roles if r[2] == :dead_end]
    fixed_allo = [r[1] for r in roles if r[2] == :allosteric]
    unknown = [r[1] for r in roles if r[2] == :unknown]
    n_unknown = length(unknown)

    all_results = AbstractMechanismSpec[]

    for reg_mask in 0:(1 << n_unknown) - 1
        de_regs = [fixed_de;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i-1)) & 1 == 0]]
        allo_regs = [fixed_allo;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i-1)) & 1 == 1]]

        # Phase 1: base mechanism pipeline
        base = _catalytic_topologies(rxn)
        base = _expand_ress_variants(base, rxn)
        base = _expand_dead_end_inhibitors(base, rxn;
            dead_end_regs=de_regs)
        base = _expand_equivalence_constraints(base, rxn)
        base = _deduplicate(base, rxn)
        append!(all_results, base)

        # Phase 2: allosteric expansion
        if !isempty(allo_regs)
            cn = catalytic_n > 0 ? catalytic_n : 1
            allo = _expand_allosteric(base, rxn;
                catalytic_n=cn, allosteric_regs=allo_regs)
            allo = _expand_tr_equivalence(allo, rxn)
            allo = _deduplicate_allosteric(allo, rxn)
            append!(all_results, allo)
        end
    end

    all_results
end
```

When `catalytic_n=0` (not explicitly oligomeric), allosteric expansion
uses `cn=1` — monomeric MWC. This subsumes general modifier and
essential activator as fitting limits (L→0 and L→∞ respectively).

## Rate Equation Derivation

### AllostericEnzymeMechanism rate equation

One function `_build_allosteric_rate_body` replaces
`_build_oligomeric_rate_body`. Always 2 conformations (T + R):

```
N_R, Q_R = _raw_symbolic_rate_polys(BaseMech)
N_T, Q_T = rename_params(N_R, Q_R, "_T")    # T-state parameters
reg_Q_R = [(1 + R_i/K_Ri_R)^m_i for each reg site]
reg_Q_T = [(1 + R_i/K_Ri_T)^m_i for each reg site]

v = E_total * CatN *
    (N_R · Q_R^(CatN-1) · ∏reg_R + L · N_T · Q_T^(CatN-1) · ∏reg_T) /
    (Q_R^CatN · ∏reg_R + L · Q_T^CatN · ∏reg_T)
```

Key properties:
- Structurally identical to current `_build_oligomeric_rate_body` with
  NConf=2
- No NConf branching — the NConf=1 code path is removed
- T/R equivalence constraints (K_S_T = K_S_R) reduce independent
  parameter count, handled by substituting T-params with R-params
  where `tr_equivalence[i] == true`

### T/R equivalence in rate equation

When `tr_equivalence[i] == true` for a binding parameter K_i:
- K_i_T is not an independent parameter
- It is set equal to K_i_R in the rate equation
- This reduces the fitted parameter count by 1 per equivalent pair
- In the generated rate equation body, K_i_T is replaced by K_i_R

### Parameter count

For `AllostericMechanismSpec`:
```
param_count = base.param_count              # base mechanism params (R-state)
            + n_T_independent_params        # T-state params (not TR-equivalent)
            + 1                             # L (allosteric constant)
            + n_reg_K_params                # K_Ri_R and K_Ri_T per reg site
```

Where:
- `n_T_independent_params` = number of base params where
  `tr_equivalence == false` (T-state has own value)
- `n_reg_K_params` = number of regulator binding constants (R-state
  + T-state for each regulator, minus TR-equivalent ones)

## compile_mechanism

`compile_mechanism` dispatches on spec type:

- `MechanismSpec` → `EnzymeMechanism` (unchanged)
- `AllostericMechanismSpec` → `AllostericEnzymeMechanism`

For `AllostericMechanismSpec`:
1. Compile `spec.base` → `EnzymeMechanism` (the base catalytic
   mechanism becomes `BaseMech` type parameter)
2. Construct `AllostericEnzymeMechanism{Mets, BaseMech, CatN, RegSites}`
   where:
   - `Mets` = metabolites from the base mechanism + regulator symbols
   - `CatN` = `spec.catalytic_n`
   - `RegSites` = tuple encoding reg sites, multiplicities, and T/R
     equivalence info

## What Gets Deleted

| Deleted | Reason |
|---------|--------|
| `OligomericEnzymeMechanism` | Replaced by `AllostericEnzymeMechanism` |
| `OligomericMechanismSpec` | Replaced by `AllostericMechanismSpec` |
| `_expand_general_modifiers` (old stage 3) | GM absorbed into MWC |
| `_expand_essential_activators` (old stage 4) | EA absorbed into MWC |
| `_deduplicate_oem` (old stage 10) | Replaced by `_deduplicate_allosteric` |
| `NConf` type parameter | Always 2 conformations |
| `_build_oligomeric_rate_body` | Replaced by `_build_allosteric_rate_body` |
| NConf=1 code paths in rate derivation | Removed |
| `_AnyMechanism` Union type | Replaced by `AbstractEnzymeMechanism` |

## What Gets Added

| Added | Purpose |
|-------|---------|
| `AbstractEnzymeMechanism` | Supertype for dispatch |
| `AllostericEnzymeMechanism` | Unified allosteric type |
| `AllostericMechanismSpec` | Enumeration spec for allosteric variants |
| `_build_allosteric_rate_body` | Rate equation (always 2 conformations) |
| `_deduplicate_allosteric` | Dedup for allosteric specs |
| `_expand_tr_equivalence` impl | Was a no-op stub, now implemented |

## What Gets Modified

| Modified | Change |
|----------|--------|
| `EnzymeMechanism` | Adds `<: AbstractEnzymeMechanism` |
| `@enzyme_mechanism` DSL | Constructs `AllostericEnzymeMechanism` instead of `OligomericEnzymeMechanism` |
| `parameters`, `metabolites`, etc. | Dispatch on `AllostericEnzymeMechanism` instead of `OligomericEnzymeMechanism` |
| `rescale_parameter_values` | Updated for new type |
| `_kcat_forward` | Updated for new type |
| `rate_equation` | Dispatches to `_build_allosteric_rate_body` |
| `rate_equation_string` | Updated for new type |
| `structural_identifiability_deficit` | Updated for new type |

## Impact on Tests

- `StageExpansionTestSpec`: remove `expected_n_general_modifier` and
  `expected_n_essential_activator` fields. Add
  `expected_n_tr_equivalence` field. Renumber stages.
- `EnumerationTestSpec`: same field changes. Expected counts will
  change (fewer total mechanisms since GM/EA variants are removed,
  but TR equivalence variants are added).
- `mechanism_definitions_for_test_enzyme_derivation.jl`:
  `OligomericEnzymeMechanism` → `AllostericEnzymeMechanism` throughout.
- `test_enzyme_derivation.jl`: type name changes; rate equation tests
  verify same mathematical structure.
- `@enzyme_mechanism` DSL tests: updated for new type name.

## Testing Strategy

Phase 1 stages (1-5) are tested as before — `StageExpansionTestSpec`
with base mechanisms through each stage independently.

Phase 2 stages (6-8) are tested via `StageExpansionTestSpec` with
allosteric regulator specs:
- Stage 6 (`_expand_allosteric`): verify set partition × multiplicity
  counts
- Stage 7 (`_expand_tr_equivalence`): verify 2^n expansion for each
  allosteric spec
- Stage 8 (`_deduplicate_allosteric`): verify dedup reduces count

End-to-end testing via `EnumerationTestSpec` verifies the complete
pipeline including regulator partitioning.

Parameter count accuracy is verified by sampling compiled mechanisms
and checking `spec.param_count == length(parameters(m))`.
