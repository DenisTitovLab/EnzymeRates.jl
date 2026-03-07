# Mechanism Enumeration Redesign

## Motivation

The current `src/mechanism_enumeration.jl` implementation is brittle and
difficult to test. Testing from `EnzymeReaction` produces too many
mechanisms to verify with hard-coded counts. The new design introduces
staged expansion functions that are independently testable, richer
regulator role modeling, and incremental parameter counting.

---

## Data Types

### `RegulatorRole` hierarchy

```julia
abstract type RegulatorRole end
struct Allosteric <: RegulatorRole end
struct DeadEnd <: RegulatorRole end
struct UnconstrainedRegulator <: RegulatorRole end
```

- `Allosteric()`: OEM expansion + essential activator + general modifier
  (single-regulator special cases)
- `DeadEnd()`: dead-end inhibitor complexes only
- `UnconstrainedRegulator()`: try all roles

### `EnzymeReaction` type parameter change

Regulators type parameter changes from a tuple of Symbols to a tuple of
`(name::Symbol, role::Symbol)` pairs:

```julia
# Old: (:R1, :R2)
# New: ((:R1, :allosteric), (:R2, :dead_end), (:R3, :unknown))
```

### `AbstractMechanismSpec` hierarchy

```julia
abstract type AbstractMechanismSpec end

struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    edges::Vector{Tuple{Int,Int}}
    n_catalytic_edges::Int
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{ParamConstraint}
    param_count::Int
end

struct OligomericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    tr_equivalence::Vector{Bool}
end
```

Stages 1-7 produce `Vector{MechanismSpec}`. Stages 8-10 produce
`Vector{OligomericMechanismSpec}`. Final output includes both.

---

## DSL

```julia
@enzyme_reaction begin
    substrates: A[C], B[N]
    products: P[C], Q[N]
    dead_end_inhibitors: I
    allosteric_regulators: R1
    regulators: R2
end
```

- `dead_end_inhibitors:` -> `DeadEnd()`
- `allosteric_regulators:` -> `Allosteric()`
- `regulators:` -> `UnconstrainedRegulator()` (try all roles)

---

## Pipeline

All stage functions have signature `(specs, reaction) -> Vector{...}`.
Forms, adjacency, etc. are computed internally from the reaction.

### Stage 1: Catalytic Topologies
`_catalytic_topologies(reaction) -> Vector{MechanismSpec}`

DFS for elementary cycles through free enzyme consuming all substrates
and releasing all products. Returns mechanisms with first isomerization
step SS (all other steps RE). `param_count` initialized here.

### Stage 2: RE/SS Assignment
`_expand_ress_variants(specs, reaction) -> Vector{MechanismSpec}`

Appends new variants with different RE/SS positions/counts. Constraints:
G >= 1, at least one SS step, G <= max_re_groups. Original specs from
stage 1 are preserved. Each SS flip: `param_count += 1`.

### Stage 3: General Modifier Expansion
`_expand_general_modifiers(specs, reaction) -> Vector{MechanismSpec}`

For regulators going through non-OEM path (all together): adds parallel
catalytic paths with regulator bound. Both E and ER catalyze. R-binding
edges always RE. Propagates RE/SS from catalytic counterparts. Appends
new variants, preserves originals.

### Stage 4: Essential Activator Expansion
`_expand_essential_activators(specs, reaction) -> Vector{MechanismSpec}`

Same constraint as stage 3 (all non-OEM regulators together). Replaces
catalytic cycle with R-bound version: only ER+S<->ESR<->ER+P, no E+S<->ES
path. E<->ER binding edge added. Appends new variants, preserves originals.

### Stage 5: Dead-End Inhibitor Expansion
`_expand_dead_end_inhibitors(specs, reaction) -> Vector{MechanismSpec}`

Adds dead-end complexes for regulators with `DeadEnd()` or
`UnconstrainedRegulator()` role. Each dead-end RE edge: `param_count += 1`.
Appends new variants, preserves originals.

### Stage 6: Equivalence Constraints
`_expand_equivalence_constraints(specs, reaction) -> Vector{MechanismSpec}`

Groups binding edges by same metabolite AND same site type. Constraints:
- Substrate-site S-binding edges can group together
- Regulator-site R-binding edges for same site can group together
- Substrate-site and regulator-site for same metabolite CANNOT group
- Two different regulatory sites for same metabolite CANNOT group

RE constraint: `param_count -= 1`. SS constraint: `param_count -= 2`.
Appends constrained variants, preserves originals.

### Stage 7: Deduplication
`_deduplicate(specs, reaction) -> Vector{MechanismSpec}`

Removes mechanisms with identical (concentration fingerprint, constraint
descriptor). Keeps the one with fewest parameters.

### Stage 8: Allosteric (OEM) Expansion
`_expand_allosteric(specs, reaction) -> Vector{OligomericMechanismSpec}`

Only for mechanisms NOT expanded in stages 3/4. Site partitioning,
multiplicity combos. `n_conf` is always 2. Appends OEM specs.

### Stage 9: T/R Equivalence
`_expand_tr_equivalence(specs, reaction) -> Vector{OligomericMechanismSpec}`

For each OEM spec, enumerates which parameter groups are T=R equivalent
vs independent. Each equivalence removes parameters. Appends new
variants, preserves originals.

### Stage 10: Post-OEM Deduplication
`_deduplicate_oem(specs, reaction) -> Vector{OligomericMechanismSpec}`

Catches T<->R mirror duplicates (same rate equation with T/R labels
swapped, L -> 1/L). May not be needed; implemented if duplicates found.

### Pipeline Routing

Regulators partition into one of two paths per mechanism:
- **Path A**: {essential_activator, general_modifier, dead_end} (stages 3+4+5)
- **Path B**: {allosteric, dead_end} (stages 5+8+9+10)

Cannot mix stages 3/4 regulators with stage 8 regulators. Both paths
are enumerated for regulators with `UnconstrainedRegulator()` role.
Results merge after stage 7.

### RE/SS Propagation Rule

RE/SS is determined at the catalytic level (stage 2) and propagated:
- Find catalytic counterpart by stripping regulator occupancy from
  both endpoint forms
- Expanded edge gets same RE/SS as its catalytic counterpart
- Regulator binding/release edges (no catalytic counterpart) are
  always RE

### Parameter Count Tracking

| Operation | Delta |
|-----------|-------|
| Base catalytic (all RE) | n_edges - n_thermo_constraints + 2 (E_total, Keq) |
| RE -> SS flip | +1 |
| Dead-end RE edge | +1 |
| Dead-end SS edge | +2 |
| RE equivalence constraint | -1 per constrained edge |
| SS equivalence constraint | -2 per constrained edge |
| General modifier edges | +n new edges (RE or SS) |
| Essential activator | net change from replacing cycle |
| OEM expansion | +L + 1 per regulator-site K |
| T/R equivalence | removes params per constrained group |

Verified by sampling ~10 compiled mechanisms per stage, calling
`parameters()`, comparing to `param_count`.

---

## `enumerate_mechanisms` API

```julia
enumerate_mechanisms(
    reaction::EnzymeReaction;
    max_re_groups::Int = 7,
    catalytic_n::Int = 2,
) -> Vector{AbstractMechanismSpec}
```

Regulator roles come from `EnzymeReaction` type parameters.

---

## Testing

### File structure

```
src/old_mechanism_enumeration.jl          # renamed from current
test/old_mechanism_definitions_...jl      # renamed from current
test/old_test_mechanism_enum_...jl        # renamed from current
test/old_reaction_definitions_...jl       # renamed from current

src/mechanism_enumeration.jl              # new implementation
test/mechanism_enumeration_test_specs.jl  # new test specs
test/test_mechanism_enumeration.jl        # new tests
```

### `StageExpansionTestSpec`

Tests individual stages starting from one hand-built `MechanismSpec`,
running through ALL stages, comparing output count at each stage to
hard-coded values.

```julia
struct StageExpansionTestSpec
    name::String
    reaction::Any
    base_mechanism::MechanismSpec

    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int
    expected_n_tr_equiv::Int
    expected_n_oem_dedup::Int
end
```

### `EnumerationTestSpec`

End-to-end from `EnzymeReaction`, comparing output count at each stage
to hard-coded values.

```julia
struct EnumerationTestSpec
    name::String
    reaction::Any

    expected_n_catalytic::Int
    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int
    expected_n_tr_equiv::Int
    expected_n_oem_dedup::Int
end
```

### Shared reaction set (both spec types)

1. Uni-Uni (no regulators)
2. Uni-Uni + 1 unconstrained regulator
3. Uni-Bi + 1 regulator
4. Bi-Bi + 2 regulators
5. Bi-Bi Ping-Pong (no regulators)
6. Bi-Bi Ping-Pong + 1 regulator

### Parameter count verification

For each stage output, sample ~10 mechanisms (deterministic seed),
compile them, call `parameters()`, verify `param_count` matches.

### Orthogonal verification helpers

Test-local combinatorial formulas to cross-check counts (not calling
`enumerate_mechanisms` internals).

### Large mechanism safety

For Bi-Bi + 2 regulators: add timeout/skip flags for stages that may
be too large without lazy implementation.

---

## Regulator Mechanics

### General Modifier (single allosteric regulator, L->small limit)
Both catalytic paths exist:
- E + S <-> ES <-> E + P
- ER + S <-> ESR <-> ER + P
- Plus R-binding: E <-> ER, ES <-> ESR (always RE)

Can act as activator or inhibitor depending on relative rates.

### Essential Activator (single allosteric regulator, L->large limit)
Only R-bound path catalyzes:
- ER + S <-> ESR <-> ER + P
- E <-> ER binding edge (always RE)
- E + S <-> ES does NOT exist

### Full Allosteric (OEM, any number of regulators)
`OligomericEnzymeMechanism` with NConf=2 (T and R states).
Catalytic mechanism duplicated for each conformation, with T/R
equivalence variants for each parameter group.

### Dead-End Inhibitor
R binds enzyme forms creating dead-end complexes (no catalytic
activity). Dead-end edges inherit RE/SS from catalytic counterpart.
