# Mechanism Enumeration Test Redesign

Restructure mechanism enumeration tests to separate test data from test
logic, use DSL macros for mechanism definitions, add combinatorial
verification of hardcoded counts, and align with the staged pipeline
design.

---

## File Structure

```
test/mechanism_enumeration_test_specs.jl  # spec types, helper, reactions, builders
test/test_mechanism_enumeration.jl        # test code only (iterates over specs)
```

`runtests.jl` includes the specs file before the test file.

---

## Helper Function

```julia
mechanism_spec_from_mechanism(mechanism::EnzymeMechanism, reaction) -> MechanismSpec
```

~15-line function that extracts edges, equilibrium steps, param constraints,
and param count from a compiled `EnzymeMechanism`'s type parameters and
returns a `MechanismSpec`. Uses `enumerate_enzyme_forms` to map species
tuples back to form indices.

---

## Spec Types

### `StageExpansionTestSpec`

Each stage runs **independently** on the same `base_mechanism` — not a
pipeline where output feeds forward. This makes each expected count
independently derivable from the base mechanism's properties.

Regulator info (dead-end vs allosteric) is derived from the reaction's
`regulator_roles()` in the test code, not stored as fields. Stage
functions continue to accept `dead_end_regs`/`allosteric_regs` kwargs
(needed because `enumerate_mechanisms` does 2^n partitioning of
`:unknown` regs and stages must know which partition they're in).

```julia
Base.@kwdef struct StageExpansionTestSpec
    name::String
    reaction::Any
    base_mechanism::EnzymeRates.MechanismSpec
    catalytic_n::Int = 0

    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int = 0
    expected_n_tr_equiv::Int = 0
    expected_n_oem_dedup::Int = 0
end
```

### `EnumerationTestSpec`

Full pipeline end-to-end, including regulator partitioning.

```julia
Base.@kwdef struct EnumerationTestSpec
    name::String
    reaction::Any
    catalytic_n::Int = 0

    expected_n_forms::Int
    expected_n_catalytic::Int
    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int = 0
    expected_n_tr_equiv::Int = 0
    expected_n_oem_dedup::Int = 0
    expected_n_total::Int
end
```

---

## Builder Pattern

Follows `mechanism_definitions_for_test_enzyme_derivation.jl` pattern:
mechanisms defined inline using `@enzyme_mechanism` macro inside `let`
blocks within builder functions.

```julia
function build_stage_expansion_specs()
    specs = StageExpansionTestSpec[]

    let  # Uni-Uni, dead-end I: catalytic-only base
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
        end
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products: P[C]
                enzymes: E, ES[C], EP[C]
            end
            steps: begin
                [E, S] <--> [ES]
                [E, P] <--> [EP]
                [ES] --> [EP]
            end
        end
        base = mechanism_spec_from_mechanism(m, rxn)
        push!(specs, StageExpansionTestSpec(
            name="Uni-Uni dead-end I, catalytic base",
            reaction=rxn,
            base_mechanism=base,
            # 3 edges: S-bind(RE), P-bind(RE), isom(SS)
            # RE/SS: 2 RE edges, 2^2-1=3 flip combos + original...
            expected_n_ress=...,
            ...
        ))
    end

    specs
end

function build_enumeration_specs()
    specs = EnumerationTestSpec[]
    # ... similar pattern ...
    specs
end
```

---

## Reaction Set (8 reactions)

Regulated reactions have two versions: `:unknown` for
`EnumerationTestSpec` (tests partitioning), explicit roles for
`StageExpansionTestSpec` (single partition, easier to reason about).

| # | Reaction | Regulators | Notes |
|---|----------|------------|-------|
| 1 | Uni-Uni | none | Simplest case, 3 forms |
| 2 | Uni-Uni + 1 reg | `:unknown` / `:dead_end` / `:allosteric` | Tests partitioning |
| 3 | Uni-Bi + 1 reg | `:unknown` / `:dead_end` / `:allosteric` | Multi-product |
| 4 | Uni-Bi + allosteric reg | `:allosteric`, catalytic_n=2 | OEM expansion |
| 5 | Bi-Bi | none | 9 topologies |
| 6 | Bi-Bi Ping-Pong | none | Distinct topology |
| 7 | Bi-Bi Ping-Pong + 1 reg | `:unknown` / `:dead_end` / `:allosteric` | Regs on ping-pong |
| 8 | Bi-Bi + 2 regs | `:unknown` / 1`:dead_end`+1`:allosteric` | Multi-regulator |

---

## StageExpansionTestSpec Instances

Multiple specs per reaction with different base mechanisms targeting
different stage behaviors.

### No-reg reactions (1, 5, 6): 1 spec each
- Catalytic-only base. Interesting: RE/SS, equivalence (for Bi-Bi
  where same metabolite binds multiple forms), dedup.

### Single-reg reactions (2, 3, 7): 3 specs each
- **Catalytic-only base, reg as `:dead_end`** — tests DE expansion,
  RE/SS propagation to dead-end edges. GM/EA pass through.
- **Catalytic-only base, reg as `:allosteric`** — tests GM and EA
  expansion. DE passes through.
- **Pre-expanded base with regulator edges, reg as `:dead_end`** —
  tests equivalence between catalytic and dead-end binding edges,
  dedup on richer graph. DE passes through (edges already present).

### Uni-Bi + allosteric, catalytic_n=2 (4): 2 specs
- **Catalytic-only base** — tests GM, EA expansion + OEM stages.
- **General-modifier base** (parallel paths present) — tests
  equivalence between E+S→ES and ER+S→ESR edges, OEM on
  pre-expanded graph.

### Bi-Bi + 2 regs (8): 3 specs
- **Catalytic-only base** — both DE and allosteric stages expand.
- **Base with dead-end edges** — allosteric expansion + equivalence
  across dead-end and catalytic edges.
- **Base with GM + dead-end edges** — equivalence and dedup on
  complex graph.

---

## EnumerationTestSpec Instances

One spec per reaction, all using `:unknown` regs where applicable.
Bi-Bi + 2 regs gets a timeout wrapper:

```julia
@testset "End-to-end: Bi-Bi + 2 regs" begin
    result = @timeout 120 begin
        collect(enumerate_mechanisms(bi_bi_two_regs_unknown))
    end
    if result === nothing
        @warn "Bi-Bi + 2 regs timed out — lazy enumeration needed"
        @test_broken false
    else
        @test length(result) == s.expected_n_total
    end
end
```

---

## Test Code Structure (`test_mechanism_enumeration.jl`)

```julia
@testset "Mechanism Enumeration Pipeline" begin

    # Stage expansion: each stage independently on base
    @testset "Stage expansion: $(s.name)" for s in STAGE_EXPANSION_SPECS
        rxn = s.reaction
        de_regs = [r[1] for r in regulator_roles(rxn)
                   if r[2] == :dead_end]
        al_regs = [r[1] for r in regulator_roles(rxn)
                   if r[2] == :allosteric]
        base = [s.base_mechanism]

        @test length(_expand_ress_variants(base, rxn)) ==
            s.expected_n_ress
        @test length(_expand_general_modifiers(base, rxn;
            allosteric_regs=al_regs)) ==
            s.expected_n_general_modifier
        @test length(_expand_essential_activators(base, rxn;
            allosteric_regs=al_regs)) ==
            s.expected_n_essential_activator
        @test length(_expand_dead_end_inhibitors(base, rxn;
            dead_end_regs=de_regs)) ==
            s.expected_n_dead_end
        @test length(_expand_equivalence_constraints(base, rxn)) ==
            s.expected_n_equivalence
        @test length(_deduplicate(base, rxn)) ==
            s.expected_n_dedup

        if s.catalytic_n > 0
            # OEM stages...
        end
    end

    # End-to-end pipeline
    @testset "End-to-end: $(s.name)" for s in ENUMERATION_SPECS
        # ... full pipeline with _run_full_pipeline_stages ...
        # ... + enumerate_mechanisms for total ...
    end

    # Property-based tests (kept from current)
    # param_count accuracy (kept from current)
    # compile_mechanism round-trip (kept from current)
    # Combinatorial cross-checks
end
```

---

## Orthogonal Verification

### Layer 1: Comments on expected values

Each expected count has an inline comment explaining the combinatorial
derivation:

```julia
# 3 edges: S-bind(RE), P-bind(RE), isom(SS).
# RE/SS flips: 2 RE edges → 2^2-1=3 non-empty subsets to flip.
# But must keep ≥1 RE group → exclude all-SS = 3 valid flips.
# Total = 1 (original) + 3 (flips) - 1 (all-SS) = 3
expected_n_ress=3,
```

### Layer 2: Executable cross-checks

A testset with independent formulas computing expected counts from
reaction properties (form count, edge count, regulator count) without
calling pipeline internals:

```julia
@testset "Combinatorial cross-checks" begin
    # Uni-Uni RE/SS: C(n_re, k) summed over valid k
    n_re = 2  # RE edges in Uni-Uni triangle
    expected = sum(binomial(n_re, k) for k in 1:n_re-1) + 1
    @test STAGE_EXPANSION_SPECS[1].expected_n_ress == expected
    # ...
end
```

---

## Preserved from Current Tests

- Property-based tests (monotonicity, dedup reduces count, etc.)
- `param_count` accuracy tests (sampled compilation + parameter count
  comparison)
- `compile_mechanism` round-trip tests
- OEM expansion property tests
