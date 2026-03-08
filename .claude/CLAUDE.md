# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

When saving new lessons learned or project knowledge, write them to this file (`.claude/CLAUDE.md`) so they are shared across all machines via git. Do not use the local auto memory for project-specific knowledge.

## Package Goal

EnzymeRates.jl identifies the best enzyme rate equation from kinetic data. Given a reaction definition and experimental rate measurements at varying concentrations, the package enumerates all biochemically valid mechanisms, fits each to data, and selects the one with fewest parameters that adequately describes the data (cross-validation). See `SPEC.md` for the full API specification.

**Primary use case**: `EnzymeReaction` + data → `IdentifyRateEquationProblem` → `identify_rate_equation()` → `IdentifyRateEquationResults`

**Secondary use cases**: manually define mechanisms via `@enzyme_mechanism` and derive/fit rate equations.

## API Design (see SPEC.md)

- **19 exported symbols** (planned): 6 types, 2 macros, 2 constants (`Full`, `Reduced`), 9 functions. Currently 16 — `IdentifyRateEquationProblem`, `IdentifyRateEquationResults`, `identify_rate_equation` are pending implementation.
- `compile_mechanism` is exported: converts `MechanismSpec` → `EnzymeMechanism` or `OligomericEnzymeMechanism`
- Enumeration internals (`SiteState`, `EnzymeFormSpec`, `MechanismSpec`, `enumerate_mechanisms`, etc.) are NOT part of the public API — accessible via `IdentifyRateEquationProblem` fields for power users
- Data tables use a `group` column (not `Article`+`Fig`) to identify measurement groups sharing the same E_total
- Cross-validation: leave-one-group-out
- Keq is always user-provided, never estimated from data

## Commands

```bash
# Run full test suite (cold — pays precompilation + JIT cost every time)
julia --project -e 'using Pkg; Pkg.test()'
```

## Workflow

- Always run tests before committing

## Code Style

- 92-character line length limit, 4-space indentation
- Prefer minimal code: inline single-use helpers, avoid unnecessary abstractions
- Remove unused features entirely — don't add parameters to disable them
- After any refactor, re-read changed files for dead code and further simplification

## Key Architecture Decisions

- `EnzymeMechanism{Species,Reactions,EquilibriumSteps,ParamConstraints}` is a singleton type encoding mechanism info in type parameters
- `OligomericEnzymeMechanism{Mets,CatalyticMech,CatN,RegSites,NConf}` represents multi-subunit MWC allosteric enzymes — see `src/types.jl` and `src/dsl.jl` for DSL syntax
- `EnzymeReaction{S,P,R}` similarly encodes reactions in types
- Each unique mechanism = unique type → affects compilation time
- `@generated` functions used for compile-time computation (metabolites, graph, stoich_matrix, rate_equation, _kcat_forward)
- `_AnyMechanism = Union{EnzymeMechanism, OligomericEnzymeMechanism}` used for shared dispatch (e.g., `rescale_parameter_values`)

### Canonical Step Form
- The `EnzymeMechanism` constructor normalizes RE steps so metabolite is always on LHS (binding direction): `[E, S] ⇌ [ES]`, never `[ES] ⇌ [E, S]`
- SS steps are NOT canonicalized (swapping kf↔kr would break analytical test formulas)
- After canonicalization, all RE metabolite K params are binding Kd (displayed as `1/K`). Non-binding RE steps (pure isomerization) retain Ka convention.
- `_binding_K_symbols` relies on this invariant: checks only for metabolite on LHS, no RHS check needed

### Regulator representation
- Regulators are `(name::Symbol, role::Symbol)` pairs in `EnzymeReaction` type parameter `R`
- `RegulatorRole` hierarchy: `Allosteric`, `DeadEnd`, `UnconstrainedRegulator` (abstract: `RegulatorRole`)
- `regulators(m)` returns a tuple of bare `Symbol` names for backward compatibility
- `regulator_roles(rxn)` returns raw `(name, role)` pairs from `EnzymeReaction`
- `@enzyme_reaction` DSL accepts `regulators:` (role=`:unknown`), `dead_end_inhibitors:` (`:dead_end`), `allosteric_regulators:` (`:allosteric`)
- Plain `Symbol` regulators auto-normalize to `(name, :unknown)` in the `EnzymeReaction` constructor
- Substrates/products are `(name, atoms)` tuples — access name via `s[1]`

### Dead-end SS/RE propagation
- Dead-end substrate/product-binding edges (e.g., ER↔ESR for S-binding) inherit RE/SS status from their catalytic counterpart (E↔ES), found by stripping regulator sites from endpoints
- Dead-end regulator-binding edges (E↔ER, ES↔ESR) remain always RE
- `_dead_end_catalytic_map` maps each dead-end edge to its unique catalytic counterpart; `_propagate_de_eq_steps!` copies the RE/SS status

### Dead-end parameter equivalence constraints
- `_find_equivalent_groups` includes dead-end edges with catalytic counterparts: when ER↔ESR binds the same non-product metabolite S as catalytic E↔ES, they join the same equiv group
- This adds constrained variants where K_S_dead_end = K_S_catalytic (fewer params: R doesn't affect S-binding)
- Dead-end edges always have lower index than catalytic, so `g[1]` is always the catalytic edge

### Mechanism enumeration staged pipeline
- `MechanismSpec` has 6 fields: `reaction, edges, n_catalytic_edges, equilibrium_steps, param_constraints, param_count`
- `OligomericMechanismSpec` has 5 fields: `base::MechanismSpec, catalytic_n, allosteric_reg_sites, allosteric_multiplicities, tr_equivalence`
- Pipeline order: `_catalytic_topologies` (stage 1) → `_expand_ress_variants` (stage 2) → `_expand_general_modifiers` (stage 3) → `_expand_essential_activators` (stage 4) → `_expand_dead_end_inhibitors` (stage 5) → `_expand_equivalence_constraints` (stage 6) → `_deduplicate` (stage 7) → `_expand_allosteric` (stage 8) → `_expand_tr_equivalence` (stage 9, currently passthrough) → `_deduplicate_oem` (stage 10)
- Stages 1-7 produce `Vector{MechanismSpec}`, stages 8-10 produce `Vector{OligomericMechanismSpec}`
- Regulator partitioning (2^n_unknown masks for unknown-role regs) happens in `enumerate_mechanisms` orchestration; stage functions take explicit `dead_end_regs`/`allosteric_regs` kwargs
- `_set_partitions` enumerates all Bell-number set partitions of allosteric regulators
- `compile_mechanism` converts `MechanismSpec` → `EnzymeMechanism` or `OligomericMechanismSpec` → `OligomericEnzymeMechanism`
- Same-site regulators share a `(1 + R1/K_R1 + R2/K_R2)^m` denominator factor

### General modifier and essential activator
- `_expand_general_modifiers`: duplicates catalytic cycle with regulator bound (parallel paths), R-binding edges always RE
- `_expand_essential_activators`: replaces catalytic cycle with R-bound version (only ER+S→ESR path), adds E→ER binding edge

## Source Layout

- `src/types.jl` — `EnzymeReaction`, `EnzymeMechanism`, `OligomericEnzymeMechanism` structs; `RegulatorRole` hierarchy (`Allosteric`, `DeadEnd`, `UnconstrainedRegulator`); `EnzymeMechanism` and `OligomericEnzymeMechanism` accessors; `regulator_roles()`; `RateEquationMode` hierarchy
- `src/dsl.jl` — `@enzyme_reaction` (supports `substrates:`, `products:`, `regulators:`, `dead_end_inhibitors:`, `allosteric_regulators:` labels) and `@enzyme_mechanism` macros (handles both `EnzymeMechanism` and `OligomericEnzymeMechanism` DSL)
- `src/sym_poly_for_rate_eq_derivation.jl` — Symbolic polynomial algebra (`Poly` type); `_rename_poly_T`, `_count_oligomeric_rate_monomials` for MWC identifiability
- `src/rate_eq_derivation.jl` — King-Altman/Cha rate equation derivation via `@generated` functions; parameters API; identifiability checks; kcat computation (`_is_ss_rate_constant`, `_kcat_components`, `_kcat_forward`) and `rescale_parameter_values`; OligomericEnzymeMechanism MWC rate equation assembly (`_build_oligomeric_rate_body`, `rate_equation_string`, `structural_identifiability_deficit`)
- `src/thermodynamic_constr_for_rate_eq_derivation.jl` — Haldane/Wegscheider thermodynamic constraints, dependent parameter elimination, `_build_rate_body` for `@generated rate_equation`
- `src/fitting.jl` — `FittingProblem`, `loss!`, `fit_rate_equation` using Optimization.jl
- `src/mechanism_enumeration.jl` — Staged pipeline for mechanism enumeration. Types: `SiteDefinition`, `EnzymeFormSpec`, `MechanismSpec` (6 fields), `OligomericMechanismSpec`, `MechanismIterator`. 10 stages: `_catalytic_topologies` → `_expand_ress_variants` → `_expand_general_modifiers` → `_expand_essential_activators` → `_expand_dead_end_inhibitors` → `_expand_equivalence_constraints` → `_deduplicate` → `_expand_allosteric` → `_expand_tr_equivalence` → `_deduplicate_oem`. `compile_mechanism` converts specs to `EnzymeMechanism`/`OligomericEnzymeMechanism`. Helpers: `_dead_end_catalytic_map`/`_propagate_de_eq_steps!`, `_set_partitions`/`_partition_mult_count`, `_concentration_fingerprint`/`_constraint_descriptor` for dedup

## Vmax Normalization (kcat factoring) — IMPLEMENTED

### Implementation
- `_kcat_forward(m, params)`: `@generated` function computing kcat analytically from polynomial structure
- `_kcat_components(M)`: extracts (num_k, den_k) candidate pairs by grouping polynomials by metabolite pattern
- `_is_ss_rate_constant(sym)`: classifies symbols as SS rate constants (lowercase `k` followed by digit)
- `rescale_parameter_values(m, params; kcat=1.0)`: public API, scales SS k's uniformly so kcat = target

### Key properties
- kcat is homogeneous degree-1 in SS k's, independent of RE K's
- Uniform k-degree in denominator guarantees v/(E_total * kcat) is scale-invariant
- For mechanisms with multiple catalytic paths (e.g., non-essential activator), kcat = max over all paths
- For OligomericEnzymeMechanism with NConf=2, kcat depends on regulator corner; returns max over 2^n_lig corners

## Known Issues

### `rate_equation` compilation limits for large mechanisms
- `rate_equation(m, conc, params)` uses `@generated` functions that derive the rate equation at compile time via King-Altman/Cha method
- For mechanisms with many enzyme forms/steps, compilation can be extremely slow, exhaust memory, or StackOverflow
- This is inherent to the type-parameter-based architecture: each unique `EnzymeMechanism` type triggers full symbolic derivation at compile time
- Workaround in tests: only the simplest mechanisms (first 10 by form count) are tested with `rate_equation`; larger mechanisms are tested only for enumeration correctness
- Future fix: `identify_rate_equation` should order candidates by `param_count_estimate` (ascending) and skip mechanisms that exceed a time/memory budget

## Testing

- Tests include Aqua (quality) and JET (static analysis)
- Don't leave profiling deps (SnoopCompile) in Project.toml — Aqua stale deps check will fail
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` defines shared mechanisms used by multiple test files — must be included before those tests
- Mechanism enumeration tests use two data-driven spec types: `StageExpansionTestSpec` (hand-built `MechanismSpec` through individual stages) and `EnumerationTestSpec` (end-to-end from `EnzymeReaction` through full pipeline). Both verify expected counts at each stage.
- Parameter count verification: per-stage sampling compiles mechanisms and checks `param_count == length(parameters(m))`
- `MechanismTestSpec` has optional `analytical_kcat_fn` field for per-mechanism kcat formula validation
- kcat/rescaling tests (scale invariance, rate proportionality, V≈1, custom target) run for ALL mechanism specs in the main `run_all_tests` loop — not in a separate file
