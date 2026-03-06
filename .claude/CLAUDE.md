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
- Regulators are plain `Symbol`s (no atom compositions) — they don't participate in chemical transformation
- `regulators(m)` returns a tuple of bare `Symbol`s for both `EnzymeMechanism` and `OligomericEnzymeMechanism`
- Substrates/products are `(name, atoms)` tuples — access name via `s[1]`

## Source Layout

- `src/types.jl` — `EnzymeReaction`, `EnzymeMechanism`, `OligomericEnzymeMechanism` structs; `EnzymeMechanism` and `OligomericEnzymeMechanism` accessors; `RateEquationMode` hierarchy
- `src/dsl.jl` — `@enzyme_reaction` and `@enzyme_mechanism` macros (handles both `EnzymeMechanism` and `OligomericEnzymeMechanism` DSL)
- `src/sym_poly_for_rate_eq_derivation.jl` — Symbolic polynomial algebra (`Poly` type); `_rename_poly_T`, `_count_oligomeric_rate_monomials` for MWC identifiability
- `src/rate_eq_derivation.jl` — King-Altman/Cha rate equation derivation via `@generated` functions; parameters API; identifiability checks; kcat computation (`_is_ss_rate_constant`, `_kcat_components`, `_kcat_forward`) and `rescale_parameter_values`; OligomericEnzymeMechanism MWC rate equation assembly (`_build_oligomeric_rate_body`, `rate_equation_string`, `structural_identifiability_deficit`)
- `src/thermodynamic_constr_for_rate_eq_derivation.jl` — Haldane/Wegscheider thermodynamic constraints, dependent parameter elimination, `_build_rate_body` for `@generated rate_equation`
- `src/fitting.jl` — `FittingProblem`, `loss!`, `fit_rate_equation` using Optimization.jl
- `src/mechanism_enumeration.jl` — `SiteState`/`EnzymeFormSpec`/`MechanismSpec`/`MechanismIterator` types, `enumerate_enzyme_forms` (all valid enzyme forms from reaction), `enumerate_mechanisms` (catalytic cycle → regulator partitioning → dead-end → RE/SS with G≤7 cap → oligomeric expansion), `compile_mechanism` (MechanismSpec → EnzymeMechanism or OligomericEnzymeMechanism)

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
- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl` defines shared reactions for mechanism enumeration tests — must be included before those tests
- Mechanism enumeration tests use data-driven `EnumerationTestSpec` approach — verification helpers (`_compute_expected_dead_end_count`, `_compute_expected_n_total`) use only public struct fields, no `EnzymeRates._*` calls
- `enumerate_mechanisms` pipeline: Catalytic → regulator partitioning (2^n_reg) → dead-end → RE/SS with G≤7 cap → oligomeric expansion (when `catalytic_n > 0`)
- `MechanismTestSpec` has optional `analytical_kcat_fn` field for per-mechanism kcat formula validation
- kcat/rescaling tests (scale invariance, rate proportionality, V≈1, custom target) run for ALL mechanism specs in the main `run_all_tests` loop — not in a separate file
