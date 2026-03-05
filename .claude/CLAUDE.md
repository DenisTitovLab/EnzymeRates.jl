# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

When saving new lessons learned or project knowledge, write them to this file (`.claude/CLAUDE.md`) so they are shared across all machines via git. Do not use the local auto memory for project-specific knowledge.

## Package Goal

EnzymeRates.jl identifies the best enzyme rate equation from kinetic data. Given a reaction definition and experimental rate measurements at varying concentrations, the package enumerates all biochemically valid mechanisms, fits each to data, and selects the one with fewest parameters that adequately describes the data (cross-validation). See `SPEC.md` for the full API specification.

**Primary use case**: `EnzymeReaction` + data → `IdentifyRateEquationProblem` → `identify_rate_equation()` → `IdentifyRateEquationResults`

**Secondary use cases**: manually define mechanisms via `@enzyme_mechanism` and derive/fit rate equations.

## API Design (see SPEC.md)

- **16 exported symbols** (planned): 5 types, 2 macros, 2 constants (`Full`, `Reduced`), 7 functions. Currently 13 — `IdentifyRateEquationProblem`, `IdentifyRateEquationResults`, `identify_rate_equation` are pending implementation.
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
- `@generated` functions used for compile-time computation (metabolites, graph, stoich_matrix, rate_equation)

### Canonical Step Form
- The `EnzymeMechanism` constructor normalizes RE steps so metabolite is always on LHS (binding direction): `[E, S] ⇌ [ES]`, never `[ES] ⇌ [E, S]`
- SS steps are NOT canonicalized (swapping kf↔kr would break analytical test formulas)
- After canonicalization, all RE metabolite K params are binding Kd (displayed as `1/K`). Non-binding RE steps (pure isomerization) retain Ka convention.
- `_binding_K_symbols` relies on this invariant: checks only for metabolite on LHS, no RHS check needed

## Source Layout

- `src/types.jl` — `EnzymeReaction`, `EnzymeMechanism`, `OligomericEnzymeMechanism` structs; `EnzymeMechanism` and `OligomericEnzymeMechanism` accessors; `RateEquationMode` hierarchy
- `src/dsl.jl` — `@enzyme_reaction` and `@enzyme_mechanism` macros (handles both `EnzymeMechanism` and `OligomericEnzymeMechanism` DSL)
- `src/sym_poly_for_rate_eq_derivation.jl` — Symbolic polynomial algebra (`Poly` type); `_rename_poly_T`, `_count_oligomeric_rate_monomials` for MWC identifiability
- `src/rate_eq_derivation.jl` — King-Altman/Cha rate equation derivation via `@generated` functions; parameters API; identifiability checks; OligomericEnzymeMechanism MWC rate equation assembly (`_build_oligomeric_rate_body`, `rate_equation_string`, `structural_identifiability_deficit`)
- `src/thermodynamic_constr_for_rate_eq_derivation.jl` — Haldane/Wegscheider thermodynamic constraints, dependent parameter elimination, `_build_rate_body` for `@generated rate_equation`
- `src/fitting.jl` — `FittingProblem`, `loss!`, `fit_rate_equation` using Optimization.jl
- `src/mechanism_enumeration.jl` — `SiteState`/`EnzymeFormSpec`/`MechanismSpec`/`PreRessEntry`/`MechanismIterator` types, `enumerate_enzyme_forms` (all valid enzyme forms from reaction), `enumerate_mechanisms` (catalytic cycle enumeration → dead-end lattice → lazy RE/SS × equivalent step constraints), `enumerate_mechanism_stages` (returns intermediate results at each pipeline stage)

## Vmax Normalization (kcat factoring)

### Problem
Different datasets often have unknown Vmax (= kcat * E_total, i.e., specific
activity). The fitting must identify binding constants and rate constant ratios
(shape parameters) independently of the unknown per-dataset scale.

### Mathematical basis
King-Altman rate equations have the structure:
```
v = E_total * N(k's, K's, concs) / D(K's, k's, concs)
```
Two key structural properties make Vmax factoring always possible:
1. **Numerator factoring**: After Haldane substitution, the numerator always
   factors as `f_fwd(k's) * (substrate_product - product_product / Keq)`.
   This holds for ALL mechanism types (ordered, random, ping-pong) because
   the Haldane relation is a thermodynamic identity.
2. **Uniform k-degree in denominator**: All denominator terms have the same
   total k-degree (G-1 for G groups), because each King-Altman spanning tree
   has exactly (G-1) edges. This guarantees that v/(E_total * kcat) is
   scale-invariant (k-degree 0).

### Design: Option B (post-fit kcat normalization)
- Fit with per-figure centering in log-space (unchanged `loss!`)
- After fitting, compute `kcat_fitted` analytically from the polynomial
  structure (ratio of numerator/denominator leading k-coefficients at
  saturating substrates)
- Normalize all SS rate constants: `k_norm = k_fitted / kcat_fitted`
  - This gives `kcat(k_norm) = 1` (kcat is homogeneous degree-1 in SS k's)
  - K's (binding constants) are unchanged (they're ratios k_f/k_r)
  - All k-ratios are unchanged (uniform scaling)
- For prediction: user provides Vmax (= kcat * E_total) as a separate
  quantity; the normalized params + Vmax fully determine the rate
- To recover true physical k values: `k_true = kcat_measured * k_norm`
  (where kcat_measured is independently measured)

### kcat properties
- kcat = f_fwd / leading_denom_k_coeff (always well-defined, even for
  ping-pong where denominator has no constant term)
- kcat is homogeneous degree-1 in SS k's, independent of RE K's
- For ordered Bi-Bi: kcat = k3f*k4f / (k3f+k4f)
- For ping-pong Bi-Bi: kcat = k2f*k4f / (k2f+k4f)
- For simple uni-uni with 1 RE + 1 SS: kcat = k2f

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
- Mechanism enumeration tests use data-driven `EnumerationTestSpec` approach via `enumerate_mechanism_stages` — verification helpers (`_compute_expected_dead_end_count`, `_compute_expected_n_total`) use only public struct fields, no `EnzymeRates._*` calls
