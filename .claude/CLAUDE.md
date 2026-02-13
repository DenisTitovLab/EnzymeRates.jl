# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run full test suite
julia --project -e 'using Pkg; Pkg.test()'

# Run a single test file (e.g., test_fitting.jl)
julia --project -e 'using EnzymeRates, Test, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_fitting.jl")'

# Run tests interactively (useful for debugging)
julia --project
# then: using EnzymeRates, Test; include("test/runtests.jl")
```

## Key Architecture Decisions

- `EnzymeMechanism{Species,Reactions,EquilibriumSteps,ParamConstraints}` is a singleton type encoding mechanism info in type parameters
- `EnzymeReaction{S,P,R}` similarly encodes reactions in types
- Each unique mechanism = unique type → affects compilation time
- `@generated` functions used for compile-time computation (metabolites, graph, stoich_matrix, rate_equation)

## Source Layout

- `src/types.jl` — `EnzymeReaction`, `EnzymeMechanism` structs, accessors, `RateEquationMode` hierarchy
- `src/dsl.jl` — `@enzyme_reaction` and `@enzyme_mechanism` macros
- `src/sym_poly_for_rate_eq_derivation.jl` — Symbolic polynomial algebra (`Poly` type) used by rate equation derivation
- `src/rate_eq_derivation.jl` — King-Altman/Cha rate equation derivation via `@generated` functions; parameters API (`parameters`, `fitted_params`, `param_count_estimate`)
- `src/rate_eq_rewriting.jl` — Haldane/Wegscheider thermodynamic constraints, dependent parameter elimination, `_build_rate_body` for `@generated rate_equation`
- `src/fitting.jl` — `FittingProblem`, `loss!`, `fit_rate_equation` using Optimization.jl
- `src/mechanism_enumeration.jl` — `enumerate_enzyme_forms` (all valid enzyme forms), `enumerate_mechanisms` (reaction graph → cycle enumeration → topology combination → dead-end lattice → RE/SS × equivalent step constraints)
- `src/mechanism_selection.jl` — `select_mechanism` (beam search + leave-one-figure-out CV + 1-SE rule)

## Performance Pattern: Function Barriers

- `enumerate_mechanisms` uses a function barrier (`_enumerate_mechanisms_impl`) with `@nospecialize` on reaction-specific args
- This prevents recompilation of the entire enumeration logic for each new reaction type
- Raw mechanism data (tuples) is collected and deduplicated BEFORE creating EnzymeMechanism types
- Result: ~5x faster first call, ~18x faster different-reaction calls

## Testing

- Tests include Aqua (quality) and JET (static analysis)
- Don't leave profiling deps (SnoopCompile) in Project.toml — Aqua stale deps check will fail
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` defines shared mechanisms used by multiple test files — must be included before those tests
