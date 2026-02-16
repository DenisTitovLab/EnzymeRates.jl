# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run full test suite (cold — pays precompilation + JIT cost every time)
julia --project -e 'using Pkg; Pkg.test()'
```

## MCP Julia REPL Server

A persistent Julia session is available via MCP (`.mcp.json`). Claude Code auto-starts it on session launch — no manual setup needed.

- **Tool**: `exec_julia` — execute Julia code in a persistent session (the only tool; no static resources)
- **Pre-loaded**: Revise, EnzymeRates, Test, Random, LinearAlgebra
- **Revise**: Source edits in `src/` are picked up automatically
- **Introspection**: Use `exec_julia` to query exported names (`names(EnzymeRates)`), method tables (`methods(f)`), docstrings (`@doc f`), type hierarchies (`subtypes`/`supertypes`), and source locations
- **include()**: Relative paths resolve from project root (e.g., `include("test/test_fitting.jl")`)
- **First start**: ~30-60s (package loading + JIT). Subsequent calls are fast.
- **Startup timeout**: If the server fails to connect, start Claude Code with `MCP_TIMEOUT=120000 claude` to allow enough time for Julia startup.
- **Server script**: `.claude/mcp_julia_server.jl`

## Key Architecture Decisions

- `EnzymeMechanism{Species,Reactions,EquilibriumSteps,ParamConstraints}` is a singleton type encoding mechanism info in type parameters
- `EnzymeReaction{S,P,R}` similarly encodes reactions in types
- Each unique mechanism = unique type → affects compilation time
- `@generated` functions used for compile-time computation (metabolites, graph, stoich_matrix, rate_equation)

## Source Layout

- `src/types.jl` — `EnzymeReaction`, `EnzymeMechanism` structs, accessors, `RateEquationMode` hierarchy
- `src/dsl.jl` — `@enzyme_reaction` and `@enzyme_mechanism` macros
- `src/sym_poly_for_rate_eq_derivation.jl` — Symbolic polynomial algebra (`Poly` type) used by rate equation derivation
- `src/rate_eq_derivation.jl` — King-Altman/Cha rate equation derivation via `@generated` functions; parameters API (`parameters`, `fitted_params`, `param_count_estimate`); identifiability checks
- `src/rate_eq_rewriting.jl` — Haldane/Wegscheider thermodynamic constraints, dependent parameter elimination, `_build_rate_body` for `@generated rate_equation`
- `src/fitting.jl` — `FittingProblem`, `loss!`, `fit_rate_equation` using Optimization.jl
- `src/mechanism_enumeration.jl` — `SiteState`/`EnzymeFormSpec`/`MechanismSpec`/`PreRessEntry`/`MechanismIterator` types, `enumerate_enzyme_forms` (all valid enzyme forms from reaction), `enumerate_mechanisms` (catalytic cycle enumeration → dead-end lattice → lazy RE/SS × equivalent step constraints), `enumerate_mechanism_stages` (returns intermediate results at each pipeline stage)

## Performance Pattern: @nospecialize

- Enumeration functions use `@nospecialize` on `EnzymeReaction` args to prevent recompilation for each reaction type
- Raw mechanism data (tuples) is collected and deduplicated BEFORE creating EnzymeMechanism types

## Mechanism Enumeration Architecture

- `MechanismIterator` is lazy — stages 1-3 (catalytic, activator, dead-end) are eager; stage 4 (RE/SS + constraints) generates `MechanismSpec` on demand via `PreRessEntry` state machine
- `length(iter)` is O(1) — precomputed via `_count_ress_variants`
- `MechanismSpec` instances from the same iterator share `form_names`/`form_atoms` vectors (no redundant copies)
- `enumerate_mechanism_stages` exposes all intermediate pipeline results for testing and inspection
- Stages 1-3 with generous `max_forms` can still be slow for multi-regulator reactions (dead-end enumeration is combinatorial)

## Dead-End Mechanism Combinatorics

For reactions with r regulators, each regulator is either an activator (part of the catalytic topology) or an inhibitor (creates dead-ends). No mixed roles: an activator does not participate in dead-end formation.

- Dead-end configs per activator config = `(2^r_inh)^n_topo`
  - `r_inh`: number of inhibitor regulators (binding site never occupied in any topology form)
  - `n_topo`: number of forms in the catalytic topology
- For Uni-Uni (n_cat=3 base catalytic forms):
  - No activator: n_topo = 3
  - Essential activator: n_topo = n_cat + 1 = 4 (E, EA, EAS, EAP)
  - Non-essential activator: n_topo = 2 × n_cat = 6 (bare + activated cycles)
- Each activator can be essential or non-essential → 2 configs per activator choice
- Total = Σ over activator configs of `(2^r_inh)^n_topo`
- Test helper `_compute_expected_dead_end_count` verifies this formula

## Testing

- Tests include Aqua (quality) and JET (static analysis)
- Don't leave profiling deps (SnoopCompile) in Project.toml — Aqua stale deps check will fail
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` defines shared mechanisms used by multiple test files — must be included before those tests
- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl` defines shared reactions for mechanism enumeration tests — must be included before those tests
- Mechanism enumeration tests use data-driven `EnumerationTestSpec` approach via `enumerate_mechanism_stages` — verification helpers (`_compute_expected_dead_end_count`, `_compute_independent_ress_count`) use only public struct fields, no `EnzymeRates._*` calls
