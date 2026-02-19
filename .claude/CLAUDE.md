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
- **Stuck REPL**: If the REPL is stuck on a long-running computation, kill it with `pkill -f mcp_julia_server` via Bash, wait a few seconds, then call `exec_julia` again — the server auto-restarts

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

## Lessons Learned

### Rate Equation Derivation
- When all RE forms are in one group (G=1), the SS isomerization step flux IS the overall rate (no sign correction needed)
- `_compute_alpha` BFS handles forward/reverse RE traversal with K parameters
- `is_k_parameter` must match both `k1f`/`k1r` patterns AND `K1`/`K2` patterns (but not `Keq`)
- K convention: binding RE steps use Kd (dissociation), non-binding (product release) use Ka (forward eq const)
- K→1/K inversion must be applied consistently in: rate expr, dep_exprs substitution, constraint strings, and test helpers
- When `poly_str` displays inverted params with multiple denominators, need parentheses: `A * B / (K1 * K2)`

### Constraint Handling
- Constraint substitution operates at raw polynomial level (before K→1/K inversion)
- When constraining K params, both must be same type (both binding or both non-binding) to avoid inversion mismatch
- In DSL constraint parsing, use `Ref` for mutable coeff and read back with `coeff_ref[]`
- HW constraint matrix must merge constrained columns into replacement columns before Gaussian elimination

### Mechanism Enumeration
- Activator/inhibitor exclusivity: a regulator is EITHER activator OR inhibitor, never both. `_enumerate_dead_end_configs` excludes reg positions occupied in any topo form
- Essential activator: only entry binding edge (bare enzyme → shadow), not all base→shadow binding edges. Gives n_cat+1 topo forms (E, EA, EAS, EAP), not 2×n_cat
- Dead-end formula `(2^r_inh)^n_topo` is only valid when `max_forms` is unconstrained
- Ping-pong (CX/NX) product release must be TWO elementary steps: isomerization (product forms on enzyme, substrate → residual) then release (product leaves). A single combined step (e.g., `E_A → E_X + P`) violates detailed balance. The intermediate form (e.g., `E_X_P`) must be in the cycle.
- `_is_valid_isomerization` uses atom balance only (`_core_atoms` equality) for ping-pong (residual) case — the strict all-subs/all-prods occupancy pattern only applies to standard isomerizations

### Type System and Compatibility
- Default `eq_steps` = all false preserves backward compatibility with 2-arg constructor
- ParamConstraints default `()` preserves backward compat with 3-arg constructor
- `_binding_K_symbols` identifies binding steps: metabolite on LHS, enzyme-only on RHS
- JET requires `::SubString` type assertions on regex captures to avoid Union{Nothing,SubString} errors

## MCP REPL + Revise Pitfall

- After `git stash`/`git stash pop` (or any bulk file operation that cycles through intermediate states), **always restart the MCP REPL** (`pkill -f mcp_julia_server`) rather than relying on `Revise.revise()`. Revise can silently retain stale compiled methods, producing wrong numerical results that look plausible. Always confirm with `Pkg.test()` in a fresh process before trusting REPL results after stash operations.
- **Never use `exit()` in the REPL** — it kills the MCP connection permanently for the session. Use `pkill -f mcp_julia_server` via Bash instead (auto-restarts on next `exec_julia` call).
- **Use Bash with timeout for potentially slow computations** (e.g., enumeration with regulators). The REPL has no easy Ctrl+C — if a call blocks, the only recovery is `pkill -f mcp_julia_server`.

## Testing

- Tests include Aqua (quality) and JET (static analysis)
- Don't leave profiling deps (SnoopCompile) in Project.toml — Aqua stale deps check will fail
- `test/mechanism_definitions_for_test_enzyme_derivation.jl` defines shared mechanisms used by multiple test files — must be included before those tests
- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl` defines shared reactions for mechanism enumeration tests — must be included before those tests
- Mechanism enumeration tests use data-driven `EnumerationTestSpec` approach via `enumerate_mechanism_stages` — verification helpers (`_compute_expected_dead_end_count`, `_compute_expected_n_total`) use only public struct fields, no `EnzymeRates._*` calls
