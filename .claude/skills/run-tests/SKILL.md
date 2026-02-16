---
name: run-tests
description: Run Julia tests via MCP REPL (fast, avoids recompilation)
---

# Run Tests

Run EnzymeRates test files using the MCP `exec_julia` tool (persistent session with Revise).

## Arguments

- Test file name(s) (e.g., `test_fitting`, `test_dsl test_types`), or `all` for the full suite.
- Omitting arguments runs all tests.

## Procedure

### 1. Run tests via `exec_julia`

**Shared definitions must always be included first.** The two definition files are:
- `mechanism_definitions_for_test_enzyme_derivation.jl` — needed by: `test_enzyme_derivation`, `test_fitting`, `test_mechanism_enum_of_enz_reaction`
- `reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl` — needed by: `test_mechanism_enum_of_enz_reaction`

Test files that do NOT need shared definitions: `test_accessors`, `test_types`, `test_dsl`, `test_aqua_jet`.

**For `all` or no arguments**, run the full suite:

```julia
include("test/runtests.jl")
```

**For specific test files**, prepend required definitions in the same `exec_julia` call. Example for `test_fitting`:

```julia
include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/test_fitting.jl")
```

### 2. Report results

Report test results (pass/fail counts) to the user. If any tests fail, show the failure details.

If the MCP server is not available, fall back to cold run:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
