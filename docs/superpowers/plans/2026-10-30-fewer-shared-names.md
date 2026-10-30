# Fewer Shared Names Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cut the non-exported functions crossing each module boundary to the counts in the spec,
by storing facts on the mechanism, moving shared graph analyses to Mechanisms, splitting
`parameters`, and renaming what the name review approved.

**Architecture:** Same branch and PR (`deep-modules-simplification`, version 0.10.0). First the
verification tools are rebuilt (they were lost with `/tmp`), then the stored facts land with an
allocation pin, then callers switch and old names go, then the moves, the API changes and the
renames, then the docs. Every commit keeps laws, parameter names, child sets and printed text
identical (census), except where a task says it changes an API.

**Tech Stack:** Julia 1.12–1.13; Test, TestEnv, Aqua, JET; Documenter + DocumenterMermaid.

**Spec:** `docs/superpowers/specs/2026-10-30-fewer-shared-names-design.md` (goals 1–10, the
interface table, Performance). Read it before any task.

## Global Constraints

- **CLAUDE.md applies in full** (TDD; smallest change; pristine output; `_testhelper_` prefix;
  ABOUTME headers; 92-character lines; evergreen comments; no coverage loss). The spec is
  Denis's permission for the changes it names.
- **Version stays `0.10.0`.** No change to `rate_equation`, `compile_mechanism` or the singleton
  types (`EnzymeMechanism{Sig}`, `AllostericEnzymeMechanism{Sig}`).
- **One Julia process at a time** (7.7 GB, no swap); `--heap-size-hint=2500M`. Scratch files go
  in `~/claude-scratchpad/fsn/` (persistent), never in `/tmp` and never in the repo.
- **Focused test run:** `julia --project --heap-size-hint=2500M -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random, REPL; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<Module>/<file>.jl")'`.
  **Full suite** once before each task's commit:
  `julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'` (pristine).
- **Census gate** (every commit that touches `src/`): the Task 0 census of the tree compares
  equal to the baseline built at `a6eb1d1b` — zero differences in printed text, fitted names
  and rate-law text.
- **Benchmark gate:** the Task 0 benchmarks, paired ABBA against the task's parent commit;
  more than 10% slower on any of them, reproduced in a second ABBA set, means stop and report
  BLOCKED. Julia 1.13.
- **Module map and interface tables:** a task that changes which names cross modules regenerates
  the map (`julia test/update_module_map.jl`) and edits the per-module interface tables in
  `docs/src/developer.md` in the same commit; the boundary test checks both.
- **Commit trailer:** `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` then
  `Claude-Session: https://claude.ai/code/session_01BFcXywjZosSMs2X8r6rZip`. Never push.
- **Names fixed by the spec and Denis:** `rapid_equilibrium_segments`, `has_turnover_cycle`,
  `has_unbalanced_cycle` (in Mechanisms); `RateConstants` (not exported);
  `dependent_parameter_values` (exported); `check_valid_rate_equation(Bool, m)`;
  `fitted_parameter_count`; `EnzymeForm`; `Metabolite`; `starting_mechanisms`;
  `rate_equation_term_count`; `test_rate_equation` keeps its name and is not exported.
- **Stored fields** (this plan's names, used by every task):
  - `Mechanism`: `modifiers::Vector{Symbol}`, `groups::Vector{Vector{Int}}`,
    `step_kinds::Vector{Symbol}`, `toward_products::BitVector`,
    `step_sides::Vector{Tuple{String, String}}`;
  - `AllostericMechanism`: `active_steps::BitVector`, `inactive_steps::BitVector` (today's
    `conformation_steps(am, :active/:inactive)`), `excluded_active::BitVector`,
    `excluded_inactive::BitVector` (per step: does its group's tag exclude that conformation);
  - a plain `Mechanism` answers the allosteric questions through `catalytic` and its callers'
    existing branches; it gets no allosteric fields.
  - Derived fields are left out of `==` and `hash`.

## Review Focus

1. **A mechanism built from the singleton (the `Sig` decode) carries the same stored facts as
   one built by the constructor** — otherwise a compiled mechanism's law reads empty groups
   (Task 1, "stored facts survive the singleton round trip").
2. **Two mechanisms equal except for a derived field cannot exist, and `unique!` still dedups
   enumerated children** (Task 1, "derived fields stay out of == and hash": equal canonical
   fields ⇒ equal and hash-equal).
3. **A user who called `parameters(m, values)` gets a clear error pointing at
   `dependent_parameter_values`** rather than a `MethodError` (Task 4).
4. **`check_valid_rate_equation(Bool, m)` never throws, for an invalid mechanism too; the
   one-argument form throws today's message** (Task 5).
5. **The rebuilt census really detects a difference** — a deliberately changed law in a scratch
   copy makes `compare.jl` fail (Task 0, negative control).

---

### Task 0: Rebuild the census and benchmark tools

**Files (outside the repo):** create `~/claude-scratchpad/fsn/census/{env/,census.jl,compare.jl}`,
`~/claude-scratchpad/fsn/bench/{bench_enum.jl,bench_guard.jl,bench_pass1.jl,pair.sh}`,
`~/claude-scratchpad/fsn/gates.md`.

**Interfaces:**
- Produces: `census.jl <repo> <outdir>` writes one TSV per reaction; `compare.jl <a> <b>` prints
  `RESULT: PASS` or the differences; each bench script prints `TOTAL <seconds>` (median of 3 runs
  after a warm-up); `pair.sh <bench> <parent-sha>` runs new, parent, parent, new (parent from a
  `git worktree` in `~/claude-scratchpad/fsn/wt`, removed afterwards).
- Scripts reach the package through one `const ER = EnzymeRates` and small helper functions, so a
  rename is re-pointed in one place (list the internals each script uses in its header).

- [ ] **Step 1: The census.** Six reactions: uni-uni, ordered bi-bi (A, B → P, Q), ping-pong
  (with a shared-site pair), ter-bi, a reaction with an allosteric activator and a competitive
  inhibitor, and one with an allosteric regulator on a second site. For each: every mechanism of
  `init_mechanisms`, then `expand_mechanisms` on a fixed sample of parents (seeded RNG, at most
  40 parents per level, two levels). Per mechanism record `sprint(show, m)`, its fitted parameter
  names and `rate_equation_string` (with coefficients), or the guard's rejection message.
- [ ] **Step 2: Negative control.** Build the baseline at `a6eb1d1b` twice → `compare` PASS; then
  hand-edit one law in a copy of the output → `compare` reports exactly that difference.
- [ ] **Step 3: Benchmarks.** `bench_enum` (init + expansion over the census sample),
  `bench_guard` (`check_valid_rate_equation` over every census mechanism, rebuilt with
  `DSL.parse_mechanism`), `bench_pass1` (identify's per-mechanism record over 2,000 census
  mechanisms). Record baselines at `a6eb1d1b` in `gates.md`.
- [ ] **Step 4: No commit** (nothing in the repo changes). Report paths and baseline numbers.

### Task 1: Store the mechanism's facts

**Files:** `src/Mechanisms/types.jl` (structs, constructors, `_unchecked_mechanism`, the `Sig`
decode, `==`/`hash`), `src/Mechanisms/accessors.jl`, move the conformation code
(`excludes`, `exclusion_components`, `conformation_steps`) from `src/Derivation/conformations.jl`
into Mechanisms as private helpers the constructor calls; `test/Mechanisms/*`.

**Interfaces:** Produces the stored fields in Global Constraints. In this task the existing
functions (`modifiers`, `groups`, `step_kind`, `runs_toward_products`, `step_sides`,
`Derivation.conformation_steps`, `Derivation.excludes`) keep working by reading the fields, so no
caller changes yet.

- [ ] **Step 1: Pin today's allocations.** At BASE, measure bytes allocated by one
  `Mechanism(substrates, products, nodes, steps)` (and the allosteric constructor) for three
  fixtures: Michaelis–Menten, the largest bi-bi `init_mechanisms` seed, and one allosteric
  mechanism from `MECHANISM_TEST_SPECS`. Write testset "constructor allocation budget" in
  `test/Mechanisms/test_reaction_graph.jl`: each fixture's `@allocated` (measured inside a
  `_testhelper_` function after a warm-up call) is `≤ ceil(Int, 1.3 * <BASE bytes>)`, the BASE
  bytes written as constants with a comment naming the commit. It passes at BASE.
- [ ] **Step 2: Failing tests:** "stored facts equal the old functions" — over every
  `MECHANISM_TEST_SPECS` mechanism and the bi-bi seeds: `m.modifiers == modifiers(m)`,
  `m.groups == groups(m)`, `m.step_kinds == [step_kind(m, k) for k in eachindex(m.steps)]`,
  likewise `toward_products`, `step_sides`, and for allosteric ones the four masks against
  `Derivation.conformation_steps`/`Derivation.excludes`; "stored facts survive the singleton round
  trip" (`Mechanism(compile_mechanism(m))` has equal fields, all of them); "derived fields stay out
  of == and hash". Run: FAIL (no such fields).
- [ ] **Step 3: Implement**; the old functions return the stored fields.
- [ ] **Step 4:** focused Mechanisms and Derivation tests, full suite, census PASS, allocation
  budget green, benchmarks (all three) within the gate. **Commit:** "Mechanisms store their
  derived facts".

### Task 2: Callers read the stored facts; delete the accessor functions

**Files:** every `src/` caller of `modifiers`, `groups`, `step_kind`, `runs_toward_products`,
`step_sides`, `conformation_steps`, `excludes` (Derivation, Enumeration, DSL); the Mechanisms and
Derivation `API` lists; tests; the module map and interface tables.

**Interfaces:** Deleted: `Mechanisms.modifiers`, `groups`, `step_kind`, `runs_toward_products`,
`step_sides` (the functions); `Derivation.conformation_steps`, `Derivation.excludes`. If
Enumeration's `_state_conformations` (a bare tag, no mechanism) cannot be served by a stored field,
keep one function for the tag rule in Mechanisms' API (name it `excludes`) and say why in the
report. The DSL stops calling `catalytic` (it reads the field it needs).

- [ ] **Step 1:** Remove the names from the `API` lists; run the boundary test: FAIL at every
  remaining caller. Convert callers and tests (tests of the deleted functions become assertions on
  the fields; no assertion about mechanism content is dropped).
- [ ] **Step 2:** regenerate the module map; update the interface tables.
- [ ] **Step 3:** full suite, census, benchmarks. **Commit:** "Read stored mechanism facts with
  dots".

### Task 3: Shared graph analyses move to Mechanisms

**Files:** `src/Derivation/segments.jl` → `src/Mechanisms/segments.jl` (the segment walk,
`re_segments` and `segment_graph` merged, `has_turnover_cycle`, `has_unbalanced_cycle`); callers in
Derivation and Enumeration; `test/Derivation/test_segments.jl` → `test/Mechanisms/`; API lists;
map and tables.

**Interfaces:** Produces `Mechanisms.rapid_equilibrium_segments(m, on, re = <m's own flags>)` →
one record with today's `re_segments` fields (`segment_of`, `segments`, `offsets`, `bottom`) plus
`segment_graph`'s `steps` (the kept steady-state steps between segments);
`Mechanisms.has_turnover_cycle(m, on)`; `Mechanisms.has_unbalanced_cycle(from, to, weight, on,
parent, offset)`. Deleted from Derivation's API: `re_segments`, `segment_graph`,
`has_turnover_cycle`, `has_unbalanced_cycle`.

- [ ] **Step 1:** Move the tests first, retargeted (`Mechanisms.rapid_equilibrium_segments` for both
  former callers' cases, asserting today's values): FAIL. Then move the code verbatim apart from the
  merge; convert callers.
- [ ] **Step 2:** map, tables, full suite, census, benchmarks. **Commit:** "Rapid-equilibrium
  segments and turnover cycles live in Mechanisms".

### Task 4: `parameters` returns names only

**Files:** `src/Derivation/rate_eq_derivation.jl`, `src/Derivation/Derivation.jl`,
`src/Fitting/fitting.jl`, `src/EnzymeRates.jl` (export), `docs/src/` pages that call
`parameters(m, values)` (the textbooks page) or name `scale_rate_constants`; tests.

**Interfaces:**
- `struct RateConstants end` in Derivation (in its `API`, not exported):
  `parameters(m, RateConstants())` → the steady-state rate constant names of both conformations
  (today's `_ss_rate_constant_names`, as a tuple in a fixed order), on concrete and compiled
  mechanisms.
- `dependent_parameter_values(m, values::NamedTuple) -> NamedTuple`: only the dependent constants,
  in name order, with today's checks; error text starts `dependent_parameter_values:`. Exported.
- `parameters(m, values::NamedTuple)` is removed; calling it errors with a message naming
  `dependent_parameter_values` (a method that throws that message is acceptable).
- Deleted: `scale_rate_constants`. Fitting scales the names `RateConstants()` returns.

- [ ] **Step 1: Failing tests:** "RateConstants lists what Fitting scales" (for every spec
  mechanism, equal as a set to the names today's `scale_rate_constants` changes); port the
  `parameters(m, values)` testsets to `merge(values, dependent_parameter_values(m, values))` with
  identical expected values; "parameters(m, values) points to dependent_parameter_values"
  (`@test_throws` with the message). FAIL.
- [ ] **Step 2:** implement; Fitting's default and reference scaling keep their tests green.
- [ ] **Step 3:** docs, map, tables, full suite, census, bench_pass1. **Commit:** "parameters
  returns names; dependent_parameter_values; RateConstants" (body states the API change).

### Task 5: Function renames

**Files:** every `src/`, `test/`, `docs/src/` use of the renamed functions; API lists; map;
tables; CLAUDE.md lines that name them.

**Interfaces:** `check_valid_rate_equation(m)` throws (today's message) and returns `true`;
`check_valid_rate_equation(Bool, m)::Bool` never throws (replaces `; throw = false`);
`parameter_count` → `fitted_parameter_count`; `base_mechanisms` → `starting_mechanisms`;
`rate_equation_complexity` → `rate_equation_term_count`.

- [ ] **Step 1: Failing test:** "check_valid_rate_equation(Bool, m) never throws" — a valid and an
  invalid fixture give `true`/`false`, `@inferred`; the one-argument form throws today's message.
  Rename the other three in the tests first. FAIL.
- [ ] **Step 2:** rename in `src/` and docs; map; tables; full suite; census. **Commit:**
  "Rename for clarity: check_valid_rate_equation(Bool, m), fitted_parameter_count,
  starting_mechanisms, rate_equation_term_count".

### Task 6: Type renames

**Files:** every use of `Node` (the mechanism type) and `Met`; API lists; docs; CLAUDE.md.

**Interfaces:** `Node` → `EnzymeForm`, `Met` → `Metabolite` (same fields). `show` output and the
census unchanged (the printed DSL never names these types).

- [ ] **Step 1:** rename in the tests first (FAIL), then in `src/`. Check that no local variable or
  unrelated identifier containing `Met`/`Node` was changed (review the diff by eye; `Metabolites`
  of reactions and `metabolites()` stay).
- [ ] **Step 2:** map; tables; full suite; census. **Commit:** "Rename Node to EnzymeForm and Met to
  Metabolite".

### Task 7: Exports, CLAUDE.md and docs

**Files:** `src/EnzymeRates.jl` (drop `test_rate_equation` from the exports; it stays in Fitting's
`API`), `docs/src/api.md`, `docs/src/developer.md`, the fitting and identify docs pages that call
`test_rate_equation` bare, `.claude/CLAUDE.md`, `ROADMAP.md`.

- [ ] **Step 1:** Update the export test (if one pins the export list) first; then the export;
  docs call `EnzymeRates.test_rate_equation`.
- [ ] **Step 2: CLAUDE.md.** Add the principle in one or two sentences to "Module boundaries":
  keep names shared across modules to a minimum, so each module is deep and can be worked on
  knowing only the other modules' interfaces. Update "API Design" (the type names, the
  construction surface) and any line naming a renamed or deleted function.
- [ ] **Step 3:** docs build with doctests (`julia --project=<docsenv> docs/make.jl`; rebuild a
  docs env in `~/claude-scratchpad/fsn/docsenv` if needed), full suite on Julia 1.13 and 1.12,
  final census, final benchmarks paired against `a6eb1d1b`. **Commit:** "test_rate_equation is
  internal; CLAUDE.md states the shared-names principle".
