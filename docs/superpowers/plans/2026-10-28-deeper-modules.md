# Deeper Modules Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make each EnzymeRates module deeper: the mechanism becomes plain data read with dots,
each function moves to the module whose job it is, two small modules (`Data`, `DSL`) join, and
the developer docs explain the result with a diagram.

**Architecture:** Same branch and PR as the first module split (`deep-modules-simplification`).
Storage-independent interface work comes first (Tasks 1–7), then the storage port with its
performance checkpoint (Task 8), then the callers move to dot notation and the old surface is
deleted (Tasks 9–14), then the docs (Task 15). Every commit keeps the census, the golden files
and printed text byte-identical, except where a task says it changes behaviour.

**Tech Stack:** Julia 1.10–1.13; Test, TestEnv, Aqua, JET; Optimization.jl; Documenter.

**Spec:** `docs/superpowers/specs/2026-10-28-deeper-modules-design.md` (§1–§9). Read it before
any task. The first split's spec and plan (`…/2026-10-27-deep-modules-simplification*`) give
background only.

## Global Constraints

- **CLAUDE.md applies in full:** TDD; the smallest change; pristine test output;
  `_testhelper_` prefix on every test helper; `ABOUTME:` two-line header on every new file;
  92-character lines, 4-space indentation; evergreen comments (no history, no "new/old/moved");
  never remove a comment unless it is false. The spec is Denis's permission for the rewrites it
  names; any other rewrite needs his.
- **Version stays `0.10.0`.** No bump in this PR.
- **One Julia process at a time** (7.7 GB, no swap). Pass `--heap-size-hint=2500M` to every run.
- **Focused test run:**
  `julia --project --heap-size-hint=2500M -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<Module>/<file>.jl")'`
  (`test/<file>.jl` for files directly in `test/`).
- **Full suite** once before each task's final commit:
  `julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`.
  Pristine: no Fail, Error, Broken or warning beyond the two Pkg manifest warnings recorded in
  `$S/dm/baseline.md`.
- **`rate_equation`** stays allocation-free and under CLAUDE.md's bound
  (`_testhelper_test_rate_equation_performance`). A change that would touch it: stop, report.
- **Census gate** (every commit): build the tree's corpus and compare with `DM_base`:
  `julia --project=$S/dm/census/env --heap-size-hint=2500M $S/dm/census/corpus.jl /home/denis.linux/.julia/dev/EnzymeRates $S/dm/census/DM_T<N> $S/dm/census/rename_phase3.toml`
  then `compare.jl DM_base DM_T<N>` must print `RESULT: PASS` with zero differences. The census
  and benchmark scripts reach internals through `ER.` names; when a task renames or moves a
  name they call, the task re-points the script (it lives in `$S`, outside the repo) and says
  so in its report.
- **Golden files stay byte-identical:** `test/fixtures/canonical_text_golden.txt`,
  `test/fixtures/phase2_init_golden.txt`, `test/reference/allosteric_golden_reference.txt`.
  `show` output of every mechanism and reaction is unchanged (spec, Non-goals).
- **Benchmark gate:** `bench_enum.jl`, `bench_guard.jl`, `bench_pass1.jl` in `$S/dm/`, each
  paired ABBA against the task's parent commit in the same session (ruling R8). More than 5%
  slower means stop: write the numbers to `$S/dm/gates.md` and report BLOCKED. Run the
  benchmark for every task that touches `src/Mechanisms`, `src/Derivation` or
  `src/Enumeration`.
- **Boundary test** (`test/test_boundaries.jl`) stays green at every commit. A task that adds
  or removes an interface name edits the module's `const API` and the test's expectations in the
  same commit.
- **Attribution.** Every commit message ends with:
  ```
  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01BFcXywjZosSMs2X8r6rZip
  ```
- **Scratchpad.** `$S` = `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/276b07d3-bafa-4735-bad9-ab925057d392/scratchpad`.
  Logs go to `$S/dm2/`, never into the repo.
- **Names fixed by the spec and Denis:** modules `Data`, `DSL`; types `RateData`, `Node`, `Step`,
  `Reactant`, `Regulator`, `IndependentParameters`, `FittedParameters`, `AllParameters`; fields
  `mwc_tags`, `effect`; functions `parameters`, `check_valid_rate_equation(m; throw = true)::Bool`,
  `rate_equation_dedup_key`, `test_rate_equation`, `parse_mechanism`, `step_sides`, `groups`,
  `step_kind` with values `:isomerization`, `:binding`, `:fused_binding`, `:chemistry`.

## Review Focus

1. **A relative fit without `reference` whose predicted rate is zero or non-finite at every data
   row.** Expect `retcode = :ReferenceFailed`, `loss = Inf`, no throw, so identify drops it
   (Task 4, "default scale fails softly").
2. **Saved identify results written before this PR** must parse back to equal mechanisms.
   `show` stays byte-identical and `parse_mechanism(sprint(show, m)) == m` for every golden
   mechanism (Task 8 and Task 13, the golden-file round trip).
3. **A macro used from a module that imports only the macro**
   (`using EnzymeRates: @enzyme_mechanism`) still expands after the macros move to `DSL`
   (Task 6, "macros work when imported alone").
4. **`parameters` agrees across representations:** `parameters(m, kind)` on a `Mechanism` and on
   its `compile_mechanism(m)` return the same tuple for all three kinds (Task 1).
5. **A mechanism written in another order** (nodes permuted, duplicate nodes, steps reversed,
   group labels renumbered) builds a `Mechanism` that is `==` and `hash`-equal to the canonical
   one (Task 13, "construction is canonical").

---

## Stage 1: Interfaces that do not depend on storage

### Task 1: One `parameters` with kinds

**Files:**
- Modify: `src/Derivation/rate_eq_derivation.jl` (parameter API), `src/Derivation/Derivation.jl`
  (`API`), `src/EnzymeRates.jl` (exports, imports), `src/Fitting/fitting.jl`,
  `src/identify_rate_equation.jl`, `docs/src/deriving/textbooks.md`, `docs/src/api.md`.
- Test: `test/Mechanisms/test_accessors.jl` (the `all_parameters` testsets),
  `test/Derivation/test_rate_eq_derivation.jl`, every test calling `fitted_params` or
  `all_parameters` (`git grep -n`).

**Interfaces:**
- Produces, in `Derivation`: `struct IndependentParameters end`, `struct FittedParameters end`,
  `struct AllParameters end`;
  - `parameters(m)` ≡ `parameters(m, IndependentParameters())` → `(fitted..., :Keq, :E_total)`;
  - `parameters(m, FittedParameters())` → today's `fitted_params(m)`, `@generated` on
    `AbstractEnzymeMechanism` subtypes (so `loss!` stays allocation-free);
  - `parameters(m, AllParameters())` → `(parameters(m)..., dependents...)`, dependents in the
    order `all_parameters` adds them (by name);
  - `parameters(m, values::NamedTuple)` → today's `all_parameters(m, values)`, same errors with
    the prefix `parameters:` in place of `all_parameters:`.
  - All four on `AbstractEnzymeMechanism`, `Mechanism` and `AllostericMechanism`.
- Deleted: `fitted_params`, `all_parameters`. Exported: the three kind types; `all_parameters`
  leaves the exports. `parameter_count` stays (ruling: Enumeration counts before it checks
  validity, `moves.jl:395–400`, and the `Mechanism` method of `FittedParameters` checks first).

- [ ] **Step 1: Write the failing tests.** In `test/Derivation/test_rate_eq_derivation.jl`, a
  testset "parameters kinds agree across representations" over every `MECHANISM_TEST_SPECS`
  mechanism `em` and `m = EnzymeRates.Mechanism(em)` (and allosteric ones via
  `AllostericMechanism`):
  ```julia
  for kind in (IndependentParameters(), FittedParameters(), AllParameters())
      @test parameters(m, kind) == parameters(em, kind)
  end
  @test parameters(em) == parameters(em, IndependentParameters())
  @test parameters(em, IndependentParameters()) ==
        (parameters(em, FittedParameters())..., :Keq, :E_total)
  @test (@allocated parameters(em, FittedParameters())) == 0
  ```
  Port every `all_parameters` testset to `parameters(m, values)` with identical assertions,
  and add `keys(parameters(m, values)) ⊆ parameters(m, AllParameters())` for a full `values`.
- [ ] **Step 2: Run the two test files.** Expected: FAIL (`IndependentParameters` undefined).
- [ ] **Step 3: Implement** in `rate_eq_derivation.jl`; replace every `fitted_params` caller
  with `parameters(m, FittedParameters())` and every `all_parameters` caller with
  `parameters(m, values)`; update `Derivation.API`, exports and the top-level `using` line;
  update the docs pages that show `all_parameters` or `fitted_params`.
- [ ] **Step 4: Run the focused files, then the full suite.** Expected: PASS, pristine.
  Re-point `$S/dm/census/corpus.jl` (`ER.fitted_params` → `ER.parameters(em,
  ER.FittedParameters())`) and the bench scripts if they call either name.
- [ ] **Step 5: Census and benchmark gates; commit.** Message: "One parameters with kinds:
  IndependentParameters, FittedParameters, AllParameters", body naming the removed exports.

### Task 2: One validity check; the dedup key

**Files:**
- Modify: `src/Derivation/guard.jl:1–30`, `src/Derivation/key.jl:150–276`,
  `src/Derivation/Derivation.jl`, every caller of `has_valid_rate_equation`,
  `check_valid_rate_equation`, `rate_equation_key` and `Mechanisms.mechanism_id`
  (`src/Enumeration/*`, `src/identify_rate_equation.jl`, `src/macros.jl`,
  `src/Derivation/admission.jl`), `src/Mechanisms/dsl.jl:649–650`, `src/Mechanisms/Mechanisms.jl`.
- Test: `test/Derivation/test_key.jl`, `test/Derivation/test_admission.jl`, the tests of
  `has_valid_rate_equation` and `mechanism_id` (`git grep -n`).

**Interfaces:**
- Produces: `check_valid_rate_equation(m; throw::Bool = true)::Bool` — with `throw = true`
  errors with the guard's message on a fault, else returns `true`; with `throw = false` returns
  whether the verdict has no fault. The keyword shadows `Base.throw` inside the body; the body
  must not call `throw`.
- Produces: `rate_equation_dedup_key(m) -> (; key, parameters, fit_key, fit_parameters, capped)`
  (renamed `rate_equation_key`). Its capped branch builds the text id itself:
  `bytes2hex(sha256(sprint(show, m)))[1:16]`, so key strings do not change.
- `mechanism_id` leaves `Mechanisms`; identify keeps a private `_mechanism_id(text)` with the
  same formula.
- Deleted: `has_valid_rate_equation`, `rate_equation_key`, `Mechanisms.mechanism_id`.

- [ ] **Step 1: Failing tests.** In `test_admission.jl`, testset "check_valid_rate_equation
  returns a Bool": for one well-posed and one rejected mechanism (take them from an existing
  rejection testset), `check_valid_rate_equation(good) === true`,
  `check_valid_rate_equation(bad; throw = false) === false`,
  `check_valid_rate_equation(good; throw = false) === true`, and `@test_throws` with today's
  message for `check_valid_rate_equation(bad)`. Add `@inferred` on both calls. In
  `test_key.jl`, rename the calls and add "a capped key keeps its text id": for a capped
  fixture already in the file, `k.key` equals the string the test pinned before the rename.
- [ ] **Step 2: Run them.** Expected: FAIL.
- [ ] **Step 3: Implement**; replace callers; delete the three names; update API lists.
- [ ] **Step 4: Focused files, full suite.** PASS, pristine. Re-point census/bench scripts.
- [ ] **Step 5: Gates; commit.** "check_valid_rate_equation returns a Bool; rate_equation_dedup_key".

### Task 3: The `Data` module

**Files:**
- Create: `src/Data/Data.jl` (module `Data`, `const API = (:RateData,)`), `src/Data/rate_data.jl`.
- Modify: `src/EnzymeRates.jl` (include `Data` first), `src/Fitting/fitting.jl` (move
  `rate_data`, `_rate_table`, `_reference`, `_check_anchors` out), `src/Fitting/Fitting.jl`,
  `src/identify_rate_equation.jl`, `test/test_boundaries.jl` (`_testhelper_MODULES`,
  `_testhelper_LAYERS`: `:Data => Symbol[]`, `:Fitting => [:Mechanisms, :Derivation, :Data]`).
- Test: create `test/Data/test_rate_data.jl` (moving the validation testsets from
  `test/Fitting/test_fitting.jl` and `test_fitting_scale.jl`); add it to `test/runtests.jl`.

**Interfaces:**
- Produces: `RateData{D<:NamedTuple}` with fields `data::D` (the validated
  `Tables.columntable`), `Keq::Float64`, `absolute::Bool`,
  `reference::Vector{Pair{NamedTuple, Float64}}` (empty in absolute mode);
  constructor `RateData(table, metabolite_names; Keq, absolute = false, reference = nothing,
  owner = "mechanism")` with today's checks and messages (`owner` is `"mechanism"` or
  `"reaction"`, as in `rate_data`).
- `FittingProblem` holds `rate_data::RateData` in place of its `data`, `Keq`, `absolute`,
  `reference` fields; `IdentifyRateEquationProblem` likewise. Their public constructors keep
  today's signatures. `loss!` stays allocation-free.
- Deleted: `Fitting.rate_data`.

- [ ] **Step 1: Failing test.** `test/Data/test_rate_data.jl` holds the moved validation
  testsets, calling `EnzymeRates.Data.RateData(...)` in place of `FittingProblem` where they
  test table checks only. Plus "RateData holds the validated table":
  `rd.Keq === 2.0`, `rd.absolute === false`, `rd.data.Rate == table.Rate`.
- [ ] **Step 2: Run it.** Expected: FAIL (`Data` undefined).
- [ ] **Step 3: Implement**, moving code verbatim where possible.
- [ ] **Step 4: Focused Fitting, Data and identify tests; boundary test; full suite.** The
  fitting allocation test (`loss!` 0 allocations) must still pass.
- [ ] **Step 5: Gates (bench_pass1 only: no Mechanisms/Derivation/Enumeration change); commit.**
  "Data module: RateData holds the validated rate table".

### Task 4: Default reference scale and `test_rate_equation`

Behaviour change and API change.

**Files:**
- Modify: `src/Fitting/fitting.jl`, `src/Fitting/Fitting.jl`, `src/Data/rate_data.jl`
  (`reference` optional in relative mode), `src/identify_rate_equation.jl:1025–1060`
  (cross-validation), `src/EnzymeRates.jl` (export `test_rate_equation`; drop
  `rescale_parameter_values`), `docs/src/` pages that show `rescale_parameter_values` or say
  `reference` is required.
- Test: `test/Fitting/test_fitting_scale.jl`, `test/Fitting/test_fitting.jl`,
  `test/Data/test_rate_data.jl`, `test/test_identify_rate_equation.jl`.

**Interfaces:**
- Relative mode with `reference === nothing` is valid. `fit_rate_equation` then multiplies the
  steady-state rate constants (`Derivation.scale_rate_constants`) by `1 / max_i |v_i|`, where
  `v_i` is `rate_equation` at data row `i`'s concentrations with `E_total = 1` and the problem's
  `Keq`. If that maximum is zero or non-finite: `loss = Inf`, `retcode = :ReferenceFailed`,
  unscaled params, no throw. With `reference`, behaviour is today's.
- Produces: `test_rate_equation(m, params::NamedTuple, table; Keq, absolute = false)::Float64`,
  the loss (`loss!`) of `m` at `params` on `table`. `params` must hold every name of
  `parameters(m, FittedParameters())` (extra names ignored); a missing one errors naming it.
  Exported.
- Identify's cross-validation becomes a private `_cross_validation_loss` in
  `identify_rate_equation.jl`: fit on the other groups (`fit_rate_equation`), score the held-out
  group (`test_rate_equation` or the same `loss!` call), same non-finite error message.
- Deleted: `rescale_parameter_values`, `cross_validation_loss`. Fitting's API:
  `(:FittingProblem, :fit_rate_equation, :test_rate_equation)`.

- [ ] **Step 1: Failing tests.** In `test_fitting_scale.jl`:
  - "default scale sets the fastest data row to 1": fit a simple mechanism in relative mode
    without `reference`; with `p = merge(fit.params, (Keq = Keq, E_total = 1.0))`,
    `maximum(abs(rate_equation(m, row_concs, p)) for row in data) ≈ 1`.
  - "default scale fails softly": data whose rows give a zero predicted rate everywhere (all
    substrate and product concentrations 0) → `fit.retcode === :ReferenceFailed`,
    `fit.loss == Inf`.
  - "test_rate_equation is the fitting loss": for a fit `fit` on problem `fp`,
    `test_rate_equation(m, fit.params, table; Keq) ≈ fit.loss` (relative, no reference — the
    loss is scale-free) and a missing parameter errors naming it.
  Change "reference is required in relative mode" into "relative mode accepts no reference".
  Delete the `rescale_parameter_values` testsets; list their names in the commit message.
- [ ] **Step 2: Run them.** Expected: FAIL.
- [ ] **Step 3: Implement.** Document in the `FittingProblem` docstring and the fitting docs page
  that default-scaled rate constants are in units of the fastest measured condition, not kcats.
- [ ] **Step 4: Focused Fitting, Data, identify tests; full suite.**
- [ ] **Step 5: Gates (bench_pass1); census (identify is not in the census, so PASS expected);
  commit.** "Default reference scale; test_rate_equation; drop rescale_parameter_values",
  body stating the API change.

### Task 5: Reaction records

**Files:**
- Modify: `src/Mechanisms/types.jl:1–135, 320–452` (metabolite and regulator types,
  `EnzymeReaction`, its accessors and `show`), `src/Mechanisms/dsl.jl` (`reaction_expr`),
  `src/Mechanisms/Mechanisms.jl` (`API`), every caller of `reactants`, `regulators`,
  `atoms`, `metabolite`, `regulator`, `reg_type`, `name` on a reaction metabolite,
  `substrates(r)`/`products(r)` on a reaction, `allowed_catalytic_multiplicities`,
  `shared_catalytic_site`, `AllostericRegulator`, `CompetitiveInhibitor`
  (`src/Enumeration/*`, `src/identify_rate_equation.jl`, `src/macros.jl`).
- Test: `test/Mechanisms/test_types.jl`, `test/Mechanisms/test_dsl.jl`, enumeration and
  identify tests that build reactions by hand.

**Interfaces:**
- Produces: `struct Reactant; name::Symbol; atoms::<builder's choice: Dict{Symbol,Int} or
  Vector{Pair{Symbol,Int}}>; end`; `struct Regulator; name::Symbol; kind::Symbol; effect::Symbol;
  end` (`kind ∈ (:allosteric, :competitive)`, `effect ∈ (:activator, :inhibitor,
  :unspecified)`, checked by the constructor with today's message wording);
  `EnzymeReaction` fields `substrates::Vector{Reactant}`, `products::Vector{Reactant}`,
  `regulators::Vector{Regulator}`, `allowed_catalytic_multiplicities::Vector{Int}`,
  `shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}`, all public; today's constructor checks
  (atom balance, names, decodability) and sort orders kept.
- Deleted: `Metabolite`, `Reactant` (abstract), `Regulator` (abstract), `Substrate`,
  `Product`, `AllostericRegulator`, `CompetitiveInhibitor`, `ReactantAtoms`,
  `DeclaredRegulator`, and the functions `reactants`, `regulators`, `atoms`, `metabolite`,
  `regulator`, `reg_type`, `allowed_catalytic_multiplicities`, `shared_catalytic_site`,
  `substrates(::EnzymeReaction)`, `products(::EnzymeReaction)`, `name(::Metabolite)`.
  `Mechanisms.API` gains `:Reactant`, `:Regulator`.

- [ ] **Step 1: Failing test.** In `test_types.jl`, "reaction fields are its declaration":
  build `@enzyme_reaction` with two substrates, a product, an `::Activator` allosteric
  regulator and a competitive inhibitor; assert `rxn.substrates[1].name`, `.atoms`,
  `rxn.regulators[1] == Regulator(:A, :allosteric, :activator)` (adapt to the fixture),
  `rxn.allowed_catalytic_multiplicities`, and `sprint(show, rxn)` equals today's text
  (copy the expected string from a run on the parent commit).
- [ ] **Step 2: Run it.** FAIL.
- [ ] **Step 3: Implement**; migrate callers and tests to fields.
- [ ] **Step 4: Focused Mechanisms, Enumeration, identify tests; full suite.**
- [ ] **Step 5: Census, benchmark (bench_enum), commit.** "EnzymeReaction holds plain records".

### Task 6: The `DSL` module

**Files:**
- Create: `src/DSL/DSL.jl` (module `DSL`; `const API = (:parse_mechanism,
  Symbol("@enzyme_reaction"), Symbol("@enzyme_mechanism"), Symbol("@allosteric_mechanism"))`),
  `src/DSL/read.jl` (the reader half of `src/Mechanisms/dsl.jl`), `src/DSL/macros.jl`
  (today's `src/macros.jl`).
- Modify: `src/Mechanisms/dsl.jl` keeps only the printer (`show` of mechanisms, `form_text`,
  `step_text` and their helpers); rename it `src/Mechanisms/text.jl`. `src/EnzymeRates.jl`
  (include `DSL` after `Derivation`; `using .DSL: @enzyme_reaction, @enzyme_mechanism,
  @allosteric_mechanism, parse_mechanism`; keep the exports), `src/identify_rate_equation.jl`
  (`_parse_mechanism` → `DSL.parse_mechanism`), `test/test_boundaries.jl` (module `:DSL`,
  layer `[:Mechanisms, :Derivation]`; the rule-1 parser must accept macro names in an `API`
  list).
- Test: move the reader tests of `test/Mechanisms/test_dsl.jl` into `test/DSL/test_dsl.jl`;
  keep printer tests in `test/Mechanisms/`. Update `test/runtests.jl`.

**Interfaces:**
- Produces: `DSL.parse_mechanism(text::AbstractString) -> Mechanism or AllostericMechanism`
  (today's `_parse_mechanism`, same checks and messages).
- `Mechanisms.API` loses `raw_mechanism`, `reaction_expr`, `check_uniquely_decodable`
  (private to `DSL` now, underscore-free names may stay or become `_`-prefixed; builder's call).
- The decodability rules for a reaction stay where the `EnzymeReaction` constructor can run
  them; only the mechanism-text check moves.

- [ ] **Step 1: Failing test.** In `test/DSL/test_dsl.jl`, "macros work when imported alone":
  ```julia
  m = Module(:_testhelper_MacroOnly)
  Core.eval(m, :(using EnzymeRates: @enzyme_mechanism, @enzyme_reaction))
  em = Core.eval(m, :(@enzyme_mechanism begin
      substrates: S
      products: P
      steps: begin
          E + S ⇌ E(S)
          E(S) <--> E(P)
          E(P) ⇌ E + P
      end
  end))
  @test em isa EnzymeRates.EnzymeMechanism
  ```
  plus the same for `@enzyme_reaction`, and "parse_mechanism reads show back": for every
  mechanism in `test/fixtures/canonical_text_golden.txt`,
  `sprint(show, EnzymeRates.DSL.parse_mechanism(text)) == text`.
- [ ] **Step 2: Run.** The round trip should already pass under the old name; the new name
  FAILs.
- [ ] **Step 3: Implement the move.**
- [ ] **Step 4: Focused DSL, Mechanisms, identify tests; boundary test; full suite.** Re-point
  `ER._parse_mechanism` in `$S/dm/census/*` and `$S/dm/bench_guard.jl`.
- [ ] **Step 5: Census, benchmarks, commit.** "DSL module: the reader and the macros".

### Task 7: `mwc_tags`

**Files:**
- Modify: every occurrence in `src/` of the field and accessor `cat_allo_states`, the
  `RegulatorySite` field and accessor `allo_states`, and the `RawSite`/`RawAllostericMechanism`
  fields of those names; docstrings that name them; `src/Mechanisms/Mechanisms.jl` (`API`).
- Test: every test using those names (`git grep -n "cat_allo_states\|allo_states"`).

**Interfaces:**
- Produces: field `mwc_tags` on `AllostericMechanism`, `RegulatorySite`, `RawSite` and
  `RawAllostericMechanism`; accessor `mwc_tags(am)`/`mwc_tags(s)` in place of
  `cat_allo_states(am)`/`allo_states(s)` (the accessors die in Task 11). Values unchanged.

- [ ] **Step 1:** Rename in the tests first. Run `test/Mechanisms/test_types.jl`: FAIL.
- [ ] **Step 2:** Rename in `src/`. Error messages that say "allo state" say "MWC tag".
- [ ] **Step 3: Full suite; census; commit.** "Rename the allosteric tags mwc_tags".

## Stage 2: The storage port

### Task 8: Nodes and steps as the mechanism's storage

Use the most capable model. This is the rewrite the spec's §3 authorizes.

**Files:**
- Modify: `src/Mechanisms/types.jl` (from `# ─── The reaction graph` to the end),
  `src/Mechanisms/accessors.jl`, `src/Mechanisms/write.jl`, `src/Mechanisms/names.jl`,
  `src/Mechanisms/text.jl`, `src/Mechanisms/Mechanisms.jl`.
- Test: `test/Mechanisms/*`, and any test reading `m.forms`, `m.edges`, `m.table`,
  `m.conformations` or `Counts` (`git grep -n`).

**Interfaces:**
- Produces (spec §3): `Node(conformation::Symbol, bound::Vector{Met}, added::Vector{Symbol},
  subtracted::Vector{Symbol})`, `Step(from::Int, to::Int, taken_up::Vector{Met},
  given_off::Vector{Met}, rapid_equilibrium::Bool, group::Int)`, `Mechanism` with fields
  `substrates`, `products`, `nodes`, `steps` (vectors sorted as today). Each `Vector{Met}` and
  `Vector{Symbol}` field is sorted. `==` and `hash` are field-wise. `Node`, `Step` join
  `Mechanisms.API`.
- Canonical order must equal today's: `isless(::Node)` compares `(conformation, bound, added,
  subtracted)` lexicographically on the sorted vectors (a prefix first), which matches today's
  `isless(::Form)`/`_seq_less`; free E first; steps by `(from, to, taken_up, given_off,
  rapid_equilibrium)`; groups numbered by first step. Orientation (`_orientation`,
  `_role_masks`, `_align_tied_groups!`) keeps its rules on the new types.
- In this task every function of `Mechanisms.API` keeps its signature and meaning. `Metabolites`
  becomes a wrapper of one sorted `Vector{Met}` with the same functions (`copies(s, name)`,
  `metabolite_names`, `modifier_copies`, `isempty`, `length`, `∪`, `∩`, `==`, `hash`); the
  accessors return it without copying. `Raw*` types stay as the constructors' input. Deleted:
  `Counts`, `Form`, `Edge`, `copies(::UInt128, i)`, `entries`, `total`, `_seq_less`.
- `Sig` follows the new storage; `compile_mechanism(Mechanism(em)) === em` and
  `Mechanism(compile_mechanism(m)) == m` hold.

- [ ] **Step 1: Pin the invariants.** Before touching `src/`, add to
  `test/Mechanisms/test_reaction_graph.jl` "storage round trips": for every mechanism in the
  two golden fixture files, `sprint(show, m)` equals its line, `Mechanism(compile_mechanism(m))
  == m`, and `hash(parse(text)) == hash(parse(text))` across two parses. Run: PASS on the
  parent (they pin today's behaviour).
- [ ] **Step 2: Write the new-storage tests:** "nodes and steps are the storage" — for
  `E + S ⇌ E(S); E(S) <--> E(P); E(P) ⇌ E + P`: `m.nodes[1] == Node(:E, Met[], Symbol[],
  Symbol[])`, `m.steps[1].taken_up == [Met(:S, false)]`, `m.steps[1].rapid_equilibrium`, and
  `length(m.nodes) == 3`. Run: FAIL.
- [ ] **Step 3: Port** the storage, canonicalization, checks, `Sig`, accessors, writers,
  names and printer. Migrate tests that read the deleted fields to the new ones.
- [ ] **Step 4: Run all `test/Mechanisms/*`, then the full suite.** PASS, pristine; golden
  files unchanged.
- [ ] **Step 5: Census.** `RESULT: PASS`, zero differences.
- [ ] **Step 6: Performance and memory checkpoint** (spec §7). Paired ABBA against the parent
  commit for all three benchmarks. Memory: a script `$S/dm2/memory.jl` that parses every
  mechanism of `DM_base/mechs.tsv` with `parse_mechanism` and prints
  `Base.summarysize(vector_of_mechanisms)`; run it on parent and port. Record both in
  `$S/dm/gates.md`. If any benchmark is more than 5% slower: profile, fix hot spots within
  the spec's design, re-measure; if it is still more than 5% slower, report BLOCKED with the
  numbers (the controller stops the storage-dependent tasks for Denis).
- [ ] **Step 7: Commit.** "Store a mechanism as nodes and steps".

## Stage 3: Read with dots, delete the read accessors

### Task 9: `step_kind`, `groups`, and Derivation reads with dots

**Files:**
- Modify: `src/Mechanisms/accessors.jl` (`step_kind` four values; `groups(m)`),
  `src/Mechanisms/Mechanisms.jl` (`API` gains `:groups`), every `src/Derivation/*.jl` file
  that calls a read accessor; every caller of `is_binding` in any module.
- Test: `test/Mechanisms/test_graph_accessors.jl`, Derivation tests as needed.

**Interfaces:**
- `step_kind(m, k)` returns `:isomerization`, `:binding` (a plain binding), `:fused_binding`
  (takes up one metabolite, gives off none, but its two nodes differ by more than it), or
  `:chemistry`. `is_binding(m, k)` callers become
  `step_kind(m, k) in (:binding, :fused_binding)`; `is_binding` is deleted.
- `groups(m)::Vector{Vector{Int}}`: the steps of each kinetic group, in step order.
- Derivation reads stored data with dots: `m.steps[k].from`, `.taken_up`, `m.nodes[x].bound`,
  `catalytic(m).substrates`, `am.mwc_tags[g]`, `site.ligands`… Multiset reads become vector
  code: `copies(v, n)` → `count(x -> x.name === n, v)`, `metabolite_names(v)` →
  `[x.name for x in v]`, `modifier_copies(v)` → `filter(x -> x.modifier, v)`.

- [ ] **Step 1: Failing tests** in `test_graph_accessors.jl`: one step of each kind written
  inline (an isomerization, a plain binding, a fused binding `E(A) + B <--> E(P, Q)`, a
  chemistry step), asserting `step_kind`; and `groups(m)` for a mechanism with a two-step group.
- [ ] **Step 2: Run.** FAIL.
- [ ] **Step 3: Implement**, then migrate Derivation file by file, running its focused tests
  after each file.
- [ ] **Step 4: Full suite.**
- [ ] **Step 5: Census, bench_guard and bench_pass1; commit.** "Derivation reads mechanisms
  with dots; step_kind has four kinds".

### Task 10: Enumeration, Fitting and identify read with dots

**Files:**
- Modify: `src/Enumeration/*.jl`, `src/Fitting/fitting.jl`, `src/identify_rate_equation.jl`,
  `src/DSL/*` — every read accessor call becomes a field read, as in Task 9.
- Test: as needed.

- [ ] **Step 1:** No new behaviour; the existing tests are the regression. Run
  `test/Enumeration/test_mechanism_enumeration.jl` on the parent and record its count.
- [ ] **Step 2: Migrate** file by file, running the focused enumeration test after each.
- [ ] **Step 3: Full suite; census; bench_enum; commit.** "Enumeration and identify read
  mechanisms with dots".

### Task 11: Delete the read accessors and `Metabolites`

**Files:**
- Modify: `src/Mechanisms/accessors.jl`, `src/Mechanisms/*.jl` internals that still call the
  accessors, `src/Mechanisms/Mechanisms.jl` (`API`), `test/test_boundaries.jl` (rule 3: drop
  `Form`, `Edge`, `Metabolites` from the banned types; `Node`, `Step`, `Mechanism`,
  `AllostericMechanism`, `RegulatorySite`, `Met`, `EnzymeReaction`, `Reactant`, `Regulator`
  fields are public).
- Test: every test calling a deleted accessor.

**Interfaces:**
- Deleted: `species`, `bound`, `conformation`, `residual`, `steps` (the function; `m.steps`
  stays), `ends`, `taken_up`, `given_off`, `is_rapid_equilibrium`, `rapid_equilibrium_steps`,
  `group`, `members`, `ngroups`, `tag`, `mwc_tags` (the accessor), `catalytic_multiplicity`,
  `regulatory_sites`, `ligands`, `multiplicity`, `Metabolites`, `copies`, `metabolite_names`,
  `modifier_copies`, `substrates`/`products` on mechanisms. `catalytic`, `metabolites`,
  `modifiers`, `step_kind`, `groups` stay.

- [ ] **Step 1:** Remove the names from `Mechanisms.API`; run the boundary test: FAIL wherever
  a caller remains. Fix callers (src and tests) until green.
- [ ] **Step 2: Delete the definitions.** Full suite.
- [ ] **Step 3: Census; benchmarks; commit.** "Mechanisms is read with dots: delete the read
  accessors".

## Stage 4: Each function where its job is

### Task 12: Derivation owns the naming and the allosteric semantics

**Files:**
- Create: `src/Derivation/names.jl` (the parameter half of `src/Mechanisms/names.jl`),
  `src/Derivation/conformations.jl` (`excludes`, `exclusion_components`,
  `conformation_steps`), `src/Derivation/segments.jl` (`present_species`, `segment_walk`,
  `re_segments`, `segment_graph`, `has_unbalanced_cycle`, `has_turnover_cycle`,
  `runs_toward_products`, `_name_index`).
- Modify: `src/Mechanisms/names.jl` keeps only what the constructor checks and the printer use
  (`_form_name`, `_sides` behind the public `step_sides`); `src/Mechanisms/types.jl` (remove
  the moved functions); both `API` lists; every caller (`Mechanisms.x` → `Derivation.x` in
  Enumeration; bare in Derivation).
- Test: move the chokepoint testset ("chokepoint: no Symbol(\"[KkVL]...\") outside
  parameter-name renderers") from `test/Mechanisms/test_types.jl` to
  `test/Derivation/test_names.jl` and point its scan at `src/Derivation/names.jl`; move the
  tests of the moved functions to `test/Derivation/`; update CLAUDE.md's chokepoint paragraph.

**Interfaces:**
- `Mechanisms.API` gains `:step_sides` (`step_sides(m, k) -> (String, String)`, today's
  `_sides`); loses `constant_columns`, `regulator_constant`, `name`, `excludes`,
  `exclusion_components`, `conformation_steps`, `has_turnover_cycle`, `has_unbalanced_cycle`,
  `re_segments`, `segment_graph`, `segment_walk`, `present_species`, `runs_toward_products`.
- `Derivation.API` gains `:conformation_steps`, `:excludes` (the `excludes(tag, c)` method;
  the `(am, k, c)` method stays as a Derivation method), `:has_turnover_cycle`,
  `:has_unbalanced_cycle`, `:re_segments`, `:segment_graph`.
- The moved functions read the mechanism with dots; behaviour unchanged.

- [ ] **Step 1:** Move the tests first (they now call `EnzymeRates.Derivation.x`); run: FAIL.
- [ ] **Step 2:** Move the code; fix callers; update both API lists and the boundary test.
- [ ] **Step 3: Focused Derivation and Enumeration tests; full suite.**
- [ ] **Step 4: Census; benchmarks; commit.** "Derivation owns the parameter names and the
  allosteric conformations".

## Stage 5: Write with nodes and steps

### Task 13: Constructors from nodes and steps; the DSL writes them

**Files:**
- Modify: `src/Mechanisms/types.jl` (constructors), `src/DSL/read.jl`,
  `src/Derivation/admission.jl` (its `with_flags`/`with_tags` calls), `src/Mechanisms/write.jl`.
- Test: `test/Mechanisms/test_reaction_graph.jl`, `test/Mechanisms/test_types.jl`,
  `test/DSL/test_dsl.jl`.

**Interfaces:**
- Produces: `Mechanism(substrates::Vector{Symbol}, products::Vector{Symbol},
  nodes::Vector{Node}, steps::Vector{Step})` — nodes in any order, duplicates allowed; steps in
  any direction, `group` as labels 1, 2, … each used; canonicalizes and runs today's
  representation checks with today's messages.
- Produces: `AllostericMechanism(substrates, products, nodes, steps, mwc_tags::Vector{Symbol},
  catalytic_multiplicity::Int, regulatory_sites::Vector{RegulatorySite})` with
  `mwc_tags[label]` per written group label; and the composition constructor
  `AllostericMechanism(catalytic::Mechanism, mwc_tags, catalytic_multiplicity, sites)`
  (today's) over canonical groups.
- The "declared modifier no step binds" check moves to `DSL` with its message unchanged.
- Admission rebuilds with the constructors (copy `m.steps` with changed `rapid_equilibrium`;
  the composition constructor for tags).
- The `RawMechanism`/`RawAllostericMechanism` constructors of `Mechanism` and
  `AllostericMechanism` remain only until Task 14 removes their last caller.

- [ ] **Step 1: Failing tests** in `test_reaction_graph.jl`, "construction is canonical": the
  Michaelis–Menten mechanism built from nodes and steps (a) as canonical, (b) with nodes
  permuted and one duplicated, (c) with every step reversed (forms and sides swapped), (d)
  with group labels renumbered; all four `==` and equal `hash`, and `sprint(show, ·)` equal.
  Plus one allosteric case through the written-parts constructor.
- [ ] **Step 2: Run.** FAIL.
- [ ] **Step 3: Implement**; switch the DSL reader and admission.
- [ ] **Step 4: Golden round trip** (Task 6's test) and the full suite.
- [ ] **Step 5: Census; bench_guard; commit.** "Mechanisms are built from nodes and steps".

### Task 14: Enumeration writes nodes and steps; delete the written form

**Files:**
- Modify: `src/Enumeration/seeds.jl`, `src/Enumeration/moves.jl`, `src/Enumeration/filters.jl`;
  `src/Mechanisms/write.jl` (deleted unless something remains), `src/Mechanisms/types.jl`
  (`Raw*` types and constructors), `src/Mechanisms/Mechanisms.jl` (`API`),
  `test/test_boundaries.jl` (delete rule 3, `_testhelper_FIELD_ALLOWED`,
  `_testhelper_MECHANISM_TYPES` and their testsets; renumber the rules in comments and in
  CLAUDE.md's "Module boundaries" paragraph).
- Test: `test/Enumeration/test_mechanism_enumeration.jl`, `test/Mechanisms/test_write.jl`,
  `test/Mechanisms/test_types.jl`, `test/Mechanisms/test_reaction_graph.jl`, and the shared
  fixtures file — every `Raw*` use.

**Interfaces:**
- Seeds build `Node`/`Step` lists and call the constructors. Moves rebuild by copying
  `m.steps`/`m.nodes` with changes. The split move's group-to-tag mapping (today in
  `regroup`) becomes a private Enumeration helper.
- Deleted: `RawForm`, `RawStep`, `RawMechanism`, `RawSite`, `RawAllostericMechanism`, the
  `Raw*(m, …)` converters, `with_flags`, `regroup`, `with_tags`, `canonical_graph` (public).
  If the dead-end move or the seed bases need an unchecked build (a check fails on their
  intermediate graph, or the checks cost more than the 5% gate), keep one private unchecked
  constructor in `Mechanisms` reached through a public non-underscore name only if rule 4
  requires it, and say why in the report.
- `Mechanisms.API` after this task is exactly spec §3: the 12 types and `step_kind`, `groups`,
  `metabolites`, `modifiers`, `catalytic`, `compile_mechanism`, `form_text`, `step_text`,
  `step_sides`. A boundary testset pins this list.

- [ ] **Step 1:** Add the boundary testset "Mechanisms API is the data model" pinning the list
  above; run: FAIL.
- [ ] **Step 2: Migrate seeds**, focused enumeration test after; then **moves and filters**.
- [ ] **Step 3: Migrate the tests' `Raw*` uses** to `@enzyme_mechanism` fixtures or
  nodes/steps (CLAUDE.md's enumeration-test rules apply: fixtures inline, exact children).
- [ ] **Step 4: Delete** the written form and the writers; delete rule 3.
- [ ] **Step 5: Full suite; census; bench_enum and bench_guard; commit.** "Enumeration writes
  nodes and steps; delete the written form".

## Stage 6: Docs and final verification

### Task 15: The module map in the developer docs

**Files:**
- Modify: `docs/src/developer.md` (rewrite the "Modules" section and its "Accessors, views and
  writers" subsection; fix every later mention of deleted names), `docs/src/api.md`,
  `.claude/CLAUDE.md` ("API Design", "Module boundaries" with the rule count, focused-test
  paths, the chokepoint location), `ROADMAP.md` (retire entries this work closes),
  `README.md` if it shows a deleted name.

**Interfaces (spec §8):** the "Modules" section holds, in this order:
1. The diagram: a text diagram in a code block of `Data`, `Mechanisms`, `Derivation`,
   `Fitting`, `Enumeration`, `DSL` and the top level, each arrow labelled with the main calls
   (e.g. `Enumeration → Derivation: check_valid_rate_equation, parameter_count`).
2. One short paragraph per module: what it is for and what it hides.
3. One table per module: every `API` name, a one-line description, its users.
4. The mechanism data model: the field tables of spec §3 and the rule "stored data with dots,
   derived data with functions".
5. One request traced through the modules: `identify_rate_equation` from the reaction to the
   selected fit.

- [ ] **Step 1: Write the section.** Check every name against the `API` constants (a name in
  the docs that is not in an `API` list, or the reverse, is a defect).
- [ ] **Step 2: Docs build** with doctests (`julia --project=docs docs/make.jl`, or the
  `$S/dm/docsenv` copy the first split used). Expected: no errors, no doctest failures.
- [ ] **Step 3: Full suite on Julia 1.12 and 1.13** (`julia +1.13 --project …`). Both pristine.
- [ ] **Step 4: Final census** against `DM_base`: PASS. **Final benchmarks** paired against
  the commit before Task 1 (`1dfad6d3`); record in `gates.md`.
- [ ] **Step 5: Commit.** "Developer docs: the module map".
