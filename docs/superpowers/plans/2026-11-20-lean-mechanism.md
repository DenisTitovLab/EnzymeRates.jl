# Lean Mechanism Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make a mechanism one lean reaction graph with opaque form names, behind smaller
`Mechanisms`, `Derivation` and `Enumeration` interfaces, with `rate =
compile_rate_equation(m)`, `rate(concs, params, Keq)` as the one way to evaluate a law.

**Architecture:** The work goes inside-out. First the interface changes that today's data
model can carry (names, `Keq`, parameters, compiling, candidates). Then each module stops
reading form contents except through one chokepoint. Only then does the data model
switch, in three steps: plain values, one type and one macro, lean forms. The suite passes
after every task.

**Tech Stack:** Julia 1.13 (default `julia`), the package's own test suite, JET, Aqua.

**Spec:** `docs/superpowers/specs/2026-11-18-lean-mechanism-design.md`. It is the binding
authority. Read the decision a task names before starting the task.

**Maps of today's code** (written 2026-11-20 at commit 56bea2e6, outside the repo; line
numbers drift as tasks land, so check a line before relying on it):
`/home/denis.linux/claude-scratchpad/lean-mechanism-spike/plan-inventory/`
`mechanisms_dsl.md`, `derivation.md`, `enumeration_identify.md`, `fitting_tests_docs.md`.
Each ends with a "Traps" section. Read the map of every module a task touches.

## Global Constraints

- **Branch.** Work on `guard-rule-excluded-steps`, in place. Stage files by path; never
  `git add -A` (`simulations/` is untracked and not part of this work). Never skip a
  hook. Never push.
- **Commit trailers.** Every commit message ends with
  `Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_015DwkaD8jqpCL89XaikN6mF`.
- **One Julia process at a time** (7.7 GB, no swap). Every Julia command runs under
  `flock /home/denis.linux/claude-scratchpad/julia.lock`.
  - Focused run of one test file:
    `flock /home/denis.linux/claude-scratchpad/julia.lock julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<Module>/<file>.jl")'`
  - Full suite, in the background, about 10 minutes:
    `flock /home/denis.linux/claude-scratchpad/julia.lock timeout 1200 julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`
  - About one suite in nine fails the `loss!` timing pin from stack placement alone
    (`claude-scratchpad/loss-timing-investigation/report.md`). Rerun once before treating
    that one failure as real.
- **Every task ends with a green full suite**, and the task's report gives the pass
  count. Test output is pristine: an expected error or warning is captured and asserted.
- **No verdict changes** (spec, "Tests"). Every pin of a guard verdict, of an admitting
  reading and of an enumerator population passes with the same counts. A pin that would
  have to change stops that part of the work: report it, do not edit the test. The
  exceptions are the spec's: some fitted-parameter names (the constants `_priority_key`
  picks, and `mod` leaving a name), `E_total` and `Keq`, and the constants of the steps
  that decision 1 turns (32 isomerizations and at most three exchanges, in about 25
  fixtures). There a constant is renamed after the other direction and a
  rapid-equilibrium constant becomes its inverse; the law is the same function.
- **Planner's ruling on order.** The stored order of forms, steps and groups will change
  for some mechanisms (forms are ordered by name, no longer by contents). A pin on the
  *order* of a list (`pp[1:62]`, `got == [a, b]`, a `_Law ==` that differs by form
  order) may be rewritten when the set and the count are unchanged. Every such pin goes
  in the report, with the reason. This is the planner's call and the maintainer reads the
  list at the end.
- **A failing test is never deleted or weakened.** A test of a removed feature
  (`E_total`, a parameter kind, a compiled type, `@allosteric_mechanism`, call notation)
  goes with the feature, and the report lists it.
- **Oracles stay independent.** A test helper that models the guard or a rate law keeps
  its own statement of the rule and never calls the package's.
- **Hash.** `hash` of a mechanism is computed from its fields, equal for every writing
  and in every process. The guard's limit test seeds its draws with it
  (`src/Derivation/guard.jl`, `hash(i, hash(mech))`), and identify and the seed search
  drop duplicates by hash.
- **The free enzyme is form 1**, in every intermediate state. A rapid-equilibrium
  segment is rooted at its smallest form index, and the derivation assumes the first
  tree is free E's.
- **Speed.** The compiled rate stays allocation-free and under
  `50 ns + 0.5 ns × denominator terms`; `Fitting.loss!` keeps its bounds; the
  compile-budget and type-stability tests pass. Stop and report before a change that
  would make the hot path allocate or slow down.
- **Module boundaries.** A task that changes a `public` statement, an export, a module
  docstring or a cross-module call updates, in the same commit: the interface tables and
  the diagram in `docs/src/developer.md`, the pinned name sets in
  `test/test_boundaries.jl`, and the module map (`flock … julia test/update_module_map.jl`).
  The other docs pages wait for Task 16, so the docs build may be red in between.
- **Public names** added by this plan are the spec's, already approved:
  `CatalyticScheme`, `Allostery`, `steps_in_active_conformation`,
  `steps_in_inactive_conformation` (Mechanisms); `compile_rate_equation` (exported),
  `parameter_units` (Derivation); `CandidateMechanism`, `candidate_mechanisms`
  (Enumeration). No other public name is added. A transitional private helper is fine.
- **Style.** `.claude/CLAUDE.md` binds: deep functions, no pass-through helpers, 92
  columns, two `ABOUTME:` lines per file, `_testhelper_` prefix under `test/`, names and
  comments describe the code as it is.
- **Stop and report** when: a guard verdict or an enumerator population would change; a
  law fails its oracle; a hand-written `MECHANISM_TEST_SPECS` entry or an identify
  test's true mechanism is rejected or no longer recovered; a speed or compile bound
  cannot be met; or the task cannot be done without a public name this plan does not
  list.
- **Version** stays 0.11.0.

## Review Focus

1. A mechanism written with no form named `E` (the free enzyme called `Efree`). Expected:
   an error from the macro and from the constructor that says the free enzyme is the
   form named `E`. Test in Task 15.
2. A kinetic group of isomerizations with one member written the other way round.
   Expected: it builds, it is `!=` the aligned mechanism, and reversing every member of
   the group gives a mechanism `==` the aligned one. Test in Task 9.
3. A `params` that still holds `Keq` or `E_total`, as a fit saved by an earlier version
   does. Expected: `rate`, `test_rate_equation` and the plot ignore the extra names and
   use the `Keq` argument. Test in Task 4.
4. `oligomeric_state:` or `regulatory_site(...)` on a mechanism with no allosteric
   states, or states on some statements only. Expected: an error in the user's terms.
   Test in Task 13.
5. A misspelt form name (`EAB` once written `EBA`). Expected: when the new form hangs
   off the graph with no way back, the error names it; the macro's docstring says that a
   misspelling can also pass. Test in Task 15.

---

### Task 1: Baseline and neutral keys

Nothing in the repository changes. The task records what later tasks must reproduce, in
a form that survives the change of notation.

**Files:**
- Create: `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/neutral_key.jl`
- Create: `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/neutral_before.jl`
- Create: `…/gate/neutral_before.jls`, `…/gate/neutral_before.out`,
  `…/gate/baseline-suite-56bea2e6.log`

**Interfaces:**
- Produces: `neutral_key(m) -> String`, defined in `neutral_key.jl` for today's
  `Mechanism` and `AllostericMechanism`. Task 15 adds the method for the lean
  `Mechanism`, and the two must give equal strings for the same mechanism.
- Produces: `neutral_before.jls`, a `Dict` with `:populations` (name => sorted vector of
  `(key, verdict)`) and `:fixtures` (file => sorted vector of keys). Later tasks add
  fixtures, so fixtures are compared per file as multisets: every key recorded here must
  still be present.

- [ ] **Step 1: Run the full suite** and save its log. Expected: 0 failures. Record the
  pass count.
- [ ] **Step 2: Write `neutral_key(m)`.** The key is the SHA-256 of a canonical text
  built from, in this order: the sorted substrates; the sorted products; the kinetic
  groups as a sorted list, each group a sorted list of steps with its allosteric state
  (`:None` when plain); the oligomeric state; the regulatory sites as a sorted list of
  `(ligands, multiplicity, states)`. A step is written without direction: its two sides
  sorted, a side being the form's name (today `Mechanisms._form_name`) and the sorted
  bare names of the metabolites free on that side, followed by its rapid-equilibrium
  flag. A group holding two or more isomerizations adds their stored directions as
  `(from name, to name)` pairs, taking the lexicographically smaller of the pattern and
  the pattern with every pair reversed. A modifier copy is written as its reactant's
  bare name.
- [ ] **Step 3: Write `neutral_before.jl`.** Populations: for each reaction of
  `gate/populations.jl` and of `orientation/check2.jl`, the seeds
  (`Enumeration.init_mechanisms`) and one expansion; for ter-ter, the 3,000 seeds
  `check2.jl` samples. For each mechanism store `neutral_key` and the guard's verdict
  (`:valid` or the fault's name). Fixtures: every `@enzyme_mechanism` and
  `@allosteric_mechanism` block under `test/` and `docs/src/`, found as
  `orientation/check3.jl` finds them, grouped by file; a block that does not evaluate is
  counted per file as `:not_evaluated`.
- [ ] **Step 4: Run it** and check the output: the population sizes equal those in
  `orientation/check2.out`, and at most 368 fixtures are `:not_evaluated`.
- [ ] **Step 5: Report.** No commit. The report gives the suite's pass count, the
  population sizes, the fixture count and the path of each file.

### Task 2: Remove `excludes`; make `segment_walk` and `present_species` private

Spec decision 5.

**Files:**
- Modify: `src/Mechanisms/types.jl` (`excludes`, the two constructors' masks),
  `src/Mechanisms/segments.jl`, `src/Mechanisms/Mechanisms.jl`
- Modify: `src/Derivation/thermo.jl`, `src/Derivation/names.jl`,
  `src/Derivation/guard.jl`, `src/Derivation/conformations.jl`,
  `src/Enumeration/moves.jl`
- Modify: `test/test_boundaries.jl`, `test/Mechanisms/test_segments.jl`,
  `docs/src/developer.md`, `.claude/CLAUDE.md` (module map, regenerated)

**Interfaces:**
- Produces: `Mechanisms` no longer has `excludes`, `segment_walk` or `present_species`
  among its public names (25 → 22). `Mechanisms._segment_walk` and
  `Mechanisms._present_species` stay for `rapid_equilibrium_segments`.

- [ ] **Step 1: Pin.** In `test/Derivation/test_admission.jl` add a testset
  "the guard's fast segments": for three mechanisms written inline (a plain ordered
  bi-bi, an allosteric uni-uni with an `:OnlyA` binding, a ping-pong), assert the
  `segment`, `segments`, `offset` and `cycle` that `Derivation._fast_segments` returns
  today, as literal values. Run it: PASS.
- [ ] **Step 2: Give the guard its own walk.** `_fast_segments(net)` in
  `src/Derivation/guard.jl` walks `net` itself and returns the same four values. The
  call to `Mechanisms.present_species` in `guard.jl` becomes Derivation's own loop.
- [ ] **Step 3: Remove `excludes`.** Each of its five callers compares the state:
  `tag === (inactive ? :OnlyA : :OnlyI)`. The constructor's two masks do the same.
- [ ] **Step 4: Rename** `segment_walk` to `_segment_walk` and `present_species` to
  `_present_species` in `Mechanisms`; update the `public` statement, the boundary pins,
  the developer page and the module map.
- [ ] **Step 5: Run** the focused files `test/Derivation/test_admission.jl`,
  `test/Mechanisms/test_segments.jl`, `test/test_boundaries.jl`, then the full suite.
- [ ] **Step 6: Commit** "Remove excludes; make the segment walk private to Mechanisms".

### Task 3: `oligomeric_state` and `allosteric_states`

Spec decision 1 (the renames). Mechanical; no behaviour changes.

**Files:**
- Modify: every file under `src/`, `test/` and `docs/src/developer.md` that names
  `catalytic_multiplicity`, `allowed_catalytic_multiplicities` or `mwc_tags` (counts per
  file in `fitting_tests_docs.md` section 3; about 54 lines in `src/`, 689 in `test/`).

**Interfaces:**
- Produces: the field `oligomeric_state` (was `catalytic_multiplicity`), the
  `EnzymeReaction` field `allowed_oligomeric_states` (was
  `allowed_catalytic_multiplicities`), the field `allosteric_states` on the allosteric
  mechanism and on `RegulatorySite` (was `mwc_tags`), and the mechanism macro's keyword
  `oligomeric_state:` (was `catalytic_multiplicity:`). `show` prints the new keyword and
  `parse_mechanism` reads it.

- [ ] **Step 1: Rename** in `src/`, with word-boundary replacement, then read every
  changed line: local variables named `mults`, `cn` or `tags` stay; docstrings and error
  messages follow the fields.
- [ ] **Step 2: Rename** in `test/` and in `docs/src/developer.md`.
- [ ] **Step 3: Regenerate** `test/reference/allosteric_golden_reference.txt` only if it
  holds the old keyword, and say so in the report with the diff's line count.
- [ ] **Step 4: Run** the full suite. Expected: the pass count of Task 1.
- [ ] **Step 5: Commit** "Rename catalytic_multiplicity to oligomeric_state and mwc_tags
  to allosteric_states".

### Task 4: `Keq` as an argument, no `E_total`, one parameter list

Spec decision 7, "The parameter interface", except `parameter_units` (Task 5) and the
compiled rate (Task 6).

**Files:**
- Modify: `src/Derivation/rate_eq_derivation.jl`, `codegen.jl`, `print.jl`, `thermo.jl`,
  `Derivation.jl`; `src/Fitting/fit_rate_equation.jl`; `src/Plotting/plot_rate_equation_on_data.jl`;
  `src/Data/rate_data.jl`; `src/identify_rate_equation.jl`; `src/EnzymeRates.jl`
- Modify: the tests and helpers that build `params` or pass `Keq`
  (`test/mechanism_definitions_for_test_enzyme_derivation.jl`, `test/Derivation/*.jl`,
  `test/Fitting/*.jl`, `test/test_identify_rate_equation.jl`, `test/test_type_stability.jl`,
  `test/test_compile_budget.jl`, `test/Plotting/`), `docs/src/developer.md`, the module map

**Interfaces:**
- Produces, until Task 6 renames the first:
  - `rate_equation(m, concs::NamedTuple, params::NamedTuple, Keq::Real)`: the rate per
    unit enzyme. `params` holds the names of `parameters(m)`; extra names are ignored.
  - `parameters(m) -> Tuple{Vararg{Symbol}}`: the constants a fit estimates, in today's
    `FittedParameters` order. No second argument.
  - `dependent_parameter_values(m, params::NamedTuple, Keq::Real) -> NamedTuple`: the
    dependent constants, by name.
  - `FittingProblem(m, table, Keq::Real; absolute = false, reference = nothing)`
  - `IdentifyRateEquationProblem(reaction, table, Keq::Real; absolute = false, reference = nothing)`
  - `test_rate_equation(m, params::NamedTuple, table, Keq::Real; absolute = false)`
  - `plot_rate_equation_on_data(m, params, table, Keq::Real; …)` and
    `plot_rate_equation_on_data(f::Function, params, table; …)`
  - `Data.RateData(table, metabolite_names, Keq::Real; absolute, reference, owner)`
  - A fit's `params` holds exactly `parameters(m)`.
- Removed: `IndependentParameters`, `FittedParameters`, `AllParameters` (exports and
  types), `E_total` everywhere.

- [ ] **Step 1: Write the failing tests** in `test/Derivation/test_rate_eq_derivation.jl`,
  testset "Keq is an argument":

```julia
m = @enzyme_mechanism begin
    substrates: S
    products: P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end
names = parameters(m)
@test :Keq ∉ names && :E_total ∉ names
params = NamedTuple{names}(Tuple(rand(MersenneTwister(1), length(names)) .+ 0.5))
concs = (S = 0.7, P = 0.3)
v = rate_equation(m, concs, params, 2.0)
# The law gives the rate per unit enzyme, and extra names in params are ignored.
@test rate_equation(m, concs, merge(params, (Keq = 99.0, E_total = 99.0)), 2.0) == v
@test rate_equation(m, (S = 0.3 / 2.0, P = 0.3), params, 2.0) ≈ 0 atol = 1e-12
dep = dependent_parameter_values(m, params, 2.0)
@test dep isa NamedTuple && isdisjoint(keys(dep), names)
@test !occursin("E_total", rate_equation_string(m))
@test occursin("Keq", rate_equation_string(m))
```

  Add to `test/Fitting/test_fitting.jl`, testset "a fit's params hold the fitted names":
  `fit = fit_rate_equation(FittingProblem(m, data, 2.0)); @test keys(fit.params) == parameters(m)`,
  with `data` the table that file's first fit uses; and
  `@test test_rate_equation(m, merge(fit.params, (Keq = 99.0,)), data, 2.0) == test_rate_equation(m, fit.params, data, 2.0)`.
- [ ] **Step 2: Run them.** Expected: FAIL (method errors).
- [ ] **Step 3: Derivation.** `Keq` is read from the argument in the generated body, not
  from `params` by name; `E_total` leaves the four sites the map lists (derivation.md
  section 3); `rate_equation_string` shows `Keq` as an input of its own. The three kinds
  go; `parameters(m)` is the fitted list.
- [ ] **Step 4: Fitting, Data, Plotting, identify.** `loss!` passes `fp.rate_data.Keq`
  as the argument and merges nothing into `params`. `_fit_keywords` in identify passes
  `Keq` positionally. The cross-validation path passes `Keq` to `test_rate_equation`.
- [ ] **Step 5: Tests and helpers.** Helpers that built `(…, Keq, E_total)` tuples build
  the fitted names and pass `Keq` beside them. The parts of
  `_testhelper_test_rescaling` and `_testhelper_name_units` that scale `E_total` go,
  and the report lists them. A test that used `AllParameters` to name the dependent
  constants uses `keys(dependent_parameter_values(m, params, Keq))`.
- [ ] **Step 6: Run** `test/Derivation/test_rate_eq_derivation.jl` (the speed and
  allocation gate is in it), `test/Fitting/test_fitting.jl`, then the full suite.
- [ ] **Step 7: Commit** "Make Keq an argument; remove E_total and the parameter kinds".

### Task 5: `parameter_units`

Spec decision 7, the table of units.

**Files:**
- Modify: `src/Derivation/rate_eq_derivation.jl`, `names.jl`, `Derivation.jl`;
  `src/Fitting/fit_rate_equation.jl`
- Test: `test/Derivation/test_names.jl`, `test/Fitting/test_fitting_scale.jl`,
  `test/test_type_stability.jl`, `test/test_boundaries.jl`

**Interfaces:**
- Consumes: `parameters(m)` (Task 4).
- Produces: `Derivation.parameter_units(m) -> (; rate::Vector{Int}, concentrations::Matrix{Int})`,
  public, not exported. `rate[i]` is the power of the rate in the units of
  `parameters(m)[i]`; `concentrations[i, j]` the power of `metabolites(m)[j]`. The rows
  are the spec's table.
- Removed from `Derivation`: `parameter_info`, `name`, `StepParameter`,
  `RegulatorParameter`, `AllostericParameter`. Fitting reads no step and no
  `step_kinds`.

- [ ] **Step 1: Pin first.** In `test/Fitting/test_fitting_scale.jl` add
  "starting values are unchanged": for every `MECHANISM_TEST_SPECS` entry, with the
  typical rate 3.0 and typical concentrations `10.0^-j` for metabolite `j`, record
  today's log centres (`Fitting._log_centres` or its equivalent; see the map) into a
  literal `Dict` keyed by the spec's name. Run: PASS.
- [ ] **Step 2: Write the failing test** of `parameter_units` in
  `test/Derivation/test_names.jl`: for the uni-uni mechanism of Task 4 with a
  steady-state chemistry step and rapid-equilibrium bindings, the units of each named
  constant equal the spec's table row; and for every `MECHANISM_TEST_SPECS` entry,
  `r^rate × Π s_j^power` reproduces the centres pinned in Step 1.
- [ ] **Step 3: Implement** `parameter_units`, computed from the steps and the
  regulatory sites. Fitting's centres, the rescaling of relative data and
  `_steady_state_constant_names` read it.
- [ ] **Step 4: Remove** the five names; update the type-stability rows that named
  them, the boundary pins, the developer page and the module map.
- [ ] **Step 4b: `fitted_parameter_count`.** Time it against `length(parameters(m))` on
  the bi-bi seeds. If the second is no slower, replace the first everywhere and remove
  the public name; otherwise keep it. Report the two times.
- [ ] **Step 5: Run** the two focused files, `test/test_type_stability.jl`, then the
  full suite.
- [ ] **Step 6: Commit** "Give Fitting parameter_units in place of parameter_info".

### Task 6: `compile_rate_equation`

Spec decision 3, on today's mechanism types.

**Files:**
- Modify: `src/Derivation/codegen.jl`, `rate_eq_derivation.jl`, `Derivation.jl`;
  `src/Fitting/fit_rate_equation.jl`; `src/Plotting/plot_rate_equation_on_data.jl`;
  `src/identify_rate_equation.jl`; `src/EnzymeRates.jl`
- Modify: every test that calls `rate_equation(` (counts per file in
  `fitting_tests_docs.md` section 3), `test/test_compile_budget.jl`,
  `test/test_type_stability.jl`, `test/test_boundaries.jl`, `docs/src/developer.md`

**Interfaces:**
- Consumes: `rate_equation(m, concs, params, Keq)` (Task 4).
- Produces: `rate = compile_rate_equation(m)`, exported. `rate(concs::NamedTuple,
  params::NamedTuple, Keq::Real)` returns the rate per unit enzyme, allocation-free.
  `parameters(rate)` and `Mechanisms.metabolites(rate)` return tuples of Symbols that
  are compile-time constants. `compile_rate_equation` accepts every mechanism form the
  package has at that commit.
- Produces: `FittingProblem` holds both the mechanism and its compiled rate.
- Removed: `rate_equation`, the function and the export. The callable's type is private
  to `Derivation`.

- [ ] **Step 1: Write the failing test** in `test/Derivation/test_rate_eq_derivation.jl`,
  "compile_rate_equation": with `m`, `params`, `concs` as in Task 4,

```julia
rate = compile_rate_equation(m)
@test rate(concs, params, 2.0) == rate_equation(m, concs, params, 2.0)
@test parameters(rate) == parameters(m)
@test metabolites(rate) == Tuple(metabolites(m))
@test (@allocated rate(concs, params, 2.0)) == 0
@test compile_rate_equation(m) === rate
```

- [ ] **Step 2: Run it.** Expected: FAIL, `compile_rate_equation` not defined.
- [ ] **Step 3: Implement** the callable and its generated call; `loss!` calls it and
  reads the two tuples from it.
- [ ] **Step 4: Migrate** every caller in `src/` and `test/`, then delete
  `rate_equation` and the first assertion of Step 1. A test that evaluates one law many
  times compiles once, outside its loop.
- [ ] **Step 5: Move the gates.** The speed and allocation gate, the compile-budget
  test and the type-stability rows measure `rate`.
- [ ] **Step 6: Run** `test/Derivation/test_rate_eq_derivation.jl`,
  `test/test_compile_budget.jl`, `test/test_type_stability.jl`, then the full suite.
- [ ] **Step 7: Commit** "Evaluate a law through compile_rate_equation".

### Task 7: Derivation reads the graph only

Spec decision 5, the first three bullets. After this task no file under
`src/Derivation/` reads a form's contents, `step_kinds`, `toward_products`,
`step_sides`, `modifiers` or `groups`, or orders forms or steps with `isless`.

**Files:**
- Modify: `src/Derivation/names.jl`, `thermo.jl`, `engine.jl`, `key.jl`, `print.jl`,
  `guard.jl`, `admission.jl`, `messages.jl`, `conformations.jl`
- Modify: `src/Mechanisms/types.jl` (one stored field, below)
- Test: `test/Derivation/test_names.jl`, `test_key.jl`, `test_rate_eq_derivation.jl`

**Interfaces:**
- Produces: today's `Mechanism` gains the stored field `enzyme_forms::Vector{Symbol}`,
  the form names (`_form_name`), one per node. It is the field the lean graph keeps.
  `Derivation` names constants from `enzyme_forms` and the bare metabolite names.
- Produces: `check_valid_rate_equation` rejects a mechanism whose constants would share
  a name (today the solve throws it later; derivation.md, "Five facts", 1).
- Changes named by the spec: a constant's name writes a modifier copy as the bare
  reactant name (`mod` leaves the metabolite part); the constants `_priority_key` picks
  as dependent may differ.

- [ ] **Step 1: Pin first.** Add `test/Derivation/test_names.jl` testset "fitted names
  before the graph-only rule" that records `parameters(m)` for every
  `MECHANISM_TEST_SPECS` entry as literal tuples. Run: PASS. This testset is the record
  of which names move; Step 5 rewrites it.
- [ ] **Step 2: Write the failing test** "constants that would share a name are
  rejected by the guard": find the existing test of the error that
  `src/Derivation/thermo.jl` throws for two constants with one name, and assert for its
  mechanism that `check_valid_rate_equation(Bool, m)` is false. If no mechanism can
  reach that error, report it and delete the check.
- [ ] **Step 3: Replace each read**, in the order of derivation.md section 5: kinetic
  groups from the steps' `group`; a step's shape (no metabolite, one taken up, an
  exchange) from its lists; names from `enzyme_forms`.
- [ ] **Step 4: The rule in `_priority_key`.** It reads the step's shape, its
  rapid-equilibrium flag, its kinetic group and its position in the stored order. Choose
  the rule that changes the fewest tuples of Step 1, and measure: the report gives the
  number of `MECHANISM_TEST_SPECS` entries and of the bi-bi seeds
  (`Enumeration.init_mechanisms`) whose `parameters` change.
- [ ] **Step 5: Rewrite the pins** that hold a moved name: the Step 1 testset, and any
  pin elsewhere. Each law still passes its oracle; no oracle changes.
- [ ] **Step 6: Run** every file in `test/Derivation/`, then the full suite.
- [ ] **Step 7: Commit** "Derive names and dependent constants from the graph alone".

### Task 8: A conformation's steps as functions

Spec decision 5, "A conformation's steps".

**Files:**
- Modify: `src/Mechanisms/types.jl`, `accessors.jl`, `Mechanisms.jl`
- Modify: the readers of the four masks in `src/Derivation/` and `src/Enumeration/`
  (mechanisms_dsl.md section 6)
- Test: `test/Mechanisms/test_types.jl`, `test/test_boundaries.jl`,
  `test/test_type_stability.jl`

**Interfaces:**
- Produces: `Mechanisms.steps_in_active_conformation(m) -> BitVector` and
  `Mechanisms.steps_in_inactive_conformation(m) -> BitVector`, one flag per step, for
  every mechanism form the package has at that commit. A plain mechanism gives all
  true for both.
- Removed: the stored `active_steps`, `inactive_steps`, `excluded_active`,
  `excluded_inactive`. A caller that needs the steps a conformation excludes by state
  builds that mask by comparing the states.

- [ ] **Step 1: Write the failing test** "a conformation's steps": for the spec's
  example (`E + S ⇌ E(S) :: OnlyA`, `E(S) <--> E(P) :: OnlyA`, `E + P ⇌ E(P) ::
  EqualAI`) the inactive result keeps only the P binding and the active result keeps
  all three; for a mechanism with an `:EqualAI` step between two cut-off forms that step
  is false in the inactive result.
- [ ] **Step 2: Run it.** Expected: FAIL, not defined.
- [ ] **Step 3: Implement** both, and migrate every reader. A caller in a loop calls
  once and keeps the result.
- [ ] **Step 4: Measure.** Time the focused run of
  `test/Enumeration/test_mechanism_enumeration.jl` and of
  `test/Derivation/test_admission.jl` before and after, two runs each. Report the four
  numbers. Stop and report if either is more than 15% slower.
- [ ] **Step 5: Run** the full suite.
- [ ] **Step 6: Commit** "Compute a conformation's steps on demand".

### Task 9: The stored direction of a step, from the graph

Spec decision 1, "Which way a step is stored".

**Files:**
- Modify: `src/Mechanisms/types.jl` (`_orientation`, `_align_tied_groups!`,
  `_canonicalize`, `_side_score`)
- Test: `test/Mechanisms/test_types.jl`; the pins of the steps that turn

**Interfaces:**
- Produces: the constructor stores each step in the direction the spec's table gives,
  reading only the steps, substrates and products. The steps of one kinetic group turn
  together, led by the member nearest the free enzyme among those the direction of the
  reaction decides.
- Reference implementation of the table, for comparison only:
  `claude-scratchpad/lean-mechanism-spike/orientation/check4.jl` (isomerizations, with
  the group rule) and `check6_exchange_local.jl` (exchanges).

- [ ] **Step 1: Write the failing tests**, testset "the stored direction of a step":
  - `F <--> E` in `E + A <--> E(A)`, `E(A) <--> E(P)`, `F + P <--> E(P)`, `F <--> E` is
    stored from `F` to `E` (today it is stored from `E`).
  - `E(S) + P ⇌ E(P) + S` is stored taking up `S`.
  - A dead-end isomerization `E(I) <--> Estar(I)` reached only through `E + I ⇌ E(I)`
    is stored from `E(I)`.
  - Review Focus 2: with the group `(Estar <--> E, Estar(A) <--> E(A))` in the mechanism
    of `test/Mechanisms/test_types.jl` that holds it, writing both members reversed
    gives a mechanism `==` the first, and writing one member reversed gives one `!=` it.
- [ ] **Step 2: Run them.** Expected: the first three FAIL.
- [ ] **Step 3: Implement the table** and the group rule. The alignment by
  conformation labels and its "misaligned" error go.
- [ ] **Step 4: Run the full suite** and sort the failures. A failure is expected only
  where a pin names a constant of a step that turned: the report lists each fixture,
  the step and the old and new name, and their number is at most 35 steps in about 25
  fixtures. Rewrite those pins. Any other failure stops the task.
- [ ] **Step 5: Check the enumerator.** Run `gate/neutral_before.jl`'s population part
  again into a new file and compare with `neutral_before.jls`: every key and verdict
  equal.
- [ ] **Step 6: Commit** "Store each step the way the reaction runs".

### Task 10: `candidate_mechanisms`

Spec decision 6, "The interface" and "Identify"; decision 7, "What identify asks".

**Files:**
- Modify: `src/Enumeration/Enumeration.jl`, `seeds.jl`, `expand.jl`, `filters.jl`
- Modify: `src/identify_rate_equation.jl`, `src/Mechanisms/types.jl` and
  `src/DSL/read.jl` (the decodability check leaves them)
- Modify: `test/Enumeration/test_mechanism_enumeration.jl`,
  `test/test_identify_rate_equation.jl`,
  `test/mechanism_definitions_for_test_enzyme_derivation.jl`, `test/test_boundaries.jl`,
  `test/test_type_stability.jl`, `docs/src/developer.md`

**Interfaces:**
- Produces, public in `Enumeration`:

```julia
struct CandidateMechanism
    mechanism          # the mechanism a caller reads
    # private fields follow in Task 11
end

candidate_mechanisms(reaction; max_param_count = typemax(Int),
                     eq_complexity_filter = typemax(Int),
                     optional_allosteric_regulators = Symbol[],
                     optional_competitive_inhibitors = Symbol[]) -> Vector{CandidateMechanism}
candidate_mechanisms(reaction, parent::CandidateMechanism;
                     max_param_count = typemax(Int),
                     eq_complexity_filter = typemax(Int)) -> Vector{CandidateMechanism}
```

  The `mechanism` field is typed with the mechanism type or types of that commit. Two
  candidates are `==` and hash alike exactly when their mechanisms are, and
  `hash(c) == hash(c.mechanism)`.
- Produces, in the shared test helpers: `_testhelper_candidate(mechanism) ->
  CandidateMechanism`. Until Task 15 it wraps the mechanism.
- Removed: `Enumeration.starting_mechanisms`, `Enumeration.expand_mechanisms`,
  identify's `_required_regulators` (moves into Enumeration), identify's calls to
  `check_valid_rate_equation` and `rate_equation_term_count`, and
  `Mechanisms.check_decodable_names` as a public name (Enumeration owns the check and
  runs it in the first method).

- [ ] **Step 1: Measure today.** For the bi-bi reaction with one competitive inhibitor
  and for the ter-ter reaction of the enumeration tests, time
  `Enumeration.starting_mechanisms` and, for the first, one expansion of the first 50
  starting mechanisms. Record the numbers.
- [ ] **Step 2: Write the failing tests**, in the enumeration test file:
  - "candidate_mechanisms: the starting set": for the uni-uni reaction, the mechanisms
    of `candidate_mechanisms(rxn)` equal as a set those `starting_mechanisms` gives
    today (write them out as the existing testset does).
  - "candidate_mechanisms: the size limits": for the bi-bi reaction,
    `candidate_mechanisms(rxn; max_param_count = 6)` equals the unlimited result
    filtered by `Derivation.fitted_parameter_count(c.mechanism) <= 6`, and the same for
    `eq_complexity_filter = 40` with `Derivation.rate_equation_term_count`.
  - "candidate_mechanisms: a parent's children": for one parent of an existing
    `expand_mechanisms` testset, the children's mechanisms equal that testset's set.
  - "candidate_mechanisms: an undecodable reaction": the reaction the existing
    decodability test uses makes the first method throw, and constructing the reaction
    no longer throws.
- [ ] **Step 3: Run them.** Expected: FAIL, not defined.
- [ ] **Step 4: Implement** the two methods. The guard and the two limits apply to what
  a method returns, never inside a move's own search. The first method throws on an
  enumerator bug; it returns no failure rows.
- [ ] **Step 5: Migrate identify.** It holds candidates, passes its four keywords
  through, fits and reports `candidate.mechanism`, and asks Derivation three things per
  candidate. Its own private record of a candidate is renamed so that neither its type
  nor its fields can be confused with `CandidateMechanism` (the map lists the struct).
  The progress log drops the two skip counts, and `initial_mechanisms.csv` has no
  failure rows from the starting set.
- [ ] **Step 6: Migrate the tests.** A test that called `expand_mechanisms` or
  `starting_mechanisms` calls the method that replaces it with the default keywords and
  compares `Set(c.mechanism for c in children) == Set(expected)`. The identify tests of
  the two skip counts and of starting-set failure rows go with those features; list
  them.
- [ ] **Step 7: Measure again** as in Step 1 and report both sets of numbers. Stop and
  report if the ter-ter starting call is more than twice as slow.
- [ ] **Step 8: Run** the enumeration and identify test files, then the full suite.
- [ ] **Step 9: Commit** "Give Enumeration one function that returns candidates".

### Task 11: The enumerator's own record of each form

Spec decision 6, "Inside Enumeration". No behaviour changes: the enumeration test file
passes without an edit.

**Files:**
- Modify: `src/Enumeration/seeds.jl`, `moves.jl`, `filters.jl`, `expand.jl`,
  `Enumeration.jl`
- Create: `src/Enumeration/record.jl`

**Interfaces:**
- Consumes: `CandidateMechanism` (Task 10).
- Produces, private to `Enumeration`, in `record.jl`:
  - a record type holding, per form and at the form's index in the mechanism: the
    conformation label, the metabolites bound productively, those bound as modifier
    copies (a pure inhibitor counts as one, as today), and the residual (substrates
    added, products subtracted);
  - `CandidateMechanism` gains its second, private field holding that record;
  - one function that builds a candidate from forms' records and steps, the only place
    in `Enumeration` that calls a `Mechanisms` constructor;
  - one function that makes the record of a mechanism written by hand, the only place
    that reads a form's contents from a `Mechanisms` type. `_testhelper_candidate` calls
    the package's path for now; Task 15 gives it its own.
- After this task no other line under `src/Enumeration/` reads `nodes`, `bound`,
  `added`, `subtracted`, `conformation`, `modifier`, `modifiers`, `groups`,
  `step_kinds`, `toward_products` or `step_sides` from a `Mechanisms` value. Step kinds
  and kinetic groups are computed inside `Enumeration` from the steps and the record.

- [ ] **Step 1: Read** enumeration_identify.md sections 4, 5 and 8 and list every read
  to replace; put the list in the report.
- [ ] **Step 2: Write `record.jl`** with the type and the two functions.
- [ ] **Step 3: Move the reads**, one file at a time, running the enumeration test file
  after each: `seeds.jl`, `moves.jl`, `filters.jl`. Whether a step binds a copy is the
  difference between its two forms' records.
- [ ] **Step 4: Check** with `grep` that the fields above appear only in `record.jl`.
- [ ] **Step 5: Run** the enumeration file, the identify file, then the full suite.
  Expected: no test edited.
- [ ] **Step 6: Commit** "Keep the enumerator's knowledge of a form in its own record".

### Task 12: Plain values; the compiled mechanism types go

Spec decision 3, first bullet, and the removals of decision 1 that concern types.

**Files:**
- Modify: `src/Mechanisms/types.jl`, `accessors.jl`, `Mechanisms.jl`;
  `src/DSL/macros.jl`, `read.jl`; `src/Derivation/codegen.jl`,
  `rate_eq_derivation.jl`, `thermo.jl`; `src/Fitting/fit_rate_equation.jl`;
  `src/Plotting/plot_rate_equation_on_data.jl`; `src/identify_rate_equation.jl`;
  `src/EnzymeRates.jl`
- Modify: tests that wrap a macro in `Mechanism(…)` or `AllostericMechanism(…)` (303
  and 242 in the enumeration file), that call `compile_mechanism`, or that name
  `EnzymeMechanism`; `test/test_compile_budget.jl`; `test/test_boundaries.jl`

**Interfaces:**
- Consumes: `compile_rate_equation` (Task 6).
- Produces: `@enzyme_mechanism` and `@allosteric_mechanism` return plain values
  (`Mechanism`, `AllostericMechanism`). `parameters(m)` and `metabolites(m)` are
  ordinary functions. The encoding of a mechanism in a type parameter lives in
  `Derivation`, behind `compile_rate_equation`. `Mechanism` is exported.
- Removed: `EnzymeMechanism`, `AllostericEnzymeMechanism`, `AbstractEnzymeMechanism`,
  `compile_mechanism`, the conversions between the two forms.

- [ ] **Step 1: Write the failing test** in `test/Mechanisms/test_types.jl`,
  "a macro returns a plain value": `m = @enzyme_mechanism …` is a
  `EnzymeRates.Mechanism`, `m == @enzyme_mechanism …` written with its steps in another
  order, and `compile_rate_equation(m) === compile_rate_equation(deepcopy(m))`.
- [ ] **Step 2: Run it.** Expected: FAIL.
- [ ] **Step 3: Move the encoding** into `Derivation` and make the macros return
  values.
- [ ] **Step 4: Migrate** `src/` and the tests; delete the wrappers.
- [ ] **Step 5: The compile budget.** `test/test_compile_budget.jl` measures the
  compiled rate and counts its specializations; restate its bounds for one mechanism
  type and report the old and new numbers.
- [ ] **Step 6: Run** the full suite.
- [ ] **Step 7: Commit** "Return plain mechanisms; compile only the rate equation".

### Task 13: One type and one macro

Spec decision 1 (the four structs) and decision 4 (one macro), with the catalytic
scheme still holding today's forms.

**Files:**
- Modify: `src/Mechanisms/types.jl`, `accessors.jl`, `segments.jl`, `text.jl`,
  `Mechanisms.jl`; `src/DSL/macros.jl`, `read.jl`, `DSL.jl`; every reader of
  `Mechanisms.catalytic`, of a plain mechanism's fields and of the allosteric fields
  in `src/` (mechanisms_dsl.md section 6: 54 calls of `catalytic`)
- Modify: every test that names `@allosteric_mechanism` or `AllostericMechanism`, or
  reads a mechanism's fields; `test/test_boundaries.jl`

**Interfaces:**
- Produces:

```julia
struct Allostery
    allosteric_states::Vector{Symbol} # one per kinetic group; empty for a plain mechanism
    oligomeric_state::Int
    regulatory_sites::Vector{RegulatorySite}
end

struct Mechanism
    catalytic::CatalyticScheme
    allosteric::Allostery
end
```

  `CatalyticScheme` holds what today's plain `Mechanism` holds. A mechanism is plain
  when `m.allosteric.allosteric_states` is empty, and its `oligomeric_state` is 1.
- Produces: `@enzyme_mechanism` reads both kinds. A mechanism is allosteric when its
  steps carry states: every statement has one or none has. `oligomeric_state:` and
  `regulatory_site(n): X::OnlyA` are read only with states. `allosteric_regulators:`
  is not read. `show` prints this form and `parse_mechanism` reads it back.
- Removed: `AllostericMechanism`, `catalytic`, `@allosteric_mechanism`.

- [ ] **Step 1: Write the failing tests** in `test/DSL/test_dsl.jl`, "one macro":
  - an allosteric mechanism written with `@enzyme_mechanism` equals the same one built
    from `Mechanism(CatalyticScheme(…), Allostery([:OnlyA, :EqualAI, :EqualAI], 2, sites))`;
  - `parse_mechanism(sprint(show, m)) == m` for it and for a plain one;
  - Review Focus 4: `oligomeric_state: 2` with no states, `regulatory_site(2): X::OnlyA`
    with no states, and states on two of three statements each throw an
    `ArgumentError` whose message names the keyword or the statement.
- [ ] **Step 2: Run them.** Expected: FAIL.
- [ ] **Step 3: Implement** the structs, the constructor's checks and the macro.
- [ ] **Step 4: Migrate** `src/`, then the tests: the macro's name, the field paths
  (`m.catalytic.steps`, `m.allosteric.allosteric_states`), and the constructors.
- [ ] **Step 5: Run** the full suite.
- [ ] **Step 6: Commit** "Hold plain and allosteric mechanisms in one type, written by
  one macro".

### Task 14: The fixture translator

Nothing in the repository changes. The task builds and proves the tool Task 15 applies.

**Files:**
- Create: `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/translate/translate.jl`
- Create: `…/translate/dry-run/` (translated copies of `test/` and `docs/src/`),
  `…/translate/report.md`

**Interfaces:**
- Consumes: today's `Mechanisms._form_name`; `gate/neutral_key.jl` (Task 1).
- Produces: `translate.jl`, run as
  `julia --project translate.jl <source dir> <output dir>`. It rewrites, inside every
  mechanism macro block and every `parse_mechanism` string, each form in call notation
  to its name (`E(B; residual = A - P)` to `EB_res_plusA_minusP`, `E(S::Mod)` to
  `ESmod`), each `X::Mod` on a free metabolite to the bare name, and leaves comments,
  spacing outside the forms and everything outside the blocks alone.
- Produces: `fixtures_after_task13.jls`, the neutral keys of the fixtures at this
  commit, a sorted vector per file.

- [ ] **Step 1: Regenerate the keys** at this commit (`neutral_key` for the one type of
  Task 13) and compare with `neutral_before.jls`: per file, every key recorded there is
  still present. A key is free of direction and of constant names, so the earlier tasks
  should have moved none.
- [ ] **Step 2: Write the translator.** It writes each isomerization in the direction
  stored at this commit when the written direction differs from it and the step shares
  a kinetic group with another isomerization; it lists those steps in `report.md`.
- [ ] **Step 3: Dry run** over `test/` and `docs/src/` into `dry-run/`.
- [ ] **Step 4: Check the dry run** without the new macro: for every fixture, the names
  on the translated steps equal, as a multiset, the names `neutral_key` wrote for that
  fixture. Report every block the translator could not rewrite (a form built from an
  interpolated variable, a block inside a string) with its file and line; Task 15
  rewrites those by hand.
- [ ] **Step 5: Report.** No commit. The report gives the number of blocks translated,
  the number left for hand translation, and the steps turned in Step 2.

### Task 15: The lean graph

Spec decision 1 (the data model), decision 4 (form names), decision 5 (what leaves
`Mechanisms`), decision 6 (names and the decodability check in `Enumeration`). This is
the task the earlier ones prepare; read all of decisions 1 and 4 to 6 first.

**Files:**
- Modify: `src/Mechanisms/types.jl`, `names.jl`, `text.jl`, `accessors.jl`,
  `segments.jl`, `Mechanisms.jl`; `src/DSL/read.jl`, `macros.jl`
- Modify: `src/Enumeration/record.jl` (its two chokepoints, the naming rule, the
  decodability check); `src/Derivation/` (the encoding of a mechanism in the compiled
  rate's type, and any read of a field that changed shape)
- Modify: every fixture under `test/` (apply the translator), the shared test helpers,
  `test/test_boundaries.jl`, `test/test_type_stability.jl`
- Modify: `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/neutral_key.jl`
  (a method for the lean type)

**Interfaces:**
- Produces: the spec's `Step` and `CatalyticScheme`:

```julia
struct Step
    from::Int                         # index into enzyme_forms
    to::Int
    taken_up::Vector{Int}             # indices into metabolites
    given_off::Vector{Int}
    rapid_equilibrium::Bool
    kinetic_group::Int
end

struct CatalyticScheme
    metabolites::Vector{Symbol}
    enzyme_forms::Vector{Symbol}      # the free enzyme, named E, first
    substrates::Vector{Int}           # indices into metabolites
    products::Vector{Int}
    steps::Vector{Step}
end
```

  `metabolites` holds substrates, then products, then modifiers, each by name. The
  other forms follow `E` by name. Steps are stored in the direction Task 9 gives and
  sorted; kinetic groups are numbered by first occurrence in that order.
- Produces: the constructor checks representation only: indices in range, a form named
  `E`, every form joined to it, distinct names, uniform kinetic groups, one known
  allosteric state per group.
- Produces: the macro reads a name on a step side that is not a declared metabolite as
  a form, used as written; each side holds exactly one form. `X::Mod` and call notation
  are errors that say what to write instead.
- Produces: `Enumeration` owns the naming rule (today's `_form_name`) and the
  decodability check. `_testhelper_candidate(mechanism)` reads the form names under
  that rule, stated in the helper itself, and errors on a name it cannot read.
- Removed: `EnzymeForm`, `Metabolite`, `form_text`, `check_conformation_label`, the
  stored facts (`modifiers`, `groups`, `step_kinds`, `toward_products`, `step_sides`),
  `_orientation`'s reads of contents, `isless` on forms. `Mechanisms` has the 15 public
  names the spec lists.

- [ ] **Step 1: Write the failing tests** in `test/Mechanisms/test_types.jl`,
  "the lean graph":
  - the spec's decision 4 example parses; its `catalytic.enzyme_forms[1] === :E`; the
    step `EA + B <--> EPQ` has `taken_up == [index of B]` and `given_off == []`;
  - two writings of it (steps in another order, the kinetic groups in another order,
    a binding written as a release) are `==` and hash alike;
  - Review Focus 1: the same mechanism with `E` renamed `Efree` throws an
    `ArgumentError` whose message contains "named `E`", from the macro and from
    `Mechanism(CatalyticScheme(…), Allostery(Symbol[], 1, RegulatorySite[]))`;
  - Review Focus 5: `EA + B <--> EPQ` with a later step `EQP ⇌ EQ + P` (a misspelt
    `EPQ`) throws an error that names `EQP` when no step joins `EQP` to `E`;
  - `E(A)` on a step side and `S::Mod` each throw an error naming the form to write.
- [ ] **Step 2: Run them.** Expected: FAIL.
- [ ] **Step 3: `Mechanisms`.** The structs, the constructor, `==`, `hash`, `isless` on
  steps if a sort needs it, `show`, `metabolites(m)`, `step_text`, the two conformation
  functions, `rapid_equilibrium_segments`, the two cycle functions, on the lean fields.
- [ ] **Step 4: `DSL`.** The macro and `parse_mechanism`.
- [ ] **Step 5: `Enumeration`.** The record is built by the enumerator and carried; the
  two chokepoints of Task 11 change, and nothing else under `src/Enumeration/` should
  need to.
- [ ] **Step 6: `Derivation`.** The compiled rate's type parameter holds the whole
  mechanism as nested tuples of Symbols, Ints and Bools; fall back to an integer
  encoding only if Julia rejects that, and say so in the report.
- [ ] **Step 7: Translate the fixtures** with `translate.jl`, translate by hand the
  blocks Task 14 listed, and rewrite the tests of removed names (`EnzymeForm(`,
  `Metabolite(`, `.nodes`, `form_text`, `_form_name`, `step_kinds`, `step_sides`,
  `toward_products`); list the tests that go with a removed feature.
- [ ] **Step 8: Prove the translation.** Add the lean method of `neutral_key`, compute
  the fixtures' keys and compare with `fixtures_after_task13.jls`: per file, every key
  recorded there is present. Compute the populations and compare with
  `neutral_before.jls`: every key and verdict equal.
- [ ] **Step 9: Run the full suite** and sort the failures by the Global Constraints:
  an order pin whose set is unchanged is rewritten and listed; a moved constant name
  from the spec's exceptions is rewritten and listed (a step the table's last row
  decides may turn again here, since the forms are now ordered by name); anything else
  stops the task.
- [ ] **Step 10: Commit** "Make the mechanism a lean graph with opaque form names". A
  series of commits is fine if the suite passes at each.

### Task 16: Documentation

Spec, "Documentation".

**Files:**
- Modify: every page under `docs/src/` that writes a mechanism or evaluates a law
  (counts per page in `fitting_tests_docs.md` section 6), `docs/src/api.md`,
  `docs/src/developer.md`, `docs/src/identify/enumeration_engine.md`,
  `docs/src/deriving/mwc_allostery.md`, `README.md`, `ROADMAP.md`,
  `docs/src/roadmap.md`, `.claude/CLAUDE.md`, the docstrings in `src/`

- [ ] **Step 1: Translate** the mechanisms in `docs/src/`, `README.md` and the
  docstrings with `translate.jl`, and the prose that names a form.
- [ ] **Step 2: Evaluating a law.** Every page that called `rate_equation` shows
  `rate = compile_rate_equation(m)` and `rate(concs, params, Keq)`. `E_total` leaves the
  text; where the rate is first defined, one sentence says the law gives the rate per
  unit enzyme. `Keq` is the last positional argument in every call shown.
- [ ] **Step 3: `api.md`** lists every export with its present signature
  (`checkdocs = :exports`).
- [ ] **Step 4: `developer.md`**: the data model, the interface tables, what
  `Enumeration` keeps in a candidate, the parameter interface. It and
  `mwc_allostery.md` say "island" for a set of forms that a conformation's steps join.
- [ ] **Step 5: `enumeration_engine.md`** describes what the search proposes, by the
  starting set and the moves, and shows `candidate_mechanisms` as the way to list it.
- [ ] **Step 6: `.claude/CLAUDE.md`**: the module map and the bullets under it; the
  Julia-style bullet on compiled types; the "`rate_equation` speed" guard; the rules
  for Derivation tests and for enumeration-engine tests, rewritten in form names
  ("write chemistry the enumerator's way").
- [ ] **Step 7: `ROADMAP.md`** and `docs/src/roadmap.md`: "Shrink the module
  interfaces" updated for the three modules; the stored-facts known issue updated with
  Task 8's and Task 17's numbers; a note that only tests reach the enumerator's code
  for conformation labels other than `E` until iso mechanisms are enumerated.
- [ ] **Step 8: Build.** `flock … julia --project=docs docs/make.jl` and
  `flock … julia --project=docs -e 'using Documenter, EnzymeRates; doctest(EnzymeRates)'`
  both finish without an error or a warning about a missing docstring.
- [ ] **Step 9: Run** the full suite (the boundary test reads the developer page).
- [ ] **Step 10: Commit** "Bring the documentation to the lean mechanism".

### Task 17: The gate, on the lean graph

Spec, "Tests", last bullet.

**Files:**
- Modify: `/home/denis.linux/claude-scratchpad/lean-mechanism-spike/gate/populations.jl`,
  `compare.jl` (the lean graph and the new names)
- Create: `…/gate/report-lean.md`

- [ ] **Step 1: Bring the gate script to the lean graph**: ten populations, 78,554
  mechanisms.
- [ ] **Step 2: Run it** and compare every verdict with `gate/report.md`'s: all equal.
  A difference stops the task.
- [ ] **Step 3: Measure** the stored-facts issue as `ROADMAP.md` describes it, and the
  time of `candidate_mechanisms` for the ter-ter reaction; put the numbers in
  `report-lean.md` and in the ROADMAP entry.
- [ ] **Step 4: Run the full suite** once more and record the pass count.
- [ ] **Step 5: Commit** the ROADMAP numbers, "Record the lean graph's gate and
  timings".

---

## After the last task

The controller runs the whole-branch review, then deletes this plan and the spec in the
branch's last commit (`.claude/CLAUDE.md`, "Workflow"). The maintainer reads, in the
final report: every pin rewritten under the order ruling, every constant renamed, every
test that went with a removed feature, and the timings of Tasks 8, 10 and 17.
