# EnzymeRates: instructions for agents

## Goal: low complexity

Complexity is whatever makes this package hard to understand and change. It comes from
dependencies (code that cannot be understood or changed alone) and from obscurity (facts a
reader needs and cannot see), and it grows a little with every careless change. Working
code is not enough: each change, a bug fix included, leaves the code it touches simpler
than it found it.

Deep modules keep complexity down. `src/` is a few modules, each with a small interface and
its complexity inside, so that using a module, or fixing a bug in another one, needs only
its interface and not how it works. The same holds for a function: a deep function does a
lot behind a simple signature.

Two patterns defeat this, and both arrive one small helper at a time:

- Functions that call each other across the package. Any change then needs knowledge of
  the whole package and breaks code that looked unrelated.
- Many small, shallow functions and modules. Each hides nothing, so there are many
  interfaces to learn and remember.

How each module works inside is on the Developer page of the docs
(<https://DenisTitovLab.github.io/EnzymeRates.jl/>, source `docs/src/developer.md`); read
the part for a module before changing its internals.

## Module map

<!-- module table, generated from src/ by test/update_module_map.jl -->
| Module | What it does | Public names (★ exported) | May call |
|:--|:--|:--|:--|
| EnzymeRates (top level) | Finds the best mechanism for a reaction from measured rates. | ★ `identify_rate_equation`, ★ `IdentifyRateEquationProblem`, ★ `IdentifyRateEquationResults` | all |
| DSL | Turns macro-written reactions and mechanisms into checked objects. | ★ `@enzyme_reaction`, ★ `@enzyme_mechanism`, `parse_mechanism` | Mechanisms, Derivation |
| Plotting | Draws a rate equation over measured rates. | ★ `plot_rate_equation_on_data` | Mechanisms, Derivation, Data |
| Enumeration | Generates candidate mechanisms for a reaction. | `CandidateMechanism`, `candidate_mechanisms` | Mechanisms, Derivation |
| Fitting | Fits a rate law to measured rates. | ★ `fit_rate_equation`, ★ `FittingProblem`, `test_rate_equation` | Mechanisms, Derivation, Data |
| Derivation | Derives a mechanism's rate equation and compiles it. | ★ `compile_rate_equation`, ★ `rate_equation_string`, ★ `parameters`, ★ `dependent_parameter_values`, `RateEquationTooLargeError`, `check_valid_rate_equation`, `fitted_parameter_count`, `parameter_units`, `rate_equation_dedup_key`, `rate_equation_term_count` | Mechanisms |
| Data | Checks and holds a table of measured rates. | `RateData` | none |
| Mechanisms | Stores reactions and mechanisms as plain data. | ★ `metabolites`, ★ `EnzymeReaction`, ★ `Mechanism`, `Allostery`, `CatalyticScheme`, `Reactant`, `Regulator`, `RegulatorySite`, `Step`, `has_nonzero_weight_cycle`, `has_turnover_cycle`, `rapid_equilibrium_segments`, `step_text`, `steps_running_in_allosteric_conformations` | none |
<!-- end of module table -->

- A module lives in `src/<Module>/`, declares its interface in the `public` statement of
  `src/<Module>/<Module>.jl`, and has its tests in `test/<Module>/`.
- A source file `<file>.jl` has its tests in `test_<file>.jl` in the mirrored folder
  (`src/<Module>/` → `test/<Module>/`, `src/` → `test/`, `ext/` → `test/ext/`): a testset
  goes in the file of the function it is about. A check of a package-wide rule is
  `test/test_<rule>.jl`.
- Each source file holds one piece of functionality whole, so that a reader finds what
  belongs together in one place. Add to a file before adding a file. Only a module's
  header includes files.
- Name a source file for what it gives the rest of the package: the exported or public
  function, type or macro it defines (`fit_rate_equation.jl`, `mechanism_type.jl`,
  `enzyme_mechanism_macro.jl`). When no such name says what the file does, name it for
  the calculation it carries out (`rate_equation_derivation.jl`,
  `rate_equation_correctness_guard.jl`). Never `utils.jl`, `helpers.jl` or `types.jl`.
- Keep calls between the files of a module few and one-way: a file does not use a file
  that uses it, directly or through other files. A helper lives in the file of its
  callers; one that two files need lives in the file that owns the idea. The file table
  on the Developer page counts the names each file takes from the others.
- Call another module's function qualified (`Derivation.fitted_parameter_count(m)`), and
  only if it is public there. Inside a module, call its own names bare.
- Only the exports (★) are promised to users. Users and tests can also reach
  `EnzymeRates.test_rate_equation`, `EnzymeRates.Enumeration.candidate_mechanisms` and
  `EnzymeRates.Derivation.RateEquationTooLargeError` (the error of a law over the term
  cap), without a promise. Code in `src/` still follows the bullet above.
- The fields of a `Mechanism` are part of the `Mechanisms` interface, read with dots
  (`m.catalytic.steps[k].from`): its `catalytic` scheme holds `metabolites`,
  `enzyme_forms`, `substrates` and `products` (indices into `metabolites`) and `steps`,
  and its `allosteric` part holds `allosteric_states`, `oligomeric_state` and
  `regulatory_sites`. A mechanism stores nothing derived from them. A form's name is
  opaque: only `E`, the free enzyme, means something to the package.
- Enumeration works on candidates (`CandidateMechanism`): a `Mechanism` and the
  enumerator's private record of each form; a caller reads only `c.mechanism`.
  `compile_rate_equation(m)` compiles a mechanism's rate equation when its value is needed.
- A data table's `group` column marks the measurements made with one amount of enzyme;
  a rate law gives the rate per unit enzyme. `Keq` is always given by the user, never
  estimated from data, and is the last positional argument wherever a function takes it.

## Design rules

- Make functions and modules deep: a simple interface over substantial work. An
  interface's cost is its names, its arguments and their complexity, its errors and its
  preconditions.
- A simple interface matters more than a simple implementation. Handle the hard case
  inside rather than passing it to callers as an argument, an option or an error.
- Each design decision lives in one module. When two modules know the same format or
  rule, merge them or move the knowledge behind one interface.
- Build only what is needed now, behind an interface general enough that it does not name
  its one current caller.
- No pass-through functions, no wrappers with one caller, no arguments that only travel
  through. Do not split a function for its length; split it when the parts can be
  understood separately.
- Extract repeated code when the snippet is substantial and the helper's signature is
  simple. Two similar lines are cheaper than one more interface.
- Keep error cases few. Define an operation so that the odd input is a normal case.
  Validate user input once, at the entry point, with a message in the user's terms; inside,
  trust the types.
- Write obvious code: a reader's first guess about what it does should be right. Prefer
  that to brevity or cleverness.
- A value is an internal constant unless users have a real reason to change it, as they do
  for optimizer settings.
- Remove an unused feature entirely rather than adding a way to disable it. Backward
  compatibility needs the maintainer's approval.

## Changing an interface

Design an interface before implementing it. For a function shared across files, write
the docstring first. If the docstring must be long, or must describe the implementation,
the design is shallow: change the design, not the docstring.

Adding a public name, adding an argument to one, or changing an exported function needs
the maintainer's approval before implementation. Sketch two or three different interfaces,
compare them from the caller's side, and bring the comparison with a recommendation. Prefer
the one where each fact has one owner, each task has one way to do it, and callers need to
know least.

Every change looks for a public name or an argument it can remove. Removing or shrinking a
non-exported shared name needs no approval; do it and say so in the report.

## Playbooks

### Feature

1. Design the interface as above.
2. Write failing tests against that interface and watch them fail.
3. Implement until they pass.
4. Re-read the changed files for dead code, shallow helpers and anything the change made
   unnecessary.

### Refactor

1. Pin the behaviour with tests before moving anything.
2. Delete first: dead code, one-caller wrappers, redundant checks. Then reshape.
3. Migrate every caller and remove the old name in the same change. No parallel old and
   new paths.
4. Re-read the changed files for dead code and further simplification.
5. Keep the result only if a reader now needs to know less. Judge by public names,
   arguments and cross-module calls first, and by line count second.
6. Give a refactor its own commit, apart from a fix, wherever the two can be separated.
   Both belong in one pull request when the refactor simplifies the code around the fix
   or grew out of its investigation.

### Bug fix

Follow `superpowers:systematic-debugging`: find the root cause before any fix, and never
patch a symptom. Reproduce the bug in the simplest failing test, watch it fail, fix one
thing at a time, and leave the touched code simpler than the smallest patch would.

### Optimization

1. Measure before changing anything, and name what limits the time.
2. Time both versions in alternating runs and report the median and the range. Check
   allocations with `@allocated`, type stability with `test/test_type_stability.jl` under
   `julia +1.13`, and compile time with `test/test_compile_budget.jl`.
3. Simplify the hot path before adding machinery to it.
4. Revert a change that does not measurably help, unless it makes the code simpler.

## Julia style

- Dispatch on types rather than branching with `isa` or tags; walking an `Expr` is the
  exception. Add a trait only where a second implementation exists.
- Write type-stable code and keep struct fields concretely typed. Use type parameters for
  that, not to specialize per mechanism: the compiled rate law that `compile_rate_equation`
  returns carries its mechanism in its type (`Sig`) because its call is generated, and
  specializing further costs compile time.
- On Julia 1.13 a closure over a `@nospecialize` argument is untyped; use `Base.Fix2`
  in its place.
- 92-character lines, 4-space indentation. Match the style of the surrounding code, and
  leave whitespace alone unless a formatter changes it.
- Where this file is silent, follow the SciML style guide
  (<https://docs.sciml.ai/SciMLStyle/stable/>). This repo's conventions win where they
  differ: one leading `_` for an internal name, plain `@testset`.

## Names, comments and docstrings

A name says what the thing does, specifically enough that a reader can follow a call
without opening its docstring. Read a name together with its arguments: `name(p, m)` is
clear, and dispatch means it needs no `parameter_` prefix. A longer name is better than a
vague one:

- `has_turnover_cycle(m)`, not `check(m)`
- `dependent_parameter_values(m, ...)`, not `compute_values(m, ...)`
- `rate_equation_dedup_key(m)`, not `key(m)`
- `_kept_forms_keep_a_binding(m)`, not `_check_forms(m)`

Follow Julia conventions: `lowercase` functions, `CamelCase` types, `!` on a function that
mutates its arguments, a leading `_` on an internal name. A name that is hard to choose
often marks a design that is not yet clear.

Names and comments describe the code as it is. They do not mention how it is implemented,
what it used to be or how it changed, so no "new", "old", "legacy", "improved" or
"wrapper".

Write a comment for what the code cannot say: why something exists, a chemical or
numerical reason, the bug that would happen without it, a trap for the next editor. Delete
a comment that restates the code, and update a comment in the same change as its code.

A docstring gives a caller what they need: arguments, result, units and errors. It does
not describe the implementation. A helper private to its file needs none when its name
and arguments make its job obvious.

A spec or plan under `docs/superpowers/` uses the same names as `src/` and the docs
pages. If it needs a term the code lacks, that term is a candidate name for the code, so
choose it as carefully as a function name.

Every code file starts with exactly two comment lines, each beginning `ABOUTME: `, saying
what the file holds.

## Guards

Tests enforce these rules. When one fails, fix the code; change the test only as part of
an approved interface change.

- **Module boundaries.** `test/test_boundaries.jl` parses `src/` and fails on a
  cross-module name that is not public, a module naming one its layer forbids, an
  underscore name that is public or used across modules, and an `export` in a submodule.
  It also fails when the table above or the interface tables, diagram and file table in
  `docs/src/developer.md` disagree with `src/`, when a test file is named for no source
  file, and when a source file has no test file. It enforces the file rules above: it
  fails on a source file named `utils.jl` or the like, or named for a function, type or
  macro it does not define; on an `include` outside a header; on a loop of files in a
  module, each taking names from the next (the uses `ROADMAP.md` lists aside); and on a
  code file without its two `ABOUTME` lines. After changing a `public` statement, an
  export, a module docstring, a cross-module call, a call between the files of a module
  or a file's `ABOUTME` lines, update the developer page in the same commit and run
  `julia test/update_module_map.jl`.
- **Type stability.** `test/test_type_stability.jl` requires that JET find no run-time
  dispatch under any exported or public function, internal callees included. The suite
  runs it on Julia 1.13 only (`julia +1.13`). A new public function needs a row there,
  with a named exemption set for any dispatch it keeps.
- **Rate-law speed.** The compiled rate law, `rate = compile_rate_equation(m)` called as
  `rate(concs, params, Keq)`, is allocation-free and takes under
  `50 ns + 0.5 ns × denominator terms` per call on every gated law: each mechanism in
  `MECHANISM_TEST_SPECS` and the larger laws that test file lists. Its time is the best of
  five batches, each at its own stack depth with its own copy of the inputs.
  `Fitting.loss!` is allocation-free, and its time per data row is at most 1.25 × the
  compiled law's time per call plus 8 `log` calls, met by any of 20 batches, each run at
  its own stack depth. `test/Derivation/test_rate_equation_derivation.jl` enforces both,
  with the times read in one child process started without coverage counters or
  `--check-bounds` (`_testhelper_speed_readings`, `test/testhelpers.jl`) and the
  allocations counted in the test process. `test/Fitting/test_fit_rate_equation.jl`
  requires that `loss!` keep its prediction loop vectorized with the law inlined: in its
  optimized IR on three laws, and in its time per row, at most 0.7 of the law's time per
  call on the 858-term law. A fit evaluates the law millions of times, so stop and discuss
  with the maintainer before a change that would make it allocate or slow down. CI timings
  vary by runner: rerun a job that misses the bound rather than loosening it.
- **Compile time.** `test/test_compile_budget.jl` bounds the first-call compile time of
  a compiled rate law and how many method instances each mechanism's fit compiles. Code
  specialized per mechanism trips it.
- **Enumeration memory.** `test/Enumeration/test_starting_mechanisms.jl` bounds the bytes
  each bi-bi starting candidate holds, 3.5% above today's, and the bytes the bi-bi starting
  set, one expansion and the ter-ter `init_mechanisms` call allocate. A field derived from
  the steps stored on a mechanism or candidate, or more work per candidate, trips it.

Propose a new guard when the maintainer corrects something a deterministic check could
have caught, when existing code breaks a rule in this file, or when you are about to add a
mechanically checkable rule here. Propose it in your report or in `ROADMAP.md`; do not
build it unasked. When a guard lands, shrink its prose rule here to a line naming the
test. A rule that needs judgment, such as depth or naming, stays prose: a test of a
proxy such as length would enforce the wrong thing.

## Tests

- Read all test and system output. It must be pristine: an expected error or warning is
  captured and asserted.
- Tests exercise real behaviour: no test of a mock, and no mocks in end-to-end tests.
  Tell the maintainer when you find an existing test of mocked behaviour.
- A failing test is never deleted or weakened to pass; raise it with the maintainer. Every
  failure you see is yours to fix or report, whoever caused it.
- Do not reduce coverage. New behaviour comes with tests of all of it.
- Every helper under `test/` has the `_testhelper_` prefix; only a one-line closure local
  to a testset or to another helper keeps a bare name.

### Derivation tests

Tests of a compiled rate law (`compile_rate_equation`) live in `test/Derivation/` and
follow two rules:

1. **Write the full mechanism down**, step by step with `@enzyme_mechanism`, in the
   testset or in a `MECHANISM_TEST_SPECS` entry. A fixture from the enumerator or a
   helper cannot be checked against its rate equation by reading. The exception is an
   aggregate check over a whole population, which says so in a comment.
2. **Compare with a separately derived rate equation** at random parameters and
   concentrations: one derived by hand or from a textbook (`analytical_rate_fn`), or from
   the mass-action steady state (`_testhelper_reference_qssa`). The package's own printed
   law or a stored output does not count. Property assertions and golden pins are welcome
   in addition, never instead.

### Enumeration-engine tests

Tests of `candidate_mechanisms`, `init_mechanisms`, `seed_mechanisms` and the expansion
moves live in `test/Enumeration/` (`test_starting_mechanisms.jl`, `test_expansion_moves.jl`,
`test_candidate_mechanism_type.jl`) and follow three rules:

1. **Write the mechanism in the testset**, inline with `@enzyme_mechanism`, even when that
   repeats one used elsewhere, and hand it to the enumerator as `_testhelper_candidate(m)`.
   The exception is an aggregate regression pin over a whole seed set, which says so in a
   comment.
2. **Assert the exact children**: `length(children) == n` and
   `Set(children) == Set(expected)`, the children read as mechanisms (`c.mechanism`,
   `_testhelper_move`, `_testhelper_children`), with every expected child written out.
   Property assertions are welcome in addition, never instead.
3. **Write chemistry the enumerator's way.** Name each form as the enumerator does (`EAB`,
   `ESSmod`, `EB_res_plusA_minusP`): `_testhelper_candidate` reads the form's contents from
   its name. Unmerged chemistry is an isomerization to the product-bound form followed by
   releases (`EA <--> EP_res_plusA_minusP`, then
   `EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P`). A merged complex is the product-bound
   form: the last substrate binds into it in one step (`EA + B <--> EPQ`) and the products
   leave by plain releases. A Theorell–Chance step takes up the substrate and gives off the
   product in one step (`EA + B <--> EQ + P`). Never write a merged complex on the
   substrate side (`EAB <--> EQ + P`): the enumerator never emits it, so it is a different
   mechanism from the one the move produces.

## Commands and workflow

```bash
# Run full test suite (cold — pays precompilation + JIT cost every time). The heap hint
# keeps the peak under the memory available next to a resident VS Code language server.
timeout 1200 julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'
```

```bash
# Focused run of one test file (warm: skips the full suite's precompile and JIT; the
# longest file takes about 75 s). TestEnv activates the test environment in place. A
# module's tests live in test/<Module>/; a file directly in test/ is named
# test/<file>.jl instead.
julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/testhelpers.jl"); include("test/<Module>/<file>.jl")'
```

- Use the focused run while iterating and before each commit, and the full suite once
  before pushing.
- The focused run skips Aqua and JET, and includes only the shared helpers file and the
  file you name. Keep helpers and fixtures that several test files use in
  `test/testhelpers.jl`.
- Run one Julia process at a time: the machine has 7.7 GB and no swap, and the full suite
  needs about 4 GB free.
- Run the full suite in the background under `timeout 1200`, wait for its exit notice and
  read the result before reporting.
- `TestEnv` is installed once in the default environment:
  `julia -e 'using Pkg; Pkg.add("TestEnv")'`.

Workflow:

- Read `ROADMAP.md` before starting. Update it in the commit that fixes, adds or retires
  an entry, and update `docs/src/roadmap.md` too when it lists the entry.
- Fix a broken thing when you find it, unless the fix changes an interface that needs
  approval. Anything else unrelated to the task goes in `ROADMAP.md`.
- Work on a branch for the task and commit often. If the tree has uncommitted changes or
  untracked files when you start, ask how to handle them.
- Never skip or disable a pre-commit hook. Run `git status` before any `git add -A`.
- Specs and plans under `docs/superpowers/` are working documents: delete a change's spec
  and plan in the pull request that implements it.
- A spec opens with the problem and the proposed change in plain language, one paragraph
  each, before any internal name appears.

## Working with the maintainer

- Say so when you do not know something or are stuck, and ask for help. Report results as
  they are, failures included.
- Give your honest technical judgment and push back on an idea you think is wrong, with
  the reason, or say that it is a gut feeling.
- Stop and ask when the decision is the maintainer's: a design choice between valid
  approaches, adding or changing a public name, rewriting a whole implementation or
  module, or anything hard to undo. Otherwise do the task, with its obvious follow-ups:
  choose, say what you chose and what it costs if wrong.
- When asked how to approach something, answer the question; do not start implementing.
- Doing it right matters more than doing it fast. Tedious, systematic work is often the
  correct solution.

### Writing to the maintainer

These rules cover chat messages, specs, plans and pull request descriptions.

- Use the names the code and the docs already have. Write `:NonequalAI`, not "own
  constants"; write `_expand_change_allo_state`, not "the later move" or "the state
  moves". Never give a second name to something that has one.
- Coin a term only for a thing the code has no name for. Define it in one sentence where
  it first appears in each message, and keep such terms to two or three per document.
- Assume the maintainer knows the public names in the module map and nothing below
  them. The first time a message names an internal function (a leading `_`), a private
  rule or a testset, say in one clause what it does and which file holds it. Say it again
  in a later message: do not rely on an earlier one being remembered.
- Explain as to a scientist who has not read the code. Start with the problem in terms of
  enzymes, mechanisms and rate equations, give one small concrete example, and only then
  name the code that is involved.
- Each message stands alone. A reader who skipped the earlier ones can follow it.
- When a question shows that a term or a function was unclear, explain it and stop using
  the unclear form. Do not answer with more new terms.
