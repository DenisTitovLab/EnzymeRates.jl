# Deep Modules and Simplification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Delete dead src code and split EnzymeRates into four deep submodules (`Mechanisms`, `Derivation`, `Enumeration`, `Fitting`) whose small interfaces a parse-only boundary test enforces.

**Architecture:** One PR on branch `deep-modules-simplification`, built by approach B of the spec:
1. Delete code.
2. Cut the submodules with mechanical moves, and add a boundary test whose API lists start as today's crossings.
3. Shrink each list to its target ("ratchet"). Each shrink turns the boundary test red and a commit turns it green.

Refactor commits keep the census, the exact child sets and the derived laws byte-identical.

**Tech Stack:** Julia 1.10–1.13; Test, TestEnv, Aqua, JET; Optimization.jl; Documenter.

**Spec:** `docs/superpowers/specs/2026-10-27-deep-modules-simplification-design.md`. Read it before any task. Its section numbers (§1–§9) are cited below.

## Global Constraints

- **CLAUDE.md applies in full:**
  - TDD;
  - the smallest change;
  - no rewrite of a working implementation without permission;
  - pristine test output;
  - `_testhelper_` prefix on every test helper;
  - `ABOUTME:` two-line header on every new file;
  - 92-character lines and 4-space indentation;
  - comments evergreen, with no history and no "new/old/moved";
  - never remove a comment unless it is false.
- **One Julia process at a time.** The machine has 7.7 GB and no swap. Pass `--heap-size-hint=2500M` to every run.
- **Focused test run:**

  ```
  julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'
  ```

  After Task 11 moves the tests, use `test/<Module>/<file>.jl`.
- **Full suite:** run it once before each task's final commit.

  ```
  julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'
  ```
- **`rate_equation`** stays allocation-free and under 120 ns for every `MECHANISM_TEST_SPECS` mechanism (`_testhelper_test_rate_equation_performance`). If a change would touch it, stop and report.
- **Census gate** (every refactor commit): the census corpus of the tree compares equal to `DM_base` (Task 0) with `compare.jl`: `RESULT: PASS`, zero differences. A commit flagged as a behaviour change records its census diff in the commit message instead.
- **Benchmark gate:** more than 5% slower than the Task 0 baseline on `bench_enum.jl`, `bench_guard.jl` or `bench_pass1.jl` means stop. Write the numbers to `$S/dm/gates.md` and leave the commit unlanded for Denis.
- **Version 0.11.0** (Task 6), the only bump.
- **Attribution.** Every commit message ends with the two attribution lines:

  ```
  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01BFcXywjZosSMs2X8r6rZip
  ```
- **Scratchpad.** `$S` = `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/276b07d3-bafa-4735-bad9-ab925057d392/scratchpad`. Census tools, corpora, benchmarks and logs live in `$S/dm/`, never in the repo.
- **Names fixed by the spec and Denis:**
  - module names `Mechanisms`, `Derivation`, `Enumeration`, `Fitting`;
  - `has_valid_rate_equation`, `check_valid_rate_equation`, `parameter_count`;
  - `rate_equation_key(m) -> (; key, parameters, fit_key, fit_parameters)`;
  - `rate_equation_complexity`, `scale_rate_constants`.

## Review Focus

These are failure modes the spec implies that no task's tests would otherwise exercise. Each line names its owning test.
1. **Relative-mode problem built without `reference`.** The user expects an error that names `reference` and says it is required unless `absolute = true` (Task 6, `test_fitting_scale.jl` "reference is required in relative mode").
2. **A `reference` anchor names a metabolite the mechanism lacks, or has a negative concentration.** This must fail at problem construction, not mid-fit (Task 6, "reference anchors are validated at construction").
3. **The documented internal paths keep working after the split.** `EnzymeRates.Mechanism`, `EnzymeRates.AllostericMechanism` and `EnzymeRates.compile_mechanism` still resolve (Task 11, `test_boundaries.jl` "documented internal paths resolve").
4. **A hand-written mechanism the guard rejects** errors from `@enzyme_mechanism` with the same message text after the macros move to the top level (Task 11: the existing message tests in `test_admission.jl` run unchanged, with qualified internals).
5. **An allosteric mechanism in relative mode with a reference anchor** is scaled with `L` and the regulator constants untouched (Task 6, "reference scaling leaves L and binding constants").

---

## Stage 0: Baselines

### Task 0: Census, benchmark and suite baselines

**Files:**
- Create: `$S/dm/census/` (copy of `~/claude-scratchpad/EnzymeRates-12f42995/phase3/census/{corpus.jl,compare.jl,census_common.jl,neutral_graph.jl,mkenv.jl,rename_phase3.toml}`).
- Create: `$S/dm/bench_enum.jl`, `$S/dm/bench_guard.jl`, `$S/dm/bench_pass1.jl`, `$S/dm/gates.md`, `$S/dm/baseline.md`.

**Interfaces:**
- Produces:
  - corpus `$S/dm/census/DM_base`;
  - env `$S/dm/census/env` (develops the repo by path);
  - three benchmark scripts, each printing one line `TOTAL <seconds>` (the median of 3 runs after one warm-up);
  - `baseline.md` with the full-suite test count and wall time.

- [ ] **Step 1: Copy the census tools and build the env.**

  ```
  julia --heap-size-hint=2500M $S/dm/census/mkenv.jl /home/denis.linux/.julia/dev/EnzymeRates $S/dm/census/env
  ```

  Expected output: `ok 0.10.0 from …/EnzymeRates/src/EnzymeRates.jl`.
- [ ] **Step 2: Build the baseline corpus.**

  ```
  julia --project=$S/dm/census/env --heap-size-hint=2500M $S/dm/census/corpus.jl /home/denis.linux/.julia/dev/EnzymeRates $S/dm/census/DM_base $S/dm/census/rename_phase3.toml
  ```

  Then build it a second time as `DM_base2`, and check that `compare.jl DM_base DM_base2` prints `RESULT: PASS`.
- [ ] **Step 3: Write the benchmarks.**
  - `bench_enum.jl`: `init_mechanisms` and `expand_mechanisms` over the census's six reactions, at the census's level-0 and level-1 parent samples, with a fixed RNG seed.
  - `bench_guard.jl`: `_is_well_posed` over every mechanism of `DM_base`'s `mechs.tsv`, rebuilt with `EnzymeRates._parse_mechanism`.
  - `bench_pass1.jl`: `_pass1_record(m, 337, 20)` over 2,000 mechanisms sampled from bi-bi and ping-pong level 1.

  Each script reaches internals through one `const ER = EnzymeRates`, so later tasks only re-point that one line. Run each and record `TOTAL` in `baseline.md`.
- [ ] **Step 4: Run the full suite** and record the test count and wall time in `baseline.md`. Commit nothing; `$S` is outside the repo.

---

## Stage 1: Deletions (spec §1)

### Task 1: Limit test keeps only its general methods

**Files:**
- Modify: `src/admission.jl`. Delete `_limit_by_aggregation` (1138–1201), `_law_on_face` (1238–1271) and the A2 branch of `_face_stage` (the `if one_law && length(s.aggregated) == 1 … end` block). Delete any helper left without callers; check each with `git grep -w`.
- Test: `test/test_admission.jl:1250–1300`.

**Interfaces:**
- Produces: `_face_limit_check(mech, net, s, face, seed)` with an unchanged signature. Its body becomes `$S/probe_lim/general.jl`'s `_face_limit_check_general`: `limit(pt) = _limit_by_interpolation(...)` and `law = _law_by_reconstruction(...)`.

- [ ] **Step 1: Change the A2 pin.** The row near `test_admission.jl:1279` that pins `_LIMIT_UNTESTED` for an A2-settled face must now expect the general method's verdict: the face runs the limit test and gives `_LIMIT_EQUAL`. Grep the file for other pins of `_LIMIT_UNTESTED` on a face that A2 settles, and change them the same way.
- [ ] **Step 2: Run `test/test_admission.jl`.** Expected: those rows FAIL.
- [ ] **Step 3: Make the deletions** listed under Files.
- [ ] **Step 4: Run `test/test_admission.jl` and `test/test_oracles.jl`.** Expected: PASS, pristine.
- [ ] **Step 5: Census and benchmark gates.** The census must PASS with zero differences: probe 10 found identical verdicts on 318,105 faces. Run `bench_guard.jl` and record its time. The spec accepts open faces going from 0.26 ms to 1.4 ms, so a guard slowdown up to the share of open faces is expected; anything above 5% goes into `gates.md` as a note, not a stop.
- [ ] **Step 6: Full suite, then commit.** Message: "Limit test: interpolation and reconstruction only".

### Task 2: Remove the fast split counter

**Files:**
- Modify:
  - `src/thermo.jl:476–605`: delete `_thermodynamic_constraints`, `_cycle_basis`, `_count_kind`, `_add_step_column!` and `_partition_independent_count`. Keep `_rref_partition` and `_integer_nullspace`, which other code uses.
  - `src/mechanism_enumeration.jl:1157–1170`: delete the `Mechanism` method of `_split_gain_test`. Make the `AllostericMechanism` method (`:1177`) take `m::Union{Mechanism, AllostericMechanism}`, building the child with `typeof(m)(_raw(...))`. Reword its docstring for both types. Update the `_expand_split_kinetic_group` docstring sentence that names `_partition_independent_count`.
- Test: the testsets in `test/` that call `_partition_independent_count`, `_count_kind`, `_cycle_basis` or `_thermodynamic_constraints` (`git grep -n` them).

**Interfaces:**
- Produces: `_split_gain_test(m::Union{Mechanism, AllostericMechanism}, units, flags_of)`, returning `gains(sel)::Bool`.

- [ ] **Step 1: Pin the plain split's children.** Run `test/test_mechanism_enumeration.jl` and record that its split testsets pass. They pin exact children, so they are the regression test.
- [ ] **Step 2: Retire the deleted functions' unit tests.** Delete the testsets that exist only to test the deleted functions. List each by name in the commit message.
- [ ] **Step 3: Make the deletions and the method merge.**
- [ ] **Step 4: Run `test_mechanism_enumeration.jl` and `test_joint_solve.jl`.** Expected: PASS.
- [ ] **Step 5: Census gate (zero differences) and `bench_enum.jl` gate.**
- [ ] **Step 6: Full suite, then commit.** Message: "Count split children on the built child".

### Task 3: Delete class 5

**Files:**
- Modify:
  - `src/admission.jl`:
    - delete the class 5 section (`# ─── Class 5 …`, 1927–~2016): `_law_branches`, `_support_verdict`, `_vanishing_face`, `_assert_support` and `_support_message`;
    - delete `_DENOMINATOR_VANISHES` from the `_Fault` enum and its docstring line (97–98);
    - delete `_guard_message`'s class 5 branch;
    - update the header docstring of `_check_well_posed` (lines 4–17).
  - `src/engine.jl:345–357`: `_with_law` stops calling `_assert_support`. Update the docstrings at 197–198 and 345–346.
  - `src/print.jl:380`, `src/codegen.jl:177`, `src/identify_rate_equation.jl:539`: docstrings that mention class 5.
- Test:
  - `test/test_admission.jl`: 1079–1082, 1269, 3100–3230;
  - `test/mechanism_definitions_for_test_enzyme_derivation.jl`: 1812–1819 (`_testhelper_support_verdict`) and 563.

**Interfaces:**
- Consumes: `rate_equation(m, concs, params)`, `_testhelper_random_independent_params` (shared fixtures), `EnzymeRates._faces`, `EnzymeRates._face_names`.
- Produces: `_testhelper_finite_on_every_face(mech; rng)`: true when `rate_equation` is finite at random positive parameters on every face of `_faces(catalytic(mech))`, with that face's metabolites at 0 and the others drawn uniformly from (0.1, 10).

- [ ] **Step 1: Add the numeric face test.** Write the helper in the shared fixtures file. Replace the testset "every spec's denominator keeps a monomial on every face" with "every spec's rate equation is finite on every face": `@test _testhelper_finite_on_every_face(spec.mechanism; rng = MersenneTwister(5))` for each `MECHANISM_TEST_SPECS` entry. Add "every sampled census mechanism is finite on every face", which runs the same check over `init_mechanisms` plus one `expand_mechanisms` round of the bi-bi and ping-pong reactions (written inline), on a 300-mechanism `MersenneTwister(6)` sample.
- [ ] **Step 2: Turn the class 5 fixtures into guard rejections.** `dead_end`, `unfed` and the allosteric dead end become `@test !EnzymeRates._is_well_posed(x)` with the class 3 fault and face they hit. Read the current fault with `_guard(x)` and pin it.
- [ ] **Step 3: Run both new tests on the current code.** Expected: PASS, since class 3 already rejects these fixtures.
- [ ] **Step 4: Delete the class 5 code** listed under Files, the testsets that call `_support_verdict`, `_vanishing_face`, `_law_branches` or `_assert_support`, and `_testhelper_support_verdict`.
- [ ] **Step 5: Run `test_admission.jl` and `test_rate_eq_derivation.jl`.** Expected: PASS.
- [ ] **Step 6: Census gate (zero differences). Full suite. Commit.** Message: "Delete class 5; pin finiteness on every face numerically".

### Task 4: Drop per-regulator multiplicities `I(1, 2)`

**Files:**
- Modify:
  - `src/dsl.jl`: delete the `I(1, 2)`/`A(4)` regulator syntax. A regulator entry with a call form is a parse error. Delete the docstring bullet at lines 28–31 and the `mults` plumbing of `_parse_reaction_block` (91–125) and `_reaction_expr` (223–235), keeping the catalytic `oligomeric_state:` multiplicities.
  - `src/types.jl`: `RegulatorMults` loses `allowed_multiplicities` and is renamed `DeclaredRegulator(regulator, reg_type)`, together with its accessor. `Base.show(::EnzymeReaction)` stops printing per-regulator lists; its docstring loses the `I(1, 2)` sentence.
  - `ROADMAP.md:140–150`: the entry says the syntax is rejected until the feature lands with it.
- Test: `test/test_dsl.jl`, `test/test_types.jl` (`git grep -n "RegulatorMults\|allowed_multiplicities\|I(1, 2)" test`).

**Interfaces:**
- Produces: `DeclaredRegulator(regulator::Regulator, reg_type::Symbol = :unspecified)`; `regulator(r)`, `reg_type(r)`; `regulators(rxn)::Vector{DeclaredRegulator}`.

- [ ] **Step 1: Write the failing test.** In `test_dsl.jl`, "a regulator multiplicity list is a parse error": `@test_throws "multiplicit" @eval @enzyme_reaction begin substrates: S[C]; products: P[C]; competitive_inhibitors: I(1, 2) end`. Match a substring of the message you write, naming `oligomeric_state:`.
- [ ] **Step 2: Run it.** Expected: FAIL (the syntax parses today).
- [ ] **Step 3: Make the deletions and the rename.** Change every test reference to `RegulatorMults` or `allowed_multiplicities` the same way.
- [ ] **Step 4: Run `test_dsl.jl` and `test_types.jl`.** Expected: PASS.
- [ ] **Step 5: Census gate (no census reaction uses the syntax). Full suite. Commit, including ROADMAP.md.**

### Task 5: Delete dead code

**Files:**
- Modify: `src/types.jl:1416–1419` (`catalytic_mechanism`, `catalytic_multiplicity` on `AllostericEnzymeMechanism`) and `src/identify_rate_equation.jl:43–47` (the unreachable `:L in mnames` check and its comment).
- Test: the lines that call either function (`git grep -n`), which go with them.

- [ ] **Step 1: Delete the code and the test lines that call it.**
- [ ] **Step 2: Focused runs of the affected files, then the full suite.**
- [ ] **Step 3: Commit,** listing the retired test lines in the message.

### Task 6: Relative mode scales to reference anchors; the peak scale goes

**Files:**
- Modify `src/fitting.jl`:
  - `FittingProblem` replaces the `peak_range` field with `reference::Vector{Pair{NamedTuple, Float64}}`.
  - `FittingProblem(mechanism, table; Keq, absolute = false, reference = nothing)`.
  - `fit_rate_equation` returns `(params, loss, retcode)`.
  - Delete `DEFAULT_PEAK_RANGE`, `_peak_range`, `LogConcentrationTurnover`, `_maximize_on_box`, `_coordinate_ascent!`, `_line_maximum!`, `_peak_turnover`, `_finite_or_minus_inf` and the no-anchor method of `rescale_parameter_values`.
  - Rewrite the ABOUTME header and docstrings.
- Modify `src/identify_rate_equation.jl`:
  - `IdentifyRateEquationProblem` replaces `peak_range` with `reference`, passed to every `FittingProblem` it builds, the cross-validation folds included. A fold reads only `fit.params`, and a failed reference still returns parameters (unscaled), so a fold never aborts on it.
  - Delete `peak_at_range_edge` from `_Fit`, `BatchEntry` rows, the CSV and `_failure_row`.
- Modify `Project.toml`: version `0.11.0`.
- Modify docs: `docs/src/fitting/normalized_vs_absolute.md`, `docs/src/fitting/tutorial.md`, `docs/src/identify/model_selection.md`, `docs/src/identify/tutorial.md`, `docs/src/identify/parallel.md`, `docs/src/getting_started.md`.
- Test: `test/test_fitting_scale.jl` (rewrite), `test/test_fitting.jl`, `test/test_identify_rate_equation.jl`, `test/test_compile_budget.jl`, `test/test_admission.jl`. That is every relative-mode `FittingProblem(` or `IdentifyRateEquationProblem(`: 80 of 85 construction sites.

**Interfaces:**
- Produces:
  - `_anchor_log_ratios(m, params, anchors) -> Union{Vector{Float64}, Nothing}`: `log(measured / predicted)` per anchor, or `nothing` when a prediction is zero, non-finite or of the wrong sign. Shared by `fit_rate_equation` and `rescale_parameter_values(m, params, kcats)`, which keeps its error text for the `nothing` case.
  - `_check_anchors(anchors, metabolites)`: errors on a non-`NamedTuple` concentration, an unknown metabolite, a negative or non-finite concentration, or a zero or non-finite turnover. It reuses the `kcats:` messages with the prefix as a parameter (`"reference: "` or `"kcats: "`).
- **Behaviour:** after the restarts, in relative mode, the fit computes `_anchor_log_ratios(fp.mechanism, merge(params, (; Keq)), fp.reference)`.
  - On `nothing`, it returns `(params, loss = Inf, retcode = :ReferenceFailed)` with `params` unscaled.
  - Otherwise it returns `_scale_rate_constants(m, params, exp(mean(ratios)))`.

- [ ] **Step 1: Write the failing tests** in `test_fitting_scale.jl`, replacing the peak testsets (19–294):
  - "reference is required in relative mode": `@test_throws "reference" FittingProblem(m, data; Keq = 2.0)`. With `absolute = true` and no reference, construction succeeds.
  - "reference anchors are validated at construction": an unknown metabolite, a negative concentration and a zero turnover each throw a message starting `reference:`.
  - "relative fit matches the reference turnover": a Michaelis–Menten mechanism written inline, with `reference = [(S = 10.0,) => 3.0]`. After the fit, `rate_equation(m, (S = 10.0, P = 0.0), merge(fit.params, (Keq = 2.0, E_total = 1.0)))` is `≈ 3.0` (rtol `1e-8`).
  - "an unlisted metabolite sits at zero in an anchor": an anchor `(A = 1.0,)` on a bi-bi evaluates with `B = P = Q = 0`, so the turnover is 0 and the fit returns `loss == Inf`, `retcode === :ReferenceFailed`.
  - "reference scaling leaves L and binding constants": an allosteric mechanism (inline, `catalytic_multiplicity: 2`, one `:OnlyA` group) changes only its `k_…` constants, against the unscaled absolute fit run at the same seed.
  - "rescale_parameter_values: measured turnovers" and "…: errors" stay as they are.
- [ ] **Step 2: Run `test_fitting_scale.jl`.** Expected: the new testsets FAIL.
- [ ] **Step 3: Implement and delete** as listed under Files and Interfaces.
- [ ] **Step 4: Give every other relative-mode construction a `reference`.** Use one row of the test's own data: `[(row's concentrations...) => row's Rate]`. A test that compared fitted values against peak-scaled numbers compares against the reference-scaled numbers. Recompute the expected value from the ground-truth parameters with the same anchor; never loosen a tolerance. Do the same in the six docs pages, and in the doctests if any.
- [ ] **Step 5: Run `test_fitting.jl`, `test_fitting_scale.jl` and `test_identify_rate_equation.jl`.** Expected: PASS. Run `test_compile_budget.jl` and keep its budgets unchanged.
- [ ] **Step 6: Full suite. Docs build** (`julia --project=docs docs/make.jl` if `docs/Project.toml` exists). Commit. Message: "Scale relative fits to reference anchors; drop the peak scale (0.11.0)", with the body stating the API break.

---

## Stage 2: Refactor merges (spec §2)

### Task 7: Admission face plumbing (scan item 3)

**Files:**
- Modify: `src/admission.jl`. `_guard` (449–468), `_faces_verdict` (786–793) and `_admission` (1911–1924) each loop over faces and return the first fault. Write one `_first_face_fault(check, faces) -> _Verdict`, where `check(i)` returns a verdict, and use it in all three.
- Inline `_admission_face` (1902–1905) into its two callers if `_first_face_fault` leaves it one line.
- Merge `_tags` (`admission.jl:743–744`) into `_solve_tags` (`thermo.jl:40–41`) only if every reader of `_tags` gives identical results with `fill(:None, ngroups(m))` for a plain mechanism. Otherwise keep both, and say why in the task report.
- Test: none new. This is behaviour-preserving, and the existing guard and admission tests pin it.

- [ ] **Step 1: Refactor.**
- [ ] **Step 2: Run `test_admission.jl`.** Expected: PASS.
- [ ] **Step 3: Census gate. `bench_guard.jl` gate. Full suite. Commit.**

### Task 8: Merge plain/allosteric method pairs (scan item 6)

**Files:**
- Modify:
  - `src/engine.jl:349–357`: one `_with_law(f, m::_AnyMech, s::_Solve = _solve(m))`.
  - `src/dsl.jl:848–855`: one `_checked_mechanism(raw::Union{RawMechanism, RawAllostericMechanism})`.
  - `src/types.jl`: the two raw constructors `Mechanism(raw)` / `AllostericMechanism(raw)`, wherever their bodies coincide apart from the type.

  Merge only where the bodies are identical apart from dispatch, and keep both docstrings' facts.

- [ ] **Step 1: Refactor. Step 2: Focused runs of `test_types.jl`, `test_dsl.jl`, `test_law.jl`. Step 3: Census gate, full suite, commit.**

### Task 9: One exclusion rule (scan item 7)

**Files:**
- Modify:
  - `src/mechanism_enumeration.jl:1924–1950`: `_state_conformations`, `_site_active_states`.
  - `src/rate_eq_derivation.jl:138–158`: `_kreg_params`, `_reg_K`.

  Each currently restates which conformations a tag excludes. Express each through `_excluded_by(tag, inactive)` (`src/types.jl:1678`).

- [ ] **Step 1: Refactor. Step 2: Focused runs of `test_mechanism_enumeration.jl` (merge-sites and regulator testsets) and `test_rate_eq_derivation.jl`. Step 3: Census gate, full suite, commit.**

---

## Stage 3: Cut the submodules (spec §3, §6, §7 step 3)

**Cut mechanics, shared by Tasks 11–14:**
- **Layout.** Each module is a folder `src/<Module>/` with `<Module>.jl`:
  - two ABOUTME lines;
  - `module <Module>`;
  - the `using`/`import` lines its files need: `using ..Mechanisms: Mechanisms` and so on for the siblings it may use, plus the package deps it uses;
  - `const API = (...)`;
  - its `include`s;
  - `end`.

  The top-level `src/EnzymeRates.jl` includes the module files in layer order and keeps the exports.
- **References.** Inside a module, its own names stay bare. Every reference to a sibling's name is qualified: `Mechanisms.catalytic(m)`, `x::Mechanisms.Mechanism`. A method that extends a sibling's function is written `function Mechanisms.name(...)`.
- **The qualifying script.** `$S/dm/qualify.jl <file> <Module> <names.txt>` prefixes each listed identifier with `<Module>.` using `Base.JuliaSyntax` tokens. It skips:
  - a token preceded by `.`;
  - a token that is a field name in `x.f`;
  - a keyword-argument name;
  - a name bound locally in the enclosing function: an argument, an assignment target, or a `for` or comprehension variable.

  Run it, then review the diff by eye. Never commit an unreviewed script rewrite.
- **Test rewrite.** `$S/dm/requalify_tests.jl <names.txt> <Module>` rewrites `EnzymeRates.<name>` in `test/**/*.jl` to `EnzymeRates.<Module>.<name>`, for each name the module defines.
- **API list.** `const API` starts as exactly the names siblings and the top level use from the module today (from `$S/api/deps.jl` rerun on the new tree), sorted.
- **Census adapter.** After each cut, update `$S/dm/census/neutral_graph.jl` and `corpus.jl` (`ER.<name>` → `ER.<Module>.<name>`), so the corpus format is unchanged.

### Task 10: The boundary test

**Files:**
- Create: `test/test_boundaries.jl`, registered in `test/runtests.jl` after `test_aqua_jet.jl`.

**Interfaces:**
- Produces:
  - `_testhelper_module_sources(root) -> Dict{Symbol, Vector{String}}`: module name to `.jl` files under `src/<Module>/`.
  - `_testhelper_api(root, mod) -> Set{Symbol}`: the `const API` tuple parsed from `src/<mod>/<mod>.jl`.
  - `_testhelper_cross_refs(file) -> Vector{Tuple{Symbol, Symbol, Int}}`: `(module, name, line)` for each `Mod.name` where `Mod` is one of the four modules.
  - `_testhelper_field_refs(file, fields) -> Vector{Tuple{Symbol, Int}}`: each `x.f` with `f ∈ fields`.
  - Constants:
    - `_testhelper_LAYERS = Dict(:Mechanisms => Symbol[], :Derivation => [:Mechanisms], :Enumeration => [:Mechanisms, :Derivation], :Fitting => [:Mechanisms, :Derivation])`;
    - `_testhelper_MECHANISM_FIELDS`: every field name of `Mechanism`, `AllostericMechanism`, `Form`, `Edge`, `Met`, `RegulatorySite` and the `Raw*` types;
    - `_testhelper_FIELD_BAN = Dict{Symbol, Set{Symbol}}()`: per module, the fields banned there. It starts empty, and Stage 5 fills it.
    - `_testhelper_FIELD_ALLOWED`: per file, the `(expression text)` entries allowed despite a matching name, such as `data.group`.

- [ ] **Step 1: Write the scanner self-tests.** On inline source strings:
  - `"Derivation.has_valid_rate_equation(m)"` gives `[(:Derivation, :has_valid_rate_equation, 1)]`;
  - `"x.edges[k].from"` gives `[(:edges, 1), (:from, 1)]` with `fields = Set([:edges, :from])`;
  - `"Base.show(io, x)"` gives no cross-ref.
- [ ] **Step 2: Write the four rules as testsets over `src/`:**
  - "every cross-module name is in its API list";
  - "modules respect the layers" (the top level, `src/*.jl`, may use any module);
  - "no mechanism field access outside Mechanisms" (per `_testhelper_FIELD_BAN`);
  - "no cross-module underscore names". This rule is skipped while a module's `API` still holds underscore names; mark it `@test_broken` per such module, and Task 21 turns it green.

  Add "documented internal paths resolve": `isdefined(EnzymeRates, :Mechanism)`, `:AllostericMechanism`, `:compile_mechanism`.
- [ ] **Step 3: Run it.** Expected:
  - the self-tests PASS;
  - the rule testsets PASS vacuously (no `src/<Module>/` yet);
  - "documented internal paths resolve" PASSES.
- [ ] **Step 4: Commit.** Message: "Add the module boundary test".

### Task 11: Cut `Mechanisms`

**Files:**
- Move (`git mv`): `src/types.jl` and `src/names.jl` to `src/Mechanisms/`. Move `src/dsl.jl`'s parser, printer and readback to `src/Mechanisms/dsl.jl`. The three macros, `_checked_mechanism` and `_compile_checked` move to the new `src/macros.jl` at the top level, since they call the guard.
- Move in:
  - `compile_mechanism` (from `src/mechanism_enumeration.jl`);
  - `_has_unbalanced_cycle` (from enumeration);
  - `_kept_members` (from `src/thermo.jl`; a naming rule);
  - `_kreg_params` and `_reg_K` (from `src/rate_eq_derivation.jl`);
  - `_raw_form`, if enumeration's use of it cannot go through `_raw`.
- Create: `src/Mechanisms/Mechanisms.jl`.
- Modify:
  - `src/EnzymeRates.jl`: include `Mechanisms/Mechanisms.jl` first. Add `using .Mechanisms: Mechanisms, Mechanism, AllostericMechanism, compile_mechanism` and re-export the public types (`EnzymeReaction`, `EnzymeMechanism`, `AllostericEnzymeMechanism`). Qualify `precompile` statements.
  - Every other `src/*.jl`: qualify references to `Mechanisms` names via `qualify.jl`.
- Create (scratchpad, used by Tasks 11–14): `$S/dm/qualify.jl` and `$S/dm/requalify_tests.jl` as described under "Cut mechanics".
- Move tests: `test_types.jl`, `test_reaction_graph.jl`, `test_dsl.jl` and `test_accessors.jl` go to `test/Mechanisms/`. Update `runtests.jl`, and the focused-run command in `.claude/CLAUDE.md` (`test/<Module>/<file>.jl`).

- [ ] **Step 1: Mechanical moves and the module file.** Set `API` to today's crossings.
- [ ] **Step 2: Qualify** the other src files and rewrite the tests (script plus review).
- [ ] **Step 3: Load check.** `julia --project -e 'using EnzymeRates'` gives no warnings and no errors.
- [ ] **Step 4: Run `test/test_boundaries.jl`.** Expected: PASS with `Mechanisms` now scanned: its references to siblings are none, and siblings' uses are in `API`.
- [ ] **Step 5: Census gate (adapter updated) and full suite.** Expected: the same test count as `baseline.md`, adjusted by Stage 1's recorded changes.
- [ ] **Step 6: Commit.** Message: "Cut the Mechanisms module".

### Task 12: Cut `Derivation`

**Files:**
- Move: `admission.jl`, `thermo.jl`, `poly.jl`, `engine.jl`, `normalize.jl`, `print.jl`, `key.jl`, `codegen.jl` and `rate_eq_derivation.jl` to `src/Derivation/`. Split `admission.jl` into `guard.jl` (classes 1–3 wiring, networks, faces), `limit_test.jl` (mod-p arithmetic, `_face_limit_check` and its helpers), `admission.jl` (interior and face condition) and `messages.jl` (`_guard_message` and the fix texts). The split is moves only; no line changes besides headers.
- Create: `src/Derivation/Derivation.jl`.
- Move tests: `test_rate_eq_derivation.jl`, `test_poly.jl`, `test_law.jl`, `test_key.jl`, `test_print.jl`, `test_oracles.jl`, `test_admission.jl`, `test_joint_solve.jl`, `test_pivot_priority_regression.jl`, `test_allosteric_collapse.jl`, `test_allosteric_golden.jl` and `allosteric_ground_truth.jl` go to `test/Derivation/`.
- **Generated code:** the bodies `rate_equation`'s generator returns resolve in `Derivation`. Check that every function the generated Expr calls is a `Derivation` binding or a `Base` name. The rate-equation performance tests and `test_compile_budget.jl` must pass with unchanged budgets.

- [ ] **Steps 1–6 as in Task 11.** Commit message: "Cut the Derivation module".

### Task 13: Cut `Fitting`

**Files:**
- Move: `src/fitting.jl` to `src/Fitting/fitting.jl`; create `src/Fitting/Fitting.jl`.
- Move tests: `test_fitting.jl` and `test_fitting_scale.jl` go to `test/Fitting/`.

- [ ] **Steps 1–6 as in Task 11.** Commit message: "Cut the Fitting module".

### Task 14: Cut `Enumeration`

**Files:**
- Move: `src/mechanism_enumeration.jl` to `src/Enumeration/`, split into:
  - `seeds.jl`: `init_mechanisms`, `seed_mechanisms`, topologies and seed variants;
  - `moves.jl`: the `_expand_*` moves and their searches;
  - `filters.jl`: the dwell gauge, the hyperbolic check, `_degenerate`, the emission and atom assertions, and `_respects_reg_type`;
  - `expand.jl`: `expand_mechanisms`.

  Moves only.
- Create: `src/Enumeration/Enumeration.jl`.
- Move tests: `test_mechanism_enumeration.jl` goes to `test/Enumeration/`. `test_identify_rate_equation.jl`, `test_compile_budget.jl` and `test_aqua_jet.jl` stay in `test/`.

- [ ] **Steps 1–6 as in Task 11.** Commit message: "Cut the Enumeration module".

---

## Stage 4: Ratchet the interfaces (spec §4)

Each task shrinks one `API` list to its target. The boundary test goes red, the code changes, and the test goes green again. Write the new target list first: that is the failing test.

### Task 15: `Fitting` interface

**Files:**
- Modify:
  - `src/Fitting/Fitting.jl`: set `API = (:FittingProblem, :cross_validation_loss, :fit_rate_equation, :rate_table, :rescale_parameter_values)`. `FittingProblem` and `rescale_parameter_values` are public names users already have.
  - `src/Fitting/fitting.jl`: rename `_rate_table` to `rate_table(table, metabolite_names, Keq) -> NamedTuple`. Add `cross_validation_loss(mechanism, data, held_out; Keq, optimizer, kwargs...) -> Float64`, which moves in `_cv_fold_loss` and `_fold_problem` from identify.
  - `src/identify_rate_equation.jl`: use only `Fitting.<API>`.
- Test: the identify LOOCV tests stay unchanged. Add `test/Fitting/test_fitting.jl` "cross_validation_loss is the held-out group's loss of the fit on the others": on a two-group inline dataset, it equals `loss!` of `fit_rate_equation` on group 1, evaluated on group 2, with the same `Random.seed!`.

**Interfaces:**
- Produces: `Fitting.cross_validation_loss`, `Fitting.rate_table`. The fold problems take the identify problem's `Keq`, `absolute` and `reference`, as Task 6 set.

- [ ] **Step 1: Write the target `API` and the new test. Step 2: Run the boundary test and the new test (FAIL). Step 3: Implement. Step 4: Run `test_identify_rate_equation.jl`, the boundary test and the Fitting tests (PASS). Step 5: Census gate, full suite, commit.**

### Task 16: `Enumeration` interface

**Files:**
- Modify:
  - `src/Enumeration/Enumeration.jl`: set `API = (:base_mechanisms, :expand_mechanisms)`.
  - `src/Enumeration/seeds.jl`: add `base_mechanisms(rxn, required_allo::Set{Symbol}, required_comp::Set{Symbol}) -> (mechanisms, failures::Vector{Tuple{Any, String}})`. It moves in identify's `_base_tier` (`identify_rate_equation.jl:840–857`), with the seed choice from identify (around line 950): `init_mechanisms` when both sets are empty, else `seed_mechanisms`.
  - `src/identify_rate_equation.jl`: call `Enumeration.base_mechanisms` and wrap each failure into `FitFailure`.
  - `src/EnzymeRates.jl` and `.claude/CLAUDE.md`: the API Design section's `EnzymeRates.init_mechanisms` path becomes `EnzymeRates.Enumeration.init_mechanisms`.
- Test: move identify's base-tier testsets into `test/Enumeration/` against `Enumeration.base_mechanisms`, with the same expected mechanisms.

- [ ] **Steps as in Task 15.** Commit message: "Enumeration exposes base_mechanisms and expand_mechanisms".

### Task 17: `Derivation` interface for identify, fitting and the macros

**Files:**
- Modify `src/Derivation/`:
  - add `has_valid_rate_equation(m) = _is_well_posed(m)` and `check_valid_rate_equation(m) = _check_well_posed(m)` (renames, with callers moved);
  - `rate_equation_key(m) -> (; key::String, parameters::Vector{Symbol}, fit_key::String, fit_parameters::Vector{Symbol})`. It derives the law once and computes `_canonical(L, id)`. For an `AllostericMechanism` with `_is_phantom(L)` whose active law holds every fitted constant but `L`, `fit_*` come from `_canonical(L.active, id)`; otherwise they equal `key` and `parameters`. `id` is `_mechanism_id(sprint(show, m))`, moved into `Derivation`;
  - `rate_equation_complexity(m) = _eq_complexity(m)`;
  - `scale_rate_constants(m, params, factor)`, replacing `_scale_rate_constants` and `_ss_rate_constant_names` for `Fitting`;
  - `fitted_params(m::Union{Mechanism, AllostericMechanism})` reads `_solve(m).fitted` instead of `compile_mechanism`.
- Modify `src/identify_rate_equation.jl`: `_pass1_record` uses only `rate_equation_complexity`, `check_valid_rate_equation`, `fitted_params`, `rate_equation_string(m; coefficients = true)` and `rate_equation_key`. `_Candidate` keeps its fields, filled from these. Keep the leaf-cap log line through a `capped` field only if it costs no extra export. Otherwise drop the log line and say so in the commit message.
- Modify: `src/macros.jl` uses `Derivation.check_valid_rate_equation`.
- Test: in `test/Derivation/test_key.jl`, add "rate_equation_key gives a phantom its own key and the active law's fit key", using an inline one-subunit L-phantom mechanism. Find one from the existing phantom tests in `test_identify_rate_equation.jl`.

- [ ] **Step 1: Target `API` and new tests. Step 2: Red. Step 3: Implement.**
- [ ] **Step 4: `bench_pass1.jl` gate.** Above 5% slower, commit the work to a side branch `gate/pass1`, write the timings to `$S/dm/gates.md`, and continue from the previous commit. Denis decides between the slowdown and a private derivation cache.
- [ ] **Step 5: Identify tests, census gate, full suite, commit.**

### Tasks 18–20: Enumeration asks only `has_valid_rate_equation` and `parameter_count`

These are three tasks, one per search, in this order: 18 the flip search (`_minimal_flips`), 19 the tagging search (`_tagging_bases`, `_admits`, `_stops_inactive_turnover`), 20 the split search (`_expand_split_kinetic_group`'s `gains`).

**Files:**
- Modify: `src/Enumeration/moves.jl`. Each search builds its candidate (`_build(_raw(m, re))`, or the tagged `AllostericMechanism`) and asks `Derivation.has_valid_rate_equation`. The gain test asks `Derivation.parameter_count` (`= _independent_param_count`, added to `Derivation`'s `API` in Task 18).
- Delete each guard-stage call's crossing from `Derivation`'s `API` once nothing uses it: `_faces_pass`, `_bottomless_segment`, `_admission`, `_NO_FAULT`, `_is_well_posed` (via the exits), `_independent_param_count` and `_has_turnover_cycle`. The last one moves to `Mechanisms` in Task 23.

- [ ] **Step 1: Shrink `Derivation.API`** by the names this search uses. Run the boundary test: red.
- [ ] **Step 2: Rewrite the search.**
- [ ] **Step 3: Census gate and `bench_enum.jl` gate.**
  - **Zero differences and under 5% slower:** run the full suite and commit as a refactor.
  - **Any census difference or a larger slowdown:** write the diff counts per reaction and level, three example mechanisms gained or lost, and the timings to `$S/dm/gates.md`. Commit the work to a side branch `gate/<search>` (not the PR branch) and skip to the next task. Denis decides on waking.

### Task 21: Underscore names leave the API lists

**Files:**
- Modify: every `API` that still holds an underscore name. Each such name is either replaced by an existing API function or moved into the calling module. Siblings of `Mechanisms` keep using its underscore names only until Stage 5; record those as the remaining exceptions in the boundary test's `@test_broken`.

- [ ] **Steps:** target lists, red, fix, green, census gate, full suite, commit.

---

## Stage 5: `Mechanisms` accessors (spec §5)

**Rule:** a derived view goes in `Mechanisms` only when at least two modules use it.

### Task 22: The accessors

**Files:**
- Create: `src/Mechanisms/accessors.jl` and `test/Mechanisms/test_graph_accessors.jl`.

**Interfaces:**
- Produces, all allocation-free unless the result is a fresh vector by nature:
  - **Species:** `species(m)::UnitRange{Int}`, `bound(m, x)::Metabolites`, `conformation(m, x)::Symbol`.
  - **Steps:** `steps(m)::UnitRange{Int}`, `ends(m, k)::Tuple{Int, Int}` in stored direction, `taken_up(m, k)::Metabolites`, `given_off(m, k)::Metabolites`, `is_rapid_equilibrium(m, k)::Bool`, and `step_kind(m, k)::Symbol`, one of `:binding`, `:isomerization` or `:chemistry` (from `_is_binding`, `_is_iso`, `_is_chemistry`).
  - **Groups:** `group(m, k)::Int`, `members(m, g)::Vector{Int}`, `ngroups(m)`.
  - **Metabolites:** `Metabolites`, an immutable wrapper of `Counts` plus the table it indexes, with `copies(s, name::Symbol)`, `names(s)::Vector{Symbol}`, `isempty`, `∪` and `∩`.
  - **Allosteric:** `catalytic(am)`, `tag(am, g)::Symbol`, `multiplicity(am)::Int`, `regulatory_sites(am)`.

  Write each signature's test first. Test every accessor on two inline mechanisms, plain random-order bi-bi and a two-conformation allosteric uni-uni, with exact expected values. Add `@allocated` checks equal to 0 for the scalar accessors.

- [ ] **Steps: tests red, implement, green, full suite, commit.**

### Task 23: Derived views (scan items 2, 5, 8, 10)

This is four commits, each with its gate.
1. **`conformation_steps(am, c::Symbol)::BitVector`** (`:active`/`:inactive`). It replaces `_active_mask`, `_inactive_mask`, `_active_law`, `_inactive_law` and `_state_mask` everywhere.
2. **`has_turnover_cycle(m, on::BitVector)::Bool`.** One implementation for `_has_turnover_cycle` (guard) and `_re_turnover_cycle` (enumeration, which passes the fast mask). **Gate (item 8):** first add tests that pin both current functions on every `MECHANISM_TEST_SPECS` mechanism and every subset mask of size ≤ 2 of five inline mechanisms; the merged function must give the same Booleans. `_has_unbalanced_cycle` becomes its private helper.
3. **`re_segments(m, on::BitVector, re::BitVector)`**, returning `(; segment_of, segments, bottom)`, and **`segment_graph(m, on)`**. One walker replaces `_segments`, `_fast_segments`'s graph part and `_bottomless_segment`. **Gate (item 10):** `bench_guard.jl` and `bench_enum.jl`; drop this commit if either is slower.
4. **`present_species(m, on)::BitVector`** replaces `_present_forms`.

- [ ] **For each: pinning tests, implement, switch callers, census gate, benchmark gate, full suite, commit.**

### Task 24: `Derivation` reads mechanisms only through accessors

**Files:**
- Modify: every file in `src/Derivation/` (about 344 field sites).
- Modify: `test/test_boundaries.jl`. Set `_testhelper_FIELD_BAN[:Derivation] = _testhelper_MECHANISM_FIELDS`. That is the failing test.

- [ ] **Steps:** ban (red), rewrite (`ends(m, k)` for `e.from`/`e.to`, `taken_up` for `e.consumed`, …), green, census gate, `bench_guard.jl` gate, rate-equation performance tests, full suite, commit.

### Task 25: `Enumeration` reads and writes through accessors

**Files:**
- Modify: every file in `src/Enumeration/` (about 341 field sites). Reads go through the accessors. Writes go through the `Raw*` types plus `with_flags(m, re)`, `regroup(m, groups)` and `with_tags(am, tags)`, which are added to `Mechanisms` here and replace `_raw`, `_build` and `_with`.
- Modify: `test/test_boundaries.jl` with `_testhelper_FIELD_BAN[:Enumeration]`.

- [ ] **Steps as Task 24,** with the `bench_enum.jl` gate.

### Task 26: Close the `Mechanisms` interface

**Files:**
- Modify:
  - `src/Mechanisms/Mechanisms.jl`: `API` becomes the accessors of Task 22, the views of Task 23, the write functions of Task 25, and the types, raw types, reaction accessors, `compile_mechanism` and the DSL entry points the macros call. Nothing else.
  - The top level, `Fitting` and identify: they read mechanisms through accessors.
  - `test/test_boundaries.jl`: ban fields in every module and in `src/*.jl`, and turn every remaining `@test_broken` into `@test`.

- [ ] **Steps:** target list (red), fix, green, census gate, full suite, commit.

---

## Stage 6: Docs and wrap-up (spec §9)

### Task 27: Docs, CLAUDE.md, ROADMAP

**Files:**
- Modify:
  - `docs/src/developer.md`: the module map, each module's `API` and the boundary test. Remove the peak-scale and class 5 text.
  - `docs/src/api.md`: if it lists `peak_range` or internals by old path.
  - `.claude/CLAUDE.md`:
    - the API Design paths;
    - the focused-test command with `test/<Module>/`;
    - one paragraph on the boundary test and `API` lists;
    - the chokepoint test's location, if it moved.
  - `ROADMAP.md`: retire the entries this PR closes (grep for the peak scale, class 5, the limit-test special cases and the split counter).

- [ ] **Step 1: Edit.**
- [ ] **Step 2: Docs build and doctests.**
- [ ] **Step 3: Full suite on Julia 1.12 and on 1.13** (`julia +1.13 --project …`). Both pristine.
- [ ] **Step 4: Final census against `DM_base`.** Expected: PASS, or exactly the differences of commits flagged as behaviour changes.
- [ ] **Step 5: Commit.**
