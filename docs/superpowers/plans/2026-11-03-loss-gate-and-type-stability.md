# Loss Gate, Type Stability and Shorter CSVs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fail a test on any `loss!` slowdown or allocation, require zero JET optimization reports on every exported and public function, tighten the `rate_equation` bound, and shorten the `rate_equation` text in the search CSVs.

**Architecture:** The `rate_equation` speed gate gains the `loss!` checks, so a law added to the gate gets both. A new `test/test_type_stability.jl` holds one table of functions and one list of exemptions, and grows module by module as each module's sites are fixed. The short law is an option on `rate_equation_string`, which `_pass1_record` calls.

**Tech Stack:** Julia 1.13.3, JET 0.13.1 (`report_opt`), Test, CSV, DataFrames.

**Spec:** `docs/superpowers/specs/2026-11-03-loss-gate-and-type-stability-design.md`. Read it first.

## Global Constraints

- Work in the worktree `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/3794dcd9-f12d-4146-9264-28ded321b71f/scratchpad/guard` (call it `$G`), branch `docs/fit-speed-record`, based on `main` at `9c194361`. Never touch `/home/denis.linux/.julia/dev/EnzymeRates`, which another session owns.
- One Julia process at a time: 7.7 GB, no swap. Pass `--heap-size-hint=2500M`. Never run tests in parallel with another agent.
- Focused run, from `$G`: `julia --project --heap-size-hint=2500M -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'`. A helper used by more than one test file goes in the fixtures file.
- Follow `.claude/CLAUDE.md`: TDD, `_testhelper_` prefix on every test helper, 92 columns, no history in comments or names, Derivation and Enumeration test rules (mechanisms written out inline).
- No change to a derived law, parameter name, child set, key or the text `rate_equation_string` prints today. The hashes in `test/test_compile_budget.jl` must not change, nor its counts rise.
- No new public name or argument other than `coefficient_values`. Do not rewrite a whole module: stop and report instead.
- Exemptions name a function, and only these reasons qualify: an error message, the solver barrier, `Expr` walking, a `DataFrame`, CSV or `pmap` boundary, lifting a value to a type.
- End every commit message with:
  ```
  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_015DwkaAMe9ZdejBj1fMJEtj
  ```

## Review Focus

- A law with `N`, `N_A`, `F_A` or a free-E inactive conformation prints a short form holding exactly the full form's kept lines, in order (Task 2's property test over all gated laws).
- A short law with a superscript name (`dA_ThreePG²`) survives a CSV write and read unchanged (Task 2's identify test compares text read back).
- On Windows and macOS the harness still finds package frames, so it identifies them by module and never by file path (Task 3's mutation test).
- An exemption that no longer matches anything fails, so the list cannot go stale and hide a later bug (Task 3's stale-exemption test).
- On Julia 1.11 and 1.12, including the Lowest compat job, the JET table skips cleanly and the `loss!` gate and `coefficient_values` tests still run (Task 3's version guard; Task 8 reads the CI job).

---

### Task 1: One `rate_equation` bound, and the `loss!` gate

**Files:**
- Modify: `test/mechanism_definitions_for_test_enzyme_derivation.jl`, `test/Derivation/test_rate_eq_derivation.jl`, `test/Fitting/test_fitting.jl`, `.claude/CLAUDE.md`

**Interfaces:**
- Produces, in the fixtures file: `_testhelper_loss_table(mets; rows = 64)`, `_testhelper_loss_problem(m; absolute = false)` (`Keq = 2.0`), `_testhelper_loss_allocations(fp) -> Int`.
- Produces, in `test_rate_eq_derivation.jl`: `_testhelper_rate_equation_time_limit(n_terms) -> seconds`, `_testhelper_loss_time_per_row(m) -> (; per_row, allowed)`.

- [ ] **Step 1: Write the failing gate.** Add to `_testhelper_test_performance` (line 567), replacing `@test t < 120e-9`, and to the large-law loop (line 1376), replacing its bound:

```julia
n = _testhelper_denominator_terms(m)
@test t < _testhelper_rate_equation_time_limit(n)
@test _testhelper_loss_allocations(_testhelper_loss_problem(m)) == 0
@test _testhelper_loss_allocations(_testhelper_loss_problem(m; absolute = true)) == 0
(; per_row, allowed) = _testhelper_loss_time_per_row(m)
@info "$(_testhelper_law_label): rate_equation $(round(1e9t; digits = 1)) ns, $n terms, " *
      "loss! $(round(1e9per_row; digits = 1)) ns per row, allowed $(round(1e9allowed; digits = 1))"
@test per_row <= allowed
```

Log for every law, with the spec's name or the large law's index as the label; this replaces the large loop's `@info`. Keep the `@inferred metabolites(fp.mechanism)` and `@inferred parameters(fp.mechanism, FittedParameters())` assertions from the testset this task removes.

- [ ] **Step 2: Implement the helpers.** `_testhelper_rate_equation_time_limit(n) = 25e-9 + 0.5e-9 * n`, replacing `_testhelper_seconds_per_denominator_term` and the comments on both old bounds. The new comment states why: 25 ns + 0.5 ns per term leaves about 5–6× headroom locally on most laws and about 2× over the worst of 168 CI timings. `_testhelper_loss_time_per_row(m)` times three things in one process on the same 64-row table, each the minimum over 100 repeats per row: `Fitting.loss!`; a `@noinline` loop calling `rate_equation` on each row; and a `@noinline` loop taking `log(abs(x))` of the `Rate` column. `allowed = t_rate + 6 * t_log`. Move `_testhelper_loss_table` and the allocation half of `_testhelper_loss_cost` from `test/Fitting/test_fitting.jl` to the fixtures file. In `test_fitting.jl`, delete the testsets "loss! on every spec mechanism…" and "loss! on the largest laws…", now covered by the gate, and adapt the column-type and fold testsets to the moved helpers.

- [ ] **Step 3: Run the Derivation and Fitting files.** Run the focused command on `Derivation/test_rate_eq_derivation.jl`, then on `Fitting/test_fitting.jl`. Expected: all pass, and the `@info` lines show every law's `per_row` below `allowed`. Report the largest `per_row / allowed`. If any law exceeds its bound, stop and report: do not change the constants.

- [ ] **Step 4: Prove the gate fails on the historical bug.** Append to `src/Mechanisms/types.jl` after the method whose body is `Tuple(metabolites(_concrete(M())))`:

```julia
function metabolites(m::AllostericEnzymeMechanism)
    extra = Symbol[]
    (Tuple(metabolites(_concrete(m)))..., extra...)
end
```

Run `Derivation/test_rate_eq_derivation.jl`. Expected: failures on allosteric laws, in `loss!` allocations and time. Record the counts, then `git checkout src/Mechanisms/types.jl`.

- [ ] **Step 5: Update `.claude/CLAUDE.md`'s section "`rate_equation` runtime perf is non-negotiable".** It states `t < 25 ns + 0.5 ns × denominator terms` for every gated law, `allocs == 0`, the `loss!` rules (`allocs == 0` in both modes on every gated law; time per row at most `rate_equation` per call plus 6 `log` calls, timed in one process) and that CI ran up to 3.5× this machine across 28 recent jobs. It still says to rerun a job that misses the bound rather than loosen it.

- [ ] **Step 6: Commit** the four files as "Bound rate_equation by law size and gate loss! on every law".

### Task 2: The short law in the CSVs

**Files:**
- Modify: `src/Derivation/print.jl`, `src/identify_rate_equation.jl`, `test/Derivation/test_print.jl`, `test/test_identify_rate_equation.jl`, `docs/src/getting_started.md`, `docs/src/developer.md`, `ROADMAP.md`

**Interfaces:**
- Produces: `rate_equation_string(m; coefficients::Bool = false, coefficient_values::Bool = true) -> String`, on the concrete and the singleton method.

- [ ] **Step 1: Write the failing tests in `test/Derivation/test_print.jl`,** following the Derivation test rules:

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
short(m) = rate_equation_string(m; coefficients = true, coefficient_values = false)
@test short(m) == "(; S, P) = concs\n" *
    "v = E_total * n_0 * (S - P / Keq) / (1 + d_S * S + d_P * P)"
@test_throws ArgumentError rate_equation_string(m; coefficient_values = false)
```

Add a written-out allosteric law with its exact short text (the `m_allo` at `test_rate_eq_derivation.jl:2573` shows the full form to start from). Add the aggregate check, with a comment saying it is one, over `MECHANISM_TEST_SPECS` and both large-law lists, on concrete and singleton forms: the short form's lines are the full coefficient form's lines, in order, minus the `params` line, the comment, and each line of a value (`n…`, `d…`, `f…`, `N…` left-hand sides); every remaining `Q_`, `D_`, `F_`, `W_` line is kept; the last line starts with `v = `.

- [ ] **Step 2: Run `Derivation/test_print.jl`.** Expected: FAIL, unknown keyword `coefficient_values`.

- [ ] **Step 3: Implement.** Add `coefficient_values::Bool = true` to both `rate_equation_string` methods and both `_law_text` methods. `coefficient_values = false` without `coefficients` throws `ArgumentError("coefficient_values = false needs coefficients = true")`. The short form is the `concs` line, then `Q_c`, `D_c`, `F_c`, `W_c` in their current order, then `v`; to take only the `concs` line, split `_header_lines` so nothing parses text. Document the option in the docstring, with a doctest, and say this text names coefficients it does not define and does not run.

- [ ] **Step 4: Run `Derivation/test_print.jl`.** Expected: PASS.

- [ ] **Step 5: Write the failing identify test.** In `test/test_identify_rate_equation.jl`, beside line 597, with `save_dir`'s files read back through `CSV.read`:

```julia
short(m) = rate_equation_string(m; coefficients = true, coefficient_values = false)
searchf = CSV.read(joinpath(save_dir, "equation_search_iteration_1.csv"), DataFrame)
row = first(filter(r -> !ismissing(r.rate_equation), searchf))
@test row.rate_equation == short(EnzymeRates.parse_mechanism(row.mechanism))
@test cvf.rate_equation[1] == short(EnzymeRates.parse_mechanism(cvf.mechanism[1]))
@test !any(startswith("# The coefficients"), split(cvf.rate_equation[1], '\n'))
```

Use the qualified name `parse_mechanism` is really reached by (`EnzymeRates.DSL.parse_mechanism` if it is not visible at the top). The test at line 597, that `best_equation.csv` holds the full form, stays unchanged. Run the file. Expected: FAIL on the new assertions.

- [ ] **Step 6: Wire it.** `_pass1_record` prints with `coefficient_values = false`. In `_cv_model_selection`, `best_equation.csv` gets a one-row copy of `cv_df` whose `rate_equation` is `Derivation.rate_equation_string(candidates[best_row_idx].mech; coefficients = true)`. `results.cv_results` and `loocv_results.csv` keep the short form. Update the `_pass1_record` docstring, and `save_dir` in `identify_rate_equation`'s docstring. Update `docs/src/getting_started.md` and `docs/src/developer.md` where they describe these forms and files.

- [ ] **Step 7: Run `test_identify_rate_equation.jl` and `test_compile_budget.jl`.** Expected: PASS, with `PASS1 method instances per fresh mechanism` still 0 and no pinned hash changed.

- [ ] **Step 8: Delete `ROADMAP.md`'s entry "Drop `rate_equation` from the search CSVs", and commit** as "Print the law without coefficient values in the search and LOOCV CSVs".

### Task 3: The type-stability harness and the `Mechanisms` rows

**Files:**
- Create: `test/test_type_stability.jl`
- Modify: `test/runtests.jl` (include it just before `test_aqua_jet.jl`)

**Interfaces:**
- Produces, in that file: `_TESTHELPER_DISPATCH_EXEMPTIONS::Vector{Pair{Symbol, String}}` (function name => reason), `_TESTHELPER_CHECKED_CALLS::Vector` of `(name::String, f, argtypes::Tuple)`, `_testhelper_keyword_call(f, kwtype, argtypes...)`, `_testhelper_dispatch_sites(f, argtypes) -> Vector{NamedTuple}` with `(; function_name::Symbol, file, line)`.

- [ ] **Step 1: Write the harness and the mutation test.** Write the fixtures inline, one plain and one allosteric mechanism with the `@enzyme_mechanism` and `@allosteric_mechanism` macros; use their types plus `EnzymeRates.Mechanism` and `EnzymeRates.AllostericMechanism` for the concrete forms. `$S/bench/recon4.jl` (`$S` is `$G/..`) has a working sweep to read, including the argument types of every function. Wrap everything in a guard:

```julia
if v"1.13" <= VERSION < v"1.14"
    using JET
    # testsets
else
    @test_skip "zero JET reports is required on Julia 1.13 only"
end
```

`_testhelper_dispatch_sites` runs `JET.report_opt(f, argtypes)` and, for each of `JET.get_reports(result)`, takes the innermost frame of `report.vst` whose `linfo.def.module` is `EnzymeRates` or a submodule (compare modules, never file paths). The function name drops closure and keyword-body decoration: `#_cv_model_selection#97` and `#_raw_step!##3` give `_cv_model_selection` and `_raw_step!`. It returns the unique sites whose name is not exempt. It must not filter by the last frame's module, since `_rref!`'s reports end in `Base`.

Testsets:
- `"$name reports no run-time dispatch"` for each row: `@test isempty(_testhelper_dispatch_sites(f, argtypes))`, printing the sites on failure.
- A mutation test: define a local function with a boxed captured variable or a call on a `Ref{Any}`, inside a throwaway `module _TestStabilityFixture`, and show that `report_opt` flags it, so the harness cannot silently report zero.
- A stale-exemption check: every name in the exemption list is reported by at least one row with the exemptions turned off (`_testhelper_dispatch_sites(f, argtypes; apply_exemptions = false)`).

- [ ] **Step 2: Add the rows that are already clean,** concrete and singleton where methods exist: `metabolites`, `compile_mechanism`, `catalytic`, `check_conformation_label`, `check_decodable_names`, `excludes`, `form_text`, `step_text`, `has_turnover_cycle`, `has_unbalanced_cycle`, `present_species`, `rapid_equilibrium_segments`, `segment_walk`, and the `Mechanism`, `AllostericMechanism`, `EnzymeMechanism` and `AllostericEnzymeMechanism` constructors.

- [ ] **Step 3: Run the file on its own, and verify.** Expected: all rows pass. Put the uninferable `metabolites` of Task 1 back, add a temporary row `loss!` on an allosteric problem, and confirm the row reports sites; revert both. Report the time and peak memory of the file (`/usr/bin/time -v`).

- [ ] **Step 4: Commit** as "Add the type-stability test with the Mechanisms functions".

### Task 4: `Derivation` sites

**Files:**
- Modify: `test/test_type_stability.jl`; as needed `src/Derivation/{engine,poly,normalize,thermo,print,rate_eq_derivation}.jl`

**Interfaces:** Consumes the table and exemption list from Task 3.

Reported sites are in `$S/bench/recon4.log`. Confirm each from the code; the spec's triage table gives the expected action.

| Site | Expected action |
|:--|:--|
| `engine.jl:35, 54` (`_spanning_trees`), `poly.jl:199, 201` (`_lowest_shift`) | Fix: a recursive or reassigning closure becomes a local-free helper or a typed `Ref` |
| `normalize.jl:64–83` (`_normalizer`) | Fix: captured variable; and `reduce(hcat, …)` loses the matrix type |
| `thermo.jl:146–160` (`_rref!`), `thermo.jl:274, 277` (`_solve`) | Fix the callers first, then re-run |
| `engine.jl:338`, `print.jl:181`, `rate_eq_derivation.jl:112, 113` | Fix if small; else exempt with a reason |
| `messages.jl:22, 142` (`_guard_message`, `_face_message`) | Exempt: error message |
| `codegen.jl:214` (`rate_equation` on a concrete mechanism) | Exempt: lifts a value to a type |

- [ ] **Step 1: Add the rows, and see them fail:** `rate_equation`, `rate_equation_string` (with and without each keyword, through `_testhelper_keyword_call`), `parameters`, `dependent_parameter_values`, `check_valid_rate_equation`, `fitted_parameter_count`, `name`, `parameter_info`, `rate_equation_term_count`, `rate_equation_dedup_key`, concrete and singleton, plain and allosteric.
- [ ] **Step 2: Fix or exempt the sites, one group at a time,** re-running the file after each group. JET does not look past a dispatch, so new sites can appear in callees: treat them the same way. On Julia 1.13 a closure over a `@nospecialize` argument is untyped; use `Base.Fix2` or a helper taking the value as an argument.
- [ ] **Step 3: Run the Derivation tests:** `test_type_stability.jl`, then `Derivation/test_print.jl`, `test_law.jl`, `test_key.jl`, `test_joint_solve.jl`, `test_oracles.jl`, `test_rate_eq_derivation.jl`, then `test_compile_budget.jl`, one at a time. Expected: all pass, no pinned hash or count changed.
- [ ] **Step 4: Commit,** separating each independent fix, ending with "Check the Derivation functions for run-time dispatch".

### Task 5: `Enumeration` sites

**Files:**
- Modify: `test/test_type_stability.jl`; as needed `src/Enumeration/{moves,seeds,expand}.jl`

**Interfaces:** Consumes the table and exemption list from Task 3.

| Site | Expected action |
|:--|:--|
| `moves.jl:294, 325` (`_unbalanced_blocks`), `seeds.jl:56, 134, 159` (`_catalytic_topologies`), `moves.jl:641–649` (`_dead_end_child`) | Fix: closures and captured variables |
| `seeds.jl:602–612, 674, 703` | Fix what a type assertion on `pmap`'s result fixes; else exempt |
| `expand.jl:24` | Fix if small; else exempt with a reason |
| `filters.jl:503` (`_assert_atom_conserving`) | Exempt: error message |

- [ ] **Step 1: Add the rows, and see them fail:** `starting_mechanisms`, `expand_mechanisms` on `Vector{Mechanism}`, `Vector{AllostericMechanism}` and `Vector{Union{Mechanism, AllostericMechanism}}`, `init_mechanisms`, `seed_mechanisms`. Add `_pass1_record` for both concrete types if its Derivation sites are clean; otherwise it moves to Task 7.
- [ ] **Step 2: Fix or exempt the sites,** re-running after each group.
- [ ] **Step 3: Run `test_type_stability.jl`, `Enumeration/test_mechanism_enumeration.jl` and `test_compile_budget.jl`.** Expected: all pass, with exact child sets and pinned seed counts unchanged.
- [ ] **Step 4: Commit,** ending with "Check the Enumeration functions for run-time dispatch".

### Task 6: `Fitting`, `Data` and `DSL` sites

**Files:**
- Modify: `test/test_type_stability.jl`; as needed `src/Fitting/fit_rate_equation.jl`, `src/Data/rate_data.jl`

**Interfaces:** Consumes the table and exemption list from Task 3.

| Site | Expected action |
|:--|:--|
| `fit_rate_equation.jl` (18 sites: `_fit_loss`, `_optimization_problem`, the `fit_rate_equation` body, `_reference_scale`, `_scale_rate_constants`, the `FittingProblem` body, `test_rate_equation`) | Exempt the solver barrier and the code on the `@nospecialize` problem around it; fix a site a small change fixes |
| `rate_data.jl:36, 64` | Fix if small; else exempt with a reason |
| `read.jl` (21 sites, in `_parse_labeled_line`, `_bare_symbols_from_values`, `_bare_name`, `_raw_site`, `_step_parts`, `_peel_tag!`, `_raw_step!`, `_raw_form`, `parse_mechanism`) | Exempt: `Expr` walking |

`_fit_loss` is the barrier, but `loss!` itself gets no exemption: a dispatch inside it is the bug this PR exists to catch.

- [ ] **Step 1: Add the rows, and see them fail:** `loss!` on a plain and an allosteric problem, `fit_rate_equation`, `FittingProblem`, `test_rate_equation`, `RateData`, `parse_mechanism`.
- [ ] **Step 2: Fix or exempt the sites.**
- [ ] **Step 3: Run `test_type_stability.jl`, `Fitting/test_fitting.jl`, `Fitting/test_fitting_scale.jl`, `Fitting/test_fit_search.jl`, `Data/test_rate_data.jl`, `DSL/test_dsl.jl`.** Expected: all pass.
- [ ] **Step 4: Commit,** ending with "Check the Fitting, Data and DSL functions for run-time dispatch".

### Task 7: `identify_rate_equation`, and completeness

**Files:**
- Modify: `test/test_type_stability.jl`; as needed `src/identify_rate_equation.jl`

**Interfaces:** Consumes the table and exemption list from Task 3.

The 56 sites are `DataFrame` columns, `CSV.write`, `pmap` results and progress text. Fix what a typed collection, a `Float64.` conversion or a type assertion fixes: `_process_batch`'s and `_expand_parents`'s reads of `pmap` results are the likely ones. Exempt the rest by function, expected to include `_write_rows_csv!`, `_rows_to_dataframe`, `_progress`, `_default_save_dir`, `_select_best_row` and `_cv_model_selection`. If typing this code would mean rewriting the orchestration, stop and report. Do not change behavior, written files or progress text.

- [ ] **Step 1: Add the rows, and see them fail:** `identify_rate_equation` through `Core.kwcall` with `optimizer` a `CMAEvolutionStrategyOpt`, `IdentifyRateEquationProblem`, `_pass1_record` if not already present.
- [ ] **Step 2: Fix or exempt the sites.**
- [ ] **Step 3: Write the completeness test.**

```julia
@testset "every public function has a row" begin
    checked = Set(f for (_, f, _) in _TESTHELPER_CHECKED_CALLS)
    for mod in (EnzymeRates, EnzymeRates.Mechanisms, EnzymeRates.Derivation,
                EnzymeRates.Enumeration, EnzymeRates.Fitting, EnzymeRates.DSL,
                EnzymeRates.Data), n in names(mod)
        startswith(string(n), "@") && continue  # macros run at parse time
        x = getfield(mod, n)
        x isa Function && x ∉ checked && @test "$mod.$n" == "has a row"
    end
end
```

Run it with a row removed to see it fail, then restore the row. A row for a keyword call stores the function itself so this set still finds it.

- [ ] **Step 4: Prove the table catches a closure bug.** Introduce a boxed closure in an enumeration function that has a row; run the file; expected: that row fails. Revert. Record the result for the PR.

- [ ] **Step 5: Run `test_type_stability.jl`, `test_identify_rate_equation.jl`, `test_compile_budget.jl`.** Expected: all pass. Report the final number of exemptions and the time and peak memory of the type-stability file.

- [ ] **Step 6: Commit** as "Check identify_rate_equation for run-time dispatch and require a row for each public function".

### Task 8: Evidence and record

Run by the controller, not a subagent, because it needs the whole branch and the census scripts outside the repository.

**Files:**
- Modify: `ROADMAP.md`, `docs/src/developer.md`
- Delete: the spec and this plan

- [ ] **Step 1: Census.** Add a key column to `claude-scratchpad/fsn/census/census.jl` and `compare.jl` (`rate_equation_dedup_key(m).key`). Run it with all `REACTIONS` at `9c194361` in a throwaway worktree and at the final commit, one at a time. `compare.jl` must print `RESULT: PASS`. Remove the worktree.
- [ ] **Step 2: Timings.** Time `init_mechanisms`, `expand_mechanisms` and `_pass1_record` on the same mechanisms at both commits with a script in `$S/bench`, and record the sizes from a small real search.
- [ ] **Step 3: Record.** `ROADMAP.md` gains an idea to pin operation counts per law, the machine-independent check on a moderate `rate_equation` slowdown, and any entry the work retired is removed. `docs/src/developer.md` describes the two gates and the short law, in "What a fit costs".
- [ ] **Step 4: Full suite.** Check memory first, then, in the background: `timeout 1200 julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`. Expected: passes with pristine output. Report the time and memory the new tests add, against the 4 GB limit.
- [ ] **Step 5: Delete the spec and plan,** review the whole branch, and open a PR with the mutation runs, census, timings, exemption list and CSV sizes in its body. Commit `8151748d` holds the wrong regulator explanation, so the body recommends a squash merge. Denis merges.
