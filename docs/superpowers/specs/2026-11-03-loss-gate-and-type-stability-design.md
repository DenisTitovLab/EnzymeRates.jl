# Guard the fitting hot path and make the core machinery type-stable

## Intent

From #32 to #83, five months, every allosteric fit ran 14–27× slower per row because
`loss!` dispatched at run time, and no test noticed. The tests timed `rate_equation` alone,
`loss!`'s only allocation test used a uni-uni mechanism, and the suite's JET test looks for
errors, not for dispatch. `loss!` is the package's hot path, so a slowdown there must fail a
test on the PR that causes it. Denis also wants type stability checked on every exported
function and every public function the modules share, to catch other bugs of this kind.

## Goals

1. A failing test catches any allocation or slowdown in `loss!`, on every mechanism the
   `rate_equation` speed gate covers.
2. JET's optimization analysis (`report_opt`) reports zero on every exported function, every
   public function of a submodule and the hot path, apart from a short list of named,
   justified exemptions.
3. Every site JET reports on them is fixed or exempted, in one PR.
4. A public function added later fails the test until it has a row in it.

## Non-goals

- Timestamps in `progress.log`, dropping `rate_equation` from the CSVs, fewer restarts, and
  the idle-worker question. Each is a separate change.
- A rewrite of `identify_rate_equation`'s orchestration to avoid `DataFrame`, `pmap` or CSV.
- A comparison of whole search runs between versions.

## Anti-goals

- Any change to a derived law, a parameter name, an enumeration child set, a key or printed
  text. A fix that changes one of these is wrong.
- A slower `rate_equation`, or a higher count in `test/test_compile_budget.jl`.
- A new public name, or a new argument to one.
- A gate that fails because a machine is slow or busy, or because a Julia release changed
  inference with no change to the package.

## Constraints

- The suite runs on one Julia at a time, on 7.7 GB, and on CI's newest Julia. JET is already a
  test dependency.
- A rewrite of a whole implementation or module needs Denis's permission first. If a fix to
  `thermo.jl`, `key.jl` or `identify_rate_equation.jl` grows to that size, stop and ask.
- **Base (Denis):** #101 merges first, since it changes the versions of JET, CSV, Optimization
  and OrdinaryDiffEqFIRK. `docs/fit-speed-record` then rebases on `main`, and the sweep below
  runs again before triage. Denis does the merge.
- **Julia version (Denis):** zero reports is required on Julia 1.13 only. On other versions
  the test uses `@test_skip`, as `test/test_compile_budget.jl` does, and moving the pin is a
  deliberate change. The `loss!` gate runs on every version.
- Follow `.claude/CLAUDE.md`: TDD, `_testhelper_` names, 92 columns, and an update of the
  Developer page and `ROADMAP.md`. Delete this spec and its plan before merge.

## Measured starting point

- `loss!` per row exceeds `rate_equation` per call by 3–8.5 ns, or 1.3–3.3 `log` calls, on all
  48 gated laws. As a ratio this runs from 1.09 to 6.1, because a tiny law's call takes 0.6 ns.
- A `loss!` check takes about 0.14 s per mechanism, 7 s for all 48.
- `report_opt` gives 0 reports on `loss!` and on `rate_equation` of a singleton mechanism, and 8
  on `loss!` with the uninferable `metabolites` put back. `@inferred loss!` passes on that bug.
- Across every function in scope it reports 141 package sites on Julia 1.13, before #101. The
  sweep takes 17.5 s in a process that peaks at 1.8 GB.
- 141 is a floor. JET does not look past a run-time dispatch, so removing one can reveal sites
  in the callee.
- 37 of the 141 show only in a frame below the package, among them `_rref!`'s four. A filter on
  the last frame's module would hide them, so exemptions must not filter that way.

## Approaches

The `loss!` gate is the same in each. For type stability:

- **A: hot path only.** Zero reports on `loss!`, `rate_equation` and `FittingProblem`. Small,
  with no noise, but it leaves enumeration and derivation unchecked.
- **B: pin report counts.** Add per-function counts that may not rise, with no code fixes.
  Cheap, but refactors and Julia upgrades shift the counts, so it raises false alarms.
- **C: clean up, then require zero.** Fix the real sites, exempt errors and barriers by name,
  and assert zero. Strongest, and it touches derivation code, possibly for little speed gain.

**Denis chose C, all in one PR,** because he wants these bugs caught across the core
machinery, not only on the hot path.

## The `loss!` gate

Decided with Denis:

- It runs on every law of the `rate_equation` speed gate: the `MECHANISM_TEST_SPECS` entries
  and the larger laws in `_testhelper_largest_enumerated_mechanisms` and
  `_testhelper_large_ter_bi_mechanisms`. It extends the gate rather than standing beside it,
  so a law added to the gate gets the check. Mechanisms written inline elsewhere in
  `test/Derivation` are out of scope.
- It asserts `allocs == 0`, in relative and absolute mode.
- It asserts time per row `≤ rate_equation per call + 6 × log per call`. All three are timed
  in one process, so machine speed cancels out. This replaces the 300 ns bound in commit
  `b79f9eb3`. A strict ratio cannot hold, since each row also takes a `log` and is centered.

Kept from `b79f9eb3`: `RateData` stores `Rate` and each named metabolite as a
`Vector{Float64}`, and `test/Fitting/test_fitting.jl` tests column types and a
cross-validation fold.

## The type-stability test

A new file, `test/test_type_stability.jl`, included from `test/runtests.jl`, runs `report_opt`
on the functions below and requires zero reports. The exemptions live in one list in that file,
each with its reason. A new dispatch site therefore fails the test until someone fixes it or
argues for an exemption in review.

| Where | Functions |
|:------|:----------|
| Exported | `rate_equation`, `rate_equation_string`, `parameters`, `dependent_parameter_values`, `metabolites`, `fit_rate_equation`, `identify_rate_equation`, `FittingProblem`, `IdentifyRateEquationProblem` |
| `Derivation` | `check_valid_rate_equation`, `fitted_parameter_count`, `name`, `parameter_info`, `rate_equation_term_count`, `rate_equation_dedup_key` |
| `Enumeration` | `starting_mechanisms`, `expand_mechanisms` |
| `Fitting` | `test_rate_equation` |
| `DSL`, `Data` | `parse_mechanism`, `RateData` |
| `Mechanisms` | `catalytic`, `check_conformation_label`, `check_decodable_names`, `compile_mechanism`, `excludes`, `form_text`, `has_turnover_cycle`, `has_unbalanced_cycle`, `present_species`, `rapid_equilibrium_segments`, `segment_walk`, `step_text`, and the `Mechanism`, `AllostericMechanism`, `EnzymeMechanism` and `AllostericEnzymeMechanism` constructors |
| Internal | `loss!`, `init_mechanisms`, `seed_mechanisms`, `_pass1_record` |

Each runs on every form its methods accept: plain and allosteric, concrete and singleton, and
`expand_mechanisms` on plain, allosteric and mixed vectors.

**Claude's calls:**

- A completeness check walks `names` of the package and each submodule, and fails on a public
  function with no row, which is how "all" stays true.
- Macros are out, because they run at parse time and share their code with `parse_mechanism`.
- The internal functions are in. `loss!` and `_pass1_record` are the workers' hot paths, and
  dispatches above `init_mechanisms` and `seed_mechanisms` would otherwise hide them.

## Triage of the 141 sites

This triage is a starting point from reading each line. The builder confirms each site and
moves it if the code says otherwise. Exemptions name a function, and only these reasons qualify:
an error message, the solver barrier, `Expr` walking, a `DataFrame`, CSV or `pmap` boundary, and
lifting a value to a type.

| Sites | What they are | Action |
|:------|:--------------|:-------|
| `identify_rate_equation.jl` (56) | `DataFrame` columns, `CSV.write`, `pmap` results, progress text | Fix what a typed collection or a type assertion fixes; exempt the rest by function |
| `read.jl` (21) | Walks `Expr`, whose `args` is a `Vector{Any}` | Exempt |
| `fit_rate_equation.jl` (23) | The barrier that makes the solver compile once per optimizer, and code on the `@nospecialize` problem around it | Exempt, except a site a small change fixes |
| `filters.jl:503`, `messages.jl:22, 142` | Building an error message | Exempt |
| `codegen.jl:214` | `rate_equation` on a concrete mechanism lifts it to a type | Exempt |
| `engine.jl:35, 54`; `moves.jl:294, 325`; `poly.jl:199, 201`; `seeds.jl:56, 134, 159` | Local closures that capture a variable they reassign, or recurse | Fix |
| `moves.jl:641–649`; `normalize.jl:64–83` | A captured variable | Fix |
| `thermo.jl:146–160` (`_rref!`), `thermo.jl:274, 277` | Probably an unknown matrix type from the caller; `reduce(hcat, …)` in `normalize.jl` loses it | Fix the callers first, then re-run |
| `seeds.jl:602–612, 674, 703` | `pmap` results and a failure message | Fix what a type assertion fixes; else exempt |
| `engine.jl:338`, `expand.jl:24`, `print.jl:181`, `rate_data.jl:36, 64`, `rate_eq_derivation.jl:112, 113` | Mixed | Fix if a small change does it; else exempt with a reason |

## Done means

- Each gate fails when the bug it guards against is put back. At the least: the uninferable
  `metabolites` fails the `loss!` gate and the type-stability test, and a boxed closure in
  enumeration fails the type-stability test. The PR reports these runs.
- A census shows zero differences in laws, keys, parameter names, child sets and printed text
  between the rebased base and the final commit. **Claude's call:** it reuses
  `claude-scratchpad/fsn/census` with all its reactions, adds a key column and stays out of the
  repository. It samples parents, and the PR states that limit.
- `rate_equation`'s gate and `test/test_compile_budget.jl` pass unchanged.
- The full suite passes on Julia 1.13, and the PR reports the time and memory the new tests add,
  with the suite's peak staying below the 4 GB the machine allows. **Claude's call:** the peak is
  held and reported on 1.13, since 1.11 and 1.12 already sit at about 4 GB without these tests.
- The PR reports enumeration and PASS 1 timings before and after, so the benefit of the cleanup
  is measured, whichever way it comes out.

## Left to the builder

- How to express exemptions in JET, within the rule on last-frame filters above.
- The order of the fixes, and how to split them into commits.
- The exact fixtures in the type-stability test.
- Whether a site in a "Fix … else exempt" row is fixed or exempted, within the rule in that row.
