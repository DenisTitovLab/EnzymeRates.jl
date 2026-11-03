# One parameter type and a data-based fit search: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix issue #97 by giving the package one `Parameter` type family, then searching each fitted parameter around the concentrations and rates in the data, and naming parameters that end at a box edge.

**Architecture:** `StepParameter` replaces `_Column` and the four private types, so a joint-solve column and a parameter are the same value. `Derivation.parameter_info(m)` lists the fitted parameters, and `Fitting` computes each start and box from them and the mechanism. `rate_equation` and `loss!` do not change.

**Tech Stack:** Julia 1.12+, Optimization.jl, Test. Module rules: `.claude/CLAUDE.md`.

**Spec:** `docs/superpowers/specs/2026-11-02-data-based-fit-search-design.md`. Read it first; it holds the type definitions, `name`, `_log_centre`, `_scale_rate_constants`, doc wording and roadmap entries that this plan refers to and does not repeat.

## Global Constraints

- `rate_equation` and `loss!` do not change, and the `rate_equation` performance tests stay untouched.
- No names are exported. `parameters()` keeps returning Symbols.
- 92-character lines, 4-space indentation, `ABOUTME` header on new files, `_testhelper_` prefix on helpers under `test/`.
- One Julia process at a time: check `pgrep julia` first. Never park on a background monitor. Machine: 7.7 GB, no swap.
- Focused test run:
  `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<Module>/<file>.jl")'`
- Full suite, run before a task's commit when the task changes `src/`:
  `julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`
- A change to an `API`, a module docstring or a cross-module call updates `docs/src/developer.md` in the same commit, and `julia test/update_module_map.jl` rewrites the diagram.
- Comments are evergreen: no "new", "replaces", "used to". Remove no comment unless it is false.
- Do not loosen, reseed or delete a failing test. Only the two removals in the spec's "Removals Denis approved" are allowed. If a seeded identify test breaks in Task 3, record the test, its output and your diagnosis in the task report, leave the test as it was, and continue.
- Do not push, open a PR or bump the version.
- End every commit message with:
  ```
  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_015DwkaCNDvp2TjTtzLpaqqh
  ```

## Review Focus

- A table with integer concentrations or integer rates: the search works as with floats (Task 3 test).
- An explicit `lb`/`ub` that excludes the data-based centre: starts clamp into the box, and the fit runs and reports `at_bound` (Task 4 test).
- A cross-validation fold whose training rows hold a metabolite only at zero, while other folds do not: the fold fits through the fallback and does not throw (Task 3 test, on a table with that fold's rows).
- An L phantom reading its source's fit: its row reports the source's `at_bound` names, never `:L` (Task 5 test).
- Mutating the vector `parameter_info` returns must not corrupt the memoized solve (Task 1 test).

---

### Task 1: `StepParameter` replaces `_Column` and the private types

No behaviour changes: every parameter keeps its name and order.

**Files:**
- Modify: `src/Derivation/Derivation.jl` (delete `_Column`; `API` adds `StepParameter`, `RegulatorParameter`, `AllostericParameter`, `parameter_info`, `name`)
- Modify: `src/Derivation/names.jl`, `thermo.jl`, `rate_eq_derivation.jl`, `engine.jl:220`, `codegen.jl:51`, and every `c.edge` on a column in `admission.jl`, `guard.jl`, `messages.jl`, `limit_test.jl` (grep; other `.edge` fields stay)
- Modify: `src/Mechanisms/types.jl:863` (error text says "rendered regulator parameter names" in place of "rendered Kreg names") and its pins in `test/Mechanisms/test_types.jl:1074`, `test/Mechanisms/test_reaction_graph.jl:214,489`
- Modify tests, call forms only: `test/Derivation/test_names.jl`, `test_joint_solve.jl`, `test_admission.jl`, `test/Mechanisms/test_types.jl`, `test/Enumeration/test_mechanism_enumeration.jl`, `test/mechanism_definitions_for_test_enzyme_derivation.jl`, and any other file `grep` finds
- Modify docs: `docs/src/developer.md` (API table, parameter section, module map), `.claude/CLAUDE.md` chokepoint section

**Interfaces:**
- Produces: `Derivation.StepParameter(step::Int, slot::Int, inactive::Bool)`, `Derivation.RegulatorParameter(ligand::Symbol, inactive::Bool)`, `Derivation.AllostericParameter()`, all `<: Derivation.Parameter` (private); `Derivation.name(p, m[, kept]) -> Symbol`, concrete mechanisms only; `Derivation.parameter_info(m::_AnyMech) -> Vector{Parameter}`, a copy of what `_Solve` stores, in `parameters(m, FittedParameters())` order, after `check_valid_rate_equation(m)`.

- [ ] **Step 1: Write the failing tests** in `test/Derivation/test_names.jl`

```julia
@testset "parameter_info lists the fitted parameters in parameters() order" begin
    for spec in MECHANISM_TEST_SPECS
        m = _testhelper_lift(spec.mechanism)
        info = ER.Derivation.parameter_info(m)
        @test [ER.Derivation.name(p, m) for p in info] ==
              collect(parameters(m, FittedParameters()))
        @test all(p -> p isa ER.Derivation.Parameter, info)
        empty!(info)
        @test !isempty(ER.Derivation.parameter_info(m))
    end
end

@testset "parameter_info points at the steps of Uni-Uni" begin
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end)
    P = ER.Derivation.StepParameter
    @test ER.Derivation.parameter_info(m) ==
          [P(1, 2, false), P(2, 2, false), P(2, 1, false)]
    @test parameters(m, FittedParameters()) == (:k_ES_to_E_P, :k_ES_to_E_S, :k_E_S_to_ES)
end
```

Add an allosteric testset on the spec's fixture (`E + S ⇌ E(S) :: NonequalAI`, `E(S) <--> E(P) :: EqualAI`, `E(P) ⇌ E + P :: EqualAI`, `regulatory_site(2): R::NonequalAI`, `allosteric_regulators: R`, `catalytic_multiplicity: 2`) asserting the list ends `RegulatorParameter(:R, false), RegulatorParameter(:R, true), AllostericParameter()` and renders `:K_A_Rreg, :K_I_Rreg, :L`.

- [ ] **Step 2: Run `test_names.jl`, confirm failure** with `UndefVarError` on `parameter_info`.

- [ ] **Step 3: Implement in `names.jl` and `thermo.jl`**
  - The types and `name` as in the spec. `_state(m, p::StepParameter)` holds the logic of today's `_column_parameter`: `:None` in a plain mechanism, `:I` when `p.inactive`, `:EqualAI` for a step in an `:EqualAI` group, else `:A`.
  - `name(p::RegulatorParameter, am)` finds the ligand's tag in `am.regulatory_sites`, and renders `K_<ligand>reg` for `:EqualAI`, else the `K_A_` or `K_I_` form. `name(::AllostericParameter, ::AllostericMechanism) = :L`. Keep the docstrings' explanation of how names read.
  - `regulator_constant(am, lig, tag, inactive)` keeps its signature and `nothing` result, rendering through `name`.
  - Delete `_edge_constants`, `_column_parameter`, `Krapid`, `Kfor`, `Krev`, `Kreg`. In `_constant_sides`, drop the fallback for a slot with no class: use `kept[(g, slot, inactive)]` and delete that docstring sentence (`_check_representation` guarantees one flag per group).
  - `_Solve` gains `info::Vector{Parameter}`: the independents' `StepParameter`s, then the regulator parameters (same order as `_regulator_names`), then `AllostericParameter()` for an allosteric mechanism. `fitted` is the names of `info`, computed as today from `names[independents]` and `name` of the rest.
  - `parameter_info(m::_AnyMech) = (check_valid_rate_equation(m); copy(_memo_solve(m).info))`.
  - Keep `fitted_parameter_count` and `_edge_names` as they are, reading `StepParameter` columns.

- [ ] **Step 4: Migrate tests** with no expected name or value changed.
  - Add `_testhelper_step_parameters(m, k, inactive = false)` to the shared fixtures file: one `StepParameter` per slot of step `k`. It replaces each `_edge_constants` call; a `state` of `:I` becomes `inactive = true`, and any other state `false`.
  - `T(1, :None)` tests become equality and inequality tests of the new types' fields.
  - `Krapid(bind, :I)`-style assertions on a plain mechanism move to an allosteric fixture with the same expected names.
  - Assertions naming a combination `src/` never builds (a K on a steady-state step, a k on a rapid-equilibrium one) move to a fixture whose step has that flag when it gives the same expected name. Delete one only when no such fixture can give that name, and list each deletion in the task report.
  - The compiled-mechanism `MethodError` checks stay.

- [ ] **Step 5: Update `docs/src/developer.md`, `.claude/CLAUDE.md` and the module map.** Run `julia test/update_module_map.jl`.

- [ ] **Step 6: Run `test_names.jl`, `test_joint_solve.jl`, then the full suite.** Expected: pass with pristine output, including `test_boundaries.jl` and the chokepoint guard. Exit code 0.

- [ ] **Step 7: Commit:** `Derivation: one Parameter family, StepParameter as the solve's column label`

### Task 2: Scale through `parameter_info`, remove `RateConstants`

**Files:**
- Modify: `src/Fitting/fit_rate_equation.jl:293-304`, `src/Derivation/rate_eq_derivation.jl` (delete `RateConstants`, `_rate_constant_names`, `_edge_names`), `src/Derivation/Derivation.jl` (`API`), `docs/src/developer.md`
- Modify: shared fixtures file (add `_testhelper_edge_names`, unchanged from `_edge_names`)
- Modify tests: `test/Derivation/test_rate_eq_derivation.jl:687,1762-1812,2590-2632`, `test/Derivation/test_oracles.jl:478`, `test_print.jl:642`, `test_law.jl:482`, `test/Mechanisms/test_types.jl:2953-2955,3043`, `test/Mechanisms/test_accessors.jl:135`, `test/Fitting/test_fitting_scale.jl`

**Interfaces:**
- Consumes: `parameter_info`, `name`, `StepParameter`.
- Produces: private `Fitting._concrete(m)`, returning a concrete mechanism from a compiled or concrete one; `_scale_rate_constants(mechanism, params, factor)`, accepting both.

- [ ] **Step 1: Rewrite the tests.** The testset at `test_rate_eq_derivation.jl:1762` becomes "parameter_info gives what Fitting scales". Keep its scaling check over every `MECHANISM_TEST_SPECS` entry on concrete and compiled mechanisms, and keep its comment's promise. Its expected table shrinks to fitted names (Denis approved this), each computed as the names of the `parameter_info` entries on steady-state steps:

```julia
expected = Dict(
    "Uni-Uni" => (:k_ES_to_E_P, :k_ES_to_E_S, :k_E_S_to_ES),
```

Fill the other four entries (`RE Uni-Uni`, `MWC Dimer [AllostericEnzymeMechanism]`, `m_OnlyA_prod`, `LDH i-state NonequalAI 5-group`) from the old table's names that appear in `parameters(m, FittedParameters())`, read from today's output. `test_types.jl:3043` becomes `@test_throws msg EnzymeRates.Derivation.parameter_info(m)`.

- [ ] **Step 2: Run the file, confirm failure** on the scaling or table assertions.

- [ ] **Step 3: Implement** `_scaled`, `_concrete` and `_scale_rate_constants` as in the spec, and delete `RateConstants`, `_rate_constant_names` and `_edge_names`. Move the helper to the fixtures file, and update its callers.

- [ ] **Step 4: Update `developer.md` and the module map.**

- [ ] **Step 5: Run the full suite.** Expected: pass, exit code 0.

- [ ] **Step 6: Commit:** `Fitting scales through parameter_info; RateConstants goes`

### Task 3: Search around the data

**Files:**
- Modify: `src/Fitting/fit_rate_equation.jl`
- Create: `test/Fitting/test_fit_search.jl`, included after `test_fitting_scale.jl` in `test/runtests.jl`

**Interfaces:**
- Consumes: `parameter_info`, `_concrete`.
- Produces, private to `Fitting`:
  - `_concentration_scales(fp) -> Dict{Symbol, Float64}`: each of `Mechanisms.metabolites(fp.mechanism)` to the geometric mean of its nonzero concentrations in `fp.rate_data.data`. An always-zero metabolite takes the geometric mean of all positive concentrations of those metabolites in the table; with none at all, 1.0.
  - `_rate_scale(fp) -> Float64`: geometric mean of `|Rate|`, which validation guarantees nonzero (`exp(mean(fp.log_abs_rates))`).
  - `_log_centre(p, m, s::Dict{Symbol, Float64}, r::Real) -> Float64` as in the spec, with `_sum_log(s, mets) = sum(log(s[x.name]) for x in mets; init = 0.0)`.
  - `_log_centres(fp) -> Vector{Float64}`: `_log_centre` of each `parameter_info` entry, in order.

- [ ] **Step 1: Write the failing tests** in `test_fit_search.jl`:
  - **"Facts match the law":** for every `MECHANISM_TEST_SPECS` entry and for inline `E + A + B <--> E(A, B)`, `E + A + B ⇌ E(A, B)`, `F + A + A ⇌ F(A, A)`, `E(P, Q) <--> E + P + Q` mechanisms and one modifier copy (`S::Mod` with `modifiers:` declared), with random `f_M` and `g`: scale each metabolite by `f_M`, `Keq` by the product of `f_M` to its stoichiometry (products over substrates), and each fitted parameter `p` by `exp(_log_centre(p, m, sf, g * r) - _log_centre(p, m, s, r))`. Then `rate_equation(scaled) ≈ g * rate_equation(original)` with `rtol = 1e-9`.
  - **"Issue #97 regression":** `E + S ⇌ E(S)`, `E(S) <--> E + P`, K = 0.2e-6 in M, noiseless; the fit recovers K within 1%, and not at `exp(-15)`.
  - **"Unit independence":** one seeded absolute fit in M, mM and µM returns the same determined parameters after unit conversion, within `1e-3` log units (throwaway runs agreed to about 1e-6). A relative fit passes the same comparison after its scale.
  - **"Both modes":** a relative and an absolute fit each recover a steady-state law.
  - **"Two-metabolite law":** one fit recovers a law with `E + A + B <--> E(A, B)`.
  - **"Always-zero metabolite":** a table with `P` always zero fits with a finite loss, and `_concentration_scales(fp)[:P]` equals the fallback's exact value. A table with no positive concentration gives 1.0. A training subset where a metabolite is zero only in those rows fits too (Review Focus).
  - **"Integer columns":** a table of `Int` concentrations and rates gives the same centres as its `Float64` copy (Review Focus).
  - **"Search box":** with a stub optimizer that records `prob.u0`, `prob.lb` and `prob.ub`, defaults are `centre .± 15` and `u0` lies within them. Use `Random.seed!` and the `Optimization` solver already used in `test_fitting.jl` for real fits.

- [ ] **Step 2: Run the file, confirm failure** with `UndefVarError` on `_log_centre`.

- [ ] **Step 3: Implement** the functions above. `lb` defaults to `_log_centres(fp) .- 15.0`, `ub` to `_log_centres(fp) .+ 15.0`. `x0 = clamp.(centre .+ randn(np) .* 2.0, lb, ub)`, where `centre = _log_centres(fp)` is computed once per fit. Update the docstring: the defaults, and that `lb` and `ub` hold absolute log values. `FittingProblem`'s docstring gains the pointer the spec asks for in Task 6.

- [ ] **Step 4: Run the file, then the full suite.** Expected: new tests pass. If a seeded test elsewhere breaks, follow Global Constraints, and pick it up in the task report.

- [ ] **Step 5: Commit:** `Fit searches around the data's concentrations and rates`

### Task 4: `at_bound` in fit results

**Files:**
- Modify: `src/Fitting/fit_rate_equation.jl`
- Modify tests: `test/Fitting/test_fitting.jl:411,433,455`, `test/Fitting/test_fitting_scale.jl:79`, `test/Fitting/test_fit_search.jl`

**Interfaces:**
- Produces: `fit_rate_equation(...) -> (; params, loss, retcode, at_bound)`, `at_bound::Tuple{Vararg{Symbol}}`, names whose best point is within 0.01 log units of its `lb` or `ub`, read before the relative-mode scale; `()` when no restart reaches a finite loss.

- [ ] **Step 1: Write the failing tests:** with a stub optimizer returning a chosen `u`:
  - `u` at a known `ub`, within 0.009, names that parameter, and at 0.011 does not. The same at `lb`, and with explicit `lb`/`ub`.
  - An explicit box excluding the centre: every `prob.u0` lies in the box, and `at_bound` names the parameter when the stub returns the clamped start (Review Focus).
  - A fit inside the box gives `()`.
  - A relative-mode fit whose scale moves a reported k outside the box keeps `at_bound == ()`.
  - A `:NoFiniteLoss` fit gives `()`.
  - A `:ReferenceFailed` fit reports its best point's `at_bound`.
  - `fitted, loss, retcode = fit_rate_equation(...)` still destructures.
  - The `keys(...)` pins gain `:at_bound` last.

- [ ] **Step 2: Run, confirm failure.**

- [ ] **Step 3: Implement** in all three returns, keeping `at_bound` last, and document the field: "at bound" does not mean wrong, because a parameter no data constrain wanders across the box.

- [ ] **Step 4: Run the Fitting test files, then the full suite.**

- [ ] **Step 5: Commit:** `Fit results name the parameters at a box edge`

### Task 5: `at_bound` through identify and the CSVs

**Files:**
- Modify: `src/identify_rate_equation.jl` (`_Fit` at 392, `fit!` at 670, row at 707-722, failure row at 478-488, `_rows_to_dataframe` at 300-317, docstring at ~168)
- Test: `test/test_identify_rate_equation.jl`

**Interfaces:**
- Consumes: `fit.at_bound`.
- Produces: `_Fit.at_bound::Vector{Bool}`, flags at the same canonical positions as `values` (`Bool[x in fit.at_bound for x in _fit_order(r)]`); rows and CSVs gain `at_bound`, a `String` of names joined with `";"`, `""` when none, `missing` on a failed row.

- [ ] **Step 1: Write the failing tests:**
  - Two mechanisms sharing one fit report their own names for the same positions: names come from each candidate's `_fit_order(c)` with the fit's flags.
  - An L phantom reads `[fit.at_bound; false]`, so `:L` never appears (Review Focus).
  - The search CSVs, `loocv_results.csv` and `best_equation.csv` each have an `at_bound` column after a short `identify_rate_equation` run with `save_dir`, and `best_equation.csv` reports the all-rows fit, not a fold's.
  - A failed mechanism's row has `missing`.

- [ ] **Step 2: Run, confirm failure.**

- [ ] **Step 3: Implement.** `_rows_to_dataframe` feeds all CSVs, so one column there covers them. `src/` has no CSV reader of these files, confirmed by `grep -n "CSV.read\|CSV.File" -r src`.

- [ ] **Step 4: Run `test_identify_rate_equation.jl`, then the full suite.**

- [ ] **Step 5: Commit:** `identify_rate_equation reports at_bound in its CSVs`

### Task 6: User docs and roadmap

**Files:**
- Modify: `docs/src/getting_started.md` near the `Keq` paragraph, the `FittingProblem` docstring, `ROADMAP.md`

- [ ] **Step 1: Add the spec's docs text verbatim.** Add the note that the box spans about 6.5 decades either side of each centre, and explain `at_bound`. Add the pointer from `FittingProblem`'s docstring. `docs/src/getting_started.md:129` keeps its destructuring.

- [ ] **Step 2: Add the two `ROADMAP.md` entries in the spec's wording** (Known issues; Decided against). If `ROADMAP.md` or `docs/src/roadmap.md` lists issue #97 or this bound, retire that entry.

- [ ] **Step 3: Run `test/test_boundaries.jl` focused and any doctest the docs file carries.** Expected: pass.

- [ ] **Step 4: Commit:** `Docs: units, the search box and at_bound; roadmap entries`

### Task 7: Generality check and clean-up

**Files:**
- Create in the scratchpad only: `one_off_generality.jl`
- Delete: the spec and this plan

- [ ] **Step 1: Run one-off fits, with no repo changes:** one with 5% multiplicative noise and one on a law of about 20 fitted parameters from `_testhelper_largest_enumerated_mechanisms`, each with today's default box (explicit `lb = fill(-15, np)`, `ub = fill(15, np)`) and with the default box, same seed and budget. Report loss and `at_bound` for both, in the task report. A worse result is a finding to report; do not change the design to chase it.

- [ ] **Step 2: Re-read changed `src/` files for dead code,** such as an unused helper or a stale docstring mention of `_Column`, `Krapid`, `RateConstants` or `_edge_constants`: `grep -rn "_Column\|Krapid\|Kfor\|Krev\|Kreg\b\|RateConstants\|_edge_constants\|_column_parameter\|fitted_parameters" src docs/src test .claude/CLAUDE.md`. Expected: no matches other than the `test_rate_eq_derivation.jl` testset wording you chose.

- [ ] **Step 3: `git rm` the spec and this plan,** as CLAUDE.md's Workflow requires.

- [ ] **Step 4: Run the full suite one last time.** Expected: pass, exit code 0.

- [ ] **Step 5: Commit:** `Remove the implemented spec and plan`
