# Deep modules and simplification: design

Main at 4134a5ff (after PRs #85–#95). One PR, version 0.11.0.

## Intent

PRs #85–#95 rebuilt the package around the reaction graph. They left two problems:
- **Code with no current use.** Several features and special cases are no longer needed.
- **No boundaries.** Every file reaches into every other file's internals, so working in one part means knowing all of it. The current coupling is about 117 function names crossing into the mechanism types, 690 direct field accesses on them outside their files, and about 50 internals crossing among derivation, enumeration, fitting and identify.

This PR deletes the dead code and splits the package into **deep modules**, after Ousterhout's *A Philosophy of Software Design*. Each module hides a lot of code behind a few functions. Someone, human or agent, working in one module needs to know only the other modules' short interfaces. A test enforces the interfaces.

## Goals

- **Deletions.** Delete the src code listed in §1.
- **Four submodules.** Split the package into `Mechanisms`, `Derivation`, `Enumeration` and `Fitting`. Identify, the macros and the exports stay in the top-level module.
- **Small interfaces.** Each module exposes a few functions (§4). The mechanism types are reached only through about 30 accessors (§5).
- **A boundary test.** A parse-only test fails on any cross-module name outside a module's API list, and on any mechanism field access outside `Mechanisms` (§6).

## Non-goals

- **No deletions from tests, docs or src comments.** The deletion goal counts src code only. Tests of deleted code are rewritten or retired, and §8 lists each one. Docstrings go with the code they document.
- **No new features**, apart from the `reference` anchors that replace the peak scale.
- **No change to `rate_equation`.** Its generated code, its speed and its zero allocations stay as they are.

## Anti-goals

Any one of these makes the work a failure, even if the tests pass:
- **A changed law or child set.** A derived rate law, a parameter name or an enumeration child set changes in a commit that claims to be a refactor.
- **A slower `rate_equation`.** It allocates or slows down; CLAUDE.md's bound holds.
- **A slowdown nobody approved.** Enumeration, the guard or identify gets measurably slower (more than 5% on the benchmark set) without Denis approving it.
- **A dead-end hand-written mechanism.** The guard rejects a hand-written mechanism without the message and fix suggestions it gives today. Writing mechanisms by hand stays first-class.
- **A flat package in disguise.** The module split leaves code reaching sibling internals through qualified names.

## Constraints

- **CLAUDE.md applies:**
  - TDD;
  - the smallest change;
  - no rewrite of a working implementation without permission;
  - pristine test output;
  - `_testhelper_` prefixes;
  - ABOUTME headers;
  - one Julia process at a time (7.7 GB RAM, no swap).
- **Julia compatibility:** 1.10 through 1.13. The `public` keyword needs 1.11, so API lists are plain constants, not `public` declarations.
- **No Julia privacy.** Julia has none, and ExplicitImports.jl treats every submodule under one package root as internal. The boundary test is the only enforcement.
- **`Core` is taken.** No submodule may be named `Core`.

## Approaches

**A. Layer by layer.**
1. Deletions.
2. Rewrite the consumers onto `Mechanisms` accessors.
3. Cut the submodules and add the boundary test last.

This is simple, but nothing enforces the boundaries during the long rewrite.

**B. Boundary test first, then ratchet (chosen).**
1. Deletions.
2. Cut the four submodules with mechanical moves only, and add the boundary test at once. Each API list starts as today's crossings and the field-access ban list starts empty, so the test passes.
3. Shrink each list towards the target in §4–§5. Each shrink turns the test red, and a commit turns it green again.

The test drives the refactor, progress shows as shorter lists, and every commit is small and checkable.

**C. Rebuild each module behind its new interface.** Rejected. It throws away working implementations, which CLAUDE.md forbids, and leaves the census pins nothing to compare against stage by stage.

Denis chose B.

## 1. Deletions

Each deletion below is its own commit. A commit that changes behaviour says so in its message.

**Limit test: keep only the general methods.**
- Delete `_limit_by_aggregation` (`src/admission.jl:1138–1201`) and `_law_on_face` (`1238–1271`).
- Delete the A2 shortcut in `_face_stage`.
- Keep Case I.
- Evidence (probe 10): verdicts identical on 318,105 faces. Open faces cost 1.4 ms in place of 0.26 ms. A2 never fired in the corpus.

**Peak scale: dropped; relative mode requires reference anchors.** Behaviour change and API break.
- `FittingProblem` and `IdentifyRateEquationProblem` take `reference = [(conc...) => turnover, ...]`. A metabolite not listed is at zero.
- With `absolute = false`, `reference` is required.
- The fit scales the steady-state rate constants through the existing anchor code (`rescale_parameter_values` with anchors).
- A zero, non-finite or wrong-sign predicted turnover at a reference point makes the fit fail. Failing means no throw: `fit_rate_equation` returns loss `Inf` with retcode `:ReferenceFailed`, so identify drops the mechanism from the beam (`_ingest!` skips non-finite losses) and never aborts.
- The cross-validation folds pass `reference` like any fit. A failed reference returns the parameters unscaled, and a fold scores only the parameters with the scale-free relative loss, so a fold never aborts on it.
- Nothing replaces `peak_at_range_edge` in the fit result, `_Fit` or the CSV.
- Delete:
  - `DEFAULT_PEAK_RANGE`, `_peak_range` and the `peak_range` keyword;
  - `LogConcentrationTurnover`, `_maximize_on_box`, `_coordinate_ascent!`, `_line_maximum!`, `_peak_turnover` and `_finite_or_minus_inf`;
  - the `peak_at_range_edge` field and CSV column;
  - the no-anchor method of `rescale_parameter_values`.

**Fast split counter.**
- Delete:
  - the plain `_split_gain_test`;
  - `_partition_independent_count`, `_thermodynamic_constraints`, `_cycle_basis`, `_add_step_column!` and `_count_kind` (`src/thermo.jl:476–605`).
- The allosteric gain test, which counts the built child with `_independent_param_count`, becomes the only gain test.
- Evidence (probe 12c): identical children, and no measurable slowdown of `expand_mechanisms`.

**Class 5.**
- Delete `_assert_support`, `_support_verdict`, `_vanishing_face`, class 5's message branch, the `_DENOMINATOR_VANISHES` fault, and the calls at `src/engine.jl:351,356`.
- Delete `_law_branches` if nothing else uses it.
- Evidence (probe 13):
  - No mechanism the guard accepts has a 0/0 face: 0 of 2,293 enumerated mechanisms and 0 of 9,396 single-group RE/SS flip variants.
  - Every class 5 test fixture also fails class 3.
  - `loss!` already scores a NaN prediction as a sign mismatch.

**Per-regulator multiplicities `I(1, 2)`.**
- Stored but never read.
- The DSL rejects the syntax with a parse error.
- Delete the `RegulatorMults` multiplicity list (`src/dsl.jl:28`, `src/types.jl:444`).
- Update the ROADMAP entry "Make per-regulator multiplicities such as `A(1, 2)` work": the syntax now returns with the feature.

**Dead code.**
- `catalytic_mechanism` and `catalytic_multiplicity` on `AllostericEnzymeMechanism` (`src/types.jl:1416–1419`).
- The unreachable L check (`src/identify_rate_equation.jl:43–47`).

## 2. Refactor merges

These are from the simplification scan (items 3, 6 and 7), each a behaviour-preserving commit:
- **Item 3. Admission face plumbing.** Drop `_admission_face`'s duplicate plumbing, share `_first_fault`, merge `_tags` with `_solve_tags`, and inline the four one-call helpers.
- **Item 6. Merge plain/allosteric pairs** where the bodies coincide: `_with_law`, `_checked_mechanism`, and the two raw constructors.
- **Item 7. One exclusion rule.** `_excluded_by` serves `_state_conformations`, `_site_active_states`, `_reg_K` and `_kreg_params`.

Four more scan items land in §5 as `Mechanisms` derived views, so they are done once:
- **Item 2:** one conformation-steps view replaces `_active_mask`, `_inactive_mask`, `_active_law`, `_inactive_law`, `_state_mask` and the admission loop.
- **Item 5:** the shared graph walks `_reached` and `_has_unbalanced_cycle`.
- **Item 8:** one turnover-cycle test replaces `_has_turnover_cycle` and `_re_turnover_cycle`. **Gate:** first add tests that pin both functions' current behaviour, including malformed input.
- **Item 10:** one RE-segment walker replaces `_segments`, `_fast_segments` and `_bottomless_segment`. **Gate:** benchmark the guard and enumeration, and drop the merge if it is slower.

## 3. Modules and layering

```
Mechanisms ← Derivation ← Fitting
                        ← Enumeration
top level (identify, macros, exports) uses all four
```

| Module | Contents |
|---|---|
| `Mechanisms` | the reaction and mechanism types, the written form (`Raw*`), canonical form, `compile_mechanism` and the `Sig` types, naming rules (`names.jl` plus `_kreg_params`, `_reg_K` and `_kept_members` moved in), DSL parsing and printing, the accessors and shared derived views (§5) |
| `Derivation` | the guard (classes 1–3, limit test, admission, messages; `admission.jl` split into four files), thermo, poly, engine, normalize, print, key, codegen, `rate_eq_derivation` |
| `Enumeration` | `mechanism_enumeration.jl` split into files for moves, seeds and filters, plus identify's base tier (`_base_tier`, moved in) |
| `Fitting` | `fitting.jl` plus identify's cross-validation fold (`_cv_fold_loss`, `_fold_problem`, moved in) |
| top level | `identify_rate_equation.jl`, the three macros (they call the guard, so they can't live in `Mechanisms`), the exports, the precompile statements |

**Layout.** Each module has a folder `src/<Module>/` containing `<Module>.jl`, which declares the module and includes its files. Tests mirror this layout in `test/<Module>/`.

**Layering is enforced:** `Mechanisms` names no sibling; `Derivation` names only `Mechanisms`; `Enumeration` and `Fitting` name only `Mechanisms` and `Derivation`. The moves needed to make the graph acyclic:
- `compile_mechanism` goes to `Mechanisms`.
- `_has_unbalanced_cycle` goes from enumeration to `Mechanisms`.
- `_raw_form` stops being called by enumeration directly.
- The naming rules move into `Mechanisms`.

## 4. Interfaces

A cross-module call is always qualified (`Derivation.has_valid_rate_equation(m)`). Code inside a module stays bare. Names in an API list carry no leading underscore; every underscore name is private to its module. Final names may change during planning; the counts may not grow without Denis.

**`Derivation`:**

| Function | For | Today |
|---|---|---|
| `has_valid_rate_equation(m) -> Bool` | Enumeration, identify | `_is_well_posed` |
| `check_valid_rate_equation(m)` (throws the guard's message) | macros, identify | `_check_well_posed` |
| `parameter_count(m)` | Enumeration's gain test | `_independent_param_count` |
| `rate_equation_key(m) -> (; key, parameters, fit_key, fit_parameters)` | identify's fit dedup | `_canonical` on the derived law, plus `_is_phantom`. `key` and `parameters` are the law's own key and its canonical parameter order. `fit_key` and `fit_parameters` name the fit the mechanism reads. They equal `key` and `parameters` except for an L phantom, where they are the active law's key and order (without `L`). Identify keeps today's phantom behaviour: its own key in the beam, the CSV and `n_params`; demotion when the source fit throws; the `L_phantom` column; the cross-validation tie-break. The plan settles how the leaf-cap log line gets `capped`. |
| `rate_equation_complexity(m)` | identify's `eq_complexity_filter` | `_eq_complexity` |
| `scale_rate_constants(m, params, factor)` | Fitting's reference scaling | `_ss_rate_constant_names` plus `_scale_rate_constants` |
| already public: `rate_equation`, `rate_equation_string`, `fitted_params`, `parameters`, `all_parameters`, `metabolites` | Fitting, identify | |

- **Identify drops `_with_law`, `_solve`, `_law_text` and `_Canonical`.** It calls the public `fitted_params(m)` and `rate_equation_string(m; coefficients = true)` on the concrete mechanism, plus `rate_equation_key`. `rate_equation_string` on a concrete mechanism compiles nothing. `fitted_params` on a concrete mechanism today goes through `compile_mechanism` (`src/rate_eq_derivation.jl:35`), which generates a type per mechanism; its concrete method changes to read the joint solve (`_solve`) instead. Each call derives the law again. **Gate:** measure identify's PASS 1 before and after. If the slowdown exceeds 5%, add a private derivation cache inside `Derivation` rather than widening the interface.
- **`has_valid_rate_equation` excludes class 5** (deleted) and runs no derivation.
- **`MAX_RATE_EQUATION_TERMS` (5,000) stays** as the hard derivation cap. Identify's filter default (337) is a separate number. Probe 14: every one of the 557 ter-bi mechanisms over 337 terms has at most 20 parameters, so `max_param_count` does not cover the filter.

**`Enumeration`:**

| Function | For |
|---|---|
| the base tier of a reaction | identify |
| the children of one mechanism | identify |

- **The base tier** is what identify fits first: `init_mechanisms` (no required regulators) or `seed_mechanisms` (with them), then the base tier. Non-degenerate seeds are kept; each degenerate seed is replaced by its non-degenerate RE→SS flip children. It returns its expansion failures beside the mechanisms.
- **The children** are today's `expand_mechanisms` on one parent.
- `init_mechanisms`, `seed_mechanisms`, `expand_mechanisms`, `_degenerate` and the hyperbolic check become private, reached by tests through qualified names.
- **Enumeration calls only `Derivation.has_valid_rate_equation` and `Derivation.parameter_count`.**
  - Inside the flip, split and tagging searches it currently runs guard stages on bare flag vectors: `_faces_pass`, `_bottomless_segment`, `_admits`, `_stops_inactive_turnover`.
  - These partial stages are not the full guard:
    - the flip search leaves out admission's face condition;
    - the tagging search leaves out class 3;
    - the split search checks only bottomless segments.

    A set that passes a partial stage blocks its supersets. So with the full `has_valid_rate_equation`, a set that fails a stage left out today becomes a failed set, and the search tries its supersets. Child sets can change, not only speed.
  - **Decision (Denis, option a):** each search builds its candidate and asks `has_valid_rate_equation` (and `parameter_count` for the gain test).
  - **Gate, one search per commit:** run the census and the speed benchmark (`expand_mechanisms` and `init_mechanisms` on the reference reactions). A census difference, or a slowdown above 5%, goes to Denis with the numbers before that commit lands. With no difference, the commit is a pure refactor.
  - **Prior evidence for the flip search** (phase 4 final review, `review4/measure.md` §7): adding the face condition to the flip predicate gained and lost 0 of 44,801 flip children, with the whole census identical. The flip move took about 2× as long; an `expand_mechanisms` round did not change. The tagging and split searches are unmeasured.
  - For plain mechanisms the flip search's stages already equal the full guard, since a flip cannot change class 1 (phase-4 fix wave, FW-09).
  - Enumeration's structural shortcuts (RE-segment count, flux-carrying groups, chain flanks, dwell gauge, hyperbolic check, degeneracy) stay private. They read the mechanism only through `Mechanisms` accessors.

**`Fitting`:**

| Function | For |
|---|---|
| `fit_rate_equation` | users, identify |
| `cross_validate(m, data; ...)` | identify: the held-out loss of each group |
| a validated rate-data constructor | `FittingProblem` and `IdentifyRateEquationProblem` both validate the table through it |

`loss!`, `_rate_table`, `_cv_fold_loss` and `_fold_problem` become private.

## 5. `Mechanisms`: the graph and its accessors

**The model.** A mechanism is a graph:
- **Nodes** are enzyme species (forms).
- **Edges** are steps. Each step carries the metabolites it takes up and gives off, its RE/SS flag and its kinetic group.
- **Allosteric mechanisms** add a tag per group, the catalytic multiplicity and the regulatory sites.

The encoding stays as it is. What gets hidden is every storage decision:
- the `Counts` bit packing;
- metabolite-table indices;
- edge direction (stored physical-forward);
- form numbering;
- canonical order;
- the edge masks.

**Draft accessors.** Forms and steps are integer handles, and a set of steps is a `BitVector` over handles. Final names and signatures are settled in the plan, from the inventory of current uses.

| Layer | Accessors |
|---|---|
| Species | `species(m)`, `bound(m, x)`, `conformation(m, x)` |
| Steps | `steps(m)`, `ends(m, k)`, `taken_up(m, k)`, `given_off(m, k)`, `is_rapid_equilibrium(m, k)`, `step_kind(m, k)` (binding, isomerization, chemistry) |
| Groups | `group(m, k)`, `members(m, g)`, `ngroups(m)` |
| Metabolites | an opaque multiset with `copies`, `names`, `isempty`, `∪`, `∩`; `substrates(m)`, `products(m)`, `modifiers(m)` |
| Allosteric | `catalytic(am)`, `tag(am, g)`, `multiplicity(am)`, `regulatory_sites(am)`, `conformation_steps(am, :active \| :inactive)` |
| Derived | `re_segments(m, steps)`, `segment_graph(m, steps)` (the SS graph between segments), `has_turnover_cycle(m, steps)`, `present_species(m, steps)` |
| Write | the written form (`RawMechanism`, `RawStep`, `RawForm`, `RawAllostericMechanism`, `RawSite`), the checked constructors, `with_flags`, `regroup`, `with_tags` |

- **Rule (agreed with Denis).** A derived view lives in `Mechanisms` only when at least two modules use it. Otherwise it stays private to its module. So the guard's networks and faces, spanning trees, the cycle basis, flux-carrying steps and the dwell gauge stay where they are.
- **What the budget of about 30 covers.** It covers the graph accessors in the table above. `Mechanisms`' API list also carries the types (`Mechanism`, `AllostericMechanism`, `EnzymeReaction`, the metabolite and `Sig` types, the `Raw*` types), the reaction accessors (`substrates(rxn)`, `regulators`, …), `compile_mechanism` and the DSL entry points the macros call. These don't count towards the 30, but the plan enumerates them and they may not grow beyond what exists today.
- **Accessors must not allocate** where today's field access doesn't. In particular, the metabolite multiset keeps the bit encoding inside.
- **Scope.** About 650 consumer sites move onto the accessors: enumeration 341, admission 221, thermo 63, engine 46, other derivation files 14.

## 6. The boundary test

`test/test_boundaries.jl` parses `src/` with `Meta.parseall`, without loading the package. It fails when:
1. **A sibling's name is outside its API list.** The check applies to `Sibling.name` in any module, and to any name the top level uses from a submodule. Each module declares its list as a constant in `<Module>.jl`, for example `const API = (:has_valid_rate_equation, :check_valid_rate_equation, ...)`.
2. **The layering is broken:** a module names a sibling it may not depend on (§3).
3. **A mechanism field is accessed outside `Mechanisms`.** The test bans `x.f` outside `src/Mechanisms/` for every field name of `Mechanism`, `AllostericMechanism`, `Form`, `Edge`, `Met`, `RegulatorySite` and the `Raw*` types. The test lists each unrelated `x.f` it allows (`data.group`, `fit.params`, …) by file.
4. **A cross-module name starts with an underscore.**

The scanner prototype is at `scratchpad/api/deps.jl` (`Meta.parseall`, definitions and references per file).

## 7. Order of work (one PR, approach B)

1. **Deletions** (§1), one commit each, with behaviour changes flagged in the message.
2. **Merges** (§2, items 3, 6 and 7).
3. **Cut the submodules.** Mechanical moves and qualified names only. Move the tests into `test/<Module>/` and rewrite their internal references to qualified names by script (about 2,239 references to 307 internals). Add the boundary test, with each API list set to today's crossings and an empty field ban list. Everything passes.
4. **Ratchet the interfaces**, one shrink per commit, test red then green:
   1. `Fitting`;
   2. `Enumeration` (the base tier moves in);
   3. `Derivation` for identify (`rate_equation_key`, public calls);
   4. `Derivation` for enumeration (`has_valid_rate_equation`/`parameter_count`), one search per commit, each behind the census and benchmark gate (§4).
5. **Ratchet `Mechanisms`.** Add the accessors and derived views, ban one field name at a time and rewrite its uses, and merge scan items 2, 5, 8 and 10 (with their gates) as derived views.
6. **Docs, ROADMAP, CLAUDE.md, version 0.11.0.**

**Census pins.** Every refactor commit (all but the flagged behaviour changes) keeps these unchanged:
- the enumeration census counts;
- the exact child sets in `test_mechanism_enumeration.jl`;
- the derived laws of `MECHANISM_TEST_SPECS`;
- the golden files;
- the `rate_equation` performance test.

## 8. Tests that change

Tests are rewritten, never deleted silently. Each change is listed here, and any not listed here goes to Denis.
- **Peak scale.** Retire the peak-search testsets in `test/test_fitting_scale.jl` (lines 19–294) with the feature. Rewrite "rescale_parameter_values: peak scale". Add reference-anchor tests: required in relative mode, unlisted metabolites at zero, failed fit on zero, non-finite or wrong-sign turnover. The two other `rescale_parameter_values` testsets stay.
- **Class 5.**
  - "every spec's denominator keeps a monomial on every face" becomes a numeric check: `rate_equation` is finite at random parameters with each face's metabolites at zero.
  - The class 5 fixtures (`test_admission.jl:1079`, `1269`, `3100–3230`) become assertions that the guard rejects them.
  - Add the same finite-on-every-face check over a sampled enumeration census.
- **Limit test.** `test_admission.jl:1279`, which pins `_LIMIT_UNTESTED` for A2, changes to the general method's verdict.
- **Fast counter.** Tests of `_partition_independent_count` and the plain gain test retire. The built-child gain test must give identical children (probe 12c).
- **`I(1, 2)`** becomes a parse-error test.
- **Qualified names.** Every test reaching an internal uses its qualified name (`EnzymeRates.Derivation._solve`).

## 9. Docs and process

- **Developer page** (`docs/src/developer.md`): add the module map, the interfaces and the boundary test, and remove the peak scale and class 5.
- **User docs:** `reference` replaces `peak_range`.
- **CLAUDE.md:**
  - the API Design section (the `EnzymeRates.init_mechanisms` path becomes qualified);
  - the focused-test command, for the new test layout;
  - a line on the boundary test.
- **ROADMAP.md:** the `I(1, 2)` entry.
- **Version:** 0.11.0, the only bump in the PR.

## Evidence

Scripts and logs are in the session scratchpad under `probe_lim/`.
- **Probe 10:** general versus special-case limit methods give identical verdicts on 318,105 faces.
- **Probe 12c:** the built-child split counter gives identical children, with no expand slowdown.
- **Probe 13:** class 5 rejects nothing the guard accepts (2,293 enumerated mechanisms and 9,396 flip variants).
- **Probe 14:** the complexity filter is independent of `max_param_count`. Ter-bi has 557 of 13,260 mechanisms over 337 terms, all with at most 20 parameters; the maximum is 756,000.
- **Probe 15:** distinct laws by key versus by printed string:

  | Reaction | Mechanisms | Distinct strings | Distinct keys |
  |---|---|---|---|
  | bi-bi | 2,997 | 2,997 | 2,643 |
  | ping-pong | 21,392 | 21,014 | 16,056 |
  | bi-bi, R, n=2 | 2,299 | 2,299 | 2,079 |

## Left to the builder

- **Accessor names and signatures**, chosen from the inventory. The count stays near 30.
- **Enumeration's two function names.** Prefer today's names where they still fit.
- **The split of `admission.jl` and `mechanism_enumeration.jl` into files.**
- **Benchmark sets** for the three gates (enumeration on built candidates, identify's PASS 1, scan item 10).
- **Generated code and precompile.** Name scoping inside `rate_equation`'s generated code (bodies resolve in `Derivation`), and qualifying the precompile statements.

## Calls I made (for Denis to check)

1. **Identify stays in the top-level module** rather than a fifth submodule: nothing calls it, so a module would hide nothing.
2. **`check_valid_rate_equation(m)`** sits beside `has_valid_rate_equation(m)`: the macros and identify need the guard's message.
3. **`scale_rate_constants(m, params, factor)`** is `Derivation`'s export to `Fitting`, so `Fitting` never learns which constants are rate constants.
4. **A validated rate-data constructor** is `Fitting`'s third export, so identify validates its table without `_rate_table`.
5. **No leading underscore on API-list names**, and the boundary test rejects any cross-module underscore name.
6. **The 5% threshold** for the benchmark gates.
7. **Scan items 2, 5, 8 and 10 land as `Mechanisms` derived views** (§5) rather than as standalone merges.
8. **A bad reference turnover gives loss `Inf`, retcode `:ReferenceFailed` and unscaled parameters**, not a throw. The cross-validation folds pass `reference` like any fit, and nothing replaces `peak_at_range_edge` (§1).
9. **The budget of about 30 covers only the graph accessors** (§5).
10. **`fitted_params` on a concrete mechanism reads the joint solve** rather than compiling (§4).
