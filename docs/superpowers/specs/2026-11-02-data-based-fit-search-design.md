# One parameter type and a data-based fit search

Fixes issue #97. One PR with two stages: unify how the package represents a parameter, then
use that representation to search fits around the data.

## Intent

`fit_rate_equation` searches every constant in `log(p) ∈ [-15, 15]` from
`x0 = 2 * randn()`, whatever the data's units. For data in M, this floor (3.06e-7) cuts off
real binding constants, and the fit gives no warning. James's NADPH K sat exactly on it.

The search should start and run near the concentrations and rates in the data, so the
choice of unit stops mattering, and a fit that ends on the box edge should say so.

A correct start needs to know what each constant is. Today that knowledge is split: private
types render names, and `Fitting` only sees the resulting Symbols. Parameters are a core
concept of the package, alongside mechanisms and `rate_equation`. They should have one
representation, from which names, fitting information and future uses are all computed.

## Goals

- One `Parameter` type family, shared across modules, replaces every other description of a
  constant. Everything else about a constant is computed from it.
- The fit recovers a constant the data determine, in M, mM or µM alike.
- Each kind of constant gets its own start and box, by dispatch on its type.
- One rule serves relative and absolute mode.
- Fit results and CSVs name each constant at a box edge.

## Non-goals

- `lb`/`ub` or a box width in `identify_rate_equation`. Fewer choices make the package
  easier to use. The parameter count varies across the beam, so per-constant lists could
  not apply anyway.
- A declared concentration unit or a physical-plausibility warning.
- Reporting that the data only limit a K from one side (see Roadmap).
- Changing what `parameters()` returns to callers: names, as Symbols.
- Parameter types inside the joint solve. `_Column` stays its private position label.

## Anti-goals

- Any change to `rate_equation` or `loss!`, including their performance. The rate law takes
  its names at compile time, as it does now.
- Hard limits at physical values. They recreate the bug for real enzymes just outside them.
- Parsing constant names to recover facts about them.
- A second family of types, or fields that repeat what the mechanism's steps hold.

## Approaches

**The search** (Denis chose the last):

| Approach | Trade-off | Outcome |
|---|---|---|
| Global box ±ln(1e9) | Two lines, but 1e-9 means 1 nM only in M; K = 0.2 nM in M failed in 20 of 20 restarts | Rejected |
| Physical ranges per class (K 1 nM–100 mM, k 0.1–1e9) | Walls are meaningful, but need known units; K = 0.2 nM, K = 300 mM and absolute k = 0.02/s failed in 20 of 20 restarts | Rejected |
| Physical ranges +3 decades, data-based starts | Fit as well as centring on data, but need known units | Rejected |
| Starts and box centred on the data | No new user input; unit-independent; needs facts about each constant | **Chosen** |

**How `Fitting` learns about a constant** (Denis chose the last):

| Approach | Trade-off | Outcome |
|---|---|---|
| Parse names in `Fitting` | No new shared names; brittle, and hides a dependence on the text format | Rejected |
| `parameter_info(m)` returning NamedTuples | One shared name; `kind` checks in place of dispatch | Rejected |
| New types beside the private ones, or a stored name | Two descriptions of one thing, or a name that duplicates the facts | Rejected |
| `Parameter{Kind}`, one parametric struct | Fewest shared names; less conventional | Rejected |
| Abstract `Parameter` with concrete subtypes that hold facts about the step | Textbook Julia, but the fields copy the step | Rejected |
| Abstract `Parameter` with concrete subtypes that point at a step | Textbook Julia; no copied facts; more shared names, justified because parameters are a core concept; every use takes the mechanism | **Chosen** |

## Design

### Stage 1: the parameter types (a refactor)

No behaviour changes in this stage. Every constant keeps its name and order.

A parameter says which step, which direction and which conformation, and nothing else. The
mechanism already holds every other fact, so the type repeats none of them. In
`src/Derivation/names.jl`:

```julia
abstract type Parameter end

struct RateConstant <: Parameter
    step::Int
    reverse::Bool
    inactive::Bool
end

struct EquilibriumConstant <: Parameter
    step::Int
    inactive::Bool
end

struct RegulatorConstant <: Parameter
    ligand::Symbol
    inactive::Bool
end

struct AllostericConstant <: Parameter end
```

These replace `Krapid`, `Kfor`, `Krev` and `Kreg`. `step` indexes the steps of the catalytic
graph. `steps[k].rapid_equilibrium` picks the type, as `_edge_constants` does today.

Facts computed from the mechanism:

| Fact | Source |
|---|---|
| Metabolites in the direction the constant runs | `steps[p.step].taken_up` and `given_off`, swapped when `reverse` |
| Reactant name of a modifier copy | `Metabolite.name` |
| Whether a K reads as a dissociation constant | `step_kinds[p.step]` is `:binding` or `:fused_binding` |
| Name | `name(p, m, kept)`, from `step_sides` and the kept member, as today |
| Conformation tag such as `:A` or `:EqualAI` | `inactive` and the step's or ligand's tag |

`name(p, m[, kept])` stays the chokepoint, for concrete mechanisms only, and derives the tag
that the private types' `state` field held. `_solve` builds `AllostericConstant()` and takes
its name through `name`, so `:L` is no longer a special case. `regulator_constant` keeps
returning a name or `nothing`, computed through `name`, because the law's variables are
keyed by name.

**`_Column` stays.** The joint solve keeps `(edge, slot, inactive)` as its private label for
an unknown, because it needs positions and no parameter types. Names and types meet it in
one place: `_column_parameter` converts a column to a `Parameter`.

New function, in `Derivation`:

```julia
fitted_parameters(m) -> Vector{Parameter}
```

It returns the fitted constants in `FittedParameters` order, for concrete and compiled
mechanisms alike: `_column_parameter` of each independent column, then the regulator
constants, then `AllostericConstant()`. A fitted column is its class's kept member, so its
`step` holds the metabolites its name describes.

`Derivation.API` changes as follows. Adds: `Parameter`, `RateConstant`,
`EquilibriumConstant`, `RegulatorConstant`, `AllostericConstant`, `fitted_parameters` and
`name`. Removes: `RateConstants` and its private helper `_rate_constant_names`, once
`_scale_rate_constants` dispatches on the type. None of these is exported, so the public API
does not change. Stage 1 moves `_scale_rate_constants` to the types, because that use is
what lets `RateConstants` go.

### Stage 2: the search

**References, from the table of the `FittingProblem`:**

- `s_M`: the geometric mean of metabolite `M`'s nonzero concentrations.
- `r`: the geometric mean of `|Rate|`, in both modes. The relative loss ignores a common k
  scale, so `r` changes nothing there: rates scaled by 1e-4 gave identical fits in the laws
  tested.

**Centres, private to `Fitting`, by dispatch:**

```julia
function _log_centre(p::RateConstant, m, s, r)
    e = Mechanisms.catalytic(m).steps[p.step]
    log(r) - _sum_log(s, p.reverse ? e.given_off : e.taken_up)
end
function _log_centre(p::EquilibriumConstant, m, s, r)
    cm = Mechanisms.catalytic(m)
    e = cm.steps[p.step]
    cm.step_kinds[p.step] in (:binding, :fused_binding) ? _sum_log(s, e.taken_up) :
        _sum_log(s, e.given_off) - _sum_log(s, e.taken_up)
end
_log_centre(p::RegulatorConstant, m, s, r) = log(s[p.ligand])
_log_centre(::AllostericConstant, m, s, r) = 0.0
```

`m` is the concrete mechanism, lifted once per fit from `fp.mechanism`. `_sum_log` reads each
`Metabolite.name`.

**Search:** `x0 = centre + 2 * randn()`, clamped to the box, and box `centre ± 15`. Each
restart draws anew, and the box width stays as it is today.

**`fit_rate_equation`:** keeps `lb` and `ub`, which still hold absolute log values. Only
their defaults become data-based. Cross-validation folds get their centres from their
training rows, because each fold builds its own `FittingProblem`.

**Relative-mode scale, done in stage 1:** `_scale_rate_constants` dispatches on the type and
keeps its behaviour. `merge` keeps the order of `params` and leaves entries that are not
fitted constants, such as `Keq` and `E_total`, unchanged; tests pass those.

```julia
_scaled(::Parameter, value, factor) = value
_scaled(::RateConstant, value, factor) = value * factor

_scale_rate_constants(m, params, factor) = merge(params, NamedTuple(
    Derivation.name(p, m) => _scaled(p, params[Derivation.name(p, m)], factor)
    for p in Derivation.fitted_parameters(m)))
```

**At-bound flag:** a fourth result field, `at_bound`, a tuple of fitted-constant names. A
name appears when the best point lies within 0.01 log units of its `lb` or `ub`. The field
reads the point before the relative-mode scale, which can move a reported k outside the box.
The flag means "at bound", never "wrong": a constant no data constrain wanders across the
box.

`at_bound` reaches the search CSVs, `loocv_results.csv` and `best_equation.csv` as an
`at_bound` column, reporting the fit on all rows, never a fold's.

**Docs:** `docs/src/getting_started.md` near the `Keq` paragraph, with a pointer from the
`FittingProblem` docstring:

> We strongly recommend concentrations in M, because `Keq` values usually assume a 1 M
> standard state. Other units work as long as all metabolites, `Keq` and any `reference`
> concentrations use the same one. K and k then carry that unit: a K in mM, a second-order
> k per mM per time. The fit searches around the concentrations in your data, so the unit
> doesn't otherwise change the constants your data determine.

The docs also note that the box spans about 6.5 decades either side of each centre, and
explain `at_bound`. `fit_rate_equation`'s docstring describes the new defaults and field.

**Project docs, each in the commit that changes the code it describes:**

- `.claude/CLAUDE.md` chokepoint section: the new types.
- `docs/src/developer.md`: the API table, the parameter types, and the module map
  (`julia test/update_module_map.jl`).
- `src/Mechanisms/types.jl`: the error text "rendered Kreg names" and the tests that pin it,
  in `test/Mechanisms/test_types.jl` and `test/Mechanisms/test_reaction_graph.jl`.

**Roadmap:** `ROADMAP.md` gains these entries:

- Known issues: a regulator measured only at saturating levels gets an arbitrary K. With an
  activator at 0 and 100× K, restarts that fit perfectly reported 0.015× to 2× the true K.
- Decided against: a global ±ln(1e9) box, physical ranges or a plausibility warning (either
  needs a declared unit), and `lb`/`ub` in `identify_rate_equation`.

## Tests (written first)

**Stage 1:**

- **Same names, same order:** `[name(p, m) for p in fitted_parameters(m)] == parameters(m,
  FittedParameters())` for every mechanism in `MECHANISM_TEST_SPECS`, plain and allosteric.
- **`test/Derivation/test_names.jl`:** about 40 lines move to the new types, with the same
  expected names. The `T(1, :None)` tests become equality and inequality tests of the new
  types' fields. The compiled-mechanism `MethodError` checks stay. Assertions that give a
  plain mechanism an `:I` state, such as `Krapid(bind, :I)`, move to an allosteric fixture
  with the same expected names, because the conformation tag now comes from the mechanism.
- **`RateConstants` tests:** the scaling assertions in `test_rate_eq_derivation.jl` and
  `test_fitting_scale.jl` read the fitted `RateConstant`s. `test_types.jl:3043` asserts that
  `fitted_parameters(em)` throws the same message.
- **Guard:** the chokepoint testset still passes, and fails on a stray renderer.
- **Permitted test edits in this stage:** call forms only: the new constructors, the
  `_edge_constants` arguments and the replacements above, in `test_names.jl`,
  `test_joint_solve.jl`, `test_types.jl`, `test_reaction_graph.jl`,
  `test_rate_eq_derivation.jl`, `test_fitting_scale.jl` and the shared fixtures file. No
  expected name or value changes, so a seeded result that shifts is a refactoring bug.

**Stage 2:**

- **Facts match the law:** for every mechanism in `MECHANISM_TEST_SPECS`, scale each
  metabolite by its own factor `f_M`, `r` by `g`, and each constant by
  `exp(_log_centre(p, s .* f, g * r) - _log_centre(p, s, r))`, and scale `Keq` by the
  product of `f_M` to its stoichiometry. `rate_equation` must scale by exactly `g`. This
  checks the type's facts and `_log_centre` against the derived law. The test also covers
  steps moving two metabolites in one direction, which a search of the specs' file did not
  find there: `E + A + B <--> E(A, B)`, `E + A + B ⇌ E(A, B)`, `F + A + A ⇌ F(A, A)` and
  `E(P, Q) <--> E + P + Q`, each written inline. Modifier copies are covered too.
- **Unit independence:** the same seeded fit in M, mM and µM returns the same determined
  constants after conversion, compared after the relative-mode scale. Throwaway runs agreed
  to about 1e-6 log units; the tolerance sits above that.
- **Issue #97 regression:** K = 0.2 µM with data in M comes back at 0.2 µM, not at 3.06e-7.
- **Both modes:** one relative and one absolute fit recover a steady-state law.
- **Two-metabolite law:** one fit recovers a law with a step moving two metabolites in one
  direction, which no throwaway run covered.
- **Always-zero metabolite:** the fit runs and uses the fallback below.
- **`at_bound`:** a stub optimizer at a known edge names that constant, explicit `lb`/`ub`
  included; a fit inside the box gives `()`; a relative-mode scale does not change it.
- **CSV column:** the search CSVs, `loocv_results.csv` and `best_equation.csv` carry it.
- **Regressions:** the `rate_equation` performance tests stay untouched. Tests that pin
  `keys(res) == (:params, :loss, :retcode)` gain the new field.

## Claude's calls, for Denis to check

- The type names and `fitted_parameters`. `AllostericConstant` holds only `L`.
- `RateConstant.reverse` in place of a slot number, and `inactive` in place of the `state`
  symbol, because the mechanism gives the tag.
- The rule that a binding's K reads as a dissociation constant appears in `name`, in
  `_joint_system` and in `_log_centre`. A third copy costs less than another shared name, and
  the "Facts match the law" test catches a disagreement.
- The new types stay unexported. Other modules reach them qualified.
- The 0.01 tolerance for `at_bound`. In the throwaway reproduction of the issue, it found
  a constant at a wall in 19 of 20 restarts. No run compared a tighter or looser value.
- A metabolite that is always zero uses the geometric mean of all positive concentrations
  in the table, which reached the true fit in 20 of 20 restarts. A table with no positive
  concentration at all uses a centre of 1.
- The flag appears in results and CSVs, with no console warning: a beam search would print
  many, and undetermined constants would trigger them harmlessly.
- `at_bound` goes last in the result, so `fitted, loss, retcode = fit_rate_equation(...)`
  in the docs keeps working.

## Left to the builder

- Carrying `at_bound` through the memoized `_Fit`, which `pmap` returns and which stores
  values in the order of `_fit_order`, so mechanisms sharing one fit report their own names,
  and through L phantoms.
- Handling fits with no finite loss, which have no best point and report `()`.
- The CSV encoding of the names, and any reader of those CSVs that must accept the column.

## Known risks

- Seeded tests that pin an identify result may change in stage 2, because starts change. If
  one breaks, bring it to Denis, and neither loosen it nor reseed it to pass.
- Evidence so far is noiseless synthetic data on laws of 5 to 15 parameters, with one seed
  per setting. Neither James's data, a 20-parameter law nor a step moving two metabolites in
  one direction has been run.
- Whether this PR bumps the version is Denis's call at PR time.
- The specs and plans under `docs/superpowers/` for this work are deleted in the PR that
  implements it.
