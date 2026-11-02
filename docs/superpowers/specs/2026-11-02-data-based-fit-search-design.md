# Data-based fit search

Fixes issue #97. A short design: a small change to one fitting function, plus one new
`Derivation` function, a result field, docs and roadmap entries.

## Intent

`fit_rate_equation` searches every constant in `log(p) ∈ [-15, 15]` from
`x0 = 2 * randn()`, whatever the data's units. For data in M, this floor (3.06e-7) cuts off
real binding constants, and the fit gives no warning. James's NADPH K sat exactly on it.

The search should start and run near the concentrations and rates in the data, so the
choice of unit stops mattering, and a fit that ends on the box edge should say so.

## Goals

- The fit recovers a constant the data determine, in M, mM or µM alike.
- A K, a first-order k and a second-order k each get their own start and box.
- One rule serves relative and absolute mode.
- Fit results and CSVs name each constant at a box edge.

## Non-goals

- `lb`/`ub` or a box width in `identify_rate_equation`. Fewer choices make the package
  easier to use. The parameter count varies across the beam, so per-constant lists could
  not apply anyway.
- A declared concentration unit or a physical-plausibility warning.
- Reporting that the data only limit a K from one side (see Roadmap).

## Anti-goals

- Any change to `rate_equation` or `loss!`, including their performance.
- Hard limits at physical values. They recreate the bug for real enzymes just outside them.

## Approaches

| Approach | Trade-off | Outcome |
|---|---|---|
| Global box ±ln(1e9) | Two lines, but 1e-9 means 1 nM only in M; K = 0.2 nM in M failed in 20 of 20 restarts | Rejected |
| Physical ranges per class (K 1 nM–100 mM, k 0.1–1e9) | Walls are meaningful, but need known units; K = 0.2 nM, K = 300 mM and absolute k = 0.02/s failed in 20 of 20 restarts | Rejected |
| Physical ranges +3 decades, data-based starts | Fit as well as centring on data, but need known units | Rejected |
| Starts and box centred on the data | No new user input; unit-independent; needs each constant's concentration order | **Chosen by Denis** |

## Design

**References, from the table of the `FittingProblem`:**

- `s_M`: the geometric mean of metabolite `M`'s nonzero concentrations.
- `r`: the geometric mean of `|Rate|`, in both modes. The relative loss ignores a common
  k scale, so `r` changes nothing there: rates scaled by 1e-4 gave identical fits in the
  laws tested.

**Centre of each fitted constant:** `log(r) * t + Σ_M order_M * log(s_M)`, where `t` is 1
for a rate constant and 0 otherwise, and `order_M` is the constant's concentration order
in `M`.

| Constant | Centre | `order_M` |
|---|---|---|
| Dissociation K (`:binding`, `:fused_binding`), regulator K | `s_M` | +1 |
| First-order k | `r` | 0 |
| Second-order k | `r / s_M` | −1 |
| Isomerization K, `L` | 1 | 0 |
| Other equilibrium K, stored direction, such as `E + A + B ⇌ E(A, B)` | product of `s` | given off +1, taken up −1 |
| Other rate constant | `r` over product of `s` | −1 per metabolite taken up |

**Search:** `x0 = centre + 2 * randn()`, clamped to the box, and box `centre ± 15`. Each
restart draws anew, and the box width stays as it is today.

**`fit_rate_equation`:** keeps `lb` and `ub`, which still hold absolute log values. Only
their defaults become data-based. Cross-validation folds get their centres from their
training rows, because each fold builds its own `FittingProblem`.

**Concentration orders:** one new `Derivation` function gives each fitted constant's orders
in `FittedParameters` order, read from `Step.taken_up` and `Step.given_off` of the joint
solve's columns, plus the regulator constants and `L`. Its name carries no underscore and
joins `Derivation.API`. `docs/src/developer.md` and the module map update in the same
commit (`julia test/update_module_map.jl`).

**At-bound flag:** a fourth result field, `at_bound`, a tuple of fitted-constant names. A
name appears when the best point lies within 0.01 log units of its `lb` or `ub`. The field
reads the point before the relative-mode scale, which can move a reported k outside the
box. Search CSVs, `loocv_results.csv` and `best_equation.csv` gain an `at_bound` column.
The flag means "at bound", never "wrong": a constant no data constrain wanders across the
box.

**Docs:** `docs/src/getting_started.md` near the `Keq` paragraph, with a pointer from the
`FittingProblem` docstring:

> We strongly recommend concentrations in M, because `Keq` values usually assume a 1 M
> standard state. Other units work as long as all metabolites, `Keq` and any `reference`
> concentrations use the same one. K and k then carry that unit: a K in mM, a second-order
> k per mM per time. The fit searches around the concentrations in your data, so the unit
> doesn't otherwise change the constants your data determine.

The docs also note that the box spans about 6.5 decades either side of each centre, and
explain `at_bound`.

**Roadmap:** `ROADMAP.md` gains these entries:

- Known issues: a regulator measured only at saturating levels gets an arbitrary K. With an
  activator at 0 and 100× K, restarts that fit perfectly reported 0.015× to 2× the true K.
- Decided against: a global ±ln(1e9) box, physical ranges or a plausibility warning (either
  needs a declared unit), and `lb`/`ub` in `identify_rate_equation`.

## Tests (written first)

- **Orders match the law:** for every mechanism in `MECHANISM_TEST_SPECS`, scale each
  metabolite by its own factor, each constant by its orders and `Keq` to match, and assert
  that `rate_equation` does not change. This checks the orders against the derived law.
  `MECHANISM_TEST_SPECS` holds a Theorell–Chance step, but a search of its file found no
  step moving two metabolites in one direction, so the test also covers
  `E + A + B <--> E(A, B)`, `E + A + B ⇌ E(A, B)`, `F + A + A ⇌ F(A, A)` and
  `E(P, Q) <--> E + P + Q`, each written inline.
- **Unit independence:** the same seeded fit in M, mM and µM returns the same determined
  constants after conversion, compared after the relative-mode scale. Throwaway runs
  agreed to about 1e-6 log units; the tolerance sits above that.
- **Issue #97 regression:** K = 0.2 µM with data in M comes back at 0.2 µM and not at
  3.06e-7.
- **Both modes:** one relative and one absolute fit recover a steady-state law.
- **Always-zero metabolite:** the fit runs and uses the fallback below.
- **`at_bound`:** a stub optimizer at a known edge names that constant, explicit `lb`/`ub`
  included; a fit inside the box gives `()`; a relative-mode scale does not change it.
- **CSV column:** the search CSVs, `loocv_results.csv` and `best_equation.csv` carry it.
- **Regressions:** the `rate_equation` performance tests stay untouched. Tests that pin
  `keys(res) == (:params, :loss, :retcode)` gain the new field.

## Claude's calls, for Denis to check

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

- Carrying `at_bound` through the memoized `_Fit`, which stores values in the order of
  `_fit_order`, so mechanisms sharing one fit report their own names.
- Handling fits with no finite loss, which have no best point and report `()`.
- The new function's name, and the CSV encoding of the names.

## Known risks

- Seeded tests that pin an identify result may change, because starts change. If one
  breaks, bring it to Denis, and neither loosen it nor reseed it to pass.
- Evidence so far is noiseless synthetic data on laws of 5 to 15 parameters, with one seed
  per setting. Neither James's data, a 20-parameter law nor a step moving two metabolites
  in one direction has been run.
- Whether this PR bumps the version is Denis's call at PR time.
