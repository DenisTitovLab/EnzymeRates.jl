# plot_rate_equation_on_data Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the exported `plot_rate_equation_on_data`, which draws a mechanism's rate
equation at given parameters over a rate table, one panel per `group`.

**Architecture:** A top-level file in `src/` turns the mechanism, parameters and table
into plain panel descriptions, with every number already in display units. A package
extension triggered by `Makie` draws those descriptions and nothing else. The exported
function lives in `src/`, calls the preparation, then hands the panels to the drawing
function whose only method is in the extension.

**Tech Stack:** Julia 1.11+, package extensions, Makie 0.24 (CairoMakie 0.15 in tests
and docs).

**Spec:** `docs/superpowers/specs/2026-11-16-plot-rate-equation-on-data-design.md`.
Read it first; this plan does not repeat its panel rules.

## Global Constraints

- Read `.claude/CLAUDE.md` before starting. TDD, 92-character lines, two `ABOUTME:`
  lines at the top of every new code file, `_testhelper_` prefix on test helpers,
  mechanisms in tests written out with the macros.
- The approved signature, not to be changed without Denis:
  `plot_rate_equation_on_data(m, params::NamedTuple, table; Keq, absolute = false,
  concentration_unit = "M", columns = 5, figure = (;), axis = (;)) -> Makie.Figure`
- No new exported or `public` name other than `plot_rate_equation_on_data`.
- Loading EnzymeRates alone must not load Makie: `Makie` goes in `[weakdeps]` only.
- Run one Julia process at a time. Focused test run:
  `julia --project -e 'using TestEnv; TestEnv.activate(); using Test, EnzymeRates, LinearAlgebra, Random; include("test/mechanism_definitions_for_test_enzyme_derivation.jl"); include("test/<file>.jl")'`
- Test output must be pristine: no warnings printed by the new tests.
- The mockup code is gone; the spec's "What a panel shows" is the full description.
  The approved look is at https://claude.ai/artifact/Hfp7p2Uk12P7QWQo6JQmxz.

## Review Focus

- A table whose `group` values are not strings (integers, symbols): titles and the
  zero-x error must still work. Test in Task 1.
- A group with one row, or with every x at zero: points, no curve, no throw. Task 1.
- A relative-mode group whose predictions are all zero: points, no curve. Task 1.
- Fewer groups than `columns`: the figure must not reserve empty columns. Task 2.
- A call made before any Makie backend is loaded: an error that names the fix, raised
  before any preparation work. Task 2.

## File Structure

- Create `src/plot_rate_equation_on_data.jl`: panel types, the preparation, the
  exported function, the drawing function's declaration.
- Create `ext/EnzymeRatesMakieExt.jl`: the drawing method.
- Create `test/test_plot_rate_equation_on_data.jl`: preparation tests and the
  no-backend error, without Makie.
- Create `test/test_plot_rate_equation_on_data_makie.jl`: drawing tests with CairoMakie.
- Modify `src/EnzymeRates.jl`, `Project.toml`, `test/runtests.jl`,
  `test/test_aqua_jet.jl`, `test/test_type_stability.jl`, `docs/src/developer.md`,
  `docs/src/api.md`, `docs/src/fitting/tutorial.md`, `docs/Project.toml`,
  `docs/src/roadmap.md`, `ROADMAP.md`.

---

### Task 1: Panel preparation

**Files:**
- Create: `src/plot_rate_equation_on_data.jl`
- Modify: `src/EnzymeRates.jl` (add `include("plot_rate_equation_on_data.jl")` after
  `identify_rate_equation.jl`)
- Create: `test/test_plot_rate_equation_on_data.jl`
- Modify: `test/runtests.jl` (include the test file after
  `test_identify_rate_equation.jl`)

**Interfaces:**
- Consumes: `Fitting.FittingProblem(m, table; Keq, absolute)` for validation and the
  compiled mechanism (`fp.mechanism`, `fp.rate_data.data`, `fp.rate_data.Keq`);
  `Derivation.rate_equation`, `Derivation.parameters(m, Derivation.FittedParameters())`,
  `Mechanisms.metabolites`. Do not use `fp.group_point_indexes`: it loses group names
  and order.
- Produces, all in module `EnzymeRates`, none exported or public:

```julia
struct RateEquationSeries
    label::String                 # varying concentrations, e.g. "2 mM, 0"; "" if none vary
    x::Vector{Float64}            # measured x, in display units
    rate::Vector{Float64}         # measured rate, divided by 10^rate_power
    curve_x::Vector{Float64}      # empty when the group has no curve
    curve_rate::Vector{Float64}   # scaled prediction, divided by 10^rate_power
end

struct RateEquationPanel
    title::String                         # group name, "_" shown as " "
    x_label::String                       # "PEP (mM)"
    rate_power::Int                       # 0 unless the largest |Rate| < 0.01
    fixed_concentrations::Vector{String}  # ["ADP 2 mM"]
    legend_title::String                  # "F16BP, Phenylalanine"; "" if none vary
    series::Vector{RateEquationSeries}
    rate_limits::Tuple{Float64, Float64}  # in the same units as `rate`
end

rate_equation_panels(m, params::NamedTuple, table; Keq, absolute::Bool,
    concentration_unit::AbstractString, xscale) -> Vector{RateEquationPanel}
```

Decisions the spec leaves open, fixed here:

- `params` is checked as `test_rate_equation` checks it, with the same error text.
- Display units: with `concentration_unit == "M"` the panel's x values are divided by
  1e-3, 1e-6 or 1e-9 by the spec's thresholds on the group's largest x, and by 1 with
  label unit `M` when that largest x is 0. With any other unit nothing is divided.
- Concentration text: `"0"` for zero; otherwise the number at four significant digits
  with no trailing zeros or trailing `.0`, a space, and the unit. Under `"M"` each
  value picks its own prefix by the same thresholds.
- `rate_power = floor(Int, log10(largest |Rate|))` when that largest is below 0.01,
  else 0.
- Curve x: with `xscale === identity`, 200 points evenly spaced from 0 to the largest
  x. Otherwise 200 points geometrically spaced from the smallest x to the largest,
  and a group with any x equal to 0 throws an `ArgumentError` whose message contains
  the group name, the x metabolite and "cannot show zero".
- Group scale: `exp(mean(log(|Rate| / |predicted|)))` over the rows whose prediction
  is finite and nonzero; 1 in absolute mode. No such row in relative mode, or a
  largest x of 0, means every series of the group has empty `curve_x` and
  `curve_rate`.
- Rate limits, with `hi`/`lo` the largest/smallest measured rate and `chi`/`clo` the
  largest/smallest finite curve value (use `hi`/`lo` when there is no curve):
  `top = hi > 0 ? max(1.12hi, min(1.6hi, 1.05chi)) : 0.0`,
  `bottom = lo < 0 ? min(1.12lo, max(1.6lo, 1.05clo)) : 0.0`.
  All-positive groups get `(-0.05top, top)`, all-negative `(bottom, -0.05bottom)`,
  mixed `(bottom, top)`.
- Series order: sorted by the tuple of varying concentrations. Metabolites listed in
  `legend_title`, in series labels and in `fixed_concentrations` follow
  `metabolites(m)` order.

- [ ] **Step 1: Write the failing tests**

In `test/test_plot_rate_equation_on_data.jl`, write a random-order bi-bi mechanism over
`A, B, P, Q` with `@enzyme_mechanism` (copy the step style of an entry in
`test/mechanism_definitions_for_test_enzyme_derivation.jl`). Pick fixed parameter values
with `_testhelper_random_independent_params` or by hand, and build tables whose `Rate`
column is `factor * rate_equation(...)`, so expected values come from `rate_equation`
called directly. Testsets and their assertions:

- `"x metabolite and series"`: group `"exp_1"` with A at `[1e-4, 2e-4, 5e-4, 1e-3]`
  for each B in `(1e-3, 2e-3)`, `Q = 5e-5`, `P = 0`, `Rate = 3 * predicted`. Assert
  one panel; `title == "exp 1"`; `x_label == "A (mM)"`; `legend_title == "B"`;
  `fixed_concentrations == ["Q 50 µM"]`; two series labelled `"1 mM"` and `"2 mM"`;
  first series `x ≈ [0.1, 0.2, 0.5, 1.0]`; `rate_power == 0`.
- `"tie goes to metabolites order"`: a group where A and B each take two distinct
  values; assert the x label names `first(metabolites(m))` among the tied ones.
- `"group scale"`: for the table above in relative mode, `curve_rate` at every
  `curve_x` equals `3 * rate_equation` at that x (converted back from display units)
  with `rtol = 1e-10`, `length(curve_x) == 200`, `curve_x[1] == 0`,
  `curve_x[end] ≈ 1.0`. With `absolute = true` the same check holds with factor 1.
- `"panels follow first appearance"`: groups appearing in order `"b"`, `"a"`, `"b"`
  give titles `["b", "a"]`.
- `"small rates"`: rates scaled so the largest is 3e-5 give `rate_power == -5` and
  `maximum(series.rate) ≈ 3.0`.
- `"rate limits"`: an all-positive group whose curve stays below the data gives
  `rate_limits == (-0.05 * 1.12hi, 1.12hi)`; a table with `Rate` at one tenth of the
  prediction scale in absolute mode gives a top of `1.6hi`; an all-negative group
  (products only, substrates zero) gives a negative bottom and a top of
  `-0.05 * bottom`.
- `"concentration units"`: with `concentration_unit = "mM"` and A values
  `[0.1, 0.2]`, `x_label == "A (mM)"` and `x == [0.1, 0.2]`; a fixed value of 2.0
  reads `"B 2 mM"`. Under `"M"`: 1.25e-3 reads `"1.25 mM"`, 1.08e-3 `"1.08 mM"`,
  5e-9 `"5 nM"`.
- `"no curve"`: a one-row group, and a group with every candidate x at zero, each
  give a panel with points and empty `curve_x`; the all-zero one has x label unit
  `M`. A relative-mode group at `Keq`-equilibrium concentrations, where every
  prediction is zero, gives empty curves and does not throw.
- `"non-string groups"`: `group = [1, 1, 2, 2]` gives titles `["1", "2"]`.
- `"log x"`: `xscale = log10` gives `curve_x[1] ≈` the smallest x and
  `curve_x[2] / curve_x[1] ≈ curve_x[3] / curve_x[2]`; a group with an x of zero
  throws `ArgumentError` matching the group name and `"cannot show zero"`.
- `"input checks"`: a table without a `group` column throws as `FittingProblem` does;
  `params` missing a fitted name throws with `"params is missing"`.

- [ ] **Step 2: Run the focused test; expect failure on undefined `rate_equation_panels`**

- [ ] **Step 3: Implement the two structs and `rate_equation_panels`**

Build one `FittingProblem` for validation and the compiled mechanism, then group rows
by `data.group` in first-appearance order. Keep helpers private (leading `_`) and few.

- [ ] **Step 4: Run the focused test; expect all pass with no printed warnings**

- [ ] **Step 5: Commit** (`Prepare the panels of a rate equation plot`)

---

### Task 2: The exported function and the Makie extension

**Files:**
- Modify: `src/plot_rate_equation_on_data.jl`, `src/EnzymeRates.jl` (export; the
  ABOUTME count becomes 20), `Project.toml`
- Create: `ext/EnzymeRatesMakieExt.jl`, `test/test_plot_rate_equation_on_data_makie.jl`
- Modify: `test/test_plot_rate_equation_on_data.jl`, `test/runtests.jl`,
  `test/test_aqua_jet.jl`, `test/test_type_stability.jl`, `docs/src/developer.md`,
  `docs/src/api.md`

**Interfaces:**
- Consumes: Task 1's `rate_equation_panels`, `RateEquationPanel`, `RateEquationSeries`.
- Produces:

```julia
# src/plot_rate_equation_on_data.jl
function draw_rate_equation_panels end   # only method is in the extension
plot_rate_equation_on_data(m, params::NamedTuple, table; Keq, absolute = false,
    concentration_unit = "M", columns = 5, figure = (;), axis = (;))

# ext/EnzymeRatesMakieExt.jl
EnzymeRates.draw_rate_equation_panels(panels::Vector{EnzymeRates.RateEquationPanel};
    columns::Integer, figure::NamedTuple, axis::NamedTuple) -> Makie.Figure
```

Decisions fixed here:

- `plot_rate_equation_on_data` first checks
  `Base.get_extension(EnzymeRates, :EnzymeRatesMakieExt) === nothing` and, if so,
  throws an `ErrorException` with exactly: ``plot_rate_equation_on_data needs a Makie
  backend: run `using CairoMakie` (or another Makie backend) first``. It then calls
  `rate_equation_panels` with `xscale = get(axis, :xscale, identity)` and passes the
  result to `draw_rate_equation_panels`.
- `Project.toml`: `[weakdeps] Makie`, `[extensions] EnzymeRatesMakieExt = "Makie"`,
  `[compat] Makie = "0.24"`, `CairoMakie = "0.15"`; `CairoMakie` in `[extras]` and the
  test target; `version = "0.10.3"`.
- `test/test_aqua_jet.jl`: add `:CairoMakie` to the `stale_deps` ignore list.
- Drawing defaults, all applied inside `Makie.with_theme` and overridden by the
  user's `figure` and `axis`: `fs = get(figure, :fontsize, 10)`; columns used
  `= min(columns, length(panels))`; figure size `(320 * columns_used, 200 * rows)`;
  figure padding 10; column gap 12, row gap 10; no grid; no top or right spine;
  titles in the regular font at `fs`; tick labels and legend text at `fs - 1`; three
  or four ticks per axis (`WilkinsonTicks(3; k_min = 3, k_max = 4)`); y label
  `"Rate"` on the first column only; lines 1.6 wide; markers of size 7 with a
  0.5-wide white stroke; colours `Makie.wong_colors()` and markers `[:circle, :rect,
  :diamond, :utriangle, :dtriangle, :pentagon, :hexagon, :star5]`, each cycling by
  series index.
- Each panel is its own `GridLayout`: the `Axis` at `[1, 1]`, and at `[1, 2]` a
  top-aligned column holding the fixed concentrations as one left-justified `Label`
  (one per line, grey) and, when `legend_title != ""`, a frameless `Legend` titled
  with it. `ylims!` from `rate_limits`. When `rate_power != 0`, a small grey
  `×10ⁿ` text in the corner on the data's side of zero (top left for positive
  rates, bottom left for negative).
- With a non-identity `xscale`, leave the x ticks to Makie.

- [ ] **Step 1: Write the failing tests**

In `test/test_plot_rate_equation_on_data.jl`, a testset `"needs a Makie backend"`
placed where no backend is loaded yet: assert
`Base.get_extension(EnzymeRates, :EnzymeRatesMakieExt) === nothing`, then
`@test_throws "needs a Makie backend" plot_rate_equation_on_data(m, params, table; Keq)`.

In `test/test_plot_rate_equation_on_data_makie.jl` (`using CairoMakie`), on a
three-group table from the Task 1 mechanism:

- `"returns a figure with one axis per group"`: result `isa Makie.Figure`; exactly
  three `Makie.Axis` among its content; their titles equal the group titles in order.
- `"grid"`: with `columns = 2` the figure size is `(640, 400)`; with the default
  `columns` and three groups it is `(960, 200)`.
- `"figure and axis keywords"`: `figure = (; size = (500, 300), fontsize = 8)` gives
  that size; `axis = (; xgridvisible = true)` sets it on every axis.
- `"log x"`: `axis = (; xscale = log10)` on a table without zero x returns a figure
  whose axes have `xscale[] === log10`; with a zero x it throws `ArgumentError`.
- `"global theme unchanged"`: `Makie.current_default_theme().fontsize[]` is equal
  before and after a call made with `figure = (; fontsize = 7)`.
- `"saves"`: `save` to a `.png` and a `.pdf` in `mktempdir()` gives non-empty files.
- `"legend and labels"`: a group with a varying and a fixed metabolite has one
  `Makie.Legend` in the figure; a group with no varying metabolite has none.

Include the Makie test file in `test/runtests.jl` directly after the preparation test
file.

- [ ] **Step 2: Run both focused test files; expect failures on the undefined function**

- [ ] **Step 3: Implement the function, its docstring, the export and the extension**

The docstring covers arguments, the result, units, how series are formed, that x is
linear unless `axis` sets `xscale`, that a log axis rejects zero x, and the errors.

- [ ] **Step 4: Run both focused test files; expect pass, pristine output**

- [ ] **Step 5: Update the records that tests guard**

- `docs/src/developer.md`: add a `plot_rate_equation_on_data` row to "The top level"
  interface table (exported ✓), add `src/plot_rate_equation_on_data.jl` and the
  extension to the top-level paragraph, then run `julia test/update_module_map.jl`.
- `docs/src/api.md`: add the function to the Functions `@docs` block.
- `test/test_type_stability.jl`: add a row for `plot_rate_equation_on_data`. It must
  give the same result whether or not CairoMakie is loaded in the process. Drawing a
  figure dispatches dynamically by nature and never runs per data row, so a new
  named exemption set with that reason is acceptable; keep its function list as
  short as the reported sites allow. Run it with `julia +1.13`.

- [ ] **Step 6: Run the focused tests for `test/test_boundaries.jl`,
  `test/test_aqua_jet.jl` and, on `julia +1.13`, `test/test_type_stability.jl`; all pass**

- [ ] **Step 7: Commit** (`Add plot_rate_equation_on_data, drawn by a Makie extension`)

---

### Task 3: Documentation and roadmap

**Files:**
- Modify: `docs/src/fitting/tutorial.md`, `docs/Project.toml`, `docs/src/roadmap.md`,
  `ROADMAP.md`

- [ ] **Step 1: Add a "Plotting the fit on the data" section to the tutorial**

Place it after "Running the fit". Reuse the page's existing `fitting` example block
names so it plots that page's own mechanism, fit and table: `using CairoMakie`, the
call, and the figure as the block's result. Two or three sentences say what a panel
shows and how the x metabolite and series are chosen, and one names `figure`, `axis`
and `concentration_unit`. Add `CairoMakie` to `docs/Project.toml` with compat `"0.15"`.

- [ ] **Step 2: Update the roadmaps**

Remove the "Plotting" section from `docs/src/roadmap.md`. In `ROADMAP.md` under
"Ideas", add one entry listing what was considered and left out: overlaying two
mechanisms, a fit summary (measured against predicted, error per group), a tick rule
for log axes, and a keyword to override the inferred x metabolite.

- [ ] **Step 3: Build the docs**

Run: `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'`
Expected: the build finishes with no warnings from the new section, and
`docs/build/fitting/tutorial/index.html` references a generated image. `deploydocs`
does nothing off CI.

- [ ] **Step 4: Commit** (`Document plot_rate_equation_on_data`)

---

### Task 4: Whole-suite verification and cleanup

- [ ] **Step 1: Run the full suite in the background and measure its peak memory**

Run: `/usr/bin/time -v timeout 1200 julia --project -e 'using Pkg; Pkg.test(julia_args=["--heap-size-hint=2500M"])'`
Expected: all tests pass. Record "Maximum resident set size" of the largest process.
If the suite is killed or fails for memory, stop and report: the fallback is to run
the Makie test file in its own process or CI job, which changes
`test/runtests.jl` and `.github/workflows/CI.yml`.

- [ ] **Step 2: Run the full suite on `julia +1.13`** the same way; expect pass.

- [ ] **Step 3: Delete the spec and this plan**

`git rm docs/superpowers/specs/2026-11-16-plot-rate-equation-on-data-design.md docs/superpowers/plans/2026-11-16-plot-rate-equation-on-data.md`

- [ ] **Step 4: Commit** (`Remove the plot design and plan now that they are implemented`)
