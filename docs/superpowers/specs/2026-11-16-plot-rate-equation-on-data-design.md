# Plotting a rate equation on data: design

## Intent

After a fit, Denis wants to see whether the rate equation follows the data, one
measurement group at a time, and to use the same figure in papers and talks. The old
package, DataDrivenEnzymeRateEqs.jl, had `plot_fit_on_data` for this. EnzymeRates has
nothing, and `docs/src/roadmap.md` lists plotting as planned.

This change adds one exported function that draws a mechanism's rate equation, at
given parameters, over a rate table: one panel per `group`.

## Goals

- One call draws every group of a table, with no per-group hand-written labels.
- The figure is good enough for a paper without rework: clean axes, readable units,
  nothing overlapping at the default size.
- Loading EnzymeRates alone never installs or loads Makie.
- The user can change figure size, fonts and axis attributes through Makie's own names.

## Non-goals

- Special handling for log x axes: no tick rule of our own, and no way to show x = 0
  points. Linear is the common case. A log axis works through `axis` when the data
  allow it (see "Log x through `axis`").
- Fit summaries (measured against predicted, error per group), error numbers on
  panels, and overlaying two mechanisms. These were mocked up and declined.
- Choosing groups, ordering panels or overriding the x metabolite by argument. The
  user filters and orders the table before the call.
- Custom colours, markers or legend placement.

## Anti-goals

- The function must not fit anything. It draws the rate equation at exactly the
  parameters given.
- It must not change the user's global Makie theme.
- It must not assume a table has columns beyond those `FittingProblem` needs. The old
  function needed `source` and `X_axis_label` columns.

## Constraints

- A new exported name and its arguments need Denis's approval. The signature below is
  approved as written; any change to it goes back to him.
- `test/test_type_stability.jl` needs a row for every exported function.
- `test/test_boundaries.jl` and the tables and module map in `docs/src/developer.md`
  must agree with `src/`.
- One Julia process at a time on the development machine (7.7 GB, no swap).
- The package supports Julia 1.11 and later, so package extensions are available.

## Approaches

Three ways to bring Makie in were considered.

1. **Package extension (chosen).** Makie is a weak dependency. EnzymeRates declares
   the function; the drawing method lives in `ext/` and loads when the user loads a
   Makie backend. Fitting on a cluster never pays for Makie. The cost is a compat bound
   on Makie to maintain, CairoMakie in the test and docs environments, and a helpful
   error when the function is called without Makie.
2. **Makie as an ordinary dependency.** The simplest code, and what the old package
   did with CairoMakie. Rejected because every user, including cluster jobs that never
   plot, would install and precompile the plotting stack.
3. **No Makie.** Return the prepared series and let the user draw. Rejected because
   the multi-panel layout is the part worth having.

Denis chose 1. The design keeps approach 3's output as an internal layer: plain-Julia
preparation inside the package, and an extension that only draws.

`MakieCore` was checked as a lighter base: Makie has not depended on it since 0.23,
so it does not apply.

## The function

```julia
plot_rate_equation_on_data(m, params::NamedTuple, table; Keq, absolute = false,
    concentration_unit = "M", columns = 5, figure = (;), axis = (;)) -> Makie.Figure
```

```julia
using EnzymeRates, CairoMakie
fig = plot_rate_equation_on_data(m, fit.params, data; Keq = 20_000.0)
save("pkm2.pdf", fig)
```

The first five arguments mean what they mean in
`test_rate_equation(m, params, table; Keq, absolute = false)`:

- `m`: a mechanism, in any form `FittingProblem` accepts.
- `params`: must hold every name of `parameters(m, FittedParameters())`; other names
  are ignored. Values need not come from a fit.
- `table`: needs `group`, `Rate` and one column per metabolite of `m`, and passes the
  same checks as `FittingProblem`. Other columns are ignored.
- `Keq`: required.
- `absolute`: `false` for relative data, `true` for absolute per-enzyme turnover.

The remaining keywords:

- `concentration_unit`: the unit of the table's concentrations. With `"M"`, labels
  use mM, µM or nM. Any other string is shown as given and the numbers are not
  converted.
- `columns`: panels per row.
- `figure`: a `NamedTuple` passed to `Makie.Figure`, for example
  `(; size = (510, 330), fontsize = 7.5)`.
- `axis`: a `NamedTuple` passed to every panel's `Makie.Axis`.

Values in `figure` and `axis` override the built-in defaults. A user can make labels
collide this way; that is accepted.

Called without a Makie backend loaded, the function fails with a message that says to
load CairoMakie or another Makie backend.

## What a panel shows

These rules were settled on mockups drawn from the PKM2 table (591 rows, 28 groups).

- **Panels.** One per distinct `group`, in order of first appearance in the table,
  filled row by row, `columns` per row. The title is the group name with underscores
  shown as spaces.
- **X metabolite.** The metabolite with the most distinct concentrations in the group.
  A tie goes to the earlier one in `metabolites(m)`. On PKM2 this matched all 28
  hand-written labels of the old table.
- **Series.** Each distinct combination of the other metabolites' concentrations is a
  series, in sorted order. A series gets one colour and one marker shape, so colour is
  never the only cue.
- **Points.** The measured `Rate` against the x metabolite.
- **Curve.** One per series: the rate equation over 200 evenly spaced x values from 0
  to the group's largest x, at the series' other concentrations, with `E_total = 1`.
- **Group scale.** With `absolute = false`, each group's curves are multiplied by one
  factor, `exp(mean(log(|Rate| / |predicted|)))` over the group's rows. This is the
  centring the fitting loss does and it changes no rate constant. With
  `absolute = true` there is no factor.
- **Legend.** To the right of the panel. Metabolites held at one nonzero concentration
  across the group are listed first as text, for example `ADP 2 mM`. Metabolites held
  at zero are not mentioned. Below that, a legend headed by the names of the
  metabolites that vary between series, with one entry per series giving their
  concentrations.
- **Y axis.** Linear, including zero. The upper limit covers the data with a small
  margin, and extends to show curves up to 1.6 × the largest data point; a curve above
  that is clipped. Groups of negative rates mirror this rule. The label `Rate` appears
  on the first column only. When the largest absolute rate is below 0.01, tick labels
  are shown as multiples of a power of ten and a small `×10ⁿ` note sits in the panel
  corner.
- **X axis.** Linear from 0, three or four ticks. The label is the metabolite and
  unit, for example `PEP (mM)`.
- **Units.** With `concentration_unit = "M"`, each panel's x axis uses the prefix that
  fits its largest x (mM at 1e-3 and above, µM at 1e-6 and above, nM below), and each
  legend value picks its own prefix. Numbers print with at most four significant
  digits and no trailing zeros.
- **Style.** No grid, no top or right spine, regular-weight titles. Every text size
  derives from the figure's `fontsize` (default 10), so one number rescales the whole
  figure. The default figure size is 320 × 200 units per panel.
- **Theme.** The defaults are applied with `Makie.with_theme` inside the call. A spike
  confirmed that the global theme is unchanged afterwards and that the returned figure
  keeps its look outside the block.

## Log x through `axis`

Denis wants `axis = (; xscale = log10)` to work when Makie allows it, that is, when no
x value is zero. Two rules make it work, and they apply whenever `axis` holds an
`xscale` other than `identity`:

- Each curve runs from the group's smallest x to its largest, not from 0, at 200
  points evenly spaced on the transformed axis, so the curve stays smooth at the low
  end.
- A group with a point at x = 0 raises an error that names the group and says a log
  axis cannot show zero.

Ticks are Makie's own for that scale, labelled in the panel's unit.

## Calls I made, for Denis to check

- **Sign mismatches.** The group scale uses absolute values over all rows, so a curve
  with the wrong sign is drawn with the wrong sign instead of being hidden. A group
  where a predicted rate is zero skips that row in the scale.
- **Mixed-sign groups.** If a group has both positive and negative rates, the y limits
  apply the rule above to each side. No PKM2 group does this.
- **More series than colours.** Colours come from Makie's seven-colour Wong palette
  and eight marker shapes, each cycling. Past seven series a colour repeats with a
  different marker.
- **Many series.** Where two metabolites change together on every row, each row
  becomes its own series, and a 20-row group draws 20 curves and 20 legend entries.
  This is accepted and not capped. The docstring says how series are formed, so the
  user can see why.
- **Groups with nothing to draw a curve from.** The call never throws on a table
  `FittingProblem` accepts. A group whose largest x is 0 (for example one row, or
  replicates with the tied x metabolite at zero) draws its points and no curve, and
  its axis uses the unprefixed unit. A relative-mode group where no row gives a
  finite, nonzero prediction also draws its points and no curve.
- **Zero x on a log axis.** The error names the group, instead of leaving Makie to
  raise its own.
- **Version.** 0.10.3: a new export is not breaking.
- **Code placement.** The preparation is a top-level file in `src/`, beside
  `identify_rate_equation.jl`, with no new submodule and no new public name beyond the
  export. It returns plain data describing each panel. The extension is
  `ext/EnzymeRatesMakieExt.jl`, triggered by `Makie`.

## Tests

- **Preparation, without Makie.** On small tables written in the test: the x
  metabolite including a tie; the series split; fixed and varying metabolites; the
  group scale against a hand calculation in relative mode and its absence in absolute
  mode; curve values against `rate_equation` called directly; unit labels for both
  kinds of `concentration_unit`; the table checks `FittingProblem` enforces; the curve
  range and the zero-x error under a log `xscale`.
- **Drawing, with CairoMakie.** The call returns a `Figure` with one axis per group in
  the expected grid; `figure` and `axis` values reach the figure and axes; the global
  theme is unchanged after the call; a PNG and a PDF save without error.
- **Without Makie.** The error message names the fix.
- Mechanisms are written out with the macros, as elsewhere in the suite.

## Docs and records

- A docstring on the function.
- A section in the fitting tutorial that shows the call and the resulting figure,
  built by Documenter, so `docs/Project.toml` gains CairoMakie.
- `docs/src/developer.md`: the export, and the extension.
- `docs/src/roadmap.md`: the "Plotting" section goes, or shrinks to what remains
  unplanned. `ROADMAP.md` gains or retires entries to match.
- This spec and the plan are deleted in the pull request that implements them.

## Left to the builder

- The shape of the internal panel description, and how the extension reaches it.
- The Makie compat bound. The spike ran on Makie 0.24.17 through CairoMakie 0.15.17.
- How `test/test_type_stability.jl` covers a function whose only method lives in an
  extension and whose drawing dispatches dynamically by nature.
- Whether `test/test_boundaries.jl` needs to know about `ext/`.
- The order of tests so that the no-Makie error is checked before CairoMakie loads.
- Exact margins, gaps and tick counts, matched by eye to the approved mockup.
- Whether the drawing tests fit in the main test process. The full suite needs about
  4 GB free of 7.7 GB under a 2500M heap hint, and CairoMakie adds to that peak.
  Measure it; if it does not fit, run the drawing tests in a separate process or CI
  job.
- What CairoMakie adds to CI time. It took 140 s to precompile once on the development
  machine and was not measured on CI.

## Reference material

The mockup code is throwaway and lives outside the repo. Its logic is described above
in full. The approved linear mockup is at
https://claude.ai/artifact/Hfp7p2Uk12P7QWQo6JQmxz (private to Denis).
