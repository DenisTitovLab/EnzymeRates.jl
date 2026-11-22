# Normalized vs absolute rate

Enzyme kinetic data is often reported on an arbitrary rate scale. The enzyme
concentration in the assay may be unknown, and figures from different
preparations may reflect different states of purity, so the absolute value of a
measured rate carries no reliable meaning. What survives is the *shape* of the
response: changing the amount of enzyme scales every rate by the same factor and
leaves that shape intact.

EnzymeRates.jl handles both cases through two `FittingProblem` keywords, `absolute` and
`reference`. With the default `absolute = false` the data is relative: the fit ignores an
arbitrary factor on each group's rates, and `reference`, one or more turnovers
measured at stated concentrations, sets the scale of the returned rate constants. Without
a `reference` the scale is set so that the largest rate per unit enzyme the mechanism
predicts over the data rows is 1.
With `absolute = true` the data is per-enzyme turnover: the fit scores absolute
magnitudes, returns the constants as fitted and takes no `reference`.
`IdentifyRateEquationProblem` takes the same keywords.

## The two modes

### Relative mode (the default)

- **Loss**: each group's log-ratios are **mean-centered before squaring**, which removes the
  arbitrary per-group scale, such as the group's enzyme amount. The loss is invariant to
  rescaling all rates in a group by the same factor.
- **Output scale**: the loss cannot see a common factor on every steady-state rate
  constant, so the fitter fixes it after the fit: the returned constants make the
  mechanism reproduce the **reference** turnovers or, without a reference, put the
  fastest measured condition at a rate of 1 (below). The factor multiplies only the rate
  constants — the lowercase-`k` parameters every steady-state step contributes, in both
  directions (`k_<from>_to_<to>`) — and leaves the binding constants `K`, `Keq`,
  the allosteric `L`, and regulatory K values unchanged.

### The reference

`reference` is optional in relative mode: a vector of anchors `concentrations =>
turnover`.

- Each anchor's concentrations are a `NamedTuple` in the units of your data, and a
  metabolite you leave out is zero.
- Its turnover is per unit enzyme, per active site: a kcat, or a specific activity in
  turnover units.
- When you know no turnover in absolute units, pass a row of your own data, its
  concentrations `=>` its `Rate`. The constants then predict rates on that row's scale.
- Turnover is per active site, so divide a kcat quoted per oligomer by its number of
  active sites first.

```julia
fp = FittingProblem(m, data, 2.0; reference = [(S = 10.0, P = 0.0) => 42.0])
```

The fit multiplies every steady-state rate constant by the geometric mean of measured over
predicted turnover at the anchors; several anchors combine in log space. A fit whose
predicted turnover at an anchor is zero, non-finite or of the other sign has no such factor,
for example a law that needs a substrate the anchor leaves at zero. Such a fit fails:
`fit_rate_equation` returns its constants unscaled, with `loss = Inf` and
`retcode = :ReferenceFailed`, and identify drops the mechanism from the search. Name every
substrate and every essential activator in an anchor so each candidate law runs there: a law
that needs an activator predicts zero turnover where the anchor leaves it at zero, and its
fit fails. In identify, a mechanism that leaves out an optional regulator ignores that
regulator's concentration in the anchors.

Without a `reference` the fit multiplies every steady-state rate constant by one factor
so that the largest absolute rate per unit enzyme the mechanism predicts over your data's
concentration rows is 1. The default-scaled rate constants are therefore in units of
the fastest measured condition, not kcats: to read them as turnovers, pass a `reference`
with measured turnovers. A fit whose largest predicted rate is zero or non-finite fails
like a failed anchor: constants unscaled, `loss = Inf`, `retcode = :ReferenceFailed`.

Relative mode is the right choice when your data reports relative rates
(arbitrary units, percent activity, or values normalized per some enzyme-concentration
estimate that you do not trust).

### Absolute mode

Pass `absolute = true`. In absolute mode:

- **Loss**: log-ratios are squared **without centering**. The absolute
  magnitude of the rate is scored, so the data must report per-enzyme turnover
  in consistent units.
- **Output**: the raw fitted parameters, never rescaled. The data fixes the absolute
  scale.

Absolute mode is the right choice when your data is already on a per-enzyme
turnover basis.

!!! note "The fitter has no enzyme-concentration input"
    The rate equation gives the rate per unit enzyme, so the loss compares your
    measured `Rate` against the per-enzyme turnover. The data table carries no
    enzyme-concentration column. Relative mode makes this moot: per-group centering
    cancels any constant factor, the enzyme amount included. Absolute mode does not, so
    the measured `Rate` must already be per-enzyme — divide rates recorded at a
    known enzyme amount by that amount before fitting.

## Scoring a fit on other data

[`EnzymeRates.test_rate_equation`](@ref) scores a fit's constants on data the fit did not
see, as cross-validation does for each held-out group. In relative mode the score does not
depend on the scale of the constants.

```@example testfit
using EnzymeRates

uni_uni = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> ES
        ES <--> E + P
    end
end

params = (k_ES_to_E_P = 3.0, k_ES_to_E_S = 1.0, k_E_S_to_ES = 4.0)
heldout = (group = ["G1", "G1", "G2"], Rate = [0.5, 0.9, 1.1],
           S = [1.0, 5.0, 3.0], P = [0.0, 0.0, 0.0])

EnzymeRates.test_rate_equation(uni_uni, params, heldout, 2.0)
```
