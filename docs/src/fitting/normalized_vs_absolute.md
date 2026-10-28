# Normalized vs absolute rate

Enzyme kinetic data is often reported on an arbitrary rate scale. The enzyme
concentration in the assay may be unknown, and figures from different
preparations may reflect different states of purity, so the absolute value of a
measured rate carries no reliable meaning. What survives is the *shape* of the
response: changing the amount of enzyme scales every rate by the same factor and
leaves that shape intact.

EnzymeRates.jl handles both cases through two `FittingProblem` keywords, `absolute` and
`reference`. With the default `absolute = false` the data is relative: the fit ignores an
arbitrary factor on each group's rates, and you pass `reference`, one or more turnovers
measured at stated concentrations, which sets the scale of the returned rate constants.
With `absolute = true` the data is per-enzyme turnover: the fit scores absolute
magnitudes, returns the constants as fitted and takes no `reference`.
`IdentifyRateEquationProblem` takes the same keywords.

## The two modes

### Relative mode (the default)

- **Loss**: each group's log-ratios are **mean-centered before squaring**,
  which removes the arbitrary per-group `E_total` scale. The loss is invariant
  to rescaling all rates in a group by the same factor.
- **Output scale**: the loss cannot see a common factor on every steady-state rate
  constant, so the fitter fixes it after the fit: the returned constants make the
  mechanism reproduce the **reference** turnovers, as `rescale_parameter_values` does
  with measured turnovers (below).

### The reference

`reference` is required in relative mode: a vector of anchors `concentrations =>
turnover`.

- Each anchor's concentrations are a `NamedTuple` in the units of your data, and a
  metabolite you leave out is zero.
- Its turnover is per unit `E_total`, the active-site concentration: a kcat, or a
  specific activity in turnover units.
- When you know no turnover in absolute units, pass a row of your own data, its
  concentrations `=>` its `Rate`. The constants then predict rates on that row's scale.

```julia
fp = FittingProblem(m, data; Keq = 2.0,
                    reference = [(S = 10.0, P = 0.0) => 42.0])
```

The fit multiplies every steady-state rate constant by the geometric mean of measured over
predicted turnover at the anchors. A fit whose predicted turnover at an anchor is zero,
non-finite or of the other sign has no such factor, for example a law that needs a
substrate the anchor leaves at zero. Such a fit fails: `fit_rate_equation` returns its
constants unscaled, with `loss = Inf` and `retcode = :ReferenceFailed`, and identify drops
the mechanism from the search. Name every substrate and every essential activator in an
anchor so each candidate law runs there: a law that needs an activator predicts zero
turnover where the anchor leaves it at zero, and its fit fails. In identify, a mechanism
that leaves out an optional regulator ignores that regulator's concentration in the
anchors.

Relative mode is the right choice when your data reports relative rates
(arbitrary units, percent activity, or values normalized per some `E_total`
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

!!! note "The fitter has no `E_total` input"
    Every prediction is evaluated at `E_total = 1`, so the loss compares your
    measured `Rate` against the per-enzyme turnover — in effect it fits
    `Rate / E_total`. The data table carries no `E_total` column, and the fitter
    ignores one if present. Relative mode makes this moot: per-group centering
    cancels any constant factor, `E_total` included. Absolute mode does not, so
    the measured `Rate` must already be per-enzyme — divide rates recorded at a
    known enzyme amount by that amount before fitting.

## Rescaling after the fit

`rescale_parameter_values` re-expresses fitted parameters on another scale. It rescales
only the rate constants — the lowercase-`k` parameters every steady-state step
contributes, in both directions (`k_<from>_to_<to>`) — and leaves the binding constants
`K`, `Keq`, `E_total`, the allosteric `L`, and regulatory K values unchanged.

The rate constants are the only parameters that carry time in their units, and
the measured rate is the only thing that supplies a time scale, so it can pin
down only what lives in them. The binding constants, `Keq`, and `E_total` are
concentrations or concentration ratios, and `L` is a pure number — none carry a
time dimension, so the rate constrains them only as far as the response shape
already does. Multiplying every rate constant by a common factor multiplies every
rate by that factor, and the K's stay put: a rapid-equilibrium step stores its binding
constant directly, so the rescaling never touches it, while a steady-state binding
constant is the ratio of its two rate constants, in which the common factor cancels.

The parameters you pass must include `Keq`, which the dependent rate constants need; a
fit's `params` holds only the fitted constants, so merge it in.

Pass `kcats`, a vector of `concentrations => turnover` pairs in the form `reference`
takes.

- Each is a measured turnover — a kcat or a specific activity in turnover units — at
  stated concentrations, in the units of your data. Metabolites you leave out are zero.
- Turnover is per unit `E_total`, the active-site concentration, so divide a kcat quoted
  per oligomer by its number of active sites first.
- One anchor makes the predicted turnover equal the measured one there; several combine
  by the geometric mean of measured over predicted turnover.
- The second return value, the spread, is the standard deviation of the anchors'
  `log(measured/predicted)`: 0 for one anchor, and large when the anchors disagree with
  the model or with each other.

A relative fit comes back scaled to its `reference` this way; rescale it to put it on
other anchors.

```@example rescale
using EnzymeRates

uni_uni = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E + P
    end
end

params = (k_ES_to_E_P = 3.0, k_ES_to_E_S = 1.0, k_E_S_to_ES = 4.0, Keq = 2.0)

# A turnover of 42 measured at S = 1000 with no product
anchored, spread = rescale_parameter_values(uni_uni, params, [(S = 1e3,) => 42.0])
anchored
```

The rescaled law reproduces the anchor:

```@example rescale
rate_equation(uni_uni, (S = 1e3, P = 0.0), merge(anchored, (E_total = 1.0,)))
```

Only the rate constants changed; `Keq` is the same in `params` and `anchored`.
