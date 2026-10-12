# Normalized vs absolute rate

Enzyme kinetic data is often reported on an arbitrary rate scale. The enzyme
concentration in the assay may be unknown, and figures from different
preparations may reflect different states of purity, so the absolute value of a
measured rate carries no reliable meaning. What survives is the *shape* of the
response: changing the amount of enzyme scales every rate by the same factor and
leaves that shape intact.

EnzymeRates.jl handles both cases through one `FittingProblem` keyword, `absolute`.
With the default `absolute = false` the data are relative: the fit ignores an
arbitrary factor on each group's rates, and the returned rate constants are scaled so
the mechanism's peak turnover is 1. With `absolute = true` the data are per-enzyme
turnover: the fit scores absolute magnitudes and returns the constants as fitted.
`IdentifyRateEquationProblem` takes the same keyword.

## The two modes

### Relative mode (the default)

- **Loss**: each group's log-ratios are **mean-centered before squaring**,
  which removes the arbitrary per-group `E_total` scale. The loss is invariant
  to rescaling all rates in a group by the same factor.
- **Output scale**: the loss cannot see a common factor on every steady-state rate
  constant, so the fitter fixes it after the fit: the returned constants make the
  mechanism's **peak turnover** equal 1. The peak turnover is the largest net forward
  rate per enzyme, `v/E_total`, over every metabolite's concentration — substrates,
  products, inhibitors and regulators — within `peak_range`. A deterministic search in
  log concentrations finds it, the range's corners included. Multiply the returned rate
  constants by a measured peak turnover to put them in absolute units, or rescale them
  to measured kcats (below).

Relative mode is the right choice when your data reports relative rates
(arbitrary units, percent activity, or values normalized per some `E_total`
estimate that you do not trust).

### Absolute mode

Pass `absolute = true`. In absolute mode:

- **Loss**: log-ratios are squared **without centering**. The absolute
  magnitude of the rate is scored, so the data must report per-enzyme turnover
  in consistent units.
- **Output**: the raw fitted parameters, never rescaled. The data fix the absolute
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

## The peak range

`peak_range = (lo, hi)` bounds the concentrations the peak search visits, in the units
of your data. The default `(1e-12, 1e3)` covers physiological concentrations whether
your data are in M, the unit [Fitting tutorial & data format](@ref) recommends, or in
mM. For data in µM, where `1e3` means 1 mM, raise the upper bound. `FittingProblem`,
`IdentifyRateEquationProblem` and `rescale_parameter_values` all take `peak_range`.

Most laws peak inside the range or saturate well below its upper edge. A law still
rising at the upper edge — by more than 10% over the range's last decade — is flagged:
the fit result's `peak_at_range_edge` is `true`, and identify's CSVs record the flag in
a column of the same name. The fit is still scaled, to the turnover at the edge; the
flag says that a wider range would change the scale.

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

- **Measured turnovers.** Pass `kcats`, a vector of `concentrations => turnover` pairs,
  each a measured turnover — a kcat or a specific activity in turnover units — at stated
  concentrations. Metabolites you leave out are zero. Concentrations use the units of
  your data. Turnover is per unit `E_total`, the active-site concentration, so divide a
  kcat quoted per oligomer by its number of active sites first. One anchor makes the
  predicted turnover equal the measured one there; several combine by the geometric mean
  of measured over predicted turnover. The second return value, the spread, is the
  standard deviation of the anchors' `log(measured/predicted)`: 0 for one anchor, and
  large when the anchors disagree with the model or with each other.
- **Peak scale.** Without `kcats`, the function applies the scale relative fits get:
  peak turnover 1 within `peak_range`. Its second return value is `nothing`.

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

The peak scale, which relative fits return:

```@example rescale
peak_scaled, _ = rescale_parameter_values(uni_uni, params)
peak_scaled
```

Only the rate constants changed; `Keq` is the same in both.
