# Normalized vs absolute rate

Enzyme kinetic data is often reported on an arbitrary rate scale. The enzyme
concentration in the assay may be unknown, and figures from different
preparations may reflect different states of purity, so the absolute value of a
measured rate carries no reliable meaning. What survives is the *shape* of the
response: changing the amount of enzyme scales every rate by the same factor and
leaves that shape intact.

EnzymeRates.jl handles both cases through one `FittingProblem` field,
`scale_k_to_kcat`. For a trustworthy per-enzyme turnover scale, set it to
`nothing`; the fit then scores absolute magnitudes. For an arbitrary scale, set
it to a number — the enzyme's kcat if you know it, or simply `1`; the fit then
becomes invariant to that arbitrary factor and pins the returned turnover to the
number you chose. That single field controls two things at once: which loss
formula the fitter uses, and how it normalizes the returned parameters.

## The two modes

### Relative mode (the default)

Pass a positive `Real` — the default is `1.0`. In relative mode:

- **Loss**: each group's log-ratios are **mean-centered before squaring**,
  which removes the arbitrary per-group `E_total` scale. The loss is invariant
  to rescaling all rates in a group by the same factor.
- **Output normalization**: the fitted SS rate constants are rescaled so that
  the forward catalytic rate constant (kcat) equals `scale_k_to_kcat`. With
  the default `1.0`, the returned equation has kcat = 1; the absolute turnover
  is recovered by multiplying by a separately measured kcat value. A custom
  target (e.g. `scale_k_to_kcat = 42.0`) gives kcat = 42 directly.

Relative mode is the right choice when your data reports relative rates
(arbitrary units, percent activity, or values normalized per some `E_total`
estimate that you do not trust to three significant figures).

### Absolute mode

Pass `nothing`. In absolute mode:

- **Loss**: log-ratios are squared **without centering**. The absolute
  magnitude of the rate is scored, so the data must report per-enzyme turnover
  in consistent units.
- **Output**: the raw fitted parameters are returned without kcat
  normalization. The data fixes the absolute scale.

Absolute mode is the right choice when your data is in units of
turnover-per-enzyme and you have reliable `E_total` measurements.

## Rescaling after the fact

`rescale_parameter_values(mechanism, params; scale_k_to_kcat)` is the public
API for kcat normalization. It rescales only the rate constants — the
lowercase-`k` parameters (`kon_…`, `koff_…`, and the steady-state
interconversion `k_…`) — while leaving the binding constants `K`, `Keq`,
`E_total`, the allosteric `L`, and regulatory K values unchanged.

The rate constants are the only parameters that carry time in their units, and
the measured rate is the only thing that supplies a time scale, so it can pin
down only what lives in them. The binding constants, `Keq`, and `E_total` are
concentrations or concentration ratios, and `L` is a pure number — none carry a
time dimension, so the rate constrains them only as far as the response shape
already does. Multiplying every rate constant by a common factor scales kcat by
that factor (kcat is homogeneous of degree one in them), and the K's stay put: a
rapid-equilibrium step stores its binding constant directly, so the rescaling
never touches it, while a steady-state binding constant is the ratio `koff/kon`,
in which the common factor cancels. One uniform rescaling therefore sets kcat to
any target while leaving every other parameter fixed.

`fit_rate_equation` calls this function internally in relative mode. You can
also call it directly on any parameter `NamedTuple` to renormalize after the
fact.

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

# These params have kcat = 3.0 (kon_P_ES is the bottleneck)
params = (koff_P_ES = 6.0, kon_P_ES = 3.0, kon_S_E = 4.0, Keq = 2.0, E_total = 1.0)

rescaled = rescale_parameter_values(uni_uni, params; scale_k_to_kcat = 1.0)
rescaled
```

Notice that only the rate constants changed (`koff_P_ES`, `kon_P_ES`,
`kon_S_E`); `Keq` and `E_total` are unchanged. Calling
`rescale_parameter_values` again on `rescaled` with `scale_k_to_kcat = 1.0`
returns the same values — the internal kcat computation on `rescaled` gives
≈ 1.
