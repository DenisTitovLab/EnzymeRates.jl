# Normalized vs absolute rate

The single `scale_k_to_kcat` field on `FittingProblem` controls two things at
once: which loss formula the fitter uses, and how it normalizes the returned
parameters.

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

A non-positive `Real` (e.g. `scale_k_to_kcat = 0.0`) raises an error at
construction.

## Confirming the default

```jldoctest
julia> using EnzymeRates

julia> m = @enzyme_mechanism begin
           substrates: S
           products:   P
           steps: begin
               E + S <--> E(S)
               E(S) <--> E + P
           end
       end;

julia> data = (group = ["G1"], Rate = [1.0], S = [1.0], P = [0.1]);

julia> FittingProblem(m, data; Keq = 2.0).scale_k_to_kcat === 1.0
true
```

## Rescaling after the fact

`rescale_parameter_values(mechanism, params; scale_k_to_kcat)` is the public
API for kcat normalization. It rescales only the SS rate constants — the
`Kon`, `Koff`, `Kfor`, `Krev` family — while leaving RE binding constants,
`Keq`, `E_total`, allosteric `L`, and regulatory K values unchanged. kcat is
homogeneous degree-1 in the SS rate constants and independent of the RE K's,
so rescaling the SS k's uniformly sets kcat to any target without disturbing
anything else.

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

Notice that only the SS rate constants changed (`koff_P_ES`, `kon_P_ES`,
`kon_S_E`); `Keq` and `E_total` are unchanged. Calling
`rescale_parameter_values` again on `rescaled` with `scale_k_to_kcat = 1.0`
returns the same values — the internal kcat computation on `rescaled` gives
≈ 1.
