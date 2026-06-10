# Identify tutorial

`identify_rate_equation` takes an `EnzymeReaction` and rate data, enumerates
all biochemically valid mechanisms, fits each to data, and returns the
simplest one that generalizes by leave-one-group-out cross-validation.

This page walks through a fast, fully runnable example. The example uses
noiseless data and a collapsed width-1 beam to finish in seconds; the full
production search uses the wider defaults (`min_beam_width=50`,
`loss_rel_threshold=2.0`, `loss_abs_threshold=0.01`, `max_param_count=20`)
and typically runs for roughly an hour.

```@setup identify_fast
using EnzymeRates
```

## A reaction and a generating mechanism

```@example identify_fast
using EnzymeRates

# The reaction: a reversible uni-uni S ⇌ P.
rxn = @enzyme_reaction begin
    substrates: S[C]
    products:   P[C]
end

# A concrete, non-degenerate uni-uni mechanism to generate the data.
generator = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E + P
    end
end

println("fitted params: ", EnzymeRates.fitted_params(generator))
println("metabolites:   ", metabolites(generator))
```

The mechanism is all-steady-state with three independent parameters:
`koff_P_ES`, `kon_P_ES`, `kon_S_E`; `Keq` is supplied by the user and
`E_total` is absorbed into the rate scale.

## Simulate noiseless data

```@example identify_fast
Keq = 2.0
true_params = (koff_P_ES = 2.5, kon_P_ES = 5.0, kon_S_E = 10.0,
               Keq = Keq, E_total = 1.0)

# Two measurement groups — the minimum for leave-one-group-out CV.
concs  = [(S = 1.0, P = 0.1), (S = 2.0, P = 0.1), (S = 5.0, P = 0.1),
          (S = 1.0, P = 0.5), (S = 2.0, P = 0.5), (S = 5.0, P = 0.5)]
groups = ["G1", "G1", "G1", "G2", "G2", "G2"]

rates = [rate_equation(generator, c, true_params) for c in concs]
data  = (group = groups, Rate = rates,
         S = [c.S for c in concs], P = [c.P for c in concs])
nothing # hide
```

The `data` table has a `:group` column, a `:Rate` column, and one column per
metabolite. Each unique `group` value becomes one cross-validation fold, so at
least two groups are required.

## The data contract

`IdentifyRateEquationProblem` validates the table at construction:

- `:group` and `:Rate` columns must be present.
- One column per substrate, product, and regulator (names match
  `metabolites(mechanism)` exactly).
- Every `Rate` must be nonzero — the loss function works in log space.
- At least two distinct `group` values are required for cross-validation.

`Keq` is a required keyword argument, always user-supplied; the package never
estimates it from data.

## Run the fast search

```@example identify_fast
using OptimizationPyCMA

prob = IdentifyRateEquationProblem(rxn, data; Keq = Keq)

results = identify_rate_equation(prob;
    optimizer          = PyCMAOpt(),
    loss_rel_threshold = 1.0,    # cutoff == best loss …
    loss_abs_threshold = 0.0,    # … no additive slack …
    min_beam_width     = 1,      # … so exactly one survivor per level.
    max_param_count    = 4,      # small cap ⇒ seconds, not hours
    n_restarts         = 3,
    maxtime            = 5.0,
    pmap_function      = map,    # serial; pass `pmap` to distribute
    show_progress      = false,
    save_dir           = mktempdir(),
)
nothing # hide
```

With `loss_rel_threshold=1.0` and `loss_abs_threshold=0.0`, the beam cutoff
equals the best loss at each parameter-count level (`_select_beam` in
`src/identify_rate_equation.jl`), keeping exactly one survivor per level.
`min_beam_width=1` reinforces this: even if the cutoff would admit zero
survivors, one is kept. Together, they collapse the search to a deterministic
single-path trace through the mechanism space.

The full production search uses the wider defaults (`min_beam_width=50`,
`loss_rel_threshold=2.0`, `loss_abs_threshold=0.01`, `max_param_count=20`,
`n_restarts=20`, `maxtime=60.0`) and explores far more candidates in
parallel via `pmap_function=pmap`.

## Read the result

`IdentifyRateEquationResults` has exactly two fields: `best` and `cv_results`.

```@example identify_fast
results.best
```

`results.best` is an `AbstractEnzymeMechanism`. Pass it to
[`rate_equation_string`](@ref) to see the symbolic rate equation:

```@example identify_fast
print(rate_equation_string(results.best))
```

The equation holds the same kinetic information as the generating mechanism
(the same simplified fractional form), which is why the search recovers it on
noiseless data.

`results.cv_results` is a `DataFrame` with one row per candidate equation that
entered cross-validation:

```@example identify_fast
first(results.cv_results, 5)
```

The schema of `cv_results` — its columns, the CV score definition, and the
two-test model-selection rule — is detailed on the [Model selection](@ref)
page.

The enumeration strategy — `init_mechanisms`, the six expansion moves, and
deduplication — is described on the [The enumeration engine](@ref) page.
