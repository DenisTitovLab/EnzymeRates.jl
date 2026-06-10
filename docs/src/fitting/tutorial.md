# Fitting tutorial & data format

This page shows how to fit a rate equation to experimental data from start to
finish. It covers the data format, constructing a `FittingProblem`, running the
multi-start optimizer, and reading the result.

## Data format

`FittingProblem` accepts any
[Tables.jl](https://github.com/JuliaData/Tables.jl)-compatible source — a
`DataFrame`, a CSV-loaded table, or a plain `NamedTuple` of equal-length
vectors. The required columns are:

| Column | Type | Description |
|--------|------|-------------|
| `group` | any | Groups rows that share the same `E_total`. Leave-one-group-out cross-validation folds on this column. |
| `Rate` | nonzero `Real` | Measured reaction rate. Must be nonzero — the loss works in log space. |
| one per metabolite | `Real` | Concentration of each metabolite. Column names must match `metabolites(mechanism)` exactly. |

Call `metabolites(mechanism)` to find which concentration columns your data
needs before constructing the problem.

```@example fitting
using EnzymeRates

uni_uni = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E + P
    end
end

metabolites(uni_uni)
```

A missing required column, a missing metabolite column, or a zero `Rate` each
raises an `ErrorException` at construction — the check runs before any fitting.

## Building the `FittingProblem`

The example below generates synthetic data by evaluating the true rate equation
on a concentration grid, then wraps it in a `FittingProblem`.

```@example fitting
# Independent fitted params for the reduced equation: koff_P_ES, kon_P_ES, kon_S_E
true_params = (koff_P_ES = 3.0, kon_P_ES = 5.0, kon_S_E = 10.0, Keq = 2.0, E_total = 1.0)

concs_list = [(S = s, P = p) for s in (0.5, 1.0, 2.0, 5.0, 10.0) for p in (0.1, 0.5)]

data = (
    group = fill("G1", length(concs_list)),
    Rate  = [rate_equation(uni_uni, c, true_params) for c in concs_list],
    S     = [c.S for c in concs_list],
    P     = [c.P for c in concs_list],
)

fp = FittingProblem(uni_uni, data; Keq = 2.0)
```

`Keq` is a required keyword argument and is always user-supplied — the package
never estimates it from data. The constructor also accepts a concrete
`Mechanism` or `AllostericMechanism` and compiles it once at construction, so
the fitting hot path pays no compilation overhead.

## Running the fit

`fit_rate_equation` requires an explicit optimizer — the base package depends
only on `Optimization.jl` and ships no solver backend. Install
`OptimizationPyCMA` (the recommended choice; needs Python and `pycma`) or
`OptimizationBBO` (a pure-Julia alternative) and pass its optimizer object. See
[Loss & optimizers](@ref) for details.

```@example fitting
using OptimizationPyCMA

result = fit_rate_equation(
    fp, PyCMAOpt();
    n_restarts = 3, maxtime = 10.0,
)

keys(result)
```

```@example fitting
result.retcode isa Symbol
```

```@example fitting
result.params
```

The fit runs `n_restarts` independent optimizations from random initial points
and returns the best. Because the optimizer is stochastic, `result.params` and
`result.loss` vary from run to run.

Check `result.retcode === :Success` to confirm the optimizer converged on its
own criteria. Any other value — `:MaxTime`, `:Failure`, `:NoFiniteLoss` — means
the fit is un-converged. When time or restart budgets are tight, `:MaxTime` is
common; increase `n_restarts` and `maxtime` for production fits.

The returned `params` reflect kcat normalization: with the default
`scale_k_to_kcat = 1.0`, the SS rate constants are rescaled so that kcat = 1.
[Normalized vs absolute rate](@ref) explains `scale_k_to_kcat`;
[Loss & optimizers](@ref) covers the loss function and how to supply your own
optimizer.
