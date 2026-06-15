# Loss & optimizers

This page explains the loss function that `fit_rate_equation` minimizes and how
to choose an optimization algorithm through Optimization.jl.

## The loss function

`loss!(x, fp)` computes a **log-ratio loss**. Parameters in `x` are in log
space — the actual rate constants are `exp.(x)`. For each of the `N` data points
it squares the log-ratio of predicted to measured rate, then sums over all
points and divides by `N`:

```
loss = (1/N) Σᵢ (log(|predicted rateᵢ|) − log(|measured rateᵢ|))²
```

### Centered vs uncentered

The loss formula switches on `fp.scale_k_to_kcat` (see
[Normalized vs absolute rate](@ref)):

- **Relative mode** (`scale_k_to_kcat` is a `Real`): each group's log-ratios
  are **mean-subtracted before squaring**. This makes the loss invariant to
  multiplying all rates in a group by the same factor.
- **Absolute mode** (`scale_k_to_kcat === nothing`): the log-ratios are
  **squared without centering**. The absolute magnitude is scored, so the
  data must be in consistent per-enzyme turnover units.

### Sign-mismatch penalty

When the predicted rate sign disagrees with the measured sign, or when the
prediction is zero, the log-ratio is replaced by a sentinel value of `10.0`
in the working buffer. After the main loop, a flat `100.0` penalty per
mismatched point is added to the total. In centered mode this prevents an
all-mismatch group from contributing zero loss (the uniform sentinel would
cancel under mean-subtraction); the post-hoc penalty keeps it positive.

!!! note "When the sign penalty can bite"
    This penalty is harmless for most enzyme kinetic data. An assay usually
    omits at least one substrate or one product, so the net rate keeps one sign
    across the whole dataset and stays clear of zero. When every substrate and
    product is present, though, the net rate can pass through zero, and the
    penalty can steer the fit toward parameters that reproduce the sign of
    near-zero rates and fit the rest of the data poorly. Watch for this when
    your measurements straddle the reaction's equilibrium.

## Optimizers

EnzymeRates fits with any optimization algorithm available through
[Optimization.jl](https://docs.sciml.ai/Optimization/stable/).
`fit_rate_equation` wraps `loss!` into an `Optimization.OptimizationFunction`,
builds an `OptimizationProblem`, and calls `Optimization.solve(prob, optimizer;
…)` with the `optimizer` you pass. EnzymeRates depends only on Optimization.jl
and ships no solver itself, so install the sub-package for the algorithm you
want and pass its optimizer object. Two pure-Julia choices cover most needs:

```julia
] add OptimizationCMAEvolutionStrategy   # CMA-ES — recommended
] add OptimizationBBO                     # BBO differential evolution
```

`CMAEvolutionStrategyOpt()` from `OptimizationCMAEvolutionStrategy` is the
recommended optimizer and the one the [`identify_rate_equation`](@ref) search
uses: CMA-ES (Covariance Matrix Adaptation Evolution Strategy) handles the
correlated, non-convex landscapes that the Haldane/Wegscheider parameter
reduction produces. `BBO_adaptive_de_rand_1_bin_radiuslimited()` from
`OptimizationBBO` is a solid alternative for cross-checking a fit against an
independent global optimizer. Any other Optimization.jl solver works as well.

## Passing solver options

`fit_rate_equation` separates two kinds of solver options:

- **Common options** are named keyword arguments: `maxtime`, `maxiters`,
  `abstol`, `reltol`, and `callback`. These are the
  [Optimization.jl common solver options](https://docs.sciml.ai/Optimization/stable/API/solve/)
  every solver understands. `maxtime` and `maxiters` are always forwarded to
  `Optimization.solve`; `abstol`, `reltol`, and `callback` are forwarded only
  when set, so each solver keeps its own default otherwise.
- **Solver-specific options** go in `solver_kwargs`, a `NamedTuple` forwarded
  verbatim to `Optimization.solve`. When a key appears in both a named common
  option and `solver_kwargs`, the `solver_kwargs` value wins.

`n_restarts` is separate from both: it sets how many independent optimizations
`fit_rate_equation` runs from different random starting points, keeping the best
result. Finding the global minimum of a rate equation — a ratio of polynomials —
is NP-hard, so any single run can settle in a local minimum; restarting from
independent points raises the chance of reaching the global optimum. In our
tests on complex rate equations with datasets of 500–1000 points,
`n_restarts = 10–20` returned the same loss on every run.
