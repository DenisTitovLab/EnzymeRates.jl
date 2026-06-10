# Loss & optimizers

This page explains the loss function that `fit_rate_equation` minimizes and the
bring-your-own-optimizer design.

## The loss function

`loss!(x, fp)` computes a **log-ratio loss**. Parameters in `x` are in log
space — the actual rate constants are `exp.(x)`. For each data point the
working value is:

```
log(|predicted rate|) − log(|measured rate|)
```

The final loss is the sum of squares of these values divided by the number of
data points.

### Centered vs uncentered

The loss formula switches on `fp.scale_k_to_kcat` (see
[Normalized vs absolute rate](@ref)):

- **Relative mode** (`scale_k_to_kcat` is a `Real`): each group's log-ratios
  are **mean-subtracted before squaring**. This removes the per-group
  `E_total` scale, so the loss is invariant to multiplying all rates in a
  group by the same factor.
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

### Zero-allocation hot path

`loss!` is allocation-free. `FittingProblem` pre-allocates a `log_ratios_buffer`
at construction, and `loss!` reuses it on every call. The loss depends on
`rate_equation`, which is itself allocation-free and sub-100 ns per call for
every mechanism in the package. These two contracts together are what makes
multi-start fitting over millions of loss evaluations practical.

## Bring your own optimizer

**EnzymeRates depends only on `Optimization.jl`; it ships no solver backend.**

`fit_rate_equation` wraps `loss!` into an `Optimization.OptimizationFunction`,
constructs an `OptimizationProblem`, and calls `Optimization.solve(prob,
optimizer; maxtime, kwargs...)`. The `optimizer` argument is whatever
Optimization.jl solver object the caller supplies. Install one of the
Optimization.jl solver sub-packages and pass its optimizer object.

Install one of these — `OptimizationPyCMA` (the recommended choice; needs
Python and `pycma`) or `OptimizationBBO` (a tested, pure-Julia alternative):

```julia
] add OptimizationPyCMA   # CMA-ES — recommended; requires Python and pycma
] add OptimizationBBO     # BBO differential evolution — pure Julia, no Python needed
```

### Recommended: `PyCMAOpt()` from `OptimizationPyCMA`

CMA-ES (Covariance Matrix Adaptation Evolution Strategy) is the recommended
optimizer for rate-equation fitting and the one the identify pipeline uses.
It handles the correlated, non-convex landscapes that arise from the
Haldane/Wegscheider parameter reduction well. It requires Python with
`pycma` installed, which `OptimizationPyCMA` arranges via `CondaPkg.jl`.

### Alternative: `BBO_adaptive_de_rand_1_bin_radiuslimited()` from `OptimizationBBO`

A pure-Julia differential-evolution optimizer with no Python dependency. A
good choice for environments where Python is unavailable.

## Example: running a fit

```@example optimizers
using EnzymeRates, OptimizationPyCMA

uni_uni = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S <--> E(S)
        E(S) <--> E + P
    end
end

true_params = (koff_P_ES = 3.0, kon_P_ES = 5.0, kon_S_E = 10.0, Keq = 2.0, E_total = 1.0)
concs_list = [(S = s, P = p) for s in (0.5, 1.0, 2.0, 5.0, 10.0) for p in (0.1, 0.5)]
data = (
    group = fill("G1", length(concs_list)),
    Rate  = [rate_equation(uni_uni, c, true_params) for c in concs_list],
    S     = [c.S for c in concs_list],
    P     = [c.P for c in concs_list],
)

fp = FittingProblem(uni_uni, data; Keq = 2.0)

result = fit_rate_equation(
    fp, PyCMAOpt();
    n_restarts = 3, maxtime = 10.0,
)

(keys(result), result.retcode isa Symbol)
```

The same call with the pure-Julia BBO optimizer (no Python needed):

```julia
using OptimizationBBO

result = fit_rate_equation(
    fp, BBO_adaptive_de_rand_1_bin_radiuslimited();
    n_restarts = 3, maxtime = 10.0,
)
```

## Passing solver-specific options

Solver keywords flow through `fit_rate_equation`'s `kwargs...` directly to
`Optimization.solve`. `maxtime` and `n_restarts` are arguments of
`fit_rate_equation` itself; solver-specific knobs such as `popsize` and
`maxiters` go in `kwargs`:

```julia
result = fit_rate_equation(
    fp, PyCMAOpt();
    n_restarts = 10, maxtime = 60.0, popsize = 50,
)
```

`n_restarts` and `maxtime` are the primary budget controls. `popsize` is
CMA-ES-specific and controls the population size per generation; larger values
slow each generation but improve coverage of high-dimensional spaces.
