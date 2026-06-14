# Solver-agnostic kwarg forwarding for fitting

Date: 2026-06-14

## Problem

`identify_rate_equation` and `fit_rate_equation` forward optimizer
keyword arguments to `Optimization.solve` in a way that only works for
CMA-ES-style solvers:

- `identify_rate_equation` hardcodes `popsize::Int = 200` and
  `verbose::Int = -9` as named kwargs and **always** splats them into
  `Optimization.solve` (via `fitting_kwargs` → `_process_batch` →
  `fit_rate_equation` → `solve`).
- `popsize` and `verbose` are **solver-specific** options (CMA-ES /
  pycma). Per the Optimization.jl docs, the only options handled by the
  wrapper layer across solvers are the *common* options `maxiters`,
  `maxtime`, `abstol`, `reltol`, `callback` — unsupported common options
  *warn*; solver-specific options forwarded to a solver that does not
  accept them are a **hard error**.
- Net effect: the pipeline only runs against solvers that happen to
  accept `popsize`/`verbose`. Any other Optimization.jl solver errors
  out. The README's own `identify_rate_equation` example would fail with
  `CMAEvolutionStrategyOpt` for exactly this reason.

`fit_rate_equation` is already partly clean — it force-passes only
`maxtime` and slurps everything else through a `kwargs...` catch-all —
but it inherits whatever `identify_rate_equation` hands it, and its
catch-all silently leaks any stray kwarg to `solve`.

`lb`/`ub` are **not** part of this problem: they go to the
`OptimizationProblem` *constructor*, not `solve`, and their length is
per-mechanism (`length(fitted_params(...))`), which is why they live at
the `fit_rate_equation` level and `identify_rate_equation` does not
expose them. They are unchanged by this work.

## Goal

Support **every** Optimization.jl solver. Separate the kwargs the
package owns an opinionated default for (and that are universally safe)
from solver-specific options, which the caller supplies explicitly and
which are forwarded verbatim.

## Background: the `OptimizationCMAEvolutionStrategy` wrapper

The recommended optimizer is moving from `PyCMAOpt` (OptimizationPyCMA)
to `CMAEvolutionStrategyOpt` (OptimizationCMAEvolutionStrategy) — pure
Julia, no Python dependency, which removes a source of HPC friction.

Reading the installed wrapper source (`OptimizationCMAEvolutionStrategy`
versions `MySJY` and `1x519` — identical behavior) and the underlying
`CMAEvolutionStrategy.minimize`:

- The wrapper's `__map_optimizer_args` has a **fixed** kwarg signature —
  `(callback, maxiters, maxtime, abstol, reltol, verbose)` — with **no
  `kwargs...` catch-all**. It maps: `maxiters → maxiter`, `maxtime →
  maxtime`, `abstol → ftol`, `verbose::Bool → logger verbosity (1 or
  0)`; it **warns and ignores** `reltol`.
- It **never forwards `popsize`** to `minimize` and has no parameter for
  it, so `solve(prob, CMAEvolutionStrategyOpt(); popsize=…)` →
  `"unsupported keyword argument popsize"`. `minimize` itself supports
  `popsize` (default `4 + floor(3·ln(N))` ≈ 8 for 5 params), but it is
  unreachable through the wrapper.
- It **hardcodes the initial step size `sigma0 = 0.1`** (the
  `minimize(_loss, u0, 0.1; …)` call), so that is not settable either.
- `verbose` is supported **only as a `Bool`** (an integer like the old
  `-9` fails the `::Bool` signature). The default `false` → verbosity 0
  → **silent by default**, so the old `verbose=-9` suppression is
  unnecessary on this solver.

Consequences for switching the recommended solver to
`CMAEvolutionStrategyOpt`:

1. The `verbose` problem disappears — fits are quiet by default.
2. `popsize`/`sigma0` robustness tuning is **inaccessible** through the
   wrapper. Fits use the default `popsize ≈ 8` and `sigma0 = 0.1`. This
   is a real reduction in per-fit global-search robustness versus
   `PyCMAOpt(popsize=200)`; it is partly compensable with more
   `n_restarts` but is not equivalent. Restoring `popsize` control
   requires an **upstream fix** to the wrapper (see Out of Scope).

This is orthogonal to the redesign: the `solver_kwargs` bag forwards
whatever the caller puts; for `CMAEvolutionStrategyOpt` there simply are
few useful knobs to forward today.

## Design

### Kwarg taxonomy (both `fit_rate_equation` and `identify_rate_equation`)

**Named, with package defaults — EnzymeRates fitting knobs:**

- `n_restarts::Int = 20` — multi-start count (an EnzymeRates concept,
  not an Optimization.jl option).
- `lb`, `ub` — log-space bounds → `OptimizationProblem` constructor.
  `fit_rate_equation` only; per-mechanism length.

**Named, with defaults — Optimization.jl common solver options** (the
full documented common set; safe across solvers — unsupported ones warn,
never error):

- `maxtime::Real = 60.0` — always forwarded (the real per-fit budget).
- `maxiters::Integer = 10_000_000` — always forwarded; effectively
  "let `maxtime` govern", preserving current behavior.
- `abstol::Union{Real,Nothing} = nothing`
- `reltol::Union{Real,Nothing} = nothing`
- `callback = nothing`

  The optional three default to `nothing` and are forwarded **only when
  set**, so each solver keeps its own default otherwise. (On
  `CMAEvolutionStrategyOpt`, `reltol` is a documented no-op — the wrapper
  warns; that is expected and acceptable for a common option.)

**Pass-through bag — solver-specific options:**

- `solver_kwargs = (;)` — a `NamedTuple` (or pairs) forwarded **verbatim**
  to `Optimization.solve`. Holds anything solver-specific (`popsize`,
  `verbose`, `sigma0`, `seed`, `multi_threading`, …) and may also override
  a common default. The caller owns matching its contents to the chosen
  solver.

### Forwarding / merge semantics

In `fit_rate_equation`, assemble the solve kwargs and merge so the bag
can override a named default without a duplicate-keyword error:

```julia
common = (; maxtime, maxiters)
abstol   === nothing || (common = (; common..., abstol))
reltol   === nothing || (common = (; common..., reltol))
callback === nothing || (common = (; common..., callback))
solve_kw = merge(common, solver_kwargs)   # solver_kwargs wins on conflict
...
sol = Optimization.solve(prob, optimizer; solve_kw...)
```

`merge((; maxtime, …), solver_kwargs)` is last-wins: a key present in
both takes the `solver_kwargs` value, so a caller can override `maxtime`
or `maxiters` through the bag without Julia raising a duplicate-keyword
error.

### Threading through `identify_rate_equation`

`identify_rate_equation` exposes the same named commons + `solver_kwargs`
and bundles them into the existing `fitting_kwargs` NamedTuple:

```julia
fitting_kwargs = (; n_restarts, maxtime, maxiters,
                    abstol, reltol, callback, solver_kwargs)
```

`fitting_kwargs` already flows down through `_beam_search` →
`_process_batch` → `fit_rate_equation` and `_cv_model_selection` →
`_loocv` → `fit_rate_equation` via their existing `kwargs...` forwarding.
Because `solver_kwargs` is a single named entry inside `fitting_kwargs`,
those intermediate functions need **no signature changes** — only
`fit_rate_equation` and `identify_rate_equation` change.

### Clean break (no back-compat)

- Remove `popsize` and `verbose` as named kwargs of
  `identify_rate_equation`.
- Remove `identify_rate_equation`'s `optim_kwargs...` catch-all and
  `fit_rate_equation`'s trailing `kwargs...` catch-all. The API becomes
  fully explicit: the named commons plus `solver_kwargs`. A stray/unknown
  kwarg now errors at the EnzymeRates boundary (clear "unsupported
  keyword argument") instead of silently leaking to `solve`.
- Existing callers passing `popsize=…`/`verbose=…` migrate to
  `solver_kwargs=(; popsize=…, verbose=…)`. `maxiters=…` stays a named
  kwarg.

## Recommended usage & docs (`CMAEvolutionStrategyOpt` only)

README and docstrings recommend `CMAEvolutionStrategyOpt` as the single
optimizer; PyCMA is dropped from the docs.

- `using OptimizationCMAEvolutionStrategy` (was `OptimizationPyCMA`).
- `fit_rate_equation(fp, CMAEvolutionStrategyOpt(); n_restarts=3,
  maxtime=5.0)` — drop `popsize=50` (unsupported by this wrapper; fits
  use the default population). No `solver_kwargs` needed for the common
  case; fits are silent by default.
- `identify_rate_equation(prob; optimizer=CMAEvolutionStrategyOpt(),
  max_param_count=10, pmap_function=map)` — now works because the
  defaults no longer inject `popsize`/`verbose`.
- Document `solver_kwargs` in both docstrings: what it is, that it is
  forwarded verbatim to `Optimization.solve`, that contents must match
  the chosen solver, and a concrete example (e.g. for a solver that
  supports it, `solver_kwargs=(; popsize=200)`).

## Testing

Goal: prove solver-agnosticism across multiple solvers, and prove the
forwarding/clean-break boundary. Add `OptimizationCMAEvolutionStrategy`
to `test/Project.toml` (and the `test` target list) so the recommended
path is exercised.

1. **Recommended-path / regression (the bug):** `fit_rate_equation` and
   `identify_rate_equation` run to success with `CMAEvolutionStrategyOpt()`
   and **default** (empty) `solver_kwargs`. This fails on `main` today
   (force-injected `popsize=200`), so it is the core regression test.
2. **Forwarding to a supporting solver:** a `PyCMAOpt()` test passing
   `solver_kwargs=(; popsize=…, verbose=-9)` succeeds — proves the bag is
   forwarded and honored where supported. (Keep OptimizationPyCMA as a
   test-only dep; the Python dep is acceptable in CI even though it is
   being dropped from user-facing docs.)
3. **Verbatim forwarding / wrapper limitation:** passing
   `solver_kwargs=(; popsize=200)` to `CMAEvolutionStrategyOpt()` raises
   the expected "unsupported keyword argument" error — proves the bag is
   forwarded unmodified and documents the wrapper gap.
4. **Clean break:** passing `popsize=…`/`verbose=…` as **top-level**
   kwargs to `identify_rate_equation` now raises "unsupported keyword
   argument" (they are no longer accepted there).
5. **Existing call-site migration:**
   - `test/test_identify_rate_equation.jl:250` `popsize=200` →
     `solver_kwargs=(; popsize=200)`.
   - `test/test_identify_rate_equation.jl:389-390`
     `maxiters=500, popsize=40, verbose=-9` →
     `maxiters=500, solver_kwargs=(; popsize=40, verbose=-9)`.
   - The `BBO_…` `fit_rate_equation` tests (`test/test_fitting.jl`) and
     all default-path identify tests need no kwarg changes; they already
     pass only `n_restarts`/`maxtime`. They additionally serve as
     non-CMA solver coverage.
6. **README test:** the in-test README block must use
   `CMAEvolutionStrategyOpt` and the migrated kwargs and stay green.

## Files touched

- `src/fitting.jl` — `fit_rate_equation` signature + solve-kwarg
  assembly; docstring.
- `src/identify_rate_equation.jl` — `identify_rate_equation` signature
  (drop `popsize`/`verbose`/`optim_kwargs`, add common kwargs +
  `solver_kwargs`); `fitting_kwargs` assembly; docstring. No changes to
  `_process_batch`/`_loocv`/`_beam_search`/`_cv_model_selection`.
- `test/Project.toml` — add `OptimizationCMAEvolutionStrategy`.
- `test/test_fitting.jl`, `test/test_identify_rate_equation.jl` — new
  solver-agnostic tests; migrate the `popsize`/`verbose` call-sites.
- `README.md` — switch to `CMAEvolutionStrategyOpt`; migrate kwargs;
  document `solver_kwargs`.

## Out of scope (tracked follow-up)

An upstream PR to `OptimizationCMAEvolutionStrategy` to forward `popsize`
and `sigma0` (and add a `kwargs...` catch-all to pass other
`CMAEvolutionStrategy.minimize` options through) so that
`solver_kwargs=(; popsize=200)` works on `CMAEvolutionStrategyOpt` once
merged. Not a deliverable of this spec.
