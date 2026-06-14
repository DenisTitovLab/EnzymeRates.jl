# Solver-agnostic kwarg forwarding — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `fit_rate_equation` / `identify_rate_equation` work with every Optimization.jl solver by separating named common solver options from a `solver_kwargs` pass-through bag, and switch the recommended/tested solver from PyCMA to CMAEvolutionStrategy.

**Architecture:** Named kwargs (`n_restarts`, `lb`, `ub`, plus the Optimization.jl common options `maxtime`, `maxiters`, `abstol`, `reltol`, `callback`) keep package defaults; everything solver-specific goes in `solver_kwargs::NamedTuple`, merged into the `Optimization.solve` call (`merge((; maxtime, maxiters, …), solver_kwargs)`, bag wins on conflict). Clean break: `popsize`/`verbose` named kwargs and the `kwargs...`/`optim_kwargs...` catch-alls are removed. Source spec: `docs/superpowers/specs/2026-06-14-solver-kwargs-passthrough-design.md`.

**Tech Stack:** Julia, Optimization.jl, OptimizationCMAEvolutionStrategy (new test dep), OptimizationBBO (kept), Test/Aqua/JET.

---

## Files

- **Modify** `src/fitting.jl:184-252` — `fit_rate_equation` signature, solve-kwarg assembly, docstring.
- **Modify** `src/identify_rate_equation.jl:112-135,179-197` — `identify_rate_equation` docstring + signature (drop `popsize`/`verbose`/`optim_kwargs`, add common kwargs + `solver_kwargs`); `fitting_kwargs` assembly. Body unchanged; `_process_batch`/`_loocv`/`_beam_search`/`_cv_model_selection` unchanged.
- **Modify** `Project.toml` — `[extras]`/`[compat]`/`[targets]`: add `OptimizationCMAEvolutionStrategy` (Task 1), remove `OptimizationPyCMA` (Task 5).
- **Modify** `test/test_identify_rate_equation.jl` — migrate all PyCMA usages to CMAEvolutionStrategy; rework the three `beam_fraction` tests; add regression + clean-break tests.
- **Modify** `test/test_fitting.jl` — add fit-level forwarding/merge/clean-break testset (BBO tests unchanged).
- **Modify** `README.md:81,118-119,124-132,167-173` — switch to `CMAEvolutionStrategyOpt`, drop `popsize`, document `solver_kwargs`.

## Running tests

**Authoritative (every task's final check):**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```
This is the full suite (Aqua + JET + compile budget + all files) — slow/cold; pays precompile each run.

**Optional fast iteration** (one-time install of TestEnv into your global env, then focused file includes):
```bash
julia -e 'using Pkg; Pkg.add("TestEnv")'   # one-time, global
```
Focused run of the identify tests:
```bash
julia --project -e '
using TestEnv; TestEnv.activate()
using EnzymeRates, Test, DataFrames, CSV, Random, Statistics, LinearAlgebra
using OptimizationCMAEvolutionStrategy, OptimizationBBO
include("test/mechanism_definitions_for_test_enzyme_derivation.jl")
include("test/test_identify_rate_equation.jl")'
```
Focused run of the fitting tests: swap the last `include` to `test/test_fitting.jl`.

---

## Task 1: Add OptimizationCMAEvolutionStrategy as a test dependency

**Files:**
- Modify: `Project.toml`

`OptimizationCMAEvolutionStrategy` UUID = `bd407f91-200f-4536-9381-e4ba712f53f8`, installed version 0.3.x. Keep `OptimizationPyCMA` for now (removed in Task 5) so the suite stays green during migration.

- [ ] **Step 1: Add the `[compat]` entry**

In `Project.toml`, in the `[compat]` block, add a line (alphabetical, between `OptimizationBBO` and `OptimizationPyCMA`):

```toml
OptimizationCMAEvolutionStrategy = "0.3"
```

- [ ] **Step 2: Add the `[extras]` entry**

In the `[extras]` block, add (after the `OptimizationBBO` line):

```toml
OptimizationCMAEvolutionStrategy = "bd407f91-200f-4536-9381-e4ba712f53f8"
```

- [ ] **Step 3: Add to the test target**

Replace the `[targets]` line:

```toml
test = ["Test", "OrdinaryDiffEqFIRK", "OptimizationBBO", "OptimizationPyCMA", "Aqua", "JET"]
```

with:

```toml
test = ["Test", "OrdinaryDiffEqFIRK", "OptimizationBBO", "OptimizationCMAEvolutionStrategy", "OptimizationPyCMA", "Aqua", "JET"]
```

- [ ] **Step 4: Resolve and run the full suite**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS (the new dep is just made available; nothing else changed). Resolution must succeed (no version conflict on `OptimizationCMAEvolutionStrategy = "0.3"`).

- [ ] **Step 5: Commit**

```bash
git add Project.toml
git commit -m "Add OptimizationCMAEvolutionStrategy as a test dependency"
```

---

## Task 2: Rewrite source + migrate identify tests to CMAEvolutionStrategy

This is the coupled core: `identify_rate_equation` force-injects `popsize=200`/`verbose=-9`, and it forwards them to `fit_rate_equation`, so both functions and the identify tests must change together to keep the suite green. The migrated tests fail on the unmodified source (proving the bug), then pass after the rewrite.

**Files:**
- Modify: `test/test_identify_rate_equation.jl`
- Modify: `src/fitting.jl:184-252`
- Modify: `src/identify_rate_equation.jl:112-135,179-197`

- [ ] **Step 1: Migrate every PyCMA usage in `test/test_identify_rate_equation.jl`**

Apply these exact replacements:

(a) Line 9 — the import:
```julia
using OptimizationPyCMA
```
→
```julia
using OptimizationCMAEvolutionStrategy
```

(b) Line 81 — the shared optimizer binding:
```julia
    pycma_opt = PyCMAOpt()
```
→
```julia
    cmaes_opt = CMAEvolutionStrategyOpt()
```

(c) Line 229 (pipeline) and line 365 (`save_dir` test) — both read `        optimizer=pycma_opt,`. Replace **both** occurrences with `        optimizer=cmaes_opt,` (use replace-all on the exact line `optimizer=pycma_opt,`).

(d) The mechanism-recovery fit (lines 247-250):
```julia
        fit_best = fit_rate_equation(
            fp_best, pycma_opt;
            n_restarts=3, maxtime=10.0,
            popsize=200)
```
→
```julia
        fit_best = fit_rate_equation(
            fp_best, cmaes_opt;
            n_restarts=3, maxtime=10.0)
```

(e) The `_loocv` per-fold-scores fit (lines 386-390):
```julia
        scores = EnzymeRates._loocv(
            m, prob;
            optimizer=PyCMAOpt(),
            n_restarts=2, maxtime=2.0,
            maxiters=500, popsize=40, verbose=-9)
```
→
```julia
        scores = EnzymeRates._loocv(
            m, prob;
            optimizer=CMAEvolutionStrategyOpt(),
            n_restarts=2, maxtime=2.0,
            maxiters=500)
```

(f) The `_loocv` loud-on-failure test (lines 416-420):
```julia
        # An unsupported optimizer kwarg makes every fold fit throw; _loocv
        # must NOT swallow it (the old behavior returned Float64[]).
        @test_throws Exception EnzymeRates._loocv(
            m, prob; optimizer=PyCMAOpt(),
            n_restarts=1, maxtime=1.0, beam_fraction=0.5)
```
→
```julia
        # An unrecognized kwarg (`beam_fraction`) is rejected by
        # `fit_rate_equation` per fold; _loocv must NOT swallow that error
        # (the old behavior returned Float64[]).
        @test_throws Exception EnzymeRates._loocv(
            m, prob; optimizer=CMAEvolutionStrategyOpt(),
            n_restarts=1, maxtime=1.0, beam_fraction=0.5)
```

(g) The all-base-fail testset — title (line 495):
```julia
@testset "beam_fraction kwarg removed: passing it errors" begin
```
→
```julia
@testset "all base fits fail: failure CSV written, then raises" begin
```
its comment (lines 496-502):
```julia
    # `beam_fraction` is not a recognized kwarg; it is forwarded to
    # `Optimization.solve`, which rejects it. Per-mechanism fit failures
    # are isolated in `_process_batch`, so an unknown optimizer kwarg
    # fails every fit; the base tier is then empty and the pipeline raises.
    # The contract under test is that passing it errors (no silent
    # acceptance), AND that the all-base-fail path persists the failure
    # rows to `initial_mechanisms.csv` before raising (for cluster debugging).
```
→
```julia
    # A `solver_kwargs` option the optimizer rejects (`popsize` is
    # unsupported by CMAEvolutionStrategy) is forwarded verbatim to
    # `Optimization.solve` and makes every fit throw. Per-mechanism fit
    # failures are isolated in `_process_batch`, so the base tier is then
    # empty and the pipeline raises. The contract under test is that the
    # all-base-fail path persists the failure rows to
    # `initial_mechanisms.csv` before raising (for cluster debugging).
```
and its call (lines 513-516):
```julia
    @test_throws ErrorException identify_rate_equation(
        prob; beam_fraction=0.5,
        optimizer=PyCMAOpt(),
        n_restarts=1, maxtime=1.0, save_dir=tmp)
```
→
```julia
    @test_throws ErrorException identify_rate_equation(
        prob; solver_kwargs=(; popsize=200),
        optimizer=CMAEvolutionStrategyOpt(),
        n_restarts=1, maxtime=1.0, save_dir=tmp)
```

(h) The `_process_batch` testset — two `optimizer=PyCMAOpt(),` lines (972, 983) become `optimizer=CMAEvolutionStrategyOpt(),`; and the failure-path block (lines 988-992):
```julia
    # config error (unknown optimizer kwarg) → every fit throws → all failures,
    # no entries; each failure carries a non-empty error string.
    fail_entries, fail_failures = EnzymeRates._process_batch(ms, prob;
        pmap_function=map, optimizer=PyCMAOpt(),
        max_param_count=20, n_restarts=1, maxtime=1.0, beam_fraction=0.5)
```
→
```julia
    # config error (solver rejects an option) → every fit throws → all
    # failures, no entries; each failure carries a non-empty error string.
    fail_entries, fail_failures = EnzymeRates._process_batch(ms, prob;
        pmap_function=map, optimizer=CMAEvolutionStrategyOpt(),
        max_param_count=20, n_restarts=1, maxtime=1.0,
        solver_kwargs=(; popsize=200))
```

- [ ] **Step 2: Add the regression test** (drives the source change)

Add this standalone testset to the end of `test/test_identify_rate_equation.jl` (after the `_process_batch` / `_ingest!` testsets, at top level):

```julia
@testset "identify runs on a solver that rejects popsize (regression)" begin
    # On the pre-fix source, identify_rate_equation force-injected
    # popsize=200, which CMAEvolutionStrategy rejects → every base fit
    # failed → ErrorException. The default path must now run cleanly on a
    # solver that does not accept popsize.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (group = ["G1", "G1", "G2", "G2"],
            Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0],
            P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    tmp = mktempdir()
    results = identify_rate_equation(prob;
        optimizer=CMAEvolutionStrategyOpt(),
        min_beam_width=1, loss_rel_threshold=1.0, loss_abs_threshold=0.0,
        max_param_count=6, n_cv_candidates=1, n_restarts=1, maxtime=1.0,
        pmap_function=map, save_dir=tmp, show_progress=false)
    @test results isa IdentifyRateEquationResults
end
```

- [ ] **Step 3: Run to verify RED on the unmodified source**

Run: `julia --project -e 'using Pkg; Pkg.test()'` (or the focused identify recipe above).
Expected: FAIL. The main pipeline (`results = identify_rate_equation(prob; … optimizer=cmaes_opt …)` near line 221) and the new regression testset throw `ErrorException` with a message containing `Every base-tier fit failed` / `unsupported keyword argument popsize` — the bug.

- [ ] **Step 4: Rewrite `fit_rate_equation`** (`src/fitting.jl`)

Replace the entire docstring + function (lines 184-252) with:

```julia
"""
    fit_rate_equation(fp::FittingProblem, optimizer;
        n_restarts=20, maxtime=60.0, maxiters=10_000_000,
        abstol=nothing, reltol=nothing, callback=nothing,
        lb=fill(-15.0, length(fitted_params(fp.mechanism))),
        ub=fill(15.0, length(fitted_params(fp.mechanism))),
        solver_kwargs=(;))

Fit rate constants by minimizing `loss!` using Optimization.jl.

Runs `n_restarts` independent optimizations from random initial points and returns
the best result.

`maxtime` and `maxiters` are Optimization.jl common solver options forwarded to
every `Optimization.solve` call. `abstol`, `reltol`, and `callback` are common
options forwarded only when set (otherwise each solver uses its own default).
`solver_kwargs` is a `NamedTuple` of solver-specific options forwarded verbatim
to `solve` (e.g. `(; popsize=200)` for a CMA-ES solver that supports it). A key
present in both a named common option and `solver_kwargs` takes the
`solver_kwargs` value. `lb`/`ub` are the log-space bounds passed to the
`OptimizationProblem`.

Rescaling is driven by `fp.scale_k_to_kcat`: when it is a `Real`, the returned
parameters are rescaled so that `_kcat_forward(mechanism, params) ≈
fp.scale_k_to_kcat`; when it is `nothing`, the raw (unrescaled) parameters are
returned (the data fixes the absolute scale).

Returns a NamedTuple `(params, loss, retcode)` where:
- `params`: fitted rate constants as a NamedTuple
- `loss`: the best loss value achieved
- `retcode`: the `Symbol` form of the best restart's `sol.retcode`. Only
  `:Success` indicates the optimizer converged on its own criteria; any other
  value (e.g. `:MaxTime` — hit the time budget; `:Failure`; `:Default` — ran but
  did not flag success) means the fit should be treated as un-converged (check
  `retcode !== :Success`).
"""
function fit_rate_equation(fp::FittingProblem, optimizer;
    n_restarts::Int=20,
    maxtime::Real=60.0,
    maxiters::Integer=10_000_000,
    abstol::Union{Real,Nothing}=nothing,
    reltol::Union{Real,Nothing}=nothing,
    callback=nothing,
    lb=fill(-15.0, length(fitted_params(fp.mechanism))),
    ub=fill(15.0, length(fitted_params(fp.mechanism))),
    solver_kwargs=(;),
)
    obj = Optimization.OptimizationFunction((x, p) -> loss!(x, p))
    np = length(fitted_params(fp.mechanism))

    # Common solver options: maxtime/maxiters always forwarded; the optional
    # ones only when set, so each solver keeps its own default otherwise.
    # solver_kwargs is forwarded verbatim and wins on key conflicts.
    common = (; maxtime, maxiters)
    abstol   === nothing || (common = (; common..., abstol))
    reltol   === nothing || (common = (; common..., reltol))
    callback === nothing || (common = (; common..., callback))
    solve_kwargs = merge(common, solver_kwargs)

    best_x = zeros(np)
    best_loss = Inf
    # Sentinel for "no restart produced a finite objective" (loss stays Inf, so
    # the fit is dropped by the beam's non-finite filter). Deliberately NOT a
    # real SciMLBase ReturnCode name — `Symbol(ReturnCode.Default) === :Default`,
    # so `:NoFiniteLoss` stays distinct from a genuine `:Default` solver return
    # (which only fires when a restart achieves a finite objective).
    best_retcode = :NoFiniteLoss

    for _ in 1:n_restarts
        x0 = clamp.(randn(np) .* 2.0, lb, ub)
        prob = Optimization.OptimizationProblem(obj, x0, fp; lb=lb, ub=ub)
        sol = Optimization.solve(prob, optimizer; solve_kwargs...)
        if sol.objective < best_loss
            best_loss = sol.objective
            best_x .= sol.u
            best_retcode = Symbol(sol.retcode)
        end
    end

    pnames = fitted_params(fp.mechanism)
    result_params = NamedTuple{pnames}(ntuple(i -> exp(best_x[i]), Val(length(pnames))))
    if fp.scale_k_to_kcat !== nothing
        full = merge(result_params, (Keq = fp.Keq, E_total = 1.0))
        rp = rescale_parameter_values(
            fp.mechanism, full; scale_k_to_kcat=fp.scale_k_to_kcat,
        )
        result_params = NamedTuple{pnames}(
            ntuple(i -> rp[pnames[i]], Val(length(pnames))),
        )
    end
    return (params = result_params, loss = best_loss, retcode = best_retcode)
end
```

- [ ] **Step 5: Rewrite `identify_rate_equation` docstring + signature** (`src/identify_rate_equation.jl`)

(a) Replace the four docstring bullets (lines 112-118):
```julia
- `maxtime::Real = 60.0`: max time per fit (seconds)
- `maxiters::Int = 10_000_000`: max iterations per
  optimizer run (forwarded to `Optimization.solve`)
- `popsize::Int = 200`: population size for optimizer
  (forwarded to `Optimization.solve`)
- `verbose::Int = -9`: optimizer verbosity
  (forwarded to `Optimization.solve`)
```
with:
```julia
- `maxtime::Real = 60.0`: max time per fit (seconds; common solver
  option, forwarded to `Optimization.solve`)
- `maxiters::Integer = 10_000_000`: max iterations per optimizer run
  (common solver option, forwarded to `Optimization.solve`)
- `abstol`/`reltol`/`callback = nothing`: Optimization.jl common solver
  options, forwarded to `Optimization.solve` only when set
- `solver_kwargs::NamedTuple = (;)`: solver-specific options forwarded
  verbatim to `Optimization.solve` (e.g. `(; popsize=200)` for a CMA-ES
  solver that supports it); the caller matches its contents to `optimizer`
```

(b) Delete the trailing "Extra kwargs" docstring bullet (lines 134-135):
```julia
- Extra kwargs are forwarded to `fit_rate_equation`
  and then to `Optimization.solve`.
```

(c) Replace the signature + `fitting_kwargs` block (lines 179-197):
```julia
    maxtime::Real = 60.0,
    maxiters::Int = 10_000_000,
    popsize::Int = 200,
    verbose::Int = -9,
    # Model selection
    n_cv_candidates::Int = 5,
    se_threshold::Float64 = 1.0,
    perm_p_threshold::Float64 = 0.16,
    # Output & parallelism
    save_dir::String = _default_save_dir(),
    show_progress::Bool = true,
    pmap_function::Function = pmap,
    # Extra fitting/optimizer kwargs
    optim_kwargs...
)
    fitting_kwargs = (;
        n_restarts, maxtime,
        maxiters, popsize, verbose,
        optim_kwargs...)
```
with:
```julia
    maxtime::Real = 60.0,
    maxiters::Integer = 10_000_000,
    abstol::Union{Real,Nothing} = nothing,
    reltol::Union{Real,Nothing} = nothing,
    callback = nothing,
    solver_kwargs = (;),
    # Model selection
    n_cv_candidates::Int = 5,
    se_threshold::Float64 = 1.0,
    perm_p_threshold::Float64 = 0.16,
    # Output & parallelism
    save_dir::String = _default_save_dir(),
    show_progress::Bool = true,
    pmap_function::Function = pmap,
)
    fitting_kwargs = (;
        n_restarts, maxtime, maxiters,
        abstol, reltol, callback, solver_kwargs)
```

- [ ] **Step 6: Run the full suite — verify GREEN**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS. The regression test, the migrated pipeline/recovery/loocv tests, and the reworked all-base-fail + `_process_batch` failure tests all pass on CMAEvolutionStrategy. (The reworked `solver_kwargs=(; popsize=200)` paths exercise the verbatim-forwarding rejection.)

- [ ] **Step 7: Commit**

```bash
git add src/fitting.jl src/identify_rate_equation.jl test/test_identify_rate_equation.jl
git commit -m "Forward solver options via named commons + solver_kwargs bag

Drop the force-injected popsize/verbose and the kwargs catch-alls;
fit/identify now run on any Optimization.jl solver. Migrate identify
tests from PyCMA to CMAEvolutionStrategy."
```

---

## Task 3: Explicit forwarding / merge / clean-break tests

The behavior is already correct after Task 2; these lock the contract with focused assertions (some overlap the migrated tests, but state the contract directly).

**Files:**
- Modify: `test/test_fitting.jl`
- Modify: `test/test_identify_rate_equation.jl`

- [ ] **Step 1: Add the fit-level testset** to `test/test_fitting.jl`

Place it immediately after the `@testset "scale_k_to_kcat normalization"` testset (after line 362), at the same nesting level, so the file's `uni_uni` fixture and `make_synthetic_data` helper are in scope:

```julia
    # ── Test: solver-option forwarding (named commons + solver_kwargs) ──
    @testset "solver kwarg forwarding" begin
        using OptimizationCMAEvolutionStrategy
        Keq_val = 2.0
        true_params = (kon_S_E = 10.0, kon_P_ES = 5.0, koff_P_ES = 1.0,
            Keq = Keq_val, E_total = 1.0)
        concs_list = [
            (S = 0.5, P = 0.1), (S = 1.0, P = 0.1), (S = 2.0, P = 0.1),
            (S = 5.0, P = 0.1), (S = 10.0, P = 0.1),
        ]
        data = make_synthetic_data(uni_uni, true_params, concs_list)
        fp = FittingProblem(uni_uni, data; Keq=Keq_val)

        # Default (empty) solver_kwargs runs on a solver that rejects popsize —
        # no solver-specific option is force-injected.
        res = fit_rate_equation(fp, CMAEvolutionStrategyOpt();
            n_restarts=1, maxtime=1.0)
        @test haskey(res, :params)
        @test res.retcode isa Symbol

        # solver_kwargs is forwarded verbatim: an option this solver rejects
        # surfaces as an error (proves the bag reaches `solve`).
        @test_throws Exception fit_rate_equation(
            fp, CMAEvolutionStrategyOpt();
            n_restarts=1, maxtime=1.0, solver_kwargs=(; popsize=200))

        # Merge semantics: the same key in both a named common option and
        # solver_kwargs does NOT raise a duplicate-keyword error (a naive
        # double-splat would); solver_kwargs wins.
        res2 = fit_rate_equation(
            fp, CMAEvolutionStrategyOpt();
            n_restarts=1, maxtime=60.0, solver_kwargs=(; maxtime=1.0))
        @test haskey(res2, :params)

        # Clean break: popsize/verbose are no longer accepted named kwargs.
        @test_throws Exception fit_rate_equation(
            fp, CMAEvolutionStrategyOpt(); n_restarts=1, maxtime=1.0, popsize=200)
        @test_throws Exception fit_rate_equation(
            fp, CMAEvolutionStrategyOpt(); n_restarts=1, maxtime=1.0, verbose=-9)
    end
```

- [ ] **Step 2: Add the identify-level clean-break testset** to `test/test_identify_rate_equation.jl`

Add as a standalone testset at the end of the file (top level, alongside the other standalone testsets):

```julia
@testset "removed kwargs error at the identify boundary" begin
    # popsize/verbose are no longer named kwargs and there is no catch-all,
    # so they are rejected immediately at the call boundary (before any
    # fitting or CSV write) — distinct from a solver-rejected solver_kwargs
    # option, which fails inside fitting.
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products: P[C]
    end
    data = (group = ["G1", "G1", "G2", "G2"],
            Rate = [0.5, 0.8, 1.0, 1.1],
            S = [1.0, 2.0, 3.0, 4.0],
            P = [0.1, 0.2, 0.3, 0.4])
    prob = IdentifyRateEquationProblem(rxn, data; Keq=10.0)
    @test_throws Exception identify_rate_equation(
        prob; popsize=200, optimizer=CMAEvolutionStrategyOpt(),
        n_restarts=1, maxtime=1.0, save_dir=mktempdir())
    @test_throws Exception identify_rate_equation(
        prob; verbose=-9, optimizer=CMAEvolutionStrategyOpt(),
        n_restarts=1, maxtime=1.0, save_dir=mktempdir())
end
```

- [ ] **Step 3: Run the full suite — verify GREEN**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS (new testsets included).

- [ ] **Step 4: Commit**

```bash
git add test/test_fitting.jl test/test_identify_rate_equation.jl
git commit -m "Test solver_kwargs forwarding, merge override, and clean break"
```

---

## Task 4: Migrate README to CMAEvolutionStrategy + document solver_kwargs

`test/test_readme_runs.jl` executes the README's non-skipped ```julia blocks; the fit block uses `PyCMAOpt()`/`popsize=50`, so the README must migrate before the PyCMA dep is removed (Task 5).

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Switch the import** (line 81)

```julia
using OptimizationPyCMA, Random
```
→
```julia
using OptimizationCMAEvolutionStrategy, Random
```

- [ ] **Step 2: Update the fit prose** (lines 118-119)

```
The fit runs `fit_rate_equation` on a `FittingProblem`, using the PyCMA
optimizer (multi-start CMA-ES) recommended for rate-equation fitting.
```
→
```
The fit runs `fit_rate_equation` on a `FittingProblem`, using the
CMAEvolutionStrategy optimizer (multi-start CMA-ES) recommended for
rate-equation fitting.
```

- [ ] **Step 3: Update the fit code block** (lines 126-127)

```julia
result = fit_rate_equation(fp, PyCMAOpt();
    n_restarts=3, maxtime=5.0, popsize=50)
```
→
```julia
result = fit_rate_equation(fp, CMAEvolutionStrategyOpt();
    n_restarts=3, maxtime=5.0)
```

- [ ] **Step 4: Document `solver_kwargs`** — add this paragraph immediately after the closing ``` of the fit code block (after line 132):

```markdown
Solver-specific options are passed through `solver_kwargs`, a `NamedTuple`
forwarded verbatim to `Optimization.solve`; the Optimization.jl common
options `maxtime`, `maxiters`, `abstol`, `reltol`, and `callback` are named
keyword arguments. Match `solver_kwargs` to your chosen optimizer — for
example `CMAEvolutionStrategyOpt` does not accept `popsize`.
```

- [ ] **Step 5: Update the identify code block** (line 170)

```julia
    optimizer=PyCMAOpt(),
```
→
```julia
    optimizer=CMAEvolutionStrategyOpt(),
```

- [ ] **Step 6: Run the full suite — verify GREEN**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS — `test/test_readme_runs.jl` runs the migrated fit block on `CMAEvolutionStrategyOpt` (the identify block stays `# README-SKIP-IN-TEST`). No remaining `PyCMAOpt`/`popsize` in any non-skipped block.

- [ ] **Step 7: Commit**

```bash
git add README.md
git commit -m "README: use CMAEvolutionStrategy and document solver_kwargs"
```

---

## Task 5: Remove OptimizationPyCMA from test dependencies

No test or README references PyCMA after Tasks 2-4. Remove it from all three `Project.toml` locations (Aqua checks that `[extras]`, `[compat]`, and `[targets]` stay consistent).

**Files:**
- Modify: `Project.toml`

- [ ] **Step 1: Remove the `[compat]` entry**

Delete the line:
```toml
OptimizationPyCMA = "1.0.0"
```

- [ ] **Step 2: Remove the `[extras]` entry**

Delete the line:
```toml
OptimizationPyCMA = "fb0822aa-1fe5-41d8-99a6-e7bf6c238d3b"
```

- [ ] **Step 3: Remove from the test target**

Replace the `[targets]` line:
```toml
test = ["Test", "OrdinaryDiffEqFIRK", "OptimizationBBO", "OptimizationCMAEvolutionStrategy", "OptimizationPyCMA", "Aqua", "JET"]
```
with:
```toml
test = ["Test", "OrdinaryDiffEqFIRK", "OptimizationBBO", "OptimizationCMAEvolutionStrategy", "Aqua", "JET"]
```

- [ ] **Step 4: Resolve and run the full suite — verify GREEN**

Run: `julia --project -e 'using Pkg; Pkg.test()'`
Expected: PASS. Aqua's project/compat/extras consistency checks pass (PyCMA gone from all three locations); no test imports `OptimizationPyCMA`.

- [ ] **Step 5: Commit**

```bash
git add Project.toml
git commit -m "Remove OptimizationPyCMA test dependency"
```

---

## Self-Review

**Spec coverage:**
- Taxonomy (named commons + `solver_kwargs`) → Task 2 Steps 4-5. ✓
- Merge semantics → Task 2 Step 4 (`merge(common, solver_kwargs)`); tested Task 3 Step 1. ✓
- Thread through identify (`fitting_kwargs`) → Task 2 Step 5(c). ✓
- Clean break (drop popsize/verbose + catch-alls) → Task 2 Steps 4-5; tested Task 3 Steps 1-2. ✓
- Recommended/test solver = CMAEvolutionStrategy only; PyCMA dropped → Tasks 1, 2(1), 4, 5. ✓
- Test plan #1 regression → Task 2 Step 2. #2 BBO second solver → existing `test_fitting.jl` BBO tests (unchanged, now also get `maxiters`). #3 verbatim-forward rejection → Task 2 Step 1(g/h) + Task 3 Step 1. #4 merge → Task 3 Step 1. #5 clean break → Task 3 Steps 1-2. #6 call-site migration → Task 2 Step 1. #7 README → Task 4. ✓
- `Project.toml` extras/compat/targets (test-target pattern, not `test/Project.toml`) → Tasks 1, 5. ✓
- Out-of-scope upstream wrapper PR → not in plan (correct). ✓

**Placeholder scan:** none — every step gives exact code/old-new strings and a run command with expected result.

**Type/name consistency:** `solver_kwargs` (NamedTuple), `cmaes_opt`/`CMAEvolutionStrategyOpt()`, `maxiters::Integer`, `abstol`/`reltol`/`callback` named identically in `fit_rate_equation`, `identify_rate_equation`, and tests. `merge(common, solver_kwargs)` matches the spec. `make_synthetic_data`/`uni_uni` reused from existing `test/test_fitting.jl` scope.
