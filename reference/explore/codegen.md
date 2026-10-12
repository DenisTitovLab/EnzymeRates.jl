# What compiles per mechanism, and what costs memory in a long identification run

## Bottom line

1. **The hypothesis does not hold for the current pipeline.** I ran each mechanism through what `identify_rate_equation` does to it: PASS1 compile and render, PASS2 fit, and one LOOCV fold. The derivation side costs about **0.08 s of the 1.21 s of compile** per mechanism and about **0.12 MB of the 2.25 MB of permanent memory**. That is every `@generated` body plus its forwarders. The fitting stack costs **85–90%**, because it is compiled again for every mechanism type.
2. **Cause.** `FittingProblem{M,D}` carries the singleton mechanism type `M`. `OptimizationProblem(obj, x0, fp; …)` (`src/fitting.jl:251`) puts `fp` into the problem type. From there everything downstream specializes on `M`:
   - `solve`, `OptimizationCache`, `OptimizationState` and `OptimizationSolution`;
   - the `_cb` and `_loss` closures in OptimizationCMAEvolutionStrategy (`…/OptimizationCMAEvolutionStrategy/MySJY/src/OptimizationCMAEvolutionStrategy.jl:83,101,111`);
   - `CMAEvolutionStrategy.minimize`, `run!`, `evaluate`, `log!`, `BasicLogger` and `Optimizer`.

   That adds up to about 370 new method instances per mechanism. LOOCV pays the same again: it builds the fold data as views (`src/identify_rate_equation.jl:892`), which is a different `D`.
3. **A type barrier at the optimizer boundary removes most of it, and I measured this.** The prototype, `fit_erased` in `expB.jl`, hands the optimizer one concrete problem type for every mechanism: `ErasedObjective(f::Function)`. Per mechanism it cut:
   - compile from 1.21 s to 0.32 s (−73%);
   - permanent memory from 2.25 MB to 0.44 MB (−80%);
   - new method instances from 744 to 94.

   The runtime cost is +1 to 60 ns on a 0.3–2 µs loss call (≤3%) and a 16-byte allocation per call (the boxed return value). End-to-end fit time did not change.
4. **Growth is linear, with no saturation:** 2.16 MB per mechanism with the package fitter, 0.41 MB with the prototype, 0.10 MB for derivation alone. The first and second halves of each run grow at the same rate.
5. **Derivation compile scales with equation size**, by about 66 ms and about 21 KB of native code per 1000 optimized IR statements. The fitting stack does not depend on size: about 0.55 s and 305 KB of native code per data type.
6. **The derivation's own Julia work is small.** It takes 1–4 ms per mechanism, which is 4% (at most 9%) of the first `rate_equation` call. The rest is inference and LLVM on the generated body. Moving derivation code out of `@generated` does not reduce codegen; only a smaller generated body would.
7. **Two side findings.**
   - (a) `string(typeof(em))`, used for the CSV `mechanism_type` column (`src/identify_rate_equation.jl:466`), compiles `show` for every new shape of the Sig tuple. That is 72 ms per mechanism, the largest PASS1 cost.
   - (b) `Base.gc_live_bytes()` rises about 36–49 KB per `_dependent_param_exprs` call, but no real memory is behind it. The cause is the `Rational{BigInt}` work. After about 40k calls the phantom total nears the 2.5 GB heap hint, and from then on every collection is a full sweep: GC time triples and wall time rises 75%.

## Setup

- Julia 1.12.5 (aarch64), one thread. EnzymeRates 0.9.0 at `277a785`, OptimizationCMAEvolutionStrategy 0.3.2, CMAEvolutionStrategy 0.2.6, Optimization 4.8.0, OptimizationBase 2.14.0, SciMLBase 2.155.1. Single process, no `pmap` workers; the per-process numbers equal per-worker numbers.
- **Population:** 74 distinct compiled types from the reaction `A[C], B[N] → P[C], Q[N]` with allosteric regulator `R`.
  - 48 base-tier mechanisms, up to 4 per King–Altman complexity class (1 to 245).
  - 8 plain children, 12 allosteric children, 6 allosteric grandchildren.
  - 5–9 fitted parameters, 6–191 denominator terms, 47–2443 optimized IR statements in `rate_equation`.
  - Three warm-up mechanisms outside the population ran first. The population then ran in a fixed shuffled order.
- **Measurement around each stage:** two full GCs, then `malloc_trim`, then these readings:
  - cumulative compile time;
  - pool live bytes (`jl_gc_pool_live_bytes`);
  - native code (`Base.jit_total_bytes`);
  - malloc in use (`mallinfo2`);
  - RSS;
  - a method-instance census over `Core.methodtable`;
  - per-entry compile attribution from `--trace-compile=stderr --trace-compile-timing`, delimited by markers.
- **Memory metric.** I use pool + native code + malloc as permanent memory, not `gc_live_bytes`, because the latter drifts (finding 7b).
- **Repeatability.** A first run of all three experiments (kept in `run1/`) reproduced the compile totals within 3% and the native-code and malloc numbers exactly.

## 1. What specializes per mechanism

| Function family | New method instances per mechanism | Evidence |
|---|---|---|
| `rate_equation`: generated body plus forwarder (`src/rate_eq_derivation.jl:413,427`) | 2 per (concs, params) NamedTuple pair; 4 in the pipeline | census; trace 49 ms/mech; 8–221 KB native code |
| `fitted_params` (`:47`) | 1. Returns a constant, so no native code | 1 ms, ~1 KB |
| `parameters` (`:38,41`); `metabolites` (`src/types.jl:1322,1336`) | 2 and ~1 | trivial; identify does not call `parameters` |
| `_kcat_forward` (`:593,639`), reached via `rescale_parameter_values` | 1–2 | 7.5 ms, ~3 KB native code |
| `rate_equation_string` forwarder | 2, plus about 10 inference-only instances of the `@nospecialize` helpers (`_dependent_param_exprs`, `_num_den_exprs`, `_concrete`, `Mechanism`, `_equation_text`, …), which get no native code | 9 ms, 20 KB pool |
| Per Sig shape, not strictly per mechanism: `show` on the Sig (`string(typeof(em))`); `collect`/`setindex_widen_up_to` in `_sig_of` | about 26 Base instances; 31 of 74 mechanisms hit a new shape, at 70–330 ms each | `show` 72 ms/mech; `_sig_of` 5 ms/mech |
| `_poly_to_expr` sort key `Tuple((string(s), e) …)` (`src/sym_poly_for_rate_eq_derivation.jl:109`): `isless`/`Order.lt` once per pair of tuple lengths | ~5 early in a run, then saturates | 11 ms/mech here, about 0.8 s per process |
| EnzymeRates fitting entry points: `FittingProblem` and its keyword constructor, `fit_rate_equation` and its body, `loss!`, `_cv_fold_loss`, `rescale_parameter_values` body, `_rows_to_dataframe` closure | ~21 per data-container type | census |
| Optimization, SciMLBase, OptimizationBase and OptimizationCMAEvolutionStrategy | 76 | census (B2) |
| CMAEvolutionStrategy: `minimize`, `Optimizer`, `BasicLogger`, `run!`, `evaluate`, `log!`, `noise_handling!`, `xbest` | 135 | trace: `evaluate` 100, `minimize` 68, `log!` 30 ms/mech |
| Base and Core on the per-mechanism problem and NamedTuple types: `merge`, `getproperty`, `NamedTuple`, `keys`, `iterate`, `Core.Typeof` | ~135 | census |
| Julia 1.12's own compiler closures (`Base.Compiler.var"#inferiterate_2arg#…"`) specialized on NamedTuple-merge types that carry each mechanism's parameter names | ~1.7 entries | 124 ms/mech; 274 trace entries with the package fitter, 6 with the prototype |

The singleton type itself is cheap: `summarysize(typeof(em))` is 0.65–1.2 KB, and `compile_mechanism` adds 9 KB of pool on average.

## 2. Measurements

**Per-mechanism totals** (mean / median / max; memory is pool + native code + malloc):

| Run | Compile (s) | Pool (KB) | Native code (KB) | Malloc (KB) | Permanent (KB) | New MIs |
|---|---|---|---|---|---|---|
| A: derivation only, all entry points | 0.079 / 0.060 / 0.501 | 50 / 31 / 443 | 28 / 19 / 272 | 43 / 30 / 546 | **121** / 78 / 1084 | 40 |
| B package fitter: PASS1 + PASS2 | 0.679 / 0.661 / 1.348 | 507 | 340 | 351 | **1197** / 1151 / 2531 | 401 |
| B package fitter: plus 1 LOOCV fold | 1.210 / 1.194 / 1.897 | 967 | 619 | 662 | **2248** / 2200 / 3758 | 744 |
| B prototype: PASS1 + PASS2 | 0.202 / 0.168 / 0.790 | 79 | 90 | 110 | **280** / 243 / 1478 | 67 |
| B prototype: plus 1 LOOCV fold | 0.324 / 0.283 / 0.910 | 123 | 154 | 162 | **439** / 412 / 1647 | 94 |

- **Ratios.** The fit alone costs 8.6× A in compile and 9.9× in memory. With the LOOCV fold it is 15.4× and 18.5×. The prototype over the package fitter is 0.27 in compile and 0.20 in memory.
- **Whole-run RSS growth over 74 mechanisms.** A +5.5 MB, package fitter +98.8 MB, prototype +18.5 MB. The GC reuses free heap pages early in a run, so RSS lags the permanent-memory sum.
- **Allosteric (n=18) against non-allosteric (n=56).** Same picture: A 0.080 vs 0.078 s; package fitter with fold 1.19 vs 1.22 s and 2.35 vs 2.22 MB; prototype 0.28 vs 0.34 s and 0.44 vs 0.44 MB.

**Per stage** (means):

| Stage | Compile (s) | Permanent (KB) | Native code (KB) | MIs |
|---|---|---|---|---|
| A1 `rate_equation` (with `parameters`, `metabolites`) | 0.053 | 91 | 20 | 15 |
| A2 `fitted_params` | 0.001 | 1 | 0 | 1 |
| A3 `_kcat_forward` | 0.008 | ~8 (median) | 3 | 1 |
| A4 `rate_equation_string` | 0.011 | 25 | 2 | 15 |
| B1 PASS1 (same in both fitters) | 0.093 | 105 | 34 | 31 |
| B2 fit, package fitter / prototype | 0.585 / 0.109 | 1092 / 177 | 305 / 56 | 370 / 36 |
| B3 LOOCV fold, package fitter / prototype | 0.531 / 0.122 | 1051 / 160 | 280 / 64 | 344 / 27 |

**Generator Julia work against compile time** (mean / median / max):

| Item | Time |
|---|---|
| `_dependent_param_exprs` | 0.22 / 0.20 / 0.44 ms |
| `_num_den_exprs` | 1.03 / 0.59 / 3.63 ms |
| `rate_equation_string`, warm | 1.6 / 1.2 / 5.2 ms |
| A1 compile | 53 / 29 / 400 ms |

The derivation work is 4% of the A1 compile (median 3%, max 9%).

**Scaling of the `rate_equation` compile with size:**

| IR statements | Denominator terms | Compile (ms) | Native code (KB) |
|---|---|---|---|
| 0–100 | 6–12 | 11 | 8 |
| 700–1300 | about 109 | 80 | 29 |
| 2000–3000 | 191 | 186 | 72 |

- The fit is 3.4 ms + 65.9 ms per 1000 statements (r = 0.71).
- The package fitter's fit compile is 0.547 s + 0.055 s per 1000 statements. Its native code does not track size (r = 0.03).
- At the `eq_complexity_filter` cap of 337 terms (about 4300 statements) I extrapolate about 0.3 s and 100–200 KB of native code. That is still below the fitting stack.

**Runtime cost of the barrier** (from `bench_erased.jl`):

| Measure | Direct | Erased |
|---|---|---|
| Loss call, 5 mechanisms (ns) | 303 / 276 / 1576 / 1934 / 1951 | 307 / 277 / 1617 / 1994 / 1979 |
| Allocation per call | 0 B | 16 B |
| CMA-ES solve, 300 iterations | 2–8 ms | same 2–8 ms |

**The `gc_live_bytes` drift** (`probe_leak2.jl`, `probe_leak3.jl`, `probe_drift_gc.jl`):

- 6000 calls of `_dependent_param_exprs` on one mechanism: `gc_live` +218 MB, pool live +0.03 MB, RSS +32 MB. The RSS growth is a one-time allocator step, not linear.
- 100k dead BigInts: `gc_live` +3.6 MB, with RSS and malloc flat.
- After about 40k derivation calls, full sweeps went from 4 of 18 collections to all of them. GC time per 10k calls rose from 1.0 to 3.2 s and wall time from 2.7 to 4.9 s, while RSS stayed at 1.0–1.15 GB.

I only tested this on 1.12.5.

## 3. Interpretation

- **Memory growth comes from the fitting stack being re-specialized per mechanism and per data-container type.** It is not derivation codegen, and it is not type proliferation of the singleton.
- **NamedTuple types** keyed on each mechanism's parameter names cost something real but secondary:
  - with the package fitter, about 135 Base/Core instances per fit and the 124 ms of compiler-closure compiles;
  - with the barrier in place, about 22 Base/Core instances and about 15 ms of `merge`, `iterate`, the `rescale_parameter_values` body and the rows closure.
- **The derivation side grows in proportion to equation size but stays small in this population.** Per mechanism: the kernel's codegen (11–186 ms), `show` on the Sig (72 ms), the `Order.lt` warm-up (saturates), and about 1 ms each for the small generated helpers.
- **Projection for one process.** Each fit with LOOCV costs about 2.2 MB now, so 5,000 fitted equations reach about 11 GB, above this 7.7 GB machine. With the barrier the same 5,000 cost about 2.2 GB. Mechanisms that only pass PASS1 (skipped or inherited) cost about 0.1 MB each.

## 4. Recommendations, ranked, with expected savings

1. **Add a type barrier at the optimizer boundary** (`src/fitting.jl:216-267`). Measured saving: −0.89 s and −1.8 MB per mechanism.
   - The prototype boxes the return value. A `Ref{Float64}` output argument or FunctionWrappers.jl would remove the 16-byte allocation per call.
   - Keep only `loss!` and the rate kernel specialized.
2. **Give LOOCV folds the same container type as PASS2.** Copy each fold's rows into plain `Vector{Float64}` columns, or keep one data matrix, instead of `view` (`src/identify_rate_equation.jl:892`).
   - With the package fitter this avoids a second full stack (−0.53 s, −1.05 MB per candidate).
   - With the barrier it still saves the 62 ms `fold` closure plus about 50 ms of `loss!` for the view type.
3. **Build `FittingProblem` without specializing on the mechanism.** Use `@nospecialize` and a data layout that does not depend on the mechanism. Estimated −30 ms per mechanism per data type.
4. **Stop calling `string(typeof(em))`** (`src/identify_rate_equation.jl:466`). Print the concrete `Mechanism` with a printer that does not specialize, which also serves goal 13. Estimated −72 ms, about −26 method instances and most of PASS1's 105 KB.
5. **Give `_poly_to_expr` a sort key of one fixed type** (`src/sym_poly_for_rate_eq_derivation.jl:104-111`). Estimated −0.8 s per process.
6. **Keep fitted values as `Vector{Float64}` plus a names vector** and pass kernel arguments by position. Estimated −20 to 40 ms and about 20 method instances per mechanism. This also lets `metabolites` and `fitted_params` stop being `@generated`, but that part saves only about 1 ms and about 1 KB each.
7. **Do not make de-generating `rate_equation` a memory priority.** Moving everything except the numeric kernel out of `@generated` saves about 10–20 ms per mechanism. Moving the derivation work out of the generator saves less than 4 ms. Only shrinking the body cuts kernel codegen, for example by sharing spanning-tree subproducts or factoring by metabolite. The realistic floor with items 1–6 done is about 0.06–0.1 s and 0.1–0.2 MB per mechanism.
8. **Remove the drift** by replacing `Rational{BigInt}` in the constraint solves with checked `Int` or modular arithmetic. The places are:
   - `src/thermodynamic_constr_for_rate_eq_derivation.jl:123,148,370,488,531`;
   - `src/rate_eq_derivation.jl:323,971`;
   - enumeration's `_partition_independent_count` and `_onlya_haldane_violation`.

   This keeps the GC's accounting honest in long runs.

## Files

All files are in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/codegen/`. No repository files were changed, and no Julia process is left running.

Scripts:
- `make_population.jl` → `pop.jls`
- `common.jl` (measurement tools)
- `expA.jl`
- `expB.jl` (argument `orig` or `erased`; contains `fit_erased` and `cv_fold_erased`)
- `bench_erased.jl`
- `probe_mi.jl`, `probe_mi2.jl`, `probe_mi3.jl`, `probe_mi4.jl`, `probe_mi5.jl`
- `probe_leak.jl`, `probe_leak2.jl`, `probe_leak3.jl`
- `probe_drift_gc.jl`
- `analyze.py`, `analyze_trace.py`, `scaling.py`

Data:
- `expA.tsv`, `expA_bymeth.tsv`, `expA_sizes.tsv`, `expA.trace`
- `expB_{orig,erased}.tsv`, `expB_{orig,erased}_bymeth.tsv`, `expB_{orig,erased}.trace`, `expB_{orig,erased}.out`
- `run1/` (the first, reproducing run)