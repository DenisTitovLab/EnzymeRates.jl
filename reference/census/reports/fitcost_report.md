# Fit-cost calibration (census track `fitcost`)

## Bottom line

1. **Fits are cheap and never come near the time limits.** I fitted each candidate with identify's own defaults. One candidate's fit (20 restarts, cold) took 0.05–4.2 s, median 0.20–0.23 s, on 48-point relative data. The longest single restart took 0.51 s against a maxtime of 600 s. Of 13,280 restarts, 13,277 stopped on CMA-ES's ftol rule (retcode Success) and 3 on stagnation (Stalled). None stopped at maxtime or maxiters.
2. **Rank-deficient candidates do not cost more.** At equal parameter count they need about 25% fewer generations (log coefficient −0.29 ± 0.06, pooled over both datasets). All 3,880 of their restarts ended on ftol, and the longest took 0.23 s.
3. **The parameter count drives the cost; law size matters little at 48 points.** Generations per restart grow as np^2.3, with a ×1.9 scatter between candidates that np, law size and rank do not explain. One generation costs the CMA-ES overhead plus λ objective calls, and the call cost is linear in law size. At 48 points the objective is only a median 40% of a generation.
4. **Fixed costs per candidate are large.** They come to about 68 ms + 71 µs × law size: identify's pass 1, building the FittingProblem, the first `loss!` compile and compilation inside `fit_rate_equation`. That is 23–29% of the summed time and 35–39% of a median candidate's time.
5. **Weighting rule for the census.** Rank-deficient candidates' share of predicted time lies between their count share and their D-terms share. In bi_bi round 1 they are 22% of the count, 13% of the predicted time and 7% of the D terms.

## What identify does per fit (from the code)

- **identify defaults.** `identify_rate_equation` requires an optimizer; its docstring recommends `CMAEvolutionStrategyOpt()`. The other defaults are `n_restarts=20`, `maxtime=600.0`, `maxiters=10_000_000`, `abstol`/`reltol`/`callback = nothing` and `solver_kwargs=(;)`.
- **The restart loop.** `fit_rate_equation` (src/fitting.jl) draws `x0 = clamp.(2·randn(np), ±15)` and uses bounds of ±15 in log space. It passes `maxtime` and `maxiters` to every `Optimization.solve` call, so **maxtime applies per restart**: a fit can take up to 20 × 600 s, although identify's docstring says "max time per fit". It keeps the best restart, then rescales the result to kcat = 1.
- **The CMA-ES wrapper** (OptimizationCMAEvolutionStrategy 0.3.13, CMAEvolutionStrategy 0.2.6):
  - Settings: sigma0 = 0.1, population λ = 4 + ⌊3 ln n⌋, and `maxiters` counts generations.
  - Stop rules: ftol 1e-11 on the population's f-range gives Success; stagnation gives Stalled; xtol is off.
  - The wrapper reports `fevals` = generations; I counted evaluations as generations × λ.
- **Dedup, workers and LOOCV.** identify fits one representative per `eq_hash` and spreads the fits over workers with `pmap`. LOOCV then refits the top 5 per parameter count once per held-out group.

## Setup

- **Code and environment.** Package branch `onlya-guard-stranded-segment` (11cac53), Julia 1.12.5. The environment is a copy of `ci_repro/env` (Optimization 5.9.1). Julia time was about 24 min over 7 runs, one Julia process at a time.
- **Pools.** I built them with `census/small/lib.jl`, which follows identify's seeds, base tier, expansion, dedup, caps and `eq_hash`.

  | reaction | round 0 | round 1 | round 2 |
  |---|---|---|---|
  | bi_bi | 239 seeds fitted; 239 eligible | all 239 expanded, 2752 fitted; 2398 eligible | 30 parents sampled, 240 fitted; 207 eligible |
  | ter_bi | 4952 fitted (206 over the complexity cap); facts on 424; 424 eligible | 25 parents, 573 fitted; 356 eligible | 15 parents, 259 fitted; 146 eligible |

  "Eligible" means the guard passes (no W), the candidate is not a duplicate key (D) or the same law as its parent (Z), and the kcat precondition holds.
- **Sample: 332 candidates.** bi_bi has 168 (rounds 0/1/2: 17/118/33). ter_bi has 164 (rounds 0/1/2: 62/62/40).
  - Selection: half farthest-point sampling in (np, log law size), half random.
  - Rank-deficient (P) = package-law rank of ∂v/∂log θ below the fitted count. There are 97 P candidates: 50 in bi_bi (all allosteric n=1) and 47 in ter_bi (9 plain, 38 allosteric). Deficit 1 for 86, 2 for 10, 3 for 1.
  - np in the sample runs from 5 to 11.
  - Law size `ops` = call nodes in the package law's numerator, denominator and dependent constants. The slice's D terms ≈ (ops + 131)/12.7 (R² 0.985).
- **Data.** Relative mode, Keq 10, 48 points in 4 groups. Each series runs from 0.05 to 10; two groups carry one product at 1. Each group is scaled by 10^(0.3·N) and each rate by noise exp(0.05·N).
  - Main references: bi_bi #235 (np 6) and ter_bi #1448 (np 7).
  - Second dataset: bi_bi #49 and ter_bi #296, with different constants and noise.
  - The reference fits reach the noise floor (best loss 0.0028 and 0.0019 against σ² = 0.0025).
- **Budget.** Identify's defaults throughout; no scaling down was needed.

## Results (main dataset; the second dataset agrees)

| reaction | class | n | np median | ops median | generations/restart median | t_fit median / mean / max (s) | mean per-candidate total (s) |
|---|---|---|---|---|---|---|---|
| bi_bi | full rank | 118 | 6 | 149 | 332 | 0.16 / 0.28 / 2.21 | 0.34 |
| bi_bi | P | 50 | 6 | 48 | 209 | 0.09 / 0.10 / 0.25 | 0.15 |
| ter_bi | full rank | 117 | 7 | 214 | 833 | 0.46 / 0.65 / 4.19 | 0.74 |
| ter_bi | P | 47 | 8 | 103 | 1587 | 0.59 / 0.57 / 1.14 | 0.64 |

- **P against full rank at equal np** (median t_fit, s):

  | reaction, np | P | full rank |
  |---|---|---|
  | bi_bi, 6 | 0.08 (n=26) | 0.11 (n=61) |
  | bi_bi, 7 | 0.10 (n=13) | 0.41 (n=43) |
  | ter_bi, 8 | 0.66 (n=23) | 0.71 (n=30) |
  | ter_bi, 9 | 0.72 (n=9) | 0.91 (n=16) |

  P is cheaper or about equal in every bin on both datasets.
- **Why rank deficiency does not slow the fit** (my reading of the code, not tested). The ftol rule looks only at the spread of loss values in the population, and there is no x-tolerance, so flat directions never delay the stop. In relative mode every law already has one exactly flat direction, the overall rate scale.
- **Loss evaluation.** One `loss!` call takes 0.34–15.6 µs (median 0.65 µs), which is 7.3 ns + 0.049 ns × ops per data point. Inside CMA-ES each evaluation costs about 2.7× that once the optimizer's overhead is included.
- **Never-saturating candidates: none exist in the pools.**
  - My randomized numeric probe flagged one: bi_bi #1858 (allosteric n=1, np 7). The structural no-Vmax test found no never-saturating candidates in census/small's run of these bi_bi rounds.
  - Evaluated directly at 4 parameter draws, #1858 plateaus as the substrates go from 1 to 1e12, in both BigFloat and Float64.
  - Its fit was ordinary: 0.54 s, 20/20 Success.
- **Garbage collection depends on the process.** I fitted the same 80 candidates in a process that had first built the pools (about 8,000 compiled mechanism types, maxrss 2.0 GB) and in a fresh process. Generations and losses were identical (80/80). The loaded process spent 2.08× as long optimizing (67.9 s vs 32.7 s), with GC at 47% of fit time against 5%. Identify's workers accumulate compiled mechanisms over a search, so real runs should fall between the two.
- **More data makes law size matter.** I fitted 12 candidates on 48 and on 480 points.
  - Generations stayed about the same (median ratio ≈ 1.0, range 0.55–2.1).
  - `loss!` time grew 10×, time per generation 3.4–9×, and fit time 3.5–18×.
  - The per-generation formula below, with `t_loss` scaled by 10, predicts the per-generation time within about 5%.
- **A measurement artifact I found and fixed.** My first instrumented restart loop recompiled for every mechanism, adding about 75 ms per fit; the first run (`data_run1/`) and `fits_fresh.tsv` carry it. The corrected runs (`fresh2`, `alt`) time `EnzymeRates.fit_rate_equation` itself, cold, and run a non-specializing copy of the loop for per-restart statistics. That copy reproduces the package's best loss in 664/664 fits.

## Cost model (pooled: 662 fits = 331 candidates × 2 datasets)

**t_candidate = F + 20 · G · g**

| term | meaning | formula | fit quality |
|---|---|---|---|
| F | fixed cost per candidate | 68 ms + 71 µs · ops | |
| g | time per generation | 3.8 µs + 0.147 µs · np² + 0.94 · λ(np) · t_loss, with λ = 4 + ⌊3 ln np⌋ | R² 0.965 |
| t_loss | one `loss!` call, 48 points | 0.35 µs + 2.38 ns · ops | R² 0.90 |
| G | generations per restart | log G = 1.63 + 2.26(±0.19) · log np − 0.29(±0.06) · P + 0.65 · ter_bi | R² 0.40, residual sd 0.63 |

- G is the median; multiply by 1.22 for the mean.
- The `ter_bi` term (×1.9) depends on the dataset: it falls to ×1.4 on the second dataset.

**Scatter.** Actual/predicted has median 0.94, IQR 0.75–1.44, 5–95% range 0.53–2.54 and extremes 0.36–5.51. Summed predicted time is within 2% of actual on the main dataset and 4% on the second. Binned by predicted-time decile, the actual medians run 0.12 → 1.01 s while the predicted bins span 0.11 → 0.77+ s. The binned table is in `logs/pooled.log`.

**Log-linear alternative** (main dataset): t_candidate ≈ 3.7e-4 s · np^2.64 · ops^0.34 · 0.85^P (R² 0.58, scatter ×1.75). With slice D terms instead: 4.95e-4 · np^2.66 · Dterms^0.43 · 0.81^P. On the second dataset the np exponent drops to 1.84, while the size exponent (0.32) and the P factor (0.87) hold. Per-candidate generation counts correlate only 0.55 (log scale) between the two datasets.

**P's share of predicted serial fit time in my pools** (48 points):

| pool round | records | P share of count | of summed np | of D terms | of model time | model total (s) |
|---|---|---|---|---|---|---|
| bi_bi r0 | 239 | 0% | 0% | 0% | 0% | 45 |
| bi_bi r1 | 2398 | 22.2% | 21.8% | 7.2% | 13.3% | 599 |
| bi_bi r2 | 207 | 49.8% | 50.6% | 36.4% | 44.8% | 50 |
| ter_bi r0 | 424 | 1.4% | 1.7% | 0.5% | 1.3% | 209 |
| ter_bi r1 | 356 | 19.1% | 18.7% | 7.2% | 10.3% | 245 |
| ter_bi r2 | 146 | 39.7% | 40.8% | 22.1% | 31.6% | 105 |

At 480 points the bi_bi round-1 P share falls to about 8%.

## Examples

Notation: `⇌` is rapid equilibrium, `<-->` is steady state, `(…)` is a kinetic group, `::Tag` is an allosteric tag. Times are the main dataset.

**Full rank, small law**
- bi_bi #49, np 5, ops 37 — 0.05 s, 112 generations/restart:
  `(E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP); (E + B ⇌ EB, EA + B ⇌ EAB, EQ + B ⇌ EBQ); (E + P ⇌ EP, EA + P ⇌ EAP, EQ + P ⇌ EPQ); (E + Q ⇌ EQ, EB + Q ⇌ EBQ, EP + Q ⇌ EPQ); EAB <--> EPQ`
- bi_bi #1479, np 6, ops 54 — 0.22 s, 693 generations/restart:
  `E + A <--> EA ::OnlyA; E + Q ⇌ EQ ::EqualAI; EA + B ⇌ EPQ ::OnlyA; EQ + P <--> EPQ ::EqualAI | n=1`

**Full rank, large law**
- bi_bi #2871, np 7, ops 2815 — 2.21 s:
  `(E + A ⇌ EA, EB + A ⇌ EPQ, EQ + A ⇌ EAQ); (E + B <--> EB, EA + B <--> EPQ); (E + P <--> EP, EB + P <--> EBP, EQ + P <--> EPQ); (E + Q <--> EQ, EA + Q <--> EAQ, EP + Q <--> EPQ); EP + B <--> EBP`
- ter_bi #5621, np 8, ops 4951 — 4.19 s, the slowest fit; its first call alone took 0.25 s:
  `(E + A ⇌ EA, EB + A ⇌ EAB, EBD + A ⇌ EPQ, ED + A ⇌ EAD); (E + B ⇌ EB, EA + B ⇌ EAB, EAD + B ⇌ EPQ, ED + B ⇌ EBD); (E + D <--> ED, EA + D <--> EAD, EAB + D <--> EPQ, EB + D <--> EBD, EP + D <--> EDP); (E + P <--> EP, EQ + P <--> EPQ); (E + Q <--> EQ, EP + Q <--> EPQ); ED + P <--> EDP`

**Rank-deficient**
- bi_bi #248, np 6, rank 5 — 0.09 s:
  `E + A ⇌ EA ::OnlyA; E + Q ⇌ EQ ::OnlyA; (EA + B ⇌ EAB, EQ + B ⇌ EBQ) ::EqualAI; (EA + P ⇌ EAP, EQ + P ⇌ EPQ) ::EqualAI; EAB <--> EPQ ::OnlyA | n=1`
- ter_bi #4884 (plain), np 8, rank 7 — 0.29 s:
  `E + D ⇌ ED; E + P ⇌ EP; EAD <--> EQ_res_+A_+D_-Q; EB_res_+A_+D_-Q <--> EP; ED + A ⇌ EAD; E_res_+A_+D_-Q + B ⇌ EB_res_+A_+D_-Q; E_res_+A_+D_-Q + Q ⇌ EQ_res_+A_+D_-Q`
- ter_bi #5760 (plain), np 9, rank 8 — 1.11 s, the slowest P candidate:
  `E + D ⇌ ED; E + P ⇌ EP; EABD <--> EPQ; EAD + B <--> EABD; EBD + A <--> EABD; ED + A ⇌ EAD; ED + B ⇌ EBD; (ED + Q ⇌ EDQ, EP + Q ⇌ EPQ)`

## Not checked

- An end-to-end `identify_rate_equation` run: `pmap` and worker overhead, expansion, CSV writing and LOOCV fold fits.
- Other reactions (uni_uni, bi_bi_pp, ter_ter, pyruvate carboxylase), reactions with regulators, absolute mode (which has no flat scale direction), optimizers other than CMA-ES, and data beyond 480 points.
- Candidates in classes W, D and Z were excluded from the sample. They are fitted today, and their cost should follow the same per-generation mechanics, but I did not measure it.
- No hand-written never-saturating laws, and saturation of #1858 was not proven for all parameter values.
- Only two heap states for garbage collection.
- The machine is shared aarch64 hardware, so timings carry some noise.

## Files

All under `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/fitcost/`:

- **Raw tables (corrected runs):**
  - `data/fits_fresh2.tsv`, `data/restarts_fresh2.tsv` — main dataset, one row per fit and per restart.
  - `data/fits_alt.tsv`, `data/restarts_alt.tsv` — second dataset.
- **Raw tables (other experiments):**
  - `data/datasize.tsv` — 48 vs 480 points.
  - `data/cma_overhead_*.tsv` — CMA-ES overhead with a trivial objective.
  - `data_run1/fits.tsv` — loaded-heap run.
  - `data/fits_fresh.tsv` — superseded run with the recompiling loop.
- **Sample, pools and data:**
  - `data/selection.tsv`, `data/selection_big.tsv` — the sample with step lists.
  - `data/population_selection_big.tsv` — eligible pool records.
  - `data/data_*.tsv` — the synthetic datasets.
- **Analysis logs:** `logs/pooled.log`, `logs/analyze_fresh2.log`, `logs/analyze_alt.log`, `logs/final_fresh2.log`, `logs/examples_fresh2.log`, `logs/check.log`.
- **Run logs:** `logs/all.log`, `logs/select_big.log`, `logs/fit_fresh2.log`, `logs/fit_alt.log`.
- **Scripts:** `fitlib.jl`, `select.jl`, `fit.jl`, `all.jl`, `check.jl`, `analyze.py`, `examples.py`, `final.py`, `pooled.py`.