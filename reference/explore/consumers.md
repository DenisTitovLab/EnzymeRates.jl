# Consumers of the derivation: fitting, identification, CSV output and dedup

I worked read-only and ran no Julia. Where a claim is my reasoning rather than something I read in the code, it says so. The HPC numbers come from `docs/hpc_results/`, which git ignores (`.gitignore:10`). Every one of those runs is older than v0.6.0 (`8662a6b`, 2026-09-22). They therefore use the old parameter names (`kon_ADP_E`) and older enumeration code.

File abbreviations used below:
- **rd** `src/rate_eq_derivation.jl`
- **tc** `src/thermodynamic_constr_for_rate_eq_derivation.jl`
- **sp** `src/sym_poly_for_rate_eq_derivation.jl`
- **ty** `src/types.jl`
- **id** `src/identify_rate_equation.jl`
- **fit** `src/fitting.jl`
- **en** `src/mechanism_enumeration.jl`
- **tid / trd / tf** `test/test_identify_rate_equation.jl`, `test/test_rate_eq_derivation.jl`, `test/test_fitting.jl`

## 1. Which derivation functions the consumers call

| API | Defined | Called from | What the caller relies on |
|---|---|---|---|
| `compile_mechanism` | en:388-389, which builds `EnzymeMechanism{_sig_of(m)}()` (ty:1130) or the allosteric form (ty:1201-1209) | fit:99; id:342, 461, 491, 846, 855; convenience forwarders rd:62, 64, 419, 492 | The singleton type is the identity for code generation and for the CSV `string(typeof(em))`. It is also `results.best` (id:855, 875). |
| `fitted_params` (`@generated`) | rd:47-50; the generator runs `_dependent_param_exprs(M)` | fit:121, 223-224, 227; id:462-463 | A compile-time-constant tuple of names. Its order is the optimizer's coordinate order (fit:128, 260), the order of user-supplied `lb`/`ub`, and the order `_cv_fold_loss` reads `fit.params` back in (id:895-896). Its length is the beam's `n_params` and is checked against `max_param_count`. Plain mechanisms sort it by name (tc:548-552). Allosteric mechanisms list catalytic names sorted, then active-state `Kreg`s, then inactive-state `Kreg`s, then `L` (rd:1003-1004). |
| `metabolites` (`@generated`) | ty:1322-1327, 1336-1339 | fit:50 (required columns), fit:122; rd:383, 514 | A constant tuple, so the concentration NamedTuple in `loss!` is type-stable (comment at ty:1323-1325). Regulators declared but bound by no step are left out (ty:1054-1066). |
| `rate_equation(m, concs, params)`, Reduced mode (`@generated`) | rd:412-413, 427-433 | fit:138 | Reads both NamedTuples by name (rd:359-362, 381-385). Needs `Keq` and `E_total` in `params`. Returns a signed Float64; the loss uses the sign (fit:139-145). Must not allocate. |
| `rate_equation_string`, Reduced mode | rd:486-492, 526-537 | id:464, for every mechanism not seen before | Deterministic text. It feeds `eq_hash` and the CSV `rate_equation` column. |
| `_rate_eq_dedup_key` | id:262-268; uses `ANNOTATION_SUBSTITUTED` (rd:16) | id:467 | The exact string format (section 2). |
| `rescale_parameter_values` | rd:736-746 | fit:261-265 | Calls `_kcat_forward` (`@generated`, rd:593-606 and 639-726) on the independent parameters plus `Keq`. Uses `_ss_rate_constant_names` (rd:550-555) to decide which names to scale. |
| `_eq_complexity` | rd:319-332 | id:460, before compiling | The V×τ term count from the graph alone; the same number guards derivation (`_assert_derivable`, rd:340-347). For an allosteric mechanism it counts only the active-state graph (rd:332). |
| `==` / `hash` of `Mechanism`, `AllostericMechanism` | ty:971-982 | id:448-452 (`seen`, which stores only the `UInt64`); id:608-616 (`parent_of` Dict); `unique!` at id:646 and 732-734; en:1980-1987 | Structural identity. |
| Exceptions thrown by derivation | e.g. tc:539-541, rd:191-194, rd:619-620 | caught at id:468-470 and 495-497 | Derivation must throw rather than hang, so failures become `FitFailure` rows. |

**Not used by fitting or identification:**
- `parameters`.
- Full mode (`Full` appears only in rd, ty and the export list).
- `_independent_param_count`, which only enumeration uses (en:1130-1131).

**Name conventions the consumers depend on:**
- `Keq` and `E_total` are reserved (fit:129, 263), and so is `L` (rd:716, 1004, 1045).
- `_ss_rate_constant_names` collects the rendered `Kfor`/`Krev` names, so in practice it selects the lowercase `k_` names (ty:1398-1403).
- State tags are `A_`, `I_` or none (ty:1355-1358); regulator constants are `K_<tag><ligand>reg` (ty:1410-1411).
- The CSV gets one column per parameter name (id:230, 246-248).

**Public result:** `IdentifyRateEquationResults.best` is the singleton type (id:58-61). Users pass it to `rate_equation_string` (README.md:65-66).

## 2. How `eq_hash` deduplication works

**Where the key is used:**
- `_process_batch` first drops structures already seen by `hash(m)` (id:448-452).
- On the workers it then applies the complexity cap, compiles, applies the parameter cap, renders the string and hashes it (id:458-471).
- One representative per key not yet in `memo` is fit (id:474-479, 489-498). Every row with that key copies the fit, including the rescaled parameters (id:504-535). `memo` lasts for the whole search (id:687-690).
- `_offer_cv!` keeps the best `n_cv_candidates` distinct keys per parameter count for cross-validation (id:572-585).
- The frontier still keeps every structurally distinct mechanism, duplicates included (id:550-565).

**What the key is.** It is `hash` of the `rate_equation_string` text after two kinds of lines are removed (id:262-268):
- lines starting with `#`, which are the `# Wegscheider constraints:` / `# Haldane constraints:` headers;
- lines containing `"  (substituted into v)"`.

What remains comes from `_equation_text` (rd:505-519):
- the params line: independent names plus `Keq, E_total`;
- the concentrations line;
- one `sym = rhs` line per dependent parameter, sorted by name;
- the `v = E_total * (num) / (den)` line.

The key depends on several canonical orderings that exist for this purpose: independent names sorted (tc:548-552), power-product factors sorted (sp:216-221), monomials sorted (sp:104-111), and dependent lines sorted (rd:500-503, which calls this "load-bearing").

**What it catches.** Only equations whose text is identical. The main case is a kinetic group that was split but that a single-symbol Wegscheider tie folds back together (`_build_wegscheider_rename_map`, rd:89-103). The folded constant disappears from `v`, and its "substituted" line is stripped. `_testhelper_dedup_twins` pins this (tid:1199-1300). It cannot merge different functions:
- Identical text means the same `v`, the same independent parameters and the same dependent definitions.
- Fits are copied by name, which needs identical names in identical order (tid:1235-1236).

**What it misses:**
1. **Renamings.** Each constant is named after its kinetic group's representative step (`_group_rep`, tc:109-110) and its allosteric tag. When grouping or tags differ, the same constant gets a different name.
   - tid:1306-1357 pins one LDH pair: "eq_hash does not see through which tied name survives".
   - In allosteric mechanisms, `:EqualAI` renders a bare name while `:NonequalAI` renders `K_A_…`/`K_I_…`. So a split that thermodynamics ties back still renders differently (my reading of ty:1355-1358 and rd:953-1005).
2. **Reparameterizations and uncancelled common factors** (ROADMAP:34-41).

**Evidence from the HPC runs (pre-v0.6):**
- **Exact-key hits are rare.** Inherited rows were 0.4–14.7% of all rows, summed from each `progress.log`: LDH 07/17 4.1%, LDH 09/17 3.2% with inhibitors and 14.7% without, PFKP 0.4%, PGK 07/12 3.1%, PGK 09/17 2.8% and 4.5%.
- **Equal-loss duplicates are common.**
  - Share of distinct keys whose training loss matches another key's to 8 significant digits at the same count:
    - LDH 09/17 without inhibitors: 61.7%.
    - PGK 09/17 with inhibitors: 84.1%.
    - PGK 09/17 without inhibitors: 61.7%.
  - That is the size of the ROADMAP's "66% under renaming" figure (ROADMAP:73).
  - Equal loss is not proof of an equal function. Among those keys, about 10% have exactly the same set of fitted values as another member, which is the signature of a pure renaming. Another 28–31% have a K or L at the ±15 log bound, which means the fit collapsed to a smaller nested model.
- **The cross-validation pool fills with copies.** In each run's 45 LOOCV rows, 7–21 rows share both loss and `cv_score` with another row at the same count while having different keys.
  - A concrete case is in `ldh_hpc_results/2026_09_17_results_no_comp_inh/loocv_results.csv` at n_params = 9: all five slots hold one model.
  - The five keys carry identical values under different names. 0.000329 is called `K_Pyruvate_ENADH`, `K_A_Pyruvate_ENADH` or `K_A_Pyruvate_E`. 0.000112 is called `K_I_Pyruvate_E` or `K_Pyruvate_E`. The other seven parameters match exactly.

**Costs beyond wasted fits:**
- Renamed copies take cross-validation slots, which defeats the "distinct equations" purpose of `n_cv_candidates`.
- They also use up the `min_beam_width` budget, because `_select_count!` ranks raw losses (id:299-305).
- Each copy is compiled separately.

**What a change in format, naming or dependent choice would do:**
- **Any format change moves every hash.** That includes removing `1 ^ 2`, normalizing the denominator, or renaming residuals (ROADMAP:98-104). This is harmless: nothing reads keys back across runs. The only aggregate pin is the 239-class count over bi-bi `init_mechanisms` (tid:717, 734).
- **Grouping-independent names plus a dependent choice determined by the function (goal 5)** would merge renaming classes.
- **Normalizing to a leading 1 (goal 6)** would merge forms that differ by a common factor.
- **A rename-invariant key breaks the copy-by-name contract** (id:515, 529) unless each row also carries a map from the representative's names to its own.
- **Display and evaluation use different constraint solves.**
  - For plain mechanisms, the string's dependent lines come from a solve with no rename (rd:529), while `rate_equation` uses the renamed solve (rd:431-432). I think they differ only in the folded-tie lines; I have not checked this.
  - Allosteric strings never show folded ties (rd:534-537).
- **The key must give the same value on every worker,** because workers compute it and the main process compares it. Today that relies on `hash(::String)` being unseeded (my inference). A key built from polynomials must avoid hashing that depends on `objectid` or Dict iteration order.

## 3. What the identify run writes to CSV

**Files:**
- `initial_mechanisms.csv` and `equation_search_iteration_N.csv` (id:713-714)
- `loocv_results.csv` and `best_equation.csv` (id:872-873)
- `progress.log` (id:361-370)

**Columns** (from `_rows_to_dataframe`, id:228-250):
- `n_params`, `parent_n_params`, `loss`
- `mechanism_type`, `parent_mechanism_type`
- `rate_equation`, `retcode`, `error`
- `eq_hash` (16 hex characters, id:530)
- `fit_inherited`
- then one column per fitted-parameter name across all rows, sorted, with `missing` where a row lacks it.
- The cross-validation files add `cv_score`, `cv_score_se` and one `cv_fold_<group>` per group (id:850-861).
- The docs list these columns in `docs/src/identify/model_selection.md:134-151`.

**How the mechanism is written:**
- `mechanism_type` is `string(typeof(em))` (id:466). If compiling fails, it falls back to `"EnzymeRates.Mechanism"` (id:341-345).
- `parent_mechanism_type` repeats the parent row's full string (id:524-525).

**Example.** This is from `pgk_hpc_results/2026_09_17_results_w_comp_inh/best_equation.csv`; the whole string is 3,347 characters.
```
AllostericEnzymeMechanism{EnzymeMechanism{(((((:Substrate, :ADP), ((:C, 10), (:H, 15), (:N, 5), (:O, 10), (:P, 2))), ((:Product, :ATP), ((:C, 10), (:H, 16), (:N, 5), (:O, 13), (:P, 3))), …, (((:CompetitiveInhibitor, :ADP), (1,)), …), (1,)), (((((), :E, ((), ())), (((:Substrate, :ADP),), :E, ((), ())), (:Substrate, :ADP), false), …
… (1, (:EqualAI, :EqualAI, :OnlyA, :OnlyA, :EqualAI, :EqualAI, :EqualAI, :OnlyA, :EqualAI)), ()}
```
The current format differs: since `e59ebd6` (2026-10-05) a step stores tuples of consumed and released metabolites. `test/test_pivot_priority_regression.jl:17` hard-codes a string in the current format.

**Size.** In PGK 09/17 iteration 3 (4,949 rows, 32.7 MB, 100 columns):
- `mechanism_type` averages 2,721 characters, `parent_mechanism_type` 2,515 and `rate_equation` 1,052.
- The two type columns therefore make up about 79% of the file.
- The selected row's `rate_equation` is 29,071 characters.

**Does anything read these files back?**
- The package source only writes them (`CSV.write` at id:219, 872-873).
- The tests round-trip `mechanism_type` with `Core.eval(EnzymeRates, Meta.parse(…))` (tid:257, 433-434).
- Denis's `pgk_hpc_results/analyze_fit_results.jl:7-21` reads only `n_params`, `cv_score` and `cv_fold_*`.
- Loading old CSVs was declined (ROADMAP:162-163).

**Readability problems:**
- **The mechanism string is not readable.** It carries atom formulas, metabolite-type tags, an empty `((), ())` residual on every form, and `true`/`false` for rapid-equilibrium versus steady-state. You cannot read or grep the steps, kinetic groups or allosteric tags.
- **Every child row repeats its parent's full string.**
- **Round-tripping means eval'ing CSV text, and old files no longer decode.** The 09/17 strings store a single-metabolite `consumed`, which `_mechanism_from_sig` (ty:1074-1096) would reject (my inference).
- **`show` is not a faithful listing.**
  - `show(::EnzymeMechanism)` (ty:1255-1265) prints a readable flat list of steps but drops kinetic grouping.
  - The allosteric `show` (ty:1269-1285) keeps groups and tags.
  - `Mechanism` just dumps its fields (ty:913).
- **The equation text has artifacts.**
  - `1 ^ 2` appears because `_reg_site_expr` returns the literal `1` for a site with no ligand in that conformation, and `_power_expr` then raises it (rd:1030-1040, 1145-1146).
  - Inactive-state terms print the metabolite first (`S * k_I_ES_to_EP / …`), because only active-state symbols count as parameters when sorting (rd:1129-1131; fixture trd:2012-2025).
  - Residual names such as `k_EA_to_EP_res_+A_-P` are not valid identifiers, in the text or in CSV headers (trd:1577-1595).

## 4. What compiles once per mechanism, where, and why

**Parallelism.**
- The package uses only `Distributed.pmap` (EnzymeRates.jl:29; id:458, 489, 598, 845; en:1996). `Threads` does not appear in `src`.
- `pmap` has no worker affinity. One mechanism's render (pass 1), fit (pass 2) and cross-validation folds run on arbitrary workers, and each worker compiles whatever it touches.
- The cross-validation grid spreads one candidate's folds over up to one worker per group (id:842-847).
- With no workers, everything compiles in the main process (`docs/src/identify/parallel.md:17`).

**Compiled per mechanism type M:**

| Item | Where | Specialized on | Note |
|---|---|---|---|
| The singleton type | ty:1130, 1201-1209 | the Sig (about 2.7 KB of text in PGK) | Kept for the life of the process; an allosteric mechanism also creates its catalytic `EnzymeMechanism` type |
| `fitted_params`, `parameters`, `metabolites` | rd:38-50; ty:1322-1339 | M | Tiny bodies; the generators run the constraint solve or the Sig decode |
| `rate_equation` and its 3-argument forwarder | rd:412-433 | M, plus the concentration and parameter NamedTuple types | The largest generated body; grows with V×τ |
| `_kcat_forward` | rd:593-606, 639-726 | M, plus the parameter NamedTuple type | Called once per fit; the allosteric version emits 2^(number of ligands) regulator corners per pattern (rd:676-719) |
| `rescale_parameter_values` | rd:736-746 | the parameter NamedTuple type (`m` is `@nospecialize`) | Shared by mechanisms with the same names |
| `FittingProblem` constructor and type | fit:48-59 | M and the data type D | D is the full columns in pass 2 but `SubArray` views in cross-validation (id:891-893), so it compiles twice per candidate |
| `loss!` | fit:119 | x, M, D | x is a `SubArray` from CMA-ES (`CMAEvolutionStrategy/src/optimizer.jl:200-211`) and a `Vector` in the cross-validation test loss (id:896) |
| `fit_rate_equation` (with its keyword body), `_cv_fold_loss` | fit:216; id:884 | M, D | |
| The Optimization/CMA-ES layer | see below | M, D | All of `solve` → `init` → `__solve` → `minimize` recompiles per mechanism |

How the Optimization/CMA-ES layer picks up M:
- `OptimizationProblem` takes `P = FittingProblem{M,D}`. The cache stores `ReInitCache(prob.u0, prob.p)` (OptimizationBase 2.14.0 `src/cache.jl:8-23, 36-39`; 5.7.0 has the same at lines 50-69, 85-88 and 182).
- The CMA-ES `__solve` builds `_loss` and `_cb` closures over that cache. This is OptimizationCMAEvolutionStrategy 0.3.13 `src:97, 117-131`, and versions 0.3.2–0.3.12 are the same.
- `minimize`, `run!`, `evaluate` and `BasicLogger{F}` specialize on those closures (`CMAEvolutionStrategy.jl:32-51, 114-118`; `log.jl:2-9`).

**Barriers already in place.**
- The functions that turn a singleton back into a concrete mechanism, and the derivation helpers, take `@nospecialize` (ty:1074, 1139, 1173, 1187-1188; rd:300, 353, 379, 505, 521-537, 550, 737, 1007, 1123).
- The pipeline holds plain `Union{Mechanism, AllostericMechanism}` values and only compiles inside the `pmap` closures.
- There is no barrier between `FittingProblem` and the solver.

**The derivation repeats within one worker** (counted from the code; on another worker the count starts again):
- Pass 1:
  - The `fitted_params` generator runs 2 constraint solves.
  - `rate_equation_string` runs 5 more solves and one King–Altman derivation: two through `_dependent_param_exprs`, one with no rename, and two that build rename maps.
- Pass 2:
  - The `rate_equation` generator runs 3 solves and one King–Altman.
  - The `_kcat_forward` generator runs 3 solves and one King–Altman.
- That totals about 13 cycle-basis solves and 3 King–Altman derivations per fitted mechanism on one worker.
- Allosteric mechanisms double the King–Altman work per call. `_dependent_param_exprs(am)` also re-derives both conformations whenever a rename exists (rd:984-993).

**What I expect dominates (my estimate; the measuring agent should confirm):** two costs per fitted mechanism per worker.
- (a) The Optimization/CMA-ES layer specialized on `FittingProblem{M,D}`. This is a fixed cost made of many method instances, and it is paid twice for cross-validation candidates.
- (b) Compiling the generated `rate_equation` body, which wins for large equations.
- `_kcat_forward` comes third. `fitted_params` and `metabolites` are negligible as code.

Useful measurements:
- Time the first and second call of `fit_rate_equation(…; maxiters=1)` on a new mechanism. Repeat on a renamed twin that has the same names.
- Track the change in `Base.jit_total_bytes()` and `Sys.maxrss()` per new mechanism.
- Split SnoopCompile inference time between the SciML packages, CMAEvolutionStrategy and EnzymeRates.

## 5. What the fitter needs from a rate equation

- **Two ordered name lists that are constant at compile time:** N fitted parameters and K metabolites (fit:121-137). The inner loop itself needs no names; it only maps `x[i]` to parameter i and data column j to metabolite j.
- **An allocation-free evaluator** of the signed rate at fixed `Keq` and `E_total = 1.0` (fit:129), which computes the dependent parameters internally.
- **Log-space parameters:**
  - Each parameter is `exp(x[i])` (fit:128).
  - Every parameter kind gets the default box [-15, 15] (fit:223-224).
  - Starting points are `clamp.(randn·2)` (fit:250).
  - `lb`/`ub` are positional, in `fitted_params` order.
- **Rescaling in relative mode** needs kcat as a function of the independent parameters plus `Keq`, and the set of steady-state rate-constant names (fit:261-265).
  - This only works if `v` is homogeneous of degree 1 in those constants (`docs/src/fitting/normalized_vs_absolute.md:72-77`). Each steady-state step enters a cycle row as `+c·log kf − c·log kr`, so every dependent k also scales by one factor (my reading of tc:272-283).
  - It follows that the per-group-centered loss (fit:158-169) is flat along a common scale in every relative-mode fit (my inference).
- **Data:** `group`, `Rate` and one column per `metabolites(m)` (fit:70-90). There is no `E_total` column.
- **Cross-validation** refits each fold and scores it with `loss!` at `log.(fit.params)`, in NamedTuple order (id:894-896).
- **The identify search also needs:**
  - `n_params`;
  - a deterministic equation key;
  - text for the CSV;
  - a size bound computable before derivation;
  - derivation errors thrown as exceptions, not hangs.

## 6. Constraints and opportunities for the redesign

**Load-bearing; change only on purpose and with Denis's agreement:**
1. **Allocation and speed tests.** `rate_equation` called with NamedTuples allocates nothing and runs under 120 ns (trd:460-476, 881-887). `loss!` allocates nothing (tf:201-217, 234-244).
2. **The `fitted_params` order.** It ties together optimizer coordinates, `lb`/`ub`, `fit.params` and `_cv_fold_loss`. Its length is the beam's parameter count.
3. **The `eq_hash` contract.** An equal key must mean the same function and the same names in the same order.
4. **Rescaling** by name class, and the degree-1 homogeneity it assumes.
5. **Keys that are identical on every process.**
6. **Tests tied to formats:**
   - the Sig round-trip tests and the hard-coded Sig fixture;
   - the byte-identical allosteric string (trd:2012-2025);
   - the 239-class count;
   - the documented CSV columns.

**Free to change.** Nothing outside the tests reads these back:
- the Sig layout;
- the string format;
- how the key is built;
- how kcat is computed;
- whether `fitted_params` and `metabolites` are `@generated`, as long as the inner loop stays type-stable.

**Opportunities, in order of expected payoff:**
1. **Keep M out of the solver's types.**
   - Keep one compiled numeric kernel per mechanism behind a zero-allocation barrier, for example a positional `v(c::NTuple{K}, θ::NTuple{N}, Keq)` held through FunctionWrappers or stored with its type erased.
   - The Optimization/CMA-ES layer would then compile once per data type, which goes straight at the uncollectable-code memory problem.
   - A plain dynamic dispatch per `loss!` call would box one Float64 and fail the zero-allocation test (my inference).
2. **Derive each mechanism once per process.**
   - Cache a concrete `derive(m)` result (polynomials, names, dependents, kcat groups) on each worker, and have the generators and the string read from it.
   - That removes about 11 constraint solves and 2 King–Altman runs per fitted mechanism.
   - Also copy the cross-validation fold columns into plain `Vector`s so D matches pass 2. The copy is tiny next to a second compile of the solver layer.
3. **Compute kcat without generated code.** Either evaluate the grouped polynomials numerically once per fit, or anchor kcat at the assay concentrations as ROADMAP:90-97 plans.
4. **Build keys from the polynomials, plus a second key that sees through renaming.**
   - For example, bucket by the concentration monomials of the numerator and denominator plus the parameter count, then check for an exact bijection between names.
   - Use it to keep cross-validation slots and the beam budget distinct.
   - Carry the name map so a fit can be copied across names.
5. **Write a readable CSV.**
   - Use a canonical listing in the DSL's own syntax, one line per kinetic group with tags and regulatory sites, plus a short stable `mechanism_id`.
   - Use a `parent_id` column instead of repeating the parent's string.
   - Parse it back with a parser shared with the DSL instead of `eval`.
6. **Removing `reaction` from `Mechanism` is easy here.** The identify code already has `prob.reaction` wherever one is needed (`compile_mechanism`, `FittingProblem(::Mechanism)`, `_degenerate`). `hash` and `==` would also stop rehashing the whole reaction for every `seen`, `unique!` and `parent_of` lookup.

**Other things I noticed, all latent and none verified by running:**
- **`loss!` is not thread-safe.** It writes the shared `log_ratios_buffer` (fit:120, 141-145), so CMA-ES `multi_threading=true` passed through `solver_kwargs` would race (`optimizer.jl:203-208`). The Float64 buffer also rules out optimizers that need automatic differentiation.
- **A metabolite name can shadow a parameter.** `_rate_body` unpacks `params` and then `concs` into one scope (rd:381-385). A metabolite named `L`, `Keq` or `E_total` would therefore silently replace that parameter. I found no reserved-name check in ty or `src/dsl.jl`.
- **`eq_complexity_filter` undercounts allosteric mechanisms.** It measures only the active-state graph (rd:332), while the generated MWC body evaluates both conformations plus the powers.