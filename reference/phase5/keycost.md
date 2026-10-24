# Key, complexity filter and codegen (R5.5, R5.6, R5.7; engine.md Q8, Q15; key.md §5 items 6, 7, 8, 11, 12)

Paths are relative to `scratchpad/phase45/research5/keycost/` unless they start with
`scratchpad/`. Inputs were regenerated with `graph_mechanism_wt` at 7205ba4 (`env/`, census mkenv
recipe plus CMA-ES and CSV; `logs/mkenv.log`), Julia 1.12.5 unless stated. Slice = the stock
`scratchpad/slice/src/`; the key is its `canonical` with concentrations labelled by name
(`klib.jl` `canonical_named`). Every Julia run went through `phase45/jl_research.sh`; peak RSS
2.0 GB (`logs/export5_terter.log`).

## 0. Answer

- **R5.5:** with a regulator and the census reactions (two subunits), counting both conformations
  drops 0 candidates in every round: 2 rounds as asked, and 5-6 rounds sampled uniformly or
  largest-first (`logs/r55_table.log:2-27`).
- The rule bites only on one-subunit allosteric laws:
  - sampled uniformly, it drops 0.05-0.48% of each round (`logs/r55_table.log:29-31, 38-40, 48-49`);
  - walked largest-first, it drops 10-25% of rounds 3-4, and 55% of bi-bi round 4 (`:35-36, 44-45, 53-54`).
- The dropped laws cost what the plain laws the filter already drops cost:
  - 0.092 ns per V×τ unit, against 0.091 for admitted plain laws (`logs/r55_cost_ratio.log:1, 3`);
  - recommendation: the sum, default unchanged (F1).
- **R5.6 passes.** On 1,390 ter-ter laws up to 3,438 key terms (V×τ ≤ 4,896), no law hit the leaf cap.
  - Every law needed exactly 1 leaf (`logs/kt_terter_summary.log:3-9`).
  - Key time is 1.21% of a round's fit time with the filter at 5,000, and 0.34% at the default 337 (`logs/key_vs_fit.log:6-7`).
- **Item 12:** at the 100,000-leaf cap, `canonical` throws:
  - toy symmetric laws hit the cap after 4.7-17 s (`logs/leafcap.log:4, 7`);
  - proposal: cap 1,000 leaves, then a per-mechanism key (F3).
- **R5.7:** 0 allocations on 578 timed laws (`logs/c_calib.log:1`).
  - Worst local times per denominator term: 0.108 ns on the gate laws, 0.151 ns on any law with D ≥ 300 (`:2, 5-10`).
  - Specs: at most 13.1 ns, against 41.6 ns for the package (`:11`; `data/x5_small.tsv`, spec rows).
  - Provisional c = 1.0 ns per term (F8); the PR's CI sets the final value.
- **Q15 / item 7:** the gate should count D after substitution.
  - It equals today's count on 67 of 67 plain laws, so the pins [396, 396, 396] and [589, 799, 858] stay (`logs/dterms_vs_pkg.log:1`).
- **Item 6:** a 16-hex key reads back as a number from 2.1e-3 of one-row CSVs (`best_equation.csv` is one row); a 64-hex key from none (F4).
- **Item 8:** test with a subprocess plus pinned key digests. Body digests vary with the Julia version: 2 of 78 differ between 1.12 and 1.13 (F5).
- **Item 11:** for q > 2, emit `x^(one(x)/q)`.
  - It costs the same as today's Float64 exponent and is exact in BigFloat: error 1.5e-77, against 8.5e-17 (`logs/roots.log:1-2`) (F6).

## 1. R5.5 and Q8: the filter over both conformations

Method (`pops5.jl`):
- Round 0 is identify's base tier: `seed_mechanisms` with R required, or `init_mechanisms` for bi-bi without a regulator, then `_base_tier`.
- Round k expands parents drawn from round k−1 within today's filter (V×τ(active) ≤ 337).
  - `uniform` draws parents at random (the census method).
  - `largest` takes the parents with the largest V×τ.
  - The quota is 1,500 parents per round (round 1: all parents).
- For each candidate:
  - V×τ(active) is `_eq_complexity(m)`;
  - V×τ(inactive) is `_eq_complexity(_inactive_groups(am))`;
  - "same law" means every step is kept and every group is `:EqualAI`.
- Runs: 20-101 s each, RSS ≤ 1.52 GB (`logs/pops5_*.log`, last lines).

Rules counted, each over the candidates today's filter admits:
- sum: V×τA + V×τI > 337;
- distinct: the inactive part counted only when its law differs;
- max: the larger of the two;
- sum > 674: the default retuned to 2 × 337.

| Reaction, mode | Rounds | Newly over, per round (share of admitted) | Source |
|---|---|---|---|
| bi-bi + R, n = 2 (asked) | 0-5 uniform, 0-6 largest | 0 every round; the largest admitted sum is 56 and 40 | r55_table.log:2-14 |
| ping-pong + R, n = 2 (asked) | 0-5 uniform, 0-6 largest | 0 every round; the largest sum is 50 and 72 | r55_table.log:15-27 |
| bi-bi + R, n = 1 | uniform 1-3 | 24 (0.09%), 38 (0.30%), 59 (0.48%) | r55_table.log:29-31 |
| bi-bi + R, n = 1 | largest 2-4 | 230 (3.3%), 1,028 (17.4%), 1,256 (25.1%) | r55_table.log:34-36 |
| ping-pong + R, n = 1 | uniform 1-3 | 12 (0.05%), 33 (0.26%), 24 (0.20%) | r55_table.log:38-40 |
| ping-pong + R, n = 1 | largest 2-4 | 120 (1.9%), 553 (10.4%), 920 (19.8%) | r55_table.log:43-45 |
| census bi-bi, no regulator | uniform 2-3 | 19 (0.05%), 165 (0.42%) | r55_table.log:48-49 |
| census bi-bi, no regulator | largest 2-4 | 48 (0.2%), 3,528 (15.0%), 8,899 (55.3%) | r55_table.log:52-54 |

- **Why the n = 2 rounds never grow.** Every regulated seed of the census reactions has two subunits. Inferred, not run: the hyperbolic filter blocks the flips that would grow V×τ there (`graph_mechanism_wt/src/mechanism_enumeration.jl:899`).
  - The largest admitted V×τ(active) stays at 9-36 through 6 rounds (`r55_table.log:2-27`).
  - Seed counts: 2,112 and 1,856 (`logs/pops5_bibiR_uniform.log:2`, `pops5_pingpong_uniform.log:2`).
- "distinct" equals "sum" in every row: no newly-over candidate has an inactive law identical to its active one (`r55_table.log`, columns 4-5).
- "max" drops 0, except 192 in bi-bi largest round 4, where V×τI > V×τA (`r55_table.log:54`). It
  undercounts, because the generated code evaluates both conformations.
- Retuning to "sum > 674" admits 68-4,920 plain candidates per round that today's default drops (`r55_table.log:47-53`, column 6).
- Cost check on the slice (`r55_sample.py` → `data/kt_r55.tsv`): 160 newly-over allosteric laws, 60 allosteric laws under the threshold, 60 plain laws in (200, 337] and 60 in (337, 720].
  - Newly over: V×τ sum median 351 (338-719); 1,096 operations; 33.5 ns median, 74.5 ns max; 11 ms key (`logs/r55_cost.log:4`).
  - Plain laws in (200, 337]: 864 operations, 27.1 ns (`:2`).
  - Plain laws over 337: 1,208 operations, 34.8 ns (`:5`).
  - Per V×τ unit: 0.092 ns and 3.04 operations for the allosteric laws, against 0.091 ns and 2.71 for the plain laws (`logs/r55_cost_ratio.log:1, 3`). Sharing locals between conformations does not make the sum overstate cost.
  - 0 allocations in all 340 (`logs/r55_cost.log`, last column).
- The package's V×τ per conformation equals the slice's `complexity` with `inactive_mask` in 1,808 of 1,808 rows (`logs/cx_pkg_vs_slice.log:1-3`).

## 2. R5.6 and item 12: key cost on ter-ter laws

Inputs (`export5.jl terter`, `logs/export5_terter.log`):
- `init_mechanisms(ter_ter)` gives 250,855 mechanisms.
- Their V×τ bins: ≤ 50: 162,067; 51-337: 50,562; 338-1,000: 12,588; 1,001-2,000: 3,696; 2,001-5,000: 2,754; over 5,000: 19,188 (`:3-8`).
- Sampled: 150, 150, 150, 200 and 300 per bin, plus the 40 largest at V×τ 4,896, plus 400 children from round 1 of 25 parents (274 allosteric; `:9`).
- `keytime.jl` derives the slice law and times `canonical_named`. A key over 5 ms is timed twice and keeps the minimum; both runs must give the same key.

| Key-part terms | Laws | Key ms: median / p90 / max | Law ms, median | Key ÷ law, median | Leaves, max |
|---|---|---|---|---|---|
| ≤ 50 | 334 | 0.13 / 0.27 / 2.4 | 0.2 | 0.81 | 1 |
| 51-337 | 245 | 1.70 / 3.23 / 5.0 | 0.6 | 2.40 | 1 |
| 338-1,000 | 241 | 16.1 / 24.2 / 32.1 | 7.7 | 2.14 | 1 |
| 1,001-2,000 | 324 | 46.5 / 90.0 / 115 | 26.2 | 2.27 | 1 |
| 2,001-5,000 | 246 | 124 / 212 / 351 | 38.8 | 3.53 | 1 |

(`logs/kt_terter_summary.log:3-7`)

- 0 errors, 0 leaf-cap hits. Every law took 1 leaf and 1 search node: refinement alone separates all constants (`:1, 9`).
- Over all 1,968 rows timed (1,800 distinct laws, including specs, gate laws, the R5.5 sample and allosteric laws), every law took 1 leaf (`logs/compile_and_rows.log:3`).
- Key time grows as key-part terms^1.63; law time as terms^1.54 (`kt_terter_summary.log:8`).
- On large laws the key costs 2-3.5× the law (median, `:5-7`).
- Slowest key: 351 ms (round-1 child, 3,363 key terms; `:11`).
- Allosteric children: median 0.17 ms, max 159 ms (`:10`).
- **Fit time, measured** (`fittime.jl`; `logs/fittime_1_4.log:2-5`, `fittime_5_8.log:2-5`, `fittime_9_12.log:2-5`):
  - 12 ter-ter mechanisms: 2 per bin and 2 allosteric children.
  - The fit is identify's: 20 CMA-ES restarts, relative loss, 48 synthetic points from the true model.
  - Per candidate:
    - 0.21 s in the ≤ 50 bin;
    - 0.39 s in 51-337;
    - 0.57 s in 338-1,000;
    - 0.76 s in 1,001-2,000;
    - 1.33 s in 2,001-5,000 (fixed cost plus 20 restarts; `logs/key_vs_fit.log:1-5`).
  - These are easy fits: the data come from the true model. The census critic gives 80 ms fixed plus 0.1-1.5 s of
    optimization (`scratchpad/census/reports/critic.md:31`), so real shares are lower.
- Key ÷ fit per candidate (mean): 0.1%, 0.7%, 2.8%, 5.0% and 7.9% by bin; the worst single candidate is 21.9% (`key_vs_fit.log:1-5`).
- **A round** (all 231,667 init mechanisms ≤ 5,000, weighted by bin): key 0.23 h against 18.6 h of serial fits, **1.21%**. With the default 337: 3.0 min against 14.8 h, **0.34%** (`:6-7`).
- The key runs after the complexity filter (PASS1 checks complexity first, `graph_mechanism_wt/src/identify_rate_equation.jl:507`), so a filtered candidate pays nothing.
- **Item 12, the leaf cap.** `canonical` throws `"canonical labelling: over 100000 leaves"` (`scratchpad/slice/src/key.jl:102`). Toy laws (`leafcap.jl`):
  - 9 symmetric constants (9! leaves) hit the cap after 4.65 s;
  - 17 tied pairs (2^17) after 17.2 s;
  - 16 pairs (65,536 leaves) finish in 10.6 s (`logs/leafcap.log:3-7`).
  - Inferred, not run: on a 3,000-term law each leaf costs up to the 1-leaf key time (≤ 0.36 s), so reaching the cap can take hours.
  - The cap bounds leaves, not time.

## 3. R5.7, Q15 and item 7: codegen locals, the gate's term count and c

Inputs:
- the 42 specs, the six gate laws (`test_rate_eq_derivation.jl:1253-1302`) and the 30 largest ter-bi seeds (`export5.jl small`);
- 160 large ter-ter laws (`data/x5_tt_code.tsv`) and the 340 R5.5 laws.

The time per call is the minimum of 5 × 10,000 calls through a `@noinline` barrier (`klib.jl` `time_call`).

| Laws | Allocations | Max ns | Worst ns per D term | Source |
|---|---|---|---|---|
| 42 specs | 0 | 13.1 (the package: 41.6) | – | c_calib.log:11; x5_small.tsv |
| bi-bi gate laws, D 396 × 3 | 0 | 42.7 (the package: 61-63) | 0.108 | c_calib.log:5-7 |
| ter-bi gate laws, D 589 / 799 / 858 | 0 | 68.9 (the package: 87-128) | 0.080 | c_calib.log:8-10 |
| 30 largest ter-bi seeds, D 589-858 | 0 | 75.2 | 0.103 | kt_small.tsv |
| 160 ter-ter laws, D up to 3,236 | 0 | 169.7 (9 laws over 120 ns, D 2,017-3,009) | 0.091 | kt_tt_code_summary.log:14-15; compile_and_rows.log:2 |
| all 578 laws timed | 0 | 169.7 | 0.151 (an allosteric bi-bi law, D 446) | c_calib.log:1-3 |

- Compile plus first call on the ter-ter set: median 0.17 s, max 0.44 s; the program and body take 0.02 s median (`compile_and_rows.log:1`).
- **Q15 / item 7, which term count.**
  - *D after substitution* (lowest terms):
    - equals the package's `_raw_symbolic_rate_polys` count on 67 of 67 plain laws, so the pins [396, 396, 396] and [589, 799, 858] stay (`logs/dterms_vs_pkg.log:1`);
    - ns grows as D^0.75 on the small set and as D^0.86 on ter-ter (`logs/kt_small_summary.log:18`; `kt_tt_code_summary.log:18`). Per-term cost falls with size, so a bound calibrated at 396 terms holds above it.
    - On the mixed R5.5 sample, D tracks ns poorly: slope 0.33, against 0.78 for operations (`logs/kt_r55_summary.log:18`). Those laws all run under 75 ns, inside the 120-ns floor (`:14`).
  - *V×τ*: worst fit; on ter-ter the log-log slope is 0.04 and ns per V×τ spans 7.9× (`:16, 18`).
  - *Body operations*: spans 3.4×, against 4.1× for D (`:15, 17`). Slightly tighter, but the test would have to walk the body; D is the spec's own wording.
  - For an MWC law, D = DA + DI, with a flip-form-only inactive conformation counting 1 (`klib.jl` `den_terms`).
- **c by phase 1's rule.** The rule: twice the worst CI time per term, rounded up to one significant figure (`scratchpad/phase1_logs/pr_body_draft.md:43-45`).
  - Local worst on the gate set: 0.108 ns. On local readings alone the rule gives 0.3.
  - Today's CI-to-local ratios, from the package on the bi-bi gate laws:
    - local: 63.4/396 = 0.160 ns per term (`data/x5_small.tsv`, gate rows);
    - CI: 0.313 on ubuntu, 0.459 on windows, 0.732 on macOS (`graph_mechanism_wt/test/test_rate_eq_derivation.jl:1309-1310`);
    - ratios: 1.96, 2.87 and 4.58.
  - Predicted worst CI reading for the new code: 0.108 × 4.58 = 0.49 ns, so the rule gives **c = 1.0**.
  - The same rule on today's readings would give 1.5; add414a's 0.8 cleared the readings rather than doubling them.
  - 0 local laws exceed max(120 ns, c·D) at any c ≥ 0.2 (`c_calib.log:13-17`).
- **Gate harness finding.** The package helper's inlined loop measures 4 of 42 specs at under 0.5 ns, because the call is hoisted out of the loop (`c_calib.log:11`). Those specs pass the gate without being measured (F9).
- **Julia 1.10, the compat floor.** The slice fails on 5 of 78 laws (allosteric specs without sites).
  - `AlloLaw`'s parametric constructor gets an untyped empty `map` over sites (`data/kt_small_110.tsv`, rows with "error").
  - The other 73 match 1.12's keys and bodies (`logs/xversion.log:5`).
  - Phase-5 code must give the site list an explicit element type and run its tests on 1.10.

## 4. Items 6, 8 and 11

- **Item 6, key length** (`csvkey.jl`; CSV.jl 0.10.17):
  - 200,000 one-row CSVs per length, written and read back.
  - 16 hex digits:
    - 305 values read back as Float64 (1.53e-3) and 112 as Int64 (5.6e-4); the total, 2.09e-3, matches the analytic 2.06e-3 (`logs/csvkey.log:2, 6`);
    - Parsers.jl takes `e` and `f` as exponent markers, so `7170f82284841813` reads back as `Inf` (`:4`);
    - an all-digit key loses its leading zeros;
    - 64 hex digits: 0 of 200,000 (`:5`).
  - A column with any hex letter stays String (two-row example, `:8`), so the risk sits in one-row and few-row files: `best_equation.csv` is always one row (`identify_rate_equation.jl:921`).
  - False merges at n keys: 64-bit 2.7e-8 at 1e6 keys and 2.7e-4 at 1e8; 256-bit 4e-66 (`:9-11`).
  - Phase 3's 16-hex `mechanism_id` and `parent_id` (R16) carry the same 2.1e-3 hazard in one-row files.
- **Item 8, second-process test.**
  - Within 1.12, a second process gives identical keys and bodies for 78 of 78 laws (`logs/xversion.log:1`).
  - Through a subprocess (as `test_compile_budget.jl` spawns them): 48 of 48 (`:4`).
  - 1.13 gives identical keys for 78 of 78, but **bodies differ for 2 of 78** (`:2, 6-7`).
    - Cause (`logs/bodydiff_112.log` against `bodydiff_113.log`): dependent locals are numbered in order of first use.
    - That order follows `Dict` iteration over a polynomial's terms (`slice/src/codegen.jl:47-55`, `poly.jl:152`), and Julia 1.13 hashes differently.
    - The operation counts are equal; only the local names and order change.
  - Cost: package load 0.9 s per subprocess; keys and bodies of 48 laws in a fresh slice process 9.8 s, Julia start included (`logs/subproc.log:1-3`).
- **Item 11, roots for q > 2** (`roots.jl`; 100,000 Float64 values, 2,000 BigFloat-256 values):
  - `x^(1/q)` with a Float64 exponent (the slice, `codegen.jl:94`): 9.05 ns.
    - BigFloat error: 8.5e-17 for q = 3, 5.1e-17 for q = 5, 4.3e-17 for q = 6 (`logs/roots.log:1, 9, 12`). That fails a 1e-20 oracle.
  - `x^(one(x)/q)`: the same 9.06-9.10 ns, BigFloat error ≤ 1.5e-77 (`:2, 6, 10, 13`).
  - `cbrt`: 3.14 ns (`:3`). `sqrt(sqrt(x))`: 0.98 ns (`:7`). Float64 errors stay within 1.3 eps for every form.
  - No regenerated law needs a fractional exponent: 0 of 1,800 have an exponent denominator ≥ 2 (`kt_terter_summary.log:9`; `kt_tt_code_summary.log:7`; `kt_small_summary.log:7`; `kt_r55_summary.log:7`). Only hand-written fixtures (the slice's `cube_root`) reach q > 2.

## 5. Proposed rulings (ruling — why — cost if wrong)

- **F1 (Q8, R5.5).** `eq_complexity_filter` compares V×τ(active) + V×τ(inactive) with the threshold. This is the cap's own measure through one function, and a flip-form-only inactive conformation counts 0.
  - The default stays 337. The docstring says that a one-subunit K-type variant of the 336-term random-order bi-bi counts 672.
  - — Why:
    - 0 newly over in every round of the asked populations;
    - 0.05-0.48% per round when sampled uniformly;
    - the dropped laws cost what dropped plain laws cost (section 1);
    - "distinct" equals "sum"; "max" undercounts; retuning admits up to 4,920 plain laws per round.
  - — Cost if wrong: a beam walking toward large one-subunit allosteric laws loses 10-55% of deep rounds' candidates (V×τA 169-337). The user raises the keyword.
- **F2 (R5.6).** The key ships as prototyped (refinement, then individualization) with concentrations labelled by name. It is computed only for candidates that pass the filter.
  - — 1.21% of a round's fit time at the 5,000 cap, 0.34% at 337.
  - — Per candidate it reaches 8-22% above V×τ 2,000. If rounds full of 5,000-term laws matter, the next step is a cheaper encode.
- **F3 (item 12).** Lower the leaf cap to 1,000. At the cap the candidate gets a per-mechanism key (SHA-256 of a marker plus its `mechanism_id`), is fitted on its own, and the log names it; it is not a failure row.
  - — No real law needs a second leaf (1,800 of 1,800). A failure row would drop a valid candidate.
  - — The cap bounds leaves, not time; at 100,000 it could take hours on a large law.
  - — Cost if wrong: a law with over 1,000 symmetric labellings loses dedup, at one extra fit per twin. Automorphism pruning (about 30 lines) is the exact alternative, not needed today.
- **F4 (item 6).** The `key` column holds the full 64-hex SHA-256.
  - — Why: 16-hex values read back as numbers from 2.1e-3 of one-row CSVs; 64-hex from none; the spec says SHA-256.
  - — Cost: 48 more characters per row, against a median ter-ter law of 2,630 characters.
  - Open for Denis: give `mechanism_id` and `parent_id` (R16) 64 digits too, or document reading them with `types = String`.
- **F5 (item 8).** Two parts:
  - one subprocess in `test_compile_budget.jl`'s pattern computes keys and body digests of the 42 specs and 6 gate laws, which must equal the parent's (about 10 s);
  - key digests of about 5 specs are pinned, since keys agree across 1.10, 1.12 and 1.13.
  - Body digests are not pinned. Optionally, number the dependent locals by solve column, which makes bodies version-independent (about 3 lines).
  - — Why: the spec's success criterion 3 asks for identical keys and code in a second process; the 2 of 78 body differences on 1.13 would fail pinned bodies.
  - — Cost: about 10 s of suite time.
- **F6 (item 11).** Codegen emits `sqrt` for q = 2 and `x^(one(x)/q)` for q > 2.
  - — Exact in BigFloat at Float64's cost; one code path.
  - — None: `cbrt` would save 6 ns on laws no enumeration produces.
- **F7 (Q15, item 7).** The size gate counts D after substitution, with an MWC law counting DA + DI.
  - — It equals today's count on 67 of 67 laws, so the pins stay; it is the spec's wording; V×τ tracks ns worst.
  - — Cost if wrong: a law whose operation count outgrows D gets a looser bound. MWC laws track ns less well through D (slope 0.33 on the R5.5 sample), but the 120-ns floor covers them today.
- **F8 (R5.7).** c follows phase 1's rule on the phase-5 PR's CI: two runs on three runners, twice the worst reading, rounded up to one significant figure.
  - The provisional value for the first CI run is 1.0 ns per term.
  - — Local worst 0.108 times today's worst CI-to-local ratio of 4.58, then doubled. Local readings alone give 0.3.
  - — Too loose misses only a regression above 2× on 600+-term laws, which the 0-allocation gate partly covers; too tight flakes, as 120 ns does today.
- **F9 (gate harness, needs Denis: CLAUDE.md names this helper).** `_testhelper_test_rate_equation_performance` times its calls through a `@noinline` barrier.
  - — 4 of 42 specs are hoisted to under 0.5 ns, so the gate checks nothing for them.
  - — Cost: about 2 lines, and readings rise by the call cost of about 1-2 ns.

## 6. On phase 3's types (S3_design.md sections 1, 7)

- **Filter.**
  - The slice's `complexity(m, mask)` (`slice/src/engine.jl:207-220`) ports line for line:
    - `m.edges[mask]`;
    - `e.from`/`e.to` for `a`/`b`;
    - `e.rapid_equilibrium` for `re`;
    - segments from `_segments(m, mask, _flags(m))`.
  - Allosteric: `complexity(catalytic(am), trues) + complexity(catalytic(am), _inactive_mask(am))`.
  - The package's `_inactive_groups` and the slice's `inactive_mask` agree on 1,808 of 1,808 laws, so `_inactive_mask(am)` (S3_design.md:982) gives these counts.
- **Key.** `canonical_named` reads only the law (`key_parts`, `vars.names`, `vars.conc`) and names from `metabolites(m)`.
  - It is independent of the mechanism type.
  - It sits on law type C (law.md L1) unchanged.
- **Codegen.**
  - `state!` reads `m.edges`, `bare_law`, `sol.deps` and the inactive mask; these become `m.edges[k]`, `_inactive_mask(am)` and the structured solve (law.md L4).
  - `program`, `horner`, `emit` and `body` are type-free.
  - F5's renumbering iterates `sol.deps` in column order.

## 7. Limits and reproduction

- Everything is on the slice. The package port reruns the key cost and the gate on its own law type and CI.
- The R5.5 sampling brackets the beam but is not identify's loss-based selection; max_param_count is not applied.
- Fit times are 12 easy synthetic fits: true model, 48 points, the package's `rate_equation`. Real fits are longer, so key shares are lower.
- Reproduce with `jl_research.sh`:
  - `pops5.jl <rxn> <mode> <quotas…>`, then `r55_table.py`, `r55_sample.py`, `r55_cost.py`;
  - `export5.jl small|terter`, `run_chunks.sh <in> <out> key|code <prefix> [+ver]`, `kt_summary.py`;
  - `fittime.jl`, `leafcap.jl`, `roots.jl`, `csvkey.jl`, `subproc.jl`, `bodydiff.jl`.
