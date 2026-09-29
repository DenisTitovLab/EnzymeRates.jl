# Cross-validation of engine A, engine B and the package export

**Verdict.** The three sources agree everywhere they can be compared. No engine bug turned up, so no engine was changed and no CHANGELOG was written. No package discrepancy turned up either.

| comparison | compared | disagreements |
|---|---|---|
| engine A vs engine B on the curated suite (fitted, 4-5 rank paths, supports, exact law values) | 69 valid + 8 invalid | 0 |
| engine A vs engine B on random fuzz mechanisms (same checks) | 384 valid + 116 invalid | 0 |
| package pkg_fitted vs engine A fitted | 12,556 | 0 |
| package pkg_fitted vs engine B fitted | 12,556 | 0 |
| package fitted-name set is a valid free set | 12,556 | 0 |
| package rank (pkg_rank.jsonl) vs engine A rank | 250 | 0 |
| engine A fast rank vs engine A symbolic rank | 12,556 | 0 |
| engine A fast rank vs engine B numeric rank | 12,556 | 0 |

## 1. Suite (cv/suite/, 77 mechanisms)

The suite holds:
- engine A's validation set: 9 suite cases, 14 extras and 4 invalid;
- engine B's validation set: 9 suite cases, 5 valid extras and 4 invalid;
- 16 new hand-written mechanisms (n01-n11, listed below);
- 16 export picks across R1-R6, levels 0-2 and all three moves, including 6 ping-pong and 7 inhibitor mechanisms, all with a pkg_rank.

**Method.** Both engines run with the symbolic law. compare.py then does the following:
1. It builds one random rational point of all group parameters that satisfies engine A's constraints. Free values are perfect M-th powers, M = lcm(D_A, L_B), so every root either law takes stays rational.
2. It checks engine B's own dependent relations exactly at that point.
3. It evaluates v exactly at 5 random rational concentration points in four ways: A's symbolic law, A's linear-solve evaluator over Q, B's brute-force RE-limit point law (full mass action with symbolic Lambda), and B's symbolic law. All four must be equal as rationals.

A negative control shows the comparison can fail: scaling one free parameter, or one dependent parameter alone, breaks the equality.

Ranks come from four independent paths:
- A symbolic, exact over Q;
- A fast, over GF(2^61-1);
- B numeric, mpmath at 150 digits;
- the exact coefficient-map rank of B's law. This path is skipped only for the sqrt-exponent case.

**Result.**
- All 69 valid mechanisms agree on fitted count, every rank path, supports and exact law values.
- Engine B's dependent relations hold at engine A's point in 69 of 69.
- All 8 invalid mechanisms raise in both engines with the same error class: infinite rate 2, stoichiometry 2, Keq forced 2, zero rate 1, disconnected 1.
- Engine B's smallest singular-value gap is 10^54.7.

New cases (fitted/rank, identical in every path):

| case | fitted/rank | note |
|---|---|---|
| n01 ordered bi-bi, all SS, two central complexes | 9/7 | Segel IX-87 is 7/7 |
| n02, n02b merged X with mixed RE/SS | 5/5, 5/5 | |
| n03 Theorell-Chance with RE outer bindings | 3/3 | |
| n04 ping-pong, six steps, all SS | 11/6 | the four-step Segel IX-140 is 7/6 |
| n05 ping-pong RE seed | 7/6 | |
| n05b ping-pong RE seed + dead ends | 9/8 | |
| n06 all-SS random bi-bi + one central isomerization | 17/15 | 15/15 without the isomerization |
| n07 random bi-bi with merged central complex | 9/9 | |
| n08 uni-uni RE with S dead-end at E | 4/3 | K_cat and K_dead enter only as kf*K_cat and K_cat+K_dead |
| n08b uni-uni SS with S dead-end at E | 6/3 | |
| n09a / n09b context shared vs split | 5/5 vs 6/6 | the split adds one identifiable parameter |
| n09c split inside a Wegscheider cycle | 5/5 | the cycle forces Ka1 = Ka2, so the split adds no freedom |
| n10 RE Wegscheider tie through a shared group | 5/5 | K3 = K2 is forced |
| n11 SS Wegscheider tie through a shared group | 8/8 | |

## 2. Random differential test

fuzz_gen.py produces random mechanisms with these features:
- uni-uni, sequential bi-bi (one or two routes sharing forms) and ping-pong;
- Theorell-Chance and merged-complex folding, and central isomerizations;
- 0-3 dead-ends, some binding two copies in one step;
- parallel steps, random step directions, random RE/SS types and group sharing.

The "wild" set also merges random same-type groups, which produces forced Keq, fractional exponents and cross-metabolite ties.

| set | total | valid | valid, all checks agree | valid, with phantoms | both raise, same class | one raises |
|---|---|---|---|---|---|---|
| fuzz | 300 | 235 | 235 | 209 | 65 (infinite rate) | 0 |
| fuzzw | 200 | 149 | 149 | 98 | 48 infinite rate + 3 Keq forced | 0 |

Engine B's smallest gap is 10^52.0 on fuzz and 10^49.3 on fuzzw, with no ambiguous rank.

## 3. Export (12,556 mechanisms)

- **Engine A fast rank.** It ran on all 12,556 in 14.7 s total, 1.1 ms mean and 6.3 ms max per mechanism.
- **Engine A symbolic path.** It ran on all 12,556 and matches the fast rank everywhere.
- **Engine B.** Its fitted count and numeric rank ran on all 12,556 at 0.28 s mean per mechanism. The smallest gap was 10^52.9, with 0 ambiguous ranks and 0 errors.
- **Fitted counts.** pkg_fitted equals engine A and engine B on 12,556 of 12,556.
- **Package free sets.** pkg_freeset_check.py maps every package fitted-parameter name (K_, kon_/koff_, k_a_to_b, Kiso_, <met>inh) to one group parameter. In 12,556 of 12,556, the non-fitted parameters are all fixed by the fitted ones and Keq, so the package's fitted set is a valid free parameterization everywhere. In 8,458 it is identical to engine A's pivot choice.
- **Package ranks.** The package's rank equals engine A's on 250 of 250. In that sample 46 mechanisms have phantoms, 62 phantoms in total, and both engines confirm all of them.
- **Hand check of R4_00061**, the package's own level-0 phantom (ping-pong with an RE second chemistry step, fitted 6). With the RE segment weights times B:
  - D = B + K1 AB + K2(1+1/K5) BQ + K1K3 ABP + K2/(K5K6) Q + K2K3/(K5K6) PQ + K1K2K3/(K5K6) APQ, and N = kf4 K1 (AB - PQ/Keq).
  - The identifiable combinations are K1, K3, kf4, K2(1+1/K5) and K2/(K5K6), so the rank is 5, as all sources report.
  - Engine A's supports and its expression kr4 = K1K5K6kf4/(K2K3Keq) match this derivation.

The full disagreement table is results/export_disagreements.tsv. It covers any fitted count, any rank or the free-set check disagreeing, and it has no rows (header only).

## 4. Fixes made

None, because no comparison disagreed. I also re-read both engines' core math:
- engine A: Cha weights, King-Altman in-arborescence orientation, telescoped omega weights, the reverse-mode gradient's handling of the normalization row, and the rank[dcoef|coef]-1 kernel argument;
- engine B: cycle de-duplication, unit-pivot parameterization, bordered Bareiss solve and the Lambda-limit coefficient ratio.

I found no defect.

## 5. Remaining discrepancies

None. What this cross-validation does not establish:
- **Package rank beyond the sample.** The package rank is known only on the 250 sampled mechanisms, one random draw each (smallest gap 10^44.6). On the rest, the valid free set proves the package's fitted parameters reparameterize the engines' parameters, but the package rate law itself was not evaluated (that needs Julia).
- **Common specification.** The two engines share no code but were written from the same specification.
- **Random-point ranks.** All ranks are exact or numeric at random points (Schwartz-Zippel), not symbolic proofs of generic rank.

The "move" label records the move that first discovered a mechanism, so a mechanism reachable by two moves counts under the first one. For example, R5 level-2 re_to_ss equals R4's 842.

## 6. Phantom table

See the phantom_table field. Tables split by move and by ping-pong vs sequential are in results/tables.md. Some of those rows:
- R4 L2: re_to_ss 842 / 295 with phantoms / 479 phantoms; split 546 / 88 / 100.
- R6 L2: dead_end 3718 / 585 / 611.
- Ping-pong vs sequential: R4 L2 171 / 77 / 105 vs 1217 / 306 / 474.

## Reproducing (from cv/, Python only)

build_suite.py; run_engine.py a|b [dir]; compare.py [dir]; negative_controls.py; fuzz_gen.py 300 20260928 fuzz; fuzz_gen.py 200 777 fuzzw wild; run_export_a.py; run_export_a_sym.py; run_export_b.py 3; pkg_freeset_check.py; tables.py; summarize.py
