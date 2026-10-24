Denis, here is the census for the small and medium reactions. It covers the six uni-uni reactions (plain, plus the ones with an inhibitor, oligomeric state or allosteric regulators), uni-bi, bi-bi and ping-pong bi-bi. The three regulator reactions run allosteric expansion with regulator moves. No fits were run. "Fitted" means the representative identify would fit.

## Bottom line

- **Today's filters remove very little.**
  - Complexity cap: 7 per bi-bi reaction (6 in round 1, 1 in round 2).
  - Parameter cap: 0. Compile errors: 0.
  - eq_hash: 6 in uni-bi round 2, 22 in bi-bi, 28 in ping-pong.
  - Seed screen: 7 of 264 ping-pong seeds, all condition C (a reactant drops out), replaced by 21 flip children. It removes nothing elsewhere.
- **W (your guard)** rejects only `:OnlyA` tag sets in the 60-parent rounds. Ping-pong: 331/3,220 in round 1 (10.3%) and 93/685 in round 2 (13.6%). All of them are islands on ping-pong residual forms.
  - Among all fitted candidates there is no v = 0, no RE turnover and no 0/0 on a face.
  - Face splits appear only at depth 2, in the larger supplementary runs: 24 of 17,290 bi-bi round-2 candidates (0.14%) and 7 of 10,967 ping-pong (0.06%).
- **D (same rename-proof key as an earlier fit, which eq_hash misses)** is 13–27% of fitted in rounds 1–2 of uni-bi, bi-bi and ping-pong. It is 0 in every uni-uni reaction.
  - About 90% or more of D candidates are `to_allosteric` siblings that differ only in the tags of groups the inactive conformation does not contain.
  - The rest are split children that keep the parent's law.
- **P (rank-deficient)** is 27–30% of bi-bi and ping-pong fits in round 1, and 43–46% in round 2.
  - It is dominated by the one-subunit L phantom.
  - In a stratified sample, all 30 L-phantom candidates have the same law as a smaller candidate in the run, usually their parent.
- **N (never saturates):** 0 of every fitted candidate in every run.
- **Z (same law as its parent):** 0.6–2.9%, all of them also D.
- **OK** (none of the above), bi-bi: 68% in round 1, 43–48% in round 2. Ping-pong: 62% in round 1, 41–42% in round 2.

## How it ran

- **Pipeline, in identify's order.**
  - Seeds come from `init_mechanisms`, or from `seed_mechanisms` when the reaction declares regulators (every declared regulator required, identify's default).
  - Then `_base_tier`, dedup by structural hash against everything seen, complexity cap 337, compile, parameter cap 20, `rate_equation_string` and eq_hash. One representative is kept per eq_hash not fitted earlier in the run.
  - Expansion repeats `expand_mechanisms` move by move, so each child carries its move label. It matched `expand_mechanisms` on up to 20 parents per reaction.
- **Rounds.**
  - Round 1 expands every round-0 entry; identify would select parents by loss.
  - Round 2 expands a seeded sample of at most 60 round-1 entries.
  - Supplements: round 2 from all 2,752 bi-bi parents, and from 1,200 of the 3,220 ping-pong parents. I stopped the full ping-pong run because memory was running out (peak RSS 4.3 GB on the bi-bi run).
- **W.**
  - RE turnover, using the guard probe's `re_turnover`.
  - v ≡ 0 and 0/0 on a face, using the support test on the combined law, generalized to every face.
  - Bottomless segments, splits and traps, using the guard's `face_check` on the full graph and on the inactive conformation.
  - `:OnlyA` admission, using `Slice.admit`.
  - Faces: all products at zero, all substrates at zero, and each single reactant at zero.
- **D.** The slice's `canonical` key, with concentrations labelled by name rather than position (the skeptic's note 1).
- **P.**
  - The rank of ∂v/∂log θ of the package's own law (what identify fits) at np+6 random points.
  - Gradients are exact forward-mode through the package's generated Exprs. Ambiguous gaps were redone in BigFloat.
  - Cross-check: wherever the package and slice parameter counts agree (about 41,000 rank comparisons over all runs), the ranks agree too.
  - This is absolute mode. In the default relative mode every law loses one more direction (the overall scale), which the kcat rescale handles. I did not count it.
- **N.** Condition V (`_has_vmax`) in both directions.
  - A numeric probe (v at 10⁹ against 10⁶) flagged 14 extra supplementary candidates.
  - Re-run in BigFloat, all 14 saturate, so the probe was wrong (`logs/unbounded.log`).
- **Validation.** The classifiers fire on the guard probe's known cases (`logs/validate.log`):
  - 34 of 34 flagged enumerator children (splits) and the hand-built LDH-like trap;
  - both RE-turnover cases and both v ≡ 0 cases;
  - 4 mutants whose D vanishes on a face;
  - a uni-uni law that never saturates.

## Tables (main census; % of fitted)

Columns: gen = candidates after first-parent dedup; seen = already produced earlier in the run; "any" = in at least one class.

| reaction | round | gen | seen | cx>337 | eq_hash | **fitted** | W | D | P | N | Z | any | OK |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| uni_uni_rxn | 0/1/2 | 1/3/10 | 0 | 0 | 0 | **1/3/10** | 0 | 0 | 0/3/10 | 0 | 0 | 0/100%/100% | 1/0/0 |
| uni_uni_with_reg | 0/1/2 | 1/7/26 | 0 | 0 | 0 | **1/7/26** | 0 | 0 | 0/7/26 | 0 | 0 | 0/100%/100% | 1/0/0 |
| uni_uni_allo (n=2) | 0/1/2 | 1/3/10 | 0 | 0 | 0 | **1/3/10** | 0 | 0 | 0/0/2 | 0 | 0 | 0/0/20% | 1/3/8 |
| uni_uni_allo_reg | 0/1/2 | 8/30/47 | 0 | 0 | 0 | **8/30/47** | 0 | 0 | 0/4/22 | 0 | 0 | 0/13%/47% | 8/26/25 |
| uni_uni_allo_2reg | 0/1/2 | 16/92/224 | 0/0/10 | 0 | 0 | **16/92/214** | 0 | 0 | 0/8/56 | 0 | 0 | 0/9%/26% | 16/84/158 |
| uni_uni_reg_and_inhibitor | 0/1/2 | 16/66/114 | 0 | 0 | 0 | **16/66/114** | 0 | 0 | 0/8/47 | 0 | 0 | 0/12%/41% | 16/58/67 |
| uni_bi | 0 | 7 | 0 | 0 | 0 | **7** | 0 | 0 | 0 | 0 | 0 | 0 | 7 |
| uni_bi | 1 | 44 | 0 | 0 | 0 | **44** | 0 | 6 (13.6%) | 25 (56.8%) | 0 | 0 | 56.8% | 19 |
| uni_bi | 2 (all 44 parents) | 164 | 0 | 0 | 6 | **158** | 0 | 28 (17.7%) | 124 (78.5%) | 0 | 2 | 78.5% | 34 |
| bi_bi | 0 | 239 | 0 | 0 | 0 | **239** | 0 | 0 | 0 | 0 | 0 | 0 | 239 |
| bi_bi | 1 | 2,758 | 0 | 6 | 0 | **2,752** | 0 | 354 (12.9%) | 736 (26.7%) | 0 | 0 | 886 (32.2%) | 1,866 |
| bi_bi | 2 (60 parents) | 663 | 2 | 1 | 22 | **638** | 0 | 104 (16.3%) | 283 (44.4%) | 0 | 8 | 335 (52.5%) | 303 |
| bi_bi_pp | 0 | 278 | 0 | 0 | 0 | **278** | 0 | 0 | 5 (1.8%) | 0 | 0 | 1.8% | 273 |
| bi_bi_pp | 1 | 3,226 | 0 | 6 | 0 | **3,220** | 331 (10.3%) | 563 (17.5%) | 952 (29.6%) | 0 | 0 | 1,241 (38.5%) | 1,979 |
| bi_bi_pp | 2 (60 parents) | 716 | 2 | 1 | 28 | **685** | 93 (13.6%) | 154 (22.5%) | 295 (43.1%) | 0 | 4 | 399 (58.2%) | 286 |

**Supplements, round 2:**

| reaction | parents | gen | seen | cx>337 | eq_hash | **fitted** | W | D | P | Z | any | OK |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bi_bi | all 2,752 | 17,754 | 28 | 58 | 378 | **17,290** | 80 (0.5%) | 4,638 (26.8%) | 7,680 (44.4%) | 502 (2.9%) | 9,804 (56.7%) | 7,486 (43.3%) |
| bi_bi_pp | 1,200 | 11,429 | 20 | 46 | 396 | **10,967** | 1,505 (13.7%) | 2,844 (25.9%) | 5,054 (46.1%) | 178 (1.6%) | 6,464 (58.9%) | 4,503 (41.1%) |

**Overlaps:**
- bi_bi round 1: P only 532, D+P 204, D only 150.
- bi_bi_pp round 1: P 556, D+P 204, D 150, W+D+P 129, W+D 80, W+P 63, W 59.

**Weighting by fit cost.** The "any" candidates carry this share of each round's fitted-parameter sum and D-term sum:

| run | parameters | D terms |
|---|---|---|
| bi_bi r1 | 31.8% | 12.1% |
| bi_bi r2 (60 parents) | 53.0% | 29.3% |
| bi_bi r2 (all parents) | 56.4% | 33.2% |
| bi_bi_pp r1 | 39.8% | 17.4% |
| bi_bi_pp r2 (60 parents) | 59.6% | 35.6% |
| bi_bi_pp r2 (1,200 parents) | 60.0% | 36.2% |
| uni_bi r1 | 56.8% | 46.4% |
| uni_bi r2 | 78.5% | 66.5% |

P and D laws are small. Median D terms in bi_bi round 1: P 11, D 10, OK 18 (p90 125). Per-candidate `:np` and `:dterms` are in the data files.

## Classes in detail

**W.**
- Ping-pong: `:OnlyA` islands only, 424 in the main runs and 1,813 in the supplement, all on residual forms. Example #3048 (r1, `to_allosteric`): `E + A <--> EP_res ::OnlyA; (E + Q <--> EQ, EQ + Q <--> EBQ_res) ::OnlyA; (EQ + B ⇌ EBQ, E_res + B ⇌ EQ) ::OnlyA; E_res + P ⇌ EP_res ::OnlyA | n=1`. Admission message: "E(; residual = A - P) reach the free enzyme only through :OnlyA steps that are not bindings toward them".
- Bi-bi full round 2: G1 rejections the package's own guard passes.
  - 48 "arc inside a component", e.g. #17966: `(E + A <--> EA, EB + A <--> EPQ) ::OnlyA; (E + B ⇌ EB, EA + B ⇌ EPQ) ::OnlyA; E + Q ⇌ EQ ::EqualAI; (EA + P <--> EAP, EQ + P <--> EPQ) ::OnlyA; EB + Q ⇌ EBQ ::EqualAI; EQ + B ⇌ EBQ ::EqualAI | n=1`.
  - 16 "arcs close a cycle", e.g. #17968.
- Face split, all from `re_to_ss`, e.g. #14463 (split at S=0 and at A=0): `E + A ⇌ EA; (E + P <--> EP, EQ + P <--> EPQ); (E + Q <--> EQ, EP + Q <--> EPQ); EA + P ⇌ EAP; (EA + B <--> EPQ, EQ + B <--> EBQ); EP + A ⇌ EAP`.

**D.**
- Cause: of all D candidates, uni_bi 32/34, bi_bi 450/458, pp 713/717, bi_bi full 4,490/4,992 and pp 1,200-parent 3,229/3,407 differ from the earlier fit only in the tags of groups outside the inactive conformation. The remainder are the Z split children.
- Example (uni_bi): #13 `E + Q ⇌ EQ ::OnlyA; E + S ⇌ ES ::EqualAI; EQ + P ⇌ EPQ ::OnlyA; ES <--> EPQ ::OnlyA | n=1` has the same key as #9, which tags `EQ + P ⇌ EPQ ::EqualAI`.
- Example (pp): #291, #295 and #299 all share #287's key.
- Spot check: 90 pairs (30 each in uni-bi, bi-bi, pp) evaluate identically under the key's constant correspondence (worst relative difference 0).
- Reverse direction: 0 of the 852 eq_hash duplicates have a different key from their representative.

**P.**
- Who is affected:
  - bi_bi: 1,015 of 1,019 P are allosteric at n=1. Of the 2,389 fitted n=1 candidates, 1,015 are P.
  - bi_bi_pp: 1,218 of 1,252 P are allosteric at n=1.
  - The oligomeric_state-2 uni-uni reactions have no L phantom.
- Deficit units by cause, bi_bi full round 2 (greedy attribution, a column is dropped when the rest span it): L phantom (n=1) 6,829, flux-free SS step 844, NonequalAI copy 812, chain flank 256, other 75. The package alone adds 192 more, from the ROADMAP's known issues 1 and 2.
- RE isomerization gauge: 4 ping-pong seed flips in round 0, e.g. #273 `(E + A ⇌ EA, EP_res + A ⇌ EAP_res); E + Q <--> EQ; (EA + P ⇌ EAP, E_res + P ⇌ EP_res); EA <--> EP_res; EB_res ⇌ EQ; E_res + B ⇌ EB_res`.
- Same-law sample. I took 131 P candidates, at most 6 per sub-class and reaction, and tested equality of coefficient images both ways: Levenberg–Marquardt on log-coefficients, at 2 points each way, residual below 1e-8.

| primary sub-class | sampled | same law as a smaller candidate | a strictly larger family | no smaller candidate with the same support |
|---|---|---|---|---|
| L phantom | 30 | 30 | 0 | 0 |
| chain flank | 29 | 27 | 2 | 0 |
| NonequalAI copy | 11 | 11 | 0 | 0 |
| other | 11 | 10 | 1 | 0 |
| flux-free SS step | 32 | 9 | 18 | 5 |
| RE isomerization gauge | 6 | 1 | 3 | 2 |
| package-only phantom | 12 | 4 | 8 | 0 |

- Examples of "same law":
  - uni_uni #2 `E + P ⇌ EP ::OnlyA; E + S ⇌ ES ::EqualAI; ES <--> EP ::OnlyA | n=1` (np 4, rank 3) has the seed's law.
  - bi_bi #340 has the law of its parent #6.
  - bi_bi #3052 (chain flank) has the law of its parent #345.
- Example of a larger family (flux-free SS step): bi_bi #3483 `E + A ⇌ EA ::OnlyA; (E + P ⇌ EP, EQ + P ⇌ EPQ) ::OnlyA; (E + Q <--> EQ, EP + Q <--> EPQ) ::NonequalAI; (EA + B <--> EPQ, EQ + B <--> EBQ) ::OnlyA | n=1`.

**N.** None.
- The 7 seeds the seed screen removes are all condition C, not V. Example: `E + A ⇌ EA; E + Q ⇌ EQ; EA <--> EP_res; EB_res ⇌ EQ; E_res + B ⇌ EB_res; E_res + P ⇌ EP_res`.
- Condition C never occurs in any fitted candidate.
- The kcat precondition fails nowhere, so no fit would end as a FitFailure through the rescale.

**Z.** Split children of allosteric parents. The package counts one more parameter (known issue 2), but the law equals the parent's. Example: pp #3618, a split of #587.

## Cost per candidate

Bi-bi full run, 20,659 candidates, median / p95.

| step | median | p95 |
|---|---|---|
| identify's own pre-fit step (compile, `fitted_params`, render, eq_hash) | 1.35 ms | 4.3 ms |
| slice conversion | 31 µs | |
| slice solve + law | 0.23 ms | 1.0 ms |
| key | 0.18 ms | 1.9 ms |
| W: support test, all faces | 28 µs | 55 µs |
| W: RE turnover + face split/trap, all faces + inactive conformation | 31 µs | 43 µs |
| W: admission | 4 µs | |
| N: condition V both ways + C | 15 µs | |
| P: rank on the slice law (exact BigFloat log-gradients) | 1.9 ms | 18.7 ms |
| P: rank on the package law | 0.31 ms | 5.5 ms |

- The package-law rank also needs about 9.6 ms to rebuild the package's Exprs.
- Maxima reach 0.8–0.9 s, probably GC pauses (I did not check).
- The same-law test is far costlier: about 0.4 s per sampled candidate, plus a search for comparators.
- I did not measure fit cost.

## What looks worth filtering

1. **D: yes.**
   - It is 13–27% of uni-bi, bi-bi and ping-pong fits, it is exact, and it costs about 0.2 ms.
   - Two places to do it: dedupe by key in identify, or emit one tag set per inactive-equivalence class in `_expand_to_allosteric`, which removes about 90% or more at the source.
2. **W `:OnlyA`: already decided.** It is 10–14% of ping-pong fits and costs 4 µs. The face checks (0–0.14%) and v = 0 / RE turnover (0) cost almost nothing but buy nothing in enumeration.
3. **P: worth acting on, but not by dropping candidates blindly.**
   - L-phantom, chain-flank, NonequalAI-copy and "other" deficits were almost always the same law as a smaller candidate: 78 of 81 sampled. Their fits are wasted.
   - Flux-free SS, RE-gauge and package-only deficits usually carry a larger family with one phantom parameter.
   - So skip or count by rank (about 2 ms), but keep the n=1 tag children as parents. Some of their children are new families (the "larger family" cases above). I did not check how many.
4. **N and Z: no.** N is 0, and every Z is already D.

## Not checked

- No fits: no losses, no beam selection, no data-dependent failures.
- No mass-action oracle on census candidates.
- Faces with modifiers at zero.
- MWC laws with n > 1 outside the uni-uni regulator reactions.
- Ter-bi, ter-ter and the larger reactions.
- The candidates `_seed_candidate_screen` rejects inside `init_mechanisms`.
- Ranks were taken at one parameter draw per candidate.
- The same-law verdicts are numerical, on a sample.
- Cause attribution for P is greedy.

Julia time was about 35 minutes. Every run went through the lock, and neither the package nor `slice/src` changed.

Everything is in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/small/`:
- **Scripts:** `lib.jl`, `census.jl`, `report.jl`, `tables.jl`, `p_sample.jl`, `validate.jl`, `sanity.jl`, `wexamples.jl`, `unbounded.jl`
- **Data:** `data/<reaction>.jls` (records) and `data/<reaction>.mechs.jls`; supplements `data/bi_bi_rxn.r2_100000.*`, `data/bi_bi_pp_rxn.r2_1200.*`
- **Logs** (in `logs/`):
  - `tables.log`: funnels, classes, weights, D causes, admission kinds, degenerate seeds
  - `report_all.log` and `report_supp.log`: sub-classes, examples, costs
  - `p_sample.log`: same-law sample and D spot check
  - `validate.log`, `sanity_all.log`, `wexamples.log`, `unbounded.log`
  - `run_all.log`, `run_r2all.log`, `run_pp_r2_1200.log`