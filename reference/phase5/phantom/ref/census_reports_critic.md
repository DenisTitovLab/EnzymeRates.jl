Denis, here is the synthesis. Three filters are worth adding: dropping islands, deduplicating tag sets, and dropping one-subunit L-phantom tag variants. Together they remove 21–47% of an expansion round's predicted fit time at 48 data points, and 14–39% at 480. Seed rounds lose under 2%. The other classes are too small to matter.

## Share of fit time by class

**Method.**
- I applied the fitcost track's pooled model to every fitted census candidate: t = F + 20·G·g, using the mean of G.
- I converted slice D terms to package law size with a log-log fit on the fitcost pools (R² 0.984).
- Each column excludes the columns to its left, in the order the filters would run. Shares are of predicted serial fit time per round.
- I exported the saved records under the lock, about 5 s per run; nothing was re-run.

| reaction, round | W | D, not W | L phantom, not W or D | other P | total | total, 480 points | census count share |
|---|---|---|---|---|---|---|---|
| bi-bi r1 | 0 | 9.3 | 11.6 | 0.5 | 21.4 | 13.9 | 32.2 |
| bi-bi r2 (all 2,752 parents) | 0.8* | 20.7 | 15.1 | 6.2 | 42.9 | 33.8 | 56.7 |
| ping-pong r1 | 12.2 | 7.7 | 9.6 | 1.2 | 30.8 | 22.3 | 38.5 |
| ping-pong r2 (1,200 parents) | 17.3 | 12.7 | 13.7 | 5.4 | 49.1 | 39.0 | 58.9 |
| ter-bi r1 | 9.2 | 14.2 | 7.5 | 1.2 | 32.1 | 26.7 | 37.3 |
| ter-bi r2 | 8.4 | 17.3 | 7.3 | 5.9 | 39.0 | 32.6 | 49.7 |
| ter-ter r1 | 18.2 | 18.0 | 2.3 | 4.3 | 42.9 | 36.6 | 47.2 |
| PC r1 | 31.5 | 13.7 | 1.4 | 5.0 | 51.6 | 45.3 | 56.5 |

\*The bi-bi W candidates are G1 and face cases; there are no islands in bi-bi. Elsewhere W is almost entirely islands.

- **Z and N.** Z is 2% or less and sits inside D. N is 0 everywhere.
- **Uni-uni regulator reactions.** They spend 8–46% of round time on chain-flank and flux-free P, but each round takes under a minute serially.
- **Scale.** Bi-bi round 2 from all parents is 1.6 h serial at 48 points and 6.9 h at 480.
- **Cost weighting.** Every share lands between the census's two proxies. For bi-bi round 2: 56.4% by parameter count, 42.9% by the model, 33.2% by D terms.

## Filters worth adding

A fitted candidate costs about 80 ms of fixed overhead (building the fitting problem, the first `loss!` compile, the fitter compile) plus 0.1–1.5 s of optimization. Every check below costs under 3% of that.

1. **Duplicate tag sets (D).** 90–99% of D are `to_allosteric` siblings that differ only in the tags of groups the inactive conformation cannot reach. The rest are Z.
   - **Now, in the package.** In `_expand_to_allosteric`, emit one K-type subset per pattern of inactive-reachable groups and their tags. It is exact: `_inactive_groups` keeps only steps between forms reachable from the free enzyme, so those tags change names only. The tested pairs have the same package law, 72 of 72. It costs microseconds.
   - **Reachability.** With all 2,752 bi-bi parents expanded, none of the 12,652 new round-2 laws was reachable only through a dropped sibling.
   - **Ping-pong is inconclusive.** Up to 44 keys are lost (7 full rank), all from `change_allo_state`. But only 1,200 of 3,220 parents were expanded, so the other paths may simply be missing.
   - **Later, with the rewrite.** Use the rename-proof key as identify's fit memo: duplicates inherit the fit and stay in the frontier. This costs 0.4–0.6 ms median, has no reachability risk, and also catches Z and plain renames. `_offer_cv!` would need to deduplicate by the key too.
   - **Order does not matter.** W and D never split a pair: in all six runs, both members of every D pair are W or neither is.
2. **Islands (`Slice.admit`).** 4 µs or less per candidate; 8–32% of round time in ping-pong and the ter reactions.
   - **Reachability is safe.** In ping-pong, every round-2 law found only through an island is itself W.
   - **Validity depends on the `:OnlyA` semantics (B1).** If `:OnlyA` steps are absent from the inactive conformation, as in today's package, the law is exact. If they vanish with a fixed ratio, as in your B1 note, the steady state depends on a hidden inactive constant. It did so in 24/24 islands tested in the small reactions and 48/48 in the large ones, and the package law is then one special case.
   - **Under your answer 2.** Dropping islands in enumeration is fine either way. Making them an error for hand-written mechanisms holds only under the second reading.
   - **Unchecked.** 30–40% of islands are neither D nor P, which is 2.5–9.9% of round time. Nobody tested whether any other candidate has their law.
3. **One-subunit L phantoms (P).** After W and D they take 9.6–21.9% of round time in bi-bi and ping-pong, but only 1.4–7.5% in the ter reactions, where most are already W or D.
   - **Same law elsewhere.** The census found a smaller candidate with the same law in 30/30 and 47/47 samples. The verifiers found the parent's law in 85/99 and 17/18.
   - **Cost of not emitting them** (the ROADMAP's planned item). In the complete bi-bi depth-2 run this loses at most 4 of 7,486 new full-rank laws, all `change_allo_state` children. It also loses 1,184 L-phantom keys that most likely duplicate smaller laws.
   - **It needs an exact detector.** Candidates are the slice rank (18 µs to 1.9 ms) or the 3.5 ms modular rank. The structural test has not been derived.
   - **Alternative.** Skip the fit and inherit the parent's loss. This keeps every law but mis-scores up to 14% of these candidates (one verifier's parent-only comparison found 14 of 99 whose law is not the parent's).

## Not worth a filter
- **N:** 0 everywhere.
- **Z:** 2% or less, all inside D. Today's package fits a strictly larger law for 36 of 38 tested Z children; that is ROADMAP known issue 2, so fix the derivation instead.
- **Face splits, traps and G1 cases:** 0.8% of time or less, so this is a correctness question only. The trap verdict is unreliable: 16 of 16 were false alarms.
- **v ≡ 0, RE turnover and 0/0 on a face:** 0 in enumerated output.
- **Other P** (flux-free steady-state steps, RE isomerization gauge, chain flank, NonequalAI copy, package-only): 0.3–7.9% of time. They often carry a larger family (18 of 32 flux-free samples), so count them by rank instead of dropping them.

## Overclaims, contradictions, unverified numbers
1. **Corrections the verifiers confirmed.**
   - D and Z are exact only for the slice law. All 694 np+1 pairs are the Z children.
   - The key misses same-law duplicates: in uni-bi round 2, 122 candidates share an earlier fit's law against 28 counted as D.
   - "Face splits all from `re_to_ss`" is wrong: 21 of 31 are.
   - The ping-pong supplement's "islands only" leaves out 18 G1 rejections.
   - The census says 852 eq_hash duplicates; the data hold 830.
   - In the large census, all 16 trap verdicts are false alarms.
   - "Package count exceeds slice count in 250" is 913 over all fitted children. The "one more dimension" is one to four.
   - The L-phantom share is 13.8–24.7%, not 14–27%.
   - In "0 of 23 gauge/other", 22 were never compared with anything, and PC (1,1955) is a counterexample.
   - D after the island rule is 16.3–24.8% of expansion fits, not 22–40%.
2. **Neither census weighted by fit cost.** Their "fit-cost weights" are sums of parameters and of D terms; the table above is the cost-model version.
3. **The census filter costs are compared with the wrong baseline.** They compare against identify's first pass (1.35 ms). What a filter saves is the full fixed cost (about 80 ms) plus the fit.
4. **"W already decided" holds for enumeration only.** For hand-written mechanisms it now depends on answer 2 and B1. Neither census separated one-step islands from multi-step ones, the distinction your B1 note draws.
5. **The rank timings are unverified and disagree 100-fold:** 1.9 ms in the small census (BigFloat) against 18 µs in the large one (Float64). Neither figure matters against a fit.
6. **The fitcost track had no verifier.**
   - Its coefficients reproduce when I re-run `pooled.py`, but that is a re-run, not a check.
   - Its sample was 48-point synthetic data, np 5–11, CMA-ES only, and it excluded W, D and Z candidates.
   - Generation counts correlate only 0.55 between its two datasets.
   - I applied its ter-bi factor to ter-ter and pyruvate carboxylase without calibration.
   - Its `maxtime` finding is correct: the time limit applies to each restart (`src/fitting.jl:266`), while identify's docstring says "max time per fit" (`src/identify_rate_equation.jl:118`).
7. **My numbers are bounds and model estimates.**
   - The reachability figures are upper bounds: the census credits each child to its first parent only.
   - Law size is extrapolated for laws above 5,244 operations.
   - The 480-point column assumes generation counts do not change.

## Still uncertain
- **Loss-based selection.** Both censuses expanded random parents. In identify, duplicates tie in loss, so they enter the beam together. `_offer_cv!` deduplicates by eq_hash only, so D siblings can fill several of the five cross-validation slots per parameter count. That cost is not in any number here.
- **Coverage gaps.** Depth 3 and beyond; round 2 for ter-ter and pyruvate carboxylase; pyruvate dehydrogenase; MWC with more than one subunit outside uni-uni; ping-pong reachability.
- **Real data sizes.** Shares fall 5–10 points from 48 to 480 points.

Files are in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/critic/`:
- `export_small.jl`, `export_large.jl`: export the census records
- `small.tsv`, `large.tsv`: the exported records
- `shares.py`, `shares2.py`: cost weighting
- `reach.py`, `reach2.py`, `reach3.py`: reachability bounds

The package and `slice/src` are unchanged.