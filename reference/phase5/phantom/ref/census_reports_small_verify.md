Denis, the census numbers hold. Every table cell recomputes exactly from the saved records, and one row (uni-bi, round 2) reproduces from scratch. Three findings change how to read it:

- **D and Z are exact only for the correct law.** For 694 D pairs, today's package fits a model with one more identifiable parameter (ROADMAP known issue 2). The rewrite's law would make them equal.
- **The key misses many same-law duplicates.** In uni-bi round 2, 122 of 158 candidates share their package-law family with an earlier fit; D counts 28. This includes 8–11% of the "OK" class.
- **Under your answer 2, W splits in two.** For islands and G1 cycles (2,257 candidates, including all 424 W in the main runs), the package law is the steady state of the package's own ODE. They fail only the stricter `:OnlyA`-limit requirement.

## Verdicts

| class | verdict | evidence |
|---|---|---|
| W | CONFIRMED membership; 2 report statements wrong; one sub-check incomplete | admission: 0 disagreements in 34,125 candidates; faces: 31/31 flagged fail, 0 false alarms outside G1 |
| D | CONFIRMED for the correct law; PARTIAL for today's package | same package law in 72/75 random pairs; np+1 pairs: package law larger, correct law equal (30/30) |
| P | CONFIRMED | rank 260/260 sampled plus 158/158 in the row; 81/92 of the census's same-law pairs hold on package laws |
| N | CONFIRMED | 0 of 1,888 laws grow without bound |
| Z | CONFIRMED for the correct law; REFUTED for today's package | 0/38 Z children have their parent's package law (36 strictly larger) |

## Methods (none use the census code)
- **Laws:** the compiled package `rate_equation`, with complex-step derivatives. I used Float64 first and Complex{BigFloat} (256 bits) where the gap was ambiguous.
- **Same law:** I fitted one law's log-parameters to the other's rates (Levenberg–Marquardt) at 2·np+12 random points, both ways, at two parameter draws. "Same" means both residuals are below 1e-8.
- **Oracle** (`moracle.jl`, written for this check):
  - Random constants consistent with thermodynamics and the kinetic groups, built in exact rationals; BigFloat solves at 512 bits.
  - Allosteric mechanisms at one subunit: the active network plus the package's inactive network, joined by a fast flip at the free form.
  - Faces: the RE limit exactly on the face (Λ = 1e30 and 1e40) under two RE-speed patterns (10^[−3,3]), against the interior limit (ε = 1e-25, Λ = 1e60). I also computed the long-time state reached from the free enzyme.
- **Admission:** re-implemented from its definition on the package's own types.

## W
- **Admission.** My code agrees on all 34,125 fitted allosteric candidate-records in 7 runs (the supplements repeat rounds 0–1).
  - Islands: 424 in ping-pong, 1,813 in ping-pong with 1,200 parents.
  - G1: 48 inside a component and 16 cycles in bi-bi with all parents; 14 and 4 in ping-pong with 1,200 parents.
- **`:OnlyA`-limit oracle.** Every form and step is present in the inactive state. `:OnlyA` plain bindings have their affinity divided by e^(60·e). Other `:OnlyA` steps keep a free inactive ratio, and their steady-state rates are scaled by 1e-24.
  - Islands: the rate depends on the free constant in 24/24. The spread between two draws at a positive point runs from 2e-4 to 0.91 (median 0.40).
  - The island-populated laws lie outside the package's family in 20/20 (best fit residual 2e-3 to 1.0). So the package law is one special case, not a reparametrization.
  - G1: no consistent family in 28/28.
  - Admitted: 84/84 have a unique limit equal to the package model to 1e-12. One of them, #19519, needed exponents taken from component heights.
- **Faces.**
  - The 31 flagged candidates all fail. In the 21 plain splits the law equals the interior limit to 1e-16. The exact-face limit depends on RE speeds and differs from the law by 1.4e-3 to 0.71 in relative terms.
  - Unflagged: 0 of 104 fail, and 1,608 of 1,609 agree in a larger sample. The one disagreement, bi-bi #19907, is a G1-inside candidate.
  - From the free enzyme, all 62 G1-inside candidates fail at the all-products face and at one single-product face, by 0.5% to 84%. The census flags 10 of them. All 20 G1 cycles pass.
  - Inhibitor at zero (not in the census): 230 candidates tested, 0 fail.
- **The package law matches the deletion-model steady state.** The deletion model is the package's own: non-`:OnlyA` steps on the forms reachable from the free form, a flip at the free form, and the enzyme starting free. Fitting the package law to oracle rates at positive points:
  - Admitted: 40/40.
  - G1: 20/20.
  - Islands: 20/20, though the 10-parameter #13220 converged on only 1 of 3 further draws.
- **Examples.** `E(res)` abbreviates `E(; residual = A − P)`.
  - Island, ping-pong #3273: `(E + A ⇌ EA, EQ + A ⇌ EAQ) ::EqualAI; (E + Q ⇌ EQ, EA + Q ⇌ EAQ) ::OnlyA; EA <--> EP(res) ::OnlyA; (EB(res) + P ⇌ EBP(res), E(res) + P ⇌ EP(res)) ::OnlyA; EB(res) ⇌ EQ ::OnlyA; (EP(res) + B <--> EBP(res), E(res) + B <--> EB(res)) ::EqualAI | n=1`. The package model gives −0.0207; the populated variant gives −0.0138 or −0.0108, depending on the free constant.
  - G1-inside with a face failure the census missed, bi-bi #18045: `(E + A <--> EA, EB + A <--> EPQ) ::OnlyA; (E + B <--> EB, EA + B <--> EPQ) ::OnlyA; (E + Q <--> EQ, EB + Q <--> EBQ) ::EqualAI; (EA + P ⇌ EAP, EQ + P ⇌ EPQ) ::OnlyA; EQ + B <--> EBQ ::EqualAI | n=1`. In the inactive state E(B) reaches E only through Q-bound forms. At Q = 0 it is cut off but keeps its weight in the law, which then misses the exact-face steady state by up to 46%.
  - Plain split, bi-bi #14463: `E + A ⇌ EA; (E + P <--> EP, EQ + P <--> EPQ); (E + Q <--> EQ, EP + Q <--> EPQ); EA + P ⇌ EAP; (EA + B <--> EPQ, EQ + B <--> EBQ); EP + A ⇌ EAP`. At all substrates zero, with the package's own constants, the law gives −0.0119211. The exact-face limit gives −0.0119197 or −0.0119131 under the two speed patterns.

## D
- **Cause, from the package types.** Siblings that differ only in tags: 32/34, 450/458, 713/717, 4,490/4,992 and 3,229/3,407. These match the census exactly.
- **Package laws, 75 random pairs.**
  - 72 are the same family. 6 of those needed 40 optimizer starts.
  - 70 of the 75 print different parameter names, which is why eq_hash misses them.
  - The other 3 are np+1 pairs.
- **np+1 pairs: 694 across the runs** (2, 8, 4, 502 and 178 in uni-bi, bi-bi, ping-pong and the two supplements).
  - In a sample of 30, the later candidate's package law is strictly larger in 28 and different at rank+1 in 2.
  - Its correct law (oracle) lies inside the earlier candidate's package family in 30/30, with residual ≤ 9e-12.
  - Example: uni-bi #213 `E + P ⇌ EP ::EqualAI; E + S <--> EPQ ::OnlyA; (E + Q <--> EQ, EP + Q <--> EPQ) ::OnlyA; EQ + P ⇌ EPQ ::EqualAI` (np 6, rank 5) against #50 `(E + P ⇌ EP, EQ + P ⇌ EPQ) ::EqualAI; E + S <--> EPQ ::OnlyA; (E + Q <--> EQ, EP + Q <--> EPQ) ::OnlyA` (np 5, rank 4).
- **Same law, different key (missed by the census).**
  - Uni-bi round 2: 122 of 158 candidates share an earlier fit's law, against 28 D. Of the 96 extra, 92 are P.
  - The other 4 have full rank and sit in the census's OK class. Three are certain and one probable (one of three draws stalled at 1.5e-5). Example: #127 `(E + P <--> EP, EQ + P <--> EPQ) ::OnlyA; (E + Q ⇌ EQ, EP + Q ⇌ EPQ) ::OnlyA; E + S ⇌ ES ::OnlyA; ES <--> EPQ ::OnlyA` has the law of #123, which differs only in `E + S ⇌ ES ::EqualAI`.
  - Round 1, census-OK, compared only with the parent and earlier same-rank siblings: 9 of 80 in bi-bi and 6 of 80 in ping-pong. This is a lower bound. 14 of the 15 are `to_allosteric` children with a plain sibling's law. Example: bi-bi #409 `E + A ⇌ EA ::OnlyA; E + P ⇌ EP ::EqualAI; (EA + B ⇌ EAB, EP + B ⇌ EBP) ::OnlyA; EAB <--> EPQ ::OnlyA; EP + Q ⇌ EPQ ::OnlyA | n=1` has the law of #395 `E + A ⇌ EA; E + P ⇌ EP; (EA + B ⇌ EAB, EP + B ⇌ EBP); EAB <--> EPQ; EP + Q <--> EPQ`.

## P
- **Ranks.** All 260 sampled candidates agree with the census (150 P and 110 not, in 11 runs).
  - A first Float64 pass disagreed on 7. All 28 ambiguous cases then agreed in BigFloat over 3 draws (kept ≥ 5e-8, dropped ≤ 3e-75).
- **The census's 92 "same law" pairs, re-tested on package laws:**
  - 66 are the same law.
  - 15 lie strictly inside the comparator at equal rank, so the fit is still wasted. Example: uni-uni #5 `E + P <--> EP ::OnlyA; E + S ⇌ ES ::EqualAI; ES <--> EP ::OnlyA | n=1` against the seed `E + P ⇌ EP; E + S ⇌ ES; ES <--> EP`.
  - 9 are strictly larger. Four of these are package-only, known issues 1 and 2, e.g. bi-bi #1663 `E + A <--> EA ::OnlyA; E + P ⇌ EP ::EqualAI; (EA + Q <--> EAQ, EP + Q <--> EPQ) ::EqualAI; EA + B ⇌ EPQ ::OnlyA | n=1`.
  - 2 differ; these may be optimizer misses.
  - The census's named examples (#4, #340, #3052) are confirmed.
- **Against the parent only (my own sample):** 85 of 99 L-phantom candidates are the same law, against 11 of 33 chain-flank and 0 of 16 flux-free SS candidates.

## N
- I scanned 1,888 laws for growth along joint and single-reactant rays, both directions, at 3 draws:
  - every fitted candidate in the uni-uni and uni-bi runs;
  - 250 random fits in each of bi-bi, ping-pong and the two supplements;
  - the 14 hits of the census's numeric probe.
- None grows. One slope of 0.25 (bi-bi #2349) saturates when followed to 1e30 in BigFloat.

## Z
- Every Z candidate is also D, in every run.
- On package laws, 36 of 38 Z children are strictly larger than their parent and 2 differ. So "the law equals the parent's" is true only for the correct law.

## Row recomputed: uni-bi, round 2
- I re-ran the pipeline with my own code and the public `expand_mechanisms`: 164 generated, 6 eq_hash duplicates, 158 fitted. These are the same candidates, in the same order, with the same parents.
- W 0, N 0 and P 124 all match the census. A first Float64 pass gave rank 4 for #119; BigFloat gives 6, as the census has.
- D (by law) is 122 against the census's 28. 26 of the 28 agree; #202 and #213 are np+1.
- Z (by law) is 101 against 2.
- The candidates in any class come to 128 against 124, so OK is 30 against 34.

## Your answer 2
- **Consistent with answer 2:** the 31 face splits and 62 G1-inside candidates. Their law is not the ODE steady state where initial-rate data sit, with all products at zero.
- **Matching the package's own ODE:** the islands (2,237) and G1 cycles (20). Their law is the deletion-model steady state at every point and face I tested. They are rejected only because the full two-conformation model, with the `:OnlyA` steps vanishing, has no unique limit (islands) or no consistent one (G1 cycles).
- Whether to keep this second group depends on which ODE you mean, so that is your call.

## Also missed or misstated
- The census's face check looks only at rapid-equilibrium segments. It misses inactive forms that are cut off through steady-state steps, which is 52 of the 62 G1-inside candidates.
  - It changes no count, because admission rejects these first.
  - In all 62, the offending binding shares its kinetic group with a fused binding. That is why the package's own guard lets them through.
- "Face split, all from `re_to_ss`" is wrong. 21 of 31 are; 10 are `to_allosteric` G1-inside children.
- "Ping-pong: islands only, 1,813 in the supplement" leaves out 18 G1 rejections there.
- `faces_of` labels the all-products face "P=0", the same as P alone. The "W faces" tallies merge the two.
- The report counts "852 eq_hash duplicates"; the data hold 830. The zero still holds.

## Not checked
- Laws with more than one subunit under the oracle: I used one subunit, with regulators at zero.
- Ter-bi, ter-ter and larger reactions; fits and losses; the census's cost figures.
- Same-law verdicts are numerical, so "different" may be an optimizer miss. The OK-duplicate figure is a lower bound.
- My `:OnlyA` family is one construction of the slice's limit.
- The census's P cause attribution and its 7 "no comparator" verdicts.

All Julia runs went through the lock with hard timeouts, about 27 minutes in total. The package (git clean) and `slice/src` are unchanged.

Everything is in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/small/verify/`:
- **Libraries:** `vlib.jl`, `moracle.jl`
- **Scripts:** `v1_tables.jl` through `v18_g1groups.jl`
- **Logs** (in `logs/`):
  - tables: `v1_tables.log`
  - row: `v7_row.log`, `v15_followup.log`
  - W: `v3_W.log`, `v9_facesample.log`, `v11_g1.log`, `v10_fit.log`, `v16_refit.log`, `v8_extra.log`, `v18_g1groups.log`
  - D and Z: `v5_DZ.log`, `v12_Dpkg.log`, `v17_okdups.log`
  - P: `v4_P.log`, `v4b_rank.log`, `v14_psample.log`
  - N: `v6_N.log`, `v13_growth.log`