Denis, the census holds up. My own driver reproduced its pipeline and every number in its W and D rows. Each class passed an independent test. There is one real error, inside W: all 16 "trap" verdicts on the all-products or all-substrates face are false alarms. The 17 face-flagged mechanisms are still W, because each fails on a single-reactant face. Several numbers in the census's prose are also off; they are listed below.

## Verdicts

| Class | Verdict | Independent test |
|---|---|---|
| **W, island rule** (1,273) | **CONFIRMED** | **Oracle with the `:OnlyA` steps vanishing:** <br>– 48 of 48 flagged give a different rate depending on how the steps vanish. <br>– 48 of 48 unflagged give one limit, equal to the law. |
| **W, faces** (17) | **CONFIRMED** as a class; the trap verdicts are **REFUTED** | **Exact-face oracle on every face of all 12,061 fitted laws:** <br>– Exactly these 17 fail. <br>– The failure is always on the single-reactant face the "split" verdict names, never on the full face a "trap" verdict names. |
| **D** (2,046) | **CONFIRMED** | **Law values under the key's renaming:** <br>– All 2,046 pairs agree to 2e-76. <br>– 36 pairs checked by the oracle agree to 2e-28. |
| **P** | **CONFIRMED** | **Rank at new points (512-bit) and from the oracle's Jacobian:** <br>– 48 of 48 flagged are deficient. <br>– 70 of 70 unflagged are full rank. |
| **N** (0) | **CONFIRMED** | Exact leading exponents of all 12,061 laws along 200+ growth directions: none grows without bound. |
| **Z** (18) | **CONFIRMED** | 18 of 18 have exactly the parent's law; the oracle agrees to 5e-27. |

The "oracle" is a BigFloat mass-action steady state I wrote in `vlib.jl`. Rapid-equilibrium steps run at speed Λ × a random factor. Allosteric mechanisms use the n=1 network with both conformations. It shares no code with the slice's law derivation. Two caveats:
- Its constants come from the slice's `all_parameters`; the inactive constants outside the law are solved from the slice's cycle rows.
- Every set of constants satisfied all cycle rows (residual ≤ 1e-153).

At one random positive point, the slice law equals the oracle for all 12,061 fitted candidates (largest difference 3e-17 at Λ = 1e22). None diverges and none has more than one steady state.

## One row, recomputed end to end
`v5_e2e.jl` is my own driver built on package calls, with the census's seeds and sampling. All seven columns match the census exactly in these rows:
- Input
- Seed screen (65/179, 95/243, 174/430)
- Complexity filter
- eq_hash dedup (8)
- Fitted (allosteric)

The fitted sets match by structural hash: 12,061 of 12,061.

On that recomputed set:
- **W row: 0 | 105 | 39 | 8 | 434 | 7 | 697.** This is my own island rule, written on the package mechanisms rather than the slice, plus the face results from the oracle.
- **D row: 0 | 288 | 172 | 0 | 786 | 0 | 800.**

Both are exact.

## Evidence

**W, island rule**
- My rebuilt island rule agrees with `admit` on all 4,791 fitted allosteric records (1,273 rejected, 3,518 accepted). Every rejection carries the island message.
- **Vanishing test:**
  - Every inactive form is kept and the `:OnlyA` steps vanish with thermodynamically consistent constants: a binding by kon → 0, any other step by both rates → 0.
  - Components the free enzyme cannot reach through `:OnlyA` bindings get an energy offset: above the free enzyme (A, and A2 with other free constants), level with it (B), or below it (C).
- **Flagged, 48 sampled (12 per reaction-round):**
  - A and A2 equal the law to 1e-8 or better.
  - B differs from the law by 2.6e-8 to 0.99.
  - C drives the rate below 1e-30.
  - So the law the slice and today's package derive is just the limit in which the stranded forms vanish.
- **Unflagged, 48 sampled:** all patterns equal the law to 2.5e-10.
- **Example, ter-bi (1,1157):**
  - Steps: `E + P <--> E(; residual = A + D - Q) + B [OnlyA]; (E + A <--> E(A), E(Q; residual = A + D - Q) + A <--> E(A, Q; residual = A + D - Q)) [OnlyA]; (E(; residual = A + D - Q) + Q <--> E(Q; residual = A + D - Q), E(A) + Q <--> E(A, Q), E(Q; residual = A + D - Q) + Q <--> E(A, D, Q)) [OnlyA]; (E(A) + D ⇌ E(Q; residual = A + D - Q), E(Q; residual = A + D - Q) + D ⇌ E(D, Q; residual = A + D - Q), E(A, Q) + D ⇌ E(A, D, Q)) [OnlyA]`
  - Law and A: 0.00188262. B: 0.00184748. B with other vanishing speeds: 0.0017858. C: 8e-41.
- **Example, PC (1,1057):** law −0.0886, B −0.00186.

**W, faces**
- For the 17 flagged, I ran 3 parameter points × 20 speed orderings in 10^[−6, 6].
- **16 records with a full-face trap verdict:** on that face the law equals the oracle to 5e-23. The trap is fed only from a form that itself holds the zeroed product, so it is never entered. The round-4 probe had already flagged this weakness.
- **Their single-reactant "split" face:** the limit is unique (spread over orderings ≤ 4e-23) but differs from the law by 3e-5 to 0.95.
- **Only ter-bi (2,460) depends on the speeds:** the spread is up to 28% on substrates=0 and up to 0.6% on A=0.
- **Example, ter-ter (0,473):**
  - Steps: `(E + B ⇌ E(B), E(R) + B ⇌ E(B, R), E(A, R) + B ⇌ E(A, B, R)); (E + P ⇌ E(P), E(R) + P ⇌ E(P, R), E(A, R) + P ⇌ E(A, P, R)); (E + R ⇌ E(R), E(B) + R ⇌ E(B, R), E(P) + R ⇌ E(P, R), E(A, B) + R ⇌ E(A, B, R), E(A, P) + R ⇌ E(A, P, R)); (E(B) + A <--> E(A, B), E(P) + A <--> E(A, P), E(R) + A <--> E(A, R), E(B, R) + A <--> E(A, B, R), E(P, R) + A <--> E(A, P, R)); (E(B) + Q <--> E(B, Q), E(P, R) + Q <--> E(P, Q, R)); (E(P) + D ⇌ E(D, P), E(A, B) + D ⇌ E(P, Q, R))`
  - The census says "trap" at P=0: the law is correct on the all-products face (to 2e-25).
  - The census says "split" at R=0: the law is off by 6e-3 to 0.33 there.
- **Side issue:** in ter-bi and ter-ter, "P=0" names both the all-products face and the single-product face. That is why (0,652) lists "P=0" twice.

**D**
- **Renaming-blind fingerprint.** I evaluated every fitted law with all its fitted constants set to one value, which no renaming can change. Every key duplicate also matches by fingerprint.
- **507 records match an earlier fingerprint but not its key.** I fitted 86 of these pairs both ways with Levenberg–Marquardt. All 86 are different laws: the worse direction stalls at a residual of 1.2e-3 or more.
- **Pair structure.** 2,025 pairs are same-round allosteric siblings with identical steps whose tags differ only in groups outside both inactive conformations:

  | Groups whose tag differs | Pairs |
  |---|---|
  | 1 | 1,369 |
  | 2 | 508 |
  | 3 | 132 |
  | 4 | 16 |

  The census says "a group"; it is one to four. The other 21 pairs are the 18 Z and three plain sibling pairs.

**P**
- **New points, central differences in 512-bit arithmetic:**
  - Flagged: dropped directions ≤ 1e-81, kept ≥ 8.5e-9.
  - 21 borderline records (census weakest direction < 1e-6): 19 are full rank (kept ≥ 3.4e-8) and 2 are deficient, matching the census.
  - 16 records deficient only in the package law: the slice law is full rank.
- **Oracle Jacobian:** 36 records, rank equal to the slice rank in all 36 (kept ≥ 3.7e-5, dropped ≤ 3e-29).
- **Package law in BigFloat:** 37 records, rank equal to the census package rank in all 37.
- **Sub-kinds by deleting one column:** L adds no rank in 15 of 15 L-phantom samples; the RE-isomerization constant adds none in 15 of 15 gauge samples.
- **Same law as the parent, tested both ways** (the census tested one way):
  - 17 of 18 L-phantom children have the parent's law.
  - The exception is ter-bi (1,1255), which is also W_admit. Its parent reproduces the child at only 5 of 12 random child points, even with 51 starts.
  - All 14 admitted L-phantom children in the sample have the parent's law.

**N**
- **Exact leading-exponent test, all 12,061 laws:**
  - Substrates grow with products at zero.
  - Products grow with substrates at zero.
  - Each axis, pair, the diagonal and 200 random directions, plus the substrate axes with products held positive.
  - No direction is unbounded.
- **Numeric scans:** 80 of 15,864 scans (300 candidates per reaction) kept rising up to t = 1e30. All 80 have an exact exponent gap ≤ 0 and level off or decay by t = 1e240.
- **Oracle scans:** all 36 candidates level off.
- **The 334 degenerate seeds:** today's law grows without bound in none.

**Z**
- **Unflagged:** 18 children whose key differs from their parent's, fitted both ways:
  - 17 have a strictly larger law than the parent.
  - PC (1,1955), a plain split child with an RE-isomerization gauge, has exactly its parent (0,2251)'s law. The key keeps it apart because of one extra phantom constant; P catches it.

## What the census got wrong or missed
1. **Trap false alarms, 16 of 16.** "Findings only on single-reactant faces: 0" is wrong in substance: 16 of the 17 face-W mechanisms fail only on a single-reactant face. A face check built on this trap test would reject mechanisms whose full-face law is right.
2. **"The package count exceeds the slice count in 250 fitted children."** 250 is the count within P. Over all fitted children it is 913: 109, 82, 387 and 335 per expansion round.
3. **"One more identifiable dimension" for the 902 Xlaw.**
   - It is one more in 775, two in 120, three in 5 and four in 2.
   - Where only the package law is deficient, the oracle agrees with the slice. Ter-bi (1,639) `E + B ⇌ E(B) [OnlyA]; E + Q <--> E(Q) [EqualAI]; (E(B) + A <--> E(A, B), E(Q) + A <--> E(A, Q)) [EqualAI]; E(Q) + P <--> E(A, B) + D [OnlyA]` has slice rank 7/7, oracle rank 7, package rank 8/9.
4. **L-phantom share.** It is 13.8–24.7%, not 14–27%; 27% is all of P in ter-bi round 2.
   - There are 1,048 L-phantom records, not 1,047.
   - There are 63 RE-isomerization gauges (57 without an L phantom), not 60.
5. **"Gauge and other: 0 of 23 have one."** 22 of those 23 were never compared with any candidate. PC (1,1955) is a counterexample.
6. **Cheapest-first arithmetic.** 841 of the 2,046 D are also W_admit. After the island rule, the key removes 1,205 of the expansion fits, not 22–40%:

   | | ter-bi R1 | ter-bi R2 | ter-ter R1 | PC R1 |
   |---|---|---|---|---|
   | D not also W_admit | 211 (16.3%) | 142 (24.8%) | 517 (21.9%) | 335 (16.8%) |

7. **The class union counts candidates deficient only in today's package law.** With the rewrite's classes alone, the share of expansion fits drops:

   | | ter-bi R1 | ter-bi R2 | ter-ter R1 | PC R1 |
   |---|---|---|---|---|
   | Census union | 37.3% | 49.7% | 47.2% | 56.5% |
   | Rewrite's classes only | 36.6% | 47.7% | 45.3% | 55.0% |

8. **"14 flux-free" is right in substance.** Each has a :NonequalAI steady-state group and no chemistry in the inactive conformation, even though the census's own log labels them `ss_local`.
9. **No duplicates hidden by sampling (ter-bi only).** No child repeats the key of an unsampled base mechanism; the full base tier has 4,952 distinct keys.

## Not checked
- Fits and the cost of each check; pyruvate dehydrogenase; corner faces; n > 1 or regulatory sites (none occur in this census).
- Keys of the unsampled ter-ter and PC seeds.
- Whether the slice's parameters cover every thermodynamically allowed constant set. The oracle uses the slice's parameters.
- The face scan over all 12,061 laws used one parameter point and one speed ordering per face, with two more orderings only where the law disagreed.
- At t = 60 the steepest vanishing patterns ran out of precision, so those verdicts rest on t = 30 and 45.
- A "different" verdict from fitting could in principle be an optimizer failure.

## Files
Everything is under `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/large/verify`, with one log per script in `logs/`:
- **`scripts/`:**
  - Libraries: `vlib.jl`, `vfit.jl`
  - All-candidate and W checks: `v0_positive.jl`, `v1_w.jl`, `v1b_w.jl`, `v1c_wall.jl`
  - D and Z: `v2_dz.jl`, `v2b_dz.jl`, `v2c_fp_pc.jl`
  - P: `v3_p.jl`, `v3b_samelaw.jl`, `v3c_case.jl`, `v8_attr.jl`
  - N: `v4_n.jl`, `v4b_n.jl`
  - End to end and recounts: `v5_e2e.jl`, `v5b_compare.jl`, `v6_counts.jl`, `v7_fullkey.jl`
  - Queries: `q1.jl`–`q3.jl`
- **Ignore part (a) of `logs/v2_dz.log`.** My bug gave L a random value in one member of each pair; `logs/v2b_dz.log` is the corrected run.
- **`data/`:** `v0_positive.jls`, `v2_splits_*.jls`, `v5_*.jls`.

Julia time was about 45 minutes, every run through the lock. The package checkout (`git status` is clean) and `slice/src` are untouched.