# Track t3b — Which unmerged central isomerizations to keep (Denis's points 10 and 12), measured on the R4/R5 populations

All numbers below come from the package's own enumerator output (export at `ident/export/`, package HEAD 96cdf23) and from engine A. Nothing in the repository was modified. Julia was not run. Scripts and data: `ident/t3b_point12/`.

## 0. Findings in one page

1. **MERGE always costs exactly 2 fitted parameters** (EXACT: 4,412 valid merges, R4 1,464 + R5 2,948). It is valid (no all-RE catalytic cycle) for 1,464 of the 1,866 SS isomerizations in R4. The 402 invalid ones include all 62 seeds: a seed's isomerization is the only SS step on a catalytic cycle.
2. **The concentration terms stay the same only when both ends of the isomerization are steady-state-only.** All 24 "allSS" merges keep both supports. Only 100 of 1,440 "someRE" merges do; 1,340 lose terms. Point 10's "often the same conc terms" holds for allSS neighbourhoods and fails for the rest.
3. **Exact theory.** When the merge is *clean* (no kinetic group at X1/X2 also covers steps elsewhere), MERGE is exactly the parent with the isomerization in rapid equilibrium, relabelled (Lemma 1: PROVED; explicit parameter maps checked on all 740 clean merges, worst mismatch 1.1e-42). The result therefore lies in the closure of the parent's family (Lemma 2). If, in addition, both ends are *series* forms (exactly one other step, SS), the families are **exactly equal over positive parameters** (Theorem 4: PROVED; exact sympy check). The rank then stays the same, and 2 phantoms disappear.
4. **Point 12, answered on the data.** At depth ≤ 2 only 24/1,866 SS isomerizations (R4) have allSS ends: 22 have both ends series, and 10 of those are clean. The 10 clean parents are exactly their merged forms, and the merged forms are exactly the 5-parameter seed laws they descend from (3/3 checked both ways < 1e-35). These parents are pointless: the seed's law plus 2 phantoms. The 12 obstructed ones are closure-equal to their own level-1 parents, so they are redundant too. Their merged forms are 5/5 laws with the seed's supports that the seed does not fit (12/12, residuals ≥ 2e-2). **All 1,842 other SS isomerizations should be kept unmerged.** Their merge changes the law (Δrank −1 to −5).
5. **Rules.** R1 (drop allSS parents) and R2 (replace allSS parents by merge) each lose 2 families: the ping-pong mechanisms R4_01650 (8/8) and R4_01788 (8/7), whose ends are not series. R2c (replace only clean allSS) loses nothing, but its 10 additions duplicate seed laws. **Series rule** (drop a parent whose SS isomerization has two series ends): removes 22 mechanisms and 33 phantoms, loses 0 families. R3 (add every merge as a child) adds 1,464 distinct mechanisms (+80%) to R4, 1,400 of them with new (supports, rank) signatures, 172 phantoms, 0 losses. 220 of them have only 4 fitted parameters, one fewer than any seed.
6. **Rank policies.** P1 (reject a flip/split child whose rank does not exceed its parent's) removes 328/1,819 (R4) and 422/4,715 (R5), all phantom, and 85% of R4's phantoms (616 → 92). It removes no full-rank mechanism. 320/328 removed families are provably kept by a kept equal-rank ancestor. 8 R4 mechanisms (rank 8) are reachable only through a rank-neutral split; they stay unreachable at depth 3, and no kept relative fits them. P3 (fit only full-rank mechanisms) loses families (108–115 R4 candidates; 10/10 tested not contained). P2 re-tiers 418 R4 mechanisms and removes nothing.
7. **Cost.** Engine A's modular rank takes 1.6 ms (median; p95 4.7 ms, max 7.4 ms) per R4/R5 mechanism in pure Python: 3.6 s for all of R4 and 7.7 s for all of R5. The package spends 1,098 s and 3,643 s compiling the same populations. A rank test costs about 0.3% of the current per-mechanism cost, and P1 would skip compiling 18% (R4) / 9% (R5). **A rank filter would speed enumeration up, not slow it down**, but only with a modular rank. The package's BigFloat SVD rank (≈0.6 s per mechanism) would not do.
8. **Tool caveat.** family_fit (8-start local LM) misses proven containments often. On 42 clean pairs whose merged law provably lies in the closure of the parent's family, it returned < 1e-20 for 19, 1e-20..1e-3 for 13, and ≥ 1e-3 for 10. A family_fit residual ≥ 1e-3 is therefore not evidence of non-containment. The "incomparable" outcomes for obstructed merges are open.

## 1. Setup and definitions

**Populations.** R4 = bi-bi `A[CX] + B[N] ⇌ P[C] + Q[NX]` (ping-pong-capable), levels 0/1/2 = 62/369/1,388 mechanisms. R5 = R4 + competitive inhibitor I, 62/719/3,934. Levels 1–2 are the children of `_expand_re_to_ss`, `_expand_split_kinetic_group` and `_expand_add_dead_end_regulator`. Structurally, R5 contains every R4 mechanism. fitted/rank values are engine A's (exact rank over GF(2^61−1); they agree with the package on every mechanism cross-checked).

**Complete move graph** (needed for P1/P3). I rebuilt every parent→child edge. Flip and split come from the Python re-implementation in `t7_classes/moves.py`, with one fix I added: the package excludes groups binding a regulator from the flip units, and the original omitted that exclusion (254 extra R5 children before the fix). Dead-end edges are the export's provenance, which is complete because the dead-end move produced 0 duplicates. **Validation (EXACT):** per level and per move, the raw child counts equal the package's (R4 0→1: 255/114; 1→2: 1,522/624; R5 0→1: 255/114/350; 1→2: 2,961/1,731/2,153). The set of new children equals the exported next level. Every exported first parent is among the produced parents. Edges: R4 2,515, R5 7,564.

**MERGE(X1,X2)** of a ligand-free SS step s: X1→X2. The two forms become one form "X1=X2". Every other step keeps its metabolites, type and group. s and its group are deleted. *Groups holding other steps:* only in R5, where 18 SS chemistry steps (9 groups) share a group with their I-bound mirror copy (ping-pong with I on the modified enzyme). I merge every step of that group; the two I-binding RE steps then become identical (same ends, same metabolite, same group) and are collapsed into one. All of these are obstructed. **ELIM(X)** for a merged form with exactly two steps Y→X (L1, R1) and X→Z (L2, R2): delete X and both steps, and add Y→Z consuming L1+L2 and releasing R1+R2 in a new group, typed SS if either input was SS, else RE. **P_RE** = the parent with s turned RE.

**Neighbourhood of s:** *allSS* when every other step at X1 and X2 is SS; otherwise *someRE*. **Series end:** an end X_i whose only steps are s and exactly one other step, which is SS. **Clean:** every kinetic group (other than s's) holding a step at X1 or X2 holds only steps at that same end, all with that end on the same side (all "from" or all "to"). No other step joins X1 and X2. Otherwise the merge is **obstructed**. The typical obstruction is a binding group shared with an abortive dead-end binding elsewhere (e.g. B binding EA→EAB and EQ→EBQ in one group).

**Identifiable** = rank = fitted. **Phantoms** = fitted − rank. **fam(M)** = the set of rate laws v(concentrations) over positive parameters (Keq and E_total fixed). **cl** = Zariski closure. **family_fit(src,tgt)**: worst-over-5-laws best residual of fitting tgt to src's random laws; ~0 means fam(src) ⊂ fam(tgt). "fit P to M's laws" = family_fit(M, P). Classes: < 1e-20, 1e-20..1e-3, ≥ 1e-3, timeout (600 s per direction), memory (LM drove parameters to ~e^(10^k) and mpmath exhausted a 1.8 GB cap; that is a divergence towards a limit).

## 2. Proofs

**Lemma 1 (clean MERGE = RE limit, relabelled).** If the merge of s is clean, fam(M) = fam(P_RE).
*Proof.* Let K be s's constant in P_RE ([X2] = K[X1]), c1 = 1/(1+K), c2 = K/(1+K). Map P_RE parameters to M parameters group by group. A group whose steps have end X_i as "to": RE K_g' = K_g/c_i; SS kr_g' = c_i·kr_g. A group with X_i as "from": RE K_g' = c_i·K_g; SS kf_g' = c_i·kf_g. Other groups are unchanged. Clean ⇔ the map is well defined per group. Take P_RE's steady state and set [X] = [X1]+[X2], all other forms unchanged. Every mass-action term of every step of M then equals the corresponding term of P_RE, because c_i[X] = [X_i]. Every RE relation holds. The segment balance equations are the same equations, since s is internal to a segment of P_RE. Total enzyme is unchanged. By uniqueness of the steady state this is M's steady state, and v (net A uptake, which s does not touch) is equal. Thermodynamics: a cycle of M through X that enters at X_i and leaves at X_j picks up log(c_j/c_i). That is 0 for i = j and ±log K for i ≠ j, which equals s's contribution to the corresponding P_RE cycle. So the constraints map onto each other. Conversely, given M parameters, choose any K > 0 and invert the map. ∎
*Check (EXACT):* on all 740 clean valid merges of R4 ∪ R5, both directions, with 3 random laws × 6 points in 50-digit arithmetic: the mapped parameters satisfy the target's constraints and the laws agree. Worst relative mismatch 1.1e-42 (`explicit_map.py`).

**Lemma 2 (RE limit is a limit).** If s joins two different RE segments of P, then fam(P_RE) ⊂ cl(fam(P)).
*Proof.* Put kf_s = λκK, kr_s = λκ. Each King–Altman spanning tree of the segment graph uses s at most once, so N and D are affine in λ. The λ-coefficients are the trees through s, i.e. the trees of the graph with s's ends contracted, which gives the Cha law of P_RE. So v_P → v_{P_RE} as λ → ∞. Every valid merge has s between different segments: were X1 and X2 in one segment of P, the RE path between them would carry the overall reaction, and P_RE, and with it M, would have an all-RE catalytic cycle. Remark: the same argument shows MERGE is invalid iff s is the only SS step on some catalytic cycle of P. ∎

**Lemma 3 (equal dimension).** Both closures are irreducible, with dimension equal to the generic Jacobian rank. If fam(M) ⊂ cl(fam(P)) and rank(M) = rank(P), then cl(fam(M)) = cl(fam(P)). ∎

**Theorem 4 (both-series clean merges are exact).** If s is clean and both ends are series forms, then fam(M) = fam(P) over positive parameters.
*Proof.* Such a P contains a chain Y +L1 ⇌ X1 ⇌ X2 ⇌ Z +R2 (three SS steps, own groups). The rest of the mechanism sees it only through the flux J = A·L1·[Y] − B·R2·[Z] and the chain's enzyme content [X1]+[X2] = C·L1·[Y] + D·R2·[Z]. The cycle constraints see only A/B. For M's two-step chain, (k1,k−1,k2,k−2) ↦ (A,B,C,D) = (k1k2, k−1k−2, k1, k−2)/(k−1+k2) is a bijection onto R+^4, with inverse k2 = A/C, k−1 = B/D, k1 = C(k−1+k2), k−2 = D(k−1+k2). Hence P ⊂ M. For M ⊂ P, the three-step chain (a1,b1,s,s',a2,b2) gives A = a1a2s/Σ, B = b1b2s'/Σ, C = a1(a2+s+s')/Σ, D = b2(b1+s+s')/Σ, with Σ = a2b1 + a2s + b1s'. Given a target, take s = s' = σ > max(k2,k−1), a2 = 2σk2/(σ−k2) and b1 = 2σk−1/(σ−k−1), so that A/C = k2 and B/D = k−1. Then a1 and b2 are fixed by C and D, all positive. ∎ *Check (EXACT):* sympy derives both chains' (A,B,C,D). The construction reproduces 200 random rational two-step chains exactly with positive rates, and the inverse map reproduces the three-step quadruple symbolically (`gadget.py`). The uni-uni all-SS chain (R1_00004, 5/3) merges to 3/3 with the same supports. That is Denis's point 4, and it is a special case.

**Lemma 5 (ELIM is a limit when unshared).** If X's two steps are in their own groups and X is the bound form on both sides (true for every merged central complex here: EA+B→X, EQ+P→X, or E+A→X, E_res+P→X), then fam(ELIM) ⊂ cl(fam(M)). Take kr_in = λα and kf_out = λβ with λ → ∞. Then [X] → 0 and J → (kf_in·β·L1[Y] − α·kr_out·R2[Z])/(α+β), which is a fused step with arbitrary positive (kf, kr). An RE step at X is handled the same way (K_in → 0 and k_out → ∞ with the product fixed). ∎

**Lemma 6 (moves are monotone).** For a flip child C of P, P = RE limit of C (Lemma 2), so P ⊂ cl(C). For a split child, P = C with the two parts equal, so P ⊂ C. For a dead-end child, P is the limit with no inhibitor binding, so P ⊂ cl(C). Hence rank(C) ≥ rank(P). If rank(C) = rank(P), then cl(C) = cl(P) by Lemma 3. *Check (EXACT):* 0 rank drops on all 2,515 R4 and 7,564 R5 edges.

## 3. Part A — MERGE and ELIM on every SS isomerization of R4 and R5

### A1. Validity (R4; R5 in brackets)
SS isomerizations: 1,866 [4,811]. Valid MERGE: 1,464 [2,948]. Invalid (all-RE catalytic cycle): 402 [1,863]. Invalid ⇔ s is the only SS step on some catalytic cycle (PROVED, remark to Lemma 2). By level (R4): level 0: 62/62 invalid (every seed); level 1: 114/376; level 2: 226/1,428. Group merges (R5 only): 18 steps in 9 two-step groups, 12 valid, 6 invalid.

| rx | nbhd | parent | class | SS isos | valid | invalid |
|---|---|---|---|---|---|---|
| R4 | allSS | ident | ping-pong | 1 | 1 | 0 |
| R4 | allSS | phantom | ping-pong | 3 | 3 | 0 |
| R4 | allSS | phantom | sequential | 20 | 20 | 0 |
| R4 | someRE | ident | ping-pong | 160 | 142 | 18 |
| R4 | someRE | ident | sequential | 1,272 | 1,046 | 226 |
| R4 | someRE | phantom | ping-pong | 108 | 90 | 18 |
| R4 | someRE | phantom | sequential | 302 | 162 | 140 |
| R5 | allSS | (same 24) | | 24 | 24 | 0 |
| R5 | someRE | ident | ping-pong | 484 | 328 | 156 |
| R5 | someRE | ident | sequential | 3,741 | 2,226 | 1,515 |
| R5 | someRE | phantom | ping-pong | 196 | 144 | 52 |
| R5 | someRE | phantom | sequential | 366 | 226 | 140 |

### A2. Δfitted, Δrank, identifiability (valid merges)
| rx | nbhd | parent | class | n | Δfitted | Δrank | P ident | M ident | P phantoms | M phantoms |
|---|---|---|---|---|---|---|---|---|---|---|
| R4 | allSS | ident | ping-pong | 1 | −2:1 | −2:1 | 1 | 1 | 0 | 0 |
| R4 | allSS | phantom | ping-pong | 3 | −2:3 | −2:1 0:2 | 0 | 1 | 6 | 2 |
| R4 | allSS | phantom | sequential | 20 | −2:20 | −1:12 0:8 | 0 | 20 | 28 | 0 |
| R4 | someRE | ident | ping-pong | 142 | −2:142 | −2:142 | 142 | 142 | 0 | 0 |
| R4 | someRE | ident | sequential | 1,046 | −2:1,046 | −4:8 −3:24 −2:1,014 | 1,046 | 1,014 | 0 | 40 |
| R4 | someRE | phantom | ping-pong | 90 | −2:90 | −2:41 −1:45 0:4 | 0 | 36 | 109 | 56 |
| R4 | someRE | phantom | sequential | 162 | −2:162 | −2:54 −1:108 | 0 | 88 | 182 | 74 |
| R5 | allSS | (the same 24 as R4) | | 24 | −2:24 | −2:2 −1:12 0:10 | 1 | 22 | 34 | 2 |
| R5 | someRE | ident | ping-pong | 328 | −2 | −2:328 | 328 | 328 | 0 | 0 |
| R5 | someRE | ident | sequential | 2,226 | −2 | −5:8 −4:32 −3:40 −2:2,146 | 2,226 | 2,146 | 0 | 128 |
| R5 | someRE | phantom | ping-pong | 144 | −2 | −2:67 −1:73 0:4 | 0 | 58 | 169 | 88 |
| R5 | someRE | phantom | sequential | 226 | −2 | −2:54 −1:172 | 0 | 152 | 246 | 74 |

Δfitted = −2 in every valid merge (EXACT). A merge never gains rank (Δrank ≤ 0, EXACT). It can create phantoms: 32 identifiable R4 parents give phantom merges (Δrank −3/−4).

### A3. Supports (R4 [R5])
| nbhd | num, den | n |
|---|---|---|
| allSS | same, same | 24 [24] |
| someRE | same, same | 100 [100] |
| someRE | same, den loses terms | 1,272 [2,680] |
| someRE | num and den lose terms | 56 [116] |
| someRE | num gains, den loses and gains | 12 [28] |

### A4. Structure predicts the rank change (R4; R5 identical in pattern)
| sharing | series ends | Δrank: n |
|---|---|---|
| clean | both | 0: 10 |
| clean | one | 0: 4 (ping-pong), −1: 114 |
| clean | none | −1: 23, −2: 131 |
| obstructed | both | −1: 12 |
| obstructed | one | −1: 8, −2: 276 |
| obstructed | none | −1: 8, −2: 846, −3: 24, −4: 8 |

Δrank = 0 occurs in 14 of the 1,464 valid R4 merges (EXACT): all 10 clean both-series merges (8 sequential, 2 ping-pong) and 4 clean one-series ping-pong merges. No obstructed merge and no other sequential merge keeps the rank. As a rule of thumb, each clean series end removes one phantom.

### A5. Family relation on a stratified sample
88 pairs: all 24 allSS, plus 8 per someRE cell of parent-ident × sequential/ping-pong × clean/obstructed. Cells were drawn from R4 ∪ R5 with seed 20260928.

| nbhd, class, sharing | n | Δrank | proved relation | fit P to M's laws | fit M to P's laws |
|---|---|---|---|---|---|
| allSS seq clean | 8 | 0 | **EQUAL** (Thm 4) | <1e-20: 8 | <1e-20: 2, 1e-20..1e-3: 6; recheck (24 starts, new laws): <1e-36: 8/8 |
| allSS pp clean | 2 | 0 | **EQUAL** (Thm 4) | <1e-20: 2 | <1e-20: 1, timeout 1; recheck <1e-35: 2/2 |
| allSS seq obstructed | 12 | −1 | not EQUAL (rank) | <1e-20: 4 (RESULT ⊂ PARENT, SAMPLED), 1e-20..1e-3: 4, ≥1e-3: 4 | ≥1e-3: 12 |
| allSS pp obstructed | 2 | −2 | not EQUAL (rank) | ≥1e-3: 2 | ≥1e-3: 2 |
| someRE pp clean | 16 | 0:1, −1:6, −2:9 | M ⊂ cl(P) (L1+L2), strict unless Δrank 0 | <1e-20: 4, small: 8, ≥1e-3: 4 | ≥1e-3: 14, timeout 2 |
| someRE seq clean | 16 | −1:8, −2:8 | M ⊊ cl(P) | <1e-20: 5, small: 5, ≥1e-3: 6 | ≥1e-3: 15, memory 1 |
| someRE pp obstructed | 16 | −1:5, −2:11 | not EQUAL (rank) | <1e-20: 1, small: 5, ≥1e-3: 10 | ≥1e-3: 16 |
| someRE seq obstructed | 16 | −1:4, −2:10, −3:2 | not EQUAL (rank) | <1e-20: 1, small: 4, ≥1e-3: 11 | ≥1e-3: 15, timeout 1 |

For the whole populations (not only the sample), the relation follows exactly from Lemmas 1–4 and the exact ranks. R4: 10 EQUAL (clean both-series); 4 closure-EQUAL (clean, Δrank 0, ping-pong one-series); 268 RESULT ⊊ cl(PARENT) (clean, Δrank < 0); 1,182 obstructed, all with Δrank < 0, hence not EQUAL. Whether an obstructed M lies in cl(P) is open.

**Positive-control calibration** of family_fit: on the 42 clean sample pairs (M ⊂ cl(P) proved), fit P to M's laws gave < 1e-20: 19, 1e-20..1e-3: 13, ≥ 1e-3: 10. On the 10 exact-EQUAL pairs, fit M to P's laws needed 24 starts to reach < 1e-35 in 8 cases. So the "≥ 1e-3 in both directions" rows for obstructed pairs are not evidence of incomparability.

**Uni-uni (Denis's point 4, R1):** R1_00004, E+S→ES→EP→E+P all SS (5 fitted / rank 3), merges to E+S→X→E+P, 3/3, same supports. EQUAL: Theorem 4, and fits < 1e-39 on 5/5 laws (P ⊂ M) and < 1e-37 on 4/5 (M ⊂ P; the fifth 4e-12 is LM stalling on the redundant target). The level-1 uni-uni parents (4/3) merge to 2/2 with a lost S or P term (clean one-series, Δrank −1).

### A6. Examples (step lists; → is the stored direction; RE/SS; g = group)
* **Clean both-series, R4_00481 (level 2, 7/5):** E+A⇌EA RE g1; E+Q⇌EQ RE g2; EA+B→EAB SS g3; EAB→EPQ SS g4; EQ+P→EPQ SS g5. MERGE → EA+B→X SS; EQ+P→X SS (5/5). Supports num {AB, PQ}, den {1, A, AB, PQ, Q} for both. EQUAL to P (Thm 4) and EQUAL to seed R4_00004 (E⇌EA⇌EAB RE, EAB→EPQ SS, EPQ⇌EQ⇌E RE, 5/5): fits 2e-38 / 5e-39. The same holds for R4_00467 ↔ R4_00003 (1e-38 / 2e-39) and R4_00500 ↔ R4_00006 (6e-39 / 2e-36). Closure equality with the seed also follows from Lemma 6 along the equal-rank chain seed (5/5) → level 1 (6/5) → P (7/5).
* **Obstructed both-series, R4_00441 (7/6):** as above, but B also binds EQ→EBQ in g3 and P also binds EA→EAP in g4 (abortive complexes sharing the groups). M is 5/5 with the same supports as seed R4_00001, but the seed does not fit M's laws (0.1) and M does not fit the seed's (LM diverged). P itself is closure-equal to its level-1 parents R4_00065/66 (6/6) by Lemma 6.
* **allSS but not series, R4_01650 (ping-pong, 8/8, identifiable):** E+A→EA SS and EP_res+A→EAP_res SS (g1, substrate inhibition by A on F·P); EA+P→EAP and E_res+P→EP_res SS (g3); EA→EP_res SS; the rest RE. M is 6/6 with identical supports: a more constrained law (point 10's case). The parent's 8-dimensional family is lost if the parent is replaced.
* **someRE, level 1, R4_00079 (6/6):** E+A→EA SS; everything else RE (E+Q, EA+B, EQ+P); EAB→EPQ SS. M (E+A→EA the only SS step, X in rapid equilibrium) is 4/4:
  v = K3·kf1·(Keq·A·B − P·Q) / (K3·Keq·B + K2·K3·Keq·B·Q + K2·K5·Keq·P·Q + K2·K3·K5·Keq·B·P·Q).
  This is the rate-limiting-A-binding law (Cha). The parent's denominator has {1, A, AB, B, BPQ, BQ, PQ, Q}; the merge keeps only {B, BQ, PQ, BPQ} (divided by B: 1 + Q + PQ/B + PQ), so there is no A term. 220 distinct merges of R4 level-1 parents are 4-parameter (196 identifiable); no seed has fewer than 5.

### A7. ELIM of merged forms with exactly two steps
R4: 366 of 1,464 merged forms have degree 2 (R5: 734); all ELIMs are valid, 258 identifiable (R5 534).

| nbhd | fused | steps at X in shared groups | Δfitted vs M | Δrank vs M | n (R4) |
|---|---|---|---|---|---|
| allSS | SS | 0 | −2 | −2 | 10 (9 identifiable Theorell–Chance) |
| allSS | SS | 1 | 0 | −1 | 8 |
| allSS | SS | 2 | +2 | 0 | 4 |
| someRE | RE | 0 | −1 | −1 (84), −2 (8) | 92 |
| someRE | RE | 1 | 0 | 0 | 40 |
| someRE | RE | 2 | +1 | +1 | 12 |
| someRE | SS | 0 | −1 | −1 | 104 |
| someRE | SS | 1 | 0 / +1 | 0 | 72 |
| someRE | SS | 2 | +2 | +1 | 24 |

The denominator always loses terms (the central-complex terms). Example: R4_00481's M → EA+B→EQ+P SS, with E⇌EA and E⇌EQ RE: 3/3, den {1, A, Q}. With unshared steps, ELIM ⊊ cl(M) (Lemma 5 + Δrank < 0). With shared groups the shared groups keep their other members, so ELIM *adds* parameters (0 to +2). Fits on 6 such pairs: ≥ 1e-3 both ways (open, see calibration). *Other type choice:* making an RE-fused step SS instead gives +1 fitted and +1 rank (136 cases) or +2 rank (8), with different supports, in all 144 R4 cases: a different, larger model.

## 4. Part B — rules for point 12 (R4; R5 identical in every count except R3)

| rule | removed | added (distinct) | net | phantom mechs / phantoms removed | phantoms added | families lost | new families |
|---|---|---|---|---|---|---|---|
| R1 keep SS iso unmerged only if an end has an RE step (drop allSS) | 24 | 0 | −24 | 23 / 34 | 0 | 2: R4_01650 (rank 8), R4_01788 (rank 7) | 0 |
| R2 merge-and-replace allSS | 24 | 24 | 0 | 23 / 34 | 2 | same 2 | 14: 12 obstructed merges (not fitted by the seed, SAMPLED) + 2 new signatures |
| R2c replace clean allSS only | 10 | 10 | 0 | 10 / 21 | 1 | 0 (Thm 4) | 0: all 10 equal seed laws |
| **series rule**: drop a parent whose SS iso has both ends series | 22 | 0 (or 12 obstructed merges) | −22 | 22 / 33 | 0 | 0 (clean: Thm 4; obstructed: closure-equal to kept level-1 parent, Lemma 6) | 0 (12 if merges added) |
| R3 add every valid merge as a child | 0 | 1,464 [R5 2,942] | +1,464 (+80%) [+2,942, +62%] | 0 | 172 in 162 mechs [292 in 242] | 0 | 1,400 new (supports, rank) signatures [2,870] |
| R3s add merges for someRE only | 0 | 1,440 [2,918] | +1,440 | 0 | 170 | 0 | 1,398 [2,868] |

"Families lost" = removed parents with no kept equal-rank ancestor (Lemma 6) whose merge is not EQUAL. The 2 R1/R2 losses are exact: their merges drop the rank by 2, and no kept relative has their rank. Signature checks on R3: of 42 R4 merges whose (supports, rank) signature is already in the population, 8 were fitted both ways against the match. 3 are EQUAL: R4_00429's merge ≡ R4_00062 (the 6/5 phantom ping-pong seed, so the merge gives it an identifiable 5/5 form), R4_01794 ≡ R4_00421, R4_01814 ≡ R4_00428. 1 is one-way (K ⊂ M). 4 are unfitted.

## 5. Part C — rank-based policies and cost (depth ≤ 2, complete move graph)

| policy | R4 removed | of which full rank | phantom mechs (phantoms) left | families |
|---|---|---|---|---|
| none | 0 | – | 418 (616) | – |
| **P1** reject flip/split child with rank ≤ parent rank | 328 (24 at level 1, 304 at level 2) | 0 | 90 (92) | 304 closure-equal to a kept parent, 16 to a kept ancestor (PROVED, Lemma 6); **8 lost**: R4_01761..64 and R4_01782..85 (9/8). Their only route is a rank-neutral split of ping-pong seed R4_00059/60 (P binding at E_res vs EB_res), then a flip. They are unreachable at depth 3 too (graph extended to 3,627 new level-3 mechanisms). Fits against kept same-support rank-8 relatives: 2/2 not contained (2e-2 to 2e-1). |
| P1 on R5 | 422 | 0 | 145 (147) | 398 + 16 kept; 8 lost (the same pattern) |
| P2 charge rank, not fitted | 0 | – | 418 (616) | none removed; 418 R4 mechanisms (23%) drop tier: by 1: 225, by 2: 188, by 3: 5 (R5: 567) |
| P3a fit only full rank, no expansion of phantom mechs | 423 (incl. 2 phantom ping-pong seeds) | 5 (descendants of seed R4_00061) | 0 | 308 kept by equal-rank relatives; 115 candidates; 10 tested against kept same/superset-support relatives: 0/10 contained (fits 3e-2 to 4e-1) |
| P3b fit only full rank, expand everything | 418 | 0 | 0 | 310 kept by equal-rank relatives; 108 candidates (70 have a kept parent/sibling with the same supports and rank ≥; 38 have none) |

All P1 removals are phantom by construction: fitted(C) > fitted(P) ≥ rank(P) = rank(C). P1 removes the same set whether or not dead-end edges are also tested. After P1 the remaining R4 phantom mechanisms are 2 seeds, 9 at level 1 and 79 at level 2 (mostly 8/7). P3 variants discard mechanisms that are the only depth-≤2 carriers of some families. For example, R4_00062's 5-dimensional family has no full-rank representative in the population; the MERGE of R4_00429 would supply one.

**Cost (measured, this machine under load).** Engine A fast rank per mechanism: median 1.63 ms, p95 4.67 ms, max 6.0 ms (R4); 1.52 / 3.18 / 7.4 ms (R5). Total 3.6 s (R4) and 7.7 s (R5); 0 disagreements with the stored ranks. Package compile + fitted_params + show: 1,098 s (R4, 0.60 s per mechanism) and 3,643 s (R5, 0.77 s). Package BigFloat rank: 155 s for 250 mechanisms. So P1 costs ~0.3% of today's per-mechanism time, runs before compile_mechanism, and removes 18% (R4) / 9% (R5) of compiles. P2 and P3 cost the same rank call. A Julia port of the GF(p) rank (King–Altman/Cha linear solve per random point, reverse-mode dv/dθ) is needed. The rank never needs the symbolic law.

## 6. Answers
* **Point 10 (merge regardless of identifiability).** Correct that the merge is worth offering regardless of identifiability: it always removes exactly 2 parameters and never loses a family (the parent is kept). But it changes the concentration terms in 1,340 of 1,440 someRE cases. When clean, it is exactly the RE limit of the central isomerization, a new, smaller model (e.g. rate-limiting substrate binding). Identical terms occur only when both ends are steady-state-only (24/24) or in 100 someRE cases. On identifiable parents the merge is identifiable in 1,157/1,189 (R4). 32 merges of identifiable parents are phantom.
* **Point 12 (which unmerged to keep).** Keep every SS central isomerization unless both of its ends are series forms (one other step each, SS). Those parents add nothing: when clean they are exactly the merged law (Thm 4), which at depth ≤ 2 is the seed's law. When obstructed they are closure-equal to their level-1 parents. The allSS criterion (R1/R2) is slightly too broad: it drops two ping-pong families (ranks 8 and 7, R4_01650 and R4_01788) whose ends are not series. Denis's "fully SS ordered E=EA=EAB=EPQ=EQ=E is pointless" is right. Theorem 4 applies to every clean series pair around the central complex (each step in its own group), so the all-SS ordered bi-bi equals its one-complex merge exactly (9 → 7 fitted, as in the shared facts). The same holds for point 6: in the six-step all-SS ping-pong both chemistry steps have two series ends, so two applications of Theorem 4 give the textbook four-step form exactly (11 → 9 → 7 fitted, rank 6 throughout; the remaining phantom is the ping-pong one, not a merge issue). No rank filter is needed for this rule. On this population P1 also removes all 22 series-rule parents, since each has the rank of its level-1 parents.
* **Rank filter speed.** Measured: not a slowdown (see Cost). P1 captures far more (328 vs 22 mechanisms) and has 8 known losses at depth 2.

## 7. Open questions
1. For obstructed merges (1,182 in R4), is M ⊂ cl(P)? family_fit cannot decide (calibration above); an exact closure test is needed.
2. P1's 8 losses: allowing a rank-neutral split to be extended by a flip (a cross-move minimal-gaining set) would recover them; that design is not worked out.
3. The 220 four-parameter merges of level-1 parents sit below the seeds; should they be seeds (point 9)?
4. For clean one-series merges (Δrank −1), is there a variant (merge + flipping the RE side) that keeps the rank? Not tested.
5. Obstructed allSS merges have the seed's supports and rank but are not fitted by the seed in either direction. Are these genuinely different 5-dimensional laws, or a positivity boundary? Open.
6. All results are depth ≤ 2. The all-SS random-order bi-bi (point 7) is not in these populations. Its merge is obstructed when A/B binding groups are shared across contexts, so Theorem 4 does not cover it.

## 8. Files (all under ident/t3b_point12/)
- mlib.py: MERGE, ELIM and clean-flag operations.
- run_merge.py → merges.jsonl; tab_merge.py and features.py → merge_rows.json; make_tables.py → tables_A.txt.
- explicit_map.py → explicit_map.json; gadget.py.
- ff_merge.py → ff_merge.jsonl; ff_final.py → ff_final.txt; ff_batch2.py → ff_batch2.jsonl.
- novelty.py → novelty.json, novelty_pairs.json; rules.py → rules.json.
- edges.py → edges_R4.json, edges_R5.json; policies.py → policies.json.
- lost.py and lost_search.py → lost_classes.json, lost_candidates.json, lost_search.json.
- depth3.py → depth3.json; time_rank.py → time_rank.json; elim_alt.py → elim_alt.json.
- uniuni_merges.txt, uniuni_ff.json.
