# Skeptic verdicts: Points 10 and 12: which unmerged isomerizations to keep

Independent recheck of the track report `t3b_point12_census_report.md` with separate code.

## Summary

Denis, I rechecked ten of the track's claims with my own code: my own MERGE, a fitted-count routine that matches the package on all 6,534 mechanisms, a rebuilt R4 move graph that reproduces the export exactly, engine B ranks, supports and laws, and sympy. **All the counts reproduce.** That covers MERGE (Δfitted −2 on 1,464 R4 and 2,936 R5 merges, validity, supports, novelty), the series rule (22 mechanisms, 33 phantoms), P1 (328 removed, 616 → 92 phantoms, no full-rank mechanism removed), P3 (5 unreachable full-rank mechanisms), the two ping-pong families lost under R1/R2, and the cost of a modular rank (1.3 ms median, 0 disagreements with engine B).

**The main correction is one the track missed.** Two mechanisms with equal rank, one inside the other's closure, do not have the same family over positive parameters.
- **P1 loses real models.** R4_00413 and R4_00419 are rank-neutral splits of two ping-pong seeds, and P1 removes them. Their positive families are strictly larger than their seeds' (PROVED from exact laws). The split lifts an inequality that holds in the seed: in effect it allows cooperative binding on the modified enzyme. No kept rank-6 mechanism covers those laws. In a sample of 25 mechanisms P1 removes, 8 had random laws their kept equal-rank ancestors could not fit. So P1 and P3 lose more than the track says, and 'non-identifiable, hence useless' is not safe as a premise.
- **The series rule still holds at depth ≤ 2, but for a different reason.** The 10 clean cases are exactly equal to their seed's law. I proved a stronger lemma: they equal the mechanism with both flanking steps in rapid equilibrium, which here is always a seed. In the 12 obstructed cases, neither level-1 parent covers the mechanism on its own. The union of both parents and the seed does (by hand for R4_00441; sampled for four more).
- **Beyond depth 2, the rule needs a guarantee.** It is lossless only if the enumerator also emits those parents or the merge. The safest version replaces such a mechanism with its merge.

The claim "P1 loses 8 mechanisms" is right in effect. Four of the eight do have kept equal-rank relatives through split edges that only appear when you expand to depth 3, but those relatives fail to fit their laws.

I did not recheck R5's P1/P3 counts or the 12 R5 merges of shared-group isomerizations. Everything is in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/verify_t3b_point12_census/:
- vlib.py, moves.py, graph.py: MERGE, fitted count, and the rebuilt move graph (graph_R4_d3.json)
- census.py, run_b.py, analyze_merges.py: merge census and engine B runs (b_R4_*.jsonl, b_R5_*.jsonl)
- p1.py, p1b.py: P1 and P3 simulation
- series.py: series-rule check
- chain_check.py: sympy check of the chain lemmas
- law413.py, ineq.py: the separating inequalities
- ffu.py, p1_cover.py: union fits (p1_cover.jsonl, log_ffu_series.txt)
- time_rank.py: rank timing

## Verdicts

### C1: confirmed

Reproduced with my own MERGE and my own cycle-constraint fitted count, which matches pkg_fitted on all 1,819 R4 and 4,715 R5 mechanisms, plus engine B ranks. R4: 1,464/1,464 valid merges have Δfitted = −2 (EXACT). Δrank is 0 in 14, −1 in 165, −2 in 1,253, −3 in 24 and −4 in 8; none is positive. R5: 2,936/2,936 single-step merges have Δfitted −2 and Δrank ≤ 0. I did not recheck the track's 12 R5 merges of shared-group isomerizations.

### C2: confirmed

My own RE union-find and cycle-uptake check. R4 has 1,866 SS isomerizations, 1,464 valid and 402 invalid: level 0 62/62, level 1 114/376, level 2 226/1,428. R5 has 4,811 SS isomerizations: 1,857 invalid, 2,936 valid, and 18 in shared groups. The 'iff' is sound. Contracting s maps an all-RE cycle through the merged form onto an RE path X1→X2 in P, which together with s is a catalytic cycle whose only SS step is s. The converse holds the same way. No valid merge creates a bottomless RE segment (0/1,464 R4, 0/2,936 R5).

### C6: confirmed

chain_check.py: sympy re-derives the three-step chain (A,B,C,D). Its construction reproduces 300 random rational two-step quadruples exactly with positive rates. The result is stronger than claimed. An RE-flanked chain Y⇌X1 (RE), X1→X2 (SS), X2⇌Z (RE) gives (A,B,C,D) = (s·Ka, s'·Kb, Ka, Kb), which is also onto R+^4. So a clean both-series C equals C_ab exactly, where C_ab is C with both flank groups turned RE. For all 10 clean R4 cases C_ab is a seed (series.py). The 10 therefore equal their seed laws exactly over positive parameters. The track's C8 had only 3/3 sampled plus closure equality.

### C10: confirmed

Conclusion reproduced. R4_01650 is 8/8 and merges to 6/6; R4_01788 is 8/7 and merges to 6/5 (engine B). The track's argument ('no kept relative has equal rank') does not by itself rule out non-relatives or higher-rank containers, so I checked further. R4_01650: no R4 mechanism of rank ≥ 8 has superset supports, even shifted. R4_01788: the rank-7 superset-support candidates cannot be closure-equal because their generic supports are larger. The only rank-8 candidate that R1 keeps, R4_01720, fits R4_01788's laws at 0.21–0.56 on 3 laws. The other 22 allSS parents are not lost (see C11).

### C11: overstated

The counts reproduce exactly: 22 isomerizations in 22 mechanisms, 33 phantoms, 10 clean. For all 22, C_ab is a seed. Clean cases are exact (see C6). The obstructed 12 are where the track's proof fails. Its argument, 'closure-equal to a level-1 parent, so nothing lost', is invalid over positive parameters. R4_00441 by hand: the phantom is u + ρw = 1 + ρ. Parent R4_00065 only reaches u < 1, parent R4_00066 only u > 1, the seed only u = 1, and C covers all of (0, 1+ρ). Union fits on 5 obstructed cases (00441, 00456, 00524, 00754, 00814) agree. Each law is fitted by one level-1 parent (~1e-39); where I have per-parent residuals, the other parent fails (3e-4 to 1.5e-1). The union of both parents and the seed covered every law. So at depth ≤ 2 in R4 nothing is lost, because both single-flank parents and C_ab exist. That is PROVED for the 10 clean cases and for R4_00441, and SAMPLED for the other 4 tested. As a general rule it is lossless only if C_ab (and, when obstructed, both single-flank parents) are also emitted.

### C12: confirmed

My merges plus engine B for R4. 1,464 merges, all distinct, none structurally equal to a population mechanism. 1,400 have a (supports, rank) signature new to the population. 162 phantom mechanisms carry 172 phantoms. 220 merges of level-1 parents have 4 fitted parameters, 196 of them identifiable, and the minimum seed count is 5. R5 single-step merges: 2,936 distinct, 2,864 new signatures, 242 mechanisms with 292 phantoms. That fits the track's 2,942 / 2,870 once group merges are included. Supports (C3) also reproduce: allSS 24/24 unchanged; someRE 100 unchanged, 1,272 lose denominator terms only, 56 lose numerator and denominator terms, 12 change in the numerator.

### C15: overstated

The counts reproduce on my own rebuilt R4 move graph, which I validated against the export: raw 255/114 and 1,522/624, new sets equal the export, 3,627 new at level 3. With engine B ranks: 328 removed (24 L1, 304 L2), 0 full rank, 90 phantom mechanisms left, phantoms 616 → 92. Also 0 rank drops on all 2,663 in-population edges (C13). But the track shows only that the 320 removed families are closure-equal to a kept ancestor, not that the ancestor contains them over positive parameters, and here that fails. PROVED: R4_00413 and R4_00419 (7/6, rank-neutral splits of seeds R4_00059 and R4_00060, removed by P1) have strictly larger positive families. From engine B's coprime laws, normalized by the B coefficient: seed 00059 has a_BQ·a_PQ − a_BPQ·a_Q = K2²K4/(K5K6) > 0 always, while in 00413 it equals K2²K8(K5K7 − K6 + K7)/(K5²K7²), negative in 37% of random points. Likewise seed 00060 has a_AB·a_BQ − a_ABQ = K1K2/K4 > 0, while 00419 gives K2K6(K2K5 + K2 − K3K5)/(K3K5). No kept rank-6 mechanism shares their supports. SAMPLED: in 8 of 25 stratified P1-removed mechanisms, at least 1 of 4 random laws was not reached by the union of kept equal-rank ancestors (residuals 5e-3 to 2.7e-1). So P1 loses more than 8 families. I did not recheck R5.

### C16: confirmed

R4 part reproduced. R4_01761–64 and R4_01782–85 (9/8) are reachable only through the P1-rejected splits R4_00413 and R4_00419. They stay unreachable at depth 3. One nuance: my depth-3 graph has lateral split edges from kept 8/8 mechanisms into 4 of them (R4_01756→01762, 01759→01763, 01770→01774→… i.e. 01770→01783 and 01774→01784). By the track's own closure criterion those 4 would count as kept, so its loss criterion is inconsistent. Over positive parameters they are still lost: R4_01756 failed on 4 of 6 laws of R4_01762 (1e-4 to 7e-2). R5 not rechecked.

### C17: confirmed

My graph gives the same results. P3a removes 423, and its 5 unreachable full-rank mechanisms are R4_00423, 01789, 01798, 01803 and 01804, as stated. For P3b I count 104 phantoms with no kept equal-rank relative using depth-3 edges; the track has 108 with depth-2 edges. 'P3 loses families' is now PROVED rather than sampled: R4_00413 and R4_00419 are phantom, P3 drops them, and their positive families are not contained in any kept rank-6 mechanism (see C15).

### C19: confirmed

idkit.rank on all 1,819 R4 mechanisms: median 1.33 ms, p95 2.0 ms, max 3.1 ms, total 2.4 s. 0 disagreements with the engine B ranks. Package compile + fitted_params cost is 1,098 s / 1,819 ≈ 0.60 s per mechanism (export README). The ratio is about 0.2%.

## New issues

- Positive families differ from their closures, and every closure-based 'no loss' argument in this track misses it. Equal rank plus P ⊂ cl(C) (Lemmas 3 and 6) does not give fam(C) ⊆ fam(P). A rank-neutral ('pure phantom') child can lift an inequality that holds on the parent's positive family, for example allowing cooperativity between B and P binding on F in ping-pong. PROVED for R4_00413 ⊋ R4_00059 and R4_00419 ⊋ R4_00060 (ineq.py); SAMPLED for 8 of 25 P1-removed mechanisms. So P1 (and P3) remove genuine models, not only redundant ones, and Denis's premise 'non-identifiable, hence useless' needs this caveat.
- Series rule, obstructed half: a single parent does not cover C. The union of both single-flank parents (C with one flank RE) and C_ab covers it. So the rule is lossless only if all three are emitted. They are at depth ≤ 2 in R4, but nothing guarantees it at deeper levels.
- Series rule, clean half: fam(C) = fam(C_ab) exactly, where C_ab is C with both flank groups RE. This lemma is new and PROVED: an RE-flanked chain realizes all of R+^4. The rule is lossless only if C_ab or the merge M is emitted. At depth ≤ 2 in R4, C_ab is always a seed (22/22). Beyond that, the minimal-gaining-set flip logic may never emit C_ab. A safer rule replaces C by its merge M, which has the same family with 2 fewer parameters.
- Depth-3 expansion adds lateral split edges between level-2 mechanisms (e.g. R4_01756 → R4_01762). The track's depth-≤2 edge set therefore undercounts relatives. Minor count shifts follow: P3b has 104, not 108, phantoms with no equal-rank kept relative.
- Not rechecked: the 12 R5 merges of shared-group isomerizations (all obstructed), and all R5 P1/P3 counts, which would need the dead-end move reimplemented.
