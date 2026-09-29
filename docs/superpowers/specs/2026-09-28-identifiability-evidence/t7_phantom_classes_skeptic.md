# Skeptic verdicts: Phantom classes of flip and split children

Independent recheck of the track report `t7_phantom_classes_report.md` with separate code.

## Summary

Denis, I checked 10 claims with my own code. I rebuilt the flip and split moves from the Julia source, wrote my own segment-level zero-flux test, used engine B for ranks and laws, and verified parameter maps exactly with sympy.

**Confirmed:**
- **C3, C4, C5 (zero-flux class):** 188 mechanisms with 336 directions in R4 and the same in R5. Reverting the zero-flux groups gives exactly the same rate law (376/376 exact checks per reaction). Every one of the 148 flip cases reverts to its own parent.
- **C12:** the Zariski-closure argument holds.
- **C13:** both positivity gaps (R4_00847 and R4_00490) reproduce independently.
- **C18:** the reachability losses reproduce exactly: 64 identifiable mechanisms lost at level 3 and 272 at level 4 (280 under the reject-phantom rule).

**Overstated:**
- **C19 (the lossless zero-flux completion).** The track's numbers reproduce exactly, but only for its own implementation. That implementation runs a separate search for minimal supersets, capped at 2 extra groups. The recommended rule says to make a zero-flux flip set count as failing inside the package's existing minimal-set search. Implemented that way, it gives a graph identical to "reject zero-flux children", which loses 64 identifiable mechanisms at level 3 and 272 at level 4. There is also a gap at the level of families: 72 dropped split children at level 4 have reverted forms that the moves never reach. Emitting that reverted form (same law, fewer parameters) instead of dropping the child closes the gap.
- **C9 (merging the central complex).** Merging keeps the family only when the steps that enter the block have private groups. When a group is shared, the rank drops: R4_00441 goes from fitted 7, rank 6 to fitted 5, rank 5, and this happens in 12 of 22 such cases in the export.
- **C10 (deleting a parallel channel).** The same privacy condition applies.
- **C20 (repair table).** Its numbers match, but the statement that only reverting repairs AB is false. Fusing the SS chemistry step with the neighbouring RE step into one SS step (for example EAB → EQ + P) keeps the rank with one fewer parameter in 119 of 147 AB mechanisms. It gave the same family in both directions for 8 of 8 pairs I tested.

**Rule 5 (ping-pong seeds).** Fixing seeds R4_00061/62 removes the phantom from only 35 of the 60 affected mechanisms; the other 25 come from other seeds. Dropping the two seeds would lose their two families, since nothing else reachable matches them.

Everything is in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/verify_t7_phantom_classes/. The main scripts are vlib.py (moves and zero-flux test), graph.py and compare.py (reachability), zlaw.py, sbcheck.py, pccheck.py, richeck.py, posgap.py, abcheck.py and abelim.py.

## Verdicts

### C3: confirmed

Proof reviewed: Hill cycle-flux on Cha's segment-level system holds, and so does the lemma that an edge lies on a catalytic cycle iff its biconnected block's cycle space has sigma != 0 (cycles through one edge span a 2-connected block). Independent check with my own block-based zero-flux test (vlib.zero_flux_groups) and engine B. Reverting the Z groups keeps the rank and lowers fitted by #Z in 188/188 R4 and 188/188 R5 mechanisms. The law is unchanged under K=kf/kr at random consistent points (engine-B laws, 80 digits): 376/376 in R4 and 376/376 in R5, worst relative difference 8.9e-80. Negative control (K=1.01 kf/kr on R4_00470) gives 1.1e-4. R6 sample: 36 Z mechanisms, revert keeps the rank (engine A) and lowers fitted by #Z in 36/36. Scripts: zcheck.py, zlaw.py, negctl.py, zR6.py.

### C4: confirmed

EXACT with my own classifier plus engine-B ranks. R4: 336 Z directions in 188 mechanisms. Flips give 296 directions in 148 mechanisms; splits give 40 directions in 40 mechanisms, and every split's zero-flux step has form degrees (1,3). 174 mechanisms are pure Z and 14 are mixed. R5 counts are identical. zcheck.py -> zrecs_R4.json, zrecs_R5.json.

### C5: confirmed

EXACT over all 148 Z flip children. Reverting the zero-flux groups gives exactly the recorded parent (canonical key) in 148/148. The parent's SS chemistry step has both ends in one RE segment in 148/148. My re-implementation of _flux_carrying_groups (parent form-graph biconnected blocks) flags every flipped group as flux-carrying in 148/148 (zcheck.py).

### C9: overstated

The lemma holds when the block's two entry steps have private groups. The two-way maps from the lemma, checked on engine-B laws at 80 digits (3 points each way), are exact on all 12 private degree-2 blocks: R1_00004, R2_00005 and 10 in R4 (sbcheck.py). At depth <=4, 165/165 such instances keep the rank. 'Entered from outside by exactly two steps' is not sufficient. Counterexample R4_00441: EA+B->EAB shares its group with EQ+B->EBQ, and EQ+P->EPQ shares with EA+P->EAP. MERGE of EAB->EPQ goes from fitted/rank 7/6 to 5/5, so it is not EQUAL. This happens in 12 of 22 such blocks in the R4 export and in 96 of 261 degree-2 instances at depth <=4 (repairs.py, deeprep.py).

### C10: overstated

Confirmed for private channel groups. The explicit map kr_e'=kr_e+kr_f, kf_e'=kf_e(kr_e+kr_f)/kr_e gives law identity (engine B, 80 digits) for every PC-pattern mechanism in the R4 export: 74 mechanisms, 148 deletions, 296/296 points. Rank is kept and fitted drops by 1 in all (pccheck.py). At depth <=4, 2016/2016 private deletions keep the rank with -1 fitted. The general statement needs privacy: with shared channel groups (all inside Z-bearing nodes), 24 depth-4 deletions raise both fitted and rank by 1, and 40 change neither (deeprep.py, deeprep2.py).

### C12: confirmed

Proof reviewed. A flip parent is the Lambda->infinity limit, and a split parent is a specialization, so the parent family lies in the closure of the child family. The closures of positive-parameter images, which are connected real-analytic sets, are irreducible and have dimension equal to the generic rank. Equal dimension therefore forces equal closures, and the converse is immediate. It says nothing about positive families; C13 shows that gap.

### C13: confirmed

Independent check with engine-B laws (posgap.py). For parents R4_00203 (d_AP*d_B - d_P*d_AB) and R4_00085 (d_PQ*d_A - d_Q*d_AP), the invariant is a ratio of polynomials whose coefficients are all positive, so it is positive everywhere. The children take negative values: R4_00847 at 991/2000 random rational points, R4_00490 at 756/2000. Explicit negative points are given in posgap.py output.

### C18: confirmed

EXACT with my own flip/split re-implementation, written from the Julia source, not moves.py. It reproduces export levels 1-2 (369/1388, 0 missing, 0 extra) and the track's graph (L3 5446 nodes; L4 12448 nodes, 29146 edges). Reject Z children: L3 4426 nodes, 64 identifiable lost (all fitted 9); L4 272 lost (240 at 10, 32 at 11). No rank gain: L3 3998 nodes, the identical 64/272 sets. Reject phantom: L3 3697 nodes, 64 lost at L3 and 280 at L4. Nuance: 8 of the 280 phantom-rule losses stay reachable when only Z nodes are removed, so 'every loss runs through a zero-flux intermediate' holds for the Z and rank rules only (graph.py, compare.py).

### C19: overstated

My independent implementation of the track's variant reproduces its numbers exactly: L3 6420 nodes (5567 identifiable, 853 phantom), 1808 new identifiable (1238/506/64 at 10/11/12 fitted), 0 Z nodes. L4 has 12634 identifiable, 2090 phantom and 0 lost. With the rank filter, L3 has 5726 nodes, 159 phantom and 0 lost. Three caveats. (a) The recommended wording, 'count the set as failing so the minimal-set search extends it', does not produce this variant. Plugged into the package's Apriori _minimal_gaining_sets, it yields a graph identical node for node to 'reject zero-flux children' (4426 nodes at L3), which loses 64 identifiable at L3 and 272 at L4. The track's code (moves.expand_flip_zf_complete) instead runs a separate minimal-superset search per failing set. (b) The counts depend on its max_extra=2 cap: uncapped gives 6476 nodes and 1864 new identifiable at L3, and 15012 nodes at L4. (c) 'Loses nothing' covers identifiable nodes only. 72 dropped Z split children at depth 4 have Z-reverts that appear nowhere in the whole move closure (50,039 mechanisms, exhausted at depth 12), and no variant node has an exactly equal parameter space (samespace.py, cover2.py, deep.py).

### C20: overstated

Engine-B recomputation reproduces most of the table. Z revert keeps the rank 188/188. Deleting a channel keeps the rank in all 74 PC-pattern mechanisms. MERGE of the RE isomerization keeps the rank in 60 mechanisms; 56 of them are exact reparametrizations with private groups (richeck.py, 112/112 points), and all 97 identifiable instances lose rank. MERGE keeps the rank on the track's SB ids. For AB (all 147 ids), some single-group revert keeps the rank in 147/147, MERGE of the SS isomerization in 0/147, and ELIM of a two-step all-SS transit in 0. But under the shared-context ELIM definition (fused step is SS if either step is SS), ELIM of the absorber keeps the rank with -1 fitted in 119/147 AB mechanisms. Example: EAB->EPQ (SS) plus EQ+P<->EPQ (RE) fused into EAB->EQ+P (SS). Exact coefficient solves at random points show EQUAL in both directions for 8/8 tested pairs, including R4_00490 (abcheck.py, abelim.py); R4_00081 is EQUAL by a hand map. So 'AB: only revert' and 'MERGE and ELIM always shrink the family' are false for ELIM in general.

## New issues

- Rule 1 as worded is lossy. Making a zero-flux flip set 'fail' inside the package's Apriori _minimal_gaining_sets gives a graph identical, node for node, to rejecting zero-flux children: 4426 nodes at L3, with 64 identifiable mechanisms lost at L3 and 272 at L4. Apriori only extends a failing set when all its subsets also failed, and a sibling cut has usually gained already. The lossless variant needs a separate minimal-superset search for each failing minimal set. The track capped that search at 2 extra units; uncapped gives 6476/15012 nodes at L3/L4 and also loses 0.
- Family coverage under zero-flux completion. At depth 4 the variant drops 72 Z split children. Their Z-reverts are unreachable anywhere in the full R4 move closure (50,039 mechanisms), and no variant node has an exactly equal parameter space. For 7 sampled 8/8 reverts, an 8-parameter same-support node contains the law at a random point. For the eight 9/8 reverts, only 4 have a containing node at the test point. Alternative: emit the Z-revert (EQUAL, fewer parameters) instead of dropping the child, together with completion. Every dropped node then has an exact twin at depth <=4, and the variant gains 112 nodes: 104 identifiable, including 64 identifiable 8-parameter mechanisms that the baseline offers only as 9-fitted Z nodes (graph.py variant zrev).
- The SB rule's predicate ('ligand-free SS step whose two forms are both top transits') is too loose. MERGE loses 1 rank in 12 R4 export mechanisms whose entry groups are shared, e.g. R4_00441 (7/6 -> 5/5). At depth <=4, on nodes free of zero-flux groups, 72 shared degree-2 instances and 312 higher-degree instances lose 1-2 ranks. Only 'each transit has exactly the iso step plus one entry step, and the entry groups are private' is safe (165/165 at depth <=4). The predicate also misses the 4 ping-pong SB-class mechanisms (R4_01811/14/17/18), where MERGE of EA->EP_res (only one side is a top transit) keeps the rank at -2 fitted.
- Rule 5 (ping-pong seeds 61/62) covers only 35 of the 60 RI-phantom mechanisms (those where MERGE of the RE isomerization keeps the rank). The other 25 descend only from seeds R4_00056-60, e.g. R4_00408 (flip of 59) and R4_01685 (split of 56). Dropping seeds 61/62 would lose their rank-5 families: within depth 4, no node reachable from any other seed has the same monomial support.
- AB class: ELIM of the absorber complex keeps the rank with -1 fitted in 119/147 AB mechanisms, and is EQUAL in 8/8 tested pairs. It fuses the SS chemistry step with the neighbouring RE binding/release step, e.g. EAB -> EQ + P. This is exactly the Theorell-Chance-style fused step of Denis's points 8/11. It needs a Step whose end forms differ by more than the released metabolite. The track's 'no single operation other than revert' holds only for ELIM of all-SS transits.
- family_fit calibration. A provably EQUAL pair (identical per-step parameter spaces after Wegscheider ties; bijective map verified by sympy simplify == 0) gave family_fit residuals of 1.1e-4 and 2.2e-2 (5 laws, 8 starts). Single-law residuals up to ~0.2 in the track's SAMPLED tables (C16) therefore cannot count as evidence of non-equality. Exact coefficient solving with sympy, allowing parametric solutions and positivity sampling, was reliable here (positive control passes).
- The rank-gain rule (optional rule 3) drops 572 non-Z phantom nodes at L3 (1257 at L4) that have no exact twin, including R4_00847, which has a proved positivity gap. '0 identifiable lost' is true, but family losslessness is not established for this rule.
