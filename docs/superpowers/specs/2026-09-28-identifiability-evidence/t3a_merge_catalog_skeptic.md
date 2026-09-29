# Skeptic verdicts: Points 3-8 and 10: merge and Theorell-Chance catalog

Independent recheck of the track report `t3a_merge_catalog_report.md` with separate code.

## Summary

Denis, I checked 10 claims by recomputing them with my own code. I wrote an exact-rational steady-state evaluator with an exact RE limit, an exact rank via forward-mode derivatives, my own chain, merge and inverse maps, and a sympy rederivation of Table L, and cross-checked against engine B. Seven claims reproduce exactly. L1 is right only under a condition it leaves out. D2 is right, but a 'limit only' remark in the report is wrong for 20 ping-pong cases. X2's counts are exact, but its example (and P5b) cites the wrong mechanism.

**Confirmed exactly:**
- **Chain lemma, standard orientation:** holds on all 227 pairs, in both directions.
- **X1:** the export-wide phantom removal numbers all match. R4 removes 159 of 616, R5 297 of 771, R6 178 of 975. Mechanisms made identifiable: 127, 259 and 159. No mismatches.
- **X2:** duplicate counts match: R4 108, R5 188, R6 80.
- **RE merge (D4a):** exact on 207/207 random bi-bi cases, and none of the 207 parents is identifiable.
- **Package ping-pong seeds (D4b):** R4_00056–60 are 6/6 and their merges 5/5, so the merge loses a dimension and is not equal.
- **Ping-pong (P6):** six-step 63/63 non-identifiable. The gauge direction checks out: the law fixes only 1/kf3 + 1/kf6. The textbook family splits into RSSS / SSSR / RSSR, and I showed this with exact parameter maps, which is stronger than the track's fits.
- **Random bi-bi (P7c):** all counts match: 256/81, 49 invalid merges, rank drops 176/30/1, monomials kept in 56, 39 of 80.
- **Theorell–Chance (P8) and counts (A3):** parameter counts and the lost monomials match on every ordered, uni-uni and ping-pong case.

**Problems found:**
1. **Rule 1 is too broad.** The chain lemma only works when the two flank steps bind ligands into the chain and release them out of it (or carry none). The rule's text never says this, and the track's chain finder doesn't check it. Counterexample: iso uni-uni, E+S→ES, ES→F+P, F→E, all SS. The rule would merge F with E. The parent is 5/4 and has a P·S term; the "representative" is 3/3 and loses it, so the iso family disappears. The export has no such chains, so X1 and X2 are unaffected. RE merges don't need the condition: they stay exact wherever the ligands sit.
2. **Wrong mechanism ID in P5b and §6.** R4_00467 belongs to R4_00003's family (the seed with dead ends at EAQ). The seven-parameter RSSSR child of R4_00004 is R4_00481.
3. **"Limit only" is wrong in 20 ping-pong cases.** There the plain merge keeps rank 6 and gives a proper piece of the parent's family at finite parameters. This doesn't change any rule.

**Implementation caveats:**
- **Merged forms can't go back into expansion as the package stands.** Merged forms are exactly the folded-chemistry forms `_assert_chemistry_is_iso` rejects before expansion (src/mechanism_enumeration.jl:163-172, called at 2487). For now, canonicalize only for dedup and parameter counting.
- **Canonicalizing before expansion can make families unreachable.** A dead end bound to only one of the two chain forms can't be reached from the merged form. The current moves never produce this case (0 found in R3, R5, R6).

Scripts and data are in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/verify_t3a_merge_catalog/ (vk.py, chainmaps.py, table_l.py, verify_chain_own.py and .json, chain_export_own.py and chain_export_own_R*.json, dup_own.py, re_merge_own.py, pp_seeds.py, pp_verify.py, random_ss.py and .jsonl, ordered_catalog.py, misc_checks.py, plain_merge_d2.py, iso_chain.py, onesided.py, iso_counterexample.json, deadend_chain.py).

## Verdicts

### L1: overstated

Reproduced for the setting the proof uses. I rederived Table L symbolically with sympy (table_l.py): all 7 unmerged rows and merged SS match. I wrote my own chain maps (chainmaps.py) and checked all 227 parent/representative pairs of uni-uni, ordered and ping-pong in both directions (verify_chain_own.py). Every pair satisfied Haldane exactly after mapping, gave identical laws at exact rational points, had the predicted fitted drop and had equal exact rank: 0 failures. BUT the claim and Rule 1 leave out a condition the proof silently assumes: s1 must only bind ligands into X1 and s2 must only release ligands out of X2 (or the flank must be ligand-free). Counterexample (iso_chain.py, own rank plus engine B). Iso uni-uni all SS: E+S->ES, ES->F+P, F->E. F-E is ligand-free, F and E each have exactly one other step, and every step owns its group, so it passes the stated rule as an SSS chain. The parent is 5/4 with den {1,P,S,P*S}. The class-representative merge is 3/3 with den {1,P,S}: P*S is lost and the rank drops, so it is NOT EQUAL. A one-sided case also fails: E+S<->ES RE, ES<->F+P RE, F->F2 SS, F2<->E RE is 4/3, but its class representative is 4/4, a different family. With the orientation condition added, the lemma holds everywhere I tested, including all qualifying export chains (X1).

### A3: confirmed

Own fitted counts (vk.py). MERGE: -1 on all 207 random RE-isomerization merges and -2 on all 207 valid random SS-isomerization merges. The chain maps confirmed every predicted drop. ELIM on the 15 merged ordered types (ordered_catalog.py): -2 exactly when both flanks are SS (RSSR 5->3, RSSS 6->4, SSSR 6->4, SSSS 7->5), otherwise -1, and the fused step is RE only for RR flanks (RRRS, SRRR, SRRS). Uni-uni SS,SS gives 3->1. Ping-pong X1 and X2 ELIM give 7->5. Every invalid MERGE or ELIM (49 random, 1 ordered, 3 ping-pong double merges) was an all-RE catalytic cycle, never a Keq-forcing or stoichiometry failure, consistent with A1 and A2.

### D4a: confirmed

re_merge_own.py uses my own rescaling map (K*p_i, kr/p_i, p1=1+K0, p2=(1+K0)/K0) and an inverse with a random K0. On all 207 valid random bi-bi RE-isomerization assignments: Haldane holds exactly on the mapped points, laws are identical at exact rational points in both directions, fitted drops by 1, the exact rank is unchanged, and 0/207 parents are identifiable. Unlike the SS chain rule, this needs no condition on where the flank ligands sit: the iso uni-uni free-enzyme RE isomerization F<->E (types SSR, RSR, SRR) also merges exactly (4->3, 3->2, 3->2).

### D4b: confirmed

Own counts and exact rank (pp_seeds.py): R4_00056, 57, 58, 59 and 60 are 6/6 as parents and 5/5 after merging EB_res=EQ, so the rank drops and the merge is not EQUAL. R4_00061 and 62 are 6/5 -> 5/5. Wording caveat: 'groups straddling X1 and X2' fits only some seeds. In R4_00059 the shared groups link steps at X1 (EB_res+P, E_res+B->EB_res) to steps elsewhere. In R4_00060 they link steps at X2 (EQ+A->EAQ, E+Q) to E+A->EA and EA+Q. Rule 2's actual test ('every step incident to X1 or X2 owns its group') is the right one and excludes all five. The NEITHER relation for R4_00057 was not checked.

### D2: confirmed

plain_merge_d2.py checks all uni-uni, ordered and ping-pong chain cases. The rank is kept for every RE isomerization (RRR 17, RRS 21, SRR 21, SRS 21) and for SSS (21). RSR: 4 invalid, 17 lose rank. RSS/SSR: 22 lose rank. Each rank-keeping case is shown EQUAL by exact maps (D4a, L1). Nuance: in 20 ping-pong cases (one chain RSS or SSR, the other chain FULL), the plain merge KEEPS rank 6 because the parent carries the ping-pong gauge phantom. These are still not EQUAL, by P6c's sign invariant: (C,FULL) laws have D_AB*D_P - D_AP*D_B > 0. But the result there is an attained proper piece (R strictly inside P at finite parameters), not the 'limit only' relation that report §7(ii) asserts for one-RE-flank merges.

### X1: confirmed

chain_export_own.py uses my own chain detection, class-representative merge, fitted counts and exact rank, plus exact forward-map law identity with Haldane on every qualifying mechanism. Phantoms removed: R1 4/4, R2 6/6, R3 8/16, R4 159/616, R5 297/771, R6 178/975. Mechanisms made identifiable: 3/3, 5/5, 3/10, 127/418, 259/567, 159/906. Mechanisms with a qualifying chain: 4, 7, 11, 262, 607, 660. 0 fitted, law or Haldane mismatches. My single-draw rank disagreed on 6 R5/R6 mechanisms; with 5 random draws all 6 agree, so those were unlucky points.

### X2: confirmed

dup_own.py groups the export by canonical form after iterative class-representative merging. It counts mechanisms whose group holds a cheaper member: R1 3, R2 5, R3 7, R4 108, R5 188, R6 80. No group has tied minima. The counts match exactly. The example ID is wrong, though (the same error is in P5b and §6). R4_00467 is the RSSSR child of R4_00003 (the seed with dead ends at EAQ, parent R4_00076), and its family equals R4_00003's. The RSSSR child of R4_00004 is R4_00481 (7/5, parent R4_00081), which does fall in R4_00004's group with R4_00081 and R4_00082 (6/5).

### P6c: confirmed

My check is stronger than the track's family fits: exact constructive maps (pp_verify.py). For 20 random textbook SSSS laws I computed the sign b1*c2 - a2*d1 and moved along the gauge to d1=0 or c2=0. I then built explicit RSSS (12 laws) or SSSR (7 laws) parameters: Haldane holds exactly and the laws are identical at exact points, 0 failures. Conversely, 10 random RSSS laws are exact textbook laws (moving epsilon < 0) and all have sign > 0. So the pieces are disjoint and each lies inside the textbook family. Supporting results: the six-step ping-pong has 63 valid assignments, none identifiable. The merged menu has 15 valid, and only SSSS (7/6) is non-identifiable. The effective D formula matches the law exactly. The gauge move leaves v bit-identical and keeps 1/kf3+1/kf6, 1/kr1+1/kr4, kr3/kf3 and kr6/kf6 while kf3 changes (P6a and P6b confirmed).

### P7c: confirmed

random_ss.py, with my own exact rank and engine B rank in agreement, and monomials from engine B. SS isomerization: 256 valid assignments, 81 identifiable. The plain merge is invalid in 49 (all-RE cycle) and valid in 207, all with -2 fitted. Delta rank: 2 in 176, 1 in 30, 0 in 1 (RRSSSSSRR). Monomials are unchanged in 56 (Delta rank 2/1/0: 39/16/1). All 80 identifiable parents with a valid merge give identifiable merges, and 39 of those keep the monomials. Fully SS is 15/15 with 18 numerator and 48 denominator monomials; the merge is 13/13 with the same 18 and 48. P7e is also confirmed: RRSSSSSRR equals the combined-port effective law exactly (5 parameter points x 4 concentration points), and the three counts 7/7, 11/7 and 9/7 are reproduced.

### P8: confirmed

ordered_catalog.py, engine B monomials. On all 15 merged ordered types, ELIM loses exactly the claimed sets and gains none. SSSS, SSRS, SRSS: {ABP, BPQ}. RRSR: {AB}. RSRR: {PQ}. RSSR: {AB, PQ}. SSSR: {AB, BPQ}. RSSS: {ABP, PQ}. TC SSS is 5/5. Uni-uni loses {S,P} (1/1, v=kf(S-P/Keq)). Ping-pong: X1 loses {AP}, X2 loses {BQ}, and fusing both halves gives 3/3 with den {A,B,P,Q} and no constant term. Counts are -2 exactly for SS/SS flanks. I did not prove the general support formula independently; only these cases were reproduced.

## New issues

- Rule 1 (canonicalize every chain) needs an orientation condition that its text and claim L1 leave out: each flank step's ligands must sit on the outer form, meaning s1 only binds ligands into X1 and s2 only releases ligands out of X2, or the flank is ligand-free. Without it the rule merges free-enzyme isomerizations of iso mechanisms. Iso uni-uni E+S->ES, ES->F+P, F->E (all SS, own groups) is 5/4 with den {1,P,S,P*S}. The rule's 'representative' E+S->ES->E+P is 3/3 with den {1,P,S}. The P*S family is lost, and the chain lemma would wrongly predict 2 phantoms where there is 1. The track's chain_export.py find_chains has no such check. The export has no inner-ligand chains (all export law identities held), so the X1/X2 numbers stand, but the rule as a general enumerator rule is unsafe for iso or conformational mechanisms. Scripts: verify_t3a_merge_catalog/iso_chain.py, onesided.py; mechanisms in iso_counterexample.json.
- The RE merge (Rule 2 / D4a) needs no orientation condition: RE merges of F<->E in iso uni-uni (SSR, RSR, SRR) are exact both ways. So an RE isomerization with exclusive incident groups is a phantom regardless of where the ligands sit, but an SS isomerization chain is canonicalizable only under the outer-ligand condition.
- Implementation conflict: for the package's ordered and ping-pong chains, the class representative is exactly the 'folded chemistry' form that the expansion path rejects. _assert_chemistry_is_iso is at src/mechanism_enumeration.jl:163-172 and is called before expansion at line 2487, and CLAUDE.md enumeration-test rule 3 says the moves misread folded forms. An ordered merged X=(EAB=EPQ) also has no single bound set for the current Species model. So Rules 1 and 2 cannot be dropped in as 'canonicalize after each move' without Species/Step and move changes. The workable interim form is to expand from the unmerged form and canonicalize only for dedup and parameter counting.
- Reachability caveat for canonicalizing before expansion: a dead end bound to only one chain form is not reachable from the merged X. Uni-uni RSR + I bound to ES only is 4/4 with den {1,P,S,I*S}; the merged SS,SS + I bound to X is 4/4 with den {1,P,S,I*P,I*S}, a different family. Under the current dead-end move this never occurs: 0 such children in R3, R5 and R6, because competitive inhibitors bind only free-site forms, and dead ends at ping-pong EA or EP_res share K across contexts. It would become a loss if a move ever targets a single central complex.
- The claim P5b and the §6 example cite the wrong ID. R4_00467 (7/5) is the RSSSR grandchild of R4_00003, the seed with dead ends at EAQ, via R4_00076, and its family equals R4_00003's. The RSSSR child of the plain ordered seed R4_00004 is R4_00481 (7/5, via R4_00081).
- Report §7(ii) says a one-RE-flank plain merge is 'limit only'. That is wrong in 20 ping-pong cases (one chain RSS or SSR, the other chain FULL): the plain merge keeps rank 6 and gives an attained proper piece of the parent, the same pieces as RSSS/SSSR in P6c. This does not change Rule 3, but it means the plain merge is itself a direct way to generate Rule 6's identifiable pieces.
- Partial data on the open question about shared-flank chains, from engine-A results in cv/results/export_a.jsonl. R4_00001 is 5/5. Flipping one shared flank group to SS keeps it identifiable (R4_00063, 64, 65, 66 are all 6/6). Flipping both flank groups creates phantoms (R4_00441 is 7/6; split children R4_00442 and R4_00445 are 8/7). So shared-flank chains do produce chain-like phantoms that the lemma does not cover.
- D4b wording: the non-exact RE merges in R4_00059 and R4_00060 are not caused by groups straddling X1 and X2. In each, a group joins a step at one of the two forms to a step elsewhere: R4_00060's EQ+A->EAQ shares its group with E+A->EA. Rule 2's 'incident steps own their groups' test is correct; only the explanation is off.
