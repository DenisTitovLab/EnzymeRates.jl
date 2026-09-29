# Skeptic verdicts: Point 9: new seeds

Independent recheck of the track report `t4_new_seeds_report.md` with separate code.

## Summary

Denis, I checked 10 of the track's claims with my own code and every one holds. I rebuilt all 1,323 seeds and MERGE/ELIM variants myself and computed exact degeneracy flags, fitted parameter counts and ranks with engine B; the track used engine A. I also added an exact Gröbner test that decides whether one mechanism's rate law is exactly a member of another's family.

**The 10 claims:** C4, C6, C7, C8, C10, C11, C12, C17, C18 and C22 are all confirmed.
- Two are stronger than the track showed. The 12 decorated ordered seeds are exactly NEITHER against their merged twin; the track had only sampled this with family_fit. The 8 ordered/random duplicate pairs are exactly EQUAL.
- The proposed 179-seed bi-bi set checks out. Every member is fitted = rank = 5 and non-degenerate. No EQUAL duplicates remain among the 23 pairs that share a rate-law support.

**Where the recommendations go too far:**
1. **The V/C/I check wrongly rejects some clean mechanisms.** It is exact only when no substrate or product also binds as a dead-end inhibitor. When one does, the V (Vmax) part rejects clean, identifiable mechanisms. In merged R6 variants it rejected 98 of 704 clean ones, and in merged R3 variants 2 of 4. The smallest example is E+S→X (SS), X⇌E+P (RE), E+S⇌ESinh (RE). Its rate law is (kf·S − kr·P/K)/(1 + S/Ki + P/K) with 3 of 3 parameters identifiable, yet V rejects it. This only happens when a merge is applied to a mechanism that already has such dead-end copies, as a merge-anywhere move would do.
2. **Rules 1 and 4 together cut off identifiable ping-pong families.** Rule 1 rejects degenerate mechanisms in every move; rule 4 drops today's ping-pong seeds. Today's moves never leave a seed's topology, so all 132 clean ping-pong mechanisms at levels 1–2 descend only from the 7 degenerate ping-pong seeds. With those seeds gone, none of the 132 is reached. 83 of them are identifiable, and none equals any merged variant without splits.
3. **The I check turns out to be unnecessary.** Suppose validity already requires every SS group to carry nonzero steady-state flux (your point 1). Then V+C alone gives the same verdict as V+C+I on everything I tested. The only 4 mechanisms that need I also have SS groups with zero flux.

Scripts and data are in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/verify_t4_new_seeds/; verification_summary.json holds the key numbers. No repository file was modified and no Julia was run.

## Verdicts

### C4: confirmed

I rebuilt all 248 single-SS-group MERGE variants with my own MERGE code and took exact flags from engine B's exact rational-point law (two parameter points, always the same result). All 248 are degenerate. In every one the flipped ligand is unbounded (b flag), and each also either has no Vmax or has another ligand that is not needed. 24 have phantoms (16 with one, 8 with two). One caveat: the proof in §5.2 has a gap. It says all forms lie in one RE segment, but that is false when the flipped group also holds an abortive step: in R4_00001 with only {B} SS, EBQ is its own segment. The EXACT result still holds.

### C6: confirmed

Engine B on my 372 two-group MERGE variants: 130 are clean, and every one has fitted = rank = 5. By class: 60 ordered/ordered (3 per seed), 56 ordered/random (2 per seed), 0 random/random, and 14 ping-pong (the {A,Q} and {B,P} sets). The 8 degenerate two-group variants with phantoms have rank 4. Files: my_results.jsonl, detail.py.

### C7: confirmed

I rederived the proof. The steady state gives X = (k1·L1·Y1 + k−2·L2·Y2)/(k−1 + k2) and J = (k1·k2·L1·Y1 − k−1·k−2·L2·Y2)/(k−1 + k2). These equal X1 + X2 and kf·X1 − kr·X2 of the RE–SS–RE box under the stated map, and the inverse map is positive. Exact check (lump_verify.py): the rational-function identity N_seed·D_var − N_var·D_seed ≡ 0 in A, B, P, Q holds at 3 random rational parameter points for all 8 seeds whose two X-adjacent groups are single-step. It fails for all 12 seeds where an adjacent group also holds an abortive step.

### C8: confirmed

I strengthened this with an exact test instead of family_fit (contain.py). It solves the coefficient-matching equations with a Gröbner basis over Q, on 2 random laws in each direction. The twins of seeds 00003, 00004, 00006, 00007, 00018, 00022, 00025 and 00026 each have a unique positive finite solution both ways: EQUAL. For the other 12 ordered seeds the basis is [1] in both directions, so there is no finite complex solution: NEITHER at finite parameters, exact per law rather than SAMPLED. Hand check for R4_00005: the variant's normalized coefficients satisfy (AP·Q)/(A·PQ) = 1 + Ka1·kf3/(Ka2·Keq·kf4) > 1, while the seed requires this ratio to equal 1.

### C10: confirmed

tc_check.py builds both 3-form Theorell–Chance cycles symbolically in sympy. Under the claimed map (b1, s1, b2, s2, u) = (t, r1, t_rev, r2, a1), the difference of the two laws is exactly 0. The inverse map (a2 = u·b1·s2/(s1·Keq·b2)) also gives 0. Gröbner (3 laws each way): 00004-TC ≡ 00026-TC and 00007-TC ≡ 00022-TC are EQUAL; 00004-TC vs 00007-TC has no solution.

### C11: confirmed

My ELIM variants, checked with engine B and exact flags. Every ELIM variant with one or two SS groups is degenerate. The all-SS triple (outer substrate, TC step, outer product) is clean for all 20 ordered/ordered seeds. It has fitted = rank = 5 for exactly 00003, 00004, 00006, 00007, 00018, 00022, 00025 and 00026; of the other 12, 8 have 6 and 4 have 7. For ping-pong, 4 clean half-TC variants are 5/5 (00057-X1, 00061-X2, 00062-X1, 00062-X2), and both-halves ELIM is always degenerate. Groups left holding only abortive steps fail my flux check whenever they are SS.

### C12: confirmed

Engine B plus exact flags. All 55 sequential seeds are clean with 5/5. All 7 ping-pong seeds carry the flags a:fwd:B, f:fwd:B, a:rev:Q and f:rev:Q. R4_00061 and R4_00062 are 6/5, the other five 6/6. The symbolic zero-product law of both R4_00062 and R4_00056 is Ka1·kf·A/(1 + Ka1·A), which also confirms C13 (H1).

### C17: confirmed

I wrote my own V/C/I implementation from the report text and compared it with my exact flags. 0 disagreements on: 1,163 valid seeds and MERGE/ELIM variants, the 1,819 R4 export mechanisms, 15 ordered ter-ter MERGE pairs, 1,460 R5 merged variants with I dead-ends, and 2,231 R6 export mechanisms (all of L0 and L1, plus 400 sampled from L2). The claim as scoped to its 3,056 mechanisms therefore reproduces. It does not generalize, though: V falsely rejects 100 clean, identifiable merged variants that carry substrate or product dead-end copies (see new_issues).

### C18: confirmed

reach.py on r4_results.jsonl (engine B plus exact flags): 93 of 1,819 are degenerate, all ping-pong. That is 7/7 seeds, 26/47 at level 1 and 60/171 at level 2; 47 of the 93 have phantoms. Side finding confirmed with an independent criterion (an SS group whose net steady-state flux is zero at a generic point): 188 mechanisms fail it, all at level 2 and all with phantoms.

### C22: confirmed

All 179 IDs map to variants I built independently, and all are valid, fitted = rank = 5 (engine B) and clean (exact flags). The mix is 55 seeds + 114 two-group MERGE + 6 TC + 4 ping-pong half-TC. They have 156 distinct supports. The 23 pairs that share a support are all non-EQUAL by Gröbner (basis [1] both ways), so no duplicates remain. The 18 clean 5/5 variants left out are exactly the 8 lumping twins, one member of each of the 8 C9 pairs, and 2 TC duplicates. I proved the C9 pairs EQUAL exactly (unique positive Gröbner solution both ways, 2 laws each), where the track had only SAMPLED them.

## New issues

- The V check falsely rejects clean, identifiable mechanisms when a substrate or product also binds as a dead-end inhibitor (R3/R6-type reactions), so 'rule OK ⟺ non-degenerate' and 'no component ever falsely rejects' do not generalize. Smallest example, R3_00005 merged with SS={S}: E+S→X (SS), X⇌E+P (RE), E+S⇌ESinh (RE). Its law v = (kf·S − kr·P/K)/(1 + S/Ki + P/K) is clean and 3/3, but V_fwd rejects it: the dead-end copy caps the SS binding flux, and V cannot see that. Counts: R3 merged, 2 of 4 clean variants rejected. R6 merged (88 level-1 dead-end parents, 1,320 variants), 98 of 704 clean variants rejected (52 by V_fwd, 46 by V_rev). All 100 are identifiable: 12 have one SS group (3/3 or 5/5), 88 have two (6/6). R5 (inhibitor I only): 0 of 1,460 disagree. Today's R6 export: 0 disagreements, because V only fails in merged representations. With today's moves V is monotone (RE→SS flips and pendant dead-ends cannot make it fail), so this bites only when a MERGE is applied to a mechanism that already carries S/P dead-end copies, e.g. the merge-anywhere move of points 3 and 10. Files: inh_results_R6.jsonl, vfalse_results.jsonl, v_counterexample.py.
- Rules 1 and 4 together ('reject degenerate in every move' and 'drop the ping-pong seeds') lose reachable identifiable families; 'loses: Nothing' is wrong. Moves keep the step set, only add SS flags and only refine groups, so all 132 clean ping-pong mechanisms at levels 1–2 of today's R4 enumeration descend only from the 7 degenerate ping-pong seeds (21 of them have no other potential parent at all). 83 of the 132 are identifiable: 16 at level 1 with 7/7 (e.g. R4_00386, 00388, 00389, 00395–00397) and 67 at level 2. None of them shares a (numerator support, denominator support, rank) with any of my 1,101 valid merged/ELIM variants, so no merged family without splits is EQUAL to them. Only the plain seed's children R4_00428/00429/00430 (7/6) have merged 6/6 counterparts; by Gröbner, a random law of each lies in two 3-SS merged families. Possible fixes: keep degenerate seeds as non-fitted expansion parents, or accept the loss explicitly. Files: r4_lost.json, lost_twins_18.json.
- The I (visibility) component is redundant once validity uses a semantic flux check (every SS group carries nonzero net steady-state flux, which is Denis's point 1). On every flux-valid mechanism I tested (1,163 own variants, 1,631 R4 export, 1,466 R3/R5 merged, 1,320 R6 merged), V+C alone gives the same verdict as V+C+I. The only 4 mechanisms that need I (R4_01732, R4_01744, R4_01765, R4_01781, which are the track's 'V+C has 4 false accepts') all have zero-flux SS groups. So with the point-1 flux rule, the rule can be just V+C.
- The C4 proof (§5.2) is incomplete. When the flipped group also holds an abortive step, the dead-end form is its own segment rather than a self-loop inside one segment (e.g. R4_00001 with SS={B}: EBQ). The 248/248 EXACT result still holds.
- Ordered/random seeds also have an inclusion-minimal clean three-group SS set with 6 parameters, one per seed and 28 in total (of the form {outer substrate, P, Q} or its mirror). They sit alongside the 28 random/random ones. The seed rule sends both to level 1, but the report mentions only the random/random ones.
- The C8 NEITHER verdicts and the C9 EQUAL verdicts can be settled exactly without family_fit. My contain.py does Gröbner coefficient matching on a random rational law and runs in about 0.1 s per pair for bi-bi. It avoids both family_fit failure modes that C23 describes.
