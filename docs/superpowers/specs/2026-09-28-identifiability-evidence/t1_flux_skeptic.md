# Skeptic verdicts: Point 1: zero-flux SS groups

Independent recheck of the track report `t1_flux_report.md` with separate code.

## Summary

Denis, I rederived Track T1 with my own code and no Julia. I did not import any track scripts. My tools:
- an exact GF(2^61-1) full mass-action solver that takes the RE limit by interpolating polynomials in Lambda, with no Cha segments;
- my own potential-based parameter sampler and fitted count, which matches both the package and engine B on all 12,556 mechanisms;
- a zero-flux predicate defined directly as 'lies on a simple cycle of nonzero first-substrate weight';
- the package screen written from its definition;
- my own port of _expand_re_to_ss, which reproduces the raw flip-children counts for R4 (1,777), R5 (3,216) and R6 (7,533) with 0 misses.

Everything the track reports about the export reproduces exactly:
- 25,850 SS steps, 1,726 zero-flux, 0 predicate mismatches; all flux-carrying steps are one-way.
- 739 zero-flux groups in 412 mechanisms (188 distinct structures).
- The package-screen confusion is 19,862 / 654 / 85 / 0, and the screen equals the all-SS predicate on all 67,235 groups.
- Converting a zero-flux group to RE keeps v exactly equal (scaling kf and kr together leaves v unchanged) and lowers fitted by exactly 1 per group; the twin is always enumerated (for flips, the parent itself).
- Flux status: flips are monotone, splits and dead-end additions leave it unchanged, MERGE and ELIM behave as claimed (16 flux losses in 11,332 single merges; 0 changes in 40,058 ELIMs).
- The flip rule removes 148 children and adds 0 at levels 0-1 in R4, R5 and R6.
- The 64 level-3 mechanisms are reached at level 4.
- Phantom totals (2,388 to 1,619; R4 616 to 266) and the 743 of 819 mixed-only mechanisms that are identifiable are confirmed with engine B ranks.

So the core rule is sound: reject a mechanism with any steady-state group that has no flux-carrying step, checked group by group, never step by step.

Three corrections:

1. **The current package screen breaks for your proposed mechanisms.** Keeping `_flux_carrying_groups` as the flip pre-filter only works while every turnover cycle contains an isomerization step. Theorell–Chance, merged central complex and ELIM mechanisms have none, so the screen flags nothing and the flip move produces no children. My counterexamples get 0 flip children with the package screen and 2 with the correct one. The fix is to use the new cycle-weight test on the all-steady-state version of the mechanism; on the export it gives the same flags as today's screen, with less code.

2. **"Extend" is required, not optional.** On the R5 level-2 parents, the extend variant emits 8 valid pair flips (for example R5_04329, flipping g8 and g9) that the drop variant never produces from any R5 level-2 flip. Today the package reaches them only one move later, through a zero-flux intermediate.

3. **The "check only the flipped groups" shortcut needs every move to keep the invariant.** It assumes no parent already has a zero-flux group. That holds for flip, split (with the new filter) and dead-end, but the allosteric moves have not been analyzed and a single-step SS merge breaks it. Checking every steady-state group of every child costs the same and removes that dependency.

Files are in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/verify_t1_flux/:
- `vlib.py`: the solver, predicate and package screen.
- `flipport.py` and `flipport2.py`: the flip-move port.
- `vflux.jsonl`: per-mechanism results.
- `summary.json`: all the numbers above.
- The check scripts: `analyze.py`, `convert_check.py`, `moves.py`, `scaling.py`, `merge_elim_v.py`, `check_flipport2.py`, `deeper.py`, `added.py`, `added2.py`, `level3b.py`, `r6twins.py`, `screen_merged.py`.

## Verdicts

### C1: confirmed

My own method: exact GF(2^61-1) mass action with all forms kept. The RE limit comes from interpolating det and det*J as polynomials in Lambda (no Cha segments). Parameters are sampled from a potential parametrization with integer-nullspace group ties. For all 412 zero-group mechanisms: scaling kr (K fixed) of every zero group leaves v exactly unchanged at 3 random points (412/412); the RE-converted mechanism gives the same v at the same K (412/412); my potential-space count gives fitted drop = #zero groups (412/412) and single-group drop 1 (739/739). The RE twin has engine-B rank equal to the original's (412/412). Rank equality follows from v_Z = v_C composed with the projection kf/kr.

### C2: confirmed

Exact GF(p) finite-scaling test on every SS group of all 12,556 mechanisms. Scaling (kf,kr) jointly by a random factor changes v for all 19,862 flux-carrying groups and leaves it invariant for all 739 zero-flux groups (0 exceptions). The general two-way/grouped case stays a CONJECTURE, as the track says.

### C3: confirmed

Checked the potential/Thomson proof: J_e = c_e(psi_a-psi_b), the block restriction of a circulation is a circulation, and the sum of c_e*dpsi^2 = 0. Exact: all 1,726 steps my predicate calls zero have J = 0 in the exact Lambda limit at 2 random GF(p) points.

### C4: confirmed

My predicate is defined directly as 'lies on a simple cycle of nonzero first-substrate weight', by DFS over simple cycles with no block algorithm. It agrees with the independent exact fluxes on all 25,850 SS steps (1,726 zero, 24,124 nonzero, 0 mismatches). All 24,124 flux-carrying steps are one-way (0 two-way) by my own cycle enumeration. The grouped two-way case is correctly left as a CONJECTURE; it could only make the rule miss a rejection, never wrongly reject.

### C5: overstated

Reproduced on the export exactly. Using the definition (a group shares a simple cycle with an iso step), the SS-group confusion is TT 19,862 / TF 654 / FF 85 / FT 0, and the screen equals the all-SS cycle predicate on all 67,235 groups of the 10,071 distinct mechanisms. The general soundness ('never flags a flux-carrying group as zero') only holds while every unbalanced cycle contains an isomerization step, which is true for the current single-metabolite Step model. It fails for Theorell-Chance / merged-central-complex / ELIM mechanisms. Counterexample: E+A<=>EA (RE), EA+B->EQ+P (SS), E+Q<=>EQ (RE). All steps carry exact flux, the screen flags no group, and the flip move gets 0 units.

### C6: confirmed

Independent exact flux gives 412 mechanisms, 188 distinct structures (group-relabel-invariant key) and 739 zero groups. Strata are identical: R4/R5 L2 re_to_ss 148 (296 groups), split 40 (40); R6 sample re_to_ss 31 (62), split 5 (5); 0 at levels 0/1, 0 from dead_end, 0 in R1-R3. 654 of the 739 are package-flagged.

### C7: confirmed

Mapped each re_to_ss child's steps to its parent's by (from,to,consumed,released) and used exact fluxes: 4,172 parent SS steps, 0 flux->zero and 0 zero->flux. All 654 zero groups in flip children consist only of newly flipped steps.

### C8: confirmed

4,097 parent steps across split children, 0 status changes. 85 split-created zero groups, each a single step split off a flux-carrying parent group: EAP 23, EBQ 20, EBP 18, EAQ 16, EAP_res 4, EBQ_res 4.

### C9: confirmed

11,846 parent steps in dead_end children, 0 changes. All 557 new SS steps match a parent step of the same group with the inhibitor tag stripped, and none differs in flux status.

### C10: confirmed

My own MERGE with exact flux over the 10,071 distinct mechanisms: 11,332 single-step merges (6,396 give an infinite rate), 16 flux->zero, all SS (R5_04414 and similar), 0 on RE merges, 0 zero->flux, 0 predicate mismatches. Whole-group merges: 39 valid, 27 infinite, 0 flux->zero.

### C11: confirmed

My own ELIM (fused step typed SS if either input is SS, in a new group, TC steps allowed): 40,058 applications, 0 status changes, 0 predicate-vs-exact mismatches.

### C12: confirmed

RE conversion followed by a key lookup finds a twin for 412/412. 327 flip cases equal their own parent; 408 have a twin in the same reaction. The 4 R6 split cases (R6_05800/07918/12554/15966) have twins R4_00655/00761/01004/01179. My flip port applied to their R6 level-1 copies (R6_00369/00510/00850/01061) generates each twin (4/4). Exact rational maps are replaced by exact GF(p) identity checks.

### C13: confirmed

My own Python port of _expand_re_to_ss reproduces the raw flip-children counts exactly for R4 (1,777), R5 (3,216 = 255+2,961) and R6 (7,533 = 255+7,278 including the unexported R6 level 2), with 0 misses. On R4, R5 and R6 level-0/1 parents both the extend and drop variants remove exactly 148 and add 0. This holds only at level 0/1 and for R4 level 2; see new_issues for R5 level-2 parents, where the variants differ.

### C14: confirmed

Flip-only simulation on R4: from 1,388 level-2 parents there are 3,211 new level-3 flip children, 2,555 flux-valid, and 64 of those are reached only via zero-flux level-2 parents. Starting from the 1,200 valid level-≤2 mechanisms with rule (extend) flips, all 64 appear at level 4. The other 2,491 valid level-3 children are still reached at level 3.

### C15: confirmed

819 mechanisms have zero-flux steps only inside flux-carrying groups. 743 of them are identifiable by engine B (rank_b == fitted_b; 0 ambiguous ranks); the track used engine A.

### C17: confirmed

With engine B ranks: all 12,556→12,144 mechanisms, phantoms 2,388→1,619, mechanisms with phantoms 1,909→1,497. R4: 1,819→1,631, 616→266, 418→230. R4 L2: 1,388→1,200, 579→229, 383→195.

### C16: unverified

Not rerun. The SAMPLED family_fit is superseded by the exact GF(p) identity v_Z = v_C∘pi, which holds on 412/412 (C1/C12).

## New issues

- The recommendation to keep _flux_carrying_groups as the unit pre-filter is wrong for mechanisms whose turnover cycles have no isomerization step: Theorell-Chance, merged central complex, and ELIM products (Denis's points 3-8/11). The screen only flags groups that share a block with an iso step, so it flags nothing in such mechanisms and the flip move gets no units. Counterexamples (exact flux: every step carries flux):
  - TC: E+A<=>EA (RE g1), EA+B->EQ+P (SS g2), E+Q<=>EQ (RE g3). The package-logic port gives 0 flip children; the all-SS cycle-weight screen gives 2 (flip g1, flip g3).
  - E=EA=(X)=EQ=E: E+A<=>EA (RE), EA+B->X (SS), EQ+P->X (SS), E+Q<=>EQ (RE). Again 0 vs 2 children.
  Fix: replace the pre-filter with the new cycle-weight helper run on the all-SS mechanism. It equals the package screen on all 67,235 exported groups, and it is less code.
- 'Extend' versus 'drop' in the flip rule is not a free choice (the track's open question 3), and 'adds no children' holds only at levels 0-1 and for R4.
  - On the complete R5 level-2 parents (3,934), extend emits 14,580 children and drop emits 14,572.
  - Extend adds 8 valid pair flips that the package does not emit today, e.g. R5_04329 flipping g8 {EIinh_res+B->EBIinh_res} together with g9 {EQ+B->EBQ, E_res+B->EB_res}. Also R5_04334, 04336, 04338, 04400, 04403, 04687 and 04689, all ping-pong split children with I.
  - Today these are reached only one move later, through a zero-flux bridge (flip g9 alone, whose g9 carries zero flux).
  - Under drop, no flip of any R5 level-2 mechanism produces them.
  The rule must use extend (count the set as failed).
- Checking only the flipped groups relies on an invariant: no parent has a zero-flux SS group. That holds only if every move preserves it. Flip, split (with the filter) and dead-end are verified. Not analyzed: the four allosteric moves (to_allosteric, change_allo_state, merge_regulatory_sites, promote_catalytic_to_onlya) and allosteric A/I graphs. A single-step SS merge breaks the invariant (16 cases on the export). Suggest running the O(steps+forms) check on all SS groups of every emitted child, or at least a test asserting the invariant over all children.
- The R6 numbers in the track come from the level-2 sample. On R6's complete level-2 flip set (7,278 raw children, rebuilt by the port), the flip rule removes 148 distinct children, the same structures as in R4, not 31.
- Supporting point, with no change needed: split and dead-end children of a zero-flux mechanism are always zero-flux themselves (step statuses are preserved), so rejecting a zero-flux mechanism can only cost flip descendants. For R4 that is the 64 level-3 mechanisms, all reached at level 4 under the extend rule. Beam pruning of the extra intermediate is the only remaining risk.
