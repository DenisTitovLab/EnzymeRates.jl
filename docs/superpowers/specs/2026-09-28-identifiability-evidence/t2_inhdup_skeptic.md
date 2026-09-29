# Skeptic verdicts: Point 2: substrates and products as inhibitors

Independent recheck of the track report `t2_inhdup_report.md` with separate code.

## Summary

Denis, the track's core numbers hold up. The recommended rule ("a dead-end copy must create a new complex") does what it claims on R6 up to level 2. It is not yet a complete guard, for two reasons.

**What I reproduced independently.** I wrote my own port of the dead-end move from the Julia source and my own twin, SD and WSD tests, and used engine B where the track used engine A.
- The port matches the package: every exported child, and every raw and distinct count.
- On the exhaustive R6 regeneration (36,844 events): 3,082 phantoms, all of them all-twin; 0 of 26,292 non-twin events create a phantom; the rule removes 10,552 events (7,470 of them identifiable).
- Engine B agrees on all 7,631 exported dead-end children and on a stratified sample of 215 events.
- The non-productive-binding theorem is correct: its parameter maps make the laws identical symbolically on U1, R6_00075 and R6_00111 (a case with a mirror step).
- The foreign-inhibitor theorem holds (R5 up to level 3: dPh is never positive).
- Uni-uni copies are always phantoms (hand proof of family equality).
- The wider-family certificates for B10 and R6_03189 reproduce.
- The split-sibling counts reproduce; the "+0" is relative to the child, not the parent.
- WSD is exact on the 31,088 events whose parents bind only by RE.

**Where the recommendation needs changes.**
1. **Split move.** It can split a kept copy's group into a part whose only site duplicates an existing complex. That recreates the non-productive-binding phantom: 42 of 184 exported split children of rule-kept parents (engine B), about 285 across all of level 2 by extrapolation. The twin test must be a check on every move's children, or copy groups must never be split.
2. **Composition key.** It misses a complex that matches an existing form only through an RE isomerization. At level 3 there are 3 counterexamples to "a new complex never creates a phantom", and one has a lineage the rule keeps (R6_28303 + Q\*, 8/8 → 9/8). Also treating a site as a twin when it matches an existing form's (RE segment, net uptake) catches every phantom at levels ≤ 2 and 3, for 4 and 10 extra identifiable losses. This matters for your merge moves too.
3. **Rule 2 as written.** Its combinatorial text has a false positive: hand example H1 is identifiable at 7/7 with both engines, yet the text rejects it. The track's actual code is correct.
4. **Rule 3.** WSD misses 8 level-3 phantoms whose parents bind only by RE, so it cannot be the only filter.

Everything is in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/verify_t2_inhdup/ (key_numbers.json summarises the counts; hand/ has the counterexample mechanisms).

## Verdicts

### C13: confirmed

I wrote my own port of `_expand_add_dead_end_regulator_native` from src/mechanism_enumeration.jl:1828-1960, without using deadend_port.py. It reproduces the package exactly. Raw/distinct child counts: R2 1/2 and 1/2; R3 2/6 and 2/5; R5 350/2,153; R6 1,400/35,444 raw and 1,400/25,344 distinct. Every exported dead-end child (R2 3, R3 7, R5 2,503, R6 5,118) is regenerated from its recorded parent, with 0 misses. Files: /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/verify_t2_inhdup/vport.py and check_port.py.

### C5: refuted

The claim holds on its stated finite set. On my own regeneration (engine A), 0 of 26,292 non-twin R6 events up to level 2 are phantoms. Engine B agrees on all 5,118 exported R6 children: (dPh, all-twin) = (1,T) 436, (0,T) 966, (0,F) 3,708, (-1,F) 8. It also agrees on a stratified sample of 215 events.

As a general rule it is false. I applied the package's own move to the 4,169 exported level-2 mechanisms (level 3, 70,687 events). Three children have a new complex yet create a phantom; engine B confirms all three.

The cleanest is R6_28303 + Q* at {E, EPinh}:
- E·Q* duplicates EQ.
- EPinh·Q* has composition {P,Q}, which is new.
- Parent 8/8, child 9/8.

R6_28303's lineage (init seed R6_00058 → split → non-twin dead end) passes the rule, so this is a real miss. Cause: EPinh·Q* has the same RE segment and net uptake (P+Q) as EBP_res, because EB_res ⇌ EQ is an RE isomerization. A composition key cannot see that. Files: hand/L3_*.json and ev_R6_L3.jsonl.

### C12: confirmed

Reproduced exactly with my own twin predicate and port (engine A ranks):
- Rejecting all-twin children removes 10,552 of 36,844 events: 3,082 phantom and 7,470 identifiable.
- No phantom-creating dead-end event is left.
- At seed level it removes 360 events (124 phantom), and every one of the 248 seed×ligand pairs keeps at least 2 children.
- R3 loses all its children; R2 and R5 are unchanged.

The guarantee holds only in the stated scope (R6, levels ≤ 2). The next two new_issues show it fails once the enumeration goes further (split route, level 3).

### C3: confirmed

I checked the proof step by step. The Lemma-1 gauge with μ = log(K_g/(K_g+K*)) on C gives K_g' = K_g + K*, and constants leaving C are scaled by ρ.

I then substituted the map into engine-B symbolic laws, and the child law minus the mapped parent law cancels to exactly 0 in three cases:
- U1 vs U0;
- R6_00075 vs R6_00001 (two sites, SS step leaving C);
- R6_00111 vs R6_00003 (two sites plus an RE mirror in the Q group).

The inverse map exists for any K* < K_g', so the families are EQUAL.

My own implementation of the SD predicate selects exactly 2,392 R6 events, all with dPh = +1 (0 FP, 690 FN); in R3 it selects 6 of 8 events.

Caveat: this confirms the μ-form of the theorem. The combinatorial one-line version written into rule 2 is weaker (see new_issues).

### C2: confirmed

The I = 0 slice argument is sound. The child's fitted count is always the parent's +1: I checked this on all 36,844 R6 events, 2,503 R5 events and 70,687 level-3 events. Null(child) therefore injects into null(parent).

One implicit assumption: ∂v/∂K* must not vanish identically. The move never builds a complete inhibitor-bound layer, so this holds in practice.

Data: R5 levels 1–2 give dPh ∈ {0,-1} on 2,503 events (engine A, and engine B on all exported children). R5 level 3 gives dPh ∈ {0,-1} on 8,680 events. R2 gives 3/3 with dPh = 0.

### C15: confirmed

All 8 R3 events (7 distinct children) have dPh = +1 with engines A and B, both on the exported set and in my sample.

Hand proof of family equality for S binding by SS (U3: E+S⇌ES SS, ES⇌EP SS, EP⇌E+P RE):
- Law: v = (αS − (α/Keq)P)/(1 + γS + δP), with α = kf1·kf2/(kr1+kf2), γ = kf1/(kr1+kf2) and δ = K3 + α/(kr1·Keq).
- This map is onto the whole positive (α, γ, δ) orthant: take kf2 = α/γ, kr1 large, then K3 = δ − α/(kr1·Keq) > 0.
- Adding E·S* only shifts γ by K*, so the child's family is the same orthant.

The move can place S* or P* only at E, because sites must carry neither S nor P, so every copy duplicates ES or EP.

### C10: confirmed

Recomputed from engine-B laws, with the denominator normalized to constant term 1 and gcd(N, D) = 1 checked:
- R6_00179: c_AB − c_A·c_B = Ka_1·Ka_4, always > 0.
- R6_03189: Ka_1·Ka_5 − Ka_2·Ka_3, which takes either sign.
- B9 (rebuilt by me): c_PQ − c_A·c_B/Keq = Ka_4·Ka_5 > 0.
- B10: Ka_4·Ka_5·(kf_1 − Ka_9·kr_3)/kf_1, which takes either sign. This equals the track's expression once kr_1 is substituted.

So each child reaches laws that lie outside the closure of its parent's family.

File: verify_t2_inhdup/cert.py.

### C4: confirmed

The proof holds: after splitting groups into the classes leave C / enter C / SS inside C / rest, every group is uniform for the single-L-group gauge.

The explicit map R1 → R3 (K1' = K1 + K*, K11 = K1, K12 = K2·K1/(K1+K*)) makes the engine-B law difference cancel to 0.

Counts, recomputed independently with my own scope test and sibling construction (engine-A fitted counts):
- in scope: 8,572;
- identifiable in scope: 5,900 (4,632 / 1,156 / 112);
- non-SD phantoms in scope: 412 (384 / 28);
- SD events in scope: 2,260;
- out of scope: 410 phantom (278 once the 132 SD events are excluded) and 1,570 identifiable.

The '+0/+1/+2' counts are measured against the CHILD's parameter count, not the parent's; against the parent they are +1/+2/+3.

### C7: overstated

On the stated set it holds. My own WSD reproduces dPh exactly on all 31,088 R6 events (levels ≤ 2) whose parent is init, split or dead_end. My WSD uses group rows plus one row per (RE segment, net uptake) class, with engine-B consistent rational parameters over GF(p). On re_to_ss parents it gives TP 262, FN 224, FP 0, matching the track's C8.

The general phrasing ('when every binding step of the parent is RE, WSD predicts dPh exactly') fails at level 3. Among RE-binding parents it misses 8 phantom events, all confirmed by engine B. Examples: R6_28155 + B* at {E_res}, 9/8 vs 8/8; R6_29723 + P* at {EA, E_res}. All 8 are ping-pong children where the gauge crosses the F segment's shared B/P groups. In 6 of them the parent's WSD nullity is 0, so a finer class definition cannot catch them either.

### C18: confirmed

Timed on 300 random exported R6 dead-end children while another agent's Python jobs were running:
- engine-A Model plus rank(): 1.03 ms per mechanism, consistent with the claimed 1.2 ms;
- composition-twin lookups for every form: 0.016 ms per mechanism, a negligible set lookup as claimed.

I did not time the track's own SD or WSD code. My WSD, which builds engine-B parameters, runs about 1 ms per mechanism.

## New issues

- The split move brings back the phantom class the rule removes. `_expand_split_kinetic_group` can split the copy's single kinetic group by context, leaving a group whose only site is a twin. Example: R6_00117 (R6_00003 + B* at {EA, EQ}, identifiable 6/6, kept by the rule because EQ·B* is a new complex) splits into R6_02495, where EA→EAB* sits alone in group 4 and duplicates EAB. Engine B gives 7/6. Among the exported level-2 split children of dead-end parents the rule keeps, 42 of 184 are phantoms (engine B), and every one has a copy group whose sites are all twins. Scaling from the uniform 4,169-of-28,304 sample gives roughly 285 such level-2 mechanisms (SAMPLED estimate). Fix: apply the twin test as a mechanism-level invariant to every move's children (reject a child in which some copy-binding group has only twin sites), or stop the split move from splitting copy groups.
- The composition key misses complexes that match an existing form's concentration dependence only through a ligand-free RE step. Examples are a ping-pong second chemistry step written as RE (EB_res ⇌ EQ), and, in future, merged central complexes. At level 3 (dead-end children of the 4,169 exported level-2 mechanisms), 3 non-twin children are phantoms (engine A and B). One of them, R6_28303 + Q* at {E, EPinh} (8/8 → 9/8), has a lineage the rule keeps. Fix: also count a site as a twin when X·L* has the same (RE segment, net uptake) as an existing form, and reject if every site is a twin under either key. This union catches all 5,737 level-3 phantom events and all 3,082 at levels ≤ 2. It costs 10 extra identifiable events at level 3 and 4 at levels ≤ 2. The uptake key alone misses the SS-binding twins (48 at levels ≤ 2, 210 at level 3), so both keys are needed.
- Rule 2's combinatorial description is not the predicate the track evaluated, and read literally it has a false positive. Hand example H1 (in verify_t2_inhdup/hand/):
  - E+A⇌EA (g1), EB+A⇌EAB (g2), E+B⇌EB (g3);
  - EA+B⇌EAB and EQ+B⇌EBQ share g4;
  - EAB→EPQ (SS), EQ+P⇌EPQ, E+Q⇌EQ;
  - plus A* at {E, EB}, with an RE mirror in g3.
  Every literal condition holds: each L-binding group is site-sourced; g4 has one step inside C and one outside, so 'none cross'; the mirror is RE. Yet the child is identifiable: 7/7 with both engines. The track's own sd.py (a null-space test with compensation rows) correctly returns False. To implement rule 2 as a cheap combinatorial scan, require that a group's steps inside, leaving or entering C all attach to XL of the same L-group, or else use the null-space version. H1 is reachable by two context splits, although R6 at levels ≤ 2 contains no such case (my literal and refined versions agree on all 36,844 events).
- Rule 3 (WSD) is not exact beyond level 2: it misses 8 phantom level-3 events whose parents bind only by RE (C7). All 8 are all-twin, so rule 1 and the exact modular rank still catch them. WSD should not be used as the only filter on the claim that it 'loses nothing and is exact for RE-binding parents'.
- At level 3, 12 all-twin dead-end children have dPh = −1: the copy breaks a gauge the parent already had. The composition rule discards them although they are more identifiable than their parents. None occur at levels ≤ 2. Minor.
