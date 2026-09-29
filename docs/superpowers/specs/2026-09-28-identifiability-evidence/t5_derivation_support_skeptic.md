# Skeptic verdicts: Point 11: derivation support

Independent recheck of the track report `t5_derivation_support_report.md` with separate code.

## Summary

Denis, I checked 10 of the T5 claims with my own code: a BigFloat mass-action steady state (RE steps at rate scale 1e60, thermodynamically consistent parameter values) compared with rate_equation on 612 mechanisms. The track's facts reproduce. One design conclusion does not hold.

**Confirmed (6 of 10):**
- **C1:** Theorell–Chance loses P in the DSL.
- **C3:** an RE step that releases its metabolite hits a BoundsError.
- **C4:** wrong stoichiometry, then a thermodynamic error.
- **C8:** the same fused step written two ways gives mechanisms that are not ==.
- **C11:** `_kcat_forward` emits `max()` with no arguments.
- **C12:** the modular rank prototype matches engines A and B on all 12,534 mechanisms. I re-ran it on 431 R4 mechanisms at about 3 ms each, and all agreed.

**C7 is confirmed and wider than stated.** The spurious "no reaction cut" hits 15 of 255 merged random bi-bi flag patterns. One of them is the merged form of the random RE seed.

**Overstated (3):**
- **C5** ("no orientation works under both flags, so a merge move can't be bolted onto the current Step"). It holds only when the merged complex is kept as the substrate-side species (EA + B → EAB, then a fused release EAB → EQ + P). Keeping the product-side species instead (EA + B → E(P,Q), then the pure release E(P,Q) → EQ + P) needs no fused release. That form derives exactly today:
  - uni-uni: 3 of 3 RE/SS combinations;
  - ordered bi-bi: 15 of 15;
  - ping-pong: 15 of 15;
  - random bi-bi: 161 of 176 patterns with a finite rate. The other 15 are the C7 cut bug.

  Fitted counts are unchanged (Segel 7, ping-pong 7). The atom check passes, and the dead-end move keeps all 16 patterns. So MERGE does not require the Step redesign. The redesign is needed only for Theorell–Chance and fused steps with two metabolites. There is also a corner case where the "PROVED" claim fails outright: a fused release with no bound metabolites on either side works under both flags.
- **C6:** "exactly −v" holds only when the cut contains nothing but reversed steps. A cut that mixes a reversed and a forward step gives a law that is silently wrong and not a sign flip.
- **C9:** every listed failure reproduces for the substrate-side form. In the product-side form, only the `is_iso`-based checks still block: the flux-carrying test (so RE→SS gives no children), the allosteric moves, and the chemistry assertion for covalent ping-pong steps.

**New hazards:**
- A fused binding and a pure binding of the same metabolite at the same form get the same parameter name, silently (kon_B_EA appears twice in fitted_params).
- 31 of 255 merged random patterns are rejected at construction by the bottomless-segment check.

Everything is in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/verify_t5_derivation_support/:
- scripts: masscheck.jl, run_cases.jl, run_random.jl, corner.jl, probes.jl, rank_repro.jl, cmp_modrank.py
- results: results_all.jsonl, results_rndPS.jsonl, results_rndSSo1.jsonl, probes.log, probes2.log, rank_repro.tsv

## Verdicts

### C1: confirmed

Reproduced in my own probe (probes.log). The DSL line `E(A) + B <--> E(Q) + P` is stored as EA -> EQ with bm=B, SS, and P is dropped. fitted_params then errors with 'Cycle 1 produces metabolite change not proportional to net reaction'. `E + A + B <--> E(A, B)` fails at macro expansion with 'step side has more than one metabolite term (A, B)'. fieldnames(Step) = (:from_species, :to_species, :bound_metabolite, :is_equilibrium). By reading, dsl.jl:744-745 takes the LHS metabolite, else the RHS one.

### C3: confirmed

Method: my own BigFloat mass-action steady state (RE steps at rate scale 1e60, thermodynamically consistent potentials), compared with rate_equation.
- U3, U8, B5 and P4 crash with 'BoundsError: attempt to access 0-element Vector{Symbol} at index [1]'. Their fitted counts are 2, 2, 6, 6, which equal a hand count of (#params − 1 Haldane) for these single-cycle mechanisms.
- Every substrate-side pattern with an RE fused release crashes: ordered 7/7, uni-uni 1/1, ping-pong 11/11, random bi-bi 130/130. This includes one where the RE release closes an RE cycle.
- Root cause read at rate_eq_derivation.jl:295/305 (m_l[1] indexed unconditionally).

### C4: confirmed

U4, UB2 and B4 reproduce the error 'Cycle 1 produces metabolite change not proportional to net reaction'.
- Stored steps: U4 is E->ES bm=P SS; UB2 is EA->EQ bm=P RE; B4 is EQ->EAB bm=P SS.
- Exhaustive check: the ordered fused release written `E(Q) + P op E(A, B)` fails in exactly the 8 of 15 RE/SS patterns where that step is SS. The uni-uni analogue fails 2/2.

### C5: overstated

The claim holds when the fused release's source form carries bound metabolites. Checked over all RE/SS combinations:
- Ordered: `E(A,B) op E(Q)+P` fails iff the step is RE (7/15); `E(Q)+P op E(A,B)` fails iff it is SS (8/15).
- Same split for uni-uni and ping-pong.

The universal 'PROVED' statement is false:
- A fused release with no bound metabolite on either side works under both flags when written `E + Q <--> E(; residual = A - P)`. This is uni-bi ping-pong, F -> E + Q: exact under SS (fitted 5) and RE (fitted 4). Written `F <--> E + Q`, it fails under both.

The design conclusion built on C5 ('a merge move cannot be bolted onto the current Step') is refuted. Keep the product-side species as the merged complex, e.g. E(A) + B -> E(P,Q) then a pure release E(P,Q) -> E(Q) + P. That mechanism has no fused release, and it derives exactly with today's Step under every RE/SS combination tested:
- uni-uni 3/3;
- ordered bi-bi 15/15;
- ping-pong 15/15;
- random bi-bi 161/176 finite, non-bottomless patterns; the other 15 are the C7 cut bug.

Worst relative error was 1.2e-56 against the mass-action rate.

### C6: overstated

Reproduced the sign flip:
- U6c: relative error vs −v is 1.5e-60; vs +v it is 2.0.
- B6c: vs −v 4.6e-59.

'Exactly −v' holds only when the cut consists solely of reversed steps. Counterexample rndPS_rev4_SSSS: merged random bi-bi, all SS, with `E(P, Q) <--> E(B) + A` written reversed. The chosen 2-step A cut mixes one forward and one reversed step. The law is neither +v nor −v (relative error 1.75 vs v, 1.84 vs −v), so the error is silent and not detectable by sign. The same mechanism with step 3 reversed happened to pick a clean cut and is exact.

### C7: confirmed

RB2 reproduces 'rate_equation: no rapid-equilibrium-consistent reaction cut'. My mass-action rate is finite: v at scales 1e40 and 1e50 agree to 1e-5.

The bug is broader than one case:
- Product-side merged random bi-bi: 15/255 patterns have a finite rate and fail this way. All 15 have an all-SS produce-X or consume-X cut.
- Substrate-side: 8/255 fail this way.
- The failures include RRSSRRRR (all RE except the two fused substrate steps), the merged form of the random RE seed.

Root cause confirmed by reading: central_forms is built only from `r.typ === :chem` endpoints (rate_eq_derivation.jl:572-575).

### C8: confirmed

U5 (`E + S <--> E(P); E(P) <--> E + P`) vs U6 (`E(P) <--> E + S; E + P <--> E(P)`):
- Mechanism == is false and typeof == is false.
- Fitted names are (koff_S_E, kon_P_E, kon_S_E) vs (koff_S_EP, kon_P_E, kon_S_EP), exactly as claimed.

The same holds for a product-side fused binding written reversed: B6 vs B6b are not ==, with names kon_B_EA vs kon_B_EPQ. Constructor code (types.jl:160-170) swaps only when bm is in bound(from) and not in bound(to).

### C9: overstated

Every listed fact reproduces for the substrate-side merged form (probes2.log):
- Atom error text is identical: 'atom-non-conserving step EAB → EQ (bound P): Δatoms Dict(:C => -1) ≠ Dict(:C => 1)'.
- Flux flags are all 0; RE→SS gives 0 children vs 4 for the unmerged ordered RE seed, and 0 vs 2 for uni-uni.
- expand_mechanisms errors.
- Dead-end under R6 gives 12 vs 16 children; the 4 lost are exactly the patterns that include EQ.

The headline 'does not accept merged mechanisms' is false for the product-side merged form, e.g. E+A⇌EA, EA+B⇌EPQ, EPQ<-->EQ+P, EQ⇌E+Q:
- the atom check passes;
- _assert_chemistry_is_iso passes (non-covalent);
- expand_mechanisms returns 9 children (merged uni-uni: 1);
- dead-end gives all 16 patterns, identical to the unmerged seed.

What remains broken there:
- _flux_carrying_groups uses is_iso, so RE→SS still gives 0 children;
- the allosteric moves use is_iso (read at mechanism_enumeration.jl:2007/2199/2240, thermo:565);
- _assert_chemistry_is_iso rejects covalent ping-pong fused bindings such as E + A -> E(P; residual = A - P).

### C11: confirmed

ER._kcat_forward raises 'MethodError: no method matching max()' on U4b and U6c. It returns numbers on U1 (2.0) and on product-side uni-uni RS (12.0). By reading (rate_eq_derivation.jl:919-921), an empty candidate list emits Expr(:call, :max) with no arguments.

### C12: confirmed

My own comparison script (cmp_modrank.py) against the track's TSVs:
- The TSVs have 1819/4715/6000 rows.
- Rank and fitted equal engine A (cv/results/export_a.jsonl) on 12,534/12,534.
- Rank also equals engine B (export_b.jsonl, rank_b) on 12,534/12,534.
- Phantom counts are 418/567/906; median t_rank is 3.5/4.5/6.2 ms.

I also re-ran modrank.jl myself on my own enumeration of R4 levels 0-1 (431 mechanisms = 62+369):
- All 431 agree with engines A/B. I matched by printed text plus fitted count, because 106 export texts are shared by two ids that differ only in grouping.
- Median 2.9 ms, p90 5.1 ms, max 0.10 s.

Caveats:
- The prototype reads the package's own num/den and Haldane dependents, so it is a fast rank of the package's law, not an independent derivation.
- I did not re-measure the compile-time comparison (C13).

## New issues

- The track concludes the merge needs a Step redesign, but it only tested keeping the SUBSTRATE-side species. That representation needs fused releases, and those fail under one flag or the other (C5). Keeping the PRODUCT-side species needs no fused release: EA + B -> E(P,Q) is a fused binding with B consumed, and E(P,Q) -> EQ + P is a pure release. It derives exactly today under every RE/SS combination (verified against my mass-action steady state): uni-uni 3/3, ordered 15/15, ping-pong 15/15, random bi-bi 161/176 finite non-bottomless patterns (the other 15 are C7). Fitted counts equal the substrate-side ones wherever both derive, e.g. Segel 7 and ping-pong 7. Parameter names come out natural (kon_B_EA, kon_P_EQ). So MERGE can ship without redesigning Step. The redesign is needed only for Theorell–Chance and ELIM, which put two free metabolites in one step. The track's recommendation to keep the substrate-side species is the choice that fails today.
- The spurious 'no rapid-equilibrium-consistent reaction cut' (C7) is a family, not one case: 15/255 product-side and 8/255 substrate-side merged random bi-bi flag patterns have finite rates but fail. It includes RRSSRRRR (RE bindings to E, SS fused substrate steps, RE releases), the merged form of the random RE seed, i.e. the child a MERGE of that seed would produce. Every failing product-side pattern has an all-SS produce-X or consume-X cut on the merged form, so the numerator-cut fix is needed in either representation.
- 31/255 merged random bi-bi patterns (same count in both representations) are rejected at Mechanism construction by the bottomless-RE-segment check, for example segment {EA, EPQ, EB} or {EP, EPQ, EQ}. The mass-action rate is finite. This is existing policy, not a merge bug, but a MERGE move on random parents will produce such children, and the move must skip or handle them.
- Parameter names can collide silently. A fused binding and a pure binding of the same metabolite from the same form, both SS (EA + B <--> E(P,Q) and EA + B <--> E(A,B)), both render as kon_B_EA/koff_B_EA. fitted_params then lists kon_B_EA and koff_B_EA twice, with no error. Any fused-step naming keyed on (metabolite, from-form), which is what the product-side representation uses today, needs a uniqueness assertion. The track's redesign avoids this by naming fused steps by form pair.
- A reversed fused substrate step inside a multi-step cut gives a law that is neither +v nor −v (rndPS_rev4_SSSS: relative error 1.75 vs v, 1.84 vs −v). Checking the sign of v at saturating substrate therefore cannot detect it. Whatever representation is chosen, the Step constructor must orient fused substrate steps canonically (substrate consumed). Today it swaps only pure bindings.
- The C5 'no orientation works for both flags' result has an exception. A fused release whose endpoints have no bound metabolites, such as uni-bi ping-pong F -> E + Q, derives under both flags when stored as E + Q -> F: exact under SS (5 fitted) and RE (4 fitted). _step_sides branch 3 is skipped when both bound lists are empty, so branch 5 marks the metabolite as consumed. The same step written F -> E + Q fails under both flags.
- Under the product-side representation, three enumeration blockers remain:
  - _flux_carrying_groups (is_iso), so RE→SS gives 0 children;
  - the allosteric moves' is_iso chemistry detection;
  - _assert_chemistry_is_iso for covalent (ping-pong) fused bindings such as E + A -> E(P; residual = A - P).
  The atom check and the dead-end move already work (16/16 patterns).
