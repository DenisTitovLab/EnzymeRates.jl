Denis, the short answer: letting mechanisms through and leaving the rest to the loss works for what the enumerator makes. It fails for three classes that only hand-written mechanisms produce. The enumerator does make one more class, face splits and traps: their laws look plausible, but their value at zero products or zero substrates is ambiguous.

- **Enumerator output** had no v ≡ 0 law and no 0/0 face: 0 of 4,612 parents and 0 of 12,967 flip/split children. The seed screen excludes RE turnover, and the move that turns RE steps steady-state cannot create it (argued, not computed). Face splits and traps make up 0.15–0.4% of level-1 children and 0.6–0.8% of level-2 children.
- **The loss catches** v ≡ 0 (a flat 100, or 200 in absolute mode) and product activation (104.8 over 8 points). Today both still run a full fit, then crash in the kcat rescale.
- **The loss misses four classes:**
  - RE turnover beside a steady-state route. 453/flip2 prints an ordinary saturating law, but the oracle's rate grows with RE speed.
  - Laws that are 0/0 on a face, whenever the data carry products (loss 0.013–0.19 at P = Q = 0.01).
  - Laws that never saturate, fitted to data that don't saturate either.
  - Face splits and traps.
- **kcat = 0 laws are not nonsense.** They are complete substrate inhibition. Keep them.

## E2/E9: kcat classes and the gauge

**Answer**
- **kcat = 0 mechanisms are feasible.** v > 0 holds exactly at every positive substrate point for all 284 such laws, the allosteric (MWC) ones included, and for 221 of 221 in a sample of parts 5–11. v falls to 0 only as one "culprit" substrate grows. On the fitted full law, 127 of the 257 plain laws have two or three such substrates, not "occasionally two".
- **How common each class is** (all counts are on the slice's laws, not the package's):

  | Class | Where and how often |
  |---|---|
  | kcat = 0 | 7.4% of ter-ter init (parts 1–4); about 8.4% in a sample of parts 5–11; 3.7% of expand_terter; 6 of 4,786 core laws, all MWC; no plain bi-bi |
  | Infinite | no enumerated law; 13 hand-built or adversarial; 120 of 807 SS→RE or mixed toggles with a forward rate; 0 of 398 RE→SS toggles |
  | No unique limit | about 2.0% of ter-ter (parts 1–4; 2.8% over all parts per the page); 8 of 1,512 expand_terter; PFK |
  | Rate above kcat | at least 3.1% of laws with a finite kcat (a lower bound) |

  The catalytic-core kcat also misdescribes how the fitted law saturates in 12.9% of modifier-free laws.
- **A gauge over a finite reference set is always defined when Q ≢ 0,** as long as the set holds one strictly positive point away from equilibrium. Q and D have positive coefficients, so v = H·Q/D cannot vanish there. It fails only at face points: 0 for specs/24 and adv/14, 0/0 for fuzz/262. The probe's "no failures" could not have failed on its population. The real failures are v ≡ 0, which gives a gauge of 0 (68 of 875 toggles), and laws that never saturate.
- **E9's no-anchor gauge (the peak over all concentrations) is principled but not simple.** It does handle kcat = 0, overshoot and direction dependence. But for a hyperbolic law the peak lies at infinite substrate, so it needs limits or an expanding search, and a fixed box underestimates it. It is infinite for laws that never saturate.

**Uncertain**
- Parts 5–11 were only sampled.
- MWC traps and limits were checked only on the slice's own law.
- Nothing ran on the package. Today's `_kcat_forward` is not a limit: it takes the maximum over product-free monomial patterns found in both numerator and denominator (src/rate_eq_derivation.jl:610–640). Reading the code, for complete substrate inhibition it likely returns the apparent kcat, not 0, and throws only when no pattern matches. This has not been run.

**Recommendation**
- Never reject kcat = 0.
- Adopt E9's anchors and drop the saturating_core option.
- Make the no-anchor scale always succeed. When the peak is unbounded, keep the fit and leave the constants unscaled with a flag. In relative mode the loss cannot see the scale, so the scale must never turn a fit into a FitFailure. This also removes the E2/E9 conflict listed below.

## E2: can one law-based test replace the guards?

**Notation.** The slice writes v = E_t·H·Q/D, so Q is the positive numerator factor. The intended tests are (a) D ≢ 0, (b) Q ≢ 0, and (c) D ≢ 0 on both faces (all products zero; all substrates zero). Running the face test on Q, as the question words it, would reject product activation, which E2 allows.

**Answer: no.**
- **What the law tests catch:** every bottomless segment; all-RE mechanisms, through (b); 1,027 accepted v ≡ 0 laws; 22 accepted laws that are NaN on a face, with limits that depend on direction. (a) never fired.
- **What they miss:**
  - RE turnover beside a steady-state route. The package's structural check catches all 1,654 cases; the slice has no such check.
  - A declared reactant that no step uses.
  - :OnlyA admission.
  - Face splits and traps.
  - Non-unique steady states on a face (2 fuzz records; nothing catches these).
  - The corner where products and modifiers are both zero. fuzz/262 (p2_fuzz_dump_flip/4) is 0/0 there but passed (c); all 27 of (c)'s hits were mutants. This is my inference from both tracks.
- **What they over-reject:** 6 bottomless mutants whose 0/0 is removable (probe only). Running (b) per conformation would reject 5,072 of 6,555 allosteric records, so it must run on the combined law.

**Uncertain**
- That (b) is equivalent to "a steady-state step closes a turnover cycle".
- The face check's 6 false alarms and 2 misses.
- Faces with only one product, or with modifiers, at zero.
- Each verdict rests on one parameter point per record.
- Ter-ter children and children of the dead-end and tagging moves were not checked.

**Recommendation.** Take the page's minimal_guard with two changes:
1. Reject v ≡ 0 using (b) on the combined law. It never fires on enumerator output; today such a law costs a full fit, crashes the rescale and gives any gauge a value of 0.
2. Also run (c) with every modifier at zero.

Keep stoichiometry and unused reactants, RE turnover (port the package's check to the slice) and B1 admission as structural checks. (c) can replace the bottomless check.

Allow face splits and traps and document them, but fix the page's description:
- **Size:** with RE speeds spread over 10^±3, the gap is median 2.3%, 90th percentile 47%, maximum 314%. On a face with only one product at zero it reaches 80% for the hand-built LDH-like case.
- **Reading:** for traps, the law is the exact RE limit when the RE step feeding the trap is the slowest. For splits, it lies strictly between the limits of the different speed orderings.

Non-saturating laws (condition V) and laws where a reactant drops out (condition C) can go to the loss once the scale always succeeds. No track counted how many seed variants `_degenerate` removes today, so the cost in extra fits is unknown.

## E10: flip at the free form E

**Answer**
- **It is exact.** An independent oracle covered 12 shapes, n = 1–3 and 208 points; the worst error was 3.2e-69. With a single flip form, the flip edge is a bridge, so it carries no net flux and only L matters.
- **What it fixes:** iso mechanisms, where free forms are joined by a steady-state step and PR #86 rejects them. For free forms inside E's RE segment, PR #86's law is already the same family (L ↦ L·(K_I/K_A)^n, 16 of 16 comparisons); only L's label changes.
- **No effect on enumerated mechanisms.** 9,308 distinct mechanisms all have a single free form, and 300 + 223 laws came out identical.
- **It brings three new problems:**
  1. Choosing the form by the name E rejects documented names: the E_c fixture errors in the suite, and swapping names changes the law by 21–100%.
  2. It admits product-activated laws, on which `_kcat_forward` throws (24 of 24 runs).
  3. It admits R2–R4 and S1 by deleting the stranded forms, which silently decides B1.
- **The :OnlyA island cases:**
  - **S1:** deletion is the limit under every reading.
  - **R3 and R4:** if :OnlyA means a kinetic block, the limit is the EqualAI or NonequalAI twin for every ε (R4 is off from deletion by 1.4e-2 to 5.2e-1). B1's recommended rule rejects them; the prototype admits them.
  - **R2:** deletion is the limit at every P > 0. At P = 0 it depends on how fast the :OnlyA steps vanish — the same face ambiguity as E2's traps. This is my inference.
  - **Reading B's justification is refuted.** For `F + P <--> E(P) :: OnlyA`, deletion corresponds to K_I → 0, while K_I → ∞ drives the rate to 0. Its rapid-equilibrium variant passes B3, so B1's island rule must check orientation. That is untested.

**Uncertain:** regulatory sites with two free forms; n > 3; whether swapping names changes the law family, not just its values; fitting and identification on the newly admitted laws.

**Recommendation: not in PR #86.** PR #86's rejections produced no wrong law in these tracks. Do it in the rewrite under E10 store_node, keeping flipE's alpha factor (which places the flip at the designated form inside its RE segment) and its oracle tests:
- Within one RE segment, pick the flip form by content (B3); across segments, require an explicit DSL line; never pick by the name E.
- E10's "no derivation change later" holds only if the rewrite flips at the stored node, with that weighting, from day one.
- Land it after B1 is decided and after E9's scale handles product activation.

## Overclaims and contradictions
1. **kcat probe: "the other 2,566 plain laws never overshoot".** Its own verify2.log has adv_pkg/12 at 12,300×. The verifier adds p2_adv_pkg/11 and two branched hand-built cases with no abortive complex (42,960× and 3,093×). The peak equals kcat only when every turnover passes one complex reached by substrate bindings alone, with no substrate-dependent step after it. Abortive complexes and parallel routes both break that.
2. **kcat verifier's proof.** "v = k_chem·[E(A,B)] ≤ k_chem·E_t = kcat" fails when product release is steady state, because then kcat < k_chem. The conclusion still holds: every turnover passes E(A,B), and finite substrate only adds time upstream of it.
3. **Guard probe calls face laws "wrong" and cites LDH's skipped ODE test.** The verifier corrected both: a trap's law is one speed ordering's exact limit, and every allosteric spec skips the ODE test.
4. **Guard verifier's non-saturating example is condition V.** The 7 of 264 ping-pong seeds that `_degenerate` keeps from being fitted are condition C: a reactant drops out at zero products (src/mechanism_enumeration.jl:764). Those seeds are still expanded.
5. **flipE probe.** It claimed reading B is the existing K_I → ∞ convention, which the verifier refuted. Its oracle reused the package's constant names. Its "adopt now" missed the kcat throw.
6. **Decision page:**
   - E9 says kcat is "below the peak for steady-state random-order laws". That is false for the textbook random-order bi-bi (0 of 300).
   - E2 understates the size of the face gap.
   - E10 says "exact at one subunit"; it is now shown exact up to n = 3.
   - E2's minimal_guard allows v ≡ 0 and laws that never saturate, but E9's no-anchor peak is 0 or infinite for exactly those.

## Unverified numbers
- The census counts for parts 1–4 (96,218 / 7,116 / 1,900 / 5,844) and the 1,612/288 split.
- Overshoot ratios "up to 3,706".
- The guard's timings (about 10 µs per test against about 1 ms per derivation).
- The face check's 286/288 recall; the verifier confirmed it only on plain enumerator children (20/20 and 31/31).
- Agreement between the bare and reduced laws, and `admit`'s 834 rejections.
- The 81 :OnlyA ping-pong mechanisms (taken from B1).
- The test-suite counts, which come from a log only.