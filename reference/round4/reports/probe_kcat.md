## kcat saturation census: findings

The hypothesis holds. All 284 laws with kcat = 0 run forward: the rate is positive at every finite substrate concentration with products at zero. It falls to 0 only when one substrate grows, occasionally two. For the 257 plain laws the enzyme is trapped in an abortive complex that holds that substrate plus a product. Allosteric cases trap a different complex. The reference-rate gauge was defined for every law I tested.

About 30 minutes of Julia, one process at a time under the lock. Everything is in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/`, called `K/` below.

**Method in brief:**
- **Derivable** means `check_boundary`, `admit`, `solve` and `rate_law` all succeed.
- **kcat class** comes from `S.kcat_candidates`. A law counts as "zero" when every candidate's numerator is identically zero. I also classified the whole law (`kcat_terms` with every edge), called "full law" below.
- **Parameters** are thermodynamically consistent: `random_values` from `integrate/lib.jl` draws independents and Keq log-uniform in 1e-2..1e2, and dependents follow from the solve. The oracle gets the same constants through `step_constants`/`pop_constants`.
- **Numeric law** agrees with the slice's BigFloat `law_value` to 1.8e-15. At zero products it agrees with the BigFloat mass-action and MWC oracles to 1.3e-15.
- **Ter-ter sample**: `kcat_zero.jl`'s 3,000 laws, reproduced exactly; it gives the same 222 zeros. I added 300 random laws from the 6,468 "depends on how A, B, D grow" laws.

### Task 1: kcat class per set

| set | records | derivable | finite | 0 | infinite | no unique limit | full law: finite / 0 / no unique |
|---|---|---|---|---|---|---|---|
| 27 core sets | 4,917 | 4,786 | 4,766 | 6 (allo_bibi 4, AD3b, AD3d) | 13 (hand-built or adversarial) | 1 (PFK) | 4,462 / 311 / 2 |
| expand_terter | 1,926 | 1,512 | 1,448 | 56 | 0 | 8 | 1,159 / 242 / 111 |
| init_terter, kcat_zero sample | 3,000 | 3,000 | 2,778 | 222 | 0 | 0 (excluded by how the sample was drawn) | 2,162 / 637 / 201 |
| init_terter, "depends" sample | 300 | 300 | 0 | 0 | 0 | 300 | 0 / 14 / 286 |
| init_terter parts 1–4 (full census) | — | 96,218 | 87,202 | 7,116 (7.4%) | 0 | 1,900 | 79,642 / 12,960 / 3,616 |

- **No "other" errors** in the slice's kcat. The full law has one (fuzz/262: 0/0 at zero products).
- **Census agrees with the pipeline**: its "no unique limit" verdicts match the pipeline rows on all 96,218 laws.
- **Census stopped after part 4.** Parts 5–11 (about 135,000 laws) would have taken another 60–80 minutes, and another track was waiting on the lock. Those parts are covered only by the samples.
- **Not derivable**: 131 core records (130 rejected by admission, 1 cycle error) and 414 expand_terter records (226 admission, 188 spanning-tree cap).
- Logs: `K/logs/sum_all.log`, `analyze_*.log`, `sum_census.log`.

**Side finding: kcat ignores some substrate inhibition.** `catalytic_core` drops abortive complexes that sit on a leaf of the graph. When it does, kcat is finite but the full law goes to 0 at saturation. That happens in 906 of 8,992 finite-kcat laws (init_bibi 32/239, init_bibi_pp 38/264, expand_bibi 103/684, p1_enum_inh 124/1041) and in 5,844 laws of census parts 1–4. Example, init_bibi/57: kcat 18.98, full-law limit 0, highest rate found 0.246.
```
E + A ⇌ E(A);  E + Q <--> E(Q);  E(A) + P ⇌ E(A, P) [dead-end];  E(A) + B <--> E(P, Q)
E(Q) + B <--> E(B, Q) [dead-end];  E(Q) + P ⇌ E(P, Q)
```
So kcat = 0 is the case where the same kind of complex sits on a cycle, which the filter keeps. When it sits on a leaf, kcat is finite but the mechanism never reaches it.

### Task 2: the 284 laws with kcat = 0

**The rate is positive everywhere.**
- All 284 laws give v > 0 at all 20 random substrate points (1e-4..1e4) at each of 3 parameter points.
- For the 257 plain laws this is exact, not sampled: Q at zero products keeps at least one term, and every coefficient of Q and D is positive.
- The "kcat = 0 means the mechanism cannot catalyze" reading is therefore wrong. For example, allo_bibi/1023 reaches up to 0.26.

**Which substrate drives v to 0.** I compared the degree in each substrate of numerator and denominator at zero products. It is exact for the plain laws.
- In every plain zero law at least one substrate, grown alone, sends v to 0 like 1/X.
- On the law kcat uses (dead-end branches removed), that is one substrate in 240 of 257 laws and two in 17.
- Census parts 1–4: one substrate in 6,288 laws, two in 828, out of 7,116.
- No substrate ever makes v grow without bound.
- The numeric growth scans never contradict the exact degrees. In 67 of 1,188 cases, growing by 1e8 was not enough to show the fall.

**What holds the enzyme (plain laws).** I compared steady-state form populations from the oracle at 1e7 and 1e8 times the substrate.
- For the substrate that drives kcat to 0: all 672 traps are a complex holding that substrate plus a product or a residual.
- At zero products, their only exit is releasing that substrate. Their other exit binds a product, and that product rebinding is what closes the cycle and keeps the complex in the catalytic core.
- The other 280 traps belong to substrates that inhibit only the full law. They are the same kind of complex on a leaf, with no product-binding exit.
- One of the 257 plain laws had no trap found at any point (asymptote beyond 1e8).

**Allosteric laws are different (27 laws).** They are still dead ends at zero products, but the trapped complex holds no product.
- AD3b and AD3d trap a second copy of the substrate in the active conformation, `E(S, S::Mod)`.
- allo_bibi (4) and 21 expand_terter laws trap the substrate complex in the inactive conformation, where chemistry is OnlyA. Example: inactive-state `E(A, B)`, with v falling as B⁻² (n = 2).

Example 1, init_terter.part7/4757 (the same mechanism as the slice's example). Exact degree gaps for A, B, D are (−1, 0, 0). Growing A alone gives −1 per decade, B and D alone give 0. The trap is E(A, P), whose product exit is `E(A, P) + R ⇌ E(A, P, R)`.
```
E + B ⇌ E(B)  [g1];  E + D ⇌ E(D)  [g2];  E + P <--> E(P)  [g3];  E(B) + D ⇌ E(B, D)  [g2]
E(D) + B ⇌ E(B, D)  [g1];  E(D) + R ⇌ E(D, R)  [g4, dead-end];  E(P) + A <--> E(A, P)  [g5]
E(P) + R ⇌ E(P, R)  [g4];  E(B, D) + A <--> E(P, Q, R)  [g5];  E(A, P) + R ⇌ E(A, P, R)  [g4]
E(P, R) + A <--> E(A, P, R)  [g5];  E(P, R) + Q ⇌ E(P, Q, R)  [g6]
```
Example 2, init_terter.part7/5384. The trap is E(A, P), with exit `E(A, P) + R ⇌ E(A, P, R)`. v ranges 1.4e-6 to 13.5 over the random points.
```
E + B <--> E(B);  E + D ⇌ E(D);  E + P <--> E(P);  E(B) + D ⇌ E(B, D);  E(D) + B <--> E(B, D)
E(P) + A ⇌ E(A, P);  E(P) + R ⇌ E(P, R);  E(B, D) + A <--> E(P, R) + Q;  E(A, P) + R ⇌ E(A, P, R)
E(P, R) + A ⇌ E(A, P, R)
```
Example 3, allo_bibi/1023 (MWC n = 2). Growing A alone gives 0 per decade, B gives −2. The trap is inactive-state E(A, B), which leaves only by releasing B.
```
E + A <--> E(A) [EqualAI];  E + P ⇌ E(P) [OnlyA];  E(A) + B ⇌ E(A, B) [EqualAI]
E(P) + Q ⇌ E(P, Q) [EqualAI];  E(A, B) <--> E(P, Q) [OnlyA]
```
**"No unique limit" laws (census, 1,900).** In 1,612 of them one substrate alone drives v to 0 on the core law, and growing a second substrate relieves it. In the other 288 every substrate alone saturates at a finite rate, and the limit depends only on their ratio.

Logs: `check_traps.log`, `check_drives.log`, `sum_census.log`.

### Task 3: does the rate ever exceed kcat?

Yes. In 278 of 8,992 finite-kcat laws (3.1%), the highest rate found exceeds 1.01 × kcat:

| population | finite-kcat laws | exceed kcat | share |
|---|---|---|---|
| core sets | 4,766 | 162 | 3.4% |
| expand_terter | 1,448 | 66 | 4.6% |
| init_terter sample | 2,778 | 50 | 1.8% |

- Within the core sets, rates run up to 16% (allo_uniuni_inh) and 12.8% (allo_unibi).
- Ratios reach 3,706. In 271 of the 278 the overshoot appears with modifiers held at zero.
- I re-checked at 5 fresh parameter points per law. 250 of the 278 overshoot again. The oracle confirms the rate at the maximum to 7.3e-16. With all substrates near 1e12 the rate returns to kcat (median ratio 1).
- These counts are lower bounds: the search used 3 parameter points per law and a 1e-4..1e4 box.

Which mechanisms overshoot:
- **Plain laws:** only when a steady-state step touches an abortive complex on a non-turnover cycle, 130 of 2,539 such laws. The other 2,566 plain laws never overshoot, including 542 with steady-state random binding order and no abortive complex.
- **MWC, n = 1:** 30 of 678, all with such steady-state abortive cycles; 0 of 876 without.
- **MWC, n = 2:** 115 of 1,915 where the substrate binds the two conformations differently, by tags or by a steady-state step. The remaining 3 are hand-built P5a cases.
- **MWC, n = 4:** 0 of 5.

**Steady-state bi-bi.** I built the textbook random-order steady-state bi-bi, with every step in its own group, and tried 300 parameter points:
- At fixed B, the rate as a function of A peaks above its A → ∞ value in 91 of 300. This is the classic non-monotonic behaviour at fixed co-substrate, and Denis's statement is right in that sense.
- But the rate over all A and B never exceeds kcat (0 of 300). When every substrate saturates, the enzyme sits in the ternary complex.
- Among 112 enumerated random-order laws with no steady-state abortive cycle, none overshoots at 30 points each with a wider search. The only hit was an adversarial law, adv_pkg/12.

Examples:
- expand_bibi/388: kcat 0.084, highest rate 1.34 at A = 0.45, B = 1e4. Steps: `E + A <--> E(A); E + P <--> E(P); E(A) + B <--> E(P) + Q; E(A) + P <--> E(A, P); E(P) + A <--> E(A, P)`
- allo_uniuni_inh/51 (n = 2): kcat 4.6e-9, highest rate 8.3e-6 at S = 0.035. Steps: `E + P ⇌ E(P) [OnlyA]; E + S <--> E(S) [NonequalAI]; E(S) <--> E(P) [OnlyA]`
- p1_enum_inh/1001: random-order steady-state bi-bi with abortive E(A, Q) and E(B, P). kcat 0.0146, highest rate 4.58 at A = 1e4, B = 815.

Logs: `sum_all.log`, `verify.log`, `verify2.log`, `classify_over.log`, `classify_allo.log`, `ferdinand.log`.

### Task 4: the reference-rate gauge

No failures. I took 3 parameter points and 50 concentration vectors per law, with 10 of the vectors at zero products. Across all 9,598 analysed laws and the 96,218 census laws there was no NaN or Inf, and max |v| was finite and positive at every parameter point.

There is one caveat about the reference conditions. Three hand-built laws give v = 0, or 0/0, when products and modifiers are both zero:
- **specs/24:** R is an essential activator, so v = 0 without it. `catalytic_core` also labels `E + R ⇌ E(R)` as a dead-end branch.
- **adv/14:** the inhibitor I is needed to recycle F.
- **fuzz/262:** v is 0/0 at zero products and zero I.

So a reference set needs at least one point with products or modifiers present. Logs: `check_p0m0.log`, `inspect_core.log`, `show_two.log`.

### Not checked
- Census parts 5–11 of init_terter.
- Exact (symbolic) limits for the allosteric laws; I only ran numeric growth scans on them.
- Any law built by the package itself; everything here is the slice's laws.
- A complete search for the highest rate: 3 parameter points per law within a 1e-4..1e4 box, plus the re-checks above.
- The trap populations for MWC laws come from my reimplementation of the two-conformation oracle (`mwc_pops`, written from `v_flip`). I validated only the rate it gives, not the populations themselves.

Scripts are in `K/scripts/` (`klib.jl`, `analyze.jl`, `summarize.jl`, `census_terter.jl`, `verify.jl`). Per-law rows are in `K/data/rows_*.jls` and `K/data/census_terter.jls`.