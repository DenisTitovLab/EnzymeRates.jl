## kcat track: verification of the probe's report

**Bottom line for the decision.** The kcat evidence supports Denis's "let them through" for enumerated mechanisms and for the enumerator's own RE→SS children. It does not support dropping checks for mechanisms where RE steps close a turnover cycle, which in practice means user-written ones.

- **Enumerated laws and RE→SS children:**
  - All 284 kcat = 0 laws run forward. This is now exact for all 284, including the 27 MWC laws the probe only sampled, and holds on 221 more laws from parts 5–11.
  - The gauge cannot fail on these laws.
  - Among 398 RE→SS toggles, none gave v ≡ 0 and none failed to saturate.
- **SS→RE toggles:**
  - 64 of 267 derive v ≡ 0. The slice's current checks accept all 64.
  - 72 of the remaining 203 never saturate.
  - A few derive finite laws that my oracle contradicts. Those laws would not produce a huge loss.

I used my own scripts in `V/`, which is `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat/verify`. About 7 minutes of Julia ran, one process at a time under the lock.

The independent route is my own BigFloat mass-action oracle, written on the slice's `Mechanism` with RE steps at speed Λ, plus my own exact degree analysis. Against the law it agreed to 1.4e-15 on 283 plain laws (`smoke.log`), and its rate at equilibrium was zero to 7e-19.

---

### Claim 1: all 284 kcat = 0 laws run forward; exact for the 257 plain laws → CONFIRMED

- **Class counts reproduced exactly** (`summarize_zero.log`):

  | set | derivable | finite | 0 | infinite | no unique limit |
  |---|---|---|---|---|---|
  | core | 4,786 | 4,766 | 6 | 13 | 1 |
  | expand_terter | 1,512 | 1,448 | 56 | 0 | 8 |
  | ter-ter sample | 3,000 | 2,778 | 222 | 0 | 0 |

  Non-derivable counts also match: core 130 admission and 1 cycle error; expand_terter 226 admission and 188 term cap.
- **Exact positivity, 284 of 284:** every coefficient of every component is positive, and the numerator and denominator both keep terms at P = 0. So v > 0 at every positive substrate point. This includes the 27 MWC laws, which I checked compositionally.
  - My oracle matches the law at P = 0 on all 257 plain laws (worst 5.6e-16).
  - All 284 have a negative gap under uniform scaling.
- **Parts 5–11, which the census skipped:** I sampled 2,800 laws (400 per part). 221 have kcat = 0 (weighted share 8.4%; 2.8–14.5% per part), and 221 of 221 are exactly v > 0 (`census_sample.log`).
- **Example, init_terter.part7/4757.** Exact gaps for (A, B, D) are (−1, 0, 0); the trap is E(A, P).
  `E + B ⇌ E(B); E + D ⇌ E(D); E + P <--> E(P); E(B) + D ⇌ E(B, D); E(D) + B ⇌ E(B, D); E(D) + R ⇌ E(D, R); E(P) + A <--> E(A, P); E(P) + R ⇌ E(P, R); E(B, D) + A <--> E(P, Q, R); E(A, P) + R ⇌ E(A, P, R); E(P, R) + A <--> E(A, P, R); E(P, R) + Q ⇌ E(P, Q, R)`
- **Example, allo_bibi/1023 (MWC, n = 2).** Exact v > 0; gaps (A 0, B −2).
  `E + A <--> E(A) [EqualAI]; E + P ⇌ E(P) [OnlyA]; E(A) + B ⇌ E(A, B) [EqualAI]; E(P) + Q ⇌ E(P, Q) [EqualAI]; E(A, B) <--> E(P, Q) [OnlyA]`
- **Scope limit.** This does not carry over to SS→RE mechanisms: 68 of the 87 kcat = 0 mutants below have v ≡ 0. For those, Denis's "cannot catalyze" reading is literally true.

### Claim 2: kcat = 0 means one substrate traps the enzyme in a complex with a product; no substrate makes v unbounded → CONFIRMED, one wording caveat

- **Exact per-substrate gaps on the core law:** one culprit substrate in 240 laws, two in 17, identical to the probe. No gap is above 0 in the 284 laws or in the 221-law parts 5–11 sample.
- **Caveat on "occasionally two".** That holds for the core law only. On the full law, the one actually fitted, 130 plain laws have one culprit, 115 have two and 12 have three. In the parts 5–11 sample, 110 of 221 have two or more.
- **Traps, from my oracle at X × 1e9.** There were 396 traps across the 257 plain laws.
  - 274 belong to core culprits; 271 of these are clean: one complex holding X plus a product or residual, whose exits either release X or need a product. Example, expand_terter/613: the trap is E(D, R), with exits `E(R) + D ⇌ E(D, R)` and `E(D, R) + P <--> E(D, P, R)`.
  - 122 belong to substrates that inhibit only the full law (leaf dead ends); 117 of these are clean, with "release X" as the only exit.
  - The 8 unclean traps all resolve to a single such complex at X × 1e12–1e15 (`trap_check.log`). Example: expand_terter/1697 at D × 1e12 holds 100% in E(D, R).
- **Side claim: kcat ignores some substrate inhibition (906 of 8,992) → CONFIRMED.** I got exactly 906 with my own test, and the per-set counts match (init_bibi 32, init_bibi_pp 38, expand_bibi 103, p1_enum_inh 124).
  - Among modifier-free finite-kcat laws, the full law's saturating rate also depends on the direction of growth in 297 more laws.
  - In 2 more it sits below kcat (p2_fuzz_dump_flip/3: 0.52 against kcat 1.19).
  - So the core kcat misdescribes the full law in 1,146 of 8,879 modifier-free finite-kcat laws (12.9%) (`summarize_zero2.log`).

### Claim 3: the reference-rate gauge was defined for every law → PARTIAL

- **True for the tested population, but it could not have failed there.**
  - Every derivable law of the 9,843 records has all coefficients positive (0 exceptions) and Q ≢ 0. So any vector of positive concentrations gives a finite, nonzero v.
  - Only points with products and modifiers both at zero can give 0 or 0/0. That happens in exactly the probe's three laws: specs/24, p2_adv_pkg/14 and p2_fuzz_dump_flip/4.
- **It fails on the population the decision is about.** I made 878 mutants by toggling the RE/SS flag of one or two whole kinetic groups in 300 plain and 150 allosteric enumerated parents, and derived them without the checks.
  - 68 of 875 derived laws are v ≡ 0, so the gauge is 0. All 68 pass the slice's `check_boundary` and `admit`; the slice lacks the package's "no finite rate" check.
  - 120 of the 807 laws with a forward rate never saturate, against 13 of 4,786 core laws. A peak-turnover gauge, the E9 option, is infinite for these. For 70 of the 77 plain ones my oracle agrees with the law, so this behaviour is real, not an artifact.
- **Split by toggle direction** (`summarize_mutants3.log`):

  | direction | derived | v ≡ 0 | never saturating |
  |---|---|---|---|
  | RE→SS (the enumerator's move) | 398 | 0 | 0 of 398 |
  | SS→RE | 267 | 64 | 72 of 203 |
  | mixed | 210 | 4 | 48 of 206 |

- **Example of v ≡ 0, expand_bibi/120** (one group toggled each way). The chemistry is RE, so the cycle through E(B, P) turns over on RE steps alone. The oracle gives −3.8e27 at RE speed 1e30, growing ×1e10 from speed 1e20.
  `E + B ⇌ E(B); E + P <--> E(P); E(B) + A ⇌ E(A, B); E(B) + P ⇌ E(B, P); E(P) + B ⇌ E(B, P); E(P) + Q ⇌ E(P, Q); E(A, B) ⇌ E(P, Q)`

### Claim 4: v exceeds kcat in 278 laws; for plain laws only via an abortive complex; the textbook SS random-order bi-bi never does → PARTIAL

- **The overshoot is real.** Re-searching 60 random flagged laws at 5 fresh points each, 58 overshoot again. My oracle confirms the rate at the argmax to 5.2e-16 (`overshoot.log`).
- **278 is a lower bound.** 7 of 300 unflagged laws in the abortive class overshoot under a wider search (constants 1e-3..1e3, box 1e-8..1e8), by up to 268× (expand_bibi/600). My class table reproduces the probe's exactly: 868/843/542/313/2,539 laws, with all 130 flagged plain laws in the last class.
- **Textbook bi-bi: confirmed.** Over 300 parameter points the best sup/kcat is 1.000. Theory guarantees it: at P = Q = 0, v = k_chem·[E(A, B)] ≤ k_chem·E_t, which equals kcat.
- **"Only via an abortive complex" is false as a rule.** When saturation closes a faster first-order branch, no abortive complex is needed. Both hand-built cases below have no abortive form, and my oracle agrees at the argmax (8e-17 and 3.5e-16).
  - Mixed ping-pong and sequential: sup/kcat 42,960 (kcat 0.0025, sup 107 at A = 1e8, B = 95).
    `E + A ⇌ E(A); E(A) + B <--> E(A, B); E(A, B) <--> E(P, Q); E(Q) + P ⇌ E(P, Q); E + Q ⇌ E(Q); E(A) <--> E(P; +A −P); E(; +A −P) + P ⇌ E(P; +A −P); E(; +A −P) + B <--> E(B; +A −P); E(B; +A −P) <--> E(Q)`
  - Isomer branch at E(A): sup/kcat 3,093.
    `E + A ⇌ E(A); E(A) + B <--> E(A, B); E(A, B) <--> E(P, Q); E(Q) + P ⇌ E(P, Q); E + Q ⇌ E(Q); E(A) <--> F(A); F(A) + B <--> F(A, B); F(A, B) <--> E(P, Q)`
  - Two adversarial core laws in the probe's own "never overshoot" classes do overshoot: p1_adv_pkg/12 by 12,300× and p2_adv_pkg/11 by 1.032×. The probe's own `verify2.log` found adv_pkg/12, which contradicts its sentence "the other 2,566 plain laws never overshoot".
  - So Denis's E9 intuition holds for branched mechanisms, not for the textbook one.
- **MWC claims hold** (`overshoot2.log`): n = 1 without an abortive cycle, 0 of 300; n = 2 without conformation-dependent binding, 2 of 300, both hand-built P5a cases.

### What the probe missed that would change its conclusion

1. **Its population cannot answer the question.** It covered only laws that pass the slice's checks, drawn from enumerated or hand-built sets. The guards Denis wants to drop mainly stop RE steps that close turnover cycles, and the probe never saw those. See Claim 3 for the counts.
2. **"Huge loss" fails for laws that are finite but wrong.** These laws are smooth and positive, so they fit data like any other law and would be reported under the wrong mechanism.
   - Of 431 plain mutants with a forward rate, 14 disagree with the oracle.
   - The 2 worst are RE turnover running beside an SS route. p1_enum_inh/742 (two groups toggled, one each way) gives a law of 3.6e-4 against an oracle of 4.9e27, growing ×1e10 with RE speed.
   - 3 of the 158 plain RE→SS children are wrong on the forward face, so the enumerator's own move can produce them. expand_bibi/696 with one group made SS gives a face law of 0.0384 against an oracle of 0.0218, and changing the RE speed pattern moves the oracle by 0.9%:
     `E + Q ⇌ E(Q); E + A <--> E(P; +A −P); E(; +A −P) + B <--> E(Q); E(; +A −P) + P <--> E(P; +A −P); E(Q) + A ⇌ E(A, Q); E(Q) + P <--> E(B, P; +A −P); E(P; +A −P) + Q ⇌ E(A, Q); E(P; +A −P) + B <--> E(B, P; +A −P)`
   - The guard probe found the same classes at a larger scale.
3. **The gauge's definition matters.** Fifty random positive vectors cannot fail on a positive-coefficient law. The E9 peak over all concentrations is a different quantity:
   - It is infinite for the 13 never-saturating core laws and for 15% of SS→RE or mixed toggles.
   - For hyperbolic laws its peak lies at infinite substrate, so a search in a finite box underestimates it.

### Not checked

- The trap contents of the MWC laws, and the MWC laws against an independent oracle; I used the slice's law there.
- Whether the 16 admission-rejected mutants, or the 356 admission-rejected records, give correct laws when derived without `admit`.
- The package itself: every result is the slice's law, and I did not run the package's checks on my mutants.
- A recount of census parts 1–4, the 300 "no unique limit" laws, and the gauge on the 96,218 census laws (covered only by the structural argument).
- My mutants are random toggles of whole groups, not the enumerator's filtered flip move. I added no deletions or dead ends, and the oracle checks used one parameter point per mutant.

**Files.** Scripts are in `V/scripts/` (`vlib.jl`, `zero_laws.jl`, `mutants.jl`, `overshoot.jl`, `overshoot2.jl`, `overshoot_lib.jl`, `trap_check.jl`, `census_sample.jl`, plus the `summarize_*.jl` scripts). Logs are in `V/logs/`, and per-law rows are in `V/data/zero_rows.jls` and `V/data/mutant_rows.jls`.