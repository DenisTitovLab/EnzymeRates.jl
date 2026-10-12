Denis, none of the four claims I tested is refuted, but two corrections change how the decision reads. I used three independent routes: an oracle I wrote from scratch (BigFloat mass-action steady state, no slice or package code), the package's own laws with their constants passed to that oracle, and the slice.

- **Face class is larger than reported.** The probe's numbers understate how far the face value can move.
- **"Wrong" is too strong for traps.** On the full face, the law is the exact RE limit for one ordering of the RE speeds.
- **The existing screens are missing.** The probe never mentions the package's no-Vmax and condition-C screens, which are the "no saturation" checks your question is about.

## Verdicts

### 1. RE turnover beside an SS route: the slice gives a finite, wrong law and no law test catches it — CONFIRMED
- **Oracle:** v(Λ=1e40) / v(Λ=1e20) = 1e20 on 5 of 5 random constant sets, for both cases below. There is no finite RE limit.
  - `[S → P]`: `E + S ⇌ E(S)`, `E(S) ⇌ E(P)`, `E + P ⇌ E(P)`, `E(S) <--> F(P)`, `F(P) <--> E + P`
  - `mut/expand_bibi/453/flip2`: `E + A ⇌ E(A)`, `E + P ⇌ E(P)`, `(E + Q ⇌ E(Q), E(P) + Q ⇌ E(P, Q))`, `E(A) + B ⇌ E(P, Q)`, `E(Q) + P <--> E(P, Q)`
- **Slice:** `check_boundary` and `admit` pass and `derive` succeeds, with Q ≢ 0. At random points the rates are 0.013 / 0.023 / 0.014 and −0.019 / −0.40 / 0.16.
- **A fit would not flag it.** 453/flip2 prints as an ordinary saturating law: `v = E_total * (A * B - P * Q / Keq) * k_EPQ_to_EQ_P * Keq / (K_EP_to_E_P * K_EQ_to_E_Q) / (1 + A / K_EA_to_E_A + P / K_EP_to_E_P + Q / K_EQ_to_E_Q + A * B * Keq / (K_EP_to_E_P * K_EQ_to_E_Q))`.
- **Package:** it rejects 453/flip2 with "the rapid-equilibrium steps close a catalytic cycle … no finite rate".
- **All-RE case:** for `E + S ⇌ E(S)`, `E(S) ⇌ E(P)`, `E(P) ⇌ E + P`, the slice prints `* () /`, so v ≡ 0, but the true rate is infinite.
- **Counts reproduce the probe's 1,654:** N≡0 in 959 to_re, 56 flip2 and 636 all_re records (635 mutants plus the skeptic's case), plus the 3 flip2 records with a finite law.
- Logs: `logs/v1_oracle.log`, `logs/v4_slice.log`, `logs/v3_pkg.log`.

### 2. Face split/trap: emitted by the enumerator; the law is finite but is not the face RE limit — PARTIAL (the facts hold; the reading needs correcting)
**Confirmed:**
- **Emission:** all 34 flagged children come out of the public `expand_mechanisms`, not only `_expand_re_to_ss`.
  - `_degenerate` is false for all 34.
  - None of their 15 plain bi-bi parents is an `init_mechanisms` seed.
- **Gap on every plain flagged child.** I ran the oracle on the package's own constants; the package law matches it at positive points to about 1e-15. All 30 plain flagged children show a gap between the law and the exact-face RE limit on the flagged face. On the other face the gap is ≤3e-16.
- **The probe's face check agrees with my oracle on plain enumerator output:**
  - Level-1 children: 20 flagged, 20 gaps, out of 5,353. No misses, no false alarms.
  - Level-2 children of 800 random level-1 children: 31 flagged, 31 gaps, out of 3,798.
  - The rate grows with depth: 0.37% at level 1, 0.82% at level 2.

**Corrections:**
- **Magnitude.** The probe drew RE speed factors from 10^[−1,1]; I used 10^[−3,3]. Over the 30 plain flagged children × 2 draws, max |face − law| / |law| has median 2.3%, p90 47% and max 314%.
  - Hand LDH-like case: up to 83%.
  - Spec 40 (active conformation as a plain network): up to 61%.
  - Spec 14: up to 56%.
- **Interpretation.** I enumerated every extreme speed ordering (each RE step at 1e±6) over 3 draws.
  - **Traps, full face:** the law equals the exact-face RE limit of one ordering (the trap-feeding step slowest) to 1e−12–1e−15. This holds for the hand LDH-like case, spec 40 active, spec 14 and child 9.
  - **The split (child 1):** the law lies strictly inside the range; the nearest ordering is 3e−4 to 2.5e−2 away.
  - **Spec 40 at NAD = 0 alone:** also inside the range, 5e−5 to 6e−4 away.
- **The LDH ODE skip is not evidence.** The skip is not specific to LDH: all 11 allosteric specs set `run_ode_test=false`, because `_testhelper_build_ode_rhs` takes only `m::EnzymeMechanism`. Spec 14 does run the ODE check, but at concentrations 0.1–10, so it never samples the face.

**Examples:**
- **Split, child 1:** `(E + A <--> E(A), E(B) + A <--> E(P, Q))`, `(E + B <--> E(B), E(A) + B <--> E(P, Q))`, `E + Q ⇌ E(Q)`, `E(B) + Q ⇌ E(B, Q)`, `E(Q) + B ⇌ E(B, Q)`, `E(Q) + P <--> E(P, Q)`. At P = Q = 0 the law gives 0.0149124; the exact-face RE limit ranges 0.0148788–0.0165063.
- **Trap, hand LDH-like:** `E + A ⇌ E(A)`, `E(A) + B <--> E(A, B)`, `E(A, B) <--> E(P, Q)`, `E(P) + Q ⇌ E(P, Q)`, `E(Q) + P ⇌ E(P, Q)`, `E + P <--> E(P)`, `E + Q ⇌ E(Q)`. Law 0.563988; exact face 0.0964–0.563984.
- **Trap, spec 40 active network (13 steps):** law 0.690457; exact face 0.2715–0.690432.

Logs: `logs/v2_pkg.log`, `logs/v3_pkg.log`, `logs/v7_orderings.log`, `logs/v8_recall.log`, `logs/v5_census.log`.

### 3. N ≡ 0: 1,027 mechanisms the package accepts and prints as v = 0 — CONFIRMED
- **Recount:** exactly 1,027 of 2,119 del_cat mutants have Q ≡ 0, matching the probe's 618 + 409.
  - Enumerated parents: 0 of 4,612.
  - Their flip and split children: 0 of 12,967 (21 more hit the term cap).
- **Package:** for `mut/init_bibi/26/del_cat` (`E + B ⇌ E(B)`, `E + P ⇌ E(P)`, `E(B) + A ⇌ E(A, B)`, `E(A, B) <--> E(P, Q)`) it prints `v = E_total * (0) / (...)`. `_degenerate` is false.
- **Not in the probe:** the fitter turns these into a FitFailure after a wasted fit.
  - Today's `loss!` is a flat 100.0 at 5 of 5 random points in relative mode, and 200.0 in absolute mode.
  - `rescale_parameter_values` then throws "_kcat_forward: no kcat components…", which identify records as a FitFailure.

Logs: `logs/v4_slice.log`, `logs/v3_pkg.log`.

### 4. (c) D = 0 on a face: the package gives NaN, the limit depends on direction, and the enumerator never makes one — CONFIRMED
- **Package:** NaN on the face for the hand `ss_dead_end_at_face` case and 4 of the probe's records (p1_enum_inh/196 and /710, init_bibi_pp/224 and /198).
- **Direction dependence:** limits approached from three directions are 0.00173 / 0.0650 / 0.0644 for init_bibi_pp/224, and 0.218 / 0.0078 / 0.224 for the hand case.
- **Exact face:** the RE limit is unique and equals 0 in all 5. In init_bibi_pp/224/del_cat, E(P) has no exit at Q = 0, so all the enzyme ends up there. Its steps are `E + A <--> E(A)`, `E + B <--> E(B)`, `E + Q <--> E(Q)`, `E(A) + B <--> E(P, Q)`, `E(B) + Q <--> E(B, Q)`, `E(B) + A <--> E(P, Q)`, `E(P) + Q <--> E(P, Q)`, `E(Q) + B <--> E(B, Q)`, `E(Q) + P ⇌ E(P, Q)`.
- **Recount:** 21 del_cat mutants (11 forward face, 10 reverse) plus 6 bottomless. None among enumerated parents or children.
- **Not in the probe:** `loss!` is 100.0 on initial-rate data (P = Q = 0) at 5 of 5 points, but 0.013–0.19 at the same points with P = Q = 0.01.

## What the probe missed that changes the conclusion
1. **The existing degeneracy screens.** The probe never discusses `_has_vmax` and `_chemistry_equilibrates_both_sides`, which `_degenerate` uses in `_base_tier` and which `_seed_candidate_screen` also runs. Its law tests cannot see these cases.
   - **Example:** `E + S <--> E(S)`, `E(S) ⇌ E(P)`, `E + P ⇌ E(P)`. The package derives `v = E_total * (k_E_S_to_ES * S - k_ES_to_E_S * P / (K_EP_to_E_P * K_ES_to_EP)) / (1 + …)`, which never saturates: at P = 0, v = 1, 1e2, 1e4, 1e6 as S goes from 1 to 1e6.
   - `_degenerate` is true. The law passes (b) and (c) and has no RE turnover. By my reading it also passes the face check, but I did not run that check on it.
   - **Today's effect:** 7 of the 264 bi-bi ping-pong mechanisms from `init_mechanisms` are degenerate (all condition C) and are never fitted; 0 of 239 bi-bi and 0 of 7 uni-bi are.
   - **"Huge loss" is not guaranteed** for a non-saturating law, because it fits sub-saturating data well.
2. **The "let the fitter sort it out" premise splits by class.**
   - **(b) and (c):** it holds, because both already score a loss of 100 per point. Neither class ever comes from the enumerator, so their guards only matter for DSL input.
   - **RE turnover beside an SS route:** it fails, because the law looks normal. This class comes only from DSL input and must stay a structural check; the package already has one.
   - **The face class:** it fails too, and the enumerator does produce it. Because the law is a valid RE limit for traps, rejecting the class is a modelling choice rather than a correctness fix. It would also reject spec 40, one of the three "LDH i-state" fixture mechanisms, and spec 14.
3. **One-product faces carry the gap too,** which matters for product-inhibition data:
   - hand LDH-like at Q = 0 alone: up to 80%
   - spec 40 active at NAD = 0 alone: up to 60%
   - child 1 at Q = 0 alone: up to 9.5%

   No plain enumerated mechanism had a gap only on a one-reactant face (0 of 1,792 parents, 5,353 level-1 children and 2,178 level-2 children). So checking full faces finds every case, but the effect on data is larger than the full-face numbers suggest.
4. **A reference-rate gauge needs (b) and (c) as well.** v ≡ 0 gives a zero gauge, and a gauge point on a (c) face gives NaN.

## Not checked
- **Allosteric mechanisms with my oracle:** I only built plain networks, and spec 40 only as its active conformation. Its n = 4 law is finite on the face (0.0155, 0.0078 and 4.21 at three draws).
- **RE turnover in enumerator output:** I argue it cannot occur, since the seed screen rejects it, flips only remove RE steps and splits keep flags. I did not compute it.
- **Other probe claims:**
  - the 6 removable bottomless cases
  - the 6 false alarms
  - the 2 non-unique fuzz cases
  - timings
  - bare-vs-reduced agreement
  - the claimed equivalence of (b) with a steady-state turnover cycle
  - `admit`'s 834 rejections
- **Other children and faces:** children of the dead-end and tagging moves, ter-ter children, and faces with modifiers or regulators at zero.
- **Sampling:** each verdict rests on 1–5 parameter draws and 3–20 speed patterns per mechanism.

All Julia runs went through the lock with a hard timeout, about 20 minutes in total. Nothing changed in the package checkout (`git status` is clean) or in `slice/src`.

Files are in `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/guard/verify/`:
- `voracle.jl`
- `v1_cases.jl`, `v1_oracle.jl`
- `vpkglib.jl`, `v2_pkg.jl`, `v3_pkg.jl`
- `v4_slice.jl`
- `v5_parts.jl`, `v5_census.jl`
- `v6_degenerate.jl`, `v7_orderings.jl`, `v8_recall.jl`
- `logs/v1_oracle.log` … `logs/v8_recall.log`