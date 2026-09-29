# Point 9: what replaces the seeds' central isomerization (bi-bi R4: A[CX]+B[N] ⇌ P[C]+Q[NX])

Track `t4_seeds`. Scratch directory: `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/t4_seeds/`. No repository file was modified and no Julia was run.

## 0. Answer in brief

- **MERGE** removes each seed's chemistry isomerization. What remains is 4 binding groups (one per metabolite), all RE, so the rate is infinite. With k of those groups SS, the fitted count is **3 + k**. The seed has 5 (sequential) or 6 (ping-pong).
- **k = 1** (the first RE→SS flip) gives 4 parameters, *fewer* than the seed. All 248 are valid and **all 248 are degenerate** (PROVED + EXACT): the flipped ligand never saturates.
- **k = 2** gives 5 parameters: the *same* as a sequential seed, one fewer than a ping-pong seed. **130 of 372 are clean, all with rank 5.** 8 of these are the seed itself (lumping theorem, PROVED); the other 122 are new. After dedup they give **114 new families**. None of them equals any R4 mechanism emitted at levels 0–2.
- **k = 3** gives 6 parameters: all 248 clean, no phantoms. Some are exact twins of today's level-1 children (lumping).
- **Random/random**: no clean k = 2 variant exists. The first clean merged variant has 6 parameters, one *more* than the seed.
- **ELIM (Theorell–Chance)**: clean only with all 3 groups SS. That gives 5 parameters (the seed count) for 8 ordered seeds. The two binding orders give one family (PROVED).
- **Today's 7 ping-pong seeds are all degenerate**: at zero products the rate does not need B. The merged ping-pong variants with two SS groups have 5 parameters and give the classic ping-pong law.
- A **cheap structural rule (V + C + I)** matches the exact degeneracy flags on **3,056 mechanisms with 0 disagreements**. It flags 93 of today's 1,819 R4 mechanisms (all ping-pong).
- **Proposed bi-bi init: 179 seeds**, all 5-parameter, clean and identifiable, against 62 today (§8). Ordered topologies gain n + m − 2 new seeds each; this was checked on uni-uni through ter-ter.

## 1. Setup and definitions

**Mechanisms.** These use the shared JSON format. Every package step binds or releases one metabolite, or is an isomerization. A group is a set of steps sharing constants. RE groups have one association constant K; SS groups have kf and kr, in the written direction. Keq and Et are known. Group labels below use the metabolite letter, with ˣ when the group has a step touching the merged central form. F = E' = the enzyme carrying the covalent residual.

**MERGE(X1, X2).** Delete the ligand-free step X1⇌X2 and its group, and re-attach every other step to one form X. In the seeds, every chemistry group holds one step. **ELIM(X)**: if X has exactly two steps Y→X and X→Z, fuse them into one step Y→Z carrying both ligands, in a new group ("TC").

**Validity (Part B)**
1. Finite rate: no RE cycle with nonzero net stoichiometry.
2. No bottomless RE segment: the package's `_bottomless_re_segment`, re-implemented. It never fired, and it flags none of the 1,819 exported R4 mechanisms, as the package requires.
3. Every SS group carries flux: some step of the group lies on a simple cycle of the RE-contracted segment graph whose net stoichiometry is nonzero. A self-loop counts through its own stoichiometry. I checked this by brute-force enumeration of the simple cycles through each edge.

**Degeneracy flags (Part C).** Computed from the exact reduced law of engine A.
- **a**: the forward initial rate (P = Q = 0) does not depend on some substrate. **a′**: the same for the reverse rate and a product.
- **b**: at P = Q = 0 the rate grows without bound in one substrate (numerator degree > denominator degree). **b′**: the reverse analog.
- **e**: no forward Vmax, meaning v(tA, tB, 0, 0) is unbounded as t → ∞. **e′**: the reverse analog.
- **f**: the forward rate at zero products stays > 0 when a substrate is set to 0 (that substrate is not needed). **f′**: the reverse analog.
- **c**: no A·B (or P·Q) monomial in the denominator.
- **d** (informational only): the forward numerator has monomials other than A·B, such as the A²B and AB² terms of steady-state random order.

**Degenerate** means any of a, a′, b, b′, e, e′, f, f′. On all 1,149 seeds and variants, a ⊆ f, while b and c always come with a, e or f (EXACT). So the flags that matter are **e** (no Vmax) and **f** (a ligand not needed).

**Relations.** Two families are EQUAL when each reproduces the other at finite parameters. That needs an explicit map, or family_fit < 1e-20 on 5 laws in both directions with best-fit |log10 θ| ≤ 8.
- *var ⊂ cl(seed)*: reached only as parameters diverge.
- *NEITHER(fit)*: some law is not fit (≥1e-3).
- *NEITHER(supp)*: neither support contains the other. Equal families have equal generic supports and rank, so containment is impossible.

**Caveat on family_fit (EXACT).**
- A worst residual below 1e-20 can be a *limit*. Example: R4_00004|MERGE|{Pˣ} → seed fits to 3e-41, but only with K5 = 10^-1702. I therefore record the best-fit parameters (`ff_params.py`).
- A residual of about 1e-2 can be an *optimizer failure*. Pairs that are exactly EQUAL by §5.1 failed on 1–2 of 5 laws with 8 starts and passed with 40.
- idkit's LM can raise MemoryError on runaway steps. I wrapped `_mp_eval` so that |log θ| > 1e4 counts as a rejected step.
- For targets with phantoms, a large |log θ| is not diagnostic, because the fit can slide along the unidentifiable direction.

## 2. Part A: the 62 seeds

EXACT: each seed has **4 binding groups**, one per metabolite, each holding all catalytic and abortive steps of that metabolite. Chemistry is one SS isomerization EAB→EPQ in the 55 sequential seeds (fitted 5 = 4 K + kf + kr − 1). In the 7 ping-pong seeds it is EA→EP' (SS) plus EB'⇌EQ (RE) (fitted 6).

| class (substrate order / product order) | seeds | abortive decorations |
|---|---|---|
| A-first / P-off-first (E+A, EA+B, EPQ→EQ+P, EQ→E+Q) | 00001–00005 | EAP+EBQ, EBQ, EAQ, –, EAP |
| A-first / Q-off-first | 00006–00010 | EAP, –, EAQ+EBP, EBP, EAQ |
| A-first / random | 00011–00017 | 7 combos of {EAP,EBQ} or {EAQ,EBP} |
| B-first / P-off-first | 00018–00022 | EBQ, EAQ+EBP, EBP, EAQ, – |
| B-first / Q-off-first | 00023–00027 | EAP+EBQ, EBQ, EBP, –, EAP |
| B-first / random | 00028–00034 | 7 combos |
| random / P-off-first | 00035–00041 | 7 combos |
| random / Q-off-first | 00042–00048 | 7 combos |
| random / random | 00049–00055 | 7 combos |
| ping-pong (A, P off, B, Q off) | 00056–00062 | EA·P, F'P·A, EQ·B, FB·Q, EQ·A, FB·P combos; 00062 none |

Abortive steps share the metabolite's catalytic group. Note that EAQ and EBP are sometimes alternative catalytic routes (for example EQ + A → EAQ → EA + Q), not dead ends.

## 3. Part B: valid SS sets

- **MERGE**: all 15 nonempty subsets are valid for all 62 seeds (930 variants, EXACT). The only invalid set is the empty one. Finite rate for any nonempty set is PROVED: every catalytic cycle has net uptake of each metabolite, so it contains a step of that metabolite's single group. The **minimal sets are the 4 single groups**; the **sets one larger are the 6 pairs**.
- **ELIM**: applies to the 20 ordered/ordered seeds, and to ping-pong half X1 in 00057 and 00062 and half X2 in 00061 and 00062. Groups: outer substrate, outer product, TC, plus any group left holding only an abortive step. Such a leftover group can never be SS, because it carries no flux (point 1). **Minimal**: the 3 singletons. **One larger**: the 3 pairs. The triple is also valid. Ping-pong with both halves eliminated: minimal {TC1}, {TC2}; one larger {TC1, TC2}.

## 4. Part C: results

### 4.1 By number of SS groups (all 62 seeds)

| base | k | fitted − seed | n | clean | phantoms | relation to seed |
|---|---|---|---|---|---|---|
| MERGE seq | 1 | −1 | 220 | 0 | 24 | 55 ⊂cl(seed); rest NEITHER |
| MERGE seq | 2 | 0 | 330 | 116 | 8 (all degenerate) | 8 EQUAL (all clean) |
| MERGE seq | 3 | +1 | 220 | 220 | 0 | 18 seed⊂cl(var) |
| MERGE seq | 4 | +2 | 55 | 55 | 0 | 6 seed⊂cl(var) |
| MERGE pp | 1 / 2 / 3 / 4 | −2 / −1 / 0 / +1 | 28 / 42 / 28 / 7 | 0 / 14 / 28 / 7 | 0 / 0 / 0 / 1 | {A,P} EQUAL to seed (degenerate) |
| ELIM seq | 1 / 2 / 3 | −2…0 / −1…+1 / 0…+2 | 60 / 60 / 20 | 0 / 0 / 20 | 8 / 0 / 0 | NEITHER, 8 ⊂cl |
| ELIM pp | 1 / 2 / 3 | −4…−3 / −3…−2 / −1 | 8 / 7 / 2 (+14 extra for 00062) | 0 / 0 / 2 (+2) | 0 | – |

All 41 variants with phantoms are degenerate except R4_00062|MERGE|{A,B,P,Q} (the textbook four-step ping-pong, 7/6).

### 4.2 Two-group MERGE variants by structural type (EXACT)

| class | {s_adj, p_adj} | {s_out, p_adj} | {s_adj, p_out} | {s_out, p_out} | same side |
|---|---|---|---|---|---|
| ord/ord (20) | 20 clean (8 = seed) | 20 clean | 20 clean | 20 degenerate (f, f′) | 40 degenerate (e / e′) |
| ord/rand + rand/ord (28) | 56 clean (8 seeds: the two are EQUAL to each other) | – | – | 28 pairs using the ordered side's outer group: degenerate | 56 degenerate (e / e′) |
| rand/rand (7) | 28 degenerate (f) | – | – | – | 14 degenerate (e / e′) |
| ping-pong (7) | {A,Q}, {B,P}: 14 clean; {A,P}, {B,Q}: 14 degenerate (f) | | | | {A,B}, {P,Q}: 14 degenerate (e / e′) |

In a random side both groups touch X, so every pair there is of the adj type.

### 4.3 Worked example: ordered A-first / P-off-first, no abortives (R4_00004)

Seed: g1 E+A⇌EA (RE), g3 EA+B⇌EAB (RE), g4 EAB→EPQ (SS), g5 EQ+P⇌EPQ (RE), g2 E+Q⇌EQ (RE). Fitted 5, rank 5. Den {1, A, AB, PQ, Q}. After MERGE: E+A⇌EA, EA+B⇌X, X⇌EQ+P, EQ⇌E+Q.

| base | SS groups | fitted | rank | denominator | flags | rule | relation to seed |
|---|---|---|---|---|---|---|---|
| MERGE | {A} | 4 | 4 | B, BPQ, BQ, PQ | e f:B f′:PQ b:A | no | NEITHER(supp) |
| MERGE | {Q} | 4 | 4 | AB, ABP, AP, P | e′ f:AB f′:P b′:Q | no | NEITHER(supp) |
| MERGE | {Bˣ} | 4 | 4 | 1, A, PQ, Q | e b:B | no | ⊂cl(seed) |
| MERGE | {Pˣ} | 4 | 4 | 1, A, AB, Q | e′ b′:P | no | ⊂cl(seed) (fits only at K5→0) |
| MERGE | {A,Q} | 5 | 5 | AB, ABP, AP, B, BPQ, BQ, P, PQ | f:B f′:P | no | NEITHER |
| MERGE | {A,Bˣ} | 5 | 5 | 1, A, B, BPQ, BQ, PQ, Q | e | no | NEITHER |
| MERGE | **{A,Pˣ}** | 5 | 5 | 1, A, AB, B, BPQ, BQ, PQ, Q | clean | ok | NEITHER (new) |
| MERGE | **{Bˣ,Q}** | 5 | 5 | 1, A, AB, ABP, AP, P, PQ, Q | clean | ok | NEITHER (new) |
| MERGE | {Pˣ,Q} | 5 | 5 | 1, A, AB, ABP, AP, P, Q | e′ | no | NEITHER |
| MERGE | {Bˣ,Pˣ} | 5 | 5 | 1, A, AB, PQ, Q | clean | ok | **EQUAL** (lumping) |
| MERGE | 3 groups (4 variants) | 6 | 6 | … | clean | ok | {A,Bˣ,Pˣ} and {Bˣ,Pˣ,Q} are lumping twins of level-1 R4_00079 / R4_00080 |
| MERGE | all four | 7 | 7 | 11 terms | clean | ok | twin of level-2 R4_00476 (exact check) |
| ELIM | {A}, {Q}, {TC} | 3 | 3 | … | degenerate | no | NEITHER / ⊂cl |
| ELIM | {A,Q}, {A,TC}, {Q,TC} | 4 | 4 | … | f / e / e′ | no | NEITHER |
| ELIM | **{A,Q,TC}** (TC) | 5 | 5 | 1, A, AB, AP, B, BQ, P, PQ, Q | clean | ok | NEITHER (new) |

Laws at zero products (engine A, confirmed by engine B):
- **{A,Pˣ}**: v₀ = V·A·B / (Kia·Kb + Kb·A + Ka·B + A·B). Here V = k_P,off, Kb = 1/K_B,assoc, Ka = V/k_A,on, Kia = k_A,off/k_A,on. This is the full steady-state ordered forward law. The reverse law has the RE-ordered form 1 + Q + PQ.
- **{Bˣ,Q}**: the mirror image, with the RE-ordered forward form (1, A, AB) and the full steady-state reverse law.
- **TC**: V = k_Q,off, Ka = V/k_A,on, Kb = V/k_TC, Kia = k_A,off/k_A,on.
- **{A,Q}**: v₀ = c₁A/(c₂A + c₃), with no B.

### 4.4 Ping-pong without abortives (R4_00062; seed 6 fitted, rank 5, degenerate)

After MERGE: E+A⇌X1 (A), X1⇌F+P (P), F+B⇌X2 (B), X2⇌E+Q (Q).

| SS groups | fitted | rank | denominator | flags | relation |
|---|---|---|---|---|---|
| single groups | 4 | 4 | … | degenerate | – |
| **{A,Q}** (PP-AQ) | 5 | 5 | A, AB, AP, B, BQ, P, PQ, Q | clean | NEITHER (new) |
| **{B,P}** (PP-BP) | 5 | 5 | same | clean | NEITHER (new); NEITHER vs PP-AQ |
| {A,P} | 5 | 5 | AB, B, BQ, PQ, Q | f:B f′:Q | **EQUAL to the seed** (5 laws, |log10 θ| ≤ 6) |
| {B,Q} | 5 | 5 | A, AB, AP, P, PQ | f:A f′:P | NEITHER |
| {A,B}, {P,Q} | 5 | 5 | … | e / e′ | NEITHER |
| 3 groups (4 variants) | 6 | 6 | full ping-pong support | clean | each lies inside R4_00429 (5/5 laws); R4_00429 is reached only in part |
| all four | 7 | **6** | full ping-pong support | clean | lumping twin of level-1 R4_00429 (both chemistries SS) |
| half-TC: E+A→F+P (TC1) SS, B SS, Q SS | 5 | 5 | A, AB, B, BQ, P, PQ, Q | clean | new |
| half-TC mirror: A SS, P SS, F+B→E+Q (TC2) SS | 5 | 5 | A, AB, AP, B, P, PQ, Q | clean | new |
| both halves TC | ≤3 | | | always e / e′ | – |

PP-AQ at zero products: v₀ = V·A·B/(Kb·A + Ka·B + A·B), the classic ping-pong form.

### 4.5 Random sides

- **random / P-off-first (R4_00041)**: clean pairs are {Aˣ,Pˣ} and {Bˣ,Pˣ}. Their numerators carry AB² or A²B terms (the usual steady-state random-order terms). {A,Q} and {B,Q} fail f; {A,B} fails e (no forward Vmax although AB terms exist); {P,Q} fails e′.
- **A-first / random products (R4_00017)**: clean pairs are {Bˣ,Pˣ} and {Bˣ,Qˣ}.
- **random / random (R4_00055)**: every pair is degenerate. The clean minimum is 3 groups, 6 parameters: {A,B,P}, {A,B,Q}, {A,P,Q}, {B,P,Q}. No R4 mechanism with 6 or 7 parameters has a support containing any of them, so they are new at 6.

### 4.6 Theorell–Chance duplicates (PROVED)

T1: E+A⇌EA (a1, r1), EA+B⇌EQ+P (t, Haldane-dependent reverse), EQ⇌E+Q (Q release r2, Q binding a2).
T2: E+B⇌EB (b1, s1), EB+A⇌EP+Q (u), EP⇌E+P (P release s2, P binding b2).

Dividing each law by its constant term gives 11 coefficient ratios. They match exactly under b1 = t, s1 = r1, u = a1, s2 = r2, b2 = t·a1·r2/(Keq·a2·r1) = T1's reverse TC rate. The map is positive and bijective, so **the two binding orders are one family**. family_fit agrees (≤2e-37 both ways). The same holds for the A-first/Q-off-first vs B-first/P-off-first pair (≤7e-38, SAMPLED). The 8 TC variants at 5 parameters therefore give 6 distinct families.

### 4.7 Reachability in today's enumeration

**EQUAL-reachability.** Equal families must share (rank, num support, den support). Among the 122 new clean two-group variants, the only R4 mechanisms sharing a signature are their own seeds (12 dead-end ordered and 8 ordered/random seeds). Those pairs fit NEITHER (≥4e-3). So **no new 5-parameter family is emitted today at any level 0–2**.

**⊂ at +1 parameter.** Plain-seed representatives were fit against every support-compatible R4 mechanism with 5 or 6 parameters (160 pairs, two-stage, best-fit parameters recorded):

| candidate | result |
|---|---|
| OrdA-P {A,Pˣ} | only limits: 5 laws ~0 at |log10 θ| 26–3735 in R4_00079 / R4_00151; others ≥4e-2 |
| OrdB-Q {Bˣ,Q} | only limits (LM stalls at 2e-12, |log10 θ| 10–12) |
| R4_00017 {Bˣ,Pˣ} | only limits (≥4e-12, |log10 θ| 9–16) |
| R4_00041 {Aˣ,Pˣ} | only limits (|log10 θ| ≥ 26) |
| TC | partly inside R4_00335 / R4_00336 (children of random seed 00049 with the A or B group SS). The exact Gröbner map is unique and rational but gives a negative parameter on 2 of 4 random TC laws (EXACT) |
| PP-AQ, PP-BP | inside no ≤6-parameter R4 mechanism (21 candidates, all ≥1e-3). PP-BP lies inside the 7-parameter R4_00429 at finite parameters (5 laws) |

## 5. Theorems used

### 5.1 Lumping theorem (PROVED)

Let M contain X1 ⇌ X2 (SS, kf, kr), where X1's only other step is Y1 + L1 ⇌ X1 (RE, dissociation constant K1) and X2's only other step is X2 ⇌ Y2 + L2 (RE, dissociation constant K2). Each of these three steps is alone in its group. Let M′ replace X1 and X2 by one form X with Y1 + L1 ⇌ X (SS, k1, k−1) and X ⇌ Y2 + L2 (SS, k2 release, k−2 binding).

*Proof.* The steady state of X gives
X = (k1·L1·Y1 + k−2·L2·Y2)/(k−1 + k2)
and the net flux through the box
J = (k1·k2·L1·Y1 − k−1·k−2·L2·Y2)/(k−1 + k2).
In M, X1 + X2 = L1·Y1/K1 + L2·Y2/K2 and J = kf·L1·Y1/K1 − kr·L2·Y2/K2. Under kf = k2, kr = k−1, K1 = (k−1 + k2)/k1, K2 = (k−1 + k2)/k−2, the box mass and J are the same functions of (Y1, Y2, L1, L2). The rest of M sees the box only through J and through enzyme conservation, so every steady state and v coincide. For cycles through the box, (1/K1)·(kf/kr)·K2 = (k1/k−1)·(k2/k−2), so the Haldane rows coincide. The inverse map k2 = kf, k−1 = kr, k1 = (kf + kr)/K1, k−2 = (kf + kr)/K2 is positive. ∎

Exact rational checks (v difference 0 at random points, with the box's Haldane relation satisfied) on three pairs: R4_00004 vs merged {Bˣ,Pˣ}, R4_00079 vs merged {A,Bˣ,Pˣ}, R4_00476 vs merged all-SS. The theorem **fails** if a box group also holds an abortive step. That is the case for 12 ordered seeds, e.g. R4_00001: fits 0.13–0.47 and 1e-3–0.13.

Consequence: a merged form whose two adjacent groups are both SS and single-step duplicates today's RE-in / SS-iso / RE-out box. The new families are exactly those where **one adjacent step of X stays RE**, i.e. the chemistry is at equilibrium with one binding. Today's representation, which always keeps X1⇌X2, cannot express that.

### 5.2 A single SS group is degenerate (PROVED)

Removing one metabolite's steps leaves every catalytic cycle a connected path, so all forms lie in one RE segment and the SS steps are self-loops. Take the direction in which M binds. Every M-bound catalytic form is RE-connected back to E through releases of the other side's ligands, so it carries such a factor and vanishes when that side is zero. What survives are M-free forms whose weights do not involve M. So v = kf·[M]·(M-independent fraction): flag b or b′. EXACT on 248/248.

### 5.3 Vmax rule, sufficiency (PROVED)

Suppose a set F of SS steps whose forward rate involves no substrate cuts every catalytic cycle. Then at zero products v ≤ Σ_{e∈F} kf_e·Et, because net turnover must cross F. So a forward Vmax exists. The converse was EXACT on all 3,056 mechanisms.

### 5.4 Invisibility (PROVED, sound direction)

If every binding step of S starts from a form that vanishes at zero products, all S-binding rates vanish there, so the zero-product rate is S-independent.

## 6. Part D: today's seeds and hypotheses

- **Seeds**: the 55 sequential seeds are clean and identifiable. All 7 ping-pong seeds are degenerate (f: B not needed, f′: Q not needed); 00061 and 00062 also have one phantom.
- **H1: PROVED for all 7.** The only SS step is EA→EP'. Every other form lies in one RE segment reached from E through Q binding, P binding or abortive bindings, and carries P or Q. At P = Q = 0 only E and EA survive, so v = kf·K_A·A/(1 + K_A·A), with no B.
- **H2**: PROVED (§5.1, with exactly your map) when the B and P groups hold only the central-complex steps (8/20 ordered seeds). REFUTED for the 12 seeds whose flipped groups share an abortive step.
- **H3: PROVED** (§4.3): a new identifiable 5-parameter family with the full SS ordered forward law. It is not reachable today except as a limit of 6-parameter mechanisms.
- **H4**: PROVED for the plain mechanism. Segment {EA, X, EQ} drains to EQ at P = 0, so B changes only the EA:X ratio of vanishing forms. With an EQ·B abortive (R4_00001, R4_00002), B affects the rate only as an inhibitor, v₀ = A/(A·B·c1 + A·c2 + c3), so B is still not needed (f). Refuted as stated, degenerate in all 5 seeds.

## 7. A cheap structural rule for degeneracy

- **V (Vmax)**: make every SS step that binds a substrate RE; no RE cycle with nonzero stoichiometry may appear. Do the same for products.
- **C (crossing)**: no chemistry node (a merged form, a set of forms joined by RE isomerizations, or an RE TC step) may be left by both an RE step releasing a substrate and an RE step releasing a product.
- **I (visibility)**: every substrate has a binding step whose substrate-free form carries no product in its segment weight (the package's `_re_segment_extras`). Every product likewise with respect to substrates.

Rule OK ⟺ not degenerate on 1,149 seeds and variants, 14 extra ping-pong ELIM variants, 1,819 R4 export mechanisms, 95 ordered MERGE variants (uni-uni to ter-ter) and 41 ELIM variants. That is **3,056 mechanisms with 0 disagreements** (EXACT). No component, and no pair of components, suffices alone. On today's R4 export the rule removes 93 mechanisms, all ping-pong: 7/7 seeds, 26/47 level-1 and 60/171 level-2 ping-pong mechanisms; 47 of the 93 have phantoms.

## 8. Part E: answer to point 9 and the seed rule

1. **Count.** After MERGE, a single RE→SS flip gives *fewer* parameters, but it is always degenerate. Two flips give the *same* count (sequential) or *fewer* (ping-pong: 5 vs 6). Random/random needs three flips, which is *more* (6). TC needs all 3 groups SS and then equals the seed count.
2. **Rule.**
   - MERGE each chemistry step.
   - Emit every inclusion-minimal whole-group SS set that is valid and passes V/C/I.
   - Also ELIM each central form that has exactly two steps, with the same filter.
   - Apply V/C/I to today's seeds as well; this drops the ping-pong seeds.
   - Drop the lumping twin (both adjacent groups SS and single-step).
   - Put in init those with count ≤ the RE seed count; send the rest (e.g. random/random) to level 1.
   - Generalization (EXACT): an ordered topology with n substrates and m products has n + m − 1 clean pairs {s_i, p_j} with i = n or j = 1. One is the lumping twin, so n + m − 2 are new (uni-bi 1, bi-bi 2, ter-ter 4). Ordered ter-ter ELIM adds 4 clean TC variants at the seed count.
3. **Proposed bi-bi seed set** (all fitted = rank = 5, non-degenerate):

| class | today's seeds kept | new MERGE (2 groups) | new TC | total |
|---|---|---|---|---|
| ord/ord | 20 | 52 (60 − 8 lumping twins) | 6 (8 − 2 TC duplicates) | 78 |
| ord/rand + rand/ord | 28 | 48 (56 − 8 EQUAL pairs) | – | 76 |
| rand/rand | 7 | 0 (minimum 6 parameters) | – | 7 |
| ping-pong | 0 (7 degenerate seeds dropped) | 14 | 4 half-TC | 18 |
| **total** | 55 | 114 | 10 | **179** |

Your "strange equations" guess in point 9 is confirmed: {A,B} and {P,Q} lack the AB / PQ terms and have no Vmax, and singles never saturate. My recommendation is not to fit them at all, since V/C/I removes them exactly and at no cost.

## 9. Side finding (point 1)

My flux criterion is stricter than the package's `_flux_carrying_groups`. Mine asks for a cycle with nonzero stoichiometry in the RE-contracted graph; the package asks for a biconnected block containing a chemistry step. 188 level-2 R4 mechanisms fail mine, and all 188 have phantoms.

## Files (in t4_seeds/)

- **Inputs and flags**: variants.jsonl (all variants, with laws and mechanisms), extra_variants.jsonl, allflags.jsonl, final_variants.jsonl (relations and rule verdicts).
- **R4 export flags**: r4_flags.jsonl.
- **Candidate sets**: rule_candidates.json, proposed_init.json, sig_matches.json.
- **Engine B cross-check**: engine_b_check.jsonl.
- **family_fit results**: ff_cache*.jsonl; ff_more_starts.jsonl (40-start reruns).
- **Scripts**: common.py (MERGE, ELIM, validity), variants.py, rule3.py (V/C/I), generalize*.py, lump_check.py, tc_map*.py, assemble.py, tables.py.