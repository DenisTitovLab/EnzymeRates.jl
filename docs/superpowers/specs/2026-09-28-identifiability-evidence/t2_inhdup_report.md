# Substrates and products as dead-end inhibitors (Denis's point 2)

Track t2. All claims carry PROVED / EXACT / SAMPLED / CONJECTURE labels. Numbers come from engine A (fast GF(p) rank, symbolic laws) unless stated otherwise. Engine B agrees on every mechanism I cross-checked: 22 hand examples and 40 regenerated children.

## 0. Answer in plain words

A substrate or product declared as a competitive inhibitor is modelled as a copy L* that binds with L's own concentration. Such a copy is non-identifiable only when it forms no new enzyme complex: every complex X·L* it creates already exists in the mechanism with the same bound composition. The existing complex is the productive XL, an init abortive complex, or a complex another copy built.

Such a copy is a second, dead orientation of an existing complex ("non-productive binding"), and steady-state data see only the total. Denis's uni-uni example is the textbook case: K_S and K_S,inh enter only as K_S + K_S,inh, with kcat rescaled, and the child's rate laws are exactly the parent's.

Copies that create at least one new complex never add a phantom (EXACT, 0 of 26,292 events). Examples are E·B in ordered bi-bi, abortive EQ·A, and F·A in ping-pong. Foreign inhibitors can never add a phantom (PROVED).

Among "no-new-complex" children, 29% are phantoms. The rest are identifiable only because a shared kinetic group pins the productive complex. For those in the proof's scope the child's family lies inside a kinetic-group split of the parent (PROVED).

## 1. Definitions

- **Reactant copy L\***: a `CompetitiveInhibitor` whose name equals a substrate or product L. It binds with concentration [L].
  - **Sites** I are the forms X it binds. All site steps X → X·L\* are RE and share one new group with association constant K\*, so [X·L\*] = K\*·L·[X].
  - **Mirror steps**: every parent step between two sites is copied between their copy forms. The copy keeps the parent step's group and type.
- **Composition** of a form: the multiset of bound metabolite names (copy flag ignored), plus conformation and covalent residual.
  - Site X is a **twin site** if composition(X)+L equals the composition of a form already in the parent.
  - A child is **all-twin** if every site is a twin: it creates no new complex.
- **Constants** follow engine A. An RE group has association-type K with [to]Π[released] = K[from]Π[consumed]. An SS group has kf and kr in the written direction.
- **Phantoms** = fitted − rank. **dPh** = phantoms(child) − phantoms(parent). A child *creates a phantom* when dPh ≥ 1.
- **Productive complex** of site X: the parent form XL reached from X by a non-copy step binding L. C = {XL : X ∈ I}.

## 2. Part A: what the package does

`_expand_add_dead_end_regulator_native` works per declared, not-yet-bound CompetitiveInhibitor, and per pattern: non-empty competing substrates × non-empty competing products × any existing inhibitors.

- **Sites** are the sources of binding steps of competing metabolites (copy steps included), filtered two ways:
  - A site must not bind all substrates or all products. Copies count as bound (form "EAinh" has bound set {A}).
  - A site must carry no competing metabolite.
- **Deduplication** is by the sorted site list.
- **When L competes with itself**, the sites are every form where L binds productively. The move then builds non-productive binding by construction.
- **In `seed_mechanisms`**, required competitive inhibitors are added by this move, so phantom children become seeds.

init's `_expand_substrate_product_dead_ends` adds f+m only when that name is not a catalytic form (`de_name in cat_forms && continue`), so it always creates a new composition. It adds no phantom at level 0 (EXACT):

- In R4 level 0, 51 of 62 seeds carry abortive complexes.
- The only phantom seeds are R4_00061 and R4_00062, the ping-pong seeds with an RE second chemistry step.

Those abortive complexes are what later copies duplicate. Example: R6_00075 places B\* at {EA, EQ}, duplicating the productive EAB and the abortive EBQ. Both steps are in group 3, and the child has dPh = +1.

**Python port** (`deadend_port.py`). It reproduces the move exactly (EXACT):

- Raw child counts: R2 1/2, R3 2/6, R5 350/2,153, R6 1,400/35,444. R6 level 2 has 25,344 distinct children, which equals the index.
- Every exported dead-end child is regenerated from its recorded parent.

This lets me cover all 25,344 R6 level-2 dead-end children, not only the 3,718 exported ones.

## 3. Part B: theory

**Lemma 1 (dwell gauge, PROVED).** Let μ assign a number to every form, and require μ to be uniform on each kinetic group:

- RE group: μ_from − μ_to is the same for all its steps.
- SS group: μ_from is the same for all steps, and so is μ_to.

Transform the constants: log K_g += μ_a − μ_b (RE); log kf_g += μ_a and log kr_g += μ_b (SS). Then:

1. Every step flux is unchanged. The forward rate kf·e^{μ_a}·e^{−μ_a}[a] is unchanged, and so is the reverse. The RE relations hold with the new K.
2. Around every cycle the μ-changes telescope to 0, so all Haldane and Wegscheider relations hold.
3. The steady state maps [F] → e^{−μ_F}[F]; it is unique up to scale because the form graph is connected.

Hence v/E_t is multiplied by Σ[F]/Σe^{−μ_F}[F]. The law is invariant iff Σ_F(1 − e^{−μ_F})[F] ≡ 0 in the concentrations. Infinitesimally the condition is Σ μ_F[F] ≡ 0, and every such μ ≠ 0 is a phantom direction.

**Theorem 1 (foreign inhibitor, PROVED).** A copy of a ligand I whose concentration occurs nowhere else has dPh ≤ 0.

Proof. At I = 0 every copy form has weight 0 and every mirror carries no flux, so v_child(θ, K\*; x, I = 0) = v_parent(θ; x). A child null direction (δθ, δK\*) therefore needs δθ ∈ null(parent). If δθ = 0, then δK\* ≠ 0 is impossible, because ∂v/∂K\* ≠ 0 at I > 0. So null(child) injects into null(parent). ∎

This agrees with the data: R2 (3/3) and R5 (2,503/2,503) have dPh ∈ {0, −1}. The argument fails for reactant copies because [L] cannot be set to 0.

**Theorem 2 (non-productive binding, PROVED).** Assume:

- (a) every site X has an RE parent step X → XL binding L, with C ∩ I = ∅;
- (b) μ = c_g on the XL of L-binding group g, s on the copy forms and 0 elsewhere satisfies the group rows for generic constants. Combinatorially: each L-binding group has only site-sourced steps; every other group touching C is uniform with respect to C's boundary; no SS mirror step exists.

Then the child's family equals the parent's, and dPh ≥ 1.

Proof. Child → parent: set ρ_g = K_g/(K_g + K\*) and apply Lemma 1 with μ = log ρ_g on XL and 0 elsewhere, then drop the copy forms. This gives:

- K'_g = K_g + K\*;
- constants of steps leaving C multiplied by ρ_g (RE: K·ρ; SS from C: kf·ρ);
- constants of steps entering C: RE K/ρ, SS kr·ρ.

Then [XL]' = (K_g + K\*)L[X], which is the child's [XL] + [X·L\*]. Every other form and every SS flux is unchanged, because copy forms touch no SS step. So v/E_t is equal.

Parent → child: pick any K\* ∈ (0, min_g K'_g), set K_g = K'_g − K\*, and invert the scalings. ∎

Denis's example (U1) is: E+S⇌ES (RE, K_S), ES→EP (SS, k), E+P⇌EP (RE, K_P), plus E+S⇌E·S\* (K\*).

- Rate law: v/E_t = k·K_S(S − P/Keq)/(1 + (K_S + K\*)S + K_P·P).
- Parent map: K_S' = K_S + K\*, k' = k·K_S/(K_S + K\*).
- Fitted/rank: 4/3 against the parent's 3/3.

**Theorem 3 (split sibling, PROVED).** Scope:

- all sites RE-bind L through one kinetic group;
- no SS mirror;
- C ∩ I = ∅.

Let S be the parent with every group split into four classes: steps leaving C, steps entering C, SS steps inside C, and the rest. Then family(child) ⊆ family(S), and S = parent exactly when Theorem 2 applies. Proof: in S the Theorem-2 map is admissible. ∎

**Observations (EXACT).**

- (i) A child with at least one non-twin site never creates a phantom: 0 of 26,292 R6 events, plus all of R2, R3 and R5.
- (ii) Every phantom-creating child keeps its parent's monomial support (all 7,631 exported children).
- (iii) The "within-segment dwell" test (WSD) predicts dPh exactly on all 31,088 R6 events whose parent has only RE binding steps. WSD uses group rows plus one row per (RE segment, ligand content) class, with Cha weights at random consistent parameters.

**Hand examples** (fitted/rank; engine A = engine B on all 22 checked). A\* denotes a copy of A.

| id | mechanism | parent | child | why |
|---|---|---|---|---|
| U1 | uni-uni RE seed + E·S\* | 3/3 | 4/3 | twin ES, Thm 2 |
| U2 | + E·P\* | 3/3 | 4/3 | product, same |
| U4 | S binding SS + E·S\* | 4/3 | 5/3 | cross-segment dwell; family = parent (fits < 1e-39) |
| U6 | merged E+S⇌X⇌E+P (SS) + E·S\* | 3/3 | 4/3 | saturated family; equal (fits < 1e-37) |
| B1 | ordered RE (E+A⇌EA, EA+B⇌EAB, EAB→EPQ SS, EQ+P⇌EPQ, E+Q⇌EQ) + E·A\* | 5/5 | 6/5 | twin EA |
| B2/B3/B4 | + EA·B\* / E·Q\* / EQ·P\* | 5/5 | 6/5 | twins EAB / EQ / EPQ |
| B5 | + E·B\* | 5/5 | 6/6 | new complex (classic dead-end EB) |
| B6 | + EQ·A\* | 5/5 | 6/6 | new complex (abortive EQA) |
| B7 | + {E·A\*, EQ·A\*}, RE mirror | 5/5 | 6/6 | one new complex |
| B10 | ordered, A binding SS, + E·A\* | 6/6 | 7/6 | twin via SS; family wider than B9 (certificate) |
| B12→B13 | B0 + abortive EAQ (K_A, K_Q shared) + E·A\* | 5/5 | 6/6 | tie-only; inside its 6/6 split sibling |
| B14 | B12 + {E·A\*, EQ·A\*} + RE mirror | 5/5 | 6/5 | gauge admissible |
| R1 | random RE, shared K_A, K_B + E·A\* | 5/5 | 6/6 | tie-only; R1 ⊂ R3 (explicit map, proper) |
| R2 | + {E·A\*, EB·A\*} + mirror in the K_B group | 5/5 | 6/5 | gauge admissible |
| R4 | random with A and B groups split (R3, 6/6) + E·A\* | 6/6 | 7/6 | split frees the gauge |
| T1 | Theorell–Chance (E+A⇌EA, EA+B→EQ+P SS, E+Q⇌EQ) + E·A\* | 3/3 | 4/3 | twin |
| T2 | TC + EA·B\* | 3/3 | 4/4 | new complex EAB |
| T5 | TC all SS + E·A\* | 5/5 | 6/6 | SS twin, identifiable |
| P1 | ping-pong RE seed + F·B\* | 7/6 | 8/6 | twin FB |
| P2/P3 | + F·A\* / E·B\* | 7/6 | 8/7 | classic dead-ends, new complexes |

**Answers to the sub-questions of Part B.**

- **RE vs SS productive binding.**
  - An RE twin with an admissible gauge is always a phantom, and the child family equals the parent's.
  - An RE twin whose gauge is blocked by a shared group is identifiable and lies inside S.
  - An SS twin is a phantom in 48 of 360 R6 events. It is always a phantom in uni-uni, where the family is saturated. It can be wider than its parent (B10).
- **Copy where L does not bind catalytically.** It is a phantom only if X·L\* reproduces a complex built by another route. Example: EP\*·Q\* has the same composition as the productive EPQ; 176 such events appear once twin sites are filtered.
- **Products** behave exactly like substrates: at seed level, 31 of 350 events are phantoms for each of A, B, P and Q.
- **Mirrors.** RE mirrors keep the gauge. SS mirrors make copy forms carry flux and fall outside the proofs: 204 all-RE events have them, 116 phantom and 88 identifiable.
- **Group sharing and Wegscheider ties** decide between phantom and tie-only; compare B13/B14 and R1/R2/R4. Haldane never obstructs (Lemma 1).
- **The gauge lens**, taken literally (Theorem 2 / SD), is compatible iff every group touching XL crosses C's boundary uniformly and the mirrors are RE. It has 0 false positives and recall 78%.
- **King–Altman-level phantoms** also exist. In R6_00065 the law obeys z = (y − x)/Keq, with x = c_AB/c_A, y = c_BQ/c_Q and z = (c_PQ − c_Q c_AP/c_A)/c_A (PROVED symbolically). B\* at {EA, EQ} shifts x and y by the same K\*, so the rank stays at 6.

## 4. Part C: predicates against exact rank

**Exported children (7,631)**

- R2: 3 children, 0 phantoms.
- R5: 2,503 children (6 with dPh = −1), 0 phantoms. Every predicate is negative.
- R3: 7 children, all dPh = +1. All-twin 7; SD 5; WSD 5; exact dwell 7.
- R6: 5,118 children (436 phantom):

| predicate | TP | FP | FN |
|---|---|---|---|
| all-twin | 436 | 966 | 0 |
| SD | 343 | 0 | 93 |
| WSD | 393 | 0 | 43 |
| exact dwell | 396 | 0 | 40 |

**Exhaustive R6 regeneration (36,844 events, 26,744 distinct children).**

- dPh: +1 = 3,082 (8.4%), 0 = 33,686, −1 = 76.
- Phantom share by parent move: init 124/1,400; split 406/2,856; dead_end 2,066/26,832; re_to_ss 486/5,756.

| predicate | TP | FP | FN | TN |
|---|---|---|---|---|
| all-twin (no new complex) | 3,082 | 7,470 | 0 | 26,292 |
| SD (Theorem 2) | 2,392 | 0 | 690 | 33,762 |
| WSD | 2,858 | 0 | 224 | 33,762 |
| exact dwell | 2,886 | 0 | 196 | 33,762 |
| all sites bind L by RE | 2,922 | 6,078 | 160 | 27,684 |

**Mismatches explained.**

- WSD is exact for parents whose binding steps are all RE.
- All 224 WSD misses have re_to_ss parents: 28 are cross-segment dwell gauges, and 196 are King–Altman-level phantoms like the R6_01882 invariant above.
- The 76 events with dPh = −1 are copies whose mirrors break a parent gauge; none of them is all-twin.

## 5. Part D: families

All-twin events split into proof classes:

- **P1: SD, 2,392 phantom events.** EQUAL by Theorem 2. Fits child→parent: 45/45 below 1e-35 (1 MemoryError). Hand U1 and B1 are below 1e-39 in both directions. Rejecting these loses nothing.
- **P2: in-scope all-RE phantoms with S ≠ parent, 412 events.** Child ⊆ S (Theorem 3). S has the same parameter count in 384 of them and one more in 28, and S itself carries a phantom.
  - The child can be wider than its parent. Certificate for R6_03189 against R6_00179: on the parent family c_AB − c_A·c_B = K1·K4 > 0; the child has K1·K5 − K2·K3, of either sign.
- **P3–P5: 118 all-RE out of scope, 48 SS-twin and 112 mixed-twin phantoms.** There is no proof here.
  - B10 against B9: c_PQ − c_A·c_B/Keq = K4·K5 > 0 on B9, but K4·K5 − K2·K9·kf3/(Keq·kr1) on B10 (certificate).
  - Child→parent fits for the WSD/XSEG/NONDWELL classes: 57 of 93 below 1e-20, 35 between 1e-4 and 0.36. The nonzero values are optimizer upper bounds; two of them are confirmed by the certificates.
- **I1: 5,900 identifiable in-scope (tie-only) events.** Child ⊆ S (Theorem 3). S has +0 parameters in 4,632, +1 in 1,156 and +2 in 112.
  - Child→S fits: 26 of 30 below 1e-36.
  - The 4 nonzero ones (3e-3 to 1e-2 with 8 starts) are proven contained, so they are optimizer shortfalls. B13→S dropped from 5e-3 to below 1e-40 with 40 starts.
  - S→child fits of 0.14–0.35 show the children are proper subsets.
- **I2–I4: 178 + 312 + 1,080 identifiable events.** No proof.

What an identifiable representation of the wider phantom families (B10, R6_03189) would be: I found none in the current step language. That a Theorell–Chance-style merge gives one is CONJECTURE.

## 6. Part E: rule

**Recommended: a dead-end copy must create a new complex.**

```julia
composition(sp, extra...) = (sort!([String.(name.(bound(sp))); String.(extra)]),
                             conformation(sp), residual(sp))
existing = Set(composition(sp) for sp in values(form_sp))    # once per mechanism
# after `active` is built:
all(f -> composition(form_sp[f], reg_name) in existing, active) && continue
```

**Effect** (EXACT):

- R6: removes 10,552 of 36,844 events (28.6%) and leaves 0 phantom-creating dead-end events.
- Seed level: removes 124 phantom and 236 tie-only events; every seed × ligand keeps at least 2 children.
- R3: removes all 7 children (all phantoms).
- R2/R5: unchanged.
- A per-site variant (drop twin sites) also leaves 0 phantoms, but replaces about 21k children with about 13k new reduced-site ones.

**What it loses**: the tie-only identifiable children, and at most 278 phantom children whose family can be wider than the parent's.

**Conservative alternatives**, run only on all-twin children:

| test | catches | loses | cost (Python, per mechanism) |
|---|---|---|---|
| SD, combinatorial | 78% | nothing | 0.03 ms |
| WSD | 93% | nothing identifiable | 0.4 ms |
| modular rank | 100% | nothing identifiable | 1.2 ms |

For comparison, the package's fitted_params costs 0.5–1.7 s per mechanism.

## 7. Open questions for Denis

1. In uni-uni, a required S or P inhibitor has no identifiable placement, so it leaves zero seeds. Should the requirement be dropped with a warning, or should the move learn second-site complexes?
2. Keep or drop the tie-only children?
3. Merged forms (point 3) carry two compositions; the twin test must check both.
4. Allosteric mechanisms are untested.
5. About 1.5% of events rest on local-optimizer fits.

Scripts and data are in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/t2_inhdup/.