# Merge and Theorell–Chance catalog (Denis's points 3–8 and 10)

Track t3a. Everything below was rederived in this track; the earlier spec was not used as evidence.
Scripts and data: `/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/t3a_merge_catalog/` (file list in §9).

## 0. Short answers

* **MERGE and ELIM keep thermodynamic consistency** (PROVED, §2). With every step in its own group,
  MERGE removes exactly 1 fitted parameter when the merged step is RE and exactly 2 when it is SS
  (EXACT: all 634 valid MERGEs in the catalog). ELIM removes 2 when both input steps are SS and 1
  otherwise.
* **The chain lemma (§3, PROVED)** decides the uni-uni, ordered and ping-pong cases completely.
  Call X1–X2 a *chain* when the isomerization joins two forms that each have exactly one other
  step (binding into X1, release from X2). The rest of the mechanism sees the chain only through
  four numbers (a, b, c, d): flux J = a·[Y]L − b·[Z]R and occupancy [X1]+[X2] = c·[Y]L + d·[Z]R.
  So MERGE, with the flank types set by the chain's class, gives exactly the parent's family. It
  removes exactly the chain's phantoms: 2 for SSS, 1 for SSR, RSS, SRS, RRS, SRR and RRR, and 0
  for RSR. A plain MERGE that keeps the flank types is EQUAL only when the isomerization is RE or
  all three steps are SS.
* **Point 4:** fully SS uni-uni is 5/3, and E+S⇌(ES=EP)⇌E+P is 3/3 with an EQUAL family (explicit
  maps both ways). All 7 valid uni-uni assignments fall into just 3 families. The package's 4 R1
  mechanisms are one family.
* **Point 5:** fully SS ordered bi-bi is 9/7, and E⇌EA⇌(EAB=EPQ)⇌EQ⇌E is 7/7 and EQUAL.
  Theorell–Chance is 5/5 and lacks exactly A·B·P and B·P·Q. It is a limit of the merged family,
  never reached at finite parameters. Worse, the package's RE seed (5/5) has RE→SS children at
  6/5, 6/5 and 7/5 (R4_00081, R4_00082, R4_00467) whose family is **identical** to the seed's.
* **Point 6:** the six-step all-SS ping-pong is 11/6; the textbook four-step is 7/6 and EQUAL to
  it. The textbook's invisible direction trades kf3 against kf6 at fixed 1/kf3 + 1/kf6 (kf3 is
  (EA=FP)→F+P, kf6 is (FB=EQ)→E+Q), with kr3/kf3, kr6/kf6 and the four half-reaction transfer
  coefficients held fixed. So the law fixes the total forward exit time but not its split between
  the two half-reactions.
  **No mechanism in this space is both identifiable and EQUAL to the textbook** (PROVED). Its
  family splits exactly into two identifiable 6/6 families, RSSS and SSSR, plus a 5/5 boundary
  (RSSR). The sign of D_AB·D_P − D_AP·D_B decides which piece a law belongs to (10/10 laws
  confirmed).
* **Point 7:** for fully SS random bi-bi (15/15), MERGE gives 13/13 with *identical* monomials
  (18 in the numerator, 48 in the denominator), and thermodynamics holds. But the family is **not**
  EQUAL: the result is the parent's RE-isomerization limit, with a rank 2 lower. Over all 256
  random assignments with an SS isomerization: 49 plain MERGEs are invalid and 207 give −2 fitted.
  Their rank drops by 2 in 176, by 1 in 30 and by 0 in 1 (the 1 is EQUAL, PROVED). Monomials are
  unchanged in 56 of the 207. All 207 valid random assignments with an RE isomerization carry at
  least one phantom, and merging them is exact.
* **Point 8:** ELIM of the merged central complex removes 2 parameters when both flanks are SS and
  1 otherwise. It always changes the law: it removes exactly the monomials that only the central
  complex's occupancy produced (PROVED). The result is a limit of the merged family.
* **Point 10:** yes. Among identifiable parents, 39 of 80 random parents with an SS isomerization
  (and a valid merge) give an identifiable merge with 2 fewer parameters and identical monomials.
  In the ordered topology, 3 of the 4 identifiable parents have a variant with 1 fewer parameter
  and identical monomials. All of these are strict sub-families, not EQUAL.
* **Rules D(i)–(v)**, proved or refuted in §7: (iii), the RE sandwich, holds for chains (a
  bijection with the same count) and fails for forms with more than two steps. (iv), MERGE of an RE
  isomerization, is exact whenever no step touching X1 or X2 shares its group with another step.
  Such an RE isomerization is therefore always a phantom. With shared groups it fails: package
  seeds R4_00056–60 are identifiable 6/6 and their
  merges (5/5) are not EQUAL (NEITHER for R4_00057, SAMPLED).

## 1. Definitions and conventions

**Mechanism.** Forms (nodes) and steps (edges) in the shared JSON format. A step `from→to` consumes
`consumed` and releases `released`. RE steps have one constant K, with [to]·Π[released] =
K·[from]·Π[consumed]. SS steps have kf and kr. Each step is in its own group unless stated.
`fitted` = number of parameters − rank of the cycle constraints (Keq and E_t known). `rank` = generic
rank of ∂v/∂log θ (engine A, exact). `phantoms` = fitted − rank. "f/r" means fitted/rank.

**Type strings.** One letter per step, in the order below (R = RE, S = SS).

| topology | steps, in type-string order |
|---|---|
| uni-uni | s1 E+S→ES, s2 ES→EP, s3 EP→E+P |
| ordered bi-bi | s1 E+A→EA, s2 EA+B→EAB, s3 EAB→EPQ, s4 EPQ→EQ+P, s5 EQ→E+Q |
| ping-pong | s1 E+A→EA, s2 EA→FP, s3 FP→F+P, s4 F+B→FB, s5 FB→EQ, s6 EQ→E+Q |
| random bi-bi | s1 E+A→EA, s2 E+B→EB, s3 EA+B→EAB, s4 EB+A→EAB, s5 EAB→EPQ, s6 EPQ→EP+Q, s7 EPQ→EQ+P, s8 EP→E+P, s9 EQ→E+Q |

After a MERGE the merged step is dropped from the string: ordered "RSSR" is E+A (R), EA+B→X (S),
X→EQ+P (S), EQ→E+Q (R), with X = (EAB=EPQ). After an ELIM the fused step is appended: ordered TC
"RRS" is E+A (R), EQ→E+Q (R), EA+B→EQ+P (S).

**Family.** F(M) is the set of rate laws v(concentrations; Keq) over all positive parameter values
that satisfy M's cycle constraints. cl F is its closure: laws reached as limits when parameters
diverge. Relations between a parent P and a result R:

* EQUAL: F(R) = F(P).
* R ⊂ P: F(R) ⊊ F(P).
* R ⊂ cl(P), limit only: every law of R is a limit of laws of P, and none is attained at finite
  parameters.
* NEITHER: neither family contains the other.

**Evidence labels.** PROVED means the proof is here. EXACT means an exact computation over a named
finite set. SAMPLED means family fits (engine A's mpmath Levenberg–Marquardt at 40 digits, 5 laws,
8 starts), each reported with the largest |log10 θ| of the best fit. CONJECTURE means unproved.

## 2. Part A: MERGE and ELIM

**Implementation** (`transforms.py`):

* `merge(m, X1, X2)` contracts the single ligand-free step X1–X2 into one form. Every other step
  keeps its metabolites, type and group. It refuses when the merged step's group has other steps
  or when another step also joins X1 and X2. No catalog case hit either refusal; for the
  package export see "Shared groups on the merged step itself" below.
* `elim(m, X)` orients X's two steps as Y→X (L1, R1) and X→Z (L2, R2) and replaces them with one
  step Y→Z that consumes L1+L2 and releases R1+R2, in a new group. The step is SS if either input
  was SS, else RE; `fused_type` overrides this. It refuses when an input's group has other steps.

**Theorem A1 (MERGE is thermodynamically consistent). PROVED.** Let G' = G/e0 be the form graph
with e0 = X1X2 contracted (X1 ≠ X2, e0 the only step joining them).

1. *Cycles correspond one to one.* Restricting a cycle vector to E∖{e0} maps the cycle space of
   G linearly onto that of G'.
   - It is well defined: the boundary at X is the sum of the boundaries at X1 and X2, and e0's
     contributions cancel.
   - It is injective: a cycle supported on e0 alone is zero, because e0 is not a loop.
   - The dimensions match: |E|−|V|+1 = (|E|−1)−(|V|−1)+1.

   So every cycle z' of G' is the image of exactly one cycle z of G: insert e0 wherever z' enters
   X through a step of X1 and leaves through a step of X2 (or the reverse).
2. *Same stoichiometry.* e0 carries no ligand, so z and z' have the same net stoichiometry and
   the same turnover number n.
3. *Same product of equilibrium constants.* Let K0 = [X2]/[X1] at equilibrium, and set p1 = 1+K0,
   p2 = (1+K0)/K0. Multiply the equilibrium constant of every step that enters X from the X_i side
   by p_i, and divide that of every step that leaves X from the X_i side by p_i. Concretely:
   RE K' = K·p_i or K/p_i; SS kr' = kr/p_i for a step into X_i, kf' = kf/p_i for a step out of
   X_i. A pass through X that enters at X1 and leaves at X2 then picks up p1/p2 = K0, which is
   exactly e0's factor in the lifted cycle. Entering and leaving on the same side picks up 1. So
   every cycle of G' has the same product of constants as its preimage, hence equal to Keq^n.

   The map carries every thermodynamically valid parent point to a valid merged point, so the
   merged constraints are consistent (they never force Keq).

**Theorem A2 (ELIM is thermodynamically consistent). PROVED.** X has degree 2, so every cycle
through X uses s1 and s2 with the same multiplicity (flow conservation at X). Replacing (s1, s2)
by the fused step f is therefore an isomorphism of cycle spaces. Since σ_f = σ1 + σ2, the
stoichiometry is unchanged. Setting K_f = K1·K2 (for SS, kf_f/kr_f = kf1kf2/(kr1kr2)) keeps every
cycle's product. If Y = Z the fused step is a self-loop whose own cycle has n = 1. In uni-uni,
ELIM gives E+S→E+P with v = kf·(S − P/Keq), 1/1. An RE self-loop is an all-RE catalytic cycle
and is rejected.

**Parameter counts. PROVED; EXACT on the catalog.** With own groups, the constraint rows are the
cycle vectors written in per-step coordinates. Both maps above are injective on cycle spaces, so
the constraint rank is unchanged and the fitted count changes by the change in raw parameters.

| operation | inputs | raw params | Δfitted |
|---|---|---|---|
| MERGE | RE isomerization | −1 | **−1** |
| MERGE | SS isomerization | −2 | **−2** |
| ELIM (default type) | SS+SS → SS | 4→2 | **−2** |
| ELIM (default type) | RE+SS → SS | 3→2 | **−1** |
| ELIM (default type) | RE+RE → RE | 2→1 | **−1** |
| ELIM (other choice) | RE+SS → RE | 3→1 | −2 (changes RE connectivity) |
| ELIM (other choice) | SS+SS → RE | 4→1 | −3 (changes RE connectivity) |
| ELIM (other choice) | RE+RE → SS | 2→2 | 0 (gains monomials) |

The MERGE counts were checked on every valid result: MERGE with an RE isomerization gives −1 in
225/225 cases (uni-uni 3, ordered 15, random 207). With an SS isomerization it gives −2 in 225/225
(uni-uni 3, ordered 15, random 207). Ping-pong MERGE1 gives −1/−2 in 31/31 of each type, and
MERGE12 gives −2/−3/−4 as the sum in 60/60. ELIM counts match the table in every row of §4.1–4.2.

With shared groups, the constraint rank must be recomputed and MERGE (iv) can fail (§7).

**Shared groups on the merged step itself.** Deleting one step of a multi-step group removes no
parameter, so `merge` refuses this case. It never occurs in the four catalog topologies. In the
package export it occurs for 133 of 14,120 isomerization steps (R5 and R6 only, EXACT). In each of
them the same isomerization appears a second time on an inhibitor-bound form with the same constant
(for example R5_00731: EB_res⇌EQ and EBIinh_res⇌EIinhQ share an RE group). Merging both copies
together would be the natural extension; it is not analyzed here.

## 3. The chain lemma (the tool that settles the textbook topologies)

**Setting.** X1 and X2 are joined by a ligand-free step s0. X1 has exactly one other step s1, and
X2 exactly one other step s2. Each of s1, s0 and s2 is alone in its group. Orient the chain as

Y + L --s1--> X1 --s0--> X2 --s2--> Z + R    (Y = Z is allowed, as in uni-uni).

Write u = [Y]·ΠL and w = [Z]·ΠR. The chain type is the three letters t(s1)t(s0)t(s2). For SS
steps, (f_i, r_i) = (kf, kr) in the written direction. For RE steps, K1 = [X1]/u, K0 = [X2]/[X1]
and K2 = w/[X2].

**Lemma.**

(a) At steady state (RE steps in the RE limit), [X1]+[X2] = c·u + d·w and the chain flux is
J = a·u − b·w, with (a, b, c, d) given in Table L. For RRR the whole chain is one RE segment:
w = Kc·u and [X1]+[X2] = cs·u.

(b) The rate law depends on the chain's constants only through (a, b, c, d) (or (Kc, cs)).

(c) Each unmerged type maps onto a class, and each merged Y→X→Z type is a bijection onto one:

| class | unmerged types | image | merged type (bijective) |
|---|---|---|---|
| FULL | SSS, RSR, SSR, RSS, SRS | all of (0,∞)⁴ | SS |
| C | RRS | {d = 0} | RS |
| D | SRR | {c = 0} | SR |
| RSTOR | RRR | (Kc, cs) | RR |

(d) The cycle constraints see the chain only through a/b (Kc for RRR), which equals the product
of the chain's equilibrium constants.

**Table L** (Σ = f0f2 + r1f2 + r1r0):

| chain | a | b | c | d |
|---|---|---|---|---|
| SSS | f1f0f2/Σ | r1r0r2/Σ | f1(f0+r0+f2)/Σ | r2(r1+f0+r0)/Σ |
| RSR | K1f0 | r0/K2 | K1 | 1/K2 |
| SSR (s = r1+f0) | f1f0/s | r1r0/(K2s) | f1/s | (r0+s)/(K2s) |
| RSS (s = r0+f2) | K1f0f2/s | r0r2/s | K1(s+f0)/s | r2/s |
| SRS (σ = r1+f2K0) | f1f2K0/σ | r1r2/σ | f1(1+K0)/σ | r2(1+K0)/σ |
| RRS | K1K0f2 | r2 | K1(1+K0) | 0 |
| SRR | f1 | r1/(K0K2) | 0 | (1+K0)/(K0K2) |
| RRR | Kc = K1K0K2 | | cs = K1(1+K0) | |
| merged SS (s = r1+f2) | f1f2/s | r1r2/s | f1/s | r2/s |
| merged RS | K1f2 | r2 | K1 | 0 |
| merged SR | f1 | r1/K2 | 0 | 1/K2 |
| merged RR | Kc = K1K2 | | cs = K1 | |

**Proof.**

(a) X1 and X2 touch no other step, so their balance equations (or RE relations) involve only u,
w, [X1] and [X2]. Solving them gives Table L. For SSS, for example: f1u + r0[X2] = (r1+f0)[X1] and
f0[X1] + r2w = (r0+f2)[X2]. In Cha's formulation an RE s1 puts X1 in Y's segment with weight K1u,
and the same algebra holds when the SS step is internal to a segment.

(b) At steady state, every other form's balance meets the chain only through J: J leaves Y via s1
and enters Z via s2. Enzyme conservation meets it only through [X1]+[X2]. The rate v counts the
chain's turnover only through J. Two chains with equal (a, b, c, d) therefore give the same linear
system and the same law.

(c) The images follow from the formulas. Onto-ness is shown by explicit preimages, where
α = a/c and β = b/d:

* merged SS: f2 = α, r1 = β, f1 = c(α+β), r2 = d(α+β).
* SSS: pick t > max(α, β). Set f0 = r0 = t, f2 = 2tα/(t−α), r1 = 2tβ/(t−β), f1 = cΣ/(2t+f2) and
  r2 = dΣ/(r1+2t).
* RSR: K1 = c, f0 = α, r0 = β, K2 = 1/d.
* SSR: f0 = α. Pick r0 = u > β, then r1 = β(u+f0)/(u−β), f1 = c(r1+f0) and
  K2 = (r0+r1+f0)/(d(r1+f0)).
* RSS: r0 = β. Pick f0 = u > α, then f2 = α(r0+u)/(u−α), r2 = d(r0+f2) and
  K1 = c(r0+f2)/(r0+f2+f0).
* SRS: pick K0. Then f2 = α(1+K0)/K0, r1 = β(1+K0), f1 = cσ/(1+K0) and r2 = dσ/(1+K0).
* RRS: pick K0. Then K1 = c/(1+K0), f2 = a/(K0K1), r2 = b.
* SRR: pick K0. Then K2 = (1+K0)/(K0d), f1 = a, r1 = b·K0K2.
* RRR: pick K0. Then K1 = cs/(1+K0) and K2 = Kc/(K1K0).

The free picks (t, u, K0) are the fibres, that is, the phantoms.

(d) a/b equals the product of the chain's constants in every row. ∎

**Corollaries (PROVED).**

1. *Class-representative MERGE.* Replace the chain by Y→X→Z with flank types FULL→(SS,SS),
   C→(RE,SS), D→(SS,RE), RSTOR→(RE,RE), and keep every other step. The result is EQUAL to the
   parent, with the same rank. Fitted drops by the chain's phantoms: SSS 2; SSR, RSS, SRS, RRS,
   SRR and RRR 1; RSR 0.
2. *Plain MERGE* (flank types kept) is EQUAL iff the class is kept: the isomerization is RE, or
   the chain is SSS. For SSR→SR and RSS→RS it drops c or d, giving a boundary of FULL: R ⊂ cl(P).
   For RSR→RR the chain collapses into one RE segment: the result is invalid (an all-RE catalytic
   cycle) or the fast-chain limit of the parent (both merged steps driven to RE), again R ⊂ cl(P).
3. *Useless RE→SS flips.* Flipping a flank of an RSR chain to SS (giving RSS, SSR or SSS) never
   changes the family. It only adds 1 or 2 phantoms.
4. *Closed-form rank* for these topologies (EXACT on all 101 valid assignments, not proved in
   general): rank = Σ(parameters of non-chain steps) + Σ(class
   dimensions: FULL 4, C 3, D 3, RSTOR 2) − 1 (Haldane), minus 1 more when both ping-pong chains
   are FULL (the gauge of §5.3).

**Verification.**

* EXACT, `verify_chain.py`: all 227 (parent, representative) pairs (uni-uni 7; ordered 31;
  ping-pong 63 × {chain 1, chain 2, both}). Every pair has the predicted fitted drop, equal rank,
  and identical numerator and denominator monomials. The forward and backward maps give
  bit-identical v in exact rational arithmetic (6 random points per pair). Each mapped point
  satisfies the other mechanism's cycle constraints exactly.
* Negative controls fail as expected (`negctl.py`): a perturbed map breaks both identity and
  Haldane, and SSR has no preimage in the plain-merged SR.
* EXACT: the closed-form rank matches engine A on 101/101 valid assignments (`rank_formula.py`).

## 4. Part B: catalog

### 4.1 Uni-uni (8 assignments; RRR is invalid, an all-RE catalytic cycle)

"R⊂cl(P)" means a limit only (PROVED where monomials are lost, by the positivity argument in
§7(ii)). ELIM(rep) is the Theorell–Chance form of the class representative.

| parent | fitted/rank | den | chain type(s) → class | types-kept MERGE | Δden | relation | class rep | rep fitted/rank | ELIM(rep) | ELIM fitted/rank | ELIM Δden vs parent |
|---|---|---|---|---|---|---|---|---|---|---|---|
| RRS | 3/2 | 2 | RRS→C | RS 2/2 | same | EQUAL | RS | 2/2 | S | 1/1 | −{S} |
| RSR | 3/3 | 3 | RSR→FULL | RR invalid | – | invalid | SS | 3/3 | S | 1/1 | −{P,S} |
| RSS | 4/3 | 3 | RSS→FULL | RS 2/2 | −{P} | R⊂cl(P) | SS | 3/3 | S | 1/1 | −{P,S} |
| SRR | 3/2 | 2 | SRR→D | SR 2/2 | same | EQUAL | SR | 2/2 | S | 1/1 | −{P} |
| SRS | 4/3 | 3 | SRS→FULL | SS 3/3 | same | EQUAL | SS | 3/3 | S | 1/1 | −{P,S} |
| SSR | 4/3 | 3 | SSR→FULL | SR 2/2 | −{S} | R⊂cl(P) | SS | 3/3 | S | 1/1 | −{P,S} |
| SSS | 5/3 | 3 | SSS→FULL | SS 3/3 | same | EQUAL | SS | 3/3 | S | 1/1 | −{P,S} |

The 7 valid assignments give only 3 families: FULL = reversible Michaelis–Menten
V(S−P/Keq)/(1+S/Ks+P/Kp) (RSR, RSS, SSR, SRS, SSS); C, with no P term (RRS); D, with no S term
(SRR). The package's R1 mechanisms (RSR seed, RSS, SSR, SSS) are therefore one family: its 3 R1
children add nothing.

### 4.2 Ordered bi-bi (32 assignments; RRRRR invalid)

| parent | fitted/rank | den | chain type(s) → class | types-kept MERGE | Δden | relation | class rep | rep fitted/rank | ELIM(rep) | ELIM fitted/rank | ELIM Δden vs parent |
|---|---|---|---|---|---|---|---|---|---|---|---|
| RRRRS | 5/4 | 4 | RRR→RSTOR | RRRS 4/4 | same | EQUAL | RRRS | 4/4 | RSR | 3/3 | −{A*B*P} |
| RRRSR | 5/4 | 4 | RRS→C | RRSR 4/4 | same | EQUAL | RRSR | 4/4 | RRS | 3/3 | −{A*B} |
| RRRSS | 6/5 | 7 | RRS→C | RRSS 5/5 | same | EQUAL | RRSS | 5/5 | RSS | 4/4 | −{A*B*P} |
| RRSRR | 5/5 | 5 | RSR→FULL | RRRR invalid | – | invalid | RSSR | 5/5 | RRS | 3/3 | −{A*B,P*Q} |
| RRSRS | 6/6 | 8 | RSR→FULL | RRRS 4/4 | −{1,A,P*Q,Q} | R⊂cl(P) | RSSS | 6/6 | RSS | 4/4 | −{A*B*P,P*Q} |
| RRSSR | 6/5 | 5 | RSS→FULL | RRSR 4/4 | −{P*Q} | R⊂cl(P) | RSSR | 5/5 | RRS | 3/3 | −{A*B,P*Q} |
| RRSSS | 7/6 | 8 | RSS→FULL | RRSS 5/5 | −{P*Q} | R⊂cl(P) | RSSS | 6/6 | RSS | 4/4 | −{A*B*P,P*Q} |
| RSRRR | 5/4 | 4 | SRR→D | RSRR 4/4 | same | EQUAL | RSRR | 4/4 | RRS | 3/3 | −{P*Q} |
| RSRRS | 6/5 | 8 | SRR→D | RSRS 5/5 | same | EQUAL | RSRS | 5/5 | RSS | 4/4 | −{A*B*P,P*Q} |
| RSRSR | 6/5 | 5 | SRS→FULL | RSSR 5/5 | same | EQUAL | RSSR | 5/5 | RRS | 3/3 | −{A*B,P*Q} |
| RSRSS | 7/6 | 8 | SRS→FULL | RSSS 6/6 | same | EQUAL | RSSS | 6/6 | RSS | 4/4 | −{A*B*P,P*Q} |
| RSSRR | 6/5 | 5 | SSR→FULL | RSRR 4/4 | −{A*B} | R⊂cl(P) | RSSR | 5/5 | RRS | 3/3 | −{A*B,P*Q} |
| RSSRS | 7/6 | 8 | SSR→FULL | RSRS 5/5 | same | R⊂cl(P) | RSSS | 6/6 | RSS | 4/4 | −{A*B*P,P*Q} |
| RSSSR | 7/5 | 5 | SSS→FULL | RSSR 5/5 | same | EQUAL | RSSR | 5/5 | RRS | 3/3 | −{A*B,P*Q} |
| RSSSS | 8/6 | 8 | SSS→FULL | RSSS 6/6 | same | EQUAL | RSSS | 6/6 | RSS | 4/4 | −{A*B*P,P*Q} |
| SRRRR | 5/4 | 4 | RRR→RSTOR | SRRR 4/4 | same | EQUAL | SRRR | 4/4 | SRR | 3/3 | −{B*P*Q} |
| SRRRS | 6/5 | 8 | RRR→RSTOR | SRRS 5/5 | same | EQUAL | SRRS | 5/5 | SSR | 4/4 | −{A*B*P,B*P*Q} |
| SRRSR | 6/5 | 8 | RRS→C | SRSR 5/5 | same | EQUAL | SRSR | 5/5 | SRS | 4/4 | −{A*B,B*P*Q} |
| SRRSS | 7/6 | 11 | RRS→C | SRSS 6/6 | same | EQUAL | SRSS | 6/6 | SSS | 5/5 | −{A*B*P,B*P*Q} |
| SRSRR | 6/6 | 8 | RSR→FULL | SRRR 4/4 | −{1,A,A*B,Q} | R⊂cl(P) | SSSR | 6/6 | SRS | 4/4 | −{A*B,B*P*Q} |
| SRSRS | 7/7 | 11 | RSR→FULL | SRRS 5/5 | −{1,A,Q} | R⊂cl(P) | SSSS | 7/7 | SSS | 5/5 | −{A*B*P,B*P*Q} |
| SRSSR | 7/6 | 8 | RSS→FULL | SRSR 5/5 | same | R⊂cl(P) | SSSR | 6/6 | SRS | 4/4 | −{A*B,B*P*Q} |
| SRSSS | 8/7 | 11 | RSS→FULL | SRSS 6/6 | same | R⊂cl(P) | SSSS | 7/7 | SSS | 5/5 | −{A*B*P,B*P*Q} |
| SSRRR | 6/5 | 7 | SRR→D | SSRR 5/5 | same | EQUAL | SSRR | 5/5 | SRS | 4/4 | −{B*P*Q} |
| SSRRS | 7/6 | 11 | SRR→D | SSRS 6/6 | same | EQUAL | SSRS | 6/6 | SSS | 5/5 | −{A*B*P,B*P*Q} |
| SSRSR | 7/6 | 8 | SRS→FULL | SSSR 6/6 | same | EQUAL | SSSR | 6/6 | SRS | 4/4 | −{A*B,B*P*Q} |
| SSRSS | 8/7 | 11 | SRS→FULL | SSSS 7/7 | same | EQUAL | SSSS | 7/7 | SSS | 5/5 | −{A*B*P,B*P*Q} |
| SSSRR | 7/6 | 8 | SSR→FULL | SSRR 5/5 | −{A*B} | R⊂cl(P) | SSSR | 6/6 | SRS | 4/4 | −{A*B,B*P*Q} |
| SSSRS | 8/7 | 11 | SSR→FULL | SSRS 6/6 | same | R⊂cl(P) | SSSS | 7/7 | SSS | 5/5 | −{A*B*P,B*P*Q} |
| SSSSR | 8/6 | 8 | SSS→FULL | SSSR 6/6 | same | EQUAL | SSSR | 6/6 | SRS | 4/4 | −{A*B,B*P*Q} |
| SSSSS | 9/7 | 11 | SSS→FULL | SSSS 7/7 | same | EQUAL | SSSS | 7/7 | SSS | 5/5 | −{A*B*P,B*P*Q} |

Menu of the merged topology E⇌EA⇌X⇌EQ⇌E (15 valid of 16, all identifiable), by type string:

| type | f/r | den |
|---|---|---|
| SSSS | 7/7 | 11 |
| SSSR | 6/6 | 8 |
| SSRS | 6/6 | 11 |
| SRSS | 6/6 | 11 |
| RSSS | 6/6 | 8 |
| RSSR | 5/5 | 5 |
| RSRS | 5/5 | 8 |
| SRSR | 5/5 | 8 |
| SRRS | 5/5 | 8 |
| RRSS | 5/5 | 7 |
| SSRR | 5/5 | 7 |
| RRSR | 4/4 | 4 |
| RSRR | 4/4 | 4 |
| RRRS | 4/4 | 4 |
| SRRR | 4/4 | 4 |

Theorell–Chance menu (7 valid of 8, all identifiable): SSS 5/5 (9), SSR 4/4 (6), SRS 4/4 (6),
RSS 4/4 (6), RRS 3/3 (3), RSR 3/3 (3), SRR 3/3 (3). Laws are in `catalog_ordered.json`. For
example, TC SSS has den {1, A, B, P, Q, AB, AP, BQ, PQ}. The merged SSSS adds ABP and BPQ.

Identifiable ordered parents are exactly the 4 with an RSR central chain: RRSRR 5/5 (the package
seed R4_00004), RRSRS 6/6, SRSRR 6/6 and SRSRS 7/7. Every other assignment has 1 or 2 chain
phantoms.

### 4.3 Ping-pong (64 assignments; RRRRRR invalid)

Parent f/r by chain types. Chain 1 = (s1 s2 s3) is E+A→EA→FP→F+P; chain 2 = (s4 s5 s6) is
F+B→FB→EQ→E+Q.

| chain1 \\ chain2 | SSS (FULL) | RSR (FULL) | SSR (FULL) | RSS (FULL) | SRS (FULL) | RRS (C) | SRR (D) | RRR (RSTOR) |
|---|---|---|---|---|---|---|---|---|
| SSS (FULL) | 11/6 | 9/6 | 10/6 | 10/6 | 10/6 | 9/6 | 9/6 | 8/5 |
| RSR (FULL) | 9/6 | 7/6 | 8/6 | 8/6 | 8/6 | 7/6 | 7/6 | 6/5 |
| SSR (FULL) | 10/6 | 8/6 | 9/6 | 9/6 | 9/6 | 8/6 | 8/6 | 7/5 |
| RSS (FULL) | 10/6 | 8/6 | 9/6 | 9/6 | 9/6 | 8/6 | 8/6 | 7/5 |
| SRS (FULL) | 10/6 | 8/6 | 9/6 | 9/6 | 9/6 | 8/6 | 8/6 | 7/5 |
| RRS (C) | 9/6 | 7/6 | 8/6 | 8/6 | 8/6 | 7/5 | 7/5 | 6/4 |
| SRR (D) | 9/6 | 7/6 | 8/6 | 8/6 | 8/6 | 7/5 | 7/5 | 6/4 |
| RRR (RSTOR) | 8/5 | 6/5 | 7/5 | 7/5 | 7/5 | 6/4 | 6/4 | invalid |

**No six-step ping-pong assignment is identifiable** (EXACT, 63/63). Each assignment's
class-representative double MERGE is EQUAL to it (EXACT maps, §3).

Merged (textbook) menu E⇌X1⇌F⇌X2⇌E, types (s1 s3 s4 s6):

| type | f/r | den |
|---|---|---|
| SSSS | 7/6 | 8 |
| RSSS | 6/6 | 8 |
| SRSS | 6/6 | 8 |
| SSRS | 6/6 | 8 |
| SSSR | 6/6 | 8 |
| RSSR | 5/5 | 8 |
| SRRS | 5/5 | 8 |
| RSRS | 5/5 | 7 |
| SRSR | 5/5 | 7 |
| RRSS | 5/5 | 5 |
| SSRR | 5/5 | 5 |
| RRRS | 4/4 | 4 |
| RRSR | 4/4 | 4 |
| RSRR | 4/4 | 4 |
| SRRR | 4/4 | 4 |

Only SSSS is non-identifiable.

ELIM menus:

* One half fused (E+A→F+P, other half merged; type string = the other half's two flanks, then the
  fused step): SSS 5/5 (den 7), RSS 4/4 (6), SRS 4/4 (6), SSR 4/4 (4), RRS, RSR, SRR 3/3. Fusing
  F+B→E+Q instead gives the same menu.
* Both halves fused, E+A→F+P and F+B→E+Q: SS 3/3, v = (k1k2·AB − k₋1k₋2·PQ)/(k1A + k2B + k₋1P
  + k₋2Q), no constant term. RS and SR are 2/2.

ELIM of (EA=FP) removes exactly {AP}; ELIM of (FB=EQ) removes exactly {BQ}; both remove
{AB, AP, BQ, PQ}.

### 4.4 Random bi-bi (512 assignments, each step its own group; engine A symbolic on all 512)

* 463 are valid; the other 49 are all-RE catalytic cycles.
* **Isomerization s5 RE: 207 valid, 0 identifiable** (phantoms 1:80, 2:82, 3:41, 4:4). The plain
  MERGE is EXACT in all 207: −1 fitted, same rank, same monomials, and bit-identical v both ways
  under the Theorem A1 map with Haldane preserved (`verify_re_merge.py`).
* **Isomerization s5 SS: 256 valid, 81 identifiable.** The plain MERGE is invalid in 49 (all-RE
  cycle) and valid in 207:
  - Δfitted is −2 in all 207.
  - Δrank is −2 in 176, −1 in 30 and 0 in 1.
  - Numerator and denominator monomials are unchanged in 56 (39 with Δrank −2, 16 with −1, 1 with
    0).
  - Of the 80 identifiable parents with a valid merge, the merge is identifiable in all 80 and
    keeps the monomials in 39.
* Fully SS: SSSSSSSSS 15/15, 18 numerator and 48 denominator monomials. The merged SSSSSSSS is
  13/13 with identical monomials.
* The RE-seed analog with own groups, RRRRSRRRR, is 7/7 with den {1, A, B, AB, P, Q, PQ}. Its plain
  MERGE is invalid.

**Package grouped random seed R4_00055** (A, B, P and Q bindings each share one K across their two
contexts; s5 SS): 5/5. The plain MERGE is invalid (all-RE catalytic cycle). A merge with SS flanks
would have to split the shared binding groups, which raises the count to 9 fitted. So MERGE has
nothing to offer this seed.

### 4.5 Package ping-pong seeds (R4 level 0, the RE second chemistry EB_res ⇌ EQ)

| seed | parent f/r | MERGE(EB_res, EQ) f/r | den | groups at X1, X2 | relation |
|---|---|---|---|---|---|
| R4_00062 (the six-step RSRRRR) | 6/5 | 5/5 | same | exclusive | EQUAL (PROVED iv; SAMPLED <1e-37 both ways) |
| R4_00061 | 6/5 | 5/5 | same | exclusive | EQUAL (PROVED iv; SAMPLED <1e-40 both ways) |
| R4_00057 | 6/6 | 5/5 | same | Q and B binding groups shared between a step at EQ and a step at EB_res | NEITHER (SAMPLED 0.037 / 0.014, parameters moderate) |
| R4_00056, 58, 59, 60 | 6/6 | 5/5 | same | shared, as in 57 | rank drops, so not EQUAL (EXACT) |

In the package's first chemistry (the SS step EA→EP_res), the plain MERGE is invalid (all-RE
catalytic cycle) for all 7 ping-pong seeds.

### 4.6 Family fits (SAMPLED)

"P⊆R" means fitting the result to parent laws; "R⊆P" means fitting the parent to result laws.
The entry is the largest per-law residual over 5 laws, with the largest |log10 θ| of the best
fits. A residual near 0 with moderate |log10 θ| means the law is attained. A residual near 0 only
with |log10 θ| ≫ 5 means the law is approached as a limit.

| # | pair (parent vs result) | P ⊆ R: max residual (max abs log10 θ) | R ⊆ P: max residual (max abs log10 θ) | theory |
|---|---|---|---|---|
| 0 | uni SSS vs merged SS | 2.7e-39 (3.2) | 1.2e-39 (41.9) | EQUAL (chain lemma) |
| 1 | ord RRSRR vs merged RSSR | 6.2e-41 (3.2) | 4.4e-39 (2.5) | EQUAL (RE sandwich) |
| 2 | ord SSSSS vs merged SSSS | 1.4e-36 (2.1) | 1.1e-37 (1.6) | EQUAL (chain lemma) |
| 3 | ord SRSRS vs merged SSSS | 1.6e-37 (2.7) | 2.7e-36 (2.7) | EQUAL (chain lemma) |
| 5 | uni SSR vs kept SR | 6.5e-01 (6283774791092589.0) | 3.5e-36 (2163275409490.4) | R ⊂ cl(P), limit only |
| 6 | ord RSSRS vs kept RSRS | 5.8e-02 (2.5) | 4.1e-06 (12.9) | R ⊂ cl(P), limit only |
| 7 | ord RRSRS vs kept RRRS | 8.4e-01 (338.5) | 6.9e-08 (13.5) | R ⊂ cl(P), limit only (chain driven to RE) |
| 8 | ord merged SSSS vs TC SSS | 4.4e-01 (33.4) | 3.7e-12 (13.6) | R ⊂ cl(P), limit only |
| 9 | ord merged RSSR vs TC RRS | 7.4e-01 (16.4) | 8.2e-13 (238460419.1) | R ⊂ cl(P), limit only |
| 10 | uni merged SS vs ELIM S | 8.2e-01 (2.0) | 1.3e-12 (2253.6) | R ⊂ cl(P), limit only |
| 11 | pp textbook SSSS vs RSSS | 1.7e-03 (4401479281.2) | 6.8e-38 (1.9) | R ⊂ P (piece) |
| 12 | pp textbook SSSS vs SSSR | 7.2e-02 (16.5) | stopped unfinished | R ⊂ P (piece) |
| 13 | pp textbook SSSS vs SRSS | 1.4e-02 (16.6) | 5.9e-39 (606.3) | R ⊂ P (piece) |
| 14 | pp textbook SSSS vs SSRS | 1.1e-03 (15.5) | 3.5e-37 (1.8) | R ⊂ P (piece) |
| 17 | rand RRSSSRRSS vs merged RRSSRRSS | 2.0e-01 (989744335.3) | 6.0e-02 (731.8) | R ⊂ cl(P), rank −1 |
| 20 | rand SSSSRSSSS vs merged SSSSSSSS | 4.5e-02 (90041671.9) | stopped unfinished | EQUAL (iv) |
| 22 | ord RRSRS vs merged RSRS | 2.5e-01 (17.5) | 4.1e-12 (12.7) | R ⊂ cl(P), rank −1 |
| 23 | ord SRSRS vs merged SRSS | 1.5e-01 (21643597.9) | 2.7e-12 (14.2) | R ⊂ cl(P), rank −1 |
| 24 | ord RRSRR vs merged RRSR | 4.5e-01 (17.8) | 3.8e-27 (67883646124.8) | R ⊂ cl(P) (loses PQ) |
| 25 | pkg R4_00057 vs MERGE RE iso | 3.7e-02 (1.7) | 1.4e-02 (1.5) | not EQUAL (rank 6 vs 5) |
| 26 | pkg R4_00061 vs MERGE RE iso | 1.7e-40 (2.3) | 6.3e-41 (1.9) | EQUAL (iv) |
| 27 | pkg R4_00062 vs MERGE RE iso | 2.9e-40 (1.5) | 6.4e-38 (1.4) | EQUAL (iv) |

Rows 11–14 are per-law piece tests: the P ⊆ R maximum is large because some of the 5 textbook
laws lie outside that piece. The laws inside fit to ≤3e-35 (§5.3).

**Caveats.**

* For targets with 11 or more parameters the local fitter is unreliable. For SSSSRSSSS (14/13)
  versus its RE-merge (13/13), which is EQUAL by EXACT maps, it returned 0.045 (row 20). Row 17's
  R ⊆ P value (0.06) is likewise an optimizer failure: theory says R ⊂ cl(P). Random bi-bi
  relations here therefore rest on proofs and ranks, not on fits.
* The six-step versus textbook fit (job 4, target 7 or 11 parameters with a gauge direction) and
  row 12's R ⊆ P were stopped unfinished after more than 15 minutes. Both are covered by EXACT maps
  (§3) and by rows 11, 13 and 14.

## 5. Part C: the points

### 5.1 Point 4 (uni-uni)

* Fully SS E+S⇌ES⇌EP⇌E+P is **5/3** (engine A; EXACT). The chain is SSS (FULL, 2 phantoms).
* E+S⇌(ES=EP)⇌E+P with both steps SS is **3/3**. Its law is
  v = kf1·kf3·kr1·(Keq·S − P)/(Keq·kf1·kr1·S + Keq·kf3·kr1 + Keq·kr1² + kf1·kf3·P), i.e. reversible
  Michaelis–Menten.
* The two are **EQUAL**: PROVED by the chain lemma with explicit maps both ways (§3), EXACT
  (bit-identical laws), and SAMPLED (<3e-39 both ways, |log10 θ| ≤ 3.2 except one law at 41.9
  along a phantom direction of the 5-parameter parent).

Denis's statement holds: 5 → 3 fitted with an identical family. The same 3/3 family is also reached
by the RE seed RSR (E+S⇌ES RE, ES→EP SS, EP⇌E+P RE), which the package emits as R1 level 0. That
seed is EQUAL to the merged SS with the same count (the RE sandwich, §7(iii)). So the package's R1
already has an identifiable representative. Its RE→SS children RSS, SSR and SSS (4/3, 4/3, 5/3)
are exact duplicates of the seed's family.

### 5.2 Point 5 (ordered bi-bi)

* Fully SS E⇌EA⇌EAB⇌EPQ⇌EQ⇌E is **9/7**, and E⇌EA⇌(EAB=EPQ)⇌EQ⇌E (SSSS) is **7/7**. They are
  **EQUAL**, with identical 11 denominator monomials: PROVED, EXACT, and SAMPLED (<1.5e-36 both
  ways, |log10 θ| ≤ 2.1).
* Theorell–Chance E⇌EA⇌EQ⇌E (EA+B→EQ+P) is **5/5**. It lacks exactly {ABP, BPQ} and has den {1, A,
  B, P, Q, AB, AP, BQ, PQ}.
* TC ⊂ cl(merged), a limit only. PROVED: set c, d → 0 in the chain (§7(v)). SAMPLED: 3.7e-12 with
  |log10 θ| ≈ 11–14. In the other direction, merged ⊄ TC (0.44).

**Extra finding.** The package's ordered RE seed RRSRR (5/5; E+A RE, EA+B RE, EAB→EPQ SS, EPQ⇌EQ+P
RE, EQ⇌E+Q RE) is EQUAL to the merged RSSR (5/5, same count). That is, E+A RE, EA+B→X SS, X→EQ+P SS,
EQ→E+Q RE. PROVED, and SAMPLED at <5e-39. Its RE→SS children flip a flank of the RSR chain and are
**exactly the same family**:

* RRSSR and RSSRR are 6/5.
* RSSSR is 7/5.
* In the export these are R4_00081, R4_00082 (level 1) and R4_00467 (level 2).

(§6 gives the export-wide count.)

### 5.3 Point 6 (ping-pong)

**Six-step to textbook.** SSSSSS is **11/6**. Its double MERGE, the textbook E⇌(EA=FP)⇌F⇌(FB=EQ)⇌E
(SSSS), is **7/6**. They are **EQUAL** (PROVED by the chain lemma on each half, and EXACT). So the
enumerator's six-step all-SS form has five more fitted parameters than its identifiable content
(11 − 6). Merging removes four of them; the fifth is the textbook's own phantom.

**The textbook's invisible direction** (PROVED, EXACT). Each half is a merged SS chain with
(a_i, b_i, c_i, d_i). Solving the two-node cycle E⇌F gives

v = (a1a2·AB − b1b2·PQ) / D,

D = a1A + a2B + b1P + b2Q + (a2c1 + a1c2)·AB + (b1c1 + a1d1)·AP + (a2d2 + b2c2)·BQ + (b1d2 + b2d1)·PQ.

This formula matches engine A's textbook law exactly, 15/15 rational points
(`pp_effective_check.py`). With b2 set by Haldane, the coefficient map has rank 6 on 7 free
parameters (`pp_gauge.py`). Its kernel is

δ(c1, d1, c2, d2) = ε·(1, −b1/a1, −a2/a1, b2/a1).

This moves flux-proportional occupancy between the two central complexes: δ[X1] = +ε·v/a1 and
δ[X2] = −ε·v/a1, so total enzyme is unchanged. In rate constants (the merged-chain formulas give
c1/a1 = 1/kf3, c2/a2 = 1/kf6, c1b1/a1 + d1 = kr3/kf3, c2b2/a2 + d2 = kr6/kf6):

* The law depends on kf3 and kf6 only through **1/kf3 + 1/kf6**, the forward kcat relation
  1/kcat = 1/kf3 + 1/kf6. It also fixes kr3/kf3, kr6/kf6 and a1, b1, a2, b2.
* Likewise only 1/kr1 + 1/kr4 is fixed.
* So which half-reaction limits turnover is invisible.

Moving kf3 and kf6 along this direction, with compensating kf1, kr1, kf4, kr4, kr3 and kr6, left v
bit-identical at 21/21 exact rational points (`pp_gauge_check.py`).

**Is any mechanism in this space identifiable and EQUAL to the textbook? No** (PROVED).

1. Rank 6 needs two chains with class dimensions summing to 7 or 8 (§3, corollary 4). ELIM (class
   dimension 2) caps the rank at 5.
2. FULL+FULL has rank 6 but at least 7 fitted: 7/6 at best.
3. The identifiable rank-6 mechanisms are therefore one FULL chain plus one C or D chain, e.g.
   merged RSSS, SSSR, SRSS, SSRS (6/6), or the same class pairs with an RSR half left unmerged. Each has one of
   c1, d1, c2, d2 equal to 0.
4. Along the gauge line, a textbook law can be moved until d1 or c2 reaches 0 (ε > 0), or until
   c1 or d2 reaches 0 (ε < 0).
5. Hence F(textbook) = F(RSSS) ⊔ F(SSSR) ⊔ F(RSSR), sorted by the sign of

   D_AB·D_P − D_AP·D_B = a1(b1c2 − a2d1): > 0 in RSSS, < 0 in SSSR, = 0 in RSSR (5/5).

   Likewise F(textbook) = F(SSRS) ⊔ F(SRSS) ⊔ F(SRRS), sorted by D_AB·D_Q − D_BQ·D_A. Each piece is
   a proper subset of the textbook family.

SAMPLED check (`pp_pieces.py`): 10/10 random textbook laws matched the prediction.

* The two predicted candidates fitted to ≤1e-35 with |log10 θ| ≤ 3.2.
* The two non-predicted candidates stayed at ≥8e-5.

The 5-law family fits agree: the RSSS and SSSR pieces took laws {0,1,3,4} and {2}; SSRS took
{0,1,3} and SRSS took {2,4}. Each law landed in exactly one piece of each pair.

**Practical consequence.** The textbook 7/6 has no identifiable equivalent. Either count it as 6,
or emit the pair RSSS and SSSR (6/6 each), which together with the 5/5 RSSR covers exactly its
family.

### 5.4 Point 7 (random bi-bi)

**Thermodynamics:** consistent (Theorem A1, PROVED for any MERGE).

**Concentration terms:** for fully SS, identical (EXACT: 18 numerator and 48 denominator monomials).
Fitted drops 15 → 13 (EXACT).

**Family: not EQUAL.** Plain MERGE of an SS isomerization is the same thing as switching the
isomerization to RE and then applying the exact RE merge (PROVED: Theorem A1 with (iv)). So the
result is the parent's RE-isomerization limit, R ⊂ cl(P). The fully SS parent is identifiable at
15/15, so the 13/13 result is a strict 13-dimensional sub-family.

Across the 256 SS-isomerization assignments (§4.4):

* 49 plain merges are invalid.
* Of the 207 valid ones, 206 lose rank (176 by 2, 30 by 1).
* One is EQUAL: RRSSSSSRR (11/7) → RRSSSSRR (9/7).

That EQUAL case is PROVED. EA, EB, EP and EQ all sit in E's RE segment, so the two substrate ports
carry the same monomial A·B·[E] and the two product ports P·Q·[E]. The hub is then a chain with
combined ports. RRRRSRRRR (7/7), RRSSSSSRR (11/7) and RRSSSSRR (9/7) all share the law

v = (a·AB − b·PQ)/(1 + K1A + K2B + P/K8 + Q/K9 + c·AB + d·PQ).

This was checked exactly at 36/36 points (`random_effective.py`). All three have the full
7-dimensional orthant as their family, so all are EQUAL. The RE seed with own groups (7/7) is the
identifiable member.

**Answer.** Denis's intuition about monomials holds for fully SS random order, and thermodynamics
holds. But the move is rank-reducing (point 10), not a pure win. In random order, the pure wins
are merges of **RE** isomerizations: every one of the 207 is a phantom removal.

### 5.5 Point 8 (Theorell–Chance ELIM)

**Counts** (EXACT on §4.1–4.3): −2 fitted relative to the merged form when both flanks are SS; −1
when one or both flanks are RE (the fused type follows the default rule).

**Monomials** (PROVED, §7(v)): the removed denominator monomials are exactly those that only the
central complex's occupancy c·u·W'_Y + d·w·W'_Z produced. Examples:

* ordered all-SS: {ABP, BPQ};
* ordered with RE flanks: {AB} or {PQ};
* uni-uni: {S, P};
* ping-pong halves: {AP} or {BQ}.

Checked independently with a matrix-tree computation on 6 all-SS cases (`elim_monomials.py`).
ELIM never adds monomials under the default type rule (EXACT on the catalog).

**Relation:** ELIM ⊂ cl(merged), limit only (PROVED; SAMPLED 3.7e-12 and 8e-13 at |log10 θ| of
11 to 2e8). Denis's claim is right: −2 parameters and a changed equation whenever the flanks are
SS.

The "other" fused type matters:

* RE+SS→RE can be invalid. Ordered RRSR→RRR is an all-RE cycle.
* SS+SS→RE (−3) and RE+RE→SS (0, and it gains monomials) change which forms share an RE segment.
  They are not a storage-removal limit.

### 5.6 Point 10 (param-reducing merges of identifiable mechanisms)

Yes. Such moves exist; every one of them is a strict sub-family (the rank drops), and about half keep all concentration terms (EXACT):

* **Random, SS isomerization.** 39 of the 80 identifiable parents with a valid merge give an
  identifiable result with −2 fitted and identical monomials (for example RSSSSSSSS 14/14 →
  RSSSSSSS 12/12). The other 41 lose some monomials.
* **Ordered.** The 4 identifiable parents give:

  | parent | same-monomial variants |
  |---|---|
  | RRSRS 6/6 | RSRS 5/5 |
  | SRSRR 6/6 | SRSR 5/5 |
  | SRSRS 7/7 | SRSS 6/6, SSRS 6/6 |
  | RRSRR 5/5 | none (RRSR 4/4 loses PQ, RSRR 4/4 loses AB) |

  All are rank −1 sub-families. SAMPLED: RRSRS → RSRS gives 0.25 for P ⊆ R and 4.1e-12 for R ⊆ P only at
  |log10 θ| ≈ 11–13 (a limit); SRSRS → SRSS gives 0.15 and 2.7e-12 at |log10 θ| ≈ 10–14.
* **TC.** ELIM always gives −2 (SS flanks) or −1, with fewer monomials.

**Chain topologies.** Every parameter-reducing move that keeps the family removes a phantom. So
for an identifiable chain parent, a parameter-reducing move necessarily shrinks the family
(PROVED: EQUAL ⇒ same rank).

## 6. Where the package already loses

**Chain lemma on the whole export** (EXACT, `chain_export.py`). A chain qualifies when all three of
its steps are alone in their groups.

| reaction | mechanisms | with ≥1 qualifying chain | engine-A phantoms | phantoms removed by class-rep MERGE | mechanisms made identifiable | representative mismatches |
|---|---|---|---|---|---|---|
| R1 | 4 | 4 | 4 | 4 | 3 of 3 | 0 |
| R2 | 7 | 7 | 6 | 6 | 5 of 5 | 0 |
| R3 | 11 | 11 | 16 | 8 | 3 of 10 | 0 |
| R4 | 1819 | 262 | 616 | 159 | 127 of 418 | 0 |
| R5 | 4715 | 607 | 771 | 297 | 259 of 567 | 0 |
| R6 | 6000 | 660 | 975 | 178 | 159 of 906 | 0 |

In every qualifying mechanism, the class representative has fitted = parent − predicted and
rank = parent rank (fast rank). A further 180, 488 and 665 mechanisms (R4, R5, R6) have a chain
whose flank group is shared with a dead-end binding, for example R4_00001, where EA+B⇌EAB shares
its K with EQ+B⇌EBQ. The lemma does not apply to them.

**Provable duplicate families** (EXACT, `dup_families.py`: group the export by class
representative):

| reaction | exported mechanisms equal in family to another exported mechanism with fewer parameters |
|---|---|
| R1 | 3 of 4 |
| R2 | 5 of 7 |
| R3 | 7 of 11 |
| R4 | 108 of 1,819 |
| R5 | 188 of 4,715 |
| R6 | 80 of 6,000 |

Ranks are equal within every group. Example: R4_00004 (seed, 5/5) ≡ R4_00081 ≡ R4_00082 (6/5) ≡
R4_00467 (7/5).

## 7. Part D: general rules

**(i) When does MERGE of an SS isomerization keep the concentration monomials?**

* *Chains* (PROVED). Write D_P = D_rest(a,b) + c·u·W_Y + d·w·W_Z, where W_Y and W_Z are the
  King–Altman weights of Y and Z in the network with the chain replaced by the edge (a, b). Every
  term has positive coefficients, so no cancellation occurs and the supports do not depend on
  parameter values.
  - The class-representative MERGE keeps every monomial.
  - The plain MERGE that turns SSR into SR (or RSS into RS) loses exactly
    supp(u·W_Y) ∖ (supp D_rest ∪ supp(w·W_Z)) (or its mirror).
  - That set can be empty, so the monomials survive while the family shrinks: ordered RSSRS → RSRS,
    SRSSR → SRSR, SRSSS → SRSS and SSSRS → SSRS (4 of 11 non-EQUAL rows in §4.2) keep all monomials.
  - Or it can be non-empty: RRSSR → RRSR loses PQ.
* *Non-chain* (random): EXACT, 56 of 207 valid plain merges keep all monomials. The one EQUAL case
  also keeps them.

**Rule:** identical monomials are necessary for EQUAL but not sufficient. Counterexample: ordered
RSSRS (7/6) → RSRS (5/5) has identical monomials, yet it is R ⊂ cl(P). SAMPLED: P ⊆ R 0.058; R ⊆ P
only as a limit (≤4e-6 at |log10 θ| ≈ 10–13).

**(ii) When does it keep the family?**

* *Chains* (PROVED): the plain MERGE keeps the family iff the isomerization is RE, or both flanks
  are SS as well (SSS).
  - With exactly one RE flank, R ⊂ cl(P), limit only. Where a monomial is lost this is PROVED: the
    parent law is c·(AB·Keq − PQ)/D with D > 0 on the positive orthant, so D is coprime to the
    irreducible driving force. Supports are therefore the same at every positive parameter point,
    and a law missing a monomial is never attained. Where no monomial is lost (e.g. RSSRS → RSRS),
    "limit only" is SAMPLED: the fits reach the result's laws only at |log10 θ| ≈ 10–13.
  - With both flanks RE (the RE sandwich), the result is invalid, or the fast-chain limit (PROVED
    via the RE limit of both merged steps): ordered RRSRS → RRRS loses {1, A, Q, PQ}; SAMPLED,
    R ⊆ P is reached only at |log10 θ| ≈ 11–13.
  - The family-preserving move is always the class-representative MERGE.
* *General* (PROVED): plain MERGE of an SS isomerization equals the RE-isomerization version of
  the parent, then merged. Hence R ⊂ cl(P) always, with EQUAL only if the rank does not drop.
  EXACT in random: 1 of 207, a chain-like hub.
  - A sufficient condition for EQUAL: the hub has combined ports (all of X1's outside neighbors
    lie in one RE segment and carry the same monomial, and likewise for X2), and all hub-incident
    steps are SS.

**(iii) The "RE sandwich" hypothesis.**

* *True for chains* (PROVED): RSR maps bijectively onto FULL via (K1, f0, r0, K2) ↔ (c, a/c, b/d,
  1/d), and so does the merged SS. So Y⇌X1 RE, X1→X2 SS, X2⇌Z RE is EQUAL to X with an SS entry
  and SS exit, with the **same** count. EXACT (uni-uni RSR, ordered RRSRR, RRSRS, SRSRR, SRSRS)
  and SAMPLED.
* The sandwich's plain merge (RE, RE) is wrong: invalid, or only the fast-chain limit.
* *False when X1 or X2 has more than two steps* (EXACT). Random SSRRSRRSS (E+A and E+B SS; entries
  RE; isomerization SS; exits RE; EP→E+P and EQ→E+Q SS) is 11/11. X with SS entries and exits
  (SSSSSSSS) is 13/13. Both are identifiable with different dimensions, so they are not EQUAL.
* Even when it is EQUAL for more than two steps, it costs parameters. The RE seed RRRRSRRRR (7/7)
  equals RRSSSSRR (9/7), PROVED in §5.4, so the SS-flank version carries 2 phantoms.

**(iv) Is MERGE of a ligand-free RE step always exact (−1, same family)?**

* *Exact* (PROVED, Theorem A1 map) whenever every step incident to X1 or X2 is alone in its group.
  More generally, it is exact whenever every group's steps receive the same rescaling factor p.
* *Corollary* (PROVED): such an RE isomerization is **always a phantom**. EXACT confirmations:
  - 207/207 random RE-isomerization assignments, maps checked;
  - all RE-isomerization chain cases;
  - package seeds R4_00061 and R4_00062.
* *Not exact* when:
  1. A group straddles. Package seeds R4_00056–R4_00060 bind Q to E and to EB_res, and B to EQ and
     to E_res, with shared K. The rescaling would need K·p2 = K/p1, i.e. (1+K0)²/K0 = 1, which is
     impossible. These parents are identifiable (6/6), and the merge (5/5) gives **NEITHER** for
     R4_00057 (SAMPLED 0.037 / 0.014, parameters moderate).
  2. Another step also joins X1 and X2 (MERGE would create a self-loop).
  3. The merged step's own group has other steps (undefined; 133 of 14,120 exported
     isomerizations, all in R5/R6, see §2).

**(v) Which monomials does ELIM remove?** (PROVED, with ELIM applied to the class representative.)

* With both flanks SS, D_merged = D' + c·L1·W'_Y + d·R2·W'_Z. Here D' and W' are the King–Altman
  quantities of the ELIM network, whose fused edge carries the same (a, b).
* So ELIM removes exactly (L1·supp W'_Y ∪ R2·supp W'_Z) ∖ supp D'. With a C, D or RSTOR flank, only
  the corresponding term is present.
* EXACT on 6 all-SS cases via an independent matrix-tree computation: ordered {ABP, BPQ};
  uni-uni {S, P}; ping-pong X1 {AP}; X2 {BQ}; X2 after X1 {AB, BQ, PQ}; six-step with half 2
  merged {BQ}.
* Mixed-type cases match the same rule in every row of §4.2. The default fused type never adds a
  monomial (EXACT, catalog).

## 8. What this suggests for the enumerator (for Denis to decide)

1. **Canonicalize every qualifying chain** (degree-2 X1 and X2, exclusive groups) to its
   class-representative MERGE. This loses nothing (EQUAL, PROVED), needs no rank computation (an
   O(steps) degree check), and removes every chain phantom and the §6 duplicates.
   - Equivalently: emit catalytic chains only in merged form, and let RE→SS act on merged flanks
     (4 flank types) instead of the 8 unmerged chain types.
2. **Merge every RE ligand-free step whose neighbours own their groups:** −1, always a phantom.
   Keep it when groups straddle, as in R4_00056–60, which are identifiable.
3. **Theorell–Chance** needs Steps with ligands on both sides, which the package lacks. The TC
   family is a new family, not a duplicate. Ordered TC with RE outer steps (RRS 3/3) and the RE
   merged variants RRSR/RSRR (4/4) have **fewer** parameters than the RE seed (5/5), so they
   would belong in `init_mechanisms` (point 9).
4. **Random order:** keep the unmerged SS-isomerization forms; they are genuinely larger families.
   Offer their merge as a Δ = −2 move (point 10).
5. **Ping-pong textbook:** no identifiable equivalent exists (§5.3).

## 9. Files (all in the track directory)

| file | contents |
|---|---|
| `transforms.py` | MERGE, ELIM, topology builders |
| `catalog.py` → `catalog_{uniuni,ordered,pingpong,random}.json` | per-assignment f/r, monomials, laws, MERGE and ELIM results |
| `chainlemma.py` | Table L maps, inverses, exact evaluator |
| `verify_chain.py` → `verify_chain.json` | 227 exact pair checks |
| `negctl.py` | negative controls |
| `rank_formula.py` | closed-form rank, 101/101 |
| `verify_re_merge.py` | (iv), 207/207 exact |
| `random_effective.py` | 36/36 exact |
| `pp_gauge.py`, `pp_gauge_check.py`, `pp_effective_check.py` | ping-pong gauge, symbolic and exact |
| `pp_pieces.py` → `pp_pieces.json` | per-law piece test, 10/10 |
| `elim_monomials.py` | (v) |
| `point10.py` | ordered variants of identifiable parents |
| `random_menu.py` | random merged-menu candidates |
| `chain_export.py` → `chain_export.json` | export-wide chain phantoms |
| `dup_families.py` | provable duplicates |
| `famfits.py`, `famfits_run.py`, `fitx.py` → `famfits_*.jsonl`, `famfits_run_*.jsonl` | family fits |
| `tables.py` → `table_*.txt` | the tables above |
| `catalog_summary.json` | compact machine-readable catalog (all 616 assignments with results) |