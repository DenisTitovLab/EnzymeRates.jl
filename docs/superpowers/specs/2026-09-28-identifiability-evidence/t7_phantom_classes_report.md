# Why flip and split children are non-identifiable (track t7_classes)

Everything below was rederived from scratch. The old spec document was not used as evidence. Labels: PROVED (proof given here), EXACT (exact computation over a named finite set), SAMPLED (family fits), CONJECTURE.

## 0. Data and methods

**Populations.** Every exported mechanism with engine-A rank < fitted:

| reaction | phantom mechanisms | phantom directions |
|---|---|---|
| R1 | 3 | 4 |
| R2 | 5 | 6 |
| R4 | 418 | 616 |
| R5 | 567 | 771 |

**Exact kernels.** For each mechanism I evaluated the Jacobian of v with respect to the log free parameters exactly in rationals (engine A's `_point` with the `_Exact` field). I used random rational parameters and fitted+6 random rational concentration points, and took the null space over Q. Each kernel vector was mapped to all group parameters with d log θ_p = Σ_j coef[p][j] d log θ_free,j. This makes the directions independent of the Haldane pivot choice. I repeated this at three parameter points; the intersection of the three kernels gives the point-independent directions. Code: `klib.py`, `analyze_all.py`.

**Peeling classifier.** Repeatedly apply the first operation that keeps the rank and lowers fitted, in this order:
1. revert segment-level zero-flux SS groups;
2. MERGE an RE isomerization;
3. MERGE an SS isomerization;
4. delete one of two parallel channels;
5. ELIM a two-step SS transit;
6. revert one SS group.

The operation that removes a direction names its class. Code: `peel.py`, `classify.py`. A class counts as *new* in a child when the child has more of it than its parent; otherwise it is *inherited*.

**Move re-implementation.** `moves.py` re-implements `_expand_re_to_ss` and `_expand_split_kinetic_group` in Python. It reproduces the R4 export exactly (EXACT):
- level 1: 369/369 mechanisms;
- level 2: 1388/1388 mechanisms;
- recorded first parents: 0 mismatches.

With it I built the full R4 move graph, with every parent of every node:
- to level 3: 5,446 nodes;
- to level 4: 12,448 nodes and 29,146 edges.

**Family fits.** `idkit.family_fit_details` (8 starts, 40 digits); `ff_retry.py` replays a law exactly and adds starts.

*Calibration.* In the Z2 stratum, where the child's law is provably the parent's, 5 of 42 laws gave residuals between 4e-3 and 0.28 with 8 starts. A replay drove one of them to 6e-41. So an isolated nonzero residual is inconclusive. Where it matters I used exact coefficient algebra instead.

## 1. Definitions

- **Segment.** A connected component of the RE steps. The **segment graph** has one node per segment and one edge per SS step (a self-loop when both ends lie in one segment). An SS edge's stoichiometry is σ_e = uptake of the step + u(from) − u(to), where u is the uptake potential along RE paths inside the segment.
- **Catalytic cycle.** A cycle, or self-loop, of the segment graph with σ ≠ 0.
- **Zero-flux SS group.** An SS group none of whose steps lies on a catalytic cycle of the segment graph.
- **Top transit T.** A form that is a singleton segment (all its steps are SS) and whose every exit consumes nothing: T is left only by releases or isomerizations. Examples: EAB, EPQ, ES, EP.
- **Parallel channels.** Two SS steps that join the same two segments and whose 2-cycle is non-catalytic.
- **Operations.** MERGE and ELIM as in the shared context. revert(g) turns SS group g into RE with K = kf/kr. DEL(e) deletes step e.

## 2. Classes (part A)

**Kernel structure (EXACT, all 993 mechanisms).** The point-independent part of the kernel equals span{e_kf + e_kr : g a zero-flux SS group}. Every other phantom direction varies with the parameter point: its fibre is curved, not a scaling.

**Peeling coverage (EXACT).** Single rank-preserving operations remove:
- R4: 613 of 616 directions;
- R5: 765 of 771;
- R1+R2: 10 of 10.

The rest are class X.

**Directions per class** (new/inherited counted relative to the recorded parent):

| class | R4 directions (mechanisms) | new in R4 by move | inherited in R4 | R5 directions | R1 / R2 |
|---|---|---|---|---|---|
| Z zero-flux channel | 336 (188) | flip 296, split 40 | 0 | 336 | – |
| AB absorbed top transit | 147 (147) | flip 100, split 30 | 17 | 227 (80 more inherited through dead-end) | 2 / 4 (2 inherited through dead-end) |
| RI RE isomerization, same monomial | 60 (54) | seed 2, flip 3, split 7 | 48 | 132 | – |
| PC parallel channels | 42 (42) | flip 42 | 0 | 42 | – |
| SB series block of top transits | 28 (14) | flip 28 | 0 | 28 | 2 / 2 |
| X two-segment ping-pong | 3 (3) | flip 1, split 2 | 0 | 6 | – |

- **Pure Z mechanisms** (R4 and R5 alike, EXACT): 174 (144 flip, 30 split); 14 more mix Z with another class.
- **Dead-end move** (EXACT): it never creates a phantom. All 2 R2 and 141 R5 dead-end phantom children have Δfitted = Δrank = 1 and inherit their parent's phantoms.

### 2.1 Symbolic check of every orbit (PROVED)

`orbits.py` checks with sympy that the rate law is unchanged under each stated orbit:

| class | mechanism | orbit |
|---|---|---|
| Z | R4_00470 | (kf1, kr1) → λ(kf1, kr1) |
| AB | R4_00081 | the transit/absorber orbit below |
| PC | R4_00586 | kf6 += t·K3, kf7 −= t·K2 |
| RI | R4_00062 | K4 → λK4, with K2(1+1/K4) and K2/(K4K5) held fixed |
| X | R4_00429 | the affine orbit in 2.7 |

### 2.2 Z: zero-flux SS channel

**Pattern.** An SS group on no catalytic cycle of the segment graph. Two structural subtypes:

**Z1, dead-end (40 directions, all from splits).** A split separates a dead-end abortive binding step from its catalytic twin inside an SS group. In all 40 cases the step's forms have degrees (1, 3). Example R4_00442, a split of R4_00065 (fitted 8, rank 7):
- E+A ⇌ EA (RE)
- E+Q ⇌ EQ (RE)
- EA+B → EAB (SS, g3)
- EQ+B → EBQ (SS, g6; zero flux)
- EA+P ⇌ EAP and EQ+P ⇌ EPQ (RE, shared g4)
- EAB → EPQ (SS)

Kernel: (kf6, kr6) = (1, 1).

**Z2, short circuit (296 directions in 148 mechanisms, all from flips of a minimal pair of groups).** In all 148 cases the parent's chemistry step is already a self-loop of one RE segment: catalysis closes through an RE abortive route (EXACT). The flip pair cuts off a piece of that segment (E, EA, EQ, EAQ, or two-to-three-form pieces). The piece then hangs on two SS steps whose segment 2-cycle is non-catalytic. Example R4_00470, a flip of R4_00078 (parent 6/6, child 8/6):
- E+A → EA (SS g1)
- E+Q → EQ (SS g2)
- EA+B ⇌ EAB (RE)
- EA+Q ⇌ EAQ and EQ+A ⇌ EAQ (RE)
- EQ+P ⇌ EPQ (RE)
- EAB → EPQ (SS)

The route EA ⇌ EAQ ⇌ EQ closes catalysis without E. The kernel is (1,1) on g1 plus (1,1) on g2. kf2 does not occur in the law at all.

**Invisibility (PROVED).** Cha's reduction turns the steady state into a linear King–Altman system on segments. Its effective edge rates are kf·w(from)/W_s and kr·w(to)/W_t. Around any segment cycle, the product of forward over reverse rates is the product of step equilibrium constants times RE weight ratios. For a non-catalytic cycle this product is 1 (Wegscheider). By Hill's theorem the net cycle flux J_κ = (Π₊ − Π₋)·Σ_κ/Σ is therefore zero for every non-catalytic cycle. A step's net flux is the signed sum of the fluxes of the cycles through it, so a step on no catalytic cycle carries zero flux in every steady state. Then kf[X]Π[consumed] = kr[Y]Π[released], and the steady state solves the RE system with K = kf/kr. Uniqueness makes the two laws identical. So v depends on (kf, kr) only through kf/kr: one phantom along (1,1), and revert gives the same law.

**Computational check (EXACT).** Reverting the zero-flux groups:
- preserves the rank and lowers fitted by their number: 188/188 in R4 and in R5;
- returns exactly the parent: 148/148 flip children;
- gives an exported mechanism: 40/40 split children.

**Why the package lets Z2 through.** `_flux_carrying_groups` flags a group when it shares a biconnected block with a chemistry step in the *parent's form graph*. All Z2 groups pass that test. They are zero-flux only in the *child's segment graph*.

### 2.3 AB: absorbed top transit

**Pattern.** A flip of the binding step that forms a central complex makes it a top transit T with two private steps: a binding entry from s1 and a chemistry step to s2. s2 lies in an RE segment, and its weight is set by a constant that nothing else uses (the absorber). A split creates the same pattern when it un-shares the absorber's constant.

**Kernel support (EXACT).** 86 of the class's elementary vectors are exactly {T's SS binding group, T's SS iso group, one RE binding group}.

**Example R4_00081** (flip of R4_00004 = rapid-equilibrium ordered bi-bi, 5/5; child 6/5). Steps:
- E+A ⇌ EA (RE, K1)
- E+Q ⇌ EQ (RE, K2)
- EA+B → EAB (SS, kf3, kr3)
- EAB → EPQ (SS, kf4, kr4 dependent)
- EQ+P ⇌ EPQ (RE, K5)

Laws (engine A, association constants, Den = kf4 + kr3):

child: v = [K1·kf3·kf4/Den]·(AB − PQ/Keq) / (1 + K1·A + K2·Q + (K1·kf3/Den)·AB + (K2·K5 + K1·kf3·kf4/(Keq·kr3·Den))·PQ)

parent: v = K1·K3·kf4·(AB − PQ/Keq) / (1 + K1·A + K2·Q + K1·K3·AB + K2·K5·PQ)

Map child → parent: K3 = kf3/Den and K5' = K5 + K1·kf3·kf4/(K2·Keq·kr3·Den); all positive.

Map parent → child: choose any kr3 > K1·K3·kf4/(K2·K5·Keq), then set kf3 = K3·(kf4 + kr3) and K5_child = K5 − K1·K3·kf4/(K2·Keq·kr3).

So the families are EQUAL (PROVED), and the free kr3 is the phantom. **Invisibility:** the PQ coefficient is the sum of EPQ's own weight and the transit's population fed back through the reverse chemistry step, and nothing else in the law tells the two apart.

**When sharing blocks the absorber.** If the absorber's constant is shared, the phantom disappears. Examples:
- R4_00072: the B group also forms EBQ, so the flip gains rank.
- R4_00085: the P group also forms EAP.

**Structural predicate (EXACT over 1819 R4 mechanisms).** "Two-step top transit with private groups whose iso neighbour is a private RE leaf" flags 104 mechanisms: 100 of the 147 AB-bearing ones plus 4 without. It misses 47, mostly because of group sharing on the entry side.

### 2.4 PC: parallel channels

**Pattern.** A flip set makes both random-order release channels (EP+Q → EPQ and EQ+P → EPQ) or both binding channels (EA+B → EAB and EB+A → EAB) SS. The central complex T becomes a singleton top transit joined to one segment by two channels whose 2-cycle is non-catalytic. Kernel support: exactly the two channel groups (30/30 elementary vectors). Example R4_00586: parent 6/6, child 8/7.

**Proof (PROVED).** The 2-cycle is non-catalytic, so both channels imply the same equilibrium ratio ρ between [T] and [B] (Wegscheider). Each channel's net flux is then J_e = β_e·(ρ[B] − [T]), where β_e is its exit rate from T. Exits from a top transit are concentration-free, so the whole steady state depends on β_e + β_f only.

**Repair.** DEL of one channel gives an EQUAL family. Explicit check for R4_01411: the original law equals the DEL law under kf7' = kf7 + K1·kf5/K2 and kr6 = K1·kf5·kf6/(K3·K4·Keq·kr5), and the reverse map is positive (symbolic check). 30 of the 42 PC children still gain one rank over their parent.

### 2.5 SB: series block of two top transits

**Pattern.** Both central complexes are top transits joined by the SS chemistry step. Bi-bi: EA+B and EQ+P both SS. Uni-uni: R1_00004, all SS (5 fitted, rank 3). Kernel support: the two SS bindings plus the iso (36 elementary vectors).

**Lemma (PROVED).** Take a block X1 ⇌ X2 (iso rates k, k') with entries α1·x1·[s1] → X1 and α2·x2·[s2] → X2, and concentration-free exits β1, β2. Its determinant is Δ = β1β2 + β1k' + β2k. Solving the 2×2 steady state gives:
- [X1] + [X2] = (α1·x1·[s1]·(β2 + k + k') + α2·x2·[s2]·(β1 + k + k'))/Δ;
- net transfer s1 → s2 = (α1·β2·k·x1·[s1] − α2·β1·k'·x2·[s2])/Δ.

So the rest of the mechanism sees only four numbers (d1, d2, f, r), tied by Haldane. A merged single form realizes any positive such tuple: β2' = f/d1, β1' = r/d2, α_i' = d_i·(β1' + β2'). Conversely, for any k > f/d1 and k' > r/d2, the block realizes the tuple with β2 = (f/d1)(k + k')/(k − f/d1), β1 = (r/d2)(k + k')/(k' − r/d2) and α_i = d_i·Δ/(β_j + k + k'). That leaves a 2-dimensional fibre (k, k'): 2 phantoms, and MERGE is EQUAL.

**R1 corollary (PROVED).** All four R1 mechanisms have the same family: v = a·(S − P/Keq)/(1 + d_S·S + d_P·P) with any positive a, d_S, d_P. Their fitted counts are 3, 4, 4 and 5.

### 2.6 RI: RE isomerization between same-monomial forms

**Pattern.** Ping-pong seeds R4_00061/62 (6/5) have an RE second chemistry step EB_res ⇌ EQ. EQ and EB_res carry the same monomial Q. The denominator (after multiplying by B) is B + K1·AB + K2(1 + 1/K4)·QB + (K2/(K4K5))·Q + (K2K6/(K4K5))·PQ, and the numerator is K1·kf3·(AB − PQ/Keq)/B. So K2, K4 and K5 enter only through two combinations.

Of the 60 R4 directions, 48 are inherited from these seeds. The 10 new ones appear when a split un-shares the neighbour groups (the identifiable seeds R4_00059/60) or when a flip lets the RE step's constant be absorbed. MERGE (EB_res = EQ) keeps the rank in 54/54 mechanisms (EXACT).

### 2.7 X: two-segment ping-pong

**Pattern.** Segment 1 = {E, EA, EQ} with W1 = 1 + K1·A + K2·Q. Segment 2 = {E_res, EB_res, EP_res} with W2 = 1 + K5·B + K6·P. The two are joined only by the SS chemistry steps EA → EP_res and EB_res → EQ. Examples: R4_00429 (7/6), R4_01758 and R4_01777 (8/7).

**Proof (PROVED).** Write a = kf3·K1, q = kr4·K2, b = kf4·K5, p = kr3·K6. Then:
- D = (b·B + p·P)·W1 + (a·A + q·Q)·W2;
- N = ab·AB − pq·PQ;
- the Haldane constraint is ab/(pq) = Keq.

Under (K1, K2) += t·(a, q) and (K5, K6) −= t·(b, p), with a, q, b, p held fixed, every cross coefficient is unchanged; for example AB: b·K1 + a·K5. N and the Haldane constraint are unchanged too. The kernel vector matches: K1 moves against kf3, K2 against kr4, K5 against kf4, K6 against kr3.

No single revert, MERGE or ELIM keeps the rank (0/3). This phantom is intrinsic to the two-segment ping-pong.

## 3. Lossless rejection (part B)

### 3.1 Criterion (PROVED)

A flip's parent is the RE limit of the child (kf, kr → ∞ at fixed ratio). A split's parent is a specialization of the child (both parts' constants set equal). In both cases the parent's family lies in the closure of the child's family. Both families are images of irreducible parameter spaces, so their Zariski closures are irreducible and have dimension equal to the generic rank. Therefore rank(child) = rank(parent) exactly when the Zariski closures are equal. Equal closures do **not** imply equal families over positive parameters (section 3.4).

### 3.2 Rank change against the recorded parent (EXACT, R4 children with new phantoms)

| class, move | Δrank = 0 | Δrank > 0 |
|---|---|---|
| Z flip | 148 | 0 |
| Z split | 8 (together with AB) | 32 |
| AB flip | 98 (12 of them together with PC) | 2 |
| AB split | 30 | 0 |
| PC flip | 12 | 30 |
| SB flip | 13 | 1 |
| RI flip / split | 3 / 7 | 0 |
| X flip / split | 0 / 2 | 1 / 0 |

49 more children only inherit their phantoms (28 flips, 21 splits, all Δrank > 0).

### 3.3 Family fits, child → parent (SAMPLED; 3 laws, 8 starts; at least 15 pairs per stratum, all if fewer)

| stratum | pairs | pairs with all laws < 1e-20 | laws < 1e-20 | largest residual |
|---|---|---|---|---|
| Z2 flip, Δrank 0 (proved equal) | 15 (1 timed out) | 9 | 37/42 | 0.28 (all optimizer misses) |
| AB flip, Δrank 0 | 15 | 11 | 40/45 | 0.19 |
| AB+PC flip, Δrank 0 | 12 | 11 | 35/36 | 0.085 (retry: 3.6e-37) |
| SB flip, Δrank 0 | 15 | 14 | 42/45 | 0.02 (R4_01814, ping-pong, persists with 16 extra starts) |
| RI flip, Δrank 0 | 3 | 3 | 9/9 | 6e-37 |
| X split, Δrank 0 | 2 | 2 | 6/6 | 3e-38 |
| AB split, Δrank 0 | 15 | 1 | 28/45 | 0.25 |
| AB+Z split, Δrank 0 | 8 | 0 | 6/24 | 0.11 |
| RI split, Δrank 0 | 7 | 2 | 14/21 | 0.12 |
| PC flip, Δrank 1 | 15 | 0 | 0/45 | min-of-max 1.2e-3 |
| Z1 split, Δrank 1 | 15 | 0 | 0/45 | min-of-max 7.9e-3 |
| AB, SB, X flips, Δrank 1 | 4 | 0 | 0/12 | min-of-max 0.32 |

Retries with 40 extra starts:
- R4_01751 (AB flip, ping-pong): stays at 0.14 and 0.083.
- R4_00847: stays at 0.058; section 3.4 proves the gap exact.

### 3.4 Positivity gaps (PROVED)

**R4_00847 (AB flip) versus parent R4_00203.** The child has E+B ⇌ EB, E+P ⇌ EP, EAB → EPQ (SS), EB+A → EAB with EP+A → EAP (SS, shared) and EP+Q → EPQ (SS). The parent is the same with EP+Q ⇌ EPQ at RE. Both laws have numerator monomials {AB, PQ}, denominator monomials {1, B, P, AB, AP, PQ} and rank 6.
- For every positive parameter point the parent has d_AP·d_B − d_P·d_AB = K2³K5²Keq²kr3²/(kf3·(K1·kf4 + K2·K5·Keq·kr3)) > 0.
- The child at K1 = 3, K2 = 1, kf3 = 2, kr3 = 1, kf4 = 2, kr4 = 1, kf5 = 1, Keq = 2 has d_B = 3, d_P = 1, d_AP = 2, d_AB = 38/3, so d_AP·d_B − d_P·d_AB = −20/3.

The child therefore reaches laws the parent cannot reach.

**R4_00490 (AB split of R4_00085).**
- Parent: d_PQ − d_Q·d_AP/d_A = K1·kf3·kf5/(Keq·kr3·(kf5 + kr3)) > 0.
- Child: d_PQ − d_Q·d_AP/d_A = the same term + K2·(K6 − K4), which is negative for K4 ≫ K6.

The child is still EQUAL to its exported unflipped split sibling, E+A ⇌ EA, E+Q ⇌ EQ, EA+B ⇌ EAB, EA+P ⇌ EAP (K4), EQ+P ⇌ EPQ (K6), EAB → EPQ (6 fitted). The sibling's coefficients (K1, K2, K1K4, K1K3, K2K6, K1K3·kf) cover the whole positive orthant. The child covers it too: K1 = d_A, K2 = d_Q, K4 = d_AP/d_A, kf5 = n/d_AB; pick kr3 with n/(Keq·kr3) < d_PQ, then K6 = (d_PQ − n/(Keq·kr3))/K2 and kf3 = d_AB·(kf5 + kr3)/K1.

### 3.5 Descendants (EXACT, full R4 move graph, every parent of every node)

| rule | level-3 nodes | identifiable | phantom | identifiable lost within level 3 | lost within level 4 |
|---|---|---|---|---|---|
| none | 5446 | 3759 | 1687 | 0 | 0 |
| drop children without rank gain | 3998 | 3695 | 303 | 64 | 272 (the 64 are reached at level 4) |
| reject phantom children | 3697 | 3695 | 2 | 64 | 280 |
| reject zero-flux children | 4426 | 3695 | 731 | 64 | 272 |
| zero-flux completion of flip sets; drop zero-flux split children | 6420 | 5567 | 853 | 0 | 0 (level-4 variant: 12634 identifiable, 2090 phantom) |
| zero-flux completion + drop children without rank gain | 5726 | 5567 | 159 | 0 | not built |

- The 64 lost level-3 mechanisms all have 9 parameters (and rank 9). Each is reachable only through a level-2 zero-flux intermediate, e.g. R4_00575 (8/6), then one more flip turns the zero-flux groups flux-carrying.
- Every loss under the rank-based rules goes through such an intermediate.
- Zero-flux completion also reaches 1808 identifiable 10–12-parameter mechanisms by level 3.

## 4. Repairs (part C)

### 4.1 Which operations keep the rank (EXACT; mechanisms with at least one rank-preserving instance / mechanisms where the operation applies; class = first peel class)

| class (n) | revert SS group | MERGE SS iso | MERGE RE iso | ELIM | delete one channel |
|---|---|---|---|---|---|
| Z (188) | 188/188 | 0/188 | 6/16 | 0/108 | 140/148 |
| AB (124) | 124/124 | 0/124 | 0/10 | 0/120 | 0/20 |
| PC (42) | 12/42 | 0/42 | – | – | 42/42 |
| SB (15) | 15/15 | 15/15 | 0/1 | 15/15 | – |
| RI (54) | 14/54 | 1/54 | 54/54 | 1/15 | 0/8 |
| X (3) | 0/3 | 0/3 | – | – | – |

### 4.2 EQUAL tests (family fits both ways, 5 laws each, 8 starts; plus proofs)

| class | pairs | result | proof |
|---|---|---|---|
| Z | 6 | 5/5 both ways on 5 pairs; the sixth has one 6e-2 law | Hill argument (2.2) |
| AB (revert) | 8 | 5/5 both ways on R1_00002, R1_00003, R4_00733; 1–3 nonzero laws on R4_00737, R4_00788, R4_01383, R4_01704, R4_01707 | R4_00081 EQUAL (2.3) |
| PC (delete one channel) | 6 | orig → repaired has 1–3 nonzero laws on 5 of 6 pairs | R4_01411 EQUAL (2.4), so these are optimizer misses |
| SB (MERGE) | 8 | repaired → orig 5/5 on every pair; orig → repaired 5/5 on R4_00481 and R4_00839, 4/5 (6e-6, 2e-5, 2e-5) on three others; R4_01692 timed out | two-form block lemma (2.5) |
| RI (MERGE) | 6 | 5/5 both ways on 4 pairs; one 5e-3 or 3e-2 law on the other 2 | orbit verified (2.1) |

**Representability (EXACT).** None of the 67 MERGE repairs and none of the 42 channel deletions is an exported mechanism. All 188 Z repairs and 111 of 118 revert repairs are. A merged complex holds two different bound sets, which one `Species` cannot express.

## 5. Summary table (part D)

| class | R4 directions (mechanisms) | origin | child ⊆ parent? | lossless rejection | exact repair (EQUAL) | simple rule |
|---|---|---|---|---|---|---|
| Z1 dead-end | 40 (40) | split of an SS group holding a dead-end abortive step and its catalytic twin | no (Δrank 1 in 32) | yes: the zero-flux-reverted mechanism is exported (40/40) | revert the group | exact segment-level zero-flux test; drop the child |
| Z2 short circuit | 296 (148) | flip pair that isolates a piece of a segment whose chemistry is a self-loop (RE abortive return path) | equal to the parent (148/148) | for that parent yes; descendants no (64/272 delayed) unless zero-flux completion | revert, which gives the parent | exact test plus completion; lossless to level 4 |
| AB absorbed transit | 147 (147) | flip of the binding step that forms a central complex; split that un-shares the absorber | same Zariski closure in 128/130; positive family: mostly for flips, not for R4_00847 or most splits | usually, via the parent or the exported unflipped sibling; not guaranteed | revert the binding group; MERGE and ELIM change the family | none exact; local predicate 100/147 (4 false positives); exact test = rank |
| PC parallel channels | 42 (42) | flip of both random release (or binding) channels of a central complex | 12 yes; 30 carry a real rank gain | no for those 30 | delete one channel (not emitted today) | exact predicate (non-catalytic 2-cycle into a top transit) |
| SB series block | 28 (14) + R1/R2 4 | flip of the second binding around the central step (both complexes transits) | yes (13/14) | yes | MERGE the central step (−2) | exact predicate (SS iso between two top transits); needs a merged-complex representation |
| RI RE isomerization | 60 (54) | ping-pong seeds R4_00061/62; splits that un-share neighbours | yes (10/10 new) | yes | MERGE the RE step | seed-level fix |
| X two-segment ping-pong | 3 (3) | flip of the RE second chemistry step; splits | 2 yes, 1 no | – | none | unfixable by simple operations |

**What remains unfixable by any simple rule:**
- **X**, which is intrinsic.
- **AB**: an exact test needs the absorber analysis or a rank comparison, and even equal rank can hide a positivity gap.
- **PC children that gain rank**: they need a topology-changing variant (delete a channel).
- **Positivity differences between same-rank families**, which no rank-based rule sees.

## 6. Files (all in /tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/8633b908-eb4b-4a04-bd9a-bfa7049838d2/scratchpad/ident/t7_classes/)

- **Scripts:** klib.py (exact kernels, segment-level zero-flux test), ops.py, analyze_all.py, stage2.py, stage3.py, peel.py, classify.py, tfeat.py, moves.py (flip and split re-implementation plus variants), validate_moves.py, build_graph.py, reach.py, build_variant.py, orbits.py, posimage.py, repair_table.py, pred_ab.py, mono_all.py, count_table.py, ffrun.py, ff_retry.py, make_jobs_B.py, make_jobs_C.py.
- **Data:** r4.jsonl, r5.jsonl, small.jsonl, s2_*.jsonl, s3_*.jsonl, peel_*.jsonl, classes_*.jsonl, graph_R4_L3.json, graph_R4_L4.json (+ _reach.json), variant_*.json, ff_B.jsonl, ff_C.jsonl, ff_retry*.jsonl, jobs_C_meta.json, repair_table.txt, mono_all_R4.json, r4_edges.json.
