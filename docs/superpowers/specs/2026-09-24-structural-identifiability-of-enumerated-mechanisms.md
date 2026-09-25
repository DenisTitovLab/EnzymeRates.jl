# Structural identifiability of enumerated catalytic mechanisms

Date: 2026-09-24 (fusion analysis 2026-09-25). Status: survey, catalog of causes and verified
fixes; no design decided.

## Goal and scope

Denis's goal: the enumeration engine emits only **structurally identifiable** rate equations.
A rate equation is structurally identifiable when every fitted parameter is determined by ideal
steady-state rate data, with all substrates and products at arbitrary positive concentrations.
Practical identifiability depends on the data actually collected; it is out of scope here and
could become separate, data-driven functionality. Only catalytic mechanisms (`Mechanism`) are in
scope: an MWC wrapper is identifiable by construction, apart from the catalytic mechanism inside
it, which this analysis covers.

Terms used below:

- **Fitted parameters**: `fitted_params(m)`, the constants left after the Haldane and
  Wegscheider constraints remove dependent ones. This is the count the beam charges.
- **Rank**: the structural rank of `∂ ln v / ∂ ln θ`, the number of independent directions the
  rate law can see.
- **Phantom**: a fitted direction the rate law cannot see; phantoms = fitted − rank.
- **Family**: the set of rate laws a mechanism reaches over all positive parameter values. A
  fix is **family-preserving** when the fixed mechanism reaches every rate law of the original.

## Summary

- **Non-identifiability is common.** Of the 37,264 catalytic bi-bi mechanisms the moves reach
  within the beam's caps (up to 14 fitted parameters), 14,873 (39.9%) carry at least one
  phantom, rising from 0% at 5 parameters to 50% at 14. Seeds are almost always identifiable;
  phantoms enter through RE→SS flips (40–45% of sampled flip children) and splits (33–58% of
  sampled split children). Every uni-uni flip child has a phantom.
- **Every phantom has a known structural cause.** All 559 phantom directions in the bi-bi and
  uni-uni profiles fall into ten classes, each with an exact explanation. Five classes account
  for 95% of them: zero-flux SS groups (EQ), Briggs–Haldane transit forms (BH), serial transit
  chains (CHAIN), RE-chemistry lumping in ping-pong (PAIR) and the dwell-time gauge (CUT). PGK
  adds one more family, duplicate inhibitor binding modes, because it lists its own substrates
  and products as competitive inhibitors.
- **Every class has a family-preserving fix**, proved exact for all but one (HYBRID-a, where the
  evidence is numerical) and confirmed with the package's own derivation on representative
  mechanisms. Fixes take two forms: turning a specific SS group back into RE (or fusing two
  central complexes), or adding one dimensionless constraint row. PGK's duplicate modes are the
  exception: rejecting the offending move usually works, and no general fix is known for the
  rest.
- **Rejecting moves is not enough.** On the complete bi-bi closure (37,264 mechanisms), a
  gate that rejects children whose rank does not rise never cuts off a full-rank mechanism, but
  it loses roughly 0.7–5.5% of the rate-law families and still emits non-identifiable
  mechanisms; emitting only full-rank mechanisms loses about 10%. A rank check alone cannot
  validate a fix either: several edits restore full rank and still lose rate laws.
- **Fusing steps is exact only in specific places (§5).** Merging consecutive central complexes
  and merging the ping-pong RE half-reaction keep every rate law and remove phantoms
  (Michaelis–Menten all SS 5 → 3 parameters, ordered bi-bi 9 → 7, ping-pong 11 → 7). Eliminating
  a central complex (Theorell–Chance) or merging an SS isomerization between branch points (random
  bi-bi) gives a different, simpler rate law. Even the textbook ping-pong keeps one unidentifiable
  parameter (6 identifiable of 7). The package derives some fused steps correctly but returns
  silently wrong rate laws for others; those bugs come first.
- **Canonicalization reaches the goal.** Expanding as today but emitting each mechanism's
  canonical identifiable form (its RE-reduced core plus exact slice rows) removes every phantom,
  loses almost no family, and fits about a third fewer mechanisms. It needs new machinery: SS→RE
  group conversion, dependent slice rows, class detection, merge normalization, and some
  canonical forms the moves never build today.

## 1. Method

**Mechanisms.** Catalytic mechanisms were generated with the enumerator's own code:
`init_mechanisms` seeds, then the catalytic moves `_expand_re_to_ss`,
`_expand_split_kinetic_group` and `_expand_add_dead_end_regulator` (allosteric moves excluded).
Five reactions: uni-uni S ⇌ P; uni-uni with a competitive inhibitor I; a generic bi-bi
A + B ⇌ P + Q whose atoms allow a ping-pong residual; that bi-bi with an inhibitor I; and PGK
(BPG + ADP ⇌ 3PG + ATP, with its shared catalytic sites and all four metabolites declared
competitive inhibitors). Uni-uni was enumerated exhaustively. The other reactions were explored
by breadth-first search through two expansion levels plus random walks from the seeds, up to 14
fitted parameters, and then sampled at 60–100 mechanisms per parameter count. Mechanisms above
the beam's `eq_complexity_filter` (V×τ > 337) were excluded from all headline numbers, since the
beam never fits them.

**Rank.** For each mechanism the package's `rate_equation` was differentiated in BigFloat
(256-bit, central differences) with respect to the log of every fitted parameter, at 3n random
concentration points with all metabolites present, parameters and concentrations drawn within
10^±1 of each other; rank is the count of singular values above 1e-9 of the largest, maximized
over up to three draws. Controls matched their known ranks: RE random bi-bi 5 of 5,
context-shared iso 7 of 7, forward mnemonic with context-split binding 12 of 13, context-split
flip-flop 15 of 16, fully SS random bi-bi 15 of 15. Where a phantom exists, the gap between kept
and dropped singular values is about 10¹². Concentrations far below the binding constants fake
rank loss; the first control run caught this.

**Independent checks.** Two independent exact engines written in Python (Cha segments with a
matrix-tree steady state, and forward-mode differentiation modulo a large prime) reproduce the
package's fitted count and rank on all 1,530 profiled bi-bi and uni-uni rows; one of them ranks a
mechanism in about 20 ms, which allowed ranking the whole generated bi-bi population. Fixes were
checked for family preservation by fitting the fixed mechanism to exact rates generated from
random parameter draws of the original. The fits were checked against known answers:
self-reproduction reaches 10⁻¹⁵, the RE Michaelis–Menten parent and its SS-binding flip child
reproduce each other at 10⁻¹⁴, and RE random bi-bi fails to reproduce SS random bi-bi at 0.02.

## 2. Extent

Mechanisms with at least one phantom, by reaction (profiled samples, beam-fittable only):

| Reaction | Profiled | With phantoms | Phantoms per mechanism | Ping-pong | Not ping-pong |
|---|---|---|---|---|---|
| uni-uni | 4 | 3/4 (75%) | 1.00 | — | 3/4 (75%) |
| uni-uni + inhibitor | 8 | 6/8 (75%) | 1.00 | — | 6/8 (75%) |
| bi-bi | 870 | 331/870 (38%) | 0.63 | 134/209 (64%) | 197/661 (29%) |
| bi-bi + inhibitor | 595 | 189/595 (31%) | 0.56 | 55/121 (45%) | 134/474 (28%) |
| PGK | 572 | 248/572 (43%) | 0.79 | 67/101 (66%) | 181/471 (38%) |

The complete bi-bi closure (all 37,264 fittable mechanisms up to 14 parameters, ranked with the
fast independent engine; §6) gives 39.9%, matching the sample's 38%, but it rises more slowly
with parameter count than the sample does: 6% at 6 parameters, 25% at 8, 35% at 10, 44% at 12
and 50% at 14. The random walks over-represent phantom mechanisms at middle counts. Ping-pong
mechanisms are about twice as likely to carry phantoms as sequential ones.

By fitted-parameter count:

| Fitted params | bi-bi | bi-bi + inhibitor | PGK |
|---|---|---|---|
| 5 | 0/55 (0%) | 0/55 (0%) | 0/32 (0%) |
| 6 | 10/100 (10%) | 4/60 (6%) | 5/60 (8%) |
| 7 | 9/100 (9%) | 6/60 (10%) | 16/60 (26%) |
| 8 | 38/100 (38%) | 9/60 (15%) | 22/60 (36%) |
| 9 | 38/97 (39%) | 12/60 (20%) | 25/60 (41%) |
| 10 | 52/96 (54%) | 24/60 (40%) | 28/60 (46%) |
| 11 | 51/92 (55%) | 29/60 (48%) | 36/60 (60%) |
| 12 | 49/87 (56%) | 35/60 (58%) | 37/60 (61%) |
| 13 | 40/69 (57%) | 36/60 (60%) | 36/60 (60%) |
| 14 | 44/74 (59%) | 34/60 (56%) | 43/60 (71%) |

By the move that created the mechanism:

| Last move | bi-bi | bi-bi + inhibitor | PGK |
|---|---|---|---|
| seed | 2/61 (3%) | 0/57 (0%) | 0/32 (0%) |
| re_to_ss | 227/501 (45%) | 114/282 (40%) | 49/108 (45%) |
| split | 102/308 (33%) | 71/198 (35%) | 101/172 (58%) |
| dead_end_reg | — | 4/58 (6%) | 98/260 (37%) |

Seeds are identifiable except for 2 bi-bi ping-pong seeds (4 of 69 in the whole bi-bi
population), all of the PAIR class. In uni-uni, the enumerator builds only four catalytic
mechanisms: Michaelis–Menten with RE binding (identifiable) and its three RE→SS flip children
(all non-identifiable: a Briggs–Haldane transit form when one binding is SS, a transit chain
when both are). In PGK, dead-end inhibitor moves create phantoms far more often than in the
generic bi-bi with an inhibitor (37% against 6%), through the duplicate binding modes described
in §3.

## 3. Why constants become invisible

**The common mechanism.** Almost every class is a reweighting of enzyme forms: the population of
a form X is multiplied by (1 + ε_X) while every constant leaving X is divided by the same
factor. The one-way fluxes out of X are unchanged, so if the weighted populations still sum to
the same total at every steady state, the rate law cannot tell the difference. That happens in a
few structural situations: a step that carries no net flux, a ligand-free cut between two
blocks of forms, a fixed RE ratio between two forms, or an intermediate that the rate law sees
only through its inflow and outflow. Sharing a constant between kinetic groups blocks a
reweighting whenever the shared steps would have to move by different factors. That is why the
enumerator's context-shared seeds are almost always identifiable, and why splits expose
phantoms.

Counting the rate law's coefficients (its distinct concentration monomials) explains only 37 of
331 phantom bi-bi mechanisms. Most phantoms are algebraic dependences among coefficients, not
missing terms.

Classes, with the phantom directions they explain in the uni-uni and bi-bi profiles (559
directions in 340 mechanisms; one mechanism can carry several classes):

| Class | What makes the constants invisible | Directions | Family-preserving fix | Status |
|---|---|---|---|---|
| EQ | An SS group carries no net flux, so it sits at equilibrium | 181 | Make the group RE (on a flip child: reject the flip), or pin kon = 1 | Proved; package-checked |
| BH | Briggs–Haldane: a transit form is seen only through its inflow and outflow | 133 | Make the transit form's binding entry RE | Proved for Michaelis–Menten; package-checked |
| CHAIN | Consecutive transit forms, e.g. SS binding, SS chemistry, SS release | 131 | Make each transit form's entry RE, working inward from the chain's end; or merge the central complexes (§5), an equivalent exact fix | Proved; package-checked |
| PAIR | Two forms joined by an RE chemistry step trade mass freely (ping-pong) | 42 | Merge the two forms (§5), or slice Kiso = 1 | Proved (merge: when no shared group touches the pair); package-checked |
| CUT | Dwell-time gauge across a ligand-free SS cut (ping-pong with both half-reactions SS) | 41 | Slice: product of the crossing steps' K = 1; for ping-pong also a union of three identifiable forms split by a sign invariant (§5) | Proved; package-checked |
| CHAIN-b | A chain one of whose entry groups also holds a dead-end step | 9 | Make the chain's chemistry RE (fuse the two central complexes), or slice k_chem² = koff_in·koff_out | Orbit proved; package-checked |
| PAR | Two SS bindings enter the same complex in parallel | 9 | With BH: make both entries RE. Alone: slice koff₁ = koff₂ | Orbit proved; package-checked |
| HYBRID-a | A BH-type collision next to a context-shared binding square (mostly ping-pong) | 9 | Make the transit form's binding entry RE, confirmed by a rank check | Package-checked on 2; not proved |
| RING (HYBRID-b, JOINT) | Ring gauge: in a ping-pong ring, the transit blocks on opposite arcs carry equal flux, so four exit lifetimes shift together | 4 | Slice: product of the two rising exit rates = product of the two falling ones (plus Kiso = 1 when an RE-iso pair is involved). JOINT alone: rejecting it is lossless if both of its RE parents stay in the beam | Proved; invariance checked with the package's `rate_equation` |
| TRI | Fully shared random bi-bi with one SS substrate and one SS product binding | 0 (4 in the population) | Slice k_f·k_r = koff_B·koff_Q | Proved; package-checked |

"Proved" means an exact orbit argument shows the fix crosses or reaches every orbit of
equivalent parameter sets. "Package-checked" means that on two or three representative
mechanisms the fixed mechanism was built with the package, derived at full rank, and reproduced
every tested rate law of the original (worst loss 4.3e-11), while at least one plausible but
wrong edit failed. Package checks are evidence, not proof: with only a few random rate laws per
test, a fix that loses a small region can pass. One lossy pin (k_chem = koff_in for CHAIN-b)
passed 6 of 6 laws yet loses 13% of rate laws in the same parameter box. Every exactness claim
therefore rests on an orbit proof.

### Details by class

**EQ: zero-flux SS group.** A step that lies on no catalytic cycle of the segment graph (RE
segments contracted to nodes) carries zero net flux at steady state, so it satisfies its own
equilibrium relation and only K = koff/kon enters the rate law; the common scale of kon and
koff is invisible. Two shapes occur: a dead-end or abortive SS binding (a bridge of the form
graph), and an SS route running parallel to an RE route between the same segments while
catalysis runs elsewhere. It is the largest class. It comes from RE→SS flips of groups whose
steps are all off catalytic cycles, and from splits that separate the zero-flux context of a
shared group. The segment count test in `_expand_re_to_ss` does not catch it, because the
dead-end form becomes its own segment. Seeds never contain it.

**BH: Briggs–Haldane transit form.** A transit form is an enzyme form entered by SS binding and
left only by SS steps that release a ligand or isomerize. Its population is its total inflow
divided by its total exit rate, and the part driven by a ligand-free neighbour merges with that
neighbour's population, so one of its four port quantities is invisible. Canonical case:
Michaelis–Menten with SS substrate binding, whose family is exactly the RE parent's family
(proved). Every uni-uni flip child with one SS binding is BH.

**CHAIN: serial transit forms.** Two or more transit forms in a row, typically
E(A,B) <--> E(P,Q) with the entry binding, the chemistry and the exit release all SS. Only the
block's inflows and outflows are visible. Uni-uni with every step SS is Cleland's classic case:
rank 3 of 5, the same family as the RE seed. An ordered bi-bi with all steps SS has rank 7 of 7
with one central complex, 7 of 9 with two and 7 of 11 with three, and the two-complex family
equals the one-complex family.

**PAIR: RE-chemistry lumping.** In a ping-pong seed, the RE half-reaction E(A) ⇌ F(P) fixes the
ratio of the two forms at every concentration, so mass can move between them with each
attached constant absorbing the change. Kiso sweeps (0, ∞) along the orbit, so the slice
Kiso = 1 meets every orbit once. It is the only class present in seeds.

**CUT: dwell-time gauge.** Deleting every ligand-free SS step disconnects the forms into two
blocks (in bi-bi: ping-pong with both half-reactions SS). The rate fixes the total round-trip
time between the blocks but not its split (`2026-09-23-iso-kinetic-cooperativity-findings.md`,
§9). The slice "product of
the crossing steps' equilibrium constants = 1" meets every orbit exactly once.

**CHAIN-b.** A two-complex chain whose entry group also contains a dead-end step. The dead end
pins the entry's K, so reverting a binding to RE loses rate laws, but every orbit ends where the
chemistry becomes infinitely fast: making the chemistry RE (fusing the central complexes) is
exact, and so is the slice k_chem² = koff_in·koff_out. The package accepts the RE-chemistry
sequential mechanism, but the enumerator never generates one.

**PAR: parallel entries.** Two SS bindings in different groups land on the same complex (for
example E(A) + B and E(B) + A into E(A,B)); only their summed exit and entry are visible.

**HYBRID-a.** A context-shared random-order binding square blocks the simple BH pattern. It is
usually ping-pong (every bi-bi instance), but PGK has two sequential instances. It behaves like
BH: making the transit form's binding entry RE reproduced every tested rate law, but no proof
exists, and only a rank check licenses the conversion.

**RING (HYBRID-b and JOINT).** In a ping-pong ring, the transit blocks on either side of a hub
form carry the same steady-state flux, so their port masses can be traded without changing v:
two exit lifetimes rise by λ while two fall by λ, and every other constant follows at a fixed
ratio (proved on all four profiled mechanisms; the package's `rate_equation` confirms the
invariance to 10⁻⁷⁵). Each orbit end is where one lifetime reaches zero, and which end depends
on the parameter point, so no single SS→RE conversion covers the family: each conversion
contains only the rate laws on one side of a point-dependent inequality. The ratio
h = (product of the rising lifetimes)/(product of the falling ones) runs from 0 to ∞ along
every orbit, so the log-linear row "product of two exit rates = product of the other two" is
crossed exactly once. HYBRID-b always co-occurs with PAIR and needs Kiso = 1 as well. JOINT's
family is exactly the union of its two RE parents' families (P binding RE and B binding RE), so
rejecting JOINT loses nothing as long as both parents stay in the beam; HYBRID-b has no such
cover. On the generated populations, an automatic rule (pick the four SS lifetimes that move at
equal and opposite speed along the null vector) restored full rank on all 29 JOINT and 9
HYBRID-b mechanisms.

**TRI.** A fully context-shared random bi-bi with one SS substrate binding and one SS product
binding forms a triangle of three RE segments closed by the chemistry; the dimensionless slice
k_f·k_r = koff_B·koff_Q is exact. It does not occur in the profiled sample; the population holds
four such mechanisms, all relabelings of one another.

### PGK: duplicate inhibitor binding modes (INH-DUP)

PGK declares its own substrates and products as competitive inhibitors, so the dead-end move
can add a second binding mode of the same ligand to forms where it already binds catalytically.
The inhibitor-bound form (for example E(ATP::Inh)) then carries the same concentration monomial
as the catalytic form (E(ATP)), and the rate law sees only 1/K_ATP + 1/K_ATP,inh. The catalytic
constant's other appearances are absorbed by a constant that appears nowhere else: the
Haldane-dependent reverse chemistry, a downstream binding constant, or SS constants. Proved on
representatives: the rate law depends on the two constants only through their reciprocal sum.

The catalytic core (the mechanism with inhibitor forms removed) carries the same classes as
bi-bi. The new family accounts for most of PGK's phantoms:

| Subclass | Profile directions | Population directions | What it is |
|---|---|---|---|
| DUP-RE | 178 | 671 | Both duplicate forms sit in the same RE segment |
| DUP-T | 55 | 191 | The catalytic partner is an SS transit complex |
| DUP-BH | 27 | 108 | A duplicate combined with a Briggs–Haldane collision |
| DUP-PAR | 8 | 51 | The duplicates leave through parallel SS steps into the same segment |
| INH-EQ | 3 | 7 | An inhibitor complex opens an RE bypass, so a catalytic SS step carries no flux |

Of 2,999 generated PGK mechanisms, 1,067 carry phantoms (1,836 directions): 807 directions from
the catalytic-core classes, 1,029 from INH-DUP, and one unexplained direction (mechanism
15013533285179965847). An inhibitor mode whose forms share no concentration monomial with
another form never adds a phantom (0 of 508 such modes), but duplication alone is not sufficient
(209 of 429 fully duplicated modes add one).

Fixes, and their limits:

- **Rejecting the dead-end child** whose rank rises by less than its parameter count kept
  every tested rate law on 149 of 164 population edges, and on 92 of 92 when it is the first
  inhibitor mode on a purely catalytic parent. The evidence is weak (two sampled laws per edge);
  for the first-mode DUP-RE case the orbit argument supports it, and the other first-mode
  subclasses are undetermined. The losses occur when an earlier inhibitor mode is already
  present.
- **A dimensionless slice** fixed the one lossy representative tested, 3267908235929210394
  (K_ATP_E·K_ThreePGinh_EATP = K_ATPinh_E·K_ThreePG_EATP; proved to cross every orbit once).
  Rejecting that child keeps only the rate laws on one side of an inequality, about half of a
  symmetric parameter box. No general slice rule for INH-DUP is known yet.
- **Rejecting splits or flips that create INH-DUP** is often lossy (at least 21 of 71 splits and
  9 of 23 flips).
- **No local graph rule** separates phantom-creating inhibitor modes from useful ones; the best
  catches 120 of 164 but also fires on 49 identifiable modes. An exact rank check (about 10 ms
  per mechanism) is the reliable detector.
- **A modeling question for Denis.** Forbidding dead-end modes of a reaction's own substrates
  and products at forms where the same ligand binds catalytically would remove most DUP-RE, but
  it would also drop 1,145 identifiable modes that enlarge the family (for example
  1086605419517338646, where a context-dependent ADP affinity is thermodynamically consistent
  only through the dead-end mode).

## 4. Routes to identifiability

Every fix is one of three kinds.

1. **Avoid.** Do not emit the mechanism, or emit a canonical representative of its family.
   This is exact only when the rejected mechanism's family lies inside the family of something
   that is emitted.
2. **Quotient.** Replace the mechanism with an equivalent one that has fewer constants: turn a
   zero-flux or transit group into RE, fuse two central complexes into one (Cleland's net rate
   constants), or merge two forms joined by RE chemistry (§5 gives the exact conditions).
   Thermodynamics carries over exactly: a fused or RE step's equilibrium constant is the product
   of the constants it replaces, and the Haldane and Wegscheider relations apply to it like any
   other step.
3. **Slice.** Keep the mechanism and add one dependent-parameter row per phantom, chosen so
   that every orbit of equivalent parameter sets crosses it exactly once. The row changes only
   which parameters are fitted, not `rate_equation`.

**When a quotient is exact.** Every test so far fits one rule (argued, not proved): turning an
SS group into RE keeps every rate law when every phantom orbit ends at that group's RE limit
while all other constants stay finite. That holds for EQ, BH, CHAIN, CHAIN-b and HYBRID-a. Where
an orbit ends depends on the parameter point for CUT, TRI, PAR-only and RING, and for PAIR the
orbit ends where a form's population goes to zero; those classes need a slice.

**What does not work.**

- **Pinning a single constant.** It is exact only when that constant sweeps (0, ∞) along every
  orbit (the kon scale of a zero-flux group, Kiso in PAIR). Elsewhere it loses rate laws: in
  Michaelis–Menten with SS binding, koff = kcat misses exactly the laws with V1/V2 ≤ 1.
- **Ties the grouping can express.** A tie either leaves the phantom in place or cuts off part
  of the family.
- **Free Cleland constants.** Reparameterizing into free kinetic constants over-covers: about a
  third of positive iso coefficient sets are not iso rate laws.
- **Trusting rank alone.** Several edits restore full rank and still lose rate laws (for
  example, making the chemistry port RE in HYBRID-a). Family preservation must be checked
  separately.

How each class maps onto the three kinds of fix:

| Kind | Classes where it is exact | Where it loses rate laws |
|---|---|---|
| Avoid: reject the child | EQ; BH, CHAIN and PAR+BH flips (47/47 tested flip pairs kept every rate law); JOINT, when both of its RE parents are also emitted | CHAIN-b, deep HYBRID, TRI and PAIR→CUT flips; JOINT when a parent is missing; many splits (33 of 56 tested same-rank split pairs reach rate laws their parent cannot) |
| Quotient: SS group → RE, fuse or merge | EQ, BH, CHAIN (entry RE or merge), PAR+BH, CHAIN-b (chemistry), HYBRID-a (entry), PAIR (merge) | Converting the wrong group: every class has a plausible conversion that fails; merging an SS isomerization between branch points |
| Slice: one dimensionless row | EQ (kon = 1), PAIR (Kiso = 1), CUT (crossing product), PAR-only (koff₁ = koff₂), CHAIN-b (k_chem² = koff_in·koff_out), TRI (k_f·k_r = koff_B·koff_Q), RING (product of two exit rates = product of the other two) | A row the orbit does not cross (it deletes a real direction), or any row applied where no phantom exists |

A slice that involves a Haldane-dependent constant carries a ±log Keq term in the package's
fitted coordinates; whether it does depends on which constant the constraint solve makes
dependent in that mechanism, not on the class. None of the slices can be written as a mechanism
edit: adding them means a new kind of dependent row in the thermodynamic constraint solve.

## 5. Fusing steps and merging central complexes

Denis proposed removing intermediates that have a single route in and a single route out, and
merging central complexes such as E(A,B) and E(P,Q). Two operations make this precise:

- **ELIM(X)** removes a form X with exactly two steps Y ⇌ X ⇌ Z and replaces them with one
  fused step Y ⇌ Z that carries both steps' free metabolites.
- **MERGE(X₁, X₂)** contracts a ligand-free step X₁ ⇌ X₂ (typically chemistry between central
  complexes) into a single form.

**Thermodynamics carries over exactly (proved).** A fused step's equilibrium constant is the
product of the two it replaces, and every cycle of the original maps to a cycle of the fused
mechanism with the same product of constants and the same net turnover. Merging a ligand-free
step keeps the cycle structure. Both results are thermodynamically consistent. Fusing two SS
steps, or merging an SS step, removes two fitted parameters; the RE cases remove one.

**When fusion keeps every rate law (proved).** Take a pass-through form X between two SS steps,
Y + L₁ → X → Z + R₂, where leaving X binds nothing. At steady state X's population is
α·L₁·[Y] + β·R₂·[Z], and the rest of the mechanism sees exactly the fused step. The original
mechanism is therefore the fused one plus two free "dead-end masses", α·L₁·[Y] and β·R₂·[Z],
and the constants map one-to-one. Fusion keeps the family exactly when both masses can be
absorbed elsewhere in the mechanism: when a mass carries no ligand, when an RE partner of Y has
the same content, or when another pass-through form is entered from Y by the same ligand (after
merging parallel duplicate steps and RE pairs). Denis's rule needs two refinements:

- The form must bind nothing when it is left; otherwise the fused step's rate depends on a
  concentration and is not mass action.
- Neither of its steps may belong to a context-shared kinetic group; otherwise the fused step
  needs fresh constants and the fitted count rises.

Checking only that the dropped monomials still appear elsewhere is not enough (50
counterexamples), nor is adding a check that rank is unchanged (20 ping-pong counterexamples).

**Two more exact rules (proved).**

- **Merging a ligand-free RE step is always exact** and removes exactly one phantom, provided no
  context-shared group touches the merged pair. This settles the PAIR class: merging the
  ping-pong RE half-reaction is now a proved fix, not an argued one.
- **The RE sandwich.** Y + L ⇌ X₁ (RE), X₁ <--> X₂ (SS), X₂ ⇌ Z + R (RE) has the same fitted
  count and the same family as Y + L <--> X <--> Z + R with a single merged complex and both
  steps SS. The enumerator's RE ordered seed (5 of 5) is therefore the same model as the merged
  central complex with SS second binding and SS first release.

**Denis's examples, computed exactly** (fitted / rank):

| Mechanism | Fitted / rank | Relation |
|---|---|---|
| Michaelis–Menten, all three steps SS | 5 / 3 | |
| Michaelis–Menten, merged central complex E + S ⇌ X ⇌ E + P (SS) | 3 / 3 | same family as the all-SS and the RE forms |
| Ordered bi-bi, all SS, E(A,B) and E(P,Q) separate | 9 / 7 | |
| Ordered bi-bi, merged central complex | 7 / 7 | same family, same 11-term denominator |
| Theorell–Chance, E(A) + B ⇌ E(Q) + P | 5 / 5 | different family: lacks the ABP and BPQ terms; the limit of the ordered family as the central complex empties |
| Ordered ter-bi, all SS | 11 / 9 | merged central complex 9 / 9, same family |
| Ping-pong, six steps, all SS | 11 / 6 | |
| Ping-pong, textbook four steps | 7 / 6 | same family as the six-step form; one phantom remains |
| Random bi-bi, all SS | 15 / 15 | |
| Random bi-bi, merged central complex | 13 / 13 | same 48-term denominator but a different, lower-dimensional family on the boundary (whether it lies inside the SS family is unproven) |

**Why the textbook ping-pong has 6, not 7, identifiable parameters.** All 8 denominator
coefficients, and the 8 Cleland constants built from them (Vf, Vb, K_mA, K_mB, K_mP, K_mQ,
K_iA, K_iQ), can be measured; Vmax is already one of them, because numerator and denominator
share a common scale. For every ping-pong mechanism the coefficients obey two identities (proved
symbolically; ρᵢ = k₋ᵢ/kᵢ):

1. c_A·c_B = Keq·c_P·c_Q, the ping-pong Haldane relation among kinetic constants.
2. Normalized by the numerator, c_A = x₃ + ρ₃x₄, c_B = x₁ + ρ₁x₂, c_AB = x₂ + x₄ and
   c_PQ/(ρ₂ρ₄) = ρ₃x₁ + ρ₁x₃ with xᵢ = 1/kᵢ; the fourth equals ρ₁·(first) + ρ₃·(second) −
   ρ₁ρ₃·(third).

So 8 − 2 = 6 independent numbers against 7 fitted constants. The invisible direction shifts
1/k₂ up and 1/k₄ down by the same amount (with k₁ and k₃ compensating), leaving
Vf = 1/(1/k₂ + 1/k₄) unchanged: the data fix the total time of the two product-release
half-reactions, not its split. No ELIM, MERGE or RE conversion removes this phantom. Two exact
options exist: the slice kon(A→X₁)·kon(Q→X₂) = kon(P→X₁)·kon(B→X₂), which gives 6 of 6, or
the union of three identifiable four-step forms selected by the sign of
σ = d_AP·d_B − d_AB·d_P (A binding RE when σ < 0, Q release RE when σ > 0, both when σ = 0);
a second, equivalent decomposition uses σ′ = d_BQ·d_A − d_AB·d_Q. Today the enumerator reaches
the textbook family at 7 fitted parameters (RE bindings with both chemistry steps SS, a CUT
mechanism), and the six-step all-SS form at 11.

**Denis's point 5: merge the central complex whenever a mechanism has more than one SS step.**
Not valid in general. Merging an RE isomerization is always exact, but merging an SS
isomerization keeps the family only in special cases: when an end is a pass-through whose other
mass can be absorbed (ordered sides), or through the RE sandwich. In exact scans of mechanisms
with per-step constants, 372 of 429 SS merges lost rate laws. On the complete bi-bi closure,
merging lowers the rank in 96% of the mechanisms it changes. With a guard that keeps rank, it
fixes 888 of 14,873 non-identifiable mechanisms (6.0%), almost all of them PAIR or CHAIN, and
removes 9.8% of phantom directions.

**Combined with canonicalization.** Applying the guarded merge after the RE-reduced core (§6)
makes 1,070 of the 3,002 cores that still carry phantoms fully identifiable (1,034 of the 1,134
PAIR cores), cutting the directions that need a slice from 3,144 to 2,036. For CHAIN, merging
the chain and making its entries RE are equivalent exact fixes with the same count; the merge
produces Cleland's single central complex, which the enumerator does not build today.

**Denis's point 9: reducing identifiable mechanisms.** Merging reduces 98.7% of the full-rank
closure mechanisms by at least two parameters (elimination reduces 18%), and almost all results
are identifiable. Because the original was identifiable, the reduced mechanism always has a
different family (proved by dimension). In every certified case it lies on the original's
boundary but outside it, as Theorell–Chance does for ordered bi-bi, and three quarters have no
counterpart among the reachable mechanisms. These are simpler nested models, not
reparameterizations: useful candidates for parsimony, with genuinely different rate laws.

**Denis's points 7 and 8: new initial mechanisms.** From the 69 bi-bi seeds, merging or
eliminating the central complex and choosing which other steps become SS gives 4,686 variants,
all thermodynamically consistent. 1,554 have no more parameters than their seed, and 1,244 of
those are identifiable, but most are degenerate: they lack the AB or PQ denominator term, or
their rate does not vanish as a substrate goes to zero. Specifically:

- **Duplicates.** Merging the central complex with SS entry and exit gives exactly the RE seed's
  family (proved), so these add nothing.
- **Fewer parameters.** A merged complex with a single SS step has 4 parameters and lacks the PQ
  or the AB term; an RE Theorell–Chance mechanism has 3. These are limits of the seed family,
  not seed-like models.
- **Same count, new, well-behaved.** For ordered bi-bi: the merged complex with A binding and P
  release SS, the merged complex with B binding and Q release SS, and the all-SS
  Theorell–Chance mechanism. All are identifiable 5-parameter families that the moves cannot
  reach today. They lie on the boundary of reachable families.
- **Random topology.** No variant with 5 or fewer parameters behaves like the seeds.
- **Ping-pong.** The 5-parameter merges with A binding and Q release SS (or B binding and P
  release SS) have the textbook 8-term denominator, but they are one-dimension-smaller faces of
  the textbook family, which the enumerator already reaches at 7 parameters. The merge with SS
  steps around the first merged form equals the family of the plain PAIR seed, making it an
  identifiable 5-parameter replacement for two of the four PAIR seeds.
- **A flaw in today's ping-pong seeds.** With no products present, every RE ping-pong seed gives
  v = V·A/(K + A), independent of B: the RE second half-reaction pulls all enzyme back to E
  instantly.

**What the package supports today (checked under the package's derivation against independent
mass-action rate laws).** Correct: fused chemistry-plus-release written in the release direction
(twelve test fixtures use it, including the Segel ping-pong), fused binding-plus-chemistry
written in the binding direction, and the merged ordered central complex. Broken:

- **Silently wrong results.**
  - A fused binding written in the reverse direction with an RE release returns −v.
  - The DSL keeps one metabolite of a Theorell–Chance step and silently drops the other
    (`dsl.jl:744-745`).
  - A fused binding and an ordinary binding from the same form get the same parameter names, so
    their constants are silently tied.
  - `_compute_numerator` (`rate_eq_derivation.jl:549-551`) leaves steps at mixed
    substrate/product complexes out of the reaction cut. Whenever such a step carries catalytic
    flux, the numerator is wrong without an error. This also affects unmerged hand-written
    mechanisms with RE chemistry and an SS abortive binding (package 1.249 against exact
    1.841). Today's enumerated mechanisms are unaffected, since each contains a single SS
    chemistry step.
- **Crashes.** An RE fused release raises a BoundsError in `_compute_alpha`
  (`rate_eq_derivation.jl:305`).
- **Representation.** `Step` stores one free metabolite and not which side it is on, so
  Theorell–Chance steps cannot be written, and the side of a fused step is guessed from its
  shape. Fused steps are stored as written, so the same step written in each direction gives two
  different mechanisms and names.
- **Enumerator.** The moves recognize chemistry only as a ligand-free isomerization:
  `_flux_carrying_groups` gives a parent with fused chemistry no flips at all, the atom check
  rejects fused releases, and the dead-end move places inhibitor bindings on the wrong form for
  a fused release.

These must be fixed before the enumerator emits any fused or merged mechanism; the first group
already affects hand-written mechanisms.

## 6. Enumerator-level policies

The sampled profile cannot measure reachability, so the policies were compared on complete
closures: every mechanism the catalytic moves reach from the 69 bi-bi seeds within the beam's
caps. Bi-bi up to 14 fitted parameters has 37,264 mechanisms, 14,873 of them (39.9%) with
phantoms; bi-bi with an inhibitor up to 10 parameters has 88,076, 18,740 (21.3%) with phantoms.
Every mechanism was ranked with the fast exact engine. The sampled population had overstated the
reachability losses of rank gates; the closures show none.

| Policy | Fits (bi-bi) | Still non-identifiable | Full-rank mechanisms cut off | Families lost (estimated) | New machinery |
|---|---|---|---|---|---|
| P0 today | 37,264 | 14,873 | 0 | 0 | none |
| P4 count rank on the beam axis | 37,264 | 14,873 | 0 | 0 | exact rank; does not meet the goal |
| P2 reject re_to_ss children whose rank does not rise | 30,359 | 7,968 | 0 | ≈ 0.7% | exact rank |
| P1 reject re_to_ss and split children whose rank does not rise | 25,785 | 3,394 | 0 | ≈ 5.5% | exact rank |
| G1 emit only full-rank mechanisms | 22,301 | 0 | 90 (0.4%) | ≈ 9.6% | exact rank |
| P5 expand as today, emit each mechanism's canonical identifiable form (RE-reduced core plus slices), deduplicated | ≤ 25,333 | 0 | 0 | ≈ 0 (see below) | SS→RE group conversion, slice rows, class detection, merge normalization, some canonical forms outside today's move space |

Findings behind the table:

- **A gate that requires rank to rise never cuts off a full-rank mechanism.** Every move adds
  at least one parameter, so a full-rank child always gains rank over its parent:
  rank(child) = n(child) > n(parent) ≥ rank(parent). It is lost only if all of its parents are
  rejected, which never happened in either closure. The stricter gate "rank rises by the full
  parameter increase" does cut off full-rank mechanisms (66 in bi-bi, 8 with the inhibitor), by
  removing their non-identifiable ancestors. What a gate loses is families: a rejected child's
  family often does not lie inside its parent's. Gates also leave mechanisms whose rank rose,
  but not to full.
- **Emitting only full-rank mechanisms meets the goal but loses about a tenth of the families,**
  mostly split children and slice-class phantoms whose rate laws no identifiable mechanism in
  the move space reaches.
- **Canonicalization meets the goal and keeps almost every family.** Replacing a mechanism with
  its core (zero-flux and absorbable transit groups made RE) kept every tested rate law (60 of 60
  new pairs, 33 of 33 earlier), and the remaining phantoms take a slice. The core has the same
  best achievable loss as the mechanism it replaces, because the mechanism's family lies inside
  the core's and the core is the mechanism's RE limit. The comparison first counted 82 bi-bi
  mechanisms (0.2%) in classes without a fix; the ring-gauge and PAR-variant slices (§3), found
  in parallel, restore full rank on all of them, but their exactness is proved only on five
  representatives.
- **Cores must be merge-normalized before deduplication.** 360 of the bi-bi cores (972 with
  the inhibitor) are ordinary enumerated mechanisms with one ligand's RE group written in two
  pieces whose constants are provably equal; merging them identifies the families. Without
  this step those mechanisms look like lost or new families, and P5's emitted count includes up
  to 360 duplicates. After merging, 872 bi-bi cores (1,324 with the inhibitor) remain outside
  the space the moves reach today; the `Mechanism` constructor accepted every one.
- **Cost.** Exact rank takes 21–125 ms per mechanism in pure Python (median 47 ms), and the core
  0.2 ms, against about 290 ms for `compile_mechanism` plus `fitted_params` and 710 ms for
  `rate_equation_string` in the package: about 4–10% of the work the beam already does per
  mechanism.

Family losses were estimated from small samples (110 bi-bi and 36 inhibitor children) drawn
with a skewed law sampler, testing only parents as covering mechanisms, and then corrected for
merge-equivalent cores; treat them as rough. The bi-bi-with-inhibitor losses were not
re-estimated after the correction.

## 7. Open questions

**Design choices for Denis.**

- **Which policy.** Requiring rank to rise on re_to_ss only is cheap and loses little (≈ 0.7% of
  families) but leaves about half of today's phantom mechanisms; gating splits too leaves about
  a quarter and loses ≈ 5.5%. Emitting only full-rank mechanisms meets the goal but loses about
  a tenth of the families. Canonicalization meets the
  goal with almost no loss but needs SS→RE group conversion, dependent slice rows in the
  constraint solve, class detection and merge normalization, and it emits some mechanisms the
  moves never build today.
- **Substrates and products declared as inhibitors.** Allowing their dead-end modes at forms
  where the same ligand binds catalytically is the main source of PGK's phantoms; forbidding it
  would also drop identifiable modes that enlarge the family.
- **The beam's complexity axis.** Until every emitted mechanism is identifiable, fitted counts
  overcharge phantom mechanisms in the beam, the parsimony filter and the LOOCV tie-break.

- **New initial mechanisms.** Which of the same-count variants (merged ordered complexes with
  A binding and P release SS, or B binding and Q release SS; all-SS Theorell–Chance) should seed
  the search, whether the degenerate fewer-parameter variants belong there too, and whether the
  RE ping-pong seeds (B-independent at zero products) should be replaced by merged forms.
- **Parsimony candidates.** Merging or eliminating steps of identifiable mechanisms yields
  simpler, different rate laws that the moves never reach; whether the enumerator should offer
  them as moves.
- **Step representation.** Supporting fused steps properly means storing which side each free
  metabolite is on (and allowing one on each side for Theorell–Chance), canonicalizing their
  direction, and giving them distinct parameter names; it touches `Step`, `_step_sides`,
  `_compute_alpha`, `_compute_numerator`, naming and the moves.
- **Package bugs.** The silently wrong results listed in §5 affect hand-written mechanisms today
  and are independent of any enumeration change.

**Proofs still missing.**

- The orbit-endpoint rule for when an SS→RE conversion is exact, in general.
- The HYBRID-a conversion (evidence only).
- Combined slices for mechanisms carrying CHAIN-b and PAR-only together (92 bi-bi cores, 12 of
  them needing three slices).
- The ring and PAR-variant slices beyond the five representatives where they are proved; on the
  closures they are known only to restore rank.
- Lossless rejection of the first dead-end inhibitor mode for the DUP-T, DUP-BH and DUP-PAR
  subclasses.
- A general slice rule for INH-DUP where rejection loses rate laws.
- Whether the merged random bi-bi family lies inside the SS random family or only on its
  boundary, and whether the absorbability conditions for fusion are necessary in general.
- Whether the package-checked CHAIN-b representatives (chemistry made RE) avoid the numerator
  bug, which bites when an SS step at a mixed substrate/product complex carries flux.
- True family-loss rates of the rank-guarded merge for BH+CUT and HYBRID-a mechanisms; the
  sampled estimates are upper bounds, since the fitter reported false losses in classes proved
  exact.

**Scope limits.**

- "Identifiable" here means locally identifiable (full rank). Discrete symmetries remain: for
  example, swapping the roles of chemistry and relaxation in an iso mechanism gives the same
  rate law, so a fit can be two-to-one even at full rank. Uniqueness on each slice was probed
  numerically only.
- The profile covers uni-uni, bi-bi and PGK. Not analyzed: ter-ter reactions;
  multi-conformation mechanisms such as iso and flip-flop (see
  `2026-09-23-iso-kinetic-cooperativity-findings.md`); mechanisms with three or more gauge
  blocks; and the package's default inhibitor-required subgraph. The bi-bi-with-inhibitor
  closure stops at 10 parameters because of memory.
- Family-inclusion evidence comes from a few sampled rate laws per test, and the Python sampler
  was skewed; loss fractions should be re-estimated with unbiased samples in fitted coordinates
  before they are relied on.

## Reproducing

The profile, catalog and verification scripts were throwaway session work and are not in the
repository. Everything above can be reproduced from this description: generate mechanisms with
the enumerator's catalytic moves, rank them with the BigFloat method in §1, and test fixes by
fitting the fixed mechanism to exact rates of the original.
