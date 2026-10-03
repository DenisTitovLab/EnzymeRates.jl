# Exact filters: zero-flux groups and duplicate inhibitor copies (sub-project B)

Date: 2026-09-30. Status: design agreed with Denis section by section on 2026-09-30. Part of
the A + B + C plan in `2026-09-28-identifiable-enumeration-findings-and-plan.md`, whose points 1
and 2 hold the evidence cited here; the track reports `t1_flux_report.md` and
`t2_inhdup_report.md` and their skeptic verdicts are in `2026-09-28-identifiability-evidence/`.
Built on branch `step-explicit-stoichiometry` after A, for the one pull request that carries
A, B and C.

## Why

Two classes of mechanism the enumerator emits carry parameters that steady-state data cannot
see, and both have exact structural tests.

- **A steady-state (SS) kinetic group with no flux-carrying step** adds exactly one phantom: the
  rate depends on its two constants only through their ratio, and making the group RE gives the
  same family with one parameter fewer (proved; exact on all 412 affected exported mechanisms).
  Today's screen, `_flux_carrying_groups`, is sound but runs on the parent's uncontracted graph
  before the flip, keyed on `is_iso`, and the split move never calls it. In bi-bi at depth 2 it
  lets through 148 flip children, each its parent plus two phantoms, and the split move adds 40
  children that isolate an abortive binding in its own SS group: 188 mechanisms and 350 of the
  616 phantoms. The screen also flags nothing in a mechanism with no isomerization step, so a
  Theorell–Chance or merged parent gets no flip children; C needs the correct test.
- **A dead-end copy of a substrate or product that creates no new complex** is a second, dead
  orientation of an existing complex. In Denis's uni-uni example the copy's constant enters only
  as K_S + K_S,inh, with kcat rescaled, and the child's family equals the parent's (proved). At
  depth 2 in R6, 3,082 of 36,844 dead-end events create a phantom, every one of them a copy whose
  complexes all duplicate existing forms. The split move can recreate the class by isolating a
  duplicate-only part of a copy's group (42 of 184 exported cases), and composition alone misses
  a complex that matches an existing form only through a rapid-equilibrium isomerization (3 cases
  at level 3).

## Goals

1. Every mechanism the moves emit has, in every SS kinetic group, a step that carries net
   steady-state flux, judged on the A-state catalytic graph.
2. Every mechanism the moves emit has no redundant copy group: every kinetic group that binds a
   competitive inhibitor has, in some conformational state where it binds, a complex that no
   productive form duplicates, or a dwell gauge that the mechanism's kinetic groups block
   (section 2).
3. Rule 1 loses no family: a flip set that fails it is extended, and a split part that fails it
   is reverted to rapid equilibrium, so the representative without the zero-flux phantom is
   emitted (it may still carry a phantom of another class, such as the chain class C merges).
4. Both tests are structural, exact, O(steps + forms), free of numerics, and read a step's
   metabolite lists, so they hold for the fused and Theorell–Chance steps C introduces.
5. A required regulator that no mechanism can bind under the rules fails loudly.

## Non-goals

- `:NonequalAI` SS groups whose I-state copy carries no flux (a dead-inactive I-state has no
  chemistry, so every SS group there sits at equilibrium and only kf_I/kr_I is visible). Denis
  (2026-09-30): the RE/SS type belongs to the step and both conformations share it, so an
  I-state-only conversion is not a mechanism. The class stays, recorded below.
- Substrate inhibition at a uni-uni product complex (E·P·S): the dead-end move excludes every
  site that carries all products, so uni-uni cannot express it. Predates B.
- Chemistry detection through `is_iso` in `_hyperbolic_catalysis`, `_expand_to_allosteric`,
  `_expand_change_allo_state` and `_assert_chemistry_is_iso` (C).
- MERGE and ELIM moves; rank on the beam's complexity axis (deferred in the findings).
- The `Mechanism` constructor does not enforce either rule. Both are emission policies of the
  moves: a hand-written mechanism with a zero-flux SS group or a redundant copy still derives.

## Decisions (Denis, 2026-09-30)

- A split part with no flux-carrying step is **reverted** to rapid equilibrium, not dropped, and
  the gain test runs on the reverted child. At depth 4, 72 dropped split children have reverted
  forms that no other move sequence reaches.
- The flux test runs on the **A-state projection** only (see Non-goals).
- The **plain copy rule**: a copy must create a new complex, judged by composition and by
  (RE segment, offsets). It drops the tie-only children, which are identifiable only because a
  shared group blocks the absorption; a non-productive orientation of an existing complex is not
  a distinct hypothesis.
- `seed_mechanisms` **errors** when no mechanism binds every required regulator.
- Controller rulings during execution, for Denis to confirm: the twin test is judged per
  conformational state (section 2, "Conformational states"); the mechanism-level invariant
  judges twins against copy-free forms while the dead-end move's placement rule judges against
  all forms (withdrawn by Denis on 2026-10-02, below); `_expand_change_allo_state` filters its
  children by the invariant; a reverted split part may keep a phantom of the chain class
  (Goal 3).
- Denis (2026-10-02): the rule's priority is to keep every identifiable mechanism and tolerate a
  few percent of non-identifiable ones; one reading of the twin test (productive complexes only,
  per conformation, free enzyme always present) replaces the two readings; the split's
  all-forms test of 2026-10-01 is withdrawn.
- Denis (2026-10-02): competition in the dead-end move is decided per site, not per name. A
  substrate declared as a competitive inhibitor is a copy that binds its own dead-end site, and
  the reactant's catalytic site and the copy's dead-end site are different sites by definition:
  a form carrying only the copy of M is a site for an inhibitor competing with M, and a copy does
  not count toward the capacity test (section 5).
- Denis (2026-10-02): option 3 everywhere — a copy group is rejected only when every complex
  duplicates a productive complex and the dwell gauge of Theorem 2 is consistent; Case 3
  placements return.
- Controller ruling (2026-10-02), adopting the gauge wave's fix: a twin is a productive complex
  in the copy complex's RE segment with its offsets (Theorem 2's hypothesis (a)); a composition
  match formed at steady state is not one, since four such rejections had full rank.
  `:EqualAI` groups, and an `:EqualAI` copy's K*, are tied across conformations.
- Review fixes (2026-10-02), for Denis to confirm: the copy's complexes carry one gauge factor s,
  never their twins' factors, as Theorem 2's gauge has it. Their twins' factors rejected two
  identifiable mechanisms with a second copy bound at a complex (7 fitted, rank 7). With one
  factor, an RE mirror between two complexes has ratio 1 and no longer blocks the gauge: 140 R6
  level-2 mechanisms are newly rejected, every one with a phantom. The copy's own group is
  checked like any other group.

## Terms

- **Flux-carrying step**: an SS step whose net steady-state flux is not identically zero over the
  parameters. **Zero-flux group**: an SS group none of whose steps carries flux. A zero-flux step
  inside a group with a flux-carrying step costs nothing (819 exported mechanisms, 743 of them
  identifiable), so the rule is group-level.
- **Copy**: the `CompetitiveInhibitor` species that stands for a substrate or product also
  declared as a competitive inhibitor; it binds with the reactant's concentration. **Copy step**:
  a pure binding whose metabolite is a `CompetitiveInhibitor`, a reactant's copy or a foreign
  inhibitor alike; its **site** is `from_species` and its **complex** `to_species`. **Copy
  group**: a kinetic group of copy steps. A foreign inhibitor's complexes are always new, so
  rule 2 constrains only reactant copies, and one predicate serves both.
- **Composition** of a form: the multiset of its bound metabolites' names, ignoring the copy flag,
  with its conformation and residual.
- **Offsets**: `_re_segment_extras` gives, for every form, how many more of each metabolite it
  carries than the lowest form of its RE segment; within a segment two forms with equal offsets
  have proportional weights.
- **Twin site**: a copy step whose complex lies in a productive form's RE segment with equal
  offsets, so the two weights are proportional. A **productive** form carries no competitive
  inhibitor.

## Design

### 1. The flux predicate

`_flux_carrying_steps(groups, reaction)` returns one `BitVector` per group, and
`_flux_carrying_groups(groups, reaction)` reduces each to "some step carries flux". A
convenience method takes a mechanism and passes `steps(m)` and `reaction(m)`; an
`AllostericMechanism` is measured on `_state_mechanism(am, :A)`, as `_re_segment_count` is.

1. RE segments and offsets come from `_re_segment_extras(groups)`, which `_bottomless_re_segment`
   already computes on the same groups.
2. Every SS step is an edge between its two forms' segments. Its weight is
   Σ_x ρ_x·(σ_e[x] + u_a[x] − u_b[x]), where σ_e[x] is the copies of x the step consumes minus the
   copies it releases, u is the offset vector, and ρ_x is +1 for a substrate name of the reaction,
   −1 for a product name, 0 otherwise. Around any cycle of the segment graph the offsets
   telescope and the uptakes sum to n·ρ for n net turnovers, so the weights sum to n·Σρ²: zero
   exactly when the cycle is balanced. (The track's "first substrate" component gives the same
   verdict; this form needs no choice of metabolite.)
3. A self-loop, an SS step whose two forms share a segment, carries flux iff its weight is
   nonzero. Every other edge carries flux iff its biconnected block (`_edge_blocks`, which
   leaves self-loops unlabelled) holds an unbalanced cycle: one BFS spanning tree per block gives
   each vertex a potential, and the block is unbalanced when some edge's weight is not the
   potential difference of its ends. By the theta argument every edge of such a block lies on an
   unbalanced cycle.

The zero-flux direction is proved for any parameters and grouping; the flux-carrying direction
for one-way steps, which every flux-carrying step of the export is. The test matches exact flux
on all 25,850 exported SS steps and equals today's screen on all 67,235 exported groups. It
reads the metabolite lists, so a Theorell–Chance step EA + B → EQ + P is an unbalanced edge like
any other.

### 2. The copy predicate: all twins and a consistent gauge

A copy group is **redundant** when, in every conformational state where the copy binds, (i)
every complex it forms has a productive twin, a complex whose weight is proportional to its own
(the group is all-twin), and (ii) the dwell gauge of track 2's Theorem 2 is consistent, over the
states together (the test below). The rate law then has a direction, Lemma 1's dwell gauge,
along which the copy's constant moves and nothing observable changes: the copy adds a phantom.
Within Theorem 2's scope the gauge is a finite merge: the copy's constant enters the rate only
as K_g + K* beside the existing binding's, and the child's family is its parent's. No move
emits a redundant copy group (`_redundant_copy_groups`). A group that is all-twin but whose
gauge fails is kept: Case 3 (the shared A group that pins K_A at E(Q)), the other shared-group
cases and H1 are identifiable or unproven, and the rule's priority is never to drop an
identifiable mechanism. The gauge is the proof.

**Twins.** A copy step is a twin site when its complex has a productive twin: a form that
carries no competitive inhibitor, lies in the complex's rapid-equilibrium segment and has the
complex's offsets (`_re_segment_extras`), so that rapid equilibrium makes the two weights
proportional (`_productive_twin`, which returns the twin, preferring the complex's composition
when several forms qualify). This is Theorem 2's hypothesis (a), and the merge below is exact
only under it. The dead-end move evaluates it for a candidate site before the child exists: the
candidate complex lies in the site's segment and has the site's offsets plus one of the copy's
name, since the copy step is RE. The offsets reach the ping-pong copy that matches a form across
EB_res ⇌ EQ, where no form shares the complex's composition. A form of the complex's
composition in another segment, formed by a steady-state binding, is not a twin: its weight is
not a constant multiple of the complex's, so the merge is no reparameterization. Taken as a
twin it gave four identifiable rejections at depth 2 of R6 (an SS A binding with A* at E and a
random-order product side: 7 fitted, rank 7). Only productive forms are twin sources: Theorem 2
concerns productive twins, and a complex that duplicates only another copy's complex shares its
weight with that copy's constant, which the sites that pin the other copy may keep separable.

#### The gauge test

The test is track 2's Lemma 1 with Theorem 2's factors (`_gauge_rescaling`). Give every form a
factor σ: ρ_c on a twin of class c, one factor s on every complex of the copy, and 1 on every
other form. Multiplying an RE group's K by σ(from)/σ(to), and an SS group's kf by σ(from) and its
kr by σ(to), divides each form's weight by its factor and leaves every flux as it was. A twin's
weight is a constant multiple of its complexes', so the factors can keep each twin's total with
its complexes, and with it the rate law, unchanged to first order while K* moves. A kinetic
group shares its constants, so the rescaling must be one for all its steps. The gauge exists iff
that holds for every group, the copy's own included, with the factors treated as generic numbers
that are equal only where the structure forces them equal. Within Theorem 2's scope (every
complex a twin formed from its site by an RE binding, no step at a complex but the copy's own
and RE mirrors) the gauge is the merge of each copy complex X_i·L* into its twin T_i: T_i's
weight grows by the factor 1/ρ_i, with ρ_i = w(T_i)/(w(T_i) + w(X_i·L*)), every constant of a
step that leaves T_i is multiplied by ρ_i, every constant of a step that enters T_i is divided by
it, and the twin's binding group takes K_g' = K_g + K*.

- **Scale factors.** σ(F) = ρ_class(F) for a twin F; σ = s on every complex of the copy, one
  factor for all of them; σ(F) = 1 for every other form. The copy's own group, at 1 on each
  site, holds all its complexes at one factor. A complex never takes its twin's factor: its
  weight moves opposite to the twin's (in the merge its share of the twin's new weight is
  K*/(K_g + K*) = 1 − ρ_g). Labelled with its twin's factor, a second copy bound at a complex
  and at a twin (P* at E(A*) and at E(A, Q) on B14) read one ratio where there are two, and two
  identifiable mechanisms (7 fitted, rank 7) were rejected.
- **Classes.** Two twins share a class iff they are formed from their sites X_i and X_j by RE
  bindings of L in the same kinetic group g, each the twin of one complex (then
  ρ = K_g/(K_g + K*) for both). Every other twin has a class of its own: a twin reached from its
  site another way (across an RE isomerization), one formed by an L binding in another group, a
  twin of two complexes.
- **Consistency**, per kinetic group, the copy's own included:
  - an RE group needs one ratio σ(from)/σ(to) over its steps (its K becomes K·σ(from)/σ(to)).
    The L-binding group g of a twin has steps (1, ρ_g), and its K becomes K_g/ρ_g = K_g + K*,
    Theorem 2's K_g'. A step of g binding L at a form that is not a site gives (1, 1) and breaks
    the group: Case 3.
  - an SS group needs one σ(from) and one σ(to) (kf scales by σ(from), kr by σ(to)).
  - a mirror step joins two complexes, both at s. An RE mirror has ratio 1, as its parent step
    between two sites has, so it does not block the gauge whatever the twins' classes. An SS
    mirror carries flux between the complexes and needs s beside its parent step's 1: never
    consistent.
- **Bookkeeping.** σ is a class id, 0 for the factor 1 and the copy's own group's id for s. An
  RE step's key is (σ(from), σ(to)), or (0, 0) when the two are equal; an SS step's key is the
  pair itself; a group is consistent when its steps share one key. No parameters, no random
  points.

The cases, ranked with the test oracle (fitted/rank):

| Case | Mechanism | Gauge | Verdict | Fitted/rank |
|---|---|---|---|---|
| B1 | ordered bi-bi, A* at E | E(A) is formed by the A group alone | redundant | 6/5 |
| Case 3 | A group {E + A, E(Q) + A} with abortive E(A, Q), A* at E | A group: 1/ρ at E, 1 at E(Q) | kept | 6/6 |
| B14 | Case 3's parent, A* at {E, E(Q)}, E + Q ⇌ E(Q) mirrored | both twins share the A group's ρ; the Q group has ratio 1 throughout | redundant | 6/5 |
| H1 | A* at {E, E(B)}, twins from two A groups, shared B group | B group: ρ₁/ρ₂ beside 1 | kept | 7/7 |
| Case 1 | uni-uni, S* at E | S group alone | redundant | 4/3 |
| ping-pong | RE second chemistry, abortive E(B, P; res), P* at E, Q* at E(P*) | offsets twin E(B, P; res), touched by one single-step group | redundant | 9/7 (parent 8/7) |
| allosteric | S and chemistry `:OnlyA`, `:EqualAI` S* at E | new in the inactive state | kept | 5/5 |
| SS binding | A binds E at steady state, A* at E, random-order product side | E(A) is in another segment: no twin | kept | 7/7 (parent 6/6) |
| tie | A binding `:NonequalAI`, B binding `:EqualAI`, chemistry `:NonequalAI`, A* at E | each state consistent; the B group needs K_B·ρ_A in one, K_B·ρ_I in the other | kept | rank one above the parent's |
| two classes | A in two groups, A* at {E, E(B)}, E + B ⇌ E(B) mirrored | the RE mirror between the complexes has ratio 1 | redundant | 7/6 |
| B14 split | B14 with the A bindings and E(A) + Q ⇌ E(A, Q) in groups of their own | twins of two classes; the RE Q mirror has ratio 1 | redundant | 7/6 |
| SS mirror | the same with the Q group at steady state | the mirror needs s beside 1 | kept | 8/6 |
| P* at a complex and a twin | B14, P* at {E(A*), E(A, Q)} | P* group: s beside ρ | kept | 7/7 |
| P* at a complex and its twin | B14, P* at {E(A*), E(A)} | P* group: s beside ρ | kept | 7/6 |
| P* at two complexes | B14 split, P* at {E(A*), E(A*, Q)} | both P* steps leave a complex at s | A* redundant | 8/7 |

**Conformational states.** A copy binds the active state always and the inactive state unless
its tag is `:OnlyA`. The test runs on the graph of each state where the copy binds; the inactive
graph is `_state_mechanism(am, :I)`, which prunes `:OnlyA` groups and the forms they strand. A
site that graph lacks binds nothing there, and the free enzyme is always present, since every
conformation holds it. A copy is redundant iff it is all-twin in every state where it binds and
the gauge is consistent over the states together. In a state where the copy binds nothing every
factor is 1. The states form one system: a kinetic group tagged `:EqualAI` has one set of
constants in both, so it must take the same rescaling in each; a shared class is one number in
both states only when its binding group and the copy are both `:EqualAI` (one K_g, one K*), the
copy's factor s is one number when the copy is `:EqualAI` (one K*), and every other class is a
number per state. So an `:OnlyA` copy whose twin is left by an `:EqualAI` group is not redundant
(the group is rescaled in the active state only), nor is a copy whose twin is formed by a
`:NonequalAI` binding and left by an `:EqualAI` one. An `:EqualAI` copy of S at E, in a
mechanism whose S binding is `:OnlyA`, duplicates E(S) in the active state but is the only
S-bound form in the inactive one; its constant is visible through that state's S-dependence
(rank rises by one), and the placement stands. This is the route by which substrate inhibition
enters an allosteric mechanism. With S binding `:NonequalAI` the copy is redundant in both
states and is skipped.

**Which moves can break the rule.** A flip only cuts segments and turns RE groups to SS: twins
can disappear but never appear (a twin formed by a flipped binding leaves its site's segment),
and an SS group's condition implies an RE group's, so a flip cannot complete a gauge whose
twins it keeps. A dead-end addition adds only forms that carry the new copy, none of them
productive, and steps to the groups an older copy's gauge reads, so older groups keep their
status, and the move tests the new copy's group alone (section 5). A split
can complete a gauge, by separating a binding that forms or leaves a twin from the bindings
elsewhere in its group that blocked it (Theorem 3's split sibling), so the split filters its
children as well as its units (section 4). `_expand_change_allo_state` changes the tags the
per-state test reads (relaxing an `:OnlyA` binding to `:NonequalAI` brings its complex into the
inactive state), so its children are filtered; the other allosteric moves change neither steps
nor the tags of existing groups. The parent assertion of section 6 covers hand-written input.

Edges left as they are (for Denis to decide):

- Two copies that are twins only of each other (Q* only at E(A) and A* only at E(Q): 7 fitted,
  rank 6; after two splits, 9 fitted, rank 8) are emitted, the tolerated family. The rule has no
  proof against them, and rejecting them would also reject identifiable splits whose other copy
  is pinned. The same holds for a second copy whose twin's group also binds at the first copy's
  form through a mirror (Q* at E beside A* at {E, E(Q)}, A* at E beside Q* at {E, E(A)}): the
  mirror's ratio 1 sits beside the twin step's, the gauge fails, and the placement is emitted
  with one phantom (7 fitted, rank 6). The regression record measures the fraction at level 3.
- A `:NonequalAI` copy group has one constant per state, and the test runs on each state's graph
  with that state's constant. A copy redundant in the active state alone, new in the inactive
  one, keeps one phantom (its active constant) and is emitted; a "redundant in either state"
  rule would reject identifiable groups as well (521 of 17,483 sampled allosteric mechanisms hold
  an all-twin-in-one-state group, with mixed ranks). The same principle as the I-state class
  under Non-goals: a constant that one conformation cannot see is not a reason to change the
  group's tag.
- A copy whose only match is a complex of the same composition formed at steady state is kept
  even when it is a phantom (U4, B10 in the track: the uni-uni and ordered copies at E with S or
  A bound at steady state). No exact merge exists there, B10's family is wider than its
  parent's, and four such copies at depth 2 of R6 have full rank.
- A second copy bound at a complex of the first and at that complex's own twin (P* at E(A*)
  and E(A) on B14) pairs E(A*, P*) with E(A, P*) and keeps a phantom (7 fitted, rank 6) that
  the test does not reach: it gives the complex its factor s and the twin ρ, so the P* group
  reads two ratios. E(A, P*) is no productive twin, and the test rescales only the copy's own
  complexes. Emitted, the tolerated family.
- A mechanism reachable only through a mechanism with a redundant copy group (a split of a
  redundant copy group; two copies each redundant alone) is not enumerated, though the rule on
  its own groups admits it; the regression record ranks all of them.

### 3. The flip move (`_expand_re_to_ss`)

- **Units**: all-RE groups binding no regulator, as today, that carry flux in the all-SS copy of
  the parent's groups. A flux-carrying step stays flux-carrying in every flip descendant (the
  lift lemma), so a group that is zero-flux with every step SS can never flip usefully; the
  all-SS graph is the exact pre-filter. This replaces the `is_iso` block test and gives
  Theorell–Chance and merged parents their flip children.
- **Gain test**: `gains(sel)` requires, in addition to the segment-count rise and the bottomless
  check, that every flipped group carries flux in `flipped_groups(sel)`. A zero-flux set counts
  as failed, so `_minimal_gaining_sets` extends it to its minimal flux-carrying supersets. The
  skeptic showed "extend" is required: 8 valid pair flips of R5 level-2 parents are reached no
  other way. Only the flipped groups are tested; the lift lemma keeps every other group's
  status, given a parent that satisfies rule 1 (section 6).
- Cost: one `_re_segment_extras` and one block decomposition per candidate, the same order as
  the two probes already in `gains`; the plan may share the extras with the bottomless probe.
- `_hyperbolic_catalysis` filtering of emitted children is unchanged.

### 4. The split move (`_expand_split_kinetic_group`)

- Per-step flux flags are computed once per parent: a split moves no edge, and reverting a
  zero-flux part contracts edges of a balanced block, which changes no cycle's net
  stoichiometry, so the parent's flux flags hold for every candidate. Twin status is judged on
  the candidate's parts against productive complexes (section 2).
- **Revert**: for a unit of an SS group, a part with no flux-carrying step is rebuilt with every
  step at rapid equilibrium. A parent that satisfies rule 1 has a flux-carrying step in one part,
  so at most one part reverts. The reverted child has the same family as the raw split child with
  one parameter fewer (point 1's conversion theorem).
- **Gain test on the reverted child.** `_partition_independent_count`'s counter takes a per-step
  kind alongside the group labels, so a reverted step is counted as a binding or isomerization
  `K` on the same cycle basis; the allosteric path builds the child and needs no change. A
  reverted child whose new constant the Wegscheider ties pull back to the parent's count is
  absorbed like any tied split, so no reparameterization of the parent is emitted. A candidate
  whose reverted groups leave a bottomless RE segment counts as failed, as in the flip; the
  construction of a rejected child is never attempted.
- **Copy groups**: a bipartition of a copy group in which a part would be redundant as a group of
  its own in the child (all-twin, `_twin_only`, checked first on the parent; then the gauge on
  the child that splits only that group) is not a unit. Every superset of such a unit keeps the
  part, and splitting other groups only relaxes the part's gauge. Excluding the unit loses only
  children whose part is provably a phantom.
- **Children**: a split of any group can complete a kept copy's gauge, by separating a binding
  that forms or leaves a twin from the bindings elsewhere in its group that blocked it (Case 3
  with its A and Q groups split is redundant: 7 fitted, rank 6). A child with a redundant copy
  group is not emitted; by Theorem 2 its family is that of the same split of the mechanism
  without the copy.
- The partner search is unchanged. It ignores groups holding an SS step, so a reverted RE part
  whose constant a further split could free is not extended within one move; the same child is
  reachable by splitting the partner first.
- A reverted child that equals an existing mechanism (all 40 depth-2 cases equal a level-2 flip
  child) is removed by the existing structural deduplication.

### 5. The dead-end move (`_expand_add_dead_end_regulator_native`)

Competition is decided per site (Denis, 2026-10-02): a copy occupies its own dead-end site, so a
form carrying only the copy of M is a legitimate site for an inhibitor competing with M, and
copies do not count toward the all-substrates or all-products capacity test. The move's target
sites for competition with a reactant are the forms where it binds productively.

After the pattern's `active` sites are chosen, a pattern whose every site duplicates a
productive complex (a twin, section 2) is a candidate for rejection, and it is skipped when the
copy's dwell gauge is consistent in the child (the candidate's sites, their twins and the
mirrors the move adds; section 2). On a parent that obeys the rule only the new copy's group
can be redundant. For a `Mechanism` the gauge runs on the candidate's groups with the parent's
twins before the child is built, since the copy's forms are not productive and join no two
segments, so every parent form keeps its segment and offsets; an allosteric child is built and
tested whole. A placement whose every site is a twin but whose gauge fails,
Case 3 among them, is emitted. The same kernel builds the required-regulator seeds in
`seed_mechanisms`, so required copies obey the rule. Effect of the all-twin skip at depth 2 in
R6: 10,552 of 36,844 events are all-twin, every phantom-creating one among them, with 7,470
identifiable tie-only children; the gauge returns those whose gauge fails, and the regression
record measures the populations. Every bi-bi seed keeps at least two children per copied
ligand.

### 6. `expand_mechanisms` and `seed_mechanisms`

- `expand_mechanisms` asserts both rules on each parent, next to `_assert_chemistry_is_iso`, with
  errors that name the offending group's steps (in canonical group order, like the group-rule
  errors). This is what lets the flip test only its flipped groups, and it makes a move that
  emits a violator, such as C's merge, fail at the next expansion rather than silently
  propagate. Three allosteric moves copy or pass through `steps` and leave the tags of existing
  groups, so they preserve both rules; the split and `_expand_change_allo_state` filter their
  children (section 2).
- `seed_mechanisms` errors when it finds no seed. The message names the required competitive
  inhibitors and allosteric regulators it could not place and the keywords that make a regulator
  optional (`optional_competitive_inhibitors`, `optional_allosteric_regulators`). Today the beam
  returns an empty result with no message. The first case this catches: a uni-uni reaction with
  S or P declared as a competitive inhibitor, whose only admissible site is free E and whose copy
  complex duplicates ES or EP, formed by a binding group that nothing else shares, so the gauge
  is consistent. Whether C's merged seeds give that copy a placement is C's decision.

### 7. Unchanged

`init_mechanisms` (its abortive complexes always have new compositions and add no phantom at
level 0), `_assert_mechanism_invariants`, the `Mechanism` and `AllostericMechanism`
constructors, every derivation function, `_minimal_gaining_sets`, `_context_bipartitions`,
`_edge_blocks`, `_re_segment_extras`.

## Test changes

Every new test follows the CLAUDE.md enumeration-test rules: fixtures inline in the testset,
exact child sets and counts, `_testhelper_` prefix for file-level helpers.

New tests:

- The flux predicate: today's six cases keep their flags; Theorell–Chance (every step carries
  flux) and a merged central complex (E + A ⇌ EA, EA + B → X, EQ + P → X, E + Q ⇌ EQ; every
  step carries flux), both of which today's screen flags empty; a child with an RE shunt, the
  E1 pattern (two context-split binding groups flipped together isolate free E, and the
  abortive route carries the turnover), where today's screen says flux and the new test says
  zero; an unbalanced self-loop (chemistry inside a segment) and a balanced one.
- The flip on the E1 parent: the exact child set, in which the zero-flux pair is extended.
- The flip on a Theorell–Chance parent: two children.
- The split on the E2 parent (mixed SS group {EA + B → EAB, EQ + B → EBQ}): the exact child set
  with EQ + B ⇌ EBQ reverted. The plan searches small parents for one whose reverted child the
  ties absorb and, when it finds one, pins that child's absence.
- The dead-end move: Case 1 (uni-uni + S: no child), Case 2 (ordered + B at E: kept), Case 3
  (shared A group with abortive EAQ, A at E: kept, its gauge fails), a copy at {E, EQ} kept
  because EQ·A is new, and the ping-pong (segment, offsets) case.
- The split that would isolate EA → EAB·B\* from its copy group: the exact child set without it.
- The gauge's factors (`_redundant_copy_groups`): a second copy bound at a complex and at a
  twin (kept, 7 fitted, rank 7, a valid parent), at a complex and its own twin (kept, 7 and 6),
  at two complexes whose twins have two classes (redundant, 8 and 7); an RE mirror between
  complexes of two classes, in two mechanisms (redundant, 7 and 6 each); an SS mirror (kept).
- `seed_mechanisms`: uni-uni with S required errors with the regulator's name in the message;
  bi-bi with A required returns seeds.
- A population test: R4 through `expand_mechanisms` to depth 2 (about half a second), asserting
  both rules on every child and pinning 62, 369 and 1,200 distinct mechanisms per level under
  the convention that gives 62, 369, 1,388 today.
- Rank-oracle spot checks (`_testhelper_identifiable_rank`) on the hand cases: each rejected
  child has its parent's rank; each kept child has full rank.

Existing tests whose expectations change (Denis, 2026-09-30: rewrite, do not delete):

- Four fixtures built from a uni-uni with S declared as its own inhibitor, whose only dead-end
  child is Case 1: "Substrate-as-dead-end-inhibitor overlap" under `_expand_re_to_ss`
  (`test/test_mechanism_enumeration.jl:1870`), its allosteric counterpart (1922),
  "keeps inhibitor bindings RE" (1963), and the dead-end overlap test (2655), which asserts the
  phantom child and its +1 parameter. Each keeps its intent on a fixture with an admissible
  placement (a bi-bi copy that creates a new complex), and the dead-end test asserts zero
  children for the uni-uni case.
- The `_hyperbolic_catalysis` tests with A as its own inhibitor (5748–5766): placements whose
  every site is a twin disappear; the child counts and denominator degrees are re-derived by
  hand from the rule.
- Every exact-children flip or split test whose parent is a split child is re-derived by hand
  from the rules; the rank oracle checks each changed expectation.
- The six `_flux_carrying_groups` tests call the new signature.

## Verification

Every change is made test-first.

1. **Regression record**, like A's. A scratch script enumerates, on a detached worktree at
   a2a02b1 and on the branch tip, `init_mechanisms` plus two levels of the three catalytic moves
   on R1–R6 (R6 level 2 in full or as A's fixed-seed sample) and the allosteric sample of A's
   record, keyed by steps with explicit lists. It asserts: every mechanism removed violates rule
   1 or rule 2 under the new predicates; every kept and every new mechanism satisfies both; each
   removed zero-flux mechanism's RE-converted twin is present; the rank oracle on every removed
   R4 zero-flux mechanism gives fitted minus rank at least its zero-flux group count; derivation
   strings are unchanged on a random sample of kept mechanisms (B touches no derivation code).
   The record reports per-level counts for every reaction and the level-3 phantom fraction.
   Added mechanisms are reported with their phantom fraction. Two oracle checks on R6 levels
   1–2 (`MersenneTwister(20260930)`): of 100 kept mechanisms with an all-twin copy group, the
   record reports how many have full rank (the test is conservative, so some may not); every
   one of 100 mechanisms removed for a redundant copy group must have fitted − rank ≥ 1, and a
   single exception makes the gauge test unsound.
2. **Timing**: enumeration of R6 levels 0–2 on both commits; the design's limit is a 20% rise.
   The ter-ter worst-case split timing test in the suite stays green.
3. **Full suite** before every commit; the branch's version rises once over `main`, in A
   (0.8.0); B adds no bump.

## Risks

- **Search cost in the flip.** A zero-flux set is extended rather than dropped, so Apriori visits
  its supersets. The all-SS pre-filter removes units that can never gain, and the timing check
  bounds the rest.
- **Hand-derived expectations.** Rewritten exact-children tests are derived from the rules by
  hand; each is cross-checked with the rank oracle and, for the R4 population, with the track's
  counts.
- **Fewer seeds.** A required copy of a substrate or product now needs a new complex; bi-bi keeps
  at least two children per seed and ligand, uni-uni has none and errors. C's merged seeds may
  change the uni-uni answer.

## Documentation

- `docs/src/developer.md`: the paragraph on the two refinement moves (flux rule, revert,
  copy rule).
- `docs/src/identify/enumeration_engine.md`: the flip's "two kinds of group never flip"; the
  split section (revert); the dead-end section (new-complex rule); "Dead-end branches stay at
  equilibrium" restated as the flux rule; the sentence that a substrate declared as a dead-end
  inhibitor is how substrate inhibition enters an allosteric mechanism, qualified by the rule;
  the "never a provable copy of its parent" list gains the zero-flux set and the duplicate copy.
- The findings document's section "B. Exact filters" points here and records the four decisions.
