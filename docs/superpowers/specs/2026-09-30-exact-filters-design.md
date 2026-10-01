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
2. Every mechanism the moves emit has, in every kinetic group that binds a competitive inhibitor,
   a binding that creates a complex no other form duplicates.
3. Rule 1 loses no family: a flip set that fails it is extended, and a split part that fails it
   is reverted to rapid equilibrium, so the representative without the zero-flux phantom is
   emitted (it may still carry a phantom of another class, such as the chain class C merges).
4. Both tests are structural, exact, O(steps + forms), free of numerics, and read a step's
   metabolite lists, so they hold for the fused and Theorell–Chance steps C introduces.
5. A required regulator that no mechanism can bind under the rules fails loudly.

## Non-goals

- The gauge test that would keep the identifiable "tie-only" copies (section 2).
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
  moves: a hand-written mechanism with a zero-flux SS group or a duplicate copy still derives.

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
  all forms (section 2, "Two readings"); `_expand_change_allo_state` filters its children by
  the invariant; a reverted split part may keep a phantom of the chain class (Goal 3).

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
- **Twin site**: a copy step whose complex has another form's composition, or lies in another
  form's RE segment with equal offsets.

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

### 2. The twin predicate

A copy step is a twin site when its complex duplicates another form of the mechanism under
either key: the composition key, or the (segment, offsets) key with the complex's segment and
offsets read from `_re_segment_extras`. The dead-end move evaluates it for a candidate site
before the child exists: the candidate complex has the site's composition plus the copy's name,
lies in the site's segment, and has the site's offsets plus one of the copy's name, since the
copy step is RE. Both keys are needed: composition alone misses the ping-pong copy that matches
a form across EB_res ⇌ EQ, and offsets alone miss a copy of a metabolite bound at steady state
(48 cases at depth 2). Their union catches every phantom-creating dead-end event at depths ≤ 2
and 3 of R6, for 14 identifiable events beyond the composition key alone.

A copy group satisfies rule 2 when one of its steps is not a twin site.

Edges left as they are (found during the final review, 2026-09-30; for Denis to decide):

- A `:NonequalAI` copy group has one constant per state, and the rule asks only that the copy
  be a twin in every state where it binds. A copy that is a twin in the active state alone
  keeps one phantom (its active constant) and is emitted; a "twin in either state" rule would
  reject identifiable groups as well (521 of 17,483 sampled allosteric mechanisms hold such a
  group, with mixed ranks). The same principle as the I-state class under Non-goals: a
  constant that one conformation cannot see is not a reason to change the group's tag.
- A copy placed at free E in a parent whose inactive graph has no steps (every catalytic group
  `:OnlyA`) is judged absent from the inactive state and rejected at placement when its
  complex duplicates E(S) in the active one. Counting free E as always present would emit 107
  such placements in the sampled allosteric population, 62 of them identifiable and 45 with a
  phantom. The conservative verdict stands until Denis decides.
- The mechanism-level invariant passes a hand-written mechanism whose two copies are twins only
  of each other (Q* only at E(A) and A* only at E(Q): 7 fitted, rank 6). No move builds it: the
  dead-end placement and the split's part filter both judge against all forms.

Conformational states (added 2026-09-30 during implementation, for Denis to confirm): a copy
in an allosteric mechanism is a twin only when it duplicates a form in every state where it
binds. A copy binds the active state always and the inactive state unless its tag is `:OnlyA`;
a site that the inactive state's graph lacks (`_state_mechanism(am, :I)` prunes `:OnlyA` groups
and the forms they strand) binds nothing there. So an `:EqualAI` copy of S at E, in a mechanism
whose S binding is `:OnlyA`, duplicates E(S) in the active state but is the only S-bound form
in the inactive one; its constant is visible through that state's S-dependence (rank rises by
one), and the placement stands. This is the route by which substrate inhibition enters an
allosteric mechanism. With S binding `:NonequalAI` the copy duplicates E(S) in both states and
is skipped.

Two readings of the rule (settled 2026-09-30 during implementation, for Denis to confirm). The
dead-end move judges a new copy's sites against every form the parent has, other copies'
complexes included: that is the per-placement rule the findings measured (no phantom among the
26,292 placements that create a new complex). The split's part filter uses the same all-forms
reading: once two copies are each split down to complexes of one composition, their two
constants enter the law through one coefficient (Denis, 2026-10-01). The mechanism-level
invariant that `_expand_change_allo_state` and the parent assertion enforce judges a copy
group's sites against copy-free forms only, the forms bound to no competitive inhibitor. The two
differ when a later copy's complex has the composition of an older copy's, as E(Q, A*) has that
of E(A, Q*): the older group is then all-twin by composition, yet the two constants enter the
law through K_A·K_Q* + K_Q·K_A*, and the sites that pin the other copy keep them separable, so
such a mechanism is generically identifiable and is not rejected.

Which moves can break rule 2: a flip only cuts segments and changes no composition, so a twin
can disappear but never appear; a dead-end addition adds only forms that carry the new copy,
none of them copy-free, so older groups keep their status under the invariant;
`_expand_change_allo_state` changes the tags the per-state test reads (relaxing an `:OnlyA`
binding to `:NonequalAI` brings its complex into the inactive state and can make a kept copy a
twin in both states), so its children are filtered by the invariant; the other allosteric moves
change neither steps nor the tags of existing groups. The split move, the dead-end move and
`_expand_change_allo_state` enforce the rule; the parent assertion of section 6 covers
hand-written input.

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
  the candidate's parts against every form (section 2, "Two readings").
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
- **Copy groups**: a bipartition of a copy group in which a part has only twin sites
  (`_copy_twin_test` against every form, in every conformational state where the copy binds) is
  not a unit. Every superset of such a unit recreates the same duplicate-only group, so excluding
  the unit loses nothing.
- The partner search is unchanged. It ignores groups holding an SS step, so a reverted RE part
  whose constant a further split could free is not extended within one move; the same child is
  reachable by splitting the partner first.
- A reverted child that equals an existing mechanism (all 40 depth-2 cases equal a level-2 flip
  child) is removed by the existing structural deduplication.

### 5. The dead-end move (`_expand_add_dead_end_regulator_native`)

After the pattern's `active` sites are chosen, a pattern whose every site is a twin (section 2)
is skipped. The same kernel builds the required-regulator seeds in `seed_mechanisms`, so
required copies obey the rule. Effect at depth 2 in R6: 10,552 of 36,844 events removed, every
phantom-creating one among them, 7,470 identifiable tie-only children with them; every bi-bi
seed keeps at least two children per copied ligand.

### 6. `expand_mechanisms` and `seed_mechanisms`

- `expand_mechanisms` asserts both rules on each parent, next to `_assert_chemistry_is_iso`, with
  errors that name the offending group's steps (in canonical group order, like the group-rule
  errors). This is what lets the flip test only its flipped groups, and it makes a move that
  emits a violator, such as C's merge, fail at the next expansion rather than silently
  propagate. Three allosteric moves copy or pass through `steps` and leave the tags of existing
  groups, so they preserve both rules; `_expand_change_allo_state` filters its children
  (section 2).
- `seed_mechanisms` errors when it finds no seed. The message names the required competitive
  inhibitors and allosteric regulators it could not place and the keywords that make a regulator
  optional (`optional_competitive_inhibitors`, `optional_allosteric_regulators`). Today the beam
  returns an empty result with no message. The first case this catches: a uni-uni reaction with
  S or P declared as a competitive inhibitor, whose only admissible site is free E and whose copy
  complex duplicates ES or EP. Whether C's merged seeds give that copy a placement is C's
  decision.

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
  (shared A group with abortive EAQ, A at E: rejected), a copy at {E, EQ} kept because EQ·A is
  new, and the ping-pong (segment, offsets) case.
- The split that would isolate EA → EAB·B\* from its copy group: the exact child set without it.
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
   The record reports per-level counts for every reaction.
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
