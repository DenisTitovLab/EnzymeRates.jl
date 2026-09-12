# Design: group-set RE→SS flips and context-bipartition splits

Date: 2026-09-11
Status: approved design, measured, not implemented
Branch: `refactor-re-to-ss-steps-to-only-keep-mechanisms-with-more-params`
Supersedes the RE→SS design on branch `spec-wegscheider-reachability-second-pass`
(`2026-07-23-re-to-ss-segment-splitting-flips-design.md`), whose generator this
document keeps and extends. The diagnosis in
`2026-07-23-split-canonicalization-reachability-bug.md` stands; its fix
directions are replaced by this document.

## Summary

Two expansion moves waste compute and, in one case, lose models.

The **split move** gives one step its own binding constant. In random-order
mechanisms a Wegscheider cycle usually forces the new constant back to the old
one, so `_canonical_mechanism` undoes the split and the move drops it. Because the
tied split never enters the frontier, the second split that would free a real
parameter is never taken. Random-order mechanisms can never reach independent
binding constants. Every split from a random-order bi-bi seed is dropped today.

The **RE→SS move** flips one kinetic group at a time. Once splits have separated
a metabolite's steps, a single-group flip can leave the rate equation unchanged.
Such children are fit anyway and enter the beam one parameter count above their
parent with a phantom parameter.

This design replaces both moves and deletes structural canonicalization:

1. The RE→SS move flips **whole kinetic groups**, in the **smallest sets of
   groups** whose joint flip raises the rapid-equilibrium segment count. A group
   holding no flux-carrying step never flips.
2. The split move divides one group by **binding context**: the steps whose
   source form already carries another ligand Y against the rest, one
   bipartition per Y. It emits the **smallest sets of simultaneous
   bipartitions** whose independent-parameter count rises, and the constraint
   solver's "count unchanged" verdict is a proof that the child's equation is
   the parent's.
3. `_canonical_mechanism` and `_merge_tied_kinetic_groups` are deleted. The
   compile-time equation dedup, which already exists and is tested, is the only
   dedup.

Enumeration stays structural and deterministic. No numerical identifiability
test runs inside the search; one lives in the test suite as a validator.

Measured on the 55 bi-bi seeds: the flip move emits 220 children, identical to
today; the split move emits 102; the richest seed reaches all 13 independent
groups and 9 independent parameters, which today's search never reaches; and
the closure of both moves reaches every one of the 1,434 segmentations of that
seed's 9 forms. On the largest ter-ter random-order seed (55 steps) the flip
move emits 6 children and the split search tries 8,878 candidate sets and
emits 12.

## Decisions

Denis made these choices during design, in this order.

- **Invariant.** No child of the flip or split move is emitted that is
  *provably* a reparameterization of its parent, and every rejection is a
  proof. Some no-op flips survive (see *What is accepted*); none is lost.
- **No numerical filter in enumeration.** A rank test could reject a real gain
  through a numerical wobble. It is confined to tests.
- **Flip whole groups, then minimal sets.** Denis first chose a per-metabolite
  rule, then an unrestricted segment cut. Measurement showed the unrestricted
  cut emits 53 children per parent on the richest bi-bi seed and tens of
  thousands on ter-ter random order, almost all of them cutting through a
  kinetic group. He then settled on whole-group flips with pair-and-larger
  sets when a single group cannot raise the segment count. Measurement showed
  the flip-plus-split closure reaches the same 1,434 segmentations as the
  unrestricted cut, at 4 to 7 children per parent instead of 53.
- **Split by binding context, not any bipartition.** Arbitrary bipartitions
  give a twelve-step group 2,047 ways to split and a million pair candidates
  on ter-ter. The context rule, "A's affinity may depend on which other ligand
  is already bound", gives at most one bipartition per other ligand and
  encodes the biochemical reason for the partition.
- **No level cap on the split search.** A cap could silently reintroduce the
  reachability loss for reactions with more substrates than were measured. The
  search terminates on its own (one bipartition per group at most), and the
  cost is paid by a faster predicate, not by a heuristic.
- **Predicate cost is a design requirement.** The cycle basis is computed once
  per parent, because a split moves no edges.

## The invariant, stated honestly

A child is emitted only if no proof shows its equation equals its parent's.
The proofs available:

| Move | Rejection | Why it is a proof |
|------|-----------|-------------------|
| split | independent-parameter count unchanged | the new constant is tied by a cycle; the derivation substitutes it away, so the rendered equation is the parent's |
| flip | the group set does not raise the RE segment count | every flipped step keeps both endpoints in one segment, lands on the diagonal of the segment matrix, and never reaches the equation |
| flip | a group holds no flux-carrying step | the group's steps lie in binding-only pendant blocks with zero net flux at steady state; only their equilibrium ratio can enter the equation, which the parent already has |

The reverse claim, that every emitted child gains, is **not** made. Uni-uni is
the clean counterexample: steady-state and rapid-equilibrium binding give the
same three-parameter rate law, so every flip of a uni-uni seed is a no-op. The
July measurement found 52 of 195 segment-raising single-group flips rank-flat
on split bi-bi parents. These children are emitted, fit once, and lose to their
parent on parsimony.

## Flux-carrying steps

A step is **flux-carrying** if it lies on a cycle of the step graph that contains
a chemistry step. Chemistry steps are the isomerization steps (`is_iso`).

Computation: take the undirected multigraph whose vertices are enzyme forms and
whose edges are all steps of the mechanism, RE and SS alike. Find its
biconnected blocks. A step is flux-carrying iff its block contains an
isomerization step. A bridge is its own block, so a dead-end leaf is never
flux-carrying. Two edges lie on a common cycle iff they share a block, so the
test is exact for the definition.

Why the definition is right: a binding-only cycle satisfies detailed balance
(the Wegscheider constraints the derivation enforces), so its circulation at
steady state is zero. A step whose every cycle is binding-only therefore carries
zero net flux, its forward and reverse rate enter the equation only through
their ratio, and flipping it to SS adds a phantom parameter and nothing else.
A step on a cycle with chemistry carries the reaction's net flux, and its two
rate constants enter separately.

A group that mixes a dead-end leaf with catalytic steps may flip whole: the leaf
shares the group's `kon` and `koff`, which the catalytic steps identify, so no
phantom appears. Only a group with **no** flux-carrying step is a proven no-op,
and after a split such a lone dead-end group can exist.

Flux-carrying-ness depends only on the step graph, not on RE/SS flags. No move
removes a catalytic step, and the moves that add steps add dead-end regulator
bindings and their mirrors, so it is recomputed per mechanism and never cached.

Competitive-inhibitor bindings stay RE by the existing rule. An inhibitor bound
to two forms that a catalytic step connects gets a mirror step, which puts the
inhibitor binding on a cycle with chemistry, so that rule is a modeling choice
("inhibitor binding is fast") rather than a consequence of this definition. The
docs say so.

## The flip move

`_expand_re_to_ss(m)` emits one child per minimal segment-raising set of
eligible groups.

**Eligible groups** are all-RE, bind no regulator, and hold at least one
flux-carrying step. Eligibility replaces today's mirror-class units
(`_re_to_ss_flip_units`): a catalytic step and its inhibitor mirror share a
group until a split separates them, and after that they are separate units
that a set can flip together when the segment count demands it.

**A child** flips every step of every group in the set to SS
(`_flip_group_to_ss`, unchanged). Reaction, allo-state tags, multiplicity, and
regulatory sites pass through unchanged.

**Minimal sets** are found Apriori-style, as July designed. Level 1 tries each
eligible group alone; a group that raises the RE segment count is emitted, one
that does not seeds level 2. Level `j` tries a `j`-set only if every
`(j−1)`-subset failed. The loop ends when a level fails nothing. The predicate
is `_compute_re_groups`, a union-find at microseconds per probe. On the 55
bi-bi seeds every single group raises the count (220 probes, 220 children); on
their split children the move emits 6.7 children per mechanism with pairs
included.

Why minimal sets suffice: a larger set with the same segment count adds only
segment-flat steps, so it is a reparameterization of the smaller set; a larger
set with a higher count is reached by expanding the smaller child.

For `AllostericMechanism` the segment count is taken on the A-state projection
(`_state_allo_mechanism(am, :A)`), which holds every group, and the
flux-carrying blocks on the shared step graph. See *Risks* for the I-state
projection.

## The split move

`_expand_split_kinetic_group(m)` emits minimal gaining sets of context
bipartitions.

**A context bipartition** of a group binding metabolite X, for another ligand
Y: the steps whose source form has Y bound, against the rest. Y ranges over
every ligand that appears on some source form of the group, substrates,
products, and regulators alike, except X itself. Both parts must be nonempty,
and two ligands that induce the same partition give one bipartition. For an
isomerization group Y ranges the same way, so a chemistry step and its
inhibitor-bound mirror split by the inhibitor. A group of twelve steps on
ter-ter random order has four context bipartitions instead of 2,047.

Where Y is a competitive inhibitor this is the catalytic-step-versus-mirror
split, so today's ability to separate a mirror from its counterpart is a case
of the same rule. For `AllostericMechanism` both parts inherit the group's
catalytic allo-state tag, as the carved part does today.

**The predicate** is the independent-parameter count from the thermodynamic
constraint solve on the concrete mechanism, without compiling. A child is
accepted iff its count exceeds the parent's.

**Predicate cost.** `_dependent_param_exprs_kernel` rebuilds the cycle basis
of the step graph (`_integer_nullspace` over `BigInt`) on every call, at 0.18 s
and heavy allocation on a 55-step mechanism. A split moves no edges, so the
cycle basis is the same for every candidate of one parent; only the
group-merge map changes. The move computes the basis once per parent and tests
each candidate by the rank of the merged constraint matrix. The plan factors
the kernel so the existing type-dispatching path and this path share the
basis computation, and measures the worst ter-ter seed against a time budget.
The measured search on that seed is 8,878 candidates; at the current solve
cost that is 26 minutes, and the target is seconds.

**Minimal sets** are found Apriori-style over `(group, bipartition)` pairs
from distinct groups, exhaustively.

1. Level 1: every context bipartition of every group with two or more steps.
   Accepted candidates are emitted; rejected ones seed level 2.
2. Level `j`: extend each rejected `(j−1)`-set by one bipartition of a group not
   yet in the set. Keep a `j`-set only if every `(j−1)`-subset was rejected at
   the previous level.
3. Stop when a level rejects nothing or produces no candidates. A set holds at
   most one bipartition per group, so the levels end at the number of
   splittable groups and the candidate count is bounded by the product over
   groups of (1 + bipartitions). No level cap: measurement shows random-order
   bi-bi and ter-ter gain at level 2 and never above, but that is two data
   points, and a cap set too low would silently reintroduce the reachability
   loss for a larger reaction.

**Sound pruning at level 2 and above.** A rejected bipartition's new constant is
tied through an RE cycle that passes through the carved steps. Breaking the tie
requires changing a group on that cycle, and every group on it lies in the same
RE segment as the carved steps. So the partners tried for a rejected set are the
RE groups in the RE segments of its carved steps. Groups holding an SS step are
never partners; measured on bi-bi, every split of such a group gains alone
(652/652) and every absorbed split carves an all-RE group (2064/2064).

**Fallback if the fast predicate is still too slow** on the worst ter-ter
seed: tie-guided extension. A rejected set is extended only by bipartitions of
groups named in the dependency expression that absorbed its new constant. This
is principled but unmeasured; it is designed, not built, and it replaces
nothing above.

Why minimal sets suffice: a larger set with the same count is a
reparameterization of the smaller one, and a larger set with a higher count is
reached by expanding the smaller child. Every partition reachable by context
bipartitions whose count exceeds the seed's is reached by a chain of minimal
gaining sets.

## What is deleted, what is untouched

Deleted, with their tests (listed under *Test surface*):

- `_canonical_mechanism`, both overloads
- `_merge_tied_kinetic_groups`, both overloads
- `_split_one_step` (replaced by the context-bipartition builder)
- `_re_to_ss_flip_units` and `_step_core` (replaced by group eligibility plus
  minimal sets)

Untouched:

- The compile-time dedup by `eq_hash` in `_compile_batch`, `_fit_batch`,
  `_offer_cv!`, and the LOOCV key.
- The structural `fitted` set that expands each structure once.
- `_select_beam`, `_select_count!`, and bucketing by `length(fitted_params)`.
- `init_mechanisms` and `seed_mechanisms`. Seeds are all-RE with one SS
  chemistry step, which is the coarsest model and the right start.
- `_flip_group_to_ss`.
- `rate_equation` derivation and its 0-allocation / sub-120 ns contract. This
  change is enumeration-time only.

## What is accepted

- **No-op flips survive.** Each costs one fit and one slot in a beam bucket
  above its true dimension, where it ties its parent on loss and loses on
  parsimony. Its descendants are structures the parent's descendants also
  reach, and the `fitted` set fits each structure once, so the waste is
  linear in the number of no-op children, not exponential.
- **Phantoms from no-op flips.** These are the July `fitted_params`
  over-count defect. The no-ops this design removes by proof (segment-flat,
  pendant) were 38 of the 90 measured phantoms. Whether the rest shrink is
  measured, not promised.
- **Kernel over-count on splits.** If a parent already carries a phantom, the
  kernel may accept a split that gains nothing. That is the safe direction.
- **Partitions outside the context rule are unreachable.** "A binds EB with
  its own constant and all other forms share" is reached only as an
  intersection of contexts after a second split, and some partitions are not
  reached at all. This is the documented modeling choice.
- **The reachable space is large.** Fixing the split bug exposes it: 1,434
  segmentations and about 13,000 structures in closure from the richest bi-bi
  seed, against 16 today. What the beam fits is set by children per parent,
  which stays at today's 4 flips plus about 2 splits at the seeds.
- **The renaming-duplicate problem** (~66 % of distinct-hash equations are
  renamings `eq_hash` misses) is untouched and unrelated.

## Measurements

Done before implementation, with throwaway scripts against the real
`Mechanism` API (session scratchpad, `option1.jl`, `option1b.jl`,
`segcount.jl`, `option2.jl`, `terter3.jl`). Structural only; no compiles.

| # | Question | Result |
|---|----------|--------|
| M2 | Flip children at bi-bi seeds, new vs today | 220 vs 220 |
| M2 | Segmentations reachable, richest bi-bi seed | 1,434 under flip + split, equal to every segmentation any edge assignment reaches |
| M3 | Split candidates, bi-bi seeds | 102 children, 426 solves, 4 s |
| M3 | Split candidates, 55-step ter-ter random order | 24 singles (0 gain), 240 pairs (12 gain), 8,614 at levels 3 to 6 (0 gain); 26 min at 0.18 s per solve |
| M3 | Split candidates, 47- and 39-step ter-ter | 0 singles gain, 11 and 10 pairs gain, 816 triples 0 gain |
| M6 | Richest bi-bi seed reaches 13 groups, 9 independent parameters | yes, in a closure of 16 structures |

Remaining, run during the plan. None is a suite test; the suite asserts exact
properties only.

| # | Question | Population | Pass criterion |
|---|----------|------------|----------------|
| M1 | Do the proof-based rejections ever reject a gain? | bi-bi seeds + depth-2 expansion, plus an allosteric bi-bi | every rejected split has the parent's `eq_hash` (exact); every ineligible-group flip applied by hand has the parent's rank (numerical) |
| M3′ | Fast predicate on the worst ter-ter seed | the 55-step seed, full search | seconds, not minutes; reported |
| M4 | Residual no-op flips | flip children at depth 1–2, bi-bi and an allosteric reaction | reported: fraction rank-flat |
| M5 | Phantoms | same children | reported: `length(fitted_params)` minus rank |

The numerical rank for M1, M4, and M5 is finite differences of `rate_equation`
over the fitted parameters at 60 concentration points, maximum over three
parameter draws, as in July. Numerical wobble there fails a measurement, never
loses a model.

## Test surface

**Delete** (they test deleted functions; listed for Denis's approval per the
"never delete a test" rule):

- `test/test_types.jl` "`_merge_tied_kinetic_groups` re-merges tied splits"
- `test/test_identify_rate_equation.jl` the `_canonical_mechanism` fixed-point
  and non-convergence tests
- `test/test_mechanism_enumeration.jl` "expand_mechanisms output is canonical"

**Update:**

- Exact split-move child counts in `test_mechanism_enumeration.jl`: children are
  now context-bipartition sets, and random-order seeds emit children.
- Exact RE→SS child counts at seeds stay at their current values (measured
  220 over the 55 bi-bi seeds); the mirror-lock fixture keeps its per-child
  property, gains a non-vacuity guard, and its premise changes: a mirror flips
  with its counterpart because they share a group or because the set demands
  it, not because of a mirror class.
- The six-named-moves pin compares `expand_mechanisms` against the same
  functions, so it survives as is.
- `test/fixtures/phase2_init_golden.txt` and `test_compile_budget.jl` are
  untouched; `init_mechanisms` does not change.

**Add:**

- Random-order bi-bi reaches the fully independent 9-fitted-parameter form
  through the split closure (the headline bug).
- A group of four or more steps is split into two parts of two or more by a
  context bipartition.
- Every emitted split child has a strictly higher kernel count than its
  parent, and no emitted set has an emitted proper subset (minimality).
- Every rejected split at level 1 compiles to the parent's `eq_hash` (the
  proof, tested exactly).
- The fast predicate agrees with the full kernel count on every candidate of
  the bi-bi seeds and on the split children.
- Every flip child has more RE segments than its parent; no emitted set has an
  emitted proper subset; a group with no flux-carrying step is never flipped.
  Checked at expansion depth two, allosteric included.
- The reachability test from July's design: a mechanism with two groups that
  each leave the segment count unchanged but split a segment together emits
  the two-group child.
- A flux-carrying oracle test on hand-built graphs: a dead-end leaf, a pendant
  inhibitor block, a mirror on the main cycle, a ping-pong second half.
- The uni-uni flips are emitted (documented no-op, not rejected).
- The 55-step ter-ter random-order seed's split expansion finishes under a
  time budget.

## Documentation

Three pages change. Apply the elements-of-style skill.

**`docs/src/identify/enumeration_engine.md`.** Rewrite moves 1 and 2 to
describe the new behavior: whole-group flips in minimal segment-raising sets,
flux-carrying steps; context bipartitions, minimal sets, the count test. Add a
**Modeling choices** section that states, in one paragraph each, with the
reason:

- Steady-state detail enters by flipping whole kinetic groups, the smallest set
  of groups that separates a new rapid-equilibrium segment. Which steps can
  flip together is set by how their constants are shared, so sharing structure
  is refined first and steady-state detail follows it.
- A group whose steps carry no net flux at steady state (a dead-end branch on
  its own) never flips, because the equation cannot tell the difference.
- Competitive-inhibitor binding stays at rapid equilibrium.
- A binding constant is split by context: the affinity for a metabolite may
  depend on which other ligand is already bound. Splits are never arbitrary
  partitions, a group is never split three ways in one move, and groups never
  merge.
- A child is never a provable copy of its parent; a few equation-identical
  children survive (uni-uni) and cost one fit each.

**`docs/src/identify/combinatorics.md`.** Rewrite for a biology audience.
Detach it from the engine: the page answers "why are there so many candidate
mechanisms for even a simple enzyme?" as a list of independent ways the count
multiplies, each with a concrete bi-bi number and a one-sentence biological
meaning. No summation formulas, no Stirling or Bell numbers, no ordered Bell
numbers; multiplication of small counts only. Sections: binding order; which
forms share an affinity (context splits); equilibrium versus steady-state
binding (the number of ways to divide the enzyme forms into equilibrated
regions, 1,434 for the random-order bi-bi with dead ends); dead-end complexes;
competitive inhibitors; allosteric states. Close with why the search is
filtered. Delete the "Known limitation" admonition. Correct the ter-ter seed
count (measured 35,665, the page says 42,220).

**`docs/src/developer.md`.** In *Enumeration engine architecture*, replace the
canonicalization story with the invariant, the proofs, the once-per-parent
cycle basis, and where the rank oracle lives (tests only).

## Non-goals

- The pivot-priority defect in `_assemble_constraints`; own PR, moves golden
  fixtures.
- Fixing `fitted_params` over-count beyond what the new flip move removes.
- Renaming duplicates that `eq_hash` misses.
- Beam-budget changes (Direction B′'s `fit_inherited` exclusion). With
  provable no-ops gone, duplicates in the frontier should be rare; revisit
  only if M4 says otherwise.
- A merge move, splits into more than two parts, or partitions outside the
  context rule.
- Any change to `rate_equation` or its performance contract.

## Risks

- **Fast predicate correctness.** Reusing the cycle basis and testing rank
  must agree with the full kernel count on every candidate. A dedicated test
  compares the two on the bi-bi population. If they disagree, the full kernel
  is the reference and the fast path is wrong.
- **Fast predicate still too slow.** Unknown until M3′. The tie-guided
  extension is the fallback and is designed, not built. A level cap is not a
  fallback.
- **Allosteric flux-carrying.** The I-state projection drops `:OnlyA` groups,
  so a group flux-carrying in the shared graph may sit in a pendant region of
  the I state. Flipping it then adds a phantom in the I-state polynomial. M4
  and M5 run on an allosteric reaction to size this; if it matters,
  eligibility is required in both projections.
- **Waste from surviving no-ops.** Bounded and measured (M4), not eliminated.
  If it is large on a real LDH or PFKP run, the numerical oracle can be
  reconsidered with exact rational arithmetic, which Denis declined for now.
- **Context rule coverage.** Two data points say pairs of context splits
  suffice on random order. A four-substrate reaction is unmeasured; the
  exhaustive search covers it at whatever cost the fast predicate allows.
