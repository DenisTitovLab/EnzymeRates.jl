# Design: RE→SS segment cuts and bipartition splits

Date: 2026-09-11
Status: approved design, not implemented
Branch: `refactor-re-to-ss-steps-to-only-keep-mechanisms-with-more-params`
Supersedes the RE→SS design on branch `spec-wegscheider-reachability-second-pass`
(`2026-07-23-re-to-ss-segment-splitting-flips-design.md`). The diagnosis in
`2026-07-23-split-canonicalization-reachability-bug.md` stands and is the
starting point here; its fix directions are replaced by this document.

## Summary

Two expansion moves waste compute and, in one case, lose models.

The **split move** gives one step its own binding constant. In random-order
mechanisms a Wegscheider cycle usually forces the new constant back to the old
one, so `_canonical_mechanism` undoes the split and the move drops it. Because the
tied split never enters the frontier, the second split that would free a real
parameter is never taken. Random-order mechanisms can never reach independent
binding constants. Every split from a random-order bi-bi seed is dropped today.

The **RE→SS move** flips one kinetic group at a time. Many flips leave the rate
equation unchanged. They are fit anyway, and they enter the beam one parameter
count above their parent with a phantom parameter.

This design replaces both moves and deletes structural canonicalization:

1. The RE→SS move is a **segment cut**: it divides one rapid-equilibrium
   segment into two connected halves and turns every rapid-equilibrium edge
   across the divide into a steady-state edge. Each child has exactly one more
   segment than its parent. A cut that only detaches a region carrying no net
   flux at steady state (a dead-end branch) is never made.
2. A split divides one group into **two parts, any way**, and the move emits the
   **smallest sets of simultaneous splits** whose independent-parameter count
   rises. The thermodynamic constraint solver decides, and its "count unchanged"
   verdict is a proof that the child's equation is the parent's.
3. `_canonical_mechanism` and `_merge_tied_kinetic_groups` are deleted. The
   compile-time equation dedup, which already exists and is tested, is the only
   dedup.

Enumeration stays structural and deterministic. No numerical identifiability
test runs inside the search; one lives in the test suite as a validator.

## Decisions

Denis made these choices during the design conversation.

- **Invariant.** No child of the flip or split move is emitted that is
  *provably* a reparameterization of its parent, and every rejection is a
  proof. Some no-op flips survive (see *What is accepted*); none is lost.
- **Segment cuts, any connected cut.** The RE→SS move is defined by its effect
  on the segment graph, which is all the equation sees: one more segment per
  child, never a cut that only detaches a dead-end region. Denis first chose a
  per-metabolite restriction (a metabolite binds at equilibrium everywhere or at
  steady state everywhere) and then replaced it with the unrestricted cut,
  because the per-metabolite rule cannot express models such as "free enzyme
  binds every substrate slowly, the ternary complex forms at equilibrium." The
  plan measures the size of the cut space; the per-metabolite rule is the
  documented fallback if the measurements are ugly.
- **Bipartition splits.** A split may divide a group of four steps 2 + 2, not
  only 1 + 3. The plan measures the candidate count on ter-ter before the
  search shape is final.
- **No numerical filter in enumeration.** A rank test could reject a real gain
  through a numerical wobble. It is confined to tests.

## The invariant, stated honestly

A child is emitted only if no proof shows its equation equals its parent's.
The proofs available:

| Move | Rejection | Why it is a proof |
|------|-----------|-------------------|
| split | independent-parameter count unchanged | the new constant is tied by a cycle; the derivation substitutes it away, so the rendered equation is the parent's |
| cut | one half of the cut holds no form on a chemistry cycle | that half is a union of binding-only pendant blocks; it carries zero net flux at steady state, sits at equilibrium with the other half, and only its equilibrium ratios can enter the equation |

Two further no-ops are excluded by construction rather than by test: an SS
step whose endpoints share a segment never enters the equation (it lands on the
diagonal of the segment matrix), and a cut adds only crossing edges, so every
new SS step joins two different segments.

The reverse claim, that every emitted child gains, is **not** made. Uni-uni is
the clean counterexample: steady-state and rapid-equilibrium binding give the
same three-parameter rate law, so every cut of a uni-uni seed is a no-op. The
July measurement found 52 of 195 segment-raising flips rank-flat on split bi-bi
parents. These children are emitted, fit once, and lose to their parent on
parsimony.

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

A **chemistry-cycle form** is an enzyme form incident to a flux-carrying step.
The cut move uses forms, not steps: a region of a segment with no
chemistry-cycle form is a union of pendant blocks and sits at equilibrium with
the rest of the segment however its edges are flagged.

Flux-carrying-ness depends only on the step graph, not on RE/SS flags. No move
removes a catalytic step, and the moves that add steps add dead-end regulator
bindings and their mirrors, so it is recomputed per mechanism and never cached.

Competitive-inhibitor bindings stay RE by the existing rule. An inhibitor bound
to two forms that a catalytic step connects gets a mirror step, which puts the
inhibitor binding on a cycle with chemistry, so that rule is a modeling choice
("inhibitor binding is fast") rather than a consequence of this definition. The
docs say so.

## The cut move

`_expand_re_to_ss(m)` emits one child per admissible cut of one RE segment.

**A cut** of a segment divides its forms into two nonempty halves, each
connected by the segment's RE edges. Enumerate the subsets that contain the
segment's root form (`_segment_root`, the form with the fewest bound
metabolites), keep those whose complement is nonempty and connected. Order the
cuts by the sorted names of the non-root half for deterministic output.
Segments with one form have no cut.

**Admissible** means both halves contain at least one chemistry-cycle form. A
cut failing this detaches a pendant equilibrium region and is a proven no-op.
The cut that isolates a dead-end leaf is the common case excluded here.

**The child** rebuilds every RE edge crossing the divide with
`is_equilibrium = false`. Crossing edges are grouped by their parent group: the
crossing members of one parent group become one new SS group with `kon`/`koff`,
and the non-crossing members stay RE in the parent group with its dissociation
constant. A group therefore splits as a side effect when the cut crosses only
some of its members. Grouping by parent group rather than by metabolite keeps
earlier splits intact and keeps the seed's one-constant-per-metabolite grouping
for a cut that crosses a whole group. Reaction, allo-state tags, multiplicity,
and regulatory sites pass through unchanged, as today. Splitting a group that
carries a catalytic allo-state tag gives both parts the tag.

Each child has exactly one more RE segment than its parent, by construction.
A test asserts it anyway.

Reachability: models with several cuts are reached by chains of single cuts.
The random-order square reaches "A binding at steady state everywhere" through
the cut {E, EB} | {EA, EAB}, and "free enzyme slow to bind anything" through
{E} | {EA, EB, EAB}, which flips one A edge and one B edge together. The
per-group flip of today reaches neither cleanly: it emits the segment-flat
single-edge flips as no-ops and needs a second no-op to close the cut.

For `AllostericMechanism` the segments and cuts are computed on the A-state
projection (`_state_allo_mechanism(am, :A)`), which holds every group, and the
flux-carrying blocks on the shared step graph. See *Risks* for the I-state
projection.

## The split move

`_expand_split_kinetic_group(m)` emits minimal gaining sets of bipartitions.

**A bipartition** divides one group of `n ≥ 2` steps into two nonempty parts.
There are `2^(n−1) − 1` of them; enumerate the subsets that contain the group's
first step to avoid counting each bipartition twice. For `n ≤ 3` these are
exactly today's single carves. For `AllostericMechanism`, both parts inherit
the group's catalytic allo-state tag, as the carved part does today.

**The predicate** is the independent-parameter count from the thermodynamic
constraint solve on the concrete mechanism, without compiling: the Wegscheider
rename map plus `_dependent_param_exprs_kernel` for `Mechanism`, and the combined
state solve for `AllostericMechanism`. The plan factors the existing
type-dispatching `_dependent_param_exprs` so the concrete path and the compiled
path share one function. A child is accepted iff its count exceeds the parent's.

**Minimal sets** are found Apriori-style. A candidate is a set of bipartitions
from distinct groups.

1. Level 1: every bipartition of every group with two or more steps. Accepted
   candidates are emitted; rejected ones seed level 2.
2. Level `j`: extend each rejected `(j−1)`-set by one bipartition of a group not
   yet in the set. Keep a `j`-set only if every `(j−1)`-subset was rejected at
   the previous level. Test, emit the accepted, keep the rejected for level
   `j+1`.
3. Stop when a level rejects nothing.

**Sound pruning at level 2 and above.** A rejected bipartition's new constant is
tied through an RE cycle that passes through the carved steps. Breaking the tie
requires changing a group on that cycle, and every group on it lies in the same
RE segment as the carved steps. So the partners tried for a rejected set are the
RE groups in the RE segments of its carved steps. Groups holding an SS step are
never partners; measured on bi-bi, every split of such a group gains alone
(652/652) and every absorbed split carves an all-RE group (2064/2064).

**Cost guard.** The ter-ter random-order seeds carry groups of four and more
steps, so level-2 and level-3 candidate counts can be large. The plan measures
them before the search shape is fixed. If they prove unaffordable, tie-guided
pruning is added: the kernel's dependency expression for the absorbed constant
names the groups on the tying cycle, and only those are tried. No silent cap;
if a budget is ever needed it logs.

Why minimal sets suffice: a larger set with the same count is a
reparameterization of the smaller one, and a larger set with a higher count is
reached by expanding the smaller child. Every partition reachable by
bipartitions whose count exceeds the seed's is reached by a chain of minimal
gaining sets.

## What is deleted, what is untouched

Deleted, with their tests (listed under *Test surface*):

- `_canonical_mechanism`, both overloads
- `_merge_tied_kinetic_groups`, both overloads
- `_split_one_step` (replaced by the bipartition builder)
- `_re_to_ss_flip_units`, `_step_core`, and `_flip_group_to_ss` (replaced by
  the cut enumeration; a cut crosses a step and its mirror together whenever
  their endpoints fall on the same sides, so no separate mirror lock is needed)

Untouched:

- The compile-time dedup by `eq_hash` in `_compile_batch`, `_fit_batch`,
  `_offer_cv!`, and the LOOCV key.
- The structural `fitted` set that expands each structure once.
- `_select_beam`, `_select_count!`, and bucketing by `length(fitted_params)`.
- `init_mechanisms` and `seed_mechanisms`. Seeds are all-RE with one SS
  chemistry step, which is the coarsest model and the right start.
- `rate_equation` derivation and its 0-allocation / sub-120 ns contract. This
  change is enumeration-time only.

## What is accepted

- **No-op cuts survive.** Each costs one fit and one slot in a beam bucket
  above its true dimension, where it ties its parent on loss and loses on
  parsimony. Its descendants are structures the parent's descendants also
  reach, and the `fitted` set fits each structure once, so the waste is
  linear in the number of no-op children, not exponential.
- **Phantoms from no-op cuts.** These are the July `fitted_params`
  over-count defect. The no-ops this design removes by proof (segment-flat,
  pendant) were 38 of the 90 measured phantoms. Whether the rest shrink is
  measured, not promised.
- **The cut space is larger than 2^metabolites.** Random-order bi-bi with
  products has on the order of 25 single cuts of the seed and a few hundred
  segmentations in all, against 16 under a per-metabolite rule. Denis accepted
  this on the grounds that it is no larger than what the per-group flip reaches
  today once no-ops are removed; M2 checks that.
- **Kernel over-count on splits.** If a parent already carries a phantom, the
  kernel may accept a split that gains nothing. That is the safe direction.
- **The renaming-duplicate problem** (~66 % of distinct-hash equations are
  renamings `eq_hash` misses) is untouched and unrelated.

## Measurements

These run during the plan, with scripts kept in the scratchpad, and their
numbers go into the docs. None is a suite test; the suite asserts exact
properties only.

| # | Question | Population | Pass criterion |
|---|----------|------------|----------------|
| M1 | Do the proof-based rejections ever reject a gain? | bi-bi seeds + depth-2 expansion, plus an allosteric bi-bi | every rejected split has the parent's `eq_hash` (exact); every inadmissible cut applied by hand has the parent's rank (numerical) |
| M2 | Size of the cut space against today's flips | 55 bi-bi seeds closed under the move, old flip vs new cut | reported: reachable structures and distinct equations, before and after; the new count should not exceed the old |
| M3 | Split candidate count | ter-ter random-order seeds, all levels | reported per level; decides whether tie-guided pruning is built |
| M4 | Residual no-op cuts | cut children at depth 1–2, bi-bi and an allosteric reaction | reported: fraction rank-flat |
| M5 | Phantoms | same children | reported: `length(fitted_params)` minus rank |
| M6 | Reachability | random-order bi-bi | the 13-group, 9-parameter fully independent form is reached by the split closure; a 2 + 2 partition of a four-step group is reached |

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
  now raw bipartition sets, and random-order seeds emit children.
- Exact RE→SS child counts and the mirror-lock fixture: one child per
  admissible cut, plus a non-vacuity guard on the loop. A cut crosses a
  catalytic step and its inhibitor mirror together whenever both endpoints'
  mirrors sit on the same sides, and the fixture is re-derived rather than
  assumed.
- The six-named-moves pin compares `expand_mechanisms` against the same
  functions, so it survives as is.
- `test/fixtures/phase2_init_golden.txt` and `test_compile_budget.jl` are
  untouched; `init_mechanisms` does not change.

**Add:**

- Random-order bi-bi reaches the fully independent 9-fitted-parameter form
  through the split closure (the headline bug).
- A four-step group is split 2 + 2.
- Every emitted split child has a strictly higher kernel count than its
  parent, and no emitted set has an emitted proper subset (minimality).
- Every rejected split at level 1 compiles to the parent's `eq_hash` (the
  proof, tested exactly).
- Every cut child has exactly one more RE segment than its parent; both halves
  of every emitted cut hold a chemistry-cycle form; a cut that isolates a
  dead-end leaf is never emitted; every new SS step joins two different
  segments. Checked at expansion depth two, allosteric included.
- The random-order square emits both {E, EB} | {EA, EAB} and
  {E} | {EA, EB, EAB}.
- A cut that crosses part of a group splits the group and leaves the rest RE.
- A flux-carrying oracle test on hand-built graphs: a dead-end leaf, a pendant
  inhibitor block, a mirror on the main cycle, a ping-pong second half.
- The uni-uni cuts are emitted (documented no-op, not rejected).
- A ter-ter random-order seed's split expansion finishes under a time budget.

## Documentation

Three pages change. Apply the elements-of-style skill.

**`docs/src/identify/enumeration_engine.md`.** Rewrite moves 1 and 2 to
describe the new behavior: the segment cut, chemistry-cycle forms, one more
segment per child; bipartitions, minimal sets, the count test. Add a
**Modeling choices** section that states, in one paragraph each, with the
reason:

- Steady-state detail enters one segment cut at a time. A rapid-equilibrium
  segment is divided into two connected halves and the binding steps across the
  divide become steady state; one enzyme form may bind a metabolite at
  equilibrium while another binds it at steady state.
- Regions that carry no net flux at steady state (dead-end branches) are never
  cut off, because the equation cannot tell the difference.
- Competitive-inhibitor binding stays at rapid equilibrium.
- A split divides one group into two parts; a group is never split three ways
  in one move, and groups never merge.
- A child is never a provable copy of its parent; a few equation-identical
  children survive (uni-uni) and cost one fit each.

**`docs/src/identify/combinatorics.md`.** Rewrite for a biology audience.
Detach it from the engine: the page answers "why are there so many candidate
mechanisms for even a simple enzyme?" as a list of independent ways the count
multiplies, each with a concrete bi-bi number and a one-sentence biological
meaning. No summation formulas, no Stirling or Bell numbers, no ordered Bell
numbers; multiplication of small counts only. Sections: binding order; which
forms share an affinity (the split factor, now partitions); equilibrium versus
steady-state binding (now `2^metabolites`, not `2^steps`); dead-end complexes;
competitive inhibitors; allosteric states. The equilibrium-versus-steady-state
section counts the ways to divide the enzyme forms into equilibrated regions,
with the bi-bi number from M2. Close with why the search is filtered. Delete
the "Known limitation" admonition. Use the M2 and M3 numbers.

**`docs/src/developer.md`.** In *Enumeration engine architecture*, replace the
canonicalization story with the invariant, the two proofs, and where the rank
oracle lives (tests only).

## Non-goals

- The pivot-priority defect in `_assemble_constraints`; own PR, moves golden
  fixtures.
- Fixing `fitted_params` over-count beyond what the new cut move removes.
- Renaming duplicates that `eq_hash` misses.
- Beam-budget changes (Direction B′'s `fit_inherited` exclusion). With
  provable no-ops gone, duplicates in the frontier should be rare; revisit
  only if M4 says otherwise.
- A merge move, or splits into more than two parts.
- Any change to `rate_equation` or its performance contract.

## Risks

- **Cut-space size.** A segment of `n` forms has up to `2^(n−1) − 1` cuts
  before the connectivity and admissibility filters, and ter-ter random order
  has segments of a dozen forms. Enumeration is a subset walk with two
  connectivity checks per subset, cheap at this size, but the number of emitted
  children per parent is unmeasured until M2. The per-metabolite rule is the
  documented fallback.
- **Ter-ter split cost.** Unknown until M3. The tie-guided pruning is the
  fallback and is designed, not built.
- **Allosteric flux-carrying.** The I-state projection drops `:OnlyA` groups,
  so a form on a chemistry cycle in the shared graph may sit in a pendant
  region of the I state. A cut admissible in the A state can then detach a
  pendant I-state region and add a phantom in the I-state polynomial. M4 and M5
  run on an allosteric reaction to size this; if it matters, admissibility is
  required in both projections.
- **Waste from surviving no-ops.** Bounded and measured (M4), not eliminated.
  If it is large on a real LDH or PFKP run, the numerical oracle can be
  reconsidered with exact rational arithmetic, which Denis declined for now.
