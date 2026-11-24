# Enumerate every feasible allosteric state: design

Date: 2026-11-24. Branch: `allosteric-state-enumeration`.

## Intent

The search should reach every allosteric version of a mechanism that the package accepts.
Today it reaches about a quarter of them. PDH and OGDH are ping-pong and allosteric, and
the versions they most plausibly need are among the missing: a K-type whose inactive
conformation runs the same chemistry with different affinities, and an inactive
conformation that stops one half-reaction.

The cause is two rules in the allosteric moves, not chemistry and not the derivation. This
design replaces both rules with one that proposes every state and lets the joint solve and
the guard decide.

## Terms

An allosteric mechanism gives each kinetic group one state. This document calls a group
**stopped** when it is `:OnlyA`, **shared** when it is `:EqualAI`, and says it has **its
own constants** when it is `:NonequalAI`. A **version** of a plain mechanism is one choice
of state for every group, not all shared. A **stop set** is a set of groups stopped with
every other group shared. An **own-constants set** is a set of groups given their own
constants.

## What the census established

A census ran on the 278 starting mechanisms of the bi-bi reaction
`A[CX] + B[N] → P[C] + Q[NX]`: 39 ping-pong and 239 sequential. These are the base tier,
what `candidate_mechanisms(rxn)` returns with no regulator declared: the seeds of
`init_mechanisms` with the 7 degenerate ping-pong seeds replaced by their 21 flip children.
The census gave every group every state, at one catalytic subunit with no regulatory site,
and asked the enumerator's exit test (`_enumerator_keeps`) about each version. Scripts and
outputs are on the maintainer's machine in `claude-scratchpad/allosteric-state-census/`.

**Today's moves reach a quarter of the kept versions.** On the 39 ping-pong parents the
exit test keeps 9,054 of 16,512 versions. The tagging, relaxing and regulator moves reach
2,169. Counted by chemistry state:

| Chemistry in the inactive conformation | Kept | Reached |
|:--|--:|--:|
| All stopped | 1,343 | 1,343 |
| All with its own constants | 738 | 370 |
| Part stopped, the rest shared | 2,241 | 456 |
| Part stopped, the rest with its own constants | 2,241 | 0 |
| All shared | 699 | 0 |
| Part shared, part with its own constants | 1,612 | 0 |
| All three mixed | 180 | 0 |

**The two rules behind it.** `_tagging_bases` starts every allosteric mechanism with
enough chemistry stopped to end the inactive conformation's turnover, and takes the
every-chemistry set alone whenever the guard admits it. `_expand_change_allo_state` then
gives all chemistry groups their own constants together, and only while no binding is
stopped. So chemistry is all stopped or all independent, never shared and never per
half-reaction.

**The missed versions are sound.** On the six-step ping-pong with both chemistry steps at
steady state, 25 versions were checked three ways: the joint solve for thermodynamic
consistency, the test suite's oracles for the law against the two-conformation mass-action
limit at interior points and on every face, and `_testhelper_identifiable_rank`. All 25
pass the first two.

**One binding cannot differ alone.** A conformation that runs the whole cycle must satisfy
the same `Keq`. When one binding has its own constant and every other group is shared, the
joint solve returns `K_I = K_A`: the conformations are identical and `L` cancels. The
fitted count shows it: the version fits the plain count plus one, as an all-shared version
would. Two bindings with their own constants fit one more, and thermodynamics ties the
second to the first (`K_I_P = K_A_P · K_I_A / K_A_A` when A and P differ).

**The exit test reads only the stops.** In all 32,706 versions that combine a stop with
own constants (13,390 ping-pong, 19,316 sequential), the exit test's answer equals its
answer on the stop set alone. No version without a stop is rejected (4,465 sequential, and
every ping-pong one). So every kept version is a kept stop set plus own constants. This
holds on the base tier. It fails on the degenerate ping-pong seeds, where the guard rejects
both half-reactions with their own constants; the search never fits or expands those seeds.

**Ten "all chemistry stopped" versions have no limit.** On 10 of the 39 ping-pong parents
the guard rejects the version that stops both half-reactions and nothing else (3 `_FILLED`,
7 `_BRIDGE`). The ε-path oracle finds no limit under any consistent reading, at zero B and
at zero A and B. These need the conformations to flip at the covalent form, which is
outside this design.

## Goals

- Every version the exit test keeps, and whose fitted count shows that its conformations
  differ, is reachable from its plain mechanism.
- Every child has more fitted constants than its parent.
- The allosteric moves hold fewer rules than today: no tagging bases, no relax-all
  chemistry rule, no clause about stopped bindings.

## Non-goals

- A conformational flip at a second free form (the ten rejected versions above).
- Making the guard admit more.
- A catalytic group tagged `:OnlyI`. The enumerator does not build one today and will not.
- Proofs that remove entries whose new constants the data cannot determine (see "Known
  costs").

## Anti-goals

- A move that changes a mechanism without adding a fitted constant.
- A bound on the entries that would hide a feasible version.
- A numerical identifiability test in `src/`.

## Constraints

- No public name changes. `candidate_mechanisms` returns more mechanisms, so
  `identify_rate_equation` can select differently on existing data.
- `Enumeration` calls only public names of `Derivation`. The gain test uses
  `Derivation.fitted_parameter_count`, as the split move does.
- The enumeration-engine test rules in `.claude/CLAUDE.md` apply: mechanisms written inline,
  exact children, every expected child written out.

## Approaches

**A, picked: one propose-and-test rule, with every stop set at entry.** The move that makes
a plain mechanism allosteric emits every stop set and every smallest own-constants set that
gains. The later move only gives further shared groups their own constants. The maintainer
picked it because every child then has more fitted constants than its parent, which is the
engine's principle. It costs about twice as many children when a mechanism becomes
allosteric.

**A with a stop move, set aside.** The entry would emit only the smallest sets, and a later
move would stop a shared group. A plain parent would get as many children as today. A stop
adds no fitted constant, so the move breaks the engine's principle, and a stop child
arriving at a parameter count whose beam floor is spent enters only by its loss.

**B, set aside: keep today's moves and add the missing transitions** (un-stop one chemistry
group, relax one chemistry group, relax a pair of bindings). It changes less code. A
shared-chemistry K-type would sit four or more moves from the plain mechanism, and nothing
would guarantee that every version is reachable.

**C, set aside: emit every kept version as a child of the plain mechanism.** It needs no
later move. It costs about 177 fits per ping-pong parent where A costs about 40.

## Design

### The entry move: `_expand_to_allosteric(c, rxn)`

For each oligomeric state the reaction allows, the move emits three kinds of child, each
with `c`'s catalytic scheme. The rule that a parent failing `_hyperbolic_catalysis` gets
children at one subunit only stays as it is, and so does the drop of a child that fails
`_hyperbolic_where_required`.

1. **Stop entries.** One child per nonempty set of kinetic groups: those groups `:OnlyA`,
   every other group `:EqualAI`, no regulatory site. The move emits all `2^g − 1` sets for
   `g` groups and tests none; the exits drop those `_enumerator_keeps` refuses, as they drop
   today's K-types.
2. **Own-constants entries.** One child per smallest set of groups, by inclusion, that
   gains when tagged `:NonequalAI` with every other group `:EqualAI`, no regulatory site. A
   set gains when its fitted count exceeds the plain mechanism's by more than one, the one
   being `L`. `_minimal_gaining_sets` finds the sets; a set that fails is extended by any
   other group.
3. **Regulator pairing.** Each child of kind 1 or 2 whose set holds chemistry groups only
   (`_chemistry_groups`) is emitted again with each declared allosteric regulator at a new
   site of the mechanism's multiplicity, once `:OnlyA` and once `:OnlyI`.

The maintainer set the pairing rule for stop entries: a version that changes only chemistry
shows nothing in a substrate or product titration that the plain mechanism lacks, so
without a regulator it can drop out of the beam before the regulator move reaches it. A
version that stops a binding changes the titration and survives on its own. Extending the
pairing to own-constants entries that change only chemistry is the author's call, by the
same argument; it awaits the maintainer's check.

### The later move: `_expand_change_allo_state(c)`

- **Catalytic groups.** One child per smallest set of `:EqualAI` groups, by inclusion, that
  raises the fitted count over `c`'s when tagged `:NonequalAI`. `_minimal_gaining_sets`
  finds the sets, as in the entry move; the two moves share that search.
- **Regulatory ligands.** One child per ligand relaxed to `:NonequalAI`, as today.
- **Filters.** The drops for a redundant copy group and for a non-hyperbolic law above one
  subunit stay.

A stopped group never changes state. A version with a group at its own constants where
another version stops it is reached from the entry that leaves the group shared.

### Why this reaches every version

Call a version **irreducible** when the exit test keeps it and making any one of its
own-constants groups shared lowers the fitted count. A version that is not irreducible
holds a group whose constants thermodynamics fixes to the shared ones, so it has the law of
the version without that group. The two moves should reach exactly the irreducible
versions.

Take an irreducible version with stop set `S` and own-constants set `N`. The exit test
reads only the stops, so `S` alone is kept, and it is an entry (or `S` is empty). From
there the groups of `N` are added in smallest gaining sets. What a group gains never falls
as other groups get their own constants: a group's own constants are new unknowns, the
thermodynamic rows that bind them are fixed, and the rank those rows add shrinks as more
columns are free. So a set that gained when added still gains in the finished version, and
while groups of `N` remain to add, some smallest set of them gains.

Both premises are measured, not proved. The exit test read only the stops in 32,706 of
32,706 versions on the bi-bi base tier. The argument about gains assumes that
`fitted_parameter_count` is the number of unknowns less the rank of the constraint rows. A
test keeps both (see "Tests").

### Seeds

Two changes, both needed so that a required regulator, the default, hides no version.

**Rule 1 of `_is_seed_node`** changes from "tags nothing `:NonequalAI`" to "no regulatory
ligand is `:NonequalAI`". `seed_mechanisms` runs the entry move and the regulator move,
never the later move, so a seed's catalytic `:NonequalAI` groups can only be an
own-constants entry's. Without the change every shared-chemistry K-type would be missing,
since no later move removes a stop.

**The seed search starts from the plain base tier**, not from `init_mechanisms`. Today it
runs the entry move on the 7 degenerate ping-pong seeds and flips the allosteric results
afterwards (`_base_tier`). A version of one of the 21 flip children is then a starting
mechanism only if the same states pass the exit test on the degenerate seed, which refuses
each one-half stop of the six-step ping-pong and both halves with their own constants. With
a regulator required no plain parent is ever fitted, so those versions are unreachable.
Measured with today's moves and one required regulator at one and two subunits: the
starting set grows from 6,428 to 7,126 and loses nothing; the six-step ping-pong whose Q
release is at steady state goes from 12 versions to 80. `candidate_mechanisms(reaction)`
still applies `_base_tier` to the seeds, for a dead-end regulator's sake.

This second change is the author's call. It answers a hole the design review found, and it
fixes the same hole in today's one-half bases.

### What goes

`_tagging_bases`, `_valid_tagging`, `_stops_inactive_turnover`, the K-type and V-type
branches of the entry move, the relax-all-chemistry rule and its "no `:OnlyA` binding"
clause, and the relaxation of a stopped group.

### A worked case

The six-step ping-pong with both chemistry steps at steady state has six groups: the
bindings of A and Q at free E, of B and P at the covalent form, and the two half-reactions.

- **Stop entries:** 63 emitted. The census parent of this shape keeps 42.
- **Own-constants entries:** eight. Each half-reaction alone (two), and six pairs of
  bindings: A with P, A with B, Q with B, Q with P, A with Q, B with P. A single binding
  gains nothing and is extended to its pairs.
- **Regulator pairing:** the three stop sets that hold only half-reactions (the first, the
  second, both) and the two single half-reactions with their own constants, each with each
  regulator in two states.

## Tests

Written first and watched to fail, under the enumeration-engine rules.

- **Entry move, exact children**, on fixtures small enough to write every child out: the
  uni-uni seed, the three-step ping-pong (three groups, seven stop sets) and one bi-bi
  sequential mechanism. One testset pins the regulator pairing, with a reaction that
  declares a regulator.
- **Own-constants entries on the six-step ping-pong**: the eight sets above, and that no
  single binding is among them.
- **Later move, exact children**: from a one-half stop, from a shared-chemistry K-type, and
  from a parent whose every group is stopped or independent (no child).
- **Reachability**, an aggregate check that says so in a comment: on the bi-bi base tier,
  or a sample if the full set is slow, the closure of the two moves equals the irreducible
  versions, found by brute force over every version; and the exit test's answer on a
  version equals its answer on the version's stop set. The test runs on the base tier, not
  on the degenerate seeds.
- **Seeds**: with a required regulator, the starting set holds, on the six-step ping-pong
  with both chemistry steps at steady state (one of the 21 flip children), a one-half stop
  and a shared-chemistry K-type, each with the regulator bound. The aggregate pin in
  "seed_mechanisms drops a seed the guard rejects" gets its new count.
- **Derivation**: one shared-chemistry K-type ping-pong, written out, compared with the
  mass-action steady state, in `test/Derivation/`.
- **Removed helpers**: the testsets in `test/test_oracles.jl` and
  `test/Enumeration/test_expansion_moves.jl` that call `_tagging_bases`, `_valid_tagging`
  or `_stops_inactive_turnover` are rewritten against the moves' children.

About 50 existing testsets pin the two moves' children (29 the entry move, 21 the later
move, in `test/Enumeration/test_expansion_moves.jl`), and `test/test_oracles.jl` builds
populations from the entry move's children. Their expected sets change. A six-group fixture
now has 63 stop entries per oligomeric state, too many to write out. How these testsets are
carried over is the maintainer's decision (see "Open for the maintainer"); whatever the
answer, a new expected set is reviewed against this design, not regenerated from the new
output.

## Guards and docs

- **Allocation bound.** `test/Enumeration/test_starting_mechanisms.jl` bounds the bytes one
  expansion of the bi-bi starting set allocates. The entry move builds more children, so
  the bound will trip. Re-pinning it needs the maintainer's approval.
- **Engine page** (`docs/src/identify/enumeration_engine.md`): moves 4 and 6 are rewritten
  around stop sets and own-constants sets, with the thermodynamic reason a binding cannot
  differ alone.
- **Developer page** (`docs/src/developer.md`): the sentence on where the moves run the
  guard, and the enumeration section where it names bases.
- **MWC page** (`docs/src/deriving/mwc_allostery.md`): the sentence "as the search's
  allosteric bases do".
- **`ROADMAP.md`**: the Known issue "An allosteric ping-pong whose covalent forms empty only
  by binding B cannot stop both half-reactions" shrinks to what the census found: ten
  starting mechanisms, no limit under any reading, fixed only by a flip at the covalent
  form. Its three sub-entries go. The Planned entry on the hyperbolic drops loses its
  mention of bases. `docs/src/roadmap.md` keeps its entry.

## Search cost

Children per plain parent at one oligomeric state, without a regulator, after the exits
(medians):

| Population | Today | Stop entries | Own-constants entries |
|:--|--:|--:|--:|
| Bi-bi sequential (239) | 7 | 11 | 3 |
| Bi-bi ping-pong (39) | 13 | 32 | 8 |
| Ter-bi (60 sampled of 5,158) | 12 | 23 | not measured |

- **Growth.** Stop entries about double with each added group: 13, 23, 53 and 90 kept at
  four, five, six and seven groups on ter-bi. Today's K-types double with each binding
  group.
- **Guard time.** All stop sets of the 60 ter-bi parents took about a second.
- **Regulator pairing.** A sequential parent has one chemistry-only stop set (median), a
  ping-pong parent three. Each adds two children per declared regulator.

## Known costs

- **Entries the data cannot determine at one subunit.** The gain test is the fitted count.
  On the six-step ping-pong, the pairs A with Q and B with P fit two constants more than the
  plain mechanism and have its rank at one subunit, one more at two. The L-phantom
  certificate does not cover them, since their inactive conformation turns over, so they are
  fitted. `ROADMAP.md` gets an entry with the count over the bi-bi starting mechanisms.
- **Optional regulators on entries that stop a binding.** These are emitted bare. At one
  subunit an L phantom among them has its parent's loss and enters the beam only through
  the width floor, as today's K-types do.
- **Selections can change** on data already analysed, since more mechanisms compete.

## Interaction with the iso-enumeration design

`docs/superpowers/specs/2026-11-24-iso-enumeration-design.md` (branch `iso-enumeration`)
states its allosteric behaviour in terms of tagging bases. Under this design the same
decisions read: the flip is not a chemistry group, so a stop set or own-constants set that
holds it gets no regulator pairing, and the later move gives it its own constants like any
group. Whichever branch lands second rewrites that section and rebases
`expansion_moves.jl`.

## Left to the builder

- How the two moves share the own-constants search, and whether its extension needs a
  narrower rule than "any other shared group" once measured on ter-ter parents.
- The order of the children, as long as it does not depend on the number of workers.
- The fixtures for the exact-children tests, within the test rules.
- Whether the reachability test runs on all 278 bi-bi starting mechanisms or a sample.
- A recovery check during implementation: data generated from a shared-chemistry K-type
  ping-pong and from an "all chemistry stopped" mechanism with a regulator, each recovered
  by `identify_rate_equation`. It need not become a suite test.

## Open for the maintainer

1. How the existing testsets are carried over. The author's proposal: a testset that pins a
   removed rule (a tagging base, the every-chemistry set, the relax-all rule) is deleted
   and named in the pull request, its place taken by the new tests of the rule that
   replaces it. A testset about another aspect (the hyperbolic drop, a multiplicity, one
   regulator) keeps its assertion and moves to a fixture of at most four groups, or to one
   oligomeric state, so that every child is still written out.
2. The regulator pairing covers own-constants entries that change only chemistry (author's
   call above). Confirm or narrow to stop entries.
3. Approval to re-pin the allocation bound of one bi-bi expansion.
4. The version: results change and no public name does.
