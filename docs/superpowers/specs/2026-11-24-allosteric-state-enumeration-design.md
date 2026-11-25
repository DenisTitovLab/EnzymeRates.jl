# Enumerate every feasible allosteric state: design

Date: 2026-11-25. Branch: `allosteric-state-enumeration`.

## Intent

The search should reach every allosteric version of a mechanism that the package accepts.
Today it reaches about a quarter of them. PDH and OGDH are ping-pong and allosteric, and
the versions they most plausibly need are among the missing: a K-type whose inactive
conformation runs the same chemistry with different affinities, and an inactive
conformation that stops one half-reaction.

The cause is two rules in the allosteric moves, not chemistry and not the derivation. This
design replaces both. A mechanism becomes allosteric by stopping any set of groups in the
inactive conformation, and a later move gives groups their own constants there. The guard
decides which stops are valid and the fitted count decides which own constants are real.

The maintainer added a second aim: a competitive inhibitor may bind the inactive
conformation alone. Today it binds both conformations or the active one alone.

## Terms

An allosteric mechanism gives each kinetic group one state. This document calls a group
**stopped** when it is `:OnlyA`, **shared** when it is `:EqualAI`, and says it has **its
own constants** when it is `:NonequalAI`. An **inhibitor group** is a kinetic group whose
steps bind a competitive inhibitor, and it alone may be **inactive-only**, `:OnlyI`. A
**version** of a plain mechanism is one choice of state for every group, not all shared. A
**stop set** is a set of groups stopped with every other group shared. An **own-constants
set** is a set of groups given their own constants. A set **gains** when giving its groups
their own constants raises the fitted count, and a **smallest gaining set** has no subset
that gains. An **entry** is a child of the move that makes a plain mechanism allosteric.

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
2,169. Counted by law, they reach 941 of 4,609. By chemistry state, in versions:

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

**A stop is not always a valid mechanism.** On the bi-bi mechanism that binds A then B,
the guard refuses an inactive conformation that binds neither and keeps the chemistry
(`_BRIDGE`): stopping a substrate requires stopping the chemistry too. The same two
bindings with their own constants and the chemistry shared are kept. Of the smallest
own-constants sets that add a fitted constant, the exit test refuses the stop of 414 of
953 on the sequential parents and of 98 of 222 on the ping-pong ones.

**Ten "all chemistry stopped" versions have no limit.** On 10 of the 39 ping-pong parents
the guard rejects the version that stops both half-reactions and nothing else (3 `_FILLED`,
7 `_BRIDGE`). The ε-path oracle finds no limit under any consistent reading, at zero B and
at zero A and B. In all ten the first half-reaction is at steady state and the second
half's chemistry is a rapid-equilibrium step. Each keeps the stop of its first half. These
need the conformations to flip at the covalent form, which is outside this design.

**Today's later move makes children that add no constant, and none costs a fit.** Over two
rounds of `_expand_change_allo_state` from today's entries, 3,687 of 25,509 kept children
on the base tier add no fitted constant, and 2,907 of 35,395 on starting mechanisms with
two regulators, separate or merged into one site. Each is one binding relaxed alone, and
each has its parent's law key, so it shares its parent's fit. Every relaxation of a
regulatory ligand adds one constant (17,331 of 17,331): a regulatory site's constants enter
no thermodynamic cycle.

**An inhibitor group can be inactive-only.** Two censuses ran on starting mechanisms of the
same reaction with a competitive inhibitor `I`. The first gave every group every state on
42 of the 1,548 mechanisms that bind `I`. The second took 61 of the 17,060 that bind `I`
and a copy of `A`, and gave them every entry of the design below.

- The exit test reads only the stops and the inactive-only groups: its answer on a version
  equals its answer with the own-constants groups shared, in 29,832 of 29,832 versions.
- An inactive-only inhibitor adds no constant when the stops cut every form it binds off
  the inactive conformation. It then binds nowhere: 24 such entries, each at the plain
  mechanism's fitted count. Every other inactive-only entry adds one constant.
- Giving an inactive-only inhibitor its own constants is always kept and always adds one
  constant (765 of 765).
- A reactant's copy made inactive-only can be redundant (`_redundant_copy_groups`): 14 of
  3,892 kept entries with one inactive-only group. No entry without one is (0 of 3,804).

## Goals

- Every version the exit test keeps, each of whose own-constants groups adds to the fitted
  count, is reachable from its plain mechanism.
- Every child of the two moves has more fitted constants than its parent.
- The allosteric moves hold fewer rules than today: no tagging bases, no relax-all
  chemistry rule, no clause about stopped bindings.

## Non-goals

- A conformational flip at a second free form (the ten rejected versions above).
- Making the guard admit more.
- `:OnlyI` on a group that binds a substrate or product at its catalytic site, abortive
  complexes included. Only an inhibitor group takes it.
- Proofs that remove versions whose new constants the data cannot determine (see "Known
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

**A, picked: every stop set at entry, own constants later.** The move that makes a plain
mechanism allosteric emits every stop set. The later move gives groups their own
constants, from shared or from stopped, in the smallest sets that gain. The maintainer
picked it because every child then has more fitted constants than its parent, which is the
engine's principle. It costs about twice as many children when a mechanism becomes
allosteric.

**A with a stop move, set aside.** The entry would emit only the smallest stop sets, and a
later move would stop a shared group. A plain parent would get as many children as today.
A stop adds no fitted constant, so the move breaks the engine's principle, and a stop
child arriving at a parameter count whose beam floor is spent enters only by its loss.

**B, set aside: keep today's moves and add the missing transitions** (un-stop one chemistry
group, relax one chemistry group, relax a pair of bindings). It changes less code. A
shared-chemistry K-type would sit four or more moves from the plain mechanism, and nothing
would guarantee that every version is reachable.

**C, set aside: emit every kept version as a child of the plain mechanism.** It needs no
later move. It gives a ping-pong parent about 150 children, with 117 distinct laws, where
A gives about 35.

### Own constants at entry

A stop entry adds one fitted constant, `L`. An own-constants entry adds two. The
maintainer asked whether the entry could be stops alone. Three ways, measured on the bi-bi
base tier at one subunit:

| Entry | Misses | Entries adding two constants, per parent (median) |
|:--|:--|--:|
| Stops alone | 140 of 4,574 ping-pong laws; 566 of 14,999 sequential versions | 0 |
| Stops, and own constants where the stop is refused | nothing | 2 sequential, 3 ping-pong |
| Stops, and every smallest gaining set | nothing | 3 sequential, 8 ping-pong |

The maintainer picked the middle row. Stops alone misses versions with no stop, mostly the
shared-chemistry K-type of a sequential mechanism, on 227 of the 239 sequential parents and
28 of the 39 ping-pong ones. The bottom row adds a route that does not pass through a stop
entry (see "Known costs") for one more entry per sequential parent and five per ping-pong
parent.

### Regulator pairing

The maintainer set it: an entry that stops only chemistry is also emitted with each
declared regulator. Pairing own-constants versions too was set aside. It would put an
inactive conformation with slower chemistry and a regulator in the starting set, at one or
two more starting mechanisms per parent, regulator state and oligomeric state. That
version is one later move from the paired stop entry.

## Design

### The entry move: `_expand_to_allosteric(c, rxn)`

For each oligomeric state the reaction allows, the move emits three kinds of child, each
with `c`'s catalytic scheme. The rule that a parent failing `_hyperbolic_catalysis` gets
children at one subunit only stays as it is, and so does the drop of a child that fails
`_hyperbolic_where_required`.

1. **Stop entries.** One child per nonempty set of kinetic groups: those groups `:OnlyA`,
   every other group `:EqualAI`, no regulatory site. The move emits all `2^g − 1` sets for
   `g` groups and tests none; the exits drop those `_enumerator_keeps` refuses, as they drop
   today's K-types. Each inhibitor group has a third choice, `:OnlyI`, so with `k` inhibitor
   groups the move emits `2^(g−k) · 3^k − 1` children; one with no stop and an inactive-only
   inhibitor is an entry too. The move drops a child with an inactive-only group in two
   cases: no step of the group runs in the inactive conformation
   (`Mechanisms.steps_running_in_allosteric_conformations`), so the inhibitor binds
   nowhere; or the child holds a redundant copy group (`_redundant_copy_groups`).
2. **Own-constants entries, where the stop is refused.** Take each smallest set of
   groups, by inclusion, that gains when tagged `:NonequalAI` with every other group
   `:EqualAI`. A set gains when its fitted count exceeds the plain mechanism's by more than
   one, the one being `L`. `_minimal_gaining_sets` finds the sets; a set that fails is
   extended by any other group. The move emits the set with its own constants, no
   regulatory site, only when `_enumerator_keeps` refuses the same set stopped. Where the
   stop is kept, it is a kind-1 entry and the later move reaches the own-constants version
   from it.
3. **Regulator pairing.** Each stop entry whose set holds chemistry groups only
   (`_chemistry_groups`), with no inactive-only group, is emitted a second time with each
   declared allosteric regulator at a new site of the mechanism's multiplicity, once
   `:OnlyA` and once `:OnlyI`. The bare entry stays.

Kind 2 exists because a stop is not always a valid mechanism (see the census). The K-type
whose inactive conformation binds A and B with their own affinities and runs the shared
chemistry has no stop to grow from, so it enters directly.

The reason for kind 3: a version that stops only chemistry shows nothing in a substrate or
product titration that the plain mechanism lacks, so without a regulator it can drop out
of the beam before the regulator move reaches it. A version that stops a binding changes
the titration and survives on its own. No paired child holds a `:NonequalAI` group.

### The later move: `_expand_change_allo_state(c)`

- **Catalytic groups.** One child per smallest set of groups not yet `:NonequalAI`, by
  inclusion, that raises the fitted count over `c`'s when tagged `:NonequalAI`. The groups
  are `:EqualAI`, `:OnlyA` or `:OnlyI`, bindings and chemistry alike.
  `_minimal_gaining_sets` finds the sets, as in the entry move; the two moves share that
  search. A set that holds a stopped or inactive-only group changes what the exit test
  reads, so its child is counted only if `_enumerator_keeps` keeps it.
- **Regulatory ligands.** One child per ligand relaxed to `:NonequalAI`, as today.
- **Filters.** The drops for a redundant copy group and for a non-hyperbolic law above one
  subunit stay.

A group never returns to `:EqualAI`, `:OnlyA` or `:OnlyI`. A stop is the limit of own
constants as the inactive affinity or rate goes to zero, and an inactive-only inhibitor is
the limit as its active affinity does. Where the inactive conformation binds a ligand
weakly, the stop entry fits better than the shared one at the same fitted count, the beam
may hold only the stop entry, and the step from `:OnlyA` lets that entry grow into the true
version.

### Why this reaches every version

Call a version **irreducible** when the exit test keeps it, each of its inactive-only
groups binds in the inactive conformation, and making any one of its own-constants groups
shared lowers the fitted count. A version that is not irreducible
holds a group whose own constants add nothing: thermodynamics fixes them to the shared
ones, or the inactive conformation never reaches the group's forms and they do not show.
It has the law of the version with that group shared. The two moves should reach every
irreducible version.

Take an irreducible version with stop set `S` and own-constants set `N`. Read an
inactive-only group as one more member of `S`: the exit test reads it as it reads a stop.

- **`S` is not empty.** The exit test reads only the stops, so `S` alone is kept and is an
  entry. From there the groups of `N` are added in smallest gaining sets. What a group
  gains never falls as other groups get their own constants: a group's own constants are
  new unknowns, the thermodynamic rows that bind them are fixed, and the rank those rows
  add shrinks as more columns are free. So a set that gained when added still gains in the
  finished version, and while groups of `N` remain to add, some smallest set of them gains.
- **`S` is empty.** `N` holds a smallest gaining set `T`. If the exits refuse `T` stopped,
  `T` with its own constants is a kind-2 entry. If they keep it, it is a stop entry, and
  the later move takes its groups to their own constants. The rest of `N` follows as
  above.

The premises are measured, not proved. The exit test read only the stops in 32,706 of
32,706 versions on the bi-bi base tier. The argument about gains assumes that
`fitted_parameter_count` is the number of unknowns less the rank of the constraint rows.
That a kept stop of a smallest gaining set always grows into that set's own-constants
version is the third premise; it held on every parent measured. A test keeps all three
(see "Tests").

The conclusion is measured too, from the census's fitted counts and law keys at one
subunit. On the 39 ping-pong parents the two moves reach all 5,861 irreducible versions,
which hold 4,574 distinct laws, and 52 further versions that each share a law, and so a
fit, with one of them. On the 239 sequential parents they reach all 14,999. At two
subunits, on the 14 ping-pong parents that allow them, they reach all 2,214. On the 42
parents with an inhibitor they reach all 14,898, of which 3,719 hold the inhibitor
inactive-only.

### Other moves on a parent with an inactive-only group

The flip, split, dead-end and regulator moves run on such a parent as on any other. Two
places assume that the whole catalytic graph is the active conformation's, which an
inactive-only group ends:

- **The dead-end move** adds its inhibitor group shared, as today. On an allosteric parent
  it runs `_redundant_copy_groups` only when the copy has a twin at each binding, with the
  twins read on the whole graph (`_dead_end_child`). On a parent with an inactive-only
  group it runs the test without that shortcut.
- **`_assert_emission_rules`** reads a steady-state group's flux on the whole graph. That
  stays right for a step the active conformation runs: a cycle through an inhibitor's
  complexes follows mirrored steps, and replacing each by the step it mirrors gives a
  cycle of the active conformation. This is an argument, not a measurement. Its docstring
  loses the sentence that no move tags a catalytic group `:OnlyI`.

### Seeds

Two changes, both needed so that a required regulator, the default, hides no version.

**Rule 1 of `_is_seed_node`** changes from "tags nothing `:NonequalAI`" to "no regulatory
ligand is `:NonequalAI`". `seed_mechanisms` runs the entry move and the regulator move,
never the later move, so a seed's catalytic `:NonequalAI` groups can only be a kind-2
entry's. Without the change a required regulator would hide every kind-2
entry.

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

The second change is the author's call. It answers a hole the design review found, and
today's one-half bases have the same hole.

### What goes

`_tagging_bases`, `_valid_tagging`, `_stops_inactive_turnover`, the K-type and V-type
branches of the entry move, and the later move's two rules for catalytic groups: one
binding at a time with no gain test, and all chemistry together while no binding is
`:OnlyA`.

### A worked case

The six-step ping-pong with both chemistry steps at steady state has six groups: the
bindings of A and Q at free E, of B and P at the covalent form, and the two half-reactions.

- **Stop entries:** 63 emitted. The census parent of this shape keeps 42.
- **Smallest gaining sets:** eight. Each half-reaction alone (two), and six pairs of
  bindings: A with P, A with B, Q with B, Q with P, A with Q, B with P. A single binding
  gains nothing and is extended to its pairs.
- **Own-constants entries:** three of the eight, the pairs whose stop the exits refuse: A
  with B, Q with P, B with P. The other five are stop entries, and the later move takes
  each to its own constants. A stopped half-reaction goes alone. The stop of A and P goes
  as a pair, since A alone beside a stopped P is refused.
- **Regulator pairing:** the three stop sets that hold only half-reactions (the first, the
  second, both), each with each regulator in two states.

## Tests

Written first and watched to fail, under the enumeration-engine rules.

- **Entry move, exact children**, on fixtures small enough to write every child out: the
  uni-uni seed, the three-step ping-pong (three groups, seven stop sets) and one bi-bi
  sequential mechanism. One testset pins the regulator pairing, with a reaction that
  declares a regulator.
- **Own-constants entries on the six-step ping-pong**: the three sets above, and that no
  single binding and no set with a kept stop is among them.
- **Later move, exact children**: from a one-half stop (the half to its own constants among
  them), from the stop of A and P (both to their own constants in one child), from a
  shared-chemistry K-type, and from a parent whose every group has its own constants (no
  child).
- **Reachability**, an aggregate check that says so in a comment: on the bi-bi base tier
  and on a sample of the starting mechanisms with a competitive inhibitor, or on samples
  of both if the full sets are slow, the closure of the two moves holds every irreducible
  version, found by brute force over every version; every child has more fitted constants
  than its parent; and the exit test's answer on a version equals its answer on the
  version's stop set. The test runs on the base tier, not on the degenerate seeds.
- **Inactive-only inhibitor**, exact children: the entry move on the uni-uni seed with a
  competitive inhibitor bound at free E (four groups); a fixture whose stops cut the
  inhibitor's only site off the inactive conformation, with no inactive-only child for
  those stop sets; a fixture with a reactant's copy whose inactive-only child is dropped
  as redundant; the later move taking an inactive-only inhibitor to its own constants; and
  the dead-end move adding a second inhibitor to a parent that holds an inactive-only one.
  One mechanism with an inactive-only inhibitor, written out, is compared with the
  mass-action steady state in `test/Derivation/`, unless the suite already holds one.
- **Seeds**: with a required regulator, the starting set holds, on the six-step ping-pong
  with both chemistry steps at steady state (one of the 21 flip children), a one-half stop
  and the kind-2 entry of A with B, each with the regulator bound. The aggregate pin in
  "seed_mechanisms drops a seed the guard rejects" gets its new count.
- **Derivation**: one shared-chemistry K-type ping-pong, written out, compared with the
  mass-action steady state, in `test/Derivation/`.
- **Removed helpers**: the testsets in `test/test_oracles.jl` and
  `test/Enumeration/test_expansion_moves.jl` that call `_tagging_bases`, `_valid_tagging`
  or `_stops_inactive_turnover` are rewritten against the moves' children.

About 50 existing testsets pin the two moves' children (29 the entry move, 21 the later
move, in `test/Enumeration/test_expansion_moves.jl`), and `test/test_oracles.jl` builds
populations from the entry move's children. Their expected sets change. A six-group fixture
now has 63 stop entries per oligomeric state, too many to write out. The maintainer
approved how they are carried over:

- A testset that pins a removed rule (a tagging base, the every-chemistry set, the
  relax-all rule) is deleted and named in the pull request. The new tests of the rule that
  replaces it take its place.
- A testset about another aspect (the hyperbolic drop, a multiplicity, one regulator)
  keeps its assertion and moves to a fixture of at most four groups, or to one oligomeric
  state, so that every child is still written out.

A new expected set is reviewed against this design, not regenerated from the new output.

## Guards and docs

- **Allocation bound.** `test/Enumeration/test_starting_mechanisms.jl` bounds the bytes one
  expansion of the bi-bi starting set allocates. The entry move builds more children, so
  the bound will trip. The maintainer approved re-pinning it; the pull request states the
  bytes before and after.
- **Version.** 0.11.1: results change and no exported function does.
- **Engine page** (`docs/src/identify/enumeration_engine.md`): moves 4 and 6 are rewritten
  around stop sets and own-constants sets, with the thermodynamic reason a binding cannot
  differ alone and the three states an inhibitor group can enter with.
- **Developer page** (`docs/src/developer.md`): the sentence on where the moves run the
  guard, and the enumeration section where it names bases.
- **MWC page** (`docs/src/deriving/mwc_allostery.md`): the sentence "as the search's
  allosteric bases do".
- **`ROADMAP.md`**: the Known issue "An allosteric ping-pong whose covalent forms empty only
  by binding B cannot stop both half-reactions" shrinks to what the census found: ten
  starting mechanisms, no limit under any reading, fixed only by a flip at the covalent
  form. Its three sub-entries go. `docs/src/roadmap.md` keeps its entry.
- **`ROADMAP.md`, other entries on these moves.**
  - Planned, "Drop `change_allo_state` relaxations that free no parameter": removed. The
    later move emits only sets that gain, and a regulatory ligand's relaxation always adds
    one constant (see the census).
  - Planned, "Decide whether to keep the hyperbolic drops above one subunit": decided, the
    drops stay (see "Open for the maintainer"). The entry moves to "Decided against" with
    the reason and the populations measured here.
  - Ideas, "Apply the chain-flank exclusion at one catalytic subunit": stays, with what
    this work measured. On versions of the 15 base-tier parents that hold a steady-state
    chain, 126 of 128 flank flips at one subunit add no identifiable constant, and all 34
    whose parent is an L phantom. Two flips of a parent with a stop and an own-constants
    group add one, so the exclusion cannot cover every one-subunit mechanism.
  - Known issue, "The guard rejects some allosteric mechanisms whose law is the limit":
    its count of ter-bi starting mechanisms with no allosteric child is measured again.
  - Planned, new: let the dead-end move give its inhibitor any of the three states on an
    allosteric parent (see "Open for the maintainer").
  - Known issues, new: one for each of the first, second and fourth items under "Known
    costs".

## Search cost

Children per plain parent at one oligomeric state, without a regulator, after the exits
(medians):

| Population | Today | Stop entries | Own-constants entries |
|:--|--:|--:|--:|
| Bi-bi sequential (239) | 7 | 11 | 2 |
| Bi-bi ping-pong (39) | 13 | 32 | 3 |
| Ter-bi (60 sampled of 5,158) | 12 | 23 | not measured |

- **Growth.** Stop entries about double with each added group: 13, 23, 53 and 90 kept at
  four, five, six and seven groups on ter-bi. Today's K-types double with each binding
  group.
- **Guard time.** All stop sets of the 60 ter-bi parents took about a second.
- **Regulator pairing.** A sequential parent has one chemistry-only stop set (median), a
  ping-pong parent three. Each adds two children per declared regulator.
- **Later move.** Over the 5,861 ping-pong versions it emits about 16,000 children, a
  median of three per mechanism.
- **Inhibitors.** A parent with one inhibitor group keeps 25 stop entries and 12 with the
  inhibitor inactive-only (medians of 42 parents). With two inhibitor groups it keeps 47
  and 60 (61 parents). A starting set that requires an inhibitor and an allosteric
  regulator grows by the same factor.

## Known costs

- **Versions the data cannot determine at one subunit.** The gain test is the fitted count.
  On the six-step ping-pong, the pairs A with Q and B with P with their own constants fit
  two constants more than the plain mechanism and have its rank at one subunit, one more at
  two. The L-phantom certificate does not cover them, since their inactive conformation
  turns over, so they are fitted. The `ROADMAP.md` entry gives the count over the bi-bi
  starting mechanisms.
- **Inactive-only entries the data cannot determine at one subunit.** Of 122 sampled
  entries with the inhibitor inactive-only, 47 fit one constant more than their rank at
  one subunit, and the L-phantom certificate covers 4 of them: `L` and the inhibitor's
  constant enter the law only as a ratio, and at `L = 0` the inhibitor leaves the law. So
  43 are fitted for nothing. The same 122 with the inhibitor active-only: 28, of which 19
  are certified. At two subunits 3 of 122 are undetermined.
- **Optional regulators on entries that stop a binding.** These are emitted bare. At one
  subunit an L phantom among them has its parent's loss and enters the beam only through
  the width floor, as today's K-types do.
- **Own-constants versions behind a phantom stop.** Where a smallest gaining set's stop is
  kept, its own-constants version is reached only through that stop entry. At one subunit
  the stop entry is an L phantom for 227 of the 539 such sets on the sequential parents
  and 11 of the 124 on the ping-pong ones: it has its parent's loss and enters the beam
  only through the width floor. Above one subunit there are no phantoms. Emitting every
  smallest gaining set at entry would remove this cost (see "Own constants at entry").
- **A steady-state group that runs only beside an inactive-only inhibitor.** A split can
  put a mirrored step in a group of its own, and a flip can then take it to steady state.
  With the inhibitor inactive-only that group runs in the inactive conformation alone, and
  where the stops end that conformation's turnover its two rate constants enter as a
  ratio. The flux test reads the whole graph and does not see it.
- **Selections can change** on data already analysed, since more mechanisms compete.

## Interaction with the iso-enumeration design

`docs/superpowers/specs/2026-11-24-iso-enumeration-design.md` (branch `iso-enumeration`)
states its allosteric behaviour in terms of tagging bases. Under this design the same
decisions read: the flip is not a chemistry group, so a stop set that holds it gets no
regulator pairing, and the later move gives it its own constants like any group.
Whichever branch lands second rewrites that section and rebases `expansion_moves.jl`.

## Left to the builder

- How the two moves share the search for smallest gaining sets, and whether its extension
  needs a narrower rule than "any other group" once measured on ter-ter parents.
- Where the later move asks `_enumerator_keeps` about a set that holds a stopped group,
  before it reads the fitted count.
- The order of the children, as long as it does not depend on the number of workers.
- The fixtures for the exact-children tests, within the test rules.
- Whether the reachability test runs on all 278 bi-bi starting mechanisms or a sample.
- A recovery check during implementation: data generated from a shared-chemistry K-type
  ping-pong and from an "all chemistry stopped" mechanism with a regulator, each recovered
  by `identify_rate_equation`. It need not become a suite test.

## Open for the maintainer

The maintainer decided the approach, the entry rule, the regulator pairing, how the
existing testsets are carried over, the allocation bound, the version, and that an
inhibitor group can be inactive-only. Three calls in this document are the author's and
wait for a yes or no:

- **The hyperbolic drops stay.** Above one subunit both moves drop a child whose
  conformations' laws are not hyperbolic. The drops remove nothing the exits keep: on
  hyperbolic parents every one of 13,678 kept base-tier versions passes the test, and
  every one of 7,846 kept versions and entries on parents with inhibitors. But the exits
  refuse the failing candidates as bridges, and the Known issue on the guard proposes
  admitting more bridges. Removing the drops would rest a modeling rule on that
  conservatism, to save two filter calls and two testsets.
- **The dead-end move keeps adding its inhibitor shared.** Every version is reachable
  without more: add the inhibitor to the plain mechanism, then enter. An inhibitor added
  after entry reaches active-only or inactive-only binding only as a limit of its own
  constants, at one constant more. Giving the move all three states would triple its
  children on allosteric parents and put a second place in the engine that chooses
  exclusive states.
- **The seed search starts from the plain base tier** (see "Seeds").
