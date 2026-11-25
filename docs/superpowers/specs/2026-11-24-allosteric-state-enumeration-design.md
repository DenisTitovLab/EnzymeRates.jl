# Enumerate every feasible allosteric state: design

Date: 2026-11-25. Branch: `allosteric-state-enumeration`.

## The problem

An allosteric enzyme has two conformations, active and inactive. A mechanism says, for
each group of steps that share rate constants, what the inactive conformation does with
them: it does not run them, it runs them with the active conformation's constants, or it
runs them with constants of its own. The package accepts many such combinations for one
catalytic scheme. The search that `identify_rate_equation` runs builds about a quarter of
them.

The missing ones matter. PDH and OGDH are ping-pong and allosteric. An inactive
conformation that runs the same chemistry and binds the substrates more weakly is never
built for a ping-pong mechanism. Neither is one that stops only one of the two
half-reactions.

## The change

A mechanism becomes allosteric by choosing any set of groups that the inactive
conformation does not run. A second step then gives groups constants of their own in the
inactive conformation, a few groups at a time, and only where that adds a constant a fit
can estimate. Every combination the package accepts is then built.

The maintainer added four things along the way:

- A competitive inhibitor may bind the inactive conformation alone, whether it is added
  before or after the mechanism becomes allosteric.
- The rule that skips a pointless change of a binding step to steady state applies to
  allosteric mechanisms as it does to plain ones.
- The rule that an enzyme with several subunits has a hyperbolic catalytic scheme is
  applied where the enumerator returns its mechanisms, not inside three moves.
- When a reaction requires an allosteric regulator, the starting mechanisms are built
  from the same plain mechanisms the search otherwise starts from.

## Names used in this document

### Allosteric states

A **kinetic group** is a set of steps that share their constants. An allosteric mechanism
stores one state per kinetic group (`m.allosteric.allosteric_states`):

| State | The inactive conformation | The active conformation |
|:--|:--|:--|
| `:EqualAI` | runs the group with the same constants | runs it |
| `:OnlyA` | does not run the group | runs it |
| `:NonequalAI` | runs the group with constants of its own | runs it |
| `:OnlyI` | runs the group | does not run it |

A plain mechanism has no states. `L` is the equilibrium constant between the two
conformations, one more fitted constant.

### Functions

All but the ones marked public are internal to `Enumeration`. The moves are in
`src/Enumeration/expansion_moves.jl`; the docs page "The enumeration engine"
(`docs/src/identify/enumeration_engine.md`) numbers them.

| Function | What it does |
|:--|:--|
| `candidate_mechanisms(reaction)` | Public. Returns the starting mechanisms of a reaction. |
| `candidate_mechanisms(reaction, parent)` | Public. Returns the children of `parent`: it runs every move below and filters the results. |
| `_expand_to_allosteric(c, rxn)` | Move 4. Makes a plain mechanism allosteric: each child has `c`'s steps and a state for every kinetic group. |
| `_expand_change_allo_state(c)` | Move 6. Changes the state of kinetic groups, or of ligands at a regulatory site, of an allosteric mechanism. |
| `_expand_add_dead_end_regulator(c, rxn)` | Move 3. Adds a competitive inhibitor: new dead-end complexes, in one new kinetic group. |
| `_expand_add_allosteric_regulator(c, rxn)` | Move 5. Adds an allosteric regulator at a regulatory site. |
| `_expand_re_to_ss(c)` | Move 1. Changes kinetic groups from rapid equilibrium to steady state. |
| `_expand_split_kinetic_group(c)` | Move 2. Splits a kinetic group so that its parts get separate constants. |
| `_enumerator_keeps(m)` | The test every returned mechanism passes: `Derivation.check_valid_rate_equation(Bool, m)`, which the docs call the well-posedness guard, and `_kept_forms_keep_a_binding(m)`. In `candidate_mechanism_type.jl`. |
| `Derivation.fitted_parameter_count(m)` | Public. The number of constants a fit of `m` estimates. This document calls it the **fitted count**. |
| `Derivation.rate_equation_dedup_key(m)` | Public. A key that two mechanisms share when they have one rate equation. The search fits one mechanism per key. |
| `_minimal_gaining_sets(gains, partners)` | Finds the smallest sets, by inclusion, that pass a test. Moves 1 and 2 use it today. |
| `_chemistry_groups(c)` | Says which kinetic groups hold chemistry: an isomerization, a fused binding or a Theorell–Chance step. |

Other internal functions are introduced where they first appear.

### Three terms this document adds

- A **version** of a plain mechanism is one choice of state for every kinetic group, not
  all `:EqualAI`.
- A set of kinetic groups **gains** when tagging its groups `:NonequalAI` raises the
  fitted count. A **minimal gaining set** has no subset that gains.
- A version is **irreducible** when no group's `:NonequalAI` state is idle; "Why every
  version is reached" defines it.

## What was measured

The measurements ran on the 278 starting mechanisms of the bi-bi reaction
`A[CX] + B[N] → P[C] + Q[NX]`: 39 ping-pong and 239 sequential. This is what
`candidate_mechanisms(rxn)` returns when no regulator is declared. It comes from
`init_mechanisms`, which builds the simplest mechanisms, and `_base_tier`, which replaces
each degenerate one by its children under move 1 (both in
`src/Enumeration/starting_mechanisms.jl`). A mechanism is degenerate (`_degenerate`) when
its rate equation ignores a substrate or never saturates; 7 ping-pong mechanisms are, and
they have 21 such children. This document calls the 278 the **base tier**, after the
function.

Each measurement gave every kinetic group every state, at one catalytic subunit with no
regulatory site, and asked `_enumerator_keeps` about each version. Scripts and outputs are
on the maintainer's machine in `claude-scratchpad/allosteric-state-census/`.

**Today's moves reach a quarter of the kept versions.** On the 39 ping-pong mechanisms
`_enumerator_keeps` keeps 9,054 of 16,512 versions. Moves 4, 5 and 6 reach 2,169. Counted
by rate equation, they reach 941 of 4,609. By the state of the chemistry groups:

| Chemistry groups | Kept | Reached |
|:--|--:|--:|
| All `:OnlyA` | 1,343 | 1,343 |
| All `:NonequalAI` | 738 | 370 |
| Some `:OnlyA`, the rest `:EqualAI` | 2,241 | 456 |
| Some `:OnlyA`, the rest `:NonequalAI` | 2,241 | 0 |
| All `:EqualAI` | 699 | 0 |
| Some `:EqualAI`, the rest `:NonequalAI` | 1,612 | 0 |
| All three states | 180 | 0 |

**The two rules behind it.** `_expand_to_allosteric` calls `_tagging_bases`, which starts
every allosteric mechanism with enough chemistry groups `:OnlyA` to end the inactive
conformation's turnover. It takes the set of all chemistry groups alone whenever
`_enumerator_keeps` accepts it, and smaller sets only when it does not.
`_expand_change_allo_state` then tags all chemistry groups `:NonequalAI` together, and
only while no binding group is `:OnlyA`. So a chemistry group is `:EqualAI` only where the
set of all chemistry groups was refused, and no chemistry group becomes `:NonequalAI`
without the others.

**The missed versions are sound.** On the six-step ping-pong with both chemistry steps at
steady state, 25 versions were checked three ways: the solve for thermodynamic
consistency, the test suite's oracles, which compare the rate equation with the
mass-action model of the two conformations, and `_testhelper_identifiable_rank`
(`test/testhelpers.jl`), which counts the constants that data determine. All 25 pass the
first two.

**One binding cannot be `:NonequalAI` alone.** A conformation that runs the whole cycle
must give the same `Keq`. When one binding group is `:NonequalAI` and every other group is
`:EqualAI`, thermodynamics returns `K_I = K_A`: the conformations are identical and `L`
cancels. The fitted count shows it: the version fits the plain count plus one, as an
all-`:EqualAI` version would. Two binding groups tagged `:NonequalAI` fit one more, and
thermodynamics ties the second to the first (`K_I_P = K_A_P · K_I_A / K_A_A` when A and P
differ).

**`_enumerator_keeps` reads only the `:OnlyA` groups.** In all 32,706 versions that
combine `:OnlyA` with `:NonequalAI` (13,390 ping-pong, 19,316 sequential), its answer
equals its answer with every `:NonequalAI` group made `:EqualAI`. No version without an
`:OnlyA` group is rejected (4,465 sequential, and every ping-pong one). So every kept
version is a kept set of `:OnlyA` groups plus `:NonequalAI` groups. This holds on the base
tier. It fails on the 7 degenerate ping-pong mechanisms. `identify_rate_equation` never
fits those; today only `seed_mechanisms` expands them, and "Starting mechanisms" ends
that.

**A set of `:OnlyA` groups is not always a valid mechanism.** On the bi-bi mechanism that
binds A then B, `check_valid_rate_equation` refuses an inactive conformation that binds
neither and keeps the chemistry (its `_BRIDGE` rejection): tagging a substrate's binding
`:OnlyA` requires tagging the chemistry `:OnlyA` too. The same two binding groups tagged
`:NonequalAI`, with the chemistry `:EqualAI`, are kept. Of the minimal gaining sets,
`_enumerator_keeps` refuses the set tagged `:OnlyA` for 414 of 953 on the sequential
mechanisms and for 98 of 222 on the ping-pong ones.

**Ten versions with all chemistry `:OnlyA` have no limit.** On 10 of the 39 ping-pong
mechanisms `check_valid_rate_equation` rejects the version that tags both half-reactions
`:OnlyA` and nothing else (3 as `_FILLED`, 7 as `_BRIDGE`). The test suite's oracle finds
no mass-action limit under any reading, at zero B and at zero A and B. In all ten the
first half-reaction is at steady state and the second half's chemistry is a
rapid-equilibrium step. Each keeps the version with its first half `:OnlyA`. These need
the conformations to interconvert at the covalent form, which is outside this design.

**`_expand_change_allo_state` makes children that add no constant today, and none costs a
fit.** Over two rounds of the move from today's children of `_expand_to_allosteric`, 3,687
of 25,509 kept children on the base tier add no fitted constant, and 2,907 of 35,395 on
starting mechanisms with two regulators, at separate sites or at one. Each is one binding
group tagged `:NonequalAI` alone, and each has its parent's dedup key, so it shares its
parent's fit. Tagging a regulatory ligand `:NonequalAI` always added one constant (17,331
of 17,331): a regulatory site's constants enter no thermodynamic cycle.

**A competitive inhibitor's group can be `:OnlyI`.** `_expand_add_dead_end_regulator`
puts an inhibitor's binding steps in one kinetic group; this document calls it the
inhibitor's group. Two measurements ran on starting mechanisms of the same reaction with a
competitive inhibitor `I`. The first gave every group every state on 42 of the 1,548
mechanisms that bind `I`. The second took 61 of the 17,060 that bind `I` and a second copy
of `A` as an inhibitor, and gave them every kind (a) child of the design below.

- `_enumerator_keeps` reads only the `:OnlyA` and `:OnlyI` groups: its answer on a version
  equals its answer with every `:NonequalAI` group made `:EqualAI`, in 29,832 of 29,832
  versions.
- An `:OnlyI` inhibitor adds no constant when the `:OnlyA` groups leave the inactive
  conformation none of the enzyme forms the inhibitor binds. It then binds nowhere: 24
  such versions, each at the plain mechanism's fitted count. Every other version with an
  `:OnlyI` inhibitor and no `:NonequalAI` group adds one constant.
- Tagging an `:OnlyI` inhibitor `:NonequalAI` is always kept and always adds one constant
  (765 of 765).
- A reactant bound a second time as an inhibitor can become redundant once its group is
  `:OnlyI`: 14 of 3,892 kept versions with one `:OnlyI` group, and 0 of 3,804 without
  one. `_redundant_copy_groups` (`candidate_mechanism_type.jl`) finds such groups. An
  inhibitor is redundant when each of its complexes has a constant multiple of the weight
  of a form the mechanism already has, so that data cannot separate the inhibitor's
  constant. No mechanism the enumerator returns may hold one.

**An inhibitor added to an allosteric mechanism can take three states.**
`_expand_add_dead_end_regulator` tags the new group `:EqualAI` today. It ran on 241 of
today's children of `_expand_to_allosteric`, at one and two subunits, from 47 of the 278
base-tier mechanisms, and gave 1,316 children. Each was rebuilt with the inhibitor's group
`:OnlyA` and `:OnlyI`.

- All 1,316 `:OnlyA` children are kept, and each adds one constant, as the `:EqualAI`
  child does.
- 1,261 `:OnlyI` children are kept, bind in the inactive conformation and add one
  constant. The other 55 bind nowhere and add none.
- In those 55 the `:EqualAI` and the `:OnlyA` child have one dedup key: the inactive
  conformation has none of the forms the inhibitor binds. In the other 1,261 the three
  children have three keys.
- No child has a redundant inhibitor.

**Changing a chain's flank to steady state seldom shows in an allosteric mechanism.** A
chain is a steady-state isomerization, alone in its kinetic group, between two enzyme
forms that each have one other step, a binding. Those two bindings are its flanks.
`_chain_flank_groups` (`expansion_moves.jl`) returns them for a plain mechanism, and
`_expand_re_to_ss` skips them: the docs prove that a flank at steady state adds constants
and no new rate equation ("Apparent constants on a steady-state chain",
`docs/src/deriving/re_vs_ss.md`). For an allosteric mechanism the function returns
nothing, so flanks are changed. Of the 278 base-tier mechanisms 15 hold a chain. On four
versions of each, the fitted count was compared with the number of constants data
determine. At one subunit 126 of 128 children that change only flanks add no such
constant; the two that add one belong to a version with an `:OnlyA` and a `:NonequalAI`
group. At two subunits 76 of 96 add none and 20 add one.

## Goals

- Every version that `_enumerator_keeps` keeps, that has no redundant inhibitor, and each
  of whose `:NonequalAI` groups adds to the fitted count, is reachable from its plain
  mechanism.
- Every child of `_expand_to_allosteric`, `_expand_change_allo_state` and
  `_expand_add_dead_end_regulator` has more fitted constants than its parent.
- A competitive inhibitor's group can be `:OnlyI`.
- A required allosteric regulator hides no version that an optional one would reach.
- The allosteric moves hold fewer rules than today: no `_tagging_bases`, no rule that all
  chemistry groups change together, no clause about `:OnlyA` binding groups.
- `_expand_re_to_ss` skips chain flanks in plain and allosteric mechanisms alike.

## Non-goals

- Conformations that interconvert at a second free enzyme form (the ten rejected versions
  above).
- Making `check_valid_rate_equation` accept more. The maintainer left that for a change of
  its own.
- Removing allosteric mechanisms at one subunit. The maintainer keeps them (see "Known
  costs").
- `:OnlyI` on a group that binds a substrate or product at its catalytic site, abortive
  complexes included. Only a competitive inhibitor's group takes it.
- Proofs that remove versions whose new constants data cannot determine (see "Known
  costs").

## Anti-goals

- A child of the three moves named in "Goals" with no more fitted constants than its
  parent.
- A bound on the children of `_expand_to_allosteric` that would hide a feasible version.
- A numerical identifiability test in `src/`.

## Constraints

- No public name changes. `candidate_mechanisms` returns more mechanisms, so
  `identify_rate_equation` can select differently on existing data.
- `Enumeration` calls only public names of `Derivation`. The test for a gain uses
  `Derivation.fitted_parameter_count`, as `_expand_split_kinetic_group` does.
- The enumeration-engine test rules in `.claude/CLAUDE.md` apply: mechanisms written inline,
  exact children, every expected child written out.

## Approaches

**A, picked: every set of `:OnlyA` groups from `_expand_to_allosteric`, `:NonequalAI` from
`_expand_change_allo_state`.** `_expand_to_allosteric` emits a child for every set of
groups tagged `:OnlyA`. `_expand_change_allo_state` tags groups `:NonequalAI`, from
`:EqualAI` or from `:OnlyA`, in the minimal gaining sets. The maintainer picked it because
every child then has more fitted constants than its parent, as the children of every
other move but the merge of regulatory sites have. A plain mechanism gets about twice as
many allosteric children.

**A with a move that tags a group `:OnlyA` later, set aside.** `_expand_to_allosteric`
would emit only the smallest `:OnlyA` sets, and a later move would tag an `:EqualAI` group
`:OnlyA`. A plain mechanism would get as many children as today. But tagging a group
`:OnlyA` adds no fitted constant.

**B, set aside: keep today's moves and add the missing changes of state** (one chemistry
group from `:OnlyA` to `:EqualAI`, one chemistry group to `:NonequalAI`, a pair of binding
groups to `:NonequalAI`). It changes less code. A mechanism whose inactive conformation
runs the same chemistry with other affinities would sit four or more moves from the plain
mechanism, and nothing would guarantee that every version is reachable.

**C, set aside: emit every kept version as a child of the plain mechanism.** It needs no
second move. It gives a ping-pong mechanism about 150 children, with 117 distinct rate
equations, where A gives about 35.

### `:NonequalAI` children of `_expand_to_allosteric`

A child with only `:OnlyA` and `:EqualAI` groups adds one fitted constant, `L`. A child
with `:NonequalAI` groups adds two. The maintainer asked whether `_expand_to_allosteric`
could emit `:OnlyA` sets alone. Three ways, measured on the base tier at one subunit:

| `_expand_to_allosteric` emits | Misses | Children adding two constants, per parent (median) |
|:--|:--|--:|
| `:OnlyA` sets alone | 140 of 4,574 ping-pong rate equations; 566 of 14,999 sequential versions | 0 |
| `:OnlyA` sets, and a minimal gaining set as `:NonequalAI` where it is refused as `:OnlyA` | nothing | 2 sequential, 3 ping-pong |
| `:OnlyA` sets, and every minimal gaining set as `:NonequalAI` | nothing | 3 sequential, 8 ping-pong |

The maintainer picked the middle row. The top row misses versions with no `:OnlyA` group,
mostly the sequential mechanism whose inactive conformation runs the same chemistry with
other affinities, on 227 of the 239 sequential mechanisms and 28 of the 39 ping-pong ones.
The bottom row adds a second route to some versions (see "Known costs") for one more child
per sequential mechanism and five per ping-pong one.

### Children with a regulator

The maintainer set it: a child of `_expand_to_allosteric` that tags only chemistry groups
`:OnlyA` is also emitted with each declared allosteric regulator. A version with chemistry
groups `:NonequalAI` and a regulator is not emitted this way. It is one
`_expand_change_allo_state` step from a child that is.

## Design

### `_expand_to_allosteric(c, rxn)`

For each oligomeric state the reaction allows, the move emits three kinds of child, each
with `c`'s catalytic scheme.

**(a) Every set of `:OnlyA` groups.** One child per nonempty set of kinetic groups: those
groups `:OnlyA`, every other group `:EqualAI`, no regulatory site. The move emits all
`2^g − 1` sets for `g` groups and does not call `_enumerator_keeps` on them;
`candidate_mechanisms(reaction, parent)` drops those it refuses, as it does with every
move's children today.

A competitive inhibitor's group has a third choice, `:OnlyI`. With `k` such groups the
move emits `2^(g−k) · 3^k − 1` children. One with no `:OnlyA` group and an `:OnlyI`
inhibitor is among them. The move drops a child with an `:OnlyI` group in two cases:

- No step of the group runs in the inactive conformation, so the inhibitor binds nowhere.
  `Mechanisms.steps_running_in_allosteric_conformations(m)`, a public function, says
  which steps each conformation runs.
- `_redundant_copy_groups` finds a redundant inhibitor in the child.

Nothing stands in for a dropped child. No mechanism holds an inhibitor that cannot bind,
whether the inhibitor is required or optional. The same `:OnlyA` set is still emitted with
the inhibitor's group `:EqualAI` and `:OnlyA`. The same `:OnlyA` set without the inhibitor
is a child of the plain mechanism without it, which the search holds when the inhibitor is
optional, and `_expand_add_dead_end_regulator` adds the inhibitor to that child as to any
parent.

**(b) A minimal gaining set tagged `:NonequalAI`, where the same set tagged `:OnlyA` is
not returned.** Take each minimal gaining set, with every other group `:EqualAI` and no
regulatory site. Here a set gains when its fitted count exceeds the plain mechanism's by
more than one, the one being `L`, and its child would be returned (see "When a child would
be returned"). `_minimal_gaining_sets` finds the sets: a set that does not gain is
extended by any other group. The move emits the set tagged `:NonequalAI` only when the
same set tagged `:OnlyA`, at the same oligomeric state, would not be returned. Where it
would, that is a kind (a) child, and `_expand_change_allo_state` reaches the `:NonequalAI`
version from it.

Kind (b) exists because a set of `:OnlyA` groups is not always a valid mechanism (see
"What was measured"). The inactive conformation that binds A and B with affinities of its
own and runs the same chemistry has no `:OnlyA` child to grow from, so it is emitted
directly.

**(c) A kind (a) child with a regulator.** Each kind (a) child whose `:OnlyA` groups are
all chemistry groups (`_chemistry_groups`), and that has no `:OnlyI` group, is emitted a
second time with each declared allosteric regulator at a new regulatory site, once with
the regulator `:OnlyA` and once `:OnlyI`. The child without the regulator stays.

The reason for kind (c): a version that tags only chemistry `:OnlyA` shows nothing in a
substrate or product titration that the plain mechanism lacks. Without a regulator it can
drop out of the search before `_expand_add_allosteric_regulator` reaches it. A version
with an `:OnlyA` binding group changes the titration and survives on its own. No kind (c)
child holds a `:NonequalAI` group.

A parent that fails `_hyperbolic_catalysis` gets children at one subunit only, as today
(see "The hyperbolic rule").

### `_expand_change_allo_state(c)`

- **Kinetic groups.** One child per minimal gaining set among the groups not yet
  `:NonequalAI`: its groups are tagged `:NonequalAI`. The groups may be `:EqualAI`,
  `:OnlyA` or `:OnlyI`, binding and chemistry groups alike. Here a set gains when its
  child has a higher fitted count than `c` and would be returned.
  `_minimal_gaining_sets` finds the sets, as in `_expand_to_allosteric`; the two moves
  share that search.
- **Regulatory ligands.** One child per ligand tagged `:NonequalAI`, kept only when its
  fitted count exceeds `c`'s. The maintainer asked for the check on every change of
  state. It never fired in the measurements (17,331 of 17,331 add one constant), so it
  covers cases they did not sample.

No group leaves `:NonequalAI`. `:OnlyA` is the limit of `:NonequalAI` as the inactive
conformation's affinity or rate goes to zero, and `:OnlyI` is the limit as the active
conformation's affinity does. Where the inactive conformation binds a ligand weakly, the
`:OnlyA` child fits better than the `:EqualAI` one at the same fitted count, the search
may keep only the `:OnlyA` child, and the step from `:OnlyA` to `:NonequalAI` lets it grow
into the true version.

### When a child would be returned

Both searches for minimal gaining sets ask whether a set's child would be returned by
`candidate_mechanisms(reaction, parent)`. It would when all three hold:

- `_enumerator_keeps` keeps it;
- it passes `_hyperbolic_where_required` (see "The hyperbolic rule");
- `_redundant_copy_groups` finds no redundant inhibitor in it.

A set whose child fails one of them does not gain, and the search extends it by any other
group. So a child that would be dropped never hides a larger set. A set that holds only
`:EqualAI` groups cannot change the first answer, by the measurement above; the builder
may skip that call for such a set.

### Why every version is reached

Call a version **irreducible** when it would be returned, each of its `:OnlyI` groups
binds in the inactive conformation, and making any one of its `:NonequalAI` groups
`:EqualAI` lowers the fitted count. A version that fails only the last condition holds a
`:NonequalAI` group that adds nothing: thermodynamics fixes its constants to the active
conformation's, or the inactive conformation never reaches the group's enzyme forms. It
has the rate equation of the version with that group `:EqualAI`. The two moves should
reach every irreducible version.

Take an irreducible version. Let `S` be its `:OnlyA` and `:OnlyI` groups and `N` its
`:NonequalAI` groups.

- **`S` is not empty.** `_enumerator_keeps` reads only `S`, so the version with `S` as it
  is and every other group `:EqualAI` is kept and is a kind (a) child. From there the
  groups of `N` are added in minimal gaining sets. What a group gains never falls as other
  groups become `:NonequalAI`: the constants a `:NonequalAI` group adds are new unknowns,
  the thermodynamic equations that bind them are fixed, and they bind less as more
  constants are free. So a
  set that gained when added still gains in the finished version, and while groups of `N`
  remain to add, some minimal gaining set of them exists.
- **`S` is empty.** `N` holds a minimal gaining set `T`. If `T` tagged `:OnlyA` would not
  be returned, `T` tagged `:NonequalAI` is a kind (b) child. If it would, that is a kind
  (a) child, and `_expand_change_allo_state` tags its groups `:NonequalAI`. The rest of
  `N` follows as above.

The premises are measured, not proved. `_enumerator_keeps` read only the `:OnlyA` groups
in 32,706 of 32,706 versions on the base tier. The argument about gains assumes that
`fitted_parameter_count` is the number of unknowns less the rank of the thermodynamic
equations. The third premise is that a kind (a) child of a minimal gaining set always
grows into that set's `:NonequalAI` version; it held on every mechanism measured. The
reachability test checks the first premise and the conclusion on every run (see "Tests").

The conclusion is measured too, from fitted counts and dedup keys at one subunit. On the
39 ping-pong mechanisms the two moves reach all 5,861 irreducible versions, which hold
4,574 distinct rate equations, and 52 further versions that each share a rate equation,
and so a fit, with one of them. On the 239 sequential mechanisms they reach all 14,999. At
two subunits, on the 14 ping-pong mechanisms that allow them, they reach all 2,214. On the
42 mechanisms with an inhibitor they reach all 14,898, of which 3,719 hold the inhibitor's
group `:OnlyI`.

Two cases were not measured. The argument takes the hyperbolic rule and the redundant
inhibitor to behave as `_enumerator_keeps` does, that is, to depend on `S` alone. None of
the 42 mechanisms has a redundant inhibitor in any version, and reach was not measured on
mechanisms with a reactant bound as an inhibitor, where 14 kind (a) children are
redundant. The reachability test runs on a sample of those too (see "Tests"). If it finds
an unreached version, the builder brings it to the maintainer.

### `_expand_add_dead_end_regulator(c, rxn)` on an allosteric parent

The move places the inhibitor on a set of enzyme forms, its placement. The maintainer
decided that on an allosteric parent the new group takes each state that adds one
constant. For each placement the move emits three children, the inhibitor's group
`:EqualAI`, `:OnlyA` and `:OnlyI`, where today it emits the `:EqualAI` one. It drops a
child in the two cases of `_expand_to_allosteric`: its `:OnlyI` group runs no step in the
inactive conformation, or it has a redundant inhibitor.

Where the inactive conformation has none of the forms the inhibitor binds, the `:EqualAI`
and the `:OnlyA` child have one rate equation and share a fit. Both stay: they part once
`_expand_change_allo_state` tags an `:OnlyA` group `:NonequalAI` and the inactive
conformation reaches such a form.

`_expand_to_allosteric` and this move then treat an inhibitor's group alike. Take a plain
mechanism, an inhibitor placement and a nonempty set of the mechanism's groups to tag
`:OnlyA`. Adding the inhibitor and then calling `_expand_to_allosteric` gives the same
three mechanisms as calling `_expand_to_allosteric` and then adding the inhibitor. A
version whose only group that is not `:EqualAI` is the inhibitor's comes from the first
order alone, since no allosteric mechanism has every group `:EqualAI`.

One detail of the code changes. `_dead_end_child` builds each child. For an allosteric
child it assumes that the active conformation runs every step, and calls
`_redundant_copy_groups` only when a cheaper test on the whole catalytic scheme says it
might matter. An `:OnlyI` group, in the parent or in the child, ends that assumption. The
move calls `_redundant_copy_groups` on every allosteric child; the builder may keep the
cheaper test where no group is `:OnlyI`.

### Other moves on a parent with an `:OnlyI` group

`_expand_re_to_ss`, `_expand_split_kinetic_group` and `_expand_add_allosteric_regulator`
run on such a parent as on any other.

`_assert_emission_rules` (`expansion_moves.jl`) checks every parent before its moves run.
One of its checks is that each steady-state kinetic group has a step that carries net
flux, read on the whole catalytic scheme. That stays right for a step the active
conformation runs: a cycle through an inhibitor's complexes follows copies of ordinary
steps, and replacing each copy by its original gives a cycle of the active conformation.
This is an argument, not a measurement. The function's docstring says that no move tags a
catalytic group `:OnlyI` and names the moves that filter their children; both sentences
are rewritten.

### The hyperbolic rule

A rate equation is hyperbolic when every substrate and product concentration appears to
at most the first power. Above one subunit the two-conformation model itself raises
concentrations to a power, and data cannot tell that from a catalytic scheme that does
the same. So above one subunit the enumerator returns only mechanisms in which each
conformation's rate equation is hyperbolic. `_hyperbolic_catalysis(c)` tests the two
conformations, and `_hyperbolic_where_required(c)` is true at one subunit and equals that
test above it (both in `candidate_mechanism_type.jl`).

Three moves filter their children with `_hyperbolic_where_required` today:

- `_expand_re_to_ss`, which can change the catalytic scheme's rate equation;
- `_expand_to_allosteric` and `_expand_change_allo_state`, which cannot change the scheme
  but can change one conformation's rate equation. With some groups `:OnlyA`, the inactive
  conformation runs part of the scheme, and that part can carry a power the whole scheme
  lacks. A ping-pong with one half-reaction `:OnlyA` gets `A²` in its inactive rate
  equation.

In the last two moves the filter removes nothing that `_enumerator_keeps` keeps today (see
"Guards and docs" for the counts), but only because `check_valid_rate_equation` refuses
those children for another reason.

The filter moves to the two places where `_enumerator_keeps` already filters what the
enumerator returns: `candidate_mechanisms(reaction, parent)` and `seed_mechanisms`, the
search for starting mechanisms described under "Starting mechanisms". The three moves lose
their filters, and the rule no longer rests on `check_valid_rate_equation`.

- **No child changes.** Each filter ran on the children a move had already built, and the
  other moves do not change the test's answer.
- **The two searches for minimal gaining sets still read the rule** (see "When a child
  would be returned"): tagging one more group `:NonequalAI` can make the inactive
  conformation's rate equation hyperbolic again, so a set that fails must be extended.
  `_expand_re_to_ss` has no such case: more steady-state steps never restore a hyperbolic
  rate equation, as its docstring says.
- **`_expand_to_allosteric` keeps its shortcut.** A parent that fails
  `_hyperbolic_catalysis` gets children at one subunit only, since above one subunit every
  child would fail.
- **`_base_tier` calls `_expand_re_to_ss`.** After the change under "Starting mechanisms"
  it does so on plain mechanisms only, which the rule does not concern.
- **Measured.** On today's starting mechanisms with a required regulator, none of the
  2,218 above one subunit fails the rule.

The maintainer chose this over removing the filter from `_expand_to_allosteric` and
`_expand_change_allo_state`. That would change no child today either, but the rule would
then rest on `check_valid_rate_equation` refusing those children.

### `_expand_re_to_ss(c)`: one rule for chain flanks

`_chain_flank_groups` loses its branch for allosteric mechanisms, so `_expand_re_to_ss`
skips a chain's flanks in every mechanism ("What was measured" defines chain and flank).
The chain is read on the catalytic scheme, whatever the states of its groups. For a plain
mechanism the docs prove that a flank at steady state adds no rate equation. For an
allosteric one the rule is a modeling choice, like the hyperbolic rule: most such children
add no constant that data determine, and the few that do leave the search (see "Known
costs"). After the change, the hyperbolic rule is the only check of the catalytic scheme
that asks whether a mechanism is allosteric.

### Starting mechanisms

When a reaction requires a regulator, which is the default for a declared one,
`candidate_mechanisms(reaction)` calls `seed_mechanisms`
(`src/Enumeration/starting_mechanisms.jl`). It applies `_expand_to_allosteric`,
`_expand_add_allosteric_regulator` and `_expand_add_dead_end_regulator` repeatedly, never
`_expand_change_allo_state`, and returns the mechanisms that bind every required
regulator. `_is_seed_node` decides which intermediate mechanisms it keeps expanding. Two
changes keep a required regulator from hiding versions.

**Rule 1 of `_is_seed_node`** changes from "no group or ligand is `:NonequalAI`" to "no
regulatory ligand is `:NonequalAI`". A starting mechanism's `:NonequalAI` kinetic groups
can then only be those of a kind (b) child. Without the change a required regulator would
hide every kind (b) child.

**`seed_mechanisms` starts from the base tier**, where it starts from `init_mechanisms`
today. Today it runs `_expand_to_allosteric` on the 7 degenerate ping-pong mechanisms, and
`candidate_mechanisms(reaction)` applies `_base_tier` to the allosteric results
afterwards. A version of one of the 21 base-tier ping-pong mechanisms that come from them
is then a starting mechanism only if `_enumerator_keeps` accepts the same states on the
degenerate mechanism. On the six-step ping-pong it refuses each version with one
half-reaction `:OnlyA`, and both halves `:NonequalAI`. With a regulator required no plain
mechanism is fitted, so those versions are unreachable. Measured with today's moves and
one required regulator at one and two subunits: the starting mechanisms grow from 6,428 to
7,126 and none is lost; the six-step ping-pong whose Q release is at steady state goes
from 12 versions to 80.

So `_base_tier` runs once on the result of `init_mechanisms`, before `seed_mechanisms`
expands anything, whatever the reaction requires, and not on the results of
`seed_mechanisms`. The three moves that `seed_mechanisms` applies leave every step at
rapid equilibrium or steady state as it was, so none of its results is degenerate. With a
required inhibitor and no regulator the starting mechanisms are the same 1,548 either way.

The maintainer approved the second change.

### What goes

`_tagging_bases`, `_valid_tagging` and `_stops_inactive_turnover` (helpers of
`_expand_to_allosteric`), the two branches of `_expand_to_allosteric` that build today's
children, and the two rules of `_expand_change_allo_state` for kinetic groups: one binding
group at a time with no test for a gain, and all chemistry groups together while no
binding group is `:OnlyA`. Also the `_hyperbolic_where_required` filter in three moves,
the allosteric branch of `_chain_flank_groups`, and the call of `_base_tier` on the
results of `seed_mechanisms`.

### A worked case

The six-step ping-pong with both chemistry steps at steady state has six kinetic groups:
the bindings of A and Q at free E, the bindings of B and P at the covalent form, and the
two half-reactions.

- **Kind (a) children:** 63 emitted. `_enumerator_keeps` keeps 42 on the measured
  mechanism of this shape.
- **Minimal gaining sets:** eight. Each half-reaction alone (two), and six pairs of
  binding groups: A with P, A with B, Q with B, Q with P, A with Q, B with P. A single
  binding group does not gain and is extended to its pairs.
- **Kind (b) children:** three of the eight, the pairs that `_enumerator_keeps` refuses
  when tagged `:OnlyA`: A with B, Q with P, B with P. The other five are kind (a)
  children, and `_expand_change_allo_state` tags each `:NonequalAI`. An `:OnlyA`
  half-reaction changes alone. A and P, both `:OnlyA`, change together, since A
  `:NonequalAI` beside P `:OnlyA` is refused.
- **Kind (c) children:** the three kind (a) children whose `:OnlyA` groups are only
  half-reactions (the first, the second, both), each with each regulator in two states.

## Tests

Written first and watched to fail, under the enumeration-engine rules of
`.claude/CLAUDE.md`.

- **`_expand_to_allosteric`, exact children**, on mechanisms small enough to write every
  child out: the uni-uni mechanism, the three-step ping-pong (three groups, seven `:OnlyA`
  sets) and one bi-bi sequential mechanism. At least one of them has a kind (b) child.
  One testset pins kind (c), with a reaction that declares a regulator.
- **Kind (b) children on the six-step ping-pong**, in addition to the exact-children
  testsets: the three sets above are among the children, and no single binding group and
  no set kept as `:OnlyA` is tagged `:NonequalAI` in any child.
- **`_expand_change_allo_state`, exact children**: from one half-reaction `:OnlyA` (that
  half tagged `:NonequalAI` among them), from A and P `:OnlyA` (both tagged `:NonequalAI`
  in one child), from binding groups `:NonequalAI` with chemistry `:EqualAI`, and from a
  parent whose every group is `:NonequalAI` (no child).
- **Reachability**, an aggregate check that says so in a comment: on the base tier, on a
  sample of the starting mechanisms with a competitive inhibitor, and on a sample of those
  with a reactant bound as an inhibitor, the two moves reach every irreducible version,
  found by trying every version; every child has more fitted constants than its parent;
  and the answer of `_enumerator_keeps` on a version equals its answer with every
  `:NonequalAI` group made `:EqualAI`. The builder samples the base tier too if the full
  set is slow. The test runs on the base tier, not on the degenerate mechanisms. The check
  that every child has more fitted constants than its parent also runs on a sample of
  mechanisms with two regulators, at separate sites and at one. That sample is the only
  test of the check on regulatory ligands, which no known mechanism makes fail.
- **`:OnlyI` inhibitor**, exact children: `_expand_to_allosteric` on the uni-uni mechanism
  with a competitive inhibitor bound at free E (four groups); a mechanism whose `:OnlyA`
  groups leave the inactive conformation without the inhibitor's only form, with no
  `:OnlyI` child for those sets; a mechanism with a reactant bound as an inhibitor whose
  `:OnlyI` child is dropped as redundant; `_expand_change_allo_state` tagging an `:OnlyI`
  inhibitor `:NonequalAI`; and `_expand_add_dead_end_regulator` adding a second inhibitor
  to a parent that holds an `:OnlyI` one. One mechanism with an `:OnlyI` inhibitor,
  written out, is compared with the mass-action steady state in `test/Derivation/`, unless
  the suite already holds one.
- **`_expand_add_dead_end_regulator` on an allosteric parent**, exact children: a small
  parent with one placement gives the three children; a parent whose `:OnlyA` groups
  leave the inactive conformation without the placement's forms gives the `:EqualAI` and
  the `:OnlyA` child alone. One testset takes a plain mechanism, one placement and one
  `:OnlyA` set, and checks that the two orders give the same three mechanisms.
- **Hyperbolic rule**: the testsets that pin a dropped child of `_expand_re_to_ss` pin it
  through `candidate_mechanisms(reaction, parent)`. The two that pin a dropped child of
  `_expand_to_allosteric` and `_expand_change_allo_state` cannot show the rule there,
  since `check_valid_rate_equation` refuses the same child. They assert that
  `_hyperbolic_where_required` refuses the written-out child, and that
  `candidate_mechanisms(reaction, parent)` does not return it.
- **Chain flanks**, exact children: `_expand_re_to_ss` on an allosteric parent with a
  chain, at one subunit and at two, has no child that changes a flank. Testsets that pin
  such a child get their new expected sets.
- **Starting mechanisms**: with a required regulator they hold, on the six-step ping-pong
  with both chemistry steps at steady state, a version with one half-reaction `:OnlyA` and
  the kind (b) child of A with B, each with the regulator bound. The testset
  "seed_mechanisms drops a seed the guard rejects"
  (`test/Enumeration/test_starting_mechanisms.jl`) pins a count and asserts that named
  mechanisms are starting mechanisms. The count is re-pinned. Each named mechanism is
  checked against the new starting mechanisms; one that is no longer among them goes to
  the maintainer before the assertion changes.
- **Derivation**: one ping-pong with binding groups `:NonequalAI` and chemistry
  `:EqualAI`, written out, compared with the mass-action steady state, in
  `test/Derivation/`.
- **Removed helpers**: the testsets in `test/test_oracles.jl` and
  `test/Enumeration/test_expansion_moves.jl` that call `_tagging_bases`, `_valid_tagging`
  or `_stops_inactive_turnover` are rewritten against the moves' children.

About 50 existing testsets pin the children of the two moves (29 for
`_expand_to_allosteric`, 21 for `_expand_change_allo_state`, in
`test/Enumeration/test_expansion_moves.jl`), and `test/test_oracles.jl` builds populations
from the children of `_expand_to_allosteric`. Their expected sets change. A mechanism with
six groups now has 63 kind (a) children per oligomeric state, too many to write out. The
maintainer approved how they are carried over:

- A testset that pins a removed rule (a result of `_tagging_bases`, the set of all
  chemistry groups, all chemistry groups changing together) is deleted and named in the
  pull request. The new tests of the rule that replaces it take its place.
- A testset about another aspect (the hyperbolic rule, an oligomeric state, one regulator)
  keeps its assertion and moves to a mechanism of at most four groups, or to one
  oligomeric state, so that every child is still written out.

A new expected set is reviewed against this design, not regenerated from the new output.

## Guards and docs

- **Allocation bound.** `test/Enumeration/test_starting_mechanisms.jl` bounds the bytes one
  expansion of the bi-bi starting mechanisms allocates. `_expand_to_allosteric` builds
  more children, so the bound will trip. The maintainer approved re-pinning it; the pull
  request states the bytes before and after.
- **Version.** 0.11.1: results change and no exported function does.
- **Engine page** (`docs/src/identify/enumeration_engine.md`):
  - Moves 4 and 6 are rewritten around sets of `:OnlyA` groups and minimal gaining sets,
    with the thermodynamic reason one binding group cannot be `:NonequalAI` alone and the
    three states of an inhibitor's group.
  - Move 3 gains the three states on an allosteric parent.
  - Move 1 says that chain flanks are skipped in allosteric mechanisms too, and that
    there the rule is a choice and not a proof.
  - The hyperbolic rule is stated once, for every move.
  - "Starting from the required regulators" says that the search starts from the base
    tier and that a kinetic group of a starting mechanism can be `:NonequalAI`.
- **`docs/src/deriving/re_vs_ss.md`**, "Apparent constants on a steady-state chain": the
  same sentence about allosteric mechanisms as under move 1.
- **Developer page** (`docs/src/developer.md`): the sentence on where the moves call
  `check_valid_rate_equation`, the sentence that three moves drop a child that fails the
  hyperbolic test, and the enumeration section where it names bases. Where a call between
  the files of `Enumeration` moves, `julia test/update_module_map.jl` regenerates the
  file table in the same commit.
- **MWC page** (`docs/src/deriving/mwc_allostery.md`): the sentence "as the search's
  allosteric bases do".
- **Docstrings** of every function this design changes, among them `seed_mechanisms`,
  `_dead_end_child`, `_assert_emission_rules`, `_chain_flank_groups` and
  `_expand_re_to_ss`.
- **`ROADMAP.md`**: the Known issue "An allosteric ping-pong whose covalent forms empty only
  by binding B cannot stop both half-reactions" shrinks to what was measured: ten starting
  mechanisms, no limit under any reading, fixed only by conformations that interconvert at
  the covalent form. Its three sub-entries go. `docs/src/roadmap.md` keeps its entry.
- **`ROADMAP.md`, other entries on these moves.**
  - Planned, "Drop `change_allo_state` relaxations that free no parameter": removed.
    `_expand_change_allo_state` keeps only children that raise the fitted count. The
    entry on search cost that lists "zero-gain `change_allo_state` children" as unmeasured
    loses that item.
  - Planned, "Decide whether to keep the hyperbolic drops above one subunit": removed.
    The rule stays and is applied where the enumerator returns its mechanisms. Measured
    here: every one of 13,678 kept base-tier versions on hyperbolic mechanisms passes
    `_hyperbolic_catalysis` at two subunits, and every one of 7,846 on mechanisms with
    inhibitors, but only because `check_valid_rate_equation` refuses the failing ones,
    which the Known issue below proposes to change.
  - Ideas, "Apply the chain-flank exclusion at one catalytic subunit": removed. Flanks are
    now skipped at every subunit count.
  - Ideas, new: replace allosteric mechanisms at one subunit by a regulator that binds
    forms of the catalytic scheme (see "Known costs").
  - Known issue, "The guard rejects some allosteric mechanisms whose law is the limit":
    its count of ter-bi starting mechanisms with no allosteric child is measured again.
  - The entry that names `_stops_inactive_turnover` as a caller of `has_turnover_cycle`
    loses that caller.
  - Known issues, new: one for each of items 1, 2, 4, 5 and 6 under "Known costs".
- **This document** and its plan are deleted in the pull request that implements them.

## Search cost

Children of `_expand_to_allosteric` per plain mechanism at one oligomeric state, without a
regulator, that `_enumerator_keeps` keeps (medians):

| Mechanisms | Today | Kind (a) | Kind (b) |
|:--|--:|--:|--:|
| Bi-bi sequential (239) | 7 | 11 | 2 |
| Bi-bi ping-pong (39) | 13 | 32 | 3 |
| Ter-bi (60 sampled of 5,158) | 12 | 23 | not measured |

- **Growth.** Kind (a) children about double with each added kinetic group: 13, 23, 53 and
  90 kept at four, five, six and seven groups on ter-bi. Today's children double with each
  binding group.
- **Time in `check_valid_rate_equation`.** All kind (a) children of the 60 ter-bi
  mechanisms took about a second.
- **Kind (c).** A sequential mechanism has one kind (a) child with only chemistry groups
  `:OnlyA` (median), a ping-pong mechanism three. Each adds two children per declared
  regulator.
- **`_expand_change_allo_state`.** Over the 5,861 ping-pong versions it emits about 16,000
  children, a median of three per mechanism.
- **Inhibitors.** A mechanism with one inhibitor keeps 25 kind (a) children without an
  `:OnlyI` group and 12 with one (medians of 42 mechanisms). With two inhibitors it keeps
  47 and 60 (61 mechanisms). Starting mechanisms that require an inhibitor and an
  allosteric regulator grow by about the same factor; that was not measured.
- **`_expand_add_dead_end_regulator`.** Its children on an allosteric parent triple: a
  median of 6 becomes 18 on today's children of `_expand_to_allosteric` with one inhibitor
  declared. All but the 55 pairs above are distinct rate equations.

## Known costs

Several items concern an **L phantom**. `rate_equation_dedup_key` holds a test, the
L-phantom certificate, that a one-subunit mechanism's `L` can be absorbed into its other
constants. A mechanism that passes is not fitted: it reads the fit of its rate equation at
`L = 0`. The test requires that the inactive conformation turn nothing over and that each
term of its denominator be a term of the active conformation's.

1. **Versions that data cannot determine at one subunit.** A set gains by the fitted
   count. On the six-step ping-pong, A with Q, and B with P, tagged `:NonequalAI`, fit two
   constants more than the plain mechanism; at one subunit data determine no more
   constants than the plain mechanism's, and at two subunits one more. The certificate
   does not cover them, since their inactive conformation turns over, so they are fitted.
   How many versions of the base tier this covers was not counted.
2. **`:OnlyI` inhibitors that data cannot determine at one subunit.** Of 122 sampled
   kind (a) children with the inhibitor's group `:OnlyI`, 47 fit one constant more than
   data determine at one subunit. With no `:OnlyA` group, in 4 of 5 sampled, the rate
   equation is the plain mechanism's with its inhibitor constant replaced by
   `K·(1 + L)/L`. The certificate covers 4 of the 47: the inhibitor's terms are in the
   inactive conformation's denominator alone. So 43 are fitted with a constant data cannot
   determine. The same 122 with the inhibitor's group `:OnlyA`: 28, of which 19 pass the
   certificate. At two subunits 3 of 122 are undetermined.
3. **Optional regulators on kind (a) children with an `:OnlyA` binding group.** These
   get no kind (c) child. At one subunit an L phantom among them has its parent's loss, so
   it enters the beam of `identify_rate_equation` only through the beam's minimum width,
   as today's children do.
4. **`:NonequalAI` versions reached only through an L phantom.** Where a minimal gaining
   set is kept as `:OnlyA`, its `:NonequalAI` version is reached only through that kind
   (a) child. At one subunit that child is an L phantom for 227 of the 539 such sets on
   the sequential mechanisms and 11 of the 124 on the ping-pong ones. Above one subunit
   there are no L phantoms. Emitting every minimal gaining set from `_expand_to_allosteric`
   would remove this cost (the bottom row of the table under "Approaches").
5. **Flank changes that show.** In the sample, 20 of 96 children that change only flanks
   to steady state at two subunits, and 2 of 128 at one, add a constant that data
   determine. `_expand_re_to_ss` no longer makes them. So a slow binding beside a slow
   isomerization is missing from the allosteric search where the two conformations treat
   the chain differently.
6. **A steady-state group that runs only beside an `:OnlyI` inhibitor.**
   `_expand_split_kinetic_group` can give a step between two of the inhibitor's complexes
   a group of its own, and `_expand_re_to_ss` can then change it to steady state. With the
   inhibitor's group `:OnlyI` that step runs in the inactive conformation alone, and where
   that conformation turns nothing over its two rate constants enter as a ratio. The flux
   check of `_assert_emission_rules` reads the whole scheme and does not see it.
7. **One subunit.** The items above that name one subunit share a cause: there the
   equilibrium between the conformations only reweights the enzyme forms, so the rate
   equation is a plain one written with `L` and two sets of constants. The maintainer
   keeps one subunit. It is the only way a monomer binds an allosteric regulator, and the
   only allosteric version of a scheme that is not hyperbolic. Its replacement would be a
   regulator that binds forms of the catalytic scheme and leaves them active with
   constants of their own, a family of mechanisms that needs a design of its own.
8. **Selections can change** on data already analysed, since more mechanisms compete.

## Interaction with the iso-enumeration design

`docs/superpowers/specs/2026-11-24-iso-enumeration-design.md` (branch `iso-enumeration`)
adds a step that isomerizes the free enzyme. Three of its statements meet this design:

- It states its allosteric behaviour in terms of `_tagging_bases`. Here the same decisions
  read: that step's group is not a chemistry group, so a kind (a) child that tags it
  `:OnlyA` gets no kind (c) child, and `_expand_change_allo_state` tags it `:NonequalAI`
  like any group.
- It relies on `_expand_to_allosteric` dropping a parent that is not hyperbolic above one
  subunit. That shortcut stays.
- It describes the order of the starting mechanisms through `seed_mechanisms` and
  `_base_tier`, and says that a degenerate ping-pong mechanism stays in the starting set
  to be expanded. Here `_base_tier` runs before `seed_mechanisms`.

Whichever branch lands second rewrites those passages and rebases `expansion_moves.jl` and
`starting_mechanisms.jl`.

## Left to the builder

- How the two moves share the search for minimal gaining sets, and whether extending a
  set by "any other group" needs a narrower rule once measured on ter-ter mechanisms.
- The order of the three tests under "When a child would be returned" and of the fitted
  count, by cost.
- Whether `seed_mechanisms` calls `_base_tier` itself or receives the base tier from
  `candidate_mechanisms(reaction)`.
- Where in `seed_mechanisms` the hyperbolic rule runs: on each intermediate mechanism or
  on the results.
- Whether `_expand_to_allosteric` and `_expand_add_dead_end_regulator` share the code that
  gives an inhibitor's group its three states and drops the two cases.
- The order of the children, as long as it does not depend on the number of workers.
- The mechanisms for the exact-children tests, within the test rules.
- A recovery check during implementation: data generated from a ping-pong with binding
  groups `:NonequalAI` and chemistry `:EqualAI`, and from a mechanism with all chemistry
  `:OnlyA` and a regulator, each recovered by `identify_rate_equation`. It need not become
  a suite test.

## Decided by the maintainer

- Approach A.
- `_expand_to_allosteric` emits kinds (a), (b) and (c) as described.
- `_expand_change_allo_state` tags groups `:NonequalAI` from `:EqualAI`, `:OnlyA` and
  `:OnlyI`, and keeps a child only when its fitted count rises, for kinetic groups and
  regulatory ligands.
- A competitive inhibitor's group can be `:OnlyI`.
- `_expand_add_dead_end_regulator` gives its inhibitor three states on an allosteric
  parent.
- `_expand_re_to_ss` skips chain flanks in allosteric mechanisms too.
- The hyperbolic rule is applied where the enumerator returns its mechanisms.
- `seed_mechanisms` starts from the base tier.
- Allosteric mechanisms at one subunit stay.
- `check_valid_rate_equation` is left for another change.
- How the existing testsets are carried over, the re-pinned allocation bound, and version
  0.11.1.

One point the maintainer has not seen: the two searches for minimal gaining sets read the
hyperbolic rule and the redundant-inhibitor test, so that a child they would drop never
hides a larger set (see "When a child would be returned"). It follows from the decision
above and changes no child today.
