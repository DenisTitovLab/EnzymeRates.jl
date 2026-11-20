# One lean `Mechanism`: design

Date: 2026-11-18. Status: approved by the maintainer on 2026-11-19. Stage 1, the guard's
rule (Sequence, item 1), is built in 0.11.0; the later stages are not.

## Why

Too many names cross between `Mechanisms`, `Enumeration` and `Derivation`. `Mechanisms`
has 25 public names, four mechanism types, an abstract type, conversions both ways and a
set of stored facts that only some callers read. A reader who wants to change the
derivation must learn enzyme forms, residuals, conformation labels, step kinds and the
integer encoding of a compiled type, although the rate equation depends on none of them.

The rate equation depends only on a reaction graph: which forms exist, which steps join
them, what each step takes up and gives off, which steps are fast and which share
constants. This design makes that graph the whole mechanism.

## Words used here

- **Conformation:** the active or the inactive shape of the enzyme (MWC).
- **Allosteric state** of a kinetic group: one of `:OnlyA`, `:OnlyI`, `:EqualAI`,
  `:NonequalAI`, written after `::` in the macro.
- **Excluded step:** a step that does not run in a conformation because of its group's
  allosteric state: an `:OnlyA` step in the inactive conformation, an `:OnlyI` step in
  the active one.
- **Cut-off form:** a form that the steps still running in a conformation do not join to
  free E.

## Goals

1. One mechanism type, `Mechanism`, a plain value, for plain and allosteric mechanisms.
2. The mechanism stores a lean reaction graph and nothing derived from it.
3. `Derivation` reads only that graph: no form contents and no step kinds.
4. Form names are opaque Symbols.
5. The guard still admits a mechanism only when its law is the unique limit of
   mass-action kinetics, and its rule for excluded steps is read from the graph alone.

## Non-goals

Each of these gets its own design (see "Sequence"):

- What Enumeration needs to know about a form, and whether it can read it from the graph.
- Derivation's interface to its callers: identify's questions, an Identification module,
  a units function for Fitting, and whether Enumeration should ask `Derivation` about the
  rate law instead of predicting it from the graph.
- Dropping `E_total` and reshaping the parameter kinds.

Also outside this design: new enumeration moves, iso mechanisms, a flip at a covalent
form, identify's keyword arguments.

## Anti-goals

- Old and new paths side by side. Every caller migrates in the change that removes the old
  name.
- A reader for the old macro notation or for mechanism text saved by earlier versions.
- Code specialized per mechanism outside the compiled rate equation.
- A weakened oracle test (see "Tests").

## Constraints

- The compiled rate equation stays allocation-free and under
  `50 ns + 0.5 ns × denominator terms`; `Fitting.loss!` keeps its bounds. The compile-budget
  and type-stability tests keep passing.
- One Julia process at a time; the full suite needs about 4 GB.
- Adding or changing a public name needs the maintainer's approval. This document is that
  request for the names it lists.
- The module map in `.claude/CLAUDE.md` and `docs/src/developer.md` is regenerated with
  `julia test/update_module_map.jl` in the commit that changes a `public` statement.

## Approaches considered

### What the mechanism stores

1. **Lean graph only (chosen).** The one thing `Derivation` reads from form contents
   today is whether an excluded step is a plain binding. A rule read from the graph
   replaces that (decision 2).
2. **Lean graph plus one Bool per step** ("this uptake also does the chemistry"). Keeps
   today's verdicts exactly. Rejected: the Bool encodes a naming convention, not a fact.
   A merged complex is one node that stands for both `E(A, B)` and `E(P, Q)`, and today's
   verdict depends on which name the enumerator gives it.
3. **Keep contents in `Mechanisms`.** Rejected: it is today's design.

### How an excluded step is read

1. **By form contents (today).** A plain binding loses its affinity; anything else is
   stopped. Needs the contents.
2. **Every one-metabolite uptake between separated forms loses its affinity.** Rejected:
   on the test population it rejects 328 mechanisms admitted today, where the bound form
   is still joined to free E.
3. **An uptake into a cut-off form loses its affinity; every other excluded step is
   stopped.** Loses nothing admitted today and adds 2,352, all with laws equal to the
   numerical limit. Rejected as it stands: it lets a rapid-equilibrium step be the slowest
   step in the other conformation, which today's guard also does for a step that is not a
   plain binding.
4. **As 3, and a step keeps its kind in both conformations (chosen).** The maintainer's
   principle: rapid equilibrium and steady state differ more than 100-fold, and a step
   does not jump from the fastest kind in one conformation to slower than every
   steady-state step in the other. A steady-state uptake into a cut-off form keeps its
   kind either way, so the guard may read it as stopped too, and tries three fixed
   choices (added after the gate and the enumeration tests, on the maintainer's
   decisions). Built, and kept as the first three readings of approach 5.
5. **Energy levels (chosen).** Every cut-off group of forms lies higher in energy; an
   excluded step between two groups loses its uphill rate in either direction, or stops
   if steady-state; the guard searches for a reading that passes its checks (decision 2).
   It admits everything approach 4 does, most allosteric ping-pongs that stop both
   half-reactions, and about a thousand mechanisms more in the populations measured. It
   needs one face check approach 4 does not. The maintainer chose it after the spikes,
   knowing the evidence for that check is numerical.

### How a mechanism is compiled

1. **`Mechanism{Sig}` with fields, `Sig = Nothing` until compiled.** One type. Costs a
   rule at every public entry point: drop `Sig` before internal work.
2. **Plain `Mechanism` and a separate compiled rate equation (chosen).** The pattern of
   ModelingToolkit and Catalyst: the model is a value, and code is generated when a problem
   is built. One concrete type everywhere inside the package. Chosen because the maintainer
   expects no user to need a fast rate function outside a fit: a user who wants the law in
   another model copies the printed equation.
3. **Typed only.** Rejected on measurement: a type costs 1.2 KB that is never freed (a
   value costs 1.8 KB and is freed), a type parameter cannot hold vectors, so every read
   decodes (72 µs today), and Enumeration would keep a value form anyway.
4. **Interpret the law, never compile** (as SymbolicRegression does). Rejected for this
   refactor: it rewrites the code generator.

## Decisions

### 1. The data model

All in `Mechanisms`. Field names are part of the interface and are read with dots.

```julia
struct Step
    from::Int                         # index into enzyme_forms
    to::Int
    taken_up::Vector{Int}             # indices into metabolites
    given_off::Vector{Int}
    rapid_equilibrium::Bool
    kinetic_group::Int                # steps of one kinetic group share their constants
end

struct CatalyticScheme
    metabolites::Vector{Symbol}
    enzyme_forms::Vector{Symbol}      # opaque names; the free enzyme first
    substrates::Vector{Int}           # indices into metabolites
    products::Vector{Int}
    steps::Vector{Step}
end

struct Allostery
    allosteric_states::Vector{Symbol} # one per kinetic group; empty for a plain mechanism
    catalytic_multiplicity::Int       # 1 for a plain mechanism
    regulatory_sites::Vector{RegulatorySite}
end

struct Mechanism
    catalytic::CatalyticScheme
    allosteric::Allostery
end
```

- A mechanism is plain when `allosteric.allosteric_states` is empty.
- `RegulatorySite` names its ligands' states `allosteric_states` too, in place of
  `mwc_tags`: both hold the same four values.
- A metabolite that is neither a substrate nor a product is a modifier. A reactant that
  binds as a modifier (substrate inhibition) is a step that takes up that reactant; the
  graph has no separate modifier copy.
- The free enzyme is the first form. Nothing inside the package reads its name; the macro
  decides which written form it is (decision 4).
- Two writings of one mechanism (step order, written direction, group order, order of the
  forms) are `==` and hash alike. Form names are part of a mechanism's identity. A
  one-metabolite step is stored as the uptake.
- The constructor checks representation only: indices in range, every form joined to the
  free enzyme, distinct names, uniform kinetic groups, one known allosteric state per
  group. It does not run the guard, and it has no rule about which groups may be `:OnlyI`
  (decision 2).

Removed: `AllostericMechanism`, `EnzymeMechanism`, `AllostericEnzymeMechanism`,
`AbstractEnzymeMechanism`, `EnzymeForm`, `Metabolite`, `catalytic`, `compile_mechanism`,
`check_conformation_label`, `excludes`, the stored facts (`modifiers`, `groups`,
`step_kinds`, `toward_products`, `step_sides`, the four step masks) and the integer
encoding of names.

### 2. The guard's rule for excluded steps

#### What a conformation keeps

In a conformation that excludes steps:

1. An excluded step does not run.
2. A form stays only if the steps that still run join it to free E. The other forms are
   cut off and hold no enzyme of that conformation.
3. A step between cut-off forms does not run either, whatever its allosteric state.

The conformation's law is the law of what is left. Example, in the inactive conformation:

```
E + S ⇌ ES      :: OnlyA      excluded
ES <--> EP      :: OnlyA      excluded
E + P ⇌ EP      :: EqualAI    runs; E and EP stay, ES is cut off
```

#### How an excluded step is switched off

The guard asks whether that law is still the limit when the excluded steps are nearly
absent, not absent. The science it assumes: an excluded step has a very high or very low
equilibrium constant, or very low rate constants, and every thermodynamic constraint still
holds in each conformation. The guard reads this from the graph, not from the form
contents, and a step keeps its kind (rapid equilibrium or steady state) in both
conformations.

The picture is energy. The steps a conformation keeps join the forms into groups. Free E's
group holds the law's forms. Every other group is cut off: in that conformation it lies
higher in energy, and so holds no enzyme. Then:

| Excluded step | How it is switched off |
|:--|:--|
| Between two groups | It loses its uphill rate and keeps its downhill rate: the far side is higher. A rapid-equilibrium step stays fast. |
| Between two groups, steady state | It may be stopped instead: both rates vanish. |
| Inside one group, steady state | It is stopped. |
| Inside one group, rapid equilibrium | Fault. Its two sides have one energy, so thermodynamics pins its equilibrium constant, and it could be switched off only by turning slow. |

Which group lies above which is not given. A **reading** chooses, for each excluded step
between two groups, which side is higher, or, for a steady-state step, that it stops. A
reading is consistent when the steps that lose, each pointing uphill:

1. close no cycle, so energies exist;
2. never point into free E's group, since nothing lies below it; and
3. lead from free E's group to every group, so every cut-off group lies above it and has a
   step by which it empties.

The law is the same under every reading, since the law drops the excluded steps. The guard
admits the mechanism when some consistent reading passes the face checks below.

Why a step loses its uphill rate: the law assumes the cut-off forms are empty, and a step
into a higher group that kept its uphill rate would fill it. "Loses its affinity" is the
case in which the step is a binding and its bound side is the higher one.

Why a rapid-equilibrium step may lose in either direction but never stop: the maintainer's
principle. Rapid equilibrium and steady state differ more than 100-fold, and turning a
rapid-equilibrium step into the slowest step of the other conformation is not a reasonable
approximation. A mechanism that needs it is written with that step steady-state.

Why thermodynamics holds: each group has one energy level, and every excluded step's
equilibrium constant shifts by the difference between its two ends' levels. The shifts
cancel around every cycle.

Four examples, each in the inactive conformation, the unlisted steps `:EqualAI`:

| Mechanism | Verdict |
|:--|:--|
| `E + S ⇌ ES :: OnlyA`, `ES <--> EP :: OnlyA`, `E + P ⇌ EP` | Admitted. `ES` is cut off: S binding loses its uptake, and the chemistry loses its rate into `ES` or stops. |
| `E + S ⇌ ES :: OnlyA`, `ES <--> EP`, `E + P ⇌ EP` | Fault. `ES` stays in free E's group, made from P, so S binding is inside one group and a rapid-equilibrium step cannot be stopped. |
| `E + S <--> ES :: OnlyA`, `ES <--> EP`, `E + P ⇌ EP` | Admitted. The steady-state binding is inside one group and is stopped. Enumeration does not build this one (decision 5). |
| Ping-pong `E + A ⇌ EA`, `EA <--> FP :: OnlyA`, `F + P ⇌ FP`, `F + B ⇌ FB`, `FB ⇌ EQ :: OnlyA`, `E + Q ⇌ EQ` (F the modified enzyme) | Admitted. `F`, `FP` and `FB` are one cut-off group: `EA → FP` and `EQ → FB` vanish, their reverses keep their rates, and the group empties through both halves. |

#### The search for a reading

- The guard first tries the three readings built in stage 1's earlier tasks (Evidence
  describes them), in which only a step that takes up one metabolite into a cut-off form
  loses, toward its bound side. They are cheap and decide most mechanisms, and trying
  them first keeps every mechanism the built guard admits.
- If none passes, it searches the other consistent readings: each conformation's choices
  against the three graph conditions above, then the surviving pairs, one choice per
  conformation, against the face checks, in a fixed order, stopping at the first that
  passes. Every reading, the three built ones included, gets the same face checks.
- The search is skipped, whole, for a mechanism that is too large: when a conformation
  has more choices than one internal constant, or the pairs that survive the graph
  conditions number more than a second (the spike's were 20,000 and 5,000, which kept the
  slowest guard call near 0.1 s; the builder sets them from the timings). It never stops
  part-way, so a verdict does not depend on the order of the steps, and a mechanism and
  its mirror image get the same one. A mechanism too large to search that the three built
  readings reject is rejected. It is never wrongly admitted.
- The guard can report, to the tests, the reading it admitted with.
- A rejection's message describes a reading the guard tried and never says that a step
  cannot be switched off when some reading switches it off. A rapid-equilibrium step
  inside one group is reported first, since no reading changes it.

#### The face checks

With a reading chosen, on every face (the interior, and each set of metabolites at zero
that the guard checks today) class 3 and admission's face condition run as built. They
read only which transitions vanish and which are fast, so they take any reading. Two
things are added.

**No form fills faster than it empties (needed for soundness).** In the built readings the
direction a losing step keeps is always a release, which runs whatever is present. Under
this rule a step can keep its uptake, which dies on the face where its metabolite is zero.
Its vanishing direction then runs one way into the higher form. The check: on a face, a
rapid-equilibrium transition that vanishes and is alive while its reverse, which does not
vanish, is dead is refused when

- the enzyme reaches its source on that face, and alive rapid-equilibrium transitions
  (vanishing or not) lead from the law's forms to that source, so the flow into it grows
  with the speed of the fast steps; and
- no alive rapid-equilibrium transitions that do not vanish lead from its target to a law
  form, so the target empties only through a slow step.

The argument: a form's occupancy is its inflow over its outflow. Past a slow transition
the inflow is bounded by a rate that vanishes. Through fast transitions alone it grows
with their speed, and then only an outflow that grows as fast keeps the form empty in
every order of the two limits. A steady-state transition filling one way is covered by the
built stranded-form check. This argument is not a proof that the face condition with the
check is sufficient for every reading; the evidence below is numerical, and the gate
(Tests) is built to be hard to pass by accident.

**Groups that cannot matter are set aside (removes over-cautious rejections).** The face
condition refuses a mechanism when vanishing transitions open a faster route between two
of the law's forms and it cannot prove the route harmless. When it refuses a face as a
bridge, the guard runs it once more on that face with some groups set aside, and the face
passes if either run does. (A face refused for a stranded form is not run again: the
stranded form's group cannot drain, so it is never set aside, and the second run would
refuse the face again.) A group is the cut-off forms the enzyme reaches on the face that
transitions join to one another, with every alive transition at them. It is set aside
when, on this face:

- the law's forms it touches lie within one set that the law's own rapid-equilibrium
  steps, alive in both directions, hold in equilibrium (one touching form counts);
- its forms drain to the law's forms however fast the rapid-equilibrium steps are. Take
  the group's alive transitions that do not vanish. Its fast classes are the sets of its
  forms that the fast ones among them cannot leave (a form with no fast way out is a class
  of its own). Every fast class must lead to a law form by slow transitions out of fast
  classes, each followed along fast transitions to the classes or law forms they reach. A
  slow transition out of a form that has a fast way out does not count: its share vanishes
  as the fast steps speed up. This is the chain class 3 builds on the law's network,
  built on the group; and
- each of its transition pairs is safe: it runs both ways; or only its direction that does
  not vanish is alive, so it only empties its form; or it joins two cut-off forms, only
  its vanishing direction is alive, and that direction is slow, or its source is not one
  that alive rapid-equilibrium transitions reach from the law's forms, or its target
  reaches a law form along alive rapid-equilibrium transitions that do not vanish. A pair
  at a law form with only its vanishing direction alive is not safe.

A cut-off form the enzyme never reaches on the face, which empties only into a group set
aside, goes with it: left alone it would read as a second place the enzyme could settle.

The argument: a second route between two forms already in equilibrium changes nothing,
and with one touching form there is no route between two different forms. The group's own
occupancy vanishes because every form drains at a rate that does not shrink as the fast
steps speed up, and none is filled faster than it empties. Class 2 has already refused any cycle of fast transitions that turns the reaction over, so
a cycle inside the group carries a flux that vanishes with its occupancy.

A wider version was measured and is not built: a group touching forms that are in
equilibrium for another reason (a dead-end branch, or a part of the law's network with no
cycle that turns the reaction over), in detailed balance with the law's steady state. It
admitted nothing in any population.

#### Against the built guard

| Case | Built guard (three tries) | This rule |
|:--|:--|:--|
| Rapid-equilibrium step inside one group | Fault | Fault |
| Steady-state step inside one group | Stopped | Stopped |
| Step that takes up one metabolite into a cut-off form | Loses its affinity; steady state may stop | The same, as one reading among several |
| Rapid-equilibrium step between two groups that is not such an uptake | Fault | Loses its uphill rate, in either direction |
| Steady-state step between two groups that is not such an uptake | Stopped | Loses its uphill rate, in either direction, or stopped |
| A direction kept on a face where its metabolite is zero | Cannot arise | The filling check |
| Cut-off forms touching the law's forms within one equilibrated set | Face condition may refuse | Set aside |

Nothing the built guard admits is rejected. Its readings are tried first; the filling
check never fires on them, since the direction they keep is always a release; and setting
groups aside can only pass a face that failed.

#### `:OnlyI` on any group

Today the constructor accepts `:OnlyI` only on a group that binds modifiers alone, which
it reads from form contents. That check is removed and nothing replaces it (the
maintainer's decision): any group may be `:OnlyI` as long as the active conformation
still catalyzes the reaction, which the guard already checks (`_NO_ACTIVE_TURNOVER`). A
mechanism with two reaction paths may give one to the inactive conformation alone.

The two states are then treated alike but for that one rule. A mechanism and its mirror
image (every `:OnlyA` and `:OnlyI` swapped) can both be admitted when both conformations
catalyze; they are one model with `L` inverted. Only hand-written mechanisms meet this:
no enumeration move makes a chemistry step `:OnlyI`.

#### Evidence

From throwaway spikes that patched today's guard and the test suite's oracles in a
running session. "Law = limit" is the suite's own numerical test on every face: the
ε-path oracle (the two-conformation mass-action network with every excluded step scaled
by powers of ε, along six directions and five orders of the two limits) and the RE-speed
oracle against the package's law.

| Population | Admitted today | Of those, rejected by the rule | Newly admitted | Newly admitted, law = limit |
|:--|--:|--:|--:|:--|
| Every way of making groups `:OnlyA` in the plain init mechanisms of uni-uni, bi-bi and ping-pong (9,616) | 6,748 | 258 | 1,030 | all 1,030 |
| The 2,584 distinct allosteric mechanisms `test/Derivation/test_admission.jl` puts through the guard | 2,224 | 151 | 24 | all 24 |
| 42,475 sampled mechanisms with an `:OnlyI` group (bi-bi and ping-pong, init and one expansion round) | 480 | 30 | 1,602 | 400 of 400 sampled |

- Every mechanism the rule rejects though it is admitted today has an excluded
  rapid-equilibrium step that would have to stop. In the first two populations that step
  binds and reacts at once in 350 of the 409, and is an isomerization in nearly all the
  rest. They leave on purpose.
- Sampled from the mechanisms admitted today and by the rule, the law equals the limit
  under the rule's readings: 200 of 200.
- With the excluded steps absent rather than slow, the law equals the steady state for
  9,567 of the 9,616 in the first population; the other 49 have a rapid-equilibrium fault
  that every rule rejects.
- Eight hand-written mechanisms with two reaction paths, one of them `:OnlyI`: the six
  whose excluded steps are steady-state (two isomerization paths; an exchange step beside
  a ternary-complex path) are admitted with laws equal to the limit, and so are their
  mirror images. The two written with rapid-equilibrium bindings are rejected with their
  mirrors, as the rule says.
- 746 sampled mirror pairs in which both conformations catalyze get the same verdict in
  every pair.
- The enumerator's output for the census reactions, with its own filter (decision 5):
  uni-uni unchanged; bi-bi loses 12 of 2,060 children; ping-pong loses 38 of 1,593 in round
  one and 316 of 7,014 in round two, and gains 4 and 17.

The gate (see Tests) reran these on the guard as first built, with a sample of 2,000
ter-ter taggings added. That guard lets only a step that takes up one metabolite into a
cut-off form lose, toward its bound side, and reads the steady-state ones in three fixed
ways, tried in order: all lose their affinity; only those that lead away from free E do
(their bound side's group is more such steps from free E's than their unbound side's);
none does. "The new fault" below is a rapid-equilibrium step that guard cannot switch off.

- The first population came out exactly as in the table. In the `:OnlyI` sample every
  lost mechanism has the new fault, and 400 of 400 sampled newly admitted ones have laws
  equal to the limit.
- Ter-ter, 1,202 admitted before: 26 lost, 169 newly admitted, all 169 with laws equal to
  the limit. Of the 26, 22 have the new fault and 4 are rejected as a bridge. Each of the
  4 holds an excluded step that binds and reacts at once into a cut-off form, which the
  old guard read as stopped.
- In 2 of the 4 that step is steady-state. Reading every such step as stopped (the third
  of those readings) admits both, and one ter-ter mechanism more, each with a law equal to
  the limit on every face. It changes no verdict in the other two populations.
- In the other 2 that step is rapid equilibrium, so the old guard admitted them only by
  turning a fast step slow. They stay rejected, as a bridge and not with the new fault
  (the maintainer's decision). For one of them the numerical limit does depend on the
  path.
- With the first and third readings alone, the enumeration tests then found 64
  mechanisms, admitted before and outside the gate's populations, that the guard rejected:
  of the 9,878 K-type taggings of the bi-bi split children, the ones in which one excluded
  steady-state binding must lose its affinity and another must stop. The second reading
  admits all 64, and 5 `:OnlyI` and 2 ter-ter mechanisms more, each of the 71 with a law
  equal to the limit on every face (the maintainer's decision to add it).
- No single reading replaces the three. The second alone rejects 81 ter-ter mechanisms the
  first admits; the first and second without the third lose one of the two steady-state
  ter-ter mechanisms.
- With the three tries, every mechanism admitted before and now rejected, in all four
  populations, has the new fault or is one of the two rapid-equilibrium ter-ter ones.
- Trying every mixture of losing and stopping for those steps admits 2 mechanisms more of
  the 65,000, both ter-ter and both rejected before this work (entries 54057 and 54341 of
  the gate's populations), each with a law equal to the limit. The search of decision 2
  covers them.

The energy-level rule was measured in further spikes (November 2026), which fed readings
to the built guard's own class-2 and face checks and, separately, to the suite's oracles.
The harsh test is the ε-path oracle with corner directions; the mild one, with standard
directions, misses the failure below.

| Population | Admitted by the built guard | More under the rule | Of those, set aside rule | Oracles agree |
|:--|--:|--:|--:|:--|
| 25 ping-pong init mechanisms with a covalent form, chemistry `:OnlyA` and every step `:OnlyA` (50) | 0 | 40 | 0 | 40 of 40 (38 harsh) |
| `:OnlyA` taggings of uni-uni, bi-bi and ping-pong inits (9,616) | 7,520 | 180 | 35 | 180 of 180 (154 harsh) |
| Sampled taggings with an `:OnlyI` group (43,311) | 2,200 | 34 | 0 | 34 of 34 (10 harsh) |
| Sampled ter-ter taggings (2,000) | 1,350 | 254 | 23 | 143 of 143 checked (10 harsh) |
| K-type taggings of the bi-bi split children (9,878) | 7,474 | 600 | 0 | 120 of 120 checked (30 harsh) |

- Without the filling check the rule admits laws that are not the limit: at least 14 of
  the 253 it then admits in the second population, each with a limit at zero B that
  depends on how the steps vanish. The mild test passed them; the harsh test did not.
- A first form of the filling check, which refused every one-way fast filling, was sound
  but refused 13 of 42 ping-pongs, most of them valid: where the source is fed through a
  slow step the flow vanishes. The form above is the second.
- Of the three ping-pongs with both chemistry steps `:OnlyA` that stay rejected, none has
  a limit under any reading: their covalent forms empty only by binding B.
- None of the mechanisms the guard as first built dropped on purpose comes back (the 542
  of the gate's first four populations that were admitted before this work and held a
  rapid-equilibrium step the three readings could not switch off). In 126 that step lies
  inside one group. In the other 416 it lies between two groups, so the rule lets it lose
  a direction: 392 are bindings that react at once, which lose the release, 23 are
  isomerizations, and the face checks refuse every consistent reading (383 by the filling
  check, 33 as a bridge). They refuse rightly: of 80 sampled, none has a limit under any
  consistent reading. Other mechanisms with such a step between two groups, never admitted
  before, are admitted now, the six-step ping-pong among them. An observation, reported by
  the gate, not a condition.
- One fixed reading (levels by breadth-first distance) admits 6,252 of the 7,520, so a
  search is needed. Every reading cost 0.15 to 1.2 s per population and at most 91 ms for
  one mechanism; 53 of 644 ter-ter mechanisms searched exceeded the spike's cap.
- What stays rejected though valid, from 120 sampled of the 832 mechanisms the second
  population still loses at the face checks: 2 as a bridge; and 11 to 16 by the filling
  check whose rate equals the law's although the cut-off forms hold most of the enzyme on
  that face (the inactive conformation's kept chemistry has the active one's constants).
  No simple rule found separates the last kind from the wrong ones; they stay rejected.
- The set-aside rule's drain condition was first weaker: each form of the group had only
  to reach a law form along alive transitions that do not vanish. Building it showed that
  form unsound. In four hand-written test mechanisms (the R2funnel family of
  `test/Derivation/test_admission.jl` and the enumeration tests), read with one step
  stopped, the group's only such way out leaves a form that a fast step empties elsewhere,
  so the group drains ever more slowly as the fast steps speed up; the harsh test finds a
  limit at zero P that depends on the direction. The populations above held no such
  mechanism. With the condition as stated above the four are no longer set aside under
  that reading, and the search admits each under another reading (the stopped step loses
  its uphill rate instead), with a law equal to the limit.
- The 673 distinct allosteric mechanisms written in the test files, 201 of them open to a
  reading: the search with the filling check admits 55 more than the built guard, oracles
  agree on 55 of 55 (51 harsh). The set-aside rule adds 10, oracles agree on 65 of 65 (60
  harsh). It still adds the same 35 and 23 in the second and fourth populations of the
  table, oracles agree on 35 of 35 (25 harsh) and 23 of 23 (mild).
- The safety condition decides alone on 215 faces of 40 small mechanisms under searched
  readings, and on none under the built ones. Each is a stopped steady-state step at a law
  form that runs one way along its vanishing direction because its uptake is dead on the
  face. Of the twelve smallest, the harsh test finds a limit that depends on the direction
  in 4 and a limit equal to the law in 8, which stay rejected.

The gate ran on the guard as built, against the guard as first built with its three
readings. It added five populations to the five of the table: `:OnlyA` taggings of bi-bi
with a competitive inhibitor, of ter-bi and of bi-ter (3,000 each), the allosteric
children of one expansion round of the census reactions (1,990), and the 2,709 distinct
allosteric mechanisms `test/Derivation/test_admission.jl` puts through the guard.

| Population | Admitted before | Newly admitted | Lost |
|:--|--:|--:|--:|
| Ping-pong taggings (50) | 0 | 40 | 0 |
| `:OnlyA` taggings of uni-uni, bi-bi and ping-pong inits (9,616) | 7,520 | 180 | 0 |
| Taggings with an `:OnlyI` group (43,311) | 2,200 | 34 | 0 |
| Ter-ter taggings (2,000) | 1,350 | 253 | 0 |
| K-type taggings of the bi-bi split children (9,878) | 7,474 | 600 | 0 |
| Bi-bi with a competitive inhibitor (3,000) | 2,503 | 0 | 0 |
| Ter-bi (3,000) | 2,126 | 160 | 0 |
| Bi-ter (3,000) | 2,145 | 158 | 0 |
| One expansion round (1,990) | 1,939 | 51 | 0 |
| The admission tests' mechanisms (2,709) | 2,260 | 99 | 0 |

- Every newly admitted mechanism of 8 steps or fewer agrees with both oracles under the
  harsh test: 591 of 591. Of the larger ones, 598 of 598 sampled agree under the mild
  test, and a spaced 100 of the ter-ter ones agree under the harsh test.
- With the filling check switched off, 18 of 250 mechanisms of the second population fail
  the harsh test. With the set-aside rule switched off the guard admits 125 fewer and no
  other mechanism changes sides.
- At two catalytic subunits, the concerted network's steady state equals the law for 40 of
  40 sampled newly admitted mechanisms.
- The ter-ter count is one below the spike's. The spike set aside a form no step feeds when
  any of its alive transitions led into a set-aside group; the rule as stated wants all of
  them to. The entry concerned, one of the two accepted bridges, converges under the mild
  test and stays rejected.
- The third built reading changes no verdict in any population once the search exists: the
  search finds a reading for each of the 4 mechanisms it admits, and it decides nothing for
  the 441 too large to search.
- The slowest guard call is 0.3 s (a ter-bi mechanism of 25 steps with 2,157 pairs), above
  the spike's 0.1 s; every population still takes under 2 s. 33 mechanisms are too large to
  search and have no consistent built reading.
- What stays rejected though valid, from 40 sampled per population of the mechanisms of 8
  steps or fewer refused by a face check: 9 of 40 in the second population (an estimated
  127 of 566), 7 of 40 ter-bi, 1 of 40 bi-ter, 4 of 8 ping-pong taggings, none in the
  `:OnlyI`, ter-ter, split and inhibitor populations. Nearly all are ping-pongs that the
  filling check refuses under the reading the oracles accept, at the second substrate's
  binding to the covalent form; one of the 9 is refused as a bridge instead.

The spike scripts and logs are kept at
`/home/denis.linux/claude-scratchpad/lean-mechanism-spike/`, outside the repository. They
are a reference for the comparison, not code to keep.

### 3. Compiling

- `Mechanism` has no type parameter. The macro, Enumeration and `identify_rate_equation`
  return plain mechanisms.
- `Derivation` owns code generation. Its one function for it is exported by the package
  (the maintainer's decision):

  ```julia
  eq = compile_rate_equation(m)   # compiles once
  v = eq(concs, params)           # fast, no allocation
  ```

  `eq` is a callable object whose type carries the graph (`Sig`), so its call is generated
  per mechanism as `rate_equation` is today. `Sig` holds the graph as nested tuples of
  Symbols, Ints and Bools; the builder falls back to an integer encoding only if Julia
  rejects that as a type parameter. No dependency is added.
- `rate_equation(m, concs, params)` is removed altogether (the maintainer's decision): one
  way to evaluate a law, and no slow path. Users who want the law in other code copy
  `rate_equation_string`.
- `FittingProblem(m, data)` compiles inside and keeps both the mechanism and `eq`. Fitting
  reads the parameter names as a compile-time constant from `eq` (`parameters(eq, kind)`);
  `parameters(m, kind)` on a mechanism is an ordinary function.
- Plotting compiles once per plot. Tests evaluate a law through `compile_rate_equation`.
- The speed, allocation, compile-budget and type-stability tests move to `eq`.

Interface compared, as `.claude/CLAUDE.md` asks:

| Sketch | Caller writes | Trade-off |
|:--|:--|:--|
| A function returning a callable (chosen) | `eq = compile_rate_equation(m)`, then `eq(concs, params)` | One name, one way to evaluate. |
| A compiled object passed to `rate_equation` | `rate_equation(eq, concs, params)` | Keeps a second name, and invites a slow `rate_equation(m, …)` beside it. |
| Fitting builds the type itself | `Derivation.CompiledRateEquation{Derivation.sig(m)}()` | Exposes the encoding to Fitting. |

### 4. One macro

`@allosteric_mechanism` is removed. `@enzyme_mechanism` reads both kinds:

```julia
@enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    modifiers: I
    catalytic_multiplicity: 2              # only with allosteric states
    steps: begin
        E + A ⇌ EA                         :: EqualAI
        EA + B <--> EPQ                    :: OnlyA
        EQ + P ⇌ EPQ                       :: EqualAI
        E + Q ⇌ EQ                         :: EqualAI
        (E + I ⇌ EI, EA + I ⇌ EAI)         :: OnlyI
    end
    regulatory_site(2): X::OnlyA, Y::OnlyI # only with allosteric states
end
```

- A name on a step side that is not a declared metabolite is a form name, used as written.
  Each side holds exactly one form.
- **The free enzyme is the form named `E`.** This is the only place a form's name carries
  meaning. The macro's docstring states it, and a mechanism with no form named `E` gets an
  error that states it.
- A mechanism is allosteric when its steps carry allosteric states: every statement has
  one or none has.
- `X::Mod` and the call notation `E(A, B; residual = A - P)` are gone. A step that binds a
  reactant as a modifier simply takes it up: `ES + S ⇌ ESS`.
- The macro checks what it reads (names, one form per side, `E` present, states), builds
  the mechanism, then runs the guard, whose message names steps by their text.
- A misspelt form name makes a new form. The connectivity check and the guard catch some
  such errors, not all.
- `show` prints a mechanism as this macro writes it, and `parse_mechanism` reads that text
  back. Mechanism text saved by earlier versions does not parse.

### 5. Who knows what

- **`Derivation`** computes from the graph, in internal functions, what it reads today
  from stored facts: the kinetic groups, what a step is by its shape (no metabolite, one
  taken up, or an exchange), and the constant names from the form and metabolite names
  (`K_EAB_to_EA_B`, `k_ES_to_EP`), built at the last step of derivation. It rejects a
  mechanism whose constants would share a name.
- **The dependent constants.** `_priority_key` reads step kinds and form contents to pick
  which constants are dependent. It will read the step's shape and canonical order
  instead. The law is unchanged, but some fitted-parameter names change; the pull request
  lists which pins moved.
- **A conformation's steps.** `Derivation` and `Enumeration` both read the steps that run
  in a conformation (today the stored `active_steps` and `inactive_steps`). Two public
  functions in `Mechanisms` compute them from the allosteric states and the graph:
  `steps_in_active_conformation(m)` and `steps_in_inactive_conformation(m)`, each a
  `BitVector`. The alternatives were storing the masks (today's stored facts) or letting
  each module derive them (two owners of "what a conformation keeps").
- **What Enumeration does not build** (the maintainer's decision). The guard admits a
  mechanism in which a form stays although every binding into it is excluded, such as the
  third example of decision 2: it is correct, so a user may write it. Enumeration does not
  propose it. It keeps a child only if, in each conformation, every form that stays and
  has bindings into it has one that still runs.
- **Provisional: how Enumeration gets form contents.** Until the Enumeration design says
  otherwise, Enumeration owns the naming rule for the forms it builds (today's
  `_form_name`: `E`, `EAB`, `EPmodQ`, `EB_res_plusA_minusP`) and the check that names
  decode one way (`check_decodable_names`, run when enumeration starts, not in
  `@enzyme_reaction`), and reads contents back from form names where its public functions
  take a `Mechanism`. A mechanism whose form names are not the enumerator's then cannot be
  expanded. If that design finds that the moves can read what they need from the graph,
  the naming rule, the check and the limit all go.
- **`Mechanisms`** keeps only what two modules use. `segment_walk` and `present_species`,
  which only `Derivation` calls, move into `Derivation`.
- **Deferred to the Derivation-interface design.** Enumeration calls
  `steps_in_active_conformation`, `steps_in_inactive_conformation`,
  `rapid_equilibrium_segments` and `has_turnover_cycle` only to predict properties of the
  rate law from the graph (a maximal rate, first order in each metabolite, an invisible
  constant, which flips add a parameter, inactive turnover). Whether it should ask
  `Derivation` instead, which would make all four internal to `Derivation`, is settled
  there, with timings.

Public names of `Mechanisms` after the change, from 25: `Mechanism`, `CatalyticScheme`,
`Allostery`, `Step`, `RegulatorySite`, `EnzymeReaction`, `Reactant`, `Regulator`,
`metabolites`, `step_text`, `steps_in_active_conformation`,
`steps_in_inactive_conformation`, `rapid_equilibrium_segments`, `has_turnover_cycle`,
`has_unbalanced_cycle`: 15. The builder removes `has_unbalanced_cycle` if each module can
own its use; the Derivation-interface design may remove four more.

Exports after the change: `Mechanism` and `compile_rate_equation` join;
`EnzymeMechanism`, `AllostericEnzymeMechanism`, `@allosteric_mechanism` and
`rate_equation` leave.

## Tests

- **Pin first.** Before the data model moves, pin what must not change: the law of every
  `MECHANISM_TEST_SPECS` entry against its oracle, and the enumerator's populations.
- **Fixtures** translate mechanically from call notation to form names, by a script that
  applies today's `_form_name` (`E(B; residual = A - P)` becomes `EB_res_plusA_minusP`).
  The enumeration rule in `.claude/CLAUDE.md` ("write chemistry the enumerator's way") is
  rewritten in form names.
- **The oracles stay independent.** The test helpers that model an excluded step
  (`_testhelper_network`, `_testhelper_exclusion`, `_testhelper_epsilon_span`,
  `_testhelper_corner_directions`, `_testhelper_lp_existence`) get their own statement of
  the table in decision 2; they do not call the package's.
- **Verdicts that change.** The guard's promise does not change: it admits a mechanism
  only when its law is the unique limit of mass-action kinetics. A test may change only in
  these ways, each listed in the pull request:
  - a fixture rejected today is admitted, and the test then asserts, through the
    independent oracles, that its law equals the numerical limit on every face;
  - a fixture admitted today is rejected because an excluded rapid-equilibrium step would
    have to stop, and the test then asserts that fault;
  - a fixture rejected today is rejected for another fault, or with another message;
  - the constructor's `:OnlyI` error becomes the guard's verdict.

  The checks that stop a law which is not the limit keep their force: a limit that depends
  on how fast an isomerization vanishes, a pumped basin, a bridge, a second closed class,
  a trap, a split, the conservative rejections. Several of their hand-built fixtures hold
  an excluded rapid-equilibrium step, so the guard now rejects those fixtures earlier, for
  that step. Each such check therefore also gets a fixture that passes the earlier checks
  (the same mechanism with that step steady-state), on which the guard returns the
  face-level fault and both oracles agree on every face. A check that no mechanism can
  reach any more is reported to the maintainer, not left untested in silence.
  (An earlier draft said these tests pass unchanged; that was measured before the
  principle that a step keeps its kind was added.)
- **The guard rule is a gate.** Stage 1 lands and passes before any fixture is
  translated. The work stops and returns to the maintainer if a mechanism admitted today
  is rejected for any reason but an excluded rapid-equilibrium step that today's guard
  read as stopped, or if a newly admitted law is not the numerical limit. Such a step
  shows as the new fault where it cannot lose its affinity, and as a bridge in the two
  ter-ter mechanisms of decision 2's evidence, where it can. The data model leaves no
  room for a per-step flag.
- **The gate for the energy-level rule.** It lands only when all of these hold, each
  reported in the pull request:
  - **Nothing lost.** No mechanism the built three-try guard admits is rejected, in every
    population below.
  - **Every small newly admitted mechanism passes the harsh test.** Each mechanism of 8
    steps or fewer that the rule admits and the built guard does not agrees with both
    oracles on every face, the ε-path oracle with corner directions.
  - **Larger ones.** A spaced sample of at least 150 per population under the mild test,
    and one long run of the harsh test on at least 100 newly admitted ter-ter mechanisms.
  - **Wider populations than the spikes'.** The five of decision 2's evidence, and: the
    `:OnlyA` taggings of a bi-bi reaction with a competitive inhibitor; of a ter-bi and a
    bi-ter reaction; one expansion round of the census reactions (the enumerator's own
    output); and the distinct allosteric mechanisms `test/Derivation/test_admission.jl`
    puts through the guard.
  - **The failure mode, exhaustively.** Every mechanism whose admitting reading (the one
    the guard reports) has a vanishing transition alive on a face while its reverse is
    dead is listed, and every such mechanism of 8 steps or fewer passes the harsh test
    under that reading. Every oracle check in the gate uses the admitting reading.
  - **The checks one at a time.** With the filling check off, the harsh test must reject at
    least 14 mechanisms of the 9,616, entries 6580, 6595 and 8524 of that population among
    them (it proves the gate can see the failure). With the set-aside rule off, the
    verdicts of everything else are unchanged.
  - **More than one subunit.** The verdict does not depend on `catalytic_multiplicity`
    (already pinned), and a sample of newly admitted mechanisms agrees with the suite's
    protomer-network oracle at two subunits where that oracle applies.
  The work stops and returns to the maintainer if a newly admitted mechanism fails an
  oracle, or a mechanism the built guard admits is lost.
- **Permanent tests for the rule.** The suite keeps: the six-step ping-pong with both
  chemistry steps `:OnlyA`, admitted, its law against the limit; a mechanism the filling
  check refuses and its near-twin that it passes, each against the harsh oracle; a
  mechanism the set-aside rule admits and one it must not; a sampled population test with
  corner directions on small mechanisms; and the tests' own statement of the rule, which
  does not call the package's.
- **What counts as checked.** The comparisons of decision 2 are rerun once on the
  implemented rule and reported in the pull request, with a sample of at least 2,000
  ter-ter mechanisms and the 9,878 K-type taggings of the bi-bi split children added. The
  suite keeps a small permanent version: a sampled
  population test that the only mechanisms lost are the ones the rule names, and the
  oracle agreement on a sample of the newly admitted ones.

## Documentation

The assumption behind an excluded step is explained to users, not only to developers (the
maintainer's request). It lands with stage 1, when the rule changes:

- **A user-facing section, "How an excluded step is switched off",** on the docs page that
  introduces allosteric mechanisms. It gives, in the user's terms: what `:OnlyA` and
  `:OnlyI` assume (a very high dissociation constant or very low rate constants, with
  every thermodynamic constraint kept in each conformation); that a step keeps its kind,
  rapid equilibrium or steady state, in both conformations; the energy picture and the
  four-row table of decision 2; that the package looks for a choice of energies under
  which the rate equation is the limit, and accepts the mechanism when one exists; the
  four worked examples, the ping-pong among them; what the error means when a
  rapid-equilibrium step lies inside one group, with the two fixes (make the step
  steady-state, or exclude the steps that keep its two sides joined); and that the last
  check is validated numerically, with what was tested.
- **The mechanism macro's docstring** states the rule in two or three sentences and links
  to that section.
- **The guard's error message** for that fault names the step and says the same two fixes.
- **`docs/src/developer.md`** carries the full rule and how the guard checks it.

## Left to the builder

- The canonical order of forms and steps, and the orientation of a step that exchanges
  several metabolites.
- The replacement rule in `_priority_key`.
- The name and message of the fault for an excluded rapid-equilibrium step that would
  have to stop, what happens to `_NO_LIMIT`, and the rewritten guard section of
  `docs/src/developer.md`.

## Sequence

1. This design is approved. Its first stage is planned and built alone: the guard's rule
   on today's data model. It removes `Derivation`'s last read of form contents, touches
   no Enumeration internals, and is the gate. Stage 1 was built first with three fixed
   readings (reviewed, gate clean); the energy-level rule is then built on top of it, in
   the same branch, before the enumeration tests are brought to the rule.
2. The Enumeration design: what its moves need to know about a form.
3. The Derivation-interface design.
4. The rest of this design is planned and built: the lean graph, the one type and the one
   macro with the fixtures, then compiling inside `FittingProblem`.

The `E_total` design can come anywhere in that list.

## Consequences for the maintainer to know

- **Search results change.** The enumerator stops proposing mechanisms in which an
  excluded rapid-equilibrium step would have to stop (the numbers in decision 2), so
  `identify_rate_equation` can return a different best mechanism on data where one of
  those won.
- **Allosteric ping-pong.** `ROADMAP.md` plans "let an allosteric ping-pong stop both
  half-reactions in the inactive conformation" and says admission rejects it. Under the
  energy-level rule most are admitted (40 of 50 measured), with no second flip. The
  tagging move already prefers the child that stops all chemistry when the guard admits
  it, so the enumerator's ping-pong children change without new enumerator code. The
  entry is rewritten as part of the work.
- **The guard's last check rests on numerical evidence.** The filling check and the
  set-aside rule passed every test run, but no proof shows the face checks sufficient for
  every reading. A wrong law the guard admitted would not be caught by a fit: it is a
  smooth function that can fit well, with a wrong mechanism label and unreliable
  extrapolation toward the face where it fails. The gate in Tests is the protection.
- **A check that nothing reaches.** With the built guard the tagging move's hyperbolic
  drop above one subunit no longer drops anything. It is re-examined once the energy-level
  rule is in, and kept or removed on the maintainer's decision.
- **Valid mechanisms still rejected.** Those above the search cap that only an untried
  reading would admit; a few the face condition refuses as a bridge; and those the filling
  check refuses although their rate equals the law's. `ROADMAP.md` records them.
- **Fitted-parameter names** change for some mechanisms (decision 5).
- **Saved mechanism text** from earlier versions does not parse (decision 4).
- **The decodability check** moves from `@enzyme_reaction` to the start of enumeration
  (decision 5, provisional).
- **Version.** The change is breaking: 0.11.0, bumped once.
- **`ROADMAP.md`.** The stored-facts known issue is measured again and retired or updated;
  "Shrink the module interfaces" is updated for `Mechanisms`.
