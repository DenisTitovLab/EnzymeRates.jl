# One lean `Mechanism`: design

Date: 2026-11-18. Approved by the maintainer on 2026-11-19; decisions 6 and 7 were added
with him on 2026-11-20, and he answered the builder check's questions the same day.
Stage 1, the guard's rule (decision 2), is built in 0.11.0. This document now specifies
what is left: decisions 1 and 3 to 7. All of it lands in the same pull request as stage 1.

## Why

Too many names cross between `Mechanisms`, `Enumeration` and `Derivation`. `Mechanisms`
has 25 public names, four mechanism types, an abstract type, conversions both ways and a
set of stored facts that only some callers read. A reader who wants to change the
derivation must learn enzyme forms, residuals, conformation labels, step kinds and the
integer encoding of a compiled type, although the rate equation depends on none of them.

The rate equation depends only on a reaction graph: which forms exist, which steps join
them, what each step takes up and gives off, which steps are fast and which share
constants. This design makes that graph the mechanism. An allosteric mechanism adds its
allosteric states, its oligomeric state and its regulatory sites.

## Words used here

- **Conformation:** the active or the inactive shape of the enzyme (MWC). A
  **conformation label** is something else: the prefix of a form in the enumerator's
  record (`E`, `Estar`, `F`), which tells apart forms of the enzyme that an isomerization
  joins.
- **Kinetic group:** steps that share their constants.
- **Allosteric state** of a kinetic group, or of a regulatory site's ligand: one of
  `:OnlyA`, `:OnlyI`, `:EqualAI`, `:NonequalAI`, written after `::` in the macro. The
  `oligomeric_state` is unrelated: the number of subunits.
- **Excluded step:** a step that does not run in a conformation because of its group's
  allosteric state: an `:OnlyA` step in the inactive conformation, an `:OnlyI` step in
  the active one.
- **Cut-off form:** a form that the steps still running in a conformation do not join to
  free E.
- **Island:** in a conformation, a set of forms joined to each other by the steps the
  conformation does not exclude (decision 2).

## Goals

1. One mechanism type, `Mechanism`, a plain value, for plain and allosteric mechanisms.
2. The mechanism stores a lean reaction graph and nothing derived from it.
3. `Derivation` reads only that graph: no form contents and no step kinds.
4. Form names are opaque Symbols.
5. The guard still admits a mechanism only when its law is the unique limit of
   mass-action kinetics, and its rule for excluded steps is read from the graph alone.

## Non-goals

Outside this design: new enumeration moves, iso mechanisms, a flip at a covalent form,
and identify's keyword arguments other than `Keq` (decision 7).

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

1. **Lean graph only (chosen).** The guard read one thing from form contents: whether an
   excluded step is a plain binding. The rule of decision 2, now built, reads the graph.
   The constant names and `_priority_key` still read contents; decision 5 replaces those
   reads.
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
4. **As 3, and a step keeps its kind in both conformations (built first).** The maintainer's
   principle: rapid equilibrium and steady state differ more than 100-fold, and a step
   does not jump from the fastest kind in one conformation to slower than every
   steady-state step in the other. A steady-state uptake into a cut-off form keeps its
   kind either way, so the guard may read it as stopped too, and tries three fixed
   choices (added after the gate and the enumeration tests, on the maintainer's
   decisions). Built, and kept as the first three readings of approach 5.
5. **Energy levels (chosen).** Every cut-off island of forms lies higher in energy; an
   excluded step between two islands loses its uphill rate in either direction, or stops
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
   is built. One concrete type everywhere inside the package. A user who wants a fast
   rate function outside a fit calls `compile_rate_equation` (decision 3).
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
    oligomeric_state::Int             # subunits that change conformation together
    regulatory_sites::Vector{RegulatorySite}
end

struct Mechanism
    catalytic::CatalyticScheme
    allosteric::Allostery
end
```

- A mechanism is plain when `allosteric.allosteric_states` is empty; its
  `oligomeric_state` is then 1.
- `oligomeric_state` replaces `catalytic_multiplicity` (the maintainer, 2026-11-20), in
  the field and in the mechanism macro. `@enzyme_reaction` already writes
  `oligomeric_state:`; the reaction's field `allowed_catalytic_multiplicities` becomes
  `allowed_oligomeric_states`. Several catalytic sites on one subunit would belong to the
  catalytic scheme. A regulatory site keeps its `multiplicity`.
- `RegulatorySite` names its ligands' states `allosteric_states` too, in place of
  `mwc_tags`: both hold the same four values.
- `CatalyticScheme.metabolites` holds the metabolites the steps take up or give off. One
  that is neither a substrate nor a product is a modifier. A regulatory site's ligands
  are named in the site, and are listed there only if a step also binds them.
  `metabolites(m)` lists both, the sites' ligands last, as today.
- A reactant that binds as a modifier (substrate inhibition) is a step that takes up that
  reactant; the graph has no separate modifier copy.
- The free enzyme is the form named `E`, stored first. The constructor requires it, so
  that a mechanism built directly prints as text the macro reads back (decision 4).
  Nothing else in the package reads a form's name.
- Two writings of one mechanism are `==` and hash alike. Step order, group order, the
  order of the forms and the written direction of a step do not matter, with one
  exception: the direction of an isomerization relative to the others of its kinetic
  group (below). Form names are part of a mechanism's identity.
- The constructor checks representation only: indices in range, a form named `E`, every
  form joined to it, distinct names, uniform kinetic groups, one known allosteric state
  per group. It does not run the guard, which alone decides which groups may be
  `:OnlyI`.

#### Which way a step is stored

A form is only a name, so the constructor orients each step from the graph (the
maintainer's rule, 2026-11-20). The stored direction names the step's constants and
decides whether a rapid-equilibrium isomerization's constant is `K` or `1/K`.

The first row that decides gives the direction:

| Step | Stored direction (`from` to `to`) |
|:--|:--|
| Only takes up, or only gives off, metabolites | As the uptake. |
| Takes up some and gives off others | The way the reaction runs. Count the substrates the step takes up and the products it gives off, less the substrates it gives off and the products it takes up. The direction in which the count is positive is stored. |
| The count is zero (one modifier for another), or the step is an isomerization | The way the reaction runs: along the shortest cycle through the step that turns the reaction over. |
| No such cycle, or one each way of equal length | Away from the free enzyme: `from` is the end fewer steps from it. |
| Both ends equally far | `from` is the end whose form comes first in `enzyme_forms`. |

**The steps of one kinetic group turn together.** The steps of a group share their
constants. Steps that take up or give off metabolites all take up and give off the same
ones, and that fixes how the members correspond. For isomerizations it is a fact about
the mechanism that the graph does not hold. In `(Estar <--> E, EstarA <--> EA)` the
enzyme's isomerization has one pair of constants at both forms, yet the reaction runs
through the two steps in opposite directions. The table orients one member and the others
follow it: steps with metabolites by their lists, isomerizations by the directions they
were written in relative to it. That member is the one nearest the free enzyme among
those the direction of the reaction decides, or among all of them when it decides none
(the maintainer's choice, 2026-11-20). So turning a whole group as written changes
nothing, and turning one isomerization against its group writes a different mechanism,
one that shares other constants. Nothing can check that a group is written aligned: the
writer aligns it, as the enumerator and the translated fixtures do.

Checked on today's data model
(`claude-scratchpad/lean-mechanism-spike/orientation/`): the table gives the stored
direction of all 3,566 isomerizations in the seeds and first expansions of four reactions
and in 3,000 ter-ter seeds. Of the 1,301 isomerizations in the mechanisms the tests
write, 32 turn, in at most 25 fixtures. Each is alone in its kinetic group, and is an
isomerization of the enzyme itself (`E <--> F`) or a step off the cycle, which today
falls back to the order of the forms. All 46 that share a group keep today's direction;
in three groups the table alone would turn the members apart, which is why a group turns
together.

Of the steps that take up some metabolites and give off others, the second row gives the
stored direction of all 913 enumerated ones. Of the 121 the tests write, it gives today's
direction for 100 and turns one (`E(S) + P ⇌ E(P) + S`, stored today against the
reaction). The other 20 exchange one modifier for another and fall to the later rows,
where at most two turn.

Removed from the public names of `Mechanisms`: `AllostericMechanism`, `EnzymeMechanism`,
`AllostericEnzymeMechanism`, `AbstractEnzymeMechanism`, `EnzymeForm`, `Metabolite`,
`catalytic`, `compile_mechanism`, `check_conformation_label`, `form_text` (a form is its
name) and, by decisions 5 and 6, `excludes`, `present_species`, `segment_walk` and
`check_decodable_names`. Also removed: the stored facts (`modifiers`, `groups`,
`step_kinds`, `toward_products`, `step_sides`, the four step masks) and today's integer
encoding of a mechanism in a type parameter.

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

The picture is energy. The steps a conformation does not exclude join the forms into
**islands**: sets of forms joined to each other by such steps. An island is not a kinetic
group and not a rapid-equilibrium segment. Free E's island holds the law's forms. Every
other island is cut off: in that conformation it lies higher in energy, and so holds no
enzyme. Then:

| Excluded step | How it is switched off |
|:--|:--|
| Between two islands | It loses its uphill rate and keeps its downhill rate, toward the side the reading puts lower. A rapid-equilibrium step stays fast. |
| Between two islands, steady state | It may be stopped instead: both rates vanish. |
| Inside one island, steady state | It is stopped. |
| Inside one island, rapid equilibrium | Fault. Its two sides have one energy, so thermodynamics pins its equilibrium constant, and it could be switched off only by turning slow. |

Which island lies above which is not given. A **reading** chooses, for each excluded step
between two islands, which side is higher, or, for a steady-state step, that it stops. A
reading is consistent when the steps that lose, each pointing uphill:

1. close no cycle, so energies exist;
2. never point into free E's island, since nothing lies below it; and
3. lead from free E's island to every island, so every cut-off island lies above it and
   has a step by which it empties.

The law is the same under every reading, since the law drops the excluded steps. The guard
admits the mechanism when some consistent reading passes the guard's face checks
(`docs/src/deriving/mwc_allostery.md`).

Why a step loses its uphill rate: the law assumes the cut-off forms are empty, and a step
into a higher island that kept its uphill rate would fill it. "Loses its affinity" is the
case in which the step is a binding and its bound side is the higher one.

Why a rapid-equilibrium step may lose in either direction but never stop: the maintainer's
principle. Rapid equilibrium and steady state differ more than 100-fold, and turning a
rapid-equilibrium step into the slowest step of the other conformation is not a reasonable
approximation. A mechanism that needs it is written with that step steady-state.

Why thermodynamics holds: each island has one energy level, and every excluded step's
equilibrium constant shifts by the difference between its two ends' levels. The shifts
cancel around every cycle.

Four examples, each in the inactive conformation, the unlisted steps `:EqualAI`:

| Mechanism | Verdict |
|:--|:--|
| `E + S ⇌ ES :: OnlyA`, `ES <--> EP :: OnlyA`, `E + P ⇌ EP` | Admitted. `ES` is cut off: S binding loses its uptake, and the chemistry loses its rate into `ES` or stops. |
| `E + S ⇌ ES :: OnlyA`, `ES <--> EP`, `E + P ⇌ EP` | Fault. `ES` stays in free E's island, made from P, so S binding is inside one island and a rapid-equilibrium step cannot be stopped. |
| `E + S <--> ES :: OnlyA`, `ES <--> EP`, `E + P ⇌ EP` | Admitted. The steady-state binding is inside one island and is stopped. Enumeration does not build this one (decision 5). |
| Ping-pong `E + A ⇌ EA`, `EA <--> FP :: OnlyA`, `F + P ⇌ FP`, `F + B ⇌ FB`, `FB ⇌ EQ :: OnlyA`, `E + Q ⇌ EQ` (F the modified enzyme) | Admitted. `F`, `FP` and `FB` are one cut-off island: `EA → FP` and `EQ → FB` vanish, their reverses keep their rates, and the island empties through both halves. |

#### Built

This decision is built (0.11.0). The search for a reading, the face checks, the
comparison with the earlier guard and the evidence are in `docs/src/developer.md` (the
guard), `docs/src/deriving/mwc_allostery.md` ("How an excluded step is switched off") and
`ROADMAP.md`; the gate's numbers are in
`claude-scratchpad/lean-mechanism-spike/gate/report.md` on the maintainer's machine. What
is left to build must not change a verdict (see "Tests").

### 3. Compiling

- `Mechanism` has no type parameter. The macro and `identify_rate_equation` return plain
  mechanisms; Enumeration returns candidates that each hold one (decision 6).
- `Derivation` owns code generation. Its one function for it is exported by the package
  (the maintainer's decision):

  ```julia
  rate = compile_rate_equation(m)   # compiles once
  v = rate(concs, params, Keq)      # fast, no allocation
  ```

  `Keq` is its own argument and `E_total` is gone (decision 7).

  `rate` is a callable object whose type carries the whole mechanism (`Sig`: the
  catalytic scheme and the allostery), so its call is generated per mechanism as
  `rate_equation` is today. `Sig` holds it as nested tuples of Symbols, Ints and Bools;
  the builder falls back to an integer encoding only if Julia rejects that as a type
  parameter. No dependency is added.
- `rate_equation(m, concs, params)` is removed altogether (the maintainer's decision): one
  way to evaluate a law, and no slow path. Users who want the law in other code copy
  `rate_equation_string`.
- `FittingProblem(m, data, Keq)` compiles inside and keeps both the mechanism and `rate`.
  Fitting reads the parameter names and the metabolite names as compile-time constants
  from `rate` (`parameters(rate)`, `metabolites(rate)`), as `loss!` needs. `Derivation`
  gives the two functions a method for the compiled rate; they are methods, not new
  names. `parameters(m)` and `metabolites(m)` on a mechanism are ordinary functions.
- Plotting compiles once per plot. Tests evaluate a law through `compile_rate_equation`.
- The speed, allocation, compile-budget and type-stability tests move to `rate`.

Interface compared, as `.claude/CLAUDE.md` asks:

| Sketch | Caller writes | Trade-off |
|:--|:--|:--|
| A function returning a callable (chosen) | `rate = compile_rate_equation(m)`, then `rate(concs, params, Keq)` | One name, one way to evaluate. |
| A compiled object passed to `rate_equation` | `rate_equation(rate, concs, params)` | Keeps a second name, and invites a slow `rate_equation(m, …)` beside it. |
| Fitting builds the type itself | `Derivation.CompiledRateEquation{Derivation.sig(m)}()` | Exposes the encoding to Fitting. |

### 4. One macro

`@allosteric_mechanism` is removed. `@enzyme_mechanism` reads both kinds:

```julia
@enzyme_mechanism begin
    substrates: A, B
    products: P, Q
    modifiers: I
    oligomeric_state: 2                    # only with allosteric states
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
- **The free enzyme is the form named `E`.** This is the only meaning a form's name
  carries, and the constructor requires it (decision 1). The macro's docstring states it,
  and a mechanism with no form named `E` gets an error that states it.
- **`allosteric_regulators:` is gone from the mechanism macro.** A regulatory site names
  its ligands, so the line repeated them. `modifiers:` stays: a name on a step side must
  be told from a form name. `@enzyme_reaction` keeps its `allosteric_regulators:`.
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
  (`K_EAB_to_EA_B`, `k_ES_to_EP`), built at the last step of derivation.
  `check_valid_rate_equation` rejects a mechanism whose constants would share a name.
- **`mod` leaves the metabolite part of a constant's name.** Today a step that binds a
  reactant's modifier copy writes the metabolite as `Smod` in its constants' names. The
  lean step takes up the reactant, so that part is the bare `S`. The form's name still
  tells the two bindings apart (`ESmod` beside `ES`), so no two constants share a name.
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
- **How Enumeration gets form contents.** It keeps its own record of them beside each
  mechanism it builds (decision 6). Nothing in `src/` reads a form's contents from its
  name.
- **`Mechanisms`** keeps only what two modules use. `present_species` becomes private:
  it is five lines, and `Derivation` writes its own. `segment_walk` becomes private too
  (the maintainer, 2026-11-20): `rapid_equilibrium_segments` keeps it, and the guard in
  `Derivation`, its one outside caller, walks its own network of both conformations
  itself.
- **`Mechanisms.excludes(state, conformation)` is removed** (the maintainer,
  2026-11-20). It restates what the states' names say: `:OnlyI` is absent from the active
  conformation and `:OnlyA` from the inactive. Each of its five callers compares the
  state itself (`state === :OnlyI`). `Derivation` and `Enumeration` both need the steps a
  conformation excludes by their state, beside the steps that still run, and each builds
  that mask with the same comparison. The states stay Symbols.
- **Enumeration predicts the rate law from the graph itself** (the maintainer's decision,
  2026-11-20; decision 7). It calls `steps_in_active_conformation`,
  `steps_in_inactive_conformation`, `rapid_equilibrium_segments` and `has_turnover_cycle`
  for that, so the four stay public in `Mechanisms`.

Public names of `Mechanisms` after the change, from 25: `Mechanism`, `CatalyticScheme`,
`Allostery`, `Step`, `RegulatorySite`, `EnzymeReaction`, `Reactant`, `Regulator`,
`metabolites`, `step_text`, `steps_in_active_conformation`,
`steps_in_inactive_conformation`, `rapid_equilibrium_segments`, `has_turnover_cycle`,
`has_unbalanced_cycle`: 15. The last stays because Enumeration's seeds and Derivation's
admission both call it. Decision 1 lists the 14 names that go.

Exports after the change: `Mechanism` and `compile_rate_equation` join.
`EnzymeMechanism`, `AllostericEnzymeMechanism`, `@allosteric_mechanism`, `rate_equation`,
`IndependentParameters`, `FittedParameters` and `AllParameters` are removed, not only
unexported.

### 6. Enumeration's interface: candidates

Decided with the maintainer on 2026-11-20. His principle: the enumerator relies only on
the arguments of `identify_rate_equation`, and everything else about it is internal.
`identify_rate_equation` is its only caller in `src/`, and it expands only mechanisms the
enumerator built.

#### What the moves read from a form

| Where | What it reads |
|:--|:--|
| Seeds | Builds every form from its bound metabolites and covalent residual |
| Dead-end move | Which metabolites a form holds, to decide where an inhibitor binds |
| Split move | Divides a kinetic group by what else the form holds, and by its residual |
| Redundant-inhibitor filter | Pairs an inhibitor complex with a productive form of the same composition |
| Atom-conservation check | Bound metabolites plus residual |

The moves also read whether a metabolite binds as a reactant's modifier copy. The lean
graph holds none of this.

#### Approaches considered

1. **Candidates carry their contents (chosen).** The enumerator returns its own value: the
   lean `Mechanism` plus its record of each form. The search proposes exactly what it
   proposes today. Costs one public type and a test helper.
2. **Form names carry the contents.** The enumerator returns plain mechanisms and parses
   names such as `EB_res_plusA_minusP`. Rejected: names then carry meaning inside
   `Enumeration`, a parser and a decodability check stay in `src/`, and a misnamed form
   silently changes what a move does.
3. **A table from mechanism to contents inside the enumerator.** Rejected: hidden state,
   and identify expands parents on worker processes.
4. **Recompute the contents from the graph.** Rejected: the graph does not always give
   them. A merged complex is reached by binding A then B and by binding Q then P, so the
   graph says it holds both pairs, while the enumerator places inhibitors by the product
   pair. The populations would change, for a generality nobody needs.

#### The interface

```julia
struct CandidateMechanism
    mechanism::Mechanism       # the only field a caller reads
    # the enumerator's record of each form: private
end

candidate_mechanisms(reaction; max_param_count, eq_complexity_filter,
                     optional_allosteric_regulators,
                     optional_competitive_inhibitors) -> Vector{CandidateMechanism}
candidate_mechanisms(reaction, parent::CandidateMechanism; max_param_count,
                     eq_complexity_filter) -> Vector{CandidateMechanism}
```

Every keyword is a keyword of `identify_rate_equation`, passed through unchanged. Each
has a default that limits nothing: no size limit, and no optional regulator. Identify
owns the defaults a search uses (20 and 337) and always passes all four.

- The first method returns the mechanisms the search fits first (today's
  `starting_mechanisms`). It works out which regulators are required from the reaction
  and the two optional lists, which identify does today (`_required_regulators`).
- The second returns every child of `parent` by the expansion moves (today's
  `expand_mechanisms` on a one-element vector).
- **Every candidate returned is one identify may fit** (the maintainer's decision): it
  has passed the guard, as today, and it is within both size limits, which identify
  applies today. A mechanism is dropped when its law has more than `eq_complexity_filter`
  terms (`Derivation.rate_equation_term_count`) or more than `max_param_count` fitted
  constants (`Derivation.fitted_parameter_count`). The limits apply to the mechanisms
  returned, not inside a move's own search, so the search proposes what it proposes
  today: a mechanism over a limit was never fitted or expanded.
- Two candidates are `==` and hash alike when their mechanisms are. Form names are built
  from the contents, so equal mechanisms hold equal records.
- Errors: the reaction has no catalytic cycle; its metabolite names could give two forms
  one name (the first method checks this); no mechanism binds every required regulator;
  an internal check fails (atom conservation, the emission rules). Identify catches the
  error of one parent's expansion and records it, as today.
- The function and the type replace `starting_mechanisms` and `expand_mechanisms` as the
  public names of `Enumeration`, and those two are gone. `init_mechanisms`,
  `seed_mechanisms` and the moves stay internal, and the tests call them qualified as
  they do today. A test that called `expand_mechanisms` or `starting_mechanisms` calls
  the method that replaces it with the default keywords, so nothing is dropped for size.

Interface compared, as `.claude/CLAUDE.md` asks:

| Sketch | Caller writes | Trade-off |
|:--|:--|:--|
| One function, two methods (chosen) | `candidate_mechanisms(rxn; optional…)`, then `candidate_mechanisms(rxn, parent)` | One function. Identify passes the user's keyword arguments through and knows nothing of seeds, tiers or moves. |
| Two functions on candidates | `starting_mechanisms(rxn, required_allo, required_comp)`, `expand_mechanisms(parents, rxn)` | Today's shape. Identify computes the required sets and wraps each parent in a vector. |
| A candidate that holds its reaction | `candidate_mechanisms(rxn; optional…)`, then `children(parent)` | A shorter second call, but two names, and every candidate sent to a worker carries the reaction. |

#### Inside Enumeration

- **The record of a form:** its conformation label, the metabolites it binds
  productively, those it binds as modifier copies, and its covalent residual (substrates
  added, products subtracted). The type that holds it is private. The enumerator builds
  the label `E` only; the label stays (the maintainer, 2026-11-20) because iso mechanisms
  will need it, and the moves keep the behaviour their tests pin for other labels, such
  as the split of a group by conformation. `ROADMAP.md` says that only tests reach that
  code until iso mechanisms are enumerated.
- **Names.** Enumeration owns the naming rule for the forms it builds (today's
  `_form_name`: `E`, `EAB`, `EPmodQ`, `EB_res_plusA_minusP`). A name is for display and
  identity. Enumeration also owns the check that a reaction's metabolite names cannot
  give two forms one name (today's `check_decodable_names`), run when the first method is
  called. The check leaves the `EnzymeReaction` constructor and the mechanism macro,
  which run it today: a written mechanism needs only distinct form names.
- **Modifier copies.** In the lean graph a step that binds a reactant's modifier copy
  takes up the reactant; that it is a copy is in the record of the form it makes.
- **Step kinds and kinetic groups** (a plain binding, a binding that also does the
  chemistry, an isomerization, an exchange; the steps of each group) are computed inside
  Enumeration from the graph and its record.
- Its predictions of the rate law from the graph stay in Enumeration (decision 7).

#### Identify

Identify holds candidates. It fits, deduplicates and reports `candidate.mechanism`, hands
a candidate back to get its children, and returns plain mechanisms to the user.
`_required_regulators` moves into Enumeration. Identify's own private record of a
candidate (`_Candidate`, whose `mechanism` field is text) gets a name and field names that
cannot be confused with `CandidateMechanism`.

The enumerator keeps proposing mechanisms that differ only in the allosteric state of a
kinetic group whose steps a conformation has cut off anyway (354 of the 2,125 allosteric
candidates in the bi-bi seeds and their first expansion). They share one law, identify
fits it once, and their children may differ (the maintainer, 2026-11-20). `==` and `hash`
tell them apart.

#### Consequences

- **The starting call throws on an enumerator bug** (confirmed by the maintainer,
  2026-11-20). Today `starting_mechanisms` returns `(mechanisms, failures)`, and a
  degenerate seed whose flip expansion throws becomes a failure row while the search goes
  on without it. With one return shape for both methods, such an error stops the run at
  its start, seconds in. A parent's expansion error mid-search is still recorded and the
  run goes on.
- **`docs/src/identify/enumeration_engine.md`** keeps describing what the search
  proposes, organized by the starting set and the moves. It stops presenting
  `init_mechanisms`, `seed_mechanisms` and `expand_mechanisms` as functions to call, and
  shows `candidate_mechanisms` as the way to list what the search would propose.

### 7. Derivation's interface, `Keq` and `E_total`

Decided with the maintainer on 2026-11-20.

#### Who asks `Derivation` what, today

| Caller | What it asks | How |
|:--|:--|:--|
| Fitting | The law's value, the fitted names, a starting scale for each parameter, which constants are rate constants | 8 names. It also reads step fields and `step_kinds` to work out each parameter's units. |
| Identify | Per candidate: term count, validity, fitted names, the law's text, the dedup key | 5 functions and a kind type |
| Enumeration | Validity and the fitted-parameter count | 2 functions |
| Plotting, the macro | The law's value and names; validity | 3 names and a kind type |

#### Enumeration keeps predicting the rate law from the graph

Its filters predict six properties of a law without deriving it: a maximal rate, chemistry
in equilibrium with both sides, first order in each reactant, an invisible inhibitor
constant, which steps carry flux, inactive turnover. They stay in `Enumeration` (the
maintainer's decision): only the enumerator needs them, they are its search policy, and
asking `Derivation` would mean deriving every candidate's law before screening it. The
graph functions both modules use stay public in `Mechanisms` (decision 5).

#### The parameter interface

```julia
parameters(m)                                # the constants a fit estimates
parameter_units(m)                           # their units, for Fitting
dependent_parameter_values(m, params, Keq)
rate = compile_rate_equation(m)
rate(concs, params, Keq)                     # the rate per unit enzyme
```

- **`Keq` is its own argument** (the maintainer's design). It is a property of the
  reaction that is given, never fitted, and the signature says so. `parameters(m)` is then
  exactly the names of `params`. It depends on the conditions the data were measured in,
  so the user supplies it with the data.
- **`E_total` is gone.** The law returns the rate per unit enzyme; a user who wants an
  absolute rate multiplies by the enzyme concentration. Inside the package it was always
  1. It leaves the parameters, the printed equation and the docs.
- **One parameter list.** `IndependentParameters`, `FittedParameters` and `AllParameters`
  are removed. The dependent constants' names are the keys of
  `dependent_parameter_values`, and the printed equation shows their relations.
- **A fit's `params`** holds exactly `parameters(m)`; today it also carries `Keq`.
- **`rate_equation_string`** shows `Keq` as an input of its own.
- **`Keq` is always the last positional argument** (the maintainer, 2026-11-20), and
  nowhere a keyword: `rate(concs, params, Keq)`,
  `dependent_parameter_values(m, params, Keq)`,
  `test_rate_equation(m, params, data, Keq; absolute)`,
  `plot_rate_equation_on_data(m, params, data, Keq; ...)`,
  `FittingProblem(m, data, Keq; absolute, reference)`,
  `IdentifyRateEquationProblem(reaction, data, Keq; ...)` and
  `Data.RateData(table, metabolite_names, Keq; ...)`. The test and the plot read it from
  `params` today, and the other three take it as a keyword. The plot's method for a
  user's own rate function takes no `Keq`:
  `plot_rate_equation_on_data(f, params, data; ...)`.

`parameter_units(m)` returns each fitted constant's units as integer powers:

```julia
units = parameter_units(m)
units.rate            # Vector{Int}, one entry per name of parameters(m)
units.concentrations  # Matrix{Int}: rows follow parameters(m), columns metabolites(m)
```

| Constant | Power of the rate | Powers of concentration |
|:--|--:|:--|
| Rate constant of a steady-state step | 1 | −1 for each metabolite the step takes up in that direction |
| Dissociation constant of a rapid-equilibrium binding | 0 | +1 for the metabolite that binds |
| Other rapid-equilibrium constant | 0 | +1 for each metabolite given off, −1 for each taken up |
| Regulator's constant | 0 | +1 for its ligand |
| `L` | 0 | none |

The rate is per unit enzyme, in the units of the data's `Rate` column. Fitting starts each
constant at `r^rate × Π s_x^power` (`r` the data's typical rate, `s_x` each metabolite's
typical concentration) and, for relative data, rescales the constants whose power of the
rate is 1. `parameter_units` is public in `Derivation` and not exported. It replaces
`parameter_info`, `name`, `StepParameter`, `RegulatorParameter` and `AllostericParameter`,
and Fitting's reads of the steps.

Interface compared, as `.claude/CLAUDE.md` asks:

| Sketch | Caller writes | Trade-off |
|:--|:--|:--|
| `Keq` as an argument, one list (chosen) | `rate(concs, params, Keq)`; `parameters(m)` | `parameters(m)` is exactly what `params` holds, and the signature shows what is not fitted. |
| `Keq` inside `params`, one list | `rate(concs, params)`; `parameters(m)` | One argument fewer, but `parameters(m)` no longer lists everything `params` must hold. |
| Two kinds | `parameters(m)` for the law's inputs, `parameters(m, FittedParameters())` for a fit's | Keeps an export and a second meaning of `parameters`. |
| Today | three kinds, `E_total` merged in by every caller | |

#### What identify asks

Per candidate, identify asks `Derivation` three things: the fitted names (`parameters`),
the law's text (`rate_equation_string`) and the dedup key (`rate_equation_dedup_key`).
Today it asks five. The guard and the two size limits are the enumerator's (decision 6),
so identify's calls to `check_valid_rate_equation` and `rate_equation_term_count` go.

- **The progress log loses two counts** (accepted by the maintainer): how many children
  of a batch were skipped for each limit. A batch whose children are all too large shows
  no children.
- **The dedup key keeps its shape.** Hiding the `L`-phantom rule inside `Derivation`
  would change which mechanism is fitted first for a key, and so could change
  selections.
- **No `Identification` module.** Moving `identify_rate_equation.jl` into a module would
  remove no name and no argument from any interface.

#### Public names of `Derivation` after the change

From 16: `compile_rate_equation`, `rate_equation_string`, `parameters`,
`dependent_parameter_values` (these four exported), `parameter_units`,
`check_valid_rate_equation`, `fitted_parameter_count`, `rate_equation_dedup_key`,
`rate_equation_term_count`: 9.

## Tests

- **Pin first.** Before the data model moves, pin what must not change: the law of every
  `MECHANISM_TEST_SPECS` entry against its oracle, and the enumerator's populations.
- **Fixtures** translate mechanically from call notation to form names, by a script that
  applies today's `_form_name` (`E(B; residual = A - P)` becomes `EB_res_plusA_minusP`).
  It writes each step in the direction stored today, so every kinetic group stays aligned.
  The enumeration rule in `.claude/CLAUDE.md` ("write chemistry the enumerator's way") is
  rewritten in form names.
- **Enumeration tests keep their rules** (the maintainer's decision, 2026-11-20): every
  parent and every expected child is written inline with the macro, and the children,
  from the move or the method under test, are asserted exactly. A test helper,
  `_testhelper_candidate(mechanism)`, turns a written parent into a candidate by reading
  its form names under the enumerator's naming rule, and errors on a name it cannot read. It states the rule itself and does not call the
  package's. Children are compared by their mechanisms:
  `Set(c.mechanism for c in children) == Set(expected)`.
- **The oracles stay independent.** The test helpers that model an excluded step
  (`_testhelper_network`, `_testhelper_exclusion`, `_testhelper_epsilon_span`,
  `_testhelper_corner_directions`, `_testhelper_lp_existence`) keep their own statement
  of the table in decision 2 when they move to the lean graph; they do not call the
  package's.
- **No verdict changes.** The data model changes and the guard's rule does not. Every pin
  of a guard verdict, of an admitting reading and of an enumerator population passes
  after the fixtures are translated, with the same counts. A pin that would have to
  change stops the work and goes to the maintainer. The exceptions are named in this
  document: some fitted-parameter names (decision 5: the constants `_priority_key` picks,
  and `mod` leaving a name), `E_total` and `Keq` (decision 7), and the constants of the
  steps that decision 1 turns: 32 isomerizations and at most three exchanges, in about 25
  fixtures across eight test files. There a constant is renamed after the other direction
  and a rapid-equilibrium constant becomes its inverse; the law is the same function.
- **The gate script still runs.** `claude-scratchpad/lean-mechanism-spike/gate/` compares
  guard verdicts across ten populations. It is brought to the lean graph and rerun once
  at the end: every verdict must equal the one stage 1 gave.

## Documentation

- **The macro's docstring and the docs pages that write mechanisms** use the one macro
  and form names. Every `@allosteric_mechanism` and every call-notation form in
  `docs/src/`, `README.md` and the docstrings is translated, and the doctests pass.
- **Evaluating a law.** Every page that calls `rate_equation(m, concs, params)` shows
  `rate = compile_rate_equation(m)` and `rate(concs, params, Keq)`. `E_total` leaves the
  text, with one sentence where the rate is first defined: the law gives the rate per
  unit enzyme.
- **`docs/src/developer.md`:** the data model, the module tables and the diagram
  (regenerated by `julia test/update_module_map.jl`), what `Enumeration` keeps in a
  candidate, and the parameter interface.
- **`docs/src/identify/enumeration_engine.md`** as decision 6 says.
- **Islands.** `docs/src/developer.md` and `docs/src/deriving/mwc_allostery.md` call the
  sets of forms that a conformation's steps join "groups", which reads as kinetic groups.
  Both say "islands", as decision 2 does.
- **`.claude/CLAUDE.md`:** the module map and the bullets under it (the stored fields,
  the names reachable without a promise, `compile_mechanism`, the shared `E_total`); the
  Julia-style bullet on `EnzymeMechanism{Sig}`; the "`rate_equation` speed" guard; the
  rules for Derivation tests and enumeration-engine tests, which name
  `@allosteric_mechanism`, `expand_mechanisms` and the call notation.
- **`ROADMAP.md`:** the stored-facts known issue is measured again and retired or
  updated; "Shrink the module interfaces" is updated for `Mechanisms`, `Derivation` and
  `Enumeration`; a note says that only tests reach the enumerator's code for other
  conformation labels until iso mechanisms are enumerated (decision 6).

## Left to the builder

- The canonical order of forms, steps and kinetic groups. It may not depend on the
  writing, since the last row of the table of stored directions reads it.
- Whether `Sig` needs the integer encoding (decision 3).
- How the shortest turnover cycle through a step is found, and how a tie
  between two members of a group equally near the free enzyme is broken; neither may
  depend on the writing.
- The replacement rule in `_priority_key`.
- The private type that holds a candidate's record of its forms, and whether the moves
  work on it or on the lean graph beside it.
- The order of the two size checks and the guard at the enumerator's exits.
- Whether `fitted_parameter_count` folds into `length(parameters(m))`, by timing.

## Sequence

1. Stage 1, the guard's rule on today's data model: built.
2. The Enumeration design (decision 6) and the Derivation-interface design with `Keq` and
   `E_total` (decision 7): decided on 2026-11-20.
3. The rest is planned and built: the lean graph, the one type, the one macro with the
   fixtures, candidates in `Enumeration`, compiling inside `FittingProblem`, `Keq` as an
   argument, no `E_total`, the parameter interface. The plan chooses the order. Decision
   1 removes the compiled mechanism types and decision 3 supplies their replacement, so
   the two cannot be built apart; every commit passes the suite.

## Consequences for the maintainer to know

- **Stage 1's consequences** (search results, allosteric ping-pong, the numerical
  evidence, valid mechanisms still rejected) are recorded in `ROADMAP.md` and the docs.
- **Breaking changes for users:** `rate_equation(m, concs, params)` becomes
  `compile_rate_equation(m)(concs, params, Keq)`; `E_total` is gone, so a law returns the
  rate per unit enzyme; `parameters(m)` lists the fitted constants only, and the three
  kind types are gone; a fit's `params` no longer holds `Keq`;
  `dependent_parameter_values` takes `Keq` as an argument;
  `plot_rate_equation_on_data`, `test_rate_equation`, `FittingProblem` and
  `IdentifyRateEquationProblem` take `Keq` as their last positional argument;
  `catalytic_multiplicity` is `oligomeric_state`, `allowed_catalytic_multiplicities` is
  `allowed_oligomeric_states` and `mwc_tags` is `allosteric_states`; a mechanism is
  written with form names (`EA`, not `E(A)`), its free enzyme is named `E`, and `X::Mod`
  and `allosteric_regulators:` leave the mechanism macro; `@allosteric_mechanism`,
  `EnzymeMechanism` and `AllostericEnzymeMechanism` are gone.
- **The progress log** no longer counts the children skipped for each size limit
  (decision 7).
- **An enumerator bug in the starting set** stops a run at its start (decision 6).
- **Fitted-parameter names** change for some mechanisms (decision 5), and the constants
  of the steps decision 1 turns are named after the other direction.
- **A kinetic group of isomerizations** is aligned by its writer; today the constructor
  aligns it by the conformation labels (decision 1).
- **Saved mechanism text** from earlier versions does not parse (decision 4).
- **The decodability check** moves from the `EnzymeReaction` constructor and the
  mechanism macro to the start of enumeration (decision 6).
- **All of the refactor lands in one pull request** (the maintainer's decision,
  2026-11-20): stage 1 and everything this document specifies.
- **Version.** The change is breaking: 0.11.0, bumped once.
