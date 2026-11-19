# One lean `Mechanism`: design

Date: 2026-11-18. Status: approved by the maintainer on 2026-11-19. Nothing is built yet.

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
   steady-state step in the other.

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
absent, not absent. The science it assumes is unchanged: an excluded step has a very high
dissociation constant, or very low rate constants, and every thermodynamic constraint
still holds in each conformation. What changes is how the guard chooses between the two
for a step. Today it reads the form contents. Under this rule it reads the graph, and a
step keeps its kind (rapid equilibrium or steady state) in both conformations:

| Excluded step | How it is switched off |
|:--|:--|
| Rapid equilibrium, takes up one metabolite into a cut-off form | It stays fast and loses its affinity ("K very high"). |
| Rapid equilibrium, anything else | Fault. Thermodynamics pins its equilibrium constant, so it could be switched off only by turning slow. |
| Steady state, takes up one metabolite into a cut-off form | It loses its affinity: the uptake rate vanishes and the release keeps its rate. |
| Steady state, anything else | It is stopped: both rates vanish, and its equilibrium constant stays what thermodynamics requires. |

"Into a cut-off form" means: the bound side is cut off, and the steps that still run do
not join the two sides.

Why a step into a cut-off form must lose its affinity: the law assumes the cut-off form is
empty. With both rates low and the affinity unchanged, the form would keep its equilibrium
share however slow the step, and the law that drops it would be wrong.

Why thermodynamics holds: each set of forms that the running steps keep joined has one
energy level. A step that loses its affinity raises the level of its bound side, and every
excluded step's equilibrium constant shifts by the difference between its two ends'
levels. The shifts cancel around every cycle.

Three examples, each in the inactive conformation, the unlisted steps `:EqualAI`:

| Mechanism | Verdict |
|:--|:--|
| `E + S ⇌ ES :: OnlyA`, `ES <--> EP :: OnlyA`, `E + P ⇌ EP` | Admitted. `ES` is cut off: S binding loses its affinity, the chemistry is stopped. |
| `E + S ⇌ ES :: OnlyA`, `ES <--> EP`, `E + P ⇌ EP` | Fault. `ES` stays, made from P, so thermodynamics pins the affinity of S, and a rapid-equilibrium step cannot be stopped. |
| `E + S <--> ES :: OnlyA`, `ES <--> EP`, `E + P ⇌ EP` | Admitted. The steady-state binding is stopped: S neither binds nor leaves, and `ES` exists, made from P. Enumeration does not build this one (decision 5). |

#### The checks

With each excluded step read that way, the guard's checks run as today: the steps that
lose their affinity close no directed cycle, free E reaches every cut-off part along them,
the active conformation still catalyzes the reaction, and the face condition holds on
every face. Against today's guard:

| Case | Today | This rule |
|:--|:--|:--|
| Rapid-equilibrium plain binding whose two sides stay joined | Fault (`_NO_LIMIT`) | Fault |
| Steady-state plain binding whose two sides stay joined | Fault (`_NO_LIMIT`) | Admitted: stopped |
| Rapid-equilibrium step that binds and reacts at once, or isomerizes | Admitted: fast in one conformation, stopped in the other | Fault |
| Step that binds and reacts at once into a cut-off form | Stopped, so the cut-off form is an island | Loses its affinity |

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

Not checked: ter-ter; and an excluded rapid-equilibrium isomerization that could stay
fast with its equilibrium pushed to the side that remains, which the rule rejects for now.

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

  The tests of a law that is not the limit pass unchanged: a limit that depends on how
  fast an isomerization vanishes, a pumped basin, a bridge, a second closed class, a trap,
  a split, the conservative rejections. They did in the spike.
- **The guard rule is a gate.** Stage 1 lands and passes before any fixture is
  translated. The work stops and returns to the maintainer if a mechanism admitted today
  is rejected for any reason but an excluded rapid-equilibrium step that would have to
  stop, or if a newly admitted law is not the numerical limit. The data model leaves no
  room for a per-step flag.
- **What counts as checked.** The comparisons of decision 2 are rerun once on the
  implemented rule and reported in the pull request, with a sample of at least 2,000
  ter-ter mechanisms added. The suite keeps a small permanent version: a sampled
  population test that the only mechanisms lost are the ones the rule names, and the
  oracle agreement on a sample of the newly admitted ones.

## Documentation

The assumption behind an excluded step is explained to users, not only to developers (the
maintainer's request). It lands with stage 1, when the rule changes:

- **A user-facing section, "How an excluded step is switched off",** on the docs page that
  introduces allosteric mechanisms. It gives, in the user's terms: what `:OnlyA` and
  `:OnlyI` assume (a very high dissociation constant or very low rate constants, with
  every thermodynamic constraint kept in each conformation); that a step keeps its kind,
  rapid equilibrium or steady state, in both conformations; the four-row table of
  decision 2; the three worked examples; and what the error means when a
  rapid-equilibrium step cannot be switched off, with the two fixes (make the step
  steady-state, or exclude the steps that keep its bound form reachable).
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
   no Enumeration internals, and is the gate.
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
  rule some of these are admitted, because an excluded first half-reaction that takes up
  one metabolite into cut-off forms loses its affinity. The entry is re-examined and
  rewritten as part of the work.
- **Fitted-parameter names** change for some mechanisms (decision 5).
- **Saved mechanism text** from earlier versions does not parse (decision 4).
- **The decodability check** moves from `@enzyme_reaction` to the start of enumeration
  (decision 5, provisional).
- **Version.** The change is breaking: 0.11.0, bumped once.
- **`ROADMAP.md`.** The stored-facts known issue is measured again and retired or updated;
  "Shrink the module interfaces" is updated for `Mechanisms`.
