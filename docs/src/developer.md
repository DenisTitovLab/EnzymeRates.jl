# Developer / Architecture

This page describes how EnzymeRates works inside, for maintainers; the [API Reference](@ref)
lists the exported names. A mechanism is a plain value, `Mechanism`, everywhere; only its
compiled rate law has a type of its own (see [Derivation architecture](@ref)).

## Modules

EnzymeRates is seven submodules and a top level. Each submodule lives in `src/<Module>/`:
`<Module>.jl` declares the module, includes its files and lists its interface in a `public`
statement. A module's public names are its interface to the other modules; only the exports
are promised to users. Its tests live in `test/<Module>/`, a file `test_<file>.jl` for each
source file `<file>.jl`. The diagram below shows every exported name and every name one
module takes from another, and the tables list every module's public names;
`test/test_boundaries.jl` fails when they disagree with `src/`. The diagram is generated:
each box's sentence is the first sentence of its module's docstring, and
`julia test/update_module_map.jl` rewrites it, together with the table of files below and
the table of modules in `.claude/CLAUDE.md`, which the same test checks.

### The module map

```mermaid
%%{init: {"theme": "base", "flowchart": {"curve": "basis", "nodeSpacing": 40, "rankSpacing": 70, "htmlLabels": true},
  "themeVariables": {"fontFamily": "Helvetica, Arial, sans-serif", "fontSize": "13px",
  "primaryColor": "#ffffff", "primaryBorderColor": "#666666", "primaryTextColor": "#222222",
  "lineColor": "#777777", "edgeLabelBackground": "#ffffff", "clusterBkg": "#eaf2ee",
  "clusterBorder": "#666666", "titleColor": "#1f4f42"}}}%%
flowchart TB
  top["<b><big>EnzymeRates</big></b><br/>Finds the best mechanism for a reaction from measured rates.<br/><br/>Defines:<br/><b>★#nbsp;identify_rate_equation()</b><br/><b>★#nbsp;IdentifyRateEquationProblem</b><br/><b>★#nbsp;IdentifyRateEquationResults</b>"]
  dsl["<b><big>DSL</big></b><br/>Turns macro-written reactions and mechanisms<br/>into checked objects."]
  plot["<b><big>Plotting</big></b><br/>Draws a rate equation over measured rates."]
  enum["<b><big>Enumeration</big></b><br/>Generates candidate mechanisms for a reaction."]
  fit["<b><big>Fitting</big></b><br/>Fits a rate law to measured rates."]
  der["<b><big>Derivation</big></b><br/>Derives a mechanism's rate equation and compiles it."]
  data["<b><big>Data</big></b><br/>Checks and holds a table of measured rates."]
  mech["<b><big>Mechanisms</big></b><br/>Stores reactions and mechanisms as plain data.<br/><br/>Types every module but Data reads:<br/><b>★#nbsp;EnzymeReaction</b><br/><b>★#nbsp;Mechanism</b><br/>Allostery<br/>CatalyticScheme<br/>Reactant<br/>Regulator<br/>RegulatorySite<br/>Step"]
  top -->|"<b>★#nbsp;@enzyme_reaction</b><br/><b>★#nbsp;@enzyme_mechanism</b>"| dsl
  top -->|"<b>★#nbsp;plot_rate_equation_on_data()</b>"| plot
  top -->|"CandidateMechanism<br/>candidate_mechanisms()"| enum
  top -->|"<b>★#nbsp;fit_rate_equation()</b><br/><b>★#nbsp;FittingProblem</b><br/>test_rate_equation()"| fit
  top -->|"<b>★#nbsp;compile_rate_equation()</b><br/><b>★#nbsp;rate_equation_string()</b><br/><b>★#nbsp;parameters()</b><br/><b>★#nbsp;dependent_parameter_values()</b><br/>rate_equation_dedup_key()"| der
  top -->|"RateData"| data
  top -->|"<b>★#nbsp;metabolites()</b>"| mech
  dsl -->|"check_valid_rate_equation()"| der
  plot -->|"<b>★#nbsp;compile_rate_equation()</b><br/><b>★#nbsp;parameters()</b>"| der
  plot -->|"RateData"| data
  plot -->|"<b>★#nbsp;metabolites()</b>"| mech
  enum -->|"check_valid_rate_equation()<br/>fitted_parameter_count()<br/>rate_equation_term_count()"| der
  enum -->|"has_nonzero_weight_cycle()<br/>has_turnover_cycle()<br/>rapid_equilibrium_segments()<br/>step_text()<br/>steps_running_in_allosteric_conformations()"| mech
  fit -->|"<b>★#nbsp;compile_rate_equation()</b><br/><b>★#nbsp;parameters()</b><br/>parameter_units()"| der
  fit -->|"RateData"| data
  fit -->|"<b>★#nbsp;metabolites()</b>"| mech
  der -->|"<b>★#nbsp;metabolites()</b><br/>has_nonzero_weight_cycle()<br/>has_turnover_cycle()<br/>rapid_equilibrium_segments()<br/>step_text()<br/>steps_running_in_allosteric_conformations()"| mech
  classDef mod fill:#eaf2ee,stroke:#666666,color:#222222;
  class top,dsl,plot,enum,fit,der,data,mech mod;
```

★ bold: exported by EnzymeRates. Arrows carry every name one module takes from another;
the Mechanisms types are listed in its box.

An arrow points from a module to a module whose names it uses. A submodule names only the
modules its arrows reach; `Mechanisms` and `Data` name none. Inside a module, code calls
its own names bare; a call into another module is qualified:
`Derivation.check_valid_rate_equation(m)`.

### What each module does

**`Data`** checks a rate table once and holds it. A `RateData` keeps the table's columns
with `Keq`, the `absolute` flag and the reference anchors; `FittingProblem` and
`IdentifyRateEquationProblem` each hold one. `Data` hides the column and anchor checks.

**`Mechanisms`** is the data model. It holds reactions (`EnzymeReaction`) and mechanisms
(`Mechanism`: a `CatalyticScheme` and its `Allostery`) as plain values whose fields any
module reads with dots. Its constructors put a mechanism in canonical form and check it,
and `show` prints a mechanism as the macro block that builds it. It hides canonical form
and the representation checks. A mechanism stores its graph and nothing derived from it;
`Mechanisms` answers the graph analyses that two modules read: the steps each
conformation runs, rapid-equilibrium segments, turnover cycles and cycles of nonzero
weight.

**`Derivation`** turns a mechanism into its rate law. It decides whether the law is
defined (the guard), solves the thermodynamic constraints, names the parameters, and
derives the law with its generated code, printed text and key. It hides the King–Altman/Cha
engine, the polynomials, the joint solve and a memo that derives each mechanism once for
several questions.

**`Fitting`** fits a rate law to data. It hides the loss, the multi-start restarts, the
scaling to the reference and the call into Optimization.jl.

**`Enumeration`** proposes candidate mechanisms for a reaction (`candidate_mechanisms`):
the starting candidates identify fits first, and the children of each candidate it keeps,
each within identify's two size limits. It hides the seeds, the expansion moves, the
structural filters that prune the search and the check of the reaction's names.

**`DSL`** reads the macro language: `@enzyme_reaction`, `@enzyme_mechanism`, which
writes plain and allosteric mechanisms alike, and `parse_mechanism`, which reads a
mechanism's printed text back. It hides the parser. `@enzyme_mechanism` runs the guard,
so the DSL cannot live in `Mechanisms`.

**`Plotting`** draws a rate equation over the rates it was fitted to, with
`plot_rate_equation_on_data`. It hides the preparation of each panel: the choice of x
axis, the series, the curves and their scaling, and the labels. The package extension
`ext/EnzymeRatesMakieExt.jl`, which Julia loads with Makie, draws the panels, so loading
EnzymeRates alone never loads Makie. The extension defines
`draw_rate_equation_panels(panels; columns, figure, axis)`, which
`plot_rate_equation_on_data` calls through `Base.get_extension` with the panels it prepared.

**The top level** (`src/EnzymeRates.jl`, `src/identify_rate_equation.jl`) holds
`identify_rate_equation` and the exports. It hides the beam search, the one fit per rate
law, cross-validation and the CSV output.

### Module docstrings

The first sentence of a module's docstring is the sentence in its box of the diagram.

```@docs
EnzymeRates
EnzymeRates.Data
EnzymeRates.Mechanisms
EnzymeRates.Derivation
EnzymeRates.Fitting
EnzymeRates.Enumeration
EnzymeRates.DSL
EnzymeRates.Plotting
```

### Interfaces

Each table lists a module's public names; the top level's lists the names it defines and
exports. "Used by" names the modules that call each name, and ✓ marks a name that
`using EnzymeRates` brings in. The tests call more.

**`Data`**

| Name | What it does | Used by | Exported |
|:-----|:-------------|:--------|:---------|
| `RateData` | A checked rate table with `Keq` and its anchors | Fitting, Plotting, top level | |

**`Mechanisms`**

| Name | What it does | Used by | Exported |
|:-----|:-------------|:--------|:---------|
| `EnzymeReaction` | A reaction and its regulators | DSL, Enumeration, top level | ✓ |
| `Reactant` | A substrate or product with its atoms | DSL | |
| `Regulator` | A declared regulator with its kind and effect | DSL | |
| `Mechanism` | A catalytic scheme and its allostery, in | DSL, Derivation, | ✓ |
| | canonical form | Enumeration, Fitting, | |
| | | Plotting, top level | |
| `CatalyticScheme` | A reaction graph: metabolites, forms and steps | DSL, Derivation, | |
| | | Enumeration | |
| `Allostery` | A mechanism's allosteric states, its oligomeric | DSL, Derivation, | |
| | state and its sites | Enumeration | |
| `Step` | A step between two enzyme forms | DSL, Derivation, Enumeration | |
| `RegulatorySite` | A regulatory site | DSL, Derivation, Enumeration | |
| `metabolites` | The rate law's concentration names | Derivation, Fitting, | ✓ |
| | | Plotting, top level | |
| `step_text` | A step as the macro writes it: `E + A ⇌ EA` | Derivation, Enumeration | |
| `steps_running_in_allosteric_conformations` | The steps each conformation runs | Derivation, Enumeration | |
| `rapid_equilibrium_segments` | The rapid-equilibrium segments of a set of steps, and the | Derivation, Enumeration | |
| | steady-state steps between them | | |
| `has_turnover_cycle` | Whether a cycle of steps turns the reaction over | Derivation, Enumeration | |
| `has_nonzero_weight_cycle` | Whether the steps a mask keeps close a cycle whose weights do | Derivation, Enumeration | |
| | not sum to zero | | |

**`Derivation`**

| Name | What it does | Used by | Exported |
|:-----|:-------------|:--------|:---------|
| `compile_rate_equation` | The rate law as a function, allocation-free per call | Fitting, Plotting, top level | ✓ |
| `rate_equation_string` | The rate law as text | top level | ✓ |
| `parameters` | The names of the constants a fit estimates | Fitting, Plotting, top level | ✓ |
| `dependent_parameter_values` | The dependent constants' values from the fitted ones and `Keq` | top level | ✓ |
| `parameter_units` | Each fitted constant's units, as powers of the rate and concentrations | Fitting | |
| `check_valid_rate_equation` | Whether the law is defined | DSL, Enumeration | |
| `fitted_parameter_count` | The number of fitted constants, from the joint solve | Enumeration | |
| `rate_equation_term_count` | The law's term count, V×τ, before derivation | Enumeration | |
| `rate_equation_dedup_key` | The law's key: one key, one fit | top level | |
| `RateEquationTooLargeError` | The error a law over the term cap throws | | |

**`Fitting`**

| Name | What it does | Used by | Exported |
|:-----|:-------------|:--------|:---------|
| `FittingProblem` | A mechanism, its compiled rate law and its `RateData` | top level | ✓ |
| `fit_rate_equation` | Fits a problem with an optimizer | top level | ✓ |
| `test_rate_equation` | The loss of given parameters on a rate table | top level | |

**`Enumeration`**

| Name | What it does | Used by | Exported |
|:-----|:-------------|:--------|:---------|
| `CandidateMechanism` | A mechanism the enumerator proposes, read as `c.mechanism` | top level | |
| `candidate_mechanisms` | The starting candidates of a reaction, or a candidate's children | top level | |

**`DSL`**

| Name | What it does | Used by | Exported |
|:-----|:-------------|:--------|:---------|
| `@enzyme_reaction` | Builds an `EnzymeReaction` | top level | ✓ |
| `@enzyme_mechanism` | Builds and checks a mechanism, plain or allosteric | top level | ✓ |
| `parse_mechanism` | Reads a mechanism's printed text back | | |

**`Plotting`**

| Name | What it does | Used by | Exported |
|:-----|:-------------|:--------|:---------|
| `plot_rate_equation_on_data` | Draws a rate equation over its data, a panel per group | top level | ✓ |

**The top level**

| Name | What it does | Used by | Exported |
|:-----|:-------------|:--------|:---------|
| `IdentifyRateEquationProblem` | A reaction with its checked `RateData` | | ✓ |
| `IdentifyRateEquationResults` | The selected mechanism and the CV table | | ✓ |
| `identify_rate_equation` | Searches mechanisms, selects one by CV | | ✓ |

### Files

Three rules govern the source files of a module:

1. A file holds one piece of functionality whole, so that a reader finds what belongs
   together in one place.
2. A file is named for what it gives the rest of the package: the exported or public
   function, type or macro it defines (`fit_rate_equation.jl`, `mechanism_type.jl`,
   `enzyme_mechanism_macro.jl`), or, when no such name says what the file does, the
   calculation it carries out (`rate_equation_derivation.jl`,
   `rate_equation_correctness_guard.jl`). A name then leads to its file, which `utils.jl`,
   `helpers.jl` or `types.jl` never do.
3. Calls between the files of a module are few and one-way: a file does not use a file
   that uses it, directly or through other files. A reader then understands a file from
   the files it uses, and a change to a file can break only the files that use it.

The table below is generated from `src/`; "What it holds" is each file's two `ABOUTME`
lines. In the row of `starting_mechanisms.jl`, `expansion_moves` (6) under "Takes names
from" says that the file uses six names that `expansion_moves.jl` defines. Two pairs of
`Derivation` files use each other, which joins three files in loops; `ROADMAP.md` lists
both pairs under "Known issues", with the use in each that would have to go.

```@raw html
<!-- file table, generated from src/ by test/update_module_map.jl -->
```

| Module | File | What it holds | Public names it defines | Takes names from |
|:--|:--|:--|:--|:--|
| top level | `identify_rate_equation.jl` | identify\_rate\_equation, its problem and its results: the beam search over the candidates, one fit per rate-law key, selection by cross-validation, the CSVs. | ★ `identify_rate_equation`, ★ `IdentifyRateEquationProblem`, ★ `IdentifyRateEquationResults` |  |
| DSL | `enzyme_reaction_macro.jl` | @enzyme\_reaction, which builds an EnzymeReaction, and the reader of its block: reactants with their atom brackets, regulators, oligomeric states, shared sites. | ★ `@enzyme_reaction` |  |
| DSL | `enzyme_mechanism_macro.jl` | @enzyme\_mechanism and parse\_mechanism, which build a mechanism from its text and check it with the guard, and the reader of a mechanism block as written. | ★ `@enzyme_mechanism`, `parse_mechanism` |  |
| Plotting | `plot_rate_equation_on_data.jl` | plot\_rate\_equation\_on\_data and the plain description of each measurement group's panel it prepares; ext/EnzymeRatesMakieExt.jl draws the panels with Makie. | ★ `plot_rate_equation_on_data` |  |
| Enumeration | `candidate_mechanism_type.jl` | CandidateMechanism, a mechanism with the enumerator's record of each form, the structural tests of a candidate, and the helpers the moves and the seeds share. | `CandidateMechanism` |  |
| Enumeration | `expansion_moves.jl` | A candidate's children, candidate\_mechanisms(reaction, parent), and the moves that make them: RE→SS flips, group splits, regulators, states, site merges. | `candidate_mechanisms` | `candidate_mechanism_type` (23) |
| Enumeration | `starting_mechanisms.jl` | The starting mechanisms of a reaction, candidate\_mechanisms(reaction): catalytic topologies, dead-end decorations, merged and Theorell–Chance variants, seeds. | `candidate_mechanisms` | `candidate_mechanism_type` (22), `expansion_moves` (6) |
| Fitting | `fit_rate_equation.jl` | Fitting rate equations to experimental data and scaling them to a reference scale: FittingProblem, the loss, multi-start fits and test\_rate\_equation. | ★ `fit_rate_equation`, ★ `FittingProblem`, `test_rate_equation` |  |
| Derivation | `thermodynamic_constraints.jl` | A mechanism's constants under the joint thermodynamic solve: their types and names, the solve, parameters, their units and the dependent constants' values. | ★ `parameters`, ★ `dependent_parameter_values`, `fitted_parameter_count`, `parameter_units` | `rate_equation_correctness_guard` (1) |
| Derivation | `rate_equation_derivation.jl` | A mechanism's rate law: sparse polynomials, the law on edge masks and the MWC law, its code (compile\_rate\_equation), its text and the unit rule it prints by. | ★ `compile_rate_equation`, ★ `rate_equation_string`, ★ `parameters`, `RateEquationTooLargeError`, `rate_equation_term_count` | `thermodynamic_constraints` (5), `rate_equation_correctness_guard` (3) |
| Derivation | `rate_equation_dedup_key.jl` | The canonical key of a rate law, by colour refinement and individualization, and the rescaling certificate that flags a one-subunit L phantom. | `rate_equation_dedup_key` | `thermodynamic_constraints` (1), `rate_equation_derivation` (7), `rate_equation_correctness_guard` (1) |
| Derivation | `rate_equation_correctness_guard.jl` | The guard against a wrong rate law, check\_valid\_rate\_equation: its verdicts, networks and faces, classes 1 to 3 with the limit test, and its error messages. | `check_valid_rate_equation` | `thermodynamic_constraints` (1), `allosteric_rate_equation_correctness_guard` (8) |
| Derivation | `allosteric_rate_equation_correctness_guard.jl` | Admission, the guard's check that an allosteric law is the limit as the excluded steps are switched off: the readings, its interior and face conditions, messages. |  | `rate_equation_correctness_guard` (31) |
| Data | `rate_data_type.jl` | RateData: a rate table validated against a metabolite list, with its Keq and the reference anchors that set a relative fit's scale. | `RateData` |  |
| Mechanisms | `enzyme_reaction_type.jl` | EnzymeReaction, a reaction as a value: its reactants with their atoms, its regulators, the checks of metabolite names, and its text as the macro block. | ★ `EnzymeReaction`, `Reactant`, `Regulator` |  |
| Mechanisms | `mechanism_type.jl` | Mechanism, a reaction graph in canonical form and its allostery, with their checks, identity and text, and what the graph gives: steps, segments, cycles. | ★ `metabolites`, ★ `Mechanism`, `Allostery`, `CatalyticScheme`, `RegulatorySite`, `Step`, `has_nonzero_weight_cycle`, `has_turnover_cycle`, `rapid_equilibrium_segments`, `step_text`, `steps_running_in_allosteric_conformations` | `enzyme_reaction_type` (5) |

```@raw html
<!-- end of file table -->
```

### The mechanism data model

Stored data is read with dots; derived data is a function. A catalytic scheme is a list
of metabolite names, a list of enzyme-form names and a list of steps. A metabolite, a form
or a step is an integer index into `c.metabolites`, `c.enzyme_forms` or `c.steps`, and a
set of steps is a `BitVector` over `c.steps`. Code reads `m.catalytic.enzyme_forms[f]`,
`m.catalytic.steps[k].from` or `m.allosteric.allosteric_states[g]` directly, and calls a
function only for a fact the fields do not store, such as `step_text(m.catalytic, k)`. A
mechanism holds no reaction, and `==` and `hash` read its fields alone, alike in every
process.

| Type | Field | Holds |
|:-----|:------|:------|
| `Step` | `from::Int`, `to::Int` | its two forms, in its stored direction |
| | `taken_up::Vector{Int}` | the metabolites it takes up going `from → to` |
| | `given_off::Vector{Int}` | the metabolites it gives off going `from → to` |
| | `rapid_equilibrium::Bool` | whether it is at rapid equilibrium |
| | `kinetic_group::Int` | its kinetic group |
| `Mechanism` | `catalytic::CatalyticScheme` | the catalytic scheme |
| | `allosteric::Allostery` | its allostery; a plain mechanism's has no states |
| `CatalyticScheme` | `metabolites::Vector{Symbol}` | the substrates, the products, then |
| | | the modifiers, each group by name |
| | `enzyme_forms::Vector{Symbol}` | free E first, then the other forms by name |
| | `substrates::Vector{Int}` | the substrates' indices into `metabolites` |
| | `products::Vector{Int}` | the products' indices into `metabolites` |
| | `steps::Vector{Step}` | each step in the direction the graph gives, sorted |
| `Allostery` | `allosteric_states::Vector{Symbol}` | one per kinetic group: `:OnlyA`, |
| | | `:OnlyI`, `:EqualAI` or `:NonequalAI`; |
| | | empty for a plain mechanism |
| | `oligomeric_state::Int` | the number of subunits that change conformation together |
| | `regulatory_sites::Vector{RegulatorySite}` | the regulatory sites |
| `RegulatorySite` | `ligands::Vector{Symbol}` | the ligands that bind the site |
| | `multiplicity::Int` | the number of such sites |
| | `allosteric_states::Vector{Symbol}` | one allosteric state per ligand |
| `Reactant` | `name::Symbol` | the substrate's or product's name |
| | `atoms::Vector{Pair{Symbol,Int}}` | its atoms, `element => count` |
| `Regulator` | `name::Symbol` | the regulator's name |
| | `kind::Symbol` | `:allosteric` or `:competitive` |
| | `effect::Symbol` | `:activator`, `:inhibitor` or `:unspecified` |
| `EnzymeReaction` | `substrates::Vector{Reactant}` | the substrates |
| | `products::Vector{Reactant}` | the products |
| | `regulators::Vector{Regulator}` | the declared regulators |
| | `allowed_oligomeric_states::Vector{Int}` | the subunit counts to enumerate |
| | `shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}` | pairs that share one site |

A form's name is opaque: the scheme stores no form contents, and only `E`, the free enzyme,
means something. So a reactant bound where it does not react is the reactant taken up, and
only the forms' names tell the two bindings apart (`ESSmod` beside `ES`). Two different
sides can also render one constant name, as a form named `ES_I` and the side `ES + I` do;
the guard rejects such a mechanism (`_shared_constant_names`). The derivation and the guard
tell a step's shape from its lists alone. A binding takes up one metabolite and gives off
none, and an isomerization takes up and gives off nothing. A Theorell–Chance step takes up
one metabolite and gives off another (`k_EA_B_to_EQ_P` for `EA + B <--> EQ + P`). The
derivation reads every binding alike. So a fused binding, straight into the product-bound
form (`EA + B ⇌ EPQ`), has a dissociation constant at rapid equilibrium (`K_EPQ_to_EA_B`),
and a binding and a release rate constant at steady state.

The constructors check the representation and store the mechanism in canonical form. The
docstrings of `CatalyticScheme` and `Mechanism` say what each takes and throws, `Mechanism`
checking the `Allostery` it is given, and `_orientation` gives the rules for a step's stored
direction. The canonical form is load-bearing. Deduplication by `==` and `hash` works only
because construction is canonical, and every constant is named after the stored sides of its
class's kept member, so do not relax or reorder it. The steps of a kinetic group share their
constants, so they turn together with the group's lead (`_turn_groups_together!`), and an
isomerization keeps the direction it was written in relative to the lead. In
`(Estar <--> E, EstarA <--> EA)` the reaction runs through the two steps in opposite
directions, so nothing in the graph says how the members correspond: the writer aligns the
group, as the enumerator does, and turning one isomerization against its group writes
another mechanism. The steps of a group take up and give off the same metabolites and carry
one flag (`_check_uniform_groups`), so a fused and a plain binding of one metabolite may
share a group, and a Theorell–Chance step cannot share one with a binding.

To change a mechanism, copy `m.catalytic.enzyme_forms` and `m.catalytic.steps` with new
fields and call the constructors. Never mutate a built mechanism's vectors in place: its
canonical form, hash, compiled Sig and memo entries depend on them.

### One request through the modules

`identify_rate_equation` crosses every module on its way from a reaction to the selected
fit:

1. **The problem.** `@enzyme_reaction` (`DSL`) builds the `EnzymeReaction`
   (`Mechanisms`). `IdentifyRateEquationProblem(reaction, table, Keq)` (top level) checks
   the table into a `Data.RateData`.
2. **The starting mechanisms.** `Enumeration.candidate_mechanisms(reaction; ...)`, given
   identify's size limits and optional regulators, writes each seed as forms and steps,
   and the `Mechanism` constructor puts it in canonical form. It returns only candidates
   whose rate law is defined (`Derivation.check_valid_rate_equation(Bool, m)`) and within
   the size limits (`Derivation.rate_equation_term_count`,
   `Derivation.fitted_parameter_count`).
3. **One fit per law.** The top level asks `Derivation` three questions about each new
   candidate's mechanism: `parameters(m)` names its fitted constants,
   `rate_equation_string` prints the law, and `rate_equation_dedup_key` gives its key.
   Mechanisms with one key share one law, so one mechanism per key is fitted.
4. **The fit.** `FittingProblem` (`Fitting`) compiles that mechanism's rate law once
   (`Derivation.compile_rate_equation`), and `fit_rate_equation` fits it. The loss calls
   the compiled law on every data row at every evaluation.
5. **The beam.** The top level keeps the laws whose loss is close to the best at their
   parameter count, and `Enumeration.candidate_mechanisms(reaction, parent; ...)` proposes
   each candidate's children. Steps 3 to 5 repeat until no new law qualifies.
6. **The selection.** The top level cross-validates the best candidates at each
   parameter count: it refits each candidate once per data group, holding that group out
   (`fit_rate_equation`), and scores the fit on the held-out group (`test_rate_equation`).
   The 1-SE rule picks the result. The CSVs store each mechanism as its printed text
   (`show`), which `DSL.parse_mechanism` can read back.

### The boundary test

`test/test_boundaries.jl` parses `src/` with `Meta.parseall`, without loading the
package, and fails on:

1. a name of another module, as `Mod.name` or imported by `using`, that is not public in
   that module;
2. a module naming a module it may not depend on;
3. an underscore name in a `public` statement or used across modules;
4. an `export` in a module, which would let `using ..Mod` bring in bare names.

To use another module's function, it must be public there, under a name without an
underscore. The test also parses this page and fails when an interface table under
"Modules" disagrees with `src/`, naming the module and the names to add or remove.

The test also checks the test files' names. A test file is `test_<file>.jl` for the one
source file `<file>.jl` whose functions its testsets are about, in the mirrored folder:
`src/<Module>/` → `test/<Module>/`, `src/` → `test/`, `ext/` → `test/ext/`. A test file
directly in `test/` may instead check a package-wide rule (`test_boundaries.jl`,
`test_type_stability.jl`, `test_compile_budget.jl`, `test_aqua_jet.jl`, `test_oracles.jl`).
The test fails on a test file that names no source file; on a file in a module's test
folder or `test/ext/` that is not named `test_<file>.jl` or sits in a folder below it; on a
folder of `test/` that mirrors no source folder (`ext`, `fixtures` and `reference` aside);
and on a source file other than a module's header that has no test file.

The test also checks the source files, module headers aside, against the rules under
"Files". It fails when a file calls `include`, which only a header does; when a file's name
says nothing of what it holds (`utils.jl`, `helpers.jl`, `types.jl` and the like); when a
file is named for a public name of its module, or an export at the top level, that it does
not define; when a file `<x>_type.jl` defines no type whose name is `<x>` in snake case, or
a file `<x>_macro.jl` no macro `@<x>`; when a file of a module takes names from a file that
uses it, directly or through other files, once the uses the test lists as known exceptions
are set aside, naming the loop's files in order and the names each takes from the next; when
a listed use is not needed to break a loop, because it is gone or because the other listed
uses already break every loop through it (of uses that break the same loops, the test keeps
the last in sorted order and flags the rest); and when the table of files differs from the
one generated from `src/`. A file takes a name from another file of its module when that
file defines the name at top level, the file does not, and the file's code uses it; a name
in a string, comment, docstring or quote, a field after a dot, and a keyword or named-tuple
field name are no use. Last, the test fails on a `.jl` file under `src/`, `ext/`, `test/` or
`docs/` that does not start with exactly two `# ABOUTME: ` lines.

## Derivation architecture

A mechanism's rate equation is derived once, at compile time, by the
King–Altman/Cha method in `@generated` functions
(`src/Derivation/rate_equation_derivation.jl`). To make that possible,
`compile_rate_equation(m)` returns a **singleton callable**, `_CompiledRateEquation{Sig}`,
whose one type parameter, `Sig`, encodes the whole mechanism in canonical form, its
catalytic scheme and its allostery, as nested tuples of Symbols, Ints and Bools (`_sig_of`),
so the order or direction in which the steps were written does not change the resulting
type. The compiler specializes the derivation per type, so the symbolic King–Altman work
runs once per mechanism and session rather than at every call. The mechanism itself is a
plain value; only its compiled law carries a type of its own.

This is the most important architectural decision in the package, and it is deliberate.
Moving the derivation to compile time leaves the compiled law,
`rate = compile_rate_equation(m)` called as `rate(conc, params, Keq)`, as straight-line
arithmetic, written in Horner form (`_horner`), that must be **allocation-free and fast**;
the speed gate ([What a fit costs](@ref)) holds every gated law to that.
`rate_equation_string` (`src/Derivation/rate_equation_derivation.jl`) prints the law rather
than the code: each denominator term over the unit rule's normalizer, dimensionless whenever
the law has one. The compiled law's speed is the binding constraint on the whole package:
the fitter is a multi-start, global, gradient-free optimizer that evaluates the compiled law
millions of times per fit, and a single rate equation can take minutes to fit, so any
per-call allocation or microsecond-scale overhead would make fitting — and therefore
`identify_rate_equation`, which fits thousands of candidates — impractical. `loss!` is held
to the same standard: `FittingProblem` pre-allocates its `log_ratios_buffer` once and
`loss!` reuses it, so the inner optimization loop allocates nothing. The cost is the flip
side of the benefit: each unique compiled law triggers a full symbolic derivation, and a
mechanism with many enzyme forms or steps can be slow to compile, exhaust memory, or
`StackOverflow`.

The encoding of a mechanism in a type parameter lives in `Derivation`, behind
`compile_rate_equation`: `_sig_of` writes the Sig, and the generated methods decode it
(`_derivation_input`, `_mechanism_from_sig`) through the `Mechanisms` constructors, which
check it, so a mechanism survives the round trip unchanged. Equal mechanisms give the same
callable. Only the `@generated` methods do work per compiled law: the callable's call
runs the derivation, `parameters(rate)` runs its constraint solve, and `metabolites(rate)`
reads the names off the decoded mechanism. The decoder reads the Sig at run time, and the
package image holds it, so it compiles once for all mechanisms. Each compiled law then
pays for its generated bodies and the derivation they run, and for the methods that
specialize on it: Base's `text/plain` display and `FittingProblem` and `loss!` below. A
compiled law prints one line (`compiled rate equation`), whatever its Sig. A mechanism
prints as the macro block that builds it (`show`, `src/Mechanisms/mechanism_type.jl`);
`DSL.parse_mechanism` reads that text back through the macro's walker and checks, never
evaluating it, and identify's CSVs store it with a short id (`mechanism_id`).
`FittingProblem` holds the mechanism and its compiled law, and `loss!` specializes on the
law's type on purpose: it then reads `parameters(rate)` and `metabolites(rate)` as
constants and inlines the law into the loop that computes the predictions. The optimizer
never sees that type.
`fit_rate_equation` takes its `FittingProblem` `@nospecialize` and hands
`Optimization.solve` a `FitObjective`, which holds `x -> loss!(x, fp)` behind a
`FunctionWrapper`, so every mechanism gives the solver one problem type: `solve` and the
solver compile once per optimizer, not once per mechanism.

### Parameters

The parameter interface is five calls, all in `Derivation`:

```julia
parameters(m)                                # the constants a fit estimates
parameter_units(m)                           # their units, for Fitting
dependent_parameter_values(m, params, Keq)   # the constants the law computes
rate = compile_rate_equation(m)
rate(concs, params, Keq)                     # the rate per unit enzyme
```

`parameters(m)` is exactly the names a fit's `params` holds: `Keq` is given, never fitted,
and is the last positional argument wherever it is taken, never a keyword. The law has no
enzyme concentration; a caller who wants an absolute rate multiplies by it.
`parameters(rate)` and `metabolites(rate)` give the compiled law's names as compile-time
constants, which `loss!` needs.

One type family in `src/Derivation/thermodynamic_constraints.jl`, the union `_Parameter`,
describes every constant a fit estimates, and the joint solve lists them in `parameters(m)`
order (`_Solve.info`):

| Type | What it points at |
|:-----|:------------------|
| `_StepParameter(step, slot, inactive)` | Slot `slot` of step `step` of the catalytic graph, in the inactive conformation when `inactive` |
| `_RegulatorParameter(ligand, inactive)` | The constant of a regulator ligand, in the active or inactive conformation |
| `_AllostericParameter()` | The allosteric constant `L` |

A value holds no other fact and means nothing without its mechanism. At slot 1 a
`_StepParameter` is a rapid-equilibrium step's equilibrium constant or a steady-state step's
rate constant in its stored direction; at slot 2, a steady-state step's rate constant in the
reverse direction. The step's flag says which. A K on a steady-state step or a k on a
rapid-equilibrium one constructs but has no name: `_name` errors on it. `_name(p, m)`
renders the name from the step's sides and the kept member of its kinetic class, and
derives the conformation's prefix (`A_`, `I_` or none) from `inactive` and the group's or
ligand's allosteric state.
`parameters` returns these names as Symbols, and the rate law takes them at compile time.

`parameter_units(m)` gives each constant's units, read from its step or ligand. Fitting
centres its search on these units and scales the constants whose units hold the rate, so
it reads no step itself.

### The joint thermodynamic solve

`_solve` (`src/Derivation/thermodynamic_constraints.jl`) builds one integer system per
mechanism over the constants of its steps, from its thermodynamic cycles and the ties
within its kinetic classes (`_joint_system`), and reduces it exactly
(`_reduce_joint_system`); its pivot columns are the dependent constants. The order of the
columns, read from the graph alone, decides which constants are fitted, never the law
(`_priority_key`). `_solve` returns the dependents with the independents and the fitted
names as a `_Solve`, which the engine, `parameters`, code generation and the printer read.
Each kinetic class is named after its kept member (`_kept_members`) through the
`_name(p, m)` chokepoint, so the naming rule lives in
`src/Derivation/thermodynamic_constraints.jl`. A testset in
`test/Derivation/test_thermodynamic_constraints.jl` walks `src/` and fails on a
`Symbol("K…")`, `Symbol("k…")`, `Symbol("V…")` or `Symbol("L…")` literal outside a
parameter-name renderer, and on a renderer outside that file. The engine (`_law`,
`src/Derivation/rate_equation_derivation.jl`) derives each conformation's law on a mask of
the edges it keeps (`steps_running_in_allosteric_conformations(am)`; every edge of a plain
mechanism). Every edge keeps its stored direction, so both conformations share the
mechanism's orientation and names.

### The law, its code, its text and its key

`src/Derivation/rate_equation_derivation.jl` holds sparse polynomials over integer-indexed
variables (`_Poly`), with rational exponents for constants and checked `Int64` coefficients:
under the term cap no coefficient leaves `Int64`, so an overflow is an internal error.
`_law(m, s, mask)` (`src/Derivation/rate_equation_derivation.jl`) derives one conformation's
law, `_Law`, as `v = H·Q/D`. It contracts each rapid-equilibrium segment to a node rooted at
its first form (free E roots its own), sums the spanning trees of the segment graph, builds
the numerator from each steady-state edge's net flux, substitutes the solve, divides out
`H = Π[S] − Π[P]/Keq` exactly and brings `Q` and `D` to lowest terms in the concentrations.
An MWC law (`_AlloLaw`) holds both conformations' laws, their free-E weights over a common
factor and the regulatory sites' polynomials. `rate_equation_term_count` gives the law's
size before any polynomial is built; the term cap, which `_law` checks, and identify's
`eq_complexity_filter` read it. `_Poly`, `_Law` and `_AlloLaw` are mutable structs with
`const` fields, never an immutable wrapper of one reference: Julia 1.12 and 1.13 can hoist
such a wrapper's box out of a loop that merges it with another type, without a write
barrier.

`_rate_body` (`src/Derivation/rate_equation_derivation.jl`) emits the law as straight-line
locals: the dependents in solve-column order, each segment's weight and tree sums, `Q`, `H`
and the MWC terms, all in greedy Horner form with binary `+` and `*`. Each input is read by
name into a positional local, so a `NamedTuple` in any field order works, and reciprocals
and roots (`sqrt`, or `x^(one(x)/q)`) are computed once.

`rate_equation_string` divides the law by the unit rule's normalizer (`_normalizer`; the
derivation tutorial's [Normalization](@ref) states the rule for users), so every denominator
term is dimensionless whenever the law has one. With `coefficients = true` it prints the law
in named coefficients, and with `coefficient_values = false` without the lines that define
them. `identify_rate_equation` writes that short form to the search and LOOCV CSVs and the
full coefficient form to `best_equation.csv`.

`_canonical` (`src/Derivation/rate_equation_dedup_key.jl`) labels a law canonically: each
part to lowest terms, the fitted constants coloured by refinement on the term–constant
incidence and individualized while colours tie, the concentrations by name. The SHA-256 of
its text is the law's key, by which identify fits one mechanism per law and reads that fit
through each mechanism's own names (the canonical order). At 100 labellings the search
stops, and the law takes its mechanism's own key. The L-phantom certificate
(`_phantom_weights`) solves exactly for a rescaling of the other constants that absorbs `L`
in a one-subunit variant whose inactive conformation turns nothing over, and it alone
decides an L phantom (`_is_phantom`). `rate_equation_dedup_key(m)` gives identify a
mechanism's key and the key of the fit it reads, and `_process_batch` says when an L
phantom reads that fit.

Identify asks one mechanism three questions in a row: `parameters`,
`rate_equation_string` and `rate_equation_dedup_key`. Each needs the guard verdict, the
joint solve or the law, so a private memo (`_Memo`, `src/Derivation/Derivation.jl`) keeps
those of the last mechanism asked about, and each is derived once for all three.

### Reading `Base.gc_live_bytes`

`BigInt` arithmetic raises `Base.gc_live_bytes` by bytes that no live object holds, so
the count overstates what a derivation retains. Class 3's model point (`_model_point`,
`src/Derivation/rate_equation_correctness_guard.jl`; an integer null space in
`Rational{BigInt}`) raises it by 101 KiB on the six-step bi-bi ping-pong at two
conformations, and a 6×6 determinant outside the package does the same in `BigInt` and
`Rational{BigInt}`, not in `Rational{Int}`. The joint solve adds nothing: it runs in
`Rational{Int}` on 10,226 enumerated bi-bi and ter-bi mechanisms without overflow. The count
holds no memory and raises no peak: a process whose count reaches 2.8 GiB peaks at 1.2 GB
resident. Measure retained memory by the resident set.

## Well-posedness guard

The guard (`src/Derivation/rate_equation_correctness_guard.jl`, with
`allosteric_rate_equation_correctness_guard.jl` beside it) decides whether a mechanism's
rate law is the unique fast-equilibrium limit of mass-action kinetics at interior
concentrations and on every **face**, each set of zero concentrations where data sit
(`_faces`). `check_valid_rate_equation` states its conditions. `_guard(m)` reads the
reaction graph through masks and never constructs a mechanism. It returns an isbits
`_Verdict`, never a `Union` that a loop reads: Julia 1.12 and 1.13 can hoist the box of a
`Union` value out of a loop without a write barrier, which corrupts the heap.

The guard reads three networks (`_Network`). Two are of one catalytic subunit whose
conformations flip at free E. The **law network** (`_law_network`) holds the steps each
conformation runs. The **ε-network** (`_eps_network`) holds every step of both
conformations, each excluded step switched off as a **reading** says (`_Reading`). The third
holds every step once, each fast. The guard first checks on the mechanism that no two
constants share a name (`_shared_constant_names`). Then, before any class, it checks on the
third network that every cycle balances each metabolite (`_unbalanced_cycle`). Form names
are opaque, so a step that leaves out a metabolite, writes one on the wrong side or names
the wrong form breaks only a cycle, and the joint solve has no law for that cycle; the
message states the cycle's net change. The checks that follow are classes 1 to 3, which
`_Fault` names, and **admission**. Admission is the guard's check that an allosteric law is
the limit as the excluded steps are switched off, in the interior and on each face.
`_guard_with_reading` and `_verdict_under` give the order of all the checks, and
`_check_face` the steps of class 3.

**Readings.** The steps a conformation keeps join the forms into **islands**
(`exclusion_components`), which are neither kinetic groups nor rapid-equilibrium segments.
Every island but free E's is cut off: it lies higher in energy and holds no enzyme of that
conformation. A reading switches off each excluded step (`_SwitchOff`): a step between two
islands loses the rate of its uphill direction, or at steady state may stop, and a
steady-state step inside an island stops. The guard accepts only a consistent reading, one
under which this picture holds; `_interior_verdict` states the conditions and the reason for
each. A rapid-equilibrium step inside an island cannot be switched off
(`_FAST_STEP_PINNED`): its two sides lie at one energy, so thermodynamics fixes its
equilibrium constant, and a step keeps its kind in both conformations. The law is the same
under every reading, since it drops the excluded steps. Each of the guard's faults but three
holds alike under every consistent reading (`_FAULTS_A_READING_CHANGES` gives the case
analysis). The three are faults of the face condition, which read which transitions vanish,
and only they tell consistent readings apart.

**The search.** `_guard_with_reading` admits a mechanism under the first consistent reading
that passes every check; its docstring gives the order in which it tries readings and the
order of the faults a rejection reports. It tries the three built readings first
(`_built_reading`, `_SlowUptakeReading`): they are cheap and decide most mechanisms. Since
only three of the face condition's faults tell consistent readings apart, any other fault
under the first consistent reading ends the search. Then it searches the consistent
readings (`_consistent_columns`) against admission's face condition alone, since class 3
reads no reading and has passed once on the law network. Above a fixed number of readings
(`_MAX_CHOICES_PER_CONFORMATION`, `_MAX_READING_PAIRS`) it searches none rather than some,
so whether it admits never depends on the order of the steps. A rejection reports the
verdict of the first consistent reading tried, so its message never says that a step
cannot be switched off when another reading switches it off.

**The face condition.** Admission's face condition (`_face_condition`) reads which
transitions of the ε-network are alive on a face and which vanish; its docstring lists the
four faults it refuses, in order, and the condition behind each. Two parts of admission on a
face rest on arguments, not proofs. The filling check is one of the refusals of
`_face_condition`. The set-aside rule is applied by `_admission_face` to a face that
`_face_condition` refuses as a bridge.

The filling check (`_FILLED`) looks at a fast transition alive only in its vanishing
direction that the law's forms feed through fast transitions. It refuses one whose target
has no fast way back to the law's forms that does not vanish. A form's weight is its inflow
over its outflow. Past a slow transition the inflow is bounded by a rate that vanishes;
through fast transitions alone it grows with their speed, and then only an outflow that
grows as fast keeps the form empty in every order of the two limits. A step that loses a
direction can keep its uptake, which dies on the face where its metabolite is zero and
leaves the vanishing direction to fill the form one way. A built reading never does this:
the direction it keeps is a release, alive on every face. A steady-state transition that
fills one way is the case of a stranded basin (`_PUMPED`): fed forms outside the law network
that reach it along no alive transitions that do not vanish.

The set-aside rule (`_without_set_aside_groups` says which groups it sets aside) rests on
this argument: such a group cannot change the limit. The law's forms it touches lie in one
set that the law's own fast steps hold in equilibrium, and a second route between two forms
already in equilibrium changes nothing; with one touching form there is no route between
two forms. The group's own weight vanishes: every form drains at a rate that does not
shrink as the fast steps speed up, none fills faster than it empties, and class 2 has
refused any fast cycle in it that turns the reaction over. A slow way out of a form that
the fast steps carry elsewhere does not count as draining, since its share vanishes as
they speed up.

Neither argument proves that the face condition with these checks suffices under every
reading; the evidence is numerical. A gate ran the guard on 78,554 mechanisms in ten
populations: taggings of uni-uni, bi-bi, ping-pong, ter-bi, bi-ter and ter-ter mechanisms,
with and without an inhibitor or an `:OnlyI` group, one expansion round of the search, and
the test suite's allosteric mechanisms. The guard admits 1,575 of them only through the
search or the set-aside rule. The ε-path and RE-speed oracles follow the two-conformation
mass-action steady state as the rapid-equilibrium steps speed up and the excluded steps
vanish. They ran on the 591 of those with 8 steps or fewer and on samples of the rest. They
agree with the guard on all 591 along corner directions, on 598 larger ones along standard
directions, and on 100 ter-ter ones along corner directions. At two catalytic subunits the
steady state equals the law for 40 of 40 sampled. With the filling check off, the guard
admits 96 more mechanisms of 8 steps or fewer among the `:OnlyA` taggings of the uni-uni,
bi-bi and ping-pong starting mechanisms, and 18 of them fail along corner directions.
`test/Derivation/test_allosteric_rate_equation_correctness_guard.jl` keeps a sampled
population test with corner directions ("sampled taggings: the energy-level rule agrees with
the oracles").

**Faults.** `_Fault` lists every fault a verdict can carry, by the check that finds it.
`_admission_message` writes the text of admission's faults, with each fix whose variant
passes the guard.

On a mechanism the guard rejects, the macro, the parser and the derivation entries throw;
`init_mechanisms`, `seed_mechanisms` and both methods of `Enumeration.candidate_mechanisms`
drop it. The flip move runs the guard on each candidate it builds, the split move on each
candidate whose count rises, and the tagging and relaxing moves on each set of kinetic
groups tagged `:NonequalAI` whose count rises, the tagging move also on that set tagged
`:OnlyA`. The constructors and the Sig decoder never run it, so code that builds a
mechanism inside the package never trips a check meant for the mechanisms it returns.

## Enumeration engine architecture

The enumeration works on candidates (`CandidateMechanism`,
`src/Enumeration/candidate_mechanism_type.jl`): a mechanism (see
[The mechanism data model](@ref)) and the enumerator's record of each form (`_FormRecord`).
The enumerator owns the naming rule (`_form_name`): it names each form from its record
(`EAB`, `ESSmod`, `E_res_plusA_minusP`). Its moves read a form's contents from the record
alone. They tell a plain binding from a fused one, straight into the product-bound form, by
the records of the step's two forms (`_is_plain_binding`). They tell a step that binds a
modifier copy from one that binds the reactant by the copies the two records hold
(`_copy_change`). The seeds write their topologies, dead ends and bases as records and
steps; the moves copy a parent's records and steps with new flags, groups or bound
inhibitors (`_regrouped`, `_with_flags`, `_dead_end_graph`). One function, `_candidate`,
names the forms and builds every candidate, the only place in `Enumeration` that calls a
mechanism constructor. The moves that tag a mechanism or change its regulatory sites keep
the parent's catalytic scheme, so the constructor runs only the state and site checks.

The enumerator refuses a reaction whose names spell one string two ways, such as `NAD`, `P`
and `NADP` (`_check_decodable_names`, run by `candidate_mechanisms`). Its exits and its
flip, split, tagging and relaxing searches drop a mechanism that `_enumerator_keeps`
refuses. It refuses one the guard rejects (see [Well-posedness guard](@ref)). It also
refuses an allosteric mechanism in which a conformation keeps a form, other than free E,
that only bindings it excludes enter (`_kept_forms_keep_a_binding`). The guard admits the
second when its law is the limit, so a user may write it, but the enumerator does not
propose it.

Mechanisms and candidates are ordinary value types to avoid excessive precompilation
costs. The enumeration builds, expands, and deduplicates many thousands of candidate
mechanisms (see [The enumeration engine](@ref)). If each candidate were a singleton
parametric type, merely constructing it would trigger compiler specialization — the same
precompilation cost the derivation pays — multiplied across thousands of mechanisms, which
would be prohibitive. Carrying mechanism structure as runtime values instead means
enumeration and deduplication cost no compilation at all. Only the candidates the search
actually fits have their rate laws compiled, one at a time, through
`compile_rate_equation`.

The moves never emit a child that is provably a reparameterization of its parent. Each
proof sits in the docstring of the function that applies it; this list says which:

- The split move (`_expand_split_kinetic_group`) accepts a set of splits only if the
  fitted-parameter count rises and the guard accepts the built child. The count comes from
  the joint solve (`fitted_parameter_count`), which derives no rate polynomial. A split
  never changes admission's interior verdict, since both parts keep the group's tag.
- The relaxing move (`_expand_change_allo_state`) tags a set of kinetic groups
  `:NonequalAI` only if the fitted-parameter count rises and the enumerator would return
  the child, and a regulatory ligand only if the count rises.
- The flip move (`_expand_re_to_ss`) flips whole groups and accepts a set only if the
  rapid-equilibrium segment count rises and every flipped group carries net flux in the
  child. It never flips a flank of a steady-state chain (`_chain_flank_groups`; see
  [Apparent constants on a steady-state chain](@ref)). In a plain mechanism that is a proof;
  in an allosteric one it is a modeling choice.
- Flux is decided on the graph of rapid-equilibrium segments: a steady-state step carries
  flux exactly when its biconnected block holds a cycle that runs net reaction
  (`_flux_carrying_edges`). A steady-state group with no such step exposes only the ratio
  of its constants, so a flip set that leaves one counts as failed and is extended, and a
  split part with none is emitted at rapid equilibrium.
- A dead-end copy of a substrate or product is redundant when each of its complexes has a
  productive twin and the dwell gauge moves its constant along a direction the rate law
  cannot see (`_redundant_copy_groups`, `_productive_twin`, `_gauge_rescaling`). The
  dead-end and tagging moves drop a child that holds a redundant copy group or a group
  that binds a competitive inhibitor in neither conformation (`_inhibitors_bind`), the
  split drops every child that holds a redundant copy group, and the searches of the
  tagging and relaxing moves count a set whose child holds one as failed.

`candidate_mechanisms` asserts both rules, flux and no redundant copy, on every parent
(`_assert_emission_rules`). The flip, split, tagging and relaxing moves share one
minimal-set search (`_minimal_gaining_sets`). Duplicate laws that survive these proofs
share one fit, an L phantom's key aside: identify fits one mechanism per key
(`_canonical`, `src/Derivation/rate_equation_dedup_key.jl`). A numerical identifiability
rank exists only in the test suite, as an oracle for the proofs and the L-phantom
certificate; nothing in `src/` estimates identifiability numerically.

`init_mechanisms` appends to the seeds their merged and Theorell–Chance variants
(`_seed_variants`, plain `Mechanism`s only), whose rules [The enumeration engine](@ref)
states. A candidate's tests run cheapest first. `_seed_candidate_screen` runs three of
them on index arrays, with a weighted union-find (`has_nonzero_weight_cycle`), and the flux
test also reads only the flags. Only a candidate that passes those four is built and goes
through the guard (`_minimal_flips`). `_degenerate`, the failure of condition V or C on the
active-state graph, marks the starting mechanisms that `candidate_mechanisms`
(`_base_tier`) replaces by their flip children that are not degenerate, before
`seed_mechanisms` binds any required regulator. Both methods of
`candidate_mechanisms` return a `CandidateMechanism` for each mechanism within the two size
limits, which they apply to what they return and never inside a move's own search.

Over more than one catalytic subunit a conformational equilibrium raises every
catalytic-site binding to a power. So a mechanism above one subunit must keep a hyperbolic
catalytic scheme (`_requires_hyperbolic_catalysis`), one whose King–Altman denominator has
degree at most 1 in every substrate and product; `_hyperbolic_catalysis` tests that.
`_expand_to_allosteric` promotes a parent that fails the test at oligomeric state 1 only.
`_enumerator_returns` applies the rule where the enumerator returns its mechanisms,
`candidate_mechanisms(reaction, parent)` and `seed_mechanisms`, whichever move built them;
the searches of the tagging and relaxing moves for `:NonequalAI` sets read it too, since
tagging one more group can make a conformation's law hyperbolic again. The flip search does
not: more steady-state steps never restore a hyperbolic law. The test reads the degree
from the rapid-equilibrium segment graph without deriving the law, and the
`_hyperbolic_catalysis` docstring proves that reading exact. It decides "degree at
most 1" with reachability checks (`_all_reach`) rather than a max-arborescence solve. The
testset "_hyperbolic_catalysis matches the derived denominator"
(`test/Enumeration/test_candidate_mechanism_type.jl`) checks it against the derived
denominator over the uni-bi enumeration and the ping-pong-capable bi-bi enumeration, which
contains every bi-bi mechanism. A future conformational type (KNF, mnemonic, slow
isomerization) adds one `_requires_hyperbolic_catalysis` method and calls
`_hyperbolic_catalysis` in its promotion move.

`test/Enumeration/test_starting_mechanisms.jl` bounds the memory the enumeration holds
and allocates: the bytes each candidate of the bi-bi starting set holds
(`Base.summarysize`), the bytes the bi-bi starting set and one expansion allocate, and the
bytes the ter-ter `init_mechanisms` call allocates, each for a known number of candidates.
The bytes held are the same in every process, on Julia 1.11.9 and 1.13.4 and on every CI
platform, so their bound sits only 3.5% above them: one derived `Int` vector stored on
`Mechanism` trips it, as would any field derived from the steps on `Mechanism`, `Step` or a
candidate's record. The allocation bounds sit 25% above the larger of the values on the two
Julia versions, which differ by about 10%, and the expansion's 3.5% above its value on
1.13.4; more work per candidate trips them. Time and peak resident memory depend on the
machine and its collector, so the ter-ter testset logs both and bounds only the time
outside garbage collection. On the macOS CI runner the ter-ter call spent 13.9 to 116.1 s
in garbage collection in three runs of nearly the same code; `ROADMAP.md` lists these runs
under the slow macOS tests, whose cause is not known.

## Optimization algorithm architecture

Fitting reaches solvers only through Optimization.jl; the package ships no solver of
its own. `fit_rate_equation` (`src/Fitting/fit_rate_equation.jl`) builds an
`OptimizationProblem` whose objective, `_fit_loss`, scores each point with `loss!`
through a `FitObjective` (see [Derivation architecture](@ref)), and calls
`Optimization.solve` with whatever optimizer the caller passes. This gives the
package access to the global, gradient-free optimizers that non-convex
rate-equation fitting needs — CMA-ES (`OptimizationCMAEvolutionStrategy`) and
BBO differential evolution (`OptimizationBBO`) are the tested choices — without
taking on a solver dependency or locking users into one algorithm. Optimization.jl
common options (`maxtime`, `maxiters`, …) are named keyword arguments, and
solver-specific options pass through a `solver_kwargs` named tuple; see
[Loss & optimizers](@ref) for the user-facing view.

### What a fit costs

A fit costs `n_restarts` × evaluations per restart × table rows × time per row, plus
compiling the law once and the solver's own work. The tables below come from 118 laws
sampled from the 2026-11-02 searches, 59 from each, each fitted once, cold, with 20 restarts
on one aarch64 core on Julia 1.13.3. A benchmark copy of `_fit_loss` counted the
evaluations, and `loss!` was timed alone on parameter sets from that fit:

| Median (range) | PGK, 386 rows | PFKP, 486 rows |
|:---------------|--------------:|---------------:|
| Evaluations per restart | 2,830 (1,160–22,900) | 6,230 (2,820–30,100) |
| Evaluations per fit | 56,700 (23,300–458,000) | 125,000 (56,500–601,000) |
| `loss!` per call | 3.3 µs (2.8–18.3) | 6.3 µs (5.1–8.9) |
| `loss!` per row | 8.5 ns (7.2–47) | 13 ns (10–18) |
| Warm time per restart | 0.014 s (0.004–0.27) | 0.050 s (0.020–0.34) |
| Whole fit | 0.38 s (0.18–5.6) | 1.12 s (0.51–7.0) |
| Warm fit ÷ (evaluations × `loss!`) | 1.28 (1.08–1.54) | 1.25 (1.12–1.43) |
| Compiling the law | 0.11 s (0.07–0.22) | 0.12 s (0.12–0.19) |
| Compile share of a fit | 28% (3–54) | 11% (3–23) |

Over all the sampled time, which large laws dominate, a fit divides as follows:

| Search | Mean fit | `loss!` | Compiling the law | Solver and the rest |
|:-------|---------:|--------:|------------------:|--------------------:|
| PGK | 1.05 s | 75% | 11.5% | 13.5% |
| PFKP | 1.41 s | 71% | 9% | 20% |

`loss!` costs about the same inside a fit as alone, so a warm fit is close to evaluations ×
`loss!` time plus 1.5–1.8 µs per evaluation. A profile of warm fits puts the objective at
77–92% of their time, CMA-ES's eigendecomposition (`MEigen`, in CMAEvolutionStrategy) at
5–14% and garbage collection at 1–5%. The objective's wrapper costs too little to measure,
and the package's own work around the solver stays near 1%. Compiling takes 0.10–0.19 s per
law whatever its size, so it reaches 30–50% of a fit on PGK's laws of 6–8 parameters and
under 5% on the largest laws. `test/Fitting/test_fit_rate_equation.jl` bounds a fit's time
by its counted `loss!` calls.

The tables and the profile above measured a `loss!` that called the law once per row, inside
one loop. The present `loss!` runs three loops over the rows. The first writes each row's
prediction to `log_ratios_buffer` with the law inlined; it holds no branch and no call, so
on aarch64 the compiler runs it over two rows at a time in vector registers. The second
turns each prediction into its log-ratio or the penalty, and the third sums the squares. The
loss is the same to the bit: 194,880 comparisons with `===` found no difference, over the 48
gated laws, seven tables of 1 to 500 rows, both modes and 290 parameter points each, among
them points far outside the box where predictions are zero, infinite or NaN. Per row of a
500-row table whose rate signs match the law, on one aarch64 core on Julia 1.13.4, the three
loops take this share of the one loop's time (median over the gated laws of each class,
three alternating processes per version):

| Laws (denominator terms) | Count | One loop, ns/row | Three loops, ns/row | Share |
|:-------------------------|------:|-----------------:|--------------------:|------:|
| Plain, up to 16 | 22 | 4.1–8.4 | 3.8–4.8 | 0.72 |
| Allosteric, up to 16 | 7 | 5.5–11.4 | 4.2–6.1 | 0.68 |
| Plain, 25–96 | 9 | 8.3–17.1 | 4.5–8.0 | 0.55 |
| Allosteric, 18–22 | 3 | 7.7–10.5 | 4.7–5.7 | 0.61 |
| Allosteric, 284 (LDH) | 1 | 33.4 | 17.3 | 0.52 |
| Plain, 396–858 | 6 | 46–96 | 20–37 | 0.44 |

Warm fits with the fit-cost test's setup (500 rows, 5 restarts, CMA-ES; medians of 2
processes × 3 rounds per version) take 0.60–0.67 of the one-loop time on the test's two
laws, PK, PFK-1, HK and Segel's ordered bi-bi, 0.55–0.56 at 96 and 284 terms, 0.48 at 396,
0.39 at 858 and 0.88 on Uni-Uni, with the same evaluations and the same fitted constants.
With `loss!` that fast, the solver's own cost per evaluation weighs more: on the fit-cost
test's two laws a warm fit takes 1.37–1.40 times evaluations × `loss!`, against 1.21–1.22
with one loop. Inlining costs compile time: a fresh mechanism's first `loss!`, law call and
short fit take 4–10% longer at 12–16 terms, 22% at 96 terms and 77–117% at 396–858 terms
(medians of 2 processes). Fits on x86 are not measured; the guard below has read `loss!`
against the law there. The test jobs of CI and Downgrade run with coverage counters
(`julia-runtest`'s default), and under them the loop ran one row at a time on the three
laws checked locally.
They also run with `--check-bounds=yes`, `julia-runtest`'s default too, under which no
`@inbounds` holds; on Julia 1.13 the loop stays scalar with only its inner `@inbounds`
removed.

On the cluster a fit costs more than on one core. The PFKP and PGK searches of 2026-11-02
ran 119,453 and 15,700 fits in 62 and 13.3 minutes, at 25 and 40 core-seconds per fit.
Core-seconds are wall time × workers ÷ fits, so they include startup, enumeration,
cross-validation and idle workers. In `sacct`, the average worker of those runs used the
CPU for 42:47 of 61:18 on PFKP (70%) and 3:16 of 13:01 on PGK (25%), which comes to 17 and
10 CPU-seconds per fit. Workers connected within 15 s, and the master's CPU accounts for at
most 8.5 of PFKP's 18.5 idle minutes and 2.1 of PGK's 9.8. The rest is not measured,
because `progress.log` carries no timestamps. How those 17 and 10 CPU-seconds divide
between the fit itself, slower or contended cores, enumeration, cross-validation and
message passing is also not measured. A fit alone takes 1.41 s on PFKP's laws and 1.05 s
on PGK's on average on one core (the second table above).

Two choices in the generated code (`_rate_body`) carry much of the law's speed. They were
measured in `loss!` per row on the PFKP table (486 rows) and the PGK table (386 rows), on
one aarch64 core on Julia 1.13.3. The PFKP figures are one allosteric law each, of 13 and 20
parameters. The PGK figures are means over 60 laws drawn from those whose loss no other fit
shares, 6 allosteric and 54 plain. Allosteric laws made up 55% of the search's fits, but an
allosteric law often repeats the loss of a plain one.

- The Horner form (`_horner`) removes multiplications. Against the expanded sums it takes
  PGK's allosteric laws from 88 to 37 ns per row and its plain laws from 90 to 29. The
  PFKP laws gain nothing: 39 ns per row either way, and 77 against 76. The code emits no
  `muladd`, and Julia does not contract `x * y + z` by itself, so a law executes no fused
  multiply-adds.
- An MWC law writes each integer power out as multiplications. This was measured later,
  after other changes had taken the two PFKP laws to 27 and 64 ns per row. Julia's general
  integer power (`pow_body`) costs about 2 ns a call, and those laws would make 7 and 20
  such calls per evaluation. Written out, they take 11 ns per row in place of 27, and 20 in
  place of 64.

None of the figures of the two bullets above or of the next paragraph is today's: the first
table above gives the one-loop `loss!` a median of 8.5 ns per row on PGK's sampled laws and
13 ns on PFKP's, and the present `loss!` takes a median 0.62 (0.37–0.98) of the one-loop
time per row over the 48 gated laws, on aarch64 and on 500-row tables whose rate signs agree
with the law.

The law's own time cannot show a slow `loss!`, and `@inferred loss!` passes on a slow path,
since the result is a `Float64` either way. `loss!` reads `parameters(rate)` and
`metabolites(rate)` as compile-time constants (see [Parameters](@ref)): their generated
methods return tuple literals. It builds `NamedTuple` types from those names and their
lengths, one for the parameters and one for each row's concentrations. If either name list
were a tuple whose length the compiler cannot infer, those types would be built, and the law
called, by run-time dispatch, so a plain mechanism would pay it too.
`_testhelper_test_performance` checks both lists with `@inferred`. The cost was measured
when the `metabolites` method of an allosteric mechanism splatted a `Vector` into its
result, and the compiler inferred `Tuple{Symbol, Symbol, Symbol, Symbol, Vararg{Symbol}}`
for the measured law's mechanism. Each row then dispatched `ntuple`, the law, `sign` and
`log` at run time, and plain mechanisms did not pay it. On the 13-parameter PFKP law `loss!`
took about 1,000 ns per row and allocated 778 KB per call, against about 41 ns per call for
the law. In the PFKP search measured with that dispatch, every fit was allosteric. In the
PGK search, allosteric fits were 55% of the fits and, in a random sample measured with that
dispatch, 97% of the fit time.

So the speed gate in `test/Derivation/test_rate_equation_derivation.jl` times `loss!` as
well as the compiled law. It runs on every gated law: each mechanism of
`MECHANISM_TEST_SPECS`, the three largest enumerated bi-bi laws (396 denominator terms each)
and three ter-bi laws (589, 799 and 858 terms). `.claude/CLAUDE.md` ("Guards") states the
bounds. `_testhelper_test_rate_equation_performance` times the compiled law, and
`_testhelper_loss_time_per_row` times `loss!`, both in `test/testhelpers.jl`.

The test process does not time them. CI runs the suite with coverage counters on every line
of the package and with `--check-bounds=yes`, under which no `@inbounds` holds, so in the
test process the gate would time code no user runs: under coverage the first loop of `loss!`
runs one row at a time, and the placement effects on the law read about twice as large
(`ROADMAP.md`, "Where their data fall in memory"). The first test that asks for a time
starts one child process from the bare Julia executable, without those flags
(`_testhelper_speed_readings`). The child includes `test/testhelpers.jl`, times every gated
law and prints the readings, and the tests compare them with their bounds. The allocation
checks stay in the test process, so CI runs them under coverage and bounds checks. Locally
(aarch64) the child takes about 15 s: 0.6 s to load the package, 2.6 s to include
`testhelpers.jl` and 12 s to time the 48 laws and read `loss!` for the guard below. Each
margin has a reason:

- The compiled law's bound has a large intercept because small laws take 1-8 ns on an
  aarch64 development machine and 5-38 ns on CI under coverage counters, varying by runner.
  It catches only a large slowdown of a small law; its part per denominator term binds on
  large laws. Where the law's frame and its inputs fall in a 4 KB page can slow every call
  of a process. So each of the five batches has its own stack depth and its own copy of the
  inputs. In ten child processes on aarch64 no law took more than 0.22 of its bound.
- `loss!` is timed against the compiled law and `log` in one process, so machine speed
  cancels out. The bound holds if any batch meets it. Each batch times `loss!` on a new
  `FittingProblem`, since where a problem's arrays fall on the heap can slow `loss!` up to
  about 2× on aarch64, and one stack frame deeper than the last, since one scan of stack
  depths slowed HK 1.50–1.56× at three adjacent depths. x86 and 16 KB pages are not
  measured. `_testhelper_loss_time_per_row` gives the reasons for the factor on the law's
  time and for the `log` calls. In ten child processes on aarch64 `loss!` took 0.11–0.30 of
  the time it is allowed.

Both bounds were set from timings taken under coverage, and stay until CI's runners have
been read without it (`ROADMAP.md`, Planned).

Neither bound sees `loss!` lose its vectorized first loop. With the loop `loss!` takes
0.40–0.46 of the 858-term law's time per row, and with the law called once per row about
1.1, while the gate allows 1.25 times the law and more. So the testset "loss! runs its
predictions in a vectorized loop with the law inlined"
(`test/Fitting/test_fit_rate_equation.jl`) reads two more things from the child:

- `loss!`'s optimized IR (`code_llvm`) on Segel's ordered bi-bi, PFK-1 and the random-order
  bi-bi (`_testhelper_loss_ir_mechanisms`) must hold a division of a vector of `Float64`
  inside a block that LLVM's loop vectorizer labels `vector.body`, and no call to the
  compiled law. The prediction loop is the only loop of `loss!` that divides. The SLP
  vectorizer also puts vector divisions in straight-line code, so only those in a
  `vector.body` block count; the sum of squares in absolute mode is a vectorized loop with
  no division. A loop with a call never vectorizes, so the count of calls only names the
  cause. Uni-uni would not serve: Julia inlines its law without `@inline`, and its loop
  vectorizes without the inner `@inbounds`.
- On the 858-term law, `loss!`'s time per row must be at most 0.7 of the law's time per
  call, each the best over 20 batches at 20 stack depths with the law timed at each batch's
  depth (`_testhelper_loss_over_law`). The gate's statistic, the best batch's ratio with the
  law at one depth, let slower variants pass a bound of 0.75 times the law plus two `log`
  calls.

Run on copies of the package on aarch64, the testset fails each of four variants of `loss!`.
With the law called once per row in one loop, or without `@inline` on the law, all three
tests fail (1.03–1.10 of the law). Without the inner `@inbounds`, or with the sign test and
`log` back in the first loop, `loss!` reads 0.77–0.79, only 10–13% over the bound, and the
IR check is the reliable catch. The testset passes the present code and an `ifelse` in the
second loop, and the spike that designed it passed `@simd ivdep` and the loop without
`@simd`, which changes neither its vector body nor its time, on Julia 1.13 and 1.11. The
present code's ratio read at most 0.412 in the spike's 36 processes, 0.59 of the bound, and
0.397–0.414 in ten child processes.

On CI's runners, read in one run, the loop vectorized on all three laws in every job, and
the ratio was 0.16 on Ubuntu and 0.19 on Windows (x86-64), 0.24 on Ubuntu with Julia 1.11
and 0.42 on macOS (aarch64).

What it cannot see:

- the speed of any other law's `loss!`: only the 858-term law is timed, and the three laws
  read in IR are checked for vectorization only;
- removing `@simd`, and on Julia 1.11 removing the inner `@inbounds`, neither of which costs
  anything there;
- the second loop and the sums, about a tenth of `loss!` on the 858-term law, so doubling
  them would move the ratio from about 0.41 to 0.45 (inferred);
- the heap-placement effect, which the best of 20 batches ignores on purpose;
- a later loop that divides and vectorizes, such as an inlined vectorized `log`, which
  would pass the IR check as the prediction loop does (inferred).

`test/Fitting/test_fit_rate_equation.jl` checks `loss!` on each column type a table can
bring and on a cross-validation fold. `RateData` stores `Rate` and each metabolite as a
`Vector{Float64}`, because an abstract column type costs a run-time dispatch on each read.

`test/test_type_stability.jl` requires JET's optimization analysis to report nothing on
every exported and public function and on the workers' internal functions, on Julia 1.13
only, since inference changes between versions. Each row reads as a call and says whether it
must be free of dispatch (`@_testhelper_no_dispatch`) or may dispatch only for named reasons
(`@_testhelper_dispatches_in`), and a row that stops needing a reason fails. A public
function with no row fails the test; a public type does not, so deleting the row of
`FittingProblem`, `RateData` or `Mechanism` would pass. Five limits apply:

- JET does not look past a run-time dispatch or into a `pmap` closure, so the bodies the
  workers run have rows of their own.
- JET gives a report in a shared callee once per row, so an exempt function that reaches it
  first hides it from other functions.
- The test drops a report below a call that cannot compile, which mostly describes code that
  never runs at those types, and a report through `Base.slowlock`, which is not package
  code. A `@nospecialize` body is therefore checked only at the concrete types its rows
  name, and at the abstract signature that actually compiles, which has rows of its own
  (`fit_rate_equation` on an abstract `FittingProblem`, `_scale_rate_constants` on a
  mechanism of type `Any`).
- JET misses a dynamic call written directly in the body of a comprehension or generator
  whose result is untyped, such as `[(r[] + x)::Float64 for x in v]` with `r::Ref{Any}`. It
  does report the same call with a typed result (`Float64[...]`), in an explicit loop, or in
  a callee the comprehension calls.
- JET misses a single dynamic call in the body of a `@generated` function, so the allocation
  and time bounds above stay the hard guard on the compiled law.
