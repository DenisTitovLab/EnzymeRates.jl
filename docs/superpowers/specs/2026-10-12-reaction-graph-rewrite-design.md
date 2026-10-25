# Reaction-graph rewrite: design

Status: approved by Denis on 2026-10-12, with three reviews folded in afterwards; section 18
lists the author's proposals that implementation follows unless Denis objects. This spec is
the written record of the decisions made between 2026-10-09 and 2026-10-12: decision cards
A1–E11, the brainstorming sessions and the chat-decision ledger. Appendix A maps each card
to the section that carries it. Implementation follows one plan per phase (section 15), each
written with the writing-plans workflow and executed test-first with subagent-driven
development.

## 1. Summary

EnzymeRates derives a rate law for each enzyme mechanism, fits it to data and searches
mechanism space for the law that fits best. The derivation grew around per-step symbols,
a second rename pass and a separate inactive-state mechanism, which made results depend on
how a mechanism was written (fitted names and dedup hashes changed under regrouping) and
left known bugs (phantom names, fold hazards, ill-defined `:OnlyA` cases, a saturating kcat
that is 0 or infinite for valid mechanisms).

The rewrite derives every law from the reaction graph: one thermodynamic solve over the
cycle space, both allosteric conformations from one graph, `:OnlyA` and `:OnlyI` tags as
limits along thermodynamically consistent paths, and one well-posedness guard. It replaces
the mechanism types, the DSL, the parameter API, the printed laws, the dedup key and the
fitting scale. A vertical-slice prototype of most of this design derived 237,965 laws from
257,698 mechanisms, and every law matched an exact oracle; Appendix B lists the parts it did
not prototype.

## 2. Goals, scope and success criteria

### Goals

- **Correct laws.** Every derived law equals the steady state of mass-action kinetics in the
  limit the mechanism states (fast rapid-equilibrium steps, vanishing `:OnlyA` or `:OnlyI`
  steps), and satisfies detailed balance.
- **Deterministic results.** A mechanism's law, fitted names, parameter order and dedup key
  depend only on the mechanism, never on step order, step direction, group labels or
  process.
- **Readable output.** Printed laws use independent parameters and a Cleland-style
  coefficient form, and their denominator terms are dimensionless whenever the unit rule
  finds a normalizer (section 9). `show` prints the DSL.
- **One rule per question.** One well-posedness guard decides which mechanisms the package
  accepts; admission, one of its checks, decides which allosteric tag sets it accepts; one
  scale rule reports fitted rate constants.
- **Fast evaluation.** `rate_equation` stays allocation-free and within the speed gates
  (section 11).
- **Fast construction.** The `Mechanism` constructor builds a mechanism in about 2 µs (the
  slice's median was 1.78 µs); a benchmark test gates it from phase 3.

### Out of scope

- Enumerating iso mechanisms (the DSL accepts them; no move generates them).
- Hill's cycle-flux numerator and block factorization (C1).
- A textbook layer in the package, such as a coefficient map or Cleland constants (C5).
- Generating fewer duplicates at the source (C7); the key and the fit cache handle them.
- The iso dwell-time-gauge phantom and two-conformation kinetic cooperativity (ROADMAP).
- The rapid-equilibrium isomerization gauge of 5 enumerated ping-pong mechanisms (a D1
  finding that B9 counts as part of the declared model; proposal, section 18) and the
  flux-free steady-state step that keeps both rate constants fitted (ROADMAP known issue
  3, allowed by B8); the rank analysis tool reports both.
- Using the identifiable rank in parameter counts or model selection (B9); the one
  enumeration use is the L-phantom flag in section 13.

### Success criteria

Each criterion is a test in the suite.

1. **Exactness.** On populations regenerated in the tests (section 14) and on the
   adversarial fixtures (R2, the AD set, the RE/SS mutants and the trap specs), every law
   matches the exact mass-action oracle, plain and two-conformation; rapid-equilibrium laws
   match the RE-speed oracle as the RE speeds grow, including on zero-concentration faces;
   `:OnlyA` and `:OnlyI` laws match the ε-path oracle. Section 14 states each oracle's
   tolerance.
2. **Detailed balance.** At equilibrium concentrations the net rate is zero, per
   conformation and combined.
3. **Determinism.** Rewriting a mechanism (step order, step direction, group labels)
   changes no law, fitted name or key; a same-model regrouping changes no fitted name or
   key; renaming fitted constants changes no key, while mechanisms that differ only in
   which metabolite binds (I versus J) get different keys; a second process gives
   identical keys and generated code.
4. **No phantom names or fold hazards.** Over the regenerated bi-bi and ping-pong
   allosteric populations, no `:EqualAI` constant takes a second fitted name, and
   `length(fitted_params(m))` equals the declared model's dimension, which the test counts
   on its own: the constants in the law after group ties, less the rank of the cycle rows
   over them, both built by the test from the mechanism (ROADMAP known issues 1 and 2).
5. **Speed and compile cost.** Section 11's gates and the constructor benchmark hold.

## 3. Architecture

```
@enzyme_mechanism / @allosteric_mechanism        (DSL; parser reads show() output back)
        │
        ▼
Mechanism / AllostericMechanism                   (canonical value types; representation checks)
        │
        ▼
well-posedness guard + admission                  (DSL boundary, enumeration, derivation entry)
        │
        ▼
joint thermodynamic solve                         (cycle space; both conformations; dependents)
        │
        ▼
law (engine) + face support check                 (the guard's check on the derived law)
        │   ──► normalization and printing        (rate_equation_string, coefficient form)
        │   ──► dedup key                         (identify's fit cache and beam)
        │   ──► code generation                   (rate_equation, @generated over the Sig)
        ▼
fitting (FittingProblem, fit_rate_equation) ──► scale (peak turnover; kcat rescale after fit)
        │
        ▼
identify (enumeration moves, beam, cross-validation)
```

Planned source layout (proposal, section 18), one concern per file: `types.jl` (value
types, canonical form, representation checks, Sig), `names.jl` (the parameter-naming
chokepoint), `dsl.jl` (macros, parser, `show`), `poly.jl` (exact polynomials), `thermo.jl`
(joint solve), `admission.jl` (guard and admission), `engine.jl` (law), `normalize.jl` and
`print.jl`, `key.jl`, `codegen.jl`, `fitting.jl` (with the scale),
`mechanism_enumeration.jl`, `identify_rate_equation.jl`. The plans may merge small files.

The constructors enforce representation invariants only (proposal, section 18). The guard
and admission run at the DSL boundary (macros and parser), on every enumerated child and at
the derivation entry, never inside internal constructions; the face support test, which
needs the law, runs when the law is derived. Internal code that builds a sub-graph reads
the graph directly instead of constructing a mechanism, so a check meant for user input
cannot fire on it (the class of bug PR #87 fixed).

## 4. Mechanism types (phase 3)

### Mechanism

`Mechanism` is the reaction graph:

- substrate and product names;
- a metabolite table; each entry is a name and whether it binds as a modifier, so a
  reactant bound as a modifier (`P::Mod`) is its own entry;
- conformation labels;
- enzyme forms (conformation label, bound metabolites, residual), in canonical order with
  free E as form 1;
- edges (two form indices, consumed and released metabolites, an RE flag), stored
  canonically;
- the kinetic group of each edge.

`AllostericMechanism` holds a catalytic `Mechanism` plus one tag per kinetic group, the
catalytic multiplicity and the regulatory sites (composition). A future model with its own
derivation, such as KNF, wraps the same graph in its own type.

Invariants:

- **Canonical form.** Two writings of one mechanism produce equal fields. Step direction
  follows today's orientation tiers; ties fall to form index (C9).
- **Content identity.** `==` and `hash` compare content and are stable across processes.
- **Free E.** Every mechanism contains a free form whose conformation label is `E`
  (nothing bound, no residual); it is form 1 (E10). Phase 3 renames the enumeration test
  fixture whose free form is `E_c` (`test/test_mechanism_enumeration.jl`).
- **Limits.** At most 32 metabolite-table entries and 15 copies of one entry on a form or
  step side, with errors that say so (C9).
- **No `EnzymeReaction` field.** Atom checks take the reaction's atom table from the
  caller. `EnzymeReaction` stays as the enumeration input.

### Constructor checks

The constructor rejects:

- missing or overlapping substrates and products; a residual that adds a product or
  subtracts a substrate;
- an undeclared name, or a name bound as a modifier without being declared one;
- a self-loop, a metabolite on both sides of a step, the same step written twice;
- group members that, in their physical-forward direction, take up or give off different
  metabolites or differ in RE/SS flag; or members whose direction the chemistry does not
  fix (orientation tiers 2–3) and that no orientation aligns to one change of conformation
  label. Members whose direction the chemistry fixes align by chemistry, whatever their
  change of conformation label (B7);
- two distinct forms that render one name (`E(Ac, CoA)` and `E(AcCoA)`, or `E(P::Mod)`
  and `E(Pmod)`), since every constant is named after the forms its step joins;
- steps that do not form one connected graph;
- a binding that changes a metabolite's role (`E + P::Mod ⇌ E(P)`);
- a mechanism without free E;
- reserved names (`E_total`, `Keq`, `L`, `N` and the printed local names) and `_` in a
  metabolite name (E5);
- a declared reactant that no step uses;
- limits exceeded.

### Compiled types and the Sig

`EnzymeMechanism{Sig}` and `AllostericEnzymeMechanism{Sig}` carry the mechanism as their
only type parameter: one homogeneous integer tuple that decodes back to the mechanism
(`decode(sig(m)...) == m`). The generator decodes the Sig by splatting and derives the law
(C3), so each mechanism compiles one method instance and equal mechanisms share it.

## 5. DSL, `show` and CSV (phase 3)

### Two macros, one vocabulary

`@enzyme_mechanism` builds a `Mechanism`; `@allosteric_mechanism` builds an
`AllostericMechanism`. Each returns its compiled type (`EnzymeMechanism{Sig}`,
`AllostericEnzymeMechanism{Sig}`), as today, so `rate_equation` on macro output takes the
allocation-free path (proposal, section 18). They share the labels `substrates:`,
`products:`, `modifiers:` and `steps:` and the step syntax. Only `@allosteric_mechanism`
takes `catalytic_multiplicity:`, `allosteric_regulators:`, `regulatory_site(n):` lines and
a tag on every step or group.

```julia
@allosteric_mechanism begin
    substrates: A, B
    products: P, Q
    modifiers: I, P                    # bind the catalytic graph without reacting
    allosteric_regulators: X, Y        # bind regulatory sites
    catalytic_multiplicity: 2
    steps: begin
        E + A ⇌ E(A)                              :: EqualAI
        E(A) + B ⇌ E(A, B)                        :: NonequalAI
        E(A, B) <--> E(P, Q)                      :: OnlyA
        E(P, Q) ⇌ E(Q) + P                        :: EqualAI
        E(Q) ⇌ E + Q                              :: EqualAI
        (E + I ⇌ E(I), E(A) + I ⇌ E(A, I))        :: EqualAI
        E(Q) + P::Mod ⇌ E(P::Mod, Q)              :: EqualAI
    end
    regulatory_site(2): X::OnlyA, Y::OnlyI
end
```

- `modifiers:` lists names that bind the catalytic graph without reacting. `::Mod` marks
  only a reactant bound in that role. Constant names: `K_EI_to_E_I` for a modifier-only
  name, `K_EPmodQ_to_EQ_Pmod` for a reactant's modifier copy. Residuals print as
  `_res_plusA_minusP`.
- Every declared allosteric regulator appears in exactly one `regulatory_site(n):` line,
  which carries the site's multiplicity and each ligand's tag. There are no default sites.
- Kept: steps written in either direction, kinetic groups as parenthesized tuples, bare
  step tags, residuals (`E(P; residual = A - P)`), conformation labels such as `Eprime`,
  `Estar` or `F` for forms other than free E.
- Removed: `catalytic_steps:` (now `steps:`), `catalytic_inhibitors:` and the plain
  macro's `regulators:` (now `modifiers:`), `::Inh` (now `::Mod`), the parenthesized
  single-step tag form, and the split allosteric declaration with implicit default sites.

`@enzyme_reaction` (enumeration input) keeps `substrates:` and `products:` with atoms,
`competitive_inhibitors:`, `allosteric_regulators:` (with `::Activator`/`::Inhibitor` and
per-regulator multiplicities such as `Z(1, 2)`), `oligomeric_state:` and
`shared_catalytic_site:`. It drops `dead_end_inhibitors:` (now `competitive_inhibitors:`),
`allowed_catalytic_multiplicities:` (now `oligomeric_state:`, which takes an Int or a
tuple such as `(1, 2, 4)`, so the search still spans several multiplicities; proposal,
section 18) and multi-argument atom brackets (`B[N, P]`). It rejects the same reserved
names and `_` in metabolite names as the mechanism constructors (E5), so a bad name fails
at the reaction, not inside enumeration.

### Display and CSV

- `show(m)` prints exactly the DSL text that builds `m`, with the macro for its type. The
  compiled types print the same text.
- An internal parser reads that text back by walking the parsed expression, never with
  `eval`, and returns a mechanism equal to `m`.
- identify's CSV stores the DSL text, a short stable `mechanism_id` and `parent_id`, the
  parent's `mechanism_id`, in place of today's type strings (about 2,700 characters, up to
  79% of a results file). `mechanism_id` is the first 16 hex digits of the SHA-256 of the
  DSL text, so it is stable across processes and Julia versions; two mechanisms with one id
  are an error. The key (section 10) replaces the `eq_hash` column. The CSV keeps the
  `fit_inherited` column and adds the fit source's `mechanism_id`, an L-phantom flag
  (section 13) and a peak-at-range-edge flag (section 12). Proposed column names:
  `mechanism`, `mechanism_id`, `parent_id`, `key`, `fit_inherited`, `fit_source_id`,
  `L_phantom` and `peak_at_range_edge` (section 18).
- The CSV's `rate_equation` column holds each law's coefficient form
  (`rate_equation_string(m; coefficients = true)`, section 9), since users browse the laws
  there and skim their terms; `rate_equation_string` prints the full law by default.

## 6. Thermodynamics, names and the parameter API (phase 2)

### One joint solve

- **Columns.** One bare constant per step slot: a steady-state step has a forward and a
  reverse rate constant; a rapid-equilibrium step has one equilibrium constant (a binding's
  K is its dissociation constant).
- **Allosteric columns.** The inactive conformation adds no column for an `:EqualAI` group
  (it reads the active column) and one per slot of every member of any other group: the
  members of a `:NonequalAI` group are tied by group rows, and each member of an `:OnlyA`
  group is its own class there (B2). The `:OnlyA` columns lie outside the law and follow
  the limits of section 7. An `:OnlyI` group mirrors this with the conformations
  swapped: its finite constants are its inactive columns, tied by group rows, and its
  active columns, each member its own class, vanish (phase 4, when `:OnlyI` catalytic
  groups become legal).
- **Rows.** The fundamental cycles of the graph, in each conformation: a Haldane row
  (`= log Keq` times the cycle's multiple of the net reaction) or a Wegscheider row
  (`= 0`). Inactive cycle rows cover the whole graph, islands included. Group rows tie the
  members of a kinetic group slot by slot, aligned forward with forward (B7): each further
  member equals the first.
- **Solve.** One exact rational row reduction (checked `Int64`, retried in `BigInt` on
  overflow) chooses the dependents and expresses each as a power product of independents
  and `Keq`. Rational exponents such as `Keq^(1/2)` are accepted (B6). There is no second
  rename pass. A cycle whose net change is no multiple of the net reaction, and rows that
  reduce to 0 = log Keq, are errors, as today.

### Fixed priority

The solve orders columns by these keys in turn; the earlier a column, the sooner it
becomes a dependent:

1. the block: columns outside the law first (inactive copies on pruned forms, vanishing
   `:OnlyA` and `:OnlyI` copies), then the other inactive columns, then the active ones;
2. the step's kind: isomerizations (steps that exchange no metabolite, on-enzyme chemistry
   included), then other chemistry, then plain bindings;
3. steps whose two forms both hold a metabolite before the others (steps that touch free E
   or a form that holds only a residual);
4. steady-state before rapid-equilibrium;
5. the constant that runs against the reaction direction first (E7): the reverse of the
   physical-forward direction, except that for a product binding the uptake runs against.

Ties fall to a content key built from base names: each form's conformation label, bound
(name, role) pairs and residual name sets, and the step's consumed and released (name,
role) pairs; never rendered names. The smaller key is kept (C8). Phase 2 computes this on
today's types, with the same content key and with C9's form-index orientation for the steps
whose direction today's Tier 3 decides by name; phase 3 computes it on the graph type with
identical results, so phase 3 cannot change the dependents. Renaming a metabolite can still
change which constant is dependent (C8; 34 enumerated mechanisms), and the docs say so
(section 16).

### Names

- A rate constant is named after the direction it drives (`k_EA_to_EP`); a binding's K
  after its release (`K_EA_to_E_A`); any other K after its forward direction.
- Each kinetic class is named after its member column that comes last in elimination
  order, which is the one the solve keeps whenever it keeps one (E3), so a name depends on
  the model alone.
- `A_`/`I_` state tags appear except on `:EqualAI` groups and plain mechanisms;
  `:EqualAI` regulator ligands print bare (`K_Rreg`, B10). Fixed names: `L`, `Keq`,
  `E_total`.
- Each vanishing `:OnlyA` or `:OnlyI` copy is named after its own step. These names occur
  only inside the solve; no law contains them, and admission messages name the step by its
  DSL text.
- All `Parameter → Symbol` rendering goes through the naming chokepoint, as today.

### Parameter API

- `parameters(m)`: the independent constants plus `Keq` and `E_total`. Exported, as today.
- `fitted_params(m)`: the independent constants. Internal, as today.
- `all_parameters(m, params)`: `params` merged with every dependent step constant, computed
  from the independent values and `Keq` in `params`. Exported, which makes it new public
  API (proposal, section 18); it replaces the test helper `_testhelper_compute_all_params`.
- Lists keep today's order (E4: alphabetical, which is today's order): the catalytic
  constants sorted by name in Julia's string order (`isless`, byte order, so `K_…` precedes
  `k_…`), then the regulator constants (active, then inactive, by site and ligand), then
  `L`; `parameters` appends `Keq` and `E_total`.
- Removed (B5): the `Full`/`Reduced` modes and exports, their types, the three `Full`
  methods of `parameters`, `rate_equation` and `rate_equation_string`, and the Full-mode
  section of the textbooks page.
- Flux-free steady-state steps stay allowed (B8).

## 7. Allosteric model and admission (phases 2 and 4)

### Model

- The conformations interconvert only at free E (form 1), with `L = [E_I]/[E_A]` (E10),
  from phase 4 on. Each conformation runs its catalytic cycle on the same graph; the rate
  combines them with today's MWC formula (`docs/src/deriving/mwc_allostery.md`), exact at
  one catalytic subunit and the documented commit-per-cycle model above that.
- Catalytic groups take `:EqualAI`, `:NonequalAI` or `:OnlyA`. `:OnlyI` is allowed only on
  groups whose steps bind a modifier (a `modifiers:` name or a reactant's `::Mod` copy),
  never on a step that takes up or releases a substrate or product as a reactant; phase 4
  lifts today's ban on `:OnlyI` catalytic groups for these. Regulatory-site ligands take
  any of the four tags.
- Phase 2 derives both conformations from the one graph with today's semantics (deletion
  and pruning) and today's guard. Phase 4 brings the flip at free E, the limits and the
  admission below; until then PR #86's one-flip-segment check stays.

### Limits (B1, B2)

In the conformation a tag excludes:

- a plain binding loses its affinity: K → ∞, or for a steady-state binding kon → 0 with
  koff fixed;
- any other step (an isomerization, a fused binding, a Theorell–Chance step) is blocked:
  both rate constants go to zero with their ratio set by the cycle rows; a
  rapid-equilibrium member keeps the K the rows set, and its speed vanishes;
- each member of an `:OnlyA` (or `:OnlyI`) group vanishes at its own rate; the rows still
  bind them along the path. Per-member rates change which tag sets admission accepts,
  never a derived law.

Mixed groups (a plain and a fused binding of one metabolite) stay allowed.

### Admission

Admission replaces `_onlya_haldane_violation`; it is one of the guard's checks (section 8).
For each excluding conformation it checks three things exactly (the slice's admission took
about 1 µs); a tagged step here is one whose tag excludes that conformation:

1. **Existence.** A thermodynamically consistent path to the limit exists: a direction of
   vanishing exponents that satisfies every cycle row with the right signs. The graph test
   ships. Take the components of the forms that the conformation's other steps join, and
   make each tagged plain binding an arc from its unbound side's component to its bound
   side's. A path exists exactly when no arc stays inside one component and the arcs close
   no directed cycle. The test is provably equivalent to the exact cone (LP) test: with
   every other constant fixed, a direction's equilibrium-ratio exponents are differences of
   component potentials, and only the arcs need a sign. It agreed with an exact LP on all
   2,995 allosteric records, and the tests keep the LP as an oracle.
2. **Uniqueness.** The limit does not depend on how fast the tagged constants vanish.
   Structurally this holds when free E's component reaches every component along the arcs,
   so every island lies behind a tagged plain binding, on its bound side; the limit is then
   deletion plus pruning. It must also hold on every face where data sit (section 8), which
   rejects cases such as R2 below. Phase 4's face condition (B v3) rejects R2 and is sound
   against the ε-path oracle; conservative on R-B1 and R-B2, listed: it also rejects
   equilibrium bridges (R-B1) and stranded basins fed only through vanishing steps (R-B2)
   whose limit exists, and ROADMAP lists both.
3. **Turnover.** The active conformation still turns the reaction over.

A tag set that fails is rejected: an error naming the step by its DSL text for a
hand-written mechanism, dropped by enumeration. Expected changes against today: AD1 and
AD11 flip to rejected, AD15 flips to admitted, and AD14 and AD14b stay admitted. These
hand-built fixtures live on the reference branch (`reference/proto/skeptic_P5a-onlyA/`:
`skp_build.jl`, `skp_build3.jl`, `skp_build4.jl`), and phase 4 writes each inline in its
tests. R2 (`E + S ⇌ E(S) :: OnlyA`, `E(S) <--> E(P)`, `F + P ⇌ E(P)`,
`F <--> E :: OnlyA`) is rejected: at zero product its inactive limit depends on how fast
the isomerization vanishes.

**Allosteric ping-pong.** Today's tagging move tags every chemistry group `:OnlyA`. On a
ping-pong mechanism that stops both half-reactions in the inactive conformation, and the
covalent forms between them become an island reached through chemistry alone: its
weight, and with it the rate, depends on how fast each half stops. The probes counted 81
such tag sets in their samples: 16 in a bi-bi expansion and 65 among the tagging move's
`:OnlyA` children on a ping-pong-capable bi-bi and a ter-bi reaction (recipe and sampling
seeds in that folder's `skp_export.jl`). Every K-type tag set builds on its V-type one,
and the change-tag move relaxes chemistry only from a parent with no `:OnlyA` binding, so
dropping these would leave identify no allosteric ping-pong at all. Denis requires
allosteric ping-pong (2026-10-13): PDH and OGDH are ping-pong and allosteric. Phase 4
therefore changes the tagging move. It keeps today's V-type child when admission accepts
it; otherwise it emits each smallest set of chemistry groups that stops every catalytic
cycle of the inactive conformation and passes admission. A ping-pong mechanism gets one
V-type child per half-reaction, the other half running in both conformations, as PDH's
phosphorylation stops E1 and leaves E2 and E3 running. K-type children build on each.
A mass-action check at one catalytic subunit, with the flip at free E, ε = 1e-60 and
512-bit arithmetic, backs this on the interior and every single-reactant face: with one
half `:OnlyA`, alone or with that half's substrate binding, the rate equals the deletion
law to 1e-57 for every vanishing speed; with both halves, it moves by up to 100% with the
speeds. Phase 4 adds that check to the oracle tests and confirms every new child against
the ε-path and RE-speed oracles before its plan. Stopping both halves needs the
conformations to flip at the covalent form as well, with an equilibrium constant of its
own; that belongs with iso mechanisms (ROADMAP).

### Thermodynamic consistency

1. **By construction.** Every finite constant, in the law or not, is a column of the joint
   solve and satisfies every cycle row at every parameter point.
2. **By admission.** Every vanishing constant is the limit of a path that satisfies every
   row at every ε.
3. **By test.** The detailed-balance, ε-path and RE-speed oracles (section 14) check both.

## 8. Well-posedness guard (phase 4)

A mechanism is well defined when its rate law is the unique fast-RE limit of mass-action
kinetics wherever data sit. The package fits every well-defined mechanism and rejects the
rest (E2): an error with the cause and the fix for a hand-written mechanism, a drop in
enumeration.

### Rejected

1. **No catalytic cycle.** No cycle turns the reaction over, so v ≡ 0. Structural.
2. **RE turnover.** A cycle of rapid-equilibrium steps turns the reaction over, so the rate
   is infinite. Structural; today's "no finite rate" check, ported.
3. **A law that is not the unique fast-RE limit on a face.** The faces are every product
   zero, every substrate zero and each single reactant zero; each of these, and the
   interior, also with every modifier zero (proposal, section 18). Per rapid-equilibrium
   segment and face, the segment's lowest forms are those whose weight within the segment
   carries, for each of the face's metabolites, the segment's lowest power of it, so that
   weight stays positive on the face. A structural check finds:
   - a bottomless segment: no lowest form;
   - a split: the lowest forms are joined only through forms that are not lowest;
   - a trap: a form that is not lowest, is fed on the face (by a step that still runs
     there, from a form the face populates), and whose fast exits (rapid-equilibrium
     directions that take up none of the face's metabolites) all fail to reach the lowest
     forms. A slow exit does not clear a trap. For a trap the law equals the limit when
     the trap-feeding step is the slowest RE step, so the error names that step and
     suggests making it steady-state.
4. **`:OnlyA` or `:OnlyI` tags with no limit.** Admission (section 7).
5. **A denominator that vanishes on a face**, modifiers at zero included. This support
   test runs on the derived law, when the law is derived, and catches the steady-state
   dead ends the structural check misses.

### Fitted

Substrate inhibition and abortive complexes, flux-free steady-state steps, product
activation, and hand-written laws that never saturate. Identify keeps today's `_degenerate`
screen: it never fits a degenerate mechanism (no maximal rate in one direction, or
chemistry in rapid equilibrium with both sides) but expands it through its flip children.

### Phase 4's first task

Make the face check exact against the RE-speed oracle. Admission's face condition is
sound; conservative on R-B1 and R-B2, listed: it rejects equilibrium bridges (R-B1) and
stranded basins fed only through vanishing steps (R-B2) that the ε-path oracle admits, and
ROADMAP lists both. The probes found false alarms (6 on
mutants; 16 of 16 trap verdicts in the large census, each a trap fed only from a form that
itself holds a zeroed product) and 2 misses (non-unique face steady states). The guard
ships only after the check matches the oracle on populations, mutants and hand-built
cases, R2's active conformation among them as a case with no trap. Then the two specs with
traps ("Numerator: redundant SS-bind", "LDH i-state NonequalAI 6-group") become rejection
tests, each with a well-defined variant that keeps what the spec was written to test. An
independent oracle already confirmed both traps on the full face, with gaps up to 56% and
61% (the LDH spec checked on its active conformation; `round4/reports/verify_guard.md`).

## 9. Engine, numerator, normalization and printing (phase 5)

### Engine (C1)

- Rapid-equilibrium segments contract to nodes; each form's weight is its segment root's
  weight times the equilibrium ratios along a content-ordered tree.
- The denominator sums, over segments, the segment's weight times its spanning trees of the
  segment graph (King–Altman/Cha).
- Exact polynomial arithmetic: checked `Int64`. Under the cap no coefficient leaves
  `Int64` (proof in the PR), so an overflow is an internal error.
- The size cap keeps today's measure (at most 5,000 spanning-tree products); for allosteric
  laws it applies to both conformations together.

### Numerator

The numerator is the rate at which the first substrate is consumed, summed over
steady-state steps, as today: each step contributes its net flux times the copies of that
substrate it effectively moves, including copies carried by rapid-equilibrium bindings on
either side (the RE-weight exponents). It needs no cut and ignores written direction.
After the solve is substituted, the numerator is divided exactly by the Haldane factor
`H = Π[S] − Π[P]/Keq`, so every law prints as `H` times a positive polynomial.

### Normalization: the unit rule (C4)

All denominator terms share one unit, each concentration counting as a dimension of its
own. Dividing numerator and denominator by a normalizer with that unit makes every
denominator term dimensionless; rule 4 covers a law that has no such normalizer:

1. **Constant term.** If the denominator has one, divide by it. Every term is then
   dimensionless; in a rapid-equilibrium law whose only form that binds nothing is free E,
   the constant term is free E's weight, and the denominator sums the enzyme forms'
   abundances relative to free E (Michaelis–Menten: `1 + S/K_m`).
2. **Cleland's normalizer.** Otherwise, when the denominator holds the term `ΠS/X` for
   every substrate `X` and the term `ΠS`, `ΠS` being the product of the `n` distinct
   substrates, divide by `N = Π_X c(ΠS/X) / c(ΠS)^(n−1)`, `c(m)` being the coefficient of
   the term `m` (ping-pong bi-bi: `K_a·K_b`). Each substrate term then reads as its
   concentrations over their Michaelis constants, `A·B·C/(K_a·K_b·K_c)`. A one-substrate law
   never takes this rule: its `ΠS/A` is the constant term.
3. **Virtual constant term.** Otherwise solve the units equation: weights `λ_m` over the
   denominator terms with `Σ λ_m = 1` and `Σ λ_m · a_m = 0`, `a_m` being a term's
   concentration-exponent vector; the normalizer is `Π c_m^λ_m`.
   Fractional weights are allowed and print as roots. To make the choice unique, terms are
   taken in a fixed order (lowest total degree first, then fewest non-substrate
   metabolites, then metabolite order: substrates, products, then the other names, each by
   name); the first affinely independent ones form the basis, and λ is the unique solution
   on that basis.
4. **No normalizer.** When no such weights exist, no product of the coefficients has the
   denominator's unit: some weighting of the concentrations gives every term of the
   denominator the same nonzero degree, as a rapid-equilibrium Theorell–Chance step
   `E(A) + B ⇌ E(Q) + P` does by making every term carry exactly one of `B` and `P`. The
   law then prints as derived, and its terms keep their units.

When the normalizer has a factor that is not a monomial, or an integer to a fractional
power, that factor prints once as `N = …`. An allosteric law normalizes each conformation's
polynomials by that conformation's own rule. The rule divides a conformation's numerator,
partition function and free-E weight alike, so the MWC combination, which reads their
ratios, is unchanged.

### Printing (E8, C10)

- `rate_equation_string` prints v in the independent parameters. Each Haldane and
  Wegscheider relation prints as a reference line.
- Next to it, the coefficient form (`coefficients = true`), in the order it evaluates: the
  `params` and `concs` lines; a comment line; `N` when it prints (`N_A` and `N_I` in an
  allosteric law); one line per coefficient in the fitted constants; in an allosteric law,
  the polynomials `Q_c`, `D_c` and `F_c` in those coefficients and the site factors `W_c`;
  and last `v = E_total·n_0·(A·B − P·Q/Keq) / (1 + d_A·A + d_B·B + d_A_B·A·B + …)`. The
  reference lines print only in the full law. A coefficient is named by the monomial it
  multiplies, a power as a superscript (`d_A²_B`; proposal, section 18).
- Allosteric laws print each conformation's polynomials, then the MWC combination.
- No `1 ^ n` or `L * 1` artifacts.

## 10. Dedup key and identify (phase 5)

### Key (C6)

1. Bring the law to lowest terms.
2. Label the fitted constants canonically: refine colours on the term–constant incidence
   and, while colours tie, individualize each constant of the first tied class in turn and
   keep the smallest text. The labelling is a canonical form, so the key does not depend
   on constant names, and it lists the constants in canonical order.
3. Label concentrations by metabolite name. For an allosteric law, the key's parts also
   hold each conformation's free-E weight, which the MWC combination reads.
4. Hash the canonical text with SHA-256 (the `SHA` standard library joins the package's
   dependencies).

Two mechanisms share a key exactly when their laws are identical up to renaming fitted
constants. Mechanisms that differ only in which metabolite binds (I versus J) get
different keys, and renaming a metabolite may change the key. Laws equal only after a
reparametrization (one constant standing for a product of two, for example) get different
keys; the census found such pairs, and they are fitted separately.

### identify

- The key replaces `eq_hash`: identify fits one representative per key and stores its fit
  by canonical position, the order of step 2 (C6). Every mechanism with that key reads the
  fit through its own labelling, so its CSV row carries its own parameter names; a test
  fits two same-key mechanisms whose fitted names differ. The cross-validation candidate
  list (`_offer_cv!`) deduplicates by key.
- One-subunit L phantoms (section 13) inherit the fit of the mechanism that holds the key of
  their law at L = 0 (often the parent) when that law holds every fitted constant but L and
  that fit was made in an earlier iteration or is made in their own, and are flagged; that
  mechanism has another key, so section 13 defines how the fit maps onto the variant.
  Otherwise, or when that fit throws, a phantom is fitted itself, unflagged (section 13).
- The beam counts each key once, so duplicates cannot crowd out distinct candidates, and a
  sweep expands every member it meets of a selected key, so no law reachable through one of
  them is lost; only an L phantom's key can be refused and later selected (section 13).

## 11. Code generation and speed gates (phases 0, 1 and 5)

- `rate_equation` is `@generated` over the compiled types. The generator decodes the Sig,
  solves, derives the law and emits a Horner body (C2), with locals for each dependent, the
  rapid-equilibrium segment sums, the Haldane factor and the MWC terms. Reciprocals and
  roots are computed once. Phase 0 (PR #85) already emits the Horner body, and dependents
  are locals today; phase 5 adds the segment, Haldane and MWC locals.
- The body allocates nothing and reads inputs by name, so a `NamedTuple` in any field order
  works. `fitted_params`, `parameters` and `metabolites` return constant tuples.
- **Gates (E1):** zero allocations for every tested law; at most 120 ns per call for
  `MECHANISM_TEST_SPECS`; for the three largest enumerated bi-bi laws
  (`_testhelper_largest_enumerated_mechanisms` from PR #85, 396 denominator terms each)
  and larger laws, a bound proportional to the number of denominator terms,
  `max(120 ns, c × terms)`. Its constant `c` is calibrated on CI hardware in phase 1
  (0.8 ns per term), on PR #85's Horner code, and re-calibrated in phase 5 once the locals
  land. The compile cost per fresh mechanism stays within `test_compile_budget`'s limits.

## 12. Fitting scale (phase 1)

In relative mode the loss cannot see a common factor on every steady-state rate constant,
so the fitter fixes it after the fit (E9).

- **API.**

  ```julia
  fp  = FittingProblem(m, data; Keq = 2.0)                    # relative (default)
  fp  = FittingProblem(m, data; Keq = 2.0, peak_range = (1e-9, 1e5))
  fp  = FittingProblem(m, data; Keq = 2.0, absolute = true)   # per-enzyme turnover data
  fit = fit_rate_equation(fp, optimizer)

  kcats = [(A = 10.0, B = 5.0) => 42.0,                       # conditions => turnover
           (A = 20.0, B = 5.0) => 44.0]
  p, spread = rescale_parameter_values(m, merge(fit.params, (; Keq = 2.0)), kcats)
  ```

  `IdentifyRateEquationProblem` takes `absolute` and `peak_range` too, and no kcat
  keyword. `scale_k_to_kcat` and `_kcat_forward` are removed (a breaking change).
- **Default scale.** In relative mode, `fit_rate_equation` and `identify_rate_equation`
  return constants scaled so the mechanism's peak forward rate per enzyme is 1. The peak is
  the largest net forward rate per enzyme (v/E_total) over every metabolite's
  concentration, products included, within `peak_range` (default `(1e-12, 1e3)`, in the
  data's concentration units; the units are a proposal, section 18); a deterministic
  multi-start search in log concentrations, range corners included, finds it. Every law is
  scaled. A law still rising at the range's upper edge, by more than 10% over the range's
  last decade, is flagged:
  the fit result gains a `peak_at_range_edge` field, and identify's CSV a column of the
  same name (tolerance and names proposed, section 18). The scale never turns a fit into a
  failure. Absolute mode never rescales.
- **Rescaling.** `rescale_parameter_values(m, params, kcats; peak_range = (1e-12, 1e3))`
  re-expresses a fit; `params` must include `Keq`, which the dependents need (proposal,
  section 18). With `kcats`, one or more measured turnovers (kcat or specific activity in
  turnover units) at stated concentrations, unlisted metabolites zero, it rescales to
  them: several anchors combine by the geometric mean of measured over predicted turnover,
  the log-space least-squares scale. An anchor at which the model predicts zero or
  opposite-signed turnover is an error. Without `kcats`, it applies the peak scale over
  `peak_range`. It returns `(params, spread)`: the rescaled constants and the standard
  deviation of the anchors' log(measured/predicted) turnover, 0 for one anchor and
  `nothing` without `kcats` (proposal, section 18).
- **Units.** Every function that takes `peak_range` documents that its default covers
  physiological concentrations in M, which the fitting tutorial recommends, or in mM; for
  data in µM, where 1e3 means 1 mM, raise the upper bound (proposal, section 18).
- **Absolute mode** keeps its loss and never rescales; only its selector changes, from
  `scale_k_to_kcat = nothing` to `absolute = true`.
- **Docstrings.** `identify_rate_equation`'s docstring calls `maxtime` the "max time per
  fit" (`src/identify_rate_equation.jl`); it applies to each restart, as
  `fit_rate_equation`'s docstring already says. Phase 1 fixes it and rewrites the
  docstrings of `FittingProblem`, `fit_rate_equation` and `IdentifyRateEquationProblem` for
  the new keywords.

## 13. Enumeration (phases 3, 4 and 5)

- Phase 3 ports to the graph `Mechanism`, with unchanged children, every move (the
  rapid-equilibrium to steady-state flip, the split, the dead-end binding, and the
  allosteric moves: tagging, adding a regulator, changing a tag, merging regulatory sites),
  seed construction (`init_mechanisms`, `seed_mechanisms`) and the emission rules
  (`_assert_emission_rules`). The slice prototyped only the flip, split and dead-end moves.
  Children change in phase 2, where the joint solve fixes known issues 1 and 2 and so
  changes the allosteric split gains, and in phase 4, where admission and the guard
  replace `_onlya_haldane_violation` and the tagging move replaces a V-type tag set that
  admission rejects by the smallest admissible ones (section 7).
- Every child passes the constructor, the well-posedness guard and, if allosteric,
  admission before identify sees it. Identify keeps today's `_degenerate` screen
  (section 8).
- No move tags a catalytic group `:OnlyI`; the regulator moves keep emitting `:OnlyI`
  site tags.
- **Islands that no tagged plain binding explains** are rejected by admission (phase 4): an
  island whose component free E's component does not reach along the arcs (section 7),
  such as the covalent island of a ping-pong reached through chemistry alone. The tagging
  move no longer emits that tag set (section 7).
- **Duplicate tag sets** (siblings that differ only in the tags of groups the inactive
  conformation never reaches) stay in the frontier and share one fit through the key
  (phase 5).
- **One-subunit L phantoms** stay in the frontier, inherit the fit of the mechanism that
  holds the key of their law at L = 0, their source (often the parent), when that fit was
  made in an earlier iteration or is made in their own, and are flagged (phase 5). A
  variant is an L phantom when an exact rescaling certificate in `src/` (`key.jl`) shows
  that its family is the family of its law at L = 0. The certificate alone decides: it
  needs one catalytic subunit and an inactive conformation that turns nothing over,
  whatever that conformation binds, and nothing screens a variant first; identifiable ranks
  are test oracles. As L → 0 the variant's law becomes the source's, up to the names of
  its constants, so its reported fit is the source's constants under the variant's names
  with L = 0, and the source's loss. In the beam and in cross-validation it counts its
  own `length(fitted_params)` (B9). Cross-validation does not refit it either: on every
  fold its best fit equals the source's, so it takes the source's fold losses. Its CSV row
  records the fit source's `mechanism_id` and the L-phantom flag (section 5). A phantom
  is fitted as an ordinary candidate, unflagged, when no mechanism of the search has fitted
  its law at L = 0 by the end of the phantom's own iteration, when that fit throws, or when
  that law lacks one of the phantom's fitted constants other than L.

## 14. Testing

- **Oracles: reuse and extend.** The suite already has the ODE steady-state check on every
  plain spec (`_testhelper_ode_steady_state_flux`), the analytical textbook laws in
  `MECHANISM_TEST_SPECS`, the two-conformation mass-action gates in
  `allosteric_ground_truth.jl` (with the oligomer oracle for 1–3 protomers), and
  detailed-balance tests. They all stay; the ODE check remains an independent cross-check.
  The rewrite adds:
  - one exact linear steady-state solver (rational or BigFloat) that builds its network
    from any mechanism, so new gates need no hand-written networks;
  - the RE-speed oracle: the law is the limit of the mass-action steady state as RE speeds
    grow, at interior points and on faces, under several speed patterns (today's ODE
    checks sample concentrations of 0.1–10 only);
  - the ε-path oracle: moving `:OnlyA` or `:OnlyI` constants along the certified direction
    converges to the law;
  - detailed balance at equilibrium concentrations over whole populations;
  - the exact LP of admission's existence check, as an oracle for the graph test.

  The reference branch (Appendix B) holds the probes' versions of these oracles. The tests
  use these tolerances (proposal, section 18): the exact solver agrees with the law to
  1e-20 relative in 256-bit BigFloat (the slice's worst error was 7e-31); the RE-speed
  oracle runs every RE step at Λ times a random factor, log-uniform over 10^±3 (the
  verifier's spread; the probe's `round4/guard/lib.jl` used 10^±1), at Λ = 1e10 and 1e20
  with two speed patterns, and the law agrees at Λ = 1e20 within 1e-8 relative; along the
  ε-path the error shrinks in proportion to ε, from ε = 1e-3 to 1e-12.
- **Populations** are regenerated inside the tests from `init_mechanisms` and
  `expand_mechanisms`, sampled to keep the population tests within about one minute of the
  full suite (proposal, section 18), never loaded from scratch files. Over them the tests
  check oracle agreement, determinism under rewrites, regroupings and renames, and key
  stability.
- **Textbook bridge (C5)** stays in the tests, rewritten by role in phase 2, since it needs
  `all_parameters`: each test reads its written steps and takes the package constant of the
  matching stored step in the written direction. Phase 3 removes the DSL's test-only
  source-order feature.
- **Enumeration tests** keep CLAUDE.md's rules: inline fixtures and exact children.
- **Pins.** Each phase updates the pins of the behavior it changes, in the same PR. Before
  changing a test that checks only a string, the phase adds an oracle value check. Tests
  pin `length(fitted_params(m))` (B9).
- **Gates.** Section 11's speed, allocation and compile-budget tests, and the constructor
  benchmark (section 2).
- **By phase.** Phase 2 adds the exact steady-state solver, detailed balance over
  populations and the rewritten bridge; phase 4 adds the RE-speed and ε-path oracles, the
  LP oracle, the face tests, R2 and the AD fixtures; phase 5 adds the key's stability tests
  and the test that two same-key mechanisms read one fit under their own names.

## 15. Phases

Each phase is one PR, green with the full suite and with a test for every behavior it
decides, proposing a version bump for Denis to decide (breaking changes are acceptable
before 1.0). Before phase 1 (D1):

- PRs #85 (Horner code generation), #86 (one flip segment, ROADMAP bug list) and #87
  (`:OnlyA` guard crash) merge;
- the fix for the ter-bi expansion abort (one oversized split candidate aborts a parent's
  whole expansion) lands as its own PR; it can stay minimal, since the abort comes through
  the rename pass that phase 2 removes;
- the ROADMAP is brought in line with this spec (section 16), in PR #86 before it merges or
  in a PR of its own.

| Phase | Contents | Gate or dependency |
|---|---|---|
| 0 | Horner code generation (PR #85) | merges before phase 1 |
| 1 | Fitting scale; speed-gate calibration on CI (sections 11, 12) | the pre-phase PRs |
| 2 | Joint solve, names, parameter API, single-graph allosteric derivation, on today's types; exact and detailed-balance oracles; textbook bridge (sections 6, 7, 14) | none |
| 3 | Graph `Mechanism`, compiled types and Sig, DSL, `show`, parser, CSV, constructor checks, free E; every move, seed construction, today's derivation, fitting and identify ported with unchanged children; source-order feature removed; constructor benchmark (sections 4, 5, 13, 14) | GC-crash probe (E6) |
| 4 | Admission, flip at E, `:OnlyI` rule, well-posedness guard, enumeration checks; RE-speed, ε-path and LP oracles, face tests (sections 7, 8, 13, 14) | face check exact first; admission's face condition sound, conservative on R-B1 and R-B2, listed |
| 5 | Engine, numerator, normalization, printing, key, identify's key cache and beam, L-phantom flag, code-generation locals, speed-gate re-calibration (sections 9, 10, 11, 13) | none |

Phase 3 moves today's constructor-side checks (the bottomless segment,
`_onlya_haldane_violation`, PR #86's one-flip-segment check) into one boundary function
that the macros, the parser, the enumerator and the derivation entry call; phase 4 replaces
its contents with the guard and admission (proposal, section 18). Fitted names change
twice: in phase 2 (E3, E7, B10) and again in phase 3 (modifier and residual names), so
identify's CSV parameter columns change in both. Each phase updates the docs and ROADMAP
entries it affects (section 16).

## 16. Docs and ROADMAP

- **Docs.** The docs build runs the doctests and checks that every export has a docstring
  (`docs/make.jl`), so each phase rewrites the pages it makes false:
  - phase 1: `fitting/normalized_vs_absolute.md` (default scale, rescaling after the fit,
    absolute mode), `fitting/tutorial.md`, `fitting/loss_and_optimizers.md`,
    `developer.md` (`_kcat_forward`), and the HPC scripts
    (`docs/hpc_results/*/identify_*.jl`, which pass `scale_k_to_kcat`);
  - phase 2: `deriving/thermodynamic_constraints.md` (the fixed priority, and that
    renaming a metabolite can change which constant is dependent, C8),
    `deriving/textbooks.md` (Full-mode section removed), and the "Canonical Step Form"
    section of CLAUDE.md, whose claim that the reduction picks dependents by step order
    becomes false;
  - phase 3: the macro docstrings in `src/dsl.jl`, rendered on `api.md` (two macros,
    `modifiers:`, `::Mod`, `regulatory_site` lines, free E), `deriving/regulators.md`
    (`regulators:`, `::Inh`), `deriving/ping_pong.md` (residual names), the DSL examples
    of `deriving/mwc_allostery.md` and `identify/tutorial.md`, the `::Inh` text of
    `developer.md`, the CSV format in `identify/tutorial.md` and
    `identify/model_selection.md`, and the "API Design" section of CLAUDE.md;
  - phase 4: `deriving/mwc_allostery.md` (flip at E, limits, admission) and the `:OnlyI`
    rule in the macro docstrings;
  - phase 5: `deriving/textbooks.md` and `deriving/ping_pong.md` (normalization,
    coefficient form), the identify pages (`identify/tutorial.md`,
    `identify/model_selection.md`, `identify/enumeration_engine.md`: key, inherited fits,
    L phantoms) and `developer.md` (architecture).

  Every phase also updates any other page that shows a name, a printed law or an API the
  phase changes, `README.md`'s quickstart and `getting_started.md` among them.
- **ROADMAP.** Before phase 1 the ROADMAP is brought in line with this spec (section 15):
  PR #86's Planned graph-rewrite entry, which still describes one macro, kcat as the
  analytic limit and an older phase order, becomes a pointer to this spec; the `:OnlyA`
  entry ("Derive both conformations from one catalytic rate equation") flips at free E
  instead of at one flip segment; the Decided-against kcat-anchor entry goes, since E9
  adopts anchors after the fit; the Planned L-phantom entry keeps and flags the variants
  instead of dropping them (section 13); and entries for the `maxtime` docstring and the
  ter-bi abort are added, which phase 1 and the ter-bi fix retire. Each phase then retires
  or rewrites the entries it settles, as CLAUDE.md requires:

  | Phase | ROADMAP entries it retires or rewrites |
  |---|---|
  | 1 | the `maxtime` docstring; `_kcat_forward` reports "no kcat components" (#86) |
  | 2 | known issues 1 and 2; Full mode lists rapid-equilibrium constants (#86); splitting a tied kinetic group changes fitted names and `eq_hash` (#86: the names here, `eq_hash` in phase 5); Planned "Derive both conformations from one catalytic rate equation" (its `:OnlyA` half moves to phase 4); Planned "Build the Wegscheider rename into the constraint columns" |
  | 3 | Planned residual names; Idea "Encode a mechanism's Sig without compiling per mix of step shapes"; Decided against "Have the `Mechanism` constructors check step metabolites against the reaction" (moot: the graph type holds no reaction) |
  | 4 | the `:OnlyA` guard misjudges mixed groups (#86); the `:OnlyA` half of the conformations entry |
  | 5 | the `eq_hash` half of the regrouping entry; Planned L-phantom entry; Idea "Count identifiable parameters with an exact modular rank" (rewritten: the L-phantom test lands in `src/`, other uses stay ideas) |

  Known issue 3 and PR #86's rapid-equilibrium isomerization entry stay (section 2); the
  `eq_complexity_filter` and GC live-byte entries are open (section 18).

## 17. Risks and how each is handled

| Risk | Handling |
|---|---|
| The slice's unexplained "GC error (probable corruption)" in a workload that deserializes mechanisms | A probe of at most half a day before phase 3 (E6): rebuild instead of deserializing, run on Julia 1.10 and on x86, report upstream if Julia is at fault |
| The face check or admission's face condition rejects valid mechanisms or misses invalid ones | Phase 4 makes the face check exact against the RE-speed oracle and admission's face condition sound against the ε-path oracle; conservative on R-B1 and R-B2, listed in ROADMAP |
| The L-phantom detector costs or errs | Phase 5 measures it; oracle tests compare flagged variants with the parent's law |
| The peak search misses the peak | Deterministic starts, range corners, tests against laws with known peaks (substrate inhibition, product-assisted release) |
| Test-suite time grows | Population tests are sampled to about one minute of the full suite (section 14) |
| Parts the slice never prototyped (Appendix B) behave unlike its measurements | Each phase's plan budgets their validation and tests them against the oracles |
| Probe artifacts vanish with the session | This spec records every decision; the reference branch archives the reports, prototype code and decision records (Appendix B) |

## 18. Open for Denis

Denis has not confirmed the author's proposals below. Implementation follows each one
unless Denis objects; the sections cited mark them "proposal".

1. **Exports (section 6).** `parameters` and `all_parameters` are exported, and
   `fitted_params` stays internal. `parameters` and `fitted_params` keep today's status;
   exporting `all_parameters` adds public API that B5 described as a test aid.
2. **Constructor and guard split (sections 3, 15).** The constructors enforce
   representation invariants only; the guard and admission run at the DSL boundary, on
   every enumerated child and at the derivation entry. Phase 3 moves today's
   constructor-side checks into one boundary function, which phase 4 fills with the guard.
   Every user path reaches the DSL or the derivation entry, so no input skips a check.
3. **Source layout (section 3).** One concern per file, as listed; the plans may merge
   small files.
4. **Peak scale (section 12).** Denis set the default `peak_range = (1e-12, 1e3)`. The
   proposals: the range is in the data's concentration units, and the docstrings say the
   default covers physiological concentrations in M, the unit the fitting tutorial
   recommends, or in mM. A law is flagged when it still rises by more than 10% over the
   range's last decade; for a hyperbolic law that means an asymptote more than about 1%
   above the range's peak, so ordinary saturating laws stay unflagged.
5. **`rescale_parameter_values` (section 12).** `params` must include `Keq`, and the
   function returns `(params, spread)`, the spread being the standard deviation of the
   anchors' log(measured/predicted).
6. **`oligomeric_state:` (section 5).** It takes an Int or a tuple, which keeps the search
   over several catalytic multiplicities that `allowed_catalytic_multiplicities:` offers
   today. The alternative drops that search.
7. **Macro results (section 5).** Both macros return the compiled types, as today, which
   keeps `rate_equation` on macro output allocation-free and the docs' `m isa
   EnzymeMechanism` true. The alternative returns the value types and changes the exports
   and the docs.
8. **Modifier-zero faces (section 8).** E2 listed reactant faces only. The guard also
   checks each face, and the interior, with every modifier zero, because most data carry no
   modifier and the support test already covers those faces.
9. **Names and fields (sections 5, 9, 12).** The CSV columns `mechanism`, `mechanism_id`
   (16 hex digits of the SHA-256 of the DSL text), `parent_id`, `key`, `fit_inherited`,
   `fit_source_id`, `L_phantom` and `peak_at_range_edge`; the fit result's
   `peak_at_range_edge` field; coefficient names for a power, such as `d_A²_B`.
10. **Test budget and tolerances (section 14).** Population tests take about one minute of
    the full suite. The exact solver agrees to 1e-20 relative; the RE-speed oracle draws
    speed factors log-uniform over 10^±3, runs at Λ = 1e10 and 1e20 with two patterns and
    agrees within 1e-8 at Λ = 1e20; the ε-path error falls in proportion to ε from 1e-3
    to 1e-12.
11. **Two ROADMAP entries (section 16).** The size cap counts both conformations
    (section 9), so `eq_complexity_filter` counts them too from phase 5, which retires its
    known issue. Phase 2 re-measures the GC live-byte growth, which its checked `Int64`
    solve may remove, and retires that entry if it is gone.
12. **Rapid-equilibrium isomerization gauge (section 2).** D1's chosen option fixes the
    probes' bugs in the redesign phases, and its list includes this gauge in 5 enumerated
    ping-pong mechanisms; B9 keeps the identifiable rank an analysis tool and counts the
    gauge as part of the declared model. The spec follows B9: the gauge stays out of
    scope and the rank analysis tool reports it. The alternative gives it a failing test
    and a fix in the phase that rewrites the solve, as D1's option describes.
13. **Allosteric ping-pong (section 7).** Denis requires that ping-pong mechanisms stay
    allosteric. The proposal: the tagging move replaces a V-type tag set that admission
    rejects by its smallest admissible ones, which stops one half-reaction per child on a
    ping-pong; stopping both halves waits for iso mechanisms. The alternative expands
    rejected parents without fitting them, so the change-tag move still reaches
    `:NonequalAI` chemistry, but then no ping-pong with a blocked inactive conformation is
    ever fitted.

### Review proposals kept as decided

- The consistency review proposed eliminating isomerizations, then plain bindings, then
  other chemistry. The spec keeps the probed order approved with the joint solve:
  isomerizations, then other chemistry, then plain bindings (section 6).
- The consistency review proposed rejecting the two trap specs only if the exact face
  check and the RE-speed oracle find them ill-defined. The spec keeps the decision to
  reject both (E2); an independent oracle already confirmed their gaps (section 8).
- The codebase review proposed splitting phase 3 into the DSL, `show` and parser on
  today's types, then the type swap, since the tests hold 222 `catalytic_steps:` and 162
  `regulators:` labels, 283 lines with `::Inh` and about 370 direct `Mechanism(` and
  `Step(` calls. The spec keeps phase 3 as one PR, as decided.

## Appendix A: Decision index

| Card | Decision | Section |
|---|---|---|
| A1 | Proceed: the redesign is for correctness, determinism and speed, not size | 2 |
| A2, E11 | Phased PRs in the critic's order; the tests-only phase was dropped in brainstorming | 15 |
| A3 | Vertical-slice probe run | Appendix B |
| B1 | `:OnlyA` bindings by affinity, chemistry by kinetic block; island rule; cone test | 7, 13 |
| B2 | Per-step classes in the excluding conformation; mixed groups allowed | 6, 7 |
| B3 | Superseded by E10 | 7 |
| B4 | Superseded by E9 | 12 |
| B5 | No all-parameters law; `all_parameters` replaces the test helper | 6 |
| B6 | Rational exponents accepted | 6, 9 |
| B7 | Tie forward with forward | 4, 6 |
| B8 | Flux-free steady-state steps allowed | 6, 8 |
| B9 | Identifiable rank is an analysis tool; one enumeration use (L phantoms); tests pin `fitted_params` lengths | 2, 13, 14 |
| B10 | `:EqualAI` regulator constants bare | 6 |
| C1 | Lean engine | 2, 9 |
| C2 | Horner with segment, Haldane and MWC locals | 11 |
| C3 | Mechanism as the only type parameter; derive in the generator | 4, 11 |
| C4 | Constant term, else virtual constant term; refined into the unit rule | 9 |
| C5 | No textbook layer in the package; textbook bridge in the tests, by role | 2, 14 |
| C6 | Rename-proof key; concentrations labelled by name; fits stored by canonical position | 10 |
| C7 | No source-level duplicate suppression | 2, 13 |
| C8 | Priority ties broken by content names, and documented | 6, 16 |
| C9 | Encoding limits with clear errors; orientation ties by form index | 4 |
| C10 | Coefficient form named by monomial | 9 |
| D1 | Cheap fixes first (PRs #85, #86, #87) | 2, 15 |
| E1 | 120 ns gate plus a size-scaled bound, timed on CI hardware from the first PR | 11 |
| E2 | Fit what is well defined, reject what is not | 8 |
| E3 | Names after the kept member | 6 |
| E4 | Alphabetical parameter order | 6 |
| E5 | Reserved names; no `_` in metabolite names | 4, 5 |
| E6 | GC-crash probe before phase 3 | 15, 17 |
| E7 | "Against" means against the reaction direction | 6 |
| E8 | Printed laws in independent parameters | 9 |
| E9 | Relative fits scaled to peak forward rate 1 within `peak_range` (default `(1e-12, 1e3)`), a keyword wherever the scale is applied; measured kcats applied after the fit | 12 |
| E10 | Flip at free E; free E mandatory in every mechanism | 4, 7 |

Decisions without a card, from the brainstorming sessions and the chat-decision ledger: one
umbrella spec, then one plan per phase executed test-first with subagent-driven
development; the slice is reference only, and the oracles the tests need move into
`test/`; two macros with one vocabulary and composed types; `modifiers:` and `::Mod`;
`:OnlyI` limited to modifier bindings, with admission mirrored in the active conformation,
which must still turn over; the face check made exact before the guard ships; today's
numerator, divided exactly by the Haldane factor; the unit rule; kcat rescaling after the
fit; islands that no tagged plain binding explains rejected by admission; duplicates and L
phantoms kept with inherited fits, and no separate PR that deduplicates tag sets in
today's package; the beam counts keys once; the `maxtime` docstring fixed in the fitting
phase.

## Appendix B: Evidence

The branch `graph-rewrite-reference` (commit `b259716`) archives the material behind this
spec under `reference/`: the vertical-slice code and its interface notes (`slice/src`,
`slice/INTERFACES.md`), every probe, verification and census report (`slice/reports`,
`probes`, `round4/reports`, `census`, `peaksearch`, `followups`, `explore`), the oracles
the tests will reuse (`round4/guard/lib.jl`, `round4/flipE/scripts/lib.jl`), and the
decision records (`decisions/`: the card page, Denis's choices and notes, the replies and
the chat-decision ledger). `reference/README.md` maps each folder to this spec's sections.
The branch is reference only and is never merged.

Key results:

- **Slice.** 237,965 laws from 257,698 mechanisms (all 250,855 ter-ter seeds included)
  matched their oracles. The other 19,733 were refused: 19,376 ter-ter laws over the
  5,000-term size cap, the rest by the input checks or admission. Worst error 7e-31 on the
  core sets (4,917 records in 27 sets: enumerated populations, the package's test specs and
  the probes' adversarial and hand-built cases). Dependents identical across 24,580
  rewritten inputs and 12,608 same-model regroupings. Admission's existence test equalled
  the exact LP on 2,995 of 2,995 allosteric records. Code size about 2,045 → 1,980 lines;
  derivation 926 → 1,110. The largest bi-bi laws ran in about 42 ns; 343 ter-ter laws ran
  over 120 ns.
- **Not prototyped.** The slice left out E3's naming (it named classes by today's rule),
  E4's order (it used content order), E9's scale (it used the analytic kcat), admission's
  face condition, the key's name labels (it labelled concentrations by position), the size
  cap over both conformations (it capped each apart), four of the seven moves (tagging,
  adding a regulator, changing a tag, merging sites), seed construction and the L-phantom
  detector.
- **kcat = 0 laws.** All 284 are complete substrate inhibition, exactly: positive rate at
  every finite substrate concentration, falling to 0 only as one substrate grows.
- **Peak search.** The textbook steady-state random-order bi-bi never exceeds its
  both-saturated rate with no product present (proved two ways); with a product present it
  can (Q rescuing a slow E(P) release, 125.6x in a simple example).
- **Well-posedness.** Among 13,645 records (8,626 of them mutants of enumerated parents):
  1,654 records with RE turnover, all mutants but one (the slice's law v ≡ 0, or in 3
  records finite and wrong); 1,027 mutants with a deleted catalytic step whose accepted
  law is v ≡ 0; 21 mutants whose law is NaN on a face (22 with a hand-built case outside
  these records); face splits and traps with gaps of 2.3% median, 47% at the 90th
  percentile and 314% at most, and up to 80% with one product at zero. Enumerated output
  never produced v ≡ 0, RE turnover or a 0/0 face.
- **Flip at E.** Exact against an independent oracle (12 shapes, 1 to 3 subunits, worst
  3.2e-69); no change for 9,308 enumerated mechanisms.
- **Census.** Islands that admission rejects take 8–32% of a round's fit time in ping-pong
  and the ter reactions; duplicate tag sets 8–21%; one-subunit L phantoms 10–15% in bi-bi
  and ping-pong (`census/reports/critic.md`, each share net of the classes before it).

## Appendix C: Glossary

- **Anchor.** A measured turnover at stated concentrations, used to rescale a relative fit.
- **Certified direction.** The direction of vanishing exponents that admission's existence
  check finds; the ε-path oracle moves the vanishing constants along it.
- **Chemistry.** In the limit rules, every step that is not a plain binding: an
  isomerization, a fused binding or a Theorell–Chance step.
- **Conformation.** One of the two MWC conformations of an allosteric mechanism, active or
  inactive; each runs the same graph. Not a form's conformation label.
- **Conformation label.** The label of an enzyme form (`E`, `F`, `Eprime`, `Estar`); free E
  carries `E`. A step that changes it makes a conformational change.
- **Face.** A set of concentration points where some metabolites are exactly zero, such as
  every product zero (initial-rate data).
- **Fused binding.** A binding into a merged complex, which also does chemistry
  (`E(A) + B <--> E(P, Q)`).
- **Group labels.** The identifiers that assign steps to kinetic groups; renumbering them
  changes nothing.
- **Island.** A form free E cannot reach in a conformation once that conformation's
  excluded steps are removed. Admission accepts an island only behind a tagged plain
  binding, on its bound side (section 7).
- **Isomerization.** A step that takes up and gives off nothing: a conformational change or
  on-enzyme chemistry (`E(A, B) <--> E(P, Q)`).
- **Iso mechanism.** A mechanism whose turnover passes through a free form with another
  conformation label, such as Cleland's iso uni-uni (`E + S ⇌ E(S)`, `E(S) <--> F + P`,
  `F <--> E`).
- **Kinetic group, class.** Steps that share constants. A class is the columns of one group
  slot, tied by group rows; it takes the name of its column that comes last in elimination
  order, the one the solve keeps when it keeps any (E3).
- **L phantom.** An allosteric variant in which L enters the law only through a combination
  other constants absorb, so the data cannot determine it.
- **Model.** The family of rate laws a mechanism gives over all values of its constants.
- **Mutant.** A mechanism a probe derived from an enumerated parent by changing RE/SS flags
  (one group either way, two groups, or every step to RE) or by deleting one catalytic step.
- **Peak turnover.** The largest net forward rate per enzyme over every metabolite's
  concentration within `peak_range`, products included.
- **Physical-forward.** The direction in which a step is stored and read: a binding as
  uptake; any other step in the direction that turns substrate atoms into product atoms
  (orientation tier 1), else from the form where products leave toward the form where
  substrates bind (tier 2), else by form index (tier 3, C9).
- **Plain binding.** A step that takes up one metabolite and gives off none, its bound form
  holding the other form's metabolites plus that one, with the same residual; the
  conformation label may change.
- **Same-model regrouping.** Splitting or merging a kinetic group whose members the cycle
  rows already tie; the law does not change.
- **Segment.** The forms that rapid-equilibrium steps join: a connected component of the
  graph of RE steps.
