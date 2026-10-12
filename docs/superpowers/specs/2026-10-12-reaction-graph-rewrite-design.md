# Reaction-graph rewrite: design

Status: draft for review, 2026-10-12. This spec is the written record of the decisions made
between 2026-10-09 and 2026-10-12: decision cards A1–E11 and the brainstorming sessions.
Appendix A maps each card to the section that carries it. Implementation follows one plan
per phase (section 15), each written with the writing-plans workflow and executed test-first.

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
fitting scale. A vertical-slice prototype of this design derived 237,965 laws from 257,698
mechanisms, and every law matched an exact oracle (Appendix B).

## 2. Goals, scope and success criteria

### Goals

- **Correct laws.** Every derived law equals the steady state of mass-action kinetics in the
  limit the mechanism states (fast rapid-equilibrium steps, vanishing `:OnlyA` or `:OnlyI`
  steps), and satisfies detailed balance.
- **Deterministic results.** A mechanism's law, fitted names, parameter order and dedup key
  depend only on the mechanism, never on step order, step direction, group labels or
  process.
- **Readable output.** Printed laws use independent parameters, dimensionless denominator
  terms and a Cleland-style coefficient form. `show` prints the DSL.
- **One rule per question.** One well-posedness guard decides which mechanisms are
  admitted; one admission rule decides which allosteric tag sets are; one scale rule
  reports fitted rate constants.
- **Fast evaluation.** `rate_equation` stays allocation-free and within the speed gates
  (section 11).

### Out of scope

- Enumerating iso mechanisms (the DSL accepts them; no move generates them).
- Hill's cycle-flux numerator and block factorization (C1).
- Generating fewer duplicates at the source (C7); the key and the fit cache handle them.
- The iso dwell-time-gauge phantom and two-conformation kinetic cooperativity (ROADMAP).
- Using the identifiable rank in parameter counts or model selection (B9); the one
  enumeration use is the L-phantom flag in section 13.

### Success criteria

Each criterion is a test in the suite.

1. **Exactness.** On populations regenerated in the tests (section 14) and on adversarial
   fixtures, every law matches an exact mass-action oracle, plain and two-conformation;
   rapid-equilibrium laws match the RE-speed oracle as the RE speeds grow, including on
   zero-concentration faces; `:OnlyA` and `:OnlyI` laws match the ε-path oracle.
2. **Detailed balance.** At equilibrium concentrations the net rate is zero, per
   conformation and combined.
3. **Determinism.** Rewriting a mechanism (step order, step direction, group labels)
   changes no law, fitted name or key; a same-model regrouping changes no fitted name;
   renaming metabolites changes no key; a second process gives identical keys and
   generated code.
4. **Known issues 1 and 2 are gone.** No phantom names and no fold hazards.
5. **Speed and compile cost.** Section 11's gates hold.
6. **Coverage.** Every decided behavior has a test; every phase PR is green with the full
   suite.

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
law (engine) ──► normalization and printing       (rate_equation_string, coefficient form)
        │   ──► dedup key                         (identify's fit cache and beam)
        │   ──► code generation                   (rate_equation, @generated over the Sig)
        ▼
fitting (FittingProblem, fit_rate_equation) ──► scale (peak turnover; kcat rescale after fit)
        │
        ▼
identify (enumeration moves, beam, cross-validation)
```

Planned source layout, one concern per file: `types.jl` (value types, canonical form,
representation checks, Sig), `names.jl` (the parameter-naming chokepoint), `dsl.jl`
(macros, parser, `show`), `poly.jl` (exact polynomials), `thermo.jl` (joint solve),
`admission.jl` (guard and admission), `engine.jl` (law), `normalize.jl` and `print.jl`,
`key.jl`, `codegen.jl`, `fitting.jl` (with the scale), `mechanism_enumeration.jl`,
`identify_rate_equation.jl`. The plans may merge small files.

The constructors enforce representation invariants only. The guard and admission run at
the DSL boundary (macros and parser), on every enumerated child and at the derivation
entry, never inside internal constructions. Internal code that builds a sub-graph reads
the graph directly instead of constructing a mechanism, so a check meant for user input
cannot fire on it (the class of bug PR #87 fixed).

## 4. Mechanism types (phase 3)

### Mechanism

`Mechanism` is the reaction graph:

- substrate and product names;
- a metabolite table; each entry is a name and whether it binds as a modifier, so a
  reactant bound as a modifier (`P::Mod`) is its own entry;
- conformation labels;
- enzyme forms (bound metabolites, residual, conformation), in canonical order with free E
  as form 1;
- edges (two form indices, consumed and released metabolites, an RE flag), stored
  canonically;
- the kinetic group of each edge.

`AllostericMechanism` holds a catalytic `Mechanism` plus one tag per kinetic group, the
catalytic multiplicity and the regulatory sites (composition). A future model with its own
derivation, such as KNF, wraps the same graph in its own type.

Invariants:

- **Canonical form.** Equal mechanisms have equal fields however they were written.
  Step direction follows today's orientation tiers; ties fall to form index (C9).
- **Content identity.** `==` and `hash` compare content and are stable across processes.
- **Free E.** Every mechanism contains a free form whose conformation label is `E`
  (nothing bound, no residual); it is form 1 (E10).
- **Limits.** At most 32 metabolite-table entries and 15 copies of one entry on a form or
  step side, with errors that say so (C9).
- **No `EnzymeReaction` field.** Atom checks take the reaction's atom table from the
  caller. `EnzymeReaction` stays as the enumeration input.

### Constructor checks

The constructor rejects:

- missing or overlapping substrates and products; a residual that adds a product or
  subtracts a substrate;
- an undeclared name, or a name bound as a modifier without being declared one;
- a self-loop, a metabolite on both sides of a step, one reaction in two steps;
- group members that, in their physical-forward direction, take up or give off different
  metabolites, differ in RE/SS flag or make different conformational changes;
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
`AllostericMechanism`. They share the labels `substrates:`, `products:`, `modifiers:` and
`steps:` and the step syntax. Only `@allosteric_mechanism` takes `catalytic_multiplicity:`,
`allosteric_regulators:`, `regulatory_site(n):` lines and a tag on every step or group.

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
  step tags, residuals (`E(Q; residual = A - P)`), conformation labels such as `Eprime`,
  `Estar` or `F` for forms other than free E.
- Removed: `catalytic_steps:` (now `steps:`), `catalytic_inhibitors:` and the plain
  macro's `regulators:` (now `modifiers:`), the parenthesized single-step tag form, and the
  split allosteric declaration with implicit default sites.

`@enzyme_reaction` (enumeration input) keeps `substrates:` and `products:` with atoms,
`competitive_inhibitors:`, `allosteric_regulators:` (with `::Activator`/`::Inhibitor` and
per-regulator multiplicities such as `Z(1, 2)`), `oligomeric_state:` and
`shared_catalytic_site:`. It drops `dead_end_inhibitors:` (now `competitive_inhibitors:`),
`allowed_catalytic_multiplicities:` (now `oligomeric_state:`) and multi-argument atom
brackets (`B[N, P]`).

### Display and CSV

- `show(m)` prints exactly the DSL text that builds `m`, with the macro that built it. The
  compiled types print the same text.
- An internal parser reads that text back by walking the parsed expression, never with
  `eval`, and returns a mechanism equal to `m`.
- identify's CSV stores the DSL text, a short stable `mechanism_id` (from the content hash)
  and the parent's `parent_id`, in place of today's type strings (about 2,700 characters,
  79% of a results file). It keeps recording whether a row's fit is inherited, and from
  which `mechanism_id`.

## 6. Thermodynamics, names and the parameter API (phase 2)

### One joint solve

- **Columns.** One bare constant per step slot: a steady-state step has a forward and a
  reverse rate constant; a rapid-equilibrium step has one equilibrium constant (a binding's
  K is its dissociation constant). Members of a kinetic group share columns slot by slot,
  aligned forward with forward (B7).
- **Allosteric columns.** The inactive conformation adds none for an `:EqualAI` group (it
  reads the active column), one per `:NonequalAI` group, and one vanishing copy per member
  of an `:OnlyA` group (B2). An `:OnlyI` group's finite constants belong to the inactive
  conformation, and the active conformation adds one vanishing copy per member.
- **Rows.** The fundamental cycles of the graph, in each conformation: a Haldane row
  (`= log Keq` times the cycle's multiple of the net reaction) or a Wegscheider row
  (`= 0`). Inactive cycle rows cover the whole graph, islands included.
- **Solve.** One exact rational row reduction (checked `Int64`, retried in `BigInt` on
  overflow) chooses the dependents and expresses each as a power product of independents
  and `Keq`. Rational exponents such as `Keq^(1/2)` are accepted (B6). There is no second
  rename pass.

### Fixed priority

Columns are ordered by structure only:

1. the step's kind: isomerization, plain binding, or other chemistry;
2. whether all the step's forms are occupied;
3. steady-state before rapid-equilibrium;
4. the constant that runs against the reaction direction first (E7): the reverse of the
   physical-forward direction, except that for a product binding the uptake runs against.

Columns outside the law (inactive copies on pruned forms, vanishing `:OnlyA` and `:OnlyI`
copies) come first so they are eliminated first. Ties fall to canonical content: the
step's forms and metabolites by name (C8). Phase 2 computes this on today's types and phase 3
on the graph type with identical results, so phase 3 cannot change the dependents.

### Names

- A rate constant is named after the direction it drives (`k_EA_to_EP`); a binding's K
  after its release (`K_EA_to_E_A`); any other K after its forward direction.
- Each kinetic class is named after the member the solve keeps (E3), so a name depends on
  the model alone.
- `A_`/`I_` state tags appear except on `:EqualAI` groups and plain mechanisms;
  `:EqualAI` regulatory sites print bare (`K_Rreg`, B10). Fixed names: `L`, `Keq`,
  `E_total`.
- Each vanishing `:OnlyA` or `:OnlyI` copy is named after its own step. These names occur
  only inside the solve and in admission messages; no law contains them.
- All `Parameter → Symbol` rendering goes through the naming chokepoint, as today.

### Parameter API

- `parameters(m)`: the independent constants plus `Keq` and `E_total`. Exported.
- `fitted_params(m)`: the independent constants. Internal, as today.
- `all_parameters(m, params)`: every step constant, dependents computed from the
  independent values in `params`. Exported; it replaces the test helper
  `_testhelper_compute_all_params`.
- Every list is alphabetical (E4).
- Removed (B5): the `Full`/`Reduced` modes and exports, their types, the three `Full`
  methods of `parameters`, `rate_equation` and `rate_equation_string`, and the Full-mode
  section of the textbooks page.
- Flux-free steady-state steps stay allowed (B8).

## 7. Allosteric model and admission (phases 2 and 4)

### Model

- The conformations interconvert only at free E (form 1), with `L = [E_I]/[E_A]` (E10).
  Each conformation runs its catalytic cycle on the same graph; the rate combines them with
  today's MWC formula (`docs/src/deriving/mwc_allostery.md`), exact at one catalytic
  subunit and the documented commit-per-cycle model above that.
- Catalytic groups take `:EqualAI`, `:NonequalAI` or `:OnlyA`. `:OnlyI` is allowed only on
  groups whose steps bind a modifier (a `modifiers:` name or a reactant's `::Mod` copy),
  never on a step that takes up or releases a substrate or product as a reactant.
  Regulatory-site ligands take any of the four tags.
- Phase 2 derives both conformations from the one graph with today's semantics (deletion
  and pruning) and today's guard. Phase 4 brings the limits and admission below; until
  then PR #86's one-flip-segment check stays.

### Limits (B1, B2)

In the conformation a tag excludes:

- a binding's affinity goes to zero: K → ∞, or for a steady-state binding kon → 0 with
  koff fixed;
- a chemistry step is blocked: both rate constants go to zero with their ratio set by the
  cycle rows;
- each member of an `:OnlyA` (or `:OnlyI`) group vanishes at its own rate; the rows still
  bind them along the path. This changes only which tag sets admission accepts, never a
  derived law.

Mixed groups (a plain and a fused binding of one metabolite) stay allowed.

### Admission

Admission replaces `_onlya_haldane_violation`. For each excluding conformation it checks,
exactly and in microseconds:

1. **Existence.** A thermodynamically consistent path to the limit exists: a direction of
   vanishing exponents that satisfies every cycle row with the right signs. The slice's
   graph test (no tagged binding inside a component, no directed cycle of them) equals the
   exact LP on all 2,995 records it was compared on.
2. **Uniqueness.** The limit does not depend on how fast the tagged constants vanish. This
   holds when every form free E cannot reach in that conformation (an island) lies behind
   a tagged binding, on its bound side; the limit is then deletion plus pruning. It must
   hold on every face where data sit (section 8), which rejects cases such as R2 below.
3. **Turnover.** The active conformation still turns the reaction over.

A tag set that fails is rejected: an error naming the step for a hand-written mechanism,
dropped by enumeration. Expected changes against today: AD1 and AD11 rejected; AD14, AD14b
and AD15 admitted; the 81 enumerated `:OnlyA` ping-pong tag sets rejected. R2
(`E + S ⇌ E(S) :: OnlyA`, `E(S) <--> E(P)`, `F + P ⇌ E(P)`, `F <--> E :: OnlyA`) is
rejected: at zero product its inactive limit depends on how fast the isomerization
vanishes.

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

1. **No steady-state turnover.** No catalytic cycle runs through a steady-state step, so
   v ≡ 0. Structural.
2. **RE turnover.** A cycle of rapid-equilibrium steps turns the reaction over, so the rate
   is infinite. Structural; today's "no finite rate" check, ported.
3. **A law that is not the unique fast-RE limit on a face.** The faces are: every product
   zero; every substrate zero; each single reactant zero. Per rapid-equilibrium segment and
   face, a structural check finds:
   - a bottomless segment: no form free of the face's metabolites;
   - a split: the segment's lowest forms are joined only through forms holding the face's
     metabolites;
   - a trap: a form fed on the face whose fast exits all need the face's metabolites. For a
     trap the law equals the limit when the trap-feeding step is the slowest RE step, so
     the error names that step and suggests making it steady-state.

   A support test on the derived law (the denominator does not vanish on any face, with
   modifiers at zero too) backs these up for steady-state dead ends.
4. **`:OnlyA` or `:OnlyI` tags with no limit.** Admission (section 7).

### Fitted

Substrate inhibition and abortive complexes, flux-free steady-state steps, product
activation, and hand-written laws that never saturate. Enumeration keeps not emitting
never-saturating mechanisms (today's `_degenerate` screen).

### Phase 4's first task

Make the face check exact against the RE-speed oracle. The probes found false alarms (6 on
mutants; 16 of 16 trap verdicts in the large census) and 2 misses (non-unique face steady
states). The guard ships only after the check matches the oracle on populations, mutants
and hand-built cases. Then the two specs with traps ("Numerator: redundant SS-bind",
"LDH i-state NonequalAI 6-group") become rejection tests, each with a well-defined variant
that keeps what the spec was written to test.

## 9. Engine, numerator, normalization and printing (phase 5)

### Engine (C1)

- Rapid-equilibrium segments contract to nodes; each form's weight is its segment root's
  weight times the equilibrium ratios along a content-ordered tree.
- The denominator sums, over segments, the segment's weight times its spanning trees of the
  segment graph (King–Altman/Cha).
- Exact polynomial arithmetic: checked `Int64`, retried in `BigInt` on overflow.
- The size cap keeps today's measure (at most 5,000 spanning-tree products); for allosteric
  laws it applies to both conformations together.

### Numerator

The numerator is the rate at which the first substrate is consumed, summed over
steady-state steps, as today: each step contributes its net flux times the copies of that
substrate it effectively moves, including copies carried by rapid-equilibrium bindings on
either side (the RE-weight exponents). It needs no cut and ignores written direction.
After the solve is substituted, the numerator is divided exactly by the Haldane factor
`H = Π[S] − Π[P]/Keq`, so every law prints as `H` times a positive polynomial.

### Normalization: the unit rule

All denominator terms share one unit. Dividing numerator and denominator by a normalizer
with that unit makes every denominator term dimensionless:

1. **Constant term.** If the denominator has one, divide by it. Each term is then an enzyme
   form's abundance relative to free E (Michaelis–Menten: `1 + S/K_m`).
2. **Virtual constant term.** Otherwise solve the units equation: weights `λ_m` over the
   denominator terms with `Σ λ_m = 1` and `Σ λ_m · a_m = 0`, `a_m` being a term's
   concentration-exponent vector; the normalizer is `Π c_m^λ_m`. Cleland's
   `N = Π_X c(ΠS/X) / c(ΠS)^(n−1)` is one solution (ping-pong bi-bi: `K_a·K_b`).
   Fractional weights are allowed and print as roots. To make the choice unique, terms are
   taken in a fixed order (total degree, then fewest non-substrate metabolites, then
   metabolite order); the first affinely independent ones form the basis, and λ is the
   unique solution on that basis.
3. **Ratio form.** If no solution exists, the law depends on some concentrations only
   through their ratio (a rapid-equilibrium Theorell–Chance step `E(A) + B ⇌ E(Q) + P`
   makes every term carry B or P). Divide by the first term in the same order; terms such
   as `(P/B)·c` are dimensionless.

When the normalizer is a polynomial, it prints once as `N = …`. An allosteric law
normalizes each conformation's polynomials by that conformation's own rule.

### Printing (E8, C10)

- `rate_equation_string` prints v in the independent parameters. Each Haldane and
  Wegscheider relation prints as a reference line.
- Next to it, the coefficient form:
  `v = E_total·n_0·(A·B − P·Q/Keq) / (1 + d_A·A + d_B·B + d_A_B·A·B + …)`, then one legend
  line per coefficient in the fitted constants. A coefficient is named by the monomial it
  multiplies.
- Allosteric laws print each conformation's polynomials, then the MWC combination.
- No `1 ^ n` or `L * 1` artifacts.

## 10. Dedup key and identify (phase 5)

### Key (C6)

1. Bring the law to lowest terms.
2. Label each fitted constant by its structural role (colour refinement), so renaming
   constants or metabolites changes nothing.
3. Label concentrations by metabolite name; mark free E.
4. Hash with SHA-256.

Two mechanisms share a key when their laws are identical up to renaming fitted constants.
Laws equal only after a reparametrization (one constant standing for a product of two, for
example) get different keys; the census found such pairs, and they are fitted separately.

### identify

- The key replaces `eq_hash`: identify fits one representative per key and copies the fit
  to every mechanism with that key (today's memo, with a better key). The cross-validation
  candidate list deduplicates by key.
- One-subunit L phantoms (section 13) inherit the parent's fit and are flagged.
- The beam counts each key once, so duplicates cannot crowd out distinct candidates, and
  expands every member of a selected key, so no law becomes unreachable.

## 11. Code generation and speed gates (phase 0 and phase 5)

- `rate_equation` is `@generated` over the compiled types. The generator decodes the Sig,
  solves, derives the law and emits a Horner body (C2), with locals for each dependent, the
  rapid-equilibrium segment sums, the Haldane factor and the MWC terms. Reciprocals and
  roots are computed once. Phase 0 (PR #85) already emits the Horner body; phase 5 adds
  the locals.
- The body allocates nothing and reads inputs by name, so a `NamedTuple` in any field order
  works. `fitted_params`, `parameters` and `metabolites` return constant tuples.
- **Gates (E1):** zero allocations for every tested law; at most 120 ns per call for
  `MECHANISM_TEST_SPECS` and the three largest enumerated bi-bi laws; above that size, a
  bound proportional to the number of denominator terms, its constant calibrated on CI
  hardware in phase 5. The compile cost per fresh mechanism stays within
  `test_compile_budget`'s limits.

## 12. Fitting scale (phase 1)

In relative mode the loss cannot see a common factor on every steady-state rate constant,
so the fitter fixes it after the fit (E9).

- **API.**

  ```julia
  fp  = FittingProblem(m, data; Keq = 2.0)                    # relative (default)
  fp  = FittingProblem(m, data; Keq = 2.0, absolute = true)   # per-enzyme turnover data
  fit = fit_rate_equation(fp, optimizer)

  kcats = [(A = 10.0, B = 5.0) => 42.0,                       # conditions => turnover
           (A = 20.0, B = 5.0) => 44.0]
  p = rescale_parameter_values(m, fit.params, kcats)
  ```

  identify takes `absolute` and no kcat keyword. `scale_k_to_kcat` and `_kcat_forward`
  are removed (a breaking change).
- **Default scale.** A relative fit returns constants scaled so the mechanism's peak forward
  rate per enzyme over all concentrations, products included, is 1. A deterministic
  multi-start search in log concentrations, each metabolite from 1e-12 to 1e12 in the
  data's units, box corners included, finds the peak; it underestimates a saturating peak
  by about 1e-12 relative.
- **Unscaled fits.** A law still rising at the box edge (one that never saturates) is
  returned unscaled and flagged. The scale never turns a fit into a failure.
- **Rescaling.** `rescale_parameter_values(m, params, kcats)` takes one or more measured
  turnovers (kcat or specific activity in turnover units) at stated concentrations;
  unlisted metabolites are zero. Several anchors combine by the geometric mean of measured
  over predicted turnover, the log-space least-squares scale; the result reports their
  spread. An anchor at which the model predicts zero or opposite-signed turnover is an
  error. Without `kcats`, the function applies the default scale.
- **Absolute mode** is unchanged.
- **Docstring fix.** `maxtime` applies to each restart, not to the whole fit.

## 13. Enumeration (phases 4 and 5)

- The moves (rapid-equilibrium to steady-state flip, split, dead-end, allosteric tagging)
  are ported to the graph `Mechanism`. Their children match today's except where known
  issues 1 and 2 change today's counts.
- Every child passes the constructor, the well-posedness guard and, if allosteric,
  admission before identify sees it. Today's `_degenerate` screen stays.
- The tagging move emits no `:OnlyI` tags.
- **Islands** are rejected by admission (phase 4).
- **Duplicate tag sets** (siblings that differ only in the tags of groups the inactive
  conformation never reaches) stay in the frontier and share one fit through the key
  (phase 5).
- **One-subunit L phantoms** stay in the frontier, inherit the parent's fit and are
  flagged (phase 5). A variant is an L phantom when its inactive conformation binds nothing
  and does not catalyze, and the identifiable rank shows L does not change the law; its law
  is then the parent's, with constants rescaled. A variant whose inactive conformation binds
  something is fitted.

## 14. Testing

- **Oracles move into `test/`:**
  - the exact mass-action steady state, plain and two-conformation, in rational or
    BigFloat arithmetic;
  - the RE-speed oracle: the law is the limit of the mass-action steady state as RE speeds
    grow, at interior points and on faces, under several speed patterns;
  - the ε-path oracle: moving `:OnlyA` or `:OnlyI` constants along the certified direction
    converges to the law;
  - detailed balance at equilibrium concentrations.
- **Populations** are regenerated inside the tests from `init_mechanisms` and
  `expand_mechanisms`, sampled to keep the suite within its time budget, never loaded from
  scratch files. Over them the tests check oracle agreement, determinism under rewrites,
  regroupings and renames, and key stability.
- **Textbook bridge (C5)** stays in the tests, rewritten by role: each test reads its
  written steps and takes the package constant of the matching stored step in the written
  direction, through `all_parameters`. The DSL's test-only source-order feature goes.
- **Enumeration tests** keep CLAUDE.md's rules: inline fixtures and exact children.
- **Pins.** Each phase updates the pins of the behavior it changes, in the same PR. Before
  changing a test that checks only a string, the phase adds an oracle value check.
- **Gates.** Section 11's speed, allocation and compile-budget tests.

## 15. Phases

Each phase is one PR, green with the full suite, proposing a version bump for Denis to
decide (breaking changes are acceptable before 1.0). Before phase 1: PRs #85 (Horner code
generation), #86 (one flip segment, ROADMAP bug list) and #87 (`:OnlyA` guard crash) merge,
and the fix for the ter-bi expansion abort (one oversized split candidate aborts a parent's
whole expansion) lands as its own PR.

| Phase | Contents | Gate or dependency |
|---|---|---|
| 0 | Horner code generation (PR #85) | open |
| 1 | Fitting scale (section 12) | none |
| 2 | Joint solve, names, parameter API, single-graph allosteric derivation, on today's types (sections 6, 7) | none |
| 3 | Graph `Mechanism`, compiled types and Sig, DSL, `show`, parser, CSV, constructor checks, free E (sections 4, 5) | GC-crash probe (E6) |
| 4 | Admission, flip at E, `:OnlyI` rule, well-posedness guard, enumeration checks (sections 7, 8, 13) | face check exact first |
| 5 | Engine, numerator, normalization, printing, key, identify's key cache and beam, code-generation locals, scaled speed gate (sections 9, 10, 11, 13) | none |

Each phase updates the docs and ROADMAP entries it affects (section 16).

## 16. Docs and ROADMAP

- **Rewritten pages:** the DSL reference (two macros, `modifiers:`, `::Mod`,
  `regulatory_site` lines, free E, the `:OnlyI` rule); `deriving/mwc_allostery.md` (flip at
  E, limits, admission); `deriving/textbooks.md` (Full-mode section removed, normalization);
  `fitting/normalized_vs_absolute.md` (default scale, rescaling after the fit, absolute
  mode); the identify pages (CSV format, key, inherited fits); `developer.md`
  (architecture).
- **ROADMAP:** each phase retires the entries it fixes, as CLAUDE.md requires. New entries
  until fixed: the `maxtime` docstring (phase 1) and the ter-bi expansion abort (pre-phase
  PR).

## 17. Risks and how each is handled

| Risk | Handling |
|---|---|
| The slice's unexplained "GC error (probable corruption)" in a workload that deserializes mechanisms | A probe of at most half a day before phase 3 (E6): rebuild instead of deserializing, run on Julia 1.10 and on x86, report upstream if Julia is at fault |
| The face check rejects valid mechanisms or misses invalid ones | Phase 4 starts by making it exact against the RE-speed oracle |
| The L-phantom detector costs or errs | Phase 5 measures it; oracle tests compare flagged variants with the parent's law |
| The peak search misses the peak | Deterministic starts, box corners, tests against laws with known peaks (substrate inhibition, product-assisted release) |
| Test-suite time grows | Population tests are sampled to a stated budget |
| Probe artifacts vanish with the session | This spec records every decision; probe code is reference only |

## Appendix A: Decision index

| Card | Decision | Section |
|---|---|---|
| A1 | Proceed; the derivation stays about the same size | 1 |
| A2, E11 | Phased PRs in the critic's order; the tests-only phase was dropped in brainstorming | 15 |
| A3 | Vertical-slice probe run | Appendix B |
| B1 | `:OnlyA` bindings by affinity, chemistry by kinetic block; island rule; cone test | 7 |
| B2 | Per-step classes in the excluding conformation; mixed groups allowed | 6, 7 |
| B3 | Superseded by E10 | 7 |
| B4 | Superseded by E9 | 12 |
| B5 | No all-parameters law; `all_parameters` replaces the test helper | 6 |
| B6 | Rational exponents accepted | 6, 9 |
| B7 | Tie forward with forward | 6 |
| B8 | Flux-free steady-state steps allowed | 6, 8 |
| B9 | Identifiable rank is an analysis tool; one enumeration use (L phantoms) | 2, 13 |
| B10 | `:EqualAI` regulator constants bare | 6 |
| C1 | Lean engine | 9 |
| C2 | Horner with segment, Haldane and MWC locals | 11 |
| C3 | Mechanism as the only type parameter; derive in the generator | 4, 11 |
| C4 | Constant term, else virtual constant term; refined into the unit rule | 9 |
| C5 | Textbook bridge in the tests, by role | 14 |
| C6 | Rename-proof key; concentrations labelled by name | 10 |
| C7 | No source-level duplicate suppression | 2, 13 |
| C8 | Priority ties broken by content names | 6 |
| C9 | Encoding limits with clear errors; orientation ties by form index | 4 |
| C10 | Coefficient form named by monomial | 9 |
| D1 | Cheap fixes first (PRs #85, #86, #87) | 15 |
| E1 | 120 ns gate plus a size-scaled bound | 11 |
| E2 | Fit what is well defined, reject what is not | 8 |
| E3 | Names after the kept member | 6 |
| E4 | Alphabetical parameter order | 6 |
| E5 | Reserved names; no `_` in metabolite names | 4 |
| E6 | GC-crash probe before phase 3 | 17 |
| E7 | "Against" means against the reaction direction | 6 |
| E8 | Printed laws in independent parameters | 9 |
| E9 | Peak scale over all concentrations; measured kcats applied after the fit | 12 |
| E10 | Flip at free E; free E mandatory in every mechanism | 4, 7 |

Brainstorming decisions without a card: two macros with one vocabulary and composed types;
`modifiers:` and `::Mod`; `:OnlyI` limited to modifier bindings; `parameters` and
`all_parameters` exported, `fitted_params` internal; the unit rule; kcat rescaling after the
fit; duplicates and L phantoms kept with inherited fits; the beam counts keys once.

## Appendix B: Evidence

From the vertical slice and the verification probes, in the session scratchpad (reference
only):

- **Slice.** 237,965 laws from 257,698 mechanisms (all 250,855 ter-ter seeds included)
  matched their oracles; worst error 7e-31 on the core sets. Dependents identical across
  24,580 rewritten inputs and 12,608 same-model regroupings. Admission equalled the exact LP
  on 2,995 of 2,995 allosteric records. Code size about 2,045 → 1,980 lines; derivation
  926 → 1,110. The largest bi-bi laws ran in about 42 ns; 343 ter-ter laws ran over 120 ns.
- **kcat = 0 laws.** All 284 are complete substrate inhibition, exactly: positive rate at
  every finite substrate concentration, falling to 0 only as one substrate grows.
- **Peak search.** The textbook steady-state random-order bi-bi never exceeds its
  both-saturated rate with no product present (proved two ways); with a product present it
  can (Q rescuing a slow E(P) release, 125.6x in a simple example).
- **Well-posedness.** Among 13,645 records (8,626 RE/SS mutants): 1,654 with RE turnover
  (the slice's law v ≡ 0 or finite and wrong); 1,027 accepted laws with v ≡ 0; 22 NaN on a
  face; face splits and traps with gaps of 2.3% median, 47% at the 90th percentile and
  314% at most, and up to 80% with one product at zero. Enumerated output never produced
  v ≡ 0, RE turnover or a 0/0 face.
- **Flip at E.** Exact against an independent oracle (12 shapes, 1 to 3 subunits, worst
  3.2e-69); no change for 9,308 enumerated mechanisms.
- **Census.** Islands take 8–32% of a round's fit time in ping-pong and the ter reactions;
  duplicate tag sets 9–21%; one-subunit L phantoms 10–22% in bi-bi and ping-pong.

## Appendix C: Glossary

- **Face.** A set of concentration points where some metabolites are exactly zero, such as
  every product zero (initial-rate data).
- **Flip form.** The form at which the two conformations interconvert: free E.
- **Island.** A form free E cannot reach in a conformation once that conformation's
  excluded steps are removed.
- **Kinetic group, class.** Steps that share constants; a class is one shared constant in
  the joint solve.
- **L phantom.** An allosteric variant in which L enters the law only through a combination
  other constants absorb, so the data cannot determine it.
- **Anchor.** A measured turnover at stated concentrations, used to rescale a relative fit.
- **Peak turnover.** The largest forward rate per enzyme over all concentrations.
