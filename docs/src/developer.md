# Developer / Architecture

This page is for maintainers and contributors. It documents the internal
architecture of EnzymeRates.jl, corrected against the source as it stands on
`main`. Public API lives on the [API Reference](@ref) page.

---

## Concrete vs. singleton mechanism types

The package uses two parallel mechanism representations.

**Concrete types** (`Mechanism`, `AllostericMechanism`) are the working
representation for the enumeration pipeline. `Mechanism` (`src/types.jl`)
has exactly two fields:

- `reaction::EnzymeReaction` — the reaction the mechanism catalyzes.
- `steps::Vector{Vector{Step}}` — kinetic groups. The outer vector is one
  entry per group; each inner vector holds the `Step`s that share that
  group's kinetic parameters.

`Step` has four fields: `from_species::Species`, `to_species::Species`,
`bound_metabolite::Union{Metabolite,Nothing}`, and `is_equilibrium::Bool`.
The entire enumeration pipeline — `init_mechanisms`, `expand_mechanisms`, and
`_dedup_flat!` — operates end to end on these concrete structs. There is no
separate intermediate representation.

**Singleton types** (`EnzymeMechanism{Sig}`, `AllostericEnzymeMechanism{...}`)
are used only by the `@generated` rate-equation derivation. A singleton type
carries all mechanism data as a Julia type parameter, which lets the compiler
specialize the derivation at compile time.

`EnzymeMechanism{Sig}` (`src/types.jl`) encodes one `Mechanism` as the tuple
`(reaction_sig, steps_sig)` produced by `_sig_of(::Mechanism)`. The `Sig` is
purely structural: two mechanisms that differ only in the order their steps
were written collapse to the same `EnzymeMechanism` type, because the
`Mechanism` constructor canonicalizes step and group order before `_sig_of`
runs.

The boundary converters are:

- `EnzymeMechanism(m::Mechanism)` — lifts to the singleton type. It first
  calls `_drop_unbound_regulators(m)` to remove any regulators declared on the
  reaction but not bound by any step (e.g., a dead-end inhibitor before any
  expansion move binds it), so they neither appear in `regulators(em)` nor
  contribute a parameter.
- `Mechanism(em::EnzymeMechanism{Sig})` — lifts back via `_mechanism_from_sig`
  (`src/types.jl`).
- `compile_mechanism` (internal, not exported) is the single entry point for
  both families: `compile_mechanism(m::Mechanism) = EnzymeMechanism(m)` and
  `compile_mechanism(am::AllostericMechanism) = AllostericEnzymeMechanism(am)`
  (`src/mechanism_enumeration.jl`). It is accessible as
  `EnzymeRates.compile_mechanism`.

### Three type parameters for `AllostericEnzymeMechanism`

`AllostericEnzymeMechanism{CatalyticMech, CatSites, RegSites}` (`src/types.jl`)
uses three type parameters rather than folding everything into one value-tuple
`Sig` as `EnzymeMechanism{Sig}` does. The reason is a Julia constraint: the
first slot is a `DataType` (a concrete `EnzymeMechanism` subtype), and Julia
rejects a `DataType` in the value-tuple position of a type parameter. The three
slots are therefore kept separate.

---

## The `name(p, m)` parameter-naming chokepoint

All `Parameter → Symbol` rendering in `src/` flows through one function family:
`name(p::Parameter, m)` at `src/types.jl`. Routing all parameter-name
production through one place means any name-scheme change is a single-function
edit.

The `Parameter` subtypes divide into three groups:

- **Step-governed** (`Kd`, `Kiso`, `Kon`, `Koff`, `Kfor`, `Krev`) — carry a
  `step::Step` and a `state::Symbol`.
- **Regulator-site governed** (`Kreg`) — carries a `RegulatorySite`, an
  `AllostericRegulator`, and a `state::Symbol`.
- **Stateless scalars** (`Keq`, `Etot`, `Lallo`) — carry no step or state.

Three private helpers do the rendering:

- `_state_tag(state)` maps the allosteric state token to a name prefix: `:A →
  "A_"`, `:I → "I_"`, `:EqualAI → ""`, `:None → ""`.
- `_render_binding(prefix, rep, state)` names a binding parameter as metabolite
  + pre-binding form — for example, binding K for substrate S at free enzyme E
  renders as `:K_S_E`; the inactive-state binding K renders as `:K_I_S_E`.
- `_render_iso(prefix, from, to, state)` names an iso parameter by directed
  species pair — for example, `:k_ES_to_EP`.

### AST-walker enforcement

The chokepoint is enforced by an AST-walker test in `test/test_types.jl`
(the `"chokepoint: no Symbol(\"[KkVL]...\") outside parameter-name renderers"`
testset). The walker uses `_walk_violations!`, which parses each `src/*.jl`
file with `Meta.parseall`, identifies a chokepoint method body via
`_is_chokepoint_def` (a `name` method dispatching on a `Parameter` subtype
value), and fails the build if any `Symbol("K…")`/`Symbol("k…")`/
`Symbol("V…")`/`Symbol("L…")` literal appears outside a chokepoint body.
This is an AST walk, not a regex search; string-building patterns such as
`Symbol("K_", …)` are caught correctly.

---

## Canonical Step Form

Canonicalization is load-bearing, not cosmetic. The Haldane/Wegscheider
reduction (in `src/thermodynamic_constr_for_rate_eq_derivation.jl`) chooses
which parameters are dependent by step order. Without canonical step order, the
reduced rate equation and `fitted_params` would differ between two mechanically
identical mechanisms written in different step orderings.

**Binding steps** (RE and SS) are canonicalized in the `Step` constructor
(`src/types.jl`): the bound metabolite is always on `to_species`, so the free
enzyme and free metabolite are on the from-side. `E + S ⇌ ES` stores as
written; `ES ⇌ E + S` stores reversed. Product-release `EP ⇌ E + P` therefore
stores as `E + P → EP` — a product-release step and a product-binding step are
the same representation.

**Iso steps** (RE and SS) are canonicalized to the physical-forward direction
by `_canonical_iso_direction` (`src/types.jl`), called for every group by
`_canonicalize_iso_groups` inside the `Mechanism`/`AllostericMechanism`
constructors. The function applies three tiers:

1. **Tier 1 (atom-balance)** — the form with more substrates bound (or fewer
   products bound) is `from`. This suffices for most mechanisms.
2. **Tier 2 (binding-graph context)** — when Tier 1 ties, `_entry_kind`
   (`src/types.jl`) classifies each endpoint by whether it appears as the free
   side of substrate-binding steps, product-binding steps, both, or neither.
   The function examines ALL binding steps (RE and SS), not just RE: the
   substrate-entry / product-exit property is a chemistry fact independent of
   whether the binding step is rapid-equilibrium. A product-only form → a
   substrate-only form is treated as forward.
3. **Tier 3 (lexical fallback)** — lexicographic comparison of form names,
   so the result is always deterministic.

**Step and group order** is canonicalized by `_canonical_group_order!`
(`src/types.jl`) using `_step_canonical_key`. Steps within a group are sorted
first, and then groups are sorted by the canonical key of each group's first
step. Because the constructor always canonicalizes, `_dedup_flat!` (which
`unique!` implements) is a pure non-mutating `==`/`hash` comparison with no
mutation needed.

---

## The `@generated` compile-budget Known Issue

`rate_equation(m, conc, params)` and the supporting functions use `@generated`
functions that derive the rate equation at compile time by the King–Altman/Cha
method. Each unique `EnzymeMechanism` type triggers a full symbolic derivation.
For mechanisms with many enzyme forms or steps, this can be very slow, exhaust
memory, or `StackOverflow`. `identify_rate_equation` caps expansion via
`max_param_count`, bounding search depth, but not per-mechanism compile cost.
One very large mechanism can therefore be slow even with a low `max_param_count`.

Once compiled, `rate_equation` must be allocation-free and sub-100 ns per call.
This constraint is enforced by `test_rate_equation_performance` in
`test/test_rate_eq_derivation.jl` (asserts `allocs == 0` and `t < 100e-9`
for every mechanism in the shared test fixtures). The fitter calls
`rate_equation` millions of times per cross-validation fold; any per-call
allocation or microsecond-scale overhead makes the package unusable in
practice. If a change would force `rate_equation` to allocate or slow down,
stop and discuss before implementing.

`loss!` is held to the same standard. `FittingProblem` pre-allocates the
`log_ratios_buffer` once at construction, and `loss!` reuses it on every call,
so the fitter's inner loop allocates nothing. Combined with the `rate_equation`
contract above, this is what makes multi-start fitting over millions of loss
evaluations practical.

---

## Source layout

A file-by-file map for maintainers.

- **`src/types.jl`** — Concrete structs `EnzymeReaction`, `Mechanism`,
  `AllostericMechanism`, `Step`, `Species`, `Residual`, `ReactantAtoms`,
  `RegulatorMults`, `RegulatorySite`, `Metabolite` family (`Substrate`,
  `Product`, `CompetitiveInhibitor`, `AllostericRegulator`), `Parameter` family
  (`Kd`/`Kiso`/`Kon`/`Koff`/`Kfor`/`Krev`/`Kreg`/`Keq`/`Etot`/`Lallo`);
  singleton derivation types `EnzymeMechanism{Sig}` and
  `AllostericEnzymeMechanism{...}` with `_sig_of` / `_mechanism_from_sig`
  converters; the `EnzymeMechanism(m::Mechanism)` lift and
  `AllostericEnzymeMechanism(cm, cat_sites, reg_sites)` constructor with
  allosteric-state validation; struct accessors; the `name(p, m)` chokepoint;
  `RateEquationMode` hierarchy.

- **`src/dsl.jl`** — `@enzyme_reaction` (supports atom brackets);
  `@enzyme_mechanism` and `@allosteric_mechanism` accepting the
  decomposed-`Species` call grammar `E(S)` / `E(A, B)`. A ligand or free
  metabolite may carry a `::Inh` role tag (`G6P::Inh`) to bind in its
  `CompetitiveInhibitor` role. The opaque bound-form bare-`Symbol` grammar
  (`:ES`, `:E_S`) is rejected at parse with a migration error.

- **`src/sym_poly_for_rate_eq_derivation.jl`** — Laurent symbolic polynomial
  algebra (`POLY`, `MONO` with possibly-negative exponents, `_poly_to_expr`
  rendering negative exponents as denominators); `_rename_symbols`,
  `_zero_symbols_in_poly` for MWC allosteric-state-driven substitution;
  `_reduce_conc_lowest_terms` (the concentration-GCD: shifts each
  concentration's exponent so its minimum across numerator and denominator is
  zero, so no concentration ever sits in a denominator).

- **`src/rate_eq_derivation.jl`** — King–Altman/Cha rate equation derivation
  via `@generated` functions. Each rapid-equilibrium segment is referenced to
  its free enzyme, and the Cha α/numerator/denominator are computed directly as
  fractional/Laurent `POLY`s (no common-denominator linearization), then
  `_reduce_conc_lowest_terms` clears any residual concentration denominator.
  The parameters API; kcat computation (`_ss_rate_constant_names`,
  `_kcat_groups_from_polys`); `rescale_parameter_values`;
  `AllostericEnzymeMechanism` MWC rate equation assembly
  (`_build_allosteric_rate_body`, `rate_equation_string`); helpers for
  allosteric symbol selection, renaming, and dependent-parameter assignments.
  Parameter-symbol rendering goes through `name(p, m)` — no direct
  `Symbol("K…")` literals.

- **`src/thermodynamic_constr_for_rate_eq_derivation.jl`** — Haldane/Wegscheider
  thermodynamic constraints; `_dependent_param_exprs` builds the kinetic-group
  merge map up front and applies column merging before Gaussian elimination.

- **`src/fitting.jl`** — `FittingProblem`, `loss!`, and `fit_rate_equation`
  using Optimization.jl. `FittingProblem` carries
  `scale_k_to_kcat::Union{Real,Nothing}`: a `Real` selects relative data
  (per-group-centered loss) plus rescaling of SS k's to that kcat; `nothing`
  selects absolute per-enzyme turnover (uncentered loss, no rescale).
  `fit_rate_equation` returns `(params, loss, retcode)`, where `retcode` is
  the best restart's `Symbol(sol.retcode)`.

- **`src/identify_rate_equation.jl`** — `IdentifyRateEquationProblem` and
  `identify_rate_equation`, which run the advancing-target beam search over
  actual fitted-parameter count with leave-one-group-out cross-validation.
  `_process_batch` fuses compile + cap + fit per worker into `BatchEntry`s;
  `_ingest!` maintains the frontier and bounded CV pool; equation dedup uses
  the comment-stripped `_rate_eq_dedup_key`. `save_dir` is mandatory, writing
  `initial_mechanisms.csv` and `equation_search_iteration_N.csv`. A
  `FitFailure` captures the exception text of any compile/fit that throws.

- **`src/mechanism_enumeration.jl`** — Building blocks: `init_mechanisms`,
  `expand_mechanisms`, and `_dedup_flat!`. Native expansion moves:
  `_expand_re_to_ss`, `_expand_split_kinetic_group`,
  `_expand_add_dead_end_regulator`, `_expand_to_allosteric`,
  `_expand_add_allosteric_regulator`, `_expand_change_allo_state` — each
  dispatches on `Mechanism` / `AllostericMechanism`. The pipeline builds
  decomposed `Mechanism` / `AllostericMechanism` directly from `Step` /
  `Species`; there is no intermediate working representation.

---

## Internal API

```@docs
EnzymeRates.Mechanism
EnzymeRates.Step
EnzymeRates.init_mechanisms
EnzymeRates.compile_mechanism
```
