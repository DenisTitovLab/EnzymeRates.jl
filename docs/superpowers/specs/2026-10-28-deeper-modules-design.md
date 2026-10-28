# Deeper modules: design

Branch `deep-modules-simplification` at 05315f78, on top of the first module split
(`2026-10-27-deep-modules-simplification-design.md`). Same PR, version 0.10.0 (Denis: PRs #85–#95 and this PR ship together as 0.10.0).

## Intent

The first split gave the package four modules and a boundary test. Their interfaces are still
wide: `Mechanisms` exposes 78 names, `Derivation` 11, `Fitting` 5. Most of the `Mechanisms`
names exist because one rule (outside `Mechanisms`, never read a mechanism's fields) meets a
storage that cannot be read without decoding: metabolites packed four bits each into a
`UInt128`, indexed by a per-mechanism table.

This work makes each module deeper. The mechanism becomes plain data that any module reads
with dot notation, functions remain only where they compute something, and each function
lives in the module whose job it is. The goal is an interface someone can learn in minutes,
written up in the developer docs with a diagram.

## Goals

- **Mechanisms is plain data.** A mechanism is a list of `Node`s (enzyme forms) and `Step`s
  (mechanism steps) with public, documented fields. The written form (`Raw*`) and the packed
  storage (`Counts`, `Form`, `Edge`, `Metabolites`) are deleted. About 78 names become about
  21: 12 types and 9 functions (§3).
- **Each function lives where its job is.** Allosteric law semantics, the turnover-cycle test
  and the segment walks move to `Derivation`; the parameter naming moves to `Derivation`; the
  DSL reader and the macros move to a new `DSL` module; the rate table moves to a new `Data`
  module (§2).
- **Smaller Derivation and Fitting interfaces**, with related functions merged (§4, §5).
- **A developer-docs page section** that explains the module structure with a diagram and a
  one-line description of every interface name (§8).

## Non-goals

- **No change to `rate_equation`.** Its generated code, speed and zero allocations stay.
- **No change to what enumeration produces.** Child sets, laws and parameter names stay as
  they are; the census proves it.
- **No change to printed text.** `show` of every mechanism and reaction stays byte-identical,
  so canonical order is today's, saved identify results still parse, and
  `test/fixtures/canonical_text_golden.txt` stays as it is. (Claude's call: N2 must
  reproduce today's `isless(::Form)` order anyway for the parameter names.)
- **No new features** beyond the default reference scale and `test_rate_equation` (§5).

## Anti-goals

Any one of these makes the work a failure, even with green tests:
- A derived law, parameter name or enumeration child set changes.
- `rate_equation` allocates or slows down (CLAUDE.md's bound).
- Enumeration, the guard or identify gets slower than the gate Denis sets (§7) without his
  approval.
- A hand-written mechanism the guard rejects loses today's message or fix suggestions.
- A module reaches a sibling's internals.

## Constraints

- CLAUDE.md applies in full: TDD, smallest change, no rewrite without permission (this
  document is Denis's permission for the rewrites it names), pristine output, `_testhelper_`
  prefixes, ABOUTME headers, one Julia process at a time.
- Julia 1.10 through 1.13. API lists stay plain constants, not `public` declarations.
- Canonical form stays load-bearing: deduplication by `==`/`hash` and every parameter name
  rely on it (CLAUDE.md, "Canonical Step Form").

## Approaches

**Mechanism storage.**
- *N1, read-only views over today's packed storage.* `nodes(m)` and `steps(m)` return lazy
  `Node`/`Step` records; the `Raw*` write side stays. No performance risk, but two models
  remain (read and write) and the packing stays.
- *N2, one model (chosen).* `Node` and `Step` hold metabolites as sorted `Vector{Met}`; the
  constructor takes them in any order and canonicalizes. The spike (§7) measured hash 0.87×,
  reads 0.83×, build 1.46×, node sort 1.55×, 5.6× more allocations per mechanism, and
  estimated enumeration at 3–10% slower.

**Reading the model.** Accessor functions per fact (today), or dot notation on public fields
(chosen). Accessors inline to the same code and protect nothing (`steps(m)` and `m.steps` hand
out the same vector), and Julia's `getproperty` can later turn a field into a computed
property, so fields lock nothing in. Rule: stored data is read with dots; derived data is a
function.

**Text.** A separate printing module was rejected: the `Mechanisms` constructors build their
error messages from `form_text`/`step_text`, so printing must sit in or below `Mechanisms`,
and Julia keeps `show` with its type. A reaction DSL module apart from the mechanism DSL was
rejected as too shallow.

## 1. Layering

```
              Data ──────────────┐
                                 ▼
Mechanisms ◄── Derivation ◄── Fitting
                    ▲
                    ├──────── Enumeration
                    └──────── DSL
top level (identify, exports) uses all six
```

A module names only the modules its arrows reach: `Mechanisms` and `Data` name none,
`Derivation` names `Mechanisms`, and `Fitting`, `Enumeration` and `DSL` name `Mechanisms` and
`Derivation`; `Fitting` also names `Data`. The boundary test's layering rule encodes this.

## 2. What moves

| From | To | What |
|---|---|---|
| `Mechanisms` | `Derivation` | `conformation_steps`, `excludes`, `exclusion_components` (private), `has_turnover_cycle`, `has_unbalanced_cycle`, `re_segments`, `segment_graph`, `segment_walk` (private), `present_species` (private), `runs_toward_products` (private) |
| `Mechanisms` | `Derivation` | the parameter naming: `constant_columns`, `regulator_constant`, the `name(p, m)` chokepoint, `_priority_key`, `_ordered_columns`, `_kept_members`, all private there; the chokepoint test moves to `test/Derivation/`. The form-name rendering they build on stays in `Mechanisms`, whose constructor checks need it (distinct form names, step messages), and is reached through `step_sides` (§3) |
| `Mechanisms` + `macros.jl` | `DSL` | the reader (`raw_mechanism`, `reaction_expr`, `check_uniquely_decodable`), the three macros and `_parse_mechanism` |
| `Fitting` | `Data` | the rate-table and reference validation (`rate_data`, `_rate_table`, `_reference`, `_check_anchors`) |
| `Fitting` | top level | `cross_validation_loss`, private to identify, rebuilt on `fit_rate_equation` and `test_rate_equation` |
| `Mechanisms` | top level | `mechanism_id`, private to identify |

The printer (`show`, `form_text`, `step_text`) stays in `Mechanisms`.

## 3. Mechanisms

**Types and their public fields.**

| Type | Fields |
|---|---|
| `Met` | `name::Symbol`, `modifier::Bool` |
| `Node` | `conformation::Symbol`, `bound::Vector{Met}`, `added::Vector{Symbol}`, `subtracted::Vector{Symbol}` (the residual) |
| `Step` | `from::Int`, `to::Int` (node indices), `taken_up::Vector{Met}`, `given_off::Vector{Met}`, `rapid_equilibrium::Bool`, `group::Int` |
| `Mechanism` | `substrates::Vector{Symbol}`, `products::Vector{Symbol}`, `nodes::Vector{Node}`, `steps::Vector{Step}` |
| `RegulatorySite` | `ligands::Vector{Symbol}`, `multiplicity::Int`, `mwc_tags::Vector{Symbol}` |
| `AllostericMechanism` | `catalytic::Mechanism`, `mwc_tags::Vector{Symbol}` (one per kinetic group), `catalytic_multiplicity::Int`, `regulatory_sites::Vector{RegulatorySite}` |
| `Reactant` | `name::Symbol`, `atoms` (element => count) |
| `Regulator` | `name::Symbol`, `kind::Symbol` (`:allosteric` or `:competitive`), `effect::Symbol` (`:activator`, `:inhibitor` or `:unspecified`; today's `reg_type`, which the enumeration filter `_respects_reg_type` and `show` read) |
| `EnzymeReaction` | `substrates::Vector{Reactant}`, `products::Vector{Reactant}`, `regulators::Vector{Regulator}`, `allowed_catalytic_multiplicities::Vector{Int}`, `shared_catalytic_site::Vector{Tuple{Symbol,Symbol}}` |
| `EnzymeMechanism{Sig}`, `AllostericEnzymeMechanism{Sig}`, `AbstractEnzymeMechanism` | unchanged |

`mwc_tags` replaces `cat_allo_states` and `allo_states`; the values stay `:EqualAI`,
`:NonequalAI`, `:OnlyA`, `:OnlyI`. `Reactant` and `Regulator` replace `Substrate`, `Product`,
`ReactantAtoms`, `DeclaredRegulator`, `AllostericRegulator` and `CompetitiveInhibitor`.

**Functions** (each computes something):

| Function | What it does |
|---|---|
| `step_kind(m, k)` | `:isomerization`, `:binding`, `:fused_binding` or `:chemistry`; replaces `step_kind` and `is_binding` (`is_binding` ≡ `:binding` or `:fused_binding`) |
| `groups(m)` | the steps of each kinetic group |
| `metabolites(m)` | the concentration names in rate-law order; also on the compiled types; exported |
| `modifiers(m)` | every name bound as a modifier |
| `catalytic(m)` | `am.catalytic`, or `m` itself for a plain mechanism |
| `compile_mechanism(m)` | lifts to `EnzymeMechanism{Sig}` / `AllostericEnzymeMechanism{Sig}` |
| `form_text(m, x)`, `step_text(m, k)` | the DSL text of a node or a step |
| `step_sides(m, k)` | the two side names every constant of step `k` is named after (`("EA_B", "EQ_P")`); today's `_sides` |

That is 9 functions. The "declared modifier that no step binds" check moves from the
constructor to `DSL`, with its message unchanged, since the constructor no longer takes a
modifier list.

**Construction.** `Mechanism(substrates, products, nodes, steps)` takes nodes in any order,
duplicates allowed, and steps in any direction, with `group` as labels 1, 2, …; it
canonicalizes (deduplicates nodes, orients and orders steps, renumbers groups) and runs
today's representation checks. `AllostericMechanism` has two constructors: from the written
parts with `mwc_tags[label]` per written group label (what `RawAllostericMechanism` does
today), and the composition `AllostericMechanism(catalytic, mwc_tags, multiplicity, sites)`
over a canonical catalytic mechanism. Canonicalization is private: `canonical_graph`, the
writers (`with_flags`, `regroup`, `with_tags`) and the `Raw*(m, …)` converters are deleted.
A caller rebuilds a mechanism by copying `m.steps` with changed fields and calling the
constructor. The split move's group-to-tag mapping (today inside `regroup`) becomes an
Enumeration helper. If the dead-end move and the seed bases cannot use the checked
constructor (a check fails on their intermediate graphs, or the checks cost more than the
gate allows), the builder keeps an unchecked constructor private to `Mechanisms` and says why
in the plan.

**Sig.** The encoding follows the new storage; its layout stays internal.
`compile_mechanism(Mechanism(em)) === em` keeps holding.

## 4. Derivation

| Name | Change |
|---|---|
| `parameters(m)` | ≡ `parameters(m, IndependentParameters())`: fitted constants, `Keq`, `E_total` (what `rate_equation` takes) |
| `parameters(m, FittedParameters())` | replaces `fitted_params`; stays `@generated` on the compiled types so `loss!` stays allocation-free |
| `parameters(m, AllParameters())` | the names above plus the dependent constants |
| `parameters(m, values::NamedTuple)` | replaces `all_parameters`: `values` plus every dependent constant's value |
| `check_valid_rate_equation(m; throw = true)::Bool` | replaces `has_valid_rate_equation` (`throw = false`) and today's throwing check |
| `rate_equation_dedup_key(m)` | renames `rate_equation_key`; its capped fallback computes the same text id itself (`bytes2hex(sha256(text))[1:16]` of the printed mechanism), so key strings do not change |
| `parameter_count(m)` | folds into `length(parameters(m, FittedParameters()))` if Enumeration's call sites (`moves.jl:395–400`, which count before checking validity) allow it; otherwise stays |
| `rate_equation`, `rate_equation_string`, `rate_equation_complexity`, `scale_rate_constants` | unchanged |
| moved in (§2) | `conformation_steps`, `excludes(tag, c)`, `has_turnover_cycle`, `has_unbalanced_cycle`, `re_segments`, `segment_graph` |

`IndependentParameters`, `FittedParameters` and `AllParameters` are exported.

## 5. Data and Fitting

**Data.** `RateData(table, metabolite_names; Keq, absolute = false, reference = nothing)`
holds the validated column table, `Keq`, `absolute` and `reference`, with today's checks.
`FittingProblem` and `IdentifyRateEquationProblem` each hold one instead of repeating those
fields.

**Default reference.** In relative mode `reference` becomes optional. Without it,
`fit_rate_equation` scales the steady-state rate constants (`scale_rate_constants`) so the
largest |predicted rate| over the data's concentration rows, at `E_total = 1`, equals 1. With
it, behaviour is today's. Absolute mode is unchanged; `reference` with `absolute = true` still
errors. The docs say plainly that default-scaled rate constants are in units of the fastest
measured condition, not kcats.

**`test_rate_equation(m, params, table; Keq, absolute = false) -> Float64`**: the fitting loss
of `m` at `params` on `table`, exported. Identify's cross-validation becomes a private
"fit on the other groups, test on the held-out one".

**Deleted:** `rescale_parameter_values` (and its tests), `cross_validation_loss`,
`rate_data`. Fitting's API: `FittingProblem`, `fit_rate_equation`, `test_rate_equation`.

## 6. DSL

The `DSL` module holds the reader for `@enzyme_reaction`, `@enzyme_mechanism` and
`@allosteric_mechanism`, the macros themselves, the decodability check, and
`parse_mechanism(text)` (today's `_parse_mechanism`, which identify uses to read saved
mechanism text back). It builds `Node`s and `Step`s and calls the constructors, then
`Derivation.check_valid_rate_equation`. Its API: the three macros and `parse_mechanism`. The
top level re-exports the macros. The macros' hygiene (today they resolve names in
`EnzymeRates`) must keep working in a module that imports only the macro.

## 7. Performance gate

The spike (`/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/276b07d3-bafa-4735-bad9-ab925057d392/scratchpad/n2/`) put the storage-bound share of
`bench_enum` at about 24% and estimated N2 at +3–10% on enumeration. **Gate: to be set by
Denis** (the first split used ≤5% paired). It applies to `bench_enum`, `bench_guard` and
`bench_pass1` (in the scratchpad's `dm/`), each measured paired against its parent commit in
the same session (ruling R8 of the first split). It is checked first right after the storage
port, before any later step builds on it; a miss stops the work for Denis's call.
`rate_equation`'s bound is unchanged and absolute.

Memory is measured at the same checkpoint: the retained size of the census corpus (80,524
mechanisms) before and after. The spike's estimate is 1.9 kB → 3.0 kB per mechanism, about
+90 MB for the corpus. The full suite must still run in the machine's 7.7 GB.

The census gate stays: zero differences against `DM_base` in laws, parameter names and child
sets, and printed text (§Non-goals).

## 8. Developer docs

Rewrite the "Modules" section of `docs/src/developer.md` for a newcomer:
- **A diagram** of the six modules and the top level, arrows labelled with the main calls
  (for example `Enumeration → Derivation: check_valid_rate_equation, parameters`).
- **One short paragraph per module**: what it is for and what it hides.
- **One table per module** listing every API name, a one-line description and its users.
- **The mechanism data model**: the `Node`/`Step`/`Mechanism` field tables from §3 and the
  rule "stored data with dots, derived data with functions".
- **One request traced through the modules**: `identify_rate_equation` from the reaction to
  the selected fit.

Also update `docs/src/api.md`, CLAUDE.md's "Module boundaries" paragraph, and ROADMAP entries
this work closes.

## 9. Boundary test

- Rules 1 (API lists), 2 (layering, now six modules), 4 (no underscore names across modules)
  and 5 (no `export` in a submodule) stay.
- Rule 3 (the field ban) and `_testhelper_FIELD_ALLOWED` are deleted: every field a sibling
  can reach belongs to an API type whose fields are documented.

## Left to the builder

- Task order, with the storage port measured before the steps that build on it.
- The internal names of helpers.
- Whether `Reactant.atoms` is a `Dict` or a vector of pairs.
- How tests of deleted accessors are rewritten (to fields), within the no-coverage-loss rule.

## Calls made by Claude, for Denis to check

- The `Data` module and its `RateData` type (the names).
- `test_rate_equation` and the three parameter-kind types are exported.
- The diagram is a text diagram in a code block, which renders on GitHub and in Documenter and
  diffs cleanly.
- `Mechanism`'s constructor takes steps with integer node indices into a node list that may
  hold duplicates, so there is one `Step` type for reading and writing.
- `parameter_count` stays only if Enumeration needs its unchecked count.
- `Regulator` keeps today's activator/inhibitor declaration as a third field, `effect`.
- The form-name rendering stays in `Mechanisms` behind one new function, `step_sides`,
  rather than being copied into `Derivation` or moving the constructor's distinct-name check.
- Printed text stays byte-identical (§Non-goals).
