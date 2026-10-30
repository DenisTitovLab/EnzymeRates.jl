# Fewer shared names: design

Branch `deep-modules-simplification` at 111661fc, after the deeper-modules work
(`2026-10-28-deeper-modules-design.md`). Same PR, version 0.10.0.

## Intent

Each module should be understandable from its own interface alone: someone working in one
module should need to know almost nothing about another. The deeper-modules work cut the
interfaces down; several boundaries still carry 6–9 functions. This round cuts every boundary
to a handful of non-exported functions, by storing facts on the mechanism (read with dots), by
moving functions to the module that uses them, and by deleting names that duplicate another.

## Counting rule

- Count the **non-exported functions** one module calls from another, per edge.
- Free: types; names `EnzymeRates` exports; names the top level only re-exports; standard
  protocol overloads (`==`, `hash`, `isless`, `show` of a whole value, iteration).
- Counted: overloads that act as a module's own API (a computed `getproperty`, or `show` used to
  produce message text).
- Target: 1–2 per edge. Accepted exceptions, by Denis: EnzymeRates → Derivation (3), and
  Derivation → Mechanisms and Enumeration → Mechanisms (see §2).

## Goals

1. **Mechanisms stores its facts.** The constructor computes them once and stores them as real
   fields, read with dots (no `getproperty` overloads):
   - `Mechanism`: `modifiers`, `groups` (the steps of each kinetic group);
   - `Step`: `kind` (`:isomerization`, `:binding`, `:fused_binding`, `:chemistry`), whether it
     runs toward products as stored, and its two side names (the names its constants are
     built from);
   - `AllostericMechanism`: the steps of each conformation (today's `conformation_steps`), and
     for each kinetic group and each site ligand the conformations it exists in (today's
     `excludes`). The tag rule itself stays private in Mechanisms.
   These replace the functions `modifiers`, `groups`, `step_kind`, `runs_toward_products`,
   `step_sides`, `conformation_steps` and `excludes`. The per-step facts are stored as
   vectors on the `Mechanism` (one entry per step), not as new `Step` fields, so the
   six-argument `Step` stays the writing form. Every construction path computes them, the
   singleton lift included. They are functions of the canonical fields, so `==` and `hash`
   leave them out and deduplication costs nothing extra. Where a call site needs the MWC tag
   rule on a bare tag (Enumeration's `_state_conformations`) and no stored field can serve it,
   `excludes` stays in the API and the report says why.
2. **Graph analyses both modules need live in Mechanisms.** `segments(m, on, re)` (today's
   `re_segments` and `segment_graph`, merged into one function whose record serves both
   callers) and `has_turnover_cycle(m, on)` move from Derivation, and so does
   `has_unbalanced_cycle`, which Derivation's admission runs on its law network and the seed
   screen on index arrays (Denis: keep it shared rather than duplicate the union-find).
3. **The name checks stay where they are.** `check_decodable_names` and
   `check_conformation_label` stay in Mechanisms' API for the DSL: running the decodability
   check in the constructor would cost every enumerated child about a microsecond, and an
   enumerated mechanism takes its names from a reaction already checked. The DSL reads
   `catalytic` facts through fields instead of calling `catalytic`.
4. **`compile_mechanism` and the singleton types stay as they are** (Denis: "let's not touch
   the singleton mechanisms for now").
5. **`scale_rate_constants` leaves Derivation's interface.** A new kind,
   `parameters(m, Derivation.RateConstants())`, returns the steady-state rate constant names
   (both conformations, the inactive `I_` names included, as `scale_rate_constants` scales
   today) for concrete and compiled mechanisms; Fitting scales them itself. `RateConstants` is
   not exported (Denis); as a type it does not count against the edge.
6. **`test_rate_equation` is no longer exported.** It stays in Fitting's `API` for identify
   and stays documented as internal but usable (`EnzymeRates.test_rate_equation`).
7. **Enumeration keeps its filters** (`_degenerate`, `_hyperbolic_catalysis`,
   `_redundant_copy_groups`, the gauge test, `_has_vmax`): they are search heuristics built on
   graph analyses, not derivation machinery. Derivation keeps "does this mechanism produce a
   valid rate law" (`check_valid_rate_equation`, `parameter_count`).
8. **CLAUDE.md gains the principle** (succinct): keep names shared across modules to a minimum,
   so each module is deep and can be worked on knowing only the other modules' interfaces.

## The interfaces after this round

| Edge | Non-exported functions |
|---|---|
| EnzymeRates → Enumeration | `base_mechanisms`, `expand_mechanisms` |
| EnzymeRates → Derivation | `check_valid_rate_equation`, `rate_equation_complexity`, `rate_equation_dedup_key` |
| EnzymeRates → Fitting | `test_rate_equation` |
| EnzymeRates → Mechanisms | `compile_mechanism` |
| EnzymeRates → Data, DSL | none |
| DSL → Derivation | `check_valid_rate_equation` |
| DSL → Mechanisms | `compile_mechanism`, `check_decodable_names`, `check_conformation_label` |
| Fitting → Derivation | none |
| Fitting → Mechanisms | `compile_mechanism` |
| Enumeration → Derivation | `check_valid_rate_equation`, `parameter_count` |
| Enumeration → Mechanisms | `catalytic`, `step_text`, `form_text`, `segments`, `has_turnover_cycle`, `has_unbalanced_cycle` |
| Derivation → Mechanisms | `catalytic`, `step_text`, `form_text`, `segments`, `has_turnover_cycle`, `has_unbalanced_cycle`, `compile_mechanism` |

Denis accepted these counts, with `compile_mechanism` and `has_unbalanced_cycle` kept.

`catalytic` stays a Mechanisms function (Denis: a plain/allosteric design issue for later).
`step_text` and `form_text` stay functions (a `Step` holds node indices, so it cannot print
itself; a `show(::Node)` would count the same as `form_text`). `metabolites` stays the exported
function. Mechanisms' functions: `catalytic`, `step_text`, `form_text`, `segments`,
`has_turnover_cycle`, `has_unbalanced_cycle`, `compile_mechanism`, `check_decodable_names`,
`check_conformation_label`, plus the exported `metabolites`.

## Non-goals and anti-goals

- No change to derived laws, parameter names, enumeration child sets or printed text: a census
  with zero differences proves it.
- No change to `rate_equation` or the singleton types.
- The constructor must not get meaningfully slower or allocate meaningfully more (§Performance).

## Performance

- Today (111661fc, Julia 1.13.3): `Mechanism(substrates, products, nodes, steps)` takes 1.3 µs
  and 7.8 kB for Michaelis–Menten, 4.3 µs and 22 kB for a 9-node, 13-step bi-bi seed. It cannot
  be allocation-free (a mechanism is vectors; copying them alone allocates 1.2 kB).
- **A deterministic allocation test** pins the bytes allocated per construction for a set of
  fixtures (Michaelis–Menten, a bi-bi seed, an allosteric mechanism) to at most today's value
  plus a stated margin (the plan sets it from a measurement of the stored fields' sizes).
- **Benchmarks** (enumeration, the guard, identify's pass 1) stay within 10% of 111661fc,
  paired, on Julia 1.13 (the allocation pin is measured there too). A miss stops the work for
  Denis's call.

## Verification tooling

The census and benchmark scripts were lost with `/tmp`. The plan's first task rebuilds a
compact census in `~/claude-scratchpad`: for the six census reactions, the init mechanisms and
the children of a fixed sample of parents, recording each mechanism's printed text, fitted
parameter names and rate-law text; a compare script reports any difference. The baseline is
built at 111661fc (whose behaviour equals the earlier census baseline). The three benchmarks are
rebuilt the same way.

## Docs

The generated module map updates itself (`julia test/update_module_map.jl`); the developer
page's interface tables, the API page and CLAUDE.md's "API Design" and "Module boundaries"
paragraphs follow the code.

## Left to the builder

Field names for the stored facts; how the DSL re-prefixes constructor errors; task order (the
census tooling first, the stored fields before the moves that use them).
