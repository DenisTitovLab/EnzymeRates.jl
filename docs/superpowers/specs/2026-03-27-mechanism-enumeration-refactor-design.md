# Mechanism Enumeration Refactor Design

## Goal

Replace the current beam enumeration pipeline (`beam_enumeration.jl`) and old staged pipeline (`mechanism_enumeration.jl`) with a single, simpler mechanism enumeration module. The new design grows mechanisms incrementally by parameter count using composable building blocks rather than a monolithic pipeline.

## File Changes

### Renames (preserve old code)
- `src/mechanism_enumeration.jl` → `src/old_mechanism_enumeration.jl`
- `src/beam_enumeration.jl` → `src/old_beam_enumeration.jl`
- `test/test_beam_enumeration.jl` → `test/old_test_beam_enumeration.jl`

### New files
- `src/mechanism_enumeration.jl` — self-contained: types, init, expansion moves, dedup, type constructors
- `test/test_mechanism_enumeration.jl` — unit tests per move + integration tests per reaction

The new `mechanism_enumeration.jl` must be fully self-contained, copying any needed helpers from old files rather than depending on them.

After implementation, update CLAUDE.md architecture notes (mechanism enumeration pipeline section, `MechanismSpec` field descriptions, etc.) to reflect the new design.

## Type Changes

### `EnzymeReaction`: add `oligomeric_state` type parameter

`EnzymeReaction` is a singleton type (`struct EnzymeReaction{S,P,R} end`) with all data in type parameters. Add `oligomeric_state` as a 4th type parameter: `EnzymeReaction{S,P,R,N}` where `N::Int` (default 1 for monomer). This preserves the singleton property and compatibility with `@generated` functions.

This value is used as the subunit count for allosteric mechanisms — both catalytic subunits and regulator site multiplicities. It is user-provided, not enumerated. The `@enzyme_reaction` DSL should accept an optional `oligomeric_state:` label.

### Existing types carried forward

- **`AbstractMechanismSpec`** — abstract supertype for `MechanismSpec` and `AllostericMechanismSpec`
- **`StepSpec`** — `(reactants, products, is_equilibrium)`, unchanged
- **`ParamConstraint`** — `Tuple{Symbol, Int, Vector{Tuple{Symbol, Int}}}`, unchanged
- **`MechanismSpec`** — 4 fields: `reaction`, `steps`, `param_constraints`, `max_estimated_param_count`
- **`AllostericMechanismSpec`** — same fields as current. `catalytic_n` comes from `reaction.oligomeric_state`. All regulator site multiplicities equal `catalytic_n`.

### `compile_mechanism` removal

`compile_mechanism` is currently exported but should not have been. Remove the export and replace all call sites with the new type constructors `EnzymeMechanism(spec)` and `AllostericEnzymeMechanism(spec)`. No backward compatibility alias.

## Architecture: No Monolithic Pipeline

There is no `enumerate_mechanisms` function. The module provides three building blocks:

1. **`init_mechanisms(reaction)`** — produces the initial pool
2. **`expand_mechanisms(specs, reaction)`** — produces +1 and +2 expansions
3. **`dedup!(cache)`** — removes structural duplicates

The caller owns the loop and the cache. This allows `identify_rate_equation` (future) to control early stopping, fitting budget, and logging.

### Caller loop pattern

```julia
cache = Dict{Int, Vector{AbstractMechanismSpec}}()

init_specs = init_mechanisms(reaction)
min_pc = first(init_specs).max_estimated_param_count
cache[min_pc] = init_specs
dedup!(cache)

for pc in min_pc:max_params
    level = pop!(cache, pc, AbstractMechanismSpec[])
    isempty(level) && (isempty(cache) ? break : continue)

    # Use this level (compile, fit, evaluate, collect, etc.)
    # ...

    new_specs = expand_mechanisms(level, reaction)
    for (target_pc, specs) in new_specs
        target_pc > max_params && continue
        append!(get!(cache, target_pc, AbstractMechanismSpec[]), specs)
    end
    dedup!(cache)
end
```

**Termination**: stops when cache is empty (natural exhaustion) or `pc > max_params` (safety cap). Both conditions apply.

## `init_mechanisms(reaction) → Vector{MechanismSpec}`

Produces all mechanisms at the minimum parameter count.

### Algorithm

1. Enumerate all catalytic topologies via backtracking (same algorithm as current `_catalytic_topologies`)
2. For each topology:
   - Assign exactly 1 SS step (first isomerization), rest RE
   - Add all equivalence constraints (all K's for same metabolite constrained equal)
   - Find all substrate/product dead-end opportunities (forms where a substrate or product can bind off-cycle, subject to mixed-binding constraint)
   - Enumerate all 2^n subsets of dead-end opportunities, each producing a distinct mechanism
   - For each subset, add mirror steps where both endpoints of a catalytic step have dead-end extensions with the same metabolite (mirror inherits RE/SS from catalytic step)
   - All dead-end K's constrained to catalytic counterpart

### Param count invariant

Every mechanism from `init_mechanisms` has:

```
max_estimated_param_count = n_substrates + n_products + 3
```

Where +3 = kcat + Keq + E_total. Topology complexity (extra forms, steps, thermo constraints) cancels out. Dead-end extensions with constrained K add zero. Mirror steps create cycles in the enzyme form graph whose Wegscheider conditions are automatically satisfied by the equivalence constraints (e.g., if K_R is constrained equal on forms E and ES, the cycle E→ES→ESR→ER→E has a trivially redundant Wegscheider condition), contributing zero additional parameters.

### Form naming

- Catalytic forms: `_form_name(bound_subs, bound_prods, has_residual)` — sorts metabolites, produces e.g. `E_A_B`, `Estar_P`
- Dead-end forms: `_dead_end_form_name(base_bound, added_met)` — sorts all bound metabolites
- Regulator forms: use `__reg` suffix to avoid collision when metabolite is both substrate/product and regulator (e.g., `X__reg1`)

## Parameter count: estimated vs. actual

`MechanismSpec.max_estimated_param_count` is an upper bound computed cheaply during enumeration. It assumes each expansion move adds exactly +1 or +2 parameters. However, removing an equivalence constraint can un-trivialize a Wegscheider thermodynamic condition, making the actual gain +0 instead of +1. For example, in a bi-bi random mechanism with K_A constrained equal on two forms, the Wegscheider condition K1*K3 = K2*K4 is redundant; removing the constraint makes it active, consuming the freed parameter.

The enumeration pipeline does not attempt to compute exact parameter counts — that complexity is deferred to the fitting pipeline. When a mechanism is compiled via `EnzymeMechanism(spec)`, `parameters()` returns the true parameter list accounting for all thermodynamic constraints. The fitting pipeline organizes mechanisms by their true `length(parameters(m))`, not by `max_estimated_param_count`.

Consequence: the enumeration cache is keyed by `max_estimated_param_count`, so mechanisms are processed in approximately-correct order (simpler before more complex, since the estimate is an upper bound). A few mechanisms may appear at a slightly higher level than their true param count warrants — this is acceptable.

## `expand_mechanisms(specs, reaction) → Dict{Int, Vector{AbstractMechanismSpec}}`

Applies all expansion moves to input specs. Returns results grouped by target `max_estimated_param_count`.

### Move 1: RE→SS conversion (+1)

Convert one RE step to SS. Each eligible RE step produces one new mechanism.

- Skip RE steps whose K is involved in an equivalence constraint (breaking constraints is a separate move)
- Dead-end steps that mirror this catalytic step inherit the new SS status
- Mechanisms with all SS steps yield no new mechanisms from this move
- After Move 2 removes a constraint at level n, the previously-constrained RE step becomes eligible for Move 1 at level n+1

### Move 2: Remove equivalence constraint (+1)

Remove one equivalence constraint, making previously-shared K (or kf/kr) parameters independent.

- Each removable constraint produces one new mechanism
- Substrate/product K constraints and regulator K constraints are independent — removing one does not affect the other (they bind to different sites)
- Mechanisms with no constraints yield nothing from this move

### Move 3: Add dead-end regulator (+1)

Add a new regulator (not yet in the mechanism) as a dead-end inhibitor.

- Applicable to regulators with `:unknown` or `:dead_end` role
- Find all eligible catalytic forms where the regulator can bind
- Enumerate all non-empty subsets of eligible forms (2^n - 1 variants), each at +1
- Each variant adds k binding steps (one per form in the subset) but only 1 free K parameter because all K_R are constrained equal, hence +1 regardless of subset size
- Add mirror steps where regulator binds both endpoints of a catalytic step (mirror inherits RE/SS)
- When metabolite is both substrate/product and regulator: regulator CAN bind to forms where the same metabolite is already bound as substrate/product (different binding site). Use `__reg` suffix in form naming.

### Move 4: Add allosteric regulator (+1)

Add a new regulator to an already-allosteric mechanism's T/R state binding.

- Only applicable to `AllostericMechanismSpec`
- Applicable to regulators with `:unknown` or `:allosteric` role
- Three flavors, each +1:
  - `r_only`: regulator binds only R-state (K_R only)
  - `t_only`: regulator binds only T-state (K_T only)
  - `tr_equiv`: regulator binds both states equally (K_T = K_R, one K parameter)
- Site assignment: same site as an existing regulator (shared denominator factor) or new site (separate factor)
- When metabolite is both substrate/product and allosteric regulator: substrate/product TR equivalences and regulator TR equivalences are independent

### Move 5: Remove TR equivalence (+1)

Remove one T/R equivalence from an `AllostericMechanismSpec`, making T-state and R-state parameters independent.

- Only applicable to `AllostericMechanismSpec`
- For metabolite K's: TR-equivalent K_T = K_R becomes independent K_T and K_R (+1)
- For catalytic rate constants (small k, SS steps): TR-equivalent kf_T = kf_R becomes independent (+1)
- When same metabolite is substrate/product and regulator: removing TR-equiv for substrate/product K is separate from removing TR-equiv for regulator K — each produces a distinct mechanism
- Ping-pong mechanisms with multiple SS steps: TR equivalences for different SS steps can be removed in stages (one per level)
- Specs with no TR equivalences left yield nothing from this move

### Move 6: Allosteric conversion (+2)

Convert a non-allosteric `MechanismSpec` to `AllostericMechanismSpec`.

- Only applicable to `MechanismSpec` (non-allosteric)
- +2 = L (allosteric constant) + one metabolite differentiation
- The one differentiated metabolite can be any substrate, product, or regulator with `:unknown` or `:allosteric` role. Regulators already added as dead-end inhibitors in the current mechanism are excluded. Regulators not yet in the mechanism are eligible (allosteric conversion introduces them). The differentiated metabolite is made either `r_only` or `t_only`.
- All other metabolites start TR-equivalent (K_T = K_R)
- All non-binding SS catalytic steps start TR-equivalent (kf_T = kf_R)
- `oligomeric_state` read from `reaction`
- Already-allosteric specs yield nothing from this move
- Rationale for minimum +2: all-TR-equiv is biochemically pointless (L cancels out, allostery has no effect). At least one metabolite must differentiate T from R for the allosteric model to be meaningful.

## Deduplication

### Canonical form

Canonicalize `MechanismSpec` / `AllostericMechanismSpec` for structural equality comparison:

- Steps sorted by canonical key: `(sorted(reactants), sorted(products), is_equilibrium)`
- Param constraints sorted by target symbol
- Form names are already canonical (metabolites sorted during construction)

### `dedup!(cache)`

- Canonicalize all specs in each max_estimated_param_count bucket
- Remove duplicates by structural equality within each bucket
- No semantic/fingerprint-based dedup (may result in some kinetically-equivalent mechanisms with different structures surviving — acceptable, address later if needed)
- Called once per loop iteration after merging expansion results into cache

## Type Constructors

Replace `compile_mechanism` with Julian type constructors:

- **`EnzymeMechanism(spec::MechanismSpec)`** → `EnzymeMechanism`
- **`AllostericEnzymeMechanism(spec::AllostericMechanismSpec)`** → `AllostericEnzymeMechanism`

These convert specs to the compiled type-parameter-encoded mechanism types used by `rate_equation` and other `@generated` functions.

## Test Plan

All functions are internal (not public API). Tests use `EnzymeRates.` prefix to access internals.

### Test helpers

- **`mechanism_spec_from_mechanism(m, rxn)`** — copied from old tests. Converts compiled `EnzymeMechanism` back to `MechanismSpec`. Its correctness is verified by the round-trip test: define mechanism with `@enzyme_mechanism`, convert to `MechanismSpec`, compile back, verify `=== original`.
- **`enumerate_all(reaction; max_params)`** — simple loop (as shown in caller loop pattern above) that collects all levels into a `Dict{Int, Vector}`.

### Unit tests: `init_mechanisms`

**Catalytic topology counts & round-trip** (adapted from old Stage 1 tests):
- Define each expected topology with `@enzyme_mechanism`
- Uni-uni: 1 topology, verify `EnzymeMechanism(topo) === m_uu`
- Uni-bi: 3 topologies, verify each matches hand-defined mechanism
- Bi-bi: 9 topologies
- Bi-bi ping-pong: 10 topologies
- All have exactly 1 SS step

**Param count invariant**:
- For each reaction type, every output has `max_estimated_param_count == n_substrates + n_products + 3`

**Dead-end saturation** (adapted from old Stage 3a tests):
- Uni-uni: no dead-end forms possible (passthrough)
- Bi-bi random topology: 4 unique dead-end forms → 16 variants per topology
- Uni-bi ordered: no mixed forms possible (passthrough)
- Bi-bi ping-pong: 3 dead-end forms → 8 variants

**Round-trip**:
- Define mechanism with `@enzyme_mechanism`, convert to `MechanismSpec` via `mechanism_spec_from_mechanism`, compile back with `EnzymeMechanism(spec)`, verify `=== original`

### Unit tests: expansion moves

Each move tested with:
- Mechanisms that have the intended behavior (positive cases)
- Mechanisms that should yield nothing from this move (negative cases)
- Edge cases

**Move 1: RE→SS**:
- Mechanism with multiple RE steps → one new spec per eligible RE step, each +1
- Mechanism with constrained RE steps → those skipped
- Mechanism with ALL SS steps → yields nothing
- Mechanism with only 1 RE step → converts correctly

**Move 2: Remove constraint**:
- Mechanism with K equivalence constraints → removing one yields +1
- Mechanism with no constraints → yields nothing
- Edge case: same metabolite as substrate/product AND regulator → constraints are independent, removing one doesn't affect the other, no error

**Move 3: Add dead-end regulator**:
- Uni-uni + new regulator: all eligible form subsets enumerated, each +1
- Mechanism with regulator binding both endpoints → mirror step with inherited RE/SS
- Mechanism with no available regulators → yields nothing
- Edge case: metabolite is both substrate and regulator → regulator CAN bind to forms where same metabolite is bound as substrate (different site), produces correct dead-end forms with `__reg` suffix

**Move 4: Add allosteric regulator**:
- Already-allosteric spec + second regulator: r_only, t_only, tr_equiv × site options, each +1
- Non-allosteric spec → yields nothing

**Move 5: Remove TR equivalence**:
- Allosteric spec with TR-equivalent metabolite K's → removing one yields +1
- Allosteric spec with TR-equivalent catalytic k's (SS rate constants) → removing one yields +1
- Ping-pong mechanism: multiple SS steps with TR-equivalent k's → removable in stages across levels
- Edge case: same metabolite as substrate and regulator → removing substrate TR-equiv is separate from removing regulator TR-equiv, both appear as distinct results
- Spec with no TR equivalences left → yields nothing

**Move 6: Allosteric conversion**:
- Non-allosteric spec → variants for each substrate/product/eligible-regulator × {r_only, t_only}, all +2
- Already-allosteric spec → yields nothing
- Verify `oligomeric_state` comes from reaction

### Unit tests: dedup

- Same mechanism with steps in different order → canonicalized, one removed
- Two genuinely different mechanisms at same max_estimated_param_count → both preserved
- Idempotent: dedup twice = dedup once

### Integration tests

For each test reaction (uni-uni, uni-bi, bi-bi, bi-bi with regulators, bi-bi ping-pong, reactions with allosteric regulators):
- Run the full `enumerate_all` loop
- Verify monotonically increasing parameter counts across levels
- Verify every mechanism compiles and `length(parameters(m)) <= max_estimated_param_count`
- Verify no duplicate fingerprints within a level (regression check)
- Verify expected count at each level (golden values, established during implementation)
