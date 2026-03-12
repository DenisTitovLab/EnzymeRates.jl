# Mechanism Enumeration Refactor: Constructive Builder with Simplified Representation

## Goal

Simplify `src/mechanism_enumeration.jl` (~1253 lines) by replacing the
"enumerate all forms, build adjacency, search for cycles" paradigm with direct
constructive building of mechanisms. Target: eliminate ~400-500 lines of
infrastructure code (`enumerate_enzyme_forms`, `_classify_edge`,
`_build_adjacency`, `SiteDefinition`, `EnzymeFormSpec`, etc.) while making the
remaining code more readable.

**Out of scope:** Allosteric stages 6-8 (~100 lines) stay as-is.

**Requires adaptation:** Dedup infrastructure (`_concentration_fingerprint`,
`_spanning_arborescence_monomials`, `_compute_re_partition`, ~145 lines) keeps
the same algorithms but its interface must change to work with the new
`MechanismSpec` representation. Specifically, it currently consumes form
indices + adjacency dict; it will need to consume step-based metabolite info
and form name mappings instead. The core math (spanning arborescences,
union-find) is unchanged.

## New Representation

### MechanismSpec

```julia
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{Tuple{Vector{Symbol}, Vector{Symbol}}}
    n_catalytic_steps::Int
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{ParamConstraint}
    param_count::Int
end
```

Each step is in **canonical binding direction** — metabolite always on LHS:
- Substrate binding: `([:E, :S], [:ES])`
- Product binding (catalytically a release): `([:E, :P], [:EP])`
- Isomerization: `([:EAB], [:EPQ])`

Everything needed by later stages is derivable from steps:
- **Form names**: unique first elements from each side of all steps
- **Edge metabolite**: `length(step[1]) == 2 ? step[1][2] : nothing`
- **Binding direction**: metabolite is always on LHS (canonical form)
- **Form atoms**: trace from E, add substrate atoms on binding, subtract
  product atoms on release

### Regulatory Dummy Names

When a metabolite binds at a regulatory site, use a dummy name with a
collision-resistant separator (e.g., `:X__reg1`, `:X__reg2` — double
underscore) to distinguish it from catalytic-site binding of the same
metabolite. This means:
- Equivalence grouping just groups by metabolite name — no `site_type` field
- `compile_mechanism` strips the `__regN` suffix to recover the real name
- Same metabolite at two different regulatory sites: `:X__reg1`, `:X__reg2`
- Same metabolite as substrate and regulator: `:X` (catalytic), `:X__reg1`
  (regulatory)

### AllostericMechanismSpec

Unchanged — wraps the new `MechanismSpec`.

## Eliminated Code

The following types and functions are no longer needed:
- `SiteDefinition` — per-site binding info, replaced by step metabolites
- `EnzymeFormSpec` — per-site occupancy vectors, replaced by form names
- `enumerate_enzyme_forms` — combinatorial form enumeration
- `_classify_edge` — pairwise form comparison to find valid steps
- `_build_adjacency` — O(n²) adjacency from all form pairs
- `_is_binding_direction` — occupancy comparison for binding direction
- `_find_dead_end` — form lookup by modified occupancy
- `_dead_end_catalytic_map` — mapping dead-end edges to catalytic counterparts
- `_propagate_de_eq_steps!` — copying RE/SS from catalytic to dead-end edges

## Stage 1: Constructive Catalytic Topology Builder

Build catalytic paths by direct construction using backtracking search.

### Algorithm

From free enzyme E, recursively try all valid next steps. Track accumulated
atoms on the enzyme, bound substrates, and released products.

| Enzyme state | Allowed next steps |
|---|---|
| Free enzyme (no residual, no substrates) | Bind any substrate |
| Substrates bound (no residual) | Bind another substrate, OR ping-pong isomerize (if product atoms ⊆ accumulated), OR final isomerize (if all substrates bound) |
| Residual only (E*) | Bind any remaining substrate |
| Residual + substrate(s) bound | **Must** isomerize (→ product + new residual, or → product + free enzyme) |
| Products bound (post-final-isomerize) | Release any bound product |

### Ping-Pong Rules

- Ping-pong isomerization can happen with one OR multiple substrates bound
- Check: any unbound product's atoms ⊆ currently accumulated enzyme atoms
- Two-step process: (1) isomerization (rearranges atoms, product becomes bound,
  residual stays on enzyme), (2) release product
- Multiple ping-pong steps allowed per path (reactions beyond Bi-Bi)
- After binding a substrate to a residual intermediate, must isomerize next

### Dead-End Paths

If a recursive branch has no valid next steps and hasn't returned to E, it
produces no paths. Standard backtracking — no special abort logic needed.

### Multi-Cycle Topologies

1. Enumerate all valid catalytic paths (complete cycles E → ... → E)
2. Take all power-set combinations of paths (each subset of 1 or more paths)
3. For each combination, collect all forms and edges from the included paths
4. Each combination = one catalytic topology (e.g., random-order = union of
   ordered paths with different substrate binding orderings)
5. Deduplicate by form set (different path subsets may yield the same topology)

Non-essential activator is handled by allosteric mechanisms (stage 6), not
multi-cycle catalytic topologies. Multi-cycle unions serve random-order
binding variants (different substrate/product orderings sharing forms).

## Stage 2: Dead-End Expansion

One pipeline stage, two internal functions for independent testability:
- `_expand_substrate_product_dead_ends(specs, reaction)` — substrate/product
  dead-end complexes
- `_expand_regulator_dead_ends(specs, reaction; dead_end_regs)` — regulator
  dead-end complexes

Both called from `_expand_dead_end(specs, reaction; dead_end_regs)`.

All dead-end binding follows the same constraint: **only bind to catalytic
forms where neither all substrates nor all products are bound.**

All dead-end binding edges are **always RE**.

**Note on RE/SS:** Stage 2 creates mirror edges structurally (topology only).
All edges created here start as RE. Stage 3 later assigns RE/SS to catalytic
edges via bitmask, and mirrored edges inherit from their catalytic counterpart
at that time.

### 2a. Substrate/Product Dead-Ends (`_expand_substrate_product_dead_ends`)

- For each eligible catalytic form F, each substrate S (or product P) not
  bound to F:
  - F+S must NOT already be a catalytic form
  - Validity check: resulting form must not have (all substrates + any product)
    or (all products + any substrate)
  - Create dead-end form, add binding edge (always RE)
- Create mirror edges for catalytic edges where both endpoints can have the
  dead-end metabolite added and remain valid

### 2b. Regulator Dead-Ends (`_expand_regulator_dead_ends`)

- For each eligible catalytic form F, each dead-end regulator R → create form
  FR, add binding edge F ↔ FR (always RE)
- Use dummy metabolite name (`:R__reg1`) for the binding step
- Create mirror edges: if F₁ ↔ F₂ exists catalytically, and both F₁R and
  F₂R are valid, add F₁R ↔ F₂R

### 2c. Combinations

`_expand_dead_end` applies 2a and 2b in a single combinatorial expansion:
for each eligible catalytic form, enumerate the power set over all available
dead-end metabolites (substrate dead-ends from 2a + regulator dead-ends from
2b). Dead-end forms can have multiple dead-end metabolites simultaneously
(e.g., EAPR = substrate A + product P + regulator R on free enzyme). Each
combination must satisfy the validity constraint.

## Stage 3: RE/SS Assignment

- The RE/SS bitmask enumerates over **catalytic edges only**
- Dead-end mirrored edges **inherit** RE/SS from their catalytic counterpart
  (not in the bitmask — derived after each catalytic assignment)
- Dead-end binding edges (regulator and substrate/product) are always RE
  (not in the bitmask — fixed)
- At least 1 step must be SS (all-RE = pure equilibrium, no rate equation)
- Any step type can be SS (not restricted to isomerization)
- Filter: RE-connected groups ≤ `max_re_groups` (1 group is valid)
- Compute `param_count = n_re + 2*n_ss - n_thermo_constraints + 2`
  (the +2 accounts for kcat_forward and Keq)

## Stage 4: Equivalence Constraints

Group edges by metabolite name. Edges binding the same metabolite (same name,
after dummy-name disambiguation) with the same RE/SS status can be constrained
to share parameters.

- With dummy regulatory names, grouping by metabolite name alone correctly
  separates catalytic-site from regulatory-site binding, and separates
  different regulatory sites from each other
- For each subset of equivalence groups, generate a variant with those
  constraints active
- Each active constraint reduces `param_count` (by 1 for RE, by 2 for SS,
  per constrained pair)

## Stage 5: Deduplication

Same algorithm: concentration fingerprint + constraint descriptor as dedup
key, keep mechanism with fewest parameters. The interface adapts to the new
representation:
- Build a form-name-to-index mapping from the steps
- Build adjacency (metabolite per edge) from step metabolites
- Derive binding direction from canonical step format (metabolite on LHS)
- Feed these into `_compute_re_partition`, `_concentration_fingerprint`
  (core algorithms unchanged)

## Stages 6-8: Allosteric

Unchanged. Operate on the new `MechanismSpec` via the `base` field of
`AllostericMechanismSpec`.

## compile_mechanism

Reads directly from `MechanismSpec`:
- **Form names**: collected from steps (first element of each side)
- **Form atoms**: computed by BFS from free enzyme E through the step graph.
  Uses `spec.reaction` to determine substrate/product atom definitions.
  At each binding step, if the metabolite is a substrate, add its atoms;
  if a product, subtract its atoms. Isomerization steps (no metabolite)
  don't change atoms. For branching topologies (random-order), BFS handles
  multiple paths — each form's atoms are the same regardless of path taken
  (atom conservation guarantees this).
- **Reactions**: steps are already in the canonical binding-direction format,
  convert `Vector{Symbol}` to `Tuple{Symbol,...}`
- **Equilibrium steps and param constraints**: read directly from spec

Strips dummy regulatory suffixes (`__regN`) from metabolite names when
producing the final `EnzymeMechanism` / `AllostericEnzymeMechanism`.

## Removed Types

- `EnumerationStage` (abstract type + subtypes `Catalytic`, `WithDeadEnd`,
  `FullEnumeration`): used only for the `stage` kwarg in
  `enumerate_mechanisms` to return early. Replace with simple keyword
  argument (e.g., `stop_after::Int=8`), or remove if not needed by tests.
- `MechanismIterator`: lazy wrapper around `Iterators.flatten`. Keep or
  replace with direct return of `Vector{AbstractMechanismSpec}` —
  implementation decision, not spec-critical.

## Pipeline Order

Renumbered to reflect the new execution order (dead-end expansion now
precedes RE/SS assignment so mirrored edges exist before RE/SS is assigned):

```
Stage 1: _catalytic_topologies (constructive builder)
Stage 2: _expand_dead_end
  ├─ 2a: _expand_substrate_product_dead_ends
  └─ 2b: _expand_regulator_dead_ends
Stage 3: _expand_ress_variants (RE/SS bitmask over catalytic edges, inherit to mirrors)
Stage 4: _expand_equivalence_constraints (group by metabolite name)
Stage 5: _deduplicate (concentration fingerprint)
Stage 6: _expand_allosteric (MWC variants)
Stage 7: _expand_tr_equivalence (T/R parameter equivalence)
Stage 8: _deduplicate_allosteric (remove T↔R mirrors)
```

## TDD: Bug Fixes and Test Mapping

The test file `test/test_mechanism_enumeration.jl` documents all known bugs
via `@test_broken`. The refactor should fix these using TDD — flip
`@test_broken` to `@test` as each stage is implemented. Below maps each bug
to the implementation step that fixes it.

### Bugs fixed naturally by the refactor

These bugs exist because of design flaws in the current search-based
architecture. The constructive builder eliminates the root causes.

| Bug | Test line(s) | Root cause | Fixed by |
|-----|-------------|------------|----------|
| Ter-Ter OOMs in `_catalytic_topologies` | 268 | O(n²) form enumeration + adjacency for 3+3 reactions | Stage 1: constructive builder only builds valid paths |
| RE/SS only toggles binding edges, not isomerization | 281, 302, 316 | `_find_first_isomerization` pins isom to SS, only iterates "other" edges | Stage 3: bitmask over ALL catalytic edges, any can be SS |
| Regulators bind to fully-occupied forms | 391, 409, 427, 439 | `_expand_dead_end_inhibitors` doesn't check occupancy eligibility | Stage 2b: explicit "neither all subs nor all prods" check |
| Bi-Bi dedup undercounts (3 instead of 5) | 626 | Cascading from isom toggle bug — fewer single-SS variants exist | Cascade fix from Stage 3 isom toggle fix |

### Bugs requiring explicit fixes

These are independent of the representation change and need targeted fixes.

| Bug | Test line(s) | Root cause | Fix approach |
|-----|-------------|------------|--------------|
| T/R mirror dedup broken | 851 | `_allosteric_canonical_key` doesn't map complementary TR-equiv sets to same key | Use `min(tr_set, complement_set)` as canonical form |
| param_count off by 1 (36/126 constrained Bi-Bi) | 1157 | Formula doesn't account for equivalence constraints making Wegscheider constraints redundant | Detect when param constraint reduces an independent cycle, adjust n_thermo accordingly |

### New features with test expectations

These `@test_broken` lines document expected behavior for features that
don't exist yet. The refactor implements them.

| Feature | Test line(s) | Implemented by |
|---------|-------------|----------------|
| Substrate/product dead-end expansion | 349, 358, 367 | Stage 2a: `_expand_substrate_product_dead_ends` |
| Same metabolite as substrate and regulator | 533 | Dummy name approach (`:X__reg1`) |
| Substrate/product dead-end equivalence | 540 | Stage 4: equivalence grouping over all edges including dead-end |

### TDD implementation order

Work through stages sequentially. For each stage:
1. Flip relevant `@test_broken` to `@test` (they should fail)
2. Implement the stage
3. Run tests to confirm they pass
4. Commit

```
Step 1: Stage 1 — Constructive catalytic topology builder
  Flip: line 268 (Ter-Ter)
  Verify: existing Stage 1 tests still pass with new representation
  Verify: compile_mechanism round-trips still work

Step 2: Stage 2a — Substrate/product dead-end expansion
  Flip: lines 349, 358, 367
  New tests for _expand_substrate_product_dead_ends

Step 3: Stage 2b — Regulator dead-end expansion
  Flip: lines 391, 409, 427, 439
  Verify: existing Stage 3 tests pass with eligibility fix

Step 4: Stage 3 — RE/SS assignment
  Flip: lines 281, 302, 316
  Cascade fix: line 626 (dedup undercount)
  Verify: existing Stage 2 tests pass with all-edge bitmask

Step 5: Stage 4 — Equivalence constraints
  Flip: lines 533, 540
  Verify: existing Stage 4 tests pass

Step 6: Stage 5 — Deduplication (interface adaptation)
  Verify: existing Stage 5 tests pass with new representation

Step 7: Stage 8 — T/R mirror dedup fix
  Flip: line 851
  Fix _allosteric_canonical_key independently

Step 8: param_count formula fix
  Flip: line 1157
  Fix Wegscheider/equivalence interaction

Step 9: End-to-end verification
  Re-run all end-to-end pipeline tests
  Update expected counts if bug fixes change them
  (e.g., Uni-Uni + unknown reg count may change with eligibility fix)
```

### Note on end-to-end count changes

Several end-to-end tests (lines 1038-1088) assert specific mechanism counts
(e.g., Bi-Bi = 207, Uni-Uni + unknown reg = 21). These counts will change
because:
- Bug fixes (eligibility, isom toggle) change which mechanisms are generated
- New substrate/product dead-end expansion adds more mechanisms
- param_count formula fix changes dedup winners

These counts should be re-derived after the refactor and updated in tests.

## Expected Outcome

- Significant reduction in code size (current: ~1253 lines)
- Eliminated: `SiteDefinition`, `EnzymeFormSpec`, `enumerate_enzyme_forms`,
  `_classify_edge`, `_build_adjacency`, `_is_binding_direction`,
  `_find_dead_end`, `_dead_end_catalytic_map`, `_propagate_de_eq_steps!`
- Simplified: all stage functions (no adjacency recomputation)
- New: constructive catalytic path builder, substrate/product dead-end support
- 6 bugs fixed, 3 new features implemented
- All `@test_broken` in test_mechanism_enumeration.jl flipped to `@test`
- Conceptual simplification: no "enumerate then search" paradigm, all
  information flows forward through construction
