# Plan: Refactor mechanism_enumeration.jl to ~200 lines

## Context

`src/mechanism_enumeration.jl` is 1272 lines with over-engineered algorithms: 150-line edge classifier, 300-line DFS/BFS cycle enumerator, 200+ lines of helper functions. The core algorithms (catalytic cycles, activator/dead-end/RE-SS enumeration) can be expressed far more compactly using a precomputed adjacency matrix and structured walks, eliminating the need for `edge_class`, `_find_form`, permutation DFS, and ping-pong DFS entirely.

## Key Design Decisions (from interview)

1. **EnzymeFormSpec** → singleton type `EnzymeFormSpec{Name, Sites} <: AbstractEnzymeForm`
2. **SiteState** → removed; site data encoded in `Sites` type param, accessed via `getproperty` returning NamedTuples
3. **Both `enumerate_enzyme_forms` and `adjacency_matrix`** → `@generated` functions
4. **MechanismSpec** → stores `edges::Vector{Tuple{Int,Int}}` instead of `reactions`; `reactions()` accessor computes reaction tuples on demand
5. **`edge_class`** → removed entirely; adjacency matrix replaces all edge classification
6. **Catalytic cycles** → structured walk on adjacency matrix (not generic cycle enumeration)
7. **Tests** → adapted to new API; all enumeration counts must match exactly

## What Gets Eliminated (~900 lines)

| Current code | Lines | Replacement |
|---|---|---|
| `edge_class` + `_released_metabolite` + `_residual_metabolite` + `_is_valid_isomerization` + `_core_atoms` | ~150 | `_can_react()` (~15 lines) in `_compute_adjacency` |
| `_find_form` + `_build_standard_form_set` | ~70 | Structured walk on adjacency matrix (~20 lines) |
| `_enumerate_pingpong_form_sets!` + `_pingpong_dfs!` + `_release_prods_dfs!` | ~150 | **Eliminated** — ping-pong produces same cycle form sets as standard walk |
| `_permutations` | 12 | Inline 4-line version |
| `_derive_edges` + `_edges_to_reactions` + `_spec_to_edges` + `_used_forms` + `_used_form_count` | ~90 | ~10 lines total (edges stored directly, trivial accessors) |
| `_find_shadow_form` + `_find_dead_end_form` + `_reg_site_positions` | ~50 | Simplified site comparison using extracted form data (~15 lines) |
| `SiteState` struct | 10 | NamedTuples via `getproperty` |

## Implementation: File Structure (~200 lines)

### Section 1: Types (~25 lines)

```
abstract type AbstractEnzymeForm end
struct EnzymeFormSpec{Name, Sites} <: AbstractEnzymeForm end
Base.getproperty(::EnzymeFormSpec{N,S}, field) — returns Name or NamedTuple sites

const ParamConstraint = Tuple{Symbol, Int, Vector{Tuple{Symbol, Int}}}

struct MechanismSpec — reaction, forms, form_atoms, edges, equilibrium_steps, param_constraints
struct MechanismIterator — specs::Vector{MechanismSpec}, iterate, length
```

**`getproperty` implementation**: `f.sites` returns `map(Sites) do s ... end` producing a Tuple of NamedTuples with fields `(metabolite, index, atoms, role, full_atoms)`. Atoms are converted from type-param tuples `((:C,1),)` to runtime `[:C => 1]` vectors. This preserves test compatibility: `f.sites[k].metabolite == :A` works.

### Section 2: Compile-time helpers (~45 lines)

**`_compute_forms(S, P, R)`** (~30 lines): Called at compile time inside `@generated`. Logic identical to current `enumerate_enzyme_forms` but operates on type params and returns `Vector{Tuple{Symbol, Tuple}}` of `(name, sites_data)` pairs.

- Computes ping-pong residuals via bitmask over products (~8 lines)
- Builds per-site options with (metabolite, index, full_atoms, role) (~8 lines)
- Cartesian product with exclusion filter (~10 lines)
- Builds name symbol and sites tuple (~4 lines)

**`_compute_adjacency(form_data)`** (~15 lines): Builds `Matrix{Bool}` from form data. Two forms are adjacent if:
- They differ at exactly 1 site: one has atoms=nothing, other has atoms≠nothing (binding/release)
- They differ at exactly 1 site: both have atoms≠nothing but different content, and atoms exist in product list (ping-pong release/binding)
- They differ at ≥2 core sites and represent valid isomerization (all-subs-bound ↔ all-prods-bound, atom conservation)

This replaces the entire 150-line `edge_class` infrastructure with a ~15-line `_can_react(site_a, site_b)` check.

### Section 3: @generated functions (~15 lines)

```julia
@generated enumerate_enzyme_forms(::EnzymeReaction{S,P,R})
    → builds tuple of EnzymeFormSpec{Name, Sites}() singletons

@generated adjacency_matrix(::EnzymeReaction{S,P,R})
    → returns precomputed Matrix{Bool}

n_sites(reaction) — inline sum over substrates/products/regulators
```

### Section 4: Catalytic cycle enumeration (~30 lines)

**`_enumerate_only_catalytic_mechanisms(forms, reaction; max_forms)`** — entry point.

**`_enumerate_catalytic_form_sets(adj, forms, sub_names, prod_names)`** (~20 lines):
1. Build occupancy lookup: `Dict{Tuple{UInt,UInt}, Int}` mapping `(sub_bitmask, prod_bitmask)` → form index. Only core forms (no reg sites occupied, no extra sites, no residuals).
2. For each `sub_perm × prod_perm`: structured walk on adjacency matrix. From free enzyme, bind subs one at a time (lookup intermediate forms), find isomerization neighbor, release prods one at a time. Return `Set{Int}` or `nothing`.
3. `_combine_form_sets(cycles)`: BFS union of cycle pairs (same 8-line algorithm as current).

**Key insight**: Ping-pong DFS eliminated. Research confirmed ping-pong reactions produce identical catalytic cycle form sets to standard reactions (residual forms only appear via dead-end binding).

**Helpers** (inlined or minimal):
- `_permutations(v)` — 4-line recursive
- `_add_unique!(sets, new)` — 3 lines
- `_build_spec(form_set, adj, forms, reaction)` — derive edges from adj matrix for given form set, build MechanismSpec

### Section 5: Activator configs (~20 lines)

**`_generate_activator_configs(spec, forms, reaction)`**:
1. For each regulator, find shadow pairs: form index pairs `(base, shadow)` where shadow differs from base only by having the regulator site occupied. Linear scan comparing site tuples.
2. Three options per regulator: absent (no edges), non-essential (shadow edges + base kept), essential (shadow edges, shadowed base edges removed).
3. Cartesian product across regulators. Build MechanismSpec per combo.

Shadow detection replaces `_find_shadow_form` (~20 lines) with ~5-line site comparison using extracted form site data.

### Section 6: Dead-end configs (~20 lines)

**`_enumerate_dead_end_configs(spec, forms; max_forms)`**:
1. Per topology form: find empty regulatory site positions (site.role == :reg && site.atoms === nothing).
2. Enumerate regulator subsets via bitmask. Box closure: for each chosen subset, verify all sub-subsets have corresponding forms (bitmask iteration).
3. Cartesian product with max_forms budget. Build specs with topology edges + binding-only edges between all forms.

`_find_dead_end_form` inlined: for base form + set of positions to fill, linear scan for matching form.

### Section 7: RE/SS + constraints (~20 lines)

**`_enumerate_ress_and_constraints(spec, forms)`**: Structurally identical to current. Bitmask over edges (exclude all-RE). Find equivalent binding groups. Per valid group: constrained or unconstrained option. Build MechanismSpec.

**`_find_equivalent_groups(edges, forms)`**: ~10 lines. For each binding edge, identify which metabolite/site it binds (site comparison). Group edges binding same metabolite at same site index.

### Section 8: Pipeline + conversion (~25 lines)

**`enumerate_mechanisms(reaction; max_forms)`** — chains stages 4→5→6→7, returns `MechanismIterator`. (~10 lines)

**`_used_forms(spec)`** — `Set(spec.forms[i] for i in Set(Iterators.flatten(spec.edges)))` (1 line)

**`_used_form_count(spec)`** — `length(Set(Iterators.flatten(spec.edges)))` (1 line)

**`reactions(spec, forms)`** (~8 lines) — accessor computing reaction tuples from edges. For each `(a, b)` edge: compare sites to determine binding/release/isomerization, construct `([LHS...], [RHS...])` tuple.

**`EnzymeMechanism(spec)`** (~8 lines) — calls `enumerate_enzyme_forms(spec.reaction)` to get forms, calls `reactions(spec, forms)`, builds species/reactions/eq_steps/constraints tuples, delegates to constructor.

## Line Budget

| Section | Lines |
|---|---|
| 1. Types + accessors | 25 |
| 2. Compile-time helpers | 45 |
| 3. @generated functions | 15 |
| 4. Catalytic cycles | 30 |
| 5. Activator configs | 20 |
| 6. Dead-end configs | 20 |
| 7. RE/SS + constraints | 20 |
| 8. Pipeline + conversion | 25 |
| **Total** | **200** |

## Files to Modify

1. **`src/mechanism_enumeration.jl`** — complete rewrite (1272 → ~200 lines)
2. **`src/EnzymeRates.jl`** — update exports: remove `SiteState`, add `AbstractEnzymeForm`, `adjacency_matrix`
3. **`test/test_mechanism_enum_of_enz_reaction.jl`** — adapt tests:
   - `spec.reactions` → `spec.edges` (lines 319, 383, 399, 405, 451, 539)
   - `EnzymeRates._spec_to_edges(spec, forms)` → `spec.edges` (lines 411, 412, 419)
   - `edge_class` tests → `adjacency_matrix` tests (section 7)
   - `forms` is now a Tuple not Vector — generators/iteration still work, but type annotations may need updating
4. **`test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`** — adapt test helpers:
   - `EnzymeRates._spec_to_edges(spec, forms)` → `spec.edges` (lines 42, 111)
   - `EnzymeRates._reg_site_positions(forms[fi])` → inline or keep as internal function
   - `EnzymeRates._find_dead_end_form(fi, positions, forms)` → inline or keep
   - `EnzymeRates._find_equivalent_groups(edges, forms)` → keep (still exists)

## Implementation Order

1. Write new `mechanism_enumeration.jl` with all sections
2. Update `EnzymeRates.jl` exports
3. Adapt test files
4. Run tests via MCP REPL, fix issues iteratively
5. Run full test suite, verify all 8 enumeration spec counts match
6. Profile compilation time for 7-site reaction (must be < 1s cold, < 0.1s warm)

## Verification

```julia
# In MCP REPL:
include("test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl")
include("test/test_mechanism_enum_of_enz_reaction.jl")
```

Critical checks:
- All 8 `ENUMERATION_TEST_SPECS` pass with exact counts
- `EnzymeMechanism(spec)` construction works for all specs
- `rate_equation_string` works for converted mechanisms
- Compilation time test passes (7-site reaction < 1s cold)
- `wc -l src/mechanism_enumeration.jl` ≤ 200

## Risks

1. **Ping-pong cycle elimination**: Research confirmed ping-pong Bi-Bi produces same 9 catalytic topologies as standard Bi-Bi. If a future reaction depends on true ping-pong cycles, this would need revisiting. Mitigated by test coverage.
2. **Line target**: 200 lines is tight. If proper documentation pushes to 210-220, I'll flag it and let you decide what to cut.
3. **@generated compilation cost**: 184-form reaction creates a 184-element tuple literal. The existing test requires < 1s cold. Will profile and fall back to Vector if needed.
