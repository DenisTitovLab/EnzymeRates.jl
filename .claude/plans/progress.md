# Mechanism Enumeration Refactor — Progress Summary

## Status: Phase 1 Complete (all tests passing)

All 636 tests pass including all 8 `ENUMERATION_TEST_SPECS`.

## What Was Done

### Files Changed
1. **`src/mechanism_enumeration.jl`** — rewritten (1534 → 1020 lines, 33% reduction)
2. **`src/EnzymeRates.jl`** — removed 7 mechanism enumeration exports (now internal)
3. **`test/test_mechanism_enum_of_enz_reaction.jl`** — adapted to new `stage=` API
4. **`test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`** — adapted helpers for edge-based `MechanismSpec`

### Architectural Changes
- **`SiteState` struct → NamedTuple** with 4 fields: `(metabolite, atoms, role, full_atoms)`. Dropped `index` field (multi-site support removed).
- **`MechanismSpec` now stores `edges::Vector{Tuple{Int,Int}}`** (form index pairs) instead of `reactions::Vector{Tuple{Vector{Symbol}, Vector{Symbol}}}`. This eliminated `_spec_to_edges`, `_build_spec_from_edges`, and the `shared_names`/`shared_atoms` plumbing.
- **Precomputed adjacency graph** via `_build_adjacency(forms) → Dict{Tuple{Int,Int}, EdgeInfo}`. The `EdgeInfo` NamedTuple stores `(type::Symbol, metabolite::Union{Nothing,Symbol})`. `_classify_edge` replaces old `edge_class`.
- **Stage parameter API**: `enumerate_mechanisms(rxn; stage=Catalytic())` replaces `enumerate_mechanism_stages(rxn)`. Stage types: `Catalytic`, `WithActivator`, `WithDeadEnd`, `Full`.
- **Lazy pipeline**: stages 3-4 use `Iterators.flatten(Iterators.map(...))`. `MechanismIterator` wraps the lazy iterator with O(1) `length()`.
- **`EnzymeMechanism(spec)`** now calls `enumerate_enzyme_forms` + `_build_adjacency` internally to reconstruct forms/edges for conversion.

### What Was Dropped
- `PreRessEntry` struct and custom state machine iterator
- `_spec_to_edges`, `_edges_to_reactions`, `_build_spec_from_edges` conversion functions
- `shared_names`/`shared_atoms` parameter threading
- `_used_forms`/`_used_form_count` (replaced by `_used_form_set`)
- `n_sites` function (inlined as `length(substrates) + length(products) + length(regulators)`)
- `_free_enzyme_index` (inlined)
- Multi-site loops (`for i in 2:nsites(s)`)
- `enumerate_mechanism_stages` (replaced by `stage=` parameter)

### What Was Kept (and why)
- **Permutation-based standard cycle enumeration** (`_build_standard_form_set`, `_permutations`) — the unified DFS approach from the original plan was attempted but doesn't work. A graph-walking DFS finds cycles that are topologically valid but not biochemically valid (e.g., product-only cycles with no substrate involvement, or PP cycles that don't complete full catalytic turnovers). The permutation-based approach guarantees stoichiometric completeness.
- **Constrained PP DFS** (`_pingpong_dfs!`, `_release_prods_dfs!`) — needed for correct ping-pong cycle enumeration. The `n_subs_bound > n_prods_released` guard is essential.
- **`_is_pure_topology` filter** — still needed to reject mixed standard+PP topologies.
- **`_is_valid_isomerization` and `_core_atoms`** — complex but essential for correct PP isomerization detection.

## Phase 2 Complete (all tests passing)

All 636 tests pass. File reduced from 1020 → 892 lines (13% further reduction).

### Phase 2 Changes

**Form index dictionary** (`_build_form_index`, `_lookup_form`): O(1) form lookup by site-atoms tuple, replacing O(n) linear scans. Enabled simplification of 3 lookup functions:
- `_find_form` (31 lines) → `_lookup_form` (20 lines, shared by all callers)
- `_find_shadow_form` (18 lines) → `_find_shadow` (9 lines)
- `_find_dead_end_form` (12 lines) → `_find_dead_end` (7 lines)

**`_is_valid_isomerization` simplified** (52 → 29 lines): Single-pass loop collects diffs, then standard case uses `a_sub_occ != a_prod_occ` check (safe because `enumerate_enzyme_forms` guarantees all-subs-full ↔ no-prods-occupied invariant).

**`_classify_edge` simplified** (39 → 36 lines): Merged binding/release paths for sub/prod and reg sites. Inlined `_residual_metabolite` (removed as separate function).

**`_is_pure_topology` rewritten** (35 → 14 lines): Replaced imperative boolean tracking with declarative helper lambdas (`_is_residual`, `_is_free`, `_all_subs_full`).

**Atom summing unified**: `_core_atoms` and `_form_atoms` now share `_sum_atoms` helper (21 → 11 lines).

**Other cleanups**:
- Inlined `_enumerate_pingpong_form_sets!` into `_catalytic_topologies`
- Removed unused `adj` parameter from `_expand_activators`
- Merged sub/prod/reg site-data loops in `enumerate_enzyme_forms`
- Simplified dead-end edge construction
- Simplified activator Cartesian product edge merging

### Failed approach: Cycle pre-separation
Attempted to separate std/PP cycles before `_combine_form_sets` to eliminate `_is_pure_topology`. Failed because PP cycles validly combine with each other through `_combine_form_sets` (Bi-Bi PP: 12 vs expected 10). `_is_pure_topology` is needed to reject invalid combinations.

## Phase 3 Complete (all tests passing)

All 636 tests pass. File reduced from 892 → 750 lines (16% further reduction).

### Phase 3 Changes

**Key algorithmic insight**: Verified empirically across all 8 test specs that `_pingpong_dfs!` finds ALL standard cycles (standard cycles ⊆ PP DFS cycles). This enabled removing `_permutations` (12 lines), `_build_standard_form_set` (26 lines), and `_add_unique!` (6 lines). `_catalytic_topologies` now calls only `_pingpong_dfs!`.

**Extracted `_atom_residual` helper**: Shared between `enumerate_enzyme_forms` (residual computation) and `_pingpong_dfs!` Option 3, eliminating duplicated atom subtraction logic.

**Precomputed `prod_full` dict**: Built once in `_catalytic_topologies` and passed to `_pingpong_dfs!`, replacing repeated per-form lookups.

**Inlined `_derive_edges`** into `_catalytic_topologies`: Edges are now built directly from the form set + adjacency dict.

**Inlined `_edges_to_reactions`** into `EnzymeMechanism(spec)`: Simplified `enzs_t` generator and release branch in edge classification.

**Simplified `_expand_activators`**: Dict comprehension for shadow map, `filter!` instead of manual loop, more compact mirror/mirrored construction.

**Simplified `_dead_end_configs`**: Removed `seen` dedup set and `valid` flag — enzyme forms are a Cartesian product of per-site options, so if individual reg bindings exist, all subset combinations also exist.

**Simplified `_find_equivalent_groups`**: Uses `adj` metabolite info directly instead of re-scanning form sites.

**Simplified `_build_constraints`**: Unified RE/SS cases with prefix/suffix pattern.

**Simplified `_is_valid_isomerization`**: Removed redundant `a_sub_occ != a_prod_occ` check — guaranteed by `enumerate_enzyme_forms`' exclusion filter.

**Simplified `_is_pure_topology`**: Collected residual forms once, then checked conditions.

### Remaining functions (algorithmically necessary)

1. **`_pingpong_dfs!`** (~85 lines): 3-option DFS. Inherently complex.
2. **`enumerate_enzyme_forms`** (~75 lines): Residual computation + Cartesian product.
3. **`_dead_end_configs`** (~65 lines): Two-level bitmask loops + edge construction.
4. **`_expand_activators`** (~55 lines): Shadow pair logic + Cartesian product.

## Key Insight: Why Unified DFS Doesn't Work

The original plan proposed a single DFS on the adjacency graph to find all catalytic cycles. This was attempted and found to produce incorrect results:

- For Uni-Bi (A[C2] → P1[C] + P2[C]): DFS found 9 cycles instead of expected 3
- Spurious cycles include product-only rings (E ↔ E_P1 ↔ E_P1_P2 ↔ E_P2 ↔ E) that don't consume any substrate
- The `_is_pure_topology` filter catches some but not all spurious cycles
- Root cause: form set adjacency is undirected and doesn't encode catalytic directionality
- A stoichiometry-based post-filter was considered but requires tracking cycle paths (not just form sets), adding complexity without saving lines
