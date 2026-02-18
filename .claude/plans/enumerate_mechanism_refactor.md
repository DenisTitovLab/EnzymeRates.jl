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

## Remaining Work (Phase 2 — further line reduction)

The file is 1020 lines vs the 200-line target. The gap is primarily in:

1. **Cycle enumeration** (~250 lines): `_find_form`, `_permutations`, `_build_standard_form_set`, `_enumerate_pingpong_form_sets!`, `_pingpong_dfs!`, `_release_prods_dfs!`. These are algorithmically necessary. Minor cleanup possible but not dramatic.

2. **Edge classification** (~130 lines): `_is_valid_isomerization` alone is 50 lines due to the dual standard/PP logic. `_classify_edge` is 35 lines. Hard to reduce without losing correctness.

3. **Dead-end enumeration** (~80 lines): Already fairly compact. The Cartesian product with budget is inherently complex.

4. **Activator configs** (~75 lines): Shadow pair logic is compact but has many moving parts.

Potential further simplifications:
- Inline `_permutations` using `Combinatorics.permutations` (but adds dependency)
- Merge `_pingpong_dfs!` and `_release_prods_dfs!` (they share some structure)
- Simplify `_is_valid_isomerization` if PP case handling can be unified
- Review whether `_find_form` can be replaced by direct adjacency lookups

## Key Insight: Why Unified DFS Doesn't Work

The original plan proposed a single DFS on the adjacency graph to find all catalytic cycles. This was attempted and found to produce incorrect results:

- For Uni-Bi (A[C2] → P1[C] + P2[C]): DFS found 9 cycles instead of expected 3
- Spurious cycles include product-only rings (E ↔ E_P1 ↔ E_P1_P2 ↔ E_P2 ↔ E) that don't consume any substrate
- The `_is_pure_topology` filter catches some but not all spurious cycles
- Root cause: form set adjacency is undirected and doesn't encode catalytic directionality
- A stoichiometry-based post-filter was considered but requires tracking cycle paths (not just form sets), adding complexity without saving lines
