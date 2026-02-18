# Plan: Refactor mechanism_enumeration.jl to ~200 lines

## Context

`src/mechanism_enumeration.jl` is 1534 lines. The file is over-engineered:
- 139-line edge classifier split across 4 functions (`edge_class`,
  `_is_valid_isomerization`, `_residual_metabolite`, `_core_atoms`)
- 200+ lines for two separate cycle enumerators (standard + ping-pong DFS)
  that do fundamentally the same thing (walk the form graph)
- 90+ lines for a custom lazy iterator (`MechanismIterator`, `PreRessEntry`,
  state machine) when `Iterators.flatten(map(...))` suffices
- Redundant data round-trips: `_spec_to_edges` converts reaction tuples back
  to edge indices, `_edges_to_reactions` does the reverse
- `shared_names`/`shared_atoms` plumbing threaded through every function
- 42-line `_is_pure_topology` post-filter that can be eliminated by not
  mixing cycle types in the first place

## Approach

Replace the current architecture with:
1. **Precomputed adjacency graph** with full edge metadata — three clearly
   named rule functions (`_is_binding`, `_is_release`, `_is_isomerization`)
2. **Single unified DFS** for all catalytic cycles (standard AND ping-pong)
   — the graph structure naturally constrains the DFS to valid cycles only
3. **Edge-based MechanismSpec** (4 fields) — eliminates all conversion functions
4. **Lazy `Iterators.flatten(map(...))` pipeline** — eliminates custom iterator
5. **Drop multi-site metabolite support** (untested, never used)
6. **SiteState struct → NamedTuple** (drop `index` field entirely)

## Key Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Form types | Runtime structs | No compilation cost; singleton types add complexity for no benefit |
| SiteState | NamedTuple `(metabolite, atoms, role, full_atoms)` | 4 fields, no `index` (multi-site dropped) |
| MechanismSpec fields | `reaction, edges, eq_steps, constraints` | No `form_names`/`form_atoms` — reconstructed from `spec.reaction` when needed |
| Edge classification | Three named functions + `_build_adjacency` | User wants clear, readable rules in one place |
| Cycle enumeration | Single unified DFS on adjacency graph | Standard cycles = degenerate case where no residuals appear; PP cycles found naturally by same DFS |
| Cycle combination | Separation rule: standard combine (BFS), PP standalone | Investigation confirmed no valid mixed topologies exist; eliminates `_is_pure_topology` |
| Stage API | Type hierarchy: `Catalytic`, `WithActivator`, `WithDeadEnd`, `Full` | `enumerate_mechanisms(rxn; stage=Full())` replaces `enumerate_mechanism_stages` |
| Laziness | Materialize stages 1-2, lazy stages 3-4 via flatmap | Dead-end × RE/SS can produce billions of mechanisms |
| Iterator wrapper | None — return plain Julia iterators | `IdentifyRateEquationProblem` can wrap later |
| Exports | Nothing — all internal per SPEC.md | |
| Permutations | Not needed — unified DFS finds all cycles by graph walk | Eliminates Combinatorics.jl dependency |
| Multi-site (`S[C, 2]`) | Dropped | Untested, never used; users list metabolite multiple times instead |

## What Gets Eliminated (with line counts)

| Current code | Lines | Replacement |
|---|---|---|
| `SiteState` struct + docstring | 21 | NamedTuple (0 lines) |
| `edge_class` + `_residual_metabolite` + `_is_valid_isomerization` + `_core_atoms` | 139 | `_is_binding` + `_is_release` + `_is_isomerization` + `_build_adjacency` (~32 lines) |
| `_find_form` + `_build_standard_form_set` | 68 | Unified DFS walks adjacency directly |
| `_enumerate_pingpong_form_sets!` + `_pingpong_dfs!` + `_release_prods_dfs!` | 160 | Unified DFS (~20 lines total with standard cycles) |
| `_is_pure_topology` | 42 | Eliminated — separation rule |
| `_combine_form_sets` (general BFS) | 20 | Simplified to standard-only combination (~8 lines) |
| `_permutations` | 12 | Not needed — DFS finds all cycles |
| `_derive_edges` + `_edges_to_reactions` + `_spec_to_edges` | 47 | Eliminated — specs store edges directly |
| `_used_forms` + `_used_form_count` | 11 | Inlined: `Set(Iterators.flatten(spec.edges))` |
| `_build_spec_from_edges` + `shared_names`/`shared_atoms` plumbing | 16+ | Eliminated |
| `MechanismIterator` + `PreRessEntry` + state machine | 90 | `Iterators.flatten(map(...))` (~15 lines) |
| `enumerate_mechanism_stages` | 35 | `enumerate_mechanisms(; stage=...)` |
| `n_sites` | 12 | Inlined as metabolite count |
| Multi-site loops (`for i in 2:nsites(s)`) | 6 | Dropped |
| `_free_enzyme_index` | 3 | Inlined |
| `_form_atoms` helper | 10 | Kept (needed for conversion) |

## File Structure: mechanism_enumeration.jl (~200 lines)

### Section 1: Types + Constants (~15 lines)

```julia
const EdgeInfo = @NamedTuple{type::Symbol, metabolite::Union{Nothing,Symbol}}
const ParamConstraint = Tuple{Symbol, Int, Vector{Tuple{Symbol, Int}}}

struct EnzymeFormSpec
    name::Symbol
    sites::Vector{@NamedTuple{metabolite::Symbol,
        atoms::Union{Nothing,Vector{Pair{Symbol,Int}}},
        role::Symbol, full_atoms::Vector{Pair{Symbol,Int}}}}
end

struct MechanismSpec
    reaction::Any
    edges::Vector{Tuple{Int,Int}}
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{ParamConstraint}
end

abstract type EnumerationStage end
struct Catalytic    <: EnumerationStage end
struct WithActivator <: EnumerationStage end
struct WithDeadEnd   <: EnumerationStage end
struct Full          <: EnumerationStage end
```

### Section 2: Edge Classification + Adjacency (~32 lines)

Three short named rule functions, all in one place, clearly readable:

**`_is_binding(sa, sb)`** (~3 lines): `sa.atoms === nothing && sb.atoms !== nothing`
→ returns `sa.metabolite` or `nothing`

**`_is_release(sa, sb, form)`** (~7 lines): `sa.atoms !== nothing && sb.atoms === nothing`
Handles standard release (full_atoms → empty) and residual release (residual
atoms → empty: find which product's `full_atoms` matches the residual content).
→ returns metabolite or `nothing`

**`_is_isomerization(fa, fb, diffs)`** (~10 lines): All diffs at `:sub`/`:prod`
sites. Standard case: all core sites differ, all-subs ↔ all-prods pattern.
PP case: has residual content, atom balance conserved via `_core_atoms`.
→ returns `Bool`

**`_core_atoms(form)`** (~4 lines): Dict summing atoms across `:sub`/`:prod` sites.

**`_build_adjacency(forms)`** (~8 lines): For each pair, compute diffs. Single
diff → try `_is_binding`/`_is_release`. Multiple diffs → try `_is_isomerization`.
Store `EdgeInfo` in `Dict{Tuple{Int,Int}, EdgeInfo}`.

### Section 3: Form Enumeration (~35 lines)

**`enumerate_enzyme_forms(reaction)`**: No multi-site. Each metabolite gets
exactly one site. Steps:
1. Compute ping-pong residuals: for each substrate, enumerate which product
   atom subsets can be subtracted leaving a non-empty, non-full residual (~12 lines)
2. Build per-site options: `[(nothing, "0"), (full_atoms, name), residuals...]` (~6 lines)
3. Cartesian product with exclusion filter: reject forms where all subs full
   + any prod occupied, or vice versa (~12 lines)
4. Build `EnzymeFormSpec` per valid combo (~5 lines)

### Section 4: Catalytic Cycle Enumeration (~25 lines)

**`_find_cycles(adj, forms)`** (~5 lines): Identify core forms (reg sites all
empty), find free enzyme index, call DFS, split results into standard vs PP.

**`_cycle_dfs!(cycles, adj, core, free_idx, path, visited)`** (~12 lines):
Start at free enzyme. At each form, try all neighbors in adjacency:
- Skip non-core forms
- If neighbor is free enzyme and path ≥ 3: record cycle form set
- If neighbor not yet visited: recurse

This naturally finds BOTH standard and PP cycles — the adjacency graph's
isomerization rules (all core sites differ + atom balance) prevent invalid
partial cycles. Standard cycles are the special case where no residual forms
are visited.

**`_combine_standard_cycles(std_cycles)`** (~8 lines): BFS union of pairs.
PP cycles are appended as individuals (no combination across types).

### Section 5: Activator Configs (~22 lines)

**`_activator_configs(spec, adj, forms, reaction)`**:
For each regulator, scan adjacency for shadow pairs (topo form → form with
reg also bound). Three options per reg: absent, non-essential (shadow cycle +
base kept), essential (entry binding + shadow cycle, base edges removed).
Cartesian product across regulators. Build `MechanismSpec` per combo.

### Section 6: Dead-End Configs (~22 lines)

**`_dead_end_configs(spec, adj, forms; max_forms)`**:
Identify inhibitor positions (reg sites never occupied in topo). Per topo form:
enumerate regulator subsets via bitmask, verify sub-subsets have forms in
adjacency. Cartesian product with `max_forms` budget. Build `MechanismSpec`
with topo edges + binding edges.

### Section 7: RE/SS Lazy Generator (~22 lines)

**`_ress_variants(spec, adj, forms)`**: Returns a lazy generator.
Compute equivalent binding groups from adjacency metadata. For each valid
RE/SS mask (not all-RE), for each constraint mask over valid equiv groups,
yield a `MechanismSpec` with appropriate `eq_steps` and `param_constraints`.

**`_find_equivalent_groups(edges, adj)`** (~10 lines): Group binding edges
by (metabolite, site) using adjacency metadata.

### Section 8: Pipeline + Conversion (~27 lines)

**`enumerate_mechanisms(reaction; stage=Full(), max_forms=...)`** (~15 lines):
```julia
forms = enumerate_enzyme_forms(reaction)
adj = _build_adjacency(forms)

catalytic = _catalytic_topologies(forms, adj, reaction; max_forms)
stage isa Catalytic && return catalytic

with_act = _expand_activators(catalytic, adj, forms, reaction; max_forms)
stage isa WithActivator && return with_act

de_iter = Iterators.flatten(
    Iterators.map(s -> _dead_end_configs(s, adj, forms; max_forms), with_act))
stage isa WithDeadEnd && return de_iter

Iterators.flatten(
    Iterators.map(s -> _ress_variants(s, adj, forms), de_iter))
```

Default `max_forms`: `3 * (length(substrates(rxn)) + length(products(rxn)) + length(regulators(rxn)))` — inline, no `n_sites` function.

**`EnzymeMechanism(spec::MechanismSpec)`** (~12 lines):
Reconstruct forms + adjacency from `spec.reaction`. Convert edges to reaction
tuples using adjacency metadata. Build species/reactions/eq_steps/constraints.

## Line Budget

| Section | Lines |
|---|---|
| 1. Types + constants | 15 |
| 2. Edge classification + adjacency | 32 |
| 3. Form enumeration | 35 |
| 4. Catalytic cycles (unified DFS + combination) | 25 |
| 5. Activator configs | 22 |
| 6. Dead-end configs | 22 |
| 7. RE/SS lazy generator + equiv groups | 22 |
| 8. Pipeline + conversion | 27 |
| **Total** | **~200** |

## Files to Modify

1. **`src/mechanism_enumeration.jl`** — complete rewrite (1534 → ~200 lines)
2. **`src/EnzymeRates.jl`** — remove all mechanism enumeration exports:
   `SiteState`, `EnzymeFormSpec`, `MechanismSpec`, `enumerate_mechanisms`,
   `enumerate_mechanism_stages`, `n_sites`, `enumerate_enzyme_forms`
3. **`test/test_mechanism_enum_of_enz_reaction.jl`** — rewrite tests:
   - `enumerate_mechanism_stages(rxn)` → `enumerate_mechanisms(rxn; stage=X)`
   - `spec.reactions` → `spec.edges`
   - Use `collect()` + `length()` for counting at each stage
   - Drop `EnzymeRates._spec_to_edges` calls
4. **`test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`** —
   update test helpers:
   - `_spec_to_form_edges` → use `spec.edges` directly
   - `_find_equiv_groups` → update for adjacency-based edge metadata
   - `_compute_expected_dead_end_count` → update for new spec structure
   - `_compute_expected_n_total` → update for new spec structure

## Implementation Order

1. Write new `src/mechanism_enumeration.jl` (all 8 sections)
2. Remove old exports from `src/EnzymeRates.jl`
3. Adapt test files to new API
4. Run tests via MCP REPL, fix issues iteratively
5. Run full test suite (`Pkg.test()`) to verify all 8 test specs match
6. Verify: `wc -l src/mechanism_enumeration.jl` ≤ 210

## Verification

All 8 `ENUMERATION_TEST_SPECS` must pass with exact counts at every stage:

```julia
@testset "Pipeline: $(spec.name)" for spec in ENUMERATION_TEST_SPECS
    for (stage, expected) in [
        (Catalytic(), spec.expected_n_catalytic),
        (WithActivator(), spec.expected_n_cat_with_act),
        # WithDeadEnd and Full use collect()
    ]
        result = enumerate_mechanisms(spec.reaction;
            stage=stage, max_forms=spec.max_forms)
        @test length(result) == expected
    end
    # ...
end
```

Additional checks:
- `EnzymeMechanism(spec)` construction works for catalytic specs
- `rate_equation_string` works for converted mechanisms
- Lazy iteration doesn't OOM (verify with Bi-Bi + 2 regulators)
- `wc -l src/mechanism_enumeration.jl` ≤ 210

## Risks

1. **Unified DFS correctness**: Must reproduce exactly the same cycle counts
   as the current separate standard + PP enumerators. The current PP DFS was
   recently debugged (commits c255f23, c80052c). Mitigation: all 8 test specs
   must match exactly.
2. **DFS performance**: General DFS explores more paths than permutation-based
   standard cycle enumeration. For reactions with many core forms, this could
   be slower. Mitigation: `max_forms` limits cycle size; profile Bi-Bi + reg.
3. **Separation rule**: Verified only on Bi-Bi PP (10 topologies). If a future
   reaction produces valid mixed topologies, this would miss them. Mitigation:
   all 8 test specs pass; can add debug assertion if desired.
4. **Lazy stages 3-4**: `length()` not available on the final iterator.
   `IdentifyRateEquationProblem` will compute length separately.
5. **Line target**: 200 is tight. Edge classification + PP isomerization rules
   are inherently complex. May reach 210-220 if readability demands more space.
