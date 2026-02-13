# Implementation Plan: `enumerate_mechanisms`

## Summary of Design Decisions (from interview)

| Decision | Choice |
|----------|--------|
| Topology algorithm | Cycle-first: find all valid 1× simple cycles through E, then combine subsets |
| Lazy vs eager | Eager topology enumeration, lazy Cartesian product over (topology × dead-ends × RE/SS × constraints) |
| Regulator pruning | No pruning — all forms enter reaction graph, stoichiometry constraints filter naturally |
| Cycle combining | Shared forms & edges allowed between cycles |
| Stoichiometry check | Post-hoc validation (find cycles first, then filter by stoichiometry) |
| Emergent cycle validation | Enumerate all simple cycles in the (small) union graph via DFS |
| Dead-end limits | Global max_forms cap only, no per-form depth cap |
| Neighbor engine | Fully independent implementation |
| MechanismSpec shape | Flat runtime struct with form names + reaction tuples; convert to nested tuples on EnzymeMechanism() |
| Cycle direction | Directed from start (each undirected edge → 2 directed edges) |
| Edge direction model | Each undirected reaction yields 2 directed edges; traversal direction determines consumed/produced metabolites |
| Topology dedup | Sorted form set + edge set hash |
| Dead-end coupling | Accept topology-dependence; recompute dead-end lattice per topology |
| Ping-pong edges | Generalized release: site loses atoms, difference matched to known metabolite by exact atom difference |
| Equivalent steps | Same metabolite + same site index + different enzyme state (renamed from "symmetric steps") |
| Validation on convert | Always validate via existing EnzymeMechanism constructor |
| Cycle combo limit | No limit on number of cycles combined |
| max_forms default | 3 × n_sites |
| MechanismSpec display | No display — data struct only |
| File structure | Single file: src/mechanism_enumeration.jl |
| Test file | Extend existing test/test_mechanism_enumeration.jl |
| Forms bug | Fix copy(sites) bug as part of this PR |

---

## Phase 0: Bug Fix

**Fix the shared `sites` vector bug in `enumerate_enzyme_forms`** (line 129 of `src/mechanism_enumeration.jl`).

Currently `push!(forms, EnzymeFormSpec(..., sites))` reuses the same `sites` vector for all forms. Fix: `copy(sites)` in the push call.

---

## Phase 1: Reaction Graph Construction

Build a directed reaction graph from the output of `enumerate_enzyme_forms`.

### Data Structures

```julia
struct ReactionEdge
    from::Int          # index into forms vector
    to::Int            # index into forms vector
    metabolite::Union{Nothing, Symbol}  # nothing for isomerization
    edge_type::Symbol  # :binding, :release, :isomerization
end
```

The reaction graph is stored as `Vector{ReactionEdge}` plus the `Vector{EnzymeFormSpec}` from Phase 0.

### Edge Construction Algorithm

For every pair of forms (F1, F2) where F1 ≠ F2, check if a valid elementary reaction exists. Since each undirected reaction produces 2 directed edges, we process each unordered pair once and emit both directions.

**Binding edge** (F1 + M → F2): F1 and F2 differ at exactly one site `i`. Site `i` is unoccupied (`nothing`) in F1 and occupied (has atoms) in F2. All other sites identical. The metabolite M is identified by `F2.sites[i].metabolite`. Emit:
- Forward: `ReactionEdge(f1_idx, f2_idx, M, :binding)`
- Reverse: `ReactionEdge(f2_idx, f1_idx, M, :release)`

**Generalized release edge** (ping-pong): F1 and F2 differ at exactly one site `i`. Site `i` is occupied in both F1 and F2, but F1 has MORE atoms than F2 (F2 could be `nothing` for full release, or partial atoms for ping-pong). The difference `F1.atoms - F2.atoms` must exactly match some known metabolite M's atoms. Emit:
- Forward: `ReactionEdge(f1_idx, f2_idx, M, :release)`
- Reverse: `ReactionEdge(f2_idx, f1_idx, M, :binding)`

Note: The "standard" binding/release case (site goes empty↔occupied) is a special case of this. We can unify: for every pair differing at exactly one site, compute the atom difference and try to match a metabolite.

**Isomerization edge** (F1 → F2): F1 and F2 differ only at core catalytic sites (index-1 sites of substrates and products). All non-catalytic sites (extra sites with index ≥ 2, regulator sites) must be identical. Total atoms across all differing core sites must be equal on both sides. F1 ≠ F2. Both directions are valid:
- Forward: `ReactionEdge(f1_idx, f2_idx, nothing, :isomerization)`
- Reverse: `ReactionEdge(f2_idx, f1_idx, nothing, :isomerization)`

**Spec update**: The isomerization rule explicitly includes regulator sites as non-catalytic sites that must be identical between F1 and F2. This naturally allows isomerization between regulator-bound forms (e.g., EAS → EAP where regulator A is bound identically in both).

### Complexity

O(n² × s) where n = number of forms, s = number of sites per form. For glucokinase with ~128 forms and 7 sites, this is ~115k pair checks — fast.

### Helper: Metabolite Atom Lookup

Build a `Dict{Vector{Pair{Symbol,Int}}, Symbol}` mapping sorted atom vectors to metabolite names. Used to match atom differences to metabolites during release edge detection.

---

## Phase 2: Cycle Enumeration

Find all valid simple directed cycles through the free enzyme E in the reaction graph.

### Algorithm: DFS from E

1. Start DFS from the free enzyme form index.
2. Maintain a visited set and current path (as a vector of `(form_idx, edge_idx)` pairs).
3. At each node, explore all outgoing directed edges to unvisited nodes (or back to E if path length ≥ 2).
4. When we return to E, record the cycle (the sequence of directed edges).
5. Backtrack and continue.

### Stoichiometry Validation

For each found cycle, compute net stoichiometry:
- For each edge in the cycle:
  - `:binding` with metabolite M → consumed: M (stoich −1)
  - `:release` with metabolite M → produced: M (stoich +1)
  - `:isomerization` → no metabolite change (stoich 0)
- Sum up per-metabolite stoichiometry across all edges.
- Classify:
  - **1× cycle**: net stoichiometry matches −1 for each substrate, +1 for each product, 0 for each regulator
  - **0× cycle**: net stoichiometry is 0 for all metabolites (futile)
  - **Invalid**: anything else → discard

Keep only **1× cycles**. (0× cycles are emergent from topology combinations, not enumerated directly.)

### Output

`Vector{Vector{ReactionEdge}}` — each inner vector is a directed cycle of edges.

---

## Phase 3: Topology Combination

Combine subsets of 1× cycles into valid topologies.

### Algorithm

1. Start with the set of all 1× cycles from Phase 2.
2. Enumerate all non-empty subsets of cycles (2^k − 1 subsets for k cycles).
3. For each subset:
   a. Compute the **union topology**: set of all forms and edges across all cycles in the subset.
   b. Check form count ≤ max_forms.
   c. **Validate emergent cycles**: enumerate ALL simple cycles through E in the union graph (DFS). Each must be 1× or 0×. If any cycle has invalid stoichiometry (2×, −1×, etc.), reject.
   d. Canonicalize: sorted (form_set, edge_set) tuple.
   e. Deduplicate against previously seen canonical forms.
4. Store valid, unique topologies.

### Optimization: Incremental Combination

Instead of enumerating all 2^k subsets (which could be huge for large k):
- Sort cycles by form count (smallest first).
- Build topologies incrementally: start with single cycles, then try adding each additional cycle.
- Use BFS/DFS over the "add a cycle" operation.
- Prune: if adding a cycle exceeds max_forms, skip.
- Prune: if the union already has an invalid emergent cycle, skip further additions.
- Dedup at each level using the canonical hash.

This avoids enumerating all 2^k subsets when most are invalid.

### Emergent Cycle Validation (DFS)

For a union graph with ~15 forms, enumerate all simple cycles through E via DFS. For each cycle, compute stoichiometry and check it's 1× or 0×. This is fast for small graphs.

### Output

`Vector{Topology}` where:

```julia
struct Topology
    form_indices::Vector{Int}       # indices into forms vector
    edges::Vector{ReactionEdge}     # directed edges in this topology
end
```

---

## Phase 4: Dead-End Lattice Enumeration

For each topology, compute the dead-end attachment options.

### Per-Cycle-Form Dead-End Computation

For each form F in the topology:
1. Find all forms G in the full forms list such that:
   - G differs from F by exactly one additional site occupied (any metabolite)
   - G is NOT in the topology's form set
   - The binding edge F→G (or G→F as release) exists in the full reaction graph
2. These are the **direct dead-end children** of F.
3. Recursively, each dead-end child can have its own children (dead-end chains): forms that differ by one more binding step, also not in the topology.
4. The result is a **dead-end tree** rooted at F (actually a DAG if multiple paths lead to the same dead-end).

### Downward-Closed Subsets

A valid dead-end configuration for F is any downward-closed subset of F's dead-end lattice:
- ∅ (no dead-ends for this form)
- Any single direct child
- Any direct child + its children (if we include a grandchild, we must include the intermediate)
- Etc.

The count of downward-closed subsets of a DAG is computed via DP. For a tree of depth d with branching factor b, the count is product over nodes of (1 + count_of_subtree_configs).

**Simplification**: In practice, dead-end lattices are small (typically 1-3 levels deep, 1-5 branches). The downward-closed subset count per form is manageable.

### Total Dead-End Multiplier

```
dead_end_multiplier(topology) = product(n_dead_end_configs(F) for F in topology.forms)
```

This is precomputed per topology and stored alongside it.

### max_forms Enforcement

Dead-end forms count toward the global max_forms limit. Since dead-end options across forms are independent, we need to enforce this constraint during iteration (when materializing a specific dead-end combination), not during counting.

**Approach**: For each topology, compute the dead-end lattice per form. During iteration, when selecting a specific dead-end configuration (via modular arithmetic indexing), check that total forms (cycle + dead-ends) ≤ max_forms. Skip configurations that exceed the limit.

Alternative: Precompute only configurations that satisfy max_forms. This is more complex but avoids wasted iteration steps. Given that most topologies are small, the skip approach is simpler and sufficient.

---

## Phase 5: RE/SS Assignment Enumeration

For a topology with N edges (steps), enumerate all 2^N − 1 valid RE/SS assignments (excluding all-RE).

This is a simple bit-vector enumeration. For each topology, N is typically 3-10, giving 7 to 1023 assignments.

### Indexing

For a given assignment index `i` (1 to 2^N − 1), convert to a Bool tuple:
- bit j of i → step j is RE (true) if bit is 1, SS (false) if bit is 0
- Skip i = 2^N − 1 (all bits set = all RE)

---

## Phase 6: Equivalent Step Detection & Constraint Enumeration

### Detection

Two steps are **equivalent** if:
- Both are binding reactions (or both release — but since every undirected edge has both, we focus on the binding direction)
- Both bind the same metabolite M
- Both bind at the same site index (e.g., site 1 of metabolite I)
- They differ only in the enzyme state they bind to (e.g., E+I→EI vs ES+I→ESI)

Group all such steps into equivalence classes. Each class with ≥ 2 members generates a constraint option.

### Constraint Variants

For M equivalence groups, generate 2^M variants:
- For each group, independently choose constrained or unconstrained.
- Constrained: all steps in the group share the same rate/equilibrium parameters.

### Constraint Representation

For a group of equivalent steps {step_i, step_j, ...}:
- If all steps are RE: constraint is `K_j = K_i` (all equal to the first)
- If all steps are SS: constraints are `k_jf = k_if, k_jr = k_ir`
- Mixed RE/SS within a group: the constraint only applies between steps of the same type (or we skip constraint for mixed groups)

**Note**: Equivalent step constraints interact with RE/SS assignments. A group of steps may be all-RE in one assignment and mixed in another. The constraint enumeration must be done per (topology × RE/SS assignment) combination, not independently.

**Revised factoring**: Constraint options depend on RE/SS assignment. So the iteration order is:
1. For each topology:
2.   For each RE/SS assignment:
3.     Detect equivalent step groups (considering the current RE/SS assignment)
4.     For each constraint variant (2^M options):
5.       For each dead-end configuration:
6.         Yield MechanismSpec

This means we can't fully factor constraints independently. But since equivalent step detection is cheap (scan edges for same metabolite + same site index), this is fine.

---

## Phase 7: MechanismSpec Struct & Iterator

### MechanismSpec

```julia
struct MechanismSpec
    reaction::Any                          # the source EnzymeReaction (untyped to avoid specialization)
    forms::Vector{Symbol}                  # enzyme form names in this mechanism
    form_atoms::Vector{Vector{Pair{Symbol,Int}}}  # atoms for each form
    reactions::Vector{Tuple{Vector{Symbol}, Vector{Symbol}}}  # (lhs, rhs) per step
    equilibrium_steps::Vector{Bool}        # RE/SS per step
    param_constraints::Vector{Tuple{Symbol, Int, Vector{Tuple{Symbol,Int}}}}  # (target, coeff, factors)
end
```

### MechanismIterator

```julia
struct MechanismIterator
    # Precomputed data (from reaction graph + forms)
    all_forms::Vector{EnzymeFormSpec}
    reaction_graph::Vector{ReactionEdge}
    metabolite_atoms::Dict{Symbol, Vector{Pair{Symbol,Int}}}
    reaction::Any  # the source EnzymeReaction

    # Eagerly enumerated topologies
    topologies::Vector{Topology}

    # Per-topology precomputed data (dead-end lattices, lazily computed)
    # Computed on first access to length() or during iteration
end
```

### Iterator Protocol

```julia
function Base.iterate(iter::MechanismIterator, state=nothing)
    # state tracks: (topology_idx, ress_idx, constraint_idx, deadend_idx)
    # Uses modular arithmetic to advance through the Cartesian product
    # Skips combinations where total forms > max_forms
    # Materializes MechanismSpec for each valid combination
end

function Base.length(iter::MechanismIterator)
    # Sum over topologies of:
    #   sum over RE/SS assignments of:
    #     constraint_options(topology, ress) × dead_end_options(topology)
    # Dead-end options are topology-dependent, computed lazily and cached
end

Base.eltype(::Type{MechanismIterator}) = MechanismSpec
```

---

## Phase 8: EnzymeMechanism(spec::MechanismSpec) Conversion

### Conversion Logic

1. Build the `species` tuple from the MechanismSpec:
   - Substrates, products, regulators: from `spec.reaction`
   - Enzyme forms: from `spec.forms` + `spec.form_atoms`
2. Build the `reactions` tuple from `spec.reactions`.
3. Build the `equilibrium_steps` tuple from `spec.equilibrium_steps`.
4. Build the `param_constraints` tuple from `spec.param_constraints`.
5. Call `EnzymeMechanism(species, reactions, eq_steps, constraints)` — uses the existing constructor with full validation.

### Location

Add to `src/types.jl` as a new constructor method:

```julia
function EnzymeMechanism(spec::MechanismSpec)
    # Convert flat vectors to nested tuples
    # Call existing constructor
end
```

---

## Phase 9: Exports & Integration

### src/EnzymeRates.jl

Add exports:
```julia
export enumerate_mechanisms, MechanismSpec
```

### Update CLAUDE.md

Add `src/mechanism_enumeration.jl` description update noting mechanism topology enumeration.

---

## Phase 10: Tests

All tests go in `test/test_mechanism_enumeration.jl`, extending the existing file.

### Test 1: Uni-Uni Topology Count
- 2 topologies: 2-step (E→ES→E+P, no isomerization) and 3-step (E→ES→EP→E+P)
- Per topology: RE/SS options (2^2−1=3 for 2-step, 2^3−1=7 for 3-step)
- No regulators → no dead-ends, no equivalent steps
- Total: 3 + 7 = 10 mechanisms
- Verify `length(iter) == 10`

### Test 2: Bi-Bi Topologies
- Must include ordered Bi-Bi and random-order Bi-Bi
- Verify all specs convert to valid EnzymeMechanism

### Test 3: Regulators — Dead-End Inhibitor
- Uni-Uni + regulator I: dead-end complexes EI, ESI, etc.
- Verify mechanisms with and without dead-ends are present

### Test 3b: Activator Mechanism
- Uni-Uni + regulator A: must include topology where A is in the cycle
- E→EA→EAS→EA→E is a valid 1× cycle (A binds and releases, S→P)
- Verify this topology exists in the output

### Test 3c: Product Inhibition / Abortive Complexes
- Bi-Bi: must include dead-end forms like EP, EAQ
- Verify these attach to cycle forms via valid reactions

### Test 4: All Mechanisms Constructible
- For Uni-Uni: convert every MechanismSpec to EnzymeMechanism, verify no errors

### Test 5: Random-Order Stoichiometry
- Random-order Bi-Bi: verify per-cycle stoichiometry is 1× even though total edges > minimal cycle

### Test 6: Count Verification
- Manually count expected mechanisms for Uni-Uni, verify length(iter) matches

### Test 7: Compilation Guard
- Glucokinase (7 sites + regulators): cold < 30s, warm < 5s
- mechanism count < 10^6

### Test 8: max_forms Limit
- Verify smaller max_forms → fewer or equal mechanisms
- Verify all forms count ≤ max_forms

### Test 9: Equivalent Step Constraints
- Uni-Uni + regulator with same metabolite binding at same site index to E and ES
- Verify both constrained and unconstrained variants exist

---

## Implementation Order

1. **Phase 0**: Fix `copy(sites)` bug
2. **Phase 1**: Reaction graph construction + unit tests for edge detection
3. **Phase 2**: Cycle enumeration (DFS) + stoichiometry filter
4. **Phase 3**: Topology combination (incremental, with dedup and validation)
5. **Phase 4**: Dead-end lattice computation
6. **Phase 5**: RE/SS bit-vector enumeration
7. **Phase 6**: Equivalent step detection + constraint variants
8. **Phase 7**: MechanismSpec struct + MechanismIterator with iterate/length
9. **Phase 8**: EnzymeMechanism(spec) conversion in types.jl
10. **Phase 9**: Exports
11. **Phase 10**: Full test suite

Phases 1-3 are the algorithmic core and highest risk. Phases 4-6 are combinatorial but straightforward. Phases 7-9 are plumbing. Phase 10 validates everything.

---

## Risk Areas

1. **Cycle explosion**: For large reaction graphs, the number of simple cycles through E could be very large. Mitigation: incremental topology building with max_forms pruning eliminates most combinations early.

2. **Emergent cycle validation**: Enumerating all simple cycles in a union graph is worst-case exponential. Mitigation: enzyme mechanism graphs are sparse and small (typically < 15 forms in a topology).

3. **Ping-pong edge detection**: Matching atom differences to metabolites requires exact atom arithmetic. Must handle edge cases where multiple metabolites have the same atom content (unlikely but possible with degenerate reactions).

4. **Dead-end / max_forms interaction**: Skipping dead-end configurations that exceed max_forms during iteration may waste index space (length() overcounts). Alternative: precompute valid configurations. Decision: accept overcounting for simplicity, document that length() is an upper bound.

Actually, on reflection, `length()` should be exact. We should precompute the valid dead-end configuration count per topology respecting max_forms. This means enumerating downward-closed subsets up to a size budget. This is more complex but gives an accurate count.

**Revised approach for dead-ends + max_forms**: For each topology with C cycle forms, compute the dead-end lattice per form. Then enumerate valid combinations where total dead-end forms across all cycle forms ≤ (max_forms - C). This is a constrained Cartesian product — the counts per form are no longer independent. We'll need a DP or enumeration approach that respects the global budget.

**Simplification**: If the total possible dead-ends per form is small (typically 0-5), the constrained enumeration is fast. Use a recursive enumeration that tracks remaining budget as it assigns dead-end configs to each form.
