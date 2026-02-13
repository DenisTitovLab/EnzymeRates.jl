# Plan: Enumerate All Valid EnzymeMechanisms from EnzymeFormSpecs

## Context

Given an `EnzymeReaction`, `enumerate_enzyme_forms` produces all possible enzyme form specs (e.g., ~23 core forms for hexokinase). We need a new function `enumerate_mechanisms` that produces all valid `EnzymeMechanism` topologies (with SS/RE assignments) as lightweight tuples. The user can then materialize specific mechanisms on demand.

**Key constraint**: Return tuples, not `EnzymeMechanism` types, to avoid compiling 10^6 distinct types.

## Files to Modify

1. **`src/mechanism_enumeration.jl`** — Add all new functions (~300 lines)
2. **`src/EnzymeRates.jl`** — Add exports
3. **`test/test_mechanism_enumeration.jl`** — Add tests

## Step 1: Fix `enumerate_enzyme_forms` sites vector reuse bug

`src/mechanism_enumeration.jl:129` — The `sites` vector is allocated once (line 108) and reused across iterations. All `EnzymeFormSpec` objects share the same vector reference.

Fix: `push!(forms, EnzymeFormSpec(..., copy(sites)))`.

## Step 2: Add internal function `_enumerate_forms_with_atoms`

Returns `Vector{Tuple{Symbol, Tuple{Vararg{Tuple{Symbol,Int}}}}}` — form name + total atom content as a tuple suitable for `EnzymeMechanism` species format.

Computed inline during the same Cartesian product loop as `enumerate_enzyme_forms`, accumulating `total_atoms = Dict{Symbol,Int}` from site contents per form, then converting to sorted tuple.

Also returns per-form `site_contents` (which sites are occupied) and `is_regulator_site` flags — needed for the core/regulator decomposition.

## Step 3: Compute all elementary edges

Function `_compute_edges(form_atoms, met_atoms_dict) → Vector{NamedTuple}`

For each pair of forms `(i, j)` where `i < j`, compute atom delta `atoms_j - atoms_i`:
- **Empty delta**: isomerization (same total atoms, different site arrangement). Edge: `(i, j, metabolite=nothing)`.
- **Delta matches a known metabolite M** (all entries positive): M binds going `i → j`. Edge: `(i, j, metabolite=M, binding_from=i)`.
- **Negative delta matches metabolite M**: M binds going `j → i`. Edge: `(i, j, metabolite=M, binding_from=j)`.
- **Otherwise**: not an elementary step, skip.

Met atoms dict includes substrates, products, AND regulators. Built from `EnzymeReaction` type params.

## Step 4: Topology enumeration — grow-and-forbid DFS

Core algorithm. Enumerates all connected edge subsets containing free enzyme that satisfy stoichiometry parity constraints.

### Stoichiometry constraints per metabolite M with `n_M` edges:
- Substrate (target = -1): `n_M` odd, `n_M ≥ 1`
- Product (target = +1): `n_M` odd, `n_M ≥ 1`
- Regulator (target = 0): `n_M` even (0 = regulator not used)

### DFS algorithm:
```
grow(start_edge_idx):
    if stoichiometry_valid(met_counts): record topology
    if length(included) >= max_steps: return

    for e_idx in start_edge_idx:n_edges:
        edge = edges[e_idx]
        # Only add edges incident to current subgraph (guarantees connectivity)
        if neither endpoint in current forms: continue

        # Add edge, update form_refcount and met_counts
        push!(included, e_idx)

        # Pruning: can remaining edges (idx > e_idx) fix parity for all metabolites?
        if can_complete(): grow(e_idx + 1)

        # Backtrack
        pop!(included)
```

Start: `forms = {free_enzyme}`, `included = []`, `grow(1)`.

### Pruning `_can_complete`:
For each metabolite M with current count `c_M`:
- Count remaining available edges for M (index > current)
- If `c_M` has wrong parity and no remaining edges can fix it → prune
- If M is substrate/product with `c_M == 0` and no remaining edges → prune

## Step 5: Assign edge orientations

For each valid topology, determine reaction directions. Each metabolite edge can be oriented as binding (metabolite on LHS, contributes -1 to net) or release (metabolite on RHS, contributes +1).

For each metabolite M with `n_M` edges and target `t_M`:
- `n_binding = (n_M - t_M) / 2` (guaranteed integer by parity check)
- Choose WHICH `n_binding` edges are binding: `C(n_M, n_binding)` choices
- All choices produce equivalent rate equations → pick canonical: first `n_binding` edges in index order as binding

For isomerization: canonical direction = lower-index form on LHS.

Build `(lhs, rhs)` tuples:
- Binding: `((enzyme_form, metabolite), (enzyme_form_with_met,))`
- Release: `((enzyme_form_with_met,), (enzyme_form, metabolite))`
- Isomerization: `((form_lo,), (form_hi,))`

## Step 6: Enumerate SS/RE assignments

For each oriented topology with `n` steps: iterate `mask` from `0` to `2^n - 2`:
- `eq_steps[i] = ((mask >> (i-1)) & 1) == 1`
- mask = 0: all SS; mask = 2^n - 1: all RE (skipped)

## Step 7: Build output specs

Each mechanism spec is a `NamedTuple{(:species, :reactions, :eq_steps)}`:
- `species = (subs_tuple, prods_tuple, regs_tuple, enzymes_tuple)` matching `EnzymeMechanism` format
  - `subs_tuple/prods_tuple` from reaction: `((name, atoms), ...)`
  - `regs_tuple`: only include regulators that appear in mechanism edges
  - `enzymes_tuple`: `((form_name, form_atoms), ...)` for forms used in this mechanism
- `reactions`: tuple of `(lhs, rhs)` from step 5, sorted canonically
- `eq_steps`: tuple of `Bool` from step 6

Deduplicate via `Set` on the full spec tuple.

## Step 8: Public API

```julia
function enumerate_mechanisms(reaction::EnzymeReaction; max_steps::Union{Nothing,Int}=nothing)
    # Function barrier: extract type params to runtime values
    _enumerate_mechanisms_impl(reaction; max_steps)
end

@nospecialize function _enumerate_mechanisms_impl(reaction; max_steps)
    # All enumeration logic here
end

function materialize_mechanism(spec::NamedTuple)
    EnzymeMechanism(spec.species, spec.reactions, spec.eq_steps)
end
```

## Step 9: Exports

Add to `src/EnzymeRates.jl`:
```julia
export enumerate_mechanisms, materialize_mechanism
```

## Step 10: Tests

Add to `test/test_mechanism_enumeration.jl`:

1. **Uni-Uni hand-verification**: 3 forms, 3 edges → enumerate all mechanisms, verify count matches hand calculation, verify all materialize without error
2. **Bi-Bi with max_steps**: verify mechanisms found, all materialize
3. **max_steps monotonicity**: more steps allowed → more mechanisms
4. **No duplicates**: all specs unique
5. **Rate equations derivable**: materialized mechanisms produce valid rate equation strings
6. **Regulator handling**: mechanisms with regulators have balanced edges and materialize correctly
7. **Compilation time**: subprocess test for hexokinase, cold < 5s, warm < 2s (with max_steps=6 or similar)

## Verification

After implementation:
1. `julia --project -e 'using Pkg; Pkg.test()'` — all tests pass
2. Manual check: enumerate Uni-Uni mechanisms and verify they match known textbook mechanisms
3. Hexokinase enumeration completes in reasonable time with count < 10^6

## Key Design Decisions

### Regulators
The `EnzymeMechanism` constructor requires net stoichiometry = 0 for regulators. For dead-end inhibitor complexes, this means regulator edges must come in balanced pairs (binding + release orientations). Mechanisms can also simply exclude regulator-bound forms (no inhibition modeled). The enumeration handles both cases: regulator edges must have even count (0 = not used, 2+ = balanced pairs).

### Orientation Equivalence
Different orientation assignments for the same physical edge set produce the same rate equation (swapping forward/reverse rate constants). Therefore, we pick ONE canonical orientation per topology rather than enumerating all equivalent orientations.

### Compilation Strategy
- Return lightweight `NamedTuple` specs, not `EnzymeMechanism` types
- Use `@nospecialize` function barrier for enumeration implementation
- `materialize_mechanism(spec)` creates the actual type on demand
- Follows existing pattern from `enumerate_enzyme_forms`

### Algorithm Complexity
For hexokinase (~23 core forms, ~57 core edges):
- DFS with stoichiometry pruning is the core algorithm
- `max_steps` parameter limits mechanism size and search depth
- Core/regulator decomposition can be added as optimization if needed
