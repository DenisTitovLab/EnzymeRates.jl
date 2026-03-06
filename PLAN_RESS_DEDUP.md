# Implementation Plan: RE/SS Deduplication (Phase 1 + Phase 2)

## Overview

**Goal:** Reduce redundant mechanism enumeration by deduplicating RE/SS variants that produce observationally indistinguishable rate equations.

**Files to modify:**
- `src/mechanism_enumeration.jl` — core enumeration logic
- `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl` — test helpers + expected counts
- `test/test_mechanism_enum_of_enz_reaction.jl` — new validation tests

**Key constraint:** `_count_ress_variants` (production) and `_compute_expected_n_total` (test helper) must implement identical dedup logic, since tests verify they agree.

## Background

The King-Altman/Cha rate equation `v(concs) = E_total * N(concs) / D(concs)` is a rational function of metabolite concentrations. Its **monomial fingerprint** — which concentration monomials (like 1, S, P, S*P, A*B, etc.) appear in N and D — is fully determined by the RE group partition combined with the metabolite labeling of edges.

Two RE/SS assignments that yield the same monomial fingerprint produce **observationally indistinguishable** rate equations from steady-state data. The one with more SS steps just has extra non-identifiable microscopic parameters (more unknowns mapping to the same apparent kinetic coefficients), making fitting harder with zero gain in model expressiveness.

### Phase 1: Partition-based dedup
Skip SS masks that produce an RE group partition already seen (cheap hash check). Catches cases where cycles in the enzyme form graph provide redundant RE paths.

### Phase 2: Fingerprint-based dedup
Skip SS masks whose concentration-monomial fingerprint matches one already seen (lightweight set computation). Strictly stronger than Phase 1 — catches cases where different partitions yield the same polynomial form (e.g., all partitions of a linear mechanism).

Both phases keep the **first** (fewest SS steps) variant for each duplicate group, guaranteed by iterating masks in order of ascending `count_ones`.

---

## Step 1 — Partition helper + Phase 1 dedup

**Can be done by one agent. All changes are local to the two enumeration files.**

### 1.1 Add `_compute_re_partition` in `src/mechanism_enumeration.jl`

Place after `_compute_re_group_count` (~line 492). Same union-find logic but returns the actual partition instead of just the count:

```julia
"""
    _compute_re_partition(edges, eq_steps) -> Vector{Vector{Int}}

Compute RE group partition: connected components when only RE edges
are considered. Returns sorted vector of sorted form-index vectors.
"""
function _compute_re_partition(edges, eq_steps)
    form_indices = collect(Set(Iterators.flatten(edges)))
    parent = Dict(i => i for i in form_indices)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    for (idx, (a, b)) in enumerate(edges)
        eq_steps[idx] || continue
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end
    groups = Dict{Int, Vector{Int}}()
    for i in form_indices
        r = find(i)
        push!(get!(groups, r, Int[]), i)
    end
    sort!([sort!(v) for v in values(groups)])
end
```

The length of the returned vector equals G. Keep `_compute_re_group_count` for backward compat where only G is needed.

Also add a matching `_compute_re_partition_test` in the test file (independent reimplementation, mirrors the existing `_compute_re_group_count_test` pattern).

### 1.2 Modify `_count_ress_variants`

Current code (lines 596-624) iterates `ss_mask in 0:(1 << n_other) - 1`. Change to:

1. **Sort masks by popcount**: `sort!(collect(0:(1 << n_other) - 1); by=count_ones)`. This ensures the first mask producing each partition has the fewest SS steps (most parsimonious).
2. **Track seen partitions**: `seen = Set{Vector{Vector{Int}}}()`.
3. After computing `eq_steps` and checking G bounds, compute `partition = _compute_re_partition(edges, eq_steps)`.
4. `partition in seen && continue`; otherwise `push!(seen, partition)`.
5. Rest of counting logic unchanged.

```julia
function _count_ress_variants(spec, adj, forms; max_re_groups::Int=7)
    edges = spec.edges
    n = length(edges)
    n == 0 && return 0
    n_cat = spec.n_catalytic_edges
    iso_idx = _find_first_isomerization(edges, adj)
    equiv_groups = _find_equivalent_groups(edges, adj, forms, n_cat)
    other_indices = [i for i in 1:n_cat if i != iso_idx]
    n_other = length(other_indices)
    masks = sort!(collect(0:(1 << n_other) - 1); by=count_ones)
    seen = Set{Vector{Vector{Int}}}()
    total = 0
    for ss_mask in masks
        eq_steps = fill(true, n)
        eq_steps[iso_idx] = false
        for (bit, idx) in enumerate(other_indices)
            (ss_mask >> (bit - 1)) & 1 == 1 && (eq_steps[idx] = false)
        end
        G = _compute_re_group_count(edges, eq_steps)
        (G < 2 || G > max_re_groups) && continue
        partition = _compute_re_partition(edges, eq_steps)
        partition in seen && continue
        push!(seen, partition)
        valid_groups = [g for g in equiv_groups
            if all(eq_steps[s] == eq_steps[g[1]] for s in g)]
        total += 1 << length(valid_groups)
    end
    total
end
```

### 1.3 Modify `_ress_variants`

Same changes as 1.2 but for the lazy iterator version (lines 550-588). The `seen` set is captured by the closure in `Iterators.flatmap` and populated as the iterator is consumed:

```julia
function _ress_variants(spec, adj, forms; max_re_groups::Int=7)
    edges = spec.edges
    n = length(edges)
    n_cat = spec.n_catalytic_edges
    iso_idx = _find_first_isomerization(edges, adj)
    equiv_groups = _find_equivalent_groups(edges, adj, forms, n_cat)
    other_indices = [i for i in 1:n_cat if i != iso_idx]
    n_other = length(other_indices)
    masks = sort!(collect(0:(1 << n_other) - 1); by=count_ones)
    seen = Set{Vector{Vector{Int}}}()
    Iterators.flatmap(masks) do ss_mask
        eq_steps = fill(true, n)
        eq_steps[iso_idx] = false
        for (bit, idx) in enumerate(other_indices)
            (ss_mask >> (bit - 1)) & 1 == 1 && (eq_steps[idx] = false)
        end
        G = _compute_re_group_count(edges, eq_steps)
        (G < 2 || G > max_re_groups) && return ()
        partition = _compute_re_partition(edges, eq_steps)
        partition in seen && return ()
        push!(seen, partition)
        # ... generate constraint variants (unchanged) ...
    end
end
```

### 1.4 Update test helper `_compute_expected_n_total`

In `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl` (lines 73-102), apply identical changes:

1. Add a test-local `_compute_re_partition_test` (independent reimplementation, like the existing `_compute_re_group_count_test`).
2. Sort masks by popcount.
3. Track seen partitions, skip duplicates.

### 1.5 Compute new expected counts and update test specs

Run the modified code to get new `expected_n_total` values for all 7 test specs:

```bash
julia --project -e '
    using EnzymeRates
    include("test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl")
    for spec in ENUMERATION_TEST_SPECS
        final = EnzymeRates.enumerate_mechanisms(spec.reaction)
        println("$(spec.name): n_total=$(length(final))")
    end
'
```

Update `expected_n_total` in each `EnumerationTestSpec`. Note: `expected_n_forms`, `expected_n_catalytic`, `expected_n_cat_de` are UNCHANGED (dedup only affects the RE/SS stage).

### 1.6 Run full test suite

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

## Step 2 — Concentration fingerprint function (Phase 2 foundation)

**Can be done by a separate agent in parallel with Step 1.** This step only ADDS a new function — it doesn't modify any existing code.

### 2.1 Add `_concentration_fingerprint` in `src/mechanism_enumeration.jl`

Place after Phase 1 functions. This computes the set of metabolite-only monomials that would appear in the rate equation denominator, using only the topology (edges, eq_steps, forms, adj) — no full King-Altman derivation.

**Type:** A concentration monomial is `Vector{Pair{Symbol,Int}}` (sorted by symbol), representing a product of metabolite concentrations like `A^1 * B^1` = `[:A => 1, :B => 1]`. The fingerprint is a `Set` of these.

**Algorithm (5 sub-steps):**

#### 2.1a — Compute RE partition and form-to-group mapping

```julia
partition = _compute_re_partition(edges, eq_steps)
G = length(partition)
group_set = [Set(g) for g in partition]
form_to_group = Dict(i => g for (g, grp) in enumerate(partition) for i in grp)
```

#### 2.1b — Compute alpha concentration monomials per form

For each RE group, BFS from the reference form (first in sorted group) through RE edges. At each step, if the edge involves a metabolite binding, add that metabolite to the concentration monomial.

```julia
alpha_conc = Dict{Int, Vector{Pair{Symbol,Int}}}()
for (g, group) in enumerate(partition)
    ref = group[1]
    alpha_conc[ref] = Pair{Symbol,Int}[]  # empty = constant 1
    queue = [ref]
    while !isempty(queue)
        cur = popfirst!(queue)
        for (idx, (a, b)) in enumerate(edges)
            eq_steps[idx] || continue
            neighbor = (a == cur && b in group_set[g]) ? b :
                       (b == cur && a in group_set[g]) ? a : nothing
            (neighbor === nothing || haskey(alpha_conc, neighbor)) && continue
            info = adj[minmax(a, b)]
            parent_mono = alpha_conc[cur]
            child_mono = if info.metabolite !== nothing
                is_binding_direction = (info.type == :binding && a == cur) ||
                                      (info.type == :release && b == cur)
                is_binding_direction ? _add_met(parent_mono, info.metabolite) : parent_mono
            else
                copy(parent_mono)
            end
            alpha_conc[neighbor] = child_mono
            push!(queue, neighbor)
        end
    end
end
```

Helper `_add_met(mono, met)` returns a new sorted `Vector{Pair{Symbol,Int}}` with `met`'s exponent incremented by 1.

#### 2.1c — Compute sigma concentration monomials per group

```julia
sigma_conc = [Set(alpha_conc[i] for i in group) for (_, group) in enumerate(partition)]
```

#### 2.1d — Compute rate matrix concentration monomials R_conc[g1,g2]

For each SS step connecting form `i` (group g1) to form `j` (group g2), the concentration monomials come from the metabolites on each side of the step combined with the alpha of the source form.

Edge classification direction: the adjacency `adj` stores `(a,b)` with `a < b`. The `info.type` is `:binding` if going `a->b` adds a metabolite, `:release` if going `a->b` removes one. Cross-reference with `_classify_edge` (lines 143-194).

For the rate matrix:
- Binding edge `E + S <=> ES` (type=:binding, met=S): forward `a->b` has S as concentration, reverse `b->a` has no metabolite
- Release edge `EP <=> E + P` (type=:release, met=P): forward `a->b` has no metabolite, reverse `b->a` has P as concentration
- Isomerization: neither direction has metabolite

```julia
R_conc = Dict{Tuple{Int,Int}, Set{Vector{Pair{Symbol,Int}}}}()
for (idx, (a, b)) in enumerate(edges)
    eq_steps[idx] && continue  # only SS steps
    info = adj[minmax(a, b)]
    g1, g2 = form_to_group[a], form_to_group[b]
    g1 == g2 && continue  # intra-group SS step

    # Forward: a->b
    fwd_met = (info.type == :binding) ? info.metabolite : nothing
    fwd_mono = fwd_met !== nothing ?
        _add_met(alpha_conc[a], fwd_met) : copy(alpha_conc[a])
    push!(get!(R_conc, (g1, g2), Set{Vector{Pair{Symbol,Int}}}()), fwd_mono)

    # Reverse: b->a
    rev_met = (info.type == :release) ? info.metabolite : nothing
    rev_mono = rev_met !== nothing ?
        _add_met(alpha_conc[b], rev_met) : copy(alpha_conc[b])
    push!(get!(R_conc, (g2, g1), Set{Vector{Pair{Symbol,Int}}}()), rev_mono)
end
```

The implementing agent should cross-reference with `_raw_symbolic_rate_polys` (rate_eq_derivation.jl lines 382-408) which uses `_split_reaction_side` to get metabolites, to ensure the fingerprint matches the actual rate equation.

#### 2.1e — Compute cofactor concentration monomials D_conc[g] and combine

For each root `g`, enumerate all spanning arborescences of the G-node directed group graph rooted at g. An arborescence rooted at g is a selection of exactly one incoming edge `(parent -> child)` for each non-root node, such that following parents from any node leads to g.

For G <= 7, use recursive enumeration. The Matrix-Tree theorem guarantees all cofactor terms are non-negative, so no cancellation occurs — the monomial set is the union over all spanning trees.

**Simplifications for common cases:**
- G=2: `D_conc[1] = R_conc[2,1]`, `D_conc[2] = R_conc[1,2]`. No tree enumeration needed.
- G=3: 3 possible arborescences per root. Enumerate directly.

**Final fingerprint:**
```julia
fingerprint = Set{Vector{Pair{Symbol,Int}}}()
for g in 1:G
    D_g = _spanning_arborescences(G, R_conc, g)
    for s in sigma_conc[g], d in D_g
        push!(fingerprint, _mono_mul_conc(s, d))
    end
end
return fingerprint
```

### 2.2 Add test-local `_concentration_fingerprint_test`

In the test file, add an independent reimplementation for verification. It should use only public `MechanismSpec` fields and `EnzymeFormSpec` data.

### 2.3 Add unit tests for the fingerprint function

In `test/test_mechanism_enum_of_enz_reaction.jl`, add a new `@testset "Concentration fingerprint"`:

1. **Uni-uni 3-step:** Verify all three G=2 partitions (`{E,ES}|{EP}`, `{ES,EP}|{E}`, `{E,EP}|{ES}`) produce the same fingerprint `{[], [S=>1], [P=>1]}`. This confirms the analysis from our brainstorming.

2. **Ordered bi-uni 4-step:** Verify the G=2 partitions with different sigma content produce DIFFERENT fingerprints.

3. **Symmetry check:** For partitions Y `{EA,EAB}|{EP}|{E}` and Z `{E,EP}|{EA}|{EAB}` from the 4-step mechanism, verify they produce the SAME fingerprint `{[], [A=>1], [B=>1], [P=>1], [A=>1,B=>1], [B=>1,P=>1]}`.

These tests validate the fingerprint computation against hand-computed results.

---

## Step 3 — Integrate Phase 2 into enumeration

**Depends on Steps 1 and 2 being complete.**

### 3.1 Replace partition dedup with fingerprint dedup

In `_ress_variants` and `_count_ress_variants`, change the `seen` set from `Set{Vector{Vector{Int}}}` (partition) to `Set{Set{Vector{Pair{Symbol,Int}}}}` (fingerprint).

Option: keep the partition check as a fast pre-filter before the more expensive fingerprint computation:

```julia
# Fast pre-filter: skip if partition already seen
partition = _compute_re_partition(edges, eq_steps)
partition in seen_partitions && return ()
push!(seen_partitions, partition)

# Slower full check: skip if fingerprint matches a different partition
fp = _concentration_fingerprint(edges, eq_steps, forms, adj)
fp in seen_fingerprints && return ()
push!(seen_fingerprints, fp)
```

Profile to see if the two-level check is worth the complexity. If fingerprint computation is fast enough (likely, given G <= 7), just use fingerprint alone.

### 3.2 Update test helper

Apply matching changes to `_compute_expected_n_total` in the test file.

### 3.3 Compute new expected counts

Run the code again to get Phase-2-reduced counts. Update `expected_n_total` in all test specs.

### 3.4 Add validation test: same-fingerprint mechanisms produce same rate equation

In `test/test_mechanism_enum_of_enz_reaction.jl`, add:

```julia
@testset "Fingerprint-equivalent rate equations" begin
    # For uni-uni, compile all 3 G=2 partitions (before dedup) and
    # verify they produce the same rate equation polynomial form.
    rxn = @enzyme_reaction begin
        substrates:S[C]
        products:P[C]
    end
    forms = EnzymeRates.enumerate_enzyme_forms(rxn)
    adj = EnzymeRates._build_adjacency(forms)
    # Manually construct 3 MechanismSpecs with different SS assignments
    # ... compile each, call rate_equation_string, extract monomial sets,
    # verify they match.
end
```

This test proves the fingerprint correctly identifies equivalent models by checking against the actual symbolic rate equation.

### 3.5 Run full test suite

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

## Agent Assignment

| Agent | Steps | Dependencies | Notes |
|-------|-------|-------------|-------|
| **A** | 1.1-1.6 | None | Phase 1: partition dedup in prod + test code, new expected counts, all tests pass |
| **B** | 2.1-2.3 | None (parallel with A) | Phase 2 foundation: fingerprint function + unit tests. Does NOT modify existing enumeration code. |
| **C** | 3.1-3.5 | A and B complete | Integration: replace partition with fingerprint dedup, update counts, validation test |

---

## What NOT to change

- `enumerate_enzyme_forms` — unchanged (no dedup at form level)
- `_catalytic_topologies` — unchanged (catalytic cycle enumeration is orthogonal)
- `_expand_inhibitors` — unchanged (dead-end expansion is orthogonal)
- `_expand_oligomeric_variants` — unchanged (OEM expansion wraps RE/SS results)
- `compile_mechanism`, `EnzymeMechanism(spec)` — unchanged
- `expected_n_forms`, `expected_n_catalytic`, `expected_n_cat_de` in test specs — unchanged (dedup only affects RE/SS stage)
