# Mechanism Enumeration Redesign — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reimplement mechanism enumeration as a staged pipeline with richer regulator roles, incremental parameter counting, and independently testable expansion stages.

**Architecture:** Each expansion stage is a pure function `Vector{MechanismSpec} -> Vector{MechanismSpec}`. Stages compose into a pipeline orchestrated by `enumerate_mechanisms`. Types split into `MechanismSpec` (stages 1-7) and `OligomericMechanismSpec` (stages 8-10). `EnzymeReaction` gains regulator role info in type parameters.

**Tech Stack:** Julia, existing EnzymeRates.jl infrastructure (types.jl, dsl.jl, rate_eq_derivation.jl)

**Design doc:** `docs/plans/2026-03-07-mechanism-enumeration-redesign.md`

---

## Task 0: Rename old files

Preserve the old implementation while developing the new one.

**Files:**
- Rename: `src/mechanism_enumeration.jl` -> `src/old_mechanism_enumeration.jl`
- Rename: `test/test_mechanism_enum_of_enz_reaction.jl` -> `test/old_test_mechanism_enum_of_enz_reaction.jl`
- Rename: `test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl` -> `test/old_reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl`
- Modify: `src/EnzymeRates.jl:35` — change `include("mechanism_enumeration.jl")` to `include("old_mechanism_enumeration.jl")`
- Modify: `test/runtests.jl:7` — change `include("reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl")` to `include("old_reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl")`
- Modify: `test/runtests.jl:17` — change `include("test_mechanism_enum_of_enz_reaction.jl")` to `include("old_test_mechanism_enum_of_enz_reaction.jl")`

**Step 1: Rename files**

```bash
cd /home/denis.linux/.julia/dev/EnzymeRates
git mv src/mechanism_enumeration.jl src/old_mechanism_enumeration.jl
git mv test/test_mechanism_enum_of_enz_reaction.jl test/old_test_mechanism_enum_of_enz_reaction.jl
git mv test/reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl test/old_reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl
```

**Step 2: Update includes in `src/EnzymeRates.jl`**

Change line 35 from:
```julia
include("mechanism_enumeration.jl")
```
to:
```julia
include("old_mechanism_enumeration.jl")
```

**Step 3: Update includes in `test/runtests.jl`**

Change line 7 from:
```julia
include("reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl")
```
to:
```julia
include("old_reaction_definitions_for_test_mechanism_enum_of_enz_reaction.jl")
```

Change line 17 from:
```julia
include("test_mechanism_enum_of_enz_reaction.jl")
```
to:
```julia
include("old_test_mechanism_enum_of_enz_reaction.jl")
```

**Step 4: Verify old tests still pass**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass (only file renames, no code changes).

**Step 5: Commit**

```bash
git add -A && git commit -m "Rename old mechanism enumeration files with old_ prefix"
```

---

## Task 1: Add `RegulatorRole` types and update `EnzymeReaction`

Update the type system to support regulator roles.

**Files:**
- Modify: `src/types.jl` — add `RegulatorRole` hierarchy, change `EnzymeReaction` constructor and `regulators()` accessor
- Modify: `src/dsl.jl` — update `@enzyme_reaction` to parse new DSL labels
- Modify: `test/test_dsl.jl` — add tests for new DSL syntax
- Modify: `test/test_types.jl` — add tests for `RegulatorRole` types

**Step 1: Add `RegulatorRole` types to `src/types.jl`**

Add after the `EnzymeReaction` struct definition (after line 14):

```julia
"""Regulator role in mechanism enumeration."""
abstract type RegulatorRole end

"""Allosteric regulator: OEM expansion + essential activator +
general modifier special cases."""
struct Allosteric <: RegulatorRole end

"""Dead-end inhibitor: creates dead-end complexes only."""
struct DeadEnd <: RegulatorRole end

"""Unconstrained: try all roles (allosteric + dead-end)."""
struct UnconstrainedRegulator <: RegulatorRole end
```

**Step 2: Update `EnzymeReaction` constructor in `src/types.jl`**

The `Regulators` type parameter changes from a tuple of Symbols to a tuple of `(Symbol, Symbol)` pairs where the second element is the role:

```julia
# Old: EnzymeReaction{subs, prods, (:R1, :R2)}()
# New: EnzymeReaction{subs, prods, ((:R1, :unknown), (:R2, :dead_end))}()
```

Update the `EnzymeReaction` constructor (lines 16-31) to accept regulators as either:
- A tuple of `Symbol`s (backward compatibility) — converted to `(name, :unknown)` pairs
- A tuple of `(Symbol, Symbol)` pairs (new format)

```julia
function EnzymeReaction(subs::Tuple, prods::Tuple, regs::Tuple=())
    isempty(subs) && error("Substrates must not be empty")
    isempty(prods) && error("Products must not be empty")
    subs_names = [s[1] for s in subs]
    prods_names = [s[1] for s in prods]
    length(subs_names) != length(Set(subs_names)) &&
        error("Duplicate substrate names")
    length(prods_names) != length(Set(prods_names)) &&
        error("Duplicate product names")
    # Normalize regulators to (name, role) pairs
    normalized_regs = if !isempty(regs) && regs[1] isa Symbol
        Tuple((r, :unknown) for r in regs)
    else
        regs
    end
    for r in normalized_regs
        r isa Tuple{Symbol, Symbol} ||
            error("Regulators must be (Symbol, Symbol) pairs, got $r")
    end
    reg_names = [r[1] for r in normalized_regs]
    length(reg_names) != length(Set(reg_names)) &&
        error("Duplicate regulator names")
    subs = _sort_species(subs)
    prods = _sort_species(prods)
    sorted_regs = Tuple(sort(collect(normalized_regs); by=first))
    EnzymeReaction{subs, prods, sorted_regs}()
end
```

**Step 3: Update `regulators()` accessor**

The `regulators` function on `EnzymeReaction` (line 409 of types.jl) returns just the names (backward compatible):

```julia
regulators(::EnzymeReaction{S,P,R}) where {S,P,R} =
    Tuple(r[1] for r in R)
```

Add a new accessor for roles:

```julia
"""Return regulator (name, role) pairs."""
regulator_roles(::EnzymeReaction{S,P,R}) where {S,P,R} = R
```

**Step 4: Update `@enzyme_reaction` DSL in `src/dsl.jl`**

Update `_parse_labeled_block` call in the `@enzyme_reaction` macro (line 83) to accept new labels:

```julia
macro enzyme_reaction(block)
    parsed = _parse_labeled_block(block,
        Set([:substrates, :products, :regulators,
             :dead_end_inhibitors, :allosteric_regulators]))
    haskey(parsed, :substrates) || error("substrates not specified")
    haskey(parsed, :products) || error("products not specified")
    # Build regulator list with roles
    regs = Expr(:tuple)
    for (label, role_sym) in [
        (:regulators, :unknown),
        (:dead_end_inhibitors, :dead_end),
        (:allosteric_regulators, :allosteric),
    ]
        if haskey(parsed, label)
            syms = _regulator_tuple_to_symbols(parsed[label])
            for s in syms.args
                push!(regs.args,
                    Expr(:tuple, s, QuoteNode(role_sym)))
            end
        end
    end
    return esc(:(EnzymeReaction($subs, $prods, $regs)))
end
```

Note: The existing `_regulator_tuple_to_symbols` function already extracts Symbol names from parsed species tuples, so it can be reused.

**Step 5: Write tests for new DSL in `test/test_dsl.jl`**

Add to the `@enzyme_reaction` testset:

```julia
# New regulator role syntax
spec_roles = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    dead_end_inhibitors: I
    allosteric_regulators: A
    regulators: R
end
@test spec_roles isa EnzymeReaction
@test Set(EnzymeRates.regulators(spec_roles)) == Set([:I, :A, :R])
roles = EnzymeRates.regulator_roles(spec_roles)
@test length(roles) == 3
role_dict = Dict(r[1] => r[2] for r in roles)
@test role_dict[:I] == :dead_end
@test role_dict[:A] == :allosteric
@test role_dict[:R] == :unknown

# Backward compatibility: plain regulators
spec_plain = @enzyme_reaction begin
    substrates: S[C]
    products: P[C]
    regulators: R1, R2
end
roles_plain = EnzymeRates.regulator_roles(spec_plain)
@test all(r[2] == :unknown for r in roles_plain)
@test Set(r[1] for r in roles_plain) == Set([:R1, :R2])
```

**Step 6: Run tests**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all tests pass. Existing tests that use `regulators()` should still work since it returns just names.

**Step 7: Commit**

```bash
git commit -m "Add RegulatorRole types and update EnzymeReaction for regulator roles"
```

---

## Task 2: New `MechanismSpec` and `OligomericMechanismSpec` types

Create the new type hierarchy for mechanism specs.

**Files:**
- Create: `src/mechanism_enumeration.jl` — new file with types only (stages added in later tasks)
- Modify: `src/EnzymeRates.jl` — add `include("mechanism_enumeration.jl")` after old include

**Step 1: Create `src/mechanism_enumeration.jl` with types**

```julia
# ─── New Mechanism Enumeration Types ───────────────────────────

abstract type AbstractMechanismSpec end

struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    edges::Vector{Tuple{Int,Int}}
    n_catalytic_edges::Int
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{ParamConstraint}
    param_count::Int
end

# Convenience constructors
function MechanismSpec(reaction, edges, equilibrium_steps;
                       param_count::Int=0)
    MechanismSpec(reaction, edges, length(edges),
        equilibrium_steps, ParamConstraint[], param_count)
end

struct OligomericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    tr_equivalence::Vector{Bool}
end
```

**Step 2: Add include to `src/EnzymeRates.jl`**

After the `include("old_mechanism_enumeration.jl")` line, add:

```julia
include("mechanism_enumeration.jl")
```

**Step 3: Run tests to verify no conflicts**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `MechanismSpec` is now redefined (conflicts with old file). We need to handle this.

The old file defines `MechanismSpec` too. Since both files are included, we'll get a redefinition error. Two options:
- a) Wrap old types in a module: `module OldEnumeration ... end`
- b) Name the new type differently temporarily

**Better approach:** Don't include the new file yet. Only include it when we're ready to replace the old one. For now, put the new types in the test file for testing, and move them into `src/mechanism_enumeration.jl` when we swap old for new.

**Revised Step 2:** Don't modify `src/EnzymeRates.jl` yet. Create the file but don't include it. We'll develop and test using direct `include` from the test file, then swap at the end.

**Step 3: Verify tests still pass (no changes to includes)**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

**Step 4: Commit**

```bash
git commit -m "Add new MechanismSpec and OligomericMechanismSpec type definitions"
```

---

## Task 3: Shared infrastructure — form enumeration, adjacency, helpers

Port the reusable functions from old implementation to new file. These functions are shared between old and new: `enumerate_enzyme_forms`, `_build_adjacency`, `_classify_edge`, `_is_binding_direction`, `_compute_re_partition`, concentration fingerprint functions.

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add shared helper functions

**Step 1: Add helper functions**

These are copied from `src/old_mechanism_enumeration.jl` with minimal changes. The key functions to port:

- `_classify_edge` (lines 141-196) — unchanged
- `_build_adjacency` (lines 205-214) — unchanged
- `_is_binding_direction` (lines 223-230) — unchanged
- `_compute_re_partition` (lines 497-518) — unchanged
- `_concentration_fingerprint` (lines 572-652) — unchanged
- `_constraint_descriptor` (lines 661-671) — unchanged
- `_dead_end_catalytic_map` (lines 695-707) — unchanged
- `_propagate_de_eq_steps!` (lines 716-721) — unchanged
- `_find_equivalent_groups` (lines 733-757) — update for new equivalence rules
- `_build_constraints` (lines 762-772) — unchanged
- `_set_partitions` (lines 998-1013) — unchanged
- `_partition_mult_count` (lines 1022-1032) — unchanged

Note: `enumerate_enzyme_forms` (lines 241-307), `SiteDefinition`, `EnzymeFormSpec` are already in the old file and shared. Since both old and new files are included, these don't need duplicating. But since we're building the new file to eventually replace the old one, we should include all needed functions.

**Important:** Since both old and new files will eventually coexist temporarily, and many helper functions are identical, we should either:
- (a) Keep helpers in old file and reference them (fragile)
- (b) Move shared helpers to a separate file
- (c) Duplicate in new file (will replace old later)

Go with (c) — duplicate in new file. When we swap, old file gets deleted and new file is self-contained.

Add all the above functions to `src/mechanism_enumeration.jl`. Also add `SiteDefinition`, `EnzymeFormSpec`, `enumerate_enzyme_forms` since the new file needs to be self-contained when it replaces the old.

**BUT:** Since we can't have two definitions of `SiteDefinition` etc. while both files are included, we need a different approach.

**Revised approach:** Develop the new `mechanism_enumeration.jl` as a complete replacement. Don't include it until ready. Test it via a separate test script that loads it directly. When ready, swap the include.

**Step 1: Copy shared types and functions into new file**

Copy from old file into new `src/mechanism_enumeration.jl`:
- `SiteDefinition` struct
- `EnzymeFormSpec` struct
- `ParamConstraint` const (already in old file, but this is used by new `MechanismSpec`)
- `enumerate_enzyme_forms`
- All `_classify_edge`, `_build_adjacency`, `_is_binding_direction` helpers
- All `_compute_re_partition`, fingerprint, constraint helpers
- `_dead_end_catalytic_map`, `_propagate_de_eq_steps!`
- `_find_equivalent_groups` (with updated site-type-aware grouping)
- `_build_constraints`
- `_set_partitions`, `_partition_mult_count`

**Step 2: Update `_find_equivalent_groups` for new equivalence rules**

The key change from the design: equivalence groups must only contain edges binding the same metabolite AND same site type. Update to track site role alongside metabolite:

```julia
function _find_equivalent_groups(edges, adj, site_defs, forms,
    n_catalytic_edges, de_cat_map=nothing)
    # Group by (metabolite, site_role) instead of just metabolite
    groups = Dict{Tuple{Symbol, Symbol}, Vector{Int}}()
    product_metabolites = Set(
        sd.metabolite for sd in site_defs if sd.role == :prod)
    for i in 1:n_catalytic_edges
        (a, b) = edges[i]
        met = get(adj, minmax(a, b), missing)
        !ismissing(met) && met !== nothing &&
            met ∉ product_metabolites || continue
        # Determine site role from the differing site
        site_role = _edge_site_role(edges[i], site_defs, forms)
        push!(get!(groups, (met, site_role), Int[]), i)
    end
    if de_cat_map !== nothing
        for (di, cat_idx) in enumerate(de_cat_map)
            cat_idx === nothing && continue
            edge_idx = n_catalytic_edges + di
            met = adj[minmax(edges[edge_idx]...)]
            met === nothing && continue
            met in product_metabolites && continue
            site_role = _edge_site_role(
                edges[edge_idx], site_defs, forms)
            push!(get!(groups, (met, site_role), Int[]),
                edge_idx)
        end
    end
    sort!([sort(v) for v in values(groups) if length(v) >= 2];
        by=first)
end

"""Return the site role (:sub, :prod, :reg) for a binding edge."""
function _edge_site_role(edge, site_defs, forms)
    (a, b) = edge
    for k in eachindex(site_defs)
        forms[a].occupancy[k] != forms[b].occupancy[k] &&
            return site_defs[k].role
    end
    :iso  # isomerization (shouldn't reach here for binding edges)
end
```

**Step 3: Commit**

```bash
git commit -m "Add shared infrastructure to new mechanism_enumeration.jl"
```

---

## Task 4: Stage 1 — Catalytic Topologies

Port catalytic topology enumeration with the new `MechanismSpec` format and `param_count` initialization.

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_catalytic_topologies`

**Step 1: Implement `_catalytic_topologies`**

Port from old file (lines 317-410) with these changes:
- Signature: `_catalytic_topologies(reaction) -> Vector{MechanismSpec}`
- Compute `site_defs`, `forms`, `adj` internally from `reaction`
- Initialize `param_count` for each topology:
  - Count edges in the topology
  - Compute thermodynamic constraints (1 Haldane + n_independent_cycles - 1 Wegscheider)
  - `param_count = n_edges - n_thermo_constraints + 2` (E_total + Keq)
  - But wait: the initial topology has first isomerization SS, all others RE. So SS edges contribute 2 params, RE edges contribute 1. More precisely: `param_count = n_RE_edges * 1 + n_SS_edges * 2 - n_thermo_constraints + 2`
  - For the initial state (1 SS isomerization, rest RE): `param_count = (n_edges - 1) + 2 - n_thermo + 2 = n_edges + 3 - n_thermo`
- Set `equilibrium_steps`: first isomerization step is `false` (SS), all others `true` (RE)

```julia
function _catalytic_topologies(
    @nospecialize(reaction::EnzymeReaction)
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)
    max_forms = length(forms)
    # ... (DFS code from old implementation, adapted) ...

    # For each topology, set default RE/SS and compute param_count
    result = MechanismSpec[]
    for form_set in combined
        edges = [(a, b) for ((a, b), _) in adj
                 if a ∈ form_set && b ∈ form_set]
        # Find first isomerization
        iso_idx = 0
        for (i, (a, b)) in enumerate(edges)
            adj[minmax(a, b)] === nothing && (iso_idx = i; break)
        end
        iso_idx == 0 && (iso_idx = 1)
        eq_steps = fill(true, length(edges))
        eq_steps[iso_idx] = false

        n_edges = length(edges)
        n_forms = length(form_set)
        # Haldane: 1 per independent cycle
        # For simple cycle: n_independent_cycles = n_edges - n_forms + 1
        n_independent_cycles = n_edges - n_forms + 1
        n_thermo = n_independent_cycles  # 1 Haldane per cycle
        # RE edges: n_edges - 1, SS edges: 1
        n_re = n_edges - 1
        n_ss = 1
        param_count = n_re + 2 * n_ss - n_thermo + 2  # +2 for E_total, Keq

        push!(result, MechanismSpec(reaction, edges, length(edges),
            eq_steps, ParamConstraint[], param_count))
    end
    result
end
```

**Step 2: Commit**

```bash
git commit -m "Add Stage 1: catalytic topology enumeration"
```

---

## Task 5: Stage 2 — RE/SS Assignment

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_expand_ress_variants`

**Step 1: Implement `_expand_ress_variants`**

```julia
function _expand_ress_variants(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    max_re_groups::Int=7,
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)
    result = MechanismSpec[]
    # Always include originals (they already have valid RE/SS from stage 1)
    append!(result, specs)

    for spec in specs
        edges = spec.edges
        n = length(edges)
        n == 0 && continue

        # Find all edge indices (any can be SS)
        # Constraint: at least 1 SS step, G >= 1, G <= max_re_groups
        for ss_mask in 1:(1 << n) - 1  # skip 0 (no SS steps)
            eq_steps = [((ss_mask >> (i-1)) & 1) == 0
                        for i in 1:n]
            all(eq_steps) && continue  # no SS steps

            # Skip if identical to the original spec's eq_steps
            eq_steps == spec.equilibrium_steps && continue

            partition = _compute_re_partition(edges, eq_steps)
            G = length(partition)
            G > max_re_groups && continue

            # Compute param_count delta from base
            n_ss_new = count(!, eq_steps)
            n_ss_orig = count(!, spec.equilibrium_steps)
            delta = n_ss_new - n_ss_orig  # each RE->SS adds 1 param

            push!(result, MechanismSpec(reaction, edges,
                spec.n_catalytic_edges, eq_steps,
                ParamConstraint[], spec.param_count + delta))
        end
    end
    result
end
```

**Step 2: Commit**

```bash
git commit -m "Add Stage 2: RE/SS assignment expansion"
```

---

## Task 6: Stage 3 — General Modifier Expansion

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_expand_general_modifiers`

**Step 1: Implement `_expand_general_modifiers`**

For each regulator going through the general modifier path: duplicate catalytic cycle with R bound (parallel paths), add R-binding edges (always RE). Propagate RE/SS from catalytic counterparts.

```julia
function _expand_general_modifiers(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)
    form_lookup = Dict(
        ntuple(k -> f.occupancy[k], length(site_defs)) => i
        for (i, f) in enumerate(forms))

    # Get allosteric regulators (those that could be general modifiers)
    allosteric_regs = [r[1] for r in regulator_roles(reaction)
                       if r[2] in (:allosteric, :unknown)]
    isempty(allosteric_regs) && return specs

    result = copy(specs)  # preserve originals

    for spec in specs
        # For each subset of allosteric regulators to add as general modifiers
        # (only add all together, not mixed with OEM — per design constraint)
        # Create modifier variant with all allosteric regs as general modifiers
        _add_general_modifier_variant!(
            result, spec, allosteric_regs,
            site_defs, forms, adj, form_lookup)
    end
    result
end
```

The helper `_add_general_modifier_variant!` needs to:
1. Find forms where regulators are bound (ER, ESR, etc.)
2. Add edges between those forms (duplicating the catalytic cycle structure)
3. Add R-binding edges (E<->ER, ES<->ESR) as always RE
4. Propagate RE/SS from catalytic counterparts
5. Compute param_count delta

**Step 2: Commit**

```bash
git commit -m "Add Stage 3: general modifier expansion"
```

---

## Task 7: Stage 4 — Essential Activator Expansion

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_expand_essential_activators`

**Step 1: Implement `_expand_essential_activators`**

For each regulator going through essential activator path: replace catalytic cycle with R-bound version. Only ER+S<->ESR<->ER+P path exists. E<->ER binding edge added.

The key difference from general modifier: the original catalytic edges (E+S<->ES) are removed, replaced with R-bound versions (ER+S<->ESR).

```julia
function _expand_essential_activators(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)
    form_lookup = Dict(
        ntuple(k -> f.occupancy[k], length(site_defs)) => i
        for (i, f) in enumerate(forms))

    allosteric_regs = [r[1] for r in regulator_roles(reaction)
                       if r[2] in (:allosteric, :unknown)]
    isempty(allosteric_regs) && return specs

    result = copy(specs)  # preserve originals

    for spec in specs
        _add_essential_activator_variant!(
            result, spec, allosteric_regs,
            site_defs, forms, adj, form_lookup)
    end
    result
end
```

The helper `_add_essential_activator_variant!` needs to:
1. Find R-bound counterparts for each catalytic form
2. Build new edge list: R-bound catalytic edges + E<->ER binding
3. No E+S<->ES edges
4. Propagate RE/SS from original catalytic edges
5. Compute param_count

**Step 2: Commit**

```bash
git commit -m "Add Stage 4: essential activator expansion"
```

---

## Task 8: Stage 5 — Dead-End Inhibitor Expansion

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_expand_dead_end_inhibitors`

**Step 1: Implement `_expand_dead_end_inhibitors`**

Port from old `_expand_inhibitors` (lines 431-487) with changes:
- Signature: `_expand_dead_end_inhibitors(specs, reaction) -> Vector{MechanismSpec}`
- Compute forms/adj internally
- Only consider regulators with `DeadEnd()` or `UnconstrainedRegulator()` role
- Track `param_count` increment: +1 per dead-end RE edge, +2 per dead-end SS edge
- Preserve originals in output
- Dead-end edges inherit RE/SS from catalytic counterparts via `_propagate_de_eq_steps!`

```julia
function _expand_dead_end_inhibitors(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)
    form_lookup = Dict(
        ntuple(k -> f.occupancy[k], length(site_defs)) => i
        for (i, f) in enumerate(forms))
    max_forms = length(forms)

    dead_end_regs = [r[1] for r in regulator_roles(reaction)
                     if r[2] in (:dead_end, :unknown)]
    isempty(dead_end_regs) && return specs

    result = copy(specs)  # preserve originals

    for spec in specs
        # ... port logic from old _expand_inhibitors ...
        # For each dead-end mask combination:
        #   1. Find dead-end forms
        #   2. Add edges
        #   3. Propagate RE/SS
        #   4. Compute param_count delta
        #   5. Push new MechanismSpec
    end
    result
end
```

**Step 2: Commit**

```bash
git commit -m "Add Stage 5: dead-end inhibitor expansion"
```

---

## Task 9: Stage 6 — Equivalence Constraints

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_expand_equivalence_constraints`

**Step 1: Implement `_expand_equivalence_constraints`**

Port from old `_ress_variants` constraint logic with the updated equivalence grouping rules:
- Same metabolite AND same site type can group
- Substrate-site and regulator-site for same metabolite CANNOT group
- Two different regulatory sites for same metabolite CANNOT group

```julia
function _expand_equivalence_constraints(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)

    result = MechanismSpec[]

    for spec in specs
        edges = spec.edges
        n_cat = spec.n_catalytic_edges
        de_cat_map = _dead_end_catalytic_map(
            edges, n_cat, site_defs, forms)
        equiv_groups = _find_equivalent_groups(
            edges, adj, site_defs, forms, n_cat, de_cat_map)

        # Filter to groups where all edges share same RE/SS status
        valid_groups = [g for g in equiv_groups
            if all(spec.equilibrium_steps[s] ==
                   spec.equilibrium_steps[g[1]] for s in g)]

        # Always include the unconstrained original
        push!(result, spec)

        # Enumerate constraint masks (skip 0 = no constraints)
        for mask in 1:(1 << length(valid_groups)) - 1
            constraints = _build_constraints(
                valid_groups, spec.equilibrium_steps, mask, edges)
            # Compute param_count delta
            delta = 0
            for (gi, g) in enumerate(valid_groups)
                (mask >> (gi-1)) & 1 == 1 || continue
                n_constrained = length(g) - 1
                if spec.equilibrium_steps[g[1]]  # RE
                    delta -= n_constrained
                else  # SS
                    delta -= 2 * n_constrained
                end
            end
            push!(result, MechanismSpec(spec.reaction, edges,
                n_cat, spec.equilibrium_steps,
                constraints, spec.param_count + delta))
        end
    end
    result
end
```

**Step 2: Commit**

```bash
git commit -m "Add Stage 6: equivalence constraint expansion"
```

---

## Task 10: Stage 7 — Deduplication

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_deduplicate`

**Step 1: Implement `_deduplicate`**

Port fingerprint-based dedup from old `_ress_variants`. Uses (concentration fingerprint, constraint descriptor) as dedup key. Keeps mechanism with fewest parameters.

```julia
function _deduplicate(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    isempty(specs) && return specs
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)

    best = Dict{_DedupKey, MechanismSpec}()
    for spec in specs
        edges = spec.edges
        eq_steps = spec.equilibrium_steps
        partition = _compute_re_partition(edges, eq_steps)
        fp = _concentration_fingerprint(
            edges, eq_steps, site_defs, forms, adj, partition)

        n_cat = spec.n_catalytic_edges
        de_cat_map = _dead_end_catalytic_map(
            edges, n_cat, site_defs, forms)
        equiv_groups = _find_equivalent_groups(
            edges, adj, site_defs, forms, n_cat, de_cat_map)
        valid_groups = [g for g in equiv_groups
            if all(eq_steps[s] == eq_steps[g[1]] for s in g)]

        # Compute constraint mask from param_constraints
        constraint_mask = _constraints_to_mask(
            spec.param_constraints, valid_groups, eq_steps, edges)
        desc = _constraint_descriptor(
            edges, adj, eq_steps, valid_groups, constraint_mask)

        key = (fp, desc)
        if !haskey(best, key) ||
                spec.param_count < best[key].param_count
            best[key] = spec
        end
    end
    collect(values(best))
end
```

Note: Need to implement `_constraints_to_mask` helper that reverses the constraint → mask mapping.

**Step 2: Commit**

```bash
git commit -m "Add Stage 7: deduplication by concentration fingerprint"
```

---

## Task 11: Stage 8 — Allosteric (OEM) Expansion

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_expand_allosteric`

**Step 1: Implement `_expand_allosteric`**

Port from old `_expand_oligomeric_variants` and `_set_partitions`. Only applied to mechanisms that did NOT go through stages 3/4 (general modifier / essential activator).

```julia
function _expand_allosteric(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    catalytic_n::Int=2,
)
    allosteric_regs = [r[1] for r in regulator_roles(reaction)
                       if r[2] in (:allosteric, :unknown)]
    isempty(allosteric_regs) && return OligomericMechanismSpec[]

    result = OligomericMechanismSpec[]
    partitions = _set_partitions(allosteric_regs)

    for spec in specs
        # Skip specs that were expanded via general modifier / essential activator
        # (detected by checking if spec has regulator-bound forms in catalytic edges)
        _was_modifier_expanded(spec, reaction) && continue

        for partition in partitions
            n_groups = length(partition)
            for combo in Iterators.product(
                    ntuple(_ -> 1:catalytic_n, n_groups)...)
                push!(result, OligomericMechanismSpec(
                    spec, catalytic_n, partition,
                    collect(combo), Bool[]))
            end
        end
    end
    result
end
```

**Step 2: Commit**

```bash
git commit -m "Add Stage 8: allosteric OEM expansion"
```

---

## Task 12: Stage 9 — T/R Equivalence

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_expand_tr_equivalence`

**Step 1: Implement `_expand_tr_equivalence`**

For each OEM spec, enumerate which parameter groups are T=R equivalent vs independent. Each parameter group corresponds to a binding step or catalytic step that has mirror parameters in T and R states.

```julia
function _expand_tr_equivalence(
    specs::Vector{OligomericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = OligomericMechanismSpec[]

    for spec in specs
        base = spec.base
        # Identify parameter groups that can be T=R constrained
        # Each edge in the catalytic mechanism = one parameter group
        n_groups = length(base.edges)

        # Enumerate all T/R equivalence masks
        for mask in 0:(1 << n_groups) - 1
            tr_equiv = [((mask >> (i-1)) & 1) == 1
                        for i in 1:n_groups]
            push!(result, OligomericMechanismSpec(
                base, spec.catalytic_n,
                spec.allosteric_reg_sites,
                spec.allosteric_multiplicities,
                tr_equiv))
        end
    end
    result
end
```

**Step 2: Commit**

```bash
git commit -m "Add Stage 9: T/R equivalence expansion"
```

---

## Task 13: Stage 10 — Post-OEM Deduplication

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `_deduplicate_oem`

**Step 1: Implement `_deduplicate_oem`**

Check for T<->R mirror duplicates. Two OEM specs are mirrors if swapping T/R labels (and L -> 1/L) produces the same rate equation.

```julia
function _deduplicate_oem(
    specs::Vector{OligomericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    # For now, implement a simple canonical form:
    # Sort T/R equivalence vector and compare
    # Full mirror detection may be added later if needed
    seen = Dict{Any, OligomericMechanismSpec}()
    for spec in specs
        key = _oem_canonical_key(spec)
        if !haskey(seen, key)
            seen[key] = spec
        end
    end
    collect(values(seen))
end
```

**Step 2: Commit**

```bash
git commit -m "Add Stage 10: post-OEM deduplication"
```

---

## Task 14: Pipeline orchestration — `enumerate_mechanisms`

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `enumerate_mechanisms` function

**Step 1: Implement pipeline routing**

```julia
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    max_re_groups::Int=7,
    catalytic_n::Int=2,
)
    # Stage 1: Catalytic topologies
    catalytic = _catalytic_topologies(reaction)

    # Stage 2: RE/SS assignment
    with_ress = _expand_ress_variants(
        catalytic, reaction; max_re_groups)

    # Route regulators into paths
    regs = regulator_roles(reaction)
    has_allosteric = any(r -> r[2] in (:allosteric, :unknown), regs)
    has_dead_end = any(r -> r[2] in (:dead_end, :unknown), regs)

    # Path A: general modifier + essential activator
    path_a = with_ress
    if has_allosteric
        path_a = _expand_general_modifiers(path_a, reaction)
        path_a = _expand_essential_activators(path_a, reaction)
    end

    # Stage 5: Dead-end inhibitors (both paths)
    if has_dead_end
        path_a = _expand_dead_end_inhibitors(path_a, reaction)
    end

    # Stage 6: Equivalence constraints
    path_a = _expand_equivalence_constraints(path_a, reaction)

    # Stage 7: Dedup
    em_specs = _deduplicate(path_a, reaction)

    # Path B: Allosteric OEM expansion (only for non-modifier specs)
    oem_specs = OligomericMechanismSpec[]
    if has_allosteric && catalytic_n > 0
        oem_specs = _expand_allosteric(
            em_specs, reaction; catalytic_n)
        oem_specs = _expand_tr_equivalence(oem_specs, reaction)
        oem_specs = _deduplicate_oem(oem_specs, reaction)
    end

    # Merge results
    AbstractMechanismSpec[em_specs; oem_specs]
end
```

**Step 2: Commit**

```bash
git commit -m "Add enumerate_mechanisms pipeline orchestration"
```

---

## Task 15: Update `compile_mechanism` for new types

**Files:**
- Modify: `src/mechanism_enumeration.jl` — add `compile_mechanism` methods for new types

**Step 1: Add `compile_mechanism(::MechanismSpec)` and `compile_mechanism(::OligomericMechanismSpec)`**

Port from old file (lines 927-989). The `MechanismSpec` version is largely identical. The `OligomericMechanismSpec` version uses the `base` field.

```julia
function compile_mechanism(spec::MechanismSpec)
    # ... same as old EnzymeMechanism(spec::MechanismSpec) ...
end

function compile_mechanism(spec::OligomericMechanismSpec)
    cm = compile_mechanism(spec.base)
    rxn = spec.base.reaction
    mets = Tuple(vcat(
        [s[1] for s in substrates(rxn)],
        [p[1] for p in products(rxn)],
        collect(regulators(rxn)),
    ))
    reg_sites = Tuple(
        (Tuple(group), mult) for (group, mult) in zip(
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities))
    OligomericEnzymeMechanism{
        mets, typeof(cm), spec.catalytic_n,
        reg_sites, 2,  # n_conf always 2
    }()
end
```

**Step 2: Commit**

```bash
git commit -m "Add compile_mechanism for new MechanismSpec types"
```

---

## Task 16: Test specs and stage tests

**Files:**
- Create: `test/mechanism_enumeration_test_specs.jl` — test spec definitions
- Create: `test/test_mechanism_enumeration.jl` — stage-by-stage and end-to-end tests

**Step 1: Create `test/mechanism_enumeration_test_specs.jl`**

```julia
# Test specifications for new mechanism enumeration pipeline

struct StageExpansionTestSpec
    name::String
    reaction::Any
    base_mechanism::Any  # MechanismSpec (from new types)

    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int
    expected_n_tr_equiv::Int
    expected_n_oem_dedup::Int
end

struct EnumerationTestSpec2  # "2" to avoid conflict with old
    name::String
    reaction::Any

    expected_n_catalytic::Int
    expected_n_ress::Int
    expected_n_general_modifier::Int
    expected_n_essential_activator::Int
    expected_n_dead_end::Int
    expected_n_equivalence::Int
    expected_n_dedup::Int
    expected_n_allosteric::Int
    expected_n_tr_equiv::Int
    expected_n_oem_dedup::Int
end

function build_stage_expansion_test_specs()
    specs = StageExpansionTestSpec[]

    # Build specs for each reaction using hand-built base mechanisms
    # Expected counts will be filled in after running the pipeline
    # once and verifying by hand / orthogonal computation

    # 1. Uni-Uni (no regulators)
    # 2. Uni-Uni + 1 unconstrained regulator
    # 3. Uni-Bi + 1 regulator
    # 4. Bi-Bi + 2 regulators
    # 5. Bi-Bi Ping-Pong (no regulators)
    # 6. Bi-Bi Ping-Pong + 1 regulator

    # ... (detailed specs with hard-coded expected counts) ...

    return specs
end

function build_enumeration_test_specs2()
    specs = EnumerationTestSpec2[]
    # ... same reactions as above ...
    return specs
end
```

**Step 2: Create `test/test_mechanism_enumeration.jl`**

```julia
@testset "New Mechanism Enumeration" begin
    @testset "Stage Expansion: $(spec.name)" for spec in STAGE_TEST_SPECS
        # Stage 2: RE/SS
        ress = EnzymeRates._expand_ress_variants(
            [spec.base_mechanism], spec.reaction)
        @test length(ress) == spec.expected_n_ress

        # Stage 3: General modifier
        gm = EnzymeRates._expand_general_modifiers(ress, spec.reaction)
        @test length(gm) == spec.expected_n_general_modifier

        # ... etc for each stage ...

        # Parameter count verification (sample 10)
        rng = Random.MersenneTwister(42)
        sample = rand(rng, ress, min(10, length(ress)))
        for s in sample
            m = compile_mechanism(s)
            expected_params = length(parameters(m))
            @test s.param_count == expected_params
        end
    end

    @testset "End-to-End: $(spec.name)" for spec in ENUM_TEST_SPECS2
        result = enumerate_mechanisms(spec.reaction)
        # ... verify counts at each stage ...
    end
end
```

**Step 3: Commit**

```bash
git commit -m "Add test specs and tests for new mechanism enumeration"
```

---

## Task 17: Swap old for new

Replace old implementation with new one.

**Files:**
- Modify: `src/EnzymeRates.jl` — remove old include, keep new include
- Delete: `src/old_mechanism_enumeration.jl` (or keep for reference)
- Modify: `test/runtests.jl` — include new test files

**Step 1: Update `src/EnzymeRates.jl`**

Remove the `include("old_mechanism_enumeration.jl")` line. The `include("mechanism_enumeration.jl")` line should already be present.

**Step 2: Update `test/runtests.jl`**

Add includes for new test spec and test files. Remove or comment out old test includes.

**Step 3: Run full test suite**

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Expected: all new tests pass. Old tests are no longer included.

**Step 4: Commit**

```bash
git commit -m "Swap old mechanism enumeration for new staged pipeline"
```

---

## Task 18: Update CLAUDE.md

**Files:**
- Modify: `.claude/CLAUDE.md` — update architecture docs for new pipeline

Update the following sections:
- Source Layout — update `src/mechanism_enumeration.jl` description
- Key Architecture Decisions — add staged pipeline, regulator roles
- Testing — update test file descriptions

**Step 1: Update CLAUDE.md**

**Step 2: Commit**

```bash
git commit -m "Update CLAUDE.md for new mechanism enumeration architecture"
```

---

## Parallelization Notes

Tasks that can be executed in parallel (no dependencies between them):

- **Tasks 4, 5, 6, 7, 8** (stages 1-5) can be developed in parallel after Task 3 (shared infrastructure) is complete, since each stage function is independent. However, they must be composed sequentially for testing.

- **Tasks 9, 10** (equivalence constraints + dedup) depend on stages 1-5 being complete.

- **Tasks 11, 12, 13** (OEM stages 8-10) can be developed in parallel with tasks 9-10 since they operate on different types.

- **Task 14** (pipeline orchestration) depends on all stages being complete.

- **Task 15** (compile_mechanism) is independent of pipeline stages and can be done in parallel with tasks 4-13.

- **Task 16** (tests) should be developed incrementally alongside each stage task, but the spec values can only be finalized after the implementation is working.

Recommended execution order:
1. Task 0 (rename) — first
2. Task 1 (RegulatorRole + EnzymeReaction) — after Task 0
3. Task 2 (new types) — parallel with Task 1
4. Task 3 (shared infrastructure) — after Tasks 1 + 2
5. Tasks 4-8 (stages 1-5) — sequentially after Task 3
6. Tasks 9-10 (stages 6-7) — after Task 8
7. Tasks 11-13 (stages 8-10) — after Task 10
8. Task 14 (pipeline) — after all stages
9. Task 15 (compile_mechanism) — parallel with tasks 4-13
10. Task 16 (tests) — after Task 14
11. Task 17 (swap) — after Task 16
12. Task 18 (docs) — after Task 17
