# ─── Data Types ─────────────────────────────────────────────

"""
    SiteState

State of a single binding site on an enzyme form.

# Fields
- `metabolite::Symbol`: which metabolite binds here (e.g., `:S`, `:P`, `:R`).
- `atoms::Union{Nothing, Vector{Pair{Symbol,Int}}}`: atoms currently at
  this site (`nothing` = empty).
- `role::Symbol`: `:sub`, `:prod`, or `:reg`.
- `full_atoms::Vector{Pair{Symbol,Int}}`: atoms when the metabolite is
  fully bound (always populated, even when `atoms === nothing`).
"""
const SiteState = @NamedTuple{
    metabolite::Symbol,
    atoms::Union{Nothing, Vector{Pair{Symbol,Int}}},
    role::Symbol,
    full_atoms::Vector{Pair{Symbol,Int}},
}

"""
    EnzymeFormSpec

A named enzyme form with its binding site states. For example, `E_S_0`
has substrate S bound and product site empty.

# Fields
- `name::Symbol`: display name (e.g., `:E_S_0`).
- `sites::Vector{SiteState}`: one entry per binding site.
"""
struct EnzymeFormSpec
    name::Symbol
    sites::Vector{SiteState}
end

"""
Constraint on kinetic parameters: target parameter equals a linear
combination of source parameters.
Format: `(target_sym, coeff, [(src_sym, src_coeff), ...])`.
"""
const ParamConstraint = Tuple{Symbol, Int, Vector{Tuple{Symbol, Int}}}

"""
    MechanismSpec

A concrete mechanism topology: a reaction, a graph of enzyme forms
(as edge pairs into the `forms` vector), equilibrium step flags, and
parameter constraints between equivalent steps.

# Fields
- `reaction`: the `EnzymeReaction` this mechanism belongs to.
- `edges::Vector{Tuple{Int,Int}}`: edges between form indices.
- `equilibrium_steps::Vector{Bool}`: `true` for rapid-equilibrium steps.
- `param_constraints::Vector{ParamConstraint}`: equivalence constraints.
"""
struct MechanismSpec
    reaction::Any
    edges::Vector{Tuple{Int,Int}}
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{ParamConstraint}
end
MechanismSpec(reaction, edges) =
    MechanismSpec(reaction, edges, fill(false, length(edges)), ParamConstraint[])

"""
    EnumerationStage

Controls which pipeline stage `enumerate_mechanisms` stops at:
- `Catalytic()`: catalytic cycle topologies only.
- `WithActivator()`: add activator (essential/non-essential) variants.
- `WithDeadEnd()`: add dead-end inhibitor configurations.
- `FullEnumeration()`: add RE/SS + parameter constraint variants (lazy iterator).
"""
abstract type EnumerationStage end
struct Catalytic    <: EnumerationStage end
struct WithActivator <: EnumerationStage end
struct WithDeadEnd   <: EnumerationStage end
struct FullEnumeration <: EnumerationStage end

"""
    MechanismIterator

Lazy iterator over `MechanismSpec` instances. Wraps an inner iterator
with a precomputed O(1) total count.
"""
struct MechanismIterator
    inner::Any
    total::Int
end

Base.eltype(::Type{MechanismIterator}) = MechanismSpec
Base.IteratorSize(::Type{MechanismIterator}) = Base.HasLength()
Base.length(iter::MechanismIterator) = iter.total
Base.iterate(iter::MechanismIterator, s...) = iterate(iter.inner, s...)

# ─── Edge Classification + Adjacency ─────────────────────────

const EdgeInfo = @NamedTuple{type::Symbol, metabolite::Union{Nothing,Symbol}}

"""
    _classify_edge(form_a, form_b) → EdgeInfo or nothing

Classify the elementary step between two enzyme forms by comparing their
binding sites. Returns `nothing` if no valid step exists.

- **Single-site diff**: binding (empty→occupied) or release (occupied→empty).
- **Multi-site diff**: isomerization if substrate+product sites change
  with conserved atom balance (includes ping-pong residual transfers).
"""
function _classify_edge(form_a::EnzymeFormSpec, form_b::EnzymeFormSpec)
    diffs = [k for k in eachindex(form_a.sites)
             if form_a.sites[k].atoms != form_b.sites[k].atoms]
    isempty(diffs) && return nothing
    if length(diffs) == 1
        k = diffs[1]
        site_a, site_b = form_a.sites[k], form_b.sites[k]
        if site_a.atoms === nothing
            # Binding: substrate/product must bind with full atoms
            site_a.role in (:sub, :prod) && site_b.atoms != site_a.full_atoms &&
                return nothing
            return (type=:binding, metabolite=site_a.metabolite)
        elseif site_b.atoms === nothing
            # Release: determine which metabolite is leaving.
            # Full atoms or regulator → release own metabolite.
            # Partial (residual) → find the product whose full_atoms match.
            met = if site_a.atoms == site_a.full_atoms || site_a.role == :reg
                site_a.metabolite
            else
                idx = findfirst(
                    s -> s.role == :prod && s.full_atoms == site_a.atoms,
                    form_a.sites,
                )
                idx !== nothing ? form_a.sites[idx].metabolite : nothing
            end
            return met !== nothing ? (type=:release, metabolite=met) : nothing
        end
        return nothing
    end
    # Isomerization: multiple sub/prod sites change with conserved atoms
    all(k -> form_a.sites[k].role in (:sub, :prod), diffs) || return nothing
    any(k -> form_a.sites[k].role == :sub, diffs) &&
        any(k -> form_a.sites[k].role == :prod, diffs) || return nothing
    has_residual = any(k -> form_a.sites[k].role == :sub && any(
        x -> x !== nothing && x != form_a.sites[k].full_atoms,
        (form_a.sites[k].atoms, form_b.sites[k].atoms)), diffs)
    (has_residual || length(diffs) == count(
        s -> s.role in (:sub, :prod), form_a.sites)) || return nothing
    # Atom balance: net change across differing sites must be zero
    delta = Dict{Symbol,Int}()
    for k in diffs
        for (v, sign) in ((form_a.sites[k].atoms, 1),
                          (form_b.sites[k].atoms, -1))
            v === nothing && continue
            for (a, c) in v
                delta[a] = get(delta, a, 0) + sign * c
            end
        end
    end
    all(iszero, values(delta)) || return nothing
    (type=:isomerization, metabolite=nothing)
end

"""
    _build_adjacency(forms) → Dict{Tuple{Int,Int}, EdgeInfo}

Build undirected adjacency from all valid elementary steps between forms.
Keys are `(i, j)` with `i < j`.
"""
function _build_adjacency(forms::Vector{EnzymeFormSpec})
    Dict{Tuple{Int,Int},EdgeInfo}(
        (i, j) => info
        for i in eachindex(forms) for j in (i+1):length(forms)
        for info in (_classify_edge(forms[i], forms[j]),)
        if info !== nothing)
end

# ─── Enzyme Form Enumeration ─────────────────────────────────

"""
    enumerate_enzyme_forms(reaction::EnzymeReaction)

Enumerate all possible enzyme forms for the given reaction.

Each form is a Cartesian-product combination of per-site occupancy states:
empty, fully bound, or (for substrates in ping-pong mechanisms) partially
bound after product atom removal.
"""
function enumerate_enzyme_forms(reaction::EnzymeReaction{S,P,R}) where {S,P,R}
    sorted_atoms(spec) = sort([a => c for (a, c) in spec[2]]; by=first)
    atoms_label(v) = join(
        string(s) * (c > 1 ? string(c) : "") for (s, c) in v)
    # Product atom dicts for computing ping-pong residuals
    product_atom_dicts = [Dict{Symbol,Int}(a => c for (a, c) in p[2])
                          for p in P if !isempty(p[2])]
    metabolites = Symbol[]
    full_atom_lists = Vector{Pair{Symbol,Int}}[]
    roles = Symbol[]
    per_site_options =
        Vector{Tuple{Union{Nothing,Vector{Pair{Symbol,Int}}},String}}[]
    for (group, role) in ((S, :sub), (P, :prod), (R, :reg)), spec in group
        fa = sorted_atoms(spec)
        push!(metabolites, spec[1])
        push!(full_atom_lists, fa)
        push!(roles, role)
        site_options = [(nothing, "0"), (fa, string(spec[1]))]
        # Ping-pong residuals: subtract product atoms from substrate
        # (e.g., substrate A[CX] minus product P[C] → residual [X])
        if role == :sub && !isempty(spec[2]) && !isempty(product_atom_dicts)
            substrate_atoms = Dict{Symbol,Int}(fa)
            for m in 1:(1 << length(product_atom_dicts)) - 1
                pa = reduce(mergewith(+), (product_atom_dicts[i]
                    for i in 1:length(product_atom_dicts)
                    if (m >> (i-1)) & 1 == 1))
                all(get(substrate_atoms, a, 0) >= c
                    for (a, c) in pa) || continue
                residual = sort(
                    [a => c - get(pa, a, 0)
                     for (a, c) in substrate_atoms
                     if c > get(pa, a, 0)]; by=first)
                !isempty(residual) && residual != fa &&
                    !any(o -> o[1] == residual, site_options) &&
                    push!(site_options, (residual, atoms_label(residual)))
            end
        end
        push!(per_site_options, site_options)
    end
    # Cartesian product with exclusion filter:
    # Skip unreachable forms where all substrates are fully loaded
    # while any product is present (or vice versa for reverse direction).
    n = length(metabolites)
    forms = EnzymeFormSpec[]
    for combo in Iterators.product(per_site_options...)
        all_subs_full = all(
            i -> roles[i] != :sub || combo[i][1] == full_atom_lists[i],
            1:n)
        all_prods_full = all(
            i -> roles[i] != :prod || combo[i][1] == full_atom_lists[i],
            1:n)
        any_sub_occupied = any(
            i -> roles[i] == :sub && combo[i][1] !== nothing, 1:n)
        any_prod_occupied = any(
            i -> roles[i] == :prod && combo[i][1] !== nothing, 1:n)
        ((all_subs_full && any_prod_occupied) ||
            (all_prods_full && any_sub_occupied)) && continue
        sites = [(metabolite=metabolites[i], atoms=combo[i][1],
            role=roles[i], full_atoms=full_atom_lists[i]) for i in 1:n]
        push!(forms, EnzymeFormSpec(
            Symbol("E_" * join((c[2] for c in combo), "_")), sites))
    end
    forms
end

# ─── Catalytic Cycle Enumeration ──────────────────────────────

"""
    _catalytic_topologies(forms, adj, reaction; max_forms) → Vector{MechanismSpec}

Enumerate catalytic topologies: elementary cycles through the free enzyme
form that consume all substrates and release all products.

DFS finds elementary cycles, then power-set union generates all valid
multi-cycle topologies (e.g., random-order product release = union of
two sequential-release cycles).
"""
function _catalytic_topologies(
    forms::Vector{EnzymeFormSpec},
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    @nospecialize(reaction::EnzymeReaction);
    max_forms::Int,
)
    sub_set = Set(s[1] for s in substrates(reaction))
    prod_set = Set(p[1] for p in products(reaction))
    free = findfirst(f -> all(s -> s.atoms === nothing, f.sites), forms)
    free === nothing && return MechanismSpec[]
    # Restrict to forms without regulator binding (catalytic-only)
    cat_forms = Set(i for (i, f) in enumerate(forms)
        if all(s -> s.role != :reg || s.atoms === nothing, f.sites))
    cat_adj = Dict((a, b) => info for ((a, b), info) in adj
                   if a ∈ cat_forms && b ∈ cat_forms)
    # Form predicates for cycle validation
    _has_residual(fi) = any(
        s -> s.role == :sub && s.atoms !== nothing &&
             s.atoms != s.full_atoms, forms[fi].sites)
    _is_pure_intermediate(fi) = all(
        s -> (s.role != :sub || s.atoms != s.full_atoms) &&
             (s.role != :prod || s.atoms === nothing), forms[fi].sites)
    _all_substrates_full(fi) = all(
        s -> s.role != :sub || s.atoms == s.full_atoms, forms[fi].sites)
    # DFS for elementary cycles through the free enzyme form.
    # Tracks substrate binding and product release to ensure completeness.
    cycles = Set{Set{Int}}()
    function dfs(cur, path, visited, bound, released, has_residual)
        if cur == free && length(path) > 1 &&
                bound == sub_set && released == prod_set
            push!(cycles, Set(path))
            return
        end
        for ((a, b), info) in cat_adj
            neighbor = a == cur ? b : b == cur ? a : nothing
            neighbor === nothing && continue
            neighbor ∈ visited && neighbor != free && continue
            # Derive directed edge type from undirected adjacency:
            # binding when traversing in stored (a→b) direction
            edge_type = info.type == :isomerization ? :isomerization :
                ((a == cur) == (info.type == :binding)) ? :binding :
                :release
            met = info.metabolite
            if edge_type == :binding
                met ∈ sub_set && met ∉ bound && !has_residual || continue
                push!(bound, met)
            elseif edge_type == :release
                met ∈ prod_set && met ∉ released || continue
                push!(released, met)
            else
                neighbor ∉ visited || continue
            end
            push!(path, neighbor)
            push!(visited, neighbor)
            next_residual = edge_type == :release ? false :
                edge_type == :isomerization ? _has_residual(neighbor) :
                has_residual
            dfs(neighbor, path, visited, bound, released, next_residual)
            delete!(visited, neighbor)
            pop!(path)
            edge_type == :binding && delete!(bound, met)
            edge_type == :release && delete!(released, met)
        end
    end
    dfs(free, [free], Set([free]), Set{Symbol}(), Set{Symbol}(), false)
    # Combine elementary cycles via power-set union, then filter by
    # max_forms and ping-pong validity
    n_cycles = length(cycles)
    n_cycles == 0 && return MechanismSpec[]
    cycle_list = collect(cycles)
    combined = unique!([union((cycle_list[i] for i in 1:n_cycles
        if (m >> (i - 1)) & 1 == 1)...) for m in 1:(1 << n_cycles) - 1])
    filter!(combined) do form_set
        length(form_set) > max_forms && return false
        residual_forms = [fi for fi in form_set if _has_residual(fi)]
        isempty(residual_forms) && return true
        any(_is_pure_intermediate, residual_forms) &&
            !any(_all_substrates_full, form_set)
    end
    [MechanismSpec(reaction,
        [(a, b) for ((a, b), _) in adj if a ∈ form_set && b ∈ form_set])
     for form_set in combined]
end

# ─── Regulator Expansion ─────────────────────────────────────

"""Find dead-end form index: `base` form with additional regulators bound
at `occupied_positions`. Returns form index or `nothing`."""
function _find_dead_end(form_lookup, base::EnzymeFormSpec, occupied_positions)
    key = Tuple(k in occupied_positions ? s.full_atoms : s.atoms
                for (k, s) in enumerate(base.sites))
    get(form_lookup, key, nothing)
end

"""Build per-regulator activator options for a given catalytic topology.

Returns a vector of `(added_edges, removed_edges)` tuples representing:
1. Regulator absent (no edges added or removed).
2. Non-essential activator (all topology forms mirrored with reg bound).
3. Essential activator (only entry edge from free enzyme to reg-bound).
"""
function _build_activator_options(reg, topo, forms, form_lookup, edges)
    opts = [(Tuple{Int,Int}[], Set{Tuple{Int,Int}}())]
    pos = findfirst(s -> s.metabolite == reg && s.role == :reg &&
        s.atoms === nothing, forms[first(topo)].sites)
    pos === nothing && return opts
    # Shadow pairs: (base_form, reg-bound_form) for each topology form
    shadow_pairs = [
        (base_idx, _find_dead_end(form_lookup, forms[base_idx], (pos,)))
        for base_idx in sort(collect(topo))]
    filter!(((_, shadow),) -> shadow !== nothing, shadow_pairs)
    length(shadow_pairs) < 2 && return opts
    shadow_map = Dict(shadow_pairs)
    # Mirror the topology edges into the activated layer
    mirrored_edges = [(a, b) for (a, b) in edges
                      if haskey(shadow_map, a) && haskey(shadow_map, b)]
    mirror_edges = [minmax(shadow_map[a], shadow_map[b])
                    for (a, b) in mirrored_edges]
    bind_edges = [minmax(b, s) for (b, s) in shadow_pairs]
    # Non-essential: all base↔shadow binding + mirrored topology
    push!(opts, ([bind_edges; mirror_edges], Set{Tuple{Int,Int}}()))
    # Essential: only the entry edge (free enzyme → activated), removing
    # original topology edges that are replaced by the shadow cycle
    entry = findfirst(((base_idx, _),) -> all(
        s -> s.role == :reg || s.atoms === nothing,
        forms[base_idx].sites), shadow_pairs)
    if entry !== nothing
        push!(opts,
            ([minmax(shadow_pairs[entry]...); mirror_edges],
             Set(mirrored_edges)))
    end
    opts
end

"""
    _expand_activators(specs, forms, adj, form_lookup, reaction; max_forms)

Expand catalytic topologies with activator regulator variants.

For each regulator, three options are considered via `_build_activator_options`:
absent, non-essential activator, or essential activator.
Returns all valid Cartesian-product combinations across regulators.
"""
function _expand_activators(
    specs::Vector{MechanismSpec},
    forms::Vector{EnzymeFormSpec},
    form_lookup,
    @nospecialize(reaction::EnzymeReaction);
    max_forms::Int,
)
    reg_names = [r[1] for r in regulators(reaction)]
    isempty(reg_names) && return specs
    result = MechanismSpec[]
    for spec in specs
        topo = Set(Iterators.flatten(spec.edges))
        per_reg = map(reg_names) do reg
            _build_activator_options(
                reg, topo, forms, form_lookup, spec.edges)
        end
        for combo in Iterators.product(per_reg...)
            removals = union((r for (_, r) in combo)...)
            merged = [e for e in spec.edges if e ∉ removals]
            for (add, _) in combo
                append!(merged, add)
            end
            topo_nodes = Set(Iterators.flatten(merged))
            length(topo_nodes) > max_forms && continue
            push!(result, MechanismSpec(reaction, merged))
        end
    end
    result
end

"""
    _expand_inhibitors(specs, forms, adj, form_lookup; max_forms)

Expand mechanism topologies with dead-end inhibitor configurations.

For each topology, regulators not already acting as activators are treated
as inhibitors. Each topology form independently gets a bitmask of which
inhibitors bind there, generating `(2^n_inhibitors)^n_topo_forms` configs.
"""
function _expand_inhibitors(
    specs::Vector{MechanismSpec},
    forms::Vector{EnzymeFormSpec},
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    form_lookup;
    max_forms::Int,
)
    result = MechanismSpec[]
    for spec in specs
        topo_nodes = Set(Iterators.flatten(spec.edges))
        topo_sorted = sort!(collect(topo_nodes))
        budget = max_forms - length(topo_sorted)
        # Activator positions: reg sites occupied in any topology form
        activator_positions = Set(
            k for fi in topo_sorted
            for (k, s) in enumerate(forms[fi].sites)
            if s.role == :reg && s.atoms !== nothing)
        # Inhibitor positions: unoccupied reg sites not used as activators
        inhibitor_positions = [
            k for (k, s) in enumerate(forms[topo_sorted[1]].sites)
            if s.role == :reg && s.atoms === nothing &&
               k ∉ activator_positions]
        n_inh = length(inhibitor_positions)
        existing = Set(minmax(e...) for e in spec.edges)
        # Precompute dead-end forms for each (topo_form_index, inh_mask)
        dead_end_lookup = Dict(
            (ti, mask) => fi2
            for (ti, fi) in enumerate(topo_sorted)
            for mask in 1:(1 << n_inh) - 1
            for fi2 in (_find_dead_end(form_lookup, forms[fi],
                Tuple(inhibitor_positions[j] for j in 1:n_inh
                      if (mask >> (j - 1)) & 1 == 1)),)
            if fi2 !== nothing && fi2 ∉ topo_nodes)
        # Each topology form gets an independent inhibitor binding mask
        for dead_end_masks in Iterators.product(
                ntuple(_ -> 0:(1 << n_inh) - 1, length(topo_sorted))...)
            dead_end_forms = Set(
                fi2 for ((ti, m), fi2) in dead_end_lookup
                if m & dead_end_masks[ti] == m)
            length(dead_end_forms) > budget && continue
            all_forms = union(topo_nodes, dead_end_forms)
            new_edges = [
                (a, b) for ((a, b), info) in adj
                if a ∈ all_forms && b ∈ all_forms &&
                   (a, b) ∉ existing && info.type == :binding]
            push!(result, MechanismSpec(
                spec.reaction, [spec.edges; new_edges]))
        end
    end
    result
end

# ─── RE/SS + Constraint Lazy Generator ────────────────────────

"""
    _find_equivalent_groups(edges, adj, forms) → Vector{Vector{Int}}

Find groups of binding edges that involve the same non-product metabolite.
Edges within a group can be constrained to share kinetic parameters.
"""
function _find_equivalent_groups(edges, adj, forms)
    groups = Dict{Symbol, Vector{Int}}()
    product_metabolites = Set(
        s.metabolite for s in forms[1].sites if s.role == :prod)
    for (i, (a, b)) in enumerate(edges)
        info = get(adj, minmax(a, b), nothing)
        info !== nothing && info.type in (:binding, :release) &&
            info.metabolite ∉ product_metabolites &&
            push!(get!(groups, info.metabolite, Int[]), i)
    end
    sort!([sort(v) for v in values(groups) if length(v) >= 2]; by=first)
end

"""
    _ress_variants(spec, adj, forms)

Generate RE/SS (rapid-equilibrium / steady-state) + parameter constraint
variants for a mechanism.

Iterates over `2^n - 1` RE/SS masks (excluding all-RE), and for each
valid mask, over constraint combinations on equivalent step groups.
"""
function _ress_variants(spec, adj, forms)
    edges = spec.edges
    n = length(edges)
    equiv_groups = _find_equivalent_groups(edges, adj, forms)
    Iterators.flatmap(0:((1 << n) - 2)) do re_mask
        eq_steps = Bool[(re_mask >> (i-1)) & 1 == 1 for i in 1:n]
        # Valid groups: all edges in group share the same RE/SS assignment
        valid_groups = [g for g in equiv_groups
            if all(eq_steps[s] == eq_steps[g[1]] for s in g)]
        Iterators.map(0:((1 << length(valid_groups)) - 1)) do constraint_mask
            constraints = ParamConstraint[
                (Symbol("$p$(g[j])$s"), 1, [(Symbol("$p$(g[1])$s"), 1)])
                for (gi, g) in enumerate(valid_groups)
                if (constraint_mask >> (gi-1)) & 1 == 1
                for (p, ss) in ((eq_steps[g[1]] ?
                    ("K", ("",)) : ("k", ("f", "r"))),)
                for j in 2:length(g) for s in ss]
            MechanismSpec(spec.reaction, edges, eq_steps, constraints)
        end
    end
end

# ─── Pipeline ─────────────────────────────────────────────────

"""
    enumerate_mechanisms(reaction; stage=FullEnumeration(), max_forms=...)

Enumerate valid mechanism topologies for the given reaction.

Pipeline stages:
1. `Catalytic()` — catalytic cycles through the free enzyme.
2. `WithActivator()` — activator variants (essential/non-essential).
3. `WithDeadEnd()` — dead-end inhibitor configurations.
4. `FullEnumeration()` — RE/SS + parameter constraint variants (lazy `MechanismIterator`).
"""
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    stage::EnumerationStage=FullEnumeration(),
    max_forms::Int=3 * (length(substrates(reaction)) +
                        length(products(reaction)) +
                        length(regulators(reaction))),
)
    forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(forms)
    form_lookup = Dict(Tuple(s.atoms for s in f.sites) => i
                       for (i, f) in enumerate(forms))

    catalytic = _catalytic_topologies(forms, adj, reaction; max_forms)
    stage isa Catalytic && return catalytic

    with_activators = _expand_activators(
        catalytic, forms, form_lookup, reaction; max_forms)
    stage isa WithActivator && return with_activators

    with_inhibitors = _expand_inhibitors(
        with_activators, forms, adj, form_lookup; max_forms)
    stage isa WithDeadEnd && return with_inhibitors

    # FullEnumeration: compute RE/SS variant count, wrap in lazy iterator.
    # Count formula: 2^(n - Σgᵢ) × ∏(2^gᵢ + 2) - 2^k
    total = sum(with_inhibitors; init=0) do s
        n = length(s.edges)
        n == 0 && return 0
        equiv_groups = _find_equivalent_groups(s.edges, adj, forms)
        n_groups = length(equiv_groups)
        total_group_size = sum(length, equiv_groups; init=0)
        (1 << (n - total_group_size)) *
            prod(1 << length(g) + 2
                 for g in equiv_groups; init=1) -
            (1 << n_groups)
    end
    inner = Iterators.flatmap(
        s -> _ress_variants(s, adj, forms), with_inhibitors)
    MechanismIterator(inner, total)
end

# ─── MechanismSpec → EnzymeMechanism Conversion ──────────────

"""
    EnzymeMechanism(spec::MechanismSpec)

Convert a `MechanismSpec` to a type-parameterized `EnzymeMechanism`.
"""
function EnzymeMechanism(spec::MechanismSpec)
    rxn = spec.reaction
    forms = enumerate_enzyme_forms(rxn)
    adj = _build_adjacency(forms)
    used = Set(Iterators.flatten(spec.edges))
    _name_atoms(g) = Tuple((n, a) for (n, a, _) in g)
    _form_atoms(sites) = sort!([a => c for (a, c) in reduce(mergewith(+),
        (Dict(s.atoms) for s in sites if s.atoms !== nothing);
        init=Dict{Symbol,Int}())]; by=first)
    species = (
        _name_atoms(substrates(rxn)),
        _name_atoms(products(rxn)),
        _name_atoms(regulators(rxn)),
        Tuple((forms[i].name, Tuple(Tuple.(_form_atoms(forms[i].sites))))
            for i in eachindex(forms) if i ∈ used))
    reactions = Tuple(let info = adj[minmax(a, b)]
        ((info.type == :binding ? (forms[a].name, info.metabolite) : (forms[a].name,)),
         (info.type == :release ? (forms[b].name, info.metabolite) : (forms[b].name,)))
    end for (a, b) in spec.edges)
    EnzymeMechanism(species, reactions, Tuple(spec.equilibrium_steps),
        Tuple((t, c, Tuple(Tuple.(f))) for (t, c, f) in spec.param_constraints))
end
