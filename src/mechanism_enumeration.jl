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
(as edge pairs into the `forms` vector), equilibrium step flags,
parameter constraints between equivalent steps, regulator partition
info, and oligomeric expansion parameters.

# Fields
- `reaction`: the `EnzymeReaction` this mechanism belongs to.
- `edges::Vector{Tuple{Int,Int}}`: edges between form indices.
  Catalytic edges come first (`edges[1:n_catalytic_edges]`),
  followed by dead-end binding edges.
- `n_catalytic_edges::Int`: number of catalytic edges at the start
  of `edges`. Dead-end edges (`edges[n_catalytic_edges+1:end]`)
  inherit their RE/SS status from the catalytic edge connecting the
  same forms with regulator sites stripped. Regulator-binding
  dead-end edges (no catalytic counterpart) remain always RE.
- `equilibrium_steps::Vector{Bool}`: `true` for rapid-equilibrium steps.
- `param_constraints::Vector{ParamConstraint}`: equivalence constraints.
- `dead_end_regulators::Vector{Symbol}`: regulators treated as dead-end
  inhibitors in this partition.
- `allosteric_regulators::Vector{Symbol}`: regulators treated as
  allosteric (for OligomericEnzymeMechanism) in this partition.
- `catalytic_n::Int`: 0 = EnzymeMechanism, >0 = OligomericEnzymeMechanism
  with this many catalytic sites.
- `n_conf::Int`: number of conformational states (1 or 2).
- `allosteric_multiplicities::Vector{Int}`: multiplicity for each
  allosteric regulator (one entry per `allosteric_regulators` element).
"""
struct MechanismSpec
    reaction::Any
    edges::Vector{Tuple{Int,Int}}
    n_catalytic_edges::Int
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{ParamConstraint}
    dead_end_regulators::Vector{Symbol}
    allosteric_regulators::Vector{Symbol}
    catalytic_n::Int
    n_conf::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
end
MechanismSpec(reaction, edges) =
    MechanismSpec(reaction, edges, length(edges),
        fill(false, length(edges)),
        ParamConstraint[], Symbol[], Symbol[], 0, 1,
        Vector{Symbol}[], Int[])
MechanismSpec(reaction, edges, eq_steps, constraints) =
    MechanismSpec(reaction, edges, length(edges), eq_steps,
        constraints, Symbol[], Symbol[], 0, 1,
        Vector{Symbol}[], Int[])
MechanismSpec(reaction, edges, n_cat, eq_steps, constraints,
              dead_end_regs, allosteric_regs) =
    MechanismSpec(reaction, edges, n_cat, eq_steps,
        constraints, dead_end_regs, allosteric_regs, 0, 1,
        Vector{Symbol}[], Int[])

"""
    EnumerationStage

Controls which pipeline stage `enumerate_mechanisms` stops at:
- `Catalytic()`: catalytic cycle topologies only.
- `WithDeadEnd()`: add dead-end inhibitor configurations.
- `FullEnumeration()`: add RE/SS + parameter constraint variants (lazy iterator).
"""
abstract type EnumerationStage end
struct Catalytic    <: EnumerationStage end
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
    # Substrates and products: (name, atoms) tuples
    for (group, role) in ((S, :sub), (P, :prod)), spec in group
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
    # Regulators: plain Symbols (no atoms)
    for name in R
        push!(metabolites, name)
        push!(full_atom_lists, Pair{Symbol,Int}[])
        push!(roles, :reg)
        push!(per_site_options,
            [(nothing, "0"), (Pair{Symbol,Int}[], string(name))])
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

"""
    _expand_inhibitors(specs, forms, adj, form_lookup;
                       max_forms, dead_end_regs, allosteric_regs)

Expand mechanism topologies with dead-end inhibitor configurations.

Only regulators in `dead_end_regs` are used for dead-end expansion.
Each topology form independently gets a bitmask of which dead-end
regulators bind there, generating `(2^n_dead_end)^n_topo_forms`
configs. The resulting MechanismSpecs carry partition info.
"""
function _expand_inhibitors(
    specs::Vector{MechanismSpec},
    forms::Vector{EnzymeFormSpec},
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    form_lookup;
    max_forms::Int,
    dead_end_regs::AbstractVector{Symbol}=Symbol[],
    allosteric_regs::AbstractVector{Symbol}=Symbol[],
)
    result = MechanismSpec[]
    for spec in specs
        topo_nodes = Set(Iterators.flatten(spec.edges))
        topo_sorted = sort!(collect(topo_nodes))
        budget = max_forms - length(topo_sorted)
        # Only use regulator positions matching dead-end regulators
        reg_positions = [
            k for (k, s) in enumerate(forms[topo_sorted[1]].sites)
            if s.role == :reg && s.atoms === nothing &&
               s.metabolite in dead_end_regs]
        n_inh = length(reg_positions)
        existing = Set(minmax(e...) for e in spec.edges)
        # Precompute dead-end forms for each (topo_form_index, inh_mask)
        dead_end_lookup = Dict(
            (ti, mask) => fi2
            for (ti, fi) in enumerate(topo_sorted)
            for mask in 1:(1 << n_inh) - 1
            for fi2 in (_find_dead_end(form_lookup, forms[fi],
                Tuple(reg_positions[j] for j in 1:n_inh
                      if (mask >> (j - 1)) & 1 == 1)),)
            if fi2 !== nothing && fi2 ∉ topo_nodes)
        # Each topology form gets an independent inhibitor binding mask
        for dead_end_masks in Iterators.product(
                ntuple(_ -> 0:(1 << n_inh) - 1,
                    length(topo_sorted))...)
            dead_end_forms = Set(
                fi2 for ((ti, m), fi2) in dead_end_lookup
                if m & dead_end_masks[ti] == m)
            length(dead_end_forms) > budget && continue
            all_forms = union(topo_nodes, dead_end_forms)
            new_edges = [
                (a, b) for ((a, b), info) in adj
                if a ∈ all_forms && b ∈ all_forms &&
                   (a, b) ∉ existing && info.type == :binding]
            all_edges = [spec.edges; new_edges]
            push!(result, MechanismSpec(
                spec.reaction, all_edges,
                length(spec.edges),
                fill(false, length(all_edges)),
                ParamConstraint[],
                dead_end_regs, allosteric_regs,
                0, 1, Vector{Symbol}[], Int[]))
        end
    end
    result
end

# ─── RE Group Count + Isomerization Detection ─────────────────

"""
    _compute_re_group_count(edges, eq_steps) → Int

Compute G, the number of RE groups (connected components when only
RE edges are considered). Each form starts as its own group; RE
edges merge the groups of their endpoints.
"""
function _compute_re_group_count(edges, eq_steps)
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
        eq_steps[idx] || continue  # only RE steps merge
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end
    length(Set(find(i) for i in form_indices))
end

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

# ─── Concentration Fingerprint ─────────────────────────────────

"""Increment the exponent of `met` in a sorted monomial."""
function _add_met(mono::MONO, met::Symbol)::MONO
    result = copy(mono)
    idx = findfirst(p -> p.first == met, result)
    if idx !== nothing
        result[idx] = met => result[idx].second + 1
    else
        push!(result, met => 1)
        sort!(result; by=first)
    end
    result
end

"""
Enumerate all spanning arborescences of a G-node directed graph
rooted at `root`, returning the set of concentration monomials.
Each arborescence selects one incoming edge per non-root node.
"""
function _spanning_arborescence_monomials(
    G::Int, R_conc::Dict{Tuple{Int,Int}, Set{MONO}}, root::Int,
)
    result = Set{MONO}()
    non_root = [g for g in 1:G if g != root]
    isempty(non_root) && return Set{MONO}([MONO()])
    # Recursive enumeration: for each non-root node, pick an
    # incoming edge and accumulate the product of edge monomials
    function enumerate_trees!(idx, current_mono)
        if idx > length(non_root)
            push!(result, current_mono)
            return
        end
        node = non_root[idx]
        # King-Altman: cofactor for root = sum of in-arborescences.
        # Each non-root node has an edge FROM itself TO a destination
        # (its parent, closer to root).
        for dst in 1:G
            dst == node && continue
            edge_monos = get(R_conc, (node, dst), nothing)
            edge_monos === nothing && continue
            for emono in edge_monos
                enumerate_trees!(idx + 1, _mono_mul(current_mono, emono))
            end
        end
    end
    enumerate_trees!(1, MONO())
    result
end

"""
    _concentration_fingerprint(edges, eq_steps, forms, adj)
    → Set{MONO}

Compute the set of concentration monomials that appear in the
rate equation denominator, using only topology (no King-Altman).
"""
function _concentration_fingerprint(edges, eq_steps, forms, adj)
    partition = _compute_re_partition(edges, eq_steps)
    G = length(partition)
    group_set = [Set(g) for g in partition]
    form_to_group = Dict(
        i => g for (g, grp) in enumerate(partition) for i in grp)

    # BFS alpha concentration monomials per form
    alpha_conc = Dict{Int, MONO}()
    for (g, group) in enumerate(partition)
        ref = group[1]
        alpha_conc[ref] = MONO()
        queue = [ref]
        while !isempty(queue)
            cur = popfirst!(queue)
            for (idx, (a, b)) in enumerate(edges)
                eq_steps[idx] || continue
                neighbor = if a == cur && b in group_set[g]
                    b
                elseif b == cur && a in group_set[g]
                    a
                else
                    nothing
                end
                (neighbor === nothing || haskey(alpha_conc, neighbor)) &&
                    continue
                info = adj[minmax(a, b)]
                parent_mono = alpha_conc[cur]
                child_mono = if info.metabolite !== nothing
                    is_binding = (info.type == :binding && a == cur) ||
                                 (info.type == :release && b == cur)
                    is_binding ? _add_met(parent_mono, info.metabolite) :
                                 copy(parent_mono)
                else
                    copy(parent_mono)
                end
                alpha_conc[neighbor] = child_mono
                push!(queue, neighbor)
            end
        end
    end

    # Sigma: set of alpha monomials per group
    sigma_conc = [Set{MONO}(alpha_conc[i] for i in group)
                  for group in partition]

    # R_conc: rate matrix concentration monomials between groups
    R_conc = Dict{Tuple{Int,Int}, Set{MONO}}()
    for (idx, (a, b)) in enumerate(edges)
        eq_steps[idx] && continue  # only SS steps
        info = adj[minmax(a, b)]
        g1, g2 = form_to_group[a], form_to_group[b]
        g1 == g2 && continue  # intra-group

        # Forward: a→b
        fwd_met = (info.type == :binding) ? info.metabolite : nothing
        fwd_mono = fwd_met !== nothing ?
            _add_met(alpha_conc[a], fwd_met) : copy(alpha_conc[a])
        push!(get!(R_conc, (g1, g2), Set{MONO}()), fwd_mono)

        # Reverse: b→a
        rev_met = (info.type == :release) ? info.metabolite : nothing
        rev_mono = rev_met !== nothing ?
            _add_met(alpha_conc[b], rev_met) : copy(alpha_conc[b])
        push!(get!(R_conc, (g2, g1), Set{MONO}()), rev_mono)
    end

    # Combine: fingerprint = union over all groups g of
    #   sigma_conc[g] * D_conc[g]
    fingerprint = Set{MONO}()
    for g in 1:G
        D_g = _spanning_arborescence_monomials(G, R_conc, g)
        for s in sigma_conc[g], d in D_g
            push!(fingerprint, _mono_mul(s, d))
        end
    end
    fingerprint
end

"""
    _constraint_descriptor(edges, adj, eq_steps, valid_groups,
                           constraint_mask)
    → Set{Tuple{Symbol, Symbol}}

Compute the constraint descriptor: the set of (metabolite, mode)
pairs for each constrained equivalence group.
"""
function _constraint_descriptor(edges, adj, eq_steps, valid_groups,
                                constraint_mask)
    descriptor = Set{Tuple{Symbol, Symbol}}()
    for (gi, g) in enumerate(valid_groups)
        (constraint_mask >> (gi - 1)) & 1 == 1 || continue
        met = adj[minmax(edges[g[1]]...)].metabolite
        mode = eq_steps[g[1]] ? :RE : :SS
        push!(descriptor, (met, mode))
    end
    descriptor
end

"""
    _find_first_isomerization(edges, adj) → Int

Find index of the first isomerization edge in canonical edge order.
Falls back to index 1 if none found.
"""
function _find_first_isomerization(edges, adj)
    for (i, (a, b)) in enumerate(edges)
        info = adj[minmax(a, b)]
        info.type == :isomerization && return i
    end
    return 1
end

# ─── Dead-End SS/RE Propagation ───────────────────────────────

"""
    _dead_end_catalytic_map(edges, n_cat, forms)
    → Vector{Union{Nothing, Int}}

For each dead-end edge (index > n_cat), return the index of the
catalytic edge connecting the same forms with regulator sites
stripped, or `nothing` for regulator-binding edges (always RE).
"""
function _dead_end_catalytic_map(edges, n_cat, forms)
    _strip_regs(fi) = Tuple(
        s.role == :reg ? nothing : s.atoms for s in forms[fi].sites)
    _key(a, b) = let sa = _strip_regs(a), sb = _strip_regs(b)
        hash(sa) <= hash(sb) ? (sa, sb) : (sb, sa)
    end
    cat_stripped = Dict(
        _key(a, b) => i
        for (i, (a, b)) in enumerate(edges) if i <= n_cat)
    [get(cat_stripped, _key(a, b), nothing)
     for (a, b) in @view edges[n_cat+1:end]]
end

"""
    _propagate_de_eq_steps!(eq_steps, n_cat, de_cat_map)

Set dead-end edge eq_steps by copying from their catalytic
counterpart. R-binding dead-end edges (`nothing` in map) keep
their initialized value (RE).
"""
function _propagate_de_eq_steps!(eq_steps, n_cat, de_cat_map)
    for (di, cat_idx) in enumerate(de_cat_map)
        cat_idx === nothing && continue
        eq_steps[n_cat + di] = eq_steps[cat_idx]
    end
end

# ─── RE/SS + Constraint Lazy Generator ────────────────────────

"""
    _find_equivalent_groups(edges, adj, forms, n_catalytic_edges)

Find groups of catalytic binding edges that involve the same
non-product metabolite. Only catalytic edges (first `n_catalytic_edges`)
are considered — dead-end edges don't participate in parameter
constraint enumeration.
"""
function _find_equivalent_groups(edges, adj, forms,
    n_catalytic_edges, de_cat_map=nothing)
    groups = Dict{Symbol, Vector{Int}}()
    product_metabolites = Set(
        s.metabolite for s in forms[1].sites if s.role == :prod)
    for i in 1:n_catalytic_edges
        (a, b) = edges[i]
        info = get(adj, minmax(a, b), nothing)
        info !== nothing && info.type in (:binding, :release) &&
            info.metabolite ∉ product_metabolites &&
            push!(get!(groups, info.metabolite, Int[]), i)
    end
    if de_cat_map !== nothing
        for (di, cat_idx) in enumerate(de_cat_map)
            cat_idx === nothing && continue
            edge_idx = n_catalytic_edges + di
            info = adj[minmax(edges[edge_idx]...)]
            info.metabolite === nothing && continue
            info.metabolite in product_metabolites && continue
            push!(get!(groups, info.metabolite, Int[]), edge_idx)
        end
    end
    sort!([sort(v) for v in values(groups) if length(v) >= 2];
        by=first)
end

"""
Build parameter constraints from valid_groups and constraint_mask.
Extracted helper for `_ress_variants`.
"""
function _build_constraints(valid_groups, eq_steps, constraint_mask,
                            edges)
    ParamConstraint[
        (Symbol("$p$(g[j])$s"), 1,
            [(Symbol("$p$(g[1])$s"), 1)])
        for (gi, g) in enumerate(valid_groups)
        if (constraint_mask >> (gi-1)) & 1 == 1
        for (p, ss) in ((eq_steps[g[1]] ?
            ("K", ("",)) : ("k", ("f", "r"))),)
        for j in 2:length(g) for s in ss]
end

const _DedupKey = Tuple{Set{MONO}, Set{Tuple{Symbol,Symbol}}}

"""
    _ress_variants(spec, adj, forms; max_re_groups=7)

Generate RE/SS + parameter constraint variants for a mechanism,
deduplicated by (concentration fingerprint, constraint descriptor).

Two-pass: first collects all variants and keeps the one with fewest
SS steps per dedup key, then emits the kept variants.
"""
function _ress_variants(spec, adj, forms; max_re_groups::Int=7)
    edges = spec.edges
    n = length(edges)
    n_cat = spec.n_catalytic_edges
    iso_idx = _find_first_isomerization(edges, adj)
    de_cat_map = _dead_end_catalytic_map(edges, n_cat, forms)
    equiv_groups = _find_equivalent_groups(edges, adj, forms,
        n_cat, de_cat_map)
    other_indices = [i for i in 1:n_cat if i != iso_idx]
    n_other = length(other_indices)

    # Sort masks by popcount so first seen has fewest SS steps
    masks = sort!(collect(0:(1 << n_other) - 1); by=count_ones)

    best = Dict{_DedupKey, MechanismSpec}()
    for ss_mask in masks
        eq_steps = fill(true, n)
        eq_steps[iso_idx] = false
        for (bit, idx) in enumerate(other_indices)
            (ss_mask >> (bit - 1)) & 1 == 1 &&
                (eq_steps[idx] = false)
        end
        _propagate_de_eq_steps!(eq_steps, n_cat, de_cat_map)
        G = _compute_re_group_count(edges, eq_steps)
        (G < 2 || G > max_re_groups) && continue
        fp = _concentration_fingerprint(edges, eq_steps, forms, adj)
        valid_groups = [g for g in equiv_groups
            if all(eq_steps[s] == eq_steps[g[1]] for s in g)]
        for constraint_mask in 0:(1 << length(valid_groups)) - 1
            desc = _constraint_descriptor(
                edges, adj, eq_steps, valid_groups, constraint_mask)
            key = (fp, desc)
            haskey(best, key) && continue
            constraints = _build_constraints(
                valid_groups, eq_steps, constraint_mask, edges)
            best[key] = MechanismSpec(
                spec.reaction, edges, n_cat, eq_steps,
                constraints, spec.dead_end_regulators,
                spec.allosteric_regulators, 0, 1,
                Vector{Symbol}[], Int[])
        end
    end
    values(best)
end

"""
    _count_ress_variants(spec, adj, forms; max_re_groups=7) → Int

Count deduplicated RE/SS + constraint variants.
Same dedup logic as `_ress_variants` but returns only the count.
"""
function _count_ress_variants(
    spec, adj, forms; max_re_groups::Int=7,
)
    edges = spec.edges
    n = length(edges)
    n == 0 && return 0
    n_cat = spec.n_catalytic_edges
    iso_idx = _find_first_isomerization(edges, adj)
    de_cat_map = _dead_end_catalytic_map(edges, n_cat, forms)
    equiv_groups = _find_equivalent_groups(edges, adj, forms,
        n_cat, de_cat_map)
    other_indices = [i for i in 1:n_cat if i != iso_idx]
    n_other = length(other_indices)

    seen = Set{_DedupKey}()
    for ss_mask in 0:(1 << n_other) - 1
        eq_steps = fill(true, n)
        eq_steps[iso_idx] = false
        for (bit, idx) in enumerate(other_indices)
            (ss_mask >> (bit - 1)) & 1 == 1 &&
                (eq_steps[idx] = false)
        end
        _propagate_de_eq_steps!(eq_steps, n_cat, de_cat_map)
        G = _compute_re_group_count(edges, eq_steps)
        (G < 2 || G > max_re_groups) && continue
        fp = _concentration_fingerprint(edges, eq_steps, forms, adj)
        valid_groups = [g for g in equiv_groups
            if all(eq_steps[s] == eq_steps[g[1]]
                   for s in g)]
        for constraint_mask in 0:(1 << length(valid_groups)) - 1
            desc = _constraint_descriptor(
                edges, adj, eq_steps, valid_groups,
                constraint_mask)
            push!(seen, (fp, desc))
        end
    end
    length(seen)
end

# ─── Pipeline ─────────────────────────────────────────────────

"""
    enumerate_mechanisms(reaction; stage, max_re_groups, catalytic_n)

Enumerate valid mechanism topologies for the given reaction.

Pipeline stages:
1. `Catalytic()` — catalytic cycles through the free enzyme.
2. `WithDeadEnd()` — dead-end inhibitor configurations for all
   regulator partitions (each regulator → dead-end or allosteric).
3. `FullEnumeration()` — RE/SS + parameter constraint variants
   (lazy `MechanismIterator`).

Regulator partitioning: for reactions with n regulators, all 2^n
partitions of regulators into {dead-end, allosteric} are enumerated.
Catalytic topologies (stage 1) are computed once and reused across
all partitions.

When `catalytic_n > 0`, also produces `OligomericEnzymeMechanism`
candidates (NConf=2) for each RE/SS variant. Each allosteric
regulator gets multiplicity ∈ 1:catalytic_n; all multiplicities are
enumerated via Cartesian product.

# Keywords
- `max_re_groups::Int=7`: maximum number of RE groups (King-Altman
  rate matrix dimension). Mechanisms with G > max_re_groups are
  filtered out during RE/SS expansion.
- `catalytic_n::Int=0`: when > 0, also enumerate
  `OligomericEnzymeMechanism` candidates with this many catalytic
  subunits and NConf=2.
"""
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    stage::EnumerationStage=FullEnumeration(),
    max_re_groups::Int=7,
    catalytic_n::Int=0,
)
    forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(forms)
    form_lookup = Dict(Tuple(s.atoms for s in f.sites) => i
                       for (i, f) in enumerate(forms))
    max_forms = length(forms)

    catalytic = _catalytic_topologies(forms, adj, reaction; max_forms)
    stage isa Catalytic && return catalytic

    # Regulator partitioning: enumerate all 2^n_reg partitions
    regs = collect(Symbol, regulators(reaction))
    n_reg = length(regs)

    all_inhibitor_specs = MechanismSpec[]
    for reg_mask in 0:(1 << n_reg) - 1
        dead_end_regs = Symbol[
            regs[i] for i in 1:n_reg
            if (reg_mask >> (i - 1)) & 1 == 0]
        allosteric_regs = Symbol[
            regs[i] for i in 1:n_reg
            if (reg_mask >> (i - 1)) & 1 == 1]
        partition_specs = _expand_inhibitors(
            catalytic, forms, adj, form_lookup;
            max_forms, dead_end_regs, allosteric_regs)
        append!(all_inhibitor_specs, partition_specs)
    end
    stage isa WithDeadEnd && return all_inhibitor_specs

    # FullEnumeration: count RE/SS variants (brute-force with G cap),
    # wrap in lazy iterator.
    em_total = sum(all_inhibitor_specs; init=0) do s
        _count_ress_variants(s, adj, forms; max_re_groups)
    end

    if catalytic_n > 0
        oem_total = sum(all_inhibitor_specs; init=0) do s
            _oligomeric_count(
                s,
                _count_ress_variants(s, adj, forms;
                    max_re_groups),
                catalytic_n)
        end
        total = em_total + oem_total
        inner = Iterators.flatmap(all_inhibitor_specs) do s
            ress = _ress_variants(s, adj, forms;
                max_re_groups)
            Iterators.flatmap(ress) do em_spec
                Iterators.flatten((
                    (em_spec,),
                    _expand_oligomeric_variants(
                        em_spec, catalytic_n),
                ))
            end
        end
    else
        total = em_total
        inner = Iterators.flatmap(
            s -> _ress_variants(s, adj, forms;
                max_re_groups),
            all_inhibitor_specs)
    end
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
    _form_atoms(sites) = sort!([a => c for (a, c) in reduce(mergewith(+),
        (Dict(s.atoms) for s in sites if s.atoms !== nothing);
        init=Dict{Symbol,Int}())]; by=first)
    species = (
        substrates(rxn),
        products(rxn),
        regulators(rxn),
        Tuple((forms[i].name, Tuple(Tuple.(_form_atoms(forms[i].sites))))
            for i in eachindex(forms) if i ∈ used))
    reactions = Tuple(let info = adj[minmax(a, b)]
        ((info.type == :binding ? (forms[a].name, info.metabolite) : (forms[a].name,)),
         (info.type == :release ? (forms[b].name, info.metabolite) : (forms[b].name,)))
    end for (a, b) in spec.edges)
    EnzymeMechanism(species, reactions, Tuple(spec.equilibrium_steps),
        Tuple((t, c, Tuple(Tuple.(f))) for (t, c, f) in spec.param_constraints))
end

"""
    compile_mechanism(spec::MechanismSpec)

Convert a `MechanismSpec` to its mechanism type:
- `catalytic_n == 0` → `EnzymeMechanism`
- `catalytic_n > 0` → `OligomericEnzymeMechanism` with `NConf=n_conf`
"""
function compile_mechanism(spec::MechanismSpec)
    cm = EnzymeMechanism(spec)
    spec.catalytic_n == 0 && return cm
    rxn = spec.reaction
    mets = Tuple(vcat(
        [s[1] for s in substrates(rxn)],
        [p[1] for p in products(rxn)],
        collect(regulators(rxn)),
    ))
    reg_sites = if isempty(spec.allosteric_reg_sites)
        Tuple(
            ((reg,), mult) for (reg, mult) in zip(
                spec.allosteric_regulators,
                spec.allosteric_multiplicities))
    else
        Tuple(
            (Tuple(group), mult) for (group, mult) in zip(
                spec.allosteric_reg_sites,
                spec.allosteric_multiplicities))
    end
    OligomericEnzymeMechanism{
        mets, typeof(cm), spec.catalytic_n,
        reg_sites, spec.n_conf,
    }()
end

# ─── Set Partitions + Oligomeric Expansion ─────────────────────

"""
    _set_partitions(elements::Vector{Symbol})

Enumerate all set partitions of `elements` (Bell number partitions).
Returns `Vector{Vector{Vector{Symbol}}}` — each partition is a list
of groups, each group is a list of symbols.
"""
function _set_partitions(elements::Vector{Symbol})
    n = length(elements)
    n == 0 && return [Vector{Symbol}[]]
    n == 1 && return [Vector{Symbol}[elements]]
    result = Vector{Vector{Vector{Symbol}}}()
    for partition in _set_partitions(elements[1:end-1])
        last_elem = elements[end]
        for i in eachindex(partition)
            new_part = [copy(g) for g in partition]
            push!(new_part[i], last_elem)
            push!(result, new_part)
        end
        push!(result, [partition; [Symbol[last_elem]]])
    end
    result
end

"""
    _partition_mult_count(k, N) → Int

Count OEM multiplicity variants across all set partitions of k
regulators: sum_{g=1}^{k} S(k,g) * N^g where S(k,g) is the
Stirling number of the second kind.
"""
function _partition_mult_count(k::Int, N::Int)
    k == 0 && return 1
    S = zeros(Int, k, k)
    S[1, 1] = 1
    for n in 2:k
        for g in 1:n
            S[n, g] = (g > 1 ? S[n-1, g-1] : 0) + g * S[n-1, g]
        end
    end
    sum(S[k, g] * N^g for g in 1:k)
end

"""
    _oligomeric_count(spec, ress_count, catalytic_n) → Int

Count OligomericEnzymeMechanism variants for one dead-end spec.
For k allosteric regulators with `catalytic_n=N`, enumerates all
set partitions of regulators into site groups, then multiplicity
combos per partition: `ress_count * _partition_mult_count(k, N)`.
"""
function _oligomeric_count(spec, ress_count, catalytic_n)
    k = length(spec.allosteric_regulators)
    ress_count * _partition_mult_count(k, catalytic_n)
end

"""
    _expand_oligomeric_variants(em_spec, catalytic_n)

For a single EnzymeMechanism RE/SS spec, generate the corresponding
OligomericEnzymeMechanism specs by enumerating set partitions of
allosteric regulators into site groups, then multiplicity combos
(each ∈ 1:catalytic_n) per site group.
"""
function _expand_oligomeric_variants(em_spec, catalytic_n)
    regs = em_spec.allosteric_regulators
    k = length(regs)
    if k == 0
        oem = MechanismSpec(
            em_spec.reaction, em_spec.edges,
            em_spec.n_catalytic_edges,
            em_spec.equilibrium_steps,
            em_spec.param_constraints,
            em_spec.dead_end_regulators,
            em_spec.allosteric_regulators,
            catalytic_n, 2, Vector{Symbol}[], Int[])
        return (oem,)
    end
    partitions = _set_partitions(regs)
    Iterators.flatmap(partitions) do partition
        n_groups = length(partition)
        Iterators.map(
            Iterators.product(
                ntuple(_ -> 1:catalytic_n, n_groups)...)
        ) do combo
            MechanismSpec(
                em_spec.reaction, em_spec.edges,
                em_spec.n_catalytic_edges,
                em_spec.equilibrium_steps,
                em_spec.param_constraints,
                em_spec.dead_end_regulators,
                em_spec.allosteric_regulators,
                catalytic_n, 2, partition, collect(combo))
        end
    end
end
