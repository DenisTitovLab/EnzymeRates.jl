# ─── Data Types ─────────────────────────────────────────────

"""
    SiteDefinition

Per-reaction binding site definition (shared across all enzyme forms).

# Fields
- `metabolite::Symbol`: which metabolite binds here (e.g., `:S`, `:P`, `:R`).
- `role::Symbol`: `:sub`, `:prod`, or `:reg`.
- `full_atoms::Vector{Pair{Symbol,Int}}`: atoms when fully bound.
"""
struct SiteDefinition
    metabolite::Symbol
    role::Symbol
    full_atoms::Vector{Pair{Symbol,Int}}
end

"""
    EnzymeFormSpec

A named enzyme form with its per-site occupancy states.

# Fields
- `name::Symbol`: display name (e.g., `:E_S_0`).
- `occupancy::Vector{Union{Nothing, Vector{Pair{Symbol,Int}}}}`:
  one entry per site — `nothing` = empty, otherwise the bound atoms.
"""
struct EnzymeFormSpec
    name::Symbol
    occupancy::Vector{Union{Nothing, Vector{Pair{Symbol,Int}}}}
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

# Adjacency values: Union{Nothing, Symbol}
#   nothing = isomerization (no metabolite exchanged)
#   Symbol  = metabolite exchanged in this step

"""
    _classify_edge(site_defs, form_a, form_b)

Classify the elementary step between two enzyme forms.
Returns `missing` if no valid step exists, `nothing` for
isomerization, or the metabolite `Symbol` for binding/release.
"""
function _classify_edge(
    site_defs::Vector{SiteDefinition},
    form_a::EnzymeFormSpec,
    form_b::EnzymeFormSpec,
)
    diffs = [k for k in eachindex(site_defs)
             if form_a.occupancy[k] != form_b.occupancy[k]]
    isempty(diffs) && return missing
    if length(diffs) == 1
        k = diffs[1]
        sd = site_defs[k]
        occ_a, occ_b = form_a.occupancy[k], form_b.occupancy[k]
        if occ_a === nothing
            # Binding: substrate/product must bind with full atoms
            sd.role in (:sub, :prod) && occ_b != sd.full_atoms &&
                return missing
            return sd.metabolite
        elseif occ_b === nothing
            # Release: determine which metabolite is leaving
            met = if occ_a == sd.full_atoms || sd.role == :reg
                sd.metabolite
            else
                idx = findfirst(
                    j -> site_defs[j].role == :prod &&
                         site_defs[j].full_atoms == occ_a,
                    eachindex(site_defs),
                )
                idx !== nothing ? site_defs[idx].metabolite : nothing
            end
            return met !== nothing ? met : missing
        end
        return missing
    end
    # Isomerization: multiple sub/prod sites change with conserved atoms
    all(k -> site_defs[k].role in (:sub, :prod), diffs) || return missing
    any(k -> site_defs[k].role == :sub, diffs) &&
        any(k -> site_defs[k].role == :prod, diffs) || return missing
    has_residual = any(k -> site_defs[k].role == :sub && any(
        x -> x !== nothing && x != site_defs[k].full_atoms,
        (form_a.occupancy[k], form_b.occupancy[k])), diffs)
    (has_residual || length(diffs) == count(
        sd -> sd.role in (:sub, :prod), site_defs)) || return missing
    # Atom balance: net change across differing sites must be zero
    delta = Dict{Symbol,Int}()
    for k in diffs
        for (v, sign) in ((form_a.occupancy[k], 1),
                          (form_b.occupancy[k], -1))
            v === nothing && continue
            for (a, c) in v
                delta[a] = get(delta, a, 0) + sign * c
            end
        end
    end
    all(iszero, values(delta)) || return missing
    nothing  # isomerization: no metabolite
end

"""
    _build_adjacency(site_defs, forms)

Build undirected adjacency from all valid elementary steps between
forms. Keys are `(i, j)` with `i < j`. Values are the metabolite
`Symbol` (for binding/release) or `nothing` (for isomerization).
"""
function _build_adjacency(
    site_defs::Vector{SiteDefinition},
    forms::Vector{EnzymeFormSpec},
)
    Dict{Tuple{Int,Int}, Union{Nothing, Symbol}}(
        (i, j) => met
        for i in eachindex(forms) for j in (i+1):length(forms)
        for met in (_classify_edge(site_defs, forms[i], forms[j]),)
        if !ismissing(met))
end

"""
    _is_binding_direction(forms, from, to)

Return true if traversing from form `from` to form `to` binds a
metabolite (first differing site goes empty→occupied).
Only valid for non-isomerization edges (single-site diff).
"""
function _is_binding_direction(forms::Vector{EnzymeFormSpec},
                               from::Int, to::Int)
    for k in eachindex(forms[from].occupancy)
        forms[from].occupancy[k] == forms[to].occupancy[k] && continue
        return forms[from].occupancy[k] === nothing
    end
    error("unreachable: no differing site")
end

# ─── Enzyme Form Enumeration ─────────────────────────────────

"""
    enumerate_enzyme_forms(reaction::EnzymeReaction)

Enumerate all possible enzyme forms for the given reaction.
Returns `(site_defs, forms)` where `site_defs` is a vector of
per-site definitions and `forms` is all valid occupancy combos.
"""
function enumerate_enzyme_forms(reaction::EnzymeReaction{S,P,R}) where {S,P,R}
    sorted_atoms(spec) = sort([a => c for (a, c) in spec[2]]; by=first)
    atoms_label(v) = join(
        string(s) * (c > 1 ? string(c) : "") for (s, c) in v)
    # Product atom dicts for computing ping-pong residuals
    product_atom_dicts = [Dict{Symbol,Int}(a => c for (a, c) in p[2])
                          for p in P if !isempty(p[2])]
    site_defs = SiteDefinition[]
    per_site_options =
        Vector{Tuple{Union{Nothing,Vector{Pair{Symbol,Int}}},String}}[]
    # Substrates and products: (name, atoms) tuples
    for (group, role) in ((S, :sub), (P, :prod)), spec in group
        fa = sorted_atoms(spec)
        push!(site_defs, SiteDefinition(spec[1], role, fa))
        site_options = [(nothing, "0"), (fa, string(spec[1]))]
        # Ping-pong residuals: subtract product atoms from substrate
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
        push!(site_defs, SiteDefinition(name, :reg, Pair{Symbol,Int}[]))
        push!(per_site_options,
            [(nothing, "0"), (Pair{Symbol,Int}[], string(name))])
    end
    # Cartesian product with exclusion filter:
    # Skip unreachable forms where all substrates are fully loaded
    # while any product is present (or vice versa for reverse direction).
    n = length(site_defs)
    forms = EnzymeFormSpec[]
    for combo in Iterators.product(per_site_options...)
        all_subs_full = all(
            i -> site_defs[i].role != :sub ||
                 combo[i][1] == site_defs[i].full_atoms, 1:n)
        all_prods_full = all(
            i -> site_defs[i].role != :prod ||
                 combo[i][1] == site_defs[i].full_atoms, 1:n)
        any_sub_occupied = any(
            i -> site_defs[i].role == :sub && combo[i][1] !== nothing,
            1:n)
        any_prod_occupied = any(
            i -> site_defs[i].role == :prod && combo[i][1] !== nothing,
            1:n)
        ((all_subs_full && any_prod_occupied) ||
            (all_prods_full && any_sub_occupied)) && continue
        push!(forms, EnzymeFormSpec(
            Symbol("E_" * join((c[2] for c in combo), "_")),
            [combo[i][1] for i in 1:n]))
    end
    (site_defs, forms)
end

# ─── Catalytic Cycle Enumeration ──────────────────────────────

"""
    _catalytic_topologies(site_defs, forms, adj, reaction; max_forms)

Enumerate catalytic topologies: elementary cycles through the free
enzyme form that consume all substrates and release all products.
"""
function _catalytic_topologies(
    site_defs::Vector{SiteDefinition},
    forms::Vector{EnzymeFormSpec},
    adj::Dict{Tuple{Int,Int}, Union{Nothing, Symbol}},
    @nospecialize(reaction::EnzymeReaction);
    max_forms::Int,
)
    sub_set = Set(s[1] for s in substrates(reaction))
    prod_set = Set(p[1] for p in products(reaction))
    free = findfirst(
        f -> all(o === nothing for o in f.occupancy), forms)
    free === nothing && return MechanismSpec[]
    # Restrict to forms without regulator binding (catalytic-only)
    cat_forms = Set(i for (i, f) in enumerate(forms)
        if all(k -> site_defs[k].role != :reg ||
                    f.occupancy[k] === nothing,
               eachindex(site_defs)))
    cat_adj = Dict(
        (a, b) => met for ((a, b), met) in adj
        if a ∈ cat_forms && b ∈ cat_forms)
    # Form predicates for cycle validation
    _has_residual(fi) = any(
        k -> site_defs[k].role == :sub &&
             forms[fi].occupancy[k] !== nothing &&
             forms[fi].occupancy[k] != site_defs[k].full_atoms,
        eachindex(site_defs))
    _is_pure_intermediate(fi) = all(
        k -> (site_defs[k].role != :sub ||
              forms[fi].occupancy[k] != site_defs[k].full_atoms) &&
             (site_defs[k].role != :prod ||
              forms[fi].occupancy[k] === nothing),
        eachindex(site_defs))
    _all_substrates_full(fi) = all(
        k -> site_defs[k].role != :sub ||
             forms[fi].occupancy[k] == site_defs[k].full_atoms,
        eachindex(site_defs))
    # DFS for elementary cycles through the free enzyme form.
    cycles = Set{Set{Int}}()
    function dfs(cur, path, visited, bound, released, has_res)
        if cur == free && length(path) > 1 &&
                bound == sub_set && released == prod_set
            push!(cycles, Set(path))
            return
        end
        for ((a, b), met) in cat_adj
            neighbor = a == cur ? b : b == cur ? a : nothing
            neighbor === nothing && continue
            neighbor ∈ visited && neighbor != free && continue
            # Derive directed edge semantics from form occupancy
            binding = met !== nothing &&
                _is_binding_direction(forms, cur, neighbor)
            releasing = met !== nothing && !binding
            if met === nothing  # isomerization
                neighbor ∉ visited || continue
            elseif binding
                met ∈ sub_set && met ∉ bound && !has_res ||
                    continue
                push!(bound, met)
            else  # release
                met ∈ prod_set && met ∉ released || continue
                push!(released, met)
            end
            push!(path, neighbor)
            push!(visited, neighbor)
            next_res = releasing ? false :
                met === nothing ? _has_residual(neighbor) : has_res
            dfs(neighbor, path, visited, bound, released,
                next_res)
            delete!(visited, neighbor)
            pop!(path)
            binding && delete!(bound, met)
            releasing && delete!(released, met)
        end
    end
    dfs(free, [free], Set([free]), Set{Symbol}(), Set{Symbol}(),
        false)
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

"""Find dead-end form index: `base` form with additional regulators
bound at `occupied_positions`. Returns form index or `nothing`."""
function _find_dead_end(form_lookup, site_defs, base::EnzymeFormSpec,
                        occupied_positions)
    key = ntuple(length(site_defs)) do k
        k in occupied_positions ?
            site_defs[k].full_atoms : base.occupancy[k]
    end
    get(form_lookup, key, nothing)
end

"""
    _expand_inhibitors(specs, site_defs, forms, adj, form_lookup;
                       max_forms, dead_end_regs, allosteric_regs)

Expand mechanism topologies with dead-end inhibitor configurations.
"""
function _expand_inhibitors(
    specs::Vector{MechanismSpec},
    site_defs::Vector{SiteDefinition},
    forms::Vector{EnzymeFormSpec},
    adj::Dict{Tuple{Int,Int}, Union{Nothing, Symbol}},
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
        reg_positions = [
            k for (k, sd) in enumerate(site_defs)
            if sd.role == :reg &&
               forms[topo_sorted[1]].occupancy[k] === nothing &&
               sd.metabolite in dead_end_regs]
        n_inh = length(reg_positions)
        existing = Set(minmax(e...) for e in spec.edges)
        # Precompute dead-end forms for each (topo_form_idx, inh_mask)
        dead_end_lookup = Dict(
            (ti, mask) => fi2
            for (ti, fi) in enumerate(topo_sorted)
            for mask in 1:(1 << n_inh) - 1
            for fi2 in (_find_dead_end(form_lookup, site_defs,
                forms[fi],
                Tuple(reg_positions[j] for j in 1:n_inh
                      if (mask >> (j - 1)) & 1 == 1)),)
            if fi2 !== nothing && fi2 ∉ topo_nodes)
        # Each topology form gets an independent inhibitor mask
        for dead_end_masks in Iterators.product(
                ntuple(_ -> 0:(1 << n_inh) - 1,
                    length(topo_sorted))...)
            dead_end_forms = Set(
                fi2 for ((ti, m), fi2) in dead_end_lookup
                if m & dead_end_masks[ti] == m)
            length(dead_end_forms) > budget && continue
            all_forms = union(topo_nodes, dead_end_forms)
            new_edges = [
                (a, b) for ((a, b), met) in adj
                if a ∈ all_forms && b ∈ all_forms &&
                   (a, b) ∉ existing && met !== nothing]
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

# ─── RE Partition (union-find) ─────────────────────────────────

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
"""
function _spanning_arborescence_monomials(
    G::Int, R_conc::Dict{Tuple{Int,Int}, Set{MONO}}, root::Int,
)
    result = Set{MONO}()
    non_root = [g for g in 1:G if g != root]
    isempty(non_root) && return Set{MONO}([MONO()])
    function enumerate_trees!(idx, current_mono)
        if idx > length(non_root)
            push!(result, current_mono)
            return
        end
        node = non_root[idx]
        for dst in 1:G
            dst == node && continue
            edge_monos = get(R_conc, (node, dst), nothing)
            edge_monos === nothing && continue
            for emono in edge_monos
                enumerate_trees!(
                    idx + 1, _mono_mul(current_mono, emono))
            end
        end
    end
    enumerate_trees!(1, MONO())
    result
end

"""
    _concentration_fingerprint(edges, eq_steps, site_defs, forms,
                               adj, partition) → Set{MONO}

Compute the set of concentration monomials that appear in the
rate equation denominator, using only topology (no King-Altman).
"""
function _concentration_fingerprint(
    edges, eq_steps, site_defs, forms, adj, partition,
)
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
                (neighbor === nothing ||
                    haskey(alpha_conc, neighbor)) && continue
                met = adj[minmax(a, b)]
                parent_mono = alpha_conc[cur]
                child_mono = if met !== nothing
                    binding = _is_binding_direction(
                        forms, cur, neighbor)
                    binding ? _add_met(parent_mono, met) :
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
        met = adj[minmax(a, b)]
        g1, g2 = form_to_group[a], form_to_group[b]
        g1 == g2 && continue  # intra-group

        # Forward: a→b
        fwd_met = (met !== nothing &&
            _is_binding_direction(forms, a, b)) ? met : nothing
        fwd_mono = fwd_met !== nothing ?
            _add_met(alpha_conc[a], fwd_met) : copy(alpha_conc[a])
        push!(get!(R_conc, (g1, g2), Set{MONO}()), fwd_mono)

        # Reverse: b→a
        rev_met = (met !== nothing &&
            _is_binding_direction(forms, b, a)) ? met : nothing
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

Compute the constraint descriptor: the set of (metabolite, mode)
pairs for each constrained equivalence group.
"""
function _constraint_descriptor(edges, adj, eq_steps, valid_groups,
                                constraint_mask)
    descriptor = Set{Tuple{Symbol, Symbol}}()
    for (gi, g) in enumerate(valid_groups)
        (constraint_mask >> (gi - 1)) & 1 == 1 || continue
        met = adj[minmax(edges[g[1]]...)]
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
        adj[minmax(a, b)] === nothing && return i
    end
    return 1
end

# ─── Dead-End SS/RE Propagation ───────────────────────────────

"""
    _dead_end_catalytic_map(edges, n_cat, site_defs, forms)

For each dead-end edge (index > n_cat), return the index of the
catalytic edge connecting the same forms with regulator sites
stripped, or `nothing` for regulator-binding edges (always RE).
"""
function _dead_end_catalytic_map(edges, n_cat, site_defs, forms)
    _strip_regs(fi) = Tuple(
        site_defs[k].role == :reg ? nothing : forms[fi].occupancy[k]
        for k in eachindex(site_defs))
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

# ─── RE/SS + Constraint Variants ──────────────────────────────

"""
    _find_equivalent_groups(edges, adj, site_defs, forms,
                            n_catalytic_edges, [de_cat_map])

Find groups of binding edges that involve the same non-product
metabolite. Catalytic edges + dead-end edges with catalytic
counterparts are considered.
"""
function _find_equivalent_groups(edges, adj, site_defs, forms,
    n_catalytic_edges, de_cat_map=nothing)
    groups = Dict{Symbol, Vector{Int}}()
    product_metabolites = Set(
        sd.metabolite for sd in site_defs if sd.role == :prod)
    for i in 1:n_catalytic_edges
        (a, b) = edges[i]
        met = get(adj, minmax(a, b), missing)
        !ismissing(met) && met !== nothing &&
            met ∉ product_metabolites &&
            push!(get!(groups, met, Int[]), i)
    end
    if de_cat_map !== nothing
        for (di, cat_idx) in enumerate(de_cat_map)
            cat_idx === nothing && continue
            edge_idx = n_catalytic_edges + di
            met = adj[minmax(edges[edge_idx]...)]
            met === nothing && continue
            met in product_metabolites && continue
            push!(get!(groups, met, Int[]), edge_idx)
        end
    end
    sort!([sort(v) for v in values(groups) if length(v) >= 2];
        by=first)
end

"""
Build parameter constraints from valid_groups and constraint_mask.
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
    _ress_variants(spec, adj, site_defs, forms; max_re_groups=7)

Generate RE/SS + parameter constraint variants for a mechanism,
deduplicated by (concentration fingerprint, constraint descriptor).
Returns a collected `Vector{MechanismSpec}`.
"""
function _ress_variants(spec, adj, site_defs, forms;
                        max_re_groups::Int=7)
    edges = spec.edges
    n = length(edges)
    n == 0 && return MechanismSpec[]
    n_cat = spec.n_catalytic_edges
    iso_idx = _find_first_isomerization(edges, adj)
    de_cat_map = _dead_end_catalytic_map(
        edges, n_cat, site_defs, forms)
    equiv_groups = _find_equivalent_groups(edges, adj, site_defs,
        forms, n_cat, de_cat_map)
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
        partition = _compute_re_partition(edges, eq_steps)
        G = length(partition)
        (G < 2 || G > max_re_groups) && continue
        fp = _concentration_fingerprint(
            edges, eq_steps, site_defs, forms, adj, partition)
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
    collect(values(best))
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

When `catalytic_n > 0`, also produces `OligomericEnzymeMechanism`
candidates (NConf=2) for each RE/SS variant.

# Keywords
- `max_re_groups::Int=7`: maximum RE groups (G ≤ max_re_groups).
- `catalytic_n::Int=0`: when > 0, also enumerate
  `OligomericEnzymeMechanism` candidates with NConf=2.
"""
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    stage::EnumerationStage=FullEnumeration(),
    max_re_groups::Int=7,
    catalytic_n::Int=0,
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)
    form_lookup = Dict(
        ntuple(k -> f.occupancy[k], length(site_defs)) => i
        for (i, f) in enumerate(forms))
    max_forms = length(forms)

    catalytic = _catalytic_topologies(
        site_defs, forms, adj, reaction; max_forms)
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
            catalytic, site_defs, forms, adj, form_lookup;
            max_forms, dead_end_regs, allosteric_regs)
        append!(all_inhibitor_specs, partition_specs)
    end
    stage isa WithDeadEnd && return all_inhibitor_specs

    # Collect all RE/SS variants eagerly (each _ress_variants
    # already materializes a Dict for dedup, so no extra cost)
    all_ress = [_ress_variants(s, adj, site_defs, forms;
                    max_re_groups)
                for s in all_inhibitor_specs]
    em_total = sum(length, all_ress)

    if catalytic_n > 0
        oem_total = sum(
            zip(all_inhibitor_specs, all_ress)) do (s, ress)
            _oligomeric_count(s, length(ress), catalytic_n)
        end
        total = em_total + oem_total
        inner = Iterators.flatmap(all_ress) do ress
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
        inner = Iterators.flatten(all_ress)
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
    site_defs, forms = enumerate_enzyme_forms(rxn)
    adj = _build_adjacency(site_defs, forms)
    used = Set(Iterators.flatten(spec.edges))
    _form_atoms(occ) = sort!([a => c for (a, c) in reduce(
        mergewith(+),
        (Dict(a) for a in occ if a !== nothing);
        init=Dict{Symbol,Int}())]; by=first)
    species = (
        substrates(rxn),
        products(rxn),
        regulators(rxn),
        Tuple((forms[i].name,
               Tuple(Tuple.(_form_atoms(forms[i].occupancy))))
            for i in eachindex(forms) if i ∈ used))
    reactions = Tuple(let met = adj[minmax(a, b)]
        if met === nothing  # isomerization
            ((forms[a].name,), (forms[b].name,))
        elseif _is_binding_direction(forms, a, b)
            ((forms[a].name, met), (forms[b].name,))
        else
            ((forms[a].name,), (forms[b].name, met))
        end
    end for (a, b) in spec.edges)
    EnzymeMechanism(species, reactions,
        Tuple(spec.equilibrium_steps),
        Tuple((t, c, Tuple(Tuple.(f)))
              for (t, c, f) in spec.param_constraints))
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
"""
function _oligomeric_count(spec, ress_count, catalytic_n)
    k = length(spec.allosteric_regulators)
    ress_count * _partition_mult_count(k, catalytic_n)
end

"""
    _expand_oligomeric_variants(em_spec, catalytic_n)

Generate OligomericEnzymeMechanism specs by enumerating set
partitions of allosteric regulators into site groups, then
multiplicity combos (each ∈ 1:catalytic_n) per site group.
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
