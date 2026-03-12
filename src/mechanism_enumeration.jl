# ─── New Mechanism Enumeration (Staged Pipeline) ──────────────

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

# ─── Mechanism Spec Types ──────────────────────────────────────

abstract type AbstractMechanismSpec end

"""Elementary step in canonical binding direction (metabolite on LHS)."""
struct StepSpec
    reactants::Vector{Symbol}   # [:E, :S] or [:EAB]
    products::Vector{Symbol}    # [:ES] or [:EPQ]
    is_equilibrium::Bool
end

"""
    MechanismSpec <: AbstractMechanismSpec

Represents a monomeric enzyme mechanism specification in the
staged enumeration pipeline. Steps are `StepSpec` values with
inline form names and equilibrium status.
"""
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{StepSpec}
    param_constraints::Vector{ParamConstraint}
    param_count::Int
end

"""
    AllostericMechanismSpec <: AbstractMechanismSpec

Represents an allosteric enzyme mechanism built from a
base `MechanismSpec` plus allosteric site and multiplicity info.
`tr_equiv_metabolites` lists metabolites with K_T = K_R.
"""
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    tr_equiv_metabolites::Vector{Symbol}
end

# ─── StepSpec Helpers ──────────────────────────────────────────

"""Return the metabolite for a step, or nothing for isomerization."""
step_metabolite(s::StepSpec) =
    length(s.reactants) == 2 ? s.reactants[2] : nothing

"""Return (from_form, to_form) for a step."""
step_forms(s::StepSpec) = (s.reactants[1], s.products[1])

"""Collect all unique form names from steps."""
function all_form_names(spec::MechanismSpec)
    forms = Set{Symbol}()
    for s in spec.steps
        push!(forms, s.reactants[1])
        push!(forms, s.products[1])
    end
    forms
end

function all_form_names(steps::Vector{StepSpec})
    forms = Set{Symbol}()
    for s in steps
        push!(forms, s.reactants[1])
        push!(forms, s.products[1])
    end
    forms
end

# ─── Enumeration Stage Types (used by enumerate_mechanisms) ────

abstract type EnumerationStage end
struct Catalytic       <: EnumerationStage end
struct WithDeadEnd     <: EnumerationStage end
struct FullEnumeration <: EnumerationStage end

struct MechanismIterator
    inner::Any
    total::Int
end

Base.eltype(::Type{MechanismIterator}) = AbstractMechanismSpec
Base.IteratorSize(::Type{MechanismIterator}) = Base.HasLength()
Base.length(iter::MechanismIterator) = iter.total
Base.iterate(iter::MechanismIterator, s...) = iterate(iter.inner, s...)

# ─── Edge Classification + Adjacency ─────────────────────────

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
            sd.role in (:sub, :prod) && occ_b != sd.full_atoms &&
                return missing
            return sd.metabolite
        elseif occ_b === nothing
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
    nothing  # isomerization
end

"""
    _build_adjacency(site_defs, forms)

Build undirected adjacency from all valid elementary steps.
Keys are `(i, j)` with `i < j`.
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

Return true if traversing from→to binds a metabolite.
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
Returns `(site_defs, forms)`.
"""
function enumerate_enzyme_forms(
    reaction::EnzymeReaction{S,P,R},
) where {S,P,R}
    sorted_atoms(spec) = sort([a => c for (a, c) in spec[2]]; by=first)
    atoms_label(v) = join(
        string(s) * (c > 1 ? string(c) : "") for (s, c) in v)
    product_atom_dicts = [Dict{Symbol,Int}(a => c for (a, c) in p[2])
                          for p in P if !isempty(p[2])]
    site_defs = SiteDefinition[]
    per_site_options =
        Vector{Tuple{Union{Nothing,Vector{Pair{Symbol,Int}}},String}}[]
    for (group, role) in ((S, :sub), (P, :prod)), spec in group
        fa = sorted_atoms(spec)
        push!(site_defs, SiteDefinition(spec[1], role, fa))
        site_options = [(nothing, "0"), (fa, string(spec[1]))]
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
    # Regulators: (name, role) pairs — extract just the name
    for reg in R
        rname = reg isa Symbol ? reg : reg[1]
        push!(site_defs, SiteDefinition(rname, :reg, Pair{Symbol,Int}[]))
        push!(per_site_options,
            [(nothing, "0"), (Pair{Symbol,Int}[], string(rname))])
    end
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

# ─── Dead-End Helpers ────────────────────────────────────────

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
counterpart.
"""
function _propagate_de_eq_steps!(eq_steps, n_cat, de_cat_map)
    for (di, cat_idx) in enumerate(de_cat_map)
        cat_idx === nothing && continue
        eq_steps[n_cat + di] = eq_steps[cat_idx]
    end
end

# ─── RE Partition (union-find) ─────────────────────────────────

"""
    _compute_re_partition(edges, eq_steps) -> Vector{Vector{Int}}

Connected components when only RE edges are considered.
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

Concentration monomials in the rate equation denominator.
"""
function _concentration_fingerprint(
    edges, eq_steps, site_defs, forms, adj, partition,
)
    G = length(partition)
    group_set = [Set(g) for g in partition]
    form_to_group = Dict(
        i => g for (g, grp) in enumerate(partition) for i in grp)

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

    sigma_conc = [Set{MONO}(alpha_conc[i] for i in group)
                  for group in partition]

    R_conc = Dict{Tuple{Int,Int}, Set{MONO}}()
    for (idx, (a, b)) in enumerate(edges)
        eq_steps[idx] && continue
        met = adj[minmax(a, b)]
        g1, g2 = form_to_group[a], form_to_group[b]
        g1 == g2 && continue

        fwd_met = (met !== nothing &&
            _is_binding_direction(forms, a, b)) ? met : nothing
        fwd_mono = fwd_met !== nothing ?
            _add_met(alpha_conc[a], fwd_met) : copy(alpha_conc[a])
        push!(get!(R_conc, (g1, g2), Set{MONO}()), fwd_mono)

        rev_met = (met !== nothing &&
            _is_binding_direction(forms, b, a)) ? met : nothing
        rev_mono = rev_met !== nothing ?
            _add_met(alpha_conc[b], rev_met) : copy(alpha_conc[b])
        push!(get!(R_conc, (g2, g1), Set{MONO}()), rev_mono)
    end

    fingerprint = Set{MONO}()
    for g in 1:G
        D_g = _spanning_arborescence_monomials(G, R_conc, g)
        for s in sigma_conc[g], d in D_g
            push!(fingerprint, _mono_mul(s, d))
        end
    end
    fingerprint
end

# ─── Equivalence Groups + Constraints ──────────────────────────

"""
    _constraint_descriptor(edges, adj, eq_steps, valid_groups,
                           constraint_mask)

Constraint descriptor: set of (metabolite, mode) pairs for
each constrained equivalence group.
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

Find index of the first isomerization edge. Falls back to 1.
"""
function _find_first_isomerization(edges, adj)
    for (i, (a, b)) in enumerate(edges)
        adj[minmax(a, b)] === nothing && return i
    end
    return 1
end

"""
    _find_equivalent_groups(edges, adj, site_defs, forms,
                            n_catalytic_edges, [de_cat_map])

Find groups of binding edges involving the same non-product
metabolite AND same site type.
"""
function _find_equivalent_groups(edges, adj, site_defs, forms,
    n_catalytic_edges, de_cat_map=nothing)
    groups = Dict{Tuple{Symbol,Symbol}, Vector{Int}}()
    product_metabolites = Set(
        sd.metabolite for sd in site_defs if sd.role == :prod)
    for i in 1:n_catalytic_edges
        (a, b) = edges[i]
        met = get(adj, minmax(a, b), missing)
        !ismissing(met) && met !== nothing &&
            met ∉ product_metabolites || continue
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
    :iso
end

"""Build parameter constraints from valid_groups and constraint_mask."""
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

# ─── Set Partitions + Allosteric Helpers ─────────────────────

"""
    _set_partitions(elements::Vector{Symbol})

Enumerate all set partitions (Bell number partitions).
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

Count allosteric multiplicity variants: sum_{g=1}^{k} S(k,g) * N^g.
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

# ─── Stage 1: Catalytic Topologies ───────────────────────────

"""
    _catalytic_topologies(reaction) -> Vector{MechanismSpec}

Enumerate catalytic cycle topologies with initial RE/SS assignment
and param_count.
"""
function _catalytic_topologies(
    @nospecialize(reaction::EnzymeReaction),
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)
    max_forms = length(forms)

    sub_set = Set(s[1] for s in substrates(reaction))
    prod_set = Set(p[1] for p in products(reaction))
    free = findfirst(
        f -> all(o === nothing for o in f.occupancy), forms)
    free === nothing && return MechanismSpec[]
    cat_forms = Set(i for (i, f) in enumerate(forms)
        if all(k -> site_defs[k].role != :reg ||
                    f.occupancy[k] === nothing,
               eachindex(site_defs)))
    cat_adj = Dict(
        (a, b) => met for ((a, b), met) in adj
        if a ∈ cat_forms && b ∈ cat_forms)
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
            binding = met !== nothing &&
                _is_binding_direction(forms, cur, neighbor)
            releasing = met !== nothing && !binding
            if met === nothing
                neighbor ∉ visited || continue
            elseif binding
                met ∈ sub_set && met ∉ bound && !has_res ||
                    continue
                push!(bound, met)
            else
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

    result = MechanismSpec[]
    for form_set in combined
        edges = [(a, b) for ((a, b), _) in adj
                 if a ∈ form_set && b ∈ form_set]
        n_edges = length(edges)
        iso_idx = _find_first_isomerization(edges, adj)
        eq_steps = fill(true, n_edges)
        eq_steps[iso_idx] = false

        n_forms = length(form_set)
        n_independent_cycles = n_edges - n_forms + 1
        n_thermo = n_independent_cycles
        n_re = n_edges - 1  # all except the one SS
        n_ss = 1
        param_count = n_re + 2 * n_ss - n_thermo + 2

        push!(result, MechanismSpec(reaction, edges, n_edges,
            eq_steps, ParamConstraint[], param_count))
    end
    result
end

# ─── Stage 2: RE/SS Assignment ───────────────────────────────

"""
    _expand_ress_variants(specs, reaction; max_re_groups=7)
        -> Vector{MechanismSpec}

Enumerate RE/SS assignment combinations for catalytic edges.
The first isomerization edge is always SS. Dead-end edges
inherit RE/SS from catalytic counterparts.
"""
function _expand_ress_variants(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    max_re_groups::Int=7,
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)

    result = MechanismSpec[]
    for spec in specs
        edges = spec.edges
        n = length(edges)
        n == 0 && continue
        n_cat = spec.n_catalytic_edges
        iso_idx = _find_first_isomerization(edges, adj)
        de_cat_map = _dead_end_catalytic_map(
            edges, n_cat, site_defs, forms)
        other_indices = [i for i in 1:n_cat if i != iso_idx]
        n_other = length(other_indices)

        # Sort masks by popcount so fewer SS steps come first
        masks = sort!(collect(0:(1 << n_other) - 1); by=count_ones)

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

            # Compute param_count
            n_re = count(eq_steps)
            n_ss = n - n_re
            n_forms = length(Set(Iterators.flatten(edges)))
            n_independent_cycles = n - n_forms + 1
            n_thermo = n_independent_cycles
            param_count = n_re + 2 * n_ss - n_thermo + 2

            push!(result, MechanismSpec(reaction, edges, n_cat,
                eq_steps, ParamConstraint[], param_count))
        end
    end
    result
end

# ─── Stage 3: Dead-End Inhibitor Expansion ────────────────────

"""
    _expand_dead_end_inhibitors(specs, reaction;
        dead_end_regs) -> Vector{MechanismSpec}

Expand mechanism topologies with dead-end inhibitor configurations.
Only regulators in `dead_end_regs` create dead-end complexes.
"""
function _expand_dead_end_inhibitors(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    dead_end_regs::Vector{Symbol}=Symbol[],
)
    site_defs, forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(site_defs, forms)
    form_lookup = Dict(
        ntuple(k -> f.occupancy[k], length(site_defs)) => i
        for (i, f) in enumerate(forms))
    max_forms = length(forms)

    isempty(dead_end_regs) && return specs

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
        dead_end_lookup = Dict(
            (ti, mask) => fi2
            for (ti, fi) in enumerate(topo_sorted)
            for mask in 1:(1 << n_inh) - 1
            for fi2 in (_find_dead_end(form_lookup, site_defs,
                forms[fi],
                Tuple(reg_positions[j] for j in 1:n_inh
                      if (mask >> (j - 1)) & 1 == 1)),)
            if fi2 !== nothing && fi2 ∉ topo_nodes)
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
            n_new = length(all_edges)

            eq_steps = fill(true, n_new)
            for i in 1:length(spec.edges)
                eq_steps[i] = spec.equilibrium_steps[i]
            end
            n_cat = spec.n_catalytic_edges
            de_cat_map = _dead_end_catalytic_map(
                all_edges, n_cat, site_defs, forms)
            _propagate_de_eq_steps!(eq_steps, n_cat, de_cat_map)

            n_re = count(eq_steps)
            n_ss = n_new - n_re
            n_all_forms = length(all_forms)
            n_independent_cycles = n_new - n_all_forms + 1
            n_thermo = n_independent_cycles
            param_count = n_re + 2 * n_ss - n_thermo + 2

            push!(result, MechanismSpec(reaction, all_edges,
                n_cat, eq_steps, ParamConstraint[], param_count))
        end
    end
    result
end

# ─── Stage 4: Equivalence Constraint Expansion ────────────────

"""
    _expand_equivalence_constraints(specs, reaction)
        -> Vector{MechanismSpec}

For each spec, enumerate equivalence constraint masks.
"""
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
        valid_groups = [g for g in equiv_groups
            if all(spec.equilibrium_steps[s] ==
                   spec.equilibrium_steps[g[1]] for s in g)]

        for constraint_mask in 0:(1 << length(valid_groups)) - 1
            constraints = _build_constraints(
                valid_groups, spec.equilibrium_steps,
                constraint_mask, edges)
            # Compute param_count delta from constraints
            delta = 0
            for (gi, g) in enumerate(valid_groups)
                (constraint_mask >> (gi-1)) & 1 == 1 || continue
                n_constrained = length(g) - 1
                if spec.equilibrium_steps[g[1]]  # RE
                    delta -= n_constrained
                else  # SS
                    delta -= 2 * n_constrained
                end
            end
            push!(result, MechanismSpec(spec.reaction, edges,
                n_cat, spec.equilibrium_steps, constraints,
                spec.param_count + delta))
        end
    end
    result
end

# ─── Stage 5: Deduplication ────────────────────────────────────

"""
    _deduplicate(specs, reaction) -> Vector{MechanismSpec}

Deduplicate by (concentration fingerprint, constraint descriptor).
Keeps mechanism with fewest parameters.
"""
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

"""Reverse-map param_constraints to a bitmask over valid_groups."""
function _constraints_to_mask(constraints, valid_groups, eq_steps,
                              edges)
    mask = 0
    constrained_edge_indices = Set{Int}()
    for (target, _, srcs) in constraints
        # Parse edge index from target symbol name
        m = match(r"[kK](\d+)", string(target))
        m !== nothing && m[1] !== nothing &&
            push!(constrained_edge_indices, parse(Int, m[1]::SubString))
    end
    for (gi, g) in enumerate(valid_groups)
        if any(idx in constrained_edge_indices for idx in g[2:end])
            mask |= (1 << (gi - 1))
        end
    end
    mask
end

# ─── Allosteric Expansion ─────────────────────────────────────

"""
    _expand_allosteric(specs, reaction;
        catalytic_n, allosteric_regs)
        -> Vector{AllostericMechanismSpec}

Expand monomeric specs into allosteric (MWC) variants.
"""
function _expand_allosteric(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    catalytic_n::Int=2,
    allosteric_regs::Vector{Symbol}=Symbol[],
)
    result = AllostericMechanismSpec[]
    if isempty(allosteric_regs)
        # No allosteric regulators: catalytic subunits only
        for spec in specs
            push!(result, AllostericMechanismSpec(
                spec, catalytic_n,
                Vector{Symbol}[], Int[], Symbol[]))
        end
        return result
    end
    partitions = _set_partitions(allosteric_regs)
    for spec in specs
        for partition in partitions
            n_groups = length(partition)
            for combo in Iterators.product(
                    ntuple(_ -> 1:catalytic_n, n_groups)...)
                push!(result, AllostericMechanismSpec(
                    spec, catalytic_n, partition,
                    collect(combo), Symbol[]))
            end
        end
    end
    result
end

# ─── T/R Equivalence ─────────────────────────────────────────

"""
    _expand_tr_equivalence(specs, reaction)
        -> Vector{AllostericMechanismSpec}

Enumerate T/R parameter equivalence variants. For each metabolite
with a T-state parameter, K_T can equal K_R (fewer params) or be
independent. Produces 2^n variants per input spec where n is the
number of metabolites with T-state binding parameters.
"""
function _expand_tr_equivalence(
    specs::Vector{AllostericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = AllostericMechanismSpec[]
    site_defs, forms = EnzymeRates.enumerate_enzyme_forms(reaction)
    adj = EnzymeRates._build_adjacency(site_defs, forms)

    for spec in specs
        # Collect metabolites with T-state params
        t_mets = Symbol[]

        # 1. All metabolites in RE binding edges (catalytic + dead-end)
        for (i, (ei, ej)) in enumerate(spec.base.edges)
            spec.base.equilibrium_steps[i] || continue
            key = minmax(ei, ej)
            met = get(adj, key, nothing)
            met !== nothing && met ∉ t_mets && push!(t_mets, met)
        end

        # 2. Regulator ligands
        for site in spec.allosteric_reg_sites
            for lig in site
                lig ∉ t_mets && push!(t_mets, lig)
            end
        end

        n = length(t_mets)
        for mask in 0:(1 << n) - 1
            equiv = Symbol[t_mets[i] for i in 1:n
                          if ((mask >> (i-1)) & 1) == 1]
            push!(result, AllostericMechanismSpec(
                spec.base, spec.catalytic_n,
                spec.allosteric_reg_sites,
                spec.allosteric_multiplicities,
                equiv))
        end
    end
    result
end

# ─── Post-Allosteric Deduplication ────────────────────────────

"""
    _deduplicate_allosteric(specs, reaction) -> Vector{AllostericMechanismSpec}

Remove T<->R mirror duplicates.
"""
function _deduplicate_allosteric(
    specs::Vector{AllostericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    seen = Dict{Any, AllostericMechanismSpec}()
    for spec in specs
        key = _allosteric_canonical_key(spec)
        if !haskey(seen, key)
            seen[key] = spec
        end
    end
    collect(values(seen))
end

"""Canonical key for allosteric dedup: includes sorted TR equiv metabolites."""
function _allosteric_canonical_key(spec::AllostericMechanismSpec)
    base_key = (spec.base.edges, spec.base.equilibrium_steps,
                spec.base.param_constraints)
    # Sort reg sites and multiplicities together
    pairs = collect(zip(spec.allosteric_reg_sites,
                        spec.allosteric_multiplicities))
    sort!(pairs)
    sorted_sites = [p[1] for p in pairs]
    sorted_mults = [p[2] for p in pairs]
    (base_key, spec.catalytic_n, sorted_sites, sorted_mults,
     sort(spec.tr_equiv_metabolites))
end

# ─── MechanismSpec → EnzymeMechanism Conversion ──────────────

"""
    compile_mechanism(spec::MechanismSpec)

Convert a `MechanismSpec` to an `EnzymeMechanism`.
"""
function compile_mechanism(spec::MechanismSpec)
    _compile_enzyme_mechanism(spec)
end

"""Construct EnzymeMechanism from MechanismSpec (backward compat)."""
EnzymeMechanism(spec::MechanismSpec) = _compile_enzyme_mechanism(spec)

function _compile_enzyme_mechanism(spec::MechanismSpec)
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
        if met === nothing
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

function compile_mechanism(spec::AllostericMechanismSpec)
    cm = compile_mechanism(spec.base)
    cat_mets = metabolites(cm)

    # Build Metabolites tuple (catalytic + regulatory)
    reg_syms = Symbol[]
    for site in spec.allosteric_reg_sites
        for s in site
            s in reg_syms || s in cat_mets || push!(reg_syms, s)
        end
    end
    mets = (cat_mets..., reg_syms...)

    # Build CatSites: (catalytic_metabolites, multiplicity, tr_equiv_mets)
    cat_tr = Tuple(m for m in cat_mets
                   if m in spec.tr_equiv_metabolites)
    cat_sites = (cat_mets, spec.catalytic_n, cat_tr)

    # Build RegSites with TR equivalence info
    reg_sites = Tuple(
        (Tuple(group), mult,
         Tuple(lig for lig in group
               if lig in spec.tr_equiv_metabolites))
        for (group, mult) in zip(
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities))

    AllostericEnzymeMechanism{mets, typeof(cm), cat_sites, reg_sites}()
end

# ─── Pipeline Orchestration ──────────────────────────────────

"""
    enumerate_mechanisms(reaction; max_re_groups=7, catalytic_n=0)

Enumerate valid mechanism topologies for the given reaction
using a staged pipeline.
"""
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    stage::EnumerationStage=FullEnumeration(),
    max_re_groups::Int=7,
    catalytic_n::Int=0,
)
    # Stage 1: Catalytic topologies
    catalytic = _catalytic_topologies(reaction)
    stage isa Catalytic && return catalytic

    # Regulator partitioning: fixed-role regs + 2^n partitions
    # of unknown-role regs
    roles = regulator_roles(reaction)
    fixed_dead_end = Symbol[r[1] for r in roles if r[2] == :dead_end]
    fixed_allosteric = Symbol[r[1] for r in roles
                              if r[2] == :allosteric]
    unknown = Symbol[r[1] for r in roles if r[2] == :unknown]
    n_unknown = length(unknown)

    # WithDeadEnd stage: run dead-end on catalytic topologies
    # directly (for backward compat with tests that need large
    # mechanisms without running the full pipeline)
    if stage isa WithDeadEnd
        all_de = MechanismSpec[]
        for reg_mask in 0:(1 << n_unknown) - 1
            de_regs = Symbol[fixed_dead_end;
                [unknown[i] for i in 1:n_unknown
                 if (reg_mask >> (i - 1)) & 1 == 0]]
            append!(all_de, _expand_dead_end_inhibitors(
                catalytic, reaction; dead_end_regs=de_regs))
        end
        return all_de
    end

    all_base = MechanismSpec[]
    all_allosteric = AllostericMechanismSpec[]

    for reg_mask in 0:(1 << n_unknown) - 1
        de_regs = Symbol[fixed_dead_end;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 0]]
        allo_regs = Symbol[fixed_allosteric;
            [unknown[i] for i in 1:n_unknown
             if (reg_mask >> (i - 1)) & 1 == 1]]

        # Phase 1: base mechanism pipeline (chained)
        base = _expand_ress_variants(
            catalytic, reaction; max_re_groups)
        base = _expand_dead_end_inhibitors(
            base, reaction; dead_end_regs=de_regs)
        base = _expand_equivalence_constraints(base, reaction)
        base = _deduplicate(base, reaction)
        append!(all_base, base)

        # Phase 2: allosteric expansion (independent)
        if !isempty(allo_regs)
            cn = catalytic_n > 0 ? catalytic_n : 1
            allo = _expand_allosteric(base, reaction;
                catalytic_n=cn, allosteric_regs=allo_regs)
            allo = _expand_tr_equivalence(allo, reaction)
            allo = _deduplicate_allosteric(allo, reaction)
            append!(all_allosteric, allo)
        end
    end

    total = length(all_base) + length(all_allosteric)
    inner = Iterators.flatten((all_base, all_allosteric))
    MechanismIterator(inner, total)
end
