# ─── Enzyme Form Enumeration ─────────────────────────────────────────────────
#
# Given an EnzymeReaction with per-metabolite max binding sites,
# enumerate all possible enzyme forms via unified Cartesian product
# over per-site options (standard + ping-pong residuals).

"""
    SiteState

State of a single binding site on the enzyme.

- `metabolite`: which metabolite this site is for
- `index`: site number (1, 2, ...)
- `atoms`: atoms present (`nothing` = unoccupied, empty vector = occupied with no-atom species)
"""
struct SiteState
    metabolite::Symbol
    index::Int
    atoms::Union{Nothing, Vector{Pair{Symbol,Int}}}
end

"""
    EnzymeFormSpec

Specification of an enzyme form with named binding sites.

- `name`: canonical name encoding the site state vector (e.g., `:E_S_0`)
- `sites`: ordered vector of `SiteState`, one per binding site
"""
struct EnzymeFormSpec
    name::Symbol
    sites::Vector{SiteState}
end

"""
    enumerate_enzyme_forms(reaction::EnzymeReaction)

Enumerate all possible enzyme forms for the given reaction.

Each binding site is distinguishable. Forms include standard forms (each site
independently empty or fully occupied) and ping-pong intermediates (sites with
partial residual atom content after product release).

**Exclusion rule**: Forms where all substrate sites are fully bound while any
product site is occupied (or vice versa) are excluded as physically impossible.
Ping-pong residuals count as "occupied" but not "fully bound."

Site order: core substrates, core products, extra substrate sites, extra product
sites, regulator sites.

Returns a `Vector{EnzymeFormSpec}`.
"""
function enumerate_enzyme_forms(reaction::EnzymeReaction{S,P,R}) where {S,P,R}
    fatoms(spec) = sort([a => c for (a, c) in spec[2]]; by=first)
    nsites(spec) = length(spec) >= 3 ? spec[3] : 1
    astr(v) = join(string(s) * (c > 1 ? string(c) : "") for (s, c) in v)

    # 1. Ping-pong residuals: partial atom contents remaining after product release
    prod_atoms_list = [Dict{Symbol,Int}(a => c for (a, c) in p[2])
                       for p in P if !isempty(p[2])]
    residuals = Dict{Symbol, Vector{Vector{Pair{Symbol,Int}}}}()
    for sub_spec in S
        isempty(sub_spec[2]) && continue
        isempty(prod_atoms_list) && continue
        sub_atoms = Dict{Symbol,Int}(a => c for (a, c) in sub_spec[2])
        sub_residuals = Vector{Pair{Symbol,Int}}[]
        for mask in 1:(2^length(prod_atoms_list) - 1)
            combined = reduce(mergewith(+),
                (prod_atoms_list[i] for i in 1:length(prod_atoms_list) if (mask >> (i-1)) & 1 == 1))
            all(get(sub_atoms, atom, 0) >= count for (atom, count) in combined) || continue
            residual = Dict{Symbol,Int}(atom => count - get(combined, atom, 0)
                                         for (atom, count) in sub_atoms
                                         if count > get(combined, atom, 0))
            if !isempty(residual) && residual != sub_atoms
                r = sort([k => v for (k, v) in residual]; by=first)
                r ∉ sub_residuals && push!(sub_residuals, r)
            end
        end
        !isempty(sub_residuals) && (residuals[sub_spec[1]] = sub_residuals)
    end

    # 2. Build per-site data with pre-computed (content, label) options
    mets = Symbol[]
    idxs = Int[]
    fulls = Vector{Pair{Symbol,Int}}[]
    roles = Symbol[]  # :sub, :prod, :reg
    opts = Vector{Tuple{Union{Nothing, Vector{Pair{Symbol,Int}}}, String}}[]
    function add!(met, idx, full, role)
        push!(mets, met)
        push!(idxs, idx)
        push!(fulls, full)
        push!(roles, role)
        site_opts = Tuple{Union{Nothing, Vector{Pair{Symbol,Int}}}, String}[
            (nothing, "0"), (full, string(met))]
        for r in get(residuals, met, Vector{Pair{Symbol,Int}}[])
            push!(site_opts, (r, astr(r)))
        end
        push!(opts, site_opts)
    end
    # Site order: core substrates, core products, extra substrate sites, extra product sites, regulators
    for s in S;  add!(s[1], 1, fatoms(s), :sub);  end
    for p in P;  add!(p[1], 1, fatoms(p), :prod); end
    for s in S, i in 2:nsites(s);  add!(s[1], i, fatoms(s), :sub);  end
    for p in P, i in 2:nsites(p);  add!(p[1], i, fatoms(p), :prod); end
    for r in R, i in 1:nsites(r);  add!(r[1], i, fatoms(r), :reg);  end

    # 3. Enumerate Cartesian product with exclusion filter
    #    Uses flat modular-arithmetic indexing to avoid Iterators.product splat,
    #    which would create 2^n specialized tuple types and explode compilation.
    n = length(mets)
    forms = EnzymeFormSpec[]
    total = prod(length(o) for o in opts)
    name_parts = Vector{String}(undef, n)
    sites = Vector{SiteState}(undef, n)
    for idx in 0:total-1
        rem = idx
        all_sub_full = true
        all_prod_full = true
        any_sub_occ = false
        any_prod_occ = false
        for i in n:-1:1
            content, label = opts[i][rem % length(opts[i]) + 1]
            rem ÷= length(opts[i])
            name_parts[i] = label
            sites[i] = SiteState(mets[i], idxs[i],
                                 content === nothing ? nothing : copy(content))
            if roles[i] == :sub
                content != fulls[i] && (all_sub_full = false)
                content !== nothing && (any_sub_occ = true)
            elseif roles[i] == :prod
                content != fulls[i] && (all_prod_full = false)
                content !== nothing && (any_prod_occ = true)
            end
        end
        (all_sub_full && any_prod_occ) && continue
        (all_prod_full && any_sub_occ) && continue
        push!(forms, EnzymeFormSpec(Symbol("E_" * join(name_parts, "_")), copy(sites)))
    end

    return forms
end

# ─── Mechanism Enumeration ───────────────────────────────────────────────────
#
# Given an EnzymeReaction, enumerate all valid mechanism topologies
# (catalytic cycles + dead-ends + RE/SS assignments + equivalent step constraints).

"""
    ReactionEdge

A directed elementary reaction between two enzyme forms.
"""
struct ReactionEdge
    from::Int                           # index into forms vector
    to::Int                             # index into forms vector
    metabolite::Union{Nothing, Symbol}  # nothing for isomerization
    edge_type::Symbol                   # :binding, :release, :isomerization
end

"""
    Topology

A valid catalytic topology: a connected subgraph of the reaction graph
containing the free enzyme, with all cycles being 1× or 0× stoichiometry.
"""
struct Topology
    form_indices::Vector{Int}
    edges::Vector{ReactionEdge}
end

"""
    MechanismSpec

Lightweight runtime description of a mechanism (not type-parameterized).
Convert to `EnzymeMechanism` via `EnzymeMechanism(spec)`.
"""
struct MechanismSpec
    reaction::Any
    forms::Vector{Symbol}
    form_atoms::Vector{Vector{Pair{Symbol,Int}}}
    reactions::Vector{Tuple{Vector{Symbol}, Vector{Symbol}}}
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{Tuple{Symbol, Int, Vector{Tuple{Symbol,Int}}}}
end

"""
    DeadEndTree

Dead-end forms reachable from a cycle form, organized as a tree for
downward-closed subset enumeration.
"""
struct DeadEndTree
    form_idx::Int                # index into all_forms
    children::Vector{DeadEndTree}
end

"""
    MechanismIterator

Lazy iterator over all valid mechanisms for a reaction. Eagerly enumerates
topologies; lazily iterates over (dead-end × RE/SS × constraint) combinations.
"""
struct MechanismIterator
    all_forms::Vector{EnzymeFormSpec}
    reaction_graph::Vector{ReactionEdge}
    topologies::Vector{Topology}
    reaction::Any
    max_forms::Int
    # Per-topology precomputed data
    dead_end_trees::Vector{Vector{Vector{DeadEndTree}}}  # [topo][cycle_form] -> trees
    adjacency::Vector{Vector{Int}}  # outgoing edges per form (indices into reaction_graph)
end

Base.eltype(::Type{MechanismIterator}) = MechanismSpec
Base.IteratorSize(::Type{MechanismIterator}) = Base.HasLength()

"""
    n_sites(reaction::EnzymeReaction)

Total number of binding sites across all substrates, products, and regulators.
"""
function n_sites(@nospecialize(reaction::EnzymeReaction))
    count = 0
    for spec in Iterators.flatten((substrates(reaction), products(reaction), regulators(reaction)))
        count += length(spec) >= 3 ? spec[3] : 1
    end
    count
end

# ─── Shared Helpers ──────────────────────────────────────────────────────────

"""Number of catalytic (substrate + product) site positions in a form's sites vector.
Sites at positions > n_catalytic are regulator sites."""
function _n_catalytic_sites(@nospecialize(reaction::EnzymeReaction))
    sum(length(s) >= 3 ? s[3] : 1 for s in substrates(reaction)) +
    sum(length(p) >= 3 ? p[3] : 1 for p in products(reaction))
end

"""Undirected edge key for deduplication: canonical (min, max, metabolite) triple."""
_undirected_key(e::ReactionEdge) = (min(e.from, e.to), max(e.from, e.to), e.metabolite)

"""Build expected net stoichiometry from a reaction: -1 per substrate, +1 per product, 0 per regulator."""
function _expected_stoichiometry(@nospecialize(reaction::EnzymeReaction))
    stoich = Dict{Symbol,Int}()
    for s in substrates(reaction); stoich[s[1]] = get(stoich, s[1], 0) - 1; end
    for p in products(reaction); stoich[p[1]] = get(stoich, p[1], 0) + 1; end
    for r in regulators(reaction); stoich[r[1]] = get(stoich, r[1], 0); end
    stoich
end

"""Find the index of the free enzyme form (all sites empty)."""
function _free_enzyme_index(forms::Vector{EnzymeFormSpec})
    findfirst(f -> all(s.atoms === nothing for s in f.sites), forms)
end

"""
    _max_cycle_forms(forms, reaction)

Upper bound on the number of forms in any valid 1× catalytic cycle.

Formula: `n_only_sub_patterns + n_only_prod_patterns + 1 + 2 * n_regulators`

- `n_only_sub_patterns`: distinct catalytic-site fingerprints with only substrate positions occupied
- `n_only_prod_patterns`: same for product positions
- `+1`: free enzyme
- `+2 * n_reg`: each regulator could be an essential activator (bind + unbind = 2 extra forms)

The bound is safe because every form in a valid 1× cycle must have either only substrates
or only products at catalytic positions (mixed sub+prod violates stoichiometry).
"""
function _max_cycle_forms(forms::Vector{EnzymeFormSpec}, @nospecialize(reaction::EnzymeReaction))
    n_cat = _n_catalytic_sites(reaction)
    n_reg = length(regulators(reaction))
    sub_names = Set(s[1] for s in substrates(reaction))
    is_sub = [forms[1].sites[k].metabolite in sub_names for k in 1:n_cat]

    only_sub = Set{Vector{Union{Nothing, Vector{Pair{Symbol,Int}}}}}()
    only_prod = Set{Vector{Union{Nothing, Vector{Pair{Symbol,Int}}}}}()

    for form in forms
        fp = Vector{Union{Nothing, Vector{Pair{Symbol,Int}}}}([form.sites[k].atoms for k in 1:n_cat])
        has_sub = any(fp[k] !== nothing for k in 1:n_cat if is_sub[k])
        has_prod = any(fp[k] !== nothing for k in 1:n_cat if !is_sub[k])
        if has_sub && !has_prod
            push!(only_sub, fp)
        elseif has_prod && !has_sub
            push!(only_prod, fp)
        end
    end
    length(only_sub) + length(only_prod) + 1 + 2 * n_reg
end

# ─── Reaction Graph Construction ─────────────────────────────────────────────

"""Build lookup from sorted atom vectors to metabolite symbols."""
function _build_met_atoms_lookup(@nospecialize(reaction::EnzymeReaction))
    lookup = Dict{Vector{Pair{Symbol,Int}}, Symbol}()
    for spec in Iterators.flatten((substrates(reaction), products(reaction), regulators(reaction)))
        atoms = sort([a => c for (a, c) in spec[2]]; by=first)
        !isempty(atoms) && (lookup[atoms] = spec[1])
    end
    lookup
end

"""
Compute the atom content of a form at its core catalytic sites.
Core sites = index-1 sites at positions 1:n_catalytic (substrates and products only).
"""
function _core_atoms(form::EnzymeFormSpec, n_catalytic::Int)
    atoms = Dict{Symbol,Int}()
    for k in 1:n_catalytic
        site = form.sites[k]
        site.index != 1 && continue
        site.atoms === nothing && continue
        for (a, c) in site.atoms
            atoms[a] = get(atoms, a, 0) + c
        end
    end
    atoms
end

"""
Build the directed reaction graph from enzyme forms.

Returns a Vector{ReactionEdge} containing all valid elementary reactions.
"""
function _build_reaction_graph(forms::Vector{EnzymeFormSpec}, @nospecialize(reaction::EnzymeReaction))
    met_lookup = _build_met_atoms_lookup(reaction)
    n = length(forms)
    nsites = length(forms[1].sites)
    edges = ReactionEdge[]

    # Full atom content per metabolite — used to distinguish standard binding from
    # ping-pong residuals.  Only standard binding (empty ↔ fully occupied) is valid;
    # transitions between different occupancy levels are handled by the ping-pong case.
    met_full_atoms = Dict{Symbol, Vector{Pair{Symbol,Int}}}()
    for spec in Iterators.flatten((substrates(reaction), products(reaction), regulators(reaction)))
        met_full_atoms[spec[1]] = sort([a => c for (a, c) in spec[2]]; by=first)
    end

    # Identify which sites are core substrate/product (index 1, non-regulator)
    sub_names = Set(s[1] for s in substrates(reaction))
    prod_names = Set(p[1] for p in products(reaction))

    # Number of catalytic site positions (substrates + products).
    # Sites at positions > n_catalytic are regulator sites.
    n_catalytic = _n_catalytic_sites(reaction)

    for i in 1:n, j in (i+1):n
        fi, fj = forms[i], forms[j]
        diff_idx = 0
        n_diff = 0
        for k in 1:nsites
            if fi.sites[k].atoms != fj.sites[k].atoms
                n_diff += 1
                n_diff > 1 && break
                diff_idx = k
            end
        end

        if n_diff == 1
            # Single-site difference → binding/release (standard or ping-pong)
            si, sj = fi.sites[diff_idx], fj.sites[diff_idx]
            if si.atoms === nothing && sj.atoms !== nothing
                # i is unoccupied, j is occupied → i + M → j (binding)
                # Only valid when j has full metabolite atoms (not a ping-pong residual)
                if sj.atoms == met_full_atoms[sj.metabolite]
                    push!(edges, ReactionEdge(i, j, sj.metabolite, :binding))
                    push!(edges, ReactionEdge(j, i, sj.metabolite, :release))
                end
            elseif si.atoms !== nothing && sj.atoms === nothing
                # i is occupied, j is unoccupied → i → j + M (release)
                # Only valid when i has full metabolite atoms (not a ping-pong residual)
                if si.atoms == met_full_atoms[si.metabolite]
                    push!(edges, ReactionEdge(i, j, si.metabolite, :release))
                    push!(edges, ReactionEdge(j, i, si.metabolite, :binding))
                end
            elseif si.atoms !== nothing && sj.atoms !== nothing
                # Both occupied with different atoms → ping-pong partial release
                ai = Dict{Symbol,Int}(a => c for (a, c) in si.atoms)
                aj = Dict{Symbol,Int}(a => c for (a, c) in sj.atoms)
                # Try i→j: atom decrease (release direction)
                diff_ij = Dict{Symbol,Int}()
                for (a, c) in ai
                    d = c - get(aj, a, 0)
                    d != 0 && (diff_ij[a] = d)
                end
                for (a, c) in aj
                    haskey(ai, a) && continue
                    diff_ij[a] = -c
                end
                if all(v > 0 for v in values(diff_ij))
                    # i has more atoms → i releases metabolite to get j
                    diff_sorted = sort([a => c for (a, c) in diff_ij]; by=first)
                    met = get(met_lookup, diff_sorted, nothing)
                    if met !== nothing
                        push!(edges, ReactionEdge(i, j, met, :release))
                        push!(edges, ReactionEdge(j, i, met, :binding))
                    end
                elseif all(v < 0 for v in values(diff_ij))
                    # j has more atoms → j releases metabolite to get i
                    diff_sorted = sort([a => -c for (a, c) in diff_ij]; by=first)
                    met = get(met_lookup, diff_sorted, nothing)
                    if met !== nothing
                        push!(edges, ReactionEdge(j, i, met, :release))
                        push!(edges, ReactionEdge(i, j, met, :binding))
                    end
                end
            end
        elseif n_diff >= 2
            # Check for isomerization: all diffs must be at catalytic site positions
            # (positions 1:n_catalytic), with index 1 and sub/prod metabolite.
            # Regulator sites (positions > n_catalytic) must not change.
            all_core = true
            for k in 1:nsites
                fi.sites[k].atoms == fj.sites[k].atoms && continue
                if k > n_catalytic
                    all_core = false
                    break
                end
                site = fi.sites[k]
                if site.index != 1 || !(site.metabolite in sub_names || site.metabolite in prod_names)
                    all_core = false
                    break
                end
            end
            if all_core
                # Check isomerization rule: F1 has all sub sites occupied + all prod sites empty,
                # F2 has all sub sites empty + all prod sites occupied (or vice versa).
                # Only check catalytic site positions (1:n_catalytic) — regulator sites
                # may share metabolite names with products but are not part of catalysis.
                i_sub_occ = true; i_prod_empty = true
                j_sub_occ = true; j_prod_empty = true
                i_sub_empty = true; i_prod_occ = true
                j_sub_empty = true; j_prod_occ = true
                for k in 1:n_catalytic
                    site = fi.sites[k]
                    site.index != 1 && continue
                    if site.metabolite in sub_names
                        fi.sites[k].atoms === nothing && (i_sub_occ = false)
                        fi.sites[k].atoms !== nothing && (i_sub_empty = false)
                        fj.sites[k].atoms === nothing && (j_sub_occ = false)
                        fj.sites[k].atoms !== nothing && (j_sub_empty = false)
                    elseif site.metabolite in prod_names
                        fi.sites[k].atoms !== nothing && (i_prod_empty = false)
                        fi.sites[k].atoms === nothing && (i_prod_occ = false)
                        fj.sites[k].atoms !== nothing && (j_prod_empty = false)
                        fj.sites[k].atoms === nothing && (j_prod_occ = false)
                    end
                end
                valid_iso = (i_sub_occ && i_prod_empty && j_sub_empty && j_prod_occ) ||
                            (i_sub_empty && i_prod_occ && j_sub_occ && j_prod_empty)
                if valid_iso
                    # Verify atom conservation across core catalytic sites only
                    atoms_i = _core_atoms(fi, n_catalytic)
                    atoms_j = _core_atoms(fj, n_catalytic)
                    if atoms_i == atoms_j
                        push!(edges, ReactionEdge(i, j, nothing, :isomerization))
                        push!(edges, ReactionEdge(j, i, nothing, :isomerization))
                    end
                end
            end
        end
    end
    edges
end

"""Build adjacency list: for each form, list of outgoing edge indices."""
function _build_adjacency(forms::Vector{EnzymeFormSpec}, edges::Vector{ReactionEdge})
    adj = [Int[] for _ in 1:length(forms)]
    for (idx, e) in enumerate(edges)
        push!(adj[e.from], idx)
    end
    adj
end

# ─── Cycle Enumeration ───────────────────────────────────────────────────────

"""
Find all simple directed 1× cycles through the free enzyme.

A 1× cycle has net stoichiometry: -1 for each substrate, +1 for each product, 0 for regulators.
Cycles are limited to at most `max_cycle_forms` forms to prevent combinatorial explosion
when regulators create many enzyme forms.
"""
function _find_valid_cycles(forms::Vector{EnzymeFormSpec}, edges::Vector{ReactionEdge},
                            adj::Vector{Vector{Int}}, max_cycle_forms::Int,
                            @nospecialize(reaction::EnzymeReaction))
    free_idx = _free_enzyme_index(forms)
    free_idx === nothing && return Vector{ReactionEdge}[]

    expected_stoich = _expected_stoichiometry(reaction)

    cycles = Vector{ReactionEdge}[]
    visited = falses(length(forms))
    path_edges = ReactionEdge[]

    function dfs(node::Int)
        # path_edges has N edges → N+1 forms visited (including free_idx).
        # Closing back to free_idx adds one edge but no new form.
        # So a cycle with K forms has K edges, and path_edges has K-1 before closing.
        length(path_edges) >= max_cycle_forms && return
        for ei in adj[node]
            e = edges[ei]
            if e.to == free_idx && length(path_edges) >= 2
                push!(path_edges, e)
                stoich = _cycle_stoichiometry(path_edges)
                _is_1x_stoich(stoich, expected_stoich) && push!(cycles, copy(path_edges))
                pop!(path_edges)
            elseif !visited[e.to] && e.to != free_idx
                visited[e.to] = true
                push!(path_edges, e)
                dfs(e.to)
                pop!(path_edges)
                visited[e.to] = false
            end
        end
    end

    visited[free_idx] = true
    dfs(free_idx)
    cycles
end

# ─── Topology Combination ────────────────────────────────────────────────────

"""Compute cycle stoichiometry: returns Dict{Symbol,Int} of net metabolite changes."""
function _cycle_stoichiometry(cycle_edges::Vector{ReactionEdge})
    stoich = Dict{Symbol,Int}()
    for e in cycle_edges
        e.metabolite === nothing && continue
        if e.edge_type == :binding
            stoich[e.metabolite] = get(stoich, e.metabolite, 0) - 1
        elseif e.edge_type == :release
            stoich[e.metabolite] = get(stoich, e.metabolite, 0) + 1
        end
    end
    stoich
end

"""Check if stoichiometry exactly matches expected (1× cycle)."""
function _is_1x_stoich(stoich::Dict{Symbol,Int}, expected::Dict{Symbol,Int})
    for (met, exp) in expected
        get(stoich, met, 0) != exp && return false
    end
    for (met, _) in stoich
        haskey(expected, met) || return false
    end
    true
end

"""Check if stoichiometry is 1× (matches expected) or 0× (all zero)."""
function _is_valid_stoich(stoich::Dict{Symbol,Int}, expected::Dict{Symbol,Int})
    all_zero = all(v == 0 for v in values(stoich)) && all(get(stoich, m, 0) == 0 for m in keys(expected))
    all_zero && return true
    _is_1x_stoich(stoich, expected)
end

"""
Find all simple cycles through start_node in a subgraph defined by form_set and edge_list.
Returns true if all cycles have valid (1× or 0×) stoichiometry, false if any invalid cycle found.
"""
function _validate_all_cycles(start_node::Int, form_set::Set{Int},
                              sub_adj::Dict{Int, Vector{ReactionEdge}},
                              expected::Dict{Symbol,Int})
    visited = Set{Int}()
    path_edges = ReactionEdge[]

    function dfs(node::Int)
        for e in get(sub_adj, node, ReactionEdge[])
            e.to ∉ form_set && continue
            if e.to == start_node && length(path_edges) >= 2
                push!(path_edges, e)
                stoich = _cycle_stoichiometry(path_edges)
                valid = _is_valid_stoich(stoich, expected)
                pop!(path_edges)
                valid || return false
            elseif e.to ∉ visited && e.to != start_node
                push!(visited, e.to)
                push!(path_edges, e)
                result = dfs(e.to)
                result || return false
                pop!(path_edges)
                delete!(visited, e.to)
            end
        end
        return true
    end

    push!(visited, start_node)
    dfs(start_node)
end

"""
Combine 1× cycles into valid topologies via incremental BFS.
"""
function _combine_cycles(cycles::Vector{Vector{ReactionEdge}},
                         forms::Vector{EnzymeFormSpec},
                         all_edges::Vector{ReactionEdge},
                         max_forms::Int,
                         @nospecialize(reaction::EnzymeReaction))
    isempty(cycles) && return Topology[]

    expected_stoich = _expected_stoichiometry(reaction)
    free_idx = _free_enzyme_index(forms)::Int

    # Extract (form_set, undirected_edge_set) for each cycle
    UEdge = Tuple{Int,Int,Union{Nothing,Symbol}}
    cycle_data = map(cycles) do cycle
        fset = Set{Int}()
        uset = Set{UEdge}()
        for e in cycle
            push!(fset, e.from)
            push!(fset, e.to)
            push!(uset, _undirected_key(e))
        end
        (fset, uset)
    end

    # Canonical key for dedup: sorted form set + sorted undirected edge set
    function canonical_key(fset, uset)
        fs = sort!(collect(fset))
        us = sort!(collect(uset))
        (fs, us)
    end

    seen = Set{UInt64}()
    topologies = Topology[]

    # BFS queue: each entry is (form_set, undirected_edge_set, directed_edges, max_cycle_idx_used)
    queue = Vector{Tuple{Set{Int}, Set{UEdge}, Vector{ReactionEdge}, Int}}()

    # Seed with individual cycles
    for (ci, (fset, uset)) in enumerate(cycle_data)
        length(fset) > max_forms && continue
        key = canonical_key(fset, uset)
        h = hash(key)
        h in seen && continue
        push!(seen, h)

        dir_edges = _directed_edges_from_cycle(cycles[ci])
        push!(topologies, Topology(sort!(collect(fset)), dir_edges))
        push!(queue, (fset, uset, dir_edges, ci))
    end

    # BFS: try adding more cycles
    while !isempty(queue)
        fset, uset, dir_edges, max_ci = popfirst!(queue)
        for ci in (max_ci+1):length(cycles)
            cfset, cuset = cycle_data[ci]
            new_fset = union(fset, cfset)
            length(new_fset) > max_forms && continue
            new_uset = union(uset, cuset)
            key = canonical_key(new_fset, new_uset)
            h = hash(key)
            h in seen && continue
            push!(seen, h)

            # Build subgraph adjacency for cycle validation (both directions)
            sub_adj = Dict{Int, Vector{ReactionEdge}}()
            for e in all_edges
                _undirected_key(e) ∉ new_uset && continue
                e.from ∉ new_fset && continue
                e.to ∉ new_fset && continue
                push!(get!(sub_adj, e.from, ReactionEdge[]), e)
            end
            valid = _validate_all_cycles(free_idx, new_fset, sub_adj, expected_stoich)
            valid || continue

            new_dir_edges = _merge_cycle_edges(dir_edges, cycles[ci])
            push!(topologies, Topology(sort!(collect(new_fset)), new_dir_edges))
            push!(queue, (new_fset, new_uset, new_dir_edges, ci))
        end
    end

    topologies
end

"""Collect ONE directed edge per undirected key from a cycle's edges."""
function _directed_edges_from_cycle(cycle::Vector{ReactionEdge})
    seen = Set{Tuple{Int,Int,Union{Nothing,Symbol}}}()
    result = ReactionEdge[]
    for e in cycle
        ukey = _undirected_key(e)
        ukey in seen && continue
        push!(seen, ukey)
        push!(result, e)
    end
    result
end

"""Merge directed edges from a new cycle into existing directed edge list.
New edges (by undirected key) get the cycle's direction; existing edges keep their direction."""
function _merge_cycle_edges(existing::Vector{ReactionEdge}, cycle::Vector{ReactionEdge})
    existing_keys = Set(_undirected_key(e) for e in existing)
    result = copy(existing)
    for e in cycle
        ukey = _undirected_key(e)
        ukey in existing_keys && continue
        push!(existing_keys, ukey)
        push!(result, e)
    end
    result
end

# ─── Dead-End Computation ────────────────────────────────────────────────────

"""Return true if the binding edge from `from` to `to` adds a catalytic metabolite at its
catalytic site (index 1). Such edges should not form dead-ends — dead-end inhibition
requires binding at a non-catalytic (allosteric) site."""
function _is_catalytic_site_binding(from::EnzymeFormSpec, to::EnzymeFormSpec,
                                     catalytic_mets::Set{Symbol})
    for (sf, st) in zip(from.sites, to.sites)
        sf.atoms === st.atoms && continue  # unchanged site
        # st is the newly occupied site (binding edge: from→to adds a metabolite)
        return st.metabolite in catalytic_mets && st.index == 1
    end
    return false
end

"""
Build dead-end trees for each cycle form in a topology.
Returns Vector{Vector{DeadEndTree}} — one vector of trees per cycle form.
"""
function _build_dead_end_trees(topo::Topology, all_forms::Vector{EnzymeFormSpec},
                               adj::Vector{Vector{Int}}, edges::Vector{ReactionEdge},
                               catalytic_mets::Set{Symbol})
    cycle_set = Set(topo.form_indices)
    trees_per_form = Vector{Vector{DeadEndTree}}()

    for fi in topo.form_indices
        # Find binding edges from fi to forms NOT in cycle
        roots = DeadEndTree[]
        for ei in adj[fi]
            e = edges[ei]
            e.edge_type == :binding || continue
            e.to in cycle_set && continue
            _is_catalytic_site_binding(all_forms[e.from], all_forms[e.to], catalytic_mets) && continue
            tree = _build_dead_end_subtree(e.to, cycle_set, all_forms, adj, edges,
                                            catalytic_mets, Set{Int}())
            push!(roots, tree)
        end
        push!(trees_per_form, roots)
    end
    trees_per_form
end

"""Recursively build a dead-end subtree rooted at form_idx."""
function _build_dead_end_subtree(form_idx::Int, excluded::Set{Int},
                                  all_forms::Vector{EnzymeFormSpec},
                                  adj::Vector{Vector{Int}}, edges::Vector{ReactionEdge},
                                  catalytic_mets::Set{Symbol}, visited::Set{Int})
    push!(visited, form_idx)
    children = DeadEndTree[]
    for ei in adj[form_idx]
        e = edges[ei]
        e.edge_type == :binding || continue
        e.to in excluded && continue
        e.to in visited && continue
        _is_catalytic_site_binding(all_forms[e.from], all_forms[e.to], catalytic_mets) && continue
        child = _build_dead_end_subtree(e.to, excluded, all_forms, adj, edges,
                                         catalytic_mets, visited)
        push!(children, child)
    end
    delete!(visited, form_idx)
    DeadEndTree(form_idx, children)
end

"""
Enumerate all valid dead-end configurations for an entire topology,
respecting global max_forms budget.

Returns Vector{Vector{Int}} where each inner vector is the set of dead-end form indices to add.
"""
function _enumerate_topology_dead_ends(topo::Topology, dead_end_trees::Vector{Vector{DeadEndTree}},
                                        max_forms::Int)
    budget = max_forms - length(topo.form_indices)
    budget < 0 && return Vector{Int}[]

    # For each cycle form, enumerate its valid dead-end configs (respecting budget)
    # Then take constrained Cartesian product across cycle forms
    per_form_options = [_enumerate_single_form_dead_ends(trees, budget) for trees in dead_end_trees]

    # Constrained Cartesian product: total dead-end forms across all cycle forms ≤ budget, no duplicates
    result = Vector{Int}[]
    _cartesian_dead_ends(per_form_options, 1, Int[], budget, Set{Int}(), result)
    result
end

"""
Return all downward-closed subsets of a dead-end tree (including empty set).
Each subset is a Vector{Int} of form indices.
"""
function _dc_subsets(tree::DeadEndTree)
    # DC subsets that include root = {root} × Cartesian product of child DC subsets
    child_sublists = Vector{Vector{Int}}[_dc_subsets(c) for c in tree.children]
    including_root = Vector{Int}[]
    _cartesian_dc(child_sublists, 1, [tree.form_idx], including_root)
    pushfirst!(including_root, Int[])  # empty subset = don't include this root
    including_root
end

"""Cartesian product of DC subset lists from each child subtree."""
function _cartesian_dc(sublists::Vector{Vector{Vector{Int}}}, idx::Int,
                        current::Vector{Int}, results::Vector{Vector{Int}})
    if idx > length(sublists)
        push!(results, copy(current))
        return
    end
    for subset in sublists[idx]
        append!(current, subset)
        _cartesian_dc(sublists, idx + 1, current, results)
        resize!(current, length(current) - length(subset))
    end
end

"""Enumerate dead-end configs for a single cycle form's trees, up to budget forms."""
function _enumerate_single_form_dead_ends(trees::Vector{DeadEndTree}, budget::Int)
    per_root = Vector{Vector{Int}}[_dc_subsets(tree) for tree in trees]
    configs = Vector{Int}[]
    _cartesian_dc_budget(per_root, 1, Int[], budget, configs)
    configs
end

"""Cartesian product of per-root DC subsets, respecting total form budget."""
function _cartesian_dc_budget(per_root::Vector{Vector{Vector{Int}}}, idx::Int,
                               current::Vector{Int}, budget::Int,
                               configs::Vector{Vector{Int}})
    if idx > length(per_root)
        push!(configs, copy(current))
        return
    end
    for subset in per_root[idx]
        length(subset) > budget && continue
        append!(current, subset)
        _cartesian_dc_budget(per_root, idx + 1, current, budget - length(subset), configs)
        resize!(current, length(current) - length(subset))
    end
end

"""Constrained Cartesian product across cycle forms, respecting total budget and no duplicates."""
function _cartesian_dead_ends(per_form::Vector{Vector{Vector{Int}}}, form_idx::Int,
                               current::Vector{Int}, budget::Int,
                               used::Set{Int},
                               result::Vector{Vector{Int}})
    if form_idx > length(per_form)
        push!(result, copy(current))
        return
    end
    for config in per_form[form_idx]
        length(config) > budget && continue
        # Skip configs that include forms already used by other cycle forms
        any(f in used for f in config) && continue
        for f in config; push!(used, f); end
        append!(current, config)
        _cartesian_dead_ends(per_form, form_idx + 1, current, budget - length(config), used, result)
        resize!(current, length(current) - length(config))
        for f in config; delete!(used, f); end
    end
end

# ─── Equivalent Step Detection ───────────────────────────────────────────────

"""
Find groups of equivalent steps in a mechanism.

Two steps are equivalent if both are binding edges for the same metabolite
at the same site index, differing only in enzyme state.

Returns Vector{Vector{Int}} where each inner vector is a group of step indices (1-based).
"""
function _find_equivalent_groups(step_edges::Vector{ReactionEdge},
                                  forms::Vector{EnzymeFormSpec})
    # Map each step to (metabolite, site_index) for binding edges
    binding_key = Dict{Tuple{Symbol, Int}, Vector{Int}}()
    for (i, e) in enumerate(step_edges)
        e.edge_type == :binding || continue
        # Find which site differs between from and to
        f_from = forms[e.from]
        f_to = forms[e.to]
        for k in 1:length(f_from.sites)
            if f_from.sites[k].atoms === nothing && f_to.sites[k].atoms !== nothing
                key = (f_to.sites[k].metabolite, f_to.sites[k].index)
                push!(get!(binding_key, key, Int[]), i)
                break
            end
        end
    end
    groups = Vector{Int}[]
    for (_, indices) in binding_key
        length(indices) >= 2 && push!(groups, sort(indices))
    end
    sort!(groups; by=first)
    groups
end

# ─── MechanismSpec Materialization ───────────────────────────────────────────

"""
Build a MechanismSpec from a topology + dead-end config + RE/SS assignment + constraint choice.
"""
function _build_mechanism_spec(topo::Topology, dead_end_forms::Vector{Int},
                                eq_steps::Vector{Bool},
                                equiv_groups::Vector{Vector{Int}},
                                constraint_mask::Int,
                                all_forms::Vector{EnzymeFormSpec},
                                all_edges::Vector{ReactionEdge},
                                adj::Vector{Vector{Int}},
                                @nospecialize(reaction))
    form_indices = copy(topo.form_indices)
    append!(form_indices, dead_end_forms)

    step_edges = _get_step_edges(topo, dead_end_forms, all_forms, all_edges, adj)

    form_names = [all_forms[fi].name for fi in form_indices]
    form_atoms_list = [_form_atoms(all_forms[fi].sites) for fi in form_indices]

    rxn_tuples = Tuple{Vector{Symbol}, Vector{Symbol}}[]
    for e in step_edges
        from_name = all_forms[e.from].name
        to_name = all_forms[e.to].name
        if e.edge_type == :binding
            lhs = e.metabolite === nothing ? [from_name] : [from_name, e.metabolite]
            rhs = [to_name]
        elseif e.edge_type == :release
            lhs = [from_name]
            rhs = e.metabolite === nothing ? [to_name] : [to_name, e.metabolite]
        else  # isomerization
            lhs = [from_name]
            rhs = [to_name]
        end
        push!(rxn_tuples, (lhs, rhs))
    end

    # Build param constraints from equivalent groups + constraint mask
    constraints = Tuple{Symbol, Int, Vector{Tuple{Symbol,Int}}}[]
    for (gi, group) in enumerate(equiv_groups)
        ((constraint_mask >> (gi - 1)) & 1) == 0 && continue
        # All steps in group should have same RE/SS assignment
        first_step = group[1]
        is_re = eq_steps[first_step]
        # Check all steps in group have same RE/SS
        all_same = all(eq_steps[s] == is_re for s in group)
        all_same || continue

        if is_re
            # Constrain K_j = K_i for j > 1 in group
            for j in 2:length(group)
                push!(constraints, (Symbol("K$(group[j])"), 1,
                                   [(Symbol("K$(group[1])"), 1)]))
            end
        else
            # Constrain k_jf = k_if, k_jr = k_ir
            for j in 2:length(group)
                push!(constraints, (Symbol("k$(group[j])f"), 1,
                                   [(Symbol("k$(group[1])f"), 1)]))
                push!(constraints, (Symbol("k$(group[j])r"), 1,
                                   [(Symbol("k$(group[1])r"), 1)]))
            end
        end
    end

    MechanismSpec(reaction, form_names, form_atoms_list, rxn_tuples, eq_steps, constraints)
end

"""Compute total atom content for an enzyme form from its sites."""
function _form_atoms(sites::Vector{SiteState})
    atoms = Dict{Symbol,Int}()
    for site in sites
        site.atoms === nothing && continue
        for (a, c) in site.atoms
            atoms[a] = get(atoms, a, 0) + c
        end
    end
    sort([a => c for (a, c) in atoms]; by=first)
end

# ─── Iterator Implementation ─────────────────────────────────────────────────

function Base.length(iter::MechanismIterator)
    total = 0
    for (ti, topo) in enumerate(iter.topologies)
        dead_end_configs = _enumerate_topology_dead_ends(topo, iter.dead_end_trees[ti], iter.max_forms)
        for de_config in dead_end_configs
            step_edges = _get_step_edges(topo, de_config, iter.all_forms,
                                        iter.reaction_graph, iter.adjacency)
            n_steps = length(step_edges)
            n_ress = (1 << n_steps) - 1

            # For each RE/SS, compute equivalent groups and constraint variants
            for ress_idx in 1:n_ress
                eq_steps = _ress_index_to_vec(ress_idx, n_steps)
                equiv_groups = _find_equivalent_groups(step_edges, iter.all_forms)
                valid_groups = [g for g in equiv_groups if all(eq_steps[s] == eq_steps[g[1]] for s in g)]
                n_constraints = 1 << length(valid_groups)
                total += n_constraints
            end
        end
    end
    total
end

"""Convert RE/SS index (1-based) to Bool vector. Bit=1 means RE (true).
idx ranges from 1 to 2^n - 1, mapping to bit patterns 0 through 2^n - 2
(skipping all-RE = 2^n - 1)."""
function _ress_index_to_vec(idx::Int, n_steps::Int)
    bits = idx - 1  # 0-based: idx=1 → bits=0 (all SS), idx=2 → bits=1, etc.
    Bool[((bits >> (i-1)) & 1) == 1 for i in 1:n_steps]
end

"""
Get canonical step edges for a topology + dead-end config.

Each undirected step is represented once, preferring binding direction for
bind/release pairs, and lower-index-first for isomerization.
Searches the full reaction graph (`all_edges`) for canonical directions.
"""
function _get_step_edges(topo::Topology, de_forms::Vector{Int},
                          all_forms::Vector{EnzymeFormSpec},
                          all_edges::Vector{ReactionEdge},
                          adj::Vector{Vector{Int}})
    form_set = Set(topo.form_indices)
    for fi in de_forms; push!(form_set, fi); end

    step_edges = ReactionEdge[]
    seen_pairs = Set{Tuple{Int,Int}}()

    # Topology edges: keep cycle-traversal direction (preserves net stoichiometry)
    for e in topo.edges
        pair = minmax(e.from, e.to)
        pair in seen_pairs && continue
        push!(seen_pairs, pair)
        push!(step_edges, e)
    end

    # Dead-end edges: find binding edges from mechanism forms to dead-end forms
    for de_idx in de_forms
        for fi in Iterators.flatten((topo.form_indices, de_forms))
            fi == de_idx && continue
            fi ∉ form_set && continue
            for ei in adj[fi]
                e = all_edges[ei]
                e.to == de_idx && e.edge_type == :binding || continue
                pair = minmax(e.from, e.to)
                pair in seen_pairs && continue
                push!(seen_pairs, pair)
                push!(step_edges, e)
                break
            end
        end
    end

    step_edges
end


function Base.iterate(iter::MechanismIterator, state::Union{Nothing, Any}=nothing)
    if state === nothing
        state = (1, 1, 1, 1, _precompute_topo_data(iter, 1)...)
    end

    topo_idx, de_idx, ress_idx, constraint_idx, dead_end_configs, step_edges_cache = state

    while topo_idx <= length(iter.topologies)
        topo = iter.topologies[topo_idx]

        while de_idx <= length(dead_end_configs)
            de_config = dead_end_configs[de_idx]
            cached_step_edges = get!(step_edges_cache, de_idx) do
                _get_step_edges(topo, de_config, iter.all_forms,
                               iter.reaction_graph, iter.adjacency)
            end
            n_steps = length(cached_step_edges)
            n_ress = (1 << n_steps) - 1

            while ress_idx <= n_ress
                eq_steps_vec = _ress_index_to_vec(ress_idx, n_steps)
                equiv_groups = _find_equivalent_groups(cached_step_edges, iter.all_forms)
                valid_groups = [g for g in equiv_groups
                                if all(eq_steps_vec[s] == eq_steps_vec[g[1]] for s in g)]
                n_constraints = 1 << length(valid_groups)

                while constraint_idx <= n_constraints
                    cmask = constraint_idx - 1
                    spec = _build_mechanism_spec(topo, de_config, eq_steps_vec,
                                                 valid_groups, cmask,
                                                 iter.all_forms, iter.reaction_graph,
                                                 iter.adjacency, iter.reaction)

                    # Advance to next state before returning
                    next_state = _advance_state(iter, topo_idx, de_idx, ress_idx,
                                                constraint_idx + 1, n_constraints, n_ress,
                                                dead_end_configs, step_edges_cache)
                    return (spec, next_state)
                end
                constraint_idx = 1
                ress_idx += 1
            end
            ress_idx = 1
            de_idx += 1
        end
        de_idx = 1
        topo_idx += 1
        if topo_idx <= length(iter.topologies)
            dead_end_configs, step_edges_cache = _precompute_topo_data(iter, topo_idx)
        end
    end
    nothing
end

"""Compute next iterator state after yielding, cascading overflows upward."""
function _advance_state(iter, topo_idx, de_idx, ress_idx, constraint_idx,
                        n_constraints, n_ress, dead_end_configs, step_edges_cache)
    if constraint_idx <= n_constraints
        return (topo_idx, de_idx, ress_idx, constraint_idx,
                dead_end_configs, step_edges_cache)
    end
    # Overflow constraint → advance RE/SS
    ress_idx += 1
    if ress_idx <= n_ress
        return (topo_idx, de_idx, ress_idx, 1,
                dead_end_configs, step_edges_cache)
    end
    # Overflow RE/SS → advance dead-end config
    de_idx += 1
    if de_idx <= length(dead_end_configs)
        return (topo_idx, de_idx, 1, 1,
                dead_end_configs, step_edges_cache)
    end
    # Overflow dead-end → advance topology
    topo_idx += 1
    if topo_idx <= length(iter.topologies)
        new_precomp = _precompute_topo_data(iter, topo_idx)
        return (topo_idx, 1, 1, 1, new_precomp...)
    end
    # Past end
    return (topo_idx, 1, 1, 1,
            Vector{Int}[], Dict{Int, Vector{ReactionEdge}}())
end

function _precompute_topo_data(iter::MechanismIterator, topo_idx::Int)
    topo = iter.topologies[topo_idx]
    dead_end_configs = _enumerate_topology_dead_ends(topo, iter.dead_end_trees[topo_idx], iter.max_forms)
    step_edges_cache = Dict{Int, Vector{ReactionEdge}}()
    (dead_end_configs, step_edges_cache)
end

# ─── Main Entry Point ────────────────────────────────────────────────────────

"""
    enumerate_mechanisms(reaction::EnzymeReaction; max_forms = 3 * n_sites(reaction))

Enumerate all valid mechanism topologies for the given reaction.

Returns a `MechanismIterator` that lazily yields `MechanismSpec` structs.
Each spec can be converted to an `EnzymeMechanism` via `EnzymeMechanism(spec)`.

The enumeration covers:
1. Catalytic topologies (valid cycle combinations through free enzyme)
2. Dead-end attachments (substrate/product/regulator-bound forms branching off cycles)
3. RE/SS assignments (rapid-equilibrium vs steady-state per step)
4. Equivalent step constraints (shared parameters for equivalent binding steps)
"""
function enumerate_mechanisms(@nospecialize(reaction::EnzymeReaction);
                              max_forms::Int = 3 * n_sites(reaction))
    forms = enumerate_enzyme_forms(reaction)
    edges = _build_reaction_graph(forms, reaction)
    adj = _build_adjacency(forms, edges)
    max_cf = min(_max_cycle_forms(forms, reaction), max_forms)
    cycles = _find_valid_cycles(forms, edges, adj, max_cf, reaction)
    topologies = _combine_cycles(cycles, forms, edges, max_forms, reaction)

    # Build set of catalytic metabolites (substrates + products) — dead-end inhibition
    # at catalytic sites (index 1) is disallowed; only non-catalytic sites (index 2+) allowed
    catalytic_mets = Set{Symbol}()
    for s in substrates(reaction); push!(catalytic_mets, s[1]); end
    for p in products(reaction); push!(catalytic_mets, p[1]); end

    # Precompute dead-end trees per topology
    dead_end_trees = [_build_dead_end_trees(topo, forms, adj, edges, catalytic_mets)
                      for topo in topologies]

    MechanismIterator(forms, edges, topologies, reaction, max_forms, dead_end_trees, adj)
end

# ─── MechanismSpec → EnzymeMechanism Conversion ─────────────────────────────

"""
    EnzymeMechanism(spec::MechanismSpec)

Convert a lightweight `MechanismSpec` to a type-parameterized `EnzymeMechanism`.
Delegates to the standard constructor for full validation.
"""
function EnzymeMechanism(spec::MechanismSpec)
    rxn = spec.reaction
    subs_t = Tuple((name, atoms) for (name, atoms, _) in substrates(rxn))
    prods_t = Tuple((name, atoms) for (name, atoms, _) in products(rxn))
    regs_t = Tuple((name, atoms) for (name, atoms, _) in regulators(rxn))
    enzs_t = Tuple((spec.forms[i], Tuple(Tuple.(spec.form_atoms[i]))) for i in 1:length(spec.forms))
    species = (subs_t, prods_t, regs_t, enzs_t)

    reactions = Tuple(
        (Tuple(r[1]), Tuple(r[2])) for r in spec.reactions
    )

    eq_steps = Tuple(spec.equilibrium_steps)

    constraints = Tuple(
        (target, coeff, Tuple(Tuple.(factors)))
        for (target, coeff, factors) in spec.param_constraints
    )

    EnzymeMechanism(species, reactions, eq_steps, constraints)
end
