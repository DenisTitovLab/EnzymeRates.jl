# ─── Trait Types ──────────────────────────────────────────────────

abstract type AbstractEdgeClass end
struct MustExist <: AbstractEdgeClass end
struct CouldExist <: AbstractEdgeClass end
struct Forbidden <: AbstractEdgeClass end

abstract type AbstractSiteOccupancy end
struct EmptySite <: AbstractSiteOccupancy end
struct OccupiedSite <: AbstractSiteOccupancy end
struct ResidualSite <: AbstractSiteOccupancy end

# ─── Enriched SiteState ───────────────────────────────────────────

"""
    SiteState

State of a single binding site on the enzyme.

- `metabolite`: which metabolite this site is for
- `index`: site number (1, 2, ...)
- `atoms`: atoms present (`nothing` = unoccupied,
  empty vector = occupied with no-atom species)
- `role`: `:sub`, `:prod`, or `:reg`
- `full_atoms`: the full atom vector for this site's metabolite
"""
struct SiteState
    metabolite::Symbol
    index::Int
    atoms::Union{Nothing, Vector{Pair{Symbol,Int}}}
    role::Symbol
    full_atoms::Vector{Pair{Symbol,Int}}
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

const ParamConstraint =
    Tuple{Symbol, Int, Vector{Tuple{Symbol, Int}}}

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
    param_constraints::Vector{ParamConstraint}
end

"""
    MechanismIterator

Lazy iterator over all valid mechanisms for a reaction.
"""
struct MechanismIterator
    specs::Vector{MechanismSpec}
end

Base.eltype(::Type{MechanismIterator}) = MechanismSpec
Base.IteratorSize(::Type{MechanismIterator}) = Base.HasLength()
Base.length(iter::MechanismIterator) = length(iter.specs)
Base.iterate(iter::MechanismIterator, state=1) =
    state > length(iter.specs) ? nothing : (iter.specs[state], state + 1)

# ─── Core Helpers ─────────────────────────────────────────────────

"""
    site_occupancy(site::SiteState) → AbstractSiteOccupancy

Determine occupancy of a single site:
- `nothing` atoms → `EmptySite()`
- atoms == full_atoms → `OccupiedSite()`
- otherwise → `ResidualSite()` (ping-pong partial)
"""
function site_occupancy(site::SiteState)
    site.atoms === nothing && return EmptySite()
    site.atoms == site.full_atoms && return OccupiedSite()
    return ResidualSite()
end

"""
    edge_class(form_a, form_b) → (AbstractEdgeClass, metabolite, edge_type)

Classify the edge between two enzyme forms based on site occupancy diffs.

Returns `(Forbidden(), nothing, :none)` if no valid edge exists.
"""
function edge_class(fa::EnzymeFormSpec, fb::EnzymeFormSpec)
    nsites = length(fa.sites)
    # Find differing sites
    diff_positions = Int[]
    for k in 1:nsites
        if fa.sites[k].atoms != fb.sites[k].atoms
            push!(diff_positions, k)
        end
    end

    isempty(diff_positions) && return (Forbidden(), nothing, :none)

    if length(diff_positions) == 1
        k = diff_positions[1]
        sa, sb = fa.sites[k], fb.sites[k]
        oa, ob = site_occupancy(sa), site_occupancy(sb)
        is_core = sa.role in (:sub, :prod) && sa.index == 1

        if is_core
            # Empty↔Occupied: binding/release of metabolite
            if oa isa EmptySite && ob isa OccupiedSite
                return (MustExist(), sa.metabolite, :binding)
            elseif oa isa OccupiedSite && ob isa EmptySite
                return (MustExist(), sa.metabolite, :release)
            # Occupied↔Residual: ping-pong product release
            elseif oa isa OccupiedSite && ob isa ResidualSite
                # Find which product was released by computing atom diff
                met = _released_metabolite(sa.atoms, sb.atoms, fa)
                met === nothing && return (Forbidden(), nothing, :none)
                return (MustExist(), met, :release)
            elseif oa isa ResidualSite && ob isa OccupiedSite
                met = _released_metabolite(sb.atoms, sa.atoms, fa)
                met === nothing && return (Forbidden(), nothing, :none)
                return (MustExist(), met, :binding)
            # Residual↔Empty: release of remaining atoms
            elseif oa isa ResidualSite && ob isa EmptySite
                met = _residual_metabolite(sa.atoms, fa)
                met === nothing && return (Forbidden(), nothing, :none)
                return (MustExist(), met, :release)
            elseif oa isa EmptySite && ob isa ResidualSite
                return (Forbidden(), nothing, :none)
            else
                return (Forbidden(), nothing, :none)
            end
        else
            # Regulatory or extra site
            if oa isa EmptySite && ob isa OccupiedSite
                return (CouldExist(), sa.metabolite, :binding)
            elseif oa isa OccupiedSite && ob isa EmptySite
                return (CouldExist(), sa.metabolite, :release)
            else
                return (Forbidden(), nothing, :none)
            end
        end
    end

    # ≥2 diffs: check for isomerization (all diffs at core sites)
    all_core = all(diff_positions) do k
        fa.sites[k].role in (:sub, :prod) && fa.sites[k].index == 1
    end
    !all_core && return (Forbidden(), nothing, :none)

    # Check isomerization pattern
    if _is_valid_isomerization(fa, fb)
        return (MustExist(), nothing, :isomerization)
    end
    return (Forbidden(), nothing, :none)
end

"""
    _released_metabolite(more, fewer, form)

Find the metabolite released when atoms decrease
from `more` to `fewer` at a site.
"""
function _released_metabolite(
    more::Vector{Pair{Symbol,Int}},
    fewer::Vector{Pair{Symbol,Int}},
    form::EnzymeFormSpec,
)
    da = Dict{Symbol,Int}(a => c for (a, c) in more)
    db = Dict{Symbol,Int}(a => c for (a, c) in fewer)
    diff = Dict{Symbol,Int}()
    for (a, c) in da
        d = c - get(db, a, 0)
        d != 0 && (diff[a] = d)
    end
    for (a, c) in db
        haskey(da, a) && continue
        diff[a] = -c
    end
    all(v > 0 for v in values(diff)) || return nothing
    diff_sorted = sort([a => c for (a, c) in diff]; by=first)
    # Match against product full_atoms
    for s in form.sites
        if s.role == :prod && s.index == 1 && s.full_atoms == diff_sorted
            return s.metabolite
        end
    end
    nothing
end

"""Find the metabolite whose atoms match a residual content exactly."""
function _residual_metabolite(atoms::Vector{Pair{Symbol,Int}}, form::EnzymeFormSpec)
    # The residual atoms at this site match some product's atoms
    for s in form.sites
        s.role == :prod && s.index == 1 && s.full_atoms == atoms && return s.metabolite
    end
    nothing
end

"""Check if two forms represent a valid isomerization (all-subs ↔ all-prods)."""
function _is_valid_isomerization(fa::EnzymeFormSpec, fb::EnzymeFormSpec)
    # Classify core site occupancy pattern for a form
    function _core_pattern(form)
        sub_occ = true;  sub_empty = true
        prod_occ = true; prod_empty = true
        for s in form.sites
            s.role in (:sub, :prod) && s.index == 1 || continue
            if s.role == :sub
                s.atoms === nothing && (sub_occ = false)
                s.atoms !== nothing && (sub_empty = false)
            else
                s.atoms === nothing && (prod_occ = false)
                s.atoms !== nothing && (prod_empty = false)
            end
        end
        (; sub_occ, sub_empty, prod_occ, prod_empty)
    end

    a = _core_pattern(fa)
    b = _core_pattern(fb)
    # One form has all-subs-occupied + all-prods-empty,
    # the other has all-subs-empty + all-prods-occupied
    valid = (a.sub_occ && a.prod_empty &&
             b.sub_empty && b.prod_occ) ||
            (a.sub_empty && a.prod_occ &&
             b.sub_occ && b.prod_empty)
    !valid && return false

    _core_atoms(fa) == _core_atoms(fb)
end

"""Compute total atoms across core catalytic sites (role ∈ (:sub,:prod), index == 1)."""
function _core_atoms(form::EnzymeFormSpec)
    atoms = Dict{Symbol,Int}()
    for s in form.sites
        s.role in (:sub, :prod) && s.index == 1 || continue
        s.atoms === nothing && continue
        for (a, c) in s.atoms
            atoms[a] = get(atoms, a, 0) + c
        end
    end
    atoms
end

# ─── Enzyme Form Enumeration ──────────────────────────────────────

"""
    enumerate_enzyme_forms(reaction::EnzymeReaction)

Enumerate all possible enzyme forms for the given reaction.

Each binding site is distinguishable. Forms include standard forms (each site
independently empty or fully occupied) and ping-pong intermediates (sites with
partial residual atom content after product release).

**Exclusion rule**: Forms where all substrate sites are fully bound while any
product site is occupied (or vice versa) are excluded as physically impossible.

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
                (prod_atoms_list[i]
                 for i in 1:length(prod_atoms_list)
                 if (mask >> (i-1)) & 1 == 1))
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
    roles = Symbol[]
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
    for s in S
        add!(s[1], 1, fatoms(s), :sub)
    end
    for p in P
        add!(p[1], 1, fatoms(p), :prod)
    end
    for s in S, i in 2:nsites(s)
        add!(s[1], i, fatoms(s), :sub)
    end
    for p in P, i in 2:nsites(p)
        add!(p[1], i, fatoms(p), :prod)
    end
    for r in R, i in 1:nsites(r)
        add!(r[1], i, fatoms(r), :reg)
    end

    # 3. Enumerate Cartesian product with exclusion filter
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
                                 content === nothing ? nothing : copy(content),
                                 roles[i], fulls[i])
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

# ─── n_sites ──────────────────────────────────────────────────────

"""
    n_sites(reaction::EnzymeReaction)

Total number of binding sites across all substrates, products, and regulators.
"""
function n_sites(@nospecialize(reaction::EnzymeReaction))
    count = 0
    all_specs = Iterators.flatten((
        substrates(reaction),
        products(reaction),
        regulators(reaction),
    ))
    for spec in all_specs
        count += length(spec) >= 3 ? spec[3] : 1
    end
    count
end

# ─── Catalytic Cycle Construction ─────────────────────────────────

"""Find the index of the free enzyme form (all sites empty)."""
function _free_enzyme_index(forms::Vector{EnzymeFormSpec})
    findfirst(f -> all(s.atoms === nothing for s in f.sites), forms)
end

"""
    _find_form(forms, occupied_subs, occupied_prods, residual_subs) → index or nothing

Find form index where core sites match the given occupancy pattern.
Regulatory and extra sites must be empty.
"""
function _find_form(
    forms::Vector{EnzymeFormSpec},
    occupied_subs::Set{Symbol},
    occupied_prods::Set{Symbol},
    residual_subs::Dict{Symbol,Vector{Pair{Symbol,Int}}},
)
    function _matches_pattern(f)
        for s in f.sites
            if s.role in (:sub, :prod) && s.index == 1
                expected = _expected_atoms(
                    s, occupied_subs, occupied_prods, residual_subs,
                )
                s.atoms != expected && return false
            elseif (s.role == :reg || s.index > 1)
                s.atoms !== nothing && return false
            end
        end
        return true
    end

    findfirst(_matches_pattern, forms)
end

"""Expected atom content for a core site given occupancy pattern."""
function _expected_atoms(s::SiteState, occupied_subs, occupied_prods, residual_subs)
    if s.role == :sub
        if haskey(residual_subs, s.metabolite)
            return residual_subs[s.metabolite]
        elseif s.metabolite in occupied_subs
            return s.full_atoms
        end
    elseif s.role == :prod
        if s.metabolite in occupied_prods
            return s.full_atoms
        end
    end
    return nothing
end

"""
    _enumerate_catalytic_cycles(forms, reaction)

Construct catalytic cycles directly by permuting substrate binding and product
release orders. Returns Vector{Vector{Tuple{Int,Int}}} where each cycle is
a sequence of (from_form_idx, to_form_idx) edges.
"""
function _enumerate_catalytic_cycles(forms::Vector{EnzymeFormSpec},
                                      @nospecialize(reaction::EnzymeReaction))
    sub_names = [s[1] for s in substrates(reaction)]
    prod_names = [p[1] for p in products(reaction)]
    free_idx = _free_enzyme_index(forms)
    free_idx === nothing && return Vector{Vector{Tuple{Int,Int}}}()

    cycles = Vector{Tuple{Int,Int}}[]

    # Standard cycles: permutation of substrate binding × permutation of product release
    for sub_perm in _permutations(sub_names)
        for prod_perm in _permutations(prod_names)
            cycle = _build_standard_cycle(forms, free_idx, sub_perm, prod_perm)
            cycle !== nothing && _add_unique_cycle!(cycles, cycle)
        end
    end

    # Ping-pong cycles: interleaved substrate binding and product release
    _enumerate_pingpong_cycles!(cycles, forms, free_idx, sub_names, prod_names)

    cycles
end

"""Build a standard (non-ping-pong) catalytic cycle."""
function _build_standard_cycle(forms::Vector{EnzymeFormSpec}, free_idx::Int,
                                sub_perm::Vector{Symbol}, prod_perm::Vector{Symbol})
    sequence = [free_idx]
    occupied_subs = Set{Symbol}()
    empty_residual = Dict{Symbol,Vector{Pair{Symbol,Int}}}()

    # Bind substrates one at a time
    for sub in sub_perm
        push!(occupied_subs, sub)
        fi = _find_form(forms, occupied_subs, Set{Symbol}(), empty_residual)
        fi === nothing && return nothing
        push!(sequence, fi)
    end

    # Isomerize: all subs → all prods
    all_prods = Set{Symbol}(prod_perm)
    fi = _find_form(forms, Set{Symbol}(), all_prods, empty_residual)
    fi === nothing && return nothing

    # Verify isomerization is valid
    ec, _, _ = edge_class(forms[sequence[end]], forms[fi])
    ec isa MustExist || return nothing
    push!(sequence, fi)

    # Release products one at a time
    remaining_prods = Set{Symbol}(prod_perm)
    for prod in prod_perm
        delete!(remaining_prods, prod)
        if isempty(remaining_prods)
            # Last product release returns to free enzyme
            fi = free_idx
        else
            fi = _find_form(forms, Set{Symbol}(), remaining_prods, empty_residual)
            fi === nothing && return nothing
        end
        push!(sequence, fi)
    end

    # Convert sequence to edges
    edges = Tuple{Int,Int}[(sequence[i], sequence[i+1]) for i in 1:length(sequence)-1]
    edges
end

"""Enumerate ping-pong catalytic cycles with interleaved binding/release."""
function _enumerate_pingpong_cycles!(cycles::Vector{Vector{Tuple{Int,Int}}},
                                      forms::Vector{EnzymeFormSpec},
                                      free_idx::Int,
                                      sub_names::Vector{Symbol},
                                      prod_names::Vector{Symbol})
    # For ping-pong, we need residual states. Check if any exist.
    has_residuals = any(forms) do f
        any(f.sites) do s
            s.role == :sub && s.index == 1 &&
                site_occupancy(s) isa ResidualSite
        end
    end
    !has_residuals && return

    # Enumerate all valid interleavings of bind-substrate and release-product
    # where products are released from substrate sites (creating residuals)
    _pingpong_dfs!(cycles, forms, free_idx, sub_names, prod_names,
                   [free_idx], Set{Symbol}(), Set{Symbol}(),
                   Dict{Symbol,Vector{Pair{Symbol,Int}}}(),
                   0, 0)
end

"""DFS to enumerate ping-pong interleavings.

`n_subs_bound` and `n_prods_released` track totals across the entire history
(unlike `bound_subs` which is cleared on isomerization).
"""
function _pingpong_dfs!(cycles, forms, free_idx, all_subs, all_prods,
                         sequence, bound_subs, released_prods, residual_state,
                         n_subs_bound::Int, n_prods_released::Int)
    # If all substrates bound and all products released, we should be back at free enzyme
    if n_subs_bound == length(all_subs) && n_prods_released == length(all_prods)
        if sequence[end] == free_idx
            edges = Tuple{Int,Int}[
                (sequence[i], sequence[i+1])
                for i in 1:length(sequence)-1
            ]
            _add_unique_cycle!(cycles, edges)
        end
        return
    end

    # Option 1: Bind next substrate (if any remain)
    for sub in all_subs
        sub in bound_subs && continue
        push!(bound_subs, sub)
        fi = _find_form(forms, bound_subs, Set{Symbol}(), residual_state)
        if fi !== nothing
            ec, _, _ = edge_class(forms[sequence[end]], forms[fi])
            if ec isa MustExist
                push!(sequence, fi)
                _pingpong_dfs!(cycles, forms, free_idx, all_subs, all_prods,
                               sequence, bound_subs, released_prods, residual_state,
                               n_subs_bound + 1, n_prods_released)
                pop!(sequence)
            end
        end
        delete!(bound_subs, sub)
    end

    # Option 2: Isomerize (if all subs bound and no prods released yet from this batch)
    if length(bound_subs) == length(all_subs) && isempty(residual_state)
        remaining = setdiff(Set{Symbol}(all_prods), released_prods)
        if !isempty(remaining)
            empty_residual = Dict{Symbol,Vector{Pair{Symbol,Int}}}()
            fi = _find_form(
                forms, Set{Symbol}(), remaining, empty_residual,
            )
            if fi !== nothing
                ec, _, _ = edge_class(forms[sequence[end]], forms[fi])
                if ec isa MustExist
                    push!(sequence, fi)
                    old_bound = copy(bound_subs)
                    empty!(bound_subs)
                    _release_prods_dfs!(
                        cycles, forms, free_idx,
                        sequence, released_prods, remaining,
                    )
                    union!(bound_subs, old_bound)
                    pop!(sequence)
                end
            end
        end
    end

    # Option 3: Release a product from a substrate site (ping-pong)
    # This creates a residual state at the substrate site
    if !isempty(bound_subs)
        current_form = forms[sequence[end]]
        for prod in all_prods
            prod in released_prods && continue
            for s in current_form.sites
                s.role == :sub && s.index == 1 && s.atoms !== nothing || continue
                # Find product atoms for this product
                prod_atoms = nothing
                for ps in current_form.sites
                    if ps.role == :prod && ps.index == 1 &&
                            ps.metabolite == prod
                        prod_atoms = ps.full_atoms
                        break
                    end
                end
                prod_atoms === nothing && continue
                # Compute residual by subtracting product atoms from substrate atoms
                sa = Dict{Symbol,Int}(a => c for (a, c) in s.atoms)
                pa = Dict{Symbol,Int}(a => c for (a, c) in prod_atoms)
                all(get(sa, a, 0) >= c for (a, c) in pa) || continue
                res = Dict{Symbol,Int}(
                    a => c - get(pa, a, 0)
                    for (a, c) in sa if c > get(pa, a, 0)
                )
                isempty(res) && continue
                res_vec = sort([a => c for (a, c) in res]; by=first)
                res_vec == s.atoms && continue  # no change

                new_residual = copy(residual_state)
                new_residual[s.metabolite] = res_vec
                fi = _find_form(forms, bound_subs, Set{Symbol}(), new_residual)
                if fi !== nothing
                    ec, _, _ = edge_class(forms[sequence[end]], forms[fi])
                    if ec isa MustExist
                        push!(sequence, fi)
                        push!(released_prods, prod)
                        old_residual = copy(residual_state)
                        merge!(residual_state, new_residual)
                        _pingpong_dfs!(
                            cycles, forms, free_idx,
                            all_subs, all_prods,
                            sequence, bound_subs,
                            released_prods, residual_state,
                            n_subs_bound, n_prods_released + 1,
                        )
                        merge!(residual_state, old_residual)
                        for k in keys(new_residual)
                            haskey(old_residual, k) || delete!(residual_state, k)
                        end
                        delete!(released_prods, prod)
                        pop!(sequence)
                    end
                end
            end
        end
    end
end

"""Release remaining products after isomerization, trying all permutations."""
function _release_prods_dfs!(cycles, forms, free_idx,
                              sequence, released_prods, remaining_prods)
    if isempty(remaining_prods)
        if sequence[end] == free_idx
            edges = Tuple{Int,Int}[
                (sequence[i], sequence[i+1])
                for i in 1:length(sequence)-1
            ]
            _add_unique_cycle!(cycles, edges)
        end
        return
    end

    for prod in collect(remaining_prods)
        delete!(remaining_prods, prod)
        push!(released_prods, prod)
        if isempty(remaining_prods)
            fi = free_idx
        else
            empty_residual = Dict{Symbol,Vector{Pair{Symbol,Int}}}()
            fi = _find_form(
                forms, Set{Symbol}(), remaining_prods,
                empty_residual,
            )
        end
        if fi !== nothing
            ec, _, _ = edge_class(forms[sequence[end]], forms[fi])
            if ec isa MustExist
                push!(sequence, fi)
                _release_prods_dfs!(cycles, forms, free_idx,
                                    sequence, released_prods, remaining_prods)
                pop!(sequence)
            end
        end
        push!(remaining_prods, prod)
        delete!(released_prods, prod)
    end
end

"""Add cycle if not already present (compare as sorted edge sets)."""
function _add_unique_cycle!(
    cycles::Vector{Vector{Tuple{Int,Int}}},
    cycle::Vector{Tuple{Int,Int}},
)
    key = _topo_key(cycle)
    for existing in cycles
        _topo_key(existing) == key && return
    end
    push!(cycles, cycle)
end

"""Generate all permutations of a vector."""
function _permutations(v::Vector{T}) where T
    length(v) <= 1 && return [copy(v)]
    result = Vector{T}[]
    for i in 1:length(v)
        rest = [v[j] for j in 1:length(v) if j != i]
        for perm in _permutations(rest)
            pushfirst!(perm, v[i])
            push!(result, perm)
        end
    end
    result
end

# ─── Cycle Combination ───────────────────────────────────────────

"""
    _combine_cycles(cycles, forms, max_forms) → Vector{Vector{Tuple{Int,Int}}}

Combine individual catalytic cycles into multi-cycle topologies.
Each topology is the union of edges from compatible cycles.
Two cycles are compatible if every emergent cycle in the
union has valid (1x or 0x) stoichiometry.
"""
function _combine_cycles(cycles::Vector{Vector{Tuple{Int,Int}}},
                          forms::Vector{EnzymeFormSpec},
                          max_forms::Int)
    isempty(cycles) && return Vector{Vector{Tuple{Int,Int}}}()

    # Start with individual cycles as topologies
    topologies = [copy(c) for c in cycles if _edge_form_count(c) <= max_forms]

    # BFS: try combining pairs
    queue = [(copy(t), i) for (i, t) in enumerate(topologies)]
    seen = Set{UInt64}([hash(_topo_key(t)) for t in topologies])

    while !isempty(queue)
        topo, max_ci = popfirst!(queue)
        for ci in (max_ci+1):length(cycles)
            merged = _merge_edges(topo, cycles[ci])
            _edge_form_count(merged) > max_forms && continue
            key = hash(_topo_key(merged))
            key in seen && continue
            push!(seen, key)

            # Validate: all emergent cycles must have valid stoichiometry
            _validate_topology(merged, forms) || continue

            push!(topologies, merged)
            push!(queue, (merged, ci))
        end
    end

    topologies
end

"""
    _enumerate_only_catalytic_mechanisms(forms, reaction; max_forms)

Enumerate all catalytic topologies (single and multi-cycle) and return as
`Vector{MechanismSpec}`. Each spec has all forms, all-SS steps, and no
constraints. Merges `_enumerate_catalytic_cycles` + `_combine_cycles`.
"""
function _enumerate_only_catalytic_mechanisms(
    forms::Vector{EnzymeFormSpec},
    @nospecialize(reaction::EnzymeReaction);
    max_forms::Int=3 * n_sites(reaction),
)
    cycles = _enumerate_catalytic_cycles(forms, reaction)
    topologies = _combine_cycles(cycles, forms, max_forms)
    [_build_spec_from_edges(topo, forms, reaction) for topo in topologies]
end

"""Count unique forms referenced by an edge list."""
_edge_form_count(edges::Vector{Tuple{Int,Int}}) =
    length(Set(Iterators.flatten(edges)))

"""Canonical key for a topology: sorted undirected edge set."""
function _topo_key(topo::Vector{Tuple{Int,Int}})
    sort([(min(a,b), max(a,b)) for (a,b) in topo])
end

"""Merge two edge lists, deduplicating by undirected key."""
function _merge_edges(a::Vector{Tuple{Int,Int}}, b::Vector{Tuple{Int,Int}})
    existing = Set((min(x,y), max(x,y)) for (x,y) in a)
    result = copy(a)
    for (x,y) in b
        key = (min(x,y), max(x,y))
        key in existing && continue
        push!(existing, key)
        push!(result, (x,y))
    end
    result
end

"""
Validate that all simple cycles through the free enzyme in a topology
have valid stoichiometry (1x catalytic or 0x futile, never partial).
"""
function _validate_topology(topo::Vector{Tuple{Int,Int}}, forms::Vector{EnzymeFormSpec})
    form_set = Set(Iterators.flatten(topo))
    free_idx = _free_enzyme_index(forms)
    free_idx === nothing && return false
    free_idx ∉ form_set && return false

    # Build adjacency for the subgraph
    adj = Dict{Int, Vector{Int}}()
    for (a, b) in topo
        push!(get!(adj, a, Int[]), b)
        push!(get!(adj, b, Int[]), a)
    end

    # DFS to enumerate all simple cycles through free_idx and validate each
    visited = Set{Int}(free_idx)
    path = Int[free_idx]

    function dfs(node::Int)
        for next in get(adj, node, Int[])
            next ∉ form_set && continue
            if next == free_idx && length(path) >= 3
                stoich = _path_stoichiometry(path, forms)
                _is_valid_stoich(stoich) || return false
            elseif next ∉ visited && next != free_idx
                push!(visited, next)
                push!(path, next)
                dfs(next) || return false
                pop!(path)
                delete!(visited, next)
            end
        end
        return true
    end

    dfs(free_idx)
end

"""Compute net stoichiometry along a path (closing back to first node)."""
function _path_stoichiometry(path::Vector{Int}, forms::Vector{EnzymeFormSpec})
    stoich = Dict{Symbol,Int}()
    for i in 1:length(path)
        j = i == length(path) ? 1 : i + 1
        ec, met, etype = edge_class(forms[path[i]], forms[path[j]])
        met === nothing && continue
        if etype == :binding
            stoich[met] = get(stoich, met, 0) - 1
        elseif etype == :release
            stoich[met] = get(stoich, met, 0) + 1
        end
    end
    stoich
end

"""
Check if stoichiometry is valid: all-zero (futile) or has
both consumption and production.
"""
function _is_valid_stoich(stoich::Dict{Symbol,Int})
    all(v == 0 for v in values(stoich)) && return true
    any(v < 0 for v in values(stoich)) && any(v > 0 for v in values(stoich))
end

# ─── Activator Shadow Cycles ──────────────────────────────────────

"""
    _generate_activator_configs(topo, forms, reaction) → Vector{Vector{Tuple{Int,Int}}}

For each regulator, generate activator shadow cycle variants:
  1. Absent — no activator edges
  2. Non-essential — shadow cycle added, base cycle retained
  3. Essential — shadow cycle added, shadowed base edges removed

Returns a list of topology variants including the original.
"""
function _generate_activator_configs(topo::Vector{Tuple{Int,Int}},
                                      forms::Vector{EnzymeFormSpec},
                                      @nospecialize(reaction::EnzymeReaction))
    reg_names = [r[1] for r in regulators(reaction)]
    isempty(reg_names) && return [topo]

    topo_forms = Set(Iterators.flatten(topo))

    # Each option: (edges_to_add, base_edges_to_remove)
    OptionT = Tuple{Vector{Tuple{Int,Int}}, Set{Tuple{Int,Int}}}
    per_reg_options = Vector{Vector{OptionT}}()
    for reg in reg_names
        options = OptionT[]
        # Option 1: absent
        push!(options, (Tuple{Int,Int}[], Set{Tuple{Int,Int}}()))

        # Find shadow forms: for each cycle form, find the
        # form with this regulator also bound
        shadow_pairs = Tuple{Int,Int}[]  # (base_form_idx, shadow_form_idx)
        for base_idx in sort(collect(topo_forms))
            shadow_idx = _find_shadow_form(forms, base_idx, reg)
            shadow_idx !== nothing && push!(shadow_pairs, (base_idx, shadow_idx))
        end

        if length(shadow_pairs) >= 2
            shadow_map = Dict(bi => si for (bi, si) in shadow_pairs)

            # Shadow cycle edges: mirror base edges in shadow forms
            shadow_cycle = Tuple{Int,Int}[]
            shadowed_base_edges = Set{Tuple{Int,Int}}()
            for (a, b) in topo
                sa = get(shadow_map, a, nothing)
                sb = get(shadow_map, b, nothing)
                if sa !== nothing && sb !== nothing
                    push!(shadow_cycle, (sa, sb))
                    push!(shadowed_base_edges, (a, b))
                end
            end

            # Option 2: non-essential (full shadow, keep base)
            full_shadow_edges = Tuple{Int,Int}[
                shadow_pairs; shadow_cycle
            ]
            push!(options, (full_shadow_edges, Set{Tuple{Int,Int}}()))

            # Option 3: essential (full shadow, remove shadowed base)
            push!(options, (full_shadow_edges, shadowed_base_edges))
        end
        push!(per_reg_options, options)
    end

    # Cartesian product of per-regulator options
    results = Vector{Tuple{Int,Int}}[]
    for combo in Iterators.product(per_reg_options...)
        all_remove = Set{Tuple{Int,Int}}()
        for (_, remove) in combo
            union!(all_remove, remove)
        end
        merged = [e for e in topo if e ∉ all_remove]
        for (add, _) in combo
            append!(merged, add)
        end
        push!(results, merged)
    end
    results
end

"""Find the shadow form: same as base but with regulator also bound."""
function _find_shadow_form(forms::Vector{EnzymeFormSpec}, base_idx::Int, reg::Symbol)
    base = forms[base_idx]
    for (i, f) in enumerate(forms)
        i == base_idx && continue
        # All non-target sites must match; target reg site must gain atoms
        found_reg = false
        all_ok = true
        for (sa, sb) in zip(base.sites, f.sites)
            if sa.metabolite == reg && sa.role == :reg &&
                    sa.atoms === nothing && sb.atoms !== nothing
                found_reg = true
            elseif sa.atoms != sb.atoms
                all_ok = false
                break
            end
        end
        all_ok && found_reg && return i
    end
    nothing
end

"""
    _generate_activator_configs(spec, forms, reaction) → Vector{MechanismSpec}

MechanismSpec interface: returns the input spec plus specs with activator
shadow edges. If no regulators, returns `[spec]`.
"""
function _generate_activator_configs(
    spec::MechanismSpec,
    forms::Vector{EnzymeFormSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    topo = _spec_to_edges(spec, forms)
    edge_variants = _generate_activator_configs(topo, forms, reaction)
    [_build_spec_from_edges(variant, forms, reaction) for variant in edge_variants]
end

# ─── Dead-End Enumeration ─────────────────────────────────────────

"""
    _reg_site_positions(form::EnzymeFormSpec) → Vector{Int}

Site positions that are regulatory and currently empty (available for binding).
"""
function _reg_site_positions(form::EnzymeFormSpec)
    [k for k in eachindex(form.sites)
     if form.sites[k].role == :reg && form.sites[k].atoms === nothing]
end

"""
    _find_dead_end_form(base_idx, occupied_positions, forms) → Int or nothing

Find the form matching `forms[base_idx]` at all sites except the specified
regulatory positions, which must be occupied.
"""
function _find_dead_end_form(base_idx::Int, occupied_positions,
                              forms::Vector{EnzymeFormSpec})
    base = forms[base_idx]
    findfirst(forms) do fj
        all(1:length(base.sites)) do k
            if k in occupied_positions
                fj.sites[k].atoms !== nothing
            else
                base.sites[k].atoms == fj.sites[k].atoms
            end
        end
    end
end

"""
    _enumerate_dead_end_configs(topo, forms, max_forms) → Vector{Vector{Int}}

Thermodynamic box rule: for each topology form, choose which regulators
bind. Choosing multiple regulators forces the multi-regulator form
(box closure). Cartesian product across topology forms with budget.
"""
function _enumerate_dead_end_configs(topo::Vector{Tuple{Int,Int}},
                                      forms::Vector{EnzymeFormSpec},
                                      max_forms::Int)
    topo_forms_set = Set(Iterators.flatten(topo))
    topo_forms = sort(collect(topo_forms_set))
    n_topo = length(topo_forms)
    budget = max_forms - n_topo
    budget < 0 && return Vector{Int}[]

    # Per topology form: compute all regulator-subset options
    per_form_options = Vector{Vector{Int}}[]
    for fi in topo_forms
        reg_sites = _reg_site_positions(forms[fi])
        n_reg = length(reg_sites)
        seen_options = Set{Vector{Int}}()
        # Empty option (no dead-end binding) is always valid
        options = [Int[]]
        push!(seen_options, Int[])

        for subset_mask in 1:((1 << n_reg) - 1)
            chosen = [reg_sites[k] for k in 1:n_reg
                      if (subset_mask >> (k - 1)) & 1 == 1]
            # Powerset closure: all non-empty subsets of chosen
            # regulators must have corresponding forms
            de_forms = Int[]
            valid = true
            for sub_mask in 1:((1 << length(chosen)) - 1)
                positions = [chosen[k]
                             for k in 1:length(chosen)
                             if (sub_mask >> (k - 1)) & 1 == 1]
                form_idx = _find_dead_end_form(
                    fi, positions, forms,
                )
                if form_idx === nothing
                    valid = false
                    break
                end
                form_idx in topo_forms_set && continue
                push!(de_forms, form_idx)
            end
            !valid && continue
            sort!(de_forms)
            if de_forms ∉ seen_options
                push!(seen_options, de_forms)
                push!(options, de_forms)
            end
        end
        push!(per_form_options, options)
    end

    # Cartesian product across topology forms, filtered by budget
    configs = Vector{Int}[]
    _dead_end_cartesian!(configs, per_form_options, 1, Int[], budget)
    configs
end

"""Cartesian product of dead-end options respecting max_forms budget."""
function _dead_end_cartesian!(configs, options, idx, current, budget)
    if idx > length(options)
        push!(configs, copy(current))
        return
    end
    for option in options[idx]
        length(option) > budget && continue
        append!(current, option)
        _dead_end_cartesian!(
            configs, options, idx + 1,
            current, budget - length(option),
        )
        resize!(current, length(current) - length(option))
    end
end

"""
    _enumerate_dead_end_configs(spec, forms; max_forms) → Vector{MechanismSpec}

MechanismSpec interface: returns the input spec (no dead-ends) plus specs
with dead-end binding steps appended. Always returns ≥ 1 element.
"""
function _enumerate_dead_end_configs(
    spec::MechanismSpec,
    forms::Vector{EnzymeFormSpec};
    max_forms::Int=3 * length(forms),
)
    topo = _spec_to_edges(spec, forms)
    de_configs = _enumerate_dead_end_configs(topo, forms, max_forms)
    [
        _build_spec_from_edges(
            _get_step_edges(topo, de, forms),
            forms, spec.reaction,
        )
        for de in de_configs
    ]
end

# ─── SS/RE + Constraint Enumeration ───────────────────────────────

"""
    _enumerate_ress_and_constraints(step_edges, forms)

Enumerate all valid RE/SS assignments and equivalent step constraints.
"""
function _enumerate_ress_and_constraints(n_steps::Int,
                                          step_edges::Vector{Tuple{Int,Int}},
                                          forms::Vector{EnzymeFormSpec})
    equiv_groups = _find_equivalent_groups(step_edges, forms)
    results = Tuple{Vector{Bool}, Vector{ParamConstraint}}[]

    # Iterate all non-zero bitmasks (at least one step must be SS, i.e. not all RE)
    for re_mask in 0:((1 << n_steps) - 2)
        eq_steps = Bool[(re_mask >> (i-1)) & 1 == 1 for i in 1:n_steps]

        # Check equivalent groups have consistent RE/SS
        valid_groups = Vector{Int}[]
        for group in equiv_groups
            first_re = eq_steps[group[1]]
            if all(eq_steps[s] == first_re for s in group)
                push!(valid_groups, group)
            end
        end

        # For each valid group, enumerate constrained/unconstrained
        n_vg = length(valid_groups)
        for cmask in 0:((1 << n_vg) - 1)
            constraints = ParamConstraint[]
            for (gi, group) in enumerate(valid_groups)
                ((cmask >> (gi - 1)) & 1) == 0 && continue
                is_re = eq_steps[group[1]]
                if is_re
                    for j in 2:length(group)
                        push!(constraints, (Symbol("K$(group[j])"), 1,
                                           [(Symbol("K$(group[1])"), 1)]))
                    end
                else
                    for j in 2:length(group)
                        push!(constraints, (Symbol("k$(group[j])f"), 1,
                                           [(Symbol("k$(group[1])f"), 1)]))
                        push!(constraints, (Symbol("k$(group[j])r"), 1,
                                           [(Symbol("k$(group[1])r"), 1)]))
                    end
                end
            end
            push!(results, (eq_steps, constraints))
        end
    end

    results
end

"""
Find groups of equivalent steps. Two steps are equivalent if both are binding
edges for the same metabolite at the same site index.
"""
function _find_equivalent_groups(step_edges::Vector{Tuple{Int,Int}},
                                  forms::Vector{EnzymeFormSpec})
    binding_key = Dict{Tuple{Symbol,Int}, Vector{Int}}()
    for (i, (a, b)) in enumerate(step_edges)
        ec, met, etype = edge_class(forms[a], forms[b])
        etype == :binding || continue
        # Find which site differs
        for k in 1:length(forms[a].sites)
            if forms[a].sites[k].atoms === nothing && forms[b].sites[k].atoms !== nothing
                key = (forms[b].sites[k].metabolite, forms[b].sites[k].index)
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

"""
    _enumerate_ress_and_constraints(spec, forms) → Vector{MechanismSpec}

MechanismSpec interface: returns all valid RE/SS × constraint variants
(including the input's all-SS variant).
"""
function _enumerate_ress_and_constraints(
    spec::MechanismSpec,
    forms::Vector{EnzymeFormSpec},
)
    edges = _spec_to_edges(spec, forms)
    n = length(edges)
    configs = _enumerate_ress_and_constraints(n, edges, forms)
    [
        _build_spec_from_edges(edges, forms, spec.reaction, eq_steps, constraints)
        for (eq_steps, constraints) in configs
    ]
end

# ─── MechanismSpec Construction ───────────────────────────────────

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

"""Convert edge list (form index pairs) to reaction tuples."""
function _edges_to_reactions(
    edges::Vector{Tuple{Int,Int}},
    forms::Vector{EnzymeFormSpec},
)
    map(edges) do (a, b)
        _, met, etype = edge_class(forms[a], forms[b])
        from_name = forms[a].name
        to_name = forms[b].name
        if etype == :binding
            ([from_name, met], [to_name])
        elseif etype == :release
            ([from_name], [to_name, met])
        else  # isomerization
            ([from_name], [to_name])
        end
    end
end

"""Extract form index edge pairs from a MechanismSpec's reactions."""
function _spec_to_edges(
    spec::MechanismSpec,
    forms::Vector{EnzymeFormSpec},
)
    name_to_idx = Dict(f.name => i for (i, f) in enumerate(forms))
    edges = Tuple{Int,Int}[]
    for (lhs, rhs) in spec.reactions
        from_idx = nothing
        to_idx = nothing
        for sym in lhs
            idx = get(name_to_idx, sym, nothing)
            idx !== nothing && (from_idx = idx)
        end
        for sym in rhs
            idx = get(name_to_idx, sym, nothing)
            idx !== nothing && (to_idx = idx)
        end
        if from_idx !== nothing && to_idx !== nothing
            push!(edges, (from_idx, to_idx))
        end
    end
    edges
end

"""Build a MechanismSpec from edges, including ALL forms."""
function _build_spec_from_edges(
    edges::Vector{Tuple{Int,Int}},
    forms::Vector{EnzymeFormSpec},
    @nospecialize(reaction),
    eq_steps::Vector{Bool}=fill(false, length(edges)),
    constraints::Vector{ParamConstraint}=ParamConstraint[],
)
    all_form_names = [f.name for f in forms]
    all_form_atoms = [_form_atoms(f.sites) for f in forms]
    rxn_tuples = _edges_to_reactions(edges, forms)
    MechanismSpec(
        reaction, all_form_names, all_form_atoms,
        rxn_tuples, eq_steps, constraints,
    )
end

"""Set of enzyme form names referenced in a spec's reactions."""
function _used_forms(spec::MechanismSpec)
    form_set = Set(spec.forms)
    used = Set{Symbol}()
    for (lhs, rhs) in spec.reactions
        for sym in Iterators.flatten((lhs, rhs))
            sym in form_set && push!(used, sym)
        end
    end
    used
end

"""Count enzyme forms referenced in a spec's reactions."""
_used_form_count(spec::MechanismSpec) = length(_used_forms(spec))

"""Get canonical step edges for a topology + dead-end config.

Wires ALL valid binding edges between all forms in the mechanism
(topology ∪ dead-end). This creates thermodynamically complete graphs;
Wegscheider constraints handle dependent parameters at the rate equation
derivation stage.
"""
function _get_step_edges(topo_edges::Vector{Tuple{Int,Int}},
                          de_forms::Vector{Int},
                          forms::Vector{EnzymeFormSpec})
    topo_forms = Set(Iterators.flatten(topo_edges))
    all_forms = sort(collect(union(topo_forms, Set(de_forms))))

    seen = Set{Tuple{Int,Int}}()
    step_edges = Tuple{Int,Int}[]

    # Topology edges (directed as given — preserves cycle direction)
    for (a, b) in topo_edges
        key = (min(a, b), max(a, b))
        key in seen && continue
        push!(seen, key)
        push!(step_edges, (a, b))
    end

    # All valid binding edges between all forms in the mechanism
    for i in 1:length(all_forms), j in (i + 1):length(all_forms)
        fi, fj = all_forms[i], all_forms[j]
        key = (fi, fj)  # fi < fj since all_forms is sorted
        key in seen && continue
        ec, _, etype = edge_class(forms[fi], forms[fj])
        (ec isa CouldExist || ec isa MustExist) && etype == :binding || continue
        push!(seen, key)
        push!(step_edges, (fi, fj))
    end

    step_edges
end

# ─── Main Entry Point ────────────────────────────────────────────

"""
    enumerate_mechanisms(reaction::EnzymeReaction; max_forms = 3 * n_sites(reaction))

Enumerate all valid mechanism topologies for the given reaction.

Returns a `MechanismIterator` that lazily yields `MechanismSpec` structs.
Each spec can be converted to an `EnzymeMechanism` via `EnzymeMechanism(spec)`.

The enumeration covers:
1. Catalytic cycles (constructive permutation-based)
2. Multi-cycle topologies (cycle combination with stoichiometry validation)
3. Activator shadow cycles (regulators participating in catalysis)
4. Dead-end complexes (inhibitors, abortive complexes, product inhibition)
5. RE/SS assignments (rapid-equilibrium vs steady-state per step)
6. Equivalent step constraints (shared parameters for equivalent binding steps)
"""
function enumerate_mechanisms(@nospecialize(reaction::EnzymeReaction);
                              max_forms::Int = 3 * n_sites(reaction))
    forms = enumerate_enzyme_forms(reaction)

    # Stage 1+2: Catalytic cycles + combination
    specs = _enumerate_only_catalytic_mechanisms(forms, reaction; max_forms)

    # Stage 3: Activator shadow cycles
    specs = MechanismSpec[
        s for spec in specs
        for s in _generate_activator_configs(spec, forms, reaction)
    ]
    filter!(s -> _used_form_count(s) <= max_forms, specs)

    # Stage 4: Dead-end complexes
    specs = MechanismSpec[
        s for spec in specs
        for s in _enumerate_dead_end_configs(spec, forms; max_forms)
    ]

    # Stage 5+6: RE/SS assignments + equivalent step constraints
    specs = MechanismSpec[
        s for spec in specs
        for s in _enumerate_ress_and_constraints(spec, forms)
    ]

    MechanismIterator(specs)
end

# ─── MechanismSpec → EnzymeMechanism Conversion ──────────────────

"""
    EnzymeMechanism(spec::MechanismSpec)

Convert a lightweight `MechanismSpec` to a type-parameterized `EnzymeMechanism`.
Delegates to the standard constructor for full validation.
"""
function EnzymeMechanism(spec::MechanismSpec)
    rxn = spec.reaction
    drop_sites(t) = Tuple((n, a) for (n, a, _) in t)
    subs_t = drop_sites(substrates(rxn))
    prods_t = drop_sites(products(rxn))
    regs_t = drop_sites(regulators(rxn))

    used = _used_forms(spec)
    enzs_t = Tuple(
        (spec.forms[i], Tuple(Tuple.(spec.form_atoms[i])))
        for i in eachindex(spec.forms)
        if spec.forms[i] in used
    )
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
