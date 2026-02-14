# ─── Trait Types ──────────────────────────────────────────────────────────────

abstract type EdgeClass end
struct MustExist <: EdgeClass end
struct CouldExist <: EdgeClass end
struct Forbidden <: EdgeClass end

abstract type SiteOccupancy end
struct EmptySite <: SiteOccupancy end
struct OccupiedSite <: SiteOccupancy end
struct ResidualSite <: SiteOccupancy end

# ─── Enriched SiteState ──────────────────────────────────────────────────────

"""
    SiteState

State of a single binding site on the enzyme.

- `metabolite`: which metabolite this site is for
- `index`: site number (1, 2, ...)
- `atoms`: atoms present (`nothing` = unoccupied, empty vector = occupied with no-atom species)
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

# ─── Core Helpers ─────────────────────────────────────────────────────────────

"""
    site_occupancy(site::SiteState) → SiteOccupancy

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
    edge_class(form_a, form_b) → (EdgeClass, metabolite, edge_type)

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
            if (oa isa EmptySite && ob isa OccupiedSite)
                return (CouldExist(), sa.metabolite, :binding)
            elseif (oa isa OccupiedSite && ob isa EmptySite)
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

"""Find the metabolite released when atoms decrease from `more` to `fewer` at a site."""
function _released_metabolite(more::Vector{Pair{Symbol,Int}}, fewer::Vector{Pair{Symbol,Int}},
                               form::EnzymeFormSpec)
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
        s.role == :prod && s.index == 1 && s.full_atoms == diff_sorted && return s.metabolite
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
    # One form must have all sub core sites occupied + all prod core sites empty,
    # and the other must have the reverse
    a_sub_occ, a_prod_empty = true, true
    b_sub_occ, b_prod_empty = true, true
    a_sub_empty, a_prod_occ = true, true
    b_sub_empty, b_prod_occ = true, true

    for k in 1:length(fa.sites)
        s = fa.sites[k]
        s.role in (:sub, :prod) && s.index == 1 || continue
        if s.role == :sub
            fa.sites[k].atoms === nothing && (a_sub_occ = false)
            fa.sites[k].atoms !== nothing && (a_sub_empty = false)
            fb.sites[k].atoms === nothing && (b_sub_occ = false)
            fb.sites[k].atoms !== nothing && (b_sub_empty = false)
        else  # :prod
            fa.sites[k].atoms !== nothing && (a_prod_empty = false)
            fa.sites[k].atoms === nothing && (a_prod_occ = false)
            fb.sites[k].atoms !== nothing && (b_prod_empty = false)
            fb.sites[k].atoms === nothing && (b_prod_occ = false)
        end
    end

    valid = (a_sub_occ && a_prod_empty && b_sub_empty && b_prod_occ) ||
            (a_sub_empty && a_prod_occ && b_sub_occ && b_prod_empty)
    !valid && return false

    # Verify atom conservation across core sites
    atoms_a = _core_atoms(fa)
    atoms_b = _core_atoms(fb)
    atoms_a == atoms_b
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

# ─── Enzyme Form Enumeration ─────────────────────────────────────────────────

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
    for s in S;  add!(s[1], 1, fatoms(s), :sub);  end
    for p in P;  add!(p[1], 1, fatoms(p), :prod); end
    for s in S, i in 2:nsites(s);  add!(s[1], i, fatoms(s), :sub);  end
    for p in P, i in 2:nsites(p);  add!(p[1], i, fatoms(p), :prod); end
    for r in R, i in 1:nsites(r);  add!(r[1], i, fatoms(r), :reg);  end

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

# ─── n_sites ──────────────────────────────────────────────────────────────────

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

# ─── Catalytic Cycle Construction ─────────────────────────────────────────────

"""Find the index of the free enzyme form (all sites empty)."""
function _free_enzyme_index(forms::Vector{EnzymeFormSpec})
    findfirst(f -> all(s.atoms === nothing for s in f.sites), forms)
end

"""
    _find_form(forms, occupied_subs, occupied_prods, residual_subs) → index or nothing

Find form index where core sites match the given occupancy pattern.
Regulatory and extra sites must be empty.
"""
function _find_form(forms::Vector{EnzymeFormSpec}, occupied_subs::Set{Symbol},
                    occupied_prods::Set{Symbol}, residual_subs::Dict{Symbol,Vector{Pair{Symbol,Int}}})
    for (i, f) in enumerate(forms)
        match = true
        for s in f.sites
            s.role in (:sub, :prod) && s.index == 1 || continue
            if s.role == :sub
                if haskey(residual_subs, s.metabolite)
                    s.atoms != residual_subs[s.metabolite] && (match = false; break)
                elseif s.metabolite in occupied_subs
                    s.atoms != s.full_atoms && (match = false; break)
                else
                    s.atoms !== nothing && (match = false; break)
                end
            else  # :prod
                if s.metabolite in occupied_prods
                    s.atoms != s.full_atoms && (match = false; break)
                else
                    s.atoms !== nothing && (match = false; break)
                end
            end
        end
        # Regulatory and extra sites must be empty for catalytic cycle forms
        if match
            for s in f.sites
                (s.role == :reg || s.index > 1) && s.atoms !== nothing && (match = false; break)
            end
        end
        match && return i
    end
    nothing
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
        any(s -> s.role == :sub && s.index == 1 && site_occupancy(s) isa ResidualSite, f.sites)
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
            edges = Tuple{Int,Int}[(sequence[i], sequence[i+1]) for i in 1:length(sequence)-1]
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
            fi = _find_form(forms, Set{Symbol}(), remaining, Dict{Symbol,Vector{Pair{Symbol,Int}}}())
            if fi !== nothing
                ec, _, _ = edge_class(forms[sequence[end]], forms[fi])
                if ec isa MustExist
                    push!(sequence, fi)
                    old_bound = copy(bound_subs)
                    empty!(bound_subs)
                    _release_prods_dfs!(cycles, forms, free_idx,
                                        sequence, released_prods, remaining)
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
                    ps.role == :prod && ps.index == 1 && ps.metabolite == prod && (prod_atoms = ps.full_atoms; break)
                end
                prod_atoms === nothing && continue
                # Compute residual by subtracting product atoms from substrate atoms
                sa = Dict{Symbol,Int}(a => c for (a, c) in s.atoms)
                pa = Dict{Symbol,Int}(a => c for (a, c) in prod_atoms)
                all(get(sa, a, 0) >= c for (a, c) in pa) || continue
                res = Dict{Symbol,Int}(a => c - get(pa, a, 0) for (a, c) in sa if c > get(pa, a, 0))
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
                        _pingpong_dfs!(cycles, forms, free_idx, all_subs, all_prods,
                                       sequence, bound_subs, released_prods, residual_state,
                                       n_subs_bound, n_prods_released + 1)
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
            edges = Tuple{Int,Int}[(sequence[i], sequence[i+1]) for i in 1:length(sequence)-1]
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
            fi = _find_form(forms, Set{Symbol}(), remaining_prods, Dict{Symbol,Vector{Pair{Symbol,Int}}}())
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
function _add_unique_cycle!(cycles::Vector{Vector{Tuple{Int,Int}}}, cycle::Vector{Tuple{Int,Int}})
    key = sort([(min(a,b), max(a,b)) for (a,b) in cycle])
    for existing in cycles
        existing_key = sort([(min(a,b), max(a,b)) for (a,b) in existing])
        existing_key == key && return
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

# ─── Cycle Combination ────────────────────────────────────────────────────────

"""
    _combine_cycles(cycles, forms, max_forms) → Vector{Vector{Tuple{Int,Int}}}

Combine individual catalytic cycles into multi-cycle topologies.
Each topology is the union of edges from compatible cycles.
Two cycles are compatible if every emergent cycle in the union has valid (1× or 0×) stoichiometry.
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

"""Check if stoichiometry is valid: all-zero (futile) or has both consumption and production."""
function _is_valid_stoich(stoich::Dict{Symbol,Int})
    all(v == 0 for v in values(stoich)) && return true
    any(v < 0 for v in values(stoich)) && any(v > 0 for v in values(stoich))
end

# ─── Activator Shadow Cycles ─────────────────────────────────────────────────

"""
    _generate_activator_configs(topo, forms, reaction) → Vector{Vector{Tuple{Int,Int}}}

For each regulator, generate activator shadow cycle variants.
Returns a list of topology variants including the original.
"""
function _generate_activator_configs(topo::Vector{Tuple{Int,Int}},
                                      forms::Vector{EnzymeFormSpec},
                                      @nospecialize(reaction::EnzymeReaction))
    reg_names = [r[1] for r in regulators(reaction)]
    isempty(reg_names) && return [topo]

    topo_forms = Set(Iterators.flatten(topo))

    # For each regulator, compute configs: absent, full_shadow, free_only_shadow
    per_reg_options = Vector{Vector{Vector{Tuple{Int,Int}}}}()
    for reg in reg_names
        options = Vector{Tuple{Int,Int}}[]
        push!(options, Tuple{Int,Int}[])  # absent: no extra edges

        # Find shadow forms: for each cycle form, find the form with this regulator also bound
        shadow_pairs = Tuple{Int,Int}[]  # (base_form_idx, shadow_form_idx)
        for base_idx in sort(collect(topo_forms))
            shadow_idx = _find_shadow_form(forms, base_idx, reg)
            shadow_idx !== nothing && push!(shadow_pairs, (base_idx, shadow_idx))
        end

        if length(shadow_pairs) >= 2
            # Full shadow: all base forms connect to their shadows, shadow cycle mirrors base
            full_shadow_edges = Tuple{Int,Int}[]
            # Connection edges: base ↔ shadow
            for (bi, si) in shadow_pairs
                push!(full_shadow_edges, (bi, si))
            end
            # Shadow cycle edges: mirror the base cycle edges among shadow forms
            shadow_map = Dict(bi => si for (bi, si) in shadow_pairs)
            for (a, b) in topo
                sa = get(shadow_map, a, nothing)
                sb = get(shadow_map, b, nothing)
                sa !== nothing && sb !== nothing && push!(full_shadow_edges, (sa, sb))
            end
            push!(options, full_shadow_edges)

            # Free-enzyme-only connection: only connect at free enzyme
            free_idx = _free_enzyme_index(forms)
            if free_idx !== nothing && free_idx in topo_forms
                free_shadow = nothing
                for (bi, si) in shadow_pairs
                    bi == free_idx && (free_shadow = si; break)
                end
                if free_shadow !== nothing
                    free_only_edges = Tuple{Int,Int}[(free_idx, free_shadow)]
                    # Shadow cycle edges (same as above)
                    for (a, b) in topo
                        sa = get(shadow_map, a, nothing)
                        sb = get(shadow_map, b, nothing)
                        sa !== nothing && sb !== nothing && push!(free_only_edges, (sa, sb))
                    end
                    push!(options, free_only_edges)
                end
            end
        end
        push!(per_reg_options, options)
    end

    # Cartesian product of per-regulator options
    results = Vector{Tuple{Int,Int}}[]
    for combo in Iterators.product(per_reg_options...)
        merged = copy(topo)
        for option in combo
            append!(merged, option)
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
        match = true
        found_reg = false
        for (sa, sb) in zip(base.sites, f.sites)
            if sa.metabolite == reg && sa.role == :reg && sa.atoms === nothing && sb.atoms !== nothing
                found_reg = true
            elseif sa.atoms != sb.atoms
                match = false
                break
            end
        end
        match && found_reg && return i
    end
    nothing
end

# ─── Dead-End Enumeration ─────────────────────────────────────────────────────

"""
    _enumerate_dead_end_configs(topo, forms, max_forms) → Vector{Vector{Int}}

For each form in the topology, find CouldExist edges to other forms (dead-ends).
Enumerate all subsets of dead-end forms, respecting max_forms budget.
"""
function _enumerate_dead_end_configs(topo::Vector{Tuple{Int,Int}},
                                      forms::Vector{EnzymeFormSpec},
                                      max_forms::Int)
    topo_forms = Set(Iterators.flatten(topo))
    n_topo = length(topo_forms)
    budget = max_forms - n_topo
    budget < 0 && return Vector{Int}[]

    # Find all candidate dead-end forms
    candidates = Set{Int}()
    for fi in topo_forms
        for (j, fj) in enumerate(forms)
            j in topo_forms && continue
            ec, _, _ = edge_class(forms[fi], fj)
            ec isa CouldExist && push!(candidates, j)
        end
    end

    # Also check dead-ends of dead-ends (2 levels deep for regulatory chains)
    level2 = Set{Int}()
    for ci in candidates
        for (j, fj) in enumerate(forms)
            j in topo_forms && continue
            j in candidates && continue
            ec, _, _ = edge_class(forms[ci], fj)
            ec isa CouldExist && push!(level2, j)
        end
    end
    union!(candidates, level2)

    candidates_vec = sort(collect(candidates))

    # Enumerate subsets up to budget size
    configs = Vector{Int}[]
    push!(configs, Int[])  # empty = no dead-ends
    _enumerate_subsets!(configs, candidates_vec, 1, Int[], budget, topo_forms, forms)
    configs
end

"""Enumerate valid dead-end subsets: each dead-end must be reachable from the topology
or from another included dead-end via a CouldExist edge."""
function _enumerate_subsets!(configs, candidates, idx, current, budget, topo_forms, forms)
    budget <= 0 && return
    for i in idx:length(candidates)
        c = candidates[i]
        # Check that c is reachable: connected to topo or to already-included dead-end
        reachable = false
        for fi in topo_forms
            ec, _, _ = edge_class(forms[fi], forms[c])
            ec isa CouldExist && (reachable = true; break)
        end
        if !reachable
            for fi in current
                ec, _, _ = edge_class(forms[fi], forms[c])
                ec isa CouldExist && (reachable = true; break)
            end
        end
        !reachable && continue

        push!(current, c)
        push!(configs, copy(current))
        _enumerate_subsets!(configs, candidates, i + 1, current, budget - 1, topo_forms, forms)
        pop!(current)
    end
end

# ─── SS/RE + Constraint Enumeration ──────────────────────────────────────────

"""
    _enumerate_ress_and_constraints(step_edges, forms) → Vector{Tuple{Vector{Bool}, constraints}}

Enumerate all valid RE/SS assignments and equivalent step constraints.
"""
function _enumerate_ress_and_constraints(n_steps::Int,
                                          step_edges::Vector{Tuple{Int,Int}},
                                          forms::Vector{EnzymeFormSpec})
    equiv_groups = _find_equivalent_groups(step_edges, forms)
    results = Tuple{Vector{Bool}, Vector{Tuple{Symbol,Int,Vector{Tuple{Symbol,Int}}}}}[]

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
            constraints = Tuple{Symbol,Int,Vector{Tuple{Symbol,Int}}}[]
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

# ─── MechanismSpec Construction ───────────────────────────────────────────────

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

"""
Build a MechanismSpec from topology edges, dead-end forms, RE/SS assignment, and constraints.
"""
function _build_mechanism_spec(topo_edges::Vector{Tuple{Int,Int}},
                                dead_end_forms::Vector{Int},
                                eq_steps::Vector{Bool},
                                constraints::Vector{Tuple{Symbol,Int,Vector{Tuple{Symbol,Int}}}},
                                forms::Vector{EnzymeFormSpec},
                                @nospecialize(reaction))
    # Collect all form indices
    topo_form_set = Set(Iterators.flatten(topo_edges))
    all_form_indices = sort(collect(union(topo_form_set, Set(dead_end_forms))))

    # Get step edges: topo edges + dead-end connection edges
    step_edges = _get_step_edges(topo_edges, dead_end_forms, forms)

    form_names = [forms[fi].name for fi in all_form_indices]
    form_atoms_list = [_form_atoms(forms[fi].sites) for fi in all_form_indices]

    rxn_tuples = Tuple{Vector{Symbol}, Vector{Symbol}}[]
    for (a, b) in step_edges
        ec, met, etype = edge_class(forms[a], forms[b])
        from_name = forms[a].name
        to_name = forms[b].name
        if etype == :binding
            lhs = met === nothing ? [from_name] : [from_name, met]
            rhs = [to_name]
        elseif etype == :release
            lhs = [from_name]
            rhs = met === nothing ? [to_name] : [to_name, met]
        else  # isomerization
            lhs = [from_name]
            rhs = [to_name]
        end
        push!(rxn_tuples, (lhs, rhs))
    end

    MechanismSpec(reaction, form_names, form_atoms_list, rxn_tuples, eq_steps, constraints)
end

"""Get canonical step edges for a topology + dead-end config."""
function _get_step_edges(topo_edges::Vector{Tuple{Int,Int}},
                          de_forms::Vector{Int},
                          forms::Vector{EnzymeFormSpec})
    topo_forms = Set(Iterators.flatten(topo_edges))
    all_forms_set = union(topo_forms, Set(de_forms))

    seen = Set{Tuple{Int,Int}}()
    step_edges = Tuple{Int,Int}[]

    # Topology edges (directed as given — preserves cycle direction)
    for (a, b) in topo_edges
        key = (min(a, b), max(a, b))
        key in seen && continue
        push!(seen, key)
        push!(step_edges, (a, b))
    end

    # Dead-end edges: find binding direction from mechanism to dead-end
    for de_idx in de_forms
        for fi in sort(collect(all_forms_set))
            fi == de_idx && continue
            ec, _, etype = edge_class(forms[fi], forms[de_idx])
            (ec isa CouldExist || ec isa MustExist) && etype == :binding || continue
            key = (min(fi, de_idx), max(fi, de_idx))
            key in seen && continue
            push!(seen, key)
            push!(step_edges, (fi, de_idx))
            break
        end
    end

    step_edges
end

# ─── Main Entry Point ────────────────────────────────────────────────────────

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

    # Stage 1: Catalytic cycles
    cycles = _enumerate_catalytic_cycles(forms, reaction)

    # Stage 2: Multi-cycle topologies
    topologies = _combine_cycles(cycles, forms, max_forms)

    # Stages 3-6: For each topology, generate activator configs, dead-ends, RE/SS, constraints
    specs = MechanismSpec[]
    for topo in topologies
        # Stage 3: Activator shadow cycles
        activated_topos = _generate_activator_configs(topo, forms, reaction)

        for act_topo in activated_topos
            _edge_form_count(act_topo) > max_forms && continue

            # Stage 4: Dead-end enumeration
            de_configs = _enumerate_dead_end_configs(act_topo, forms, max_forms)

            for de_config in de_configs
                # Get step edges for this configuration
                step_edges = _get_step_edges(act_topo, de_config, forms)
                n_steps = length(step_edges)

                # Stage 5-6: RE/SS + constraints
                ress_configs = _enumerate_ress_and_constraints(n_steps, step_edges, forms)

                for (eq_steps, constraints) in ress_configs
                    spec = _build_mechanism_spec(act_topo, de_config, eq_steps,
                                                  constraints, forms, reaction)
                    push!(specs, spec)
                end
            end
        end
    end

    MechanismIterator(specs)
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
    enzs_t = Tuple((spec.forms[i], Tuple(Tuple.(spec.form_atoms[i]))) for i in eachindex(spec.forms))
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
