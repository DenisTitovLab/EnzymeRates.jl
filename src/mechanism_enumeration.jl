# ─── Data Types ─────────────────────────────────────────────

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
    PreRessEntry

Pre-computed data for one dead-end topology, used by the lazy iterator
to generate RE/SS + constraint variants on demand.
"""
struct PreRessEntry
    rxn_tuples::Vector{Tuple{Vector{Symbol}, Vector{Symbol}}}
    edges::Vector{Tuple{Int,Int}}
    equiv_groups::Vector{Vector{Int}}
    ress_count::Int
end

"""
    MechanismIterator

Lazy iterator over all valid mechanisms for a reaction.
Stages 1-3 (catalytic, activator, dead-end) are eagerly materialized.
Stage 4 (RE/SS + constraints) is generated lazily during iteration.
"""
struct MechanismIterator
    reaction::Any
    form_names::Vector{Symbol}
    form_atoms::Vector{Vector{Pair{Symbol,Int}}}
    pre_ress::Vector{PreRessEntry}
    total::Int
end

Base.eltype(::Type{MechanismIterator}) = MechanismSpec
Base.IteratorSize(::Type{MechanismIterator}) = Base.HasLength()
Base.length(iter::MechanismIterator) = iter.total

function Base.iterate(iter::MechanismIterator)
    isempty(iter.pre_ress) && return nothing
    _produce_ress_spec(iter, 1, 0, 0)
end

function Base.iterate(iter::MechanismIterator, state)
    state === nothing && return nothing
    _produce_ress_spec(iter, state...)
end

"""Build the MechanismSpec at state (entry_idx, re_mask, cmask) and advance."""
function _produce_ress_spec(iter::MechanismIterator,
                             entry_idx::Int, re_mask::Int, cmask::Int)
    entry = iter.pre_ress[entry_idx]
    n = length(entry.edges)
    eq_steps = Bool[(re_mask >> (i - 1)) & 1 == 1 for i in 1:n]

    # Build constraints from valid equiv groups and cmask
    constraints = ParamConstraint[]
    n_valid = 0
    for group in entry.equiv_groups
        first_re = eq_steps[group[1]]
        all(eq_steps[s] == first_re for s in group) || continue
        n_valid += 1
        ((cmask >> (n_valid - 1)) & 1) == 0 && continue
        if first_re
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

    spec = MechanismSpec(
        iter.reaction, iter.form_names, iter.form_atoms,
        entry.rxn_tuples, eq_steps, constraints,
    )
    next = _advance_ress_state(iter, entry_idx, re_mask, cmask, n_valid)
    return (spec, next)
end

"""Advance to the next valid (entry_idx, re_mask, cmask) state."""
function _advance_ress_state(iter::MechanismIterator, entry_idx::Int,
                              re_mask::Int, cmask::Int, n_valid::Int)
    entry = iter.pre_ress[entry_idx]
    n = length(entry.edges)

    # Next cmask within same re_mask
    if cmask < (1 << n_valid) - 1
        return (entry_idx, re_mask, cmask + 1)
    end

    # Next re_mask within same entry
    next_re = re_mask + 1
    if next_re <= (1 << n) - 2
        return (entry_idx, next_re, 0)
    end

    # Next entry
    next_entry = entry_idx + 1
    next_entry > length(iter.pre_ress) && return nothing
    return (next_entry, 0, 0)
end

"""Count total RE/SS + constraint variants for given edges and equiv groups."""
function _count_ress_variants(n_edges::Int,
                               equiv_groups::Vector{Vector{Int}})
    n_edges == 0 && return 0
    count = 0
    for re_mask in 0:((1 << n_edges) - 2)
        n_valid = 0
        for group in equiv_groups
            first_re = (re_mask >> (group[1] - 1)) & 1
            consistent = all(
                ((re_mask >> (s - 1)) & 1) == first_re for s in group)
            consistent && (n_valid += 1)
        end
        count += 1 << n_valid
    end
    count
end

# ─── Core Helpers ──────────────────────────────────────────────

"""
    edge_class(form_a, form_b) → (is_valid, metabolite, edge_type)

Classify the edge between two enzyme forms based on site occupancy diffs.

Returns `(false, nothing, :none)` if no valid edge exists.
Returns `(true, metabolite, edge_type)` otherwise, where edge_type is
`:binding`, `:release`, or `:isomerization`.
"""
function edge_class(fa::EnzymeFormSpec, fb::EnzymeFormSpec)
    nsites = length(fa.sites)
    diff_positions = Int[]
    for k in 1:nsites
        if fa.sites[k].atoms != fb.sites[k].atoms
            push!(diff_positions, k)
        end
    end

    isempty(diff_positions) && return (false, nothing, :none)

    if length(diff_positions) == 1
        k = diff_positions[1]
        sa, sb = fa.sites[k], fb.sites[k]
        is_core = sa.role in (:sub, :prod) && sa.index == 1

        if is_core
            a_empty = sa.atoms === nothing
            b_empty = sb.atoms === nothing
            a_full = sa.atoms == sa.full_atoms
            b_full = sb.atoms == sb.full_atoms
            a_residual = !a_empty && !a_full
            b_residual = !b_empty && !b_full

            # Empty↔Occupied: binding/release of metabolite
            if a_empty && b_full
                return (true, sa.metabolite, :binding)
            elseif a_full && b_empty
                return (true, sa.metabolite, :release)
            # Residual↔Empty: release of remaining atoms
            elseif a_residual && b_empty
                met = _residual_metabolite(sa.atoms, fa)
                met === nothing && return (false, nothing, :none)
                return (true, met, :release)
            else
                return (false, nothing, :none)
            end
        else
            # Regulatory or extra site
            if sa.atoms === nothing && sb.atoms !== nothing
                return (true, sa.metabolite, :binding)
            elseif sa.atoms !== nothing && sb.atoms === nothing
                return (true, sa.metabolite, :release)
            else
                return (false, nothing, :none)
            end
        end
    end

    # ≥2 diffs: check for isomerization (all diffs at core sites)
    all_core = all(diff_positions) do k
        fa.sites[k].role in (:sub, :prod) && fa.sites[k].index == 1
    end
    !all_core && return (false, nothing, :none)

    if _is_valid_isomerization(fa, fb)
        return (true, nothing, :isomerization)
    end
    return (false, nothing, :none)
end

"""Find the metabolite whose atoms match a residual content exactly."""
function _residual_metabolite(atoms::Vector{Pair{Symbol,Int}}, form::EnzymeFormSpec)
    for s in form.sites
        s.role == :prod && s.index == 1 && s.full_atoms == atoms && return s.metabolite
    end
    nothing
end

"""
Check if two forms represent a valid isomerization.

Standard case: ALL core sites differ (all-subs ↔ all-prods).
Non-differing substrate sites must be empty in the standard case.

Ping-pong case: a substrate site has a residual atom content (partial
transformation). Non-differing substrate sites may be occupied (e.g.,
B stays bound while A undergoes partial transformation). Only atom
balance is required.

Total atom balance is always verified.
"""
function _is_valid_isomerization(fa::EnzymeFormSpec, fb::EnzymeFormSpec)
    a_sub_occ = true;  a_sub_empty = true
    a_prod_occ = true; a_prod_empty = true
    b_sub_occ = true;  b_sub_empty = true
    b_prod_occ = true; b_prod_empty = true
    has_sub_diff = false; has_prod_diff = false
    has_residual = false
    non_diff_sub_occupied = false

    for k in eachindex(fa.sites)
        s = fa.sites[k]
        s.role in (:sub, :prod) && s.index == 1 || continue
        if fa.sites[k].atoms == fb.sites[k].atoms
            if s.role == :sub && fa.sites[k].atoms !== nothing
                non_diff_sub_occupied = true
            end
            continue
        end
        if s.role == :sub
            has_sub_diff = true
            a, b = fa.sites[k].atoms, fb.sites[k].atoms
            if (a !== nothing && a != s.full_atoms) ||
               (b !== nothing && b != s.full_atoms)
                has_residual = true
            end
            a === nothing && (a_sub_occ = false)
            a !== nothing && (a_sub_empty = false)
            b === nothing && (b_sub_occ = false)
            b !== nothing && (b_sub_empty = false)
        else
            has_prod_diff = true
            fa.sites[k].atoms === nothing && (a_prod_occ = false)
            fa.sites[k].atoms !== nothing && (a_prod_empty = false)
            fb.sites[k].atoms === nothing && (b_prod_occ = false)
            fb.sites[k].atoms !== nothing && (b_prod_empty = false)
        end
    end
    (!has_sub_diff || !has_prod_diff) && return false

    if has_residual
        # Ping-pong isomerization: atom balance is sufficient.
        # Non-differing occupied subs are allowed (e.g., B stays bound
        # while A undergoes partial transformation to product).
        return _core_atoms(fa) == _core_atoms(fb)
    end

    # Standard isomerization: ALL core sites must differ, and
    # non-differing substrate sites must be empty.
    non_diff_sub_occupied && return false
    for k in eachindex(fa.sites)
        s = fa.sites[k]
        s.role in (:sub, :prod) && s.index == 1 || continue
        fa.sites[k].atoms == fb.sites[k].atoms && return false
    end

    valid = (a_sub_occ && a_prod_empty &&
             b_sub_empty && b_prod_occ) ||
            (a_sub_empty && a_prod_occ &&
             b_sub_occ && b_prod_empty)
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

# ─── Enzyme Form Enumeration ───────────────────────────────────

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

# ─── n_sites ───────────────────────────────────────────────────

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

# ─── Form Set Helpers ──────────────────────────────────────────

"""Find the index of the free enzyme form (all sites empty)."""
function _free_enzyme_index(forms::Vector{EnzymeFormSpec})
    findfirst(f -> all(s.atoms === nothing for s in f.sites), forms)
end

"""
    _derive_edges(forms, form_set) → Vector{Tuple{Int,Int}}

Given a set of form indices, find all valid edges between them.
"""
function _derive_edges(forms::Vector{EnzymeFormSpec}, form_set::Set{Int})
    sorted = sort(collect(form_set))
    edges = Tuple{Int,Int}[]
    for i in 1:length(sorted), j in (i + 1):length(sorted)
        fi, fj = sorted[i], sorted[j]
        is_valid, _, _ = edge_class(forms[fi], forms[fj])
        is_valid && push!(edges, (fi, fj))
    end
    edges
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

# ─── Catalytic Form Set Enumeration ────────────────────────────

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
    function _matches(f)
        for s in f.sites
            if s.role in (:sub, :prod) && s.index == 1
                if s.role == :sub
                    if haskey(residual_subs, s.metabolite)
                        s.atoms != residual_subs[s.metabolite] && return false
                    elseif s.metabolite in occupied_subs
                        s.atoms != s.full_atoms && return false
                    else
                        s.atoms !== nothing && return false
                    end
                else  # :prod
                    if s.metabolite in occupied_prods
                        s.atoms != s.full_atoms && return false
                    else
                        s.atoms !== nothing && return false
                    end
                end
            elseif s.role == :reg || s.index > 1
                s.atoms !== nothing && return false
            end
        end
        return true
    end
    findfirst(_matches, forms)
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

"""
    _enumerate_catalytic_form_sets(forms, reaction) → Vector{Set{Int}}

Enumerate all valid catalytic form sets by constructing minimal cycles
(standard and ping-pong) then combining into multi-cycle unions.
"""
function _enumerate_catalytic_form_sets(
    forms::Vector{EnzymeFormSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    sub_names = [s[1] for s in substrates(reaction)]
    prod_names = [p[1] for p in products(reaction)]
    free_idx = _free_enzyme_index(forms)
    free_idx === nothing && return Set{Int}[]

    cycles = Set{Int}[]
    empty_residual = Dict{Symbol,Vector{Pair{Symbol,Int}}}()

    # Standard cycles: permutation of binding × permutation of release
    for sub_perm in _permutations(sub_names)
        for prod_perm in _permutations(prod_names)
            fs = _build_standard_form_set(
                forms, free_idx, sub_perm, prod_perm, empty_residual,
            )
            fs !== nothing && _add_unique!(cycles, fs)
        end
    end

    # Ping-pong cycles
    _enumerate_pingpong_form_sets!(
        cycles, forms, free_idx, sub_names, prod_names,
    )

    # Multi-cycle: enumerate all 2^n unions of individual cycles
    combined = _combine_form_sets(cycles)
    filter!(fs -> _is_pure_topology(fs, forms), combined)
    combined
end

"""Build a standard catalytic cycle as a form set."""
function _build_standard_form_set(
    forms::Vector{EnzymeFormSpec},
    free_idx::Int,
    sub_perm::Vector{Symbol},
    prod_perm::Vector{Symbol},
    residual::Dict{Symbol,Vector{Pair{Symbol,Int}}},
)
    form_set = Set{Int}(free_idx)
    occupied_subs = Set{Symbol}()

    # Bind substrates one at a time
    all_subs_fi = free_idx
    for sub in sub_perm
        push!(occupied_subs, sub)
        fi = _find_form(forms, occupied_subs, Set{Symbol}(), residual)
        fi === nothing && return nothing
        push!(form_set, fi)
        all_subs_fi = fi
    end

    # Isomerize: all subs → all prods
    all_prods = Set{Symbol}(prod_perm)
    fi = _find_form(forms, Set{Symbol}(), all_prods, residual)
    fi === nothing && return nothing
    is_valid, _, _ = edge_class(forms[all_subs_fi], forms[fi])
    is_valid || return nothing
    push!(form_set, fi)

    # Release products one at a time
    remaining = Set{Symbol}(prod_perm)
    for prod in prod_perm
        delete!(remaining, prod)
        if isempty(remaining)
            push!(form_set, free_idx)
        else
            fi = _find_form(forms, Set{Symbol}(), remaining, residual)
            fi === nothing && return nothing
            push!(form_set, fi)
        end
    end

    form_set
end

"""Enumerate ping-pong catalytic form sets."""
function _enumerate_pingpong_form_sets!(
    cycles::Vector{Set{Int}},
    forms::Vector{EnzymeFormSpec},
    free_idx::Int,
    sub_names::Vector{Symbol},
    prod_names::Vector{Symbol},
)
    # Check if any residual forms exist
    has_residuals = any(forms) do f
        any(f.sites) do s
            s.role == :sub && s.index == 1 &&
                s.atoms !== nothing && s.atoms != s.full_atoms
        end
    end
    !has_residuals && return

    # Use sequence-based DFS (like old code) to correctly track paths,
    # then convert completed sequences to form sets
    _pingpong_dfs!(cycles, forms, free_idx, sub_names, prod_names,
                   [free_idx], Set{Symbol}(), Set{Symbol}(),
                   Dict{Symbol,Vector{Pair{Symbol,Int}}}(),
                   0, 0)
end

"""DFS to enumerate ping-pong interleavings.

Tracks form sequences (not sets) for correct backtracking.
Completed cycles are converted to form sets and deduplicated.
"""
function _pingpong_dfs!(cycles, forms, free_idx, all_subs, all_prods,
                         sequence::Vector{Int},
                         bound_subs, released_prods, residual_state,
                         n_subs_bound::Int, n_prods_released::Int)
    if n_subs_bound == length(all_subs) && n_prods_released == length(all_prods)
        if sequence[end] == free_idx
            _add_unique!(cycles, Set{Int}(sequence))
        end
        return
    end

    # Option 1: Bind next substrate
    for sub in all_subs
        sub in bound_subs && continue
        push!(bound_subs, sub)
        fi = _find_form(forms, bound_subs, Set{Symbol}(), residual_state)
        if fi !== nothing
            is_valid, _, _ = edge_class(forms[sequence[end]], forms[fi])
            if is_valid
                push!(sequence, fi)
                _pingpong_dfs!(cycles, forms, free_idx, all_subs, all_prods,
                               sequence, bound_subs, released_prods,
                               residual_state,
                               n_subs_bound + 1, n_prods_released)
                pop!(sequence)
            end
        end
        delete!(bound_subs, sub)
    end

    # Option 2: Isomerize (all subs bound → remaining prods)
    # Require a new substrate binding since the last ping-pong release:
    # n_subs_bound > n_prods_released ensures this (blocks e.g. Uni-Bi
    # where A[C2]→P1[C]+P2[C] has only 1 substrate, so after releasing
    # P1 from the A-site no second substrate can bind before isomerization)
    if length(bound_subs) == length(all_subs) && n_subs_bound > n_prods_released
        remaining = setdiff(Set{Symbol}(all_prods), released_prods)
        if !isempty(remaining)
            empty_residual = Dict{Symbol,Vector{Pair{Symbol,Int}}}()
            fi = _find_form(forms, Set{Symbol}(), remaining, empty_residual)
            if fi !== nothing
                is_valid, _, _ = edge_class(forms[sequence[end]], forms[fi])
                if is_valid
                    push!(sequence, fi)
                    old_bound = copy(bound_subs)
                    empty!(bound_subs)
                    _release_prods_dfs!(
                        cycles, forms, free_idx, sequence,
                        released_prods, remaining,
                    )
                    union!(bound_subs, old_bound)
                    pop!(sequence)
                end
            end
        end
    end

    # Option 3: Ping-pong isomerization + product release (two steps)
    # Step 1: isomerize (product appears on enzyme, substrate → residual)
    # Step 2: release product from enzyme
    if !isempty(bound_subs)
        current_form = forms[sequence[end]]
        for prod in all_prods
            prod in released_prods && continue
            for s in current_form.sites
                s.role == :sub && s.index == 1 && s.atoms !== nothing || continue
                prod_atoms = nothing
                for ps in current_form.sites
                    if ps.role == :prod && ps.index == 1 && ps.metabolite == prod
                        prod_atoms = ps.full_atoms
                        break
                    end
                end
                prod_atoms === nothing && continue
                sa = Dict{Symbol,Int}(a => c for (a, c) in s.atoms)
                pa = Dict{Symbol,Int}(a => c for (a, c) in prod_atoms)
                all(get(sa, a, 0) >= c for (a, c) in pa) || continue
                res = Dict{Symbol,Int}(
                    a => c - get(pa, a, 0)
                    for (a, c) in sa if c > get(pa, a, 0)
                )
                isempty(res) && continue
                res_vec = sort([a => c for (a, c) in res]; by=first)
                res_vec == s.atoms && continue

                new_residual = copy(residual_state)
                new_residual[s.metabolite] = res_vec

                # Step 1: Find intermediate form (product on enzyme +
                # residual on substrate site)
                inter_fi = _find_form(
                    forms, bound_subs, Set{Symbol}([prod]), new_residual,
                )
                inter_fi === nothing && continue
                is_isom, _, _ = edge_class(
                    forms[sequence[end]], forms[inter_fi],
                )
                is_isom || continue

                # Step 2: Find post-release form (product gone)
                post_fi = _find_form(
                    forms, bound_subs, Set{Symbol}(), new_residual,
                )
                post_fi === nothing && continue
                is_rel, _, _ = edge_class(
                    forms[inter_fi], forms[post_fi],
                )
                is_rel || continue

                push!(sequence, inter_fi)
                push!(sequence, post_fi)
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
                pop!(sequence)
            end
        end
    end
end

"""Release remaining products after isomerization, trying all permutations."""
function _release_prods_dfs!(cycles, forms, free_idx,
                              sequence::Vector{Int},
                              released_prods, remaining_prods)
    if isempty(remaining_prods)
        if sequence[end] == free_idx
            _add_unique!(cycles, Set{Int}(sequence))
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
            fi = _find_form(forms, Set{Symbol}(), remaining_prods, empty_residual)
        end
        if fi !== nothing
            is_valid, _, _ = edge_class(forms[sequence[end]], forms[fi])
            if is_valid
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

"""Add form set if not already present."""
function _add_unique!(sets::Vector{Set{Int}}, new_set::Set{Int})
    for existing in sets
        existing == new_set && return
    end
    push!(sets, copy(new_set))
end

"""Combine individual form set cycles into multi-cycle unions."""
function _combine_form_sets(cycles::Vector{Set{Int}})
    isempty(cycles) && return Set{Int}[]

    result = copy(cycles)
    seen = Set{Set{Int}}(cycles)

    # BFS: try combining pairs
    queue = [(c, i) for (i, c) in enumerate(cycles)]
    while !isempty(queue)
        fs, max_ci = popfirst!(queue)
        for ci in (max_ci + 1):length(cycles)
            merged = union(fs, cycles[ci])
            merged in seen && continue
            push!(seen, merged)
            push!(result, merged)
            push!(queue, (merged, ci))
        end
    end

    result
end

"""
    _is_pure_topology(form_set, forms) → Bool

Check whether a catalytic topology is pure sequential or pure ping-pong.

A topology is accepted if EITHER:
- **Pure sequential**: no residual forms at all (standard ternary-complex path)
- **Pure ping-pong**: has a free enzyme intermediate (only residuals, all other
  core sites empty) AND does NOT have the all-substrates-fully-bound form

Mixed topologies (combining both types, or having residuals without a free
intermediate) are rejected as biochemically implausible.
"""
function _is_pure_topology(form_set::Set{Int}, forms::Vector{EnzymeFormSpec})
    has_residual = false
    has_free_intermediate = false
    has_all_subs_full = false

    for fi in form_set
        f = forms[fi]
        form_has_residual = false
        form_is_free_intermediate = true
        form_all_subs_full = true

        for s in f.sites
            s.index == 1 || continue
            s.role == :reg && continue

            if s.role == :sub
                if s.atoms !== nothing && s.atoms != s.full_atoms
                    form_has_residual = true
                    form_all_subs_full = false
                elseif s.atoms === nothing
                    form_all_subs_full = false
                else  # s.atoms == s.full_atoms
                    form_is_free_intermediate = false
                end
            elseif s.role == :prod
                if s.atoms !== nothing
                    form_is_free_intermediate = false
                end
            end
        end

        if form_has_residual
            has_residual = true
            form_is_free_intermediate && (has_free_intermediate = true)
        end

        form_all_subs_full && (has_all_subs_full = true)
    end

    !has_residual && return true
    has_free_intermediate && !has_all_subs_full && return true
    return false
end

# ─── Catalytic Mechanism Construction ─────────────────────────

"""
    _enumerate_only_catalytic_mechanisms(forms, reaction; max_forms)

Enumerate all catalytic topologies (single and multi-cycle) and return as
`Vector{MechanismSpec}`. Each spec has all-SS steps and no constraints.
"""
function _enumerate_only_catalytic_mechanisms(
    forms::Vector{EnzymeFormSpec},
    @nospecialize(reaction::EnzymeReaction);
    max_forms::Int=3 * n_sites(reaction),
    shared_names::Union{Nothing,Vector{Symbol}}=nothing,
    shared_atoms::Union{Nothing,Vector{Vector{Pair{Symbol,Int}}}}=nothing,
)
    form_sets = _enumerate_catalytic_form_sets(forms, reaction)
    filter!(fs -> length(fs) <= max_forms, form_sets)
    [_build_spec_from_edges(
        _derive_edges(forms, fs), forms, reaction;
        shared_names, shared_atoms)
     for fs in form_sets]
end

"""Build a MechanismSpec from explicit edges, including ALL forms."""
function _build_spec_from_edges(
    edges::Vector{Tuple{Int,Int}},
    forms::Vector{EnzymeFormSpec},
    @nospecialize(reaction),
    eq_steps::Vector{Bool}=fill(false, length(edges)),
    constraints::Vector{ParamConstraint}=ParamConstraint[];
    shared_names::Union{Nothing,Vector{Symbol}}=nothing,
    shared_atoms::Union{Nothing,Vector{Vector{Pair{Symbol,Int}}}}=nothing,
)
    fn = shared_names !== nothing ? shared_names :
        [f.name for f in forms]
    fa = shared_atoms !== nothing ? shared_atoms :
        [_form_atoms(f.sites) for f in forms]
    rxn_tuples = _edges_to_reactions(edges, forms)
    MechanismSpec(reaction, fn, fa, rxn_tuples, eq_steps, constraints)
end

# ─── Activator Variants ───────────────────────────────────────

"""
    _generate_activator_configs(spec, forms, reaction) → Vector{MechanismSpec}

For each regulator, generate activator shadow variants:
  1. Absent — no activator edges
  2. Non-essential — shadow cycle added, base cycle retained
  3. Essential — shadow cycle added, shadowed base edges removed

Returns specs including the original.
"""
function _generate_activator_configs(
    spec::MechanismSpec,
    forms::Vector{EnzymeFormSpec},
    @nospecialize(reaction::EnzymeReaction);
    shared_names::Union{Nothing,Vector{Symbol}}=nothing,
    shared_atoms::Union{Nothing,Vector{Vector{Pair{Symbol,Int}}}}=nothing,
)
    reg_names = [r[1] for r in regulators(reaction)]
    isempty(reg_names) && return [spec]

    topo = _spec_to_edges(spec, forms)
    topo_forms = Set(Iterators.flatten(topo))

    # Each option: (edges_to_add, base_edges_to_remove)
    OptionT = Tuple{Vector{Tuple{Int,Int}}, Set{Tuple{Int,Int}}}
    per_reg_options = Vector{Vector{OptionT}}()

    for reg in reg_names
        options = OptionT[]
        # Option 1: absent
        push!(options, (Tuple{Int,Int}[], Set{Tuple{Int,Int}}()))

        # Find shadow forms: for each cycle form, the form with reg also bound
        shadow_pairs = Tuple{Int,Int}[]
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

            # Full shadow edges: binding between base↔shadow + shadow cycle
            full_shadow_edges = Tuple{Int,Int}[shadow_pairs; shadow_cycle]

            # Option 2: non-essential (full shadow, keep base)
            push!(options, (full_shadow_edges, Set{Tuple{Int,Int}}()))

            # Option 3: essential (entry binding + shadow cycle, remove base)
            # Only the bare enzyme connects to the shadow cycle; other base
            # forms (ES, EP, …) are not part of the essential topology.
            entry_idx = findfirst(shadow_pairs) do (bi, _)
                all(s -> s.role == :reg || s.atoms === nothing,
                    forms[bi].sites)
            end
            if entry_idx !== nothing
                essential_edges = Tuple{Int,Int}[
                    shadow_pairs[entry_idx]; shadow_cycle]
                push!(options, (essential_edges, shadowed_base_edges))
            end
        end
        push!(per_reg_options, options)
    end

    # Cartesian product of per-regulator options
    results = MechanismSpec[]
    for combo in Iterators.product(per_reg_options...)
        all_remove = Set{Tuple{Int,Int}}()
        for (_, remove) in combo
            union!(all_remove, remove)
        end
        merged = [e for e in topo if e ∉ all_remove]
        for (add, _) in combo
            append!(merged, add)
        end
        push!(results, _build_spec_from_edges(
            merged, forms, reaction;
            shared_names, shared_atoms))
    end
    results
end

"""Find the shadow form: same as base but with regulator also bound."""
function _find_shadow_form(forms::Vector{EnzymeFormSpec}, base_idx::Int, reg::Symbol)
    base = forms[base_idx]
    for (i, f) in enumerate(forms)
        i == base_idx && continue
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

# ─── Dead-End Enumeration ─────────────────────────────────────

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
    _enumerate_dead_end_configs(spec, forms; max_forms) → Vector{MechanismSpec}

Enumerate dead-end binding configurations for a mechanism spec.
Thermodynamic box rule: choosing multiple regulators forces multi-regulator forms.
"""
function _enumerate_dead_end_configs(
    spec::MechanismSpec,
    forms::Vector{EnzymeFormSpec};
    max_forms::Int=3 * length(forms),
    shared_names::Union{Nothing,Vector{Symbol}}=nothing,
    shared_atoms::Union{Nothing,Vector{Vector{Pair{Symbol,Int}}}}=nothing,
)
    topo = _spec_to_edges(spec, forms)
    topo_forms_set = Set(Iterators.flatten(topo))
    topo_forms = sort(collect(topo_forms_set))
    n_topo = length(topo_forms)
    budget = max_forms - n_topo
    budget < 0 && return MechanismSpec[]

    # Determine activator regulator positions (occupied in any topo form).
    # A regulator is either an activator or an inhibitor, never both.
    activator_positions = Set{Int}()
    for fi in topo_forms
        for k in eachindex(forms[fi].sites)
            if forms[fi].sites[k].role == :reg &&
               forms[fi].sites[k].atoms !== nothing
                push!(activator_positions, k)
            end
        end
    end

    # Per topology form: compute all regulator-subset options
    per_form_options = Vector{Vector{Int}}[]
    for fi in topo_forms
        reg_sites = [k for k in _reg_site_positions(forms[fi])
                     if k ∉ activator_positions]
        n_reg = length(reg_sites)
        seen_options = Set{Vector{Int}}()
        options = [Int[]]
        push!(seen_options, Int[])

        for subset_mask in 1:((1 << n_reg) - 1)
            chosen = [reg_sites[k] for k in 1:n_reg
                      if (subset_mask >> (k - 1)) & 1 == 1]
            de_forms = Int[]
            valid = true
            for sub_mask in 1:((1 << length(chosen)) - 1)
                positions = [chosen[k]
                             for k in 1:length(chosen)
                             if (sub_mask >> (k - 1)) & 1 == 1]
                form_idx = _find_dead_end_form(fi, positions, forms)
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

    # Cartesian product with budget, build specs
    de_configs = Vector{Int}[]
    _dead_end_cartesian!(de_configs, per_form_options, 1, Int[], budget)

    [
        _build_spec_from_edges(
            _topo_plus_binding_edges(topo, de, forms),
            forms, spec.reaction;
            shared_names, shared_atoms,
        )
        for de in de_configs
    ]
end

"""Get edges for a topology + dead-end config.

Includes all topology edges plus binding-only edges between
all forms in the mechanism. Isomerization/release edges between
dead-end forms are not added (those belong to activator configs).
"""
function _topo_plus_binding_edges(
    topo_edges::Vector{Tuple{Int,Int}},
    de_forms::Vector{Int},
    forms::Vector{EnzymeFormSpec},
)
    topo_forms = Set(Iterators.flatten(topo_edges))
    all_forms = sort(collect(union(topo_forms, Set(de_forms))))

    seen = Set{Tuple{Int,Int}}()
    edges = Tuple{Int,Int}[]

    # Topology edges (directed as given)
    for (a, b) in topo_edges
        key = minmax(a, b)
        key in seen && continue
        push!(seen, key)
        push!(edges, (a, b))
    end

    # Additional binding-only edges between all forms
    for i in 1:length(all_forms), j in (i + 1):length(all_forms)
        fi, fj = all_forms[i], all_forms[j]
        (fi, fj) in seen && continue
        is_valid, _, etype = edge_class(forms[fi], forms[fj])
        is_valid && etype == :binding || continue
        push!(seen, (fi, fj))
        push!(edges, (fi, fj))
    end

    edges
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

# ─── SS/RE + Constraint Enumeration ───────────────────────────

"""
    _enumerate_ress_and_constraints(spec, forms) → Vector{MechanismSpec}

Enumerate all valid RE/SS assignments and equivalent step constraints.
"""
function _enumerate_ress_and_constraints(
    spec::MechanismSpec,
    forms::Vector{EnzymeFormSpec},
)
    edges = _spec_to_edges(spec, forms)
    n = length(edges)
    equiv_groups = _find_equivalent_groups(edges, forms)
    results = MechanismSpec[]

    all_form_names = [f.name for f in forms]
    all_form_atoms = [_form_atoms(f.sites) for f in forms]
    rxn_tuples = _edges_to_reactions(edges, forms)

    # Iterate all non-zero bitmasks (at least one step must be SS)
    for re_mask in 0:((1 << n) - 2)
        eq_steps = Bool[(re_mask >> (i-1)) & 1 == 1 for i in 1:n]

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
            push!(results, MechanismSpec(
                spec.reaction, all_form_names, all_form_atoms,
                rxn_tuples, eq_steps, constraints,
            ))
        end
    end

    results
end

"""
Find groups of equivalent steps. Two steps are equivalent if both involve
the same substrate or regulator metabolite at the same site index.

Direction-independent: works regardless of edge ordering in `step_edges`,
since it checks which form has the occupied site rather than relying on
`edge_class` direction. Product sites are excluded (those edges are
catalytic release steps, not binding events).
"""
function _find_equivalent_groups(step_edges::Vector{Tuple{Int,Int}},
                                  forms::Vector{EnzymeFormSpec})
    binding_key = Dict{Tuple{Symbol,Int}, Vector{Int}}()
    for (i, (a, b)) in enumerate(step_edges)
        _, _, etype = edge_class(forms[a], forms[b])
        etype in (:binding, :release) || continue
        for k in 1:length(forms[a].sites)
            a_occ = forms[a].sites[k].atoms !== nothing
            b_occ = forms[b].sites[k].atoms !== nothing
            a_occ == b_occ && continue
            # Identify the bound form's site info
            site = a_occ ? forms[a].sites[k] : forms[b].sites[k]
            site.role == :prod && break
            key = (site.metabolite, site.index)
            push!(get!(binding_key, key, Int[]), i)
            break
        end
    end
    groups = Vector{Int}[]
    for (_, indices) in binding_key
        length(indices) >= 2 && push!(groups, sort(indices))
    end
    sort!(groups; by=first)
    groups
end

# ─── Main Entry Point ─────────────────────────────────────────

"""
    enumerate_mechanisms(reaction::EnzymeReaction; max_forms = 3 * n_sites(reaction))

Enumerate all valid mechanism topologies for the given reaction.

Returns a `MechanismIterator` that lazily yields `MechanismSpec` structs.
Each spec can be converted to an `EnzymeMechanism` via `EnzymeMechanism(spec)`.

The enumeration covers:
1. Catalytic form sets (standard and ping-pong cycles, multi-cycle unions)
2. Activator shadow variants (regulators participating in catalysis)
3. Dead-end complexes (inhibitors, abortive complexes, product inhibition)
4. RE/SS assignments (rapid-equilibrium vs steady-state per step)
5. Equivalent step constraints (shared parameters for equivalent binding steps)
"""
function enumerate_mechanisms(@nospecialize(reaction::EnzymeReaction);
                              max_forms::Int = 3 * n_sites(reaction))
    stages = enumerate_mechanism_stages(reaction; max_forms)
    stages.final
end

"""
    enumerate_mechanism_stages(reaction; max_forms) → NamedTuple

Run the enumeration pipeline and return intermediate results at each stage:

- `forms`: all enzyme forms
- `catalytic`: catalytic topologies (stage 1)
- `with_activator`: after activator configs (stage 2)
- `with_dead_end`: after dead-end configs (stage 3)
- `final`: lazy `MechanismIterator` for RE/SS + constraints (stage 4)
"""
function enumerate_mechanism_stages(@nospecialize(reaction::EnzymeReaction);
                                     max_forms::Int = 3 * n_sites(reaction))
    forms = enumerate_enzyme_forms(reaction)
    form_names = [f.name for f in forms]
    form_atoms = [_form_atoms(f.sites) for f in forms]

    catalytic = _enumerate_only_catalytic_mechanisms(
        forms, reaction; max_forms,
        shared_names=form_names, shared_atoms=form_atoms)

    with_activator = MechanismSpec[
        s for spec in catalytic
        for s in _generate_activator_configs(
            spec, forms, reaction;
            shared_names=form_names, shared_atoms=form_atoms)]
    filter!(s -> _used_form_count(s) <= max_forms, with_activator)

    with_dead_end = MechanismSpec[
        s for spec in with_activator
        for s in _enumerate_dead_end_configs(
            spec, forms; max_forms,
            shared_names=form_names, shared_atoms=form_atoms)]

    pre_ress = PreRessEntry[]
    for spec in with_dead_end
        edges = _spec_to_edges(spec, forms)
        equiv_groups = _find_equivalent_groups(edges, forms)
        rxn_tuples = _edges_to_reactions(edges, forms)
        rc = _count_ress_variants(length(edges), equiv_groups)
        push!(pre_ress, PreRessEntry(rxn_tuples, edges, equiv_groups, rc))
    end
    total = sum(e.ress_count for e in pre_ress; init=0)
    final = MechanismIterator(
        reaction, form_names, form_atoms, pre_ress, total)

    (; forms, catalytic, with_activator, with_dead_end, final)
end

# ─── MechanismSpec → EnzymeMechanism Conversion ──────────────

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
