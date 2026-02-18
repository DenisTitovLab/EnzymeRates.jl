# ─── Data Types ─────────────────────────────────────────────

const SiteState = @NamedTuple{
    metabolite::Symbol,
    atoms::Union{Nothing, Vector{Pair{Symbol,Int}}},
    role::Symbol,
    full_atoms::Vector{Pair{Symbol,Int}},
}

struct EnzymeFormSpec
    name::Symbol
    sites::Vector{SiteState}
end

const ParamConstraint = Tuple{Symbol, Int, Vector{Tuple{Symbol, Int}}}

struct MechanismSpec
    reaction::Any
    edges::Vector{Tuple{Int,Int}}
    equilibrium_steps::Vector{Bool}
    param_constraints::Vector{ParamConstraint}
end

abstract type EnumerationStage end
struct Catalytic    <: EnumerationStage end
struct WithActivator <: EnumerationStage end
struct WithDeadEnd   <: EnumerationStage end
struct Full          <: EnumerationStage end

struct MechanismIterator
    inner::Any
    total::Int
end

Base.eltype(::Type{MechanismIterator}) = MechanismSpec
Base.IteratorSize(::Type{MechanismIterator}) = Base.HasLength()
Base.length(iter::MechanismIterator) = iter.total
Base.iterate(iter::MechanismIterator) = iterate(iter.inner)
Base.iterate(iter::MechanismIterator, state) = iterate(iter.inner, state)

# ─── Form Index ─────────────────────────────────────────────

"""Build Dict mapping site-atoms tuple → form index for O(1) lookup."""
_build_form_index(forms) =
    Dict(Tuple(s.atoms for s in f.sites) => i for (i, f) in enumerate(forms))

"""Look up form index by occupancy pattern."""
function _lookup_form(fidx, template, occupied_subs, occupied_prods,
                      residual=nothing)
    key = ntuple(length(template)) do i
        s = template[i]
        if s.role == :reg
            nothing
        elseif s.role == :sub
            if residual !== nothing && haskey(residual, s.metabolite)
                residual[s.metabolite]
            elseif s.metabolite in occupied_subs
                s.full_atoms
            else
                nothing
            end
        else # :prod
            s.metabolite in occupied_prods ? s.full_atoms : nothing
        end
    end
    get(fidx, key, nothing)
end

# ─── Edge Classification + Adjacency ─────────────────────────

const EdgeInfo = @NamedTuple{type::Symbol, metabolite::Union{Nothing,Symbol}}

"""Sum atoms across sites, optionally filtering by role."""
function _sum_atoms(sites, filter_roles=nothing)
    atoms = Dict{Symbol,Int}()
    for s in sites
        filter_roles !== nothing && !(s.role in filter_roles) && continue
        s.atoms === nothing && continue
        for (a, c) in s.atoms; atoms[a] = get(atoms, a, 0) + c; end
    end
    atoms
end
_core_atoms(form::EnzymeFormSpec) = _sum_atoms(form.sites, (:sub, :prod))
_form_atoms(sites) = sort([a => c for (a, c) in _sum_atoms(sites)]; by=first)

"""Check if two forms represent a valid isomerization."""
function _is_valid_isomerization(fa::EnzymeFormSpec, fb::EnzymeFormSpec)
    has_sub = false; has_prod = false; has_residual = false
    n_core = 0; n_diff = 0
    for k in eachindex(fa.sites)
        s = fa.sites[k]
        s.role in (:sub, :prod) || continue
        n_core += 1
        fa.sites[k].atoms == fb.sites[k].atoms && continue
        n_diff += 1
        s.role == :sub ? (has_sub = true) : (has_prod = true)
        if s.role == :sub
            a, b = fa.sites[k].atoms, fb.sites[k].atoms
            if (a !== nothing && a != s.full_atoms) ||
               (b !== nothing && b != s.full_atoms)
                has_residual = true
            end
        end
    end
    (!has_sub || !has_prod) && return false
    has_residual && return _core_atoms(fa) == _core_atoms(fb)
    # Standard case: ALL sub/prod sites must differ, with all-full↔all-empty
    # pattern. Since no residual and all sites differ, each is full or empty.
    n_diff != n_core && return false
    a_sub_occ = all(fa.sites[k].atoms !== nothing
        for k in eachindex(fa.sites) if fa.sites[k].role == :sub)
    a_prod_occ = all(fa.sites[k].atoms !== nothing
        for k in eachindex(fa.sites) if fa.sites[k].role == :prod)
    a_sub_occ != a_prod_occ && _core_atoms(fa) == _core_atoms(fb)
end

"""Classify edge between two enzyme forms → EdgeInfo or nothing."""
function _classify_edge(fa::EnzymeFormSpec, fb::EnzymeFormSpec)
    diffs = Int[]
    for k in eachindex(fa.sites)
        fa.sites[k].atoms != fb.sites[k].atoms && push!(diffs, k)
    end
    isempty(diffs) && return nothing

    if length(diffs) == 1
        k = diffs[1]
        sa, sb = fa.sites[k], fb.sites[k]
        if sa.atoms === nothing && sb.atoms !== nothing
            (sa.role in (:sub, :prod) && sb.atoms != sa.full_atoms) &&
                return nothing
            return (type=:binding, metabolite=sa.metabolite)
        end
        if sb.atoms === nothing && sa.atoms !== nothing
            met = if sa.atoms == sa.full_atoms || sa.role == :reg
                sa.metabolite
            else
                m = nothing
                for s in fa.sites
                    s.role == :prod && s.full_atoms == sa.atoms &&
                        (m = s.metabolite; break)
                end
                m
            end
            met !== nothing && return (type=:release, metabolite=met)
        end
        return nothing
    end

    all(k -> fa.sites[k].role in (:sub, :prod), diffs) || return nothing
    _is_valid_isomerization(fa, fb) &&
        return (type=:isomerization, metabolite=nothing)
    nothing
end

"""Build adjacency dict from enzyme forms."""
function _build_adjacency(forms::Vector{EnzymeFormSpec})
    adj = Dict{Tuple{Int,Int}, EdgeInfo}()
    for i in 1:length(forms), j in (i+1):length(forms)
        info = _classify_edge(forms[i], forms[j])
        info !== nothing && (adj[(i, j)] = info)
    end
    adj
end

# ─── Enzyme Form Enumeration ─────────────────────────────────

"""
    enumerate_enzyme_forms(reaction::EnzymeReaction)

Enumerate all possible enzyme forms for the given reaction.
"""
function enumerate_enzyme_forms(reaction::EnzymeReaction{S,P,R}) where {S,P,R}
    fatoms(spec) = sort([a => c for (a, c) in spec[2]]; by=first)
    astr(v) = join(string(s) * (c > 1 ? string(c) : "") for (s, c) in v)

    # Ping-pong residuals
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
            all(get(sub_atoms, atom, 0) >= count
                for (atom, count) in combined) || continue
            residual = Dict{Symbol,Int}(
                atom => count - get(combined, atom, 0)
                for (atom, count) in sub_atoms
                if count > get(combined, atom, 0))
            if !isempty(residual) && residual != sub_atoms
                r = sort([k => v for (k, v) in residual]; by=first)
                r ∉ sub_residuals && push!(sub_residuals, r)
            end
        end
        !isempty(sub_residuals) && (residuals[sub_spec[1]] = sub_residuals)
    end

    # Build per-site data
    mets = Symbol[]; fulls = Vector{Pair{Symbol,Int}}[]
    roles = Symbol[]
    OptT = Tuple{Union{Nothing, Vector{Pair{Symbol,Int}}}, String}
    opts = Vector{OptT}[]
    for (group, role) in ((S, :sub), (P, :prod), (R, :reg))
        for spec in group
            push!(mets, spec[1]); push!(fulls, fatoms(spec))
            push!(roles, role)
            site_opts = OptT[(nothing, "0"),
                             (fatoms(spec), string(spec[1]))]
            for r in get(residuals, spec[1], Vector{Pair{Symbol,Int}}[])
                push!(site_opts, (r, astr(r)))
            end
            push!(opts, site_opts)
        end
    end

    # Cartesian product with exclusion filter
    n = length(mets)
    forms = EnzymeFormSpec[]
    total = prod(length(o) for o in opts)
    name_parts = Vector{String}(undef, n)
    sites_buf = Vector{SiteState}(undef, n)
    for idx in 0:total-1
        rem = idx
        all_sub_full = true; all_prod_full = true
        any_sub_occ = false; any_prod_occ = false
        for i in n:-1:1
            content, label = opts[i][rem % length(opts[i]) + 1]
            rem ÷= length(opts[i])
            name_parts[i] = label
            sites_buf[i] = (
                metabolite=mets[i],
                atoms=content === nothing ? nothing : copy(content),
                role=roles[i], full_atoms=fulls[i],
            )
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
        push!(forms, EnzymeFormSpec(
            Symbol("E_" * join(name_parts, "_")), copy(sites_buf)))
    end
    forms
end

# ─── Catalytic Cycle Enumeration ──────────────────────────────

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

"""Build a standard catalytic cycle as a form set."""
function _build_standard_form_set(fidx, template, adj, free_idx,
                                   sub_perm, prod_perm)
    form_set = Set{Int}(free_idx)
    occupied_subs = Set{Symbol}()
    prev = free_idx
    for sub in sub_perm
        push!(occupied_subs, sub)
        fi = _lookup_form(fidx, template, occupied_subs, Set{Symbol}())
        fi === nothing && return nothing
        push!(form_set, fi)
        prev = fi
    end
    all_prods = Set{Symbol}(prod_perm)
    fi = _lookup_form(fidx, template, Set{Symbol}(), all_prods)
    fi === nothing && return nothing
    haskey(adj, minmax(prev, fi)) || return nothing
    push!(form_set, fi)
    for i in eachindex(prod_perm)
        remaining = Set{Symbol}(@view prod_perm[i+1:end])
        fi = isempty(remaining) ? free_idx :
            _lookup_form(fidx, template, Set{Symbol}(), remaining)
        fi === nothing && return nothing
        push!(form_set, fi)
    end
    form_set
end

function _pingpong_dfs!(cycles, forms, fidx, template, adj, free_idx,
                         all_subs, all_prods, sequence, bound_subs,
                         released_prods, residual_state,
                         n_subs_bound, n_prods_released)
    if n_subs_bound == length(all_subs) &&
            n_prods_released == length(all_prods)
        sequence[end] == free_idx &&
            _add_unique!(cycles, Set{Int}(sequence))
        return
    end

    # Option 1: Bind next substrate
    for sub in all_subs
        sub in bound_subs && continue
        push!(bound_subs, sub)
        fi = _lookup_form(fidx, template, bound_subs, Set{Symbol}(),
                          residual_state)
        if fi !== nothing && haskey(adj, minmax(sequence[end], fi))
            push!(sequence, fi)
            _pingpong_dfs!(cycles, forms, fidx, template, adj, free_idx,
                           all_subs, all_prods, sequence, bound_subs,
                           released_prods, residual_state,
                           n_subs_bound + 1, n_prods_released)
            pop!(sequence)
        end
        delete!(bound_subs, sub)
    end

    # Option 2: Standard isomerization → release remaining products
    if length(bound_subs) == length(all_subs) &&
            n_subs_bound > n_prods_released
        remaining = setdiff(Set{Symbol}(all_prods), released_prods)
        if !isempty(remaining)
            fi = _lookup_form(fidx, template, Set{Symbol}(), remaining)
            if fi !== nothing && haskey(adj, minmax(sequence[end], fi))
                push!(sequence, fi)
                _release_prods_dfs!(cycles, fidx, template, adj, free_idx,
                                    sequence, remaining)
                pop!(sequence)
            end
        end
    end

    # Option 3: Ping-pong isomerization + product release (two steps)
    if !isempty(bound_subs)
        current_form = forms[sequence[end]]
        for prod in all_prods
            prod in released_prods && continue
            for s in current_form.sites
                s.role == :sub && s.atoms !== nothing || continue
                prod_atoms = nothing
                for ps in current_form.sites
                    if ps.role == :prod && ps.metabolite == prod
                        prod_atoms = ps.full_atoms; break
                    end
                end
                prod_atoms === nothing && continue
                sa = Dict{Symbol,Int}(a => c for (a, c) in s.atoms)
                pa = Dict{Symbol,Int}(a => c for (a, c) in prod_atoms)
                all(get(sa, a, 0) >= c for (a, c) in pa) || continue
                res = Dict{Symbol,Int}(
                    a => c - get(pa, a, 0)
                    for (a, c) in sa if c > get(pa, a, 0))
                isempty(res) && continue
                res_vec = sort([a => c for (a, c) in res]; by=first)
                res_vec == s.atoms && continue

                new_residual = copy(residual_state)
                new_residual[s.metabolite] = res_vec

                inter = _lookup_form(fidx, template, bound_subs,
                                     Set{Symbol}([prod]), new_residual)
                inter === nothing && continue
                haskey(adj, minmax(sequence[end], inter)) || continue

                post = _lookup_form(fidx, template, bound_subs,
                                    Set{Symbol}(), new_residual)
                post === nothing && continue
                haskey(adj, minmax(inter, post)) || continue

                push!(sequence, inter, post)
                push!(released_prods, prod)
                old_residual = copy(residual_state)
                merge!(residual_state, new_residual)
                _pingpong_dfs!(cycles, forms, fidx, template, adj,
                               free_idx, all_subs, all_prods, sequence,
                               bound_subs, released_prods, residual_state,
                               n_subs_bound, n_prods_released + 1)
                merge!(residual_state, old_residual)
                for k in keys(new_residual)
                    haskey(old_residual, k) || delete!(residual_state, k)
                end
                delete!(released_prods, prod)
                pop!(sequence); pop!(sequence)
            end
        end
    end
end

"""Release remaining products after isomerization."""
function _release_prods_dfs!(cycles, fidx, template, adj, free_idx,
                              sequence, remaining)
    if isempty(remaining)
        sequence[end] == free_idx &&
            _add_unique!(cycles, Set{Int}(sequence))
        return
    end
    for prod in collect(remaining)
        delete!(remaining, prod)
        fi = isempty(remaining) ? free_idx :
            _lookup_form(fidx, template, Set{Symbol}(), remaining)
        if fi !== nothing && haskey(adj, minmax(sequence[end], fi))
            push!(sequence, fi)
            _release_prods_dfs!(cycles, fidx, template, adj, free_idx,
                                sequence, remaining)
            pop!(sequence)
        end
        push!(remaining, prod)
    end
end

function _add_unique!(sets::Vector{Set{Int}}, new_set::Set{Int})
    for existing in sets
        existing == new_set && return
    end
    push!(sets, new_set)
end

"""Check whether a form set is pure sequential or pure ping-pong."""
function _is_pure_topology(form_set::Set{Int}, forms::Vector{EnzymeFormSpec})
    _is_residual(s) = s.role == :sub && s.atoms !== nothing &&
        s.atoms != s.full_atoms
    _is_free(f) = !any(s -> (s.role == :sub && s.atoms == s.full_atoms) ||
        (s.role == :prod && s.atoms !== nothing), f.sites)
    _all_subs_full(f) = all(s -> s.role != :sub || s.atoms == s.full_atoms,
        f.sites)
    has_residual = any(fi -> any(_is_residual, forms[fi].sites), form_set)
    !has_residual && return true
    has_free_inter = any(form_set) do fi
        any(_is_residual, forms[fi].sites) && _is_free(forms[fi])
    end
    has_free_inter && !any(fi -> _all_subs_full(forms[fi]), form_set)
end

"""Combine individual cycles into multi-cycle unions via BFS."""
function _combine_form_sets(cycles::Vector{Set{Int}})
    isempty(cycles) && return Set{Int}[]
    result = copy(cycles)
    seen = Set{Set{Int}}(cycles)
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

"""Derive edges from adjacency for a set of form indices."""
function _derive_edges(adj::Dict{Tuple{Int,Int}, EdgeInfo}, form_set)
    sorted = sort(collect(form_set))
    edges = Tuple{Int,Int}[]
    for i in 1:length(sorted), j in (i+1):length(sorted)
        haskey(adj, (sorted[i], sorted[j])) &&
            push!(edges, (sorted[i], sorted[j]))
    end
    edges
end

"""Enumerate catalytic topologies → Vector{MechanismSpec}."""
function _catalytic_topologies(
    forms::Vector{EnzymeFormSpec},
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    fidx, template,
    @nospecialize(reaction::EnzymeReaction);
    max_forms::Int,
)
    sub_names = [s[1] for s in substrates(reaction)]
    prod_names = [p[1] for p in products(reaction)]
    free_idx = findfirst(
        f -> all(s -> s.atoms === nothing, f.sites), forms)
    free_idx === nothing && return MechanismSpec[]

    cycles = Set{Int}[]

    for sub_perm in _permutations(sub_names)
        for prod_perm in _permutations(prod_names)
            fs = _build_standard_form_set(
                fidx, template, adj, free_idx, sub_perm, prod_perm)
            fs !== nothing && _add_unique!(cycles, fs)
        end
    end

    # Ping-pong cycles
    has_residuals = any(forms) do f
        any(f.sites) do s
            s.role == :sub && s.atoms !== nothing && s.atoms != s.full_atoms
        end
    end
    if has_residuals
        _pingpong_dfs!(cycles, forms, fidx, template, adj, free_idx,
                       sub_names, prod_names, [free_idx], Set{Symbol}(),
                       Set{Symbol}(),
                       Dict{Symbol,Vector{Pair{Symbol,Int}}}(), 0, 0)
    end

    combined = _combine_form_sets(cycles)
    filter!(fs -> _is_pure_topology(fs, forms) && length(fs) <= max_forms,
            combined)

    map(combined) do fs
        edges = _derive_edges(adj, fs)
        MechanismSpec(reaction, edges,
                      fill(false, length(edges)), ParamConstraint[])
    end
end

# ─── Activator Configs ────────────────────────────────────────

"""Find shadow form: base with regulator also bound."""
function _find_shadow(fidx, base::EnzymeFormSpec, reg::Symbol)
    any(s -> s.metabolite == reg && s.role == :reg && s.atoms === nothing,
        base.sites) || return nothing
    key = Tuple(
        s.metabolite == reg && s.role == :reg ? s.full_atoms : s.atoms
        for s in base.sites
    )
    get(fidx, key, nothing)
end

"""Set of form indices in a spec's edges."""
_used_form_set(spec::MechanismSpec) = Set(Iterators.flatten(spec.edges))

"""Generate activator variants for specs."""
function _expand_activators(
    specs::Vector{MechanismSpec},
    forms::Vector{EnzymeFormSpec},
    fidx,
    @nospecialize(reaction::EnzymeReaction);
    max_forms::Int,
)
    reg_names = [r[1] for r in regulators(reaction)]
    isempty(reg_names) && return specs

    result = MechanismSpec[]
    for spec in specs
        topo_forms = _used_form_set(spec)
        OptionT = Tuple{Vector{Tuple{Int,Int}}, Set{Tuple{Int,Int}}}
        per_reg = Vector{Vector{OptionT}}()

        for reg in reg_names
            options = OptionT[]
            push!(options, (Tuple{Int,Int}[], Set{Tuple{Int,Int}}()))

            shadow_pairs = Tuple{Int,Int}[]
            for bi in sort(collect(topo_forms))
                si = _find_shadow(fidx, forms[bi], reg)
                si !== nothing && push!(shadow_pairs, (bi, si))
            end

            if length(shadow_pairs) >= 2
                shadow_map = Dict(bi => si for (bi, si) in shadow_pairs)
                shadow_cycle = Tuple{Int,Int}[]
                shadowed_base = Set{Tuple{Int,Int}}()
                for (a, b) in spec.edges
                    sa = get(shadow_map, a, nothing)
                    sb = get(shadow_map, b, nothing)
                    if sa !== nothing && sb !== nothing
                        push!(shadow_cycle,
                              (min(sa, sb), max(sa, sb)))
                        push!(shadowed_base, (a, b))
                    end
                end
                binding = [(min(b, s), max(b, s))
                           for (b, s) in shadow_pairs]
                full_shadow = [binding; shadow_cycle]

                # Non-essential
                push!(options, (full_shadow, Set{Tuple{Int,Int}}()))

                # Essential
                entry_idx = findfirst(shadow_pairs) do (bi, _)
                    all(s -> s.role == :reg || s.atoms === nothing,
                        forms[bi].sites)
                end
                if entry_idx !== nothing
                    bp = shadow_pairs[entry_idx]
                    essential = [(min(bp...), max(bp...));
                                 shadow_cycle]
                    push!(options, (essential, shadowed_base))
                end
            end
            push!(per_reg, options)
        end

        for combo in Iterators.product(per_reg...)
            removals = union((r for (_, r) in combo)...)
            merged = [e for e in spec.edges if e ∉ removals]
            for (add, _) in combo; append!(merged, add); end
            new_spec = MechanismSpec(reaction, merged,
                fill(false, length(merged)), ParamConstraint[])
            length(_used_form_set(new_spec)) <= max_forms &&
                push!(result, new_spec)
        end
    end
    result
end

# ─── Dead-End Configs ─────────────────────────────────────────

"""Find dead-end form: base + specified reg positions occupied."""
function _find_dead_end(fidx, base::EnzymeFormSpec, occ_positions)
    key = Tuple(
        k in occ_positions ? s.full_atoms : s.atoms
        for (k, s) in enumerate(base.sites)
    )
    get(fidx, key, nothing)
end

"""Enumerate dead-end binding configurations."""
function _dead_end_configs(
    spec::MechanismSpec,
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    forms::Vector{EnzymeFormSpec},
    fidx;
    max_forms::Int,
)
    topo_set = _used_form_set(spec)
    topo_forms = sort(collect(topo_set))
    budget = max_forms - length(topo_forms)
    budget < 0 && return MechanismSpec[]

    # Activator positions: reg sites occupied in any topo form
    act_pos = Set{Int}()
    for fi in topo_forms, k in eachindex(forms[fi].sites)
        forms[fi].sites[k].role == :reg &&
            forms[fi].sites[k].atoms !== nothing &&
            push!(act_pos, k)
    end

    # Per topo form: compute regulator-subset options
    per_form_opts = Vector{Vector{Int}}[]
    for fi in topo_forms
        reg_sites = [k for k in eachindex(forms[fi].sites)
                     if forms[fi].sites[k].role == :reg &&
                        forms[fi].sites[k].atoms === nothing &&
                        k ∉ act_pos]
        n_reg = length(reg_sites)
        seen = Set{Vector{Int}}()
        options = [Int[]]
        push!(seen, Int[])

        for mask in 1:((1 << n_reg) - 1)
            chosen = [reg_sites[k] for k in 1:n_reg
                      if (mask >> (k-1)) & 1 == 1]
            de_forms = Int[]
            valid = true
            for sub_mask in 1:((1 << length(chosen)) - 1)
                positions = [chosen[k] for k in 1:length(chosen)
                             if (sub_mask >> (k-1)) & 1 == 1]
                fi2 = _find_dead_end(fidx, forms[fi], positions)
                if fi2 === nothing
                    valid = false; break
                end
                fi2 in topo_set && continue
                push!(de_forms, fi2)
            end
            !valid && continue
            sort!(de_forms)
            if de_forms ∉ seen
                push!(seen, de_forms)
                push!(options, de_forms)
            end
        end
        push!(per_form_opts, options)
    end

    de_configs = Vector{Int}[]
    _dead_end_cartesian!(de_configs, per_form_opts, 1, Int[], budget)

    map(de_configs) do de
        existing = Set(minmax(e...) for e in spec.edges)
        edges = copy(spec.edges)
        all_sorted = sort(collect(union(topo_set, Set(de))))
        for i in 1:length(all_sorted), j in (i+1):length(all_sorted)
            fi, fj = all_sorted[i], all_sorted[j]
            (fi, fj) in existing && continue
            info = get(adj, (fi, fj), nothing)
            info !== nothing && info.type == :binding || continue
            push!(existing, (fi, fj)); push!(edges, (fi, fj))
        end
        MechanismSpec(spec.reaction, edges,
                      fill(false, length(edges)), ParamConstraint[])
    end
end

function _dead_end_cartesian!(configs, options, idx, current, budget)
    if idx > length(options)
        push!(configs, copy(current)); return
    end
    for option in options[idx]
        length(option) > budget && continue
        append!(current, option)
        _dead_end_cartesian!(configs, options, idx + 1,
                             current, budget - length(option))
        resize!(current, length(current) - length(option))
    end
end

# ─── RE/SS + Constraint Lazy Generator ────────────────────────

"""Find equivalent binding step groups from edges."""
function _find_equivalent_groups(
    edges::Vector{Tuple{Int,Int}},
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    forms::Vector{EnzymeFormSpec},
)
    binding_key = Dict{Symbol, Vector{Int}}()
    for (i, (a, b)) in enumerate(edges)
        info = get(adj, minmax(a, b), nothing)
        info === nothing && continue
        info.type in (:binding, :release) || continue
        for k in eachindex(forms[a].sites)
            forms[a].sites[k].atoms == forms[b].sites[k].atoms && continue
            site = forms[a].sites[k].atoms !== nothing ?
                forms[a].sites[k] : forms[b].sites[k]
            site.role == :prod && break
            push!(get!(binding_key, site.metabolite, Int[]), i)
            break
        end
    end
    groups = [sort(v) for v in values(binding_key) if length(v) >= 2]
    sort!(groups; by=first)
end

"""Count RE/SS + constraint variants."""
function _count_ress_variants(n_edges::Int, equiv_groups)
    n_edges == 0 && return 0
    count = 0
    for re_mask in 0:((1 << n_edges) - 2)
        n_valid = 0
        for group in equiv_groups
            first_re = (re_mask >> (group[1] - 1)) & 1
            all(((re_mask >> (s - 1)) & 1) == first_re
                for s in group) && (n_valid += 1)
        end
        count += 1 << n_valid
    end
    count
end

"""Generate RE/SS + constraint variants lazily."""
function _ress_variants(spec, adj, forms)
    edges = spec.edges
    n = length(edges)
    equiv_groups = _find_equivalent_groups(edges, adj, forms)

    Iterators.flatten(
        Iterators.map(0:((1 << n) - 2)) do re_mask
            eq_steps = Bool[(re_mask >> (i-1)) & 1 == 1 for i in 1:n]
            valid_groups = [g for g in equiv_groups
                if all(eq_steps[s] == eq_steps[g[1]] for s in g)]
            n_vg = length(valid_groups)

            (MechanismSpec(spec.reaction, edges, eq_steps,
                _build_constraints(cmask, eq_steps, valid_groups))
             for cmask in 0:((1 << n_vg) - 1))
        end
    )
end

function _build_constraints(cmask, eq_steps, valid_groups)
    constraints = ParamConstraint[]
    for (gi, group) in enumerate(valid_groups)
        ((cmask >> (gi-1)) & 1) == 0 && continue
        if eq_steps[group[1]]
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
    constraints
end

# ─── Pipeline ─────────────────────────────────────────────────

"""
    enumerate_mechanisms(reaction; stage=Full(), max_forms=...)

Enumerate valid mechanism topologies for the given reaction.
"""
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    stage::EnumerationStage=Full(),
    max_forms::Int=3 * (length(substrates(reaction)) +
                        length(products(reaction)) +
                        length(regulators(reaction))),
)
    forms = enumerate_enzyme_forms(reaction)
    adj = _build_adjacency(forms)
    fidx = _build_form_index(forms)
    template = forms[1].sites

    catalytic = _catalytic_topologies(
        forms, adj, fidx, template, reaction; max_forms)
    stage isa Catalytic && return catalytic

    with_act = _expand_activators(
        catalytic, forms, fidx, reaction; max_forms)
    stage isa WithActivator && return with_act

    de_iter = Iterators.flatten(Iterators.map(
        s -> _dead_end_configs(s, adj, forms, fidx; max_forms), with_act))
    stage isa WithDeadEnd && return de_iter

    # Full: materialize dead-end, compute count, wrap lazy RE/SS
    de_specs = collect(de_iter)
    total = sum(de_specs; init=0) do s
        eg = _find_equivalent_groups(s.edges, adj, forms)
        _count_ress_variants(length(s.edges), eg)
    end
    inner = Iterators.flatten(Iterators.map(
        s -> _ress_variants(s, adj, forms), de_specs))
    MechanismIterator(inner, total)
end

# ─── MechanismSpec → EnzymeMechanism Conversion ──────────────

"""Convert edges to reaction tuples using adjacency."""
function _edges_to_reactions(
    edges::Vector{Tuple{Int,Int}},
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    forms::Vector{EnzymeFormSpec},
)
    map(edges) do (a, b)
        info = adj[minmax(a, b)]
        fa, fb = forms[a].name, forms[b].name
        if info.type == :binding
            ([fa, info.metabolite], [fb])
        elseif info.type == :release
            ([fa], [fb, info.metabolite])
        else
            ([fa], [fb])
        end
    end
end

"""
    EnzymeMechanism(spec::MechanismSpec)

Convert a `MechanismSpec` to a type-parameterized `EnzymeMechanism`.
"""
function EnzymeMechanism(spec::MechanismSpec)
    rxn = spec.reaction
    forms = enumerate_enzyme_forms(rxn)
    adj = _build_adjacency(forms)

    used = _used_form_set(spec)
    subs_t = Tuple((n, a) for (n, a, _) in substrates(rxn))
    prods_t = Tuple((n, a) for (n, a, _) in products(rxn))
    regs_t = Tuple((n, a) for (n, a, _) in regulators(rxn))
    enzs_t = Tuple(
        (forms[i].name, Tuple(Tuple.(fa)))
        for (i, fa) in (
            (i, _form_atoms(forms[i].sites))
            for i in 1:length(forms))
        if i in used
    )
    species = (subs_t, prods_t, regs_t, enzs_t)

    rxn_tuples = _edges_to_reactions(spec.edges, adj, forms)
    reactions = Tuple((Tuple(r[1]), Tuple(r[2])) for r in rxn_tuples)
    eq_steps = Tuple(spec.equilibrium_steps)
    constraints = Tuple(
        (target, coeff, Tuple(Tuple.(factors)))
        for (target, coeff, factors) in spec.param_constraints
    )
    EnzymeMechanism(species, reactions, eq_steps, constraints)
end
