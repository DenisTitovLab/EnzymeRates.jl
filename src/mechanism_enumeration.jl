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

"""Compute residual atoms after subtracting `to_remove` from `sub_atoms`.
Returns sorted Vector{Pair} or nothing if subtraction is invalid/trivial."""
function _atom_residual(sub_atoms, to_remove)
    sa = Dict{Symbol,Int}(sub_atoms)
    pa = Dict{Symbol,Int}(to_remove)
    all(get(sa, a, 0) >= c for (a, c) in pa) || return nothing
    res = Dict{Symbol,Int}(
        a => c - get(pa, a, 0) for (a, c) in sa if c > get(pa, a, 0))
    isempty(res) && return nothing
    r = sort([a => c for (a, c) in res]; by=first)
    r == sub_atoms ? nothing : r
end

"""Check if two forms represent a valid isomerization."""
function _is_valid_isomerization(fa::EnzymeFormSpec, fb::EnzymeFormSpec)
    diffs = [(k, fa.sites[k].role) for k in eachindex(fa.sites)
             if fa.sites[k].role in (:sub, :prod) &&
                fa.sites[k].atoms != fb.sites[k].atoms]
    any(d -> d[2] == :sub, diffs) &&
        any(d -> d[2] == :prod, diffs) || return false
    has_residual = any(diffs) do (k, role)
        role == :sub && any(x -> x !== nothing && x != fa.sites[k].full_atoms,
                           (fa.sites[k].atoms, fb.sites[k].atoms))
    end
    n_core = count(s -> s.role in (:sub, :prod), fa.sites)
    (has_residual || length(diffs) == n_core) &&
        _core_atoms(fa) == _core_atoms(fb)
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
                idx = findfirst(s -> s.role == :prod &&
                    s.full_atoms == sa.atoms, fa.sites)
                idx !== nothing ? fa.sites[idx].metabolite : nothing
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
        sub_atoms_vec = fatoms(sub_spec)
        sub_residuals = Vector{Pair{Symbol,Int}}[]
        for mask in 1:(2^length(prod_atoms_list) - 1)
            combined = reduce(mergewith(+),
                (prod_atoms_list[i]
                 for i in 1:length(prod_atoms_list)
                 if (mask >> (i-1)) & 1 == 1))
            r = _atom_residual(sub_atoms_vec, combined)
            r !== nothing && r ∉ sub_residuals && push!(sub_residuals, r)
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

function _pingpong_dfs!(cycles, forms, fidx, template, adj, free_idx,
                         all_subs, all_prods, prod_full, sequence,
                         bound_subs, released_prods, residual_state,
                         n_subs_bound, n_prods_released)
    if n_subs_bound == length(all_subs) &&
            n_prods_released == length(all_prods)
        sequence[end] == free_idx &&
            (s = Set{Int}(sequence); s ∉ cycles && push!(cycles, s))
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
                           all_subs, all_prods, prod_full, sequence,
                           bound_subs, released_prods, residual_state,
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
            pa = get(prod_full, prod, nothing)
            pa === nothing && continue
            for s in current_form.sites
                s.role == :sub && s.atoms !== nothing || continue
                res_vec = _atom_residual(s.atoms, pa)
                res_vec === nothing && continue

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
                _pingpong_dfs!(cycles, forms, fidx, template, adj,
                               free_idx, all_subs, all_prods, prod_full,
                               sequence, bound_subs, released_prods,
                               new_residual,
                               n_subs_bound, n_prods_released + 1)
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
            (s = Set{Int}(sequence); s ∉ cycles && push!(cycles, s))
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

"""Check whether a form set is pure sequential or pure ping-pong."""
function _is_pure_topology(form_set::Set{Int}, forms::Vector{EnzymeFormSpec})
    _has_res(f) = any(s -> s.role == :sub && s.atoms !== nothing &&
        s.atoms != s.full_atoms, f.sites)
    res_forms = [fi for fi in form_set if _has_res(forms[fi])]
    isempty(res_forms) && return true
    has_free_res = any(res_forms) do fi
        !any(s -> (s.role == :sub && s.atoms == s.full_atoms) ||
            (s.role == :prod && s.atoms !== nothing), forms[fi].sites)
    end
    has_free_res && !any(fi -> all(
        s -> s.role != :sub || s.atoms == s.full_atoms,
        forms[fi].sites), form_set)
end

"""Combine individual cycles into all unique multi-cycle unions."""
function _combine_form_sets(cycles::Vector{Set{Int}})
    isempty(cycles) && return Set{Int}[]
    n = length(cycles)
    result = Set{Int}[]
    seen = Set{Set{Int}}()
    for mask in 1:(1 << n) - 1
        merged = union((cycles[i] for i in 1:n
                        if (mask >> (i-1)) & 1 == 1)...)
        merged in seen && continue
        push!(seen, merged); push!(result, merged)
    end
    result
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

    prod_full = Dict(s.metabolite => s.full_atoms
                     for s in template if s.role == :prod)
    cycles = Set{Int}[]
    _pingpong_dfs!(cycles, forms, fidx, template, adj, free_idx,
                   sub_names, prod_names, prod_full, [free_idx],
                   Set{Symbol}(), Set{Symbol}(),
                   Dict{Symbol,Vector{Pair{Symbol,Int}}}(), 0, 0)

    combined = _combine_form_sets(cycles)
    filter!(fs -> _is_pure_topology(fs, forms) && length(fs) <= max_forms,
            combined)

    map(combined) do fs
        sorted = sort(collect(fs))
        edges = [(sorted[i], sorted[j])
                 for i in 1:length(sorted) for j in (i+1):length(sorted)
                 if haskey(adj, (sorted[i], sorted[j]))]
        MechanismSpec(reaction, edges,
                      fill(false, length(edges)), ParamConstraint[])
    end
end

# ─── Activator Configs ────────────────────────────────────────

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

    OptionT = Tuple{Vector{Tuple{Int,Int}}, Set{Tuple{Int,Int}}}
    result = MechanismSpec[]
    for spec in specs
        topo_forms = _used_form_set(spec)
        per_reg = Vector{Vector{OptionT}}()

        for reg in reg_names
            options = OptionT[(Tuple{Int,Int}[], Set{Tuple{Int,Int}}())]
            shadow_pairs = Tuple{Int,Int}[]
            for bi in sort(collect(topo_forms))
                pos = findfirst(s -> s.metabolite == reg &&
                    s.role == :reg && s.atoms === nothing,
                    forms[bi].sites)
                pos === nothing && continue
                si = _find_dead_end(fidx, forms[bi], (pos,))
                si !== nothing && push!(shadow_pairs, (bi, si))
            end

            if length(shadow_pairs) >= 2
                smap = Dict(shadow_pairs)
                mirror = [minmax(smap[a], smap[b])
                          for (a, b) in spec.edges
                          if haskey(smap, a) && haskey(smap, b)]
                mirrored = Set((a, b) for (a, b) in spec.edges
                               if haskey(smap, a) && haskey(smap, b))
                binding = [minmax(b, s) for (b, s) in shadow_pairs]

                # Non-essential: all binding + mirror cycle
                push!(options, ([binding; mirror], Set{Tuple{Int,Int}}()))

                # Essential: entry binding + mirror, remove mirrored base
                entry = findfirst(((bi, _),) -> all(
                    s -> s.role == :reg || s.atoms === nothing,
                    forms[bi].sites), shadow_pairs)
                if entry !== nothing
                    bp = shadow_pairs[entry]
                    push!(options,
                          ([minmax(bp...); mirror], mirrored))
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

    # Per topo form: compute regulator-subset options.
    # Enzyme forms are a Cartesian product of per-site options, so if
    # individual reg bindings exist, all subset combinations also exist.
    per_form_opts = Vector{Vector{Int}}[]
    for fi in topo_forms
        reg_sites = [k for k in eachindex(forms[fi].sites)
                     if forms[fi].sites[k].role == :reg &&
                        forms[fi].sites[k].atoms === nothing &&
                        k ∉ act_pos]
        n_reg = length(reg_sites)
        options = [Int[]]
        for mask in 1:((1 << n_reg) - 1)
            chosen = [reg_sites[k] for k in 1:n_reg
                      if (mask >> (k-1)) & 1 == 1]
            de_forms = Int[]
            for sub_mask in 1:((1 << length(chosen)) - 1)
                positions = [chosen[k] for k in 1:length(chosen)
                             if (sub_mask >> (k-1)) & 1 == 1]
                fi2 = _find_dead_end(fidx, forms[fi], positions)
                fi2 !== nothing && fi2 ∉ topo_set && push!(de_forms, fi2)
            end
            push!(options, sort!(de_forms))
        end
        push!(per_form_opts, options)
    end

    results = MechanismSpec[]
    for combo in Iterators.product(per_form_opts...)
        de = vcat(combo...)
        length(de) > budget && continue
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
        push!(results, MechanismSpec(spec.reaction, edges,
                      fill(false, length(edges)), ParamConstraint[]))
    end
    results
end

# ─── RE/SS + Constraint Lazy Generator ────────────────────────

"""Find equivalent binding step groups from edges."""
function _find_equivalent_groups(
    edges::Vector{Tuple{Int,Int}},
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    forms::Vector{EnzymeFormSpec},
)
    template = forms[1].sites
    binding_key = Dict{Symbol, Vector{Int}}()
    for (i, (a, b)) in enumerate(edges)
        info = get(adj, minmax(a, b), nothing)
        info === nothing && continue
        info.type in (:binding, :release) || continue
        any(s -> s.metabolite == info.metabolite && s.role == :prod,
            template) && continue
        push!(get!(binding_key, info.metabolite, Int[]), i)
    end
    groups = [sort(v) for v in values(binding_key) if length(v) >= 2]
    sort!(groups; by=first)
end

"""Count RE/SS + constraint variants using closed-form formula.

For `n` edges and `k` non-overlapping equiv groups of sizes g₁,...,gₖ:
  f(n; g₁,...,gₖ) = 2^(n - Σgᵢ) × ∏(2^gᵢ + 2) - 2^k
"""
function _count_ress_variants(n_edges::Int, equiv_groups)
    n_edges == 0 && return 0
    k = length(equiv_groups)
    sum_g = sum(length, equiv_groups; init=0)
    (1 << (n_edges - sum_g)) *
        prod(1 << length(g) + 2 for g in equiv_groups; init=1) - (1 << k)
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
        pfx, sfxs = eq_steps[group[1]] ? ("K", ("",)) : ("k", ("f", "r"))
        for j in 2:length(group), sfx in sfxs
            push!(constraints, (Symbol("$pfx$(group[j])$sfx"), 1,
                               [(Symbol("$pfx$(group[1])$sfx"), 1)]))
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
        (forms[i].name, Tuple(Tuple.(_form_atoms(forms[i].sites))))
        for i in 1:length(forms) if i in used)
    species = (subs_t, prods_t, regs_t, enzs_t)

    reactions = Tuple(
        let info = adj[minmax(a, b)]
            lhs = info.type == :binding ?
                (forms[a].name, info.metabolite) : (forms[a].name,)
            rhs = info.type == :release ?
                (forms[b].name, info.metabolite) : (forms[b].name,)
            (lhs, rhs)
        end
        for (a, b) in spec.edges)
    eq_steps = Tuple(spec.equilibrium_steps)
    constraints = Tuple(
        (target, coeff, Tuple(Tuple.(factors)))
        for (target, coeff, factors) in spec.param_constraints)
    EnzymeMechanism(species, reactions, eq_steps, constraints)
end
