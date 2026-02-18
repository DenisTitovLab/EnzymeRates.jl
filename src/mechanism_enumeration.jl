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
function _lookup_form(fidx, template, occ_subs, occ_prods, residual=nothing)
    key = ntuple(length(template)) do i
        s = template[i]
        s.role == :reg && return nothing
        if s.role == :sub
            r = residual !== nothing ? get(residual, s.metabolite, nothing) :
                nothing
            return r !== nothing ? r :
                s.metabolite in occ_subs ? s.full_atoms : nothing
        end
        s.metabolite in occ_prods ? s.full_atoms : nothing
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

"""Classify edge between two enzyme forms → EdgeInfo or nothing."""
function _classify_edge(fa::EnzymeFormSpec, fb::EnzymeFormSpec)
    diffs = [k for k in eachindex(fa.sites)
             if fa.sites[k].atoms != fb.sites[k].atoms]
    isempty(diffs) && return nothing
    if length(diffs) == 1
        k = diffs[1]; sa, sb = fa.sites[k], fb.sites[k]
        if sa.atoms === nothing && sb.atoms !== nothing
            (sa.role in (:sub, :prod) && sb.atoms != sa.full_atoms) &&
                return nothing
            return (type=:binding, metabolite=sa.metabolite)
        elseif sb.atoms === nothing && sa.atoms !== nothing
            met = sa.atoms == sa.full_atoms || sa.role == :reg ?
                sa.metabolite :
                let idx = findfirst(s -> s.role == :prod &&
                        s.full_atoms == sa.atoms, fa.sites)
                    idx !== nothing ? fa.sites[idx].metabolite : nothing
                end
            return met !== nothing ? (type=:release, metabolite=met) :
                nothing
        end
        return nothing
    end
    all(k -> fa.sites[k].role in (:sub, :prod), diffs) || return nothing
    any(k -> fa.sites[k].role == :sub, diffs) &&
        any(k -> fa.sites[k].role == :prod, diffs) || return nothing
    has_res = any(k -> fa.sites[k].role == :sub && any(
        x -> x !== nothing && x != fa.sites[k].full_atoms,
        (fa.sites[k].atoms, fb.sites[k].atoms)), diffs)
    (has_res || length(diffs) == count(
        s -> s.role in (:sub, :prod), fa.sites)) &&
        _core_atoms(fa) == _core_atoms(fb) &&
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
    pa_list = [Dict{Symbol,Int}(a => c for (a, c) in p[2])
               for p in P if !isempty(p[2])]
    residuals = Dict{Symbol,Vector{Vector{Pair{Symbol,Int}}}}()
    for ss in S
        (isempty(ss[2]) || isempty(pa_list)) && continue
        sav = fatoms(ss); sr = Vector{Pair{Symbol,Int}}[]
        for m in 1:(1 << length(pa_list)) - 1
            comb = reduce(mergewith(+), (pa_list[i]
                for i in 1:length(pa_list) if (m >> (i - 1)) & 1 == 1))
            r = _atom_residual(sav, comb)
            r !== nothing && r ∉ sr && push!(sr, r)
        end
        !isempty(sr) && (residuals[ss[1]] = sr)
    end
    # Build per-site options
    mets = Symbol[]; fulls = Vector{Pair{Symbol,Int}}[]; roles = Symbol[]
    OptT = Tuple{Union{Nothing,Vector{Pair{Symbol,Int}}},String}
    opts = Vector{OptT}[]
    for (group, role) in ((S, :sub), (P, :prod), (R, :reg)), spec in group
        push!(mets, spec[1]); push!(fulls, fatoms(spec))
        push!(roles, role)
        so = OptT[(nothing, "0"), (fatoms(spec), string(spec[1]))]
        for r in get(residuals, spec[1], Vector{Pair{Symbol,Int}}[])
            push!(so, (r, astr(r)))
        end
        push!(opts, so)
    end
    # Cartesian product with exclusion filter
    n = length(mets); forms = EnzymeFormSpec[]
    for combo in Iterators.product(opts...)
        asf = apf = true; aso = apo = false
        for i in 1:n
            c = combo[i][1]
            if roles[i] == :sub
                c != fulls[i] && (asf = false)
                c !== nothing && (aso = true)
            elseif roles[i] == :prod
                c != fulls[i] && (apf = false)
                c !== nothing && (apo = true)
            end
        end
        ((asf && apo) || (apf && aso)) && continue
        sites = [(metabolite=mets[i],
            atoms=combo[i][1] === nothing ? nothing : copy(combo[i][1]),
            role=roles[i], full_atoms=fulls[i]) for i in 1:n]
        push!(forms, EnzymeFormSpec(
            Symbol("E_" * join((c[2] for c in combo), "_")), sites))
    end
    forms
end

# ─── Catalytic Cycle Enumeration ──────────────────────────────

"""Enumerate catalytic topologies → Vector{MechanismSpec}."""
function _catalytic_topologies(
    forms::Vector{EnzymeFormSpec},
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    fidx, template,
    @nospecialize(reaction::EnzymeReaction);
    max_forms::Int,
)
    sub_ns = [s[1] for s in substrates(reaction)]
    prod_ns = [p[1] for p in products(reaction)]
    free = findfirst(f -> all(s -> s.atoms === nothing, f.sites), forms)
    free === nothing && return MechanismSpec[]
    prod_full = Dict(s.metabolite => s.full_atoms
                     for s in template if s.role == :prod)

    cycles = Set{Int}[]
    seq = [free]; bound = Set{Symbol}(); released = Set{Symbol}()
    residual = Dict{Symbol,Vector{Pair{Symbol,Int}}}()
    lf(s, p, r=residual) = _lookup_form(fidx, template, s, p, r)

    function _try(fi, f)
        fi !== nothing && haskey(adj, minmax(seq[end], fi)) || return
        push!(seq, fi); f(); pop!(seq)
    end

    function _release()
        rem = setdiff(Set(prod_ns), released)
        if isempty(rem)
            seq[end] == free &&
                (s = Set(seq); s ∉ cycles && push!(cycles, s))
            return
        end
        for p in collect(rem)
            delete!(rem, p); push!(released, p)
            _try(isempty(rem) ? free :
                lf(Set{Symbol}(), rem, nothing), _release)
            delete!(released, p); push!(rem, p)
        end
    end

    function _dfs()
        if length(bound) == length(sub_ns) &&
                length(released) == length(prod_ns)
            seq[end] == free &&
                (s = Set(seq); s ∉ cycles && push!(cycles, s))
            return
        end
        for sub in sub_ns  # bind substrate
            sub in bound && continue
            push!(bound, sub)
            _try(lf(bound, Set{Symbol}()), _dfs)
            delete!(bound, sub)
        end
        if length(bound) == length(sub_ns) &&
                length(bound) > length(released)  # isomerize + release
            _try(lf(Set{Symbol}(),
                setdiff(Set(prod_ns), released), nothing), _release)
        end
        isempty(bound) && return  # PP isomerization + release
        cf = forms[seq[end]]
        for prod in prod_ns
            prod in released && continue
            pa = get(prod_full, prod, nothing)
            pa === nothing && continue
            for s in cf.sites
                s.role == :sub && s.atoms !== nothing || continue
                rv = _atom_residual(s.atoms, pa)
                rv === nothing && continue
                old = get(residual, s.metabolite, nothing)
                residual[s.metabolite] = rv
                inter = lf(bound, Set([prod]))
                if inter !== nothing &&
                        haskey(adj, minmax(seq[end], inter))
                    post = lf(bound, Set{Symbol}())
                    if post !== nothing &&
                            haskey(adj, minmax(inter, post))
                        push!(seq, inter, post)
                        push!(released, prod)
                        _dfs()
                        delete!(released, prod)
                        pop!(seq); pop!(seq)
                    end
                end
                old === nothing ? delete!(residual, s.metabolite) :
                    (residual[s.metabolite] = old)
            end
        end
    end
    _dfs()

    # Combine cycles into unions, filter, convert to MechanismSpec
    n = length(cycles)
    n == 0 && return MechanismSpec[]
    combined = unique!([union((cycles[i] for i in 1:n
        if (m >> (i - 1)) & 1 == 1)...) for m in 1:(1 << n) - 1])
    _hr(f) = any(s -> s.role == :sub && s.atoms !== nothing &&
        s.atoms != s.full_atoms, f.sites)
    filter!(combined) do fs
        length(fs) > max_forms && return false
        res = [fi for fi in fs if _hr(forms[fi])]
        isempty(res) && return true
        any(fi -> !any(s -> (s.role == :sub &&
            s.atoms == s.full_atoms) || (s.role == :prod &&
            s.atoms !== nothing), forms[fi].sites), res) &&
        !any(fi -> all(s -> s.role != :sub ||
            s.atoms == s.full_atoms, forms[fi].sites), fs)
    end
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
    OptT = Tuple{Vector{Tuple{Int,Int}}, Set{Tuple{Int,Int}}}
    result = MechanismSpec[]
    for spec in specs
        topo = _used_form_set(spec)
        per_reg = [begin
            opts = OptT[(Tuple{Int,Int}[], Set{Tuple{Int,Int}}())]
            spairs = [(bi, _find_dead_end(fidx, forms[bi], (pos,)))
                for bi in sort(collect(topo))
                for pos in (findfirst(s -> s.metabolite == reg &&
                    s.role == :reg && s.atoms === nothing,
                    forms[bi].sites),) if pos !== nothing]
            filter!(((_, si),) -> si !== nothing, spairs)
            if length(spairs) >= 2
                sm = Dict(spairs)
                mirr = [minmax(sm[a], sm[b]) for (a, b) in spec.edges
                        if haskey(sm, a) && haskey(sm, b)]
                mirrd = Set((a, b) for (a, b) in spec.edges
                            if haskey(sm, a) && haskey(sm, b))
                bind = [minmax(b, s) for (b, s) in spairs]
                push!(opts, ([bind; mirr], Set{Tuple{Int,Int}}()))
                entry = findfirst(((bi, _),) -> all(
                    s -> s.role == :reg || s.atoms === nothing,
                    forms[bi].sites), spairs)
                entry !== nothing && push!(opts,
                    ([minmax(spairs[entry]...); mirr], mirrd))
            end
            opts
        end for reg in reg_names]
        for combo in Iterators.product(per_reg...)
            removals = union((r for (_, r) in combo)...)
            merged = [e for e in spec.edges if e ∉ removals]
            for (add, _) in combo; append!(merged, add); end
            ns = MechanismSpec(reaction, merged,
                fill(false, length(merged)), ParamConstraint[])
            length(_used_form_set(ns)) <= max_forms && push!(result, ns)
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
    topo = _used_form_set(spec)
    budget = max_forms - length(topo)
    budget < 0 && return MechanismSpec[]
    act_pos = Set(k for fi in topo for k in eachindex(forms[fi].sites)
        if forms[fi].sites[k].role == :reg &&
           forms[fi].sites[k].atoms !== nothing)
    pairs = [(fi, k) for fi in sort(collect(topo))
             for k in eachindex(forms[fi].sites)
             if forms[fi].sites[k].role == :reg &&
                forms[fi].sites[k].atoms === nothing && k ∉ act_pos]
    results = MechanismSpec[]
    for mask in 0:((1 << length(pairs)) - 1)
        active = Dict{Int,Vector{Int}}()
        for (idx, (fi, k)) in enumerate(pairs)
            ((mask >> (idx - 1)) & 1) == 1 &&
                push!(get!(active, fi, Int[]), k)
        end
        de = Int[]
        for (fi, positions) in active,
                sm in 1:((1 << length(positions)) - 1)
            pos = Tuple(positions[j] for j in 1:length(positions)
                        if (sm >> (j - 1)) & 1 == 1)
            fi2 = _find_dead_end(fidx, forms[fi], pos)
            fi2 !== nothing && fi2 ∉ topo && fi2 ∉ de && push!(de, fi2)
        end
        length(de) > budget && continue
        existing = Set(minmax(e...) for e in spec.edges)
        all_fi = sort(collect(union(topo, de)))
        new_edges = [(fi, fj) for fi in all_fi for fj in all_fi
            if fi < fj && (fi, fj) ∉ existing &&
               haskey(adj, (fi, fj)) && adj[(fi, fj)].type == :binding]
        push!(results, MechanismSpec(spec.reaction,
            [spec.edges; new_edges],
            fill(false, length(spec.edges) + length(new_edges)),
            ParamConstraint[]))
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
