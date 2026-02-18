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
MechanismSpec(reaction, edges) =
    MechanismSpec(reaction, edges, fill(false, length(edges)), ParamConstraint[])

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
_form_atoms(sites) = sort([a => c for (a, c) in _sum_atoms(sites)]; by=first)

"""Compute residual atoms after subtracting `to_remove` from `sub_atoms`.
Returns sorted Vector{Pair} or nothing if subtraction is invalid/trivial."""
function _atom_residual(sub_atoms, to_remove)
    sa, pa = Dict{Symbol,Int}(sub_atoms), Dict{Symbol,Int}(to_remove)
    all(get(sa, a, 0) >= c for (a, c) in pa) || return nothing
    r = sort([a => c - get(pa, a, 0) for (a, c) in sa
              if c > get(pa, a, 0)]; by=first)
    isempty(r) || r == sub_atoms ? nothing : r
end

"""Classify edge between two enzyme forms → EdgeInfo or nothing."""
function _classify_edge(fa::EnzymeFormSpec, fb::EnzymeFormSpec)
    diffs = [k for k in eachindex(fa.sites)
             if fa.sites[k].atoms != fb.sites[k].atoms]
    isempty(diffs) && return nothing
    if length(diffs) == 1
        k = diffs[1]; sa, sb = fa.sites[k], fb.sites[k]
        if sa.atoms === nothing && sb.atoms !== nothing
            sa.role in (:sub, :prod) && sb.atoms != sa.full_atoms &&
                return nothing
            return (type=:binding, metabolite=sa.metabolite)
        end
        if sb.atoms === nothing && sa.atoms !== nothing
            met = sa.atoms == sa.full_atoms || sa.role == :reg ?
                sa.metabolite :
                let i = findfirst(s -> s.role == :prod &&
                        s.full_atoms == sa.atoms, fa.sites)
                    i !== nothing ? fa.sites[i].metabolite : nothing end
            return met !== nothing ? (type=:release, metabolite=met) : nothing
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
        _sum_atoms(fa.sites, (:sub, :prod)) ==
            _sum_atoms(fb.sites, (:sub, :prod)) &&
        return (type=:isomerization, metabolite=nothing)
    nothing
end

"""Build adjacency dict from enzyme forms."""
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
        asf = all(i -> roles[i] != :sub || combo[i][1] == fulls[i], 1:n)
        apf = all(i -> roles[i] != :prod || combo[i][1] == fulls[i], 1:n)
        aso = any(i -> roles[i] == :sub && combo[i][1] !== nothing, 1:n)
        apo = any(i -> roles[i] == :prod && combo[i][1] !== nothing, 1:n)
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

"""Build directed adjacency list from undirected edge dict.

Reverses binding↔release for the back direction; isomerization is symmetric."""
function _build_directed_adj(adj)
    nbrs = Dict{Int, Vector{Tuple{Int, Symbol, Union{Nothing, Symbol}}}}()
    for ((a, b), info) in adj
        push!(get!(nbrs, a, []), (b, info.type, info.metabolite))
        rtype = info.type == :binding ? :release :
                info.type == :release ? :binding : :isomerization
        push!(get!(nbrs, b, []), (a, rtype, info.metabolite))
    end
    nbrs
end

"""Enumerate catalytic topologies → Vector{MechanismSpec}.

Uses a graph-walk DFS over the directed adjacency list. Substrate binding,
product release, and isomerization (including ping-pong) are all handled
uniformly as edge traversals—no special-case code per reaction type."""
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
    # Restrict to forms with all regulator sites empty (catalytic topology)
    cat_forms = Set(i for (i, f) in enumerate(forms)
        if all(s -> s.role != :reg || s.atoms === nothing, f.sites))
    nbrs = _build_directed_adj(adj)

    cycles = Set{Set{Int}}()
    path = [free]; visited = Set{Int}([free])
    bound = Set{Symbol}(); released = Set{Symbol}()
    # After PP isomerization (producing a product + residual), force release before
    # any further binding. This prevents non-ping-pong interleaving of bind/release.
    must_release = false

    function step!(next)
        push!(path, next); push!(visited, next)
        dfs()
        delete!(visited, next); pop!(path)
    end
    function dfs()
        cur = path[end]
        if cur == free && length(path) > 1 &&
                bound == sub_set && released == prod_set
            push!(cycles, Set(path))
            return
        end
        for (next, etype, met) in get(nbrs, cur, ())
            next ∈ cat_forms || continue
            next ∈ visited && next != free && continue
            if etype == :binding && met ∈ sub_set && met ∉ bound &&
                    !must_release
                push!(bound, met); step!(next); delete!(bound, met)
            elseif etype == :release && met ∈ prod_set && met ∉ released
                old_mr = must_release; must_release = false
                push!(released, met); step!(next); delete!(released, met)
                must_release = old_mr
            elseif etype == :isomerization && next ∉ visited
                has_res = any(s -> s.role == :sub && s.atoms !== nothing &&
                    s.atoms != s.full_atoms, forms[next].sites)
                old_mr = must_release; must_release = has_res
                step!(next)
                must_release = old_mr
            end
        end
    end
    dfs()

    # Combine elementary cycles via power-set union, filter, convert
    n = length(cycles)
    n == 0 && return MechanismSpec[]
    cvec = collect(cycles)
    combined = unique!([union((cvec[i] for i in 1:n
        if (m >> (i - 1)) & 1 == 1)...) for m in 1:(1 << n) - 1])
    _has_residual(f) = any(s -> s.role == :sub && s.atoms !== nothing &&
        s.atoms != s.full_atoms, f.sites)
    filter!(combined) do fs
        length(fs) > max_forms && return false
        res = [fi for fi in fs if _has_residual(forms[fi])]
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
        MechanismSpec(reaction, edges)
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
            pos = findfirst(s -> s.metabolite == reg && s.role == :reg &&
                s.atoms === nothing, forms[first(topo)].sites)
            spairs = pos === nothing ? typeof((0, 0))[] :
                [(bi, _find_dead_end(fidx, forms[bi], (pos,)))
                 for bi in sort(collect(topo))]
            filter!(((_, si),) -> si !== nothing, spairs)
            if length(spairs) >= 2
                sm = Dict(spairs)
                me = [(a, b) for (a, b) in spec.edges
                      if haskey(sm, a) && haskey(sm, b)]
                mirr = [minmax(sm[a], sm[b]) for (a, b) in me]
                bind = [minmax(b, s) for (b, s) in spairs]
                push!(opts, ([bind; mirr], Set{Tuple{Int,Int}}()))
                entry = findfirst(((bi, _),) -> all(
                    s -> s.role == :reg || s.atoms === nothing,
                    forms[bi].sites), spairs)
                entry !== nothing && push!(opts,
                    ([minmax(spairs[entry]...); mirr], Set(me)))
            end
            opts
        end for reg in reg_names]
        for combo in Iterators.product(per_reg...)
            removals = union((r for (_, r) in combo)...)
            merged = [e for e in spec.edges if e ∉ removals]
            for (add, _) in combo; append!(merged, add); end
            ns = MechanismSpec(reaction, merged)
            length(_used_form_set(ns)) <= max_forms && push!(result, ns)
        end
    end
    result
end

# ─── Dead-End Configs ─────────────────────────────────────────

"""Find dead-end form: base + specified reg positions occupied."""
function _find_dead_end(fidx, base::EnzymeFormSpec, occ_positions)
    key = Tuple(k in occ_positions ? s.full_atoms : s.atoms
                for (k, s) in enumerate(base.sites))
    get(fidx, key, nothing)
end

"""Enumerate dead-end binding configurations.

Precomputes a lookup mapping (topology_form, inhibitor_mask) → dead-end form
index, then uses subset filtering instead of online subset enumeration."""
function _dead_end_configs(
    spec::MechanismSpec,
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    forms::Vector{EnzymeFormSpec},
    fidx;
    max_forms::Int,
)
    sorted_topo = sort(collect(_used_form_set(spec)))
    budget = max_forms - length(sorted_topo)
    budget < 0 && return MechanismSpec[]
    act_pos = Set(k for fi in sorted_topo for (k, s) in enumerate(forms[fi].sites)
                  if s.role == :reg && s.atoms !== nothing)
    inh = [k for (k, s) in enumerate(forms[sorted_topo[1]].sites)
           if s.role == :reg && s.atoms === nothing && k ∉ act_pos]
    r = length(inh)
    topo_set = Set(sorted_topo)
    existing = Set(minmax(e...) for e in spec.edges)
    # Precompute: for each topology form, map inhibitor mask → dead-end form
    de_lookup = [Dict(mask => fi2
        for mask in 1:(1 << r) - 1
        for fi2 in (_find_dead_end(fidx, forms[fi],
            Tuple(inh[j] for j in 1:r if (mask >> (j - 1)) & 1 == 1)),)
        if fi2 !== nothing && fi2 ∉ topo_set) for fi in sorted_topo]
    results = MechanismSpec[]
    for combo in Iterators.product(
            ntuple(_ -> 0:(1 << r) - 1, length(sorted_topo))...)
        de = Set{Int}()
        for (ti, fi_mask) in enumerate(combo)
            for (m, fi2) in de_lookup[ti]
                m & fi_mask == m && push!(de, fi2)
            end
        end
        length(de) > budget && continue
        all_fi = sort([sorted_topo; collect(de)])
        new_edges = [(fi, fj) for fi in all_fi for fj in all_fi
            if fi < fj && (fi, fj) ∉ existing &&
               haskey(adj, (fi, fj)) && adj[(fi, fj)].type == :binding]
        push!(results, MechanismSpec(spec.reaction,
            [spec.edges; new_edges]))
    end
    results
end

# ─── RE/SS + Constraint Lazy Generator ────────────────────────

"""Find equivalent binding step groups from edges."""
function _find_equivalent_groups(edges, adj, forms)
    groups = Dict{Symbol, Vector{Int}}()
    for (i, (a, b)) in enumerate(edges)
        info = get(adj, minmax(a, b), nothing)
        info === nothing && continue
        info.type in (:binding, :release) || continue
        any(s -> s.metabolite == info.metabolite && s.role == :prod,
            forms[1].sites) && continue
        push!(get!(groups, info.metabolite, Int[]), i)
    end
    sort!([sort(v) for v in values(groups) if length(v) >= 2]; by=first)
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
    edges = spec.edges; n = length(edges)
    equiv_groups = _find_equivalent_groups(edges, adj, forms)
    Iterators.flatmap(0:((1 << n) - 2)) do re_mask
        eq_steps = Bool[(re_mask >> (i-1)) & 1 == 1 for i in 1:n]
        vg = [g for g in equiv_groups
            if all(eq_steps[s] == eq_steps[g[1]] for s in g)]

        Iterators.map(0:((1 << length(vg)) - 1)) do cmask
            pc = ParamConstraint[]
            for (gi, g) in enumerate(vg)
                ((cmask >> (gi-1)) & 1) == 0 && continue
                p, ss = eq_steps[g[1]] ? ("K", ("",)) : ("k", ("f", "r"))
                for j in 2:length(g), s in ss
                    push!(pc, (Symbol("$p$(g[j])$s"), 1,
                               [(Symbol("$p$(g[1])$s"), 1)]))
                end
            end
            MechanismSpec(spec.reaction, edges, eq_steps, pc)
        end
    end
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

    catalytic = _catalytic_topologies(forms, adj, reaction; max_forms)
    stage isa Catalytic && return catalytic

    with_act = _expand_activators(
        catalytic, forms, fidx, reaction; max_forms)
    stage isa WithActivator && return with_act

    de_iter = Iterators.flatmap(
        s -> _dead_end_configs(s, adj, forms, fidx; max_forms), with_act)
    stage isa WithDeadEnd && return de_iter

    # Full: materialize dead-end, compute count, wrap lazy RE/SS
    de_specs = collect(de_iter)
    total = sum(de_specs; init=0) do s
        eg = _find_equivalent_groups(s.edges, adj, forms)
        _count_ress_variants(length(s.edges), eg)
    end
    inner = Iterators.flatmap(
        s -> _ress_variants(s, adj, forms), de_specs)
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
    _na(g) = Tuple((n, a) for (n, a, _) in g)
    species = (_na(substrates(rxn)), _na(products(rxn)), _na(regulators(rxn)),
        Tuple((forms[i].name, Tuple(Tuple.(_form_atoms(forms[i].sites))))
            for i in 1:length(forms) if i in used))
    reactions = Tuple(let i = adj[minmax(a, b)]
        ((i.type == :binding ? (forms[a].name, i.metabolite) : (forms[a].name,)),
         (i.type == :release ? (forms[b].name, i.metabolite) : (forms[b].name,)))
    end for (a, b) in spec.edges)
    eq_steps = Tuple(spec.equilibrium_steps)
    constraints = Tuple(
        (target, coeff, Tuple(Tuple.(factors)))
        for (target, coeff, factors) in spec.param_constraints)
    EnzymeMechanism(species, reactions, eq_steps, constraints)
end
