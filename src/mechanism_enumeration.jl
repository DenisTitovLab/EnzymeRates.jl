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
Base.iterate(iter::MechanismIterator, s...) = iterate(iter.inner, s...)

# ─── Edge Classification + Adjacency ─────────────────────────

const EdgeInfo = @NamedTuple{type::Symbol, metabolite::Union{Nothing,Symbol}}

"""Classify edge between two enzyme forms → EdgeInfo or nothing.

Single-site diff → binding (empty→occupied) or release (occupied→empty).
Multi-site diff → isomerization if sub+prod sites change with atom balance."""
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
    # Isomerization: multiple sub/prod sites change with conserved atoms
    all(k -> fa.sites[k].role in (:sub, :prod), diffs) || return nothing
    any(k -> fa.sites[k].role == :sub, diffs) &&
        any(k -> fa.sites[k].role == :prod, diffs) || return nothing
    has_res = any(k -> fa.sites[k].role == :sub && any(
        x -> x !== nothing && x != fa.sites[k].full_atoms,
        (fa.sites[k].atoms, fb.sites[k].atoms)), diffs)
    (has_res || length(diffs) == count(
        s -> s.role in (:sub, :prod), fa.sites)) || return nothing
    # Atom balance: net change across differing sites must be zero
    delta = Dict{Symbol,Int}()
    for k in diffs, (v, sign) in ((fa.sites[k].atoms, 1), (fb.sites[k].atoms, -1))
        v === nothing && continue
        for (a, c) in v; delta[a] = get(delta, a, 0) + sign * c; end
    end
    all(iszero, values(delta)) || return nothing
    (type=:isomerization, metabolite=nothing)
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
    # Ping-pong residuals: subtract product atom combos from substrate atoms
    pa_list = [Dict{Symbol,Int}(a => c for (a, c) in p[2])
               for p in P if !isempty(p[2])]
    residuals = Dict{Symbol,Vector{Vector{Pair{Symbol,Int}}}}()
    for ss in S
        (isempty(ss[2]) || isempty(pa_list)) && continue
        sav = fatoms(ss)
        sa = Dict{Symbol,Int}(sav); sr = Vector{Pair{Symbol,Int}}[]
        for m in 1:(1 << length(pa_list)) - 1
            pa = reduce(mergewith(+), (pa_list[i]
                for i in 1:length(pa_list) if (m >> (i - 1)) & 1 == 1))
            all(get(sa, a, 0) >= c for (a, c) in pa) || continue
            r = sort([a => c - get(pa, a, 0) for (a, c) in sa
                      if c > get(pa, a, 0)]; by=first)
            !isempty(r) && r != sav && r ∉ sr && push!(sr, r)
        end
        !isempty(sr) && (residuals[ss[1]] = sr)
    end
    # Build per-site options
    mets = Symbol[]; fulls = Vector{Pair{Symbol,Int}}[]; roles = Symbol[]
    opts = Vector{Tuple{Union{Nothing,Vector{Pair{Symbol,Int}}},String}}[]
    for (group, role) in ((S, :sub), (P, :prod), (R, :reg)), spec in group
        push!(mets, spec[1]); push!(fulls, fatoms(spec))
        push!(roles, role)
        so = [(nothing, "0"), (fatoms(spec), string(spec[1]))]
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

"""Enumerate catalytic topologies → Vector{MechanismSpec}.

DFS finds elementary cycles through the free enzyme form, then power-set
union generates all valid multi-cycle topologies. Edge direction is computed
inline from the undirected adjacency—no separate directed adjacency needed."""
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
    cat_forms = Set(i for (i, f) in enumerate(forms)
        if all(s -> s.role != :reg || s.atoms === nothing, f.sites))
    cat_adj = Dict((a, b) => info for ((a, b), info) in adj
                   if a ∈ cat_forms && b ∈ cat_forms)
    _has_res(fi) = any(s -> s.role == :sub && s.atoms !== nothing &&
        s.atoms != s.full_atoms, forms[fi].sites)
    _pure_inter(fi) = all(s -> (s.role != :sub || s.atoms != s.full_atoms) &&
        (s.role != :prod || s.atoms === nothing), forms[fi].sites)
    _subs_full(fi) = all(s -> s.role != :sub ||
        s.atoms == s.full_atoms, forms[fi].sites)
    # DFS for elementary cycles through free form. Edge direction is derived
    # inline: (a==cur)==(info.type==:binding) means forward=binding.
    cycles = Set{Set{Int}}()
    function dfs(cur, path, vis, bound, released, mr)
        if cur == free && length(path) > 1 &&
                bound == sub_set && released == prod_set
            push!(cycles, Set(path)); return
        end
        for ((a, b), info) in cat_adj
            nxt = a == cur ? b : b == cur ? a : nothing
            nxt === nothing && continue
            nxt ∈ vis && nxt != free && continue
            et = info.type == :isomerization ? :isomerization :
                ((a == cur) == (info.type == :binding)) ? :binding : :release
            met = info.metabolite
            if et == :binding
                met ∈ sub_set && met ∉ bound && !mr || continue
                push!(bound, met)
            elseif et == :release
                met ∈ prod_set && met ∉ released || continue
                push!(released, met)
            else
                nxt ∉ vis || continue
            end
            push!(path, nxt); push!(vis, nxt)
            dfs(nxt, path, vis, bound, released,
                et == :release ? false : et == :isomerization ? _has_res(nxt) : mr)
            delete!(vis, nxt); pop!(path)
            et == :binding && delete!(bound, met)
            et == :release && delete!(released, met)
        end
    end
    dfs(free, [free], Set([free]), Set{Symbol}(), Set{Symbol}(), false)
    # Combine elementary cycles via power-set union, filter, convert
    n = length(cycles); n == 0 && return MechanismSpec[]
    cvec = collect(cycles)
    combined = unique!([union((cvec[i] for i in 1:n
        if (m >> (i - 1)) & 1 == 1)...) for m in 1:(1 << n) - 1])
    filter!(combined) do fs
        length(fs) > max_forms && return false
        res = [fi for fi in fs if _has_res(fi)]
        isempty(res) && return true
        any(_pure_inter, res) && !any(_subs_full, fs)
    end
    [MechanismSpec(reaction,
        [(a, b) for ((a, b), _) in adj if a ∈ fs && b ∈ fs])
     for fs in combined]
end

# ─── Regulator Expansion (Activators + Dead-Ends) ────────────

"""Generate activator variants: absent, non-essential, or essential per regulator."""
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
        topo = Set(Iterators.flatten(spec.edges))
        per_reg = map(reg_names) do reg
            opts = [(Tuple{Int,Int}[], Set{Tuple{Int,Int}}())]
            pos = findfirst(s -> s.metabolite == reg && s.role == :reg &&
                s.atoms === nothing, forms[first(topo)].sites)
            pos === nothing && return opts
            spairs = [(bi, _find_dead_end(fidx, forms[bi], (pos,)))
                      for bi in sort(collect(topo))]
            filter!(((_, si),) -> si !== nothing, spairs)
            length(spairs) < 2 && return opts
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
            opts
        end
        for combo in Iterators.product(per_reg...)
            removals = union((r for (_, r) in combo)...)
            merged = [e for e in spec.edges if e ∉ removals]
            for (add, _) in combo; append!(merged, add); end
            length(Set(Iterators.flatten(merged))) <= max_forms &&
                push!(result, MechanismSpec(reaction, merged))
        end
    end
    result
end

"""Find dead-end form: base + specified reg positions occupied."""
function _find_dead_end(fidx, base::EnzymeFormSpec, occ_positions)
    key = Tuple(k in occ_positions ? s.full_atoms : s.atoms
                for (k, s) in enumerate(base.sites))
    get(fidx, key, nothing)
end

"""Enumerate dead-end binding configurations.

Precomputes a flat lookup mapping (topo_index, inhibitor_mask) → dead-end form
index, then uses subset filtering via a one-pass set comprehension."""
function _dead_end_configs(
    spec::MechanismSpec,
    adj::Dict{Tuple{Int,Int}, EdgeInfo},
    forms::Vector{EnzymeFormSpec},
    fidx;
    max_forms::Int,
)
    topo = sort(collect(Set(Iterators.flatten(spec.edges))))
    budget = max_forms - length(topo)
    budget < 0 && return MechanismSpec[]
    act_pos = Set(k for fi in topo for (k, s) in enumerate(forms[fi].sites)
                  if s.role == :reg && s.atoms !== nothing)
    inh = [k for (k, s) in enumerate(forms[topo[1]].sites)
           if s.role == :reg && s.atoms === nothing && k ∉ act_pos]
    r = length(inh)
    topo_set = Set(topo)
    existing = Set(minmax(e...) for e in spec.edges)
    # Flat lookup: (topo_index, inhibitor_mask) → dead-end form index
    de_forms = Dict((ti, mask) => fi2
        for (ti, fi) in enumerate(topo) for mask in 1:(1 << r) - 1
        for fi2 in (_find_dead_end(fidx, forms[fi],
            Tuple(inh[j] for j in 1:r if (mask >> (j - 1)) & 1 == 1)),)
        if fi2 !== nothing && fi2 ∉ topo_set)
    results = MechanismSpec[]
    for combo in Iterators.product(
            ntuple(_ -> 0:(1 << r) - 1, length(topo))...)
        de = Set(fi2 for ((ti, m), fi2) in de_forms
                 if m & combo[ti] == m)
        length(de) > budget && continue
        all_forms = union(topo_set, de)
        new_edges = [(a, b) for ((a, b), info) in adj
            if a ∈ all_forms && b ∈ all_forms &&
               (a, b) ∉ existing && info.type == :binding]
        push!(results, MechanismSpec(spec.reaction,
            [spec.edges; new_edges]))
    end
    results
end

# ─── RE/SS + Constraint Lazy Generator ────────────────────────

"""Find equivalent binding step groups from edges."""
function _find_equivalent_groups(edges, adj, forms)
    groups = Dict{Symbol, Vector{Int}}()
    prod_mets = Set(s.metabolite for s in forms[1].sites if s.role == :prod)
    for (i, (a, b)) in enumerate(edges)
        info = get(adj, minmax(a, b), nothing)
        info !== nothing && info.type in (:binding, :release) &&
            info.metabolite ∉ prod_mets &&
            push!(get!(groups, info.metabolite, Int[]), i)
    end
    sort!([sort(v) for v in values(groups) if length(v) >= 2]; by=first)
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
    fidx = Dict(Tuple(s.atoms for s in f.sites) => i
                for (i, f) in enumerate(forms))

    catalytic = _catalytic_topologies(forms, adj, reaction; max_forms)
    stage isa Catalytic && return catalytic

    with_act = _expand_activators(
        catalytic, forms, fidx, reaction; max_forms)
    stage isa WithActivator && return with_act

    de_iter = Iterators.flatmap(
        s -> _dead_end_configs(s, adj, forms, fidx; max_forms), with_act)
    stage isa WithDeadEnd && return de_iter

    # Full: materialize dead-end, compute RE/SS variant count, wrap lazy
    # Count formula: 2^(n - Σgᵢ) × ∏(2^gᵢ + 2) - 2^k
    de_specs = collect(de_iter)
    total = sum(de_specs; init=0) do s
        n = length(s.edges); n == 0 && return 0
        eg = _find_equivalent_groups(s.edges, adj, forms)
        k = length(eg); sum_g = sum(length, eg; init=0)
        (1 << (n - sum_g)) *
            prod(1 << length(g) + 2 for g in eg; init=1) - (1 << k)
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
    used = Set(Iterators.flatten(spec.edges))
    _na(g) = Tuple((n, a) for (n, a, _) in g)
    _atoms(sites) = sort!([a => c for (a, c) in reduce(mergewith(+),
        (Dict(s.atoms) for s in sites if s.atoms !== nothing);
        init=Dict{Symbol,Int}())]; by=first)
    species = (_na(substrates(rxn)), _na(products(rxn)), _na(regulators(rxn)),
        Tuple((forms[i].name, Tuple(Tuple.(_atoms(forms[i].sites))))
            for i in eachindex(forms) if i ∈ used))
    reactions = Tuple(let info = adj[minmax(a, b)]
        ((info.type == :binding ? (forms[a].name, info.metabolite) : (forms[a].name,)),
         (info.type == :release ? (forms[b].name, info.metabolite) : (forms[b].name,)))
    end for (a, b) in spec.edges)
    EnzymeMechanism(species, reactions, Tuple(spec.equilibrium_steps),
        Tuple((t, c, Tuple(Tuple.(f))) for (t, c, f) in spec.param_constraints))
end
