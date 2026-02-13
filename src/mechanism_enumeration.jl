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
    mets = Symbol[]; idxs = Int[]; fulls = Vector{Pair{Symbol,Int}}[]
    roles = Symbol[]  # :sub, :prod, :reg
    opts = Vector{Tuple{Union{Nothing, Vector{Pair{Symbol,Int}}}, String}}[]
    function add!(met, idx, full, role)
        push!(mets, met); push!(idxs, idx); push!(fulls, full); push!(roles, role)
        site_opts = Tuple{Union{Nothing, Vector{Pair{Symbol,Int}}}, String}[
            (nothing, "0"), (full, string(met))]
        for r in get(residuals, met, Vector{Pair{Symbol,Int}}[])
            push!(site_opts, (r, astr(r)))
        end
        push!(opts, site_opts)
    end
    for s in S; add!(s[1], 1, fatoms(s), :sub); end
    for p in P; add!(p[1], 1, fatoms(p), :prod); end
    for s in S; for i in 2:nsites(s); add!(s[1], i, fatoms(s), :sub); end; end
    for p in P; for i in 2:nsites(p); add!(p[1], i, fatoms(p), :prod); end; end
    for r in R; for i in 1:nsites(r); add!(r[1], i, fatoms(r), :reg); end; end

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
        all_sub_full = true; all_prod_full = true
        any_sub_occ = false; any_prod_occ = false
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
        push!(forms, EnzymeFormSpec(Symbol("E_" * join(name_parts, "_")), sites))
    end

    return forms
end
