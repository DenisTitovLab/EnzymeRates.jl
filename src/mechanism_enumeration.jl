# ─── Enzyme Form Enumeration ─────────────────────────────────────────────────
#
# Given an EnzymeReaction with per-metabolite max binding sites (occupancy),
# enumerate all possible enzyme forms. Each binding site is distinguishable.

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

function Base.show(io::IO, f::EnzymeFormSpec)
    print(io, "EnzymeFormSpec(", f.name, ")")
end

"""Total atoms of an enzyme form, computed by summing across all occupied sites."""
function total_atoms(form::EnzymeFormSpec)
    result = Dict{Symbol, Int}()
    for site in form.sites
        site.atoms === nothing && continue
        for (atom, count) in site.atoms
            result[atom] = get(result, atom, 0) + count
        end
    end
    sort([k => v for (k, v) in result]; by=first)
end

# ─── Internal Helpers ────────────────────────────────────────────────────────

"""Template for a binding site: metabolite, index, and full atom content when occupied."""
struct _SiteTemplate
    metabolite::Symbol
    index::Int
    full_atoms::Vector{Pair{Symbol,Int}}
end

"""
Build the ordered list of site templates for a reaction.

Site order:
1. Core sites: 1st site for each substrate (sorted), then 1st site for each product (sorted)
2. Extra sites: 2nd+ sites for substrates (sorted), then 2nd+ sites for products (sorted)
3. Regulator sites: all regulator sites (sorted)
"""
function _build_site_list(S, P, R)
    templates = _SiteTemplate[]

    _make_atoms(spec) = sort([a => c for (a, c) in spec[2]]; by=first)
    _max_sites(spec) = length(spec) >= 3 ? spec[3] : 1

    # Core sites: 1st site per substrate, then 1st site per product
    for s in S
        push!(templates, _SiteTemplate(s[1], 1, _make_atoms(s)))
    end
    for p in P
        push!(templates, _SiteTemplate(p[1], 1, _make_atoms(p)))
    end

    # Extra sites: 2nd+ for substrates, then products
    for s in S
        atoms = _make_atoms(s)
        for i in 2:_max_sites(s)
            push!(templates, _SiteTemplate(s[1], i, atoms))
        end
    end
    for p in P
        atoms = _make_atoms(p)
        for i in 2:_max_sites(p)
            push!(templates, _SiteTemplate(p[1], i, atoms))
        end
    end

    # Regulator sites
    for r in R
        atoms = _make_atoms(r)
        for i in 1:_max_sites(r)
            push!(templates, _SiteTemplate(r[1], i, atoms))
        end
    end

    return templates
end

"""Format atoms as a string like `CX` or `C2N`."""
function _atoms_to_string(atoms::Vector{Pair{Symbol,Int}})
    io = IOBuffer()
    for (sym, count) in sort(atoms; by=first)
        print(io, sym)
        count > 1 && print(io, count)
    end
    String(take!(io))
end

"""Compute the canonical name for an enzyme form from its site contents."""
function _form_name(templates, contents)
    parts = ["E"]
    for (t, content) in zip(templates, contents)
        if content === nothing
            push!(parts, "0")
        elseif content == t.full_atoms
            push!(parts, string(t.metabolite))
        else
            push!(parts, _atoms_to_string(content))
        end
    end
    Symbol(join(parts, "_"))
end

# ─── Standard Form Enumeration ───────────────────────────────────────────────

"""Enumerate all standard forms (each site either empty or fully occupied by its metabolite)."""
function _enumerate_standard_forms(templates, max_total_bound)
    n = length(templates)
    forms = EnzymeFormSpec[]

    for mask in 0:(2^n - 1)
        count_ones(mask) > max_total_bound && continue

        contents = Vector{Union{Nothing, Vector{Pair{Symbol,Int}}}}(undef, n)
        sites = Vector{SiteState}(undef, n)
        for i in 1:n
            occupied = (mask >> (i-1)) & 1 == 1
            c = occupied ? copy(templates[i].full_atoms) : nothing
            contents[i] = c
            sites[i] = SiteState(templates[i].metabolite, templates[i].index, c)
        end

        name = _form_name(templates, contents)
        push!(forms, EnzymeFormSpec(name, sites))
    end

    return forms
end

# ─── Ping-Pong Residual Computation ─────────────────────────────────────────

"""
Compute valid ping-pong residual atom contents for substrate sites.

For each substrate, find partial atom contents that could remain in its site
after releasing one or more products (whose combined atoms are a proper subset
of the substrate's atoms).

Returns `Dict{Symbol, Vector{Vector{Pair{Symbol,Int}}}}` mapping substrate name
to a list of valid residual atom vectors.
"""
function _compute_ping_pong_residuals(S, P)
    residuals = Dict{Symbol, Vector{Vector{Pair{Symbol,Int}}}}()

    for sub_spec in S
        sub_name = sub_spec[1]
        sub_atoms_raw = sub_spec[2]
        isempty(sub_atoms_raw) && continue
        sub_atoms = Dict{Symbol,Int}(a => c for (a, c) in sub_atoms_raw)

        prod_atoms_list = Dict{Symbol,Int}[]
        for prod_spec in P
            prod_atoms_raw = prod_spec[2]
            isempty(prod_atoms_raw) && continue
            push!(prod_atoms_list, Dict{Symbol,Int}(a => c for (a, c) in prod_atoms_raw))
        end
        isempty(prod_atoms_list) && continue

        sub_residuals = Vector{Pair{Symbol,Int}}[]

        for mask in 1:(2^length(prod_atoms_list) - 1)
            combined = Dict{Symbol,Int}()
            for (i, pa) in enumerate(prod_atoms_list)
                if (mask >> (i-1)) & 1 == 1
                    for (atom, count) in pa
                        combined[atom] = get(combined, atom, 0) + count
                    end
                end
            end

            # Check combined atoms are a subset of substrate atoms
            valid = true
            residual = Dict{Symbol,Int}()
            for (atom, count) in sub_atoms
                diff = count - get(combined, atom, 0)
                if diff < 0
                    valid = false
                    break
                elseif diff > 0
                    residual[atom] = diff
                end
            end
            if valid
                for atom in keys(combined)
                    if !haskey(sub_atoms, atom)
                        valid = false
                        break
                    end
                end
            end

            # Residual must be non-empty (not fully released) and differ from full substrate
            if valid && !isempty(residual) && residual != sub_atoms
                r = sort([k => v for (k, v) in residual]; by=first)
                r ∉ sub_residuals && push!(sub_residuals, r)
            end
        end

        if !isempty(sub_residuals)
            residuals[sub_name] = sub_residuals
        end
    end

    return residuals
end

# ─── Ping-Pong Form Enumeration ─────────────────────────────────────────────

"""
Enumerate ping-pong intermediate forms.

For each site that can have a partial (residual) atom content, enumerate all
combinations of standard occupancy on the other sites.
"""
function _enumerate_ping_pong_forms(templates, residuals, max_total_bound)
    n = length(templates)
    forms = EnzymeFormSpec[]
    seen_names = Set{Symbol}()

    for (pp_idx, t) in enumerate(templates)
        haskey(residuals, t.metabolite) || continue

        other_indices = [j for j in 1:n if j != pp_idx]
        n_other = length(other_indices)

        for residual_atoms in residuals[t.metabolite]
            for mask in 0:(2^n_other - 1)
                n_bound = count_ones(mask) + 1  # +1 for the partial site
                n_bound > max_total_bound && continue

                contents = Vector{Union{Nothing, Vector{Pair{Symbol,Int}}}}(undef, n)
                sites = Vector{SiteState}(undef, n)

                contents[pp_idx] = copy(residual_atoms)
                sites[pp_idx] = SiteState(t.metabolite, t.index, copy(residual_atoms))

                for (k, j) in enumerate(other_indices)
                    occupied = (mask >> (k-1)) & 1 == 1
                    c = occupied ? copy(templates[j].full_atoms) : nothing
                    contents[j] = c
                    sites[j] = SiteState(templates[j].metabolite, templates[j].index, c)
                end

                name = _form_name(templates, contents)
                if name ∉ seen_names
                    push!(seen_names, name)
                    push!(forms, EnzymeFormSpec(name, sites))
                end
            end
        end
    end

    return forms
end

# ─── Main API ────────────────────────────────────────────────────────────────

"""
    enumerate_enzyme_forms(reaction::EnzymeReaction; max_total_bound=4, max_ping_pong_intermediates=2)

Enumerate all possible enzyme forms for the given reaction.

Each binding site is distinguishable (site 1 and site 2 for the same metabolite
are distinct). Forms include:

- **Standard forms**: each site independently empty or occupied by its full metabolite
- **Ping-pong intermediates**: sites with partial atom content (residual atoms
  remaining after product release)

# Arguments
- `reaction`: The enzyme reaction specification (with per-metabolite max binding sites)
- `max_total_bound`: Maximum number of simultaneously occupied sites (default 4)
- `max_ping_pong_intermediates`: Maximum number of sites that can simultaneously
  have partial (residual) atom content (default 2). Set to 0 to disable ping-pong
  intermediate enumeration entirely.

# Returns
A `Vector{EnzymeFormSpec}` of all valid enzyme forms.
"""
function enumerate_enzyme_forms(reaction::EnzymeReaction{S,P,R};
    max_total_bound::Int = 4,
    max_ping_pong_intermediates::Int = 2) where {S,P,R}

    templates = _build_site_list(S, P, R)

    # Standard forms
    forms = _enumerate_standard_forms(templates, max_total_bound)
    seen_names = Set{Symbol}(f.name for f in forms)

    # Ping-pong intermediate forms
    if max_ping_pong_intermediates > 0
        has_atoms = any(!isempty(s[2]) for group in (S, P) for s in group)
        if !has_atoms
            n_subs = length(S)
            n_prods = length(P)
            if n_subs >= 2 || n_prods >= 2
                @warn "Reaction has no atom annotations; ping-pong intermediates cannot be determined. " *
                      "Add atom annotations (e.g., A[CX]) to enable ping-pong intermediate enumeration."
            end
        else
            pp_residuals = _compute_ping_pong_residuals(S, P)
            if !isempty(pp_residuals)
                pp_forms = _enumerate_ping_pong_forms(templates, pp_residuals, max_total_bound)
                for f in pp_forms
                    if f.name ∉ seen_names
                        push!(seen_names, f.name)
                        push!(forms, f)
                    end
                end
            end
        end
    end

    return forms
end
