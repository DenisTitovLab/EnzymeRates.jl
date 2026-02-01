using Graphs

_atoms_dict(atoms) = Dict{Symbol,Int}(atom => count for (atom, count) in atoms)

function _species_from_group(group, role)
    [Species(name, role, _atoms_dict(atoms)) for (name, atoms) in group]
end

"""Return substrates (with stoichiometric multiplicity) for the mechanism."""
function substrates(::EnzymeMechanism{SpeciesT}) where {SpeciesT}
    _species_from_group(SpeciesT[1], metabolite)
end

"""Return products (with stoichiometric multiplicity) for the mechanism."""
function products(::EnzymeMechanism{SpeciesT}) where {SpeciesT}
    _species_from_group(SpeciesT[2], metabolite)
end

"""Return regulators (with multiplicity as defined) for the mechanism."""
function regulators(::EnzymeMechanism{SpeciesT}) where {SpeciesT}
    _species_from_group(SpeciesT[3], metabolite)
end

"""Return all distinct enzyme forms in the mechanism."""
function enzyme_forms(::EnzymeMechanism{SpeciesT}) where {SpeciesT}
    enzs = SpeciesT[4]
    [Species(name, enzyme, _atoms_dict(atoms)) for (name, atoms) in enzs]
end

"""Return all distinct metabolites in the mechanism."""
function metabolites(::EnzymeMechanism{SpeciesT}) where {SpeciesT}
    subs, prods, regs = SpeciesT[1:3]
    seen = Set{Symbol}()
    mets = Species[]
    for group in (subs, prods, regs)
        for (name, atoms) in group
            if name ∉ seen
                push!(seen, name)
                push!(mets, Species(name, metabolite, _atoms_dict(atoms)))
            end
        end
    end
    mets
end

"""Number of distinct enzyme states."""
function n_states(::EnzymeMechanism{SpeciesT}) where {SpeciesT}
    length(SpeciesT[4])
end

"""
Build a directed graph of enzyme-form connectivity.
Returns (graph, forms) where forms[i] is the Species for node i.
"""
function graph(m::EnzymeMechanism{SpeciesT, Reactions}) where {SpeciesT, Reactions}
    forms = enzyme_forms(m)
    enz_names = Tuple(s.name for s in forms)
    name_to_idx = Dict(n => i for (i, n) in enumerate(enz_names))
    enz_set = Set(enz_names)
    g = SimpleDiGraph(length(forms))
    for (lhs, rhs) in Reactions
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        add_edge!(g, name_to_idx[e_lhs], name_to_idx[e_rhs])
        add_edge!(g, name_to_idx[e_rhs], name_to_idx[e_lhs])
    end
    g, forms
end

"""
Stoichiometry matrix: rows = metabolites, columns = steps.
Positive = produced, negative = consumed.
"""
function stoich_matrix(m::EnzymeMechanism{SpeciesT, Reactions}) where {SpeciesT, Reactions}
    mets = metabolites(m)
    met_idx = Dict(s.name => i for (i, s) in enumerate(mets))
    enz_names = Set(s.name for s in enzyme_forms(m))
    S = zeros(Int, length(mets), length(Reactions))
    for (step_j, (lhs, rhs)) in enumerate(Reactions)
        for s in lhs
            s in enz_names && continue
            S[met_idx[s], step_j] -= 1
        end
        for s in rhs
            s in enz_names && continue
            S[met_idx[s], step_j] += 1
        end
    end
    S
end

"""
Default parameter grouping: steps that bind/release the same metabolite share parameters.
Returns a vector of vectors of step indices.
"""
function param_groups(m::EnzymeMechanism{SpeciesT, Reactions}) where {SpeciesT, Reactions}
    groups = Dict{Set{Symbol}, Vector{Int}}()
    enz_names = Set(s.name for s in enzyme_forms(m))
    for (step_i, (lhs, rhs)) in enumerate(Reactions)
        mets_in_step = Set{Symbol}()
        for s in lhs
            s in enz_names || push!(mets_in_step, s)
        end
        for s in rhs
            s in enz_names || push!(mets_in_step, s)
        end
        push!(get!(groups, mets_in_step, Int[]), step_i)
    end
    collect(values(groups))
end

function param_groups(m::EnzymeMechanism, overrides::Dict)
    param_groups(m)
end

"""
Reconstruct the raw steps as `Vector{Pair{Vector{Species}, Vector{Species}}}`.
"""
function steps(::EnzymeMechanism{SpeciesT, Reactions}) where {SpeciesT, Reactions}
    subs, prods, regs, enzs = SpeciesT

    met_species = Dict{Symbol, Species}()
    for (name, atoms) in (subs..., prods..., regs...)
        if !haskey(met_species, name)
            met_species[name] = Species(name, metabolite, _atoms_dict(atoms))
        end
    end

    enz_species = Dict{Symbol, Species}()
    for (name, atoms) in enzs
        enz_species[name] = Species(name, enzyme, _atoms_dict(atoms))
    end

    result = Pair{Vector{Species}, Vector{Species}}[]
    for (lhs, rhs) in Reactions
        lhs_vec = Species[]
        rhs_vec = Species[]
        for s in lhs
            push!(lhs_vec, haskey(enz_species, s) ? enz_species[s] : met_species[s])
        end
        for s in rhs
            push!(rhs_vec, haskey(enz_species, s) ? enz_species[s] : met_species[s])
        end
        push!(result, lhs_vec => rhs_vec)
    end
    result
end

"""Number of steps in the mechanism."""
function n_steps(::EnzymeMechanism{SpeciesT, Reactions}) where {SpeciesT, Reactions}
    length(Reactions)
end
