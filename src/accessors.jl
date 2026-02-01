using Graphs

"""Return substrates (with stoichiometric multiplicity)."""
substrates(::EnzymeMechanism{Species}) where {Species} = Species[1]
substrates(::Type{<:EnzymeMechanism{Species}}) where {Species} = Species[1]
substrates(::EnzymeReaction{S,P,R}) where {S,P,R} = S

"""Return products (with stoichiometric multiplicity)."""
products(::EnzymeMechanism{Species}) where {Species} = Species[2]
products(::Type{<:EnzymeMechanism{Species}}) where {Species} = Species[2]
products(::EnzymeReaction{S,P,R}) where {S,P,R} = P

"""Return regulators."""
regulators(::EnzymeMechanism{Species}) where {Species} = Species[3]
regulators(::Type{<:EnzymeMechanism{Species}}) where {Species} = Species[3]
regulators(::EnzymeReaction{S,P,R}) where {S,P,R} = R

"""Return all enzyme forms as a tuple of (name, atoms)."""
enzyme_forms(::EnzymeMechanism{Species}) where {Species} = Species[4]
enzyme_forms(::Type{<:EnzymeMechanism{Species}}) where {Species} = Species[4]

"""Return unique metabolites as a tuple of (name, atoms)."""
function metabolites(::EnzymeMechanism{Species}) where {Species}
    subs, prods, regs = Species[1:3]
    seen = Set{Symbol}()
    mets = Tuple{Symbol,Any}[]
    for group in (subs, prods, regs)
        for (name, atoms) in group
            if name ∉ seen
                push!(seen, name)
                push!(mets, (name, atoms))
            end
        end
    end
    Tuple(mets)
end

"""Number of distinct enzyme states."""
n_states(::EnzymeMechanism{Species}) where {Species} = length(Species[4])

"""Number of steps in the mechanism."""
n_steps(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions} = length(Reactions)

"""Return the reactions tuple directly."""
reactions(::EnzymeMechanism{Species, R}) where {Species, R} = R
reactions(::Type{<:EnzymeMechanism{Species, R}}) where {Species, R} = R

"""
Build a directed graph of enzyme-form connectivity.
Returns (graph, enzyme_forms_tuple).
"""
function graph(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions}
    enzs = Species[4]
    enz_names = Tuple(e[1] for e in enzs)
    name_to_idx = Dict(n => i for (i, n) in enumerate(enz_names))
    enz_set = Set(enz_names)
    g = SimpleDiGraph(length(enzs))
    for (lhs, rhs) in Reactions
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        add_edge!(g, name_to_idx[e_lhs], name_to_idx[e_rhs])
        add_edge!(g, name_to_idx[e_rhs], name_to_idx[e_lhs])
    end
    g, enzs
end

"""
Stoichiometry matrix: rows = metabolites, columns = steps.
Positive = produced, negative = consumed.
"""
function stoich_matrix(m::EnzymeMechanism{Species, Reactions}) where {Species, Reactions}
    mets = metabolites(m)
    met_idx = Dict(m[1] => i for (i, m) in enumerate(mets))
    enz_names = Set(e[1] for e in enzyme_forms(m))
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
function param_groups(m::EnzymeMechanism{Species, Reactions}) where {Species, Reactions}
    groups = Dict{Set{Symbol}, Vector{Int}}()
    enz_names = Set(e[1] for e in enzyme_forms(m))
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

param_groups(m::EnzymeMechanism, overrides::Dict) = param_groups(m)

"""Return rate constant names as a tuple of Symbols, e.g. `(:k1f, :k1r, :k2f, :k2r)`."""
parameters(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions} =
    ntuple(i -> Symbol("k", (i+1)÷2, isodd(i) ? "f" : "r"), 2 * length(Reactions))
