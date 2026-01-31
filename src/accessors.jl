using Graphs

"""Return all distinct enzyme forms in the mechanism."""
function enzyme_forms(::EnzymeMechanism{N, Steps, FormNames, MetAtoms}) where {N, Steps, FormNames, MetAtoms}
    [Species(name, enzyme) for name in FormNames]
end

"""Return all distinct metabolites in the mechanism."""
function metabolites(::EnzymeMechanism{N, Steps, FormNames, MetAtoms}) where {N, Steps, FormNames, MetAtoms}
    [Species(name, metabolite, Dict{Symbol,Int}(atom => count for (atom, count) in atoms))
     for (name, atoms) in MetAtoms]
end

"""Number of distinct enzyme states."""
n_states(::EnzymeMechanism{N}) where {N} = N

"""
Build a directed graph of enzyme-form connectivity.
Returns (graph, forms) where forms[i] is the Species for node i.
"""
function graph(m::EnzymeMechanism{N, Steps}) where {N, Steps}
    forms = enzyme_forms(m)
    g = SimpleDiGraph(N)
    for (i, j, kf, kr, met_f, met_r) in Steps
        add_edge!(g, i, j)
        add_edge!(g, j, i)
    end
    g, forms
end

"""
Stoichiometry matrix: rows = metabolites, columns = steps.
Positive = produced, negative = consumed.
"""
function stoich_matrix(m::EnzymeMechanism{N, Steps, FormNames, MetAtoms}) where {N, Steps, FormNames, MetAtoms}
    mets = metabolites(m)
    met_idx = Dict(s.name => i for (i, s) in enumerate(mets))
    S = zeros(Int, length(mets), length(Steps))
    for (step_j, (i, j, kf, kr, met_f, met_r)) in enumerate(Steps)
        met_f !== nothing && (S[met_idx[met_f], step_j] -= 1)
        met_r !== nothing && (S[met_idx[met_r], step_j] += 1)
    end
    S
end

"""
Default parameter grouping: steps that bind/release the same metabolite share parameters.
Returns a vector of vectors of step indices.
"""
function param_groups(m::EnzymeMechanism{N, Steps}) where {N, Steps}
    groups = Dict{Set{Symbol}, Vector{Int}}()
    for (step_i, (i, j, kf, kr, met_f, met_r)) in enumerate(Steps)
        mets_in_step = Set{Symbol}()
        met_f !== nothing && push!(mets_in_step, met_f)
        met_r !== nothing && push!(mets_in_step, met_r)
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
function steps(::EnzymeMechanism{N, Steps, FormNames, MetAtoms}) where {N, Steps, FormNames, MetAtoms}
    met_species = Dict{Symbol, Species}()
    for (name, atoms) in MetAtoms
        met_species[name] = Species(name, metabolite, Dict{Symbol,Int}(a => c for (a, c) in atoms))
    end
    result = Pair{Vector{Species}, Vector{Species}}[]
    for (i, j, kf, kr, met_f, met_r) in Steps
        lhs = Species[Species(FormNames[i], enzyme)]
        met_f !== nothing && push!(lhs, met_species[met_f])
        rhs = Species[Species(FormNames[j], enzyme)]
        met_r !== nothing && push!(rhs, met_species[met_r])
        push!(result, lhs => rhs)
    end
    result
end

"""Number of steps in the mechanism."""
n_steps(::EnzymeMechanism{N, Steps}) where {N, Steps} = length(Steps)
