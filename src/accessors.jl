using Graphs

"""Return all distinct enzyme forms in the mechanism."""
function enzyme_forms(m::EnzymeMechanism)
    forms = Species[]
    seen = Set{Symbol}()
    for (lhs, rhs) in m.steps
        for s in vcat(lhs, rhs)
            if s.role == enzyme && s.name ∉ seen
                push!(seen, s.name)
                push!(forms, s)
            end
        end
    end
    forms
end

"""Return all distinct metabolites in the mechanism."""
function metabolites(m::EnzymeMechanism)
    mets = Species[]
    seen = Set{Symbol}()
    for (lhs, rhs) in m.steps
        for s in vcat(lhs, rhs)
            if s.role == metabolite && s.name ∉ seen
                push!(seen, s.name)
                push!(mets, s)
            end
        end
    end
    mets
end

"""Number of distinct enzyme states."""
n_states(m::EnzymeMechanism) = length(enzyme_forms(m))

"""
Build a directed graph of enzyme-form connectivity.
Returns (graph, forms) where forms[i] is the Species for node i.
"""
function graph(m::EnzymeMechanism)
    forms = enzyme_forms(m)
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))
    n = length(forms)
    g = SimpleDiGraph(n)
    for (lhs, rhs) in m.steps
        e_lhs = [s for s in lhs if s.role == enzyme]
        e_rhs = [s for s in rhs if s.role == enzyme]
        length(e_lhs) == 1 && length(e_rhs) == 1 || error("Each step must have exactly one enzyme form on each side")
        i = name_to_idx[e_lhs[1].name]
        j = name_to_idx[e_rhs[1].name]
        add_edge!(g, i, j)
        add_edge!(g, j, i)
    end
    g, forms
end

"""
Stoichiometry matrix: rows = metabolites, columns = steps.
Positive = produced, negative = consumed.
"""
function stoich_matrix(m::EnzymeMechanism)
    mets = metabolites(m)
    met_idx = Dict(s.name => i for (i, s) in enumerate(mets))
    S = zeros(Int, length(mets), length(m.steps))
    for (j, (lhs, rhs)) in enumerate(m.steps)
        for s in lhs
            s.role == metabolite && (S[met_idx[s.name], j] -= 1)
        end
        for s in rhs
            s.role == metabolite && (S[met_idx[s.name], j] += 1)
        end
    end
    S
end

"""
Default parameter grouping: steps that bind/release the same metabolite share parameters.
Returns a vector of vectors of step indices.
"""
function param_groups(m::EnzymeMechanism)
    # Group by the set of metabolites involved in each step
    groups = Dict{Set{Symbol}, Vector{Int}}()
    for (i, (lhs, rhs)) in enumerate(m.steps)
        mets_in_step = Set{Symbol}()
        for s in vcat(lhs, rhs)
            s.role == metabolite && push!(mets_in_step, s.name)
        end
        key = mets_in_step
        push!(get!(groups, key, Int[]), i)
    end
    collect(values(groups))
end

function param_groups(m::EnzymeMechanism, overrides::Dict)
    # Allow custom grouping
    base = param_groups(m)
    # overrides could map step indices to group labels
    base
end
