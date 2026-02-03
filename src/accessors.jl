using Graphs

"""Return substrates (with stoichiometric multiplicity)."""
substrates(::EnzymeMechanism{Species}) where {Species} = Species[1]
substrates(::EnzymeReaction{S,P,R}) where {S,P,R} = S

"""Return products (with stoichiometric multiplicity)."""
products(::EnzymeMechanism{Species}) where {Species} = Species[2]
products(::EnzymeReaction{S,P,R}) where {S,P,R} = P

"""Return regulators."""
regulators(::EnzymeMechanism{Species}) where {Species} = Species[3]
regulators(::EnzymeReaction{S,P,R}) where {S,P,R} = R

"""Return all enzyme forms as a tuple of (name, atoms)."""
enzyme_forms(::EnzymeMechanism{Species}) where {Species} = Species[4]

"""Compile-time helper: collect unique metabolites from Species type parameter."""
function _unique_metabolites(Species)
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
    return mets
end

"""Return unique metabolites as a tuple of (name, atoms)."""
@generated function metabolites(::EnzymeMechanism{Species}) where {Species}
    return Tuple(_unique_metabolites(Species))
end

"""Number of distinct enzyme states."""
n_states(::EnzymeMechanism{Species}) where {Species} = length(Species[4])

"""Number of steps in the mechanism."""
n_steps(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions} = length(Reactions)

"""Return the reactions tuple directly."""
reactions(::EnzymeMechanism{Species, R}) where {Species, R} = R

"""
Build a directed graph of enzyme-form connectivity.
Returns (graph, enzyme_forms_tuple).
"""
@generated function graph(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions}
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
    return g, enzs
end

"""
Stoichiometry matrix: rows = metabolites, columns = steps.
Positive = produced, negative = consumed.
"""
@generated function stoich_matrix(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions}
    mets = _unique_metabolites(Species)
    met_idx = Dict(m[1] => i for (i, m) in enumerate(mets))
    enz_names = Set(e[1] for e in Species[4])
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
    return S
end

# ─── Mode-dispatched parameters ──────────────────────────────────────────────

"""
    parameters(m::EnzymeMechanism, [mode])

Return the parameter names required for the given mode as a tuple of Symbols.

# Modes
- `IdentifiableHaldaneWegscheider` (default): identifiable combinations + Keq + E_total
- `HaldaneWegscheider`: independent k's + Keq + E_total
- `Raw`: all 2N k's + E_total

# Examples
```julia
parameters(m)                          # identifiable + Keq + E_total
parameters(m, Raw)                     # all k's + E_total
parameters(m, HaldaneWegscheider)      # independent k's + Keq + E_total
```
"""
function parameters end

# Default: IdentifiableHaldaneWegscheider mode
parameters(m::EnzymeMechanism) = parameters(m, IdentifiableHaldaneWegscheider)

# Raw mode: all 2N k-parameters + E_total
@generated function parameters(::EnzymeMechanism{Species, Reactions}, ::RawMode) where {Species, Reactions}
    ks = ntuple(i -> Symbol("k", (i+1)÷2, isodd(i) ? "f" : "r"), 2 * length(Reactions))
    return (ks..., :E_total)
end

# HaldaneWegscheider mode: independent k's + Keq + E_total
@generated function parameters(::M, ::HaldaneWegscheiderMode) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    return (indep..., :Keq, :E_total)
end

# IdentifiableHaldaneWegscheider mode: identifiable combinations + Keq + E_total
# (Defined in identifiability.jl after the required functions are available)

# ─── Legacy API (for backwards compatibility) ────────────────────────────────

"""Return all rate constant names (all 2N k's). Equivalent to `parameters(m, Raw)` without E_total."""
@generated function all_parameters(::EnzymeMechanism{Species, Reactions}) where {Species, Reactions}
    return ntuple(i -> Symbol("k", (i+1)÷2, isodd(i) ? "f" : "r"), 2 * length(Reactions))
end

"""Return only independent parameter names (excludes dependent k's, Keq, E_total)."""
@generated function independent_parameters(::M) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    return indep
end

"""Return dependent parameters as a tuple of `(symbol, expression_string)` pairs."""
@generated function dependent_parameters(::M) where {M <: EnzymeMechanism}
    dep_exprs, _ = _dependent_param_exprs(M)
    pairs = Tuple{Symbol, String}[]
    for (sym, expr) in sort(collect(dep_exprs); by=first)
        push!(pairs, (sym, string(expr)))
    end
    return Tuple(pairs)
end
