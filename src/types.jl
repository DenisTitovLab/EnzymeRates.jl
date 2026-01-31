@enum SpeciesRole enzyme metabolite

struct Species
    name::Symbol
    role::SpeciesRole
    atoms::Dict{Symbol,Int}
end

Species(name::Symbol, role::SpeciesRole) = Species(name, role, Dict{Symbol,Int}())

function Base.:(==)(a::Species, b::Species)
    a.name == b.name && a.role == b.role && a.atoms == b.atoms
end

function Base.hash(a::Species, h::UInt)
    hash(a.name, hash(a.role, hash(a.atoms, h)))
end

Base.show(io::IO, s::Species) = print(io, s.name)

struct ReactionSpec
    substrates::Vector{Species}
    products::Vector{Species}
    regulators::Vector{Species}
end

ReactionSpec(s, p) = ReactionSpec(s, p, Species[])

"""
    AbstractEnzymeMechanism

Abstract supertype for enzyme mechanisms. Used as element type in collections
of mechanisms with different type parameters (e.g. from `enumerate_mechanisms`).
"""
abstract type AbstractEnzymeMechanism end

"""
    EnzymeMechanism{N, Steps, FormNames, MetAtoms}

Singleton type encoding an enzyme mechanism entirely in type parameters for
compile-time rate equation generation.

- `N::Int`: number of enzyme forms
- `Steps`: tuple of tuples `(i, j, kf, kr, met_f_or_nothing, met_r_or_nothing)`
- `FormNames`: tuple of `Symbol`s naming the enzyme forms in order
- `MetAtoms`: tuple of `(name, ((atom, count), ...))` for each metabolite
"""
struct EnzymeMechanism{N, Steps, FormNames, MetAtoms} <: AbstractEnzymeMechanism end

"""
    EnzymeMechanism(steps::Vector{Pair{Vector{Species},Vector{Species}}})

Construct an `EnzymeMechanism` from a vector of elementary steps, encoding all
mechanism data into type parameters.
"""
function EnzymeMechanism(steps::Vector{Pair{Vector{Species},Vector{Species}}})
    # Extract enzyme forms in discovery order
    forms = Species[]
    seen = Set{Symbol}()
    for (lhs, rhs) in steps
        for s in vcat(lhs, rhs)
            if s.role == enzyme && s.name ∉ seen
                push!(seen, s.name)
                push!(forms, s)
            end
        end
    end
    n = length(forms)
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))
    form_names = Tuple(s.name for s in forms)

    # Extract metabolites with atoms
    mets = Species[]
    met_seen = Set{Symbol}()
    for (lhs, rhs) in steps
        for s in vcat(lhs, rhs)
            if s.role == metabolite && s.name ∉ met_seen
                push!(met_seen, s.name)
                push!(mets, s)
            end
        end
    end
    met_atoms = Tuple(
        (s.name, Tuple(Tuple(p) for p in sort!(collect(s.atoms); by=first)))
        for s in mets
    )

    # Encode steps
    step_tuples = map(enumerate(steps)) do (step_idx, (lhs, rhs))
        e_lhs = first(s for s in lhs if s.role == enzyme)
        e_rhs = first(s for s in rhs if s.role == enzyme)
        i = name_to_idx[e_lhs.name]
        j = name_to_idx[e_rhs.name]
        m_lhs = [s for s in lhs if s.role == metabolite]
        m_rhs = [s for s in rhs if s.role == metabolite]
        kf = Symbol("k$(step_idx)f")
        kr = Symbol("k$(step_idx)r")
        met_f = isempty(m_lhs) ? nothing : m_lhs[1].name
        met_r = isempty(m_rhs) ? nothing : m_rhs[1].name
        (i, j, kf, kr, met_f, met_r)
    end

    EnzymeMechanism{n, Tuple(step_tuples), form_names, met_atoms}()
end
