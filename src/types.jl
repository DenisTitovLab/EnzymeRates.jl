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

struct EnzymeMechanism
    steps::Vector{Pair{Vector{Species},Vector{Species}}}
end

"""
    TypedMechanism{N, Steps}

Type-level encoding of an enzyme mechanism for use with `@generated` functions.
- `N`: number of enzyme forms
- `Steps`: tuple of tuples `(i, j, kf, kr, met_f_or_nothing, met_r_or_nothing)`
"""
struct TypedMechanism{N, Steps} end

"""
    typed_mechanism(m::EnzymeMechanism) -> TypedMechanism{N, Steps}

Convert an `EnzymeMechanism` to a `TypedMechanism` with structure encoded in type parameters.
"""
function typed_mechanism(m::EnzymeMechanism)
    forms = enzyme_forms(m)
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))
    n = length(forms)

    step_tuples = map(enumerate(m.steps)) do (step_idx, (lhs, rhs))
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

    TypedMechanism{n, Tuple(step_tuples)}()
end
