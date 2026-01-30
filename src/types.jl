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
