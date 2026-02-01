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
    EnzymeMechanism{Species,Reactions}

Singleton type encoding an enzyme mechanism in type parameters.

- `Species`: `(substrates, products, regulators, enzyme_species)` where each entry is a tuple
  of `(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
- `Reactions`: tuple of `(lhs, rhs)` where each side is a tuple of species `Symbol`s.
"""
struct EnzymeMechanism{Species, Reactions} <: AbstractEnzymeMechanism end

_atoms_tuple(atoms::Dict{Symbol,Int}) =
    Tuple((a, c) for (a, c) in sort!(collect(atoms); by=first))

function _collect_metabolites(steps)
    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    met_order = Symbol[]
    for (lhs, rhs) in steps
        for s in vcat(lhs, rhs)
            s.role == metabolite || continue
            if !haskey(met_atoms, s.name)
                met_atoms[s.name] = copy(s.atoms)
                push!(met_order, s.name)
            else
                met_atoms[s.name] == s.atoms || error("Inconsistent atoms for metabolite $(s.name)")
            end
        end
    end
    return met_atoms, met_order
end

function _collect_enzyme_names(steps)
    names = Symbol[]
    seen = Set{Symbol}()
    for (lhs, rhs) in steps
        for s in vcat(lhs, rhs)
            s.role == enzyme || continue
            if s.name ∉ seen
                push!(seen, s.name)
                push!(names, s.name)
            end
        end
    end
    return names
end

function _validate_elementary_steps(steps)
    for (idx, (lhs, rhs)) in enumerate(steps)
        lhs_enz = [s for s in lhs if s.role == enzyme]
        rhs_enz = [s for s in rhs if s.role == enzyme]
        length(lhs_enz) == 1 || error("Step $idx: lhs must contain exactly one enzyme form")
        length(rhs_enz) == 1 || error("Step $idx: rhs must contain exactly one enzyme form")

        lhs_mets = [s for s in lhs if s.role == metabolite]
        rhs_mets = [s for s in rhs if s.role == metabolite]
        length(lhs_mets) <= 1 || error("Step $idx: lhs has more than one metabolite")
        length(rhs_mets) <= 1 || error("Step $idx: rhs has more than one metabolite")
    end
end

function _build_reactions_tuple(steps)
    reactions = map(steps) do (lhs, rhs)
        e_lhs = [s for s in lhs if s.role == enzyme][1]
        e_rhs = [s for s in rhs if s.role == enzyme][1]
        m_lhs = [s for s in lhs if s.role == metabolite]
        m_rhs = [s for s in rhs if s.role == metabolite]
        lhs_syms = isempty(m_lhs) ? (e_lhs.name,) : (e_lhs.name, m_lhs[1].name)
        rhs_syms = isempty(m_rhs) ? (e_rhs.name,) : (e_rhs.name, m_rhs[1].name)
        (lhs_syms, rhs_syms)
    end
    return Tuple(reactions)
end

function _net_stoich(reactions, enzyme_set::Set{Symbol})
    net = Dict{Symbol,Int}()
    for (lhs, rhs) in reactions
        for s in lhs
            s in enzyme_set && continue
            net[s] = get(net, s, 0) - 1
        end
        for s in rhs
            s in enzyme_set && continue
            net[s] = get(net, s, 0) + 1
        end
    end
    return net
end

function _spec_metabolites(spec::ReactionSpec)
    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    function add!(s::Species)
        if !haskey(met_atoms, s.name)
            met_atoms[s.name] = copy(s.atoms)
        else
            met_atoms[s.name] == s.atoms || error("Inconsistent atoms for metabolite $(s.name) in spec")
        end
    end
    for s in spec.substrates
        add!(s)
    end
    for s in spec.products
        add!(s)
    end
    for s in spec.regulators
        add!(s)
    end
    return met_atoms
end

function _spec_stoich(spec::ReactionSpec)
    stoich = Dict{Symbol,Int}()
    for s in spec.substrates
        stoich[s.name] = get(stoich, s.name, 0) - 1
    end
    for s in spec.products
        stoich[s.name] = get(stoich, s.name, 0) + 1
    end
    for s in spec.regulators
        stoich[s.name] = get(stoich, s.name, 0) + 0
    end
    return stoich
end

function _infer_enzyme_atoms(enzyme_names::Vector{Symbol}, reactions, met_atoms::Dict{Symbol,Dict{Symbol,Int}})
    root = :E
    root in enzyme_names || error("Free enzyme :E not found in enzyme forms")
    enzyme_set = Set(enzyme_names)

    enzyme_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    enzyme_atoms[root] = Dict{Symbol,Int}()

    visited = Set{Symbol}([root])
    queue = [root]

    while !isempty(queue)
        current = popfirst!(queue)
        for (lhs, rhs) in reactions
            e_lhs = first(s for s in lhs if s in enzyme_set)
            e_rhs = first(s for s in rhs if s in enzyme_set)
            m_lhs = nothing
            m_rhs = nothing
            for s in lhs
                s in enzyme_set && continue
                m_lhs = s
            end
            for s in rhs
                s in enzyme_set && continue
                m_rhs = s
            end

            for (from, to, consumed, produced) in (
                (e_lhs, e_rhs, m_lhs, m_rhs),
                (e_rhs, e_lhs, m_rhs, m_lhs),
            )
                from == current || continue
                new_atoms = copy(enzyme_atoms[from])
                if consumed !== nothing
                    for (atom, count) in met_atoms[consumed]
                        new_atoms[atom] = get(new_atoms, atom, 0) + count
                    end
                end
                if produced !== nothing
                    for (atom, count) in met_atoms[produced]
                        new_atoms[atom] = get(new_atoms, atom, 0) - count
                    end
                end
                filter!(p -> p.second != 0, new_atoms)
                if to in visited
                    new_atoms == enzyme_atoms[to] ||
                        error("Inconsistent enzyme atom assignment for $(to)")
                else
                    enzyme_atoms[to] = new_atoms
                    push!(visited, to)
                    push!(queue, to)
                end
            end
        end
    end

    length(visited) == length(enzyme_names) ||
        error("Enzyme forms are disconnected from :E")
    return enzyme_atoms
end

function _check_atomic_conservation(reactions, enzyme_atoms::Dict{Symbol,Dict{Symbol,Int}}, met_atoms::Dict{Symbol,Dict{Symbol,Int}})
    enzyme_set = Set(keys(enzyme_atoms))
    for (lhs, rhs) in reactions
        lhs_atoms = Dict{Symbol,Int}()
        rhs_atoms = Dict{Symbol,Int}()
        for s in lhs
            atoms_to_add = s in enzyme_set ? enzyme_atoms[s] : met_atoms[s]
            for (atom, count) in atoms_to_add
                lhs_atoms[atom] = get(lhs_atoms, atom, 0) + count
            end
        end
        for s in rhs
            atoms_to_add = s in enzyme_set ? enzyme_atoms[s] : met_atoms[s]
            for (atom, count) in atoms_to_add
                rhs_atoms[atom] = get(rhs_atoms, atom, 0) + count
            end
        end
        filter!(p -> p.second != 0, lhs_atoms)
        filter!(p -> p.second != 0, rhs_atoms)
        lhs_atoms == rhs_atoms || return false
    end
    return true
end

"""
    EnzymeMechanism(steps::Vector{Pair{Vector{Species},Vector{Species}}})

Construct an `EnzymeMechanism` from a vector of elementary steps. Substrates,
products, and regulators are inferred from net stoichiometry.
"""
function EnzymeMechanism(steps::Vector{Pair{Vector{Species},Vector{Species}}})
    _validate_elementary_steps(steps)

    enzyme_names = _collect_enzyme_names(steps)
    met_atoms, met_order = _collect_metabolites(steps)
    reactions = _build_reactions_tuple(steps)

    enzyme_atoms = _infer_enzyme_atoms(enzyme_names, reactions, met_atoms)
    _check_atomic_conservation(reactions, enzyme_atoms, met_atoms) ||
        error("Atomic conservation failed")

    enzyme_set = Set(enzyme_names)
    net = _net_stoich(reactions, enzyme_set)

    subs = Species[]
    prods = Species[]
    regs = Species[]
    for name in met_order
        coeff = get(net, name, 0)
        if coeff < 0
            for _ in 1:(-coeff)
                push!(subs, Species(name, metabolite, met_atoms[name]))
            end
        elseif coeff > 0
            for _ in 1:coeff
                push!(prods, Species(name, metabolite, met_atoms[name]))
            end
        else
            push!(regs, Species(name, metabolite, met_atoms[name]))
        end
    end

    subs_t = Tuple((s.name, _atoms_tuple(s.atoms)) for s in subs)
    prods_t = Tuple((s.name, _atoms_tuple(s.atoms)) for s in prods)
    regs_t = Tuple((s.name, _atoms_tuple(s.atoms)) for s in regs)
    enzs_t = Tuple((name, _atoms_tuple(enzyme_atoms[name])) for name in enzyme_names)

    EnzymeMechanism{(subs_t, prods_t, regs_t, enzs_t), reactions}()
end

"""
    EnzymeMechanism(spec::ReactionSpec, steps::Vector{Pair{Vector{Species},Vector{Species}}})

Construct an `EnzymeMechanism` from a `ReactionSpec` and elementary steps.
"""
function EnzymeMechanism(spec::ReactionSpec, steps::Vector{Pair{Vector{Species},Vector{Species}}})
    _validate_elementary_steps(steps)

    enzyme_names = _collect_enzyme_names(steps)
    spec_met_atoms = _spec_metabolites(spec)
    reactions = _build_reactions_tuple(steps)

    # Ensure all metabolites in steps are in spec and atoms match
    met_atoms, _ = _collect_metabolites(steps)
    for (name, atoms) in met_atoms
        haskey(spec_met_atoms, name) || error("Metabolite $(name) not in ReactionSpec")
        spec_met_atoms[name] == atoms || error("Inconsistent atoms for metabolite $(name)")
    end

    enzyme_atoms = _infer_enzyme_atoms(enzyme_names, reactions, spec_met_atoms)
    _check_atomic_conservation(reactions, enzyme_atoms, spec_met_atoms) ||
        error("Atomic conservation failed")

    enzyme_set = Set(enzyme_names)
    net = _net_stoich(reactions, enzyme_set)
    expected = _spec_stoich(spec)

    for (name, coeff) in expected
        get(net, name, 0) == coeff ||
            error("Net stoichiometry mismatch for $(name)")
    end
    for (name, coeff) in net
        haskey(expected, name) || error("Metabolite $(name) not in ReactionSpec")
        expected[name] == coeff || error("Net stoichiometry mismatch for $(name)")
    end

    subs_t = Tuple((s.name, _atoms_tuple(s.atoms)) for s in spec.substrates)
    prods_t = Tuple((s.name, _atoms_tuple(s.atoms)) for s in spec.products)
    regs_t = Tuple((s.name, _atoms_tuple(s.atoms)) for s in spec.regulators)
    enzs_t = Tuple((name, _atoms_tuple(enzyme_atoms[name])) for name in enzyme_names)

    EnzymeMechanism{(subs_t, prods_t, regs_t, enzs_t), reactions}()
end
