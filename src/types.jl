"""
    AbstractEnzymeReaction

Abstract supertype for enzyme reactions with different type parameters.
"""
abstract type AbstractEnzymeReaction end

"""
    EnzymeReaction{Substrates, Products, Regulators}

Singleton type encoding an enzyme reaction specification in type parameters.

Each of `Substrates`, `Products`, `Regulators` is a tuple of
`(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
"""
struct EnzymeReaction{Substrates, Products, Regulators} <: AbstractEnzymeReaction end

function EnzymeReaction(subs::Tuple, prods::Tuple, regs::Tuple=())
    subs_names = Set(s[1] for s in subs)
    prods_names = Set(s[1] for s in prods)
    regs_names = Set(s[1] for s in regs)
    for name in regs_names
        name in subs_names && error("Regulator $(name) also listed as substrate")
        name in prods_names && error("Regulator $(name) also listed as product")
    end
    EnzymeReaction{subs, prods, regs}()
end

"""
    AbstractEnzymeMechanism

Abstract supertype for enzyme mechanisms. Used as element type in collections
of mechanisms with different type parameters.
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

"""
    EnzymeMechanism(species::Tuple, reactions::Tuple)

Construct an `EnzymeMechanism` from explicit species and reaction tuples.

- `species` must be `(substrates, products, regulators, enzymes)` where each entry is a tuple
  of `(name::Symbol, atoms::Tuple{Vararg{Tuple{Symbol,Int}}})`.
- `reactions` is a tuple of `(lhs, rhs)`; each side is a tuple of symbols.
"""
function EnzymeMechanism(species::Tuple, reactions::Tuple)
    length(species) == 4 || error("species must be (substrates, products, regulators, enzymes)")
    subs, prods, regs, enzs = species
    subs isa Tuple || error("substrates must be a tuple of (name, atoms)")
    prods isa Tuple || error("products must be a tuple of (name, atoms)")
    regs isa Tuple || error("regulators must be a tuple of (name, atoms)")
    enzs isa Tuple || error("enzymes must be a tuple of (name, atoms)")

    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    function add_met(name, atoms)
        atoms_dict = Dict{Symbol,Int}(a => c for (a, c) in atoms)
        if haskey(met_atoms, name)
            met_atoms[name] == atoms_dict || error("Inconsistent atoms for metabolite $(name)")
        else
            met_atoms[name] = atoms_dict
        end
    end
    for (name, atoms) in subs
        add_met(name, atoms)
    end
    for (name, atoms) in prods
        add_met(name, atoms)
    end
    for (name, atoms) in regs
        add_met(name, atoms)
    end

    subs_names = [name for (name, _) in subs]
    prods_names = [name for (name, _) in prods]
    regs_names = [name for (name, _) in regs]
    regs_set = Set(regs_names)
    for name in subs_names
        name in regs_set && error("Regulator $(name) also listed as substrate")
    end
    for name in prods_names
        name in regs_set && error("Regulator $(name) also listed as product")
    end

    enzyme_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    for (name, atoms) in enzs
        haskey(enzyme_atoms, name) && error("Duplicate enzyme species $(name)")
        enzyme_atoms[name] = Dict{Symbol,Int}(a => c for (a, c) in atoms)
    end
    isempty(enzyme_atoms) && error("No enzyme species defined")

    for name in keys(enzyme_atoms)
        haskey(met_atoms, name) && error("Species $(name) defined as both enzyme and metabolite")
    end

    expected = Dict{Symbol,Int}()
    for name in subs_names
        expected[name] = get(expected, name, 0) - 1
    end
    for name in prods_names
        expected[name] = get(expected, name, 0) + 1
    end
    for name in regs_names
        expected[name] = get(expected, name, 0) + 0
    end
    for name in regs_names
        get(expected, name, 0) == 0 || error("Regulator $(name) has nonzero net stoichiometry")
    end

    enzyme_set = Set(keys(enzyme_atoms))
    metabolite_set = Set(keys(met_atoms))
    net = Dict{Symbol,Int}()

    for (step_idx, reaction) in enumerate(reactions)
        reaction isa Tuple || error("Reaction $(step_idx) is not a tuple")
        length(reaction) == 2 || error("Reaction $(step_idx) must be (lhs, rhs)")
        lhs, rhs = reaction
        lhs isa Tuple || error("Reaction $(step_idx) lhs must be a tuple")
        rhs isa Tuple || error("Reaction $(step_idx) rhs must be a tuple")

        lhs_enz = 0
        rhs_enz = 0
        lhs_mets = 0
        rhs_mets = 0

        lhs_atoms = Dict{Symbol,Int}()
        rhs_atoms = Dict{Symbol,Int}()

        for s in lhs
            if s in enzyme_set
                lhs_enz += 1
                for (atom, count) in enzyme_atoms[s]
                    lhs_atoms[atom] = get(lhs_atoms, atom, 0) + count
                end
            elseif s in metabolite_set
                lhs_mets += 1
                for (atom, count) in met_atoms[s]
                    lhs_atoms[atom] = get(lhs_atoms, atom, 0) + count
                end
                net[s] = get(net, s, 0) - 1
            else
                error("Reaction $(step_idx) uses unknown species $(s)")
            end
        end
        for s in rhs
            if s in enzyme_set
                rhs_enz += 1
                for (atom, count) in enzyme_atoms[s]
                    rhs_atoms[atom] = get(rhs_atoms, atom, 0) + count
                end
            elseif s in metabolite_set
                rhs_mets += 1
                for (atom, count) in met_atoms[s]
                    rhs_atoms[atom] = get(rhs_atoms, atom, 0) + count
                end
                net[s] = get(net, s, 0) + 1
            else
                error("Reaction $(step_idx) uses unknown species $(s)")
            end
        end

        lhs_enz == 1 || error("Reaction $(step_idx) lhs must contain exactly one enzyme form")
        rhs_enz == 1 || error("Reaction $(step_idx) rhs must contain exactly one enzyme form")
        lhs_mets <= 1 || error("Reaction $(step_idx) lhs has more than one metabolite")
        rhs_mets <= 1 || error("Reaction $(step_idx) rhs has more than one metabolite")

        filter!(p -> p.second != 0, lhs_atoms)
        filter!(p -> p.second != 0, rhs_atoms)
        lhs_atoms == rhs_atoms || error("Atomic conservation failed at step $(step_idx)")
    end

    for (name, coeff) in expected
        net_coeff = get(net, name, 0)
        if coeff == 0
            net_coeff == 0 || error("Regulator $(name) has nonzero net stoichiometry")
        else
            net_coeff == 0 && error("Net stoichiometry mismatch for $(name)")
            sign(net_coeff) == sign(coeff) || error("Net stoichiometry mismatch for $(name)")
            abs(net_coeff) % abs(coeff) == 0 || error("Net stoichiometry mismatch for $(name)")
        end
    end
    for (name, _) in net
        haskey(expected, name) || error("Metabolite $(name) not in species tuple")
    end

    EnzymeMechanism{species, reactions}()
end
