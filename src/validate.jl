"""
Check atomic conservation and regulator net-zero stoichiometry.
"""
function validate(m::EnzymeMechanism)
    subs, prods, regs, enzs = _species_data(m)
    reactions = _reactions_data(m)

    enzyme_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    for (name, atoms) in enzs
        enzyme_atoms[name] = Dict{Symbol,Int}(a => c for (a, c) in atoms)
    end
    enzyme_set = Set(keys(enzyme_atoms))

    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    for (name, atoms) in (subs..., prods..., regs...)
        if !haskey(met_atoms, name)
            met_atoms[name] = Dict{Symbol,Int}(a => c for (a, c) in atoms)
        end
    end

    for (lhs, rhs) in reactions
        lhs_enz = [s for s in lhs if s in enzyme_set]
        rhs_enz = [s for s in rhs if s in enzyme_set]
        length(lhs_enz) == 1 || return false
        length(rhs_enz) == 1 || return false

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
    reg_names = Set(s[1] for s in regs)
    for name in reg_names
        get(net, name, 0) == 0 || return false
    end

    return true
end
