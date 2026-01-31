"""
Check atomic conservation for each step of the mechanism.
For each step, sum atoms on LHS must equal sum atoms on RHS.

Enzyme forms carry atoms implicitly: an enzyme form's atomic content
is inferred from the metabolites bound to it. We track this by computing
the "enzyme atom balance" — atoms that must be carried by the enzyme form.
"""
function validate(m::EnzymeMechanism)
    raw = steps(m)
    forms = enzyme_forms(m)
    n = length(forms)
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))

    enzyme_atoms = Dict{Symbol, Dict{Symbol,Int}}()

    root_name = forms[1].name
    enzyme_atoms[root_name] = Dict{Symbol,Int}()

    visited = Set{Symbol}([root_name])
    queue = [root_name]

    while !isempty(queue)
        current = popfirst!(queue)
        for (lhs, rhs) in raw
            e_lhs = [s for s in lhs if s.role == enzyme]
            e_rhs = [s for s in rhs if s.role == enzyme]
            length(e_lhs) == 1 && length(e_rhs) == 1 || continue

            src = e_lhs[1].name
            dst = e_rhs[1].name

            for (from, to, consumed, produced) in [
                (src, dst, lhs, rhs),
                (dst, src, rhs, lhs)
            ]
                if from == current && to ∉ visited
                    new_atoms = copy(enzyme_atoms[from])
                    for s in consumed
                        s.role == metabolite || continue
                        for (atom, count) in s.atoms
                            new_atoms[atom] = get(new_atoms, atom, 0) + count
                        end
                    end
                    for s in produced
                        s.role == metabolite || continue
                        for (atom, count) in s.atoms
                            new_atoms[atom] = get(new_atoms, atom, 0) - count
                        end
                    end
                    filter!(p -> p.second != 0, new_atoms)
                    enzyme_atoms[to] = new_atoms
                    push!(visited, to)
                    push!(queue, to)
                end
            end
        end
    end

    for (lhs, rhs) in raw
        lhs_atoms = Dict{Symbol,Int}()
        for s in lhs
            atoms_to_add = s.role == enzyme ? get(enzyme_atoms, s.name, Dict{Symbol,Int}()) : s.atoms
            for (atom, count) in atoms_to_add
                lhs_atoms[atom] = get(lhs_atoms, atom, 0) + count
            end
        end
        rhs_atoms = Dict{Symbol,Int}()
        for s in rhs
            atoms_to_add = s.role == enzyme ? get(enzyme_atoms, s.name, Dict{Symbol,Int}()) : s.atoms
            for (atom, count) in atoms_to_add
                rhs_atoms[atom] = get(rhs_atoms, atom, 0) + count
            end
        end
        filter!(p -> p.second != 0, lhs_atoms)
        filter!(p -> p.second != 0, rhs_atoms)
        if lhs_atoms != rhs_atoms
            return false
        end
    end

    return true
end
