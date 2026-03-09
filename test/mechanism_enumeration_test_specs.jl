# ABOUTME: Test specifications for mechanism enumeration pipeline
# ABOUTME: Defines spec types, reactions, helper functions, and builder functions

using Random

# ── Helper: convert EnzymeMechanism → MechanismSpec ──────────

"""
Convert a compiled EnzymeMechanism back to a MechanismSpec.
Matches mechanism species to enumerated form indices by walking
the reaction graph using the adjacency from enumerate_enzyme_forms.
"""
function mechanism_spec_from_mechanism(
    m::EnzymeMechanism, @nospecialize(rxn::EnzymeReaction);
    n_catalytic_edges::Int=0)
    site_defs, forms = EnzymeRates.enumerate_enzyme_forms(rxn)
    adj = EnzymeRates._build_adjacency(site_defs, forms)

    mech_forms = EnzymeRates.enzyme_forms(m)
    mf_names = Set(n for (n, _) in mech_forms)

    # Free enzyme (no atoms) → form index 1
    name_to_idx = Dict{Symbol,Int}()
    for (mf_name, mf_atoms) in mech_forms
        if isempty(mf_atoms)
            name_to_idx[mf_name] = 1
        end
    end

    # Iteratively resolve unknown forms via known ones + reactions
    rxns = EnzymeRates.reactions(m)
    for _ in 1:length(rxns)
        for (lhs, rhs) in rxns
            el = [s for s in lhs if s ∈ mf_names]
            er = [s for s in rhs if s ∈ mf_names]
            length(el) == 1 && length(er) == 1 || continue
            # Determine known/unknown sides
            known, unknown = if haskey(name_to_idx, el[1]) &&
                    !haskey(name_to_idx, er[1])
                el[1], er[1]
            elseif haskey(name_to_idx, er[1]) &&
                    !haskey(name_to_idx, el[1])
                er[1], el[1]
            else
                continue
            end
            kidx = name_to_idx[known]
            met = [s for s in Iterators.flatten((lhs, rhs))
                   if s ∉ mf_names]
            met_sym = isempty(met) ? nothing : met[1]
            for j in 1:length(forms)
                j == kidx && continue
                key = minmax(kidx, j)
                haskey(adj, key) || continue
                adj[key] == met_sym || continue
                name_to_idx[unknown] = j
                break
            end
        end
        length(name_to_idx) == length(mech_forms) && break
    end

    length(name_to_idx) == length(mech_forms) || error(
        "Could not resolve all enzyme forms: " *
        "resolved=$(keys(name_to_idx))")

    eq_steps_tuple = EnzymeRates.equilibrium_steps(m)
    pc = EnzymeRates.param_constraints(m)

    edges = Tuple{Int,Int}[]
    for (lhs, rhs) in rxns
        enz_lhs = [s for s in lhs if s ∈ mf_names]
        enz_rhs = [s for s in rhs if s ∈ mf_names]
        length(enz_lhs) == 1 && length(enz_rhs) == 1 ||
            error("Expected exactly 1 enzyme on each side")
        push!(edges, (name_to_idx[enz_lhs[1]],
                      name_to_idx[enz_rhs[1]]))
    end

    eq_steps = collect(Bool, eq_steps_tuple)
    constraints = [
        (t, c, [(s, sc) for (s, sc) in f])
        for (t, c, f) in pc
    ]

    n_cat = n_catalytic_edges > 0 ?
        n_catalytic_edges : length(edges)
    EnzymeRates.MechanismSpec(
        rxn, edges, n_cat, eq_steps,
        constraints, length(parameters(m)))
end
