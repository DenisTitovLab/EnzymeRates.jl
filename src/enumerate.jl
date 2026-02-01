using Graphs

"""
    enumerate_mechanisms(::EnzymeReaction, n_params::Int) -> Vector{AbstractEnzymeMechanism}

Enumerate all valid enzyme mechanisms for the given reaction specification
that have exactly `n_params` independent kinetic parameters.
"""
function enumerate_mechanisms(::EnzymeReaction{S, P, R}, n_params::Int) where {S, P, R}
    # Build met_atoms lookup from type params
    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    all_met_names = Symbol[]
    for group in (S, P, R)
        for (name, atoms) in group
            if !haskey(met_atoms, name)
                met_atoms[name] = Dict{Symbol,Int}(a => c for (a, c) in atoms)
            end
            push!(all_met_names, name)
        end
    end
    unique!(all_met_names)

    # Build enzyme forms and track which metabolites are "bound"
    all_enzyme_forms = Symbol[:E]
    form_bound = Dict{Symbol, Vector{Symbol}}(:E => Symbol[])

    # Single-metabolite complexes
    for met in all_met_names
        fn = Symbol("E$(met)")
        push!(all_enzyme_forms, fn)
        form_bound[fn] = [met]
    end

    # Two-metabolite complexes (for multi-substrate)
    if length(all_met_names) > 1
        for i in 1:length(all_met_names)
            for j in (i+1):length(all_met_names)
                m1, m2 = all_met_names[i], all_met_names[j]
                fn = Symbol("E$(m1)$(m2)")
                push!(all_enzyme_forms, fn)
                form_bound[fn] = [m1, m2]
            end
        end
    end

    # Ping-pong forms (modified enzyme carrying atomic fragment)
    ping_pong_forms = _generate_ping_pong_forms(S, P, all_met_names, met_atoms)
    for (form_name, bound_mets) in ping_pong_forms
        if form_name ∉ all_enzyme_forms
            push!(all_enzyme_forms, form_name)
            form_bound[form_name] = bound_mets
        end
    end

    # Generate ALL candidate steps between pairs of enzyme forms
    candidate_steps = Pair{Vector{Symbol}, Vector{Symbol}}[]
    enzyme_set = Set(all_enzyme_forms)

    for i in 1:length(all_enzyme_forms)
        for j in 1:length(all_enzyme_forms)
            i == j && continue
            fi, fj = all_enzyme_forms[i], all_enzyme_forms[j]
            mi = form_bound[fi]
            mj = form_bound[fj]
            mi_names = Set(mi)
            mj_names = Set(mj)

            # Case 1: Simple binding
            if length(mj) == length(mi) + 1
                extra = setdiff(mj_names, mi_names)
                missing_from_j = setdiff(mi_names, mj_names)
                if length(extra) == 1 && isempty(missing_from_j)
                    met_name = first(extra)
                    met_name in keys(met_atoms) && push!(candidate_steps, [fi, met_name] => [fj])
                end
            end

            # Case 2: Catalytic release
            if length(mi) == length(mj) + 1
                lost_names = setdiff(mi_names, mj_names)
                gained_names = setdiff(mj_names, mi_names)
                if length(lost_names) == 1 && isempty(gained_names)
                    lost_name = first(lost_names)

                    # Simple unbinding
                    lost_name in keys(met_atoms) && push!(candidate_steps, [fi] => [fj, lost_name])

                    # Catalytic release: release a different metabolite
                    for alt_met in all_met_names
                        alt_met == lost_name && continue
                        push!(candidate_steps, [fi] => [fj, alt_met])
                    end
                end
            end

            # Case 3: Internal isomerization
            if mi_names == mj_names && !isempty(mi_names) && fi < fj
                push!(candidate_steps, [fi] => [fj])
            end

            # Case 4: Catalytic exchange
            if length(mi) == length(mj)
                lost = setdiff(mi_names, mj_names)
                gained = setdiff(mj_names, mi_names)
                if length(lost) == 1 && length(gained) == 1
                    met_in = first(gained)
                    met_out = first(lost)
                    if met_in in keys(met_atoms) && met_out in keys(met_atoms)
                        push!(candidate_steps, [fi, met_in] => [fj, met_out])
                    end
                end
            end
        end
    end

    unique!(candidate_steps)

    valid_mechanisms = AbstractEnzymeMechanism[]
    _enumerate_subsets!(valid_mechanisms, candidate_steps, S, P, R, met_atoms, enzyme_set, n_params)
    valid_mechanisms
end

"""Generate ping-pong enzyme forms that carry atomic fragments."""
function _generate_ping_pong_forms(subs_t, prods_t, all_met_names, met_atoms)
    forms = Tuple{Symbol, Vector{Symbol}}[]

    for (sub_name, sub_atoms_t) in subs_t
        sub_atoms = met_atoms[sub_name]
        for (prod_name, prod_atoms_t) in prods_t
            prod_atoms = met_atoms[prod_name]
            transferred = Dict{Symbol,Int}()
            for (atom, count) in sub_atoms
                diff = count - get(prod_atoms, atom, 0)
                if diff > 0
                    transferred[atom] = diff
                end
            end
            if !isempty(transferred)
                push!(forms, (:F, Symbol[]))
                for met in all_met_names
                    if met != sub_name && met != prod_name
                        push!(forms, (Symbol("F$(met)"), [met]))
                    end
                end
                break
            end
        end
    end

    forms
end

function _atoms_tuple_from_dict(atoms::Dict{Symbol,Int})
    Tuple((a, c) for (a, c) in sort!(collect(atoms); by=first))
end

function _collect_enzyme_names(raw_steps, enzyme_set)
    names = Symbol[]
    seen = Set{Symbol}()
    for (lhs, rhs) in raw_steps
        for s in vcat(lhs, rhs)
            s in enzyme_set || continue
            if s ∉ seen
                push!(seen, s)
                push!(names, s)
            end
        end
    end
    names
end

function _build_reactions_tuple(raw_steps, enzyme_set)
    rxns = map(raw_steps) do (lhs, rhs)
        e_lhs = first(s for s in lhs if s in enzyme_set)
        e_rhs = first(s for s in rhs if s in enzyme_set)
        m_lhs = [s for s in lhs if s ∉ enzyme_set]
        m_rhs = [s for s in rhs if s ∉ enzyme_set]
        lhs_syms = isempty(m_lhs) ? (e_lhs,) : (e_lhs, m_lhs[1])
        rhs_syms = isempty(m_rhs) ? (e_rhs,) : (e_rhs, m_rhs[1])
        (lhs_syms, rhs_syms)
    end
    return Tuple(rxns)
end

function _infer_enzyme_atoms(enzyme_names, reactions, met_atoms::Dict{Symbol,Dict{Symbol,Int}})
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
                    new_atoms == enzyme_atoms[to] || return nothing
                else
                    enzyme_atoms[to] = new_atoms
                    push!(visited, to)
                    push!(queue, to)
                end
            end
        end
    end

    length(visited) == length(enzyme_names) || return nothing
    return enzyme_atoms
end

"""Enumerate subsets of candidate steps that form valid mechanisms."""
function _enumerate_subsets!(results, candidate_steps, subs_t, prods_t, regs_t, met_atoms, enzyme_set, target_n_params)
    n_steps = length(candidate_steps)
    max_steps = min(n_steps, 8)

    for size in 2:max_steps
        _combinations!(results, candidate_steps, subs_t, prods_t, regs_t, met_atoms, enzyme_set, target_n_params, size, 1, Pair{Vector{Symbol},Vector{Symbol}}[])
    end
end

function _combinations!(results, candidates, subs_t, prods_t, regs_t, met_atoms, enzyme_set, target_n_params, size, start, current)
    if length(current) == size
        _check_and_add!(results, current, subs_t, prods_t, regs_t, met_atoms, enzyme_set, target_n_params)
        return
    end
    for i in start:length(candidates)
        push!(current, candidates[i])
        _combinations!(results, candidates, subs_t, prods_t, regs_t, met_atoms, enzyme_set, target_n_params, size, i + 1, current)
        pop!(current)
    end
end

"""Check validity of raw steps and add to results if valid."""
function _check_and_add!(results, raw_steps, subs_t, prods_t, regs_t, met_atoms, enzyme_set, target_n_params)
    # Extract enzyme forms present in these steps
    enz_names = Symbol[]
    seen = Set{Symbol}()
    for (lhs, rhs) in raw_steps
        for s in vcat(lhs, rhs)
            if s in enzyme_set && s ∉ seen
                push!(seen, s)
                push!(enz_names, s)
            end
        end
    end

    :E in enz_names || return

    # Build graph
    n = length(enz_names)
    name_to_idx = Dict(s => i for (i, s) in enumerate(enz_names))
    g = SimpleDiGraph(n)
    for (lhs, rhs) in raw_steps
        e_lhs = [s for s in lhs if s in enzyme_set]
        e_rhs = [s for s in rhs if s in enzyme_set]
        length(e_lhs) == 1 && length(e_rhs) == 1 || continue
        i = name_to_idx[e_lhs[1]]
        j = name_to_idx[e_rhs[1]]
        add_edge!(g, i, j)
        add_edge!(g, j, i)
    end
    is_connected(g) || return

    # Stoichiometry check
    met_names = Symbol[]
    met_seen = Set{Symbol}()
    for (lhs, rhs) in raw_steps
        for s in vcat(lhs, rhs)
            if s ∉ enzyme_set && s ∉ met_seen
                push!(met_seen, s)
                push!(met_names, s)
            end
        end
    end
    met_idx = Dict(s => i for (i, s) in enumerate(met_names))
    SM = zeros(Int, length(met_names), length(raw_steps))
    for (j, (lhs, rhs)) in enumerate(raw_steps)
        for s in lhs
            s ∉ enzyme_set && (SM[met_idx[s], j] -= 1)
        end
        for s in rhs
            s ∉ enzyme_set && (SM[met_idx[s], j] += 1)
        end
    end
    net = vec(sum(SM, dims=2))

    for (name, _) in subs_t
        idx = findfirst(==(name), met_names)
        idx === nothing && return
        net[idx] < 0 || return
    end
    for (name, _) in prods_t
        idx = findfirst(==(name), met_names)
        idx === nothing && return
        net[idx] > 0 || return
    end
    for (name, _) in regs_t
        idx = findfirst(==(name), met_names)
        idx === nothing && return
        net[idx] == 0 || return
    end

    # Independent params count
    s = length(raw_steps)
    edges_set = Set{Set{Symbol}}()
    for (lhs, rhs) in raw_steps
        e_lhs = first(sp for sp in lhs if sp in enzyme_set)
        e_rhs = first(sp for sp in rhs if sp in enzyme_set)
        push!(edges_set, Set([e_lhs, e_rhs]))
    end
    n_unique_edges = length(edges_set)
    n_parallel_extra = s - n_unique_edges
    n_graph_cycles = n_unique_edges - n + 1
    n_cycles = n_graph_cycles + n_parallel_extra
    if n_cycles == 0
        n_cycles = 1
    end
    n_indep = 2 * s - n_cycles
    n_indep == target_n_params || return

    local_enzyme_set = Set(enz_names)
    rxns = _build_reactions_tuple(raw_steps, local_enzyme_set)
    enzyme_atoms = _infer_enzyme_atoms(enz_names, rxns, met_atoms)
    enzyme_atoms === nothing && return

    enzs_t = Tuple((name, _atoms_tuple_from_dict(enzyme_atoms[name])) for name in enz_names)

    try
        m = EnzymeMechanism((subs_t, prods_t, regs_t, enzs_t), rxns)
        push!(results, m)
    catch
        return
    end
end
