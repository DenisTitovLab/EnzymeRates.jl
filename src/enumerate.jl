using Graphs

"""
    enumerate_mechanisms(spec::ReactionSpec, n_params::Int) -> Vector{AbstractEnzymeMechanism}

Enumerate all valid enzyme mechanisms for the given reaction specification
that have exactly `n_params` independent kinetic parameters.
"""
function enumerate_mechanisms(spec::ReactionSpec, n_params::Int)
    all_metabolites = vcat(spec.substrates, spec.products, spec.regulators)

    E = Species(:E, enzyme)

    # Build enzyme forms and track which metabolites are "bound"
    all_enzyme_forms = Species[E]
    form_bound = Dict{Symbol, Vector{Species}}(:E => Species[])

    # Single-metabolite complexes
    for met in all_metabolites
        fn = Symbol("E$(met.name)")
        f = Species(fn, enzyme)
        push!(all_enzyme_forms, f)
        form_bound[fn] = [met]
    end

    # Two-metabolite complexes (for multi-substrate)
    if length(all_metabolites) > 1
        for i in 1:length(all_metabolites)
            for j in (i+1):length(all_metabolites)
                m1, m2 = all_metabolites[i], all_metabolites[j]
                fn = Symbol("E$(m1.name)$(m2.name)")
                f = Species(fn, enzyme)
                push!(all_enzyme_forms, f)
                form_bound[fn] = [m1, m2]
            end
        end
    end

    # Ping-pong forms (modified enzyme carrying atomic fragment)
    ping_pong_forms = _generate_ping_pong_forms(spec, all_metabolites)
    for (form, bound_mets) in ping_pong_forms
        if form.name ∉ [f.name for f in all_enzyme_forms]
            push!(all_enzyme_forms, form)
            form_bound[form.name] = bound_mets
        end
    end

    # Generate ALL candidate steps between pairs of enzyme forms.
    candidate_steps = Pair{Vector{Species}, Vector{Species}}[]

    for i in 1:length(all_enzyme_forms)
        for j in 1:length(all_enzyme_forms)
            i == j && continue
            fi, fj = all_enzyme_forms[i], all_enzyme_forms[j]
            mi = form_bound[fi.name]
            mj = form_bound[fj.name]
            mi_names = Set(s.name for s in mi)
            mj_names = Set(s.name for s in mj)

            # Case 1: Simple binding
            if length(mj) == length(mi) + 1
                extra = setdiff(mj_names, mi_names)
                missing_from_j = setdiff(mi_names, mj_names)
                if length(extra) == 1 && isempty(missing_from_j)
                    met_name = first(extra)
                    met = _find_met(all_metabolites, met_name)
                    met !== nothing && push!(candidate_steps, [fi, met] => [fj])
                end
            end

            # Case 2: Catalytic release
            if length(mi) == length(mj) + 1
                lost_names = setdiff(mi_names, mj_names)
                gained_names = setdiff(mj_names, mi_names)
                if length(lost_names) == 1 && isempty(gained_names)
                    lost_name = first(lost_names)

                    # Simple unbinding
                    met = _find_met(all_metabolites, lost_name)
                    met !== nothing && push!(candidate_steps, [fi] => [fj, met])

                    # Catalytic release: release a different metabolite
                    for alt_met in all_metabolites
                        alt_met.name == lost_name && continue
                        push!(candidate_steps, [fi] => [fj, alt_met])
                    end
                end
            end

            # Case 3: Internal isomerization
            if mi_names == mj_names && !isempty(mi_names) && fi.name < fj.name
                push!(candidate_steps, [fi] => [fj])
            end

            # Case 4: Catalytic exchange
            if length(mi) == length(mj)
                lost = setdiff(mi_names, mj_names)
                gained = setdiff(mj_names, mi_names)
                if length(lost) == 1 && length(gained) == 1
                    met_in = _find_met(all_metabolites, first(gained))
                    met_out = _find_met(all_metabolites, first(lost))
                    if met_in !== nothing && met_out !== nothing
                        push!(candidate_steps, [fi, met_in] => [fj, met_out])
                    end
                end
            end
        end
    end

    # Deduplicate candidate steps
    unique!(candidate_steps)

    # Enumerate valid mechanisms — work on raw steps, only construct EnzymeMechanism for valid ones
    valid_mechanisms = AbstractEnzymeMechanism[]
    _enumerate_subsets!(valid_mechanisms, candidate_steps, spec, n_params)
    valid_mechanisms
end

function _find_met(all_metabolites, name)
    idx = findfirst(s -> s.name == name, all_metabolites)
    idx === nothing ? nothing : all_metabolites[idx]
end

"""Generate ping-pong enzyme forms that carry atomic fragments."""
function _generate_ping_pong_forms(spec::ReactionSpec, all_metabolites::Vector{Species})
    forms = Tuple{Species, Vector{Species}}[]

    for sub in spec.substrates
        for prod in spec.products
            transferred = Dict{Symbol,Int}()
            for (atom, count) in sub.atoms
                diff = count - get(prod.atoms, atom, 0)
                if diff > 0
                    transferred[atom] = diff
                end
            end
            if !isempty(transferred)
                form = Species(:F, enzyme)
                push!(forms, (form, Species[]))
                for met in all_metabolites
                    if met.name != sub.name && met.name != prod.name
                        push!(forms, (Species(Symbol("F$(met.name)"), enzyme), [met]))
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

function _collect_enzyme_names_from_steps(raw_steps)
    names = Symbol[]
    seen = Set{Symbol}()
    for (lhs, rhs) in raw_steps
        for s in vcat(lhs, rhs)
            s.role == enzyme || continue
            if s.name ∉ seen
                push!(seen, s.name)
                push!(names, s.name)
            end
        end
    end
    names
end

function _build_reactions_tuple_from_steps(raw_steps)
    reactions = map(raw_steps) do (lhs, rhs)
        e_lhs = [s for s in lhs if s.role == enzyme][1].name
        e_rhs = [s for s in rhs if s.role == enzyme][1].name
        m_lhs = [s for s in lhs if s.role == metabolite]
        m_rhs = [s for s in rhs if s.role == metabolite]
        lhs_syms = isempty(m_lhs) ? (e_lhs,) : (e_lhs, m_lhs[1].name)
        rhs_syms = isempty(m_rhs) ? (e_rhs,) : (e_rhs, m_rhs[1].name)
        (lhs_syms, rhs_syms)
    end
    return Tuple(reactions)
end

function _infer_enzyme_atoms_from_steps(enzyme_names, reactions, met_atoms::Dict{Symbol,Dict{Symbol,Int}})
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
function _enumerate_subsets!(results, candidate_steps, spec, target_n_params)
    n_steps = length(candidate_steps)
    max_steps = min(n_steps, 8)

    for size in 2:max_steps
        _combinations!(results, candidate_steps, spec, target_n_params, size, 1, Pair{Vector{Species},Vector{Species}}[])
    end
end

function _combinations!(results, candidates, spec, target_n_params, size, start, current)
    if length(current) == size
        _check_and_add_raw!(results, current, spec, target_n_params)
        return
    end
    for i in start:length(candidates)
        push!(current, candidates[i])
        _combinations!(results, candidates, spec, target_n_params, size, i + 1, current)
        pop!(current)
    end
end

"""Check validity of raw steps and add to results if valid."""
function _check_and_add_raw!(results, raw_steps, spec, target_n_params)
    # Extract enzyme forms
    forms = Species[]
    seen = Set{Symbol}()
    for (lhs, rhs) in raw_steps
        for s in vcat(lhs, rhs)
            if s.role == enzyme && s.name ∉ seen
                push!(seen, s.name)
                push!(forms, s)
            end
        end
    end

    any(f -> f.name == :E, forms) || return

    # Build graph
    n = length(forms)
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))
    g = SimpleDiGraph(n)
    for (lhs, rhs) in raw_steps
        e_lhs = [s for s in lhs if s.role == enzyme]
        e_rhs = [s for s in rhs if s.role == enzyme]
        length(e_lhs) == 1 && length(e_rhs) == 1 || continue
        i = name_to_idx[e_lhs[1].name]
        j = name_to_idx[e_rhs[1].name]
        add_edge!(g, i, j)
        add_edge!(g, j, i)
    end
    is_connected(g) || return

    # Stoichiometry check
    mets = Species[]
    met_seen = Set{Symbol}()
    for (lhs, rhs) in raw_steps
        for s in vcat(lhs, rhs)
            if s.role == metabolite && s.name ∉ met_seen
                push!(met_seen, s.name)
                push!(mets, s)
            end
        end
    end
    met_idx = Dict(s.name => i for (i, s) in enumerate(mets))
    SM = zeros(Int, length(mets), length(raw_steps))
    for (j, (lhs, rhs)) in enumerate(raw_steps)
        for s in lhs
            s.role == metabolite && (SM[met_idx[s.name], j] -= 1)
        end
        for s in rhs
            s.role == metabolite && (SM[met_idx[s.name], j] += 1)
        end
    end
    net = vec(sum(SM, dims=2))

    for sub in spec.substrates
        idx = findfirst(s -> s.name == sub.name, mets)
        idx === nothing && return
        net[idx] < 0 || return
    end
    for prod in spec.products
        idx = findfirst(s -> s.name == prod.name, mets)
        idx === nothing && return
        net[idx] > 0 || return
    end
    for reg in spec.regulators
        idx = findfirst(s -> s.name == reg.name, mets)
        idx === nothing && return
        net[idx] == 0 || return
    end

    # Independent params count (on raw steps)
    s = length(raw_steps)
    edges_set = Set{Set{Symbol}}()
    for (lhs, rhs) in raw_steps
        e_lhs = [sp for sp in lhs if sp.role == enzyme][1].name
        e_rhs = [sp for sp in rhs if sp.role == enzyme][1].name
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

    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    for s in spec.substrates
        met_atoms[s.name] = copy(s.atoms)
    end
    for s in spec.products
        met_atoms[s.name] = copy(s.atoms)
    end
    for s in spec.regulators
        met_atoms[s.name] = copy(s.atoms)
    end

    enzyme_names = _collect_enzyme_names_from_steps(raw_steps)
    reactions = _build_reactions_tuple_from_steps(raw_steps)
    enzyme_atoms = _infer_enzyme_atoms_from_steps(enzyme_names, reactions, met_atoms)
    enzyme_atoms === nothing && return

    subs_t = Tuple((s.name, _atoms_tuple_from_dict(s.atoms)) for s in spec.substrates)
    prods_t = Tuple((s.name, _atoms_tuple_from_dict(s.atoms)) for s in spec.products)
    regs_t = Tuple((s.name, _atoms_tuple_from_dict(s.atoms)) for s in spec.regulators)
    enzs_t = Tuple((name, _atoms_tuple_from_dict(enzyme_atoms[name])) for name in enzyme_names)

    try
        m = EnzymeMechanism((subs_t, prods_t, regs_t, enzs_t), reactions)
        push!(results, m)
    catch
        return
    end
end
