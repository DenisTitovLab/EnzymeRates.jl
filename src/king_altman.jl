using LinearAlgebra

"""
Build the King-Altman rate matrix for the mechanism.
Returns (rate_matrix_expr, forms, met_names) where rate_matrix_expr[i,j]
gives the pseudo-first-order rate from state i to state j as a symbolic expression.

Each step has forward rate k{step_idx}f and reverse rate k{step_idx}r.
Binding steps multiply by metabolite concentration.
"""
function _build_rate_info(m::EnzymeMechanism)
    forms = enzyme_forms(m)
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))
    n = length(forms)

    # Collect all directed edges with their rate expressions
    # edges[i][j] = list of (rate_constant_name, metabolite_name_or_nothing)
    edges = Dict{Tuple{Int,Int}, Vector{Tuple{Symbol,Union{Symbol,Nothing}}}}()

    for (step_idx, (lhs, rhs)) in enumerate(m.steps)
        e_lhs = [s for s in lhs if s.role == enzyme][1]
        e_rhs = [s for s in rhs if s.role == enzyme][1]
        i = name_to_idx[e_lhs.name]
        j = name_to_idx[e_rhs.name]

        m_lhs = [s for s in lhs if s.role == metabolite]
        m_rhs = [s for s in rhs if s.role == metabolite]

        # Forward: i → j, rate = k{step_idx}f * [metabolites on lhs]
        kf = Symbol("k$(step_idx)f")
        met_f = isempty(m_lhs) ? nothing : m_lhs[1].name
        push!(get!(edges, (i, j), []), (kf, met_f))

        # Reverse: j → i, rate = k{step_idx}r * [metabolites on rhs]
        kr = Symbol("k$(step_idx)r")
        met_r = isempty(m_rhs) ? nothing : m_rhs[1].name
        push!(get!(edges, (j, i), []), (kr, met_r))
    end

    return edges, forms, name_to_idx
end

"""
Enumerate all spanning trees of an undirected graph rooted at a given node,
directed toward that root. Returns list of vectors of (from, to) edges.
Uses Kirchhoff's theorem approach via direct enumeration.
"""
function _spanning_trees_toward(n::Int, all_edges::Vector{Tuple{Int,Int}}, root::Int)
    # We need spanning trees of the complete undirected graph
    # directed toward root. Each tree has n-1 edges, one from each non-root node.
    # For each non-root node, exactly one outgoing edge.

    # Build adjacency: for each node, what edges go out from it (directed toward root means
    # each non-root node picks one neighbor to point to)
    # Actually: a spanning arborescence rooted at `root` is a set of n-1 directed edges
    # such that every node except root has exactly one outgoing edge,
    # and following edges from any node leads to root.

    non_root = [i for i in 1:n if i != root]
    # For each non-root node, list possible targets
    adj = Dict{Int, Vector{Int}}()
    for node in non_root
        adj[node] = Int[]
    end
    for (i, j) in all_edges
        if i != root && i in keys(adj)
            push!(adj[i], j)
        end
    end
    # Remove duplicates
    for node in non_root
        adj[node] = unique(adj[node])
    end

    # Enumerate by recursion: pick edge for each non-root node
    trees = Vector{Vector{Tuple{Int,Int}}}()
    _enum_trees!(trees, non_root, adj, root, 1, Tuple{Int,Int}[], n)
    trees
end

function _enum_trees!(result, non_root, adj, root, idx, current_edges, n)
    if idx > length(non_root)
        # Check if it's a valid arborescence: all nodes reach root
        if _reaches_root(current_edges, non_root, root, n)
            push!(result, copy(current_edges))
        end
        return
    end
    node = non_root[idx]
    for target in adj[node]
        push!(current_edges, (node, target))
        _enum_trees!(result, non_root, adj, root, idx + 1, current_edges, n)
        pop!(current_edges)
    end
end

function _reaches_root(edges, non_root, root, n)
    # Build successor map
    succ = Dict{Int,Int}()
    for (i, j) in edges
        succ[i] = j
    end
    for node in non_root
        visited = Set{Int}()
        cur = node
        while cur != root
            cur in visited && return false
            push!(visited, cur)
            haskey(succ, cur) || return false
            cur = succ[cur]
        end
    end
    true
end

"""
Precompute all King-Altman data needed by both `rate_function` and `rate_equation_string`.

Returns `(edges, trees_by_root, step_info, n)` where `step_info` is a vector of
`(i, j, kf_name, kr_name, met_lhs, met_rhs)` tuples describing each mechanism step
(enzyme form indices i→j, rate constant symbols, and metabolite names or nothing).
"""
function _king_altman_data(m::EnzymeMechanism)
    edges, forms, name_to_idx = _build_rate_info(m)
    n = length(forms)

    all_dir_edges = collect(keys(edges))

    trees_by_root = Dict{Int, Vector{Vector{Tuple{Int,Int}}}}()
    for root in 1:n
        trees_by_root[root] = _spanning_trees_toward(n, all_dir_edges, root)
    end

    # Collect per-step info: enzyme form indices + rate constants + metabolites
    step_info = Tuple{Int, Int, Symbol, Symbol, Union{Symbol,Nothing}, Union{Symbol,Nothing}}[]
    for (step_idx, (lhs, rhs)) in enumerate(m.steps)
        e_lhs = [s for s in lhs if s.role == enzyme][1]
        e_rhs = [s for s in rhs if s.role == enzyme][1]
        i = name_to_idx[e_lhs.name]
        j = name_to_idx[e_rhs.name]
        m_lhs = [s for s in lhs if s.role == metabolite]
        m_rhs = [s for s in rhs if s.role == metabolite]
        kf = Symbol("k$(step_idx)f")
        kr = Symbol("k$(step_idx)r")
        met_f = isempty(m_lhs) ? nothing : m_lhs[1].name
        met_r = isempty(m_rhs) ? nothing : m_rhs[1].name
        push!(step_info, (i, j, kf, kr, met_f, met_r))
    end

    return edges, trees_by_root, step_info, n
end

function _rate_term_str(k_name, met_name)
    met_name === nothing ? string(k_name) : "$(k_name)*$(met_name)"
end

"""
Compute the King-Altman rate equation symbolically and return a compiled function.
The function signature is: f(params::NamedTuple, concs::NamedTuple) -> Float64
where params has keys like :k1f, :k1r, :k2f, :k2r, ...
and concs has metabolite concentration keys like :S, :P, ...
Plus :E_total for total enzyme concentration.
"""
function rate_function(m::EnzymeMechanism)
    edges, trees_by_root, step_info, n = _king_altman_data(m)

    let edges=edges, trees_by_root=trees_by_root, n=n, step_info=step_info
        function(params, concs)
            E_total = haskey(concs, :E_total) ? concs[:E_total] : 1.0

            # Compute D[root] for each root (cofactors)
            D = zeros(n)
            for root in 1:n
                for tree in trees_by_root[root]
                    prod_val = 1.0
                    for (i, j) in tree
                        edge_rate = 0.0
                        for (k_name, met_name) in edges[(i, j)]
                            r = params[k_name]
                            if met_name !== nothing
                                r *= concs[met_name]
                            end
                            edge_rate += r
                        end
                        prod_val *= edge_rate
                    end
                    D[root] += prod_val
                end
            end

            D_total = sum(D)

            # Numerator: flux through step 1 (net consumption of first substrate)
            # v = E_total * (R_fwd * D[i] - R_rev * D[j]) / D_total
            i, j, kf, kr, met_f, met_r = step_info[1]
            rf = Float64(params[kf])
            if met_f !== nothing
                rf *= concs[met_f]
            end
            rr = Float64(params[kr])
            if met_r !== nothing
                rr *= concs[met_r]
            end

            E_total * (rf * D[i] - rr * D[j]) / D_total
        end
    end
end

"""
Return a string representation of the rate equation.
"""
function _denom_tree_str(trees, edges)
    parts = String[]
    for tree in trees
        term_parts = String[]
        for (i, j) in tree
            edge_terms = [_rate_term_str(k, met) for (k, met) in edges[(i, j)]]
            if length(edge_terms) == 1
                push!(term_parts, edge_terms[1])
            else
                push!(term_parts, "($(join(edge_terms, " + ")))")
            end
        end
        push!(parts, join(term_parts, "*"))
    end
    parts
end

function rate_equation_string(m::EnzymeMechanism)
    edges, trees_by_root, step_info, n = _king_altman_data(m)

    # Build D[root] strings
    D_strs = String[]
    for root in 1:n
        parts = _denom_tree_str(trees_by_root[root], edges)
        push!(D_strs, join(parts, " + "))
    end

    # Numerator: flux through step 1
    i, j, kf, kr, met_f, met_r = step_info[1]
    fwd_str = _rate_term_str(kf, met_f)
    rev_str = _rate_term_str(kr, met_r)
    Di_str = D_strs[i]
    Dj_str = D_strs[j]
    num_str = "E_total * ($fwd_str * ($Di_str) - $rev_str * ($Dj_str))"

    # Denominator: sum of all D[root]
    all_denom_parts = String[]
    for root in 1:n
        append!(all_denom_parts, _denom_tree_str(trees_by_root[root], edges))
    end

    "($num_str) / ($(join(all_denom_parts, " + ")))"
end
