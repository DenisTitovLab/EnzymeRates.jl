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
Compute the King-Altman rate equation symbolically and return a compiled function.
The function signature is: f(params::NamedTuple, concs::NamedTuple) -> Float64
where params has keys like :k1f, :k1r, :k2f, :k2r, ...
and concs has metabolite concentration keys like :S, :P, ...
Plus :E_total for total enzyme concentration.
"""
function rate_function(m::EnzymeMechanism)
    edges, forms, name_to_idx = _build_rate_info(m)
    n = length(forms)

    # Build list of all directed edges for spanning tree enumeration
    all_dir_edges = Tuple{Int,Int}[]
    for (i, j) in keys(edges)
        push!(all_dir_edges, (i, j))
    end

    # For each root node, enumerate spanning trees toward that root
    # D[root] = sum of products of edge rates for all spanning trees toward root
    tree_products = Dict{Int, Vector{Vector{Tuple{Tuple{Int,Int}, Tuple{Symbol, Union{Symbol,Nothing}}}}}}()

    for root in 1:n
        trees = _spanning_trees_toward(n, all_dir_edges, root)
        tree_products[root] = Vector{Vector{Tuple{Tuple{Int,Int}, Tuple{Symbol, Union{Symbol,Nothing}}}}}()
        for tree in trees
            # For each edge in the tree, get the rate contributions
            # Multiple rate constants can contribute to same directed edge (parallel edges)
            # The combined rate for edge i→j is sum of individual rates
            # But in spanning tree, we use the combined edge rate
            # Actually, King-Altman with combined edges: the rate for edge i→j is
            # the SUM of all individual rate terms for that direction
            # So each spanning tree contributes a product of combined edge rates
            # But expanding: it's a sum of products where each factor picks one
            # individual rate term from each edge.

            # For simplicity with combined edges, we'll handle it in the numeric function
            edge_info = Tuple{Tuple{Int,Int}, Tuple{Symbol, Union{Symbol,Nothing}}}[]
            for (i, j) in tree
                # This edge uses the combined rate for i→j
                push!(edge_info, ((i, j), (:_combined, nothing)))
            end
            push!(tree_products[root], edge_info)
        end
    end

    # Build the rate function as a closure
    # Combined edge rate for (i,j): sum of k*[met] for all rate terms on that edge
    function _eval_edge_rate(i, j, params, concs, edges)
        rate = 0.0
        for (k_name, met_name) in edges[(i, j)]
            r = params[k_name]
            if met_name !== nothing
                r *= concs[met_name]
            end
            rate += r
        end
        rate
    end

    # Precompute the spanning trees structure
    trees_by_root = Dict{Int, Vector{Vector{Tuple{Int,Int}}}}()
    for root in 1:n
        trees_by_root[root] = _spanning_trees_toward(n, all_dir_edges, root)
    end

    # The numerator: use any step to compute net flux
    # v = sum over all steps of (forward_flux - reverse_flux) ... but actually
    # for King-Altman: numerator = E_total * (prod_all_forward - prod_all_reverse)
    # where the products are over all steps in the cycle

    # Collect all forward and reverse rate constants and their metabolites
    fwd_rates = Tuple{Symbol, Union{Symbol,Nothing}}[]  # (k_name, met_or_nothing)
    rev_rates = Tuple{Symbol, Union{Symbol,Nothing}}[]
    for (step_idx, (lhs, rhs)) in enumerate(m.steps)
        m_lhs = [s for s in lhs if s.role == metabolite]
        m_rhs = [s for s in rhs if s.role == metabolite]
        push!(fwd_rates, (Symbol("k$(step_idx)f"), isempty(m_lhs) ? nothing : m_lhs[1].name))
        push!(rev_rates, (Symbol("k$(step_idx)r"), isempty(m_rhs) ? nothing : m_rhs[1].name))
    end

    let edges=edges, trees_by_root=trees_by_root, n=n, fwd_rates=fwd_rates, rev_rates=rev_rates
        function(params, concs)
            E_total = haskey(concs, :E_total) ? concs[:E_total] : 1.0

            # Compute denominator: sum over all roots of sum of spanning tree products
            denom = 0.0
            for root in 1:n
                for tree in trees_by_root[root]
                    prod_val = 1.0
                    for (i, j) in tree
                        # Combined edge rate
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
                    denom += prod_val
                end
            end

            # Compute numerator
            num_fwd = 1.0
            for (k_name, met_name) in fwd_rates
                num_fwd *= params[k_name]
                if met_name !== nothing
                    num_fwd *= concs[met_name]
                end
            end
            num_rev = 1.0
            for (k_name, met_name) in rev_rates
                num_rev *= params[k_name]
                if met_name !== nothing
                    num_rev *= concs[met_name]
                end
            end

            E_total * (num_fwd - num_rev) / denom
        end
    end
end

"""
Return a string representation of the rate equation.
"""
function rate_equation_string(m::EnzymeMechanism)
    edges, forms, name_to_idx = _build_rate_info(m)
    n = length(forms)

    all_dir_edges = collect(keys(edges))

    trees_by_root = Dict{Int, Vector{Vector{Tuple{Int,Int}}}}()
    for root in 1:n
        trees_by_root[root] = _spanning_trees_toward(n, all_dir_edges, root)
    end

    # Build numerator string
    function _rate_term_str(k_name, met_name)
        if met_name === nothing
            return string(k_name)
        else
            return "$(k_name)*$(met_name)"
        end
    end

    # Numerator
    fwd_terms = String[]
    rev_terms = String[]
    for (step_idx, (lhs, rhs)) in enumerate(m.steps)
        m_lhs = [s for s in lhs if s.role == metabolite]
        m_rhs = [s for s in rhs if s.role == metabolite]
        push!(fwd_terms, _rate_term_str(Symbol("k$(step_idx)f"), isempty(m_lhs) ? nothing : m_lhs[1].name))
        push!(rev_terms, _rate_term_str(Symbol("k$(step_idx)r"), isempty(m_rhs) ? nothing : m_rhs[1].name))
    end

    num_str = "E_total * ($(join(fwd_terms, "*")) - $(join(rev_terms, "*")))"

    # Denominator
    denom_parts = String[]
    for root in 1:n
        for tree in trees_by_root[root]
            term_parts = String[]
            for (i, j) in tree
                edge_terms = String[]
                for (k_name, met_name) in edges[(i, j)]
                    push!(edge_terms, _rate_term_str(k_name, met_name))
                end
                if length(edge_terms) == 1
                    push!(term_parts, edge_terms[1])
                else
                    push!(term_parts, "($(join(edge_terms, " + ")))")
                end
            end
            push!(denom_parts, join(term_parts, "*"))
        end
    end

    "($num_str) / ($(join(denom_parts, " + ")))"
end
