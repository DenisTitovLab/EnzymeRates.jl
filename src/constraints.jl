using Graphs

"""
Count independent thermodynamic cycles (Wegscheider conditions).
For a mechanism with N states and S steps, the number of independent cycles
is S - N + 1 (assuming connected graph).
Each cycle gives one Wegscheider constraint.
Additionally, the overall Haldane relation gives one constraint
(relating forward/reverse rate constants to Keq).

Returns the number of independent kinetic parameters.
Each reversible step has 2 rate constants (forward + reverse).
Total raw params = 2 * n_steps.
Constraints = n_cycles (Wegscheider, which includes Haldane).
Independent params = 2 * n_steps - n_cycles.
"""
function n_independent_params(m::EnzymeMechanism)
    n = n_states(m)
    s = length(m.steps)
    # Number of independent cycles in the undirected enzyme-form graph
    # For a connected graph: cycles = edges - nodes + 1
    # But note: multiple steps can connect the same pair of enzyme forms
    # We count unique undirected edges in the enzyme graph
    g, forms = graph(m)
    # The undirected graph has ne(g)/2 edges (since we add both directions)
    # But actually our graph may have parallel edges collapsed by SimpleDiGraph
    # Count unique undirected connections from steps
    edges_set = Set{Set{Symbol}}()
    for (lhs, rhs) in m.steps
        e_lhs = [s for s in lhs if s.role == enzyme][1].name
        e_rhs = [s for s in rhs if s.role == enzyme][1].name
        push!(edges_set, Set([e_lhs, e_rhs]))
    end
    n_unique_edges = length(edges_set)

    # Steps beyond unique edges create additional cycles
    # (parallel edges between same pair of states)
    n_parallel_extra = s - n_unique_edges

    # Graph cycles from topology
    n_graph_cycles = n_unique_edges - n + 1

    # Total independent cycles (each gives one constraint)
    n_cycles = n_graph_cycles + n_parallel_extra

    # But for a simple linear chain (n_cycles = 0), there's still the Haldane relation
    # Actually, Haldane IS one of the cycle constraints when the cycle is the
    # overall reaction loop. For a linear chain with no topological cycles,
    # Haldane gives 1 constraint.
    # For mechanisms with cycles, Haldane is included in the cycle constraints.

    if n_cycles == 0
        # Linear chain: Haldane relation gives 1 constraint
        n_cycles = 1
    end

    2 * s - n_cycles
end
