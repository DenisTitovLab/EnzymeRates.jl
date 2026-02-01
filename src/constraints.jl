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
    reactions = _reactions_data(m)
    enz_names = Set(s.name for s in enzyme_forms(m))
    n = length(enz_names)
    s = length(reactions)

    edges_set = Set{Set{Symbol}}()
    for (lhs, rhs) in reactions
        e_lhs = first(s for s in lhs if s in enz_names)
        e_rhs = first(s for s in rhs if s in enz_names)
        push!(edges_set, Set([e_lhs, e_rhs]))
    end
    n_unique_edges = length(edges_set)

    n_parallel_extra = s - n_unique_edges
    n_graph_cycles = n_unique_edges - n + 1
    n_cycles = n_graph_cycles + n_parallel_extra

    if n_cycles == 0
        n_cycles = 1
    end

    2 * s - n_cycles
end
