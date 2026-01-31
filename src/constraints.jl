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
function n_independent_params(m::EnzymeMechanism{N, Steps, FormNames}) where {N, Steps, FormNames}
    s = length(Steps)
    g, _ = graph(m)

    edges_set = Set{Set{Symbol}}()
    for (i, j, kf, kr, met_f, met_r) in Steps
        push!(edges_set, Set([FormNames[i], FormNames[j]]))
    end
    n_unique_edges = length(edges_set)

    n_parallel_extra = s - n_unique_edges
    n_graph_cycles = n_unique_edges - N + 1
    n_cycles = n_graph_cycles + n_parallel_extra

    if n_cycles == 0
        n_cycles = 1
    end

    2 * s - n_cycles
end
