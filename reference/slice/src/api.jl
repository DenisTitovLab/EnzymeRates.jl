# ABOUTME: The parameter API over the joint solve: the parameters a rate law takes, the
# ABOUTME: fitted subset, and every constant computed from the independent values.

"""The name of each column of `sol`, as `constant_names` renders it in its state."""
function column_names(m::Mechanism, sol::Solution)
    active = constant_names(m)
    inactive = is_allosteric(m) ? constant_names(m; inactive = true) : active
    [(c.inactive ? inactive : active)[c.edge][c.slot] for c in sol.cols]
end

"""
    fitted_params(m) -> Tuple

The independent constants that reach the law, in parameter order (active, then inactive;
by edge and slot), then for an allosteric mechanism its regulator constants
(`regulator_names` order) and `:L`.
"""
function fitted_params(m::Mechanism)
    sol = solve(m)
    out = column_names(m, sol)[sol.indep]
    if is_allosteric(m)
        append!(out, (n for (_, ns) in regulator_names(m) for n in ns))
        push!(out, :L)
    end
    Tuple(out)
end

"""The parameters of `m`'s rate law: `fitted_params(m)`, then `:Keq` and `:E_total`."""
parameters(m::Mechanism) = (fitted_params(m)..., :Keq, :E_total)

"""
    all_parameters(m, params) -> NamedTuple

`params` (the values of `parameters(m)`, by name) with every dependent constant added,
computed from the independent values and `Keq` (`x^(p/q)` for a rational exponent): one
entry per name, so every per-step constant of the law reads its value by its name.
"""
function all_parameters(m::Mechanism, params::NamedTuple)
    sol = solve(m)
    names = column_names(m, sol)
    value(d) = prod((pow(params[names[i]], e) for (i, e) in d.terms);
                    init = pow(params.Keq, d.keq))
    merge(params, (; (names[d.col] => value(d) for d in sol.deps
                      if d.label !== :kinetic_group)...))
end
