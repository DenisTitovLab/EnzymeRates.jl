# ABOUTME: Helpers shared by the engine-stage check scripts: law values, a package sample's
# ABOUTME: independents, units of a law's variables, and evaluation of printed laws.
include(joinpath(@__DIR__, "..", "solve", "lib.jl"))

"""`v/E_total` of law `L` at concentrations `conc` (by name) and parameters `pars` (by name),
in BigFloat."""
function lawval(L, conc, pars)
    vals = [L.vars.conc[i] ? BigFloat(conc[n]) : BigFloat(getfield(pars, n))
            for (i, n) in enumerate(L.vars.names)]
    S.evaluate(L.H, vals) * S.evaluate(L.Q, vals) / S.evaluate(L.D, vals)
end

"""BigFloat values of every package name at sample `s`: fitted, dependents from their
monomials, Pass-2 folds."""
function pkg_values(p, s)
    v = Dict{Symbol, BigFloat}(n => BigFloat(x) for (n, x) in s.params)
    for (d, mono) in p.dependent_monomials
        v[d] = prod((v[x]^BigFloat(e) for (x, e) in mono); init = big(1))
    end
    foreach(((a, t),) -> v[a] = v[t], p.pass2_rename)
    v
end

"""My parameters (independents, Keq, E_total = 1) read from the package's per-step values
`vals` of record `r` (`pkg_values`)."""
function my_params(r, sol, vals)
    m = r.c.mech
    names = S.column_names(m, sol)
    x = S.edge_constants(r.c, [BigFloat[vals[n] for n in ns] for ns in r.p.step_names])
    (; (names[i] => x[sol.cols[i].edge][sol.cols[i].slot] for i in sol.indep)...,
     Keq = vals[:Keq], E_total = big(1))
end

"""Units `(concentration, time)` of each variable of law `L` of `m` under solve `sol`."""
function var_units(m, L, sol)
    fwd, names = S.orient(m), S.column_names(m, sol)
    u = Dict{Int, Tuple{Rational{Int}, Rational{Int}}}()
    for (i, n) in enumerate(L.vars.names)
        L.vars.conc[i] && (u[i] = (1, 0))
        n === :Keq && (u[i] = (length(m.products) - length(m.substrates), 0))
    end
    for c in sol.indep
        col, o = sol.cols[c], S.oriented(m, sol.cols[c].edge, fwd)
        u[L.vars.index[names[c]]] = !o.re ?
            (-(col.slot == 1 ? S.nibsum(o.cons) : S.nibsum(o.rel)), -1) :
            S.is_binding(o) ? (1, 0) : (S.nibsum(o.rel) - S.nibsum(o.cons), 0)
    end
    u
end

"""Units of a monomial under the variable units `u`."""
mono_units(mo, u) = reduce((a, b) -> a .+ b, (e .* u[v] for (v, e) in mo);
                           init = (0 // 1, 0 // 1))

"""`v` of a printed law at `pars` and `conc` (NamedTuples): every line run in order, except
that the coefficient form's `v` line (`legend_first`) runs after its legend."""
function eval_text(text, pars, conc; legend_first = false)
    mod = Module()
    Core.eval(mod, :(params = $pars))
    Core.eval(mod, :(concs = $conc))
    foreach(x -> Core.eval(mod, :($x = $(getfield(pars, x)))), keys(pars))
    foreach(x -> Core.eval(mod, :($x = $(getfield(conc, x)))), keys(conc))
    lines = split(text, "\n")
    legend_first && (lines = [lines[2:end]; lines[1]])
    for l in lines
        startswith(l, "#") || Core.eval(mod, Meta.parse(l))
    end
    Core.eval(mod, :v)
end
