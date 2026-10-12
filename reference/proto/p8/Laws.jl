# ABOUTME: Final rate laws of Pop.PMech mechanisms as Canon.Law values (probe P8): the redesign's
# ABOUTME: (Engine law + Thermo joint solve; P5b's MWC model) and the package's own exported laws.
const _P8_COMMON = joinpath(@__DIR__, "..", "common")
isdefined(Main, :Pop) || include(joinpath(_P8_COMMON, "Population.jl"))
isdefined(Main, :Poly) || include(joinpath(_P8_COMMON, "Poly.jl"))
isdefined(Main, :Engine) || include(joinpath(_P8_COMMON, "Engine.jl"))
isdefined(Main, :Thermo) || include(joinpath(_P8_COMMON, "Thermo.jl"))
isdefined(Main, :Oracle) || include(joinpath(_P8_COMMON, "Oracle.jl"))
isdefined(Main, :Allo) || include(joinpath(@__DIR__, "vendor", "Allo_p5b.jl"))
isdefined(Main, :Canon) || include(joinpath(@__DIR__, "Canon.jl"))

module Laws

import ..Pop, ..Poly, ..Engine, ..Thermo, ..Oracle, ..Allo, ..Canon

"""A final law: `law` over the fitted names and `Keq` (dependents substituted), `alllaw` over every
name with folds applied but each other dependent kept as a symbol, `defs` those dependents'
definitions (`name => [symbol => exponent]`, `:Keq` included), the kind (`:plain` v = num/den,
or `:mwc`), the catalytic and site multiplicities of an MWC law, and what its numeric checks need
(the Thermo solution of a plain redesign law, the P5b model of an MWC redesign law)."""
struct Built
    law::Canon.Law
    alllaw::Canon.Law
    defs::Vector{Pair{Symbol, Vector{Pair{Symbol, Rational{Int}}}}}
    kind::Symbol
    n::Int
    mults::Vector{Int}
    extra::Any
end

"""Whether a definition is a fold: one other name with exponent 1 and no `Keq`."""
isfold(rhs) = length(rhs) == 1 && rhs[1].first !== :Keq && rhs[1].second == 1

"""The law with every fold of `defs` applied and every other dependent kept as a parameter, and
those dependents' definitions."""
function all_names_law(law::Canon.Law, defs)
    folds = Dict(d => rhs for (d, rhs) in defs if isfold(rhs))
    keep = [d => rhs for (d, rhs) in defs if !isfold(rhs)]
    comps = [substitute(collect(p), s -> get(folds, s, nothing)) for p in law.comps]
    present = Set(s for c in comps for m in keys(c) for (s, _) in m)
    used = [d => rhs for (d, rhs) in keep if d in present]
    params = sort!(unique([law.params; first.(used)]))
    Canon.Law(comps, law.fixed, law.units, params,
              string(law.header, "|I=", length(law.params))), used
end

toint(q) = (isinteger(q) || error("fractional exponent $q"); Int(q))

"""Terms `monomial => coefficient` (monomials as `symbol => exponent`) with each symbol `s` for
which `sub(s)` returns `symbol => exponent` pairs replaced by them, to the term's power."""
function substitute(terms, sub)
    out = Canon.Poly()
    for (m, c) in terms
        d = Dict{Symbol, Int}()
        for (s, e) in m
            r = sub(s)
            if r === nothing
                d[s] = get(d, s, 0) + e
            else
                for (s2, e2) in r
                    d[s2] = get(d, s2, 0) + toint(e2 * e)
                end
            end
        end
        k = sort!([s => e for (s, e) in d if e != 0]; by = first)
        v = get(out, k, big(0) // 1) + Rational{BigInt}(c)
        iszero(v) ? delete!(out, k) : (out[k] = v)
    end
    out
end

"""The terms of an Engine/Poly polynomial with variables as their names."""
symterms(p::Poly.P, V::Poly.Vars) = [[V.names[v] => e for (v, e) in mo] => c for (mo, c) in p.t]

"""Error unless every name of `law` is fixed or a declared parameter."""
function check_names(law::Canon.Law, id)
    extra = setdiff(Set(s for c in law.comps for m in keys(c) for (s, _) in m),
                    law.fixed, law.params)
    isempty(extra) || error("$id: names outside the parameters: $(collect(extra))")
    law
end

# ─── Plain laws ────────────────────────────────────────────────────────

"""The redesign's plain law: the Engine's bare law with every step constant written through
Thermo's joint solve as `Keq^q · Π independent^e` (folds included); parameters are the
content-keyed independent column names."""
function plain_redesign(m::Pop.PMech)
    L = Engine.bare_rate_law(m; check = false)
    sol = Thermo.solve(m)
    cname(i) = Symbol(sol.cols[i].name)
    expr = Dict{Int, Vector{Pair{Symbol, Rational{Int}}}}(i => [cname(i) => 1 // 1]
                                                          for i in sol.indep)
    for d in sol.deps
        expr[d.col] = [[cname(f) => e for (f, e) in d.terms];
                       d.keq == 0 ? Pair{Symbol, Rational{Int}}[] : [:Keq => d.keq]]
    end
    at = Dict((c.j, c.r) => i for (i, c) in enumerate(sol.cols))
    sub = Dict{Symbol, Vector{Pair{Symbol, Rational{Int}}}}()
    for (j, vs) in enumerate(L.stepvars)
        c = sol.cs[j]
        name(v) = L.vars.names[v]
        if c.re
            e = c.flipped && !c.bindK ? -1 : 1
            sub[name(vs[1])] = [s => x * e for (s, x) in expr[at[(j, 1)]]]
        else
            sub[name(vs[1])] = expr[at[(j, c.flipped ? 2 : 1)]]
            sub[name(vs[2])] = expr[at[(j, c.flipped ? 1 : 2)]]
        end
    end
    f(s) = get(sub, s, nothing)
    num, den = (substitute(symterms(p, L.vars), f) for p in (L.num, L.den))
    law = Canon.Law([num, den], Set([m.metabolites; :Keq]), [[1, 2]],
                    sort!([cname(i) for i in sol.indep]), "plain")
    # the same law over every column name: each step constant is its column (to the RE flip)
    defs = [cname(d.col) => [[cname(f) => e for (f, e) in d.terms];
                             d.keq == 0 ? Pair{Symbol, Rational{Int}}[] : [:Keq => d.keq]]
            for d in sol.deps]
    raw = Dict{Symbol, Vector{Pair{Symbol, Rational{Int}}}}()
    for (j, vs) in enumerate(L.stepvars)
        c = sol.cs[j]
        if c.re
            raw[L.vars.names[vs[1]]] = [cname(at[(j, 1)]) => (c.flipped && !c.bindK ? -1 : 1) // 1]
        else
            raw[L.vars.names[vs[1]]] = [cname(at[(j, c.flipped ? 2 : 1)]) => 1 // 1]
            raw[L.vars.names[vs[2]]] = [cname(at[(j, c.flipped ? 1 : 2)]) => 1 // 1]
        end
    end
    rawlaw = Canon.Law([substitute(symterms(p, L.vars), s -> get(raw, s, nothing))
                        for p in (L.num, L.den)], law.fixed, law.units, law.params, "plain")
    alllaw, used = all_names_law(rawlaw, defs)
    Built(check_names(law, m.id), check_names(alllaw, m.id), used, :plain, 1, Int[], sol)
end

"""The package's plain law: its Reduced polynomials (Pass-2 folds applied) with every
dependent replaced by its monomial in the fitted names and `Keq`."""
function plain_today(m::Pop.PMech)
    deps = Dict(m.dependent_monomials)
    f(s) = get(deps, s, nothing)
    law = Canon.Law([substitute(m.pkg_num, f), substitute(m.pkg_den, f)],
                    Set([m.metabolites; :Keq]), [[1, 2]], copy(m.fitted_params), "plain")
    raw = Canon.Law([substitute(m.pkg_num, s -> nothing), substitute(m.pkg_den, s -> nothing)],
                    law.fixed, law.units, copy(m.fitted_params), "plain")
    alllaw, used = all_names_law(raw, m.dependent_monomials)
    Built(check_names(law, m.id), check_names(alllaw, m.id), used, :plain, 1, Int[], nothing)
end

# ─── MWC laws ──────────────────────────────────────────────────────────

"""Regulator-site factor `1 + Σ [lig]/K` of every site in each conformation, as terms, with
the site multiplicities; `regname(st, lig, tag)` names a ligand's constant or returns nothing."""
function site_terms(m::Pop.PMech, regname)
    comps, mults = Vector{Vector{Pair{Vector{Pair{Symbol, Int}}, Rational{Int}}}}(), Int[]
    for site in m.sites
        for st in (:A, :I)
            t = [Pair{Symbol, Int}[] => 1 // 1]
            for (lig, tag) in zip(site.ligands, site.states)
                K = regname(st, lig, tag)
                K === nothing || push!(t, sort!([lig => 1, K => -1]; by = first) => 1 // 1)
            end
            push!(comps, t)
        end
        push!(mults, site.multiplicity)
    end
    comps, mults
end

"""An MWC law from each state's `(num, den, dfree)` terms and the site terms: components
`[N_A, d_A, D_A, N_I, d_I, D_I, W_A1, W_I1, …]`, each state one normalization unit and each
site factor its own; `L` is fixed (it is never renamed)."""
function mwc_law(m::Pop.PMech, stateA, stateI, sites, mults, params, f)
    comps = [substitute(p, f) for p in [stateA[1], stateA[3], stateA[2],
                                        stateI[1], stateI[3], stateI[2], sites...]]
    units = [[[1, 2, 3], [4, 5, 6]]; [[i] for i in 7:length(comps)]]
    Canon.Law(comps, Set([m.metabolites; m.regulators; :Keq; :L]), units,
              filter(!=(:L), params), string("mwc|n=", m.cat_multiplicity, "|sites=", mults))
end

"""The package's name for a regulator constant (`_reg_K`)."""
pkg_regname(st, lig, tag) = tag === (st === :I ? :OnlyA : :OnlyI) ? nothing :
                            Symbol("K_", st === :I && tag !== :EqualAI ? "I_" : "A_", lig, "reg")

"""The package's MWC law: its exported per-state polynomials with every dependent replaced by
its monomial in the fitted names and `Keq`, and its regulator-site factors."""
function mwc_today(m::Pop.PMech)
    deps = Dict(m.dependent_monomials)
    f(s) = get(deps, s, nothing)
    sites, mults = site_terms(m, pkg_regname)
    law = mwc_law(m, m.state_polys_A, m.state_polys_I, sites, mults, m.fitted_params, f)
    raw = mwc_law(m, m.state_polys_A, m.state_polys_I, sites, mults, m.fitted_params, s -> nothing)
    alllaw, used = all_names_law(raw, m.dependent_monomials)
    Built(check_names(law, m.id), check_names(alllaw, m.id), used, :mwc, m.cat_multiplicity,
          mults, nothing)
end

"""The redesign's MWC law (P5b's single-graph model): its per-state laws with every name
written through the joint solve in the fitted names and `Keq`."""
function mwc_redesign(m::Pop.PMech)
    md = Allo.derive(m)
    nm = md.names
    fitted = Set(nm.fitted)
    function f(s)
        (s in fitted || s === :Keq || s in m.metabolites) && return nothing
        haskey(nm.folds, s) && return [nm.folds[s] => 1 // 1]
        q, ts = nm.value[s]
        q == 0 ? ts : [ts; :Keq => q]
    end
    st(P) = [symterms(p, md.T) for p in P]
    sites, mults = site_terms(m, Allo.regname)
    law = mwc_law(m, st(md.A), st(md.I), sites, mults, nm.fitted, f)
    raw = mwc_law(m, st(md.A), st(md.I), sites, mults, nm.fitted, s -> nothing)
    defs = [d => (q == 0 ? ts : [ts; :Keq => q]) for (d, (q, ts)) in nm.value
            if !(d in fitted) && ts != [d => 1]]
    alllaw, used = all_names_law(raw, [[d => [t => 1 // 1] for (d, t) in nm.folds]; defs])
    Built(check_names(law, m.id), check_names(alllaw, m.id), used, :mwc, m.cat_multiplicity,
          mults, md)
end

# ─── Numbers ───────────────────────────────────────────────────────────

const PREC = 512

"""Value of a polynomial at `vals` (BigFloat)."""
evalpoly(p, vals) = sum((BigFloat(c) * prod((vals[s]^e for (s, e) in m); init = big(1))
                         for (m, c) in p); init = big(0))

"""`v / E_total` of a built law at `vals` (every parameter, `Keq`, `L` and concentration)."""
function value(b::Built, vals)
    setprecision(BigFloat, PREC) do
        e = [evalpoly(p, vals) for p in b.law.comps]
        b.kind === :plain && return e[1] / e[2]
        nA, qA, nI, qI = e[1] / e[2], e[3] / e[2], e[4] / e[5], e[6] / e[5]
        WA = prod((e[5 + 2k]^b.mults[k] for k in eachindex(b.mults)); init = big(1))
        WI = prod((e[6 + 2k]^b.mults[k] for k in eachindex(b.mults)); init = big(1))
        n, L = b.n, vals[:L]
        (nA * qA^(n - 1) * WA + L * nI * qI^(n - 1) * WI) / (qA^n * WA + L * qI^n * WI)
    end
end

"""`v / E_total` of a redesign law by an independent route: the mass-action oracle at the bare
constants of the joint solve (plain), or P5b's own MWC evaluation (MWC)."""
function reference_value(b::Built, m::Pop.PMech, vals, concs)
    if b.kind === :plain
        sol = b.extra
        idx = Dict(Symbol(sol.cols[i].name) => i for i in sol.indep)
        x = setprecision(BigFloat, PREC) do
            Thermo.all_values(sol, Dict(idx[p] => vals[p] for p in keys(idx)), vals[:Keq])
        end
        return Oracle.steady_state_v(m, Thermo.stored_constants(sol, x), concs)
    end
    setprecision(BigFloat, PREC) do
        fit = Dict{Symbol, BigFloat}(s => vals[s] for s in b.extra.names.fitted)
        fit[:Keq] = vals[:Keq]
        Allo.rate(b.extra, Allo.all_values(b.extra, fit), concs)
    end
end

"""`v / E_total` of the package's law by the oracle (plain): fitted values, dependents from
their monomials, Pass-2 absorbed names at their targets' values."""
function oracle_today(m::Pop.PMech, vals, concs)
    allv = Dict{Symbol, BigFloat}(vals)
    setprecision(BigFloat, PREC) do
        for (d, mono) in m.dependent_monomials
            allv[d] = prod((BigFloat(vals[s])^BigFloat(e) for (s, e) in mono); init = big(1))
        end
    end
    Oracle.steady_state_v(m, Oracle.constants_from_names(m, allv; rename = Dict(m.pass2_rename)),
                          concs)
end

end # module
