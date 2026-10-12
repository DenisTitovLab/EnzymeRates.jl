# ABOUTME: Validates the single-graph allosteric derivation (Allo.jl) against the package's exported
# ABOUTME: outputs, a two-conformation mass-action oracle at n = 1, and finite-difference ranks.
using Random, Printf, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "Allo.jl"))
include(joinpath(@__DIR__, "..", "common", "Oracle.jl"))

const POPS = isempty(ARGS) ? ["specs", "handwritten", "allo_uniuni", "allo_bibi", "expand_bibi"] :
             ARGS
const RANK_ALL = get(ENV, "RANK_ALL", "1") == "1"
big5(x) = setprecision(() -> BigFloat(x), BigFloat, 512)
relerr(a, b) = abs(a - b) / max(abs(a), abs(b), big(10)^-300)

# ─── The package's law, rebuilt from its exported pieces ────────────────

"""The package's values of every name at fitted values `fit` (with `:Keq`): its dependents from
their monomials."""
function pkg_values(m, fit::AbstractDict)
    v = Dict{Symbol, BigFloat}(k => big5(x) for (k, x) in fit)
    for (d, mono) in m.dependent_monomials
        v[d] = prod((v[s]^BigFloat(e) for (s, e) in mono); init = big(1))
    end
    v
end

"""The package's `v / E_t` from its exported per-state polynomials (normalized per state; the
package's three renderings are equal to this) and regulator factors."""
function pkg_rate(m, vals, concs; n = m.cat_multiplicity, regs = true)
    all = merge(Dict{Symbol, BigFloat}(k => big5(x) for (k, x) in concs), vals)
    ev(p) = Oracle.eval_poly(p, all)
    nA, qA, dA = ev.(m.state_polys_A)
    nI, qI, dI = ev.(m.state_polys_I)
    WA = regs ? Allo.reg_factor(m, :A, vals, concs) : big(1)
    WI = regs ? Allo.reg_factor(m, :I, vals, concs) : big(1)
    a, b, c, d, L = nA / dA, qA / dA, nI / dI, qI / dI, vals[:L]
    (a * b^(n - 1) * WA + L * c * d^(n - 1) * WI) / (b^n * WA + L * d^n * WI)
end

"""The bare stored constants `(state, step, slot) => value` the package's values imply: active
names from `step_names_A`, inactive ones from the names its pruned graph renders, each through
that state's Pass-2 rename."""
function pkg_bare(m, vals)
    rA, rI = Dict(m.pass2_rename_A), Dict(m.pass2_rename_I)
    x = Dict{Allo.Col, BigFloat}()
    for (j, s) in enumerate(m.steps)
        for r in 1:Allo.nslots(s)
            n = m.step_names_A[j][r]; x[(:A, j, r)] = vals[get(rA, n, n)]
        end
        m.inactive_kept[j] || continue
        val(n) = vals[get(rI, n, n)]
        ns = m.step_names_I_pkg[j]
        if !m.inactive_reversed[j]
            for r in 1:Allo.nslots(s); x[(:I, j, r)] = val(ns[r]); end
        elseif s.re
            x[(:I, j, 1)] = Pop.is_binding(s) ? val(ns[1]) : 1 / val(ns[1])
        else
            x[(:I, j, 1)] = val(ns[2]); x[(:I, j, 2)] = val(ns[1])
        end
    end
    x
end

"""The package name of each of my columns (for reading package values off my bare point)."""
function pkg_colname(m, (st, j, r)::Allo.Col)
    if st === :A
        n = m.step_names_A[j][r]; get(Dict(m.pass2_rename_A), n, n)
    else
        ns = m.step_names_I_pkg[j]
        n = m.inactive_reversed[j] ? (m.steps[j].re ? ns[1] : ns[3 - r]) : ns[r]
        get(Dict(m.pass2_rename_I), n, n)
    end
end

"""The package's normalization case from its exported free-enzyme weights."""
function pkg_norm_case(m)
    dA, dI = m.state_polys_A[3], m.state_polys_I[3]
    dA == dI && return :same
    mf(d) = length(d) == 1 && !any(s in m.metabolites for (s, _) in d[1].first)
    mf(dA) && mf(dI) ? :divide : :cross
end

"""Whether a weight carries a concentration (`:met`), else a monomial (`:mono`) or a sum
(`:poly`)."""
dkind(d, mets) = any(s in mets for (mo, _) in d for (s, _) in mo) ? :met :
                 length(d) == 1 ? :mono : :poly

"""The package's fold hazards: a pair folded in opposite directions in the two states (shared
symbols, or twin names), and a symbol folded in one state that the other state still uses."""
function pkg_fold_issues(m)
    rA, rI = Dict(m.pass2_rename_A), Dict(m.pass2_rename_I)
    syms(ps) = Set(s for p in ps for (mo, _) in p for (s, _) in mo)
    sA, sI = syms(m.state_polys_A), syms(m.state_polys_I)
    twin = Dict(m.step_names_A[j][r] => m.step_names_I_pkg[j][r]
                for j in eachindex(m.steps) if m.inactive_kept[j]
                for r in 1:Allo.nslots(m.steps[j]))
    opp = any(get(rI, y, nothing) == x for (x, y) in rA) ||
          any(haskey(twin, x) && haskey(twin, y) && get(rI, twin[y], nothing) == twin[x]
              for (x, y) in rA)
    one = any(x in sI for x in keys(rA)) || any(x in sA for x in keys(rI))
    opp, one
end

"""The constants the mechanism declares: every group's active names, the inactive names of
its non-`:OnlyA` groups, the regulator constants and `L`."""
function declared(m)
    g = Pop.group_of_step(m)
    d = Set{Symbol}([:L])
    for j in eachindex(m.steps)
        union!(d, m.step_names_A[j])
        m.group_tags[g[j]] === :OnlyA || union!(d, m.step_names_I[j])
    end
    for site in m.sites, (l, t) in zip(site.ligands, site.states), st in (:A, :I)
        n = Allo.regname(st, l, t); n === nothing || push!(d, n)
    end
    d
end

include(joinpath(@__DIR__, "validate_lib.jl"))

# ─── Finite-difference identifiable rank ─────────────────────────────────

"""Rank of `J` by Gaussian elimination with complete pivoting: pivots above `tol` times the
first; also the smallest kept and the largest dropped pivot, relative to the first."""
function rank_gap(J; tol = big(10)^-35)
    A = copy(J); m, n = size(A); piv = BigFloat[]
    for k in 1:min(m, n)
        _, idx = findmax(abs.(A[k:m, k:n]))
        i, j = idx[1] + k - 1, idx[2] + k - 1
        A[[k, i], :] = A[[i, k], :]; A[:, [k, j]] = A[:, [j, k]]
        push!(piv, abs(A[k, k]))
        iszero(A[k, k]) && break
        for r in (k + 1):m
            f = A[r, k] / A[k, k]
            A[r, k:n] .-= f .* A[k, k:n]
        end
    end
    rk = count(p -> p > tol * piv[1], piv)
    rk, Float64(piv[rk] / piv[1]), rk < length(piv) ? Float64(piv[rk + 1] / piv[1]) : 0.0
end

"""Identifiable rank of `f(θ, concs)` in the log of the parameters `θ0` over `npts` random
concentration points: central differences in 256-bit BigFloat."""
function fd_rank(f, θ0::Vector{BigFloat}, mets, rng)
    setprecision(BigFloat, 256) do
        p = length(θ0); npts = 2p + 6; h = big(10)^-25
        J = zeros(BigFloat, npts, p)
        for i in 1:npts
            c = Dict(x => BigFloat(10)^(4 * big(rand(rng)) - 2) for x in mets)
            for k in 1:p
                θp = copy(θ0); θp[k] *= exp(h)
                θm = copy(θ0); θm[k] *= exp(-h)
                J[i, k] = (f(θp, c) - f(θm, c)) / (2h)
            end
        end
        rank_gap(J)
    end
end

logu(rng) = big5(10)^(4 * big5(rand(rng)) - 2)

# ─── Text checks ─────────────────────────────────────────────────────────

const ARTIFACT = r"(?<![\w.])1 \^ \d|\^ 1(?![\d.])|(?<![\w.^])1 \* |\* 1(?![\d.\w])|L \* 0"

"""A polynomial with every name replaced by its tag-free role name (`rolekey`)."""
function stripped(p, T)
    S = Allo.Vars()
    Allo.subst(p, v -> Allo.mono(Allo.var!(S, Symbol(Allo.rolekey(T.names[v])[1])) => 1)), S
end

"""The package's printed term order of a polynomial (its `_poly_to_expr` sort with the active
state's names as parameters), each term written with tag-free names."""
function pkg_term_order(p, params)
    key(x) = (sum((e for (s, e) in x.first if s ∉ params); init = 0), x.second < 0,
              Tuple((string(s), e) for (s, e) in x.first))
    [join([string(Allo.rolekey(s)[1], "^", e) for (s, e) in
           sort(x.first; by = sp -> (sp.first in params ? 0 : 1, string(sp.first)))], "*")
     for x in sort(p; by = key)]
end

# ─── One mechanism ───────────────────────────────────────────────────────

function check(m::Pop.PMech, rng)
    r = Dict{Symbol, Any}(:id => m.id, :source => m.source, :n => m.cat_multiplicity,
                          :tags => m.group_tags, :pkg_count => length(m.fitted_params),
                          :pkg_fitted => m.fitted_params)
    t0 = time()
    md = try
        Allo.derive(m)
    catch e
        r[:error] = sprint(showerror, e); return r
    end
    r[:derive_s] = time() - t0
    nm = md.names
    r[:keepI_match] = md.keepI == m.inactive_kept
    r[:count] = length(nm.fitted); r[:fitted] = nm.fitted
    decl = declared(m)
    r[:undeclared_mine] = setdiff(nm.fitted, decl)
    r[:undeclared_pkg] = setdiff(m.fitted_params, decl)
    r[:folds] = length(nm.folds); r[:collapses] = count(==(:collapse), values(nm.tie_kind))
    r[:ndeps] = length(nm.deps)
    r[:re_subst] = count(c -> m.steps[c[2]].re && haskey(nm.deps, Allo.colname(m, c)), md.cols)
    r[:nonint] = any(!isinteger(e) for (q, ts) in md.expr for (_, e) in ts) ||
                 any(!isinteger(q) for (q, _) in md.expr)
    # (b) mine: a folded or dependent name never reaches the law or the fitted list
    lawsyms = Set(md.T.names[v] for p in (md.A..., md.I...) for mo in keys(p.t) for (v, _) in mo)
    r[:dep_in_law] = count(n -> n in lawsyms, keys(nm.folds)) +
                     count(n -> n in nm.fitted, [collect(keys(nm.folds)); collect(keys(nm.deps))])
    twin = Dict(Allo.colname(m, (:A, j, k)) => Allo.colname(m, (:I, j, k))
                for j in eachindex(m.steps) if md.keepI[j] for k in 1:Allo.nslots(m.steps[j]))
    r[:opp_mine] = count(((x, y),) -> haskey(twin, x) && haskey(twin, y) && twin[x] != x &&
                         get(nm.folds, twin[y], nothing) == twin[x], collect(nm.folds))
    r[:pkg_opp], r[:pkg_one] = pkg_fold_issues(m)
    r[:norm_pkg] = pkg_norm_case(m); r[:norm_mine] = Allo.norm_case(md)
    mets = Set(m.metabolites)
    r[:dkind_pkg] = (dkind(m.state_polys_A[3], mets), dkind(m.state_polys_I[3], mets))
    conc_vars = Set(md.T.index[x] for x in m.metabolites if haskey(md.T.index, x))
    mk(d) = any(v in conc_vars for mo in keys(d.t) for (v, _) in mo) ? :met :
            Allo.nterms(d) == 1 ? :mono : :poly
    r[:dkind_mine] = (mk(md.A[3]), mk(md.I[3]))
    r[:terms] = (Allo.nterms(md.A[1]), Allo.nterms(md.A[2]), Allo.nterms(md.I[1]), Allo.nterms(md.I[2]))

    setprecision(BigFloat, 512) do
        # (d) the package's samples: its bare point against my constraints, then the rates
        cons, vs, vb, pk = 0.0, 0.0, 0.0, 0.0
        for s in m.samples
            pv = pkg_values(m, Dict(s.params))
            x = pkg_bare(m, pv)
            for (i, (q, ts)) in enumerate(md.expr)
                i in md.indep && continue
                pred = pv[:Keq]^BigFloat(q) * prod((x[md.cols[f]]^BigFloat(e) for (f, e) in ts); init = big(1))
                cons = max(cons, Float64(relerr(pred, x[md.cols[i]])))
            end
            fit = Dict{Symbol, BigFloat}(:Keq => pv[:Keq], :L => pv[:L])
            for i in md.indep; fit[Allo.colname(m, md.cols[i])] = x[md.cols[i]]; end
            for n in nm.fitted; n === :L || haskey(fit, n) || (fit[n] = pv[n]); end
            c = Dict(s.concs)
            vm = Allo.rate(md, Allo.all_values(md, fit), c)
            vp = pkg_rate(m, pv, c)
            vs = max(vs, Float64(relerr(vm, big5(s.v))))
            vb = max(vb, Float64(relerr(vm, vp)))
            pk = max(pk, Float64(relerr(vp, big5(s.v))))
        end
        r[:cons_err], r[:v_vs_sample], r[:v_vs_pkgbig], r[:pkg_rebuild] = cons, vs, vb, pk

        # (d1) my random point read into the package's names: does its law reproduce mine?
        spread, vr, bare = 0.0, 0.0, 0.0
        for _ in 1:2
            fit = Dict{Symbol, BigFloat}(n => logu(rng) for n in nm.fitted)
            fit[:Keq] = logu(rng)
            vals = Allo.all_values(md, fit)
            got = Dict{Symbol, Vector{BigFloat}}()
            for c in md.cols
                push!(get!(got, pkg_colname(m, c), BigFloat[]), vals[Allo.colname(m, c)])
            end
            pf = Dict{Symbol, BigFloat}(:Keq => fit[:Keq])
            for n in m.fitted_params
                if haskey(got, n)
                    g = got[n]; pf[n] = g[1]
                    spread = max(spread, maximum(Float64(relerr(y, g[1])) for y in g))
                else
                    pf[n] = vals[n]                     # regulators and L share my names
                end
            end
            pv = pkg_values(m, pf)
            xb = pkg_bare(m, pv)
            bare = max(bare, maximum(Float64(relerr(xb[c], vals[Allo.colname(m, c)])) for c in md.cols))
            c = Oracle.concs_at_ratio(m, rng, 10 * logu(rng))
            vr = max(vr, Float64(relerr(Allo.rate(md, vals, c), pkg_rate(m, pv, c))))
        end
        r[:restrict_spread], r[:restrict_v], r[:restrict_bare] = spread, vr, bare

        # detailed balance of my law and the package's, and the n = 1 oracle
        db, dbp, orc = 0.0, 0.0, 0.0
        for _ in 1:2
            fit = Dict{Symbol, BigFloat}(n => logu(rng) for n in nm.fitted); fit[:Keq] = logu(rng)
            vals = Allo.all_values(md, fit)
            c1 = Oracle.concs_at_ratio(m, rng, fit[:Keq])
            c2 = copy(c1); c2[m.products[1]] /= 10
            db = max(db, Float64(abs(Allo.rate(md, vals, c1)) / abs(Allo.rate(md, vals, c2))))
            pf = Dict{Symbol, BigFloat}(n => logu(rng) for n in m.fitted_params); pf[:Keq] = fit[:Keq]
            pv = pkg_values(m, pf)
            dbp = max(dbp, Float64(abs(pkg_rate(m, pv, c1)) / abs(pkg_rate(m, pv, c2))))
            pm, k = combined(md, vals)
            c = Oracle.concs_at_ratio(m, rng, 10 * logu(rng))
            orc = max(orc, Float64(relerr(Oracle.steady_state_v(pm, k, c),
                                          Allo.rate(md, vals, c; n = 1, regs = false))))
        end
        r[:db_mine], r[:db_pkg], r[:oracle_n1] = db, dbp, orc
    end

    # (f) the printed law
    txt = Allo.law_text(md)
    r[:text] = txt
    r[:artifacts_mine] = length(collect(eachmatch(ARTIFACT, txt)))
    r[:artifacts_pkg] = length(collect(eachmatch(ARTIFACT, m.rate_equation_string)))
    isconc(s) = s in m.metabolites
    sym, same_txt = 0, 0
    for k in 1:2
        (sA, SA), (sI, SI) = stripped(md.A[k], md.T), stripped(md.I[k], md.T)
        Allo.Poly.symbolic(sA, SA) == Allo.Poly.symbolic(sI, SI) || continue
        sym += 1
        strip_txt(p) = replace(Allo.poly_text(p, md.T, isconc), r"\b([kK])_[AI]_" => s"\1_")
        same_txt += strip_txt(md.A[k]) == strip_txt(md.I[k])
    end
    r[:sym_pairs], r[:sym_same_text] = sym, same_txt
    params = Set(n for ns in m.step_names_A for n in ns)
    psym, psame = 0, 0
    for k in 1:2
        a, b = m.state_polys_A[k], m.state_polys_I[k]
        key(p) = sort([sort([Allo.rolekey(s)[1] => e for (s, e) in mo]) => c for (mo, c) in p])
        key(a) == key(b) || continue
        psym += 1
        psame += pkg_term_order(a, params) == pkg_term_order(b, params)
    end
    r[:pkg_sym_pairs], r[:pkg_sym_same_order] = psym, psame

    # (c) identifiable ranks where the counts differ (and everywhere when RANK_ALL)
    if RANK_ALL || r[:count] != r[:pkg_count]
        θ = [logu(rng) for _ in nm.fitted]; keq = logu(rng)
        fm(θv, c) = Allo.rate(md, Allo.all_values(md, merge(Dict(nm.fitted .=> θv), Dict(:Keq => keq))), c)
        r[:rank_mine] = fd_rank(fm, θ, m.metabolites, rng)
        θp = [logu(rng) for _ in m.fitted_params]
        fp(θv, c) = pkg_rate(m, pkg_values(m, merge(Dict(m.fitted_params .=> θv), Dict(:Keq => keq))), c)
        r[:rank_pkg] = fd_rank(fp, θp, m.metabolites, rng)
    end
    r
end

# ─── Run ─────────────────────────────────────────────────────────────────

results = Dict{String, Vector{Dict{Symbol, Any}}}()
for name in POPS
    pop = filter(m -> m.allosteric, Pop.load_population(name))
    rng = Xoshiro(20261009)
    t = @elapsed rs = [check(m, rng) for m in pop]
    results[name] = rs
    @printf("%s: %d mechanisms in %.1f s (%d errors)\n", name, length(pop), t,
            count(r -> haskey(r, :error), rs))
    flush(stdout)
    serialize(joinpath(@__DIR__, "data", "results_$(name).jls"), rs)
end
