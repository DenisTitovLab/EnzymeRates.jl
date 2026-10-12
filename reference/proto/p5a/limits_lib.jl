# ABOUTME: The :OnlyA limit study (P5a): paths along consistent constants, the oracle sweep in ε,
# ABOUTME: the symbolic leading face, deletion, and the dependence on the choice of independents.
include(joinpath(@__DIR__, "pkgvals.jl"))
using Printf, Random, Serialization

const ORDERS = Any[:pkg, :rev, :afirst, (:perm, 1), (:perm, 2), (:perm, 3)]
const EPS = setprecision(() -> [BigFloat(10)^(-k) for k in (4, 6, 8, 10, 12, 14, 20, 30)],
                         BigFloat, MWC.PREC)
logu(rng) = MWC.big1024(10.0^(4 * rand(rng) - 2))
relerr(a, b) = (d = max(abs(a), abs(b)); d == 0 ? 0.0 : Float64(abs(a - b) / d))

"""A rejected uni-uni twin carries every group :NonequalAI; restore the requested tags from
its `source` ("uni:<kinds>:tS,tC,tP")."""
function retag(m)
    startswith(m.source, "uni:") || return m
    tS, tC, tP = Symbol.(split(split(m.source, ":")[3], ","))
    tags = map(m.groups) do g
        s = m.steps[g[1]]
        any(x -> x.name === :S, s.consumed) ? tS : any(x -> x.name === :P, s.consumed) ? tP : tC
    end
    d = Dict(f => getfield(m, f) for f in fieldnames(Pop.PMech))
    d[:group_tags] = tags
    Pop.PMech(; d...)
end

"""Random regulator constants per state (an :EqualAI ligand shares its A value)."""
function draw_reg(m, rng)
    A = Dict{Tuple{Int, Symbol}, Any}(); I = Dict{Tuple{Int, Symbol}, Any}()
    for (si, site) in enumerate(m.sites), (lig, t) in zip(site.ligands, site.states)
        A[(si, lig)] = t === :OnlyI ? nothing : logu(rng)
        I[(si, lig)] = t === :OnlyA ? nothing : t === :EqualAI ? A[(si, lig)] : logu(rng)
    end
    A, I
end
regnames(m, st) = Dict((si, lig) => (st === :A ? t === :OnlyI : t === :OnlyA) ? nothing :
                       Symbol("Kreg", st === :I && t === :EqualAI ? :A : st, "_", si, "_", lig)
                       for (si, site) in enumerate(m.sites) for (lig, t) in zip(site.ligands, site.states))

"""Fast-pair scale that keeps every RE step far faster than every SS step."""
function lambda_for(J, x)
    maxss, minre = big(1), big(1)
    for (i, v) in x
        c = J.cols[i]
        if J.cs[c.st][c.j].re
            minre = min(minre, v, 1 / v)
        else
            maxss = max(maxss, v)
        end
    end
    setprecision(() -> big(10)^40 * maxss * 10^4 / minre, BigFloat, MWC.PREC)
end

"""Oracle rate along the path at `eps` (all inactive steps present)."""
function v_path(J, hat, a, eps, keq, L, rA, rI, concs)
    m = J.m
    x = MWC.path_values(J, hat, a, eps, keq)
    MWC.v_factorized(m, MWC.stored(J, x, :A), MWC.stored(J, x, :I), J.re[:I],
                     trues(length(m.steps)), 1:length(m.species), L, rA, rI, concs;
                     n = m.cat_multiplicity, Λ = lambda_for(J, x))
end

"""Oracle rate with today's semantics: :OnlyA steps deleted from the inactive state and the
forms free enzyme can no longer reach pruned (constants at ε = 1)."""
function v_deletion(J, hat, keq, L, rA, rI, concs)
    m = J.m
    x = MWC.all_values(J, hat, keq)
    active = [t !== :OnlyA for t in J.tags]
    MWC.v_factorized(m, MWC.stored(J, x, :A), MWC.stored(J, x, :I), [s.re for s in m.steps],
                     active, MWC.reachable(m, active), L, rA, rI, concs;
                     n = m.cat_multiplicity, Λ = lambda_for(J, x))
end

"""The combined law of `J` at a hat point, built once: the expressions, the inactive
state's polynomials, the values by table variable, and the table index of every
independent column."""
function sym_setup(J, hat, keq, L, rA, rI, concs)
    m = J.m
    tab = MWC.Table(Poly.Vars())
    lawA = MWC.state_law(J, :A, tab); lawI = MWC.state_law(J, :I, tab)
    nA, nI = regnames(m, :A), regnames(m, :I)
    WA = MWC.W_poly(m, nA, tab); WI = MWC.W_poly(m, nI, tab)
    Lv = Poly.var!(tab.T, :L)
    num, den = MWC.combine(lawA, lawI, WA, WI, Poly.term(Poly.P{Int}, Poly.mono(Lv => 1)),
                           m.cat_multiplicity)
    T = tab.T
    vals = Vector{BigFloat}(undef, length(T.names))
    vidx = Dict{Int, Int}()
    setprecision(BigFloat, MWC.PREC) do
        for i in J.indep
            k = get(T.index, Symbol(J.cols[i].name), 0)
            k == 0 && continue
            vidx[i] = k; vals[k] = hat[i]
        end
        vals[T.index[:Keq]] = keq; vals[Lv] = L
        for x in m.metabolites
            haskey(T.index, x) && (vals[T.index[x]] = concs[x])
        end
        for (nm, rv) in ((nA, rA), (nI, rI)), (key, s) in nm
            s === nothing || (vals[T.index[s]] = rv[key])
        end
    end
    (num = num, den = den, lawI = lawI, vals = vals, vidx = vidx)
end

"""The leading face of a built law along exponents `a`: the limit value at the hat point,
the orders of the inactive N, Q, D, and whether a sum cancelled on its face."""
function sym_eval(S, a)
    dir = Dict{Int, MWC.Q}(k => a[i] for (i, k) in S.vidx)
    res = MWC.lead_value(S.num, S.den, [dir], S.vals)
    ordI = [MWC.initial(p, [dir])[1] for p in S.lawI]
    (v = res.v, on = res.on, od = res.od, cancel = res.cancel, ordI = ordI)
end

sym_limit(J, a, hat, keq, L, rA, rI, concs) = sym_eval(sym_setup(J, hat, keq, L, rA, rI, concs), a)

"""Relative difference against a scale: |a − b| / max(|a|, |b|, |scale|)."""
reldiff(a, b, scale) = (d = max(abs(a), abs(b), abs(scale)); d == 0 ? 0.0 : Float64(abs(a - b) / d))

"""Every quantity of the study for one mechanism, convention and draw.

Two path semantics are compared. The RREF path scales the joint solve's independents by
the convention and lets the dependents follow their rows (so it depends on which constants
are independent). The cone path imposes the convention on every moved column, holds every
other constant fixed, and enumerates the feasible directions (`MWC.cone_directions`)."""
function study(m, conv; seed = 1, orders = ORDERS, sweep = true, ncone = 12)
    rng = Random.Xoshiro(seed)
    J = MWC.joint_solve(m; order = :pkg)
    aint = MWC.intended(J, conv)
    a = MWC.realized(J, aint)
    partner, wrong = MWC.path_report(J, aint, a)
    hat = MWC.random_hat(J, rng)
    keq, L = logu(rng), logu(rng)
    rA, rI = draw_reg(m, rng)
    concs = Dict(x => logu(rng) for x in m.metabolites)
    ceq = Oracle.concs_at_ratio(m, rng, keq)
    c10 = Oracle.concs_at_ratio(m, rng, keq / 10)
    out = Dict{Symbol, Any}(:id => m.id, :conv => conv, :admitted => m.ok,
                            :partner => [J.cols[i].name for i in partner],
                            :wrong => [J.cols[i].name for i in wrong],
                            :n => m.cat_multiplicity)
    v1 = v_path(J, hat, a, big(1), keq, L, rA, rI, concs)          # ε = 1: the scale
    out[:v1] = v1
    if sweep
        vs = [v_path(J, hat, a, e, keq, L, rA, rI, concs) for e in EPS]
        cons = maximum(abs(v_path(J, hat, a, e, keq, L, rA, rI, ceq)) /
                       max(abs(v_path(J, hat, a, e, keq, L, rA, rI, c10)), big(1e-300))
                       for e in EPS)
        out[:v_eps] = Float64.(vs)
        out[:conv_err] = [reldiff(v, vs[end], v1) for v in vs]
        out[:consistency] = Float64(cons)
        out[:v_lim_oracle] = vs[end]
    end
    S = sym_setup(J, hat, keq, L, rA, rI, concs)
    s = sym_eval(S, a)
    out[:v_lim_sym] = s.v
    out[:ord_num], out[:ord_den], out[:cancel], out[:ordI] = s.on, s.od, s.cancel, s.ordI
    sweep && (out[:sym_vs_oracle] = reldiff(s.v, out[:v_lim_oracle], v1))
    out[:v_del] = v_deletion(J, hat, keq, L, rA, rI, concs)
    out[:lim_vs_del] = reldiff(s.v, out[:v_del], 0)
    # the same bare point under other choices of independents (RREF-path semantics)
    θ = MWC.all_values(J, hat, keq)
    key(Jx, i) = (Jx.cols[i].st, Jx.cols[i].j, Jx.cols[i].r)
    θk = Dict(key(J, i) => v for (i, v) in θ)
    per = Dict{Any, Any}()
    for o in orders
        o === :pkg && continue
        Jo = MWC.joint_solve(m; order = o)
        ho = Dict(i => θk[key(Jo, i)] for i in Jo.indep)
        ao = MWC.realized(Jo, MWC.intended(Jo, conv))
        po, wo = MWC.path_report(Jo, MWC.intended(Jo, conv), ao)
        so = sym_limit(Jo, ao, ho, keq, L, rA, rI, concs)
        per[o] = (v = so.v, partner = [Jo.cols[i].name for i in po],
                  wrong = [Jo.cols[i].name for i in wo])
    end
    out[:per_order] = per
    out[:order_spread] = maximum((reldiff(p.v, s.v, 0) for p in values(per)); init = 0.0)
    # cone semantics: every feasible direction that moves only the :OnlyA columns
    dirs, d = MWC.cone_directions(J, conv)
    out[:cone_dim], out[:cone_n] = d, length(dirs)
    out[:rref_in_cone] = isempty(partner) && isempty(wrong)
    if !isempty(dirs)
        pick = unique([1; round.(Int, range(1, length(dirs); length = min(ncone, length(dirs))))])
        res = map(pick) do k
            ak = dirs[k]
            sk = sym_eval(S, ak)
            (a = [Int(ak[i]) for i in sort!(collect(keys(MWC.cone_signs(J, conv))))],
             v = sk.v, ordI = sk.ordI,
             seq = sk.ordI[2] !== nothing && sk.ordI[3] !== nothing && sk.ordI[2] < sk.ordI[3])
        end
        out[:cone] = res
        out[:cone_v] = res[1].v
        out[:cone_spread] = maximum(reldiff(r.v, res[1].v, 0) for r in res)
        out[:cone_seq] = any(r -> r.seq, res)
        out[:cone_vs_del] = maximum(reldiff(r.v, out[:v_del], 0) for r in res)
        # the oracle along the first and last picked directions
        epsk(ak) = setprecision(() -> BigFloat(10)^(-max(6, 30 ÷ max(1, Int(maximum(abs, ak))))),
                                BigFloat, MWC.PREC)
        out[:cone_oracle_err] = maximum(reldiff(
            v_path(J, hat, dirs[k], epsk(dirs[k]), keq, L, rA, rI, concs), r.v, v1)
            for (k, r) in zip(pick, res) if k in (pick[1], pick[end]))
    end
    out
end
