# ABOUTME: Helpers of the kcat saturation census: a slice law evaluated numerically at one
# ABOUTME: parameter point, kcat classes, growth scans, trap analysis, sup search, gauge check.
include("/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/slice/integrate/lib.jl")
using Random, Printf, Serialization, Statistics, LinearAlgebra

const KDIR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/round4/kcat"
const KDATA = joinpath(KDIR, "data")

# ─── numeric polynomials in the concentrations ───────────────────────────────────────────

"""A polynomial in the concentrations with Float64 coefficients: term `t` is
`coef[t] * Π c[idx[t][j]]^ex[t][j]`."""
struct NPoly
    coef::Vector{Float64}
    idx::Vector{Vector{Int}}
    ex::Vector{Vector{Int}}
end

"""`p` with its constants evaluated at `xv` (indexed by variable; concentration entries are
ignored)."""
function NPoly(p::S.Poly, vars::S.Vars, xv::Vector{Float64})
    coef, idx, ex = Float64[], Vector{Int}[], Vector{Int}[]
    for (mo, cp) in S.split_conc(p, vars)
        all(e -> isinteger(e), last.(mo)) || error("non-integer concentration exponent")
        push!(coef, S.evaluate(cp, xv))
        push!(idx, [v for (v, _) in mo])
        push!(ex, [Int(e) for (_, e) in mo])
    end
    NPoly(coef, idx, ex)
end

function ev(p::NPoly, c::AbstractVector{T}) where {T}
    s = zero(T)
    @inbounds for t in eachindex(p.coef)
        x = T(p.coef[t])
        for j in eachindex(p.idx[t])
            x *= c[p.idx[t][j]]^p.ex[t][j]
        end
        s += x
    end
    s
end

"""A law at one parameter point: `v/E_total` as a function of the concentrations
(`metabolites(m)` order). Plain: `P = [Q, D]`; MWC: `[QA, DA, QI, DI, fA, fI, WA_s..., WI_s...]`."""
struct NLaw
    nc::Int
    subs::Vector{Int}
    prods::Vector{Int}
    mods::Vector{Int}
    keq::Float64
    allo::Bool
    P::Vector{NPoly}
    mult::Vector{Int}
    n::Int
    Lval::Float64
end

"""Indices of the substrates, products and modifiers (names that are no reactant) in
`metabolites(m)`."""
function roles(m)
    names = S.metabolites(m)
    ([i for i in eachindex(names) if names[i] in m.substrates],
     [i for i in eachindex(names) if names[i] in m.products],
     [i for i in eachindex(names) if !S.is_reactant(m, names[i])])
end

"""Variable values of `vars` at the parameter values `x` (by name), 1.0 for concentrations."""
values_of(vars::S.Vars, x) = Float64[vars.conc[i] ? 1.0 : Float64(x[n])
                                     for (i, n) in enumerate(vars.names)]

function NLaw(m, L, x)
    subs, prods, mods = roles(m)
    nc = length(S.metabolites(m))
    vars = L.vars
    xv = values_of(vars, x)
    np(p) = NPoly(p, vars, xv)
    if L isa S.Law
        return NLaw(nc, subs, prods, mods, Float64(x[:Keq]), false, [np(L.Q), np(L.D)],
                    Int[], 1, 0.0)
    end
    QI, DI = S.inactive_terms(L)
    P = [np(L.active.Q), np(L.active.D), np(QI), np(DI), np(L.fA), np(L.fI)]
    append!(P, [np(s[1]) for s in L.sites])
    append!(P, [np(s[2]) for s in L.sites])
    NLaw(nc, subs, prods, mods, Float64(x[:Keq]), true, P, [s[3] for s in L.sites], L.n,
         Float64(x[:L]))
end

"""`v/E_total` of `nl` at the concentrations `c`."""
function rate(nl::NLaw, c::AbstractVector{T}) where {T}
    H = prod(c[i] for i in nl.subs) - prod(c[i] for i in nl.prods) / T(nl.keq)
    nl.allo || return H * ev(nl.P[1], c) / ev(nl.P[2], c)
    QA, DA, QI, DI, fA, fI = (ev(nl.P[i], c) for i in 1:6)
    ns = length(nl.mult)
    WA = prod((ev(nl.P[6 + s], c)^nl.mult[s] for s in 1:ns); init = one(T))
    WI = prod((ev(nl.P[6 + ns + s], c)^nl.mult[s] for s in 1:ns); init = one(T))
    n, L = nl.n, T(nl.Lval)
    H * (fI^n * QA * DA^(n - 1) * WA + L * fA^n * QI * DI^(n - 1) * WI) /
        (fI^n * DA^n * WA + L * fA^n * DI^n * WI)
end

# ─── derivation and kcat classes ─────────────────────────────────────────────────────────

"""`(ok, reason, sol, L)`: whether the slice derives `m`'s law (boundary checks, admission,
the solve and the law, as `derive` runs them before printing and code)."""
function derivable(m)
    try
        S.check_boundary(m)
        rej = S.admit(m)
        rej === nothing || return (false, "admit: " * first(rej, 80), nothing, nothing)
        sol = S.solve(m)
        L = S.rate_law(m, sol)
        (true, "", sol, L)
    catch err
        (false, "error: " * first(sprint(showerror, err), 80), nothing, nothing)
    end
end

"""The kcat class of `m`: `(:finite | :zero | :infinite | :nounique | :other, msg, vars,
terms)`; `:zero` when every candidate's numerator vanishes identically."""
function kcat_class(m, sol)
    vars, terms = try
        S.kcat_candidates(m, sol)
    catch err
        msg = err isa ErrorException ? err.msg : sprint(showerror, err)
        cls = occursin("never saturates", msg) ? :infinite :
              occursin("depends on how", msg) ? :nounique : :other
        return (cls, first(msg, 160), nothing, nothing)
    end
    (all(t -> iszero(t[1]), terms) ? :zero : :finite, "", vars, terms)
end

"""kcat at parameter values `x` from the candidates `(vars, terms)`."""
function kcat_value(vars, terms, x)
    xv = values_of(vars, x)
    maximum(S.evaluate(n, xv) / S.evaluate(d, xv) for (n, d) in terms)
end

"""Random consistent parameter values of `m` (log-uniform in [1e-2, 1e2]; dependents follow
from the independents and Keq): `x` by name (BigFloat) and every column's value `cv`."""
params_at(m, sol, rng) = random_values(m, sol, rng)

# ─── structure ───────────────────────────────────────────────────────────────────────────

"""Whether form `i` holds a product (bound, either role) or a residual."""
holds_product(m, i) = (f = m.forms[i]; f.res != 0 ||
                       any(x -> x.name in m.products, S.mets_of(m, f.bound)))

"""Copies of `X` (substrate or its modifier copy) form `i` holds."""
holds(m, i, X) = count(x -> x.name === X, S.mets_of(m, m.forms[i].bound))

"""
Structural features of `m` restricted to the edges in `mask`: SS and RE edge counts,
ping-pong (a form with a residual), and the support of the non-turnover cycle space (cycles
whose net change is no reaction): `ss_nt` SS edges on it, split by whether both end forms
hold no product (`ss_nt_sub`, alternative substrate-binding orders) or not (`ss_nt_prod`).
"""
function features(m, mask = trues(length(m.edges)))
    fwd = S.orient(m)
    zs, cs = S.cycles(m, fwd)
    E = length(m.edges)
    basis = Vector{Int}[]
    ref = findfirst(!=(0), cs)
    for (z, c) in zip(zs, cs)
        if c == 0
            push!(basis, z)
        elseif ref !== nothing && z !== zs[ref]
            push!(basis, cs[ref] .* z .- c .* zs[ref])
        end
    end
    support = falses(E)
    foreach(z -> support .|= (z .!= 0), basis)
    support .&= mask
    ss_nt = [k for k in 1:E if support[k] && !m.edges[k].re]
    sub = [k for k in ss_nt if !holds_product(m, m.edges[k].a) && !holds_product(m, m.edges[k].b)]
    (nss = count(k -> mask[k] && !m.edges[k].re, 1:E), nre = count(k -> mask[k] && m.edges[k].re, 1:E),
     pingpong = any(k -> mask[k] && (m.forms[m.edges[k].a].res != 0 || m.forms[m.edges[k].b].res != 0), 1:E),
     nt = count(support), ss_nt = length(ss_nt), ss_nt_sub = length(sub),
     ss_nt_prod = length(ss_nt) - length(sub))
end

"""The steps of `m` as the DSL writes them, physical-forward, with group and core flags."""
function step_lines(m; core = S.catalytic_core(m))
    fwd = S.orient(m)
    [string(S.step_text(m, S.oriented(m, k, fwd)), "   [g", m.group[k],
            isempty(m.tags) ? "" : " " * String(m.tags[m.group[k]]),
            core[k] ? "" : ", dead-end branch", "]") for k in eachindex(m.edges)]
end

# ─── symbolic limits at zero products (plain laws) ───────────────────────────────────────

"""`p` with every term holding one of the concentration variables `zero_vars` dropped."""
drop_vars(p::S.Poly{T}, zero_vars) where {T} =
    S.Poly{T}(filter(t -> !any(v -> v.first in zero_vars, t.first), p.terms))

"""Highest exponent of variable `v` in `p` (-1 for the zero polynomial)."""
deg(p::S.Poly, v) = isempty(p.terms) ? -1 : Int(maximum(mo -> S.mexp(mo, v), keys(p.terms)))

"""For a plain law `L` at zero products and modifiers: whether `Q` keeps a term (then
`v > 0` at every positive substrate concentration, all coefficients being positive), the
signs check, and per substrate `deg_X(ΠS·Q) - deg_X(D)` (negative: `v → 0` as X grows)."""
function zero_product_limits(m, L::S.Law)
    subs, prods, mods = roles(m)
    Q0, D0 = drop_vars(L.Q, [prods; mods]), drop_vars(L.D, [prods; mods])
    pos = all(>(0), values(L.Q.terms)) && all(>(0), values(L.D.terms))
    (qkeeps = !isempty(Q0.terms), dkeeps = !isempty(D0.terms), positive = pos,
     gap = [deg(Q0, i) + 1 - deg(D0, i) for i in subs])
end

# ─── numeric scans ───────────────────────────────────────────────────────────────────────

"""log-uniform value in [10^lo, 10^hi]."""
logu10(rng, lo, hi) = 10.0^(lo + (hi - lo) * rand(rng))

"""Growth scans at zero products and modifiers: from a random base point (log-uniform in
[1e-3, 1e3]), grow each substrate alone, then all together, by 10^0..10^8. Returns the base
point and per scan the nine rates."""
function growth_scans(nl::NLaw, rng)
    base = zeros(nl.nc)
    for i in nl.subs
        base[i] = logu10(rng, -3, 3)
    end
    scan(ids) = [rate(nl, (c = copy(base); c[ids] .*= 10.0^k; c)) for k in 0:8]
    base, [[scan([i]) for i in nl.subs]; [scan(nl.subs)]]
end

"""Minimizes `f` from `x0` by Nelder–Mead (box `[lo, hi]` by clamping), `iters` steps."""
function nelder_mead(f, x0, step, lo, hi; iters = 300)
    d = length(x0)
    cl(x) = clamp.(x, lo, hi)
    pts = [cl(x0)]
    for j in 1:d
        x = copy(x0); x[j] += step; push!(pts, cl(x))
    end
    fs = f.(pts)
    for _ in 1:iters
        o = sortperm(fs); pts, fs = pts[o], fs[o]
        cen = sum(pts[1:d]) / d
        xr = cl(cen .+ (cen .- pts[end])); fr = f(xr)
        if fr < fs[1]
            xe = cl(cen .+ 2 .* (cen .- pts[end])); fe = f(xe)
            fe < fr ? (pts[end] = xe; fs[end] = fe) : (pts[end] = xr; fs[end] = fr)
        elseif fr < fs[end - 1]
            pts[end] = xr; fs[end] = fr
        else
            xc = cl(cen .+ 0.5 .* (pts[end] .- cen)); fc = f(xc)
            if fc < fs[end]
                pts[end] = xc; fs[end] = fc
            else
                for j in 2:(d + 1)
                    pts[j] = cl(pts[1] .+ 0.5 .* (pts[j] .- pts[1])); fs[j] = f(pts[j])
                end
            end
        end
        maximum(abs.(fs .- fs[1])) <= 1e-13 * abs(fs[1]) && break
    end
    i = argmin(fs)
    pts[i], fs[i]
end

"""The largest `v/E_total` found at zero products over substrates and (with `mods`)
modifiers in [1e-4, 1e4] (log space; without `mods` the modifiers stay at zero): `nrand`
random points, then Nelder–Mead from the best three. Returns `(sup, argmax concentrations)`."""
function sup_search(nl::NLaw, rng; nrand = 0, mods = true, lo = -4.0, hi = 4.0)
    vary = mods ? [nl.subs; nl.mods] : nl.subs
    d = length(vary)
    nrand = nrand > 0 ? nrand : 150 * d
    c = zeros(nl.nc)
    function f(y)
        for (j, i) in enumerate(vary)
            c[i] = 10.0^y[j]
        end
        v = rate(nl, c)
        isfinite(v) ? -v : Inf
    end
    ys = [lo .+ (hi - lo) .* rand(rng, d) for _ in 1:nrand]
    fs = f.(ys)
    best = sortperm(fs)[1:min(3, nrand)]
    out = [nelder_mead(f, ys[b], 0.5, lo, hi) for b in best]
    k = argmin(last.(out))
    y = first(out[k])
    -last(out[k]), Dict(zip(vary, 10.0 .^ y))
end

"""Rates at 50 random concentration vectors (every concentration log-uniform in
[1e-3, 1e3]; the first `nzero` with the products at zero)."""
function gauge_rates(nl::NLaw, rng; npts = 50, nzero = 10)
    map(1:npts) do i
        c = [logu10(rng, -3, 3) for _ in 1:nl.nc]
        i <= nzero && (c[nl.prods] .= 0.0)
        rate(nl, c)
    end
end
