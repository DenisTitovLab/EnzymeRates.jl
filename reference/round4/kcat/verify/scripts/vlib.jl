# ABOUTME: Shared helpers of the kcat-claims verifier: record loading, exact law analysis at
# ABOUTME: zero products, fast numeric laws, a mass-action oracle and hand-built mechanisms.

const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
const VD = joinpath(SP, "round4", "kcat", "verify")
include(joinpath(SP, "proto", "common", "Population.jl"))
include(joinpath(SP, "slice", "src", "Slice.jl"))
const S = Slice
using Serialization, Random, Printf, LinearAlgebra, Statistics
setprecision(BigFloat, 512)

const CDATA = joinpath(SP, "slice", "core", "data")
const IDATA = joinpath(SP, "slice", "integrate", "data")
const KDATA = joinpath(SP, "round4", "kcat", "data")       # the probe's data (ids only)
const CORE = [String(split(l, '\t')[1])
              for l in readlines(joinpath(CDATA, "MANIFEST.tsv"))[2:end]]

# ─── records ─────────────────────────────────────────────────────────────────────────────

"""Every record of a converted core set as `(set, index, m)`."""
core_records(set) = [(set = set, index = c.index, m = c.mech)
                     for c in deserialize(joinpath(CDATA, set * ".jls"))]

"""Records of a ter-ter Pop export (`pop_<set>.jls`), all or the indices `idx`."""
function pop_records(set, idx = nothing)
    ps = deserialize(joinpath(IDATA, "pop_$set.jls"))
    out = NamedTuple[]
    for i in (idx === nothing ? eachindex(ps) : idx)
        m = try
            S.Converted(ps[i], i).mech
        catch
            nothing
        end
        m === nothing || push!(out, (set = set, index = i, m = m))
    end
    ps = nothing
    GC.gc()
    out
end

"""The steps of `m` in DSL text, physical-forward, with group labels."""
function steps_text(m)
    fwd = S.orient(m)
    [S.step_text(m, S.oriented(m, k, fwd)) * "  [g$(m.group[k])" *
     (isempty(m.tags) ? "" : " " * String(m.tags[m.group[k]])) * "]" for k in eachindex(m.edges)]
end

# ─── derivation under test ───────────────────────────────────────────────────────────────

"""`(ok, why, sol, L)` with the slice's checks (`check_boundary`, `admit`) when `checks`."""
function derive_law(m; checks = true)
    try
        if checks
            S.check_boundary(m)
            r = S.admit(m)
            r === nothing || return (false, "admit", nothing, nothing)
        end
        sol = S.solve(m)
        (true, "", sol, S.rate_law(m, sol))
    catch e
        (false, "error: " * first(replace(sprint(showerror, e), '\n' => ' '), 120), nothing,
         nothing)
    end
end

"""The slice's kcat class: `(:finite | :zero | :infinite | :nounique | :other, vars, terms)`."""
function kcat_class(m, sol)
    try
        vars, terms = S.kcat_candidates(m, sol)
        (all(t -> iszero(t[1]), terms) ? :zero : :finite, vars, terms)
    catch e
        msg = sprint(showerror, e)
        (occursin("never saturates", msg) ? :infinite :
         occursin("depends on how", msg) ? :nounique : :other, nothing, nothing)
    end
end

# ─── exact analysis of a law ─────────────────────────────────────────────────────────────

"""Concentration indices of substrates, products and modifiers (non-reactants)."""
function roles(m)
    names = S.metabolites(m)
    (names, [i for (i, n) in enumerate(names) if n in m.substrates],
     [i for (i, n) in enumerate(names) if n in m.products],
     [i for (i, n) in enumerate(names) if !S.is_reactant(m, n)])
end

"""`p` without the terms that hold one of the variables `zs` (those set to zero)."""
drop(p::S.Poly{T}, zs) where {T} =
    S.Poly{T}(Dict{S.Mono, T}(mo => c for (mo, c) in p.terms if !any(q -> q.first in zs, mo)))

"""Highest exponent of variable `v` in `p`; `nothing` for the zero polynomial."""
degv(p::S.Poly, v) = isempty(p.terms) ? nothing :
    maximum(mo -> sum((e for (w, e) in mo if w == v); init = 0 // 1), keys(p.terms))
"""Highest total exponent of the variables `vs` in `p`; `nothing` for zero."""
tdeg(p::S.Poly, vs) = isempty(p.terms) ? nothing :
    maximum(mo -> sum((e for (w, e) in mo if w in vs); init = 0 // 1), keys(p.terms))
allpos(p::S.Poly) = all(>(0), values(p.terms))

"""Degree bookkeeping with `nothing` for the zero polynomial: products add, sums take the
largest (all coefficients positive, so no cancellation)."""
dmul(xs...) = any(isnothing, xs) ? nothing : sum(xs)
dmax(xs...) = (ys = filter(!isnothing, collect(xs)); isempty(ys) ? nothing : maximum(ys))

"""
Exact zero-product facts of a law (plain or MWC), products and modifiers at zero:
`pos` (every coefficient of every component positive), `num0` and `den0` (numerator and
denominator keep a term), `gap[i]` = deg_X(num) - deg_X(den) for each substrate X alone
growing, `ugap` the same for all substrates growing together (total degree).
"""
function zero_facts(m, L)
    names, subs, prods, mods = roles(m)
    zs = [prods; mods]
    smult = [count(==(names[i]), m.substrates) for i in subs]       # exponent of X in ΠS
    if L isa S.Law
        Q0, D0 = drop(L.Q, zs), drop(L.D, zs)
        pos = allpos(L.Q) && allpos(L.D)
        dN(f) = dmul(f(Q0), 0)
        num_deg(v) = dmul(degv(Q0, v), smult[findfirst(==(v), subs)])
        gap = [isnothing(num_deg(v)) || isnothing(degv(D0, v)) ? nothing :
               num_deg(v) - degv(D0, v) for v in subs]
        ug = isempty(Q0.terms) || isempty(D0.terms) ? nothing :
             tdeg(Q0, subs) + sum(smult) - tdeg(D0, subs)
        return (pos = pos, num0 = !isempty(Q0.terms), den0 = !isempty(D0.terms), gap = gap,
                ugap = ug)
    end
    QI, DI = S.inactive_terms(L)
    comps = Dict(:QA => L.active.Q, :DA => L.active.D, :QI => QI, :DI => DI, :fA => L.fA,
                 :fI => L.fI)
    pos = all(allpos(p) || isempty(p.terms) for p in values(comps)) &&
          all(allpos(s[1]) && allpos(s[2]) for s in L.sites)
    z = Dict(k => drop(p, zs) for (k, p) in comps)
    n = L.n
    # site factors with the non-reactant ligands at zero (1 unless a ligand is a reactant)
    WA = [(drop(s[1], zs), s[3]) for s in L.sites]
    WI = [(drop(s[2], zs), s[3]) for s in L.sites]
    sc(x, k) = isnothing(x) ? nothing : k * x
    function degs(f)
        wa = dmul(0 // 1, (sc(f(w), k) for (w, k) in WA)...)
        wi = dmul(0 // 1, (sc(f(w), k) for (w, k) in WI)...)
        b1n = dmul(sc(f(z[:fI]), n), f(z[:QA]), n == 1 ? 0 // 1 : sc(f(z[:DA]), n - 1), wa)
        b2n = dmul(sc(f(z[:fA]), n), f(z[:QI]), n == 1 ? 0 // 1 : sc(f(z[:DI]), n - 1), wi)
        b1d = dmul(sc(f(z[:fI]), n), sc(f(z[:DA]), n), wa)
        b2d = dmul(sc(f(z[:fA]), n), sc(f(z[:DI]), n), wi)
        dmax(b1n, b2n), dmax(b1d, b2d)
    end
    nz(p) = isempty(p.terms) ? nothing : 0 // 1
    num0, den0 = degs(nz)
    gap = map(enumerate(subs)) do (j, v)
        a, b = degs(p -> degv(p, v))
        isnothing(a) || isnothing(b) ? nothing : a + smult[j] - b
    end
    a, b = degs(p -> tdeg(p, subs))
    (pos = pos, num0 = !isnothing(num0), den0 = !isnothing(den0), gap = gap,
     ugap = isnothing(a) || isnothing(b) ? nothing : a + sum(smult) - b)
end

# ─── numeric laws ────────────────────────────────────────────────────────────────────────

"""A polynomial in the concentrations with Float64 coefficients at one parameter point."""
struct NP
    c::Vector{Float64}
    v::Vector{Vector{Int}}
    e::Vector{Vector{Int}}
end
function NP(p::S.Poly, vars::S.Vars, pv::Vector{Float64})
    acc = Dict{Vector{Pair{Int, Int}}, Float64}()
    for (mo, c) in p.terms
        cm, val = Pair{Int, Int}[], Float64(c)
        for (v, e) in mo
            if vars.conc[v]
                push!(cm, v => Int(e))
            else
                val *= denominator(e) == 1 ? pv[v]^numerator(e) : pv[v]^Float64(e)
            end
        end
        acc[cm] = get(acc, cm, 0.0) + val
    end
    ks = collect(keys(acc))
    NP([acc[k] for k in ks], [first.(k) for k in ks], [last.(k) for k in ks])
end
function (p::NP)(x::AbstractVector{T}) where {T}
    s = zero(T)
    @inbounds for t in eachindex(p.c)
        y = T(p.c[t])
        for j in eachindex(p.v[t])
            y *= x[p.v[t][j]]^p.e[t][j]
        end
        s += y
    end
    s
end

"""Random parameter values by name: every non-concentration variable of `vars`, log-uniform
in [10^lo, 10^hi]."""
function draw_params(vars::S.Vars, rng; lo = -2.0, hi = 2.0)
    Dict{Symbol, Float64}(n => 10.0^(lo + (hi - lo) * rand(rng))
                          for (i, n) in enumerate(vars.names) if !vars.conc[i])
end
pvec(vars::S.Vars, pd) = Float64[vars.conc[i] ? 1.0 : pd[n] for (i, n) in enumerate(vars.names)]

"""A law at one parameter point: `v/E_total` of concentrations in `metabolites(m)` order."""
struct NL
    allo::Bool
    P::Vector{NP}            # plain: H, Q, D; MWC: H, QA, DA, QI, DI, fA, fI, WA..., WI...
    mult::Vector{Int}
    n::Int
    Lv::Float64
end
function NL(L, pd)
    if L isa S.Law
        pv = pvec(L.vars, pd)
        return NL(false, [NP(p, L.vars, pv) for p in (L.H, L.Q, L.D)], Int[], 1, 0.0)
    end
    pv = pvec(L.vars, pd)
    QI, DI = S.inactive_terms(L)
    P = [NP(p, L.vars, pv) for p in (L.H, L.active.Q, L.active.D, QI, DI, L.fA, L.fI)]
    append!(P, [NP(s[1], L.vars, pv) for s in L.sites])
    append!(P, [NP(s[2], L.vars, pv) for s in L.sites])
    NL(true, P, [s[3] for s in L.sites], L.n, pd[:L])
end
function rate(nl::NL, c::AbstractVector{T}) where {T}
    H = nl.P[1](c)
    nl.allo || return H * nl.P[2](c) / nl.P[3](c)
    QA, DA, QI, DI, fA, fI = (nl.P[i](c) for i in 2:7)
    ns = length(nl.mult)
    WA = prod((nl.P[7 + s](c)^nl.mult[s] for s in 1:ns); init = one(T))
    WI = prod((nl.P[7 + ns + s](c)^nl.mult[s] for s in 1:ns); init = one(T))
    n, L = nl.n, T(nl.Lv)
    H * (fI^n * QA * DA^(n - 1) * WA + L * fA^n * QI * DI^(n - 1) * WI) /
        (fI^n * DA^n * WA + L * fA^n * DI^n * WI)
end

"""kcat at parameter values `pd` from the candidates `(vars, terms)` (BigFloat)."""
function kcat_value(vars, terms, pd)
    x = BigFloat[vars.conc[i] ? big(1) : big(pd[n]) for (i, n) in enumerate(vars.names)]
    Float64(maximum(myeval(n, x) / myeval(d, x) for (n, d) in terms))
end

"""The value of a slice polynomial at `x` (by variable), rational exponents included."""
function myeval(p::S.Poly, x::Vector{T}) where {T}
    s = zero(T)
    for (mo, c) in p.terms
        t = T(c)
        for (v, e) in mo
            t *= denominator(e) == 1 ? x[v]^numerator(e) : x[v]^(T(numerator(e)) / denominator(e))
        end
        s += t
    end
    s
end

# ─── mass-action oracle (plain mechanisms) ───────────────────────────────────────────────

"""Every per-edge constant (`constant_names` slot order) of plain `m` at independent values
`pd` (by name, with :Keq), dependents from the solve's relations, in BigFloat."""
function edge_constants(m, sol, pd)
    names = S.column_names(m, sol)
    cv = Dict{Int, BigFloat}(i => big(pd[names[i]]) for i in sol.indep)
    bp(x, e) = denominator(e) == 1 ? x^numerator(e) : x^(big(numerator(e)) / denominator(e))
    for d in sol.deps
        cv[d.col] = prod((bp(cv[f], e) for (f, e) in d.terms); init = bp(big(pd[:Keq]), d.keq))
    end
    [[cv[sol.index[S.Col(k, s, false)]] for s in 1:(e.re ? 1 : 2)]
     for (k, e) in enumerate(m.edges)]
end

"""
The full mass-action steady state of plain `m` with each RE step running at speed `Λ`
(times `w[k]`): returns `(v, x)`, `v` the net consumption of the first substrate per
enzyme (per net reaction) and `x` the form populations (sum 1). `conc` by
`metabolites(m)` index; a modifier copy reads its reactant's concentration.
"""
function oracle(m, kv, conc::AbstractVector; Λ = big(10)^30, w = nothing)
    names = S.metabolites(m)
    ci = Dict(n => i for (i, n) in enumerate(names))
    fwd = S.orient(m)
    F = length(m.forms)
    W = zeros(BigFloat, F, F)
    s1 = m.substrates[1]
    flows = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
    for (k, e) in enumerate(m.edges)
        o = S.oriented(m, k, fwd)
        cons, rel = S.mets_of(m, o.cons), S.mets_of(m, o.rel)
        cc = prod((big(conc[ci[x.name]]) for x in cons); init = big(1))
        cr = prod((big(conc[ci[x.name]]) for x in rel); init = big(1))
        if !o.re
            kf, kr = kv[k][1], kv[k][2]
        else
            K, r = kv[k][1], Λ * (w === nothing ? 1 : w[k])
            kf, kr = S.is_binding(o) ? (r, r * K) : (r * K, r)
        end
        a, b = Int(o.a), Int(o.b)
        rab, rba = kf * cc, kr * cr
        W[b, a] += rab; W[a, a] -= rab; W[a, b] += rba; W[b, b] -= rba
        st = count(x -> x.name === s1 && !x.modifier, cons) -
             count(x -> x.name === s1 && !x.modifier, rel)
        push!(flows, (a, b, rab, rba, st))
    end
    A = copy(W)
    A[1, :] .= 1
    rhs = zeros(BigFloat, F)
    rhs[1] = 1
    x = A \ rhs
    v = sum(st * (rab * x[a] - rba * x[b]) for (a, b, rab, rba, st) in flows) /
        count(==(s1), m.substrates)
    v, x
end

# ─── hand-built mechanisms from DSL-like text ────────────────────────────────────────────

"""A form written `E`, `E(A, B)`, `F(B; +A -P)` or `E(S::Mod)`."""
function parse_form(t)
    t = strip(t)
    i = findfirst('(', t)
    i === nothing && return S.RawForm(Symbol(t), S.Met[], Symbol[], Symbol[])
    conf, inner = Symbol(t[1:i-1]), t[i+1:end-1]
    parts = split(inner, ';')
    bound = S.Met[]
    for x in split(parts[1], ',')
        x = strip(x)
        isempty(x) && continue
        push!(bound, endswith(x, "::Mod") ? S.Met(Symbol(x[1:end-5]), true) :
                     S.Met(Symbol(x), false))
    end
    added, sub = Symbol[], Symbol[]
    if length(parts) > 1
        for x in split(strip(parts[2]))
            x[1] == '+' ? push!(added, Symbol(x[2:end])) : push!(sub, Symbol(x[2:end]))
        end
    end
    S.RawForm(conf, bound, added, sub)
end

isform(t, confs) = (u = strip(t); any(c -> u == c || startswith(u, c * "("), confs))

"""
A plain mechanism from step strings `"E(A) + B <--> E(A, B)"` (`<-->` steady state, `<=>`
rapid equilibrium), optionally `" [g3]"` for the kinetic group (default: one per step).
"""
function mech(subs, prods, lines; confs = ("E", "F", "G"))
    steps, group = S.RawStep[], Int[]
    for (j, l) in enumerate(lines)
        g = j
        mt = match(r"\[g(\d+)\]", l)
        if mt !== nothing
            g = parse(Int, mt[1])
            l = replace(l, mt.match => "")
        end
        re = occursin("<=>", l)
        lhs, rhs = strip.(split(l, re ? "<=>" : "<-->"))
        sideL, sideR = split_side(lhs), split_side(rhs)
        fL = only(filter(t -> isform(t, confs), sideL))
        fR = only(filter(t -> isform(t, confs), sideR))
        met(x) = endswith(x, "::Mod") ? S.Met(Symbol(x[1:end-5]), true) : S.Met(Symbol(x), false)
        cons = [met(strip(t)) for t in sideL if !isform(t, confs)]
        rel = [met(strip(t)) for t in sideR if !isform(t, confs)]
        push!(steps, S.RawStep(parse_form(fL), parse_form(fR), cons, rel, re))
        push!(group, g)
    end
    S.Mechanism(Symbol.(subs), Symbol.(prods), steps, relabel_groups(group))
end
relabel_groups(g) = (u = unique(g); [findfirst(==(x), u) for x in g])

"""Splits one side of a step at the top-level `+` signs (not inside parentheses)."""
function split_side(s)
    out, depth, cur = String[], 0, IOBuffer()
    for ch in s
        if ch == '(' ; depth += 1; elseif ch == ')'; depth -= 1; end
        if ch == '+' && depth == 0
            push!(out, strip(String(take!(cur))))
        else
            write(cur, ch)
        end
    end
    push!(out, strip(String(take!(cur))))
    filter(!isempty, out)
end
