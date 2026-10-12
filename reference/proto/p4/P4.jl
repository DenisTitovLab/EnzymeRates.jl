# ABOUTME: P4 probe library (module P4): joint-solved laws with readable names, constant-term criteria,
# ABOUTME: coefficient maps, units, the Hill factor, rational functions and Cleland coefficient ratios.
isdefined(Main, :Pop) || include(joinpath(@__DIR__, "..", "common", "Population.jl"))
isdefined(Main, :Poly) || include(joinpath(@__DIR__, "..", "common", "Poly.jl"))
isdefined(Main, :Engine) || include(joinpath(@__DIR__, "..", "common", "Engine.jl"))
isdefined(Main, :Thermo) || include(joinpath(@__DIR__, "..", "common", "Thermo.jl"))

module P4

import ..Pop, ..Engine, ..Thermo
using ..Poly: Poly, Vars, var!, Mono, mono, mmul, mpow, mexp, P, term, subst, lowest_shift,
              nterms, scale, expand, addto!

const PI = P{Int}

toint(q) = (denominator(q) == 1 || error("fractional exponent $q"); Int(numerator(q)))

# ─── Solves ────────────────────────────────────────────────────────────

"""Thermo's joint solve with a choice of cycle rows: `:all` (every cycle, IndependentParameters),
`:re` (cycles of the RE subgraph only: AllParameters with the RE closure) or `:none` (kinetic-group
rows only: AllParameters as today's Full)."""
function solve(m::Pop.PMech; cycles::Symbol = :all, front = nothing)
    front === nothing && cycles === :all && return Thermo.solve(m)
    front === nothing && cycles === :none && return Thermo.solve(m; mode = :all)
    T = Thermo
    cs = T.canonical_steps(m)
    cols = T.columns(cs)
    front === nothing || (cols = [filter(c -> c.name in front, cols); filter(c -> !(c.name in front), cols)])
    colidx = Dict((c.j, c.r) => i for (i, c) in enumerate(cols))
    rows, rhs = Vector{Int}[], Int[]
    if cycles === :all
        deltas = T.step_delta(cs)
        for z in T.fundamental_cycles(length(m.species), cs)
            push!(rows, T.cycle_row(cs, colidx, z)); push!(rhs, T.cycle_multiple(m, deltas, z))
        end
    elseif cycles === :re
        re = findall(c -> c.re, cs)
        for zs in T.fundamental_cycles(length(m.species), cs[re])
            z = zeros(Int, length(cs)); z[re] .= zs
            push!(rows, T.cycle_row(cs, colidx, z)); push!(rhs, 0)
        end
    end
    for r in T.group_rows(cs, colidx, T.group_alignment(m, cs, m.groups))
        push!(rows, r); push!(rhs, 0)
    end
    n = length(cols)
    M = zeros(T.Q, length(rows), n + 1)
    for (i, r) in enumerate(rows); M[i, 1:n] .= r; M[i, n + 1] = rhs[i]; end
    piv, big = T.rref!(M)
    indep = setdiff(1:n, piv)
    (n + 1) in piv && error("$(m.id): the constraints reduce to 0 = log Keq")
    deps = [T.Dep(p, M[r, n + 1], [f => -M[r, f] for f in indep if M[r, f] != 0],
                  M[r, n + 1] != 0 ? :haldane : :wegscheider, false) for (r, p) in enumerate(piv)]
    T.Solution(m, cs, cols, indep, deps, length(rows), 0, 0, big)
end

# ─── Laws over readable names ──────────────────────────────────────────

"""A rate law `v/E_t = num/den` over named variables: concentrations (`conc`), independent
constants named by the package's group names in the package's stored convention, and `Keq`
(`keq`, 0 when absent). `units[v]` is (concentration exponent, time exponent). `stepmono[j]`
holds each step constant in Engine's convention (`[K]` or `[kf, kr]`, stored direction) as a
monomial; `colvar[i]` the variable of independent column `i`, `sgn[i]` its exponent in that
column's value."""
struct NLaw
    m::Pop.PMech
    sol::Thermo.Solution
    vars::Vars
    conc::Dict{Symbol, Int}
    keq::Int
    num::PI
    den::PI
    units::Dict{Int, NTuple{2, Int}}
    stepmono::Vector{Vector{Mono}}
    colvar::Dict{Int, Int}
    colsgn::Dict{Int, Int}
end

"""The law of `m` with the solve of `cycles` substituted (see `solve`)."""
function law(m::Pop.PMech; cycles::Symbol = :all, front = nothing,
             L = Engine.bare_rate_law(m; check = false))
    sol = solve(m; cycles, front)
    W = Vars()
    conc = Dict(x => var!(W, x) for x in m.metabolites)
    units = Dict{Int, NTuple{2, Int}}(v => (1, 0) for v in values(conc))
    kv = cycles === :all ? var!(W, :Keq) : 0
    kv == 0 || (units[kv] = (length(m.products) - length(m.substrates), 0))
    sgn(j) = (sol.cs[j].re && !sol.cs[j].bindK && sol.cs[j].flipped) ? -1 : 1
    slot(c) = sol.cs[c.j].re ? 1 : (sol.cs[c.j].flipped ? 3 - c.r : c.r)
    colmono, colvar, colsgn = Dict{Int, Mono}(), Dict{Int, Int}(), Dict{Int, Int}()
    for i in sol.indep
        c = sol.cols[i]; s = m.steps[c.j]
        name = m.step_names[c.j][slot(c)]
        haskey(W.index, name) && error("$(m.id): $(name) names two independents")
        v = var!(W, name)
        units[v] = s.re ? (Pop.is_binding(s) ? (1, 0) : (length(s.released) - length(s.consumed), 0)) :
                   (-(slot(c) == 1 ? length(s.consumed) : length(s.released)), -1)
        colmono[i] = mono(v => sgn(c.j)); colvar[i] = v; colsgn[i] = sgn(c.j)
    end
    for d in sol.deps
        d.keq != 0 && kv == 0 && error("$(m.id): Keq in a dependent of a Keq-free solve")
        mo = d.keq == 0 ? Mono() : mono(kv => toint(d.keq))
        for (f, e) in d.terms
            mo = mmul(mo, mpow(colmono[f], toint(e)))
        end
        colmono[d.col] = mo
    end
    colidx = Dict((c.j, c.r) => i for (i, c) in enumerate(sol.cols))
    stepmono = map(enumerate(m.steps)) do (j, s)
        s.re ? [mpow(colmono[colidx[(j, 1)]], sgn(j))] :
        [colmono[colidx[(j, sol.cs[j].flipped ? 2 : 1)]], colmono[colidx[(j, sol.cs[j].flipped ? 1 : 2)]]]
    end
    f = Dict{Int, Mono}()
    for (j, vs) in enumerate(L.stepvars), (r, v) in enumerate(vs)
        f[v] = stepmono[j][r]
    end
    for v in L.concvars
        f[v] = mono(conc[L.vars.names[v]] => 1)
    end
    num, den = subst(L.num, v -> f[v]), subst(L.den, v -> f[v])
    shift = lowest_shift((num, den), values(conc))
    NLaw(m, sol, W, conc, kv, num * shift, den * shift, units, stepmono, colvar, colsgn)
end

# ─── Coefficient maps ──────────────────────────────────────────────────

"""`p` split by concentration monomial: concentration monomial => parameter polynomial."""
function cmap(p::PI, concvars)
    cv = Set(concvars)
    out = Dict{Mono, PI}()
    for (mo, c) in p.t
        cm = Mono([q for q in mo if q.first in cv])
        pm = Mono([q for q in mo if !(q.first in cv)])
        addto!(get!(() -> PI(), out, cm).t, pm, c)
    end
    out
end
cmap(L::NLaw, p::PI) = cmap(p, values(L.conc))

"""Concentration monomial of the metabolites `xs` (each once)."""
cmono(L::NLaw, xs) = foldl(mmul, (mono(L.conc[x] => 1) for x in xs); init = Mono())

degree(mo::Mono) = sum((e for (_, e) in mo); init = 0)

# ─── Units ─────────────────────────────────────────────────────────────

unit(L::NLaw, mo::Mono) = foldl(((a, b), (v, e)) -> (a + e * L.units[v][1], b + e * L.units[v][2]),
                                mo; init = (0, 0))
"""The set of units of the terms of `p`."""
units_of(L::NLaw, p::PI) = Set(unit(L, mo) for mo in keys(p.t))

"""Whether every denominator term has one unit and every numerator term that unit times 1/time."""
function homogeneous(L::NLaw)
    ud, un = units_of(L, L.den), units_of(L, L.num)
    length(ud) == 1 && length(un) == 1 && only(un) == only(ud) .+ (0, -1)
end

# ─── Constant-term criteria ────────────────────────────────────────────

"""The terminal strongly connected components of the digraph on `1:n` with edges `a => b`."""
function terminal_sccs(n, es)
    adj = [Int[] for _ in 1:n]
    for (a, b) in es
        a == b || push!(adj[a], b)
    end
    reach = map(1:n) do s
        seen = falses(n); seen[s] = true; st = [s]
        while !isempty(st)
            u = pop!(st)
            for w in adj[u]
                seen[w] || (seen[w] = true; push!(st, w))
            end
        end
        seen
    end
    comps, done = Vector{Int}[], falses(n)
    for i in 1:n
        done[i] && continue
        c = [j for j in 1:n if reach[i][j] && reach[j][i]]
        done[c] .= true; push!(comps, c)
    end
    [c for c in comps if all(j -> !reach[c[1]][j] || j in c, 1:n)]
end

"""Directions of step `s` that take up no metabolite, as form pairs."""
free_dirs(s::Pop.PStep) = vcat(isempty(s.consumed) ? [s.from => s.to] : Pair{Int, Int}[],
                               isempty(s.released) ? [s.to => s.from] : Pair{Int, Int}[])

"""RE segment of each form (numbered 1:G) and the segments' forms."""
function segments(m::Pop.PMech)
    n = length(m.species)
    p = collect(1:n)
    root(i) = p[i] == i ? i : (p[i] = root(p[i]))
    for s in m.steps
        s.re && (p[root(s.from)] = root(s.to))
    end
    roots = unique(root(i) for i in 1:n)
    seg = [findfirst(==(root(i)), roots) for i in 1:n]
    seg, [findall(==(g), seg) for g in eachindex(roots)]
end

"""Terminal-SCC predictions of a constant term:
- `c0`: the full form graph, every step (RE included) a pair of directed edges;
- `c1`: the Cha-contracted graph, the metabolite-free directions of SS steps between segments;
- `c2`: `c1` with a segment's out-edge kept only from a form that every form of its segment
  reaches by metabolite-free RE directions (an exit form), and no prediction of a constant term
  when some segment has no exit form.
Each entry is the number of terminal SCCs (`c2` = 0 when a segment has no exit form), plus
whether the free enzyme's form/segment is terminal."""
function criteria(m::Pop.PMech)
    n = length(m.species)
    seg, segs = segments(m)
    G = length(segs)
    t0 = terminal_sccs(n, [d for s in m.steps for d in free_dirs(s)])
    ss_dirs = [d for s in m.steps if !s.re for d in free_dirs(s) if seg[d.first] != seg[d.second]]
    t1 = terminal_sccs(G, [seg[a] => seg[b] for (a, b) in ss_dirs])
    exits = map(segs) do F
        es = [d for s in m.steps if s.re for d in free_dirs(s) if seg[d.first] == seg[F[1]]]
        loc = Dict(f => i for (i, f) in enumerate(F))
        tl = terminal_sccs(length(F), [loc[a] => loc[b] for (a, b) in es])
        length(tl) == 1 ? Set(F[i] for i in only(tl)) : Set{Int}()
    end
    t2 = any(isempty, exits) ? Vector{Int}[] :
         terminal_sccs(G, [seg[a] => seg[b] for (a, b) in ss_dirs if a in exits[seg[a]]])
    # c3: a segment's zero-concentration limit is its equilibrium distribution: the forms of
    # least net uptake (α exponent vector componentwise minimal); none when no form attains the
    # componentwise minimum (an exchange step such as an RE Theorell-Chance step)
    up = uptake(m, seg, segs)
    mins = map(F -> Dict(x => minimum(get(up[i], x, 0) for i in F) for x in m.metabolites), segs)
    low = [Set(i for i in F if all(get(up[i], x, 0) == mins[g][x] for x in m.metabolites))
           for (g, F) in enumerate(segs)]
    t3 = any(isempty, low) ? Vector{Int}[] :
         terminal_sccs(G, [seg[a] => seg[b] for (a, b) in ss_dirs if a in low[seg[a]]])
    ifree = findfirst(Pop.is_free, m.species)
    (c0 = length(t0), c1 = length(t1), c2 = length(t2), c3 = length(t3),
     free_terminal_c3 = ifree !== nothing && length(t3) == 1 && seg[ifree] in only(t3),
     segments = G, exitless = count(isempty, exits), no_minimal = count(isempty, low),
     nonzero_min = count(d -> any(!=(0), values(d)), mins))
end

"""Net metabolite uptake (by concentration name) of each form from its segment's root along
RE steps: the concentration exponents of the form's Cha weight α."""
function uptake(m::Pop.PMech, seg, segs)
    hasres(sp) = !isempty(sp.res_added) || !isempty(sp.res_subtracted)
    up = Vector{Union{Nothing, Dict{Symbol, Int}}}(nothing, length(m.species))
    for F in segs
        r = argmin(i -> (length(m.species[i].bound), hasres(m.species[i]) ? 1 : 0,
                         string(m.species[i].name)), F)
        up[r] = Dict{Symbol, Int}(); queue = [r]
        while !isempty(queue)
            u = popfirst!(queue)
            for s in m.steps
                s.re || continue
                for (a, b, sg) in ((s.from, s.to, 1), (s.to, s.from, -1))
                    (a == u && up[b] === nothing) || continue
                    d = copy(up[a])
                    for x in s.consumed; d[x.name] = get(d, x.name, 0) + sg; end
                    for x in s.released; d[x.name] = get(d, x.name, 0) - sg; end
                    up[b] = d; push!(queue, b)
                end
            end
        end
    end
    up
end

"""The Engine law's denominator before its lowest-terms shift, and the shift per metabolite."""
function den_unshifted(L::Engine.Law)
    σ, F, A = L.stats[:sigma], L.stats[:factors], L.stats[:attach]
    d = sum(σ[g] * expand(PI, [F[b][A[b][g]] for b in eachindex(F)]) for g in eachindex(σ))
    minexp(p, v) = minimum(mexp(t, v) for t in keys(p.t))
    shift = Dict(L.vars.names[v] => minexp(L.den, v) - minexp(d, v) for v in L.concvars)
    d * mono((v => shift[L.vars.names[v]] for v in L.concvars)...) == L.den ||
        error("unshifted denominator does not reproduce the law")
    d, shift
end

# ─── Hill factor ───────────────────────────────────────────────────────

"""Exact quotient of `p` by `1 − x` for a monomial `x` holding `kv` (Keq) to the power −1, by
chain division (P1's `divide_one_minus`); `nothing` when `1 − x` does not divide `p`."""
function divide_one_minus(p::PI, x::Mono, kv::Int)
    chains = Dict{Mono, Vector{Pair{Int, Int}}}()
    for (t, c) in p.t
        e = mexp(t, kv)
        push!(get!(chains, mmul(t, mpow(x, e)), Pair{Int, Int}[]), -e => c)
    end
    q = PI()
    for (key, terms) in chains
        sort!(terms; by = first)
        acc = 0
        for (k, (i, c)) in enumerate(terms)
            acc += c
            k == length(terms) && break
            for j in i:(terms[k + 1].first - 1)
                addto!(q.t, mmul(key, mpow(x, j)), acc)
            end
        end
        iszero(acc) || return nothing
    end
    q
end

"""`num = (ΠS − ΠP/Keq)·Q`: `Q` (or `nothing`), whether it is a polynomial in the
concentrations, has positive coefficients, and is free of Keq."""
function hill(L::NLaw)
    S, Pr = cmono(L, L.m.substrates), cmono(L, L.m.products)
    q = divide_one_minus(L.num, mmul(mmul(Pr, mpow(S, -1)), mono(L.keq => -1)), L.keq)
    q === nothing && return (Q = nothing, poly = false, positive = false, keq_free = false)
    Q = q * mpow(S, -1)
    (Q = Q, poly = all(mexp(t, v) >= 0 for t in keys(Q.t) for v in values(L.conc)),
     positive = all(>(0), values(Q.t)), keq_free = all(mexp(t, L.keq) == 0 for t in keys(Q.t)))
end

"""Class of a dependent column: `(re, kind, against)` (kind 3 iso, 2 chemistry, 1 binding;
`against` = the SS constant drives against the physical forward direction)."""
depclass(sol, d) = (c = sol.cols[d.col]; (sol.cs[c.j].re, sol.cs[c.j].kind, c.prio[4] == 1))

# ─── Printing ──────────────────────────────────────────────────────────

function fmt_mono(mo::Mono, W::Vars)
    nu, de = String[], String[]
    for (v, e) in sort(mo; by = q -> string(W.names[q.first]))
        s = string(W.names[v])
        push!(e > 0 ? nu : de, abs(e) == 1 ? s : string(s, "^", abs(e)))
    end
    n = isempty(nu) ? "1" : join(nu, "*")
    isempty(de) ? n : string(n, "/", length(de) == 1 ? de[1] : "(" * join(de, "*") * ")")
end

function fmt_poly(p::PI, W::Vars)
    iszero(p) && return "0"
    ts = sort(collect(p.t); by = t -> (degree(t[1]), fmt_mono(t[1], W)))
    s = join([(c == 1 ? "" : c == -1 ? "-" : string(c, "*")) * fmt_mono(mo, W) for (mo, c) in ts],
             " + ")
    replace(s, "+ -" => "- ")
end

# ─── Exact division and display ratios ─────────────────────────────────

"""Exact quotient `p ÷ q` of Laurent polynomials, or `nothing` when `q` does not divide `p` over
the integers (both shifted to polynomials, then lex-order division)."""
function exact_div(p::PI, q::PI)
    iszero(q) && error("division by zero")
    vs = sort!(unique(v for x in (p, q) for t in keys(x.t) for (v, _) in t))
    sp, sq = lowest_shift((p,), vs), lowest_shift((q,), vs)
    p, q = p * sp, q * sq
    key(t) = [mexp(t, v) for v in vs]
    lq = argmax(key, collect(keys(q.t))); cq = q.t[lq]
    r, out = p, PI()
    while !iszero(r)
        lr = argmax(key, collect(keys(r.t))); cr = r.t[lr]
        t = mmul(lr, mpow(lq, -1))
        (rem(cr, cq) == 0 && all(e >= 0 for (_, e) in t)) || return nothing
        addto!(out.t, t, div(cr, cq))
        r = r - term(PI, t, div(cr, cq)) * q
    end
    out * mmul(sq, mpow(sp, -1))
end

"""Factors of `p` that are sums of single variables, found as groups of variables that never
share a term (exact when `p` is a monomial times a product of such linear forms)."""
function linear_factors(p::PI)
    vs = unique(v for t in keys(p.t) for (v, _) in t)
    together = Set(minmax(a, b) for t in keys(p.t) for (a, _) in t for (b, _) in t if a != b)
    comps = Vector{Int}[]
    for v in vs
        k = findfirst(c -> all(w -> !(minmax(v, w) in together), c), comps)
        k === nothing ? push!(comps, [v]) : push!(comps[k], v)
    end
    [f for f in (PI(Dict(mono(v => 1) => 1 for v in c)) for c in comps if length(c) > 1)
     if exact_div(p, f) !== nothing]
end

"""`a/b` with common monomial content and common linear factors cancelled."""
function reduce_ratio(a::PI, b::PI)
    vs = unique(v for x in (a, b) for t in keys(x.t) for (v, _) in t)
    s = lowest_shift((a, b), vs)
    a, b = a * s, b * s
    q = exact_div(a, b)
    q !== nothing && return q, one(PI)
    for f in linear_factors(b)
        qa = exact_div(a, f)
        qa === nothing || (a = qa; b = exact_div(b, f))
    end
    a, b
end

"""`p` printed as monomial content times its linear factors times the remaining factor."""
function fmt_factored(p::PI, W::Vars)
    vs = unique(v for t in keys(p.t) for (v, _) in t)
    s = lowest_shift((p,), vs)
    c = mpow(s, -1)
    rest, facs = p * s, PI[]
    for f in linear_factors(rest)
        q = exact_div(rest, f)
        q === nothing || (push!(facs, f); rest = q)
    end
    parts = String[]
    isempty(c) || push!(parts, fmt_mono(c, W))
    append!(parts, ["(" * fmt_poly(f, W) * ")" for f in facs])
    rest == one(PI) || push!(parts, nterms(rest) == 1 ? fmt_poly(rest, W) : "(" * fmt_poly(rest, W) * ")")
    isempty(parts) ? "1" : join(parts, "*")
end

# ─── Rational functions ────────────────────────────────────────────────

"""A rational function `n/d` of Laurent polynomials, kept with no common monomial factor."""
struct RF
    n::PI
    d::PI
end
function rf(n::PI, d::PI)
    iszero(d) && error("division by zero")
    iszero(n) && return RF(n, one(PI))
    vs = unique(v for p in (n, d) for t in keys(p.t) for (v, _) in t)
    s = lowest_shift((n, d), vs)
    RF(n * s, d * s)
end
RF(p::PI) = rf(p, one(PI))
RF(x::Integer) = RF(scale(one(PI), x), one(PI))
RF(x::Rational) = RF(scale(one(PI), numerator(x)), scale(one(PI), denominator(x)))
RF(x::AbstractFloat) = RF(Rational{Int}(x))
Base.:+(a::RF, b::RF) = a.d == b.d ? rf(a.n + b.n, a.d) : rf(a.n * b.d + b.n * a.d, a.d * b.d)
Base.:-(a::RF) = RF(-a.n, a.d)
Base.:-(a::RF, b::RF) = a + (-b)
Base.:*(a::RF, b::RF) = rf(a.n * b.n, a.d * b.d)
Base.:/(a::RF, b::RF) = rf(a.n * b.d, a.d * b.n)
Base.:^(a::RF, k::Integer) = k >= 0 ? foldl(*, fill(a, k); init = RF(1)) : RF(1) / a^(-k)
for op in (:+, :-, :*, :/)
    @eval Base.$op(a::RF, b::Number) = $op(a, RF(b))
    @eval Base.$op(a::Number, b::RF) = $op(RF(a), b)
end
"""Whether two rational functions are equal (cross-multiplied, exact)."""
same(a::RF, b::RF) = a.n * b.d == b.n * a.d

# ─── Cleland constants ─────────────────────────────────────────────────

"""Cleland's kinetic constants of a bi-bi law as coefficient ratios (Cleland 1963; the ordered
Bi Bi ratios as in Segel p. 564 and DynaFit TN-2015-03). `A, B` are the substrates in binding
order and `P, Q` the products in release order. An inhibition constant is c(T)/c(T·X) for the
first context T, in a fixed list, whose two terms are both present: Kia (∅, then P), Kib (PQ,
then ∅), Kip (AB, then ∅), Kiq (∅, then B). `nothing` marks a constant the law does not
define."""
function cleland_bibi(L::NLaw; A, B, P, Q)
    N, D = cmap(L, L.num), cmap(L, L.den)
    c(xs...) = get(D, cmono(L, xs), PI())
    n(xs...) = get(N, cmono(L, xs), PI())
    r(a, b) = iszero(b) ? nothing : iszero(a) ? RF(0) : RF(a) / RF(b)
    ki(pairs...) = (for (t, tx) in pairs
                        (iszero(c(t...)) || iszero(c(tx...))) || return r(c(t...), c(tx...))
                    end; nothing)
    (V1 = r(n(A, B), c(A, B)), V2 = r(-n(P, Q), c(P, Q)), Ka = r(c(B), c(A, B)),
     Kb = r(c(A), c(A, B)), Kp = r(c(Q), c(P, Q)), Kq = r(c(P), c(P, Q)),
     Kia = ki(((), (A,)), ((P,), (A, P))), Kib = ki(((P, Q), (B, P, Q)), ((), (B,))),
     Kip = ki(((A, B), (A, B, P)), ((), (P,))), Kiq = ki(((), (Q,)), ((B,), (B, Q))),
     Keq = r(n(A, B), -n(P, Q)))
end

"""Cleland's uni-uni constants: V1, V2, Km of A (Ka) and of P (Kp), Keq."""
function cleland_uniuni(L::NLaw; A, P)
    N, D = cmap(L, L.num), cmap(L, L.den)
    c(xs...) = get(D, cmono(L, xs), PI())
    n(xs...) = get(N, cmono(L, xs), PI())
    r(a, b) = iszero(b) ? nothing : iszero(a) ? RF(0) : RF(a) / RF(b)
    (V1 = r(n(A), c(A)), V2 = r(-n(P), c(P)), Ka = r(c(), c(A)), Kp = r(c(), c(P)),
     Keq = r(n(A), -n(P)))
end

end # module
