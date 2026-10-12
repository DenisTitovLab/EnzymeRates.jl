# ABOUTME: P7 probe (module P7): positional rate-law IRs of Pop.PMech records, flat (the package's
# ABOUTME: expanded polynomials) and factored (Haldane factor, segment and block factors, MWC nodes).
const _COMMON = joinpath(@__DIR__, "..", "common")
isdefined(Main, :Pop) || include(joinpath(_COMMON, "Population.jl"))
isdefined(Main, :Poly) || include(joinpath(_COMMON, "Poly.jl"))
isdefined(Main, :Engine) || include(joinpath(_COMMON, "Engine.jl"))
isdefined(Main, :EngineLean) || include(joinpath(@__DIR__, "..", "p1_engine", "EngineLean.jl"))

module P7

import ..Pop, ..Poly, ..Engine, ..EngineLean
using ..Poly: Vars, var!, Mono, mono, mmul, mpow, mexp, P, term, subst, lowest_shift, nterms

const PI = P{Int}

"""`coef · Π var^exp · Π node^pow`."""
struct Term
    coef::Int
    mono::Mono
    facs::Vector{Pair{Int, Int}}
end

"""A rate law over positional variables: parameters `1:np` (fitted…, Keq, E_total), then the
dependents (`deps[i]` defines variable `np + i` as a power product of parameters), then the
concentrations. `v = x[scale] · node[num] / node[den]`; a node is a sum of `Term`s and may refer
to earlier nodes."""
struct Law
    vars::Vars
    np::Int
    deps::Vector{Mono}
    nc::Int
    nodes::Vector{Vector{Term}}
    num::Int
    den::Int
    scale::Int
end

# ─── Variables ────────────────────────────────────────────────────────────

"""The variable table of `m`: fitted parameters, Keq, E_total, dependents (by name), concentrations."""
function law_vars(m::Pop.PMech)
    V = Vars()
    foreach(s -> var!(V, s), [m.fitted_params; :Keq; :E_total])
    np = length(V.names)
    foreach(d -> var!(V, first(d)), m.dependent_monomials)
    foreach(s -> var!(V, s), m.metabolites)
    V, np
end
isconc(V, np, nd, v) = v > np + nd
dep_monos(m, V) = Mono[mono((V.index[s] => Int(e) for (s, e) in pp)...) for (_, pp) in m.dependent_monomials]

"""A `Pop.PPoly` over the names of `V` (integer coefficients)."""
topoly(pp, V) = PI(Dict{Mono, Int}(mono((V.index[s] => e for (s, e) in mn)...) => Int(c) for (mn, c) in pp))

"""`p` with every dependent variable replaced by its power product (Laurent in the parameters)."""
substdeps(p::PI, np, deps) = subst(p, v -> np < v <= np + length(deps) ? deps[v - np] : nothing)

node_of(p::PI) = sort!([Term(c, mn, Pair{Int, Int}[]) for (mn, c) in p.t]; by = t -> (t.mono, t.coef))

# ─── Flat laws: the package's expanded polynomials ───────────────────────

"""The flat law of `m`: the package's own expanded polynomials, dependents kept as variables.
An allosteric mechanism assembles the MWC sum from its two states' polynomials as
`_num_den_exprs` does (state-polynomial nodes are shared between the numerator and denominator)."""
function flat_law(m::Pop.PMech)
    V, np = law_vars(m)
    deps = dep_monos(m, V)
    m.allosteric || return Law(V, np, deps, length(m.metabolites),
                               [node_of(topoly(m.pkg_num, V)), node_of(topoly(m.pkg_den, V))], 1, 2,
                               V.index[:E_total])
    mwc_law(m, V, np, deps, nothing)
end

# ─── MWC assembly ─────────────────────────────────────────────────────────

"""Regulatory-site factor of `site` in state `st` (`:A` or `:I`) as a polynomial, `nothing` when
the site holds no ligand in that state (its factor is 1)."""
function reg_factor(site::Pop.PSite, st::Symbol, V)
    p = PI(Dict{Mono, Int}(Mono() => 1))
    for (lig, tag) in zip(site.ligands, site.states)
        tag === (st === :I ? :OnlyA : :OnlyI) && continue
        K = Symbol("K_", st === :I && tag !== :EqualAI ? "I" : "A", "_", lig, "reg")
        p = p + term(PI, mono(V.index[lig] => 1, V.index[K] => -1))
    end
    nterms(p) == 1 ? nothing : p
end

"""The MWC law of an allosteric `m`. `fac === nothing` keeps each state numerator as one node;
otherwise `fac(N)` returns `(H, Q)` with `N == H·Q` (the Haldane factor, shared) or `nothing`."""
function mwc_law(m::Pop.PMech, V, np, deps, fac)
    n = m.cat_multiplicity
    NA, QA, DA = (topoly(p, V) for p in m.state_polys_A)
    NI, QI, DI = (topoly(p, V) for p in m.state_polys_I)
    concs = Set(np + length(deps) + 1:length(V.names))
    monofree(p) = nterms(p) == 1 && !any(v in concs for (v, _) in first(keys(p.t)))
    cross = false
    if DA == DI
    elseif monofree(DA) && monofree(DI)
        ia, ii = mpow(first(keys(DA.t)), -1), mpow(first(keys(DI.t)), -1)
        all(==(1), values(DA.t)) && all(==(1), values(DI.t)) || error("dfree coefficient")
        NA, QA, NI, QI = NA * ia, QA * ia, NI * ii, QI * ii
    else
        cross = true
    end
    nodes = Vector{Term}[]
    add!(p) = (push!(nodes, node_of(p)); length(nodes))
    L = V.index[:L]
    fA = fac === nothing ? nothing : fac(NA)
    fI = fac === nothing || iszero(NI) ? nothing : fac(NI)
    shareH = fA !== nothing && (iszero(NI) || fI !== nothing)
    H = shareH ? add!(fA[1]) : 0
    nA = shareH ? add!(fA[2]) : add!(NA)
    qA = add!(QA)
    nI = iszero(NI) ? 0 : shareH ? add!(fI[2]) : add!(NI)
    qI = add!(QI)
    rA = [add!(r) for r in (reg_factor(s, :A, V) for s in m.sites) if r !== nothing]
    rI = [add!(r) for r in (reg_factor(s, :I, V) for s in m.sites) if r !== nothing]
    mA = [m.sites[i].multiplicity for i in eachindex(m.sites) if reg_factor(m.sites[i], :A, V) !== nothing]
    mI = [m.sites[i].multiplicity for i in eachindex(m.sites) if reg_factor(m.sites[i], :I, V) !== nothing]
    dA = cross ? add!(DA) : 0
    dI = cross ? add!(DI) : 0
    wA = [rA .=> mA; cross ? [dI => n] : Pair{Int, Int}[]]       # factors of the A-state terms
    wI = [rI .=> mI; cross ? [dA => n] : Pair{Int, Int}[]]
    pw(q, k) = k == 0 ? Pair{Int, Int}[] : [q => k]
    numt = [Term(1, Mono(), [nA => 1; pw(qA, n - 1); wA])]
    iszero(NI) || push!(numt, Term(1, mono(L => 1), [nI => 1; pw(qI, n - 1); wI]))
    if shareH
        push!(nodes, numt)
        numt = [Term(1, Mono(), [H => 1, length(nodes) => 1])]
    end
    push!(nodes, numt)
    num = length(nodes)
    push!(nodes, [Term(1, Mono(), [qA => n; wA]), Term(1, mono(L => 1), [qI => n; wI])])
    Law(V, np, deps, length(m.metabolites), nodes, num, length(nodes), V.index[:E_total])
end

# ─── Factoring tools ──────────────────────────────────────────────────────

"""`(ΠS − ΠP/Keq, N/(ΠS − ΠP/Keq))` for a numerator `N` over the parameters, Keq and the
concentrations (dependents substituted), or `nothing` when the Haldane factor does not divide."""
function haldane_split(N::PI, m::Pop.PMech, V)
    kv = V.index[:Keq]
    S = mono((V.index[s] => 1 for s in m.substrates)...)
    Pm = mono((V.index[s] => 1 for s in m.products)...)
    x = mmul(mmul(Pm, mpow(S, -1)), mono(kv => -1))             # Γ/Keq
    Q = EngineLean.divide_one_minus(N, x, kv)
    Q === nothing && return nothing
    H = term(PI, S) - term(PI, mmul(Pm, mono(kv => -1)))
    H, Q * mpow(S, -1)
end

"""Total order on Laurent monomials (lexicographic on exponent vectors)."""
function mless(a::Mono, b::Mono)
    i, j = 1, 1
    while i <= length(a) || j <= length(b)
        va = i <= length(a) ? a[i].first : typemax(Int)
        vb = j <= length(b) ? b[j].first : typemax(Int)
        v = min(va, vb)
        ea = va == v ? a[i].second : 0
        eb = vb == v ? b[j].second : 0
        ea != eb && return ea < eb
        va == v && (i += 1); vb == v && (j += 1)
    end
    false
end

"""Exponent range of every variable over the terms of `p`."""
function ranges(p::PI)
    r = Dict{Int, Tuple{Int, Int}}()
    vs = Set(v for mn in keys(p.t) for (v, _) in mn)
    for v in vs
        es = [mexp(mn, v) for mn in keys(p.t)]
        r[v] = (minimum(es), maximum(es))
    end
    r
end

"""Exact quotient `p / f` of Laurent polynomials, `nothing` when `f` does not divide `p`. Every
quotient term must lie in the exponent box `[min_v p − min_v f, max_v p − max_v f]`, which bounds
the loop."""
function exactdiv(p::PI, f::PI)
    rp, rf = ranges(p), ranges(f)
    lead(q) = foldl((a, b) -> mless(a, b) ? b : a, keys(q.t))
    lf = lead(f); cf = f.t[lf]
    q, r = PI(), p
    while !iszero(r)
        lr = lead(r)
        rem(r.t[lr], cf) == 0 || return nothing
        t = mmul(lr, mpow(lf, -1))
        for v in union(keys(rp), keys(rf))
            lo = get(rp, v, (0, 0))[1] - get(rf, v, (0, 0))[1]
            hi = get(rp, v, (0, 0))[2] - get(rf, v, (0, 0))[2]
            lo <= mexp(t, v) <= hi || return nothing
        end
        c = div(r.t[lr], cf)
        q = q + term(PI, t, c)
        r = r - Poly.scale(f, c) * t
    end
    q
end

"""`p` with each concentration's lowest exponent brought to zero, and the monomial taken out."""
function lift(p::PI, concs)
    s = lowest_shift((p,), concs)
    p * s, mpow(s, -1)
end

# ─── Factored laws ────────────────────────────────────────────────────────

"""The factored law of `m`. Plain: `den = Σ_g σ_g Π_b F_b` (segment partition functions and
block tree sums from P1's engine, renamed to the package's names; `blocks = false` keeps each
segment's tree sum `D_g = Π_b F_b` expanded as one factor) or, with `expand_den`, the
package's expanded (merged) denominator; `num = H · Q · Π extracted factors` (`H` the Haldane
factor, `Q` positive, segment and block factors divided out exactly). Allosteric: the MWC law with
the Haldane factor taken out of the numerator sum (its state denominators are always expanded)."""
function factored_law(m::Pop.PMech; expand_den = false, blocks = true)
    V, np = law_vars(m)
    deps = dep_monos(m, V)
    sub(p) = substdeps(p, np, deps)
    m.allosteric && return mwc_law(m, V, np, deps, N -> haldane_split(sub(N), m, V))
    concs = collect(np + length(deps) + 1:length(V.names))
    E = Engine.bare_rate_law(m; check = false)
    ren = Dict(m.pass2_rename)
    vmap = Dict{Int, Mono}()
    for (j, vs) in enumerate(E.stepvars), (i, v) in enumerate(vs)
        s = m.step_names[j][i]
        vmap[v] = mono(V.index[get(ren, s, s)] => 1)
    end
    for x in m.metabolites; vmap[E.vars.index[x]] = mono(V.index[x] => 1); end
    rn(p) = subst(p, v -> vmap[v])
    num, den = rn(E.num), rn(E.den)                  # in lowest terms already
    topoly(m.pkg_num, V) == num && topoly(m.pkg_den, V) == den || error("$(m.id): engine != package")
    st = E.stats
    bf(g) = [rn(st[:factors][b][st[:attach][b][g]]) for b in eachindex(st[:factors])]
    # with `blocks = false` a segment's tree sum D_g is one (expanded) factor
    gfactors(g) = [rn(st[:sigma][g]); blocks ? bf(g) : [foldl(*, bf(g); init = one(PI))]]
    # the engine's lowest-terms monomial: den = shift · Σ_g σ_g Π_b F_b (the factors are unshifted)
    raw = sum(foldl(*, gfactors(g)) for g in eachindex(st[:sigma]))
    lead(q) = foldl((a, b) -> mless(a, b) ? b : a, keys(q.t))
    shift = mmul(lead(den), mpow(lead(raw), -1))
    raw * shift == den || error("$(m.id): factored denominator != engine denominator")

    nodes = Vector{Term}[]
    ids = Dict{Dict{Mono, Int}, Int}()          # factor polynomial (terms) => node
    nodeof!(q) = get!(() -> (push!(nodes, node_of(q)); length(nodes)), ids, q.t)
    cands = PI[]                                # distinct multi-term factors, lifted
    for g in eachindex(st[:sigma]), p in gfactors(g)
        q, _ = lift(p, concs)
        nterms(q) > 1 && !any(==(q), cands) && push!(cands, q)
    end
    dterms = Term[]
    expand_den || for g in eachindex(st[:sigma])
        facs, mn, coef = Pair{Int, Int}[], shift, 1
        for p in gfactors(g)            # a monomial factor folds into the term, others are nodes
            q, c = lift(p, concs)
            mn = mmul(mn, c)
            if nterms(q) == 1
                k, cq = first(q.t); mn = mmul(mn, k); coef *= cq
            else
                push!(facs, nodeof!(q) => 1)
            end
        end
        push!(dterms, Term(coef, mn, sort!(facs)))
    end
    # Numerator: the Haldane factor, then every factor that divides the quotient exactly.
    sp = haldane_split(sub(num), m, V)
    nfacs = Pair{Int, Int}[]
    if sp === nothing
        push!(nodes, node_of(num)); push!(nfacs, length(nodes) => 1); nmono = Mono(); ncoef = 1
    else
        H, Q = sp
        push!(nodes, node_of(H)); push!(nfacs, length(nodes) => 1)
        for f in cands
            while (q = exactdiv(Q, sub(f))) !== nothing
                Q = q; push!(nfacs, nodeof!(f) => 1)
            end
        end
        Q, c = lift(Q, concs)
        any(<(0), values(Q.t)) && error("$(m.id): negative quotient term")
        nmono, ncoef = c, 1
        if nterms(Q) == 1
            k, ncoef = first(Q.t); nmono = mmul(nmono, k)
        else
            push!(nodes, node_of(Q)); push!(nfacs, length(nodes) => 1)
        end
    end
    push!(nodes, [Term(ncoef, nmono, sort!(nfacs))])
    push!(nodes, expand_den ? node_of(den) : sort!(dterms; by = t -> (t.facs, t.mono)))
    Law(V, np, deps, length(m.metabolites), nodes, length(nodes) - 1, length(nodes), V.index[:E_total])
end

# ─── Checks and statistics ────────────────────────────────────────────────

"""Node `k` of `law` expanded to a polynomial (dependents kept as variables)."""
function expand_node(law::Law, k::Int, memo = Dict{Int, PI}())
    haskey(memo, k) && return memo[k]
    acc = PI()
    for t in law.nodes[k]
        x = term(PI, t.mono, t.coef)
        for (f, pw) in t.facs, _ in 1:pw
            x = x * expand_node(law, f, memo)
        end
        acc = acc + x
    end
    memo[k] = acc
end

"""Value of `law` at `vals` (BigFloat; dependents computed from the parameters)."""
function value(law::Law, vals::AbstractVector)
    x = copy(vals)
    for (i, d) in enumerate(law.deps)
        x[law.np + i] = prod((x[v]^e for (v, e) in d); init = one(eltype(x)))
    end
    memo = Dict{Int, eltype(x)}()
    nv(k) = get!(memo, k) do
        sum((t.coef * prod((x[v]^e for (v, e) in t.mono); init = one(eltype(x))) *
             prod((nv(f)^pw for (f, pw) in t.facs); init = one(eltype(x))) for t in law.nodes[k]);
            init = zero(eltype(x)))
    end
    x[law.scale] * nv(law.num) / nv(law.den)
end

"""Positional values of a package sample: parameters, dependents (left 1, computed by `value`),
concentrations."""
function sample_vals(law::Law, s::Pop.PSample; T = BigFloat)
    x = ones(T, length(law.vars.names))
    for (k, v) in s.params; x[law.vars.index[k]] = v; end
    for (k, v) in s.concs; x[law.vars.index[k]] = v; end
    x
end

"""Arithmetic operations of the flat evaluation of a node (multiplications, additions, divisions
by a per-term denominator) for size statistics."""
function node_terms(law::Law)
    sum(length(n) for n in law.nodes)
end

end # module
