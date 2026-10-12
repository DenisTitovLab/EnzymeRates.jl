# ABOUTME: Rename-proof key of a rate law (module Canon, probe P8): lowest-terms normalization,
# ABOUTME: canonical parameter labeling (color refinement plus search) and SHA-256 of the relabeled law.
module Canon

using SHA, Random

"""A monomial: `name => exponent` pairs sorted by name, no zero exponent."""
const Mono = Vector{Pair{Symbol, Int}}
"""A polynomial: monomial => nonzero coefficient."""
const Poly = Dict{Mono, Rational{BigInt}}

"""
A rate law to key. `comps` are its component polynomials in a fixed role order (a plain
law: num, den; an MWC law: each state's num, den and free-enzyme weight, then each
regulatory site's factor per state). `fixed` names (concentrations, `Keq`) are never
renamed; every other name is a parameter. `units[u]` lists the components that share one
normalization factor; its last component fixes the sign. `params` are the declared fitted
names (a name absent from every component still counts) and `header` the law's scalar
structure (multiplicities).
"""
struct Law
    comps::Vector{Poly}
    fixed::Set{Symbol}
    units::Vector{Vector{Int}}
    params::Vector{Symbol}
    header::String
end

"""
The components with each unit brought to lowest terms: multiplied by the monomial that
takes every variable's lowest exponent over the unit's terms to zero, scaled to coprime
integer coefficients, with a positive coefficient sum in the unit's last component. Equal
functions get equal components when no unit has a non-monomial common factor.
"""
function normalize(law::Law)
    out = [Dict{Mono, BigInt}() for _ in law.comps]
    for u in law.units
        nterm = sum(length(law.comps[i]) for i in u)
        lo, seen = Dict{Symbol, Int}(), Dict{Symbol, Int}()
        for i in u, m in keys(law.comps[i]), (s, e) in m
            lo[s] = min(get(lo, s, e), e); seen[s] = get(seen, s, 0) + 1
        end
        shift = Dict(s => -(seen[s] < nterm ? min(e, 0) : e) for (s, e) in lo)
        scale = foldl(lcm, (denominator(c) for i in u for c in values(law.comps[i]));
                      init = big(1))
        g = foldl(gcd, (numerator(c * scale) for i in u for c in values(law.comps[i]));
                  init = big(0))
        sgn = sign(sum(values(law.comps[u[end]]); init = big(0) // 1))
        sgn == 0 && error("normalize: the sign component sums to zero")
        for i in u, (m, c) in law.comps[i]
            d = Dict(m)
            for (s, e) in shift
                d[s] = get(d, s, 0) + e
            end
            out[i][sort!([s => e for (s, e) in d if e != 0]; by = first)] =
                sgn * numerator(c * scale) ÷ g
        end
    end
    out
end

"""The rank of each element of `xs` among its distinct values (1 = smallest)."""
function ranks(xs)
    u = sort!(unique(xs))
    d = Dict(x => i for (i, x) in enumerate(u))
    [d[x] for x in xs]
end

"""Color refinement of parameters and terms over the incidences `inc[t] = [(p, e)]` and
`occ[p] = [(t, e)]` to the coarsest stable ordered partition below the given colors."""
function refine(pcol, tcol, inc, occ)
    pcol, tcol = ranks(pcol), ranks(tcol)
    np, nt = maximum(pcol; init = 0), maximum(tcol; init = 0)
    while true
        tcol = ranks([(tcol[t], sort!([(pcol[p], e) for (p, e) in inc[t]]))
                      for t in eachindex(inc)])
        pcol = ranks([(pcol[p], sort!([(tcol[t], e) for (t, e) in occ[p]]))
                      for p in eachindex(occ)])
        np2, nt2 = maximum(pcol; init = 0), maximum(tcol; init = 0)
        (np2, nt2) == (np, nt) && return pcol, tcol
        np, nt = np2, nt2
    end
end

"""The law relabeled by the discrete coloring `pcol` (parameter `p` becomes `x<pcol[p]>`), as
sorted lines of text."""
encode(pcol, inc, labels) =
    join(sort!([string(labels[t], "|", join(("x$(pcol[p])^$e" for (p, e) in
                                             sort(inc[t]; by = x -> pcol[x[1]])), "*"))
                for t in eachindex(inc)]), "\n")

"""Search state: the smallest leaf encoding, its parameter order, leaves visited, leaves that
reached the smallest encoding, and the leaf cap."""
mutable struct Search
    best::Union{Nothing, String}
    order::Vector{Int}
    leaves::Int
    hits::Int
    cap::Int
end

"""Individualization-refinement: refine; at a discrete coloring encode a leaf; otherwise
individualize each member of the first non-singleton cell in turn. The smallest leaf
encoding over the whole tree is the canonical form."""
function search!(S::Search, pcol, tcol, inc, occ, labels)
    pcol, tcol = refine(pcol, tcol, inc, occ)
    if allunique(pcol)
        S.leaves += 1
        S.leaves > S.cap && error("canonical labeling: more than $(S.cap) leaves")
        enc = encode(pcol, inc, labels)
        if S.best === nothing || enc < S.best
            S.best, S.order, S.hits = enc, sortperm(pcol), 1
        elseif enc == S.best
            S.hits += 1
        end
        return
    end
    c = minimum(x for x in pcol if count(==(x), pcol) > 1)
    for p in findall(==(c), pcol)
        q = 2 .* pcol
        q[p] -= 1
        search!(S, q, tcol, inc, occ, labels)
    end
end

"""A law's canonical form: its parameters in canonical order (absent declared ones last, by
name), the canonical text, its SHA-256 (the key), the leaves searched and how many reached
the canonical text (more than one: the law has a parameter symmetry)."""
struct Canonical
    order::Vector{Symbol}
    encoding::String
    key::String
    leaves::Int
    hits::Int
end

"""The canonical form of `law` (`cap` bounds the search leaves)."""
function canonical(law::Law; cap = 100_000)
    comps = normalize(law)
    present = sort!(unique(s for c in comps for m in keys(c) for (s, _) in m
                           if !(s in law.fixed)))
    absent = sort!(setdiff(law.params, present))
    pidx = Dict(p => i for (i, p) in enumerate(present))
    labels, inc = String[], Vector{Tuple{Int, Int}}[]
    for (ci, c) in enumerate(comps), (m, coef) in c
        push!(labels, string(ci, "|", coef, "|",
                             join(("$s^$e" for (s, e) in m if s in law.fixed), "*")))
        push!(inc, [(pidx[s], e) for (s, e) in m if !(s in law.fixed)])
    end
    occ = [Tuple{Int, Int}[] for _ in present]
    for (t, ps) in enumerate(inc), (p, e) in ps
        push!(occ[p], (t, e))
    end
    S = Search(nothing, Int[], 0, 0, cap)
    search!(S, ones(Int, length(present)), ranks(labels), inc, occ, labels)
    enc = string(law.header, "|P=", length(law.params), "|absent=", length(absent), "\n",
                 something(S.best, ""))
    Canonical([present[S.order]; absent], enc, bytes2hex(sha256(enc)), S.leaves, S.hits)
end

"""The name map that copies a fit from the representative `r` to the member `m` of one key:
each of `r`'s parameters => the member's parameter in the same canonical position."""
name_map(r::Canonical, m::Canonical) = Dict(zip(r.order, m.order))

"""Whether `law` with its parameters renamed by `map` equals `other`, component by component,
after both are normalized (an exact polynomial identity)."""
function same_under(map, law::Law, other::Law)
    law.header == other.header && length(law.params) == length(other.params) || return false
    a, b = normalize(law), normalize(other)
    ren(m) = sort!([get(map, s, s) => e for (s, e) in m]; by = first)
    length(a) == length(b) && all(Dict(ren(m) => c for (m, c) in a[i]) == b[i] for i in eachindex(a))
end

"""`law` with its parameters renamed by `map` (names missing from `map` kept)."""
rename(law::Law, map) =
    Law([Poly(sort!([get(map, s, s) => e for (s, e) in m]; by = first) => c for (m, c) in p)
         for p in law.comps], law.fixed, law.units, [get(map, s, s) for s in law.params],
        law.header)

# ─── Diagnostic: non-monomial common factors ───────────────────────────

const PRIME = 2_147_483_647

mulmod(a, b) = Int(mod(Int128(a) * Int128(b), PRIME))

"""Value of a monomial modulo `PRIME` at `vals`, leaving out the variable `x`: returns the
coefficient factor and `x`'s exponent."""
function mono_mod(m, vals, x)
    f, ex = 1, 0
    for (s, e) in m
        s === x && (ex = e; continue)
        b = e >= 0 ? vals[s] : invmod(vals[s], PRIME)
        f = mulmod(f, powermod(b, abs(e), PRIME))
    end
    f, ex
end

"""Coefficients (low degree first) modulo `PRIME` of a polynomial as a univariate Laurent
polynomial in `x`, shifted to start at degree 0, at the values `vals` of the other names."""
function univariate(p, x, vals)
    t = Dict{Int, Int}()
    for (m, c) in p
        f, ex = mono_mod(m, vals, x)
        t[ex] = mod(get(t, ex, 0) + mulmod(Int(mod(c, PRIME)), f), PRIME)
    end
    filter!(kv -> kv.second != 0, t)
    isempty(t) && return Int[]
    lo, hi = extrema(keys(t))
    [get(t, d, 0) for d in lo:hi]
end

trim!(a) = (while !isempty(a) && a[end] == 0; pop!(a); end; a)

"""Greatest common divisor of two univariate polynomials modulo `PRIME` (monic or empty)."""
function gcd_mod(a, b)
    a, b = trim!(copy(a)), trim!(copy(b))
    while !isempty(b)
        while length(a) >= length(b)
            f = mulmod(a[end], invmod(b[end], PRIME))
            k = length(a) - length(b)
            for i in eachindex(b)
                a[i + k] = mod(a[i + k] - mulmod(f, b[i]), PRIME)
            end
            trim!(a)
            isempty(a) && break
        end
        a, b = b, a
    end
    isempty(a) ? a : (inv = invmod(a[end], PRIME); [mulmod(c, inv) for c in a])
end

"""Variables `x` for which the components of some unit of `law` keep a common factor of
positive degree in `x` after every other name takes a random value modulo `PRIME`: a
non-monomial common factor (a false alarm has probability about degree/PRIME)."""
function common_factor_vars(law::Law, rng::AbstractRNG)
    comps = normalize(law)
    out = Symbol[]
    for u in law.units
        length(u) < 2 && continue
        ps = [comps[i] for i in u if !isempty(comps[i])]
        length(ps) < 2 && continue
        vars = sort!(unique(s for p in ps for m in keys(p) for (s, _) in m))
        for x in vars
            vals = Dict(s => rand(rng, 2:PRIME - 1) for s in vars if s !== x)
            g = foldl(gcd_mod, (univariate(p, x, vals) for p in ps))
            length(g) > 1 && push!(out, x)
        end
    end
    unique!(out)
end

end # module
