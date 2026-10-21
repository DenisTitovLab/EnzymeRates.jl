# ABOUTME: Exact sparse polynomials over integer-indexed variables: concentrations with
# ABOUTME: integer exponents, constants with rational exponents, checked Int64 coefficients.

"""
A variable table: variable `i` is `names[i]`, a concentration when `conc[i]`; `index` maps
each name to its variable. Two tables are equal when they name the same variables in the
same order.
"""
struct _Vars
    names::Vector{Symbol}
    conc::BitVector
    index::Dict{Symbol, Int}
end
_Vars() = _Vars(Symbol[], falses(0), Dict{Symbol, Int}())
Base.:(==)(a::_Vars, b::_Vars) = a.names == b.names && a.conc == b.conc
Base.hash(t::_Vars, h::UInt) = hash(t.conc, hash(t.names, hash(:_Vars, h)))

"""The index of variable `s` in `t`, added when new; `conc` marks a concentration. Errors
when `s` is already a variable of the other kind."""
function _var!(t::_Vars, s::Symbol; conc::Bool = false)
    i = get(t.index, s, 0)
    i == 0 || return t.conc[i] == conc ? i :
        error("internal error: $s is a concentration in one use and a constant in another")
    push!(t.names, s)
    push!(t.conc, conc)
    t.index[s] = length(t.names)
end

"""An exponent: an integer for a concentration, a rational for a constant."""
const _Exp = Rational{Int}

"""A monomial: `variable => exponent` pairs sorted by variable, no zero exponent. It is
never mutated once built, since polynomials share their monomials as `Dict` keys."""
const _Mono = Vector{Pair{Int, _Exp}}

"""The product of two monomials (a merge of the sorted pairs)."""
function _mmul(a::_Mono, b::_Mono)
    r = _Mono()
    i, j = 1, 1
    while i <= length(a) || j <= length(b)
        if j > length(b) || (i <= length(a) && a[i].first < b[j].first)
            push!(r, a[i])
            i += 1
        elseif i > length(a) || b[j].first < a[i].first
            push!(r, b[j])
            j += 1
        else
            e = a[i].second + b[j].second
            iszero(e) || push!(r, a[i].first => e)
            i += 1
            j += 1
        end
    end
    r
end

"""`a` to the power `e`."""
_mpow(a::_Mono, e::Union{Integer, _Exp}) =
    iszero(e) ? _Mono() : _Mono([v => x * e for (v, x) in a])

"""The monomial of the `variable => exponent` pairs `ps`; a variable may repeat."""
_monomial(ps::Pair{Int}...) =
    foldl(_mmul, (_Mono([v => _Exp(e)]) for (v, e) in ps if e != 0); init = _Mono())

"""The exponent of variable `v` in `a`."""
function _mexp(a::_Mono, v::Int)
    i = findfirst(p -> p.first == v, a)
    i === nothing ? zero(_Exp) : a[i].second
end

"""
A polynomial: monomial => nonzero coefficient. Every sum and product is checked, so a
coefficient that leaves `Int64` throws `OverflowError` instead of wrapping. It is mutable,
with a `const` field, so a loop that merges it with another type into a `Union` holds a
reference, never an inline box the compiler may hoist. Two polynomials are equal when
their terms are.
"""
mutable struct _Poly
    const terms::Dict{_Mono, Int}
end
_Poly() = _Poly(Dict{_Mono, Int}())

"""Add `c·m` into `p`, dropping a zero sum; returns `p`."""
function _addto!(p::_Poly, m::_Mono, c::Integer)
    v = Base.checked_add(get(p.terms, m, 0), Int(c))
    iszero(v) ? delete!(p.terms, m) : (p.terms[m] = v)
    p
end

"""The one-term polynomial `c·m`."""
_term(m::_Mono, c::Integer = 1) = _addto!(_Poly(), m, c)

Base.one(::Type{_Poly}) = _term(_Mono())
Base.zero(::Type{_Poly}) = _Poly()
Base.iszero(p::_Poly) = isempty(p.terms)
Base.:(==)(p::_Poly, q::_Poly) = p.terms == q.terms
Base.hash(p::_Poly, h::UInt) = hash(p.terms, hash(:_Poly, h))

"""The number of terms of `p`."""
_nterms(p::_Poly) = length(p.terms)

function Base.:+(p::_Poly, q::_Poly)
    r = _Poly(copy(p.terms))
    for (m, c) in q.terms
        _addto!(r, m, c)
    end
    r
end
Base.:-(p::_Poly) = _scale(p, -1)
Base.:-(p::_Poly, q::_Poly) = p + (-q)

function Base.:*(p::_Poly, q::_Poly)
    r = _Poly()
    for (m1, c1) in p.terms, (m2, c2) in q.terms
        _addto!(r, _mmul(m1, m2), Base.checked_mul(c1, c2))
    end
    r
end
# Multiplying by a monomial maps distinct monomials to distinct ones, so no terms merge.
Base.:*(p::_Poly, m::_Mono) =
    _Poly(Dict{_Mono, Int}(_mmul(k, m) => c for (k, c) in p.terms))

"""`c·p` for an integer `c`."""
_scale(p::_Poly, c::Integer) = iszero(c) ? _Poly() :
    _Poly(Dict{_Mono, Int}(m => Base.checked_mul(x, Int(c)) for (m, x) in p.terms))

"""
`p` with each variable `v` for which `f(v)` returns a monomial (a power product) replaced
by it, `v^e` by `f(v)^e`; `f(v) === nothing` keeps `v`.
"""
function _subst(f, p::_Poly)
    r = _Poly()
    for (m, c) in p.terms
        out = _Mono()
        for (v, e) in m
            s = f(v)
            out = _mmul(out, s === nothing ? _Mono([v => e]) : _mpow(s, e))
        end
        _addto!(r, out, c)
    end
    r
end

"""
    _divexact(p, d) -> p / d

The exact quotient of `p` by a two-term `d` (the Haldane factor); an `ArgumentError` when
`d` does not divide `p`. With `d = c_u·u + c_w·w` and `t = w/u`, the terms of `p/u` fall
into chains `key·t^k` over integer `k`; dividing a chain by `c_u + c_w·t` is the recurrence
`b_k = (a_k − c_w·b_{k−1}) / c_u`, which must end at zero.
"""
function _divexact(p::_Poly, d::_Poly)
    _nterms(d) == 2 || throw(ArgumentError("_divexact divides by a two-term polynomial"))
    (u, cu), (w, cw) = collect(d.terms)
    t = _mmul(w, _mpow(u, -1))
    v, et = first(t)
    chains = Dict{_Mono, Vector{Pair{Int, Int}}}()
    for (m, c) in p.terms
        m = _mmul(m, _mpow(u, -1))
        k = floor(Int, _mexp(m, v) / et)
        push!(get!(chains, _mmul(m, _mpow(t, -k)), Pair{Int, Int}[]), k => c)
    end
    q = _Poly()
    for (key, chain) in chains
        sort!(chain; by = first)
        b, j = 0, 1
        for k in chain[1].first:chain[end].first
            a = 0
            if chain[j].first == k
                a = chain[j].second
                j += 1
            end
            r = Base.checked_sub(a, Base.checked_mul(cw, b))
            rem(r, cu) == 0 ||
                throw(ArgumentError("_divexact: the quotient is not integral"))
            b = div(r, cu)
            k == chain[end].first || iszero(b) || _addto!(q, _mmul(key, _mpow(t, k)), b)
        end
        iszero(b) || throw(ArgumentError("_divexact: the divisor does not divide"))
    end
    q
end

"""
    _lowest_shift(ps, counted = every variable) -> _Mono

The monomial in the `counted` variables that, multiplied into every polynomial of `ps`,
brings each such variable's lowest exponent over all their terms to zero, a term without
the variable counting as exponent 0: the shift to lowest terms.
"""
function _lowest_shift(ps, counted = _ -> true)
    low, seen, n = Dict{Int, _Exp}(), Dict{Int, Int}(), 0
    for p in ps, m in keys(p.terms)
        n += 1
        for (v, e) in m
            counted(v) || continue
            low[v] = min(get(low, v, e), e)
            seen[v] = get(seen, v, 0) + 1
        end
    end
    _monomial((v => -(seen[v] < n ? min(e, 0) : e) for (v, e) in low)...)
end

"""The shift to lowest terms in the concentrations of `t` alone (`_lowest_shift`)."""
_conc_shift(ps, t::_Vars) = _lowest_shift(ps, v -> t.conc[v])

"""`p` as `concentration monomial => its coefficient`, a polynomial in the constants."""
function _split_conc(p::_Poly, t::_Vars)
    out = Dict{_Mono, _Poly}()
    for (m, c) in p.terms
        _addto!(get!(_Poly, out, filter(q -> t.conc[q.first], m)),
                filter(q -> !t.conc[q.first], m), c)
    end
    out
end

"""
The nonzero polynomial `p` as `c·mu·rest`: its integer content `c` (the coefficients'
gcd, positive, or a one-term polynomial's own coefficient), its monomial content `mu`
(each variable's lowest exponent over the terms, a term without the variable counting as
exponent 0) and the rest, `one(_Poly)` for a one-term polynomial.
"""
function _content(p::_Poly)
    c = foldl(gcd, values(p.terms))
    mu = _monomial((v => minimum(mo -> _mexp(mo, v), keys(p.terms))
                    for v in unique(first(t) for mo in keys(p.terms) for t in mo))...)
    inv_mu = _mpow(mu, -1)
    c, mu, _Poly(Dict{_Mono, Int}(_mmul(mo, inv_mu) => div(x, c) for (mo, x) in p.terms))
end
