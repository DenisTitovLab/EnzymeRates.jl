# ABOUTME: Sparse Laurent polynomials over integer-indexed variables with a symbol table (module Poly).
# ABOUTME: Monomials are sorted var => exponent vectors; Int or Rational{Int} coefficients, overflow-checked.
module Poly

"""Symbol table of a polynomial ring: variable `i` is `names[i]`."""
struct Vars
    names::Vector{Symbol}
    index::Dict{Symbol, Int}
end
Vars() = Vars(Symbol[], Dict{Symbol, Int}())

"""Index of variable `s` in `t`, added when new."""
var!(t::Vars, s::Symbol) = get!(() -> (push!(t.names, s); length(t.names)), t.index, s)

"""A Laurent monomial: `variable => exponent` pairs sorted by variable, no zero exponent."""
const Mono = Vector{Pair{Int, Int}}

"""The monomial of the given pairs (distinct variables)."""
mono(ps...) = sort!(Mono([p for p in ps if p.second != 0]); by = first)

"""Product of two monomials (a merge of the sorted pairs)."""
function mmul(a::Mono, b::Mono)
    r = Mono()
    i, j = 1, 1
    while i <= length(a) || j <= length(b)
        if j > length(b) || (i <= length(a) && a[i].first < b[j].first)
            push!(r, a[i]); i += 1
        elseif i > length(a) || b[j].first < a[i].first
            push!(r, b[j]); j += 1
        else
            e = a[i].second + b[j].second
            e == 0 || push!(r, a[i].first => e)
            i += 1; j += 1
        end
    end
    r
end

"""`a` to the integer power `k` (negative allowed)."""
mpow(a::Mono, k::Integer) = k == 0 ? Mono() : Mono([v => e * k for (v, e) in a])

"""Exponent of variable `v` in `a`."""
mexp(a::Mono, v::Int) = (i = findfirst(p -> p.first == v, a); i === nothing ? 0 : a[i].second)

"""A polynomial: monomial => nonzero coefficient."""
struct P{T <: Union{Int, Rational{Int}}}
    t::Dict{Mono, T}
end
P{T}() where {T} = P{T}(Dict{Mono, T}())
"""The one-term polynomial `c·m`."""
term(::Type{P{T}}, m::Mono, c = one(T)) where {T} = P{T}(Dict{Mono, T}(m => c))
Base.one(::Type{P{T}}) where {T} = term(P{T}, Mono())
Base.zero(::Type{P{T}}) where {T} = P{T}()

chk_add(a::Int, b::Int) = Base.checked_add(a, b)
chk_add(a, b) = a + b                      # Rational{Int} arithmetic is checked
chk_mul(a::Int, b::Int) = Base.checked_mul(a, b)
chk_mul(a, b) = a * b

"""Add `c·m` into the term dictionary `d`, dropping a zero sum."""
function addto!(d::Dict{Mono, T}, m::Mono, c) where {T}
    v = chk_add(get(d, m, zero(T)), T(c))
    iszero(v) ? delete!(d, m) : (d[m] = v)
    d
end

function Base.:+(p::P{T}, q::P{T}) where {T}
    d = copy(p.t)
    for (m, c) in q.t
        addto!(d, m, c)
    end
    P{T}(d)
end
Base.:-(p::P{T}) where {T} = P{T}(Dict{Mono, T}(m => -c for (m, c) in p.t))
Base.:-(p::P{T}, q::P{T}) where {T} = p + (-q)

function Base.:*(p::P{T}, q::P{T}) where {T}
    d = Dict{Mono, T}()
    for (m1, c1) in p.t, (m2, c2) in q.t
        addto!(d, mmul(m1, m2), chk_mul(c1, c2))
    end
    P{T}(d)
end
# Multiplying by a monomial maps distinct monomials to distinct monomials.
Base.:*(p::P{T}, m::Mono) where {T} = P{T}(Dict{Mono, T}(mmul(k, m) => c for (k, c) in p.t))

"""`c·p` for a scalar `c`."""
scale(p::P{T}, c) where {T} =
    iszero(c) ? P{T}() : P{T}(Dict{Mono, T}(m => chk_mul(x, T(c)) for (m, x) in p.t))

Base.:(==)(p::P, q::P) = p.t == q.t
Base.iszero(p::P) = isempty(p.t)
nterms(p::P) = length(p.t)

"""Expansion of a product of factors."""
expand(::Type{P{T}}, factors) where {T} = foldl(*, factors; init = one(P{T}))

"""`p` with each variable `v` for which `f(v)` returns a monomial replaced by that monomial
(the Laurent power follows the exponent); `f(v) === nothing` keeps `v`."""
function subst(p::P{T}, f) where {T}
    d = Dict{Mono, T}()
    for (m, c) in p.t
        r = Mono()
        for (v, e) in m
            s = f(v)
            r = mmul(r, s === nothing ? Mono([v => e]) : mpow(s, e))
        end
        addto!(d, r, c)
    end
    P{T}(d)
end

"""The monomial that brings the lowest exponent of each variable of `vs` over all terms of
the polynomials `ps` to zero (a variable absent from a term has exponent 0 there)."""
function lowest_shift(ps, vs)
    mins = Dict{Int, Int}(v => 0 for v in vs)
    first_term = true
    for p in ps, m in keys(p.t)
        for v in vs
            e = mexp(m, v)
            mins[v] = first_term ? e : min(mins[v], e)
        end
        first_term = false
    end
    mono((v => -e for (v, e) in mins)...)
end

"""`p` with symbol-named monomials sorted by symbol, for comparison with other engines."""
symbolic(p::P{T}, t::Vars) where {T} =
    Dict{Vector{Pair{Symbol, Int}}, Rational{Int}}(
        sort!([t.names[v] => e for (v, e) in m]; by = first) => c for (m, c) in p.t)

"""Value of `p` at `vals` (indexed by variable)."""
evaluate(p::P, vals::AbstractVector) =
    sum((c * prod((vals[v]^e for (v, e) in m); init = one(eltype(vals)))
         for (m, c) in p.t); init = zero(eltype(vals)))

end # module
