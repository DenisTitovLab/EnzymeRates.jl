# ABOUTME: Variant F (synth/mkF.py) of the slice file below. Variant E (verifier): the slice's polynomials with Poly a mutable struct (const field); concentrations with integer
# ABOUTME: exponents, constants with rational exponents, checked Int or BigInt coefficients.

"""A variable table: variable `i` is `names[i]`, a concentration when `conc[i]`."""
struct Vars
    names::Vector{Symbol}
    conc::BitVector
    index::Dict{Symbol, Int}
end
Vars() = Vars(Symbol[], falses(0), Dict{Symbol, Int}())

"""The index of variable `s` in `t`, added when new; `conc` marks a concentration."""
function var!(t::Vars, s::Symbol; conc::Bool = false)
    i = get(t.index, s, 0)
    i == 0 || return t.conc[i] == conc ? i : error("$s is a concentration in one use only")
    push!(t.names, s)
    push!(t.conc, conc)
    t.index[s] = length(t.names)
end

"""An exponent: an integer for a concentration, a rational for a constant."""
const Exp = Rational{Int}

"""A monomial: `variable => exponent` pairs sorted by variable, no zero exponent."""
const Mono = Vector{Pair{Int, Exp}}

"""The product of two monomials (a merge of the sorted pairs)."""
function mmul(a::Mono, b::Mono)
    r = Mono()
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
mpow(a::Mono, e) = iszero(e) ? Mono() : Mono([v => x * e for (v, x) in a])

"""The monomial of `variable => exponent` pairs (a variable may repeat)."""
mono(ps::Pair...) = foldl(mmul, (Mono([v => Exp(e)]) for (v, e) in ps if e != 0);
                          init = Mono())

"""The exponent of variable `v` in `a`."""
function mexp(a::Mono, v)
    i = findfirst(p -> p.first == v, a)
    i === nothing ? zero(Exp) : a[i].second
end

"""A polynomial: monomial => nonzero coefficient, an overflow-checked `Int` or a
`BigInt`."""
mutable struct Poly{T <: Integer}
    const terms::Dict{Mono, T}
end
Poly{T}() where {T} = Poly{T}(Dict{Mono, T}())
"""The one-term polynomial `c·m`."""
term(::Type{T}, m::Mono, c = one(T)) where {T} = addto!(Poly{T}(), m, c)
Base.one(::Type{Poly{T}}) where {T} = term(T, Mono())
Base.zero(::Type{Poly{T}}) where {T} = Poly{T}()
Base.:(==)(p::Poly, q::Poly) = p.terms == q.terms
Base.hash(p::Poly, h::UInt) = hash(p.terms, hash(:Poly, h))
Base.iszero(p::Poly) = isempty(p.terms)
nterms(p::Poly) = length(p.terms)
Base.convert(::Type{Poly{BigInt}}, p::Poly) =
    Poly{BigInt}(Dict{Mono, BigInt}(m => c for (m, c) in p.terms))

add(a::Int, b::Int) = Base.checked_add(a, b)
add(a, b) = a + b
mul(a::Int, b::Int) = Base.checked_mul(a, b)
mul(a, b) = a * b

"""Add `c·m` into `p`, dropping a zero sum; returns `p`."""
function addto!(p::Poly{T}, m::Mono, c) where {T}
    v = add(get(p.terms, m, zero(T)), T(c))
    iszero(v) ? delete!(p.terms, m) : (p.terms[m] = v)
    p
end

function Base.:+(p::Poly{T}, q::Poly{T}) where {T}
    r = Poly{T}(copy(p.terms))
    foreach(((m, c),) -> addto!(r, m, c), q.terms)
    r
end
Base.:-(p::Poly{T}) where {T} = scale(p, -one(T))
Base.:-(p::Poly{T}, q::Poly{T}) where {T} = p + (-q)

function Base.:*(p::Poly{T}, q::Poly{T}) where {T}
    r = Poly{T}()
    for (m1, c1) in p.terms, (m2, c2) in q.terms
        addto!(r, mmul(m1, m2), mul(c1, c2))
    end
    r
end
Base.:*(p::Poly{T}, m::Mono) where {T} =
    Poly{T}(Dict{Mono, T}(mmul(k, m) => c for (k, c) in p.terms))

"""`c·p` for an integer `c`."""
scale(p::Poly{T}, c) where {T} =
    iszero(c) ? Poly{T}() : Poly{T}(Dict{Mono, T}(m => mul(x, T(c)) for (m, x) in p.terms))

"""`p` with each variable `v` for which `f(v)` returns a monomial (a power product) replaced
by it, `v^e` by `f(v)^e`; `f(v) === nothing` keeps `v`."""
function subst(p::Poly{T}, f) where {T}
    r = Poly{T}()
    for (m, c) in p.terms
        out = Mono()
        for (v, e) in m
            s = f(v)
            out = mmul(out, s === nothing ? Mono([v => e]) : mpow(s, e))
        end
        addto!(r, out, c)
    end
    r
end

"""`x^e` for a rational `e`: integer powers exactly, half powers through `sqrt`."""
function pow(x::T, e::Exp) where {T}
    n, d = numerator(e), denominator(e)
    d == 1 ? x^n : d == 2 ? sqrt(x)^n : x^(T(n) / d)
end

"""The value of `p` at `vals` (indexed by variable) in the element type of `vals`."""
evaluate(p::Poly, vals::AbstractVector{T}) where {T} =
    sum((T(c) * prod((pow(vals[v], e) for (v, e) in m); init = one(T))
         for (m, c) in p.terms); init = zero(T))

"""
    divexact(p, d) -> p / d

The exact quotient of `p` by a two-term `d` (the Haldane factor), or an error when `d`
does not divide `p`. With `d = c_u·u + c_w·w` and `t = w/u`, the terms of `p/u` fall into
chains `key·t^k` over integer `k`; dividing a chain by `c_u + c_w·t` is the recurrence
`b_k = (a_k - c_w·b_{k-1}) / c_u`, which must end at zero.
"""
function divexact(p::Poly{T}, d::Poly{T}) where {T}
    nterms(d) == 2 || throw(ArgumentError("divexact divides by a two-term polynomial"))
    (u, cu), (w, cw) = collect(d.terms)
    t = mmul(w, mpow(u, -1))
    v, et = first(t)
    chains = Dict{Mono, Vector{Pair{Int, T}}}()
    for (m, c) in p.terms
        m = mmul(m, mpow(u, -1))
        k = floor(Int, mexp(m, v) / et)
        push!(get!(chains, mmul(m, mpow(t, -k)), Pair{Int, T}[]), k => c)
    end
    q = Poly{T}()
    for (key, chain) in chains
        sort!(chain; by = first)
        b, j = zero(T), 1
        for k in chain[1].first:chain[end].first
            a = zero(T)
            if chain[j].first == k
                a = chain[j].second
                j += 1
            end
            r = add(a, -mul(cw, b))
            rem(r, cu) == 0 ||
                throw(ArgumentError("divexact: the quotient is not integral"))
            b = div(r, cu)
            k == chain[end].first || iszero(b) || addto!(q, mmul(key, mpow(t, k)), b)
        end
        iszero(b) || throw(ArgumentError("divexact: the divisor does not divide"))
    end
    q
end

"""
    lowest_shift(ps, counted = every variable) -> Mono

The monomial in the `counted` variables that, multiplied into every polynomial of `ps`,
brings each such variable's lowest exponent over all their terms to zero: the shift to
lowest terms.
"""
function lowest_shift(ps, counted = v -> true)
    low, seen, n = Dict{Int, Exp}(), Dict{Int, Int}(), 0
    for p in ps, m in keys(p.terms)
        n += 1
        for (v, e) in m
            counted(v) || continue
            low[v] = min(get(low, v, e), e)
            seen[v] = get(seen, v, 0) + 1
        end
    end
    mono((v => -(seen[v] < n ? min(e, 0) : e) for (v, e) in low)...)
end

"""The shift to lowest terms in the concentrations of `t` alone (`lowest_shift`)."""
conc_shift(ps, t::Vars) = lowest_shift(ps, v -> t.conc[v])

"""`p` as `concentration monomial => its coefficient`, a polynomial in the constants."""
function split_conc(p::Poly{T}, t::Vars) where {T}
    out = Dict{Mono, Poly{T}}()
    for (m, c) in p.terms
        addto!(get!(Poly{T}, out, filter(q -> t.conc[q.first], m)),
               filter(q -> !t.conc[q.first], m), c)
    end
    out
end
