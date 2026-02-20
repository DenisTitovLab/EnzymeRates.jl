# Lightweight symbolic polynomial type for compile-time rate equation derivation.
# All computation happens on POLY values. Conversion to Expr happens once at the end.

const MONO = Vector{Pair{Symbol,Int}}
const POLY = Dict{MONO, Int}

_mono(pairs...) = sort!(MONO(collect(pairs)); by=first)

poly_zero() = POLY()
poly_one() = POLY(_mono() => 1)
poly_const(n::Integer) = n == 0 ? POLY() : POLY(_mono() => Int(n))
poly_sym(s::Symbol) = POLY(_mono(s => 1) => 1)

function _poly_addop(a::POLY, b::POLY, sign::Int)
    r = copy(a)
    for (k, v) in b; r[k] = get(r, k, 0) + sign * v; end
    filter!(p -> p.second != 0, r)
end
poly_add(a::POLY, b::POLY) = _poly_addop(a, b, 1)
poly_sub(a::POLY, b::POLY) = _poly_addop(a, b, -1)
poly_neg(a::POLY) = POLY(k => -v for (k, v) in a)

function poly_mul(a::POLY, b::POLY)
    r = POLY()
    for (k1, v1) in a, (k2, v2) in b
        k = _mono_mul(k1, k2)
        r[k] = get(r, k, 0) + v1 * v2
    end
    filter!(p -> p.second != 0, r)
end

function _mono_mul(a::MONO, b::MONO)
    d = Dict{Symbol,Int}()
    for (s, e) in a; d[s] = get(d, s, 0) + e; end
    for (s, e) in b; d[s] = get(d, s, 0) + e; end
    filter!(p -> p.second != 0, d)
    sort!(MONO(collect(d)); by=first)
end

# Cofactor determinant expansion for symbolic matrices
function sym_det(M::Matrix{POLY}, n::Int)
    n == 0 && return poly_one()
    n == 1 && return M[1,1]
    result = poly_zero()
    for j in 1:n
        isempty(M[1,j]) && continue
        minor = Matrix{POLY}(undef, n-1, n-1)
        for r in 2:n, c in 1:n
            c < j && (minor[r-1, c] = M[r, c])
            c > j && (minor[r-1, c-1] = M[r, c])
        end
        cofactor = sym_det(minor, n-1)
        term = poly_mul(M[1,j], cofactor)
        result = iseven(j-1) ? poly_add(result, term) : poly_sub(result, term)
    end
    result
end

# Convert POLY num/den to a Julia Expr for @generated function bodies (bare symbols)
function to_rate_expr(
    num::POLY, den::POLY,
    param_syms::Set{Symbol}, conc_syms::Set{Symbol},
    inverted_params::Set{Symbol}=Set{Symbol}(),
)
    num_expr = _poly_to_expr(num, param_syms, conc_syms, inverted_params)
    den_expr = _poly_to_expr(den, param_syms, conc_syms, inverted_params)
    :(E_total * ($num_expr) / ($den_expr))
end

function _poly_to_expr(p::POLY, param_syms::Set{Symbol}, conc_syms::Set{Symbol},
                       inverted_params::Set{Symbol}=Set{Symbol}())
    isempty(p) && return 0
    pos, neg = Any[], Any[]
    sorted = sort(
        collect(p);
        by=x -> (
            sum(e for (s,e) in x[1] if s ∉ param_syms; init=0),
            x[2] < 0,
        ),
    )
    for (mono, coeff) in sorted
        nf, df = Any[], Any[]
        abs(coeff) != 1 && push!(nf, abs(coeff))
        sorted_mono = sort(
            mono;
            by=sp -> (sp.first in param_syms ? 0 : 1, string(sp.first)),
        )
        for (s, e) in sorted_mono
            if s in inverted_params
                tgt, ex = e > 0 ? (df, e) : (nf, -e)
            else
                tgt, ex = nf, e
            end
            ex != 0 && push!(tgt, ex == 1 ? s : :($s ^ $ex))
        end
        num_part = isempty(nf) ? abs(coeff) : _nest_binary(:*, nf)
        term = isempty(df) ? num_part : :($num_part / $(_nest_binary(:*, df)))
        coeff > 0 ? push!(pos, term) : push!(neg, term)
    end
    _combine_terms(pos, neg)
end

function _combine_terms(pos::Vector{Any}, neg::Vector{Any})
    pe = isempty(pos) ? nothing : _nest_binary(:+, pos)
    ne = isempty(neg) ? nothing : _nest_binary(:+, neg)
    if pe !== nothing && ne !== nothing
        :($pe - $ne)
    elseif pe !== nothing
        pe
    elseif ne !== nothing
        :(- $ne)
    else
        0
    end
end

"""Build balanced binary tree: +(+(a,b), +(c,d)) — O(log N) depth, zero-alloc runtime."""
function _nest_binary(op::Symbol, terms::Vector{Any})
    n = length(terms)
    if n == 1
        terms[1]
    elseif n == 2
        Expr(:call, op, terms[1], terms[2])
    else
        mid = n >> 1
        Expr(:call, op,
            _nest_binary(op, terms[1:mid]),
            _nest_binary(op, terms[mid+1:end]))
    end
end

"""Check if a symbol is a rate/equilibrium parameter (k or K), not Keq or E_total."""
function is_k_parameter(sym::Symbol)
    sym in (:Keq, :E_total) && return false
    s = string(sym)
    startswith(s, "k") || (startswith(s, "K") && length(s) > 1 && isdigit(s[2]))
end

# Build power expression for constraint substitution: Keq^a * prod(k_i^exp_i)
function build_power_expr(keq_exp::Rational, factors)
    function _pa(sym, exp)
        if exp == 1
            sym
        elseif exp == -1
            :(1 / $sym)
        elseif denominator(exp) == 1
            if Int(exp) > 0
                :($sym ^ $(Int(exp)))
            else
                :(1 / $sym ^ $(Int(-exp)))
            end
        else
            :($sym ^ $(Float64(exp)))
        end
    end
    terms = Any[]
    keq_exp != 0 && push!(terms, _pa(:Keq, keq_exp))
    for (sym, exp) in factors
        exp != 0 && push!(terms, _pa(sym, exp))
    end
    if isempty(terms)
        :(1)
    elseif length(terms) == 1
        terms[1]
    else
        Expr(:call, :*, terms...)
    end
end

# Substitute symbols in an Expr tree (bare symbol matching)
function substitute_params_expr(expr, subs::Dict{Symbol, Expr})
    if expr isa Symbol
        get(subs, expr, expr)
    elseif expr isa Expr
        args = Any[
            substitute_params_expr(a, subs) for a in expr.args
        ]
        Expr(expr.head, args...)
    else
        expr
    end
end

# ─── Parameter constraint substitution in POLY ─────────────

"""
Substitute `target` symbol in polynomial `p` with
`coeff * prod(sym^exp for (sym,exp) in replacement)`.
"""
function _substitute_sym_in_poly(p::POLY, target::Symbol, coeff::Int, replacement)
    result = POLY()
    for (mono, val) in p
        idx = findfirst(pair -> pair.first == target, mono)
        if idx === nothing
            result[mono] = get(result, mono, 0) + val
        else
            e = mono[idx].second
            base = MONO([pair for (i, pair) in enumerate(mono) if i != idx])
            repl = sort!(MONO([sym => exp * e for (sym, exp) in replacement]); by=first)
            final = _mono_mul(base, repl)
            result[final] = get(result, final, 0) + val * coeff^e
        end
    end
    filter!(p -> p.second != 0, result)
end

"""Apply all parameter constraints sequentially to a polynomial."""
function _apply_param_constraints(p::POLY, constraints)
    for (target, coeff, factors) in constraints
        p = _substitute_sym_in_poly(p, target, coeff, factors)
    end
    p
end
