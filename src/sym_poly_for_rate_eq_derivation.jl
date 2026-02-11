# Lightweight symbolic polynomial type for compile-time rate equation derivation.
# All computation happens on Poly values. Conversion to Expr happens once at the end.

const Mono = Vector{Pair{Symbol,Int}}
const Poly = Dict{Mono, Int}

_mono(pairs...) = sort!(Mono(collect(pairs)); by=first)

poly_zero() = Poly()
poly_one() = Poly(_mono() => 1)
poly_const(n::Integer) = n == 0 ? Poly() : Poly(_mono() => Int(n))
poly_sym(s::Symbol) = Poly(_mono(s => 1) => 1)

function _poly_addop(a::Poly, b::Poly, sign::Int)
    r = copy(a)
    for (k, v) in b; r[k] = get(r, k, 0) + sign * v; end
    filter!(p -> p.second != 0, r)
end
poly_add(a::Poly, b::Poly) = _poly_addop(a, b, 1)
poly_sub(a::Poly, b::Poly) = _poly_addop(a, b, -1)
poly_neg(a::Poly) = Poly(k => -v for (k, v) in a)

function poly_mul(a::Poly, b::Poly)
    r = Poly()
    for (k1, v1) in a, (k2, v2) in b
        k = _mono_mul(k1, k2)
        r[k] = get(r, k, 0) + v1 * v2
    end
    filter!(p -> p.second != 0, r)
end

function _mono_mul(a::Mono, b::Mono)
    d = Dict{Symbol,Int}()
    for (s, e) in a; d[s] = get(d, s, 0) + e; end
    for (s, e) in b; d[s] = get(d, s, 0) + e; end
    filter!(p -> p.second != 0, d)
    sort!(Mono(collect(d)); by=first)
end

# Cofactor determinant expansion for symbolic matrices
function sym_det(M::Matrix{Poly}, n::Int)
    n == 0 && return poly_one()
    n == 1 && return M[1,1]
    result = poly_zero()
    for j in 1:n
        isempty(M[1,j]) && continue
        minor = Matrix{Poly}(undef, n-1, n-1)
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

# Convert Poly num/den to a Julia Expr for @generated function bodies (bare symbols)
function to_rate_expr(num::Poly, den::Poly, param_syms::Set{Symbol}, conc_syms::Set{Symbol},
                      inverted_params::Set{Symbol}=Set{Symbol}())
    num_expr = _poly_to_expr(num, param_syms, conc_syms, inverted_params)
    den_expr = _poly_to_expr(den, param_syms, conc_syms, inverted_params)
    _binarize(:(E_total * ($num_expr) / ($den_expr)))
end

"""Convert n-ary +/* calls to left-folded binary for efficient codegen (avoids vararg dispatch)."""
_binarize(x) = x
function _binarize(ex::Expr)
    args = Any[_binarize(a) for a in ex.args]
    if ex.head == :call && length(args) > 3 && args[1] in (:+, :*)
        foldl((a, b) -> Expr(:call, args[1], a, b), args[2:end])
    else
        Expr(ex.head, args...)
    end
end

function _poly_to_expr(p::Poly, param_syms::Set{Symbol}, conc_syms::Set{Symbol},
                       inverted_params::Set{Symbol}=Set{Symbol}())
    isempty(p) && return 0
    pos, neg = Any[], Any[]
    for (mono, coeff) in sort(collect(p); by=x -> (sum(e for (s,e) in x[1] if s ∉ param_syms; init=0), x[2] < 0))
        nf, df = Any[], Any[]
        abs(coeff) != 1 && push!(nf, abs(coeff))
        for (s, e) in sort(mono; by=sp -> (sp.first in param_syms ? 0 : 1, string(sp.first)))
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
    pe !== nothing && ne !== nothing ? :($pe - $ne) :
    pe !== nothing ? pe : ne !== nothing ? :(- $ne) : 0
end

"""Build n-ary operator call: +(a, b, c, d) — prints as `a + b + c + d`."""
_nest_binary(op::Symbol, terms::Vector{Any}) = length(terms) == 1 ? terms[1] : Expr(:call, op, terms...)

"""Check if a symbol is a rate/equilibrium parameter (k or K), not Keq or E_total."""
function is_k_parameter(sym::Symbol)
    sym in (:Keq, :E_total) && return false
    s = string(sym)
    startswith(s, "k") || (startswith(s, "K") && length(s) > 1 && isdigit(s[2]))
end

# Build power expression for constraint substitution: Keq^a * prod(k_i^exp_i)
function build_power_expr(keq_exp::Rational, factors)
    _pa(sym, exp) = exp == 1 ? sym : exp == -1 ? :(1 / $sym) :
        denominator(exp) == 1 ? (Int(exp) > 0 ? :($sym ^ $(Int(exp))) : :(1 / $sym ^ $(Int(-exp)))) :
        :($sym ^ $(Float64(exp)))
    terms = Any[]
    keq_exp != 0 && push!(terms, _pa(:Keq, keq_exp))
    for (sym, exp) in factors; exp != 0 && push!(terms, _pa(sym, exp)); end
    isempty(terms) ? :(1) : length(terms) == 1 ? terms[1] : Expr(:call, :*, terms...)
end

# Substitute symbols in an Expr tree (bare symbol matching)
substitute_params_expr(expr, subs::Dict{Symbol, Expr}) =
    expr isa Symbol ? get(subs, expr, expr) :
    expr isa Expr ? Expr(expr.head, Any[substitute_params_expr(a, subs) for a in expr.args]...) : expr

# ─── Parameter constraint substitution in Poly ───────────────────────────────

"""
Substitute `target` symbol in polynomial `p` with `coeff * prod(sym^exp for (sym,exp) in replacement)`.
"""
function _substitute_sym_in_poly(p::Poly, target::Symbol, coeff::Int, replacement)
    result = Poly()
    for (mono, val) in p
        idx = findfirst(pair -> pair.first == target, mono)
        if idx === nothing
            result[mono] = get(result, mono, 0) + val
        else
            e = mono[idx].second
            base = Mono([pair for (i, pair) in enumerate(mono) if i != idx])
            repl = sort!(Mono([sym => exp * e for (sym, exp) in replacement]); by=first)
            final = _mono_mul(base, repl)
            result[final] = get(result, final, 0) + val * coeff^e
        end
    end
    filter!(p -> p.second != 0, result)
end

"""Apply all parameter constraints sequentially to a polynomial."""
function _apply_param_constraints(p::Poly, constraints)
    for (target, coeff, factors) in constraints
        p = _substitute_sym_in_poly(p, target, coeff, factors)
    end
    p
end
