# Lightweight symbolic polynomial type for compile-time rate equation derivation.
# All computation happens on Poly values. Conversion to Expr happens once at the end.

const Mono = Vector{Pair{Symbol,Int}}
const Poly = Dict{Mono, Int}

_mono(pairs...) = sort!(Mono(collect(pairs)); by=first)

poly_zero() = Poly()
poly_one() = Poly(_mono() => 1)
poly_const(n::Integer) = n == 0 ? Poly() : Poly(_mono() => Int(n))
poly_sym(s::Symbol) = Poly(_mono(s => 1) => 1)

function poly_add(a::Poly, b::Poly)
    r = copy(a)
    for (k, v) in b
        r[k] = get(r, k, 0) + v
    end
    filter!(p -> p.second != 0, r)
end

function poly_sub(a::Poly, b::Poly)
    r = copy(a)
    for (k, v) in b
        r[k] = get(r, k, 0) - v
    end
    filter!(p -> p.second != 0, r)
end

function poly_neg(a::Poly)
    Poly(k => -v for (k, v) in a)
end

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

poly_iszero(p::Poly) = isempty(p)

# Cofactor determinant expansion for symbolic matrices
function sym_det(M::Matrix{Poly}, n::Int)
    n == 0 && return poly_one()
    n == 1 && return M[1,1]
    result = poly_zero()
    for j in 1:n
        poly_iszero(M[1,j]) && continue
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

# Convert Poly num/den to a Julia Expr for @generated function bodies
function to_rate_expr(num::Poly, den::Poly, param_syms::Set{Symbol}, conc_syms::Set{Symbol})
    num_expr = _poly_to_expr(num, param_syms, conc_syms)
    den_expr = _poly_to_expr(den, param_syms, conc_syms)
    et = :(params.E_total)
    :($et * ($num_expr) / ($den_expr))
end

function _poly_to_expr(p::Poly, param_syms::Set{Symbol}, conc_syms::Set{Symbol})
    poly_iszero(p) && return 0
    # Separate positive and negative terms for proper a - b expressions
    pos_terms = Any[]
    neg_terms = Any[]
    for (mono, coeff) in sort(collect(p); by=x -> (_mono_sort_key(x[1], param_syms), x[2] < 0))
        factors = Any[]
        abs(coeff) != 1 && push!(factors, abs(coeff))
        for (s, e) in sort(mono; by=sp -> (sp.first in param_syms ? 0 : 1, string(sp.first)))
            accessor = s in param_syms ? :(params.$s) : :(concs.$s)
            if e == 1
                push!(factors, accessor)
            else
                push!(factors, :($accessor ^ $e))
            end
        end
        term = isempty(factors) ? abs(coeff) :
               length(factors) == 1 ? factors[1] :
               _nest_binary(:*, factors)
        coeff > 0 ? push!(pos_terms, term) : push!(neg_terms, term)
    end
    # Build expression: pos_sum - neg_sum using nested binary ops
    pos_expr = isempty(pos_terms) ? nothing :
               length(pos_terms) == 1 ? pos_terms[1] :
               _nest_binary(:+, pos_terms)
    neg_expr = isempty(neg_terms) ? nothing :
               length(neg_terms) == 1 ? neg_terms[1] :
               _nest_binary(:+, neg_terms)
    if pos_expr !== nothing && neg_expr !== nothing
        return :($pos_expr - $neg_expr)
    elseif pos_expr !== nothing
        return pos_expr
    elseif neg_expr !== nothing
        return :(- $neg_expr)
    else
        return 0
    end
end

"""Build nested binary operator calls: +(a, +(b, +(c, d))) to avoid N-ary allocation."""
function _nest_binary(op::Symbol, terms::Vector{Any})
    length(terms) == 1 && return terms[1]
    result = terms[end]
    for i in (length(terms)-1):-1:1
        result = Expr(:call, op, terms[i], result)
    end
    result
end

function _mono_sort_key(mono::Mono, param_syms::Set{Symbol})
    # Sort: constant terms first, then by number of concentration factors
    n_conc = sum(e for (s, e) in mono if s ∉ param_syms; init=0)
    n_conc
end

# Pretty-print a polynomial (bare symbols, no params./concs. prefix)
function poly_str(p::Poly, inverted_params::Set{Symbol}=Set{Symbol}())
    isempty(p) && return "0"
    ts = sort(collect(p); by=x -> (x[2] < 0, _str_sort_key(x[1])))
    parts = String[]
    for (i, (mono, coeff)) in enumerate(ts)
        num_syms = String[]
        den_syms = String[]
        for (s, e) in sort(mono; by=sp -> (_sym_sort_priority(sp.first), string(sp.first)))
            if s in inverted_params
                # Inverted: K appears as 1/K, so positive exponent -> denominator
                if e > 0
                    e == 1 ? push!(den_syms, string(s)) : push!(den_syms, "$(s)^$e")
                elseif e < 0
                    (-e) == 1 ? push!(num_syms, string(s)) : push!(num_syms, "$(s)^$(-e)")
                end
            else
                if e == 1; push!(num_syms, string(s))
                else push!(num_syms, "$(s)^$e"); end
            end
        end
        abs_c = abs(coeff)
        # Build numerator part
        num_parts = String[]
        abs_c != 1 && push!(num_parts, "$abs_c")
        append!(num_parts, num_syms)
        num_str = isempty(num_parts) ? "$abs_c" : join(num_parts, " * ")
        # Build full monomial string
        m = if isempty(den_syms)
            num_str
        else
            den_str = length(den_syms) == 1 ? den_syms[1] : "($(join(den_syms, " * ")))"
            "$num_str / $den_str"
        end
        if i == 1
            push!(parts, coeff < 0 ? "-$m" : m)
        else
            push!(parts, coeff < 0 ? " - $m" : " + $m")
        end
    end
    join(parts)
end

function _sym_sort_priority(s::Symbol)
    is_k_parameter(s) ? 0 : 1
end

function _str_sort_key(mono::Mono)
    syms = sort([string(s) for (s, _) in mono])
    syms
end

"""
    is_k_parameter(sym::Symbol) -> Bool

Check if a symbol is a rate/equilibrium parameter (k or K), not Keq or E_total.
"""
function is_k_parameter(sym::Symbol)
    s = string(sym)
    sym == :Keq && return false
    sym == :E_total && return false
    startswith(s, "k") && return true
    startswith(s, "K") && length(s) > 1 && isdigit(s[2]) && return true
    return false
end

# Build power expression for constraint substitution: Keq^a * prod(k_i^exp_i)
function build_power_expr(keq_exp::Rational, factors::Vector{Tuple{Symbol, Rational{BigInt}}})
    terms = Expr[]
    keq_exp != 0 && push!(terms, _power_accessor(:Keq, keq_exp))
    for (sym, exp) in factors
        exp == 0 && continue
        push!(terms, _power_accessor(sym, exp))
    end
    isempty(terms) && return :(1)
    length(terms) == 1 && return terms[1]
    Expr(:call, :*, terms...)
end

function _power_accessor(sym::Symbol, exp::Rational)
    base = :(params.$sym)
    exp == 1 && return base
    exp == -1 && return :(1 / $base)
    if denominator(exp) == 1
        e = Int(exp)
        return e > 0 ? :($base ^ $e) : :(1 / $base ^ $(-e))
    end
    :($base ^ $(Float64(exp)))
end

# Substitute dependent params in an Expr tree (for HaldaneWegscheider mode)
function substitute_params_expr(expr, subs::Dict{Symbol, Expr})
    if expr isa Expr && expr.head == :. && length(expr.args) == 2 &&
       expr.args[1] == :params && expr.args[2] isa QuoteNode
        sym = expr.args[2].value
        return get(subs, sym, expr)
    elseif expr isa Expr
        return Expr(expr.head, Any[substitute_params_expr(a, subs) for a in expr.args]...)
    else
        return expr
    end
end
