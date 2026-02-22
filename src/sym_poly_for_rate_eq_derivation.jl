# Lightweight symbolic polynomial type for compile-time rate equation derivation.
# All computation happens on POLY values. Conversion to Expr happens once at the end.

# Maximum raw polynomial terms allowed in a rate equation.
# Equations exceeding this limit would take too long to compile
# via @generated functions and are unlikely to be useful.
const MAX_RATE_EQUATION_TERMS = 5000

const MONO = Vector{Pair{Symbol,Int}}
const POLY = Dict{MONO, Rational{Int}}

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

function _mono_div(a::MONO, b::MONO)
    d = Dict{Symbol,Int}()
    for (s, e) in a; d[s] = get(d, s, 0) + e; end
    for (s, e) in b; d[s] = get(d, s, 0) - e; end
    filter!(p -> p.second != 0, d)
    sort!(MONO(collect(d)); by=first)
end

"""Divide POLY by a single-term POLY (exact division, assumes divisibility)."""
function _poly_div_mono(p::POLY, divisor::POLY)::POLY
    m = first(keys(divisor))
    POLY(_mono_div(k, m) => v for (k, v) in p)
end

# Cofactor determinant expansion for symbolic matrices.
# Checks intermediate term count against MAX_RATE_EQUATION_TERMS to
# abort early for mechanisms whose rate equations would be too large.
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
        if length(result) > MAX_RATE_EQUATION_TERMS
            error(
                "Rate equation for this mechanism has more than " *
                "$MAX_RATE_EQUATION_TERMS polynomial terms " *
                "(limit: $MAX_RATE_EQUATION_TERMS). Equations " *
                "this large take a very long time to compile " *
                "and are unlikely to be practically useful " *
                "for parameter fitting.",
            )
        end
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
        abs_c = abs(coeff)
        cn = Int(numerator(abs_c))
        cd = Int(denominator(abs_c))
        cn != 1 && push!(nf, cn)
        cd != 1 && push!(df, cd)
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
        num_part = isempty(nf) ? 1 : _nest_binary(:*, nf)
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
function substitute_params_expr(expr, subs::AbstractDict)
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
function _substitute_sym_in_poly(p::POLY, target::Symbol, coeff, replacement)
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

"""Apply all parameter constraints sequentially to a polynomial.
When `binding_Ks` is provided, constraints between two binding K parameters
use the reciprocal coefficient (1/c instead of c) to correct for the K→1/K
inversion that happens later in the expression builder."""
function _apply_param_constraints(
    p::POLY, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
)
    for (target, coeff, factors) in constraints
        is_binding_to_binding = !isempty(binding_Ks) &&
            target ∈ binding_Ks &&
            all(f -> f[1] ∈ binding_Ks, factors)
        c = is_binding_to_binding ? 1 // coeff : coeff
        p = _substitute_sym_in_poly(p, target, c, factors)
    end
    p
end

# ─── Factored denominator types ──────────────────────────────

"""Product of POLY factors with integer exponents: prod(factors[i]^exponents[i])."""
struct FactoredPoly
    factors::Vector{POLY}
    exponents::Vector{Int}
end

"""Sum of coefficient * FactoredPoly terms: sum(coefficients[i] * products[i])."""
struct FactoredSigma
    coefficients::Vector{POLY}
    products::Vector{FactoredPoly}
end

"""One term of the factored denominator: sigma[g] * D[g]."""
struct DenomTerm
    sigma::FactoredSigma
    cofactor::POLY
end

"""Wrap a plain POLY as a degenerate DenomTerm (1 sub-group, 1 factor, exp 1)."""
function unfactored_denom_term(sigma_num::POLY, cofactor::POLY)
    DenomTerm(
        FactoredSigma([poly_one()], [FactoredPoly([sigma_num], [1])]),
        cofactor,
    )
end

# ─── Factored Expr generation ──────────────────────────────────

"""Convert a FactoredPoly to Expr: (f1)^e1 * (f2)^e2 * ..."""
function _factored_poly_to_expr(
    fp::FactoredPoly,
    param_syms::Set{Symbol},
    conc_syms::Set{Symbol},
    inverted_params::Set{Symbol},
)
    terms = Any[]
    for (f, e) in zip(fp.factors, fp.exponents)
        f_expr = _poly_to_expr(f, param_syms, conc_syms, inverted_params)
        if e == 1
            push!(terms, f_expr)
        else
            push!(terms, :(($f_expr) ^ $e))
        end
    end
    isempty(terms) ? 1 : _nest_binary(:*, terms)
end

"""Convert a FactoredSigma to Expr: sum of coeff * factored_poly."""
function _factored_sigma_to_expr(
    fs::FactoredSigma,
    param_syms::Set{Symbol},
    conc_syms::Set{Symbol},
    inverted_params::Set{Symbol},
)
    terms = Any[]
    for (coeff, fp) in zip(fs.coefficients, fs.products)
        fp_expr = _factored_poly_to_expr(
            fp, param_syms, conc_syms, inverted_params,
        )
        if coeff == poly_one()
            push!(terms, fp_expr)
        elseif fp_expr isa Integer && fp_expr == 1
            c_expr = _poly_to_expr(
                coeff, param_syms, conc_syms, inverted_params,
            )
            push!(terms, c_expr)
        else
            c_expr = _poly_to_expr(
                coeff, param_syms, conc_syms, inverted_params,
            )
            push!(terms, :($c_expr * $fp_expr))
        end
    end
    _nest_binary(:+, terms)
end

"""Convert Vector{DenomTerm} to a single denominator Expr."""
function _denom_terms_to_expr(
    terms::Vector{DenomTerm},
    param_syms::Set{Symbol},
    conc_syms::Set{Symbol},
    inverted_params::Set{Symbol},
)
    exprs = Any[]
    for dt in terms
        s_expr = _factored_sigma_to_expr(
            dt.sigma, param_syms, conc_syms, inverted_params,
        )
        if dt.cofactor == poly_one()
            push!(exprs, s_expr)
        else
            c_expr = _poly_to_expr(
                dt.cofactor, param_syms, conc_syms, inverted_params,
            )
            push!(exprs, :($s_expr * $c_expr))
        end
    end
    _nest_binary(:+, exprs)
end

"""Build rate Expr from numerator POLY and factored denominator terms."""
function to_rate_expr(
    num::POLY, denom_terms::Vector{DenomTerm},
    param_syms::Set{Symbol}, conc_syms::Set{Symbol},
    inverted_params::Set{Symbol}=Set{Symbol}(),
)
    num_expr = _poly_to_expr(num, param_syms, conc_syms, inverted_params)
    den_expr = _denom_terms_to_expr(
        denom_terms, param_syms, conc_syms, inverted_params,
    )
    :(E_total * ($num_expr) / ($den_expr))
end

# ─── Constraint application for factored types ────────────────

function _apply_param_constraints(
    fp::FactoredPoly, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
)
    FactoredPoly(
        [_apply_param_constraints(f, constraints; binding_Ks) for f in fp.factors],
        copy(fp.exponents),
    )
end

function _apply_param_constraints(
    fs::FactoredSigma, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
)
    FactoredSigma(
        [_apply_param_constraints(c, constraints; binding_Ks)
         for c in fs.coefficients],
        [_apply_param_constraints(fp, constraints; binding_Ks)
         for fp in fs.products],
    )
end

function _apply_param_constraints(
    dt::DenomTerm, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
)
    DenomTerm(
        _apply_param_constraints(dt.sigma, constraints; binding_Ks),
        _apply_param_constraints(dt.cofactor, constraints; binding_Ks),
    )
end

# ─── Expansion and estimation helpers ─────────────────────────

"""Expand a FactoredPoly to a flat POLY by multiplying out factors^exponents."""
function _expand_factored_poly(fp::FactoredPoly)::POLY
    result = poly_one()
    for (f, e) in zip(fp.factors, fp.exponents)
        p = f
        for _ in 2:e
            p = poly_mul(p, f)
        end
        result = poly_mul(result, p)
        if length(result) > MAX_RATE_EQUATION_TERMS
            error(
                "Rate equation for this mechanism has more than " *
                "$MAX_RATE_EQUATION_TERMS polynomial terms " *
                "(limit: $MAX_RATE_EQUATION_TERMS). Equations " *
                "this large take a very long time to compile " *
                "and are unlikely to be practically useful " *
                "for parameter fitting.",
            )
        end
    end
    result
end

"""Expand a FactoredSigma to a flat POLY."""
function _expand_factored_sigma(fs::FactoredSigma)::POLY
    result = poly_zero()
    for (coeff, fp) in zip(fs.coefficients, fs.products)
        result = poly_add(result, poly_mul(coeff, _expand_factored_poly(fp)))
    end
    result
end

"""Expand Vector{DenomTerm} to a single flat POLY denominator."""
function _expand_to_poly(terms::Vector{DenomTerm})::POLY
    result = poly_zero()
    for dt in terms
        sigma_poly = _expand_factored_sigma(dt.sigma)
        result = poly_add(result, poly_mul(sigma_poly, dt.cofactor))
    end
    result
end

"""
Estimate the expanded term count for factored denominator terms.
Upper bound: accounts for factor lengths and exponents but not cancellation.
"""
function _estimate_expanded_term_count(terms::Vector{DenomTerm})::Int
    total = 0
    for dt in terms
        sigma_count = 0
        for (coeff, fp) in zip(dt.sigma.coefficients, dt.sigma.products)
            fp_count = prod(
                length(f)^e for (f, e) in zip(fp.factors, fp.exponents)
            )
            sigma_count += length(coeff) * fp_count
        end
        total += sigma_count * max(length(dt.cofactor), 1)
    end
    total
end
