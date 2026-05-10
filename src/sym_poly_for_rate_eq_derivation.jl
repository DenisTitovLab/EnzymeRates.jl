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

function _mono_op(a::MONO, b::MONO, sign::Int)
    d = Dict{Symbol,Int}()
    for (s, e) in a; d[s] = get(d, s, 0) + e; end
    for (s, e) in b; d[s] = get(d, s, 0) + sign * e; end
    filter!(p -> p.second != 0, d)
    sort!(MONO(collect(d)); by=first)
end
_mono_mul(a::MONO, b::MONO) = _mono_op(a, b, 1)
_mono_div(a::MONO, b::MONO) = _mono_op(a, b, -1)

"""Raise a POLY to a non-negative integer power via repeated multiplication."""
function _poly_power(p::POLY, n::Int)
    n == 0 && return poly_one()
    result = poly_one()
    for _ in 1:n
        result = poly_mul(result, p)
    end
    result
end

"""Divide POLY by a single-term POLY (exact division, assumes divisibility)."""
function _poly_div_mono(p::POLY, divisor::POLY)::POLY
    m = first(keys(divisor))
    POLY(_mono_div(k, m) => v for (k, v) in p)
end

"""
Try exact polynomial division `dividend ÷ divisor`.
Returns the quotient POLY, or `nothing` if division has a remainder.
Requires the divisor to have a constant term of 1 (binding polynomials).
Uses metabolite-degree layering to guarantee termination with negative exponents.
"""
function _try_poly_exact_div(dividend::POLY, divisor::POLY)
    isempty(divisor) && return nothing
    const_mono = MONO()
    haskey(divisor, const_mono) && divisor[const_mono] == 1 ||
        return nothing

    # Identify metabolite symbols (positive-exponent symbols in divisor)
    met_syms = Set{Symbol}()
    for (dm, _) in divisor
        for (s, e) in dm
            e > 0 && push!(met_syms, s)
        end
    end
    isempty(met_syms) && return copy(dividend)

    met_degree(mono) = sum(e for (s, e) in mono if s ∈ met_syms; init=0)

    # Non-constant part of divisor (BigInt for overflow safety)
    T = Rational{BigInt}
    σ_prime = [(k, T(v)) for (k, v) in divisor if k != const_mono]

    # Group dividend by metabolite degree
    by_deg = Dict{Int,Dict{MONO,T}}()
    for (m, c) in dividend
        d = met_degree(m)
        dd = get!(by_deg, d, Dict{MONO,T}())
        dd[m] = T(c)
    end

    # Process degrees low→high: degree-d Q terms are the remainder after
    # subtracting contributions from lower-degree Q terms via σ'
    quot = Dict{MONO,T}()
    for d in sort!(collect(keys(by_deg)))
        remaining = copy(by_deg[d])
        for (qm, qc) in quot
            for (dm, dc) in σ_prime
                pm = _mono_mul(qm, dm)
                met_degree(pm) == d || continue
                remaining[pm] = get(remaining, pm, T(0)) - qc * dc
            end
        end
        filter!(p -> p.second != 0, remaining)
        merge!(quot, remaining)
    end

    result = POLY(k => Rational{Int}(v) for (k, v) in quot if v != 0)
    poly_mul(result, divisor) == dividend ? result : nothing
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

# Convert POLY to a Julia Expr for @generated function bodies (bare symbols).
# `inverted_params` is retained as an optional parameter for the surviving
# allosteric-branch callers (`_factored_sigma_to_expr`, etc.) and the 4-arg
# `_poly_to_expr` overload that wraps a POLY as `FactoredSigma`. Commit 2
# rewrites those callers to consume flat POLYs and drops `inverted_params`.
function _poly_to_expr(p::POLY, param_syms::Set{Symbol}, conc_syms::Set{Symbol},
                       inverted_params::Set{Symbol}=Set{Symbol}())
    isempty(p) && return 0
    pos, neg = Any[], Any[]
    sorted = sort(
        collect(p);
        by=x -> (
            sum(e for (s,e) in x[1] if s ∉ param_syms; init=0),
            x[2] < 0,
            Tuple(string(s) for (s, _) in x[1]),
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
                tgt, ex = e > 0 ? (nf, e) : (df, -e)
            end
            ex != 0 && push!(tgt, ex == 1 ? s : :($s ^ $ex))
        end
        num_part = isempty(nf) ? 1 : _nest_binary(:*, nf)
        term = isempty(df) ? num_part : :($num_part / $(_nest_binary(:*, df)))
        coeff > 0 ? push!(pos, term) : push!(neg, term)
    end
    pe = isempty(pos) ? nothing : _nest_binary(:+, pos)
    ne = isempty(neg) ? nothing : _nest_binary(:+, neg)
    pe !== nothing && ne !== nothing && return :($pe - $ne)
    pe !== nothing && return pe
    ne !== nothing && return :(- $ne)
    return 0
end

"""Build flat n-ary expression: +(a, b, c, d). Single term returns unwrapped."""
function _nest_binary(op::Symbol, terms::Vector{Any})
    length(terms) == 1 ? terms[1] : Expr(:call, op, terms...)
end

"""Operator precedence for Expr→String conversion."""
_op_prec(op::Symbol) =
    op in (:+, :-) ? 1 : op in (:*, :/) ? 2 : op == :^ ? 3 : 0

"""Wrap in parens if inner expression has lower precedence."""
function _str_paren(expr, threshold; right=false)
    s = _expr_to_string(expr)
    expr isa Expr && expr.head == :call || return s
    ip = _op_prec(expr.args[1])
    (right ? ip <= threshold : ip < threshold) && ip > 0 ?
        "($s)" : s
end

"""
Precedence-aware Expr→String conversion for rate equations.
Avoids unnecessary parentheses that Julia's `string()` adds.
"""
function _expr_to_string(x)
    x isa Union{Number, Symbol} && return string(x)
    x isa Expr && x.head == :call || return string(x)
    op, args = x.args[1], @view(x.args[2:end])
    # Unary minus
    op == :- && length(args) == 1 &&
        return "-$(_str_paren(args[1], 1))"
    if op == :+
        return join((_expr_to_string(a) for a in args), " + ")
    elseif op == :-
        return "$(_expr_to_string(args[1])) - " *
               _str_paren(args[2], 1; right=true)
    elseif op == :*
        parts = [begin
            s = _expr_to_string(a)
            need_parens = a isa Expr && a.head == :call &&
                (_op_prec(a.args[1]) < 2 || a.args[1] == :/)
            need_parens ? "($s)" : s
        end for a in args]
        return join(parts, " * ")
    elseif op == :/
        return "$(_str_paren(args[1], 2)) / " *
               _str_paren(args[2], 2; right=true)
    elseif op == :^
        return "$(_str_paren(args[1], 3; right=true)) ^ " *
               _expr_to_string(args[2])
    else
        return string(x)
    end
end

"""Check if a symbol is a rate/equilibrium parameter (k or K), not Keq or E_total."""
function is_k_parameter(sym::Symbol)
    sym in (:Keq, :E_total) && return false
    s = string(sym)
    startswith(s, "k") || (startswith(s, "K") && length(s) > 1 && isdigit(s[2]))
end

"""
    build_power_expr(keq_exp::Rational, factors)

Build a power-product expression for thermodynamic-constraint
substitution of the form `Keq^keq_exp * prod(k_i^exp_i)` from
`factors`, an iterable of `(k_i, exp_i)` pairs. Used by the
constraint solver to materialize the Keq×rate-constant power product
that replaces a dependent rate constant.
"""
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

"""Check if an expression references any symbol in the given set."""
function _expr_references_any(expr, syms::Set{Symbol})
    if expr isa Symbol
        return expr ∈ syms
    elseif expr isa Expr
        return any(_expr_references_any(a, syms) for a in expr.args)
    end
    false
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

# ─── Symbol renaming in POLY ───────────────────────────────

"""
Rename symbols in a polynomial. `rename_map` is a `Dict{Symbol, Symbol}`;
absent keys are left unchanged. Used to alias non-representative kinetic-group
parameter symbols to their representative (e.g., `K2 → K1` when steps 1 and 2
share a kinetic group).
"""
function _rename_symbols(p::POLY, rename_map::AbstractDict{Symbol, Symbol})
    isempty(rename_map) && return p
    result = POLY()
    for (mono, val) in p
        new_mono = sort!(
            MONO([get(rename_map, s, s) => e for (s, e) in mono]);
            by=first,
        )
        # Combine like-monomial entries by exponent merging
        combined = Dict{Symbol, Int}()
        for (s, e) in new_mono
            combined[s] = get(combined, s, 0) + e
        end
        filter!(p -> p.second != 0, combined)
        canon = sort!(MONO(collect(combined)); by=first)
        result[canon] = get(result, canon, 0) + val
    end
    filter!(p -> p.second != 0, result)
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
    inverted_params::Set{Symbol}=Set{Symbol}(),
)
    # Sort factors: monomials (1 term) first, then by term count descending
    order = sortperm(
        collect(zip(fp.factors, fp.exponents));
        by=((f, _),) -> (
            length(f) == 1 ? 0 : 1,
            -length(f),
            Tuple(sort!([string(s) for (mono, _) in f for (s, _) in mono])),
        ),
    )
    terms = map(order) do i
        f_expr = _poly_to_expr(
            fp.factors[i], param_syms, conc_syms, inverted_params,
        )
        fp.exponents[i] == 1 ? f_expr : :(($f_expr) ^ $(fp.exponents[i]))
    end
    isempty(terms) ? 1 : _nest_binary(:*, Any[terms...])
end

"""Convert a FactoredSigma to Expr: sum of coeff * factored_poly."""
function _factored_sigma_to_expr(
    fs::FactoredSigma,
    param_syms::Set{Symbol},
    conc_syms::Set{Symbol},
    inverted_params::Set{Symbol}=Set{Symbol}(),
)
    terms = map(fs.coefficients, fs.products) do coeff, fp
        fp_expr = _factored_poly_to_expr(
            fp, param_syms, conc_syms, inverted_params,
        )
        if coeff == poly_one()
            fp_expr
        else
            c_expr = _poly_to_expr(
                coeff, param_syms, conc_syms, inverted_params,
            )
            fp_expr isa Integer && fp_expr == 1 ? c_expr :
                :($c_expr * $fp_expr)
        end
    end
    _nest_binary(:+, Any[terms...])
end

"""Convert Vector{DenomTerm} to a single denominator Expr."""
function _denom_terms_to_expr(
    dts::Vector{DenomTerm},
    param_syms::Set{Symbol},
    conc_syms::Set{Symbol},
    inverted_params::Set{Symbol}=Set{Symbol}(),
)
    exprs = map(dts) do dt
        s = _factored_sigma_to_expr(
            dt.sigma, param_syms, conc_syms, inverted_params,
        )
        if dt.cofactor == poly_one()
            s
        else
            c = _poly_to_expr(
                dt.cofactor, param_syms, conc_syms, inverted_params,
            )
            :($s * $c)
        end
    end
    _nest_binary(:+, Any[exprs...])
end

"""Build rate Expr from numerator (POLY or FactoredSigma) and factored denominator."""
function to_rate_expr(
    num::Union{POLY, FactoredSigma},
    denom_terms::Vector{DenomTerm},
    param_syms::Set{Symbol}, conc_syms::Set{Symbol},
    inverted_params::Set{Symbol}=Set{Symbol}(),
)
    num_expr = num isa POLY ?
        _poly_to_expr(num, param_syms, conc_syms, inverted_params) :
        _factored_sigma_to_expr(
            num, param_syms, conc_syms, inverted_params,
        )
    den_expr = _denom_terms_to_expr(
        denom_terms, param_syms, conc_syms, inverted_params,
    )
    :(E_total * ($num_expr) / ($den_expr))
end

# ─── Symbol renaming for factored types ───────────────────────

function _rename_symbols(fp::FactoredPoly, rename_map::AbstractDict{Symbol, Symbol})
    FactoredPoly(
        [_rename_symbols(f, rename_map) for f in fp.factors],
        copy(fp.exponents),
    )
end

function _rename_symbols(fs::FactoredSigma, rename_map::AbstractDict{Symbol, Symbol})
    FactoredSigma(
        [_rename_symbols(c, rename_map) for c in fs.coefficients],
        [_rename_symbols(fp, rename_map) for fp in fs.products],
    )
end

function _rename_symbols(dt::DenomTerm, rename_map::AbstractDict{Symbol, Symbol})
    DenomTerm(
        _rename_symbols(dt.sigma, rename_map),
        _rename_symbols(dt.cofactor, rename_map),
    )
end

# ─── Expansion and estimation helpers ─────────────────────────

"""Expand a FactoredPoly to a flat POLY by multiplying out factors^exponents."""
function _expand_factored_poly(fp::FactoredPoly)::POLY
    result = poly_one()
    for (f, e) in zip(fp.factors, fp.exponents)
        result = poly_mul(result, _poly_power(f, e))
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

# ─── AllostericEnzymeMechanism POLY helpers ──────────────────────

"""Remove monomials containing any of the given symbols from a POLY."""
function _zero_symbols_in_poly(p::POLY, sym_set::Set{Symbol})
    isempty(sym_set) && return p
    result = POLY()
    for (mono, coeff) in p
        has_sym = any(s ∈ sym_set for (s, _) in mono)
        has_sym || (result[mono] = coeff)
    end
    result
end
