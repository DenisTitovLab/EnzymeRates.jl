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

# Convert POLY to a Julia Expr for @generated function bodies (bare symbols)
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
    # Sort factors: monomials (1 term) first, then by term count descending
    order = sortperm(
        collect(zip(fp.factors, fp.exponents));
        by=((f, _),) -> (
            length(f) == 1 ? 0 : 1,
            -length(f),
            Tuple(sort!([string(s) for (mono, _) in f for (s, _) in mono])),
        ),
    )
    terms = Any[]
    for i in order
        f, e = fp.factors[i], fp.exponents[i]
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

"""Build rate Expr from factored numerator and factored denominator terms."""
function to_rate_expr(
    num::FactoredSigma, denom_terms::Vector{DenomTerm},
    param_syms::Set{Symbol}, conc_syms::Set{Symbol},
    inverted_params::Set{Symbol}=Set{Symbol}(),
)
    num_expr = _factored_sigma_to_expr(
        num, param_syms, conc_syms, inverted_params,
    )
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

# ─── OligomericEnzymeMechanism POLY helpers ──────────────────────

"""Rename all K/k symbols (not Keq) in a POLY with _T suffix."""
function _rename_poly_T(p::POLY)
    POLY(
        sort!(MONO([
            (is_k_parameter(s) && s != :Keq ? _rename_params_T(s) : s) => e
            for (s, e) in mono
        ]); by=first) => coeff
        for (mono, coeff) in p
    )
end

"""
Count distinct concentration monomials in the full oligomeric rate numerator and
denominator. Treats all K/k/L symbols as parameters (strips them from monomials).
Returns `(n_num, n_denom)`.
"""
function _count_oligomeric_rate_monomials(CM, CatN, RS, NConf)
    num_fs, denom_terms = _raw_symbolic_rate_polys(CM)
    N_cat_R = _expand_factored_sigma(num_fs)
    Q_cat_R = _expand_to_poly(denom_terms)

    # Build reg site partition polynomials (1 + ligand_sym for each ligand)
    reg_Q_R = POLY[
        reduce(poly_add, (poly_add(poly_one(), poly_sym(lig)) for lig in ligs))
        for (ligs, _) in RS
    ]

    function num_poly_for_conf(N_cat, Q_cat, reg_Qs, L_factor)
        n_term = poly_mul(N_cat, _poly_power(Q_cat, CatN - 1))
        for (idx, (_, n_reg)) in enumerate(RS)
            n_reg == CatN || continue
            n_term = poly_mul(n_term, _poly_power(reg_Qs[idx], n_reg))
        end
        L_factor === nothing ? n_term : poly_mul(poly_sym(L_factor), n_term)
    end

    function den_poly_for_conf(Q_cat, reg_Qs, L_factor)
        d_term = _poly_power(Q_cat, CatN)
        for (idx, (_, n_reg)) in enumerate(RS)
            d_term = poly_mul(d_term, _poly_power(reg_Qs[idx], n_reg))
        end
        L_factor === nothing ? d_term : poly_mul(poly_sym(L_factor), d_term)
    end

    full_num = num_poly_for_conf(N_cat_R, Q_cat_R, reg_Q_R, nothing)
    full_den = den_poly_for_conf(Q_cat_R, reg_Q_R, nothing)

    if NConf == 2
        N_cat_T = _rename_poly_T(N_cat_R)
        Q_cat_T = _rename_poly_T(Q_cat_R)
        reg_Q_T = POLY[_rename_poly_T(q) for q in reg_Q_R]

        full_num = poly_add(full_num, num_poly_for_conf(N_cat_T, Q_cat_T, reg_Q_T, :L))
        full_den = poly_add(full_den, den_poly_for_conf(Q_cat_T, reg_Q_T, :L))
    end

    # Count distinct concentration monomials (strip all k/K/L params)
    conc_mono(mono) = sort!(MONO([
        s => e for (s, e) in mono
        if !is_k_parameter(s) && s != :E_total && s != :Keq && s != :L
    ]); by=first)

    n_num = length(unique(conc_mono(k) for (k, _) in full_num))
    n_denom = length(unique(conc_mono(k) for (k, _) in full_den))
    n_num, n_denom
end
