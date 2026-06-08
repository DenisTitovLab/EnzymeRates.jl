# ABOUTME: Lightweight symbolic polynomial type (POLY = Dict{MONO, Rational{Int}})
# ABOUTME: for compile-time rate equation derivation.

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

"""
Reduce num/den to lowest terms over concentrations only. For each concentration
symbol, shift its exponent in every monomial of num ∪ den so the minimum exponent
across all those monomials is 0 — clears any 1/conc coupling without ever
touching a parameter symbol. A concentration absent from a monomial counts as
exponent 0, so the shift is the identity unless the symbol appears in *every*
monomial (min > 0) or with a negative exponent somewhere (min < 0). Identity for
sequential mechanisms (constant term present); for a 1/conc-coupled mechanism it
clears the coupling, yielding the standard division-free form.
"""
function _reduce_conc_lowest_terms(num::POLY, den::POLY, conc_set::Set{Symbol})
    isempty(conc_set) && return num, den
    mins = Dict{Symbol,Int}()
    for p in (num, den), mono in keys(p)
        present = Dict{Symbol,Int}(s => e for (s, e) in mono if s in conc_set)
        for s in conc_set
            e = get(present, s, 0)
            mins[s] = haskey(mins, s) ? min(mins[s], e) : e
        end
    end
    filter!(p -> p.second != 0, mins)
    isempty(mins) && return num, den
    function shift(p)
        out = POLY()
        for (mono, v) in p
            md = Dict{Symbol,Int}(mono)
            for (s, mn) in mins
                md[s] = get(md, s, 0) - mn
            end
            filter!(pr -> pr.second != 0, md)
            out[sort!(MONO(collect(md)); by=first)] = v
        end
        out
    end
    shift(num), shift(den)
end

"""
Cofactor determinant expansion for symbolic matrices. Checks
intermediate term count against `MAX_RATE_EQUATION_TERMS` to abort
early for mechanisms whose rate equations would be too large.
"""
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

"""Convert `POLY` to a Julia `Expr` for `@generated` function bodies (bare symbols)."""
function _poly_to_expr(p::POLY, param_syms::Set{Symbol}, conc_syms::Set{Symbol})
    isempty(p) && return 0
    pos, neg = Any[], Any[]
    sorted = sort(
        collect(p);
        by=x -> (
            sum(e for (s,e) in x[1] if s ∉ param_syms; init=0),
            x[2] < 0,
            Tuple((string(s), e) for (s, e) in x[1]),
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
            tgt, ex = e > 0 ? (nf, e) : (df, -e)
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

"""
Build a balanced binary `+`/`*` tree so every emitted call has exactly two
operands. Required for zero-allocation `rate_equation` runtime: Julia inlines
binary `+(::Float64, ::Float64)` into fused scalar arithmetic, but falls back
to a varargs path that boxes the operand tuple once the chain exceeds ~30
terms. See `test_rate_equation_performance` for the contract this enforces.
"""
function _nest_binary(op::Symbol, terms::Vector{Any})
    n = length(terms)
    n == 1 && return terms[1]
    n == 2 && return Expr(:call, op, terms[1], terms[2])
    mid = n >> 1
    Expr(:call, op,
        _nest_binary(op, terms[1:mid]),
        _nest_binary(op, terms[mid+1:end]))
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

"""Substitute symbols in an `Expr` tree (bare symbol matching)."""
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
absent keys are left unchanged. Used by the allosteric derivation to rename
A-state symbols to their I-state counterparts when building the inactive-
state polynomial (e.g., `:K_A_ATP_E → :K_I_ATP_E`).
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

"""
Set of symbols that survive the inactive-state masking applied by
`_zero_symbols_in_poly(_, a_only_syms)`. A symbol `s` survives iff it
appears in at least one monomial of `num_A ∪ den_A` that contains
NO symbol from `a_only_syms` (those monomials are the ones that remain
non-zero in the inactive-state polynomial).

Used by the AllostericEnzymeMechanism dep-exprs filter to avoid declaring
`:NonequalAI` I-state parameters that never appear in the rate equation
body — a phantom-parameter case where `p` is only present in active-state
monomials that get zeroed when constructing the inactive-state polynomial.
"""
function _i_state_surviving_syms(num_A::POLY, den_A::POLY, a_only_syms::Set{Symbol})
    surviving = Set{Symbol}()
    for poly in (num_A, den_A)
        for mono in keys(poly)
            any(s ∈ a_only_syms for (s, _) in mono) && continue
            for (s, _) in mono
                push!(surviving, s)
            end
        end
    end
    surviving
end
