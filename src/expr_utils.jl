"""
Utilities for working with Julia expression trees.

Provides common operations for pattern matching and transforming expressions
that reference `params.X` and `concs.X` accessors used throughout the rate
equation code generation.
"""

# ─── Pattern matching predicates ─────────────────────────────────────────────

"""
    is_param_accessor(expr) → Bool

Check if `expr` is a `params.X` accessor pattern.
"""
function is_param_accessor(expr)
    expr isa Expr && expr.head == :. && length(expr.args) == 2 &&
    expr.args[1] == :params && expr.args[2] isa QuoteNode
end

"""
    is_conc_accessor(expr) → Bool

Check if `expr` is a `concs.X` accessor pattern.
"""
function is_conc_accessor(expr)
    expr isa Expr && expr.head == :. && length(expr.args) == 2 &&
    expr.args[1] == :concs && expr.args[2] isa QuoteNode
end

"""
    get_accessor_symbol(expr) → Symbol

Extract the symbol from a `params.X` or `concs.X` accessor.
Assumes `is_param_accessor(expr)` or `is_conc_accessor(expr)` is true.
"""
function get_accessor_symbol(expr)
    expr.args[2].value
end

"""
    is_k_parameter(sym::Symbol) → Bool

Check if a symbol is a k-parameter (rate constant), not Keq or other params.
"""
function is_k_parameter(sym::Symbol)
    s = string(sym)
    startswith(s, "k") && sym != :Keq
end

"""
    is_call_expr(expr, op::Symbol) → Bool

Check if `expr` is a function call with operator `op`.
"""
function is_call_expr(expr, op::Symbol)
    expr isa Expr && expr.head == :call && length(expr.args) >= 1 && expr.args[1] == op
end

"""
    is_arithmetic_call(expr) → Bool

Check if `expr` is an arithmetic operation (+, -, *, /, ^).
"""
function is_arithmetic_call(expr)
    expr isa Expr && expr.head == :call && length(expr.args) >= 2 &&
    expr.args[1] in (:+, :-, :*, :/, :^)
end

"""
    get_call_op(expr) → Symbol

Get the operator from a call expression. Assumes `expr.head == :call`.
"""
get_call_op(expr) = expr.args[1]

"""
    get_call_args(expr) → Vector

Get the arguments from a call expression (excluding the operator).
"""
get_call_args(expr) = expr.args[2:end]

# ─── Expression construction helpers ─────────────────────────────────────────

"""
    make_param_accessor(sym::Symbol) → Expr

Create a `params.sym` accessor expression.
"""
make_param_accessor(sym::Symbol) = :(params.$sym)

"""
    make_conc_accessor(sym::Symbol) → Expr

Create a `concs.sym` accessor expression.
"""
make_conc_accessor(sym::Symbol) = :(concs.$sym)

"""
    make_product(terms) → Expr

Create a product expression from a list of terms.
Returns `1` for empty list, the single term for length 1, or `*(terms...)` otherwise.
"""
function make_product(terms)
    isempty(terms) && return :(1)
    length(terms) == 1 && return terms[1]
    return Expr(:call, :*, terms...)
end

"""
    make_sum(terms) → Expr

Create a sum expression from a list of terms.
Returns `0` for empty list, the single term for length 1, or `+(terms...)` otherwise.
"""
function make_sum(terms)
    isempty(terms) && return :(0)
    length(terms) == 1 && return terms[1]
    return Expr(:call, :+, terms...)
end

"""
    make_division(num, denom) → Expr

Create a division expression `num / denom`.
"""
make_division(num, denom) = :($num / $denom)

"""
    make_power(base, exp::Integer) → Expr

Create a power expression with simplification for common cases.
"""
function make_power(base, exp::Integer)
    exp == 0 && return :(1)
    exp == 1 && return base
    exp == -1 && return :(1 / $base)
    exp > 0 && return :($base ^ $exp)
    return :(1 / $base ^ $(-exp))
end

"""
    make_power(base, exp::Rational) → Expr

Create a power expression for rational exponents.
"""
function make_power(base, exp::Rational)
    exp == 0 && return :(1)
    exp == 1 && return base
    exp == -1 && return :(1 / $base)
    if denominator(exp) == 1
        return make_power(base, Int(numerator(exp)))
    end
    # General rational exponent - use Float64
    return :($base ^ $(Float64(exp)))
end

# ─── Expression transformation ───────────────────────────────────────────────

"""
    map_expr(f, expr)

Recursively apply `f` to all subexpressions, bottom-up.
`f` receives each node after its children have been transformed.
"""
function map_expr(f, expr)
    if expr isa Expr
        new_args = [map_expr(f, a) for a in expr.args]
        return f(Expr(expr.head, new_args...))
    else
        return f(expr)
    end
end

"""
    transform_expr(should_transform, transformer, expr)

Transform expression nodes that match `should_transform` predicate using `transformer`.
Non-matching nodes are recursively processed but not transformed.
"""
function transform_expr(should_transform::Function, transformer::Function, expr)
    if should_transform(expr)
        return transformer(expr)
    elseif expr isa Expr
        new_args = Any[expr.args[1]]  # Keep head/operator
        for a in expr.args[2:end]
            push!(new_args, transform_expr(should_transform, transformer, a))
        end
        return Expr(expr.head, new_args...)
    else
        return expr
    end
end

"""
    substitute_params(expr, subs::Dict{Symbol, Expr})

Recursively substitute `params.dep_sym` with its expression from `subs`.
"""
function substitute_params(expr, subs::Dict{Symbol, Expr})
    if is_param_accessor(expr)
        sym = get_accessor_symbol(expr)
        return get(subs, sym, expr)
    elseif expr isa Expr
        new_args = Any[substitute_params(a, subs) for a in expr.args]
        return Expr(expr.head, new_args...)
    else
        return expr
    end
end
