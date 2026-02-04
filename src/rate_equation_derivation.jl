using Combinatorics: permutations

# ─── Expression Utilities ─────────────────────────────────────────────────────
# (from expr_utils.jl)

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

# ─── Parameters API ───────────────────────────────────────────────────────────

"""
    parameters(m::EnzymeMechanism, [mode])

Return the parameter names required for the given mode as a tuple of Symbols.

# Modes
- `HaldaneWegscheider` (default): independent k's + Keq + E_total
- `Raw`: all 2N k's + E_total

# Examples
```julia
parameters(m)                          # independent k's + Keq + E_total
parameters(m, Raw)                     # all k's + E_total
parameters(m, HaldaneWegscheider)      # independent k's + Keq + E_total
```
"""
function parameters end

# Default: HaldaneWegscheider mode
parameters(m::EnzymeMechanism) = parameters(m, HaldaneWegscheider)

# Raw mode: all 2N k-parameters + E_total
@generated function parameters(::EnzymeMechanism{Species, Reactions}, ::RawMode) where {Species, Reactions}
    ks = ntuple(i -> Symbol("k", (i+1)÷2, isodd(i) ? "f" : "r"), 2 * length(Reactions))
    return (ks..., :E_total)
end

# HaldaneWegscheider mode: independent k's + Keq + E_total
@generated function parameters(::M, ::HaldaneWegscheiderMode) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    return (indep..., :Keq, :E_total)
end

# ─── Raw Rate Equation Derivation (King-Altman) ───────────────────────────────

"""
    _raw_symbolic_rate_expr(::Type{<:EnzymeMechanism})

Build a symbolic `Expr` for the QSSA rate equation using Laplacian cofactor
determinants expanded via the Leibniz formula, using all k parameters (no substitution).
"""
function _raw_symbolic_rate_expr(M::Type{<:EnzymeMechanism})
    m = M()
    subs = substrates(m)
    prods = products(m)
    enzs = enzyme_forms(m)
    rxns = reactions(m)

    isempty(subs) && error("No substrates defined")
    ref_name = subs[1][1]
    nu_ref = 0
    for (name, _) in subs
        name == ref_name && (nu_ref -= 1)
    end
    for (name, _) in prods
        name == ref_name && (nu_ref += 1)
    end
    nu_ref == 0 && error("Reference substrate has zero net stoichiometry")

    enz_names = Tuple(e[1] for e in enzs)
    enz_set = Set(enz_names)
    N = length(enzs)

    # 1. Build symbolic rate matrix R[i,j] as Expr (or 0 meaning absent)
    R = Matrix{Any}(fill(0, N, N))

    for (idx, (lhs, rhs)) in enumerate(rxns)
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        i = findfirst(==(e_lhs), enz_names)
        j = findfirst(==(e_rhs), enz_names)
        kf = Symbol("k$(idx)f")
        kr = Symbol("k$(idx)r")
        met_f = isempty(m_lhs) ? nothing : first(m_lhs)
        met_r = isempty(m_rhs) ? nothing : first(m_rhs)

        # Forward: i → j
        fwd = met_f === nothing ? make_param_accessor(kf) :
              :($(make_param_accessor(kf)) * $(make_conc_accessor(met_f)))
        R[i, j] = R[i, j] == 0 ? fwd : :($(R[i, j]) + $fwd)

        # Reverse: j → i
        rev = met_r === nothing ? make_param_accessor(kr) :
              :($(make_param_accessor(kr)) * $(make_conc_accessor(met_r)))
        R[j, i] = R[j, i] == 0 ? rev : :($(R[j, i]) + $rev)
    end

    # 2. Build symbolic Laplacian
    L = Matrix{Any}(fill(0, N, N))
    for i in 1:N
        diag_terms = Any[]
        for j in 1:N
            i == j && continue
            if R[i, j] != 0
                L[i, j] = :(-$(R[i, j]))
                push!(diag_terms, R[i, j])
            end
        end
        if !isempty(diag_terms)
            L[i, i] = length(diag_terms) == 1 ? diag_terms[1] : Expr(:call, :+, diag_terms...)
        end
    end

    # 3. Cofactor determinants via Leibniz formula with sparsity pruning
    @assert N <= 12 "Leibniz expansion is O(N!); N=$N exceeds the safety limit of 12"
    D = Vector{Any}(undef, N)
    for root in 1:N
        # Delete row root and col root
        rows = [r for r in 1:N if r != root]
        cols = [c for c in 1:N if c != root]
        n_sub = N - 1

        if n_sub == 0
            D[root] = 1
            continue
        end

        terms = Any[]
        for perm in permutations(1:n_sub)
            sign = _perm_sign(perm)

            # Product of L_sub[k, perm[k]] = L[rows[k], cols[perm[k]]]
            factors = Any[]
            all_nonzero = true
            for k in 1:n_sub
                entry = L[rows[k], cols[perm[k]]]
                if entry == 0
                    all_nonzero = false
                    break
                end
                push!(factors, entry)
            end
            !all_nonzero && continue

            prod_expr = length(factors) == 1 ? factors[1] : Expr(:call, :*, factors...)
            if sign == -1
                prod_expr = :(-$prod_expr)
            end
            push!(terms, prod_expr)
        end

        if isempty(terms)
            D[root] = 0
        elseif length(terms) == 1
            D[root] = terms[1]
        else
            D[root] = Expr(:call, :+, terms...)
        end
    end

    # 4. Denominator: sum of all cofactors
    nonzero_D = [D[k] for k in 1:N if D[k] != 0]
    denom = length(nonzero_D) == 1 ? nonzero_D[1] : Expr(:call, :+, nonzero_D...)

    # E_total
    et_expr = make_param_accessor(:E_total)

    # 5. Net consumption of reference substrate
    terms = Any[]
    for (idx, (lhs, rhs)) in enumerate(rxns)
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        i = findfirst(==(e_lhs), enz_names)
        j = findfirst(==(e_rhs), enz_names)
        kf = Symbol("k$(idx)f")
        kr = Symbol("k$(idx)r")
        met_f = isempty(m_lhs) ? nothing : first(m_lhs)
        met_r = isempty(m_rhs) ? nothing : first(m_rhs)

        rf_expr = met_f === nothing ? make_param_accessor(kf) :
                  :($(make_param_accessor(kf)) * $(make_conc_accessor(met_f)))
        rr_expr = met_r === nothing ? make_param_accessor(kr) :
                  :($(make_param_accessor(kr)) * $(make_conc_accessor(met_r)))
        flux = :($et_expr * ($rf_expr * $(D[i]) - $rr_expr * $(D[j])) / $denom)
        if met_f === ref_name
            push!(terms, flux)
        elseif met_r === ref_name
            push!(terms, :(-$flux))
        end
    end

    net_expr = isempty(terms) ? 0 :
               length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...)
    abs_nu = abs(nu_ref)
    abs_nu == 1 ? net_expr : :($net_expr / $abs_nu)
end

"""
    _symbolic_rate_expr(::Type{<:EnzymeMechanism})

Build a symbolic `Expr` for the QSSA rate equation with dependent parameters
substituted by expressions in terms of independent params + Keq.
"""
function _symbolic_rate_expr(M::Type{<:EnzymeMechanism})
    raw_expr = _raw_symbolic_rate_expr(M)
    dep_exprs, _ = _dependent_param_exprs(M)
    if !isempty(dep_exprs)
        raw_expr = substitute_params(raw_expr, dep_exprs)
    end
    return raw_expr
end

"""Compute sign of a permutation (as +1 or -1) using cycle decomposition."""
function _perm_sign(perm)
    n = length(perm)
    visited = falses(n)
    sign = 1
    for i in 1:n
        visited[i] && continue
        visited[i] = true
        cycle_len = 1
        j = perm[i]
        while j != i
            visited[j] = true
            cycle_len += 1
            j = perm[j]
        end
        if iseven(cycle_len)
            sign = -sign
        end
    end
    sign
end

# ─── Mode-dispatched rate_equation ────────────────────────────────────────────

"""
    rate_equation(m::EnzymeMechanism, params, concs, [mode])

Compute the QSSA steady-state rate (net consumption of the first substrate,
normalized by its stoichiometric coefficient). The body is generated at
compile time as a single arithmetic expression with no allocations, loops,
or matrix ops.

# Modes
- `HaldaneWegscheider` (default): Uses independent k-parameters + Keq (Haldane/Wegscheider constraints)
- `Raw`: Uses all 2N microscopic rate constants (no constraints)

# Parameters
The `params` NamedTuple must contain the parameters appropriate for the mode:
- `Raw`: all k's (k1f, k1r, k2f, k2r, ...) + E_total
- `HaldaneWegscheider`: independent k's + Keq + E_total
"""
function rate_equation end

# Default: HaldaneWegscheider mode
rate_equation(m::EnzymeMechanism, params, concs) =
    rate_equation(m, params, concs, HaldaneWegscheider)

# Raw mode: all 2N k-parameters
@generated function rate_equation(
    m::M, params::NamedTuple, concs::NamedTuple, ::RawMode
) where {M <: EnzymeMechanism}
    _raw_symbolic_rate_expr(M)
end

# HaldaneWegscheider mode: independent k-parameters + Keq
@generated function rate_equation(
    m::M, params::NamedTuple, concs::NamedTuple, ::HaldaneWegscheiderMode
) where {M <: EnzymeMechanism}
    _symbolic_rate_expr(M)
end

# ─── Polynomial representation for pretty-printing ────────────────────────────
#
# A monomial is represented as a sorted Vector{Symbol} (with repetition for
# powers), and a polynomial is a Dict mapping monomials to integer coefficients.
# Converting the Expr tree to this form automatically expands products over sums,
# cancels opposite terms (coefficients sum to zero), and sorts factors within
# each monomial.

const _Poly = Dict{Vector{Symbol},Int}

"""Sort key for symbols inside a monomial: k-constants before metabolites, then alphabetical."""
_monomial_sort_key(s) = (startswith(string(s), "k") ? 0 : 1, string(s))

"""
    _pmul(a, b) → _Poly

Multiply two polynomials by distributing every term pair.
"""
function _pmul(a::_Poly, b::_Poly)
    r = _Poly()
    for (k1, v1) in a, (k2, v2) in b
        k = sort([k1; k2]; by=_monomial_sort_key)
        r[k] = get(r, k, 0) + v1 * v2
    end
    filter!(p -> p[2] != 0, r)
end

"""
    _padd(a, b) → _Poly

Add two polynomials by merging coefficients. Zero-coefficient terms are dropped.
"""
function _padd(a::_Poly, b::_Poly)
    r = copy(a)
    for (k, v) in b
        r[k] = get(r, k, 0) + v
    end
    filter!(p -> p[2] != 0, r)
end

"""
    _to_poly(e) → _Poly

Recursively convert an `Expr` (that does not contain `/`) into a polynomial,
stripping `params.` and `concs.` prefixes so that only bare symbol names remain.
"""
function _to_poly(e)
    # Leaf: params.X or concs.X accessor
    if is_param_accessor(e) || is_conc_accessor(e)
        return _Poly([get_accessor_symbol(e)] => 1)
    end
    # Non-call or literal: treat as scalar constant
    if !(e isa Expr && e.head == :call)
        return _Poly(Symbol[] => (e isa Integer ? e : 1))
    end
    op = get_call_op(e)
    args = get_call_args(e)
    if op == :*
        return reduce(_pmul, _to_poly.(args))
    elseif op == :+
        return reduce(_padd, _to_poly.(args))
    elseif op == :- && length(args) == 1
        # Unary minus: negate all coefficients
        return _Poly(k => -v for (k, v) in _to_poly(args[1]))
    elseif op == :-
        # Binary minus: args[1] - args[2]
        return _padd(_to_poly(args[1]), _Poly(k => -v for (k, v) in _to_poly(args[2])))
    else
        return _Poly(Symbol[] => 1)
    end
end

"""
    _strip_div(e)

Recursively strip all `/` call nodes, keeping only the numerator side.
"""
function _strip_div(e)
    if is_call_expr(e, :/)
        return _strip_div(get_call_args(e)[1])
    elseif e isa Expr && e.head == :call
        return Expr(:call, e.args[1], _strip_div.(e.args[2:end])...)
    end
    return e
end

"""
    _find_denom(e)

Walk the expression tree and return the first denominator found inside a `/` call,
or `nothing` if no division is present.
"""
function _find_denom(e)
    if is_call_expr(e, :/)
        return get_call_args(e)[2]
    elseif e isa Expr && e.head == :call
        for arg in get_call_args(e)
            d = _find_denom(arg)
            d !== nothing && return d
        end
    end
    return nothing
end

"""
    _poly_str(p) → String

Pretty-print a polynomial: positive terms first, then negative, joined with +/-.
"""
function _poly_str(p::_Poly)
    isempty(p) && return "0"
    ts = sort(collect(p); by=x -> (x[2] < 0, x[1]))
    parts = String[]
    for (i, (k, v)) in enumerate(ts)
        m = if isempty(k)
            "$(abs(v))"
        elseif abs(v) == 1
            join(k, " * ")
        else
            "$(abs(v)) * " * join(k, " * ")
        end
        if i == 1
            push!(parts, v < 0 ? "-$m" : m)
        else
            push!(parts, v < 0 ? " - $m" : " + $m")
        end
    end
    return join(parts)
end

# ─── String Representation ────────────────────────────────────────────────────

"""
    rate_equation_string(m::EnzymeMechanism, [mode])

Return a string representation of the rate equation.

# Modes
- `HaldaneWegscheider` (default): Shows Haldane/Wegscheider constraints and independent parameters
- `Raw`: Shows raw equation with all 2N k-parameters
"""
function rate_equation_string end

# Default: HaldaneWegscheider mode
rate_equation_string(m::EnzymeMechanism) = rate_equation_string(m, HaldaneWegscheider)

"""
    _format_rate_equation_core(expr) → String

Format a rate expression as "v = E_total * (numerator) / (denominator)".
This is the common formatting logic shared by all modes.
"""
function _format_rate_equation_core(expr)
    num_poly = _Poly(filter(!=(:E_total), k) => v for (k, v) in _to_poly(_strip_div(expr)))
    denom_expr = _find_denom(expr)
    denom_poly = denom_expr === nothing ? _Poly(Symbol[] => 1) : _to_poly(denom_expr)
    "v = E_total * ($(_poly_str(num_poly))) / ($(_poly_str(denom_poly)))"
end

# Raw mode: all 2N k-parameters, no constraints
function rate_equation_string(::M, ::RawMode) where {M<:EnzymeMechanism}
    expr = _raw_symbolic_rate_expr(M)
    _format_rate_equation_core(expr)
end

# HaldaneWegscheider mode: with Haldane/Wegscheider constraints prepended
function rate_equation_string(::M, ::HaldaneWegscheiderMode) where {M<:EnzymeMechanism}
    expr = _raw_symbolic_rate_expr(M)
    eq = _format_rate_equation_core(expr)
    constraints = _constraint_expr_strings(M)
    if isempty(constraints)
        return eq
    else
        return join(constraints, "\n") * "\n\n" * eq
    end
end

# ─── Structural Identifiability ───────────────────────────────────────────────

"""
    _count_rate_monomials(M::Type{<:EnzymeMechanism}) → (n_num, n_denom)

Count distinct metabolite monomials in numerator and denominator of rate equation.

Returns the number of unique monomials (terms with distinct metabolite combinations,
ignoring k-parameters, Keq, and E_total) in each polynomial.

Note: Uses the raw symbolic expression (before Haldane substitution) to count monomials,
because Haldane substitution introduces Keq which doesn't change the structural form.
"""
function _count_rate_monomials(M::Type{<:EnzymeMechanism})
    # Use raw expression to avoid Keq appearing in monomials
    expr = _raw_symbolic_rate_expr(M)
    num_poly = _to_poly(_strip_div(expr))
    denom_expr = _find_denom(expr)
    denom_poly = denom_expr === nothing ? _Poly(Symbol[] => 1) : _to_poly(denom_expr)

    # Extract unique metabolite monomials (filter out k-parameters, Keq, and E_total)
    function metabolite_monomial(k::Vector{Symbol})
        sort(filter(s -> !is_k_parameter(s) && s != :E_total && s != :Keq, k))
    end

    num_monomials = unique([metabolite_monomial(k) for (k, v) in num_poly])
    denom_monomials = unique([metabolite_monomial(k) for (k, v) in denom_poly])

    return length(num_monomials), length(denom_monomials)
end

"""
    structural_identifiability_deficit(m::EnzymeMechanism) → Int

Number of excess parameters beyond what is structurally identifiable from kinetic data.

The deficit is computed as:
    n_k - (n_num + n_denom - 2)

where:
- n_k = number of independent k-parameters (after Haldane/Wegscheider constraints)
- n_num = number of distinct metabolite monomials in numerator
- n_denom = number of distinct metabolite monomials in denominator

The `-2` accounts for:
- `-1` from Haldane constraining the ratio of numerator coefficients (a_forward/a_reverse = Keq)
- `-1` from scaling freedom (can normalize one denominator coefficient)

Returns:
- `0` if exactly identifiable (parameters match identifiable coefficients)
- `> 0` if underdetermined (more parameters than identifiable coefficients)
- `< 0` if overdetermined (fewer parameters than identifiable coefficients)

# Examples
```julia
m_uu = make_uni_uni()  # E + S ⇌ ES ⇌ E + P
structural_identifiability_deficit(m_uu)  # Returns 0 (exactly identifiable)

m_3step = make_three_step_isomerization()  # E + S ⇌ ES ⇌ ES' ⇌ E + P
structural_identifiability_deficit(m_3step)  # Returns 2 (underdetermined)
```
"""
@generated function structural_identifiability_deficit(::M) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    n_k = length(indep)  # After Haldane/Wegscheider on k's

    n_num, n_denom = _count_rate_monomials(M)

    # Structurally identifiable coefficients:
    # - Numerator: n_num - 1 (Haldane constrains ratio a_forward/a_reverse = Keq)
    # - Denominator: n_denom - 1 (scaling freedom allows normalizing one coefficient)
    n_identifiable = (n_num - 1) + (n_denom - 1)

    return n_k - n_identifiable
end

"""
    is_identifiable(m::EnzymeMechanism) → Bool

Check if all mechanism parameters can be uniquely determined from steady-state kinetic data.
Returns `true` if structurally identifiable, `false` otherwise.

A mechanism is structurally identifiable when the number of independent k-parameters
is at most the number of structurally identifiable coefficients (deficit ≤ 0).

# Examples
```julia
m_uu = make_uni_uni()  # E + S ⇌ ES ⇌ E + P
is_identifiable(m_uu)  # Returns true

m_3step = make_three_step_isomerization()  # E + S ⇌ ES ⇌ ES' ⇌ E + P
is_identifiable(m_3step)  # Returns false
```
"""
is_identifiable(m::EnzymeMechanism) = structural_identifiability_deficit(m) <= 0
