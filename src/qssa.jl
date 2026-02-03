using Combinatorics: permutations

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

# ─── Mode-dispatched rate_equation ───────────────────────────────────────────

"""
    rate_equation(m::EnzymeMechanism, params, concs, [mode])

Compute the QSSA steady-state rate (net consumption of the first substrate,
normalized by its stoichiometric coefficient). The body is generated at
compile time as a single arithmetic expression with no allocations, loops,
or matrix ops.

# Modes
- `IdentifiableHaldaneWegscheider` (default): Uses identifiable parameter combinations + Keq
- `HaldaneWegscheider`: Uses independent k-parameters + Keq (Haldane/Wegscheider constraints)
- `Raw`: Uses all 2N microscopic rate constants (no constraints)

# Parameters
The `params` NamedTuple must contain the parameters appropriate for the mode:
- `Raw`: all k's (k1f, k1r, k2f, k2r, ...) + E_total
- `HaldaneWegscheider`: independent k's + Keq + E_total
- `IdentifiableHaldaneWegscheider`: identifiable combinations + Keq + E_total
"""
function rate_equation end

# Default: IdentifiableHaldaneWegscheider mode
rate_equation(m::EnzymeMechanism, params, concs) =
    rate_equation(m, params, concs, IdentifiableHaldaneWegscheider)

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

# IdentifiableHaldaneWegscheider mode: identifiable combinations + Keq
# (Defined in identifiability.jl after the required functions are available)

# ─── Polynomial representation for pretty-printing ──────────────────────────
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

# ─── Consolidated rate_equation_string (Option D) ────────────────────────────

"""
    rate_equation_string(m::EnzymeMechanism, [mode])

Return a string representation of the rate equation.

# Modes
- `IdentifiableHaldaneWegscheider` (default): Shows identifiable parameters and constraints
- `HaldaneWegscheider`: Shows Haldane/Wegscheider constraints and independent parameters
- `Raw`: Shows raw equation with all 2N k-parameters
"""
function rate_equation_string end

# Default: IdentifiableHaldaneWegscheider mode
rate_equation_string(m::EnzymeMechanism) = rate_equation_string(m, IdentifiableHaldaneWegscheider)

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

# IdentifiableHaldaneWegscheider mode: defined in identifiability.jl

"""
    constraint_strings(m::EnzymeMechanism)

Return a tuple of human-readable constraint strings for the dependent parameters.
"""
@generated function constraint_strings(::M) where {M<:EnzymeMechanism}
    strs = _constraint_expr_strings(M)
    return Tuple(strs)
end
