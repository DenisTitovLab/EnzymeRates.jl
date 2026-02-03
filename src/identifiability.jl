using LinearAlgebra: dot, I

"""
Structural identifiability analysis for enzyme rate equations.

Detects scaling symmetries in rate equations after Haldane/Wegscheider substitution.
Parameters that can be arbitrarily scaled while preserving the rate equation are
non-identifiable; only certain combinations of these parameters are identifiable.

The approach finds parameter transformations `k_i → k_i * c^θ_i` that leave
the rate equation N/D invariant. This is done by:
1. Extracting monomials from the rate equation (after dependent param substitution)
2. Building the monomial exponent matrix A where A[m,k] = exponent of parameter k in monomial m
3. For rate invariance, all monomials must scale identically: A * θ = α * 1
4. Computing the null space of the difference matrix (rows of A differing from first row)
5. Null space vectors = non-identifiable scaling directions

Known limitations:
- Misses discrete/local identifiability (finitely many equivalent parameter sets)
- Doesn't detect permutation symmetries in symmetric mechanisms
- Assumes all concentrations can be varied independently
"""

# ─── Monomial extraction from symbolic expressions ───────────────────────────

"""
    _extract_monomials_from_poly(expr) → Set{Dict{Symbol,Int}}

Extract all monomials from a polynomial expression (no division).
Each monomial is represented as Dict{Symbol,Int} mapping k-parameter symbols to exponents.
Only extracts k-parameters (symbols starting with 'k'), ignoring concentrations and other params.
"""
function _extract_monomials_from_poly(expr)
    monomials = Set{Dict{Symbol,Int}}()
    _collect_monomials!(monomials, expr, Dict{Symbol,Int}())
    return monomials
end

"""
    _collect_monomials!(monomials, expr, current)

Recursively collect monomials from an expression tree into `monomials` set.
`current` tracks the k-parameter exponents for the monomial being built.
"""
function _collect_monomials!(monomials::Set{Dict{Symbol,Int}}, expr, current::Dict{Symbol,Int})
    # Handle params.X accessor
    if is_param_accessor(expr)
        sym = get_accessor_symbol(expr)
        if is_k_parameter(sym)
            new_current = copy(current)
            new_current[sym] = get(new_current, sym, 0) + 1
            push!(monomials, new_current)
        else
            # Non-k parameter (Keq, E_total) - treat as constant
            push!(monomials, copy(current))
        end
        return
    end

    # Handle concs.X accessor - concentrations are external, treat as constant
    if is_conc_accessor(expr)
        push!(monomials, copy(current))
        return
    end

    # Handle numeric literals and non-expressions
    if !(expr isa Expr) || expr.head != :call
        push!(monomials, copy(current))
        return
    end

    op = get_call_op(expr)
    args = get_call_args(expr)

    if op == :+
        # Addition: collect from each term separately
        for arg in args
            _collect_monomials!(monomials, arg, current)
        end
    elseif op == :-
        if length(args) == 1
            # Unary minus: doesn't change monomials
            _collect_monomials!(monomials, args[1], current)
        else
            # Binary minus: collect from both sides
            _collect_monomials!(monomials, args[1], current)
            _collect_monomials!(monomials, args[2], current)
        end
    elseif op == :*
        # Multiplication: combine all factors
        _collect_product_monomials!(monomials, args, current)
    elseif op == :^
        # Power: x^n where n is positive integer
        base, exp = args
        if exp isa Integer && exp > 0
            repeated_args = fill(base, exp)
            _collect_product_monomials!(monomials, repeated_args, current)
        else
            push!(monomials, copy(current))
        end
    else
        # Unknown operation - treat as constant
        push!(monomials, copy(current))
    end
end

"""
    _collect_product_monomials!(monomials, factors, current)

Collect monomials from a product of factors using Cartesian product of factor monomials.
"""
function _collect_product_monomials!(monomials::Set{Dict{Symbol,Int}}, factors, current::Dict{Symbol,Int})
    if isempty(factors)
        push!(monomials, copy(current))
        return
    end

    # Get monomials from each factor
    factor_monomials = [_extract_monomials_from_poly(f) for f in factors]

    # Compute Cartesian product
    result_monomials = [Dict{Symbol,Int}()]
    for fms in factor_monomials
        new_results = Dict{Symbol,Int}[]
        for rm in result_monomials
            for fm in fms
                combined = copy(rm)
                for (sym, exp) in fm
                    combined[sym] = get(combined, sym, 0) + exp
                end
                push!(new_results, combined)
            end
        end
        result_monomials = new_results
    end

    # Combine with current monomial and add to set
    for rm in result_monomials
        combined = copy(current)
        for (sym, exp) in rm
            combined[sym] = get(combined, sym, 0) + exp
        end
        push!(monomials, combined)
    end
end

"""
    _extract_all_monomials(M::Type{<:EnzymeMechanism}) → Set{Dict{Symbol,Int}}

Extract all unique monomials (in terms of k-parameters) from the rate equation
after Haldane/Wegscheider substitution.
"""
function _extract_all_monomials(M::Type{<:EnzymeMechanism})
    expr = _symbolic_rate_expr(M)

    # Separate numerator and denominator
    num_expr = _strip_div(expr)
    denom_expr = _find_denom(expr)

    all_monomials = Set{Dict{Symbol,Int}}()

    union!(all_monomials, _extract_monomials_from_poly(num_expr))

    if denom_expr !== nothing
        union!(all_monomials, _extract_monomials_from_poly(denom_expr))
    end

    # Remove empty monomials (constant terms have no k-parameters)
    filter!(m -> !isempty(m), all_monomials)

    return all_monomials
end

# ─── Core identifiability analysis ───────────────────────────────────────────

"""
    _monomial_exponent_matrix(M::Type{<:EnzymeMechanism}) → (A, k_params, all_monomials)

Build matrix A where A[m,k] = exponent of k-parameter in monomial m.
Only independent k-parameters (after Haldane/Wegscheider) appear in the matrix.

Returns:
- `A`: Integer matrix of shape (n_monomials, n_params)
- `k_params`: Sorted vector of independent k-parameter symbols
- `all_monomials`: Vector of monomial dictionaries
"""
function _monomial_exponent_matrix(M::Type{<:EnzymeMechanism})
    # Get independent k-parameters (excludes dependent k's, Keq, E_total)
    _, indep = _dependent_param_exprs(M)
    k_params = sort!(Symbol[s for s in indep])

    # Extract all monomials from the rate expression
    monomials_set = _extract_all_monomials(M)
    all_monomials = collect(monomials_set)

    # Build exponent matrix
    n_monomials = length(all_monomials)
    n_params = length(k_params)

    A = zeros(Int, n_monomials, n_params)
    for (m_idx, mono) in enumerate(all_monomials)
        for (k_idx, k) in enumerate(k_params)
            A[m_idx, k_idx] = get(mono, k, 0)
        end
    end

    return A, k_params, all_monomials
end

"""
    _scaling_null_space(M::Type{<:EnzymeMechanism}) → Matrix{Int}

Compute the null space of the scaling constraint matrix.

For a rate equation to be invariant under scaling k_i → k_i * c^θ_i,
all monomials must scale by the same factor. If monomial m has exponent
A[m,k] for parameter k, then monomial m scales as c^(Σ_k A[m,k]*θ_k).
For invariance, this sum must be equal for all monomials, which means
(A[i,:] - A[j,:]) · θ = 0 for all pairs (i,j).

Returns a matrix where each column is a non-identifiable scaling direction.
An empty matrix (0 columns) means fully identifiable.
"""
function _scaling_null_space(M::Type{<:EnzymeMechanism})
    A, k_params, _ = _monomial_exponent_matrix(M)

    n_monomials, n_params = size(A)

    # Edge cases: fully identifiable
    (n_params == 0 || n_monomials <= 1) && return zeros(Int, n_params, 0)

    # Build difference matrix: D[i,:] = A[i+1,:] - A[1,:]
    D = zeros(Int, n_monomials - 1, n_params)
    for i in 1:(n_monomials - 1)
        for k in 1:n_params
            D[i, k] = A[i + 1, k] - A[1, k]
        end
    end

    return _integer_nullspace(D)
end

# ─── Option C: Simplified identifiable combinations (split into sub-functions) ─

"""
    _classify_parameters(NS::Matrix{Int}, n_params::Int) → (individually_ident, involved_in_nonident)

Classify parameters into individually identifiable (not in any null space vector)
and involved in non-identifiability (appear in at least one null space vector).

Returns:
- `individually_ident`: BitVector where true means parameter is individually identifiable
- `involved_in_nonident`: BitVector where true means parameter is involved in non-identifiability
"""
function _classify_parameters(NS::Matrix{Int}, n_params::Int)
    n_null = size(NS, 2)

    individually_ident = trues(n_params)
    involved_in_nonident = falses(n_params)

    for j in 1:n_null
        for i in 1:n_params
            if NS[i, j] != 0
                individually_ident[i] = false
                involved_in_nonident[i] = true
            end
        end
    end

    return individually_ident, involved_in_nonident
end

"""
    _find_identifiable_basis(NS_involved::Matrix{Int}, n_involved::Int) → Matrix{Int}

Find a basis for the identifiable subspace orthogonal to the null space.
Uses integer Gram-Schmidt style orthogonalization.

Returns a matrix where each column is an identifiable combination exponent vector
(in terms of the involved parameters only).
"""
function _find_identifiable_basis(NS_involved::Matrix{Int}, n_involved::Int)
    n_null = size(NS_involved, 2)
    n_ident = n_involved - n_null

    n_ident <= 0 && return zeros(Int, n_involved, 0)

    # Strategy: Find n_ident vectors that, together with NS_involved columns,
    # span the full space. We do this by finding vectors orthogonal to NS_involved.

    # Use rational arithmetic for exact computation
    NS_rat = Matrix{Rational{BigInt}}(NS_involved)

    # Find the row space of NS (which vectors span)
    # The orthogonal complement gives us identifiable directions

    # Build augmented matrix [NS' ; I] and find its null space
    # Actually, simpler: find vectors v such that NS' * v = 0
    # These are orthogonal to all null space vectors

    # Use the approach: extend NS columns to full basis, complement is identifiable
    result = zeros(Int, n_involved, n_ident)

    # Project out the null space directions from standard basis vectors
    # Keep those that remain linearly independent

    # Start with identity matrix columns
    candidates = Matrix{Rational{BigInt}}(I, n_involved, n_involved)

    # Gram-Schmidt against null space vectors
    for j in 1:n_null
        ns_vec = NS_rat[:, j]
        ns_norm_sq = dot(ns_vec, ns_vec)
        ns_norm_sq == 0 && continue
        for i in 1:n_involved
            # Project out the null space component
            proj_coeff = dot(candidates[:, i], ns_vec) / ns_norm_sq
            candidates[:, i] .-= proj_coeff .* ns_vec
        end
    end

    # Find n_ident linearly independent columns from candidates
    selected = Int[]
    for i in 1:n_involved
        length(selected) >= n_ident && break
        # Check if candidates[:,i] is linearly independent of selected
        col = candidates[:, i]
        all(col .== 0) && continue

        if isempty(selected)
            push!(selected, i)
        else
            # Check linear independence
            prev_cols = candidates[:, selected]
            # Simple check: if adding this column increases rank
            aug = hcat(prev_cols, col)
            if _rational_rank(aug) > length(selected)
                push!(selected, i)
            end
        end
    end

    # Convert selected columns to integers
    for (k, col_idx) in enumerate(selected)
        col = candidates[:, col_idx]
        # Clear denominators
        denoms = [denominator(x) for x in col if x != 0]
        if !isempty(denoms)
            scale = lcm(denoms...)
            col = col .* scale
        end
        # Convert to Int and normalize
        int_col = Int.(round.(col))
        g = gcd(filter(!=(0), abs.(int_col))...)
        g > 0 && (int_col .÷= g)
        # Make first nonzero positive
        first_nz = findfirst(!=(0), int_col)
        if first_nz !== nothing && int_col[first_nz] < 0
            int_col .*= -1
        end
        result[:, k] = int_col
    end

    return result
end

"""Compute rank of a rational matrix."""
function _rational_rank(A::Matrix{<:Rational})
    m, n = size(A)
    work = copy(A)
    rank = 0
    col = 1
    for row in 1:m
        col > n && break
        # Find pivot
        piv = 0
        for r in row:m
            if work[r, col] != 0
                piv = r
                break
            end
        end
        if piv == 0
            col += 1
            continue
        end
        work[row, :], work[piv, :] = work[piv, :], work[row, :]
        work[row, :] ./= work[row, col]
        for r in 1:m
            r == row && continue
            if work[r, col] != 0
                work[r, :] .-= work[r, col] .* work[row, :]
            end
        end
        rank += 1
        col += 1
    end
    return rank
end

"""
    _name_combination(exponents::AbstractVector{Int}, params::Vector{Symbol}) → Symbol

Generate a name for a parameter combination based on its exponents.
Products use underscore separator, ratios use "_over_".

Examples:
- [1, 1, 0] with [:k1f, :k2f, :k3f] → :k1f_k2f
- [1, 0, -1] with [:k1f, :k2f, :k3r] → :k1f_over_k3r
"""
function _name_combination(exponents::AbstractVector{Int}, params::Vector{Symbol})
    num_parts = Symbol[]
    denom_parts = Symbol[]

    for (i, exp) in enumerate(exponents)
        if exp > 0
            for _ in 1:exp
                push!(num_parts, params[i])
            end
        elseif exp < 0
            for _ in 1:(-exp)
                push!(denom_parts, params[i])
            end
        end
    end

    if isempty(num_parts) && isempty(denom_parts)
        return :_const
    elseif isempty(denom_parts)
        return Symbol(join(num_parts, "_"))
    elseif isempty(num_parts)
        return Symbol("1_over_", join(denom_parts, "_"))
    else
        return Symbol(join(num_parts, "_"), "_over_", join(denom_parts, "_"))
    end
end

"""
    _build_combination_expr(exponents::AbstractVector{Int}, params::Vector{Symbol}) → Expr

Build an expression for a parameter combination.
"""
function _build_combination_expr(exponents::AbstractVector{Int}, params::Vector{Symbol})
    num_terms = Expr[]
    denom_terms = Expr[]

    for (i, exp) in enumerate(exponents)
        if exp > 0
            for _ in 1:exp
                push!(num_terms, make_param_accessor(params[i]))
            end
        elseif exp < 0
            for _ in 1:(-exp)
                push!(denom_terms, make_param_accessor(params[i]))
            end
        end
    end

    num_expr = make_product(num_terms)
    denom_expr = make_product(denom_terms)

    if isempty(denom_terms)
        return num_expr
    elseif isempty(num_terms)
        return make_division(:(1), denom_expr)
    else
        return make_division(num_expr, denom_expr)
    end
end

"""
    _identifiable_combinations_from_nullspace(NS, k_params) → Vector{Tuple{Symbol,Expr,Vector{Int}}}

Convert null space to identifiable parameter combinations.
Returns triples of (name, expression, exponent_vector).

This is the main entry point that uses the sub-functions above.
"""
function _identifiable_combinations_from_nullspace(NS::Matrix{Int}, k_params::Vector{Symbol})
    n_params = length(k_params)
    n_null = size(NS, 2)

    # If no null space, all parameters are individually identifiable
    if n_null == 0
        combinations = Tuple{Symbol,Expr,Vector{Int}}[]
        for (i, k) in enumerate(k_params)
            exp_vec = zeros(Int, n_params)
            exp_vec[i] = 1
            push!(combinations, (k, make_param_accessor(k), exp_vec))
        end
        return combinations
    end

    # Classify parameters
    individually_ident, involved_in_nonident = _classify_parameters(NS, n_params)

    combinations = Tuple{Symbol,Expr,Vector{Int}}[]

    # Add individually identifiable parameters
    for (i, k) in enumerate(k_params)
        if individually_ident[i]
            exp_vec = zeros(Int, n_params)
            exp_vec[i] = 1
            push!(combinations, (k, make_param_accessor(k), exp_vec))
        end
    end

    # Find identifiable combinations for involved parameters
    involved_idx = findall(involved_in_nonident)
    n_involved = length(involved_idx)

    if n_involved == 0
        return combinations
    end

    # Extract null space for involved parameters only
    NS_involved = NS[involved_idx, :]

    # Find identifiable basis
    ident_basis = _find_identifiable_basis(NS_involved, n_involved)

    # Convert basis vectors to named combinations
    involved_params = k_params[involved_idx]
    for col_idx in 1:size(ident_basis, 2)
        basis_vec = ident_basis[:, col_idx]

        # Expand to full parameter space
        exp_vec = zeros(Int, n_params)
        for (j, idx) in enumerate(involved_idx)
            exp_vec[idx] = basis_vec[j]
        end

        name = _name_combination(basis_vec, involved_params)
        expr = _build_combination_expr(basis_vec, involved_params)
        push!(combinations, (name, expr, exp_vec))
    end

    return combinations
end

# ─── Identifiable basis for expression transformation ────────────────────────

"""
    _identifiable_basis(M::Type{<:EnzymeMechanism})

Return (basis_names, basis_matrix, k_params) where:
- basis_names: Vector{Symbol} of identifiable parameter names
- basis_matrix: Matrix{Int} where column i is the exponent vector for basis_names[i]
- k_params: Vector{Symbol} of original k-parameter names (sorted)
"""
function _identifiable_basis(M::Type{<:EnzymeMechanism})
    NS = _scaling_null_space(M)
    _, k_params, _ = _monomial_exponent_matrix(M)

    combinations = _identifiable_combinations_from_nullspace(NS, k_params)

    basis_names = Symbol[c[1] for c in combinations]

    n_params = length(k_params)
    n_ident = length(combinations)
    if n_ident == 0
        basis_matrix = zeros(Int, n_params, 0)
    else
        basis_matrix = zeros(Int, n_params, n_ident)
        for (j, c) in enumerate(combinations)
            basis_matrix[:, j] = c[3]
        end
    end

    return basis_names, basis_matrix, k_params
end

# ─── Expression transformation to identifiable form ──────────────────────────

"""
    _exponent_to_identifiable_expr(exp_vec, basis_names, basis_matrix)

Convert an exponent vector (in terms of original k-params) to an expression
using identifiable parameter names.

Uses exact rational arithmetic to solve basis_matrix * c = exp_vec.
"""
function _exponent_to_identifiable_expr(exp_vec::Vector{Int}, basis_names::Vector{Symbol}, basis_matrix::Matrix{Int})
    n_params, n_ident = size(basis_matrix)

    # Solve basis_matrix * c = exp_vec using rational arithmetic
    B = Matrix{Rational{BigInt}}(basis_matrix)
    e = Vector{Rational{BigInt}}(exp_vec)

    # Linear solve - convert result back to Rational to avoid BigFloat
    c_raw = B \ e
    c = Vector{Rational{BigInt}}(c_raw)

    # Build expression from coefficients
    terms = Any[]
    for (i, name) in enumerate(basis_names)
        ci = c[i]
        ci == 0 && continue

        base = make_param_accessor(name)
        if ci == 1
            push!(terms, base)
        elseif ci == -1
            push!(terms, make_division(:(1), base))
        elseif denominator(ci) == 1
            exp_int = Int(numerator(ci))
            push!(terms, make_power(base, exp_int))
        else
            # Fractional exponent (shouldn't happen for well-formed rate equations)
            push!(terms, :($base ^ $(Float64(ci))))
        end
    end

    return make_product(terms)
end

"""
    _transform_to_identifiable(expr, basis_names, basis_matrix, k_to_idx, n_params)

Transform a rate expression to use identifiable parameter names.
"""
function _transform_to_identifiable(expr, basis_names::Vector{Symbol}, basis_matrix::Matrix{Int},
                                     k_to_idx::Dict{Symbol,Int}, n_params::Int)
    # Handle params.X accessor for k-parameters
    if is_param_accessor(expr)
        sym = get_accessor_symbol(expr)
        if haskey(k_to_idx, sym)
            exp_vec = zeros(Int, n_params)
            exp_vec[k_to_idx[sym]] = 1
            return _exponent_to_identifiable_expr(exp_vec, basis_names, basis_matrix)
        else
            return expr  # Keq, E_total - keep as-is
        end
    end

    # Non-expression - return as-is
    expr isa Expr || return expr

    # Handle multiplication specially - combine k-parameters
    if is_call_expr(expr, :*)
        return _transform_product(get_call_args(expr), basis_names, basis_matrix, k_to_idx, n_params)
    end

    # For other expressions, recurse on arguments
    new_args = Any[expr.args[1]]  # Keep the operator
    for a in expr.args[2:end]
        push!(new_args, _transform_to_identifiable(a, basis_names, basis_matrix, k_to_idx, n_params))
    end
    return Expr(expr.head, new_args...)
end

"""
    _transform_product(factors, basis_names, basis_matrix, k_to_idx, n_params)

Transform a product expression, combining k-parameters into identifiable form.
"""
function _transform_product(factors, basis_names::Vector{Symbol}, basis_matrix::Matrix{Int},
                            k_to_idx::Dict{Symbol,Int}, n_params::Int)
    # Collect k-parameter exponents and non-k factors separately
    exp_vec = zeros(Int, n_params)
    other_factors = Any[]

    for f in factors
        k_exp = _extract_k_exponents(f, k_to_idx, n_params)
        if k_exp !== nothing
            exp_vec .+= k_exp
        else
            # Recursively transform non-k factors
            push!(other_factors, _transform_to_identifiable(f, basis_names, basis_matrix, k_to_idx, n_params))
        end
    end

    # Convert k-parameter product to identifiable form
    if any(exp_vec .!= 0)
        k_expr = _exponent_to_identifiable_expr(exp_vec, basis_names, basis_matrix)
        pushfirst!(other_factors, k_expr)
    end

    return make_product(other_factors)
end

"""
    _extract_k_exponents(expr, k_to_idx, n_params)

Extract k-parameter exponents from an expression that is a pure k-product.
Returns the exponent vector, or `nothing` if expr contains non-k factors.
"""
function _extract_k_exponents(expr, k_to_idx::Dict{Symbol,Int}, n_params::Int)
    # Single k-parameter
    if is_param_accessor(expr)
        sym = get_accessor_symbol(expr)
        if haskey(k_to_idx, sym)
            exp_vec = zeros(Int, n_params)
            exp_vec[k_to_idx[sym]] = 1
            return exp_vec
        end
        return nothing  # Non-k param
    end

    # Product of k-parameters
    if is_call_expr(expr, :*)
        total_exp = zeros(Int, n_params)
        for f in get_call_args(expr)
            sub_exp = _extract_k_exponents(f, k_to_idx, n_params)
            if sub_exp === nothing
                return nothing  # Contains non-k factors
            end
            total_exp .+= sub_exp
        end
        return total_exp
    end

    return nothing
end

"""
    _identifiable_symbolic_rate_expr(M::Type{<:EnzymeMechanism})

Build a symbolic Expr for the rate equation using identifiable parameter combinations.
"""
function _identifiable_symbolic_rate_expr(M::Type{<:EnzymeMechanism})
    # Get HaldaneWegscheider expression
    hw_expr = _symbolic_rate_expr(M)

    # Check if fully identifiable
    NS = _scaling_null_space(M)
    if size(NS, 2) == 0
        return hw_expr  # Already in identifiable form
    end

    # Build identifiable basis and transform
    basis_names, basis_matrix, k_params = _identifiable_basis(M)
    k_to_idx = Dict(k => i for (i, k) in enumerate(k_params))
    n_params = length(k_params)

    return _transform_to_identifiable(hw_expr, basis_names, basis_matrix, k_to_idx, n_params)
end

# ─── Mode-dispatched functions for IdentifiableHaldaneWegscheider ────────────

@generated function rate_equation(
    m::M, params::NamedTuple, concs::NamedTuple, ::IdentifiableHaldaneWegscheiderMode
) where {M <: EnzymeMechanism}
    _identifiable_symbolic_rate_expr(M)
end

@generated function parameters(::M, ::IdentifiableHaldaneWegscheiderMode) where {M <: EnzymeMechanism}
    NS = _scaling_null_space(M)
    _, k_params, _ = _monomial_exponent_matrix(M)
    combinations = _identifiable_combinations_from_nullspace(NS, k_params)
    ident_names = Tuple(c[1] for c in combinations)
    return (ident_names..., :Keq, :E_total)
end

# ─── Mode-dispatched rate_equation_string for IdentifiableHaldaneWegscheider ─

function rate_equation_string(::M, ::IdentifiableHaldaneWegscheiderMode) where {M<:EnzymeMechanism}
    NS = _scaling_null_space(M)
    _, k_params, _ = _monomial_exponent_matrix(M)

    # Get HaldaneWegscheider string as base
    hw_str = rate_equation_string(M(), HaldaneWegscheider)

    if size(NS, 2) == 0
        return hw_str  # Fully identifiable
    end

    # Build identifiability annotation
    combinations = _identifiable_combinations_from_nullspace(NS, k_params)

    lines = String[]
    push!(lines, "# Identifiable parameter combinations:")
    for (name, expr, _) in combinations
        s = string(expr)
        s = replace(s, "params." => "")
        push!(lines, "#   $name = $s")
    end

    push!(lines, "#")
    push!(lines, "# Non-identifiable scaling directions:")
    for j in 1:size(NS, 2)
        parts = String[]
        for (i, k) in enumerate(k_params)
            if NS[i, j] != 0
                exp = NS[i, j]
                if exp == 1
                    push!(parts, "$k -> $k * c")
                elseif exp == -1
                    push!(parts, "$k -> $k / c")
                elseif exp > 0
                    push!(parts, "$k -> $k * c^$exp")
                else
                    push!(parts, "$k -> $k / c^$(-exp)")
                end
            end
        end
        push!(lines, "#   Direction $j: " * join(parts, ", "))
    end

    push!(lines, "")
    push!(lines, hw_str)

    return join(lines, "\n")
end

# ─── Public API ──────────────────────────────────────────────────────────────

"""
    is_identifiable(m::EnzymeMechanism) → Bool

Return true if all independent parameters are structurally identifiable
(no non-trivial scaling symmetries exist).
"""
@generated function is_identifiable(::M) where {M<:EnzymeMechanism}
    NS = _scaling_null_space(M)
    return size(NS, 2) == 0
end

"""
    non_identifiable_directions(m::EnzymeMechanism) → Matrix{Int}

Return matrix where each column is a non-identifiable scaling direction.
Column i means: parameters can scale as k_j -> k_j * c^(matrix[j,i]).

The rows correspond to parameters in sorted order (use `independent_parameters(m)`).
An empty matrix (0 columns) indicates all parameters are identifiable.
"""
@generated function non_identifiable_directions(::M) where {M<:EnzymeMechanism}
    NS = _scaling_null_space(M)
    return NS
end

"""
    identifiable_combinations(m::EnzymeMechanism) → Vector{Tuple{Symbol, String}}

Return pairs of (combination_name, expression_string) for all identifiable
parameter combinations.
"""
@generated function identifiable_combinations(::M) where {M<:EnzymeMechanism}
    NS = _scaling_null_space(M)
    _, k_params, _ = _monomial_exponent_matrix(M)
    combinations = _identifiable_combinations_from_nullspace(NS, k_params)
    result = Tuple{Symbol,String}[]
    for (name, expr, _) in combinations
        s = string(expr)
        s = replace(s, "params." => "")
        push!(result, (name, s))
    end
    return result
end
