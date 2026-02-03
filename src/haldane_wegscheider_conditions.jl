"""
Haldane and Wegscheider thermodynamic constraints for enzyme mechanisms.

Identifies thermodynamic cycles in the mechanism graph (via null-space of the
enzyme incidence matrix), classifies them as Haldane (net reaction) or
Wegscheider (internal loop), and performs Gaussian elimination to express a
minimal set of dependent rate constants in terms of the independent ones and
the equilibrium constant Keq.
"""

"""
    _split_reaction_side(side, enz_set) → (enzyme_sym, met_syms)

Partition a reaction side (tuple of Symbols) into the single enzyme-form
symbol and a vector of metabolite symbols.
"""
function _split_reaction_side(side, enz_set)
    enzyme_sym = first(s for s in side if s in enz_set)
    met_syms = Symbol[s for s in side if s ∉ enz_set]
    return enzyme_sym, met_syms
end

"""
    _enzyme_incidence_matrix(M::Type{<:EnzymeMechanism})

Build the enzyme incidence matrix B (N_states × N_steps).
B[i,j] = -1 if step j consumes enzyme state i (LHS), +1 if step j produces it (RHS).
"""
function _enzyme_incidence_matrix(M::Type{<:EnzymeMechanism})
    m = M()
    enzs = enzyme_forms(m)
    rxns = reactions(m)
    enz_names = Tuple(e[1] for e in enzs)
    enz_set = Set(enz_names)
    N = length(enzs)
    nsteps = length(rxns)
    B = zeros(Int, N, nsteps)
    for (j, (lhs, rhs)) in enumerate(rxns)
        e_lhs, _ = _split_reaction_side(lhs, enz_set)
        e_rhs, _ = _split_reaction_side(rhs, enz_set)
        i_lhs = findfirst(==(e_lhs), enz_names)
        i_rhs = findfirst(==(e_rhs), enz_names)
        B[i_lhs, j] -= 1
        B[i_rhs, j] += 1
    end
    return B
end

"""
    _integer_nullspace(A::Matrix{Int})

Compute an integer basis for the null space of integer matrix A using
row echelon form with exact arithmetic (Rational).
Returns a Matrix{Int} where each column is a basis vector, or an empty 0-column matrix.

Note: `LinearAlgebra.nullspace` operates in Float64, which cannot guarantee the
exact integer results needed for symbolic constraint expressions.
"""
function _integer_nullspace(A::Matrix{Int})
    m, n = size(A)
    # Work in rationals for exact arithmetic
    R = Matrix{Rational{BigInt}}(A)
    pivot_cols = Int[]
    row = 1
    for col in 1:n
        # Find pivot row
        piv = 0
        for r in row:m
            if R[r, col] != 0
                piv = r
                break
            end
        end
        piv == 0 && continue
        # Swap rows
        R[row, :], R[piv, :] = R[piv, :], R[row, :]
        # Scale pivot row
        R[row, :] ./= R[row, col]
        # Eliminate column
        for r in 1:m
            r == row && continue
            if R[r, col] != 0
                R[r, :] .-= R[r, col] .* R[row, :]
            end
        end
        push!(pivot_cols, col)
        row += 1
    end
    # Free variables
    free_cols = setdiff(1:n, pivot_cols)
    nullity = length(free_cols)
    nullity == 0 && return zeros(Int, n, 0)
    # Build null space basis
    NS = zeros(Rational{BigInt}, n, nullity)
    for (k, fc) in enumerate(free_cols)
        NS[fc, k] = 1
        for (r, pc) in enumerate(pivot_cols)
            NS[pc, k] = -R[r, fc]
        end
    end
    # Convert to integers by clearing denominators per column
    result = zeros(Int, n, nullity)
    for k in 1:nullity
        lcm_val = lcm(denominator.(NS[:, k])...)
        for i in 1:n
            result[i, k] = Int(NS[i, k] * lcm_val)
        end
        # Normalize: make first nonzero entry positive, divide by GCD
        g = gcd(abs.(result[:, k])...)
        g > 0 && (result[:, k] .÷= g)
        first_nz = findfirst(!=(0), result[:, k])
        if first_nz !== nothing && result[first_nz, k] < 0
            result[:, k] .*= -1
        end
    end
    return result
end

"""
    _thermodynamic_constraints(M::Type{<:EnzymeMechanism})

Find all thermodynamic constraints (Haldane + Wegscheider).
Returns `(C, rhs_coeffs)`:
- `C`: n_constraints × N_steps integer matrix. Each row is a cycle vector over steps.
- `rhs_coeffs`: integer vector of length n_constraints.
  For Haldane rows: `sum_j C[i,j] * (log(k_jf) - log(k_jr)) = rhs_coeffs[i] * log(Keq)`.
  For Wegscheider rows: `rhs_coeffs[i] = 0`.
"""
function _thermodynamic_constraints(M::Type{<:EnzymeMechanism})
    B = _enzyme_incidence_matrix(M)
    NS = _integer_nullspace(B)
    n_constraints = size(NS, 2)
    n_constraints == 0 && return zeros(Int, 0, size(B, 2)), Int[]

    S = stoich_matrix(M())
    # Overall net stoichiometry vector (from species lists)
    m = M()
    subs = substrates(m)
    prods = products(m)
    mets = metabolites(m)
    met_names = [mt[1] for mt in mets]
    nu_net = zeros(Int, length(met_names))
    for (name, _) in subs
        idx = findfirst(==(name), met_names)
        nu_net[idx] -= 1
    end
    for (name, _) in prods
        idx = findfirst(==(name), met_names)
        nu_net[idx] += 1
    end

    C = NS'  # n_constraints × N_steps (rows are cycle vectors)
    rhs_coeffs = zeros(Int, n_constraints)
    for i in 1:n_constraints
        cycle = C[i, :]
        nu_cycle = S * cycle  # net metabolite change for this cycle
        # Check if nu_cycle is proportional to nu_net
        if all(nu_cycle .== 0)
            rhs_coeffs[i] = 0  # Wegscheider
        else
            # Find proportionality constant c: nu_cycle = c * nu_net
            c = nothing
            valid = true
            for j in eachindex(nu_cycle)
                if nu_net[j] == 0
                    if nu_cycle[j] != 0
                        valid = false
                        break
                    end
                else
                    c_j = nu_cycle[j] // nu_net[j]
                    if c === nothing
                        c = c_j
                    elseif c_j != c
                        valid = false
                        break
                    end
                end
            end
            if valid && c !== nothing && denominator(c) == 1
                rhs_coeffs[i] = Int(c)
            else
                error("Cycle $i produces metabolite change not proportional to net reaction")
            end
        end
    end
    return C, rhs_coeffs
end

"""
    _free_enzyme_binding_steps(M::Type{<:EnzymeMechanism})

Return indices of steps where free enzyme binds/releases a metabolite.
These steps are preferred to keep independent (deprioritized as pivots).
"""
function _free_enzyme_binding_steps(M::Type{<:EnzymeMechanism})
    m = M()
    enzs = enzyme_forms(m)
    rxns = reactions(m)
    enz_set = Set(e[1] for e in enzs)

    # Free enzyme forms: those with empty atoms
    free_enz_set = Set(e[1] for e in enzs if isempty(e[2]))

    binding_steps = Int[]
    for (j, (lhs, rhs)) in enumerate(rxns)
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        has_met = !isempty(m_lhs) || !isempty(m_rhs)
        involves_free = e_lhs in free_enz_set || e_rhs in free_enz_set
        if has_met && involves_free
            push!(binding_steps, j)
        end
    end
    return binding_steps
end

# ─── Sub-functions for _dependent_param_exprs ───────────────────────────────

"""
    _build_log_constraint_matrix(C, rhs_coeffs, nsteps)

Expand the cycle matrix C (over steps) into a constraint matrix A over individual
rate constants (2*nsteps columns: kf and kr for each step).

Returns `(A, rhs)` where A is n_constraints × 2*nsteps Rational matrix and
rhs is the Rational right-hand-side vector.
"""
function _build_log_constraint_matrix(C, rhs_coeffs, nsteps)
    n_constraints = size(C, 1)
    A = zeros(Rational{BigInt}, n_constraints, 2 * nsteps)
    rhs = Rational{BigInt}.(rhs_coeffs)
    for i in 1:n_constraints
        for j in 1:nsteps
            A[i, 2j-1] = C[i, j]    # coefficient for log(k_jf)
            A[i, 2j] = -C[i, j]     # coefficient for log(k_jr)
        end
    end
    return A, rhs
end

"""
    _pivot_priority(nsteps, free_binding, step_has_met_lhs, step_has_met_rhs)

Compute a priority vector over 2*nsteps rate-constant columns.  Higher priority
means the variable is more preferred as a pivot (i.e. more likely to become dependent).

Tiers (from highest to lowest priority):
  - Tier 3 (base 20): internal isomerisation steps (no metabolite on either side)
  - Tier 2 (base 10): metabolite-involving steps that do NOT touch free enzyme
  - Tier 1 (base  0): free-enzyme binding/release steps (keep independent)

Within each tier, +1 bonus for the direction (kf/kr) that does NOT face a metabolite,
so that the more physically interpretable binding rate stays independent.
"""
function _pivot_priority(nsteps, free_binding, step_has_met_lhs, step_has_met_rhs)
    priority = zeros(Int, 2 * nsteps)
    for j in 1:nsteps
        if j in free_binding
            base = 0   # Tier 1: free-enzyme binding — keep independent
        elseif !(step_has_met_lhs[j] || step_has_met_rhs[j])
            base = 20  # Tier 3: internal isomerisation — prefer to eliminate
        else
            base = 10  # Tier 2: metabolite step, not free-enzyme
        end
        # kf (2j-1): bonus if no metabolite on LHS (binding rate is more interpretable)
        priority[2j-1] = base + (step_has_met_lhs[j] ? 0 : 1)
        # kr (2j):   bonus if no metabolite on RHS
        priority[2j] = base + (step_has_met_rhs[j] ? 0 : 1)
    end
    return priority
end

"""
    _pivoted_gaussian_elimination(A, rhs, priority)

Perform Gaussian elimination on the constraint matrix A with column pivoting
guided by `priority`.  Returns `(work_A, work_rhs, pivot_cols)`.
"""
function _pivoted_gaussian_elimination(A, rhs, priority)
    n_constraints = size(A, 1)
    ncols = size(A, 2)
    pivot_cols = Int[]
    work_A = copy(A)
    work_rhs = copy(rhs)
    for i in 1:n_constraints
        # Find best pivot column (nonzero in current row, highest priority, not already used)
        best_col = 0
        best_pri = -1
        for col in 1:ncols
            col in pivot_cols && continue
            work_A[i, col] == 0 && continue
            if priority[col] > best_pri
                best_pri = priority[col]
                best_col = col
            end
        end
        best_col == 0 && error("Degenerate constraint matrix at row $i")
        push!(pivot_cols, best_col)

        # Scale row so pivot = 1
        piv_val = work_A[i, best_col]
        work_A[i, :] ./= piv_val
        work_rhs[i] /= piv_val

        # Eliminate this column from all other rows
        for r in 1:n_constraints
            r == i && continue
            if work_A[r, best_col] != 0
                factor = work_A[r, best_col]
                work_A[r, :] .-= factor .* work_A[i, :]
                work_rhs[r] -= factor * work_rhs[i]
            end
        end
    end
    return work_A, work_rhs, pivot_cols
end

"""
    _exprs_from_elimination(work_A, work_rhs, pivot_cols, nsteps)

Read off dependent-parameter expressions from the reduced row echelon form.
Returns `dep_exprs::Dict{Symbol, Expr}`.
"""
function _exprs_from_elimination(work_A, work_rhs, pivot_cols, nsteps)
    dep_exprs = Dict{Symbol, Expr}()
    for (i, pcol) in enumerate(pivot_cols)
        step_idx = (pcol + 1) ÷ 2
        is_forward = isodd(pcol)
        dep_sym = Symbol("k$(step_idx)", is_forward ? "f" : "r")

        keq_exp = work_rhs[i]
        factors = Tuple{Symbol, Rational{BigInt}}[]

        for c in 1:2*nsteps
            c == pcol && continue
            work_A[i, c] == 0 && continue
            coeff = work_A[i, c]
            s_idx = (c + 1) ÷ 2
            is_fwd = isodd(c)
            k_sym = Symbol("k$(s_idx)", is_fwd ? "f" : "r")
            push!(factors, (k_sym, -coeff))
        end

        dep_exprs[dep_sym] = _build_power_expr(keq_exp, factors)
    end
    return dep_exprs
end

# ─── Main entry point ────────────────────────────────────────────────────────

"""
    _dependent_param_exprs(M::Type{<:EnzymeMechanism})

Select dependent parameters and build substitution expressions.
Returns `(dep_exprs, indep_params)`:
- `dep_exprs`: Dict mapping dependent k Symbol to Expr using `params.independent_k` and `params.Keq`
- `indep_params`: tuple of independent parameter Symbols (not including Keq, E_total)
"""
function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
    C, rhs_coeffs = _thermodynamic_constraints(M)
    n_constraints = size(C, 1)
    nsteps = size(C, 2)

    # If no constraints, all params are independent
    if n_constraints == 0
        all_k = Symbol[]
        for i in 1:nsteps
            push!(all_k, Symbol("k$(i)f"))
            push!(all_k, Symbol("k$(i)r"))
        end
        return Dict{Symbol, Expr}(), Tuple(all_k)
    end

    A, rhs = _build_log_constraint_matrix(C, rhs_coeffs, nsteps)

    free_binding = Set(_free_enzyme_binding_steps(M))

    # Determine which steps have metabolites on LHS or RHS
    m = M()
    rxns = reactions(m)
    enzs = enzyme_forms(m)
    enz_set = Set(e[1] for e in enzs)
    step_has_met_lhs = Bool[]
    step_has_met_rhs = Bool[]
    for (lhs, rhs_side) in rxns
        push!(step_has_met_lhs, any(s ∉ enz_set for s in lhs))
        push!(step_has_met_rhs, any(s ∉ enz_set for s in rhs_side))
    end

    priority = _pivot_priority(nsteps, free_binding, step_has_met_lhs, step_has_met_rhs)
    work_A, work_rhs, pivot_cols = _pivoted_gaussian_elimination(A, rhs, priority)
    dep_exprs = _exprs_from_elimination(work_A, work_rhs, pivot_cols, nsteps)

    dep_set = Set(keys(dep_exprs))

    # Independent params: all k's not in dep_set
    indep = Symbol[]
    for j in 1:nsteps
        kf = Symbol("k$(j)f")
        kr = Symbol("k$(j)r")
        kf ∉ dep_set && push!(indep, kf)
        kr ∉ dep_set && push!(indep, kr)
    end

    return dep_exprs, Tuple(indep)
end

"""Build an Expr for `params.Keq^keq_exp * prod(params.k_i^exp_i)` with rational exponents."""
function _build_power_expr(keq_exp::Rational, factors::Vector{Tuple{Symbol, Rational{BigInt}}})
    terms = Expr[]

    if keq_exp != 0
        push!(terms, _power_term(:Keq, keq_exp))
    end

    for (sym, exp) in factors
        exp == 0 && continue
        push!(terms, _power_term(sym, exp))
    end

    isempty(terms) && return :(1)
    length(terms) == 1 && return terms[1]
    return Expr(:call, :*, terms...)
end

"""
    _power_term(sym, exp)

Build an Expr for `params.sym ^ exp`, simplifying for common integer exponents.
"""
function _power_term(sym::Symbol, exp::Rational)
    base = :(params.$sym)
    if exp == 1
        return base
    elseif exp == -1
        return :(1 / $base)
    elseif denominator(exp) == 1
        e = Int(exp)
        if e > 0
            return :($base ^ $e)
        else
            return :(1 / $base ^ $(-e))
        end
    else
        # General rational exponent
        return :($base ^ $(Float64(exp)))
    end
end
