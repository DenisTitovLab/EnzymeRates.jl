"""
Haldane and Wegscheider thermodynamic constraints for enzyme mechanisms.

Identifies thermodynamic cycles in the mechanism graph (via null-space of the
enzyme incidence matrix), classifies them as Haldane (net reaction) or
Wegscheider (internal loop), and performs Gaussian elimination to express a
minimal set of dependent rate constants in terms of the independent ones and
the equilibrium constant Keq.
"""

function _split_reaction_side(side, enz_set)
    enzyme_sym = first(s for s in side if s in enz_set)
    met_syms = Symbol[s for s in side if s ∉ enz_set]
    return enzyme_sym, met_syms
end

function _enzyme_incidence_matrix(M::Type{<:EnzymeMechanism})
    m = M()
    enzs = enzyme_forms(m)
    rxns = reactions(m)
    enz_names = Tuple(e[1] for e in enzs)
    enz_set = Set(enz_names)
    N = length(enzs)
    B = zeros(Int, N, length(rxns))
    for (j, (lhs, rhs)) in enumerate(rxns)
        e_lhs, _ = _split_reaction_side(lhs, enz_set)
        e_rhs, _ = _split_reaction_side(rhs, enz_set)
        B[findfirst(==(e_lhs), enz_names), j] -= 1
        B[findfirst(==(e_rhs), enz_names), j] += 1
    end
    return B
end

function _integer_nullspace(A::Matrix{Int})
    m, n = size(A)
    R = Matrix{Rational{BigInt}}(A)
    pivot_cols = Int[]
    row = 1
    for col in 1:n
        piv = 0
        for r in row:m
            R[r, col] != 0 && (piv = r; break)
        end
        piv == 0 && continue
        R[row, :], R[piv, :] = R[piv, :], R[row, :]
        R[row, :] ./= R[row, col]
        for r in 1:m
            r == row && continue
            R[r, col] != 0 && (R[r, :] .-= R[r, col] .* R[row, :])
        end
        push!(pivot_cols, col)
        row += 1
    end
    free_cols = setdiff(1:n, pivot_cols)
    nullity = length(free_cols)
    nullity == 0 && return zeros(Int, n, 0)
    NS = zeros(Rational{BigInt}, n, nullity)
    for (k, fc) in enumerate(free_cols)
        NS[fc, k] = 1
        for (r, pc) in enumerate(pivot_cols)
            NS[pc, k] = -R[r, fc]
        end
    end
    result = zeros(Int, n, nullity)
    for k in 1:nullity
        lcm_val = lcm(denominator.(NS[:, k])...)
        for i in 1:n
            result[i, k] = Int(NS[i, k] * lcm_val)
        end
        g = gcd(abs.(result[:, k])...)
        g > 0 && (result[:, k] .÷= g)
        first_nz = findfirst(!=(0), result[:, k])
        first_nz !== nothing && result[first_nz, k] < 0 && (result[:, k] .*= -1)
    end
    return result
end

function _thermodynamic_constraints(M::Type{<:EnzymeMechanism})
    B = _enzyme_incidence_matrix(M)
    NS = _integer_nullspace(B)
    n_constraints = size(NS, 2)
    n_constraints == 0 && return zeros(Int, 0, size(B, 2)), Int[]

    S = stoich_matrix(M())
    m = M()
    mets = metabolites(m)
    met_names = [mt[1] for mt in mets]
    nu_net = zeros(Int, length(met_names))
    for (name, _) in substrates(m)
        nu_net[findfirst(==(name), met_names)] -= 1
    end
    for (name, _) in products(m)
        nu_net[findfirst(==(name), met_names)] += 1
    end

    C = NS'
    rhs_coeffs = zeros(Int, n_constraints)
    for i in 1:n_constraints
        nu_cycle = S * C[i, :]
        if all(nu_cycle .== 0)
            rhs_coeffs[i] = 0
        else
            c = nothing
            valid = true
            for j in eachindex(nu_cycle)
                if nu_net[j] == 0
                    nu_cycle[j] != 0 && (valid = false; break)
                else
                    c_j = nu_cycle[j] // nu_net[j]
                    if c === nothing; c = c_j
                    elseif c_j != c; valid = false; break; end
                end
            end
            valid && c !== nothing && denominator(c) == 1 ? (rhs_coeffs[i] = Int(c)) :
                error("Cycle $i produces metabolite change not proportional to net reaction")
        end
    end
    return C, rhs_coeffs
end

"""
    _dependent_param_exprs(M::Type{<:EnzymeMechanism})

Select dependent parameters and build substitution expressions.
Returns `(dep_exprs, indep_params)`.
"""
function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
    C, rhs_coeffs = _thermodynamic_constraints(M)
    n_constraints = size(C, 1)
    nsteps = size(C, 2)
    m = M()
    eq_steps = equilibrium_steps(m)

    all_params = Symbol[]
    for i in 1:nsteps
        if eq_steps[i]; push!(all_params, Symbol("K$i"))
        else push!(all_params, Symbol("k$(i)f")); push!(all_params, Symbol("k$(i)r")); end
    end
    n_constraints == 0 && return Dict{Symbol, Expr}(), Tuple(all_params)

    # Build log-space constraint matrix: 1 col per RE step (K_i), 2 per SS step (k_jf, k_jr)
    n_vars = sum(is_re ? 1 : 2 for is_re in eq_steps)
    A = zeros(Rational{BigInt}, n_constraints, n_vars)
    rhs = Rational{BigInt}.(rhs_coeffs)
    step_cols = Vector{Any}(undef, nsteps)
    col = 1
    for (idx, is_re) in enumerate(eq_steps)
        if is_re; step_cols[idx] = col; col += 1
        else step_cols[idx] = (col, col + 1); col += 2; end
    end
    for i in 1:n_constraints, j in 1:nsteps
        C[i, j] == 0 && continue
        if eq_steps[j]; A[i, step_cols[j]] = C[i, j]
        else kf_col, kr_col = step_cols[j]; A[i, kf_col] = C[i, j]; A[i, kr_col] = -C[i, j]; end
    end

    # Pivot priority: internal isomerizations > metabolite steps > free-enzyme binding
    rxns = reactions(m)
    enzs = enzyme_forms(m)
    enz_set = Set(e[1] for e in enzs)
    free_enz_set = Set(e[1] for e in enzs if isempty(e[2]))
    priority = zeros(Int, n_vars)
    for j in 1:nsteps
        _, m_lhs = _split_reaction_side(rxns[j][1], enz_set)
        e_lhs_sym, _ = _split_reaction_side(rxns[j][1], enz_set)
        e_rhs_sym, m_rhs = _split_reaction_side(rxns[j][2], enz_set)
        has_lhs = !isempty(m_lhs); has_rhs = !isempty(m_rhs)
        is_free = e_lhs_sym in free_enz_set || e_rhs_sym in free_enz_set
        base = (has_lhs || has_rhs) ? (is_free ? 0 : 10) : 20
        if eq_steps[j]; priority[step_cols[j]] = base
        else kf_col, kr_col = step_cols[j]
            priority[kf_col] = base + (has_lhs ? 0 : 1)
            priority[kr_col] = base + (has_rhs ? 0 : 1)
        end
    end

    # Gaussian elimination with priority pivoting
    pivot_cols = Int[]
    work_A, work_rhs = copy(A), copy(rhs)
    for i in 1:n_constraints
        best_col, best_pri = 0, -1
        for c in 1:n_vars
            c in pivot_cols && continue; work_A[i, c] == 0 && continue
            priority[c] > best_pri && (best_pri = priority[c]; best_col = c)
        end
        best_col == 0 && error("Degenerate constraint matrix at row $i")
        push!(pivot_cols, best_col)
        pv = work_A[i, best_col]; work_A[i, :] ./= pv; work_rhs[i] /= pv
        for r in 1:n_constraints
            r == i && continue
            work_A[r, best_col] != 0 && (f = work_A[r, best_col]; work_A[r, :] .-= f .* work_A[i, :]; work_rhs[r] -= f * work_rhs[i])
        end
    end

    # Build column→symbol map and extract dependent expressions
    col_to_sym = Dict{Int, Symbol}()
    for j in 1:nsteps
        if eq_steps[j]; col_to_sym[step_cols[j]] = Symbol("K$j")
        else kf_col, kr_col = step_cols[j]; col_to_sym[kf_col] = Symbol("k$(j)f"); col_to_sym[kr_col] = Symbol("k$(j)r"); end
    end
    dep_exprs = Dict{Symbol, Expr}()
    for (i, pcol) in enumerate(pivot_cols)
        factors = Tuple{Symbol, Rational{BigInt}}[]
        for c in 1:n_vars
            c == pcol && continue; work_A[i, c] == 0 && continue
            push!(factors, (col_to_sym[c], -work_A[i, c]))
        end
        dep_exprs[col_to_sym[pcol]] = build_power_expr(work_rhs[i], factors)
    end

    dep_set = Set(keys(dep_exprs))
    indep = Symbol[p for p in all_params if p ∉ dep_set]
    return dep_exprs, Tuple(indep)
end

function _constraint_expr_strings(M::Type{<:EnzymeMechanism})
    dep_exprs, _ = _dependent_param_exprs(M)
    isempty(dep_exprs) && return String[]
    result = String[]
    for (dep_sym, expr) in sort(collect(dep_exprs); by=p -> string(p[1]))
        s = replace(string(expr), "params." => "")
        push!(result, "$dep_sym = $s")
    end
    result
end
