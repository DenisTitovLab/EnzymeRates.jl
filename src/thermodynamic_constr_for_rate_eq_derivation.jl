"""
Haldane and Wegscheider thermodynamic constraints for enzyme mechanisms,
plus preamble building helpers for @generated rate equation bodies.

Identifies thermodynamic cycles in the mechanism graph (via null-space of the
enzyme incidence matrix), classifies them as Haldane (net reaction) or
Wegscheider (internal loop), and performs Gaussian elimination to express a
minimal set of dependent rate constants in terms of the independent ones and
the equilibrium constant Keq.
"""

# ─── Shared Helpers ──────────────────────────────────────────────

"""Collect raw parameter symbols (K_i for RE, k_if/k_ir for SS) for each step."""
function _raw_param_symbols(eq_steps)
    ps = Symbol[]
    for (i, re) in enumerate(eq_steps)
        if re
            push!(ps, Symbol("K$i"))
        else
            push!(ps, Symbol("k$(i)f"))
            push!(ps, Symbol("k$(i)r"))
        end
    end
    ps
end

# ─── Thermodynamic Constraint Infrastructure ─────────────────────

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
        piv = findfirst(r -> R[r, col] != 0, row:m)
        piv === nothing && continue
        piv += row - 1
        R[row, :], R[piv, :] = R[piv, :], R[row, :]
        R[row, :] ./= R[row, col]
        for r in 1:m
            r != row && R[r, col] != 0 && (R[r, :] .-= R[r, col] .* R[row, :])
        end
        push!(pivot_cols, col); row += 1
    end
    free_cols = setdiff(1:n, pivot_cols)
    isempty(free_cols) && return zeros(Int, n, 0)
    NS = zeros(Rational{BigInt}, n, length(free_cols))
    for (k, fc) in enumerate(free_cols)
        NS[fc, k] = 1
        for (r, pc) in enumerate(pivot_cols); NS[pc, k] = -R[r, fc]; end
    end
    result = zeros(Int, n, length(free_cols))
    for k in axes(result, 2)
        col = @view NS[:, k]
        l = lcm(denominator.(col)...)
        result[:, k] .= Int.(col .* l)
        g = gcd(abs.(result[:, k])...)
        g > 0 && (result[:, k] .÷= g)
        fnz = findfirst(!=(0), result[:, k])
        fnz !== nothing && result[fnz, k] < 0 && (result[:, k] .*= -1)
    end
    result
end

function _thermodynamic_constraints(M::Type{<:EnzymeMechanism})
    B = _enzyme_incidence_matrix(M)
    NS = _integer_nullspace(B)
    nc = size(NS, 2)
    nc == 0 && return zeros(Int, 0, size(B, 2)), Int[]
    m = M()
    S = stoich_matrix(m)
    met_names = collect(metabolites(m))
    nu_net = zeros(Int, length(met_names))
    for (name, _) in substrates(m); nu_net[findfirst(==(name), met_names)] -= 1; end
    for (name, _) in products(m); nu_net[findfirst(==(name), met_names)] += 1; end
    C = NS'
    rhs_coeffs = [_classify_cycle(S * C[i, :], nu_net, i) for i in 1:nc]
    return C, rhs_coeffs
end

function _classify_cycle(nu_cycle, nu_net, i)
    all(nu_cycle .== 0) && return 0
    c = nothing
    for j in eachindex(nu_cycle)
        if nu_net[j] == 0
            nu_cycle[j] != 0 && error(
                "Cycle $i produces metabolite " *
                "change not proportional to " *
                "net reaction"
            )
        else
            c_j = nu_cycle[j] // nu_net[j]
            if c === nothing
                c = c_j
            elseif c_j != c
                error(
                    "Cycle $i produces metabolite " *
                    "change not proportional to " *
                    "net reaction"
                )
            end
        end
    end
    err = "Cycle $i produces metabolite change " *
          "not proportional to net reaction"
    c !== nothing && denominator(c) == 1 ? Int(c) : error(err)
end

"""Add log-space contribution for a param, merging constrained params into replacements.
When `param` is user-constrained (target = coeff * prod(src^exp)), its column is merged
into the replacement columns and `log(coeff)` is tracked as a constant on the RHS.
The constraint coefficient is stored as a Rational (Ka domain) to handle binding-to-binding
inversion: Kd constraint K3=4K1 becomes Ka constraint K3=K1/4, stored as coeff=1//4."""
function _add_log_contrib!(A, rhs_const, row, param, coeff, sym_col, csub)
    if haskey(csub, param)
        c, factors = csub[param]
        # c is Rational{BigInt} (Ka-domain coefficient).
        # log(c) = log(numerator) - log(denominator).
        # LHS constant: coeff * log(c). Moving to RHS: subtract coeff * log(c).
        p, q = Int(numerator(c)), Int(denominator(c))
        if p != 1
            rhs_const[row][p] = get(rhs_const[row], p, Rational{BigInt}(0)) - coeff
        end
        if q != 1
            rhs_const[row][q] = get(rhs_const[row], q, Rational{BigInt}(0)) + coeff
        end
        for (src, exp) in factors
            _add_log_contrib!(A, rhs_const, row, src, coeff * exp, sym_col, csub)
        end
    else
        A[row, sym_col[param]] += coeff
    end
end

"""Multiply an expression by the constant factor `prod(base^exp)` from `const_dict`."""
function _apply_const_factor(expr, const_dict)
    rational_part = Rational{BigInt}(1)
    symbolic_parts = Tuple{Int, Rational{BigInt}}[]
    for (base, exp) in const_dict
        exp == 0 && continue
        if denominator(exp) == 1
            rational_part *= Rational{BigInt}(base) ^ Int(exp)
        else
            push!(symbolic_parts, (base, exp))
        end
    end
    terms = Any[]
    if rational_part != 1
        if denominator(rational_part) == 1
            push!(terms, Int(rational_part))
        else
            n, d = Int(numerator(rational_part)), Int(denominator(rational_part))
            push!(terms, :($n // $d))
        end
    end
    for (base, exp) in symbolic_parts
        push!(terms, denominator(exp) == 1 ?
            :($base ^ $(Int(exp))) : :($base ^ $(Float64(exp))))
    end
    isempty(terms) && return expr
    const_expr = length(terms) == 1 ? terms[1] : Expr(:call, :*, terms...)
    expr == 1 && return const_expr
    Expr(:call, :*, const_expr, expr)
end

"""
    _dependent_param_exprs(M::Type{<:EnzymeMechanism})

Select dependent parameters and build substitution expressions.
Returns `(dep_exprs, indep_params)` where constrained params are excluded from both.

User constraints are merged into the Haldane/Wegscheider matrix via column merging,
with `log(coeff)` tracked as constant contributions through Gaussian elimination.
This preserves coupling between user constraints and thermodynamic constraints.
"""
function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
    C, rhs_coeffs = _thermodynamic_constraints(M)
    nc = size(C, 1)
    nsteps = size(C, 2)
    m = M()
    eq_steps = equilibrium_steps(m)
    all_raw = _raw_param_symbols(eq_steps)
    constraints = param_constraints(m)
    constrained_set = Set(c[1] for c in constraints)
    unconstrained = Symbol[p for p in all_raw if p ∉ constrained_set]
    nc == 0 && return Dict{Symbol, Union{Symbol, Expr}}(), Tuple(unconstrained)

    # Build log-space substitution from user constraints in Ka domain.
    # User constraints are in Kd domain; for binding-to-binding constraints,
    # the coefficient inverts: Kd K3=4K1 → Ka K3=(1/4)K1.
    # This matches _apply_param_constraints in sym_poly_for_rate_eq_derivation.jl.
    binding_Ks = Set(_binding_K_symbols(M))
    csub = Dict{Symbol, Tuple{Rational{BigInt}, Vector{Tuple{Symbol, Rational{BigInt}}}}}()
    for (target, coeff, factors) in constraints
        is_b2b = target ∈ binding_Ks && all(f -> f[1] ∈ binding_Ks, factors)
        ka_coeff = is_b2b ? Rational{BigInt}(1, coeff) : Rational{BigInt}(coeff)
        csub[target] = (ka_coeff, [(sym, Rational{BigInt}(exp)) for (sym, exp) in factors])
    end

    # Variable list excludes constrained params (merged into replacements)
    all_params = unconstrained
    sym_col = Dict(p => i for (i, p) in enumerate(all_params))
    n_vars = length(all_params)

    A = zeros(Rational{BigInt}, nc, n_vars)
    rhs = Rational{BigInt}.(rhs_coeffs)
    # Track constant contributions: rhs_const[i][coeff_val] = exponent
    # Effective RHS: rhs[i]*log(Keq) + sum(exp*log(val) for (val,exp) in rhs_const[i])
    rhs_const = [Dict{Int, Rational{BigInt}}() for _ in 1:nc]

    for i in 1:nc, j in 1:nsteps
        C[i, j] == 0 && continue
        if eq_steps[j]
            _add_log_contrib!(A, rhs_const, i, Symbol("K$j"),
                              Rational{BigInt}(C[i, j]), sym_col, csub)
        else
            _add_log_contrib!(A, rhs_const, i, Symbol("k$(j)f"),
                              Rational{BigInt}(C[i, j]), sym_col, csub)
            _add_log_contrib!(A, rhs_const, i, Symbol("k$(j)r"),
                              Rational{BigInt}(-C[i, j]), sym_col, csub)
        end
    end

    # Pivot priority: internal isomerizations > metabolite steps > free-enzyme binding
    rxns = reactions(m)
    enzs = enzyme_forms(m)
    enz_set = Set(e[1] for e in enzs)
    free_enz_set = Set(e[1] for e in enzs if isempty(e[2]))
    priority = zeros(Int, n_vars)
    for (j, (lhs, rhs_rxn)) in enumerate(rxns)
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs_rxn, enz_set)
        has_met = !isempty(m_lhs) || !isempty(m_rhs)
        is_free = e_lhs in free_enz_set || e_rhs in free_enz_set
        if !has_met
            base = 20
        elseif is_free
            base = 0
        else
            base = 10
        end
        if eq_steps[j]
            s = Symbol("K$j")
            haskey(sym_col, s) && (priority[sym_col[s]] = (is_free && has_met) ? -1 : base)
        else
            for (suffix, offset) in (("f", 0), ("r", 1))
                s = Symbol("k$(j)$suffix")
                haskey(sym_col, s) && (priority[sym_col[s]] = base + offset)
            end
        end
    end

    # Gaussian elimination with priority pivoting
    pivot_entries = Tuple{Int, Int}[]
    pivot_col_set = Set{Int}()
    wA, wrhs = copy(A), copy(rhs)
    for i in 1:nc
        best_col, best_pri = 0, -1
        for c in 1:n_vars
            c in pivot_col_set && continue
            wA[i, c] == 0 && continue
            if priority[c] > best_pri
                best_pri = priority[c]
                best_col = c
            end
        end
        if best_col == 0
            is_zero = wrhs[i] == 0 && all(e == 0 for (_, e) in rhs_const[i])
            is_zero && continue  # redundant constraint (0 = 0)
            error(
                "Thermodynamically contradictory mechanism: " *
                "constraint row $i reduces to " *
                "0 = $(wrhs[i]) * log(Keq)"
            )
        end
        push!(pivot_entries, (i, best_col))
        push!(pivot_col_set, best_col)
        pv = wA[i, best_col]
        wA[i, :] ./= pv
        wrhs[i] /= pv
        for (c, e) in rhs_const[i]
            rhs_const[i][c] = e / pv
        end
        for r in 1:nc
            if r != i && wA[r, best_col] != 0
                f = wA[r, best_col]
                wA[r, :] .-= f .* wA[i, :]
                wrhs[r] -= f * wrhs[i]
                for (c, e) in rhs_const[i]
                    rhs_const[r][c] = get(rhs_const[r], c, Rational{BigInt}(0)) - f * e
                end
            end
        end
    end

    dep_exprs = Dict{Symbol, Union{Symbol, Expr}}()
    for (prow, pcol) in pivot_entries
        factors = [
            (all_params[c], -wA[prow, c])
            for c in 1:n_vars
            if c != pcol && wA[prow, c] != 0
        ]
        base_expr = build_power_expr(wrhs[prow], factors)
        dep_exprs[all_params[pcol]] = _apply_const_factor(base_expr, rhs_const[prow])
    end
    dep_set = Set(keys(dep_exprs))
    return dep_exprs, Tuple(p for p in all_params if p ∉ dep_set)
end

"""Apply K→1/K inversion to Haldane dep_exprs.
When a dependent param is itself a binding K, its RHS is wrapped in `inv_fn`
to compensate for the implicit LHS inversion (Ka→Kd)."""
function _apply_kd_inversion(dep_exprs, M::Type{<:EnzymeMechanism}, inv_fn)
    binding_Ks = Set(_binding_K_symbols(M))
    isempty(binding_Ks) && return dep_exprs
    inv_subs = Dict(K => inv_fn(K) for K in binding_Ks)
    Dict(
        k => begin
            rhs = substitute_params_expr(v, inv_subs)
            k in binding_Ks ? inv_fn(rhs) : rhs
        end
        for (k, v) in dep_exprs
    )
end

function _constraint_expr_strings(M::Type{<:EnzymeMechanism})
    m = M()
    lines = String[]
    for (target, coeff, factors) in param_constraints(m)
        push!(lines, _user_constraint_to_string(target, coeff, factors))
    end
    dep_exprs, _ = _dependent_param_exprs(M)
    if !isempty(dep_exprs)
        dep_exprs = _apply_kd_inversion(dep_exprs, M, K -> :(1 / $K))
        for (sym, expr) in sort(collect(dep_exprs); by=p -> string(p[1]))
            push!(lines, "$sym = $(string(expr))")
        end
    end
    lines
end

# ─── Preamble Building Helpers ───────────────────────────────────

"""Build destructuring Expr: (; a, b, c) = source"""
function _destructuring_expr(syms, source::Symbol)
    Expr(:(=), Expr(:tuple, Expr(:parameters, syms...)), source)
end

"""
Collect sorted raw parameter symbols (k1f, k1r, K1, ...,
E_total) for a mechanism, excluding constrained params.
"""
function _sorted_raw_param_symbols(M::Type{<:EnzymeMechanism})
    m = M()
    constrained = Set(c[1] for c in param_constraints(m))
    raw = _raw_param_symbols(equilibrium_steps(m))
    Tuple(
        p for p in (raw..., :E_total) if p ∉ constrained
    )
end

"""Full mode: destructure all params + concs, then raw expr."""
function _build_rate_body(M, ::Type{FullMode})
    expr, all_params, conc_syms = _raw_rate_expr_and_symbols(M)
    Expr(:block,
        _destructuring_expr(all_params, :params),
        _destructuring_expr(conc_syms, :concs),
        expr)
end

"""Reduced mode: destructure indep params + concs, define dep params, then raw expr."""
function _build_rate_body(M, ::Type{ReducedMode})
    expr, _, conc_syms = _raw_rate_expr_and_symbols(M)
    dep_exprs, indep = _dependent_param_exprs(M)
    dep_exprs = _apply_kd_inversion(dep_exprs, M, K -> :(inv($K)))
    hw_params = (indep..., :Keq, :E_total)
    assignments = [Expr(:(=), sym, dep_exprs[sym])
                   for (sym, _) in sort(collect(dep_exprs); by=first)]
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(conc_syms, :concs),
        assignments...,
        expr)
end
