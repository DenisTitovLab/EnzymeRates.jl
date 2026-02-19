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

"""
    _dependent_param_exprs(M::Type{<:EnzymeMechanism})

Select dependent parameters and build substitution expressions.
Returns `(dep_exprs, indep_params)`.
"""
function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
    C, rhs_coeffs = _thermodynamic_constraints(M)
    nc = size(C, 1)
    nsteps = size(C, 2)
    m = M()
    eq_steps = equilibrium_steps(m)
    constrained = Set(c[1] for c in param_constraints(m))
    all_params = [p for p in _raw_param_symbols(eq_steps) if p ∉ constrained]
    nc == 0 && return Dict{Symbol, Expr}(), Tuple(all_params)

    # Constraint substitution map: target → factors (for log-space column merging)
    csub = Dict(target => factors for (target, _, factors) in param_constraints(m))

    # Build log-space constraint matrix directly for non-constrained params only.
    # Constrained params are substituted inline: if K3=K1, K3's column merges into K1's.
    sym_col = Dict(p => i for (i, p) in enumerate(all_params))
    n_vars = length(all_params)
    A = zeros(Rational{BigInt}, nc, n_vars)
    rhs = Rational{BigInt}.(rhs_coeffs)
    for i in 1:nc, j in 1:nsteps
        C[i, j] == 0 && continue
        pairs = if eq_steps[j]
            [(Symbol("K$j"), 1)]
        else
            [(Symbol("k$(j)f"), 1), (Symbol("k$(j)r"), -1)]
        end
        for (sym, sign) in pairs
            targets = if haskey(csub, sym)
                [(s, sign * e) for (s, e) in csub[sym]]
            else
                [(sym, sign)]
            end
            for (s, sgn) in targets
                haskey(sym_col, s) && (A[i, sym_col[s]] += C[i, j] * sgn)
            end
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
    pivot_cols = Int[]
    wA, wrhs = copy(A), copy(rhs)
    for i in 1:nc
        best_col, best_pri = 0, -1
        for c in 1:n_vars
            c in pivot_cols && continue
            wA[i, c] == 0 && continue
            if priority[c] > best_pri
                best_pri = priority[c]
                best_col = c
            end
        end
        best_col == 0 && error(
            "Degenerate constraint matrix at row $i"
        )
        push!(pivot_cols, best_col)
        pv = wA[i, best_col]
        wA[i, :] ./= pv
        wrhs[i] /= pv
        for r in 1:nc
            if r != i && wA[r, best_col] != 0
                f = wA[r, best_col]
                wA[r, :] .-= f .* wA[i, :]
                wrhs[r] -= f * wrhs[i]
            end
        end
    end

    dep_exprs = Dict{Symbol, Expr}()
    for (i, pcol) in enumerate(pivot_cols)
        factors = [
            (all_params[c], -wA[i, c])
            for c in 1:n_vars
            if c != pcol && wA[i, c] != 0
        ]
        dep_exprs[all_params[pcol]] = build_power_expr(wrhs[i], factors)
    end
    dep_set = Set(keys(dep_exprs))
    return dep_exprs, Tuple(p for p in all_params if p ∉ dep_set)
end

function _constraint_expr_strings(M::Type{<:EnzymeMechanism})
    m = M()
    lines = String[]
    # User constraints first
    for (target, coeff, factors) in param_constraints(m)
        push!(lines, _user_constraint_to_string(target, coeff, factors))
    end
    # Then HW constraints
    dep_exprs, _ = _dependent_param_exprs(M)
    if !isempty(dep_exprs)
        subs = Dict(K => :(1 / $K) for K in _binding_K_symbols(M))
        for (sym, expr) in sort(collect(dep_exprs); by=p -> string(p[1]))
            disp = isempty(subs) ? expr :
                substitute_params_expr(expr, subs)
            push!(lines, "$sym = $(string(disp))")
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

"""Collect sorted concentration symbols for a mechanism."""
function _sorted_conc_symbols(M::Type{<:EnzymeMechanism})
    metabolites(M())
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
    # Apply K→1/K to dep_exprs for Kd convention
    binding_Ks = Set(_binding_K_symbols(M))
    if !isempty(binding_Ks)
        inv_subs = Dict(K => :(inv($K)) for K in binding_Ks)
        dep_exprs = Dict(k => substitute_params_expr(v, inv_subs) for (k,v) in dep_exprs)
    end
    hw_params = (indep..., :Keq, :E_total)
    assignments = [Expr(:(=), sym, dep_exprs[sym])
                   for (sym, _) in sort(collect(dep_exprs); by=first)]
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(conc_syms, :concs),
        assignments...,
        expr)
end
