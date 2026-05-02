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

"""
Collect raw parameter symbols (K_i for RE, k_if/k_ir for SS) for the
representative step of each kinetic group, in step-order.
"""
function _raw_param_symbols(m::EnzymeMechanism)
    eq = equilibrium_steps(m)
    ps = Symbol[]
    for g in kinetic_groups(m)
        rep = first(steps_in_group(m, g))
        if eq[rep]
            push!(ps, Symbol("K$rep"))
        else
            push!(ps, Symbol("k$(rep)f"))
            push!(ps, Symbol("k$(rep)r"))
        end
    end
    ps
end

# ─── Thermodynamic Constraint Infrastructure ─────────────────────

function _enzyme_incidence_matrix(enz_names, enz_set, rxns)
    N = length(enz_names)
    B = zeros(Int, N, length(rxns))
    for (j, (lhs, rhs, _, _)) in enumerate(rxns)
        e_lhs, _ = _split_reaction_side(lhs, enz_set)
        e_rhs, _ = _split_reaction_side(rhs, enz_set)
        B[findfirst(==(e_lhs), enz_names), j] -= 1
        B[findfirst(==(e_rhs), enz_names), j] += 1
    end
    return B
end

function _enzyme_incidence_matrix(M::Type{<:EnzymeMechanism})
    m = M()
    enz_names = enzyme_forms(m)
    _enzyme_incidence_matrix(
        enz_names, Set(enz_names), reactions(m),
    )
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

function _thermodynamic_constraints(
    enz_names, enz_set, rxns, stoich_mat,
    met_names, subs_species, prods_species,
)
    B = _enzyme_incidence_matrix(enz_names, enz_set, rxns)
    NS = _integer_nullspace(B)
    nc = size(NS, 2)
    nc == 0 && return zeros(Int, 0, size(B, 2)), Int[]
    nu_net = zeros(Int, length(met_names))
    for name in subs_species
        nu_net[findfirst(==(name), met_names)] -= 1
    end
    for name in prods_species
        nu_net[findfirst(==(name), met_names)] += 1
    end
    C = NS'
    rhs_coeffs = [
        _classify_cycle(stoich_mat * C[i, :], nu_net, i)
        for i in 1:nc
    ]
    return C, rhs_coeffs
end

function _thermodynamic_constraints(M::Type{<:EnzymeMechanism})
    m = M()
    enz_names = enzyme_forms(m)
    _thermodynamic_constraints(
        enz_names, Set(enz_names), reactions(m),
        stoich_matrix(m)[metabolite_row_range(m), :], collect(metabolites(m)),
        substrates(m), products(m),
    )
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
    _dependent_param_exprs(M::Type{<:EnzymeMechanism}) → (dep_exprs, indep_params)

Select dependent parameters and build substitution expressions for the
Haldane / Wegscheider thermodynamic constraints. Steps in the same
kinetic group share parameters: their cycle-incidence columns are merged
into the representative step's column before Gaussian elimination, so
`dep_exprs` and `indep_params` are keyed only on representatives.
"""
function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
    m = M()
    rxns = reactions(m)
    eq_steps = equilibrium_steps(m)
    enz_names = enzyme_forms(m)
    enz_set = Set(enz_names)

    # Free enzyme is any form that's never the RHS of a canonical RE binding
    # step `[F, met...] ⇌ [F_bound]`. SS steps don't determine binding state
    # (their direction isn't canonicalized). Used for pivot priority.
    free_enz_set = Set{Symbol}(enz_names)
    for (lhs, rhs, is_eq, _) in rxns
        is_eq || continue
        _, m_l = _split_reaction_side(lhs, enz_set)
        e_r, m_r = _split_reaction_side(rhs, enz_set)
        if !isempty(m_l) && isempty(m_r)
            delete!(free_enz_set, e_r)
        end
    end

    C, rhs_coeffs = _thermodynamic_constraints(M)
    all_params = _raw_param_symbols(m)
    nc = size(C, 1)
    nsteps = size(C, 2)
    nc == 0 && return (Dict{Symbol, Union{Symbol, Expr}}(),
                       Tuple(all_params))

    sym_col = Dict(p => i for (i, p) in enumerate(all_params))
    n_vars = length(all_params)

    # Build per-step → representative-step alias map
    rename = _build_kinetic_rename_map(m)

    # Translate cycle-incidence columns into the merged-parameter A matrix.
    # Non-representative steps' columns are folded into their representative
    # via the rename map; this is mathematically equivalent to a kinetic-group
    # equality constraint (K_idx = K_rep, k_idx_f = k_rep_f, ...).
    A = zeros(Rational{BigInt}, nc, n_vars)
    rhs = Rational{BigInt}.(rhs_coeffs)
    for i in 1:nc, j in 1:nsteps
        C[i, j] == 0 && continue
        if eq_steps[j]
            sym = Symbol("K$j")
            sym = get(rename, sym, sym)
            A[i, sym_col[sym]] += C[i, j]
        else
            kf = Symbol("k$(j)f"); kr = Symbol("k$(j)r")
            kf = get(rename, kf, kf); kr = get(rename, kr, kr)
            A[i, sym_col[kf]] += C[i, j]
            A[i, sym_col[kr]] -= C[i, j]
        end
    end

    # Pivot priority: internal isomerizations > metabolite steps
    #                 > free-enzyme binding. Inherits from the representative
    #                 step (first in the kinetic group).
    priority = zeros(Int, n_vars)
    for (j, (lhs, rhs_rxn, _, _)) in enumerate(rxns)
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs_rxn, enz_set)
        has_met = !isempty(m_lhs) || !isempty(m_rhs)
        is_free = (e_lhs in free_enz_set ||
                   e_rhs in free_enz_set)
        if !has_met
            base = 20
        elseif is_free
            base = 0
        else
            base = 10
        end
        if eq_steps[j]
            s = Symbol("K$j")
            haskey(sym_col, s) && (priority[sym_col[s]] =
                (is_free && has_met) ? -1 : base)
        else
            for (suffix, offset) in (("f", 0), ("r", 1))
                s = Symbol("k$(j)$suffix")
                haskey(sym_col, s) &&
                    (priority[sym_col[s]] = base + offset)
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
            wrhs[i] == 0 && continue  # redundant constraint (0 = 0)
            error(
                "Thermodynamically contradictory mechanism: " *
                "constraint row $i reduces to " *
                "0 = $(wrhs[i]) * log(Keq)")
        end
        push!(pivot_entries, (i, best_col))
        push!(pivot_col_set, best_col)
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

    dep_exprs = Dict{Symbol, Union{Symbol, Expr}}()
    for (prow, pcol) in pivot_entries
        factors = [
            (all_params[c], -wA[prow, c])
            for c in 1:n_vars
            if c != pcol && wA[prow, c] != 0
        ]
        dep_exprs[all_params[pcol]] = build_power_expr(wrhs[prow], factors)
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
    lines = String[]
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
Collect raw parameter symbols (one K or k_f/k_r per kinetic group) plus
`E_total`, in step order.
"""
function _sorted_raw_param_symbols(M::Type{<:EnzymeMechanism})
    Tuple((_raw_param_symbols(M())..., :E_total))
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
