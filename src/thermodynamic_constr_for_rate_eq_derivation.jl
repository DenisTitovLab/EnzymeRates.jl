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
representative step of each kinetic group, in step-order. Routes
Symbol production through the `name(p::Parameter, m)` chokepoint via
`_enumerate_parameters_full`.
"""
function _raw_param_symbols(m::Mechanism)
    Symbol[name(p, m) for p in _enumerate_parameters_full(m)]
end

_raw_param_symbols(m::EnzymeMechanism) = _raw_param_symbols(Mechanism(m))

"""
For each step in `m` (in flat-iteration order), yield the Parameter
instances that govern that step: `[Kd|Kiso]` for an RE step, `[Kon|Kfor,
Koff|Krev]` for an SS step. Each Parameter is anchored on the original
step (not the rep), so `name(p, m)` renders to the rep's structural
Symbol via the value-context chokepoint, collapsing kinetic-group members.
"""
function _step_parameters(m::Mechanism)
    out = Vector{Vector{Parameter}}()
    for (s, _) in _flat_steps(m)
        params = is_equilibrium(s) ?
            Parameter[is_binding(s) ? Kd(s, :None) : Kiso(s, :None)] :
            Parameter[is_binding(s) ? Kon(s, :None)  : Kfor(s, :None),
                      is_binding(s) ? Koff(s, :None) : Krev(s, :None)]
        push!(out, params)
    end
    out
end

# ─── Structural primacy: free-enzyme set + step priority ─────────

"""
Set of enzyme-form names that are NOT the RHS of any canonical RE binding
step `F + met… ⇌ F_bound`. SS steps don't determine binding state (their
direction isn't canonicalized). Shared by the kinetic-group name
representative and the Haldane elimination pivot.
"""
function _free_enz_set(m::Union{Mechanism, AllostericMechanism})
    em = compile_mechanism(m)
    enz_names = enzyme_forms(em)
    enz_set = Set(enz_names)
    free_enz_set = Set{Symbol}(enz_names)
    for (lhs, rhs, is_eq, _) in reactions(em)
        is_eq || continue
        _, m_l = _split_reaction_side(lhs, enz_set)
        e_r, m_r = _split_reaction_side(rhs, enz_set)
        if !isempty(m_l) && isempty(m_r)
            delete!(free_enz_set, e_r)
        end
    end
    free_enz_set
end

"""
Structural primacy base score for a step (lower = more primary / less
eliminable). Free-enzyme RE binding (-1) < free-enzyme SS binding (0) <
non-free metabolite step (10) < internal isomerization (20). Shared by the
kinetic-group name representative (argmin) and the Haldane elimination pivot
(argmax, which adds a +0/+1 forward/reverse offset per rate constant).
"""
function _step_priority(s::Step, free_enz_set::Set{Symbol})
    has_met = is_binding(s)
    is_free = (name(from_species(s)) in free_enz_set) ||
              (name(to_species(s))   in free_enz_set)
    is_equilibrium(s) && has_met && is_free && return -1
    return !has_met ? 20 : is_free ? 0 : 10
end

"""
Total lexical tiebreak for two distinct steps in the same kinetic group:
species pair + bound metabolite + RE/SS flag.
"""
_step_lex_key(s::Step) =
    (String(name(from_species(s))), String(name(to_species(s))),
     String(bound_metabolite(s) === nothing ? "" : name(bound_metabolite(s))),
     is_equilibrium(s))

"""
Kinetic-group naming representative: the structurally-primary step
(`argmin _step_priority`), with a deterministic lexical tiebreak.
"""
_group_rep(group::Vector{Step}, free_enz_set::Set{Symbol}) =
    argmin(s -> (_step_priority(s, free_enz_set), _step_lex_key(s)), group)

# ─── Thermodynamic Constraint Infrastructure ─────────────────────

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
    m = M()
    enz_names = enzyme_forms(m)
    enz_set = Set(enz_names)
    rxns = reactions(m)
    stoich_mat = stoich_matrix(m)[metabolite_row_range(m), :]
    met_names = collect(metabolites(m))
    subs_species = substrates(m)
    prods_species = products(m)

    # Enzyme incidence matrix: rows = enzyme forms, cols = steps;
    # entries are −1 / +1 depending on which side an enzyme form
    # appears on.
    N = length(enz_names)
    B = zeros(Int, N, length(rxns))
    for (j, (lhs, rhs, _, _)) in enumerate(rxns)
        e_lhs, _ = _split_reaction_side(lhs, enz_set)
        e_rhs, _ = _split_reaction_side(rhs, enz_set)
        B[findfirst(==(e_lhs), enz_names), j] -= 1
        B[findfirst(==(e_rhs), enz_names), j] += 1
    end

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

    # Classify each null-space cycle as Haldane (proportional to the
    # net reaction → contributes log(Keq)) or Wegscheider (closed
    # cycle, zero net change). Errors on cycles that touch metabolites
    # but aren't proportional to the net reaction.
    function classify_cycle(nu_cycle, i)
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

    C = NS'
    rhs_coeffs = [classify_cycle(stoich_mat * C[i, :], i) for i in 1:nc]
    return C, rhs_coeffs
end

"""
    _dependent_param_exprs(M::Type{<:EnzymeMechanism}) → (dep_exprs, indep_params)

Select dependent parameters and build substitution expressions for the
Haldane / Wegscheider thermodynamic constraints. Steps in the same
kinetic group share parameters: their cycle-incidence columns are merged
into the representative step's column before Gaussian elimination, so
`dep_exprs` and `indep_params` are keyed only on representatives.

The wrapper calls `_build_wegscheider_rename_map(M)` to obtain the
rename map for absorbed single-symbol Wegscheider RE ties and forwards
to `_dependent_param_exprs_kernel`.
"""
function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
    rename = _build_wegscheider_rename_map(M)
    dep_exprs, indep = _dependent_param_exprs_kernel(M, rename)
    # Filter Pass-2-absorbed symbols out of indep. Pass 2 of
    # `_build_wegscheider_rename_map` adds entries like `K_P_E => K_S_E`
    # when a Wegscheider tie collapses two binding-K group reps to the
    # same name. After the merge, the absorbed symbol doesn't appear in
    # the v polynomial — its column has been folded into the target.
    # But the absorbed symbol is still a kinetic-group rep in the
    # mechanism, so `_raw_param_symbols` emits it and the kernel keeps it
    # in `indep`. Without this filter, `fitted_params` exposes a fittable
    # dummy dimension that doesn't affect the loss, and finite-restart
    # convergence suffers (the same rate equation can land at noticeably
    # different fitted losses depending on which absorbed symbol got
    # the dummy slot).
    indep = Tuple(p for p in indep if get(rename, p, p) == p)
    return dep_exprs, indep
end

"""
Gaussian-elimination kernel underlying `_dependent_param_exprs`. Takes
the rename map as a parameter so callers can supply either the
user-defined kinetic-group rename (Pass 1 only) or the full rename map
that also absorbs single-symbol Wegscheider RE ties (Pass 1 + Pass 2).

Pass 2 of `_build_wegscheider_rename_map` calls this kernel with the
Pass-1-only rename to discover which single-symbol ties to absorb; the
display path in `rate_equation_string` likewise calls it with the
Pass-1-only rename to keep absorbed ties visible under the
`# Wegscheider constraints:` section.
"""
function _dependent_param_exprs_kernel(
    M::Type{<:EnzymeMechanism},
    rename::AbstractDict{Symbol, Symbol},
)
    m = M()
    mech = Mechanism(m)
    rxns = reactions(m)
    eq_steps = equilibrium_steps(m)
    enz_names = enzyme_forms(m)
    enz_set = Set(enz_names)

    free_enz_set = _free_enz_set(mech)

    C, rhs_coeffs = _thermodynamic_constraints(M)
    all_params = _raw_param_symbols(mech)
    nc = size(C, 1)
    nsteps = size(C, 2)
    nc == 0 && return (Dict{Symbol, Union{Symbol, Expr}}(),
                       Tuple(all_params))

    sym_col = Dict(p => i for (i, p) in enumerate(all_params))
    n_vars = length(all_params)

    # For each step's source-index, the Parameter(s) governing it.
    # `name(p, mech)` renders to the rep-renamed Symbol (Pass-1
    # kinetic-group rename is folded into the chokepoint); `rename` then
    # applies any Pass-2 single-symbol Wegscheider ties on top.
    step_params = _step_parameters(mech)
    step_name(p::Parameter) = get(rename, name(p, mech), name(p, mech))

    # Translate cycle-incidence columns into the merged-parameter A matrix.
    # Non-representative steps' columns are folded into their representative
    # via the chokepoint; this is mathematically equivalent to a kinetic-group
    # equality constraint (K_idx = K_rep, k_idx_f = k_rep_f, ...).
    #
    # Binding K's are Kd in the polynomial; cycle products use 1/Kd, so
    # binding-K column entries get a sign flip on top of the cycle incidence.
    binding_K_set = Set{Symbol}()
    for (j, (lhs, _, _, _)) in enumerate(rxns)
        eq_steps[j] || continue
        any(s ∉ enz_set for s in lhs) || continue
        push!(binding_K_set, step_name(step_params[j][1]))
    end

    A = zeros(Rational{BigInt}, nc, n_vars)
    rhs = Rational{BigInt}.(rhs_coeffs)
    for i in 1:nc, j in 1:nsteps
        C[i, j] == 0 && continue
        if eq_steps[j]
            sym = step_name(step_params[j][1])
            sign_factor = sym in binding_K_set ? -1 : 1
            A[i, sym_col[sym]] += sign_factor * C[i, j]
        else
            kf = step_name(step_params[j][1])
            kr = step_name(step_params[j][2])
            A[i, sym_col[kf]] += C[i, j]
            A[i, sym_col[kr]] -= C[i, j]
        end
    end

    # Pivot priority: internal isomerizations > metabolite steps
    #                 > free-enzyme binding (argmax picks most eliminable).
    #                 Same `_step_priority` scoring feeds the group-naming rep
    #                 (argmin). SS steps add a +0/+1 forward/reverse offset.
    priority = zeros(Int, n_vars)
    for j in 1:nsteps
        step = step_params[j][1].step
        base = _step_priority(step, free_enz_set)
        if eq_steps[j]
            s = step_name(step_params[j][1])
            haskey(sym_col, s) && (priority[sym_col[s]] = base)
        else
            for (offset, p) in enumerate(step_params[j])
                s = step_name(p)
                haskey(sym_col, s) &&
                    (priority[sym_col[s]] = base + offset - 1)
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
    hw_params = (indep..., :Keq, :E_total)
    assignments = [Expr(:(=), sym, dep_exprs[sym])
                   for (sym, _) in sort(collect(dep_exprs); by=first)]
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(conc_syms, :concs),
        assignments...,
        expr)
end
