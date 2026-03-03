"""
Rate equation derivation and identifiability for OligomericEnzymeMechanism.

Implements the MWC (Monod-Wyman-Changeux) rate equation assembly for multi-site,
multi-conformation allosteric enzymes. Each catalytic subunit contributes
independently; regulatory sites are enzyme-level (n_reg < CatN) or per-subunit
(n_reg = CatN) as specified by the `RegSites` type parameter.

MWC rate formula (per conformation c, summed over conformations):
  num = CatN * sum_c( L_c * N_cat_c * Q_cat_c^(CatN-1)
            * prod(Q_reg_i_c^n_reg_i for i with n_reg_i == CatN) )
  den = sum_c( L_c * Q_cat_c^CatN * prod(Q_reg_i_c^n_reg_i) )
  v = E_total * num / den

Regulatory sites with n_reg_i < CatN (enzyme-level binding) appear only in
the denominator. Sites with n_reg_i == CatN (per-subunit binding) appear in
both numerator and denominator.
"""

# ─── Parameter naming ───────────────────────────────────────────

"""Append `_T` suffix to any K or k parameter symbol (not Keq, E_total)."""
function _rename_params_T(sym::Symbol)
    is_k_parameter(sym) ? Symbol(string(sym) * "_T") : sym
end

"""Name for a regulatory site parameter: K_{ligand}_reg{i} or K_{ligand}_T_reg{i}."""
function _reg_param_name(ligand::Symbol, site_idx::Int, T_state::Bool)
    T_state ? Symbol("K_$(ligand)_T_reg$(site_idx)") :
              Symbol("K_$(ligand)_reg$(site_idx)")
end

# ─── Accessors ──────────────────────────────────────────────────

"""Delegate structural accessors to the CatalyticMech singleton."""
n_states(::OligomericEnzymeMechanism{M,CM,N,RS,NC}) where {M,CM,N,RS,NC} =
    n_states(CM())
n_steps(::OligomericEnzymeMechanism{M,CM,N,RS,NC}) where {M,CM,N,RS,NC} =
    n_steps(CM())
equilibrium_steps(::OligomericEnzymeMechanism{M,CM,N,RS,NC}) where {M,CM,N,RS,NC} =
    equilibrium_steps(CM())
substrates(::OligomericEnzymeMechanism{M,CM,N,RS,NC}) where {M,CM,N,RS,NC} =
    substrates(CM())
products(::OligomericEnzymeMechanism{M,CM,N,RS,NC}) where {M,CM,N,RS,NC} =
    products(CM())
param_constraints(::OligomericEnzymeMechanism) = ()

"""Return all metabolite names (catalytic + regulatory) from the Metabolites type param."""
metabolites(::OligomericEnzymeMechanism{Mets,CM,N,RS,NC}) where {Mets,CM,N,RS,NC} =
    Tuple(m[1] for m in Mets)

# ─── Binding K symbols ──────────────────────────────────────────

"""
Return all binding (Kd-convention) K symbols: R-state, T-state, and reg site params.
Regulatory site K params are always Kd (dissociation constants).
"""
function _binding_K_symbols(
    ::Type{OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}},
) where {Mets,CM,CatN,RS,NConf}
    r_ks = Tuple(_binding_K_symbols(CM))
    t_ks = NConf == 2 ? Tuple(_rename_params_T(K) for K in r_ks) : ()
    reg_ks_r = Tuple(
        _reg_param_name(lig, i, false)
        for (i, (ligs, _)) in enumerate(RS) for lig in ligs
    )
    reg_ks_t = NConf == 2 ? Tuple(
        _reg_param_name(lig, i, true)
        for (i, (ligs, _)) in enumerate(RS) for lig in ligs
    ) : ()
    (r_ks..., t_ks..., reg_ks_r..., reg_ks_t...)
end

# ─── Dependent parameter expressions ────────────────────────────

"""
Collect all k/K symbols from an Expr/Symbol tree.
Used to build T-state renaming substitution.
"""
function _k_syms_in_expr(expr)
    syms = Set{Symbol}()
    function collect!(e)
        if e isa Symbol
            is_k_parameter(e) && push!(syms, e)
        elseif e isa Expr
            for a in e.args; collect!(a); end
        end
    end
    collect!(expr)
    syms
end

"""Rename all k/K symbols (excluding Keq) in an expression with _T suffix."""
function _rename_expr_params_T(expr)
    k_syms = _k_syms_in_expr(expr)
    subs = Dict(s => _rename_params_T(s) for s in k_syms if s != :Keq)
    isempty(subs) ? expr : substitute_params_expr(expr, subs)
end

"""
    _dependent_param_exprs(M::Type{<:OligomericEnzymeMechanism})

Return `(dep_exprs, indep_params)` for an OligomericEnzymeMechanism.

For NConf=1: delegates to CatalyticMech; adds reg site params to indep.
For NConf=2: duplicates R-state analysis with _T suffix for T-state;
             adds reg site params (R and T) and L to indep.
"""
function _dependent_param_exprs(
    ::Type{OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}},
) where {Mets,CM,CatN,RS,NConf}
    dep_R, indep_R = _dependent_param_exprs(CM)

    reg_params_r = [
        _reg_param_name(lig, i, false)
        for (i, (ligs, _)) in enumerate(RS) for lig in ligs
    ]

    if NConf == 1
        return dep_R, (indep_R..., reg_params_r...)
    end

    # NConf == 2: add T-state params
    dep_T = Dict{Symbol, Union{Symbol, Expr}}(
        _rename_params_T(k) => _rename_expr_params_T(v)
        for (k, v) in dep_R
    )
    indep_T = Tuple(_rename_params_T(p) for p in indep_R)

    reg_params_t = [
        _reg_param_name(lig, i, true)
        for (i, (ligs, _)) in enumerate(RS) for lig in ligs
    ]

    merged_dep = merge(dep_R, dep_T)
    merged_indep = (indep_R..., indep_T..., reg_params_r..., reg_params_t..., :L)
    return merged_dep, merged_indep
end

# ─── Parameters API ─────────────────────────────────────────────

parameters(m::OligomericEnzymeMechanism) = parameters(m, Reduced)

function parameters(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}, ::ReducedMode,
) where {Mets,CM,CatN,RS,NConf}
    M = OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

function fitted_params(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf},
) where {Mets,CM,CatN,RS,NConf}
    M = OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}
    _, indep = _dependent_param_exprs(M)
    indep
end

# ─── Rate body building helpers ──────────────────────────────────

"""Build the regulatory site partition function expression: 1 + lig/K_lig_reg_i + ..."""
function _reg_site_expr(ligs, site_idx::Int, T_state::Bool)
    terms = Any[1]
    for lig in ligs
        K_sym = _reg_param_name(lig, site_idx, T_state)
        push!(terms, :($(lig) / $K_sym))
    end
    _nest_binary(:+, terms)
end

"""Raise an expression to an integer power (returns 1 for n=0, expr for n=1)."""
function _power_expr(expr, n::Int)
    n == 0 && return 1
    n == 1 && return expr
    :(($expr)^$n)
end

"""Build the MWC rate equation body as an Expr block."""
function _build_oligomeric_rate_body(Mets, CM, CatN, RS, NConf)
    M_type = OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}

    # Get symbolic rate polynomials from CatalyticMech
    num_fs, denom_terms = _raw_symbolic_rate_polys(CM)

    # Build param/conc symbol sets for expression generation
    m_cat = CM()
    cat_constr = Set(c[1] for c in param_constraints(m_cat))
    cat_params = Set{Symbol}(
        p for p in _raw_param_symbols(equilibrium_steps(m_cat)) if p ∉ cat_constr
    )
    cat_mets = Set{Symbol}(metabolites(m_cat))
    binding_Ks_r = Set(_binding_K_symbols(CM))

    # R-state expressions
    N_R_expr = _factored_sigma_to_expr(num_fs, cat_params, cat_mets, binding_Ks_r)
    Q_R_expr = _denom_terms_to_expr(denom_terms, cat_params, cat_mets, binding_Ks_r)

    # Dep param assignments (R-state): k3r = inv(K2)*K1*k3f*inv(Keq), etc.
    dep_R, _ = _dependent_param_exprs(CM)
    dep_R_kd = _apply_kd_inversion(dep_R, CM, K -> :(inv($K)))
    r_assignments = Expr[
        Expr(:(=), sym, dep_R_kd[sym])
        for (sym, _) in sort(collect(dep_R_kd); by=first)
    ]

    # T-state expressions (NConf=2 only)
    N_T_expr, Q_T_expr = nothing, nothing
    t_assignments = Expr[]
    if NConf == 2
        T_subs = Dict(K => _rename_params_T(K) for K in cat_params if is_k_parameter(K))
        N_T_expr = substitute_params_expr(N_R_expr, T_subs)
        Q_T_expr = substitute_params_expr(Q_R_expr, T_subs)

        t_assignments = Expr[
            Expr(:(=), _rename_params_T(sym), _rename_expr_params_T(dep_R_kd[sym]))
            for (sym, _) in sort(collect(dep_R_kd); by=first)
        ]
    end

    # Regulatory site expressions
    reg_Q_R = Any[_reg_site_expr(ligs, i, false) for (i, (ligs, _)) in enumerate(RS)]
    reg_Q_T = NConf == 2 ?
        Any[_reg_site_expr(ligs, i, true) for (i, (ligs, _)) in enumerate(RS)] : Any[]

    # Assemble numerator and denominator
    #   Per-subunit reg sites (n_reg == CatN) go in BOTH num and den.
    #   Enzyme-level reg sites (n_reg < CatN) go in denominator only.
    function num_reg_factors(reg_Q)
        factors = Any[]
        for (idx, (_, n_reg)) in enumerate(RS)
            n_reg == CatN || continue
            push!(factors, _power_expr(reg_Q[idx], n_reg))
        end
        factors
    end

    function den_reg_factors(reg_Q)
        factors = Any[]
        for (idx, (_, n_reg)) in enumerate(RS)
            push!(factors, _power_expr(reg_Q[idx], n_reg))
        end
        factors
    end

    # Build per-conformation terms
    function num_term(N_expr, Q_expr, reg_Q)
        factors = Any[N_expr]
        CatN > 1 && push!(factors, _power_expr(Q_expr, CatN - 1))
        append!(factors, num_reg_factors(reg_Q))
        _nest_binary(:*, factors)
    end

    function den_term(Q_expr, reg_Q)
        factors = Any[_power_expr(Q_expr, CatN)]
        append!(factors, den_reg_factors(reg_Q))
        _nest_binary(:*, factors)
    end

    _num_term_R = num_term(N_R_expr, Q_R_expr, reg_Q_R)
    _den_term_R = den_term(Q_R_expr, reg_Q_R)

    full_num, full_den = if NConf == 1
        :($(CatN) * $(_num_term_R)),
        _den_term_R
    else
        _num_term_T = num_term(N_T_expr, Q_T_expr, reg_Q_T)
        _den_term_T = den_term(Q_T_expr, reg_Q_T)
        :($(CatN) * ($(_num_term_R) + L * $(_num_term_T))),
        :($(_den_term_R) + L * $(_den_term_T))
    end

    rate_expr = :(E_total * ($full_num) / ($full_den))

    # Build destructuring lines
    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq, :E_total)
    all_mets = Tuple(m[1] for m in Mets)

    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(all_mets, :concs),
        r_assignments...,
        t_assignments...,
        :(return $rate_expr))
end

# ─── Rate equation dispatch ──────────────────────────────────────

function rate_equation(m::OligomericEnzymeMechanism, concs, params)
    rate_equation(m, concs, params, Reduced)
end

@generated function rate_equation(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf},
    concs::NamedTuple, params::NamedTuple, ::ReducedMode,
) where {Mets,CM,CatN,RS,NConf}
    _build_oligomeric_rate_body(Mets, CM, CatN, RS, NConf)
end

# ─── String representation ───────────────────────────────────────

rate_equation_string(m::OligomericEnzymeMechanism) = rate_equation_string(m, Reduced)

function rate_equation_string(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}, ::ReducedMode,
) where {Mets,CM,CatN,RS,NConf}
    M = OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    all_mets = Tuple(m[1] for m in Mets)

    # Dep param assignment strings (R-state)
    dep_R, _ = _dependent_param_exprs(CM)
    dep_R_kd = _apply_kd_inversion(dep_R, CM, K -> :(1 / $K))
    dep_lines = [
        "$sym = $(_expr_to_string(dep_R_kd[sym]))"
        for (sym, _) in sort(collect(dep_R_kd); by=first)
    ]

    # T-state dep param assignments (NConf=2)
    t_dep_lines = if NConf == 2
        [
            "$(_rename_params_T(sym)) = $(_expr_to_string(_rename_expr_params_T(dep_R_kd[sym])))"
            for (sym, _) in sort(collect(dep_R_kd); by=first)
        ]
    else
        String[]
    end

    # Build the v = line using the same logic as _build_oligomeric_rate_body.
    # Extract the num and den sub-expressions from the generated body so we can
    # format them with explicit parentheses: "v = E_total * (...) / (...)".
    # (Without this, _expr_to_string formats "E_total * CatN * N / D" due to
    # left-to-right * association, which drops the opening paren the test expects.)
    body_expr = _build_oligomeric_rate_body(Mets, CM, CatN, RS, NConf)
    ret_expr = body_expr.args[end]       # :(return E_total * num / den)
    rate_v_expr = ret_expr.args[1]       # Expr(:call, :/, Expr(:call, :*, E_total, num), den)
    full_num_expr = rate_v_expr.args[2].args[3]   # num factor (after E_total in the *)
    full_den_expr = rate_v_expr.args[3]            # denominator

    v_line = "v = E_total * ($(_expr_to_string(full_num_expr))) / ($(_expr_to_string(full_den_expr)))"

    join([
        "(; $(join(hw_params, ", "))) = params",
        "(; $(join(all_mets, ", "))) = concs",
        dep_lines...,
        t_dep_lines...,
        v_line,
    ], "\n")
end

# ─── Identifiability ─────────────────────────────────────────────

"""
Raise a POLY to an integer power via repeated multiplication.
Returns poly_one() for n=0.
"""
function _poly_power(p::POLY, n::Int)
    n == 0 && return poly_one()
    result = poly_one()
    for _ in 1:n
        result = poly_mul(result, p)
    end
    result
end

"""Rename all K/k symbols (not Keq) in a POLY with _T suffix."""
function _rename_poly_T(p::POLY)
    POLY(
        sort!(MONO([
            (is_k_parameter(s) && s != :Keq ? _rename_params_T(s) : s) => e
            for (s, e) in mono
        ]); by=first) => coeff
        for (mono, coeff) in p
    )
end

"""
Count distinct concentration monomials in the full oligomeric rate numerator and
denominator. Treats all K/k/L symbols as parameters (strips them from monomials).
Returns `(n_num, n_denom)`.
"""
function _count_oligomeric_rate_monomials(CM, CatN, RS, NConf)
    num_fs, denom_terms = _raw_symbolic_rate_polys(CM)
    N_cat_R = _expand_factored_sigma(num_fs)
    Q_cat_R = _expand_to_poly(denom_terms)

    # Build reg site partition polynomials (1 + ligand_sym for each ligand)
    reg_Q_R = POLY[
        reduce(poly_add, (poly_add(poly_one(), poly_sym(lig)) for lig in ligs))
        for (ligs, _) in RS
    ]

    # Assemble for NConf conformations
    # NConf=1: single conformation (no L)
    # NConf=2: R-state + L * T-state

    function num_poly_for_conf(N_cat, Q_cat, reg_Qs, L_factor)
        n_term = poly_mul(N_cat, _poly_power(Q_cat, CatN - 1))
        for (idx, (_, n_reg)) in enumerate(RS)
            n_reg == CatN || continue
            n_term = poly_mul(n_term, _poly_power(reg_Qs[idx], n_reg))
        end
        L_factor === nothing ? n_term : poly_mul(poly_sym(L_factor), n_term)
    end

    function den_poly_for_conf(Q_cat, reg_Qs, L_factor)
        d_term = _poly_power(Q_cat, CatN)
        for (idx, (_, n_reg)) in enumerate(RS)
            d_term = poly_mul(d_term, _poly_power(reg_Qs[idx], n_reg))
        end
        L_factor === nothing ? d_term : poly_mul(poly_sym(L_factor), d_term)
    end

    full_num = num_poly_for_conf(N_cat_R, Q_cat_R, reg_Q_R, nothing)
    full_den = den_poly_for_conf(Q_cat_R, reg_Q_R, nothing)

    if NConf == 2
        N_cat_T = _rename_poly_T(N_cat_R)
        Q_cat_T = _rename_poly_T(Q_cat_R)
        reg_Q_T = POLY[_rename_poly_T(q) for q in reg_Q_R]

        full_num = poly_add(full_num, num_poly_for_conf(N_cat_T, Q_cat_T, reg_Q_T, :L))
        full_den = poly_add(full_den, den_poly_for_conf(Q_cat_T, reg_Q_T, :L))
    end

    # Count distinct concentration monomials (strip all k/K/L params)
    reg_ligands = Set{Symbol}(lig for (ligs, _) in RS for lig in ligs)
    cat_mets = Set{Symbol}(metabolites(CM()))

    function conc_mono(mono)
        sort!(MONO([
            s => e for (s, e) in mono
            if !is_k_parameter(s) && s != :E_total && s != :Keq && s != :L
        ]); by=first)
    end

    n_num = length(unique(conc_mono(k) for (k, _) in full_num))
    n_denom = length(unique(conc_mono(k) for (k, _) in full_den))
    n_num, n_denom
end

@generated function structural_identifiability_deficit(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf},
) where {Mets,CM,CatN,RS,NConf}
    M = OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}
    _, indep = _dependent_param_exprs(M)
    n_k = length(indep)
    n_num, n_denom = _count_oligomeric_rate_monomials(CM, CatN, RS, NConf)
    n_k - (n_num - 1) - (n_denom - 1)
end
