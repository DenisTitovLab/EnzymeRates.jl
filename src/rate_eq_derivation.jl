# ─── Parameters API ─────────────────────────────────────────

const _AnyMechanism = Union{EnzymeMechanism, OligomericEnzymeMechanism}

"""
    parameters(m::EnzymeMechanism, [mode])
    parameters(m::OligomericEnzymeMechanism, [mode])

Return the parameter names required for the given mode as a tuple of Symbols.

# Modes
- `Reduced` (default): independent k's + Keq + E_total
- `Full`: all 2N k's + E_total (EnzymeMechanism only)
"""
function parameters end

parameters(m::_AnyMechanism) = parameters(m, Reduced)

@generated function parameters(
    ::EnzymeMechanism{Species, Reactions, EqSteps, PC},
    ::FullMode,
) where {Species, Reactions, EqSteps, PC}
    constrained = Set(c[1] for c in PC)
    Tuple(p for p in (_raw_param_symbols(EqSteps)..., :E_total) if p ∉ constrained)
end

@generated function parameters(::M, ::ReducedMode) where {M <: _AnyMechanism}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

"""Independent rate constant names for fitting (excludes Keq, E_total)."""
@generated function fitted_params(::M) where {M <: _AnyMechanism}
    _, indep = _dependent_param_exprs(M)
    indep
end

# ─── RE Group Helpers ───────────────────────────────────────

function _split_reaction_side(side, enz_set)
    enzyme_sym = first(s for s in side if s in enz_set)
    met_syms = Symbol[s for s in side if s ∉ enz_set]
    return enzyme_sym, met_syms
end

"""Compute RE-connected groups via union-find. Returns `(groups, form_to_group)`."""
function _compute_re_groups(enz_names, enz_set, rxns, eq_steps)
    N = length(enz_names)
    parent = collect(1:N)
    function find(x)
        while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end
        x
    end
    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] || continue
        e_lhs, _ = _split_reaction_side(lhs, enz_set)
        e_rhs, _ = _split_reaction_side(rhs, enz_set)
        ra, rb = find(findfirst(==(e_lhs), enz_names)), find(findfirst(==(e_rhs), enz_names))
        ra != rb && (parent[ra] = rb)
    end
    root_to_group = Dict{Int, Int}()
    groups = Vector{Vector{Int}}()
    form_to_group = zeros(Int, N)
    for i in 1:N
        r = find(i)
        g = get!(root_to_group, r) do; push!(groups, Int[]); length(groups) end
        push!(groups[g], i); form_to_group[i] = g
    end
    groups, form_to_group
end

"""
Compute alpha factors (relative concentrations within RE groups) as POLY values.
Returns `(alpha, sigma)` where alpha[i] is a (num::POLY, den::POLY) pair
and sigma[g] is a (num::POLY, den::POLY) pair.
"""
function _compute_alpha(enz_names, enz_set, rxns, eq_steps, groups)
    N = length(enz_names)
    mp = mets -> isempty(mets) ? poly_one() : reduce(poly_mul, poly_sym.(mets))
    alpha_num = Vector{POLY}(fill(poly_one(), N))
    alpha_den = Vector{POLY}(fill(poly_one(), N))

    for group in groups
        length(group) == 1 && continue
        visited = Set{Int}([group[1]])
        queue = [group[1]]
        while !isempty(queue)
            cur = popfirst!(queue)
            for (idx, (lhs, rhs)) in enumerate(rxns)
                eq_steps[idx] || continue
                e_l, m_l = _split_reaction_side(lhs, enz_set)
                e_r, m_r = _split_reaction_side(rhs, enz_set)
                i_f = findfirst(==(e_l), enz_names)
                j_f = findfirst(==(e_r), enz_names)
                K = poly_sym(Symbol("K$idx"))
                if i_f == cur && j_f ∉ visited
                    alpha_num[j_f] = poly_mul(poly_mul(alpha_num[cur], K), mp(m_l))
                    alpha_den[j_f] = poly_mul(alpha_den[cur], mp(m_r))
                    push!(visited, j_f); push!(queue, j_f)
                elseif j_f == cur && i_f ∉ visited
                    alpha_num[i_f] = poly_mul(alpha_num[cur], mp(m_r))
                    alpha_den[i_f] = poly_mul(poly_mul(alpha_den[cur], K), mp(m_l))
                    push!(visited, i_f); push!(queue, i_f)
                end
            end
        end
    end

    # Compute sigma per group: sum of alpha_i with cleared denominators
    sigma_num = Vector{POLY}(undef, length(groups))
    sigma_den = Vector{POLY}(undef, length(groups))
    for (g, group) in enumerate(groups)
        if length(group) == 1
            sigma_num[g] = sigma_den[g] = poly_one()
        else
            sigma_den[g] = reduce(poly_mul, alpha_den[i] for i in group)
            sigma_num[g] = reduce(poly_add,
                poly_mul(alpha_num[i],
                    reduce(poly_mul, (alpha_den[j] for j in group if j != i);
                        init=poly_one()))
                for i in group)
        end
    end
    alpha_num, alpha_den, sigma_num, sigma_den
end

# ─── Algebraic Regulator Factoring ─────────────────────────────

"""
Try algebraic polynomial factoring on a constrained sigma.
Tries each metabolite present in sigma: splits terms by metabolite presence,
divides the metabolite-containing terms by K*met, and checks whether
the result cleanly divides sigma into `base * (... + K*met)`.
Returns a `FactoredSigma` or `nothing`.
"""
function _try_algebraic_factor_sigma(
    sigma::POLY, rxns, eq_steps, enz_set, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
)
    # Classify symbols present in sigma
    all_syms = Set(s for (mono, _) in sigma for (s, _) in mono)
    K_syms = Set{Symbol}(
        s for s in all_syms
        if startswith(string(s), "K") && length(string(s)) > 1 &&
           isdigit(string(s)[2])
    )
    met_syms = sort!(Symbol[
        s for s in all_syms
        if s ∉ K_syms && !startswith(string(s), "k") &&
           s != :Keq && s != :E_total
    ]; by=string)

    # 1:1 constraint resolution map for finding canonical K symbols
    resolve = Dict{Symbol, Symbol}(
        target => factors[1][1]
        for (target, coeff, factors) in constraints
        if coeff == 1 && length(factors) == 1 && factors[1][2] == 1
    )

    # Try each metabolite, pick best factoring by score
    best = nothing
    best_score = (length(sigma), 0, 0)
    for met in met_syms
        # Find canonical binding K for this metabolite
        K_R = nothing
        for (idx, (lhs, rhs_rxn)) in enumerate(rxns)
            eq_steps[idx] || continue
            _, m_lhs = _split_reaction_side(lhs, enz_set)
            _, m_rhs = _split_reaction_side(rhs_rxn, enz_set)
            ((met ∈ m_lhs) ⊻ (met ∈ m_rhs)) || continue
            K_R = Symbol("K$idx")
            seen = Set{Symbol}()
            while haskey(resolve, K_R) && K_R ∉ seen
                push!(seen, K_R); K_R = resolve[K_R]
            end
            break
        end
        K_R === nothing && continue
        K_R ∈ K_syms || continue

        # Split sigma into terms with/without this metabolite
        without_R, with_R = POLY(), POLY()
        for (mono, coeff) in sigma
            (any(s == met for (s, _) in mono) ? with_R : without_R)[mono] = coeff
        end
        isempty(with_R) && continue

        # Divide each with_R term by K_R * met
        K_R_met = _mono(K_R => 1, met => 1)
        base_R = POLY()
        valid = true
        for (mono, coeff) in with_R
            r_exp = 0; k_exp = 0
            for (s, e) in mono
                s == met && (r_exp = e)
                s == K_R && (k_exp = e)
            end
            if r_exp < 1 || k_exp < 1; valid = false; break; end
            base_R[_mono_div(mono, K_R_met)] = coeff
        end
        (!valid || base_R == poly_one()) && continue

        # Check if base_R ⊆ without_R (subset means remainder factoring possible)
        is_subset = all(
            haskey(without_R, m) && without_R[m] == c for (m, c) in base_R
        )

        coeffs = POLY[]
        products = FactoredPoly[]
        remainder = is_subset ? poly_sub(without_R, base_R) : poly_zero()
        unfactored_count = is_subset ? length(remainder) : length(without_R)
        if is_subset
            one_plus_KR = poly_add(
                poly_one(), POLY(_mono(K_R => 1, met => 1) => 1),
            )
            # Try to absorb remainder: if remainder = q * base_R,
            # then sigma = base_R * (q + 1 + K*met)
            combined = one_plus_KR
            absorbed = isempty(remainder)
            if !absorbed
                q = _try_poly_exact_div(remainder, base_R)
                if q !== nothing
                    combined = poly_add(one_plus_KR, q)
                    absorbed = true
                end
            end
            if absorbed
                # Recurse on base_R for multi-site factoring
                inner = _try_algebraic_factor_sigma(
                    base_R, rxns, eq_steps, enz_set, constraints;
                    binding_Ks,
                )
                factors, exps = if inner !== nothing &&
                       length(inner.coefficients) == 1 &&
                       inner.coefficients[1] == poly_one()
                    fp = inner.products[1]
                    [fp.factors; combined], [fp.exponents; 1]
                else
                    [base_R, combined], [1, 1]
                end
                push!(coeffs, poly_one())
                push!(products, FactoredPoly(factors, exps))
                unfactored_count = 0
                remainder = poly_zero()
            else
                push!(coeffs, remainder)
                push!(products, FactoredPoly([poly_one()], [1]))
                push!(coeffs, poly_one())
                push!(products, FactoredPoly([base_R, one_plus_KR], [1, 1]))
            end
        else
            K_R_met_poly = POLY(_mono(K_R => 1, met => 1) => Rational{Int}(1))
            if !isempty(without_R)
                push!(coeffs, without_R)
                push!(products, FactoredPoly([poly_one()], [1]))
            end
            push!(coeffs, K_R_met_poly)
            push!(products, FactoredPoly([base_R], [1]))
        end
        # Monomial shape overlap tiebreaker (ignoring rate constants)
        _shape(m) = sort!([s => e for (s, e) in m if !startswith(string(s), "k")])
        uf_poly = is_subset ? remainder : without_R
        shapes_u = Set(_shape(m) for (m, _) in uf_poly)
        shapes_b = Set(_shape(m) for (m, _) in base_R)
        overlap = length(shapes_u ∩ shapes_b)
        score = (unfactored_count, -length(base_R), -overlap)
        if score < best_score
            best = FactoredSigma(coeffs, products)
            best_score = score
        end
    end
    best
end

"""
Try to express polynomial `p` as `Q^n` for some polynomial Q and integer n >= 2.
Returns `FactoredSigma` with a single term `FactoredPoly([Q], [n])`, or `nothing`.
"""
function _try_poly_power(p::POLY)
    length(p) < 3 && return nothing
    get(p, MONO(), 0) == 1 || return nothing
    max_e = maximum((e for (mono, _) in p for (_, e) in mono), init=0)
    max_e < 2 && return nothing
    for n in max_e:-1:2
        root_terms = [mono for (mono, c) in p
                      if !isempty(mono) && c == 1 &&
                         all(e % n == 0 for (_, e) in mono)]
        isempty(root_terms) && continue
        Q = POLY(MONO() => Rational{Int}(1))
        for mono in root_terms
            Q[MONO([s => div(e, n) for (s, e) in mono])] = 1
        end
        if _poly_power(Q, n) == p
            return FactoredSigma([poly_one()], [FactoredPoly([Q], [n])])
        end
    end
    nothing
end

"""Build rate poly for one SS step direction, clearing alpha denominators within group."""
function _ss_contrib(k_poly, mets, i_form, alpha_num, alpha_den, group)
    r = isempty(mets) ? k_poly : poly_mul(k_poly, reduce(poly_mul, poly_sym.(mets)))
    r = poly_mul(r, alpha_num[i_form])
    for k in group
        k == i_form && continue
        r = poly_mul(r, alpha_den[k])
    end
    r
end

"""
Build substitution pairs merging Haldane-derived parameters that have
identical expressions after user constraints. E.g., if k8r and k7r both
resolve to `k7f / (K1 * Keq)`, returns `[:k8r => :k7r]`.
"""
function _haldane_equality_substitutions(M::Type{<:EnzymeMechanism})
    dep_exprs, _ = _dependent_param_exprs(M)
    length(dep_exprs) < 2 && return Pair{Symbol,Symbol}[]
    # Sort by step number for stable canonical choice (lowest first)
    _step_num(s) = let m = match(r"\d+", string(s))
        m === nothing ? 0 : something(tryparse(Int, m.match::SubString), 0)
    end
    sorted = sort(collect(dep_exprs); by=p -> _step_num(p[1]))
    subs = Pair{Symbol,Symbol}[]
    expr_to_canonical = Dict{Any,Symbol}()
    for (sym, expr) in sorted
        canon = get!(expr_to_canonical, expr, sym)
        canon !== sym && push!(subs, sym => canon)
    end
    subs
end

# ─── Raw Rate Equation Derivation (Unified Cha / King-Altman) ───

"""
Try to factor a polynomial: power detection → algebraic factoring → trivial wrap.
When `check_benefit` is true, algebraic factoring is only used if it reduces
the display term count compared to the unfactored polynomial.
"""
function _factor_poly(
    p::POLY, rxns, eq_steps, enz_set, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
    check_benefit::Bool=false,
)
    result = _try_poly_power(p)
    result !== nothing && return result
    afs = _try_algebraic_factor_sigma(
        p, rxns, eq_steps, enz_set, constraints; binding_Ks,
    )
    if afs !== nothing && _expand_factored_sigma(afs) == p
        if !check_benefit
            return afs
        end
        n = sum(
            (length(fp.factors) == 1 && fp.factors[1] == poly_one()) ?
                length(c) : 1
            for (c, fp) in zip(afs.coefficients, afs.products)
        )
        n < length(p) && return afs
    end
    FactoredSigma([poly_one()], [FactoredPoly([p], [1])])
end

"""Build raw numerator POLY and factored denominator terms for the rate equation."""
function _raw_symbolic_rate_polys(M::Type{<:EnzymeMechanism})
    m = M()
    subs_species = substrates(m)
    prods_species = products(m)
    enzs = enzyme_forms(m)
    rxns = reactions(m)
    eq_steps = equilibrium_steps(m)

    enz_names = Tuple(e[1] for e in enzs)
    enz_set = Set(enz_names)
    groups, form_to_group = _compute_re_groups(enz_names, enz_set, rxns, eq_steps)
    alpha_num, alpha_den, sigma_num, sigma_den = _compute_alpha(
        enz_names, enz_set, rxns, eq_steps, groups,
    )
    G = length(groups)

    # Build rate matrix R[g1,g2] with alpha denominators cleared
    R = [poly_zero() for _ in 1:G, _ in 1:G]

    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] && continue
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        i_form = findfirst(==(e_lhs), enz_names)
        j_form = findfirst(==(e_rhs), enz_names)
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        kf_poly = poly_sym(Symbol("k$(idx)f"))
        kr_poly = poly_sym(Symbol("k$(idx)r"))

        R[g1, g2] = poly_add(
            R[g1, g2],
            _ss_contrib(
                kf_poly, m_lhs, i_form,
                alpha_num, alpha_den, groups[g1],
            ),
        )
        R[g2, g1] = poly_add(
            R[g2, g1],
            _ss_contrib(
                kr_poly, m_rhs, j_form,
                alpha_num, alpha_den, groups[g2],
            ),
        )
    end

    # Build Laplacian and cofactor determinants
    L = [i == j ? poly_zero() : poly_neg(R[i,j]) for i in 1:G, j in 1:G]
    for i in 1:G
        L[i, i] = reduce(
            poly_add,
            R[i, j] for j in 1:G if j != i;
            init=poly_zero(),
        )
    end
    D = [begin
        idx = [r for r in 1:G if r != root]
        isempty(idx) ? poly_one() : sym_det(L[idx, idx], G - 1)
    end for root in 1:G]

    # Factor sigma for each RE group; try poly_power → algebraic → unfactored
    pc = param_constraints(m)
    binding_Ks = Set{Symbol}(
        Symbol("K$i") for (i, (lhs, rhs)) in enumerate(rxns)
        if eq_steps[i] && any(s ∉ enz_set for s in lhs) &&
           all(s ∈ enz_set for s in rhs)
    )
    normalize = G == 1 && sigma_den[1] != poly_one()
    denom_terms = DenomTerm[]
    for g in 1:G
        raw_sigma = if normalize
            reduce(
                poly_add,
                (_poly_div_mono(alpha_num[i], alpha_den[i])
                 for i in groups[g]),
            )
        else
            sigma_num[g]
        end
        csigma = _apply_param_constraints(raw_sigma, pc; binding_Ks)
        push!(denom_terms, DenomTerm(
            _factor_poly(csigma, rxns, eq_steps, enz_set, pc; binding_Ks),
            D[g],
        ))
    end

    # Numerator: net flux through SS steps
    num, nu_ref = _compute_numerator(
        rxns, eq_steps, enz_names, enz_set,
        alpha_num, alpha_den, form_to_group, groups,
        D, subs_species, prods_species,
    )
    normalize && (num = _poly_div_mono(num, sigma_den[1]))

    abs_nu = abs(nu_ref)
    if abs_nu != 1
        for (i, dt) in enumerate(denom_terms)
            denom_terms[i] = DenomTerm(
                dt.sigma,
                poly_mul(dt.cofactor, poly_const(abs_nu)),
            )
        end
    end

    # Apply user-defined parameter constraints
    num = _apply_param_constraints(num, pc; binding_Ks)
    denom_terms = [_apply_param_constraints(dt, pc; binding_Ks)
                   for dt in denom_terms]

    # Merge Haldane-derived equal parameters (e.g., k8r→k7r when both
    # resolve to the same thermodynamic expression after user constraints)
    haldane_subs = _haldane_equality_substitutions(M)
    if !isempty(haldane_subs)
        hc = [(t, 1, [(c, 1)]) for (t, c) in haldane_subs]
        num = _apply_param_constraints(num, hc)
        denom_terms = [_apply_param_constraints(dt, hc)
                       for dt in denom_terms]
    end

    n_terms = length(num) + _estimate_expanded_term_count(denom_terms)
    if n_terms > MAX_RATE_EQUATION_TERMS
        error(
            "Rate equation for this mechanism has $n_terms polynomial " *
            "terms (limit: $MAX_RATE_EQUATION_TERMS). Equations this " *
            "large take a very long time to compile and are unlikely " *
            "to be practically useful for parameter fitting.",
        )
    end

    # Factor numerator (only if it reduces display terms)
    num_fs = _factor_poly(
        num, rxns, eq_steps, enz_set, pc;
        binding_Ks, check_benefit=true,
    )

    num_fs, denom_terms
end

"""
Compute the numerator polynomial by selecting an appropriate metabolite
to track through SS steps. Returns `(num::POLY, nu_ref::Int)`.
"""
function _compute_numerator(
    rxns, eq_steps, enz_names, enz_set,
    alpha_num, alpha_den, form_to_group, groups,
    D, subs_species, prods_species,
)
    ref_name = subs_species[1][1]
    nu_ref = (count(s -> s[1] == ref_name, prods_species) -
              count(s -> s[1] == ref_name, subs_species))

    # Classify metabolites into SS vs RE step sets
    ss_mets, re_mets = Set{Symbol}(), Set{Symbol}()
    for (idx, (lhs, rhs)) in enumerate(rxns)
        _, m_lhs = _split_reaction_side(lhs, enz_set)
        _, m_rhs = _split_reaction_side(rhs, enz_set)
        target = eq_steps[idx] ? re_mets : ss_mets
        for met in m_lhs; push!(target, met); end
        for met in m_rhs; push!(target, met); end
    end

    # Choose tracking metabolite: prefer ref in SS-only, then alternate SS met
    met_name = nothing  # nothing = sum all SS fluxes (G=1 fallback)
    nu_met = nu_ref
    if ref_name ∈ ss_mets && ref_name ∉ re_mets
        met_name = ref_name
    elseif !isempty(ss_mets)
        all_mets = Dict{Symbol, Int}()
        for (n, _) in subs_species; all_mets[n] = get(all_mets, n, 0) - 1; end
        for (n, _) in prods_species; all_mets[n] = get(all_mets, n, 0) + 1; end
        ss_only = setdiff(ss_mets, re_mets)
        search = isempty(ss_only) ? ss_mets : ss_only
        met_name = something(
            iterate(m for m in search if get(all_mets, m, 0) != 0),
            (first(ss_mets),),
        )[1]
        nu_met = get(all_mets, met_name, 0)
    end

    # Compute flux through SS steps
    result = poly_zero()
    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] && continue
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        i_form = findfirst(==(e_lhs), enz_names)
        j_form = findfirst(==(e_rhs), enz_names)
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        rf = _ss_contrib(
            poly_sym(Symbol("k$(idx)f")), m_lhs, i_form,
            alpha_num, alpha_den, groups[g1],
        )
        rr = _ss_contrib(
            poly_sym(Symbol("k$(idx)r")), m_rhs, j_form,
            alpha_num, alpha_den, groups[g2],
        )
        flux = poly_sub(poly_mul(rf, D[g1]), poly_mul(rr, D[g2]))
        if met_name === nothing
            result = poly_add(result, flux)
        else
            met_f = isempty(m_lhs) ? nothing : first(m_lhs)
            met_r = isempty(m_rhs) ? nothing : first(m_rhs)
            (met_f !== met_name && met_r !== met_name) && continue
            result = met_f === met_name ?
                poly_add(result, flux) : poly_sub(result, flux)
        end
    end

    # Adjust for stoichiometric ratio between tracked and reference metabolite
    if nu_met != 0 && nu_met != nu_ref
        ratio = nu_ref // nu_met
        if ratio == -1
            result = poly_neg(result)
        elseif ratio != 1
            error("Non-unit stoichiometric ratio not supported")
        end
    end

    result, nu_ref
end

# ─── Expr generation from POLY ──────────────────────────────

"""
Identify K symbols for binding RE steps (where K should be Kd, not Ka).
Canonical form invariant: all RE metabolite steps have metabolite on LHS,
so a binding step is simply any RE step with a non-enzyme species on LHS.
"""
function _binding_K_symbols(M::Type{<:EnzymeMechanism})
    m = M()
    rxns = reactions(m)
    eq = equilibrium_steps(m)
    enz_set = Set(e[1] for e in enzyme_forms(m))
    [Symbol("K$i") for (i, (lhs, _)) in enumerate(rxns)
     if eq[i] && any(s ∉ enz_set for s in lhs)]
end

"""
Compute the raw rate expression (bare symbols) and sorted parameter/concentration symbols.
Returns `(expr, all_params, sorted_concs)`.
"""
function _raw_rate_expr_and_symbols(M::Type{<:EnzymeMechanism})
    num, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    constrained = Set(c[1] for c in param_constraints(m))
    raw_ps = _raw_param_symbols(equilibrium_steps(m))
    param_syms = Set{Symbol}(
        p for p in raw_ps if p ∉ constrained
    )
    conc_syms = Set{Symbol}(metabolites(m))
    inv_set = Set(K for K in _binding_K_symbols(M) if K ∉ constrained)
    expr = to_rate_expr(num, denom_terms, param_syms, conc_syms, inv_set)
    all_params = _sorted_raw_param_symbols(M)
    return expr, all_params, metabolites(m)
end

# ─── Mode-dispatched rate_equation ────────────────────────────

"""
    rate_equation(m::EnzymeMechanism, concs, params, [mode])

Compute the QSSA steady-state rate. The body is generated at compile time
as a single arithmetic expression with no allocations, loops, or matrix ops.
"""
function rate_equation end

rate_equation(m::_AnyMechanism, concs, params) = rate_equation(m, concs, params, Reduced)

@generated function rate_equation(
    m::M, concs::NamedTuple, params::NamedTuple, ::FullMode,
) where {M <: EnzymeMechanism}
    _build_rate_body(M, FullMode)
end

@generated function rate_equation(
    m::M, concs::NamedTuple, params::NamedTuple,
    ::ReducedMode,
) where {M <: EnzymeMechanism}
    _build_rate_body(M, ReducedMode)
end

# ─── String Representation ────────────────────────────────────

"""
    rate_equation_string(m, [mode])

Return a string representation of the rate equation.
"""
function rate_equation_string end

rate_equation_string(m::_AnyMechanism) = rate_equation_string(m, Reduced)

"""Build the `v = E_total * (num) / (den)` line from the raw symbolic rate polys."""
function _rate_v_line(M::Type{<:EnzymeMechanism})
    num_fs, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    constrained = Set(c[1] for c in param_constraints(m))
    raw_ps = _raw_param_symbols(equilibrium_steps(m))
    ps = Set{Symbol}(p for p in raw_ps if p ∉ constrained)
    cs = Set{Symbol}(metabolites(m))
    inv = Set(K for K in _binding_K_symbols(M) if K ∉ constrained)
    num_str = _expr_to_string(
        _factored_sigma_to_expr(num_fs, ps, cs, inv),
    )
    den_str = _expr_to_string(
        _denom_terms_to_expr(denom_terms, ps, cs, inv),
    )
    "v = E_total * ($num_str) / ($den_str)"
end

function rate_equation_string(::M, ::FullMode) where {M<:EnzymeMechanism}
    lines = ["(; $(join(_sorted_raw_param_symbols(M), ", "))) = params",
             "(; $(join(metabolites(M()), ", "))) = concs"]
    for (target, coeff, factors) in param_constraints(M())
        push!(lines, _user_constraint_to_string(target, coeff, factors))
    end
    push!(lines, _rate_v_line(M))
    join(lines, "\n")
end

function rate_equation_string(::M, ::ReducedMode) where {M<:EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    join(["(; $(join((indep..., :Keq, :E_total), ", "))) = params",
          "(; $(join(metabolites(M()), ", "))) = concs",
          _constraint_expr_strings(M)...,
          _rate_v_line(M)], "\n")
end

# ─── Structural Identifiability ───────────────────────────────

"""
    structural_identifiability_deficit(m::EnzymeMechanism) → Int

Number of excess parameters beyond what is structurally identifiable from kinetic data.
"""
@generated function structural_identifiability_deficit(::M) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    n_k = length(indep)
    num_fs, denom_terms = _raw_symbolic_rate_polys(M)
    num = _expand_factored_sigma(num_fs)
    den = _expand_to_poly(denom_terms)
    mm(mono) = sort([
        s => e for (s, e) in mono
        if !is_k_parameter(s) && s != :E_total && s != :Keq
    ])
    n_num = length(unique(mm(k) for (k, _) in num))
    n_denom = length(unique(mm(k) for (k, _) in den))
    n_k - (n_num - 1) - (n_denom - 1)
end

# ═══════════════════════════════════════════════════════════════════
# OligomericEnzymeMechanism rate equations (MWC)
#
# MWC rate formula (per conformation c, summed over conformations):
#   num = CatN * sum_c( L_c * N_cat_c * Q_cat_c^(CatN-1)
#             * prod(Q_reg_i_c^n_reg_i for i with n_reg_i == CatN) )
#   den = sum_c( L_c * Q_cat_c^CatN * prod(Q_reg_i_c^n_reg_i) )
#   v = E_total * num / den
#
# Regulatory sites with n_reg_i < CatN appear only in the denominator.
# Sites with n_reg_i == CatN appear in both numerator and denominator.
# ═══════════════════════════════════════════════════════════════════

# ─── Parameter naming ────────────────────────────────────────────

"""Append `_T` suffix to any K or k parameter symbol (not Keq, E_total)."""
function _rename_params_T(sym::Symbol)
    is_k_parameter(sym) ? Symbol(string(sym) * "_T") : sym
end

"""Build T-state parameter substitution dictionary from a collection of parameter symbols."""
_build_T_subs(params) = Dict(p => _rename_params_T(p) for p in params if is_k_parameter(p))

"""Name for a regulatory site parameter: K_{ligand}_reg{i} or K_{ligand}_T_reg{i}."""
function _reg_param_name(ligand::Symbol, site_idx::Int, T_state::Bool)
    T_state ? Symbol("K_$(ligand)_T_reg$(site_idx)") :
              Symbol("K_$(ligand)_reg$(site_idx)")
end

# ─── Binding K symbols ───────────────────────────────────────────

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

# ─── Dependent parameter expressions ─────────────────────────────

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

    # NConf == 2: rename all k/K params (not Keq) with _T suffix
    T_subs = _build_T_subs(Iterators.flatten([keys(dep_R), indep_R]))
    dep_T = Dict{Symbol, Union{Symbol, Expr}}(
        _rename_params_T(k) => substitute_params_expr(v, T_subs)
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

# parameters and fitted_params for OligomericEnzymeMechanism are handled
# by the unified _AnyMechanism methods at the top of this file.

# ─── Rate body building helpers ───────────────────────────────────

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

"""
Build R-state and T-state dep param assignment Exprs.
Returns `(r_assignments::Vector{Expr}, t_assignments::Vector{Expr})`.
Shared by `_build_oligomeric_rate_body` and `rate_equation_string`.
"""
function _oligomeric_dep_assignments(
    CM::Type{<:EnzymeMechanism}, NConf, inv_fn,
)
    dep_R, indep_R = _dependent_param_exprs(CM)
    dep_R_kd = _apply_kd_inversion(dep_R, CM, inv_fn)
    sorted_deps = sort(collect(dep_R_kd); by=first)

    r_assignments = Expr[
        Expr(:(=), sym, dep_R_kd[sym]) for (sym, _) in sorted_deps
    ]

    t_assignments = if NConf == 2
        T_subs = _build_T_subs(Iterators.flatten([keys(dep_R), indep_R]))
        Expr[
            Expr(:(=), _rename_params_T(sym),
                substitute_params_expr(dep_R_kd[sym], T_subs))
            for (sym, _) in sorted_deps
        ]
    else
        Expr[]
    end

    return r_assignments, t_assignments
end

"""
Assemble the MWC numerator and denominator Exprs.
Returns `(full_num, full_den)` where the numerator already includes the `CatN` factor.
Shared by `_build_oligomeric_rate_body` and `rate_equation_string`.
"""
function _oligomeric_num_den_exprs(CM, CatN, RS, NConf)
    num_fs, denom_terms = _raw_symbolic_rate_polys(CM)
    m_cat = CM()
    cat_constr = Set(c[1] for c in param_constraints(m_cat))
    cat_params = Set{Symbol}(
        p for p in _raw_param_symbols(equilibrium_steps(m_cat)) if p ∉ cat_constr
    )
    cat_mets = Set{Symbol}(metabolites(m_cat))
    binding_Ks_r = Set(_binding_K_symbols(CM))

    N_R = _factored_sigma_to_expr(num_fs, cat_params, cat_mets, binding_Ks_r)
    Q_R = _denom_terms_to_expr(denom_terms, cat_params, cat_mets, binding_Ks_r)
    reg_Q_R = Any[_reg_site_expr(ligs, i, false) for (i, (ligs, _)) in enumerate(RS)]

    # Per-subunit reg sites (n_reg == CatN) go in both num and den.
    # Enzyme-level reg sites (n_reg < CatN) go in denominator only.
    function make_num_term(N, Q, reg_Qs)
        factors = Any[N]
        CatN > 1 && push!(factors, _power_expr(Q, CatN - 1))
        for (idx, (_, n_reg)) in enumerate(RS)
            n_reg == CatN || continue
            push!(factors, _power_expr(reg_Qs[idx], n_reg))
        end
        _nest_binary(:*, factors)
    end

    function make_den_term(Q, reg_Qs)
        factors = Any[_power_expr(Q, CatN)]
        for (idx, (_, n_reg)) in enumerate(RS)
            push!(factors, _power_expr(reg_Qs[idx], n_reg))
        end
        _nest_binary(:*, factors)
    end

    num_R = make_num_term(N_R, Q_R, reg_Q_R)
    den_R = make_den_term(Q_R, reg_Q_R)

    NConf == 1 && return :($(CatN) * $(num_R)), den_R

    # NConf == 2: T-state
    T_subs = _build_T_subs(cat_params)
    N_T = substitute_params_expr(N_R, T_subs)
    Q_T = substitute_params_expr(Q_R, T_subs)
    reg_Q_T = Any[_reg_site_expr(ligs, i, true) for (i, (ligs, _)) in enumerate(RS)]

    num_T = make_num_term(N_T, Q_T, reg_Q_T)
    den_T = make_den_term(Q_T, reg_Q_T)

    :($(CatN) * ($(num_R) + L * $(num_T))), :($(den_R) + L * $(den_T))
end

"""Build the MWC rate equation body as an Expr block."""
function _build_oligomeric_rate_body(Mets, CM, CatN, RS, NConf)
    M_type = OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}
    full_num, full_den = _oligomeric_num_den_exprs(CM, CatN, RS, NConf)
    rate_expr = :(E_total * ($full_num) / ($full_den))

    r_assignments, t_assignments = _oligomeric_dep_assignments(
        CM, NConf, K -> :(inv($K)),
    )

    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq, :E_total)

    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(Mets, :concs),
        r_assignments...,
        t_assignments...,
        :(return $rate_expr))
end

# ─── Rate equation dispatch ───────────────────────────────────────

@generated function rate_equation(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf},
    concs::NamedTuple, params::NamedTuple, ::ReducedMode,
) where {Mets,CM,CatN,RS,NConf}
    _build_oligomeric_rate_body(Mets, CM, CatN, RS, NConf)
end

# ─── String representation ────────────────────────────────────────

function rate_equation_string(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}, ::ReducedMode,
) where {Mets,CM,CatN,RS,NConf}
    M = OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)

    r_assignments, t_assignments = _oligomeric_dep_assignments(
        CM, NConf, K -> :(1 / $K),
    )
    dep_lines = ["$(a.args[1]) = $(_expr_to_string(a.args[2]))" for a in r_assignments]
    t_dep_lines = ["$(a.args[1]) = $(_expr_to_string(a.args[2]))" for a in t_assignments]

    full_num, full_den = _oligomeric_num_den_exprs(CM, CatN, RS, NConf)
    v_line = "v = E_total * ($(_expr_to_string(full_num))) / ($(_expr_to_string(full_den)))"

    join([
        "(; $(join(hw_params, ", "))) = params",
        "(; $(join(Mets, ", "))) = concs",
        dep_lines...,
        t_dep_lines...,
        v_line,
    ], "\n")
end

# ─── Structural Identifiability ───────────────────────────────────

@generated function structural_identifiability_deficit(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf},
) where {Mets,CM,CatN,RS,NConf}
    M = OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}
    _, indep = _dependent_param_exprs(M)
    n_k = length(indep)
    n_num, n_denom = _count_oligomeric_rate_monomials(CM, CatN, RS, NConf)
    n_k - (n_num - 1) - (n_denom - 1)
end

