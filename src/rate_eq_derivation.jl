# ─── Parameters API ─────────────────────────────────────────

const _AnyMechanism = AbstractEnzymeMechanism

"""
    parameters(m::EnzymeMechanism, [mode])
    parameters(m::AllostericEnzymeMechanism, [mode])

Return the parameter names required for the given mode as a tuple of Symbols.

# Modes
- `Reduced` (default): independent k's + Keq + E_total
- `Full`: all 2N k's + E_total (EnzymeMechanism only)
"""
function parameters end

parameters(m::_AnyMechanism) = parameters(m, Reduced)

@generated function parameters(::M, ::FullMode) where {M <: EnzymeMechanism}
    Tuple((_raw_param_symbols(M())..., :E_total))
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

"""
Build a renaming map from non-representative step parameter symbols to the
representative step's parameter symbols. Used to alias `K2` → `K1` (etc.) when
steps 1 and 2 share a kinetic group.
"""
function _build_kinetic_rename_map(m::EnzymeMechanism)
    rename = Dict{Symbol, Symbol}()
    eq = equilibrium_steps(m)
    for g in kinetic_groups(m)
        idxs = steps_in_group(m, g)
        length(idxs) == 1 && continue
        rep = first(idxs)
        for idx in idxs
            idx == rep && continue
            if eq[idx]
                rename[Symbol("K$idx")] = Symbol("K$rep")
            else
                rename[Symbol("k$(idx)f")] = Symbol("k$(rep)f")
                rename[Symbol("k$(idx)r")] = Symbol("k$(rep)r")
            end
        end
    end
    rename
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
    for (idx, (lhs, rhs, _, _)) in enumerate(rxns)
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
            for (idx, (lhs, rhs, _, _)) in enumerate(rxns)
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
Try algebraic polynomial factoring on a sigma.
Tries each metabolite present in sigma: splits terms by metabolite presence,
divides the metabolite-containing terms by K*met, and checks whether
the result cleanly divides sigma into `base * (... + K*met)`.
`rename_map` aliases non-representative kinetic-group K symbols to their
representative — used to find the canonical binding K for each metabolite.
Returns a `FactoredSigma` or `nothing`.
"""
function _try_algebraic_factor_sigma(
    sigma::POLY, rxns, eq_steps, enz_set, rename_map;
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

    # Try each metabolite, pick best factoring by score
    best = nothing
    best_score = (length(sigma), 0, 0)
    for met in met_syms
        # Find canonical binding K for this metabolite
        K_R = nothing
        for (idx, (lhs, rhs_rxn, _, _)) in enumerate(rxns)
            eq_steps[idx] || continue
            _, m_lhs = _split_reaction_side(lhs, enz_set)
            _, m_rhs = _split_reaction_side(rhs_rxn, enz_set)
            ((met ∈ m_lhs) ⊻ (met ∈ m_rhs)) || continue
            K_R = Symbol("K$idx")
            K_R = get(rename_map, K_R, K_R)
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
                    base_R, rxns, eq_steps, enz_set, rename_map;
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
function _haldane_equality_substitutions(dep_exprs)
    length(dep_exprs) < 2 && return Pair{Symbol,Symbol}[]
    # Sort by step number for stable canonical choice (lowest first)
    _step_num(s) = let m = match(r"\d+", string(s))
        m === nothing ? 0 :
            something(tryparse(Int, m.match::SubString), 0)
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
    p::POLY, rxns, eq_steps, enz_set, rename_map;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
    check_benefit::Bool=false,
)
    result = _try_poly_power(p)
    result !== nothing && return result
    afs = _try_algebraic_factor_sigma(
        p, rxns, eq_steps, enz_set, rename_map; binding_Ks,
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

"""Build raw numerator POLY and factored denominator terms for the rate equation.
`rename_map` aliases non-representative kinetic-group parameter symbols to their
representative (one K or k_f/k_r per group); `dep_exprs` is the Haldane solution."""
function _raw_symbolic_rate_polys(
    subs_species, prods_species, enz_names, enz_set,
    rxns, eq_steps, rename_map, dep_exprs,
)
    groups, form_to_group = _compute_re_groups(
        enz_names, enz_set, rxns, eq_steps,
    )
    alpha_num, alpha_den, sigma_num, sigma_den =
        _compute_alpha(
            enz_names, enz_set, rxns, eq_steps, groups,
        )
    G = length(groups)

    # Build rate matrix R[g1,g2] with alpha denominators cleared
    R = [poly_zero() for _ in 1:G, _ in 1:G]

    for (idx, (lhs, rhs, _, _)) in enumerate(rxns)
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
    L = [i == j ? poly_zero() : poly_neg(R[i,j])
         for i in 1:G, j in 1:G]
    for i in 1:G
        L[i, i] = reduce(
            poly_add,
            R[i, j] for j in 1:G if j != i;
            init=poly_zero(),
        )
    end
    D = [begin
        idx = [r for r in 1:G if r != root]
        isempty(idx) ? poly_one() :
            sym_det(L[idx, idx], G - 1)
    end for root in 1:G]

    # Factor sigma for each RE group
    binding_Ks = Set{Symbol}(
        Symbol("K$i")
        for (i, (lhs, rhs, _, _)) in enumerate(rxns)
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
        csigma = _rename_symbols(raw_sigma, rename_map)
        push!(denom_terms, DenomTerm(
            _factor_poly(
                csigma, rxns, eq_steps, enz_set, rename_map;
                binding_Ks,
            ),
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

    # Apply kinetic-group renaming (K2 → K1 etc.) to numerator
    num = _rename_symbols(num, rename_map)
    denom_terms = [_rename_symbols(dt, rename_map) for dt in denom_terms]

    # Merge Haldane-derived equal parameters (e.g., k8r→k7r when
    # both resolve to the same thermodynamic expression)
    haldane_subs = _haldane_equality_substitutions(dep_exprs)
    if !isempty(haldane_subs)
        hsub_map = Dict{Symbol, Symbol}(haldane_subs)
        num = _rename_symbols(num, hsub_map)
        denom_terms = [_rename_symbols(dt, hsub_map)
                       for dt in denom_terms]
    end

    n_terms = (length(num) +
               _estimate_expanded_term_count(denom_terms))
    if n_terms > MAX_RATE_EQUATION_TERMS
        error(
            "Rate equation for this mechanism has " *
            "$n_terms polynomial terms " *
            "(limit: $MAX_RATE_EQUATION_TERMS). Equations " *
            "this large take a very long time to compile " *
            "and are unlikely to be practically useful " *
            "for parameter fitting.",
        )
    end

    # Factor numerator (only if it reduces display terms)
    num_fs = _factor_poly(
        num, rxns, eq_steps, enz_set, rename_map;
        binding_Ks, check_benefit=true,
    )

    num_fs, denom_terms
end

function _raw_symbolic_rate_polys(M::Type{<:EnzymeMechanism})
    m = M()
    enz_names = enzyme_forms(m)
    enz_set = Set(enz_names)
    rxns = reactions(m)
    eq_steps = equilibrium_steps(m)
    rename_map = _build_kinetic_rename_map(m)
    dep_exprs, _ = _dependent_param_exprs(M)
    _raw_symbolic_rate_polys(
        substrates(m), products(m), enz_names, enz_set,
        rxns, eq_steps, rename_map, dep_exprs,
    )
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
    ref_name = subs_species[1]
    nu_ref = (count(==(ref_name), prods_species) -
              count(==(ref_name), subs_species))

    # Classify metabolites into SS vs RE step sets
    ss_mets, re_mets = Set{Symbol}(), Set{Symbol}()
    for (idx, (lhs, rhs, _, _)) in enumerate(rxns)
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
        for n in subs_species; all_mets[n] = get(all_mets, n, 0) - 1; end
        for n in prods_species; all_mets[n] = get(all_mets, n, 0) + 1; end
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
    for (idx, (lhs, rhs, _, _)) in enumerate(rxns)
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
The returned symbols use kinetic-group representatives (one per RE binding group).
"""
function _binding_K_symbols(M::Type{<:EnzymeMechanism})
    m = M()
    rxns = reactions(m)
    eq_steps = equilibrium_steps(m)
    enz_set = Set(enzyme_forms(m))
    rename_map = _build_kinetic_rename_map(m)
    seen = Set{Symbol}()
    syms = Symbol[]
    for (i, (lhs, _, _, _)) in enumerate(rxns)
        eq_steps[i] || continue
        any(s ∉ enz_set for s in lhs) || continue
        K = Symbol("K$i")
        K = get(rename_map, K, K)
        K in seen && continue
        push!(seen, K)
        push!(syms, K)
    end
    syms
end

"""
Compute the raw rate expression (bare symbols) and sorted parameter/concentration symbols.
Returns `(expr, all_params, sorted_concs)`.
"""
function _raw_rate_expr_and_symbols(M::Type{<:EnzymeMechanism})
    num, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    param_syms = Set{Symbol}(_raw_param_symbols(m))
    conc_syms = Set{Symbol}(metabolites(m))
    inv_set = Set(_binding_K_symbols(M))
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
    ps = Set{Symbol}(_raw_param_symbols(m))
    cs = Set{Symbol}(metabolites(m))
    inv = Set(_binding_K_symbols(M))
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

# ─── kcat Computation Helpers ──────────────────────────────────

"""Check if a symbol is a steady-state rate constant (lowercase k followed by digit)."""
function _is_ss_rate_constant(sym::Symbol)
    s = string(sym)
    length(s) > 1 && s[1] == 'k' && isdigit(s[2])
end

"""Group `num` and `den` POLYs by metabolite monomial pattern.
Returns `(num_groups, den_groups)` where each value is a POLY of k-monomials
sharing the same met-monomial. Reverse (negative-coefficient) terms are dropped
from the numerator. Used by `_kcat_components` and `_kcat_forward` to compare
saturating metabolite patterns across R/T states."""
function _kcat_groups_from_polys(num::POLY, den::POLY)
    function split_mono(mono::MONO)
        k_mono = MONO()
        met_mono = MONO()
        for (s, e) in mono
            if is_k_parameter(s) || s == :Keq
                push!(k_mono, s => e)
            elseif s != :E_total
                push!(met_mono, s => e)
            end
        end
        sort!(k_mono; by=first), sort!(met_mono; by=first)
    end

    num_groups = Dict{MONO, POLY}()
    for (mono, coeff) in num
        coeff > 0 || continue
        k_part, met_part = split_mono(mono)
        p = get!(num_groups, met_part, POLY())
        p[k_part] = get(p, k_part, Rational{Int}(0)) + coeff
    end
    den_groups = Dict{MONO, POLY}()
    for (mono, coeff) in den
        k_part, met_part = split_mono(mono)
        p = get!(den_groups, met_part, POLY())
        p[k_part] = get(p, k_part, Rational{Int}(0)) + coeff
    end
    num_groups, den_groups
end

"""
Compute kcat candidate components from the rate equation polynomial structure.
Returns `Vector{Tuple{Any, Any}}` of (num_k_expr, den_k_expr) pairs.
kcat = max(nk/dk) over all candidates (evaluated at runtime).

Groups the expanded numerator (forward terms only, positive coefficients) and
denominator by metabolite pattern. For each matching pair, builds k-only
expressions. K's (equilibrium constants) cancel between num and den at
matching metabolite levels.

Multiple candidates arise for mechanisms with alternative catalytic pathways
(e.g., non-essential activator with/without activator bound).
"""
function _kcat_components(M::Type{<:EnzymeMechanism})
    num_fs, denom_terms = _raw_symbolic_rate_polys(M)
    num = _expand_factored_sigma(num_fs)
    den = _expand_to_poly(denom_terms)
    num_groups, den_groups = _kcat_groups_from_polys(num, den)

    # Build kcat candidates: for each forward numerator metabolite group
    # with a matching denominator group, create (num_k_expr, den_k_expr)
    empty_set = Set{Symbol}()
    components = Tuple{Any, Any}[]
    for (met_key, num_k) in sort!(collect(num_groups); by=first)
        den_k = get(den_groups, met_key, nothing)
        den_k === nothing && continue
        num_expr = _poly_to_expr(num_k, empty_set, empty_set)
        den_expr = _poly_to_expr(den_k, empty_set, empty_set)
        push!(components, (num_expr, den_expr))
    end

    components
end

"""
    _kcat_forward(m::EnzymeMechanism, params) → Float64

Compute kcat (forward) analytically from the polynomial structure.
kcat is the maximum rate at saturating substrates, zero products,
and E_total=1. For mechanisms with multiple catalytic paths
(e.g., non-essential activator), returns the max over all paths.

Uses `_kcat_components` to get (num_k_expr, den_k_expr) candidates,
then evaluates max(nk/dk) at runtime parameter values.
"""
@generated function _kcat_forward(
    ::M, params::NamedTuple,
) where {M <: EnzymeMechanism}
    components = _kcat_components(M)
    dep_exprs, indep = _dependent_param_exprs(M)
    dep_exprs = _apply_kd_inversion(dep_exprs, M, K -> :(inv($K)))
    hw_params = (indep..., :Keq)
    assignments = [Expr(:(=), sym, dep_exprs[sym])
                   for (sym, _) in sort(collect(dep_exprs); by=first)]
    # Apply Kd inversion to component expressions: raw polys use Ka,
    # but params store Kd for binding K's
    binding_Ks = Set(_binding_K_symbols(M))
    kd_subs = Dict(K => :(inv($K)) for K in binding_Ks)
    candidates = [
        :($(substitute_params_expr(nk, kd_subs)) /
          $(substitute_params_expr(dk, kd_subs)))
        for (nk, dk) in components
    ]
    result = length(candidates) == 1 ?
        candidates[1] : Expr(:call, :max, candidates...)
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        assignments...,
        :(return $result))
end

"""
    _kcat_forward(m::AllostericEnzymeMechanism, params) → Float64

Compute kcat (forward) for an MWC allosteric enzyme.

kcat is the maximum rate at saturating substrates, zero products,
E_total=1, over all regulator concentration corners (each regulator
either 0 or saturating).

With regulatory sites: regulators shift R/T balance,
so kcat depends on the regulator corner. We enumerate all 2^n_lig
corners and return the max.
"""
@generated function _kcat_forward(
    ::AllostericEnzymeMechanism{CM,CS,RS},
    params::NamedTuple,
) where {CM,CS,RS}
    M_type = AllostericEnzymeMechanism{CM,CS,RS}
    m = M_type()
    CatN = catalytic_multiplicity(m)

    # Build R-state and T-state polynomials separately so the saturating
    # metabolite pattern can be matched across conformations. T-state zeros
    # `:OnlyR` k symbols at the polynomial level — patterns that only exist
    # via an `:OnlyR` step (e.g. an `:OnlyR` substrate-binding K) drop out
    # of T-state, so their `(num_k_T, den_k_T)` are 0 and the T-state
    # contribution at saturation vanishes.
    num_fs, denom_terms = _raw_symbolic_rate_polys(CM)
    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)
    rename_T = _T_rename(m)
    num_R_poly = _zero_symbols_in_poly(_expand_factored_sigma(num_fs), t_only_syms)
    den_R_poly = _zero_symbols_in_poly(_expand_to_poly(denom_terms), t_only_syms)
    num_T_poly = _rename_symbols(
        _zero_symbols_in_poly(_expand_factored_sigma(num_fs), r_only_syms),
        rename_T)
    den_T_poly = _rename_symbols(
        _zero_symbols_in_poly(_expand_to_poly(denom_terms), r_only_syms),
        rename_T)
    num_R_groups, den_R_groups = _kcat_groups_from_polys(num_R_poly, den_R_poly)
    num_T_groups, den_T_groups = _kcat_groups_from_polys(num_T_poly, den_T_poly)

    # Choose the saturating R-state met pattern (single component for
    # mechanisms exercised here; assert keeps that constraint visible).
    r_keys = sort!([k for k in keys(num_R_groups) if haskey(den_R_groups, k)])
    @assert length(r_keys) == 1 "Catalytic mechanism should have exactly 1 kcat component"
    met_key = r_keys[1]
    empty_set = Set{Symbol}()
    raw_num_k_R = _poly_to_expr(num_R_groups[met_key], empty_set, empty_set)
    raw_den_k_R = _poly_to_expr(den_R_groups[met_key], empty_set, empty_set)

    # Apply Kd inversion: raw polys use Ka convention, params use Kd
    binding_Ks_R = Set(_binding_K_symbols(CM))
    binding_Ks_T = Set(get(rename_T, K, K) for K in binding_Ks_R)
    num_k_R_expr = substitute_params_expr(
        raw_num_k_R, Dict(K => :(inv($K)) for K in binding_Ks_R))
    den_k_R_expr = substitute_params_expr(
        raw_den_k_R, Dict(K => :(inv($K)) for K in binding_Ks_R))

    # T-state at the same metabolite pattern. Missing → set 0 (the saturating
    # pattern is unreachable in T-state, so T contributes neither flux nor
    # enzyme mass at saturation along that pattern).
    kd_subs_T = Dict(K => :(inv($K)) for K in binding_Ks_T)
    num_T_p = get(num_T_groups, met_key, nothing)
    den_T_p = get(den_T_groups, met_key, nothing)
    num_k_T_expr = num_T_p === nothing ? 0 :
        substitute_params_expr(
            _poly_to_expr(num_T_p, empty_set, empty_set), kd_subs_T)
    den_k_T_expr = den_T_p === nothing ? 0 :
        substitute_params_expr(
            _poly_to_expr(den_T_p, empty_set, empty_set), kd_subs_T)
    t_pattern_dead = den_T_p === nothing

    # Build dependent param assignments for R-state
    r_assignments, t_assignments_ =
        _build_dep_assignments(M_type, K -> :(inv($K)))

    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq)

    t_state_dead = _t_state_dead(m)

    # A_c = num_k_c * den_k_c^(CatN-1), B_c = den_k_c^CatN
    A_R = CatN == 1 ? num_k_R_expr :
        :($(num_k_R_expr) * $(den_k_R_expr)^$(CatN - 1))
    B_R = :($(den_k_R_expr)^$(CatN))

    # T-state catalytic A/B. `A_T` is zero whenever the catalytic cycle is
    # broken; `B_T` additionally drops to 0 if the saturating pattern is
    # absent in T-state (its mass is lower-order in subs and vanishes
    # relative to R-state at saturation).
    if t_pattern_dead
        A_T = 0
        B_T = 0
    else
        A_T = t_state_dead ? 0 :
              CatN == 1 ? num_k_T_expr :
                          :($(num_k_T_expr) * $(den_k_T_expr)^$(CatN - 1))
        B_T = :($(den_k_T_expr)^$(CatN))
    end

    # Skip T-state dep-assignments when T-state is dead (they may reference
    # zero-valued params and cause Inf via 1/K substitutions).
    t_assignments = (t_state_dead || t_pattern_dead) ? Expr[] : t_assignments_

    if isempty(RS)
        # No regulatory sites: single kcat value
        result = :($(CatN) * ($(A_R) + L * $(A_T)) /
                   ($(B_R) + L * $(B_T)))
        return Expr(:block,
            _destructuring_expr(hw_params, :params),
            r_assignments...,
            t_assignments...,
            :(return $result))
    end

    # Enumerate all unique ligands across reg sites
    all_ligs = Symbol[]
    for entry in RS
        for lig in entry[1]
            lig in all_ligs || push!(all_ligs, lig)
        end
    end
    n_ligs = length(all_ligs)
    lig_idx = Dict(lig => i - 1 for (i, lig) in enumerate(all_ligs))

    # For each corner (2^n_ligs), compute kcat expression
    # At corner: each lig is either 0 or ∞.
    # At ∞: Q_reg_c_i → sum(inv(K_j_c_i) for saturating j)
    # At 0: Q_reg_c_i = 1
    corner_exprs = Any[]
    for mask in 0:(2^n_ligs - 1)
        W_R_factors = Any[]
        W_T_factors = Any[]
        for site_idx in eachindex(RS)
            ligs = regulatory_site_ligands(m, site_idx)
            n_reg = regulatory_site_multiplicity(m, site_idx)
            sat_terms_R = Any[]
            sat_terms_T = Any[]
            for lig in ligs
                if (mask >> lig_idx[lig]) & 1 == 1
                    tag = regulatory_ligand_tag(m, site_idx, lig)
                    if tag != :OnlyT
                        push!(sat_terms_R,
                              :(inv($(_reg_param_name(lig, site_idx, false)))))
                    end
                    if tag != :OnlyR
                        K_T = _reg_param_name(
                            lig, site_idx, tag in (:NonequalRT, :OnlyT))
                        push!(sat_terms_T, :(inv($K_T)))
                    end
                end
            end
            if !isempty(sat_terms_R)
                q_R = length(sat_terms_R) == 1 ?
                    sat_terms_R[1] :
                    _nest_binary(:+, sat_terms_R)
                push!(W_R_factors, _power_expr(q_R, n_reg))
            end
            if !isempty(sat_terms_T)
                q_T = length(sat_terms_T) == 1 ?
                    sat_terms_T[1] :
                    _nest_binary(:+, sat_terms_T)
                push!(W_T_factors, _power_expr(q_T, n_reg))
            end
        end
        # Build kcat at this corner. Empty W_R_factors / W_T_factors mean
        # no saturating regulator at that conformation — they default to 1.
        if isempty(W_R_factors) && isempty(W_T_factors)
            # Zero regulators corner
            kcat_expr = :($(CatN) * ($(A_R) + L * $(A_T)) /
                          ($(B_R) + L * $(B_T)))
        else
            W_R = isempty(W_R_factors) ? 1 :
                length(W_R_factors) == 1 ?
                    W_R_factors[1] :
                    _nest_binary(:*, W_R_factors)
            W_T = isempty(W_T_factors) ? 1 :
                length(W_T_factors) == 1 ?
                    W_T_factors[1] :
                    _nest_binary(:*, W_T_factors)
            kcat_expr = :($(CatN) *
                ($(A_R) * $(W_R) + L * $(A_T) * $(W_T)) /
                ($(B_R) * $(W_R) + L * $(B_T) * $(W_T)))
        end
        push!(corner_exprs, kcat_expr)
    end

    result = length(corner_exprs) == 1 ?
        corner_exprs[1] : Expr(:call, :max, corner_exprs...)
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        r_assignments...,
        t_assignments...,
        :(return $result))
end

# ─── Public API: rescale_parameter_values ──────────────────────────

"""
    rescale_parameter_values(m, params::NamedTuple; kcat=1.0)

Rescale SS rate constants so that `_kcat_forward(m, result) ≈ kcat`.
Non-SS parameters (K's, Keq, E_total, L, regulatory K's) are unchanged.
"""
function rescale_parameter_values(
    m::_AnyMechanism, params::NamedTuple; kcat=1.0,
)
    kcat_current = _kcat_forward(m, params)
    scale = kcat / kcat_current
    NamedTuple{keys(params)}(Tuple(
        _is_ss_rate_constant(k) ? v * scale : v
        for (k, v) in zip(keys(params), values(params))
    ))
end

# ═══════════════════════════════════════════════════════════════════
# AllostericEnzymeMechanism rate equations (MWC)
#
# MWC rate formula (per conformation c, summed over conformations).
# Let cat_n = catalytic_multiplicity(m):
#   num = cat_n * sum_c( L_c * N_cat_c * Q_cat_c^(cat_n - 1)
#             * prod(Q_reg_i_c^n_reg_i for i with n_reg_i == cat_n) )
#   den = sum_c( L_c * Q_cat_c^cat_n * prod(Q_reg_i_c^n_reg_i) )
#   v = E_total * num / den
#
# Regulatory sites with n_reg_i < cat_n appear only in the denominator.
# Sites with n_reg_i == cat_n appear in both numerator and denominator.
# ═══════════════════════════════════════════════════════════════════

# ─── Parameter naming ────────────────────────────────────────────

"""Append `_T` suffix to any K or k parameter symbol (not Keq, E_total)."""
function _rename_params_T(sym::Symbol)
    is_k_parameter(sym) ? Symbol(string(sym) * "_T") : sym
end

"""Name for a regulatory site parameter: K_{ligand}_reg{i} or K_{ligand}_T_reg{i}."""
function _reg_param_name(ligand::Symbol, site_idx::Int, T_state::Bool)
    T_state ? Symbol("K_$(ligand)_T_reg$(site_idx)") :
              Symbol("K_$(ligand)_reg$(site_idx)")
end

# ─── Allosteric kinetic-group symbol helpers ─────────────────────

"""Parameter symbols (representative) for the kinetic group of step `idx`."""
function _group_param_symbols(cm::EnzymeMechanism, idx::Int)
    rep = first(steps_in_group(cm, kinetic_group(cm, idx)))
    equilibrium_steps(cm)[rep] ?
        (Symbol("K$rep"),) :
        (Symbol("k$(rep)f"), Symbol("k$(rep)r"))
end

"""Catalytic-cycle parameter symbols zeroed in the R-state (`:OnlyT` groups)."""
function _onlyT_syms(m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    syms = Set{Symbol}()
    for g in kinetic_groups(cm)
        group_tag(m, g) == :OnlyT || continue
        rep = first(steps_in_group(cm, g))
        for s in _group_param_symbols(cm, rep); push!(syms, s); end
    end
    syms
end

"""
The T-state catalytic cycle cannot close — and therefore both forward
and reverse net flux vanish — when any `:OnlyR` kinetic group's
representative step either has no metabolite (an iso step) or binds a
substrate. Both cases break a step the cycle traverses; binding
equilibrium for the rest forces N_T = 0 in steady state. Used by both
`rate_equation` (via `_allosteric_num_den_exprs`) and `_kcat_forward`.
"""
function _t_state_dead(m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    enz_set = Set(enzyme_forms(cm))
    sub_set = Set(substrates(cm))
    rxns = reactions(cm)
    for g in kinetic_groups(cm)
        group_tag(m, g) == :OnlyR || continue
        rep = first(steps_in_group(cm, g))
        _, mets_lhs = _split_reaction_side(rxns[rep][1], enz_set)
        if isempty(mets_lhs) || any(met in sub_set for met in mets_lhs)
            return true
        end
    end
    false
end

"""Catalytic-cycle parameter symbols zeroed in the T-state (`:OnlyR` groups)."""
function _onlyR_syms(m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    syms = Set{Symbol}()
    for g in kinetic_groups(cm)
        group_tag(m, g) == :OnlyR || continue
        rep = first(steps_in_group(cm, g))
        for s in _group_param_symbols(cm, rep); push!(syms, s); end
    end
    syms
end

"""
R→T rename map for groups whose T-state symbol differs from their R-state
symbol — `:NonequalRT` (independent K_R, K_T) and `:OnlyT` (T-state-only,
no R-state counterpart so the T-state symbol takes the canonical _T suffix).
"""
function _T_rename(m::AllostericEnzymeMechanism)
    cm = catalytic_mechanism(m)
    rename = Dict{Symbol, Symbol}()
    for g in kinetic_groups(cm)
        tag = group_tag(m, g)
        (tag == :NonequalRT || tag == :OnlyT) || continue
        rep = first(steps_in_group(cm, g))
        for s in _group_param_symbols(cm, rep)
            rename[s] = _rename_params_T(s)
        end
    end
    rename
end

# ─── Binding K symbols ───────────────────────────────────────────

"""
Return all binding (Kd-convention) K symbols: R-state, T-state, and reg site params.
Regulatory site K params are always Kd. `:EqualRT` and `:OnlyR` groups have no
T-state K (they share the R-state symbol or are absent in T).
"""
function _binding_K_symbols(
    ::Type{AllostericEnzymeMechanism{CM,CS,RS}},
) where {CM,CS,RS}
    m = AllostericEnzymeMechanism{CM,CS,RS}()
    cm = catalytic_mechanism(m)
    r_ks = Tuple(_binding_K_symbols(CM))
    t_ks = Symbol[]
    for K in r_ks
        idx = parse(Int, string(K)[2:end])
        tag = group_tag(m, kinetic_group(cm, idx))
        tag == :NonequalRT && push!(t_ks, _rename_params_T(K))
    end
    reg_ks_r = Symbol[]
    reg_ks_t = Symbol[]
    for (i, entry) in enumerate(RS)
        for lig in entry[1]
            tag = regulatory_ligand_tag(m, i, lig)
            tag == :OnlyT || push!(reg_ks_r, _reg_param_name(lig, i, false))
            tag in (:NonequalRT, :OnlyT) &&
                push!(reg_ks_t, _reg_param_name(lig, i, true))
        end
    end
    (r_ks..., t_ks..., reg_ks_r..., reg_ks_t...)
end

# ─── Dependent parameter expressions ─────────────────────────────

"""
    _dependent_param_exprs(M::Type{<:AllostericEnzymeMechanism})

Return `(dep_exprs, indep_params)` for an AllostericEnzymeMechanism.
R-state expressions come from `_dependent_param_exprs(CM)`; T-state entries
copy R with `_T` rename, with per-group filtering (`:OnlyR` zeroed in T,
`:OnlyT` zeroed in R, `:EqualRT` shares R symbol).
Adds reg site params and L to indep.
"""
function _dependent_param_exprs(
    ::Type{AllostericEnzymeMechanism{CM,CS,RS}},
) where {CM,CS,RS}
    m = AllostericEnzymeMechanism{CM,CS,RS}()
    dep_R_all, indep_R_all = _dependent_param_exprs(CM)
    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)

    # Filter R-state: drop dep entries that reference :OnlyT symbols
    # (those equal zero in R-state); drop :OnlyT independents.
    dep_R = Dict{Symbol, Union{Symbol, Expr}}()
    for (k, v) in dep_R_all
        _expr_references_any(v, t_only_syms) && continue
        dep_R[k] = v
    end
    indep_R = Symbol[p for p in indep_R_all if p ∉ t_only_syms]

    # Build T-state symbol substitution for non-`:EqualRT` groups
    rename_T = _T_rename(m)
    T_subs = Dict{Symbol, Symbol}(rename_T)

    dep_T = Dict{Symbol, Union{Symbol, Expr}}()
    indep_T_list = Symbol[]
    for (k, v) in dep_R_all
        _expr_references_any(v, r_only_syms) && continue
        t_k = get(rename_T, k, k)
        # `:EqualRT` (k unchanged) is already in dep_R; skip duplicate.
        t_k == k && continue
        dep_T[t_k] = substitute_params_expr(v, T_subs)
    end
    for p in indep_R_all
        p ∈ r_only_syms && continue
        if !haskey(rename_T, p)
            # `:EqualRT` independent: its T-state mirror equals the R-state
            # symbol. Add p_T = p as a dep, do not duplicate as indep.
            dep_T[_rename_params_T(p)] = p
        else
            # `:NonequalRT` and `:OnlyT` get a distinct T-state independent.
            push!(indep_T_list, _rename_params_T(p))
        end
    end

    # Reg-site parameters
    reg_params_r = Symbol[]
    reg_params_t_indep = Symbol[]
    reg_params_t_dep = Pair{Symbol, Symbol}[]
    for (i, entry) in enumerate(RS)
        for lig in entry[1]
            tag = regulatory_ligand_tag(m, i, lig)
            K_R = _reg_param_name(lig, i, false)
            K_T = _reg_param_name(lig, i, true)
            tag == :OnlyT || push!(reg_params_r, K_R)
            if tag == :EqualRT
                push!(reg_params_t_dep, K_T => K_R)
            elseif tag == :NonequalRT || tag == :OnlyT
                push!(reg_params_t_indep, K_T)
            end
        end
    end
    for (k, v) in reg_params_t_dep
        dep_T[k] = v
    end

    merged_dep = merge(dep_R, dep_T)
    merged_indep = (indep_R..., indep_T_list...,
                    reg_params_r..., reg_params_t_indep..., :L)
    return merged_dep, merged_indep
end

# parameters and fitted_params for AllostericEnzymeMechanism are handled
# by the unified _AnyMechanism methods at the top of this file.

# ─── Rate body building helpers ───────────────────────────────────

"""Build the regulatory site partition function expression: 1 + lig/K_lig_reg_i + ...
Skips ligands absent from the given conformation (`:OnlyR` in T-state, `:OnlyT`
in R-state). Uses the R-state K symbol when the ligand tag is `:EqualRT`."""
function _reg_site_expr(m::AllostericEnzymeMechanism, site_idx::Int, T_state::Bool)
    terms = Any[1]
    for lig in regulatory_site_ligands(m, site_idx)
        tag = regulatory_ligand_tag(m, site_idx, lig)
        if T_state
            tag == :OnlyR && continue
        else
            tag == :OnlyT && continue
        end
        use_T_name = T_state && tag in (:NonequalRT, :OnlyT)
        K_sym = _reg_param_name(lig, site_idx, use_T_name)
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
Build R-state and T-state dep-param assignment Exprs.
Returns `(r_assignments::Vector{Expr}, t_assignments::Vector{Expr})`.
Shared by `_build_allosteric_rate_body` and `rate_equation_string`.
"""
function _build_dep_assignments(
    M_type::Type{<:AllostericEnzymeMechanism}, inv_fn,
)
    m = M_type()
    CM = typeof(catalytic_mechanism(m))

    dep_R, indep_R = _dependent_param_exprs(CM)
    dep_R_kd = _apply_kd_inversion(dep_R, CM, inv_fn)
    sorted_deps = sort(collect(dep_R_kd); by=first)

    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)
    rename_T = _T_rename(m)
    T_subs = Dict{Symbol, Symbol}(rename_T)

    # R-state dependent param assignments
    # Set to zero if expression references :OnlyT symbols
    r_assignments = Expr[]
    for (sym, expr_kd) in sorted_deps
        if _expr_references_any(expr_kd, t_only_syms)
            push!(r_assignments, Expr(:(=), sym, 0))
        else
            push!(r_assignments, Expr(:(=), sym, expr_kd))
        end
    end

    t_assignments = Expr[]

    # `:EqualRT` independent catalytic params: p_T = p_R
    # (must come before dep assignments that reference them)
    for p in indep_R
        p ∈ r_only_syms && continue
        haskey(rename_T, p) && continue  # NonequalRT — handled below
        push!(t_assignments, Expr(:(=), _rename_params_T(p), p))
    end

    # `:EqualRT` reg params: K_T_reg = K_R_reg
    for i in eachindex(regulatory_sites(m))
        for lig in regulatory_site_ligands(m, i)
            regulatory_ligand_tag(m, i, lig) == :EqualRT || continue
            push!(t_assignments,
                  Expr(:(=), _reg_param_name(lig, i, true),
                             _reg_param_name(lig, i, false)))
        end
    end

    # T-state dependent param assignments
    # Set to zero if expression references :OnlyR symbols.
    # Skip `:EqualRT` keys (R-state assignment already covers them).
    for (sym, expr_kd) in sorted_deps
        t_sym = get(rename_T, sym, sym)
        t_sym == sym && continue
        if _expr_references_any(expr_kd, r_only_syms)
            push!(t_assignments, Expr(:(=), t_sym, 0))
        else
            push!(t_assignments, Expr(:(=), t_sym,
                substitute_params_expr(expr_kd, T_subs)))
        end
    end

    return r_assignments, t_assignments
end

"""Convert a flat POLY to an Expr using the standard parameter/concentration display."""
function _poly_to_expr(p::POLY, param_syms, conc_syms, inv_set)
    fs = FactoredSigma([poly_one()], [FactoredPoly([p], [1])])
    _factored_sigma_to_expr(fs, param_syms, conc_syms, inv_set)
end

"""
Assemble the MWC numerator and denominator Exprs.
Returns `(full_num, full_den)` where the numerator already includes the
`catalytic_multiplicity` factor.
"""
function _allosteric_num_den_exprs(M_type::Type{<:AllostericEnzymeMechanism})
    m = M_type()
    CM = typeof(catalytic_mechanism(m))
    CatN = catalytic_multiplicity(m)
    RS = regulatory_sites(m)

    num_fs, denom_terms = _raw_symbolic_rate_polys(CM)
    cat_params = Set{Symbol}(_raw_param_symbols(CM()))
    cat_mets = Set{Symbol}(metabolites(CM()))
    binding_Ks_r = Set(_binding_K_symbols(CM))

    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)
    rename_T = _T_rename(m)
    # T-state binding K set: renamed counterparts of R-state binding K's
    # (`:NonequalRT` / `:OnlyT` groups get T-suffixed; `:EqualRT` groups
    # pass through unchanged).
    binding_Ks_t = Set(get(rename_T, K, K) for K in binding_Ks_r)

    # R-state catalytic Exprs.
    # Use factored form when no `:OnlyT` zeroing is needed (preserves nice
    # display). Otherwise expand → zero → flat-poly Expr.
    if isempty(t_only_syms)
        N_R = _factored_sigma_to_expr(num_fs, cat_params, cat_mets, binding_Ks_r)
        Q_R = _denom_terms_to_expr(denom_terms, cat_params, cat_mets, binding_Ks_r)
    else
        num_r_poly = _zero_symbols_in_poly(_expand_factored_sigma(num_fs), t_only_syms)
        den_r_poly = _zero_symbols_in_poly(_expand_to_poly(denom_terms), t_only_syms)
        N_R = _poly_to_expr(num_r_poly, cat_params, cat_mets, binding_Ks_r)
        Q_R = _poly_to_expr(den_r_poly, cat_params, cat_mets, binding_Ks_r)
    end

    # T-state catalytic Exprs.
    # Zero `:OnlyR` symbols at POLY level, then rename `:NonequalRT` / `:OnlyT`
    # symbols to T-suffixed counterparts. `:EqualRT` symbols pass through
    # unchanged (R-state binding) and resolve through dep-param assignments.
    # When the T-state cycle is broken, force N_T = 0: the Cha framework
    # otherwise produces a non-physical reverse flux from products that
    # have nowhere to go.
    if isempty(r_only_syms)
        N_T = substitute_params_expr(
            _factored_sigma_to_expr(num_fs, cat_params, cat_mets, binding_Ks_r),
            rename_T)
        Q_T = substitute_params_expr(
            _denom_terms_to_expr(denom_terms, cat_params, cat_mets, binding_Ks_r),
            rename_T)
    else
        num_t_poly = _rename_symbols(
            _zero_symbols_in_poly(_expand_factored_sigma(num_fs), r_only_syms),
            rename_T)
        den_t_poly = _rename_symbols(
            _zero_symbols_in_poly(_expand_to_poly(denom_terms), r_only_syms),
            rename_T)
        N_T = _t_state_dead(m) ? 0 :
              _poly_to_expr(num_t_poly, cat_params, cat_mets, binding_Ks_t)
        Q_T = _poly_to_expr(den_t_poly, cat_params, cat_mets, binding_Ks_t)
    end

    reg_Q_R = Any[_reg_site_expr(m, i, false) for i in eachindex(RS)]
    reg_Q_T = Any[_reg_site_expr(m, i, true) for i in eachindex(RS)]

    function make_num_term(N, Q, reg_Qs)
        factors = Any[N]
        CatN > 1 && push!(factors, _power_expr(Q, CatN - 1))
        for i in eachindex(RS)
            n_reg = regulatory_site_multiplicity(m, i)
            n_reg == CatN || continue
            push!(factors, _power_expr(reg_Qs[i], n_reg))
        end
        _nest_binary(:*, factors)
    end

    function make_den_term(Q, reg_Qs)
        factors = Any[_power_expr(Q, CatN)]
        for i in eachindex(RS)
            push!(factors, _power_expr(reg_Qs[i],
                                       regulatory_site_multiplicity(m, i)))
        end
        _nest_binary(:*, factors)
    end

    num_R = make_num_term(N_R, Q_R, reg_Q_R)
    den_R = make_den_term(Q_R, reg_Q_R)
    num_T = make_num_term(N_T, Q_T, reg_Q_T)
    den_T = make_den_term(Q_T, reg_Q_T)

    :($(CatN) * ($(num_R) + L * $(num_T))), :($(den_R) + L * $(den_T))
end

"""Build the MWC rate equation body as an Expr block."""
function _build_allosteric_rate_body(M_type::Type{<:AllostericEnzymeMechanism})
    full_num, full_den = _allosteric_num_den_exprs(M_type)
    rate_expr = :(E_total * ($full_num) / ($full_den))

    r_assignments, t_assignments = _build_dep_assignments(M_type, K -> :(inv($K)))

    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(M_type())

    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(mets, :concs),
        r_assignments...,
        t_assignments...,
        :(return $rate_expr))
end

# ─── Rate equation dispatch ───────────────────────────────────────

@generated function rate_equation(
    ::AllostericEnzymeMechanism{CM,CS,RS},
    concs::NamedTuple, params::NamedTuple, ::ReducedMode,
) where {CM,CS,RS}
    _build_allosteric_rate_body(AllostericEnzymeMechanism{CM,CS,RS})
end

# ─── String representation ────────────────────────────────────────

function rate_equation_string(
    ::AllostericEnzymeMechanism{CM,CS,RS}, ::ReducedMode,
) where {CM,CS,RS}
    M = AllostericEnzymeMechanism{CM,CS,RS}
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(M())

    r_assignments, t_assignments = _build_dep_assignments(M, K -> :(1 / $K))
    dep_lines = ["$(a.args[1]) = $(_expr_to_string(a.args[2]))" for a in r_assignments]
    t_dep_lines = ["$(a.args[1]) = $(_expr_to_string(a.args[2]))" for a in t_assignments]

    full_num, full_den = _allosteric_num_den_exprs(M)
    v_line = "v = E_total * ($(_expr_to_string(full_num))) / ($(_expr_to_string(full_den)))"

    join([
        "(; $(join(hw_params, ", "))) = params",
        "(; $(join(mets, ", "))) = concs",
        dep_lines...,
        t_dep_lines...,
        v_line,
    ], "\n")
end

# ─── Structural Identifiability ───────────────────────────────────

@generated function structural_identifiability_deficit(
    ::AllostericEnzymeMechanism{CM,CS,RS},
) where {CM,CS,RS}
    M = AllostericEnzymeMechanism{CM,CS,RS}
    _, indep = _dependent_param_exprs(M)
    n_k = length(indep)
    n_num, n_denom = _count_allosteric_rate_monomials(M)
    n_k - (n_num - 1) - (n_denom - 1)
end

"""
Count distinct concentration monomials in the full allosteric rate numerator
and denominator. Builds R/T catalytic and reg-site polys via the same POLY-level
substitutions used by `rate_equation`, then strips parameter symbols and counts.
Returns `(n_num, n_denom)`.
"""
function _count_allosteric_rate_monomials(M_type::Type{<:AllostericEnzymeMechanism})
    m = M_type()
    CM = typeof(catalytic_mechanism(m))
    CatN = catalytic_multiplicity(m)
    RS = regulatory_sites(m)

    num_fs, denom_terms = _raw_symbolic_rate_polys(CM)
    N_cat_base = _expand_factored_sigma(num_fs)
    Q_cat_base = _expand_to_poly(denom_terms)

    t_only_syms = _onlyT_syms(m)
    r_only_syms = _onlyR_syms(m)
    rename_T = _T_rename(m)

    N_cat_R = _zero_symbols_in_poly(N_cat_base, t_only_syms)
    Q_cat_R = _zero_symbols_in_poly(Q_cat_base, t_only_syms)
    N_cat_T = _rename_symbols(_zero_symbols_in_poly(N_cat_base, r_only_syms), rename_T)
    Q_cat_T = _rename_symbols(_zero_symbols_in_poly(Q_cat_base, r_only_syms), rename_T)

    function reg_q_poly(i, T_state)
        ligs_filtered = Symbol[]
        for lig in regulatory_site_ligands(m, i)
            tag = regulatory_ligand_tag(m, i, lig)
            if T_state
                tag == :OnlyR && continue
            else
                tag == :OnlyT && continue
            end
            push!(ligs_filtered, lig)
        end
        isempty(ligs_filtered) ? poly_one() :
            reduce(poly_add, (poly_add(poly_one(), poly_sym(lig))
                              for lig in ligs_filtered))
    end
    reg_Q_R = POLY[reg_q_poly(i, false) for i in eachindex(RS)]
    reg_Q_T = POLY[reg_q_poly(i, true) for i in eachindex(RS)]

    function num_poly_for_conf(N_cat, Q_cat, reg_Qs, L_factor)
        n_term = poly_mul(N_cat, _poly_power(Q_cat, CatN - 1))
        for i in eachindex(RS)
            n_reg = regulatory_site_multiplicity(m, i)
            n_reg == CatN || continue
            n_term = poly_mul(n_term, _poly_power(reg_Qs[i], n_reg))
        end
        L_factor === nothing ? n_term : poly_mul(poly_sym(L_factor), n_term)
    end

    function den_poly_for_conf(Q_cat, reg_Qs, L_factor)
        d_term = _poly_power(Q_cat, CatN)
        for i in eachindex(RS)
            d_term = poly_mul(d_term, _poly_power(reg_Qs[i],
                                       regulatory_site_multiplicity(m, i)))
        end
        L_factor === nothing ? d_term : poly_mul(poly_sym(L_factor), d_term)
    end

    full_num = poly_add(
        num_poly_for_conf(N_cat_R, Q_cat_R, reg_Q_R, nothing),
        num_poly_for_conf(N_cat_T, Q_cat_T, reg_Q_T, :L),
    )
    full_den = poly_add(
        den_poly_for_conf(Q_cat_R, reg_Q_R, nothing),
        den_poly_for_conf(Q_cat_T, reg_Q_T, :L),
    )

    conc_mono(mono) = sort!(MONO([
        s => e for (s, e) in mono
        if !is_k_parameter(s) && s != :E_total && s != :Keq && s != :L
    ]); by=first)

    n_num = length(unique(conc_mono(k) for (k, _) in full_num))
    n_denom = length(unique(conc_mono(k) for (k, _) in full_den))
    n_num, n_denom
end

