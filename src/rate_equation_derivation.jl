# ─── Parameters API ───────────────────────────────────────────────────────────

"""
    parameters(m::EnzymeMechanism, [mode])

Return the parameter names required for the given mode as a tuple of Symbols.

# Modes
- `HaldaneWegscheider` (default): independent k's + Keq + E_total
- `Raw`: all 2N k's + E_total
"""
function parameters end

parameters(m::EnzymeMechanism) = parameters(m, HaldaneWegscheider)

@generated function parameters(::EnzymeMechanism{Species, Reactions, EqSteps}, ::RawMode) where {Species, Reactions, EqSteps}
    params = Symbol[]
    for (i, is_re) in enumerate(EqSteps)
        if is_re; push!(params, Symbol("K$i"))
        else push!(params, Symbol("k$(i)f")); push!(params, Symbol("k$(i)r")); end
    end
    (Tuple(params)..., :E_total)
end

@generated function parameters(::M, ::HaldaneWegscheiderMode) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

# ─── RE Group Helpers ─────────────────────────────────────────────────────────

"""Compute RE-connected groups via union-find. Returns `(groups, form_to_group)`."""
function _compute_re_groups(enz_names, enz_set, rxns, eq_steps)
    N = length(enz_names)
    parent = collect(1:N)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    function union!(a, b)
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end

    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] || continue
        i = findfirst(==(first(s for s in lhs if s in enz_set)), enz_names)
        j = findfirst(==(first(s for s in rhs if s in enz_set)), enz_names)
        union!(i, j)
    end

    root_to_group = Dict{Int, Int}()
    groups = Vector{Vector{Int}}()
    form_to_group = zeros(Int, N)
    for i in 1:N
        r = find(i)
        if !haskey(root_to_group, r)
            push!(groups, Int[]); root_to_group[r] = length(groups)
        end
        g = root_to_group[r]; push!(groups[g], i); form_to_group[i] = g
    end
    groups, form_to_group
end

"""
Compute alpha factors (relative concentrations within RE groups) as Poly values.
Returns `(alpha, sigma)` where alpha[i] is a (num::Poly, den::Poly) pair
and sigma[g] is a (num::Poly, den::Poly) pair.
"""
function _compute_alpha(enz_names, enz_set, rxns, eq_steps, groups)
    N = length(enz_names)
    one_p = poly_one()
    alpha_num = Vector{Poly}(fill(one_p, N))
    alpha_den = Vector{Poly}(fill(one_p, N))

    for group in groups
        length(group) == 1 && continue
        ref = group[1]
        visited = Set{Int}([ref])
        queue = [ref]

        while !isempty(queue)
            current = popfirst!(queue)
            for (idx, (lhs, rhs)) in enumerate(rxns)
                eq_steps[idx] || continue
                i_form = findfirst(==(first(s for s in lhs if s in enz_set)), enz_names)
                j_form = findfirst(==(first(s for s in rhs if s in enz_set)), enz_names)
                m_lhs = [s for s in lhs if s ∉ enz_set]
                m_rhs = [s for s in rhs if s ∉ enz_set]
                K_poly = poly_sym(Symbol("K$idx"))

                if i_form == current && j_form ∉ visited
                    # Forward: alpha_j = alpha_i * K * prod(met_lhs) / prod(met_rhs)
                    num = poly_mul(poly_mul(alpha_num[current], K_poly), _met_product(m_lhs))
                    den = poly_mul(alpha_den[current], _met_product(m_rhs))
                    alpha_num[j_form] = num; alpha_den[j_form] = den
                    push!(visited, j_form); push!(queue, j_form)
                elseif j_form == current && i_form ∉ visited
                    # Reverse: alpha_i = alpha_j * prod(met_rhs) / (K * prod(met_lhs))
                    num = poly_mul(alpha_num[current], _met_product(m_rhs))
                    den = poly_mul(poly_mul(alpha_den[current], K_poly), _met_product(m_lhs))
                    alpha_num[i_form] = num; alpha_den[i_form] = den
                    push!(visited, i_form); push!(queue, i_form)
                end
            end
        end
    end

    # Compute sigma per group: sum of alpha_i. Clear denominators by using common den.
    sigma_num = Vector{Poly}(undef, length(groups))
    sigma_den = Vector{Poly}(undef, length(groups))
    for (g, group) in enumerate(groups)
        if length(group) == 1
            sigma_num[g] = poly_one(); sigma_den[g] = poly_one()
        else
            # Common denominator = product of all alpha_den in group
            common_den = poly_one()
            for i in group; common_den = poly_mul(common_den, alpha_den[i]); end
            total_num = poly_zero()
            for i in group
                # scale_i = common_den / alpha_den[i] * alpha_num[i]
                other_dens = poly_one()
                for j in group; j == i && continue; other_dens = poly_mul(other_dens, alpha_den[j]); end
                total_num = poly_add(total_num, poly_mul(other_dens, alpha_num[i]))
            end
            sigma_num[g] = total_num; sigma_den[g] = common_den
        end
    end

    alpha_num, alpha_den, sigma_num, sigma_den
end

_met_product(mets) = isempty(mets) ? poly_one() : reduce(poly_mul, poly_sym.(mets))

# ─── Raw Rate Equation Derivation (Unified Cha / King-Altman) ────────────────

"""Build raw numerator and denominator Poly for the rate equation."""
function _raw_symbolic_rate_polys(M::Type{<:EnzymeMechanism})
    m = M()
    subs_species = substrates(m)
    prods_species = products(m)
    enzs = enzyme_forms(m)
    rxns = reactions(m)
    eq_steps = equilibrium_steps(m)

    enz_names = Tuple(e[1] for e in enzs)
    enz_set = Set(enz_names)
    N = length(enzs)

    groups, form_to_group = _compute_re_groups(enz_names, enz_set, rxns, eq_steps)
    alpha_num, alpha_den, sigma_num, sigma_den = _compute_alpha(enz_names, enz_set, rxns, eq_steps, groups)
    G = length(groups)

    # Build rate matrix R[g1,g2] over groups as Poly values
    # To handle alpha = num/den, we multiply through by the group's sigma denominator.
    # R_eff[g1,g2] = R[g1,g2] * sigma_den[g1]  (clearing the denominator)
    R = Matrix{Poly}(undef, G, G)
    for i in 1:G, j in 1:G; R[i,j] = poly_zero(); end

    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] && continue
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        i_form = findfirst(==(e_lhs), enz_names)
        j_form = findfirst(==(e_rhs), enz_names)
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        kf_poly = poly_sym(Symbol("k$(idx)f"))
        kr_poly = poly_sym(Symbol("k$(idx)r"))

        # Forward: g1 → g2, rate = kf * met_lhs * (alpha_num[i] / alpha_den[i])
        # After clearing sigma_den[g1], multiply by the other alpha_dens in the group
        fwd = poly_mul(kf_poly, _met_product(m_lhs))
        fwd = poly_mul(fwd, alpha_num[i_form])
        # Clear alpha_den[i_form] relative to sigma_den[g1] = prod of all alpha_den in group
        for k in groups[g1]; k == i_form && continue; fwd = poly_mul(fwd, alpha_den[k]); end
        R[g1, g2] = poly_add(R[g1, g2], fwd)

        # Reverse: g2 → g1
        rev = poly_mul(kr_poly, _met_product(m_rhs))
        rev = poly_mul(rev, alpha_num[j_form])
        for k in groups[g2]; k == j_form && continue; rev = poly_mul(rev, alpha_den[k]); end
        R[g2, g1] = poly_add(R[g2, g1], rev)
    end

    # Build Laplacian
    L = Matrix{Poly}(undef, G, G)
    for i in 1:G
        diag = poly_zero()
        for j in 1:G
            if i == j; L[i,j] = poly_zero()
            else L[i,j] = poly_neg(R[i,j]); diag = poly_add(diag, R[i,j]); end
        end
        L[i,i] = diag
    end

    # Cofactor determinants
    D = Vector{Poly}(undef, G)
    for root in 1:G
        n_sub = G - 1
        if n_sub == 0; D[root] = poly_one(); continue; end
        rows = [r for r in 1:G if r != root]
        cols = [c for c in 1:G if c != root]
        sub = Matrix{Poly}(undef, n_sub, n_sub)
        for (ri, r) in enumerate(rows), (ci, c) in enumerate(cols)
            sub[ri, ci] = L[r, c]
        end
        D[root] = sym_det(sub, n_sub)
    end

    # Denominator: sum of sigma_num[g] * D[g]  (sigma_den factors are in the R matrix)
    denom = poly_zero()
    for g in 1:G
        denom = poly_add(denom, poly_mul(sigma_num[g], D[g]))
    end

    # Numerator: net flux through any SS step
    ref_name = subs_species[1][1]
    nu_ref = 0
    for (name, _) in subs_species; name == ref_name && (nu_ref -= 1); end
    for (name, _) in prods_species; name == ref_name && (nu_ref += 1); end

    num = _compute_numerator(rxns, eq_steps, enz_names, enz_set, alpha_num, alpha_den,
                             form_to_group, groups, D, ref_name, nu_ref, subs_species, prods_species)

    abs_nu = abs(nu_ref)
    if abs_nu != 1
        denom = poly_mul(denom, poly_const(abs_nu))
    end

    num, denom
end

"""Compute the numerator polynomial by tracking flux of an appropriate metabolite through SS steps."""
function _compute_numerator(rxns, eq_steps, enz_names, enz_set, alpha_num, alpha_den,
                            form_to_group, groups, D, ref_name, nu_ref, subs_species, prods_species)
    # Find a metabolite that appears in SS steps for flux counting
    # Check if ref_name appears in SS and/or RE steps
    ref_in_ss, ref_in_re = false, false
    for (idx, (lhs, rhs)) in enumerate(rxns)
        _, m_lhs = _split_reaction_side(lhs, enz_set)
        _, m_rhs = _split_reaction_side(rhs, enz_set)
        met_f = isempty(m_lhs) ? nothing : first(m_lhs)
        met_r = isempty(m_rhs) ? nothing : first(m_rhs)
        if met_f === ref_name || met_r === ref_name
            eq_steps[idx] ? (ref_in_re = true) : (ref_in_ss = true)
        end
    end

    if ref_in_ss && !ref_in_re
        return _flux_numerator(rxns, eq_steps, enz_names, enz_set, alpha_num, alpha_den,
                               form_to_group, groups, D, ref_name)
    end

    # Fallback: find alternate metabolite in SS steps
    all_mets = Dict{Symbol, Int}()
    for (name, _) in subs_species; all_mets[name] = get(all_mets, name, 0) - 1; end
    for (name, _) in prods_species; all_mets[name] = get(all_mets, name, 0) + 1; end

    ss_mets, re_mets = Set{Symbol}(), Set{Symbol}()
    for (idx, (lhs, rhs)) in enumerate(rxns)
        _, m_lhs = _split_reaction_side(lhs, enz_set)
        _, m_rhs = _split_reaction_side(rhs, enz_set)
        target = eq_steps[idx] ? re_mets : ss_mets
        for met in m_lhs; push!(target, met); end
        for met in m_rhs; push!(target, met); end
    end

    if !isempty(ss_mets)
        ss_only = setdiff(ss_mets, re_mets)
        search = isempty(ss_only) ? ss_mets : ss_only
        alt_name = nothing
        for met in search
            haskey(all_mets, met) && all_mets[met] != 0 && (alt_name = met; break)
        end
        alt_name === nothing && (alt_name = first(ss_mets))
        nu_alt = get(all_mets, alt_name, 0)

        alt_num = _flux_numerator(rxns, eq_steps, enz_names, enz_set, alpha_num, alpha_den,
                                  form_to_group, groups, D, alt_name)
        if nu_alt != 0
            ratio = nu_ref // nu_alt
            if ratio == 1; return alt_num
            elseif ratio == -1; return poly_neg(alt_num)
            else; error("Non-unit stoichiometric ratio not supported"); end
        end
        return alt_num
    end

    # All metabolites in RE steps only — use flux through any SS step
    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] && continue
        e_lhs, _ = _split_reaction_side(lhs, enz_set)
        e_rhs, _ = _split_reaction_side(rhs, enz_set)
        i_form = findfirst(==(e_lhs), enz_names)
        j_form = findfirst(==(e_rhs), enz_names)
        g1, g2 = form_to_group[i_form], form_to_group[j_form]

        rf = poly_sym(Symbol("k$(idx)f"))
        rf = poly_mul(rf, alpha_num[i_form])
        for k in groups[g1]; k == i_form && continue; rf = poly_mul(rf, alpha_den[k]); end

        rr = poly_sym(Symbol("k$(idx)r"))
        rr = poly_mul(rr, alpha_num[j_form])
        for k in groups[g2]; k == j_form && continue; rr = poly_mul(rr, alpha_den[k]); end

        return poly_sub(poly_mul(rf, D[g1]), poly_mul(rr, D[g2]))
    end
    error("No SS steps found for numerator computation")
end

"""Compute net flux of `met_name` through SS steps as a Poly."""
function _flux_numerator(rxns, eq_steps, enz_names, enz_set, alpha_num, alpha_den,
                         form_to_group, groups, D, met_name)
    result = poly_zero()
    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] && continue
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        met_f = isempty(m_lhs) ? nothing : first(m_lhs)
        met_r = isempty(m_rhs) ? nothing : first(m_rhs)
        (met_f !== met_name && met_r !== met_name) && continue

        i_form = findfirst(==(e_lhs), enz_names)
        j_form = findfirst(==(e_rhs), enz_names)
        g1, g2 = form_to_group[i_form], form_to_group[j_form]

        rf = poly_mul(poly_sym(Symbol("k$(idx)f")), _met_product(m_lhs))
        rf = poly_mul(rf, alpha_num[i_form])
        for k in groups[g1]; k == i_form && continue; rf = poly_mul(rf, alpha_den[k]); end

        rr = poly_mul(poly_sym(Symbol("k$(idx)r")), _met_product(m_rhs))
        rr = poly_mul(rr, alpha_num[j_form])
        for k in groups[g2]; k == j_form && continue; rr = poly_mul(rr, alpha_den[k]); end

        flux = poly_sub(poly_mul(rf, D[g1]), poly_mul(rr, D[g2]))
        if met_f === met_name
            result = poly_add(result, flux)
        elseif met_r === met_name
            result = poly_sub(result, flux)
        end
    end
    result
end

# ─── Expr generation from Poly ───────────────────────────────────────────────

function _build_param_conc_sets(M::Type{<:EnzymeMechanism})
    m = M()
    eq_steps = equilibrium_steps(m)
    param_syms = Set{Symbol}()
    for (i, is_re) in enumerate(eq_steps)
        if is_re; push!(param_syms, Symbol("K$i"))
        else push!(param_syms, Symbol("k$(i)f")); push!(param_syms, Symbol("k$(i)r")); end
    end
    conc_syms = Set{Symbol}(mt[1] for mt in metabolites(m))
    param_syms, conc_syms
end

function _raw_symbolic_rate_expr(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    param_syms, conc_syms = _build_param_conc_sets(M)
    to_rate_expr(num, den, param_syms, conc_syms)
end

function _symbolic_rate_expr(M::Type{<:EnzymeMechanism})
    raw_expr = _raw_symbolic_rate_expr(M)
    dep_exprs, _ = _dependent_param_exprs(M)
    isempty(dep_exprs) ? raw_expr : substitute_params_expr(raw_expr, dep_exprs)
end

# ─── Mode-dispatched rate_equation ────────────────────────────────────────────

"""
    rate_equation(m::EnzymeMechanism, params, concs, [mode])

Compute the QSSA steady-state rate. The body is generated at compile time
as a single arithmetic expression with no allocations, loops, or matrix ops.
"""
function rate_equation end

rate_equation(m::EnzymeMechanism, params, concs) = rate_equation(m, params, concs, HaldaneWegscheider)

@generated function rate_equation(m::M, params::NamedTuple, concs::NamedTuple, ::RawMode) where {M <: EnzymeMechanism}
    _raw_symbolic_rate_expr(M)
end

@generated function rate_equation(m::M, params::NamedTuple, concs::NamedTuple, ::HaldaneWegscheiderMode) where {M <: EnzymeMechanism}
    _symbolic_rate_expr(M)
end

# ─── String Representation ────────────────────────────────────────────────────

"""
    rate_equation_string(m::EnzymeMechanism, [mode])

Return a string representation of the rate equation.
"""
function rate_equation_string end

rate_equation_string(m::EnzymeMechanism) = rate_equation_string(m, HaldaneWegscheider)

function _format_rate_equation_core(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    "v = E_total * ($(poly_str(num))) / ($(poly_str(den)))"
end

function rate_equation_string(::M, ::RawMode) where {M<:EnzymeMechanism}
    _format_rate_equation_core(M)
end

function rate_equation_string(::M, ::HaldaneWegscheiderMode) where {M<:EnzymeMechanism}
    eq = _format_rate_equation_core(M)
    constraints = _constraint_expr_strings(M)
    isempty(constraints) ? eq : join(constraints, "\n") * "\n\n" * eq
end

# ─── Structural Identifiability ───────────────────────────────────────────────

function _count_rate_monomials(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    function met_mono(mono::Mono)
        sort([s => e for (s, e) in mono if !is_k_parameter(s) && s != :E_total && s != :Keq])
    end
    n_num = length(unique([met_mono(k) for (k, _) in num]))
    n_den = length(unique([met_mono(k) for (k, _) in den]))
    n_num, n_den
end

"""
    structural_identifiability_deficit(m::EnzymeMechanism) → Int

Number of excess parameters beyond what is structurally identifiable from kinetic data.
"""
@generated function structural_identifiability_deficit(::M) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    n_k = length(indep)
    n_num, n_denom = _count_rate_monomials(M)
    n_k - (n_num - 1) - (n_denom - 1)
end

"""
    is_identifiable(m::EnzymeMechanism) → Bool

Check if all mechanism parameters can be uniquely determined from steady-state kinetic data.
"""
is_identifiable(m::EnzymeMechanism) = structural_identifiability_deficit(m) <= 0
