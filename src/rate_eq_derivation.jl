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
    (_raw_param_symbols(EqSteps)..., :E_total)
end

@generated function parameters(::M, ::HaldaneWegscheiderMode) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

# ─── RE Group Helpers ─────────────────────────────────────────────────────────

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
        e_lhs, _ = _split_reaction_side(lhs, enz_set)
        e_rhs, _ = _split_reaction_side(rhs, enz_set)
        union!(findfirst(==(e_lhs), enz_names), findfirst(==(e_rhs), enz_names))
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
    mp = mets -> isempty(mets) ? poly_one() : reduce(poly_mul, poly_sym.(mets))
    alpha_num = Vector{Poly}(fill(poly_one(), N))
    alpha_den = Vector{Poly}(fill(poly_one(), N))

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
    sigma_num = Vector{Poly}(undef, length(groups))
    sigma_den = Vector{Poly}(undef, length(groups))
    for (g, group) in enumerate(groups)
        if length(group) == 1
            sigma_num[g] = poly_one(); sigma_den[g] = poly_one()
        else
            sigma_den[g] = reduce(poly_mul, alpha_den[i] for i in group)
            sigma_num[g] = reduce(poly_add,
                poly_mul(reduce(poly_mul, (alpha_den[j] for j in group if j != i); init=poly_one()), alpha_num[i])
                for i in group)
        end
    end
    alpha_num, alpha_den, sigma_num, sigma_den
end

"""Build rate poly for one SS step direction, clearing alpha denominators within group."""
function _ss_contrib(k_poly, mets, i_form, alpha_num, alpha_den, group)
    r = isempty(mets) ? k_poly : poly_mul(k_poly, reduce(poly_mul, poly_sym.(mets)))
    r = poly_mul(r, alpha_num[i_form])
    for k in group; k == i_form && continue; r = poly_mul(r, alpha_den[k]); end
    r
end

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
    groups, form_to_group = _compute_re_groups(enz_names, enz_set, rxns, eq_steps)
    alpha_num, alpha_den, sigma_num, sigma_den = _compute_alpha(enz_names, enz_set, rxns, eq_steps, groups)
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

        R[g1, g2] = poly_add(R[g1, g2], _ss_contrib(kf_poly, m_lhs, i_form, alpha_num, alpha_den, groups[g1]))
        R[g2, g1] = poly_add(R[g2, g1], _ss_contrib(kr_poly, m_rhs, j_form, alpha_num, alpha_den, groups[g2]))
    end

    # Build Laplacian and cofactor determinants
    L = [i == j ? poly_zero() : poly_neg(R[i,j]) for i in 1:G, j in 1:G]
    for i in 1:G; L[i,i] = reduce(poly_add, R[i,j] for j in 1:G if j != i; init=poly_zero()); end
    D = [begin
        idx = [r for r in 1:G if r != root]
        isempty(idx) ? poly_one() : sym_det(L[idx, idx], G - 1)
    end for root in 1:G]

    denom = reduce(poly_add, poly_mul(sigma_num[g], D[g]) for g in 1:G)

    # Numerator: net flux through any SS step
    ref_name = subs_species[1][1]
    nu_ref = count(s -> s[1] == ref_name, prods_species) - count(s -> s[1] == ref_name, subs_species)

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
    # Classify metabolites into SS/RE step sets, and check ref_name
    ref_in_ss, ref_in_re = false, false
    ss_mets, re_mets = Set{Symbol}(), Set{Symbol}()
    for (idx, (lhs, rhs)) in enumerate(rxns)
        _, m_lhs = _split_reaction_side(lhs, enz_set)
        _, m_rhs = _split_reaction_side(rhs, enz_set)
        target = eq_steps[idx] ? re_mets : ss_mets
        for met in m_lhs; push!(target, met); end
        for met in m_rhs; push!(target, met); end
        met_f = isempty(m_lhs) ? nothing : first(m_lhs)
        met_r = isempty(m_rhs) ? nothing : first(m_rhs)
        if met_f === ref_name || met_r === ref_name
            eq_steps[idx] ? (ref_in_re = true) : (ref_in_ss = true)
        end
    end

    flux = (name) -> _flux_numerator(rxns, eq_steps, enz_names, enz_set, alpha_num, alpha_den,
                                     form_to_group, groups, D, name)
    ref_in_ss && !ref_in_re && return flux(ref_name)

    # Fallback: find alternate metabolite in SS steps
    all_mets = Dict{Symbol, Int}()
    for (name, _) in subs_species; all_mets[name] = get(all_mets, name, 0) - 1; end
    for (name, _) in prods_species; all_mets[name] = get(all_mets, name, 0) + 1; end

    if !isempty(ss_mets)
        ss_only = setdiff(ss_mets, re_mets)
        search = isempty(ss_only) ? ss_mets : ss_only
        alt_name = something(iterate(met for met in search if get(all_mets, met, 0) != 0),
                             (first(ss_mets),))[1]
        nu_alt = get(all_mets, alt_name, 0)
        alt_num = flux(alt_name)
        nu_alt == 0 && return alt_num
        ratio = nu_ref // nu_alt
        ratio == 1 ? (return alt_num) : ratio == -1 ? (return poly_neg(alt_num)) :
            error("Non-unit stoichiometric ratio not supported")
    end

    # All metabolites in RE steps only — use flux through any SS step
    return _flux_numerator(rxns, eq_steps, enz_names, enz_set, alpha_num, alpha_den,
                           form_to_group, groups, D)
end

"""Compute net flux through SS steps. If `met_name` is nothing, return raw flux of first SS step."""
function _flux_numerator(rxns, eq_steps, enz_names, enz_set, alpha_num, alpha_den,
                         form_to_group, groups, D, met_name::Union{Symbol,Nothing}=nothing)
    result = poly_zero()
    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] && continue
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        i_form = findfirst(==(e_lhs), enz_names)
        j_form = findfirst(==(e_rhs), enz_names)
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        rf = _ss_contrib(poly_sym(Symbol("k$(idx)f")), m_lhs, i_form, alpha_num, alpha_den, groups[g1])
        rr = _ss_contrib(poly_sym(Symbol("k$(idx)r")), m_rhs, j_form, alpha_num, alpha_den, groups[g2])
        flux = poly_sub(poly_mul(rf, D[g1]), poly_mul(rr, D[g2]))
        met_name === nothing && return flux
        met_f = isempty(m_lhs) ? nothing : first(m_lhs)
        met_r = isempty(m_rhs) ? nothing : first(m_rhs)
        (met_f !== met_name && met_r !== met_name) && continue
        result = met_f === met_name ? poly_add(result, flux) : poly_sub(result, flux)
    end
    result
end

# ─── Expr generation from Poly ───────────────────────────────────────────────

"""
Identify K symbols for binding RE steps (where K should be Kd, not Ka).
A binding step has: metabolite(s) on LHS, only enzyme forms on RHS.
"""
function _binding_K_symbols(M::Type{<:EnzymeMechanism})
    m = M()
    rxns = reactions(m)
    eq = equilibrium_steps(m)
    enz_set = Set(e[1] for e in enzyme_forms(m))
    [Symbol("K$i") for (i, (lhs, rhs)) in enumerate(rxns)
     if eq[i] && any(s ∉ enz_set for s in lhs) && all(s ∈ enz_set for s in rhs)]
end

"""
Compute the raw rate expression (bare symbols) and sorted parameter/concentration symbols.
Returns `(expr, all_params, sorted_concs)`.
"""
function _raw_rate_expr_and_symbols(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    m = M()
    param_syms = Set{Symbol}(_raw_param_symbols(equilibrium_steps(m)))
    conc_syms = Set{Symbol}(mt[1] for mt in metabolites(m))
    inv_set = Set(_binding_K_symbols(M))
    expr = to_rate_expr(num, den, param_syms, conc_syms, inv_set)
    all_params = _sorted_raw_param_symbols(M)
    sorted_concs = _sorted_conc_symbols(M)
    return expr, all_params, sorted_concs
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
    _build_rate_body(M, RawMode)
end

@generated function rate_equation(m::M, params::NamedTuple, concs::NamedTuple, ::HaldaneWegscheiderMode) where {M <: EnzymeMechanism}
    _build_rate_body(M, HaldaneWegscheiderMode)
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
    m = M()
    ps = Set{Symbol}(_raw_param_symbols(equilibrium_steps(m)))
    cs = Set{Symbol}(mt[1] for mt in metabolites(m))
    inv = Set(_binding_K_symbols(M))
    "v = E_total * ($(string(_poly_to_expr(num, ps, cs, inv)))) / ($(string(_poly_to_expr(den, ps, cs, inv))))"
end

function rate_equation_string(::M, ::RawMode) where {M<:EnzymeMechanism}
    join(["(; $(join(_sorted_raw_param_symbols(M), ", "))) = params",
          "(; $(join(_sorted_conc_symbols(M), ", "))) = concs",
          _format_rate_equation_core(M)], "\n")
end

function rate_equation_string(::M, ::HaldaneWegscheiderMode) where {M<:EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    join(["(; $(join((indep..., :Keq, :E_total), ", "))) = params",
          "(; $(join(_sorted_conc_symbols(M), ", "))) = concs",
          _constraint_expr_strings(M)...,
          _format_rate_equation_core(M)], "\n")
end

# ─── Structural Identifiability ───────────────────────────────────────────────

function _count_rate_monomials(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    mm(mono) = sort([s => e for (s, e) in mono if !is_k_parameter(s) && s != :E_total && s != :Keq])
    length(unique(mm(k) for (k, _) in num)), length(unique(mm(k) for (k, _) in den))
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
