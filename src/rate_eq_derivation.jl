# ─── Parameters API ─────────────────────────────────────────

"""
    parameters(m::EnzymeMechanism, [mode])

Return the parameter names required for the given mode as a tuple of Symbols.

# Modes
- `Reduced` (default): independent k's + Keq + E_total
- `Full`: all 2N k's + E_total
"""
function parameters end

parameters(m::EnzymeMechanism) = parameters(m, Reduced)

@generated function parameters(
    ::EnzymeMechanism{Species, Reactions, EqSteps, PC},
    ::FullMode,
) where {Species, Reactions, EqSteps, PC}
    constrained = Set(c[1] for c in PC)
    Tuple(p for p in (_raw_param_symbols(EqSteps)..., :E_total) if p ∉ constrained)
end

@generated function parameters(
    ::EnzymeMechanism{Sp, Rx, Eq, PC}, ::ReducedMode,
) where {Sp, Rx, Eq, PC}
    M = EnzymeMechanism{Sp, Rx, Eq, PC}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

"""Independent rate constant names for fitting (excludes Keq, E_total)."""
@generated function fitted_params(
    ::EnzymeMechanism{Sp, Rx, Eq, PC},
) where {Sp, Rx, Eq, PC}
    M = EnzymeMechanism{Sp, Rx, Eq, PC}
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
            sigma_num[g] = poly_one(); sigma_den[g] = poly_one()
        else
            sigma_den[g] = reduce(poly_mul, alpha_den[i] for i in group)
            sigma_num[g] = reduce(
                poly_add,
                (
                    poly_mul(
                        reduce(
                            poly_mul,
                            (alpha_den[j] for j in group if j != i);
                            init=poly_one(),
                        ),
                        alpha_num[i],
                    )
                    for i in group
                ),
            )
        end
    end
    alpha_num, alpha_den, sigma_num, sigma_den
end

# ─── Sigma Factoring Detection ─────────────────────────────────

"""
Split an RE group into conformational sub-groups connected by binding
steps only. Returns `(sub_groups, sub_group_coeffs)` where each
sub-group is a `Vector{Int}` of form indices and `sub_group_coeffs[k]`
is the POLY coefficient for sub-group k relative to sub-group 1.
"""
function _split_conformational_subgroups(
    group, enz_names, enz_set, rxns, eq_steps,
)
    group_set = Set(group)
    N = length(group)
    global_to_local = Dict(g => i for (i, g) in enumerate(group))

    parent = collect(1:N)
    function find(x)
        while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end
        x
    end
    uf_union!(a, b) = let ra = find(a), rb = find(b)
        ra != rb && (parent[ra] = rb)
    end

    iso_steps = Int[]
    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] || continue
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        i_f = findfirst(==(e_lhs), enz_names)
        j_f = findfirst(==(e_rhs), enz_names)
        (i_f ∈ group_set && j_f ∈ group_set) || continue
        has_met_lhs = !isempty(m_lhs)
        has_met_rhs = !isempty(m_rhs)
        if has_met_lhs ⊻ has_met_rhs
            uf_union!(global_to_local[i_f], global_to_local[j_f])
        elseif !has_met_lhs && !has_met_rhs
            push!(iso_steps, idx)
        end
    end

    comp_map = Dict{Int, Int}()
    sub_groups = Vector{Vector{Int}}()
    for (local_i, global_i) in enumerate(group)
        r = find(local_i)
        if !haskey(comp_map, r)
            push!(sub_groups, Int[])
            comp_map[r] = length(sub_groups)
        end
        push!(sub_groups[comp_map[r]], global_i)
    end

    length(sub_groups) == 1 && return sub_groups, [poly_one()]

    form_to_sg = Dict{Int, Int}()
    for (sg_idx, sg) in enumerate(sub_groups)
        for f in sg; form_to_sg[f] = sg_idx; end
    end

    n_sg = length(sub_groups)
    coeffs = Vector{Union{POLY, Nothing}}(fill(nothing, n_sg))
    coeffs[1] = poly_one()
    visited_sg = Set{Int}([1])
    queue = [1]

    while !isempty(queue)
        cur_sg = popfirst!(queue)
        for idx in iso_steps
            lhs, rhs = rxns[idx]
            e_lhs, _ = _split_reaction_side(lhs, enz_set)
            e_rhs, _ = _split_reaction_side(rhs, enz_set)
            i_f = findfirst(==(e_lhs), enz_names)
            j_f = findfirst(==(e_rhs), enz_names)
            (i_f ∈ group_set && j_f ∈ group_set) || continue
            sg_i = form_to_sg[i_f]
            sg_j = form_to_sg[j_f]
            K = poly_sym(Symbol("K$idx"))
            if sg_i == cur_sg && sg_j ∉ visited_sg
                coeffs[sg_j] = poly_mul(coeffs[sg_i]::POLY, K)
                push!(visited_sg, sg_j); push!(queue, sg_j)
            elseif sg_j == cur_sg && sg_i ∉ visited_sg
                # Reverse traversal requires 1/K — not representable as POLY.
                # Return nothing to signal factoring failure.
                return sub_groups, nothing
            end
        end
    end

    any(isnothing, coeffs) && return sub_groups, nothing
    return sub_groups, POLY[c for c in coeffs]
end

"""
Compute the metabolite binding state for each form in a sub-group via
BFS through binding RE steps. Returns `Dict{Int, Dict{Symbol, Int}}`
mapping form index → {metabolite => count}.
"""
function _metabolite_binding_states(
    subgroup, enz_names, enz_set, rxns, eq_steps,
)
    sg_set = Set(subgroup)
    ref = first(subgroup)
    states = Dict{Int, Dict{Symbol, Int}}(ref => Dict{Symbol, Int}())
    queue = [ref]

    while !isempty(queue)
        cur = popfirst!(queue)
        for (idx, (lhs, rhs)) in enumerate(rxns)
            eq_steps[idx] || continue
            e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
            e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
            i_f = findfirst(==(e_lhs), enz_names)
            j_f = findfirst(==(e_rhs), enz_names)
            (i_f ∈ sg_set && j_f ∈ sg_set) || continue
            has_met_lhs = !isempty(m_lhs)
            has_met_rhs = !isempty(m_rhs)
            (has_met_lhs ⊻ has_met_rhs) || continue

            met = has_met_lhs ? first(m_lhs) : first(m_rhs)
            # Determine less-bound and more-bound forms
            if has_met_lhs
                less, more = i_f, j_f  # [E, M] ⇌ [EM]
            else
                less, more = j_f, i_f  # [EM] ⇌ [E, M]
            end

            if cur == less && !haskey(states, more)
                s = copy(states[cur])
                s[met] = get(s, met, 0) + 1
                states[more] = s
                push!(queue, more)
            elseif cur == more && !haskey(states, less)
                s = copy(states[cur])
                s[met] = get(s, met, 0) - 1
                s[met] == 0 && delete!(s, met)
                states[less] = s
                push!(queue, less)
            end
        end
    end
    states
end

"""
Partition metabolites into independent binding sites. Two metabolites
share a site if they never co-occur in any form. Returns a vector of
vectors of metabolite symbols, one per site.
"""
function _partition_metabolite_sites(
    binding_states::Dict{Int, Dict{Symbol, Int}},
)
    all_mets = Set{Symbol}()
    for (_, bs) in binding_states
        for m in keys(bs); push!(all_mets, m); end
    end
    mets = sort!(collect(all_mets))
    isempty(mets) && return Vector{Symbol}[]

    # Build co-occurrence: do two metabolites ever appear together?
    cooccur = Dict{Tuple{Symbol,Symbol}, Bool}()
    for m1 in mets, m2 in mets
        m1 >= m2 && continue
        cooccur[(m1, m2)] = false
    end
    for (_, bs) in binding_states
        present = [m for m in mets if get(bs, m, 0) > 0]
        for i in eachindex(present), j in (i+1):length(present)
            a, b = minmax(present[i], present[j])
            cooccur[(a, b)] = true
        end
    end

    # Union-find: merge metabolites that NEVER co-occur
    # (they compete for the same site)
    met_idx = Dict(m => i for (i, m) in enumerate(mets))
    uf_parent = collect(1:length(mets))
    function find(x)
        while uf_parent[x] != x
            uf_parent[x] = uf_parent[uf_parent[x]]; x = uf_parent[x]
        end
        x
    end
    uf_union!(a, b) = let ra = find(a), rb = find(b)
        ra != rb && (uf_parent[ra] = rb)
    end

    for m1 in mets, m2 in mets
        m1 >= m2 && continue
        # Share a site if they never co-occur AND there's substitution evidence
        if !cooccur[(minmax(m1, m2))]
            # Check substitution: exists a form with m1 that has a neighbor with m2
            has_sub = false
            for (_, bs) in binding_states
                if get(bs, m1, 0) > 0 && get(bs, m2, 0) == 0
                    # Check if a related form has m2 instead of m1
                    for (_, bs2) in binding_states
                        if get(bs2, m2, 0) > 0 && get(bs2, m1, 0) == 0
                            # Same total binding count on this site
                            other_same = true
                            for m3 in mets
                                m3 == m1 && continue
                                m3 == m2 && continue
                                get(bs, m3, 0) != get(bs2, m3, 0) && (other_same = false; break)
                            end
                            if other_same && (
                                get(bs, m1, 0) == get(bs2, m2, 0)
                            )
                                has_sub = true; break
                            end
                        end
                    end
                    has_sub && break
                end
            end
            has_sub && uf_union!(met_idx[m1], met_idx[m2])
        end
    end

    groups = Dict{Int, Vector{Symbol}}()
    for (i, m) in enumerate(mets)
        r = find(i)
        haskey(groups, r) || (groups[r] = Symbol[])
        push!(groups[r], m)
    end
    collect(values(groups))
end

"""
Check if the forms in a sub-group equal the Cartesian product of
per-site binding states. Returns `(is_product, per_site_states)`.
"""
function _check_cartesian_product(
    subgroup, binding_states, sites,
)
    # Collect per-site state tuples (occupancy count per metabolite at that site)
    per_site_states = Vector{Vector{Dict{Symbol, Int}}}()
    for site_mets in sites
        states_set = Set{Dict{Symbol, Int}}()
        for f in subgroup
            bs = binding_states[f]
            site_state = Dict{Symbol, Int}(
                m => get(bs, m, 0) for m in site_mets if get(bs, m, 0) > 0
            )
            push!(states_set, site_state)
        end
        push!(per_site_states, sort!(collect(states_set); by=repr))
    end

    # Generate Cartesian product and check bijection with sub-group forms
    expected = Set{Dict{Symbol, Int}}()
    function _cart_product(site_idx, acc)
        if site_idx > length(per_site_states)
            push!(expected, copy(acc))
            return
        end
        for ss in per_site_states[site_idx]
            merged = copy(acc)
            for (m, c) in ss
                merged[m] = get(merged, m, 0) + c
            end
            _cart_product(site_idx + 1, merged)
        end
    end
    _cart_product(1, Dict{Symbol, Int}())

    actual = Set{Dict{Symbol, Int}}()
    for f in subgroup
        bs = binding_states[f]
        push!(actual, Dict(m => c for (m, c) in bs if c > 0))
    end

    return expected == actual, per_site_states
end

"""
Check K consistency: all binding RE steps for the same metabolite
within a sub-group must resolve to the same canonical K (after
constraint resolution). Returns `(consistent, met_to_K)`.

When `sigma_K_syms` is provided, steps whose canonical K is not in
the sigma are skipped — they correspond to BFS paths not taken.
"""
function _check_k_consistency(
    subgroup, enz_names, enz_set, rxns, eq_steps,
    binding_states, sites, constraints;
    sigma_K_syms::Union{Nothing, Set{Symbol}}=nothing,
)
    sg_set = Set(subgroup)
    # Build constraint resolution map: chase target → replacement
    resolve = Dict{Symbol, Symbol}()
    for (target, coeff, factors) in constraints
        # Only simple constraints (coeff=1, single factor with exp=1) are compatible
        coeff != 1 && return false, Dict{Symbol, Symbol}()
        length(factors) != 1 && return false, Dict{Symbol, Symbol}()
        _, exp = factors[1]
        exp != 1 && return false, Dict{Symbol, Symbol}()
        resolve[target] = factors[1][1]
    end
    # Chase to canonical form
    function canonical(sym)
        visited = Set{Symbol}()
        while haskey(resolve, sym) && sym ∉ visited
            push!(visited, sym)
            sym = resolve[sym]
        end
        sym
    end

    # Map metabolite → canonical K
    met_to_K = Dict{Symbol, Symbol}()
    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] || continue
        e_lhs, m_lhs = _split_reaction_side(lhs, enz_set)
        e_rhs, m_rhs = _split_reaction_side(rhs, enz_set)
        i_f = findfirst(==(e_lhs), enz_names)
        j_f = findfirst(==(e_rhs), enz_names)
        (i_f ∈ sg_set && j_f ∈ sg_set) || continue
        has_met_lhs = !isempty(m_lhs)
        has_met_rhs = !isempty(m_rhs)
        (has_met_lhs ⊻ has_met_rhs) || continue

        met = has_met_lhs ? first(m_lhs) : first(m_rhs)
        K_sym = canonical(Symbol("K$idx"))

        # Skip steps whose canonical K doesn't appear in the sigma
        # (BFS used a different path for this metabolite)
        if sigma_K_syms !== nothing && K_sym ∉ sigma_K_syms
            continue
        end

        if haskey(met_to_K, met)
            met_to_K[met] != K_sym && return false, Dict{Symbol, Symbol}()
        else
            met_to_K[met] = K_sym
        end
    end
    true, met_to_K
end

"""
Build per-site POLY factors for a sub-group. Returns a `FactoredPoly`
with one factor per site and detected exponents for repeated factors.

When `alpha_ratio` is provided, each factor is the direct sum of
alpha_ratio values (normalized, may have negative K exponents).
Otherwise falls back to cleared-denominator computation.
"""
function _build_site_factors(
    subgroup, enz_names, enz_set, rxns, eq_steps,
    alpha_num, alpha_den, sites, met_to_K, binding_states;
    alpha_ratio=nothing,
)
    factors = POLY[]
    for site_mets in sites
        site_met_set = Set(site_mets)
        # Collect forms relevant to this site:
        # reference (all sites empty) + forms with ONLY this site occupied
        site_forms = Int[]
        for f in subgroup
            bs = binding_states[f]
            other = any(c > 0 for (m, c) in bs if m ∉ site_met_set)
            other && continue
            push!(site_forms, f)
        end

        if alpha_ratio !== nothing
            factor = reduce(poly_add, (alpha_ratio[fi] for fi in site_forms))
        else
            factor = poly_zero()
            for fi in site_forms
                term = alpha_num[fi]
                for fj in site_forms
                    fi == fj && continue
                    term = poly_mul(term, alpha_den[fj])
                end
                factor = poly_add(factor, term)
            end
        end
        push!(factors, factor)
    end

    # Detect repeated factors: merge identical POLY values
    unique_factors = POLY[]
    exponents = Int[]
    for f in factors
        idx = findfirst(==(f), unique_factors)
        if idx !== nothing
            exponents[idx] += 1
        else
            push!(unique_factors, f)
            push!(exponents, 1)
        end
    end

    FactoredPoly(unique_factors, exponents)
end

"""
Try to factor sigma for an RE group into a `FactoredSigma`.
Returns `nothing` if factoring is not possible.

When `alpha_ratio` is provided (normalized alpha = alpha_num / alpha_den),
produces normalized factors without denominator-clearing artifacts.
"""
function _try_factor_sigma(
    group, enz_names, enz_set, rxns, eq_steps,
    alpha_num, alpha_den, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
    alpha_ratio::Union{Nothing, Vector{POLY}}=nothing,
)
    use_ratio = alpha_ratio !== nothing

    # Step 1: Split into conformational sub-groups
    sub_groups, coeffs = _split_conformational_subgroups(
        group, enz_names, enz_set, rxns, eq_steps,
    )
    coeffs === nothing && return nothing

    factored_products = FactoredPoly[]
    final_coeffs = POLY[]

    for (sg_idx, sg) in enumerate(sub_groups)
        # Sub-group coefficient: iso-step weight.
        # Without alpha_ratio, also multiply by other-sub-group denominators.
        if use_ratio
            sg_coeff = coeffs[sg_idx]
        else
            other_den = poly_one()
            sg_set = Set(sg)
            for j in group
                j ∈ sg_set && continue
                other_den = poly_mul(other_den, alpha_den[j])
            end
            sg_coeff = poly_mul(coeffs[sg_idx], other_den)
        end

        # Step 2: Compute binding states
        states = _metabolite_binding_states(sg, enz_names, enz_set, rxns, eq_steps)
        length(states) != length(sg) && return nothing

        # Step 3: Partition metabolite sites
        sites = _partition_metabolite_sites(states)
        isempty(sites) && length(sg) == 1 && begin
            # Single form, no metabolites → trivial factor
            push!(factored_products, FactoredPoly([poly_one()], [1]))
            push!(final_coeffs, sg_coeff)
            continue
        end

        # Step 4: Check Cartesian product
        is_product, per_site_states = _check_cartesian_product(sg, states, sites)
        is_product || return nothing

        # Step 5: Check K consistency
        # Compute constrained sub-group sigma to extract which K symbols appear
        if use_ratio
            sg_sigma = reduce(poly_add, (alpha_ratio[i] for i in sg))
        else
            sg_sigma = reduce(
                poly_add,
                (poly_mul(
                    alpha_num[i],
                    reduce(
                        poly_mul,
                        (alpha_den[j] for j in sg if j != i);
                        init=poly_one(),
                    ),
                ) for i in sg),
            )
        end
        sg_sigma = _apply_param_constraints(sg_sigma, constraints; binding_Ks)
        sigma_K_syms = Set{Symbol}(
            s for (mono, _) in sg_sigma for (s, _) in mono
            if startswith(string(s), "K") && length(string(s)) > 1 &&
               isdigit(string(s)[2])
        )
        consistent, met_to_K = _check_k_consistency(
            sg, enz_names, enz_set, rxns, eq_steps, states, sites,
            constraints; sigma_K_syms,
        )
        consistent || return nothing

        # Step 6: Build per-site factors
        fp = _build_site_factors(
            sg, enz_names, enz_set, rxns, eq_steps,
            alpha_num, alpha_den, sites, met_to_K, states;
            alpha_ratio,
        )

        # Divide out ref_form alpha from factors to avoid double-counting.
        ref_alpha = use_ratio ? alpha_ratio[first(sg)] : alpha_num[first(sg)]
        if ref_alpha != poly_one()
            fp = FactoredPoly(
                [_poly_div_mono(f, ref_alpha) for f in fp.factors],
                copy(fp.exponents),
            )
        end

        push!(factored_products, fp)
        push!(final_coeffs, sg_coeff)
    end

    # Verify: expanded factored sigma must match the expected sigma
    fs = FactoredSigma(final_coeffs, factored_products)
    est = _estimate_expanded_term_count([DenomTerm(fs, poly_one())])
    if est <= MAX_RATE_EQUATION_TERMS
        if use_ratio
            sigma_expected = reduce(
                poly_add, (alpha_ratio[i] for i in group),
            )
        else
            sigma_expected = reduce(
                poly_add,
                (poly_mul(
                    alpha_num[i],
                    reduce(
                        poly_mul,
                        (alpha_den[j] for j in group if j != i);
                        init=poly_one(),
                    ),
                ) for i in group),
            )
        end
        sigma_expected = _apply_param_constraints(
            sigma_expected, constraints; binding_Ks,
        )
        sigma_actual = _expand_factored_sigma(fs)
        sigma_expected != sigma_actual && return nothing
    end
    fs
end

# ─── Algebraic Regulator Factoring (Fallback) ─────────────────

"""Find the canonical K symbol for a metabolite's binding step."""
function _find_met_binding_K(met, rxns, eq_steps, enz_set, constraints)
    # Build constraint resolution chain
    resolve = Dict{Symbol, Symbol}()
    for (target, coeff, factors) in constraints
        coeff == 1 && length(factors) == 1 && factors[1][2] == 1 || continue
        resolve[target] = factors[1][1]
    end
    function canonical(sym)
        visited = Set{Symbol}()
        while haskey(resolve, sym) && sym ∉ visited
            push!(visited, sym)
            sym = resolve[sym]
        end
        sym
    end
    for (idx, (lhs, rhs)) in enumerate(rxns)
        eq_steps[idx] || continue
        _, m_lhs = _split_reaction_side(lhs, enz_set)
        _, m_rhs = _split_reaction_side(rhs, enz_set)
        has_met_lhs = met ∈ m_lhs
        has_met_rhs = met ∈ m_rhs
        (has_met_lhs ⊻ has_met_rhs) || continue
        return canonical(Symbol("K$idx"))
    end
    nothing
end

"""
Count monomial "shapes" shared between `p1` and `p2`, where shapes are
monomials with rate-constant symbols (k-prefixed) removed. High overlap
means the two polynomials have parallel structure differing only in rate
constants — a sign of clean factoring.
"""
function _mono_shape_overlap(p1::POLY, p2::POLY)
    _shape(mono) = sort!([s => e for (s, e) in mono if !startswith(string(s), "k")])
    shapes1 = Set(_shape(mono) for (mono, _) in p1)
    shapes2 = Set(_shape(mono) for (mono, _) in p2)
    length(shapes1 ∩ shapes2)
end

"""
Try to factor `sigma` by metabolite `met` with binding K symbol `K_R`.
Returns `(base_R, without_R, is_subset, remainder)` or `nothing`.
"""
function _try_factor_by_met(sigma::POLY, met::Symbol, K_R::Symbol)
    with_R = POLY()
    without_R = POLY()
    for (mono, coeff) in sigma
        if any(s == met for (s, _) in mono)
            with_R[mono] = coeff
        else
            without_R[mono] = coeff
        end
    end
    isempty(with_R) && return nothing

    # Divide each term in with_R by K_R * met
    K_R_met = _mono(K_R => 1, met => 1)
    base_R = POLY()
    for (mono, coeff) in with_R
        r_exp = 0; k_exp = 0
        for (s, e) in mono
            s == met && (r_exp = e)
            s == K_R && (k_exp = e)
        end
        (r_exp >= 1 && k_exp >= 1) || return nothing
        base_R[_mono_div(mono, K_R_met)] = coeff
    end

    is_subset = all(
        haskey(without_R, m) && without_R[m] == c for (m, c) in base_R
    )
    remainder = is_subset ? poly_sub(without_R, base_R) : nothing
    (base_R, without_R, is_subset, remainder)
end

"""
Try algebraic polynomial factoring on a constrained sigma.
Returns a `FactoredSigma` or `nothing`.
"""
function _try_algebraic_factor_sigma(
    sigma::POLY, rxns, eq_steps, enz_set, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
)
    # Collect metabolites present in sigma
    K_syms = Set{Symbol}(
        s for (mono, _) in sigma for (s, _) in mono
        if startswith(string(s), "K") && length(string(s)) > 1 &&
           isdigit(string(s)[2])
    )
    k_syms = Set{Symbol}(
        s for (mono, _) in sigma for (s, _) in mono
        if startswith(string(s), "k")
    )
    met_syms = Symbol[
        s for s in Set(
            s for (mono, _) in sigma for (s, _) in mono
        ) if s ∉ K_syms && s ∉ k_syms && s != :Keq && s != :E_total
    ]
    sort!(met_syms; by=string)

    # Try each metabolite and pick the best factoring.
    # Score tuple (lower is better):
    #   1. unfactored terms — fewer remainder/without_R terms is better
    #   2. -product_terms — more terms in product factor is better
    #   3. -shape_overlap — prefer extractions where without_R and base_R
    #      share the same monomial shapes (ignoring rate constants)
    best = nothing
    best_score = (length(sigma), 0, 0)
    for met in met_syms
        K_R = _find_met_binding_K(met, rxns, eq_steps, enz_set, constraints)
        K_R === nothing && continue
        K_R ∈ K_syms || continue
        result = _try_factor_by_met(sigma, met, K_R)
        result === nothing && continue
        base_R, without_R, is_subset, remainder = result
        unfactored = is_subset ? remainder::POLY : without_R
        overlap = _mono_shape_overlap(unfactored, base_R)

        coeffs = POLY[]
        products = FactoredPoly[]
        if is_subset
            # sigma = remainder + base_R * (1 + K_R*met)
            one_plus_KR = poly_add(
                poly_one(), POLY(_mono(K_R => 1, met => 1) => 1),
            )
            if !isempty(remainder)
                push!(coeffs, remainder)
                push!(products, FactoredPoly([poly_one()], [1]))
            end
            push!(coeffs, poly_one())
            push!(products, FactoredPoly([base_R, one_plus_KR], [1, 1]))
        else
            # sigma = without_R + K_R*met * base_R
            K_R_met = POLY(_mono(K_R => 1, met => 1) => Rational{Int}(1))
            if !isempty(without_R)
                push!(coeffs, without_R)
                push!(products, FactoredPoly([poly_one()], [1]))
            end
            push!(coeffs, K_R_met)
            push!(products, FactoredPoly([base_R], [1]))
        end
        score = (length(unfactored), -length(base_R), -overlap)
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
    # Check constant term is 1
    const_mono = MONO()
    haskey(p, const_mono) && p[const_mono] == 1 || return nothing

    # Find candidate squared terms: coefficient 1, all even exponents
    sq_terms = MONO[]
    for (mono, coeff) in p
        isempty(mono) && continue
        coeff == 1 || continue
        all(iseven(e) for (_, e) in mono) || continue
        push!(sq_terms, mono)
    end
    isempty(sq_terms) && return nothing

    # Form candidate Q = 1 + sum(sqrt_terms)
    Q = POLY(MONO() => Rational{Int}(1))
    for mono in sq_terms
        sqrt_mono = MONO([s => div(e, 2) for (s, e) in mono])
        Q[sqrt_mono] = 1
    end

    # Try Q^2 first
    Q2 = poly_mul(Q, Q)
    if Q2 == p
        return FactoredSigma(
            [poly_one()], [FactoredPoly([Q], [2])],
        )
    end

    # Try Q^3
    Q3 = poly_mul(Q2, Q)
    if Q3 == p
        return FactoredSigma(
            [poly_one()], [FactoredPoly([Q], [3])],
        )
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

"""Sum polynomial entries from `polys` at indices `group`."""
function _sum_group_polys(polys::Vector{POLY}, group)
    reduce(poly_add, (polys[i] for i in group))
end

# ─── Raw Rate Equation Derivation (Unified Cha / King-Altman) ───

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

    # Compute normalized alpha_ratio for G=1 when sigma_den is non-trivial
    pc = param_constraints(m)
    binding_Ks = Set{Symbol}(
        Symbol("K$i") for (i, (lhs, rhs)) in enumerate(rxns)
        if eq_steps[i] && any(s ∉ enz_set for s in lhs) &&
           all(s ∈ enz_set for s in rhs)
    )
    normalize = G == 1 && sigma_den[1] != poly_one()
    ar = if normalize
        N = length(enz_names)
        ratio = Vector{POLY}(undef, N)
        for i in 1:N
            ratio[i] = _poly_div_mono(alpha_num[i], alpha_den[i])
        end
        ratio
    else
        nothing
    end

    # Try to factor sigma for each RE group; fall back to unfactored
    denom_terms = DenomTerm[]
    for g in 1:G
        fs = _try_factor_sigma(
            groups[g], enz_names, enz_set, rxns, eq_steps,
            alpha_num, alpha_den, pc; binding_Ks, alpha_ratio=ar,
        )
        if fs !== nothing
            push!(denom_terms, DenomTerm(fs, D[g]))
        else
            # Algebraic fallback: factor constrained sigma by metabolite
            raw_sigma = if normalize
                _sum_group_polys(ar::Vector{POLY}, groups[g])
            else
                sigma_num[g]
            end
            csigma = _apply_param_constraints(raw_sigma, pc; binding_Ks)
            # Try polynomial power detection first (e.g., (1+S/K)^2)
            pfs = _try_poly_power(csigma)
            if pfs !== nothing
                push!(denom_terms, DenomTerm(pfs, D[g]))
            else
                afs = _try_algebraic_factor_sigma(
                    csigma, rxns, eq_steps, enz_set, pc; binding_Ks,
                )
                if afs !== nothing
                    push!(denom_terms, DenomTerm(afs, D[g]))
                else
                    push!(denom_terms, unfactored_denom_term(raw_sigma, D[g]))
                end
            end
        end
    end

    # Numerator: net flux through any SS step
    ref_name = subs_species[1][1]
    nu_ref = (count(s -> s[1] == ref_name, prods_species) -
              count(s -> s[1] == ref_name, subs_species))

    num = _compute_numerator(
        rxns, eq_steps, enz_names, enz_set,
        alpha_num, alpha_den, form_to_group, groups,
        D, ref_name, nu_ref, subs_species, prods_species,
    )

    # Normalize numerator by sigma_den when G=1
    if normalize
        num = _poly_div_mono(num, sigma_den[1])
    end

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

    n_terms = length(num) + _estimate_expanded_term_count(denom_terms)
    if n_terms > MAX_RATE_EQUATION_TERMS
        error(
            "Rate equation for this mechanism has $n_terms polynomial " *
            "terms (limit: $MAX_RATE_EQUATION_TERMS). Equations this " *
            "large take a very long time to compile and are unlikely " *
            "to be practically useful for parameter fitting.",
        )
    end

    num, denom_terms
end

"""
Compute the numerator polynomial by tracking flux of an
appropriate metabolite through SS steps.
"""
function _compute_numerator(
    rxns, eq_steps, enz_names, enz_set,
    alpha_num, alpha_den, form_to_group, groups,
    D, ref_name, nu_ref, subs_species, prods_species,
)
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

    flux = (name) -> _flux_numerator(
        rxns, eq_steps, enz_names, enz_set,
        alpha_num, alpha_den,
        form_to_group, groups, D, name,
    )
    ref_in_ss && !ref_in_re && return flux(ref_name)

    # Fallback: find alternate metabolite in SS steps
    all_mets = Dict{Symbol, Int}()
    for (name, _) in subs_species; all_mets[name] = get(all_mets, name, 0) - 1; end
    for (name, _) in prods_species; all_mets[name] = get(all_mets, name, 0) + 1; end

    if !isempty(ss_mets)
        ss_only = setdiff(ss_mets, re_mets)
        search = isempty(ss_only) ? ss_mets : ss_only
        alt_name = something(
            iterate(
                met for met in search
                if get(all_mets, met, 0) != 0
            ),
            (first(ss_mets),),
        )[1]
        nu_alt = get(all_mets, alt_name, 0)
        alt_num = flux(alt_name)
        nu_alt == 0 && return alt_num
        ratio = nu_ref // nu_alt
        if ratio == 1
            return alt_num
        elseif ratio == -1
            return poly_neg(alt_num)
        else
            error("Non-unit stoichiometric ratio not supported")
        end
    end

    # All metabolites in RE steps only — flux through any SS step
    return _flux_numerator(
        rxns, eq_steps, enz_names, enz_set,
        alpha_num, alpha_den,
        form_to_group, groups, D,
    )
end

"""
Compute net flux through SS steps. If `met_name` is nothing,
sum raw flux of all SS steps (G=1 case).
"""
function _flux_numerator(
    rxns, eq_steps, enz_names, enz_set,
    alpha_num, alpha_den, form_to_group, groups,
    D, met_name::Union{Symbol,Nothing}=nothing,
)
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
            continue
        end
        met_f = isempty(m_lhs) ? nothing : first(m_lhs)
        met_r = isempty(m_rhs) ? nothing : first(m_rhs)
        (met_f !== met_name && met_r !== met_name) && continue
        result = met_f === met_name ? poly_add(result, flux) : poly_sub(result, flux)
    end
    result
end

# ─── Expr generation from POLY ──────────────────────────────

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
    sorted_concs = _sorted_conc_symbols(M)
    return expr, all_params, sorted_concs
end

# ─── Mode-dispatched rate_equation ────────────────────────────

"""
    rate_equation(m::EnzymeMechanism, concs, params, [mode])

Compute the QSSA steady-state rate. The body is generated at compile time
as a single arithmetic expression with no allocations, loops, or matrix ops.
"""
function rate_equation end

function rate_equation(m::EnzymeMechanism, concs, params)
    rate_equation(m, concs, params, Reduced)
end

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
    rate_equation_string(m::EnzymeMechanism, [mode])

Return a string representation of the rate equation.
"""
function rate_equation_string end

rate_equation_string(m::EnzymeMechanism) = rate_equation_string(m, Reduced)

function _format_rate_equation_core(M::Type{<:EnzymeMechanism})
    num, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    constrained = Set(c[1] for c in param_constraints(m))
    raw_ps = _raw_param_symbols(equilibrium_steps(m))
    ps = Set{Symbol}(p for p in raw_ps if p ∉ constrained)
    cs = Set{Symbol}(metabolites(m))
    inv = Set(K for K in _binding_K_symbols(M) if K ∉ constrained)
    # Try algebraic factoring on the numerator (only if it reduces terms)
    rxns = reactions(m)
    eq_steps = equilibrium_steps(m)
    enz_set = Set(e[1] for e in enzyme_forms(m))
    pc = param_constraints(m)
    binding_Ks = Set(K for K in _binding_K_symbols(M) if K ∉ constrained)
    num_fs = _try_algebraic_factor_sigma(
        num, rxns, eq_steps, enz_set, pc; binding_Ks,
    )
    use_factored_num = if num_fs !== nothing
        # Count display terms: trivial products expand, non-trivial count as 1
        n = 0
        for (c, p) in zip(num_fs.coefficients, num_fs.products)
            trivial = length(p.factors) == 1 && p.factors[1] == poly_one()
            n += trivial ? length(c) : 1
        end
        n < length(num)
    else
        false
    end
    num_expr = if use_factored_num
        _factored_sigma_to_expr(num_fs::FactoredSigma, ps, cs, inv)
    else
        _poly_to_expr(num, ps, cs, inv)
    end
    num_str = _expr_to_string(num_expr)
    den_str = _expr_to_string(_denom_terms_to_expr(denom_terms, ps, cs, inv))
    "v = E_total * ($num_str) / ($den_str)"
end

function rate_equation_string(::M, ::FullMode) where {M<:EnzymeMechanism}
    lines = ["(; $(join(_sorted_raw_param_symbols(M), ", "))) = params",
             "(; $(join(_sorted_conc_symbols(M), ", "))) = concs"]
    pc = param_constraints(M())
    for (target, coeff, factors) in pc
        push!(lines, _user_constraint_to_string(target, coeff, factors))
    end
    push!(lines, _format_rate_equation_core(M))
    join(lines, "\n")
end

function rate_equation_string(::M, ::ReducedMode) where {M<:EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    join(["(; $(join((indep..., :Keq, :E_total), ", "))) = params",
          "(; $(join(_sorted_conc_symbols(M), ", "))) = concs",
          _constraint_expr_strings(M)...,
          _format_rate_equation_core(M)], "\n")
end

# ─── Structural Identifiability ───────────────────────────────

function _count_rate_monomials(M::Type{<:EnzymeMechanism})
    num, denom_terms = _raw_symbolic_rate_polys(M)
    den = _expand_to_poly(denom_terms)
    mm(mono) = sort([
        s => e for (s, e) in mono
        if !is_k_parameter(s) && s != :E_total && s != :Keq
    ])
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

