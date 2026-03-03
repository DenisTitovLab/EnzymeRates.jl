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
            e_rhs, _ = _split_reaction_side(rhs, enz_set)
            i_f = findfirst(==(e_lhs), enz_names)
            j_f = findfirst(==(e_rhs), enz_names)
            (i_f ∈ sg_set && j_f ∈ sg_set) || continue
            # Canonical form: RE binding steps have metabolite on LHS.
            # Skip RE isomerizations (no metabolite on either side).
            isempty(m_lhs) && continue

            met = first(m_lhs)
            # LHS = [E, M], RHS = [EM]: less-bound on LHS, more-bound on RHS
            less, more = i_f, j_f

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
        # Canonical form: RE binding steps have metabolite on LHS.
        # Skip RE isomerizations (no metabolite).
        isempty(m_lhs) && continue
        e_rhs, _ = _split_reaction_side(rhs, enz_set)
        i_f = findfirst(==(e_lhs), enz_names)
        j_f = findfirst(==(e_rhs), enz_names)
        (i_f ∈ sg_set && j_f ∈ sg_set) || continue

        met = first(m_lhs)
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
Try to detect Q^n structure directly from per-subunit alpha values.
For a homodimer with competitive metabolites on each of n identical
subunits, forms are multisets and `_partition_metabolite_sites` fails.
This function constructs Q from the per-subunit forms (total binding
≤ max_bind_per_subunit) and verifies Q^n == full sigma.

`sg_sigma` is the pre-computed constrained sub-group sigma polynomial.
`alpha_vals` maps form index → POLY alpha value (ratio or numerator).

When some forms are "dead-end" (only bind metabolites not present in
any multi-copy form), they are separated as additive remainder terms.
Returns `(FactoredPoly, remainder::POLY)` or `nothing`.
"""
function _try_subunit_power(
    sg, binding_states, alpha_vals, sg_sigma, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
)
    # Find metabolites that appear with count > 1 (multi-subunit metabolites)
    max_count = Dict{Symbol, Int}()
    for f in sg
        for (met, cnt) in binding_states[f]
            max_count[met] = max(get(max_count, met, 0), cnt)
        end
    end
    isempty(max_count) && return nothing

    # Identify "core" metabolites (count ≥ 2) vs "dead-end only" metabolites
    core_mets = Set{Symbol}(m for (m, c) in max_count if c >= 2)
    isempty(core_mets) && return nothing
    n = maximum(max_count[m] for m in core_mets)
    n < 2 && return nothing

    # Partition forms: core (only core metabolites) vs dead-end (has non-core)
    core_forms = Int[]
    deadend_forms = Int[]
    for f in sg
        bs = binding_states[f]
        has_noncore = any(c > 0 for (m, c) in bs if m ∉ core_mets)
        if has_noncore
            push!(deadend_forms, f)
        else
            push!(core_forms, f)
        end
    end
    length(core_forms) < 3 && return nothing

    # Try each candidate max_bind_per_subunit on core forms
    for max_bind in 1:(n - 1)
        n % max_bind != 0 && continue
        sub_n = n ÷ max_bind

        # Collect per-subunit forms, deduplicated by binding state.
        # Symmetric positions (e.g., E_0S and E_S0) have identical
        # binding states and constrained alpha — keep one representative.
        seen_states = Dict{Dict{Symbol,Int}, Int}()
        subunit_forms = Int[]
        for f in core_forms
            total = sum(values(binding_states[f]); init=0)
            total > max_bind && continue
            bs = Dict(m => c for (m, c) in binding_states[f] if c > 0)
            haskey(seen_states, bs) && continue
            seen_states[bs] = f
            push!(subunit_forms, f)
        end
        length(subunit_forms) < 2 && continue

        # Normalize alpha values by dividing out the reference form's alpha
        # (conformational coefficient like K37 for T-state sub-groups)
        ref_alpha = alpha_vals[first(subunit_forms)]
        if ref_alpha != poly_one()
            normed = f -> _poly_div_mono(alpha_vals[f], ref_alpha)
        else
            normed = f -> alpha_vals[f]
        end

        Q = reduce(poly_add, (normed(f) for f in subunit_forms))
        Q = _apply_param_constraints(Q, constraints; binding_Ks)

        Qn = Q
        for _ in 2:sub_n
            Qn = poly_mul(Qn, Q)
        end

        # Compute core sigma (normalized) and verify Q^n matches
        core_sigma = reduce(
            poly_add, (normed(f) for f in core_forms),
        )
        core_sigma = _apply_param_constraints(
            core_sigma, constraints; binding_Ks,
        )
        Qn != core_sigma && continue

        # Compute dead-end remainder (normalized)
        if isempty(deadend_forms)
            return FactoredPoly([Q], [sub_n]), poly_zero()
        end
        remainder = reduce(
            poly_add, (normed(f) for f in deadend_forms),
        )
        remainder = _apply_param_constraints(
            remainder, constraints; binding_Ks,
        )

        # Check if remainder = Qn * factor (multiplicative dead-end)
        factor = _try_poly_exact_div(remainder, Qn)
        if factor !== nothing
            cofactor = poly_add(poly_one(), factor)
            return FactoredPoly([Q, cofactor], [sub_n, 1]), poly_zero()
        end

        # Additive remainder: sigma = Qn + remainder
        return FactoredPoly([Q], [sub_n]), remainder
    end
    nothing
end

"""
Try to decompose polynomial Q into multiplicative per-site factors
using `_try_algebraic_factor_sigma` iteratively.
Returns a `FactoredPoly` with per-site factors, or the original
single-factor `FactoredPoly` if decomposition fails.
"""
function _decompose_product(
    fp::FactoredPoly, rxns, eq_steps, enz_set, constraints;
    binding_Ks::Set{Symbol}=Set{Symbol}(),
)
    length(fp.factors) != 1 && return fp
    Q = fp.factors[1]
    n = fp.exponents[1]
    length(Q) < 3 && return fp

    factors = POLY[]
    remaining = Q
    while true
        afs = _try_algebraic_factor_sigma(
            remaining, rxns, eq_steps, enz_set, constraints;
            binding_Ks,
        )
        afs === nothing && break
        # Check if it's a pure multiplicative factoring (one term, coeff=1)
        length(afs.coefficients) != 1 && break
        afs.coefficients[1] != poly_one() && break
        prod = afs.products[1]
        length(prod.factors) < 2 && break
        # Extract the smaller factors, keep the largest as remaining
        for f in prod.factors[1:end-1]
            push!(factors, f)
        end
        remaining = prod.factors[end]
    end
    push!(factors, remaining)
    length(factors) == 1 && return fp
    FactoredPoly(factors, fill(n, length(factors)))
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

        # Compute constrained sub-group sigma (used by step 5 + fallback)
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

        # Try structural factoring (steps 2-6)
        fp = nothing
        states = _metabolite_binding_states(
            sg, enz_names, enz_set, rxns, eq_steps,
        )
        if length(states) == length(sg)
            sites = _partition_metabolite_sites(states)
            if isempty(sites) && length(sg) == 1
                # Single form, no metabolites → trivial factor
                push!(factored_products, FactoredPoly([poly_one()], [1]))
                push!(final_coeffs, sg_coeff)
                continue
            end
            is_product, _ = _check_cartesian_product(sg, states, sites)
            if is_product
                sigma_K_syms = Set{Symbol}(
                    s for (mono, _) in sg_sigma for (s, _) in mono
                    if startswith(string(s), "K") &&
                       length(string(s)) > 1 &&
                       isdigit(string(s)[2])
                )
                consistent, met_to_K = _check_k_consistency(
                    sg, enz_names, enz_set, rxns, eq_steps,
                    states, sites, constraints; sigma_K_syms,
                )
                if consistent
                    fp = _build_site_factors(
                        sg, enz_names, enz_set, rxns, eq_steps,
                        alpha_num, alpha_den, sites, met_to_K,
                        states; alpha_ratio,
                    )
                    ref_alpha = use_ratio ?
                        alpha_ratio[first(sg)] :
                        alpha_num[first(sg)]
                    if ref_alpha != poly_one()
                        fp = FactoredPoly(
                            [_poly_div_mono(f, ref_alpha)
                             for f in fp.factors],
                            copy(fp.exponents),
                        )
                    end
                end
            end
        end

        # Reject trivial structural factoring (single factor, exponent 1)
        # when multi-subunit power might apply (max binding ≥ 2)
        if fp !== nothing && length(fp.factors) == 1 && fp.exponents[1] == 1
            max_bind = 0
            for f in sg
                total = sum(values(states[f]); init=0)
                max_bind = max(max_bind, total)
            end
            max_bind >= 2 && (fp = nothing)
        end

        # Try direct subunit power construction (Q^n from per-subunit alphas)
        if fp === nothing
            alpha_vals = use_ratio ? alpha_ratio : alpha_num
            sp_result = _try_subunit_power(
                sg, states, alpha_vals, sg_sigma, constraints;
                binding_Ks,
            )
            if sp_result !== nothing
                sp_fp, sp_remainder = sp_result
                sp_fp = _decompose_product(
                    sp_fp, rxns, eq_steps, enz_set, constraints;
                    binding_Ks,
                )
                if sp_remainder == poly_zero()
                    fp = sp_fp
                else
                    # Dead-end remainder: wrap as additive FactoredSigma
                    # and handle at this level by constructing the full
                    # factored sigma directly
                    push!(factored_products, sp_fp)
                    push!(final_coeffs, sg_coeff)
                    # Add remainder as a separate unfactored term
                    push!(
                        factored_products,
                        FactoredPoly([poly_one()], [1]),
                    )
                    push!(final_coeffs, poly_mul(sg_coeff, sp_remainder))
                    continue
                end
            end
        end

        # Fallback: try poly_power on constrained subgroup sigma
        if fp === nothing
            sg_sigma_div = sg_sigma
            coeff_poly = _apply_param_constraints(
                coeffs[sg_idx], constraints; binding_Ks,
            )
            if coeff_poly != poly_one()
                sg_sigma_div = _poly_div_mono(sg_sigma_div, coeff_poly)
            end
            pfs = _try_poly_power(sg_sigma_div)
            if pfs !== nothing
                fp = pfs.products[1]
                fp = _decompose_product(
                    fp, rxns, eq_steps, enz_set, constraints; binding_Ks,
                )
            end
        end

        fp === nothing && return nothing
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
    const_mono = MONO()
    haskey(p, const_mono) && p[const_mono] == 1 || return nothing

    # max exponent of any symbol determines max possible power
    max_e = 0
    for (mono, _) in p
        for (_, e) in mono
            max_e = max(max_e, e)
        end
    end
    max_e < 2 && return nothing

    # Try each candidate power n from max_e down to 2
    for n in max_e:-1:2
        # Find terms whose exponents are all divisible by n with coeff 1
        root_terms = MONO[]
        for (mono, coeff) in p
            isempty(mono) && continue
            coeff == 1 || continue
            all(e % n == 0 for (_, e) in mono) || continue
            push!(root_terms, mono)
        end
        isempty(root_terms) && continue

        # Form candidate Q = 1 + sum(nth-root terms)
        Q = POLY(MONO() => Rational{Int}(1))
        for mono in root_terms
            Q[MONO([s => div(e, n) for (s, e) in mono])] = 1
        end

        # Compute Q^n and verify
        Qn = Q
        for _ in 2:n
            Qn = poly_mul(Qn, Q)
        end
        if Qn == p
            return FactoredSigma(
                [poly_one()], [FactoredPoly([Q], [n])],
            )
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

"""Sum polynomial entries from `polys` at indices `group`."""
function _sum_group_polys(polys::Vector{POLY}, group)
    reduce(poly_add, (polys[i] for i in group))
end

"""
Build substitution pairs merging Haldane-derived parameters that have
identical expressions after user constraints. E.g., if k8r and k7r both
resolve to `k7f / (K1 * Keq)`, returns `[:k8r => :k7r]`.
"""
function _haldane_equality_substitutions(M::Type{<:EnzymeMechanism})
    dep_exprs, _ = _dependent_param_exprs(M)
    length(dep_exprs) < 2 && return Pair{Symbol,Symbol}[]
    # Group by expression identity
    groups = Dict{Any,Vector{Symbol}}()
    for (sym, expr) in dep_exprs
        pushed = false
        for (key, vec) in groups
            if key == expr
                push!(vec, sym)
                pushed = true
                break
            end
        end
        pushed || (groups[expr] = [sym])
    end
    # Build substitution pairs: later symbols → first in group
    # Sort by step number (e.g., k7r=7, k10r=10) for stable canonical choice
    _step_num(s) = something(tryparse(Int, match(r"\d+", string(s)).match), 0)
    subs = Pair{Symbol,Symbol}[]
    for (_, syms) in groups
        length(syms) < 2 && continue
        sort!(syms; by=_step_num)
        canonical = first(syms)
        for s in @view syms[2:end]
            push!(subs, s => canonical)
        end
    end
    subs
end

function _is_K_sym(s::Symbol)
    str = string(s)
    length(str) > 1 && str[1] == 'K' && isdigit(str[2])
end

"""
Try factoring numerator by dividing it with sigma factors from the
factored denominator. For each DenomTerm sigma factor σ_i, attempts to
find quotient Q_i such that `num = Σ Q_i * σ_i`. Works by partitioning
numerator terms by their binding-K symbols (each conformational subgroup
has distinct K params).
"""
function _try_factor_num_by_denom_sigmas(num::POLY, denom_terms)
    # Collect unique sigma base factors from denominator
    sigma_factors = POLY[]
    for dt in denom_terms
        for fp in dt.sigma.products
            for f in fp.factors
                f != poly_one() && f ∉ sigma_factors &&
                    push!(sigma_factors, f)
            end
        end
    end
    isempty(sigma_factors) && return nothing

    # Partition numerator into sub-polynomials by K symbols
    # (single-factor case: all terms go to one group)
    K_sets = [
        Set{Symbol}(s for (mono, _) in σ for (s, _) in mono if _is_K_sym(s))
        for σ in sigma_factors
    ]
    sub_polys = [POLY() for _ in sigma_factors]
    for (mono, coeff) in num
        mono_Ks = Set{Symbol}(s for (s, _) in mono if _is_K_sym(s))
        assigned = false
        for (i, Ks) in enumerate(K_sets)
            if !isempty(intersect(mono_Ks, Ks))
                sub_polys[i][mono] = coeff
                assigned = true
                break
            end
        end
        assigned || return nothing
    end

    coefficients = POLY[]
    products = FactoredPoly[]
    for (i, σ) in enumerate(sigma_factors)
        isempty(sub_polys[i]) && continue
        Q = _try_poly_exact_div(sub_polys[i], σ)
        Q === nothing && return nothing
        # Extract common factor from Q: 2*k7f*S/K1 → 2 * (k7f*S/K1)
        g, Q_reduced = _extract_poly_common_factor(Q)
        push!(coefficients, g)
        push!(products, FactoredPoly([Q_reduced, σ], [1, 1]))
    end

    # Verify round-trip
    expanded = reduce(
        poly_add,
        (poly_mul(c, _expand_factored_poly(p))
         for (c, p) in zip(coefficients, products)),
    )
    expanded != num && return nothing
    FactoredSigma(coefficients, products)
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
                fp = _decompose_product(
                    pfs.products[1], rxns, eq_steps, enz_set, pc;
                    binding_Ks,
                )
                pfs = FactoredSigma([poly_one()], [fp])
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

    # Wrap numerator in FactoredSigma (try factoring, fall back to trivial)
    num_fs = _try_poly_power(num)
    if num_fs !== nothing
        fp = _decompose_product(
            num_fs.products[1], rxns, eq_steps, enz_set, pc;
            binding_Ks,
        )
        num_fs = FactoredSigma([poly_one()], [fp])
    end
    if num_fs === nothing
        afs = _try_algebraic_factor_sigma(
            num, rxns, eq_steps, enz_set, pc; binding_Ks,
        )
        if afs !== nothing
            # Only use if factoring reduces display terms
            n = 0
            for (c, p) in zip(afs.coefficients, afs.products)
                trivial = length(p.factors) == 1 &&
                    p.factors[1] == poly_one()
                n += trivial ? length(c) : 1
            end
            if n < length(num)
                num_fs = afs
            end
        end
    end
    # Try denom-guided factoring: divide numerator by denominator sigma factors
    dfs = _try_factor_num_by_denom_sigmas(num, denom_terms)
    if dfs !== nothing
        n_dfs = sum(
            (length(fp.factors) == 1 && fp.factors[1] == poly_one()) ?
                length(c) : 1
            for (c, fp) in zip(dfs.coefficients, dfs.products)
        )
        n_cur = num_fs === nothing ? length(num) : sum(
            (length(fp.factors) == 1 && fp.factors[1] == poly_one()) ?
                length(c) : 1
            for (c, fp) in zip(num_fs.coefficients, num_fs.products)
        )
        if n_dfs < n_cur
            num_fs = dfs
        end
    end
    if num_fs === nothing
        num_fs = FactoredSigma(
            [poly_one()], [FactoredPoly([num], [1])],
        )
    end

    num_fs, denom_terms
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
    num_fs, denom_terms = _raw_symbolic_rate_polys(M)
    m = M()
    constrained = Set(c[1] for c in param_constraints(m))
    raw_ps = _raw_param_symbols(equilibrium_steps(m))
    ps = Set{Symbol}(p for p in raw_ps if p ∉ constrained)
    cs = Set{Symbol}(metabolites(m))
    inv = Set(K for K in _binding_K_symbols(M) if K ∉ constrained)
    # Numerator is already factored by _raw_symbolic_rate_polys
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
    num_fs, denom_terms = _raw_symbolic_rate_polys(M)
    num = _expand_factored_sigma(num_fs)
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
    T_subs = Dict(
        p => _rename_params_T(p)
        for p in Iterators.flatten([keys(dep_R), indep_R])
        if is_k_parameter(p) && p != :Keq
    )
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

# ─── Parameters API ──────────────────────────────────────────────

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
    T_subs = Dict(K => _rename_params_T(K) for K in cat_params if is_k_parameter(K))
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

    # This call only accepts Type{<:EnzymeMechanism}, which narrows CM's type for JET
    # so the subsequent _apply_kd_inversion call type-checks. Results go to cat_params.
    num_fs_body, _ = _raw_symbolic_rate_polys(CM)
    m_cat = CM()
    cat_constr = Set(c[1] for c in param_constraints(m_cat))
    cat_params = Set{Symbol}(
        p for p in _raw_param_symbols(equilibrium_steps(m_cat)) if p ∉ cat_constr
    )

    # Dep param assignments (R-state): k3r = inv(K2)*K1*k3f*inv(Keq), etc.
    dep_R, _ = _dependent_param_exprs(CM)
    dep_R_kd = _apply_kd_inversion(dep_R, CM, K -> :(inv($K)))
    r_assignments = Expr[
        Expr(:(=), sym, dep_R_kd[sym])
        for (sym, _) in sort(collect(dep_R_kd); by=first)
    ]

    # T-state dep param assignments (NConf=2 only)
    t_assignments = Expr[]
    if NConf == 2
        T_subs = Dict(K => _rename_params_T(K) for K in cat_params if is_k_parameter(K))
        t_assignments = Expr[
            Expr(:(=), _rename_params_T(sym), substitute_params_expr(dep_R_kd[sym], T_subs))
            for (sym, _) in sort(collect(dep_R_kd); by=first)
        ]
    end

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

function rate_equation(m::OligomericEnzymeMechanism, concs, params)
    rate_equation(m, concs, params, Reduced)
end

@generated function rate_equation(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf},
    concs::NamedTuple, params::NamedTuple, ::ReducedMode,
) where {Mets,CM,CatN,RS,NConf}
    _build_oligomeric_rate_body(Mets, CM, CatN, RS, NConf)
end

# ─── String representation ────────────────────────────────────────

rate_equation_string(m::OligomericEnzymeMechanism) = rate_equation_string(m, Reduced)

function rate_equation_string(
    ::OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}, ::ReducedMode,
) where {Mets,CM,CatN,RS,NConf}
    M = OligomericEnzymeMechanism{Mets,CM,CatN,RS,NConf}
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)

    # Dep param assignment strings (R-state)
    dep_R, indep_R = _dependent_param_exprs(CM)
    dep_R_kd = _apply_kd_inversion(dep_R, CM, K -> :(1 / $K))
    dep_lines = [
        "$sym = $(_expr_to_string(dep_R_kd[sym]))"
        for (sym, _) in sort(collect(dep_R_kd); by=first)
    ]

    # T-state dep param assignments (NConf=2)
    t_dep_lines = if NConf == 2
        T_subs = Dict(
            p => _rename_params_T(p)
            for p in Iterators.flatten([keys(dep_R), indep_R])
            if is_k_parameter(p) && p != :Keq
        )
        [
            "$(_rename_params_T(sym)) = $(_expr_to_string(substitute_params_expr(dep_R_kd[sym], T_subs)))"
            for (sym, _) in sort(collect(dep_R_kd); by=first)
        ]
    else
        String[]
    end

    # Build the v = line using _oligomeric_num_den_exprs directly.
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

