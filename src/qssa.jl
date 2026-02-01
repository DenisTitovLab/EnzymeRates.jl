using Combinatorics: permutations

function _steps_from_species_reactions(species_data, reactions_data)
    enzs = species_data[4]
    enz_names = Tuple(e[1] for e in enzs)
    steps = map(enumerate(reactions_data)) do (step_idx, (lhs, rhs))
        e_lhs = nothing
        e_rhs = nothing
        m_lhs = nothing
        m_rhs = nothing
        for s in lhs
            if s in enz_names
                e_lhs = s
            else
                m_lhs = s
            end
        end
        for s in rhs
            if s in enz_names
                e_rhs = s
            else
                m_rhs = s
            end
        end
        i = findfirst(==(e_lhs), enz_names)
        j = findfirst(==(e_rhs), enz_names)
        kf = Symbol("k$(step_idx)f")
        kr = Symbol("k$(step_idx)r")
        (i, j, kf, kr, m_lhs, m_rhs)
    end
    return Tuple(steps)
end

"""
    _symbolic_rate_expr(Species, Reactions, PNames, CNames)

Build a symbolic `Expr` for the QSSA rate equation using Laplacian cofactor
determinants expanded via the Leibniz formula. The rate is defined as the
net consumption of the first substrate, normalized by its stoichiometry.
"""
function _symbolic_rate_expr(species_data, reactions_data, PNames, CNames)
    subs, prods, _, enzs = species_data
    isempty(subs) && error("No substrates defined")
    ref_name = subs[1][1]
    nu_ref = 0
    for (name, _) in subs
        name == ref_name && (nu_ref -= 1)
    end
    for (name, _) in prods
        name == ref_name && (nu_ref += 1)
    end
    nu_ref == 0 && error("Reference substrate has zero net stoichiometry")

    Steps = _steps_from_species_reactions(species_data, reactions_data)
    N = length(enzs)

    # 1. Build symbolic rate matrix R[i,j] as Expr (or 0 meaning absent)
    R = Matrix{Any}(fill(0, N, N))

    for (i, j, kf, kr, met_f, met_r) in Steps
        # Forward: i → j
        fwd = met_f === nothing ? :(params.$kf) : :(params.$kf * concs.$met_f)
        R[i, j] = R[i, j] == 0 ? fwd : :($(R[i, j]) + $fwd)

        # Reverse: j → i
        rev = met_r === nothing ? :(params.$kr) : :(params.$kr * concs.$met_r)
        R[j, i] = R[j, i] == 0 ? rev : :($(R[j, i]) + $rev)
    end

    # 2. Build symbolic Laplacian
    L = Matrix{Any}(fill(0, N, N))
    for i in 1:N
        diag_terms = Any[]
        for j in 1:N
            i == j && continue
            if R[i, j] != 0
                L[i, j] = :(-$(R[i, j]))
                push!(diag_terms, R[i, j])
            end
        end
        if !isempty(diag_terms)
            L[i, i] = length(diag_terms) == 1 ? diag_terms[1] : Expr(:call, :+, diag_terms...)
        end
    end

    # 3. Cofactor determinants via Leibniz formula with sparsity pruning
    D = Vector{Any}(undef, N)
    for root in 1:N
        # Delete row root and col root
        rows = [r for r in 1:N if r != root]
        cols = [c for c in 1:N if c != root]
        n_sub = N - 1

        if n_sub == 0
            D[root] = 1
            continue
        end

        terms = Any[]
        for perm in permutations(1:n_sub)
            # Compute sign of permutation
            sign = _perm_sign(perm)

            # Product of L_sub[k, perm[k]] = L[rows[k], cols[perm[k]]]
            factors = Any[]
            all_nonzero = true
            for k in 1:n_sub
                entry = L[rows[k], cols[perm[k]]]
                if entry == 0
                    all_nonzero = false
                    break
                end
                push!(factors, entry)
            end
            !all_nonzero && continue

            prod_expr = length(factors) == 1 ? factors[1] : Expr(:call, :*, factors...)
            if sign == -1
                prod_expr = :(-$prod_expr)
            end
            push!(terms, prod_expr)
        end

        if isempty(terms)
            D[root] = 0
        elseif length(terms) == 1
            D[root] = terms[1]
        else
            D[root] = Expr(:call, :+, terms...)
        end
    end

    # 4. Denominator: sum of all cofactors
    nonzero_D = [D[k] for k in 1:N if D[k] != 0]
    denom = length(nonzero_D) == 1 ? nonzero_D[1] : Expr(:call, :+, nonzero_D...)

    # E_total
    has_etotal = :E_total in PNames
    et_expr = has_etotal ? :(params.E_total) : 1.0

    # 5. Net consumption of reference substrate
    terms = Any[]
    for (i, j, kf, kr, met_f, met_r) in Steps
        rf_expr = met_f === nothing ? :(params.$kf) : :(params.$kf * concs.$met_f)
        rr_expr = met_r === nothing ? :(params.$kr) : :(params.$kr * concs.$met_r)
        flux = :($et_expr * ($rf_expr * $(D[i]) - $rr_expr * $(D[j])) / $denom)
        if met_f === ref_name
            push!(terms, flux)
        elseif met_r === ref_name
            push!(terms, :(-$flux))
        end
    end

    net_expr = isempty(terms) ? 0 :
               length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...)
    abs_nu = abs(nu_ref)
    abs_nu == 1 ? net_expr : :($net_expr / $abs_nu)
end

"""Compute sign of a permutation (as +1 or -1)."""
function _perm_sign(perm)
    n = length(perm)
    visited = falses(n)
    sign = 1
    for i in 1:n
        visited[i] && continue
        visited[i] = true
        cycle_len = 1
        j = perm[i]
        while j != i
            visited[j] = true
            cycle_len += 1
            j = perm[j]
        end
        if iseven(cycle_len)
            sign = -sign
        end
    end
    sign
end

"""
    rate_equation(::EnzymeMechanism{Species,Reactions}, params::NamedTuple, concs::NamedTuple)

Compute the QSSA steady-state rate (net consumption of the first substrate,
normalized by its stoichiometric coefficient). The body is generated at
compile time as a single arithmetic expression with no allocations, loops,
or matrix ops.
"""
@generated function rate_equation(
    ::EnzymeMechanism{SpeciesT, Reactions},
    params::NamedTuple{PNames},
    concs::NamedTuple{CNames}
) where {SpeciesT, Reactions, PNames, CNames}
    _symbolic_rate_expr(SpeciesT, Reactions, PNames, CNames)
end

"""
    rate_equation_string(m::EnzymeMechanism)

Return a string representation of the rate equation, matching what `rate_equation` computes.
"""
function rate_equation_string(m::EnzymeMechanism)
    _rate_equation_string(m)
end

function _rate_equation_string(::EnzymeMechanism{SpeciesT, Reactions}) where {SpeciesT, Reactions}
    Steps = _steps_from_species_reactions(SpeciesT, Reactions)

    met_names = Symbol[]
    param_names = Symbol[]
    for (i, j, kf, kr, met_f, met_r) in Steps
        met_f !== nothing && met_f ∉ met_names && push!(met_names, met_f)
        met_r !== nothing && met_r ∉ met_names && push!(met_names, met_r)
        push!(param_names, kf)
        push!(param_names, kr)
    end

    PNames = tuple(param_names..., :E_total)
    CNames = tuple(met_names...)
    expr = _symbolic_rate_expr(SpeciesT, Reactions, PNames, CNames)

    # Convert to string, then strip `params.` and `concs.` prefixes
    s = string(expr)
    s = replace(s, "params." => "")
    s = replace(s, "concs." => "")
    s
end
