using Combinatorics: permutations

"""
    _symbolic_rate_expr(::Type{<:EnzymeMechanism})

Build a symbolic `Expr` for the QSSA rate equation using Laplacian cofactor
determinants expanded via the Leibniz formula. The rate is defined as the
net consumption of the first substrate, normalized by its stoichiometry.
"""
function _symbolic_rate_expr(M::Type{<:EnzymeMechanism})
    m = M()
    subs = substrates(m)
    prods = products(m)
    enzs = enzyme_forms(m)
    rxns = reactions(m)

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

    enz_names = Tuple(e[1] for e in enzs)
    N = length(enzs)

    # 1. Build symbolic rate matrix R[i,j] as Expr (or 0 meaning absent)
    R = Matrix{Any}(fill(0, N, N))

    for (idx, (lhs, rhs)) in enumerate(rxns)
        e_lhs = first(s for s in lhs if s in enz_names)
        e_rhs = first(s for s in rhs if s in enz_names)
        i = findfirst(==(e_lhs), enz_names)
        j = findfirst(==(e_rhs), enz_names)
        kf = Symbol("k$(idx)f")
        kr = Symbol("k$(idx)r")
        m_it = iterate(s for s in lhs if s ∉ enz_names)
        met_f = m_it === nothing ? nothing : m_it[1]
        m_it = iterate(s for s in rhs if s ∉ enz_names)
        met_r = m_it === nothing ? nothing : m_it[1]

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
    et_expr = :(params.E_total)

    # 5. Net consumption of reference substrate
    terms = Any[]
    for (idx, (lhs, rhs)) in enumerate(rxns)
        e_lhs = first(s for s in lhs if s in enz_names)
        e_rhs = first(s for s in rhs if s in enz_names)
        i = findfirst(==(e_lhs), enz_names)
        j = findfirst(==(e_rhs), enz_names)
        kf = Symbol("k$(idx)f")
        kr = Symbol("k$(idx)r")
        m_it = iterate(s for s in lhs if s ∉ enz_names)
        met_f = m_it === nothing ? nothing : m_it[1]
        m_it = iterate(s for s in rhs if s ∉ enz_names)
        met_r = m_it === nothing ? nothing : m_it[1]

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
    m::M, params::NamedTuple, concs::NamedTuple
) where {M <: EnzymeMechanism}
    _symbolic_rate_expr(M)
end

# Polynomial representation: monomial (sorted Symbol vector) → integer coefficient.
# Converting the Expr tree to this form automatically expands products over sums,
# cancels opposite terms (coefficients sum to zero), and sorts factors within each monomial.
const _Poly = Dict{Vector{Symbol},Int}
# Sort key: k-constants before metabolites, alphabetically within each group
_sk(s) = (startswith(string(s), "k") ? 0 : 1, string(s))
# Multiply two polynomials: distribute every term pair, concatenate and sort monomials
_pmul(a::_Poly, b::_Poly) = (r = _Poly(); for (k1,v1) in a, (k2,v2) in b; k = sort([k1;k2]; by=_sk); r[k] = get(r,k,0) + v1*v2; end; filter!(p -> p[2] != 0, r))
# Add two polynomials: merge coefficients, drop zeros (this is where cancellation happens)
_padd(a::_Poly, b::_Poly) = (r = copy(a); for (k,v) in b; r[k] = get(r,k,0) + v; end; filter!(p -> p[2] != 0, r))

# Recursively convert an Expr (without /) into a polynomial, stripping params./concs. prefixes
function _to_poly(e)
    e isa Expr && e.head == :. && return _Poly([e.args[2].value] => 1)
    (e isa Expr && e.head == :call) || return _Poly(Symbol[] => (e isa Integer ? e : 1))
    op, a = e.args[1], e.args[2:end]
    op == :* ? reduce(_pmul, _to_poly.(a)) : op == :+ ? reduce(_padd, _to_poly.(a)) :
    op == :- && length(a) == 1 ? _Poly(k => -v for (k,v) in _to_poly(a[1])) :
    op == :- ? _padd(_to_poly(a[1]), _Poly(k => -v for (k,v) in _to_poly(a[2]))) : _Poly(Symbol[] => 1)
end

# Split numerator and denominator: _strip_div removes all / nodes (keeping numerator sides),
# _find_denom finds the first denominator in the tree (all fractions share the same one)
_strip_div(e) = (e isa Expr && e.head == :call) ? (e.args[1] == :/ ? _strip_div(e.args[2]) :
    Expr(:call, e.args[1], _strip_div.(e.args[2:end])...)) : e
_find_denom(e) = (e isa Expr && e.head == :call) ? (e.args[1] == :/ ? e.args[3] :
    foldl((r, a) -> r !== nothing ? r : _find_denom(a), e.args[2:end]; init=nothing)) : nothing

# Pretty-print a polynomial: positive terms first, then negative, joined with + / -
function _poly_str(p::_Poly)
    isempty(p) && return "0"
    ts = sort(collect(p); by=x -> (x[2] < 0, x[1]))
    join([begin m = isempty(k) ? "$(abs(v))" : abs(v) == 1 ? join(k, " * ") : "$(abs(v)) * " * join(k, " * ")
        i == 1 ? (v < 0 ? "-$m" : m) : (v < 0 ? " - $m" : " + $m") end for (i,(k,v)) in enumerate(ts)])
end

"""
    rate_equation_string(m::EnzymeMechanism)

Return a string representation of the rate equation, matching what `rate_equation` computes.
"""
function rate_equation_string(::M) where {M<:EnzymeMechanism}
    expr = _symbolic_rate_expr(M)
    np = _Poly(filter(!=(:E_total), k) => v for (k,v) in _to_poly(_strip_div(expr)))
    "E_total * ($(_poly_str(np))) / ($(_poly_str(_to_poly(_find_denom(expr)))))"
end
