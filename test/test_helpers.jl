using LinearAlgebra
using Random

"""
Independent reference: compute King-Altman rate using Laplacian cofactor method.
"""
function reference_king_altman(m::EnzymeMechanism, params::NamedTuple, concs::NamedTuple; E_total=1.0)
    forms = enzyme_forms(m)
    n = length(forms)
    name_to_idx = Dict(s.name => i for (i, s) in enumerate(forms))

    # Build rate matrix R[i,j] = pseudo-first-order rate from i to j
    R = zeros(n, n)
    for (step_idx, (lhs, rhs)) in enumerate(m.steps)
        e_lhs = [s for s in lhs if s.role == enzyme][1]
        e_rhs = [s for s in rhs if s.role == enzyme][1]
        i = name_to_idx[e_lhs.name]
        j = name_to_idx[e_rhs.name]

        m_lhs = [s for s in lhs if s.role == metabolite]
        m_rhs = [s for s in rhs if s.role == metabolite]

        kf = params[Symbol("k$(step_idx)f")]
        kr = params[Symbol("k$(step_idx)r")]

        rf = kf
        if !isempty(m_lhs)
            rf *= concs[m_lhs[1].name]
        end
        R[i, j] += rf

        rr = kr
        if !isempty(m_rhs)
            rr *= concs[m_rhs[1].name]
        end
        R[j, i] += rr
    end

    # Build Laplacian
    L = zeros(n, n)
    for i in 1:n
        for j in 1:n
            if i != j
                L[i, j] = -R[i, j]
                L[i, i] += R[i, j]
            end
        end
    end

    # Cofactors
    D = zeros(n)
    for i in 1:n
        rows = [r for r in 1:n if r != i]
        cols = [c for c in 1:n if c != i]
        D[i] = det(L[rows, cols])
    end

    D_total = sum(D)
    E_conc = D ./ D_total .* E_total

    # Compute rate using first step
    lhs, rhs = m.steps[1]
    e_lhs = [s for s in lhs if s.role == enzyme][1]
    e_rhs = [s for s in rhs if s.role == enzyme][1]
    i = name_to_idx[e_lhs.name]
    j = name_to_idx[e_rhs.name]

    m_lhs = [s for s in lhs if s.role == metabolite]
    m_rhs = [s for s in rhs if s.role == metabolite]

    kf = params[:k1f]
    kr = params[:k1r]

    rf = kf * (isempty(m_lhs) ? 1.0 : concs[m_lhs[1].name])
    rr = kr * (isempty(m_rhs) ? 1.0 : concs[m_rhs[1].name])

    v = rf * E_conc[i] - rr * E_conc[j]
    return v
end

"""
    _unicyclic_denominator(f, r)

King–Altman denominator for a unicyclic network with forward rates `f[i]`
and reverse rates `r[i]` around the cycle (1-indexed; cyclic).
"""
function _unicyclic_denominator(f::AbstractVector, r::AbstractVector)
    n = length(f)
    @assert length(r) == n
    T = promote_type(eltype(f), eltype(r))
    D = zero(T)

    for i in 1:n
        Ti = zero(T)
        for m in 0:n-1
            pr = one(T)
            for j in 0:m-1
                pr *= r[mod1(i + j, n)]
            end
            pf = one(T)
            for j in (m+1):(n-1)
                pf *= f[mod1(i + j, n)]
            end
            Ti += pr * pf
        end
        D += Ti
    end

    return D
end

"""
    _unicyclic_flux(f, r, Et) -> v

Closed-form QSSA flux for a unicyclic enzyme-state cycle.
v = Et*(prod(f)-prod(r)) / sum(tree_weights)
"""
function _unicyclic_flux(f::AbstractVector, r::AbstractVector, Et)
    pf = prod(f)
    pr = prod(r)
    D  = _unicyclic_denominator(f, r)
    return Et * (pf - pr) / D
end


function random_params_concs(m::EnzymeMechanism, met_names::Vector{Symbol}; rng=Random.default_rng())
    n_steps = length(m.steps)
    param_keys = Symbol[]
    param_vals = Float64[]
    for i in 1:n_steps
        push!(param_keys, Symbol("k$(i)f"))
        push!(param_vals, 0.1 + 9.9 * rand(rng))
        push!(param_keys, Symbol("k$(i)r"))
        push!(param_vals, 0.1 + 9.9 * rand(rng))
    end
    params = NamedTuple{Tuple(param_keys)}(Tuple(param_vals))

    conc_vals = [0.1 + 9.9 * rand(rng) for _ in met_names]
    concs = NamedTuple{Tuple(met_names)}(Tuple(conc_vals))

    return params, concs
end
