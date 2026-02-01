using LinearAlgebra
using Random

"""
Independent reference: compute QSSA rate using Laplacian cofactor method.
Works directly with EnzymeMechanism type parameters.
"""
function reference_qssa(m::EnzymeMechanism{Species, Reactions}, params::NamedTuple, concs::NamedTuple) where {Species, Reactions}
    enzs = enzyme_forms(m)
    n = length(enzs)
    enz_names = Tuple(e[1] for e in enzs)
    name_to_idx = Dict(nm => i for (i, nm) in enumerate(enz_names))
    enz_set = Set(enz_names)

    ref_name, nu_ref = _reference_metabolite(m)

    # Build rate matrix R[i,j] = pseudo-first-order rate from i to j
    R = zeros(n, n)
    for (step_idx, (lhs, rhs)) in enumerate(Reactions)
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        i = name_to_idx[e_lhs]
        j = name_to_idx[e_rhs]

        m_lhs = [s for s in lhs if s ∉ enz_set]
        m_rhs = [s for s in rhs if s ∉ enz_set]

        kf = params[Symbol("k$(step_idx)f")]
        kr = params[Symbol("k$(step_idx)r")]

        rf = isempty(m_lhs) ? kf : kf * concs[m_lhs[1]]
        R[i, j] += rf

        rr = isempty(m_rhs) ? kr : kr * concs[m_rhs[1]]
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
    E_conc = D ./ D_total .* params.E_total

    # Compute net consumption of reference substrate
    v = 0.0
    for (step_idx, (lhs, rhs)) in enumerate(Reactions)
        e_lhs = first(s for s in lhs if s in enz_set)
        e_rhs = first(s for s in rhs if s in enz_set)
        i = name_to_idx[e_lhs]
        j = name_to_idx[e_rhs]

        m_lhs = [s for s in lhs if s ∉ enz_set]
        m_rhs = [s for s in rhs if s ∉ enz_set]

        kf = params[Symbol("k$(step_idx)f")]
        kr = params[Symbol("k$(step_idx)r")]

        rf = kf * (isempty(m_lhs) ? 1.0 : concs[m_lhs[1]])
        rr = kr * (isempty(m_rhs) ? 1.0 : concs[m_rhs[1]])

        flux = rf * E_conc[i] - rr * E_conc[j]
        if !isempty(m_lhs) && m_lhs[1] == ref_name
            v += flux
        elseif !isempty(m_rhs) && m_rhs[1] == ref_name
            v -= flux
        end
    end

    return v / abs(nu_ref)
end

"""
    _unicyclic_denominator(f, r)

QSSA denominator for a unicyclic network with forward rates `f[i]`
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


function random_params_concs(m, met_names::Vector{Symbol}; rng=Random.default_rng())
    ns = n_steps(m)
    param_keys = Symbol[]
    param_vals = Float64[]
    for i in 1:ns
        push!(param_keys, Symbol("k$(i)f"))
        push!(param_vals, 0.1 + 9.9 * rand(rng))
        push!(param_keys, Symbol("k$(i)r"))
        push!(param_vals, 0.1 + 9.9 * rand(rng))
    end
    push!(param_keys, :E_total)
    push!(param_vals, 0.1 + 9.9 * rand(rng))
    params = NamedTuple{Tuple(param_keys)}(Tuple(param_vals))

    conc_vals = [0.1 + 9.9 * rand(rng) for _ in met_names]
    concs = NamedTuple{Tuple(met_names)}(Tuple(conc_vals))

    return params, concs
end

"""
Test that `rate_equation` is non-allocating and fast for the given mechanism.
Must be a standalone function to avoid @testset closure boxing.
"""
function test_rate_equation_performance(m, params, concs)
    rate_equation(m, params, concs) # warmup/compile
    allocs = @allocated rate_equation(m, params, concs)
    t = @elapsed for _ in 1:10_000; rate_equation(m, params, concs); end
    return allocs, t / 10_000
end

function _reference_metabolite(m)
    subs = substrates(m)
    isempty(subs) && error("No substrate found in mechanism")
    name = subs[1][1]
    coeff = -count(s -> s[1] == name, subs)
    return name, coeff
end
