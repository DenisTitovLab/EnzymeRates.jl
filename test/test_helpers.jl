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

# Safe getter: missing concentration -> 0
_getc(c::NamedTuple, sym::Symbol) = haskey(c, sym) ? getfield(c, sym) : 0.0

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

"""
    rate_seq_uniuni(p, c) -> v
Ordered sequential uni–uni:
E + S <-> ES
ES <-> E + P
"""
function rate_seq_uniuni(p::NamedTuple, c::NamedTuple)
    S = _getc(c, :S); P = _getc(c, :P)
    f = [
        p.k1f * S,   # bind S
        p.k2f        # release P forward
    ]
    r = [
        p.k1r,       # unbind S
        p.k2r * P    # bind P (reverse of release)
    ]
    return _unicyclic_flux(f, r, p.Et)
end

"""
    rate_seq_biuni(p, c) -> v

Ordered sequential bi–uni:
E + S1 <-> ES1
ES1 + S2 <-> ES1S2
ES1S2 <-> EP1          (chem step)
EP1 <-> E + P1
"""
function rate_seq_biuni(p::NamedTuple, c::NamedTuple)
    S1 = _getc(c, :S1); S2 = _getc(c, :S2); P1 = _getc(c, :P1)

    f = [
        p.k1f * S1,   # bind S1
        p.k2f * S2,   # bind S2
        p.k3f,        # chem forward
        p.k4f         # release P1 forward
    ]
    r = [
        p.k1r,        # unbind S1
        p.k2r,        # unbind S2
        p.k3r,        # chem reverse
        p.k4r * P1    # bind P1 (reverse of release)
    ]

    return _unicyclic_flux(f, r, p.Et)
end

"""
    rate_seq_bibi(p, c) -> v

Ordered sequential bi–bi:
E + S1 <-> ES1
ES1 + S2 <-> ES1S2
ES1S2 <-> EP1P2            (chem step)
EP1P2 <-> EP2 + P1         (release P1)
EP2 <-> E + P2             (release P2)
"""
function rate_seq_bibi(p::NamedTuple, c::NamedTuple)
    S1 = _getc(c, :S1); S2 = _getc(c, :S2)
    P1 = _getc(c, :P1); P2 = _getc(c, :P2)

    f = [
        p.k1f * S1,
        p.k2f * S2,
        p.k3f,
        p.k4f,
        p.k5f
    ]
    r = [
        p.k1r,
        p.k2r,
        p.k3r,
        p.k4r * P1,
        p.k5r * P2
    ]

    return _unicyclic_flux(f, r, p.Et)
end

"""
    rate_seq_terbi(p, c) -> v

Ordered sequential ter–bi:
E + S1 <-> ES1
ES1 + S2 <-> ES1S2
ES1S2 + S3 <-> ES1S2S3
ES1S2S3 <-> EP1P2            (chem step)
EP1P2 <-> EP2 + P1           (release P1)
EP2 <-> E + P2               (release P2)
"""
function rate_seq_terbi(p::NamedTuple, c::NamedTuple)
    S1 = _getc(c, :S1); S2 = _getc(c, :S2); S3 = _getc(c, :S3)
    P1 = _getc(c, :P1); P2 = _getc(c, :P2)

    f = [
        p.k1f * S1,
        p.k2f * S2,
        p.k3f * S3,
        p.k4f,
        p.k5f,
        p.k6f
    ]
    r = [
        p.k1r,
        p.k2r,
        p.k3r,
        p.k4r,
        p.k5r * P1,
        p.k6r * P2
    ]

    return _unicyclic_flux(f, r, p.Et)
end

"""
    rate_seq_terter(p, c) -> v

Ordered sequential ter–ter:
E + S1 <-> ES1
ES1 + S2 <-> ES1S2
ES1S2 + S3 <-> ES1S2S3
ES1S2S3 <-> EP1P2P3           (chem step)
EP1P2P3 <-> EP2P3 + P1        (release P1)
EP2P3 <-> EP3 + P2            (release P2)
EP3 <-> E + P3                (release P3)
"""
function rate_seq_terter(p::NamedTuple, c::NamedTuple)
    S1 = _getc(c, :S1); S2 = _getc(c, :S2); S3 = _getc(c, :S3)
    P1 = _getc(c, :P1); P2 = _getc(c, :P2); P3 = _getc(c, :P3)

    f = [
        p.k1f * S1,
        p.k2f * S2,
        p.k3f * S3,
        p.k4f,
        p.k5f,
        p.k6f,
        p.k7f
    ]
    r = [
        p.k1r,
        p.k2r,
        p.k3r,
        p.k4r,
        p.k5r * P1,
        p.k6r * P2,
        p.k7r * P3
    ]

    return _unicyclic_flux(f, r, p.Et)
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
