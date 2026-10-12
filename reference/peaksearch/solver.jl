# ABOUTME: Independent mass-action steady-state solver for enzyme mechanisms (GTH state reduction)
# ABOUTME: with a free-energy parameterization that makes every cycle thermodynamically consistent.

module PeakSolver

using LinearAlgebra

export Mech, mech, ntheta, nss, rates, steady_state, turnover, netprod, rate, vsat_dirs,
       log10Ka, LIGANDS

const LIGANDS = (:A, :B, :P, :Q)
lig(s::Symbol) = findfirst(==(s), LIGANDS)

struct Step
    from::Int
    to::Int
    cons::Vector{Int}   # ligands taken up in the forward direction
    rel::Vector{Int}    # ligands given off in the forward direction
    re::Bool            # rapid equilibrium
end

struct Mech
    name::String
    forms::Vector{Symbol}   # forms[1] is free enzyme E (free energy 0)
    steps::Vector{Step}
    cut::Vector{Int}        # steps whose net forward fluxes sum to the turnover rate
end

"""
    mech(name, forms, stepdefs, cut)

`stepdefs` are tuples `(from, consumed, to, released, re)`, e.g.
`(:EA, [:B], :EAB, Symbol[], false)` for `E(A) + B <--> E(A, B)`.
"""
function mech(name, forms, stepdefs, cut)
    fi = Dict(f => i for (i, f) in enumerate(forms))
    steps = [Step(fi[a], fi[b], Int[lig(x) for x in c], Int[lig(x) for x in r], re)
             for (a, c, b, r, re) in stepdefs]
    Mech(name, forms, steps, cut)
end

nss(m::Mech) = count(s -> !s.re, m.steps)
ntheta(m::Mech) = nss(m) + length(m.forms) - 1

"""
    log10Ka(m, θ, log10Keq, i, T)

log10 of step i's equilibrium constant (forward/reverse, concentration units as written).
Free energies (log10 units): forms from θ (E = 0), ligands A = B = P = 0, Q = -log10 Keq, so
every cycle's product of equilibrium constants is 1 (futile) or Keq (catalytic).
"""
function log10Ka(m::Mech, θ, log10Keq, i, ::Type{T}) where {T}
    s = m.steps[i]
    n = nss(m)
    g(f) = f == 1 ? zero(T) : T(θ[n + f - 1])
    gl(l) = l == 4 ? -T(log10Keq) : zero(T)
    g(s.from) + sum((gl(l) for l in s.cons); init = zero(T)) -
        g(s.to) - sum((gl(l) for l in s.rel); init = zero(T))
end

"""
    rates(m, θ, log10Keq; T = Float64, Λ = 1e60)

Forward and reverse rate constants of every step, and the total violation (log10 units) of
the box 1e-6 ≤ k ≤ 1e6 on every steady-state rate constant (1e-12 ≤ K ≤ 1e12 for rapid-
equilibrium steps, which run at rate scale Λ). θ = [log10 kf of each steady-state step in
step order..., log10 free energy of forms[2:end]...].
"""
function rates(m::Mech, θ, log10Keq; T = Float64, Λ = T(1e60))
    ns = length(m.steps)
    kf = zeros(T, ns)
    kr = zeros(T, ns)
    viol = zero(T)
    j = 0
    for i in 1:ns
        lKa = log10Ka(m, θ, log10Keq, i, T)
        if m.steps[i].re
            kf[i] = Λ * exp10(lKa / 2)
            kr[i] = Λ * exp10(-lKa / 2)
            viol += max(zero(T), abs(lKa) - 12)
        else
            j += 1
            lkf = T(θ[j])
            lkr = lkf - lKa
            viol += max(zero(T), abs(lkf) - 6) + max(zero(T), abs(lkr) - 6)
            kf[i] = exp10(lkf)
            kr[i] = exp10(lkr)
        end
    end
    kf, kr, viol
end

"""
Stationary distribution of a CTMC from its off-diagonal rate matrix Q (overwritten) by the
Grassmann–Taksar–Heyman state reduction, which never subtracts and so keeps componentwise
relative accuracy however stiff the rates are.
"""
function gth!(Q::Matrix{T}) where {T}
    n = size(Q, 1)
    for k in n:-1:2
        s = zero(T)
        @inbounds for j in 1:k-1
            s += Q[k, j]
        end
        s > 0 || error("GTH: state $k has no exit to lower-index states")
        @inbounds for i in 1:k-1
            Q[i, k] /= s
        end
        @inbounds for j in 1:k-1
            qkj = Q[k, j]
            iszero(qkj) && continue
            for i in 1:k-1
                i == j && continue
                Q[i, j] += Q[i, k] * qkj
            end
        end
    end
    π = zeros(T, n)
    π[1] = one(T)
    @inbounds for j in 2:n
        acc = zero(T)
        for i in 1:j-1
            acc += π[i] * Q[i, j]
        end
        π[j] = acc
    end
    π ./= sum(π)
end

fwd(m, kf, c, i) = (r = kf[i]; for l in m.steps[i].cons; r *= c[l]; end; r)
rev(m, kr, c, i) = (r = kr[i]; for l in m.steps[i].rel; r *= c[l]; end; r)

function steady_state(m::Mech, kf::Vector{T}, kr::Vector{T}, c) where {T}
    n = length(m.forms)
    Q = zeros(T, n, n)
    for (i, s) in enumerate(m.steps)
        Q[s.from, s.to] += fwd(m, kf, c, i)
        Q[s.to, s.from] += rev(m, kr, c, i)
    end
    gth!(Q)
end

stepnet(m, kf, kr, c, π, i) =
    fwd(m, kf, c, i) * π[m.steps[i].from] - rev(m, kr, c, i) * π[m.steps[i].to]

"""Net production rate of ligand l (per unit total enzyme)."""
function netprod(m::Mech, kf, kr, c, π, l)
    v = zero(eltype(π))
    for (i, s) in enumerate(m.steps)
        ν = count(==(l), s.rel) - count(==(l), s.cons)
        ν == 0 && continue
        v += ν * stepnet(m, kf, kr, c, π, i)
    end
    v
end

"""
Turnover rate (A + B -> P + Q per unit enzyme). Each substrate's net uptake, each product's net
release and the net flux through the cut steps all equal it; the candidate with the smallest
gross (forward plus reverse) flux loses the least to cancellation. With P = 0 or Q = 0 that
candidate is the product's release, a sum without subtraction.
"""
function turnover(m::Mech, kf, kr, c, π)
    T = eltype(π)
    best_v, best_g = zero(T), T(Inf)
    for l in 1:4
        v, g = zero(T), zero(T)
        for (i, s) in enumerate(m.steps)
            ν = count(==(l), s.rel) - count(==(l), s.cons)
            ν == 0 && continue
            f = fwd(m, kf, c, i) * π[s.from]
            r = rev(m, kr, c, i) * π[s.to]
            v += (l <= 2 ? -ν : ν) * (f - r)
            g += abs(ν) * (f + r)
        end
        g < best_g && ((best_v, best_g) = (v, g))
    end
    v, g = zero(T), zero(T)
    for i in m.cut
        f = fwd(m, kf, c, i) * π[m.steps[i].from]
        r = rev(m, kr, c, i) * π[m.steps[i].to]
        v += f - r
        g += f + r
    end
    g < best_g ? v : best_v
end

rate(m::Mech, kf, kr, c) = turnover(m, kf, kr, c, steady_state(m, kf, kr, c))

"""
Rates at A, B -> infinity with no products along three directions: (A, B) =
(10^a, 10^b) for (a, b) in `dirs`.
"""
function vsat_dirs(m::Mech, kf::Vector{T}, kr::Vector{T};
                   dirs = ((48, 48), (48, 24), (24, 48))) where {T}
    [rate(m, kf, kr, (exp10(T(a)), exp10(T(b)), zero(T), zero(T))) for (a, b) in dirs]
end

end
