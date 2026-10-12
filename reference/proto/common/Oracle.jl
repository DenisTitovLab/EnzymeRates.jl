# ABOUTME: Numeric steady-state oracle (module Oracle) for plain Pop.PMech mechanisms, in BigFloat.
# ABOUTME: Solves the full mass-action network with RE steps as fast pairs; stdlib only.
isdefined(Main, :Pop) || include(joinpath(@__DIR__, "Population.jl"))

module Oracle

using LinearAlgebra, Random
import ..Pop

const PREC = 512
const LAMBDA_EXP = 40

big512(x) = setprecision(() -> BigFloat(x), BigFloat, PREC)

"""
The generator of `m`'s mass-action network at constants `k` and concentrations `concs`,
with the normalization row Σx = 1 in place of its last row, and each step's
`(ω, from, to, forward rate, reverse rate)`, ω being the copies of `substrates[1]` the
step takes up minus those it gives off (counted by name). Conventions as in
`steady_state_v`. Call inside `setprecision(BigFloat, PREC)`.
"""
function _network(m::Pop.PMech, k::AbstractVector, concs::AbstractDict, Λ)
    Λb = Λ === nothing ? big(10)^LAMBDA_EXP : BigFloat(Λ)
    n = length(m.species)
    M = zeros(BigFloat, n, n)
    c(x) = BigFloat(concs[x.name])
    x1 = m.substrates[1]
    steps = Tuple{Int, Int, Int, BigFloat, BigFloat}[]
    for (j, s) in enumerate(m.steps)
        if s.re
            K = BigFloat(k[j][1])
            f, r = Pop.is_binding(s) ? (Λb, Λb * K) : (Λb, Λb / K)
        else
            f, r = BigFloat(k[j][1]), BigFloat(k[j][2])
        end
        rf = f * prod((c(x) for x in s.consumed); init = big(1))
        rr = r * prod((c(x) for x in s.released); init = big(1))
        a, b = s.from, s.to
        M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
        ω = count(x -> x.name === x1, s.consumed) - count(x -> x.name === x1, s.released)
        push!(steps, (ω, a, b, rf, rr))
    end
    M[n, :] .= 1
    M, steps
end

"""Steady-state form populations of `m` (summing to 1), under the conventions of
`steady_state_v`."""
function steady_state_populations(m::Pop.PMech, k::AbstractVector, concs::AbstractDict;
                                  Λ = nothing)
    setprecision(BigFloat, PREC) do
        M, _ = _network(m, k, concs, Λ)
        rhs = zeros(BigFloat, size(M, 1)); rhs[end] = 1
        M \ rhs
    end
end

"""
    steady_state_v(m, k, concs; Λ = 10^40) -> BigFloat

Net consumption rate of the first substrate (`m.substrates[1]`) per unit enzyme at steady
state. `k[j]` holds step `j`'s constants in the package's convention: `[kf, kr]` for a
steady-state step, in its stored direction (`kf` times the consumed concentrations drives
from → to, `kr` times the released ones drives to → from); `[K]` for a rapid-equilibrium
step, realized as the fast pair `Λ` / `Λ·K` for a binding (one metabolite consumed, none
released; `K` a dissociation constant, `[to]/[from] = [M]/K`) and `Λ` / `Λ/K` for any other
step (`K = [to]·Π[released] / ([from]·Π[consumed])`). `concs` maps metabolite names to
concentrations; an inhibitor copy reads its reactant's name. Works in 512-bit BigFloat.
"""
function steady_state_v(m::Pop.PMech, k::AbstractVector, concs::AbstractDict;
                        Λ = nothing)
    setprecision(BigFloat, PREC) do
        M, steps = _network(m, k, concs, Λ)
        rhs = zeros(BigFloat, size(M, 1)); rhs[end] = 1
        x = M \ rhs
        sum((ω * (rf * x[a] - rr * x[b]) for (ω, a, b, rf, rr) in steps if ω != 0);
            init = big(0))
    end
end

"""Log-uniform BigFloat in [10^lo, 10^hi]."""
logu(rng, lo = -2, hi = 2) = setprecision(() -> big(10)^(lo + (hi - lo) * big(rand(rng))),
                                          BigFloat, PREC)

"""
    consistent_constants(m, rng; Keq) -> (k, G, mu)

Bare per-step constants (each step its own, no kinetic groups) that satisfy every Haldane
and Wegscheider relation: random free energies `G` for the forms (uniform in [-3, 3]),
standard potentials `mu` for the metabolites (by name; substrates and products shifted so
`Σmu[S] − Σmu[P] = log Keq`), the association constant of step j
`Ka = exp(G[from] + Σmu[consumed] − G[to] − Σmu[released])`, and a random log-uniform
prefactor. A steady-state step gets `[kf, kf/Ka]`; a rapid-equilibrium binding `[1/Ka]`
(its dissociation constant), any other rapid-equilibrium step `[Ka]`.
"""
function consistent_constants(m::Pop.PMech, rng::AbstractRNG; Keq)
    setprecision(BigFloat, PREC) do
        G = [6 * big(rand(rng)) - 3 for _ in m.species]
        names = unique([m.substrates; m.products; m.inhibitors; m.regulators;
                        [x.name for s in m.steps for x in [s.consumed; s.released]]])
        mu = Dict{Symbol, BigFloat}(x => 6 * big(rand(rng)) - 3 for x in names)
        gap = log(BigFloat(Keq)) - (sum(mu[x] for x in m.substrates) -
                                    sum(mu[x] for x in m.products))
        for x in m.products
            mu[x] -= gap / length(m.products)
        end
        k = Vector{BigFloat}[]
        for s in m.steps
            Ka = exp(G[s.from] + sum((mu[x.name] for x in s.consumed); init = big(0)) -
                     G[s.to] - sum((mu[x.name] for x in s.released); init = big(0)))
            if s.re
                push!(k, [Pop.is_binding(s) ? 1 / Ka : Ka])
            else
                kf = logu(rng)
                push!(k, [kf, kf / Ka])
            end
        end
        k, G, mu
    end
end

"""Random concentrations (log-uniform in [1e-2, 1e2]) for every metabolite of `m`, with the
first product rescaled so the mass-action ratio `Π[P]/Π[S]` equals `gamma`."""
function concs_at_ratio(m::Pop.PMech, rng::AbstractRNG, gamma)
    c = Dict{Symbol, BigFloat}(x => logu(rng) for x in m.metabolites)
    setprecision(BigFloat, PREC) do
        Γ = prod(c[x] for x in m.products) / prod(c[x] for x in m.substrates)
        c[m.products[1]] *= BigFloat(gamma) / Γ
    end
    c
end

"""Per-step constants from named values: step j takes `values[name]` for each of its
names `m.step_names[j]`, a name folded by `rename` (`absorbed => target`) taking its
target's value."""
function constants_from_names(m::Pop.PMech, values::AbstractDict;
                              rename = Dict{Symbol, Symbol}(), names = m.step_names)
    [BigFloat[big512(values[get(rename, n, n)]) for n in ns] for ns in names]
end

"""Random constants that respect the kinetic groups: one log-uniform value per distinct
name in `m.step_names`, shared by every step that renders it."""
function group_constants(m::Pop.PMech, rng::AbstractRNG)
    vals = Dict{Symbol, BigFloat}()
    for ns in m.step_names, n in ns
        haskey(vals, n) || (vals[n] = logu(rng))
    end
    constants_from_names(m, vals), vals
end

"""Value of a `Pop.PPoly` at `values` (symbols missing from `values` raise)."""
function eval_poly(p::Pop.PPoly, values::AbstractDict)
    setprecision(BigFloat, PREC) do
        sum((BigFloat(numerator(c)) / denominator(c) *
             prod((big512(values[s])^e for (s, e) in mono); init = big(1))
             for (mono, c) in p); init = big(0))
    end
end

end # module
