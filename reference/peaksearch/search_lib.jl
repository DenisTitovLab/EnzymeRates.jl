# ABOUTME: Search machinery for rate peaks: the joint objective R = v(c) / v_sat over (rate constants,
# ABOUTME: concentrations), an IPOP-CMA-ES optimizer, a feasible pattern-search polish, inner sups.

using Random, LinearAlgebra, Printf

isdefined(@__MODULE__, :PeakSolver) || include(joinpath(@__DIR__, "solver.jl"))
using .PeakSolver
isdefined(@__MODULE__, :MECHS) || include(joinpath(@__DIR__, "mechs.jl"))

const NCONC = Dict(:zero => 2, :Q => 3, :P => 3, :PQ => 4)

function concs(setting, lc, ::Type{T} = Float64) where {T}
    A = exp10(T(lc[1])); B = exp10(T(lc[2]))
    setting === :zero && return (A, B, zero(T), zero(T))
    setting === :Q && return (A, B, zero(T), exp10(T(lc[3])))
    setting === :P && return (A, B, exp10(T(lc[3])), zero(T))
    (A, B, exp10(T(lc[3])), exp10(T(lc[4])))
end

"""
R = v(c) / max over saturation directions of v(A, B -> inf; P = Q = 0), the violation of
the rate-constant box and of the concentration box [clo, chi] (log10), and the pieces.
"""
function evalR(m, x, setting, lK; clo = -8.0, chi = 8.0, T = Float64)
    nθ = ntheta(m)
    kf, kr, viol = rates(m, view(x, 1:nθ), lK; T = T)
    lc = view(x, nθ+1:length(x))
    for l in lc
        viol += max(0.0, l - chi) + max(0.0, clo - l)
    end
    v = rate(m, kf, kr, concs(setting, lc, T))
    vs = vsat_dirs(m, kf, kr)
    vsb = maximum(vs)
    (R = v / vsb, viol = viol, v = v, vsat = vsb, vsdirs = vs)
end

mutable struct Tracker
    R::Float64
    x::Vector{Float64}
    nevals::Int
end

function make_objective(m, setting, lK, tr::Tracker; kw...)
    x -> begin
        tr.nevals += 1
        maximum(abs, x) > 60 && return 1e12 + sum(abs, x)
        r = evalR(m, x, setting, lK; kw...)
        isfinite(r.R) || return 1e11
        if r.viol == 0 && r.R > tr.R
            tr.R = r.R
            tr.x = collect(x)
        end
        score = r.R > 0 ? log10(r.R) : -300 + r.R / (1 + abs(r.R))
        -score + 100 * r.viol
    end
end

"""Standard (μ/μ_w, λ)-CMA-ES minimizing F from x0 with step size σ0."""
function cmaes(F, x0, σ0, rng; λ = 4 + floor(Int, 3log(length(x0))), maxevals = 10^5,
               tolx = 1e-10)
    n = length(x0)
    μ = λ ÷ 2
    w = [log(μ + 0.5) - log(i) for i in 1:μ]
    w ./= sum(w)
    μeff = 1 / sum(abs2, w)
    cc = (4 + μeff / n) / (n + 4 + 2μeff / n)
    cs = (μeff + 2) / (n + μeff + 5)
    c1 = 2 / ((n + 1.3)^2 + μeff)
    cμ = min(1 - c1, 2(μeff - 2 + 1 / μeff) / ((n + 2)^2 + μeff))
    damps = 1 + 2max(0, sqrt((μeff - 1) / (n + 1)) - 1) + cs
    chiN = sqrt(n) * (1 - 1 / (4n) + 1 / (21n^2))
    xm = collect(Float64, x0)
    σ = σ0
    pc = zeros(n); ps = zeros(n)
    B = Matrix{Float64}(I, n, n); D = ones(n)
    C = Matrix{Float64}(I, n, n); invsqrtC = Matrix{Float64}(I, n, n)
    X = zeros(n, λ); Y = zeros(n, λ); f = zeros(λ)
    eigeneval = 0; counteval = 0; gen = 0
    hist = Float64[]
    window = 10 + ceil(Int, 30n / λ)
    while counteval < maxevals
        gen += 1
        for k in 1:λ
            Y[:, k] = B * (D .* randn(rng, n))
            X[:, k] = xm .+ σ .* Y[:, k]
            f[k] = F(view(X, :, k))
        end
        counteval += λ
        idx = sortperm(f)
        push!(hist, f[idx[1]])
        yw = zeros(n)
        for i in 1:μ
            yw .+= w[i] .* view(Y, :, idx[i])
        end
        xm .+= σ .* yw
        ps .= (1 - cs) .* ps .+ sqrt(cs * (2 - cs) * μeff) .* (invsqrtC * yw)
        hsig = norm(ps) / sqrt(1 - (1 - cs)^(2gen)) / chiN < 1.4 + 2 / (n + 1)
        pc .= (1 - cc) .* pc .+ (hsig ? sqrt(cc * (2 - cc) * μeff) : 0.0) .* yw
        Cμ = zeros(n, n)
        for i in 1:μ
            y = view(Y, :, idx[i])
            Cμ .+= w[i] .* (y * y')
        end
        C .= (1 - c1 - cμ) .* C .+ c1 .* (pc * pc' .+ (hsig ? 0.0 : cc * (2 - cc)) .* C) .+
             cμ .* Cμ
        σ *= exp((cs / damps) * (norm(ps) / chiN - 1))
        if counteval - eigeneval > λ / (c1 + cμ) / n / 10
            eigeneval = counteval
            C .= (C .+ C') ./ 2
            E = eigen(Symmetric(C))
            D .= sqrt.(max.(E.values, 1e-30))
            B .= E.vectors
            invsqrtC .= B * Diagonal(1 ./ D) * B'
        end
        σ * maximum(D) < tolx && break
        σ > 1e3 && break
        maximum(D) / minimum(D) > 1e7 && break
        gen > window && minimum(@view hist[1:gen-window]) - minimum(hist) < 1e-12 && break
    end
    counteval
end

"""Coordinate pattern search on x maximizing R among feasible points only."""
function polish(m, setting, lK, x; steps = (1.0, 0.3, 0.1, 0.03, 0.01, 1e-3, 1e-4, 1e-5),
                maxsweeps = 200, kw...)
    Rof(y) = (r = evalR(m, y, setting, lK; kw...); r.viol == 0 && isfinite(r.R) ? r.R : -Inf)
    x = copy(x)
    fx = Rof(x)
    for h in steps
        for _ in 1:maxsweeps
            improved = false
            for i in eachindex(x), s in (-h, h)
                y = copy(x)
                y[i] += s
                fy = Rof(y)
                if fy > fx
                    x, fx = y, fy
                    improved = true
                end
            end
            improved || break
        end
    end
    x, fx
end

"""
IPOP-CMA-ES restarts from random starting points; returns the best feasible point found
(after polish) and the total evaluation count.
"""
function search(m, setting, lK; seed = 1, restarts = 7, budget = 1_500_000, σ0 = 2.0,
                clo = -8.0, chi = 8.0)
    rng = MersenneTwister(seed)
    nθ = ntheta(m)
    n = nθ + NCONC[setting]
    tr = Tracker(-Inf, Float64[], 0)
    F = make_objective(m, setting, lK, tr; clo = clo, chi = chi)
    λ0 = 4 + floor(Int, 3log(n))
    used = 0
    for r in 1:restarts
        used >= budget && break
        x0 = vcat(3 .* (2 .* rand(rng, nθ) .- 1),
                  clo .+ (chi - clo) .* rand(rng, NCONC[setting]))
        used += cmaes(F, x0, σ0, rng; λ = λ0 * 2^(r - 1), maxevals = budget - used)
    end
    isempty(tr.x) && return (R = -Inf, x = Float64[], evals = tr.nevals)
    x, R = polish(m, setting, lK, tr.x; clo = clo, chi = chi)
    (R = R, x = x, evals = tr.nevals)
end

"""Sup over concentrations of v / v_sat for fixed rate constants: grid, then pattern search."""
function inner_sup(m, θ, lK, setting; clo = -8.0, chi = 8.0)
    nc = NCONC[setting]
    h = nc == 2 ? 0.25 : nc == 3 ? 0.5 : 1.0
    grid = collect(clo:h:chi)
    kf, kr, _ = rates(m, θ, lK)
    vsb = maximum(vsat_dirs(m, kf, kr))
    best, bestc = -Inf, zeros(nc)
    for lc in Iterators.product(ntuple(_ -> grid, nc)...)
        R = rate(m, kf, kr, concs(setting, lc)) / vsb
        R > best && ((best, bestc) = (R, collect(lc)))
    end
    for hh in (h / 2, h / 8, h / 32, 1e-3, 1e-4)
        for _ in 1:200
            improved = false
            for i in 1:nc, s in (-hh, hh)
                y = copy(bestc)
                y[i] = clamp(y[i] + s, clo, chi)
                R = rate(m, kf, kr, concs(setting, y)) / vsb
                R > best && ((best, bestc, improved) = (R, y, true))
            end
            improved || break
        end
    end
    best, bestc
end

"""
θ from rate constants: `lkf[i]` (log10 kf) for every steady-state step i, `lkr[i]` for the
steps that fix the free energies (a spanning tree grown from E); the other reverse rate
constants follow from thermodynamics.
"""
function theta_from(m, lK, lkf::Dict{Int, Float64}, lkr::Dict{Int, Float64})
    n = nss(m)
    θ = zeros(ntheta(m))
    known = falses(length(m.forms)); known[1] = true
    g = zeros(length(m.forms))
    gl(l) = l == 4 ? -lK : 0.0
    progress = true
    while progress
        progress = false
        for (i, s) in enumerate(m.steps)
            haskey(lkr, i) || continue
            lKa = lkf[i] - lkr[i]
            sc = sum((gl(l) for l in s.cons); init = 0.0)
            sr = sum((gl(l) for l in s.rel); init = 0.0)
            if known[s.from] && !known[s.to]
                g[s.to] = g[s.from] + sc - sr - lKa
                known[s.to] = true; progress = true
            elseif known[s.to] && !known[s.from]
                g[s.from] = g[s.to] + sr - sc + lKa
                known[s.from] = true; progress = true
            end
        end
    end
    all(known) || error("theta_from: the given reverse constants do not span the forms")
    j = 0
    for (i, s) in enumerate(m.steps)
        s.re && continue
        j += 1
        θ[j] = lkf[i]
    end
    θ[n+1:end] .= g[2:end]
    θ
end

function describe(m, θ, lK)
    kf, kr, viol = rates(m, θ, lK)
    io = IOBuffer()
    for (i, s) in enumerate(m.steps)
        lhs = join([String(m.forms[s.from]); [String(LIGANDS[l]) for l in s.cons]], " + ")
        rhs = join([String(m.forms[s.to]); [String(LIGANDS[l]) for l in s.rel]], " + ")
        if s.re
            @printf(io, "      %-22s RE  K = %.3e\n", "$lhs <-> $rhs", kf[i] / kr[i])
        else
            @printf(io, "      %-22s kf = %.3e  kr = %.3e\n", "$lhs <-> $rhs", kf[i], kr[i])
        end
    end
    @printf(io, "      box violation = %.2e\n", viol)
    String(take!(io))
end
