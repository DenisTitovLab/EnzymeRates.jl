# ABOUTME: Estimates duplicates a rename-proof key cannot see: pairs of laws with one concentration
# ABOUTME: support but different keys where one law reproduces the other at fitted values (LM fit).
include(joinpath(@__DIR__, "Laws.jl"))
include(joinpath(@__DIR__, "Factor.jl"))
using Random, LinearAlgebra, Printf, Statistics, Serialization

const POPS = ARGS[2:end]
const VARIANT = Symbol(ARGS[1])
const ms = reduce(vcat, [Pop.load_population(p) for p in POPS])
const rng = Xoshiro(11)

"""A law compiled for Float64 evaluation over the variables `vars` (fitted names first)."""
struct Fast
    comps::Vector{Vector{Tuple{Float64, Vector{Tuple{Int, Int}}}}}
    kind::Symbol
    n::Int
    mults::Vector{Int}
    np::Int
    vars::Vector{Symbol}
end

function compile(b::Laws.Built, law::Canon.Law, m)
    vars = [law.params; :Keq; :L; sort(unique([m.metabolites; m.regulators]))]
    at = Dict(v => i for (i, v) in enumerate(vars))
    comps = [[(Float64(c), [(at[s], e) for (s, e) in mo]) for (mo, c) in p] for p in Canon.normalize(law)]
    Fast(comps, b.kind, b.n, b.mults, length(law.params), vars)
end

function value(f::Fast, x)
    e = [sum((c * prod((x[i]^k for (i, k) in mo); init = 1.0) for (c, mo) in p); init = 0.0)
         for p in f.comps]
    f.kind === :plain && return e[1] / e[2]
    nA, qA, nI, qI = e[1] / e[2], e[3] / e[2], e[4] / e[5], e[6] / e[5]
    WA = prod((e[5 + 2k]^f.mults[k] for k in eachindex(f.mults)); init = 1.0)
    WI = prod((e[6 + 2k]^f.mults[k] for k in eachindex(f.mults)); init = 1.0)
    n, L = f.n, x[f.np + 2]
    (nA * qA^(n - 1) * WA + L * nI * qI^(n - 1) * WI) / (qA^n * WA + L * qI^n * WI)
end

"""Concentration support of each component (the sorted concentration monomials present)."""
support(law::Canon.Law, m) = [sort!(unique([[x for x in mo if x.first in Set([m.metabolites; m.regulators])]
                                            for mo in keys(p)])) for p in Canon.normalize(law)]

"""Whether law `fb` reproduces law `fa` at random fitted values of `fa`: Levenberg–Marquardt in log
space over `fb`'s fitted names from several starts, on 4·P concentration points. Returns the best
max relative residual and the largest |log value| at that fit."""
function mimics(fa::Fast, fb::Fast; starts = 8, iters = 150)
    P = fb.np
    xa = [exp(4 * rand(rng) - 2) for _ in 1:fa.np]
    keq, L = exp(4 * rand(rng) - 2), exp(4 * rand(rng) - 2)
    nc = length(fa.vars) - fa.np - 2
    pts = [[exp(4 * rand(rng) - 2) for _ in 1:nc] for _ in 1:4P]
    y = [value(fa, [xa; keq; L; c]) for c in pts]
    res(u) = [(value(fb, [exp.(u); keq; L; c]) - yi) / abs(yi) for (c, yi) in zip(pts, y)]
    best, bestu = Inf, Float64[]
    for _ in 1:starts
        u = 4 .* rand(rng, P) .- 2
        r = res(u); μ = 1e-3
        for _ in 1:iters
            J = hcat([(res(u .+ 1e-7 .* (1:P .== k)) .- r) ./ 1e-7 for k in 1:P]...)
            all(isfinite, J) && all(isfinite, r) || break
            A = J' * J; g = J' * r
            δ = -((A + μ * Diagonal(diag(A) .+ 1e-12)) \ g)
            r2 = res(u .+ δ)
            if all(isfinite, r2) && sum(abs2, r2) < sum(abs2, r)
                u, r, μ = u .+ δ, r2, μ / 3
            else
                μ *= 4
            end
            maximum(abs, r) < 1e-12 && break
        end
        mx = all(isfinite, r) ? maximum(abs, r) : Inf
        if mx < best
            best, bestu = mx, copy(u)
        end
    end
    # the fit must also hold at fresh concentrations it never saw
    hold = [[exp(6 * rand(rng) - 3) for _ in 1:nc] for _ in 1:20]
    out = maximum(abs((value(fb, [exp.(bestu); keq; L; c]) - value(fa, [xa; keq; L; c])) /
                      value(fa, [xa; keq; L; c])) for c in hold)
    max(best, out), maximum(abs, bestu; init = 0.0)
end

d = deserialize(joinpath(@__DIR__, "data", "$(join(POPS, "+")).jls"))
key = Dict(r[:id] => r[(VARIANT, :k2c, :key)] for r in d.recs if haskey(r, (VARIANT, :k2c, :key)))
fast, sup = Dict{Int, Fast}(), Dict{Int, Any}()
for (i, m) in enumerate(ms)
    haskey(key, m.id) || continue
    b = VARIANT === :today ? (m.allosteric ? Laws.mwc_today(m) : Laws.plain_today(m)) :
        (m.allosteric ? Laws.mwc_redesign(m) : Laws.plain_redesign(m))
    law = first(Factor.cancel(b.law))
    fast[i] = compile(b, law, m)
    sup[i] = (b.kind, length(law.params), support(law, m))
end
buckets = Dict{Any, Vector{Int}}()
for (i, s) in sup
    push!(get!(buckets, s, Int[]), i)
end
pairs = Tuple{Int, Int}[]
for v in values(buckets)
    reps = Dict{String, Int}()
    for i in sort(v)
        get!(reps, key[ms[i].id], i)
    end
    r = sort(collect(values(reps)))
    for x in eachindex(r), y in x+1:length(r)
        push!(pairs, (r[x], r[y]))
    end
end
println(join(POPS, "+"), " / ", VARIANT, ": ", length(fast), " laws, ", length(buckets),
        " support buckets, ", length(pairs), " cross-key pairs in a shared bucket")
sample = shuffle(rng, pairs)[1:min(end, parse(Int, get(ENV, "NPAIRS", "60")))]
hits = 0
t = @elapsed for (a, b) in sample
    r1, u1 = mimics(fast[a], fast[b])
    r1 < 1e-8 || continue
    r2, u2 = mimics(fast[a], fast[b])          # a second random point of a
    r3, u3 = mimics(fast[b], fast[a])          # and the other direction
    same = r2 < 1e-8 && r3 < 1e-8
    global hits += same
    @printf("   %s ~ %s: residuals %.1e %.1e %.1e, max |log x| %.1f %.1f %.1f%s\n", ms[a].id, ms[b].id,
            r1, r2, r3, u1, u2, u3, same ? "  SAME MODEL (both directions)" : "")
end
@printf("sampled %d pairs in %.0f s: %d reproduce each other both ways\n", length(sample), t, hits)
