# ABOUTME: Independent verifier library for the rate-peak claims: exact rational steady states,
# ABOUTME: King–Altman spanning-tree rates without subtraction, and cycle-product thermodynamic checks.

module VL

using LinearAlgebra, Printf

const LIG = (:A, :B, :P, :Q)
ligidx(s::Symbol) = (i = findfirst(==(s), LIG); i === nothing ? error("ligand $s") : i)

struct Edge
    u::Int
    w::Int
    cons::Vector{Int}   # ligands taken up going u -> w
    rel::Vector{Int}    # ligands given off going u -> w
end

struct Mech
    name::String
    states::Vector{Symbol}
    edges::Vector{Edge}
end

function mech(name, states, defs)
    id = Dict(s => i for (i, s) in enumerate(states))
    Mech(name, collect(Symbol, states),
         [Edge(id[a], id[b], Int[ligidx(x) for x in c], Int[ligidx(x) for x in r])
          for (a, c, b, r) in defs])
end

const N0 = Symbol[]
const BIND = [(:E, [:A], :EA, N0), (:E, [:B], :EB, N0), (:EA, [:B], :EAB, N0),
              (:EB, [:A], :EAB, N0)]
const REL = [(:EPQ, N0, :EQ, [:P]), (:EPQ, N0, :EP, [:Q]), (:EQ, N0, :E, [:Q]),
             (:EP, N0, :E, [:P])]

# Step and state orders follow the search agent's mechs.jl so its θ vectors decode unchanged.
const MECHS = Dict(
    "T" => mech("textbook", [:E, :EA, :EB, :EAB, :EPQ, :EP, :EQ],
                vcat(BIND, [(:EAB, N0, :EPQ, N0)], REL)),
    "c" => mech("isomerization", [:E, :EA, :EB, :EAB, :EstarAB, :EPQ, :EP, :EQ],
                vcat(BIND, [(:EAB, N0, :EstarAB, N0), (:EstarAB, N0, :EPQ, N0)], REL)),
    "d1" => mech("parallel_sharedEPQ", [:E, :EA, :EB, :EAB, :EstarAB, :EPQ, :EP, :EQ],
                 vcat([(:E, [:A], :EA, N0), (:E, [:B], :EB, N0), (:EA, [:B], :EAB, N0),
                       (:EB, [:A], :EstarAB, N0), (:EAB, N0, :EPQ, N0),
                       (:EstarAB, N0, :EPQ, N0)], REL)),
    "e" => mech("deadend_EAQ", [:E, :EA, :EB, :EAB, :EPQ, :EP, :EQ, :EAQ],
                vcat(BIND, [(:EAB, N0, :EPQ, N0)], REL, [(:EQ, [:A], :EAQ, N0)])),
    "f" => mech("abortive_releasing", [:E, :EA, :EB, :EAB, :EPQ, :EP, :EQ, :EAQ],
                vcat(BIND, [(:EAB, N0, :EPQ, N0)], REL,
                     [(:EQ, [:A], :EAQ, N0), (:EAQ, N0, :EA, [:Q])])),
    "g" => mech("ordered_release", [:E, :EA, :EB, :EAB, :EPQ, :EQ],
                vcat(BIND, [(:EAB, N0, :EPQ, N0), (:EPQ, N0, :EQ, [:P]), (:EQ, N0, :E, [:Q])])),
    "a" => mech("lumped", [:E, :EA, :EB, :EAB], vcat(BIND, [(:EAB, N0, :E, [:P, :Q])])),
    "d2" => mech("parallel_lumped", [:E, :EA, :EB, :EAB, :EstarAB],
                 [(:E, [:A], :EA, N0), (:E, [:B], :EB, N0), (:EA, [:B], :EAB, N0),
                  (:EB, [:A], :EstarAB, N0), (:EAB, N0, :E, [:P, :Q]),
                  (:EstarAB, N0, :E, [:P, :Q])]),
)

nparams(m::Mech) = length(m.edges) + length(m.states) - 1

"""
Rate constants (BigFloat) from the search agent's θ: θ[i] = log10 kf of step i, θ[ns + j] =
log10 free energy of state j + 1 (E = 0); ligand energies A = B = P = 0, Q = -log10 Keq;
log10 K_i = g(from) + Σ g(taken up) - g(to) - Σ g(given off), and kr = kf / K.
"""
function decode(m::Mech, θ, lK)
    ns = length(m.edges)
    length(θ) == nparams(m) || error("θ has $(length(θ)) entries, expected $(nparams(m))")
    g(s) = s == 1 ? big(0.0) : big(θ[ns+s-1])
    gl(l) = l == 4 ? -big(lK) : big(0.0)
    kf = BigFloat[]
    kr = BigFloat[]
    for (i, e) in enumerate(m.edges)
        lKa = g(e.u) + sum((gl(l) for l in e.cons); init = big(0.0)) - g(e.w) -
              sum((gl(l) for l in e.rel); init = big(0.0))
        push!(kf, exp10(big(θ[i])))
        push!(kr, exp10(big(θ[i]) - lKa))
    end
    kf, kr
end

const QQ = Rational{BigInt}
toQ(x) = QQ(x)

"""Pseudo-first-order forward and backward rates of every step at concentrations c = (A, B, P, Q)."""
function pseudo(m::Mech, kf, kr, c)
    T = promote_type(eltype(kf), eltype(kr), typeof(c[1]))
    fr = Vector{T}(undef, length(m.edges))
    br = Vector{T}(undef, length(m.edges))
    for (i, e) in enumerate(m.edges)
        f = T(kf[i])
        for l in e.cons
            f *= c[l]
        end
        r = T(kr[i])
        for l in e.rel
            r *= c[l]
        end
        fr[i] = f
        br[i] = r
    end
    fr, br
end

"""Exact Gaussian elimination (first nonzero pivot) over the rationals; overwrites M and b."""
function gauss_solve!(M::Matrix{QQ}, b::Vector{QQ})
    n = length(b)
    for k in 1:n
        p = findfirst(i -> !iszero(M[i, k]), k:n)
        p === nothing && error("singular system")
        p += k - 1
        if p != k
            for j in 1:n
                M[k, j], M[p, j] = M[p, j], M[k, j]
            end
            b[k], b[p] = b[p], b[k]
        end
        for i in k+1:n
            iszero(M[i, k]) && continue
            f = M[i, k] / M[k, k]
            for j in k:n
                M[i, j] -= f * M[k, j]
            end
            b[i] -= f * b[k]
        end
    end
    x = zeros(QQ, n)
    for i in n:-1:1
        s = b[i]
        for j in i+1:n
            s -= M[i, j] * x[j]
        end
        x[i] = s / M[i, i]
    end
    x
end

"""Exact stationary distribution (sums to 1) of the enzyme-state chain with step rates fr, br."""
function steady_exact(m::Mech, fr::Vector{QQ}, br::Vector{QQ})
    n = length(m.states)
    M = zeros(QQ, n, n)
    for (i, e) in enumerate(m.edges)
        M[e.w, e.u] += fr[i]
        M[e.u, e.u] -= fr[i]
        M[e.u, e.w] += br[i]
        M[e.w, e.w] -= br[i]
    end
    M[n, :] .= one(QQ)
    b = zeros(QQ, n)
    b[n] = one(QQ)
    gauss_solve!(M, b)
end

"""Net production rate of ligand l per unit enzyme."""
function netprod(m::Mech, fr, br, π, l)
    s = zero(eltype(π))
    for (i, e) in enumerate(m.edges)
        ν = count(==(l), e.rel) - count(==(l), e.cons)
        ν == 0 && continue
        s += ν * (fr[i] * π[e.u] - br[i] * π[e.w])
    end
    s
end

"""
Exact turnover per unit enzyme at c (rationals). Net P and Q release and net A and B uptake
must agree exactly; the function errors otherwise.
"""
function rate_exact(m::Mech, kfQ::Vector{QQ}, krQ::Vector{QQ}, c::NTuple{4,QQ})
    fr, br = pseudo(m, kfQ, krQ, c)
    π = steady_exact(m, fr, br)
    vP = netprod(m, fr, br, π, 3)
    vQ = netprod(m, fr, br, π, 4)
    vA = -netprod(m, fr, br, π, 1)
    vB = -netprod(m, fr, br, π, 2)
    (vP == vQ && vP == vA && vP == vB) || error("exact turnover mismatch")
    vP, π
end

qpow10(k::Integer) = k >= 0 ? QQ(big(10)^k) : QQ(1, big(10)^(-k))

"""
Thermodynamic check from the rate constants alone: for each chord of a spanning tree grown from
E, the fundamental cycle's log10 Π(kf/kr) minus n log10 Keq, where n is the cycle's net uptake
of A (the net uptake must be n(A + B) - n(P + Q)). Returns (max deviation, per-cycle rows).
"""
function thermo_check(m::Mech, kf, kr, lK)
    n = length(m.states)
    φ = Vector{Union{Nothing,BigFloat}}(nothing, n)
    λ = [zeros(Int, 4) for _ in 1:n]
    φ[1] = big(0.0)
    intree = falses(length(m.edges))
    function netd(e)
        d = zeros(Int, 4)
        for l in e.cons
            d[l] += 1
        end
        for l in e.rel
            d[l] -= 1
        end
        d
    end
    changed = true
    while changed
        changed = false
        for (i, e) in enumerate(m.edges)
            lk = log10(big(kf[i])) - log10(big(kr[i]))
            if φ[e.u] !== nothing && φ[e.w] === nothing
                φ[e.w] = φ[e.u] + lk
                λ[e.w] = λ[e.u] .+ netd(e)
                intree[i] = changed = true
            elseif φ[e.w] !== nothing && φ[e.u] === nothing
                φ[e.u] = φ[e.w] - lk
                λ[e.u] = λ[e.w] .- netd(e)
                intree[i] = changed = true
            end
        end
    end
    all(x -> x !== nothing, φ) || error("graph not connected")
    worst = big(0.0)
    rows = Tuple{Int,Int,Float64,Float64}[]
    for (i, e) in enumerate(m.edges)
        intree[i] && continue
        lk = log10(big(kf[i])) - log10(big(kr[i]))
        logprod = φ[e.u] + lk - φ[e.w]
        net = λ[e.u] .+ netd(e) .- λ[e.w]
        nA = net[1]
        net == [nA, nA, -nA, -nA] || error("cycle through step $i has net uptake $net")
        dev = abs(logprod - nA * big(lK))
        push!(rows, (i, nA, Float64(logprod), Float64(dev)))
        worst = max(worst, dev)
    end
    worst, rows
end

"""Largest violation (log10 units) of lo ≤ log10 k ≤ hi over all forward and reverse constants."""
boxviol(kf, kr; lo = -6, hi = 6) =
    maximum(max(0.0, Float64(log10(big(k))) - hi, lo - Float64(log10(big(k))))
            for k in vcat(kf, kr))

# ---------------------------------------------------------------------------------------------
# King–Altman: π_r ∝ Σ over spanning trees of Π (rates of the tree's arcs directed toward r).

struct KA
    n::Int
    trees::Vector{Matrix{Int}}   # per root: ntrees × (n - 1), +e = forward arc, -e = backward
end

function combos(m, k)
    out = Vector{Vector{Int}}()
    c = collect(1:k)
    while true
        push!(out, copy(c))
        i = k
        while i >= 1 && c[i] == m - k + i
            i -= 1
        end
        i == 0 && break
        c[i] += 1
        for j in i+1:k
            c[j] = c[j-1] + 1
        end
    end
    out
end

function build_ka(m::Mech)
    n = length(m.states)
    trees = Vector{Vector{Int}}()
    for comb in combos(length(m.edges), n - 1)
        parent = collect(1:n)
        ok = true
        for e in comb
            a = m.edges[e].u
            while parent[a] != a
                a = parent[a]
            end
            b = m.edges[e].w
            while parent[b] != b
                b = parent[b]
            end
            if a == b
                ok = false
                break
            end
            parent[a] = b
        end
        ok && push!(trees, comb)
    end
    mats = Matrix{Int}[]
    for r in 1:n
        M = zeros(Int, length(trees), n - 1)
        for (t, comb) in enumerate(trees)
            adj = [Int[] for _ in 1:n]
            for e in comb
                push!(adj[m.edges[e].u], e)
                push!(adj[m.edges[e].w], e)
            end
            seen = falses(n)
            seen[r] = true
            queue = [r]
            col = 0
            while !isempty(queue)
                x = popfirst!(queue)
                for e in adj[x]
                    ed = m.edges[e]
                    y = ed.u == x ? ed.w : ed.u
                    seen[y] && continue
                    seen[y] = true
                    push!(queue, y)
                    col += 1
                    M[t, col] = ed.u == y ? e : -e
                end
            end
            col == n - 1 || error("orientation failed")
        end
        push!(mats, M)
    end
    KA(n, mats)
end

ntrees(ka::KA) = size(ka.trees[1], 1)

function ka_weights!(w::Vector{Float64}, ka::KA, fr::Vector{Float64}, br::Vector{Float64})
    @inbounds for r in 1:ka.n
        M = ka.trees[r]
        s = 0.0
        for t in 1:size(M, 1)
            p = 1.0
            for j in 1:size(M, 2)
                e = M[t, j]
                p *= e > 0 ? fr[e] : br[-e]
            end
            s += p
        end
        w[r] = s
    end
    w
end

"""
Release flux of product l (3 = P, 4 = Q) per unit enzyme from King–Altman weights, valid when
that product's concentration is zero (its release steps are then one-way, so the flux is a sum
of positive terms and needs no subtraction).
"""
function ka_release(m::Mech, w, fr, l)
    s = 0.0
    @inbounds for (i, e) in enumerate(m.edges)
        (l in e.rel) && (s += fr[i] * w[e.u])
    end
    s / sum(w)
end

# ---------------------------------------------------------------------------------------------
# Textbook (T) closed forms, step numbers as in MECHS["T"]:
# 5 EAB->EPQ, 6 EPQ->EQ+P, 7 EPQ->EP+Q, 8 EQ->E+Q, 9 EP->E+P.

"""Turnover of T at A, B -> ∞, P = 0 and product Q: 1/T(Q) with q = k9/(k9 + k-7 Q),
T(Q) = 1/k5 + (1 + k-5/k5 + k6/k8 + k7 q/k9)/(k6 + k7 q)."""
function vsat_T(kf, kr, Q = 0)
    k5, km5, k6, k7, k8, k9, km7 = kf[5], kr[5], kf[6], kf[7], kf[8], kf[9], kr[7]
    q = k9 / (k9 + km7 * Q)
    1 / (1 / k5 + (1 + km5 / k5 + k6 / k8 + k7 * q / k9) / (k6 + k7 * q))
end

"""Turnover of g (random binding, ordered release) at A, B -> ∞ and P = Q = 0:
1 / (1/k5 + (1 + k-5/k5)/k6 + 1/k7) with 6 EPQ->EQ+P, 7 EQ->E+Q."""
vsat_g(kf, kr) = 1 / (1 / kf[5] + (1 + kr[5] / kf[5]) / kf[6] + 1 / kf[7])

"""
Expected time spent in each state of `S` before leaving S (exact), starting from `start`. All
transitions out of S count as absorption. Returns the vector in the order of S.
"""
function expected_times(m::Mech, fr::Vector{QQ}, br::Vector{QQ}, S, start)
    idx = Dict(s => i for (i, s) in enumerate(S))
    k = length(S)
    G = zeros(QQ, k, k)
    for (i, e) in enumerate(m.edges)
        for (a, b, r) in ((e.u, e.w, fr[i]), (e.w, e.u, br[i]))
            iszero(r) && continue
            haskey(idx, a) || continue
            G[idx[a], idx[a]] -= r
            haskey(idx, b) && (G[idx[a], idx[b]] += r)
        end
    end
    M = Matrix{QQ}(-transpose(G))
    rhs = zeros(QQ, k)
    rhs[idx[start]] = one(QQ)
    gauss_solve!(M, rhs)
end

end
