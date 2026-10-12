# ABOUTME: Skeptic of probe P2: checks written independently of Thermo.jl's machinery (exact
# ABOUTME: constraint dimension from potentials, FD rank, step fluxes, scrambles); stdlib only.
const PROTO = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/proto"
isdefined(Main, :Oracle) || include(joinpath(PROTO, "common", "Oracle.jl"))
isdefined(Main, :Thermo) || include(joinpath(PROTO, "common", "Thermo.jl"))

module Adv
using LinearAlgebra, Random
import ..Pop, ..Oracle, ..Thermo

const QB = Rational{BigInt}

# ─── exact linear algebra ─────────────────────────────────────────────
function rrefq(A::AbstractMatrix)
    M = QB.(A); nr, nc = size(M); piv = Int[]; r = 1
    for c in 1:nc
        r > nr && break
        p = findfirst(i -> M[i, c] != 0, r:nr)
        p === nothing && continue
        p += r - 1
        M[[r, p], :] = M[[p, r], :]
        M[r, :] .//= M[r, c]
        for i in 1:nr
            (i == r || M[i, c] == 0) && continue
            M[i, :] .-= M[i, c] .* M[r, :]
        end
        push!(piv, c); r += 1
    end
    M, piv
end
rankq(A) = size(A, 1) == 0 ? 0 : length(rrefq(A)[2])
function nullq(A)
    nc = size(A, 2)
    size(A, 1) == 0 && return Matrix{QB}(I, nc, nc)
    M, piv = rrefq(A)
    free = setdiff(1:nc, piv)
    N = zeros(QB, nc, length(free))
    for (k, f) in enumerate(free)
        N[f, k] = 1
        for (i, p) in enumerate(piv); N[p, k] = -M[i, f]; end
    end
    N
end

# ─── constraint dimension from potentials (independent of cycles and of Thermo) ──────
"""Bare constants in the STORED orientation, oracle convention: RE `[K]`, SS `[kf, kr]`."""
xcols(m) = [(j, r) for (j, s) in enumerate(m.steps) for r in (s.re ? (1,) : (1, 2))]
bindinglike(s) = length(s.consumed) + length(s.released) == 1

"""Linear map θ → log x with θ = (form potentials G, metabolite potentials μ by name, one
log prefactor per SS step), and the Keq functional f(θ) = Σμ_S − Σμ_P."""
function potential_map(m)
    n = length(m.species)
    mets = unique(Symbol[m.substrates; m.products;
                         [x.name for s in m.steps for x in [s.consumed; s.released]]])
    mi = Dict(x => n + i for (i, x) in enumerate(mets))
    ss = [j for (j, s) in enumerate(m.steps) if !s.re]
    si = Dict(j => n + length(mets) + i for (i, j) in enumerate(ss))
    nθ = n + length(mets) + length(ss)
    cols = xcols(m)
    A = zeros(QB, length(cols), nθ)
    for (row, (j, r)) in enumerate(cols)
        s = m.steps[j]
        ka = zeros(QB, nθ)                      # log Ka(from → to) in the stored direction
        ka[s.from] += 1; ka[s.to] -= 1
        for x in s.consumed; ka[mi[x.name]] += 1; end
        for x in s.released; ka[mi[x.name]] -= 1; end
        if s.re
            # oracle: a stored binding's K is its dissociation constant; any other RE K is
            # Ka in the stored direction (a stored release's K is then also a Kd)
            A[row, :] .= Pop.is_binding(s) ? -ka : ka
        else
            A[row, si[j]] = 1
            r == 2 && (A[row, :] .-= ka)
        end
    end
    f = zeros(QB, nθ)
    for x in m.substrates; f[mi[x]] += 1; end
    for x in m.products; f[mi[x]] -= 1; end
    A, f, cols
end

"""Equality rows of the kinetic groups over `xcols(m)`. `rel(j)` is the lead each member is
tied to and whether it runs reversed relative to that lead in stored orientation. Default:
the package's semantics (every member tied to the first, as stored)."""
function tie_rows(m, cols; rel = j -> nothing)
    idx = Dict(c => i for (i, c) in enumerate(cols))
    rows = Vector{QB}[]
    for g in m.groups, j in g
        lead, flip = something(rel(j), (g[1], false))
        j == lead && continue
        s = m.steps[j]
        if s.re
            row = zeros(QB, length(cols)); row[idx[(j, 1)]] += 1
            row[idx[(lead, 1)]] -= (flip && !bindinglike(s)) ? -1 : 1
            push!(rows, row)
        else
            for r in 1:2
                row = zeros(QB, length(cols))
                row[idx[(j, flip ? 3 - r : r)]] += 1; row[idx[(lead, r)]] -= 1
                push!(rows, row)
            end
        end
    end
    isempty(rows) ? zeros(QB, 0, length(cols)) : permutedims(reduce(hcat, rows))
end

"""The probe's group semantics as relative stored flips: member j is tied to lead l, reversed
when (canonical flip of j) ⊻ (alignment flip) ⊻ (canonical flip of l)."""
function probe_rel(m)
    cs = Thermo.canonical_steps(m)
    al = Thermo.group_alignment(m, cs, m.groups)
    j -> (l = al[j][1]; (l, cs[j].flipped ⊻ al[j][2] ⊻ cs[l].flipped))
end

"""Independent count (the dimension of the admissible log-constants at fixed Keq), whether
the groups contradict thermodynamics (force Keq = 1), and the AllParameters count."""
function true_count(m; semantics = :package)
    A, f, cols = potential_map(m)
    T = tie_rows(m, cols; rel = semantics === :probe ? probe_rel(m) : (j -> nothing))
    TA = T * A
    contradiction = rankq(vcat(TA, permutedims(f))) == rankq(TA)
    N = nullq(vcat(TA, permutedims(f)))
    (count = rankq(A * N), contradiction = contradiction, allp = length(cols) - rankq(T))
end

"""Whether stored per-step constants `k` respect the package's ties (equal stored constants)."""
function ties_hold(m, k; rel = j -> nothing, tol = big(10)^-60)
    for g in m.groups, j in g
        lead, flip = something(rel(j), (g[1], false))
        j == lead && continue
        s = m.steps[j]
        if s.re
            a, b = k[j][1], k[lead][1]
            b = (flip && !bindinglike(s)) ? 1 / b : b
            abs(a - b) > tol * abs(b) && return false
        else
            for r in 1:2
                a, b = k[j][flip ? 3 - r : r], k[lead][r]
                abs(a - b) > tol * abs(b) && return false
            end
        end
    end
    true
end

# ─── numerics with the oracle ───────────────────────────────────────────
"""Stored constants from random independents of `sol` (and the Keq used)."""
function random_constants(sol, rng; keq = Oracle.logu(rng))
    θ = Dict(i => Oracle.logu(rng) for i in sol.indep)
    Thermo.stored_constants(sol, Thermo.all_values(sol, θ, keq)), keq, θ
end

"""Rank of the Jacobian of v(c_1..c_N) over log independents: central differences in 512-bit
BigFloat, QR with column pivoting, rows scaled by |v|. Returns (rank, smallest kept and
largest dropped |R_ii| relative to the largest), the threshold being 1e-20."""
function fd_rank(m, sol, rng; h = big(10)^-28, Λ = nothing)
    setprecision(BigFloat, 512) do
        p = length(sol.indep)
        p == 0 && return 0, BigFloat[]
        npts = 3p + 8
        keq = Oracle.logu(rng)
        θ = Dict(i => Oracle.logu(rng) for i in sol.indep)
        cs = [Dict{Symbol, BigFloat}(x => Oracle.logu(rng) for x in m.metabolites) for _ in 1:npts]
        J = zeros(BigFloat, npts, p)
        for (col, i) in enumerate(sol.indep), sg in (1, -1)
            θp = copy(θ); θp[i] *= exp(sg * h)
            k = Thermo.stored_constants(sol, Thermo.all_values(sol, θp, keq))
            for r in 1:npts
                J[r, col] += sg * Oracle.steady_state_v(m, k, cs[r]; Λ) / (2h)
            end
        end
        k0 = Thermo.stored_constants(sol, Thermo.all_values(sol, θ, keq))
        for r in 1:npts                      # row scaling by |v| (rank-invariant)
            v = Oracle.steady_state_v(m, k0, cs[r]; Λ); J[r, :] ./= abs(v)
        end
        F = qr(J, ColumnNorm())
        d = sort!(abs.(diag(F.R)); rev = true); d ./= d[1]
        r = count(>(big(10)^-20), d)
        kept = Float64(d[r]); dropped = r < p ? Float64(d[r + 1]) : 0.0
        r, kept, dropped
    end
end

"""Net steady-state flux of every SS step (`nothing` for RE) and v, in BigFloat."""
function step_fluxes(m, k, c)
    setprecision(BigFloat, 512) do
        x = Oracle.steady_state_populations(m, k, c)
        cc(v) = prod((BigFloat(c[y.name]) for y in v); init = big(1))
        J = map(enumerate(m.steps)) do (j, s)
            s.re ? nothing : k[j][1] * cc(s.consumed) * x[s.from] - k[j][2] * cc(s.released) * x[s.to]
        end
        J, Oracle.steady_state_v(m, k, c)
    end
end

"""Plain bridges of the full step graph (no contraction): removing the step disconnects
its two forms."""
function plain_bridges(m)
    n = length(m.species)
    map(eachindex(m.steps)) do e
        adj = [Int[] for _ in 1:n]
        for (j, s) in enumerate(m.steps)
            j == e && continue
            push!(adj[s.from], s.to); push!(adj[s.to], s.from)
        end
        seen = falses(n); q = [m.steps[e].from]; seen[q[1]] = true
        while !isempty(q)
            u = popfirst!(q)
            for v in adj[u]; seen[v] || (seen[v] = true; push!(q, v)); end
        end
        !seen[m.steps[e].to]
    end
end

# ─── scrambles ───────────────────────────────────────────────────────────
"""`m` with every storage choice scrambled: species order and names, bound and residual
list order, step order and direction, consumed/released order, group and member order,
and the order of the substrate/product/metabolite lists."""
function scramble(m, rng)
    n = length(m.species); π = randperm(rng, n); inv = invperm(π)
    species = [Pop.PSpecies(shuffle(rng, sp.bound), sp.conformation, shuffle(rng, sp.res_added),
                            shuffle(rng, sp.res_subtracted), Symbol("X", rand(rng, 1:10^9)))
               for sp in m.species[π]]
    flip = rand(rng, Bool, length(m.steps))
    σ = randperm(rng, length(m.steps)); σinv = invperm(σ)
    steps = map(σ) do j
        s = m.steps[j]
        c, r = shuffle(rng, s.consumed), shuffle(rng, s.released)
        flip[j] ? Pop.PStep(inv[s.to], inv[s.from], r, c, s.re) :
                  Pop.PStep(inv[s.from], inv[s.to], c, r, s.re)
    end
    groups = [shuffle(rng, [σinv[j] for j in g]) for g in shuffle(rng, m.groups)]
    Pop.PMech(; id = m.id, source = m.source, substrates = shuffle(rng, m.substrates),
              products = shuffle(rng, m.products), inhibitors = shuffle(rng, m.inhibitors),
              metabolites = shuffle(rng, m.metabolites), species, steps, groups)
end

signature(sol) = (Thermo.independent_names(sol), Thermo.dependent_table(sol),
                  sort!([Thermo.colname(sol, d.col) => (d.fold, d.label) for d in sol.deps];
                        by = first))

"""Exponent statistics of a solve: max |numerator|/denominator met, non-integer exponents."""
function exponent_stats(sol)
    nonint = 0; maxe = 0 // 1; keqs = Rational{Int}[]
    for d in sol.deps
        d.keq == 0 || push!(keqs, d.keq)
        for (_, e) in d.terms
            isinteger(e) || (nonint += 1); abs(e) > maxe && (maxe = abs(e))
        end
        isinteger(d.keq) || (nonint += 1)
    end
    (nonint = nonint, maxexp = maxe, maxint = sol.maxint, keqs = keqs)
end

end # module
