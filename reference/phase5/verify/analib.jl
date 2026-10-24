# ABOUTME: Verifier analysis library on the stock slice: BigFloat identifiability rank of a law,
# ABOUTME: the exact L-rescaling certificate, exact rational elimination, and the unit rule.
using Random, LinearAlgebra
const QQ = Rational{BigInt}

# ── value and log-gradient (θ ∂/∂θ) of polynomials at a point ─────────────────────────────
struct VG
    v::BigFloat
    g::Vector{BigFloat}
end
Base.:+(a::VG, b::VG) = VG(a.v + b.v, a.g .+ b.g)
Base.:*(a::VG, b::VG) = VG(a.v * b.v, a.v .* b.g .+ b.v .* a.g)
vgpow(a::VG, k::Int) = k == 0 ? VG(big(1.0), zero(a.g)) : VG(a.v^k, k * a.v^(k - 1) .* a.g)

"""`p` at the variable values `x` (logs `lx`), with its gradient over the columns `col[v]`."""
function vg(p, x, lx, col, ncol)
    v, g = big(0.0), zeros(BigFloat, ncol)
    for (mo, c) in p.terms
        t = BigFloat(c)
        for (u, e) in mo
            t *= denominator(e) == 1 ? x[u]^Int(numerator(e)) : exp(BigFloat(e) * lx[u])
        end
        v += t
        for (u, e) in mo
            j = get(col, u, 0)
            j == 0 || (g[j] += BigFloat(e) * t)
        end
    end
    VG(v, g)
end

"""Rank of `J` by complete pivoting: pivots above `tol` times the largest entry; also the
smallest kept and the largest dropped pivot over that entry."""
function numrank(J::Matrix{BigFloat}; tol = big(10.0)^-40)
    A = copy(J)
    m, n = size(A)
    mx = maximum(abs, A; init = big(0.0))
    mx == 0 && return 0, big(0.0), big(0.0)
    r, kept = 0, big(Inf)
    for k in 1:min(m, n)
        best, bi, bj = big(-1.0), 0, 0
        for i in k:m, j in k:n
            a = abs(A[i, j])
            a > best && ((best, bi, bj) = (a, i, j))
        end
        best <= tol * mx && return r, kept / mx, best / mx
        A[[k, bi], :] = A[[bi, k], :]
        A[:, [k, bj]] = A[:, [bj, k]]
        r += 1
        kept = min(kept, best)
        for i in (k + 1):m
            f = A[i, k] / A[k, k]
            f == 0 && continue
            A[i, k:n] .-= f .* A[k, k:n]
        end
    end
    r, kept / mx, big(0.0)
end

"""Random values over `vars`: constants and Keq log-uniform in 10^±1, concentrations 10^±1.5."""
function draw(vars, rng)
    [big(10.0)^(vars.conc[i] ? 3 * rand(rng) - 1.5 : 2 * rand(rng) - 1) for i in eachindex(vars.names)]
end

"""Log-derivative rows of `v` over the fitted names `fit`: one θ draw, `nrow` concentration
draws. `parts(x, lx, col, ncol)` returns `(num, den)` as `VG`s."""
function jac(vars, fit, parts, nrow, rng)
    col = Dict{Int, Int}()
    for (j, f) in enumerate(fit)
        i = get(vars.index, f, 0)
        i == 0 || (col[i] = j)
    end
    base = draw(vars, rng)
    J = zeros(BigFloat, nrow, length(fit))
    for r in 1:nrow
        x = copy(base)
        for i in eachindex(x)
            vars.conc[i] && (x[i] = big(10.0)^(3 * rand(rng) - 1.5))
        end
        lx = log.(x)
        num, den = parts(x, lx, col, length(fit))
        J[r, :] .= num.g ./ num.v .- den.g ./ den.v
    end
    J
end

plain_parts(L) = (x, lx, col, nc) -> (vg(L.Q, x, lx, col, nc), vg(L.D, x, lx, col, nc))

"""The MWC numerator and denominator of an `AlloLaw` (slice allosteric.jl:10-11)."""
function allo_parts(S, A)
    QI, DI = S.inactive_terms(A)
    li = A.vars.index[:L]
    function (x, lx, col, nc)
        f(p) = vg(p, x, lx, col, nc)
        g = zeros(BigFloat, nc)
        j = get(col, li, 0)
        j == 0 || (g[j] = x[li])
        Lv = VG(x[li], g)
        n = A.n
        fA, fI, qa, da, qi, di = f(A.fA), f(A.fI), f(A.active.Q), f(A.active.D), f(QI), f(DI)
        wa = foldl(*, (vgpow(f(s[1]), s[3]) for s in A.sites); init = VG(big(1.0), zeros(BigFloat, nc)))
        wi = foldl(*, (vgpow(f(s[2]), s[3]) for s in A.sites); init = VG(big(1.0), zeros(BigFloat, nc)))
        num = vgpow(fI, n) * qa * vgpow(da, n - 1) * wa + Lv * vgpow(fA, n) * qi * vgpow(di, n - 1) * wi
        den = vgpow(fI, n) * vgpow(da, n) * wa + Lv * vgpow(fA, n) * vgpow(di, n) * wi
        num, den
    end
end

"""Identifiability rank over `fit`, the larger of two θ draws (each `length(fit) + 8` rows);
also whether the draws agreed and the worst pivot gap."""
function law_rank(vars, fit, parts, rng)
    r1, k1, d1 = numrank(jac(vars, fit, parts, length(fit) + 8, rng))
    r2, k2, d2 = numrank(jac(vars, fit, parts, length(fit) + 8, rng))
    max(r1, r2), r1 == r2, min(k1, k2), max(d1, d2)
end

# ── exact rational elimination ──────────────────────────────────────────────────────────
"""Solve `M[:, 1:end-1] * w = M[:, end]` exactly; `nothing` when inconsistent."""
function solve_exact(M::Matrix{QQ})
    A = copy(M)
    m, n = size(A, 1), size(A, 2) - 1
    piv, r = Int[], 0
    for c in 1:n
        i = findfirst(k -> !iszero(A[k, c]), (r + 1):m)
        i === nothing && continue
        i += r
        r += 1
        A[[r, i], :] = A[[i, r], :]
        A[r, :] .//= A[r, c]
        for k in 1:m
            k != r && !iszero(A[k, c]) && (A[k, :] .-= A[k, c] .* A[r, :])
        end
        push!(piv, c)
        r == m && break
    end
    all(iszero, A[(r + 1):m, end]) || return nothing
    w = zeros(QQ, n)
    for (k, c) in enumerate(piv)
        w[c] = A[k, end]
    end
    w
end

"""The L-rescaling certificate of an n = 1 MWC law with QI = 0 (phantom.md section 1,
derived again): weights `w` over the constants θ (not Keq, L or concentrations) and a scalar
`s` with e·w = s on the terms of N = fI·QA·WA and of R = Y − X, e·w = s + 1 on X, where
Y = fI·DA·WA and X = fA·DI·WI ⊂ Y with equal coefficients. Returns a symbol saying why
not, or the weights."""
function certificate(S, A)
    A.n == 1 || return :n
    QI, DI = S.inactive_terms(A)
    iszero(QI) || return :catalyzes
    P = typeof(A.H)
    WA = foldl(*, (foldl(*, fill(s[1], s[3])) for s in A.sites); init = one(P))
    WI = foldl(*, (foldl(*, fill(s[2], s[3])) for s in A.sites); init = one(P))
    N, Y, X = A.fI * A.active.Q * WA, A.fI * A.active.D * WA, A.fA * DI * WI
    all(((mo, c),) -> get(Y.terms, mo, 0) == c, X.terms) || return :not_subsum
    R = Y - X
    θ = [v for v in eachindex(A.vars.names)
         if !A.vars.conc[v] && !(A.vars.names[v] in (:Keq, :L))]
    at = Dict(v => j for (j, v) in enumerate(θ))
    rows = Vector{QQ}[]
    for (p, rhs) in ((N, 0), (R, 0), (X, 1)), mo in keys(p.terms)
        e = zeros(QQ, length(θ) + 2)
        for (v, x) in mo
            j = get(at, v, 0)
            j == 0 || (e[j] += x)
        end
        e[end - 1] = -1
        e[end] = rhs
        push!(rows, e)
    end
    unique!(rows)
    M = permutedims(reduce(hcat, rows))
    w = solve_exact(M)
    w === nothing ? :no_solution : Dict(A.vars.names[θ[j]] => w[j] for j in eachindex(θ))
end
