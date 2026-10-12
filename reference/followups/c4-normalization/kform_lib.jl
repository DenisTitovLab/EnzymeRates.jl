# ABOUTME: Card C4 follow-up library: per-metabolite K-forms of a P4 law's denominator over exact
# ABOUTME: rational arithmetic (affine basis of the concentration support, γ relation tests, N and K_X).
include(joinpath(@__DIR__, "..", "..", "proto", "p4", "P4.jl"))
using Random, Printf
using .Poly: Mono, mono, mmul, mpow, mexp, P, nterms

const PI = P{Int}
const QQ = Rational{BigInt}

"""Exact value of `p` at the point `x` (indexed by variable)."""
evalq(p::PI, x::Vector{QQ}) =
    sum((QQ(c) * prod((x[v]^e for (v, e) in mo); init = one(QQ)) for (mo, c) in p.t); init = zero(QQ))

"""Random positive integer point for every variable of the law."""
randpoint(L, rng) = QQ[rand(rng, 1:10^6) for _ in L.vars.names]

"""Exponent vector of a concentration monomial over the metabolites `mets`."""
ev(L, mo::Mono, mets) = [mexp(mo, L.conc[x]) for x in mets]

"""Coefficients `λ` with Σ λ_b B[b] = v, or `nothing` when `v` is outside the span of the (linearly
independent) vectors `B`."""
function solve_span(B::Vector{Vector{QQ}}, v::Vector{QQ})
    n, k = length(v), length(B)
    M = QQ[j <= k ? B[j][i] : v[i] for i in 1:n, j in 1:(k + 1)]
    piv, r = Int[], 1
    for c in 1:k
        r > n && break
        p = findfirst(i -> M[i, c] != 0, r:n)
        p === nothing && continue
        p += r - 1
        M[[r, p], :] = M[[p, r], :]
        M[r, :] ./= M[r, c]
        for i in 1:n
            (i != r && M[i, c] != 0) && (M[i, :] .-= M[i, c] .* M[r, :])
        end
        push!(piv, c); r += 1
    end
    any(M[i, k + 1] != 0 for i in r:n) && return nothing
    λ = zeros(QQ, k)
    for (i, c) in enumerate(piv); λ[c] = M[i, k + 1]; end
    λ
end

"""Inverse of a square rational matrix (error when singular)."""
function invq(A::Matrix{QQ})
    n = size(A, 1)
    M = [A Matrix{QQ}(I_n(n))]
    for c in 1:n
        p = findfirst(i -> M[i, c] != 0, c:n); p === nothing && error("singular")
        p += c - 1
        M[[c, p], :] = M[[p, c], :]
        M[c, :] ./= M[c, c]
        for i in 1:n
            (i != c && M[i, c] != 0) && (M[i, :] .-= M[i, c] .* M[c, :])
        end
    end
    M[:, (n + 1):end]
end
I_n(n) = QQ[i == j ? 1 : 0 for i in 1:n, j in 1:n]

"""Sort key of a concentration monomial: degree, then how many non-substrate metabolites it
holds, then its text (so the constant term, then substrate monomials, come first)."""
monokey(L, mo) = (P4.degree(mo), count(((v, _),) -> !(L.vars.names[v] in L.m.substrates), mo),
                  P4.fmt_mono(mo, L.vars))

"""Affine analysis of the concentration support of `D` (concentration monomial => coefficient):
the basis (in `monokey` order), each other monomial's relation λ over the basis, and each
relation's class at random points: `:one` (γ ≡ 1), `:number` (a constant ≠ 1) or `:ratio`
(depends on the constants). γ_j = d_j / Π_b d_b^λ_b, tested as γ^q for the common denominator q."""
function kform(L, D::Dict{Mono, PI}; rng = Xoshiro(7), npts = 3)
    mets = [x for x in L.m.metabolites if any(mo -> mexp(mo, L.conc[x]) != 0, keys(D))]
    mons = sort(collect(keys(D)); by = mo -> monokey(L, mo))
    w(mo) = QQ[1; ev(L, mo, mets)]
    basis, B = Mono[], Vector{QQ}[]
    rels = Tuple{Mono, Vector{QQ}}[]
    for mo in mons
        λ = isempty(B) ? nothing : solve_span(B, w(mo))
        λ === nothing ? (push!(basis, mo); push!(B, w(mo))) : push!(rels, (mo, λ))
    end
    pts = [randpoint(L, rng) for _ in 1:npts]
    dv = [Dict(mo => evalq(D[mo], x) for mo in mons) for x in pts]
    cls = map(rels) do (mo, λ)
        q = foldl(lcm, (denominator(l) for l in λ); init = BigInt(1))
        g = [dv[i][mo]^q / prod((dv[i][b]^(Int(numerator(l * q))) for (b, l) in zip(basis, λ));
                                init = one(QQ)) for i in 1:npts]
        allequal(g) ? (g[1] == 1 ? :one : :number) : :ratio
    end
    # N and K_X: log d_b = log N − Σ_X e_bX log K_X on the basis (unique when square)
    nN = 1 + length(mets)
    N_K = nothing
    if length(basis) == nN
        A = QQ[j == 1 ? 1 : -ev(L, b, mets)[j - 1] for b in basis, j in 1:nN]
        N_K = invq(A)          # row i: exponents over the basis coefficients of log N (i = 1), log K_X
    end
    n_determined = solve_span(B, QQ[1; zeros(QQ, length(mets))]) !== nothing
    (mets = mets, basis = basis, rels = rels, cls = cls, N_K = N_K, n_determined = n_determined,
     rank = length(basis), nmons = length(mons))
end

"""A product of coefficient symbols c(m)^e as text, from exponents over the basis."""
function cprod(L, basis, ex)
    nu, de = String[], String[]
    for (b, e) in zip(basis, ex)
        e == 0 && continue
        s = "c(" * (isempty(b) ? "1" : P4.fmt_mono(b, L.vars)) * ")"
        a = abs(e)
        t = a == 1 ? s : string(s, "^", denominator(a) == 1 ? string(numerator(a)) : string(a))
        push!(e > 0 ? nu : de, t)
    end
    n = isempty(nu) ? "1" : join(nu, "·")
    isempty(de) ? n : string(n, "/", length(de) == 1 ? de[1] : "(" * join(de, "·") * ")")
end

"""The rational function Π d_b^ex_b over the basis coefficients, reduced (integer exponents only)."""
function cratio(D, basis, ex; extra = nothing)
    nu, de = one(PI), one(PI)
    extra === nothing || (nu = nu * D[extra])
    for (b, e) in zip(basis, ex)
        e == 0 && continue
        denominator(e) == 1 || return nothing
        k = Int(numerator(e))
        f = foldl(*, fill(D[b], abs(k)); init = one(PI))
        k > 0 ? (nu = nu * f) : (de = de * f)
    end
    P4.reduce_ratio(nu, de)
end

ratio_text(L, ab) = ab === nothing ? "(fractional power)" :
    (ab[2] == one(PI) ? P4.fmt_factored(ab[1], L.vars) :
     P4.fmt_factored(ab[1], L.vars) * " / (" * P4.fmt_factored(ab[2], L.vars) * ")")

byid(id) = (pi = split(id, "/"); Pop.load_population(pi[1])[parse(Int, pi[2])])
