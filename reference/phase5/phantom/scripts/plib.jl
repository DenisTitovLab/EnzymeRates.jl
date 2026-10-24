# ABOUTME: L-phantom research library on the slice types: census text to slice mechanisms, laws as
# ABOUTME: branch sums, modular and BigFloat Jacobian ranks, the rescaling certificate, keys, identity.

const SP = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
isdefined(Main, :Slice) || include(joinpath(SP, "slice", "src", "Slice.jl"))
const S = Slice
using Random, Printf, SHA, LinearAlgebra

const QQ = Rational{BigInt}

# ─── census neutral text → slice mechanism (copied from research5/unit/ulib.jl:16-98) ────

function parse_label(lab::AbstractString, reactants)
    if endswith(lab, "mod") && Symbol(lab[1:(end - 3)]) in reactants
        return S.Met(Symbol(lab[1:(end - 3)]), true)
    end
    x = Symbol(lab)
    S.Met(x, !(x in reactants))
end

function parse_form(t::AbstractString, reactants)
    i = findfirst('(', t)
    conf, inner = Symbol(t[1:(i - 1)]), t[(i + 1):(end - 1)]
    res = ""
    occursin("; residual=", inner) && ((inner, res) = split(inner, "; residual="))
    bound = isempty(inner) ? S.Met[] : [parse_label(l, reactants) for l in split(inner, ",")]
    added, subtracted = Symbol[], Symbol[]
    for r in (isempty(res) ? SubString{String}[] : split(res, "_"))
        startswith(r, "plus") ? push!(added, Symbol(r[5:end])) :
            push!(subtracted, Symbol(r[6:end]))
    end
    S.RawForm(conf, bound, added, subtracted)
end

function parse_side(t::AbstractString, reactants)
    parts = split(t, " + ")
    parse_form(parts[1], reactants), S.Met[parse_label(l, reactants) for l in parts[2:end]]
end

"""The slice mechanism of a census neutral text; `plain = true` drops the allosteric tags,
multiplicity and sites (the catalytic mechanism)."""
function parse_neutral(text::AbstractString, subs::Vector{Symbol}, prods::Vector{Symbol};
                       plain::Bool = false)
    reactants = Set([subs; prods])
    allo = startswith(text, "allosteric")
    mult, body, sitetext = 1, "", "none"
    if allo
        rest = text[length("allosteric n=") + 1:end]
        sp = findfirst(' ', rest)
        mult = parse(Int, rest[1:(sp - 1)])
        body, sitetext = split(rest[(sp + 1):end], " sites "; limit = 2)
    else
        body = text[length("plain ") + 1:end]
    end
    steps, group, tags = S.RawStep[], Int[], Symbol[]
    for (g, gt) in enumerate(split(body, " / "))
        tag = ""
        if !endswith(gt, "}")
            gt, tag = split(gt, "}:")
            gt *= "}"
        end
        allo && push!(tags, Symbol(tag))
        for st in split(gt[2:(end - 1)], " | ")
            flag, core = st[(end - 1):end], st[1:(end - 3)]
            a, b = split(core, " <-> ")
            (fa, ca), (fb, cb) = parse_side(a, reactants), parse_side(b, reactants)
            push!(steps, S.RawStep(fa, fb, ca, cb, flag == "RE"))
            push!(group, g)
        end
    end
    sites = S.Site[]
    if allo && sitetext != "none"
        for s in split(sitetext, " / ")
            mm = match(r"^site\((\d+): (.*)\)$", s)
            ligs = [split(l, "::") for l in split(mm[2], ", ")]
            push!(sites, S.Site(Symbol.(first.(ligs)), Symbol.(last.(ligs)),
                                parse(Int, mm[1])))
        end
    end
    allo && !plain ? S.Mechanism(subs, prods, steps, group; tags, multiplicity = mult, sites) :
                     S.Mechanism(subs, prods, steps, group)
end

# ─── a prime field ──────────────────────────────────────────────────────────────────────

"""An element of the prime field of order `p` (a prime below 2^63)."""
struct Fp{p} <: Number
    v::UInt64
    Fp{p}(v::UInt64, ::Val{:raw}) where {p} = new{p}(v)
end
Fp{p}(x::Integer) where {p} = Fp{p}(UInt64(mod(x, p)), Val(:raw))
Base.:+(a::Fp{p}, b::Fp{p}) where {p} = (s = a.v + b.v; Fp{p}(s >= p ? s - p : s, Val(:raw)))
Base.:-(a::Fp{p}, b::Fp{p}) where {p} = Fp{p}(a.v >= b.v ? a.v - b.v : a.v + p - b.v, Val(:raw))
Base.:-(a::Fp{p}) where {p} = Fp{p}(a.v == 0 ? a.v : p - a.v, Val(:raw))
Base.:*(a::Fp{p}, b::Fp{p}) where {p} = Fp{p}(UInt64(widemul(a.v, b.v) % p), Val(:raw))
function Base.:^(a::Fp{p}, e::Integer) where {p}
    e < 0 && return inv(a)^(-e)
    r, b = one(Fp{p}), a
    while e > 0
        isodd(e) && (r *= b)
        b *= b
        e >>= 1
    end
    r
end
Base.inv(a::Fp{p}) where {p} = (a.v == 0 && throw(DivideError()); a^(p - 2))
Base.:/(a::Fp{p}, b::Fp{p}) where {p} = a * inv(b)
Base.zero(::Type{Fp{p}}) where {p} = Fp{p}(UInt64(0), Val(:raw))
Base.one(::Type{Fp{p}}) where {p} = Fp{p}(UInt64(1), Val(:raw))
Base.iszero(a::Fp) = a.v == 0
Base.:(==)(a::Fp{p}, b::Fp{p}) where {p} = a.v == b.v
Base.convert(::Type{Fp{p}}, x::Integer) where {p} = Fp{p}(x)

"""Deterministic Miller–Rabin for 64-bit `n`."""
function isprime64(n::UInt64)
    n < 2 && return false
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
        n % q == 0 && return n == q
    end
    d, s = n - 1, 0
    while iseven(d); d >>= 1; s += 1; end
    pw(a, e) = (r = UInt128(1); b = UInt128(a % n); while e > 0; isodd(e) && (r = r * b % n);
                b = b * b % n; e >>= 1; end; r)
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
        x = pw(UInt64(a), d)
        (x == 1 || x == n - 1) && continue
        ok = false
        for _ in 1:(s - 1)
            x = x * x % n
            x == n - 1 && (ok = true; break)
        end
        ok || return false
    end
    true
end
const P1 = UInt64(2)^61 - 1
const P2 = UInt64(2)^62 - 57
@assert isprime64(P1) && isprime64(P2)
const F1 = Fp{P1}
const F2 = Fp{P2}

# ─── laws as branch sums ────────────────────────────────────────────────────────────────

"""A polynomial compiled for evaluation: per term its coefficient and `(variable,
exponent)` pairs."""
struct CPoly
    coef::Vector{BigInt}
    mono::Vector{Vector{Pair{Int, Rational{Int}}}}
end
CPoly(p::S.Poly) = CPoly([big(c) for c in values(p.terms)], [collect(m) for m in keys(p.terms)])

"""A law `v = E_total·N/D` over `vars`, `N` and `D` sums of branches, each branch a list of
`(polynomial, power)`."""
struct BLaw
    vars::S.Vars
    polys::Vector{CPoly}
    N::Vector{Vector{Tuple{Int, Int}}}
    D::Vector{Vector{Tuple{Int, Int}}}
end

function blaw(L::S.Law)
    BLaw(L.vars, CPoly.([L.H, L.Q, L.D]), [[(1, 1), (2, 1)]], [[(3, 1)]])
end

"""The MWC law as branches (`AlloLaw` docstring, slice/src/allosteric.jl:4-13)."""
function blaw(A::S.AlloLaw)
    T = eltype(values(A.H.terms))
    QI, DI = S.inactive_terms(A)
    Lp = S.term(T, S.mono(A.vars.index[:L] => 1))
    ps = S.Poly[A.H, A.fI, A.active.Q, A.active.D, Lp, A.fA, QI, DI]
    WA = Tuple{Int, Int}[]
    WI = Tuple{Int, Int}[]
    for (wa, wi, m) in A.sites
        push!(ps, wa); push!(WA, (length(ps), m))
        push!(ps, wi); push!(WI, (length(ps), m))
    end
    n = A.n
    N = [[(1, 1); (2, n); (3, 1); (4, n - 1); WA], [(1, 1); (5, 1); (6, n); (7, 1); (8, n - 1); WI]]
    D = [[(2, n); (4, n); WA], [(5, 1); (6, n); (8, n); WI]]
    BLaw(A.vars, CPoly.(ps), N, D)
end

"""The least common multiple of the exponent denominators of every non-concentration
variable in `B`."""
denlcm(B::BLaw) = foldl(lcm, (denominator(e) for p in B.polys for m in p.mono for (v, e) in m
                              if !B.vars.conc[v]); init = 1)

"""
    point(B, x, cols, Q) -> (N, D, dlogN, dlogD) per column

`B` at the variable values `x[v] = t[v]^Q` given as bases `t` (concentrations: `t = x`),
with each column's log-derivative numerator `Σ_term (e·Q)·term` for the variables `cols`.
"""
function point(B::BLaw, t::Vector{T}, cols::Vector{Int}, Q::Int) where {T}
    nc = length(cols)
    at = Dict(v => j for (j, v) in enumerate(cols))
    pv = Vector{T}(undef, length(B.polys))
    pg = Vector{Vector{T}}(undef, length(B.polys))
    cache = Dict{Tuple{Int, Int}, T}()
    pw(v, k) = T === Float64 ? t[v]^k : get!(() -> t[v]^k, cache, (v, k))
    used = falses(length(B.polys))
    for b in Iterators.flatten((B.N, B.D)), (i, k) in b
        k > 0 && (used[i] = true)
    end
    for (i, p) in enumerate(B.polys)
        used[i] || continue
        F = zero(T)
        G = zeros(T, nc)
        for (c, mo) in zip(p.coef, p.mono)
            term = T(c)
            for (v, e) in mo
                k = B.vars.conc[v] ? Int(e) : Int(e * Q)
                term *= pw(v, k)
            end
            F += term
            for (v, e) in mo
                j = get(at, v, 0)
                j == 0 || (G[j] += T(Int(e * Q)) * term)
            end
        end
        pv[i], pg[i] = F, G
    end
    function branches(bs)
        V = zero(T)
        G = zeros(T, nc)
        for b in bs
            bv = one(T)
            bg = zeros(T, nc)
            zero_branch = false
            for (i, k) in b
                k == 0 && continue
                F = pv[i]
                iszero(F) && (zero_branch = true; break)
                bv *= F^k
                bg .+= T(k) .* pg[i] ./ F
            end
            zero_branch && continue
            V += bv
            G .+= bv .* bg
        end
        V, G
    end
    N, Ng = branches(B.N)
    D, Dg = branches(B.D)
    N, D, Ng, Dg
end

"""The Jacobian of `log v` over `cols` (columns scaled by `Q`) at `nrow` random concentration
points, the constants fixed at random; `ones` appends the relative mode's scale column."""
function jacobian(::Type{T}, B::BLaw, cols, Q, nrow, rng; ones::Bool = false) where {T}
    t = Vector{T}(undef, length(B.vars.names))
    draw() = T <: Fp ? T(rand(rng, 1:typemax(Int))) : T(10)^(2 * rand(rng) - 1)
    for v in eachindex(t)
        B.vars.conc[v] && continue
        t[v] = T <: Fp ? draw() : draw()^(1 / big(Q))
    end
    J = Matrix{T}(undef, nrow, length(cols) + ones)
    for r in 1:nrow
        foreach(v -> B.vars.conc[v] && (t[v] = draw()), eachindex(t))
        N, D, Ng, Dg = point(B, t, cols, Q)
        (iszero(N) || iszero(D)) && error("zero law value at a random point")
        J[r, 1:length(cols)] .= Ng ./ N .- Dg ./ D
        ones && (J[r, end] = one(T))
    end
    J
end

"""Rank over a prime field by Gaussian elimination."""
function rank_fp(J::Matrix{T}) where {T <: Fp}
    A = copy(J)
    m, n = size(A)
    r = 0
    for c in 1:n
        p = findfirst(i -> !iszero(A[i, c]), (r + 1):m)
        p === nothing && continue
        p += r
        r += 1
        A[[r, p], :] = A[[p, r], :]
        iv = inv(A[r, c])
        for i in (r + 1):m
            iszero(A[i, c]) && continue
            f = A[i, c] * iv
            for k in c:n
                A[i, k] -= f * A[r, k]
            end
        end
        r == m && break
    end
    r
end

"""Rank of a BigFloat matrix by complete pivoting: pivots above `tol` times the first; the
smallest kept and largest dropped pivot ratios (census/small/lib.jl:471-492)."""
function rank_bf(J; tol = big(10)^-30)
    A = copy(J); m, n = size(A); piv = BigFloat[]
    for k in 1:min(m, n)
        _, idx = findmax(abs.(A[k:m, k:n]))
        i, j = idx[1] + k - 1, idx[2] + k - 1
        A[[k, i], :] = A[[i, k], :]; A[:, [k, j]] = A[:, [j, k]]
        push!(piv, abs(A[k, k]))
        iszero(A[k, k]) && break
        for r in (k + 1):m
            A[r, k:n] .-= (A[r, k] / A[k, k]) .* A[k, k:n]
        end
    end
    (isempty(piv) || iszero(piv[1]) || !isfinite(piv[1])) && return 0, 0.0, 0.0
    rk = count(q -> q > tol * piv[1], piv)
    rk, Float64(piv[rk] / piv[1]), rk < length(piv) ? Float64(piv[rk + 1] / piv[1]) : 0.0
end

"""The modular rank of `B` over the variables `cols`: the larger of one draw in each of two
prime fields (each draw a lower bound that is exact unless a random point is special)."""
function rank_mod(B::BLaw, cols, rng; ones::Bool = false)
    isempty(cols) && return 0
    Q = denlcm(B)
    nrow = length(cols) + ones + 6
    max(rank_fp(jacobian(F1, B, cols, Q, nrow, rng; ones)),
        rank_fp(jacobian(F2, B, cols, Q, nrow, rng; ones))) - ones
end

"""The BigFloat rank of `B` over `cols` (192 bits, `2n + 6` points, pivot cut 1e-30: the
census's `slice_rank`, census/small/lib.jl:494-529)."""
function rank_big(B::BLaw, cols, rng)
    isempty(cols) && return (0, 1.0, 0.0)
    setprecision(BigFloat, 192) do
        rank_bf(jacobian(BigFloat, B, cols, denlcm(B), 2 * length(cols) + 6, rng))
    end
end

"""The variable indices of `names` in `B.vars` (error when one is missing)."""
colsof(B::BLaw, names) = [B.vars.index[x] for x in names]

# ─── the key (census/small/lib.jl:262-310, concentrations by name) ──────────────────────

function canonical_named(L)
    comps, units, header = S.key_parts(L)
    vars = L.vars
    fixed(v) = vars.conc[v] || vars.names[v] in (:Keq, :L)
    terms = Tuple{Int, S.Mono, BigInt}[]
    for u in units
        s = S.lowest_shift(comps[u])
        g = foldl(gcd, (c for i in u for c in values(comps[i].terms)); init = big(0))
        append!(terms, ((i, S.mmul(mo, s), div(big(c), g)) for i in u
                        for (mo, c) in comps[i].terms))
    end
    params = [v for v in eachindex(vars.names) if !fixed(v)]
    held = Set(v for (_, mo, _) in terms for (v, _) in mo)
    present = filter(in(held), params)
    at = Dict(v => j for (j, v) in enumerate(present))
    nm(v) = vars.conc[v] ? "c_" * string(vars.names[v]) : string(vars.names[v])
    labels = [string(i, "|", c, "|", join(("$(nm(v))^$e" for (v, e) in mo if fixed(v)),
                                          "*")) for (i, mo, c) in terms]
    inc = [[(at[v], e) for (v, e) in mo if !fixed(v)] for (_, mo, _) in terms]
    occ = [Tuple{Int, S.Exp}[] for _ in present]
    for (t, ps) in enumerate(inc), (p, e) in ps
        push!(occ[p], (t, e))
    end
    best, order, leaves = nothing, Int[], 0
    function search(pcol, tcol)
        pcol, tcol = S.refine(pcol, tcol, inc, occ)
        if allunique(pcol)
            (leaves += 1) > 100_000 && error("canonical labelling: over 100000 leaves")
            text = S.encode(pcol, inc, labels)
            (best === nothing || text < best) && ((best, order) = (text, sortperm(pcol)))
            return
        end
        c = minimum(x for x in pcol if count(==(x), pcol) > 1)
        for p in findall(==(c), pcol)
            q = 2 .* pcol
            q[p] -= 1
            search(q, tcol)
        end
    end
    search(ones(Int, length(present)), S.ranks(labels))
    absent = sort!(vars.names[setdiff(params, present)])
    text = string(header, "|P=", length(params), "|absent=", length(absent), "\n", best)
    (key = bytes2hex(sha256(text))[1:16], nabsent = length(absent),
     order = [vars.names[present[order]]; absent])
end

# ─── the rescaling certificate ──────────────────────────────────────────────────────────

ppow(p::S.Poly{T}, k::Int) where {T} = foldl(*, (p for _ in 1:k); init = one(S.Poly{T}))

"""
    certificate(A) -> weights or nothing

For a one-subunit MWC law whose inactive conformation does not catalyze (`QI = 0`), write
`v = N/(Y + L·X)` with `N = H·fI·QA·WA`, `Y = fI·DA·WA`, `X = fA·DI·WI`. When every term of
`X` is a term of `Y` with the same coefficient, `Y + L·X = (Y - X) + (1 + L)·X`; a weight
`w` per constant (fitted constants and regulator constants; concentrations, `Keq` and `L`
weigh 0) and a shift `s` with `e·w = s` on every constant monomial of `N` and `Y - X` and
`e·w = s + 1` on every one of `X` give `v(θ, L) = v₀(θ ⊙ (1 + L)^w)`, `v₀ = N/Y` the
L → 0 law, for all positive `θ` and `L`: the variant's family equals the L → 0 law's.
Returns `w` (name => weight) or `nothing`.
"""
function certificate(A::S.AlloLaw)
    A.n == 1 || return nothing
    QI, DI = S.inactive_terms(A)
    iszero(QI) || return nothing
    T = eltype(values(A.H.terms))
    WA = foldl(*, (ppow(s[1], s[3]) for s in A.sites); init = one(S.Poly{T}))
    WI = foldl(*, (ppow(s[2], s[3]) for s in A.sites); init = one(S.Poly{T}))
    N = A.H * A.fI * A.active.Q * WA
    Y = A.fI * A.active.D * WA
    X = A.fA * DI * WI
    all(((m, c),) -> get(Y.terms, m, 0) == c, X.terms) || return nothing
    vars = A.vars
    θ = [v for v in eachindex(vars.names) if !vars.conc[v] && !(vars.names[v] in (:Keq, :L))]
    at = Dict(v => j for (j, v) in enumerate(θ))
    nθ = length(θ)
    rows = Dict{Vector{QQ}, QQ}()
    function add!(p, rhs)
        for m in keys(p.terms)
            e = zeros(QQ, nθ + 1)
            for (v, x) in m
                j = get(at, v, 0)
                j == 0 || (e[j] += x)
            end
            e[end] = -1
            r = get(rows, e, nothing)
            r === nothing ? (rows[e] = rhs) : r == rhs || return false
        end
        true
    end
    add!(N, QQ(0)) && add!(Y - X, QQ(0)) && add!(X, QQ(1)) || return nothing
    M = zeros(QQ, length(rows), nθ + 2)
    for (i, (e, rhs)) in enumerate(rows)
        M[i, 1:(nθ + 1)] .= e
        M[i, end] = rhs
    end
    piv = S.rref!(M, nθ + 1)
    all(iszero, M[(length(piv) + 1):end, end]) || return nothing
    w = zeros(QQ, nθ + 1)
    for (r, c) in enumerate(piv)
        w[c] = M[r, end]
    end
    Dict(vars.names[θ[j]] => w[j] for j in 1:nθ)
end

# ─── evaluation in BigFloat, by name ────────────────────────────────────────────────────

"""`v/E_total` of `B` at values by name (`vals`, BigFloat; a variable `B` holds but `vals`
lacks is an error)."""
function value(B::BLaw, vals::Dict{Symbol, BigFloat})
    Q = denlcm(B)
    t = [B.vars.conc[v] ? vals[B.vars.names[v]] : vals[B.vars.names[v]]^(1 / big(Q))
         for v in eachindex(B.vars.names)]
    N, D, _, _ = point(B, t, Int[], Q)
    N / D
end

"""`v/E_total` of `B` over a prime field at bases by name (`x = t^Q`)."""
function value_fp(B::BLaw, t::Dict{Symbol, T}, Q::Int) where {T <: Fp}
    tv = [t[x] for x in B.vars.names]
    N, D, _, _ = point(B, tv, Int[], Q)
    N / D
end
