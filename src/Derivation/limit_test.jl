# ABOUTME: Class 3's limit test mod 2^61 − 1: residue arithmetic, polynomials mod p and the
# ABOUTME: law's limit on a face (`_face_limit_check`), by interpolation or reconstruction.

# ─── Class 3: the limit test mod 2^61 − 1 ───────────────────────────────

"""The prime 2^61 − 1: the limit test computes every rate, weight and limit modulo it."""
const _P61 = (Int64(1) << 61) - 1
"""A limit that grows without bound, outside the residues 0…p − 1."""
const _INFINITE = Int64(-1)
"""A law that has no value on a face, outside the residues 0…p − 1."""
const _UNDEFINED = Int64(-2)

_mulp(a::Int64, b::Int64) = Int64(mod(widemul(a, b), _P61))
_addp(a::Int64, b::Int64) = (s = a + b; s >= _P61 ? s - _P61 : s)
_subp(a::Int64, b::Int64) = (s = a - b; s < 0 ? s + _P61 : s)
_modp(x::Integer) = Int64(mod(x, _P61))

"""`a^k` mod p; a negative `k` raises the inverse."""
function _powp(a::Int64, k::Integer)
    k < 0 && return _powp(_invp(a), -k)
    r = Int64(1)
    while k > 0
        isodd(k) && (r = _mulp(r, a))
        a = _mulp(a, a)
        k >>= 1
    end
    r
end

"""The inverse of `a` mod p."""
_invp(a::Int64) = (a == 0 && throw(DivideError()); _powp(a, _P61 - 2))

"""A deterministic stream of residues in 2…p − 2 (splitmix64), seeded by a mechanism's
content hash and a face, so every process draws the same points."""
mutable struct _ResidueStream
    state::UInt64
end

"""The next residue of `r`."""
function _residue!(r::_ResidueStream)
    r.state += 0x9e3779b97f4a7c15
    z = r.state
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    Int64(2 + (z ⊻ (z >> 31)) % UInt64(_P61 - 3))
end

"""`(A \\ B, det(A))` mod p for square `A`, by Gauss–Jordan elimination: the determinant is
the product of the pivots, its sign flipped by each row swap; `A \\ B` is an empty matrix,
and the determinant 0, when `A` is singular."""
function _solvep(A::Matrix{Int64}, B::Matrix{Int64})
    n = size(A, 1)
    M = hcat(A, B)
    d = Int64(1)
    for j in 1:n
        k = findfirst(i -> M[i, j] != 0, j:n)
        k === nothing && return Matrix{Int64}(undef, 0, 0), Int64(0)
        k += j - 1
        if k != j
            M[[j, k], :] = M[[k, j], :]
            d = _subp(Int64(0), d)
        end
        d = _mulp(d, M[j, j])
        pivot = _invp(M[j, j])
        for c in j:size(M, 2)
            M[j, c] = _mulp(M[j, c], pivot)
        end
        for i in 1:n
            (i == j || M[i, j] == 0) && continue
            f = M[i, j]
            for c in j:size(M, 2)
                M[i, c] = _subp(M[i, c], _mulp(f, M[j, c]))
            end
        end
    end
    M[:, (n + 1):end], d
end

"""The stationary distribution mod p of the chain with off-diagonal rates `R[i, j]`
(`i → j`) and one closed class covering every state, and the determinant of the system that
gives it (the generator's transpose, its last row ones); empty, and 0, when singular."""
function _stationaryp(R::Matrix{Int64})
    n = size(R, 1)
    n == 1 && return Int64[1], Int64(1)
    A = zeros(Int64, n, n)
    for i in 1:n, j in 1:n
        i == j && continue
        A[j, i] = _addp(A[j, i], R[i, j])
        A[i, i] = _subp(A[i, i], R[i, j])
    end
    A[n, :] .= 1
    b = zeros(Int64, n, 1)
    b[n] = 1
    x, d = _solvep(A, b)
    vec(x), d
end

# Polynomials mod p: coefficient vectors, constant term first.
_ptrimp!(a::Vector{Int64}) = (while length(a) > 1 && a[end] == 0; pop!(a); end; a)
_degp(a::Vector{Int64}) =
    (b = _ptrimp!(copy(a)); length(b) == 1 && b[1] == 0 ? -1 : length(b) - 1)
_paddp(a::Vector{Int64}, b::Vector{Int64}) =
    _ptrimp!([_addp(get(a, i, Int64(0)), get(b, i, Int64(0)))
              for i in 1:max(length(a), length(b))])
_psubp(a::Vector{Int64}, b::Vector{Int64}) =
    _ptrimp!([_subp(get(a, i, Int64(0)), get(b, i, Int64(0)))
              for i in 1:max(length(a), length(b))])
function _pmulp(a::Vector{Int64}, b::Vector{Int64})
    c = zeros(Int64, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] = _addp(c[i + j - 1], _mulp(a[i], b[j]))
    end
    _ptrimp!(c)
end
function _pdivremp(a::Vector{Int64}, b::Vector{Int64})
    a, b = _ptrimp!(copy(a)), _ptrimp!(copy(b))
    db = _degp(b)
    db < 0 && throw(DivideError())
    q = zeros(Int64, max(_degp(a) - db + 1, 1))
    lead = _invp(b[end])
    while _degp(a) >= db
        k = _degp(a) - db
        c = _mulp(a[end], lead)
        q[k + 1] = c
        for i in 0:db
            a[i + k + 1] = _subp(a[i + k + 1], _mulp(c, b[i + 1]))
        end
        _ptrimp!(a)
    end
    _ptrimp!(q), a
end
_pevalp(a::Vector{Int64}, x::Int64) =
    foldr((c, acc) -> _addp(c, _mulp(acc, x)), a; init = Int64(0))
"""The lowest-order term of `a` as `(order, coefficient)`; `(typemax(Int), 0)` for 0."""
_plowp(a::Vector{Int64}) =
    (k = findfirst(!=(0), a); k === nothing ? (typemax(Int), Int64(0)) : (k - 1, a[k]))

"""Newton interpolation mod p through the points `(x, y)`."""
function _interpp(x::Vector{Int64}, y::Vector{Int64})
    n = length(x)
    c = copy(y)
    for j in 2:n, i in n:-1:j
        c[i] = _mulp(_subp(c[i], c[i - 1]), _invp(_subp(x[i], x[i - j + 1])))
    end
    p = Int64[c[n]]
    for i in (n - 1):-1:1
        p = _paddp(_pmulp(p, Int64[_subp(Int64(0), x[i]), 1]), Int64[c[i]])
    end
    p
end

"""`(N, D)` with degrees at most `dmax` and `N/D` through the points: the extended
Euclidean algorithm on the interpolant (it needs at least `2 dmax + 1` points)."""
function _ratrecp(x::Vector{Int64}, y::Vector{Int64}, dmax::Int)
    M = Int64[1]
    for xi in x
        M = _pmulp(M, Int64[_subp(Int64(0), xi), 1])
    end
    r0, r1 = M, _interpp(x, y)
    t0, t1 = Int64[0], Int64[1]
    while _degp(r1) > dmax
        q, r = _pdivremp(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, _psubp(t0, _pmulp(q, t1))
    end
    r1, t1
end

"""The rational function of one variable that `f` evaluates mod p, reconstructed from
random points with a degree bound doubled until two more points agree; empty vectors past
degree 1024."""
function _reconstructp(f, rng::_ResidueStream)
    d = 8
    while d <= 1024
        xs = [_residue!(rng) for _ in 1:(2d + 3)]
        ys = [f(x) for x in xs]
        N, D = _ratrecp(xs[1:(end - 2)], ys[1:(end - 2)], d)
        all(k -> (dv = _pevalp(D, xs[k]); dv != 0 &&
                  _mulp(_pevalp(N, xs[k]), _invp(dv)) == ys[k]),
            (length(xs) - 1):length(xs)) && return N, D
        d *= 2
    end
    Int64[], Int64[]
end

"""
A random point mod p of the model of `mech`'s catalytic graph `m` for `ncopies`
conformations on face `face`: log-free energies of every form in each conformation and of
every metabolite name, an integer null-space basis of the kinetic-group ties (members
equal in each conformation that keeps them, an `:EqualAI` member equal across) raised to
random residues; from them each edge's equilibrium constant `K[k, c]` and `L`; a
steady-state speed per kinetic group and conformation (one for both of an `:EqualAI`
group), a rapid-equilibrium speed per edge and conformation, the flip's speed;
concentrations `conc` per name, zero on the face.
"""
function _model_point(mech::_AnyMech, ncopies::Int, face::UInt64, rng::_ResidueStream)
    m = Mechanisms.catalytic(mech)
    nf, names = length(Mechanisms.species(m)), _guard_names(m)
    nnames = length(names)
    nv = ncopies * nf + nnames
    excluded(k, c) = Mechanisms.excludes(mech, k, _CONFORMATIONS[c])
    equal_ai(g) = mech isa Mechanisms.AllostericMechanism &&
                  Mechanisms.tag(mech, g) === :EqualAI
    function drop!(row, k, c, dir)
        a, b = Mechanisms.ends(m, k)
        row[(c - 1) * nf + a] += dir
        row[(c - 1) * nf + b] -= dir
        up, off = Mechanisms.taken_up(m, k), Mechanisms.given_off(m, k)
        for (n, name) in enumerate(names)
            row[ncopies * nf + n] +=
                dir * (Mechanisms.copies(up, name) - Mechanisms.copies(off, name))
        end
        row
    end
    ties = Vector{Int}[]
    for g in 1:Mechanisms.ngroups(m), c in 1:ncopies
        members = Mechanisms.members(m, g)
        lead = first(members)
        excluded(lead, c) && continue
        for k in members[2:end]
            push!(ties, drop!(drop!(zeros(Int, nv), k, c, 1), lead, c, -1))
        end
        if c == 2 && equal_ai(g)
            for k in members
                push!(ties, drop!(drop!(zeros(Int, nv), k, 2, 1), k, 1, -1))
            end
        end
    end
    B = _integer_nullspace(isempty(ties) ? zeros(Int, 0, nv) :
                           permutedims(reduce(hcat, ties)))
    t = [_residue!(rng) for _ in axes(B, 2)]
    x = [reduce(_mulp, (_powp(t[j], B[i, j]) for j in axes(B, 2) if B[i, j] != 0);
                init = Int64(1)) for i in 1:nv]
    K = ones(Int64, length(Mechanisms.steps(m)), ncopies)
    xnames = view(x, (ncopies * nf + 1):nv)
    for k in Mechanisms.steps(m), c in 1:ncopies
        a, b = Mechanisms.ends(m, k)
        num = _mul_held(x[(c - 1) * nf + a], Mechanisms.taken_up(m, k), names, xnames)
        den = _mul_held(x[(c - 1) * nf + b], Mechanisms.given_off(m, k), names, xnames)
        K[k, c] = _mulp(num, _invp(den))
    end
    ss = [_residue!(rng) for _ in 1:Mechanisms.ngroups(m), _ in 1:ncopies]
    for g in 1:Mechanisms.ngroups(m)
        ncopies == 2 && equal_ai(g) && (ss[g, 2] = ss[g, 1])
    end
    re = [_residue!(rng) for _ in Mechanisms.steps(m), _ in 1:ncopies]
    flip = _residue!(rng)
    L = ncopies == 2 ? _mulp(x[1], _invp(x[nf + 1])) : Int64(1)
    conc = [_zeroes(face, n) ? Int64(0) : _residue!(rng) for n in 1:nnames]
    (; K, ss, re, flip, L, conc)
end

"""`r` times `value[n]` once for each copy of name `n` of `names` (`_guard_names`) that
`held` holds, mod p."""
function _mul_held(r::Int64, held::Mechanisms.Metabolites, names::Vector{Symbol},
                   value::AbstractVector{Int64})
    for (n, x) in enumerate(names), _ in 1:Mechanisms.copies(held, x)
        r = _mulp(r, value[n])
    end
    r
end

"""`pt` with new rapid-equilibrium and flip speeds: another speed pattern."""
_respeed_point(pt, rng::_ResidueStream) =
    merge(pt, (re = [_residue!(rng) for _ in axes(pt.re, 1), _ in axes(pt.re, 2)],
               flip = _residue!(rng)))

"""The rate mod p of every transition of `net` at the concentrations `conc` (per name):
its speed times `K` in the stored direction, its speed alone in reverse, times each
metabolite it takes up, a fast one times `Λ` unless `Λ` is 0."""
function _rates_p(m::Mechanisms.Mechanism, net::_Network, pt, conc::Vector{Int64}, Λ::Int64)
    names = _guard_names(m)
    r = zeros(Int64, length(net.from))
    for t in eachindex(net.from)
        k, c = net.edge[t], net.conformation[t]
        if k == 0
            r[t] = net.forward[t] ? _mulp(pt.flip, pt.L) : pt.flip
        else
            speed = net.fast[t] ? pt.re[k, c] : pt.ss[Mechanisms.group(m, k), c]
            r[t] = net.forward[t] ? _mulp(speed, pt.K[k, c]) : speed
            held = net.forward[t] ? Mechanisms.taken_up(m, k) : Mechanisms.given_off(m, k)
            r[t] = _mul_held(r[t], held, names, conc)
        end
        net.fast[t] && Λ != 0 && (r[t] = _mulp(r[t], Λ))
    end
    r
end

"""ρ(u) mod p: over the slow transitions out of node `u`, their rate `r` times their
uptake of `x1` plus the fall of its potential `h` along them."""
function _slow_rho(net::_Network, r::Vector{Int64}, h::Vector{Int}, x1::Int)
    rho = zeros(Int64, net.nforms * net.ncopies)
    for t in eachindex(net.from)
        (net.fast[t] || r[t] == 0) && continue
        u = net.from[t]
        rho[u] = _addp(rho[u],
                       _mulp(r[t], _modp(net.change[x1, t] + h[u] - h[net.to[t]])))
    end
    rho
end

"""Equilibrium weights mod p within each part of `parts` over the fast transitions inside
it, from its first node, at rates `r`."""
function _eq_weights_p(net::_Network, parts::Vector{Vector{Int}}, r::Vector{Int64})
    nn = net.nforms * net.ncopies
    part = zeros(Int, nn)
    for (q, P) in enumerate(parts), u in P
        part[u] = q
    end
    next = [Int[] for _ in 1:nn]
    for t in eachindex(net.from)
        u = net.from[t]
        net.fast[t] && part[u] > 0 && part[u] == part[net.to[t]] && push!(next[u], t)
    end
    w = zeros(Int64, nn)
    for P in parts
        w[P[1]] = 1
        done = falses(nn)
        done[P[1]] = true
        todo = [P[1]]
        while !isempty(todo)
            u = pop!(todo)
            for t in next[u]
                v = net.to[t]
                done[v] && continue
                back = isodd(t) ? t + 1 : t - 1
                r[t] != 0 && r[back] != 0 ||
                    error("_face_limit_check: a fast transition inside a part runs one way")
                w[v] = _mulp(w[u], _mulp(r[t], _invp(r[back])))
                done[v] = true
                push!(todo, v)
            end
        end
        all(u -> done[u], P) || error("_face_limit_check: a part is not fast-connected")
    end
    w
end

"""The face chain's limit mod p as every rapid-equilibrium speed grows, by interpolating
the rate's numerator and its tree-sum denominator in Λ (degree below the closed class's
size) on the face chain's closed class; `_INFINITE` when the rate grows with Λ."""
function _limit_by_interpolation(m::Mechanisms.Mechanism, net::_Network, s, pt,
                                 h::Vector{Int}, x1::Int, rng::_ResidueStream)
    C = s.closed[1]
    n = length(C)
    pos = zeros(Int, net.nforms * net.ncopies)
    pos[C] .= 1:n
    xs, Ds, Ns = Int64[], Int64[], Int64[]
    for _ in 1:(n + 1)
        Λ = _residue!(rng)
        r = _rates_p(m, net, pt, pt.conc, Λ)
        Q = zeros(Int64, n, n)
        for t in eachindex(net.from)
            a, b = pos[net.from[t]], pos[net.to[t]]
            s.alive[t] && r[t] != 0 && a > 0 && b > 0 && a != b || continue
            Q[a, b] = _addp(Q[a, b], r[t])
        end
        x, D = _stationaryp(Q)
        isempty(x) && error("_face_limit_check: singular at a random Λ")
        rho = _slow_rho(net, r, h, x1)
        v = reduce(_addp, (_mulp(x[pos[u]], rho[u]) for u in C); init = Int64(0))
        push!(xs, Λ)
        push!(Ds, D)
        push!(Ns, _mulp(D, v))
    end
    cD, cN = _interpp(xs[1:n], Ds[1:n]), _interpp(xs[1:n], Ns[1:n])
    _pevalp(cD, xs[end]) == Ds[end] && _pevalp(cN, xs[end]) == Ns[end] ||
        error("_face_limit_check: the degree bound failed")
    dD, dN = _degp(cD), _degp(cN)
    dN > dD && return _INFINITE
    dN < dD && return Int64(0)
    _mulp(cN[dD + 1], _invp(cD[dD + 1]))
end

"""The law mod p at concentrations `conc` with every name present: segment weights and
the slow segment chain over every segment."""
function _law_inside(m::Mechanisms.Mechanism, net::_Network, s, pt, conc::Vector{Int64},
                     h::Vector{Int}, x1::Int)
    r = _rates_p(m, net, pt, conc, Int64(0))
    w = _eq_weights_p(net, s.segments, r)
    W = [_invp(reduce(_addp, w[S]; init = Int64(0))) for S in s.segments]
    Q = zeros(Int64, length(s.segments), length(s.segments))
    for t in eachindex(net.from)
        net.fast[t] && continue
        a, b = s.segment[net.from[t]], s.segment[net.to[t]]
        a == b && continue
        Q[a, b] = _addp(Q[a, b], _mulp(_mulp(w[net.from[t]], W[a]), r[t]))
    end
    mu = first(_stationaryp(Q))
    isempty(mu) && error("_face_limit_check: singular segment chain")
    rho = _slow_rho(net, r, h, x1)
    v = Int64(0)
    for (g, S) in enumerate(s.segments), u in S
        v = _addp(v, _mulp(_mulp(mu[g], _mulp(w[u], W[g])), rho[u]))
    end
    v
end

"""The law's value mod p on the face: its limit as the zeroed concentrations go to zero
along two random directions, by rational reconstruction in their scale; `_UNDEFINED` when
a direction meets a pole or the two differ."""
function _law_by_reconstruction(m::Mechanisms.Mechanism, net::_Network, s, pt,
                                h::Vector{Int}, x1::Int, face::UInt64,
                                rng::_ResidueStream)
    found = Int64[]
    for _ in 1:2
        direction = [_residue!(rng) for _ in pt.conc]
        law_at(δ) = _law_inside(m, net, s, pt,
                                [_zeroes(face, n) ? _mulp(direction[n], δ) : pt.conc[n]
                                 for n in eachindex(pt.conc)], h, x1)
        num, den = _reconstructp(law_at, rng)
        isempty(den) && error("_face_limit_check: no reconstruction under the degree cap")
        (a, lead_num), (b, lead_den) = _plowp(num), _plowp(den)
        push!(found, a > b ? Int64(0) : a < b ? _INFINITE :
                     _mulp(lead_num, _invp(lead_den)))
    end
    found[1] == _INFINITE || found[1] != found[2] ? _UNDEFINED : found[1]
end

"""
    _face_limit_check(mech, net, s, face, seed) -> _LimitVerdict

The limit test of class 3 on face `face` of `mech`, whose law
network is `net` and face structure `s`: at two random points of the model mod 2^61 − 1,
drawn from `seed`, the fast-RE limit of the face chain by interpolation in Λ against the
law's face value by reconstruction as the zeroed concentrations vanish. On a difference a
third limit at new speeds tells `_LIMIT_DIFFER` from `_LIMIT_SPEED`. A false
`_LIMIT_EQUAL` has probability at most twice the degree over p.
"""
function _face_limit_check(mech::_AnyMech, net::_Network, s, face::UInt64,
                           seed::UInt64)
    m = Mechanisms.catalytic(mech)
    rng = _ResidueStream(seed)
    x1 = _turnover_name(m)
    # Class 2, which runs first, leaves no fast cycle that takes `x1` up: a potential.
    h = s.offset[x1, :]
    limit(pt) = _limit_by_interpolation(m, net, s, pt, h, x1, rng)
    for _ in 1:2
        pt = _model_point(mech, net.ncopies, face, rng)
        law = _law_by_reconstruction(m, net, s, pt, h, x1, face, rng)
        law == _UNDEFINED && return _LIMIT_LAW_NAN
        v = limit(pt)
        # The face chain's fast transitions are the law network's, a subset of the
        # ε-network's, where class 2, which runs first, finds no turnover cycle.
        v == _INFINITE && error("_face_limit_check: the limit grows without bound")
        v == law && continue
        return limit(_respeed_point(pt, rng)) == v ? _LIMIT_DIFFER : _LIMIT_SPEED
    end
    _LIMIT_EQUAL
end
