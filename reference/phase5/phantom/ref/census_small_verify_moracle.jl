# ABOUTME: Mass-action RE-limit oracle written for this verification: exact-rational consistent random
# ABOUTME: constants, two-conformation networks (flip at the free form), :OnlyA limit families, faces.
setprecision(BigFloat, 512)
const RQ = Rational{BigInt}

"""Reduced row echelon form in place; the pivot columns."""
function rref!(A::Matrix{RQ})
    m, n = size(A)
    piv = Int[]
    r = 1
    for c in 1:n
        r > m && break
        p = findfirst(i -> A[i, c] != 0, r:m)
        p === nothing && continue
        p += r - 1
        if p != r
            A[[r, p], :] = A[[p, r], :]
        end
        A[r, :] ./= A[r, c]
        for i in 1:m
            i != r && A[i, c] != 0 && (A[i, :] .-= A[i, c] .* A[r, :])
        end
        push!(piv, c)
        r += 1
    end
    piv
end

"""Exact solution of `A x = b`: (feasible, particular solution, nullspace basis columns)."""
function solve_exact(A::Matrix{RQ}, b::Vector{RQ})
    m, n = size(A)
    M = hcat(A, reshape(b, m, 1))
    piv = rref!(M)
    (n + 1) in piv && return (false, zeros(RQ, n), zeros(RQ, n, 0))
    x0 = zeros(RQ, n)
    for (r, c) in enumerate(piv)
        x0[c] = M[r, n + 1]
    end
    free = setdiff(1:n, piv)
    N = zeros(RQ, n, length(free))
    for (k, f) in enumerate(free)
        N[f, k] = 1
        for (r, c) in enumerate(piv)
            N[c, k] = -M[r, f]
        end
    end
    (true, x0, N)
end
rq(rng) = RQ(rand(rng, -16:16), 8)
stackrows(rows, n) = isempty(rows) ? zeros(RQ, 0, n) : permutedims(reduce(hcat, rows))
dmu(s, mu) = sum((mu[x] for x in s.cons); init = RQ(0)) - sum((mu[x] for x in s.rel); init = RQ(0))

"""Active-state constants consistent with thermodynamics and the kinetic groups: form
potentials `G`, metabolite potentials `mu` (exact rationals), each step's `lnK = G(a) - G(b)
+ Δμ`, equal within a group; one steady-state speed per group."""
function aconsts(ms, rng)
    nf, nm = length(ms.forms), length(ms.mets)
    midx = Dict(x => i for (i, x) in enumerate(ms.mets))
    function row(s)
        v = zeros(RQ, nf + nm)
        v[s.a] += 1; v[s.b] -= 1
        foreach(x -> v[nf + midx[x]] += 1, s.cons)
        foreach(x -> v[nf + midx[x]] -= 1, s.rel)
        v
    end
    rows = Vector{RQ}[]
    for g in 1:ms.ngroups
        ss = [s for s in ms.steps if s.g == g]
        for j in 2:length(ss)
            push!(rows, row(ss[1]) .- row(ss[j]))
        end
    end
    ok, x0, N = solve_exact(stackrows(rows, nf + nm), zeros(RQ, length(rows)))
    ok || error("no consistent active constants")
    x = x0 .+ N * [rq(rng) for _ in 1:size(N, 2)]
    mu = Dict(y => x[nf + i] for (i, y) in enumerate(ms.mets))
    (G = x[1:nf], mu = mu, lnK = [sum(row(s) .* x) for s in ms.steps],
     speed = [big(lu10(rng, -1, 1)) for _ in 1:ms.ngroups])
end

"""Forms the free form reaches over the steps `keep`."""
function reach_from(ms, flip, keep)
    nf = length(ms.forms)
    r = falses(nf); r[flip] = true
    changed = true
    while changed
        changed = false
        for (k, s) in enumerate(ms.steps)
            keep[k] && r[s.a] != r[s.b] && (r[s.a] = r[s.b] = true; changed = true)
        end
    end
    r
end

"""Inactive-state constants. `family = false`: the package's model, the non-:OnlyA steps on
the forms the free form reaches over them. `family = true`: every form and step, :OnlyA plain
bindings with their association constant divided by `exp(ell·e_s)` (affinity → 0), other
:OnlyA steps unconstrained (their ratio left free), `e_s` from `expo` (default 1). Returns
`nothing` when the constraints are infeasible, else the I-form set, per-step lnK (`nothing`
for an absent step) and a sampler of the free directions."""
function iconsts(ms, A, rng; family = false, ell = RQ(60), expo = nothing, pin = nothing)
    nf = length(ms.forms)
    flip = ms.frees[1]
    kept = [ms.tags[s.g] !== :OnlyA for s in ms.steps]
    inI = family ? trues(nf) : reach_from(ms, flip, kept)
    iforms = findall(inI)
    pos = Dict(f => i for (i, f) in enumerate(iforms))
    present = [inI[s.a] && inI[s.b] && (family || kept[k]) for (k, s) in enumerate(ms.steps)]
    neq = sort!(unique([ms.steps[k].g for k in findall(present)
                        if ms.tags[ms.steps[k].g] === :NonequalAI]))
    tpos = Dict(g => length(iforms) + i for (i, g) in enumerate(neq))
    nv = length(iforms) + length(neq)
    rows = Vector{RQ}[]; rhs = RQ[]
    nb = 0
    for (k, s) in enumerate(ms.steps)
        present[k] || continue
        v = zeros(RQ, nv); v[pos[s.a]] += 1; v[pos[s.b]] -= 1
        tag = ms.tags[s.g]
        if tag === :EqualAI
            push!(rows, v); push!(rhs, A.G[s.a] - A.G[s.b])
        elseif tag === :NonequalAI
            v[tpos[s.g]] -= 1; push!(rows, v); push!(rhs, -dmu(s, A.mu))
        elseif plain_binding(ms, s)
            nb += 1
            e = expo === nothing ? 1 : expo[nb]
            push!(rows, v); push!(rhs, A.G[s.a] - A.G[s.b] - ell * e)
        end
    end
    a = zeros(RQ, nv); a[pos[flip]] = 1
    push!(rows, a); push!(rhs, RQ(0))
    if pin !== nothing                  # the package model's values on its own forms and groups
        for (f, val) in pin.G
            r = zeros(RQ, nv); r[pos[f]] = 1; push!(rows, r); push!(rhs, val)
        end
        for (g, val) in pin.t
            haskey(tpos, g) || continue
            r = zeros(RQ, nv); r[tpos[g]] = 1; push!(rows, r); push!(rhs, val)
        end
    end
    ok, x0, N = solve_exact(stackrows(rows, nv), rhs)
    ok || return nothing
    draw(r) = x0 .+ N * [rq(r) for _ in 1:size(N, 2)]
    lnk(x) = [present[k] ? x[pos[s.a]] - x[pos[s.b]] + dmu(s, A.mu) : nothing
              for (k, s) in enumerate(ms.steps)]
    ispeed = [ms.tags[g] === :EqualAI ? A.speed[g] : big(lu10(rng, -1, 1)) for g in 1:ms.ngroups]
    pins(x) = (G = Dict(f => x[pos[f]] for f in iforms), t = Dict(g => x[tpos[g]] for g in neq))
    (iforms = iforms, pos = pos, present = present, kept = kept, draw = draw, lnk = lnk,
     pins = pins, speed = ispeed, nfree = size(N, 2), nbind = nb)
end

"""Terminal strongly connected classes of a directed graph."""
function nterminal(n, arcs)
    adj = [Int[] for _ in 1:n]; radj = [Int[] for _ in 1:n]
    for (a, b) in arcs
        push!(adj[a], b); push!(radj[b], a)
    end
    seen = falses(n); order = Int[]
    for s in 1:n
        seen[s] && continue
        stack = [(s, 1)]; seen[s] = true
        while !isempty(stack)
            u, i = stack[end]
            if i <= length(adj[u])
                stack[end] = (u, i + 1)
                w = adj[u][i]
                seen[w] || (seen[w] = true; push!(stack, (w, 1)))
            else
                pop!(stack); push!(order, u)
            end
        end
    end
    comp = zeros(Int, n); c = 0
    for s in reverse(order)
        comp[s] == 0 || continue
        c += 1; comp[s] = c; st = [s]
        while !isempty(st)
            u = pop!(st)
            for w in radj[u]
                comp[w] == 0 && (comp[w] = c; push!(st, w))
            end
        end
    end
    out = falses(c)
    for (a, b) in arcs
        comp[a] != comp[b] && (out[comp[a]] = true)
    end
    count(!, out)
end

"""Steady-state net consumption of the first substrate per enzyme. A plain network when `I`
is `nothing`; else the two-conformation n = 1 network joined at the free form by a fast flip
of ratio `L`, the inactive steps from `I` with lnK `lnkI` (absent where `nothing`); an absent
or :OnlyA steady-state step in the family model runs at `slow` times its speed. RE steps at
`Λ·σ[k]`. NaN when the steady state is not unique."""
function oracle_rate(ms, A, I, lnkI, L, conc, Λ, σ; slow = big(1e-24))
    nf = length(ms.forms)
    nI = I === nothing ? 0 : length(I.iforms)
    n = nf + nI
    M = zeros(BigFloat, n, n)
    arcs = Tuple{Int, Int}[]
    fl = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
    x1 = ms.subs[1]
    cc(xs) = prod((conc[x] for x in xs); init = big(1.0))
    function add!(a, b, rf, rr, w)
        M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
        rf > 0 && push!(arcs, (a, b)); rr > 0 && push!(arcs, (b, a))
        w != 0 && push!(fl, (a, b, rf, rr, w))
    end
    wt(s) = count(==(x1), s.cons) - count(==(x1), s.rel)
    for (k, s) in enumerate(ms.steps)
        base = s.re ? Λ * σ[k] : A.speed[s.g]
        h = exp(BigFloat(A.lnK[k]) / 2)
        add!(s.a, s.b, base * h * cc(s.cons), base / h * cc(s.rel), wt(s))
    end
    if I !== nothing
        for (k, s) in enumerate(ms.steps)
            lnkI[k] === nothing && continue
            base = s.re ? Λ * σ[k] : I.speed[s.g] * (I.kept[k] ? big(1) : slow)
            h = exp(BigFloat(lnkI[k]) / 2)
            add!(nf + I.pos[s.a], nf + I.pos[s.b], base * h * cc(s.cons), base / h * cc(s.rel), wt(s))
        end
        f = ms.frees[1]
        add!(f, nf + I.pos[f], Λ * σ[end] * sqrt(L), Λ * σ[end] / sqrt(L), 0)
    end
    nterminal(n, arcs) == 1 || return big(NaN)
    Mx = copy(M); Mx[n, :] .= 1
    rhs = zeros(BigFloat, n); rhs[n] = 1
    x = Mx \ rhs
    sum((w * (rf * x[a] - rr * x[b]) for (a, b, rf, rr, w) in fl); init = big(0.0))
end

"""The faces the census checks: every product at zero, every substrate at zero, each single
reactant at zero (when its side has more than one)."""
function faces(ms)
    fs = Tuple{String, Vector{Symbol}}[("all products", ms.prods), ("all substrates", ms.subs)]
    for side in (ms.subs, ms.prods)
        length(side) > 1 && foreach(x -> push!(fs, ("$x alone", [x])), side)
    end
    fs
end

rel(a, b) = (isnan(a) || isnan(b)) ? big(Inf) : abs(a - b) / max(abs(a), abs(b), big(1e-30))

"""Face test of one mechanism with one draw of consistent constants: per face, the exact-face
RE limit under two speed patterns (Λ = 1e30, 1e40) and the interior limit (the face's names at
ε·u with ε = 1e-25, Λ = 1e60) under both patterns. A face fails when the exact-face values
depend on the pattern or differ from the interior limit."""
function face_test(ms, rng; pattern = (-3, 3))
    A = aconsts(ms, rng)
    I = ms.allo ? iconsts(ms, A, rng) : nothing
    lnkI = I === nothing ? nothing : I.lnk(I.draw(rng))
    L = big(lu10(rng, -1, 1))
    conc = Dict{Symbol, BigFloat}(x => big(lu10(rng, -1, 1)) for x in ms.mets)
    nE = length(ms.steps) + 1
    σs = [[big(lu10(rng, pattern...)) for _ in 1:nE] for _ in 1:2]
    pos = oracle_rate(ms, A, I, lnkI, L, conc, big(10)^40, σs[1])
    pos2 = oracle_rate(ms, A, I, lnkI, L, conc, big(10)^40, σs[2])
    out = NamedTuple[]
    for (lab, zs) in faces(ms)
        c0 = copy(conc); foreach(z -> c0[z] = big(0), zs)
        u = Dict(z => big(lu10(rng, -1, 1)) for z in zs)
        ce = copy(conc); foreach(z -> ce[z] = big(10)^-25 * u[z], zs)
        ex = [oracle_rate(ms, A, I, lnkI, L, c0, big(10)^Λe, σ) for σ in σs for Λe in (30, 40)]
        it = [oracle_rate(ms, A, I, lnkI, L, ce, big(10)^60, σ) for σ in σs]
        conv = max(rel(ex[1], ex[2]), rel(ex[3], ex[4]))
        pat = rel(ex[2], ex[4])
        ipat = rel(it[1], it[2])
        gap = max(rel(ex[2], it[1]), rel(ex[4], it[2]))
        push!(out, (face = lab, exact = Float64.(ex[[2, 4]]), interior = Float64.(it),
                    conv = Float64(conv), pattern = Float64(pat), ipattern = Float64(ipat),
                    gap = Float64(gap), fails = pat > 1e-8 || gap > 1e-8))
    end
    (positive_pattern = Float64(rel(pos, pos2)), faces = out)
end

"""The :OnlyA limit test of an allosteric mechanism at one constant draw: the package model
(inactive state = non-:OnlyA steps on the forms the free form reaches) against the :OnlyA
family (every form and step in the inactive state, :OnlyA plain bindings' affinity divided by
e^(60·e_s), other :OnlyA steps' ratios free, steady-state ones slowed by 1e-24), pinned to the
package model on its own forms, at two draws of the family's remaining free directions. At a
positive point and on every face. `family = nothing` when no exponents e_s in 1:4 (6 tries)
make the family thermodynamically consistent."""
function onlya_test(ms, rng; tries = 8)
    A = aconsts(ms, rng)
    Ip = iconsts(ms, A, rng)
    xp = Ip.draw(rng)
    lnkp = Ip.lnk(xp)
    pin = Ip.pins(xp)
    If = nothing
    nbind = count(k -> !Ip.kept[k] && plain_binding(ms, ms.steps[k]), eachindex(ms.steps))
    for t in 1:tries
        expo = t == 1 ? nothing : rand(rng, 1:4, nbind)
        If = iconsts(ms, A, rng; family = true, expo, pin)
        If === nothing || break
    end
    L = big(lu10(rng, -1, 1))
    conc = Dict{Symbol, BigFloat}(x => big(lu10(rng, -1, 1)) for x in ms.mets)
    σ = [big(lu10(rng, -1, 1)) for _ in 1:(length(ms.steps) + 1)]
    pts = [("positive", Symbol[]); faces(ms)]
    out = NamedTuple[]
    if If === nothing
        for (lab, zs) in pts
            c = copy(conc); foreach(z -> c[z] = big(0), zs)
            push!(out, (point = lab, pkg = Float64(oracle_rate(ms, A, Ip, lnkp, L, c, big(10)^30, σ)),
                        fam = Float64[], spread = NaN, vs_pkg = NaN))
        end
        return (family = false, nfree = -1, points = out)
    end
    x1, x2 = If.draw(rng), If.draw(rng)
    l1, l2 = If.lnk(x1), If.lnk(x2)
    for (lab, zs) in pts
        c = copy(conc); foreach(z -> c[z] = big(0), zs)
        vp = oracle_rate(ms, A, Ip, lnkp, L, c, big(10)^30, σ)
        v1 = oracle_rate(ms, A, If, l1, L, c, big(10)^30, σ)
        v2 = oracle_rate(ms, A, If, l2, L, c, big(10)^30, σ)
        push!(out, (point = lab, pkg = Float64(vp), fam = Float64.([v1, v2]),
                    spread = Float64(rel(v1, v2)), vs_pkg = Float64(max(rel(v1, vp), rel(v2, vp)))))
    end
    (family = true, nfree = If.nfree, points = out)
end

# ─── the package's own constants as mass-action steps (plain mechanisms) ─────────────────
evx(x::Symbol, P) = P[x]
evx(x::Real, P) = big(x)
evx(e::Expr, P) = getfield(Base, e.args[1])((evx(a, P) for a in e.args[2:end])...)

"""Per-step lnK (natural log of kf/kr) and steady-state speeds of a plain package mechanism
for fitted values `θ` (Dict), as its law uses them, in the oracle's constant format."""
function pkg_aconsts(m::ER.Mechanism, ms, θ::Dict{Symbol, BigFloat}, keq::BigFloat)
    dep, _ = ER._dependent_param_exprs(m)
    rename = ER._build_wegscheider_rename_map(m)
    P = copy(θ); P[:Keq] = keq
    for (k, e) in dep
        P[k] = evx(e, P)
    end
    look(k) = haskey(P, k) ? P[k] : P[rename[k]]
    gs = ER.steps(m)
    lnK = Any[]; spd = Dict{Int, BigFloat}()
    for s in ms.steps
        rep = first(gs[s.g])
        if s.re
            K = look(ER.name(ER.Krapid(rep, :None), m))
            push!(lnK, ER.is_binding(s.step) ? -log(K) : log(K))
        else
            kf = look(ER.name(ER.Kfor(rep, :None), m)); kr = look(ER.name(ER.Krev(rep, :None), m))
            push!(lnK, log(kf / kr)); spd[s.g] = sqrt(kf * kr)
        end
    end
    (lnK = lnK, speed = [get(spd, g, big(1)) for g in 1:ms.ngroups])
end

"""Package law against the oracle with the package's own constants (plain mechanisms): at a
positive point and on each face, the law, the interior limit and the exact-face RE limits
under two speed patterns."""
function law_vs_oracle(m::ER.Mechanism, ms, pl, rng)
    θ = Dict{Symbol, BigFloat}(x => big(lu10(rng, -1, 1)) for x in pl.names)
    keq = big(lu10(rng, -1, 1))
    A = pkg_aconsts(m, ms, θ, keq)
    conc = Dict{Symbol, BigFloat}(x => big(lu10(rng, -1, 1)) for x in ms.mets)
    σs = [[big(lu10(rng, -3, 3)) for _ in 1:(length(ms.steps) + 1)] for _ in 1:2]
    law(c) = prate(pl, Float64[θ[x] for x in pl.names], Float64(keq),
                   Dict{Symbol, Float64}(x => Float64(v) for (x, v) in c))
    out = NamedTuple[]
    for (lab, zs) in [("positive", Symbol[]); faces(ms)]
        c0 = copy(conc); foreach(z -> c0[z] = big(0), zs)
        u = Dict(z => big(lu10(rng, -1, 1)) for z in zs)
        ce = copy(conc); foreach(z -> ce[z] = big(10)^-25 * u[z], zs)
        ex = [oracle_rate(ms, A, nothing, nothing, big(1), c0, big(10)^40, σ) for σ in σs]
        it = oracle_rate(ms, A, nothing, nothing, big(1), ce, big(10)^60, σs[1])
        vl = law(c0)
        push!(out, (point = lab, law = vl, interior = Float64(it), exact = Float64.(ex),
                    law_vs_interior = Float64(rel(big(vl), it)),
                    exact_spread = Float64(rel(ex[1], ex[2])),
                    law_vs_exact = Float64(max(rel(big(vl), ex[1]), rel(big(vl), ex[2])))))
    end
    out
end

# ─── long-time state from the free enzyme (initial condition E_A), for faces where some forms
#     become unreachable: the steady state the enzyme actually reaches ─────────────────────
"""Strongly connected component id of every node."""
function sccs(n, arcs)
    adj = [Int[] for _ in 1:n]; radj = [Int[] for _ in 1:n]
    for (a, b) in arcs
        push!(adj[a], b); push!(radj[b], a)
    end
    seen = falses(n); order = Int[]
    for s in 1:n
        seen[s] && continue
        stack = [(s, 1)]; seen[s] = true
        while !isempty(stack)
            u, i = stack[end]
            if i <= length(adj[u])
                stack[end] = (u, i + 1)
                w = adj[u][i]
                seen[w] || (seen[w] = true; push!(stack, (w, 1)))
            else
                pop!(stack); push!(order, u)
            end
        end
    end
    comp = zeros(Int, n); c = 0
    for s in reverse(order)
        comp[s] == 0 || continue
        c += 1; comp[s] = c; st = [s]
        while !isempty(st)
            u = pop!(st)
            for w in radj[u]
                comp[w] == 0 && (comp[w] = c; push!(st, w))
            end
        end
    end
    comp
end

"""The rate the enzyme reaches from the free form of the active state: the stationary state of
the terminal class reachable from it (absorption-weighted when several are), with every form
it cannot reach at zero. Also the number of terminal classes reachable from it and in all."""
function oracle_rate_free(ms, A, I, lnkI, L, conc, Λ, σ; slow = big(1e-24))
    nf = length(ms.forms)
    nI = I === nothing ? 0 : length(I.iforms)
    n = nf + nI
    M = zeros(BigFloat, n, n)
    arcs = Tuple{Int, Int}[]
    fl = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
    x1 = ms.subs[1]
    cc(xs) = prod((conc[x] for x in xs); init = big(1.0))
    function add!(a, b, rf, rr, w)
        M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
        rf > 0 && push!(arcs, (a, b)); rr > 0 && push!(arcs, (b, a))
        w != 0 && push!(fl, (a, b, rf, rr, w))
    end
    wt(s) = count(==(x1), s.cons) - count(==(x1), s.rel)
    for (k, s) in enumerate(ms.steps)
        base = s.re ? Λ * σ[k] : A.speed[s.g]
        h = exp(BigFloat(A.lnK[k]) / 2)
        add!(s.a, s.b, base * h * cc(s.cons), base / h * cc(s.rel), wt(s))
    end
    if I !== nothing
        for (k, s) in enumerate(ms.steps)
            lnkI[k] === nothing && continue
            base = s.re ? Λ * σ[k] : I.speed[s.g] * (I.kept[k] ? big(1) : slow)
            h = exp(BigFloat(lnkI[k]) / 2)
            add!(nf + I.pos[s.a], nf + I.pos[s.b], base * h * cc(s.cons), base / h * cc(s.rel), wt(s))
        end
        f = ms.frees[1]
        add!(f, nf + I.pos[f], Λ * σ[end] * sqrt(L), Λ * σ[end] / sqrt(L), 0)
    end
    comp = sccs(n, arcs)
    nc = maximum(comp)
    out = falses(nc)
    foreach(((a, b),) -> comp[a] != comp[b] && (out[comp[a]] = true), arcs)
    ntot = count(!, out)
    start = ms.frees[1]
    R = falses(n); R[start] = true; todo = [start]
    adj = [Int[] for _ in 1:n]
    foreach(((a, b),) -> push!(adj[a], b), arcs)
    while !isempty(todo)
        u = pop!(todo)
        for w in adj[u]
            R[w] || (R[w] = true; push!(todo, w))
        end
    end
    terms = unique([comp[i] for i in 1:n if R[i] && !out[comp[i]]])
    # stationary state of each reachable terminal class
    function stat(c)
        T = findall(==(c), comp)
        Q = M[T, T]
        Q[end, :] .= 1
        rhs = zeros(BigFloat, length(T)); rhs[end] = 1
        x = zeros(BigFloat, n)
        x[T] .= Q \ rhs
        x
    end
    x = if length(terms) == 1
        stat(terms[1])
    else
        # absorption probabilities from `start` into each class: h_c on the transient part of R
        tr = [i for i in 1:n if R[i] && !(comp[i] in terms)]
        tpos = Dict(t => j for (j, t) in enumerate(tr))
        xs = zeros(BigFloat, n)
        for c in terms
            # Σ_t rate(s,t) (h(t) - h(s)) = 0 on transients; h = 1 on c, 0 on other classes
            H = zeros(BigFloat, length(tr), length(tr)); b = zeros(BigFloat, length(tr))
            for (j, s) in enumerate(tr)
                for t in 1:n
                    t == s && continue
                    r = M[t, s]
                    r > 0 || continue
                    H[j, j] -= r
                    if haskey(tpos, t)
                        H[j, tpos[t]] += r
                    elseif comp[t] == c
                        b[j] -= r
                    end
                end
            end
            p = start in tr ? (H \ b)[tpos[start]] : (comp[start] == c ? big(1) : big(0))
            xs .+= p .* stat(c)
        end
        xs
    end
    v = sum((w * (rf * x[a] - rr * x[b]) for (a, b, rf, rr, w) in fl); init = big(0.0))
    (v = v, nreach = length(terms), ntotal = ntot)
end

"""Face test with the long-time state from the free enzyme: per face, the exact-face values for two
speed patterns, the interior limit, the terminal classes reachable / in all."""
function face_test_free(ms, rng; pattern = (-3, 3))
    A = aconsts(ms, rng)
    I = ms.allo ? iconsts(ms, A, rng) : nothing
    lnkI = I === nothing ? nothing : I.lnk(I.draw(rng))
    L = big(lu10(rng, -1, 1))
    conc = Dict{Symbol, BigFloat}(x => big(lu10(rng, -1, 1)) for x in ms.mets)
    nE = length(ms.steps) + 1
    σs = [[big(lu10(rng, pattern...)) for _ in 1:nE] for _ in 1:2]
    out = NamedTuple[]
    for (lab, zs) in faces(ms)
        c0 = copy(conc); foreach(z -> c0[z] = big(0), zs)
        u = Dict(z => big(lu10(rng, -1, 1)) for z in zs)
        ce = copy(conc); foreach(z -> ce[z] = big(10)^-25 * u[z], zs)
        ex = [oracle_rate_free(ms, A, I, lnkI, L, c0, big(10)^40, σ) for σ in σs]
        it = oracle_rate_free(ms, A, I, lnkI, L, ce, big(10)^60, σs[1])
        pat = rel(ex[1].v, ex[2].v)
        gap = max(rel(ex[1].v, it.v), rel(ex[2].v, it.v))
        push!(out, (face = lab, exact = Float64.([e.v for e in ex]), interior = Float64(it.v),
                    pattern = Float64(pat), gap = Float64(gap), nreach = ex[1].nreach,
                    ntotal = ex[1].ntotal, fails = pat > 1e-8 || gap > 1e-8))
    end
    out
end
