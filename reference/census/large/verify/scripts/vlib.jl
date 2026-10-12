# ABOUTME: Verification library: a BigFloat mass-action steady-state oracle for plain and n=1
# ABOUTME: allosteric slice mechanisms (pruned or with :OnlyA steps vanishing), law values, ranks.

const CEN = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad/census/large"
const VER = joinpath(CEN, "verify")
include(joinpath(CEN, "scripts", "census_lib.jl"))   # loads EnzymeRates (ER), Slice (S), helpers
include(joinpath(CEN, "scripts", "summarize.jl"))    # label! (the census's own labels)
using LinearAlgebra, Random, Printf, Serialization

setprecision(BigFloat, 512)

"""The saved census run of reaction `name`, labelled as the census labels it."""
function load_run(name)
    d = deserialize(joinpath(CEN, "data", "$(name).jls"))
    label!(d.rounds)
    d
end

fmtb(x) = isnan(x) ? "NaN" : @sprintf("%.6g", Float64(x))
lu10(rng, lo, hi) = big(10.0)^(lo + (hi - lo) * rand(rng))

# ─── the slice law's value, assembled from its parts ────────────────────────────────────

"""The values of `L.vars`: concentrations from `conc`, every other name from `vals`."""
lawx(L, vals, conc) = [L.vars.conc[i] ? big(conc[n]) : big(vals[n])
                       for (i, n) in enumerate(L.vars.names)]

"""`v/E_total` of the slice law `L` (a `Law`, or an n-subunit `AlloLaw`) at `x`."""
function law_value(L, x::AbstractVector{T}) where {T <: Real}
    ev(p) = S.evaluate(p, x)
    L isa S.Law && return ev(L.H) * ev(L.Q) / ev(L.D)
    QA, DA = ev(L.active.Q), ev(L.active.D)
    QI, DI = L.inactive === nothing ? (zero(T), one(T)) :
             (ev(L.inactive.Q), ev(L.inactive.D))
    fA, fI, Lv, n = ev(L.fA), ev(L.fI), x[L.vars.index[:L]], L.n
    WA = prod((ev(s[1])^s[3] for s in L.sites); init = one(T))
    WI = prod((ev(s[2])^s[3] for s in L.sites); init = one(T))
    num = fI^n * QA * DA^(n - 1) * WA + Lv * fA^n * QI * DI^(n - 1) * WI
    den = fI^n * DA^n * WA + Lv * fA^n * DI^n * WI
    ev(L.H) * num / den
end

# ─── every per-step constant, the law's own and the ones outside it ─────────────────────

"""
The log value of every column of `m`'s constraint system (`Slice.constraint_system`): the
law's columns read by name from `vals` (`all_parameters`), the other columns (inactive
constants outside the law) solved from the cycle rows by exact elimination, each free one
drawn uniformly in `±spread` decades. Returns `(cols, x, resid, nfree)`; `resid` is the
largest violation of a cycle row by the law's own values (0 when they are consistent).
"""
function all_columns(m, vals, rng; spread = 1.0)
    cols, A, b = S.constraint_system(m)
    nA = S.constant_names(m)
    nI = S.is_allosteric(m) ? S.constant_names(m; inactive = true) : nA
    name(c) = (c.inactive ? nI : nA)[c.edge][c.slot]
    n = length(cols)
    x = zeros(BigFloat, n)
    known = falses(n)
    for (i, c) in enumerate(cols)
        nm = name(c)
        haskey(vals, nm) || continue
        x[i] = log(big(vals[nm]))
        known[i] = true
    end
    U = findall(.!known)
    lk = log(big(vals[:Keq]))
    r = [b[i] * lk - sum((A[i, j] * x[j] for j in findall(known)); init = big(0.0))
         for i in axes(A, 1)]
    M = Rational{BigInt}.(A[:, U])
    nr, nu = size(M)
    piv, row = Int[], 1
    for col in 1:nu
        row > nr && break
        p = findfirst(i -> !iszero(M[i, col]), row:nr)
        p === nothing && continue
        p += row - 1
        M[[row, p], :] = M[[p, row], :]
        r[row], r[p] = r[p], r[row]
        f = inv(M[row, col])
        M[row, :] .*= f
        r[row] *= BigFloat(f)
        for i in 1:nr
            (i == row || iszero(M[i, col])) && continue
            g = M[i, col]
            M[i, :] .-= g .* M[row, :]
            r[i] -= BigFloat(g) * r[row]
        end
        push!(piv, col)
        row += 1
    end
    resid = maximum(abs, r[row:end]; init = big(0.0))
    free = setdiff(1:nu, piv)
    xu = zeros(BigFloat, nu)
    foreach(f -> xu[f] = spread * (2 * rand(rng) - 1) * log(big(10.0)), free)
    for (k, pc) in enumerate(piv)
        xu[pc] = r[k] - sum((BigFloat(M[k, f]) * xu[f] for f in free); init = big(0.0))
    end
    x[U] .= xu
    cols, x, resid, length(free)
end

# ─── :OnlyA components, arcs and potentials (the island rule's objects, rebuilt) ────────

"""Components of the forms joined by steps outside `:OnlyA` groups, the `:OnlyA` plain
bindings as arcs `(unbound component, bound component, edge)`, the flip form's component,
and the components no chain of arcs from it reaches."""
function onlya_structure(m)
    nf = length(m.forms)
    parent = collect(1:nf)
    top(i) = (while parent[i] != i; i = parent[i]; end; i)
    for (k, e) in enumerate(m.edges)
        m.tags[m.group[k]] === :OnlyA && continue
        parent[top(e.a)] = top(e.b)
    end
    comp = [top(i) for i in 1:nf]
    fwd = S.orient(m)
    arcs = Tuple{Int, Int, Int}[]
    for (k, e) in enumerate(m.edges)
        m.tags[m.group[k]] === :OnlyA || continue
        o = S.oriented(m, k, fwd)
        S.is_plain_binding(m, o) && push!(arcs, (comp[o.a], comp[o.b], k))
    end
    c0 = comp[1]
    seen = Set([c0])
    todo = [c0]
    while !isempty(todo)
        u = pop!(todo)
        for (a, b, _) in arcs
            a == u && !(b in seen) && (push!(seen, b); push!(todo, b))
        end
    end
    comps = unique(comp)
    (comp = comp, arcs = arcs, c0 = c0, lost = [c for c in comps if !(c in seen)],
     comps = comps)
end

"""A potential on the components: `c0` at 0, every other source of the arc graph at
`offset`, each component above its arc predecessors by 1 to 2 (longest-path layering)."""
function potentials(st, offset, rng)
    preds = Dict(c => [a for (a, b, _) in st.arcs if b == c] for c in st.comps)
    phi = Dict{Int, BigFloat}()
    remaining = copy(st.comps)
    while !isempty(remaining)
        i = findfirst(c -> all(a -> haskey(phi, a), preds[c]), remaining)
        i === nothing && error("the :OnlyA arcs close a cycle")
        c = remaining[i]
        ps = preds[c]
        phi[c] = isempty(ps) ? (c == st.c0 ? big(0.0) : big(offset)) :
                 maximum(phi[a] for a in ps) + 1 + big(rand(rng))
        deleteat!(remaining, i)
    end
    phi
end

# ─── the mass-action network and its steady state ───────────────────────────────────────

"""
The directed transitions `(from, to, rate, w)` of `m`'s mass-action network at
concentrations `conc`: every edge forward and back, rapid-equilibrium edges at speed
`Λ·speed[k]` (ratio kept), `w` the copies of `x1` the transition takes up. For an allosteric
`m` the n=1 two-conformation network, its forms `nf + i` the inactive copies, joined at the
flip form by a rapid-equilibrium step of ratio `L`. With `mode = :pruned` the inactive
conformation holds only the edges of `inactive_mask(m)`; with `mode = :vanish` it holds every
edge, the `:OnlyA` ones in their thermodynamic inactive constants, the plain bindings' rate
of taking up multiplied by `exp(t·Δ)` and every other `:OnlyA` edge's two rates by
`exp(-t·σ ± t·Δ/2)`, `Δ` the potential drop `phi[a] - phi[b]` between its components.
"""
function transitions(m, cols, x, vals, conc; Λ, speed, mode = :pruned, phi = nothing,
                     t = big(0.0), sigma = nothing, x1 = m.substrates[1])
    fwd = S.orient(m)
    nf = length(m.forms)
    col = Dict(c => i for (i, c) in enumerate(cols))
    allo = S.is_allosteric(m)
    imask = allo ? S.inactive_mask(m) : falses(length(m.edges))
    st = mode === :vanish ? onlya_structure(m) : nothing
    out = Tuple{Int, Int, BigFloat, Int}[]
    cc(p) = prod((big(conc[y.name]) for y in S.mets_of(m, p)); init = big(1.0))
    nx(p) = count(y -> y.name === x1 && !y.modifier, S.mets_of(m, p))
    for inactive in (allo ? (false, true) : (false,)), (k, e) in enumerate(m.edges)
        tag = allo ? m.tags[m.group[k]] : :EqualAI
        if inactive
            mode === :pruned && !imask[k] && continue
        end
        o = S.oriented(m, k, fwd)
        ci(s) = col[S.Col(k, s, inactive && tag !== :EqualAI)]
        if o.re
            K = exp(x[ci(1)])
            sp = Λ * speed[k]
            kf, kr = S.is_binding(o) ? (sp, sp * K) : (sp * K, sp)
        else
            kf, kr = exp(x[ci(1)]), exp(x[ci(2)])
        end
        if inactive && tag === :OnlyA && mode === :vanish
            d = phi[st.comp[o.a]] - phi[st.comp[o.b]]
            if S.is_plain_binding(m, o)
                kf *= exp(t * d)
            else
                kf *= exp(-t * sigma[k] + t * d / 2)
                kr *= exp(-t * sigma[k] - t * d / 2)
            end
        end
        off = inactive ? nf : 0
        w = nx(o.cons) - nx(o.rel)
        rf, rr = kf * cc(o.cons), kr * cc(o.rel)
        rf > 0 && push!(out, (off + o.a, off + o.b, rf, w))
        rr > 0 && push!(out, (off + o.b, off + o.a, rr, -w))
    end
    if allo
        sp = Λ * speed[end]
        push!(out, (1, nf + 1, sp * big(vals[:L]), 0), (nf + 1, 1, sp, 0))
    end
    out
end

"""The closed classes (strongly connected components with no exit) of the graph of `tr` on
the nodes `nodes` (Kosaraju)."""
function closed_classes(nodes, tr)
    at = Dict(u => i for (i, u) in enumerate(nodes))
    n = length(nodes)
    adj, radj = [Int[] for _ in 1:n], [Int[] for _ in 1:n]
    for (a, b, _, _) in tr
        push!(adj[at[a]], at[b])
        push!(radj[at[b]], at[a])
    end
    seen, order = falses(n), Int[]
    for s in 1:n
        seen[s] && continue
        stack = [(s, 1)]
        seen[s] = true
        while !isempty(stack)
            u, i = stack[end]
            if i <= length(adj[u])
                stack[end] = (u, i + 1)
                v = adj[u][i]
                seen[v] || (seen[v] = true; push!(stack, (v, 1)))
            else
                pop!(stack)
                push!(order, u)
            end
        end
    end
    comp, c = zeros(Int, n), 0
    for s in reverse(order)
        comp[s] == 0 || continue
        c += 1
        comp[s] = c
        todo = [s]
        while !isempty(todo)
            u = pop!(todo)
            for v in radj[u]
                comp[v] == 0 && (comp[v] = c; push!(todo, v))
            end
        end
    end
    exits = falses(c)
    for u in 1:n, v in adj[u]
        comp[u] != comp[v] && (exits[comp[u]] = true)
    end
    [[nodes[u] for u in 1:n if comp[u] == k] for k in 1:c if !exits[k]]
end

"""
The steady-state net consumption of `x1` per enzyme of the transitions `tr`, the enzyme
on the nodes they touch: `(v, nclasses)`, `v = NaN` when the graph has more than one closed
class (no unique steady state). Nodes outside the closed class hold no enzyme.
"""
function steady_rate(tr)
    nodes = sort!(unique!([[a for (a, _, _, _) in tr]; [b for (_, b, _, _) in tr]]))
    classes = closed_classes(nodes, tr)
    length(classes) == 1 || return big(NaN), length(classes)
    C = classes[1]
    at = Dict(u => i for (i, u) in enumerate(C))
    n = length(C)
    G = zeros(BigFloat, n, n)
    for (a, b, r, _) in tr
        haskey(at, a) && haskey(at, b) || continue
        G[at[b], at[a]] += r
        G[at[a], at[a]] -= r
    end
    G[n, :] .= 1
    rhs = zeros(BigFloat, n)
    rhs[n] = 1
    p = G \ rhs
    v = sum((w * r * p[at[a]] for (a, b, r, w) in tr if w != 0 && haskey(at, a));
            init = big(0.0))
    v, 1
end

"""Random fitted values for `m` (`10^[lo, hi]`), `Keq` and `E_total = 1`, with every
constant of the law (`all_parameters`) as a name → BigFloat dictionary."""
function random_vals(m, rng; lo = -1.0, hi = 1.0)
    names = S.fitted_params(m)
    p = (; (x => lu10(rng, lo, hi) for x in names)..., Keq = lu10(rng, -1, 1),
         E_total = big(1.0))
    vals_from(m, p)
end
vals_from(m, p) = Dict{Symbol, BigFloat}(k => big(v) for (k, v) in pairs(S.all_parameters(m, p)))

"""Random positive concentrations of `metabolites(m)` in `10^[lo, hi]`."""
random_conc(m, rng; lo = -1.0, hi = 1.0) =
    Dict{Symbol, BigFloat}(x => lu10(rng, lo, hi) for x in S.metabolites(m))

"""Random RE speed factors, one per edge and one for the flip, in `10^[-dec, dec]`."""
random_speed(m, rng; dec = 1.0) = [lu10(rng, -dec, dec) for _ in 1:(length(m.edges) + 1)]

"""The oracle rate of `m` (`transitions` then `steady_rate`)."""
oracle(m, cols, x, vals, conc; kw...) = steady_rate(transitions(m, cols, x, vals, conc; kw...))

"""The faces the census checks: every product zero, every substrate zero, each reactant."""
vfaces(m) = [("products=0", m.products); ("substrates=0", m.substrates);
             [(string(y, "=0"), [y]) for y in [m.substrates; m.products]]]

# ─── rank by pivoted Gram–Schmidt in BigFloat ───────────────────────────────────────────

"""Column-pivoted modified Gram–Schmidt of `J` (BigFloat, twice orthogonalized): the
residual norms of the pivots in order, relative to the largest column norm."""
function pivoted_norms(J::Matrix{BigFloat})
    Q = copy(J)
    nc = size(Q, 2)
    s1 = maximum(norm(view(Q, :, j)) for j in 1:nc)
    left = collect(1:nc)
    out = BigFloat[]
    basis = Vector{Vector{BigFloat}}()
    while !isempty(left)
        for j in left, _ in 1:2, q in basis
            Q[:, j] .-= dot(q, view(Q, :, j)) .* q
        end
        nrm = [norm(view(Q, :, j)) for j in left]
        i = argmax(nrm)
        push!(out, nrm[i] / s1)
        j = left[i]
        nrm[i] > 0 && push!(basis, Q[:, j] ./ nrm[i])
        deleteat!(left, i)
    end
    out
end
