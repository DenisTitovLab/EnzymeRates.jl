# ABOUTME: Skeptic's independent route: its own constraint system, consistent random constants,
# ABOUTME: steady-state and MWC network oracles, and checks of a slice law against them.
const SCR = "/tmp/claude-501/-home-denis-linux--julia-dev-EnzymeRates/12f42995-f5ba-4fda-9510-dd9b87f7d069/scratchpad"
isdefined(Main, :Oracle) || include(joinpath(SCR, "proto", "common", "Oracle.jl"))
isdefined(Main, :Slice) || include(joinpath(SCR, "slice", "src", "Slice.jl"))
isdefined(Main, :retag) || include(joinpath(SCR, "slice", "core", "retag.jl"))
using Serialization, Random, Printf, LinearAlgebra, Statistics
const S = Slice
setprecision(BigFloat, 512)

# ─── records ────────────────────────────────────────────────────────────────────────────

const CORE = joinpath(SCR, "slice", "core", "data")
const IDATA_IN = joinpath(SCR, "slice", "integrate", "data")

"""Records `(set, p, c)` of a converted core set, or of a ter-ter export (converted here)."""
function skeptic_records(set)
    if isfile(joinpath(CORE, set * ".jls"))
        row = only(r for r in (split(l, '\t') for l in readlines(joinpath(CORE, "MANIFEST.tsv")))
                   if r[1] == set)
        pops = deserialize(row[4])
        return [(set = set, p = set == "p5a_cases" ? retag(pops[c.index]) : pops[c.index], c = c)
                for c in deserialize(joinpath(CORE, set * ".jls"))]
    end
    if set == "init_terter"           # 8 random records of each part, converted alone
        out = []
        for k in 1:11
            ps = deserialize(joinpath(IDATA_IN, "pop_init_terter.part$k.jls"))
            for i in randperm(Xoshiro(k), length(ps))[1:8]
                push!(out, (set = "init_terter.part$k", p = ps[i], c = S.Converted(ps[i], i)))
            end
        end
        return out
    end
    [(set = set, p = p, c = S.Converted(p, i))
     for (i, p) in enumerate(deserialize(joinpath(IDATA_IN, "pop_$set.jls")))]
end

# ─── the independent constraint system (Pop steps, stored direction) ───────────────────

"""Signed metabolite change of Pop step `s` going from → to, by name (released − consumed)."""
function step_change(s)
    d = Dict{Symbol, Int}()
    foreach(x -> d[x.name] = get(d, x.name, 0) - 1, s.consumed)
    foreach(x -> d[x.name] = get(d, x.name, 0) + 1, s.released)
    d
end

"""Fundamental cycles of the Pop graph from a BFS tree rooted at the last species, scanning
steps in reverse order: `(z, multiple)` with `z[j]` the signed traversals of step `j` and
`multiple` the net reactions the cycle runs (error when the change is no multiple)."""
function my_cycles(p)
    n, J = length(p.species), length(p.steps)
    path = Vector{Union{Nothing, Vector{Int}}}(nothing, n)
    path[n] = zeros(Int, J)
    tree, grown = falses(J), true
    while grown
        grown = false
        for j in reverse(1:J)
            s = p.steps[j]
            a, b = path[s.from], path[s.to]
            (a === nothing) == (b === nothing) && continue
            if a !== nothing
                path[s.to] = copy(a); path[s.to][j] += 1
            else
                path[s.from] = copy(b); path[s.from][j] -= 1
            end
            tree[j] = grown = true
        end
    end
    any(isnothing, path) && error("graph not connected")
    out = Tuple{Vector{Int}, Int}[]
    for j in 1:J
        tree[j] && continue
        s = p.steps[j]
        z = path[s.from] - path[s.to]
        z[j] += 1
        tot = Dict{Symbol, Int}()
        for (i, zi) in enumerate(z), (x, d) in step_change(p.steps[i])
            tot[x] = get(tot, x, 0) + zi * d
        end
        c = -get(tot, p.substrates[1], 0)
        for (x, v) in tot
            want = x in p.substrates ? -c : x in p.products ? c : 0
            v == want || error("cycle closed by step $j changes $x by $v, not $want")
        end
        push!(out, (z, c))
    end
    out
end

"""Alignment of group member `j` with lead `j0` in stored direction (decision 2, forward
with forward): by the metabolites taken up and given off; for steps without metabolites, by
the conformational change; `:ambiguous` for an isomerization without one."""
function alignment(p, j0, j)
    mets(s, rev) = (q = (sort([(x.name, x.role) for x in s.consumed]),
                         sort([(x.name, x.role) for x in s.released]));
                    rev ? reverse(q) : q)
    cf(i) = p.species[i].conformation
    turn(s, rev) = rev ? (cf(s.to), cf(s.from)) : (cf(s.from), cf(s.to))
    turn2(s, rev) = (q = turn(s, rev); q[1] == q[2] ? nothing : q)
    l, s = p.steps[j0], p.steps[j]
    same, flip = mets(s, false) == mets(l, false), mets(s, true) == mets(l, false)
    same && !flip && return :same
    flip && !same && return :flip
    same || error("group members differ")
    # no free metabolite: forward chemistry with forward chemistry (the bound content change)
    bnd(i) = [(x.name, x.role) for x in p.species[i].bound]
    minus(a, b) = (r = copy(a); foreach(x -> (k = findfirst(==(x), r); k === nothing ||
                                              deleteat!(r, k)), b); sort(r))
    delta(f, t) = (minus(bnd(f), bnd(t)), minus(bnd(t), bnd(f)))   # (lost, gained)
    chem(s, rev) = rev ? delta(s.to, s.from) : delta(s.from, s.to)
    cl = chem(l, false)
    if cl[1] != cl[2]
        chem(s, false) == cl && return :same
        chem(s, true) == cl && return :flip
    end
    t = turn2(l, false)
    t === nothing && turn2(s, false) === nothing && return :ambiguous
    turn2(s, false) == t ? :same : turn2(s, true) == t ? :flip : error("group members differ")
end

"""
    my_system(p, tags) -> (cols, A, b, inlaw)

The joint log-linear system over per-step constants (`(j, slot, state)`, oracle convention,
stored direction): cycle rows per state (whole graph), ties per group slot by slot after
alignment (inactive: `:NonequalAI` members tied, `:OnlyA` untied, `:EqualAI` shares the active
columns), right side in multiples of log Keq. `inlaw[i]` marks columns of the law: every
active column, and inactive columns of steps the flip species reaches through non-`:OnlyA`
steps.
"""
function my_system(p, tags)
    allo = !isempty(tags)
    g = Pop.group_of_step(p)
    tag(j) = allo ? tags[g[j]] : :plain
    cols, idx = Tuple{Int, Int, Symbol}[], Dict{Tuple{Int, Int, Symbol}, Int}()
    for st in (allo ? (:A, :I) : (:A,)), (j, s) in enumerate(p.steps), r in 1:(s.re ? 1 : 2)
        if st === :I && tag(j) === :EqualAI
            idx[(j, r, :I)] = idx[(j, r, :A)]
            continue
        end
        push!(cols, (j, r, st)); idx[(j, r, st)] = length(cols)
    end
    rows, b = Vector{Rational{BigInt}}[], Rational{BigInt}[]
    function lnka!(row, j, st, w)              # w·log(forward equilibrium ratio of step j)
        s = p.steps[j]
        if s.re
            row[idx[(j, 1, st)]] += Pop.is_binding(s) ? -w : w
        else
            row[idx[(j, 1, st)]] += w; row[idx[(j, 2, st)]] -= w
        end
    end
    for st in (allo ? (:A, :I) : (:A,)), (z, c) in my_cycles(p)
        row = zeros(Rational{BigInt}, length(cols))
        foreach(j -> z[j] == 0 || lnka!(row, j, st, z[j]), eachindex(z))
        push!(rows, row); push!(b, c)
    end
    for grp in p.groups, j in grp[2:end]
        j0 = grp[1]
        al = alignment(p, j0, j)
        al === :ambiguous && (al = :same)
        for st in (allo ? (:A, :I) : (:A,))
            st === :I && tag(j) in (:OnlyA, :EqualAI) && continue
            s = p.steps[j]
            for r in 1:(s.re ? 1 : 2)
                row = zeros(Rational{BigInt}, length(cols))
                if s.re
                    row[idx[(j, 1, st)]] += 1
                    row[idx[(j0, 1, st)]] += al === :same || Pop.is_binding(s) ? -1 : 1
                else
                    row[idx[(j, r, st)]] += 1
                    row[idx[(j0, al === :same ? r : 3 - r, st)]] -= 1
                end
                push!(rows, row); push!(b, 0)
            end
        end
    end
    inlaw = trues(length(cols))
    if allo
        flips = [i for (i, sp) in enumerate(p.species) if Pop.is_free(sp)]
        flip = flips[argmin([String(p.species[i].conformation) for i in flips])]
        seen = falses(length(p.species)); seen[flip] = true
        grown = true
        while grown
            grown = false
            for (j, s) in enumerate(p.steps)
                tag(j) === :OnlyA && continue
                seen[s.from] == seen[s.to] && continue
                seen[s.from] = seen[s.to] = grown = true
            end
        end
        keep = [tag(j) !== :OnlyA && seen[s.from] for (j, s) in enumerate(p.steps)]
        inlaw = [st === :A || keep[j] for (j, _, st) in cols]
    end
    A = isempty(rows) ? zeros(Rational{BigInt}, 0, length(cols)) : permutedims(reduce(hcat, rows))
    (; cols, idx, A, b, inlaw, allo)
end

"""Exact RREF of `M` in place over its first `n` columns; the pivot columns."""
function my_rref!(M, n)
    piv, r = Int[], 1
    for col in 1:n
        r > size(M, 1) && break
        q = findfirst(i -> !iszero(M[i, col]), r:size(M, 1))
        q === nothing && continue
        q += r - 1
        M[[r, q], :] = M[[q, r], :]
        M[r, :] ./= M[r, col]
        for i in axes(M, 1)
            i != r && !iszero(M[i, col]) && (M[i, :] .-= M[i, col] .* M[r, :])
        end
        push!(piv, col); r += 1
    end
    piv
end

"""
    my_solution(sys, rng) -> (order, M, piv, nfree_law)

RREF of `[A b]` with the columns outside the law first, each block in a random order; the
law's free count is the number of free columns among the law's columns.
"""
function my_solution(sys, rng)
    n = length(sys.cols)
    out_law, in_law = shuffle(rng, findall(.!sys.inlaw)), shuffle(rng, findall(sys.inlaw))
    order = [out_law; in_law]
    M = [sys.A[:, order] sys.b]
    piv = my_rref!(M, n)
    all(iszero, M[length(piv) + 1:end, n + 1]) || error("inconsistent system")
    nfree_law = count(i -> !(i in piv), length(out_law) + 1:n)
    (; order, M, piv, nfree_law)
end

"""Random log-constants on the solution set: free columns uniform in ±2 decades, `lkeq` log
Keq; returns per-column BigFloat values (by original column index)."""
function my_sample(sys, sol, rng, lkeq)
    n = length(sys.cols)
    x = zeros(BigFloat, n)
    free = setdiff(1:n, sol.piv)
    for f in free
        x[f] = (4 * big(rand(rng)) - 2) * log(big(10))
    end
    for (r, pc) in enumerate(sol.piv)
        x[pc] = BigFloat(sol.M[r, n + 1]) * lkeq -
                sum((BigFloat(sol.M[r, f]) * x[f] for f in free if !iszero(sol.M[r, f]));
                    init = big(0))
    end
    vals = zeros(BigFloat, n)
    for (pos, i) in enumerate(sol.order)
        vals[i] = exp(x[pos])
    end
    vals
end

"""Per Pop step constants of state `st` (oracle convention) from column values."""
step_values(p, sys, vals, st) =
    [[vals[sys.idx[(j, r, st)]] for r in 1:(s.re ? 1 : 2)] for (j, s) in enumerate(p.steps)]

# ─── oracles ────────────────────────────────────────────────────────────────────────────

"""Steady state of Pop record `p` over the steps `use` and the species `forms` with
constants `k` (oracle convention): the populations and the net consumption of the first
substrate, RE steps as fast pairs at `Λ`."""
function my_network(p, k, use, forms, concs; Λ = big(10)^60)
    loc = Dict(f => i for (i, f) in enumerate(forms))
    n = length(forms)
    n == 1 && return [big(1)], big(0), loc
    M = zeros(BigFloat, n, n)
    x1 = p.substrates[1]
    fl = Tuple{Int, Int, BigFloat, BigFloat, Int}[]
    for (j, s) in enumerate(p.steps)
        use[j] || continue
        a, b = loc[s.from], loc[s.to]
        f, r = s.re ? (Pop.is_binding(s) ? (Λ, Λ * k[j][1]) : (Λ, Λ / k[j][1])) :
               (k[j][1], k[j][2])
        rf = f * prod((BigFloat(concs[x.name]) for x in s.consumed); init = big(1))
        rr = r * prod((BigFloat(concs[x.name]) for x in s.released); init = big(1))
        M[b, a] += rf; M[a, a] -= rf; M[a, b] += rr; M[b, b] -= rr
        w = count(x -> x.name === x1, s.consumed) - count(x -> x.name === x1, s.released)
        push!(fl, (a, b, rf, rr, w))
    end
    M[n, :] .= 1
    rhs = zeros(BigFloat, n); rhs[n] = 1
    x = M \ rhs
    x, sum((w * (rf * x[a] - rr * x[b]) for (a, b, rf, rr, w) in fl if w != 0); init = big(0)),
    loc
end

"""The MWC (formulation 1) rate per site: each conformation's own network (inactive: the
law's steps only, deletion), flip at the free species `flip`, regulator factors `WA`, `WI`."""
function my_mwc(p, kA, kI, keepI, flip, L, WA, WI, n, concs; Λ = big(10)^60)
    xA, JA, locA = my_network(p, kA, trues(length(p.steps)), collect(eachindex(p.species)),
                              concs; Λ)
    formsI = sort!(unique([flip; [f for (j, s) in enumerate(p.steps) if keepI[j]
                                  for f in (s.from, s.to)]]))
    xI, JI, locI = my_network(p, kI, keepI, formsI, concs; Λ)
    wA = WA / xA[locA[flip]]^n
    wI = L * WI / xI[locI[flip]]^n
    (wA * JA + wI * JI) / (wA + wI)
end

"""The flip species and the inactive law's steps of `p` (own computation)."""
function my_inactive(p, tags)
    g = Pop.group_of_step(p)
    flips = [i for (i, sp) in enumerate(p.species) if Pop.is_free(sp)]
    flip = flips[argmin([String(p.species[i].conformation) for i in flips])]
    seen = falses(length(p.species)); seen[flip] = true
    grown = true
    while grown
        grown = false
        for (j, s) in enumerate(p.steps)
            tags[g[j]] === :OnlyA && continue
            seen[s.from] == seen[s.to] && continue
            seen[s.from] = seen[s.to] = grown = true
        end
    end
    flip, [tags[g[j]] !== :OnlyA && seen[s.from] for (j, s) in enumerate(p.steps)]
end

"""Regulator factor of one state: Π_sites (1 + Σ [l]/K)^mult over the ligands it binds."""
function my_W(m, regvals, concs, st)
    w = big(1)
    for s in m.sites
        z = big(1)
        for (l, state) in zip(s.ligands, s.states)
            (st === :A && state === :OnlyI || st === :I && state === :OnlyA) && continue
            name = state === :EqualAI ? Symbol("K_", l, "reg") :
                   Symbol("K_", st === :A ? "A" : "I", "_", l, "reg")
            z += BigFloat(concs[l]) / regvals[name]
        end
        w *= z^s.multiplicity
    end
    w
end

# ─── the slice's parameters from a sample ───────────────────────────────────────────────

"""The value of slice column `col` (edge, slot, state) from per-step values `kA`/`kI`."""
function column_value(c, fwd, col, kA, kI)
    j = findfirst(==(col.edge), c.edge_of)
    k = col.inactive ? kI : kA
    e = c.mech.edges[col.edge]
    same = fwd[col.edge] != c.flipped[j]
    same && return k[j][col.slot]
    e.re ? (S.is_binding(e) ? k[j][1] : 1 / k[j][1]) : k[j][3 - col.slot]
end

"""`parameters(m)` as a NamedTuple of BigFloat from a sample (`kA`, `kI` per step), with
random regulator constants and `L`."""
function slice_params(c, sol, kA, kI, keq, rng)
    m = c.mech
    fwd = S.orient(m)
    names = S.column_names(m, sol)
    d = Pair{Symbol, BigFloat}[names[i] => column_value(c, fwd, sol.cols[i], kA, kI)
                               for i in sol.indep]
    if S.is_allosteric(m)
        for (_, ns) in S.regulator_names(m), x in ns
            push!(d, x => Oracle.logu(rng))
        end
        push!(d, :L => Oracle.logu(rng))
    end
    push!(d, :Keq => keq, :E_total => big(1))
    (; d...)
end

relerr(a, b) = Float64(abs(a - b) / max(abs(b), big(1e-300)))

# ─── kcat by the oracle ─────────────────────────────────────────────────────────────────

"""Biconnected blocks of the Pop graph (own Tarjan): block id per step."""
function my_blocks(p)
    n, J = length(p.species), length(p.steps)
    adj = [Tuple{Int, Int}[] for _ in 1:n]
    for (j, s) in enumerate(p.steps)
        push!(adj[s.from], (s.to, j)); push!(adj[s.to], (s.from, j))
    end
    disc, low, t = zeros(Int, n), zeros(Int, n), Ref(0)
    block, stack, nb = zeros(Int, J), Int[], Ref(0)
    function dfs(u, pe)
        t[] += 1; disc[u] = low[u] = t[]
        for (w, j) in adj[u]
            j == pe && continue
            if disc[w] == 0
                push!(stack, j)
                dfs(w, j)
                low[u] = min(low[u], low[w])
                if low[w] >= disc[u]
                    nb[] += 1
                    while true
                        x = pop!(stack); block[x] = nb[]; x == j && break
                    end
                end
            elseif disc[w] < disc[u]
                push!(stack, j); low[u] = min(low[u], disc[w])
            end
        end
    end
    for v in 1:n
        disc[v] == 0 && dfs(v, 0)
    end
    block
end

"""Steps outside the dead-end branches (own computation): blocks without a turnover cycle
are dropped while they meet the remaining steps at one species at most."""
function my_core(p)
    block = my_blocks(p)
    J = length(p.steps)
    x1 = p.substrates[1]
    w(s) = count(x -> x.name === x1, s.consumed) - count(x -> x.name === x1, s.released)
    turnover = Set{Int}()
    for b in unique(block)
        js = findall(==(b), block)
        pot = Dict{Int, Int}(p.steps[js[1]].from => 0)
        grown = true
        while grown
            grown = false
            for j in js
                s = p.steps[j]
                if haskey(pot, s.from) && !haskey(pot, s.to)
                    pot[s.to] = pot[s.from] - w(s); grown = true
                elseif haskey(pot, s.to) && !haskey(pot, s.from)
                    pot[s.from] = pot[s.to] + w(s); grown = true
                end
            end
        end
        any(j -> pot[p.steps[j].to] != pot[p.steps[j].from] - w(p.steps[j]), js) &&
            push!(turnover, b)
    end
    keep = trues(J)
    while true
        removed = false
        for b in unique(block[keep])
            b in turnover && continue
            js = findall(j -> keep[j] && block[j] == b, 1:J)
            vb = Set(f for j in js for f in (p.steps[j].from, p.steps[j].to))
            vrest = Set(f for j in 1:J if keep[j] && block[j] != b
                        for f in (p.steps[j].from, p.steps[j].to))
            length(intersect(vb, vrest)) <= 1 || continue
            keep[js] .= false; removed = true
        end
        removed || break
    end
    keep
end

"""The oracle's kcat candidates: per corner of the non-reactant names (each at `Y` or "zero"),
`v/E_total` with products at "zero" and substrates at `X`·u (u random in [1, 10]); `f(concs)`
evaluates a state. "Zero" is 1e-80: the rapid-equilibrium law is the limit of fast steps at a
positive concentration, and at exactly 0 a finite-rate oracle can differ (the specs/40 note);
the evaluation runs at 1024 bits for the wide range of rates."""
function oracle_corners(m, rng, f; X = big(10)^24, Y = big(10)^12)
    mods = [x for x in S.metabolites(m) if !S.is_reactant(m, x)]
    u = Dict(x => 1 + 9 * big(rand(rng)) for x in m.substrates)
    setprecision(BigFloat, 1024) do
        map(S.all_subsets(mods)) do sat
            c = Dict{Symbol, BigFloat}(x => x in m.substrates ? X * u[x] : BigFloat(10)^-80
                                       for x in S.metabolites(m))
            foreach(x -> c[x] = Y, sat)
            f(c)
        end
    end
end
