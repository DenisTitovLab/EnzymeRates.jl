# ABOUTME: Joint log-linear solve (module Thermo) of cycle and kinetic-group constraints over bare
# ABOUTME: per-step constants, with a content-only priority; P2 probe, stdlib only, checked Int64.
isdefined(Main, :Pop) || include(joinpath(@__DIR__, "Population.jl"))

"""
One exact solve over the log of every bare step constant of a plain mechanism.

**Columns.** One per bare constant, named by content so that storage order, storage
direction and grouping never enter a name:
- a steady-state (SS) step has two directed rate constants, `k:<side> -> <side>`;
- a rapid-equilibrium (RE) step that binds one metabolite (one free metabolite on one
  side, none on the other) has its dissociation constant `Kd:<form> + M <-> <form>`;
- any other RE step has `K:<side> -> <side>` in its canonical orientation (see below).
A side is a form key followed by its free metabolites; a form key is
`conformation(bound...;+residual-residual)` with `X::Inh` for an inhibitor copy.

**Rows.** The fundamental cycles of a spanning forest of the full step graph
(`Σ z_j log Ka_j = c log Keq`, c the cycle's net metabolite change over the net
reaction, an error if it is not a multiple), plus one equality row per tied constant of a
kinetic group (each member against the group's lead, the member with the smallest
content key, oriented alike). `AllParameters` keeps the group rows only.

**Priority** (a function of each step's own content only). A constant is eliminated
earlier when its key `(kind, occupied, ss, against)` is larger, ties broken by the
constant's name (the larger name is eliminated first, so the smaller one is kept):
- `kind`: 3 an isomerization (no free metabolite), 2 any other chemistry (fused binding,
  Theorell-Chance, several metabolites on a side, residual change), 1 a plain binding
  (one metabolite taken up, the bound form holds the free form's metabolites plus it);
- `occupied`: 1 when no form of the step is unbound (a binding: its free side's form
  holds a metabolite), else 0, so free-enzyme constants stay independent longest;
- `ss`: 1 for an SS rate constant, 0 for an RE constant;
- `against`: 1 for the SS rate constant that drives against the physical forward
  direction: the direction of decreasing `progress` (substrates remaining, then products
  made, then free substrates minus free products), association for a binding that
  progress leaves tied (an inhibitor), the lexicographically smaller orientation for any
  other tie.

The reduced row echelon form of `[A | rhs]` over `Rational{Int64}` (every operation
checked) with the columns in priority order gives the dependents (pivots) as power
products of the independents and `Keq`. A dependent whose right side is one independent
with exponent 1 and no `Keq` is a fold (a rename). Labels: `:kinetic_group` when the
dependent ties to a constant of its own kinetic group that is kept longer (the class
representative), else `:haldane` when its right side holds `Keq`, else `:wegscheider`.

Usage (plain mechanisms only):
```julia
include(".../proto/common/Thermo.jl")          # defines Main.Pop and Main.Thermo
sol = Thermo.solve(m)                         # or solve(m; groups = …, mode = :all)
Thermo.independent_names(sol); Thermo.dependent_table(sol)   # content-keyed names
x = Thermo.all_values(sol, Dict(i => v for i in sol.indep), Keq)   # every bare constant
k = Thermo.stored_constants(sol, x)           # per step, Oracle convention
Thermo.zero_flux_ss(m); Thermo.with_flux_free_as_re(m)       # flux-free SS steps
```
"""
module Thermo

import ..Pop

const Q = Rational{Int64}

# ─── Content keys ─────────────────────────────────────────────────────

metkey(x::Pop.PMet) = x.role === :Inh ? string(x.name, "::Inh") :
                      x.role === :Reg ? string(x.name, "::Reg") : string(x.name)

function formkey(sp::Pop.PSpecies)
    res = isempty(sp.res_added) && isempty(sp.res_subtracted) ? "" :
          ";" * join(("+" * string(r) for r in sort(sp.res_added)), "") *
          join(("-" * string(r) for r in sort(sp.res_subtracted)), "")
    string(sp.conformation, "(", join(sort(metkey.(sp.bound)), ","), res, ")")
end

sidekey(fk::AbstractString, mets) = join([fk; sort(metkey.(mets))], " + ")

"""Reaction progress of a side (a form plus free metabolites), compared lexicographically, larger
= earlier in the reaction: (substrates held or free, minus products held or free, free
substrates minus free products). The first two entries order chemistry; the third orders a
binding, a release or an exchange that converts nothing."""
function progress(sp::Pop.PSpecies, free)
    ns(v) = count(x -> x.role === :S, v); np(v) = count(x -> x.role === :P, v)
    (ns(sp.bound) + ns(free), -(np(sp.bound) + np(free)), ns(free) - np(free))
end

# ─── Canonical steps ───────────────────────────────────────────────────

"""A step in its canonical orientation `a → b` (form indices) taking up `cons` and giving off
`rel`. `flipped` says the stored step runs `b → a`. `fwd` is +1 when the physical forward
direction is `a → b`, -1 when it is `b → a`. `names` are the bare constants: `[K]` for RE,
`[k(a→b), k(b→a)]` for SS. `bindK` marks an RE binding-like step (a dissociation constant)."""
struct CStep
    a::Int
    b::Int
    cons::Vector{Pop.PMet}
    rel::Vector{Pop.PMet}
    re::Bool
    flipped::Bool
    kind::Int
    occ::Int
    fwd::Int
    bindK::Bool
    key::String
    names::Vector{String}
end

bindinglike(cons, rel) = length(cons) + length(rel) == 1

function plain_binding(m::Pop.PMech, a, b, cons, rel)
    (length(cons) == 1 && isempty(rel)) || return false
    fa, fb = m.species[a], m.species[b]
    sort(metkey.([fa.bound; cons])) == sort(metkey.(fb.bound)) &&
        sort(fa.res_added) == sort(fb.res_added) &&
        sort(fa.res_subtracted) == sort(fb.res_subtracted)
end

"""Orientation of a step from content alone: a binding-like step takes its metabolite up;
otherwise decreasing `progress` (forward); otherwise the lexicographically smaller orientation.
Returns `(a, b, cons, rel, flipped)`."""
function canonical_orientation(m::Pop.PMech, s::Pop.PStep, fk)
    fwd = (s.from, s.to, s.consumed, s.released, false)
    rev = (s.to, s.from, s.released, s.consumed, true)
    if bindinglike(s.consumed, s.released)
        return length(s.consumed) == 1 ? fwd : rev
    end
    pf = progress(m.species[s.from], s.consumed)
    pt = progress(m.species[s.to], s.released)
    pf != pt && return pf > pt ? fwd : rev
    kf = (sidekey(fk[s.from], s.consumed), sidekey(fk[s.to], s.released))
    reverse(kf) < kf ? rev : fwd
end

function canonical_steps(m::Pop.PMech)
    fk = formkey.(m.species)
    map(m.steps) do s
        a, b, cons, rel, flipped = canonical_orientation(m, s, fk)
        sa, sb = sidekey(fk[a], cons), sidekey(fk[b], rel)
        kind = isempty(cons) && isempty(rel) ? 3 : plain_binding(m, a, b, cons, rel) ? 1 : 2
        unb(i) = isempty(m.species[i].bound)
        occ = kind == 1 ? (unb(a) ? 0 : 1) : (unb(a) || unb(b) ? 0 : 1)
        pa, pb = progress(m.species[a], cons), progress(m.species[b], rel)
        fwd = pa > pb ? 1 : pa < pb ? -1 : bindinglike(cons, rel) ? 1 : (sa <= sb ? 1 : -1)
        bk = s.re && bindinglike(cons, rel)
        key = fwd == 1 ? string(sa, " -> ", sb) : string(sb, " -> ", sa)
        names = s.re ? [bk ? string("Kd:", sa, " <-> ", sb) : string("K:", sa, " -> ", sb)] :
                [string("k:", sa, " -> ", sb), string("k:", sb, " -> ", sa)]
        CStep(a, b, cons, rel, s.re, flipped, kind, occ, fwd, bk, string(s.re ? "RE " : "SS ", key),
              names)
    end
end

# ─── Columns and priority ──────────────────────────────────────────────

"""Bare constant: step `j`, slot `r` (1 = K or k(a→b), 2 = k(b→a)), its name and priority key."""
struct Col
    j::Int
    r::Int
    name::String
    prio::NTuple{4, Int}
end

function columns(cs::Vector{CStep})
    cols = Col[]
    for (j, c) in enumerate(cs)
        if c.re
            push!(cols, Col(j, 1, c.names[1], (c.kind, c.occ, 0, 0)))
        else
            # slot 1 drives a → b; it runs against the physical forward direction when fwd = -1
            push!(cols, Col(j, 1, c.names[1], (c.kind, c.occ, 1, c.fwd == 1 ? 0 : 1)))
            push!(cols, Col(j, 2, c.names[2], (c.kind, c.occ, 1, c.fwd == 1 ? 1 : 0)))
        end
    end
    # eliminated first: larger priority key, then the larger name (the smaller name is kept)
    sort!(cols; by = c -> (c.prio, c.name), rev = true)
end

# ─── Cycles ───────────────────────────────────────────────────────────

"""Fundamental cycles of a spanning forest of the step graph, as ±1 vectors over the steps
in canonical orientation, built by breadth-first search from the lowest form index."""
function fundamental_cycles(nforms::Int, cs::Vector{CStep})
    adj = [Tuple{Int, Int, Int}[] for _ in 1:nforms]         # (neighbour, step, sign)
    for (j, c) in enumerate(cs)
        push!(adj[c.a], (c.b, j, 1)); push!(adj[c.b], (c.a, j, -1))
    end
    parent = fill((0, 0, 0), nforms)                        # (parent form, step, sign)
    depth = fill(-1, nforms)
    tree = falses(length(cs))
    for r in 1:nforms
        depth[r] >= 0 && continue
        depth[r] = 0; queue = [r]
        while !isempty(queue)
            u = popfirst!(queue)
            for (v, j, sg) in adj[u]
                depth[v] >= 0 && continue
                depth[v] = depth[u] + 1; parent[v] = (u, j, sg); tree[j] = true
                push!(queue, v)
            end
        end
    end
    cycles = Vector{Int}[]
    for (j, c) in enumerate(cs)
        tree[j] && continue
        z = zeros(Int, length(cs)); z[j] = 1
        # walk b → a through the tree: x climbs from b, y climbs from a
        x, y = c.b, c.a
        while x != y
            if depth[x] >= depth[y]
                p, k, sg = parent[x]; z[k] -= sg; x = p           # traverse x → p
            else
                p, k, sg = parent[y]; z[k] += sg; y = p           # traverse p → y later
            end
        end
        push!(cycles, z)
    end
    cycles
end

"""Net metabolite release of each canonical step, by name (an inhibitor copy counts as its
reactant)."""
function step_delta(cs::Vector{CStep})
    map(cs) do c
        d = Dict{Symbol, Int}()
        for x in c.cons; d[x.name] = get(d, x.name, 0) - 1; end
        for x in c.rel;  d[x.name] = get(d, x.name, 0) + 1; end
        filter!(p -> p.second != 0, d)
    end
end

"""The multiple `c` of the net reaction (products +1, substrates -1) that cycle `z` runs;
errors when the net change is no multiple of it."""
function cycle_multiple(m::Pop.PMech, deltas, z)
    tot = Dict{Symbol, Int}()
    for (j, zj) in enumerate(z), (k, v) in deltas[j]
        zj == 0 || (tot[k] = get(tot, k, 0) + zj * v)
    end
    filter!(p -> p.second != 0, tot)
    isempty(tot) && return 0
    r = Dict{Symbol, Int}()
    for x in m.substrates; r[x] = get(r, x, 0) - 1; end
    for x in m.products;   r[x] = get(r, x, 0) + 1; end
    k0 = first(keys(r))
    c = get(tot, k0, 0) // r[k0]
    (isinteger(c) && c != 0 && all(get(tot, k, 0) == c * v for (k, v) in r) &&
     all(haskey(r, k) for k in keys(tot))) ||
        error("$(m.id): cycle net change $(tot) is no multiple of the net reaction")
    Int(c)
end

"""Row coefficients of cycle `z` over the bare constants: Σ z_j log Ka_j where Ka is the
association constant in canonical orientation (SS: k(a→b)/k(b→a); RE binding: 1/Kd; RE
other: K)."""
function cycle_row(cs, colidx, z)
    row = zeros(Int, length(colidx))
    for (j, zj) in enumerate(z)
        zj == 0 && continue
        c = cs[j]
        if c.re
            row[colidx[(j, 1)]] += c.bindK ? -zj : zj
        else
            row[colidx[(j, 1)]] += zj; row[colidx[(j, 2)]] -= zj
        end
    end
    row
end

# ─── Kinetic groups ────────────────────────────────────────────────────

"""For each step, its kinetic group's lead (the member with the smallest content key) and
whether the member must be reversed from its canonical orientation to run like the lead:
the same metabolites taken up and given off, then the same conformation change, else its
canonical orientation. `groups` holds step indices."""
function group_alignment(m::Pop.PMech, cs, groups)
    msig(c, flip) = flip ? (sort(metkey.(c.rel)), sort(metkey.(c.cons))) :
                           (sort(metkey.(c.cons)), sort(metkey.(c.rel)))
    csig(c, flip) = flip ? (m.species[c.b].conformation, m.species[c.a].conformation) :
                           (m.species[c.a].conformation, m.species[c.b].conformation)
    al = Dict{Int, Tuple{Int, Bool}}()
    for g in groups
        lead = g[argmin([cs[j].key for j in g])]
        for j in g
            cs[j].re == cs[lead].re || error("$(m.id): group mixes RE and SS steps")
            fits = [f for f in (false, true) if msig(cs[j], f) == msig(cs[lead], false)]
            isempty(fits) && error("$(m.id): group members $(lead), $(j) differ in kind")
            if length(fits) == 2
                fc = [f for f in fits if csig(cs[j], f) == csig(cs[lead], false)]
                fits = length(fc) == 1 ? fc : [false]
            end
            al[j] = (lead, only(fits))
        end
    end
    al
end

"""Equality rows of the kinetic groups: every member's constants equal its lead's, both
oriented alike (`group_alignment`)."""
function group_rows(cs, colidx, al)
    rows = Vector{Int}[]
    for (j, (lead, flip)) in sort!(collect(al); by = first)
        j == lead && continue
        c = cs[j]
        if c.re
            row = zeros(Int, length(colidx))
            row[colidx[(j, 1)]] += 1
            # an RE binding K is direction-free; any other RE K inverts on a flip
            row[colidx[(lead, 1)]] -= (flip && !c.bindK) ? -1 : 1
            push!(rows, row)
        else
            for r in 1:2
                row = zeros(Int, length(colidx))
                row[colidx[(j, flip ? 3 - r : r)]] += 1
                row[colidx[(lead, r)]] -= 1
                push!(rows, row)
            end
        end
    end
    rows
end

# ─── Exact RREF ────────────────────────────────────────────────────────

"""Reduced row echelon form of `M` in place, every operation checked. Returns the pivot
columns and the largest numerator or denominator met on the way."""
function rref!(M::Matrix{Q})
    nr, nc = size(M)
    piv = Int[]
    big = 1
    row = 1
    for col in 1:nc
        row > nr && break
        p = findfirst(r -> M[r, col] != 0, row:nr)
        p === nothing && continue
        p += row - 1
        if p != row
            for k in 1:nc; M[row, k], M[p, k] = M[p, k], M[row, k]; end
        end
        inv = 1 // M[row, col]
        for k in 1:nc
            M[row, k] *= inv
            big = max(big, abs(numerator(M[row, k])), denominator(M[row, k]))
        end
        for r in 1:nr
            (r == row || M[r, col] == 0) && continue
            f = M[r, col]
            for k in 1:nc
                M[row, k] == 0 && continue
                M[r, k] -= f * M[row, k]
                big = max(big, abs(numerator(M[r, k])), denominator(M[r, k]))
            end
        end
        push!(piv, col); row += 1
    end
    piv, big
end

# ─── The solve ─────────────────────────────────────────────────────────

"""A dependent: column `col` = Keq^keq · Π independent^exponent."""
struct Dep
    col::Int
    keq::Q
    terms::Vector{Pair{Int, Q}}
    label::Symbol
    fold::Bool
end

"""Result of one solve. `cols` in priority order; `indep` and `deps` index into it."""
struct Solution
    m::Pop.PMech
    cs::Vector{CStep}
    cols::Vector{Col}
    indep::Vector{Int}
    deps::Vector{Dep}
    ncycles::Int
    nhaldane_fund::Int
    haldane_gcd::Int
    maxint::Int
end

"""Joint solve of `m` under `groups` (default its own). `mode = :independent` uses cycle and
group rows; `mode = :all` the group rows only."""
function solve(m::Pop.PMech; groups = m.groups, mode::Symbol = :independent)
    cs = canonical_steps(m)
    cols = columns(cs)
    colidx = Dict((c.j, c.r) => i for (i, c) in enumerate(cols))
    rows = Vector{Int}[]
    rhs = Int[]
    cyc = fundamental_cycles(length(m.species), cs)
    deltas = step_delta(cs)
    mult = [cycle_multiple(m, deltas, z) for z in cyc]
    if mode === :independent
        for (z, c) in zip(cyc, mult)
            push!(rows, cycle_row(cs, colidx, z)); push!(rhs, c)
        end
    end
    al = group_alignment(m, cs, groups)
    for r in group_rows(cs, colidx, al)
        push!(rows, r); push!(rhs, 0)
    end
    n = length(cols)
    M = zeros(Q, length(rows), n + 1)
    for (i, r) in enumerate(rows)
        M[i, 1:n] .= r; M[i, n + 1] = rhs[i]
    end
    piv, big = rref!(M)
    (n + 1) in piv && error("$(m.id): the constraints reduce to 0 = log Keq")
    indep = setdiff(1:n, piv)
    # class representative of a tied constant: the member of its kinetic-group role class
    # kept longest (last in priority order)
    gos = zeros(Int, length(m.steps))
    for (k, g) in enumerate(groups), j in g; gos[j] = k; end
    rep = Dict{Tuple{Int, Int}, Int}()
    for (i, c) in enumerate(cols)
        rep[(gos[c.j], role(cs, al, c))] = i           # later columns overwrite earlier ones
    end
    # a tied class keeps at most one free member, and it is the representative
    for (i, c) in enumerate(cols)
        i in indep && rep[(gos[c.j], role(cs, al, c))] != i &&
            error("$(m.id): free member $(c.name) is not its class representative")
    end
    deps = Dep[]
    for (r, p) in enumerate(piv)
        terms = [f => -M[r, f] for f in indep if M[r, f] != 0]
        keq = M[r, n + 1]
        fold = keq == 0 && length(terms) == 1 && terms[1].second == 1
        c = cols[p]
        label = rep[(gos[c.j], role(cs, al, c))] != p ? :kinetic_group :
                keq != 0 ? :haldane : :wegscheider
        push!(deps, Dep(p, keq, terms, label, fold))
    end
    g = foldl(gcd, mult; init = 0)
    Solution(m, cs, cols, indep, deps, length(cyc), count(!=(0), mult), g, big)
end

"""Role of a constant inside its kinetic group, aligned with the group lead
(`group_alignment`): 1 for an RE K; for SS the lead's slot it is tied to."""
role(cs, al, c::Col) = cs[c.j].re ? 1 : (al[c.j][2] ? 3 - c.r : c.r)

# ─── Names, expressions and numbers ────────────────────────────────────

colname(sol::Solution, i::Int) = sol.cols[i].name

"""Every dependent as `name => (Keq exponent, [name => exponent])`, sorted by name."""
function dependent_table(sol::Solution)
    sort!([colname(sol, d.col) => (d.keq, sort!([colname(sol, f) => e for (f, e) in d.terms];
                                                  by = first)) for d in sol.deps]; by = first)
end

independent_names(sol::Solution) = sort!([colname(sol, i) for i in sol.indep])

"""Values of every bare constant (BigFloat, by column) from values of the independents
(`vals[i]` for column `i` in `sol.indep`) and Keq."""
function all_values(sol::Solution, indep_vals::AbstractDict{Int}, keq)
    x = Dict{Int, BigFloat}(i => BigFloat(indep_vals[i]) for i in sol.indep)
    for d in sol.deps
        v = BigFloat(keq)^BigFloat(d.keq)
        for (f, e) in d.terms
            v *= x[f]^BigFloat(e)
        end
        x[d.col] = v
    end
    x
end

"""Per-step constants in the stored orientation and the oracle's convention (`[K]` RE,
`[kf, kr]` SS) from values of every bare constant by column."""
function stored_constants(sol::Solution, x::AbstractDict{Int})
    idx = Dict((c.j, c.r) => i for (i, c) in enumerate(sol.cols))
    map(enumerate(sol.cs)) do (j, c)
        if c.re
            K = x[idx[(j, 1)]]
            [c.flipped && !c.bindK ? 1 / K : K]
        else
            k1, k2 = x[idx[(j, 1)]], x[idx[(j, 2)]]
            c.flipped ? [k2, k1] : [k1, k2]
        end
    end
end

"""Values of every bare constant by column from stored per-step constants (the inverse of
`stored_constants`)."""
function bare_values(sol::Solution, k::AbstractVector)
    x = Dict{Int, BigFloat}()
    for (i, c) in enumerate(sol.cols)
        s = sol.cs[c.j]; kj = k[c.j]
        x[i] = s.re ? (s.flipped && !s.bindK ? 1 / BigFloat(kj[1]) : BigFloat(kj[1])) :
               BigFloat(kj[s.flipped ? 3 - c.r : c.r])
    end
    x
end

# ─── Zero-flux steady-state steps ──────────────────────────────────────

"""SS steps that carry no net flux at steady state for any constants that satisfy the cycle
relations: an SS step whose every cycle runs no reaction. Contract every RE segment to a
node; an SS step inside one segment closes a cycle with that segment's RE path, and an SS
step between segments lies in a biconnected block of the contracted graph (Tarjan). A
step in a block, or a self-loop, whose cycles all have net change zero is flux-free: its
two constants enter the rate only as their ratio. Returns a Bool per step."""
function zero_flux_ss(m::Pop.PMech; cs = canonical_steps(m))
    n = length(m.species)
    deltas = step_delta(cs)
    # RE segments by union-find, with each form's net-release potential from its root
    par = collect(1:n)
    find(x) = (while par[x] != x; par[x] = par[par[x]]; x = par[x]; end; x)
    for c in cs
        c.re && (par[find(c.a)] = find(c.b))
    end
    phi = Vector{Union{Nothing, Dict{Symbol, Int}}}(nothing, n)
    for r in 1:n
        (find(r) == r) || continue
        phi[r] = Dict{Symbol, Int}(); queue = [r]
        while !isempty(queue)
            u = popfirst!(queue)
            for (j, c) in enumerate(cs)
                c.re || continue
                for (x, y, sg) in ((c.a, c.b, 1), (c.b, c.a, -1))
                    (x == u && phi[y] === nothing) || continue
                    d = copy(phi[u])
                    for (k, v) in deltas[j]; d[k] = get(d, k, 0) + sg * v; end
                    phi[y] = d; push!(queue, y)
                end
            end
        end
    end
    # contracted net change of SS step j taken a → b: Δ_j + φ(a) − φ(b)
    function cdelta(j)
        c = cs[j]
        d = copy(deltas[j])
        for (k, v) in phi[c.a]; d[k] = get(d, k, 0) + v; end
        for (k, v) in phi[c.b]; d[k] = get(d, k, 0) - v; end
        filter!(p -> p.second != 0, d)
    end
    ss = [j for (j, c) in enumerate(cs) if !c.re]
    zero = falses(length(cs))
    for j in ss
        find(cs[j].a) == find(cs[j].b) && (zero[j] = isempty(cdelta(j)))
    end
    # biconnected blocks over the edges between segments (Tarjan, iterative-free recursion)
    edges = [j for j in ss if find(cs[j].a) != find(cs[j].b)]
    nodes = unique(vcat([find(cs[j].a) for j in edges], [find(cs[j].b) for j in edges]))
    adj = Dict(u => Tuple{Int, Int}[] for u in nodes)        # (neighbour, edge)
    for j in edges
        u, v = find(cs[j].a), find(cs[j].b)
        push!(adj[u], (v, j)); push!(adj[v], (u, j))
    end
    disc = Dict{Int, Int}(); low = Dict{Int, Int}(); t = Ref(0)
    stack = Int[]; blocks = Vector{Vector{Int}}()
    function dfs(u, pe)
        t[] += 1; disc[u] = low[u] = t[]
        for (v, e) in adj[u]
            e == pe && continue
            if !haskey(disc, v)
                push!(stack, e)
                dfs(v, e)
                low[u] = min(low[u], low[v])
                if low[v] >= disc[u]
                    blk = Int[]
                    while true
                        x = pop!(stack); push!(blk, x); x == e && break
                    end
                    push!(blocks, blk)
                end
            elseif disc[v] < disc[u]
                push!(stack, e); low[u] = min(low[u], disc[v])
            end
        end
    end
    for u in sort(nodes)
        haskey(disc, u) || dfs(u, 0)
    end
    for blk in blocks
        # a block runs a reaction when one of its fundamental cycles has nonzero net change
        reactive = false
        bn = unique(vcat([find(cs[j].a) for j in blk], [find(cs[j].b) for j in blk]))
        loc = Dict(u => i for (i, u) in enumerate(bn))
        # spanning tree of the block by BFS; potentials ψ over metabolites
        psi = Vector{Union{Nothing, Dict{Symbol, Int}}}(nothing, length(bn))
        psi[1] = Dict{Symbol, Int}(); queue = [1]; used = Set{Int}()
        while !isempty(queue)
            u = popfirst!(queue)
            for j in blk
                a, b = loc[find(cs[j].a)], loc[find(cs[j].b)]
                for (x, y, sg) in ((a, b, 1), (b, a, -1))
                    (x == u && psi[y] === nothing) || continue
                    d = copy(psi[u])
                    for (k, v) in cdelta(j); d[k] = get(d, k, 0) + sg * v; end
                    psi[y] = d; push!(used, j); push!(queue, y)
                end
            end
        end
        for j in blk
            j in used && continue
            a, b = loc[find(cs[j].a)], loc[find(cs[j].b)]
            d = mergewith(+, copy(psi[a]), cdelta(j))
            for (k, v) in psi[b]; d[k] = get(d, k, 0) - v; end
            any(!=(0), values(d)) && (reactive = true)
        end
        reactive || (for j in blk; zero[j] = true; end)
    end
    zero
end

"""`m` with every kinetic group of flux-free SS steps (`zero_flux_ss`) made RE: their two
rate constants enter the rate only as a ratio, which an RE step's one constant carries."""
function with_flux_free_as_re(m::Pop.PMech)
    z = zero_flux_ss(m)
    conv = falses(length(m.steps))
    for g in m.groups
        all(j -> !m.steps[j].re && z[j], g) && (conv[g] .= true)
    end
    any(conv) || return m, 0
    steps = [conv[j] ? Pop.PStep(s.from, s.to, s.consumed, s.released, true) : s
             for (j, s) in enumerate(m.steps)]
    m2 = Pop.PMech(; id = m.id, source = m.source, substrates = m.substrates,
                   products = m.products, inhibitors = m.inhibitors, metabolites = m.metabolites,
                   species = m.species, steps = steps, groups = m.groups)
    m2, count(g -> conv[first(g)], m.groups)
end

end # module
