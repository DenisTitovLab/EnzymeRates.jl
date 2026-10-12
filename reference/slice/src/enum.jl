# ABOUTME: Enumeration moves on the mechanism type: the steady-state flip, the group split
# ABOUTME: and the dead-end binding, with the graph tests and the copy rule they share.

# ─── Candidates and the minimal-set search ──────────────────────────────────────────────

"""`m` with edge `k` in kinetic group `labels[k]`, group `l` tagged `tags[l]`, and edge
`k` at rapid equilibrium exactly when `re[k]`: the groups renumbered by their first edge.
No two edges share their forms and metabolites, which order them before the flag, so the
fields stay canonical; each group must keep one flag."""
function regroup(m::Mechanism, labels, tags, re)
    group, t = relabel(labels, tags, Int(maximum(labels)))
    Mechanism(m.substrates, m.products, m.mets, m.confs, m.forms,
              [Edge(e.a, e.b, e.cons, e.rel, r) for (e, r) in zip(m.edges, re)], group, t,
              m.multiplicity, m.sites)
end

"""`m` with edge `k` at rapid equilibrium exactly when `re[k]`."""
with_re(m::Mechanism, re) = regroup(m, m.group, m.tags, re)

"""All subsets of `v`, the empty one first."""
all_subsets(v) =
    [v[[isodd(s >> (i - 1)) for i in eachindex(v)]] for s in 0:(1 << length(v)) - 1]

"""
    minimal_sets(gains, partners) -> sets

Every inclusion-minimal set of units for which `gains(set)` holds, level by level: a set
is tested only when each of its subsets one smaller was tested and failed, and a failed
set grows only by the units `partners(set)` lists. Each set is sorted.
"""
function minimal_sets(gains, partners)
    out, failed = Vector{Int}[], [Int[]]
    while !isempty(failed)
        tested, next = Set(failed), Set{Vector{Int}}()
        for set in failed, u in partners(set)
            u in set && continue
            c = sort!([set; u])
            all(x -> filter(!=(x), c) in tested, c) && push!(next, c)
        end
        failed = Vector{Int}[]
        for c in sort!(collect(next))
            gains(c) ? push!(out, c) : push!(failed, c)
        end
    end
    out
end

# ─── Graph tests ─────────────────────────────────────────────────────────────────────────

"""The nibble mask of `m`'s modifier entries."""
modifier_mask(m::Mechanism) =
    sum((Pack(0xf) << (4 * (i - 1)) for (i, x) in enumerate(m.mets) if x.modifier);
        init = zero(Pack))

"""
    edge_flux(m) -> BitVector

Whether each edge carries net steady-state flux for generic constants (a rapid-equilibrium
edge never does). A steady-state edge `a → b` weighs its uptake of substrates minus
products plus the offsets (`segment_offsets`) of `a` minus those of `b`, so a cycle weighs
the turnover it runs times the reactant count. An edge within one segment carries flux when
its weight is nonzero, any other when its block of the segment graph holds a cycle of
nonzero weight (`unbalanced_blocks`).
"""
function edge_flux(m::Mechanism)
    segment, names, ex = segment_offsets(m)
    rho = [x in m.substrates ? 1 : x in m.products ? -1 : 0 for x in names]
    sgn = [rho[searchsortedfirst(names, x.name)] for x in m.mets]
    pot = ex' * rho
    up(p) = sum(sgn[j] * count_at(p, j) for j in eachindex(m.mets))
    w(e) = up(e.cons) - up(e.rel) + pot[e.a] - pot[e.b]
    flux, cross = falses(length(m.edges)), Int[]
    for (k, e) in enumerate(m.edges)
        e.re && continue
        segment[e.a] == segment[e.b] ? (flux[k] = w(e) != 0) : push!(cross, k)
    end
    block, unbalanced = unbalanced_blocks(maximum(segment),
        [(segment[m.edges[k].a], segment[m.edges[k].b]) for k in cross],
        [w(m.edges[k]) for k in cross])
    flux[cross] .= in(unbalanced).(block)
    flux
end

"""Whether some edge of each kinetic group carries flux (`edge_flux`)."""
group_flux(m::Mechanism) =
    (f = edge_flux(m); [any(f[members(m, g)]) for g in 1:ngroups(m)])

"""
The kinetic groups that flank a steady-state isomerization alone in its group, when each
of its forms has one other edge, a binding into that form alone in its group: by the chain
lemma flipping such a flank adds constants and no rate law. None in an allosteric
mechanism, where a flank's flip can show through the two conformations.
"""
function chain_flanks(m::Mechanism)
    is_allosteric(m) && return Int[]
    alone = [count(==(g), m.group) == 1 for g in m.group]
    at = [Int[] for _ in m.forms]
    for (k, e) in enumerate(m.edges)
        push!(at[e.a], k)
        push!(at[e.b], k)
    end
    function flank(k0, x)
        rest = filter(!=(k0), at[x])
        length(rest) == 1 && is_binding(m.edges[rest[1]]) && m.edges[rest[1]].b == x &&
            alone[rest[1]] ? Int(m.group[rest[1]]) : 0
    end
    out = Int[]
    for (k, e) in enumerate(m.edges)
        alone[k] && e.cons == e.rel == 0 && !e.re || continue
        f1, f2 = flank(k, e.a), flank(k, e.b)
        f1 > 0 && f2 > 0 && push!(out, f1, f2)
    end
    out
end

"""Whether every vertex of the digraph `arcs` (`(u, v, …)` on the vertices `1:n`) reaches
`root` when each arc in `fixed` is its source's only way out, so that a spanning
arborescence toward `root` holds every fixed arc."""
function reaches_all(n, arcs, root, fixed)
    pinned = Dict(a[1] => a[2] for a in fixed)
    into = [Int[] for _ in 1:n]
    for (u, v) in arcs
        get(pinned, u, v) == v && push!(into[v], u)
    end
    seen, queue = falses(n), [root]
    seen[root] = true
    while !isempty(queue)
        for u in into[popfirst!(queue)]
            seen[u] || (seen[u] = true; push!(queue, u))
        end
    end
    all(seen)
end

"""
    hyperbolic(m) -> Bool

Whether the denominator of `m`'s catalytic scheme has degree at most 1 in every substrate
and product, leaving out every edge at a form that holds a modifier. On the graph of the
rapid-equilibrium segments, a steady-state edge taken either way scores for a metabolite
the copies it takes up that way plus its source form's offset (`segment_offsets`), and a
segment scores its largest offset. The degree exceeds 1 when an edge or a segment scores
2, or some spanning arborescence toward a scoring segment holds a scoring edge, or one
holds two (`reaches_all`).
"""
function hyperbolic(m::Mechanism)
    mods = modifier_mask(m)
    mask = [(m.forms[e.a].bound | m.forms[e.b].bound) & mods == 0 for e in m.edges]
    segment, names, ex = segment_offsets(m, mask)
    on = falses(length(m.forms))
    foreach(k -> mask[k] && (on[m.edges[k].a] = on[m.edges[k].b] = true), eachindex(mask))
    ids = unique(segment[on])
    sid = [on[i] ? findfirst(==(segment[i]), ids) : 0 for i in eachindex(m.forms)]
    n = length(ids)
    arcs = Tuple{Int, Int, Int, Pack}[]
    for (k, e) in enumerate(m.edges)
        mask[k] && !e.re && sid[e.a] != sid[e.b] || continue
        push!(arcs, (sid[e.a], sid[e.b], e.a, e.cons), (sid[e.b], sid[e.a], e.b, e.rel))
    end
    for x in [m.substrates; m.products]
        r, j = searchsortedfirst(names, x), searchsortedfirst(m.mets, Met(x, false))
        score(a) = count_at(a[4], j) + ex[r, a[3]]
        hot = filter(a -> score(a) > 0, arcs)
        top = zeros(Int, n)
        foreach(i -> top[sid[i]] = max(top[sid[i]], ex[r, i]), findall(on))
        (any(a -> score(a) > 1, hot) || any(>(1), top)) && return false
        roots = findall(==(1), top)
        any(a[1] != s && reaches_all(n, arcs, s, (a,)) for a in hot for s in roots) &&
            return false
        for (p, a) in enumerate(hot), b in hot[p + 1:end]
            a[1] != b[1] && any(1:n) do s
                s != a[1] && s != b[1] && reaches_all(n, arcs, s, (a, b))
            end && return false
        end
    end
    true
end

# ─── The copy rule ──────────────────────────────────────────────────────────────────────

"""The table entry that the binding `o` takes up."""
taken(m::Mechanism, o::Edge) = findfirst(i -> count_at(o.cons, i) > 0, eachindex(m.mets))

"""The modifier entry the oriented edge `o` binds alone; 0 for any other edge."""
copy_entry(m::Mechanism, o::Edge) =
    is_binding(o) && m.mets[taken(m, o)].modifier ? taken(m, o) : 0

"""
    twin_finder(m, mask) -> (site, x) -> form or 0

The productive twin of the complex that entry `x` forms at form `site`, over the edges in
`mask` (one conformation's graph): a form holding no modifier in `site`'s rapid-equilibrium
segment whose offsets are `site`'s plus one `x` by name (`segment_offsets`), so that its
weight is a constant multiple of the complex's. The twin holding the complex's metabolites
by name, its conformation and residual when there is one, else the first; 0 when none is.
"""
function twin_finder(m::Mechanism, mask)
    segment, names, ex = segment_offsets(m, mask)
    mods = modifier_mask(m)
    held(i, extra...) = (sort!([[x.name for x in mets_of(m, m.forms[i].bound)]; extra...]),
                         m.forms[i].conf, m.forms[i].res)
    function (site, x)
        target = ex[:, site]
        target[searchsortedfirst(names, m.mets[x].name)] += 1
        twins = [i for i in eachindex(m.forms) if segment[i] == segment[site] &&
                 m.forms[i].bound & mods == 0 && view(ex, :, i) == target]
        isempty(twins) && return 0
        k = findfirst(i -> held(i) == held(site, m.mets[x].name), twins)
        twins[something(k, 1)]
    end
end

"""
    copy_gauge(m, fwd, mask, g, twin, label) -> rescaling per group, or nothing

The dwell gauge of copy group `g` over the edges in `mask` (one conformation): every
complex of the copy takes the copy's factor; each twin (`twin`, a `twin_finder` of the same
edges) takes the factor of the rapid-equilibrium binding group that forms it from its site
when it is the twin of one complex, else its own; every other form 1. `label` names a
factor class. Each group's rescaling is the factor pair of its edges' ends (none for a
rapid-equilibrium edge whose ends scale alike), and `nothing` when a complex has no twin or
a group's edges need two rescalings.
"""
function copy_gauge(m::Mechanism, fwd, mask, g, twin, label)
    kept(h) = [oriented(m, k, fwd) for k in members(m, h) if mask[k]]
    copies = kept(g)
    x = copy_entry(m, copies[1])
    x > 0 && all(o -> twin(o.a, x) > 0, copies) || return nothing
    twin_of = Dict(Int(o.b) => twin(o.a, x) for o in copies)
    former = Dict{Tuple{Int, Int}, Int}()
    for k in eachindex(m.edges)
        o = oriented(m, k, fwd)
        mask[k] && o.re && is_binding(o) &&
            m.mets[taken(m, o)] == Met(m.mets[x].name, false) &&
            (former[(o.a, o.b)] = m.group[k])
    end
    class = Dict{Int, Any}(c => label((:group, g)) for c in keys(twin_of))
    for o in copies
        t = twin_of[o.b]
        haskey(class, t) && continue
        h = count(==(t), values(twin_of)) == 1 ? get(former, (o.a, t), 0) : 0
        class[t] = label(h > 0 ? (:group, h) : (:form, t))
    end
    rescale(o) = (a = get(class, o.a, nothing); b = get(class, o.b, nothing);
                  o.re && isequal(a, b) ? (nothing, nothing) : (a, b))
    out = Dict{Int, Any}()
    for h in 1:ngroups(m)
        os = kept(h)
        isempty(os) && continue
        r = rescale(os[1])
        all(o -> isequal(rescale(o), r), os) || return nothing
        out[h] = r
    end
    out
end

"""
    redundant_copies(m) -> groups

The copy groups whose constant no rate law sees: in every conformation the copy binds,
each complex has a productive twin and the dwell gauge exists (`copy_gauge`), an `:EqualAI`
group rescaling alike in both. The inactive conformation is `inactive_mask`'s graph; a
factor class is one number in both only when its group and the copy are both `:EqualAI`.
"""
function redundant_copies(m::Mechanism)
    fwd = orient(m)
    copies = [g for g in 1:ngroups(m)
              if copy_entry(m, oriented(m, members(m, g)[1], fwd)) > 0]
    isempty(copies) && return copies
    every = trues(length(m.edges))
    twin = twin_finder(m, every)
    is_allosteric(m) ||
        return filter(g -> copy_gauge(m, fwd, every, g, twin, identity) !== nothing, copies)
    mask = inactive_mask(m)
    twin_i = twin_finder(m, mask)
    filter(copies) do g
        shared(c) = c[1] === :group && m.tags[c[2]] === m.tags[g] === :EqualAI
        label(state) = c -> shared(c) ? c : (state, c)
        ra = copy_gauge(m, fwd, every, g, twin, label(:A))
        ra === nothing && return false
        ri = any(mask[members(m, g)]) ? copy_gauge(m, fwd, mask, g, twin_i, label(:I)) :
             Dict(Int(h) => (nothing, nothing) for h in m.group[mask])
        ri !== nothing && all(keys(ra)) do h
            m.tags[h] !== :EqualAI || !haskey(ri, h) || isequal(ra[h], ri[h])
        end
    end
end

# ─── The moves ──────────────────────────────────────────────────────────────────────────

"""
    flip_children(m) -> Vector{Mechanism}

The rapid-equilibrium to steady-state move. A unit is a rapid-equilibrium kinetic group
that binds no modifier, flanks no chain (`chain_flanks`) and carries flux once every edge is
at steady state (`group_flux`). One child per inclusion-minimal set of units
(`minimal_sets`) whose flip raises the segment count, leaves no bottomless segment and
leaves every flipped group carrying flux. Over more than one catalytic subunit the child
must be `hyperbolic`.
"""
function flip_children(m::Mechanism)
    lead = [m.edges[members(m, g)[1]] for g in 1:ngroups(m)]
    flanks, mods = chain_flanks(m), modifier_mask(m)
    flux = group_flux(with_re(m, falses(length(m.edges))))
    units = [g for (g, e) in enumerate(lead)
             if e.re && flux[g] && !(g in flanks) && (e.cons | e.rel) & mods == 0]
    base = maximum(first(re_forest(m)))
    flipped(sel) =
        with_re(m, [e.re && !(g in units[sel]) for (e, g) in zip(m.edges, m.group)])
    function passes(sel)
        c = flipped(sel)
        maximum(first(re_forest(c))) > base && bottomless(c) === nothing &&
            all(group_flux(c)[units[sel]])
    end
    children = [flipped(sel) for sel in minimal_sets(passes, _ -> eachindex(units))]
    m.multiplicity > 1 ? filter(hyperbolic, children) : children
end

"""
The ways to divide kinetic group `g` in two by context: whether the form each member binds
at (its `a` form, or its `b` form when it only gives metabolites off) holds a given ligand
other than the group's own, has a given conformation, or a given residual. Each division is
`(the part holding the group's first edge, the rest)`, once.
"""
function bipartitions(m::Mechanism, fwd, g)
    mem = members(m, g)
    os = [oriented(m, k, fwd) for k in mem]
    fs = [m.forms[o.cons == 0 && o.rel != 0 ? o.b : o.a] for o in os]
    own = is_binding(os[1]) ? os[1].cons : zero(Pack)
    tests = [[f -> count_at(f.bound, j) > 0 for j in eachindex(m.mets)
              if count_at(own, j) == 0];
             [f -> f.conf == c for c in unique(f.conf for f in fs)];
             [f -> f.res == r for r in unique(f.res for f in fs)]]
    parts = Tuple{Vector{Int}, Vector{Int}}[]
    for t in tests
        with = mem[t.(fs)]
        0 < length(with) < length(mem) || continue
        rest = setdiff(mem, with)
        push!(parts, mem[1] in with ? (with, rest) : (rest, with))
    end
    unique!(parts)
end

"""
    split_children(m) -> Vector{Mechanism}

The kinetic-group split. A unit divides one group by context (`bipartitions`); a
steady-state part none of whose edges carries flux (`edge_flux` of `m`) goes to rapid
equilibrium. One child per inclusion-minimal set of units, at most one per group, whose
joint application raises the independent-parameter count (`solve`) and, when a part went to
rapid equilibrium, leaves no bottomless segment. A failed set grows only by units of
rapid-equilibrium groups that share a segment with its own rapid-equilibrium groups. Both
parts keep the group's tag; a child with a redundant copy (`redundant_copies`) is dropped.
"""
function split_children(m::Mechanism)
    fwd, flux, G = orient(m), edge_flux(m), ngroups(m)
    units = [(g, p..., [k for q in p if !m.edges[q[1]].re && !any(flux[q]) for k in q])
             for g in 1:G for p in bipartitions(m, fwd, g)]
    function child(sel)
        labels, tags, re = Int.(m.group), copy(m.tags), [e.re for e in m.edges]
        for (i, (g, _, part, rev)) in enumerate(units[sel])
            labels[part] .= G + i
            is_allosteric(m) && push!(tags, m.tags[g])
            re[rev] .= true
        end
        regroup(m, labels, tags, re)
    end
    base = length(solve(m).indep)
    function gains(sel)
        c = child(sel)
        (all(u -> isempty(units[u][4]), sel) || bottomless(c) === nothing) &&
            length(solve(c).indep) > base
    end
    segment = first(re_forest(m))
    near = [m.edges[members(m, g)[1]].re ?
            Set(segment[m.edges[k].a] for k in members(m, g)) : Set{Int}() for g in 1:G]
    function partners(sel)
        isempty(sel) && return eachindex(units)
        used = [units[u][1] for u in sel]
        touched = union((near[g] for g in used)...)
        [u for (u, (g,)) in enumerate(units)
         if !(g in used) && !isdisjoint(near[g], touched)]
    end
    filter!(c -> isempty(redundant_copies(c)),
            [child(sel) for sel in minimal_sets(gains, partners)])
end

"""
    dead_end_children(m, declared) -> Vector{Mechanism}

The dead-end move. Each name of `declared` that `m` binds neither as a modifier nor at a
regulatory site binds as a modifier at every form of a placement, at rapid equilibrium in
one new group (`:EqualAI` in an allosteric mechanism), and each edge between two
placement forms is mirrored between their complexes in its own group. A placement comes
from a choice of competing substrates and products (at least one of each) and bound
modifiers: the forms where one of them is free on an edge, holding none of them, nor every
substrate or every product. A child whose copy is redundant (`redundant_copies`) is
dropped.
"""
function dead_end_children(m::Mechanism, declared)
    bound = inhibitors(m)
    names = [x for x in declared if !(x in bound) && !(x in regulators(m))]
    isempty(names) && return Mechanism[]
    holds(i, modifier) =
        [x.name for x in mets_of(m, m.forms[i].bound) if x.modifier == modifier]
    productive = [holds(i, false) for i in eachindex(m.forms)]
    inhibiting = [holds(i, true) for i in eachindex(m.forms)]
    function free_at(x::Met)
        j = findfirst(==(x), m.mets)
        j === nothing && return Int[]
        [Int(i) for e in m.edges for (i, p) in ((e.a, e.cons), (e.b, e.rel))
         if count_at(p, j) > 0]
    end
    places = Vector{Int}[]
    for cs in all_subsets(m.substrates)[2:end], cp in all_subsets(m.products)[2:end],
        ci in all_subsets(bound)
        comp = [cs; cp]
        at = union((free_at(Met(x, false)) for x in comp)...,
                   (free_at(Met(x, true)) for x in ci)...)
        push!(places, sort!([i for i in at if !(m.substrates ⊆ productive[i]) &&
                             !(m.products ⊆ productive[i]) &&
                             isdisjoint(productive[i], comp) &&
                             isdisjoint(inhibiting[i], ci)]))
    end
    unique!(filter!(!isempty, places))
    (subs, prods, steps, group), kw = raw(m)
    form = Dict{Int, RawForm}()
    for (k, e) in enumerate(m.edges)
        form[e.a], form[e.b] = steps[k].from, steps[k].to
    end
    children = Mechanism[]
    for x in names, place in places
        y = Met(x, true)
        complex(i) = (f = form[i]; RawForm(f.conf, [f.bound; y], f.added, f.subtracted))
        mirrors = [k for (k, e) in enumerate(m.edges) if e.a in place && e.b in place]
        new = [[RawStep(complex(m.edges[k].a), complex(m.edges[k].b), steps[k].consumed,
                        steps[k].released, m.edges[k].re) for k in mirrors];
               [RawStep(form[i], complex(i), [y], Met[], true) for i in place]]
        child = Mechanism(subs, prods, [steps; new],
                          [group; group[mirrors]; fill(ngroups(m) + 1, length(place))];
                          tags = is_allosteric(m) ? [kw.tags; :EqualAI] : kw.tags,
                          kw.multiplicity, kw.sites)
        isempty(redundant_copies(child)) && push!(children, child)
    end
    children
end

"""
    expand(m, declared, atoms) -> Vector{Mechanism}

One round of expansion of `m`: the children of the three moves, each once, checked at the
boundary (`check_boundary`, and `check_atoms` with the reaction's `atoms`; a child that
fails is an error) and admitted (`admit` drops a rejected allosteric tag set).
"""
function expand(m::Mechanism, declared, atoms)
    out = unique!([flip_children(m); split_children(m); dead_end_children(m, declared)])
    foreach(c -> (check_boundary(c); check_atoms(c, atoms)), out)
    filter!(c -> admit(c) === nothing, out)
end
