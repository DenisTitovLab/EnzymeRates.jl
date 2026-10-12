# ABOUTME: The mechanism value type: content-sorted forms, isbits edges, kinetic groups,
# ABOUTME: allosteric tags and sites; structure-only constructor, orientation, checks, Sig.

# ─── Metabolites and the raw input a DSL or an enumeration move writes ──────────────────

"""
A metabolite as a form holds it or a step exchanges it: its name, and `modifier` when it
binds as a non-reacting ligand (an inhibitor `I`, or a reactant's copy `A::Mod`).
"""
struct Met
    name::Symbol
    modifier::Bool
end
Base.isless(x::Met, y::Met) = isless((x.name, x.modifier), (y.name, y.modifier))
Base.hash(x::Met, h::UInt) = hash(x.modifier, hash(x.name, hash(:Met, h)))

"""An enzyme form as written: conformation, bound metabolites and residual (substrates
added, products subtracted), each list in any order."""
struct RawForm
    conf::Symbol
    bound::Vector{Met}
    added::Vector{Symbol}
    subtracted::Vector{Symbol}
end

"""A step as written, in either direction: `from` takes up `consumed` and gives off
`released` to become `to`; `re` marks rapid equilibrium."""
struct RawStep
    from::RawForm
    to::RawForm
    consumed::Vector{Met}
    released::Vector{Met}
    re::Bool
end

"""An allosteric regulatory site: its ligands, each ligand's state and its multiplicity."""
struct Site
    ligands::Vector{Symbol}
    states::Vector{Symbol}
    multiplicity::Int
end
Base.:(==)(s::Site, t::Site) =
    s.ligands == t.ligands && s.states == t.states && s.multiplicity == t.multiplicity
Base.hash(s::Site, h::UInt) =
    hash(s.multiplicity, hash(s.states, hash(s.ligands, hash(:Site, h))))
Base.isless(s::Site, t::Site) =
    isless((s.ligands, s.multiplicity, s.states), (t.ligands, t.multiplicity, t.states))

const TAGS = (:OnlyA, :EqualAI, :NonequalAI)
const SITE_STATES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)

# ─── Packed metabolite counts ───────────────────────────────────────────────────────────

"""Counts of each metabolite-table entry, four bits per entry: at most 32 entries and 15
copies of one entry."""
const Pack = UInt128

count_at(p::Pack, i) = Int((p >> (4 * (i - 1))) & 0xf)
nibsum(p::Pack) = (s = 0; while p != 0; s += Int(p & 0xf); p >>= 4; end; s)

"""The counts of the metabolites `xs`, entries of the sorted table `mets`."""
function pack(mets::Vector{Met}, xs)
    p = zero(Pack)
    for x in xs
        i = searchsortedfirst(mets, x)
        count_at(p, i) == 15 &&
            error("a form or a step side holds $(x.name) more than 15 times; a mechanism " *
                  "holds at most 15 copies of one metabolite there")
        p += one(Pack) << (4 * (i - 1))
    end
    p
end

# ─── The mechanism ──────────────────────────────────────────────────────────────────────

"""An enzyme form: the counts of each table entry it holds (`bound`) and of its residual
(`res`: substrates added and products subtracted, the role giving the sign), and its
conformation, an index into `confs`."""
struct Form
    bound::Pack
    res::Pack
    conf::Int32
end
Base.hash(f::Form, h::UInt) = hash(f.conf, hash(f.res, hash(f.bound, hash(:Form, h))))

"""The content order of forms: fewest bound metabolites, then fewest residual entries."""
formkey(f::Form) = (nibsum(f.bound), nibsum(f.res), f.bound, f.res, f.conf)

"""A step between forms `a` and `b`: going `a → b` it takes up `cons` and gives off `rel`.
Stored with `a < b`, except a binding-like step (one metabolite in all, on one side),
which is stored taking its metabolite up."""
struct Edge
    a::Int32
    b::Int32
    cons::Pack
    rel::Pack
    re::Bool
end
Base.isless(e::Edge, f::Edge) =
    isless((e.a, e.b, e.cons, e.rel, e.re), (f.a, f.b, f.cons, f.rel, f.re))
Base.hash(e::Edge, h::UInt) =
    hash(e.re, hash(e.rel, hash(e.cons, hash(e.b, hash(e.a, hash(:Edge, h))))))
reversed(e::Edge) = Edge(e.b, e.a, e.rel, e.cons, e.re)

"""
A mechanism as a graph, plain or allosteric. `forms` are sorted by content (`formkey`) and
`edges` by `(a, b, cons, rel, re)`; `group[k]` is the kinetic group of edge `k`, numbered
by each group's first edge, and `tags[g]` its allosteric tag (none for a plain mechanism).
`mets` and `confs` are the sorted tables that forms and edges index. Build one from raw
steps with `Mechanism(substrates, products, steps, group; tags, multiplicity, sites)`.
"""
struct Mechanism
    substrates::Vector{Symbol}
    products::Vector{Symbol}
    mets::Vector{Met}
    confs::Vector{Symbol}
    forms::Vector{Form}
    edges::Vector{Edge}
    group::Vector{Int32}
    tags::Vector{Symbol}
    multiplicity::Int
    sites::Vector{Site}
end
Base.:(==)(x::Mechanism, y::Mechanism) =
    x.substrates == y.substrates && x.products == y.products && x.mets == y.mets &&
    x.confs == y.confs && x.forms == y.forms && x.edges == y.edges &&
    x.group == y.group && x.tags == y.tags && x.multiplicity == y.multiplicity &&
    x.sites == y.sites
Base.hash(x::Mechanism, h::UInt) = hash(x.sites, hash(x.multiplicity, hash(x.tags,
    hash(x.group, hash(x.edges, hash(x.forms, hash(x.confs, hash(x.mets,
    hash(x.products, hash(x.substrates, hash(:Mechanism, h)))))))))))

"""
    canonicalize(substrates, products, steps, group; tags, multiplicity, sites)
        -> (m, edge_of, flipped)

The mechanism of `steps`, written in any order and direction, where step `j` belongs to
kinetic group `group[j] ∈ 1:G` and, for an allosteric mechanism, group `g` carries
`tags[g]`. `edge_of[j]` is the edge step `j` became and `flipped[j]` whether that edge
runs opposite to the step as written. Checks structure only (see `check_boundary`).
"""
function canonicalize(substrates, products, steps::Vector{RawStep}, group::Vector{Int};
                      tags::Vector{Symbol} = Symbol[], multiplicity::Int = 1,
                      sites::Vector{Site} = Site[])
    subs, prods = sort(substrates), sort(products)
    !isempty(subs) && !isempty(prods) && allunique([subs; prods]) ||
        error("a mechanism needs substrates and products with distinct names")
    seen, named = Met[], Symbol[]
    for s in steps
        append!(seen, s.consumed)
        append!(seen, s.released)
        for f in (s.from, s.to)
            all(in(subs), f.added) && all(in(prods), f.subtracted) ||
                error("a residual adds substrates and subtracts products")
            append!(seen, f.bound)
            foreach(x -> push!(seen, Met(x, false)), f.added)
            foreach(x -> push!(seen, Met(x, false)), f.subtracted)
            push!(named, f.conf)
        end
    end
    mets, confs = sorted_unique(seen), sorted_unique(named)
    length(mets) <= 32 || error("a mechanism holds at most 32 metabolite entries (a name " *
                                "and its modifier copy count apart); this one has " *
                                "$(length(mets))")
    for x in mets
        x.modifier || x.name in subs || x.name in prods ||
            error("$(x.name) is neither a substrate nor a product, so it binds only as a " *
                  "modifier")
    end
    form(f) = Form(pack(mets, f.bound), pack(mets, (Met(x, false) for x in f.added)) +
                                        pack(mets, (Met(x, false) for x in f.subtracted)),
                   Int32(searchsortedfirst(confs, f.conf)))
    n = length(steps)
    ends = Vector{Form}(undef, 2n)
    for (j, s) in enumerate(steps)
        ends[2j - 1], ends[2j] = form(s.from), form(s.to)
    end
    forms, at = Form[], zeros(Int32, 2n)
    for t in sortperm(map(formkey, ends))
        (isempty(forms) || forms[end] != ends[t]) && push!(forms, ends[t])
        at[t] = length(forms)
    end
    written, flipped = Vector{Edge}(undef, n), falses(n)
    for (j, s) in enumerate(steps)
        a, b = at[2j - 1], at[2j]
        a == b && error("step $j joins a form to itself")
        c, r = pack(mets, s.consumed), pack(mets, s.released)
        any(i -> count_at(c, i) > 0 && count_at(r, i) > 0, eachindex(mets)) &&
            error("step $j takes up and gives off one metabolite")
        flipped[j] = nibsum(c) + nibsum(r) == 1 ? c == 0 : a > b
        written[j] = flipped[j] ? Edge(b, a, r, c, s.re) : Edge(a, b, c, r, s.re)
    end
    p = sortperm(written)
    edges = written[p]
    for k in 2:n
        e, f = edges[k - 1], edges[k]
        (e.a, e.b, e.cons, e.rel) == (f.a, f.b, f.cons, f.rel) &&
            error("steps $(p[k - 1]) and $(p[k]) are one reaction; a mechanism holds " *
                  "each reaction once")
    end
    G = maximum(group; init = 0)
    all(in(1:G), group) || error("kinetic group labels run from 1 to the number of groups")
    if isempty(tags)
        isempty(sites) && multiplicity == 1 ||
            error("a plain mechanism has no regulatory sites and multiplicity 1; tag " *
                  "every kinetic group to make it allosteric")
    else
        length(tags) == G && all(in(TAGS), tags) && multiplicity >= 1 ||
            error("an allosteric mechanism tags each kinetic group with one of $TAGS " *
                  "and has a multiplicity of at least 1")
    end
    labels, newtags = relabel(group[p], tags, G)
    csites = sort!(map(canonical_site, sites))
    allunique(l for s in csites for l in s.ligands) ||
        error("a ligand binds at one regulatory site")
    m = Mechanism(subs, prods, mets, confs, forms, edges, labels, newtags, multiplicity,
                  csites)
    check_groups(m, orient(m))
    connected(m) || error("the steps do not form one connected graph")
    m, invperm(p), flipped
end

Mechanism(substrates, products, steps::Vector{RawStep}, group::Vector{Int}; kw...) =
    first(canonicalize(substrates, products, steps, group; kw...))

"""`v` sorted, each value once."""
sorted_unique(v) = (sort!(v); [x for (i, x) in enumerate(v) if i == 1 || x != v[i - 1]])

"""Group labels renumbered by first occurrence in `labels` (one per edge, in edge order),
with `tags` carried along; an error on a label no edge holds."""
function relabel(labels, tags::Vector{Symbol}, G::Int)
    new, group, next = zeros(Int32, G), similar(labels, Int32), Int32(0)
    for (k, l) in enumerate(labels)
        new[l] == 0 && (new[l] = next += Int32(1))
        group[k] = new[l]
    end
    next == G || error("kinetic group $(findfirst(iszero, new)) holds no step")
    t = similar(tags)
    for l in eachindex(tags)
        t[new[l]] = tags[l]
    end
    group, t
end

"""`s` with its ligands sorted, each keeping its state; an error when malformed."""
function canonical_site(s::Site)
    s.multiplicity >= 1 && !isempty(s.ligands) && length(s.ligands) == length(s.states) &&
        all(in(SITE_STATES), s.states) ||
        error("a regulatory site needs ligands, one state from $SITE_STATES per ligand " *
              "and a multiplicity of at least 1")
    p = sortperm(s.ligands)
    Site(s.ligands[p], s.states[p], s.multiplicity)
end

"""Error unless the steps of each kinetic group, in their physical-forward direction, take
up and give off the same metabolites and carry one RE/SS flag."""
function check_groups(m::Mechanism, fwd)
    lead = zeros(Int, ngroups(m))
    for k in eachindex(m.edges)
        g = m.group[k]
        lead[g] == 0 && (lead[g] = k; continue)
        o, l = oriented(m, k, fwd), oriented(m, lead[g], fwd)
        (o.cons, o.rel, o.re) == (l.cons, l.rel, l.re) && continue
        error("a kinetic group holds $(step_text(m, l)) and $(step_text(m, o)), which " *
              "differ in kind or in RE/SS; the steps of a kinetic group share their " *
              "constants, so they must be the same kind of step with the same flag")
    end
end

"""Whether the edges join every form into one graph."""
function connected(m::Mechanism)
    root = collect(eachindex(m.forms))
    top(i) = (while root[i] != i; i = root[i]; end; i)
    foreach(e -> root[top(e.a)] = top(e.b), m.edges)
    all(i -> top(i) == top(1), eachindex(root))
end

# ─── Accessors ──────────────────────────────────────────────────────────────────────────

ngroups(m::Mechanism) = Int(maximum(m.group; init = 0))
members(m::Mechanism, g) = findall(==(g), m.group)
is_allosteric(m::Mechanism) = !isempty(m.tags)
is_reactant(m::Mechanism, x::Symbol) = x in m.substrates || x in m.products
"""Whether `x` is a reactant's modifier copy (a reactant name bound as a modifier)."""
is_dual(m::Mechanism, x::Met) = x.modifier && is_reactant(m, x.name)
is_free(f::Form) = f.bound == 0 && f.res == 0
is_binding(e::Edge) = e.rel == 0 && nibsum(e.cons) == 1
is_plain_binding(m::Mechanism, e::Edge) = is_binding(e) &&
    m.forms[e.b].bound == m.forms[e.a].bound + e.cons &&
    m.forms[e.b].res == m.forms[e.a].res
"""The table entries counted in `p`, each as often as counted."""
mets_of(m::Mechanism, p::Pack) =
    [m.mets[i] for i in eachindex(m.mets) for _ in 1:count_at(p, i)]
"""Names that bind as modifiers, sorted."""
inhibitors(m::Mechanism) = unique!([x.name for x in m.mets if x.modifier])
"""Ligands of the regulatory sites, sorted."""
regulators(m::Mechanism) = sort!(unique!([l for s in m.sites for l in s.ligands]))
"""The concentration names: substrates, products, then the other names sorted."""
metabolites(m::Mechanism) =
    unique!([m.substrates; m.products; sort!(unique!([inhibitors(m); regulators(m)]))])

"""`m` written back as raw input: `(substrates, products, steps, group)` and the keywords
`(; tags, multiplicity, sites)`."""
function raw(m::Mechanism)
    function form(i)
        f = m.forms[i]
        res = [x.name for x in mets_of(m, f.res)]
        RawForm(m.confs[f.conf], mets_of(m, f.bound), filter(in(m.substrates), res),
                filter(in(m.products), res))
    end
    steps = [RawStep(form(e.a), form(e.b), mets_of(m, e.cons), mets_of(m, e.rel), e.re)
             for e in m.edges]
    (m.substrates, m.products, steps, Int.(m.group)),
    (; tags = m.tags, multiplicity = m.multiplicity, sites = m.sites)
end

# ─── Physical-forward orientation ───────────────────────────────────────────────────────

"""Masks of the table entries that are substrates and products (a modifier is neither)."""
function role_masks(m::Mechanism)
    sm = pm = zero(Pack)
    for (i, x) in enumerate(m.mets)
        x.modifier && continue
        nib = Pack(0xf) << (4 * (i - 1))
        x.name in m.substrates ? (sm |= nib) : (pm |= nib)
    end
    sm, pm
end

"""
    orientation(m) -> (fwd, tier)

Whether each edge's `a → b` is its physical-forward direction, by today's tiers, and the
tier that decided it. A binding takes its metabolite up (tier 0). Otherwise the side that
holds more substrates, bound or free, then fewer products comes first (tier 1); otherwise
the step leaves a form where only products enter or leave for one where only substrates
do (tier 2); otherwise the lower form index comes first (tier 3). The members of a kinetic
group that tier 2 or 3 decides turn with the group's lead, the member that changes
conformation first, then by tier and edge: a member taking up and giving off the lead's
metabolites keeps the lead's conformational change; an error when it cannot.
"""
function orientation(m::Mechanism)
    sm, pm = role_masks(m)
    flags(p) = UInt8((p & sm != 0) + 2 * (p & pm != 0))
    ent = zeros(UInt8, length(m.forms))
    for e in m.edges
        ent[e.a] |= flags(e.cons)
        ent[e.b] |= flags(e.rel)
    end
    score(i, free) = (nibsum(m.forms[i].bound & sm) + nibsum(free & sm),
                      -nibsum(m.forms[i].bound & pm) - nibsum(free & pm))
    fwd, tier = trues(length(m.edges)), zeros(Int, length(m.edges))
    for (k, e) in enumerate(m.edges)
        is_binding(e) && continue
        sa, sb = score(e.a, e.cons), score(e.b, e.rel)
        if sa != sb
            fwd[k], tier[k] = sa > sb, 1
        elseif (ent[e.a], ent[e.b]) in ((0x02, 0x01), (0x01, 0x02))
            fwd[k], tier[k] = ent[e.a] == 0x02, 2
        else
            tier[k] = 3
        end
    end
    conf(i) = m.forms[i].conf
    turns(e) = conf(e.a) != conf(e.b)
    ties = Dict{Int32, Vector{Int}}()
    foreach(k -> tier[k] > 1 && push!(get!(ties, m.group[k], Int[]), k), eachindex(tier))
    for tied in values(ties)
        length(tied) < 2 && continue
        lead = oriented(m, argmin(k -> (!turns(oriented(m, k, fwd)), tier[k],
                                        oriented(m, k, fwd)), tied), fwd)
        for k in tied
            o = oriented(m, k, fwd)
            same = [x for x in (o, reversed(o)) if (x.cons, x.rel) == (lead.cons, lead.rel)]
            isempty(same) && continue
            fit = findfirst(x -> !turns(x) || (conf(x.a), conf(x.b)) ==
                                              (conf(lead.a), conf(lead.b)), same)
            fit === nothing &&
                error("a kinetic group holds $(step_text(m, lead)) and " *
                      "$(step_text(m, o)), which make different conformational changes; " *
                      "the steps of a kinetic group share their constants, so they must " *
                      "make the same conformational change")
            fwd[k] = same[fit].a == m.edges[k].a
        end
    end
    fwd, tier
end

orient(m::Mechanism) = first(orientation(m))

"""Edge `k` in its physical-forward direction (`fwd` from `orient`)."""
oriented(m::Mechanism, k, fwd) = fwd[k] ? m.edges[k] : reversed(m.edges[k])

# ─── Rapid-equilibrium segments ─────────────────────────────────────────────────────────

"""
    re_forest(m, mask = all edges) -> (segment, parent, order)

The rapid-equilibrium segments of the edges in `mask`. `segment[i]` numbers form `i`'s
segment in order of its root, the segment's lowest form. `parent[i]` is the edge that
reaches form `i` in a breadth-first tree from the root that scans edges in index order:
`+k` when `i` is `edges[k].b`, `-k` when it is `edges[k].a`, 0 at a root. `order` lists the
forms root first, each after its parent.
"""
function re_forest(m::Mechanism, mask = trues(length(m.edges)))
    n = length(m.forms)
    adj = [Int[] for _ in 1:n]
    for (k, e) in enumerate(m.edges)
        mask[k] && e.re || continue
        push!(adj[e.a], k)
        push!(adj[e.b], -k)
    end
    segment, parent, order = zeros(Int, n), zeros(Int, n), Int[]
    for r in 1:n
        segment[r] == 0 || continue
        segment[r] = isempty(order) ? 1 : segment[order[end]] + 1
        push!(order, r)
        q = length(order)
        while q <= length(order)
            u = order[q]
            q += 1
            for k in adj[u]
                v = k > 0 ? m.edges[k].b : m.edges[-k].a
                segment[v] == 0 || continue
                segment[v], parent[v] = segment[r], k
                push!(order, v)
            end
        end
    end
    segment, parent, order
end

# ─── Boundary checks (DSL and enumeration), separate from the constructor ───────────────

"""
    bottomless(m, mask = all edges) -> forms or nothing

The forms of a rapid-equilibrium segment whose weights all vanish when one side's
metabolites are absent (products first, then substrates), so its rate is 0/0 there;
`nothing` when no segment does. A modifier copy counts on its reactant's side.
"""
function bottomless(m::Mechanism, mask = trues(length(m.edges)))
    segment, parent, order = re_forest(m, mask)
    names = unique!([x.name for x in m.mets])
    row = [searchsortedfirst(names, x.name) for x in m.mets]
    expo = zeros(Int, length(names), length(m.forms))
    for i in order
        k = parent[i]
        k == 0 && continue
        e = m.edges[abs(k)]
        expo[:, i] .= view(expo, :, k > 0 ? e.a : e.b)
        for j in eachindex(m.mets)
            expo[row[j], i] += sign(k) * (count_at(e.cons, j) - count_at(e.rel, j))
        end
    end
    side = [x in m.substrates ? 1 : x in m.products ? 2 : 0 for x in names]
    for s in 1:maximum(segment; init = 0)
        fs = findall(==(s), segment)
        length(fs) < 2 && continue
        low = minimum(expo[:, fs]; dims = 2)
        for sd in (2, 1)
            all(i -> any(r -> side[r] == sd && expo[r, i] > low[r], eachindex(names)),
                fs) &&
                return fs
        end
    end
    nothing
end

"""The first binding edge whose bound form holds its free form's metabolites plus the
bound one by name but not by role (`E + P::Mod → E(P)`); `nothing` when none does."""
function role_change(m::Mechanism)
    byname(p) = sort!([x.name for x in mets_of(m, p)])
    for e in m.edges
        is_binding(e) && !is_plain_binding(m, e) || continue
        fa, fb = m.forms[e.a], m.forms[e.b]
        fa.res == fb.res && byname(fb.bound) == byname(fa.bound + e.cons) && return e
    end
    nothing
end

"""
    check_atoms(m, atoms)

Error unless every step conserves atoms; `atoms[x]` lists the `element => count` pairs of
reactant `x` (a name without an entry carries none). A form holds its bound metabolites'
atoms plus its residual's: added substrates' minus subtracted products'.
"""
function check_atoms(m::Mechanism, atoms)
    function tally!(out, p, sign, res = false)
        for (i, x) in enumerate(m.mets), (el, n) in get(atoms, x.name, ())
            s = res && x.name in m.products ? -sign : sign
            out[el] = get(out, el, 0) + s * count_at(p, i) * n
        end
        out
    end
    side(i, p) = tally!(tally!(tally!(Dict{Symbol, Int}(), m.forms[i].bound, 1),
                               m.forms[i].res, 1, true), p, 1)
    for e in m.edges
        a, b = side(e.a, e.cons), side(e.b, e.rel)
        all(el -> get(a, el, 0) == get(b, el, 0), union(keys(a), keys(b))) ||
            error("$(step_text(m, e)) does not conserve atoms")
    end
end

"""Error on a bottomless rapid-equilibrium segment, a binding that changes a metabolite's
role, or colliding rendered names (`check_names`)."""
function check_boundary(m::Mechanism)
    fs = bottomless(m)
    fs === nothing ||
        error("rapid-equilibrium segment {" * join((form_text(m, i) for i in fs), ", ") *
              "} has no form free of one side's metabolites, so its rate is undefined " *
              "when that side is absent (0/0); make one of its bindings steady-state")
    e = role_change(m)
    e === nothing || error("$(step_text(m, e)) changes a metabolite's role; write a " *
                           "reactant that binds as a modifier as X::Mod on both sides")
    check_names(m)
end

# ─── The Sig: a mechanism as one homogeneous integer tuple ──────────────────────────────

"""A mechanism as a singleton type; its Sig (`sig`) is the only type parameter."""
struct EnzymeMechanism{S} end

"""
    sig(m) -> NTuple{N, Int}

`m` as one homogeneous integer tuple: the counts of substrates, products, table entries,
conformations, forms, edges, tags and sites, then the multiplicity; every name as its
length and code units; each table entry's modifier flag; each form as `bound` and `res`
in 64-bit halves and its conformation; each edge as `a`, `b`, the halves of `cons` and
`rel`, and its flag; the group labels; the tag codes; each site as its multiplicity, its
ligand count and each ligand's name and state code. Equal mechanisms give equal Sigs.
"""
function sig(m::Mechanism)
    out = Int[length(m.substrates), length(m.products), length(m.mets), length(m.confs),
              length(m.forms), length(m.edges), length(m.tags), length(m.sites),
              m.multiplicity]
    name!(s) = (u = codeunits(String(s)); push!(out, length(u)); append!(out, u))
    halves!(p) = push!(out, (p >> 64) % Int, p % Int)
    foreach(name!, m.substrates)
    foreach(name!, m.products)
    foreach(x -> (name!(x.name); push!(out, x.modifier)), m.mets)
    foreach(name!, m.confs)
    foreach(f -> (halves!(f.bound); halves!(f.res); push!(out, f.conf)), m.forms)
    for e in m.edges
        push!(out, e.a, e.b)
        halves!(e.cons)
        halves!(e.rel)
        push!(out, e.re)
    end
    append!(out, m.group)
    foreach(t -> push!(out, findfirst(==(t), TAGS)), m.tags)
    for s in m.sites
        push!(out, s.multiplicity, length(s.ligands))
        for (l, st) in zip(s.ligands, s.states)
            name!(l)
            push!(out, findfirst(==(st), SITE_STATES))
        end
    end
    Tuple(out)
end

"""
    decode(S...) -> Mechanism

The mechanism whose Sig is `S`, splatted: the varargs method, and the vector it reads the
Sig through, compile once for every length.
"""
function decode(xs::Int...)
    v, at = Int[xs...], Ref(0)
    next() = v[at[] += 1]
    name() = Symbol(String(UInt8[next() for _ in 1:next()]))
    halves() = Pack(next() % UInt64) << 64 | Pack(next() % UInt64)
    ns, np, nm, nc, nf, ne, nt, nsites, multiplicity = (next() for _ in 1:9)
    subs, prods = [name() for _ in 1:ns], [name() for _ in 1:np]
    mets = [Met(name(), next() == 1) for _ in 1:nm]
    confs = [name() for _ in 1:nc]
    forms = [Form(halves(), halves(), next()) for _ in 1:nf]
    edges = [Edge(next(), next(), halves(), halves(), next() == 1) for _ in 1:ne]
    group = Int32[next() for _ in 1:ne]
    tags = Symbol[TAGS[next()] for _ in 1:nt]
    function site()
        mult, nl = next(), next()
        ls = [(name(), SITE_STATES[next()]) for _ in 1:nl]
        Site(first.(ls), last.(ls), mult)
    end
    sites = Site[site() for _ in 1:nsites]
    Mechanism(subs, prods, mets, confs, forms, edges, group, tags, multiplicity, sites)
end

"""The singleton of `m`'s compiled type."""
compile(m::Mechanism) = EnzymeMechanism{sig(m)}()

"""The mechanism of a compiled singleton; reads the Sig at run time, so it compiles once."""
mechanism(@nospecialize(em::EnzymeMechanism)) = decode(typeof(em).parameters[1]...)
