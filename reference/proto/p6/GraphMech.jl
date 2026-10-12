# ABOUTME: Prototype graph-based mechanism type (module GM, probe P6): content-sorted forms,
# ABOUTME: isbits edges, per-edge kinetic-group labels and lazily cached graph analyses; stdlib only.
module GM

using Serialization

# ─── Metabolites and the raw input a DSL or an enumeration move writes ──────────────

"""A metabolite as a form holds it or a step exchanges it: its name, and `inh` when it
binds as a competitive inhibitor (a standalone inhibitor or a reactant's `X::Inh` copy)."""
struct Met
    name::Symbol
    inh::Bool
end
Base.isless(x::Met, y::Met) = isless((x.name, x.inh), (y.name, y.inh))
Base.hash(x::Met, h::UInt) = hash(x.inh, hash(x.name, h))

"""An enzyme form as written: conformation, bound metabolites and residual (substrates
added, products subtracted), each list in any order."""
struct RawForm
    conf::Symbol
    bound::Vector{Met}
    added::Vector{Symbol}
    subtracted::Vector{Symbol}
end

"""A step as written, in either direction: `from` takes up `consumed` and gives off
`released` to become `to`."""
struct RawStep
    from::RawForm
    to::RawForm
    consumed::Vector{Met}
    released::Vector{Met}
    re::Bool
end

"""An allosteric regulatory site: ligands with their states, and its multiplicity."""
struct Site
    ligands::Vector{Symbol}
    states::Vector{Symbol}
    multiplicity::Int
end
Base.isless(s::Site, t::Site) =
    isless((s.ligands, s.multiplicity, s.states), (t.ligands, t.multiplicity, t.states))

const TAGS = (:OnlyA, :EqualAI, :NonequalAI)
const SITE_STATES = (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)

# ─── Packed metabolite counts ────────────────────────────────────────────────────────

"""Counts of each metabolite-table entry, four bits per entry (at most 32 entries and 15
copies)."""
const Pack = UInt128
const LOW = Pack(0x11111111111111111111111111111111)

count_at(p::Pack, i) = Int((p >> (4 * (i - 1))) & 0xf)
nibsum(p::Pack) = (s = 0; while p != 0; s += Int(p & 0xf); p >>= 4; end; s)
overlap(p::Pack, q::Pack) = ((p | p >> 1 | p >> 2 | p >> 3) & (q | q >> 1 | q >> 2 | q >> 3) & LOW) != 0

function pack(mets::Vector{Met}, xs)
    p, n = zero(Pack), 0
    for x in xs
        p += one(Pack) << (4 * (findfirst(==(x), mets) - 1))
        n += 1
    end
    nibsum(p) == n || error("a metabolite occurs more than 15 times in one list")
    p
end

# ─── The mechanism ──────────────────────────────────────────────────────────────────

"""An enzyme form: counts of each metabolite-table entry it holds (`bound`) and of its
residual (`res`: a substrate added or a product subtracted), and its conformation (index
into `confs`). Forms sort by this content."""
struct Form
    bound::Pack
    res::Pack
    conf::Int32
end
Base.isless(f::Form, g::Form) = isless((f.bound, f.res, f.conf), (g.bound, g.res, g.conf))
Base.hash(f::Form, h::UInt) = hash(f.conf, hash(f.res, hash(f.bound, h)))

"""A step between forms `a` and `b`: going `a → b` it takes up `cons` and gives off `rel`.
Stored with `a < b`, except a binding-like step (one metabolite in all, on one side), which
is stored taking its metabolite up."""
struct Edge
    a::Int32
    b::Int32
    cons::Pack
    rel::Pack
    re::Bool
end
Base.isless(e::Edge, f::Edge) =
    isless((e.a, e.b, e.cons, e.rel, e.re), (f.a, f.b, f.cons, f.rel, f.re))
Base.hash(e::Edge, h::UInt) = hash(e.re, hash(e.rel, hash(e.cons, hash(e.b, hash(e.a, h)))))
reversed(e::Edge) = Edge(e.b, e.a, e.rel, e.cons, e.re)

"""Graph analyses of a mechanism or of an edge mask over it: the rapid-equilibrium segment
of each form (numbered by its lowest form) with each form's metabolite exponents above the
segment's lowest (`extras[name, form]`), the spanning-forest edges, the fundamental cycles
(`+k` runs edge `k` from `a` to `b`, `-k` back) and each cycle's multiple of the net
reaction."""
struct Analysis
    segment::Vector{Int32}
    names::Vector{Symbol}
    extras::Matrix{Int}
    tree::BitVector
    cycles::Vector{Vector{Int32}}
    haldane::Vector{Int}
end

"""
A mechanism as a graph. Forms are sorted by content and edges by `(a, b, cons, rel, re)`;
`group[k]` is the kinetic group of edge `k`, numbered by each group's first edge, and
`tags[g]` its allosteric state (empty for a plain mechanism). `mets` and `confs` are the
metabolite and conformation tables the forms and edges index. `cache` holds the analyses
of the whole graph once computed; it takes no part in `==`, `hash` or serialization.
"""
struct GMech
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
    cache::Base.RefValue{Union{Nothing, Analysis}}
end

const FIELDS = (:substrates, :products, :mets, :confs, :forms, :edges, :group, :tags,
                :multiplicity, :sites)

for (T, fs) in ((Site, fieldnames(Site)), (GMech, FIELDS))
    @eval Base.:(==)(x::$T, y::$T) =
        $(foldl((l, f) -> :($l && x.$f == y.$f), fs[2:end]; init = :(x.$(fs[1]) == y.$(fs[1]))))
    @eval Base.hash(x::$T, h::UInt) =
        $(foldl((acc, f) -> :(hash(x.$f, $acc)), fs; init = :(hash($(QuoteNode(nameof(T))), h))))
end

function Serialization.serialize(s::AbstractSerializer, m::GMech)
    Serialization.serialize_type(s, GMech, false)
    for f in FIELDS; serialize(s, getfield(m, f)); end
end
Serialization.deserialize(s::AbstractSerializer, ::Type{GMech}) =
    GMech((deserialize(s) for _ in FIELDS)..., Ref{Union{Nothing, Analysis}}(nothing))

addunique!(v, x) = (x in v || push!(v, x); v)

function formkey(mets, confs, subs, prods, f::RawForm)
    all(in(subs), f.added) && all(in(prods), f.subtracted) ||
        error("a residual adds a substrate and subtracts a product")
    res = pack(mets, (Met(x, false) for x in f.added)) +
          pack(mets, (Met(x, false) for x in f.subtracted))
    Form(pack(mets, f.bound), res, Int32(findfirst(==(f.conf), confs)))
end

"""Relabel `labels` (one per edge, already in edge order) by first occurrence, and carry
`tags` along; errors on a label no edge holds."""
function relabel(labels, tags::Vector{Symbol}, ngroups::Int)
    new = zeros(Int32, ngroups)
    next = Int32(0)
    group = Vector{Int32}(undef, length(labels))
    for (k, l) in enumerate(labels)
        new[l] == 0 && (new[l] = (next += Int32(1)))
        group[k] = new[l]
    end
    next == ngroups || error("kinetic group $(findfirst(iszero, new)) holds no step")
    t = similar(tags)
    for l in eachindex(tags); t[new[l]] = tags[l]; end
    group, t
end

"""
    canonicalize(substrates, products, steps, group; tags, multiplicity, sites)
        -> (m, edge_of, flipped)

The mechanism of `steps`, written in any order and direction, with kinetic group label
`group[j] ∈ 1:G` for step `j` and, for an allosteric mechanism, `tags[g]` per label.
`edge_of[j]` is the edge step `j` became and `flipped[j]` whether it was reversed.
"""
function canonicalize(substrates::Vector{Symbol}, products::Vector{Symbol},
                      steps::Vector{RawStep}, group::Vector{Int};
                      tags::Vector{Symbol} = Symbol[], multiplicity::Int = 1,
                      sites::Vector{Site} = Site[])
    subs, prods = sort(substrates), sort(products)
    allunique(subs) && allunique(prods) && isdisjoint(subs, prods) ||
        error("substrates and products must be distinct names")
    mets, confs = Met[], Symbol[]
    for s in steps
        for f in (s.from, s.to)
            foreach(x -> addunique!(mets, x), f.bound)
            foreach(x -> addunique!(mets, Met(x, false)), f.added)
            foreach(x -> addunique!(mets, Met(x, false)), f.subtracted)
            addunique!(confs, f.conf)
        end
        foreach(x -> addunique!(mets, x), s.consumed)
        foreach(x -> addunique!(mets, x), s.released)
    end
    sort!(mets); sort!(confs)
    length(mets) <= 32 || error("more than 32 metabolite entries")
    n = length(steps)
    keys = Vector{Form}(undef, 2n)
    for (j, s) in enumerate(steps)
        keys[2j - 1] = formkey(mets, confs, subs, prods, s.from)
        keys[2j] = formkey(mets, confs, subs, prods, s.to)
    end
    forms = unique!(sort(keys))
    edges = Vector{Edge}(undef, n)
    flipped = falses(n)
    for (j, s) in enumerate(steps)
        a = Int32(searchsortedfirst(forms, keys[2j - 1]))
        b = Int32(searchsortedfirst(forms, keys[2j]))
        a == b && error("step $j joins a form to itself")
        c, r = pack(mets, s.consumed), pack(mets, s.released)
        overlap(c, r) && error("step $j takes up and gives off one metabolite")
        flip = nibsum(c) + nibsum(r) == 1 ? c == 0 : a > b
        edges[j] = flip ? Edge(b, a, r, c, s.re) : Edge(a, b, c, r, s.re)
        flipped[j] = flip
    end
    p = sortperm(edges)
    edges = edges[p]
    for k in 2:n
        e, f = edges[k - 1], edges[k]
        (e.a, e.b, e.cons, e.rel) == (f.a, f.b, f.cons, f.rel) &&
            error("one reaction appears in two steps")
    end
    G = maximum(group; init = 0)
    isempty(tags) || length(tags) == G || error("one tag per kinetic group")
    all(in(TAGS), tags) || error("a catalytic tag must be one of $TAGS")
    g, t = relabel(view(group, p), tags, G)
    lead = zeros(Int, G)
    for (k, e) in enumerate(edges)
        lead[g[k]] == 0 && (lead[g[k]] = k; continue)
        l = edges[lead[g[k]]]
        e.re == l.re && ((e.cons, e.rel) == (l.cons, l.rel) || (e.rel, e.cons) == (l.cons, l.rel)) ||
            error("kinetic group $(g[k]) mixes kinds of step or RE and SS")
    end
    sites = sort!(map(canonical_site, sites))
    allunique(l for s in sites for l in s.ligands) || error("a ligand sits at two sites")
    m = GMech(subs, prods, mets, confs, forms, edges, g, t, multiplicity, sites,
              Ref{Union{Nothing, Analysis}}(nothing))
    m, invperm(p), flipped
end

GMech(substrates, products, steps, group; kw...) =
    first(canonicalize(substrates, products, steps, group; kw...))

function canonical_site(s::Site)
    s.multiplicity >= 1 && !isempty(s.ligands) && length(s.ligands) == length(s.states) &&
        all(in(SITE_STATES), s.states) || error("malformed regulatory site")
    p = sortperm(s.ligands)
    Site(s.ligands[p], s.states[p], s.multiplicity)
end

# ─── Accessors ──────────────────────────────────────────────────────────────────────

mets_of(m::GMech, p::Pack) = [m.mets[i] for i in eachindex(m.mets) for _ in 1:count_at(p, i)]
is_free(f::Form) = f.bound == 0 && f.res == 0
is_binding(e::Edge) = e.rel == 0 && nibsum(e.cons) == 1
is_plain_binding(m::GMech, e::Edge) = is_binding(e) &&
    m.forms[e.b].bound == m.forms[e.a].bound + e.cons && m.forms[e.b].res == m.forms[e.a].res
ngroups(m::GMech) = Int(maximum(m.group; init = 0))
inhibitors(m::GMech) = unique!([x.name for x in m.mets if x.inh])
regulators(m::GMech) = unique!([l for s in m.sites for l in s.ligands])
metabolites(m::GMech) = unique!([m.substrates; m.products; inhibitors(m); regulators(m)])

# ─── Children for enumeration ───────────────────────────────────────────────────────

"""`m` written back as raw input: `(substrates, products, steps, group, kw)`."""
function raw(m::GMech)
    function form(i)
        f = m.forms[i]
        res = mets_of(m, f.res)
        RawForm(m.confs[f.conf], mets_of(m, f.bound), [x.name for x in res if x.name in m.substrates],
                [x.name for x in res if x.name in m.products])
    end
    steps = [RawStep(form(e.a), form(e.b), mets_of(m, e.cons), mets_of(m, e.rel), e.re) for e in m.edges]
    (m.substrates, m.products, steps, Int.(m.group),
     (; tags = m.tags, multiplicity = m.multiplicity, sites = m.sites))
end

"""`m` with every edge of kinetic group `g` rapid-equilibrium (`re`) or steady-state. Edge
order does not depend on the flag, so nothing is re-sorted."""
function with_re(m::GMech, g::Integer, re::Bool)
    edges = [m.group[k] == g ? Edge(e.a, e.b, e.cons, e.rel, re) : e for (k, e) in enumerate(m.edges)]
    GMech(m.substrates, m.products, m.mets, m.confs, m.forms, edges, m.group, m.tags,
          m.multiplicity, m.sites, Ref{Union{Nothing, Analysis}}(nothing))
end

"""`m` with edges `ks`, part of one kinetic group, moved to a group of their own that keeps
the group's tag."""
function with_split(m::GMech, ks)
    labels = copy(m.group)
    labels[ks] .= ngroups(m) + 1
    tags = isempty(m.tags) ? m.tags : [m.tags; m.tags[m.group[first(ks)]]]
    g, t = relabel(labels, tags, ngroups(m) + 1)
    GMech(m.substrates, m.products, m.mets, m.confs, m.forms, m.edges, g, t,
          m.multiplicity, m.sites, Ref{Union{Nothing, Analysis}}(nothing))
end

"""Masks of the table entries that are substrates and products (an inhibitor copy is
neither)."""
function role_masks(m::GMech)
    sm = pm = zero(Pack)
    for (i, x) in enumerate(m.mets)
        x.inh && continue
        x.name in m.substrates && (sm |= Pack(0xf) << (4 * (i - 1)))
        x.name in m.products && (pm |= Pack(0xf) << (4 * (i - 1)))
    end
    sm, pm
end

"""The package's rendered name of form `i` (`EAB`, `EPinh`, `E_res_+A_-P`)."""
function form_name(m::GMech, i)
    f = m.forms[i]
    s = string(m.confs[f.conf], (x.inh ? string(x.name, "inh") : x.name for x in mets_of(m, f.bound))...)
    f.res == 0 && return s
    res = mets_of(m, f.res)
    string(s, "_res", ("_+$(x.name)" for x in res if x.name in m.substrates)...,
           ("_-$(x.name)" for x in res if x.name in m.products)...)
end

"""Form `i` as the DSL writes it: `E`, `E(A, P::Inh)`, `E(B; residual = A - P)`."""
function form_text(m::GMech, i)
    f = m.forms[i]
    ms = [x.inh ? "$(x.name)::Inh" : "$(x.name)" for x in mets_of(m, f.bound)]
    res = mets_of(m, f.res)
    r = isempty(res) ? "" : "; residual = " *
        join([x.name for x in res if x.name in m.substrates], " + ") *
        join(" - $(x.name)" for x in res if x.name in m.products)
    conf = string(m.confs[f.conf])
    isempty(ms) && isempty(r) ? conf : "$conf($(join(ms, ", "))$r)"
end

function step_text(m::GMech, e::Edge)
    side(i, p) = join([form_text(m, i); [x.inh ? "$(x.name)::Inh" : "$(x.name)" for x in mets_of(m, p)]], " + ")
    string(side(e.a, e.cons), e.re ? " ⇌ " : " <--> ", side(e.b, e.rel))
end

# ─── Physical-forward orientation (naming and priority) ────────────────────────────

"""
    orient(m; tier3 = :name) -> (fwd, tier)

Whether each edge's `a → b` is its physical-forward direction, by the package's tiers: a
binding takes its metabolite up (tier 0); otherwise substrates held or free, then products
held or free, decrease (tier 1); otherwise the step leaves a form where only products enter
or leave for one where only substrates do (tier 2); otherwise the smaller rendered form
name comes first (`tier3 = :name`, the package's rule) or the smaller form index
(`tier3 = :index`). The steps of a kinetic group that tier 1 leaves tied turn after the
group's lead: the same metabolites taken up and given off, then the same conformation
change. `tier[k]` is the tier that decided edge `k`.
"""
function orient(m::GMech; tier3::Symbol = :name)
    sm, pm = role_masks(m)
    names = tier3 === :name ? [form_name(m, i) for i in eachindex(m.forms)] : String[]
    flags(p) = UInt8((p & sm != 0) + 2 * (p & pm != 0))
    ent = zeros(UInt8, length(m.forms))
    for e in m.edges
        ent[e.a] |= flags(e.cons); ent[e.b] |= flags(e.rel)
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
            fwd[k], tier[k] = tier3 === :name ? names[e.a] <= names[e.b] : true, 3
        end
    end
    conf(i) = m.forms[i].conf
    turns(e) = conf(e.a) != conf(e.b)
    oriented(k) = fwd[k] ? m.edges[k] : reversed(m.edges[k])
    metnames(p) = join((x.name for x in mets_of(m, p)), "+")
    namekey(e) = (names[e.a], names[e.b], metnames(e.cons), metnames(e.rel), e.re)
    for g in 1:ngroups(m)
        tied = [k for k in eachindex(m.edges) if m.group[k] == g && tier[k] > 1]
        length(tied) < 2 && continue
        lead = oriented(tier3 === :name ?
            argmin(k -> (!turns(oriented(k)), tier[k], namekey(oriented(k))), tied) :
            argmin(k -> (!turns(oriented(k)), tier[k], oriented(k)), tied))
        for k in tied
            same = [o for o in (oriented(k), reversed(oriented(k))) if (o.cons, o.rel) == (lead.cons, lead.rel)]
            isempty(same) && continue
            fit = findfirst(o -> !turns(o) || (conf(o.a), conf(o.b)) == (conf(lead.a), conf(lead.b)), same)
            fit === nothing && error("kinetic group $g makes two different conformation changes")
            fwd[k] = same[fit].a == m.edges[k].a
        end
    end
    fwd, tier
end

physical_forward(m::GMech; tier3::Symbol = :name) = first(orient(m; tier3))

"""`m` with every edge stored in its physical-forward direction (orientation rule 2)."""
function physical(m::GMech; tier3::Symbol = :name)
    fwd = physical_forward(m; tier3)
    edges = [f ? e : reversed(e) for (e, f) in zip(m.edges, fwd)]
    p = sortperm(edges)
    g, t = relabel(view(m.group, p), m.tags, ngroups(m))
    GMech(m.substrates, m.products, m.mets, m.confs, m.forms, edges[p], g, t, m.multiplicity,
          m.sites, Ref{Union{Nothing, Analysis}}(nothing))
end

# ─── Graph analyses on an edge mask ─────────────────────────────────────────────────

"""Adjacency of the masked edges: signed edge ids, `+k` leaving `a`, `-k` leaving `b`."""
function adjacency(m::GMech, mask::BitVector, only_re::Bool)
    adj = [Int[] for _ in m.forms]
    for (k, e) in enumerate(m.edges)
        mask[k] && (e.re || !only_re) || continue
        push!(adj[e.a], k); push!(adj[e.b], -k)
    end
    adj
end

"""The analyses of the edges in `mask` (all edges by default)."""
function analyze(m::GMech, mask::BitVector = trues(length(m.edges)))
    n = length(m.forms)
    names = unique!([x.name for x in m.mets])
    row = [searchsortedfirst(names, x.name) for x in m.mets]
    # Rapid-equilibrium segments, breadth first from each segment's lowest form.
    adj = adjacency(m, mask, true)
    segment = zeros(Int32, n)
    expo = zeros(Int, length(names), n)
    queue = Int[]
    nseg = 0
    for r in 1:n
        segment[r] == 0 || continue
        segment[r] = (nseg += 1)
        push!(queue, r)
        while !isempty(queue)
            u = popfirst!(queue)
            for sk in adj[u]
                e = m.edges[abs(sk)]
                v = sk > 0 ? e.b : e.a
                segment[v] == 0 || continue
                segment[v] = nseg
                expo[:, v] .= view(expo, :, u)
                for i in eachindex(m.mets)
                    expo[row[i], v] += sign(sk) * (count_at(e.cons, i) - count_at(e.rel, i))
                end
                push!(queue, v)
            end
        end
    end
    for s in 1:nseg, x in eachindex(names)
        low = minimum(expo[x, i] for i in 1:n if segment[i] == s)
        for i in 1:n; segment[i] == s && (expo[x, i] -= low); end
    end
    tree, cycles = cycle_basis(m, mask)
    Analysis(segment, names, expo, tree, cycles, [haldane_multiple(m, c) for c in cycles])
end

"""Spanning forest (breadth first from the lowest form) and fundamental cycles of the
masked edges."""
function cycle_basis(m::GMech, mask::BitVector)
    n = length(m.forms)
    adj = adjacency(m, mask, false)
    parent, depth = zeros(Int, n), fill(-1, n)
    tree = falses(length(m.edges))
    queue = Int[]
    for r in 1:n
        depth[r] < 0 || continue
        depth[r] = 0
        push!(queue, r)
        while !isempty(queue)
            u = popfirst!(queue)
            for sk in adj[u]
                e = m.edges[abs(sk)]
                v = sk > 0 ? e.b : e.a
                depth[v] < 0 || continue
                depth[v], parent[v], tree[abs(sk)] = depth[u] + 1, sk, true
                push!(queue, v)
            end
        end
    end
    up(v) = (e = m.edges[abs(parent[v])]; parent[v] > 0 ? e.a : e.b)
    cycles = Vector{Int32}[]
    for (k, e) in enumerate(m.edges)
        mask[k] && !tree[k] || continue
        from_b, from_a = Int32[k], Int32[]
        u, v = Int(e.b), Int(e.a)
        while u != v
            if depth[u] >= depth[v]
                push!(from_b, -parent[u]); u = up(u)
            else
                push!(from_a, parent[v]); v = up(v)
            end
        end
        push!(cycles, append!(from_b, reverse!(from_a)))
    end
    tree, cycles
end

"""The net reaction's multiple that cycle `c` runs; errors if the cycle runs anything else."""
function haldane_multiple(m::GMech, c::Vector{Int32})
    nu = Dict{Symbol, Int}()
    for sk in c, (i, x) in enumerate(m.mets)
        e = m.edges[abs(sk)]
        d = sign(sk) * (count_at(e.rel, i) - count_at(e.cons, i))
        d == 0 || (nu[x.name] = get(nu, x.name, 0) + d)
    end
    mult = -get(nu, m.substrates[1], 0)
    all(get(nu, x, 0) == -mult for x in m.substrates) &&
        all(get(nu, x, 0) == mult for x in m.products) ||
        error("a cycle runs a reaction other than a multiple of the net reaction")
    mult
end

"""The analyses of the whole graph, computed on first use and cached."""
function analysis(m::GMech)
    a = m.cache[]
    a === nothing || return a
    m.cache[] = analyze(m)
end

"""Forms of a rapid-equilibrium segment whose weights all vanish when one side's
metabolites are absent; `nothing` when no segment does."""
function bottomless(m::GMech, a::Analysis)
    side = [x in m.substrates ? 1 : x in m.products ? 2 : 0 for x in a.names]
    for s in 1:maximum(a.segment; init = 0)
        members = findall(==(s), a.segment)
        length(members) < 2 && continue
        for sd in (2, 1)
            all(i -> any(r -> side[r] == sd && a.extras[r, i] > 0, eachindex(side)), members) &&
                return members
        end
    end
    nothing
end

"""The inactive state as an edge mask: the edges of non-`:OnlyA` groups whose forms the
free enzyme reaches over those edges."""
function inactive_mask(m::GMech)
    live = BitVector([m.tags[m.group[k]] !== :OnlyA for k in eachindex(m.edges)])
    reach = falses(length(m.forms))
    for (k, e) in enumerate(m.edges)
        live[k] || continue
        is_free(m.forms[e.a]) && (reach[e.a] = true)
        is_free(m.forms[e.b]) && (reach[e.b] = true)
    end
    changed = true
    while changed
        changed = false
        for (k, e) in enumerate(m.edges)
            live[k] && reach[e.a] != reach[e.b] || continue
            reach[e.a] = reach[e.b] = changed = true
        end
    end
    BitVector([live[k] && reach[e.a] && reach[e.b] for (k, e) in enumerate(m.edges)])
end

end # module
