# ABOUTME: Converts the probes' Pop.PMech records to mechanisms, reading only their fields,
# ABOUTME: and defines the record type of the converted populations later stages load.

"""
A converted Pop record: its `id`, `source` and `index` in its population file, the
mechanism, and for each Pop step `j` the edge `edge_of[j]` it became and whether that edge
runs against the step's stored direction (`flipped[j]`).
"""
struct Converted
    id::String
    source::String
    index::Int
    mech::Mechanism
    edge_of::Vector{Int}
    flipped::BitVector
end

"""The raw input of a `Pop.PMech` record `p`: substrates, products, steps, group labels,
tags, multiplicity and sites."""
function pop_raw(p)
    met(x) = x.role in (:S, :P, :Inh) ? Met(x.name, x.role === :Inh) :
             error("$(x.name) binds in a step as $(x.role)")
    function form(i)
        sp = p.species[i]
        RawForm(sp.conformation, Met[met(x) for x in sp.bound], sp.res_added,
                sp.res_subtracted)
    end
    steps = [RawStep(form(s.from), form(s.to), Met[met(x) for x in s.consumed],
                     Met[met(x) for x in s.released], s.re) for s in p.steps]
    group = zeros(Int, length(steps))
    for (g, js) in enumerate(p.groups), j in js
        group[j] = g
    end
    (; p.substrates, p.products, steps, group, tags = copy(p.group_tags),
     multiplicity = p.cat_multiplicity,
     sites = [Site(s.ligands, s.states, s.multiplicity) for s in p.sites])
end

"""`canonicalize` of a raw input as `pop_raw` returns it."""
canonicalize(r::NamedTuple) = canonicalize(r.substrates, r.products, r.steps, r.group;
                                           r.tags, r.multiplicity, r.sites)

"""The `Converted` record of the `index`-th record `p` of a population."""
Converted(p, index::Int) = Converted(p.id, p.source, index, canonicalize(pop_raw(p))...)

"""
The map between edge `k`'s constants (physical-forward, in `constant_names` order) and Pop
step `j`'s (its stored direction, the oracle's convention), with `fwd = orient(c.mech)`:
where the step runs against the physical-forward direction, `[kf, kr]` swap and a
non-binding `K` inverts (a binding's dissociation constant stays). Its own inverse.
"""
function pop_flip(c::Converted, fwd, j, x)
    k = c.edge_of[j]
    fwd[k] != c.flipped[j] || (length(x) == 1 && is_binding(c.mech.edges[k])) ? x :
        length(x) == 2 ? reverse(x) : inv.(x)
end

"""Per Pop step constants in the oracle's convention from per-edge constants `x`."""
function pop_constants(c::Converted, x)
    fwd = orient(c.mech)
    [pop_flip(c, fwd, j, x[c.edge_of[j]]) for j in eachindex(c.edge_of)]
end

"""Per-edge constants from per Pop step constants `y` in the oracle's convention."""
function edge_constants(c::Converted, y)
    fwd, x = orient(c.mech), similar(y)
    foreach(j -> x[c.edge_of[j]] = pop_flip(c, fwd, j, y[j]), eachindex(y))
    x
end
