# ABOUTME: A mechanism as raw input (forms, steps, the whole), the canonical graph of raw
# ABOUTME: input, and the writers that rebuild a mechanism with new flags, groups or tags.

# ─── The canonical graph ────────────────────────────────────────────────

"""
    canonical_graph(raw::RawMechanism) -> Mechanism

The canonical graph `raw` writes (`_canonicalize`), with the representation checks of
`Mechanism(raw)` left out but one: a kinetic group whose tied steps no orientation aligns
has no canonical graph, so it raises the representation error (`_check_aligned`). The
seed variants read their merged and Theorell–Chance bases this way, and the dead-end move
reads a candidate before it decides to build it.
"""
function canonical_graph(raw::RawMechanism)
    c = _canonicalize(raw)
    _check_aligned(c.mechanism, c.misaligned)
    c.mechanism
end

# ─── The written form ───────────────────────────────────────────────────

"""
    RawForm(m::Mechanism, x::Int, extra::Met...)

Form `x` of `m` as written: its conformation label, its bound metabolites with `extra`
added, and its residual.
"""
function RawForm(m::Mechanism, x::Int, extra::Met...)
    n = m.nodes[x]
    RawForm(n.conformation, Met[n.bound; extra...], n.added, n.subtracted)
end

"""
    RawStep(m::Mechanism, k::Int, extra::Met...)

Step `k` of `m` as written, in its stored direction and at its own flag, each of its two
forms also holding `extra` (`RawForm`). With a modifier as `extra` it is the step's
mirror between the complexes the modifier forms at its two ends.
"""
RawStep(m::Mechanism, k::Int, extra::Met...) =
    _written_step(m, k, m.steps[k].rapid_equilibrium, extra...)

"""Step `k` of `m` as written (`RawStep`), at rapid equilibrium when `re` holds."""
function _written_step(m::Mechanism, k::Int, re::Bool, extra::Met...)
    e = m.steps[k]
    RawStep(RawForm(m, e.from, extra...), RawForm(m, e.to, extra...), e.taken_up,
            e.given_off, re)
end

"""
    RawMechanism(m::Mechanism)

`m` as written: its declared names, and each step as `RawStep(m, k)`, labelled by its
kinetic group. Every vector is fresh, and `Mechanism(RawMechanism(m)) == m`.
"""
RawMechanism(m::Mechanism) =
    _written(m, BitVector(e.rapid_equilibrium for e in m.steps), _group_labels(m))

"""
    RawAllostericMechanism(am::AllostericMechanism)

`am` as written: its catalytic graph as `RawMechanism` writes it, each group label tagged
as its group, and its multiplicity and regulatory sites. Every vector is fresh, and
`AllostericMechanism(RawAllostericMechanism(am)) == am`.
"""
RawAllostericMechanism(am::AllostericMechanism) =
    _written(am, BitVector(e.rapid_equilibrium for e in am.catalytic.steps),
             _group_labels(am.catalytic))

"""`m` as written (`RawMechanism`), step `k` at rapid equilibrium where `re[k]` holds and
labelled `groups[k]`."""
_written(m::Mechanism, re::BitVector, groups::Vector{Int}) =
    RawMechanism(copy(m.substrates), copy(m.products), modifiers(m),
                 RawStep[_written_step(m, k, re[k]) for k in eachindex(m.steps)],
                 copy(groups))

"""`am` as written (`RawAllostericMechanism`), its catalytic graph as `_written` writes it
and each label tagged as the group of its first step."""
_written(am::AllostericMechanism, re::BitVector, groups::Vector{Int}) =
    RawAllostericMechanism(_written(am.catalytic, re, groups), _label_tags(am, groups),
                           am.catalytic_multiplicity, RawSite.(am.regulatory_sites))

"""The tag of each label of `groups`, a kinetic group label per step of `am`'s catalytic
graph: the tag of the group of the label's first step."""
_label_tags(am::AllostericMechanism, groups::Vector{Int}) =
    Symbol[am.mwc_tags[am.catalytic.steps[findfirst(==(l), groups)].group]
           for l in 1:maximum(groups; init = 0)]

"""
    RawSite(s::RegulatorySite)

Regulatory site `s` as written: its ligands, multiplicity and tags, each vector fresh.
"""
RawSite(s::RegulatorySite) = RawSite(copy(s.ligands), s.multiplicity, copy(s.mwc_tags))

# ─── Writers ────────────────────────────────────────────────────────────

"""
    with_flags(m, re::BitVector)

`m` with step `k` of its catalytic graph at rapid equilibrium where `re[k]` holds, every
step in its kinetic group: `regroup` with `m`'s own groups.
"""
with_flags(m::_AnyMech, re::BitVector) = regroup(m, _group_labels(catalytic(m)), re)

"""
    regroup(m, groups::Vector{Int}, re::BitVector)

`m` with step `k` of its catalytic graph in the kinetic group labelled `groups[k]` and at
rapid equilibrium where `re[k]` holds (by default, where it already is); the labels run
from 1 and each is used. Each new
group of an allosteric mechanism takes the tag of its first step's group, so the parts of
a split group keep its tag. The rest is `m`'s. The child is built and checked by the
nodes-and-steps constructors, so it is in canonical form: a new grouping can turn a tied
step and renumber the steps.
"""
regroup(m::Mechanism, groups::Vector{Int},
        re::BitVector = BitVector(e.rapid_equilibrium for e in m.steps)) =
    Mechanism(m.substrates, m.products, m.nodes, _relabelled(m, groups, re))
function regroup(am::AllostericMechanism, groups::Vector{Int},
                 re::BitVector = BitVector(e.rapid_equilibrium for e in am.catalytic.steps))
    cm = am.catalytic
    AllostericMechanism(cm.substrates, cm.products, cm.nodes, _relabelled(cm, groups, re),
                        _label_tags(am, groups), am.catalytic_multiplicity,
                        am.regulatory_sites)
end

"""The steps of `m`, step `k` labelled `groups[k]` and at rapid equilibrium where `re[k]`
holds."""
_relabelled(m::Mechanism, groups::Vector{Int}, re::BitVector) =
    Step[Step(e.from, e.to, e.taken_up, e.given_off, re[k], groups[k])
         for (k, e) in enumerate(m.steps)]

"""
    with_tags(am::AllostericMechanism, tags::Vector{Symbol})

`am` with kinetic group `g` of its catalytic graph tagged `tags[g]`, through the composition
constructor, which checks the tags; the catalytic graph is `am`'s, shared.
"""
with_tags(am::AllostericMechanism, tags::Vector{Symbol}) =
    AllostericMechanism(am.catalytic, tags, am.catalytic_multiplicity, am.regulatory_sites)
