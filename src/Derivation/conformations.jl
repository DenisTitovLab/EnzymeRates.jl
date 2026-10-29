# ABOUTME: The allosteric conformations of a mechanism: the steps each one excludes and
# ABOUTME: keeps (`excludes`, `exclusion_components`, `conformation_steps`).

"""
    excludes(m, k, c::Symbol) -> Bool
    excludes(tag::Symbol, c::Symbol) -> Bool

Whether conformation `c` (`:active` or `:inactive`) of `m` excludes step `k` of its
catalytic graph: an `:OnlyA` group's step in the inactive conformation, an `:OnlyI`
group's in the active one. There a plain binding loses its affinity and any other step is
blocked (spec section 7). A plain mechanism excludes no step.

With a tag, whether conformation `c` excludes a step or a regulator ligand tagged `tag`:
the one rule by which the steps, the regulator constants and the moves' ligand states
read the tags.
"""
excludes(am::Mechanisms.AllostericMechanism, k::Int, c::Symbol) =
    excludes(am.mwc_tags[am.catalytic.steps[k].group], c)
excludes(::Mechanisms.Mechanism, ::Int, c::Symbol) = (_conformation_index(c); false)
excludes(tag::Symbol, c::Symbol) = tag === (_conformation_index(c) == 2 ? :OnlyA : :OnlyI)

"""
    exclusion_components(am, c::Symbol) -> Vector{Int}

Each form's component under the edges of `am`'s catalytic graph that conformation `c`
(`:active` or `:inactive`) keeps (`excludes`), every form included, labelled by its smallest
form index, so free E's label is 1. Admission reads each plain binding `c` excludes as an
arc from its unbound side's component to its bound side's (`_arcs`).
"""
function exclusion_components(am::Mechanisms.AllostericMechanism, c::Symbol)
    m = am.catalytic
    parent = collect(eachindex(m.nodes))
    top(i) = (while parent[i] != i; i = parent[i]; end; i)
    for (k, e) in enumerate(m.steps)
        excludes(am, k, c) && continue
        a, b = top(e.from), top(e.to)
        a == b || (parent[max(a, b)] = min(a, b))
    end
    [top(f) for f in eachindex(m.nodes)]
end

"""
    conformation_steps(m, c::Symbol) -> BitVector

The steps of conformation `c` (`:active` or `:inactive`) of `m`, one flag per step in
`m`'s stored direction. For an allosteric mechanism these are the steps of its catalytic
graph that `c` does not exclude (`excludes`: an `:OnlyI` group's from the active
conformation, an `:OnlyA` group's from the inactive one) between forms they join to free
E, form 1. The conformations interconvert at free E alone, so a form those steps do not
join to it holds no enzyme of that conformation, a second free form included. The law
deletes the excluded steps and prunes the forms this strands; admission decides when that
is the limit of the vanishing constants (spec section 7). A plain mechanism has only the
active conformation: every step is active and none inactive. The laws (`_law`), the joint
solve (`_ordered_columns`), the guard's law network (`_law_network`), the copy rule
(`_redundant_copy_groups`), the degeneracy test (`_degenerate`), the tagging move's
stopping test (`_stops_inactive_turnover`) and the hyperbolic test
(`_hyperbolic_catalysis`) read each conformation from it.
"""
function conformation_steps(am::Mechanisms.AllostericMechanism, c::Symbol)
    comp = exclusion_components(am, c)
    BitVector([!excludes(am, k, c) && comp[e.from] == 1
               for (k, e) in enumerate(am.catalytic.steps)])
end
conformation_steps(m::Mechanisms.Mechanism, c::Symbol) =
    (_conformation_index(c) == 1 ? trues : falses)(length(m.steps))

"""The index of conformation `c` that `excludes` reads: 1 for `:active`, 2 for
`:inactive`."""
_conformation_index(c::Symbol) =
    c === :active ? 1 : c === :inactive ? 2 :
    throw(ArgumentError("a conformation is :active or :inactive, not :$c"))
