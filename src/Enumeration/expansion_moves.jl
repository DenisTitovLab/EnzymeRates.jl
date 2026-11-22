# ABOUTME: A candidate's children, candidate_mechanisms(reaction, parent), and the moves
# ABOUTME: that make them: RE→SS flips, group splits, regulators, states, site merges.

# ─── candidate_mechanisms(reaction, parent) ─────────────────────────────

"""
    candidate_mechanisms(reaction; max_param_count = typemax(Int),
                         eq_complexity_filter = typemax(Int),
                         optional_allosteric_regulators = Symbol[],
                         optional_competitive_inhibitors = Symbol[])
        -> Vector{CandidateMechanism}
    candidate_mechanisms(reaction, parent::CandidateMechanism;
                         max_param_count = typemax(Int),
                         eq_complexity_filter = typemax(Int))
        -> Vector{CandidateMechanism}

The candidate mechanisms of the `EnzymeReaction` `reaction`. The first method returns the
mechanisms a search fits first: they bind every regulator `reaction` declares but those
the two optional lists name, and no regulator those lists name. The second returns the
children of `parent`, a candidate of `reaction`, by every expansion move: each changes
one thing, such as a step's speed, a kinetic group's split, a bound inhibitor or an
allosteric property.

Every candidate returned has a well-posed rate law (`Derivation.check_valid_rate_equation`)
with at most `eq_complexity_filter` terms (`Derivation.rate_equation_term_count`) and at
most `max_param_count` fitted constants (`Derivation.fitted_parameter_count`). The limits
drop only what a method returns: a mechanism over a limit is never returned, so it is
never expanded either.

Errors: the reaction's metabolite names spell one string two ways (such as `Ac`, `CoA`
and `AcCoA`), so that two forms would share a name (first method); the reaction has no
catalytic cycle; no mechanism binds every required regulator; an internal check of the
enumerator fails, such as atom conservation.
"""
function candidate_mechanisms(reaction::Mechanisms.EnzymeReaction,
                              parent::CandidateMechanism;
                              max_param_count::Int = typemax(Int),
                              eq_complexity_filter::Int = typemax(Int))
    _assert_emission_rules(parent)
    children = CandidateMechanism[]
    append!(children, _expand_re_to_ss(parent))
    append!(children, _expand_split_kinetic_group(parent))
    append!(children, _expand_add_dead_end_regulator(parent, reaction))
    if !_is_allosteric(parent.mechanism)
        append!(children, _expand_to_allosteric(parent, reaction))
    else
        append!(children, _expand_add_allosteric_regulator(parent, reaction),
                _expand_change_allo_state(parent), _expand_merge_regulatory_sites(parent))
    end
    filter!(c -> _respects_reg_type(c.mechanism, reaction), children)
    foreach(child -> _assert_atom_conserving(child, reaction), children)
    # The term count reads the graph alone, where the guard derives, so the term limit
    # runs first and `_within_size_limits` applies only the other.
    eq_complexity_filter == typemax(Int) || filter!(
        c -> Derivation.rate_equation_term_count(c.mechanism) <= eq_complexity_filter,
        children)
    filter!(c -> _enumerator_keeps(c.mechanism), children)
    _within_size_limits(children, max_param_count, typemax(Int))
end

# ─── The emission rules a parent must satisfy ───────────────────────────

"""
    _assert_emission_rules(c::CandidateMechanism)

Error unless the parent candidate `c` obeys the two rules every child of the moves obeys.
Every steady-state kinetic group holds an edge that carries net flux
(`_flux_carrying_groups`), read on the whole catalytic graph. That graph is the active
conformation's: the active conformation leaves out only the groups tagged `:OnlyI` and the
steps they cut off from free E, and no move tags a catalytic group `:OnlyI`. No kinetic
group binds a competitive inhibitor redundantly (`_redundant_copy_groups`). A parent must
obey both rules because a flip tests only the groups it flips. The split, the dead-end move
and `_expand_change_allo_state` filter their children; the other moves preserve both rules.
"""
function _assert_emission_rules(c::CandidateMechanism)
    cm = c.mechanism.catalytic
    members = _kinetic_groups(cm)
    # A group as its steps' DSL lines (`step_text`).
    label(g) = join((Mechanisms.step_text(cm, k) for k in members[g]), ", ")
    flux = _flux_carrying_groups(cm, _rapid_equilibrium_flags(cm))
    for (g, ks) in enumerate(members)
        cm.steps[first(ks)].rapid_equilibrium || flux[g] || error(
            "candidate_mechanisms: steady-state kinetic group {" * label(g) * "} has no " *
            "step that carries net flux, so its two constants enter the rate only as " *
            "their ratio; write the group at rapid equilibrium")
    end
    for g in _redundant_copy_groups(c)
        error("candidate_mechanisms: kinetic group {" * label(g) * "} binds a " *
              "competitive inhibitor only where the complex duplicates a productive " *
              "form, and a dwell gauge absorbs its constant into the existing binding's; " *
              "bind the inhibitor where it forms a new complex, or drop it")
    end
    nothing
end

# ─── Kinetic groups ─────────────────────────────────────────────────────

"""Whether each kinetic group of graph `m` holds some step `k` with `pred(k)`."""
function _groups_with(pred, m::Mechanisms.CatalyticScheme)
    out = falses(maximum(e -> e.kinetic_group, m.steps; init = 0))
    for k in eachindex(m.steps)
        pred(k) && (out[m.steps[k].kinetic_group] = true)
    end
    out
end

"""Whether each kinetic group of the catalytic graph of candidate `c` holds a chemistry
step: anything but a plain binding (`_is_plain_binding`). The allosteric moves tag a
chemistry group apart from a binding group."""
_chemistry_groups(c::CandidateMechanism) =
    _groups_with(k -> !_is_plain_binding(c, k), c.mechanism.catalytic)

"""
Candidate `c` with step `k` of its catalytic graph in the kinetic group labelled
`groups[k]` and at rapid equilibrium where `re[k]` holds; the labels run from 1 and each
is used. Each new group takes the allosteric state of its first step's group
(`_label_tags`), so the parts of a split group keep its state. The rest is `c`'s. The
child is built and checked by the constructors (`_candidate`), so it is in canonical form:
a new grouping can turn a tied step and renumber the steps.
"""
function _regrouped(c::CandidateMechanism, groups::Vector{Int}, re::BitVector)
    m = c.mechanism
    cm, a = m.catalytic, m.allosteric
    _candidate(cm.metabolites, c.forms, cm.substrates, cm.products,
               _regrouped_steps(cm, groups, re),
               Mechanisms.Allostery(_label_tags(m, groups), a.oligomeric_state,
                                    a.regulatory_sites))
end

"""Candidate `c` with step `k` of its catalytic graph at rapid equilibrium where `re[k]`
holds, every step in its kinetic group (`_regrouped` with `c`'s own groups)."""
_with_flags(c::CandidateMechanism, re::BitVector) =
    _regrouped(c, Int[e.kinetic_group for e in c.mechanism.catalytic.steps], re)

_regrouped_steps(m::Mechanisms.CatalyticScheme, groups::Vector{Int}, re::BitVector) =
    Mechanisms.Step[Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off, re[k], groups[k])
                    for (k, e) in enumerate(m.steps)]

"""The allosteric state of each label of `groups`, a kinetic group label per step of
`m`'s catalytic graph: the state of the group of the label's first step. A plain
mechanism's groups have no state, so neither do its labels."""
function _label_tags(m::Mechanisms.Mechanism, groups::Vector{Int})
    states = m.allosteric.allosteric_states
    isempty(states) && return Symbol[]
    Symbol[states[m.catalytic.steps[findfirst(==(l), groups)].kinetic_group]
           for l in 1:maximum(groups; init = 0)]
end

# ─── Expansion moves ────────────────────────────────────────────────────

"""
    _inhibitor_competition_patterns(sub_names, prod_names, existing_inhibitors)
        -> Vector{Tuple{Set{Symbol},Set{Symbol},Set{Symbol}}}

Every pattern of competition for a new inhibitor, as the substrates, the products and the
existing inhibitors it competes with: at least one substrate, at least one product, and any
subset of the existing inhibitors.
"""
function _inhibitor_competition_patterns(
    sub_names::Set{Symbol},
    prod_names::Set{Symbol},
    existing_inhibitors::Vector{Symbol},
)
    subs = _subsets(sort(collect(sub_names)))[2:end]
    prods = _subsets(sort(collect(prod_names)))[2:end]
    inhs = _subsets(sort(existing_inhibitors))
    [(Set(s), Set(p), Set(i)) for s in subs for p in prods for i in inhs]
end

"""
    _expand_re_to_ss(c::CandidateMechanism) → Vector{CandidateMechanism}

The flip move: one child of candidate `c` per inclusion-minimal set of flip units whose
joint flip to steady state (`_minimal_flips`) passes three tests:

- it raises the number of rapid-equilibrium segments, since a steady-state step whose two
  ends share a segment never reaches the rate law;
- every flipped group carries flux in the child;
- `_enumerator_keeps` accepts the child.

A flip unit is a kinetic group that is all rapid equilibrium, takes up or gives off no
modifier copy, carries flux with every edge steady state (`_flux_carrying_groups`) and is
no chain flank (`_chain_flank_groups`). The child keeps everything else of `c`. A parent
that `_requires_hyperbolic_catalysis` keeps only the children that pass
`_hyperbolic_catalysis`. The engine page gives the reasons for each rule ("Flip
rapid-equilibrium groups to steady state" in "The enumeration engine").

The search extends a set that fails a test. A set that leaves a flipped group without flux
gives the parent plus one phantom per such group, and its flux-carrying supersets are the
smallest identifiable relaxations. The mechanisms that are not degenerate beyond a failed
set stay reachable, through the other cut orders of the same ring, or through the extended
set when no other order gains. The hyperbolic test runs on the emitted sets, not inside the
search: there a failing set would be extended with more flips, and more steady-state steps
never restore a hyperbolic equation, so its supersets need no visit.
"""
function _expand_re_to_ss(c::CandidateMechanism)
    cm = c.mechanism.catalytic
    flanks = _chain_flank_groups(c.mechanism)
    # A step takes up or gives off a modifier copy exactly when its two forms hold
    # different copies.
    exchanges_modifier = _groups_with(cm) do k
        c.forms[cm.steps[k].from].copies != c.forms[cm.steps[k].to].copies
    end
    eligible(g) = !(g in flanks) && !exchanges_modifier[g]
    base = _re_segment_count(cm, _rapid_equilibrium_flags(cm))
    flips = _minimal_flips(re -> _re_segment_count(cm, re) > base, c, eligible)
    children = CandidateMechanism[_with_flags(c, re) for re in flips]
    filter(_hyperbolic_where_required, children)
end

"""
Kinetic groups of `m` that are a flank of a qualifying chain whose isomerization is steady
state: an isomerization X1 → X2 alone in its group, where X1 and X2 each have exactly one
other step, each a binding into that form and alone in its group. The rest of the
mechanism sees such a chain only through its flux and enzyme content, four numbers that
the form with both flanks at rapid equilibrium already covers, so flipping a flank adds
parameters and no rate law. The tests call this the chain lemma; the docs page "Rapid
equilibrium vs steady state" states it under "Apparent constants on a steady-state chain".
An allosteric mechanism has none: at more than one catalytic subunit a flank's flip can be
visible through the two conformations.
"""
function _chain_flank_groups(am::Mechanisms.Mechanism)
    _is_allosteric(am) && return Set{Int}()
    m = am.catalytic
    # The group of each edge alone in its group, 0 for any other edge.
    lone = zeros(Int, length(m.steps))
    for (g, ks) in enumerate(_kinetic_groups(m))
        length(ks) == 1 && (lone[only(ks)] = g)
    end
    at = [Int[] for _ in eachindex(m.enzyme_forms)]
    for k in eachindex(m.steps)
        (; from, to) = m.steps[k]
        push!(at[from], k); push!(at[to], k)
    end
    # The group of the flank of edge `k0` at form `x`, 0 when there is none.
    function flank(k0, x)
        rest = filter(!=(k0), at[x])
        length(rest) == 1 || return 0
        k = only(rest)
        _is_binding(m.steps[k]) && m.steps[k].to == x ? lone[k] : 0
    end
    out = Set{Int}()
    for k0 in eachindex(m.steps)
        lone[k0] != 0 && _is_iso(m.steps[k0]) &&
            !m.steps[k0].rapid_equilibrium || continue
        f1, f2 = flank.(k0, (m.steps[k0].from, m.steps[k0].to))
        f1 == 0 || f2 == 0 || union!(out, (f1, f2))
    end
    out
end

"""
    _minimal_flips(admissible, c, eligible = _ -> true) -> Vector{BitVector}

The rapid-equilibrium flags of the catalytic graph of candidate `c`, plain or allosteric,
after every inclusion-minimal flip to steady state that passes, one flag per edge for each
flip set. A unit is a rapid-equilibrium kinetic group that carries flux with every edge
steady state (`_flux_carrying_groups`) and passes `eligible`. A candidate turns its units'
edges steady state and keeps every other flag. It passes when `admissible` holds on its
flags, every flipped group carries flux in it, and `_enumerator_keeps` accepts the mechanism
it builds.
"""
function _minimal_flips(admissible, c::CandidateMechanism, eligible = _ -> true)
    cm = c.mechanism.catalytic
    flags = _rapid_equilibrium_flags(cm)
    members = _kinetic_groups(cm)
    flux = _flux_carrying_groups(cm, falses(length(cm.steps)))
    units = [g for g in eachindex(members)
             if all(k -> flags[k], members[g]) && flux[g] && eligible(g)]
    function flipped(sel)
        re = copy(flags)
        for u in units[sel], k in members[u]
            re[k] = false
        end
        re
    end
    passes(sel) = (re = flipped(sel);
                   admissible(re) && all(_flux_carrying_groups(cm, re)[units[sel]]) &&
                   _enumerator_keeps(_with_flags(c, re).mechanism))
    [flipped(sel) for sel in _minimal_gaining_sets(passes, _ -> 1:length(units))]
end

"""
    _flux_carrying_edges(m, re) -> BitVector

Whether each edge of graph `m`, edge `k` at rapid equilibrium where `re[k]` holds,
carries net steady-state flux for generic parameter values. Rapid-equilibrium edges hold
equilibrium mass and are flagged `false`; the rule that uses these flags concerns
steady-state groups only.

Contract each rapid-equilibrium segment to one vertex (`rapid_equilibrium_segments`) and
make every steady-state edge an edge between its two forms' segments, weighted by the
substrates it takes up minus the products (`_uptake_weight`), plus its first end's offsets
minus its second end's (`from`, `to`; a modifier copy counts as its reactant). Around any
cycle the offsets telescope and the uptakes sum to the net turnover times the reactant
count, so a cycle has weight zero exactly when it runs no net reaction. A self-loop (both
forms in one segment) carries flux iff its weight is nonzero. Any other edge carries flux
iff its biconnected block (`_unbalanced_blocks`) holds a cycle of nonzero weight: every edge
of such a block lies on one (join the edge to the cycle by two disjoint paths; one of the
two resulting cycles is unbalanced), while in a balanced block detailed balance holds along
every cycle and each edge's flux vanishes. The zero-flux verdict holds for any parameters
and grouping; the flux-carrying verdict for one-way steps, which every step the enumerator
emits is. The test reads each edge's metabolite lists, so it holds for fused and
Theorell–Chance steps as well.
"""
function _flux_carrying_edges(m::Mechanisms.CatalyticScheme, re::BitVector)
    s = Mechanisms.rapid_equilibrium_segments(m, trues(length(m.steps)), re)
    rho = _name_signs(m)
    offset(f) = sum(rho[n] * s.offsets[n, f] for n in eachindex(rho); init = 0)
    flags = falses(length(m.steps))
    edges = Tuple{Int, Int}[]; weights = Int[]; owner = Int[]
    for k in eachindex(m.steps)
        re[k] && continue
        (; from, to) = m.steps[k]
        weight = _uptake_weight(m, k) + offset(from) - offset(to)
        a, b = s.segment_of[from], s.segment_of[to]
        if a == b
            flags[k] = weight != 0
        else
            push!(edges, (a, b)); push!(weights, weight); push!(owner, k)
        end
    end
    block, unbalanced = _unbalanced_blocks(length(s.segments), edges, weights)
    for (j, k) in enumerate(owner)
        flags[k] = block[j] in unbalanced
    end
    flags
end

"""
    _flux_carrying_groups(m, re) -> BitVector

One flag per kinetic group of graph `m`, edge `k` at rapid equilibrium where `re[k]`
holds: whether some edge of the group carries net steady-state flux
(`_flux_carrying_edges`). A steady-state group with no such edge exposes only the ratio
of its two constants, so the moves never emit one; a zero-flux edge inside a group that
also holds a flux-carrying edge costs nothing, because the group's shared constants are
pinned by the edge that carries flux.
"""
function _flux_carrying_groups(m::Mechanisms.CatalyticScheme, re::BitVector)
    flux = _flux_carrying_edges(m, re)
    _groups_with(k -> flux[k], m)
end

"""Each metabolite's net uptake sign, in the order of `m.metabolites` (the rows of the
offsets of `rapid_equilibrium_segments`): +1 for a substrate, −1 for a product, 0 for a
modifier."""
_name_signs(m::Mechanisms.CatalyticScheme) =
    Int[n in m.substrates ? 1 : n in m.products ? -1 : 0 for n in eachindex(m.metabolites)]

"""
    _unbalanced_blocks(nv, edges, weights) -> (block, unbalanced)

The biconnected block of every edge of an undirected multigraph without self-loops,
`edges[e] = (u, v)` with vertices `1:nv` (Tarjan's edge-stack algorithm; a bridge is a
block of its own), and the set of blocks holding a cycle of nonzero weight, edge `e`
weighing `weights[e]` along its stored direction. The search gives each vertex a
potential that rises by each tree edge's weight along its stored direction. An edge's
defect, its weight minus the potential difference of its ends, is zero on a tree edge
and, on any other edge, the weight of the cycle it closes with the tree path between its
ends, a cycle of its own block. The tree edges of a block span it, so these cycles
generate the block's cycles and weight adds over them: a block holds a cycle of nonzero
weight exactly when one of its edges has a nonzero defect.
"""
function _unbalanced_blocks(nv::Int, edges::Vector{Tuple{Int, Int}}, weights::Vector{Int})
    adj = [Int[] for _ in 1:nv]
    for (e, (u, v)) in enumerate(edges)
        push!(adj[u], e); push!(adj[v], e)
    end
    disc = zeros(Int, nv); low = zeros(Int, nv); phi = zeros(Int, nv)
    block = zeros(Int, length(edges))
    stack = Int[]; clock = Ref(0); nblocks = Ref(0)
    # It takes itself as an argument because a closure that names itself is boxed.
    function visit(self, u, parent_edge)
        clock[] += 1; disc[u] = low[u] = clock[]
        for e in adj[u]
            e == parent_edge && continue
            w = edges[e][1] == u ? edges[e][2] : edges[e][1]
            if disc[w] == 0
                phi[w] = phi[u] + (edges[e][1] == u ? weights[e] : -weights[e])
                push!(stack, e)
                self(self, w, e)
                low[u] = min(low[u], low[w])
                if low[w] >= disc[u]
                    nblocks[] += 1
                    while true
                        x = pop!(stack); block[x] = nblocks[]
                        x == e && break
                    end
                end
            elseif disc[w] < disc[u]
                push!(stack, e)
                low[u] = min(low[u], disc[w])
            end
        end
    end
    for v in 1:nv
        disc[v] == 0 && visit(visit, v, 0)
    end
    block, Set{Int}(block[e] for (e, (u, v)) in enumerate(edges)
                    if weights[e] != phi[v] - phi[u])
end

_re_segment_count(m::Mechanisms.CatalyticScheme, re::BitVector) =
    length(Mechanisms.rapid_equilibrium_segments(m, trues(length(m.steps)), re).segments)

"""
    _minimal_gaining_sets(gains, partners) -> Vector{Vector{Int}}

Every minimal set of units for which `gains(set)` holds, found Apriori-style.
`partners(set)` lists the units that may extend `set`, and `partners(Int[])` the units
tried alone. `gains` runs on a set only after failing on each of its subsets one unit
smaller, so never on a superset of a gaining set. Each set is sorted; they come by size,
then in lexicographic order.
"""
function _minimal_gaining_sets(gains, partners)
    out = Vector{Int}[]
    failed = Vector{Int}[Int[]]
    while !isempty(failed)
        failed_keys = Set(failed)
        candidates = Set{Vector{Int}}()
        for set in failed, u in partners(set)
            u in set && continue
            c = sort!(vcat(set, u))
            c in candidates && continue
            all(setdiff(c, [x]) in failed_keys for x in c) || continue
            push!(candidates, c)
        end
        failed = Vector{Int}[]
        for c in sort!(collect(candidates))
            gains(c) ? push!(out, c) : push!(failed, c)
        end
    end
    out
end

"""
    _expand_split_kinetic_group(c::CandidateMechanism) → Vector{CandidateMechanism}

The split move: one child of candidate `c` per inclusion-minimal set of context
bipartitions (`_context_bipartitions`), at most one per group, that passes two tests.
Applied together (`_apply_bipartitions`), the set raises the fitted-parameter count, which
the joint solve gives without deriving a rate polynomial (`fitted_parameter_count`). And
`_enumerator_keeps` accepts the child. Each split group's second part becomes a group of its
own, with the group's allosteric state. A part of a steady-state group none of whose edges
carries flux turns rapid equilibrium (`_revert_zero_flux_parts`): the same family with one
constant fewer. The child keeps everything else of `c`. The move drops a child that holds a
redundant copy group (`_redundant_copy_groups`), whose family is that of the same split
without the copy.

A single bipartition whose new constant a Wegscheider cycle ties straight back gains
nothing, but the search extends it, since splitting another group on that cycle breaks the
tie: random-order binding needs one split per substrate before any constant frees up. It
extends a failed set only by the rapid-equilibrium groups that share a segment with one of
its split groups: a tie needs a rapid-equilibrium cycle through the split, and every group
on that cycle lies in that segment.

A copy is redundant when each of its complexes has a twin and a dwell gauge exists. A twin
is a productive form whose weight is proportional to the complex's (`_productive_twin`). A
dwell gauge is a rescaling of the constants that moves the copy's constant where the rate
law cannot see it (`_gauge_rescaling`). Bindings elsewhere in the group of a binding that
forms or leaves a twin can block the gauge. A split that separates them completes the
gauge, and the copy, which the shared group pinned before, becomes redundant. A split of a
group that forms a copy's twins can break the gauge. So the move tests each child as a
whole.
"""
function _expand_split_kinetic_group(c::CandidateMechanism)
    cm = c.mechanism.catalytic
    members = _kinetic_groups(cm)
    flags = _rapid_equilibrium_flags(cm)
    # The parent's flux serves every child, since a split moves no edge.
    flux = _flux_carrying_edges(cm, flags)
    units = Tuple{Int, Tuple{Vector{Int}, Vector{Int}}}[]
    reverted = Vector{Int}[]
    for g in eachindex(members), bp in _context_bipartitions(c, g)
        push!(units, (g, bp)); push!(reverted, _revert_zero_flux_parts(cm, bp, flux))
    end
    isempty(units) && return CandidateMechanism[]
    # The flags of a candidate: the parent's, with each selected unit's reverted edges at
    # rapid equilibrium.
    function flags_of(sel)
        re = copy(flags)
        for u in sel
            re[reverted[u]] .= true
        end
        re
    end
    base = Derivation.fitted_parameter_count(c.mechanism)
    # The child of each set the search accepts, built once and emitted as judged.
    accepted = Dict{Vector{Int}, CandidateMechanism}()
    function gains(sel)
        child = _apply_bipartitions(c, units[sel], flags_of(sel))
        Derivation.fitted_parameter_count(child.mechanism) > base &&
            _enumerator_keeps(child.mechanism) || return false
        accepted[sel] = child
        true
    end
    every = trues(length(cm.steps))
    # RE segment ids touched by each kinetic group's edges; empty for a steady-state
    # group, which lies on no RE cycle and so never partners a set that does not gain.
    segment_of = Mechanisms.rapid_equilibrium_segments(cm, every, flags).segment_of
    segments = [flags[first(ks)] ?
                    Set{Int}(segment_of[cm.steps[k].from] for k in ks) :
                    Set{Int}() for ks in members]
    partners(sel) = begin
        isempty(sel) && return 1:length(units)
        used = Set(units[u][1] for u in sel)
        touched = reduce(union, (segments[units[u][1]] for u in sel))
        [u for u in 1:length(units)
         if !(units[u][1] in used) && !isempty(intersect(segments[units[u][1]], touched))]
    end
    filter!(child -> isempty(_redundant_copy_groups(child)),
            CandidateMechanism[accepted[sel]
                               for sel in _minimal_gaining_sets(gains, partners)])
end

"""
The end of step `k` of `m` that does not carry the step's free metabolites: the form they
bind to. Its first end (`from`) when the step takes a metabolite up (every binding, plain
or fused, is stored with its metabolite taken up), its second when it only gives
metabolites off (two or more: a step giving off one is stored as the binding it
reverses), and its first for an isomerization (both lists empty).
"""
function _context_form(m::Mechanisms.CatalyticScheme, k::Int)
    (; from, to) = m.steps[k]
    isempty(m.steps[k].taken_up) && !isempty(m.steps[k].given_off) ? to : from
end

"""
    _context_bipartitions(c, g) -> Vector{Tuple{Vector{Int}, Vector{Int}}}

The ways to divide kinetic group `g` of the catalytic graph of candidate `c` in two by
binding context, each a pair of edge lists. A context value is a ligand bound at an edge's
context form (`_context_form`), the form's conformation or its covalent residual; it
divides the edges whose context form carries it from the rest. A ligand is a name bound
productively or as a modifier copy, so a reactant and its copy are two ligands, and a
modifier separates a catalytic step from its inhibitor-bound mirror. The group's own bound
ligand is no context. Both parts are nonempty, the first holds the group's first edge, each
lists its edges in edge order, and values that divide alike give one bipartition. The
bipartitions come in a fixed order, so the move's output is deterministic: by ligand
(modifiers, then products, then substrates, each by name), then by conformation (by label),
then by residual (by the names it adds, then those it subtracts).
"""
function _context_bipartitions(c::CandidateMechanism, g::Int)
    m = c.mechanism.catalytic
    ks = _kinetic_groups(m)[g]
    # The group's own bound ligand, as (name, copy), when its first edge is a binding.
    k1 = first(ks)
    own = _is_binding(m.steps[k1]) ?
        (n = m.metabolites[only(m.steps[k1].taken_up)]; (n, _copy_change(c, k1, n) > 0)) :
        nothing
    forms = [_context_form(m, k) for k in ks]
    carries(x, (name, copy)) = name in (copy ? c.forms[x].copies : c.forms[x].productive)
    ligands = unique!([(n, copy) for x in forms
                       for n in [c.forms[x].productive; c.forms[x].copies]
                       for copy in (false, true)
                       if carries(x, (n, copy)) && (n, copy) != own])
    products = m.metabolites[m.products]
    role((name, copy)) = copy ? 0 : name in products ? 1 : 2
    division(holds) = Int[k for (k, x) in zip(ks, forms) if holds(x)]
    divisions = [division(x -> carries(x, l))
                 for l in sort!(ligands; by = l -> (role(l), first(l)))]
    append!(divisions, division(x -> c.forms[x].label === label)
            for label in sort!(unique(c.forms[x].label for x in forms)))
    residual(x) = (c.forms[x].adds, c.forms[x].subtracts)
    append!(divisions, division(x -> residual(x) == r)
            for r in sort!(unique(residual(x) for x in forms)))
    parts(with) = (without = Int[k for k in ks if !(k in with)];
                   first(ks) in with ? (with, without) : (without, with))
    unique!(first, Tuple{Vector{Int}, Vector{Int}}[
        parts(with) for with in divisions if 0 < length(with) < length(ks)])
end

"""
`m` regrouped by `selection`, pairs of a kinetic group of `m` and one of its
`_context_bipartitions`, each group at most once: edge `k` in group `group[k]`, the groups
numbered as `m`'s with each selected group's second part in a group of its own right after
it.
"""
function _bipartitioned_groups(m::Mechanisms.CatalyticScheme, selection)
    second = Dict(g => last(bp) for (g, bp) in selection)
    group, n = zeros(Int, length(m.steps)), 0
    for (g, ks) in enumerate(_kinetic_groups(m))
        group[ks] .= (n += 1)
        haskey(second, g) && (group[second[g]] .= (n += 1))
    end
    group
end

"""
The child of candidate `c` with each selected group divided in two (`_bipartitioned_groups`)
and edge `k` at rapid equilibrium where `re[k]` holds, built through the checked constructor
(`_regrouped`). For an allosteric mechanism both parts inherit the group's catalytic
allosteric state: the split relaxes parameters and must not change which conformations
run a step.
"""
function _apply_bipartitions(c::CandidateMechanism, selection,
                             re::BitVector = _rapid_equilibrium_flags(
                                 c.mechanism.catalytic))
    _regrouped(c, _bipartitioned_groups(c.mechanism.catalytic, selection), re)
end

"""
The steady-state edges of the bipartition `bp` of a group of `m` that the split emits at
rapid equilibrium: those of each part none of whose edges carries flux (`flux`, the
parent's `_flux_carrying_edges`). A rapid-equilibrium group's edges carry none and are
rapid equilibrium already, so it reverts none.
"""
_revert_zero_flux_parts(m::Mechanisms.CatalyticScheme, bp, flux::BitVector) =
    Int[k for part in bp if !any(k -> flux[k], part)
        for k in part if !m.steps[k].rapid_equilibrium]

"""
The forms of the catalytic graph of candidate `c` where metabolite `name`, as its modifier
copy when `copy` holds (`_copy_change`), is free on some edge: the edge's first end when it
takes the metabolite up, its second when it gives it off, so a Theorell–Chance step EA + B
→ EQ + P frees B at EA and P at EQ. A reactant and its modifier copy share a name but
bind different sites: the reactant asks where the metabolite binds productively, the copy
where it binds its dead-end site. A metabolite `c` does not hold is free nowhere.
"""
function _forms_where_free(c::CandidateMechanism, name::Symbol, copy::Bool)
    cm = c.mechanism.catalytic
    free = falses(length(c.forms))
    n = something(findfirst(==(name), cm.metabolites), 0)
    for k in eachindex(cm.steps)
        (; from, to, taken_up, given_off) = cm.steps[k]
        # How many of `name` the step takes up and gives off in the role asked: `d` are
        # copies, the rest bind productively.
        d = _copy_change(c, k, name)
        up, off = copy ? (max(d, 0), max(-d, 0)) :
            (_copies(taken_up, n) - max(d, 0), _copies(given_off, n) - max(-d, 0))
        up > 0 && (free[from] = true)
        off > 0 && (free[to] = true)
    end
    free
end

"""
The metabolites, the form records and the steps of the dead-end move's child of candidate
`c` that binds modifier `reg` at each form of `sites`:

- `c`'s metabolites, with `reg` after them unless it is one;
- `c`'s records, with one complex after them per site, in the site's conformation and
  residual;
- `c`'s steps in their groups, with each step between two sites mirrored between their
  complexes in the step's own group, and one new group, after `c`'s, of the
  rapid-equilibrium bindings of `reg` at the sites.
"""
function _dead_end_graph(c::CandidateMechanism, reg::Symbol, sites::Vector{Int})
    m = c.mechanism.catalytic
    metabolites = reg in m.metabolites ? m.metabolites : [m.metabolites; reg]
    inhibitor = findfirst(==(reg), metabolites)
    forms = copy(c.forms)
    complex = zeros(Int, length(c.forms))
    for f in sites
        r = c.forms[f]
        push!(forms, _FormRecord(r.label, r.productive, sort!([r.copies; reg]), r.adds,
                                 r.subtracts))
        complex[f] = length(forms)
    end
    steps = copy(m.steps)
    for e in m.steps
        all(in(sites), (e.from, e.to)) || continue
        push!(steps, Mechanisms.Step(complex[e.from], complex[e.to], e.taken_up,
                                     e.given_off, e.rapid_equilibrium, e.kinetic_group))
    end
    G = maximum(e -> e.kinetic_group, m.steps; init = 0)
    for f in sites
        push!(steps, Mechanisms.Step(f, complex[f], [inhibitor], Int[], true, G + 1))
    end
    metabolites, forms, steps
end

"""
The dead-end move's child of the mechanism `m` with the metabolites `metabolites`, the form
records `forms` and the steps `steps` (`_dead_end_graph`), binding modifier `reg` in a new
group; `nothing` when the copy is redundant in it. An allosteric child tags the new group
`:EqualAI` and keeps `m`'s oligomeric state and regulatory sites.

A plain child's new group, the group of the steps that take the copy up, goes through
`_gauge_rescaling` on the child's own twins (`_productive_twin`). The twins are the same
forms as the parent's: the complexes are not productive, and their edges attach each complex
to its site's segment and join no two segments, so every form of `m` keeps its segment and
offsets. But the complexes sort in among `m`'s forms, so `m`'s form indices do not hold in
the child. Only the new group is tested: `m` holds no redundant copy group (the parent
rule), and in each of about 218,000 enumerated parents checked, the placement made no older
copy redundant (an empirical finding, not a proof). An allosteric child goes through
`_redundant_copy_groups` only when the copy has a twin in the active conformation at each
of its bindings (`_all_twin`), as redundancy needs; its bindings are the only edges that
take it up. The twins are read on the child's whole catalytic graph, which is the active
conformation's (`_assert_emission_rules` says why).
"""
function _dead_end_child(m::Mechanisms.Mechanism, metabolites::Vector{Symbol},
                         forms::Vector{_FormRecord}, steps::Vector{Mechanisms.Step},
                         reg::Symbol)
    a = m.allosteric
    states = _is_allosteric(m) ? [a.allosteric_states; :EqualAI] : Symbol[]
    child = _candidate(metabolites, forms, m.catalytic.substrates, m.catalytic.products,
                       steps,
                       Mechanisms.Allostery(states, a.oligomeric_state, a.regulatory_sites))
    cm = child.mechanism.catalytic
    twin = _productive_twin(child, trues(length(cm.steps)))
    if !_is_allosteric(m)
        binding = findfirst(k -> _copy_change(child, k, reg) > 0, eachindex(cm.steps))::Int
        return _gauge_rescaling(child, _kinetic_groups(cm), cm.steps[binding].kinetic_group,
                                twin, _one_state_class) === nothing ? child : nothing
    end
    bindings = findall(k -> _copy_change(child, k, reg) > 0, eachindex(cm.steps))
    _all_twin(child, bindings, twin) && !isempty(_redundant_copy_groups(child)) ?
        nothing : child
end

"""
    _expand_add_dead_end_regulator(c::CandidateMechanism, rxn::EnzymeReaction)
        → Vector{CandidateMechanism}

The dead-end move: the children of candidate `c` that bind one more competitive inhibitor
of `rxn`, one that `c` binds neither as a modifier copy nor at a regulatory site. A ligand
of a regulatory site belongs there, so this move leaves it out, while the regulator move
(`_expand_add_allosteric_regulator`) takes a name bound as a copy. Each competition pattern
(`_inhibitor_competition_patterns`) places the inhibitor at the forms where a metabolite it
competes with is free (`_forms_where_free`): a substrate or product on a productive step,
an existing inhibitor at its own binding. The placement leaves out each form that holds,
productively, every substrate, every product or a substrate or product the inhibitor
competes with, and each form that holds an inhibitor it competes with. A reactant declared
as a competitive inhibitor binds as a copy at a dead-end site of its own, so a form that
carries only the copy of M is a site for an inhibitor that competes with M.

Each eligible inhibitor and each distinct nonempty placement give one child
(`_dead_end_graph`). The inhibitor binds every form of the placement at rapid equilibrium,
in one new kinetic group (one new constant). Each step between two placed forms is
mirrored between their complexes, in the step's own group. The move drops a child whose
copy is redundant (`_dead_end_child`).
"""
function _expand_add_dead_end_regulator(c::CandidateMechanism,
                                        rxn::Mechanisms.EnzymeReaction)
    sub_names = Set(s.name for s in rxn.substrates)
    prod_names = Set(p.name for p in rxn.products)
    bound_regs = _bound_comp_inhibitors(c)
    excluded = union(bound_regs, _bound_allo_regs(c.mechanism))
    eligible_regs = [reg.name for reg in rxn.regulators
                     if reg.kind === :competitive && reg.name ∉ excluded]
    isempty(eligible_regs) && return CandidateMechanism[]

    # Per form, the names bound productively and the names bound as modifiers: a
    # reactant's copy sits at a dead-end site, not at its reactant's.
    productive = [Set(f.productive) for f in c.forms]
    inhibiting = [Set(f.copies) for f in c.forms]

    # The sites of each competition pattern in form order, once each in pattern order.
    placements = Vector{Int}[]
    for (comp_subs, comp_prods, comp_inhibitors) in
            _inhibitor_competition_patterns(sub_names, prod_names, collect(bound_regs))
        comp_reactants = union(comp_subs, comp_prods)
        target_forms = falses(length(c.forms))
        for met in comp_reactants
            target_forms .|= _forms_where_free(c, met, false)
        end
        for inh in comp_inhibitors
            target_forms .|= _forms_where_free(c, inh, true)
        end
        push!(placements, [f for f in findall(target_forms)
                           if !issubset(sub_names, productive[f]) &&
                              !issubset(prod_names, productive[f]) &&
                              isempty(intersect(productive[f], comp_reactants)) &&
                              isempty(intersect(inhibiting[f], comp_inhibitors))])
    end
    unique!(filter!(!isempty, placements))

    results = CandidateMechanism[]
    for reg in eligible_regs, sites in placements
        child = _dead_end_child(c.mechanism, _dead_end_graph(c, reg, sites)..., reg)
        child === nothing || push!(results, child)
    end
    results
end

"""
    _valid_tagging(c::CandidateMechanism, tags::Vector{Symbol}) -> Bool

Whether `_enumerator_keeps` accepts the plain mechanism of candidate `c` with the tags
`tags`, one per kinetic group, over one catalytic subunit with no regulatory site: tagging
bases pass the test every other candidate passes. Neither the guard's checks nor
`_kept_forms_keep_a_binding` read the oligomeric state or the sites, so every child the
tagging move builds on these tags gets the same verdict at the exits, up to the limit
test's random draw, whose seed hashes the mechanism.
"""
_valid_tagging(c::CandidateMechanism, tags::Vector{Symbol}) =
    _enumerator_keeps(_candidate(c, tags, 1, Mechanisms.RegulatorySite[]).mechanism)

"""
    _stops_inactive_turnover(c::CandidateMechanism, tags::Vector{Symbol}) -> Bool

Whether the tags `tags` leave the inactive conformation of the plain mechanism `m` of
candidate `c` no catalytic cycle: no cycle of the steps it runs
(`steps_running_in_allosteric_conformations`) turns the reaction over
(`has_turnover_cycle`).
"""
function _stops_inactive_turnover(c::CandidateMechanism, tags::Vector{Symbol})
    am = _candidate(c, tags, 1, Mechanisms.RegulatorySite[]).mechanism
    !Mechanisms.has_turnover_cycle(
        am.catalytic, Mechanisms.steps_running_in_allosteric_conformations(am).inactive)
end

"""
    _tagging_bases(c::CandidateMechanism) -> Vector{Vector{Symbol}}

The tag sets the tagging move (`_expand_to_allosteric`) builds on, one tag per kinetic group
of candidate `c`'s plain mechanism. The first candidate tags every chemistry group
(`_chemistry_groups`) `:OnlyA` and every binding group `:EqualAI`; when it passes
`_valid_tagging` it is the only base. Otherwise the bases are the smallest sets of
chemistry groups, by inclusion, that stop every catalytic cycle of the inactive
conformation (`_stops_inactive_turnover`) and pass `_valid_tagging`, each tagged `:OnlyA`
with every other group `:EqualAI`; a set that fails leaves its supersets in the search. The
bases come by size, each size in `_subsets` order. The engine page works the ping-pong
case, where one half-reaction can be a base alone ("Promote a non-allosteric to allosteric
mechanism" in "The enumeration engine").
"""
function _tagging_bases(c::CandidateMechanism)
    chem = _chemistry_groups(c)
    tagged(onlya) = [g in onlya ? :OnlyA : :EqualAI for g in eachindex(chem)]
    groups = findall(chem)
    _valid_tagging(c, tagged(groups)) && return [tagged(groups)]
    sets = _subsets(groups)
    found = Vector{Int}[]
    for k in 1:length(groups), s in sets
        length(s) == k && !any(b -> b ⊆ s, found) || continue
        t = tagged(s)
        _stops_inactive_turnover(c, t) && _valid_tagging(c, t) && push!(found, s)
    end
    [tagged(s) for s in found]
end

"""
    _expand_to_allosteric(c::CandidateMechanism, rxn::EnzymeReaction)
        → Vector{CandidateMechanism}

The tagging move: the allosteric variants of candidate `c`, whose mechanism `m` is plain,
each built on a base of `_tagging_bases` and sharing `m`'s catalytic scheme. For each
oligomeric state of `rxn`'s `allowed_oligomeric_states` in turn, the move emits the
K-types, then the V-types:

  * K-type: a nonempty set of binding groups (those outside `_chemistry_groups`) tagged
    `:OnlyA` on top of a base, and on top of the every-chemistry tag set when that set is
    no base, since a tagged binding can restore its limit. Each is emitted bare; the exits
    drop those `_enumerator_keeps` refuses. Over one subunit `L` can be a phantom, which a
    rescaling of the other constants absorbs (the `L_phantom` column in "Best mechanism
    selection" gives the rule). The move emits such variants anyway, since their children
    can carry laws no other child does. The engine page works the cases ("Promote a
    non-allosteric to allosteric mechanism" in "The enumeration engine").
  * V-type: a base paired with a declared allosteric regulator at a new site, one variant
    per `(regulator, tag)` with `tag ∈ {:OnlyA, :OnlyI}`. On the every-chemistry base the
    inactive state binds as the active state does and runs no chemistry, so over one
    subunit, without a regulator, `L` folds into `kcat` (`v = kcat/(1+L)·shape`). Above one
    subunit, or on a base that stops one half-reaction, `L` can show without a regulator,
    yet the move pairs every base with one.

The move never emits the all-`:EqualAI` tagging: its two conformations are identical, `L`
cancels, and the law is `m`'s. The variants are distinct: each differs from the others in
its tags, its regulatory site or its oligomeric state. A parent that fails
`_hyperbolic_catalysis` gets variants at oligomeric state 1 only, and above one subunit the
move drops a variant whose conformations' laws fail it (`_requires_hyperbolic_catalysis`).
Errors on an allosteric `c`.
"""
function _expand_to_allosteric(c::CandidateMechanism, rxn::Mechanisms.EnzymeReaction)
    # The variants carry states and sites of their own; an allosteric parent's would be
    # lost.
    _is_allosteric(c.mechanism) && throw(ArgumentError(
        "_expand_to_allosteric: the mechanism is allosteric; the move makes a plain one " *
        "allosteric"))
    hyperbolic = _hyperbolic_catalysis(c)
    cns = [cn for cn in rxn.allowed_oligomeric_states
           if cn == 1 || hyperbolic]
    isempty(cns) && return CandidateMechanism[]
    chem = _chemistry_groups(c)
    vtags = [is_chem ? :OnlyA : :EqualAI for is_chem in chem]
    regs = [reg.name for reg in rxn.regulators if reg.kind === :allosteric]
    bases = _tagging_bases(c)
    kbases = vtags in bases ? bases : [bases; [vtags]]
    ktags = Vector{Symbol}[]
    for base in kbases, sel in _subsets(findall(!, chem))[2:end]
        tags = copy(base)
        tags[sel] .= :OnlyA
        push!(ktags, tags)
    end
    results = CandidateMechanism[]
    for cn in cns
        for tags in ktags
            push!(results, _candidate(c, tags, cn, Mechanisms.RegulatorySite[]))
        end
        # V-type: a base paired with a declared regulator, whose binding reveals L.
        any(chem) || continue
        for base in bases, reg in regs, tag in (:OnlyA, :OnlyI)
            push!(results,
                  _candidate(c, base, cn, [Mechanisms.RegulatorySite([reg], cn, [tag])]))
        end
    end
    filter(_hyperbolic_where_required, results)
end

"""
    _expand_add_allosteric_regulator(c::CandidateMechanism, rxn::EnzymeReaction)
        → Vector{CandidateMechanism}

The children of candidate `c`, whose mechanism is allosteric, that bind one more allosteric
regulator of `rxn`, one that no regulatory site of `c` binds, one child per site and tag:

  * `:OnlyA`, `:OnlyI` or `:NonequalAI` at a new site, or at an existing site that acts on
    a conformation the tag acts on. At an all-`:OnlyA` site an `:OnlyI` ligand, or the
    reverse, gives the rate equation of the new-site child.
  * `:EqualAI` only at an existing site that holds a ligand of another tag. A site whose
    ligands are all `:EqualAI` binds both conformations alike and cancels from the rate
    equation; this rule, not the `RegulatorySite` constructor, keeps the move from building
    one.

A name `c` binds as a dead-end competitive inhibitor stays eligible: a metabolite may play
both roles, and the two have distinct parameters. A new site takes the mechanism's
oligomeric state as its multiplicity, and each child keeps the rest of `c`.
"""
function _expand_add_allosteric_regulator(c::CandidateMechanism,
                                          rxn::Mechanisms.EnzymeReaction)
    am = c.mechanism
    taken = _bound_allo_regs(am)
    sites = am.allosteric.regulatory_sites
    results = CandidateMechanism[]
    for reg in rxn.regulators
        reg.kind === :allosteric && !(reg.name in taken) || continue
        for tag in (:OnlyA, :OnlyI, :NonequalAI), site_idx in 0:length(sites)
            site_idx >= 1 && isempty(intersect(_state_conformations(tag),
                _site_active_states(sites[site_idx]))) && continue
            push!(results, _make_am_with_added_reg(c, reg.name, tag, site_idx))
        end
        for site_idx in findall(site -> any(!=(:EqualAI), site.allosteric_states),
                                sites)
            push!(results, _make_am_with_added_reg(c, reg.name, :EqualAI, site_idx))
        end
    end
    results
end

"""
The candidate identical to `c`, whose mechanism is allosteric, except the ligand
`reg::Symbol` is added at `site_idx` (0 = create a new site, 1 to the number of sites =
append to that existing site) with allosteric state `tag`. A new site takes the
mechanism's oligomeric state as its multiplicity.
"""
function _make_am_with_added_reg(c::CandidateMechanism, reg::Symbol, tag::Symbol,
                                 site_idx::Int)
    a = c.mechanism.allosteric
    sites = copy(a.regulatory_sites)
    if site_idx == 0
        push!(sites, Mechanisms.RegulatorySite([reg], a.oligomeric_state, [tag]))
    else
        s = sites[site_idx]
        sites[site_idx] =
            Mechanisms.RegulatorySite([s.ligands; reg],
                                      s.multiplicity,
                                      [s.allosteric_states; tag])
    end
    _candidate(c, a.allosteric_states, a.oligomeric_state, sites)
end

"""
    _expand_change_allo_state(c::CandidateMechanism) → Vector{CandidateMechanism}

The relaxing move: the children of candidate `c`, whose mechanism is allosteric, that tag
`:NonequalAI` what `c` tags `:EqualAI`, `:OnlyA` or `:OnlyI`, keeping the rest of `c`:

  * one child per binding group, and one per regulatory ligand;
  * one child that relaxes every chemistry group (`_chemistry_groups`) at once, so a fully
    live inactive conformation is one move from every V-type, the tagging move's one-half
    bases included. The move makes it only while no binding group is `:OnlyA`. Beside an
    `:OnlyA` binding the restored chemistry can join the binding's two sides in the
    inactive conformation and leave its affinity no path to zero. Where the inactive
    conformation binds nothing (every binding `:OnlyA`), the skipped variant is
    rate-equivalent to the form emitted directly, so no observable hypothesis is lost.

Relaxing an `:OnlyA` binding brings its complex into the inactive state. There a copy's
complex that had no twin can gain one, and with a consistent gauge the copy's constant
shows in neither state. So the move drops a child with a redundant competitive-inhibitor
group (`_redundant_copy_groups`). An `:OnlyA` abortive binding, once relaxed, brings its
complex into the inactive law, which can then carry a power. So above one subunit the move
drops a child whose conformations' laws fail `_hyperbolic_catalysis`
(`_requires_hyperbolic_catalysis`). The exits (`_enumerator_keeps`) judge the rest.
Relaxing an `:OnlyA` binding joins its two sides, which can leave another tagged binding
inside one component or close a cycle of tagged bindings, and a face can fail. A
regulatory ligand's tag does not enter the guard, since a regulator site completes no
catalytic cycle.
"""
function _expand_change_allo_state(c::CandidateMechanism)
    a = c.mechanism.allosteric
    states = a.allosteric_states
    chem = _chemistry_groups(c)
    relaxations = [[g] for g in eachindex(states) if !chem[g] && states[g] != :NonequalAI]
    any(chem .& (states .!= :NonequalAI)) && !any(.!chem .& (states .== :OnlyA)) &&
        push!(relaxations, findall(chem))
    results = CandidateMechanism[]
    for gs in relaxations
        new_states = copy(states)
        new_states[gs] .= :NonequalAI
        push!(results, _candidate(c, new_states, a.oligomeric_state, a.regulatory_sites))
    end
    for (si, site) in enumerate(a.regulatory_sites),
        li in eachindex(site.ligands)
        site.allosteric_states[li] == :NonequalAI && continue
        new_sites = copy(a.regulatory_sites)
        new_states = copy(site.allosteric_states)
        new_states[li] = :NonequalAI
        new_sites[si] = Mechanisms.RegulatorySite(site.ligands,
                                                  site.multiplicity, new_states)
        push!(results, _candidate(c, states, a.oligomeric_state, new_sites))
    end
    filter!(child -> isempty(_redundant_copy_groups(child)), results)
    filter!(_hyperbolic_where_required, results)
end

"""
    _state_conformations(state::Symbol) -> Set{Symbol}

The conformations a single ligand's allosteric `state` acts on: `:active` for
`:OnlyA`/`:EqualAI`/`:NonequalAI`, `:inactive` for `:OnlyI`/`:EqualAI`/`:NonequalAI`.
"""
function _state_conformations(state::Symbol)
    conf = Set{Symbol}()
    state === :OnlyI || push!(conf, :active)
    state === :OnlyA || push!(conf, :inactive)
    conf
end

"""
    _site_active_states(site::RegulatorySite) -> Set{Symbol}

The conformations a regulatory site's ligands act on: the union of `_state_conformations`
over its states.
"""
function _site_active_states(site::Mechanisms.RegulatorySite)
    active = Set{Symbol}()
    for st in site.allosteric_states
        union!(active, _state_conformations(st))
    end
    active
end

"""
    _expand_merge_regulatory_sites(c::CandidateMechanism) → Vector{CandidateMechanism}

The merge move: for each pair of regulatory sites of candidate `c`, the children whose two
sites become one shared site that holds both sites' ligands, so that they compete. Every
ligand keeps its own dissociation constant, so a merge adds no parameter. The shared site
takes each of these tag assignments:

  * every ligand's own tag (co-binding), unless the two sites act on disjoint conformations
    (`_site_active_states`): an all-`:OnlyA` site merged with an all-`:OnlyI` site derives
    the rate equation of the two sites kept apart;
  * each ligand in turn retagged `:EqualAI`, an antagonist that binds both conformations
    equally and acts by competing for the site.

The move's own test drops the all-`:EqualAI` assignment, which cancels from the rate
equation and which the `RegulatorySite` constructor accepts. The shared site keeps the
first site's multiplicity (the oligomeric state), and the constructor sorts its ligands by
name, so two routes to one ligand partition give `==` mechanisms. `candidate_mechanisms`
drops a child that breaks a declared regulator type (`_respects_reg_type`).
"""
function _expand_merge_regulatory_sites(c::CandidateMechanism)
    a = c.mechanism.allosteric
    sites = a.regulatory_sites
    n = length(sites)
    results = CandidateMechanism[]
    for i in 1:(n - 1), j in (i + 1):n
        ligs = vcat(sites[i].ligands, sites[j].ligands)
        base_states = vcat(sites[i].allosteric_states,
                           sites[j].allosteric_states)
        mult = sites[i].multiplicity
        others = Mechanisms.RegulatorySite[sites[k] for k in 1:n if k != i && k != j]
        redundant = isempty(intersect(_site_active_states(sites[i]),
                                       _site_active_states(sites[j])))
        assignments = redundant ? Vector{Symbol}[] : [base_states]
        for k in eachindex(base_states)
            base_states[k] == :EqualAI && continue
            retagged = copy(base_states)
            retagged[k] = :EqualAI
            all(==(:EqualAI), retagged) || push!(assignments, retagged)
        end
        for states in assignments
            merged = Mechanisms.RegulatorySite(ligs, mult, states)
            push!(results, _candidate(c, a.allosteric_states, a.oligomeric_state,
                                      vcat(others, [merged])))
        end
    end
    results
end
