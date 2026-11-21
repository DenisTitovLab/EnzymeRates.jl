# ABOUTME: The expansion moves of the enumeration (RE→SS flips, group splits, dead-end
# ABOUTME: and allosteric regulators, allosteric states, site merges) and their searches.

# ─── Kinetic groups ──────────────────────────────────────────

"""Whether each kinetic group of graph `m` holds some step `k` with `pred(k)`."""
function _groups_with(pred, m::Mechanisms.CatalyticScheme)
    out = falses(maximum(e -> e.group, m.steps; init = 0))
    for k in eachindex(m.steps)
        pred(k) && (out[m.steps[k].group] = true)
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
    _candidate(cm.substrates, cm.products, c.forms, _relabelled(cm, groups, re),
               Mechanisms.Allostery(_label_tags(m, groups), a.oligomeric_state,
                                    a.regulatory_sites))
end

"""Candidate `c` with step `k` of its catalytic graph at rapid equilibrium where `re[k]`
holds, every step in its kinetic group (`_regrouped` with `c`'s own groups)."""
_with_flags(c::CandidateMechanism, re::BitVector) =
    _regrouped(c, Int[e.group for e in c.mechanism.catalytic.steps], re)

"""The steps of `m`, step `k` labelled `groups[k]` and at rapid equilibrium where `re[k]`
holds."""
_relabelled(m::Mechanisms.CatalyticScheme, groups::Vector{Int}, re::BitVector) =
    Mechanisms.Step[Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off, re[k], groups[k])
                    for (k, e) in enumerate(m.steps)]

"""The allosteric state of each label of `groups`, a kinetic group label per step of
`m`'s catalytic graph: the state of the group of the label's first step. A plain
mechanism's groups have no state, so neither do its labels."""
function _label_tags(m::Mechanisms.Mechanism, groups::Vector{Int})
    states = m.allosteric.allosteric_states
    isempty(states) && return Symbol[]
    Symbol[states[m.catalytic.steps[findfirst(==(l), groups)].group]
           for l in 1:maximum(groups; init = 0)]
end

# ─── Expansion Moves ─────────────────────────────────────────

"""
    _inhibitor_competition_patterns(sub_names, prod_names,
        existing_inhibitors)
        → Vector{Tuple{Set{Symbol},Set{Symbol},Set{Symbol}}}

Enumerate inhibitor competition patterns: (competing_subs,
competing_prods, competing_inhibitors). Each inhibitor must
compete with ≥1 substrate and ≥1 product. Competition with
existing inhibitors is a free binary choice.
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
    _expand_re_to_ss(c::CandidateMechanism) -> Vector{CandidateMechanism}

RE→SS expansion move. A flip unit is a whole kinetic group that is all-RE, takes up
or gives off no modifier (competitive-inhibitor binding stays at rapid equilibrium by
modeling choice), and holds a flux-carrying step (`_flux_carrying_groups` with every
edge steady state; a group with none exposes only equilibrium ratios under every
assignment and would gain a phantom parameter).
A flank of a qualifying chain whose isomerization is steady state is no unit
(`_chain_flank_groups`). One child is produced per minimal set of units whose joint flip
raises the RE segment count (`_minimal_flips`): a single group when it
cuts a segment on its own, several groups when each alone is bridged by an RE
route through the others — as happens once a split has separated a
metabolite's binding steps, or a catalytic step from its inhibitor-bound
mirror. A flip that leaves the segment count unchanged adds an SS step whose
endpoints share a segment, which the rate equation never sees. A flip set whose child
`_enumerator_keeps` refuses (the guard rejects it, or a conformation keeps a form that only
bindings it excludes enter) also counts as failed, so the minimal-set search
extends it instead of emitting it; so does
a set one of whose flipped groups carries no flux in the child
(`_flux_carrying_groups` on the flipped groups): the child is its parent plus
one phantom per such group, and its flux-carrying supersets are the smallest
identifiable relaxations. The non-degenerate mechanisms beyond a set that counts as failed
stay reachable — through the other cut orders of the same ring, or through the extended
set when no other order gains. Every other flag, the declared
names and, for an allosteric mechanism, the tags, multiplicity and regulatory sites are
kept; each child goes through the checked constructor.

A parent that `_requires_hyperbolic_catalysis` (a conformational mechanism over
more than one catalytic subunit) keeps only the children whose laws pass
`_hyperbolic_catalysis` (each conformation's, for an allosteric child). The check
runs on the emitted minimal sets, not inside the gain predicate: a failing set
would otherwise be extended with more flips, and more steady-state steps never
restore a hyperbolic equation, so its supersets need no visit.
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
other step, each a binding into that form and alone in its group. By the chain lemma the
mechanism sees such a chain only through its flux and enzyme content, four numbers that
the form with both flanks at rapid equilibrium already covers, so flipping a flank adds
parameters and no rate law. An allosteric mechanism has none: at more than one catalytic
subunit a flank's flip can be visible through the two conformations.
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
it builds. The first two tests read the flags, so only a candidate that passes them is
built. The sets come from `_minimal_gaining_sets`, so a failing set is extended and a
passing one is not.
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

"""Each metabolite name's net uptake sign, the names sorted (the rows of the
offsets of `rapid_equilibrium_segments`): +1 for a substrate, −1 for a product, 0 for a
modifier-only name."""
_name_signs(m::Mechanisms.CatalyticScheme) =
    Int[n in m.substrates ? 1 : n in m.products ? -1 : 0
        for n in sort!(Mechanisms.metabolites(m))]

"""The substrates step `k` of `m` takes up net, less the products it takes up net. A
reactant's modifier copy counts as its reactant, whose concentration it shares."""
function _uptake_weight(m::Mechanisms.CatalyticScheme, k::Int)
    t, g = m.steps[k].taken_up, m.steps[k].given_off
    net(x) = _copies(t, x) - _copies(g, x)
    sum(net, m.substrates; init = 0) -
        sum(net, m.products; init = 0)
end

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

"""Number of rapid-equilibrium segments of graph `m`, edge `k` at rapid equilibrium where
`re[k]` holds (`rapid_equilibrium_segments`). The flip move measures each candidate on its
flags, so it builds only the candidates that raise the count."""
_re_segment_count(m::Mechanisms.CatalyticScheme, re::BitVector) =
    length(Mechanisms.rapid_equilibrium_segments(m, trues(length(m.steps)), re).segments)

"""
    _minimal_gaining_sets(gains, partners) -> Vector{Vector{Int}}

Every minimal set of units for which `gains(set)` holds, found Apriori-style:
level 1 tests each unit of `partners(Int[])`; a level-`j` set is tested only if
every `(j−1)`-subset was tested and failed at the previous level, so no superset
of a gaining set is ever tested. `partners(set)` lists the units allowed to
extend `set` (units already in `set` are skipped). The loop ends when a level
fails nothing. Each returned set is sorted; sets are ordered by level, then
lexicographically.
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

Kinetic-group split move. A unit is one context bipartition of one group
(`_context_bipartitions`: the affinity for a metabolite may depend on which
other ligand is already bound). One child is produced per minimal set of units,
at most one per group, whose joint application raises the independent-parameter
count and gives a child `_enumerator_keeps` accepts (`_minimal_gaining_sets`). A single
bipartition whose new constant a Wegscheider cycle ties straight back is not a model:
the cycle ties the new constant to the others, so the child's rate function is the
parent's. Such a bipartition
is kept as a seed for pairs and larger sets, because the cycle that absorbs it
is broken by also splitting another group on that cycle — random-order binding
needs one split per substrate before any constant frees up. Partners for a
rejected set are the all-RE groups sharing a rapid-equilibrium segment with any
step of the split groups: a tie needs an RE cycle through the carve, and every
group on that cycle lies in that segment. The count is the joint solve's on the built
child (`fitted_parameter_count`), which derives no rate polynomial, and a set whose
child `_enumerator_keeps` refuses also counts as failed, so the search extends it instead
of emitting it.
A part of a steady-state group none of whose edges carries flux
(`_flux_carrying_edges`, computed once on the parent, since a split moves no edge)
is emitted at rapid equilibrium (`_revert_zero_flux_parts`), the same family with
one constant fewer; the count takes the reverted part under
its new kind, so a reverted constant the Wegscheider ties pull back is absorbed
like any tied split. A split of any group can complete a
copy's gauge, by separating a binding that forms or leaves a twin from the bindings
elsewhere in its group that blocked the gauge (the split then holds the family of the
copy its shared group pinned), and a split of a group that forms a copy's
twins can break the gauge, so the copy rule is decided on each child: a child with a
redundant copy group (`_redundant_copy_groups`) is not emitted; its family is that of
the same split without the copy.
The declared names and, for an allosteric mechanism, the multiplicity and regulatory
sites are kept; both parts of a split group inherit its allosteric state. Each
child is the parent's steps regrouped (`_apply_bipartitions`), through the checked
constructor, since a regrouping can turn a tied group's member.
"""
function _expand_split_kinetic_group(c::CandidateMechanism)
    cm = c.mechanism.catalytic
    members = _kinetic_groups(cm)
    flags = _rapid_equilibrium_flags(cm)
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

Ways to divide kinetic group `g` of the catalytic graph of candidate `c` in two by binding
context, each part a list of edges: another ligand already bound, the enzyme's conformation,
or its covalent residual. Encodes "the affinity for this metabolite may depend on what else
is bound, on the conformation, or on the covalent state" — the ligand family includes a
modifier, which separates a catalytic step from its inhibitor-bound mirror. A ligand is a
metabolite name bound productively or as a modifier copy (the form's record), so a reactant
and its modifier copy are two ligands. Each context value divides the edges whose context
form (`_context_form`) carries it from the rest; the group's own bound ligand is never a
context, and a group whose forms share one conformation and one residual gets no bipartition
from those two families. Both parts are nonempty, the first part holds the group's first
edge, each part lists its edges in edge order, context values that induce the same division
give one bipartition, and the order is ligands (modifiers, then products, then substrates,
each by name), then conformations (by label), then residuals (by the names they add, then
those they subtract), for deterministic output.
"""
function _context_bipartitions(c::CandidateMechanism, g::Int)
    m = c.mechanism.catalytic
    ks = [k for k in eachindex(m.steps) if m.steps[k].group == g]
    # The group's own bound ligand, as (name, copy), when its first edge is a binding.
    k1 = first(ks)
    own = _is_binding(m.steps[k1]) ?
        (n = only(m.steps[k1].taken_up).name; (n, _copy_change(c, k1, n) > 0)) : nothing
    forms = [_context_form(m, k) for k in ks]
    carries(x, (name, copy)) = name in (copy ? c.forms[x].copies : c.forms[x].productive)
    ligands = unique!([(n, copy) for x in forms
                       for n in [c.forms[x].productive; c.forms[x].copies]
                       for copy in (false, true)
                       if carries(x, (n, copy)) && (n, copy) != own])
    role((name, copy)) = copy ? 0 : name in m.products ? 1 : 2
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
allosteric state (splitting is a parameter-relaxation move that must not change A/I
semantics).
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
takes the metabolite up, its second when it gives it off, so a Theorell–Chance step E(A) + B
→ E(Q) + P frees B at E(A) and P at E(Q). A reactant and its modifier copy share a name but
bind different sites: the reactant asks where the metabolite binds productively, the copy
where it binds its dead-end site. The dead-end move targets these forms for competition with
the metabolite. A metabolite `c` does not hold is free nowhere.
"""
function _forms_where_free(c::CandidateMechanism, name::Symbol, copy::Bool)
    steps = c.mechanism.catalytic.steps
    free = falses(length(c.forms))
    for k in eachindex(steps)
        (; from, to, taken_up, given_off) = steps[k]
        # How many of `name` the step takes up and gives off in the role asked: `d` are
        # copies, the rest bind productively.
        d = _copy_change(c, k, name)
        up, off = copy ? (max(d, 0), max(-d, 0)) :
            (_copies(taken_up, name) - max(d, 0), _copies(given_off, name) - max(-d, 0))
        up > 0 && (free[from] = true)
        off > 0 && (free[to] = true)
    end
    free
end

"""
The records of the forms and the steps of the dead-end move's child of candidate `c`: the
records and steps of `c`'s catalytic graph, each step in its group, with modifier `reg`
also bound at each form of `sites`, a complex appended in the site's conformation and
residual; each step between two sites mirrored between their complexes, in the step's own
group; and one new group, after `c`'s, of the rapid-equilibrium bindings of `reg` at the
sites.
"""
function _dead_end_graph(c::CandidateMechanism, reg::Symbol, sites::Vector{Int})
    m = c.mechanism.catalytic
    inhibitor = Mechanisms.Metabolite(reg, true)
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
                                     e.given_off, e.rapid_equilibrium, e.group))
    end
    G = maximum(e -> e.group, m.steps; init = 0)
    for f in sites
        push!(steps, Mechanisms.Step(f, complex[f], [inhibitor], Mechanisms.Metabolite[],
                                     true, G + 1))
    end
    forms, steps
end

"""
The dead-end move's candidate child of the mechanism `m` whose forms have the records
`forms` and whose steps are `steps` (`_dead_end_graph`), binding modifier `reg` in a new
group; `nothing` when the copy is redundant in it (`_redundant_copy_groups`). A `Mechanism`
child's new group, the group of the steps that take the copy up, is tested by
`_gauge_rescaling` on the child's own twins (`_productive_twin`). The twins are the child's,
not the parent's, though they are the same forms: the complexes are not productive, and
their edges attach each complex to its site's segment and join no two segments, so every
form of `m` keeps its segment and offsets; but the complexes sort in among `m`'s forms, so
`m`'s form indices do not hold in the child. Only the new group is tested: `m` holds no
redundant copy group (the parent rule), and in each of about 218,000 enumerated parents
checked, the placement made no older copy redundant (an empirical finding, not a proof). An
allosteric child takes its copy's group tagged `:EqualAI`, with `m`'s multiplicity and
regulatory sites, and runs `_redundant_copy_groups` only when every complex of the copy has
a twin in the active state (`_all_twin`), which redundancy needs. Each child is built
through the checked constructor (`_candidate`).
"""
function _dead_end_child(m::Mechanisms.Mechanism, forms::Vector{_FormRecord},
                         steps::Vector{Mechanisms.Step}, reg::Symbol)
    a = m.allosteric
    states = _is_allosteric(m) ? [a.allosteric_states; :EqualAI] : Symbol[]
    child = _candidate(m.catalytic.substrates, m.catalytic.products, forms, steps,
                       Mechanisms.Allostery(states, a.oligomeric_state, a.regulatory_sites))
    cm = child.mechanism.catalytic
    twin = _productive_twin(child, trues(length(cm.steps)))
    # A plain child's one state takes the gauge test of the new group alone; an
    # allosteric child's two states take `_redundant_copy_groups`, run only when the
    # copy has a twin at each of its bindings, the only edges that take it up.
    if !_is_allosteric(m)
        binding = findfirst(k -> _copy_change(child, k, reg) > 0, eachindex(cm.steps))::Int
        return _gauge_rescaling(child, _kinetic_groups(cm), cm.steps[binding].group,
                                twin, _one_state_class) === nothing ? child : nothing
    end
    bindings = findall(k -> _copy_change(child, k, reg) > 0, eachindex(cm.steps))
    _all_twin(child, bindings, twin) && !isempty(_redundant_copy_groups(child)) ?
        nothing : child
end

"""
    _expand_add_dead_end_regulator(c::CandidateMechanism, rxn::EnzymeReaction)
        → Vector{CandidateMechanism}

Add a dead-end regulator binding step set. For each competitive inhibitor declared in `rxn`
but not yet bound by the steps of `c`'s mechanism, enumerate inhibitor competition patterns
(S × P × existing inhibitors); for each pattern, add RE binding steps to the forms where a
competing metabolite is free, unless the form already holds a competing metabolite or holds
every substrate or every product. Competition is decided per site, not per name. A substrate
or product declared as a competitive inhibitor is a copy that binds a dead-end site of its
own, and that site and the reactant's catalytic site are different sites by definition.
Competition with a substrate or product M targets the forms where M is free on a productive
step and excludes the forms that hold M productively; competition with an existing inhibitor
I targets the forms where I binds as a modifier (its binding steps' source forms) and
excludes the forms that hold I; the capacity test counts productive bindings only. A form
that carries only the copy of M is therefore a site for an inhibitor that competes with M.
Mirror steps inherit their catalytic counterpart's kinetic group. All new binding steps for
a single regulator share one fresh trailing kinetic group (one K_R parameter). Each child is
the parent's graph with the copy added (`_dead_end_graph`). A pattern is skipped when its
child binds the new copy redundantly (`_redundant_copy_groups`), as `_dead_end_child`
decides. Ligands of an allosteric mechanism's regulatory sites are not eligible: they belong
on a regulatory site.

The caller must pass the declared `rxn`: its regulators decide which inhibitors are
eligible.
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
`_kept_forms_keep_a_binding` read the multiplicity or the sites, so every child the tagging
move builds on these tags gets the same verdict at the exits, up to the limit test's random
draw, whose seed hashes the mechanism.
"""
_valid_tagging(c::CandidateMechanism, tags::Vector{Symbol}) =
    _enumerator_keeps(_candidate(c, tags, 1, Mechanisms.RegulatorySite[]).mechanism)

"""
    _stops_inactive_turnover(c::CandidateMechanism, tags::Vector{Symbol}) -> Bool

Whether the tags `tags` leave the inactive conformation of the plain mechanism `m` of
candidate `c` no catalytic cycle: no cycle of the steps it runs
(`steps_in_inactive_conformation`) turns the reaction over (`has_turnover_cycle`).
"""
function _stops_inactive_turnover(c::CandidateMechanism, tags::Vector{Symbol})
    am = _candidate(c, tags, 1, Mechanisms.RegulatorySite[]).mechanism
    !Mechanisms.has_turnover_cycle(am.catalytic,
                                   Mechanisms.steps_in_inactive_conformation(am))
end

"""
    _tagging_bases(c::CandidateMechanism) -> Vector{Vector{Symbol}}

The tag sets the tagging move (`_expand_to_allosteric`) builds on, one tag per kinetic group
of candidate `c`'s plain mechanism (spec section 7, "Allosteric ping-pong"). The first
candidate tags every chemistry group (`_chemistry_groups`) `:OnlyA` and every binding group
`:EqualAI`; when it passes `_valid_tagging` it is the only base. Otherwise the bases are the
smallest sets of chemistry groups, by inclusion, that stop every catalytic cycle of the
inactive conformation (`_stops_inactive_turnover`) and pass `_valid_tagging`, each tagged
`:OnlyA` with every other group `:EqualAI`; a set that fails leaves its supersets in the
search. On a ping-pong mechanism the first candidate cuts the covalent forms off in the
inactive conformation. When it fails, a half-reaction's chemistry groups can be a base
alone, the other half running in both conformations, and a half spread over several
chemistry groups can give several bases. The bases come by size, each size in `_subsets`
order.
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

Convert the plain mechanism `m` of candidate `c` into allosteric variants, each built on a
base of `_tagging_bases`: the tag set with every chemistry group `:OnlyA` and every binding
group `:EqualAI` when it passes `_valid_tagging`, else each smallest set of chemistry groups
that stops the inactive conformation's catalysis and passes. A chemistry group is one
holding a chemistry step (`_chemistry_groups`): an isomerization, a fused binding or a
Theorell–Chance step; every other group is a binding group. On a ping-pong mechanism,
tagging both half-reactions cuts the covalent forms off in the inactive conformation; when
that fails, a base can tag one half alone, the other half running in both conformations, and
a half spread over several chemistry groups can give several (spec section 7, "Allosteric
ping-pong"). Two kinds are not enumerated: the all-`:EqualAI` baseline and a V-type variant
without a regulator.

  * The all-`:EqualAI` baseline is never emitted — the two conformations
    are identical, `L` cancels, and the mechanism is indistinguishable
    from `m`.
  * K-type: a nonempty set of binding groups set `:OnlyA` on top of a base, and on top
    of the every-chemistry set when that set is no base, since a tagged binding can
    restore its limit; each is emitted bare, and the exits drop those
    `_enumerator_keeps` refuses. Over more than one catalytic subunit the conformational
    equilibrium enters the rate with the multiplicity's power. Over one subunit the inactive
    conformation is a dead-end branch of the free enzyme, and `L` is a phantom in some
    K-type variants and not in others. In every K-type variant of the rapid-equilibrium
    ordered bi-bi seed a rescaling of the active state's constants absorbs it, whatever
    that conformation binds. In the merged ordered bi-bi whose B group mixes a fused and
    a plain binding it is a phantom only in the variants whose inactive conformation
    binds nothing: no rescaling absorbs an A or Q binding kept there. Binding nothing is
    not enough: in a random-order merged bi-bi, whose A group takes A up at E and at
    E(B) under one rate constant, `L` shows in every K-type variant. identify keeps every
    K-type variant in the frontier, since its children can carry laws its parent's children
    lack. A variant that `_is_phantom` (`src/Derivation/key.jl`) flags is not fitted when
    another mechanism holds the fit of its law at `L = 0`: it reads that fit, with `L = 0`.
  * V-type: a base, with no `:OnlyA` binding, paired with a declared allosteric regulator at
    a new site, one variant per `(regulator, tag)` with `tag ∈ {:OnlyA, :OnlyI}`. On the
    every-chemistry base the inactive state binds substrate and product as the active state
    does and runs no chemistry, so over one subunit, without a regulator, `L` folds entirely
    into `kcat` (`v = kcat/(1+L)·shape`) and is not observable; above one subunit it can
    show. On a one-half base the inactive conformation runs the other half-reaction, and `L`
    shows without a regulator too. Yet every base is emitted only with one. A reaction with
    no declared allosteric regulators emits no V-type.

A parent with no declared regulator ends with no child at the exits when it has no binding
group, and so no K-type, or when `_enumerator_keeps` refuses every K-type.

For each value in `rxn`'s `allowed_oligomeric_states`, the move emits
the K-types, then the V-types, at that `oligomeric_state`. Every variant
composes `m`'s graph, which it shares, with its tags, multiplicity and sites. The
variants are distinct by construction: each differs from the others in its tags, its
regulatory site or its multiplicity.

A parent whose catalytic scheme fails `_hyperbolic_catalysis` (random-order
steady-state binding, or a substrate that traps a steady-state intermediate in
an abortive complex, whose own equation carries concentration powers) emits
variants at multiplicity 1 only: above one subunit the conformational
equilibrium would add a second source of powers
(`_requires_hyperbolic_catalysis`). Above one subunit a variant is also dropped when
either conformation's law fails the test: on a base that stops one half-reaction, a
ping-pong's abortive complexes can give the inactive law a power the parent lacks.
"""
function _expand_to_allosteric(c::CandidateMechanism, rxn::Mechanisms.EnzymeReaction)
    # The variants carry states and sites of their own; an allosteric parent's would be lost.
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
    # K-type: a nonempty set of binding groups :OnlyA on top of a base, and on top of the
    # every-chemistry set when it is no base; the exits drop those `_enumerator_keeps`
    # refuses.
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

Add one allosteric regulator declared in `rxn` not yet bound by the allosteric mechanism
`am` of candidate `c` at a regulatory site. A name already bound as a catalytic-step
dead-end competitive inhibitor is still eligible: a metabolite may play both roles, and the
two render to distinct parameters. For each (new ligand, target site, tag) combination, emit
a variant:

  * target site ∈ {new site} ∪ {existing sites}
  * tag ∈ {:OnlyA, :OnlyI, :NonequalAI} for any target site
  * tag = :EqualAI only at an existing site that already has at least
    one non-`:EqualAI` ligand. A site whose ligands are all `:EqualAI`
    binds both conformations alike and cancels from the rate equation;
    this guard, not the `RegulatorySite` constructor, keeps the move from
    building one.

New sites take the mechanism's oligomeric state as their multiplicity.
The mechanism's catalytic side and regulatory_sites' multiplicities pass
through unchanged.

The regulators declared in `rxn` decide which regulators are eligible.
"""
function _expand_add_allosteric_regulator(c::CandidateMechanism,
                                          rxn::Mechanisms.EnzymeReaction)
    am = c.mechanism
    _is_allosteric(am) || throw(ArgumentError(
        "_expand_add_allosteric_regulator: the mechanism is plain; the move adds a " *
        "regulator to an allosteric one"))
    taken = _bound_allo_regs(am)
    sites = am.allosteric.regulatory_sites
    results = CandidateMechanism[]
    for reg in rxn.regulators
        reg.kind === :allosteric && !(reg.name in taken) || continue
        # Non-:EqualAI tags at any (new or existing) site.
        for tag in (:OnlyA, :OnlyI, :NonequalAI), site_idx in 0:length(sites)
            # Appending to an existing site whose conformations are disjoint
            # from the new ligand's (an all-:OnlyA site gaining an :OnlyI
            # ligand, or the reverse) reproduces the add-at-a-new-site
            # equation — redundant, so skip it.
            site_idx >= 1 && isempty(intersect(_state_conformations(tag),
                _site_active_states(sites[site_idx]))) && continue
            push!(results, _make_am_with_added_reg(c, reg.name, tag, site_idx))
        end
        # :EqualAI only at an existing site that keeps a non-:EqualAI ligand: an
        # all-:EqualAI site binds both conformations alike and cancels from the rate.
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

Relax a "constrained" tag (`:EqualAI`, `:OnlyA`, `:OnlyI`) to `:NonequalAI`. A binding
catalytic group and a regulatory ligand each relax individually — one variant per group
not already `:NonequalAI`. The chemistry groups (those holding a chemistry step,
`_chemistry_groups`) relax **together**, one variant setting every non-`:NonequalAI`
chemistry group to `:NonequalAI` at once, so a fully live inactive conformation is one
move from every V-type, the tagging move's one-half bases included. The base catalytic
steps, multiplicity, and untouched tags are preserved.

The chemistry relaxation is tried only while no binding group is `:OnlyA`. Beside an
`:OnlyA` binding the restored chemistry can join the binding's two sides in the inactive
conformation and leave its affinity no path to zero; where the inactive conformation
binds nothing (every binding `:OnlyA`), the skipped variant is rate-equivalent to the
form emitted directly, so no observable hypothesis is lost.

The test at the enumerator's exits (`_enumerator_keeps`) judges each child:
relaxing an `:OnlyA` binding joins its two sides, which can leave another tagged binding
inside one component or close a cycle of tagged bindings, and a face can fail. A
regulatory ligand's tag does not enter the guard — a regulator site completes no catalytic
cycle. Above one subunit a child is kept only when each conformation's law passes
`_hyperbolic_catalysis` (`_requires_hyperbolic_catalysis`): relaxing an `:OnlyA` abortive
binding brings its complex into the inactive law, which can then carry a power.

A child with a redundant competitive-inhibitor group (`_redundant_copy_groups`)
is dropped. The per-state test reads the tags this move changes: relaxing an
`:OnlyA` binding to `:NonequalAI` brings its complex into the inactive state,
where a copy that was new there can become its twin with a consistent gauge, and
the copy's constant then shows in neither state.
"""
function _expand_change_allo_state(c::CandidateMechanism)
    _is_allosteric(c.mechanism) || throw(ArgumentError(
        "_expand_change_allo_state: the mechanism is plain; the move relaxes the " *
        "allosteric states of an allosteric one"))
    a = c.mechanism.allosteric
    states = a.allosteric_states
    chem = _chemistry_groups(c)
    # Binding catalytic groups relax individually; the chemistry groups relax together,
    # in one variant, and only while no binding group is :OnlyA.
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

The conformations a regulatory site's ligands act on — the union of
`_state_conformations` over its states. Two sites with disjoint active states —
an all-`:OnlyA` site and an all-`:OnlyI` site — merge to a rate equation
identical to keeping them separate, so that all-keep merge is redundant.
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

Merge each unordered pair of regulatory sites into one shared site holding
both sites' ligands, at no parameter cost (Δ0 — every ligand keeps its own
dissociation constant). This lets two regulators compete at one site,
including the activator↔antagonist ambiguity. For the merged ligands,
enumerate the Δ0-valid allosteric state assignments:

  * the all-keep assignment — every ligand retains its current state
    (co-binding), skipped when the two sites act on disjoint conformations
    (an all-`:OnlyA` site merged with an all-`:OnlyI` site), since that merge
    derives to the same rate equation as keeping the sites separate (see
    `_site_active_states`);
  * each assignment retagging exactly one ligand to `:EqualAI` — the
    antagonist forms (a ligand that binds both conformations equally,
    competing for the shared site).

The all-`:EqualAI` assignment (degenerate) is dropped by this move's own
guard; the `RegulatorySite` constructor does not reject it, so that guard is
load-bearing. The merged site reuses one site's `multiplicity` (equal to
`oligomeric_state`), and the `RegulatorySite` constructor sorts its
ligands by name, so two merge routes reaching the same ligand partition produce
`==` mechanisms and dedup by `hash`. Regulator type is not enforced here;
`candidate_mechanisms` drops every child that fails `_respects_reg_type`.
"""
function _expand_merge_regulatory_sites(c::CandidateMechanism)
    _is_allosteric(c.mechanism) || throw(ArgumentError(
        "_expand_merge_regulatory_sites: the mechanism is plain; the move merges the " *
        "regulatory sites of an allosteric one"))
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
