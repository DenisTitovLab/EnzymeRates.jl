# ABOUTME: Structural filters of the enumeration: degeneracy, hyperbolic catalysis, the
# ABOUTME: dwell gauge of redundant copies, atom and emission assertions, regulator types.

"""Whether the rate of graph `m`, edge `k` at rapid equilibrium where `re[k]` holds, has a
maximum as the metabolites `names` (its substrates or its products, sorted) grow: turning
rapid equilibrium every steady-state edge that takes up or gives off one of them leaves no
rapid-equilibrium turnover cycle (condition V). Only the edges `on` flags count, so an
allosteric mechanism's test can read its active law (`_degenerate`)."""
function _has_vmax(m::Mechanisms.Mechanism, re::BitVector, names::Vector{Symbol},
                   on::BitVector = trues(length(Mechanisms.steps(m))))
    touches = BitVector([_exchanges(m, k, names) for k in Mechanisms.steps(m)])
    !Mechanisms.has_turnover_cycle(m, (re .| touches) .& on)
end

"""Whether `s` holds a metabolite named in `names`, a reactant's modifier copy under its
reactant's name."""
_holds(s::Mechanisms.Metabolites, names::Vector{Symbol}) =
    any(n -> Mechanisms.copies(s, n) > 0, names)

"""Whether step `k` of `m` takes up or gives off a metabolite named in `names`
(`_holds`)."""
_exchanges(m::Mechanisms.Mechanism, k::Int, names::Vector{Symbol}) =
    _holds(Mechanisms.taken_up(m, k), names) || _holds(Mechanisms.given_off(m, k), names)

"""Whether step `k` of `m` takes up a metabolite named in one of `subs` and `prods` (a
reactant of one side of the reaction) and gives off one named in the other, as a
Theorell–Chance step does."""
function _crosses_sides(m::Mechanisms.Mechanism, k::Int, subs::Vector{Symbol},
                        prods::Vector{Symbol})
    t, g = Mechanisms.taken_up(m, k), Mechanisms.given_off(m, k)
    _holds(t, subs) && _holds(g, prods) || _holds(t, prods) && _holds(g, subs)
end

"""
Whether some chemistry node of graph `m`, edge `k` at rapid equilibrium where `re[k]`
holds, is left both by a rapid-equilibrium edge that releases a substrate and by one that
releases a product (condition C). A chemistry node is a merged complex (a form entered by
a fused binding of a substrate) with every form joined to it by RE isomerizations, a set
of forms joined to each other by RE isomerizations, or an RE Theorell–Chance step. An RE
edge leaves a node through a node form F and releases M when traversing it away from F
gives M off: F is its `to` form and M is consumed, or F is its `from` form and M is
released. An RE Theorell–Chance node leaves through its own edge in both directions.
Then the chemistry sits in rapid equilibrium with both sides and some reactant is not
needed at zero products or substrates. A reactant's modifier copy counts as its reactant.
Only the edges `on` flags count (`_degenerate`).
"""
function _chemistry_equilibrates_both_sides(
        m::Mechanisms.Mechanism, re::BitVector,
        on::BitVector = trues(length(Mechanisms.steps(m))))
    re = re .& on
    subs, prods = Mechanisms.substrates(m), Mechanisms.products(m)
    releases(node, side) = any(Mechanisms.steps(m)) do k
        from, to = Mechanisms.ends(m, k)
        re[k] &&
            (to in node && !(from in node) && _holds(Mechanisms.taken_up(m, k), side) ||
             from in node && !(to in node) && _holds(Mechanisms.given_off(m, k), side))
    end
    # The forms joined by RE isomerizations, one node per segment; a merged complex
    # outside every RE isomerization is a node of its own.
    isos = BitVector([re[k] && Mechanisms.step_kind(m, k) === :isomerization
                      for k in Mechanisms.steps(m)])
    s = Mechanisms.re_segments(m, isos, isos)
    nodes = [s.segments; [[x] for x in _merged_complexes(m, subs, on)
                          if s.segment_of[x] == 0]]
    any(node -> releases(node, subs) && releases(node, prods), nodes) ||
        any(k -> re[k] && _crosses_sides(m, k, subs, prods), Mechanisms.steps(m))
end

"""Whether `m` is degenerate by the structural conditions V and C: no maximal rate in one
direction (`_has_vmax`), or chemistry in rapid equilibrium with both sides
(`_chemistry_equilibrates_both_sides`). Read on the catalytic graph with its own flags,
over the active conformation's law (`conformation_steps`): every edge of a plain
mechanism, and the edges its tags keep of an allosteric one. The beam expands a degenerate
seed without fitting it (`_base_tier`)."""
function _degenerate(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism})
    cm = Mechanisms.catalytic(m)
    re = Mechanisms.rapid_equilibrium_steps(cm)
    on = Mechanisms.conformation_steps(m, :active)
    !(_has_vmax(cm, re, Mechanisms.substrates(cm), on) &&
      _has_vmax(cm, re, Mechanisms.products(cm), on)) ||
        _chemistry_equilibrates_both_sides(cm, re, on)
end

"""
Whether the enumerator may only give this mechanism a catalytic scheme that
passes `_hyperbolic_catalysis`. True when a conformational equilibrium sits over
a catalytic site with more than one subunit: the equilibrium then raises every
catalytic-site binding to the power of the multiplicity. Over one subunit it
only reweights each enzyme form and adds no power of its own.
"""
_requires_hyperbolic_catalysis(::Mechanisms.Mechanism) = false
_requires_hyperbolic_catalysis(am::Mechanisms.AllostericMechanism) =
    Mechanisms.catalytic_multiplicity(am) > 1

"""Whether `m` passes `_hyperbolic_catalysis` where it `_requires_hyperbolic_catalysis`;
true where it does not."""
_hyperbolic_where_required(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}) =
    !_requires_hyperbolic_catalysis(m) || _hyperbolic_catalysis(m)

"""
    _hyperbolic_catalysis(m::Mechanism, on = trues(length(steps(m)))) -> Bool

Whether the King–Altman denominator of `m`'s catalytic scheme has degree at most
1 in every substrate and product concentration. Every binding of a substrate or
product at its catalytic site counts, abortive complexes included. Edges touching
a form that holds a modifier are left out (a mask of the catalytic-site edges): a
modifier binds a site of its own by definition, so the powers its binding adds,
including those of a reactant's modifier copy, are a separate source that
conformational mechanisms keep. Only the edges `on` flags count, and `re_segments` leaves
out the forms they do not keep, so an allosteric mechanism's test can read each
conformation's law.

The equation is a sum over rapid-equilibrium (RE) segments and spanning
arborescences of the segment graph toward each segment. A denominator term is
the root segment's weight times the weight of every tree edge, so the exponent
of `X` in a term is the most `X` any form of the root segment carries beyond
the segment's bottom form, plus, per tree edge, one if the step binds `X` in
the tree direction and the count of `X` the edge's source form carries beyond
its bottom (`re_segments`). The degree exceeds 1 exactly when one
segment or one edge scores 2 or more, or the root scores 1 and some
arborescence toward it holds a scoring edge, or some arborescence holds two
scoring edges. An arborescence toward `S` containing given edges exists iff
every segment still reaches `S` once each given edge's source keeps that edge
as its only way out (`_all_reach`).
"""
function _hyperbolic_catalysis(m::Mechanisms.Mechanism,
                               on::BitVector = trues(length(Mechanisms.steps(m))))
    on_site = on .& BitVector([!any(_holds_modifier(m, x) for x in Mechanisms.ends(m, k))
                               for k in Mechanisms.steps(m)])
    s = Mechanisms.segment_graph(m, on_site)
    # A metabolite name's row of the offsets: its position among the sorted names.
    names = sort!(Mechanisms.metabolites(m))
    name_of(n) = searchsortedfirst(names, n)
    # Directed segment-graph edges: source segment, target segment, source form, the
    # names (their rows, `name_of`) bound in that direction.
    edges = Tuple{Int, Int, Int, Vector{Int}}[]
    rows(side) = [name_of(n) for n in Mechanisms.metabolite_names(side)]
    for k in s.steps
        from, to = Mechanisms.ends(m, k)
        a, b = s.segment_of[from], s.segment_of[to]
        push!(edges, (a, b, from, rows(Mechanisms.taken_up(m, k))))
        push!(edges, (b, a, to, rows(Mechanisms.given_off(m, k))))
    end
    n = length(s.segments)
    for y in [Mechanisms.substrates(m); Mechanisms.products(m)]
        x = name_of(y)
        score(e) = count(==(x), e[4]) + s.offsets[x, e[3]]
        carrying = [e for e in edges if score(e) > 0]
        any(e -> score(e) > 1, carrying) && return false
        root_score(k) = maximum((s.offsets[x, f] for f in s.segments[k]); init=0)
        any(k -> root_score(k) > 1, 1:n) && return false
        roots = [k for k in 1:n if root_score(k) == 1]
        for e in carrying, k in roots
            k != e[1] && _all_reach(n, edges, k, (e,)) && return false
        end
        for (p, e1) in enumerate(carrying), e2 in carrying[p + 1:end]
            e1[1] == e2[1] && continue
            any(k -> k != e1[1] && k != e2[1] && _all_reach(n, edges, k, (e1, e2)),
                1:n) && return false
        end
    end
    true
end

"""
    _hyperbolic_catalysis(am::AllostericMechanism) -> Bool

Whether each conformation's law of `am` (`conformation_steps`) passes
`_hyperbolic_catalysis`. A conformation that runs part
of the catalytic graph can carry a power the whole graph lacks: a ping-pong's abortive
complexes can put A² into the inactive law once that conformation stops one
half-reaction.
"""
_hyperbolic_catalysis(am::Mechanisms.AllostericMechanism) =
    all(c -> _hyperbolic_catalysis(Mechanisms.catalytic(am),
                                   Mechanisms.conformation_steps(am, c)),
        (:active, :inactive))

"""
Whether every segment of the segment graph reaches `root` when each edge in
`fixed` is its source segment's only way out. A digraph has a spanning
arborescence toward `root` iff every vertex reaches `root`, and with the fixed
edges as their sources' only exits every such arborescence contains them.
"""
function _all_reach(n::Int, edges, root::Int, fixed)
    pinned = Dict(e[1] => e[2] for e in fixed)
    into = [Int[] for _ in 1:n]
    for (u, v, _, _) in edges
        get(pinned, u, v) == v && push!(into[v], u)
    end
    seen = falses(n)
    seen[root] = true
    queue = [root]
    while !isempty(queue)
        v = popfirst!(queue)
        for u in into[v]
            seen[u] || (seen[u] = true; push!(queue, u))
        end
    end
    all(seen)
end

"""
A class of the dwell gauge (`_gauge_rescaling`): a kinetic group's class (`group > 0`)
or a twin form's class of its own (`form > 0`), in one conformational state (`state` 1,
or 2 for the inactive one) or in both (0). Group and form indices count apart, so group
2 and form 2 are two classes.
"""
struct _GaugeClass
    group::Int
    form::Int
    state::Int
end

"""The class of the factor 1, which every form outside a copy's complexes and twins
takes."""
const _NO_RESCALING = _GaugeClass(0, 0, 0)

"""The gauge class of group `h` or form `t` (`_GaugeClass`) in a mechanism with one
conformational state."""
_one_state_class(h::Int, t::Int) = _GaugeClass(h, t, 1)

"""
    _productive_twin(m, mask) -> (site, ligand) -> Int

The productive form of the sub-graph of `m` whose edges `mask` keeps, a form that holds
no modifier, whose weight is proportional to that of the complex a modifier named
`ligand` forms at form `site` by a rapid-equilibrium step, or 0 when none is: a form of
the site's rapid-equilibrium segment with the complex's offsets (`re_segments`; the
modifier counts as its name, so a reactant's copy as the reactant), the one of the
complex's composition when several are, else the first in form order. Within a segment
rapid equilibrium fixes each weight relative to another as a constant times the
metabolites bound along the route between them, so equal offsets mean proportional
weights, and steady-state data see only their sum. The offsets reach a copy across a
rapid-equilibrium isomerization (a ping-pong second chemistry step, a conformational
isomer), where no form shares the complex's composition. A form of the same composition
in another segment, formed by a steady-state binding, is no twin: its weight is not a
constant multiple of the complex's, merging the two is no reparameterization, and the
copy's constant may stay visible.

Only productive forms are twin sources, so the complex itself never is one. A complex
that duplicates only a copy's complex shares its weight with that copy's constant, and
the sites that pin either copy may keep both constants separable, so such a match never
makes a twin. The forms are `m`'s own, so the twins of a mechanism the dead-end move
builds are found on that mechanism (`_dead_end_child`).
"""
function _productive_twin(m::Mechanisms.Mechanism, mask::BitVector)
    s = Mechanisms.re_segments(m, mask, Mechanisms.rapid_equilibrium_steps(m))
    # A metabolite name's row of the offsets: its position among the sorted names.
    names = sort!(Mechanisms.metabolites(m))
    # A form's bound names, sorted, its conformation and its residual.
    composition(f) = (Mechanisms.metabolite_names(Mechanisms.bound(m, f)),
                      Mechanisms.conformation(m, f), Mechanisms.residual(m, f))
    function (site::Int, ligand::Symbol)
        n = searchsortedfirst(names, ligand)
        n <= length(names) && names[n] === ligand || return 0
        twins = [f for f in s.segments[s.segment_of[site]]
                 if !_holds_modifier(m, f) &&
                    all(x -> s.offsets[x, f] == s.offsets[x, site] + (x == n),
                        axes(s.offsets, 1))]
        isempty(twins) && return 0
        bound_names, label, residual = composition(site)
        wanted = (sort!(push!(bound_names, ligand)), label, residual)
        k = findfirst(t -> composition(t) == wanted, twins)
        k === nothing ? minimum(twins) : twins[k]
    end
end

"""Whether form `x` of `m` holds a modifier (`modifier_copies`)."""
_holds_modifier(m::Mechanisms.Mechanism, x::Int) =
    !isempty(Mechanisms.modifier_copies(Mechanisms.bound(m, x)))

"""The name of the modifier copy that step `k` of `m` binds (`is_binding`), or `nothing`
when the step binds no modifier copy."""
function _bound_copy(m::Mechanisms.Mechanism, k::Int)
    t = Mechanisms.taken_up(m, k)
    Mechanisms.is_binding(m, k) && !isempty(Mechanisms.modifier_copies(t)) ?
        only(Mechanisms.metabolite_names(t)) : nothing
end

"""
Whether the edges `ks` of one kinetic group of `m` bind a modifier, and only at sites
whose complex has a productive twin (`twin`, a `_productive_twin` of one state's graph).
"""
function _all_twin(m::Mechanisms.Mechanism, ks::Vector{Int}, twin)
    ligand = _bound_copy(m, first(ks))
    ligand !== nothing && all(k -> twin(Mechanisms.ends(m, k)[1], ligand) != 0, ks)
end

"""
    _gauge_rescaling(m, members, g, twin, label)
        -> Union{Dict{Int, Tuple{_GaugeClass, _GaugeClass}}, Nothing}

The dwell gauge of the modifier copy bound by kinetic group `g` of `m`, in one
conformational state's graph, whose groups hold the edges `members`: the rescaling each
nonempty kinetic group needs, or `nothing` when the gauge does not exist. It does not
exist when group `g` binds no modifier or a complex it forms has no productive twin
(`_all_twin`, with `twin` a `_productive_twin` of that graph), or when some group would
need two rescalings. The gauge gives every form a factor σ: ρ_T to each twin T, one factor
s to every complex of the copy and 1 to every other form. Multiplying a rapid-equilibrium
group's K by σ(from)/σ(to), and a steady-state group's forward constant by σ(from) and its
reverse by σ(to), divides each form's weight by its factor and leaves every flux as it
was. A twin's weight is a constant multiple of its complexes', so the factors can keep
each twin's total with its complexes, and the rate law, unchanged while K* moves: the
copy's constant is a phantom. A kinetic group shares its constants, so every group, the
copy's own included, must have one ratio σ(from)/σ(to) over its edges (rapid
equilibrium), or one σ(from) and one σ(to) (steady state). The copy's group, at 1 on each
site, holds all its complexes at the one factor s, and a complex never takes its twin's
factor: its weight moves opposite to the twin's. A rapid-equilibrium mirror edge, between
two complexes, therefore has ratio 1, as its parent edge between two sites has, whatever
the twins; a steady-state mirror needs s beside its parent edge's 1, and the gauge fails.
When each site binds the copied ligand at rapid equilibrium, the gauge is the merge of
each complex into its twin, the twin's binding group taking K_h + K* in place of K_h, and
the family is that of the mechanism without the copy.

The ρ are generic numbers, equal only where the structure makes them equal: twins
formed from their sites by RE bindings of the copied ligand in one kinetic group h, each
the twin of one complex, share one factor (the merge gives it as ρ = K_h/(K_h + K*));
every other twin has its own. `label(h, t)` names the classes (`_GaugeClass`): it
receives the group of a shared class (h) or of the copy's complexes (g) with `t = 0`, or
`h = 0` and the twin `t` of a class of its own, so the caller decides which classes are
one number. A group's rescaling is the pair of its ends' classes, `_NO_RESCALING` standing
for the factor 1, and two of them for a rapid-equilibrium group whose two ends scale
alike. The test is bookkeeping over the edges, with no parameters and no numerics.
"""
function _gauge_rescaling(m::Mechanisms.Mechanism, members::Vector{Vector{Int}}, g::Int,
                          twin, label)
    _all_twin(m, members[g], twin) || return nothing
    ligand = _bound_copy(m, first(members[g]))::Symbol
    twin_of = Dict(Mechanisms.ends(m, k)[2] => twin(Mechanisms.ends(m, k)[1], ligand)
                   for k in members[g])
    # The rapid-equilibrium bindings of the ligand itself, not of a modifier copy.
    binding_group = Dict{Tuple{Int, Int}, Int}()
    for (h, ks) in enumerate(members), k in ks
        t = Mechanisms.taken_up(m, k)
        Mechanisms.is_rapid_equilibrium(m, k) && Mechanisms.is_binding(m, k) &&
            isempty(Mechanisms.modifier_copies(t)) && Mechanisms.copies(t, ligand) == 1 &&
            (binding_group[Mechanisms.ends(m, k)] = h)
    end
    class = Dict{Int, _GaugeClass}(complex => label(g, 0) for complex in keys(twin_of))
    for k in members[g]
        from, to = Mechanisms.ends(m, k)
        t = twin_of[to]
        haskey(class, t) && continue
        h = count(==(t), values(twin_of)) == 1 ? get(binding_group, (from, t), 0) : 0
        class[t] = h > 0 ? label(h, 0) : label(0, t)
    end
    σ(f) = get(class, f, _NO_RESCALING)
    rescale(k) = (a = σ(Mechanisms.ends(m, k)[1]); b = σ(Mechanisms.ends(m, k)[2]);
                  Mechanisms.is_rapid_equilibrium(m, k) && a == b ?
                      (_NO_RESCALING, _NO_RESCALING) : (a, b))
    rescaling = Dict{Int, Tuple{_GaugeClass, _GaugeClass}}()
    for (h, ks) in enumerate(members)
        isempty(ks) && continue
        r = rescale(first(ks))
        all(k -> rescale(k) == r, ks) || return nothing
        rescaling[h] = r
    end
    rescaling
end

"""
    _redundant_copy_groups(m) -> Vector{Int}

Kinetic groups of `m` that bind a modifier redundantly: in every conformational state
where the copy binds, every complex has a productive twin, a form whose weight is
proportional to the complex's (`_productive_twin`, `_all_twin`), and the dwell gauge
exists over the states together (`_gauge_rescaling`). The gauge then moves the copy's
constant along a direction the rate law cannot see, so the copy adds a phantom, and the
mechanism's family is its family without the copy (`_gauge_rescaling`). A copy group
whose complexes all have twins but whose gauge fails is not redundant: a group that forms
or leaves a twin also binds where no copy does (a shared group that pins the existing
binding), a group holds a step that leaves a complex beside one that leaves another form
(a second copy bound at the complex and at a twin), or a mirror is steady-state. The gauge
is the proof that the copy's constant is invisible, and without it the constant may be
identifiable. Two copies that are twins only of each other, neither pinned, pass the test
and carry a phantom: the rule admits them.

An allosteric copy binds the active state unless its tag is `:OnlyI` and the inactive
state unless it is `:OnlyA`. Each state's graph is its law's (`conformation_steps`): the
edges its tags keep between forms they join to free E. The copy binds in a state only at
the sites that graph keeps; in a state where it binds nothing, every factor of that state
is 1. A copy that binds in neither state is not redundant: its
constant enters neither law and is not fitted. The two states form one system.
A group tagged `:EqualAI` has one set of constants in both, so it must take the same
rescaling in each. A shared class is one number in both states only when its binding
group and the copy are both `:EqualAI` (one K_h, one K*). The copy's complexes share one
factor in both states when the copy is `:EqualAI` (one K*). Every other class is a number
per state.
"""
function _redundant_copy_groups(m::Union{Mechanisms.Mechanism,
                                          Mechanisms.AllostericMechanism})
    cm = Mechanisms.catalytic(m)
    members = Mechanisms._group_members(cm)
    any(ks -> _bound_copy(cm, first(ks)) !== nothing, members) || return Int[]
    if m isa Mechanisms.Mechanism
        twin = _productive_twin(cm, trues(length(Mechanisms.steps(cm))))
        return [g for g in eachindex(members)
                if _gauge_rescaling(cm, members, g, twin, _one_state_class) !== nothing]
    end
    am = m::Mechanisms.AllostericMechanism
    equal(h) = Mechanisms.tag(am, h) === :EqualAI
    masks = (Mechanisms.conformation_steps(am, :active),
             Mechanisms.conformation_steps(am, :inactive))
    kept = [[filter(k -> mask[k], ks) for ks in members] for mask in masks]
    twins = [_productive_twin(cm, mask) for mask in masks]
    unscaled(c) = Dict(h => (_NO_RESCALING, _NO_RESCALING)
                       for h in eachindex(kept[c]) if !isempty(kept[c][h]))
    filter(collect(eachindex(members))) do g
        all(c -> isempty(kept[c][g]), 1:2) && return false
        one_number(h) = h > 0 && equal(h) && equal(g)
        label(state) = (h, t) -> _GaugeClass(h, t, one_number(h) ? 0 : state)
        rescaling(c) = isempty(kept[c][g]) ? unscaled(c) :
                       _gauge_rescaling(cm, kept[c], g, twins[c], label(c))
        active = rescaling(1)
        active === nothing && return false
        inactive = rescaling(2)
        inactive === nothing && return false
        all(h -> !equal(h) || !haskey(inactive, h) || active[h] == inactive[h],
            keys(active))
    end
end

# ─── Atom-conservation validation ────────────────────────────

"""
Net atom multiset of the reactants named in `plus` minus those named in `minus`, with
zero counts dropped. A name may repeat; a name that names no reactant carries no atoms.
"""
function _net_atoms(reaction::Mechanisms.EnzymeReaction, plus, minus)
    acc = Dict{Symbol,Int}()
    for ra in Mechanisms.reactants(reaction)
        k = count(==(Mechanisms.name(Mechanisms.metabolite(ra))), plus) -
            count(==(Mechanisms.name(Mechanisms.metabolite(ra))), minus)
        for (a, c) in Mechanisms.atoms(ra)
            acc[a] = get(acc, a, 0) + k * c
        end
    end
    filter!(kv -> kv.second != 0, acc)
end

"""
    _assert_emission_rules(m)

The two rules every mechanism the moves emit satisfies, checked on a parent before
it is expanded. Every steady-state kinetic group holds an edge that carries net flux
(`_flux_carrying_groups` on the catalytic graph, an allosteric mechanism's active
state). No kinetic group binds a competitive inhibitor redundantly
(`_redundant_copy_groups`). A parent must obey both because a flip tests only the groups
it flips. The split, the dead-end move, and `_expand_change_allo_state` filter their
children; the other moves preserve both rules.
"""
function _assert_emission_rules(m::Union{Mechanisms.Mechanism,
                                          Mechanisms.AllostericMechanism})
    cm = Mechanisms.catalytic(m)
    members = Mechanisms._group_members(cm)
    # A group as its edges' sides (`_sides`), the names its constants take.
    label(g) = join((join(Mechanisms._sides(cm, k), " → ") for k in members[g]), ", ")
    flux = _flux_carrying_groups(cm, Mechanisms.rapid_equilibrium_steps(cm))
    for (g, ks) in enumerate(members)
        Mechanisms.is_rapid_equilibrium(cm, first(ks)) || flux[g] || error(
            "expand_mechanisms: steady-state kinetic group {" * label(g) * "} has no " *
            "step that carries net flux, so its two constants enter the rate only as " *
            "their ratio; write the group at rapid equilibrium")
    end
    for g in _redundant_copy_groups(m)
        error("expand_mechanisms: kinetic group {" * label(g) * "} binds a competitive " *
              "inhibitor only where the complex duplicates a productive form, and a " *
              "dwell gauge absorbs its constant into the existing binding's; bind the " *
              "inhibitor where it forms a new complex, or drop it")
    end
    nothing
end

"""
    _assert_atom_conserving(m, rxn::EnzymeReaction)

Assert every edge of an enumerated mechanism `m`'s catalytic graph conserves the atoms
that reaction `rxn` gives its metabolites: the atoms of its source form plus those it
takes up equal the atoms of its target form plus those it gives off, a form carrying the
atoms of its bound metabolites and of its residual's added substrates, less those of its
residual's subtracted products, each read from the reaction's inventory by name
(`_net_atoms`). Errors on the first edge that does not, group by group, naming its forms
and the metabolites it exchanges. This is the enumeration-path guardrail against
atom-non-conserving covalent intermediates; it is intentionally NOT a `Mechanism`
constructor check, since a hand-written `@enzyme_mechanism` block carries no atoms and
may fold steps.
"""
function _assert_atom_conserving(m::Union{Mechanisms.Mechanism,
                                           Mechanisms.AllostericMechanism},
                                 rxn::Mechanisms.EnzymeReaction)
    cm = Mechanisms.catalytic(m)
    names = Mechanisms.metabolite_names
    # The names counted on one side of an edge: `form`'s bound metabolites and residual
    # additions, the residual subtractions of the `other` form, and the free metabolites.
    side(form::Int, other::Int, free::Mechanisms.Metabolites) =
        [names(Mechanisms.bound(cm, form)); names(Mechanisms.residual(cm, form)[1]);
         names(Mechanisms.residual(cm, other)[2]); names(free)]
    for k in Mechanisms._flat_edges(cm)
        f, t = Mechanisms.ends(cm, k)
        taken, given = Mechanisms.taken_up(cm, k), Mechanisms.given_off(cm, k)
        diff = _net_atoms(rxn, side(t, f, given), side(f, t, taken))
        isempty(diff) && continue
        from, to = Mechanisms._form_name(cm, f), Mechanisms._form_name(cm, t)
        error("atom-non-conserving step $from → $to (consumed $(names(taken)), " *
              "released $(names(given))): " *
              "atoms(to) + atoms(released) − atoms(from) − atoms(consumed) = $diff")
    end
    nothing
end

# ─── Regulator-Type Filter ────────────────────────────────────

"""
    _respects_reg_type(m, rxn::EnzymeReaction) -> Bool

Whether every regulatory ligand of `m` respects the type `rxn` declares for it, given the
other ligands at its site. A `Mechanism` has no regulatory sites and passes. A ligand's
type is the `reg_type` of its `AllostericRegulator` entry in `rxn`, or `:unspecified`
when `rxn` has no such entry (the name is absent, or declared only as another
`Regulator` subtype). An `:unspecified` ligand always passes. A designated activator is
never `:OnlyI` and a designated inhibitor is never `:OnlyA`; the type-matching pure
state and `:NonequalAI` are always allowed. `:EqualAI` is an antagonist state and is
rejected only when another ligand at the site is a same-type designated effector, since
it would counteract that ligand's declared direction.
"""
_respects_reg_type(::Mechanisms.Mechanism, ::Mechanisms.EnzymeReaction) = true
function _respects_reg_type(am::Mechanisms.AllostericMechanism,
                            rxn::Mechanisms.EnzymeReaction)
    declared = Dict(Mechanisms.name(Mechanisms.regulator(rm)) => Mechanisms.reg_type(rm)
                    for rm in Mechanisms.regulators(rxn)
                    if Mechanisms.regulator(rm) isa Mechanisms.AllostericRegulator)
    all(Mechanisms.regulatory_sites(am)) do site
        types = [get(declared, l, :unspecified) for l in Mechanisms.ligands(site)]
        all(zip(types, Mechanisms.allo_states(site))) do (t, state)
            t === :unspecified ||
                !(t === :activator && state === :OnlyI ||
                  t === :inhibitor && state === :OnlyA ||
                  state === :EqualAI && count(==(t), types) > 1)
        end
    end
end
