# ABOUTME: Mechanism enumeration by incremental parameter count growth
# ABOUTME: Provides init_mechanisms, expand_mechanisms building blocks


# ─── Catalytic topologies ─────────────────────────────────────────────

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

# ─── Atom-conservation validation ────────────────────────────

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
    flux = _flux_carrying_groups(cm, Mechanisms._flags(cm))
    for (g, ks) in enumerate(members)
        cm.edges[first(ks)].rapid_equilibrium || flux[g] || error(
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
    names(c::Mechanisms.Counts) = Mechanisms._entry_names(cm, c)
    # The names counted on one side of an edge: `form`'s bound metabolites and residual
    # additions, the residual subtractions of the `other` form, and the free metabolites.
    side(form::Mechanisms.Form, other::Mechanisms.Form, free::Mechanisms.Counts) =
        [names(form.bound); names(form.added); names(other.subtracted); names(free)]
    for k in Mechanisms._flat_edges(cm)
        e = cm.edges[k]
        f, t = cm.forms[e.from], cm.forms[e.to]
        diff = _net_atoms(rxn, side(t, f, e.released), side(f, t, e.consumed))
        isempty(diff) && continue
        from, to = Mechanisms._form_name(cm, e.from), Mechanisms._form_name(cm, e.to)
        error("atom-non-conserving step $from → $to (consumed $(names(e.consumed)), " *
              "released $(names(e.released))): " *
              "atoms(to) + atoms(released) − atoms(from) − atoms(consumed) = $diff")
    end
    nothing
end

"""All subsets of `v` in binary-counting order (bit `i - 1` selects `v[i]`), empty first."""
_subsets(v::AbstractVector) =
    [v[[isodd(mask >> (i - 1)) for i in eachindex(v)]] for mask in 0:(1 << length(v)) - 1]

"""
All weak orderings of `items`: sequences of nonempty levels that partition `items`, each
level in the order of `items`. They come in the order of their first level in
`_subsets(items)`, then recursively in the order of the rest.
"""
function _weak_orderings(items::Vector{Symbol})
    isempty(items) && return [Vector{Symbol}[]]
    Vector{Vector{Symbol}}[[[level]; rest] for level in _subsets(items)[2:end]
                           for rest in _weak_orderings(setdiff(items, level))]
end

"""
    _catalytic_topologies(reaction) -> Vector{Vector{RawStep}}

Build catalytic cycle topologies by constructive backtracking. Each topology is a set of
raw steps forming one or more complete catalytic cycles (E -> ... -> E), every release
written as the binding it reverses, so every binding takes its metabolite up. Errors
when the reaction has none.
"""
function _catalytic_topologies(reaction::Mechanisms.EnzymeReaction)
    sub_names = Symbol[Mechanisms.name(s) for s in Mechanisms.substrates(reaction)]
    prod_names = Symbol[Mechanisms.name(p) for p in Mechanisms.products(reaction)]

    # The form on :E binding `on_subs` and `on_prods` after a route consumed `consumed`
    # and released `released`. Its covalent residual adds the consumed substrates not
    # bound and subtracts the released and bound products; it has none exactly when
    # their atoms cancel (no covalent residue remains).
    function form(consumed, released, on_subs, on_prods)
        added, subtracted = setdiff(consumed, on_subs), [released; on_prods]
        bound = Mechanisms.Met[Mechanisms.Met.(on_subs, false);
                               Mechanisms.Met.(on_prods, false)]
        isempty(_net_atoms(reaction, added, subtracted)) ? Mechanisms.RawForm(:E, bound) :
            Mechanisms.RawForm(:E, bound, added, subtracted)
    end

    # Collect all complete catalytic paths as raw step lists. A route starts at free E,
    # binds substrates, isomerizes a substrate-bound form to a product-bound one, releases
    # those products one at a time in every order, and repeats until it has consumed
    # every substrate and released every product. `cur` is the current form, `consumed`
    # and `released` are the route's history, and `on_subs` and `on_prods` are what `cur`
    # binds; the atoms on the enzyme are those of `consumed` minus those of `released`.
    paths = Vector{Mechanisms.RawStep}[]
    function walk!(cur, consumed, released, on_subs, on_prods, path)
        if isempty(cur.bound) && length(consumed) == length(sub_names) &&
                length(released) == length(prod_names)
            push!(paths, copy(path))
            return
        end
        # The step from `cur` to the form the arguments describe, taking up `takes_up`,
        # or, for a release, the binding into `cur` it reverses; then the walk on from
        # that form.
        function step!(takes_up::Vector{Mechanisms.Met}, release::Bool, con, rel, on_s,
                       on_p)
            to = form(con, rel, on_s, on_p)
            push!(path,
                  release ? Mechanisms.RawStep(to, cur, takes_up, Mechanisms.Met[], true) :
                            Mechanisms.RawStep(cur, to, takes_up, Mechanisms.Met[], true))
            walk!(to, con, rel, on_s, on_p, path)
            pop!(path)
        end
        if !isempty(on_prods)
            # Release a bound product; nothing binds until every bound product is gone.
            for p in on_prods
                step!([Mechanisms.Met(p, false)], true, consumed, [released; p], Symbol[],
                      filter(!=(p), on_prods))
            end
            return
        end
        remaining_subs = setdiff(sub_names, consumed)
        for s in remaining_subs
            step!([Mechanisms.Met(s, false)], false, [consumed; s], released, [on_subs; s],
                  Symbol[])
        end
        # C7: only a substrate-bound form isomerizes. C6: at most three substrates react.
        (isempty(on_subs) || length(on_subs) > 3) && return
        for subset in _subsets(setdiff(prod_names, released))[2:end]
            residue = _net_atoms(reaction, consumed, [released; subset])
            any(<(0), values(residue)) && continue
            # Admissible-residual rule: an isomerization leaves a covalent residue
            # exactly when substrates remain to bind. Without a residue while substrates
            # remain, the enzyme would return to apo E mid-cycle, splitting the reaction
            # into disconnected half-cycles. With a residue after the last substrate,
            # the route strands it: only a substrate-bound form isomerizes (C7).
            isempty(residue) == isempty(remaining_subs) || continue
            # C6: iso size limit, the residue counting as one product
            length(subset) + !isempty(residue) > 3 && continue
            # C8: product-only iso form
            step!(Mechanisms.Met[], false, consumed, released, Symbol[], subset)
        end
    end
    # Start from free enzyme (:E, no bound metabolites).
    walk!(Mechanisms.RawForm(:E, Mechanisms.Met[]), Symbol[], Symbol[], Symbol[], Symbol[],
          Mechanisms.RawStep[])

    # --- Group paths by isomerization pattern ---
    iso_groups = Dict{Set{Mechanisms.RawStep}, Vector{Vector{Mechanisms.RawStep}}}()
    for path in paths
        push!(get!(iso_groups, Set(filter(Mechanisms._is_iso, path)),
                   Vector{Mechanisms.RawStep}[]), path)
    end

    # --- Build topologies: union whole paths consistent with each
    # (substrate weak-ordering, product weak-ordering). Unioning complete
    # paths (rather than cherry-picking steps) keeps every form connected to
    # the catalytic complex — paths consistent with one weak ordering never
    # carry contradictory binding orders, so no dangling single-metabolite
    # forms arise. Binding history is read from the path, so ping-pong (where
    # a consumed substrate leaves the bound set) is handled correctly. A path
    # is consistent with a weak ordering when its binding (release) order never
    # steps back to an earlier level; every path binds every substrate and
    # releases every product once.
    #
    # Iterate iso_groups deterministically (smaller iso-step counts first, then by
    # their iso steps' forms in the canonical form order) so topology output order is
    # stable — a `Set`'s iteration order is not.
    level_of(wo) = Dict{Symbol,Int}(m => i for (i, level) in enumerate(wo) for m in level)
    sub_levels = level_of.(_weak_orderings(sub_names))
    prod_levels = level_of.(_weak_orderings(prod_names))
    # The names a path binds, in order, among `names`: a release counts as the binding it
    # is written as.
    route_order(path, names) =
        Symbol[s.consumed[1].name for s in path
               if !isempty(s.consumed) && s.consumed[1].name in names]
    # A form's place in the canonical form order (`_sorted_forms`): free E first, then by
    # label, bound metabolites, residual additions and residual subtractions.
    key(f::Mechanisms.RawForm) =
        (!(isempty(f.bound) && isempty(f.added) && isempty(f.subtracted) &&
           f.conformation === :E),
         f.conformation, f.bound, f.added, f.subtracted)
    form_pair(s) = (key(s.from), key(s.to))
    sorted_iso_pats = sort!(collect(keys(iso_groups));
                            by = pat -> (length(pat), sort([form_pair(s) for s in pat])))
    result = Vector{Mechanisms.RawStep}[]
    for iso_pat in sorted_iso_pats
        group = [(path, route_order(path, sub_names), route_order(path, prod_names))
                 for path in iso_groups[iso_pat]]
        seen_topos = Set{Set{Mechanisms.RawStep}}()
        for sl in sub_levels, pl in prod_levels
            topo_keys = Set{Mechanisms.RawStep}()
            for (path, binds, releases) in group
                issorted(binds; by = m -> sl[m]) && issorted(releases; by = m -> pl[m]) &&
                    union!(topo_keys, path)
            end
            (isempty(topo_keys) || topo_keys ∈ seen_topos) && continue
            push!(seen_topos, topo_keys)

            steps = sort!(collect(topo_keys);
                          by = s -> (Mechanisms._is_iso(s), form_pair(s)...))

            # The first iso step is the (single) SS step; every other step is RE.
            iso_idx = findfirst(Mechanisms._is_iso, steps)
            push!(result, Mechanisms.RawStep[
                Mechanisms.RawStep(s.from, s.to, s.consumed, s.released, i != iso_idx)
                for (i, s) in enumerate(steps)])
        end
    end
    isempty(result) && error("no catalytic cycle for $(join(sub_names, " + ")) ⇌ " *
                             "$(join(prod_names, " + ")): every route needs an " *
                             "isomerization converting more than three substrates, or " *
                             "more than three products counting a covalent residue as one")
    result
end

# ─── Dead-End Helpers ────────────────────────────────────────

"""
    _competition_patterns(sub_names, prod_names)
        → Vector{Set{Tuple{Symbol,Symbol}}}

Enumerate all bipartite competition graphs on
substrates × products where every substrate and every
product has degree ≥ 1.
"""
function _competition_patterns(
    sub_names::Set{Symbol},
    prod_names::Set{Symbol},
)
    edges = [(s, p) for s in sort(collect(sub_names)) for p in sort(collect(prod_names))]
    [Set(pat) for pat in _subsets(edges)[2:end]
     if all(s -> any(e -> e[1] == s, pat), sub_names) &&
        all(p -> any(e -> e[2] == p, pat), prod_names)]
end

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
    _expand_substrate_product_dead_ends(topos, reaction) -> Vector{Vector{RawStep}}

For each catalytic topology, enumerate substrate/product dead-end
form combinations and return each resulting mechanism as a flat
`Vector{RawStep}`. A dead-end form is a catalytic form with one more
substrate or product bound, subject to:
- The resulting form binds at least one substrate AND
  at least one product (mixed binding required)
- The resulting form doesn't have all substrates or
  all products
- A declared `shared_catalytic_site` pair is never both bound in the
  resulting form

Each competition pattern (`_competition_patterns`) keeps the dead-end forms that bind
none of its (substrate, product) pairs, and each distinct set of kept forms gives one
mechanism: the topology plus a rapid-equilibrium binding step between every two of its
forms with the same conformation and residual that differ by one bound metabolite.

Two preconditions hold for the topologies of `_catalytic_topologies`: every binding is
rapid equilibrium, and every isomerization joins a substrate-only form to a
product-only form. The first gives a binding between two dead-end forms the flag of the
topology binding it parallels. The second keeps every dead-end form off the catalytic
cycle and leaves no isomerization a dead-end parallel, so bindings are the only steps
the dead-end forms need.
"""
function _expand_substrate_product_dead_ends(
    topos::Vector{Vector{Mechanisms.RawStep}},
    reaction::Mechanisms.EnzymeReaction,
)
    sub_names = Set(Mechanisms.name(s) for s in Mechanisms.substrates(reaction))
    prod_names = Set(Mechanisms.name(p) for p in Mechanisms.products(reaction))
    mets = sort!(collect(union(sub_names, prod_names)))

    # Competition patterns depend only on the reaction, so they are computed once.
    # A declared shared catalytic site forbids its (substrate, product) pair
    # from co-occupying the catalytic site, so keep only competition patterns
    # whose forbidden-edge set contains every declared pair. The complete
    # bipartite pattern contains all edges and always survives, so the list is
    # never empty.
    patterns = filter(pat -> all(in(pat), Mechanisms.shared_catalytic_site(reaction)),
                      _competition_patterns(sub_names, prod_names))

    held(f::Mechanisms.RawForm) = Set(x.name for x in f.bound)
    result = Vector{Mechanisms.RawStep}[]
    for topo in topos
        catalytic = unique!([f for s in topo for f in (s.from, s.to)])
        dead_ends = Mechanisms.RawForm[]
        for f in catalytic, m in mets
            b = held(f)
            m in b && continue
            push!(b, m)
            subs, prods = intersect(b, sub_names), intersect(b, prod_names)
            (isempty(subs) || isempty(prods) || subs == sub_names || prods == prod_names) &&
                continue
            push!(dead_ends, Mechanisms.RawForm(f.conformation,
                                                [f.bound; Mechanisms.Met(m, false)],
                                                f.added, f.subtracted))
        end
        # Several catalytic forms can give the same dead end.
        unique!(dead_ends)
        edges = Set((s.from, s.to) for s in topo)
        seen = Set{Vector{Mechanisms.RawForm}}()
        for pattern in patterns
            active = filter(dead_ends) do d
                b = held(d)
                !any((s, p) in pattern
                     for s in intersect(b, sub_names) for p in intersect(b, prod_names))
            end
            active in seen && continue
            push!(seen, active)

            # Fully connect the enzyme-form graph: any two forms that are identical
            # except for one bound metabolite are joined by a binding step. This adds
            # each dead-end form's bindings to the forms it extends, the parallels of
            # the topology's bindings, and the bindings between two dead-end forms.
            # Added edges are RE bindings on the differing metabolite: `_seed_groups`
            # folds them into that metabolite's kinetic group, so an added edge adds no
            # free parameter whether it reaches a dead-end form or closes a
            # thermodynamically dependent cycle.
            forms = [catalytic; active]
            hs = held.(forms)
            steps = copy(topo)
            for (i, f1) in enumerate(forms), (j, f2) in enumerate(forms)
                f1.conformation == f2.conformation &&
                    (f1.added, f1.subtracted) == (f2.added, f2.subtracted) &&
                    length(hs[j]) == length(hs[i]) + 1 && issubset(hs[i], hs[j]) &&
                    !((f1, f2) in edges) || continue
                met = only(setdiff(hs[j], hs[i]))
                push!(steps, Mechanisms.RawStep(f1, f2, [Mechanisms.Met(met, false)],
                                                Mechanisms.Met[], true))
            end
            push!(result, steps)
        end
    end
    result
end

# ─── Mechanism Enumeration ───────────────────────────────────

"""
    _seed_groups(steps) -> Vector{Int}

The kinetic group label of each raw step of a seed, labels numbered in order of first
use: steps that take up and give off the same metabolites with the same RE/SS flag
share one group, since the steps of a group share their constants; every isomerization
is a group of its own.
"""
function _seed_groups(steps::Vector{Mechanisms.RawStep})
    label = Dict{Tuple{Vector{Mechanisms.Met}, Vector{Mechanisms.Met}, Bool, Int}, Int}()
    [get!(label,
          (s.consumed, s.released, s.rapid_equilibrium, Mechanisms._is_iso(s) ? j : 0),
          length(label) + 1)
     for (j, s) in enumerate(steps)]
end


# ─── Expansion Moves ─────────────────────────────────────────

"""
    _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})

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
endpoints share a segment, which the rate equation never sees. A flip set whose child the
guard's class 3 rejects on some face (`_faces_pass`) also counts as failed, so the
minimal-set search extends it instead of emitting it; so does
a set one of whose flipped groups carries no flux in the child
(`_flux_carrying_groups` on the flipped groups): the child is its parent plus
one phantom per such group, and its flux-carrying supersets are the smallest
identifiable relaxations. The non-degenerate mechanisms beyond a set that counts as failed
stay reachable — through the other cut orders of the same ring, or through the extended
set when no other order gains. Admission's face condition does not count: a set whose
allosteric child admission rejects (`_BRIDGE`) is emitted, the exit or identify's base
tier drops that child, and its supersets are not tried. Every other flag, the declared
names and, for an allosteric mechanism, the tags, multiplicity and regulatory sites are
kept; each child goes through the checked constructor.

A parent that `_requires_hyperbolic_catalysis` (a conformational mechanism over
more than one catalytic subunit) keeps only the children whose laws pass
`_hyperbolic_catalysis` (each conformation's, for an allosteric child). The check
runs on the emitted minimal sets, not inside the gain predicate: a failing set
would otherwise be extended with more flips, and more steady-state steps never
restore a hyperbolic equation, so its supersets need no visit.
"""
function _expand_re_to_ss(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism})
    cm = Mechanisms.catalytic(m)
    flanks = _chain_flank_groups(m)
    modifiers = Mechanisms._modifier_mask(cm)
    exchanges_modifier = Mechanisms._groups_with(
        k -> (cm.edges[k].consumed | cm.edges[k].released) & modifiers != 0, cm.group,
        Mechanisms.ngroups(cm))
    eligible(g) = !(g in flanks) && !exchanges_modifier[g]
    base = _re_segment_count(cm, Mechanisms._flags(cm))
    flips = _minimal_flips(re -> _re_segment_count(cm, re) > base, m, eligible)
    children = typeof(m)[Mechanisms._build(Mechanisms._raw(m, re)) for re in flips]
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
function _chain_flank_groups(m::Mechanisms.Mechanism)
    # The group of each edge alone in its group, 0 for any other edge.
    lone = zeros(Int, length(m.edges))
    for (g, ks) in enumerate(Mechanisms._group_members(m))
        length(ks) == 1 && (lone[only(ks)] = g)
    end
    at = [Int[] for _ in m.forms]
    for (k, e) in enumerate(m.edges)
        push!(at[e.from], k); push!(at[e.to], k)
    end
    # The group of the flank of edge `k0` at form `x`, 0 when there is none.
    function flank(k0, x)
        rest = filter(!=(k0), at[x])
        length(rest) == 1 || return 0
        e = m.edges[only(rest)]
        Mechanisms._is_binding(e) && e.to == x ? lone[only(rest)] : 0
    end
    out = Set{Int}()
    for (k0, e0) in enumerate(m.edges)
        lone[k0] != 0 && Mechanisms._is_iso(e0) && !e0.rapid_equilibrium || continue
        f1, f2 = flank(k0, e0.from), flank(k0, e0.to)
        f1 == 0 || f2 == 0 || union!(out, (f1, f2))
    end
    out
end
_chain_flank_groups(::Mechanisms.AllostericMechanism) = Set{Int}()

"""
    _minimal_flips(admissible, m, eligible = _ -> true) -> Vector{BitVector}

The rapid-equilibrium flags of the catalytic graph of `m`, plain or allosteric, after
every inclusion-minimal flip to steady state that passes, one flag per edge for each flip
set. A unit is a rapid-equilibrium kinetic group that carries flux with every edge steady
state (`_flux_carrying_groups`) and passes `eligible`. A candidate turns its units' edges
steady state and keeps every other flag. It passes when `admissible` holds on its flags,
every flipped group carries flux in it, and the guard's classes 2 and 3 pass on every face
(`_faces_pass`, which reads `m`'s law network under the candidate's flags). The sets come
from `_minimal_gaining_sets`, so a failing set is extended and a passing one is not.
Admission's face condition is not part of the test: a candidate whose allosteric child
admission rejects (`_BRIDGE`) passes, and its supersets are not tried. Nothing is built: a
candidate is its flags.
"""
function _minimal_flips(admissible, m::Mechanisms._AnyMech, eligible = _ -> true)
    cm = Mechanisms.catalytic(m)
    flags = Mechanisms._flags(cm)
    members = Mechanisms._group_members(cm)
    flux = _flux_carrying_groups(cm, falses(length(cm.edges)))
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
                   _faces_pass(m, re))
    [flipped(sel) for sel in _minimal_gaining_sets(passes, _ -> 1:length(units))]
end

"""
    _flux_carrying_edges(m, re) -> BitVector

Whether each edge of graph `m`, edge `k` at rapid equilibrium where `re[k]` holds,
carries net steady-state flux for generic parameter values. Rapid-equilibrium edges hold
equilibrium mass and are flagged `false`; the rule that uses these flags concerns
steady-state groups only.

Contract each rapid-equilibrium segment to one vertex (`_segments`) and make every
steady-state edge an edge between its two forms' segments, weighted by the substrates it
takes up minus the products, plus its `from` form's offsets minus its `to` form's (a
modifier copy counts as its reactant). Around any cycle the offsets telescope and the
uptakes sum to the net turnover times the reactant count, so a cycle has weight zero
exactly when it runs no net reaction. A self-loop (both forms in one segment) carries
flux iff its weight is nonzero. Any other edge carries flux iff its biconnected block
(`_unbalanced_blocks`) holds a cycle of nonzero weight: every edge of such a block lies
on one (join the edge to the cycle by two disjoint paths; one of the two resulting
cycles is unbalanced), while in a balanced block detailed balance holds along every
cycle and each edge's flux vanishes. The zero-flux verdict holds for any parameters and
grouping; the flux-carrying verdict for one-way steps, which every step the enumerator
emits is. The test reads each edge's metabolite lists, so it holds for fused and
Theorell–Chance steps as well.
"""
function _flux_carrying_edges(m::Mechanisms.Mechanism, re::BitVector)
    s = Mechanisms._segments(m, trues(length(m.edges)), re)
    name_of = Mechanisms._name_index(m)
    # The net uptake sign of each table entry and of each metabolite name.
    sign = _entry_signs(m)
    rho = zeros(Int, size(s.offsets, 1))
    for (i, n) in enumerate(name_of)
        rho[n] = sign[i]
    end
    offset(f) = sum(rho[n] * s.offsets[n, f] for n in eachindex(rho); init = 0)
    weight(e) = _uptake_weight(e, sign) + offset(e.from) - offset(e.to)
    flags = falses(length(m.edges))
    edges = Tuple{Int, Int}[]; weights = Int[]; owner = Int[]
    for (k, e) in enumerate(m.edges)
        re[k] && continue
        a, b = s.segment_of[e.from], s.segment_of[e.to]
        if a == b
            flags[k] = weight(e) != 0
        else
            push!(edges, (a, b)); push!(weights, weight(e)); push!(owner, k)
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
function _flux_carrying_groups(m::Mechanisms.Mechanism, re::BitVector)
    flux = _flux_carrying_edges(m, re)
    Mechanisms._groups_with(k -> flux[k], m.group, Mechanisms.ngroups(m))
end

"""Each table entry's net uptake sign by its name: +1 for a substrate, −1 for a product,
0 for a modifier-only name. A reactant's modifier copy counts as its reactant, whose
concentration it shares."""
_entry_signs(m::Mechanisms.Mechanism) =
    Int[insorted(x.name, m.substrates) ? 1 : insorted(x.name, m.products) ? -1 : 0
        for x in m.table]

"""The signs `sign` (`_entry_signs`) summed over the entries edge `e` takes up, minus
those summed over the entries it gives off."""
_uptake_weight(e::Mechanisms.Edge, sign::Vector{Int}) =
    sum(sign[i] for i in Mechanisms.entries(e.consumed); init = 0) -
    sum(sign[i] for i in Mechanisms.entries(e.released); init = 0)

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
    function visit(u, parent_edge)
        clock[] += 1; disc[u] = low[u] = clock[]
        for e in adj[u]
            e == parent_edge && continue
            w = edges[e][1] == u ? edges[e][2] : edges[e][1]
            if disc[w] == 0
                phi[w] = phi[u] + (edges[e][1] == u ? weights[e] : -weights[e])
                push!(stack, e)
                visit(w, e)
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
        disc[v] == 0 && visit(v, 0)
    end
    block, Set{Int}(block[e] for (e, (u, v)) in enumerate(edges)
                    if weights[e] != phi[v] - phi[u])
end

"""
`raw` with the isomerization `s0`, alone in its kinetic group, merged onto its product
side: `s0` and its group are removed, and every other step at `s0.from` moves onto
`s0.to`. The last substrate's binding into the substrate side becomes a fused binding
into the product side; the releases stay plain. The other groups keep their order.
"""
function _merge_isomerization(raw::Mechanisms.RawMechanism, s0::Mechanisms.RawStep)
    move(f::Mechanisms.RawForm) = f == s0.from ? s0.to : f
    keep = [s != s0 for s in raw.steps]
    steps = Mechanisms.RawStep[Mechanisms.RawStep(move(s.from), move(s.to), s.consumed,
                                                  s.released, s.rapid_equilibrium)
                               for s in raw.steps[keep]]
    Mechanisms.RawMechanism(raw.substrates, raw.products, raw.modifiers, steps,
                            first(Mechanisms._number_groups(raw.group[keep],
                                                            maximum(raw.group))))
end

"""
`raw` with the form `x` eliminated, or `nothing` when that is not possible: `x` must
have exactly two steps, both bindings into it, Y + L → x and Z + R → x with Y ≠ Z. They
become one Theorell–Chance step Y + L → Z + R, rapid equilibrium, in a group of its own
after the others, where Y + L is the binding of a substrate when one of the two is; their
groups keep their other steps.
"""
function _eliminate_form(raw::Mechanisms.RawMechanism, x::Mechanisms.RawForm)
    at_x = findall(s -> x in (s.from, s.to), raw.steps)
    length(at_x) == 2 && all(j -> raw.steps[j].to == x && isempty(raw.steps[j].released) &&
                                  length(raw.steps[j].consumed) == 1, at_x) ||
        return nothing
    binds_substrate(j) = (y = only(raw.steps[j].consumed);
                          !y.modifier && y.name in raw.substrates)
    entry, release = raw.steps[sort(at_x; by = j -> !binds_substrate(j))]
    entry.from == release.from && return nothing
    keep = [!(j in at_x) for j in eachindex(raw.steps)]
    steps = [raw.steps[keep]; Mechanisms.RawStep(entry.from, release.from, entry.consumed,
                                                 release.consumed, true)]
    group = [raw.group[keep]; maximum(raw.group) + 1]
    Mechanisms.RawMechanism(raw.substrates, raw.products, raw.modifiers, steps,
                            first(Mechanisms._number_groups(group, maximum(group))))
end

"""Whether a cycle of the edges of graph `m` that `fast` flags runs net turnover, which
makes the rate infinite when those edges are at rapid equilibrium. Each edge joins its two
forms with its uptake of the reaction's substrates minus its products as weight
(`_uptake_weight`); a turnover cycle has weight n·Σρ², never zero, so the flagged edges
hold one iff they close a cycle of nonzero weight (`_has_unbalanced_cycle`)."""
function _re_turnover_cycle(m::Mechanisms.Mechanism, fast::BitVector)
    sign = _entry_signs(m)
    n = length(m.forms)
    Mechanisms._has_unbalanced_cycle([e.from for e in m.edges], [e.to for e in m.edges],
                                     [_uptake_weight(e, sign) for e in m.edges], fast,
                                     zeros(Int, n), zeros(Int, n))
end

"""Whether the rate of graph `m`, edge `k` at rapid equilibrium where `re[k]` holds, has a
maximum as the metabolites `names` (its substrates or its products, sorted) grow: turning
rapid equilibrium every steady-state edge that takes up or gives off one of them leaves no
rapid-equilibrium turnover cycle (condition V). Only the edges `on` flags count, so an
allosteric mechanism's test can read its active law (`_degenerate`)."""
function _has_vmax(m::Mechanisms.Mechanism, re::BitVector, names::Vector{Symbol},
                   on::BitVector = trues(length(m.edges)))
    side = Mechanisms._names_mask(m, names)
    touches = BitVector([(e.consumed | e.released) & side != 0 for e in m.edges])
    !_re_turnover_cycle(m, (re .| touches) .& on)
end

"""The merged complexes of graph `m`, each once: the forms that a fused binding of a
substrate enters, a chemistry edge (`_is_chemistry`) that binds an entry of the mask
`subs`. Only the edges `on` flags count."""
_merged_complexes(m::Mechanisms.Mechanism, subs::Mechanisms.Counts,
                  on::BitVector = trues(length(m.edges))) =
    unique(e.to for (k, e) in enumerate(m.edges)
           if on[k] && Mechanisms._is_binding(e) && e.consumed & subs != 0 &&
              Mechanisms._is_chemistry(m, e))

"""Whether edge `e` takes up an entry of one of the masks `subs` and `prods` (a reactant
of one side of the reaction) and gives off one of the other, as a Theorell–Chance step
does."""
_crosses_sides(e::Mechanisms.Edge, subs::Mechanisms.Counts, prods::Mechanisms.Counts) =
    e.consumed & subs != 0 && e.released & prods != 0 ||
    e.consumed & prods != 0 && e.released & subs != 0

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
function _chemistry_equilibrates_both_sides(m::Mechanisms.Mechanism, re::BitVector,
                                            on::BitVector = trues(length(m.edges)))
    re = re .& on
    subs, prods = Mechanisms._names_mask(m, m.substrates),
                  Mechanisms._names_mask(m, m.products)
    releases(node, side) = any(eachindex(m.edges)) do k
        e = m.edges[k]
        re[k] && (e.to in node && !(e.from in node) && e.consumed & side != 0 ||
                  e.from in node && !(e.to in node) && e.released & side != 0)
    end
    # The forms joined by RE isomerizations, one node per segment; a merged complex
    # outside every RE isomerization is a node of its own.
    isos = BitVector([re[k] && Mechanisms._is_iso(e) for (k, e) in enumerate(m.edges)])
    s = Mechanisms._segments(m, isos, isos)
    nodes = [s.segments; [[x] for x in _merged_complexes(m, subs, on)
                          if s.segment_of[x] == 0]]
    any(node -> releases(node, subs) && releases(node, prods), nodes) ||
        any(k -> re[k] && _crosses_sides(m.edges[k], subs, prods), eachindex(m.edges))
end

"""Whether `m` is degenerate by the structural conditions V and C: no maximal rate in one
direction (`_has_vmax`), or chemistry in rapid equilibrium with both sides
(`_chemistry_equilibrates_both_sides`). Read on the catalytic graph with its own flags,
over the active conformation's law: every edge of a plain mechanism, and the edges
`_active_mask` keeps of an allosteric one. The beam expands a degenerate seed without
fitting it (`_base_tier`)."""
function _degenerate(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism})
    cm = Mechanisms.catalytic(m)
    re = Mechanisms._flags(cm)
    on = Mechanisms._active_law(m)
    !(_has_vmax(cm, re, cm.substrates, on) && _has_vmax(cm, re, cm.products, on)) ||
        _chemistry_equilibrates_both_sides(cm, re, on)
end

"""
The tests of `_seed_variants` that build nothing, for the graph `m` of one base of a
seed: a function of a candidate's rapid-equilibrium flags (one per edge of `m`, some
turned steady state) that is true when the candidate has no rapid-equilibrium turnover
cycle (`_re_turnover_cycle`), a maximal rate both ways (`_has_vmax`) and no chemistry in
equilibrium with both sides (`_chemistry_equilibrates_both_sides`). Each edge's ends,
uptake weight and reactants are read once per base, so a candidate costs a few passes
over arrays. A base holds no isomerization, so each chemistry node is a merged complex
or a Theorell–Chance step. A set of edges holds a turnover cycle iff it closes a cycle of
nonzero weight (`_has_unbalanced_cycle`), which reuses the screen's buffers. Each
maximal-rate test runs on the candidate's rapid-equilibrium edges and more, and adding
edges never removes such a cycle, so a turnover cycle of the candidate fails both tests.
"""
function _seed_candidate_screen(m::Mechanisms.Mechanism)
    any(Mechanisms._is_iso, m.edges) &&
        error("_seed_candidate_screen: the seed base holds an isomerization; a base " *
              "has every isomerization merged")
    sign = _entry_signs(m)
    subs, prods = Mechanisms._names_mask(m, m.substrates),
                  Mechanisms._names_mask(m, m.products)
    from = [e.from for e in m.edges]
    to = [e.to for e in m.edges]
    weight = [_uptake_weight(e, sign) for e in m.edges]
    touches(side) = BitVector([(e.consumed | e.released) & side != 0 for e in m.edges])
    sub_step, prod_step = touches(subs), touches(prods)
    crossing = BitVector([_crosses_sides(e, subs, prods) for e in m.edges])
    complexes = _merged_complexes(m, subs)
    leaves(x, side) = BitVector([e.to == x && e.consumed & side != 0 ||
                                 e.from == x && e.released & side != 0 for e in m.edges])
    exits = [(leaves(x, subs), leaves(x, prods)) for x in complexes]
    parent = zeros(Int, length(m.forms)); offset = zeros(Int, length(m.forms))
    turnover(on) = Mechanisms._has_unbalanced_cycle(from, to, weight, on, parent, offset)
    meets(a::BitVector, b::BitVector) = any(k -> a[k] && b[k], eachindex(a))
    kept = falses(length(m.edges))
    function screen(re::BitVector)
        kept .= re .| sub_step
        turnover(kept) && return false
        kept .= re .| prod_step
        turnover(kept) && return false
        meets(re, crossing) && return false
        !any(((es, ep),) -> meets(re, es) && meets(re, ep), exits)
    end
end

"""
The bases of the seed `m`'s variants, as raw input (`_raw`), the merged base first:
merging every isomerization of `m` onto its product side (`_merge_isomerization`) and
setting every step at rapid equilibrium gives the merged base; eliminating one merged
complex with two steps (`_eliminate_form`) gives a Theorell–Chance base, one per complex
that can be eliminated, in the order of the isomerizations.
"""
function _seed_bases(m::Mechanisms.Mechanism)
    raw = Mechanisms._raw(m)
    isos = filter(Mechanisms._is_iso, raw.steps)
    joined = foldl(_merge_isomerization, isos; init = raw)
    merged = Mechanisms.RawMechanism(joined.substrates, joined.products, joined.modifiers,
                                     Mechanisms.RawStep[
                                         Mechanisms.RawStep(s.from, s.to, s.consumed,
                                                            s.released, true)
                                         for s in joined.steps], joined.group)
    Mechanisms.RawMechanism[merged;
                            filter(!isnothing,
                                   [_eliminate_form(merged, s.to) for s in isos])]
end

"""Whether the flags `re` make the merged base `base` its seed's lumping twin: each merged
complex, forms `cx` of `base`, has two edges, both steady state and each alone in its
group."""
_lumping_twin(base::Mechanisms.Mechanism, cx::Vector{Int}, re::BitVector) = all(cx) do f
    at = findall(e -> f in (e.from, e.to), base.edges)
    length(at) == 2 && all(k -> !re[k] && count(==(base.group[k]), base.group) == 1, at)
end

"""
The merged and Theorell–Chance variants of the seed `m` of reaction `rxn`, whose atoms
each base's edges must conserve. Each base (`_seed_bases`) is read as the graph
`_canonicalize` gives it, unchecked; a base whose tied groups no orientation aligns raises
the representation error (`_check_aligned`). MERGE and ELIM conserve atoms, which each
base asserts. For each base, every inclusion-minimal set of its groups (`_minimal_flips`)
whose flip to steady state gives a valid, flux-carrying, non-degenerate candidate is a
variant: no rapid-equilibrium turnover cycle (`_re_turnover_cycle`), a maximal rate both
ways (`_has_vmax`), no chemistry in equilibrium with both sides
(`_chemistry_equilibrates_both_sides`), every steady-state group carrying flux, and the
guard's classes 2 and 3 on every face (`_faces_pass`). The first three tests read index
arrays (`_seed_candidate_screen`); the other two run only on a candidate that passes them.
A merged variant that is its seed's lumping twin (`_lumping_twin`) is skipped: the unmerged
form with rapid-equilibrium flanks has its family at the same count. A base's groups come in
canonical order, so the variants come in the order of their flip sets over those groups;
each is built and checked.
"""
function _seed_variants(m::Mechanisms.Mechanism, rxn::Mechanisms.EnzymeReaction)
    variants = Mechanisms.Mechanism[]
    for (i, b) in enumerate(_seed_bases(m))
        c = Mechanisms._canonicalize(b)
        Mechanisms._check_aligned(c.mechanism, c.misaligned)
        base = c.mechanism
        _assert_atom_conserving(base, rxn)
        flips = _minimal_flips(_seed_candidate_screen(base), base)
        if i == 1
            cx = _merged_complexes(base, Mechanisms._names_mask(base, base.substrates))
            filter!(re -> !_lumping_twin(base, cx, re), flips)
        end
        append!(variants, [Mechanisms.Mechanism(Mechanisms._raw(base, re)) for re in flips])
    end
    variants
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
    _hyperbolic_catalysis(m::Mechanism, on = trues(length(m.edges))) -> Bool

Whether the King–Altman denominator of `m`'s catalytic scheme has degree at most
1 in every substrate and product concentration. Every binding of a substrate or
product at its catalytic site counts, abortive complexes included. Edges touching
a form that holds a modifier are left out (a mask of the catalytic-site edges): a
modifier binds a site of its own by definition, so the powers its binding adds,
including those of a reactant's modifier copy, are a separate source that
conformational mechanisms keep. Only the edges `on` flags count, and `_segments` leaves
out the forms they do not keep, so an allosteric mechanism's test can read each
conformation's law.

The equation is a sum over rapid-equilibrium (RE) segments and spanning
arborescences of the segment graph toward each segment. A denominator term is
the root segment's weight times the weight of every tree edge, so the exponent
of `X` in a term is the most `X` any form of the root segment carries beyond
the segment's bottom form, plus, per tree edge, one if the step binds `X` in
the tree direction and the count of `X` the edge's source form carries beyond
its bottom (`_segments`). The degree exceeds 1 exactly when one
segment or one edge scores 2 or more, or the root scores 1 and some
arborescence toward it holds a scoring edge, or some arborescence holds two
scoring edges. An arborescence toward `S` containing given edges exists iff
every segment still reaches `S` once each given edge's source keeps that edge
as its only way out (`_all_reach`).
"""
function _hyperbolic_catalysis(m::Mechanisms.Mechanism,
                               on::BitVector = trues(length(m.edges)))
    modifiers = Mechanisms._modifier_mask(m)
    on_site = on .& BitVector([(m.forms[e.from].bound | m.forms[e.to].bound) &
                               modifiers == 0 for e in m.edges])
    s = Mechanisms._segments(m, on_site, Mechanisms._flags(m))
    name_of = Mechanisms._name_index(m)
    # Directed segment-graph edges: source segment, target segment, source form, the
    # names (`_name_index`) bound in that direction.
    edges = Tuple{Int, Int, Int, Vector{Int}}[]
    for (k, e) in enumerate(m.edges)
        on_site[k] && !e.rapid_equilibrium || continue
        a, b = s.segment_of[e.from], s.segment_of[e.to]
        a == b && continue
        push!(edges, (a, b, e.from, [name_of[i] for i in Mechanisms.entries(e.consumed)]))
        push!(edges, (b, a, e.to, [name_of[i] for i in Mechanisms.entries(e.released)]))
    end
    n = length(s.segments)
    for y in [m.substrates; m.products]
        x = name_of[Mechanisms._entry(m, Mechanisms.Met(y, false))]
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

Whether each conformation's law of `am` passes `_hyperbolic_catalysis`: the active one
(`_active_mask`) and the inactive one (`_inactive_mask`). A conformation that runs part
of the catalytic graph can carry a power the whole graph lacks: a ping-pong's abortive
complexes can put A² into the inactive law once that conformation stops one
half-reaction.
"""
_hyperbolic_catalysis(am::Mechanisms.AllostericMechanism) =
    _hyperbolic_catalysis(am.catalytic, Mechanisms._active_mask(am)) &&
    _hyperbolic_catalysis(am.catalytic, Mechanisms._inactive_mask(am))

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

"""Number of rapid-equilibrium segments of graph `m`, edge `k` at rapid equilibrium where
`re[k]` holds (`_segments`). The flip move measures each candidate on its flags, so it
builds only the candidates that pass."""
_re_segment_count(m::Mechanisms.Mechanism, re::BitVector) =
    length(Mechanisms._segments(m, trues(length(m.edges)), re).segments)

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
    _expand_split_kinetic_group(m::Mechanism) → Vector{Mechanism}
    _expand_split_kinetic_group(am::AllostericMechanism) → Vector{AllostericMechanism}

Kinetic-group split move. A unit is one context bipartition of one group
(`_context_bipartitions`: the affinity for a metabolite may depend on which
other ligand is already bound). One child is produced per minimal set of units,
at most one per group, whose joint application raises the independent-parameter
count (`_minimal_gaining_sets`). A single bipartition whose new constant a
Wegscheider cycle ties straight back is not a model: the cycle ties the new constant to the
others, so the child's rate function is the parent's. Such a bipartition
is kept as a seed for pairs and larger sets, because the cycle that absorbs it
is broken by also splitting another group on that cycle — random-order binding
needs one split per substrate before any constant frees up. Partners for a
rejected set are the all-RE groups sharing a rapid-equilibrium segment with any
step of the split groups: a tie needs an RE cycle through the carve, and every
group on that cycle lies in that segment. The count is the joint solve's on the built
child (`_split_gain_test`).
A part of a steady-state group none of whose edges carries flux
(`_flux_carrying_edges`, computed once on the parent, since a split moves no edge)
is emitted at rapid equilibrium (`_revert_zero_flux_parts`), the same family with
one constant fewer; the gain test counts the reverted part under
its new kind, so a reverted constant the Wegscheider ties pull back is absorbed
like any tied split, and a candidate whose reverted groups leave a rapid-equilibrium
segment without a bottom form counts as failed. A split of any group can complete a
copy's gauge, by separating a binding that forms or leaves a twin from the bindings
elsewhere in its group that blocked the gauge (the split then holds the family of the
copy its shared group pinned), and a split of a group that forms a copy's
twins can break the gauge, so the copy rule is decided on each child: a child with a
redundant copy group (`_redundant_copy_groups`) is not emitted; its family is that of
the same split without the copy.
The declared names and, for an allosteric mechanism, the multiplicity and regulatory
sites are kept; both parts of a split group inherit its catalytic allo-state tag. Each
child is the parent's raw input regrouped (`_apply_bipartitions`), through the checked
constructor, since a regrouping can turn a tied group's member.
"""
function _expand_split_kinetic_group(m::Union{Mechanisms.Mechanism,
                                               Mechanisms.AllostericMechanism})
    cm = Mechanisms.catalytic(m)
    members = Mechanisms._group_members(cm)
    flags = Mechanisms._flags(cm)
    flux = _flux_carrying_edges(cm, flags)
    units = Tuple{Int, Tuple{Vector{Int}, Vector{Int}}}[]
    reverted = Vector{Int}[]
    for g in eachindex(members), bp in _context_bipartitions(cm, g)
        push!(units, (g, bp)); push!(reverted, _revert_zero_flux_parts(cm, bp, flux))
    end
    isempty(units) && return typeof(m)[]
    # The flags of a candidate: the parent's, with each selected unit's reverted edges at
    # rapid equilibrium.
    function flags_of(sel)
        re = copy(flags)
        for u in sel
            re[reverted[u]] .= true
        end
        re
    end
    gain = _split_gain_test(m, units, flags_of)
    every = trues(length(cm.edges))
    gains(sel) = (all(u -> isempty(reverted[u]), sel) ||
                  isempty(_bottomless_segment(cm, every, flags_of(sel)))) && gain(sel)
    # RE segment ids touched by each kinetic group's edges; empty for a steady-state
    # group, which lies on no RE cycle and so never partners a set that does not gain.
    segment_of = Mechanisms._segments(cm, every, flags).segment_of
    segments = [flags[first(ks)] ? Set{Int}(segment_of[cm.edges[k].from] for k in ks) :
                Set{Int}() for ks in members]
    partners(sel) = begin
        isempty(sel) && return 1:length(units)
        used = Set(units[u][1] for u in sel)
        touched = reduce(union, (segments[units[u][1]] for u in sel))
        [u for u in 1:length(units)
         if !(units[u][1] in used) && !isempty(intersect(segments[units[u][1]], touched))]
    end
    sets = _minimal_gaining_sets(gains, partners)
    filter!(c -> isempty(_redundant_copy_groups(c)),
            typeof(m)[_apply_bipartitions(m, units[sel], flags_of(sel)) for sel in sets])
end

"""
The split move's gain test: `gains(sel)` is true when the child built from the selected
units, with the flags `flags_of(sel)` (`_apply_bipartitions`), has more independent
parameters than the parent by the joint solve (`_independent_param_count`), which derives
no rate polynomial. For an allosteric mechanism, both parts of a split group inherit its
tag and admission's interior reads each edge's tag, so a split never changes the interior
verdict (`_admission_interior`); the faces can change through the new group ties and the
parts it reverts to rapid equilibrium, and the exit (`expand_mechanisms`) drops such a
child.
"""
function _split_gain_test(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism},
                          units, flags_of)
    base = _independent_param_count(m)
    gains(sel) =
        _independent_param_count(_apply_bipartitions(m, units[sel], flags_of(sel))) > base
    gains
end

"""
The end of edge `e` that does not carry the edge's free metabolites: the form they
bind to. `e.from` when the edge takes a metabolite up (every binding, plain or fused,
is stored with its metabolite taken up), `e.to` when it only gives metabolites off (two
or more: an edge giving off one is stored as the binding it reverses), and `e.from` for
an isomerization (both lists empty).
"""
_context_form(e::Mechanisms.Edge) = e.consumed == 0 && e.released != 0 ? e.to : e.from

"""
    _context_bipartitions(m, g) -> Vector{Tuple{Vector{Int}, Vector{Int}}}

Ways to divide kinetic group `g` of graph `m` in two by binding context, each part a
list of edges: another ligand already bound, the enzyme's conformation, or its covalent
residual. Encodes "the affinity for this metabolite may depend on what else is bound, on
the conformation, or on the covalent state" — the ligand family includes a modifier,
which separates a catalytic step from its inhibitor-bound mirror. Each context value
divides the edges whose context form (`_context_form`) carries it from the rest; the
group's own bound entry is never a context, and a group whose forms share one
conformation and one residual gets no bipartition from those two families. Both parts
are nonempty, the first part holds the group's first edge, each part lists its edges in
edge order, context values that induce the same division give one bipartition, and the
order is ligands (modifiers, then products, then substrates, each by name), then
conformations (by label), then residuals (by the names they add, then those they
subtract), for deterministic output.
"""
function _context_bipartitions(m::Mechanisms.Mechanism, g::Int)
    ks = findall(==(g), m.group)
    own = Mechanisms._bound_entry(m.edges[first(ks)])
    forms = [m.forms[_context_form(m.edges[k])] for k in ks]
    ligands = unique!([i for f in forms for i in Mechanisms.entries(f.bound) if i != own])
    role(i) = m.table[i].modifier ? 0 : insorted(m.table[i].name, m.products) ? 1 : 2
    division(carries) = Int[k for (k, f) in zip(ks, forms) if carries(f)]
    divisions = [division(f -> Mechanisms.copies(f.bound, i) > 0)
                 for i in sort!(ligands; by = i -> (role(i), m.table[i].name))]
    append!(divisions, division(f -> f.conformation == c)
            for c in sort!(unique(f.conformation for f in forms)))
    append!(divisions, division(f -> (f.added, f.subtracted) == r)
            for r in sort!(unique((f.added, f.subtracted) for f in forms);
                           by = r -> (Mechanisms._entry_names(m, r[1]),
                                      Mechanisms._entry_names(m, r[2]))))
    parts(with) = (without = Int[k for k in ks if !(k in with)];
                   first(ks) in with ? (with, without) : (without, with))
    unique!(first, Tuple{Vector{Int}, Vector{Int}}[
        parts(with) for with in divisions if 0 < length(with) < length(ks)])
end

"""
`m` regrouped by `selection`, pairs of a kinetic group of `m` and one of its
`_context_bipartitions`, each group at most once, as `(group, origin)`: edge `k` in
group `group[k]`, the groups numbered as `m`'s with each selected group's second part
in a group of its own right after it, and the group of `m` that each new group `l` came
from, `origin[l]`.
"""
function _bipartitioned_groups(m::Mechanisms.Mechanism, selection)
    second = Dict(g => last(bp) for (g, bp) in selection)
    group, origin = similar(m.group), Int[]
    for (g, ks) in enumerate(Mechanisms._group_members(m))
        push!(origin, g)
        group[ks] .= length(origin)
        haskey(second, g) || continue
        push!(origin, g)
        group[second[g]] .= length(origin)
    end
    group, origin
end

"""
The child of `m` with each selected group divided in two (`_bipartitioned_groups`) and
edge `k` at rapid equilibrium where `re[k]` holds, built through the checked constructor.
For an allosteric mechanism both parts inherit the group's catalytic allo-state tag
(splitting is a parameter-relaxation move that must not change A/I semantics).
"""
function _apply_bipartitions(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism},
                             selection,
                             re::BitVector = Mechanisms._flags(Mechanisms.catalytic(m)))
    group, origin = _bipartitioned_groups(Mechanisms.catalytic(m), selection)
    m isa Mechanisms.Mechanism ? Mechanisms.Mechanism(Mechanisms._raw(m, re, group)) :
        Mechanisms.AllostericMechanism(
            Mechanisms._raw(m, re, group, Mechanisms.cat_allo_states(m)[origin]))
end

"""
The steady-state edges of the bipartition `bp` of a group of `m` that the split emits at
rapid equilibrium: those of each part none of whose edges carries flux (`flux`, the
parent's `_flux_carrying_edges`). A rapid-equilibrium group's edges carry none and are
rapid equilibrium already, so it reverts none.
"""
_revert_zero_flux_parts(m::Mechanisms.Mechanism, bp, flux::BitVector) =
    Int[k for part in bp if !any(k -> flux[k], part)
        for k in part if !m.edges[k].rapid_equilibrium]

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
the site's rapid-equilibrium segment with the complex's offsets (`_segments`; the
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
    s = Mechanisms._segments(m, mask, Mechanisms._flags(m))
    name_of = Mechanisms._name_index(m)
    modifiers = Mechanisms._modifier_mask(m)
    # A form's bound names, sorted, its conformation and its residual.
    composition(f) = (Mechanisms._entry_names(m, m.forms[f].bound), m.forms[f].conformation,
                      m.forms[f].added, m.forms[f].subtracted)
    function (site::Int, ligand::Symbol)
        i = findfirst(x -> x.name === ligand, m.table)
        i === nothing && return 0
        n = name_of[i]
        twins = [f for f in s.segments[s.segment_of[site]]
                 if m.forms[f].bound & modifiers == 0 &&
                    all(x -> s.offsets[x, f] == s.offsets[x, site] + (x == n),
                        axes(s.offsets, 1))]
        isempty(twins) && return 0
        names, label, added, subtracted = composition(site)
        wanted = (sort!(push!(names, ligand)), label, added, subtracted)
        k = findfirst(t -> composition(t) == wanted, twins)
        k === nothing ? minimum(twins) : twins[k]
    end
end

"""
Whether the edges `ks` of one kinetic group of `m` bind a modifier, and only at sites
whose complex has a productive twin (`twin`, a `_productive_twin` of one state's graph).
"""
function _all_twin(m::Mechanisms.Mechanism, ks::Vector{Int}, twin)
    i = Mechanisms._bound_entry(m.edges[first(ks)])
    i != 0 && m.table[i].modifier &&
        all(k -> twin(m.edges[k].from, m.table[i].name) != 0, ks)
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
    ligand = m.table[Mechanisms._bound_entry(m.edges[first(members[g])])].name
    twin_of = Dict(m.edges[k].to => twin(m.edges[k].from, ligand) for k in members[g])
    binding_group = Dict{Tuple{Int, Int}, Int}()
    for (h, ks) in enumerate(members), k in ks
        e = m.edges[k]
        i = Mechanisms._bound_entry(e)
        e.rapid_equilibrium && i != 0 && !m.table[i].modifier &&
            m.table[i].name === ligand && (binding_group[(e.from, e.to)] = h)
    end
    class = Dict{Int, _GaugeClass}(complex => label(g, 0) for complex in keys(twin_of))
    for k in members[g]
        t = twin_of[m.edges[k].to]
        haskey(class, t) && continue
        h = count(==(t), values(twin_of)) == 1 ?
            get(binding_group, (m.edges[k].from, t), 0) : 0
        class[t] = h > 0 ? label(h, 0) : label(0, t)
    end
    σ(f) = get(class, f, _NO_RESCALING)
    rescale(e) = (a = σ(e.from); b = σ(e.to);
                  e.rapid_equilibrium && a == b ? (_NO_RESCALING, _NO_RESCALING) : (a, b))
    rescaling = Dict{Int, Tuple{_GaugeClass, _GaugeClass}}()
    for (h, ks) in enumerate(members)
        isempty(ks) && continue
        r = rescale(m.edges[first(ks)])
        all(k -> rescale(m.edges[k]) == r, ks) || return nothing
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
state unless it is `:OnlyA`. Each state's graph is its law's (`_active_mask`,
`_inactive_mask`): the edges its tags keep between forms they join to free E. The copy
binds in a state only at the sites that graph keeps; in a state where it binds nothing,
every factor of that state is 1. A copy that binds in neither state is not redundant: its
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
    binds_copy(ks) = (i = Mechanisms._bound_entry(cm.edges[first(ks)]);
                      i != 0 && cm.table[i].modifier)
    any(binds_copy, members) || return Int[]
    if m isa Mechanisms.Mechanism
        twin = _productive_twin(cm, trues(length(cm.edges)))
        return [g for g in eachindex(members)
                if _gauge_rescaling(cm, members, g, twin, _one_state_class) !== nothing]
    end
    tags = Mechanisms.cat_allo_states(m)
    masks = (Mechanisms._active_mask(m), Mechanisms._inactive_mask(m))
    kept = [[filter(k -> mask[k], ks) for ks in members] for mask in masks]
    twins = [_productive_twin(cm, mask) for mask in masks]
    unscaled(c) = Dict(h => (_NO_RESCALING, _NO_RESCALING)
                       for h in eachindex(kept[c]) if !isempty(kept[c][h]))
    filter(collect(eachindex(members))) do g
        all(c -> isempty(kept[c][g]), 1:2) && return false
        one_number(h) = h > 0 && tags[h] === :EqualAI && tags[g] === :EqualAI
        label(state) = (h, t) -> _GaugeClass(h, t, one_number(h) ? 0 : state)
        rescaling(c) = isempty(kept[c][g]) ? unscaled(c) :
                       _gauge_rescaling(cm, kept[c], g, twins[c], label(c))
        active = rescaling(1)
        active === nothing && return false
        inactive = rescaling(2)
        inactive === nothing && return false
        all(h -> tags[h] !== :EqualAI || !haskey(inactive, h) || active[h] == inactive[h],
            keys(active))
    end
end

"""
The forms of graph `m` where table entry `entry` is free on some edge: the edge's `from`
when it takes the entry up, its `to` when it gives it off, so a Theorell–Chance step
E(A) + B → E(Q) + P frees B at E(A) and P at E(Q). A reactant and its modifier copy are
two entries that share a name but bind different sites: the reactant's entry asks where
the metabolite binds productively, the copy's where it binds its dead-end site. The
dead-end move targets these forms for competition with the metabolite. Entry 0, which
`_entry` gives for an entry `m` lacks, is free nowhere.
"""
function _forms_where_free(m::Mechanisms.Mechanism, entry::Int)
    free = falses(length(m.forms))
    entry == 0 && return free
    for e in m.edges
        Mechanisms.copies(e.consumed, entry) > 0 && (free[e.from] = true)
        Mechanisms.copies(e.released, entry) > 0 && (free[e.to] = true)
    end
    free
end

"""
The dead-end move's raw child of graph `m`: `m` as raw input (`_raw`) with modifier `reg`
also bound at each form of `sites`, a complex in the site's conformation and residual;
each edge between two sites mirrored between their complexes, in the edge's own group;
and one new group, after `m`'s, of the rapid-equilibrium bindings of `reg` at the sites.
`reg` joins the declared modifiers, and the child's table gains its entry.
"""
function _dead_end_raw(m::Mechanisms.Mechanism, reg::Symbol, sites::Vector{Int})
    raw = Mechanisms._raw(m)
    inhibitor = Mechanisms.Met(reg, true)
    complex = Dict(f => Mechanisms._raw_form(m, f, inhibitor) for f in sites)
    for (k, e) in enumerate(m.edges)
        haskey(complex, e.from) && haskey(complex, e.to) || continue
        push!(raw.steps, Mechanisms.RawStep(complex[e.from], complex[e.to],
                                            raw.steps[k].consumed, raw.steps[k].released,
                                            e.rapid_equilibrium))
        push!(raw.group, m.group[k])
    end
    for f in sites
        push!(raw.steps, Mechanisms.RawStep(Mechanisms._raw_form(m, f), complex[f],
                                            [inhibitor], Mechanisms.Met[], true))
        push!(raw.group, Mechanisms.ngroups(m) + 1)
    end
    sort!(unique!(push!(raw.modifiers, reg)))
    raw
end

"""
The dead-end move's child of `m` that `raw` (`_dead_end_raw`) writes, binding modifier
`reg` in a new group; `nothing` when the copy is redundant in it (`_redundant_copy_groups`).
A `Mechanism` candidate is canonicalized (`_canonicalize`) and its new group tested by
`_gauge_rescaling` on the candidate's own twins (`_productive_twin`), before any check
runs; it is built through the checked constructor only when kept. The twins are the
candidate's, not the parent's, though they are the same forms: the complexes are not
productive, and their edges attach each complex to its site's segment and join no two
segments, so every form of `m` keeps its segment and offsets; but the complexes sort in
among `m`'s forms and the copy's table entry moves every later one, so `m`'s form and entry
indices do not hold in the candidate. Only the new group is tested: `m` holds no redundant
copy group (the parent rule), and in each of about 218,000 enumerated parents checked, the
placement made no older copy redundant (an empirical finding, not a proof). An allosteric
child is built through the checked constructor, its copy's group tagged `:EqualAI`, with
`m`'s multiplicity and regulatory sites, and runs `_redundant_copy_groups` only when every
complex of the copy has a twin in the active state (`_all_twin`), which redundancy needs.
"""
function _dead_end_child(m::Mechanisms.Mechanism, raw::Mechanisms.RawMechanism, ::Symbol)
    c = Mechanisms._canonicalize(raw)
    cm = c.mechanism
    twin = _productive_twin(cm, trues(length(cm.edges)))
    _gauge_rescaling(cm, Mechanisms._group_members(cm),
                     c.group_of[Mechanisms.ngroups(m) + 1], twin,
                     _one_state_class) === nothing || return nothing
    Mechanisms.Mechanism(raw)
end
function _dead_end_child(am::Mechanisms.AllostericMechanism, raw::Mechanisms.RawMechanism,
                         reg::Symbol)
    child = Mechanisms.AllostericMechanism(
        Mechanisms.RawAllostericMechanism(raw, [am.cat_allo_states; :EqualAI],
                                          am.catalytic_multiplicity,
                                          Mechanisms._raw_sites(am)))
    cm = child.catalytic
    # The copy's own bindings, the only edges that take its entry up.
    i = Mechanisms._entry(cm, Mechanisms.Met(reg, true))
    bindings = findall(e -> Mechanisms.copies(e.consumed, i) > 0, cm.edges)
    _all_twin(cm, bindings, _productive_twin(cm, trues(length(cm.edges)))) &&
        !isempty(_redundant_copy_groups(child)) ? nothing : child
end

"""
    _expand_add_dead_end_regulator(m::Union{Mechanism, AllostericMechanism},
                                   rxn::EnzymeReaction) → Vector{typeof(m)}

Add a dead-end regulator binding step set. For each `CompetitiveInhibitor`
declared in `rxn` but not yet bound by `m`'s steps, enumerate inhibitor
competition patterns (S × P × existing inhibitors); for each pattern,
add RE binding steps to the forms where a competing metabolite is free, unless
the form already holds a competing metabolite or holds every substrate or every
product. Competition is decided per site, not per name. A substrate or product
declared as a competitive inhibitor is a copy that binds a dead-end site of its
own, and that site and the reactant's catalytic site are different sites by
definition. Competition with a substrate or product M targets the forms where
M is free on a productive step and excludes the forms that hold M productively;
competition with an existing inhibitor I targets the forms where I binds as a
modifier (its binding steps' source forms) and excludes the forms that hold I; the
capacity test counts productive bindings only. A form that carries only the copy of M
is therefore a site for an inhibitor that competes with M. Mirror steps inherit their
catalytic counterpart's kinetic group. All new binding steps for a single regulator
share one fresh trailing kinetic group (one K_R parameter). Each child is the parent's
raw input with the copy added (`_dead_end_raw`).
A pattern is skipped when its child binds the new copy redundantly
(`_redundant_copy_groups`). A `Mechanism` candidate is tested by `_gauge_rescaling`
inside `_dead_end_child`, before the child is built; an allosteric child is built and
tested directly. Ligands of an allosteric mechanism's regulatory sites are not eligible:
they belong on a regulatory site.

The caller must pass the declared `rxn`: its regulators decide which inhibitors are
eligible.
"""
function _expand_add_dead_end_regulator(m::Union{Mechanisms.Mechanism,
                                                  Mechanisms.AllostericMechanism},
                                        rxn::Mechanisms.EnzymeReaction)
    cm = Mechanisms.catalytic(m)
    sub_names = Set(Mechanisms.name(s) for s in Mechanisms.substrates(rxn))
    prod_names = Set(Mechanisms.name(p) for p in Mechanisms.products(rxn))
    modifiers = Mechanisms._modifier_mask(cm)
    bound_regs = _bound_comp_inhibitors(m)
    excluded = union(bound_regs, _bound_allo_regs(m))
    eligible_regs = [Mechanisms.name(Mechanisms.regulator(rm))
                     for rm in Mechanisms.regulators(rxn)
                     if Mechanisms.regulator(rm) isa Mechanisms.CompetitiveInhibitor &&
                        Mechanisms.name(Mechanisms.regulator(rm)) ∉ excluded]
    isempty(eligible_regs) && return typeof(m)[]

    # Per form, the names bound productively and the names bound as modifiers: a
    # reactant's copy sits at a dead-end site, not at its reactant's.
    productive = [Set(Mechanisms._entry_names(cm, f.bound & ~modifiers)) for f in cm.forms]
    inhibiting = [Set(Mechanisms._entry_names(cm, f.bound & modifiers)) for f in cm.forms]

    # The sites of each competition pattern in form order, once each in pattern order.
    placements = Vector{Int}[]
    for (comp_subs, comp_prods, comp_inhibitors) in
            _inhibitor_competition_patterns(sub_names, prod_names, collect(bound_regs))
        comp_reactants = union(comp_subs, comp_prods)
        target_forms = falses(length(cm.forms))
        for met in comp_reactants
            target_forms .|=
                _forms_where_free(cm, Mechanisms._entry(cm, Mechanisms.Met(met, false)))
        end
        for inh in comp_inhibitors
            target_forms .|=
                _forms_where_free(cm, Mechanisms._entry(cm, Mechanisms.Met(inh, true)))
        end
        push!(placements, [f for f in findall(target_forms)
                           if !issubset(sub_names, productive[f]) &&
                              !issubset(prod_names, productive[f]) &&
                              isempty(intersect(productive[f], comp_reactants)) &&
                              isempty(intersect(inhibiting[f], comp_inhibitors))])
    end
    unique!(filter!(!isempty, placements))

    results = typeof(m)[]
    for reg in eligible_regs, sites in placements
        child = _dead_end_child(m, _dead_end_raw(cm, reg, sites), reg)
        child === nothing || push!(results, child)
    end
    results
end

"""
    _admits(m::Mechanism, tags::Vector{Symbol}) -> Bool

Whether admission (`_admission`) accepts `m` with the tags `tags`, one per kinetic group: a
path to the limit of the tagged constants exists, every form free E reaches only through
tagged steps lies behind a tagged plain binding, the active conformation turns over, and
the face condition holds on every face. It judges one catalytic subunit with no regulatory
site, which admission never reads. Class 3 is not run: the guard at the enumerator's exits
drops a child that fails it.
"""
_admits(m::Mechanisms.Mechanism, tags::Vector{Symbol}) =
    _admission(Mechanisms.AllostericMechanism(m, tags, 1,
                                              Mechanisms.RegulatorySite[])).fault ===
        _NO_FAULT

"""
    _stops_inactive_turnover(m::Mechanism, tags::Vector{Symbol}) -> Bool

Whether the tags `tags` leave the inactive conformation of `m` no catalytic cycle: no
cycle of the edges its law keeps (`_inactive_mask`) turns the reaction over
(`_has_turnover_cycle`).
"""
_stops_inactive_turnover(m::Mechanisms.Mechanism, tags::Vector{Symbol}) =
    !_has_turnover_cycle(m,
                         Mechanisms._inactive_mask(Mechanisms.AllostericMechanism(
                             m, tags, 1, Mechanisms.RegulatorySite[])))

"""
    _tagging_bases(m::Mechanism) -> Vector{Vector{Symbol}}

The tag sets the tagging move (`_expand_to_allosteric`) builds on, one tag per kinetic
group of `m` (spec section 7, "Allosteric ping-pong"). The first candidate tags every
chemistry group (`_chemistry_groups`) `:OnlyA` and every binding group `:EqualAI`; when
admission accepts it (`_admits`) it is the only base. Otherwise the bases are the
smallest sets of chemistry groups, by inclusion, that stop every catalytic cycle of the
inactive conformation (`_stops_inactive_turnover`) and pass admission, each tagged
`:OnlyA` with every other group `:EqualAI`. A ping-pong mechanism thus usually gets one
base per half-reaction, the other half running in both conformations: a half whose
chemistry groups fail admission gives none, and a half spread over several chemistry
groups can give several. The bases come by size, each size in `_subsets` order.
"""
function _tagging_bases(m::Mechanisms.Mechanism)
    chem = Mechanisms._chemistry_groups(m)
    tagged(onlya) = [g in onlya ? :OnlyA : :EqualAI for g in eachindex(chem)]
    groups = findall(chem)
    _admits(m, tagged(groups)) && return [tagged(groups)]
    sets = _subsets(groups)
    found = Vector{Int}[]
    for k in 1:length(groups), s in sets
        length(s) == k && !any(b -> b ⊆ s, found) || continue
        t = tagged(s)
        _stops_inactive_turnover(m, t) && _admits(m, t) && push!(found, s)
    end
    [tagged(s) for s in found]
end

"""
    _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

Convert a non-allosteric `Mechanism` into allosteric variants, each built on a base of
`_tagging_bases`: the tag set with every chemistry group `:OnlyA` and every binding group
`:EqualAI` when admission accepts it, else each smallest set of chemistry groups that stops
the inactive conformation's catalysis and passes admission. A chemistry group is one
holding a chemistry step (`_is_chemistry`): an isomerization, a fused binding or a
Theorell–Chance step; every other group is a binding group. On a ping-pong mechanism,
tagging both half-reactions strands the covalent forms between them, so a base usually
tags one half alone and the other half runs in both conformations (spec section 7,
"Allosteric ping-pong"); a half whose chemistry groups fail admission gives no base, and a
half spread over several chemistry groups can give several. Two kinds are not enumerated:
the all-`:EqualAI` baseline and a V-type variant without a regulator.

  * The all-`:EqualAI` baseline is never emitted — the two conformations
    are identical, `L` cancels, and the mechanism is indistinguishable
    from `m`.
  * K-type: a nonempty set of binding groups set `:OnlyA` on top of a base, and on top
    of the every-chemistry set when admission rejects that set, since a tagged binding
    can restore its limit; each is emitted bare and kept when admission accepts it
    (`_admits`). Over more than one catalytic subunit the conformational equilibrium
    enters the rate with the multiplicity's power. Over one subunit the inactive
    conformation is a dead-end branch of the free enzyme, and `L` is a phantom in some
    K-type variants and not in others. In every K-type variant of the rapid-equilibrium
    ordered bi-bi seed a rescaling of the active state's constants absorbs it, whatever
    that conformation binds. In the merged ordered bi-bi whose B group mixes a fused and
    a plain binding it is a phantom only in the variants whose inactive conformation
    binds nothing: no rescaling absorbs an A or Q binding kept there. Binding nothing is
    not enough: in a random-order merged bi-bi, whose A group takes A up at E and at
    E(B) under one rate constant, `L` shows in every K-type variant. identify keeps every
    K-type variant in the frontier, since its children can carry laws its parent's children
    lack. A variant that `_is_phantom` (`src/key.jl`) flags is not fitted when another
    mechanism holds the fit of its law at `L = 0`: it reads that fit, with `L = 0`.
  * V-type: a base, with no `:OnlyA` binding, paired with a declared allosteric regulator at
    a new site, one variant per `(regulator, tag)` with `tag ∈ {:OnlyA, :OnlyI}`. On the
    every-chemistry base the inactive state binds substrate and product as the active state
    does and runs no chemistry, so without a regulator `L` folds entirely into `kcat`
    (`v = kcat/(1+L)·shape`) and is not observable. On a one-half base the inactive
    conformation runs the other half-reaction, and `L` shows without a regulator too, yet
    every base is emitted only with one. A reaction with no declared allosteric regulators
    emits no V-type.

Admission does not run class 3 of the guard: a child that fails a face is dropped where
`expand_mechanisms` returns it. A ping-pong parent with no declared regulator can end
with no child, when every K-type strands a form or a basin.

For each value in `rxn`'s `allowed_catalytic_multiplicities`, the move emits
the K-types, then the V-types, at that `catalytic_multiplicity`. Every variant
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
function _expand_to_allosteric(m::Mechanisms.Mechanism, rxn::Mechanisms.EnzymeReaction)
    hyperbolic = _hyperbolic_catalysis(m)
    cns = [cn for cn in Mechanisms.allowed_catalytic_multiplicities(rxn)
           if cn == 1 || hyperbolic]
    isempty(cns) && return Mechanisms.AllostericMechanism[]
    chem = Mechanisms._chemistry_groups(m)
    vtags = [c ? :OnlyA : :EqualAI for c in chem]
    regs = [Mechanisms.name(Mechanisms.regulator(rm)) for rm in Mechanisms.regulators(rxn)
            if Mechanisms.regulator(rm) isa Mechanisms.AllostericRegulator]
    bases = _tagging_bases(m)
    # K-type: a nonempty set of binding groups :OnlyA on top of a base, and on top of the
    # every-chemistry set when admission rejects it; admission drops the rest.
    kbases = vtags in bases ? bases : [bases; [vtags]]
    ktags = Vector{Symbol}[]
    for base in kbases, sel in _subsets(findall(!, chem))[2:end]
        tags = copy(base)
        tags[sel] .= :OnlyA
        _admits(m, tags) && push!(ktags, tags)
    end
    results = Mechanisms.AllostericMechanism[]
    for cn in cns
        for tags in ktags
            push!(results,
                  Mechanisms.AllostericMechanism(m, tags, cn, Mechanisms.RegulatorySite[]))
        end
        # V-type: a base paired with a declared regulator, whose binding reveals L.
        any(chem) || continue
        for base in bases, reg in regs, tag in (:OnlyA, :OnlyI)
            push!(results,
                  Mechanisms.AllostericMechanism(
                      m, base, cn, [Mechanisms.RegulatorySite([reg], cn, [tag])]))
        end
    end
    filter(_hyperbolic_where_required, results)
end

"""
    _expand_add_allosteric_regulator(am::AllostericMechanism,
                                     rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

Add one `AllostericRegulator` declared in `rxn` not yet bound by `am`
at a regulatory site. A name already bound as a catalytic-step dead-end
competitive inhibitor is still eligible: a metabolite may play both
roles, and the two render to distinct parameters. For each (new ligand,
target site, tag) combination, emit a variant:

  * target site ∈ {new site} ∪ {existing sites}
  * tag ∈ {:OnlyA, :OnlyI, :NonequalAI} for any target site
  * tag = :EqualAI only at an existing site that already has at least
    one non-`:EqualAI` ligand. A site whose ligands are all `:EqualAI`
    binds both conformations alike and cancels from the rate equation;
    this guard, not the `RegulatorySite` constructor, keeps the move from
    building one.

New sites inherit `am.catalytic_multiplicity` as their multiplicity.
The mechanism's catalytic side and regulatory_sites' multiplicities pass
through unchanged.

The regulators declared in `rxn` decide which regulators are eligible.
"""
function _expand_add_allosteric_regulator(
    am::Mechanisms.AllostericMechanism, rxn::Mechanisms.EnzymeReaction,
)
    taken = _bound_allo_regs(am)
    sites = Mechanisms.regulatory_sites(am)
    results = Mechanisms.AllostericMechanism[]
    for rm in Mechanisms.regulators(rxn)
        reg = Mechanisms.regulator(rm)
        reg isa Mechanisms.AllostericRegulator && !(Mechanisms.name(reg) in taken) ||
            continue
        # Non-:EqualAI tags at any (new or existing) site.
        for tag in (:OnlyA, :OnlyI, :NonequalAI), site_idx in 0:length(sites)
            # Appending to an existing site whose conformations are disjoint
            # from the new ligand's (an all-:OnlyA site gaining an :OnlyI
            # ligand, or the reverse) reproduces the add-at-a-new-site
            # equation — redundant, so skip it.
            site_idx >= 1 && isempty(intersect(_state_conformations(tag),
                _site_active_states(sites[site_idx]))) && continue
            push!(results, _make_am_with_added_reg(am, Mechanisms.name(reg), tag, site_idx))
        end
        # :EqualAI only at an existing site that keeps a non-:EqualAI ligand: an
        # all-:EqualAI site binds both conformations alike and cancels from the rate.
        for site_idx in findall(site -> any(!=(:EqualAI), Mechanisms.allo_states(site)),
                                sites)
            push!(results,
                  _make_am_with_added_reg(am, Mechanisms.name(reg), :EqualAI, site_idx))
        end
    end
    results
end

"""
Build an `AllostericMechanism` identical to `am` except the ligand
`reg::Symbol` is added at `site_idx` (0 = create a new site,
1..length(am.regulatory_sites) = append to that existing site) with
allosteric state `tag`. Multiplicity for a new site inherits
`am.catalytic_multiplicity`.
"""
function _make_am_with_added_reg(
    am::Mechanisms.AllostericMechanism, reg::Symbol, tag::Symbol, site_idx::Int,
)
    sites = copy(Mechanisms.regulatory_sites(am))
    if site_idx == 0
        push!(sites, Mechanisms.RegulatorySite([reg], Mechanisms.catalytic_multiplicity(am),
                                               [tag]))
    else
        s = sites[site_idx]
        sites[site_idx] =
            Mechanisms.RegulatorySite([Mechanisms.ligands(s); reg],
                                      Mechanisms.multiplicity(s),
                                      [Mechanisms.allo_states(s); tag])
    end
    Mechanisms._with(am; sites)
end

"""
    _expand_change_allo_state(am::AllostericMechanism)
        → Vector{AllostericMechanism}

Relax a "constrained" tag (`:EqualAI`, `:OnlyA`, `:OnlyI`) to `:NonequalAI`. A binding
catalytic group and a regulatory ligand each relax individually — one variant per group
not already `:NonequalAI`. The chemistry groups (those holding a chemistry step,
`_is_chemistry`) relax **together**, one variant setting every non-`:NonequalAI`
chemistry group to `:NonequalAI` at once, so a fully live inactive conformation is one
move from every V-type, the tagging move's one-half bases included. The base catalytic
steps, multiplicity, and untouched tags are preserved.

The chemistry relaxation is tried only while no binding group is `:OnlyA`. Beside an
`:OnlyA` binding the restored chemistry can join the binding's two sides in the inactive
conformation and leave its affinity no path to zero; where the inactive conformation
binds nothing (every binding `:OnlyA`), the skipped variant is rate-equivalent to the
form emitted directly, so no observable hypothesis is lost.

A catalytic relaxation is kept when admission accepts it (`_admits`): relaxing an
`:OnlyA` binding joins its two sides, which can leave another tagged binding inside one
component or close a cycle of tagged bindings. A regulatory ligand's tag does not enter
admission — a regulator site completes no catalytic cycle — so that branch needs no
filter. Admission does not run class 3 of the guard: a child that fails a face is dropped
where `expand_mechanisms` returns it. Above one subunit a child is kept only when each
conformation's law passes `_hyperbolic_catalysis` (`_requires_hyperbolic_catalysis`):
relaxing an `:OnlyA` abortive binding brings its complex into the inactive law, which can
then carry a power.

A child with a redundant competitive-inhibitor group (`_redundant_copy_groups`)
is dropped. The per-state test reads the tags this move changes: relaxing an
`:OnlyA` binding to `:NonequalAI` brings its complex into the inactive state,
where a copy that was new there can become its twin with a consistent gauge, and
the copy's constant then shows in neither state.
"""
function _expand_change_allo_state(am::Mechanisms.AllostericMechanism)
    states = Mechanisms.cat_allo_states(am)
    chem = Mechanisms._chemistry_groups(am.catalytic)
    # Binding catalytic groups relax individually; the chemistry groups relax together,
    # in one variant, and only while no binding group is :OnlyA.
    relaxations = [[g] for g in eachindex(states) if !chem[g] && states[g] != :NonequalAI]
    any(chem .& (states .!= :NonequalAI)) && !any(.!chem .& (states .== :OnlyA)) &&
        push!(relaxations, findall(chem))
    results = Mechanisms.AllostericMechanism[]
    for gs in relaxations
        new_states = copy(states)
        new_states[gs] .= :NonequalAI
        _admits(am.catalytic, new_states) || continue
        push!(results, Mechanisms._with(am; states = new_states))
    end
    for (si, site) in enumerate(Mechanisms.regulatory_sites(am)),
        li in eachindex(Mechanisms.ligands(site))
        Mechanisms.allo_states(site)[li] == :NonequalAI && continue
        new_sites = copy(Mechanisms.regulatory_sites(am))
        new_states = copy(Mechanisms.allo_states(site))
        new_states[li] = :NonequalAI
        new_sites[si] = Mechanisms.RegulatorySite(Mechanisms.ligands(site),
                                                  Mechanisms.multiplicity(site), new_states)
        push!(results, Mechanisms._with(am; sites = new_sites))
    end
    filter!(c -> isempty(_redundant_copy_groups(c)), results)
    filter!(_hyperbolic_where_required, results)
end

"""
    _state_conformations(state::Symbol) -> Set{Symbol}

The conformations a single ligand's allosteric `state` acts on: `:active` for
`:OnlyA`/`:EqualAI`/`:NonequalAI`, `:inactive` for `:OnlyI`/`:EqualAI`/`:NonequalAI`.
"""
function _state_conformations(state::Symbol)
    conf = Set{Symbol}()
    Mechanisms._excluded_by(state, false) || push!(conf, :active)
    Mechanisms._excluded_by(state, true) || push!(conf, :inactive)
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
    for st in Mechanisms.allo_states(site)
        union!(active, _state_conformations(st))
    end
    active
end

"""
    _expand_merge_regulatory_sites(am::AllostericMechanism)
        → Vector{AllostericMechanism}

Merge each unordered pair of regulatory sites into one shared site holding
both sites' ligands, at no parameter cost (Δ0 — every ligand keeps its own
dissociation constant). This lets two regulators compete at one site,
including the activator↔antagonist ambiguity. For the merged ligands,
enumerate the Δ0-valid allo-state assignments:

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
`catalytic_multiplicity`), and the `RegulatorySite` constructor sorts its
ligands by name, so two merge routes reaching the same ligand partition produce
`==` mechanisms and dedup by `hash`. Regulator type is not enforced here;
`expand_mechanisms` drops every child that fails `_respects_reg_type`.
"""
function _expand_merge_regulatory_sites(am::Mechanisms.AllostericMechanism)
    sites = Mechanisms.regulatory_sites(am)
    n = length(sites)
    results = Mechanisms.AllostericMechanism[]
    for i in 1:(n - 1), j in (i + 1):n
        ligs = vcat(Mechanisms.ligands(sites[i]), Mechanisms.ligands(sites[j]))
        base_states = vcat(Mechanisms.allo_states(sites[i]),
                           Mechanisms.allo_states(sites[j]))
        mult = Mechanisms.multiplicity(sites[i])
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
            push!(results, Mechanisms._with(am; sites = vcat(others, [merged])))
        end
    end
    results
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

"""
    expand_mechanisms(mechs, reaction) -> Vector{Union{Mechanism, AllostericMechanism}}

Apply all expansion moves (RE→SS, split kinetic group, add dead-end
regulator, to-allosteric, add allosteric regulator, change allo state,
merge regulatory sites) to each input mechanism and return the children as
a flat vector. Bucketing by parameter count is the caller's job, not
enumeration's. Each parent must satisfy the two emission rules
(`_assert_emission_rules`): every steady-state group carries flux and no
competitive-inhibitor group is redundant (`_redundant_copy_groups`). Every child that
passes the regulator-type filter must conserve atoms (`_assert_atom_conserving`), which
errors on the first child that does not; a child whose rate law the guard rejects
(`_is_well_posed`) is dropped.
"""
function expand_mechanisms(
    mechs::Vector{<:Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}},
    rxn::Mechanisms.EnzymeReaction)
    result = Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}[]
    for m in mechs
        _assert_emission_rules(m)
        append!(result, _expand_re_to_ss(m), _expand_split_kinetic_group(m),
                _expand_add_dead_end_regulator(m, rxn))
        if m isa Mechanisms.Mechanism
            append!(result, _expand_to_allosteric(m, rxn))
        else
            append!(result, _expand_add_allosteric_regulator(m, rxn),
                    _expand_change_allo_state(m), _expand_merge_regulatory_sites(m))
        end
    end
    filter!(c -> _respects_reg_type(c, rxn), result)
    foreach(child -> _assert_atom_conserving(child, rxn), result)
    filter!(_is_well_posed, result)
    result
end

# ─── Entry Points ────────────────────────────────────────────

"""
    init_mechanisms(reaction::EnzymeReaction) -> Vector{Mechanism}

Public entry point. Produces the starting mechanisms of the search for a
reaction as concrete `Mechanism` structs: first the seeds, then their merged
and Theorell–Chance variants. A seed is built for each catalytic topology
(`_catalytic_topologies`) and each substrate/product dead-end subset
(`_expand_substrate_product_dead_ends`), with one steady-state step, a
chemistry isomerization, and the binding steps of one metabolite at one RE/SS
flag collapsed into one kinetic group (`_seed_groups`). Dead-end
enumeration respects `shared_catalytic_site`. The variants of each seed
(`_seed_variants`) follow in the seeds' order, each once. The parameter counts
are mixed: a three-group merged variant or a decorated Theorell–Chance variant
fits more parameters than its seed, a merged ping-pong variant fewer. Errors when
the reaction has no catalytic topology. A seed whose rate law the guard rejects
(`_is_well_posed`) is dropped, and so is a variant with no turnover cycle
(`_has_turnover_cycle`): the flip search has run the guard's classes 2 and 3 on every
variant (`_faces_pass`).
"""
function init_mechanisms(r::Mechanisms.EnzymeReaction)
    subs = Symbol[Mechanisms.name(s) for s in Mechanisms.substrates(r)]
    prods = Symbol[Mechanisms.name(p) for p in Mechanisms.products(r)]
    seeds = _expand_substrate_product_dead_ends(_catalytic_topologies(r), r)
    mechs = [Mechanisms.Mechanism(Mechanisms.RawMechanism(subs, prods, Symbol[], steps,
                                                          _seed_groups(steps)))
             for steps in seeds]
    foreach(m -> _assert_atom_conserving(m, r), mechs)
    seen = Set(mechs)
    out = filter(_is_well_posed, mechs)
    for m in unique(mechs), v in _seed_variants(m, r)
        v in seen && continue
        push!(seen, v)
        _has_turnover_cycle(v, trues(length(v.edges))) && push!(out, v)
    end
    out
end

"""
    seed_mechanisms(rxn, required_allo::Set{Symbol}, required_comp::Set{Symbol})
        -> Vector{Union{Mechanism, AllostericMechanism}}

Fully-required seed set for the beam. Grows `init_mechanisms(rxn)` by a
breadth-first closure under the seed-build structure moves and retains the nodes
that bind every required regulator. The two allosteric-lifting moves
(`_expand_to_allosteric`, `_expand_add_allosteric_regulator`) run only when an
allosteric regulator is required; `_expand_add_dead_end_regulator` always runs.
A competitive-only required set therefore stays non-allosteric — the seeds are
`Mechanism`s at `base + n_required_comp`, with no `L`.

A child is enqueued (and marked visited by `hash`) only when it is a valid seed
node:

1. no `:NonequalAI` tag anywhere (cheap states only; the beam reaches
   `:NonequalAI` later via `change_allo_state`);
2. every regulatory site binds a single ligand — one site per required
   allosteric regulator. This bounds the closure: multi-ligand children (from
   adding a regulator to an existing site) fail it and are dropped before
   expansion;
3. no optional regulator is bound — every bound allosteric ligand ∈
   `required_allo` and every bound competitive inhibitor ∈ `required_comp`;
4. every ligand respects its declared type (`_respects_reg_type`).

The seeds are the valid nodes that additionally bind ALL of `required_allo` at
regulatory sites and ALL of `required_comp` as competitive-inhibitor dead ends.
Returned deduped (each node is visited once). Explored as a wave-parallel BFS:
each level's child generation is distributed via `pmap`, then deduped and
folded back into the next frontier in frontier order — the same order the
serial FIFO BFS would enqueue in — so the result is byte-identical to a serial
traversal. Errors when no mechanism binds every required regulator; a seed whose rate
law the guard rejects (`_is_well_posed`) is dropped. The search still expands such a node:
only the seeds it returns are checked.
"""
function seed_mechanisms(rxn::Mechanisms.EnzymeReaction, required_allo::Set{Symbol},
                         required_comp::Set{Symbol})
    visited = Set{UInt64}()
    seeds = Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}[]
    # Dedup + collect on the main node. Returns true when `m` is new, so the
    # caller advances only genuinely-new nodes to the next wave. Called in
    # frontier order, which equals the serial BFS enqueue order, so `visited`
    # and `seeds` end byte-identical to the FIFO version.
    consider!(m) = begin
        h = hash(m)
        h in visited && return false
        push!(visited, h)
        issubset(required_allo, _bound_allo_regs(m)) &&
            issubset(required_comp, _bound_comp_inhibitors(m)) && push!(seeds, m)
        true
    end
    frontier = Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}[
        m for m in init_mechanisms(rxn) if consider!(m)]
    while !isempty(frontier)
        # Per-node child generation is pure and independent — distribute it.
        childsets = pmap(frontier) do m
            children = Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}[]
            isempty(required_allo) || append!(children, m isa Mechanisms.Mechanism ?
                _expand_to_allosteric(m, rxn) : _expand_add_allosteric_regulator(m, rxn))
            append!(children, _expand_add_dead_end_regulator(m, rxn))
            filter!(c -> _is_seed_node(c, rxn, required_allo, required_comp), children)
        end
        next = Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}[]
        for cs in childsets, c in cs
            consider!(c) && push!(next, c)
        end
        frontier = next
    end
    if isempty(seeds)
        list(s) = isempty(s) ? "none" : join(sort!(collect(s)), ", ")
        error("seed_mechanisms: no mechanism binds every required regulator " *
              "(competitive inhibitors: " * list(required_comp) *
              "; allosteric regulators: " * list(required_allo) * "). A substrate or " *
              "product declared as a competitive inhibitor binds only where its complex " *
              "is not redundant; a uni-uni seed has no such site. Mark a regulator " *
              "optional with `optional_competitive_inhibitors` or " *
              "`optional_allosteric_regulators`, or remove its declaration")
    end
    filter!(_is_well_posed, seeds)
    seeds
end

"""Names of the allosteric regulators bound at `m`'s regulatory sites (empty for
a `Mechanism`)."""
_bound_allo_regs(::Mechanisms.Mechanism) = Set{Symbol}()
_bound_allo_regs(am::Mechanisms.AllostericMechanism) =
    Set{Symbol}(lig for site in Mechanisms.regulatory_sites(am)
                for lig in Mechanisms.ligands(site))

"""Names of the competitive inhibitors some step of `m` takes up as a modifier, collected
group by group in edge order (`_flat_edges`)."""
function _bound_comp_inhibitors(m::Union{Mechanisms.Mechanism,
                                          Mechanisms.AllostericMechanism})
    cm = Mechanisms.catalytic(m)
    modifiers = Mechanisms._modifier_mask(cm)
    Set{Symbol}(cm.table[i].name for k in Mechanisms._flat_edges(cm)
                for i in Mechanisms.entries(cm.edges[k].consumed & modifiers))
end

"""
    _is_seed_node(m, rxn::EnzymeReaction, required_allo::Set{Symbol},
                 required_comp::Set{Symbol}) -> Bool

A child worth expanding: cheap states, one ligand per regulatory site, no
optional regulator bound, and reg-type-respecting.
"""
function _is_seed_node(m::Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism},
                       rxn::Mechanisms.EnzymeReaction, required_allo::Set{Symbol},
                       required_comp::Set{Symbol})
    m isa Mechanisms.AllostericMechanism && (:NonequalAI in Mechanisms.cat_allo_states(m) ||
        any(site -> length(Mechanisms.ligands(site)) != 1 ||
                    only(Mechanisms.allo_states(site)) == :NonequalAI,
            Mechanisms.regulatory_sites(m))) && return false
    issubset(_bound_allo_regs(m), required_allo) || return false
    issubset(_bound_comp_inhibitors(m), required_comp) || return false
    _respects_reg_type(m, rxn) || return false
    true
end
