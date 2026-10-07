# ABOUTME: Mechanism enumeration by incremental parameter count growth
# ABOUTME: Provides init_mechanisms, expand_mechanisms building blocks


# ─── Catalytic topologies ─────────────────────────────────────────────

"""
Net atom multiset of the reactants named in `plus` minus those named in `minus`, with
zero counts dropped. A name may repeat; a name that names no reactant carries no atoms.
"""
function _net_atoms(reaction::EnzymeReaction, plus, minus)
    acc = Dict{Symbol,Int}()
    for ra in reactants(reaction)
        k = count(==(name(metabolite(ra))), plus) - count(==(name(metabolite(ra))), minus)
        for (a, c) in atoms(ra)
            acc[a] = get(acc, a, 0) + k * c
        end
    end
    filter!(kv -> kv.second != 0, acc)
end

# ─── Atom-conservation validation ────────────────────────────

"""
Assert one `Step` conserves atoms: the atoms of `from_species` plus those of
the consumed metabolites must equal the atoms of `to_species` plus those of the
released metabolites (an iso step leaves the atom multiset unchanged). A form
carries the atoms of its bound metabolites and of its residual's added substrates,
less those of its residual's subtracted products, read from the reaction's
inventory by name (`_net_atoms`). Errors naming the offending step.
"""
function _assert_step_atom_conserving(reaction::EnzymeReaction, s::Step)
    # The names counted on one side of `s`: `form`'s bound metabolites and residual
    # additions, the residual subtractions of the `other` form, and the free metabolites.
    side(form, other, free) = [name.(bound(form)); name.(added(residual(form)));
                               name.(subtracted(residual(other))); name.(free)]
    diff = _net_atoms(reaction, side(to_species(s), from_species(s), released(s)),
                      side(from_species(s), to_species(s), consumed(s)))
    isempty(diff) || error(
        "atom-non-conserving step $(name(from_species(s))) → $(name(to_species(s))) " *
        "(consumed $(name.(consumed(s))), released $(name.(released(s)))): " *
        "atoms(from) + atoms(consumed) − atoms(to) − atoms(released) = $diff")
    nothing
end

"""
    _assert_emission_rules(m)

The two rules every mechanism the moves emit satisfies, checked on a parent before
it is expanded. Every steady-state kinetic group holds a step that carries net flux
(`_flux_carrying_groups` on `steps(m)`, which for an allosteric mechanism is its
active-state graph). No kinetic group binds a competitive inhibitor redundantly
(`_redundant_copy_groups`). A parent must obey both because a flip tests only the groups
it flips. The split, the dead-end move, and `_expand_change_allo_state` filter their
children; the other moves preserve both rules.
"""
function _assert_emission_rules(m::Union{Mechanism, AllostericMechanism})
    label(g) = join((join(_forward_sides(s), " → ") for s in steps(m)[g]), ", ")
    flux = _flux_carrying_groups(steps(m), reaction(m))
    for (g, group) in enumerate(steps(m))
        is_equilibrium(first(group)) || flux[g] || error(
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
    _assert_atom_conserving(m::Mechanism)
    _assert_atom_conserving(am::AllostericMechanism)

Assert every step of an enumerated mechanism conserves atoms (see
`_assert_step_atom_conserving`). This is the enumeration-path guardrail
against atom-non-conserving covalent intermediates; it is intentionally
NOT a `Step` / `Mechanism` constructor check, since hand-written
`@enzyme_mechanism` fixtures use placeholder atoms and folded steps.
"""
function _assert_atom_conserving(m::Union{Mechanism, AllostericMechanism})
    for group in steps(m), s in group
        _assert_step_atom_conserving(reaction(m), s)
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
    _catalytic_topologies(reaction) -> Vector{Vector{Step}}

Build catalytic cycle topologies by constructive backtracking.
Each topology is a set of steps forming one or more complete
catalytic cycles (E -> ... -> E). Errors when the reaction has none.
"""
function _catalytic_topologies(reaction::EnzymeReaction)
    sub_names = Symbol[name(s) for s in substrates(reaction)]
    prod_names = Symbol[name(p) for p in products(reaction)]

    # The form on :E binding `on_subs` and `on_prods` after a route consumed `consumed`
    # and released `released`. Its covalent residual adds the consumed substrates not
    # bound and subtracts the released and bound products; it is `Residual()` exactly
    # when their atoms cancel (no covalent residue remains).
    function form(consumed, released, on_subs, on_prods)
        added, subtracted = setdiff(consumed, on_subs), [released; on_prods]
        res = isempty(_net_atoms(reaction, added, subtracted)) ? Residual() :
              Residual(Substrate.(added), Product.(subtracted))
        Species(Metabolite[Substrate.(on_subs); Product.(on_prods)], :E, res)
    end

    # Collect all complete catalytic paths as Step lists. A route starts at free E, binds
    # substrates, isomerizes a substrate-bound form to a product-bound one, releases those
    # products one at a time in every order, and repeats until it has consumed every
    # substrate and released every product. `cur` is the current form, `consumed` and
    # `released` are the route's history, and `on_subs` and `on_prods` are what `cur`
    # binds; the atoms on the enzyme are those of `consumed` minus those of `released`.
    paths = Vector{Step}[]
    function walk!(cur, consumed, released, on_subs, on_prods, path)
        if isempty(bound(cur)) && length(consumed) == length(sub_names) &&
                length(released) == length(prod_names)
            push!(paths, copy(path))
            return
        end
        # Steps from `cur` to the form the arguments describe, then walks on from it.
        function step!(takes_up, gives_off, con, rel, on_s, on_p)
            to = form(con, rel, on_s, on_p)
            push!(path, Step(cur, to, takes_up, gives_off, true))
            walk!(to, con, rel, on_s, on_p, path)
            pop!(path)
        end
        if !isempty(on_prods)
            # Release a bound product; nothing binds until every bound product is gone.
            for p in on_prods
                step!(Metabolite[], Metabolite[Product(p)], consumed, [released; p],
                      Symbol[], filter(!=(p), on_prods))
            end
            return
        end
        remaining_subs = setdiff(sub_names, consumed)
        for s in remaining_subs
            step!(Metabolite[Substrate(s)], Metabolite[], [consumed; s], released,
                  [on_subs; s], Symbol[])
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
            step!(Metabolite[], Metabolite[], consumed, released, Symbol[], subset)
        end
    end
    # Start from free enzyme (:E, no bound metabolites).
    walk!(Species(Metabolite[], :E), Symbol[], Symbol[], Symbol[], Symbol[], Step[])

    # --- Group paths by isomerization pattern ---
    iso_groups = Dict{Set{Step}, Vector{Vector{Step}}}()
    for path in paths
        push!(get!(iso_groups, Set(filter(is_iso, path)), Vector{Step}[]), path)
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
    # Iterate iso_groups deterministically (smaller iso-step counts first,
    # then by sorted iso-step names) so topology output order is stable —
    # `Set{Step}` hashing is not value-stable.
    level_of(wo) = Dict{Symbol,Int}(m => i for (i, level) in enumerate(wo) for m in level)
    sub_levels = level_of.(_weak_orderings(sub_names))
    prod_levels = level_of.(_weak_orderings(prod_names))
    route_order(path, ::Type{T}) where {T} =
        Symbol[name(bound_metabolite(s)::T) for s in path if bound_metabolite(s) isa T]
    sorted_iso_pats = sort!(collect(keys(iso_groups));
        by = pat -> (length(pat),
                     sort([(string(name(from_species(s))), string(name(to_species(s))))
                           for s in pat])))
    result = Vector{Step}[]
    for iso_pat in sorted_iso_pats
        group = [(path, route_order(path, Substrate), route_order(path, Product))
                 for path in iso_groups[iso_pat]]
        seen_topos = Set{Set{Step}}()
        for sl in sub_levels, pl in prod_levels
            topo_keys = Set{Step}()
            for (path, binds, releases) in group
                issorted(binds; by = m -> sl[m]) && issorted(releases; by = m -> pl[m]) &&
                    union!(topo_keys, path)
            end
            (isempty(topo_keys) || topo_keys ∈ seen_topos) && continue
            push!(seen_topos, topo_keys)

            steps = sort!(collect(topo_keys); by = s -> (is_iso(s),
                string(name(from_species(s))), string(name(to_species(s)))))

            # The first iso step is the (single) SS step; every other step is
            # RE. Rebuild each Step with that tag (Step is immutable; direction
            # is unaffected by is_equilibrium).
            iso_idx = findfirst(is_iso, steps)
            push!(result, Step[_with_equilibrium(s, i != iso_idx)
                               for (i, s) in enumerate(steps)])
        end
    end
    isempty(result) && error("no catalytic cycle for $(reaction): either a product " *
                             "shares a substrate's name, or every route needs an " *
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
    _expand_substrate_product_dead_ends(topos, reaction) -> Vector{Vector{Step}}

For each catalytic topology, enumerate substrate/product dead-end
form combinations and return each resulting mechanism as a flat
`Vector{Step}`. A dead-end form is a catalytic form with one more
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
    topos::Vector{Vector{Step}},
    reaction::EnzymeReaction,
)
    sub_names = Set(name(s) for s in substrates(reaction))
    prod_names = Set(name(p) for p in products(reaction))
    mets = sort!(collect(union(sub_names, prod_names)))

    # A declared shared catalytic site forbids its (substrate, product) pair
    # from co-occupying the catalytic site, so keep only competition patterns
    # whose forbidden-edge set contains every declared pair. The complete
    # bipartite pattern contains all edges and always survives, so the list is
    # never empty.
    patterns = filter(pat -> all(in(pat), shared_catalytic_site(reaction)),
                      _competition_patterns(sub_names, prod_names))

    role(m::Symbol) = m in sub_names ? Substrate(m) : Product(m)
    held(sp::Species) = Set(name(b) for b in bound(sp))
    result = Vector{Step}[]
    for topo in topos
        catalytic = _forms([topo])
        dead_ends = Species[]
        for sp in catalytic, m in mets
            b = held(sp)
            m in b && continue
            push!(b, m)
            subs, prods = intersect(b, sub_names), intersect(b, prod_names)
            (isempty(subs) || isempty(prods) || subs == sub_names || prods == prod_names) &&
                continue
            push!(dead_ends, Species(Metabolite[bound(sp)..., role(m)], conformation(sp),
                                     residual(sp)))
        end
        unique!(dead_ends)
        edges = Set((from_species(s), to_species(s)) for s in topo)
        seen = Set{Vector{Species}}()
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
            for (i, sp1) in enumerate(forms), (j, sp2) in enumerate(forms)
                conformation(sp1) == conformation(sp2) && residual(sp1) == residual(sp2) &&
                    length(hs[j]) == length(hs[i]) + 1 && issubset(hs[i], hs[j]) &&
                    !((sp1, sp2) in edges) || continue
                met = only(setdiff(hs[j], hs[i]))
                push!(steps, Step(sp1, sp2, Metabolite[role(met)], Metabolite[], true))
            end
            push!(result, steps)
        end
    end
    result
end

# ─── Compilation ─────────────────────────────────────────────

"""
    compile_mechanism(m::Mechanism)
    compile_mechanism(am::AllostericMechanism)

Convert a `Mechanism` to an `EnzymeMechanism`, or an
`AllostericMechanism` to an `AllostericEnzymeMechanism`.
"""
compile_mechanism(m::Mechanism) = EnzymeMechanism(m)
compile_mechanism(am::AllostericMechanism) = AllostericEnzymeMechanism(am)

# ─── Mechanism Enumeration ───────────────────────────────────

"""
    _seed_groups(steps) -> Vector{Vector{Step}}

The kinetic groups of a seed's flat step list. Steps that take up and give off the
same metabolites (`_step_kind`) with the same RE/SS flag share one group; every
isomerization is a group of its own. The `Mechanism` constructor canonicalizes group
and step order, so the order returned here does not matter.
"""
function _seed_groups(steps::Vector{Step})
    groups = Dict{Any, Vector{Step}}()
    for (i, s) in enumerate(steps)
        push!(get!(groups, is_iso(s) ? i : (_step_kind(s), is_equilibrium(s)), Step[]), s)
    end
    collect(values(groups))
end


# ─── Expansion Moves ─────────────────────────────────────────

"""
    _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})

RE→SS expansion move. A flip unit is a whole kinetic group that is all-RE, binds
no regulator (competitive-inhibitor binding stays at rapid equilibrium by
modeling choice), and holds a flux-carrying step (`_flux_carrying_groups` on the
all-steady-state graph, `_all_steady_state`; a group with none exposes only
equilibrium ratios under every assignment and would gain a phantom parameter).
A flank of a qualifying chain whose isomerization is steady state is no unit
(`_chain_flank_groups`). One child is produced per minimal set of units whose joint flip
raises the RE segment count (`_minimal_flips`): a single group when it
cuts a segment on its own, several groups when each alone is bridged by an RE
route through the others — as happens once a split has separated a
metabolite's binding steps, or a catalytic step from its inhibitor-bound
mirror. A flip that leaves the segment count unchanged adds an SS step whose
endpoints share a segment, which the rate equation never sees. A flip set
that leaves a rapid-equilibrium segment with no bottom form
(`_bottomless_re_segment`, which the constructor rejects) also counts as
failed, so the minimal-set search extends it instead of emitting it; so does
a set one of whose flipped groups carries no flux in the child
(`_flux_carrying_groups` on the flipped groups): the child is its parent plus
one phantom per such group, and its flux-carrying supersets are the smallest
identifiable relaxations. The non-degenerate mechanisms beyond these stay reachable —
through the other cut orders of the same ring, or through the extended set when
no other order gains. All other groups, the reaction, and (for allosteric) the
catalytic-allo tags, multiplicity, and regulatory sites are preserved
verbatim.

A parent that `_requires_hyperbolic_catalysis` (a conformational mechanism over
more than one catalytic subunit) keeps only the children whose catalytic scheme
passes `_hyperbolic_catalysis`. The check runs on the emitted minimal sets, not
inside the gain predicate: a failing set would otherwise be extended with more
flips, and more steady-state steps never restore a hyperbolic equation, so its
supersets need no visit.
"""
function _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})
    flanks = _chain_flank_groups(m)
    eligible(g) = !(g in flanks) &&
        !any(s -> any(x -> x isa Regulator, consumed(s)) ||
                  any(x -> x isa Regulator, released(s)), steps(m)[g])
    base = _re_segment_count(steps(m))
    flips = _minimal_flips(gs -> _re_segment_count(gs) > base, steps(m), reaction(m),
                           eligible)
    children = typeof(m)[_with(m; groups = gs) for gs in flips]
    _requires_hyperbolic_catalysis(m) ? filter(_hyperbolic_catalysis, children) : children
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
function _chain_flank_groups(m::Mechanism)
    groups = steps(m)
    lone = Dict(only(group) => g for (g, group) in enumerate(groups) if length(group) == 1)
    at = Dict{Species, Vector{Step}}()
    for group in groups, s in group, sp in (from_species(s), to_species(s))
        push!(get!(at, sp, Step[]), s)
    end
    flank(s0, x) = begin
        rest = filter(!=(s0), at[x])
        length(rest) == 1 || return nothing
        s = only(rest)
        is_binding(s) && to_species(s) == x ? get(lone, s, nothing) : nothing
    end
    out = Set{Int}()
    for s0 in keys(lone)
        is_iso(s0) && !is_equilibrium(s0) || continue
        f1, f2 = flank(s0, from_species(s0)), flank(s0, to_species(s0))
        f1 === nothing || f2 === nothing || union!(out, (f1, f2))
    end
    out
end
_chain_flank_groups(::AllostericMechanism) = Set{Int}()

"""`s` with its rapid-equilibrium flag set to `flag`; `s` itself when it has that flag, so
mechanisms built from one another share their unchanged steps."""
_with_equilibrium(s::Step, flag::Bool) = is_equilibrium(s) == flag ? s :
    Step(from_species(s), to_species(s), consumed(s), released(s), flag)

"""Every step of `groups` rebuilt at steady state. The flip move pre-filters its
units on this graph: a step that carries no flux with every step steady-state
carries none under any assignment (flips only refine segments), so a group with
no such step can never flip usefully."""
_all_steady_state(groups::Vector{Vector{Step}}) =
    [_with_equilibrium.(group, false) for group in groups]

"""
    _minimal_flips(admissible, base, rxn, eligible = _ -> true)

The groups of `base` after every inclusion-minimal flip to steady state that passes. A
unit is a rapid-equilibrium group of `base` that carries flux on the all-steady-state
graph (`_all_steady_state`, `_flux_carrying_groups`) and passes `eligible`. A candidate
takes each flipped group from that graph, so the candidates share their steady-state
steps, and every other group from `base` by reference. It passes when `admissible` holds
on its groups, it leaves no rapid-equilibrium segment without a bottom form
(`_bottomless_re_segment`) and every flipped group carries flux in it. The sets come from
`_minimal_gaining_sets`, so a failing set is extended and a passing one is not.
"""
function _minimal_flips(admissible, base::Vector{Vector{Step}}, rxn::EnzymeReaction,
                        eligible = _ -> true)
    steady = _all_steady_state(base)
    flux = _flux_carrying_groups(steady, rxn)
    units = [g for g in eachindex(base)
             if all(is_equilibrium, base[g]) && flux[g] && eligible(g)]
    flipped(sel) = (gs = copy(base); gs[units[sel]] = steady[units[sel]]; gs)
    passes(sel) = (gs = flipped(sel);
                   admissible(gs) && _bottomless_re_segment(gs) === nothing &&
                   all(_flux_carrying_groups(gs, rxn)[units[sel]]))
    [flipped(sel) for sel in _minimal_gaining_sets(passes, _ -> 1:length(units))]
end

"""
    _flux_carrying_steps(groups, rxn) -> Vector{BitVector}

One flag per step of `groups`, parallel to `groups`: whether the step carries net
steady-state flux for generic parameter values. Rapid-equilibrium steps hold
equilibrium mass and are flagged `false`; the rule that uses these flags concerns
steady-state groups only.

Contract each rapid-equilibrium segment to one vertex (`_re_segment_extras`) and
make every steady-state step an edge between its two forms' segments, weighted by
the substrates it takes up minus the products, plus its `from` form's offsets minus
its `to` form's. Around any cycle the offsets telescope and the uptakes sum to the
net turnover times the reactant count, so a cycle has weight zero exactly when it
runs no net reaction. A self-loop (both forms in one segment) carries flux iff its
weight is nonzero. Any other edge carries flux iff its biconnected block
(`_unbalanced_blocks`) holds a cycle of nonzero weight: every edge of such a block lies
on one (join the edge to the cycle by two disjoint paths; one of the two resulting
cycles is unbalanced), while in a balanced block detailed balance holds along every
cycle and each step's flux vanishes. The zero-flux verdict holds for any parameters
and grouping; the flux-carrying verdict for one-way steps, which every step the
enumerator emits is. The test reads each step's metabolite lists, so it holds for
fused and Theorell–Chance steps as well.
"""
function _flux_carrying_steps(groups::Vector{Vector{Step}}, rxn::EnzymeReaction)
    species, segments, extras, idx, seg = _re_segment_extras(groups)
    rho = _reactant_signs(rxn)
    offset(i) = sum(get(rho, x, 0) * e for (x, e) in extras[i]; init = 0)
    weight(s) = _uptake_weight(s, rho) +
                offset(idx[from_species(s)]) - offset(idx[to_species(s)])
    flags = [falses(length(group)) for group in groups]
    edges = Tuple{Int, Int}[]; weights = Int[]; owner = Tuple{Int, Int}[]
    for (g, group) in enumerate(groups), (j, s) in enumerate(group)
        is_equilibrium(s) && continue
        a, b = seg[idx[from_species(s)]], seg[idx[to_species(s)]]
        if a == b
            flags[g][j] = weight(s) != 0
        else
            push!(edges, (a, b)); push!(weights, weight(s)); push!(owner, (g, j))
        end
    end
    block, unbalanced = _unbalanced_blocks(length(segments), edges, weights)
    for (e, (g, j)) in enumerate(owner)
        flags[g][j] = block[e] in unbalanced
    end
    flags
end

"""Net uptake sign of each reactant name of `rxn`: +1 for a substrate, −1 for a product."""
function _reactant_signs(rxn::EnzymeReaction)
    rho = Dict{Symbol, Int}()
    for s in substrates(rxn); rho[name(s)] = get(rho, name(s), 0) + 1; end
    for p in products(rxn);   rho[name(p)] = get(rho, name(p), 0) - 1; end
    rho
end

"""The signs `rho` (`_reactant_signs`) summed over the metabolites `s` takes up, minus
those summed over the metabolites it gives off."""
_uptake_weight(s::Step, rho::Dict{Symbol, Int}) =
    sum(get(rho, name(x), 0) for x in consumed(s); init = 0) -
    sum(get(rho, name(x), 0) for x in released(s); init = 0)

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
    _flux_carrying_groups(groups, rxn) -> BitVector

One flag per kinetic group: whether some step of the group carries net
steady-state flux (`_flux_carrying_steps`). A steady-state group with no such step
exposes only the ratio of its two constants, so the moves never emit one; a
zero-flux step inside a group that also holds a flux-carrying step costs nothing,
because the group's shared constants are pinned by the step that carries flux.
"""
_flux_carrying_groups(groups::Vector{Vector{Step}}, rxn::EnzymeReaction) =
    BitVector([any(f) for f in _flux_carrying_steps(groups, rxn)])

"""
`groups` with the isomerization `s0`, alone in its group, merged onto its product side:
`s0` and its group are removed, and every other step at `from_species(s0)` moves onto
`to_species(s0)`. The last substrate's binding into the substrate side becomes a fused
binding into the product side; the releases stay plain.
"""
function _merge_isomerization(groups::Vector{Vector{Step}}, s0::Step)
    x1, x2 = from_species(s0), to_species(s0)
    move(sp) = sp == x1 ? x2 : sp
    moved(s) = x1 in (from_species(s), to_species(s)) ?
        Step(move(from_species(s)), move(to_species(s)), consumed(s), released(s),
             is_equilibrium(s)) : s
    [Step[moved(s) for s in group] for group in groups if group != [s0]]
end

"""
`groups` with the form `x` eliminated, or `nothing` when that is not possible: `x` must
have exactly two steps, both bindings into it, Y + L → x and Z + R → x with Y ≠ Z. They
become one Theorell–Chance step Y + L → Z + R, rapid equilibrium, in a group of its own,
where Y + L is the binding of a substrate when one of the two is; their groups keep their
other steps.
"""
function _eliminate_form(groups::Vector{Vector{Step}}, x::Species)
    at_x = [s for group in groups for s in group if x in (from_species(s), to_species(s))]
    length(at_x) == 2 && all(s -> is_binding(s) && to_species(s) == x, at_x) ||
        return nothing
    entry, release = sort(at_x; by = s -> !(bound_metabolite(s) isa Substrate))
    from_species(entry) == from_species(release) && return nothing
    fused = Step(from_species(entry), from_species(release), consumed(entry),
                 consumed(release), true)
    kept = [filter(s -> !(s in at_x), group) for group in groups]
    push!(filter!(!isempty, kept), [fused])
end

"""Whether the steps `k` with `on[k]` close a cycle of nonzero weight, step `k` joining
forms `from[k]` and `to[k]` with weight `weight[k]`: giving every form a potential that
rises by each step's weight along it then fails somewhere, which a weighted union-find in
the buffers `parent` and `offset` (one entry per form) detects."""
function _has_unbalanced_cycle(from, to, weight, on, parent, offset)
    parent .= eachindex(parent); fill!(offset, 0)
    function root(x)
        p = 0
        while parent[x] != x
            p += offset[x]; x = parent[x]
        end
        x, p
    end
    for k in eachindex(from)
        on[k] || continue
        ra, pa = root(from[k]); rb, pb = root(to[k])
        if ra == rb
            pb == pa + weight[k] || return true
        else
            parent[rb] = ra; offset[rb] = pa + weight[k] - pb
        end
    end
    false
end

"""Whether a cycle of rapid-equilibrium steps of `groups` runs net turnover, which makes the
rate infinite. The steps for which `fast` holds count as rapid equilibrium, by default
those flagged so. Each RE step joins its two forms with its uptake of the reaction's
substrates minus its products as weight (`_uptake_weight`); a turnover cycle has weight
n·Σρ², never zero, so the RE steps hold one iff they close a cycle of nonzero weight
(`_has_unbalanced_cycle`)."""
function _re_turnover_cycle(groups::Vector{Vector{Step}}, rxn::EnzymeReaction,
                            fast = is_equilibrium)
    rho = _reactant_signs(rxn)
    re = [s for group in groups for s in group if fast(s)]
    idx = Dict{Species, Int}()
    vertex(sp) = get!(idx, sp, length(idx) + 1)
    ends = [(vertex(from_species(s)), vertex(to_species(s))) for s in re]
    n = length(idx)
    _has_unbalanced_cycle(first.(ends), last.(ends), [_uptake_weight(s, rho) for s in re],
                          trues(length(re)), zeros(Int, n), zeros(Int, n))
end

"""Whether the rate has a maximum as the metabolites of `side` (`Substrate` or `Product`)
grow: turning rapid equilibrium every steady-state step that takes up or gives off one of
them leaves no rapid-equilibrium turnover cycle (condition V)."""
function _has_vmax(groups::Vector{Vector{Step}}, rxn::EnzymeReaction, side::Type)
    names = Set(name(x) for x in (side === Substrate ? substrates(rxn) : products(rxn)))
    touches(s) = _any_named(consumed(s), names) || _any_named(released(s), names)
    !_re_turnover_cycle(groups, rxn, s -> is_equilibrium(s) || touches(s))
end

"""Whether some metabolite of `ms` has its name in `names`."""
_any_named(ms, names) = any(m -> name(m) in names, ms)

"""Whether `s` is a fused binding of a substrate: a chemistry step (`_is_chemistry`) that
binds a metabolite named in `subs`. The form it enters is a merged complex."""
_fused_substrate_binding(s::Step, subs) =
    (m = bound_metabolite(s); m !== nothing && _is_chemistry(s) && name(m) in subs)

"""Whether `s` takes up a reactant of one side of the reaction (names `subs` or `prods`)
and gives off one of the other, as a Theorell–Chance step does."""
_crosses_sides(s::Step, subs, prods) =
    _any_named(consumed(s), subs) && _any_named(released(s), prods) ||
    _any_named(consumed(s), prods) && _any_named(released(s), subs)

"""
Whether some chemistry node of `groups` is left both by a rapid-equilibrium step that
releases a substrate and by one that releases a product (condition C). A chemistry
node is a merged complex (a form entered by a fused binding of a substrate) with every form
joined to it by RE isomerizations, a set of forms joined to each other by RE isomerizations,
or an RE Theorell–Chance step. An RE step leaves a node through a node form F and releases M
when traversing it away from F gives M off: F is its `to_species` and M is consumed, or F is
its `from_species` and M is released. An RE Theorell–Chance node leaves through its own step
in both directions. Then the chemistry sits in rapid equilibrium with both sides and some
reactant is not needed at zero products or substrates.
"""
function _chemistry_equilibrates_both_sides(groups::Vector{Vector{Step}},
                                            rxn::EnzymeReaction)
    subs = Set(name(s) for s in substrates(rxn))
    prods = Set(name(p) for p in products(rxn))
    re = [s for group in groups for s in group if is_equilibrium(s)]
    releases(n, names) = any(re) do t
        to_species(t) in n && !(from_species(t) in n) && _any_named(consumed(t), names) ||
            from_species(t) in n && !(to_species(t) in n) && _any_named(released(t), names)
    end
    # The forms joined by RE isomerizations, one node per segment; a merged complex
    # outside every RE isomerization is a node of its own.
    species, segments, _, idx, _ = _re_segment_extras([filter(is_iso, re)])
    merged = [to_species(s) for group in groups for s in group
              if _fused_substrate_binding(s, subs)]
    nodes = [[Set(species[k]) for k in segments];
             [Set([x]) for x in merged if !haskey(idx, x)]]
    any(n -> releases(n, subs) && releases(n, prods), nodes) ||
        any(t -> _crosses_sides(t, subs, prods), re)
end

"""Whether `m` is degenerate by the structural conditions V and C: no maximal rate in one
direction (`_has_vmax`), or chemistry in rapid equilibrium with both sides
(`_chemistry_equilibrates_both_sides`). Read on `steps(m)`, an allosteric mechanism's
active-state graph. The beam expands a degenerate seed without fitting it (`_base_tier`)."""
function _degenerate(m::Union{Mechanism, AllostericMechanism})
    groups, rxn = steps(m), reaction(m)
    !(_has_vmax(groups, rxn, Substrate) && _has_vmax(groups, rxn, Product)) ||
        _chemistry_equilibrates_both_sides(groups, rxn)
end

"""
The tests of `_seed_variants` that build no `Step`, for one base `groups` of a seed: a
function of a candidate's groups (those of `groups`, some flipped to steady state) that
reads each step's flag from them and is true when the candidate has no rapid-equilibrium
turnover cycle (`_re_turnover_cycle`), a maximal rate both ways (`_has_vmax`) and no
chemistry in equilibrium with both sides (`_chemistry_equilibrates_both_sides`). The forms
and each step's ends, uptake weight and reactants are indexed once per base, so a
candidate costs a few passes over arrays. A base holds no isomerization, so each chemistry
node is a merged complex or a Theorell–Chance step. A set of steps holds a turnover cycle
iff it closes a cycle of nonzero weight (`_has_unbalanced_cycle`), which reuses the
screen's buffers. Each maximal-rate test runs on the candidate's rapid-equilibrium steps
and more, and adding steps never removes such a cycle, so a turnover cycle of the
candidate fails both tests.
"""
function _seed_candidate_screen(groups::Vector{Vector{Step}}, rxn::EnzymeReaction)
    flat = [s for group in groups for s in group]
    any(is_iso, flat) && error("_seed_candidate_screen: a seed base holds no isomerization")
    rho = _reactant_signs(rxn)
    subs = Set(name(x) for x in substrates(rxn))
    prods = Set(name(x) for x in products(rxn))
    index = Dict{Species, Int}()
    for s in flat, sp in (from_species(s), to_species(s))
        get!(index, sp, length(index) + 1)
    end
    from = [index[from_species(s)] for s in flat]
    to = [index[to_species(s)] for s in flat]
    weight = [_uptake_weight(s, rho) for s in flat]
    touches(names) = BitVector([_any_named(consumed(s), names) ||
                                _any_named(released(s), names) for s in flat])
    sub_step, prod_step = touches(subs), touches(prods)
    crossing = BitVector([_crosses_sides(s, subs, prods) for s in flat])
    complexes = unique(to[k] for k in eachindex(flat)
                       if _fused_substrate_binding(flat[k], subs))
    leaves(x, names) = BitVector([to[k] == x && _any_named(consumed(flat[k]), names) ||
                                  from[k] == x && _any_named(released(flat[k]), names)
                                  for k in eachindex(flat)])
    exits = [(leaves(x, subs), leaves(x, prods)) for x in complexes]
    parent = zeros(Int, length(index)); offset = zeros(Int, length(index))
    turnover(on) = _has_unbalanced_cycle(from, to, weight, on, parent, offset)
    meets(a::BitVector, b::BitVector) = any(k -> a[k] && b[k], eachindex(a))
    re = falses(length(flat)); kept = falses(length(flat))
    function screen(gs::Vector{Vector{Step}})
        k = 0
        for g in gs, s in g
            re[k += 1] = is_equilibrium(s)
        end
        kept .= re .| sub_step
        turnover(kept) && return false
        kept .= re .| prod_step
        turnover(kept) && return false
        meets(re, crossing) && return false
        !any(((es, ep),) -> meets(re, es) && meets(re, ep), exits)
    end
end

"""
The merged and Theorell–Chance variants of the seed `m`. Merging every isomerization of `m`
onto its product side (`_merge_isomerization`) and setting every step at rapid equilibrium
gives the merged base; eliminating one merged complex with two steps (`_eliminate_form`)
gives a Theorell–Chance base. MERGE and ELIM conserve atoms, which each base asserts. For
each base, every inclusion-minimal set of its groups (`_minimal_flips`) whose flip to
steady state gives a valid, flux-carrying, non-degenerate candidate is a variant: no
rapid-equilibrium turnover cycle (`_re_turnover_cycle`), a maximal rate both ways
(`_has_vmax`), no chemistry in equilibrium with both sides
(`_chemistry_equilibrates_both_sides`), no bottomless segment, and every steady-state group
carrying flux. The first three tests read index arrays (`_seed_candidate_screen`); the
other two run only on a candidate that passes them. A merged variant whose every merged
complex has both its steps steady state, each alone in its group, is skipped: the
unmerged form with rapid-equilibrium flanks has its family at the same count.
"""
function _seed_variants(m::Mechanism)
    rxn = reaction(m)
    isos = [s for group in steps(m) for s in group if is_iso(s)]
    merged = foldl(_merge_isomerization, isos; init = steps(m))
    merged = [_with_equilibrium.(group, true) for group in merged]
    complexes = [to_species(s) for s in isos]
    bases = Vector{Vector{Step}}[[merged];
        filter(!isnothing, [_eliminate_form(merged, x) for x in complexes])]
    lumping_twin(gs) = all(complexes) do x
        at_x = [(s, group) for group in gs for s in group
                if x in (from_species(s), to_species(s))]
        length(at_x) == 2 &&
            all(((s, group),) -> !is_equilibrium(s) && length(group) == 1, at_x)
    end
    variants = Mechanism[]
    for (i, base) in enumerate(bases)
        for group in base, s in group
            _assert_step_atom_conserving(rxn, s)
        end
        for gs in _minimal_flips(_seed_candidate_screen(base, rxn), base, rxn)
            i == 1 && lumping_twin(gs) && continue
            push!(variants, Mechanism(rxn, gs))
        end
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
_requires_hyperbolic_catalysis(::Mechanism) = false
_requires_hyperbolic_catalysis(am::AllostericMechanism) = catalytic_multiplicity(am) > 1

"""
    _hyperbolic_catalysis(m) -> Bool

Whether the King–Altman denominator of `m`'s catalytic scheme has degree at most
1 in every substrate and product concentration. Every binding of a substrate or
product at its catalytic site counts, abortive complexes included. Steps touching
a form that carries a declared inhibitor are left out: an inhibitor binds a site
of its own by definition, so the powers its binding adds, including those of a
substrate declared as a dead-end inhibitor, are a separate source that
conformational mechanisms keep.

The equation is a sum over rapid-equilibrium (RE) segments and spanning
arborescences of the segment graph toward each segment. A denominator term is
the root segment's weight times the weight of every tree edge, so the exponent
of `X` in a term is the most `X` any form of the root segment carries beyond
the segment's bottom form, plus, per tree edge, one if the step binds `X` in
the tree direction and the count of `X` the edge's source form carries beyond
its bottom (`_re_segment_extras`). The degree exceeds 1 exactly when one
segment or one edge scores 2 or more, or the root scores 1 and some
arborescence toward it holds a scoring edge, or some arborescence holds two
scoring edges. An arborescence toward `S` containing given edges exists iff
every segment still reaches `S` once each given edge's source keeps that edge
as its only way out (`_all_reach`).
"""
function _hyperbolic_catalysis(m::Union{Mechanism, AllostericMechanism})
    on_catalytic_site(s) = !any(b -> b isa Regulator,
                                vcat(bound(from_species(s)), bound(to_species(s))))
    groups = [filter(on_catalytic_site, group) for group in steps(m)]
    species, segments, extras, idx, seg_of = _re_segment_extras(groups)
    # Directed segment-graph edges: source segment, target segment, source form,
    # metabolites bound in that direction.
    edges = Tuple{Int, Int, Int, Vector{Symbol}}[]
    for group in groups, s in group
        is_equilibrium(s) && continue
        a, b = idx[from_species(s)], idx[to_species(s)]
        seg_of[a] == seg_of[b] && continue
        push!(edges, (seg_of[a], seg_of[b], a, name.(consumed(s))))
        push!(edges, (seg_of[b], seg_of[a], b, name.(released(s))))
    end
    rxn = reaction(m)
    mets = [name(x) for side in (substrates(rxn), products(rxn)) for x in side]
    n = length(segments)
    for x in mets
        score(e) = count(==(x), e[4]) + get(extras[e[3]], x, 0)
        carrying = [e for e in edges if score(e) > 0]
        any(e -> score(e) > 1, carrying) && return false
        root_score(k) = maximum((get(extras[i], x, 0) for i in segments[k]); init=0)
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

"""Number of rapid-equilibrium segments (connected components of the RE subgraph) of
`groups`, by a union-find over the forms. The flip move measures each candidate on its
groups, so a candidate the constructors would reject is never constructed."""
function _re_segment_count(groups::Vector{Vector{Step}})
    index = Dict{Species, Int}()
    parent = Int[]
    vertex(sp) = get!(() -> (push!(parent, length(parent) + 1); length(parent)), index, sp)
    find(x) = (while parent[x] != x; x = parent[x] = parent[parent[x]]; end; x)
    for group in groups, s in group
        a, b = vertex(from_species(s)), vertex(to_species(s))
        is_equilibrium(s) && (parent[find(a)] = find(b))
    end
    count(i -> parent[i] == i, eachindex(parent))
end

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
Wegscheider cycle ties straight back is not a model: the derivation substitutes
the constant away and the child's equation is the parent's. Such a bipartition
is kept as a seed for pairs and larger sets, because the cycle that absorbs it
is broken by also splitting another group on that cycle — random-order binding
needs one split per substrate before any constant frees up. Partners for a
rejected set are the all-RE groups sharing a rapid-equilibrium segment with any
step of the split groups: a tie needs an RE cycle through the carve, and every
group on that cycle lies in that segment. The count is evaluated without building the child
for a `Mechanism` (`_partition_independent_count`, cycle basis once per parent)
and by the combined state solve on the built child for an `AllostericMechanism`.
A part of a steady-state group none of whose steps carries flux
(`_flux_carrying_steps`, computed once on the parent, since a split moves no edge)
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
The reaction and (for allosteric) multiplicity and regulatory sites are
preserved; both parts of a split group inherit its catalytic allo-state tag.
"""
function _expand_split_kinetic_group(m::Union{Mechanism, AllostericMechanism})
    groups = steps(m)
    flux = _flux_carrying_steps(groups, reaction(m))
    units = Tuple{Int, Tuple{Vector{Step}, Vector{Step}}}[]
    reverted = Bool[]
    for g in eachindex(groups), bp in _context_bipartitions(groups[g])
        parts = _revert_zero_flux_parts(groups[g], bp, flux[g])
        push!(units, (g, parts)); push!(reverted, parts != bp)
    end
    isempty(units) && return typeof(m)[]
    gain = _split_gain_test(m, units)
    gains(sel) =
        (!any(u -> reverted[u], sel) ||
         _bottomless_re_segment(_bipartitioned_groups(groups, units[sel])[1]) ===
         nothing) && gain(sel)
    # RE segment ids touched by each kinetic group's steps; empty for a group holding an
    # SS step, which lies on no RE cycle and is never a split partner.
    _, _, _, idx, seg = _re_segment_extras(groups)
    segments = [is_equilibrium(first(group)) ?
                Set{Int}(seg[idx[from_species(s)]] for s in group) : Set{Int}()
                for group in groups]
    partners(sel) = begin
        isempty(sel) && return 1:length(units)
        used = Set(units[u][1] for u in sel)
        touched = reduce(union, (segments[units[u][1]] for u in sel))
        [u for u in 1:length(units)
         if !(units[u][1] in used) &&
            !isempty(intersect(segments[units[u][1]], touched))]
    end
    sets = _minimal_gaining_sets(gains, partners)
    filter!(c -> isempty(_redundant_copy_groups(c)),
            typeof(m)[_apply_bipartitions(m, units[sel]) for sel in sets])
end

"""Gain test for the split move: `sel -> Bool`, true when applying the selected
units raises the independent-parameter count above the parent's. A unit's steps are
matched to the parent's flat steps by reaction (forms and metabolite lists), since a
reverted part carries the other flag; each step is counted under the kind it has in
the child."""
function _split_gain_test(m::Mechanism, units)
    counter = _partition_independent_count(m)
    flat = _flat_steps(m)
    reaction_key(s) = (from_species(s), to_species(s), consumed(s), released(s))
    position = Dict(reaction_key(s) => j for (j, (s, _)) in enumerate(flat))
    parent_ids = [g for (_, g) in flat]
    parent_kinds = [_count_kind(s) for (s, _) in flat]
    base = counter(parent_ids, parent_kinds)
    function gains(sel)
        ids = copy(parent_ids)
        kinds = copy(parent_kinds)
        for (i, u) in enumerate(sel), (part_number, part) in enumerate(units[u][2]),
            s in part
            j = position[reaction_key(s)]
            part_number == 2 && (ids[j] = length(steps(m)) + i)
            kinds[j] = _count_kind(s)
        end
        counter(ids, kinds) > base
    end
    gains
end

function _split_gain_test(am::AllostericMechanism, units)
    base = _independent_param_count(am)
    sel -> _independent_param_count(_apply_bipartitions(am, units[sel])) > base
end

"""
The endpoint of `s` that does not carry the step's free metabolites: the form
they bind to. `from_species(s)` when the step consumes a metabolite (every
binding, plain or fused, is stored with its metabolite consumed), `to_species(s)`
when it only releases metabolites (two or more: a step releasing one is stored
as the binding it reverses), and `from_species(s)` for an iso step (both lists
empty).
"""
_context_form(s::Step) =
    isempty(consumed(s)) && !isempty(released(s)) ? to_species(s) : from_species(s)

"""
    _context_bipartitions(group) -> Vector{Tuple{Vector{Step}, Vector{Step}}}

Ways to divide one kinetic group in two by binding context: another ligand
already bound, the enzyme's conformation, or its covalent residual. Encodes
"the affinity for this metabolite may depend on what else is bound, on the
conformation, or on the covalent state" — the ligand family includes a
competitive inhibitor, which separates a catalytic step from its inhibitor-bound
mirror. Each context value divides the steps whose context form
(`_context_form`) carries it from the rest; the group's own bound metabolite is
never a context, and a group whose forms share one conformation and one residual
gets no bipartition from those two families. Both parts are nonempty, the first
part holds the group's first step, context values that induce the same division
give one bipartition, and the order is ligands (by role then name), then
conformations (by name), then residuals, for deterministic output.
"""
function _context_bipartitions(group::Vector{Step})
    own = bound_metabolite(first(group))
    forms = [_context_form(s) for s in group]
    ligands = unique!(Metabolite[b for f in forms for b in bound(f) if b != own])
    division(carries) = Step[s for (s, f) in zip(group, forms) if carries(f)]
    divisions = [division(f -> y in bound(f))
                 for y in sort!(ligands; by = b -> (string(typeof(b)), string(name(b))))]
    append!(divisions, division(f -> conformation(f) == c)
            for c in sort!(unique(conformation(f) for f in forms); by = string))
    append!(divisions, division(f -> residual(f) == r)
            for r in sort!(unique(residual(f) for f in forms); by = string))
    parts(with) = (without = Step[s for s in group if !(s in with)];
                   first(group) in with ? (with, without) : (without, with))
    unique!(first, Tuple{Vector{Step}, Vector{Step}}[
        parts(with) for with in divisions if 0 < length(with) < length(group)])
end

"""
Replace each selected group by its two bipartition parts. `selection` pairs a
group index with one of that group's `_context_bipartitions`; each group appears
at most once. For an allosteric mechanism both parts inherit the group's
catalytic allo-state tag (splitting is a parameter-relaxation move that must not
change A/I semantics).
"""
function _apply_bipartitions(m::Union{Mechanism, AllostericMechanism}, selection)
    groups, origin = _bipartitioned_groups(steps(m), selection)
    m isa Mechanism ? _with(m; groups) :
        _with(m; groups, states = cat_allo_states(m)[origin])
end

"""
The bipartition `bp` of `group` with every part none of whose steps carries flux
(`flags`, one per step of `group`) rebuilt at rapid equilibrium. A rapid-equilibrium
group's steps are all flagged `false` (`_flux_carrying_steps`) and rebuild as themselves,
so its bipartition comes back equal to `bp`.
"""
function _revert_zero_flux_parts(group::Vector{Step}, bp, flags::BitVector)
    carries = Dict(zip(group, flags))
    map(part -> any(s -> carries[s], part) ? part : _with_equilibrium.(part, true), bp)
end

"""Groups of `groups` with each selected group replaced by its two parts, plus
the index of the original group each new group came from."""
function _bipartitioned_groups(groups::Vector{Vector{Step}}, selection)
    parts = Dict(selection)
    ([part for (g, group) in enumerate(groups) for part in get(parts, g, (group,))],
     [g for (g, group) in enumerate(groups) for _ in get(parts, g, (group,))])
end

"""
    _productive_twin(groups) -> (site, ligand) -> Union{Species, Nothing}

The productive form of `groups`, a form bound to no competitive inhibitor, whose weight
is proportional to that of the complex `ligand` forms at `site` by a rapid-equilibrium
step, or `nothing` when none is: a form of the site's rapid-equilibrium segment with the
complex's offsets (`_re_segment_extras`), the one of the complex's composition when
several are, else the first by name. Within a segment rapid equilibrium fixes each
weight relative to another as a constant times the metabolites bound along the route
between them, so equal offsets mean proportional weights, and steady-state data see only
their sum. The offsets reach a copy across a rapid-equilibrium isomerization (a
ping-pong second chemistry step, a conformational isomer), where no form shares the
complex's composition. A form of the same composition in another segment, formed by a
steady-state binding, is no twin: its weight is not a constant multiple of the
complex's, merging the two is no reparameterization, and the copy's constant may stay
visible.

Only productive forms are twin sources, so the complex itself never is one. A complex
that duplicates only a copy's complex shares its weight with that copy's constant, and
the sites that pin either copy may keep both constants separable, so such a match never
makes a twin.
"""
function _productive_twin(groups::Vector{Vector{Step}})
    species, segments, extras, idx, seg = _re_segment_extras(groups)
    productive(sp) = !any(b -> b isa CompetitiveInhibitor, bound(sp))
    # Sorted names of a form's bound metabolites with `extra` also bound (a copy counts as
    # the reactant it copies), with its conformation and residual.
    composition(sp::Species, extra::Metabolite...) =
        (sort!(vcat(Symbol[name(b) for b in bound(sp)],
                    Symbol[name(x) for x in extra])), conformation(sp), residual(sp))
    function (site::Species, ligand::Metabolite)
        i = idx[site]
        target = mergewith(+, extras[i], Dict(name(ligand) => 1))
        twins = [species[j] for j in segments[seg[i]]
                 if productive(species[j]) && extras[j] == target]
        isempty(twins) && return nothing
        wanted = composition(site, ligand)
        k = findfirst(t -> composition(t) == wanted, twins)
        k === nothing ? argmin(name, twins) : twins[k]
    end
end

"""
Whether `group`, the steps of one kinetic group, binds a competitive inhibitor, and only
at sites whose complex has a productive twin (`twin`, a `_productive_twin` of one state's
graph).
"""
function _all_twin(group::Vector{Step}, twin)
    ligand = bound_metabolite(first(group))
    ligand isa CompetitiveInhibitor &&
        all(s -> twin(from_species(s), ligand) !== nothing, group)
end

"""
    _gauge_rescaling(groups, g, twin, label) -> Union{Dict{Int, Tuple}, Nothing}

The dwell gauge of the competitive-inhibitor copy bound by kinetic group `g` of
`groups`, one conformational state's step graph: the rescaling each nonempty kinetic
group needs, or `nothing` when the gauge does not exist. It does not exist when group `g`
binds no competitive inhibitor or a complex it forms has no productive twin (`_all_twin`,
with `twin` a `_productive_twin` of `groups`), or when some group would need two
rescalings. The gauge gives every form a factor σ: ρ_T to each twin T, one factor s to
every complex of the copy and 1 to every other form. Multiplying a rapid-equilibrium
group's K by σ(from)/σ(to), and a steady-state group's forward constant by σ(from) and
its reverse by σ(to), divides each form's weight by its factor and leaves every flux
as it was. A twin's weight is a constant multiple of its complexes', so the factors
can keep each twin's total with its complexes, and the rate law, unchanged while K*
moves: the copy's constant is a phantom. A kinetic group shares its constants, so every
group, the copy's own included, must have one ratio σ(from)/σ(to) over its steps (rapid
equilibrium), or one σ(from) and one σ(to) (steady state). The copy's group, at 1 on
each site, holds all its complexes at the one factor s, and a complex never takes its
twin's factor: its weight moves opposite to the twin's. A rapid-equilibrium mirror step,
between two complexes, therefore has ratio 1, as its parent step between two sites
has, whatever the twins; a steady-state mirror needs s beside its parent step's 1, and
the gauge fails. When each site binds the copied ligand at rapid equilibrium, the gauge
is the merge of each complex into its twin, the twin's binding group taking K_h + K* in
place of K_h, and the family is that of the mechanism without the copy.

The ρ are generic numbers, equal only where the structure makes them equal: twins
formed from their sites by RE bindings of the copied ligand in one kinetic group h, each
the twin of one complex, share one factor (the merge gives it as ρ = K_h/(K_h + K*));
every other twin has its own. `label` names the classes: it
receives the group index of a shared class (h) or of the copy's complexes (g), or the
twin of a class of its own, and returns the class's symbol, so the caller decides which
classes are one number. A group's rescaling is the pair of its ends' symbols, `nothing`
standing for the factor 1, and `(nothing, nothing)` for a rapid-equilibrium group whose
two ends scale alike.
The test is bookkeeping over the steps, with no parameters and no numerics.
"""
function _gauge_rescaling(groups::Vector{Vector{Step}}, g::Int, twin, label)
    _all_twin(groups[g], twin) || return nothing
    ligand = bound_metabolite(first(groups[g]))::Metabolite
    twin_of = Dict(to_species(s) => twin(from_species(s), ligand) for s in groups[g])
    binding_group = Dict{Tuple{Species, Species}, Int}()
    for (h, group) in enumerate(groups), s in group
        bm = bound_metabolite(s)
        is_equilibrium(s) && bm isa Reactant && name(bm) == name(ligand) &&
            (binding_group[(from_species(s), to_species(s))] = h)
    end
    class = Dict{Species, Any}(complex => label(g) for complex in keys(twin_of))
    for s in groups[g]
        t = twin_of[to_species(s)]
        haskey(class, t) && continue
        h = count(==(t), values(twin_of)) == 1 ?
            get(binding_group, (from_species(s), t), 0) : 0
        class[t] = label(h > 0 ? h : t)
    end
    σ(sp) = get(class, sp, nothing)
    rescale(s) = (a = σ(from_species(s)); b = σ(to_species(s));
                  is_equilibrium(s) && isequal(a, b) ? (nothing, nothing) : (a, b))
    rescaling = Dict{Int, Tuple{Any, Any}}()
    for (h, group) in enumerate(groups)
        isempty(group) && continue
        r = rescale(first(group))
        all(s -> isequal(rescale(s), r), group) || return nothing
        rescaling[h] = r
    end
    rescaling
end

"""
    _redundant_copy_groups(m) -> Vector{Int}

Kinetic groups of `m` that bind a competitive inhibitor redundantly: in every
conformational state where the copy binds, every complex has a productive twin, a form
whose weight is proportional to the complex's (`_productive_twin`, `_all_twin`), and the
dwell gauge exists over the states together (`_gauge_rescaling`). The gauge then moves
the copy's constant along a direction the rate law cannot see, so the copy
adds a phantom, and the mechanism's family is its family without the copy
(`_gauge_rescaling`). A copy group whose complexes all have twins but whose gauge fails is
not redundant: a group that forms or leaves a twin also binds where no copy does (a shared
group that pins the existing binding), a group holds a step that leaves a complex beside
one that leaves another form (a second copy bound at the complex and at a twin), or a
mirror is steady-state. The gauge is the proof that the copy's constant is invisible, and
without it the constant may be identifiable. Two copies that are twins only of each other,
neither pinned, pass the test and carry a phantom: the rule admits them.

An allosteric copy binds the active state always and the inactive state unless its tag
is `:OnlyA`. The inactive state's graph is that of `_state_mechanism(am, :I)`: the steps
of the non-`:OnlyA` groups between forms reachable from the free enzyme
(`_inactive_groups`). The copy binds there only at the sites that graph keeps, and a
free-enzyme site is always among them. In a state where the copy binds nothing, every
factor of that state is 1. The two states form one system. A group tagged `:EqualAI`
has one set of constants in both, so it must take the same rescaling in each. A shared
class is one number in both states only when its binding group and the copy are both
`:EqualAI` (one K_h, one K*). The copy's complexes share one factor in both states when
the copy is `:EqualAI` (one K*). Every other class is a number per state.
"""
function _redundant_copy_groups(m::Union{Mechanism, AllostericMechanism})
    active = steps(m)
    any(group -> bound_metabolite(first(group)) isa CompetitiveInhibitor, active) ||
        return Int[]
    twin_active = _productive_twin(active)
    m isa Mechanism &&
        return [g for g in eachindex(active)
                if _gauge_rescaling(active, g, twin_active, identity) !== nothing]
    tags = cat_allo_states(m)
    inactive = _inactive_groups(m)
    twin_inactive = _productive_twin(inactive)
    filter(collect(eachindex(active))) do g
        one_number(c) = c isa Int && tags[c] === :EqualAI && tags[g] === :EqualAI
        label(state) = c -> one_number(c) ? c : (state, c)
        rescaling_active = _gauge_rescaling(active, g, twin_active, label(:A))
        rescaling_active === nothing && return false
        binds_inactive = !isempty(inactive[g])
        rescaling_inactive = binds_inactive ?
            _gauge_rescaling(inactive, g, twin_inactive, label(:I)) :
            Dict(h => (nothing, nothing) for h in eachindex(inactive)
                 if !isempty(inactive[h]))
        rescaling_inactive === nothing && return false
        all(h -> tags[h] !== :EqualAI || !haskey(rescaling_inactive, h) ||
                 isequal(rescaling_active[h], rescaling_inactive[h]),
            keys(rescaling_active))
    end
end

"""
Names of the forms where the metabolite named `met_name`, in the role `role`, is free on
some step: the step's `from_species` when it takes the metabolite up, its `to_species`
when it gives it off, so a Theorell–Chance step E(A) + B → E(Q) + P frees B at E(A) and
P at E(Q). `Reactant` asks where a substrate or product binds productively;
`CompetitiveInhibitor` where an inhibitor, a reactant's copy included, binds its dead-end
site. A reactant and its copy share a name but bind different sites, so the role tells
them apart. The dead-end move targets these forms for competition with the metabolite.
"""
function _forms_where_free(m::Union{Mechanism, AllostericMechanism},
                           role::Type{<:Metabolite}, met_name::Symbol)
    forms = Set{Symbol}()
    for group in steps(m), s in group,
        (form, mets) in ((from_species(s), consumed(s)), (to_species(s), released(s)))
        any(x -> x isa role && name(x) == met_name, mets) && push!(forms, name(form))
    end
    forms
end

"""
The dead-end move's child of `m`: `groups`, `m`'s groups with the new copy's mirrors
and the copy's group last, on reaction `rxn`; `nothing` when the new copy is redundant
in it (`_redundant_copy_groups`). `twin` is the `_productive_twin` of `m`'s steps, the
graph of its only (or active) state, and finds the child's twins there before the child
exists: the copy's forms are not productive, and its steps attach each complex to its
site's segment and join no two segments, so every form of `m` keeps its segment and
offsets. A `Mechanism` candidate takes the `_gauge_rescaling` path on `groups` and is
built only when kept. Only the new group is tested: `m` holds no redundant copy group
(the parent rule), and the placement cannot make an older one redundant. An allosteric
child is built, and runs `_redundant_copy_groups` only when every complex of the copy
has a twin in the active state (`_all_twin`), which redundancy needs; it tags the copy's
group `:EqualAI` and keeps `m`'s multiplicity and regulatory sites.
"""
function _dead_end_child(::Mechanism, groups::Vector{Vector{Step}}, rxn::EnzymeReaction,
                         twin)
    _gauge_rescaling(groups, length(groups), twin, identity) !== nothing && return nothing
    Mechanism(rxn, groups)
end
function _dead_end_child(am::AllostericMechanism, groups::Vector{Vector{Step}},
                         rxn::EnzymeReaction, twin)
    child = AllostericMechanism(rxn, groups, vcat(cat_allo_states(am), [:EqualAI]),
                                catalytic_multiplicity(am), copy(regulatory_sites(am)))
    _all_twin(groups[end], twin) && !isempty(_redundant_copy_groups(child)) ?
        nothing : child
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
`CompetitiveInhibitor` (its binding steps' source forms) and excludes the forms
that hold I; the capacity test counts productive bindings only. A form that
carries only the copy of M is therefore a site for an inhibitor that competes
with M. Mirror steps inherit their catalytic counterpart's `kinetic_group`. All
new binding steps for a single regulator share one fresh trailing kinetic group
(one K_R parameter).
A pattern is skipped when its child binds the new copy redundantly
(`_redundant_copy_groups`). A `Mechanism` candidate is tested by `_gauge_rescaling`
inside `_dead_end_child`, before the child is built; an allosteric child is built and
tested directly. Ligands of an allosteric mechanism's regulatory sites are not eligible:
they belong on a regulatory site.

The caller must pass the declared `rxn`: its regulators, not those of `m`'s reaction,
decide which inhibitors are eligible. The child's reaction is `rxn` itself.
"""
function _expand_add_dead_end_regulator(m::Union{Mechanism, AllostericMechanism},
                                        rxn::EnzymeReaction)
    sub_names = Set(name(s) for s in substrates(rxn))
    prod_names = Set(name(p) for p in products(rxn))
    bound_regs = Set{Symbol}(name(bm) for group in steps(m) for s in group
                             for bm in consumed(s) if bm isa Regulator)
    excluded = union(bound_regs, _bound_allo_regs(m))
    eligible_regs = [name(regulator(rm)) for rm in regulators(rxn)
                     if regulator(rm) isa CompetitiveInhibitor &&
                        name(regulator(rm)) ∉ excluded]
    isempty(eligible_regs) && return typeof(m)[]

    form_sp = Dict(name(sp) => sp for group in steps(m) for s in group
                   for sp in (from_species(s), to_species(s)))
    # Per form, the names bound productively and the names bound as competitive
    # inhibitors: a reactant's copy sits at a dead-end site, not at its reactant's.
    productive = Dict(f => Set(name(b) for b in bound(sp) if b isa Reactant)
                      for (f, sp) in form_sp)
    inhibiting = Dict(f => Set(name(b) for b in bound(sp) if b isa CompetitiveInhibitor)
                      for (f, sp) in form_sp)

    # The sites of each competition pattern, sorted, once each in pattern order.
    placements = Vector{Symbol}[]
    for (comp_subs, comp_prods, comp_inhibitors) in
            _inhibitor_competition_patterns(sub_names, prod_names, collect(bound_regs))
        comp_reactants = union(comp_subs, comp_prods)
        target_forms = Set{Symbol}()
        for met in comp_reactants
            union!(target_forms, _forms_where_free(m, Reactant, met))
        end
        for inh in comp_inhibitors
            union!(target_forms, _forms_where_free(m, CompetitiveInhibitor, inh))
        end
        push!(placements, [f for f in sort!(collect(target_forms))
                           if !issubset(sub_names, productive[f]) &&
                              !issubset(prod_names, productive[f]) &&
                              isempty(intersect(productive[f], comp_reactants)) &&
                              isempty(intersect(inhibiting[f], comp_inhibitors))])
    end
    unique!(filter!(!isempty, placements))

    results = typeof(m)[]
    twin = _productive_twin(steps(m))
    for reg_name in eligible_regs, active in placements
        inhibitor = CompetitiveInhibitor(reg_name)
        # The complex the copy forms at each site, in the site's conformation and residual.
        complex = Dict(f => Species(Metabolite[bound(form_sp[f])..., inhibitor],
                                    conformation(form_sp[f]), residual(form_sp[f]))
                       for f in active)
        mirrored(s) = haskey(complex, name(from_species(s))) &&
                      haskey(complex, name(to_species(s)))
        groups = [Step[group; [Step(complex[name(from_species(s))],
                                    complex[name(to_species(s))],
                                    consumed(s), released(s), is_equilibrium(s))
                               for s in group if mirrored(s)]]
                  for group in steps(m)]
        push!(groups, [Step(form_sp[f], complex[f], Metabolite[inhibitor], Metabolite[],
                            true) for f in active])
        child = _dead_end_child(m, groups, rxn, twin)
        child === nothing || push!(results, child)
    end
    results
end

"""
    _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

Convert a non-allosteric `Mechanism` into allosteric variants. Two kinds whose
conformational constant `L` never shows are not enumerated: the all-`:EqualAI` baseline
and a V-type variant without a regulator. An `:OnlyA` catalytic binding means the
inactive conformation cannot bind that metabolite, so it cannot complete the catalytic
cycle: every emitted
`:OnlyA` variant is **dead-inactive** — every chemistry group is `:OnlyA`,
and the inactive conformation only binds ligands. A chemistry group is one
holding a chemistry step (`_is_chemistry`): an isomerization, a fused
binding or a Theorell–Chance step; every other group is a binding group.

  * The all-`:EqualAI` baseline is never emitted — the two conformations
    are identical, `L` cancels, and the mechanism is indistinguishable
    from `m`.
  * K-type: every non-empty subset of binding groups is set `:OnlyA`,
    with every chemistry group `:OnlyA`, and each is emitted bare. Over
    more than one catalytic subunit the conformational equilibrium enters
    the rate with the multiplicity's power. Over one subunit the inactive
    conformation is a dead-end branch of the free enzyme, and `L` is a
    phantom in some K-type variants and not in others. In every K-type
    variant of the rapid-equilibrium ordered bi-bi seed a rescaling of the
    active state's constants absorbs it, whatever that conformation binds.
    In the merged ordered bi-bi whose B group mixes a fused and a plain
    binding it is a phantom only in the variants whose inactive
    conformation binds nothing: no rescaling absorbs an A or Q binding
    kept there. Binding nothing is not enough: in a random-order merged
    bi-bi, whose A group takes A up at E and at E(B) under one rate
    constant, `L` shows in every K-type variant. A subset that leaves a
    binding-only Wegscheider cycle unsatisfiable
    (`_onlya_haldane_violation`) is dropped.
  * V-type: no `:OnlyA` binding, every chemistry group `:OnlyA`. The
    inactive state binds substrate/product identically to the active
    state but cannot catalyze, so `L` folds entirely into `kcat`
    (`v = kcat/(1+L)·shape`) and is not observable; it is emitted ONLY
    paired with a declared allosteric regulator at a new site, one
    variant per `(regulator, tag)` with `tag ∈ {:OnlyA, :OnlyI}`. A
    reaction with no declared allosteric regulators emits no V-type.

For each value in `rxn`'s `allowed_catalytic_multiplicities`, the move emits
the K-types, then the V-types, at that `catalytic_multiplicity`. Catalytic
steps are reused by reference. The variants are distinct by construction: each
differs from the others in its tags, its regulatory site or its multiplicity.

A parent whose catalytic scheme fails `_hyperbolic_catalysis` (random-order
steady-state binding, or a substrate that traps a steady-state intermediate in
an abortive complex, whose own equation carries concentration powers) emits
variants at multiplicity 1 only: above one subunit the conformational
equilibrium would add a second source of powers
(`_requires_hyperbolic_catalysis`).
"""
function _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
    hyperbolic = _hyperbolic_catalysis(m)
    cns = [cn for cn in allowed_catalytic_multiplicities(rxn) if cn == 1 || hyperbolic]
    isempty(cns) && return AllostericMechanism[]
    chem = [any(_is_chemistry, group) for group in steps(m)]
    vtags = [c ? :OnlyA : :EqualAI for c in chem]
    regs = [name(regulator(rm)) for rm in regulators(rxn)
            if regulator(rm) isa AllostericRegulator]
    # K-type: every non-empty subset of binding groups :OnlyA, with every
    # chemistry group :OnlyA — a catalytically-dead inactive conformation.
    # A state that cannot bind a catalytic metabolite cannot complete the
    # cycle, so it runs no chemistry. Each is emitted bare; over one subunit
    # L can be a phantom (see the docstring). `_onlya_haldane_violation` drops
    # a subset that leaves a binding-only Wegscheider cycle unsatisfiable.
    ktags = Vector{Symbol}[]
    for sel in _subsets(findall(!, chem))[2:end]
        tags = copy(vtags)
        tags[sel] .= :OnlyA
        _onlya_haldane_violation(reaction(m), steps(m), tags) === nothing &&
            push!(ktags, tags)
    end
    results = AllostericMechanism[]
    for cn in cns
        for tags in ktags
            push!(results,
                  AllostericMechanism(reaction(m), steps(m), tags, cn, RegulatorySite[]))
        end
        # V-type: no :OnlyA binding, every chemistry group :OnlyA. The inactive
        # state binds identically but cannot catalyze, so L folds into kcat and is
        # unobservable — emit only paired with a declared regulator.
        any(chem) || continue
        for reg in regs, tag in (:OnlyA, :OnlyI)
            push!(results, AllostericMechanism(reaction(m), steps(m), vtags, cn,
                [RegulatorySite([AllostericRegulator(reg)], cn, [tag])]))
        end
    end
    results
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
The mechanism's catalytic side, regulatory_sites' multiplicities, and
reaction payload pass through unchanged.

The regulators declared in `rxn`, not those of `am`'s reaction, decide which
regulators are eligible: the reaction of an `@allosteric_mechanism` fixture
declares only the regulators its sites bind. The child keeps `am`'s reaction.
"""
function _expand_add_allosteric_regulator(
    am::AllostericMechanism, rxn::EnzymeReaction,
)
    taken = _bound_allo_regs(am)
    sites = regulatory_sites(am)
    results = AllostericMechanism[]
    for rm in regulators(rxn)
        reg = regulator(rm)
        reg isa AllostericRegulator && !(name(reg) in taken) || continue
        # Non-:EqualAI tags at any (new or existing) site.
        for tag in (:OnlyA, :OnlyI, :NonequalAI), site_idx in 0:length(sites)
            # Appending to an existing site whose conformations are disjoint
            # from the new ligand's (an all-:OnlyA site gaining an :OnlyI
            # ligand, or the reverse) reproduces the add-at-a-new-site
            # equation — redundant, so skip it.
            site_idx >= 1 && isempty(intersect(_state_conformations(tag),
                _site_active_states(sites[site_idx]))) && continue
            push!(results, _make_am_with_added_reg(am, name(reg), tag, site_idx))
        end
        # :EqualAI only at an existing site that keeps a non-:EqualAI ligand: an
        # all-:EqualAI site binds both conformations alike and cancels from the rate.
        for site_idx in findall(site -> any(!=(:EqualAI), allo_states(site)), sites)
            push!(results, _make_am_with_added_reg(am, name(reg), :EqualAI, site_idx))
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
    am::AllostericMechanism, reg::Symbol, tag::Symbol, site_idx::Int,
)
    sites = copy(regulatory_sites(am))
    lig = AllostericRegulator(reg)
    if site_idx == 0
        push!(sites, RegulatorySite([lig], catalytic_multiplicity(am), [tag]))
    else
        s = sites[site_idx]
        sites[site_idx] =
            RegulatorySite([ligands(s); lig], multiplicity(s), [allo_states(s); tag])
    end
    _with(am; sites)
end

"""
    _expand_change_allo_state(am::AllostericMechanism)
        → Vector{AllostericMechanism}

Relax a "constrained" tag (`:EqualAI`, `:OnlyA`, `:OnlyI`) to `:NonequalAI`. A binding
catalytic group and a regulatory ligand each relax individually — one variant per group
not already `:NonequalAI`. The chemistry groups (those holding a chemistry step,
`_is_chemistry`) relax **together**, one variant setting every
non-`:NonequalAI` chemistry group to `:NonequalAI` at once: inactive catalysis
is all-or-nothing, so a fully-`:NonequalAI` catalytic inactive conformation
cannot be reached by relaxing chemistry groups one at a time (each mixed
intermediate is a partial and is dropped). The base catalytic steps,
multiplicity, and untouched tags are preserved.

The chemistry relaxation is tried only while no binding group is `:OnlyA`. The
inactive conformation's catalysis must be all-or-nothing — fully dead (every
chemistry group `:OnlyA`, whether because an `:OnlyA` binding blocks the cycle or by
a dead-inactive V-type) or fully live — and the relaxation leaves every chemistry
group live. Beside an `:OnlyA` binding it would leave a partial conformation, which
strands enzyme in a covalent form (a kinetic sink) and crashes the
saturating-turnover extraction; a one-sided `:OnlyA` binding, legal only because
`k_I = 0`, would also be stranded by the restored finite `k_I`, leaving no
thermodynamic reading (`_onlya_haldane_violation`). Where the inactive conformation
binds nothing (all bindings `:OnlyA`), the skipped variant is rate-equivalent to the
fully-dead form emitted directly, so no observable hypothesis is lost.

A binding relaxation is dropped when it leaves an `:OnlyA` binding's cycle
unsatisfiable (`_onlya_haldane_violation`). It needs no all-or-nothing test: it tags
no group `:OnlyA` and leaves the chemistry tags alone, so a parent whose inactive
catalysis is fully dead or fully live stays so. Every mechanism the moves emit is one
of the two; a hand-built partial parent keeps its partial binding relaxations.

A regulatory ligand's tag is not an argument to the Haldane check — a
regulator site completes no catalytic cycle — so that branch needs no
filter.

A child with a redundant competitive-inhibitor group (`_redundant_copy_groups`)
is dropped. The per-state test reads the tags this move changes: relaxing an
`:OnlyA` binding to `:NonequalAI` brings its complex into the inactive state,
where a copy that was new there can become its twin with a consistent gauge, and
the copy's constant then shows in neither state.
"""
function _expand_change_allo_state(am::AllostericMechanism)
    cs, states = steps(am), cat_allo_states(am)
    chem = [any(_is_chemistry, group) for group in cs]
    # Binding catalytic groups relax individually. Chemistry groups relax together.
    # Inactive catalysis is all-or-nothing, so a fully-`:NonequalAI` catalytic inactive
    # conformation is unreachable by relaxing chemistry groups one at a time — each
    # mixed intermediate is a partial and is dropped. One variant sets every
    # non-`:NonequalAI` chemistry group to `:NonequalAI` at once, and only while no
    # binding group is `:OnlyA`, which would leave the inactive catalysis partial.
    relaxations = [[g] for g in eachindex(cs) if !chem[g] && states[g] != :NonequalAI]
    any(chem .& (states .!= :NonequalAI)) && !any(.!chem .& (states .== :OnlyA)) &&
        push!(relaxations, findall(chem))
    results = AllostericMechanism[]
    for gs in relaxations
        new_states = copy(states)
        new_states[gs] .= :NonequalAI
        _onlya_haldane_violation(reaction(am), cs, new_states) === nothing || continue
        push!(results, _with(am; states = new_states))
    end
    for (si, site) in enumerate(regulatory_sites(am)), li in eachindex(ligands(site))
        allo_states(site)[li] == :NonequalAI && continue
        new_sites = copy(regulatory_sites(am))
        new_states = copy(allo_states(site))
        new_states[li] = :NonequalAI
        new_sites[si] = RegulatorySite(copy(ligands(site)), multiplicity(site), new_states)
        push!(results, _with(am; sites = new_sites))
    end
    filter!(c -> isempty(_redundant_copy_groups(c)), results)
end

"""
    _state_conformations(state::Symbol) -> Set{Symbol}

The conformations a single ligand's allosteric `state` acts on: `:active` for
`:OnlyA`/`:EqualAI`/`:NonequalAI`, `:inactive` for `:OnlyI`/`:EqualAI`/`:NonequalAI`.
"""
function _state_conformations(state::Symbol)
    conf = Set{Symbol}()
    state in (:OnlyA, :EqualAI, :NonequalAI) && push!(conf, :active)
    state in (:OnlyI, :EqualAI, :NonequalAI) && push!(conf, :inactive)
    conf
end

"""
    _site_active_states(site::RegulatorySite) -> Set{Symbol}

The conformations a regulatory site's ligands act on — the union of
`_state_conformations` over its states. Two sites with disjoint active states —
an all-`:OnlyA` site and an all-`:OnlyI` site — merge to a rate equation
identical to keeping them separate, so that all-keep merge is redundant.
"""
function _site_active_states(site::RegulatorySite)
    active = Set{Symbol}()
    for st in allo_states(site)
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
load-bearing. The merged site reuses one site's
`multiplicity` (equal to `catalytic_multiplicity`) and its ligands are sorted
by name, so two merge routes reaching the same ligand partition produce `==`
mechanisms and dedup by `hash`. Regulator type is not enforced here;
`expand_mechanisms` drops every child that fails `_respects_reg_type`.
"""
function _expand_merge_regulatory_sites(am::AllostericMechanism)
    sites = regulatory_sites(am)
    n = length(sites)
    results = AllostericMechanism[]
    for i in 1:(n - 1), j in (i + 1):n
        ligs = vcat(ligands(sites[i]), ligands(sites[j]))
        base_states = vcat(allo_states(sites[i]), allo_states(sites[j]))
        perm = sortperm(ligs; by = lig -> String(name(lig)))
        ligs = ligs[perm]
        base_states = base_states[perm]
        mult = multiplicity(sites[i])
        others = RegulatorySite[sites[k] for k in 1:n if k != i && k != j]
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
            merged = RegulatorySite(copy(ligs), mult, states)
            push!(results, _with(am; sites = vcat(others, [merged])))
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
_respects_reg_type(::Mechanism, ::EnzymeReaction) = true
function _respects_reg_type(am::AllostericMechanism, rxn::EnzymeReaction)
    declared = Dict(name(regulator(rm)) => reg_type(rm) for rm in regulators(rxn)
                    if regulator(rm) isa AllostericRegulator)
    all(regulatory_sites(am)) do site
        types = [get(declared, name(l), :unspecified) for l in ligands(site)]
        all(zip(types, allo_states(site))) do (t, state)
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
competitive-inhibitor group is redundant (`_redundant_copy_groups`).
"""
function expand_mechanisms(
    mechs::Vector{<:Union{Mechanism, AllostericMechanism}},
    rxn::EnzymeReaction)
    result = Union{Mechanism, AllostericMechanism}[]
    for m in mechs
        _assert_emission_rules(m)
        append!(result, _expand_re_to_ss(m), _expand_split_kinetic_group(m),
                _expand_add_dead_end_regulator(m, rxn))
        if m isa Mechanism
            append!(result, _expand_to_allosteric(m, rxn))
        else
            append!(result, _expand_add_allosteric_regulator(m, rxn),
                    _expand_change_allo_state(m), _expand_merge_regulatory_sites(m))
        end
    end
    filter!(c -> _respects_reg_type(c, rxn), result)
    for child in result
        _assert_atom_conserving(child)
    end
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
the reaction has no catalytic topology.
"""
function init_mechanisms(r::EnzymeReaction)
    seeds = _expand_substrate_product_dead_ends(_catalytic_topologies(r), r)
    mechs = [Mechanism(r, _seed_groups(steps)) for steps in seeds]
    foreach(_assert_atom_conserving, mechs)
    seen = Set(mechs)
    out = copy(mechs)
    for m in unique(mechs), v in _seed_variants(m)
        v in seen || (push!(seen, v); push!(out, v))
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
traversal. Errors when no mechanism binds every required regulator.
"""
function seed_mechanisms(rxn::EnzymeReaction, required_allo::Set{Symbol},
                         required_comp::Set{Symbol})
    visited = Set{UInt64}()
    seeds = Union{Mechanism, AllostericMechanism}[]
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
    frontier = Union{Mechanism, AllostericMechanism}[
        m for m in init_mechanisms(rxn) if consider!(m)]
    while !isempty(frontier)
        # Per-node child generation is pure and independent — distribute it.
        childsets = pmap(frontier) do m
            children = Union{Mechanism, AllostericMechanism}[]
            isempty(required_allo) || append!(children, m isa Mechanism ?
                _expand_to_allosteric(m, rxn) : _expand_add_allosteric_regulator(m, rxn))
            append!(children, _expand_add_dead_end_regulator(m, rxn))
            filter!(c -> _is_seed_node(c, rxn, required_allo, required_comp), children)
        end
        next = Union{Mechanism, AllostericMechanism}[]
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
    seeds
end

"""Names of the allosteric regulators bound at `m`'s regulatory sites (empty for
a `Mechanism`)."""
_bound_allo_regs(::Mechanism) = Set{Symbol}()
_bound_allo_regs(am::AllostericMechanism) =
    Set{Symbol}(name(lig) for site in regulatory_sites(am) for lig in ligands(site))

"""Names of the competitive inhibitors bound at a dead-end step in `m`."""
function _bound_comp_inhibitors(m::Union{Mechanism, AllostericMechanism})
    bound = Set{Symbol}()
    for group in steps(m), s in group, bm in consumed(s)
        bm isa CompetitiveInhibitor && push!(bound, name(bm))
    end
    bound
end

"""
    _is_seed_node(m, rxn::EnzymeReaction, required_allo::Set{Symbol},
                 required_comp::Set{Symbol}) -> Bool

A child worth expanding: cheap states, one ligand per regulatory site, no
optional regulator bound, and reg-type-respecting.
"""
function _is_seed_node(m::Union{Mechanism, AllostericMechanism},
                       rxn::EnzymeReaction, required_allo::Set{Symbol},
                       required_comp::Set{Symbol})
    m isa AllostericMechanism && (:NonequalAI in cat_allo_states(m) ||
        any(site -> length(ligands(site)) != 1 || only(allo_states(site)) == :NonequalAI,
            regulatory_sites(m))) && return false
    issubset(_bound_allo_regs(m), required_allo) || return false
    issubset(_bound_comp_inhibitors(m), required_comp) || return false
    _respects_reg_type(m, rxn) || return false
    true
end
