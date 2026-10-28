# ABOUTME: The seed mechanisms of the enumeration: catalytic topologies, dead-end
# ABOUTME: decorations, merged and Theorell–Chance variants, and the entry points.

# ─── Catalytic topologies ─────────────────────────────────────────────

"""Whether raw step `s` takes up and gives off nothing: an isomerization."""
_is_iso(s::Mechanisms.RawStep) = isempty(s.consumed) && isempty(s.released)

"""Group labels renumbered 1, 2, … in order of first occurrence."""
function _numbered_groups(labels::Vector{Int})
    first_seen = unique(labels)
    [findfirst(==(l), first_seen) for l in labels]
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
        push!(get!(iso_groups, Set(filter(_is_iso, path)),
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
                          by = s -> (_is_iso(s), form_pair(s)...))

            # The first iso step is the (single) SS step; every other step is RE.
            iso_idx = findfirst(_is_iso, steps)
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
          (s.consumed, s.released, s.rapid_equilibrium, _is_iso(s) ? j : 0),
          length(label) + 1)
     for (j, s) in enumerate(steps)]
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
                            _numbered_groups(raw.group[keep]))
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
                            _numbered_groups(group))
end

"""The merged complexes of graph `m`, each once: the forms that a fused binding of a
substrate enters, a binding (`is_binding`) of kind `:chemistry` (`step_kind`) that takes up
a metabolite named in `subs`. Only the edges `on` flags count."""
_merged_complexes(m::Mechanisms.Mechanism, subs::Vector{Symbol},
                  on::BitVector = trues(length(Mechanisms.steps(m)))) =
    unique(Mechanisms.ends(m, k)[2] for k in Mechanisms.steps(m)
           if on[k] && Mechanisms.is_binding(m, k) &&
              _holds(Mechanisms.taken_up(m, k), subs) &&
              Mechanisms.step_kind(m, k) === :chemistry)

"""
The tests of `_seed_variants` that build nothing, for the graph `m` of one base of a
seed: a function of a candidate's rapid-equilibrium flags (one per edge of `m`, some
turned steady state) that is true when the candidate has no rapid-equilibrium turnover
cycle (`has_turnover_cycle`), a maximal rate both ways (`_has_vmax`) and no chemistry in
equilibrium with both sides (`_chemistry_equilibrates_both_sides`). Each edge's ends,
uptake weight and reactants are read once per base, so a candidate costs a few passes
over arrays. A base holds no isomerization, so each chemistry node is a merged complex
or a Theorell–Chance step. A set of edges holds a turnover cycle iff it closes a cycle of
nonzero weight (`has_unbalanced_cycle`), which reuses the screen's buffers. Each
maximal-rate test runs on the candidate's rapid-equilibrium edges and more, and adding
edges never removes such a cycle, so a turnover cycle of the candidate fails both tests.
"""
function _seed_candidate_screen(m::Mechanisms.Mechanism)
    ks = Mechanisms.steps(m)
    any(k -> Mechanisms.step_kind(m, k) === :isomerization, ks) &&
        error("_seed_candidate_screen: the seed base holds an isomerization; a base " *
              "has every isomerization merged")
    subs, prods = Mechanisms.substrates(m), Mechanisms.products(m)
    from = [Mechanisms.ends(m, k)[1] for k in ks]
    to = [Mechanisms.ends(m, k)[2] for k in ks]
    weight = [_uptake_weight(m, k) for k in ks]
    touches(side) = BitVector([_exchanges(m, k, side) for k in ks])
    sub_step, prod_step = touches(subs), touches(prods)
    crossing = BitVector([_crosses_sides(m, k, subs, prods) for k in ks])
    complexes = _merged_complexes(m, subs)
    leaves(x, side) = BitVector([to[k] == x && _holds(Mechanisms.taken_up(m, k), side) ||
                                 from[k] == x && _holds(Mechanisms.given_off(m, k), side)
                                 for k in ks])
    exits = [(leaves(x, subs), leaves(x, prods)) for x in complexes]
    n = length(Mechanisms.species(m))
    parent = zeros(Int, n); offset = zeros(Int, n)
    turnover(on) = Mechanisms.has_unbalanced_cycle(from, to, weight, on, parent, offset)
    meets(a::BitVector, b::BitVector) = any(k -> a[k] && b[k], eachindex(a))
    kept = falses(length(ks))
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
The bases of the seed `m`'s variants, as raw input (`RawMechanism`), the merged base first:
merging every isomerization of `m` onto its product side (`_merge_isomerization`) and
setting every step at rapid equilibrium gives the merged base; eliminating one merged
complex with two steps (`_eliminate_form`) gives a Theorell–Chance base, one per complex
that can be eliminated, in the order of the isomerizations.
"""
function _seed_bases(m::Mechanisms.Mechanism)
    raw = Mechanisms.RawMechanism(m)
    isos = filter(_is_iso, raw.steps)
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
    at = findall(k -> f in Mechanisms.ends(base, k), Mechanisms.steps(base))
    length(at) == 2 &&
        all(k -> !re[k] && length(Mechanisms.members(base, Mechanisms.group(base, k))) == 1,
            at)
end

"""
The merged and Theorell–Chance variants of the seed `m` of reaction `rxn`, whose atoms
each base's edges must conserve. Each base (`_seed_bases`) is read as its canonical graph
(`canonical_graph`), unchecked; a base whose tied groups no orientation aligns raises
the representation error. MERGE and ELIM conserve atoms, which each
base asserts. For each base, every inclusion-minimal set of its groups (`_minimal_flips`)
whose flip to steady state gives a valid, flux-carrying, non-degenerate candidate is a
variant: no rapid-equilibrium turnover cycle (`has_turnover_cycle`), a maximal rate both
ways (`_has_vmax`), no chemistry in equilibrium with both sides
(`_chemistry_equilibrates_both_sides`), every steady-state group carrying flux, and the
guard (`_guard_accepts`). The first three tests read index arrays
(`_seed_candidate_screen`); the other two run only on a candidate that passes them.
A merged variant that is its seed's lumping twin (`_lumping_twin`) is skipped: the unmerged
form with rapid-equilibrium flanks has its family at the same count. A base's groups come in
canonical order, so the variants come in the order of their flip sets over those groups;
each is built and checked.
"""
function _seed_variants(m::Mechanisms.Mechanism, rxn::Mechanisms.EnzymeReaction)
    variants = Mechanisms.Mechanism[]
    for (i, b) in enumerate(_seed_bases(m))
        base = Mechanisms.canonical_graph(b)
        _assert_atom_conserving(base, rxn)
        flips = _minimal_flips(_seed_candidate_screen(base), base)
        if i == 1
            cx = _merged_complexes(base, Mechanisms.substrates(base))
            filter!(re -> !_lumping_twin(base, cx, re), flips)
        end
        append!(variants, [Mechanisms.with_flags(base, re) for re in flips])
    end
    variants
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
(`_guard_accepts`) is dropped; the flip search has run the guard on every
variant.
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
    out = filter(_guard_accepts, mechs)
    for m in unique(mechs), v in _seed_variants(m, r)
        v in seen && continue
        push!(seen, v)
        push!(out, v)
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
law the guard rejects (`_guard_accepts`) is dropped. The search still expands such
a node: only the seeds it returns are checked.
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
    filter!(_guard_accepts, seeds)
    seeds
end

"""Names of the allosteric regulators bound at `m`'s regulatory sites (empty for
a `Mechanism`)."""
_bound_allo_regs(::Mechanisms.Mechanism) = Set{Symbol}()
_bound_allo_regs(am::Mechanisms.AllostericMechanism) =
    Set{Symbol}(lig for site in Mechanisms.regulatory_sites(am)
                for lig in Mechanisms.ligands(site))

"""Names of the competitive inhibitors some step of `m` takes up as a modifier
(`modifier_copies`)."""
function _bound_comp_inhibitors(m::Union{Mechanisms.Mechanism,
                                          Mechanisms.AllostericMechanism})
    cm = Mechanisms.catalytic(m)
    Set{Symbol}(n for k in Mechanisms.steps(cm)
                for n in Mechanisms.metabolite_names(
                    Mechanisms.modifier_copies(Mechanisms.taken_up(cm, k))))
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

"""
    base_mechanisms(rxn, required_allo::Set{Symbol}, required_comp::Set{Symbol})
        -> (mechanisms, failures::Vector{Tuple{Any, String}})

The mechanisms the beam fits first. Both sets empty means they start from
`init_mechanisms`; a non-empty set means they start from `seed_mechanisms`. Either way the
base tier fits the seeds that are not degenerate and, in place of those that are, their
flip children that are not (`_base_tier`). Each failure is a `(seed, message)` pair.
"""
function base_mechanisms(rxn::Mechanisms.EnzymeReaction, required_allo::Set{Symbol},
                         required_comp::Set{Symbol})
    seeds = (isempty(required_allo) && isempty(required_comp)) ?
        unique!(collect(init_mechanisms(rxn))) :
        unique!(collect(seed_mechanisms(rxn, required_allo, required_comp)))
    _base_tier(seeds, rxn)
end

"""
The mechanisms the beam fits first: every seed of `mechs` that is not degenerate
(`_degenerate`), and in place of each one that is, its flip children (`_expand_re_to_ss`)
that are not. A degenerate seed's law ignores a substrate or never saturates, so it is not
worth a fit, but it is the only parent of mechanisms that are not degenerate (such as the
ping-pong seeds' children). Only a flip can cure it: the other moves keep every catalytic
step's rapid-equilibrium flag, and the steps they add bind regulators, which neither
degeneracy condition reads. An expansion error is returned as a failure, one per seed,
carrying the seed and the first 200 characters of the error. The flip search emits only
children the guard accepts (`_guard_accepts`).
"""
function _base_tier(mechs::Vector, rxn::Mechanisms.EnzymeReaction)
    base = Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}[
        m for m in mechs if !_degenerate(m)]
    failures = Tuple{Any, String}[]
    for m in mechs
        _degenerate(m) || continue
        try
            for child in _expand_re_to_ss(m)
                _assert_atom_conserving(child, rxn)
                _degenerate(child) || push!(base, child)
            end
        catch e
            push!(failures, (m, first(sprint(showerror, e), 200)))
        end
    end
    unique!(base), failures
end
