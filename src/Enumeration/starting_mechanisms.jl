# ABOUTME: The starting mechanisms of a reaction, candidate_mechanisms(reaction): catalytic
# ABOUTME: topologies, dead-end decorations, merged and Theorell–Chance variants, seeds.

# ─── Catalytic topologies ───────────────────────────────────────────────

"""`steps` with their group labels renumbered 1, 2, … in order of first occurrence."""
function _numbered_groups(steps::Vector{Mechanisms.Step})
    number = Dict{Int, Int}()
    Mechanisms.Step[Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off,
                                    e.rapid_equilibrium,
                                    get!(number, e.kinetic_group, length(number) + 1))
                    for e in steps]
end

"""
The written graph of `steps` between the forms whose records are `forms`, as
`(; forms, steps)`: the records of the forms some step joins, in the order of `forms`, and
`steps` with their ends indexing them.
"""
function _trimmed(forms::Vector{_FormRecord}, steps::Vector{Mechanisms.Step})
    used = sort!(unique!([x for e in steps for x in (e.from, e.to)]))
    at = zeros(Int, length(forms))
    at[used] = eachindex(used)
    (; forms = forms[used],
     steps = Mechanisms.Step[Mechanisms.Step(at[e.from], at[e.to], e.taken_up, e.given_off,
                                             e.rapid_equilibrium, e.kinetic_group)
                             for e in steps])
end

"""
All weak orderings of `items`: sequences of nonempty levels that partition `items`, each
level in the order of `items`. They come in the order of their first level in
`_subsets(items)`, then recursively in the order of the rest.
"""
function _weak_orderings(items::Vector{Symbol})
    isempty(items) && return [Vector{Symbol}[]]
    orderings = Vector{Vector{Symbol}}[]
    # A comprehension here would hide the recursive call's type from inference.
    for level in _subsets(items)[2:end], rest in _weak_orderings(setdiff(items, level))
        push!(orderings, [[level]; rest])
    end
    orderings
end

"""
    _catalytic_topologies(reaction) -> Vector{@NamedTuple{forms, steps}}

The catalytic topologies of `reaction`, each a written graph (`_trimmed`) with each of its
forms once: steps forming one or more complete catalytic cycles (E -> ... -> E), every
release written as the binding it reverses, so every binding takes its metabolite up. A
step's lists index the reaction's substrates, then its products. Each step's group is 0,
for `_seed_groups` to label. Errors when the reaction has none.
"""
function _catalytic_topologies(reaction::Mechanisms.EnzymeReaction)
    sub_names = Symbol[s.name for s in reaction.substrates]
    prod_names = Symbol[p.name for p in reaction.products]
    reactants = [sub_names; prod_names]
    reactant_index(x::Symbol) = findfirst(==(x), reactants)::Int
    forms = _FormRecord[]
    index = Dict{_FormRecord, Int}()

    # The index in `forms` of the record of the form on :E binding `on_subs` and
    # `on_prods` after a route consumed `consumed` and released `released`, the form added
    # when new. Its covalent residual adds the consumed substrates not bound and subtracts
    # the released and bound products; it has none exactly when their atoms cancel (no
    # covalent residue remains).
    function form(consumed, released, on_subs, on_prods)
        added, subtracted = setdiff(consumed, on_subs), [released; on_prods]
        bound = sort!([on_subs; on_prods])
        n = isempty(_net_atoms(reaction, added, subtracted)) ?
            _FormRecord(:E, bound, Symbol[], Symbol[], Symbol[]) :
            _FormRecord(:E, bound, Symbol[], sort!(added), sort!(subtracted))
        get!(() -> (push!(forms, n); length(forms)), index, n)
    end

    # Collect all complete catalytic paths as step lists between `forms`. A route starts
    # at free E, binds substrates, isomerizes a substrate-bound form to a product-bound
    # one, releases those products one at a time in every order, and repeats until it has
    # consumed every substrate and released every product. `cur` is the index of the
    # current form, `consumed` and `released` are the route's history, and `on_subs` and
    # `on_prods` are what `cur` binds; the atoms on the enzyme are those of `consumed`
    # minus those of `released`. `walk!` takes itself as an argument because a closure that
    # names itself is boxed.
    paths = Vector{Mechanisms.Step}[]
    function walk!(self, cur, consumed, released, on_subs, on_prods, path)
        if isempty(forms[cur].productive) && length(consumed) == length(sub_names) &&
                length(released) == length(prod_names)
            push!(paths, copy(path))
            return
        end
        # The step from `cur` to the form the arguments describe, taking up `takes_up`,
        # or, for a release, the binding into `cur` it reverses; then the walk on from
        # that form.
        function step!(takes_up::Vector{Int}, release::Bool, con, rel, on_s, on_p)
            to = form(con, rel, on_s, on_p)
            push!(path,
                  release ? Mechanisms.Step(to, cur, takes_up, Int[], true, 0) :
                            Mechanisms.Step(cur, to, takes_up, Int[], true, 0))
            self(self, to, con, rel, on_s, on_p, path)
            pop!(path)
        end
        if !isempty(on_prods)
            # Release a bound product; nothing binds until every bound product is gone.
            for p in on_prods
                step!([reactant_index(p)], true, consumed, [released; p], Symbol[],
                      filter(!=(p), on_prods))
            end
            return
        end
        remaining_subs = setdiff(sub_names, consumed)
        for s in remaining_subs
            step!([reactant_index(s)], false, [consumed; s], released, [on_subs; s],
                  Symbol[])
        end
        # Only a substrate-bound form isomerizes (rule C7), and an isomerization turns at
        # most three substrates into at most three products (rule C6).
        (isempty(on_subs) || length(on_subs) > 3) && return
        for subset in _subsets(setdiff(prod_names, released))[2:end]
            residue = _net_atoms(reaction, consumed, [released; subset])
            any(<(0), values(residue)) && continue
            # An isomerization leaves a covalent residue exactly when substrates remain to
            # bind (the admissible-residual rule). Without a residue while substrates
            # remain, the enzyme would return to apo E mid-cycle, splitting the reaction
            # into disconnected half-cycles. With a residue after the last substrate, the
            # route strands it, since only a substrate-bound form isomerizes (C7).
            isempty(residue) == isempty(remaining_subs) || continue
            # Rule C6 on the products, the residue counting as one.
            length(subset) + !isempty(residue) > 3 && continue
            # Its product side binds those products and no substrate (rule C8).
            step!(Int[], false, consumed, released, Symbol[], subset)
        end
    end
    walk!(walk!, form(Symbol[], Symbol[], Symbol[], Symbol[]), Symbol[], Symbol[],
          Symbol[], Symbol[], Mechanisms.Step[])

    iso_groups = Dict{Set{Mechanisms.Step}, Vector{Vector{Mechanisms.Step}}}()
    for path in paths
        push!(get!(iso_groups, Set(filter(_is_iso, path)),
                   Vector{Mechanisms.Step}[]), path)
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
    # their iso steps' forms in content order) so topology output order is stable — a
    # `Set`'s iteration order is not.
    level_of(wo) = Dict{Symbol,Int}(m => i for (i, level) in enumerate(wo) for m in level)
    sub_levels = level_of.(_weak_orderings(sub_names))
    prod_levels = level_of.(_weak_orderings(prod_names))
    # The names a path binds, in order, among `names`: a release counts as the binding it
    # is written as.
    route_order(path, names) =
        Symbol[reactants[s.taken_up[1]] for s in path
               if !isempty(s.taken_up) && reactants[s.taken_up[1]] in names]
    # A form's place in content order: free E first, then by label, bound metabolites,
    # residual additions and residual subtractions. A topology binds no modifier, so its
    # bound metabolites are its productive ones.
    key(x::Int) = (f = forms[x];
                   (!(isempty(f.productive) && isempty(f.adds) && isempty(f.subtracts) &&
                      f.label === :E),
                    f.label, f.productive, f.adds, f.subtracts))
    form_pair(s) = (key(s.from), key(s.to))
    sorted_iso_pats = sort!(collect(keys(iso_groups));
                            by = pat -> (length(pat), sort([form_pair(s) for s in pat])))
    result = @NamedTuple{forms::Vector{_FormRecord}, steps::Vector{Mechanisms.Step}}[]
    for iso_pat in sorted_iso_pats
        group = [(path, route_order(path, sub_names), route_order(path, prod_names))
                 for path in iso_groups[iso_pat]]
        seen_topos = Set{Set{Mechanisms.Step}}()
        for sl in sub_levels, pl in prod_levels
            topo_keys = Set{Mechanisms.Step}()
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
            push!(result, _trimmed(forms, Mechanisms.Step[
                Mechanisms.Step(s.from, s.to, s.taken_up, s.given_off, i != iso_idx, 0)
                for (i, s) in enumerate(steps)]))
        end
    end
    isempty(result) && error("no catalytic cycle for $(join(sub_names, " + ")) ⇌ " *
                             "$(join(prod_names, " + ")): every route needs an " *
                             "isomerization converting more than three substrates, or " *
                             "more than three products counting a covalent residue as one")
    result
end

# ─── Dead-end helpers ───────────────────────────────────────────────────

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
    _expand_substrate_product_dead_ends(topos, reaction)
        -> Vector{@NamedTuple{forms, steps}}

For each catalytic topology, a written graph, the mechanisms that add a set of dead-end
forms to it, each a written graph: the topology's forms and the dead-end forms after them,
each once, and the steps between them. A dead-end form is a catalytic form with one more
substrate or product bound that

- binds at least one substrate and at least one product;
- binds neither every substrate nor every product;
- never binds both metabolites of a declared `shared_catalytic_site` pair.

Each competition pattern (`_competition_patterns`) keeps the dead-end forms that bind
none of its (substrate, product) pairs, and each distinct set of kept forms gives one
mechanism: the topology plus a rapid-equilibrium binding step between every two of its
forms with the same conformation and residual that differ by one bound metabolite. A
step's lists index the reaction's substrates, then its products, as the topologies' do.

Two preconditions hold for the topologies of `_catalytic_topologies`: every binding is
rapid equilibrium, and every isomerization joins a substrate-only form to a
product-only form. The first gives a binding between two dead-end forms the flag of the
topology binding it parallels. The second keeps every dead-end form off the catalytic
cycle and leaves no isomerization a dead-end parallel, so bindings are the only steps
the dead-end forms need.
"""
function _expand_substrate_product_dead_ends(
    topos::Vector{@NamedTuple{forms::Vector{_FormRecord}, steps::Vector{Mechanisms.Step}}},
    reaction::Mechanisms.EnzymeReaction,
)
    sub_names = Set(s.name for s in reaction.substrates)
    prod_names = Set(p.name for p in reaction.products)
    mets = sort!(collect(union(sub_names, prod_names)))
    reactants = Symbol[[s.name for s in reaction.substrates];
                       [p.name for p in reaction.products]]

    # Competition patterns depend only on the reaction, so they are computed once.
    # A declared shared catalytic site forbids its (substrate, product) pair
    # from co-occupying the catalytic site, so keep only competition patterns
    # whose forbidden-edge set contains every declared pair. The complete
    # bipartite pattern contains all edges and always survives, so the list is
    # never empty.
    patterns = filter(pat -> all(in(pat), reaction.shared_catalytic_site),
                      _competition_patterns(sub_names, prod_names))

    result = @NamedTuple{forms::Vector{_FormRecord}, steps::Vector{Mechanisms.Step}}[]
    for topo in topos
        forms = copy(topo.forms)
        index = Dict(f => x for (x, f) in enumerate(forms))
        # A topology binds no modifier, so a form's bound names are its productive ones.
        held(x::Int) = Set(forms[x].productive)
        catalytic = unique!([x for s in topo.steps for x in (s.from, s.to)])
        dead_ends = Int[]
        for x in catalytic, m in mets
            b = held(x)
            m in b && continue
            push!(b, m)
            subs, prods = intersect(b, sub_names), intersect(b, prod_names)
            (isempty(subs) || isempty(prods) || subs == sub_names || prods == prod_names) &&
                continue
            f = forms[x]
            d = _FormRecord(f.label, sort!([f.productive; m]), f.copies, f.adds,
                            f.subtracts)
            push!(dead_ends, get!(() -> (push!(forms, d); length(forms)), index, d))
        end
        # Several catalytic forms can give the same dead end.
        unique!(dead_ends)
        edges = Set((s.from, s.to) for s in topo.steps)
        seen = Set{Vector{Int}}()
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
            xs = [catalytic; active]
            hs = held.(xs)
            steps = copy(topo.steps)
            for (i, x1) in enumerate(xs), (j, x2) in enumerate(xs)
                f1, f2 = forms[x1], forms[x2]
                f1.label === f2.label &&
                    (f1.adds, f1.subtracts) == (f2.adds, f2.subtracts) &&
                    length(hs[j]) == length(hs[i]) + 1 && issubset(hs[i], hs[j]) &&
                    !((x1, x2) in edges) || continue
                met = only(setdiff(hs[j], hs[i]))
                push!(steps, Mechanisms.Step(x1, x2, [findfirst(==(met), reactants)],
                                             Int[], true, 0))
            end
            push!(result, _trimmed(forms, steps))
        end
    end
    result
end

# ─── Mechanism enumeration ──────────────────────────────────────────────

"""
    _seed_groups(steps) -> Vector{Int}

The kinetic group label of each step of a seed, labels numbered in order of first
use: steps that take up and give off the same metabolites with the same RE/SS flag
share one group, since the steps of a group share their constants; every isomerization
is a group of its own.
"""
function _seed_groups(steps::Vector{Mechanisms.Step})
    label = Dict{Tuple{Vector{Int}, Vector{Int}, Bool, Int}, Int}()
    [get!(label,
          (s.taken_up, s.given_off, s.rapid_equilibrium, _is_iso(s) ? j : 0),
          length(label) + 1)
     for (j, s) in enumerate(steps)]
end

"""
`steps` with the isomerization `s0`, alone in its kinetic group, merged onto its product
side: `s0` is removed, and every other step at `s0.from` moves onto `s0.to`, so no step
joins `s0.from`. The last substrate's binding into the substrate side becomes a fused
binding into the product side; the releases stay plain. The other steps keep their group
labels, and the label of `s0`'s group goes unused.
"""
function _merge_isomerization(steps::Vector{Mechanisms.Step}, s0::Mechanisms.Step)
    move(x::Int) = x == s0.from ? s0.to : x
    Mechanisms.Step[Mechanisms.Step(move(s.from), move(s.to), s.taken_up, s.given_off,
                                    s.rapid_equilibrium, s.kinetic_group)
                    for s in steps if s != s0]
end

"""
`steps` with node `x` eliminated, or `nothing` when that is not possible: `x` must have
exactly two steps, both bindings into it, Y + L → x and Z + R → x with Y ≠ Z. They
become one Theorell–Chance step Y + L → Z + R, rapid equilibrium, in a group of its own
labelled after the others, where Y + L is the binding of a substrate (an index of
`substrates`) when one of the two is; their groups keep their other steps. No step joins
`x` then.
"""
function _eliminate_form(steps::Vector{Mechanisms.Step}, x::Int, substrates::Vector{Int})
    at_x = findall(s -> x in (s.from, s.to), steps)
    length(at_x) == 2 && all(j -> steps[j].to == x && _is_binding(steps[j]), at_x) ||
        return nothing
    binds_substrate(j) = only(steps[j].taken_up) in substrates
    entry, release = steps[sort(at_x; by = j -> !binds_substrate(j))]
    entry.from == release.from && return nothing
    keep = [!(j in at_x) for j in eachindex(steps)]
    [steps[keep]; Mechanisms.Step(entry.from, release.from, entry.taken_up,
                                  release.taken_up, true,
                                  maximum(e -> e.kinetic_group, steps) + 1)]
end

"""
The tests of `_seed_variants` that build nothing, for the candidate `c` of one base of a
seed, a plain mechanism `m`: a function of a candidate's rapid-equilibrium flags (one per
edge of `m`, some turned steady state) that is true when the candidate has no
rapid-equilibrium turnover cycle (`has_turnover_cycle`), a maximal rate both ways
(`_has_vmax`) and no chemistry in equilibrium with both sides
(`_chemistry_equilibrates_both_sides`). Errors when the base holds an isomerization.

It reads each edge once per base, so a candidate costs a few passes over arrays. A base
holds no isomerization, so each chemistry node is a merged complex or a Theorell–Chance
step. A set of edges holds a turnover cycle iff it closes a cycle of nonzero weight
(`has_nonzero_weight_cycle`). Each maximal-rate test runs on the candidate's
rapid-equilibrium edges and more, and adding edges never removes such a cycle, so a
turnover cycle of the candidate fails both tests.
"""
function _seed_candidate_screen(c::CandidateMechanism)
    m = c.mechanism.catalytic
    ks = eachindex(m.steps)
    any(k -> _is_iso(m.steps[k]), ks) &&
        error("_seed_candidate_screen: the seed base holds an isomerization; a base " *
              "has every isomerization merged")
    subs, prods = m.substrates, m.products
    from = [m.steps[k].from for k in ks]
    to = [m.steps[k].to for k in ks]
    weight = [_uptake_weight(m, k) for k in ks]
    touches(side) = BitVector([_exchanges(m, k, side) for k in ks])
    sub_step, prod_step = touches(subs), touches(prods)
    crossing = BitVector([_crosses_sides(m, k, subs, prods) for k in ks])
    complexes = _merged_complexes(c, subs)
    leaves(x, side) = BitVector([to[k] == x && _holds(m.steps[k].taken_up, side) ||
                                 from[k] == x && _holds(m.steps[k].given_off, side)
                                 for k in ks])
    exits = [(leaves(x, subs), leaves(x, prods)) for x in complexes]
    n = length(c.forms)
    parent = zeros(Int, n); offset = zeros(Int, n)
    turnover(on) = Mechanisms.has_nonzero_weight_cycle(from, to, weight, on, parent, offset)
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
The bases of the variants of the seed candidate `c`, a plain mechanism `m`, each a written
graph of the records of `c`'s forms (`_trimmed`) with its group labels numbered 1, 2, …
(`_numbered_groups`). The merged base comes first: every isomerization of `m` merged onto
its product side (`_merge_isomerization`), every step at rapid equilibrium. Then comes one
Theorell–Chance base per merged complex that `_eliminate_form` can eliminate, in the order
of the isomerizations.
"""
function _seed_bases(c::CandidateMechanism)
    m = c.mechanism.catalytic
    isos = filter(_is_iso, m.steps)
    joined = foldl(_merge_isomerization, isos; init = m.steps)
    merged = Mechanisms.Step[Mechanisms.Step(s.from, s.to, s.taken_up, s.given_off, true,
                                             s.kinetic_group) for s in joined]
    bases = [merged]
    for s in isos
        tc = _eliminate_form(merged, s.to, m.substrates)
        tc === nothing || push!(bases, tc)
    end
    [_trimmed(c.forms, _numbered_groups(b)) for b in bases]
end

"""Whether the flags `re` make the merged base `base` its seed's lumping twin: each merged
complex, forms `cx` of `base`, has two edges, both steady state and each alone in its
group."""
function _lumping_twin(base::Mechanisms.CatalyticScheme, cx::Vector{Int}, re::BitVector)
    groups = _kinetic_groups(base)
    all(cx) do f
        at = findall(e -> f == e.from || f == e.to, base.steps)
        length(at) == 2 &&
            all(k -> !re[k] && length(groups[base.steps[k].kinetic_group]) == 1, at)
    end
end

"""
The merged and Theorell–Chance variants of the seed candidate `c`, a plain mechanism of
reaction `rxn`. A variant flips to steady state an inclusion-minimal set of the groups of
one base (`_seed_bases`, `_minimal_flips`). The flip passes `_seed_candidate_screen`,
leaves every flipped group carrying flux and gives a mechanism `_enumerator_keeps` accepts.
The variants come in the order of their flip sets over the base's groups, which are in
canonical order. The function skips a merged variant that is its seed's lumping twin
(`_lumping_twin`): the unmerged form with rapid-equilibrium flanks has its family at the
same count. Errors when a base's edges do not conserve `rxn`'s atoms
(`_assert_atom_conserving`); the merge (`_merge_isomerization`, MERGE in the tests) and the
elimination (`_eliminate_form`, ELIM) conserve them.
"""
function _seed_variants(c::CandidateMechanism, rxn::Mechanisms.EnzymeReaction)
    m = c.mechanism.catalytic
    variants = CandidateMechanism[]
    for (i, b) in enumerate(_seed_bases(c))
        base = _candidate(m.metabolites, b.forms, m.substrates, m.products, b.steps)
        _assert_atom_conserving(base, rxn)
        flips = _minimal_flips(_seed_candidate_screen(base), base)
        if i == 1
            cx = _merged_complexes(base, m.substrates)
            filter!(re -> !_lumping_twin(base.mechanism.catalytic, cx, re), flips)
        end
        append!(variants, [_with_flags(base, re) for re in flips])
    end
    variants
end

# ─── Entry points ───────────────────────────────────────────────────────

"""
    init_mechanisms(reaction::EnzymeReaction) -> Vector{CandidateMechanism}

The unregulated starting mechanisms of `reaction`, as candidates of plain mechanisms. First
come the seeds that `_enumerator_keeps` accepts. Then come the merged and Theorell–Chance
variants of every seed (`_seed_variants`), in the seeds' order, each once; the flip search
has run the same test on each of them. A seed is a catalytic topology
(`_catalytic_topologies`) with a set of dead-end forms
(`_expand_substrate_product_dead_ends`), its one steady-state step a chemistry
isomerization, and its steps grouped by `_seed_groups`. Errors when the reaction has no
catalytic topology.
"""
function init_mechanisms(r::Mechanisms.EnzymeReaction)
    # A seed's steps index the reactants, substrates first (`_catalytic_topologies`).
    reactants = Symbol[[s.name for s in r.substrates]; [p.name for p in r.products]]
    ns = length(r.substrates)
    seeds = _expand_substrate_product_dead_ends(_catalytic_topologies(r), r)
    mechs = CandidateMechanism[
        _candidate(reactants, s.forms, collect(1:ns), collect((ns + 1):length(reactants)),
                   Mechanisms.Step[Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off,
                                                   e.rapid_equilibrium, g)
                                   for (e, g) in zip(s.steps, _seed_groups(s.steps))])
        for s in seeds]
    foreach(c -> _assert_atom_conserving(c, r), mechs)
    seen = Set(mechs)
    out = filter(c -> _enumerator_keeps(c.mechanism), mechs)
    for c in unique(mechs), v in _seed_variants(c, r)
        v in seen && continue
        push!(seen, v)
        push!(out, v)
    end
    out
end

"""
    seed_mechanisms(rxn, required_allo::Set{Symbol}, required_comp::Set{Symbol})
        -> Vector{CandidateMechanism}

The starting mechanisms of `rxn` when it requires the allosteric regulators `required_allo`
and the competitive inhibitors `required_comp`: the mechanisms that bind every regulator of
`required_allo` at a regulatory site and every inhibitor of `required_comp` at a dead-end
site. The regulator-binding moves grow them from the base tier of `init_mechanisms(rxn)`
(`_base_tier`) through seed nodes alone. They keep every step's rapid-equilibrium flag, so
no seed is degenerate (`_degenerate`). The dead-end move (`_expand_add_dead_end_regulator`)
always runs; the tagging and regulator moves (`_expand_to_allosteric`,
`_expand_add_allosteric_regulator`) run only when an allosteric regulator is required, so
competitive inhibitors alone keep the seeds plain, with no `L`. A seed node
(`_is_seed_node`):

1. tags no regulatory ligand `:NonequalAI` (the beam reaches those later through
   `_expand_change_allo_state`). A kinetic group can be `:NonequalAI`: the tagging move
   tags a set of groups `:NonequalAI` only where the same set tagged `:OnlyA` is not
   returned, and this search never runs `_expand_change_allo_state`, the only other route
   to such a version;
2. binds one ligand per regulatory site, which bounds the search: a child that adds a
   regulator to an existing site fails, and the search drops it before expanding it;
3. binds no optional regulator: every bound allosteric ligand is in `required_allo` and
   every bound competitive inhibitor in `required_comp`;
4. respects every ligand's declared type (`_respects_reg_type`).

Each seed appears once, in an order that does not depend on the number of workers. Errors
when no mechanism binds every required regulator. The search expands a node that
`_enumerator_returns` refuses but never returns it.
"""
function seed_mechanisms(rxn::Mechanisms.EnzymeReaction, required_allo::Set{Symbol},
                         required_comp::Set{Symbol})
    visited = Set{UInt64}()
    seeds = CandidateMechanism[]
    # Whether `c` is new. Records it, and keeps it as a seed when it binds every required
    # regulator.
    consider!(c) = begin
        h = hash(c)
        h in visited && return false
        push!(visited, h)
        issubset(required_allo, _bound_allo_regs(c.mechanism)) &&
            issubset(required_comp, _bound_comp_inhibitors(c)) && push!(seeds, c)
        true
    end
    # From the base tier, not the seeds: `_enumerator_keeps` refuses some states on a
    # degenerate seed that it keeps on the seed's flip children, and the search would never
    # give those children these states.
    base = _base_tier(unique(init_mechanisms(rxn)), rxn)
    frontier = CandidateMechanism[c for c in base if consider!(c)]
    while !isempty(frontier)
        next = CandidateMechanism[]
        for cs in _seed_node_children(frontier, rxn, required_allo, required_comp), c in cs
            consider!(c) && push!(next, c)
        end
        frontier = next
    end
    if isempty(seeds)
        list(s) = isempty(s) ? "none" : join(sort!(collect(s)), ", ")
        error("candidate_mechanisms: no mechanism binds every required regulator " *
              "(competitive inhibitors: " * list(required_comp) *
              "; allosteric regulators: " * list(required_allo) * "). A substrate or " *
              "product declared as a competitive inhibitor binds only where its complex " *
              "is not redundant; a uni-uni seed has no such site. Mark a regulator " *
              "optional with `optional_competitive_inhibitors` or " *
              "`optional_allosteric_regulators`, or remove its declaration")
    end
    filter!(_enumerator_returns, seeds)
end

"""The children of each candidate of `frontier` that are seed nodes (`_is_seed_node`), in
frontier order. Generating one candidate's children is pure and independent of the
others', so the candidates are distributed over the workers."""
function _seed_node_children(
    frontier::Vector{CandidateMechanism}, rxn::Mechanisms.EnzymeReaction,
    required_allo::Set{Symbol}, required_comp::Set{Symbol}
)::Vector{Vector{CandidateMechanism}}
    # The return type states what inference cannot read through `pmap`.
    pmap(frontier) do c
        children = CandidateMechanism[]
        isempty(required_allo) || append!(children, _is_allosteric(c.mechanism) ?
            _expand_add_allosteric_regulator(c, rxn) : _expand_to_allosteric(c, rxn))
        append!(children, _expand_add_dead_end_regulator(c, rxn))
        filter!(child -> _is_seed_node(child, rxn, required_allo, required_comp),
                children)
    end
end

function _is_seed_node(c::CandidateMechanism, rxn::Mechanisms.EnzymeReaction,
                       required_allo::Set{Symbol}, required_comp::Set{Symbol})
    m = c.mechanism
    any(site -> length(site.ligands) != 1 || only(site.allosteric_states) == :NonequalAI,
        m.allosteric.regulatory_sites) && return false
    issubset(_bound_allo_regs(m), required_allo) || return false
    issubset(_bound_comp_inhibitors(c), required_comp) || return false
    _respects_reg_type(m, rxn) || return false
    true
end

function candidate_mechanisms(reaction::Mechanisms.EnzymeReaction;
                              max_param_count::Int = typemax(Int),
                              eq_complexity_filter::Int = typemax(Int),
                              optional_allosteric_regulators::Vector{Symbol} = Symbol[],
                              optional_competitive_inhibitors::Vector{Symbol} = Symbol[])
    _check_decodable_names(reaction)
    required_allo, required_comp = _required_regulators(
        reaction, optional_allosteric_regulators, optional_competitive_inhibitors)
    seeds = (isempty(required_allo) && isempty(required_comp)) ?
        _base_tier(unique(init_mechanisms(reaction)), reaction) :
        seed_mechanisms(reaction, required_allo, required_comp)
    _within_size_limits(seeds, max_param_count, eq_complexity_filter)
end

"""
    _required_regulators(rxn, optional_allosteric_regulators,
                         optional_competitive_inhibitors)
        -> (required_allo::Set{Symbol}, required_comp::Set{Symbol})

The regulators the starting mechanisms must bind: every allosteric regulator and every
competitive inhibitor declared in `rxn`, minus the names the caller marked optional. Both
sets empty means the search starts from the base tier of `init_mechanisms` (`_base_tier`),
which binds no regulator.
"""
function _required_regulators(rxn::Mechanisms.EnzymeReaction,
                              optional_allosteric_regulators::Vector{Symbol},
                              optional_competitive_inhibitors::Vector{Symbol})
    required(kind, optional) = setdiff(
        Set{Symbol}(reg.name for reg in rxn.regulators if reg.kind === kind), optional)
    (required(:allosteric, optional_allosteric_regulators),
     required(:competitive, optional_competitive_inhibitors))
end

"""
The plain mechanisms the search starts from: every seed of `mechs` that is not degenerate
(`_degenerate`), and in place of each one that is, its flip children (`_expand_re_to_ss`)
that are not. The search fits them first when the reaction requires no regulator, and
`seed_mechanisms` binds the required regulators to them when it does. A degenerate seed's
law ignores a substrate or never saturates, so it is not worth a fit, but it is the only
parent of mechanisms that are not degenerate (such as the ping-pong seeds' children). Only
a flip can cure it: the other moves keep every catalytic step's rapid-equilibrium flag, and
the steps they add bind regulators, which neither degeneracy condition reads. Errors when a
flip child does not conserve atoms (`_assert_atom_conserving`).
"""
function _base_tier(mechs::Vector{CandidateMechanism}, rxn::Mechanisms.EnzymeReaction)
    degenerate = _degenerate.(mechs)
    base = mechs[.!degenerate]
    for c in mechs[degenerate]
        for child in _expand_re_to_ss(c)
            _assert_atom_conserving(child, rxn)
            _degenerate(child) || push!(base, child)
        end
    end
    unique(base)
end

"""
Error unless the names of the metabolites `rxn` binds to an enzyme form spell every
string one way. The enumerator names each form it builds by joining the names of the
metabolites the form binds without a separator (`EAB`, `_form_name`), so two forms could
share a name (the `@enzyme_reaction` docstring gives an example). Each substrate and
product is its own label; a competitive inhibitor that is also a reactant binds as its
modifier copy, named `Xmod`; an allosteric regulator binds no form, so its name takes no
part. The Sardinas–Patterson test decides it, and the error spells both splits.
"""
function _check_decodable_names(rxn::Mechanisms.EnzymeReaction)
    reactants = Symbol[[r.name for r in rxn.substrates]; [r.name for r in rxn.products]]
    # Each label and the metabolite it stands for, as the error names it.
    written = Dict{String, String}()
    copy_label(n) = n in reactants ? "$(n)mod" => "$(n)mod ($n's modifier copy)" :
                                     String(n) => String(n)
    for (label, text) in [[String(n) => String(n) for n in reactants];
                          [copy_label(r.name) for r in rxn.regulators
                           if r.kind === :competitive]]
        other = get!(written, label, text)
        other == text || _spelled_twice_error([other], [text], label)
    end
    words = sort!(collect(keys(written)))
    tail(w, prefix) = w[ncodeunits(prefix) + 1:end]
    # Each entry holds a dangling suffix `d` and two splits: join(long) == join(short) * d.
    queue = [(tail(v, u), [v], [u]) for u in words for v in words
             if u != v && startswith(v, u)]
    seen = Set{String}()
    while !isempty(queue)
        d, long, short = popfirst!(queue)
        d in seen && continue
        push!(seen, d)
        haskey(written, d) && _spelled_twice_error(getindex.(Ref(written), [short; d]),
                                                   getindex.(Ref(written), long),
                                                   join(long))
        for w in words
            startswith(w, d) && push!(queue, (tail(w, d), [short; w], long))
            startswith(d, w) && push!(queue, (tail(d, w), long, [short; w]))
        end
    end
end

"""Error that the label splits `a` and `b`, each a list of the metabolites the labels stand
for, both spell `spelled`."""
function _spelled_twice_error(a::Vector{String}, b::Vector{String}, spelled::String)
    error("candidate_mechanisms: $(join(a, " + ")) and $(join(b, " + ")) both spell " *
          "$spelled: enzyme forms and their constants join the names of the metabolites " *
          "they bind without a separator, so two forms would share a name; rename a " *
          "metabolite")
end
