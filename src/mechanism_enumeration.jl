# ABOUTME: Mechanism enumeration by incremental parameter count growth
# ABOUTME: Provides init_mechanisms, expand_mechanisms building blocks


# ─── Catalytic topologies ─────────────────────────────────────────────

"""
Build a `Species` on conformation `:E` from sorted bound substrate /
product names plus the covalent `residual`. Names map to `Substrate` /
`Product` structs; the bound list is sorted by name (Species inner
constructor enforces this).
"""
function _make_species(
    bound_subs::Vector{Symbol},
    bound_prods::Vector{Symbol},
    residual::Residual,
)
    mets = Metabolite[Substrate.(bound_subs)...,
                      Product.(bound_prods)...]
    Species(mets, :E, residual)
end

"""
Covalent residual at an enzyme form: `added` = consumed substrates not
currently bound, `subtracted` = released products plus currently-bound
products. Reduced to `Residual()` exactly when the added/subtracted atom
multisets cancel (no covalent residue remains).
"""
function _residual_for(
    consumed::Vector{Symbol},
    on_subs::Vector{Symbol},
    released::Vector{Symbol},
    on_prods::Vector{Symbol},
    sub_atoms::Dict{Symbol,Dict{Symbol,Int}},
    prod_atoms::Dict{Symbol,Dict{Symbol,Int}},
)
    add_names = setdiff(consumed, on_subs)
    sub_names = vcat(released, on_prods)
    add_at = reduce(_add_atoms, [sub_atoms[s] for s in add_names];
                    init=Dict{Symbol,Int}())
    sub_at = reduce(_add_atoms, [prod_atoms[p] for p in sub_names];
                    init=Dict{Symbol,Int}())
    _nonzero_atoms(add_at) == _nonzero_atoms(sub_at) && return Residual()
    Residual(Substrate.(add_names), Product.(sub_names))
end

"""Extract atom counts as Dict{Symbol,Int} for a metabolite."""
function _atoms_dict(
    reaction::EnzymeReaction,
    met::Symbol,
)
    result = Dict{Symbol,Int}()
    for ra in reactants(reaction)
        m = metabolite(ra)
        (m isa Substrate || m isa Product) || continue
        name(m) == met || continue
        for (a, c) in atoms(ra)
            result[a] = get(result, a, 0) + c
        end
        return result
    end
    result
end

"""Check if product atoms are a subset of accumulated atoms."""
function _can_pingpong(
    accumulated::Dict{Symbol,Int},
    prod_atoms::Dict{Symbol,Int},
)
    for (a, c) in prod_atoms
        get(accumulated, a, 0) < c && return false
    end
    true
end

"""Subtract atom counts: accumulated minus product atoms."""
function _subtract_atoms(
    accumulated::Dict{Symbol,Int},
    prod_atoms::Dict{Symbol,Int},
)
    result = copy(accumulated)
    for (a, c) in prod_atoms
        result[a] -= c
        result[a] == 0 && delete!(result, a)
    end
    result
end

"""Add atom counts: accumulated plus substrate atoms."""
function _add_atoms(
    accumulated::Dict{Symbol,Int},
    sub_atoms::Dict{Symbol,Int},
)
    result = copy(accumulated)
    for (a, c) in sub_atoms
        result[a] = get(result, a, 0) + c
    end
    result
end

# ─── Atom-conservation validation ────────────────────────────

"""Add `sign * d` into the signed accumulator `acc` in place."""
function _accumulate_atoms!(acc::Dict{Symbol,Int}, d::Dict{Symbol,Int}, sign::Int)
    for (a, c) in d
        acc[a] = get(acc, a, 0) + sign * c
    end
    acc
end

"""Drop zero entries from a signed atom dict."""
_nonzero_atoms(d::Dict{Symbol,Int}) = filter(kv -> kv.second != 0, d)

"""
Net atom multiset carried by a `Species`: atoms of its bound metabolites
plus atoms of `residual.added` minus atoms of `residual.subtracted`, read
from the reaction's per-metabolite inventory via `_atoms_dict`.
"""
function _species_atoms(reaction::EnzymeReaction, sp::Species)
    acc = Dict{Symbol,Int}()
    for m in bound(sp)
        _accumulate_atoms!(acc, _atoms_dict(reaction, name(m)), 1)
    end
    for a in added(residual(sp))
        _accumulate_atoms!(acc, _atoms_dict(reaction, name(a)), 1)
    end
    for p in subtracted(residual(sp))
        _accumulate_atoms!(acc, _atoms_dict(reaction, name(p)), -1)
    end
    _nonzero_atoms(acc)
end

"""
Assert one `Step` conserves atoms: the atoms of `from_species` plus those of
the consumed metabolites must equal the atoms of `to_species` plus those of the
released metabolites (an iso step leaves the atom multiset unchanged). Errors
naming the offending step.
"""
function _assert_step_atom_conserving(reaction::EnzymeReaction, s::Step)
    diff = Dict{Symbol,Int}()
    _accumulate_atoms!(diff, _species_atoms(reaction, to_species(s)), 1)
    _accumulate_atoms!(diff, _species_atoms(reaction, from_species(s)), -1)
    for m in consumed(s); _accumulate_atoms!(diff, _atoms_dict(reaction, name(m)), -1); end
    for m in released(s); _accumulate_atoms!(diff, _atoms_dict(reaction, name(m)), 1); end
    diff = _nonzero_atoms(diff)
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
(`_flux_carrying_groups`). No kinetic group binds a competitive inhibitor redundantly
(`_redundant_copy_groups`). A parent must obey both because a flip tests only the groups
it flips. The split, the dead-end move, and `_expand_change_allo_state` filter their
children; the other moves preserve both rules.
"""
function _assert_emission_rules(m::Union{Mechanism, AllostericMechanism})
    label(g) = join((join(_forward_sides(s), " → ") for s in steps(m)[g]), ", ")
    flux = _flux_carrying_groups(m)
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

"""Generate all combinations of `k` elements from `arr`."""
function _combinations(arr, k)
    n = length(arr)
    k == 0 && return [eltype(arr)[]]
    k == 1 && return [[x] for x in arr]
    k == n && return [collect(arr)]
    result = Vector{Vector{eltype(arr)}}()
    for i in 1:n
        for rest in _combinations(arr[i+1:end], k - 1)
            push!(result, [arr[i]; rest])
        end
    end
    result
end

"""
Release products one at a time after a multi-product isomerization, then
continue backtracking. `pingpong_intermediate` is the ping-pong control flag
threaded to `backtrack!`; each released form carries the covalent
residual derived from the consumed/released history via `_residual_for`.
`residual_atoms` is the covalent residue remaining after the whole
`prod_subset` is released (passed to `backtrack!` as the enzyme's atoms).
"""
function _release_products!(
    all_paths, backtrack!,
    iso_species::Species,
    residual_atoms::Dict{Symbol,Int},
    consumed_subs::Vector{Symbol},
    released_prods::Vector{Symbol},
    prod_subset::Vector{Symbol},
    sub_atoms::Dict{Symbol,Dict{Symbol,Int}},
    prod_atoms::Dict{Symbol,Dict{Symbol,Int}},
    pingpong_intermediate::Bool,
    steps::Vector{Step},
)
    # Generate all release orderings of products
    function _release_recurse!(
        cur::Species,
        unreleased::Vector{Symbol},
        rel_so_far::Vector{Symbol},
    )
        if isempty(unreleased)
            # All products released, continue
            backtrack!(
                cur, residual_atoms,
                consumed_subs, rel_so_far,
                Symbol[], Symbol[],
                pingpong_intermediate, false, steps
            )
            return
        end
        for p in copy(unreleased)
            new_unreleased = filter(!=(p), unreleased)
            new_species = _make_species(
                Symbol[], new_unreleased,
                _residual_for(consumed_subs, Symbol[],
                              [rel_so_far; p], new_unreleased,
                              sub_atoms, prod_atoms))
            rel_step = Step(
                cur, new_species, Metabolite[], Metabolite[Product(p)], true)
            push!(steps, rel_step)
            _release_recurse!(
                new_species, new_unreleased, [rel_so_far; p])
            pop!(steps)
        end
    end

    _release_recurse!(
        iso_species, collect(prod_subset), copy(released_prods))
end

"""
Substrate bound-metabolite names in route (path) order.
"""
_binding_order(path::Vector{Step}) =
    Symbol[name(bound_metabolite(s)) for s in path if bound_metabolite(s) isa Substrate]

"""
Product bound-metabolite names in route (path) order.
"""
_release_order(path::Vector{Step}) =
    Symbol[name(bound_metabolite(s)) for s in path if bound_metabolite(s) isa Product]

"""
True iff `order` is a linearization of weak ordering `wo` (a vector of
levels): every metabolite of `wo` appears exactly once in `order`, and the
level index along `order` is non-decreasing (earlier levels strictly before
later levels; any order within a level).
"""
function _linearizes(order::Vector{Symbol}, wo::Vector{Vector{Symbol}})
    level = Dict{Symbol,Int}()
    for (i, lvl) in enumerate(wo), m in lvl
        level[m] = i
    end
    length(order) == length(level) || return false
    prev = 0
    for m in order
        haskey(level, m) || return false
        level[m] < prev && return false
        prev = level[m]
    end
    true
end

"""
    _catalytic_topologies(reaction) -> Vector{Vector{Step}}

Build catalytic cycle topologies by constructive backtracking.
Each topology is a set of steps forming one or more complete
catalytic cycles (E -> ... -> E).
"""
function _catalytic_topologies(
    reaction::EnzymeReaction,
)
    sub_names = Symbol[name(s) for s in substrates(reaction)]
    prod_names = Symbol[name(p) for p in products(reaction)]

    # Precompute atom dicts for each metabolite
    sub_atoms = Dict(
        s => _atoms_dict(reaction, s) for s in sub_names
    )
    prod_atoms = Dict(
        p => _atoms_dict(reaction, p) for p in prod_names
    )

    # C5: max simultaneously bound metabolites
    max_bound = max(length(sub_names), length(prod_names))

    # Collect all complete catalytic paths as Step lists
    all_paths = Vector{Vector{Step}}()

    # Backtracking state:
    # - cur_species: current enzyme form as a `Species`
    # - acc_atoms: atoms currently on the enzyme
    # - consumed_subs: substrates consumed so far (history)
    # - released_prods: products released so far (history)
    # - on_enzyme_subs: substrates currently bound
    # - on_enzyme_prods: products currently bound
    #     (post-final-isomerize)
    # - pingpong_intermediate: enzyme is in a ping-pong
    #     covalent-intermediate state (carries a residual)
    # - post_final: in product-release phase after final
    #     isomerization
    # - steps: path of Step accumulated so far
    function backtrack!(
        cur_species::Species,
        acc_atoms::Dict{Symbol,Int},
        consumed_subs::Vector{Symbol},
        released_prods::Vector{Symbol},
        on_enzyme_subs::Vector{Symbol},
        on_enzyme_prods::Vector{Symbol},
        # pingpong_intermediate: true when the enzyme is in a
        # ping-pong covalent-intermediate state. Selects the
        # bind-only / iso branches below (the covalent residue
        # itself is stored on the form as a Residual).
        pingpong_intermediate::Bool,
        post_final::Bool,
        steps::Vector{Step},
    )
        # Check for complete cycle
        if conformation(cur_species) === :E &&
                isempty(bound(cur_species)) && !isempty(steps)
            if Set(consumed_subs) == Set(sub_names) &&
                    Set(released_prods) == Set(prod_names)
                push!(all_paths, copy(steps))
                return
            end
        end

        remaining_subs = [
            s for s in sub_names if s ∉ consumed_subs
        ]
        remaining_prods = [
            p for p in prod_names if p ∉ released_prods
        ]

        if post_final
            # Release any currently bound product
            for p in copy(on_enzyme_prods)
                new_on_prods = filter(!=(p), on_enzyme_prods)
                new_released = [released_prods; p]
                new_species = _make_species(
                    Symbol[], new_on_prods,
                    _residual_for(consumed_subs, Symbol[],
                                  new_released, new_on_prods,
                                  sub_atoms, prod_atoms))
                step = Step(
                    cur_species, new_species,
                    Metabolite[], Metabolite[Product(p)], true)
                push!(steps, step)
                backtrack!(
                    new_species,
                    _subtract_atoms(
                        acc_atoms, prod_atoms[p]
                    ),
                    consumed_subs, new_released,
                    Symbol[], new_on_prods,
                    false, !isempty(new_on_prods),
                    steps
                )
                pop!(steps)
            end
            return
        end

        if isempty(on_enzyme_subs) && !pingpong_intermediate
            # Free enzyme: bind any remaining substrate
            for s in remaining_subs
                new_on = [on_enzyme_subs; s]
                new_consumed = [consumed_subs; s]
                new_species = _make_species(
                    new_on, Symbol[],
                    _residual_for(new_consumed, new_on,
                                  released_prods, Symbol[],
                                  sub_atoms, prod_atoms))
                step = Step(
                    cur_species, new_species,
                    Metabolite[Substrate(s)], Metabolite[], true)
                push!(steps, step)
                backtrack!(
                    new_species,
                    _add_atoms(acc_atoms, sub_atoms[s]),
                    new_consumed, released_prods,
                    new_on, Symbol[],
                    false, false, steps
                )
                pop!(steps)
            end
        elseif !isempty(on_enzyme_subs) && !pingpong_intermediate
            # Substrates bound, no residual
            # Option 1: bind another substrate (C5)
            if length(on_enzyme_subs) < max_bound
                for s in remaining_subs
                    new_on = [on_enzyme_subs; s]
                    new_consumed = [consumed_subs; s]
                    new_species = _make_species(
                        new_on, Symbol[],
                        _residual_for(new_consumed, new_on,
                                      released_prods, Symbol[],
                                      sub_atoms, prod_atoms))
                    step = Step(
                        cur_species, new_species,
                        Metabolite[Substrate(s)], Metabolite[], true)
                    push!(steps, step)
                    backtrack!(
                        new_species,
                        _add_atoms(
                            acc_atoms, sub_atoms[s]),
                        new_consumed, released_prods,
                        new_on, Symbol[],
                        false, false, steps
                    )
                    pop!(steps)
                end
            end
            # Option 2: ping-pong isomerize (C9)
            if !isempty(remaining_subs)
                for k in 1:length(remaining_prods)
                    for prod_subset in _combinations(
                        remaining_prods, k)
                        need = reduce(
                            _add_atoms,
                            [prod_atoms[p]
                             for p in prod_subset];
                            init=Dict{Symbol,Int}()
                        )
                        _can_pingpong(
                            acc_atoms, need
                        ) || continue
                        residual = _subtract_atoms(
                            acc_atoms, need
                        )
                        # Admissible-residual rule: a ping-pong
                        # continuation must form a genuine covalent
                        # residue. An empty residue means the enzyme
                        # returns to apo E mid-cycle while substrates
                        # remain unbound, splitting the reaction into
                        # disconnected half-cycles — not a valid
                        # mechanism. The final isomerization (Option 3,
                        # all substrates consumed) handles the
                        # legitimate return to apo E.
                        isempty(residual) && continue
                        n_prods_eff = k + 1
                        # C6: iso size limit
                        length(on_enzyme_subs) > 3 &&
                            continue
                        n_prods_eff > 3 && continue
                        # C8: product-only iso form
                        iso_species = _make_species(
                            Symbol[],
                            collect(prod_subset),
                            _residual_for(
                                consumed_subs, Symbol[],
                                released_prods,
                                collect(prod_subset),
                                sub_atoms, prod_atoms))
                        step = Step(
                            cur_species, iso_species,
                            Metabolite[], Metabolite[], true)
                        push!(steps, step)
                        # Release products one at a time. This
                        # ping-pong continuation carries a genuine
                        # covalent residual (the empty-residue case is
                        # filtered above), so the control bool is true.
                        _release_products!(
                            all_paths, backtrack!,
                            iso_species, residual,
                            consumed_subs, released_prods,
                            prod_subset, sub_atoms, prod_atoms,
                            true, steps
                        )
                        pop!(steps)
                    end
                end
            end
            # Option 3: final isomerize (all subs bound)
            if isempty(remaining_subs)
                all_prod_atoms = reduce(
                    _add_atoms,
                    [prod_atoms[p]
                     for p in remaining_prods];
                    init=Dict{Symbol,Int}()
                )
                # C6: iso size limit
                n_subs_react = length(on_enzyme_subs)
                n_prods_eff = length(remaining_prods)
                if n_subs_react <= 3 &&
                        n_prods_eff <= 3 &&
                        _can_pingpong(
                            acc_atoms, all_prod_atoms)
                    new_species = _make_species(
                        Symbol[],
                        copy(remaining_prods),
                        _residual_for(
                            consumed_subs, Symbol[],
                            released_prods,
                            copy(remaining_prods),
                            sub_atoms, prod_atoms))
                    step = Step(
                        cur_species, new_species,
                        Metabolite[], Metabolite[], true)
                    push!(steps, step)
                    backtrack!(
                        new_species, acc_atoms,
                        consumed_subs, released_prods,
                        Symbol[],
                        copy(remaining_prods),
                        false, true, steps
                    )
                    pop!(steps)
                end
            end
        elseif isempty(on_enzyme_subs) && pingpong_intermediate
            # C7: residual-bearing form with no subs — only bind, no iso
            for s in remaining_subs
                new_on = [s]
                new_consumed = [consumed_subs; s]
                new_species = _make_species(
                    new_on, Symbol[],
                    _residual_for(new_consumed, new_on,
                                  released_prods, Symbol[],
                                  sub_atoms, prod_atoms))
                step = Step(
                    cur_species, new_species,
                    Metabolite[Substrate(s)], Metabolite[], true)
                push!(steps, step)
                backtrack!(
                    new_species,
                    _add_atoms(acc_atoms, sub_atoms[s]),
                    new_consumed, released_prods,
                    new_on, Symbol[],
                    true, false, steps
                )
                pop!(steps)
            end
        elseif !isempty(on_enzyme_subs) && pingpong_intermediate
            # Residual + substrates bound
            # Option 1: bind another substrate (C5)
            if length(on_enzyme_subs) < max_bound
                for s in remaining_subs
                    new_on = [on_enzyme_subs; s]
                    new_consumed = [consumed_subs; s]
                    new_species = _make_species(
                        new_on, Symbol[],
                        _residual_for(new_consumed, new_on,
                                      released_prods, Symbol[],
                                      sub_atoms, prod_atoms))
                    step = Step(
                        cur_species, new_species,
                        Metabolite[Substrate(s)], Metabolite[], true)
                    push!(steps, step)
                    backtrack!(
                        new_species,
                        _add_atoms(
                            acc_atoms, sub_atoms[s]),
                        new_consumed, released_prods,
                        new_on, Symbol[],
                        true, false, steps
                    )
                    pop!(steps)
                end
            end
            # Option 2: isomerize to release products
            # (C9: multi-product release)
            for k in 1:length(remaining_prods)
                for prod_subset in _combinations(
                    remaining_prods, k)
                    need = reduce(
                        _add_atoms,
                        [prod_atoms[p]
                         for p in prod_subset];
                        init=Dict{Symbol,Int}()
                    )
                    _can_pingpong(
                        acc_atoms, need
                    ) || continue
                    residual_atoms = _subtract_atoms(
                        acc_atoms, need
                    )
                    n_prods_eff = k + (
                        isempty(residual_atoms) ? 0 : 1)
                    # C6: iso size limit
                    length(on_enzyme_subs) > 3 &&
                        continue
                    n_prods_eff > 3 && continue

                    has_more = !isempty(residual_atoms)
                    is_final = !has_more &&
                        isempty(remaining_subs) &&
                        k == length(remaining_prods)

                    if is_final
                        # Final iso: release all
                        # remaining products
                        new_species = _make_species(
                            Symbol[],
                            copy(remaining_prods),
                            _residual_for(
                                consumed_subs, Symbol[],
                                released_prods,
                                copy(remaining_prods),
                                sub_atoms, prod_atoms))
                        step = Step(
                            cur_species, new_species,
                            Metabolite[], Metabolite[], true)
                        push!(steps, step)
                        backtrack!(
                            new_species, acc_atoms,
                            consumed_subs,
                            released_prods,
                            Symbol[],
                            copy(remaining_prods),
                            false, true, steps
                        )
                        pop!(steps)
                    else
                        # Admissible-residual rule (mirrors the no-residual
                        # ping-pong branch): a non-final ping-pong iso must
                        # leave a genuine covalent residue. Without one the
                        # enzyme returns to apo E mid-cycle while substrates
                        # remain — a disconnected half-cycle.
                        has_more || continue
                        # C8: product-only iso form
                        iso_species = _make_species(
                            Symbol[],
                            collect(prod_subset),
                            _residual_for(
                                consumed_subs, Symbol[],
                                released_prods,
                                collect(prod_subset),
                                sub_atoms, prod_atoms))
                        step = Step(
                            cur_species, iso_species,
                            Metabolite[], Metabolite[], true)
                        push!(steps, step)
                        _release_products!(
                            all_paths, backtrack!,
                            iso_species, residual_atoms,
                            consumed_subs,
                            released_prods,
                            prod_subset, sub_atoms, prod_atoms,
                            has_more, steps
                        )
                        pop!(steps)
                    end
                end
            end
        end
    end

    # Start from free enzyme (:E, no bound metabolites).
    free_E = Species(Metabolite[], :E)
    backtrack!(
        free_E, Dict{Symbol,Int}(), Symbol[], Symbol[],
        Symbol[], Symbol[], false, false, Step[]
    )

    isempty(all_paths) && return Vector{Step}[]

    # Deduplicate paths by their structural step content.
    # `Step` equality / hash use canonical direction, so this
    # correctly identifies equal step multi-sets across paths.
    unique_paths = Vector{Vector{Step}}()
    seen_path_keys = Set{Set{Step}}()
    for path in all_paths
        key = Set(path)
        key ∈ seen_path_keys && continue
        push!(seen_path_keys, key)
        push!(unique_paths, path)
    end

    # --- Group paths by isomerization pattern ---
    _iso_pattern(path) = Set(s for s in path if is_iso(s))

    iso_groups = Dict{Set{Step}, Vector{Vector{Step}}}()
    for path in unique_paths
        pat = _iso_pattern(path)
        push!(get!(iso_groups, pat, Vector{Vector{Step}}()),
              path)
    end

    # --- Enumerate weak orderings within each group ---
    function _weak_orderings(items::Vector{T}) where T
        n = length(items)
        n == 0 && return [Vector{Vector{T}}()]
        n == 1 && return [[items]]
        orderings = Vector{Vector{Vector{T}}}()
        _wo_recurse!(orderings, Vector{Vector{T}}(), items)
        orderings
    end

    function _wo_recurse!(
        orderings, prefix, remaining::Vector{T},
    ) where T
        if isempty(remaining)
            push!(orderings, copy(prefix))
            return
        end
        for mask in 1:(2^length(remaining) - 1)
            level = T[]
            rest = T[]
            for (i, item) in enumerate(remaining)
                if (mask >> (i - 1)) & 1 == 1
                    push!(level, item)
                else
                    push!(rest, item)
                end
            end
            push!(prefix, sort(level))
            _wo_recurse!(orderings, prefix, rest)
            pop!(prefix)
        end
    end

    # --- Build topologies: union whole paths consistent with each
    # (substrate weak-ordering, product weak-ordering). Unioning complete
    # paths (rather than cherry-picking steps) keeps every form connected to
    # the catalytic complex — paths consistent with one weak ordering never
    # carry contradictory binding orders, so no dangling single-metabolite
    # forms arise. Binding history is read from the path, so ping-pong (where
    # a consumed substrate leaves the bound set) is handled correctly.
    #
    # Iterate iso_groups deterministically (smaller iso-step counts first,
    # then by sorted iso-step names) so topology output order is stable —
    # `Set{Step}` hashing is not value-stable.
    sorted_iso_pats = sort(collect(keys(iso_groups));
        by = pat -> (
            length(pat),
            sort([
                (string(name(from_species(s))),
                 string(name(to_species(s))))
                for s in pat])))

    result = Vector{Step}[]
    for iso_pat in sorted_iso_pats
        group_paths = iso_groups[iso_pat]

        sub_binding_mets = Set{Symbol}()
        prod_binding_mets = Set{Symbol}()
        for path in group_paths, step in path
            bm = bound_metabolite(step)
            bm === nothing && continue
            if bm isa Substrate
                push!(sub_binding_mets, name(bm))
            elseif bm isa Product
                push!(prod_binding_mets, name(bm))
            end
        end

        sub_orderings = _weak_orderings(
            sort(collect(sub_binding_mets)))
        prod_orderings = _weak_orderings(
            sort(collect(prod_binding_mets)))

        seen_topos = Set{Set{Step}}()
        for sub_ord in sub_orderings, prod_ord in prod_orderings
            topo_keys = Set{Step}()
            matched = false
            for path in group_paths
                _linearizes(_binding_order(path), sub_ord) || continue
                _linearizes(_release_order(path), prod_ord) || continue
                union!(topo_keys, path)
                matched = true
            end
            matched || continue
            topo_keys ∈ seen_topos && continue
            push!(seen_topos, topo_keys)

            steps = sort(collect(topo_keys); by=s -> (
                is_iso(s) ? 1 : 0,
                string(name(from_species(s))),
                string(name(to_species(s))),
            ))

            # The first iso step is the (single) SS step; every other step is
            # RE. Rebuild each Step with that tag (Step is immutable; direction
            # is unaffected by is_equilibrium).
            iso_idx = findfirst(is_iso, steps)
            push!(result, Step[_with_equilibrium(s, i != iso_idx)
                               for (i, s) in enumerate(steps)])
        end
    end
    result
end

# ─── Dead-End Helpers ────────────────────────────────────────

"""
    _substrate_product_dead_end_opportunities(
        form_sp, bound, cat_forms, sub_names, prod_names,
        add_metabolite)

Find (form, metabolite) dead-end opportunities for
substrates and products. `form_sp` maps each catalytic
form name to its `Species`; `add_metabolite(species, met)`
returns the `Species` with `met` added to its bound list,
used to render the candidate dead-end form name. A dead-end
is valid when:
- The form doesn't already bind all substrates or all
  products
- The metabolite isn't already bound at the form
- The resulting form isn't a catalytic form
- The result binds at least one substrate AND at least
  one product (mixed binding required)
- The result doesn't have all substrates or all products
"""
function _substrate_product_dead_end_opportunities(
    form_sp::Dict{Symbol, Species},
    bound::Dict{Symbol, Set{Symbol}},
    cat_forms::Set{Symbol},
    sub_names::Set{Symbol},
    prod_names::Set{Symbol},
    add_metabolite,
)
    all_mets = union(sub_names, prod_names)
    opportunities = Tuple{Symbol, Symbol}[]
    for f in sort(collect(cat_forms))
        haskey(bound, f) || continue
        fb = bound[f]
        fb_subs = intersect(fb, sub_names)
        fb_prods = intersect(fb, prod_names)
        # Eligible: neither all subs nor all prods
        (fb_subs == sub_names ||
            fb_prods == prod_names) && continue
        for m in sort(collect(all_mets))
            m in fb && continue
            de_name = name(add_metabolite(form_sp[f], m))
            de_name in cat_forms && continue
            new_bound = union(fb, Set([m]))
            new_subs = intersect(
                new_bound, sub_names)
            new_prods = intersect(
                new_bound, prod_names)
            # Must bind at least one of each type
            if isempty(new_subs) || isempty(new_prods)
                continue
            end
            # Must not bind all of either type
            if new_subs == sub_names ||
                    new_prods == prod_names
                continue
            end
            push!(opportunities, (f, m))
        end
    end
    opportunities
end

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
    subs = sort(collect(sub_names))
    prods = sort(collect(prod_names))
    edges = [(s, p) for s in subs for p in prods]
    n = length(edges)
    result = Set{Tuple{Symbol,Symbol}}[]
    for mask in 1:(1 << n) - 1
        pat = Set{Tuple{Symbol,Symbol}}()
        for j in 1:n
            if (mask >> (j - 1)) & 1 == 1
                push!(pat, edges[j])
            end
        end
        all(s -> any(
            p -> (s, p) in pat, prods), subs) ||
            continue
        all(p -> any(
            s -> (s, p) in pat, subs), prods) ||
            continue
        push!(result, pat)
    end
    result
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
    subs = sort(collect(sub_names))
    prods = sort(collect(prod_names))
    inhs = sort(existing_inhibitors)
    n_s = length(subs)
    n_p = length(prods)
    n_i = length(inhs)

    result = Tuple{
        Set{Symbol}, Set{Symbol}, Set{Symbol}}[]
    for s_mask in 1:(1 << n_s) - 1
        comp_subs = Set{Symbol}()
        for j in 1:n_s
            if (s_mask >> (j - 1)) & 1 == 1
                push!(comp_subs, subs[j])
            end
        end
        for p_mask in 1:(1 << n_p) - 1
            comp_prods = Set{Symbol}()
            for j in 1:n_p
                if (p_mask >> (j - 1)) & 1 == 1
                    push!(comp_prods, prods[j])
                end
            end
            for i_mask in 0:(1 << n_i) - 1
                comp_inhs = Set{Symbol}()
                for j in 1:n_i
                    if (i_mask >> (j - 1)) & 1 == 1
                        push!(comp_inhs, inhs[j])
                    end
                end
                push!(result, (
                    comp_subs, comp_prods, comp_inhs))
            end
        end
    end
    result
end

"""
    _expand_substrate_product_dead_ends(topos, reaction)
        -> Vector{Tuple{Vector{Step}, Vector{Int}}}

For each catalytic topology, enumerate substrate/product dead-end
form combinations and return each resulting mechanism as a flat
`Vector{Step}` paired with a parallel `Vector{Int}` of kinetic-group
ids. A dead-end form is created when a substrate or product binds to
a catalytic form where it doesn't normally bind, subject to:
- The resulting form is not already a catalytic form
- The resulting form binds at least one substrate AND
  at least one product (mixed binding required)
- The resulting form doesn't have all substrates or
  all products
- A declared `shared_catalytic_site` pair is never both bound in the
  resulting form
"""
function _expand_substrate_product_dead_ends(
    topos::Vector{Vector{Step}},
    reaction::EnzymeReaction,
)
    sub_names = Set(name(s) for s in substrates(reaction))
    prod_names = Set(name(p) for p in products(reaction))
    all_mets = union(sub_names, prod_names)

    # Competition patterns depend only on the reaction,
    # not the topology — compute once.
    patterns = _competition_patterns(
        sub_names, prod_names)

    # A declared shared catalytic site forbids its (substrate, product) pair
    # from co-occupying the catalytic site, so keep only competition patterns
    # whose forbidden-edge set contains every declared pair. The complete
    # bipartite pattern contains all edges and always survives, so the list is
    # never empty.
    shared = shared_catalytic_site(reaction)
    if !isempty(shared)
        patterns = filter(
            pat -> all(edge -> edge in pat, shared), patterns)
    end

    _role(m::Symbol) = m in sub_names ? Substrate(m) : Product(m)
    _add(sp::Species, m::Symbol) = Species(
        Metabolite[bound(sp)..., _role(m)],
        conformation(sp), residual(sp))

    result = Tuple{Vector{Step}, Vector{Int}}[]
    for topo in topos
        # Form name → Species and → bound-metabolite-name set.
        form_sp = Dict{Symbol, Species}()
        for s in topo
            form_sp[name(from_species(s))] = from_species(s)
            form_sp[name(to_species(s))] = to_species(s)
        end
        boundmap = Dict{Symbol, Set{Symbol}}(
            f => Set(name(b) for b in bound(sp))
            for (f, sp) in form_sp)
        cat_forms = Set(keys(form_sp))

        de_opportunities =
            _substrate_product_dead_end_opportunities(
                form_sp, boundmap, cat_forms, sub_names,
                prod_names, _add)

        # Deduplicate: multiple catalytic forms may
        # produce the same dead-end form. Group by
        # dead-end form name.
        de_forms = Dict{Symbol,
            Vector{Tuple{Symbol, Symbol}}}()
        for (f, m) in de_opportunities
            de_name = name(_add(form_sp[f], m))
            push!(get!(de_forms, de_name,
                Tuple{Symbol, Symbol}[]), (f, m))
        end
        de_form_names = sort(collect(keys(de_forms)))

        # Map each dead-end form to its bound metabolites
        de_bound = Dict{Symbol, Set{Symbol}}()
        for de_name in de_form_names
            entries = de_forms[de_name]
            f, m = first(entries)
            de_bound[de_name] = union(
                boundmap[f], Set([m]))
        end

        seen = Set{Vector{Symbol}}()

        for pattern in patterns
            # Filter dead-end forms by competition
            allowed_de = Symbol[]
            for de_name in de_form_names
                mets = de_bound[de_name]
                de_subs = intersect(mets, sub_names)
                de_prods = intersect(
                    mets, prod_names)
                has_conflict = any(
                    (s, p) in pattern
                    for s in de_subs
                    for p in de_prods)
                has_conflict || push!(
                    allowed_de, de_name)
            end

            # Dedup by form set
            allowed_de in seen && continue
            push!(seen, allowed_de)

            active_de = Set{Symbol}(allowed_de)

            # Build new steps: original topology + dead-end.
            # Each topology step is its own initial group
            # (group id = source position).
            steps = copy(topo)
            groups = collect(1:length(topo))
            next_g = length(topo) + 1

            # Add binding steps for active dead-ends.
            # Each binding step is an equivalence-eligible
            # candidate, but during initialization every binding
            # step gets its own fresh group (init_mechanisms
            # later applies same-metabolite grouping).
            for de_name in sort(collect(active_de))
                for (cat_form, met) in de_forms[de_name]
                    base = form_sp[cat_form]
                    push!(steps, Step(
                        base, _add(base, met), Metabolite[_role(met)], Metabolite[],
                        true))
                    push!(groups, next_g)
                    next_g += 1
                end
            end

            # Add mirror steps: for each catalytic
            # step, if both endpoints have dead-end
            # forms with the same metabolite, add a
            # parallel step. Mirror inherits RE/SS
            # AND the catalytic step's kinetic_group
            # (the step's source position in the topology).
            for (ci, s) in enumerate(topo)
                from = name(from_species(s))
                to = name(to_species(s))
                for de_met in sort(collect(all_mets))
                    de_met in boundmap[from] && continue
                    de_met in boundmap[to] && continue
                    from_de = name(_add(form_sp[from], de_met))
                    to_de = name(_add(form_sp[to], de_met))
                    from_de in active_de || continue
                    to_de in active_de || continue
                    push!(steps, Step(
                        _add(form_sp[from], de_met),
                        _add(form_sp[to], de_met),
                        consumed(s), released(s), is_equilibrium(s)))
                    push!(groups, ci)
                end
            end

            # Fully connect the enzyme-form graph: any two present forms that
            # are identical except for one bound metabolite must be joined by a
            # binding step. The dead-end + mirror steps above only cover
            # single-bystander cases; this fills multi-bystander gaps (e.g. a
            # dead-end form adjacent to another dead-end form). Added edges are
            # RE bindings on the differing metabolite — the equivalence grouping
            # folds them into that metabolite's kinetic group, and under rapid
            # equilibrium the extra edge is a thermodynamically-dependent cycle
            # that adds no free parameter. A no-op when no gaps exist (e.g.
            # bi-bi), so already-connected mechanisms are unaffected.
            present = Dict{Symbol, Species}()
            have_edge = Set{Tuple{Symbol, Symbol}}()
            for s in steps
                fr, to = from_species(s), to_species(s)
                present[name(fr)] = fr
                present[name(to)] = to
                push!(have_edge, (name(fr), name(to)))
                push!(have_edge, (name(to), name(fr)))
            end
            # Sort for deterministic edge/group order (Dict value order is not
            # guaranteed); matches the defensive sorting used when assembling
            # topologies above.
            forms_list = sort(collect(values(present)); by = name)
            for sp1 in forms_list, sp2 in forms_list
                conformation(sp1) == conformation(sp2) || continue
                residual(sp1) == residual(sp2) || continue
                b1 = Set(name(mb) for mb in bound(sp1))
                b2 = Set(name(mb) for mb in bound(sp2))
                (length(b2) == length(b1) + 1 && issubset(b1, b2)) || continue
                (name(sp1), name(sp2)) in have_edge && continue
                met = only(setdiff(b2, b1))
                push!(steps, Step(sp1, sp2, Metabolite[_role(met)], Metabolite[], true))
                push!(groups, next_g); next_g += 1
                push!(have_edge, (name(sp1), name(sp2)))
                push!(have_edge, (name(sp2), name(sp1)))
            end

            push!(result, (steps, groups))
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
    _to_group_list(steps, groups) -> Vector{Vector{Step}}

Partition a flat `Vector{Step}` into kinetic groups by the parallel
`groups` id vector, ordered by first occurrence of each group id in the
flat step list. (The `Mechanism` constructor then canonicalizes group and
step order, so this ordering is not load-bearing downstream.)
"""
function _to_group_list(steps::Vector{Step}, groups::Vector{Int})
    order = Int[]
    bygroup = Dict{Int, Vector{Step}}()
    for (s, g) in zip(steps, groups)
        haskey(bygroup, g) || push!(order, g)
        push!(get!(bygroup, g, Step[]), s)
    end
    [bygroup[g] for g in order]
end

"""
Reassign kinetic-group ids so steps sharing `(consumed, released, RE/SS)`
collapse into one group. Each multi-step class is assigned a fresh id;
singleton classes and iso steps keep their existing id. Operates on the
`(steps, groups)` parallel-array form and returns the merged pair.
"""
function _apply_equivalence_grouping(
    steps::Vector{Step}, groups::Vector{Int},
)
    classes = Dict{Tuple{Tuple, Tuple, Bool}, Vector{Int}}()
    for (i, s) in enumerate(steps)
        is_iso(s) && continue
        key = (Tuple(name.(consumed(s))), Tuple(name.(released(s))), is_equilibrium(s))
        push!(get!(classes, key, Int[]), i)
    end
    next_g = maximum(groups; init=0) + 1
    new_groups = copy(groups)
    for (_, idxs) in classes
        length(idxs) >= 2 || continue
        for i in idxs
            new_groups[i] = next_g
        end
        next_g += 1
    end
    (steps, new_groups)
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
raises the RE segment count (`_minimal_gaining_sets`): a single group when it
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
    flux = _flux_carrying_groups(_all_steady_state(steps(m)), reaction(m))
    flanks = _chain_flank_groups(m)
    units = [g for g in kinetic_groups(m)
             if all(is_equilibrium, steps(m)[g]) && flux[g] && !(g in flanks) &&
                !any(s -> any(x -> x isa Regulator, consumed(s)) ||
                          any(x -> x isa Regulator, released(s)), steps(m)[g])]
    flipped_groups(sel) = begin
        groups = steps(m)
        for u in sel
            groups = _flip_group_to_ss(groups, units[u])
        end
        groups
    end
    base = _re_segment_count(m)
    gains(sel) = begin
        _re_segment_count_after_flip(m, Set(units[u] for u in sel)) > base || return false
        groups = flipped_groups(sel)
        _bottomless_re_segment(groups) === nothing || return false
        flux = _flux_carrying_groups(groups, reaction(m))
        all(u -> flux[units[u]], sel)
    end
    sets = _minimal_gaining_sets(gains, _ -> 1:length(units))
    children = typeof(m)[_with_steps(m, flipped_groups(sel)) for sel in sets]
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
    group_of = Dict(s => g for (g, group) in enumerate(groups) for s in group)
    at = Dict{Species, Vector{Step}}()
    for group in groups, s in group, sp in (from_species(s), to_species(s))
        push!(get!(at, sp, Step[]), s)
    end
    flank(s0, x) = begin
        rest = filter(!=(s0), at[x])
        length(rest) == 1 || return nothing
        s = only(rest)
        is_binding(s) && to_species(s) == x && length(groups[group_of[s]]) == 1 ?
            s : nothing
    end
    out = Set{Int}()
    for group in groups
        length(group) == 1 || continue
        s0 = only(group)
        is_iso(s0) && !is_equilibrium(s0) || continue
        f1, f2 = flank(s0, from_species(s0)), flank(s0, to_species(s0))
        f1 === nothing || f2 === nothing || union!(out, (group_of[f1], group_of[f2]))
    end
    out
end
_chain_flank_groups(::AllostericMechanism) = Set{Int}()

"""
Return a fresh `Vector{Vector{Step}}` matching `groups` but with every
Step in group `g` rebuilt with `is_equilibrium=false`. All other groups
are reused by reference (Step is immutable).
"""
_flip_group_to_ss(groups::Vector{Vector{Step}}, g::Int) =
    [gi == g ? _with_equilibrium.(gr, false) : gr for (gi, gr) in enumerate(groups)]

"""`s` with its rapid-equilibrium flag set to `flag`; `s` itself when it has that flag, so
mechanisms built from one another share their unchanged steps."""
_with_equilibrium(s::Step, flag::Bool) = is_equilibrium(s) == flag ? s :
    Step(from_species(s), to_species(s), consumed(s), released(s), flag)

"""
Biconnected blocks of an undirected multigraph. `edges[e] = (u, v)` with
vertices `1:nv`. Returns the block id of every edge (Tarjan's edge-stack
algorithm); a bridge is a block of its own.
"""
function _edge_blocks(nv::Int, edges::Vector{Tuple{Int, Int}})
    adj = [Int[] for _ in 1:nv]
    for (e, (u, v)) in enumerate(edges)
        push!(adj[u], e); push!(adj[v], e)
    end
    disc = zeros(Int, nv); low = zeros(Int, nv)
    block = zeros(Int, length(edges))
    stack = Int[]; clock = Ref(0); nblocks = Ref(0)
    function visit(u, parent_edge)
        clock[] += 1; disc[u] = low[u] = clock[]
        for e in adj[u]
            e == parent_edge && continue
            w = edges[e][1] == u ? edges[e][2] : edges[e][1]
            if disc[w] == 0
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
    block
end

"""Every step of `groups` rebuilt at steady state. The flip move pre-filters its
units on this graph: a step that carries no flux with every step steady-state
carries none under any assignment (flips only refine segments), so a group with
no such step can never flip usefully."""
_all_steady_state(groups::Vector{Vector{Step}}) =
    [_with_equilibrium.(group, false) for group in groups]

"""`_re_segment_extras(groups)` with each form's index and segment: `(species, segments,
extras, index, segment_of)`, where `index[sp]` is the form's position in `species` and
`segment_of[i]` the segment holding form `i`."""
function _indexed_re_segments(groups::Vector{Vector{Step}})
    species, segments, extras = _re_segment_extras(groups)
    index = Dict(sp => i for (i, sp) in enumerate(species))
    segment_of = zeros(Int, length(species))
    for (k, members) in enumerate(segments), i in members
        segment_of[i] = k
    end
    species, segments, extras, index, segment_of
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
(`_edge_blocks`) holds a cycle of nonzero weight: every edge of such a block lies
on one (join the edge to the cycle by two disjoint paths; one of the two resulting
cycles is unbalanced), while in a balanced block detailed balance holds along every
cycle and each step's flux vanishes. The zero-flux verdict holds for any parameters
and grouping; the flux-carrying verdict for one-way steps, which every step the
enumerator emits is. The test reads each step's metabolite lists, so it holds for
fused and Theorell–Chance steps as well.
"""
function _flux_carrying_steps(groups::Vector{Vector{Step}}, rxn::EnzymeReaction)
    species, segments, extras, idx, seg = _indexed_re_segments(groups)
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

"""The block of every edge of the weighted multigraph (`_edge_blocks`) and the set of
blocks that hold a cycle of nonzero weight (`_block_balanced`)."""
function _unbalanced_blocks(nv::Int, edges::Vector{Tuple{Int, Int}}, weights::Vector{Int})
    block = _edge_blocks(nv, edges)
    members = Dict{Int, Vector{Int}}()
    for e in eachindex(edges)
        push!(get!(members, block[e], Int[]), e)
    end
    block, Set{Int}(b for (b, es) in members if !_block_balanced(edges, weights, es))
end

"""Whether every cycle of the block made of the edges `es` has weight zero. A BFS
spanning tree gives each vertex a potential, rising by an edge's weight along its
stored direction; the block is balanced iff every edge's weight equals the
potential difference of its ends. A block is connected, so one search from any of
its vertices visits all of them."""
function _block_balanced(edges, weights, es::Vector{Int})
    adj = Dict{Int, Vector{Int}}()
    for e in es
        u, v = edges[e]
        push!(get!(adj, u, Int[]), e); push!(get!(adj, v, Int[]), e)
    end
    root = edges[first(es)][1]
    phi = Dict(root => 0)
    queue = [root]
    while !isempty(queue)
        u = popfirst!(queue)
        for e in adj[u]
            a, b = edges[e]
            v = a == u ? b : a
            w = a == u ? weights[e] : -weights[e]
            if haskey(phi, v)
                phi[v] == phi[u] + w || return false
            else
                phi[v] = phi[u] + w
                push!(queue, v)
            end
        end
    end
    true
end

"""
    _flux_carrying_groups(groups, rxn) -> BitVector
    _flux_carrying_groups(m) -> BitVector

One flag per kinetic group: whether some step of the group carries net
steady-state flux (`_flux_carrying_steps`). A steady-state group with no such step
exposes only the ratio of its two constants, so the moves never emit one; a
zero-flux step inside a group that also holds a flux-carrying step costs nothing,
because the group's shared constants are pinned by the step that carries flux.
The mechanism method reads `steps(m)`, which for an allosteric mechanism is its
A-state catalytic graph.
"""
_flux_carrying_groups(groups::Vector{Vector{Step}}, rxn::EnzymeReaction) =
    BitVector([any(f) for f in _flux_carrying_steps(groups, rxn)])
_flux_carrying_groups(m::Union{Mechanism, AllostericMechanism}) =
    _flux_carrying_groups(steps(m), reaction(m))

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

"""Whether a cycle of rapid-equilibrium steps of `groups` runs net turnover, which makes the
rate infinite. Each RE step is an edge between its two forms weighted by its uptake of the
reaction's substrates minus its products (`_uptake_weight`); a turnover cycle has weight
n·Σρ², never zero, and a block of the RE graph holds one iff it is unbalanced
(`_unbalanced_blocks`)."""
function _re_turnover_cycle(groups::Vector{Vector{Step}}, rxn::EnzymeReaction)
    rho = _reactant_signs(rxn)
    idx = Dict{Species, Int}()
    vertex(sp) = get!(idx, sp, length(idx) + 1)
    edges = Tuple{Int, Int}[]; weights = Int[]
    for group in groups, s in group
        is_equilibrium(s) || continue
        push!(edges, (vertex(from_species(s)), vertex(to_species(s))))
        push!(weights, _uptake_weight(s, rho))
    end
    !isempty(last(_unbalanced_blocks(length(idx), edges, weights)))
end

"""Whether the rate has a maximum as the metabolites of `side` (`Substrate` or `Product`)
grow: turning rapid equilibrium every steady-state step that takes up or gives off one of
them leaves no rapid-equilibrium turnover cycle (condition V)."""
function _has_vmax(groups::Vector{Vector{Step}}, rxn::EnzymeReaction, side::Type)
    names = Set(name(x) for x in (side === Substrate ? substrates(rxn) : products(rxn)))
    touches(s) = _any_named(consumed(s), names) || _any_named(released(s), names)
    !_re_turnover_cycle([[touches(s) ? _with_equilibrium(s, true) : s for s in group]
                         for group in groups], rxn)
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
    re_iso = filter(is_iso, re)
    function node(start)
        out = Set{Species}(start); frontier = collect(start)
        while !isempty(frontier)
            f = pop!(frontier)
            for s in re_iso, (a, b) in ((from_species(s), to_species(s)),
                                        (to_species(s), from_species(s)))
                a == f && !(b in out) && (push!(out, b); push!(frontier, b))
            end
        end
        out
    end
    releases(n, names) = any(re) do t
        to_species(t) in n && !(from_species(t) in n) && _any_named(consumed(t), names) ||
            from_species(t) in n && !(to_species(t) in n) && _any_named(released(t), names)
    end
    merged = [to_species(s) for group in groups for s in group
              if _fused_substrate_binding(s, subs)]
    nodes = vcat([node([x]) for x in merged],
                 [node([from_species(s), to_species(s)]) for s in re_iso])
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
function of the mask `ss` of groups flipped to steady state, every other step at rapid
equilibrium, that is true when the candidate has no rapid-equilibrium turnover cycle
(`_re_turnover_cycle`), a maximal rate both ways (`_has_vmax`) and no chemistry in
equilibrium with both sides (`_chemistry_equilibrates_both_sides`). The forms and each
step's ends, uptake weight, group and reactants are indexed once per base, so a candidate
costs a few passes over arrays. A base holds no isomerization, so each chemistry node is a
merged complex or a Theorell–Chance step. A set of steps holds a turnover cycle iff giving
every form a potential that rises by each step's weight along it fails somewhere: a cycle
of nonzero weight is exactly a conflict, found here by a weighted union-find. Each
maximal-rate test runs on the candidate's rapid-equilibrium steps and more, and adding
steps never removes a conflict, so a turnover cycle of the candidate fails both tests.
"""
function _seed_candidate_screen(groups::Vector{Vector{Step}}, rxn::EnzymeReaction)
    flat = [s for group in groups for s in group]
    any(is_iso, flat) && error("_seed_candidate_screen: a seed base holds no isomerization")
    group = [g for (g, steps_g) in enumerate(groups) for _ in steps_g]
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
    function root(x)
        p = 0
        while parent[x] != x
            p += offset[x]; x = parent[x]
        end
        x, p
    end
    function turnover(on::BitVector)
        parent .= eachindex(parent); fill!(offset, 0)
        for k in eachindex(flat)
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
    meets(a::BitVector, b::BitVector) = any(k -> a[k] && b[k], eachindex(a))
    re = falses(length(flat)); kept = falses(length(flat))
    function screen(ss::AbstractVector{Bool})
        for k in eachindex(flat)
            re[k] = !ss[group[k]]
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
each base, every inclusion-minimal set of its groups whose flip to steady state gives a
valid, flux-carrying, non-degenerate candidate is a variant: no rapid-equilibrium turnover
cycle (`_re_turnover_cycle`), a maximal rate both ways (`_has_vmax`), no chemistry in
equilibrium with both sides (`_chemistry_equilibrates_both_sides`), no bottomless segment,
and every steady-state group carrying flux. The first three tests read index arrays
(`_seed_candidate_screen`); a candidate's steps are built only when it passes them. A
merged variant whose every merged complex has both its steps steady state, each alone in
its group, is skipped: the unmerged form with rapid-equilibrium flanks has its family at
the same count.
"""
function _seed_variants(m::Mechanism)
    rxn = reaction(m)
    isos = [s for group in steps(m) for s in group if is_iso(s)]
    merged = foldl(_merge_isomerization, isos; init = steps(m))
    merged = [_with_equilibrium.(group, true) for group in merged]
    complexes = [to_species(s) for s in isos]
    bases = Tuple{Vector{Vector{Step}}, Bool}[(merged, true)]
    for x in complexes
        b = _eliminate_form(merged, x)
        b === nothing || push!(bases, (b, false))
    end
    lumping_twin(gs) = all(complexes) do x
        at_x = [(s, group) for group in gs for s in group
                if x in (from_species(s), to_species(s))]
        length(at_x) == 2 &&
            all(((s, group),) -> !is_equilibrium(s) && length(group) == 1, at_x)
    end
    variants = Mechanism[]
    for (base, is_merged) in bases
        for group in base, s in group
            _assert_step_atom_conserving(rxn, s)
        end
        steady = _all_steady_state(base)
        flux = _flux_carrying_groups(steady, rxn)
        units = [g for g in eachindex(base) if flux[g]]
        screen = _seed_candidate_screen(base, rxn)
        mask = falses(length(base))
        # A candidate takes each flipped group from `steady`, so the variants of a base
        # share their steady-state steps.
        flipped(sel) = begin
            gs = copy(base)
            for u in sel
                gs[units[u]] = steady[units[u]]
            end
            gs
        end
        admissible(sel) = begin
            fill!(mask, false)
            for u in sel
                mask[units[u]] = true
            end
            screen(mask) || return false
            gs = flipped(sel)
            _bottomless_re_segment(gs) === nothing || return false
            carries = _flux_carrying_groups(gs, rxn)
            all(g -> is_equilibrium(first(gs[g])) || carries[g], eachindex(gs))
        end
        for sel in _minimal_gaining_sets(admissible, _ -> 1:length(units))
            gs = flipped(sel)
            is_merged && lumping_twin(gs) && continue
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
    groups = filter(!isempty, [filter(on_catalytic_site, group) for group in steps(m)])
    species, segments, extras, idx, seg_of = _indexed_re_segments(groups)
    # Directed segment-graph edges: source segment, target segment, source form,
    # metabolites bound in that direction.
    edges = Tuple{Int, Int, Int, Vector{Symbol}}[]
    for group in groups, s in group
        is_equilibrium(s) && continue
        m_lhs = Symbol[name(x) for x in consumed(s)]
        m_rhs = Symbol[name(x) for x in released(s)]
        a, b = idx[from_species(s)], idx[to_species(s)]
        seg_of[a] == seg_of[b] && continue
        push!(edges, (seg_of[a], seg_of[b], a, m_lhs))
        push!(edges, (seg_of[b], seg_of[a], b, m_rhs))
    end
    rxn = reaction(m)
    mets = vcat(Symbol[name(s) for s in substrates(rxn)],
                Symbol[name(p) for p in products(rxn)])
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

"""Number of rapid-equilibrium segments (connected components of the RE
subgraph). An allosteric mechanism is measured on its A-state projection, which
holds every catalytic group."""
_re_segment_count(m::Mechanism) = length(_compute_re_groups(m)[2])
_re_segment_count(am::AllostericMechanism) = _re_segment_count(_state_mechanism(am, :A))

"""RE segment count of `m` after flipping the kinetic groups in `flipped` to SS,
computed on `m`'s own steps without building the child: a union-find over the
species joined by the RE steps of every other group. `_compute_re_groups`
gives the same answer on the built child; this form is what the flip move probes
with, so a candidate the constructors would reject is never constructed."""
function _re_segment_count_after_flip(
    m::Union{Mechanism, AllostericMechanism}, flipped,
)
    species = Species[]
    for group in steps(m), s in group
        from_species(s) in species || push!(species, from_species(s))
        to_species(s) in species   || push!(species, to_species(s))
    end
    parent = collect(1:length(species))
    function find(x)
        while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end
        x
    end
    for (g, group) in enumerate(steps(m))
        g in flipped && continue
        for s in group
            is_equilibrium(s) || continue
            ra = find(findfirst(==(from_species(s)), species))
            rb = find(findfirst(==(to_species(s)),   species))
            ra != rb && (parent[ra] = rb)
        end
    end
    count(i -> find(i) == i, eachindex(species))
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
    for g in kinetic_groups(m), bp in _context_bipartitions(groups[g])
        parts = _revert_zero_flux_parts(groups[g], bp, flux[g])
        push!(units, (g, parts)); push!(reverted, parts !== bp)
    end
    isempty(units) && return typeof(m)[]
    selection(sel) = [units[u] for u in sel]
    gain = _split_gain_test(m, units)
    gains(sel) =
        (!any(u -> reverted[u], sel) ||
         _bottomless_re_segment(_bipartitioned_groups(groups, selection(sel))[1]) ===
         nothing) && gain(sel)
    segments = _group_re_segments(m)
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
            typeof(m)[_apply_bipartitions(m, selection(sel)) for sel in sets])
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
        next_id = length(steps(m))
        for u in sel
            next_id += 1
            for (part_number, part) in enumerate(units[u][2]), s in part
                j = position[reaction_key(s)]
                part_number == 2 && (ids[j] = next_id)
                kinds[j] = _count_kind(s)
            end
        end
        counter(ids, kinds) > base
    end
    gains
end

function _split_gain_test(am::AllostericMechanism, units)
    base = _independent_param_count(am)
    sel -> _independent_param_count(
        _apply_bipartitions(am, [units[u] for u in sel])) > base
end

"""RE segment ids touched by each kinetic group's steps; empty for a group holding
an SS step, which lies on no RE cycle and is never a split partner."""
function _group_re_segments(m::Union{Mechanism, AllostericMechanism})
    cm = m isa Mechanism ? m : _state_mechanism(m, :A)
    species, _, form_to_segment = _compute_re_groups(cm)
    segment(sp) = form_to_segment[findfirst(==(sp), species)]
    map(steps(m)) do group
        all(is_equilibrium, group) || return Set{Int}()
        Set{Int}(segment(from_species(s)) for s in group)
    end
end

"""
The endpoint of `s` that does not carry the step's free metabolites: the form
they bind to. `from_species(s)` when the step consumes a metabolite (every
binding, plain or fused, is stored with its metabolite consumed), `to_species(s)`
when it only releases metabolites (two or more: a step releasing one is stored
as the binding it reverses), and `from_species(s)` for an iso step (both lists
empty).
"""
function _context_form(s::Step)
    isempty(consumed(s)) || return from_species(s)
    isempty(released(s)) || return to_species(s)
    from_species(s)
end

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
    ligands = Set{Metabolite}()
    for f in forms, b in bound(f)
        own !== nothing && b == own && continue
        push!(ligands, b)
    end
    division(carries) = Step[s for (s, f) in zip(group, forms) if carries(f)]
    divisions = [division(f -> y in bound(f))
                 for y in sort!(collect(ligands);
                                by = b -> (string(typeof(b)), string(name(b))))]
    append!(divisions, division(f -> conformation(f) == c)
            for c in sort!(unique(conformation(f) for f in forms); by = string))
    append!(divisions, division(f -> residual(f) == r)
            for r in sort!(unique(residual(f) for f in forms); by = string))
    seen = Set{Vector{Step}}()
    out = Tuple{Vector{Step}, Vector{Step}}[]
    for with in divisions
        without = Step[s for s in group if !(s in with)]
        (isempty(with) || isempty(without)) && continue
        first_part, second_part = first(group) in with ? (with, without) : (without, with)
        first_part in seen && continue
        push!(seen, first_part)
        push!(out, (first_part, second_part))
    end
    out
end

"""
Replace each selected group by its two bipartition parts. `selection` pairs a
group index with one of that group's `_context_bipartitions`; each group appears
at most once. For an allosteric mechanism both parts inherit the group's
catalytic allo-state tag (splitting is a parameter-relaxation move that must not
change A/I semantics).
"""
function _apply_bipartitions(m::Mechanism, selection)
    _with_steps(m, _bipartitioned_groups(steps(m), selection)[1])
end

function _apply_bipartitions(am::AllostericMechanism, selection)
    groups, origin = _bipartitioned_groups(steps(am), selection)
    _with_steps_and_cat_states(am, groups, cat_allo_states(am)[origin])
end

"""
The bipartition `bp` of `group` with every part none of whose steps carries flux
(`flags`, one per step of `group`) rebuilt at rapid equilibrium; `bp` itself when
no part changes, and always for a rapid-equilibrium group (see `_flux_carrying_groups`).
"""
function _revert_zero_flux_parts(group::Vector{Step}, bp, flags::BitVector)
    is_equilibrium(first(group)) && return bp
    carries = Dict(s => flags[j] for (j, s) in enumerate(group))
    revert(part) = any(s -> carries[s], part) ? part :
        _with_equilibrium.(part, true)
    p1, p2 = revert(bp[1]), revert(bp[2])
    p1 === bp[1] && p2 === bp[2] ? bp : (p1, p2)
end

"""Groups of `groups` with each selected group replaced by its two parts, plus
the index of the original group each new group came from."""
function _bipartitioned_groups(groups::Vector{Vector{Step}}, selection)
    parts = Dict(g => bp for (g, bp) in selection)
    out = Vector{Vector{Step}}()
    origin = Int[]
    for (g, group) in enumerate(groups)
        if haskey(parts, g)
            push!(out, parts[g][1]); push!(origin, g)
            push!(out, parts[g][2]); push!(origin, g)
        else
            push!(out, group); push!(origin, g)
        end
    end
    out, origin
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
    species, segments, extras, idx, seg = _indexed_re_segments(groups)
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
    shared = Dict{Species, Int}()
    for t in values(twin_of)
        shared[t] = get(shared, t, 0) + 1
    end
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
        h = shared[t] == 1 ? get(binding_group, (from_species(s), t), 0) : 0
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
(`_reachable_from_free`). The copy binds there only at the sites that graph keeps, and a
free-enzyme site is always among them. In a state where the copy binds nothing, every
factor of that state is 1. The two states form one system. A group tagged `:EqualAI`
has one set of constants in both, so it must take the same rescaling in each. A shared
class is one number in both states only when its binding group and the copy are both
`:EqualAI` (one K_h, one K*). The copy's complexes share one factor in both states when
the copy is `:EqualAI` (one K*). Every other class is a number per state.
"""
function _redundant_copy_groups(m::Mechanism)
    groups = steps(m)
    any(group -> bound_metabolite(first(group)) isa CompetitiveInhibitor, groups) ||
        return Int[]
    twin = _productive_twin(groups)
    [g for g in eachindex(groups)
     if _gauge_rescaling(groups, g, twin, identity) !== nothing]
end
function _redundant_copy_groups(am::AllostericMechanism)
    active = steps(am)
    any(group -> bound_metabolite(first(group)) isa CompetitiveInhibitor, active) ||
        return Int[]
    tags = cat_allo_states(am)
    reachable = _reachable_from_free(
        [group for (g, group) in enumerate(active) if tags[g] !== :OnlyA])
    inactive = [tags[g] === :OnlyA ? Step[] :
                Step[s for s in group if name(from_species(s)) in reachable &&
                                         name(to_species(s)) in reachable]
                for (g, group) in enumerate(active)]
    twin_active, twin_inactive = _productive_twin(active), _productive_twin(inactive)
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
Extend `rxn`'s regulators with a new `CompetitiveInhibitor(name)`,
preserving every other field (reactants, allowed catalytic
multiplicities, shared_catalytic_site). `EnzymeReaction`'s inner
constructor canonicalizes the regulator order.
"""
function _add_competitive_inhibitor(rxn::EnzymeReaction, reg_name::Symbol)
    any(rm -> name(regulator(rm)) == reg_name, regulators(rxn)) && return rxn
    new_regs = copy(regulators(rxn))
    push!(new_regs, RegulatorMults(CompetitiveInhibitor(reg_name), Int[1]))
    EnzymeReaction(copy(reactants(rxn)), new_regs,
                   copy(allowed_catalytic_multiplicities(rxn));
                   shared_catalytic_site = copy(shared_catalytic_site(rxn)))
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
                                   rxn::EnzymeReaction; exclude_regs) →
        Vector{typeof(m)}

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

The caller must pass the declared `rxn` because `m.reaction` only
carries regulators already bound by its steps; not-yet-bound regulators
live exclusively in the declared reaction. The new `Mechanism`'s
reaction is `rxn` extended with the newly-bound regulator (preserving
the substrate / product / multiplicity payload).
"""
function _expand_add_dead_end_regulator(
    m::Union{Mechanism, AllostericMechanism}, rxn::EnzymeReaction;
    exclude_regs::Set{Symbol}=Set{Symbol}(),
)
    isempty(regulators(rxn)) && return typeof(m)[]
    additional_excluded = m isa AllostericMechanism ?
        Set(name(l) for site in regulatory_sites(m) for l in ligands(site)) :
        Set{Symbol}()

    sub_names = Set(name(s) for s in substrates(rxn))
    prod_names = Set(name(p) for p in products(rxn))

    existing_regs = Set{Symbol}()
    for group in steps(m), s in group, bm in consumed(s)
        bm isa Regulator && push!(existing_regs, name(bm))
    end

    eligible_regs = Symbol[]
    for rm in regulators(rxn)
        reg = regulator(rm)
        reg isa CompetitiveInhibitor || continue
        name(reg) in existing_regs && continue
        name(reg) in additional_excluded && continue
        name(reg) in exclude_regs && continue
        push!(eligible_regs, name(reg))
    end
    sort!(eligible_regs)
    isempty(eligible_regs) && return typeof(m)[]

    form_sp = Dict{Symbol, Species}()
    for group in steps(m), s in group
        form_sp[name(from_species(s))] = from_species(s)
        form_sp[name(to_species(s))] = to_species(s)
    end
    cat_forms = Set(keys(form_sp))
    # Per form, the names bound productively and the names bound as competitive
    # inhibitors: a reactant's copy sits at a dead-end site, not at its reactant's.
    productive = Dict(f => Set(name(b) for b in bound(sp) if b isa Reactant)
                      for (f, sp) in form_sp)
    inhibiting = Dict(f => Set(name(b) for b in bound(sp) if b isa CompetitiveInhibitor)
                      for (f, sp) in form_sp)

    results = typeof(m)[]

    twin = _productive_twin(steps(m))

    for reg_name in eligible_regs
        eligible_forms = Symbol[]
        for f in sort(collect(cat_forms))
            fb = productive[f]
            (intersect(fb, sub_names) == sub_names ||
                intersect(fb, prod_names) == prod_names) && continue
            push!(eligible_forms, f)
        end
        isempty(eligible_forms) && continue

        existing_inhibitors = Symbol[]
        for group in steps(m), s in group, bm in consumed(s)
            bm isa Regulator || continue
            name(bm) == reg_name && continue
            push!(existing_inhibitors, name(bm))
        end
        sort!(unique!(existing_inhibitors))

        inh_patterns = _inhibitor_competition_patterns(
            sub_names, prod_names, existing_inhibitors)
        seen = Set{Vector{Symbol}}()

        for (comp_subs, comp_prods, comp_inhibitors) in inh_patterns
            comp_reactants = union(comp_subs, comp_prods)
            target_forms = Set{Symbol}()
            for met in comp_reactants
                union!(target_forms, _forms_where_free(m, Reactant, met))
            end
            for inh in comp_inhibitors
                union!(target_forms, _forms_where_free(m, CompetitiveInhibitor, inh))
            end

            active = Symbol[]
            for f in sort(collect(target_forms))
                f in eligible_forms || continue
                isempty(intersect(productive[f], comp_reactants)) || continue
                isempty(intersect(inhibiting[f], comp_inhibitors)) || continue
                push!(active, f)
            end
            isempty(active) && continue
            active in seen && continue
            push!(seen, active)

            de_species_map = Dict{Symbol, Species}()
            reg_group_steps = Step[]
            for cf in active
                base = form_sp[cf]
                # The complex the copy forms at `base`, in its conformation and residual.
                de_species = Species(
                    Metabolite[bound(base)..., CompetitiveInhibitor(reg_name)],
                    conformation(base), residual(base))
                de_species_map[cf] = de_species
                push!(reg_group_steps, Step(
                    base, de_species,
                    Metabolite[CompetitiveInhibitor(reg_name)], Metabolite[], true))
            end

            mirror_per_group = Dict{Int, Vector{Step}}()
            for (gi, group) in enumerate(steps(m))
                for s in group
                    fn = name(from_species(s))
                    tn = name(to_species(s))
                    haskey(de_species_map, fn) || continue
                    haskey(de_species_map, tn) || continue
                    push!(get!(mirror_per_group, gi, Step[]),
                        Step(de_species_map[fn], de_species_map[tn],
                             consumed(s), released(s), is_equilibrium(s)))
                end
            end

            new_groups = Vector{Vector{Step}}()
            for (gi, group) in enumerate(steps(m))
                extended = copy(group)
                haskey(mirror_per_group, gi) &&
                    append!(extended, mirror_per_group[gi])
                push!(new_groups, extended)
            end
            push!(new_groups, reg_group_steps)

            child = _dead_end_child(m, new_groups,
                                    _add_competitive_inhibitor(rxn, reg_name), twin)
            child === nothing || push!(results, child)
        end
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

For each value in `rxn`'s `allowed_catalytic_multiplicities`, the
multiplicity becomes the variant's `catalytic_multiplicity`. Catalytic
steps are reused by reference; duplicate variants are removed.

A parent whose catalytic scheme fails `_hyperbolic_catalysis` (random-order
steady-state binding, or a substrate that traps a steady-state intermediate in
an abortive complex, whose own equation carries concentration powers) emits
variants at multiplicity 1 only: above one subunit the conformational
equilibrium would add a second source of powers
(`_requires_hyperbolic_catalysis`).
"""
function _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
    hyperbolic = _hyperbolic_catalysis(m)
    n_g = length(steps(m))
    chem = [g for g in 1:n_g if any(_is_chemistry, steps(m)[g])]
    bind = [g for g in 1:n_g if !(g in chem)]
    regs = Symbol[]
    for rm in regulators(rxn)
        reg = regulator(rm)
        reg isa AllostericRegulator && push!(regs, name(reg))
    end
    sort!(regs)
    results = AllostericMechanism[]
    for cn in allowed_catalytic_multiplicities(rxn)
        cn > 1 && !hyperbolic && continue
        # K-type: every non-empty subset of binding groups :OnlyA, with every
        # chemistry group :OnlyA — a catalytically-dead inactive conformation.
        # A state that cannot bind a catalytic metabolite cannot complete the
        # cycle, so it runs no chemistry. Each is emitted bare; over one subunit
        # L can be a phantom (see the docstring). `_onlya_haldane_violation` drops
        # a subset that leaves a binding-only Wegscheider cycle unsatisfiable.
        for mask in 1:(2^length(bind) - 1)
            tags = Symbol[:EqualAI for _ in 1:n_g]
            for g in chem
                tags[g] = :OnlyA
            end
            for (i, g) in enumerate(bind)
                (mask >> (i - 1)) & 1 == 1 && (tags[g] = :OnlyA)
            end
            _onlya_haldane_violation(rxn, steps(m), tags) === nothing || continue
            push!(results, AllostericMechanism(
                reaction(m), copy(steps(m)), tags, cn, RegulatorySite[]))
        end
        # V-type: no :OnlyA binding, every chemistry group :OnlyA. The inactive
        # state binds identically but cannot catalyze, so L folds into kcat and is
        # unobservable — emit only paired with a declared regulator.
        if !isempty(chem) && !isempty(regs)
            vtags = Symbol[:EqualAI for _ in 1:n_g]
            for g in chem
                vtags[g] = :OnlyA
            end
            am_cat = AllostericMechanism(
                reaction(m), copy(steps(m)), vtags, cn, RegulatorySite[])
            for reg in regs, tag in (:OnlyA, :OnlyI)
                push!(results, _make_am_with_added_reg(am_cat, reg, tag, 0))
            end
        end
    end
    unique!(results)
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
    one non-`:EqualAI` ligand (otherwise the `RegulatorySite`
    constructor's all-`:EqualAI` rule would reject the variant).

New sites inherit `am.catalytic_multiplicity` as their multiplicity.
The mechanism's catalytic side, regulatory_sites' multiplicities, and
reaction payload pass through unchanged.

Caller must supply `rxn` because `am.reaction` only carries regulators
already bound by its steps; not-yet-bound regulators live in the
declared reaction.
"""
function _expand_add_allosteric_regulator(
    am::AllostericMechanism, rxn::EnzymeReaction,
)
    existing_allo = Set{Symbol}()
    for site in regulatory_sites(am), lig in ligands(site)
        push!(existing_allo, name(lig))
    end

    new_regs = Symbol[]
    for rm in regulators(rxn)
        reg = regulator(rm)
        reg isa AllostericRegulator || continue
        name(reg) in existing_allo && continue
        push!(new_regs, name(reg))
    end
    sort!(new_regs)
    isempty(new_regs) && return AllostericMechanism[]

    results = AllostericMechanism[]
    for reg in new_regs
        n_sites = length(regulatory_sites(am))
        # Non-:EqualAI tags at any (new or existing) site.
        for tag in (:OnlyA, :OnlyI, :NonequalAI)
            for site_idx in 0:n_sites
                # Appending to an existing site whose conformations are disjoint
                # from the new ligand's (an all-:OnlyA site gaining an :OnlyI
                # ligand, or the reverse) reproduces the add-at-a-new-site
                # equation — redundant, so skip it.
                site_idx >= 1 && isempty(intersect(_state_conformations(tag),
                    _site_active_states(regulatory_sites(am)[site_idx]))) && continue
                push!(results,
                    _make_am_with_added_reg(am, reg, tag, site_idx))
            end
        end
        # :EqualAI at an existing site only when that site already has
        # at least one non-:EqualAI ligand (avoids the constructor's
        # all-:EqualAI single-ligand rejection / identical-cancellation).
        for site_idx in 1:n_sites
            site = regulatory_sites(am)[site_idx]
            any(st != :EqualAI for st in allo_states(site)) || continue
            push!(results,
                _make_am_with_added_reg(am, reg, :EqualAI, site_idx))
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
    new_sites = RegulatorySite[]
    if site_idx == 0
        for site in regulatory_sites(am)
            push!(new_sites, site)
        end
        push!(new_sites, RegulatorySite(
            AllostericRegulator[AllostericRegulator(reg)],
            catalytic_multiplicity(am),
            Symbol[tag]))
    else
        for (i, site) in enumerate(regulatory_sites(am))
            if i == site_idx
                new_ligs = copy(ligands(site))
                push!(new_ligs, AllostericRegulator(reg))
                new_states = copy(allo_states(site))
                push!(new_states, tag)
                push!(new_sites, RegulatorySite(
                    new_ligs, multiplicity(site), new_states))
            else
                push!(new_sites, site)
            end
        end
    end
    AllostericMechanism(reaction(am), copy(steps(am)),
                        copy(cat_allo_states(am)),
                        catalytic_multiplicity(am), new_sites)
end

"""
    _partial_onlya_catalysis(cat_steps, cat_allo_states) → Bool

True when the inactive conformation catalyzes only partially: some catalytic
group is `:OnlyA` (a dead binding or chemistry group) while some chemistry
group (one holding a chemistry step, `_is_chemistry`) is still live (not
`:OnlyA`). The inactive conformation's catalysis must be all-or-nothing — fully
dead (every chemistry group `:OnlyA`, whether because an `:OnlyA` binding blocks
the cycle or by a dead-inactive V-type) or fully live. A partial conformation
strands enzyme in a covalent form (a kinetic sink) and crashes the
saturating-turnover extraction. The enumeration moves use this to avoid
generating such a form.
"""
function _partial_onlya_catalysis(cat_steps::Vector{Vector{Step}},
                                  cat_allo_states::Vector{Symbol})
    live = any(any(_is_chemistry, cat_steps[g]) && cat_allo_states[g] !== :OnlyA
               for g in eachindex(cat_steps))
    live && any(==(:OnlyA), cat_allo_states)
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

Relaxing an `:OnlyA` chemistry group is dropped in two cases. A one-sided
`:OnlyA` binding is only legal because `k_I = 0`; restoring a finite
`k_I` strands it, leaving no thermodynamic reading
(`_onlya_haldane_violation`). More broadly, inactive catalysis must be
all-or-nothing: a relaxation that leaves the inactive conformation
catalyzing only partially — some catalytic group `:OnlyA` while a
chemistry group stays live — is dropped (`_partial_onlya_catalysis`),
because such a conformation strands enzyme in a covalent form. Relaxing
an `:OnlyA` binding, or a chemistry group of a fully-live inactive
conformation, is retained. Where the inactive conformation binds nothing
(all bindings `:OnlyA`), the dropped partial variant is rate-equivalent
to the fully-dead form emitted directly, so no observable hypothesis is
lost.

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
    results = AllostericMechanism[]
    cs = steps(am)
    chem = [g for g in eachindex(cs) if any(_is_chemistry, cs[g])]

    # Binding catalytic groups relax individually.
    for g in eachindex(cat_allo_states(am))
        g in chem && continue
        cat_allo_states(am)[g] == :NonequalAI && continue
        new_states = copy(cat_allo_states(am))
        new_states[g] = :NonequalAI
        _onlya_haldane_violation(reaction(am), cs, new_states) ===
            nothing || continue
        _partial_onlya_catalysis(cs, new_states) && continue
        push!(results, _with_cat_allo_states(am, new_states))
    end

    # Chemistry groups relax together. Inactive catalysis is all-or-nothing, so a
    # fully-`:NonequalAI` catalytic inactive conformation is unreachable by
    # relaxing chemistry groups one at a time — each mixed intermediate is a
    # partial and is dropped. One variant sets every non-`:NonequalAI` chemistry
    # group to `:NonequalAI` at once.
    if any(cat_allo_states(am)[g] != :NonequalAI for g in chem)
        new_states = copy(cat_allo_states(am))
        for g in chem
            new_states[g] = :NonequalAI
        end
        if _onlya_haldane_violation(reaction(am), cs, new_states) === nothing &&
           !_partial_onlya_catalysis(cs, new_states)
            push!(results, _with_cat_allo_states(am, new_states))
        end
    end

    for (si, site) in enumerate(regulatory_sites(am))
        for (li, _) in enumerate(ligands(site))
            allo_states(site)[li] == :NonequalAI && continue
            new_sites = copy(regulatory_sites(am))
            new_states = copy(allo_states(site))
            new_states[li] = :NonequalAI
            new_sites[si] = RegulatorySite(
                copy(ligands(site)), multiplicity(site), new_states)
            push!(results, _with_reg_sites(am, new_sites))
        end
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
`expand_mechanisms` runs every child through `_filter_by_reg_type`.
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
        for states in _merged_site_state_assignments(base_states;
                                                     drop_all_keep=redundant)
            merged = RegulatorySite(copy(ligs), mult, states)
            push!(results, _with_reg_sites(am, vcat(others, [merged])))
        end
    end
    results
end

"""
    _merged_site_state_assignments(base_states::Vector{Symbol};
                                   drop_all_keep=false) -> Vector{Vector{Symbol}}

Δ0-valid allo-state assignments for a merged site's ligands: the all-keep
assignment (omitted when `drop_all_keep`), plus each assignment retagging
exactly one non-`:EqualAI` ligand to `:EqualAI`. The all-`:EqualAI` result is
dropped. `drop_all_keep` omits the all-keep entry for a redundant merge (see
`_site_active_states`) while keeping the antagonist retags.
"""
function _merged_site_state_assignments(base_states::Vector{Symbol};
                                        drop_all_keep::Bool=false)
    assignments = Vector{Symbol}[]
    drop_all_keep || push!(assignments, copy(base_states))
    for i in eachindex(base_states)
        base_states[i] == :EqualAI && continue
        retagged = copy(base_states)
        retagged[i] = :EqualAI
        all(==(:EqualAI), retagged) && continue
        push!(assignments, retagged)
    end
    assignments
end

# ─── Regulator-Type Filter ────────────────────────────────────

"""
    _declared_reg_type(reg_name::Symbol, rxn::EnzymeReaction) -> Symbol

The declared type (`:activator`, `:inhibitor`, or `:unspecified`) of the
`AllostericRegulator` named `reg_name` in `rxn`'s `regulators`. Returns
`:unspecified` when `reg_name` names no `AllostericRegulator` entry (either
absent entirely, or present only as some other `Regulator` subtype).
"""
function _declared_reg_type(reg_name::Symbol, rxn::EnzymeReaction)
    for rm in regulators(rxn)
        reg = regulator(rm)
        reg isa AllostericRegulator && name(reg) == reg_name && return reg_type(rm)
    end
    :unspecified
end

"""
    _state_respects_reg_type(reg_name, state::Symbol, sibling_names::Vector{Symbol},
                          rxn::EnzymeReaction) -> Bool

Whether ligand `reg_name`'s allosteric `state` is consistent with its
declared regulator type, given `sibling_names` — the OTHER ligands at its
regulatory site. `:unspecified` type always passes. A designated activator
is never `:OnlyI` and a designated inhibitor is never `:OnlyA`; the
type-matching pure state and `:NonequalAI` are always allowed. `:EqualAI`
is an antagonist state and is rejected only when a sibling is a
same-type designated effector, since that would counteract the sibling's
declared direction.
"""
function _state_respects_reg_type(reg_name, state::Symbol,
                              sibling_names::Vector{Symbol}, rxn::EnzymeReaction)
    rt = _declared_reg_type(reg_name, rxn)
    rt === :unspecified && return true
    rt === :activator && state === :OnlyI && return false
    rt === :inhibitor && state === :OnlyA && return false
    if state === :EqualAI
        any(s -> _declared_reg_type(s, rxn) === rt, sibling_names) && return false
    end
    true
end

"""
    _filter_by_reg_type(mechs::Vector, rxn::EnzymeReaction) -> Vector

Keep only mechanisms whose every regulatory ligand respects its declared
type (`_state_respects_reg_type`) given its site's other ligands. A `Mechanism`
has no regulatory sites and passes trivially.
"""
_filter_by_reg_type(mechs::Vector, rxn::EnzymeReaction) =
    filter(m -> _respects_reg_type(m, rxn), mechs)

"""Whether every regulatory ligand in `m` respects its declared type (see
`_filter_by_reg_type`)."""
_respects_reg_type(::Mechanism, ::EnzymeReaction) = true
function _respects_reg_type(am::AllostericMechanism, rxn::EnzymeReaction)
    for site in regulatory_sites(am)
        lig_names = Symbol[name(lig) for lig in ligands(site)]
        for (i, lig_name) in enumerate(lig_names)
            siblings = Symbol[lig_names[j] for j in eachindex(lig_names) if j != i]
            _state_respects_reg_type(lig_name, allo_states(site)[i], siblings, rxn) ||
                return false
        end
    end
    true
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
    result = _filter_by_reg_type(result, rxn)
    for child in result
        _assert_atom_conserving(child)
    end
    result
end

# --- Dedup ---


"""
    init_mechanisms(reaction::EnzymeReaction) -> Vector{Mechanism}

Public entry point. Produces the starting mechanisms of the search for a
reaction as concrete `Mechanism` structs: first the seeds, then their merged
and Theorell–Chance variants. A seed is built for each catalytic topology
(`_catalytic_topologies`) and each substrate/product dead-end subset
(`_expand_substrate_product_dead_ends`), with one steady-state step, a
chemistry isomerization, and binding steps sharing the same `(metabolite, RE/SS)`
class collapsed into one kinetic group (`_apply_equivalence_grouping`). Dead-end
enumeration respects `shared_catalytic_site`. The variants of each seed
(`_seed_variants`) follow in the seeds' order, each once. The parameter counts
are mixed: a three-group merged variant or a decorated Theorell–Chance variant
fits more parameters than its seed, a merged ping-pong variant fewer.
"""
function init_mechanisms(r::EnzymeReaction)
    topos = _catalytic_topologies(r)
    expanded = _expand_substrate_product_dead_ends(topos, r)
    mechs = Mechanism[]
    for (steps, groups) in expanded
        merged_steps, merged_groups =
            _apply_equivalence_grouping(steps, groups)
        m = Mechanism(r, _to_group_list(merged_steps, merged_groups))
        _assert_atom_conserving(m)
        push!(mechs, m)
    end
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
4. every ligand respects its declared type (`_filter_by_reg_type`).

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
        _binds_all_required(m, required_allo, required_comp) && push!(seeds, m)
        true
    end
    frontier = Union{Mechanism, AllostericMechanism}[
        m for m in init_mechanisms(rxn) if consider!(m)]
    while !isempty(frontier)
        # Per-node child generation is pure and independent — distribute it.
        childsets = pmap(frontier) do m
            filter(c -> _is_seed_node(c, rxn, required_allo, required_comp),
                   _seed_children(m, rxn, required_allo))
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

"""
    _seed_children(m, rxn::EnzymeReaction, required_allo::Set{Symbol})
        -> Vector{Union{Mechanism, AllostericMechanism}}

Children of `m` under the seed-build structure moves. The two
allosteric-lifting moves run only when an allosteric regulator is required, so
a competitive-only required set (`required_allo` empty) stays non-allosteric —
no `L` — and seeds at `base + n_required_comp`. The dead-end move always runs.
A `Mechanism` is lifted by `_expand_to_allosteric`; an `AllostericMechanism`
gains a regulator by `_expand_add_allosteric_regulator`.
"""
function _seed_children(m::Union{Mechanism, AllostericMechanism},
                        rxn::EnzymeReaction, required_allo::Set{Symbol})
    children = Union{Mechanism, AllostericMechanism}[]
    if !isempty(required_allo)
        append!(children, m isa Mechanism ? _expand_to_allosteric(m, rxn) :
                          _expand_add_allosteric_regulator(m, rxn))
    end
    append!(children, _expand_add_dead_end_regulator(m, rxn))
    children
end

"""Whether `m` carries a `:NonequalAI` tag, on the catalytic step or a regulatory
site."""
_has_nonequalai(::Mechanism) = false
_has_nonequalai(am::AllostericMechanism) =
    any(==(:NonequalAI), cat_allo_states(am)) ||
    any(site -> any(==(:NonequalAI), allo_states(site)), regulatory_sites(am))

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
    _has_nonequalai(m) && return false
    m isa AllostericMechanism &&
        !all(site -> length(ligands(site)) == 1, regulatory_sites(m)) &&
        return false
    issubset(_bound_allo_regs(m), required_allo) || return false
    issubset(_bound_comp_inhibitors(m), required_comp) || return false
    _respects_reg_type(m, rxn) || return false
    true
end

_binds_all_required(m::Union{Mechanism, AllostericMechanism},
                    required_allo::Set{Symbol}, required_comp::Set{Symbol}) =
    issubset(required_allo, _bound_allo_regs(m)) &&
    issubset(required_comp, _bound_comp_inhibitors(m))
