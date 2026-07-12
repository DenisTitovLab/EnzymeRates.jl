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
Assert one `Step` conserves atoms: a binding step must move exactly the
bound metabolite's atoms onto the enzyme (`atoms(to) − atoms(from) ==
atoms(bound_metabolite)`); an iso step must leave the atom multiset
unchanged (`atoms(to) == atoms(from)`). Errors naming the offending step.
"""
function _assert_step_atom_conserving(reaction::EnzymeReaction, s::Step)
    diff = Dict{Symbol,Int}()
    _accumulate_atoms!(diff, _species_atoms(reaction, to_species(s)), 1)
    _accumulate_atoms!(diff, _species_atoms(reaction, from_species(s)), -1)
    diff = _nonzero_atoms(diff)
    bm = bound_metabolite(s)
    expected = bm === nothing ? Dict{Symbol,Int}() :
        _nonzero_atoms(copy(_atoms_dict(reaction, name(bm))))
    diff == expected || error(
        "atom-non-conserving step $(name(from_species(s))) → " *
        "$(name(to_species(s))) (bound " *
        "$(bm === nothing ? "—" : name(bm))): Δatoms $diff ≠ $expected")
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
                cur, new_species, Product(p), true)
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
    Symbol[name(bound_metabolite(s)) for s in path
           if is_binding(s) && bound_metabolite(s) isa Substrate]

"""
Product bound-metabolite names in route (path) order.
"""
_release_order(path::Vector{Step}) =
    Symbol[name(bound_metabolite(s)) for s in path
           if is_binding(s) && bound_metabolite(s) isa Product]

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
                    cur_species, new_species, Product(p), true)
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
                    cur_species, new_species, Substrate(s), true)
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
                        Substrate(s), true)
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
                            nothing, true)
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
                        nothing, true)
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
                    Substrate(s), true)
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
                        Substrate(s), true)
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
                            nothing, true)
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
                            nothing, true)
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
            is_binding(step) || continue
            bm = bound_metabolite(step)
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
            push!(result, Step[
                Step(from_species(s), to_species(s),
                     bound_metabolite(s), i != iso_idx)
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
                        base, _add(base, met), _role(met), true))
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
                        bound_metabolite(s), is_equilibrium(s)))
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
                push!(steps, Step(sp1, sp2, _role(met), true))
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
Reassign kinetic-group ids so binding steps sharing `(metabolite, RE/SS)`
collapse into one group. Each multi-step class is assigned a fresh id;
singleton classes and iso steps keep their existing id. Operates on the
`(steps, groups)` parallel-array form and returns the merged pair.
"""
function _apply_equivalence_grouping(
    steps::Vector{Step}, groups::Vector{Int},
)
    classes = Dict{Tuple{Symbol,Bool}, Vector{Int}}()
    for (i, s) in enumerate(steps)
        bm = bound_metabolite(s)
        bm === nothing && continue
        push!(get!(classes, (name(bm), is_equilibrium(s)), Int[]), i)
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


# ─── Expansion-Move Helpers ──────────────────────────────────

# ─── Expansion Moves ─────────────────────────────────────────

"""
Inhibitor-free core of a step: the same catalytic binding with every
`Regulator` stripped from both its species. Two steps with equal cores are
the same binding in different inhibitor contexts (mirrors) — e.g.
`E→E·Pyruvate` and `E·Pyruvateinh→E·Pyruvate·Pyruvateinh`.
"""
function _step_core(s::Step)
    strip(sp) = Species(
        Metabolite[b for b in bound(sp) if !(b isa Regulator)],
        conformation(sp), residual(sp))
    (strip(from_species(s)), strip(to_species(s)), bound_metabolite(s))
end

"""
Partition the RE→SS-eligible kinetic groups into mirror classes: connected
components of the graph where two groups are linked if they share a step core.
Eligible = all-RE and not an inhibitor binding (invariant 1). A catalytic
binding and its inhibitor-bound mirror, once a split has separated them into
different groups, land in one class and flip together (invariant 2). Same-group
mirrors and non-mirror groups each form their own singleton class, so behavior
is unchanged except where a split has separated a mirror. Classes are returned
sorted by lowest group index for deterministic move order.
"""
function _re_to_ss_flip_units(m::Union{Mechanism, AllostericMechanism})
    elig = [g for g in kinetic_groups(m)
            if all(is_equilibrium, steps(m)[g]) &&
               !any(s -> bound_metabolite(s) isa Regulator, steps(m)[g])]
    core_groups = Dict{Any, Vector{Int}}()
    for g in elig, s in steps(m)[g]
        push!(get!(core_groups, _step_core(s), Int[]), g)
    end
    parent = Dict(g => g for g in elig)
    root(x) = parent[x] == x ? x : root(parent[x])
    for gs in values(core_groups), i in 2:length(gs)
        parent[root(gs[i])] = root(gs[1])
    end
    comps = Dict{Int, Vector{Int}}()
    for g in elig
        push!(get!(comps, root(g), Int[]), g)
    end
    sort([sort(unique(c)) for c in values(comps)]; by = first)
end

"""
    _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})

Mechanism-native overload of the RE→SS expansion move. For each mirror class of
all-RE catalytic kinetic groups (`_re_to_ss_flip_units`), produce a variant with
every group in that class flipped to SS at once. Competitive-inhibitor bindings
are never flipped (RE-only), and a catalytic step flips together with its
inhibitor-bound mirror. All other groups, the reaction, and (for allosteric) the
catalytic-allo tags, multiplicity, and regulatory sites are preserved verbatim.
"""
function _expand_re_to_ss(m::Union{Mechanism, AllostericMechanism})
    results = typeof(m)[]
    for unit in _re_to_ss_flip_units(m)
        new_groups = steps(m)
        for g in unit
            new_groups = _flip_group_to_ss(new_groups, g)
        end
        push!(results, _with_steps(m, new_groups))
    end
    results
end

"""
Return a fresh `Vector{Vector{Step}}` matching `groups` but with every
Step in group `g` rebuilt with `is_equilibrium=false`. All other groups
are reused by reference (Step is immutable).
"""
function _flip_group_to_ss(groups::Vector{Vector{Step}}, g::Int)
    new_groups = Vector{Vector{Step}}()
    for (gi, gr) in enumerate(groups)
        if gi == g
            flipped = Step[
                Step(from_species(s), to_species(s),
                     bound_metabolite(s), false)
                for s in gr]
            push!(new_groups, flipped)
        else
            push!(new_groups, gr)
        end
    end
    new_groups
end

"""
    _expand_split_kinetic_group(m::Mechanism) → Vector{Mechanism}
    _expand_split_kinetic_group(am::AllostericMechanism) → Vector{AllostericMechanism}

Mechanism-native overload of the kinetic-group split move. For each
catalytic group with 2+ members, produce one variant per member where
that member is carved out into a fresh trailing group. The reaction
(and, for allosteric, multiplicity / regulatory sites) is preserved.
Catalytic allo-state tags are extended with the parent group's tag
appended (splitting is a parameter-relaxation move that MUST NOT
change A/I semantics).

Splitting a group adds a parameter, but a Wegscheider cycle often forces
that new parameter straight back to dependent, making the split a
model-space no-op that fits the parent's equation. `_canonical_mechanism`
merges such splits back, so each candidate is canonicalized and dropped
when it returns to the parent. These no-op splits dominate — up to ~2/3 of
the mechanisms in bi-bi enumeration — and, because every enumerated
mechanism is canonical, an un-dropped one would re-enter the beam as a
self-loop.
"""
function _expand_split_kinetic_group(m::Mechanism)
    results = Mechanism[]
    mc = _canonical_mechanism(m)
    for g in kinetic_groups(m)
        length(steps(m)[g]) >= 2 || continue
        for split_idx in 1:length(steps(m)[g])
            child = _canonical_mechanism(
                _with_steps(m, _split_one_step(steps(m), g, split_idx)))
            child == mc || push!(results, child)
        end
    end
    results
end

function _expand_split_kinetic_group(am::AllostericMechanism)
    results = AllostericMechanism[]
    mc = _canonical_mechanism(am)
    for g in kinetic_groups(am)
        length(steps(am)[g]) >= 2 || continue
        for split_idx in 1:length(steps(am)[g])
            new_groups = _split_one_step(steps(am), g, split_idx)
            new_states = vcat(cat_allo_states(am), [cat_allo_states(am)[g]])
            child = _canonical_mechanism(
                _with_steps_and_cat_states(am, new_groups, new_states))
            child == mc || push!(results, child)
        end
    end
    results
end

"""
Canonical kinetic-group partition: merge kinetic groups whose binding-K
representatives are single-symbol Wegscheider-tied (the relation
`_build_wegscheider_rename_map` finds), so split and merged encodings of the
same rate-equivalent graph collapse to one partition. Returns `mech.steps`
unchanged when nothing is tied.
"""
function _merge_tied_kinetic_groups(mech::Mechanism)
    rename = _build_wegscheider_rename_map(mech)
    isempty(rename) && return mech.steps
    groups = mech.steps
    step_params = _step_parameters(mech)
    # Canonical representative binding-K name per group (nothing if no binding step).
    rep = Vector{Union{Symbol, Nothing}}(nothing, length(groups))
    for (idx, (s, g)) in enumerate(_flat_steps(mech))
        is_equilibrium(s) && is_binding(s) || continue
        k = name(step_params[idx][1], mech)
        rep[g] = get(rename, k, k)
    end
    byrep = Dict{Symbol, Vector{Int}}()
    for (gi, r) in enumerate(rep)
        r === nothing && continue
        push!(get!(byrep, r, Int[]), gi)
    end
    any(length(v) > 1 for v in values(byrep)) || return mech.steps
    merged = Vector{Vector{Step}}()
    done = falses(length(groups))
    for gi in eachindex(groups)
        done[gi] && continue
        r = rep[gi]
        if r !== nothing && length(byrep[r]) > 1
            gis = byrep[r]
            push!(merged, Step[s for j in gis for s in groups[j]])
            for j in gis
                done[j] = true
            end
        else
            push!(merged, copy(groups[gi]))
            done[gi] = true
        end
    end
    merged
end

"""
Allosteric partition merge: the AllostericMechanism analog of
`_merge_tied_kinetic_groups(::Mechanism)`. Ties come from the per-state binding-K
Wegscheider relations (`_state_wegscheider_rename_map` for `:A` and `:I`). Each
catalytic group is keyed by `(tag, folded-binding-K-name(s))`; the tag is part of
the key so groups that differ in allosteric state (and therefore are not
rate-equivalent) never merge. Returns the merged
`(cat_steps, cat_allo_states)` parallel vectors, unchanged when nothing is tied.
"""
function _merge_tied_kinetic_groups(am::AllostericMechanism)
    rename_A = _state_wegscheider_rename_map(am, :A)
    rename_I = _state_wegscheider_rename_map(am, :I)
    groups = steps(am)
    tags = cat_allo_states(am)
    (isempty(rename_A) && isempty(rename_I)) && return groups, tags
    fold(d, nm) = get(d, nm, nm)
    keyof = Vector{Any}(nothing, length(groups))
    for (g, grp) in enumerate(groups)
        tag = tags[g]
        bstep = nothing
        for s in grp
            if is_equilibrium(s) && is_binding(s)
                bstep = s
                break
            end
        end
        bstep === nothing && continue
        astate = tag === :EqualAI ? :EqualAI : :A
        aK = fold(rename_A, name(Kd(bstep, astate), am))
        keyof[g] = tag === :NonequalAI ?
            (:NonequalAI, aK, fold(rename_I, name(Kd(bstep, :I), am))) :
            (tag, aK)
    end
    bykey = Dict{Any, Vector{Int}}()
    for (g, k) in enumerate(keyof)
        k === nothing && continue
        push!(get!(bykey, k, Int[]), g)
    end
    any(length(v) > 1 for v in values(bykey)) || return groups, tags
    merged_steps = Vector{Vector{Step}}()
    merged_tags = Symbol[]
    done = falses(length(groups))
    for g in eachindex(groups)
        done[g] && continue
        k = keyof[g]
        if k !== nothing && length(bykey[k]) > 1
            gis = bykey[k]
            push!(merged_steps, Step[s for j in gis for s in groups[j]])
            push!(merged_tags, tags[g])
            for j in gis
                done[j] = true
            end
        else
            push!(merged_steps, copy(groups[g]))
            push!(merged_tags, tags[g])
            done[g] = true
        end
    end
    merged_steps, merged_tags
end

"""
Canonical form of a mechanism for deduplication: the same graph with its
kinetic-group partition merged over Wegscheider-tied binding-K's. Two
graph-distinct mechanisms that reduce to the same rate function collapse to the
same canonical mechanism, so their rendered equation and `eq_hash` agree.
Applied by the split-move expansion (`_expand_split_kinetic_group`), so every
enumerated mechanism is canonical.
"""
function _canonical_mechanism(m::Mechanism; max_passes::Int = 8)
    prev = m
    for _ in 1:max_passes   # convergence is ≤2 passes in practice
        merged = Mechanism(reaction(prev), _merge_tied_kinetic_groups(prev))
        merged == prev && return merged
        prev = merged
    end
    error("_canonical_mechanism did not reach a fixed point in $max_passes " *
          "merge passes — the kinetic-group merge is not converging, a " *
          "canonicalization bug for the mechanism producing this")
end

function _canonical_mechanism(am::AllostericMechanism; max_passes::Int = 8)
    prev = am
    for _ in 1:max_passes
        cat_steps, cat_states = _merge_tied_kinetic_groups(prev)
        merged = AllostericMechanism(reaction(prev), cat_steps, cat_states,
                                     catalytic_multiplicity(prev),
                                     copy(regulatory_sites(prev)))
        merged == prev && return merged
        prev = merged
    end
    error("_canonical_mechanism did not reach a fixed point in $max_passes " *
          "merge passes — the kinetic-group merge is not converging, a " *
          "canonicalization bug for the mechanism producing this")
end

"""
Return a fresh `Vector{Vector{Step}}` matching `groups` but with the
step at `(g, split_idx)` moved into a new trailing singleton group.
Other groups are reused by reference (Step / Vector{Step} are immutable
from this caller's perspective).
"""
function _split_one_step(
    groups::Vector{Vector{Step}}, g::Int, split_idx::Int,
)
    new_groups = Vector{Vector{Step}}()
    for (gi, gr) in enumerate(groups)
        if gi == g
            remaining = Step[gr[i] for i in eachindex(gr)
                             if i != split_idx]
            push!(new_groups, remaining)
        else
            push!(new_groups, gr)
        end
    end
    push!(new_groups, Step[groups[g][split_idx]])
    new_groups
end

"""
Map each form name in `m`'s step graph to the set of bound-metabolite
names read directly from its `Species`. Returns
`Dict{form_name => Set{met_name}}`.
"""
function _bound_at_forms(m::Union{Mechanism, AllostericMechanism})
    result = Dict{Symbol, Set{Symbol}}()
    for group in steps(m), s in group
        for sp in (from_species(s), to_species(s))
            fn = name(sp)
            haskey(result, fn) ||
                (result[fn] = Set(name(b) for b in bound(sp)))
        end
    end
    result
end

"""
Return source-form names that have a binding step for the named
metabolite. The source form is the side without the metabolite (RE
binding canonicalizes the metabolite onto `to_species`).
"""
function _forms_with_binding_step_native(
    m::Union{Mechanism, AllostericMechanism}, met_name::Symbol,
)
    result = Set{Symbol}()
    for group in steps(m), s in group
        bm = bound_metabolite(s)
        bm === nothing && continue
        name(bm) == met_name || continue
        push!(result, name(from_species(s)))
    end
    result
end

"""
Extend `rxn`'s regulators with a new `CompetitiveInhibitor(name)`,
preserving every other field (reactants, allowed catalytic
multiplicities). `EnzymeReaction`'s inner constructor canonicalizes the
regulator order.
"""
function _add_competitive_inhibitor(rxn::EnzymeReaction, reg_name::Symbol)
    any(rm -> name(regulator(rm)) == reg_name, regulators(rxn)) && return rxn
    new_regs = copy(regulators(rxn))
    push!(new_regs, RegulatorMults(CompetitiveInhibitor(reg_name), Int[1]))
    EnzymeReaction(copy(reactants(rxn)), new_regs,
                   copy(allowed_catalytic_multiplicities(rxn)))
end

"""
    _expand_add_dead_end_regulator(m::Mechanism, rxn::EnzymeReaction;
                                   exclude_regs) → Vector{Mechanism}
    _expand_add_dead_end_regulator(am::AllostericMechanism,
                                   rxn::EnzymeReaction; exclude_regs) →
        Vector{AllostericMechanism}

Add a dead-end regulator binding step set. For each `CompetitiveInhibitor`
declared in `rxn` but not yet bound by `m`'s steps, enumerate inhibitor
competition patterns (S × P × existing inhibitors); for each pattern,
add RE binding steps to forms where the competing metabolite has a
binding step and the form isn't already bound by any competing
metabolite. Mirror steps inherit their catalytic counterpart's
`kinetic_group`. All new binding steps for a single regulator share one
fresh trailing kinetic group (one K_R parameter).

The caller must pass the declared `rxn` because `m.reaction` only
carries regulators already bound by its steps; not-yet-bound regulators
live exclusively in the declared reaction. The new `Mechanism`'s
reaction is `rxn` extended with the newly-bound regulator (preserving
the substrate / product / multiplicity payload).
"""
function _expand_add_dead_end_regulator(
    m::Mechanism, rxn::EnzymeReaction;
    exclude_regs::Set{Symbol}=Set{Symbol}(),
)
    raw = _expand_add_dead_end_regulator_native(
        m, rxn, Set{Symbol}();
        exclude_regs=exclude_regs)
    [Mechanism(new_reaction, new_groups)
     for (new_groups, _, new_reaction) in raw]
end

function _expand_add_dead_end_regulator(
    am::AllostericMechanism, rxn::EnzymeReaction;
    exclude_regs::Set{Symbol}=Set{Symbol}(),
)
    # Allosteric ligands are excluded — dead-end is not the right move
    # for them (they belong on a regulatory site).
    allo_ligands = Set{Symbol}()
    for site in regulatory_sites(am), lig in ligands(site)
        push!(allo_ligands, name(lig))
    end
    raw = _expand_add_dead_end_regulator_native(
        am, rxn, allo_ligands;
        exclude_regs=exclude_regs)
    [AllostericMechanism(new_reaction, new_groups,
                         vcat(cat_allo_states(am), [:EqualAI]),
                         catalytic_multiplicity(am),
                         copy(regulatory_sites(am)))
     for (new_groups, _, new_reaction) in raw]
end

"""
Shared kernel for the Mechanism / AllostericMechanism dead-end
expansion. `additional_excluded` carries allosteric ligand names that
should not be eligible (empty for plain `Mechanism`). Returns raw
`(new_groups, new_regulator_group_index, new_reaction)` tuples; the
top-level methods construct the typed result.
"""
function _expand_add_dead_end_regulator_native(
    m::Union{Mechanism, AllostericMechanism},
    rxn::EnzymeReaction,
    additional_excluded::Set{Symbol};
    exclude_regs::Set{Symbol},
)
    isempty(regulators(rxn)) &&
        return Tuple{Vector{Vector{Step}}, Int, EnzymeReaction}[]

    sub_names = Set(name(s) for s in substrates(rxn))
    prod_names = Set(name(p) for p in products(rxn))

    existing_regs = Set{Symbol}()
    for group in steps(m), s in group
        bm = bound_metabolite(s)
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
    isempty(eligible_regs) &&
        return Tuple{Vector{Vector{Step}}, Int, EnzymeReaction}[]

    form_sp = Dict{Symbol, Species}()
    for group in steps(m), s in group
        form_sp[name(from_species(s))] = from_species(s)
        form_sp[name(to_species(s))] = to_species(s)
    end
    cat_forms = Set(keys(form_sp))

    n_groups_before = length(steps(m))
    results = Tuple{Vector{Vector{Step}}, Int, EnzymeReaction}[]

    boundmap = _bound_at_forms(m)

    for reg_name in eligible_regs
        eligible_forms = Symbol[]
        for f in sort(collect(cat_forms))
            haskey(boundmap, f) || continue
            fb = boundmap[f]
            (intersect(fb, sub_names) == sub_names ||
                intersect(fb, prod_names) == prod_names) && continue
            push!(eligible_forms, f)
        end
        isempty(eligible_forms) && continue

        existing_inhibitors = Symbol[]
        for group in steps(m), s in group
            bm = bound_metabolite(s)
            bm isa Regulator || continue
            name(bm) == reg_name && continue
            push!(existing_inhibitors, name(bm))
        end
        sort!(unique!(existing_inhibitors))

        inh_patterns = _inhibitor_competition_patterns(
            sub_names, prod_names, existing_inhibitors)
        seen = Set{Vector{Symbol}}()

        for (comp_subs, comp_prods, comp_inhibitors) in inh_patterns
            target_forms = Set{Symbol}()
            for met in comp_subs
                union!(target_forms,
                       _forms_with_binding_step_native(m, met))
            end
            for met in comp_prods
                union!(target_forms,
                       _forms_with_binding_step_native(m, met))
            end
            for inh in comp_inhibitors
                union!(target_forms,
                       _forms_with_binding_step_native(m, inh))
            end

            all_competing = union(comp_subs, comp_prods, comp_inhibitors)
            active = Symbol[]
            for f in sort(collect(target_forms))
                f in eligible_forms || continue
                haskey(boundmap, f) || continue
                isempty(intersect(boundmap[f], all_competing)) || continue
                push!(active, f)
            end
            isempty(active) && continue
            active in seen && continue
            push!(seen, active)

            de_species_map = Dict{Symbol, Species}()
            reg_group_steps = Step[]
            for cf in active
                base = form_sp[cf]
                de_species = Species(
                    Metabolite[bound(base)..., CompetitiveInhibitor(reg_name)],
                    conformation(base), residual(base))
                de_species_map[cf] = de_species
                push!(reg_group_steps, Step(
                    base, de_species,
                    CompetitiveInhibitor(reg_name), true))
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
                             bound_metabolite(s), is_equilibrium(s)))
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

            new_reaction = _add_competitive_inhibitor(rxn, reg_name)
            push!(results, (new_groups, n_groups_before + 1, new_reaction))
        end
    end
    results
end

"""
    _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

Mechanism-native overload: convert a non-allosteric `Mechanism` into
allosteric variants, keeping only variants that are empirically
distinguishable from a simpler mechanism (an MWC conformational
constant `L` that has no observable effect is not enumerated):

  * The all-`:EqualAI` baseline is never emitted — the two conformations
    are identical, `L` cancels, and the mechanism is indistinguishable
    from `m`.
  * A binding group (its representative step binds a metabolite) set to
    `:OnlyA` is emitted bare: the bound metabolite's concentration
    reveals `L` (a K-type mechanism).
  * A catalytic group (its representative step is an isomerization, no
    bound metabolite) set to `:OnlyA` is emitted ONLY paired with a
    declared allosteric regulator at a new site, one variant per
    `(regulator, tag)` with `tag ∈ {:OnlyA, :OnlyI}` (a V-type
    mechanism). With no regulator bound, the inactive state binds
    substrate/product identically to the active state but cannot
    catalyze, so `L` folds entirely into `kcat`
    (`v = kcat/(1+L)·shape`) and is not observable; a reaction with no
    declared allosteric regulators emits nothing for that group.

For each value in `rxn`'s `allowed_catalytic_multiplicities`, the
multiplicity becomes the variant's `catalytic_multiplicity`. Catalytic
steps are reused by reference.
"""
function _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
    n_g = length(steps(m))
    base_tags = Symbol[:EqualAI for _ in 1:n_g]
    regs = Symbol[]
    for rm in regulators(rxn)
        reg = regulator(rm)
        reg isa AllostericRegulator && push!(regs, name(reg))
    end
    sort!(regs)
    results = AllostericMechanism[]
    for cn in allowed_catalytic_multiplicities(rxn)
        for g in 1:n_g
            new_tags = copy(base_tags)
            new_tags[g] = :OnlyA
            if is_iso(rep_step(m, g))
                am_cat = AllostericMechanism(
                    reaction(m), copy(steps(m)), new_tags, cn, RegulatorySite[])
                for reg in regs, tag in (:OnlyA, :OnlyI)
                    push!(results, _make_am_with_added_reg(am_cat, reg, tag, 0))
                end
            else
                push!(results, AllostericMechanism(
                    reaction(m), copy(steps(m)), new_tags, cn, RegulatorySite[]))
            end
        end
    end
    results
end

"""
    _expand_to_allosteric(am::AllostericMechanism, rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

`AllostericMechanism` input is already allosteric — no-op move.
Only non-allosteric mechanisms can be promoted by this move.
"""
_expand_to_allosteric(::AllostericMechanism, ::EnzymeReaction) =
    AllostericMechanism[]

"""
    _expand_add_allosteric_regulator(am::AllostericMechanism,
                                     rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

Add one `AllostericRegulator` declared in `rxn` but not yet bound by
`am` (neither at a regulatory site nor as a catalytic-step dead-end
inhibitor). For each (new ligand, target site, tag) combination, emit
a variant:

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

    existing_de = Set{Symbol}()
    for group in steps(am), s in group
        bm = bound_metabolite(s)
        bm isa Regulator && push!(existing_de, name(bm))
    end

    new_regs = Symbol[]
    for rm in regulators(rxn)
        reg = regulator(rm)
        reg isa AllostericRegulator || continue
        name(reg) in existing_allo && continue
        name(reg) in existing_de && continue
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
    _expand_add_allosteric_regulator(m::Mechanism, rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

Non-allosteric input: no-op; this move only extends allosteric mechanisms.
The dispatch shape here ensures callers walking a mixed Mechanism /
AllostericMechanism collection don't need to type-check.
"""
_expand_add_allosteric_regulator(::Mechanism, ::EnzymeReaction) =
    AllostericMechanism[]

"""
    _expand_change_allo_state(am::AllostericMechanism)
        → Vector{AllostericMechanism}

Mechanism-native overload. Relax one allo state from a "constrained"
tag (`:EqualAI`, `:OnlyA`, `:OnlyI`) to `:NonequalAI`. Variants are
emitted for each catalytic kinetic group tag and each regulatory
ligand tag that is not already `:NonequalAI`. The base catalytic
steps, multiplicity, and untouched tags are preserved.
"""
function _expand_change_allo_state(am::AllostericMechanism)
    results = AllostericMechanism[]

    for g in 1:length(cat_allo_states(am))
        cat_allo_states(am)[g] == :NonequalAI && continue
        new_states = copy(cat_allo_states(am))
        new_states[g] = :NonequalAI
        push!(results, _with_cat_allo_states(am, new_states))
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

    results
end

"""
    _expand_change_allo_state(m::Mechanism) → Vector{AllostericMechanism}

Non-allosteric input: no-op; this move only relaxes allosteric state tags.
"""
_expand_change_allo_state(::Mechanism) =
    AllostericMechanism[]

"""
    _expand_merge_regulatory_sites(am::AllostericMechanism)
        → Vector{AllostericMechanism}

Merge each unordered pair of regulatory sites into one shared site holding
both sites' ligands, at no parameter cost (Δ0 — every ligand keeps its own
dissociation constant). This lets two regulators compete at one site,
including the activator↔antagonist ambiguity. For the merged ligands,
enumerate the Δ0-valid allo-state assignments:

  * the all-keep assignment — every ligand retains its current state
    (co-binding);
  * each assignment retagging exactly one ligand to `:EqualAI` — the
    antagonist forms (a ligand that binds both conformations equally,
    competing for the shared site).

The all-`:EqualAI` assignment (degenerate) is dropped by this move's own
guard; the `RegulatorySite` constructor does not reject it, so that guard is
load-bearing. The merged site reuses one site's
`multiplicity` (equal to `catalytic_multiplicity`) and its ligands are sorted
by name, so two merge routes reaching the same ligand partition produce `==`
mechanisms and dedup by `hash`. Sign is not enforced here; `expand_mechanisms`
runs every child through `_filter_by_sign`.
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
        for states in _merged_site_state_assignments(base_states)
            merged = RegulatorySite(copy(ligs), mult, states)
            push!(results, _with_reg_sites(am, vcat(others, [merged])))
        end
    end
    results
end

"""
    _expand_merge_regulatory_sites(::Mechanism) → Vector{AllostericMechanism}

Non-allosteric input: no-op; this move only merges regulatory sites, which a
`Mechanism` has none of. Keeps callers type-uniform.
"""
_expand_merge_regulatory_sites(::Mechanism) =
    AllostericMechanism[]

# Δ0-valid allo-state assignments for a merged site's ligands: the all-keep
# assignment, plus each assignment retagging exactly one non-`:EqualAI` ligand
# to `:EqualAI`. The all-`:EqualAI` result is dropped.
function _merged_site_state_assignments(base_states::Vector{Symbol})
    assignments = Vector{Symbol}[copy(base_states)]
    for i in eachindex(base_states)
        base_states[i] == :EqualAI && continue
        retagged = copy(base_states)
        retagged[i] = :EqualAI
        all(==(:EqualAI), retagged) && continue
        push!(assignments, retagged)
    end
    assignments
end

# ─── Regulator-Sign Filter ────────────────────────────────────

"""
    _regulator_sign(reg_name::Symbol, rxn::EnzymeReaction) -> Symbol

The declared sign (`:activator`, `:inhibitor`, or `:unspecified`) of the
`AllostericRegulator` named `reg_name` in `rxn`'s `regulators`. Returns
`:unspecified` when `reg_name` names no `AllostericRegulator` entry (either
absent entirely, or present only as some other `Regulator` subtype).
"""
function _regulator_sign(reg_name::Symbol, rxn::EnzymeReaction)
    for rm in regulators(rxn)
        reg = regulator(rm)
        reg isa AllostericRegulator && name(reg) == reg_name && return sign(rm)
    end
    :unspecified
end

"""
    _state_respects_sign(reg_name, state::Symbol, sibling_names::Vector{Symbol},
                          rxn::EnzymeReaction) -> Bool

Whether ligand `reg_name`'s allosteric `state` is consistent with its
declared regulator sign, given `sibling_names` — the OTHER ligands at its
regulatory site. `:unspecified` sign always passes. A designated activator
is never `:OnlyI` and a designated inhibitor is never `:OnlyA`; the
sign-matching pure state and `:NonequalAI` are always allowed. `:EqualAI`
is an antagonist state and is rejected only when a sibling is a
same-sign designated effector, since that would counteract the sibling's
declared direction.
"""
function _state_respects_sign(reg_name, state::Symbol,
                              sibling_names::Vector{Symbol}, rxn::EnzymeReaction)
    sgn = _regulator_sign(reg_name, rxn)
    sgn === :unspecified && return true
    sgn === :activator && state === :OnlyI && return false
    sgn === :inhibitor && state === :OnlyA && return false
    if state === :EqualAI
        any(s -> _regulator_sign(s, rxn) === sgn, sibling_names) && return false
    end
    true
end

"""
    _filter_by_sign(mechs::Vector, rxn::EnzymeReaction) -> Vector

Keep only mechanisms whose every regulatory ligand respects its declared
sign (`_state_respects_sign`) given its site's other ligands. A `Mechanism`
has no regulatory sites and passes trivially.
"""
_filter_by_sign(mechs::Vector, rxn::EnzymeReaction) =
    filter(m -> _respects_sign(m, rxn), mechs)

_respects_sign(::Mechanism, ::EnzymeReaction) = true
function _respects_sign(am::AllostericMechanism, rxn::EnzymeReaction)
    for site in regulatory_sites(am)
        lig_names = Symbol[name(lig) for lig in ligands(site)]
        for (i, lig_name) in enumerate(lig_names)
            siblings = Symbol[lig_names[j] for j in eachindex(lig_names) if j != i]
            _state_respects_sign(lig_name, allo_states(site)[i], siblings, rxn) ||
                return false
        end
    end
    true
end

"""
    expand_mechanisms(mechs, reaction) -> Vector{Union{Mechanism, AllostericMechanism}}

Apply all expansion moves (RE→SS, split kinetic group, add dead-end
regulator, to-allosteric, add allosteric regulator, change allo state, merge
regulatory sites) to each input mechanism and return the children as a flat
vector. Bucketing by parameter count is the caller's job, not enumeration's.
"""
function expand_mechanisms(
    mechs::Vector{<:Union{Mechanism, AllostericMechanism}},
    rxn::EnzymeReaction)
    result = Union{Mechanism, AllostericMechanism}[]
    for m in mechs
        _add_expansions_mech!(result, m, rxn)
    end
    result = _filter_by_sign(result, rxn)
    for child in result
        _assert_atom_conserving(child)
    end
    result
end

function _add_expansions_mech!(
    result::Vector{Union{Mechanism, AllostericMechanism}},
    m::Union{Mechanism, AllostericMechanism},
    rxn::EnzymeReaction)
    append!(result, _expand_re_to_ss(m))
    append!(result, _expand_split_kinetic_group(m))
    append!(result, _expand_add_dead_end_regulator(m, rxn))
    append!(result, _expand_to_allosteric(m, rxn))
    append!(result, _expand_add_allosteric_regulator(m, rxn))
    append!(result, _expand_change_allo_state(m))
    append!(result, _expand_merge_regulatory_sites(m))
end

# --- Dedup ---


"""
    init_mechanisms(reaction::EnzymeReaction) -> Vector{Mechanism}

Public entry point. Produces all mechanisms at minimum parameter count
for a reaction as concrete `Mechanism` structs. For each catalytic
topology (`_catalytic_topologies`): 1 SS step, all substrate/product
dead-end subsets (`_expand_substrate_product_dead_ends`), with binding
steps sharing the same `(metabolite, RE/SS)` class collapsed into one
kinetic group (`_apply_equivalence_grouping`).
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
    mechs
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
4. every ligand respects its declared sign (`_filter_by_sign`).

The seeds are the valid nodes that additionally bind ALL of `required_allo` at
regulatory sites and ALL of `required_comp` as competitive-inhibitor dead ends.
Returned deduped (each node is visited once).
"""
function seed_mechanisms(rxn::EnzymeReaction, required_allo::Set{Symbol},
                         required_comp::Set{Symbol})
    visited = Set{UInt64}()
    queue = Union{Mechanism, AllostericMechanism}[]
    seeds = Union{Mechanism, AllostericMechanism}[]
    enqueue!(m) = begin
        h = hash(m)
        h in visited && return
        push!(visited, h)
        push!(queue, m)
        _binds_all_required(m, required_allo, required_comp) && push!(seeds, m)
    end
    for m in init_mechanisms(rxn)
        enqueue!(m)
    end
    while !isempty(queue)
        m = popfirst!(queue)
        for c in _seed_children(m, rxn, required_allo)
            _is_seed_node(c, rxn, required_allo, required_comp) && enqueue!(c)
        end
    end
    seeds
end

# Children of `m` under the seed-build structure moves. The two
# allosteric-lifting moves run only when an allosteric regulator is required, so
# a competitive-only required set (`required_allo` empty) stays non-allosteric —
# no `L` — and seeds at `base + n_required_comp`. The dead-end move always runs.
# Each move is a no-op on the mechanism kind it does not apply to.
function _seed_children(m::Union{Mechanism, AllostericMechanism},
                        rxn::EnzymeReaction, required_allo::Set{Symbol})
    children = Union{Mechanism, AllostericMechanism}[]
    if !isempty(required_allo)
        append!(children, _expand_to_allosteric(m, rxn))
        append!(children, _expand_add_allosteric_regulator(m, rxn))
    end
    append!(children, _expand_add_dead_end_regulator(m, rxn))
    children
end

_has_nonequalai(::Mechanism) = false
_has_nonequalai(am::AllostericMechanism) =
    any(==(:NonequalAI), cat_allo_states(am)) ||
    any(site -> any(==(:NonequalAI), allo_states(site)), regulatory_sites(am))

_bound_allo_regs(::Mechanism) = Set{Symbol}()
_bound_allo_regs(am::AllostericMechanism) =
    Set{Symbol}(name(lig) for site in regulatory_sites(am) for lig in ligands(site))

function _bound_comp_inhibitors(m::Union{Mechanism, AllostericMechanism})
    bound = Set{Symbol}()
    for group in steps(m), s in group
        bm = bound_metabolite(s)
        bm isa CompetitiveInhibitor && push!(bound, name(bm))
    end
    bound
end

# A child worth expanding: cheap states, one ligand per regulatory site, no
# optional regulator bound, and sign-respecting.
function _is_seed_node(m::Union{Mechanism, AllostericMechanism},
                       rxn::EnzymeReaction, required_allo::Set{Symbol},
                       required_comp::Set{Symbol})
    _has_nonequalai(m) && return false
    m isa AllostericMechanism &&
        !all(site -> length(ligands(site)) == 1, regulatory_sites(m)) &&
        return false
    issubset(_bound_allo_regs(m), required_allo) || return false
    issubset(_bound_comp_inhibitors(m), required_comp) || return false
    _respects_sign(m, rxn) || return false
    true
end

_binds_all_required(m::Union{Mechanism, AllostericMechanism},
                    required_allo::Set{Symbol}, required_comp::Set{Symbol}) =
    issubset(required_allo, _bound_allo_regs(m)) &&
    issubset(required_comp, _bound_comp_inhibitors(m))

"""
    _assert_mechanism_invariants(m::Mechanism) -> Nothing

Structural invariants every valid Mechanism should satisfy:
- Every group is non-empty
- Each binding step's bound_metabolite is non-nothing AND iso steps have nothing
- from_species != to_species for every step
"""
function _assert_mechanism_invariants(m::Mechanism)
    flat = collect(Iterators.flatten(steps(m)))
    isempty(flat) && error("empty steps in Mechanism")
    for g in steps(m)
        isempty(g) && error("empty kinetic group in Mechanism")
    end
    for s in flat
        if is_binding(s)
            bound_metabolite(s) === nothing &&
                error("binding step has nothing bound_metabolite")
        else
            bound_metabolite(s) === nothing ||
                error("iso step has non-nothing bound_metabolite")
        end
        from_species(s) == to_species(s) &&
            error("from_species == to_species in step $s")
    end

    # Every declared substrate/product must appear in some step. Regulators
    # are excluded — init_mechanisms declares a dead-end inhibitor that no
    # step binds yet (expand_mechanisms binds it later; _drop_unbound_regulators
    # drops it at compile time). Substrates/products are never dropped.
    appearing = Set{Symbol}()
    for s in flat
        for sp in (from_species(s), to_species(s))
            for met in bound(sp)
                push!(appearing, name(met))
            end
        end
        bm = bound_metabolite(s)
        bm === nothing || push!(appearing, name(bm))
    end
    for met in (substrates(reaction(m))..., products(reaction(m))...)
        name(met) in appearing ||
            error("declared substrate/product $(name(met)) appears in no step")
    end

    # A kinetic group of size > 1 must bind a single metabolite with a single
    # RE/SS kind (no mixing). Iso steps within such a group are ignored here;
    # the per-step loop above already enforces bound/iso consistency.
    for group in steps(m)
        length(group) == 1 && continue
        kinds = [(is_equilibrium(s), bound_metabolite(s)) for s in group
                 if bound_metabolite(s) !== nothing]
        isempty(kinds) && continue
        first_eq, first_met = kinds[1]
        for (eq, met) in kinds[2:end]
            eq == first_eq ||
                error("kinetic group mixes RE and SS binding steps")
            met == first_met ||
                error("kinetic group binds different metabolites: " *
                      "$(name(first_met)) and $(name(met))")
        end
    end
    nothing
end

function _assert_mechanism_invariants(m::AllostericMechanism)
    # Check the base catalytic-side invariants (every cat-group non-empty,
    # etc.), then the allosteric-specific invariants against the actual
    # AllostericMechanism fields.
    flat = Step[s for g in steps(m) for s in g]
    isempty(flat) && error("AllostericMechanism: empty cat_steps")
    for g in steps(m)
        isempty(g) && error("AllostericMechanism: empty catalytic kinetic group")
    end

    # cat_allo_states is one per cat group, validated by the constructor;
    # re-check defensively here:
    length(cat_allo_states(m)) == length(steps(m)) ||
        error("AllostericMechanism: cat_allo_states length " *
              "$(length(cat_allo_states(m))) ≠ cat_steps length " *
              "$(length(steps(m)))")
    valid_cat_states = (:OnlyA, :EqualAI, :NonequalAI)
    for tag in cat_allo_states(m)
        tag in valid_cat_states ||
            error("AllostericMechanism: invalid cat allo state $tag")
    end

    catalytic_multiplicity(m) ≥ 1 ||
        error("AllostericMechanism: catalytic_multiplicity " *
              "$(catalytic_multiplicity(m)) must be ≥ 1")

    # regulatory_sites is a Vector{RegulatorySite}; each site carries
    # its own ligand list + multiplicity + per-ligand allo states. The
    # constructor validates internal structure; here we only assert the
    # list is non-nothing.
    regulatory_sites(m) isa Vector{RegulatorySite} ||
        error("AllostericMechanism: regulatory_sites not Vector{RegulatorySite}")

    nothing
end
