# ABOUTME: Mechanism enumeration by incremental parameter count growth
# ABOUTME: Provides init_mechanisms, expand_mechanisms, dedup! building blocks

# ─── Mechanism Spec Types ──────────────────────────────────────

abstract type AbstractMechanismSpec end

"""
Elementary step with kinetic-group tag. Steps sharing a `kinetic_group`
share kinetic parameters (one `K` for RE groups, one `kf`/`kr` for SS).
Canonical binding direction: metabolite on LHS for RE binding steps.
"""
struct StepSpec
    reactants::Vector{Symbol}   # [:E, :S] or [:EAB]
    products::Vector{Symbol}    # [:ES] or [:EPQ]
    is_equilibrium::Bool
    kinetic_group::Int
end

Base.:(==)(a::StepSpec, b::StepSpec) =
    a.reactants == b.reactants &&
    a.products == b.products &&
    a.is_equilibrium == b.is_equilibrium &&
    a.kinetic_group == b.kinetic_group

Base.hash(s::StepSpec, h::UInt) =
    hash(s.kinetic_group, hash(s.is_equilibrium,
        hash(s.products, hash(s.reactants, h))))

"""
    MechanismSpec <: AbstractMechanismSpec

Monomeric enzyme mechanism specification. Steps are `StepSpec` values
with inline form names, equilibrium status, and kinetic-group tag.
Same-group steps share kinetic parameters.
"""
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{StepSpec}
    n_fit_params_estimate::Int
end

"""
    AllostericMechanismSpec <: AbstractMechanismSpec

Allosteric enzyme mechanism built from a base `MechanismSpec`. Each
catalytic kinetic group and each regulator-site ligand carries an
allosteric state tag indicating its R/T-state relationship. Tags:
`:OnlyR`, `:OnlyT`, `:EqualRT`, `:NonequalRT`.

`group_tags` and `reg_ligand_tags` use **dense** Dict storage — every
kinetic group present in `base.steps` has an explicit entry in
`group_tags`, and every ligand listed in `allosteric_reg_sites` has an
explicit entry in `reg_ligand_tags`. The default tag for newly
constructed specs is `:NonequalRT`, but it is stored explicitly rather
than implied by Dict absence. The constructor validates this density.

`group_tags` maps `kinetic_group → tag`.
`reg_ligand_tags` maps `ligand → tag`.
"""
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    group_tags::Dict{Int, Symbol}
    reg_ligand_tags::Dict{Symbol, Symbol}
    n_fit_params_estimate::Int

    function AllostericMechanismSpec(
        base::MechanismSpec,
        catalytic_n::Int,
        allosteric_reg_sites::Vector{Vector{Symbol}},
        allosteric_multiplicities::Vector{Int},
        group_tags::Dict{Int, Symbol},
        reg_ligand_tags::Dict{Symbol, Symbol},
        n_fit_params_estimate::Int,
    )
        used_groups = Set(s.kinetic_group for s in base.steps)
        for g in used_groups
            haskey(group_tags, g) || error(
                "AllostericMechanismSpec: kinetic_group $g present " *
                "in base.steps but missing from group_tags " *
                "(dense storage required).")
        end
        for site in allosteric_reg_sites
            for lig in site
                haskey(reg_ligand_tags, lig) || error(
                    "AllostericMechanismSpec: ligand $lig listed in " *
                    "allosteric_reg_sites but missing from " *
                    "reg_ligand_tags (dense storage required).")
            end
        end
        new(base, catalytic_n, allosteric_reg_sites,
            allosteric_multiplicities, group_tags,
            reg_ligand_tags, n_fit_params_estimate)
    end
end

# ─── StepSpec Helpers ──────────────────────────────────────────

"""Return the metabolite for a step, or nothing for isomerization."""
step_metabolite(s::StepSpec) =
    length(s.reactants) == 2 ? s.reactants[2] : nothing

"""Return (from_form, to_form) for a step."""
step_forms(s::StepSpec) = (s.reactants[1], s.products[1])

"""Collect all unique form names from steps."""
function all_form_names(spec::MechanismSpec)
    forms = Set{Symbol}()
    for s in spec.steps
        push!(forms, s.reactants[1])
        push!(forms, s.products[1])
    end
    forms
end

function all_form_names(steps::Vector{StepSpec})
    forms = Set{Symbol}()
    for s in steps
        push!(forms, s.reactants[1])
        push!(forms, s.products[1])
    end
    forms
end

"""
    _forms_with_binding_step(steps, metabolite)
        → Set{Symbol}

Return source forms that have a binding step for `metabolite`.
The source form is the one that doesn't contain the metabolite
(the reactant side of a binding step).
"""
function _forms_with_binding_step(
    steps::Vector{StepSpec}, metabolite::Symbol,
)
    result = Set{Symbol}()
    for s in steps
        step_metabolite(s) === metabolite || continue
        push!(result, s.reactants[1])
    end
    result
end

# ─── Stage 1: Catalytic Topologies ───────────────────────────

"""
Build a `Species` from sorted bound substrate / product names plus an
`is_estar` flag. Names map to `Substrate` / `Product` structs;
conformation is `:Estar` when `is_estar`, else `:E`. Bound list is
sorted by name (Species inner constructor enforces this).
"""
function _make_species(
    bound_subs::Vector{Symbol},
    bound_prods::Vector{Symbol},
    is_estar::Bool,
)
    mets = Metabolite[Substrate.(bound_subs)...,
                      Product.(bound_prods)...]
    Species(mets, is_estar ? :Estar : :E)
end

"""Extract atom counts as Dict{Symbol,Int} for a metabolite."""
function _atoms_dict(
    @nospecialize(reaction::EnzymeReactionLegacy),
    met::Symbol,
)
    result = Dict{Symbol,Int}()
    for (name, atoms) in substrates(reaction)
        name == met || continue
        for (a, c) in atoms
            result[a] = get(result, a, 0) + c
        end
        return result
    end
    for (name, atoms) in products(reaction)
        name == met || continue
        for (a, c) in atoms
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
Release products one at a time after a multi-product
isomerization, then continue backtracking.
`is_estar` controls conformation tagging (Estar prefix).
`has_residual_atoms` tracks whether residual atoms remain.
"""
function _release_products!(
    all_paths, backtrack!,
    iso_species::Species,
    residual_atoms::Dict{Symbol,Int},
    consumed_subs::Vector{Symbol},
    released_prods::Vector{Symbol},
    prod_subset::Vector{Symbol},
    prod_atoms::Dict{Symbol,Dict{Symbol,Int}},
    is_estar::Bool,
    has_residual_atoms::Bool,
    steps::Vector{Step},
)
    # Generate all release orderings of products
    function _release_recurse!(
        cur::Species,
        cur_atoms::Dict{Symbol,Int},
        unreleased::Vector{Symbol},
        rel_so_far::Vector{Symbol},
    )
        if isempty(unreleased)
            # All products released, continue
            backtrack!(
                cur, residual_atoms,
                consumed_subs, rel_so_far,
                Symbol[], Symbol[],
                is_estar, false, steps
            )
            return
        end
        for p in copy(unreleased)
            new_unreleased = filter(!=(p), unreleased)
            new_species = _make_species(
                Symbol[], new_unreleased, is_estar)
            rel_step = Step(
                cur, new_species, Product(p), true)
            push!(steps, rel_step)
            _release_recurse!(
                new_species,
                _subtract_atoms(
                    cur_atoms, prod_atoms[p]),
                new_unreleased,
                [rel_so_far; p],
            )
            pop!(steps)
        end
    end

    _release_recurse!(
        iso_species,
        _add_atoms(
            residual_atoms,
            reduce(_add_atoms,
                [prod_atoms[p] for p in prod_subset];
                init=Dict{Symbol,Int}())),
        collect(prod_subset),
        copy(released_prods),
    )
end

"""
    _catalytic_topologies(reaction) -> Vector{MechanismSpec}

Build catalytic cycle topologies by constructive backtracking.
Each topology is a set of steps forming one or more complete
catalytic cycles (E -> ... -> E).
"""
function _catalytic_topologies(
    @nospecialize(reaction::EnzymeReactionLegacy),
)
    subs = substrates(reaction)
    prods = products(reaction)
    sub_names = Symbol[s[1] for s in subs]
    prod_names = Symbol[p[1] for p in prods]

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
    # - has_residual: enzyme carries leftover atoms from
    #     ping-pong
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
        # has_residual: true when enzyme is in Estar
        # conformation (ping-pong intermediate), even if
        # no residual atoms remain (empty-residual pp).
        # Controls both conformation tagging and branch
        # selection.
        has_residual::Bool,
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
                    Symbol[], new_on_prods, false)
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

        if isempty(on_enzyme_subs) && !has_residual
            # Free enzyme: bind any remaining substrate
            for s in remaining_subs
                new_on = [on_enzyme_subs; s]
                new_consumed = [consumed_subs; s]
                new_species = _make_species(
                    new_on, Symbol[], false)
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
        elseif !isempty(on_enzyme_subs) && !has_residual
            # Substrates bound, no residual
            # Option 1: bind another substrate (C5)
            if length(on_enzyme_subs) < max_bound
                for s in remaining_subs
                    new_on = [on_enzyme_subs; s]
                    new_consumed = [consumed_subs; s]
                    new_species = _make_species(
                        new_on, Symbol[], false)
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
                        n_prods_eff = k + (
                            isempty(residual) ? 0 : 1)
                        # C6: iso size limit
                        length(on_enzyme_subs) > 3 &&
                            continue
                        n_prods_eff > 3 && continue
                        # C8: product-only iso form
                        iso_species = _make_species(
                            Symbol[],
                            collect(prod_subset),
                            true)
                        step = Step(
                            cur_species, iso_species,
                            nothing, true)
                        push!(steps, step)
                        # Release products one at a time
                        # Ping-pong from free enzyme
                        # always creates Estar
                        _release_products!(
                            all_paths, backtrack!,
                            iso_species, residual,
                            consumed_subs, released_prods,
                            prod_subset, prod_atoms,
                            true, !isempty(residual),
                            steps
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
                        copy(remaining_prods), false)
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
        elseif isempty(on_enzyme_subs) && has_residual
            # C7: Estar with no subs — only bind, no iso
            for s in remaining_subs
                new_on = [s]
                new_consumed = [consumed_subs; s]
                new_species = _make_species(
                    new_on, Symbol[], true)
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
        elseif !isempty(on_enzyme_subs) && has_residual
            # Residual + substrates bound
            # Option 1: bind another substrate (C5)
            if length(on_enzyme_subs) < max_bound
                for s in remaining_subs
                    new_on = [on_enzyme_subs; s]
                    new_consumed = [consumed_subs; s]
                    new_species = _make_species(
                        new_on, Symbol[], true)
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
                            false)
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
                        can_cont = has_more ||
                            !isempty(remaining_subs)
                        can_cont || continue
                        # C8: product-only iso form
                        iso_species = _make_species(
                            Symbol[],
                            collect(prod_subset),
                            has_more)
                        step = Step(
                            cur_species, iso_species,
                            nothing, true)
                        push!(steps, step)
                        _release_products!(
                            all_paths, backtrack!,
                            iso_species, residual_atoms,
                            consumed_subs,
                            released_prods,
                            prod_subset, prod_atoms,
                            has_more, has_more,
                            steps
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

    isempty(all_paths) && return MechanismSpec[]

    # Deduplicate paths by their structural step content.
    # `Step` equality / hash ignore source_idx and use
    # canonical direction, so this correctly identifies
    # equal step multi-sets across paths.
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

    # Select binding steps whose source species' bound mets
    # are reachable under the cumulative weak ordering.
    # `met_filter` selects which side (substrates / products)
    # this ordering covers; `bound_filter` restricts the
    # source species' bound multiset to that same side so
    # ping-pong substrate-bound iso products don't gate on
    # product orderings (or vice versa).
    function _steps_for_ordering(
        all_group_steps, ordering, bound_filter,
    )
        selected = Set{Step}()
        accessible = Set{Symbol}()
        for level in ordering
            allowed = union(accessible, Set(level))
            for met in level
                for step in all_group_steps
                    is_binding(step) || continue
                    bm = bound_metabolite(step)
                    name(bm) == met || continue
                    src = from_species(step)
                    src_mets = Symbol[
                        name(m) for m in bound(src)
                        if bound_filter(m)
                    ]
                    if all(m ∈ allowed for m in src_mets)
                        push!(selected, step)
                    end
                end
            end
            union!(accessible, level)
        end
        selected
    end

    # --- Build topologies ---
    sub_met_set = Set(sub_names)
    prod_met_set = Set(prod_names)
    is_sub(m::Metabolite) = name(m) in sub_met_set
    is_prod(m::Metabolite) = name(m) in prod_met_set

    # Iterate iso_groups in a deterministic order: prefer
    # smaller iso-step counts first (sequential before
    # ping-pong), then by sorted iso-step names. This stabilizes
    # the topology output order so test ordering invariants hold
    # — `Set{Step}` hashing is not value-stable across this
    # refactor's data-structure change, so unordered iteration
    # would reorder topologies relative to legacy.
    sorted_iso_pats = sort(collect(keys(iso_groups));
        by = pat -> (
            length(pat),
            sort([
                (string(name(from_species(s))),
                 string(name(to_species(s))))
                for s in pat])))

    result = MechanismSpec[]
    for iso_pat in sorted_iso_pats
        group_paths = iso_groups[iso_pat]
        all_group_steps = Set{Step}()
        for path in group_paths
            union!(all_group_steps, path)
        end

        # Always include all isomerization steps
        iso_steps_set = Set{Step}(
            s for s in all_group_steps if is_iso(s))

        sub_binding_mets = Set{Symbol}()
        prod_binding_mets = Set{Symbol}()
        for step in all_group_steps
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
            sub_keys = _steps_for_ordering(
                all_group_steps, sub_ord, is_sub)
            prod_keys = _steps_for_ordering(
                all_group_steps, prod_ord, is_prod)
            topo_keys = union(iso_steps_set, sub_keys,
                prod_keys)
            topo_keys ∈ seen_topos && continue
            push!(seen_topos, topo_keys)

            steps = sort(collect(topo_keys); by=s -> (
                is_iso(s) ? 1 : 0,
                string(name(from_species(s))),
                string(name(to_species(s))),
            ))

            push!(result, _mechanism_spec_from_steps(
                reaction, steps))
        end
    end
    result
end

"""
Convert a backtracker-produced `Vector{Step}` topology into a
`MechanismSpec`. Picks the first iso step as the (only) SS step
(`is_equilibrium=false`); all other steps are RE. Each step is
assigned its source position as its kinetic-group index. Computes
`n_fit_params_estimate` from the form / step counts.
"""
function _mechanism_spec_from_steps(
    reaction, steps::Vector{Step},
)
    iso_idx = findfirst(is_iso, steps)
    tagged = StepSpec[
        _stepspec_from_step(s, i != iso_idx, i)
        for (i, s) in enumerate(steps)
    ]

    form_names = Set{Symbol}()
    for s in steps
        push!(form_names, name(from_species(s)))
        push!(form_names, name(to_species(s)))
    end
    n_forms = length(form_names)
    n_steps = length(steps)
    n_cycles = n_steps - n_forms + 1
    n_re = count(s -> s.is_equilibrium, tagged)
    n_ss = n_steps - n_re
    n_fit_params_estimate = n_re + 2 * n_ss - n_cycles

    MechanismSpec(reaction, tagged, n_fit_params_estimate)
end

"""
Render a `Step` into the legacy `StepSpec` shape used by the rest
of the enumeration pipeline. Binding steps put the metabolite on
LHS (canonical binding direction). Iso steps render in the forward
catalytic direction: the species with more bound substrates is the
source; the species with more bound products is the destination.
This preserves the semantic substrate-bound → product-bound iso
direction even after `Step`'s lex-based iso canonicalization.
"""
function _stepspec_from_step(
    s::Step, is_equilibrium::Bool, kinetic_group::Int,
)
    if is_binding(s)
        from_name = name(from_species(s))
        to_name   = name(to_species(s))
        m_name    = name(bound_metabolite(s))
        StepSpec(
            [from_name, m_name], [to_name],
            is_equilibrium, kinetic_group)
    else
        n_subs_from = count(b -> b isa Substrate,
                            bound(from_species(s)))
        n_subs_to   = count(b -> b isa Substrate,
                            bound(to_species(s)))
        # Forward iso: substrate-rich side → product-rich side.
        # Ties broken by lex on species name for determinism.
        forward = (n_subs_from > n_subs_to) ||
                  (n_subs_from == n_subs_to &&
                   string(name(from_species(s))) <=
                   string(name(to_species(s))))
        a, b = forward ? (from_species(s), to_species(s)) :
                         (to_species(s), from_species(s))
        StepSpec(
            [name(a)], [name(b)],
            is_equilibrium, kinetic_group)
    end
end

# ─── Dead-End Helpers ────────────────────────────────────────

"""
    _bound_metabolites_at_forms(spec, reaction)

Map each form name to its set of bound metabolite names by
tracing binding steps from :E. Isomerization steps swap
all substrates for all products.
"""
function _bound_metabolites_at_forms(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReactionLegacy),
)
    # Collect all metabolite names (including
    # suffixed regulator dummies)
    met_set = Set{Symbol}()
    for (name, _) in substrates(reaction)
        push!(met_set, name)
    end
    for (name, _) in products(reaction)
        push!(met_set, name)
    end
    for s in spec.steps
        met = step_metabolite(s)
        met !== nothing && push!(met_set, met)
    end

    # Parse bound metabolites from form name
    # by greedily matching known metabolites
    function _parse_bound(form::Symbol)
        s = string(form)
        # Strip E_ or Estar_ prefix
        if s == "E" || s == "Estar"
            return Set{Symbol}()
        end
        body = startswith(s, "Estar_") ? s[7:end] :
               startswith(s, "E_") ? s[3:end] : s
        # Split by _ and reassemble multi-part
        # metabolite names (e.g., I1__reg)
        parts = split(body, "_")
        result = Set{Symbol}()
        i = 1
        while i <= length(parts)
            # Try longest match first
            matched = false
            for len in length(parts):-1:1
                i + len - 1 > length(parts) && continue
                candidate = Symbol(join(
                    parts[i:i+len-1], "_"))
                if candidate in met_set
                    push!(result, candidate)
                    i += len
                    matched = true
                    break
                end
            end
            if !matched
                # Unknown part — skip
                i += 1
            end
        end
        result
    end

    forms = all_form_names(spec)
    bound = Dict{Symbol, Set{Symbol}}()
    for f in forms
        bound[f] = _parse_bound(f)
    end
    bound
end

"""
    _dead_end_form_name(base_form, base_bound, added_met)

Create form name for a dead-end form: base form's bound
metabolites plus the added metabolite. Preserves the E/Estar
prefix of the base form.
"""
function _dead_end_form_name(
    base_form::Symbol,
    base_bound::Set{Symbol}, added_met::Symbol,
)
    all_mets = sort(collect(
        union(base_bound, Set([added_met]))))
    prefix = _is_estar_form(base_form) ? "Estar" : "E"
    Symbol(prefix * "_" * join(all_mets, "_"))
end

function _is_estar_form(form::Symbol)
    s = string(form)
    s == "Estar" || startswith(s, "Estar_")
end

"""
    _substrate_product_dead_end_opportunities(
        bound, cat_forms, sub_names, prod_names)

Find (form, metabolite) dead-end opportunities for
substrates and products. A dead-end is valid when:
- The form doesn't already bind all substrates or all
  products
- The metabolite isn't already bound at the form
- The resulting form isn't a catalytic form
- The result binds at least one substrate AND at least
  one product (mixed binding required)
- The result doesn't have all substrates or all products
"""
function _substrate_product_dead_end_opportunities(
    bound::Dict{Symbol, Set{Symbol}},
    cat_forms::Set{Symbol},
    sub_names::Set{Symbol},
    prod_names::Set{Symbol},
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
            de_name = _dead_end_form_name(f, fb, m)
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
    _expand_substrate_product_dead_ends(specs, reaction)
        -> Vector{MechanismSpec}

For each spec, enumerate substrate/product dead-end
form combinations. A dead-end form is created when a
substrate or product binds to a catalytic form where it
doesn't normally bind, subject to:
- The resulting form is not already a catalytic form
- The resulting form binds at least one substrate AND
  at least one product (mixed binding required)
- The resulting form doesn't have all substrates or
  all products
"""
function _expand_substrate_product_dead_ends(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReactionLegacy),
)
    sub_names = Set(s[1] for s in substrates(reaction))
    prod_names = Set(p[1] for p in products(reaction))
    all_mets = union(sub_names, prod_names)

    # Competition patterns depend only on the reaction,
    # not the topology — compute once.
    patterns = _competition_patterns(
        sub_names, prod_names)

    result = MechanismSpec[]
    for spec in specs
        bound = _bound_metabolites_at_forms(
            spec, reaction)
        cat_forms = all_form_names(spec)

        de_opportunities =
            _substrate_product_dead_end_opportunities(
                bound, cat_forms, sub_names,
                prod_names)

        # Deduplicate: multiple catalytic forms may
        # produce the same dead-end form. Group by
        # dead-end form name.
        de_forms = Dict{Symbol,
            Vector{Tuple{Symbol, Symbol}}}()
        for (f, m) in de_opportunities
            de_name = _dead_end_form_name(
                f, bound[f], m)
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
                bound[f], Set([m]))
        end

        # Look up the catalytic step's kinetic_group when
        # we want to attach a mirror step to the same group.
        cat_step_groups = Dict{Tuple{Symbol,Symbol}, Int}()
        for s in spec.steps
            cat_step_groups[step_forms(s)] = s.kinetic_group
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

            # Build new steps: original + dead-end
            new_steps = copy(spec.steps)
            next_g = maximum(
                s.kinetic_group for s in new_steps;
                init=0) + 1

            # Add binding steps for active dead-ends.
            # Each binding step is an equivalence-eligible
            # candidate, but at this stage every binding
            # step gets its own fresh group (init_mechanisms
            # later applies same-metabolite grouping).
            for de_name in sort(collect(active_de))
                entries = de_forms[de_name]
                for (cat_form, met) in entries
                    push!(new_steps, StepSpec(
                        [cat_form, met], [de_name],
                        true, next_g))
                    next_g += 1
                end
            end

            # Add mirror steps: for each catalytic
            # step, if both endpoints have dead-end
            # forms with the same metabolite, add a
            # parallel step. Mirror inherits RE/SS
            # AND the catalytic step's kinetic_group.
            for s in spec.steps
                from, to =
                    s.reactants[1], s.products[1]
                met = step_metabolite(s)
                for de_met in sort(collect(all_mets))
                    haskey(bound, from) || continue
                    haskey(bound, to) || continue
                    de_met in bound[from] && continue
                    de_met in bound[to] && continue
                    from_de = _dead_end_form_name(
                        from, bound[from], de_met)
                    to_de = _dead_end_form_name(
                        to, bound[to], de_met)
                    from_de in active_de || continue
                    to_de in active_de || continue
                    g = s.kinetic_group
                    if met !== nothing
                        push!(new_steps, StepSpec(
                            [from_de, met], [to_de],
                            s.is_equilibrium, g))
                    else
                        push!(new_steps, StepSpec(
                            [from_de], [to_de],
                            s.is_equilibrium, g))
                    end
                end
            end

            pc = _n_fit_params_estimate_from_steps(new_steps)

            push!(result, MechanismSpec(
                spec.reaction, new_steps, pc))
        end
    end
    result
end

# ─── Compilation ─────────────────────────────────────────────

"""
    compile_mechanism(spec::MechanismSpec)
    compile_mechanism(spec::AllostericMechanismSpec)
    compile_mechanism(m::Mechanism)
    compile_mechanism(am::AllostericMechanism)

Convert a `MechanismSpec` to an `EnzymeMechanism`, an
`AllostericMechanismSpec` to an `AllostericEnzymeMechanism`,
a `Mechanism` to an `EnzymeMechanism`, or an `AllostericMechanism`
to an `AllostericEnzymeMechanism`.
"""
compile_mechanism(spec::MechanismSpec) =
    EnzymeMechanism(spec)
compile_mechanism(spec::AllostericMechanismSpec) =
    AllostericEnzymeMechanism(spec)
compile_mechanism(m::Mechanism) = EnzymeMechanism(m)
compile_mechanism(am::AllostericMechanism) = AllostericEnzymeMechanism(am)

"""Strip a `__regN?` suffix from a Symbol; identity if absent."""
function _strip_reg_suffix(sym::Symbol)
    s = string(sym)
    m = match(r"^(.+)__reg\d*$", s)
    m === nothing ? sym : Symbol(m.captures[1])
end

"""
    EnzymeMechanism(spec::MechanismSpec) → EnzymeMechanism

Build the singleton `EnzymeMechanism` type. Metabolite occurrences
that carry a `__regN?` suffix are stripped to their bare reaction
name; form names matching the same pattern are stripped when no
collision would result.
"""
function EnzymeMechanism(
    spec::MechanismSpec;
    exclude_regs::Set{Symbol}=Set{Symbol}(),
)
    rxn = spec.reaction isa EnzymeReaction ?
          _to_legacy_reaction(spec.reaction) : spec.reaction
    subs = Tuple(s[1] for s in substrates(rxn))
    prods = Tuple(p[1] for p in products(rxn))
    # Allosteric regulators belong on regulatory sites, not in the
    # catalytic mechanism — auto-exclude them from the catalytic
    # `regulators` list so EnzymeMechanism can be built standalone.
    auto_exclude = Set{Symbol}(
        name for (name, role) in regulator_roles(rxn)
        if role === :allosteric)
    union!(auto_exclude, exclude_regs)
    # Build the set of names actually appearing on any step (after
    # stripping the __reg suffix used by enumeration internals).
    appears_in_steps = Set{Symbol}()
    for s in spec.steps
        for sym in Iterators.flatten((s.reactants, s.products))
            push!(appears_in_steps, _strip_reg_suffix(sym))
        end
    end
    regs = Tuple(r for r in regulators(rxn)
                 if r ∉ auto_exclude && r ∈ appears_in_steps)
    met_set = Set{Symbol}()
    for g in (subs, prods, regs); for n in g; push!(met_set, n); end; end

    # Suffixed regulator dummies (e.g. :I__reg) act as metabolites in
    # steps but compile down to their bare reaction name.
    suffixed = Set{Symbol}()
    for s in spec.steps
        for sym in Iterators.flatten((s.reactants, s.products))
            stripped = _strip_reg_suffix(sym)
            stripped == sym && continue
            stripped in met_set && push!(suffixed, sym)
        end
    end

    # Form name renaming: strip suffix when no collision; otherwise
    # keep the suffixed name.
    form_set = Set{Symbol}()
    for s in spec.steps
        for sym in Iterators.flatten((s.reactants, s.products))
            sym in met_set && continue
            sym in suffixed && continue
            push!(form_set, sym)
        end
    end
    form_name_map = Dict{Symbol,Symbol}()
    used = Set{Symbol}()
    for name in sort!(collect(form_set))
        candidate = _strip_reg_suffix(name)
        if candidate == name || candidate ∉ used
            form_name_map[name] = candidate
            push!(used, candidate)
        else
            form_name_map[name] = name
            push!(used, name)
        end
    end

    rxns = Tuple(
        (
            Tuple(get(form_name_map, x, _strip_reg_suffix(x))
                  for x in s.reactants),
            Tuple(get(form_name_map, x, _strip_reg_suffix(x))
                  for x in s.products),
            s.is_equilibrium,
            s.kinetic_group,
        )
        for s in spec.steps
    )
    EnzymeMechanism((subs, prods, regs), rxns)
end

# ─── Mechanism Enumeration ───────────────────────────────────

"""
    init_mechanisms(reaction::EnzymeReaction) -> Vector{Mechanism}
    init_mechanisms(reaction::EnzymeReactionLegacy) -> Vector{MechanismSpec}

Produce all mechanisms at minimum parameter count for a reaction.
For each catalytic topology: 1 SS step, all K's grouped equal
per metabolite, all substrate/product dead-end subsets.

Step kinetic groups are reassigned so that all binding steps
sharing the same `(metabolite, RE/SS)` class collapse into one
group; iso steps and uncollapsable bindings stay singletons.

The `EnzymeReaction` method returns concrete `Mechanism` structs;
the `EnzymeReactionLegacy` method is an internal entry point still
used by enumeration consumers (`expand_mechanisms`, `dedup!`) that
have not yet been switched to `Mechanism`.
"""
function init_mechanisms(
    @nospecialize(reaction::EnzymeReactionLegacy),
)
    topos = _catalytic_topologies(reaction)
    expanded = _expand_substrate_product_dead_ends(
        topos, reaction)

    n_s = length(substrates(reaction))
    n_p = length(products(reaction))
    floor_pc = n_s + n_p + 1

    result = MechanismSpec[]
    for spec in expanded
        push!(result, _apply_equivalence_grouping(
            spec, floor_pc))
    end
    result
end

"""
Reassign `kinetic_group` so steps sharing `(metabolite, RE/SS)`
collapse into one group. Other steps preserve their existing
group ids (singletons remain singletons).

`pc` is the maximum of the group-based count and `floor_pc`,
ensuring the upper-bound invariant `pc >= length(parameters(m))`
holds when dead-end mirror groups create non-catalytic cycles
that don't impose new Wegscheider constraints.
"""
function _apply_equivalence_grouping(
    spec::MechanismSpec, floor_pc::Int,
)
    groups = _max_equivalence_constraints(spec)
    new_steps = StepSpec[]
    for s in spec.steps
        met = step_metabolite(s)
        if met !== nothing &&
                haskey(groups, (met, s.is_equilibrium))
            g = groups[(met, s.is_equilibrium)]
            push!(new_steps, StepSpec(
                s.reactants, s.products,
                s.is_equilibrium, g))
        else
            push!(new_steps, s)
        end
    end
    pc = max(
        _n_fit_params_estimate_from_steps(new_steps),
        floor_pc)
    MechanismSpec(spec.reaction, new_steps, pc)
end

"""
Compute parameter-count estimate from a step list, counting kinetic
groups rather than individual steps. This is exact for mechanisms
whose cycles all impose independent Wegscheider constraints; dead-end
mirror cycles add 0 effective constraints, so this formula can
underestimate — callers should floor to a safe lower bound.
"""
function _n_fit_params_estimate_from_steps(steps::Vector{StepSpec})
    groups_re = Set{Int}()
    groups_ss = Set{Int}()
    for s in steps
        if s.is_equilibrium
            push!(groups_re, s.kinetic_group)
        else
            push!(groups_ss, s.kinetic_group)
        end
    end
    n_forms = length(all_form_names(steps))
    n_thermo = length(steps) - n_forms + 1
    length(groups_re) + 2 * length(groups_ss) - n_thermo
end

"""
Group catalytic-cycle binding steps that share the same
`(metabolite, is_equilibrium)` class; each multi-step class is
assigned a fresh kinetic-group integer. Returns a `Dict{(met,is_eq) → gnum}`
covering only classes with 2+ steps. Iso steps (no metabolite) and
singleton classes are not in the dict — callers leave their existing
group ids intact.
"""
function _max_equivalence_constraints(spec::MechanismSpec)
    classes = Dict{Tuple{Symbol,Bool}, Vector{Int}}()
    for (i, s) in enumerate(spec.steps)
        met = step_metabolite(s)
        met === nothing && continue
        key = (met, s.is_equilibrium)
        push!(get!(classes, key, Int[]), i)
    end

    used = Set{Int}(s.kinetic_group for s in spec.steps)
    next_g = isempty(used) ? 1 : maximum(used) + 1
    result = Dict{Tuple{Symbol,Bool}, Int}()
    for (key, idxs) in classes
        length(idxs) >= 2 || continue
        result[key] = next_g
        next_g += 1
    end
    result
end

# ─── Expansion-Move Helpers ──────────────────────────────────

"""Return the underlying `Vector{StepSpec}` for either spec type."""
_steps(s::MechanismSpec) = s.steps
_steps(s::AllostericMechanismSpec) = s.base.steps

_n_fit_params_estimate(s::AbstractMechanismSpec) = s.n_fit_params_estimate

"""
    _n_fit_params_estimate(m::Mechanism)
    _n_fit_params_estimate(am::AllostericMechanism)

Upper-bound estimate of the fittable parameter count for a `Mechanism`
or `AllostericMechanism`. Counts kinetic GROUPS (one K per RE group,
two k per SS group) and subtracts the number of independent thermodynamic
cycles (Haldane + Wegscheider) bound by the enzyme-form graph. The
allosteric variant adds the per-tag bookkeeping plus +1 for `L`.

This mirrors the spec-side `_n_fit_params_estimate_from_steps` formula
but reads from `Mechanism`'s structural fields directly.
"""
function _n_fit_params_estimate(m::Mechanism)
    n_steps = sum(length, m.steps; init = 0)
    n_re_groups = 0
    n_ss_groups = 0
    form_names = Set{Symbol}()
    for group in m.steps
        isempty(group) && continue
        is_re = is_equilibrium(first(group))
        is_re ? (n_re_groups += 1) : (n_ss_groups += 1)
        for s in group
            push!(form_names, name(from_species(s)))
            push!(form_names, name(to_species(s)))
        end
    end
    n_thermo = n_steps - length(form_names) + 1
    n_re_groups + 2 * n_ss_groups - n_thermo
end

function _n_fit_params_estimate(am::AllostericMechanism)
    base = _n_fit_params_estimate(Mechanism(am.reaction, am.cat_steps))
    # +1 for L (R/T equilibrium constant), plus per-tag bookkeeping:
    # :NonequalRT catalytic groups double; :NonequalRT regulator ligands
    # also double; per-site Kreg adds a K per ligand.
    tag_extra = 0
    for tag in am.cat_allo_states
        tag == :NonequalRT && (tag_extra += 1)
    end
    reg_extra = 0
    for site in am.regulatory_sites
        for (lig, st) in zip(site.ligands, site.allo_states)
            reg_extra += 1
            st == :NonequalRT && (reg_extra += 1)
        end
    end
    base + 1 + tag_extra + reg_extra
end

"""
Return a copy of `spec` with its steps replaced and estimate updated.
For `AllostericMechanismSpec`, all allosteric-side state is preserved.
"""
_with_steps(spec::MechanismSpec, new_steps, new_pc) =
    MechanismSpec(spec.reaction, new_steps, new_pc)

_with_steps(spec::AllostericMechanismSpec, new_steps, new_pc) =
    AllostericMechanismSpec(
        MechanismSpec(spec.base.reaction, new_steps, new_pc),
        spec.catalytic_n,
        deepcopy(spec.allosteric_reg_sites),
        copy(spec.allosteric_multiplicities),
        copy(spec.group_tags),
        copy(spec.reg_ligand_tags),
        new_pc,
    )

# ─── Expansion Moves ─────────────────────────────────────────

"""
    _expand_re_to_ss(spec) → Vector{typeof(spec)}

For each kinetic group whose members are all RE, convert the entire
group to SS (atomic per group). One K (1 param) becomes one (k_f, k_r)
pair (2 params); for `:NonequalRT` allosteric groups the T-state pair
doubles the delta. Each variant adds +1 (cheap tag) or +2 (NonequalRT).
"""
function _expand_re_to_ss(spec::AbstractMechanismSpec)
    results = typeof(spec)[]
    steps = _steps(spec)
    groups = Dict{Int, Vector{Int}}()
    for (i, s) in enumerate(steps)
        push!(get!(groups, s.kinetic_group, Int[]), i)
    end
    for (g, idxs) in groups
        all(steps[i].is_equilibrium for i in idxs) || continue
        new_steps = copy(steps)
        for i in idxs
            old = steps[i]
            new_steps[i] = StepSpec(
                old.reactants, old.products, false, g)
        end
        delta = _re_to_ss_delta(spec, g)
        push!(results, _with_steps(
            spec, new_steps,
            _n_fit_params_estimate(spec) + delta))
    end
    results
end

"""
ΔP for converting kinetic group `g` from RE to SS. +1 for plain
mechanisms or allosteric groups with a "cheap" tag (`:EqualRT`,
`:OnlyR`, `:OnlyT`); +2 for `:NonequalRT` (T-state K also splits
into kf_T, kr_T).
"""
_re_to_ss_delta(::MechanismSpec, ::Int) = 1
_re_to_ss_delta(spec::AllostericMechanismSpec, g::Int) =
    spec.group_tags[g] == :NonequalRT ? 2 : 1

"""
    _expand_re_to_ss(m::Mechanism) → Vector{Mechanism}
    _expand_re_to_ss(m::AllostericMechanism) → Vector{AllostericMechanism}

Mechanism-native overload of the RE→SS expansion move. For each
catalytic kinetic group whose members are all RE, produce a variant
with that entire group flipped to SS (atomic per group). All other
groups, the reaction, and (for allosteric) the catalytic-allo tags,
multiplicity, and regulatory sites are preserved verbatim. Each step's
`source_idx` is preserved so positional parameter naming is stable.
"""
function _expand_re_to_ss(m::Mechanism)
    results = Mechanism[]
    for g in kinetic_groups(m)
        all(is_equilibrium, m.steps[g]) || continue
        push!(results, Mechanism(reaction(m),
            _flip_group_to_ss(m.steps, g)))
    end
    results
end

function _expand_re_to_ss(m::AllostericMechanism)
    results = AllostericMechanism[]
    for g in kinetic_groups(m)
        all(is_equilibrium, m.cat_steps[g]) || continue
        push!(results, AllostericMechanism(
            reaction(m),
            _flip_group_to_ss(m.cat_steps, g),
            copy(m.cat_allo_states),
            m.catalytic_multiplicity,
            copy(m.regulatory_sites)))
    end
    results
end

"""
Return a fresh `Vector{Vector{Step}}` matching `groups` but with every
Step in group `g` rebuilt with `is_equilibrium=false`. All other groups
are reused by reference (Step is immutable). `source_idx` is preserved.
"""
function _flip_group_to_ss(groups::Vector{Vector{Step}}, g::Int)
    new_groups = Vector{Vector{Step}}()
    for (gi, gr) in enumerate(groups)
        if gi == g
            flipped = Step[
                Step(from_species(s), to_species(s),
                     bound_metabolite(s), false;
                     source_idx = s.source_idx)
                for s in gr]
            push!(new_groups, flipped)
        else
            push!(new_groups, gr)
        end
    end
    new_groups
end

"""
    _expand_split_kinetic_group(spec) → Vector{typeof(spec)}

For each kinetic group with 2+ members, split one step out into a
fresh group. Each split adds +1 (RE plain or RE cheap-tag) up to
+4 (SS `:NonequalRT`) to `n_fit_params_estimate`, doubling for SS
and again for `:NonequalRT` (T-state).
"""
function _expand_split_kinetic_group(spec::AbstractMechanismSpec)
    results = typeof(spec)[]
    steps = _steps(spec)
    groups = Dict{Int, Vector{Int}}()
    for (i, s) in enumerate(steps)
        push!(get!(groups, s.kinetic_group, Int[]), i)
    end
    used = isempty(steps) ? 0 :
        maximum(s.kinetic_group for s in steps)
    new_g = used + 1
    for (g, idxs) in groups
        length(idxs) >= 2 || continue
        for split_idx in idxs
            new_steps = copy(steps)
            old = steps[split_idx]
            new_steps[split_idx] = StepSpec(
                old.reactants, old.products,
                old.is_equilibrium, new_g)
            delta = _split_group_delta(
                spec, g, old.is_equilibrium)
            new_pc = _n_fit_params_estimate(spec) + delta
            push!(results,
                _split_with_steps(
                    spec, new_steps, new_pc, g, new_g))
        end
    end
    results
end

# Like `_with_steps` but inherits the parent group's allosteric
# tag onto the freshly-created split group. Splitting is a
# parameter-relaxation move: the new group must share R/T-state
# semantics with the parent it was carved out of.
_split_with_steps(
    spec::MechanismSpec, new_steps, new_pc, _g, _new_g,
) = _with_steps(spec, new_steps, new_pc)

function _split_with_steps(
    spec::AllostericMechanismSpec, new_steps, new_pc,
    g::Int, new_g::Int,
)
    new_tags = copy(spec.group_tags)
    new_tags[new_g] = new_tags[g]
    AllostericMechanismSpec(
        MechanismSpec(spec.base.reaction, new_steps, new_pc),
        spec.catalytic_n,
        deepcopy(spec.allosteric_reg_sites),
        copy(spec.allosteric_multiplicities),
        new_tags,
        copy(spec.reg_ligand_tags),
        new_pc,
    )
end

"""
ΔP for splitting one step out of kinetic group `g`. Base: +1 (RE)
or +2 (SS). For `:NonequalRT` allosteric groups: doubled (T-state
pair also splits).
"""
_split_group_delta(::MechanismSpec, ::Int, is_re::Bool) =
    is_re ? 1 : 2

function _split_group_delta(
    spec::AllostericMechanismSpec, g::Int, is_re::Bool,
)
    base = is_re ? 1 : 2
    tag = spec.group_tags[g]
    tag == :NonequalRT ? 2 * base : base
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
change R/T semantics). Each Step's `source_idx` is preserved.
"""
function _expand_split_kinetic_group(m::Mechanism)
    results = Mechanism[]
    for g in kinetic_groups(m)
        length(m.steps[g]) >= 2 || continue
        for split_idx in 1:length(m.steps[g])
            push!(results, Mechanism(reaction(m),
                _split_one_step(m.steps, g, split_idx)))
        end
    end
    results
end

function _expand_split_kinetic_group(am::AllostericMechanism)
    results = AllostericMechanism[]
    for g in kinetic_groups(am)
        length(am.cat_steps[g]) >= 2 || continue
        for split_idx in 1:length(am.cat_steps[g])
            new_groups = _split_one_step(am.cat_steps, g, split_idx)
            new_states = vcat(am.cat_allo_states,
                              [am.cat_allo_states[g]])
            push!(results, AllostericMechanism(
                reaction(am),
                new_groups,
                new_states,
                am.catalytic_multiplicity,
                copy(am.regulatory_sites)))
        end
    end
    results
end

"""
Return a fresh `Vector{Vector{Step}}` matching `groups` but with the
step at `(g, split_idx)` moved into a new trailing singleton group.
The split step's `source_idx` is preserved; other groups are reused
by reference (Step / Vector{Step} are immutable from this caller's
perspective).
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
    _expand_add_dead_end_regulator(spec, reaction; exclude_regs)
        → Vector{typeof(spec)}

Add a dead-end regulator binding step set. For each regulator not
yet present, enumerate inhibitor competition patterns (S × P × existing
inhibitors); for each pattern, add RE binding steps to forms where the
competing metabolite has a binding step and the form isn't already
bound by any competing metabolite. Mirror steps inherit their
catalytic counterpart's `kinetic_group`. All new binding steps for a
single regulator share one new kinetic group (one K_R parameter).

Each variant adds +1 to `n_fit_params_estimate`.
"""
function _expand_add_dead_end_regulator(
    spec::AbstractMechanismSpec,
    @nospecialize(reaction::EnzymeReactionLegacy);
    exclude_regs::Set{Symbol}=Set{Symbol}(),
)
    roles = regulator_roles(reaction)
    isempty(roles) && return typeof(spec)[]

    steps = _steps(spec)
    sub_names = Set(s[1] for s in substrates(reaction))
    prod_names = Set(p[1] for p in products(reaction))

    # Allosteric ligands (if any) are excluded — dead-end is
    # not the right move for them.
    allo_ligands = Set{Symbol}()
    if spec isa AllostericMechanismSpec
        for site in spec.allosteric_reg_sites
            for l in site
                push!(allo_ligands, l)
            end
        end
    end

    existing_mets = Set{Symbol}()
    for s in steps
        for sym in Iterators.flatten((s.reactants, s.products))
            push!(existing_mets, sym)
        end
    end

    eligible_regs = Symbol[]
    for (name, role) in roles
        (role == :unknown || role == :dead_end) || continue
        name in exclude_regs && continue
        name in allo_ligands && continue
        reg_prefix = string(name) * "__reg"
        already = any(
            contains(string(m), reg_prefix)
            for m in existing_mets)
        already && continue
        push!(eligible_regs, name)
    end
    sort!(eligible_regs)
    isempty(eligible_regs) && return typeof(spec)[]

    bound = _bound_metabolites_at_forms(
        _base_or_self(spec), reaction)
    cat_forms = all_form_names(steps)

    cat_step_groups = Dict{Tuple{Symbol,Symbol}, Int}()
    for s in steps
        cat_step_groups[step_forms(s)] = s.kinetic_group
    end

    used_g = isempty(steps) ? 0 :
        maximum(s.kinetic_group for s in steps)
    next_g = used_g + 1

    results = typeof(spec)[]

    for reg in eligible_regs
        dummy = Symbol(string(reg) * "__reg")

        eligible_forms = Symbol[]
        for f in sort(collect(cat_forms))
            haskey(bound, f) || continue
            fb = bound[f]
            (intersect(fb, sub_names) == sub_names ||
                intersect(fb, prod_names) == prod_names) &&
                continue
            push!(eligible_forms, f)
        end
        isempty(eligible_forms) && continue

        existing_inhibitors = Symbol[]
        for s in steps
            met = step_metabolite(s)
            met === nothing && continue
            s_met = string(met)
            contains(s_met, "__reg") || continue
            met === dummy && continue
            push!(existing_inhibitors, met)
        end
        sort!(unique!(existing_inhibitors))

        inh_patterns = _inhibitor_competition_patterns(
            sub_names, prod_names, existing_inhibitors)
        seen = Set{Vector{Symbol}}()

        for (comp_subs, comp_prods, comp_inhibitors) in inh_patterns
            target_forms = Set{Symbol}()
            for met in comp_subs
                union!(target_forms,
                    _forms_with_binding_step(steps, met))
            end
            for met in comp_prods
                union!(target_forms,
                    _forms_with_binding_step(steps, met))
            end
            for inh in comp_inhibitors
                union!(target_forms,
                    _forms_with_binding_step(steps, inh))
            end

            all_competing = union(
                comp_subs, comp_prods, comp_inhibitors)
            active = Symbol[]
            for f in sort(collect(target_forms))
                f in eligible_forms || continue
                haskey(bound, f) || continue
                isempty(intersect(
                    bound[f], all_competing)) || continue
                push!(active, f)
            end

            isempty(active) && continue
            active in seen && continue
            push!(seen, active)

            new_steps = copy(steps)
            de_form_map = Dict{Symbol, Symbol}()
            reg_g = next_g

            for cf in active
                de_name = _dead_end_form_name(
                    cf, bound[cf], dummy)
                de_form_map[cf] = de_name
                push!(new_steps, StepSpec(
                    [cf, dummy], [de_name], true, reg_g))
            end

            for s in steps
                from, to = step_forms(s)
                haskey(de_form_map, from) || continue
                haskey(de_form_map, to) || continue
                met = step_metabolite(s)
                from_de = de_form_map[from]
                to_de = de_form_map[to]
                g = s.kinetic_group
                if met !== nothing
                    push!(new_steps, StepSpec(
                        [from_de, met], [to_de],
                        s.is_equilibrium, g))
                else
                    push!(new_steps, StepSpec(
                        [from_de], [to_de],
                        s.is_equilibrium, g))
                end
            end

            push!(results, _dead_end_with_steps(
                spec, new_steps, reg_g))
            next_g = reg_g + 1
        end
    end
    results
end

_base_or_self(spec::MechanismSpec) = spec
_base_or_self(spec::AllostericMechanismSpec) = spec.base

"""
Build the new spec after adding a dead-end regulator's binding-step
kinetic group. For plain mechanisms: just +1. For allosteric: tag
the new group `:EqualRT` so it's one shared K (no extra T-state),
keeping the allosteric move at +1.
"""
_dead_end_with_steps(spec::MechanismSpec, new_steps, _new_g) =
    MechanismSpec(
        spec.reaction, new_steps,
        spec.n_fit_params_estimate + 1)

function _dead_end_with_steps(
    spec::AllostericMechanismSpec, new_steps, new_g::Int,
)
    new_tags = copy(spec.group_tags)
    new_tags[new_g] = :EqualRT
    AllostericMechanismSpec(
        MechanismSpec(spec.base.reaction, new_steps,
            spec.n_fit_params_estimate + 1),
        spec.catalytic_n,
        deepcopy(spec.allosteric_reg_sites),
        copy(spec.allosteric_multiplicities),
        new_tags, copy(spec.reg_ligand_tags),
        spec.n_fit_params_estimate + 1)
end

"""
    _spec_from_mechanism(m::Mechanism) → MechanismSpec

Bridge a `Mechanism` to a `MechanismSpec` for routing through spec-
based enumeration moves that have not yet been rewritten natively.
Steps are flattened (outer Vector = kinetic groups) and emitted as
StepSpec rows with the group index as `kinetic_group`. The PC field
is an estimate from the step shape — exact tracking happens at the
cache layer.
"""
function _spec_from_mechanism(m::Mechanism)
    legacy_rxn = _to_legacy_reaction(m.reaction)
    steps = StepSpec[]
    for (gi, group) in enumerate(m.steps)
        for s in group
            push!(steps,
                _stepspec_from_step(s, is_equilibrium(s), gi))
        end
    end
    pc = _n_fit_params_estimate_from_steps(steps)
    MechanismSpec(legacy_rxn, steps, pc)
end

"""
    _spec_from_mechanism(am::AllostericMechanism) → AllostericMechanismSpec

Bridge an `AllostericMechanism` to an `AllostericMechanismSpec`. The
base spec is built via the `Mechanism` overload; allo-state tags map
group index → tag; regulatory sites unpack into the dense
allosteric_reg_sites / allosteric_multiplicities / reg_ligand_tags
shape.
"""
function _spec_from_mechanism(am::AllostericMechanism)
    legacy_rxn = _to_legacy_reaction(am.reaction)
    steps = StepSpec[]
    for (gi, group) in enumerate(am.cat_steps)
        for s in group
            push!(steps,
                _stepspec_from_step(s, is_equilibrium(s), gi))
        end
    end
    pc = _n_fit_params_estimate_from_steps(steps)
    base = MechanismSpec(legacy_rxn, steps, pc)
    group_tags = Dict{Int, Symbol}(
        gi => am.cat_allo_states[gi]
        for gi in 1:length(am.cat_allo_states))
    reg_sites = Vector{Symbol}[]
    mults = Int[]
    reg_lig_tags = Dict{Symbol, Symbol}()
    for site in am.regulatory_sites
        site_ligs = [name(l) for l in ligands(site)]
        push!(reg_sites, site_ligs)
        push!(mults, multiplicity(site))
        for (i, lig) in enumerate(ligands(site))
            reg_lig_tags[name(lig)] = allo_states(site)[i]
        end
    end
    AllostericMechanismSpec(
        base, am.catalytic_multiplicity,
        reg_sites, mults, group_tags, reg_lig_tags, pc)
end

"""
Walk every step in `m` (or `am`) and parse each form's bound-metabolite
name set from the form Symbol (conformation only — `bound` is empty in
Mechanisms produced by the enumerator / legacy lift). Greedy match
against `met_set`. Used by the Mechanism-native dead-end-regulator move
to identify eligible base forms.
"""
function _bound_mets_from_form_name(form::Symbol, met_set::Set{Symbol})
    s = string(form)
    (s == "E" || s == "Estar") && return Set{Symbol}()
    body = startswith(s, "Estar_") ? s[7:end] :
           startswith(s, "E_") ? s[3:end] : s
    parts = split(body, "_")
    result = Set{Symbol}()
    i = 1
    while i <= length(parts)
        matched = false
        for len in length(parts):-1:1
            i + len - 1 > length(parts) && continue
            candidate = Symbol(join(parts[i:i+len-1], "_"))
            if candidate in met_set
                push!(result, candidate)
                i += len
                matched = true
                break
            end
        end
        matched || (i += 1)
    end
    result
end

"""
Reconstruct each form's bound-metabolite name set from the step graph
of `m`. `extra_mets` extends the known metabolite vocabulary (e.g., a
regulator declared in `rxn` but not yet bound). Returns
`Dict{form_name => Set{met_name}}`.
"""
function _bound_at_forms(m::Union{Mechanism, AllostericMechanism},
                         rxn::EnzymeReaction, extra_mets::Set{Symbol})
    met_set = Set{Symbol}()
    for ra in reactants(rxn); push!(met_set, name(metabolite(ra))); end
    for rm in regulators(rxn); push!(met_set, name(regulator(rm))); end
    union!(met_set, extra_mets)
    for group in steps(m), s in group
        bm = bound_metabolite(s)
        bm === nothing || push!(met_set, name(bm))
    end
    result = Dict{Symbol, Set{Symbol}}()
    for group in steps(m), s in group
        for sp in (from_species(s), to_species(s))
            fn = name(sp)
            haskey(result, fn) ||
                (result[fn] = _bound_mets_from_form_name(fn, met_set))
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
Build a conformation-only `Species` with the given form Symbol as the
conformation, no bound metabolites, and an empty residual. Matches the
representation produced by `_mechanism_from_legacy_sig` so dead-end /
mirror species built natively interoperate with enumerator output.
"""
_conformation_species(form::Symbol) =
    Species(Metabolite[], form, Residual())

"""
Extend `rxn`'s regulators with a new `CompetitiveInhibitor(name)`,
preserving every other field (reactants, allowed catalytic
multiplicities). `EnzymeReaction`'s inner constructor canonicalizes the
regulator order.
"""
function _add_competitive_inhibitor(rxn::EnzymeReaction, reg_name::Symbol)
    new_regs = copy(rxn.regulators)
    push!(new_regs, RegulatorMults(CompetitiveInhibitor(reg_name), Int[1]))
    EnzymeReaction(copy(rxn.reactants), new_regs,
                   copy(rxn.allowed_catalytic_multiplicities))
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
    _expand_add_dead_end_regulator_native(
        m, rxn, Set{Symbol}();
        exclude_regs=exclude_regs,
        wrap=(new_groups, new_g, new_reaction) ->
            Mechanism(new_reaction, new_groups))
end

function _expand_add_dead_end_regulator(
    am::AllostericMechanism, rxn::EnzymeReaction;
    exclude_regs::Set{Symbol}=Set{Symbol}(),
)
    # Allosteric ligands are excluded — dead-end is not the right move
    # for them (they belong on a regulatory site).
    allo_ligands = Set{Symbol}()
    for site in am.regulatory_sites, lig in ligands(site)
        push!(allo_ligands, name(lig))
    end
    _expand_add_dead_end_regulator_native(
        am, rxn, allo_ligands;
        exclude_regs=exclude_regs,
        wrap=(new_groups, new_g, new_reaction) -> AllostericMechanism(
            new_reaction, new_groups,
            vcat(am.cat_allo_states, [:EqualRT]),
            am.catalytic_multiplicity,
            copy(am.regulatory_sites)))
end

"""
Shared kernel for the Mechanism / AllostericMechanism dead-end
expansion. `additional_excluded` carries allosteric ligand names that
should not be eligible (empty for plain `Mechanism`). `wrap` constructs
the result from the assembled step groups, the new regulator's kinetic
group index, and the new reaction.
"""
function _expand_add_dead_end_regulator_native(
    m::Union{Mechanism, AllostericMechanism},
    rxn::EnzymeReaction,
    additional_excluded::Set{Symbol};
    exclude_regs::Set{Symbol},
    wrap,
)
    isempty(regulators(rxn)) && return typeof(m)[]

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
    isempty(eligible_regs) && return typeof(m)[]

    cat_forms = Set{Symbol}()
    for group in steps(m), s in group
        push!(cat_forms, name(from_species(s)))
        push!(cat_forms, name(to_species(s)))
    end

    n_groups_before = length(steps(m))
    n_steps_before = sum(length, steps(m); init = 0)
    results = typeof(m)[]

    for reg_name in eligible_regs
        # Make this regulator visible to form-name parsing so e.g.
        # parsing :E_I knows :I is a bound metabolite.
        bound = _bound_at_forms(m, rxn, Set([reg_name]))

        eligible_forms = Symbol[]
        for f in sort(collect(cat_forms))
            haskey(bound, f) || continue
            fb = bound[f]
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
                haskey(bound, f) || continue
                isempty(intersect(bound[f], all_competing)) || continue
                push!(active, f)
            end
            isempty(active) && continue
            active in seen && continue
            push!(seen, active)

            # Source-idx accounting: new steps continue past the
            # existing max so the result Mechanism's invariant
            # ("all source_idx non-zero" or "all zero") is preserved.
            next_src = n_steps_before + 1

            de_species_map = Dict{Symbol, Species}()
            reg_group_steps = Step[]
            for cf in active
                de_name = _dead_end_form_name(cf, bound[cf], reg_name)
                de_species = _conformation_species(de_name)
                de_species_map[cf] = de_species
                push!(reg_group_steps, Step(
                    _conformation_species(cf), de_species,
                    CompetitiveInhibitor(reg_name), true;
                    source_idx = next_src))
                next_src += 1
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
                             bound_metabolite(s), is_equilibrium(s);
                             source_idx = next_src))
                    next_src += 1
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
            push!(results, wrap(new_groups, n_groups_before + 1,
                                new_reaction))
        end
    end
    results
end

"""
    _expand_to_allosteric(spec, reaction)
        → Vector{AllostericMechanismSpec}

Convert a non-allosteric `MechanismSpec` to allosteric. Emits the
all-`:EqualRT` baseline plus one variant per kinetic group with that
group set to `:OnlyR`. Total: `n_groups + 1` specs.

Cost: +1 (for `L`). Other tag deltas are zero relative to the
all-`:EqualRT` baseline.
"""
function _expand_to_allosteric(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReactionLegacy),
)
    cn = oligomeric_state(reaction)
    base_pc = spec.n_fit_params_estimate

    group_info = _group_info(spec.steps)
    groups_sorted = sort!(collect(keys(group_info)))

    # All-:EqualRT baseline: cheapest non-default tag everywhere.
    base_tags = Dict{Int, Symbol}(g => :EqualRT for g in groups_sorted)

    # Each emitted variant gets its OWN base MechanismSpec wrapping a
    # copied steps vector. Sharing `spec` across variants would let
    # later in-place mutation (e.g. _canonicalize! during dedup!) on one
    # variant's base corrupt the shared steps and desynchronize the
    # other variants' group_tags from their now-renumbered groups.
    fresh_base() = MechanismSpec(spec.reaction, copy(spec.steps), base_pc)

    results = AllostericMechanismSpec[]
    push!(results, AllostericMechanismSpec(
        fresh_base(), cn,
        Vector{Symbol}[], Int[],
        copy(base_tags), Dict{Symbol, Symbol}(),
        base_pc + 1))
    for g in groups_sorted
        new_tags = copy(base_tags)
        new_tags[g] = :OnlyR
        push!(results, AllostericMechanismSpec(
            fresh_base(), cn,
            Vector{Symbol}[], Int[],
            new_tags, Dict{Symbol, Symbol}(),
            base_pc + 1))
    end
    results
end

_expand_to_allosteric(
    ::AllostericMechanismSpec, @nospecialize(::EnzymeReactionLegacy),
) = AllostericMechanismSpec[]

"""
    _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

Mechanism-native overload: convert a non-allosteric `Mechanism` into
allosteric variants. Emits the all-`:EqualRT` baseline plus one
variant per kinetic group with that group set to `:OnlyR`
(`n_groups + 1` variants total). The new mechanism inherits `rxn`'s
oligomeric state as `catalytic_multiplicity`; regulatory_sites is
empty (allosteric regulators are added later via
`_expand_add_allosteric_regulator`). Steps are reused by reference.
"""
function _expand_to_allosteric(m::Mechanism, rxn::EnzymeReaction)
    cn = only(allowed_catalytic_multiplicities(rxn))
    n_g = length(m.steps)
    base_tags = Symbol[:EqualRT for _ in 1:n_g]
    empty_sites = RegulatorySite[]
    results = AllostericMechanism[]
    push!(results, AllostericMechanism(
        m.reaction, copy(m.steps), copy(base_tags),
        cn, copy(empty_sites)))
    for g in 1:n_g
        new_tags = copy(base_tags)
        new_tags[g] = :OnlyR
        push!(results, AllostericMechanism(
            m.reaction, copy(m.steps), new_tags,
            cn, copy(empty_sites)))
    end
    results
end

"""
    _expand_to_allosteric(am::AllostericMechanism, rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

`AllostericMechanism` input is already allosteric — no-op move.
Matches the spec overload's empty-vector return.
"""
_expand_to_allosteric(::AllostericMechanism, ::EnzymeReaction) =
    AllostericMechanism[]

"""
Map kinetic_group → (is_re::Bool, iso_only::Bool). `iso_only` flags
groups whose members are all isomerization steps (no metabolite).
"""
function _group_info(steps::Vector{StepSpec})
    info = Dict{Int, Tuple{Bool, Bool}}()
    members = Dict{Int, Vector{Int}}()
    for (i, s) in enumerate(steps)
        push!(get!(members, s.kinetic_group, Int[]), i)
    end
    for (g, idxs) in members
        is_re = steps[idxs[1]].is_equilibrium
        iso_only = all(
            step_metabolite(steps[i]) === nothing for i in idxs)
        info[g] = (is_re, iso_only)
    end
    info
end

"""
ΔP when a kinetic group's allo_state changes from `from` to `to`. Group cost:
1 param for `:EqualRT`/`:OnlyR`/`:OnlyT`, 2 params for `:NonequalRT`.
RE groups have 1 K-style param per "unit"; SS groups have 2 (kf + kr).
"""
function _allo_state_delta(from::Symbol, to::Symbol, is_re::Bool)
    cost = t -> (t == :NonequalRT ? 2 : 1)
    factor = is_re ? 1 : 2
    factor * (cost(to) - cost(from))
end

"""
ΔP when a regulatory ligand's allo_state changes from `from` to `to`. Each
ligand contributes 1 binding K per state — `:NonequalRT` has K_R + K_T
(2 params), `:EqualRT`/`:OnlyR`/`:OnlyT` have 1.
"""
function _allo_lig_state_delta(from::Symbol, to::Symbol)
    cost = t -> (t == :NonequalRT ? 2 : 1)
    cost(to) - cost(from)
end

"""
    _expand_add_allosteric_regulator(spec, reaction)
        → Vector{AllostericMechanismSpec}

Add one allosteric regulator (not already present) to the mechanism,
either at a new regulatory site or at an existing one. For each
ligand × site option × tag flavor, emit a variant. Each variant adds
+1 (one K_R per reg site) plus a per-tag delta on top.
"""
function _expand_add_allosteric_regulator(
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReactionLegacy),
)
    existing_allo = Set{Symbol}()
    for site in spec.allosteric_reg_sites
        for lig in site
            push!(existing_allo, lig)
        end
    end

    existing_de = Set{Symbol}()
    for s in spec.base.steps
        for sym in Iterators.flatten((s.reactants, s.products))
            m = match(r"^(.+)__reg\d*$", string(sym))
            m !== nothing &&
                push!(existing_de, Symbol(m.captures[1]))
        end
    end

    roles = regulator_roles(reaction)
    new_regs = Symbol[]
    for (name, role) in roles
        (role == :unknown || role == :allosteric) || continue
        name in existing_allo && continue
        name in existing_de && continue
        push!(new_regs, name)
    end
    sort!(new_regs)
    isempty(new_regs) && return AllostericMechanismSpec[]

    results = AllostericMechanismSpec[]
    for reg in new_regs
        n_sites = length(spec.allosteric_reg_sites)
        # Enumerate non-:EqualRT states for any site (new or existing)
        for tag in (:OnlyR, :OnlyT, :NonequalRT)
            for site_idx in 0:n_sites
                new_sites = deepcopy(spec.allosteric_reg_sites)
                new_mults = copy(spec.allosteric_multiplicities)
                new_lig_tags = copy(spec.reg_ligand_tags)

                if site_idx == 0
                    push!(new_sites, Symbol[reg])
                    push!(new_mults, spec.catalytic_n)
                else
                    push!(new_sites[site_idx], reg)
                end
                new_lig_tags[reg] = tag

                # Cost: one K binding param at this site (+1) plus
                # a per-tag delta vs the default `:EqualRT` "free"
                # cost which already adds 1.
                delta_cost = _allo_lig_state_delta(:EqualRT, tag) + 1

                push!(results, AllostericMechanismSpec(
                    MechanismSpec(spec.base.reaction, copy(spec.base.steps),
                                  spec.base.n_fit_params_estimate),
                    spec.catalytic_n,
                    new_sites, new_mults,
                    copy(spec.group_tags), new_lig_tags,
                    spec.n_fit_params_estimate + delta_cost))
            end
        end
        # Enumerate :EqualRT only for existing sites where at least one
        # ligand is already non-:EqualRT (single-ligand or all-:EqualRT
        # site cancels identically — constructor would reject).
        for site_idx in 1:n_sites
            existing_ligs = spec.allosteric_reg_sites[site_idx]
            any(spec.reg_ligand_tags[l] != :EqualRT
                for l in existing_ligs) || continue
            new_sites = deepcopy(spec.allosteric_reg_sites)
            new_mults = copy(spec.allosteric_multiplicities)
            new_lig_tags = copy(spec.reg_ligand_tags)
            push!(new_sites[site_idx], reg)
            new_lig_tags[reg] = :EqualRT
            delta_cost = _allo_lig_state_delta(:EqualRT, :EqualRT) + 1
            push!(results, AllostericMechanismSpec(
                spec.base, spec.catalytic_n,
                new_sites, new_mults,
                copy(spec.group_tags), new_lig_tags,
                spec.n_fit_params_estimate + delta_cost))
        end
    end
    results
end

_expand_add_allosteric_regulator(
    ::MechanismSpec, @nospecialize(::EnzymeReactionLegacy),
) = AllostericMechanismSpec[]

"""
    _expand_add_allosteric_regulator(am::AllostericMechanism,
                                     rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

Add one `AllostericRegulator` declared in `rxn` but not yet bound by
`am` (neither at a regulatory site nor as a catalytic-step dead-end
inhibitor). For each (new ligand, target site, tag) combination, emit
a variant:

  * target site ∈ {new site} ∪ {existing sites}
  * tag ∈ {:OnlyR, :OnlyT, :NonequalRT} for any target site
  * tag = :EqualRT only at an existing site that already has at least
    one non-`:EqualRT` ligand (otherwise the `RegulatorySite`
    constructor's all-`:EqualRT` rule would reject the variant).

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
    for site in am.regulatory_sites, lig in ligands(site)
        push!(existing_allo, name(lig))
    end

    existing_de = Set{Symbol}()
    for group in am.cat_steps, s in group
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
        n_sites = length(am.regulatory_sites)
        # Non-:EqualRT tags at any (new or existing) site.
        for tag in (:OnlyR, :OnlyT, :NonequalRT)
            for site_idx in 0:n_sites
                push!(results,
                    _make_am_with_added_reg(am, reg, tag, site_idx))
            end
        end
        # :EqualRT at an existing site only when that site already has
        # at least one non-:EqualRT ligand (avoids the constructor's
        # all-:EqualRT single-ligand rejection / identical-cancellation).
        for site_idx in 1:n_sites
            site = am.regulatory_sites[site_idx]
            any(st != :EqualRT for st in allo_states(site)) || continue
            push!(results,
                _make_am_with_added_reg(am, reg, :EqualRT, site_idx))
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
        for site in am.regulatory_sites
            push!(new_sites, site)
        end
        push!(new_sites, RegulatorySite(
            AllostericRegulator[AllostericRegulator(reg)],
            am.catalytic_multiplicity,
            Symbol[tag]))
    else
        for (i, site) in enumerate(am.regulatory_sites)
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
    AllostericMechanism(am.reaction, copy(am.cat_steps),
                        copy(am.cat_allo_states),
                        am.catalytic_multiplicity, new_sites)
end

"""
    _expand_add_allosteric_regulator(m::Mechanism, rxn::EnzymeReaction)
        → Vector{AllostericMechanism}

Non-allosteric input: no-op (matches the spec overload's empty return).
The dispatch shape here ensures callers walking a mixed Mechanism /
AllostericMechanism collection don't need to type-check.
"""
_expand_add_allosteric_regulator(::Mechanism, ::EnzymeReaction) =
    AllostericMechanism[]

"""
    _expand_change_allo_state(spec, reaction)
        → Vector{AllostericMechanismSpec}

Relax one allo_state from a "constrained" tag (`:EqualRT`, `:OnlyR`,
`:OnlyT`) to `:NonequalRT`. Tags already at `:NonequalRT` are skipped
(no relaxation possible). Each variant adds the corresponding param
delta.

For an iso-only group already at `:OnlyR`, the `:OnlyT` direction
is forbidden by the constructor — but the move only goes to
`:NonequalRT` so this isn't a concern here.
"""
function _expand_change_allo_state(
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReactionLegacy),
)
    base_steps = spec.base.steps
    group_info = _group_info(base_steps)
    results = AllostericMechanismSpec[]

    for (g, tag) in spec.group_tags
        tag == :NonequalRT && continue
        haskey(group_info, g) || continue
        is_re, _ = group_info[g]
        delta = _allo_state_delta(tag, :NonequalRT, is_re)
        new_tags = copy(spec.group_tags)
        new_tags[g] = :NonequalRT
        push!(results, AllostericMechanismSpec(
            MechanismSpec(spec.base.reaction, copy(spec.base.steps),
                          spec.base.n_fit_params_estimate),
            spec.catalytic_n,
            deepcopy(spec.allosteric_reg_sites),
            copy(spec.allosteric_multiplicities),
            new_tags, copy(spec.reg_ligand_tags),
            spec.n_fit_params_estimate + delta))
    end

    for (lig, tag) in spec.reg_ligand_tags
        tag == :NonequalRT && continue
        delta = _allo_lig_state_delta(tag, :NonequalRT)
        new_lig_tags = copy(spec.reg_ligand_tags)
        new_lig_tags[lig] = :NonequalRT
        push!(results, AllostericMechanismSpec(
            MechanismSpec(spec.base.reaction, copy(spec.base.steps),
                          spec.base.n_fit_params_estimate),
            spec.catalytic_n,
            deepcopy(spec.allosteric_reg_sites),
            copy(spec.allosteric_multiplicities),
            copy(spec.group_tags), new_lig_tags,
            spec.n_fit_params_estimate + delta))
    end

    results
end

_expand_change_allo_state(
    ::MechanismSpec, @nospecialize(::EnzymeReactionLegacy),
) = AllostericMechanismSpec[]

"""
    _expand_change_allo_state(am::AllostericMechanism)
        → Vector{AllostericMechanism}

Mechanism-native overload. Relax one allo state from a "constrained"
tag (`:EqualRT`, `:OnlyR`, `:OnlyT`) to `:NonequalRT`. Variants are
emitted for each catalytic kinetic group tag and each regulatory
ligand tag that is not already `:NonequalRT`. The base catalytic
steps, multiplicity, and untouched tags are preserved.
"""
function _expand_change_allo_state(am::AllostericMechanism)
    results = AllostericMechanism[]

    # Catalytic-group tag relaxations: cat_allo_states[g] :…→ :NonequalRT.
    for g in 1:length(am.cat_allo_states)
        am.cat_allo_states[g] == :NonequalRT && continue
        new_states = copy(am.cat_allo_states)
        new_states[g] = :NonequalRT
        push!(results, AllostericMechanism(
            am.reaction, copy(am.cat_steps), new_states,
            am.catalytic_multiplicity,
            copy(am.regulatory_sites)))
    end

    # Regulatory-ligand tag relaxations: walk each (site, ligand) pair.
    for (si, site) in enumerate(am.regulatory_sites)
        for (li, _) in enumerate(ligands(site))
            allo_states(site)[li] == :NonequalRT && continue
            new_sites = copy(am.regulatory_sites)
            new_states = copy(allo_states(site))
            new_states[li] = :NonequalRT
            new_sites[si] = RegulatorySite(
                copy(ligands(site)),
                multiplicity(site),
                new_states)
            push!(results, AllostericMechanism(
                am.reaction, copy(am.cat_steps),
                copy(am.cat_allo_states),
                am.catalytic_multiplicity,
                new_sites))
        end
    end

    results
end

"""
    _expand_change_allo_state(m::Mechanism) → Vector{AllostericMechanism}

Non-allosteric input: no-op (matches the spec overload's empty return).
"""
_expand_change_allo_state(::Mechanism) =
    AllostericMechanism[]

"""
    AllostericEnzymeMechanism(spec::AllostericMechanismSpec) →
        AllostericEnzymeMechanism

Build the singleton allosteric type from a spec. `group_tags`
sorted by group id for canonical type identity. Each regulator
site contributes its ligand list, multiplicity, and the
non-default per-ligand tag entries.
"""
function AllostericEnzymeMechanism(spec::AllostericMechanismSpec)
    # Allosteric ligands sit at reg sites in the assembled allosteric
    # mechanism; they don't appear in the catalytic-cycle steps and
    # must not show up in the catalytic `EnzymeMechanism`'s
    # `regulators` list.
    allo_set = Set{Symbol}()
    for site in spec.allosteric_reg_sites
        for l in site; push!(allo_set, l); end
    end
    cm = EnzymeMechanism(spec.base; exclude_regs=allo_set)

    n_groups = length(unique(s.kinetic_group for s in spec.base.steps))
    cat_states = ntuple(g -> spec.group_tags[g], n_groups)

    reg_sites = ntuple(length(spec.allosteric_reg_sites)) do i
        ligs = Tuple(spec.allosteric_reg_sites[i])
        mult = spec.allosteric_multiplicities[i]
        lig_states = ntuple(k -> spec.reg_ligand_tags[ligs[k]],
                            length(ligs))
        (ligs, mult, lig_states)
    end

    AllostericEnzymeMechanism(cm, (spec.catalytic_n, cat_states), reg_sites)
end

function _push_to_dict!(
    result::Dict{Int, Vector{AbstractMechanismSpec}},
    spec::AbstractMechanismSpec)
    push!(get!(result, spec.n_fit_params_estimate,
        AbstractMechanismSpec[]), spec)
end

"""
    expand_mechanisms(specs, reaction) → Dict{Int, Vector{AbstractMechanismSpec}}

Apply all +1 and +2 expansion moves. Results grouped
by target n_fit_params_estimate.
"""
function expand_mechanisms(
    specs::Vector{<:AbstractMechanismSpec},
    @nospecialize(reaction::EnzymeReactionLegacy))
    result = Dict{Int, Vector{AbstractMechanismSpec}}()
    for spec in specs
        _add_expansions!(result, spec, reaction)
    end
    result
end

"""
    expand_mechanisms(mechs::Vector, reaction)
        → Dict{Int, Vector{Union{Mechanism, AllostericMechanism}}}

Mechanism-native overload. Applies all expansion moves to each input
mechanism (RE→SS, split kinetic group, add dead-end regulator,
to-allosteric, add allosteric regulator, change allo state) and
buckets results by their `_n_fit_params_estimate` upper bound.

Accepts a heterogeneous mix of `Mechanism` and `AllostericMechanism`
inputs — the spec-overload pipeline returns the same mix because
`_expand_to_allosteric` promotes a `Mechanism` to an
`AllostericMechanism`. The reaction argument may be either a concrete
`EnzymeReaction` or the legacy parametric form (the per-move helpers
dispatch on `EnzymeReaction`; the boundary adapter at the bottom of
this file forwards from legacy).
"""
function expand_mechanisms(
    mechs::Vector{<:Union{Mechanism, AllostericMechanism}},
    rxn::EnzymeReaction)
    result = Dict{Int, Vector{Union{Mechanism, AllostericMechanism}}}()
    for m in mechs
        _add_expansions_mech!(result, m, rxn)
    end
    result
end

# Adapter so callers holding the legacy parametric reaction form
# (e.g. `IdentifyRateEquationProblem`) can drive the Mechanism-side
# pipeline. Pulls the concrete `EnzymeReaction` from any mechanism in
# the input vector (every `Mechanism`/`AllostericMechanism` stores its
# own concrete reaction). Empty input is a no-op so no concrete
# reaction is needed for the empty result.
function expand_mechanisms(
    mechs::Vector{<:Union{Mechanism, AllostericMechanism}},
    @nospecialize(::EnzymeReactionLegacy))
    isempty(mechs) &&
        return Dict{Int,
                    Vector{Union{Mechanism, AllostericMechanism}}}()
    expand_mechanisms(mechs, first(mechs).reaction)
end

function _add_expansions_mech!(
    result::Dict{Int, Vector{Union{Mechanism, AllostericMechanism}}},
    m::Union{Mechanism, AllostericMechanism},
    rxn::EnzymeReaction)
    for s in _expand_re_to_ss(m)
        _push_mech!(result, s)
    end
    for s in _expand_split_kinetic_group(m)
        _push_mech!(result, s)
    end
    for s in _expand_add_dead_end_regulator(m, rxn)
        _push_mech!(result, s)
    end
    for s in _expand_to_allosteric(m, rxn)
        _push_mech!(result, s)
    end
    for s in _expand_add_allosteric_regulator(m, rxn)
        _push_mech!(result, s)
    end
    for s in _expand_change_allo_state(m)
        _push_mech!(result, s)
    end
end

function _push_mech!(
    result::Dict{Int, Vector{Union{Mechanism, AllostericMechanism}}},
    m::Union{Mechanism, AllostericMechanism})
    pc = _n_fit_params_estimate(m)
    push!(get!(result, pc,
               Union{Mechanism, AllostericMechanism}[]), m)
end

function _add_expansions!(
    result::Dict{Int, Vector{AbstractMechanismSpec}},
    spec::AbstractMechanismSpec,
    @nospecialize(reaction::EnzymeReactionLegacy))
    for s in _expand_re_to_ss(spec)
        _push_to_dict!(result, s)
    end
    for s in _expand_split_kinetic_group(spec)
        _push_to_dict!(result, s)
    end
    for s in _expand_add_dead_end_regulator(
            spec, reaction)
        _push_to_dict!(result, s)
    end
    for s in _expand_to_allosteric(spec, reaction)
        _push_to_dict!(result, s)
    end
    for s in _expand_add_allosteric_regulator(
            spec, reaction)
        _push_to_dict!(result, s)
    end
    for s in _expand_change_allo_state(spec, reaction)
        _push_to_dict!(result, s)
    end
end

# --- Dedup ---

"""Sort key: shape first, kinetic_group last."""
function _step_sort_key(s::StepSpec)
    (sort(s.reactants), sort(s.products),
     s.is_equilibrium, s.kinetic_group)
end

"""
Renumber kinetic groups in canonical order: assign 1 to the first
group encountered (by sorted step order), 2 to the next new group,
etc. Returns `Dict{old_group → new_group}`.
"""
function _canonical_group_renumber(steps::Vector{StepSpec})
    rename = Dict{Int, Int}()
    next_g = 0
    for s in steps
        haskey(rename, s.kinetic_group) && continue
        next_g += 1
        rename[s.kinetic_group] = next_g
    end
    rename
end

function _canonicalize!(spec::MechanismSpec)
    sort!(spec.steps, by=_step_sort_key)
    rename = _canonical_group_renumber(spec.steps)
    for (i, s) in enumerate(spec.steps)
        spec.steps[i] = StepSpec(
            s.reactants, s.products,
            s.is_equilibrium, rename[s.kinetic_group])
    end
    rename
end

function _canonicalize!(spec::AllostericMechanismSpec)
    rename = _canonicalize!(spec.base)
    new_group_tags = Dict{Int, Symbol}(
        rename[g] => t for (g, t) in spec.group_tags
    )
    empty!(spec.group_tags)
    for (g, t) in new_group_tags
        spec.group_tags[g] = t
    end
    for site in spec.allosteric_reg_sites
        sort!(site)
    end
    # Sort sites themselves (with multiplicities) by content
    if length(spec.allosteric_reg_sites) >= 2
        perm = sortperm(spec.allosteric_reg_sites)
        spec.allosteric_reg_sites .=
            spec.allosteric_reg_sites[perm]
        spec.allosteric_multiplicities .=
            spec.allosteric_multiplicities[perm]
    end
    rename
end

function _dedup_key(spec::MechanismSpec)
    Tuple(
        (Tuple(sort(s.reactants)),
         Tuple(sort(s.products)),
         s.is_equilibrium,
         s.kinetic_group)
        for s in spec.steps)
end

function _dedup_key(spec::AllostericMechanismSpec)
    base_key = _dedup_key(spec.base)
    sorted_tags = Tuple(
        sort(collect(spec.group_tags); by=first))
    sorted_lig_tags = Tuple(
        sort(collect(spec.reg_ligand_tags); by=first))
    (base_key, spec.catalytic_n,
     Tuple(Tuple.(spec.allosteric_reg_sites)),
     Tuple(spec.allosteric_multiplicities),
     sorted_tags, sorted_lig_tags)
end

"""
    dedup!(cache) → cache

Remove structural duplicates from each n_fit_params_estimate bucket
via canonical form comparison.
"""
function dedup!(
    cache::Dict{Int, Vector{AbstractMechanismSpec}})
    for (pc, specs) in cache
        for s in specs
            _canonicalize!(s)
        end
        unique!(s -> _dedup_key(s), specs)
        isempty(specs) && delete!(cache, pc)
    end
    cache
end

# ─── Mechanism-based dedup ─────────────────────────────────────────────
#
# Canonical key for a `Step` that ignores `source_idx` (presentation
# metadata). `Step`'s own `==`/`hash` already ignore `source_idx`; the
# key tuple's sole job is to give `sort!` a deterministic ordering so
# two physically-equivalent `Mechanism`s end up with identical step
# storage and therefore identical struct-based `hash` / `==`.
_step_canonical_key(s::Step) =
    (hash(s.from_species), hash(s.to_species),
     hash(s.bound_metabolite), s.is_equilibrium)

"""
Sort steps within each kinetic group by `_step_canonical_key`, then
sort the outer group vector by the canonical key of its first step.
Mutates the inner vectors of `m.steps` (and the outer vector itself)
in place. `source_idx` values stay attached to their `Step`s; only
storage order changes.
"""
function _canonicalize_mechanism!(m::Mechanism)
    for group in m.steps
        sort!(group; by = _step_canonical_key)
    end
    sort!(m.steps; by = group -> _step_canonical_key(first(group)))
    m
end

function _canonicalize_mechanism!(am::AllostericMechanism)
    # Catalytic-side step storage is canonicalized identically to
    # `Mechanism`; the regulatory side is also canonicalized so two
    # `AllostericMechanism`s differing only in site presentation order
    # collapse. Ligand order within a site is fixed by the
    # `RegulatorySite` constructor, so only the outer site vector and
    # the parallel `cat_allo_states` vector need reordering. The inner
    # sort must run BEFORE computing the outer permutation so the
    # per-group "first step" key reflects the canonical inner order.
    for group in am.cat_steps
        sort!(group; by = _step_canonical_key)
    end
    perm = sortperm(1:length(am.cat_steps);
                    by = g -> _step_canonical_key(first(am.cat_steps[g])))
    permute!(am.cat_steps, perm)
    permute!(am.cat_allo_states, perm)
    sort!(am.regulatory_sites; by = _regulatory_site_canonical_key)
    am
end

_regulatory_site_canonical_key(site::RegulatorySite) =
    (Tuple(hash(l) for l in site.ligands),
     site.multiplicity,
     Tuple(site.allo_states))

function dedup!(cache::Dict{Int, Vector{Mechanism}})
    for (pc, mechs) in cache
        for m in mechs
            _canonicalize_mechanism!(m)
        end
        unique!(mechs)
        isempty(mechs) && delete!(cache, pc)
    end
    cache
end

function dedup!(cache::Dict{Int, Vector{AllostericMechanism}})
    for (pc, mechs) in cache
        for m in mechs
            _canonicalize_mechanism!(m)
        end
        unique!(mechs)
        isempty(mechs) && delete!(cache, pc)
    end
    cache
end

"""
    dedup!(cache::Dict{Int, Vector{Union{Mechanism, AllostericMechanism}}})

Heterogeneous-bucket dedup for the `identify_rate_equation` pipeline,
which stores allosteric promotions alongside non-allosteric mechanisms
in the same param-count bucket. Canonicalizes each mechanism in place
via the type-specific `_canonicalize_mechanism!` overload, then runs
`unique!` so structurally-equivalent mechanisms collapse.
"""
function dedup!(cache::Dict{Int,
                            Vector{Union{Mechanism, AllostericMechanism}}})
    for (pc, mechs) in cache
        for m in mechs
            _canonicalize_mechanism!(m)
        end
        unique!(mechs)
        isempty(mechs) && delete!(cache, pc)
    end
    cache
end

# ─── Rate-equation canonical hash ──────────────────────────────────────

"""
Build the canonical text + name_map. Internal helper exposed
to `_canonical_rate_eq_hash_data` so callers that need to project
cached parameters across hash-equivalent specs can also retrieve
the name_map.

Strategy: walk `parameters(m, Full)` to discover every parameter
symbol the mechanism could mention (including dependents that
appear in the v-line and on constraint LHSes), scan the
rate-equation body to find each parameter's first-appearance
position, then rename them as `p_1, p_2, …` in first-appearance
order. `:E_total` is in `Full` but is excluded from renaming;
`:Keq` and metabolite names are not in `Full` and aren't renamed.

Multi-symbol Wegscheider and Haldane closure lines are kept in
the body — their LHS appears in v via a runtime substitution, so
they encode real parameterization (two mechanisms with the same
v-line but different dependent-parameter choice must hash
differently). Single-symbol `(substituted into v)` ties are
stripped because their LHS never appears in v.

`parameters(m, Full)` is defined for both `EnzymeMechanism` and
`AllostericEnzymeMechanism`. Allosteric coverage includes T-state
names, regulator-site names, and the allosteric coupling `L`
automatically.
"""
function _canonicalize_rate_eq_with_map(m::AbstractEnzymeMechanism)
    raw_body = rate_equation_string(m)
    raw_body === nothing && error(
        "rate_equation_string returned nothing for $(typeof(m))")
    body = String(raw_body)

    # Drop display-only lines so the eq_hash is invariant to the
    # section a constraint was emitted under and to the per-spec
    # choice of which absorbed symbol points at which representative:
    # - Destructure header lines (no semantic content for hashing).
    # - Section header lines like `# Wegscheider constraints:`.
    # - Single-symbol equality lines carrying the
    #   `(substituted into v)` annotation. Their LHS is folded into
    #   the v polynomial via the kinetic-group rename map, so it
    #   never appears in v itself; the line encodes only display
    #   relabeling. Two hash-equivalent mechanisms can produce the
    #   same v under different relabeling structures (e.g.,
    #   {K10=K8, K8=K4} vs {K10=K4, K8=K4}); dropping these lines
    #   makes both hash identically. Multi-symbol Wegscheider
    #   and Haldane closures stay — their LHS appears in v via a
    #   runtime substitution, so they're real parameterization info.
    # The single-symbol-equality regex restricts to true RE/SS rate-
    # constant symbols (K\d+, K\d+_T, k\d+f, k\d+r, k\d+f_T,
    # k\d+r_T) so regulator-K-to-regulator-K lines (which don't flow
    # through the polynomial-level absorption pipeline) aren't
    # dropped.
    raw_lines = split(body, '\n')

    is_destructure(ln) = occursin(r"^\s*\(; .* = (params|concs)$", ln)
    is_section_header(ln) = occursin(r"^# .+ constraints:$", ln)

    sym_pattern = "(?:K\\d+(?:_T)?|k\\d+[fr](?:_T)?)"
    annotation_escaped = replace(
        ANNOTATION_SUBSTITUTED, "(" => "\\(", ")" => "\\)")
    single_eq_re = Regex(
        "^\\s*$sym_pattern\\s*=\\s*$sym_pattern$annotation_escaped\$")
    is_single_eq(ln) = occursin(single_eq_re, ln)

    body = join(
        [ln for ln in raw_lines
         if !is_destructure(ln) && !is_section_header(ln) &&
            !is_single_eq(ln)],
        '\n')

    skip = (:E_total,)
    pnames = String[String(p) for p in parameters(m, Full)
                    if p ∉ skip]

    first_pos = Dict{String,Int}()
    for name in pnames
        rx = Regex("\\b" * name * "\\b")
        m_pos = match(rx, body)
        m_pos === nothing && continue
        first_pos[name] = m_pos.offset
    end
    appearing = collect(keys(first_pos))

    ordered = sort(appearing; by=name -> (first_pos[name], name))
    name_map = Dict(name => "p_$i"
                    for (i, name) in enumerate(ordered))

    # Substitute longest first to prevent prefix collisions
    # (e.g., rename `K1_T` before `K1`).
    for name in sort(appearing; by=length, rev=true)
        body = replace(body,
            Regex("\\b" * name * "\\b") => name_map[name])
    end

    # Re-sort multiplicative factors within each monomial run so the
    # text order matches the p_i ordering, not the original-symbol
    # ordering used by _poly_to_expr. Two mechanisms whose canonical
    # mapping pairs differ only by which raw step index (e.g., K9 vs
    # K10) plays a given role will then render with their factors at
    # matching positions. A "run" is a sequence of `<atom> * <atom>
    # * ...` terms where each atom is a bare identifier or a
    # `name ^ digits` power. The match stops at parens, +/-, or `/`.
    factor_atom = "[A-Za-z_]\\w*(?:\\s*\\^\\s*\\d+)?"
    run_re = Regex("$factor_atom(?:\\s*\\*\\s*$factor_atom)+")
    body = replace(body, run_re => _sort_run_factors)

    canonical = strip(replace(body, r"\s+" => " "))
    (canonical, name_map)
end

"""
Sort the `<atom> * <atom> * ...` factors of one multiplicative run.
The full match is split on `\\s*\\*\\s*`, each factor is keyed by
`(kind, p_index, exponent, lex_string)`, and the sorted factors are
joined back with ` * `.
"""
function _sort_run_factors(run::AbstractString)
    factors = split(strip(run), r"\s*\*\s*")
    sort!(factors; by=_factor_sort_key)
    join(factors, " * ")
end

"""Sort key for one factor of a multiplicative run."""
function _factor_sort_key(f::AbstractString)
    m = match(r"^p_(\d+)(?:\s*\^\s*(\d+))?$", f)
    if m !== nothing
        # captures[1] is `(\d+)` and the overall regex requires it to match,
        # so it can never be Nothing here. Assert the type to help JET's
        # type inference; without this JET flags `parse(Int, ::Nothing)` as
        # a potential error on the Union{Nothing, SubString} type.
        base = m.captures[1]::SubString{String}
        exp = m.captures[2]
        return (0, parse(Int, base),
                exp === nothing ? 1 : parse(Int, exp::SubString{String}),
                "")
    end
    (1, 0, 0, String(f))
end

"""
Return `(UInt64 hash, 16-char hex display string, name_map)`.
The single entry point for canonical hashing; `_canonical_rate_eq_hash`
delegates here so the canonicalizer runs once and callers that need the
name_map can retrieve it without a second pass.

Hash collision probability over 10⁴ mechanisms is ~10⁻¹² with
Julia's built-in `hash(::String)::UInt64`.
"""
function _canonical_rate_eq_hash_data(m::AbstractEnzymeMechanism)
    canonical, name_map = _canonicalize_rate_eq_with_map(m)
    h = hash(canonical)
    (h, string(h, base=16, pad=16), name_map)
end

"""
Hash a mechanism's canonicalized rate equation. Returns the
`UInt64` hash.
"""
function _canonical_rate_eq_hash(m::AbstractEnzymeMechanism)
    first(_canonical_rate_eq_hash_data(m))
end

# Adapters that accept the concrete EnzymeReaction by converting to
# the parametric EnzymeReactionLegacy form. The enumeration primitives
# below dispatch on EnzymeReactionLegacy.
#
# The public `init_mechanisms(::EnzymeReaction)` returns a
# `Vector{Mechanism}` — concrete structs from `types.jl`. The legacy
# pipeline still emits `MechanismSpec` values internally;
# `_init_mechanism_specs` exposes that internal spec list to the few
# downstream consumers (`expand_mechanisms`, `dedup!`, tests) that
# have not yet been switched to `Mechanism`.
init_mechanisms(r::EnzymeReaction) =
    [_mechanism_from_spec(spec) for spec in _init_mechanism_specs(r)]

# Mechanism-returning entry point for the EnzymeReactionLegacy form
# stored in `IdentifyRateEquationProblem.reaction`. Bridges through
# `_init_mechanism_specs` so the EnzymeReaction overload above stays
# the canonical user-facing entry.
init_mechanisms(r::EnzymeReactionLegacy, ::Type{Mechanism}) =
    [_mechanism_from_spec(spec) for spec in _init_mechanism_specs(r)]

"""
    _init_mechanism_specs(reaction) -> Vector{MechanismSpec}

Internal entry point that returns the enumerator's native
`MechanismSpec` output without the boundary conversion to
`Mechanism`. Used by `expand_mechanisms`/`dedup!` callers (and tests
that introspect spec internals) until Tasks 5.3/5.4 switch them.
Accepts both `EnzymeReaction` and `EnzymeReactionLegacy` so that
`IdentifyRateEquationProblem.reaction::EnzymeReactionLegacy` callers
work without an extra conversion.
"""
_init_mechanism_specs(r::EnzymeReaction) =
    init_mechanisms(_to_legacy_reaction(r))

_init_mechanism_specs(@nospecialize(r::EnzymeReactionLegacy)) =
    init_mechanisms(r)

"""
    _mechanism_from_spec(spec::MechanismSpec) → Mechanism

Convert an enumerator-produced `MechanismSpec` to a `Mechanism` by
routing through `EnzymeMechanism(spec)` and lifting back. Steps are
grouped by `kinetic_group` and renumbered in first-occurrence order.
"""
_mechanism_from_spec(spec::MechanismSpec) =
    Mechanism(EnzymeMechanism(spec))

_catalytic_topologies(r::EnzymeReaction) =
    _catalytic_topologies(_to_legacy_reaction(r))

_atoms_dict(r::EnzymeReaction, met::Symbol) =
    _atoms_dict(_to_legacy_reaction(r), met)

_bound_metabolites_at_forms(spec::MechanismSpec, r::EnzymeReaction) =
    _bound_metabolites_at_forms(spec, _to_legacy_reaction(r))

_expand_substrate_product_dead_ends(specs::Vector{MechanismSpec},
                                    r::EnzymeReaction) =
    _expand_substrate_product_dead_ends(specs, _to_legacy_reaction(r))

_expand_add_dead_end_regulator(spec::AbstractMechanismSpec, r::EnzymeReaction;
                               kwargs...) =
    _expand_add_dead_end_regulator(spec, _to_legacy_reaction(r); kwargs...)

_expand_to_allosteric(spec::AbstractMechanismSpec, r::EnzymeReaction) =
    _expand_to_allosteric(spec, _to_legacy_reaction(r))

_expand_add_allosteric_regulator(spec::AbstractMechanismSpec,
                                 r::EnzymeReaction) =
    _expand_add_allosteric_regulator(spec, _to_legacy_reaction(r))

_expand_change_allo_state(spec::AbstractMechanismSpec, r::EnzymeReaction) =
    _expand_change_allo_state(spec, _to_legacy_reaction(r))

expand_mechanisms(specs::Vector{<:AbstractMechanismSpec}, r::EnzymeReaction) =
    expand_mechanisms(specs, _to_legacy_reaction(r))

"""
    _assert_mechanism_invariants(m::Mechanism) -> Nothing

Structural invariants every valid Mechanism should satisfy:
- Every group is non-empty
- source_idx values are unique and dense (1 through n_steps)
- Each binding step's bound_metabolite is non-nothing AND iso steps have nothing
- from_species != to_species for every step
"""
function _assert_mechanism_invariants(m::Mechanism)
    flat = collect(Iterators.flatten(m.steps))
    isempty(flat) && error("empty steps in Mechanism")
    for g in m.steps
        isempty(g) && error("empty kinetic group in Mechanism")
    end
    src_indices = [s.source_idx for s in flat]
    sort(src_indices) == collect(1:length(flat)) ||
        error("source_idx values not dense 1..n: got $(sort(src_indices))")
    for s in flat
        if is_binding(s)
            s.bound_metabolite === nothing &&
                error("binding step has nothing bound_metabolite")
        else
            s.bound_metabolite === nothing ||
                error("iso step has non-nothing bound_metabolite")
        end
        s.from_species == s.to_species &&
            error("from_species == to_species in step $s")
    end
    nothing
end

function _assert_mechanism_invariants(m::AllostericMechanism)
    # Build a minimal Mechanism view of the catalytic side to reuse the
    # base invariant checks (every cat-group non-empty, source_idx dense
    # over flat cat steps, etc.). Then check the allosteric-specific
    # invariants against the actual AllostericMechanism fields.
    flat = Step[s for g in m.cat_steps for s in g]
    isempty(flat) && error("AllostericMechanism: empty cat_steps")
    for g in m.cat_steps
        isempty(g) && error("AllostericMechanism: empty catalytic kinetic group")
    end
    src_indices = [s.source_idx for s in flat]
    sort(src_indices) == collect(1:length(flat)) ||
        error("AllostericMechanism: cat_steps source_idx not dense 1..n")

    # cat_allo_states is one per cat group, validated by the constructor;
    # re-check defensively here:
    length(m.cat_allo_states) == length(m.cat_steps) ||
        error("AllostericMechanism: cat_allo_states length " *
              "$(length(m.cat_allo_states)) ≠ cat_steps length " *
              "$(length(m.cat_steps))")
    valid_cat_states = (:OnlyR, :EqualRT, :NonequalRT)
    for tag in m.cat_allo_states
        tag in valid_cat_states ||
            error("AllostericMechanism: invalid cat allo state $tag")
    end

    m.catalytic_multiplicity ≥ 1 ||
        error("AllostericMechanism: catalytic_multiplicity " *
              "$(m.catalytic_multiplicity) must be ≥ 1")

    # regulatory_sites is a Vector{RegulatorySite}; each site carries
    # its own ligand list + multiplicity + per-ligand allo states. The
    # constructor validates internal structure; here we only assert the
    # list is non-nothing.
    m.regulatory_sites isa Vector{RegulatorySite} ||
        error("AllostericMechanism: regulatory_sites not Vector{RegulatorySite}")

    nothing
end
