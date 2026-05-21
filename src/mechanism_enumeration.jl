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

"""Build a form name Symbol from sorted bound metabolite names."""
function _form_name(
    bound_subs::Vector{Symbol},
    bound_prods::Vector{Symbol},
    has_residual::Bool,
)
    parts = sort!(vcat(bound_subs, bound_prods))
    base = isempty(parts) ? "E" : "E_" * join(parts, "_")
    has_residual && isempty(parts) && (base = "Estar")
    has_residual && !isempty(parts) &&
        (base = "Estar_" * join(parts, "_"))
    Symbol(base)
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
`is_estar` controls form naming (Estar prefix).
`has_residual_atoms` tracks whether residual atoms remain.
"""
function _release_products!(
    all_paths, backtrack!,
    iso_form::Symbol,
    residual_atoms::Dict{Symbol,Int},
    consumed_subs::Vector{Symbol},
    released_prods::Vector{Symbol},
    prod_subset::Vector{Symbol},
    prod_atoms::Dict{Symbol,Dict{Symbol,Int}},
    is_estar::Bool,
    has_residual_atoms::Bool,
    steps::Vector{StepSpec},
)
    # Generate all release orderings of products
    function _release_recurse!(
        cur::Symbol,
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
            new_form = if isempty(new_unreleased)
                _form_name(
                    Symbol[], Symbol[], is_estar)
            else
                _form_name(
                    Symbol[], new_unreleased,
                    is_estar)
            end
            # Canonical: metabolite on LHS
            rel_step = StepSpec(
                [new_form, p], [cur], true, 0
            )
            push!(steps, rel_step)
            _release_recurse!(
                new_form,
                _subtract_atoms(
                    cur_atoms, prod_atoms[p]),
                new_unreleased,
                [rel_so_far; p],
            )
            pop!(steps)
        end
    end

    _release_recurse!(
        iso_form,
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

    # Collect all complete catalytic paths as step lists
    all_paths = Vector{Vector{StepSpec}}()

    # Backtracking state:
    # - cur_form: current enzyme form name (Symbol)
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
    # - steps: path of StepSpec accumulated so far
    function backtrack!(
        cur_form::Symbol,
        acc_atoms::Dict{Symbol,Int},
        consumed_subs::Vector{Symbol},
        released_prods::Vector{Symbol},
        on_enzyme_subs::Vector{Symbol},
        on_enzyme_prods::Vector{Symbol},
        # has_residual: true when enzyme is in Estar
        # conformation (ping-pong intermediate), even if
        # no residual atoms remain (empty-residual pp).
        # Controls both form naming (Estar prefix) and
        # branch selection.
        has_residual::Bool,
        post_final::Bool,
        steps::Vector{StepSpec},
    )
        # Check for complete cycle
        if cur_form == :E && !isempty(steps)
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
                new_form = _form_name(
                    Symbol[], new_on_prods, false
                )
                # Canonical: metabolite on LHS
                step = StepSpec(
                    [new_form, p], [cur_form], true, 0
                )
                push!(steps, step)
                backtrack!(
                    new_form,
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
                new_form = _form_name(
                    new_on, Symbol[], false
                )
                step = StepSpec(
                    [cur_form, s], [new_form], true, 0
                )
                push!(steps, step)
                backtrack!(
                    new_form,
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
                    new_form = _form_name(
                        new_on, Symbol[], false
                    )
                    step = StepSpec(
                        [cur_form, s], [new_form], true, 0
                    )
                    push!(steps, step)
                    backtrack!(
                        new_form,
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
                        iso_form = _form_name(
                            Symbol[],
                            collect(prod_subset),
                            true
                        )
                        step = StepSpec(
                            [cur_form], [iso_form], true, 0
                        )
                        push!(steps, step)
                        # Release products one at a time
                        # Ping-pong from free enzyme
                        # always creates Estar
                        _release_products!(
                            all_paths, backtrack!,
                            iso_form, residual,
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
                    new_form = _form_name(
                        Symbol[],
                        copy(remaining_prods), false
                    )
                    step = StepSpec(
                        [cur_form], [new_form], true, 0
                    )
                    push!(steps, step)
                    backtrack!(
                        new_form, acc_atoms,
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
                new_form = _form_name(
                    new_on, Symbol[], true
                )
                step = StepSpec(
                    [cur_form, s], [new_form], true, 0
                )
                push!(steps, step)
                backtrack!(
                    new_form,
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
                    new_form = _form_name(
                        new_on, Symbol[], true
                    )
                    step = StepSpec(
                        [cur_form, s], [new_form], true, 0
                    )
                    push!(steps, step)
                    backtrack!(
                        new_form,
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
                        new_form = _form_name(
                            Symbol[],
                            copy(remaining_prods),
                            false
                        )
                        step = StepSpec(
                            [cur_form], [new_form],
                            true, 0
                        )
                        push!(steps, step)
                        backtrack!(
                            new_form, acc_atoms,
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
                        iso_form = _form_name(
                            Symbol[],
                            collect(prod_subset),
                            has_more
                        )
                        step = StepSpec(
                            [cur_form], [iso_form],
                            true, 0
                        )
                        push!(steps, step)
                        _release_products!(
                            all_paths, backtrack!,
                            iso_form, residual_atoms,
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

    backtrack!(
        :E, Dict{Symbol,Int}(), Symbol[], Symbol[],
        Symbol[], Symbol[], false, false, StepSpec[]
    )

    isempty(all_paths) && return MechanismSpec[]

    # Deduplicate paths by their step content (as sets)
    StepKey = Tuple{Vector{Symbol}, Vector{Symbol}}
    _step_key(s::StepSpec) = (
        sort(s.reactants), sort(s.products)
    )::StepKey
    unique_paths = Vector{Vector{StepSpec}}()
    seen_path_keys = Set{Set{StepKey}}()
    for path in all_paths
        key = Set(_step_key(s) for s in path)
        key ∈ seen_path_keys && continue
        push!(seen_path_keys, key)
        push!(unique_paths, path)
    end

    # --- Group paths by isomerization pattern ---
    met_names = Set{Symbol}(sub_names)
    union!(met_names, prod_names)

    function _iso_pattern(path)
        Set(
            (_step_key(s) for s in path
             if length(s.reactants) == 1 &&
                length(s.products) == 1 &&
                s.reactants[1] ∉ met_names &&
                s.products[1] ∉ met_names)
        )
    end

    iso_groups = Dict{
        Set{StepKey}, Vector{Vector{StepSpec}}
    }()
    for path in unique_paths
        pat = _iso_pattern(path)
        paths_vec = get!(iso_groups, pat,
            Vector{Vector{StepSpec}}())
        push!(paths_vec, path)
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

    function _steps_for_ordering(
        all_group_steps, ordering, met_set,
    )
        selected = Set{StepKey}()
        accessible = Set{Symbol}()
        for level in ordering
            level_set = Set{Symbol}(level)
            allowed = union(accessible, level_set)
            for met in level
                for sk in all_group_steps
                    lhs_mets = [
                        s for s in sk[1] if s ∈ met_names
                    ]
                    length(lhs_mets) == 1 || continue
                    lhs_mets[1] == met || continue
                    form = [
                        s for s in sk[1] if s ∉ met_names
                    ]
                    length(form) == 1 || continue
                    form_met_parts = [
                        Symbol(p) for p in split(
                            replace(string(form[1]),
                                    "Estar" => "E"),
                            "_")
                        if Symbol(p) ∈ met_set
                    ]
                    if all(m ∈ allowed
                           for m in form_met_parts)
                        push!(selected, sk)
                    end
                end
            end
            union!(accessible, level)
        end
        selected
    end

    # --- Build topologies ---
    result = MechanismSpec[]
    for (iso_pat, group_paths) in iso_groups
        all_group_steps = Set{StepKey}()
        step_dict = Dict{StepKey, StepSpec}()
        for path in group_paths
            for s in path
                sk = _step_key(s)
                push!(all_group_steps, sk)
                step_dict[sk] = s
            end
        end

        # Always include all isomerization steps
        iso_keys = Set{StepKey}()
        for sk in all_group_steps
            lhs_mets = [s for s in sk[1] if s ∈ met_names]
            rhs_mets = [s for s in sk[2] if s ∈ met_names]
            if isempty(lhs_mets) && isempty(rhs_mets)
                push!(iso_keys, sk)
            end
        end

        sub_met_set = Set(sub_names)
        prod_met_set = Set(prod_names)

        sub_binding_mets = Set{Symbol}()
        prod_binding_mets = Set{Symbol}()
        for sk in all_group_steps
            lhs_mets = [
                s for s in sk[1] if s ∈ met_names
            ]
            length(lhs_mets) == 1 || continue
            met = lhs_mets[1]
            if met ∈ sub_met_set
                push!(sub_binding_mets, met)
            else
                push!(prod_binding_mets, met)
            end
        end

        sub_orderings = _weak_orderings(
            sort(collect(sub_binding_mets)))
        prod_orderings = _weak_orderings(
            sort(collect(prod_binding_mets)))

        seen_topos = Set{Set{StepKey}}()
        for sub_ord in sub_orderings
            for prod_ord in prod_orderings
                sub_keys = _steps_for_ordering(
                    all_group_steps, sub_ord, sub_met_set,
                )
                prod_keys = _steps_for_ordering(
                    all_group_steps, prod_ord, prod_met_set,
                )
                topo_keys = union(iso_keys, sub_keys,
                    prod_keys)
                topo_keys ∈ seen_topos && continue
                push!(seen_topos, topo_keys)

                steps = [step_dict[sk] for sk in topo_keys]
                sort!(steps; by=s -> (
                    length(s.reactants) == 1 ? 1 : 0,
                    join(sort(s.reactants), "_")
                ))

                iso_idx = findfirst(
                    s -> length(s.reactants) == 1, steps
                )
                tagged = [
                    StepSpec(
                        s.reactants, s.products,
                        i != iso_idx, i
                    )
                    for (i, s) in enumerate(steps)
                ]

                form_names = Set{Symbol}()
                for s in tagged
                    union!(form_names, s.reactants)
                    union!(form_names, s.products)
                end
                setdiff!(form_names, met_names)
                n_forms = length(form_names)
                n_steps = length(tagged)
                n_cycles = n_steps - n_forms + 1
                n_re = count(s -> s.is_equilibrium, tagged)
                n_ss = n_steps - n_re
                n_thermo = n_cycles
                n_fit_params_estimate = n_re + 2 * n_ss -
                    n_thermo

                push!(result, MechanismSpec(
                    reaction, tagged, n_fit_params_estimate
                ))
            end
        end
    end
    result
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

Convert a `MechanismSpec` to an `EnzymeMechanism`, or an
`AllostericMechanismSpec` to an `AllostericEnzymeMechanism`.
"""
compile_mechanism(spec::MechanismSpec) =
    EnzymeMechanism(spec)
compile_mechanism(spec::AllostericMechanismSpec) =
    AllostericEnzymeMechanism(spec)

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
    init_mechanisms(reaction) -> Vector{MechanismSpec}

Produce all mechanisms at minimum parameter count for a reaction.
For each catalytic topology: 1 SS step, all K's grouped equal
per metabolite, all substrate/product dead-end subsets.

Step kinetic groups are reassigned so that all binding steps
sharing the same `(metabolite, RE/SS)` class collapse into one
group; iso steps and uncollapsable bindings stay singletons.
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
# the parametric Legacy form. The enumeration primitives still dispatch
# on EnzymeReactionLegacy; Stage 5 rewrites the primitives to consume
# EnzymeReaction directly, retiring these forwarders.
init_mechanisms(r::EnzymeReaction) = init_mechanisms(_to_legacy_reaction(r))

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
