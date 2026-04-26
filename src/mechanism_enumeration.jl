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
    param_count::Int
end

"""
    AllostericMechanismSpec <: AbstractMechanismSpec

Allosteric enzyme mechanism built from a base `MechanismSpec`. Each
catalytic kinetic group and each regulator-site ligand carries a tag
indicating its R/T-state relationship. Tags: `:OnlyR`, `:OnlyT`,
`:EqualRT`, `:NonequalRT` (default for absent entries).

`group_tags` maps `kinetic_group → tag` (only non-default entries).
`reg_ligand_tags` maps `ligand → tag` (only non-default entries).
"""
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    group_tags::Dict{Int, Symbol}
    reg_ligand_tags::Dict{Symbol, Symbol}
    param_count::Int
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
    @nospecialize(reaction::EnzymeReaction),
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
    @nospecialize(reaction::EnzymeReaction),
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
                param_count = n_re + 2 * n_ss -
                    n_thermo + 2

                push!(result, MechanismSpec(
                    reaction, tagged, param_count
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
    @nospecialize(reaction::EnzymeReaction),
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
    @nospecialize(reaction::EnzymeReaction),
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

            # Compute param_count by counting kinetic groups
            groups_re = Set{Int}()
            groups_ss = Set{Int}()
            for s in new_steps
                if s.is_equilibrium
                    push!(groups_re, s.kinetic_group)
                else
                    push!(groups_ss, s.kinetic_group)
                end
            end
            n_forms = length(
                all_form_names(new_steps))
            n_thermo = length(new_steps) - n_forms + 1
            pc = length(groups_re) +
                 2 * length(groups_ss) - n_thermo + 2

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
function EnzymeMechanism(spec::MechanismSpec)
    rxn = spec.reaction
    subs = Tuple(s[1] for s in substrates(rxn))
    prods = Tuple(p[1] for p in products(rxn))
    regs = Tuple(regulators(rxn))
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
    @nospecialize(reaction::EnzymeReaction),
)
    topos = _catalytic_topologies(reaction)
    expanded = _expand_substrate_product_dead_ends(
        topos, reaction)

    n_s = length(substrates(reaction))
    n_p = length(products(reaction))
    floor_pc = n_s + n_p + 3

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
    pc = max(_param_count_from_steps(new_steps), floor_pc)
    MechanismSpec(spec.reaction, new_steps, pc)
end

"""
Compute parameter count from a step list, counting kinetic groups
rather than individual steps. This is exact for mechanisms whose
cycles all impose independent Wegscheider constraints; dead-end
mirror cycles add 0 effective constraints, so this formula can
underestimate — callers should floor to a safe lower bound.
"""
function _param_count_from_steps(steps::Vector{StepSpec})
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
    length(groups_re) + 2 * length(groups_ss) -
        n_thermo + 2
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

# Expansion moves are stubbed pending Phase 4.2 rewrite.

_expand_re_to_ss(::MechanismSpec) = MechanismSpec[]
_expand_re_to_ss(::AllostericMechanismSpec) = AllostericMechanismSpec[]

_expand_remove_constraint(::MechanismSpec) = MechanismSpec[]
_expand_remove_constraint(::AllostericMechanismSpec) = AllostericMechanismSpec[]

_expand_add_dead_end_regulator(
    ::MechanismSpec,
    @nospecialize(::EnzymeReaction);
    exclude_regs::Set{Symbol}=Set{Symbol}(),
) = MechanismSpec[]

_expand_add_dead_end_regulator(
    ::AllostericMechanismSpec,
    @nospecialize(::EnzymeReaction);
    exclude_regs::Set{Symbol}=Set{Symbol}(),
) = AllostericMechanismSpec[]

_expand_to_allosteric(
    ::MechanismSpec, @nospecialize(::EnzymeReaction),
) = AllostericMechanismSpec[]

_expand_to_allosteric(
    ::AllostericMechanismSpec, @nospecialize(::EnzymeReaction),
) = AllostericMechanismSpec[]

_expand_add_allosteric_regulator(
    ::AbstractMechanismSpec, @nospecialize(::EnzymeReaction),
) = AllostericMechanismSpec[]

_expand_remove_tr_equiv(
    ::AbstractMechanismSpec, @nospecialize(::EnzymeReaction),
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
    cm = EnzymeMechanism(spec.base)
    sorted_groups = sort(collect(spec.group_tags); by=first)
    group_tags = Tuple((g, t) for (g, t) in sorted_groups)
    cat_sites = (spec.catalytic_n, group_tags)

    reg_sites = Tuple(
        (Tuple(group),
         spec.allosteric_multiplicities[i],
         Tuple((l, spec.reg_ligand_tags[l])
               for l in group
               if haskey(spec.reg_ligand_tags, l)))
        for (i, group) in enumerate(spec.allosteric_reg_sites)
    )

    AllostericEnzymeMechanism(cm, cat_sites, reg_sites)
end

function _push_to_dict!(
    result::Dict{Int, Vector{AbstractMechanismSpec}},
    spec::AbstractMechanismSpec)
    push!(get!(result, spec.param_count,
        AbstractMechanismSpec[]), spec)
end

"""
    expand_mechanisms(specs, reaction) → Dict{Int, Vector{AbstractMechanismSpec}}

Apply all +1 and +2 expansion moves. Results grouped
by target param_count.
"""
function expand_mechanisms(
    specs::Vector{<:AbstractMechanismSpec},
    @nospecialize(reaction::EnzymeReaction))
    result = Dict{Int, Vector{AbstractMechanismSpec}}()
    for spec in specs
        _add_expansions!(result, spec, reaction)
    end
    result
end

function _add_expansions!(
    result::Dict{Int, Vector{AbstractMechanismSpec}},
    spec::AbstractMechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
    for s in _expand_re_to_ss(spec)
        _push_to_dict!(result, s)
    end
    for s in _expand_remove_constraint(spec)
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
    for s in _expand_remove_tr_equiv(spec, reaction)
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

Remove structural duplicates from each param_count bucket
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
