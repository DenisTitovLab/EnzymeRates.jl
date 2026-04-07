# ABOUTME: Mechanism enumeration by incremental parameter count growth
# ABOUTME: Provides init_mechanisms, expand_mechanisms, dedup! building blocks

# ─── Data Types ─────────────────────────────────────────────

"""
Constraint on kinetic parameters: target parameter equals a linear
combination of source parameters.
Format: `(target_sym, coeff, [(src_sym, src_coeff), ...])`.
"""
const ParamConstraint = Tuple{Symbol, Int, Vector{Tuple{Symbol, Int}}}

# ─── Mechanism Spec Types ──────────────────────────────────────

abstract type AbstractMechanismSpec end

"""Elementary step in canonical binding direction (metabolite on LHS)."""
struct StepSpec
    reactants::Vector{Symbol}   # [:E, :S] or [:EAB]
    products::Vector{Symbol}    # [:ES] or [:EPQ]
    is_equilibrium::Bool
end

Base.:(==)(a::StepSpec, b::StepSpec) =
    a.reactants == b.reactants &&
    a.products == b.products &&
    a.is_equilibrium == b.is_equilibrium

Base.hash(s::StepSpec, h::UInt) =
    hash(s.is_equilibrium,
        hash(s.products, hash(s.reactants, h)))

"""
    MechanismSpec <: AbstractMechanismSpec

Represents a monomeric enzyme mechanism specification in the
staged enumeration pipeline. Steps are `StepSpec` values with
inline form names and equilibrium status.
"""
struct MechanismSpec <: AbstractMechanismSpec
    reaction::Any
    steps::Vector{StepSpec}
    param_constraints::Vector{ParamConstraint}
    param_count::Int
end

"""
    AllostericMechanismSpec <: AbstractMechanismSpec

Represents an allosteric enzyme mechanism built from a
base `MechanismSpec` plus allosteric site and multiplicity info.
`tr_equiv_metabolites` lists metabolites with K_T = K_R.
`tr_equiv_cat_steps` lists indices of non-binding SS steps
with kf_T = kf_R (catalytic step TR equivalence).
`r_only_metabolites` lists metabolites absent from T-state (K_T = ∞).
`t_only_metabolites` lists metabolites absent from R-state (K_R = ∞).
`r_only_cat_steps` lists non-binding SS step indices where T-state
doesn't catalyze (kf_T = 0).
"""
struct AllostericMechanismSpec <: AbstractMechanismSpec
    base::MechanismSpec
    catalytic_n::Int
    allosteric_reg_sites::Vector{Vector{Symbol}}
    allosteric_multiplicities::Vector{Int}
    tr_equiv_metabolites::Vector{Symbol}
    tr_equiv_cat_steps::Vector{Int}
    r_only_metabolites::Vector{Symbol}
    t_only_metabolites::Vector{Symbol}
    r_only_cat_steps::Vector{Int}
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
                    [new_form, p], [cur_form], true
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
                    [cur_form, s], [new_form], true
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
            # Option 1: bind another substrate
            for s in remaining_subs
                new_on = [on_enzyme_subs; s]
                new_consumed = [consumed_subs; s]
                new_form = _form_name(
                    new_on, Symbol[], false
                )
                step = StepSpec(
                    [cur_form, s], [new_form], true
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
            # Option 2: ping-pong isomerize
            if !isempty(remaining_subs)
                for p in remaining_prods
                    _can_pingpong(
                        acc_atoms, prod_atoms[p]
                    ) || continue
                    residual = _subtract_atoms(
                        acc_atoms, prod_atoms[p]
                    )
                    iso_form = _form_name(
                        on_enzyme_subs, [p], true
                    )
                    step = StepSpec(
                        [cur_form], [iso_form], true
                    )
                    push!(steps, step)
                    rel_form = _form_name(
                        Symbol[], Symbol[], true
                    )
                    rel_step = StepSpec(
                        [rel_form, p], [iso_form], true
                    )
                    push!(steps, rel_step)
                    backtrack!(
                        rel_form, residual,
                        consumed_subs,
                        [released_prods; p],
                        Symbol[], Symbol[],
                        true, false, steps
                    )
                    pop!(steps)
                    pop!(steps)
                end
            end
            # Option 3: final isomerize (all subs bound)
            if isempty(remaining_subs)
                all_prod_atoms = reduce(
                    _add_atoms, values(prod_atoms);
                    init=Dict{Symbol,Int}()
                )
                if _can_pingpong(
                    acc_atoms, all_prod_atoms
                )
                    new_form = _form_name(
                        Symbol[],
                        copy(remaining_prods), false
                    )
                    step = StepSpec(
                        [cur_form], [new_form], true
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
            # Residual only (E*): bind any remaining sub
            for s in remaining_subs
                new_on = [s]
                new_consumed = [consumed_subs; s]
                new_form = _form_name(
                    new_on, Symbol[], true
                )
                step = StepSpec(
                    [cur_form, s], [new_form], true
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
            # Final isomerize: all subs consumed, release
            # remaining products
            if isempty(remaining_subs) &&
                    !isempty(remaining_prods)
                all_prod_atoms = reduce(
                    _add_atoms, (
                        prod_atoms[p]
                        for p in remaining_prods
                    );
                    init=Dict{Symbol,Int}()
                )
                if _can_pingpong(
                    acc_atoms, all_prod_atoms
                )
                    new_form = _form_name(
                        Symbol[],
                        copy(remaining_prods), false
                    )
                    step = StepSpec(
                        [cur_form], [new_form], true
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
        elseif !isempty(on_enzyme_subs) && has_residual
            # Residual + substrates bound
            # Option 1: bind another remaining substrate
            for s in remaining_subs
                new_on = [on_enzyme_subs; s]
                new_consumed = [consumed_subs; s]
                new_form = _form_name(
                    new_on, Symbol[], true
                )
                step = StepSpec(
                    [cur_form, s], [new_form], true
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
            # Option 2: isomerize to release a product
            for p in remaining_prods
                _can_pingpong(
                    acc_atoms, prod_atoms[p]
                ) || continue
                residual_atoms = _subtract_atoms(
                    acc_atoms, prod_atoms[p]
                )
                has_more = !isempty(residual_atoms)
                if has_more ||
                        !isempty(remaining_subs)
                    # Ping-pong: isomerize + release
                    iso_form = _form_name(
                        on_enzyme_subs, [p], has_more
                    )
                    step = StepSpec(
                        [cur_form], [iso_form], true
                    )
                    push!(steps, step)
                    rel_form = _form_name(
                        Symbol[], Symbol[], has_more
                    )
                    rel_step = StepSpec(
                        [rel_form, p],
                        [iso_form], true
                    )
                    push!(steps, rel_step)
                    backtrack!(
                        rel_form, residual_atoms,
                        consumed_subs,
                        [released_prods; p],
                        Symbol[], Symbol[],
                        has_more, false, steps
                    )
                    pop!(steps)
                    pop!(steps)
                end
                if !has_more &&
                        isempty(remaining_subs)
                    # Final: all subs consumed, all
                    # residual consumed — release all
                    # remaining products
                    new_form = _form_name(
                        Symbol[],
                        copy(remaining_prods), false
                    )
                    step = StepSpec(
                        [cur_form], [new_form], true
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
                        i != iso_idx
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
                    reaction, tagged,
                    ParamConstraint[], param_count
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
        n_de = length(de_form_names)

        # Map each dead-end form to its bound metabolites
        de_bound = Dict{Symbol, Set{Symbol}}()
        for de_name in de_form_names
            entries = de_forms[de_name]
            f, m = first(entries)
            de_bound[de_name] = union(
                bound[f], Set([m]))
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

            # Add binding steps for active dead-ends
            for de_name in sort(collect(active_de))
                entries = de_forms[de_name]
                for (cat_form, met) in entries
                    # [cat_form, met] → [de_name]
                    # (always RE)
                    push!(new_steps, StepSpec(
                        [cat_form, met],
                        [de_name], true))
                end
            end

            # Add mirror steps: for each catalytic
            # step, if both endpoints have dead-end
            # forms with the same metabolite, add a
            # parallel step. Mirror inherits RE/SS.
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
                    if met !== nothing
                        push!(new_steps, StepSpec(
                            [from_de, met],
                            [to_de],
                            s.is_equilibrium))
                    else
                        push!(new_steps, StepSpec(
                            [from_de], [to_de],
                            s.is_equilibrium))
                    end
                end
            end

            # Compute param_count
            n_steps = length(new_steps)
            n_re = count(
                s -> s.is_equilibrium, new_steps)
            n_ss = n_steps - n_re
            n_forms = length(
                all_form_names(new_steps))
            n_thermo = n_steps - n_forms + 1
            pc = n_re + 2 * n_ss - n_thermo + 2

            push!(result, MechanismSpec(
                spec.reaction, new_steps,
                spec.param_constraints, pc))
        end
    end
    result
end

# ─── Compilation ─────────────────────────────────────────────

"""Construct EnzymeMechanism from MechanismSpec."""
EnzymeMechanism(spec::MechanismSpec) =
    _compile_enzyme_mechanism(spec)

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

function _compile_enzyme_mechanism(spec::MechanismSpec)
    rxn = spec.reaction

    # Collect metabolite names and atom dicts
    met_atoms = Dict{Symbol,Dict{Symbol,Int}}()
    for (name, atoms) in substrates(rxn)
        met_atoms[name] = Dict{Symbol,Int}(
            a => c for (a, c) in atoms
        )
    end
    for (name, atoms) in products(rxn)
        met_atoms[name] = Dict{Symbol,Int}(
            a => c for (a, c) in atoms
        )
    end
    met_set = Set(keys(met_atoms))
    for r in regulators(rxn)
        push!(met_set, r)
        # Don't overwrite existing atom dict (e.g., substrate
        # that also appears as a regulator)
        haskey(met_atoms, r) || (met_atoms[r] = Dict{Symbol,Int}())
    end

    # Strip __regN suffixes from metabolite names
    function _clean_met(sym::Symbol)
        s = string(sym)
        m = match(r"^(.+)__reg\d*$", s)
        m !== nothing ? Symbol(m.captures[1]) : sym
    end

    # Add suffixed regulator names to met_set/met_atoms
    # so BFS recognizes them as metabolites
    for s in spec.steps
        for sym in Iterators.flatten(
                (s.reactants, s.products))
            clean = _clean_met(sym)
            if clean != sym && clean ∈ met_set
                push!(met_set, sym)
                met_atoms[sym] = met_atoms[clean]
            end
        end
    end

    # Collect form names (all symbols in steps
    # that are not metabolites)
    form_set = Set{Symbol}()
    for s in spec.steps
        for sym in s.reactants
            sym ∉ met_set && push!(form_set, sym)
        end
        for sym in s.products
            sym ∉ met_set && push!(form_set, sym)
        end
    end

    # BFS from :E to compute atoms for each form
    form_atoms = Dict{Symbol,Dict{Symbol,Int}}(
        :E => Dict{Symbol,Int}()
    )
    # Build adjacency: form -> [(neighbor, met_or_nothing,
    #   direction)] where direction = :add if metabolite
    #   binds going from form to neighbor
    adj = Dict{Symbol,
        Vector{Tuple{Symbol,Union{Nothing,Symbol},
                      Symbol}}}()
    for s in spec.steps
        from = s.reactants[1]
        to = s.products[1]
        met = length(s.reactants) == 2 ?
            s.reactants[2] : nothing
        # Canonical: met on LHS means binding
        # from→to. Reverse direction = release.
        if !haskey(adj, from)
            adj[from] = Tuple{
                Symbol,Union{Nothing,Symbol},Symbol
            }[]
        end
        if !haskey(adj, to)
            adj[to] = Tuple{
                Symbol,Union{Nothing,Symbol},Symbol
            }[]
        end
        push!(adj[from], (to, met, :add))
        push!(adj[to], (from, met, :subtract))
    end

    queue = Symbol[:E]
    while !isempty(queue)
        cur = popfirst!(queue)
        cur_atoms = form_atoms[cur]
        haskey(adj, cur) || continue
        for (nbr, met, dir) in adj[cur]
            haskey(form_atoms, nbr) && continue
            nbr_atoms = copy(cur_atoms)
            if met !== nothing
                ma = met_atoms[met]
                if dir == :add
                    for (a, c) in ma
                        nbr_atoms[a] =
                            get(nbr_atoms, a, 0) + c
                    end
                else  # :subtract
                    for (a, c) in ma
                        nbr_atoms[a] =
                            get(nbr_atoms, a, 0) - c
                        nbr_atoms[a] == 0 &&
                            delete!(nbr_atoms, a)
                    end
                end
            end
            # Isomerization: no atom change
            form_atoms[nbr] = nbr_atoms
            push!(queue, nbr)
        end
    end

    # Build form name → cleaned name mapping, preserving
    # original name when cleaning would produce a collision
    form_names = sort!(collect(form_set))
    form_name_map = Dict{Symbol,Symbol}()
    used_clean = Set{Symbol}()
    for name in form_names
        candidate = _clean_met(name)
        if candidate == name || candidate ∉ used_clean
            form_name_map[name] = candidate
            push!(used_clean, candidate)
        else
            form_name_map[name] = name
            push!(used_clean, name)
        end
    end

    enzymes = Tuple(
        (form_name_map[name], Tuple(sort!(
            [Tuple(p) for p in form_atoms[name]];
            by=first
        )))
        for name in form_names
    )

    species = (
        substrates(rxn), products(rxn),
        regulators(rxn), enzymes
    )

    reactions = Tuple(
        let r = s.reactants, p = s.products
            lhs = Tuple(
                haskey(form_name_map, x) ?
                    form_name_map[x] : _clean_met(x)
                for x in r)
            rhs = Tuple(
                haskey(form_name_map, x) ?
                    form_name_map[x] : _clean_met(x)
                for x in p)
            (lhs, rhs)
        end
        for s in spec.steps
    )

    eq_steps = Tuple(s.is_equilibrium for s in spec.steps)

    constraints = Tuple(
        (t, c, Tuple(Tuple.(f)))
        for (t, c, f) in spec.param_constraints
    )

    EnzymeMechanism(species, reactions, eq_steps,
        constraints)
end

# ─── Mechanism Enumeration ───────────────────────────────────

"""
    init_mechanisms(reaction) -> Vector{MechanismSpec}

Produce all mechanisms at minimum parameter count for a reaction.
For each catalytic topology: 1 SS step, all K's constrained equal
per metabolite, all substrate/product dead-end subsets (2^n).
"""
function init_mechanisms(
    @nospecialize(reaction::EnzymeReaction),
)
    topos = _catalytic_topologies(reaction)
    expanded = _expand_substrate_product_dead_ends(
        topos, reaction)

    n_s = length(substrates(reaction))
    n_p = length(products(reaction))
    min_pc = n_s + n_p + 3

    result = MechanismSpec[]
    for spec in expanded
        constraints = _max_equivalence_constraints(spec)
        # For specs with Estar forms (ping-pong), the
        # isomerization step adds an independent K that
        # is not reducible by equivalence constraints.
        # The +1 accounts for possible overlap between
        # equivalence and thermodynamic constraints.
        has_estar = any(all_form_names(spec.steps)) do f
            _is_estar_form(f)
        end
        pc = has_estar ?
            spec.param_count -
                length(constraints) + 1 :
            min_pc
        # Ensure upper bound invariant
        pc = max(pc, min_pc)
        push!(result, MechanismSpec(
            spec.reaction, spec.steps,
            constraints, pc))
    end
    result
end

"""
Build equivalence constraints for all groups of steps binding
the same metabolite with the same RE/SS status. Each group's
steps beyond the first are constrained to equal the first.
"""
function _max_equivalence_constraints(spec::MechanismSpec)
    # Group step indices by (metabolite, RE/SS)
    groups = Dict{Tuple{Symbol,Bool}, Vector{Int}}()
    for (i, s) in enumerate(spec.steps)
        met = step_metabolite(s)
        met === nothing && continue
        key = (met, s.is_equilibrium)
        push!(get!(groups, key, Int[]), i)
    end

    constraints = ParamConstraint[]
    for (_, g) in groups
        length(g) >= 2 || continue
        sort!(g)
        is_re = spec.steps[g[1]].is_equilibrium
        if is_re
            for j in 2:length(g)
                push!(constraints, (
                    Symbol("K$(g[j])"),
                    1,
                    [(Symbol("K$(g[1])"), 1)]
                ))
            end
        else
            for j in 2:length(g)
                for sfx in ("f", "r")
                    push!(constraints, (
                        Symbol("k$(g[j])$sfx"),
                        1,
                        [(Symbol("k$(g[1])$sfx"),
                          1)]
                    ))
                end
            end
        end
    end
    constraints
end

"""
    _step_index_from_constraint_sym(sym) -> Int or nothing

Parse the step index from a constraint parameter symbol like `K3`, `k3f`, `k3r`.
Returns nothing if the symbol doesn't match the expected pattern.
"""
function _step_index_from_constraint_sym(sym::Symbol)
    s = string(sym)
    m = match(r"^[Kk](\d+)", s)
    m === nothing && return nothing
    cap = m.captures[1]
    cap === nothing && return nothing
    parse(Int, cap)
end

"""Return the set of step indices involved in any param constraint."""
function _constrained_step_indices(constraints::Vector{ParamConstraint})
    idxs = Set{Int}()
    for (target, _, followers) in constraints
        idx = _step_index_from_constraint_sym(target)
        idx !== nothing && push!(idxs, idx)
        for (src, _) in followers
            sidx = _step_index_from_constraint_sym(src)
            sidx !== nothing && push!(idxs, sidx)
        end
    end
    idxs
end

"""
    _expand_re_to_ss(spec::MechanismSpec) → Vector{MechanismSpec}

Convert one RE step to SS. Skip constrained RE steps.
Mirror dead-end steps inherit the new SS status.
"""
function _expand_re_to_ss(spec::MechanismSpec)
    result = MechanismSpec[]
    constrained = _constrained_step_indices(spec.param_constraints)

    for (i, s) in enumerate(spec.steps)
        s.is_equilibrium || continue
        i in constrained && continue

        new_steps = [StepSpec(st.reactants, st.products, st.is_equilibrium)
                     for st in spec.steps]
        new_steps[i] = StepSpec(s.reactants, s.products, false)

        # Propagate SS to dead-end mirror steps (skip constrained steps)
        from_form, to_form = step_forms(s)
        n_mirrors = 0
        for (j, ms) in enumerate(new_steps)
            j == i && continue
            ms.is_equilibrium || continue
            j in constrained && continue
            mf, mt = step_forms(ms)
            if _is_mirror_of(mf, mt, from_form, to_form, spec.steps)
                new_steps[j] = StepSpec(ms.reactants, ms.products, false)
                n_mirrors += 1
            end
        end

        push!(result, MechanismSpec(
            spec.reaction, new_steps,
            copy(spec.param_constraints),
            spec.param_count + 1 + n_mirrors))
    end
    result
end

"""
    _expand_remove_constraint(spec::MechanismSpec) → Vector{MechanismSpec}

Remove one equivalence constraint group (+1 param for an RE K
constraint, +2 params for an SS kf/kr pair) per result.
SS kf and kr constraints for the same step are always paired and
removed together to avoid biochemically invalid mechanisms where
forward and reverse rates share a K while the other does not.
"""
function _expand_remove_constraint(spec::MechanismSpec)
    n = length(spec.param_constraints)
    # Identify kf/kr pairs: find each "f"-suffixed constraint
    # and its matching "r"-suffixed partner.
    paired = Set{Int}()
    pairs = Tuple{Int,Int}[]
    for i in 1:n
        i in paired && continue
        target_str = string(spec.param_constraints[i][1])
        endswith(target_str, "f") || continue
        kr_sym = Symbol(target_str[1:end-1] * "r")
        for j in (i+1):n
            j in paired && continue
            if spec.param_constraints[j][1] == kr_sym
                push!(pairs, (i, j))
                push!(paired, i)
                push!(paired, j)
                break
            end
        end
    end

    result = MechanismSpec[]
    # Remove each kf/kr pair together (+2 params)
    for (i, j) in pairs
        new_constraints = [
            spec.param_constraints[k]
            for k in 1:n if k != i && k != j]
        push!(result, MechanismSpec(
            spec.reaction, copy(spec.steps),
            new_constraints,
            spec.param_count + 2))
    end
    # Remove each unpaired constraint (+1 param)
    for i in 1:n
        i in paired && continue
        new_constraints = [
            spec.param_constraints[k]
            for k in 1:n if k != i]
        push!(result, MechanismSpec(
            spec.reaction, copy(spec.steps),
            new_constraints,
            spec.param_count + 1))
    end
    result
end

"""
    _is_mirror_of(mf, mt, from, to, steps) -> Bool

Check if (mf, mt) is a dead-end mirror of the catalytic step (from, to).
A mirror step connects dead-end forms that extend the catalytic endpoints
by binding the same extra metabolite.
"""
function _is_mirror_of(
    mf::Symbol, mt::Symbol,
    from::Symbol, to::Symbol,
    steps::Vector{StepSpec},
)
    # For (mf, mt) to be a mirror of (from, to):
    # there must be a binding step [from, met] → [mf] and
    # a binding step [to, met] → [mt] for the same metabolite.
    from_met = nothing
    to_met = nothing
    for s in steps
        f, t = step_forms(s)
        m = step_metabolite(s)
        m === nothing && continue
        if f == from && t == mf
            from_met = m
        elseif f == to && t == mt
            to_met = m
        end
    end
    from_met !== nothing && to_met !== nothing && from_met == to_met
end

"""
    _expand_add_dead_end_regulator(spec, reaction; exclude_regs)
        → Vector{MechanismSpec}

Add a new dead-end regulator to non-empty subsets of eligible
forms. Each variant adds +1 param (one new K, constrained equal
across all binding sites for this regulator).
"""
function _expand_add_dead_end_regulator(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction);
    exclude_regs::Set{Symbol}=Set{Symbol}(),
)
    roles = regulator_roles(reaction)
    isempty(roles) && return MechanismSpec[]

    # Find regulators not yet in mechanism
    existing_mets = Set{Symbol}()
    for s in spec.steps
        for sym in Iterators.flatten(
                (s.reactants, s.products))
            push!(existing_mets, sym)
        end
    end

    eligible_regs = Symbol[]
    for (name, role) in roles
        (role == :unknown || role == :dead_end) ||
            continue
        name in exclude_regs && continue
        reg_prefix = string(name) * "__reg"
        already = any(
            contains(string(m), reg_prefix)
            for m in existing_mets)
        already && continue
        push!(eligible_regs, name)
    end
    sort!(eligible_regs)

    isempty(eligible_regs) && return MechanismSpec[]

    sub_names = Set(s[1] for s in substrates(reaction))
    prod_names = Set(
        p[1] for p in products(reaction))
    bound = _bound_metabolites_at_forms(spec, reaction)
    cat_forms = all_form_names(spec)

    result = MechanismSpec[]

    for reg in eligible_regs
        dummy = Symbol(string(reg) * "__reg")

        # Eligible: neither all subs nor all prods bound
        eligible_forms = Symbol[]
        for f in sort(collect(cat_forms))
            haskey(bound, f) || continue
            fb = bound[f]
            (intersect(fb, sub_names) == sub_names ||
                intersect(fb, prod_names) ==
                    prod_names) && continue
            push!(eligible_forms, f)
        end

        isempty(eligible_forms) && continue

        # Find existing inhibitor dummies in steps.
        # Only metabolite dummies (not form names).
        existing_inhibitors = Symbol[]
        for s in spec.steps
            met = step_metabolite(s)
            met === nothing && continue
            s_met = string(met)
            contains(s_met, "__reg") || continue
            met === dummy && continue
            push!(existing_inhibitors, met)
        end
        sort!(unique!(existing_inhibitors))

        # Enumerate inhibitor competition patterns
        inh_patterns =
            _inhibitor_competition_patterns(
                sub_names, prod_names,
                existing_inhibitors)
        seen = Set{Vector{Symbol}}()

        for (comp_subs, comp_prods,
                comp_inhibitors) in inh_patterns
            # Find forms where competing metabolites
            # have binding steps (topology-aware)
            target_forms = Set{Symbol}()
            for met in comp_subs
                union!(target_forms,
                    _forms_with_binding_step(
                        spec.steps, met))
            end
            for met in comp_prods
                union!(target_forms,
                    _forms_with_binding_step(
                        spec.steps, met))
            end
            for inh in comp_inhibitors
                union!(target_forms,
                    _forms_with_binding_step(
                        spec.steps, inh))
            end

            # Filter: eligible AND not already
            # bound to any competing metabolite
            all_competing = union(
                comp_subs, comp_prods,
                comp_inhibitors)
            active = Symbol[]
            for f in sort(collect(target_forms))
                f in eligible_forms || continue
                haskey(bound, f) || continue
                isempty(intersect(
                    bound[f], all_competing)) ||
                    continue
                push!(active, f)
            end

            isempty(active) && continue
            active in seen && continue
            push!(seen, active)

            new_steps = copy(spec.steps)
            de_form_map = Dict{Symbol, Symbol}()

            # Add binding steps (always RE)
            binding_step_indices = Int[]
            for cf in active
                de_name = _dead_end_form_name(
                    cf, bound[cf], dummy)
                de_form_map[cf] = de_name
                push!(new_steps, StepSpec(
                    [cf, dummy], [de_name], true))
                push!(binding_step_indices,
                    length(new_steps))
            end

            # Add mirror steps for catalytic steps
            # whose both endpoints have dead-end forms
            for s in spec.steps
                from, to = step_forms(s)
                haskey(de_form_map, from) || continue
                haskey(de_form_map, to) || continue
                met = step_metabolite(s)
                from_de = de_form_map[from]
                to_de = de_form_map[to]
                if met !== nothing
                    push!(new_steps, StepSpec(
                        [from_de, met], [to_de],
                        s.is_equilibrium))
                else
                    push!(new_steps, StepSpec(
                        [from_de], [to_de],
                        s.is_equilibrium))
                end
            end

            # Equivalence constraints: all K's equal
            # for this regulator across binding sites
            new_constraints = copy(
                spec.param_constraints)
            if length(binding_step_indices) >= 2
                first_idx = binding_step_indices[1]
                for j in 2:length(
                        binding_step_indices)
                    push!(new_constraints, (
                        Symbol("K$(binding_step_indices[j])"),
                        1,
                        [(Symbol("K$(first_idx)"),
                            1)]))
                end
            end

            push!(result, MechanismSpec(
                spec.reaction, new_steps,
                new_constraints,
                spec.param_count + 1))
        end
    end
    result
end

function _expand_re_to_ss(spec::AllostericMechanismSpec)
    [_rewrap_allosteric(spec, new_base)
     for new_base in _expand_re_to_ss(spec.base)]
end

function _expand_remove_constraint(
    spec::AllostericMechanismSpec)
    [_rewrap_allosteric(spec, new_base)
     for new_base in _expand_remove_constraint(spec.base)]
end

function _expand_add_dead_end_regulator(
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction);
    exclude_regs::Set{Symbol}=Set{Symbol}())
    allo_regs = Set{Symbol}()
    for site in spec.allosteric_reg_sites
        for lig in site
            push!(allo_regs, lig)
        end
    end
    all_excluded = union(exclude_regs, allo_regs)
    [_rewrap_allosteric(spec, new_base)
     for new_base in _expand_add_dead_end_regulator(
         spec.base, reaction;
         exclude_regs=all_excluded)]
end

"""
    _valid_allosteric_differentiations(reaction, spec)

Enumerate biochemically valid T/R differentiations.
K-type: ≥1 substrate + ≥1 product absent from T-state.
V-type: all SS isomerization steps inactive in T-state.
Only r_only (T is inactive conformation).
"""
function _valid_allosteric_differentiations(
    @nospecialize(reaction::EnzymeReaction),
    spec::MechanismSpec)
    sub_names = [s[1] for s in substrates(reaction)]
    prod_names = [p[1] for p in products(reaction)]

    ss_isom = Int[]
    for (i, s) in enumerate(spec.steps)
        !s.is_equilibrium &&
            step_metabolite(s) === nothing &&
            push!(ss_isom, i)
    end

    result = @NamedTuple{
        r_only_mets::Vector{Symbol},
        r_only_cat_steps::Vector{Int}}[]

    # K-type: non-empty subsets of substrates ×
    # non-empty subsets of products, all r_only
    n_s = length(sub_names)
    n_p = length(prod_names)
    for s_mask in 1:(1 << n_s) - 1
        absent_subs = Symbol[sub_names[j]
            for j in 1:n_s
            if (s_mask >> (j - 1)) & 1 == 1]
        for p_mask in 1:(1 << n_p) - 1
            absent_prods = Symbol[prod_names[j]
                for j in 1:n_p
                if (p_mask >> (j - 1)) & 1 == 1]
            push!(result, (
                r_only_mets=Symbol[
                    absent_subs; absent_prods],
                r_only_cat_steps=Int[]))
        end
    end

    # V-type: all SS isomerization steps r_only
    if !isempty(ss_isom)
        push!(result, (r_only_mets=Symbol[],
            r_only_cat_steps=copy(ss_isom)))
    end

    result
end

"""
    _expand_to_allosteric(spec, reaction)
        → Vector{AllostericMechanismSpec}

Convert non-allosteric mechanism to allosteric (+1 param
for L). K-type: ≥1 substrate + ≥1 product absent from
T-state. V-type: all SS isomerization steps inactive in
T-state (kf_T=kr_T=0).
"""
function _expand_to_allosteric(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
    cn = oligomeric_state(reaction)
    sub_names = [s[1] for s in substrates(reaction)]
    prod_names = [p[1] for p in products(reaction)]
    all_cat = Symbol[sub_names; prod_names]

    ss_isom = Int[i for (i, s) in enumerate(spec.steps)
        if !s.is_equilibrium &&
           step_metabolite(s) === nothing]

    result = AllostericMechanismSpec[]

    for diff in _valid_allosteric_differentiations(
            reaction, spec)
        absent = Set(diff.r_only_mets)
        tr_equiv = Symbol[
            m for m in all_cat if m ∉ absent]
        tr_steps = Int[i for i in ss_isom
            if i ∉ diff.r_only_cat_steps]

        push!(result, AllostericMechanismSpec(
            spec, cn,
            Vector{Symbol}[], Int[],
            tr_equiv, tr_steps,
            diff.r_only_mets,
            Symbol[],  # no t_only metabolites
            diff.r_only_cat_steps,
            spec.param_count + 1))
    end
    result
end

function _expand_to_allosteric(
    ::AllostericMechanismSpec,
    @nospecialize(::EnzymeReaction))
    AllostericMechanismSpec[]
end

"""
    _expand_add_allosteric_regulator(spec, reaction)
        → Vector{AllostericMechanismSpec}

Add one allosteric regulator not yet in the mechanism.
Three flavors (r_only, t_only, tr_equiv) × site options
(new site or same site as each existing regulator site).
"""
function _expand_add_allosteric_regulator(
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    existing_allo = Set{Symbol}()
    for site in spec.allosteric_reg_sites
        for lig in site
            push!(existing_allo, lig)
        end
    end

    # Dead-end regulators already in base mechanism
    existing_de = Set{Symbol}()
    for s in spec.base.steps
        for sym in Iterators.flatten(
                (s.reactants, s.products))
            m = match(r"^(.+)__reg\d*$", string(sym))
            m !== nothing &&
                push!(existing_de, Symbol(m.captures[1]))
        end
    end

    roles = regulator_roles(reaction)
    new_regs = Symbol[]
    for (name, role) in roles
        (role == :unknown || role == :allosteric) ||
            continue
        name in existing_allo && continue
        name in existing_de && continue
        push!(new_regs, name)
    end
    sort!(new_regs)

    isempty(new_regs) && return AllostericMechanismSpec[]
    result = AllostericMechanismSpec[]

    for reg in new_regs
        n_sites = length(spec.allosteric_reg_sites)

        for mode in (:r_only, :t_only, :tr_equiv)
            # 0 = new site, 1..n_sites = existing site
            for site_idx in 0:n_sites
                new_sites = deepcopy(
                    spec.allosteric_reg_sites)
                new_mults = copy(
                    spec.allosteric_multiplicities)
                new_tr = copy(spec.tr_equiv_metabolites)
                new_r = copy(spec.r_only_metabolites)
                new_t = copy(spec.t_only_metabolites)

                if site_idx == 0
                    push!(new_sites, Symbol[reg])
                    push!(new_mults, spec.catalytic_n)
                else
                    push!(new_sites[site_idx], reg)
                end

                if mode == :tr_equiv
                    push!(new_tr, reg)
                elseif mode == :r_only
                    push!(new_r, reg)
                else
                    push!(new_t, reg)
                end

                push!(result, AllostericMechanismSpec(
                    spec.base, spec.catalytic_n,
                    new_sites, new_mults,
                    new_tr,
                    copy(spec.tr_equiv_cat_steps),
                    new_r, new_t,
                    copy(spec.r_only_cat_steps),
                    spec.param_count + 1))
            end
        end
    end
    result
end

function _expand_add_allosteric_regulator(
    ::MechanismSpec,
    @nospecialize(::EnzymeReaction),
)
    AllostericMechanismSpec[]
end

"""
    _tr_equiv_met_delta(met, steps, allosteric_reg_sites; param_constraints) → Int

Count how many new T-state independent params are added
when removing `met` from tr_equiv_metabolites.
RE binding steps add 1 (K_T), SS binding steps add 2
(kf_T and kr_T, both independent in T-state).
Allosteric regulators always add 1 (one K_T per reg site).
Constrained follower steps (targets of equivalence constraints)
are skipped since they share parameters with the leader step.
"""
function _tr_equiv_met_delta(
    met::Symbol, steps::Vector{StepSpec},
    allosteric_reg_sites::Vector{Vector{Symbol}}=Vector{Symbol}[];
    param_constraints::Vector{ParamConstraint}=ParamConstraint[])
    for site in allosteric_reg_sites
        met in site && return 1
    end
    follower_idxs = Set{Int}()
    for (target, _, _) in param_constraints
        idx = _step_index_from_constraint_sym(target)
        idx !== nothing && push!(follower_idxs, idx)
    end
    delta = 0
    for (idx, s) in enumerate(steps)
        step_metabolite(s) === met || continue
        idx in follower_idxs && continue
        delta += s.is_equilibrium ? 1 : 2
    end
    delta
end

"""
    _expand_remove_tr_equiv(spec, reaction)
        → Vector{AllostericMechanismSpec}

Remove one TR equivalence (metabolite or catalytic step),
making T-state and R-state parameters independent.
RE metabolite removal adds +1; SS metabolite removal adds +2
(both kf_T and kr_T become independent); catalytic step removal
adds +1.
"""
function _expand_remove_tr_equiv(
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    result = AllostericMechanismSpec[]

    for (i, met) in enumerate(spec.tr_equiv_metabolites)
        new_equiv = [spec.tr_equiv_metabolites[j]
            for j in eachindex(spec.tr_equiv_metabolites)
            if j != i]
        delta = _tr_equiv_met_delta(
            met, spec.base.steps,
            spec.allosteric_reg_sites;
            param_constraints=spec.base.param_constraints)
        push!(result, AllostericMechanismSpec(
            spec.base, spec.catalytic_n,
            deepcopy(spec.allosteric_reg_sites),
            copy(spec.allosteric_multiplicities),
            new_equiv,
            copy(spec.tr_equiv_cat_steps),
            copy(spec.r_only_metabolites),
            copy(spec.t_only_metabolites),
            copy(spec.r_only_cat_steps),
            spec.param_count + delta))
    end

    for (i, _) in enumerate(spec.tr_equiv_cat_steps)
        new_steps = [spec.tr_equiv_cat_steps[j]
            for j in eachindex(spec.tr_equiv_cat_steps)
            if j != i]
        push!(result, AllostericMechanismSpec(
            spec.base, spec.catalytic_n,
            deepcopy(spec.allosteric_reg_sites),
            copy(spec.allosteric_multiplicities),
            copy(spec.tr_equiv_metabolites),
            new_steps,
            copy(spec.r_only_metabolites),
            copy(spec.t_only_metabolites),
            copy(spec.r_only_cat_steps),
            spec.param_count + 1))
    end

    # Remove one r_only cat step → step becomes independent (+1)
    # Only when no metabolites are r_only/t_only (otherwise
    # kf_T is unidentifiable — state can't catalyze)
    if isempty(spec.r_only_metabolites) &&
            isempty(spec.t_only_metabolites)
        for (i, _) in enumerate(spec.r_only_cat_steps)
            new_r_steps = [spec.r_only_cat_steps[j]
                for j in eachindex(spec.r_only_cat_steps)
                if j != i]
            push!(result, AllostericMechanismSpec(
                spec.base, spec.catalytic_n,
                deepcopy(spec.allosteric_reg_sites),
                copy(spec.allosteric_multiplicities),
                copy(spec.tr_equiv_metabolites),
                copy(spec.tr_equiv_cat_steps),
                copy(spec.r_only_metabolites),
                copy(spec.t_only_metabolites),
                new_r_steps,
                spec.param_count + 1))
        end
    end
    result
end

function _expand_remove_tr_equiv(
    ::MechanismSpec,
    @nospecialize(::EnzymeReaction),
)
    AllostericMechanismSpec[]
end

"""Construct AllostericEnzymeMechanism from AllostericMechanismSpec."""
function AllostericEnzymeMechanism(spec::AllostericMechanismSpec)
    cm = EnzymeMechanism(spec.base)
    cat_mets = metabolites(cm)

    # Build Metabolites tuple (catalytic + regulatory)
    reg_syms = Symbol[]
    for site in spec.allosteric_reg_sites
        for s in site
            s in reg_syms || s in cat_mets ||
                push!(reg_syms, s)
        end
    end
    mets = (cat_mets..., reg_syms...)

    # Build CatSites: (catalytic_metabolites, multiplicity,
    #   tr_equiv_mets, tr_equiv_cat_steps,
    #   r_only_mets, t_only_mets, r_only_cat_steps)
    cat_tr = Tuple(m for m in cat_mets
                   if m in spec.tr_equiv_metabolites)
    cat_steps_tr = Tuple(spec.tr_equiv_cat_steps)
    cat_r_only = Tuple(m for m in cat_mets
                       if m in spec.r_only_metabolites)
    cat_t_only = Tuple(m for m in cat_mets
                       if m in spec.t_only_metabolites)
    cat_r_only_steps = Tuple(spec.r_only_cat_steps)
    cat_sites = (cat_mets, spec.catalytic_n, cat_tr,
                 cat_steps_tr, cat_r_only, cat_t_only,
                 cat_r_only_steps)

    # Build RegSites with TR equivalence and
    # r_only/t_only info
    reg_sites = Tuple(
        (Tuple(group), mult,
         Tuple(lig for lig in group
               if lig in spec.tr_equiv_metabolites),
         Tuple(lig for lig in group
               if lig in spec.r_only_metabolites),
         Tuple(lig for lig in group
               if lig in spec.t_only_metabolites))
        for (group, mult) in zip(
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities))

    AllostericEnzymeMechanism{
        mets, typeof(cm), cat_sites, reg_sites}()
end

"""
    _rewrap_allosteric(original, new_base) → AllostericMechanismSpec

Replace the base of an allosteric spec with a new base,
preserving all allosteric structure.
"""
function _rewrap_allosteric(
    original::AllostericMechanismSpec,
    new_base::MechanismSpec)
    AllostericMechanismSpec(
        new_base, original.catalytic_n,
        deepcopy(original.allosteric_reg_sites),
        copy(original.allosteric_multiplicities),
        copy(original.tr_equiv_metabolites),
        copy(original.tr_equiv_cat_steps),
        copy(original.r_only_metabolites),
        copy(original.t_only_metabolites),
        copy(original.r_only_cat_steps),
        original.param_count +
            new_base.param_count -
            original.base.param_count)
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

function _step_sort_key(s::StepSpec)
    (sort(s.reactants), sort(s.products),
     s.is_equilibrium)
end

"""Remap a constraint symbol given a step index remapping (old index → new index)."""
function _remap_constraint_sym(sym::Symbol, idx_map::Dict{Int,Int})
    s = string(sym)
    m = match(r"^([Kk])(\d+)(.*)", s)
    m === nothing && return sym
    cap_prefix = m.captures[1]
    cap_idx = m.captures[2]
    cap_suffix = m.captures[3]
    (cap_prefix === nothing || cap_idx === nothing ||
        cap_suffix === nothing) && return sym
    old_idx = parse(Int, cap_idx)
    haskey(idx_map, old_idx) || return sym
    Symbol(cap_prefix, idx_map[old_idx], cap_suffix)
end

function _canonicalize!(spec::MechanismSpec)
    # Compute permutation before sorting
    n = length(spec.steps)
    perm = sortperm(spec.steps, by=_step_sort_key)
    # Build old→new index map
    idx_map = Dict{Int,Int}(perm[new] => new for new in 1:n)
    sort!(spec.steps, by=_step_sort_key)
    # Update constraint symbols to reflect new step positions
    for i in eachindex(spec.param_constraints)
        (target, coeff, followers) = spec.param_constraints[i]
        new_target = _remap_constraint_sym(target, idx_map)
        new_followers = [
            (_remap_constraint_sym(s, idx_map), c)
            for (s, c) in followers]
        spec.param_constraints[i] = (new_target, coeff, new_followers)
    end
    sort!(spec.param_constraints, by=c -> c[1])
    idx_map
end

function _canonicalize!(spec::AllostericMechanismSpec)
    idx_map = _canonicalize!(spec.base)
    # Remap step indices that refer to base.steps positions
    map!(i -> get(idx_map, i, i), spec.tr_equiv_cat_steps,
        spec.tr_equiv_cat_steps)
    map!(i -> get(idx_map, i, i), spec.r_only_cat_steps,
        spec.r_only_cat_steps)
    sort!(spec.tr_equiv_metabolites)
    sort!(spec.tr_equiv_cat_steps)
    sort!(spec.r_only_metabolites)
    sort!(spec.t_only_metabolites)
    sort!(spec.r_only_cat_steps)
    for site in spec.allosteric_reg_sites
        sort!(site)
    end
    # Sort sites themselves (with multiplicities) by content
    if length(spec.allosteric_reg_sites) >= 2
        perm = sortperm(spec.allosteric_reg_sites)
        spec.allosteric_reg_sites .= spec.allosteric_reg_sites[perm]
        spec.allosteric_multiplicities .= spec.allosteric_multiplicities[perm]
    end
    spec
end

function _dedup_key(spec::MechanismSpec)
    steps = Tuple(
        (Tuple(sort(s.reactants)),
         Tuple(sort(s.products)),
         s.is_equilibrium)
        for s in spec.steps)
    constraints = Tuple(
        (c[1], c[2], Tuple(c[3]))
        for c in spec.param_constraints)
    (steps, constraints)
end

function _dedup_key(spec::AllostericMechanismSpec)
    base_key = _dedup_key(spec.base)
    (base_key, spec.catalytic_n,
     Tuple(Tuple.(spec.allosteric_reg_sites)),
     Tuple(spec.allosteric_multiplicities),
     Tuple(spec.tr_equiv_metabolites),
     Tuple(spec.tr_equiv_cat_steps),
     Tuple(spec.r_only_metabolites),
     Tuple(spec.t_only_metabolites),
     Tuple(spec.r_only_cat_steps))
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
