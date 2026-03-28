# ABOUTME: Mechanism enumeration by incremental parameter count growth
# ABOUTME: Provides init_mechanisms, expand_mechanisms, dedup! building blocks

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
    expected_pc = n_s + n_p + 3

    result = MechanismSpec[]
    for spec in expanded
        constraints = _max_equivalence_constraints(spec)
        push!(result, MechanismSpec(
            spec.reaction, spec.steps,
            constraints, expected_pc))
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

        # Propagate SS to dead-end mirror steps
        from_form, to_form = step_forms(s)
        for (j, ms) in enumerate(new_steps)
            j == i && continue
            ms.is_equilibrium || continue
            mf, mt = step_forms(ms)
            if _is_mirror_of(mf, mt, from_form, to_form, spec.steps)
                new_steps[j] = StepSpec(ms.reactants, ms.products, false)
            end
        end

        push!(result, MechanismSpec(
            spec.reaction, new_steps,
            copy(spec.param_constraints),
            spec.param_count + 1))
    end
    result
end

"""
    _expand_remove_constraint(spec::MechanismSpec) → Vector{MechanismSpec}

Remove one equivalence constraint (+1 estimated param).
Each removable constraint produces one new mechanism.
"""
function _expand_remove_constraint(spec::MechanismSpec)
    result = MechanismSpec[]
    for i in eachindex(spec.param_constraints)
        new_constraints = [
            spec.param_constraints[j]
            for j in eachindex(spec.param_constraints)
            if j != i]
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

    for (ri, reg) in enumerate(eligible_regs)
        dummy = Symbol(
            string(reg) * "__reg" * string(ri))

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
        n_forms = length(eligible_forms)

        # Enumerate all non-empty subsets
        for mask in 1:(1 << n_forms) - 1
            active = Symbol[]
            for (j, f) in enumerate(eligible_forms)
                if (mask >> (j - 1)) & 1 == 1
                    push!(active, f)
                end
            end

            new_steps = copy(spec.steps)
            de_form_map = Dict{Symbol, Symbol}()

            # Add binding steps (always RE)
            binding_step_indices = Int[]
            for cf in active
                de_name = _dead_end_form_name(
                    bound[cf], dummy)
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
            new_constraints = copy(
                spec.param_constraints)
            if length(binding_step_indices) >= 2
                first_idx = binding_step_indices[1]
                for j in 2:length(binding_step_indices)
                    push!(new_constraints, (
                        Symbol("K$(binding_step_indices[j])"),
                        1,
                        [(Symbol("K$(first_idx)"), 1)]))
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
