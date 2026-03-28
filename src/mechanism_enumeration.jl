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
