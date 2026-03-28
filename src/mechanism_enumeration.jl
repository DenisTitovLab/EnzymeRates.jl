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
