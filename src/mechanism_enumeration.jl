# ABOUTME: Mechanism enumeration by incremental parameter count growth
# ABOUTME: Provides init_mechanisms, expand_mechanisms, dedup! building blocks

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
