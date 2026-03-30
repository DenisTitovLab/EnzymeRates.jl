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
    _tr_equiv_met_delta(met, steps, allosteric_reg_sites) → Int

Count how many new T-state independent params are added
when removing `met` from tr_equiv_metabolites.
RE binding steps add 1 (K_T), SS binding steps add 2
(kf_T and kr_T, both independent in T-state).
Allosteric regulators always add 1 (one K_T per reg site).
"""
function _tr_equiv_met_delta(
    met::Symbol, steps::Vector{StepSpec},
    allosteric_reg_sites::Vector{Vector{Symbol}}=Vector{Symbol}[])
    for site in allosteric_reg_sites
        met in site && return 1
    end
    delta = 0
    for s in steps
        step_metabolite(s) === met || continue
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
            spec.allosteric_reg_sites)
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
    spec::MechanismSpec,
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
end

function _add_expansions!(
    result::Dict{Int, Vector{AbstractMechanismSpec}},
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction))
    for s in _expand_add_allosteric_regulator(
            spec, reaction)
        _push_to_dict!(result, s)
    end
    for s in _expand_remove_tr_equiv(spec, reaction)
        _push_to_dict!(result, s)
    end
    # Apply base moves, rewrap
    for new_base in _expand_re_to_ss(spec.base)
        _push_to_dict!(result,
            _rewrap_allosteric(spec, new_base))
    end
    for new_base in _expand_remove_constraint(spec.base)
        _push_to_dict!(result,
            _rewrap_allosteric(spec, new_base))
    end
    # Dead-end regs: exclude allosteric regs
    allo_regs = Set{Symbol}()
    for site in spec.allosteric_reg_sites
        for lig in site
            push!(allo_regs, lig)
        end
    end
    for new_base in _expand_add_dead_end_regulator(
            spec.base, reaction;
            exclude_regs=allo_regs)
        _push_to_dict!(result,
            _rewrap_allosteric(spec, new_base))
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
