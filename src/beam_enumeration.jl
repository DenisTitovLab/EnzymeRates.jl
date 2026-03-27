# ABOUTME: Beam-based mechanism enumeration with runtime parameter counting.
# ABOUTME: Avoids JIT compilation by using runtime versions of derivation functions.

# ─── Data Extraction from MechanismSpec ──────────────────────

"""Extract (enz_names, enz_set) from MechanismSpec steps."""
function _spec_enzyme_names(spec::MechanismSpec)
    forms = all_form_names(spec)
    enz_names = Tuple(sort!(collect(forms)))
    enz_names, Set(enz_names)
end

"""
Extract reactions as Vector of (lhs_tuple, rhs_tuple) from
MechanismSpec steps, matching the format of `reactions(m)`.
"""
function _spec_reactions(spec::MechanismSpec)
    [(Tuple(s.reactants), Tuple(s.products)) for s in spec.steps]
end

"""Extract equilibrium_steps as a Tuple of Bool from MechanismSpec."""
function _spec_eq_steps(spec::MechanismSpec)
    Tuple(s.is_equilibrium for s in spec.steps)
end

"""
Compute stoichiometric matrix from MechanismSpec steps.
Rows = metabolites, columns = steps.
"""
function _spec_stoich_matrix(
    spec::MechanismSpec, met_names, enz_set,
)
    S = zeros(Int, length(met_names), length(spec.steps))
    met_idx = Dict(m => i for (i, m) in enumerate(met_names))
    for (j, step) in enumerate(spec.steps)
        for s in step.reactants
            s in enz_set && continue
            S[met_idx[s], j] -= 1
        end
        for s in step.products
            s in enz_set && continue
            S[met_idx[s], j] += 1
        end
    end
    S
end

"""
Collect metabolite names from the reaction definition.
Order: substrates then products then regulators, unique.
"""
function _spec_met_names(rxn)
    seen = Set{Symbol}()
    mets = Symbol[]
    for (name, _) in substrates(rxn)
        if name ∉ seen
            push!(seen, name)
            push!(mets, name)
        end
    end
    for (name, _) in products(rxn)
        if name ∉ seen
            push!(seen, name)
            push!(mets, name)
        end
    end
    for name in regulators(rxn)
        if name ∉ seen
            push!(seen, name)
            push!(mets, name)
        end
    end
    mets
end

"""
Determine free enzyme forms from MechanismSpec steps.
A form is free if no binding step produces it (i.e., it never
appears as `products[1]` of a step with a metabolite reactant).
"""
function _spec_free_enz_set(spec::MechanismSpec, enz_set)
    bound_forms = Set{Symbol}()
    for s in spec.steps
        if step_metabolite(s) !== nothing
            push!(bound_forms, s.products[1])
        end
    end
    Set(f for f in enz_set if f ∉ bound_forms)
end

# ─── Runtime Parameter Counting ─────────────────────────────

"""
    _runtime_param_count(spec::MechanismSpec)

Count independent parameters for a mechanism without JIT
compilation. Uses runtime versions of thermodynamic constraint
and dependent parameter analysis.

Returns the total parameter count including Keq and E_total.
"""
function _runtime_param_count(spec::MechanismSpec)
    steps = spec.steps
    n_steps = length(steps)
    n_re = count(s -> s.is_equilibrium, steps)
    n_ss = n_steps - n_re
    n_forms = length(all_form_names(steps))
    n_thermo = n_steps - n_forms + 1
    n_redundant = _constraint_thermo_redundancy(spec)
    n_re + 2 * n_ss - (n_thermo - n_redundant) +
        2 - length(spec.param_constraints)
end

"""
    _constraint_thermo_redundancy(spec::MechanismSpec) → Int

Count thermodynamic constraints made redundant by parameter
equivalence constraints. When constrained steps (mirror steps)
form cycles in the enzyme form graph, each independent cycle
has a Wegscheider condition that is trivially satisfied by the
equivalence constraints, reducing the effective number of
thermodynamic constraints.

Returns the number of independent cycles in the subgraph
formed by constrained steps only.
"""
function _constraint_thermo_redundancy(spec::MechanismSpec)
    constraints = spec.param_constraints
    isempty(constraints) && return 0

    # Collect step indices that appear in constraints
    constrained_indices = Set{Int}()
    for (dep, _, refs) in constraints
        m = match(r"^[kK](\d+)", string(dep))
        if m !== nothing
            cap = m[1]::SubString
            push!(constrained_indices, parse(Int, cap))
        end
        for (ref, _) in refs
            m = match(r"^[kK](\d+)", string(ref))
            if m !== nothing
                cap = m[1]::SubString
                push!(constrained_indices, parse(Int, cap))
            end
        end
    end

    isempty(constrained_indices) && return 0

    # Build subgraph of constrained steps only. Count
    # forms (nodes) and edges in this subgraph.
    forms = Set{Symbol}()
    n_edges = 0
    for idx in constrained_indices
        idx > length(spec.steps) && continue
        s = spec.steps[idx]
        push!(forms, s.reactants[1])
        push!(forms, s.products[1])
        n_edges += 1
    end

    # Independent cycles in the constrained subgraph.
    # For a connected graph: cycles = edges - nodes + 1.
    # For multiple components: cycles = edges - nodes + components.
    n_nodes = length(forms)
    n_components = _count_components(
        forms, spec.steps, constrained_indices)
    max(0, n_edges - n_nodes + n_components)
end

"""Count connected components in the subgraph of constrained steps."""
function _count_components(
    forms::Set{Symbol},
    steps::Vector{StepSpec},
    constrained_indices::Set{Int},
)
    parent = Dict(f => f for f in forms)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    function union!(a, b)
        ra, rb = find(a), find(b)
        ra != rb && (parent[ra] = rb)
    end

    for idx in constrained_indices
        idx > length(steps) && continue
        s = steps[idx]
        union!(s.reactants[1], s.products[1])
    end

    length(Set(find(f) for f in forms))
end

# ─── Runtime Parameter Counting: AllostericMechanismSpec ───

"""
    _runtime_param_count(spec::AllostericMechanismSpec)

Count independent parameters for an allosteric mechanism
without JIT compilation. Mirrors the logic of
`_dependent_param_exprs(::Type{AllostericEnzymeMechanism})`.
"""
function _runtime_param_count(spec::AllostericMechanismSpec)
    base = spec.base
    rxn = base.reaction

    # Get base indep params (without Keq, E_total)
    enz_names, enz_set = _spec_enzyme_names(base)
    rxns = _spec_reactions(base)
    eq_steps = _spec_eq_steps(base)
    met_names = _spec_met_names(rxn)
    stoich_mat = _spec_stoich_matrix(base, met_names, enz_set)
    free_enz_set = _spec_free_enz_set(base, enz_set)
    binding_Ks = Set(_binding_K_symbols(rxns, eq_steps, enz_set))

    C, rhs_coeffs = _thermodynamic_constraints(
        enz_names, enz_set, rxns, stoich_mat,
        met_names, substrates(rxn), products(rxn),
    )
    _, indep_R = _dependent_param_exprs(
        C, rhs_coeffs, eq_steps, base.param_constraints,
        binding_Ks, rxns, enz_set, free_enz_set,
    )

    # Count TR-equivalent catalytic params among indep_R
    # (r_only/t_only catalytic params don't affect param count —
    # they only change rate equation structure, not parameter set)
    n_tr_equiv = 0
    for p in indep_R
        mode = _classify_catalytic_param(
            p, base, eq_steps, spec)
        if mode == :tr_equiv
            n_tr_equiv += 1
        end
    end

    # T-state indep = base indep minus TR-equiv params
    indep_T_count = length(indep_R) - n_tr_equiv

    # Reg R-state params: exclude t_only ligands (no R-state binding)
    n_reg_R = sum(
        count(lig -> lig ∉ spec.t_only_metabolites, site)
        for site in spec.allosteric_reg_sites; init=0)

    # Reg T-state params: exclude tr_equiv, r_only, and t_only
    n_reg_T = sum(
        count(lig -> lig ∉ spec.tr_equiv_metabolites &&
                     lig ∉ spec.r_only_metabolites &&
                     lig ∉ spec.t_only_metabolites, site)
        for site in spec.allosteric_reg_sites; init=0)

    # t_only reg params: T-state is independent
    n_reg_t_only = sum(
        count(lig -> lig ∈ spec.t_only_metabolites, site)
        for site in spec.allosteric_reg_sites; init=0)

    # Total: base_indep + indep_T + reg_R + reg_T + reg_t_only + L
    #        + Keq + E_total
    length(indep_R) + indep_T_count + n_reg_R + n_reg_T +
        n_reg_t_only + 1 + 2
end

"""Classify a catalytic parameter as :both, :tr_equiv, :r_only, or :t_only."""
function _classify_catalytic_param(
    p::Symbol, base::MechanismSpec,
    eq_steps, spec::AllostericMechanismSpec,
)
    m = match(r"^K(\d+)$", string(p))
    if m !== nothing
        cap = m.captures[1]::SubString
        idx = parse(Int, cap)
        if idx <= length(base.steps) && eq_steps[idx]
            met = step_metabolite(base.steps[idx])
            if met !== nothing
                met in spec.tr_equiv_metabolites && return :tr_equiv
                met in spec.r_only_metabolites && return :r_only
                met in spec.t_only_metabolites && return :t_only
            end
        end
        return :both
    end
    if _is_ss_rate_constant(p)
        km = match(r"^k(\d+)[fr]$", string(p))
        km === nothing && return :both
        cap = km.captures[1]::SubString
        idx = parse(Int, cap)
        met = step_metabolite(base.steps[idx])
        if met !== nothing
            met in spec.tr_equiv_metabolites && return :tr_equiv
            met in spec.r_only_metabolites && return :r_only
            met in spec.t_only_metabolites && return :t_only
        else
            idx in spec.tr_equiv_cat_steps && return :tr_equiv
            idx in spec.r_only_cat_steps && return :r_only
        end
    end
    :both
end

# ─── Kinetic Symbol Detection ───────────────────────────────

"""
    _is_kinetic_symbol(s::Symbol)

True for rate/equilibrium constant symbols like K1, k2f, k3r.
Pattern: `^[kK]\\d+[fr]?\$`.
"""
function _is_kinetic_symbol(s::Symbol)
    str = string(s)
    m = match(r"^[kK]\d+[fr]?$", str)
    m !== nothing
end

# ─── Runtime Denominator Fingerprinting ──────────────────────

"""
    _strip_to_concentration_fingerprint(denom_terms)

Expand denominator terms to a flat polynomial, then strip all
kinetic symbols from each monomial, returning a `Set{MONO}` of
concentration-only monomials.
"""
function _strip_to_concentration_fingerprint(
    denom_terms::Vector{DenomTerm},
)
    flat = _expand_to_poly(denom_terms)
    result = Set{MONO}()
    for (mono, _) in flat
        conc_mono = Pair{Symbol,Int}[
            p for p in mono if !_is_kinetic_symbol(p.first)
        ]
        push!(result, conc_mono)
    end
    result
end

"""
    expand_mechanisms_by_one_param(specs, reaction) → Vector{MechanismSpec}

Generate mechanism candidates with param_count + 1 via three moves:
RE→SS flip, remove equivalence constraint, or add dead-end binding.
"""
function expand_mechanisms_by_one_param(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = MechanismSpec[]
    for spec in specs
        _expand_re_to_ss!(result, spec)
        _expand_remove_constraint!(result, spec)
        _expand_add_dead_end!(result, spec, reaction)
    end
    result
end

"""
    expand_mechanisms_by_one_param(specs, reaction)
        → Vector{AllostericMechanismSpec}

Generate allosteric mechanism candidates with param_count + 1
via four moves: remove TR equivalence, RE→SS on base,
remove constraint on base, or add dead-end binding on base.
"""
function expand_mechanisms_by_one_param(
    specs::Vector{AllostericMechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    result = AllostericMechanismSpec[]
    for spec in specs
        _expand_remove_tr_equiv!(result, spec)
        _expand_base_moves!(result, spec, reaction)
    end
    result
end

"""Remove one TR equivalence: make one metabolite's K_T ≠ K_R,
or one catalytic step's kf_T ≠ kf_R."""
function _expand_remove_tr_equiv!(
    result::Vector{AllostericMechanismSpec},
    spec::AllostericMechanismSpec,
)
    # Remove one metabolite TR equivalence
    for (i, met) in enumerate(spec.tr_equiv_metabolites)
        new_equiv = [
            spec.tr_equiv_metabolites[j]
            for j in eachindex(spec.tr_equiv_metabolites)
            if j != i]
        push!(result, AllostericMechanismSpec(
            spec.base, spec.catalytic_n,
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities,
            new_equiv,
            copy(spec.tr_equiv_cat_steps),
            copy(spec.r_only_metabolites),
            copy(spec.t_only_metabolites),
            copy(spec.r_only_cat_steps)))
    end
    # Remove one catalytic step TR equivalence
    for (i, idx) in enumerate(spec.tr_equiv_cat_steps)
        new_cat_steps = [
            spec.tr_equiv_cat_steps[j]
            for j in eachindex(spec.tr_equiv_cat_steps)
            if j != i]
        push!(result, AllostericMechanismSpec(
            spec.base, spec.catalytic_n,
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities,
            copy(spec.tr_equiv_metabolites),
            new_cat_steps,
            copy(spec.r_only_metabolites),
            copy(spec.t_only_metabolites),
            copy(spec.r_only_cat_steps)))
    end
end

"""Apply base mechanism +1 moves (RE→SS, remove constraint,
dead-end binding) to produce AllostericMechanismSpec variants
with modified bases."""
function _expand_base_moves!(
    result::Vector{AllostericMechanismSpec},
    spec::AllostericMechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    # Get +1 expansions of the base MechanismSpec
    base_results = MechanismSpec[]
    _expand_re_to_ss!(base_results, spec.base)
    _expand_remove_constraint!(base_results, spec.base)
    _expand_add_dead_end!(base_results, spec.base, reaction)

    # Wrap each expanded base in the allosteric structure
    for new_base in base_results
        push!(result, AllostericMechanismSpec(
            new_base, spec.catalytic_n,
            spec.allosteric_reg_sites,
            spec.allosteric_multiplicities,
            copy(spec.tr_equiv_metabolites),
            copy(spec.tr_equiv_cat_steps),
            copy(spec.r_only_metabolites),
            copy(spec.t_only_metabolites),
            copy(spec.r_only_cat_steps)))
    end
end

"""
Collect all parameter symbols referenced in constraints
(both the constrained symbol and the reference symbols).
"""
function _constrained_symbols(constraints::Vector{ParamConstraint})
    syms = Set{Symbol}()
    for (dep, _, refs) in constraints
        push!(syms, dep)
        for (ref, _) in refs
            push!(syms, ref)
        end
    end
    syms
end

"""Convert each RE step to SS, producing one candidate per RE step.
Skips steps whose K parameter appears in any constraint, as flipping
them would invalidate the constraint.
"""
function _expand_re_to_ss!(result::Vector{MechanismSpec}, spec::MechanismSpec)
    constrained = _constrained_symbols(spec.param_constraints)
    for (i, step) in enumerate(spec.steps)
        step.is_equilibrium || continue
        Symbol("K$i") in constrained && continue
        new_steps = copy(spec.steps)
        new_steps[i] = StepSpec(step.reactants, step.products, false)
        candidate = MechanismSpec(
            spec.reaction, new_steps,
            copy(spec.param_constraints), 0)
        push!(result, MechanismSpec(
            spec.reaction, new_steps,
            copy(spec.param_constraints),
            _runtime_param_count(candidate)))
    end
end

"""Remove one equivalence constraint, producing one candidate per constraint."""
function _expand_remove_constraint!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
)
    for i in eachindex(spec.param_constraints)
        new_constraints = [
            spec.param_constraints[j]
            for j in eachindex(spec.param_constraints)
            if j != i]
        candidate = MechanismSpec(
            spec.reaction, copy(spec.steps),
            new_constraints, 0)
        push!(result, MechanismSpec(
            spec.reaction, copy(spec.steps),
            new_constraints,
            _runtime_param_count(candidate)))
    end
end

# ─── Dead-End Binding Helpers ─────────────────────────────────

"""
    _bound_entities(spec) → Dict{Symbol, Set{Symbol}}

For each enzyme form in the mechanism, return the set of
metabolites/regulators bound to it. Derived from step structure:
if a step binds metabolite M to form F producing form FM, then
FM has M plus everything F had.
"""
function _bound_entities(spec::MechanismSpec)
    bound = Dict{Symbol, Set{Symbol}}()
    for s in spec.steps
        for f in (s.reactants[1], s.products[1])
            haskey(bound, f) || (bound[f] = Set{Symbol}())
        end
    end

    # Propagate only through binding steps (steps with a
    # metabolite). A binding step [F, M] → [FM] means FM
    # has everything F has plus M. An unbinding step goes
    # in reverse. Isomerization steps (no metabolite) don't
    # change bound state — they transform bound molecules.
    changed = true
    while changed
        changed = false
        for s in spec.steps
            met = step_metabolite(s)
            met === nothing && continue
            from_form = s.reactants[1]
            to_form = s.products[1]

            # Forward: to_form = from_form ∪ {met}
            expected_to = union(bound[from_form], Set([met]))
            if !issubset(expected_to, bound[to_form])
                union!(bound[to_form], expected_to)
                changed = true
            end

            # Reverse: from_form = to_form \ {met}
            expected_from = setdiff(bound[to_form], Set([met]))
            if !issubset(expected_from, bound[from_form])
                union!(bound[from_form], expected_from)
                changed = true
            end
        end
    end
    bound
end

"""
    _binding_capacity(reaction) → Int

Maximum number of entities that can bind at the catalytic site.
Equal to max(n_substrates, n_products).
"""
function _binding_capacity(
    @nospecialize(reaction::EnzymeReaction),
)
    max(length(substrates(reaction)),
        length(products(reaction)))
end

# ─── Dead-End Binding Move ────────────────────────────────────

"""
    _beam_dead_end_form_name(bound, new_met) → Symbol

Generate a canonical form name for a dead-end complex.
Sorts all bound metabolites alphabetically.
"""
function _beam_dead_end_form_name(
    bound::Set{Symbol}, new_met::Symbol,
)
    all_bound = sort!(collect(union(bound, Set([new_met]))))
    Symbol("E_" * join(all_bound, "_"))
end

"""
    _catalytic_from_forms(spec) → Dict{Symbol, Symbol}

For each metabolite in the mechanism's catalytic steps, return
the enzyme form it binds to. E.g., if step is [EA, B] ⇌ [EAB],
then B's catalytic from-form is EA.
"""
function _catalytic_from_forms(spec::MechanismSpec)
    result = Dict{Symbol, Symbol}()
    for s in spec.steps
        met = step_metabolite(s)
        met === nothing && continue
        haskey(result, met) && continue
        result[met] = s.reactants[1]
    end
    result
end

"""
    _dead_end_opportunities(spec, reaction;
        catalytic_forms=nothing)
        → Vector{Tuple{Symbol, Symbol}}

Find (form, metabolite) pairs where dead-end binding is possible.
A binding is valid when:
- The form is a catalytic form (not a dead-end form from prior expansion)
- The metabolite is not already bound at the form
- The result doesn't exceed binding capacity
- The result is not already a form in the mechanism
- For substrates/products: the result must have at least one substrate
  AND at least one product bound, and must not have ALL substrates
  or ALL products bound
- For substrates/products: the metabolite's catalytic binding step
  must start from a non-free enzyme form (metabolites that bind
  free enzyme directly don't create dead-ends)

When `catalytic_forms` is provided, only those forms are eligible
as binding targets. Otherwise all forms in the spec are eligible.
"""
function _dead_end_opportunities(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction);
    catalytic_forms::Union{Nothing, Set{Symbol}}=nothing,
)
    bound = _bound_entities(spec)
    cap = _binding_capacity(reaction)
    existing_forms = all_form_names(spec)
    eligible_forms = catalytic_forms !== nothing ?
        catalytic_forms : existing_forms

    sub_names = Set(s[1] for s in substrates(reaction))
    prod_names = Set(p[1] for p in products(reaction))
    all_mets = union(sub_names, prod_names)

    # Only metabolites whose catalytic from-form is NOT the
    # free enzyme can create dead-end complexes
    cat_from = _catalytic_from_forms(spec)
    free_enz = Set(f for f in keys(bound)
        if isempty(bound[f]))
    eligible_mets = Set{Symbol}(
        m for m in all_mets
        if haskey(cat_from, m) &&
            cat_from[m] ∉ free_enz)

    roles = regulator_roles(reaction)
    de_regs = Symbol[
        r[1] for r in roles
        if r[2] == :dead_end || r[2] == :unknown]

    opportunities = Tuple{Symbol, Symbol}[]

    for form in sort!(collect(keys(bound)))
        form in eligible_forms || continue
        fb = bound[form]
        n_bound = length(fb)
        n_bound >= cap && continue

        # Substrate/product dead-end opportunities
        fb_subs = intersect(fb, sub_names)
        fb_prods = intersect(fb, prod_names)
        # Skip forms that already have all subs or all prods
        (fb_subs == sub_names ||
            fb_prods == prod_names) && continue

        for met in sort!(collect(eligible_mets))
            met in fb && continue
            de_name = _beam_dead_end_form_name(fb, met)
            de_name in existing_forms && continue

            new_bound = union(fb, Set([met]))
            new_subs = intersect(new_bound, sub_names)
            new_prods = intersect(new_bound, prod_names)
            # Must have at least one sub AND one prod
            (isempty(new_subs) || isempty(new_prods)) &&
                continue
            # Must not have ALL subs or ALL prods
            (new_subs == sub_names ||
                new_prods == prod_names) && continue

            push!(opportunities, (form, met))
        end

        # Dead-end regulator opportunities
        for reg in de_regs
            reg in fb && continue
            de_name = _beam_dead_end_form_name(fb, reg)
            de_name in existing_forms && continue
            push!(opportunities, (form, reg))
        end
    end
    opportunities
end

"""
    _expand_add_dead_end!(result, spec, reaction)

Add dead-end binding configurations at exactly +1 param.
"""
function _expand_add_dead_end!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    _expand_dead_end_at_delta!(result, spec, reaction, 1)
end

"""
    _expand_dead_end_at_delta!(result, spec, reaction, delta;
        catalytic_forms=nothing)

Shared dead-end expansion logic. Groups opportunities by
dead-end form name and enumerates subsets of dead-end forms
across all metabolites. For each subset, adds binding steps
with maximal equivalence constraints. Keeps only configurations
where computed param_count equals spec.param_count + delta.

Dead-end steps inherit RE/SS status from their catalytic
counterpart. For delta=0, equivalence constraints tie dead-end
parameters to catalytic parameters (K for RE, kf+kr for SS).

When `catalytic_forms` is provided, only those forms are eligible
as binding targets for dead-end metabolites.
"""
function _expand_dead_end_at_delta!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
    delta::Int;
    catalytic_forms::Union{Nothing, Set{Symbol}}=nothing,
)
    opps = _dead_end_opportunities(
        spec, reaction; catalytic_forms)
    isempty(opps) && return

    bound = _bound_entities(spec)

    # Group opportunities by dead-end form name
    de_forms = Dict{Symbol,
        Vector{Tuple{Symbol, Symbol}}}()
    for (form, met) in opps
        de_name = _beam_dead_end_form_name(
            bound[form], met)
        push!(get!(de_forms, de_name,
            Tuple{Symbol, Symbol}[]), (form, met))
    end
    de_form_names = sort!(collect(keys(de_forms)))
    n_de = length(de_form_names)

    # Build catalytic step lookup: met → (step_index, is_RE)
    cat_step_info = Dict{Symbol, Tuple{Int, Bool}}()
    for (i, s) in enumerate(spec.steps)
        met = step_metabolite(s)
        met === nothing && continue
        haskey(cat_step_info, met) && continue
        cat_step_info[met] = (i, s.is_equilibrium)
    end

    target_pc = spec.param_count + delta

    # Enumerate subsets of dead-end forms
    for mask in 1:(1 << n_de) - 1
        active_de = Symbol[
            de_form_names[j]
            for j in 1:n_de
            if (mask >> (j - 1)) & 1 == 1]

        new_steps = copy(spec.steps)
        new_constraints = copy(spec.param_constraints)

        # Track first step index per metabolite for
        # equivalence constraints across forms
        first_step_by_met = Dict{Symbol, Int}()

        for de_name in active_de
            entries = de_forms[de_name]
            for (cat_form, met) in entries
                # Determine RE/SS from catalytic step
                cat_idx, cat_is_re = get(
                    cat_step_info, met, (0, true))
                push!(new_steps, StepSpec(
                    [cat_form, met], [de_name],
                    cat_is_re))
                step_idx = length(new_steps)

                if haskey(first_step_by_met, met)
                    first_idx = first_step_by_met[met]
                    # Constrain to first step for same met
                    if cat_is_re
                        push!(new_constraints, (
                            Symbol("K$step_idx"), 1,
                            [(Symbol("K$first_idx"), 1)]))
                    else
                        push!(new_constraints, (
                            Symbol("k$(step_idx)f"), 1,
                            [(Symbol("k$(first_idx)f"), 1)]))
                        push!(new_constraints, (
                            Symbol("k$(step_idx)r"), 1,
                            [(Symbol("k$(first_idx)r"), 1)]))
                    end
                else
                    first_step_by_met[met] = step_idx
                end
            end
        end

        # Mirror steps for each active dead-end metabolite
        all_de_mets = Set{Symbol}()
        active_de_set = Set(active_de)
        for de_name in active_de
            for (_, met) in de_forms[de_name]
                push!(all_de_mets, met)
            end
        end
        for de_met in all_de_mets
            selected_forms = Symbol[]
            for de_name in active_de
                for (form, met) in de_forms[de_name]
                    met == de_met &&
                        push!(selected_forms, form)
                end
            end
            _add_mirror_steps!(
                new_steps, new_constraints,
                spec.steps, selected_forms,
                bound, de_met)
        end

        # For delta=0: constrain dead-end params to
        # catalytic params (only for metabolites that
        # have a catalytic binding step; regulators
        # always add a free K parameter)
        if delta == 0
            for (met, first_idx) in first_step_by_met
                haskey(cat_step_info, met) || continue
                cat_idx, cat_is_re = cat_step_info[met]
                _add_catalytic_equivalence!(
                    new_steps, new_constraints,
                    spec.steps, first_idx, met,
                    cat_is_re)
            end
        end

        candidate = MechanismSpec(
            reaction, new_steps,
            new_constraints, 0)
        pc = _runtime_param_count(candidate)
        pc == target_pc || continue

        push!(result, MechanismSpec(
            reaction, new_steps,
            new_constraints, pc))
    end
end

"""Add mirror steps for dead-end binding.

For each existing step [F1, M] ⇌ [F2], if both F1 and F2 have
dead-end extensions with the same metabolite, add a mirror step
[F1_de, M] ⇌ [F2_de] inheriting RE/SS from the original.
Equivalence constraints tie mirror params to originals.
"""
function _add_mirror_steps!(
    new_steps, new_constraints, orig_steps,
    selected_forms, bound, de_met,
)
    selected_set = Set(selected_forms)
    for (orig_idx, s) in enumerate(orig_steps)
        from = s.reactants[1]
        to = s.products[1]

        from in selected_set || continue
        to in selected_set || continue
        de_met in bound[from] && continue
        de_met in bound[to] && continue

        from_de = _beam_dead_end_form_name(
            bound[from], de_met)
        to_de = _beam_dead_end_form_name(
            bound[to], de_met)

        met = step_metabolite(s)
        if met !== nothing
            push!(new_steps, StepSpec(
                [from_de, met], [to_de],
                s.is_equilibrium))
        else
            push!(new_steps, StepSpec(
                [from_de], [to_de],
                s.is_equilibrium))
        end

        mirror_idx = length(new_steps)
        if s.is_equilibrium
            push!(new_constraints, (
                Symbol("K$mirror_idx"), 1,
                [(Symbol("K$orig_idx"), 1)]))
        else
            push!(new_constraints, (
                Symbol("k$(mirror_idx)f"), 1,
                [(Symbol("k$(orig_idx)f"), 1)]))
            push!(new_constraints, (
                Symbol("k$(mirror_idx)r"), 1,
                [(Symbol("k$(orig_idx)r"), 1)]))
        end
    end
end

"""
    _runtime_denominator_monomials(spec::MechanismSpec)

Derive the rate equation denominator at runtime and return a
concentration fingerprint (`Set{MONO}`) with kinetic symbols
stripped. Requires full symbolic derivation — slower than
`_concentration_fingerprint` but produces the exact fingerprint
after thermodynamic constraint application.
"""
function _runtime_denominator_monomials(spec::MechanismSpec)
    rxn = spec.reaction
    enz_names, enz_set = _spec_enzyme_names(spec)
    rxns = _spec_reactions(spec)
    eq_steps = _spec_eq_steps(spec)
    met_names = _spec_met_names(rxn)
    stoich_mat = _spec_stoich_matrix(spec, met_names, enz_set)
    free_enz_set = _spec_free_enz_set(spec, enz_set)
    binding_Ks = Set(_binding_K_symbols(rxns, eq_steps, enz_set))

    C, rhs_coeffs = _thermodynamic_constraints(
        enz_names, enz_set, rxns, stoich_mat,
        met_names, substrates(rxn), products(rxn),
    )
    dep_exprs, _ = _dependent_param_exprs(
        C, rhs_coeffs, eq_steps, spec.param_constraints,
        binding_Ks, rxns, enz_set, free_enz_set,
    )

    _, denom_terms = _raw_symbolic_rate_polys(
        substrates(rxn), products(rxn),
        enz_names, enz_set,
        rxns, eq_steps, spec.param_constraints,
        dep_exprs,
    )
    _strip_to_concentration_fingerprint(denom_terms)
end

# ─── Deduplication ────────────────────────────────────────────

"""
    _deduplicate_specs(specs, reaction) → Vector{MechanismSpec}

Remove duplicate mechanisms within a level. Two mechanisms are
duplicates if they have the same concentration fingerprint and
constraint descriptor. Keeps the one with lowest param_count.
"""
function _deduplicate_specs(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    isempty(specs) && return specs

    best = Dict{_DedupKey, MechanismSpec}()
    for spec in specs
        steps = spec.steps
        partition = _compute_re_partition_from_steps(steps)
        fp = _concentration_fingerprint(steps, partition)

        groups = Dict{Tuple{Symbol,Bool}, Vector{Int}}()
        for (i, s) in enumerate(steps)
            met = step_metabolite(s)
            met === nothing && continue
            key = (met, s.is_equilibrium)
            push!(get!(groups, key, Int[]), i)
        end
        valid_groups = sort!(
            [sort!(g) for (_, g) in groups
             if length(g) >= 2];
            by=first)

        constraint_mask = _constraints_to_mask(
            spec.param_constraints, valid_groups, steps)
        desc = _constraint_descriptor(
            steps, valid_groups, constraint_mask)

        dedup_key = (fp, desc)
        if !haskey(best, dedup_key) ||
                spec.param_count < best[dedup_key].param_count
            best[dedup_key] = spec
        end
    end
    collect(values(best))
end

# ─── Catalytic Equivalence for +0 Dead-End ────────────────────

"""
    _add_catalytic_equivalence!(
        new_steps, new_constraints, orig_steps,
        first_new_step_idx, met, cat_is_re)

Constrain the first new dead-end binding step's parameters to
equal the catalytic step's parameters for the same metabolite.
For RE steps, constrains K. For SS steps, constrains kf and kr.
"""
function _add_catalytic_equivalence!(
    new_steps, new_constraints,
    orig_steps, first_new_step_idx, met,
    cat_is_re::Bool,
)
    for (i, s) in enumerate(orig_steps)
        step_metabolite(s) == met || continue
        s.is_equilibrium == cat_is_re || continue
        if cat_is_re
            push!(new_constraints, (
                Symbol("K$first_new_step_idx"), 1,
                [(Symbol("K$i"), 1)]))
        else
            push!(new_constraints, (
                Symbol("k$(first_new_step_idx)f"), 1,
                [(Symbol("k$(i)f"), 1)]))
            push!(new_constraints, (
                Symbol("k$(first_new_step_idx)r"), 1,
                [(Symbol("k$(i)r"), 1)]))
        end
        return
    end
end

# ─── Expansion: Same Param Count (+0) ────────────────────────

"""Swap one RE step with one SS step, producing candidates
where param_count is unchanged."""
function _expand_ress_swap!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
)
    constrained = _constrained_symbols(spec.param_constraints)
    target_pc = spec.param_count
    for (i, si) in enumerate(spec.steps)
        for (j, sj) in enumerate(spec.steps)
            i >= j && continue
            # One must be RE, other SS
            si.is_equilibrium == sj.is_equilibrium && continue
            # Skip if constrained params would be invalidated
            re_idx = si.is_equilibrium ? i : j
            ss_idx = si.is_equilibrium ? j : i
            Symbol("K$re_idx") in constrained && continue
            Symbol("k$(ss_idx)f") in constrained && continue
            Symbol("k$(ss_idx)r") in constrained && continue

            new_steps = copy(spec.steps)
            new_steps[i] = StepSpec(
                si.reactants, si.products,
                !si.is_equilibrium)
            new_steps[j] = StepSpec(
                sj.reactants, sj.products,
                !sj.is_equilibrium)

            candidate = MechanismSpec(
                spec.reaction, new_steps,
                copy(spec.param_constraints), 0)
            pc = _runtime_param_count(candidate)
            pc == target_pc || continue
            push!(result, MechanismSpec(
                spec.reaction, new_steps,
                copy(spec.param_constraints), pc))
        end
    end
end

"""
    expand_mechanisms_same_param_count(specs, reaction)
        → Vector{MechanismSpec}

Generate +0 variants: RE↔SS swaps and dead-end configurations
that preserve param_count.

Strategy: first generate all RE↔SS swap variants (fixed point),
then apply dead-end +0 to each variant. Dead-end binding targets
are restricted to catalytic forms only (not dead-end forms from
prior expansions) to avoid combinatorial explosion.
"""
function expand_mechanisms_same_param_count(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    known = Set(
        hash((s.steps, s.param_constraints))
        for s in specs)
    to_expand = copy(specs)
    result = MechanismSpec[]

    while !isempty(to_expand)
        new_this_round = MechanismSpec[]
        for s in to_expand
            candidates = MechanismSpec[]
            _expand_dead_end_at_delta!(
                candidates, s, reaction, 0;
                catalytic_forms=all_form_names(s))
            _expand_ress_swap!(candidates, s)

            for c in candidates
                h = hash((c.steps, c.param_constraints))
                h in known && continue
                push!(known, h)
                push!(new_this_round, c)
                push!(result, c)
            end
        end
        to_expand = new_this_round
    end

    # Structural deduplication: hash-based dedup catches
    # identical specs but not structurally equivalent ones
    # (same rate equation, different step ordering)
    all_specs = vcat(specs, result)
    all_specs = _deduplicate_specs(all_specs, reaction)
    filter(s -> s ∉ specs, all_specs)
end

# ─── Expansion: Allosteric (+2 or more) ──────────────────────

"""
    expand_mechanisms_by_two_params(specs, reaction;
        max_catalytic_n=4) → Vector{AllostericMechanismSpec}

Convert base mechanisms to allosteric with minimum delta (+2):
all metabolites AND all non-binding SS steps TR-equivalent.
The +2 comes from L and K_R_reg. Generates all catalytic_n
values and regulator site partitions, keeping only specs where
`_runtime_param_count == base_param_count + 2`.

Non-TR-equivalent variants (delta +3, +4, ...) are generated
later by `_expand_remove_tr_equiv!` at subsequent beam levels.

For reactions without allosteric regulators, returns empty.
"""
function expand_mechanisms_by_two_params(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction);
    max_catalytic_n::Int=4,
)
    roles = regulator_roles(reaction)
    allo_regs = Symbol[
        r[1] for r in roles
        if r[2] == :allosteric || r[2] == :unknown]
    isempty(allo_regs) && return AllostericMechanismSpec[]

    result = AllostericMechanismSpec[]
    partitions = _set_partitions(allo_regs)

    for spec in specs
        base_pc = spec.param_count
        # All metabolites with T-state params are TR-equivalent
        t_mets = _collect_t_state_metabolites_from_spec(
            spec, allo_regs)
        # All non-binding SS steps are TR-equivalent
        all_ss_steps = _collect_nonbinding_ss_steps_from_spec(
            spec)

        for partition in partitions
            n_groups = length(partition)
            for cn in 1:max_catalytic_n
                for combo in Iterators.product(
                        ntuple(_ -> 1:cn, n_groups)...)
                    allo = AllostericMechanismSpec(
                        spec, cn, partition,
                        collect(combo), copy(t_mets),
                        copy(all_ss_steps),
                        Symbol[], Symbol[], Int[])
                    _runtime_param_count(allo) ==
                        base_pc + 2 || continue
                    push!(result, allo)
                end
            end
        end
    end
    result
end

"""
Collect metabolites that would have T-state binding
parameters for a base MechanismSpec + allosteric regulators.
"""
function _collect_t_state_metabolites_from_spec(
    spec::MechanismSpec,
    allo_regs::Vector{Symbol},
)
    t_mets = Symbol[]
    for s in spec.steps
        s.is_equilibrium || continue
        met = step_metabolite(s)
        met !== nothing && met ∉ t_mets &&
            push!(t_mets, met)
    end
    for reg in allo_regs
        reg ∉ t_mets && push!(t_mets, reg)
    end
    t_mets
end

"""Indices of non-binding SS steps in a MechanismSpec."""
function _collect_nonbinding_ss_steps_from_spec(
    spec::MechanismSpec,
)
    indices = Int[]
    for (i, s) in enumerate(spec.steps)
        if !s.is_equilibrium && step_metabolite(s) === nothing
            push!(indices, i)
        end
    end
    indices
end

# ─── Orchestrator ─────────────────────────────────────────────

"""
    enumerate_mechanisms(reaction; max_param_count=nothing,
        max_catalytic_n=4) → MechanismIterator

Enumerate all valid mechanisms for a reaction, expanding
level-by-level by param_count. At each level:
1. Merge catalytic seeds + cached +2 specs + expanded +1 specs
2. Apply expand_mechanisms_same_param_count to fixed point
3. Deduplicate
4. Yield this level's mechanisms
5. Expand by +1 and +2 for next levels
"""
function enumerate_mechanisms(
    @nospecialize(reaction::EnzymeReaction);
    max_param_count::Union{Nothing,Int}=nothing,
    max_catalytic_n::Int=4,
)
    catalytic = _catalytic_topologies(reaction)

    # Group catalytic topologies by param_count
    seeds_by_pc = Dict{Int, Vector{MechanismSpec}}()
    for spec in catalytic
        push!(get!(seeds_by_pc, spec.param_count,
            MechanismSpec[]), spec)
    end

    min_pc = minimum(keys(seeds_by_pc))
    max_pc = max_param_count

    cache = Dict{Int, Vector{AbstractMechanismSpec}}()
    all_results = AbstractMechanismSpec[]
    current_plus_one = MechanismSpec[]
    current_allo_plus_one = AllostericMechanismSpec[]

    pc = min_pc - 1
    while true
        pc += 1
        max_pc !== nothing && pc > max_pc && break
        # Assemble base MechanismSpec level
        level = MechanismSpec[]
        append!(level, get(seeds_by_pc, pc, MechanismSpec[]))
        append!(level, current_plus_one)

        # Assemble AllostericMechanismSpec level
        allo_level = copy(current_allo_plus_one)

        # Add cached specs
        cached = get(cache, pc, AbstractMechanismSpec[])
        for spec in cached
            if spec isa MechanismSpec
                push!(level, spec)
            elseif spec isa AllostericMechanismSpec
                push!(allo_level, spec)
            end
        end
        delete!(cache, pc)

        # Process allosteric specs
        if !isempty(allo_level)
            allo_deduped = _deduplicate_allosteric(
                allo_level, reaction)
            append!(all_results, allo_deduped)

            current_allo_plus_one =
                expand_mechanisms_by_one_param(
                    allo_deduped, reaction)
            # Filter to actual pc + 1
            filter!(current_allo_plus_one) do s
                _runtime_param_count(s) == pc + 1
            end
        else
            current_allo_plus_one = AllostericMechanismSpec[]
        end

        if isempty(level)
            # Still need to check termination
            has_future = _has_future_work(
                seeds_by_pc, cache,
                current_plus_one,
                current_allo_plus_one, pc, max_pc)
            has_future || break
            current_plus_one = MechanismSpec[]
            continue
        end

        # +0 expansion (fixed point handled internally)
        new_zero = expand_mechanisms_same_param_count(
            level, reaction)
        append!(level, new_zero)
        level = _deduplicate_specs(level, reaction)

        # Yield this level
        append!(all_results, level)

        # Expand +1 (base mechanisms)
        current_plus_one = expand_mechanisms_by_one_param(
            level, reaction)
        filter!(s -> s.param_count == pc + 1,
            current_plus_one)
        current_plus_one = _deduplicate_specs(
            current_plus_one, reaction)

        # Expand +2 (allosteric)
        plus_two = expand_mechanisms_by_two_params(
            level, reaction; max_catalytic_n)
        for spec in plus_two
            target_pc = _runtime_param_count(spec)
            push!(get!(cache, target_pc,
                AbstractMechanismSpec[]), spec)
        end

        # Termination check
        has_future = _has_future_work(
            seeds_by_pc, cache,
            current_plus_one,
            current_allo_plus_one, pc, max_pc)
        has_future || break
    end

    MechanismIterator(all_results, length(all_results))
end

"""Check if there's any work remaining at future levels."""
function _has_future_work(
    seeds_by_pc, cache,
    current_plus_one, current_allo_plus_one,
    pc, max_pc::Union{Nothing,Int},
)
    !isempty(current_plus_one) && return true
    !isempty(current_allo_plus_one) && return true
    for k in keys(seeds_by_pc)
        k > pc && (max_pc === nothing || k <= max_pc) &&
            return true
    end
    for k in keys(cache)
        k > pc && (max_pc === nothing || k <= max_pc) &&
            return true
    end
    false
end
