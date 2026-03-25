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
    _, indep = _dependent_param_exprs(
        C, rhs_coeffs, eq_steps, spec.param_constraints,
        binding_Ks, rxns, enz_set, free_enz_set,
    )
    # indep params + Keq + E_total
    length(indep) + 2
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
    _dead_end_opportunities(spec, reaction)
        → Vector{Tuple{Symbol, Symbol}}

Find (form, metabolite) pairs where dead-end binding is possible.
A binding is valid when:
- The metabolite is not already bound at the form
- The result doesn't exceed binding capacity
- The result is not already a form in the mechanism
"""
function _dead_end_opportunities(
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
)
    bound = _bound_entities(spec)
    cap = _binding_capacity(reaction)
    existing_forms = all_form_names(spec)

    sub_names = Set(s[1] for s in substrates(reaction))
    prod_names = Set(p[1] for p in products(reaction))
    all_mets = union(sub_names, prod_names)

    roles = regulator_roles(reaction)
    de_regs = Symbol[
        r[1] for r in roles
        if r[2] == :dead_end || r[2] == :unknown]

    opportunities = Tuple{Symbol, Symbol}[]

    for form in sort!(collect(keys(bound)))
        fb = bound[form]
        n_bound = length(fb)
        n_bound >= cap && continue

        # Substrate/product dead-end opportunities
        for met in sort!(collect(all_mets))
            met in fb && continue
            de_name = _beam_dead_end_form_name(fb, met)
            de_name in existing_forms && continue
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
    _expand_dead_end_at_delta!(result, spec, reaction, delta)

Shared dead-end expansion logic. For each metabolite/regulator,
tries binding to 1, 2, ... forms with maximal equivalence
constraints. Keeps only configurations where computed param_count
equals spec.param_count + delta.
"""
function _expand_dead_end_at_delta!(
    result::Vector{MechanismSpec},
    spec::MechanismSpec,
    @nospecialize(reaction::EnzymeReaction),
    delta::Int,
)
    opps = _dead_end_opportunities(spec, reaction)
    isempty(opps) && return

    bound = _bound_entities(spec)

    by_met = Dict{Symbol, Vector{Symbol}}()
    for (form, met) in opps
        push!(get!(by_met, met, Symbol[]), form)
    end

    target_pc = spec.param_count + delta

    for (met, forms) in by_met
        n = length(forms)
        for mask in 1:(1 << n) - 1
            selected = Symbol[
                forms[i] for i in 1:n
                if (mask >> (i - 1)) & 1 == 1]

            new_steps = copy(spec.steps)
            new_constraints = copy(spec.param_constraints)

            first_step_idx = nothing
            for (j, form) in enumerate(selected)
                de_name = _beam_dead_end_form_name(
                    bound[form], met)
                push!(new_steps, StepSpec(
                    [form, met], [de_name], true))
                step_idx = length(new_steps)

                if j == 1
                    first_step_idx = step_idx
                else
                    # Constrain K of this step to match
                    # K of first step (shared K_R)
                    push!(new_constraints, (
                        Symbol("K$step_idx"), 1,
                        [(Symbol("K$first_step_idx"), 1)]))
                end
            end

            _add_mirror_steps!(
                new_steps, new_constraints, spec.steps,
                selected, bound, met)

            if delta == 0
                _add_catalytic_equivalence!(
                    new_steps, new_constraints,
                    spec.steps, first_step_idx, met)
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
        fp = _runtime_denominator_monomials(spec)

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
        first_new_step_idx, met)

Constrain the first new dead-end binding step's K to equal the
catalytic step's K for the same metabolite. Only matches RE
catalytic steps (SS steps have kf/kr, not K).
"""
function _add_catalytic_equivalence!(
    new_steps, new_constraints,
    orig_steps, first_new_step_idx, met,
)
    for (i, s) in enumerate(orig_steps)
        s.is_equilibrium || continue
        step_metabolite(s) == met || continue
        # Constrain new dead-end K to catalytic K
        push!(new_constraints, (
            Symbol("K$first_new_step_idx"), 1,
            [(Symbol("K$i"), 1)]))
        return
    end
end

# ─── Expansion: Same Param Count (+0) ────────────────────────

"""
    expand_mechanisms_same_param_count(specs, reaction)
        → Vector{MechanismSpec}

Add dead-end configurations that result in +0 net parameter
change. Iterates to a fixed point: each pass may create new
forms that enable further +0 additions. Returns only the
newly discovered +0 variants (deduplicated).
"""
function expand_mechanisms_same_param_count(
    specs::Vector{MechanismSpec},
    @nospecialize(reaction::EnzymeReaction),
)
    all_specs = copy(specs)
    prev_count = 0
    while length(all_specs) != prev_count
        prev_count = length(all_specs)
        new_zero = MechanismSpec[]
        for spec in all_specs
            _expand_dead_end_at_delta!(
                new_zero, spec, reaction, 0)
        end
        append!(all_specs, new_zero)
        all_specs = _deduplicate_specs(
            all_specs, reaction)
    end
    filter(s -> s ∉ specs, all_specs)
end
