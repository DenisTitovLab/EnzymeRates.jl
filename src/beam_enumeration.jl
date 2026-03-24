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
