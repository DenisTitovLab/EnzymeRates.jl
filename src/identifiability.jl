"""
Structural identifiability analysis for enzyme rate equations.

NOTE: The previous scaling symmetry analysis has been removed because
QSSA enzyme rate equations are inherently scaling-identifiable due to
the topological structure of spanning trees in the King-Altman method.

Future work: Implement permutation symmetry detection for mechanisms
with equivalent pathways (e.g., symmetric binding sites).
"""

# ─── Structural identifiability deficit detection ────────────────────────────

"""
    _count_rate_monomials(M::Type{<:EnzymeMechanism}) → (n_num, n_denom)

Count distinct metabolite monomials in numerator and denominator of rate equation.

Returns the number of unique monomials (terms with distinct metabolite combinations,
ignoring k-parameters, Keq, and E_total) in each polynomial.

Note: Uses the raw symbolic expression (before Haldane substitution) to count monomials,
because Haldane substitution introduces Keq which doesn't change the structural form.
"""
function _count_rate_monomials(M::Type{<:EnzymeMechanism})
    # Use raw expression to avoid Keq appearing in monomials
    expr = _raw_symbolic_rate_expr(M)
    num_poly = _to_poly(_strip_div(expr))
    denom_expr = _find_denom(expr)
    denom_poly = denom_expr === nothing ? _Poly(Symbol[] => 1) : _to_poly(denom_expr)

    # Extract unique metabolite monomials (filter out k-parameters, Keq, and E_total)
    function metabolite_monomial(k::Vector{Symbol})
        sort(filter(s -> !is_k_parameter(s) && s != :E_total && s != :Keq, k))
    end

    num_monomials = unique([metabolite_monomial(k) for (k, v) in num_poly])
    denom_monomials = unique([metabolite_monomial(k) for (k, v) in denom_poly])

    return length(num_monomials), length(denom_monomials)
end

"""
    structural_identifiability_deficit(m::EnzymeMechanism) → Int

Number of excess parameters beyond what is structurally identifiable from kinetic data.

The deficit is computed as:
    n_k - (n_num + n_denom - 2)

where:
- n_k = number of independent k-parameters (after Haldane/Wegscheider constraints)
- n_num = number of distinct metabolite monomials in numerator
- n_denom = number of distinct metabolite monomials in denominator

The `-2` accounts for:
- `-1` from Haldane constraining the ratio of numerator coefficients (a_forward/a_reverse = Keq)
- `-1` from scaling freedom (can normalize one denominator coefficient)

Returns:
- `0` if exactly identifiable (parameters match identifiable coefficients)
- `> 0` if underdetermined (more parameters than identifiable coefficients)
- `< 0` if overdetermined (fewer parameters than identifiable coefficients)

# Examples
```julia
m_uu = make_uni_uni()  # E + S ⇌ ES ⇌ E + P
structural_identifiability_deficit(m_uu)  # Returns 0 (exactly identifiable)

m_3step = make_three_step_isomerization()  # E + S ⇌ ES ⇌ ES' ⇌ E + P
structural_identifiability_deficit(m_3step)  # Returns 2 (underdetermined)
```
"""
@generated function structural_identifiability_deficit(::M) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    n_k = length(indep)  # After Haldane/Wegscheider on k's

    n_num, n_denom = _count_rate_monomials(M)

    # Structurally identifiable coefficients:
    # - Numerator: n_num - 1 (Haldane constrains ratio a_forward/a_reverse = Keq)
    # - Denominator: n_denom - 1 (scaling freedom allows normalizing one coefficient)
    n_identifiable = (n_num - 1) + (n_denom - 1)

    return n_k - n_identifiable
end

"""
    is_structurally_identifiable(m::EnzymeMechanism) → Bool

Check if all mechanism parameters can be uniquely determined from steady-state kinetic data.
Returns `true` if structurally identifiable, `false` otherwise.

A mechanism is structurally identifiable when the number of independent k-parameters
is at most the number of structurally identifiable coefficients (deficit ≤ 0).

# Examples
```julia
m_uu = make_uni_uni()  # E + S ⇌ ES ⇌ E + P
is_structurally_identifiable(m_uu)  # Returns true

m_3step = make_three_step_isomerization()  # E + S ⇌ ES ⇌ ES' ⇌ E + P
is_structurally_identifiable(m_3step)  # Returns false
```
"""
is_structurally_identifiable(m::EnzymeMechanism) = structural_identifiability_deficit(m) <= 0

# ─── IdentifiableHaldaneWegscheider mode ─────────────────────────────────────
#
# Since scaling non-identifiability cannot occur in QSSA mechanisms,
# IdentifiableHaldaneWegscheider is equivalent to HaldaneWegscheider.

@generated function rate_equation(
    m::M, params::NamedTuple, concs::NamedTuple, ::IdentifiableHaldaneWegscheiderMode
) where {M <: EnzymeMechanism}
    _symbolic_rate_expr(M)
end

@generated function parameters(::M, ::IdentifiableHaldaneWegscheiderMode) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    return (indep..., :Keq, :E_total)
end

function rate_equation_string(m::M, ::IdentifiableHaldaneWegscheiderMode) where {M<:EnzymeMechanism}
    rate_equation_string(m, HaldaneWegscheider)
end

# ─── Placeholder for future permutation symmetry detection ───────────────────

"""
    has_permutation_symmetry(m::EnzymeMechanism) → Bool

Check if the mechanism has permutation-equivalent parameters
(e.g., symmetric binding sites where swapping pathways preserves the rate).

NOT YET IMPLEMENTED - always returns false.
"""
function has_permutation_symmetry(::EnzymeMechanism)
    return false  # TODO: implement
end

"""
    permutation_symmetric_params(m::EnzymeMechanism) → Vector{Vector{Symbol}}

Return groups of parameters that can be permuted without changing the rate equation.

NOT YET IMPLEMENTED - always returns empty vector.
"""
function permutation_symmetric_params(::EnzymeMechanism)
    return Vector{Symbol}[]  # TODO: implement
end
