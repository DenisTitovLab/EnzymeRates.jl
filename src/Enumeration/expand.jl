# ABOUTME: expand_mechanisms: every expansion move applied to each parent, with the
# ABOUTME: emission and atom assertions, the regulator-type filter and the guard.

"""
    expand_mechanisms(mechs, reaction) -> Vector{Union{Mechanism, AllostericMechanism}}

Apply all expansion moves (RE→SS, split kinetic group, add dead-end
regulator, to-allosteric, add allosteric regulator, change allo state,
merge regulatory sites) to each input mechanism and return the children as
a flat vector. Bucketing by parameter count is the caller's job, not
enumeration's. Each parent must satisfy the two emission rules
(`_assert_emission_rules`): every steady-state group carries flux and no
competitive-inhibitor group is redundant (`_redundant_copy_groups`). Every child that
passes the regulator-type filter must conserve atoms (`_assert_atom_conserving`), which
errors on the first child that does not; a child whose rate law the guard rejects
(`has_valid_rate_equation`) is dropped.
"""
function expand_mechanisms(
    mechs::Vector{<:Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}},
    rxn::Mechanisms.EnzymeReaction)
    result = Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}[]
    for m in mechs
        _assert_emission_rules(m)
        append!(result, _expand_re_to_ss(m), _expand_split_kinetic_group(m),
                _expand_add_dead_end_regulator(m, rxn))
        if m isa Mechanisms.Mechanism
            append!(result, _expand_to_allosteric(m, rxn))
        else
            append!(result, _expand_add_allosteric_regulator(m, rxn),
                    _expand_change_allo_state(m), _expand_merge_regulatory_sites(m))
        end
    end
    filter!(c -> _respects_reg_type(c, rxn), result)
    foreach(child -> _assert_atom_conserving(child, rxn), result)
    filter!(Derivation.has_valid_rate_equation, result)
    result
end
