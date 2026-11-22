# ABOUTME: A parent candidate's children: every expansion move applied to it, with the
# ABOUTME: emission and atom assertions, the regulator-type filter, the guard, the limits.

# Each child must pass the regulator-type filter (`_respects_reg_type`) and conserve atoms
# (`_assert_atom_conserving`, which errors on the first child that does not); a child that
# `_enumerator_keeps` refuses (the guard rejects it, or a conformation keeps a form that
# only bindings it excludes enter) is dropped. The parent must satisfy the two emission
# rules (`_assert_emission_rules`): every steady-state group carries flux and no
# competitive-inhibitor group is redundant (`_redundant_copy_groups`).
function candidate_mechanisms(reaction::Mechanisms.EnzymeReaction,
                              parent::CandidateMechanism;
                              max_param_count::Int = typemax(Int),
                              eq_complexity_filter::Int = typemax(Int))
    _assert_emission_rules(parent)
    children = CandidateMechanism[]
    append!(children, _expand_re_to_ss(parent))
    append!(children, _expand_split_kinetic_group(parent))
    append!(children, _expand_add_dead_end_regulator(parent, reaction))
    if !_is_allosteric(parent.mechanism)
        append!(children, _expand_to_allosteric(parent, reaction))
    else
        append!(children, _expand_add_allosteric_regulator(parent, reaction),
                _expand_change_allo_state(parent), _expand_merge_regulatory_sites(parent))
    end
    filter!(c -> _respects_reg_type(c.mechanism, reaction), children)
    foreach(child -> _assert_atom_conserving(child, reaction), children)
    # The term count reads the graph alone, where the guard derives, so the term limit
    # runs first and `_within_size_limits` applies only the other.
    eq_complexity_filter == typemax(Int) || filter!(
        c -> Derivation.rate_equation_term_count(c.mechanism) <= eq_complexity_filter,
        children)
    filter!(c -> _enumerator_keeps(c.mechanism), children)
    _within_size_limits(children, max_param_count, typemax(Int))
end
