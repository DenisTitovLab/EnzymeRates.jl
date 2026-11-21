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
    m = parent.mechanism
    _assert_emission_rules(m)
    children = Union{Mechanisms.Mechanism, Mechanisms.AllostericMechanism}[]
    # One call per move: on a `Union` parent, three moves' results in one call exceed what
    # inference splits, and it dispatches at run time.
    append!(children, _expand_re_to_ss(m))
    append!(children, _expand_split_kinetic_group(m))
    append!(children, _expand_add_dead_end_regulator(m, reaction))
    if m isa Mechanisms.Mechanism
        append!(children, _expand_to_allosteric(m, reaction))
    else
        append!(children, _expand_add_allosteric_regulator(m, reaction),
                _expand_change_allo_state(m), _expand_merge_regulatory_sites(m))
    end
    filter!(c -> _respects_reg_type(c, reaction), children)
    foreach(child -> _assert_atom_conserving(child, reaction), children)
    filter!(_enumerator_keeps, children)
    _candidates(children, max_param_count, eq_complexity_filter)
end
