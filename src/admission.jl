# ABOUTME: The well-posedness check the mechanism macros, the enumerator's exits and the
# ABOUTME: derivation entries run on a mechanism; constructors check representation only.

"""
    _check_well_posed(m)

Error unless `m`'s rate law is defined where data sit: no rapid-equilibrium segment
without a bottom form (`_assert_bottomed_re_segments`), in the catalytic graph and, for
an allosteric mechanism, in the inactive conformation's graph (`_inactive_groups`); every
free form in one rapid-equilibrium segment of the steps outside `:OnlyA` groups
(`_assert_one_flip_segment`); no `:OnlyA` binding that leaves a thermodynamic cycle
unsatisfiable (`_onlya_haldane_violation`). The mechanism macros, `init_mechanisms`,
`seed_mechanisms`, `expand_mechanisms`, identify's base tier and every derivation entry
(`_derivation_input`) call it; the constructors and the Sig decoder do not, so code that
builds a mechanism internally never trips a check meant for the mechanisms it returns.
"""
function _check_well_posed(m::Mechanism)
    _assert_bottomed_re_segments(m, steps(m))
    nothing
end
function _check_well_posed(am::AllostericMechanism)
    _assert_bottomed_re_segments(am, steps(am))
    _assert_one_flip_segment(am)
    offenders = _onlya_haldane_violation(am, steps(am), cat_allo_states(am))
    isempty(offenders) ||
        error("AllostericMechanism: ", _onlya_violation_message(am, offenders))
    _assert_bottomed_re_segments(am, _inactive_groups(am))
    nothing
end

"""
Error when a rapid-equilibrium segment of `groups`, kinetic groups over mechanism `m`'s
metabolites, has no bottom form (`_bottomless_re_segment`): its rate is 0/0 where one
side's metabolites are absent. The boundary check (`_check_well_posed`) checks a
mechanism's groups and, for an allosteric mechanism, the inactive conformation's pruned
graph (`_inactive_groups`).
"""
function _assert_bottomed_re_segments(m::_AnyMech, groups::Vector{Vector{Step}})
    segment = _bottomless_re_segment(m, groups)
    segment === nothing ||
        error("Mechanism: rapid-equilibrium segment {" *
              join(_form_label.(segment, Ref(m)), ", ") *
              "} has no form free of one side's " *
              "metabolites, so its rate is undefined when that side is absent " *
              "(0/0). Make one of the segment's binding steps steady-state.")
    nothing
end

"""
Error unless every free form of the allosteric mechanism `am`'s catalytic steps (a form
with no bound metabolite and no residual) lies in one rapid-equilibrium segment of the
steps outside `:OnlyA` groups. Under formulation 1 the conformations interconvert at a
free form. Within one such segment, flipping at any of its free forms gives the same
model with L rescaled; free forms joined any other way flip at different forms in the two
conformations or leave one stranded in the inactive conformation, a different model.
"""
function _assert_one_flip_segment(am::AllostericMechanism)
    groups, tags = steps(am), cat_allo_states(am)
    free = [f for f in _forms(groups) if isempty(bound(f)) && isempty(residual(f))]
    length(free) <= 1 && return nothing
    re = [s for (g, group) in enumerate(groups) for s in group
          if tags[g] !== :OnlyA && is_equilibrium(s)]
    segment, changed = Set([first(free)]), true
    while changed
        changed = false
        for s in re
            a, b = from_species(s), to_species(s)
            (a in segment) == (b in segment) && continue
            push!(segment, a, b)
            changed = true
        end
    end
    all(in(segment), free) || error(
        "AllostericMechanism: the free forms $(join(_form_label.(free, Ref(am)), ", ")) " *
        "are not joined by rapid-equilibrium steps outside :OnlyA groups. The " *
        "conformations interconvert at one free form, so the other free forms must be " *
        "in rapid equilibrium with it: join them by ⇌ steps outside :OnlyA groups, or " *
        "keep one free form.")
    nothing
end

"""`m`, once `_check_well_posed(m)` passes."""
_well_posed(m) = (_check_well_posed(m); m)

"""The concrete mechanism the compiled type `M` encodes (`_concrete`), checked by
`_check_well_posed`: the input of every derivation entry."""
_derivation_input(@nospecialize(M::Type{<:AbstractEnzymeMechanism})) =
    _well_posed(_concrete(M()))
