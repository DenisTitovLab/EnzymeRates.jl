# ABOUTME: The well-posedness check the mechanism macros, the enumerator's exits and the
# ABOUTME: derivation entries run on a mechanism; constructors check representation only.

"""
    _check_well_posed(m)

Error unless `m`'s rate law is defined where data sit: no rapid-equilibrium segment
without a bottom form (`_assert_bottomed_segments`), in the catalytic graph and, for an
allosteric mechanism, in the inactive conformation's graph (`_inactive_mask`); every
free form in one rapid-equilibrium segment of the edges outside `:OnlyA` groups
(`_assert_one_flip_segment`); no `:OnlyA` binding that leaves a thermodynamic cycle
unsatisfiable (`_onlya_haldane_violation`). It reads the reaction graph through masks
and builds nothing. The mechanism macros, `init_mechanisms`, `seed_mechanisms`,
`expand_mechanisms`, identify's base tier and every derivation entry
(`_derivation_input`) call it; the constructors and the Sig decoder do not, so code that
builds a mechanism internally never trips a check meant for the mechanisms it returns.
"""
function _check_well_posed(m::Mechanism)
    _assert_bottomed_segments(m, trues(length(m.edges)))
    nothing
end
function _check_well_posed(am::AllostericMechanism)
    m = am.catalytic
    _assert_bottomed_segments(m, trues(length(m.edges)))
    _assert_one_flip_segment(am)
    offenders = _onlya_haldane_violation(am)
    isempty(offenders) ||
        error("AllostericMechanism: ", _onlya_violation_message(am, offenders))
    _assert_bottomed_segments(m, _inactive_mask(am))
    nothing
end

"""
    _bottomless_segment(m, mask, re) -> Vector{Int}

The forms, breadth-first (`_segments`), of the first rapid-equilibrium segment of the
edges `mask` keeps, joined by those `re` flags, whose weights all vanish at zero
products, or all vanish at zero substrates; empty when no segment does. Each metabolite
name takes its side from `m`'s substrates and products, so a modifier copy, which shares
its reactant's name and concentration, counts on that reactant's side even in a
sub-graph that holds no edge taking up or giving off the reactant itself (the inactive
conformation's, `_inactive_mask`).

Within a segment, rapid equilibrium fixes each form's weight relative to any
other as a monomial in concentrations (one factor per RE binding step on the
path between them). Once the weights are scaled to clear every negative
exponent, a segment normally has a bottom form whose weight is free of one
side's metabolites — the form the others drain to when that side is absent. A
segment with no such form (e.g. `{E(P), E(Q), E(P, Q)}` with `E + P` and
`E + Q` at steady state: weights `Kq·P : Kp·Q : P·Q`) is empty at zero
products, and the rate there depends on how `E(P, Q)` splits between `E(P)`
and `E(Q)` — a ratio of two fast release rates the RE approximation does not
carry. The rate equation is `0/0` exactly where initial-rate data sit. A
segment whose weights vanish only at a corner mixing a substrate and a product
(a mixed abortive complex, or the ping-pong `E` / `E(; residual)` pair) is
harmless: no turnover is possible at that corner.
"""
function _bottomless_segment(m::Mechanism, mask::BitVector, re::BitVector)
    s = _segments(m, mask, re)
    name_of = _name_index(m)
    # Whether each metabolite name (`_name_index`) is a product, and whether a substrate.
    products, substrates = falses(size(s.offsets, 1)), falses(size(s.offsets, 1))
    for (i, x) in enumerate(m.table)
        products[name_of[i]] = insorted(x.name, m.products)
        substrates[name_of[i]] = insorted(x.name, m.substrates)
    end
    for segment in s.segments, on_side in (products, substrates)
        length(segment) >= 2 &&
            all(f -> any(n -> on_side[n] && s.offsets[n, f] > 0, eachindex(on_side)),
                segment) && return segment
    end
    Int[]
end

"""
Error when a rapid-equilibrium segment of the edges of `m` that `mask` keeps has no
bottom form (`_bottomless_segment`): its rate is 0/0 where one side's metabolites are
absent. The boundary check (`_check_well_posed`) checks a mechanism's graph and, for an
allosteric mechanism, the inactive conformation's (`_inactive_mask`).
"""
function _assert_bottomed_segments(m::Mechanism, mask::BitVector)
    segment = _bottomless_segment(m, mask, _flags(m))
    isempty(segment) ||
        error("Mechanism: rapid-equilibrium segment {" *
              join([_form_name(m, f) for f in segment], ", ") *
              "} has no form free of one side's " *
              "metabolites, so its rate is undefined when that side is absent " *
              "(0/0). Make one of the segment's binding steps steady-state.")
    nothing
end

"""
Error unless every free form (no bound metabolite, no residual) of the allosteric
mechanism `am`'s catalytic graph lies in one rapid-equilibrium segment of the edges
outside `:OnlyA` groups. Under formulation 1 the conformations interconvert at a free
form. Within one such segment, flipping at any of its free forms gives the same model
with L rescaled; free forms joined any other way flip at different forms in the two
conformations or leave one stranded in the inactive conformation, a different model.
"""
function _assert_one_flip_segment(am::AllostericMechanism)
    m = am.catalytic
    free = [f for f in eachindex(m.forms) if _is_free(m.forms[f])]
    length(free) <= 1 && return nothing
    segment_of = _segments(m, trues(length(m.edges)),
                           _flags(m) .& _outside_onlya(am)).segment_of
    all(f -> segment_of[f] == segment_of[first(free)], free) || error(
        "AllostericMechanism: the free forms " *
        "$(join([_form_name(m, f) for f in free], ", ")) " *
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
