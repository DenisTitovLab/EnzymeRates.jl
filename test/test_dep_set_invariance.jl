# ABOUTME: Guard for the kinetic-group naming representative.
# ABOUTME: Proves the Haldane/Wegscheider dependent-parameter CHOICE is
# ABOUTME: invariant to which step represents a kinetic group.

# The kinetic-group naming representative must not change WHICH kinetic group
# becomes the Haldane elimination's dependent parameter — only the rendered
# name of the representative.
#
# A dependent parameter names a kinetic GROUP. We identify each group
# structurally and order-independently (sorted member step-hashes), so the
# key is invariant to step order within the group — and therefore invariant
# to the rep choice. Permuting the step order inside every group (reversing
# it) changes the rep (and thus the rendered Symbol) but MUST leave the set
# of dependent groups untouched. If the dependent-parameter choice depends on
# naming/order, this set differs.

"""
Structural, order-independent identity of the dependent parameter `sym`:
`(typeof(p), group_identity, p.state)`. `group_identity` is the hash of the
sorted member step-hashes of the kinetic group the parameter governs (Kreg:
the site hash). Recovers the Parameter from the Symbol via `_param_for_symbol`,
flipping an `_I_`-state name back to its `_A_` source when needed.
"""
function _dep_struct_key(sym::Symbol, mech::EnzymeRates.Mechanism)
    p = try
        EnzymeRates._param_for_symbol(mech, sym)
    catch
        active_sym = Symbol(replace(String(sym), "_I_" => "_A_"; count=1))
        EnzymeRates._flip_to_inactive(
            EnzymeRates._param_for_symbol(mech, active_sym))
    end
    group_id = if p isa EnzymeRates.Kreg
        hash(p.site)
    else
        gh = UInt(0)
        for group in EnzymeRates.steps(mech)
            if p.step in group
                gh = hash(sort!([hash(s) for s in group]))
                break
            end
        end
        gh
    end
    return (typeof(p), group_id, p.state)
end

"""Set of structural dep-group keys for a Mechanism."""
function _dep_struct_key_set(mech::EnzymeRates.Mechanism)
    M = typeof(EnzymeRates.compile_mechanism(mech))
    dep_exprs, _ = EnzymeRates._dependent_param_exprs(M)
    return Set(_dep_struct_key(sym, mech) for sym in keys(dep_exprs))
end

@testset "dependent-param choice invariant to group-rep" begin
    for spec in MECHANISM_TEST_SPECS
        spec.mechanism isa EnzymeRates.EnzymeMechanism || continue
        @testset "$(spec.name)" begin
            mech = EnzymeRates.Mechanism(spec.mechanism)
            base_keys = _dep_struct_key_set(mech)
            # Reverse step order within every kinetic group. This flips the
            # naming rep wherever a group has >1 step, but must not change
            # which groups are dependent.
            reversed = [reverse(g) for g in EnzymeRates.steps(mech)]
            mech_rev = EnzymeRates.Mechanism(mech.reaction, reversed)
            rev_keys = _dep_struct_key_set(mech_rev)
            @test base_keys == rev_keys
        end
    end
end
