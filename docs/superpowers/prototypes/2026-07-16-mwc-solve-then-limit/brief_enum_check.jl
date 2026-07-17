# ABOUTME: Does the package let you build the P:OnlyA + active-I-catalysis trap combo?
# ABOUTME: The AllostericMechanism constructor validates via _onlya_haldane_violation; test the trap + repairs.
using EnzymeRates
const ER = EnzymeRates

# helper: try to build a uni-uni S<=>P allosteric mechanism with given (S,cat,P) tags
function try_build(Stag, cattag, Ptag)
    expr = Meta.parse("""
    ER.@allosteric_mechanism begin
        substrates: S ; products: P ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + S ⇌ E(S)     :: $Stag
            E(S) <--> E(P)   :: $cattag
            E + P ⇌ E(P)     :: $Ptag
        end
    end
    """)
    try
        Core.eval(@__MODULE__, expr)
        return (ok=true, msg="built OK")
    catch e
        s = sprint(showerror, e)
        # Constructor rejection carries "AllostericMechanism:"; anything else means
        # the mechanism CONSTRUCTED and a later (accessor) call tripped.
        if occursin("AllostericMechanism:", s)
            return (ok=false, msg=first(split(s, '\n')))
        else
            return (ok=true, msg="built OK (post-construct accessor noise: $(first(split(s,'\n'))[1:min(end,40)]))")
        end
    end
end

cases = [
    ("TRAP     (S:EqualAI,   cat:NonequalAI, P:OnlyA)", :EqualAI,    :NonequalAI, :OnlyA),
    ("TRAP     (S:NonequalAI,cat:NonequalAI, P:OnlyA)", :NonequalAI, :NonequalAI, :OnlyA),
    ("TRAP     (S:EqualAI,   cat:EqualAI,    P:OnlyA)", :EqualAI,    :EqualAI,    :OnlyA),
    ("repair-1 (S:EqualAI,   cat:OnlyA,      P:OnlyA)", :EqualAI,    :OnlyA,      :OnlyA),
    ("repair-2 (S:OnlyA,     cat:EqualAI,    P:OnlyA)", :OnlyA,      :EqualAI,    :OnlyA),
    ("valid FamilyA (S:OnlyA,cat:OnlyA,      P:EqualAI)", :OnlyA,    :OnlyA,      :EqualAI),
    ("valid no-trap (S:EqualAI,cat:NonequalAI,P:EqualAI)", :EqualAI, :NonequalAI, :EqualAI),
    ("proto-INCONSISTENT (S:OnlyA,cat:EqualAI,P:EqualAI)", :OnlyA,   :EqualAI,    :EqualAI),
    ("valid TR    (S:EqualAI,   cat:OnlyA,   P:NonequalAI)", :EqualAI,:OnlyA,     :NonequalAI),
]

println("="^92)
println("Can the package CONSTRUCT each (S,cat,P) tag combo?  (constructor validates via _onlya_haldane_violation)")
println("="^92)
for (label, S, c, P) in cases
    r = try_build(S, c, P)
    status = r.ok ? "CONSTRUCTS" : "REJECTED  "
    println(rpad(label, 50), " -> ", status, "  ", r.ok ? r.msg : r.msg[1:min(end,120)])
end
