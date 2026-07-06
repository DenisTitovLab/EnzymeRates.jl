# ABOUTME: Unit tests for _split_resolution — the honorable-split nullspace partition.
module SplitResolutionTests
using Test, EnzymeRates
const ER = EnzymeRates
const Sub=ER.Substrate; const Prd=ER.Product; const Sp=ER.Species
const St=ER.Step; const RA=ER.ReactantAtoms; const Met=ER.Metabolite
const S=Sub(:S); const P=Prd(:P); const A=Sub(:A); const B=Sub(:B)

# uni-uni: E+S->ES (bind), ES<->EP (SS catalysis), EP->E+P (release)
function uni(states)
    E=Sp(Met[],:E); ES=Sp(Met[S],:E); EP=Sp(Met[P],:E)
    rxn=ER.EnzymeReaction(RA[RA(S,[:C=>1]),RA(P,[:C=>1])], ER.RegulatorMults[], Int[2])
    steps=Vector{St}[[St(E,ES,S,true)],[St(ES,EP,nothing,false)],[St(EP,E,P,true)]]
    ER.AllostericMechanism(rxn, steps, collect(Symbol,states), 2, ER.RegulatorySite[])
end
gi(am) = [g for g in 1:length(ER.steps(am)) if ER.cat_allo_state(am,g)===:NonequalAI]

@testset "_split_resolution" begin
    # (a) single NonequalAI binding + EqualAI catalysis -> fully forbidden: NO free split.
    am = uni([:NonequalAI, :EqualAI, :EqualAI])
    r = ER._split_resolution(am)
    @test isempty(r.free)
    @test length(r.derived) == 1
    @test first(r.derived).first in gi(am)         # the S-binding group is derived (K_I=K_A)
    @test isempty(first(r.derived).second)          # derived from NO free split => K_I_S=K_A_S

    # (b) two NonequalAI bindings + EqualAI catalysis -> 1 honorable DOF.
    am2 = uni([:NonequalAI, :EqualAI, :NonequalAI])  # S-binding + P-binding NonequalAI
    r2 = ER._split_resolution(am2)
    @test length(r2.free) == 1
    @test length(r2.derived) == 1
    # the derived group's split is +1 or -1 times the free group's split (delta_P = delta_S).
    d = first(r2.derived)
    @test length(d.second) == 1
    @test d.second[1].first == r2.free[1]
    @test abs(d.second[1].second) == 1

    # (c) catalysis NonequalAI -> its reverse differs natively; the binding split is free, no collapse.
    am3 = uni([:NonequalAI, :NonequalAI, :EqualAI])
    r3 = ER._split_resolution(am3)
    @test length(r3.free) == 2      # both S-binding and catalysis keep free splits
    @test isempty(r3.derived)

    # (d) all EqualAI -> no NonequalAI groups -> empty resolution.
    r4 = ER._split_resolution(uni([:EqualAI,:EqualAI,:EqualAI]))
    @test isempty(r4.free) && isempty(r4.derived)
end
end # module
