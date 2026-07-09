# ABOUTME: Strict :EqualAI collapse regression — no silent I-name promotion; forbidden
# ABOUTME: splits collapse to K_A=K_I, honorable splits stay free, all equations thermo-consistent.
module AllostericCollapseTests
using Test, EnzymeRates, Random
const ER = EnzymeRates
const Sub=ER.Substrate; const Prd=ER.Product; const Sp=ER.Species
const St=ER.Step; const RA=ER.ReactantAtoms; const Met=ER.Metabolite
const S=Sub(:S); const P=Prd(:P)
function uni(states)
    E=Sp(Met[],:E); ES=Sp(Met[S],:E); EP=Sp(Met[P],:E)
    rxn=ER.EnzymeReaction(RA[RA(S,[:C=>1]),RA(P,[:C=>1])], ER.RegulatorMults[], Int[2])
    steps=Vector{St}[[St(E,ES,S,true)],[St(ES,EP,nothing,false)],[St(EP,E,P,true)]]
    ER.AllostericMechanism(rxn, steps, collect(Symbol,states), 2, ER.RegulatorySite[])
end
function uni_ss(states)   # all-steady-state uni-uni (bindings carry kon/koff)
    E=Sp(Met[],:E); ES=Sp(Met[S],:E); EP=Sp(Met[P],:E)
    rxn=ER.EnzymeReaction(RA[RA(S,[:C=>1]),RA(P,[:C=>1])], ER.RegulatorMults[], Int[2])
    steps=Vector{St}[[St(E,ES,S,false)],[St(ES,EP,nothing,false)],[St(EP,E,P,false)]]
    ER.AllostericMechanism(rxn, steps, collect(Symbol,states), 2, ER.RegulatorySite[])
end
function evalrate(am; seed=1, split=nothing)
    cem=ER.compile_mechanism(am); fp=ER.fitted_params(am); rng=MersenneTwister(seed)
    mets=collect(ER.metabolites(cem))
    base=[(k===:L ? 0.6 : 0.4+2rand(rng)) for k in fp]
    split!==nothing && (base=[(fp[i]===split[1] ? split[2] : base[i]) for i in 1:length(fp)])
    prm=NamedTuple{(fp...,:Keq,:E_total)}((base...,3.0,1.0))
    c=NamedTuple{Tuple(mets)}(ntuple(i->0.4+2rand(rng),length(mets)))
    v=real(ER.rate_equation(cem,c,prm))
    ec=NamedTuple{Tuple(mets)}(ntuple(i->(mets[i]===:P ? 3.0 : 1.0),length(mets)))
    veq=real(ER.rate_equation(cem,ec,prm))
    (fp, v, veq)
end

@testset "strict :EqualAI collapse" begin
    fp0,_,_ = evalrate(uni([:EqualAI,:EqualAI,:EqualAI]))     # baseline

    @testset "single NonequalAI binding + EqualAI catalysis -> full collapse" begin
        fp,v,veq = evalrate(uni([:NonequalAI,:EqualAI,:EqualAI]))
        @test isfinite(v); @test abs(veq) < 1e-8
        @test !(:K_I_S_E in fp)                 # I-twin dropped (collapsed to a mirror)
        s = ER.rate_equation_string(uni([:NonequalAI,:EqualAI,:EqualAI]))
        @test occursin("K_I_S_E=K_A_S_E", replace(s," "=>""))  # explicit mirror
        @test !occursin("k_I_", s)              # catalysis not silently un-shared
    end

    @testset "two NonequalAI bindings + EqualAI catalysis -> 1 honorable DOF" begin
        am = uni([:NonequalAI,:EqualAI,:NonequalAI])
        fp,v,veq = evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        # exactly one free I-split survives; the other is a derived mirror.
        nI = count(p->startswith(String(p),"K_I_"), fp)
        @test nI == 1
        s = replace(ER.rate_equation_string(am), " "=>"")
        @test occursin("k_I_", s) == false                    # catalysis stays shared
        # the surviving split moves the rate (identifiable)
        freeI = fp[findfirst(p->startswith(String(p),"K_I_"), fp)]
        v1 = evalrate(am; split=(freeI,1.3))[2]; v2 = evalrate(am; split=(freeI,5.0))[2]
        @test !isapprox(v1, v2)
    end

    @testset "catalysis NonequalAI -> native, no collapse, no mirror" begin
        am = uni([:NonequalAI,:NonequalAI,:EqualAI])
        fp,v,veq = evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test :K_I_S_E in fp                    # binding split free
        @test any(p->startswith(String(p),"k_I_"), fp)  # catalysis split free (native)
    end

    @testset "binding-Wegscheider single inner edge -> full collapse" begin
        # random-order bi-bi (two Wegscheider boxes); tag ONE inner box-independent
        # edge (EB+A->EAB) :NonequalAI, rest :EqualAI -> its split is forbidden -> collapses.
        A2=Sub(:A); B2=Sub(:B); Q2=Prd(:Q)
        E=Sp(Met[],:E); EA=Sp(Met[A2],:E); EB=Sp(Met[B2],:E); EAB=Sp(Met[A2,B2],:E)
        EPQ=Sp(Met[P,Q2],:E); EP=Sp(Met[P],:E); EQ=Sp(Met[Q2],:E)
        rxn=ER.EnzymeReaction(RA[RA(A2,[:C=>1]),RA(B2,[:N=>1]),RA(P,[:C=>1]),RA(Q2,[:N=>1])],
                              ER.RegulatorMults[], Int[2])
        sd=[St(E,EA,A2,true),St(E,EB,B2,true),St(EB,EAB,A2,true),St(EA,EAB,B2,true),
            St(EAB,EPQ,nothing,false),St(EP,EPQ,Q2,true),St(EQ,EPQ,P,true),
            St(E,EP,P,true),St(E,EQ,Q2,true)]
        st=fill(:EqualAI,9); st[3]=:NonequalAI
        am=ER.AllostericMechanism(rxn, Vector{St}[[s] for s in sd], st, 2, ER.RegulatorySite[])
        cem=ER.compile_mechanism(am); fp=ER.fitted_params(am)
        @test !(:K_I_A_EB in fp)                            # forbidden split collapsed
        s=replace(ER.rate_equation_string(am)," "=>"")
        @test occursin("K_I_A_EB=K_A_A_EB", s)              # explicit mirror
        rng=MersenneTwister(2)
        pv=Tuple((k===:L ? 0.6 : 0.4+2rand(rng)) for k in fp)
        prm=NamedTuple{(fp...,:Keq,:E_total)}((pv...,3.0,1.0))
        mets=collect(ER.metabolites(cem))
        ec=NamedTuple{Tuple(mets)}(ntuple(i->(mets[i] in (:A,:B) ? 1.0 : sqrt(3.0)),length(mets)))
        @test abs(real(ER.rate_equation(cem,ec,prm))) < 1e-8   # thermo-consistent
    end

    # ── Steady-state binding: affinity/speed decomposition (Option 3) ──
    # A steady-state :NonequalAI binding carries two rate constants: its affinity
    # (kon/koff) can be forbidden and collapse (deriving the reverse), while its
    # speed (the forward kon) stays free. Only the affinity collapses — never the
    # whole binding. (The steady-state Wegscheider-box case is covered by `m_ro`
    # in test_rate_eq_derivation.jl.)
    @testset "SS binding + EqualAI catalysis -> affinity collapses, speed free" begin
        am = uni_ss([:NonequalAI,:EqualAI,:EqualAI])
        fp,v,veq = evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test !(:koff_I_S_E in fp)              # affinity collapsed: reverse derived
        @test :kon_I_S_E in fp                  # speed (forward) stays free
        s = replace(ER.rate_equation_string(am), " "=>"")
        @test occursin("koff_I_S_E=", s)        # explicit reverse-rate mirror
        # the surviving speed split moves the rate (identifiable)
        v1 = evalrate(am; split=(:kon_I_S_E,1.3))[2]
        v2 = evalrate(am; split=(:kon_I_S_E,5.0))[2]
        @test !isapprox(v1, v2)
    end

    @testset "SS binding + NonequalAI catalysis -> not forbidden, stays free" begin
        am = uni_ss([:NonequalAI,:NonequalAI,:EqualAI])
        fp,v,veq = evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test (:kon_I_S_E in fp) && (:koff_I_S_E in fp)   # both free (affinity honorable)
    end

    # ── Mixed-type + Wegscheider-pivot coupling (regressions for the two review
    #    Criticals: uniform affinity sign, and indep_A-keyed collapsibility + the
    #    transitive S_I closure). Both were passing the suite while silently broken.
    @testset "mixed RE/SS coupled :NonequalAI bindings -> consistent (uniform sign)" begin
        # S-binding RE, P-release SS, both :NonequalAI, catalysis :EqualAI: the two
        # coupled affinities are *different* step types, so the constraint matrix must
        # use one uniform sign — a per-type flip inverts the coupling → nonzero flux.
        E=Sp(Met[],:E); ES=Sp(Met[S],:E); EP=Sp(Met[P],:E)
        rxn=ER.EnzymeReaction(RA[RA(S,[:C=>1]),RA(P,[:C=>1])], ER.RegulatorMults[], Int[2])
        steps=Vector{St}[[St(E,ES,S,true)],[St(ES,EP,nothing,false)],[St(EP,E,P,false)]]
        am=ER.AllostericMechanism(rxn, steps, [:NonequalAI,:EqualAI,:NonequalAI], 2,
                                  ER.RegulatorySite[])
        fp,v,veq = evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
    end

    @testset ":NonequalAI RE Wegscheider-pivot binding -> absorbed, no undefined symbol" begin
        # random-order bi-bi, all-RE bindings, SS :EqualAI catalysis; tag two box
        # edges :NonequalAI so one is the Wegscheider pivot (its Kd is derived). The
        # pivot must be absorbed (not collapsed, else a circular mirror ⇒ UndefVarError),
        # and the free split's I-symbol — reachable only through the other edge's
        # collapse mirror — must be retained (transitive S_I closure), not left undefined.
        A2=Sub(:A); B2=Sub(:B); Q2=Prd(:Q)
        E=Sp(Met[],:E); EA=Sp(Met[A2],:E); EB=Sp(Met[B2],:E); EAB=Sp(Met[A2,B2],:E)
        EPQ=Sp(Met[P,Q2],:E); EP=Sp(Met[P],:E); EQ=Sp(Met[Q2],:E)
        rxn=ER.EnzymeReaction(RA[RA(A2,[:C=>1]),RA(B2,[:N=>1]),RA(P,[:C=>1]),RA(Q2,[:N=>1])],
                              ER.RegulatorMults[], Int[2])
        sd=[St(E,EA,A2,true),St(E,EB,B2,true),St(EB,EAB,A2,true),St(EA,EAB,B2,true),
            St(EAB,EPQ,nothing,false),St(EP,EPQ,Q2,true),St(EQ,EPQ,P,true),
            St(E,EP,P,true),St(E,EQ,Q2,true)]
        st=fill(:EqualAI,9); st[2]=:NonequalAI; st[3]=:NonequalAI
        am=ER.AllostericMechanism(rxn, Vector{St}[[s] for s in sd], st, 2, ER.RegulatorySite[])
        cem=ER.compile_mechanism(am)
        fp=collect(ER.fitted_params(am)); rng=MersenneTwister(3)
        prm=NamedTuple{(fp...,:Keq,:E_total)}(
            (ntuple(i->(fp[i]===:L ? 0.6 : 0.4+2rand(rng)),length(fp))...,4.0,1.0))
        mets=collect(ER.metabolites(cem))
        ec=NamedTuple{Tuple(mets)}(ntuple(i->(mets[i] in (:A,:B) ? 1.0 : 2.0),length(mets)))
        v=real(ER.rate_equation(cem,ec,prm))         # must not throw UndefVarError
        @test isfinite(v); @test abs(v) < 1e-8
    end

    @testset "dead-I NonequalAI binding -> K_I identifiable, NOT collapsed" begin
        # I state cannot turn over (OnlyA catalysis) but binds S with its own
        # affinity: K_A_S_E and K_I_S_E are BOTH identifiable (a dead-end E_I·S is
        # in no cycle, so nothing pins K_I to K_A). HEAD over-collapses this.
        am = uni([:NonequalAI, :OnlyA, :EqualAI])
        fp, v, veq = evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test :K_I_S_E in fp                          # NOT collapsed
        v1 = evalrate(am; split=(:K_I_S_E, 1.3))[2]
        v2 = evalrate(am; split=(:K_I_S_E, 5.0))[2]
        @test !isapprox(v1, v2)                       # identifiable (moves the rate)
        s = replace(ER.rate_equation_string(am), " " => "")
        @test !occursin("K_I_S_E=K_A_S_E", s)         # no collapse mirror
    end
end
end # module
