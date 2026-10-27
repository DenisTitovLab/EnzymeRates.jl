# ABOUTME: Allosteric :NonequalAI split regression — a split with no thermodynamic freedom
# ABOUTME: collapses to K_I=K_A; an identifiable split stays free; all equations thermo-consistent.
module AllostericCollapseTests
using Test, EnzymeRates, Random
const ER = EnzymeRates
const M=ER.Mechanisms.Met; const RS=ER.Mechanisms.RawStep
# A form in conformation E that holds the substrates and products `xs`.
_testhelper_form(xs...) = ER.Mechanisms.RawForm(:E, M[M(x, false) for x in xs])
# The allosteric mechanism over substrates `subs` and products `prods`, catalytic
# multiplicity 2, whose steps are `steps`, step `j` alone in group `j`, tagged `states[j]`.
_testhelper_one_step_groups(subs, prods, steps, states) =
    ER.AllostericMechanism(ER.Mechanisms.RawAllostericMechanism(
        ER.Mechanisms.RawMechanism(subs, prods, Symbol[], steps, collect(eachindex(steps))),
        collect(Symbol, states), 2, ER.Mechanisms.RawSite[]))
# `re` flags which of the S-binding, chemistry and P-release steps are rapid
# equilibrium; a steady-state binding carries kon/koff.
function _testhelper_uni(states; re=(true,false,true))
    E=_testhelper_form(); ES=_testhelper_form(:S); EP=_testhelper_form(:P)
    steps=[RS(E,ES,[M(:S,false)],M[],re[1]),RS(ES,EP,M[],M[],re[2]),
           RS(EP,E,M[],[M(:P,false)],re[3])]
    _testhelper_one_step_groups([:S], [:P], steps, states)
end
# Random-order bi-bi with two Wegscheider boxes; `states` tags its nine one-step groups.
function _testhelper_ro_bibi(states)
    E=_testhelper_form(); EA=_testhelper_form(:A); EB=_testhelper_form(:B)
    EAB=_testhelper_form(:A,:B); EPQ=_testhelper_form(:P,:Q); EP=_testhelper_form(:P)
    EQ=_testhelper_form(:Q)
    A=M(:A,false); B=M(:B,false); P=M(:P,false); Q=M(:Q,false)
    sd=[RS(E,EA,[A],M[],true),RS(E,EB,[B],M[],true),
        RS(EB,EAB,[A],M[],true),RS(EA,EAB,[B],M[],true),
        RS(EAB,EPQ,M[],M[],false),RS(EP,EPQ,[Q],M[],true),
        RS(EQ,EPQ,[P],M[],true),RS(E,EP,[P],M[],true),
        RS(E,EQ,[Q],M[],true)]
    _testhelper_one_step_groups([:A, :B], [:P, :Q], sd, states)
end
function _testhelper_evalrate(am; seed=1, split=nothing)
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
    # baseline
    fp0,_,_ = _testhelper_evalrate(_testhelper_uni([:EqualAI,:EqualAI,:EqualAI]))

    @testset "single NonequalAI binding + EqualAI catalysis -> full collapse" begin
        fp,v,veq = _testhelper_evalrate(_testhelper_uni([:NonequalAI,:EqualAI,:EqualAI]))
        @test isfinite(v); @test abs(veq) < 1e-8
        @test !(:K_I_ES_to_E_S in fp)           # I-twin dropped (collapsed to a mirror)
        s = ER.rate_equation_string(_testhelper_uni([:NonequalAI,:EqualAI,:EqualAI]))
        @test occursin("K_I_ES_to_E_S=K_A_ES_to_E_S", replace(s," "=>""))  # explicit mirror
        @test !occursin("k_I_", s)              # catalysis not silently un-shared
    end

    @testset "two NonequalAI bindings + EqualAI catalysis -> 1 honorable DOF" begin
        am = _testhelper_uni([:NonequalAI,:EqualAI,:NonequalAI])
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        # exactly one free I-split survives; the other is a derived mirror.
        nI = count(p->startswith(String(p),"K_I_"), fp)
        @test nI == 1
        s = replace(ER.rate_equation_string(am), " "=>"")
        @test occursin("k_I_", s) == false                    # catalysis stays shared
        # the surviving split moves the rate (identifiable)
        freeI = fp[findfirst(p->startswith(String(p),"K_I_"), fp)]
        v1 = _testhelper_evalrate(am; split=(freeI,1.3))[2]
        v2 = _testhelper_evalrate(am; split=(freeI,5.0))[2]
        @test !isapprox(v1, v2)
    end

    @testset "catalysis NonequalAI -> native, no collapse, no mirror" begin
        am = _testhelper_uni([:NonequalAI,:NonequalAI,:EqualAI])
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test :K_I_ES_to_E_S in fp              # binding split free
        @test any(p->startswith(String(p),"k_I_"), fp)  # catalysis split free (native)
    end

    @testset "binding-Wegscheider single inner edge -> full collapse" begin
        # random-order bi-bi (two Wegscheider boxes); tag ONE inner box-independent
        # edge (EB+A->EAB) :NonequalAI, rest :EqualAI -> its split is forbidden -> collapses.
        st=fill(:EqualAI,9); st[3]=:NonequalAI
        am=_testhelper_ro_bibi(st)
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8             # thermo-consistent
        @test !(:K_I_EAB_to_EB_A in fp)                     # forbidden split collapsed
        dep,_=Main._testhelper_dependent_param_exprs(am)
        @test dep[:K_I_EAB_to_EB_A] == dep[:K_A_EAB_to_EB_A]  # one value in both states
    end

    # ── Steady-state binding: affinity/speed decomposition (Option 3) ──
    # A steady-state :NonequalAI binding carries two rate constants: its affinity
    # (kon/koff) can be forbidden and collapse (deriving the reverse), while its
    # speed (the forward kon) stays free. Only the affinity collapses — never the
    # whole binding. (The steady-state Wegscheider-box case is covered by `m_ro`
    # in test_rate_eq_derivation.jl.)
    @testset "SS binding + EqualAI catalysis -> affinity collapses, speed free" begin
        am = _testhelper_uni([:NonequalAI,:EqualAI,:EqualAI]; re=(false,false,false))
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test !(:k_I_ES_to_E_S in fp)           # affinity collapsed: reverse derived
        @test :k_I_E_S_to_ES in fp              # speed (forward) stays free
        s = replace(ER.rate_equation_string(am), " "=>"")
        @test occursin("k_I_ES_to_E_S=", s)     # explicit reverse-rate mirror
        # the surviving speed split moves the rate (identifiable)
        v1 = _testhelper_evalrate(am; split=(:k_I_E_S_to_ES,1.3))[2]
        v2 = _testhelper_evalrate(am; split=(:k_I_E_S_to_ES,5.0))[2]
        @test !isapprox(v1, v2)
    end

    @testset "SS binding + NonequalAI catalysis -> not forbidden, stays free" begin
        am = _testhelper_uni([:NonequalAI,:NonequalAI,:EqualAI]; re=(false,false,false))
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        # both free (affinity honorable)
        @test (:k_I_E_S_to_ES in fp) && (:k_I_ES_to_E_S in fp)
    end

    # ── Mixed-type + Wegscheider-coupled :NonequalAI bindings. Two mechanisms that
    #    were silently mis-derived before the combined solve: the constraint matrix
    #    must use one uniform binding-K sign, and an all-RE Wegscheider box with
    #    :NonequalAI edges must derive its tied binding-K with no undefined symbol.
    @testset "mixed RE/SS coupled :NonequalAI bindings -> consistent (uniform sign)" begin
        # S-binding RE, P-release SS, both :NonequalAI, catalysis :EqualAI: the two
        # coupled affinities are *different* step types, so the constraint matrix must
        # use one uniform sign — a per-type flip inverts the coupling → nonzero flux.
        am=_testhelper_uni([:NonequalAI,:EqualAI,:NonequalAI]; re=(true,false,false))
        fp,v,veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
    end

    @testset ":NonequalAI RE Wegscheider box binding -> consistent, no undefined symbol" begin
        # random-order bi-bi, all-RE bindings, SS :EqualAI catalysis; tag two box
        # edges :NonequalAI so one binding-K is Wegscheider-derived from the other.
        # The combined solve must express every derived symbol purely in free columns
        # (no circular reference ⇒ no UndefVarError) and hold detailed balance.
        st=fill(:EqualAI,9); st[2]=:NonequalAI; st[3]=:NonequalAI
        # must not throw UndefVarError
        fp,v,veq = _testhelper_evalrate(_testhelper_ro_bibi(st))
        @test isfinite(v); @test abs(veq) < 1e-8
    end

    @testset "dead-I NonequalAI binding -> K_I identifiable, NOT collapsed" begin
        # I state cannot turn over (OnlyA catalysis) but binds S with its own
        # affinity: K_A_ES_to_E_S and K_I_ES_to_E_S are BOTH identifiable (a dead-end
        # E_I·S is in no cycle, so nothing pins K_I to K_A). HEAD over-collapses this.
        am = _testhelper_uni([:NonequalAI, :OnlyA, :EqualAI])
        fp, v, veq = _testhelper_evalrate(am)
        @test isfinite(v); @test abs(veq) < 1e-8
        @test :K_I_ES_to_E_S in fp                    # NOT collapsed
        v1 = _testhelper_evalrate(am; split=(:K_I_ES_to_E_S, 1.3))[2]
        v2 = _testhelper_evalrate(am; split=(:K_I_ES_to_E_S, 5.0))[2]
        @test !isapprox(v1, v2)                       # identifiable (moves the rate)
        s = replace(ER.rate_equation_string(am), " " => "")
        @test !occursin("K_I_ES_to_E_S=K_A_ES_to_E_S", s)  # no collapse mirror
    end
end
end # module
