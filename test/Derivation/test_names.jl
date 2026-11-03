# ABOUTME: Tests of the parameter names: the name(p, m) chokepoint, the constant columns and
# ABOUTME: the regulator constants.

const ER = EnzymeRates

@testset "constant_columns: the joint solve's columns and their names" begin
    M = ER.Mechanisms
    D = ER.Derivation
    col = D.StepParameter
    # Uni-uni: E + P → E(P) and E + S → E(S) at rapid equilibrium, E(S) → E(P) steady
    # state, each its own group. The isomerization's columns come first (kind 3), its
    # reverse rate, which runs against the reaction, before its forward one; then the
    # bindings, the larger content rank first.
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    c = D.constant_columns(m)
    @test c.columns == [col(3, 2, false), col(3, 1, false), col(2, 1, false),
                        col(1, 1, false)]
    @test c.names == [:k_EP_to_ES, :k_ES_to_EP, :K_ES_to_E_S, :K_EP_to_E_P]
    @test c.in_law == trues(4)
    @test (c.flat, c.tags, c.law) == ([1, 2, 3], [:None, :None, :None], falses(3))

    # Allosteric uni-uni: E + P → E(P) :EqualAI, E + S → E(S) :NonequalAI and
    # E(S) → E(P) :OnlyA, one group each. The :NonequalAI and :OnlyA groups get inactive
    # columns; the inactive copies of the :OnlyA step lie outside the inactive law
    # (block 2) and come first, then the inactive binding (block 1), then the active
    # columns (block 0). An :OnlyA step's inactive copy takes its own name.
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: OnlyA
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end)
    a = D.constant_columns(am)
    @test a.columns == [col(3, 2, true), col(3, 1, true), col(2, 1, true),
                        col(3, 2, false), col(3, 1, false), col(2, 1, false),
                        col(1, 1, false)]
    @test a.names == [:k_I_EP_to_ES, :k_I_ES_to_EP, :K_I_ES_to_E_S, :k_A_EP_to_ES,
                      :k_A_ES_to_EP, :K_A_ES_to_E_S, :K_EP_to_E_P]
    @test a.in_law == BitVector([0, 0, 1, 1, 1, 1, 1])
    @test (a.flat, a.tags, a.law) ==
          ([1, 2, 3], [:EqualAI, :NonequalAI, :OnlyA], BitVector([1, 1, 0]))
end

@testset "regulator_constant: a regulator ligand's constant per conformation" begin
    M = ER.Mechanisms
    D = ER.Derivation
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R, X, Y, Z
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: NonequalAI
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyI
        regulatory_site(2): X::NonequalAI
        regulatory_site(2): Y::EqualAI
        regulatory_site(2): Z::OnlyA
    end)
    # Each ligand with its own tag, in the active and then the inactive conformation: an
    # :OnlyI ligand is absent from the active one and an :OnlyA ligand from the inactive
    # one; an :EqualAI ligand has one bare name in both.
    @test [D.regulator_constant(am, lig, tag, inactive)
           for (lig, tag) in ((:R, :OnlyI), (:X, :NonequalAI), (:Y, :EqualAI), (:Z, :OnlyA))
           for inactive in (false, true)] ==
          [nothing, :K_I_Rreg, :K_A_Xreg, :K_I_Xreg, :K_Yreg, :K_Yreg, :K_A_Zreg, nothing]
end

@testset "Parameter family: step, regulator and allosteric parameters" begin
    # A step parameter serves every kind of edge: a slot of an edge of a mechanism's
    # catalytic graph in one conformation.
    P = ER.Derivation.StepParameter
    p = P(1, 1, false)
    @test p isa ER.Derivation.Parameter
    @test (p.step, p.slot, p.inactive) == (1, 1, false)
    @test p == P(1, 1, false)
    @test p != P(2, 1, false) && p != P(1, 2, false) && p != P(1, 1, true)

    # A regulator constant is its ligand's in one conformation; it needs no site, since a
    # ligand sits at one site.
    R = ER.Derivation.RegulatorParameter
    kr = R(:A, false)
    @test kr isa ER.Derivation.Parameter
    @test (kr.ligand, kr.inactive) == (:A, false)
    @test kr == R(:A, false)
    @test kr != R(:A, true) && kr != R(:B, false)

    @test ER.Derivation.AllostericParameter() isa ER.Derivation.Parameter
    @test ER.Derivation.AllostericParameter() == ER.Derivation.AllostericParameter()
end

@testset "parameter_info lists the fitted parameters in parameters() order" begin
    for spec in MECHANISM_TEST_SPECS
        m = _testhelper_lift(spec.mechanism)
        info = ER.Derivation.parameter_info(m)
        @test [ER.Derivation.name(p, m) for p in info] ==
              collect(parameters(m, FittedParameters()))
        @test all(p -> p isa ER.Derivation.Parameter, info)
        empty!(info)
        @test !isempty(ER.Derivation.parameter_info(m))
    end
end

@testset "parameter_info points at the steps of Uni-Uni" begin
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
        end
    end)
    P = ER.Derivation.StepParameter
    @test ER.Derivation.parameter_info(m) ==
          [P(1, 2, false), P(2, 2, false), P(2, 1, false)]
    @test parameters(m, FittedParameters()) == (:k_ES_to_E_P, :k_ES_to_E_S, :k_E_S_to_ES)
end

@testset "parameter_info ends an allosteric mechanism with its regulators and L" begin
    m = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: EqualAI
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end)
    info = ER.Derivation.parameter_info(m)
    @test info[(end - 2):end] == [ER.Derivation.RegulatorParameter(:R, false),
                                  ER.Derivation.RegulatorParameter(:R, true),
                                  ER.Derivation.AllostericParameter()]
    @test [ER.Derivation.name(p, m) for p in info[(end - 2):end]] ==
          [:K_A_Rreg, :K_I_Rreg, :L]
    @test [ER.Derivation.name(p, m) for p in info] ==
          collect(parameters(m, FittedParameters()))
end

@testset "name(p, ::AllostericMechanism) chokepoint" begin
    aem = @allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: NonequalAI
            E(P) ⇌ E + P      :: NonequalAI
        end
        regulatory_site(1): R::NonequalAI
    end
    am = ER.AllostericMechanism(aem)

    # Edge parameters, edges of the catalytic graph picked by content.
    P = ER.Derivation.StepParameter
    rep_bind = _testhelper_edge(am, "E + S ⇌ ES")
    @test ER.Derivation.name(P(rep_bind, 1, false), am) === :K_ES_to_E_S
    @test ER.Derivation.name(P(rep_bind, 1, true), am) === :K_I_ES_to_E_S

    # An :EqualAI isomerization, steady-state and at rapid equilibrium, has no state tag.
    ss_iso = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: EqualAI
            E(P) ⇌ E + P      :: NonequalAI
        end
        regulatory_site(1): R::NonequalAI
    end)
    rep_iso = _testhelper_edge(ss_iso, "ES <--> EP")
    @test ER.Derivation.name(P(rep_iso, 1, false), ss_iso) === :k_ES_to_EP
    @test ER.Derivation.name(P(rep_iso, 2, false), ss_iso) === :k_EP_to_ES
    re_iso = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S <--> E(S)   :: NonequalAI
            E(S) ⇌ E(P)       :: EqualAI
            E(P) <--> E + P   :: NonequalAI
        end
        regulatory_site(1): R::NonequalAI
    end)
    @test ER.Derivation.name(P(_testhelper_edge(re_iso, "ES ⇌ EP"), 1, false), re_iso) ===
          :K_ES_to_EP

    lig = first(am.regulatory_sites[1].ligands)
    @test ER.Derivation.name(ER.Derivation.RegulatorParameter(lig, false), am) === :K_A_Rreg
    @test ER.Derivation.name(ER.Derivation.RegulatorParameter(lig, true), am) === :K_I_Rreg
    @test ER.Derivation.name(ER.Derivation.AllostericParameter(), am) === :L

    # The chokepoint renders on the concrete mechanism only; a compiled type
    # is lifted first.
    @test_throws MethodError ER.Derivation.name(P(rep_bind, 1, false), aem)
    @test_throws MethodError ER.Derivation.name(
        ER.Derivation.RegulatorParameter(lig, false), aem)
    @test_throws MethodError ER.Derivation.name(ER.Derivation.AllostericParameter(), aem)
end

@testset "name(p::Parameter, m) chokepoint" begin
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    # The binding of S, the isomerization and the binding of P, as edges of `m`.
    bind, iso, rel = (_testhelper_edge(m, text)
                      for text in ("E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"))

    P = ER.Derivation.StepParameter

    # Structural naming: every step constant encodes its reaction's two sides;
    # a binding K reads in the release direction, iso params in the stored one.
    @test ER.Derivation.name(P(bind, 1, false), m) === :K_ES_to_E_S
    @test ER.Derivation.name(P(iso, 1, false), m) === :k_ES_to_EP
    @test ER.Derivation.name(P(iso, 2, false), m) === :k_EP_to_ES
    @test ER.Derivation.name(P(rel, 1, false), m) === :K_EP_to_E_P

    # A binding's rate pair: the forward rate binds, the reverse rate releases
    ss_bind = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)
    k = _testhelper_edge(ss_bind, "E + S <--> ES")
    @test ER.Derivation.name(P(k, 1, false), ss_bind) === :k_E_S_to_ES
    @test ER.Derivation.name(P(k, 2, false), ss_bind) === :k_ES_to_E_S

    # I-state token on RE and SS steps
    am_ss = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: NonequalAI
            E(S) <--> E(P)    :: NonequalAI
            E(P) ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end)
    k = _testhelper_edge(am_ss, "E + S ⇌ ES")
    @test ER.Derivation.name(P(k, 1, true), am_ss) === :K_I_ES_to_E_S
    k = _testhelper_edge(am_ss, "ES <--> EP")
    @test ER.Derivation.name(P(k, 1, true), am_ss) === :k_I_ES_to_EP
    @test ER.Derivation.name(P(k, 2, true), am_ss) === :k_I_EP_to_ES

    # The RE iso's equilibrium constant, named in the stored direction
    re_iso = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> E(S)
            E(S) ⇌ E(P)
            E(P) <--> E + P
        end
    end)
    @test ER.Derivation.name(P(_testhelper_edge(re_iso, "ES ⇌ EP"), 1, false), re_iso) ===
          :K_ES_to_EP
    am_re = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products:   P
        allosteric_regulators: R
        catalytic_multiplicity: 2
        steps: begin
            E + S <--> E(S)   :: EqualAI
            E(S) ⇌ E(P)       :: NonequalAI
            E(P) <--> E + P   :: EqualAI
        end
        regulatory_site(2): R::NonequalAI
    end)
    @test ER.Derivation.name(P(_testhelper_edge(am_re, "ES ⇌ EP"), 1, true), am_re) ===
          :K_I_ES_to_EP

    # Same names resolve on the mechanism lifted back from EnzymeMechanism(m);
    # the compiled type itself is not a chokepoint argument.
    em = EnzymeMechanism(m)
    @test ER.Derivation.name(P(bind, 1, false), ER.Mechanism(em)) === :K_ES_to_E_S
    @test ER.Derivation.name(P(iso, 1, false), ER.Mechanism(em)) === :k_ES_to_EP
    @test_throws MethodError ER.Derivation.name(P(bind, 1, false), em)
end

@testset "name(p::Parameter, m) for the steps of a shared kinetic group" begin
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            (E + S ⇌ E(S), E(P) + S ⇌ E(P, S))
            E(S) <--> E(P)
            E(P) ⇌ E + P
        end
    end)

    # Both steps bind S; rep = bind. Both yield the same name.
    @test ER.Derivation.name(
        ER.Derivation.StepParameter(_testhelper_edge(m, "E + S ⇌ ES"), 1, false), m) ===
          :K_ES_to_E_S
    @test ER.Derivation.name(
        ER.Derivation.StepParameter(_testhelper_edge(m, "EP + S ⇌ EPS"), 1, false), m) ===
          :K_ES_to_E_S
end

@testset "fused steps are named by their sides" begin
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(Q) + P
            E(Q) <--> E + Q
        end
    end
    names = _testhelper_step_constant_names(m)
    @test :k_EAB_to_EQ_P in names && :k_EQ_P_to_EAB in names
    # The fused release E(A, B) → E(Q) + P is stored as the binding of P it
    # reverses, E(Q) + P → E(A, B), so its forward rate constant binds P.
    mech = ER.Mechanism(m)
    k = only(k for k in eachindex(mech.steps)
             if mech.step_kinds[k] !== :binding)
    kon, koff = _testhelper_step_parameters(mech, k)
    @test (kon.slot, koff.slot) == (1, 2)
    @test ER.Derivation.name(kon, mech) == :k_EQ_P_to_EAB &&
          ER.Derivation.name(koff, mech) == :k_EAB_to_EQ_P
    @test :k_E_A_to_EA in names && :k_EQ_to_E_Q in names
end

@testset "step constants are named by their reaction" begin
    # Michaelis–Menten, every step steady-state: each binding has one rate constant
    # per direction, named after the reaction it drives; the isomerization has k both ways.
    mm_ss = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E(P)
            E(P) <--> E + P
        end
    end
    @test issubset([:k_E_S_to_ES, :k_ES_to_E_S, :k_E_P_to_EP, :k_EP_to_E_P],
                   _testhelper_step_constant_names(mm_ss))
    # The Haldane reduction leaves the SS isomerization constant fitted.
    @test :k_ES_to_EP in ER.parameters(mm_ss)
    # An RE isomerization has one equilibrium constant, in its canonical direction.
    re_iso = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) ⇌ E(P)
            E(P) <--> E + P
        end
    end
    @test :K_ES_to_EP in _testhelper_step_constant_names(re_iso)
    # An RE fused release is stored as the binding of P it reverses, so its
    # constant is a dissociation constant, [EQ]·[P]/[EAB], named in the release
    # direction.
    re_fused = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) ⇌ E(Q) + P
            E(Q) <--> E + Q
        end
    end
    @test :K_EAB_to_EQ_P in _testhelper_step_constant_names(re_fused)
    # An RE fused binding takes a dissociation constant too, [E]·[S]/[EP], named in
    # the release direction like a plain binding's.
    re_fused_binding = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(P)
            E(P) <--> E + P
        end
    end
    @test :K_EP_to_E_S in _testhelper_step_constant_names(re_fused_binding)
    # A competitive-inhibitor copy of A and A itself bind E without colliding.
    inh = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: A
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P)
            E(P) ⇌ E + P
            E + A::Mod ⇌ E(A::Mod)
        end
    end
    @test issubset([:K_EAmod_to_E_Amod, :K_EA_to_E_A],
                   _testhelper_step_constant_names(inh))
    _testhelper_check_law_against_oracle(inh)
    # The allosteric state tag follows the prefix: a :NonequalAI binding has an
    # A-state and an I-state constant, an :EqualAI binding one shared constant.
    allo = @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)     :: NonequalAI
            E(S) <--> E(P)   :: EqualAI
            E(P) ⇌ E + P     :: EqualAI
        end
    end
    # Every constant the rate equation names, fitted or dependent.
    dep, indep = _testhelper_dependent_param_exprs(typeof(allo))
    allo_names = union(keys(dep), indep)
    @test issubset([:K_A_ES_to_E_S, :K_I_ES_to_E_S], allo_names)
    @test !(:K_ES_to_E_S in allo_names)
    @test :K_EP_to_E_P in ER.parameters(allo) && !(:K_A_EP_to_E_P in allo_names)
    @test :K_I_EP_to_E_P ∉ ER.parameters(allo)
    # A ping-pong residual form: B binds E(; residual = A - P), whose name is
    # E_res_plusA_minusP, so the release-direction K reads
    # EB_res_plusA_minusP → E_res_plusA_minusP + B.
    # E and the covalent E(; residual) share one RE segment whose weights vanish only
    # at B = Q = 0 (mixed), so the mechanism is accepted.
    pingpong = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P; residual = A - P)
            E(; residual = A - P) + P ⇌ E(P; residual = A - P)
            E(; residual = A - P) + B ⇌ E(B; residual = A - P)
            E(B; residual = A - P) ⇌ E(Q)
            E + Q ⇌ E(Q)
        end
    end
    pp_names = _testhelper_step_constant_names(pingpong)
    @test :K_EB_res_plusA_minusP_to_E_res_plusA_minusP_B in pp_names
    _testhelper_check_law_against_oracle(pingpong)
    # One name per class: a rapid-equilibrium group's K, a steady-state group's kf, kr.
    pm = ER.Mechanism(pingpong)
    @test length(pp_names) == sum(pm.steps[first(ks)].rapid_equilibrium ? 1 : 2
                                  for ks in pm.groups)
end

@testset "a modifier-only name is bare in constant names, a reactant's copy is Xmod" begin
    # I binds only as a modifier, so nothing else carries its name.
    modifier_only = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
            E + I ⇌ E(I)
        end
    end
    @test _testhelper_step_constant_names(modifier_only) ==
          Set([:K_ES_to_E_S, :k_ES_to_EP, :k_EP_to_ES, :K_EP_to_E_P, :K_EI_to_E_I])
    _testhelper_check_law_against_oracle(modifier_only)
    # P is a product that also binds as a modifier, so its modifier copy reads Pmod.
    reactant_copy = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: P
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q
            E(Q) + P::Mod ⇌ E(P::Mod, Q)
        end
    end
    @test _testhelper_step_constant_names(reactant_copy) ==
          Set([:K_EA_to_E_A, :K_EAB_to_EA_B, :k_EAB_to_EPQ, :k_EPQ_to_EAB,
               :K_EPQ_to_EQ_P, :K_EQ_to_E_Q, :K_EPmodQ_to_EQ_Pmod])
    _testhelper_check_law_against_oracle(reactant_copy)
end

@testset "a residual reads _res, then _plus per added and _minus per subtracted name" begin
    # The residual adds A and B and subtracts P; each run lists its names in order.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> E(A)
            E(A) + B <--> E(A, B)
            E(A, B) <--> E(P; residual = B + A - P)
            E(P; residual = B + A - P) <--> E(; residual = B + A - P) + P
            E(; residual = B + A - P) <--> E(Q)
            E(Q) <--> E + Q
        end
    end
    @test _testhelper_step_constant_names(m) == Set([
        :k_E_A_to_EA, :k_EA_to_E_A, :k_EA_B_to_EAB, :k_EAB_to_EA_B,
        :k_EAB_to_EP_res_plusA_plusB_minusP, :k_EP_res_plusA_plusB_minusP_to_EAB,
        :k_E_res_plusA_plusB_minusP_P_to_EP_res_plusA_plusB_minusP,
        :k_EP_res_plusA_plusB_minusP_to_E_res_plusA_plusB_minusP_P,
        :k_E_res_plusA_plusB_minusP_to_EQ, :k_EQ_to_E_res_plusA_plusB_minusP,
        :k_E_Q_to_EQ, :k_EQ_to_E_Q])
    # Every name is a Julia identifier, so a parameter tuple needs no var"…".
    @test all(Base.isidentifier, ER.parameters(m))
    _testhelper_check_law_against_oracle(m)
end

@testset "Tier 2 reads the free metabolites at both ends of a step" begin
    # F is where P leaves (the fused release E(S) → F + P) and E is where S
    # enters, so the isomerization between them runs F → E however either step
    # is written.
    forward = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            E(S) <--> F + P
            F <--> E
        end
    end
    backward = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> E(S)
            F + P <--> E(S)
            E <--> F
        end
    end
    m = ER.Mechanism(forward)
    @test ER.Mechanism(backward) == m
    k = only(k for (k, e) in enumerate(m.steps) if ER.Mechanisms._is_iso(e))
    @test ER.Mechanisms._form_name(m, m.steps[k].from) == "F" &&
          ER.Mechanisms._form_name(m, m.steps[k].to) == "E"
    kf, kr = _testhelper_step_parameters(m, k)
    @test ER.Derivation.name(kf, m) == :k_F_to_E && ER.Derivation.name(kr, m) == :k_E_to_F
end

# Chokepoint guard: no `Symbol("[KkVL]...")` literal is constructed outside
# parameter-name rendering bodies (the `name(::Parameter, m)` chokepoint).

const _CHOKEPOINT_PREFIX = r"^[KkVL][_a-zA-Z0-9]"

# Extract the function-name symbol from a signature expression.
# Handles `name(...)`, `name(...) where T`, `name(...)::Ret`, etc.
function _testhelper_sig_fn_name(sig)
    while sig isa Expr && sig.head === :where
        sig = sig.args[1]
    end
    if sig isa Expr && sig.head === :(::)
        sig = sig.args[1]
    end
    return sig isa Expr && sig.head === :call ? sig.args[1] : nothing
end

# Extract the first positional arg type-annotation as a String.
function _testhelper_sig_first_arg_str(sig)
    while sig isa Expr && sig.head === :where
        sig = sig.args[1]
    end
    if sig isa Expr && sig.head === :(::)
        sig = sig.args[1]
    end
    sig isa Expr && sig.head === :call || return ""
    args = sig.args[2:end]
    pos_args = filter(a -> !(a isa Expr && a.head === :parameters), args)
    isempty(pos_args) && return ""
    return string(pos_args[1])
end

# A method definition is a chokepoint body iff it is a `name` method
# dispatching on a Parameter subtype value.
function _testhelper_is_chokepoint_def(expr)
    expr isa Expr || return false
    sig = if expr.head === :function && length(expr.args) >= 1
        expr.args[1]
    elseif expr.head === :(=) && expr.args[1] isa Expr &&
           expr.args[1].head in (:call, :where)
        expr.args[1]
    else
        return false
    end
    fn_name = _testhelper_sig_fn_name(sig)
    fn_name === :name || return false
    arg_str = _testhelper_sig_first_arg_str(sig)
    return occursin("Parameter", arg_str)
end

# Reconstruct the string content of a `Symbol("...")` call. Supports
# both literal Strings and `:string` interpolation expressions like
# `Symbol("K\$idx")` → `Expr(:string, "K", :idx)`.
function _testhelper_symbol_call_pattern(expr)
    expr isa Expr && expr.head === :call &&
        length(expr.args) >= 2 && expr.args[1] === :Symbol || return nothing
    arg2 = expr.args[2]
    if arg2 isa String
        return arg2
    elseif arg2 isa Expr && arg2.head === :string
        # Concatenate; non-String parts become a placeholder so the
        # prefix regex can still match (e.g., "K\$idx" → "K_"). The
        # placeholder must be a character class matched by
        # `_CHOKEPOINT_PREFIX` (`[_a-zA-Z0-9]`); underscore qualifies.
        return join(p isa String ? p : "_" for p in arg2.args)
    end
    return nothing
end

function _testhelper_walk_violations!(expr, in_chokepoint::Bool, out::Vector{String})
    expr isa Expr || return
    if _testhelper_is_chokepoint_def(expr)
        for child in expr.args
            _testhelper_walk_violations!(child, true, out)
        end
    else
        pat = _testhelper_symbol_call_pattern(expr)
        if pat !== nothing && occursin(_CHOKEPOINT_PREFIX, pat) && !in_chokepoint
            push!(out, "Symbol(\"$pat\")")
        end
        for child in expr.args
            _testhelper_walk_violations!(child, in_chokepoint, out)
        end
    end
end

"""The number of `name(::Parameter, m)` method definitions in `expr`."""
_testhelper_count_chokepoint_defs(expr) =
    expr isa Expr ?
    Int(_testhelper_is_chokepoint_def(expr)) +
    sum(_testhelper_count_chokepoint_defs, expr.args; init = 0) : 0

@testset "chokepoint: no Symbol(\"[KkVL]...\") outside parameter-name renderers" begin
    src_dir = normpath(joinpath(@__DIR__, "..", "..", "src"))
    names_file = joinpath(src_dir, "Derivation", "names.jl")
    for f in sort!([joinpath(d, f) for (d, _, fs) in walkdir(src_dir) for f in fs])
        endswith(f, ".jl") || continue
        src = read(f, String)
        expr = Meta.parseall(src; filename=f)
        violations = String[]
        _testhelper_walk_violations!(expr, false, violations)
        if !isempty(violations)
            @info "chokepoint violations" file=basename(f) violations
        end
        @test isempty(violations)
        # The renderers are in src/Derivation/names.jl and nowhere else.
        @test (_testhelper_count_chokepoint_defs(expr) > 0) == (f == names_file)
    end
end

@testset "base ranks order forms and edges by content, free E among the forms" begin
    # Free E is form 1, but by content it sorts after the free form D.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> D + P
            D <--> E
        end
    end
    @test [ER.Mechanisms._form_name(m, f) for f in eachindex(m.nodes)] ==
          ["E", "D", "EA", "EP"]
    @test ER.Derivation._base_rank(m) == [2, 1, 3, 4]
    # P enters at D and A at E, so the isomerization runs from D (orientation tier 2).
    texts = [_testhelper_edge_text(m, k) for k in eachindex(m.steps)]
    @test texts == ["E + A <--> EA", "D <--> E", "D + P <--> EP", "EA <--> EP"]
    # By base rank the edges from D come first.
    @test ER.Derivation._edge_base_rank(m) == [3, 1, 2, 4]
    @test [m.step_kinds[k] === :isomerization for k in eachindex(m.steps)] ==
          [false, true, false, true]
    @test [ER.Derivation._binds_product(m, k) for k in eachindex(m.steps)] ==
          [false, false, true, false]
end
