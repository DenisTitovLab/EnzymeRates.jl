# ABOUTME: Tests for the core EnzymeRates types (the mechanisms, the reaction, the
# ABOUTME: Parameter family) and their constructors, accessors, and equality.
using Serialization

const ER = EnzymeRates

# Michaelis–Menten with rapid-equilibrium bindings and a steady-state isomerization.
const _testhelper_re_mm = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ ES
        ES <--> EP
        EP ⇌ E + P
    end
end

# The forms of S ⇌ P: E, ES and EP; its metabolites are S and P.
const _testhelper_UNIUNI_FORMS = [:E, :ES, :EP]

# The steps of S ⇌ P between `_testhelper_UNIUNI_FORMS`: S binding E, ES isomerizing to EP
# at steady state, and P binding E, the bindings at rapid equilibrium; each in group 0, for
# `_testhelper_uniuni_labelled` to label.
_testhelper_uniuni_steps() =
    (bind = ER.Mechanisms.Step(1, 2, [1], Int[], true, 0),
     iso = ER.Mechanisms.Step(2, 3, Int[], Int[], false, 0),
     rel = ER.Mechanisms.Step(1, 3, [2], Int[], true, 0))

# The uni-uni steps `steps`, step `j` in kinetic group `j`.
_testhelper_uniuni_labelled(steps) =
    [ER.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, j)
     for (j, e) in enumerate(steps)]

# The mechanism over substrate S and product P whose steps are the uni-uni `steps`, step
# `j` in kinetic group `j`.
_testhelper_uniuni_of(steps) =
    _testhelper_plain_mechanism([:S, :P], _testhelper_UNIUNI_FORMS, [1], [2],
                                _testhelper_uniuni_labelled(steps))

# The allosteric mechanism over those uni-uni `steps`, group `j` tagged `tags[j]`, with
# oligomeric state `n` and the regulatory sites `sites`.
_testhelper_tagged_uniuni(steps, tags, n, sites) = ER.Mechanism(
    ER.Mechanisms.CatalyticScheme([:S, :P], _testhelper_UNIUNI_FORMS, [1], [2],
                                  _testhelper_uniuni_labelled(steps)),
    ER.Mechanisms.Allostery(tags, n, sites))

"""Whether the macro builds mechanism `m` again from the text `show` prints for it."""
_testhelper_builds_back(m) = eval(Meta.parse(sprint(show, m))) == m

@testset "a mechanism is a graph of form names and steps by index" begin
    M = ER.Mechanisms
    # Every name on a step side that is no declared metabolite
    # is a form, and the free enzyme is the form named E.
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA                         :: EqualAI
            EA + B <--> EPQ                    :: OnlyA
            EQ + P ⇌ EPQ                       :: EqualAI
            E + Q ⇌ EQ                         :: EqualAI
            (E + I ⇌ EI, EA + I ⇌ EAI)         :: OnlyI
        end
        regulatory_site(2): X::OnlyA, Y::OnlyI
    end
    c = m.catalytic
    @test c.enzyme_forms[1] === :E
    @test c.enzyme_forms == [:E, :EA, :EAI, :EI, :EPQ, :EQ]
    @test c.metabolites == [:A, :B, :P, :Q, :I]
    @test (c.substrates, c.products) == ([1, 2], [3, 4])
    chemistry = only(e for e in c.steps if c.enzyme_forms[e.from] === :EA &&
                                         c.enzyme_forms[e.to] === :EPQ)
    @test chemistry.taken_up == [findfirst(==(:B), c.metabolites)]
    @test chemistry.given_off == Int[]
    # Another writing: the steps and the groups in another order, the names declared in
    # another order, and a binding written as the release it reverses.
    rewritten = @enzyme_mechanism begin
        substrates: B, A
        products: Q, P
        modifiers: I
        oligomeric_state: 2
        steps: begin
            (EA + I ⇌ EAI, EI ⇌ E + I)         :: OnlyI
            EQ ⇌ E + Q                         :: EqualAI
            EPQ ⇌ EQ + P                       :: EqualAI
            EA + B <--> EPQ                    :: OnlyA
            E + A ⇌ EA                         :: EqualAI
        end
        regulatory_site(2): Y::OnlyI, X::OnlyA
    end
    @test rewritten == m
    @test hash(rewritten) == hash(m)

    # The free enzyme is the form named E: renamed, the macro and the constructor both
    # say so.
    err = try
        eval(:(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                Efree + S ⇌ ES
                ES <--> Efree + P
            end
        end))
    catch e
        e
    end
    @test err isa ArgumentError && occursin("named `E`", err.msg)
    steps = [M.Step(1, 2, [1], Int[], true, 1), M.Step(2, 1, Int[], [2], false, 2)]
    err = try
        ER.Mechanism(M.CatalyticScheme([:S, :P], [:Efree, :ES], [1], [2], steps),
                     M.Allostery(Symbol[], 1, M.RegulatorySite[]))
    catch e
        e
    end
    @test err isa ArgumentError && occursin("named `E`", err.msg)

    # A misspelt form name makes a new form: here EQP for EPQ, which only the inhibitor's
    # steps join, so no step joins EQP to E.
    @test_throws r"no steps lead from free E to .*EQP" eval(:(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA
            EA + B <--> EPQ
            EPQ ⇌ EQ + P
            EQ ⇌ E + Q
            EQP + I ⇌ EQPI
        end
    end))
    # Misspelt in the release of P, EQP reaches E through EQ, so the connectivity check
    # passes it; EPQ becomes a dead end, and the guard finds no step cycle that turns the
    # reaction over. The checks catch this misspelling, but by its consequence.
    @test_throws(ErrorException("Mechanism: no cycle of steps turns substrates into " *
                                "products, so the rate is zero everywhere."),
                 eval(:(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EPQ
            EQP ⇌ EQ + P
            EQ ⇌ E + Q
        end
    end)))

    # Call notation and a modifier copy are errors that name what to write.
    @test_throws "write `ES`" eval(:(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E + P
        end
    end))
    @test_throws "write `S`" eval(:(@enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> E + P
            ES + S::Mod ⇌ ESS
        end
    end))
end

@testset "a macro returns a plain value" begin
    # The macros build ordinary values: two writings of one mechanism are equal, and
    # compiling equal mechanisms gives the same rate law.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    @test !_testhelper_is_allosteric(m)
    @test m == @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            EP ⇌ E + P
            E + S ⇌ ES
            ES <--> EP
        end
    end
    @test compile_rate_equation(m) === compile_rate_equation(deepcopy(m))
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: EqualAI
            ES <--> EP   :: OnlyA
            EP ⇌ E + P     :: EqualAI
        end
    end
    @test _testhelper_is_allosteric(am)
    @test am == @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            ES <--> EP   :: OnlyA
            EP ⇌ E + P     :: EqualAI
            E + S ⇌ ES     :: EqualAI
        end
    end
    @test compile_rate_equation(am) === compile_rate_equation(deepcopy(am))
    # The plain mechanism of the same catalytic scheme compiles another law.
    @test compile_rate_equation(am) !== compile_rate_equation(_testhelper_plain_of(am))
end

@testset "Types" begin
    @testset "Mechanism struct + accessors" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end

        @test _testhelper_substrates(m) == [:S]
        @test _testhelper_products(m) == [:P]
        @test isempty(_testhelper_regulators(m))
        @test ER.metabolites(m) == [:S, :P]
        # Steps are canonicalized at construction, so compare content
        # order-independently.
        @test Set(_testhelper_edge_text(m, k)
                  for k in eachindex(m.catalytic.steps)) ==
              Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])
        @test length(_testhelper_group_members(m.catalytic)) == 3
        @test Set(_testhelper_graph(m).forms) == Set(["E", "ES", "EP"])

        # Shared kinetic-group: two steps in group 1 (regulator R binds
        # both E and ES sharing one K).
        m2 = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: R
            steps: begin
                (E + R ⇌ ER, ES + R ⇌ ERS)
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end
        # Steps are canonicalized; the two R-binding steps share one kinetic
        # group → 5 steps, 4 groups (one 2-step group binding R).
        @test length(_testhelper_group_members(m2.catalytic)) == 4
        shared = only(ks for ks in _testhelper_group_members(m2.catalytic)
                      if length(ks) == 2)
        @test all(k -> m2.catalytic.steps[k].taken_up ==
                       [_testhelper_metabolite(m2.catalytic, :R)] &&
                       isempty(m2.catalytic.steps[k].given_off),
                  shared)
    end

    @testset "metabolites() lists substrates, then products, then regulators" begin
        # A substrate, a product and a regulator, so `metabolites` lists all three roles in
        # order — coverage the plain S/P accessor tests do not reach.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: I
            steps: begin
                E + S ⇌ ES
                ES + I ⇌ EIS
                ES <--> EP
                EP ⇌ E + P
            end
        end
        @test ER.metabolites(m) == [:S, :P, :I]
    end

    @testset "a mechanism names its substrates, products and modifiers" begin
        # Modifiers, and a reactant that a step takes up as a modifier. The concentration
        # order is substrates, products, then the modifiers, each by name, as the compiled
        # rate law lists them.
        m = @enzyme_mechanism begin
            substrates: B, A
            products: Q, P
            modifiers: J, I
            steps: begin
                E + A ⇌ EA
                EA + B ⇌ EAB
                EAB <--> EPQ
                EPQ <--> EQ + P
                EQ ⇌ E + Q
                E + J ⇌ EJ
                EA + I ⇌ EAI
                EQ + A ⇌ EAmodQ
            end
        end
        @test m.catalytic.metabolites == [:A, :B, :P, :Q, :I, :J]
        @test (m.catalytic.substrates, m.catalytic.products) == ([1, 2], [3, 4])
        @test ER.metabolites(m) == [:A, :B, :P, :Q, :I, :J]
        @test ER.metabolites(compile_rate_equation(m)) == (:A, :B, :P, :Q, :I, :J)
        # A reaction's competitive inhibitor that no step binds is no concentration of
        # its init mechanisms.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: K
        end
        on_rxn = only(_testhelper_mechanisms(ER.Enumeration.init_mechanisms(rxn)))
        @test _testhelper_regulators(on_rxn) == Symbol[]
        @test ER.metabolites(on_rxn) == [:S, :P]
        # An allosteric mechanism lists its site ligands last, site by site; a name that
        # is both a modifier and a ligand is listed once.
        am = @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: I
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES     :: EqualAI
                ES <--> EP   :: OnlyA
                EP ⇌ E + P     :: EqualAI
                E + I ⇌ EI     :: EqualAI
            end
            regulatory_site(2): Y::OnlyI
            regulatory_site(2): X::OnlyA, I::OnlyI
        end
        @test am.catalytic.metabolites == [:S, :P, :I]
        @test ER.metabolites(am) == [:S, :P, :I, :X, :Y]
        @test Tuple(ER.metabolites(am)) == ER.metabolites(compile_rate_equation(am))
    end

    @testset "Allosteric mechanism from the DSL" begin
        # Single-ligand :EqualAI reg site is allowed (degenerate but valid —
        # the enumerator won't emit it, but users may write it for teaching).
        @test _testhelper_is_allosteric(@enzyme_mechanism(begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: NonequalAI
                EP ⇌ E + P      :: NonequalAI
            end
            regulatory_site(2): I::EqualAI
        end))

        # Any catalytic group may be :OnlyI, but the active conformation cannot switch
        # off a rapid-equilibrium step whose two sides its other steps still join.
        why = " is tagged :OnlyI, but the steps that conformation keeps still join its " *
              "two sides. Thermodynamics then fixes its equilibrium constant, and a step " *
              "keeps its kind in both conformations, so it cannot be switched off."
        lead = "Mechanism: in the active conformation the rapid-equilibrium step "
        pinned(step) = lead * step * why
        @test_throws pinned("E + S ⇌ ES") @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: OnlyI
                ES <--> EP    :: NonequalAI
                EP ⇌ E + P      :: NonequalAI
            end
        end
        # Nor one that binds X and turns S into P at once: the chemistry and the X binding
        # at EP still join its sides.
        @test_throws pinned("ES + X ⇌ EPX") @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: X
            steps: begin
                E + S <--> ES          :: EqualAI
                ES <--> EP           :: EqualAI
                E + P <--> EP          :: EqualAI
                ES + X ⇌ EPX       :: OnlyI
                EP + X <--> EPX    :: EqualAI
            end
        end
        @test_throws pinned("ES + S ⇌ EPSmod") @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> ES                  :: EqualAI
                ES <--> EP                   :: EqualAI
                E + P <--> EP                  :: EqualAI
                ES + S ⇌ EPSmod     :: OnlyI
                EP + S <--> EPSmod  :: EqualAI
            end
        end
        onlyi = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: I
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> E + P   :: EqualAI
                E + I ⇌ EI      :: OnlyI
            end
        end
        @test :OnlyI in onlyi.allosteric.allosteric_states

        # Build via DSL
        m = @enzyme_mechanism begin
            substrates: S
            products:   P

            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES     :: EqualAI
                ES <--> EP   :: OnlyA
                EP ⇌ E + P     :: EqualAI
            end
            regulatory_site(2): I::OnlyI
        end
        @test m.allosteric.oligomeric_state == 2
        # Catalysis (iso step) is :OnlyA, both bindings :EqualAI. Group order is
        # canonical, so identify the tagged group by its step, not its position.
        members = _testhelper_group_members(m.catalytic)
        onlyA_g = only(g for g in eachindex(members)
                       if m.allosteric.allosteric_states[g] === :OnlyA)
        @test !ER.Enumeration._is_binding(m.catalytic.steps[first(members[onlyA_g])])
        @test all(m.allosteric.allosteric_states[g] === :EqualAI
                  for g in eachindex(members) if g != onlyA_g)
        @test only(m.allosteric.regulatory_sites).allosteric_states ==
              [:OnlyI]
    end

    @testset "a mechanism prints as the text of the macro that builds it" begin
        # The macro builds the same mechanism from the text (`_testhelper_builds_back`).
        # One step per line in stored order, a release written as the binding it
        # reverses, a form with nothing bound as its bare label.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> ES
                ES <--> E + P
            end
        end
        @test sprint(show, m) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + S <--> ES
                    E + P <--> ES
                end
            end"""
        @test _testhelper_builds_back(m)

        # A block with no steps names no form E, so the macro builds no mechanism to print.
        @test_throws(ArgumentError("Mechanism: no enzyme form is named `E`; the free " *
                                   "enzyme is the form named `E`"),
                     @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
            end
        end)

        # A branched mechanism: every step, each form with its bound metabolites.
        m_b = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> EA
                E + B <--> EB
                EA + B <--> EAB
                EB + A <--> EAB
                EAB <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
            end
        end
        @test sprint(show, m_b) == """
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    E + B <--> EB
                    E + Q <--> EQ
                    EA + B <--> EAB
                    EAB <--> EPQ
                    EB + A <--> EAB
                    EQ + P <--> EPQ
                end
            end"""
        @test _testhelper_builds_back(m_b)

        # Catalytic 3-cycle: the steps print in stored order, not chain order, and
        # rapid-equilibrium and steady-state steps keep their arrows.
        m_re = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end
        @test sprint(show, m_re) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + P ⇌ EP
                    E + S ⇌ ES
                    ES <--> EP
                end
            end"""
        @test _testhelper_builds_back(m_re)

        # Theorell–Chance: B binds and P leaves in one step, written on one line.
        m_tc = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> EA
                EA + B <--> EQ + P
                EQ <--> E + Q
            end
        end
        @test sprint(show, m_tc) == """
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    E + Q <--> EQ
                    EA + B <--> EQ + P
                end
            end"""
        @test _testhelper_builds_back(m_tc)

        # A fused binding (B binds and chemistry runs in one step) prints B on its entry
        # side, since EPQ does not name it.
        m_fb = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> EA
                EA + B <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
            end
        end
        @test sprint(show, m_fb) == """
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    E + Q <--> EQ
                    EA + B <--> EPQ
                    EQ + P <--> EPQ
                end
            end"""
        @test _testhelper_builds_back(m_fb)

        # The modifiers: line lists the metabolites that are neither substrates nor
        # products; a kinetic group of several steps prints as a tuple.
        @test !contains(sprint(show, _testhelper_re_mm), "modifiers")
        m_reg = @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: I
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                (E + I ⇌ EI, EP + I ⇌ EIP)
            end
        end
        @test sprint(show, m_reg) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                modifiers: I
                steps: begin
                    (E + I ⇌ EI, EP + I ⇌ EIP)
                    E + P ⇌ EP
                    E + S ⇌ ES
                    ES <--> EP
                end
            end"""
        @test _testhelper_builds_back(m_reg)
        m_copy = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
                EP + S ⇌ EPSmod
            end
        end
        @test sprint(show, m_copy) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + P ⇌ EP
                    E + S ⇌ ES
                    EP + S ⇌ EPSmod
                    ES <--> EP
                end
            end"""
        @test _testhelper_builds_back(m_copy)

        # A ping-pong: a form that keeps a residual prints it, with or without bound
        # metabolites.
        m_pp = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA <--> EP_res_plusA_minusP
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
                E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
                EB_res_plusA_minusP <--> EQ
                E + Q ⇌ EQ
            end
        end
        @test sprint(show, m_pp) == """
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ EA
                    E + Q ⇌ EQ
                    EA <--> EP_res_plusA_minusP
                    EB_res_plusA_minusP <--> EQ
                    E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
                    E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
                end
            end"""
        @test _testhelper_builds_back(m_pp)

        # EnzymeReaction with oligomeric_state > 1.
        rxn_oligo = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: 4
        end
        @test sprint(show, rxn_oligo) == """
            @enzyme_reaction begin
                substrates: S[C]
                products: P[C]
                oligomeric_state: 4
            end"""

        # EnzymeReaction with shared_catalytic_site: surfaced in show.
        rxn_shared = @enzyme_reaction begin
            substrates: A[C], B[C]
            products:   P[C], Q[C]
            shared_catalytic_site: (A, P)
        end
        @test sprint(show, rxn_shared) == """
            @enzyme_reaction begin
                substrates: A[C], B[C]
                products: P[C], Q[C]
                shared_catalytic_site: (A, P)
            end"""

        # An allosteric mechanism adds its regulators, its oligomeric state, a tag
        # after each kinetic group and a line per regulatory site.
        m_allo = @enzyme_mechanism begin
            substrates: F6P
            products:   F16BP
            oligomeric_state: 2
            steps: begin
                E + F6P ⇌ EF6P         :: EqualAI
                EF6P <--> EF16BP     :: EqualAI
                EF16BP ⇌ E + F16BP     :: EqualAI
            end
            regulatory_site(2): I::OnlyI
        end
        @test sprint(show, m_allo) == """
            @enzyme_mechanism begin
                substrates: F6P
                products: F16BP
                oligomeric_state: 2
                steps: begin
                    E + F16BP ⇌ EF16BP :: EqualAI
                    E + F6P ⇌ EF6P :: EqualAI
                    EF6P <--> EF16BP :: EqualAI
                end
                regulatory_site(2): I::OnlyI
            end"""
        @test _testhelper_builds_back(m_allo)

        # A regulatory site with two ligands lists both with their tags.
        m_site = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES     :: EqualAI
                ES <--> EP   :: OnlyA
                EP ⇌ E + P     :: EqualAI
            end
            regulatory_site(3): I::OnlyI, J::EqualAI
        end
        @test sprint(show, m_site) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP :: EqualAI
                    E + S ⇌ ES :: EqualAI
                    ES <--> EP :: OnlyA
                end
                regulatory_site(3): I::OnlyI, J::EqualAI
            end"""
        @test _testhelper_builds_back(m_site)
    end

    @testset "the macro rejects a rapid-equilibrium segment with no bottom form" begin
        # Random product release from EPQ at rapid equilibrium, but release from EP
        # and EQ at steady state: the RE segment {EP, EPQ, EQ} has weights
        # P·Q : Kq·P : Kp·Q, all zero at P = Q = 0. The RE approximation carries no
        # parameter for how EPQ splits between EP and EQ there, so the rate is
        # 0/0 at zero products: class 3 rejects that face (`_BOTTOMLESS`), and the RE-speed
        # oracle sees the steady state there move with the release speeds.
        random_release = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                E + P <--> EP
                E + Q <--> EQ
                EP + Q ⇌ EPQ
                EQ + P ⇌ EPQ
                ES <--> EPQ
            end
        end
        v = ER.Derivation._guard(random_release)
        @test (v.fault, v.face) == (ER.Derivation._BOTTOMLESS, 2)
        msg = ER.Derivation._guard_message(random_release)
        @test occursin("has no form whose weight keeps the segment's lowest power of " *
                       "each zeroed metabolite", msg)
        @test msg ==
              "Mechanism: at zero P and Q, the rapid-equilibrium segment " *
              "{EP, EPQ, EQ} has no form whose weight keeps the segment's lowest " *
              "power of each zeroed metabolite, so the forms' shares depend on how those " *
              "concentrations vanish and the rate is 0/0 there. Make the kinetic group " *
              "of EP + Q ⇌ EPQ steady-state."
        # The printed fix passes the guard.
        @test ER.Derivation._guard(@_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                E + P <--> EP
                E + Q <--> EQ
                EP + Q <--> EPQ
                EQ + P ⇌ EPQ
                ES <--> EPQ
            end
        end) === ER.Derivation._WELL_POSED
        @test _testhelper_respeed_verdicts(random_release; rng = Random.MersenneTwister(1),
                                           law = :limit)[2] === :pattern
        @test_throws ErrorException(msg) @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                E + P <--> EP
                E + Q <--> EQ
                EP + Q ⇌ EPQ
                EQ + P ⇌ EPQ
                ES <--> EPQ
            end
        end

        # Mirror on the substrate side: the segment {EA, EB, EAB} is 0/0 at
        # A = B = 0 (face 3).
        mirror = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> EA
                E + B <--> EB
                EA + B ⇌ EAB
                EB + A ⇌ EAB
                E + P ⇌ EP
                EAB <--> EP
            end
        end
        v = ER.Derivation._guard(mirror)
        @test (v.fault, v.face) == (ER.Derivation._BOTTOMLESS, 3)
        @test _testhelper_respeed_verdicts(mirror; rng = Random.MersenneTwister(2),
                                           law = :limit)[3] === :pattern
        @test_throws ErrorException(
            ER.Derivation._guard_message(mirror)) @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> EA
                E + B <--> EB
                EA + B ⇌ EAB
                EB + A ⇌ EAB
                E + P ⇌ EP
                EAB <--> EP
            end
        end

        # The same segment with A's modifier copy binding EB at rapid equilibrium:
        # EAmodB joins it, but every form still carries A or B (the copy's
        # concentration is A's), so the segment stays bottomless, whichever order the
        # groups are written in.
        copy_last = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> EA
                E + B <--> EB
                EA + B ⇌ EAB
                EB + A ⇌ EAB
                EB + A ⇌ EAmodB
                E + P ⇌ EP
                EAB <--> EP
            end
        end
        copy_first = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                EB + A ⇌ EAmodB
                EAB <--> EP
                E + P ⇌ EP
                EB + A ⇌ EAB
                EA + B ⇌ EAB
                E + B <--> EB
                E + A <--> EA
            end
        end
        v = ER.Derivation._guard(copy_last)
        @test (v.fault, v.face) == (ER.Derivation._BOTTOMLESS, 3)
        @test ER.Derivation._guard(copy_first) == v
        @test copy_first == copy_last
        @test occursin("segment {EA, EAB, EB, EAmodB} has no form",
                       ER.Derivation._guard_message(copy_last))
        @test_throws ErrorException(
            ER.Derivation._guard_message(copy_last)) @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> EA
                E + B <--> EB
                EA + B ⇌ EAB
                EB + A ⇌ EAB
                EB + A ⇌ EAmodB
                E + P ⇌ EP
                EAB <--> EP
            end
        end

        # Ordered release: the segment {EP, EPQ} has weights Kq : Q, so EP is its
        # bottom form and the rate is defined at zero products; every face agrees with
        # the RE-speed oracle.
        m_ordered = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES
                E + P <--> EP
                EP + Q ⇌ EPQ
                ES <--> EPQ
            end
        end
        @test all(==(:agree),
                  _testhelper_respeed_verdicts(m_ordered; rng = Random.MersenneTwister(3)))

        # Mixed abortive complex: the segment {ES, EP, EPS} has weights
        # K·S : K′·P : S·P, zero only at S = P = 0, a corner outside the guard's faces
        # where no turnover is possible anyway. Accepted.
        m_abortive = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                E + P <--> EP
                ES + P ⇌ EPS
                EP + S ⇌ EPS
                ES <--> EP
            end
        end
        @test all(==(:agree),
                  _testhelper_respeed_verdicts(m_abortive; rng = Random.MersenneTwister(4)))
    end

    @testset "a free form that holds the enzyme it receives is a second closed class" begin
        # Every step that takes up the substrate A is :OnlyA, and the conformations flip at
        # free E, so the inactive conformation keeps E and EP alone: a second free form
        # F, which binds B and joins E only through F + B <--> EB, holds no inactive
        # enzyme. In the active conformation, with every substrate absent, or A, or B, F
        # holds the enzyme it receives, together with EB while B is present, and never
        # returns it to E: a second closed class, which class 3 rejects (`_NONUNIQUE`) and
        # the RE-speed oracle sees.
        s1 = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> EA                   :: OnlyA
                EA + B <--> EAB             :: EqualAI
                EAB <--> EP                 :: OnlyA
                EP <--> E + P                   :: EqualAI
                E + A ⇌ EAmod            :: OnlyA
                EAmod + B ⇌ EAmodB      :: EqualAI
                EB + A ⇌ EAmodB      :: EqualAI
                F + B <--> EB                   :: EqualAI
            end
        end
        @test _testhelper_face_faults(s1)[3:5] == fill(ER.Derivation._NONUNIQUE, 3)
        v = ER.Derivation._guard(s1)
        @test (v.fault, v.face) == (ER.Derivation._NONUNIQUE, 3)
        @test _testhelper_respeed_verdicts(s1; rng = Random.MersenneTwister(12),
                                           law = :limit)[3:5] == fill(:nonunique, 3)
        @test_throws ErrorException(
            ER.Derivation._guard_message(s1)) @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> EA                   :: OnlyA
                EA + B <--> EAB             :: EqualAI
                EAB <--> EP                 :: OnlyA
                EP <--> E + P                   :: EqualAI
                E + A ⇌ EAmod            :: OnlyA
                EAmod + B ⇌ EAmodB      :: EqualAI
                EB + A ⇌ EAmodB      :: EqualAI
                F + B <--> EB                   :: EqualAI
            end
        end
        # Two variants keep F and pass every face: F <--> E tagged :OnlyA; or F <--> E
        # untagged with the modifier binding of A untagged, so that no :OnlyA binding has
        # its two sides joined by untagged steps.
        f_onlya = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> EA                   :: OnlyA
                EA + B <--> EAB             :: EqualAI
                EAB <--> EP                 :: OnlyA
                EP <--> E + P                   :: EqualAI
                E + A ⇌ EAmod            :: OnlyA
                EAmod + B ⇌ EAmodB      :: EqualAI
                EB + A ⇌ EAmodB      :: EqualAI
                F + B <--> EB                   :: EqualAI
                F <--> E                          :: OnlyA
            end
        end
        copy_untagged = @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> EA                   :: OnlyA
                EA + B <--> EAB             :: EqualAI
                EAB <--> EP                 :: OnlyA
                EP <--> E + P                   :: EqualAI
                E + A ⇌ EAmod            :: EqualAI
                EAmod + B ⇌ EAmodB      :: EqualAI
                EB + A ⇌ EAmodB      :: EqualAI
                F + B <--> EB                   :: EqualAI
                F <--> E                          :: EqualAI
            end
        end
        # Both oracles agree on every face: the RE-speed oracle (class 3) and the ε-path
        # (admission of the :OnlyA tags).
        for (i, variant) in enumerate((f_onlya, copy_untagged))
            @test all(_testhelper_admission_agreement(
                variant; rng = Random.MersenneTwister(12 + i)))
        end
    end

    @testset "a second free form is admitted where a reading drains it" begin
        # The conformations flip at free E alone. A second free form joined to E by a
        # steady-state step is an ordinary form, whatever its label (R1, and its twin
        # labelled D, which sorts before E). Where tags strand it under the guard's built
        # readings, which stop a tagged isomerization and cannot switch off a
        # rapid-equilibrium one, admission's verdicts under them stand: every step at E
        # :OnlyA leaves F fed at zero product with a dead way back (`_PUMPED`); a tagged
        # steady-state isomerization alone joining F to E makes F an island (`_ISLAND`);
        # a tagged rapid-equilibrium one cannot be switched off (`_FAST_STEP_PINNED`). The
        # search admits each under a reading in which the isomerization loses its rate
        # into F, so F drains into free E, the oracles agreeing on every face. Each also
        # has an admitted variant under the built readings: the S binding untagged (AD2d),
        # or the catalysis tagged in place of the isomerization.
        r1 = @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: NonequalAI
                F + P ⇌ EP      :: EqualAI
                F <--> E          :: NonequalAI
            end
        end
        r1_d = @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: NonequalAI
                D + P ⇌ EP      :: EqualAI
                D <--> E          :: NonequalAI
            end
        end
        # R1's oracle is the RE-speed one (no tag excludes a step); the law also equals the
        # exact limit at one subunit, where L means [E_I]/[E_A] for D too.
        for (i, m) in enumerate((r1, r1_d))
            @test all(==(:agree), _testhelper_respeed_verdicts(
                m; rng = Random.MersenneTwister(20 + i)))
            _testhelper_check_law_against_limit(m; seed = 20 + i)
        end
        r2 = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: OnlyA
                ES <--> EP    :: EqualAI
                F + P ⇌ EP      :: EqualAI
                F <--> E          :: OnlyA
            end
        end
        v = ER.Derivation._verdict_under(r2, ER.Derivation._built_reading(r2))
        @test (v.fault, v.face) == (ER.Derivation._PUMPED, 2)
        @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(23),
                                        face = Set([:P])) === :direction_dependent
        @test @enzyme_mechanism(begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: OnlyA
                ES <--> EP    :: EqualAI
                F + P ⇌ EP      :: EqualAI
                F <--> E          :: OnlyA
            end
        end) == r2
        ad2d = @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                F + P ⇌ EP      :: EqualAI
                F <--> E          :: OnlyA
            end
        end
        @test _testhelper_epath_verdict(ad2d;
                                        rng = Random.MersenneTwister(24)) === :converges
        r3 = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: I
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                E + P ⇌ EP      :: EqualAI
                E <--> F          :: OnlyA
                F + I ⇌ FI      :: EqualAI
            end
        end
        r4 = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: NonequalAI
                E + P ⇌ EP      :: EqualAI
                F ⇌ E             :: OnlyA
            end
        end
        @test ER.Derivation._admission_interior(
            r3, ER.Derivation._built_reading(r3)).fault === ER.Derivation._ISLAND
        @test _testhelper_epath_verdict(r3; rng = Random.MersenneTwister(26)) in
              (:direction_dependent, :wrong_limit)
        v = ER.Derivation._admission_interior(r4, ER.Derivation._built_reading(r4))
        @test (v.fault, v.conformation) == (ER.Derivation._FAST_STEP_PINNED, 2)
        @test ER.Mechanisms.step_text(r4.catalytic, v.edge) == "E ⇌ F"
        @test _testhelper_epath_verdict(r4; rng = Random.MersenneTwister(27)) === :no_path
        for (i, x) in enumerate((r2, r3, r4))
            @test ER.Derivation._guard(x) == ER.Derivation._WELL_POSED
            @test all(_testhelper_admission_agreement(
                x; rng = Random.MersenneTwister(30 + i), corners = true))
            _testhelper_check_law_against_limit(x)
        end
        @test @enzyme_mechanism(begin
            substrates: S
            products: P
            modifiers: I
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                E + P ⇌ EP      :: EqualAI
                E <--> F          :: OnlyA
                F + I ⇌ FI      :: EqualAI
            end
        end) == r3
        @test @enzyme_mechanism(begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: NonequalAI
                E + P ⇌ EP      :: EqualAI
                F ⇌ E             :: OnlyA
            end
        end) == r4
        r3_base = @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: I
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyA
                E + P ⇌ EP      :: EqualAI
                E <--> F          :: EqualAI
                F + I ⇌ FI      :: EqualAI
            end
        end
        r4_base = @enzyme_mechanism begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyA
                E + P ⇌ EP      :: EqualAI
                F ⇌ E             :: EqualAI
            end
        end
        for (i, m) in enumerate((r3_base, r4_base))
            @test _testhelper_epath_verdict(
                m; rng = Random.MersenneTwister(27 + i)) ===
                :converges
        end
        # F in E's rapid-equilibrium segment, untagged, is an ordinary form too (A5).
        @test _testhelper_is_allosteric(@enzyme_mechanism(begin
            substrates: S
            products: P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: NonequalAI
                E + P ⇌ EP      :: EqualAI
                F ⇌ E             :: NonequalAI
            end
        end))
    end

    @testset "the inactive conformation's segments need a bottom form too" begin
        # The :OnlyA P group holds the fused binding EQ + P ⇌ EstarS beside the plain
        # binding E + P ⇌ EP. In the active conformation that binding joins EP,
        # EPQ and EQ to free E at rapid equilibrium, and E is the segment's bottom
        # form. The inactive conformation lacks it and reaches those forms only through
        # the steady-state chemistry and the steady-state Q binding at E, so its segment
        # {EP, EPQ, EQ} holds a product in every form: class 3 finds it bottomless
        # at zero products. The same untagged steps join the two sides of the tagged
        # rapid-equilibrium binding E + P ⇌ EP, so the inactive conformation cannot
        # switch it off: admission rejects the tag set first (`_FAST_STEP_PINNED`), and the
        # ε-path oracle finds no path.
        am = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES                               :: EqualAI
                ES <--> EPQ                          :: EqualAI
                EPQ ⇌ EP + Q                         :: EqualAI
                EPQ ⇌ EQ + P                         :: EqualAI
                (E + P ⇌ EP, EQ + P ⇌ EstarS)        :: OnlyA
                E + Q <--> EQ                            :: EqualAI
            end
        end
        @test _testhelper_face_faults(am)[2] === ER.Derivation._BOTTOMLESS
        v = ER.Derivation._guard(am)
        @test v.fault === ER.Derivation._FAST_STEP_PINNED
        @test ER.Mechanisms.step_text(am.catalytic, v.edge) == "E + P ⇌ EP"
        @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(5)) === :no_path
        @test_throws ErrorException(
            ER.Derivation._guard_message(am)) @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ ES                               :: EqualAI
                ES <--> EPQ                          :: EqualAI
                EPQ ⇌ EP + Q                         :: EqualAI
                EPQ ⇌ EQ + P                         :: EqualAI
                (E + P ⇌ EP, EQ + P ⇌ EstarS)        :: OnlyA
                E + Q <--> EQ                            :: EqualAI
            end
        end
    end

    @testset "Allosteric mechanism validators" begin
        # A mechanism with no allosteric states is plain: one subunit and no regulatory
        # site. The filters read a plain mechanism's oligomeric state, and two plain
        # values on one scheme must be equal, so any other allostery is refused.
        c = _testhelper_re_mm.catalytic
        no_site = ER.Mechanisms.RegulatorySite[]
        plain = ER.Mechanism(c, ER.Mechanisms.Allostery(Symbol[], 1, no_site))
        @test plain == _testhelper_re_mm && hash(plain) == hash(_testhelper_re_mm)
        @test (plain.allosteric.oligomeric_state, plain.allosteric.regulatory_sites) ==
              (1, no_site)
        @test_throws(ErrorException(
            "Mechanism: a mechanism with no allosteric states has oligomeric state 1 " *
            "and no regulatory site; got oligomeric state 2 and 0 regulatory sites"),
            ER.Mechanism(c, ER.Mechanisms.Allostery(Symbol[], 2, no_site)))
        @test_throws(ErrorException(
            "Mechanism: a mechanism with no allosteric states has oligomeric state 1 " *
            "and no regulatory site; got oligomeric state 1 and 1 regulatory site"),
            ER.Mechanism(c, ER.Mechanisms.Allostery(
                Symbol[], 1, [ER.Mechanisms.RegulatorySite([:R], 1, [:OnlyA])])))

        # Wrong-length allosteric_states (4 entries for 3 kinetic groups) → error
        @test_throws ErrorException ER.Mechanism(
            _testhelper_re_mm.catalytic,
            ER.Mechanisms.Allostery([:NonequalAI, :NonequalAI, :NonequalAI, :OnlyA], 2,
                                    ER.Mechanisms.RegulatorySite[]))

        # A ligand on two regulatory sites → error: both sites would render its regulator
        # parameter name
        @test_throws(ErrorException(
            "Mechanism: allosteric regulator `A` is listed at two regulatory sites; " *
            "an allosteric regulator sits at exactly one regulatory site."),
            ER.Mechanism(
                _testhelper_re_mm.catalytic,
                ER.Mechanisms.Allostery(
                    [:EqualAI, :EqualAI, :EqualAI], 2,
                    [ER.Mechanisms.RegulatorySite([:A], n, [:OnlyA]) for n in (2, 4)])))

        # Invalid allosteric state value → error
        @test_throws ErrorException @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: NonequalAI
                EP ⇌ E + P      :: NotAState
            end
        end

        # Reg site with all-:EqualAI ligands is allowed (degenerate but valid).
        @test _testhelper_is_allosteric(@enzyme_mechanism(begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: NonequalAI
                EP ⇌ E + P      :: NonequalAI
            end
            regulatory_site(2): I::EqualAI, J::EqualAI
        end))

        # Invalid reg-site ligand allosteric state → error
        @test_throws ErrorException @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: NonequalAI
                EP ⇌ E + P      :: NonequalAI
            end
            regulatory_site(2): I::NotAState
        end
    end

    @testset "an allosteric mechanism prints every tag and every site" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES     :: NonequalAI
                ES <--> EP   :: OnlyA
                EP ⇌ E + P     :: EqualAI
            end
            regulatory_site(2): I::NonequalAI
            regulatory_site(2): J::OnlyI
        end
        s = sprint(show, m)
        @test s == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP :: EqualAI
                    E + S ⇌ ES :: NonequalAI
                    ES <--> EP :: OnlyA
                end
                regulatory_site(2): I::NonequalAI
                regulatory_site(2): J::OnlyI
            end"""
        @test _testhelper_builds_back(m)
        # Every catalytic state appears in allosteric_states line
        for g in eachindex(m.allosteric.allosteric_states)
            @test occursin(string(m.allosteric.allosteric_states[g]), s)
        end
        # No :NonequalAI ligand silently hidden from reg-site display
        for site in m.allosteric.regulatory_sites,
            (lig, state) in zip(site.ligands,
                                site.allosteric_states)
            @test occursin("$lig::$state", s)
        end
    end

    @testset "an allosteric kinetic group prints as one tagged tuple" begin
        # A mechanism with no regulatory site prints no `regulatory_site(n):` line, but
        # its oligomeric state always.
        am = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                (E + S ⇌ ES, EP + S ⇌ EPS)  :: EqualAI
                ES <--> EP                  :: EqualAI
                EP ⇌ E + P                  :: EqualAI
            end
        end
        @test repr(am) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                oligomeric_state: 2
                steps: begin
                    E + P ⇌ EP :: EqualAI
                    (E + S ⇌ ES, EP + S ⇌ EPS) :: EqualAI
                    ES <--> EP :: EqualAI
                end
            end"""
        @test _testhelper_builds_back(am)
    end

    # ─── Concrete type hierarchy ──────────────────────────────────────

    @testset "reaction fields are its declaration" begin
        rxn = @enzyme_reaction begin
            substrates: B[N], A[C2H]
            products: P[C2HN]
            competitive_inhibitors: I
            allosteric_regulators: X::Activator, Y::Inhibitor, Z
            oligomeric_state: (2, 1)
            shared_catalytic_site: (A, P)
        end
        @test [s.name for s in rxn.substrates] == [:A, :B]
        @test rxn.substrates[1].atoms == [:C => 2, :H => 1]
        @test rxn.products == [ER.Mechanisms.Reactant(:P, [:N => 1, :H => 1, :C => 2])]
        @test rxn.regulators == [ER.Mechanisms.Regulator(:I, :competitive, :unspecified),
                                 ER.Mechanisms.Regulator(:X, :allosteric, :activator),
                                 ER.Mechanisms.Regulator(:Y, :allosteric, :inhibitor),
                                 ER.Mechanisms.Regulator(:Z, :allosteric, :unspecified)]
        @test rxn.allowed_oligomeric_states == [1, 2]
        @test rxn.shared_catalytic_site == [(:A, :P)]
        @test sprint(show, rxn) ==
              "@enzyme_reaction begin\n    substrates: A[C2H], B[N]\n" *
              "    products: P[C2HN]\n    competitive_inhibitors: I\n" *
              "    allosteric_regulators: X::Activator, Y::Inhibitor, Z\n" *
              "    oligomeric_state: (1, 2)\n    shared_catalytic_site: (A, P)\nend"
    end

    @testset "Reactant and Regulator compare by content" begin
        @test ER.Mechanisms.Reactant(:X, [:C => 1]) == ER.Mechanisms.Reactant(:X, [:C => 1])
        @test ER.Mechanisms.Reactant(:X, [:C => 1]) != ER.Mechanisms.Reactant(:Y, [:C => 1])
        @test hash(ER.Mechanisms.Reactant(:X, [:C => 1])) ==
              hash(ER.Mechanisms.Reactant(:X, [:C => 1]))
        a = ER.Mechanisms.Regulator(:A, :allosteric)
        @test a == ER.Mechanisms.Regulator(:A, :allosteric, :unspecified)
        @test hash(a) == hash(ER.Mechanisms.Regulator(:A, :allosteric))
        @test a != ER.Mechanisms.Regulator(:A, :competitive)
        @test a != ER.Mechanisms.Regulator(:A, :allosteric, :activator)
    end

    @testset "RegulatorySite: validation + accessors" begin
        # A site's ligands are names.
        lig_a, lig_b = :A, :B
        site = ER.Mechanisms.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA, :NonequalAI],
        )
        @test site.ligands == [:A, :B]
        @test site.multiplicity == 4
        @test site.allosteric_states == [:OnlyA, :NonequalAI]

        # Mismatched ligand / allosteric_states length → error
        @test_throws ErrorException ER.Mechanisms.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA])

        # Multiplicity < 1 → error
        @test_throws ErrorException ER.Mechanisms.RegulatorySite(
            [lig_a], 0, [:OnlyA])

        # No ligands → error
        @test_throws ErrorException ER.Mechanisms.RegulatorySite(Symbol[], 2, Symbol[])

        # Invalid allosteric state → error
        @test_throws ErrorException ER.Mechanisms.RegulatorySite(
            [lig_a], 1, [:NotAState])

        # All four allowed states accepted; the R/T symbols are rejected.
        for st in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
            @test ER.Mechanisms.RegulatorySite([lig_a], 1, [st]) isa
                  ER.Mechanisms.RegulatorySite
        end
        for st in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT)
            @test_throws ErrorException ER.Mechanisms.RegulatorySite([lig_a], 1, [st])
        end

        # Equality / hash
        site2 = ER.Mechanisms.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA, :NonequalAI])
        @test site == site2
        @test hash(site) == hash(site2)
        # The constructor sorts the ligands by name and carries each ligand's state
        # along, so the order the ligands are given in does not matter...
        site_reordered = ER.Mechanisms.RegulatorySite(
            [lig_b, lig_a], 4, [:NonequalAI, :OnlyA])
        @test site_reordered.ligands == [lig_a, lig_b]
        @test site_reordered.allosteric_states == [:OnlyA, :NonequalAI]
        @test site == site_reordered
        @test hash(site) == hash(site_reordered)
        # ...but the pairing of ligand and state does: B::OnlyA, A::NonequalAI is
        # a different site.
        site_swapped = ER.Mechanisms.RegulatorySite(
            [lig_b, lig_a], 4, [:OnlyA, :NonequalAI])
        @test site_swapped.allosteric_states == [:NonequalAI, :OnlyA]
        @test site != site_swapped
    end


    @testset "Reactant canonicalizes atom ordering" begin
        ra1 = ER.Mechanisms.Reactant(:ATP, [:C => 10, :H => 16, :N => 5])
        ra2 = ER.Mechanisms.Reactant(:ATP, [:N => 5, :H => 16, :C => 10])
        @test ra1 == ra2
        @test hash(ra1) == hash(ra2)
        @test ra1.name === :ATP
        @test ra1.atoms == [:C => 10, :H => 16, :N => 5]

        # An element listed twice counts once with the summed count, so every writing of
        # one formula gives one value.
        ra3 = ER.Mechanisms.Reactant(:ATP, [:C => 4, :N => 5, :C => 6, :H => 16])
        @test ra3 == ra1 && hash(ra3) == hash(ra1)
        @test ra3.atoms == [:C => 10, :H => 16, :N => 5]
    end

    @testset "Reactant validation" begin
        # Mandatory atoms: empty atom list is rejected.
        @test_throws ErrorException ER.Mechanisms.Reactant(:S, Pair{Symbol,Int}[])
        # Positive counts: zero / negative rejected.
        @test_throws ErrorException ER.Mechanisms.Reactant(:S, [:C => 0])
        @test_throws ErrorException ER.Mechanisms.Reactant(:S, [:C => -1])
        # Bool is not a valid count. Use untyped vector so Bool survives to
        # the constructor (Pair{Symbol,Int}[:C => true] would convert true→1).
        @test_throws ErrorException ER.Mechanisms.Reactant(:S, [:C => true])
        # An element is what a formula in `@enzyme_reaction` writes: a capital letter and
        # any lowercase letters.
        for el in (:c, :CL, Symbol("C1"))
            @test_throws(
                ErrorException("Reactant: S has the element $el; an " *
                               "element is a capital letter and any lowercase " *
                               "letters, such as C or Cl"),
                ER.Mechanisms.Reactant(:S, [el => 1]))
        end
        # Valid construction still works.
        @test ER.Mechanisms.Reactant(:S, [:C => 6, :H => 12]).atoms == [:C => 6, :H => 12]
    end

    @testset "EnzymeReaction validation" begin
        M = ER.Mechanisms
        S = M.Reactant(:S, [:C => 1])
        P = M.Reactant(:P, [:C => 1])
        Punbal = M.Reactant(:P, [:C => 2])
        noregs = M.Regulator[]
        regS = M.Regulator(:S, :competitive)
        regA = M.Regulator(:A, :competitive)
        regA2 = M.Regulator(:A, :competitive)

        # Empty substrate set rejected.
        @test_throws ErrorException ER.EnzymeReaction(M.Reactant[], [P], noregs, Int[1])
        # Empty product set rejected.
        @test_throws ErrorException ER.EnzymeReaction([S], M.Reactant[], noregs, Int[1])
        # Duplicate product names rejected.
        @test_throws ErrorException ER.EnzymeReaction([S], [P, P], noregs, Int[1])
        # Atom imbalance rejected (S has 1 C, Punbal has 2 C).
        @test_throws ErrorException ER.EnzymeReaction([S], [Punbal], noregs, Int[1])
        # Two regulators of the same kind with the same name rejected.
        @test_throws ErrorException ER.EnzymeReaction([S], [P], [regA, regA2], Int[1])
        # A regulator MAY share a substrate/product name (`::Mod` role tag:
        # one metabolite binds as a competitive inhibitor under its real name).
        rxn_inh = ER.EnzymeReaction([S], [P], [regS], Int[1])
        @test :S in [s.name for s in rxn_inh.substrates]
        # Valid balanced reaction with a distinct-named regulator still constructs.
        rxn = ER.EnzymeReaction([S], [P], [regA], Int[1])
        @test [s.name for s in rxn.substrates] == [:S]
        @test [p.name for p in rxn.products] == [:P]
        # A name listed as both a substrate and a product is rejected, whether it
        # is the only reactant name (X → X) or one of several (A + B → A + C).
        Xs = M.Reactant(:X, [:C => 2])
        @test_throws ErrorException(
            "EnzymeReaction: X named as both a substrate and a product; " *
            "concentrations and constants are keyed by name") ER.EnzymeReaction(
            [Xs], [Xs], noregs, Int[1])
        As = M.Reactant(:A, [:C => 1])
        Bs = M.Reactant(:B, [:N => 1])
        Cp = M.Reactant(:C, [:N => 1])
        @test_throws ErrorException(
            "EnzymeReaction: A named as both a substrate and a product; " *
            "concentrations and constants are keyed by name") ER.EnzymeReaction(
            [As, Bs], [As, Cp], noregs, Int[1])
        # A reactant or a regulator takes no name a rate equation keeps for its constants
        # and the local names of its code, and no name holding `_`, which constant names
        # use to join names.
        _testhelper_reserved(n) = ErrorException(
            "EnzymeReaction: $n is a reserved name; Keq, L, N, params, concs " *
            "and v are kept for a rate equation's constants and local names")
        for n in (:Keq, :L, :N, :params, :concs, :v)
            @test_throws _testhelper_reserved(n) ER.EnzymeReaction(
                [M.Reactant(n, [:C => 1])], [P], noregs, Int[1])
            @test_throws _testhelper_reserved(n) ER.EnzymeReaction(
                [S], [P], [M.Regulator(n, :allosteric)], Int[1])
        end
        @test_throws _testhelper_reserved(:Keq) @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: Keq
        end
        _testhelper_underscore(n) = ErrorException(
            "EnzymeReaction: the metabolite name $n holds `_`; constant names join " *
            "metabolite names with `_`, so a metabolite name holds none")
        @test_throws _testhelper_underscore(:L_Ala) ER.EnzymeReaction(
            [M.Reactant(:L_Ala, [:C => 1])], [P], noregs, Int[1])
        @test_throws _testhelper_underscore(:I_X) @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: I_X
        end
        # Nor is a metabolite named E, the free enzyme form's label, which would make the
        # enumerator's forms read as the metabolite; nor does it take a name that no DSL
        # text can write.
        free_E = ErrorException("EnzymeReaction: E cannot name a metabolite; it labels " *
                                "the free enzyme form")
        @test_throws free_E ER.EnzymeReaction(
            [M.Reactant(:E, [:C => 1])], [P], noregs, Int[1])
        @test_throws free_E @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: E
        end
        for n in (Symbol("1A"), Symbol("A B"), :end, :in)
            @test_throws(ErrorException("EnzymeReaction: \"$n\" cannot name a " *
                                        "metabolite; a metabolite name is an identifier " *
                                        "that is no keyword"),
                         ER.EnzymeReaction([S], [P], [M.Regulator(n, :competitive)],
                                           Int[1]))
        end
        # A contextual keyword is a keyword only in its own construct (`abstract type`,
        # `for outer`, `var"..."`), which no DSL text writes, so it names a metabolite and
        # the reaction's text reads back.
        for n in (:abstract, :mutable, :outer, :primitive, :public, :type, :var)
            r = ER.EnzymeReaction([M.Reactant(n, [:C => 1])], [P], noregs, Int[1])
            @test eval(Meta.parse(sprint(show, r))) == r
        end
        # Julia reads some identifiers as others: a micro sign as μ, a decomposed é as one
        # character. Such a name prints as text that builds another reaction, so it is
        # rejected; a name already in the form Julia reads, Greek or with a subscript,
        # reads back. The rejected strings are written with escapes: Julia would read the
        # characters themselves as the other names.
        for (n, read) in (("\u00b5M", "\u03bcM"), ("e\u0301", "\u00e9"))
            @test_throws(ErrorException("EnzymeReaction: Julia reads the name \"$n\" as " *
                                        "\"$read\", so no DSL text can write it; use " *
                                        "\"$read\""),
                         ER.EnzymeReaction([M.Reactant(Symbol(n), [:C => 1])], [P],
                                           noregs, Int[1]))
        end
        for n in (:α, :Δ, :μM, :S₁, :é)
            r = ER.EnzymeReaction([M.Reactant(n, [:C => 1])], [P], noregs, Int[1])
            @test eval(Meta.parse(sprint(show, r))) == r
        end
        # A single name MAY be declared in BOTH regulator roles: one allosteric
        # regulator and one competitive inhibitor. The two roles render to distinct
        # parameter names, so both are kept.
        regAllo = M.Regulator(:A, :allosteric)
        rxn_dual = ER.EnzymeReaction([S], [P], [regA, regAllo], Int[1])
        @test length(rxn_dual.regulators) == 2
        @test Set(reg.kind for reg in rxn_dual.regulators) ==
              Set([:allosteric, :competitive])
        @test all(reg.name == :A for reg in rxn_dual.regulators)
        # Two allosteric regulators of the same name still rejected.
        @test_throws ErrorException ER.EnzymeReaction(
            [S], [P], [regAllo, M.Regulator(:A, :allosteric)], Int[1])
    end

    @testset "EnzymeReaction shared_catalytic_site" begin
        M = ER.Mechanisms
        A = M.Reactant(:A, [:C => 1])
        B = M.Reactant(:B, [:C => 1])
        P = M.Reactant(:P, [:C => 1])
        Q = M.Reactant(:Q, [:C => 1])
        noregs = M.Regulator[]

        # Accepts a valid (substrate, product) pair; stored normalized + sorted.
        r = ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1];
            shared_catalytic_site = [(:A, :P)])
        @test r.shared_catalytic_site == [(:A, :P)]

        # Normalizes product-first input to (substrate, product).
        r2 = ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1];
            shared_catalytic_site = [(:P, :A)])
        @test r2.shared_catalytic_site == [(:A, :P)]

        # Sorts multiple pairs.
        r3 = ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1];
            shared_catalytic_site = [(:B, :Q), (:A, :P)])
        @test r3.shared_catalytic_site == [(:A, :P), (:B, :Q)]

        # Default is empty.
        @test ER.EnzymeReaction([A], [P], noregs, Int[1]).shared_catalytic_site ==
            Tuple{Symbol,Symbol}[]

        # Rejections.
        @test_throws ErrorException ER.EnzymeReaction(  # unknown name
            [A, B], [P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :Z)])
        @test_throws ErrorException ER.EnzymeReaction(  # two substrates
            [A, B], [P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :B)])
        @test_throws ErrorException ER.EnzymeReaction(  # two products
            [A, B], [P, Q], noregs, Int[1]; shared_catalytic_site = [(:P, :Q)])
        @test_throws ErrorException ER.EnzymeReaction(  # duplicate pair
            [A, B], [P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :P), (:P, :A)])

        # A metabolite that is BOTH a substrate and a competitive inhibitor may
        # still be the substrate side of a shared pair (validation keys on the
        # substrate/product role, not the inhibitor role).
        regA = M.Regulator(:A, :competitive)
        r_dual = ER.EnzymeReaction([A, B], [P, Q], [regA], Int[1];
            shared_catalytic_site = [(:A, :P)])
        @test r_dual.shared_catalytic_site == [(:A, :P)]

        # == / hash reflect the field.
        @test r != ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1])
        @test hash(r) != hash(ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1]))
        @test r == ER.EnzymeReaction([A, B], [P, Q], noregs, Int[1];
            shared_catalytic_site = [(:A, :P)])
    end

    @testset "Regulator kind and effect" begin
        M = ER.Mechanisms
        @test M.Regulator(:X, :allosteric).effect == :unspecified      # default
        @test M.Regulator(:X, :allosteric, :activator).effect == :activator
        @test M.Regulator(:X, :allosteric) != M.Regulator(:X, :allosteric, :activator)
        @test_throws(ErrorException("Regulator: kind must be :allosteric or " *
                                    ":competitive, got bogus"),
                     M.Regulator(:X, :bogus))
        @test_throws(ErrorException("Regulator: effect must be :activator, :inhibitor, " *
                                    "or :unspecified, got bogus"),
                     M.Regulator(:X, :allosteric, :bogus))
        # A competitive inhibitor takes no type: only an allosteric regulator activates or
        # inhibits through the conformation it favors.
        for t in (:activator, :inhibitor)
            @test_throws(
                ErrorException("Regulator: X is a competitive " *
                               "inhibitor, which takes no regulator type, got " *
                               ":$t; only an allosteric regulator is an " *
                               "activator (::Activator) or an inhibitor " *
                               "(::Inhibitor)"),
                M.Regulator(:X, :competitive, t))
        end
    end

    @testset "EnzymeReaction canonicalizes reactant + regulator ordering" begin
        M = ER.Mechanisms
        r1 = EnzymeReaction(
            [M.Reactant(:B, [:C => 1]), M.Reactant(:A, [:C => 1])],
            [M.Reactant(:P, [:C => 2])],
            [M.Regulator(:Y, :allosteric), M.Regulator(:X, :allosteric)],
            [3, 1, 2],
        )
        r2 = EnzymeReaction(
            [M.Reactant(:A, [:C => 1]), M.Reactant(:B, [:C => 1])],
            [M.Reactant(:P, [:C => 2])],
            [M.Regulator(:X, :allosteric), M.Regulator(:Y, :allosteric)],
            [1, 2, 3],
        )
        @test r1 == r2
        @test hash(r1) == hash(r2)
    end

    @testset "EnzymeReaction takes one or more oligomeric states, each ≥ 1" begin
        _testhelper_reaction(mults) =
            EnzymeReaction([ER.Mechanisms.Reactant(:S, [:C => 1])],
                           [ER.Mechanisms.Reactant(:P, [:C => 1])],
                           ER.Mechanisms.Regulator[], mults)
        @test_throws(ErrorException("EnzymeReaction: the allowed oligomeric states " *
                                    "must be a nonempty list of Ints ≥ 1, got [0, 1]"),
                     _testhelper_reaction([0, 1]))
        # An empty list would leave the search no multiplicity to assign.
        @test_throws(ErrorException("EnzymeReaction: the allowed oligomeric states " *
                                    "must be a nonempty list of Ints ≥ 1, got Int64[]"),
                     _testhelper_reaction(Int[]))
        # A multiplicity listed twice is no list a DSL text writes.
        @test_throws(ErrorException("EnzymeReaction: the allowed oligomeric states list " *
                                    "a multiplicity more than once; got [2, 2]."),
                     _testhelper_reaction([2, 2]))
        # A Bool is no multiplicity: the list is a Vector{Int}, and none is converted.
        @test_throws MethodError _testhelper_reaction(Bool[true])
    end

    @testset "the free enzyme's isomerization runs the way the reaction turns over" begin
        # The Segel Iso Uni Uni case: P leaves at F and A enters at E, so the cycle that
        # turns the reaction over runs F → E, and so does the isomerization, however it is
        # written.
        s2 = @enzyme_mechanism begin
            substrates: A; products: P
            steps: begin
                E + A <--> EA; EA <--> EP; EP <--> F + P; E <--> F  # last flipped
            end
        end
        forms = s2.catalytic.enzyme_forms
        iso = only(e for e in s2.catalytic.steps
                   if ER.Mechanisms._is_iso(e) && forms[e.from] in (:E, :F))
        @test forms[iso.from] === :F  # where P leaves
        @test forms[iso.to] === :E    # where A enters
    end

    @testset "Mechanism canonicalizes step order" begin
        (; bind, iso, rel) = _testhelper_uniuni_steps()

        m = _testhelper_uniuni_of([bind, iso, rel])
        c = m.catalytic
        @test (c.metabolites, c.substrates, c.products) == ([:S, :P], [1], [2])
        flat = ER.Derivation._flat_edges(c)
        @test Set(_testhelper_edge_text(m, k) for k in flat) ==
              Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])
        @test [c.steps[k].kinetic_group for k in flat] == [1, 2, 3]
        @test length(c.steps) == 3
        # Canonical by construction: a permuted input yields identical storage.
        m_perm = _testhelper_uniuni_of([rel, bind, iso])
        @test (c.enzyme_forms, c.steps) ==
              (m_perm.catalytic.enzyme_forms, m_perm.catalytic.steps)
        @test m == m_perm && hash(m) == hash(m_perm)
        # Steps and groups sort by their forms' places, the forms by name with free E
        # first: E_res_plusA_minusP sorts after EQ (`_` after the capitals), so its steps
        # come last.
        pingpong = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA <--> EP_res_plusA_minusP
                E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP
                E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
                EB_res_plusA_minusP ⇌ EQ
                E + Q ⇌ EQ
            end
        end
        @test [_testhelper_edge_text(pingpong, k)
               for k in ER.Derivation._flat_edges(pingpong.catalytic)] ==
              ["E + A ⇌ EA",
               "E + Q ⇌ EQ",
               "EA <--> EP_res_plusA_minusP",
               "EB_res_plusA_minusP ⇌ EQ",
               "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP",
               "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP"]
        # The forms: free E, then by name.
        @test pingpong.catalytic.enzyme_forms ==
              [:E, :EA, :EB_res_plusA_minusP, :EP_res_plusA_minusP, :EQ,
               :E_res_plusA_minusP]
    end

    @testset "Allosteric mechanism (non-parametric)" begin
        (; bind, iso, rel) = _testhelper_uniuni_steps()
        sites = [ER.Mechanisms.RegulatorySite([:I], 1, [:OnlyI])]
        m = _testhelper_tagged_uniuni([bind, iso, rel], [:EqualAI, :OnlyA, :NonequalAI], 2,
                                      sites)

        @test (m.catalytic.metabolites, m.catalytic.substrates, m.catalytic.products) ==
              ([:S, :P], [1], [2])
        # Steps canonicalized; allosteric tags stay bound to their steps.
        @test Set(_testhelper_edge_text(m, k) for k in eachindex(m.catalytic.steps)) ==
              Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])
        _testhelper_state_of(text) =
            m.allosteric.allosteric_states[
                m.catalytic.steps[_testhelper_edge(m, text)].kinetic_group]
        @test _testhelper_state_of("E + S ⇌ ES") == :EqualAI
        @test _testhelper_state_of("ES <--> EP")  == :OnlyA
        @test _testhelper_state_of("E + P ⇌ EP")  == :NonequalAI
        @test m.allosteric.oligomeric_state == 2
        @test m.allosteric.regulatory_sites ==
              [ER.Mechanisms.RegulatorySite([:I], 1, [:OnlyI])]
        @test length(_testhelper_group_members(m.catalytic)) == 3

        # Written in another order, each tag beside its step.
        m2 = _testhelper_tagged_uniuni([rel, iso, bind], [:NonequalAI, :OnlyA, :EqualAI], 2,
                                       sites)
        @test m == m2
        @test hash(m) == hash(m2)
    end

    @testset "Allosteric mechanism validation errors" begin
        steps = collect(_testhelper_uniuni_steps())
        _testhelper_invalid_tag(tag) = ErrorException(
            "Mechanism: catalytic group 1 has invalid allosteric state $tag (must be " *
            "one of (:OnlyA, :OnlyI, :EqualAI, :NonequalAI))")

        # The constructor takes :OnlyI on any group and leaves the verdict to the guard:
        # the active conformation cannot switch off the rapid-equilibrium S binding,
        # whose two sides the chemistry and the P binding still join. Neither fix the
        # message offers passes: each leaves the active conformation no turnover.
        # Untagging the binding passes, but is not one of them.
        uni = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                EP ⇌ E + P
            end
        end
        k = _testhelper_edge(uni, "E + S ⇌ ES")
        onlyi_tags = fill(:EqualAI, length(_testhelper_group_members(uni.catalytic)))
        onlyi_tags[uni.catalytic.steps[k].kinetic_group] = :OnlyI
        onlyi = ER.Mechanism(uni.catalytic, ER.Mechanisms.Allostery(
            onlyi_tags, 1, ER.Mechanisms.RegulatorySite[]))
        v = ER.Derivation._guard(onlyi)
        @test (v.fault, v.conformation, v.edge) == (ER.Derivation._FAST_STEP_PINNED, 1, k)
        @test ER.Derivation._guard_message(onlyi) ==
              _testhelper_pinned_cause(onlyi, k, 1)
        @test all(t -> ER.Derivation._retag_passes(onlyi, [k], t), (:EqualAI, :NonequalAI))

        # Length mismatch
        @test_throws(ErrorException("Mechanism: allosteric_states length 2 must " *
                                    "match the number of kinetic groups 3"),
                     _testhelper_tagged_uniuni(steps, [:EqualAI, :NonequalAI], 1,
                                               ER.Mechanisms.RegulatorySite[]))

        # oligomeric_state < 1
        @test_throws(ErrorException("Mechanism: oligomeric_state must be " *
                                    "≥ 1, got 0"),
                     _testhelper_tagged_uniuni(steps, [:EqualAI, :EqualAI, :EqualAI], 0,
                                               ER.Mechanisms.RegulatorySite[]))

        # Unknown allosteric state symbol
        @test_throws(_testhelper_invalid_tag(:Bogus),
                     _testhelper_tagged_uniuni(steps, [:Bogus, :EqualAI, :EqualAI], 1,
                                               ER.Mechanisms.RegulatorySite[]))
    end

    @testset "an allosteric mechanism holds its catalytic mechanism, states and sites" begin
        am = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: NonequalAI
                ES <--> EP    :: OnlyA
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(1): A::OnlyA, B::NonequalAI
        end
        @test _testhelper_is_allosteric(am)

        # Its catalytic scheme is the plain mechanism's of the same steps.
        @test am.catalytic == _testhelper_re_mm.catalytic
        re_mm = _testhelper_re_mm.catalytic
        @test (am.catalytic.metabolites, am.catalytic.substrates, am.catalytic.products) ==
              (re_mm.metabolites, re_mm.substrates, re_mm.products)

        # allosteric_states and multiplicity extracted from CS
        @test am.allosteric.allosteric_states[1] === :EqualAI
        @test am.allosteric.allosteric_states[2] === :NonequalAI
        @test am.allosteric.allosteric_states[3] === :OnlyA
        @test am.allosteric.oligomeric_state == 2

        # Regulatory sites: ligand names
        sites = am.allosteric.regulatory_sites
        @test length(sites) == 1
        @test sites[1].ligands == [:A, :B]
        @test sites[1].multiplicity == 1
        @test sites[1].allosteric_states == [:OnlyA, :NonequalAI]
    end
end

# The exception that `f()` throws, or `nothing` when it returns.
function _testhelper_thrown(f)
    try
        f()
    catch e
        return e
    end
    nothing
end

@testset "kinetic groups hold one kind of step with one flag" begin
    @testset "binding and Theorell-Chance step in one group" begin
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> EA
                    (EA + B <--> EAB, EA + B <--> EQ + P)
                    EQ <--> E + Q
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds EA + B <--> EAB and EA + B <--> EQ + P",
                       err.msg)
    end

    @testset "RE and SS steps in one group" begin
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    (E + S ⇌ ES, EP + S <--> EPS)
                    ES <--> EP
                    EP ⇌ E + P
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds E + S ⇌ ES and EP + S <--> EPS", err.msg)
    end

    @testset "one reaction written RE and SS in one group" begin
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    (E + S ⇌ ES, E + S <--> ES)
                    ES <--> EP
                    EP ⇌ E + P
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds E + S <--> ES and E + S ⇌ ES", err.msg)
    end

    @testset "bindings of two different metabolites in one group" begin
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: A, B
                products: P
                steps: begin
                    (E + A <--> EA, EA + B <--> EAB)
                    EAB <--> E + P
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds E + A <--> EA and EA + B <--> EAB", err.msg)
    end

    @testset "accepted: context-shared bindings of one metabolite" begin
        # R binds both free E and ES with the same K (non-competitive
        # inhibitor pattern): a legitimate shared binding of one metabolite.
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: R
            steps: begin
                E + S ⇌ ES
                ES <--> EP
                E + P ⇌ EP
                (E + R ⇌ ER, ES + R ⇌ ERS)
                ER + S ⇌ ERS
            end
        end
        @test !_testhelper_is_allosteric(m)
    end

    @testset "accepted: mirrored isomerizations" begin
        # Two different conformational isomerizations, free and A-bound,
        # sharing one rate by a symmetry assumption: both are `:iso`, so
        # they are the same kind and may share a kinetic group.
        m = @enzyme_mechanism begin
            substrates: A
            products: P
            steps: begin
                E + A ⇌ EA
                (E <--> Estar, EA <--> EstarA)
                EstarA <--> EstarP
                EstarP ⇌ Estar + P
            end
        end
        @test !_testhelper_is_allosteric(m)
    end
end

@testset "Allosteric mechanism enforces the kinetic-group rules" begin
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: NonequalAI
            ES <--> EP   :: EqualAI
            EP ⇌ E + P     :: EqualAI
        end
    end
    cat = am.catalytic
    step_of(pred) = only(s for s in cat.steps if pred(s))
    bind_s = step_of(s -> s.taken_up == [_testhelper_metabolite(cat, :S)])
    bind_p = step_of(s -> s.taken_up == [_testhelper_metabolite(cat, :P)])
    iso = step_of(ER.Mechanisms._is_iso)
    # The catalytic steps `steps` between the catalytic nodes, step `k` labelled
    # `group[k]`, every group :EqualAI.
    _testhelper_rebuild(steps, group) = _testhelper_thrown() do
        ER.Mechanism(
            ER.Mechanisms.CatalyticScheme(cat.metabolites, cat.enzyme_forms,
                                          cat.substrates, cat.products,
            [ER.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off,
                                e.rapid_equilibrium, g) for (e, g) in zip(steps, group)]),
            ER.Mechanisms.Allostery(fill(:EqualAI, maximum(group)), 2,
                                    ER.Mechanisms.RegulatorySite[]))
    end
    # The two bindings in one group: different metabolites.
    err = _testhelper_rebuild([bind_s, bind_p, iso], [1, 1, 2])
    @test err isa ErrorException
    @test occursin("a kinetic group holds E + P ⇌ EP and E + S ⇌ ES", err.msg)
    # The isomerization in two groups.
    err = _testhelper_rebuild([bind_s, iso, iso, bind_p], [1, 2, 3, 4])
    @test err isa ErrorException
    @test occursin("both hold the reaction ES <--> EP", err.msg)
    # The isomerization twice in one group.
    err = _testhelper_rebuild([bind_s, iso, iso, bind_p], [1, 2, 2, 3])
    @test err isa ErrorException
    @test occursin("holds the reaction ES <--> EP twice", err.msg)
end

@testset "a competitive-inhibitor copy marks no reactant side" begin
    # P's competitive-inhibitor copy binds E. Its binding and release turn nothing over, so
    # the Iso Uni Uni isomerization still runs the way the reaction does, from F, where P
    # leaves, to E, where A enters.
    m = @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            E + P ⇌ EPmod
        end
    end
    forms = m.catalytic.enzyme_forms
    iso = only(e for e in m.catalytic.steps
               if ER.Mechanisms._is_iso(e) && forms[e.from] in (:E, :F))
    @test forms[iso.from] === :F && forms[iso.to] === :E
end

@testset "the stored direction of a step" begin
    stored(m) = Set(ER.Mechanisms.step_text(m.catalytic, k)
                    for k in eachindex(m.catalytic.steps))
    # The isomerization of the free enzyme runs F → E on the one cycle that turns the
    # reaction over, E → EA → EP → F → E, however it is written.
    iso = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> EA
            EA <--> EP
            F + P <--> EP
            F <--> E
        end
    end
    @test stored(iso) ==
          Set(["E + A <--> EA", "EA <--> EP", "F + P <--> EP", "F <--> E"])
    @test (@_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> EA
            EA <--> EP
            F + P <--> EP
            E <--> F
        end
    end) == iso
    # Two isomerizations of the free enzyme in a row, G → F → E, neither of whose forms
    # takes up or gives off a metabolite at F: the cycle orients both.
    chain = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> EA
            EA <--> EP
            G + P <--> EP
            F <--> G
            E <--> F
        end
    end
    @test stored(chain) == Set(["E + A <--> EA", "EA <--> EP", "G + P <--> EP",
                                "G <--> F", "F <--> E"])
    # A step that takes up one metabolite and gives off another runs the way the reaction
    # runs: ES + P ⇌ EP + S gives off the substrate and takes up the product, so it is
    # stored taking up S.
    exchange = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            ES + P ⇌ EP + S
        end
    end
    @test "EP + S ⇌ ES + P" in stored(exchange)
    # A dead-end isomerization lies on no cycle, so it runs away from free E, from the
    # form fewer steps from it, whichever name sorts first.
    dead_end = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            E + I ⇌ EI
            EstarI <--> EI
        end
    end
    @test "EI <--> EstarI" in stored(dead_end)
    dead_end_d = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            E + I ⇌ EI
            DI <--> EI
        end
    end
    @test "EI <--> DI" in stored(dead_end_d)
    # The enzyme's isomerization at E and at EA shares one pair of constants, and the
    # reaction runs through the two steps in opposite directions. The group turns with
    # Estar → E, its member nearest free E, and the other member keeps the direction it was
    # written in relative to it: writing both reversed changes nothing, and writing one
    # reversed writes another mechanism.
    mirrored = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ EA
            (E <--> Estar, EA <--> EstarA)
            EstarA <--> EstarP
            EstarP ⇌ Estar + P
        end
    end
    both_reversed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ EA
            (Estar <--> E, EstarA <--> EA)
            EstarA <--> EstarP
            EstarP ⇌ Estar + P
        end
    end
    one_reversed = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A ⇌ EA
            (E <--> Estar, EstarA <--> EA)
            EstarA <--> EstarP
            EstarP ⇌ Estar + P
        end
    end
    @test both_reversed == mirrored && hash(both_reversed) == hash(mirrored)
    @test one_reversed != mirrored
    @test issubset(["Estar <--> E", "EstarA <--> EA"], stored(mirrored))
    @test issubset(["Estar <--> E", "EA <--> EstarA"], stored(one_reversed))
    # The last tie is broken by the forms' names, free E first. Eb and Eab are both one
    # step from free E, and the isomerization between them lies on no cycle that turns
    # the reaction over: it runs from Eab, which comes before Eb, however it is written.
    by_name = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: b
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            E + b ⇌ Eb
            E + b ⇌ Eab
            Eb <--> Eab
        end
    end
    @test "Eab <--> Eb" in stored(by_name)
    @test (@_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: b
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            E + b ⇌ Eb
            E + b ⇌ Eab
            Eab <--> Eb
        end
    end) == by_name
    # A group's lead, among members equally near free E, is the one whose forms come first
    # by name. Neither isomerization below lies on a cycle that turns the reaction over,
    # and all four forms are one step from free E. Ea <--> Fa leads, named before
    # Ez <--> Fz; it runs from Ea, and Ez <--> Fz, written in the same direction as
    # Fa <--> Ea, runs from Fz.
    lead_by_name = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: z
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            E + z ⇌ Ez
            E + z ⇌ Fz
            E <--> Ea
            E <--> Fa
            (Ez <--> Fz, Fa <--> Ea)
        end
    end
    @test issubset(["Ea <--> Fa", "Fz <--> Ez"], stored(lead_by_name))
    @test (@_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: z
        steps: begin
            E + S <--> ES
            ES <--> EP
            EP <--> E + P
            E + z ⇌ Ez
            E + z ⇌ Fz
            E <--> Ea
            E <--> Fa
            (Ea <--> Fa, Fz <--> Ez)
        end
    end) == lead_by_name
end

@testset "the tied steps of one kinetic group turn as one" begin
    # Segel Iso Uni Uni with an inhibitor I binding E and F in one group and the
    # isomerization mirrored at the I-bound forms. The group turns with F → E, its member
    # nearest free E, which runs the way the reaction turns over, and the mirror keeps the
    # direction it is written in relative to F → E. Both run F-side to E-side, so the cycle
    # E → EI → FI → F → E constrains nothing and no line forces k_F_to_E = k_E_to_F.
    ss = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            (F <--> E, FI <--> EI)
            (E + I ⇌ EI, F + I ⇌ FI)
        end
    end
    # The conformation a form name gives: its first letter.
    label(f) = first(String(ss.catalytic.enzyme_forms[f]))
    mirror = only(ss.catalytic.steps[ks] for ks in _testhelper_group_members(ss.catalytic)
                  if length(ks) == 2 &&
                     all(k -> ER.Mechanisms._is_iso(ss.catalytic.steps[k]), ks))
    @test all(e -> label(e.from) == 'F' && label(e.to) == 'E', mirror)
    @test !occursin("Wegscheider", rate_equation_string(ss))
    # Writing both mirrored steps the other way round, in the other order, changes nothing.
    @test (@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            (EI <--> FI, E <--> F)
            (E + I ⇌ EI, F + I ⇌ FI)
        end
    end) == ss

    # At rapid equilibrium the same cycle would pin K_F_to_E to 1 on its own; the shared
    # orientation leaves it fitted.
    re = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            (F ⇌ E, FI ⇌ EI)
            (E + I ⇌ EI, F + I ⇌ FI)
        end
    end
    @test :K_F_to_E in parameters(re)
    # A group whose isomerizations join different pairs of conformations is a mechanism:
    # its members turn with the lead as written.
    @test !_testhelper_is_allosteric(@enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> E + P
            (E ⇌ F, F ⇌ G, G ⇌ E)
        end
    end)
end

# The edges of the one kinetic group of `mech`'s catalytic graph that neither binds nor
# isomerizes, and the names of what each takes up and gives off.
function _testhelper_exchanges(mech)
    cm = mech.catalytic
    group = only(ks for ks in _testhelper_group_members(cm) if !any(
        k -> ER.Mechanisms._is_iso(cm.steps[k]) || ER.Enumeration._is_binding(cm.steps[k]),
        ks))
    [(cm.metabolites[e.taken_up], cm.metabolites[e.given_off]) for e in cm.steps[group]]
end

@testset "a kinetic group ties one metabolite exchange at two conformations" begin
    # Segel Iso Uni Uni with an inhibitor I binding E and F, and J displacing I at both
    # conformations in one group. Neither exchange lies on a cycle that turns the reaction
    # over. The one at E leads, nearer free E, and both its forms are one step from free E,
    # so it runs from EI, whose name comes first, taking up J; the one at F takes up J with
    # it, however the steps are written. The allosteric twin joins F to E at rapid
    # equilibrium.
    # J binds E and F as well: without those bindings EJ and FJ are dead ends at zero
    # modifiers, a second closed class the guard rejects.
    I, J = :I, :J
    plain = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            (E + I ⇌ EI, F + I ⇌ FI)
            (E + J ⇌ EJ, F + J ⇌ FJ)
            (EI + J <--> EJ + I, FI + J <--> FJ + I)
        end
    end
    allo = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA                                   :: EqualAI
            EA <--> EP                                    :: EqualAI
            EP <--> F + P                                   :: EqualAI
            F ⇌ E                                             :: EqualAI
            (E + I ⇌ EI, F + I ⇌ FI)                      :: EqualAI
            (E + J ⇌ EJ, F + J ⇌ FJ)                      :: EqualAI
            (EI + J <--> EJ + I, FI + J <--> FJ + I)  :: NonequalAI
        end
    end
    m, am = plain, allo
    for mech in (m, am)
        @test _testhelper_exchanges(mech) == [([J], [I]), ([J], [I])]
        _testhelper_check_stored_orientation(mech)
    end
    # The written direction and order of the exchanges do not matter.
    @test (@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            (E + I ⇌ EI, F + I ⇌ FI)
            (E + J ⇌ EJ, F + J ⇌ FJ)
            (FJ + I <--> FI + J, EJ + I <--> EI + J)
        end
    end) == m
    @test (@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA                                   :: EqualAI
            EA <--> EP                                    :: EqualAI
            EP <--> F + P                                   :: EqualAI
            F ⇌ E                                             :: EqualAI
            (E + I ⇌ EI, F + I ⇌ FI)                      :: EqualAI
            (E + J ⇌ EJ, F + J ⇌ FJ)                      :: EqualAI
            (FJ + I <--> FI + J, EJ + I <--> EI + J)  :: NonequalAI
        end
    end) == am
    # The rate law derives, with the exchange constants named in the J-in direction.
    @test occursin("k_EI_J_to_EJ_I", rate_equation_string(plain))
    allo_law = rate_equation_string(allo)
    @test occursin("k_A_EI_J_to_EJ_I", allo_law)
    @test occursin("k_I_EI_J_to_EJ_I", allo_law)
    _testhelper_check_law_against_oracle(plain)
    _testhelper_check_law_against_oracle(allo)

    # Exchanges that switch conformation in opposite directions as J displaces I build;
    # without the J bindings the guard rejects them at zero I and J, for EJ.
    @test !_testhelper_is_allosteric(@_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            (E + I ⇌ EI, F + I ⇌ FI)
            (EI + J <--> FJ + I, FI + J <--> EJ + I)
        end
    end)
    # Without the J bindings the guard rejects the mechanism at zero modifiers.
    dead_ends = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            (E + I ⇌ EI, F + I ⇌ FI)
            (EI + J <--> EJ + I, FI + J <--> FJ + I)
        end
    end
    v = EnzymeRates.Derivation._guard(dead_ends)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._NONUNIQUE, 4)
    @test _testhelper_respeed_verdicts(dead_ends; rng = Random.MersenneTwister(81),
                                       law = :limit)[4] === :nonunique
end

@testset "a tied exchange turns to take up its lead's metabolites" begin
    # The exchange fixture above plus A binding EI and P binding EJ. The direction is
    # read from the steps alone, and the new bindings close no cycle through either
    # exchange, so the group still runs from EI and FI, taking up J.
    # J binds E and F as well: without those bindings EJ and FJ are dead ends at zero
    # modifiers, a second closed class the guard rejects.
    I, J = :I, :J
    plain = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            (E + I ⇌ EI, F + I ⇌ FI)
            (E + J ⇌ EJ, F + J ⇌ FJ)
            EI + A ⇌ EAI
            EJ + P ⇌ EJP
            (EI + J <--> EJ + I, FI + J <--> FJ + I)
        end
    end
    allo = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA                                   :: EqualAI
            EA <--> EP                                    :: EqualAI
            EP <--> F + P                                   :: EqualAI
            F ⇌ E                                             :: EqualAI
            (E + I ⇌ EI, F + I ⇌ FI)                      :: EqualAI
            (E + J ⇌ EJ, F + J ⇌ FJ)                      :: EqualAI
            EI + A ⇌ EAI                                :: EqualAI
            EJ + P ⇌ EJP                                :: EqualAI
            (EI + J <--> EJ + I, FI + J <--> FJ + I)  :: NonequalAI
        end
    end
    m, am = plain, allo
    for mech in (m, am)
        @test _testhelper_exchanges(mech) == [([J], [I]), ([J], [I])]
        _testhelper_check_stored_orientation(mech)
    end
    # The written direction and order of the exchanges do not matter.
    @test (@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            (E + I ⇌ EI, F + I ⇌ FI)
            (E + J ⇌ EJ, F + J ⇌ FJ)
            EI + A ⇌ EAI
            EJ + P ⇌ EJP
            (FJ + I <--> FI + J, EJ + I <--> EI + J)
        end
    end) == m
    @test (@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA                                   :: EqualAI
            EA <--> EP                                    :: EqualAI
            EP <--> F + P                                   :: EqualAI
            F ⇌ E                                             :: EqualAI
            (E + I ⇌ EI, F + I ⇌ FI)                      :: EqualAI
            (E + J ⇌ EJ, F + J ⇌ FJ)                      :: EqualAI
            EI + A ⇌ EAI                                :: EqualAI
            EJ + P ⇌ EJP                                :: EqualAI
            (FJ + I <--> FI + J, EJ + I <--> EI + J)  :: NonequalAI
        end
    end) == am
    @test occursin("k_EI_J_to_EJ_I", rate_equation_string(plain))
    allo_law = rate_equation_string(allo)
    @test occursin("k_A_EI_J_to_EJ_I", allo_law)
    @test occursin("k_I_EI_J_to_EJ_I", allo_law)
    _testhelper_check_law_against_oracle(plain)
    _testhelper_check_law_against_oracle(allo)

    # The exchange between conformations, FI + J <--> EJ + I, lies on a cycle that
    # turns the reaction over, FI → EJ → E → EA → EP → F → FI, so it runs from
    # FI, taking up J, and leads its group; the exchange at E, on no such cycle, takes up
    # J with it.
    mixed = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            (E + I ⇌ EI, F + I ⇌ FI)
            E + J ⇌ EJ
            (EI + J <--> EJ + I, FI + J <--> EJ + I)
        end
    end
    @test _testhelper_exchanges(mixed) == [([J], [I]), ([J], [I])]
    _testhelper_check_stored_orientation(mixed)
    @test (@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            (E + I ⇌ EI, F + I ⇌ FI)
            E + J ⇌ EJ
            (EJ + I <--> FI + J, EJ + I <--> EI + J)
        end
    end) == mixed
    @test occursin("k_EI_J_to_EJ_I", rate_equation_string(mixed))
    _testhelper_check_law_against_oracle(mixed)

    # Neither exchange lies on a cycle, so each runs away from free E: on its own, the
    # exchange at F runs from FJ, two steps from free E, to FI, three steps away,
    # taking up I. In one group with the exchange at E, which leads, nearer free E, and
    # runs from EI, taking up J, it turns to take up J too.
    _testhelper_lone_f(ff) = only(e for e in ff.catalytic.steps
                                  if !isempty(e.taken_up) && !isempty(e.given_off) &&
                                      startswith(String(ff.catalytic.enzyme_forms[e.from]),
                                                 "F"))
    apart = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            E + I ⇌ EI
            F + J ⇌ FJ
            EI + J <--> EJ + I
            FI + J <--> FJ + I
        end
    end
    lone = _testhelper_lone_f(apart)
    @test (apart.catalytic.metabolites[lone.taken_up],
           apart.catalytic.metabolites[lone.given_off]) == ([I], [J])
    tied = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            F <--> E
            E + I ⇌ EI
            F + J ⇌ FJ
            (EI + J <--> EJ + I, FI + J <--> FJ + I)
        end
    end
    @test _testhelper_exchanges(tied) == [([J], [I]), ([J], [I])]
    _testhelper_check_stored_orientation(tied)
end

@testset "each reaction appears once in a mechanism" begin
    # The mechanism over A and P whose steps, each its own group, are `lines`, built and
    # checked by the constructor.
    _testhelper_written(lines...) = ER.DSL._construct(ER.DSL.raw_mechanism(:(begin
        substrates: A
        products: P
        steps: $(Expr(:block, lines...))
    end)))
    # A binds E into two different forms from two groups: two reactions, whose
    # names differ in the bound form (k_E_A_to_EA, k_E_A_to_EstarA).
    m = _testhelper_written(:(E + A <--> EA), :(E + A <--> EstarA),
                            :(EA <--> EP), :(E + P <--> EP))
    names = [ER.Derivation._name(p, m) for members in _testhelper_group_members(m.catalytic)
             for p in _testhelper_step_parameters(m, first(members))]
    @test allunique(names)
    @test issubset([:k_E_A_to_EA, :k_E_A_to_EstarA], names)
    # The same binding in two groups, once written as its release: both would be
    # k_E_A_to_EA.
    err = _testhelper_thrown() do
        _testhelper_written(:(E + A <--> EA), :(EA <--> E + A), :(EA <--> EP),
                            :(E + P <--> EP))
    end
    @test err isa ErrorException
    @test occursin("both hold the reaction E + A <--> EA", err.msg)
    # The same isomerization in two groups: both would be k_EA_to_EP.
    err = _testhelper_thrown() do
        _testhelper_written(:(E + A <--> EA), :(EA <--> EP), :(EP <--> EA),
                            :(E + P <--> EP))
    end
    @test err isa ErrorException
    @test occursin("both hold the reaction EA <--> EP", err.msg)
    # The same binding written RE in one group and SS in another: one reaction
    # cannot be both fast and slow.
    err = _testhelper_thrown() do
        @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                E + S ⇌ ES
                ES ⇌ EP
                EP ⇌ E + P
            end
        end
    end
    @test err isa ErrorException
    @test occursin("both hold the reaction E + S <--> ES (SS in group 2, RE in group 3)",
                   err.msg)
    # The same binding twice in one group, once written as its release: the
    # derivation would count its edge twice.
    err = _testhelper_thrown() do
        @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                (E + S <--> ES, ES <--> E + S)
                ES <--> EP
                EP <--> E + P
            end
        end
    end
    @test err isa ErrorException
    @test occursin("kinetic group 2 holds the reaction E + S <--> ES twice", err.msg)
end

@testset "EI and EIinh, EPmod and EPinh name their constants apart" begin
    # EI and EIinh render EI and EIinh: the mechanism is accepted, and each constant
    # takes a name of its own.
    iinh = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        modifiers: I
        steps: begin
            E + Iinh ⇌ EIinh
            EIinh <--> EP
            EP ⇌ E + P
            E + I ⇌ EI
        end
    end
    @test _testhelper_step_constant_names(iinh) ==
          Set([:K_EIinh_to_E_Iinh, :k_EIinh_to_EP, :k_EP_to_EIinh, :K_EP_to_E_P,
               :K_EI_to_E_I])
    _testhelper_check_law_against_oracle(iinh)
    # So with a reactant's modifier copy: EPmod and EPinh render EPmod and EPinh.
    pinh = @enzyme_mechanism begin
        substrates: Pinh
        products: P
        steps: begin
            E + Pinh ⇌ EPinh
            EPinh <--> EP
            EP ⇌ E + P
            E + P ⇌ EPmod
        end
    end
    @test _testhelper_step_constant_names(pinh) ==
          Set([:K_EPinh_to_E_Pinh, :k_EPinh_to_EP, :k_EP_to_EPinh, :K_EP_to_E_P,
               :K_EPmod_to_E_P])
    _testhelper_check_law_against_oracle(pinh)
end

@testset "admission rejects a one-sided :OnlyA binding and admits balanced ones" begin
    # Admission's fault for `am`'s graph tagged `tags`, and the metabolite of the
    # binding it names (`nothing` when it names none).
    function _testhelper_admission_reads(am, tags)
        cm = am.catalytic
        tagged = ER.Mechanism(cm, ER.Mechanisms.Allostery(tags, 1,
                                                          ER.Mechanisms.RegulatorySite[]))
        v = ER.Derivation._verdict_under(tagged, ER.Derivation._built_reading(tagged))
        e = v.edge == 0 ? nothing : cm.steps[v.edge]
        (v.fault, e === nothing || !ER.Enumeration._is_binding(e) ? nothing :
                  cm.metabolites[only(e.taken_up)])
    end
    NO, FP = ER.Derivation._NO_FAULT, ER.Derivation._FAST_STEP_PINNED

    # Uni-uni S -> P. Tags: (S binding, chemical step, P binding). Admission's fault,
    # and the metabolite of the binding it names. Each binding is at rapid equilibrium,
    # so the inactive conformation can switch off an :OnlyA binding only by its affinity,
    # which needs its two sides apart.
    function _testhelper_uni(s_tag, cat_tag, p_tag)
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            oligomeric_state: 2
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: EqualAI
                EP ⇌ E + P      :: EqualAI
            end
        end
        firsts = _testhelper_first_steps(_testhelper_graph(m))
        # tag each group by the metabolite its step binds
        tags = [!_testhelper_is_binding(s) ? cat_tag :
                only(s.consumed) === :S ? s_tag : p_tag for s in firsts]
        _testhelper_admission_reads(m, tags)
    end

    # no :OnlyA anywhere -> admitted
    @test _testhelper_uni(:EqualAI, :EqualAI, :EqualAI) == (NO, nothing)
    # :OnlyA on the substrate only, catalysis :EqualAI: catalysis and the P binding
    # join E to ES in the inactive conformation -> the S binding cannot be switched off
    @test _testhelper_uni(:OnlyA, :EqualAI, :EqualAI) == (FP, :S)
    # :OnlyA on the product only, catalysis :EqualAI -> the P binding cannot be switched
    # off
    @test _testhelper_uni(:EqualAI, :EqualAI, :OnlyA) == (FP, :P)
    # :OnlyA on the substrate and the catalysis: ES lies behind the binding
    @test _testhelper_uni(:OnlyA, :OnlyA, :EqualAI) == (NO, nothing)
    # balanced: :OnlyA on both bindings, catalysis :EqualAI -> admitted
    @test _testhelper_uni(:OnlyA, :EqualAI, :OnlyA) == (NO, nothing)
    # V-system: :OnlyA chemical step only -> admitted
    @test _testhelper_uni(:EqualAI, :OnlyA, :EqualAI) == (NO, nothing)
    # :NonequalAI catalysis joins E to ES too -> same verdict
    @test _testhelper_uni(:OnlyA, :NonequalAI, :EqualAI) == (FP, :S)

    # Whether a binding is rapid-equilibrium or steady-state does not change
    # admission's verdict: with both bindings :OnlyA, ES and EP lie behind them.
    # Both orientations of the mix are admitted.
    mixed_uni_re_ss = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: EqualAI
            ES <--> EP    :: EqualAI
            EP <--> E + P   :: EqualAI
        end
    end
    mixed_uni_ss_re = @enzyme_mechanism begin
        substrates: S
        products:   P
        oligomeric_state: 2
        steps: begin
            E + S <--> ES   :: EqualAI
            ES <--> EP    :: EqualAI
            EP ⇌ E + P      :: EqualAI
        end
    end
    # every binding :OnlyA, catalysis :EqualAI -> balanced -> valid
    function _testhelper_both_bindings_onlya(m)
        tags = [_testhelper_is_binding(s) ? :OnlyA : :EqualAI
                for s in _testhelper_first_steps(_testhelper_graph(m))]
        _testhelper_admission_reads(m, tags)
    end
    @test _testhelper_both_bindings_onlya(mixed_uni_re_ss) == (NO, nothing)
    @test _testhelper_both_bindings_onlya(mixed_uni_ss_re) == (NO, nothing)

    # A random-order binding square: E binds A and B, and each of EA and EB binds
    # the other. A lone :OnlyA edge has both sides joined by the other three, so the
    # inactive conformation cannot switch it off; tagging P's binding too separates them.
    @testset "random-order binding square" begin
        biuni = @enzyme_mechanism begin
            substrates: A, B
            products:   P
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EA          :: EqualAI
                E + B ⇌ EB          :: EqualAI
                EA + B ⇌ EAB    :: EqualAI
                EB + A ⇌ EAB    :: EqualAI
                EAB <--> EP     :: EqualAI
                EP ⇌ E + P          :: EqualAI
            end
        end
        # Tag :OnlyA the groups named by (free form, bound metabolite); every
        # other group is :EqualAI. The chemical step's key is (:EAB, nothing).
        bu = _testhelper_graph(biuni)
        # Each group's key: the form its first step starts from, and the metabolite it
        # binds (`nothing` for the chemical step).
        bu_keys = [(Symbol(bu.forms[s.from]),
                    _testhelper_is_binding(s) ? only(s.consumed) : nothing)
                   for s in _testhelper_first_steps(bu)]
        _testhelper_tags(onlya_keys...) =
            [key in onlya_keys ? :OnlyA : :EqualAI for key in bu_keys]
        _testhelper_verdict(t) = _testhelper_admission_reads(biuni, t)

        # the square alone trips nothing: no :OnlyA anywhere -> admitted
        @test _testhelper_verdict(_testhelper_tags()) == (NO, nothing)
        # one square edge :OnlyA: the square's other three bindings join its sides
        @test _testhelper_verdict(_testhelper_tags((:E, :A))) == (FP, :A)
        # ...and the macro rejects that tagging: untagged steps join the two sides of the
        # rapid-equilibrium :OnlyA binding of A, so the inactive conformation cannot
        # switch it off (`_FAST_STEP_PINNED`).
        lone = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products:   P
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EA          :: OnlyA
                E + B ⇌ EB          :: EqualAI
                EA + B ⇌ EAB    :: EqualAI
                EB + A ⇌ EAB    :: EqualAI
                EAB <--> EP     :: EqualAI
                EP ⇌ E + P          :: EqualAI
            end
        end
        v = ER.Derivation._guard(lone)
        @test v.fault === ER.Derivation._FAST_STEP_PINNED
        @test ER.Mechanisms.step_text(lone.catalytic, v.edge) == "E + A ⇌ EA"
        @test_throws ErrorException(
            ER.Derivation._guard_message(lone)) @enzyme_mechanism begin
            substrates: A, B
            products:   P
            oligomeric_state: 2
            steps: begin
                E + A ⇌ EA          :: OnlyA
                E + B ⇌ EB          :: EqualAI
                EA + B ⇌ EAB    :: EqualAI
                EB + A ⇌ EAB    :: EqualAI
                EAB <--> EP     :: EqualAI
                EP ⇌ E + P          :: EqualAI
            end
        end
        # both A-side bindings :OnlyA: the B bindings, the catalysis and the P binding
        # still join E to EA -> they cannot be switched off
        @test _testhelper_verdict(_testhelper_tags((:E, :A), (:EB, :A)))[1] === FP
        # ...and :OnlyA on P cuts the last join -> admitted
        @test _testhelper_verdict(
            _testhelper_tags((:E, :A), (:EB, :A), (:E, :P))) == (NO, nothing)
        # :OnlyA on the chemical step leaves the square, whose lone :OnlyA edge still
        # has its sides joined -> that edge cannot be switched off
        @test _testhelper_verdict(_testhelper_tags((:EAB, nothing), (:E, :A))) ==
              (FP, :A)
    end
end

@testset "admission rejects a ter-substrate cube joining a tagged binding's two sides" begin
    # A full random-order ter binding cube with seven :OnlyA edges and the cross bindings
    # :EqualAI. In the inactive conformation the untagged bindings join forms across the
    # cube, so a tagged binding's two sides share a component: no path takes its affinity to
    # zero. Admission finds untagged steps that join the two sides of an :OnlyA binding at
    # rapid equilibrium, which the inactive conformation cannot switch off
    # (`_FAST_STEP_PINNED`).
    cube = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B, C
        products:   P
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: EqualAI
            EB + A ⇌ EAB       :: EqualAI
            EB + C ⇌ EBC       :: EqualAI
            EC + A ⇌ EAC       :: EqualAI
            EC + B ⇌ EBC       :: EqualAI
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    v = ER.Derivation._guard(cube)
    @test v.fault === ER.Derivation._FAST_STEP_PINNED
    @test_throws ErrorException(
        ER.Derivation._guard_message(cube)) @enzyme_mechanism begin
        substrates: A, B, C
        products:   P
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: EqualAI
            EB + A ⇌ EAB       :: EqualAI
            EB + C ⇌ EBC       :: EqualAI
            EC + A ⇌ EAC       :: EqualAI
            EC + B ⇌ EBC       :: EqualAI
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
end

# The rejection testset above pins admission's reject direction on a ter-substrate
# binding cube; this one pins its accept direction there. An over-rejection does not
# fail loudly: it silently shrinks the enumeration's search space by dropping valid
# mechanisms before they are fitted.
@testset "admission admits feasible ter-substrate mechanisms" begin
    # A full random-order A/B/C binding lattice, every substrate binding and the
    # chemical step :OnlyA: in the inactive conformation every form stands alone but E
    # and EP, the bindings point away from free E and close no cycle, so admission
    # accepts the cube.
    cube = @enzyme_mechanism begin
        substrates: A, B, C
        products:   P
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA             :: OnlyA
            E + B ⇌ EB             :: OnlyA
            E + C ⇌ EC             :: OnlyA
            EA + B ⇌ EAB       :: OnlyA
            EA + C ⇌ EAC       :: OnlyA
            EB + A ⇌ EAB       :: OnlyA
            EB + C ⇌ EBC       :: OnlyA
            EC + A ⇌ EAC       :: OnlyA
            EC + B ⇌ EBC       :: OnlyA
            EAB + C ⇌ EABC :: OnlyA
            EAC + B ⇌ EABC :: OnlyA
            EBC + A ⇌ EABC :: OnlyA
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    @test _testhelper_is_allosteric(cube)

    # Ordered ter with an :OnlyA B binding and an :OnlyA chemical step: EAB and
    # EABC lie behind the B binding, on its bound side, so this is admitted.
    ordered = @enzyme_mechanism begin
        substrates: A, B, C
        products:   P
        oligomeric_state: 2
        steps: begin
            E + A ⇌ EA             :: EqualAI
            EA + B ⇌ EAB       :: OnlyA
            EAB + C ⇌ EABC :: EqualAI
            EABC <--> EP     :: OnlyA
            E + P ⇌ EP             :: EqualAI
        end
    end
    @test _testhelper_is_allosteric(ordered)
end

@testset "a ping-pong with both fused releases :OnlyA, pumped as first read, passes" begin
    # EA → E_res_plusA_minusP + P runs the first half-reaction and releases P in
    # one step, EB_res_plusA_minusP → E + Q the second and releases Q. In the inactive
    # conformation the covalent forms are cut off. Under the guard's first built reading
    # the second half, read as Q taken up into EB_res_plusA_minusP, loses its affinity,
    # and the first half, whose far side EA stays, is stopped. At zero A and B the B
    # binding feeds E_res_plusA_minusP, which leaves only through the vanishing first
    # half, so admission cannot show that its weight does not depend on how fast that
    # step vanishes (`_PUMPED`), and the ε-path oracle finds a limit that depends on the
    # direction there. The search admits the ping-pong under a reading in which the first
    # half loses its rate into the covalent forms instead, both oracles agreeing on every
    # face. Either half alone :OnlyA is admitted too.
    both = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        end
    end
    v = ER.Derivation._verdict_under(both, ER.Derivation._built_reading(both))
    @test (v.fault, v.face, v.conformation) == (ER.Derivation._PUMPED, 3, 2)
    @test _testhelper_epath_verdict(both; rng = Random.MersenneTwister(6),
                                    face = Set([:A, :B])) in
          (:direction_dependent, :wrong_limit)
    @test ER.Derivation._guard(both) == ER.Derivation._WELL_POSED
    @test all(_testhelper_admission_agreement(both; rng = Random.MersenneTwister(9),
                                              corners = true))
    _testhelper_check_law_against_limit(both)
    @test @enzyme_mechanism(begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        end
    end) == both
    half1 = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: OnlyA
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: EqualAI
        end
    end
    half2 = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + A <--> EA                                       :: EqualAI
            EA <--> E_res_plusA_minusP + P                   :: EqualAI
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP :: EqualAI
            EB_res_plusA_minusP <--> E + Q                     :: OnlyA
        end
    end
    # Each face: the RE-speed oracle (class 3) and the ε-path (admission) agree with the
    # guard, which accepts it.
    for (i, half) in enumerate((half1, half2))
        @test all(_testhelper_admission_agreement(half;
                                                  rng = Random.MersenneTwister(6 + i)))
    end
end

@testset "admission admits an abortive segment an :OnlyA group strands" begin
    # The A group mixes the fused binding into EPQ with A's abortive bindings at
    # EP and EQ. Tagged :OnlyA, it strands {EAP, EAQ, EAPQ} in the
    # inactive conformation: a rapid-equilibrium segment whose every form holds a
    # product, which a mechanism may not have. The product releases still make EPQ
    # there, so no built reading switches off the rapid-equilibrium fused binding. It
    # joins EB, cut off, to free E's group, and the first consistent reading lets it
    # lose its rate into EB; at zero A and B its uptake dies, and it runs one way into
    # EB, which empties only through the slow E + B <--> EB: EB fills faster than it
    # empties (`_FILLED`), as the ε-path oracle finds under that reading.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + B <--> EB                                              :: OnlyA
            (EB + A ⇌ EPQ, EP + A ⇌ EAP, EQ + A ⇌ EAQ) :: OnlyA
            EPQ <--> EP + Q                                        :: EqualAI
            EPQ <--> EQ + P                                        :: EqualAI
            (E + P ⇌ EP, EAQ + P ⇌ EAPQ)                     :: EqualAI
            (E + Q ⇌ EQ, EAP + Q ⇌ EAPQ)                     :: EqualAI
        end
    end
    k = _testhelper_edge(am, "EB + A ⇌ EPQ")
    @test ER.Derivation._admission_interior(am, ER.Derivation._built_reading(am)) ==
          ER.Derivation._Verdict(ER.Derivation._FAST_STEP_PINNED,
                                 ER.Derivation._LIMIT_UNTESTED, 0, 2, k)
    @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(3)) === :no_path
    v = ER.Derivation._guard(am)
    @test (v.fault, _testhelper_faces(am)[v.face], v.conformation, v.edge) ==
          (ER.Derivation._FILLED, Set([:A, :B]), 2, k)
    reported = _testhelper_reading(ER.Derivation._guard_with_reading(am).reading)
    @test _testhelper_admission_agrees(
        v.fault, _testhelper_epath_verdict(am; rng = Random.MersenneTwister(4),
                                           face = Set([:A, :B]), reading = reported))
    @test ER.Derivation._guard_message(am) ==
          "Mechanism: at zero A and B, in the inactive conformation the :OnlyA " *
          "rapid-equilibrium step EB + A ⇌ EPQ runs one way into EB, at a rate " *
          "that vanishes yet grows with the rapid-equilibrium speed, and EB leaves " *
          "only through the slower step E + B <--> EB, so EB fills faster than it " *
          "empties and its weight depends on how fast the tagged step vanishes."
    # With the fused binding steady-state in a group of its own, the inactive
    # conformation stops it, and admission reads each edge: the abortive bindings reach
    # {EAP, EAQ, EAPQ} from the product forms on their unbound side only,
    # so they lose their affinity, the inactive conformation holds no enzyme there, and
    # every face passes. Both oracles agree on every face.
    fused_ss = @enzyme_mechanism begin
        substrates: A, B ; products: P, Q ; oligomeric_state: 1
        steps: begin
            E + B <--> EB                                              :: OnlyA
            EB + A <--> EPQ                                        :: OnlyA
            (EP + A ⇌ EAP, EQ + A ⇌ EAQ)                     :: OnlyA
            EPQ <--> EP + Q                                        :: EqualAI
            EPQ <--> EQ + P                                        :: EqualAI
            (E + P ⇌ EP, EAQ + P ⇌ EAPQ)                     :: EqualAI
            (E + Q ⇌ EQ, EAP + Q ⇌ EAPQ)                     :: EqualAI
        end
    end
    @test all(_testhelper_admission_agreement(fused_ss; rng = Random.MersenneTwister(4)))
end

@testset "the rational nullspace helper" begin
    R = Rational{BigInt}

    @testset "_rational_nullspace" begin
        # Zero map: every x solves M * x = 0, so the nullspace is the whole
        # 3-dimensional space.
        M = zeros(R, 2, 3)
        N = ER.Derivation._rational_nullspace(M)
        @test size(N, 2) == 3

        # Full-rank square matrix: only x = 0 solves M * x = 0.
        M = R[1 2; 3 4]
        N = ER.Derivation._rational_nullspace(M)
        @test size(N, 2) == 0

        # Rank-1 2x2 (row 2 is twice row 1): both rows reduce to x + y = 0,
        # spanned by (1, -1). n = 2, rank = 1.
        M = R[1 1; 2 2]
        N = ER.Derivation._rational_nullspace(M)
        @test M * N == zeros(R, 2, size(N, 2))
        @test size(N, 2) == 2 - 1

        # More rows than columns, rank-deficient: every row is a multiple of
        # (1, 2), i.e. x + 2y = 0, spanned by (2, -1).
        M = R[1 2; 2 4; 3 6]
        N = ER.Derivation._rational_nullspace(M)
        @test M * N == zeros(R, 3, size(N, 2))
        @test size(N, 2) == 2 - 1

        # More columns than rows: x + z = 0 and y + z = 0, spanned by
        # (-1, -1, 1).
        M = R[1 0 1; 0 1 1]
        N = ER.Derivation._rational_nullspace(M)
        @test M * N == zeros(R, 2, size(N, 2))
        @test size(N, 2) == 3 - 2
    end
end

@testset "naming a mechanism's constants never changes its identity" begin
    m1 = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B ⇌ EAB
            EAB <--> EPQ
            EPQ ⇌ EQ + P
            EQ ⇌ E + Q
        end
    end
    m2 = _testhelper_build(_testhelper_written(m1))
    am1 = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: EqualAI
            ES <--> EP   :: OnlyA
            EP ⇌ E + P     :: EqualAI
        end
        regulatory_site(2): I::OnlyI
    end
    am2 = _testhelper_build(_testhelper_written(am1))
    # Every parameter's name, which reads the mechanism and stores nothing in it.
    _testhelper_param_names(m::EnzymeRates.Mechanism) =
        !_testhelper_is_allosteric(m) ? _testhelper_edge_names(m, :None) :
        [[_testhelper_edge_names(m, state) for state in (:A, :I)];
         [[K] for state in (:A, :I) for K in _testhelper_regulator_constants(m, state)]]
    # `rebuilt` is `m` constructed again from the nodes and steps it writes back.
    function _testhelper_check_identity(m, rebuilt)
        shown = repr(m)
        names = _testhelper_param_names(m)
        @test repr(m) == shown
        @test m == rebuilt && hash(m) == hash(rebuilt)
        @test compile_rate_equation(m) === compile_rate_equation(rebuilt)
        @test names == _testhelper_param_names(rebuilt)
        io = IOBuffer(); serialize(io, m); seekstart(io)
        copied = deserialize(io)
        @test copied == m && hash(copied) == hash(m)
        @test names == _testhelper_param_names(copied)
    end
    _testhelper_check_identity(m1, m2)
    _testhelper_check_identity(am1, am2)
end

@testset "a mechanism's identity leaves out its reaction" begin
    # A macro builds its mechanism from its block alone; the same steps enumerated from a
    # declared reaction, with real atoms and regulators no step binds, are the same
    # mechanism: `==` and `hash` read the graph, tags, multiplicity and sites only.
    rxn = @enzyme_reaction begin
        substrates: S[C6H12O6]
        products: P[C6H12O6]
        competitive_inhibitors: J
        allosteric_regulators: I
        oligomeric_state: 2
    end
    m_on_rxn = only(_testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn)))
    @test m_on_rxn == _testhelper_re_mm && hash(m_on_rxn) == hash(_testhelper_re_mm)
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES     :: EqualAI
            ES <--> EP   :: OnlyA
            EP ⇌ E + P     :: EqualAI
        end
        regulatory_site(2): I::OnlyI
    end
    am_on_rxn = EnzymeRates.Mechanism(
        m_on_rxn.catalytic,
        EnzymeRates.Mechanisms.Allostery(am.allosteric.allosteric_states,
                                         am.allosteric.oligomeric_state,
                                         am.allosteric.regulatory_sites))
    @test am_on_rxn == am && hash(am_on_rxn) == hash(am)
end

@testset "every derivation entry rejects a bottomless rapid-equilibrium segment" begin
    # Random product release at rapid equilibrium between steady-state releases leaves
    # the segment {EP, EPQ, EQ} without a bottom form, so the rate is 0/0 at zero
    # products. The macro's boundary check (`check_valid_rate_equation`) rejects the
    # mechanism.
    # The constructor checks the representation only, so it builds the mechanism from
    # the walker's raw value, and the mechanism compiles; every derivation entry then
    # rejects it with the same message.
    raw = EnzymeRates.DSL.raw_mechanism(:(begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ ES
            ES <--> EPQ
            EPQ ⇌ EQ + P
            EPQ ⇌ EP + Q
            EP <--> E + P
            EQ <--> E + Q
        end
    end))
    m = EnzymeRates.DSL._construct(raw)
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._BOTTOMLESS, 2)
    msg = ErrorException(EnzymeRates.Derivation._guard_message(m))
    @test_throws msg EnzymeRates.DSL._checked_mechanism(raw)
    rate = compile_rate_equation(m)
    @test_throws msg parameters(rate)
    @test_throws msg parameters(m)
    @test_throws msg rate((S = 1.0, P = 1.0, Q = 1.0), (;), 1.0)
    @test_throws msg rate_equation_string(m)
    @test_throws msg dependent_parameter_values(m, (;), 1.0)
    @test_throws msg EnzymeRates.Derivation.parameter_units(m)
end

"""A value of each type a mechanism is built from: a step, a catalytic scheme, a mechanism,
a site and an allosteric mechanism."""
function _testhelper_value_examples()
    R = EnzymeRates
    none = Int[]
    bind = R.Mechanisms.Step(1, 2, [1], none, true, 1)
    mech = _testhelper_plain_mechanism([:A, :P, :I], [:E, :EA, :EP, :EI], [1], [2],
                       [bind, R.Mechanisms.Step(2, 3, none, none, false, 2),
                        R.Mechanisms.Step(3, 1, none, [2], true, 3),
                        R.Mechanisms.Step(1, 4, [3], none, true, 4)])
    site = R.Mechanisms.RegulatorySite([:X], 2, [:OnlyI])
    allo = R.Mechanism(mech.catalytic,
                       R.Mechanisms.Allostery([:EqualAI, :EqualAI, :EqualAI, :OnlyA], 2,
                                              [site]))
    (bind, mech.catalytic, mech, site, allo)
end

@testset "mechanism values: content equality and hash" begin
    R = EnzymeRates
    a, none = 1, Int[]
    # Steps compare by content; a step's direction, flag and group take part.
    release = R.Mechanisms.Step(2, 1, none, [a], true, 1)
    @test release == R.Mechanisms.Step(2, 1, none, [a], true, 1)
    @test release != R.Mechanisms.Step(1, 2, [a], none, true, 1)
    @test release != R.Mechanisms.Step(2, 1, none, [a], false, 1)
    @test release != R.Mechanisms.Step(2, 1, none, [a], true, 2)
    # Mechanisms, sites and allosteric mechanisms built twice are equal and hash alike;
    # each field takes part.
    examples, again = _testhelper_value_examples(), _testhelper_value_examples()
    @test all(examples .== again) && hash.(examples) == hash.(again)
    _, _, mech, site, allo = examples
    c = mech.catalytic
    @test mech != _testhelper_plain_mechanism(
        c.metabolites, c.enzyme_forms, c.substrates, c.products,
        [R.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off, true, e.kinetic_group)
         for e in c.steps])
    @test mech != _testhelper_plain_mechanism(
        c.metabolites, [:E, :EA, :EP, :EJ], c.substrates, c.products, c.steps)
    @test site != R.Mechanisms.RegulatorySite([:X], 4, [:OnlyI])
    @test site != R.Mechanisms.RegulatorySite([:X], 2, [:OnlyA])
    @test allo != R.Mechanism(
        mech.catalytic,
        R.Mechanisms.Allostery([:EqualAI, :NonequalAI, :EqualAI, :OnlyA], 2, [site]))
    @test allo != R.Mechanism(
        mech.catalytic,
        R.Mechanisms.Allostery(allo.allosteric.allosteric_states, 1, [site]))
    @test allo != R.Mechanism(
        mech.catalytic, R.Mechanisms.Allostery(allo.allosteric.allosteric_states, 2,
                                               R.Mechanisms.RegulatorySite[]))
end

@testset "a mechanism value's hash reads its type's name and its fields only" begin
    # Base hashes names, flags, numbers and vectors alike in every process, but a
    # struct's default hash mixes in its type's identity, which changes between
    # precompile sessions; so each type hashes its name, then its fields in order, all
    # of which `==` compares.
    h = UInt(17)
    for x in _testhelper_value_examples()
        @test hash(x, h) == foldl((acc, f) -> hash(getfield(x, f), acc),
                                  fieldnames(typeof(x));
                                  init = hash(nameof(typeof(x)), h))
    end
end

@testset "a regulatory-site ligand takes no reserved name and no `_`" begin
    # An allosteric regulator reaches the constructor only as a site ligand, so the
    # Mechanism constructor checks the ligands' names.
    @test_throws(ErrorException(
        "Mechanism: L is a reserved name; Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyA
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): L::OnlyI
        end)
    @test_throws(ErrorException(
        "Mechanism: the metabolite name A_R holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none"),
        @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ ES      :: EqualAI
                ES <--> EP    :: OnlyA
                EP ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): A_R::EqualAI
        end)
end

@testset "the composition constructor retags a catalytic graph" begin
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES      :: NonequalAI
            ES <--> EP    :: OnlyA
            EP ⇌ E + P      :: EqualAI
        end
        regulatory_site(2): R::OnlyI
    end
    # Canonical form numbers the groups E + P → EP, E + S → ES, ES → EP.
    _testhelper_retagged(tags) = EnzymeRates.Mechanism(
        am.catalytic, EnzymeRates.Mechanisms.Allostery(tags, am.allosteric.oligomeric_state,
                                                       am.allosteric.regulatory_sites))
    @test _testhelper_retagged([:EqualAI, :NonequalAI, :NonequalAI]) ==
          @enzyme_mechanism begin
              substrates: S
              products: P
              oligomeric_state: 2
              steps: begin
                  E + S ⇌ ES      :: NonequalAI
                  ES <--> EP    :: NonequalAI
                  EP ⇌ E + P      :: EqualAI
              end
              regulatory_site(2): R::OnlyI
          end
    @test_throws("must match the number of kinetic groups 3",
                 _testhelper_retagged([:EqualAI, :NonequalAI]))
end

@testset "a conformation's steps" begin
    M = ER.Mechanisms
    texts(am, mask) = Set(_testhelper_edge_text(am, k) for k in findall(mask))
    # The inactive conformation excludes the :OnlyA binding of S and
    # the :OnlyA chemistry, which cuts ES off; the binding of P still runs there.
    example = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES      :: OnlyA
            ES <--> EP    :: OnlyA
            E + P ⇌ EP      :: EqualAI
        end
    end
    runs = @inferred M.steps_running_in_allosteric_conformations(example)
    @test runs isa NamedTuple{(:active, :inactive), Tuple{BitVector, BitVector}}
    @test texts(example, runs.active) == Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])
    @test texts(example, runs.inactive) == Set(["E + P ⇌ EP"])
    # A plain mechanism runs every step in both conformations.
    plain = _testhelper_plain_of(example)
    @test @inferred(M.steps_running_in_allosteric_conformations(plain)) ==
          (; active = trues(3), inactive = trues(3))

    # An :EqualAI step between two cut-off forms does not run: the inactive conformation
    # excludes the :OnlyA binding of I, which cuts EI and EIP off, and with them the
    # binding of P that joins them.
    dead_end = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES          :: EqualAI
            ES <--> EP        :: EqualAI
            EP ⇌ E + P          :: EqualAI
            E + I ⇌ EI          :: OnlyA
            EI + P ⇌ EIP    :: EqualAI
        end
    end
    runs = M.steps_running_in_allosteric_conformations(dead_end)
    @test texts(dead_end, runs.active) ==
          Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP", "E + I ⇌ EI", "EI + P ⇌ EIP"])
    @test texts(dead_end, runs.inactive) == Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])

    # A ping-pong, F the modified enzyme E_res_plusA_minusP. The inactive
    # conformation excludes both :OnlyA chemistries, which cuts F, FP and FB off as one
    # island; the bindings of P and B inside it do not run either. The active
    # conformation excludes the :OnlyI binding of I, which cuts EI off.
    pingpong = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        modifiers: I
        steps: begin
            E + A ⇌ EA                                          :: EqualAI
            EA <--> EP_res_plusA_minusP                      :: OnlyA
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P    :: EqualAI
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP    :: EqualAI
            EB_res_plusA_minusP <--> EQ                      :: OnlyA
            EQ ⇌ E + Q                                          :: EqualAI
            E + I ⇌ EI                                          :: OnlyI
        end
    end
    chemistry = Set(["EA <--> EP_res_plusA_minusP", "EB_res_plusA_minusP <--> EQ"])
    island = Set(["E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP",
                  "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP"])
    runs = M.steps_running_in_allosteric_conformations(pingpong)
    @test texts(pingpong, runs.active) ==
          union(Set(["E + A ⇌ EA", "E + Q ⇌ EQ"]), chemistry, island)
    @test texts(pingpong, runs.inactive) == Set(["E + A ⇌ EA", "E + Q ⇌ EQ", "E + I ⇌ EI"])
end
