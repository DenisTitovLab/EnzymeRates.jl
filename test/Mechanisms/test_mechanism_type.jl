# ABOUTME: Tests of src/Mechanisms/mechanism_type.jl: the mechanism types, their
# ABOUTME: canonical form, checks, identity and text, segments, cycles and conformations.

# ── The mechanism types: fields, checks, identity and text ────────────────────────────

using Serialization

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

# ── Canonical form and storage ────────────────────────────────────────────────────────

using Statistics: median

const ER = EnzymeRates

"""Uni-uni with rapid-equilibrium bindings, written field by field in canonical form:
metabolites S and P, forms E, EP, ES; the bindings of P and S, then ES → EP."""
function _testhelper_graph_uniuni()
    M = ER.Mechanisms
    none = Int[]
    _testhelper_plain_mechanism(
        [:S, :P], [:E, :EP, :ES], [1], [2],
        [M.Step(1, 2, [2], none, true, 1), M.Step(1, 3, [1], none, true, 2),
         M.Step(3, 2, none, none, false, 3)])
end

"""Whether a value of type `T` hashes the same in every process: a symbol, an integer or
a flag, a vector of such values, or one of the types whose hash folds its name and its
fields (`_testhelper_hash_folds`), whose fields are such values in turn."""
function _testhelper_hash_is_content(T)
    T in (Symbol, Bool, Int) && return true
    T <: Vector && return _testhelper_hash_is_content(eltype(T))
    T in (EnzymeRates.Mechanisms.RegulatorySite, ER.Mechanisms.Step,
          ER.Mechanisms.CatalyticScheme, ER.Mechanisms.Allostery, ER.Mechanism) ||
        return false
    all(f -> _testhelper_hash_is_content(fieldtype(T, f)), fieldnames(T))
end

"""Whether `x` hashes its type's name, then its fields in declaration order."""
_testhelper_hash_folds(x, h) =
    hash(x, h) == foldl((acc, f) -> hash(getfield(x, f), acc), fieldnames(typeof(x));
                        init = hash(nameof(typeof(x)), h))

@testset "a graph hashes its content, the same in every process" begin
    # Julia hashes symbols, integers, flags and vectors of them alike in every process,
    # but a struct's default hash mixes in its type's identity, which changes between
    # precompile sessions. So each graph type hashes its name, then its fields, and its
    # fields hold only such values.
    m = _testhelper_graph_uniuni()
    c = m.catalytic
    site = EnzymeRates.Mechanisms.RegulatorySite([:R], 2, [:OnlyI])
    allostery(tags, n, sites) = ER.Mechanisms.Allostery(tags, n, sites)
    am = ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 2, [site]))
    h = UInt(17)
    for x in (c.steps[3], c, m.allosteric, m, am.allosteric, am, site)
        @test _testhelper_hash_folds(x, h)
    end
    @test all(_testhelper_hash_is_content,
              (ER.Mechanisms.Step, ER.Mechanisms.CatalyticScheme, ER.Mechanisms.Allostery,
               ER.Mechanism))
    # Each field takes part in `==` and `hash`.
    @test m == _testhelper_graph_uniuni() && hash(m) == hash(_testhelper_graph_uniuni())
    fields = fieldnames(ER.Mechanisms.CatalyticScheme)
    for f in fields
        values = Any[getfield(c, g) for g in fields]
        i = findfirst(==(f), fields)
        values[i] = f === :substrates ? [2] : f === :products ? [1] :
                    reverse(getfield(c, f))
        other = ER.Mechanisms._unchecked_scheme(values...)
        @test other != c && hash(other) != hash(c)
    end
    # The steps' groups take part too.
    regrouped = ER.Mechanisms._unchecked_scheme(c.metabolites, c.enzyme_forms, c.substrates,
        c.products,
        [ER.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, g)
         for (e, g) in zip(c.steps, [1, 3, 2])])
    @test regrouped != c && hash(regrouped) != hash(c)
    @test am == ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 2, [site]))
    @test am != ER.Mechanism(c, allostery([:EqualAI, :OnlyA, :OnlyA], 2, [site]))
    @test am != ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 4, [site]))
    @test am != ER.Mechanism(c, allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                                          EnzymeRates.Mechanisms.RegulatorySite[]))
end

@testset "a mechanism's fields are its graph" begin
    # No reaction, no cache and no fact derived from the graph: a mechanism is its
    # catalytic scheme and its allostery. The scheme holds the names of the metabolites
    # and the forms, the reactants by index and the steps with their groups; the
    # allostery the states, the oligomeric state and the sites.
    @test fieldnames(ER.Mechanism) == (:catalytic, :allosteric)
    @test fieldnames(ER.Mechanisms.CatalyticScheme) ==
          (:metabolites, :enzyme_forms, :substrates, :products, :steps)
    @test fieldnames(ER.Mechanisms.Step) ==
          (:from, :to, :taken_up, :given_off, :rapid_equilibrium, :kinetic_group)
    @test fieldnames(ER.Mechanisms.Allostery) ==
          (:allosteric_states, :oligomeric_state, :regulatory_sites)
end

@testset "construction is canonical" begin
    # Michaelis–Menten from its fields: (a) in canonical form, (b) its metabolites and
    # forms in another order, (c) every step reversed, (d) its group labels renumbered.
    # Every writing builds one mechanism, with one hash and one printed text.
    M = ER.Mechanisms
    none = Int[]
    steps = [M.Step(1, 2, [2], none, true, 1), M.Step(1, 3, [1], none, true, 2),
             M.Step(3, 2, none, none, false, 3)]
    canonical = _testhelper_plain_mechanism([:S, :P], [:E, :EP, :ES], [1], [2], steps)
    @test (canonical.catalytic.enzyme_forms, canonical.catalytic.steps) ==
          ([:E, :EP, :ES], steps)
    @test canonical == @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end
    permuted = _testhelper_plain_mechanism([:P, :S], [:ES, :E, :EP], [2], [1],
                                           [M.Step(2, 3, [1], none, true, 1),
                                            M.Step(2, 1, [2], none, true, 2),
                                            M.Step(1, 3, none, none, false, 3)])
    reversed = _testhelper_plain_mechanism([:S, :P], [:E, :EP, :ES], [1], [2],
                                           [M.Step(2, 1, none, [2], true, 1),
                                            M.Step(3, 1, none, [1], true, 2),
                                            M.Step(2, 3, none, none, false, 3)])
    relabelled_steps = [M.Step(1, 2, [2], none, true, 3), M.Step(1, 3, [1], none, true, 1),
                        M.Step(3, 2, none, none, false, 2)]
    relabelled = _testhelper_plain_mechanism([:S, :P], [:E, :EP, :ES], [1], [2],
                                             relabelled_steps)
    for m in (permuted, reversed, relabelled)
        @test m == canonical && hash(m) == hash(canonical)
        @test sprint(show, m) == sprint(show, canonical)
    end
    # A step's lists in any order: A and B taken up in one step, written as [B, A].
    _testhelper_bi_uni(up) = _testhelper_plain_mechanism(
        [:A, :B, :Q], [:E, :EAB], [1, 2], [3],
        [M.Step(1, 2, up, none, true, 1), M.Step(2, 1, none, [3], false, 2)])
    unsorted = _testhelper_bi_uni([2, 1])
    @test unsorted == _testhelper_bi_uni([1, 2])
    @test [e.taken_up for e in unsorted.catalytic.steps] == [[1, 2], [3]]
    # An allosteric mechanism from its written parts takes each tag by written group label:
    # label 1 (S's binding) is canonical group 2, label 2 (the chemistry) group 3 and
    # label 3 (P's binding) group 1.
    site = M.RegulatorySite([:R], 2, [:OnlyI])
    am = ER.Mechanism(
        ER.Mechanisms.CatalyticScheme([:S, :P], [:E, :EP, :ES], [1], [2], relabelled_steps),
        ER.Mechanisms.Allostery([:EqualAI, :OnlyA, :NonequalAI], 2, [site]))
    expected = ER.Mechanism(
        canonical.catalytic,
        ER.Mechanisms.Allostery([:NonequalAI, :EqualAI, :OnlyA], 2, [site]))
    @test am == expected && hash(am) == hash(expected)
    @test sprint(show, am) == sprint(show, expected)
end

@testset "an allosteric graph sorts its sites and checks its tags and ligands" begin
    m = _testhelper_graph_uniuni()
    tags = [:EqualAI, :EqualAI, :OnlyA]
    sites = [EnzymeRates.Mechanisms.RegulatorySite([:X], 2, [:OnlyI]),
             EnzymeRates.Mechanisms.RegulatorySite([:B, :A], 2, [:OnlyA, :EqualAI])]
    am = ER.Mechanism(m.catalytic, ER.Mechanisms.Allostery(tags, 2, sites))
    @test [s.ligands for s in am.allosteric.regulatory_sites] == [[:A, :B], [:X]]
    @test am.allosteric.regulatory_sites[1].allosteric_states == [:EqualAI, :OnlyA]
    tags[1] = :OnlyA
    @test am.allosteric.allosteric_states == [:EqualAI, :EqualAI, :OnlyA]
    @test [am.allosteric.allosteric_states[g] for g in 1:3] == [:EqualAI, :EqualAI, :OnlyA]
    @test_throws(ErrorException("Mechanism: allosteric_states length 2 must " *
                                "match the number of kinetic groups 3"),
                 ER.Mechanism(m.catalytic,
                              ER.Mechanisms.Allostery([:EqualAI, :EqualAI], 2, sites)))
    @test_throws(ErrorException("Mechanism: oligomeric_state must be ≥ " *
                                "1, got 0"),
                 ER.Mechanism(m.catalytic, ER.Mechanisms.Allostery(
                     [:EqualAI, :EqualAI, :OnlyA], 0, sites)))
    # Any group may be :OnlyI, and the guard judges the mechanism: the active
    # conformation cannot switch off the rapid-equilibrium S binding (step 2), whose two
    # sides the chemistry and the P binding still join. Neither fix the message offers
    # passes: each leaves the active conformation no turnover. Untagging the binding
    # passes, but is not one of them.
    onlyi = ER.Mechanism(m.catalytic,
                         ER.Mechanisms.Allostery([:EqualAI, :OnlyI, :OnlyA], 2, sites))
    v = ER.Derivation._guard(onlyi)
    @test (v.fault, v.conformation, v.edge) == (ER.Derivation._FAST_STEP_PINNED, 1, 2)
    @test ER.Derivation._guard_message(onlyi) == _testhelper_pinned_cause(onlyi, 2, 1)
    @test all(t -> ER.Derivation._retag_passes(onlyi, [2], t), (:EqualAI, :NonequalAI))
    @test_throws(ErrorException(
        "Mechanism: allosteric regulator `X` is listed at two regulatory sites; " *
        "an allosteric regulator sits at exactly one regulatory site."),
                 ER.Mechanism(
                     m.catalytic, ER.Mechanisms.Allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                     [sites; EnzymeRates.Mechanisms.RegulatorySite([:X], 4, [:OnlyA])])))
    @test_throws(ErrorException(
        "Mechanism: L is a reserved name; Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        ER.Mechanism(m.catalytic, ER.Mechanisms.Allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                               [EnzymeRates.Mechanisms.RegulatorySite([:L], 2, [:OnlyI])])))
    # A ligand named like a form would read as that form in the mechanism's text.
    @test_throws(ErrorException(
        "Mechanism: ligand ES is also the name of an enzyme form; a text could not tell " *
        "them apart"),
        ER.Mechanism(m.catalytic,
                     ER.Mechanisms.Allostery([:EqualAI, :EqualAI, :OnlyA], 2,
                                             [ER.Mechanisms.RegulatorySite([:ES], 2,
                                                                           [:OnlyI])])))
end

"""`raw`, a spec's raw value (`raw_mechanism`), written another way, the same mechanism:
with `order`, its forms and its steps in reverse order; with `direction`, each step from
its other end; with `labels`, its group labels numbered from the last. Its declared names
come in reverse order, the steps' lists following them. Its tags follow their labels, and
its sites and each site's ligands come in reverse order."""
function _testhelper_rewrite(raw::NamedTuple; order = true, direction = true,
                             labels = true)
    M = EnzymeRates.Mechanisms
    forms, steps, tags = raw.forms, raw.steps, raw.allosteric_states
    G, n = maximum(e -> e.kinetic_group, steps; init = 0), length(forms)
    declared = [raw.substrates; raw.products; raw.modifiers]
    redeclared = [reverse(raw.substrates); reverse(raw.products); reverse(raw.modifiers)]
    at = [findfirst(==(x), redeclared) for x in declared]
    steps = [M.Step(e.from, e.to, at[e.taken_up], at[e.given_off], e.rapid_equilibrium,
                    e.kinetic_group) for e in steps]
    direction && (steps = [M.Step(e.to, e.from, e.given_off, e.taken_up,
                                  e.rapid_equilibrium, e.kinetic_group) for e in steps])
    labels && ((steps, tags) = ([M.Step(e.from, e.to, e.taken_up, e.given_off,
                                        e.rapid_equilibrium, G + 1 - e.kinetic_group)
                                 for e in steps], reverse(tags)))
    order && ((forms, steps) = (reverse(forms),
                                reverse([M.Step(n + 1 - e.from, n + 1 - e.to, e.taken_up,
                                                e.given_off, e.rapid_equilibrium,
                                                e.kinetic_group)
                                         for e in steps])))
    (; substrates = reverse(raw.substrates), products = reverse(raw.products),
     modifiers = reverse(raw.modifiers), forms, steps, allosteric_states = tags,
     oligomeric_state = raw.oligomeric_state,
     regulatory_sites = [(; ligands = reverse(s.ligands), s.multiplicity,
                          allosteric_states = reverse(s.allosteric_states))
                         for s in reverse(raw.regulatory_sites)])
end

@testset "every writing of a mechanism gives one canonical graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        m = _testhelper_build(spec.raw)
        for (order, direction, labels) in ((true, false, false), (false, true, false),
                                           (false, false, true), (true, true, true))
            other = _testhelper_build(
                _testhelper_rewrite(spec.raw; order, direction, labels))
            @test other == m && hash(other) == hash(m)
            @test compile_rate_equation(other) === compile_rate_equation(m)
        end
    end
end

"""Whether the constructor stores each step of the mechanism block `block` turned from the
way it is written: the stored step joining its two forms with its metabolites runs from
its written second form."""
function _testhelper_flipped(block)
    raw = EnzymeRates.DSL.raw_mechanism(block)
    c = _testhelper_scheme(raw)
    declared = [raw.substrates; raw.products; raw.modifiers]
    side(f, list) = (f, sort(list))
    stored = [(side(c.enzyme_forms[s.from], c.metabolites[s.taken_up]),
               side(c.enzyme_forms[s.to], c.metabolites[s.given_off])) for s in c.steps]
    map(raw.steps) do e
        a = side(raw.forms[e.from], declared[e.taken_up])
        b = side(raw.forms[e.to], declared[e.given_off])
        (a, b) in stored && return false
        (b, a) in stored && return true
        error("no stored step joins $a and $b")
    end
end

@testset "canonical form stores each step the way the graph gives it, a group as one" begin
    # A step with one metabolite takes it up, so the fused release FP → F + P is stored
    # as F + P → FP. The isomerizations lie on the one cycle that turns the reaction
    # over, E → ES → FS → FP → F → E, and run its way.
    @test _testhelper_flipped(:(begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            FS <--> ES
            FS <--> FP
            FP <--> F + P
            F <--> E
        end
    end)) == [false, true, false, true, false]
    @test _testhelper_flipped(:(begin
        substrates: S
        products: P
        steps: begin
            ES <--> E + S
            ES <--> FS
            FP <--> FS
            F + P <--> FP
            E <--> F
        end
    end)) == [true, false, true, false, true]
    # F → E runs the way the reaction turns over and leads its kinetic group, nearer free
    # E than its mirror at the I-bound forms, which keeps the direction it is written in
    # relative to F → E. Both written the other way round, the group is stored the same
    # way; written against each other, they write another mechanism. The release of P is
    # stored as the binding it reverses.
    _testhelper_mirrored(free, bound) = :(begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> F + P
            ($free, $bound)
            (E + I ⇌ EI, F + I ⇌ FI)
        end
    end)
    c = _testhelper_mirrored(:(F <--> E), :(FI <--> EI))
    d = _testhelper_mirrored(:(E <--> F), :(EI <--> FI))
    against = _testhelper_mirrored(:(F <--> E), :(EI <--> FI))
    built(block) = _testhelper_build(EnzymeRates.DSL.raw_mechanism(block))
    @test _testhelper_flipped(c) == [false, false, true, false, false, false, false]
    @test _testhelper_flipped(d) == [false, false, true, true, true, false, false]
    @test built(c) == built(d)
    @test _testhelper_flipped(against) == [false, false, true, false, false, false, false]
    @test built(against) != built(c)
end

@testset "a graph holds any number of metabolites, and a step any number of one" begin
    M = ER.Mechanisms
    names = [[Symbol("S", i) for i in 1:32]; :P]
    # Each metabolite bound to free E into a form of its own, in a group of its own.
    _testhelper_bindings(n) = _testhelper_plain_mechanism(
        names[[1:n; 33]], [:E; [Symbol("E", x) for x in names[[1:n; 33]]]],
        collect(1:n), [n + 1],
        [M.Step(1, i + 1, [i], Int[], true, i) for i in 1:(n + 1)])
    @test length(_testhelper_bindings(31).catalytic.metabolites) == 32
    @test length(_testhelper_bindings(32).catalytic.metabolites) == 33
    for n in (15, 16)
        c = M.CatalyticScheme([:A, :P], [:E, :F], [1], [2],
                              [M.Step(1, 2, fill(1, n), [2], false, 1)])
        @test count(==(1), only(c.steps).taken_up) == n
    end
end

@testset "a step turned round swaps its ends and its lists" begin
    e = ER.Mechanisms.Step(2, 5, [1, 3], [4], true, 7)
    @test ER.Mechanisms._reversed(e) == ER.Mechanisms.Step(5, 2, [4], [1, 3], true, 7)
    @test ER.Mechanisms._reversed(ER.Mechanisms._reversed(e)) == e
end

"""The raw value of `block`, a mechanism macro's block, read by the DSL walker."""
_testhelper_walk(block) = EnzymeRates.DSL.raw_mechanism(block)

@testset "the constructors build one value from every writing of a mechanism" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        # The mechanism holds the canonical scheme, its groups numbered by first step.
        c = _testhelper_scheme(spec.raw)
        @test (x.catalytic.metabolites, x.catalytic.enzyme_forms, x.catalytic.substrates,
               x.catalytic.products) ==
              (c.metabolites, c.enzyme_forms, c.substrates, c.products)
        unlabelled(steps) = [(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium)
                             for e in steps]
        @test unlabelled(x.catalytic.steps) == unlabelled(c.steps)
        # The same partition of the steps into kinetic groups, numbered by first step.
        groups(s) = Set(Set(ks) for ks in _testhelper_group_members(s))
        @test groups(x.catalytic) == groups(c)
        @test unique(e.kinetic_group for e in x.catalytic.steps) == 1:length(groups(c))
        other = _testhelper_build(_testhelper_rewrite(spec.raw))
        @test other == x && hash(other) == hash(x)
    end
    # Each tag goes to the canonical group its label became.
    am = _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES :: EqualAI
            ES <--> EP :: OnlyA
            EP ⇌ E + P :: NonequalAI
        end
    end)))
    @test am.allosteric.allosteric_states == [:NonequalAI, :EqualAI, :OnlyA]
end

"""Whether `m`, a catalytic scheme or a mechanism, stores each of its vectors at its
length: `m` takes the memory of the value rebuilt from a copy of each vector, its steps'
lists and its sites included, which `copy` stores at its length."""
_testhelper_stored_at_length(m::ER.Mechanisms.CatalyticScheme) =
    Base.summarysize(m) == Base.summarysize(
        ER.Mechanisms._unchecked_scheme(copy(m.metabolites), copy(m.enzyme_forms),
            copy(m.substrates), copy(m.products),
            [ER.Mechanisms.Step(e.from, e.to, copy(e.taken_up), copy(e.given_off),
                                e.rapid_equilibrium, e.kinetic_group) for e in m.steps]))
_testhelper_stored_at_length(am::ER.Mechanism) =
    _testhelper_stored_at_length(am.catalytic) && Base.summarysize(am) ==
    Base.summarysize(ER.Mechanism(am.catalytic,
        ER.Mechanisms.Allostery(copy(am.allosteric.allosteric_states),
                                am.allosteric.oligomeric_state,
                                [ER.Mechanisms.RegulatorySite(copy(s.ligands),
                                                              s.multiplicity,
                                                              copy(s.allosteric_states))
                                 for s in am.allosteric.regulatory_sites])))

@testset "a mechanism stores each vector at its length" begin
    # A search keeps every mechanism it emits, so memory a vector holds past its length
    # is paid once per mechanism: the ter-ter init set holds 250,855 of them.
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_stored_at_length(_testhelper_build(spec.raw))
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
        competitive_inhibitors: I
    end
    init = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
    children = _testhelper_children(rxn, init[1:10])
    @test any(c -> _testhelper_is_allosteric(c), children)
    @test all(_testhelper_stored_at_length, init)
    @test all(_testhelper_stored_at_length, children)
end

@testset "the constructors reject what no graph represents, each with its message" begin
    # Uni-uni over S and P, E + S ⇌ ES, ES <--> EP, EP ⇌ E + P, forms E, ES and EP (1 to
    # 3), with the forms `extra` (4, 5, …) and `steps` added.
    M = ER.Mechanisms
    none = Int[]
    uniuni = [M.Step(1, 2, [1], none, true, 1), M.Step(2, 3, none, none, false, 2),
              M.Step(3, 1, none, [2], true, 3)]
    _testhelper_uniuni(extra = Symbol[], steps...; metabolites = [:S, :P],
                       substrates = [1], products = [2], forms = [:E, :ES, :EP],
                       group = [1, 2, 3]) =
        _testhelper_plain_mechanism(metabolites, [forms; extra], substrates, products,
                     [[M.Step(e.from, e.to, e.taken_up, e.given_off, e.rapid_equilibrium, g)
                       for (e, g) in zip(uniuni, group)]; steps...])
    # Uni-uni with the modifier `x`, metabolite 3, bound to free E by a step of its own.
    _testhelper_binding(x; kw...) =
        _testhelper_uniuni([Symbol(:E, x)], M.Step(1, 4, [3], none, true, 4);
                           metabolites = [:S, :P, x], kw...)
    _testhelper_allosteric(tags, n, sites) =
        ER.Mechanism(ER.Mechanisms.CatalyticScheme([:S, :P], [:E, :ES, :EP], [1], [2],
                                                   uniuni),
                     ER.Mechanisms.Allostery(tags, n, sites))
    tags = [:EqualAI, :OnlyA, :EqualAI]
    site = M.RegulatorySite
    @test _testhelper_uniuni() == _testhelper_graph_uniuni()
    # Regulatory sites.
    @test_throws(ErrorException("RegulatorySite: needs at least one ligand"),
                 _testhelper_allosteric(tags, 2, [site(Symbol[], 2, Symbol[])]))
    @test_throws(ErrorException("RegulatorySite: length(ligands)=1 must equal " *
                                "length(allosteric_states)=0"),
                 _testhelper_allosteric(tags, 2, [site([:R], 2, Symbol[])]))
    @test_throws(ErrorException("RegulatorySite: multiplicity must be ≥ 1, got 0"),
                 _testhelper_allosteric(tags, 2, [site([:R], 0, [:OnlyA])]))
    @test_throws(ErrorException("RegulatorySite: allosteric state Bogus must be one of " *
                                "(:OnlyA, :OnlyI, :EqualAI, :NonequalAI)"),
                 _testhelper_allosteric(tags, 2, [site([:R], 2, [:Bogus])]))
    @test_throws(ErrorException("RegulatorySite: allosteric regulator `R` is listed " *
                                "twice; an allosteric regulator sits at exactly one " *
                                "regulatory site."),
                 _testhelper_allosteric(tags, 2, [site([:R, :X, :R], 2,
                                                       [:OnlyA, :OnlyI, :OnlyI])]))
    @test_throws(ErrorException(
        "Mechanism: allosteric regulator `R` is listed at two regulatory sites; " *
        "an allosteric regulator sits at exactly one regulatory site."),
        _testhelper_allosteric(tags, 2, [site([:R], 2, [:OnlyA]), site([:R], 4, [:OnlyI])]))
    # A ligand takes no name a metabolite cannot take, and no form's name.
    @test_throws(ErrorException("Mechanism: \"1R\" cannot name a metabolite; a " *
                                "metabolite name is an identifier that is no keyword"),
                 _testhelper_allosteric(tags, 2, [site([Symbol("1R")], 2, [:OnlyA])]))
    @test_throws(ErrorException("Mechanism: E cannot name a metabolite; it " *
                                "labels the free enzyme form"),
                 _testhelper_allosteric(tags, 2, [site([:E], 2, [:OnlyA])]))
    @test_throws(ErrorException("Mechanism: ligand ES is also the name of an enzyme " *
                                "form; a text could not tell them apart"),
                 _testhelper_allosteric(tags, 2, [site([:ES], 2, [:OnlyA])]))
    # Tags, by written group label: label 2 is canonical group 3.
    @test_throws(ErrorException("Mechanism: allosteric_states length 2 must " *
                                "match the number of kinetic groups 3"),
                 _testhelper_allosteric(tags[1:2], 2, M.RegulatorySite[]))
    @test_throws(ErrorException("Mechanism: oligomeric_state must be ≥ " *
                                "1, got 0"),
                 _testhelper_allosteric(tags, 0, M.RegulatorySite[]))
    # A Bool is no multiplicity: each constructor takes an Int and converts no Bool.
    @test_throws MethodError site([:R], true, [:OnlyA])
    @test_throws MethodError _testhelper_allosteric(tags, true, M.RegulatorySite[])
    @test_throws MethodError ER.Mechanism(_testhelper_uniuni().catalytic,
                                          ER.Mechanisms.Allostery(tags, true,
                                                                  M.RegulatorySite[]))
    # Any group may be :OnlyI: the constructor builds the chemistry :OnlyI, and the guard
    # rejects it, since the active conformation then turns nothing over.
    onlyi = _testhelper_allosteric([:EqualAI, :OnlyI, :EqualAI], 2, M.RegulatorySite[])
    v = ER.Derivation._guard(onlyi)
    @test (v.fault, v.conformation) == (ER.Derivation._NO_ACTIVE_TURNOVER, 1)
    @test startswith(ER.Derivation._guard_message(onlyi),
                     "Mechanism: without its :OnlyI steps ES <--> EP the " *
                     "active conformation turns no substrate into product")
    # Declarations and group labels.
    @test_throws(ErrorException("Mechanism: substrates must not be empty"),
                 _testhelper_uniuni(substrates = Int[]))
    @test_throws(ErrorException("Mechanism: products must not be empty"),
                 _testhelper_uniuni(products = Int[]))
    @test_throws(ErrorException("Mechanism: duplicate substrate names"),
                 _testhelper_uniuni(metabolites = [:S, :P, :S], substrates = [1, 3]))
    @test_throws(ErrorException("Mechanism: duplicate product names"),
                 _testhelper_uniuni(metabolites = [:S, :P, :P], products = [2, 3]))
    @test_throws(ErrorException("Mechanism: S named as both a substrate and a product; " *
                                "concentrations and constants are keyed by name"),
                 _testhelper_uniuni(products = [2, 1]))
    @test_throws(ErrorException("Mechanism: duplicate modifier names"),
                 _testhelper_uniuni(metabolites = [:S, :P, :I, :I]))
    # A reactant that binds as a modifier is taken up as itself, so it is not declared
    # again as a modifier.
    @test_throws(ErrorException("Mechanism: S is declared both as a reactant and as a " *
                                "modifier; a step that binds a reactant as a modifier " *
                                "simply takes it up, as `ES + S ⇌ ESS` does"),
                 _testhelper_uniuni(metabolites = [:S, :P, :S]))
    @test_throws(ErrorException("Mechanism: each step takes one kinetic group label " *
                                "from 1 to the number of groups"),
                 _testhelper_uniuni(group = [0, 1, 2]))
    @test_throws(ErrorException("Mechanism: kinetic group 2 is empty; a kinetic group " *
                                "holds at least one step"),
                 _testhelper_uniuni(group = [1, 3, 3]))
    # A step's ends index the forms, and its lists the metabolites.
    for (from, to) in ((2, 4), (0, 1))
        @test_throws(ErrorException("Mechanism: step 4 joins forms $from and $to, but " *
                                    "there are 3 forms; a step's ends index the enzyme " *
                                    "forms"),
                     _testhelper_uniuni(Symbol[], M.Step(from, to, none, none, false, 4)))
    end
    @test_throws(ErrorException("Mechanism: step 4 takes up [3] and gives off Int64[], " *
                                "but there are 2 metabolites; a step's lists index the " *
                                "metabolites"),
                 _testhelper_uniuni([:EI], M.Step(1, 4, [3], none, true, 4)))
    # The substrates and the products index the metabolites too.
    @test_throws(ErrorException("Mechanism: substrates [3], but there are 2 metabolites; " *
                                "the substrates and products index the metabolites"),
                 _testhelper_uniuni(substrates = [3]))
    @test_throws(ErrorException("Mechanism: products [0], but there are 2 metabolites; " *
                                "the substrates and products index the metabolites"),
                 _testhelper_uniuni(products = [0]))
    # Names a metabolite cannot take: no identifier, a keyword, E, a reserved name, a
    # name holding `_`.
    _testhelper_unwritable(x) = ErrorException(
        "Mechanism: \"$x\" cannot name a metabolite; a metabolite name is an identifier " *
        "that is no keyword")
    @test_throws _testhelper_unwritable("1I") _testhelper_binding(Symbol("1I"))
    @test_throws(_testhelper_unwritable("end"),
                 _testhelper_uniuni(metabolites = [:S, :P, :end], products = [2, 3]))
    @test_throws(ErrorException("Mechanism: E cannot name a metabolite; it labels the " *
                                "free enzyme form"),
                 _testhelper_uniuni(metabolites = [:S, :P, :E], substrates = [1, 3]))
    @test_throws(ErrorException("Mechanism: Keq is a reserved name; Keq, L, N, " *
                                "params, concs and v are kept for a rate equation's " *
                                "constants and local names"),
                 _testhelper_uniuni(metabolites = [:S, :P, :Keq], substrates = [1, 3]))
    @test_throws(ErrorException("Mechanism: the metabolite name S_1 holds `_`; constant " *
                                "names join metabolite names with `_`, so a metabolite " *
                                "name holds none"),
                 _testhelper_binding(:S_1))
    # A name Julia reads as another (a micro sign as μ, a decomposed é as one character)
    # prints as text that builds another mechanism, so no DSL text can write it. The
    # strings are written with escapes: Julia would read the characters themselves as
    # the other names.
    _testhelper_read_as(who, x, read) = ErrorException(
        "$who: Julia reads the name \"$x\" as \"$read\", so no DSL text can write it; " *
        "use \"$read\"")
    @test_throws(_testhelper_read_as("Mechanism", "\u00b5M", "\u03bcM"),
                 _testhelper_binding(Symbol("\u00b5M")))
    @test_throws(_testhelper_read_as("Mechanism", "e\u0301", "\u00e9"),
                 _testhelper_uniuni(metabolites = [:S, :P, Symbol("e\u0301")],
                                    substrates = [1, 3]))
    @test_throws(_testhelper_read_as("Mechanism", "\u025b", "\u03b5"),
                 _testhelper_allosteric(tags, 2, [site([Symbol("\u025b")], 2, [:OnlyA])]))
    # Form names: no identifier, read as another, given twice, a metabolite's name.
    @test_throws(ErrorException("Mechanism: \"1E\" cannot name an enzyme form; an " *
                                "enzyme form name is an identifier that is no keyword"),
                 _testhelper_uniuni([Symbol("1E")], M.Step(1, 4, none, none, false, 4)))
    @test_throws(_testhelper_read_as("Mechanism", "E\u00b5", "E\u03bc"),
                 _testhelper_uniuni([Symbol("E\u00b5")],
                                    M.Step(1, 4, none, none, false, 4)))
    @test_throws(ErrorException("Mechanism: two enzyme forms are named ES"),
                 _testhelper_uniuni([:ES], M.Step(1, 4, none, none, false, 4)))
    @test_throws(ErrorException("Mechanism: S names both an enzyme form and a " *
                                "metabolite; a text could not tell them apart"),
                 _testhelper_uniuni([:S], M.Step(1, 4, none, none, false, 4)))
    # Steps: a self-loop, a metabolite taken up and given off.
    @test_throws(ErrorException("Step: both ends are the form ES; a step joins two forms"),
                 _testhelper_uniuni(Symbol[], M.Step(2, 2, none, none, false, 4)))
    @test_throws(ErrorException("Step: both ends are the form EA; a step joins two forms"),
                 _testhelper_plain_mechanism(
                     [:A, :P], [:E, :EA], [1], [2],
                     [M.Step(1, 2, [1], none, true, 1), M.Step(2, 1, none, [2], false, 2),
                      M.Step(2, 2, none, none, false, 3)]))
    @test_throws(ErrorException("Step EA → EQ: B both consumed and released"),
                 _testhelper_build(_testhelper_walk(:(begin
                     substrates: A, B
                     products: Q
                     steps: begin
                         E + A ⇌ EA
                         EA + B <--> EQ + B
                         EQ ⇌ E + Q
                     end
                 end))))
    # One copy of B taken up and two given off: the sides share B though their counts
    # differ. The step gives off more substrate than it takes up, so it is stored the
    # other way round.
    @test_throws(ErrorException("Step EQ → EA: B both consumed and released"),
                 _testhelper_build(_testhelper_walk(:(begin
                     substrates: A, B
                     products: Q
                     steps: begin
                         E + A ⇌ EA
                         EA + B <--> EQ + B + B
                         EQ ⇌ E + Q
                     end
                 end))))
    # A kinetic group of steps of different kinds or flags.
    for (steps, text) in (
            (:(begin
                 E + A <--> EA
                 (EA + B <--> EAB, EA + B <--> EQ + P)
                 EQ <--> E + Q
             end), "EA + B <--> EAB and EA + B <--> EQ + P"),
            (:(begin
                 (E + A ⇌ EA, EQ + A <--> EAQ)
                 EA <--> EQ + P
                 EQ <--> E + Q
             end), "E + A ⇌ EA and EQ + A <--> EAQ"))
        @test_throws(ErrorException(
            "Mechanism: a kinetic group holds $text, which differ in kind or in RE/SS; " *
            "the steps of a kinetic group share their constants, so they must be the " *
            "same kind of step with the same flag"),
            _testhelper_build(_testhelper_walk(:(begin
                substrates: A, B
                products: P, Q
                steps: $steps
            end))))
    end
    # One reaction twice: in two groups, with different flags, stored once each way, and in
    # one group.
    _testhelper_twice(steps) = _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: $steps
    end)))
    @test_throws(ErrorException(
        "Mechanism: kinetic groups 2 and 3 both hold the reaction E + S <--> ES (SS in " *
        "group 2, RE in group 3); a reaction belongs to one kinetic group, whose " *
        "constants are named after it"),
        _testhelper_twice(:(begin
            E + S <--> ES
            E + S ⇌ ES
            ES ⇌ EP
            EP ⇌ E + P
        end)))
    @test_throws(ErrorException(
        "Mechanism: kinetic groups 3 and 4 both hold the reaction ES <--> EP; a reaction " *
        "belongs to one kinetic group, whose constants are named after it"),
        _testhelper_twice(:(begin
            E + S <--> ES
            ES <--> EP
            EP <--> ES
            EP <--> E + P
        end)))
    # The tied pair and the lone step both store F → E, which runs the way the reaction
    # turns over, in groups 5 and 6.
    @test_throws(ErrorException(
        "Mechanism: kinetic groups 5 and 6 both hold the reaction F <--> E; a reaction " *
        "belongs to one kinetic group, whose constants are named after it"),
        _testhelper_build(_testhelper_walk(:(begin
            substrates: A
            products: P
            modifiers: I
            steps: begin
                E + A <--> EA
                EA <--> EP
                EP <--> F + P
                F <--> E
                (E <--> F, EI <--> FI)
                E + I ⇌ EI
                F + I ⇌ FI
                EI + P <--> EIP
                FI + A <--> FAI
            end
        end))))
    @test_throws(ErrorException(
        "Mechanism: kinetic group 2 holds the reaction E + S <--> ES twice; a mechanism " *
        "holds each reaction once"),
        _testhelper_twice(:(begin
            (E + S <--> ES, ES <--> E + S)
            ES <--> EP
            EP <--> E + P
        end)))
    # No form named E: no step at all, a free form named F, forms whose names start with
    # E but are not E, or a ping-pong whose free form is F and whose covalent forms'
    # names start with E.
    no_E = ArgumentError("Mechanism: no enzyme form is named `E`; the free enzyme is the " *
                         "form named `E`")
    @test_throws no_E _testhelper_plain_mechanism([:S, :P], Symbol[], [1], [2], M.Step[])
    @test_throws no_E _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            F + S ⇌ FS
            FS <--> FP
            FP ⇌ F + P
        end
    end)))
    @test_throws no_E _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            F + S ⇌ ES
            ES <--> EP
            EP ⇌ F + P
        end
    end)))
    @test_throws no_E _testhelper_build(_testhelper_walk(:(begin
        substrates: A, B
        products: P, Q
        steps: begin
            F + A ⇌ FA
            FA <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> FQ
            FQ ⇌ F + Q
        end
    end)))
    # Forms that no steps join to free E: an isomerization between unreached forms.
    @test_throws(ErrorException("Mechanism: the steps do not form one connected graph; " *
                                "no steps lead from free E to EstarP, EstarS"),
                 _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            EstarS <--> EstarP
        end
    end))))
    # Declared names no step exchanges: a substrate, a product, a modifier.
    @test_throws(ErrorException("Mechanism: substrate T is declared, but no step uses it"),
                 _testhelper_uniuni(metabolites = [:S, :P, :T], substrates = [1, 3]))
    @test_throws(ErrorException("Mechanism: product Q is declared, but no step uses it"),
                 _testhelper_uniuni(metabolites = [:S, :P, :Q], products = [2, 3]))
    @test_throws(ErrorException("Mechanism: modifier I is declared, but no step binds it"),
                 _testhelper_uniuni(metabolites = [:S, :P, :I]))
end

@testset "a graph written back as its forms and steps builds the same graph" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        x = _testhelper_build(spec.raw)
        raw = _testhelper_written(x)
        @test _testhelper_build(raw) == x
        # Its fields are stored as given: in canonical order, each step in its stored
        # direction, labelled by group.
        c = _testhelper_scheme(raw)
        @test c == x.catalytic
    end
end

# ── The canonical form pinned as text, and the sub-graph reads against oracles ──

"""The regression census's six reactions: uni-uni, bi-bi, ping-pong, two sites, ter-bi and
the reaction whose names hold one another as prefixes."""
function _testhelper_census_reactions()
    [
        @enzyme_reaction(begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: I
            allosteric_regulators: R
            oligomeric_state: 2
        end),
        @enzyme_reaction(begin
            substrates: A[C], B[N]
            products: P[C], Q[N]
            competitive_inhibitors: A, I
        end),
        @enzyme_reaction(begin
            substrates: A[CX], B[N]
            products: P[C], Q[NX]
            allosteric_regulators: R
            oligomeric_state: 2
            shared_catalytic_site: (A, P)
        end),
        @enzyme_reaction(begin
            substrates: A[C6H12O6]
            products: B[C6H12O6]
            allosteric_regulators: X::Activator, Y::Inhibitor
            oligomeric_state: 2
        end),
        @enzyme_reaction(begin
            substrates: A[C], B[N], D[X]
            products: P[CN], Q[X]
        end),
        @enzyme_reaction(begin
            substrates: G[CX], GT[N]
            products: Gl[C], U[NX]
            competitive_inhibitors: GTI
            shared_catalytic_site: (G, Gl)
        end)]
end

"""The mechanisms the regression census enumerates first: the init populations of its six
reactions, then the children of every init mechanism of its bi-bi reaction."""
function _testhelper_census_population()
    reactions = _testhelper_census_reactions()
    population = EnzymeRates.Mechanism[]
    foreach(r -> append!(population, _testhelper_mechanisms(
        EnzymeRates.Enumeration.init_mechanisms(r))), reactions)
    bibi = reactions[2]
    append!(population,
            _testhelper_children(bibi, _testhelper_mechanisms(
                EnzymeRates.Enumeration.init_mechanisms(bibi))))
end

"""The census's first allosteric mechanisms: the allosteric children of the init
mechanisms of its three reactions that declare an allosteric regulator."""
function _testhelper_census_allosteric()
    children = [_testhelper_children(
        r, _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(r)))
                for r in _testhelper_census_reactions()[[1, 3, 4]]]
    filter(x -> _testhelper_is_allosteric(x), reduce(vcat, children))
end

"""Graph `x`'s canonical form as text lines: its names, its forms in order, each edge in
its stored direction (`_testhelper_edge_text`) with its group (and the group's tag), then
any oligomeric state and sites."""
function _testhelper_canonical_text(x)
    m, allo = x.catalytic, _testhelper_is_allosteric(x)
    g = _testhelper_graph(x)
    lines = ["substrates: $(join(g.substrates, ", ")) | products: " *
             "$(join(g.products, ", "))" *
             (isempty(g.modifiers) ? "" : " | modifiers: $(join(g.modifiers, ", "))"),
             "forms: " * join(m.enzyme_forms, ", ")]
    for k in eachindex(m.steps)
        g = m.steps[k].kinetic_group
        push!(lines, "$k: $(_testhelper_edge_text(m, k)) | group $g" *
                     (allo ? " :: $(x.allosteric.allosteric_states[g])" : ""))
    end
    allo && push!(lines, "oligomeric_state: $(x.allosteric.oligomeric_state)")
    for s in x.allosteric.regulatory_sites
        push!(lines, "regulatory_site($(s.multiplicity)): " *
                     join(["$l::$t" for (l, t) in zip(s.ligands, s.allosteric_states)], ", "))
    end
    lines
end

"""
The inactive edges of the allosteric mechanism whose test reader (`_testhelper_graph`) is
`g`, read off `g`'s steps alone, one flag per step: no step of an `:OnlyA` group, and every
other step whose two forms the steps outside `:OnlyA` groups reach from free E.
"""
function _testhelper_inactive_reach(g)
    outside = [g.tags[s.group] !== :OnlyA for s in g.steps]
    reach = falses(length(g.forms))
    reach[1] = true
    grew = true
    while grew
        grew = false
        for (j, s) in enumerate(g.steps)
            outside[j] && reach[s.from] != reach[s.to] || continue
            reach[s.from] = reach[s.to] = true
            grew = true
        end
    end
    BitVector([outside[j] && reach[s.from] && reach[s.to] for (j, s) in enumerate(g.steps)])
end

"""
Whether the sub-graph reads of `x` (`rapid_equilibrium_segments`,
`steps_running_in_allosteric_conformations`) give what the oracles read off its test
reader's steps give (`_testhelper_re_segments`, `_testhelper_inactive_reach`): the
catalytic graph's segments, form by form and breadth-first, with each form's positive
offsets, and, for an allosteric mechanism, the inactive conformation's edges and their
segments.
"""
function _testhelper_subgraphs_agree(x)
    m, g = x.catalytic, _testhelper_graph(x)
    # Step `j` of the reader is edge `flat[j]`.
    flat = ER.Derivation._flat_edges(m)
    function _testhelper_segments_agree(mask, kept)
        s = ER.Mechanisms.rapid_equilibrium_segments(m, mask)
        segments, extras = _testhelper_re_segments(g, kept)
        present = findall(!=(0), s.segment_of)
        names = m.metabolites
        _testhelper_offsets(f) = Dict(names[n] => s.offsets[n, f]
                                      for n in eachindex(names) if s.offsets[n, f] > 0)
        s.segments == segments && [_testhelper_offsets(f) for f in present] == extras
    end
    _testhelper_segments_agree(trues(length(m.steps)), trues(length(flat))) || return false
    _testhelper_is_allosteric(x) || return true
    inactive = ER.Mechanisms.steps_running_in_allosteric_conformations(x).inactive
    kept = _testhelper_inactive_reach(g)
    inactive[flat] == kept && _testhelper_segments_agree(inactive, kept)
end

# The census's first populations (`_testhelper_census_population`,
# `_testhelper_census_allosteric`), enumerated once for the three testsets that read them
# and released after them.
let population = _testhelper_census_population(),
    allosteric = _testhelper_census_allosteric()
@testset "a graph's canonical form, pinned as text, read back from its Sig and raw" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test spec.mechanism == _testhelper_build(spec.raw)
    end
    # The specs' canonical text, a blank line between specs, pinned as a fixture.
    text = reduce((a, b) -> [a; ""; b],
                  [["# $(spec.name)";
                    _testhelper_canonical_text(_testhelper_build(spec.raw))]
                   for spec in MECHANISM_TEST_SPECS])
    fixture = joinpath(@__DIR__, "..", "fixtures", "canonical_text_golden.txt")
    if !isfile(fixture)
        mkpath(dirname(fixture)); write(fixture, join(text, "\n") * "\n")
        @warn "bootstrapped canonical text golden fixture; commit it" fixture
    end
    @test text == readlines(fixture)
    # Free E comes first, also where another form's name sorts before E.
    raw = _testhelper_walk(:(begin
        substrates: A
        products: P
        steps: begin
            E + A <--> EA
            EA <--> EP
            EP <--> D + P
            D <--> E
        end
    end))
    x = _testhelper_build(raw)
    @test x.catalytic.enzyme_forms[1:2] == [:E, :D]
    @test _testhelper_graph(x).forms[1:2] == ["E", "D"]
    # The census's first populations: each graph decodes from its Sig and builds again
    # from the forms and steps it writes back. The count: the six reactions' init
    # mechanisms, then the bi-bi reaction's children.
    @test length(population) == 5508 + 5471
    @test all(_testhelper_decode(ER.Derivation._sig_of(x)) == x for x in population)
    @test all(_testhelper_build(_testhelper_written(x)) == x for x in population)
end

@testset "a sub-graph's segments and inactive edges are what its steps give" begin
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_subgraphs_agree(spec.mechanism)
    end
    @test isempty([i for (i, x) in enumerate(population)
                   if !_testhelper_subgraphs_agree(x)])
end

@testset "a graph's DSL text reads back as the graph" begin
    # `show` prints the macro text that builds a mechanism and `DSL.parse_mechanism`
    # reads it back, for the specs, plain and allosteric, and the census's first
    # populations, plain and allosteric.
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_reads_back(spec.mechanism)
        @test ER.DSL.parse_mechanism(sprint(show, spec.mechanism)) == spec.mechanism
    end
    # The allosteric children of the uni-uni, ping-pong and two-site reactions.
    @test length(allosteric) == 5 + 885 + 5
    @test isempty([i for (i, x) in enumerate(population) if !_testhelper_reads_back(x)])
    @test isempty([i for (i, x) in enumerate(allosteric) if !_testhelper_reads_back(x)])
end
end

@testset "cycles, distances and a group's lead orient steps alike in every writing" begin
    # Isomerizations a turnover cycle orients, exchanges that run away from free E, and an
    # isomerization that keeps its direction relative to its group's lead, each case in two
    # writings, against each case's pinned canonical text. In the tied fixture neither
    # exchange lies on a cycle: each runs from its form nearer free E, EI and FI,
    # taking up J.
    ss_text = ["substrates: S | products: P",
               "forms: E, ES, F, FP, FS",
               "1: E + S <--> ES | group 1",
               "2: ES <--> FS | group 2",
               "3: F <--> E | group 3",
               "4: F + P <--> FP | group 4",
               "5: FS <--> FP | group 5"]
    tied_text = ["substrates: A | products: P | modifiers: I, J",
                 "forms: E, EA, EAI, EI, EJ, EJP, EP, F, FI, FJ",
                 "1: E + A <--> EA | group 1",
                 "2: E + I ⇌ EI | group 2",
                 "3: EA <--> EP | group 3",
                 "4: EI + A ⇌ EAI | group 4",
                 "5: EI + J <--> EJ + I | group 5",
                 "6: EJ + P ⇌ EJP | group 6",
                 "7: F <--> E | group 7",
                 "8: F + P <--> EP | group 8",
                 "9: F + I ⇌ FI | group 2",
                 "10: FI <--> EI | group 7",
                 "11: FI + J <--> FJ + I | group 5"]
    allosteric_text = ["substrates: A | products: P | modifiers: I, J",
                       "forms: E, EA, EI, EJ, EP, F, FI, FJ",
                       "1: E + A <--> EA | group 1 :: EqualAI",
                       "2: E + I ⇌ EI | group 2 :: EqualAI",
                       "3: EA <--> EP | group 3 :: EqualAI",
                       "4: EI + J <--> EJ + I | group 4 :: NonequalAI",
                       "5: F ⇌ E | group 5 :: EqualAI",
                       "6: F + P <--> EP | group 6 :: EqualAI",
                       "7: F + I ⇌ FI | group 2 :: EqualAI",
                       "8: FI + J <--> FJ + I | group 4 :: NonequalAI",
                       "oligomeric_state: 1"]
    _testhelper_tied_steps(free, bound, exchange) = :(begin
        E + A <--> EA
        EA <--> EP
        EP <--> F + P
        ($free, $bound)
        (E + I ⇌ EI, F + I ⇌ FI)
        EI + A ⇌ EAI
        EJ + P ⇌ EJP
        $exchange
    end)
    blocks = [
        :(begin
            substrates: S
            products: P
            steps: begin
                E + S <--> ES
                FS <--> ES
                FS <--> FP
                FP <--> F + P
                F <--> E
            end
        end),
        :(begin
            substrates: S
            products: P
            steps: begin
                ES <--> FS
                E + S <--> ES
                F + P <--> FP
                FP <--> FS
                E <--> F
            end
        end)]
    for (free, bound, exchange) in (
            (:(F <--> E), :(FI <--> EI),
             :((EI + J <--> EJ + I, FI + J <--> FJ + I))),
            (:(E <--> F), :(EI <--> FI),
             :((FJ + I <--> FI + J, EJ + I <--> EI + J))))
        push!(blocks, :(begin
            substrates: A
            products: P
            modifiers: I, J
            steps: $(_testhelper_tied_steps(free, bound, exchange))
        end))
    end
    texts = [ss_text, ss_text, tied_text, tied_text]
    for (raw, text) in zip(_testhelper_walk.(blocks), texts)
        @test _testhelper_canonical_text(_testhelper_build(raw)) == text
    end
    for exchange in (:((EI + J <--> EJ + I, FI + J <--> FJ + I)),
                     :((FJ + I <--> FI + J, EJ + I <--> EI + J)))
        raw = _testhelper_walk(:(begin
            substrates: A
            products: P
            modifiers: I, J
            steps: begin
                E + A <--> EA                  :: EqualAI
                EA <--> EP                   :: EqualAI
                EP <--> F + P                  :: EqualAI
                F ⇌ E                            :: EqualAI
                (E + I ⇌ EI, F + I ⇌ FI)     :: EqualAI
                $exchange                        :: NonequalAI
            end
        end))
        @test _testhelper_canonical_text(_testhelper_build(raw)) == allosteric_text
    end
end

"""Whether `x` survives its round trips: its printed text reads back as `x` and prints
again as that text, two reads of the text hash alike, and `x` decodes from the Sig of its
compiled law as itself."""
function _testhelper_round_trips(x)
    text = sprint(show, x)
    a, b = ER.DSL.parse_mechanism(text), ER.DSL.parse_mechanism(text)
    a == x && sprint(show, a) == text && hash(a) == hash(b) == hash(x) &&
        _testhelper_decode(ER.Derivation._sig_of(x)) == x
end

@testset "storage round trips" begin
    # The specs, plain and allosteric, and the bi-bi init set, the enumeration golden set
    # (`bi_bi_init_golden.txt`).
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_round_trips(spec.mechanism)
    end
    init = _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn))
    @test length(init) == 239
    @test isempty([i for (i, x) in enumerate(init) if !_testhelper_round_trips(x)])
end

"""Whether the steps each conformation of `am` runs
(`steps_running_in_allosteric_conformations`) are those its tags do not exclude
(`_testhelper_excludes`), on the forms they join to free E (`_testhelper_exclusion`). True
for a plain mechanism, whose one conformation runs every step."""
function _testhelper_conformations_hold(am::ER.Mechanism)
    cm = am.catalytic
    M = ER.Mechanisms
    (; active, inactive) = M.steps_running_in_allosteric_conformations(am)
    _testhelper_is_allosteric(am) || return active == trues(length(cm.steps))
    all(((1, active), (2, inactive))) do (c, runs)
        excluded = BitVector([_testhelper_excludes(am, k, c) for k in eachindex(cm.steps)])
        component = _testhelper_exclusion(am, c).component
        runs == BitVector([!excluded[k] && component[e.from] == 1
                           for (k, e) in enumerate(cm.steps)])
    end
end

@testset "the steps a conformation runs are what its states give" begin
    # Every construction path: the constructors (the specs' writings, the bi-bi seeds)
    # and the Sig decoder (the specs' compiled laws).
    @testset "$(spec.name)" for spec in MECHANISM_TEST_SPECS
        @test _testhelper_conformations_hold(_testhelper_build(spec.raw))
        @test _testhelper_conformations_hold(
            _testhelper_decode(ER.Derivation._sig_of(spec.mechanism)))
    end
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products: P[C], Q[N]
    end
    @test all(_testhelper_conformations_hold,
              _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(rxn)))
end

@testset "forms and steps are the storage" begin
    m = _testhelper_build(_testhelper_walk(:(begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
        end
    end)))
    c = m.catalytic
    @test (c.metabolites, c.enzyme_forms, c.substrates, c.products) ==
          ([:S, :P], [:E, :EP, :ES], [1], [2])
    # The release of P is stored as the binding it reverses, and sorts first.
    @test [e.taken_up for e in c.steps[1:2]] == [[2], [1]]
    @test (c.steps[2].from, c.steps[2].to, c.steps[2].given_off) == (1, 3, Int[])
    @test c.steps[1].rapid_equilibrium && c.steps[2].rapid_equilibrium &&
          !c.steps[3].rapid_equilibrium
    @test [s.kinetic_group for s in c.steps] == [1, 2, 3]
end

@testset "the constructor builds a spec's graph in under 30 µs, the median over specs" begin
    # Each spec's time per call is the best of five batches; the gate is on their median.
    # The spec's goal is about 2 µs and the median is about 5 µs locally; shared CI
    # runners take about 5× as long (23 µs on Ubuntu), so the gate is 30 µs.
    per_call = map(MECHANISM_TEST_SPECS) do spec
        _testhelper_build(spec.raw)
        minimum(@elapsed(for _ in 1:200
                             _testhelper_build(spec.raw)
                         end) / 200 for _ in 1:5)
    end
    @test median(per_call) < 30e-6
end

"""The bytes one construction of `x` allocates, from its own canonical fields, after a
warm-up call: `Mechanism(CatalyticScheme(metabolites, enzyme_forms, substrates, products,
steps), allostery)` with the allostery of `x`."""
function _testhelper_construction_bytes(x)
    cm, a = x.catalytic, x.allosteric
    _testhelper_construct() =
        ER.Mechanism(ER.Mechanisms.CatalyticScheme(cm.metabolites, cm.enzyme_forms,
                                                   cm.substrates, cm.products, cm.steps),
                     ER.Mechanisms.Allostery(a.allosteric_states, a.oligomeric_state,
                                             a.regulatory_sites))
    _testhelper_construct()
    @allocated _testhelper_construct()
end

@testset "constructor allocation budget" begin
    # Each construction allocates at most 1.3 times the bytes it allocated at bf49d6e1 on
    # Julia 1.13 (6,672, 22,264 and 25,912 bytes): a mechanism stores the facts it computes
    # once, and they may not make every construction in a search meaningfully dearer.
    budget(bytes) = ceil(Int, 1.3 * bytes)
    michaelis_menten = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> ES
            ES <--> E + P
        end
    end
    @test _testhelper_construction_bytes(michaelis_menten) <= budget(6_672)
    # The largest seed of the bi-bi init set (A, B → P, Q): 9 forms, 13 steps.
    random_bi_bi = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB, EP + A ⇌ EAP)
            (E + B ⇌ EB, EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            (E + P ⇌ EP, EA + P ⇌ EAP, EQ + P ⇌ EPQ)
            (E + Q ⇌ EQ, EB + Q ⇌ EBQ, EP + Q ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    @test _testhelper_construction_bytes(random_bi_bi) <= budget(22_264)
    # Hexokinase: 14 steps, two :OnlyA groups and a regulatory site.
    hk = only(s for s in MECHANISM_TEST_SPECS if s.name == "HK").mechanism
    @test _testhelper_construction_bytes(hk) <= budget(25_912)
end

# ── What the stored graph holds, and step text ────────────────────────────────────────

const ER = EnzymeRates

"""The names of the metabolites the list `list` of catalytic scheme `c` indexes, in
order."""
_testhelper_names(c, list) = c.metabolites[list]

@testset "graph accessors: random-order bi-bi" begin
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            (E + B ⇌ EB, EA + B ⇌ EAB)
            EAB <--> EPQ
            EPQ <--> EQ + P
            EQ <--> E + Q
        end
    end
    c = m.catalytic

    # The metabolites: substrates, then products, each by name. The forms: free E, then
    # by name.
    @testset "species" begin
        @test c.metabolites == [:A, :B, :P, :Q]
        @test (c.substrates, c.products) == ([1, 2], [3, 4])
        @test c.enzyme_forms == [:E, :EA, :EAB, :EB, :EPQ, :EQ]
    end

    # Step numbers and directions follow canonical form: each step is stored the way the
    # graph gives it, so a product release is stored as the binding it reverses, and the
    # steps are sorted by their stored ends.
    @testset "steps" begin
        @test length(c.steps) == 7
        @test [(e.from, e.to) for e in c.steps] ==
              [(1, 2), (1, 4), (1, 6), (2, 3), (3, 5), (4, 3), (6, 5)]
        @test [_testhelper_names(c, e.taken_up) for e in c.steps] ==
              [[:A], [:B], [:Q], [:B], Symbol[], [:A], [:P]]
        @test all(isempty(e.given_off) for e in c.steps)
        @test [e.rapid_equilibrium for e in c.steps] ==
              [true, true, false, true, false, true, false]
    end

    # Kinetic groups are numbered by their first step in step order.
    @testset "groups" begin
        @test [e.kinetic_group for e in c.steps] == [1, 2, 3, 2, 4, 1, 5]
    end
end

@testset "graph accessors: allosteric uni-uni" begin
    M = ER.Mechanisms
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
    m = am.catalytic
    @test m isa ER.Mechanisms.CatalyticScheme

    # Canonical form numbers the forms E, EP, ES and stores the steps E + P → EP,
    # E + S → ES, ES → EP; each step is its own kinetic group, in that order.
    @testset "species" begin
        @test m.enzyme_forms == [:E, :EP, :ES]
        @test m.metabolites == [:S, :P]
    end

    @testset "steps" begin
        @test [(e.from, e.to) for e in m.steps] == [(1, 2), (1, 3), (3, 2)]
        @test [_testhelper_names(m, e.taken_up) for e in m.steps] ==
              [[:P], [:S], Symbol[]]
        @test all(isempty(e.given_off) for e in m.steps)
        @test [e.rapid_equilibrium for e in m.steps] == [true, true, false]
    end

    @testset "groups" begin
        @test [e.kinetic_group for e in m.steps] == [1, 2, 3]
    end

    @testset "allosteric parts" begin
        # Each tag follows its group to the canonical group number.
        @test [am.allosteric.allosteric_states[g] for g in 1:3] ==
              [:EqualAI, :NonequalAI, :OnlyA]
        @test am.allosteric.oligomeric_state == 2
        @test am.allosteric.regulatory_sites == [M.RegulatorySite([:R], 2, [:OnlyI])]
    end
end

@testset "graph accessors: a second free form and a Theorell–Chance step" begin
    # Uni-uni through a second free form F: S binds straight into the product-bound form.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S <--> FP
            FP ⇌ F + P
            E <--> F
        end
    end
    # Canonical form numbers the forms E, F, FP. It stores the isomerization from F, where
    # only products enter or leave, to E, where only substrates do: against its writing.
    c = m.catalytic
    @test c.enzyme_forms == [:E, :F, :FP]
    @test [(e.from, e.to) for e in c.steps] == [(1, 3), (2, 1), (2, 3)]
    @test [_testhelper_names(c, e.taken_up) for e in c.steps] == [[:S], Symbol[], [:P]]

    # A Theorell–Chance step takes up B and gives off P.
    tc = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EQ + P
            EQ <--> E + Q
        end
    end
    # Canonical form numbers the forms E, EA, EQ and stores the steps E + A → EA,
    # E + Q → EQ, EA + B → EQ + P.
    c = tc.catalytic
    @test [(e.from, e.to) for e in c.steps] == [(1, 2), (1, 3), (2, 3)]
    @test (_testhelper_names(c, c.steps[3].taken_up),
           _testhelper_names(c, c.steps[3].given_off)) == ([:B], [:P])
end

@testset "graph accessors: a ping-pong's forms are names" begin
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP
            EB_res_plusA_minusP <--> E + Q
        end
    end
    # A form's residual is in its name alone: free E first, then the forms by name, `_`
    # after the capitals.
    @test m.catalytic.enzyme_forms ==
          [:E, :EA, :EB_res_plusA_minusP, :E_res_plusA_minusP]
end

@testset "graph accessors: a reactant bound as a modifier is taken up as itself" begin
    # Substrate inhibition: a second S binds ES as a modifier.
    m = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES + S ⇌ ESSmod
            ES <--> EP
            EP ⇌ E + P
        end
    end
    c = m.catalytic
    # The graph holds S once, as a substrate; the binding of the second S takes it up into
    # the form named ESSmod, which tells the two bindings apart.
    @test c.metabolites == [:S, :P]
    @test c.enzyme_forms == [:E, :EP, :ES, :ESSmod]
    @test (c.steps[4].from, c.steps[4].to) == (3, 4)
    @test c.steps[4].taken_up == c.steps[2].taken_up == [1]
    @test metabolites(m) == [:S, :P]
end

@testset "step_text: steps as the macro writes them" begin
    M = ER.Mechanisms
    m = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            (E + A ⇌ EA, EB + A ⇌ EAB)
            (E + B ⇌ EB, EA + B ⇌ EAB)
            EAB <--> EPQ
            EPQ <--> EQ + P
            EQ <--> E + Q
        end
    end
    # Each step the way the graph gives it, a release as the binding it reverses.
    @test [M.step_text(m.catalytic, k) for k in eachindex(m.catalytic.steps)] ==
          ["E + A ⇌ EA", "E + B ⇌ EB", "E + Q <--> EQ", "EA + B ⇌ EAB",
           "EAB <--> EPQ", "EB + A ⇌ EAB", "EQ + P <--> EPQ"]
    # A ping-pong: each half-reaction's release of a product is stored as its uptake.
    pingpong = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA <--> E_res_plusA_minusP + P
            E_res_plusA_minusP + B <--> EB_res_plusA_minusP
            EB_res_plusA_minusP <--> E + Q
        end
    end
    @test [M.step_text(pingpong.catalytic, k)
           for k in eachindex(pingpong.catalytic.steps)] ==
          ["E + A <--> EA", "E + Q <--> EB_res_plusA_minusP",
           "E_res_plusA_minusP + P <--> EA",
           "E_res_plusA_minusP + B <--> EB_res_plusA_minusP"]
    inhibited = @enzyme_mechanism begin
        substrates: S
        products: P
        steps: begin
            E + S ⇌ ES
            ES + S ⇌ ESSmod
            ES <--> EP
            EP ⇌ E + P
        end
    end
    @test M.step_text(inhibited.catalytic, 4) == "ES + S ⇌ ESSmod"
end

# ── Segments and cycles of a step mask ────────────────────────────────────────────────

const ER = EnzymeRates

"""Whether a cycle of the steps of graph `m` that `on` keeps takes up some metabolite net:
for some metabolite, the uptake of the kept steps is no difference of form potentials,
which holds when appending the uptake column to their incidence matrix raises its rank."""
function _testhelper_net_uptake_cycle(m, on::BitVector)
    kept = filter(k -> on[k], collect(eachindex(m.steps)))
    incidence = zeros(length(kept), length(m.enzyme_forms))
    for (r, k) in enumerate(kept)
        a, b = m.steps[k].from, m.steps[k].to
        incidence[r, a] -= 1
        incidence[r, b] += 1
    end
    metabolites = unique(reduce(vcat, [[e.taken_up; e.given_off] for e in m.steps]))
    any(metabolites) do x
        uptake = Int[count(==(x), m.steps[k].taken_up) -
                     count(==(x), m.steps[k].given_off) for k in kept]
        rank([incidence uptake]) > rank(incidence)
    end
end

@testset "turnover cycle: every subset of steps of five mechanisms" begin
    M = ER.Mechanisms
    # Each mechanism with the cycles of its steps that turn the reaction over, written as
    # step texts (`_testhelper_edge_text`). A set of steps holds a turnover cycle exactly
    # when it holds all the steps of one of these.
    cases = [
        (@enzyme_mechanism(begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES <--> E + P
            end
        end), [["E + S ⇌ ES", "E + P <--> ES"]]),
        # Uni-uni through a second free form F.
        (@enzyme_mechanism(begin
            substrates: S
            products: P
            steps: begin
                E + S <--> FP
                FP ⇌ F + P
                E <--> F
            end
        end), [["E + S <--> FP", "F + P ⇌ FP", "F <--> E"]]),
        # Random-order bi-bi: the square of the two binding orders takes up A and B once
        # each way, so it turns nothing over; each binding order with the chemistry and
        # the releases does.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                (E + A ⇌ EA, EB + A ⇌ EAB)
                (E + B ⇌ EB, EA + B ⇌ EAB)
                EAB <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
            end
        end), [["E + A ⇌ EA", "EA + B ⇌ EAB", "EAB <--> EPQ", "EQ + P <--> EPQ",
                "E + Q <--> EQ"],
               ["E + B ⇌ EB", "EB + A ⇌ EAB", "EAB <--> EPQ", "EQ + P <--> EPQ",
                "E + Q <--> EQ"]]),
        # Substrate inhibition: the dead-end binding of a second S lies on no cycle.
        (@enzyme_mechanism(begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ ES
                ES + S ⇌ ESSmod
                ES <--> EP
                EP ⇌ E + P
            end
        end), [["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"]]),
        # A Theorell–Chance step beside the merged complex EPQ: the triangle of the
        # two routes from EA to EQ takes up B and gives off P once each way.
        (@enzyme_mechanism(begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ EA
                EA + B <--> EPQ
                EPQ <--> EQ + P
                EQ <--> E + Q
                EA + B <--> EQ + P
            end
        end), [["E + A ⇌ EA", "EA + B <--> EPQ", "EQ + P <--> EPQ", "E + Q <--> EQ"],
               ["E + A ⇌ EA", "EA + B <--> EQ + P", "E + Q <--> EQ"]]),
    ]
    for (m, cycles) in cases
        n = length(eachindex(m.catalytic.steps))
        at = [[_testhelper_edge(m, text) for text in cycle] for cycle in cycles]
        for bits in 0:(2^n - 1)
            on = BitVector([isodd(bits >> (k - 1)) for k in 1:n])
            expected = any(ks -> all(k -> on[k], ks), at)
            @test (findall(on), M.has_turnover_cycle(m.catalytic, on)) ==
                  (findall(on), expected)
        end
    end
end

@testset "turnover cycle: MECHANISM_TEST_SPECS" begin
    M = ER.Mechanisms
    for spec in MECHANISM_TEST_SPECS
        mech = spec.mechanism
        m = mech.catalytic
        n = length(eachindex(m.steps))
        # Every mechanism turns the reaction over, and none at rapid equilibrium alone.
        @test (spec.name, M.has_turnover_cycle(m, trues(n))) == (spec.name, true)
        @test (spec.name, M.has_turnover_cycle(m, _testhelper_flags(m))) ==
              (spec.name, false)
        subsets = [[Int[]]; [[i] for i in 1:n]; [[i, j] for i in 1:n for j in (i + 1):n]]
        masks = [BitVector(k in s for k in 1:n) for s in subsets]
        _testhelper_is_allosteric(mech) &&
            append!(masks, M.steps_running_in_allosteric_conformations(mech))
        for on in masks
            @test (spec.name, findall(on), M.has_turnover_cycle(m, on)) ==
                  (spec.name, findall(on), _testhelper_net_uptake_cycle(m, on))
        end
    end
end

@testset "turnover cycle: malformed input" begin
    M = ER.Mechanisms
    # Hand-written graphs whose cycles do not each run the reaction a whole number of
    # times: a cycle takes up A and Q and gives off nothing, or takes up B alone. The
    # turnover test weighs each step by its uptake of the first substrate, A, so it finds
    # the first cycle and not the second. The guard rejects both graphs.
    aq = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            EAB <--> EQ + P
            EQ <--> E + Q
            EA + Q <--> EAQ
            EAQ <--> E
        end
    end
    on = BitVector(k in [_testhelper_edge(aq, t)
                         for t in ("E + A <--> EA", "EA + Q <--> EAQ", "EAQ <--> E")]
                   for k in eachindex(aq.catalytic.steps))
    @test M.has_turnover_cycle(aq.catalytic, on)
    b_only = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A <--> EA
            EA + B <--> EAB
            EAB <--> EQ + P
            EQ <--> E + Q
            E + B <--> EB
            EB <--> E
        end
    end
    on = BitVector(k in [_testhelper_edge(b_only, t)
                         for t in ("E + B <--> EB", "E <--> EB")]
                   for k in eachindex(b_only.catalytic.steps))
    @test !M.has_turnover_cycle(b_only.catalytic, on)
    # A mask with fewer flags than steps is out of bounds.
    @test_throws BoundsError M.has_turnover_cycle(aq.catalytic, trues(2))
end

@testset "segments: a steady-state step inside a segment is no step of the graph" begin
    M = ER.Mechanisms
    # Random-order bi-bi whose binding of B to free E alone is steady state: the other
    # three bindings put E, EA, EB and EAB in one rapid-equilibrium segment. The
    # guard rejects it (EB fills at zero A and B); the view reads any graph.
    m = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EB + A ⇌ EAB
            E + B <--> EB
            EA + B ⇌ EAB
            EAB <--> EPQ
            EPQ <--> EQ + P
            EQ <--> E + Q
        end
    end
    every = trues(length(eachindex(m.catalytic.steps)))
    text(k) = _testhelper_edge_text(m, k)
    g = M.rapid_equilibrium_segments(m.catalytic, every)
    # Canonical form numbers the forms E, EA, EAB, EB, EPQ, EQ.
    @test g.segments == [[1, 2, 3, 4], [5], [6]]
    # E + B <--> EB joins two forms of the first segment, so it is no step of the graph.
    @test text.(g.steps) == ["E + Q <--> EQ", "EAB <--> EPQ", "EQ + P <--> EPQ"]
end

@testset "has_nonzero_weight_cycle on plain arrays" begin
    M = ER.Mechanisms
    # The triangle 1 → 2 → 3 → 1 balances at weights 1, 1, -2 and not at 1, 1, -1; without
    # its third step it closes no cycle.
    from, to3 = [1, 2, 3], [2, 3, 1]
    buffers() = (zeros(Int, 3), zeros(Int, 3))
    @test !M.has_nonzero_weight_cycle(from, to3, [1, 1, -2], trues(3), buffers()...)
    @test M.has_nonzero_weight_cycle(from, to3, [1, 1, -1], trues(3), buffers()...)
    @test !M.has_nonzero_weight_cycle(from, to3, [1, 1, -1], BitVector([1, 1, 0]),
                                      buffers()...)
end

@testset "a sub-graph reads the forms and segments of the edges a mask keeps" begin
    # Random release of P and Q at rapid equilibrium between steady-state releases: the
    # segment {EP, EPQ, EQ} has every form holding a product.
    m = @_testhelper_unguarded @enzyme_mechanism begin
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
    end
    names(forms) = String.(m.catalytic.enzyme_forms[forms])
    @test names(eachindex(m.catalytic.enzyme_forms)) == ["E", "EP", "EPQ", "EQ", "ES"]
    # One edge per group, so the edges group by group are the edges in order.
    @test [_testhelper_edge_text(m, k) for k in ER.Derivation._flat_edges(m.catalytic)] ==
          ["E + P <--> EP", "E + Q <--> EQ", "E + S ⇌ ES", "EP + Q ⇌ EPQ", "EQ + P ⇌ EPQ",
           "ES <--> EPQ"]
    @test _testhelper_flags(m) == BitVector([0, 0, 1, 1, 1, 0])
    every = trues(length(m.catalytic.steps))
    @test ER.Mechanisms._present_species(m.catalytic, every) == trues(5)
    # Each segment breadth-first from its first form; offsets per metabolite (rows S, P,
    # Q, the metabolites in stored order) above the segment's lowest.
    s = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, every, _testhelper_flags(m))
    @test [names(seg) for seg in s.segments] == [["E", "ES"], ["EP", "EPQ", "EQ"]]
    @test s.segment_of == [1, 2, 2, 2, 1]
    @test s.offsets == [0 0 0 0 1; 0 1 1 0 0; 0 0 1 1 0]
    # The steady-state steps that join two segments, which are all of them here.
    g = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, every)
    @test g.steps == [1, 2, 6]
    @test (g.segment_of, g.segments, g.offsets) == (s.segment_of, s.segments, s.offsets)
    # Without the release of P from EPQ, EQ stands alone.
    no_p = BitVector([1, 1, 1, 1, 0, 1])
    @test ER.Mechanisms._present_species(m.catalytic, no_p) == trues(5)
    s = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, no_p, _testhelper_flags(m))
    @test [names(seg) for seg in s.segments] == [["E", "ES"], ["EP", "EPQ"], ["EQ"]]
    @test ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, no_p).steps == [1, 2, 6]
    # A mask keeps only the forms its edges touch; a form it leaves is in no segment.
    bindings = BitVector([1, 1, 0, 0, 0, 0])
    @test ER.Mechanisms._present_species(m.catalytic, bindings) ==
          BitVector([1, 1, 0, 1, 0])
    s = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, bindings,
                                                 _testhelper_flags(m))
    @test [names(seg) for seg in s.segments] == [["E"], ["EP"], ["EQ"]]
    @test s.segment_of == [1, 2, 0, 3, 0]
    @test ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, bindings).steps == [1, 2]
    # `re` picks the edges that join segments: with none, each present form stands alone.
    none = ER.Mechanisms.rapid_equilibrium_segments(m.catalytic, every, falses(6))
    @test length(none.segments) == 5
    @test none.steps == 1:6

    # The inactive conformation keeps the edges outside :OnlyA groups between forms that
    # such edges reach from free E: here EA and EAI lie behind the :OnlyA
    # binding of A.
    am = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        oligomeric_state: 1
        steps: begin
            E + A <--> EA                                   :: OnlyA
            EA <--> EP                                    :: OnlyA
            E + P ⇌ EP                                      :: EqualAI
            (EA + I ⇌ EAI, EP + I ⇌ EIP)          :: NonequalAI
        end
    end
    inactive = ER.Mechanisms.steps_running_in_allosteric_conformations(am).inactive
    @test [_testhelper_edge_text(am, k) for k in findall(inactive)] ==
          ["E + P ⇌ EP", "EP + I ⇌ EIP"]
    @test am.catalytic.enzyme_forms[
              findall(ER.Mechanisms._present_species(am.catalytic, inactive))] ==
          [:E, :EIP, :EP]
end

# ── The steps each conformation runs ──────────────────────────────────────────────────

@testset "steps_running_in_allosteric_conformations" begin
    M = ER.Mechanisms
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES          :: EqualAI
            ES <--> EP        :: OnlyA
            EP ⇌ E + P          :: EqualAI
            E + I ⇌ EI          :: OnlyI
            EI + P ⇌ EIP    :: EqualAI
        end
    end
    m = am.catalytic
    # Canonical form numbers the forms E, EI, EIP, EP, ES and stores the steps
    # E + I → EI, E + P → EP, E + S → ES, EI + P → EIP, ES → EP.
    @test [(e.from, e.to) for e in m.steps] == [(1, 2), (1, 4), (1, 5), (2, 3), (5, 4)]
    # The active conformation drops the :OnlyI binding of I, which strands EI and so
    # also the binding of P to it. The inactive conformation drops the :OnlyA chemistry
    # and keeps every form.
    @test M.steps_running_in_allosteric_conformations(am) ==
          (; active = BitVector([0, 1, 1, 0, 1]), inactive = BitVector([1, 1, 1, 1, 0]))
end

# ── A group's allosteric tag through canonical order ──────────────────────────────────

# Canonicalizing catalytic-step order in the Mechanism constructor must keep
# each catalytic kinetic group's allosteric-state tag bound to the SAME
# physical step. HK declares BOTH ATP binding and the catalytic conversion as
# :OnlyA; after construction the two :OnlyA tags must still sit on those two
# physical steps (not on some other group canonical reordering moved into their
# source slots): exactly one :OnlyA group is the ATP-binding step, and the other
# is the metabolite-free catalytic conversion.
@testset "Allosteric group tag stays bound to its step" begin
    hk = only(s for s in MECHANISM_TEST_SPECS if s.name == "HK")
    am = hk.mechanism
    cm = am.catalytic
    members = _testhelper_group_members(cm)
    onlyA_groups = [g for g in eachindex(members)
                    if am.allosteric.allosteric_states[g] === :OnlyA]
    @test length(onlyA_groups) == 2
    onlyA_firsts = [cm.steps[first(members[g])] for g in onlyA_groups]
    @test count(e -> EnzymeRates.Enumeration._is_binding(e) &&
                     cm.metabolites[only(e.taken_up)] === :ATP,
                onlyA_firsts) == 1
    @test count(!EnzymeRates.Enumeration._is_binding, onlyA_firsts) == 1
end

# ── Step text reads back as its step ──────────────────────────────────────────────────

@testset "each step a message names reads back as that step" begin
    pingpong = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA <--> EP_res_plusA_minusP
            EP_res_plusA_minusP ⇌ E_res_plusA_minusP + P
            E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP
            EB_res_plusA_minusP <--> EQ
            EQ ⇌ E + Q
        end
    end
    theorell_chance = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EQ + P
            EQ ⇌ E + Q
        end
    end
    merged = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            EA + B <--> EPQ
            EPQ ⇌ EQ + P
            EQ ⇌ E + Q
        end
    end
    copies = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ ES
            ES <--> EP
            EP ⇌ E + P
            (E + I ⇌ EI, ES + I ⇌ EIS)
            ES + P ⇌ EPmodS
        end
    end
    mechs = [[s.mechanism for s in MECHANISM_TEST_SPECS];
             [pingpong, theorell_chance, merged, copies]]
    for mech in mechs
        cm = mech.catalytic
        g = _testhelper_graph(mech)
        names = (substrates = g.substrates, products = g.products,
                 modifiers = g.modifiers,
                 regulators = unique!([l for site in mech.allosteric.regulatory_sites
                                       for l in site.ligands]))
        for (k, e) in enumerate(cm.steps)
            op, lhs, rhs = EnzymeRates.DSL._step_parts(
                Meta.parse(EnzymeRates.Mechanisms.step_text(cm, k)))
            forms = Symbol[]
            s = EnzymeRates.DSL._raw_step!(forms, op, lhs, rhs, e.kinetic_group, names)
            @test forms == cm.enzyme_forms[[e.from, e.to]] &&
                  s == EnzymeRates.Mechanisms.Step(1, 2, e.taken_up, e.given_off,
                                                   e.rapid_equilibrium, e.kinetic_group)
        end
    end
    line(m, text) = EnzymeRates.Mechanisms.step_text(m.catalytic, _testhelper_edge(m, text))
    @test line(pingpong, "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP") ==
          "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP"
    @test line(theorell_chance, "EA + B <--> EQ + P") == "EA + B <--> EQ + P"
    @test line(merged, "EA + B <--> EPQ") == "EA + B <--> EPQ"
    @test line(copies, "ES + P ⇌ EPmodS") == "ES + P ⇌ EPmodS"
    @test line(copies, "E + I ⇌ EI") == "E + I ⇌ EI"
end

# ── Canonical form and equality on written and enumerated mechanisms ──────────────────

@testset "Canonical-by-construction representation independence" begin
    # One kinetic group holding two bindings of the same metabolite S (context-
    # shared: S binds free E, and S binds A-bound E, with the same K), A's binding, the
    # chemistry and P's release. Written two ways: the statements in reverse order AND
    # the group's two steps swapped. Canonical-by-construction must collapse them to one
    # struct.
    m_orderA = @enzyme_mechanism begin
        substrates: S, A
        products: P
        steps: begin
            (E + S ⇌ ES, EA + S ⇌ EAS)
            E + A ⇌ EA
            EAS <--> EP
            EP ⇌ E + P
        end
    end
    m_orderB = @enzyme_mechanism begin
        substrates: S, A
        products: P
        steps: begin
            EP ⇌ E + P
            EAS <--> EP
            E + A ⇌ EA
            (EA + S ⇌ EAS, E + S ⇌ ES)
        end
    end
    @test m_orderA == m_orderB
    @test hash(m_orderA) == hash(m_orderB)

    # unique! collapses the duplicate orderings; the mechanisms are not mutated.
    # m_split groups the two binding steps separately, so it is structurally
    # distinct and survives alongside the collapsed orderA/orderB.
    m_split = @enzyme_mechanism begin
        substrates: S, A
        products: P
        steps: begin
            E + S ⇌ ES
            EA + S ⇌ EAS
            E + A ⇌ EA
            EAS <--> EP
            EP ⇌ E + P
        end
    end
    mechs = [m_orderA, m_split, m_orderB]
    snapshot = deepcopy(mechs)
    result = unique!(mechs)
    @test length(result) == 2
    @test all(r -> any(==(r), snapshot), result)
    @test m_orderA == snapshot[1]
    @test m_split  == snapshot[2]
    @test m_orderB == snapshot[3]
    @test (m_orderA.catalytic.enzyme_forms, m_orderA.catalytic.steps) ==
          (snapshot[1].catalytic.enzyme_forms, snapshot[1].catalytic.steps)
end

@testset "Allosteric mechanism TR equivalence" begin
# Three distinct tags on three groups: S binding (RE, :EqualAI),
# P binding (RE, :NonequalAI), iso (SS, :OnlyA). The :OnlyA sits on
# the chemical step (dropped from the Haldane graph), so the mechanism
# is Haldane-valid with no one-sided :OnlyA binding.
m_compiled = @enzyme_mechanism begin
    substrates: S
    products: P
    oligomeric_state: 2
    steps: begin
        E + S ⇌ ES             :: EqualAI
        E + P ⇌ EP             :: NonequalAI
        ES <--> EP           :: OnlyA
    end
end
# Group order is canonical; allosteric tags stay bound to their steps.
g = _testhelper_graph(m_compiled)
# The tag of the one group whose first step satisfies `pred`.
_testhelper_state_of(pred) = only(g.tags[s.group] for s in _testhelper_first_steps(g)
                                  if pred(s))
_testhelper_binds_in(names) = s -> _testhelper_is_binding(s) && only(s.consumed) in names
@test _testhelper_state_of(_testhelper_binds_in(g.substrates)) == :EqualAI
@test _testhelper_state_of(_testhelper_binds_in(g.products)) == :NonequalAI
@test _testhelper_state_of(!_testhelper_is_binding) == :OnlyA
end

# ─── canonical by construction ─────────────────────────────────────────

@testset "Allosteric mechanism — canonical by construction" begin
# Site permutation with DISTINCT multiplicities. The constructor must permute
# allosteric_states alongside the catalytic groups and produce the same
# regulatory_sites ordering (regulatory side) regardless of input site order.
am_ab = @enzyme_mechanism begin
    substrates: S; products: P
    oligomeric_state: 2
    steps: begin
        E + P ⇌ EP    :: EqualAI
        E + S ⇌ ES    :: EqualAI
        ES <--> EP  :: EqualAI
    end
    regulatory_site(2): A::OnlyA
    regulatory_site(4): B::OnlyI
end
am_ba = @enzyme_mechanism begin   # sites swapped
    substrates: S; products: P
    oligomeric_state: 2
    steps: begin
        E + P ⇌ EP    :: EqualAI
        E + S ⇌ ES    :: EqualAI
        ES <--> EP  :: EqualAI
    end
    regulatory_site(4): B::OnlyI
    regulatory_site(2): A::OnlyA
end
_testhelper_assert_mechanism_invariants(am_ab)
_testhelper_assert_mechanism_invariants(am_ba)
@test am_ab == am_ba
end

# ─── mechanism dedup (unique!) ──────────────────────────────────────────
@testset "Dedup" begin

@testset "Mechanism — expansion-path overlap: dedup actually fires" begin
    # Run two rounds of expansion on the uni-uni init seeds
    # (Mechanism path), then unique!. Assert that the flat vector shrinks,
    # proving that two different expansion paths — two moves, or one move
    # reached through two parents — produced equivalent Mechanisms.
    init_mechs = collect(
        _testhelper_mechanisms(EnzymeRates.Enumeration.init_mechanisms(uni_uni_rxn)))
    pool = unique!(_testhelper_children(uni_uni_rxn, init_mechs))
    expanded = _testhelper_children(uni_uni_rxn, pool)
    pre = length(expanded)
    unique!(expanded)
    # dedup fired: two different expansion paths produced equivalent
    # Mechanisms, so the flat vector shrank.
    @test length(expanded) < pre
end

@testset "Mechanism — permuted groups collapse via canonicalization" begin
    # Two Mechanisms with the same physics but with their outer
    # kinetic-group order arbitrarily rearranged should collapse to one
    # after unique!. This exercises the constructor's outer-group
    # sort with a non-trivial permutation, confirming that any
    # permutation of the outer Vector canonicalizes back to the same
    # struct.
    m_seed = @enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ EA
            E + Q ⇌ EQ
            (EA + B ⇌ EAB, EQ + B ⇌ EBQ)
            (EA + P ⇌ EAP, EQ + P ⇌ EPQ)
            EAB <--> EPQ
        end
    end
    _testhelper_assert_mechanism_invariants(m_seed)
    n_groups = length(_testhelper_group_members(m_seed.catalytic))
    @assert n_groups >= 3 "bi-bi init seed must have ≥3 kinetic groups " *
        "to exercise a non-trivial permutation"
    # Cyclic rotation by one, [g1, g2, …, gN] → [g2, …, gN, g1]: group g written as group
    # g - 1 and group 1 last, the steps written group by group in that order.
    label = [mod1(e.kinetic_group - 1, n_groups) for e in m_seed.catalytic.steps]
    order = sortperm(label)
    m_rotated = _testhelper_plain_mechanism(
        m_seed.catalytic.metabolites, m_seed.catalytic.enzyme_forms,
        m_seed.catalytic.substrates, m_seed.catalytic.products,
        [EnzymeRates.Mechanisms.Step(e.from, e.to, e.taken_up, e.given_off,
                                     e.rapid_equilibrium, label[k])
         for (k, e) in enumerate(m_seed.catalytic.steps)][order])
    v = EnzymeRates.Mechanism[m_seed, m_rotated]
    unique!(v)
    @test length(v) == 1
end
end

@testset "has_turnover_cycle: a balanced rapid-equilibrium cycle runs no turnover" begin
    ER = EnzymeRates
    # Random-order binding of A and B at rapid equilibrium closes the cycle
    # E → EA → EAB ← EB ← E of weight 1 + 1 − 1 − 1 = 0: it runs no net reaction.
    # With the isomerization RE as well, E → EA → EAB → EPQ → EQ → E has
    # weight 4, one turnover, in the same block.
    m = @enzyme_mechanism begin
        substrates: A, B; products: P, Q
        steps: begin
            E + A ⇌ EA
            E + B ⇌ EB
            EA + B ⇌ EAB
            EB + A ⇌ EAB
            EAB <--> EPQ
            EQ + P ⇌ EPQ
            E + Q ⇌ EQ
        end
    end
    @test !ER.Mechanisms.has_turnover_cycle(m.catalytic, _testhelper_flags(m))
    @test ER.Mechanisms.has_turnover_cycle(m.catalytic, trues(length(m.catalytic.steps)))
end
