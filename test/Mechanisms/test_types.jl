# ABOUTME: Tests for the core EnzymeRates types (the mechanisms, the reaction, the
# ABOUTME: Parameter family) and their constructors, accessors, and equality.
using Serialization

const ER = EnzymeRates

# Michaelis–Menten with rapid-equilibrium bindings and a steady-state isomerization.
const _testhelper_re_mm = @enzyme_mechanism begin
    substrates: S
    products:   P
    steps: begin
        E + S ⇌ E(S)
        E(S) <--> E(P)
        E(P) ⇌ E + P
    end
end

# The raw steps of S ⇌ P: S binding E, E(S) isomerizing to E(P) at steady state, and P
# binding E, the bindings at rapid equilibrium.
function _testhelper_uniuni_steps()
    met(x) = ER.Mechanisms.Met(x, false)
    E, ES, EP = ER.Mechanisms.RawForm(:E, ER.Mechanisms.Met[]),
        ER.Mechanisms.RawForm(:E, [met(:S)]), ER.Mechanisms.RawForm(:E, [met(:P)])
    (bind = ER.Mechanisms.RawStep(E, ES, [met(:S)], ER.Mechanisms.Met[], true),
     iso = ER.Mechanisms.RawStep(ES, EP, ER.Mechanisms.Met[], ER.Mechanisms.Met[], false),
     rel = ER.Mechanisms.RawStep(E, EP, [met(:P)], ER.Mechanisms.Met[], true))
end

# The raw mechanism over substrate S and product P whose steps are `steps`, step `j` in
# kinetic group `j`.
_testhelper_raw_uniuni(steps) = ER.Mechanisms.RawMechanism(
    [:S], [:P], Symbol[], collect(steps), collect(eachindex(steps)))

# The allosteric mechanism over those raw uni-uni `steps`, group `j` tagged `tags[j]`,
# with catalytic multiplicity `n` and the raw regulatory sites `sites`.
_testhelper_tagged_uniuni(steps, tags, n, sites) = ER.AllostericMechanism(
    ER.Mechanisms.RawAllostericMechanism(_testhelper_raw_uniuni(steps), tags, n, sites))

@testset "Types" begin
    @testset "EnzymeMechanism struct + accessors" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end

        @test _testhelper_substrates(m) == [:S]
        @test _testhelper_products(m) == [:P]
        @test isempty(_testhelper_regulators(m))
        @test ER.metabolites(m) == (:S, :P)
        # Steps are canonicalized at construction, so compare content
        # order-independently.
        mech = ER.Mechanism(m)
        @test Set(_testhelper_edge_text(mech, k) for k in eachindex(mech.edges)) ==
              Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])
        @test ER.Mechanisms.ngroups(mech) == 3
        @test Set(_testhelper_graph(m).forms) == Set(["E", "ES", "EP"])

        # Shared kinetic-group: two steps in group 1 (regulator R binds
        # both E and E(S) sharing one K).
        m2 = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: R
            steps: begin
                (E + R ⇌ E(R), E(S) + R ⇌ E(S, R))
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        # Steps are canonicalized; the two R-binding steps share one kinetic
        # group → 5 steps, 4 groups (one 2-step group binding R).
        mech2 = ER.Mechanism(m2)
        @test ER.Mechanisms.ngroups(mech2) == 4
        shared = only(ks for ks in ER.Mechanisms._group_members(mech2) if length(ks) == 2)
        @test all(k -> mech2.table[ER.Mechanisms._bound_entry(mech2.edges[k])] ==
                  ER.Mechanisms.Met(:R, true),
                  shared)
    end

    @testset "metabolites() lists substrates, then products, then regulators" begin
        # A substrate, a product and a regulator, so the lift lists all three roles in
        # order — coverage the plain S/P accessor tests do not reach.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            modifiers: I
            steps: begin
                E + S ⇌ E(S)
                E(S) + I ⇌ E(S, I)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test ER.metabolites(m) == (:S, :P, :I)
    end

    @testset "a mechanism names its substrates, products and modifiers" begin
        # A modifier-only name and a reactant's modifier copy. The concentration order
        # is substrates, products, then the other modifiers, each by name, as the
        # compiled mechanism lists them.
        em = @enzyme_mechanism begin
            substrates: B, A
            products: Q, P
            modifiers: J, A, I
            steps: begin
                E + A ⇌ E(A)
                E(A) + B ⇌ E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) ⇌ E + Q
                E + J ⇌ E(J)
                E(A) + I ⇌ E(A, I)
                E(Q) + A::Mod ⇌ E(A::Mod, Q)
            end
        end
        m = ER.Mechanism(em)
        @test ER.Mechanisms.substrates(m) == [:A, :B]
        @test ER.Mechanisms.products(m) == [:P, :Q]
        @test ER.Mechanisms.modifiers(m) == [:A, :I, :J]
        @test ER.metabolites(m) == [:A, :B, :P, :Q, :I, :J]
        @test ER.metabolites(em) == (:A, :B, :P, :Q, :I, :J)
        # A reaction's competitive inhibitor that no step binds is no concentration of
        # its init mechanisms.
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: K
        end
        on_rxn = only(ER.Enumeration.init_mechanisms(rxn))
        @test ER.Mechanisms.modifiers(on_rxn) == Symbol[]
        @test ER.metabolites(on_rxn) == [:S, :P]
        # An allosteric mechanism lists its site ligands last, site by site; a name that
        # is both a modifier and a ligand is listed once.
        am = ER.AllostericMechanism(@allosteric_mechanism begin
            substrates: S
            products: P
            modifiers: I
            allosteric_regulators: Y, X, I
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)     :: EqualAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
                E + I ⇌ E(I)     :: EqualAI
            end
            regulatory_site(2): Y::OnlyI
            regulatory_site(2): X::OnlyA, I::OnlyI
        end)
        @test (ER.Mechanisms.substrates(am), ER.Mechanisms.products(am),
               ER.Mechanisms.modifiers(am)) == ([:S], [:P], [:I])
        @test ER.metabolites(am) == [:S, :P, :I, :X, :Y]
        @test Tuple(ER.metabolites(am)) == ER.metabolites(ER.compile_mechanism(am))
    end

    @testset "AllostericEnzymeMechanism lift + DSL" begin
        # Single-ligand :EqualAI reg site is allowed (degenerate but valid —
        # the enumerator won't emit it, but users may write it for teaching).
        @test @allosteric_mechanism(begin
            substrates: S
            products:   P
            allosteric_regulators: I
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: NonequalAI
            end
            regulatory_site(2): I::EqualAI
        end) isa ER.AllostericEnzymeMechanism

        # A catalytic group may be :OnlyI only when its steps bind modifiers alone.
        @test_throws ("AllostericMechanism: E + S ⇌ E(S) is tagged :OnlyI but does more " *
                      "than bind a modifier.") @allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: OnlyI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: NonequalAI
            end
        end
        # Nor when they change what else their forms hold: that is catalysis.
        why = " Only the active conformation catalyzes, so :OnlyI fits only a group " *
              "whose steps take up or give off modifiers alone. Tag its group :NonequalAI."
        @test_throws ("AllostericMechanism: E(S) + X ⇌ E(P, X) is tagged :OnlyI but does " *
                      "more than bind a modifier." * why) @allosteric_mechanism begin
            substrates: S
            products:   P
            modifiers: X
            steps: begin
                E + S <--> E(S)          :: EqualAI
                E(S) <--> E(P)           :: EqualAI
                E + P <--> E(P)          :: EqualAI
                E(S) + X ⇌ E(P, X)       :: OnlyI
                E(P) + X <--> E(P, X)    :: EqualAI
            end
        end
        copy_msg = "AllostericMechanism: E(S) + S::Mod ⇌ E(P, S::Mod) is tagged :OnlyI " *
                   "but does more than bind a modifier." * why
        @test_throws copy_msg @allosteric_mechanism begin
            substrates: S
            products:   P
            modifiers: S
            steps: begin
                E + S <--> E(S)                  :: EqualAI
                E(S) <--> E(P)                   :: EqualAI
                E + P <--> E(P)                  :: EqualAI
                E(S) + S::Mod ⇌ E(P, S::Mod)     :: OnlyI
                E(P) + S::Mod <--> E(P, S::Mod)  :: EqualAI
            end
        end
        onlyi = @allosteric_mechanism begin
            substrates: S
            products:   P
            modifiers: I
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E + P   :: EqualAI
                E + I ⇌ E(I)      :: OnlyI
            end
        end
        @test :OnlyI in ER.Mechanisms.cat_allo_states(ER.AllostericMechanism(onlyi))

        # Build via DSL
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I

            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)     :: EqualAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
            regulatory_site(2): I::OnlyI
        end
        @test ER.AllostericMechanism(m).catalytic_multiplicity == 2
        # Catalysis (iso step) is :OnlyA, both bindings :EqualAI. Group order is
        # canonical, so identify the tagged group by its step, not its position.
        am_c = ER.AllostericMechanism(m)
        members = ER.Mechanisms._group_members(am_c.catalytic)
        onlyA_g = only(g for g in eachindex(members)
                       if ER.Mechanisms.tag(am_c, g) === :OnlyA)
        @test !ER.Mechanisms._is_binding(am_c.catalytic.edges[first(members[onlyA_g])])
        @test all(ER.Mechanisms.tag(am_c, g) === :EqualAI
                  for g in eachindex(members) if g != onlyA_g)
        @test ER.Mechanisms.allosteric_regulators(am_c) == [:I]
        @test ER.Mechanisms.allo_states(only(ER.Mechanisms.regulatory_sites(am_c))) ==
              [:OnlyI]
    end

    @testset "a mechanism prints as the text of the macro that builds it" begin
        # The text prints alike from the compiled and the concrete mechanism, and the
        # macro builds the same mechanism from it.
        _testhelper_builds_back(m) = (t = sprint(show, m);
                                      t == sprint(show, ER.Mechanisms._concrete(m)) &&
                                          eval(Meta.parse(t)) === m)
        # One step per line in stored order, a release written as the binding it
        # reverses, a form with nothing bound as its bare label.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        @test sprint(show, m) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + P <--> E(S)
                    E + S <--> E(S)
                end
            end"""
        @test _testhelper_builds_back(m)

        # A block with no steps holds no free E, so the macro builds no mechanism to print.
        @test_throws(ErrorException("Mechanism: no free enzyme form E; every mechanism " *
                                    "holds the free form E, with nothing bound and no " *
                                    "residual"),
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
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B <--> E(A, B)
                E(B) + A <--> E(A, B)
                E(A, B) <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        @test sprint(show, m_b) == """
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E + B <--> E(B)
                    E + Q <--> E(Q)
                    E(A) + B <--> E(A, B)
                    E(A, B) <--> E(P, Q)
                    E(B) + A <--> E(A, B)
                    E(Q) + P <--> E(P, Q)
                end
            end"""
        @test _testhelper_builds_back(m_b)

        # Catalytic 3-cycle: the steps print in stored order, not chain order, and
        # rapid-equilibrium and steady-state steps keep their arrows.
        m_re = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test sprint(show, m_re) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    E + P ⇌ E(P)
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                end
            end"""
        @test _testhelper_builds_back(m_re)

        # Theorell–Chance: B binds and P leaves in one step, written on one line.
        m_tc = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        @test sprint(show, m_tc) == """
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E + Q <--> E(Q)
                    E(A) + B <--> E(Q) + P
                end
            end"""
        @test _testhelper_builds_back(m_tc)

        # A fused binding (B binds and chemistry runs in one step) prints B on its entry
        # side, since E(P, Q) does not name it.
        m_fb = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(P, Q)
                E(P, Q) <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        @test sprint(show, m_fb) == """
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    E + Q <--> E(Q)
                    E(A) + B <--> E(P, Q)
                    E(Q) + P <--> E(P, Q)
                end
            end"""
        @test _testhelper_builds_back(m_fb)

        # The modifiers: line lists the names the steps bind as modifiers, a reactant's
        # copy by its name; a kinetic group of several steps prints as a tuple, a
        # reactant's copy as X::Mod.
        @test !contains(sprint(show, _testhelper_re_mm), "modifiers")
        m_reg = @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: I
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                (E + I ⇌ E(I), E(P) + I ⇌ E(I, P))
            end
        end
        @test sprint(show, m_reg) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                modifiers: I
                steps: begin
                    (E + I ⇌ E(I), E(P) + I ⇌ E(I, P))
                    E + P ⇌ E(P)
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                end
            end"""
        @test _testhelper_builds_back(m_reg)
        m_copy = @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: S
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
                E(P) + S::Mod ⇌ E(P, S::Mod)
            end
        end
        @test sprint(show, m_copy) == """
            @enzyme_mechanism begin
                substrates: S
                products: P
                modifiers: S
                steps: begin
                    E + P ⇌ E(P)
                    E + S ⇌ E(S)
                    E(P) + S::Mod ⇌ E(P, S::Mod)
                    E(S) <--> E(P)
                end
            end"""
        @test _testhelper_builds_back(m_copy)

        # A ping-pong: a form that keeps a residual prints it, with or without bound
        # metabolites.
        m_pp = @enzyme_mechanism begin
            substrates: A, B
            products: P, Q
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(P; residual = A - P)
                E(; residual = A - P) + P ⇌ E(P; residual = A - P)
                E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                E(B; residual = A - P) <--> E(Q)
                E + Q ⇌ E(Q)
            end
        end
        @test sprint(show, m_pp) == """
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A ⇌ E(A)
                    E + Q ⇌ E(Q)
                    E(; residual = A - P) + B ⇌ E(B; residual = A - P)
                    E(; residual = A - P) + P ⇌ E(P; residual = A - P)
                    E(A) <--> E(P; residual = A - P)
                    E(B; residual = A - P) <--> E(Q)
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

        # An allosteric mechanism adds its regulators, its catalytic multiplicity, a tag
        # after each kinetic group and a line per regulatory site.
        m_allo = @allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            allosteric_regulators: I
            catalytic_multiplicity: 2
            steps: begin
                E + F6P ⇌ E(F6P)         :: EqualAI
                E(F6P) <--> E(F16BP)     :: EqualAI
                E(F16BP) ⇌ E + F16BP     :: EqualAI
            end
            regulatory_site(2): I::OnlyI
        end
        @test sprint(show, m_allo) == """
            @allosteric_mechanism begin
                substrates: F6P
                products: F16BP
                allosteric_regulators: I
                catalytic_multiplicity: 2
                steps: begin
                    E + F16BP ⇌ E(F16BP) :: EqualAI
                    E + F6P ⇌ E(F6P) :: EqualAI
                    E(F6P) <--> E(F16BP) :: EqualAI
                end
                regulatory_site(2): I::OnlyI
            end"""
        @test _testhelper_builds_back(m_allo)

        # A regulatory site with two ligands lists both with their tags.
        m_site = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I, J
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)     :: EqualAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
            regulatory_site(3): I::OnlyI, J::EqualAI
        end
        @test sprint(show, m_site) == """
            @allosteric_mechanism begin
                substrates: S
                products: P
                allosteric_regulators: I, J
                catalytic_multiplicity: 2
                steps: begin
                    E + P ⇌ E(P) :: EqualAI
                    E + S ⇌ E(S) :: EqualAI
                    E(S) <--> E(P) :: OnlyA
                end
                regulatory_site(3): I::OnlyI, J::EqualAI
            end"""
        @test _testhelper_builds_back(m_site)
    end

    @testset "the macro rejects a rapid-equilibrium segment with no bottom form" begin
        # Random product release from E(P, Q) at rapid equilibrium, but release from E(P)
        # and E(Q) at steady state: the RE segment {E(P), E(P, Q), E(Q)} has weights
        # P·Q : Kq·P : Kp·Q, all zero at P = Q = 0. The RE approximation carries no
        # parameter for how E(P, Q) splits between E(P) and E(Q) there, so the rate is
        # 0/0 at zero products: class 3 rejects that face (`_BOTTOMLESS`), and the RE-speed
        # oracle sees the steady state there move with the release speeds.
        random_release = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E + P <--> E(P)
                E + Q <--> E(Q)
                E(P) + Q ⇌ E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                E(S) <--> E(P, Q)
            end
        end
        v = ER.Derivation._guard(random_release)
        @test (v.fault, v.face) == (ER.Derivation._BOTTOMLESS, 2)
        msg = ER.Derivation._guard_message(random_release, v)
        @test occursin("has no form whose weight keeps the segment's lowest power of " *
                       "each zeroed metabolite", msg)
        @test msg ==
              "Mechanism: at zero P and Q, the rapid-equilibrium segment " *
              "{E(P), E(P, Q), E(Q)} has no form whose weight keeps the segment's lowest " *
              "power of each zeroed metabolite, so the forms' shares depend on how those " *
              "concentrations vanish and the rate is 0/0 there. Make the group of " *
              "E(P) + Q ⇌ E(P, Q) steady-state."
        # The printed fix passes the guard.
        @test ER.Derivation._guard(@_testhelper_unguarded @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E + P <--> E(P)
                E + Q <--> E(Q)
                E(P) + Q <--> E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                E(S) <--> E(P, Q)
            end
        end) === ER.Derivation._WELL_POSED
        @test _testhelper_respeed_verdicts(random_release; rng = Random.MersenneTwister(1),
                                           law = :limit)[2] === :pattern
        @test_throws ErrorException(msg) @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E + P <--> E(P)
                E + Q <--> E(Q)
                E(P) + Q ⇌ E(P, Q)
                E(Q) + P ⇌ E(P, Q)
                E(S) <--> E(P, Q)
            end
        end

        # Mirror on the substrate side: the segment {E(A), E(B), E(A, B)} is 0/0 at
        # A = B = 0 (face 3).
        mirror = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B ⇌ E(A, B)
                E(B) + A ⇌ E(A, B)
                E + P ⇌ E(P)
                E(A, B) <--> E(P)
            end
        end
        v = ER.Derivation._guard(mirror)
        @test (v.fault, v.face) == (ER.Derivation._BOTTOMLESS, 3)
        @test _testhelper_respeed_verdicts(mirror; rng = Random.MersenneTwister(2),
                                           law = :limit)[3] === :pattern
        @test_throws ErrorException(
            ER.Derivation._guard_message(mirror, v)) @enzyme_mechanism begin
            substrates: A, B
            products: P
            steps: begin
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B ⇌ E(A, B)
                E(B) + A ⇌ E(A, B)
                E + P ⇌ E(P)
                E(A, B) <--> E(P)
            end
        end

        # The same segment with A's modifier copy binding E(B) at rapid equilibrium:
        # E(A::Mod, B) joins it, but every form still carries A or B (the copy's
        # concentration is A's), so the segment stays bottomless, whichever order the
        # groups are written in.
        copy_last = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P
            modifiers: A
            steps: begin
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B ⇌ E(A, B)
                E(B) + A ⇌ E(A, B)
                E(B) + A::Mod ⇌ E(A::Mod, B)
                E + P ⇌ E(P)
                E(A, B) <--> E(P)
            end
        end
        copy_first = @_testhelper_unguarded @enzyme_mechanism begin
            substrates: A, B
            products: P
            modifiers: A
            steps: begin
                E(B) + A::Mod ⇌ E(A::Mod, B)
                E(A, B) <--> E(P)
                E + P ⇌ E(P)
                E(B) + A ⇌ E(A, B)
                E(A) + B ⇌ E(A, B)
                E + B <--> E(B)
                E + A <--> E(A)
            end
        end
        v = ER.Derivation._guard(copy_last)
        @test (v.fault, v.face) == (ER.Derivation._BOTTOMLESS, 3)
        @test ER.Derivation._guard(copy_first) == v
        @test copy_first == copy_last
        for m in (copy_last, copy_first)
            segment = ER.Mechanisms.re_segments(
                m, trues(length(m.edges)), ER.Mechanisms.rapid_equilibrium_steps(m)).bottom
            @test [ER.Mechanisms._form_name(m, f) for f in segment] ==
                  ["EA", "EAB", "EB", "EAmodB"]
        end
        @test_throws ErrorException(
            ER.Derivation._guard_message(copy_last, v)) @enzyme_mechanism begin
            substrates: A, B
            products: P
            modifiers: A
            steps: begin
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B ⇌ E(A, B)
                E(B) + A ⇌ E(A, B)
                E(B) + A::Mod ⇌ E(A::Mod, B)
                E + P ⇌ E(P)
                E(A, B) <--> E(P)
            end
        end

        # Ordered release: the segment {E(P), E(P, Q)} has weights Kq : Q, so E(P) is its
        # bottom form and the rate is defined at zero products; every face agrees with
        # the RE-speed oracle.
        m_ordered = @enzyme_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)
                E + P <--> E(P)
                E(P) + Q ⇌ E(P, Q)
                E(S) <--> E(P, Q)
            end
        end
        @test m_ordered isa EnzymeMechanism
        @test all(==(:agree),
                  _testhelper_respeed_verdicts(m_ordered; rng = Random.MersenneTwister(3)))

        # Mixed abortive complex: the segment {E(S), E(P), E(S, P)} has weights
        # K·S : K′·P : S·P, zero only at S = P = 0, a corner outside the guard's faces
        # where no turnover is possible anyway. Accepted.
        m_abortive = @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E + P <--> E(P)
                E(S) + P ⇌ E(S, P)
                E(P) + S ⇌ E(S, P)
                E(S) <--> E(P)
            end
        end
        @test m_abortive isa EnzymeMechanism
        @test all(==(:agree),
                  _testhelper_respeed_verdicts(m_abortive; rng = Random.MersenneTwister(4)))
    end

    @testset "a free form that holds the enzyme it receives is a second closed class" begin
        # Every step that takes up the substrate A is :OnlyA, and the conformations flip at
        # free E, so the inactive conformation keeps E and E(P) alone: a second free form
        # F, which binds B and joins E only through F + B <--> E(B), holds no inactive
        # enzyme. In the active conformation, with every substrate absent, or A, or B, F
        # holds the enzyme it receives, together with E(B) while B is present, and never
        # returns it to E: a second closed class, which class 3 rejects (`_NONUNIQUE`) and
        # the RE-speed oracle sees.
        s1 = @_testhelper_unguarded @allosteric_mechanism begin
            substrates: A, B
            products: P
            modifiers: A
            steps: begin
                E + A <--> E(A)                   :: OnlyA
                E(A) + B <--> E(A, B)             :: EqualAI
                E(A, B) <--> E(P)                 :: OnlyA
                E(P) <--> E + P                   :: EqualAI
                E + A::Mod ⇌ E(A::Mod)            :: OnlyA
                E(A::Mod) + B ⇌ E(A::Mod, B)      :: EqualAI
                E(B) + A::Mod ⇌ E(A::Mod, B)      :: EqualAI
                F + B <--> E(B)                   :: EqualAI
            end
        end
        @test _testhelper_face_faults(s1)[3:5] == fill(ER.Derivation._NONUNIQUE, 3)
        v = ER.Derivation._guard(s1)
        @test (v.fault, v.face) == (ER.Derivation._NONUNIQUE, 3)
        @test _testhelper_respeed_verdicts(s1; rng = Random.MersenneTwister(12),
                                           law = :limit)[3:5] == fill(:nonunique, 3)
        @test_throws ErrorException(
            ER.Derivation._guard_message(s1, v)) @allosteric_mechanism begin
            substrates: A, B
            products: P
            modifiers: A
            steps: begin
                E + A <--> E(A)                   :: OnlyA
                E(A) + B <--> E(A, B)             :: EqualAI
                E(A, B) <--> E(P)                 :: OnlyA
                E(P) <--> E + P                   :: EqualAI
                E + A::Mod ⇌ E(A::Mod)            :: OnlyA
                E(A::Mod) + B ⇌ E(A::Mod, B)      :: EqualAI
                E(B) + A::Mod ⇌ E(A::Mod, B)      :: EqualAI
                F + B <--> E(B)                   :: EqualAI
            end
        end
        # Two variants keep F and pass every face: F <--> E tagged :OnlyA; or F <--> E
        # untagged with the modifier binding of A untagged, so that no :OnlyA binding has
        # its two sides joined by untagged steps.
        f_onlya = @allosteric_mechanism begin
            substrates: A, B
            products: P
            modifiers: A
            steps: begin
                E + A <--> E(A)                   :: OnlyA
                E(A) + B <--> E(A, B)             :: EqualAI
                E(A, B) <--> E(P)                 :: OnlyA
                E(P) <--> E + P                   :: EqualAI
                E + A::Mod ⇌ E(A::Mod)            :: OnlyA
                E(A::Mod) + B ⇌ E(A::Mod, B)      :: EqualAI
                E(B) + A::Mod ⇌ E(A::Mod, B)      :: EqualAI
                F + B <--> E(B)                   :: EqualAI
                F <--> E                          :: OnlyA
            end
        end
        copy_untagged = @allosteric_mechanism begin
            substrates: A, B
            products: P
            modifiers: A
            steps: begin
                E + A <--> E(A)                   :: OnlyA
                E(A) + B <--> E(A, B)             :: EqualAI
                E(A, B) <--> E(P)                 :: OnlyA
                E(P) <--> E + P                   :: EqualAI
                E + A::Mod ⇌ E(A::Mod)            :: EqualAI
                E(A::Mod) + B ⇌ E(A::Mod, B)      :: EqualAI
                E(B) + A::Mod ⇌ E(A::Mod, B)      :: EqualAI
                F + B <--> E(B)                   :: EqualAI
                F <--> E                          :: EqualAI
            end
        end
        # Both oracles agree on every face: the RE-speed oracle (class 3) and the ε-path
        # (admission of the :OnlyA tags).
        for (i, variant) in enumerate((f_onlya, copy_untagged))
            @test all(_testhelper_admission_agreement(
                ER.AllostericMechanism(variant); rng = Random.MersenneTwister(12 + i)))
        end
    end

    @testset "a second free form is admitted unless a tag strands it" begin
        # The conformations flip at free E alone. A second free form joined to E by a
        # steady-state step is an ordinary form, whatever its label (R1, and its twin
        # labelled D, which sorts before E). Where tags strand it, admission decides:
        # every step at E :OnlyA leaves F fed at zero product with a dead way back
        # (`_PUMPED`); a tagged isomerization alone joining F to E makes F an island
        # (`_ISLAND`). Each rejection has an admitted variant: the S binding untagged
        # (AD2d), or the catalysis tagged in place of the isomerization.
        r1 = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: NonequalAI
                F + P ⇌ E(P)      :: EqualAI
                F <--> E          :: NonequalAI
            end
        end
        r1_d = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: NonequalAI
                D + P ⇌ E(P)      :: EqualAI
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
        r2 = @_testhelper_unguarded @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: OnlyA
                E(S) <--> E(P)    :: EqualAI
                F + P ⇌ E(P)      :: EqualAI
                F <--> E          :: OnlyA
            end
        end
        v = ER.Derivation._guard(r2)
        @test (v.fault, v.face) == (ER.Derivation._PUMPED, 2)
        @test _testhelper_epath_verdict(r2; rng = Random.MersenneTwister(23),
                                        face = Set([:P])) === :direction_dependent
        @test_throws ErrorException(
            ER.Derivation._guard_message(r2, v)) @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: OnlyA
                E(S) <--> E(P)    :: EqualAI
                F + P ⇌ E(P)      :: EqualAI
                F <--> E          :: OnlyA
            end
        end
        ad2d = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                F + P ⇌ E(P)      :: EqualAI
                F <--> E          :: OnlyA
            end
        end
        @test _testhelper_epath_verdict(ER.AllostericMechanism(ad2d);
                                        rng = Random.MersenneTwister(24)) === :converges
        r3 = @_testhelper_unguarded @allosteric_mechanism begin
            substrates: S
            products: P
            modifiers: I
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E + P ⇌ E(P)      :: EqualAI
                E <--> F          :: OnlyA
                F + I ⇌ F(I)      :: EqualAI
            end
        end
        r4 = @_testhelper_unguarded @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: NonequalAI
                E + P ⇌ E(P)      :: EqualAI
                F ⇌ E             :: OnlyA
            end
        end
        for (i, m) in enumerate((r3, r4))
            v = ER.Derivation._guard(m)
            @test v.fault === ER.Derivation._ISLAND
            @test _testhelper_epath_verdict(m; rng = Random.MersenneTwister(25 + i)) in
                  (:direction_dependent, :wrong_limit)
        end
        msg = ER.Derivation._guard_message(r3, ER.Derivation._guard(r3))
        @test_throws ErrorException(msg) @allosteric_mechanism begin
            substrates: S
            products: P
            modifiers: I
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E + P ⇌ E(P)      :: EqualAI
                E <--> F          :: OnlyA
                F + I ⇌ F(I)      :: EqualAI
            end
        end
        msg = ER.Derivation._guard_message(r4, ER.Derivation._guard(r4))
        @test_throws ErrorException(msg) @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: NonequalAI
                E + P ⇌ E(P)      :: EqualAI
                F ⇌ E             :: OnlyA
            end
        end
        r3_base = @allosteric_mechanism begin
            substrates: S
            products: P
            modifiers: I
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: OnlyA
                E + P ⇌ E(P)      :: EqualAI
                E <--> F          :: EqualAI
                F + I ⇌ F(I)      :: EqualAI
            end
        end
        r4_base = @allosteric_mechanism begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: OnlyA
                E + P ⇌ E(P)      :: EqualAI
                F ⇌ E             :: EqualAI
            end
        end
        for (i, m) in enumerate((r3_base, r4_base))
            @test _testhelper_epath_verdict(
                ER.AllostericMechanism(m); rng = Random.MersenneTwister(27 + i)) ===
                :converges
        end
        # F in E's rapid-equilibrium segment, untagged, is an ordinary form too (A5).
        @test @allosteric_mechanism(begin
            substrates: S
            products: P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: NonequalAI
                E + P ⇌ E(P)      :: EqualAI
                F ⇌ E             :: NonequalAI
            end
        end) isa AllostericEnzymeMechanism
    end

    @testset "the inactive conformation's segments need a bottom form too" begin
        # The :OnlyA P group holds the fused binding E(Q) + P ⇌ Estar(S) beside the plain
        # binding E + P ⇌ E(P). In the active conformation that binding joins E(P),
        # E(P, Q) and E(Q) to free E at rapid equilibrium, and E is the segment's bottom
        # form. The inactive conformation lacks it and reaches those forms only through
        # the steady-state chemistry and the steady-state Q binding at E, so its segment
        # {E(P), E(P, Q), E(Q)} holds a product in every form: class 3 finds it bottomless
        # at zero products. The same untagged steps join the two sides of the tagged
        # binding E + P ⇌ E(P), so its affinity cannot vanish there: admission rejects the
        # tag set first (`_NO_LIMIT`), and the ε-path oracle finds no path.
        am = @_testhelper_unguarded @allosteric_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)                               :: EqualAI
                E(S) <--> E(P, Q)                          :: EqualAI
                E(P, Q) ⇌ E(P) + Q                         :: EqualAI
                E(P, Q) ⇌ E(Q) + P                         :: EqualAI
                (E + P ⇌ E(P), E(Q) + P ⇌ Estar(S))        :: OnlyA
                E + Q <--> E(Q)                            :: EqualAI
            end
        end
        @test _testhelper_face_faults(am)[2] === ER.Derivation._BOTTOMLESS
        v = ER.Derivation._guard(am)
        @test v.fault === ER.Derivation._NO_LIMIT
        @test ER.Mechanisms.step_text(am.catalytic, v.edge) == "E + P ⇌ E(P)"
        @test _testhelper_epath_verdict(am; rng = Random.MersenneTwister(5)) === :no_path
        @test_throws ErrorException(
            ER.Derivation._guard_message(am, v)) @allosteric_mechanism begin
            substrates: S
            products: P, Q
            steps: begin
                E + S ⇌ E(S)                               :: EqualAI
                E(S) <--> E(P, Q)                          :: EqualAI
                E(P, Q) ⇌ E(P) + Q                         :: EqualAI
                E(P, Q) ⇌ E(Q) + P                         :: EqualAI
                (E + P ⇌ E(P), E(Q) + P ⇌ Estar(S))        :: OnlyA
                E + Q <--> E(Q)                            :: EqualAI
            end
        end
    end

    @testset "AllostericEnzymeMechanism lift validators" begin
        # Wrong-length cat_allo_states (4 entries for 3 kinetic groups) → error
        cm = ER.Mechanism(_testhelper_re_mm)
        @test_throws ErrorException ER.AllostericMechanism(
            cm, [:NonequalAI, :NonequalAI, :NonequalAI, :OnlyA], 2,
            ER.Mechanisms.RegulatorySite[])

        # A ligand on two regulatory sites → error: both sites would render its Kreg name
        @test_throws(ErrorException(
            "AllostericMechanism: ligand A appears in two distinct regulatory sites; " *
            "rendered Kreg names would collide. Same-ligand-two-sites is not enumerated."),
            ER.AllostericMechanism(
                cm, [:EqualAI, :EqualAI, :EqualAI], 2,
                [ER.Mechanisms.RegulatorySite([:A], n, [:OnlyA]) for n in (2, 4)]))

        # Invalid allo state value → error
        @test_throws ErrorException @allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: NotAState
            end
        end

        # Reg site with all-:EqualAI ligands is allowed (degenerate but valid).
        @test @allosteric_mechanism(begin
            substrates: S
            products:   P
            allosteric_regulators: I, J
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: NonequalAI
            end
            regulatory_site(2): I::EqualAI, J::EqualAI
        end) isa ER.AllostericEnzymeMechanism

        # Invalid reg-site ligand allo state → error
        @test_throws ErrorException @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: NonequalAI
                E(P) ⇌ E + P      :: NonequalAI
            end
            regulatory_site(2): I::NotAState
        end
    end

    @testset "an allosteric mechanism prints every tag and every site" begin
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I, J
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)     :: NonequalAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
            regulatory_site(2): I::NonequalAI
            regulatory_site(2): J::OnlyI
        end
        s = sprint(show, m)
        @test s == """
            @allosteric_mechanism begin
                substrates: S
                products: P
                allosteric_regulators: I, J
                catalytic_multiplicity: 2
                steps: begin
                    E + P ⇌ E(P) :: EqualAI
                    E + S ⇌ E(S) :: NonequalAI
                    E(S) <--> E(P) :: OnlyA
                end
                regulatory_site(2): I::NonequalAI
                regulatory_site(2): J::OnlyI
            end"""
        am = ER.AllostericMechanism(m)
        # Every catalytic state appears in cat_allo_states line
        for g in eachindex(ER.Mechanisms.cat_allo_states(am))
            @test occursin(string(ER.Mechanisms.tag(am, g)), s)
        end
        # No :NonequalAI ligand silently hidden from reg-site display
        for site in ER.Mechanisms.regulatory_sites(am),
            (lig, state) in zip(ER.Mechanisms.ligands(site),
                                ER.Mechanisms.allo_states(site))
            @test occursin("$lig::$state", s)
        end
    end

    @testset "an allosteric kinetic group prints as one tagged tuple" begin
        # A mechanism with no regulatory site prints no allosteric_regulators: line, but
        # its catalytic multiplicity always.
        am = @allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_multiplicity: 2
            steps: begin
                (E + S ⇌ E(S), E(P) + S ⇌ E(P, S))  :: EqualAI
                E(S) <--> E(P)                      :: EqualAI
                E(P) ⇌ E + P                        :: EqualAI
            end
        end
        @test repr(am) == """
            @allosteric_mechanism begin
                substrates: S
                products: P
                catalytic_multiplicity: 2
                steps: begin
                    E + P ⇌ E(P) :: EqualAI
                    (E + S ⇌ E(S), E(P) + S ⇌ E(P, S)) :: EqualAI
                    E(S) <--> E(P) :: EqualAI
                end
            end"""
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
        @test rxn.allowed_catalytic_multiplicities == [1, 2]
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
        @test ER.Mechanisms.ligands(site) == [:A, :B]
        @test ER.Mechanisms.multiplicity(site) == 4
        @test ER.Mechanisms.allo_states(site) == [:OnlyA, :NonequalAI]

        # Mismatched ligand / allo_state length → error
        @test_throws ErrorException ER.Mechanisms.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA])

        # Multiplicity < 1 → error
        @test_throws ErrorException ER.Mechanisms.RegulatorySite(
            [lig_a], 0, [:OnlyA])

        # No ligands → error
        @test_throws ErrorException ER.Mechanisms.RegulatorySite(Symbol[], 2, Symbol[])

        # Invalid allo state → error
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
        @test ER.Mechanisms.ligands(site_reordered) == [lig_a, lig_b]
        @test ER.Mechanisms.allo_states(site_reordered) == [:OnlyA, :NonequalAI]
        @test site == site_reordered
        @test hash(site) == hash(site_reordered)
        # ...but the pairing of ligand and state does: B::OnlyA, A::NonequalAI is
        # a different site.
        site_swapped = ER.Mechanisms.RegulatorySite(
            [lig_b, lig_a], 4, [:OnlyA, :NonequalAI])
        @test ER.Mechanisms.allo_states(site_swapped) == [:NonequalAI, :OnlyA]
        @test site != site_swapped
    end


    @testset "Parameter family: edge constants and Kreg" begin
        # One RE constant and one forward/reverse rate pair serve every kind of edge, each
        # an edge of a mechanism's catalytic graph in one state.
        for T in (ER.Mechanisms.Krapid, ER.Mechanisms.Kfor, ER.Mechanisms.Krev)
            p = T(1, :None)
            @test p isa T
            @test p isa ER.Mechanisms.Parameter
            @test p == T(1, :None)
            @test p != T(1, :I) && p != T(2, :None)
        end

        # A regulator constant is its ligand's in one state; it needs no site, since a
        # ligand sits at one site.
        kr = ER.Mechanisms.Kreg(:A, :A)
        @test kr isa ER.Mechanisms.Parameter
        @test kr == ER.Mechanisms.Kreg(:A, :A)
        @test kr != ER.Mechanisms.Kreg(:A, :I) && kr != ER.Mechanisms.Kreg(:B, :A)
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
            "EnzymeReaction: $n is a reserved name; E_total, Keq, L, N, params, concs " *
            "and v are kept for a rate equation's constants and local names")
        for n in (:E_total, :Keq, :L, :N, :params, :concs, :v)
            @test_throws _testhelper_reserved(n) ER.EnzymeReaction(
                [M.Reactant(n, [:C => 1])], [P], noregs, Int[1])
            @test_throws _testhelper_reserved(n) ER.EnzymeReaction(
                [S], [P], [M.Regulator(n, :allosteric)], Int[1])
        end
        @test_throws _testhelper_reserved(:E_total) @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            allosteric_regulators: E_total
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

    @testset "EnzymeReaction takes one or more catalytic multiplicities, each ≥ 1" begin
        _testhelper_reaction(mults) =
            EnzymeReaction([ER.Mechanisms.Reactant(:S, [:C => 1])],
                           [ER.Mechanisms.Reactant(:P, [:C => 1])],
                           ER.Mechanisms.Regulator[], mults)
        @test_throws(ErrorException("EnzymeReaction: allowed_catalytic_multiplicities " *
                                    "must be a nonempty list of distinct Ints ≥ 1, got " *
                                    "[0, 1]"),
                     _testhelper_reaction([0, 1]))
        # An empty list would leave the search no multiplicity to assign.
        @test_throws(ErrorException("EnzymeReaction: allowed_catalytic_multiplicities " *
                                    "must be a nonempty list of distinct Ints ≥ 1, got " *
                                    "Int64[]"),
                     _testhelper_reaction(Int[]))
        # A multiplicity listed twice is no list a DSL text writes.
        @test_throws(ErrorException("EnzymeReaction: allowed_catalytic_multiplicities " *
                                    "must be a nonempty list of distinct Ints ≥ 1, got " *
                                    "[2, 2]"),
                     _testhelper_reaction([2, 2]))
        # A Bool is no multiplicity: the list is a Vector{Int}, and none is converted.
        @test_throws MethodError _testhelper_reaction(Bool[true])
    end

    @testset "pure conformational isomerization runs product-exit to substrate-entry" begin
        # Tier 2 (the Segel Iso Uni Uni case): the pure conformational F <--> E, tied
        # on Tier 1 and decided by entry_kind, runs F → E however it is written.
        s2 = @enzyme_mechanism begin
            substrates: A; products: P
            steps: begin
                E + A <--> E(A); E(A) <--> E(P); E(P) <--> F + P; E <--> F   # last step flipped
            end
        end
        m = ER.Mechanism(s2)
        iso = only(e for e in m.edges
                   if ER.Mechanisms._is_iso(e) &&
                       ER.Mechanisms._form_name(m, e.from) in ("E", "F"))
        @test ER.Mechanisms._form_name(m, iso.from) == "F"  # product-exit
        @test ER.Mechanisms._form_name(m, iso.to)   == "E"  # substrate-entry
    end

    @testset "Mechanism canonicalizes step order" begin
        (; bind, iso, rel) = _testhelper_uniuni_steps()

        m = ER.Mechanism(_testhelper_raw_uniuni([bind, iso, rel]))
        @test (ER.Mechanisms.substrates(m), ER.Mechanisms.products(m)) == ([:S], [:P])
        flat = ER.Mechanisms._flat_edges(m)
        @test Set(_testhelper_edge_text(m, k) for k in flat) ==
              Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])
        @test m.group[flat] == [1, 2, 3]
        @test length(m.edges) == 3
        # Canonical by construction: a permuted input yields identical storage.
        m_perm = ER.Mechanism(_testhelper_raw_uniuni([rel, bind, iso]))
        @test (m.forms, m.edges, m.group) == (m_perm.forms, m_perm.edges, m_perm.group)
        @test m == m_perm && hash(m) == hash(m_perm)
        # Steps and groups sort by content (`isless(::Edge)`), each form by its
        # content (`isless(::Form)`) with free E first. The residual form
        # E(; residual = A - P) binds nothing, so its steps come before those of E(A),
        # though its name, E_res_plusA_minusP, sorts after EA.
        pingpong = ER.Mechanism(@enzyme_mechanism begin
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
        end)
        @test [_testhelper_edge_text(pingpong, k)
               for k in ER.Mechanisms._flat_edges(pingpong)] ==
              ["E + A ⇌ EA",
               "E + Q ⇌ EQ",
               "E_res_plusA_minusP + B ⇌ EB_res_plusA_minusP",
               "E_res_plusA_minusP + P ⇌ EP_res_plusA_minusP",
               "EA <--> EP_res_plusA_minusP",
               "EB_res_plusA_minusP ⇌ EQ"]
        # The forms come in that order too: free E, then by content.
        @test [ER.Mechanisms._form_name(pingpong, f) for f in eachindex(pingpong.forms)] ==
              ["E", "E_res_plusA_minusP", "EA", "EB_res_plusA_minusP",
               "EP_res_plusA_minusP", "EQ"]
    end

    @testset "AllostericMechanism (non-parametric)" begin
        (; bind, iso, rel) = _testhelper_uniuni_steps()
        sites = [ER.Mechanisms.RawSite([:I], 1, [:OnlyI])]
        m = _testhelper_tagged_uniuni([bind, iso, rel], [:EqualAI, :OnlyA, :NonequalAI], 2,
                                      sites)

        @test (ER.Mechanisms.substrates(m), ER.Mechanisms.products(m)) == ([:S], [:P])
        # Steps canonicalized; allosteric tags stay bound to their steps.
        @test Set(_testhelper_edge_text(m, k) for k in eachindex(m.catalytic.edges)) ==
              Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])
        _testhelper_state_of(text) =
            ER.Mechanisms.tag(m, m.catalytic.group[_testhelper_edge(m, text)])
        @test _testhelper_state_of("E + S ⇌ ES") == :EqualAI
        @test _testhelper_state_of("ES <--> EP")  == :OnlyA
        @test _testhelper_state_of("E + P ⇌ EP")  == :NonequalAI
        @test ER.Mechanisms.catalytic_multiplicity(m) == 2
        @test ER.Mechanisms.regulatory_sites(m) ==
              [ER.Mechanisms.RegulatorySite([:I], 1, [:OnlyI])]
        @test ER.Mechanisms.ngroups(m.catalytic) == 3
        @test ER.Mechanisms.allosteric_regulators(m) == [:I]

        # Written in another order, each tag beside its step.
        m2 = _testhelper_tagged_uniuni([rel, iso, bind], [:NonequalAI, :OnlyA, :EqualAI], 2,
                                       sites)
        @test m == m2
        @test hash(m) == hash(m2)
    end

    @testset "AllostericMechanism validation errors" begin
        steps = collect(_testhelper_uniuni_steps())
        _testhelper_invalid_tag(tag) = ErrorException(
            "AllostericMechanism: catalytic group 1 has invalid allo state $tag (must be " *
            "one of (:OnlyA, :OnlyI, :EqualAI, :NonequalAI))")

        # :OnlyI fits only a group that binds modifiers alone: the active conformation
        # is the one that catalyzes.
        uni = ER.Mechanism(@enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end)
        onlyi_tags = fill(:EqualAI, ER.Mechanisms.ngroups(uni))
        onlyi_tags[uni.group[_testhelper_edge(uni, "E + S ⇌ ES")]] = :OnlyI
        @test_throws ErrorException(
            "AllostericMechanism: E + S ⇌ E(S) is tagged :OnlyI but does more than " *
            "bind a modifier. Only the active conformation catalyzes, so :OnlyI fits " *
            "only a group whose steps take up or give off modifiers alone. Tag its " *
            "group :NonequalAI.") ER.AllostericMechanism(uni, onlyi_tags, 1,
                                                          ER.Mechanisms.RegulatorySite[])

        # Length mismatch
        @test_throws(ErrorException("AllostericMechanism: cat_allo_states length 2 must " *
                                    "match the number of kinetic groups 3"),
                     _testhelper_tagged_uniuni(steps, [:EqualAI, :NonequalAI], 1,
                                               ER.Mechanisms.RawSite[]))

        # catalytic_multiplicity < 1
        @test_throws(ErrorException("AllostericMechanism: catalytic_multiplicity must be " *
                                    "≥ 1, got 0"),
                     _testhelper_tagged_uniuni(steps, [:EqualAI, :EqualAI, :EqualAI], 0,
                                               ER.Mechanisms.RawSite[]))

        # Unknown allo-state symbol
        @test_throws(_testhelper_invalid_tag(:Bogus),
                     _testhelper_tagged_uniuni(steps, [:Bogus, :EqualAI, :EqualAI], 1,
                                               ER.Mechanisms.RawSite[]))
    end

    @testset "AllostericMechanism(::AllostericEnzymeMechanism) converter" begin
        cm = _testhelper_re_mm
        aem = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: A, B
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: NonequalAI
                E(S) <--> E(P)    :: OnlyA
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(1): A::OnlyA, B::NonequalAI
        end

        am = ER.AllostericMechanism(aem)
        @test am isa ER.AllostericMechanism

        # Catalytic side lifted via Mechanism(CM())
        lifted = ER.Mechanism(cm)
        @test ER.Mechanisms.catalytic(am) == lifted
        @test (ER.Mechanisms.substrates(am), ER.Mechanisms.products(am),
               ER.Mechanisms.modifiers(am)) ==
              (ER.Mechanisms.substrates(lifted), ER.Mechanisms.products(lifted),
               ER.Mechanisms.modifiers(lifted))

        # cat_allo_states and multiplicity extracted from CS
        @test ER.Mechanisms.tag(am, 1) === :EqualAI
        @test ER.Mechanisms.tag(am, 2) === :NonequalAI
        @test ER.Mechanisms.tag(am, 3) === :OnlyA
        @test ER.Mechanisms.catalytic_multiplicity(am) == 2

        # Regulatory sites: ligand names
        sites = ER.Mechanisms.regulatory_sites(am)
        @test length(sites) == 1
        @test sites[1].ligands == [:A, :B]
        @test sites[1].multiplicity == 1
        @test sites[1].allo_states == [:OnlyA, :NonequalAI]

        # Idempotent: two calls give equal results
        @test ER.AllostericMechanism(aem) == am
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
        rep_bind = _testhelper_edge(am, "E + S ⇌ ES")
        @test ER.Mechanisms.name(ER.Mechanisms.Krapid(rep_bind, :None), am) === :K_ES_to_E_S
        @test ER.Mechanisms.name(ER.Mechanisms.Krapid(rep_bind, :I), am) === :K_I_ES_to_E_S

        rep_iso = _testhelper_edge(am, "ES <--> EP")
        @test ER.Mechanisms.name(ER.Mechanisms.Kfor(rep_iso, :None), am) === :k_ES_to_EP
        @test ER.Mechanisms.name(ER.Mechanisms.Krapid(rep_iso, :None), am) === :K_ES_to_EP
        @test ER.Mechanisms.name(ER.Mechanisms.Krev(rep_iso, :None), am) === :k_EP_to_ES

        lig = first(ER.Mechanisms.regulatory_sites(am)[1].ligands)
        @test ER.Mechanisms.name(ER.Mechanisms.Kreg(lig, :A), am) === :K_A_Rreg
        @test ER.Mechanisms.name(ER.Mechanisms.Kreg(lig, :I), am) === :K_I_Rreg

        # The chokepoint renders on the concrete mechanism only; a compiled type
        # is lifted first.
        @test_throws MethodError ER.Mechanisms.name(
            ER.Mechanisms.Krapid(rep_bind, :None), aem)
        @test_throws MethodError ER.Mechanisms.name(ER.Mechanisms.Kreg(lig, :A), aem)
    end

    @testset "_sig_of / _mechanism_from_sig roundtrip" begin
        # One fixture built from raw input: two substrates, two products, a
        # modifier-only name bound by a step, a covalent residual, a Theorell–Chance step,
        # and steps at rapid equilibrium and at steady state. It holds no reactant's
        # modifier copy.
        m = ER.Mechanism(ER.Mechanisms.raw_mechanism(:(begin
            substrates: A, B
            products: P, Q
            modifiers: I
            steps: begin
                E + A ⇌ E(A)
                E(A) <--> E(; residual = A - P) + P
                E(; residual = A - P) + B <--> E(Q)
                E(A) + B <--> E(Q) + P
                E + Q ⇌ E(Q)
                E + I ⇌ E(I)
            end
        end), false))

        sig = ER.Mechanisms._sig_of(m)
        @test sig isa Tuple{Vararg{Int}}
        @test ER.Mechanisms._mechanism_from_sig(sig...) == m   # roundtrip

        # The Sig is the compiled type's only parameter.
        em_type = ER.EnzymeMechanism{sig}
        @test em_type <: ER.EnzymeMechanism
        em_inst = em_type()
        @test ER.EnzymeMechanism(m) === em_inst

        # Roundtrip through the type-parameter form preserves everything.
        @test ER.Mechanism(em_inst) == m

        # A regulator that no step binds is no part of a mechanism, so the init mechanism
        # of a reaction that declares one compiles to the type of the reaction's without
        # it.
        without_u = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
        end
        with_u = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            competitive_inhibitors: U
        end
        @test ER.EnzymeMechanism(only(ER.Enumeration.init_mechanisms(with_u))) ===
              ER.EnzymeMechanism(only(ER.Enumeration.init_mechanisms(without_u)))
    end

    @testset "name(p::Parameter, m) chokepoint" begin
        m = ER.Mechanism(_testhelper_re_mm)
        # The binding of S, the isomerization and the binding of P, as edges of `m`.
        bind, iso, rel = (_testhelper_edge(m, text)
                          for text in ("E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"))

        # Structural naming: every step constant encodes its reaction's two sides;
        # a binding K reads in the release direction, iso params in the stored one.
        @test ER.Mechanisms.name(ER.Mechanisms.Krapid(bind, :None), m) === :K_ES_to_E_S
        @test ER.Mechanisms.name(ER.Mechanisms.Krapid(bind, :I),    m) === :K_I_ES_to_E_S
        @test ER.Mechanisms.name(ER.Mechanisms.Kfor(iso, :None), m) === :k_ES_to_EP
        @test ER.Mechanisms.name(ER.Mechanisms.Krev(iso, :None), m) === :k_EP_to_ES
        @test ER.Mechanisms.name(ER.Mechanisms.Krapid(rel, :None), m) === :K_EP_to_E_P

        # A binding's rate pair: the forward rate binds, the reverse rate releases
        @test ER.Mechanisms.name(ER.Mechanisms.Kfor(bind, :None), m) === :k_E_S_to_ES
        @test ER.Mechanisms.name(ER.Mechanisms.Krev(bind, :None), m) === :k_ES_to_E_S

        # I-state token on SS step
        @test ER.Mechanisms.name(ER.Mechanisms.Kfor(iso, :I), m) === :k_I_ES_to_EP
        @test ER.Mechanisms.name(ER.Mechanisms.Krev(iso, :I), m) === :k_I_EP_to_ES

        # The RE iso's equilibrium constant, named in the stored direction
        @test ER.Mechanisms.name(ER.Mechanisms.Krapid(iso, :None), m) === :K_ES_to_EP
        @test ER.Mechanisms.name(ER.Mechanisms.Krapid(iso, :I),    m) === :K_I_ES_to_EP

        # Same names resolve on the mechanism lifted back from EnzymeMechanism(m);
        # the compiled type itself is not a chokepoint argument.
        em = EnzymeMechanism(m)
        @test ER.Mechanisms.name(ER.Mechanisms.Krapid(bind, :None), ER.Mechanism(em)) ===
              :K_ES_to_E_S
        @test ER.Mechanisms.name(ER.Mechanisms.Kfor(iso, :None),
                                 ER.Mechanism(em)) === :k_ES_to_EP
        @test_throws MethodError ER.Mechanisms.name(ER.Mechanisms.Krapid(bind, :None), em)
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
        @test ER.Mechanisms.name(
            ER.Mechanisms.Krapid(_testhelper_edge(m, "E + S ⇌ ES"), :None), m) ===
              :K_ES_to_E_S
        @test ER.Mechanisms.name(
            ER.Mechanisms.Krapid(_testhelper_edge(m, "EP + S ⇌ EPS"), :None), m) ===
              :K_ES_to_E_S
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

@testset "a graph tells chemistry from a plain binding edge by edge" begin
    met(x) = ER.Mechanisms.Met(x, false)
    none = ER.Mechanisms.Met[]
    # A form of label `c` holding `xs` (names, or a `Met` for a modifier copy).
    _testhelper_form(xs...; c = :E, added = Symbol[], subtracted = Symbol[]) =
        ER.Mechanisms.RawForm(
            c, ER.Mechanisms.Met[x isa ER.Mechanisms.Met ? x : met(x) for x in xs], added,
            subtracted)
    E, EA, EAB = _testhelper_form(), _testhelper_form(:A), _testhelper_form(:A, :B)
    EPQ, EQ = _testhelper_form(:P, :Q), _testhelper_form(:Q)
    F = _testhelper_form(; added = [:A], subtracted = [:P])         # E(; residual = A - P)
    FB = _testhelper_form(:B; added = [:A], subtracted = [:P])
    A_copy = ER.Mechanisms.Met(:A, true)
    # Each written step: whether its edge is chemistry, whether it is a binding, and
    # whether it is stored turned (`missing`: left to the orientation tiers).
    cases = [
        # a plain binding, and a release stored as the binding it reverses, the label
        # free to change
        (ER.Mechanisms.RawStep(E, EA, [met(:A)], none, true), false, true, false),
        (ER.Mechanisms.RawStep(_testhelper_form(:A; c = :Estar), E, none, [met(:A)], true),
         false, true, true),
        # an isomerization
        (ER.Mechanisms.RawStep(EAB, EPQ, none, none, false), true, false, false),
        # chemistry and a release, stored as the fused binding of P it reverses
        (ER.Mechanisms.RawStep(EAB, EQ, none, [met(:P)], false), true, true, true),
        # a Theorell–Chance step
        (ER.Mechanisms.RawStep(EA, EQ, [met(:B)], [met(:P)], false), true, false, false),
        # two metabolites taken up in one step
        (ER.Mechanisms.RawStep(E, EAB, [met(:B), met(:A)], none, true), true, false,
         missing),
        # a fused binding: B taken up while chemistry runs
        (ER.Mechanisms.RawStep(EA, EPQ, [met(:B)], none, true), true, true, false),
        # a binding onto a residual form, and chemistry into it, stored as the fused
        # binding F + P → E(A)
        (ER.Mechanisms.RawStep(F, FB, [met(:B)], none, true), false, true, false),
        (ER.Mechanisms.RawStep(EA, F, none, [met(:P)], false), true, true, true),
        # the modifier copy of A binds apart from A
        (ER.Mechanisms.RawStep(E, _testhelper_form(A_copy), [A_copy], none, true), false,
         true, false),
        # an isomerization may change how often a name is bound: E(P, S) → E(P, P)
        (ER.Mechanisms.RawStep(_testhelper_form(:P, :S), _testhelper_form(:P, :P), none,
                               none, false),
         true, false, false)]
    c = ER.Mechanisms._canonicalize(ER.Mechanisms.RawMechanism(
        [:A, :B, :S], [:P, :Q], [:A], first.(cases), collect(eachindex(cases))))
    m = c.mechanism
    for (j, (_, chemistry, binding, turned)) in enumerate(cases)
        e = m.edges[c.edge_of[j]]
        @test ER.Mechanisms._is_chemistry(m, e) == chemistry
        @test ER.Mechanisms._is_binding(e) == binding
        turned === missing || @test c.flipped[j] == turned
    end
    @test ER.Mechanisms._is_iso(m.edges[c.edge_of[3]]) &&
          ER.Mechanisms._is_iso(m.edges[c.edge_of[11]])
    # The copy's binding takes up an entry of its own, into a form of its own.
    plain, copy_binding = m.edges[c.edge_of[1]], m.edges[c.edge_of[10]]
    @test m.table[ER.Mechanisms._bound_entry(copy_binding)] == A_copy &&
          m.table[ER.Mechanisms._bound_entry(plain)] == met(:A) && copy_binding.to !=
          plain.to
end

@testset "the macros reject a binding that changes a metabolite's role" begin
    _testhelper_role_change(from, to, m) = ErrorException(
        "Step $from → $to changes a metabolite's role: $to holds $from's " *
        "metabolites plus $m by name, not by role; write a modifier copy as " *
        "X::Mod on both sides")
    # The DSL spellings, with P declared as a modifier.
    for (line, from, to) in ((:(E + P::Mod ⇌ E(P)), :E, :EP),
                             (:(E(P) ⇌ E + P::Mod), :E, :EP),
                             (:(E(P) + P::Mod ⇌ E(P, P)), :EP, :EPP))
        @test_throws _testhelper_role_change(from, to, :P) eval(
            :(@enzyme_mechanism begin
                substrates: S
                products:   P
                modifiers: P
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                    $line
                end
            end))
    end
end

@testset "kinetic groups hold one kind of step with one flag" begin
    @testset "binding and Theorell-Chance step in one group" begin
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: A, B
                products: P, Q
                steps: begin
                    E + A <--> E(A)
                    (E(A) + B <--> E(A, B), E(A) + B <--> E(Q) + P)
                    E(Q) <--> E + Q
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds EA_B → EAB (SS) and EA_B → EQ_P (SS)",
                       err.msg)
    end

    @testset "RE and SS steps in one group" begin
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    (E + S ⇌ E(S), E(P) + S <--> E(P, S))
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds E_S → ES (RE) and EP_S → EPS (SS)", err.msg)
    end

    @testset "one reaction written RE and SS in one group" begin
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: S
                products: P
                steps: begin
                    (E + S ⇌ E(S), E + S <--> E(S))
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds E_S → ES (SS) and E_S → ES (RE)", err.msg)
    end

    @testset "bindings of two different metabolites in one group" begin
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: A, B
                products: P
                steps: begin
                    (E + A <--> E(A), E(A) + B <--> E(A, B))
                    E(A, B) <--> E + P
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds E_A → EA (SS) and EA_B → EAB (SS)", err.msg)
    end

    @testset "a substrate binding and its competitive-inhibitor copy in one group" begin
        # A and its modifier copy A::Mod are different metabolites, so their
        # bindings are different kinds of step and cannot share a constant.
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: A
                products: Q
                modifiers: A
                steps: begin
                    (E + A ⇌ E(A), E(Q) + A::Mod ⇌ E(A::Mod, Q))
                    E(A) <--> E(Q)
                    E(Q) <--> E + Q
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds E_A → EA (RE) and EQ_Amod → EAmodQ (RE)",
                       err.msg)
    end

    @testset "accepted: context-shared bindings of one metabolite" begin
        # R binds both free E and E(S) with the same K (non-competitive
        # inhibitor pattern): a legitimate shared binding of one metabolite.
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            modifiers: R
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E + P ⇌ E(P)
                (E + R ⇌ E(R), E(S) + R ⇌ E(S, R))
                E(R) + S ⇌ E(S, R)
            end
        end
        @test m isa EnzymeMechanism
    end

    @testset "accepted: mirrored isomerizations" begin
        # Two different conformational isomerizations, free and A-bound,
        # sharing one rate by a symmetry assumption: both are `:iso`, so
        # they are the same kind and may share a kinetic group.
        m = @enzyme_mechanism begin
            substrates: A
            products: P
            steps: begin
                E + A ⇌ E(A)
                (E <--> Estar, E(A) <--> Estar(A))
                Estar(A) <--> Estar(P)
                Estar(P) ⇌ Estar + P
            end
        end
        @test m isa EnzymeMechanism
    end
end

@testset "AllostericMechanism enforces the kinetic-group rules" begin
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)     :: NonequalAI
            E(S) <--> E(P)   :: EqualAI
            E(P) ⇌ E + P     :: EqualAI
        end
    end)
    cat = ER.Mechanisms.RawMechanism(ER.Mechanisms.catalytic(am))
    step_of(pred) = only(s for s in cat.steps if pred(s))
    bind_s = step_of(s -> s.consumed == [ER.Mechanisms.Met(:S, false)])
    bind_p = step_of(s -> s.consumed == [ER.Mechanisms.Met(:P, false)])
    iso = step_of(ER.Enumeration._is_iso)
    # The catalytic steps written as `steps`, labelled `group`, every group :EqualAI.
    _testhelper_rebuild(steps, group) = _testhelper_thrown() do
        ER.AllostericMechanism(ER.Mechanisms.RawAllostericMechanism(
            ER.Mechanisms.RawMechanism(
                cat.substrates, cat.products, cat.modifiers, steps, group),
            fill(:EqualAI, maximum(group)), 2, ER.Mechanisms.RawSite[]))
    end
    # The two bindings in one group: different metabolites.
    err = _testhelper_rebuild([bind_s, bind_p, iso], [1, 1, 2])
    @test err isa ErrorException
    @test occursin("a kinetic group holds E_P → EP (RE) and E_S → ES (RE)", err.msg)
    # The isomerization in two groups.
    err = _testhelper_rebuild([bind_s, iso, iso, bind_p], [1, 2, 3, 4])
    @test err isa ErrorException
    @test occursin("both hold the reaction ES ⇌ EP", err.msg)
    # The isomerization twice in one group.
    err = _testhelper_rebuild([bind_s, iso, iso, bind_p], [1, 2, 2, 3])
    @test err isa ErrorException
    @test occursin("holds the reaction ES ⇌ EP twice", err.msg)
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
    k = only(k for (k, e) in enumerate(mech.edges) if ER.Mechanisms._is_chemistry(mech, e))
    kon, koff = ER.Mechanisms._edge_constants(mech, k, :None)
    @test kon isa ER.Mechanisms.Kfor && koff isa ER.Mechanisms.Krev
    @test ER.Mechanisms.name(kon, mech) == :k_EQ_P_to_EAB &&
          ER.Mechanisms.name(koff, mech) == :k_EAB_to_EQ_P
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
    @test length(pp_names) == sum(pm.edges[first(ks)].rapid_equilibrium ? 1 : 2
                                  for ks in ER.Mechanisms._group_members(pm))
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
    k = only(k for (k, e) in enumerate(m.edges) if ER.Mechanisms._is_iso(e))
    @test ER.Mechanisms._form_name(m, m.edges[k].from) == "F" &&
          ER.Mechanisms._form_name(m, m.edges[k].to) == "E"
    kf, kr = ER.Mechanisms._edge_constants(m, k, :None)
    @test ER.Mechanisms.name(kf, m) == :k_F_to_E && ER.Mechanisms.name(kr, m) == :k_E_to_F
end

@testset "a competitive-inhibitor copy marks no reactant side" begin
    # P's competitive-inhibitor copy binds E. Tiers 1 and 2 read metabolite roles, so the
    # copy marks nothing at E, and the Iso Uni Uni isomerization still runs from F, where
    # P leaves, to E, where A enters.
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            E + P::Mod ⇌ E(P::Mod)
        end
    end)
    iso = only(e for e in m.edges
               if ER.Mechanisms._is_iso(e) &&
                   ER.Mechanisms._form_name(m, e.from) in ("E", "F"))
    @test ER.Mechanisms._form_name(m, iso.from) == "F" &&
          ER.Mechanisms._form_name(m, iso.to) == "E"
end

@testset "the tied steps of one kinetic group turn as one" begin
    # Segel Iso Uni Uni with an inhibitor I binding E and F in one group and the
    # isomerization mirrored at the I-bound forms. Tier 2 runs F → E; E(I) and F(I) mark
    # nothing, so on its own the mirror would fall to Tier 3 and run E(I) → F(I), and the
    # shared constant would tie F → E to the reverse of F(I) → E(I). Both run F-side to
    # E-side, so the cycle E → E(I) → F(I) → F → E constrains nothing and no line
    # forces k_F_to_E = k_E_to_F.
    ss = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            (F <--> E, F(I) <--> E(I))
            (E + I ⇌ E(I), F + I ⇌ F(I))
        end
    end
    m = ER.Mechanism(ss)
    label(f) = m.conformations[m.forms[f].conformation]
    mirror = only(m.edges[ks] for ks in ER.Mechanisms._group_members(m)
                  if length(ks) == 2 && all(k -> ER.Mechanisms._is_iso(m.edges[k]), ks))
    @test all(e -> label(e.from) == :F && label(e.to) == :E, mirror)
    @test !occursin("Wegscheider", rate_equation_string(ss))
    # The written direction and order of the mirrored steps do not matter.
    @test ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            (E(I) <--> F(I), E <--> F)
            (E + I ⇌ E(I), F + I ⇌ F(I))
        end
    end) == m

    # At rapid equilibrium the same cycle would pin K_F_to_E to 1 on its own; the shared
    # orientation leaves it fitted.
    re = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            (F ⇌ E, F(I) ⇌ E(I))
            (E + I ⇌ E(I), F + I ⇌ F(I))
        end
    end
    @test :K_F_to_E in parameters(re, FittedParameters())

    # Tied steps that change different pairs of conformations share no orientation.
    @test_throws "make different conformational changes" @enzyme_mechanism begin
        substrates: A
        products: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> E + P
            (E ⇌ F, F ⇌ G, G ⇌ E)
        end
    end
end

# The edges of the one kinetic group of `mech`'s catalytic graph that neither binds nor
# isomerizes, and what each takes up and gives off.
function _testhelper_exchanges(mech)
    cm = ER.Mechanisms.catalytic(mech)
    group = only(ks for ks in ER.Mechanisms._group_members(cm) if !any(
        k -> ER.Mechanisms._is_iso(cm.edges[k]) || ER.Mechanisms._is_binding(cm.edges[k]),
        ks))
    [(cm.table[ER.Mechanisms.entries(e.consumed)],
      cm.table[ER.Mechanisms.entries(e.released)]) for e in cm.edges[group]]
end

@testset "a kinetic group ties one metabolite exchange at two conformations" begin
    # Segel Iso Uni Uni with an inhibitor I binding E and F, and J displacing I at both
    # conformations in one group. Neither exchange changes conformation, so the group
    # stores both to take up J and give off I, however the steps are written.
    # The allosteric twin joins F to E at rapid equilibrium.
    # J binds E and F as well: without those bindings E(J) and F(J) are dead ends at zero
    # modifiers, a second closed class the guard rejects.
    I, J = ER.Mechanisms.Met(:I, true), ER.Mechanisms.Met(:J, true)
    plain = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E + J ⇌ E(J), F + J ⇌ F(J))
            (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)
        end
    end
    allo = @allosteric_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)                                   :: EqualAI
            E(A) <--> E(P)                                    :: EqualAI
            E(P) <--> F + P                                   :: EqualAI
            F ⇌ E                                             :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))                      :: EqualAI
            (E + J ⇌ E(J), F + J ⇌ F(J))                      :: EqualAI
            (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)  :: NonequalAI
        end
    end
    m, am = ER.Mechanism(plain), ER.AllostericMechanism(allo)
    for mech in (m, am)
        @test _testhelper_exchanges(mech) == [([J], [I]), ([J], [I])]
        _testhelper_check_stored_orientation(mech)
    end
    # The written direction and order of the exchanges do not matter.
    @test ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E + J ⇌ E(J), F + J ⇌ F(J))
            (F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)
        end
    end) == m
    @test ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)                                   :: EqualAI
            E(A) <--> E(P)                                    :: EqualAI
            E(P) <--> F + P                                   :: EqualAI
            F ⇌ E                                             :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))                      :: EqualAI
            (E + J ⇌ E(J), F + J ⇌ F(J))                      :: EqualAI
            (F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)  :: NonequalAI
        end
    end) == am
    # The rate law derives, with the exchange constants named in the J-in direction.
    @test occursin("k_EI_J_to_EJ_I", rate_equation_string(plain))
    allo_law = rate_equation_string(allo)
    @test occursin("k_A_EI_J_to_EJ_I", allo_law)
    @test occursin("k_I_EI_J_to_EJ_I", allo_law)
    _testhelper_check_law_against_oracle(plain)
    _testhelper_check_law_against_oracle(allo)

    # Exchanges that switch conformation in opposite directions as J displaces I share no
    # orientation.
    @test_throws "make different conformational changes" @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E(I) + J <--> F(J) + I, F(I) + J <--> E(J) + I)
        end
    end
    # Without the J bindings the guard rejects the mechanism at zero modifiers.
    dead_ends = @_testhelper_unguarded @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)
        end
    end
    v = EnzymeRates.Derivation._guard(dead_ends)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._NONUNIQUE, 4)
    @test _testhelper_respeed_verdicts(dead_ends; rng = Random.MersenneTwister(81),
                                       law = :limit)[4] === :nonunique
end

@testset "a tied exchange turns to take up its lead's metabolites" begin
    # The exchange fixture above plus A binding E(I) and P binding E(J). On its own, Tier 2
    # runs the exchange at E from E(J), where P binds, to E(I), where A binds, while Tier 3
    # runs the exchange at F from F(I) to F(J). The exchange at E leads, so the group
    # turns the one at F to take up I and give off J as well.
    # J binds E and F as well: without those bindings E(J) and F(J) are dead ends at zero
    # modifiers, a second closed class the guard rejects.
    I, J = ER.Mechanisms.Met(:I, true), ER.Mechanisms.Met(:J, true)
    plain = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E + J ⇌ E(J), F + J ⇌ F(J))
            E(I) + A ⇌ E(A, I)
            E(J) + P ⇌ E(J, P)
            (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)
        end
    end
    allo = @allosteric_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)                                   :: EqualAI
            E(A) <--> E(P)                                    :: EqualAI
            E(P) <--> F + P                                   :: EqualAI
            F ⇌ E                                             :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))                      :: EqualAI
            (E + J ⇌ E(J), F + J ⇌ F(J))                      :: EqualAI
            E(I) + A ⇌ E(A, I)                                :: EqualAI
            E(J) + P ⇌ E(J, P)                                :: EqualAI
            (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)  :: NonequalAI
        end
    end
    m, am = ER.Mechanism(plain), ER.AllostericMechanism(allo)
    for mech in (m, am)
        @test _testhelper_exchanges(mech) == [([I], [J]), ([I], [J])]
        _testhelper_check_stored_orientation(mech)
    end
    # The written direction and order of the exchanges do not matter.
    @test ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E + J ⇌ E(J), F + J ⇌ F(J))
            E(I) + A ⇌ E(A, I)
            E(J) + P ⇌ E(J, P)
            (F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)
        end
    end) == m
    @test ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)                                   :: EqualAI
            E(A) <--> E(P)                                    :: EqualAI
            E(P) <--> F + P                                   :: EqualAI
            F ⇌ E                                             :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))                      :: EqualAI
            (E + J ⇌ E(J), F + J ⇌ F(J))                      :: EqualAI
            E(I) + A ⇌ E(A, I)                                :: EqualAI
            E(J) + P ⇌ E(J, P)                                :: EqualAI
            (F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)  :: NonequalAI
        end
    end) == am
    @test occursin("k_EJ_I_to_EI_J", rate_equation_string(plain))
    allo_law = rate_equation_string(allo)
    @test occursin("k_A_EJ_I_to_EI_J", allo_law)
    @test occursin("k_I_EJ_I_to_EI_J", allo_law)
    _testhelper_check_law_against_oracle(plain)
    _testhelper_check_law_against_oracle(allo)

    # An exchange that changes conformation leads its group. On its own, Tier 3 runs
    # F(I) + J <--> E(J) + I from E(J) to F(I) and the exchange at E from E(I) to E(J);
    # the group turns the exchange at E to take up I and give off J as the lead does.
    mixed = @enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            E + J ⇌ E(J)
            (E(I) + J <--> E(J) + I, F(I) + J <--> E(J) + I)
        end
    end
    mm = ER.Mechanism(mixed)
    @test _testhelper_exchanges(mm) == [([I], [J]), ([I], [J])]
    _testhelper_check_stored_orientation(mm)
    @test ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        modifiers: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            E + J ⇌ E(J)
            (E(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)
        end
    end) == mm
    @test occursin("k_EJ_I_to_EI_J", rate_equation_string(mixed))
    _testhelper_check_law_against_oracle(mixed)
end

@testset "each reaction appears once in a mechanism" begin
    # The mechanism over A and P whose steps, each its own group, are `lines`, built and
    # checked as the constructor checks raw input.
    _testhelper_written(lines...) = ER.Mechanism(ER.Mechanisms.raw_mechanism(:(begin
        substrates: A
        products: P
        steps: $(Expr(:block, lines...))
    end), false))
    # A binds E into two different forms from two groups: two reactions, whose
    # names differ in the bound form (k_E_A_to_EA, k_E_A_to_EstarA).
    m = _testhelper_written(:(E + A <--> E(A)), :(E + A <--> Estar(A)),
                            :(E(A) <--> E(P)), :(E + P <--> E(P)))
    names = [ER.Mechanisms.name(p, m) for members in ER.Mechanisms._group_members(m)
             for p in ER.Mechanisms._edge_constants(m, first(members), :None)]
    @test allunique(names)
    @test issubset([:k_E_A_to_EA, :k_E_A_to_EstarA], names)
    # The same binding in two groups, once written as its release: both would be
    # k_E_A_to_EA.
    err = _testhelper_thrown() do
        _testhelper_written(:(E + A <--> E(A)), :(E(A) <--> E + A), :(E(A) <--> E(P)),
                            :(E + P <--> E(P)))
    end
    @test err isa ErrorException
    @test occursin("both hold the reaction E_A ⇌ EA", err.msg)
    # The same isomerization in two groups: both would be k_EA_to_EP.
    err = _testhelper_thrown() do
        _testhelper_written(:(E + A <--> E(A)), :(E(A) <--> E(P)), :(E(P) <--> E(A)),
                            :(E + P <--> E(P)))
    end
    @test err isa ErrorException
    @test occursin("both hold the reaction EA ⇌ EP", err.msg)
    # The same binding written RE in one group and SS in another: one reaction
    # cannot be both fast and slow.
    err = _testhelper_thrown() do
        @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                E + S <--> E(S)
                E + S ⇌ E(S)
                E(S) ⇌ E(P)
                E(P) ⇌ E + P
            end
        end
    end
    @test err isa ErrorException
    @test occursin("both hold the reaction E_S ⇌ ES (SS in group 2, RE in group 3)",
                   err.msg)
    # The same binding twice in one group, once written as its release: the
    # derivation would count its edge twice.
    err = _testhelper_thrown() do
        @enzyme_mechanism begin
            substrates: S
            products: P
            steps: begin
                (E + S <--> E(S), E(S) <--> E + S)
                E(S) <--> E(P)
                E(P) <--> E + P
            end
        end
    end
    @test err isa ErrorException
    @test occursin("kinetic group 2 holds the reaction E_S ⇌ ES twice", err.msg)
end

@testset "two distinct forms that render one name are rejected" begin
    # E with `star` bound and the conformation Estar with nothing bound both render
    # :Estar, so every constant named after one form would also name the other. The
    # reaction's names spell each name one way; the clash comes from a conformation label.
    err = _testhelper_thrown() do
        ER.Mechanism(ER.Mechanisms.raw_mechanism(:(begin
            substrates: star
            products: P
            steps: begin
                E + star ⇌ E(star)
                E(star) <--> E(P)
                E + P ⇌ E(P)
                E <--> Estar
            end
        end), false))
    end
    @test err isa ErrorException
    @test occursin("E(star)", err.msg)
    @test occursin(" and Estar both render", err.msg)
    @test occursin("both render the name Estar; ", err.msg)
    @test occursin("rename a metabolite", err.msg)
end

@testset "two forms whose constants would share a name are rejected" begin
    # P's modifier copy bound to E and the modifier Pmod bound to E would both name their
    # constants EPmod, so K_EPmod_to_E_Pmod would name both bindings. P::Mod and Pmod
    # spell one name, so the macro rejects them before it builds the mechanism.
    @test_throws ErrorException(
        "Mechanism: P::Mod and Pmod both spell Pmod: enzyme forms and their constants " *
        "join the names of the metabolites they bind without a separator, so two forms " *
        "would share a name; rename a metabolite") @enzyme_mechanism begin
        substrates: S
        products: P
        modifiers: P, Pmod
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P)
            E(P) ⇌ E + P
            E + P::Mod ⇌ E(P::Mod)
            E + Pmod ⇌ E(Pmod)
        end
    end
    # The allosteric macro checks its catalytic names alike: E(Ac, CoA) and E(AcCoA)
    # would both render EAcCoA.
    @test_throws ErrorException(
        "Mechanism: Ac + CoA and AcCoA both spell AcCoA: enzyme forms and their " *
        "constants join the names of the metabolites they bind without a separator, so " *
        "two forms would share a name; rename a metabolite") @allosteric_mechanism begin
        substrates: Ac, CoA
        products: AcCoA
        allosteric_regulators: X
        catalytic_multiplicity: 2
        steps: begin
            E + Ac ⇌ E(Ac)            :: EqualAI
            E(Ac) + CoA ⇌ E(Ac, CoA)  :: EqualAI
            E(Ac, CoA) <--> E(AcCoA)  :: OnlyA
            E(AcCoA) ⇌ E + AcCoA      :: EqualAI
        end
        regulatory_site(2): X::OnlyI
    end
end

@testset "E(I) and E(Iinh), E(P::Mod) and E(Pinh) name their constants apart" begin
    # E(I) and E(Iinh) render EI and EIinh: the mechanism is accepted, and each constant
    # takes a name of its own.
    iinh = @enzyme_mechanism begin
        substrates: Iinh
        products: P
        modifiers: I
        steps: begin
            E + Iinh ⇌ E(Iinh)
            E(Iinh) <--> E(P)
            E(P) ⇌ E + P
            E + I ⇌ E(I)
        end
    end
    @test _testhelper_step_constant_names(iinh) ==
          Set([:K_EIinh_to_E_Iinh, :k_EIinh_to_EP, :k_EP_to_EIinh, :K_EP_to_E_P,
               :K_EI_to_E_I])
    _testhelper_check_law_against_oracle(iinh)
    # So with a reactant's modifier copy: E(P::Mod) and E(Pinh) render EPmod and EPinh.
    pinh = @enzyme_mechanism begin
        substrates: Pinh
        products: P
        modifiers: P
        steps: begin
            E + Pinh ⇌ E(Pinh)
            E(Pinh) <--> E(P)
            E(P) ⇌ E + P
            E + P::Mod ⇌ E(P::Mod)
        end
    end
    @test _testhelper_step_constant_names(pinh) ==
          Set([:K_EPinh_to_E_Pinh, :k_EPinh_to_EP, :k_EP_to_EPinh, :K_EP_to_E_P,
               :K_EPmod_to_E_Pmod])
    _testhelper_check_law_against_oracle(pinh)
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
    return occursin(
        r"Parameter|::(Krapid|Kfor|Krev|Kreg)\b",
        arg_str)
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

@testset "chokepoint: no Symbol(\"[KkVL]...\") outside parameter-name renderers" begin
    src_dir = normpath(joinpath(@__DIR__, "..", "..", "src"))
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
    end
end

@testset "admission rejects a one-sided :OnlyA binding and admits balanced ones" begin
    # Admission's fault for `am`'s graph tagged `tags`, and the metabolite of the
    # binding it names (`nothing` when it names none).
    function _testhelper_admission_reads(am, tags)
        cm = am.catalytic
        v = ER.Derivation._admission(
            ER.AllostericMechanism(cm, tags, 1, ER.Mechanisms.RegulatorySite[]))
        e = v.edge == 0 ? nothing : cm.edges[v.edge]
        (v.fault, e === nothing || !ER.Mechanisms._is_binding(e) ? nothing :
                  cm.table[ER.Mechanisms._bound_entry(e)].name)
    end
    NO, NL = ER.Derivation._NO_FAULT, ER.Derivation._NO_LIMIT

    # Uni-uni S -> P. Tags: (S binding, chemical step, P binding). Admission's fault,
    # and the metabolite of the binding it names.
    function _testhelper_uni(s_tag, cat_tag, p_tag)
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_multiplicity: 2
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        am = ER.AllostericMechanism(m)
        firsts = _testhelper_first_steps(_testhelper_graph(am))
        # tag each group by the metabolite its step binds
        tags = [!_testhelper_is_binding(s) ? cat_tag :
                only(s.consumed) === :S ? s_tag : p_tag for s in firsts]
        _testhelper_admission_reads(am, tags)
    end

    # no :OnlyA anywhere -> admitted
    @test _testhelper_uni(:EqualAI, :EqualAI, :EqualAI) == (NO, nothing)
    # :OnlyA on the substrate only, catalysis :EqualAI: catalysis and the P binding
    # join E to E(S) in the inactive conformation -> no limit, by the S binding
    @test _testhelper_uni(:OnlyA, :EqualAI, :EqualAI) == (NL, :S)
    # :OnlyA on the product only, catalysis :EqualAI -> no limit, by the P binding
    @test _testhelper_uni(:EqualAI, :EqualAI, :OnlyA) == (NL, :P)
    # :OnlyA on the substrate and the catalysis: E(S) lies behind the binding
    @test _testhelper_uni(:OnlyA, :OnlyA, :EqualAI) == (NO, nothing)
    # balanced: :OnlyA on both bindings, catalysis :EqualAI -> admitted
    @test _testhelper_uni(:OnlyA, :EqualAI, :OnlyA) == (NO, nothing)
    # V-system: :OnlyA chemical step only -> admitted
    @test _testhelper_uni(:EqualAI, :OnlyA, :EqualAI) == (NO, nothing)
    # :NonequalAI catalysis joins E to E(S) too -> same verdict
    @test _testhelper_uni(:OnlyA, :NonequalAI, :EqualAI) == (NL, :S)

    # Whether a binding is rapid-equilibrium or steady-state does not change
    # admission's verdict: with both bindings :OnlyA, E(S) and E(P) lie behind them.
    # Both orientations of the mix are admitted.
    mixed_uni_re_ss = @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: EqualAI
            E(P) <--> E + P   :: EqualAI
        end
    end
    mixed_uni_ss_re = @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        steps: begin
            E + S <--> E(S)   :: EqualAI
            E(S) <--> E(P)    :: EqualAI
            E(P) ⇌ E + P      :: EqualAI
        end
    end
    # every binding :OnlyA, catalysis :EqualAI -> balanced -> valid
    function _testhelper_both_bindings_onlya(m)
        am = ER.AllostericMechanism(m)
        tags = [_testhelper_is_binding(s) ? :OnlyA : :EqualAI
                for s in _testhelper_first_steps(_testhelper_graph(am))]
        _testhelper_admission_reads(am, tags)
    end
    @test _testhelper_both_bindings_onlya(mixed_uni_re_ss) == (NO, nothing)
    @test _testhelper_both_bindings_onlya(mixed_uni_ss_re) == (NO, nothing)

    # A random-order binding square: E binds A and B, and each of E(A) and E(B) binds
    # the other. A lone :OnlyA edge has both sides joined by the other three, so
    # admission finds no limit; tagging P's binding too separates them.
    @testset "random-order binding square" begin
        biuni = @allosteric_mechanism begin
            substrates: A, B
            products:   P
            catalytic_multiplicity: 2
            steps: begin
                E + A ⇌ E(A)          :: EqualAI
                E + B ⇌ E(B)          :: EqualAI
                E(A) + B ⇌ E(A, B)    :: EqualAI
                E(B) + A ⇌ E(A, B)    :: EqualAI
                E(A, B) <--> E(P)     :: EqualAI
                E(P) ⇌ E + P          :: EqualAI
            end
        end
        bu_am = ER.AllostericMechanism(biuni)
        # Tag :OnlyA the groups named by (free form, bound metabolite); every
        # other group is :EqualAI. The chemical step's key is (:EAB, nothing).
        bu = _testhelper_graph(bu_am)
        # Each group's key: the form its first step starts from, and the metabolite it
        # binds (`nothing` for the chemical step).
        bu_keys = [(Symbol(bu.forms[s.from]),
                    _testhelper_is_binding(s) ? only(s.consumed) : nothing)
                   for s in _testhelper_first_steps(bu)]
        _testhelper_tags(onlya_keys...) =
            [key in onlya_keys ? :OnlyA : :EqualAI for key in bu_keys]
        _testhelper_verdict(t) = _testhelper_admission_reads(bu_am, t)

        # the square alone trips nothing: no :OnlyA anywhere -> admitted
        @test _testhelper_verdict(_testhelper_tags()) == (NO, nothing)
        # one square edge :OnlyA: the square's other three bindings join its sides
        @test _testhelper_verdict(_testhelper_tags((:E, :A))) == (NL, :A)
        # ...and the macro rejects that tagging: untagged steps join the two sides of the
        # :OnlyA binding of A, so its affinity has no limit (`_NO_LIMIT`).
        lone = @_testhelper_unguarded @allosteric_mechanism begin
            substrates: A, B
            products:   P
            catalytic_multiplicity: 2
            steps: begin
                E + A ⇌ E(A)          :: OnlyA
                E + B ⇌ E(B)          :: EqualAI
                E(A) + B ⇌ E(A, B)    :: EqualAI
                E(B) + A ⇌ E(A, B)    :: EqualAI
                E(A, B) <--> E(P)     :: EqualAI
                E(P) ⇌ E + P          :: EqualAI
            end
        end
        v = ER.Derivation._guard(lone)
        @test v.fault === ER.Derivation._NO_LIMIT
        @test ER.Mechanisms.step_text(lone.catalytic, v.edge) == "E + A ⇌ E(A)"
        @test_throws ErrorException(
            ER.Derivation._guard_message(lone, v)) @allosteric_mechanism begin
            substrates: A, B
            products:   P
            catalytic_multiplicity: 2
            steps: begin
                E + A ⇌ E(A)          :: OnlyA
                E + B ⇌ E(B)          :: EqualAI
                E(A) + B ⇌ E(A, B)    :: EqualAI
                E(B) + A ⇌ E(A, B)    :: EqualAI
                E(A, B) <--> E(P)     :: EqualAI
                E(P) ⇌ E + P          :: EqualAI
            end
        end
        # both A-side bindings :OnlyA: the B bindings, the catalysis and the P binding
        # still join E to E(A) -> no limit
        @test _testhelper_verdict(_testhelper_tags((:E, :A), (:EB, :A)))[1] === NL
        # ...and :OnlyA on P cuts the last join -> admitted
        @test _testhelper_verdict(
            _testhelper_tags((:E, :A), (:EB, :A), (:E, :P))) == (NO, nothing)
        # :OnlyA on the chemical step leaves the square, whose lone :OnlyA edge still
        # has its sides joined -> no limit, by that edge
        @test _testhelper_verdict(_testhelper_tags((:EAB, nothing), (:E, :A))) ==
              (NL, :A)
    end
end

@testset "admission rejects a ter-substrate cube joining a tagged binding's two sides" begin
    # A full random-order ter binding cube with seven :OnlyA edges and the cross bindings
    # :EqualAI. In the inactive conformation the untagged bindings join forms across the
    # cube, so a tagged binding's two sides share a component: no path takes its affinity to
    # zero. Admission finds untagged steps that join the two sides of an :OnlyA binding
    # (`_NO_LIMIT`).
    cube = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B, C
        products:   P
        catalytic_multiplicity: 2
        steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A) + C ⇌ E(A, C)       :: EqualAI
            E(B) + A ⇌ E(A, B)       :: EqualAI
            E(B) + C ⇌ E(B, C)       :: EqualAI
            E(C) + A ⇌ E(A, C)       :: EqualAI
            E(C) + B ⇌ E(B, C)       :: EqualAI
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    v = ER.Derivation._guard(cube)
    @test v.fault === ER.Derivation._NO_LIMIT
    @test_throws ErrorException(
        ER.Derivation._guard_message(cube, v)) @allosteric_mechanism begin
        substrates: A, B, C
        products:   P
        catalytic_multiplicity: 2
        steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A) + C ⇌ E(A, C)       :: EqualAI
            E(B) + A ⇌ E(A, B)       :: EqualAI
            E(B) + C ⇌ E(B, C)       :: EqualAI
            E(C) + A ⇌ E(A, C)       :: EqualAI
            E(C) + B ⇌ E(B, C)       :: EqualAI
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
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
    # and E(P), the bindings point away from free E and close no cycle, so admission
    # accepts the cube.
    cube = @allosteric_mechanism begin
        substrates: A, B, C
        products:   P
        catalytic_multiplicity: 2
        steps: begin
            E + A ⇌ E(A)             :: OnlyA
            E + B ⇌ E(B)             :: OnlyA
            E + C ⇌ E(C)             :: OnlyA
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A) + C ⇌ E(A, C)       :: OnlyA
            E(B) + A ⇌ E(A, B)       :: OnlyA
            E(B) + C ⇌ E(B, C)       :: OnlyA
            E(C) + A ⇌ E(A, C)       :: OnlyA
            E(C) + B ⇌ E(B, C)       :: OnlyA
            E(A, B) + C ⇌ E(A, B, C) :: OnlyA
            E(A, C) + B ⇌ E(A, B, C) :: OnlyA
            E(B, C) + A ⇌ E(A, B, C) :: OnlyA
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    @test cube isa ER.AllostericEnzymeMechanism

    # Ordered ter with an :OnlyA B binding and an :OnlyA chemical step: E(A, B) and
    # E(A, B, C) lie behind the B binding, on its bound side, so this is admitted.
    ordered = @allosteric_mechanism begin
        substrates: A, B, C
        products:   P
        catalytic_multiplicity: 2
        steps: begin
            E + A ⇌ E(A)             :: EqualAI
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A, B) + C ⇌ E(A, B, C) :: EqualAI
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    @test ordered isa ER.AllostericEnzymeMechanism
end

@testset "a ping-pong with both fused releases :OnlyA is an island; one half passes" begin
    # E(A) → E(; residual = A - P) + P runs the first half-reaction and releases P in
    # one step, E(B; residual = A - P) → E + Q the second and releases Q: chemistry
    # (`_is_chemistry`), though each gives off one metabolite and takes up none.
    # Admission reads each edge: in the inactive conformation the covalent forms join
    # free E only through the tagged chemistry, an island (`_ISLAND`) whose weight depends
    # on how fast the tagged constants vanish. Either half alone :OnlyA is admitted, and
    # both oracles agree on every face.
    both = @_testhelper_unguarded @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: OnlyA
        end
    end
    v = ER.Derivation._guard(both)
    @test v.fault === ER.Derivation._ISLAND
    @test _testhelper_epath_verdict(both; rng = Random.MersenneTwister(6)) in
          (:direction_dependent, :wrong_limit)
    @test_throws ErrorException(
        ER.Derivation._guard_message(both, v)) @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: OnlyA
        end
    end
    half1 = @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: EqualAI
        end
    end
    half2 = @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: EqualAI
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: OnlyA
        end
    end
    # Each face: the RE-speed oracle (class 3) and the ε-path (admission) agree with the
    # guard, which accepts it.
    for (i, half) in enumerate((half1, half2))
        @test all(_testhelper_admission_agreement(ER.AllostericMechanism(half);
                                                  rng = Random.MersenneTwister(6 + i)))
    end
end

@testset "admission admits an abortive segment an :OnlyA chemistry group strands" begin
    # The A group mixes the fused binding into E(P, Q) with A's abortive bindings at
    # E(P) and E(Q), so it is chemistry. Tagged :OnlyA, it strands
    # {E(A, P), E(A, Q), E(A, P, Q)} in the inactive conformation: a rapid-equilibrium
    # segment whose every form holds a product, which a mechanism may not have.
    am = @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        steps: begin
            E + B <--> E(B)                                              :: OnlyA
            (E(B) + A ⇌ E(P, Q), E(P) + A ⇌ E(A, P), E(Q) + A ⇌ E(A, Q)) :: OnlyA
            E(P, Q) <--> E(P) + Q                                        :: EqualAI
            E(P, Q) <--> E(Q) + P                                        :: EqualAI
            (E + P ⇌ E(P), E(A, Q) + P ⇌ E(A, P, Q))                     :: EqualAI
            (E + Q ⇌ E(Q), E(A, P) + Q ⇌ E(A, P, Q))                     :: EqualAI
        end
    end
    @test am isa ER.AllostericEnzymeMechanism
    # Admission reads each edge: the abortive bindings reach {E(A, P), E(A, Q),
    # E(A, P, Q)} from the product forms on their unbound side only, so the inactive
    # conformation holds no enzyme there and every face passes (the macro accepted it).
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
    m1 = EnzymeRates.Mechanism(@enzyme_mechanism begin
        substrates: A, B
        products: P, Q
        steps: begin
            E + A ⇌ E(A)
            E(A) + B ⇌ E(A, B)
            E(A, B) <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(Q) ⇌ E + Q
        end
    end)
    m2 = EnzymeRates.Mechanism(EnzymeRates.Mechanisms.RawMechanism(m1))
    am1 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: I
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)     :: EqualAI
            E(S) <--> E(P)   :: OnlyA
            E(P) ⇌ E + P     :: EqualAI
        end
        regulatory_site(2): I::OnlyI
    end)
    am2 = EnzymeRates.AllostericMechanism(
        EnzymeRates.Mechanisms.RawAllostericMechanism(am1))
    # Every parameter's name, which reads the mechanism and stores nothing in it.
    _testhelper_param_names(m::EnzymeRates.Mechanism) =
        EnzymeRates.Derivation._edge_names(m, :None)
    _testhelper_param_names(m::EnzymeRates.AllostericMechanism) =
        [[EnzymeRates.Derivation._edge_names(m, state) for state in (:A, :I)];
         [[K] for state in (:A, :I) for K in _testhelper_regulator_constants(m, state)]]
    # `rebuilt` is `m` constructed again from the raw input it writes back.
    function _testhelper_check_identity(m, rebuilt)
        shown = repr(m)
        names = _testhelper_param_names(m)
        @test repr(m) == shown
        @test m == rebuilt && hash(m) == hash(rebuilt)
        @test EnzymeRates.compile_mechanism(m) === EnzymeRates.compile_mechanism(rebuilt)
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
    m = EnzymeRates.Mechanism(_testhelper_re_mm)
    m_on_rxn = only(EnzymeRates.Enumeration.init_mechanisms(rxn))
    @test m_on_rxn == m && hash(m_on_rxn) == hash(m)
    am = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: I
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)     :: EqualAI
            E(S) <--> E(P)   :: OnlyA
            E(P) ⇌ E + P     :: EqualAI
        end
        regulatory_site(2): I::OnlyI
    end)
    am_on_rxn = EnzymeRates.AllostericMechanism(
        m_on_rxn, EnzymeRates.Mechanisms.cat_allo_states(am),
        EnzymeRates.Mechanisms.catalytic_multiplicity(am),
        EnzymeRates.Mechanisms.regulatory_sites(am))
    @test am_on_rxn == am && hash(am_on_rxn) == hash(am)
end

@testset "every derivation entry rejects a bottomless rapid-equilibrium segment" begin
    # Random product release at rapid equilibrium between steady-state releases leaves
    # the segment {E(P), E(P, Q), E(Q)} without a bottom form, so the rate is 0/0 at zero
    # products. The macro's boundary check (`check_valid_rate_equation`) rejects the
    # mechanism.
    # The constructor checks the representation only, so it builds the mechanism from
    # the walker's raw value, and the mechanism compiles; every derivation entry then
    # rejects it with the same message.
    raw = EnzymeRates.Mechanisms.raw_mechanism(:(begin
        substrates: S
        products: P, Q
        steps: begin
            E + S ⇌ E(S)
            E(S) <--> E(P, Q)
            E(P, Q) ⇌ E(Q) + P
            E(P, Q) ⇌ E(P) + Q
            E(P) <--> E + P
            E(Q) <--> E + Q
        end
    end), false)
    m = EnzymeRates.Mechanism(raw)
    v = EnzymeRates.Derivation._guard(m)
    @test (v.fault, v.face) == (EnzymeRates.Derivation._BOTTOMLESS, 2)
    msg = ErrorException(EnzymeRates.Derivation._guard_message(m, v))
    @test_throws msg EnzymeRates._compiled_mechanism(raw)
    em = EnzymeRates.compile_mechanism(m)
    @test_throws msg parameters(em)
    @test_throws msg parameters(em, FittedParameters())
    @test_throws msg parameters(em, AllParameters())
    @test_throws msg parameters(m, FittedParameters())
    @test_throws msg parameters(m)
    @test_throws msg rate_equation(em, (S = 1.0, P = 1.0, Q = 1.0),
                                   (Keq = 1.0, E_total = 1.0))
    @test_throws msg rate_equation_string(em)
    @test_throws msg parameters(em, (Keq = 1.0,))
    @test_throws msg parameters(m, (Keq = 1.0,))
    @test_throws msg EnzymeRates.Derivation._ss_rate_constant_names(em)
end

"""A raw-input value of each type: a metabolite, a form with a residual, a step, a
mechanism, a site and an allosteric mechanism."""
function _testhelper_raw_examples()
    R = EnzymeRates
    a, p, i = R.Mechanisms.Met(:A, false), R.Mechanisms.Met(:P, false),
        R.Mechanisms.Met(:I, true)
    free = R.Mechanisms.RawForm(:E, R.Mechanisms.Met[])
    bind = R.Mechanisms.RawStep(free, R.Mechanisms.RawForm(:E, [a]), [a],
                                R.Mechanisms.Met[], true)
    chem =
        R.Mechanisms.RawStep(R.Mechanisms.RawForm(:E, [a]), R.Mechanisms.RawForm(:E, [p]),
                             R.Mechanisms.Met[], R.Mechanisms.Met[], false)
    release = R.Mechanisms.RawStep(
        R.Mechanisms.RawForm(:E, [p]), free, R.Mechanisms.Met[], [p], true)
    dead_end = R.Mechanisms.RawStep(free, R.Mechanisms.RawForm(:E, [i]), [i],
                                    R.Mechanisms.Met[], true)
    mech = R.Mechanisms.RawMechanism([:A], [:P], [:I], [bind, chem, release, dead_end],
                                     [1, 2, 3, 4])
    site = R.Mechanisms.RawSite([:X], 2, [:OnlyI])
    allo = R.Mechanisms.RawAllostericMechanism(
        mech, [:EqualAI, :OnlyA, :EqualAI, :EqualAI], 2, [site])
    (i, R.Mechanisms.RawForm(:Estar, [a, i], [:A], [:P]), bind, mech, site, allo)
end

@testset "raw input: sorted lists, content equality and hash" begin
    R = EnzymeRates
    a, b, i = R.Mechanisms.Met(:A, false), R.Mechanisms.Met(:B, false),
        R.Mechanisms.Met(:I, true)
    # A metabolite orders by name, a modifier after the reactant of its name.
    @test sort([R.Mechanisms.Met(:P, true), R.Mechanisms.Met(:P, false),
                R.Mechanisms.Met(:A, true)]) == [
        R.Mechanisms.Met(:A, true), R.Mechanisms.Met(:P, false), R.Mechanisms.Met(:P, true)]
    @test R.Mechanisms.Met(:A, true) != a && hash(R.Mechanisms.Met(:A, true)) != hash(a)
    # A form sorts its bound metabolites and its residual's two lists, so two writings
    # of one form are equal and hash alike.
    form = R.Mechanisms.RawForm(:E, [i, b, a], [:B, :A], [:Q, :P])
    @test form.bound == [a, b, i] && form.added == [:A, :B] && form.subtracted == [:P, :Q]
    same = R.Mechanisms.RawForm(:E, [a, i, b], [:A, :B], [:P, :Q])
    @test form == same && hash(form) == hash(same)
    @test R.Mechanisms.RawForm(:E, [a]) == R.Mechanisms.RawForm(:E, [a], Symbol[], Symbol[])
    @test R.Mechanisms.RawForm(:Estar, [a]) != R.Mechanisms.RawForm(:E, [a])
    # A step sorts its lists and keeps its written direction: a release stays a release.
    free, ea, eab = R.Mechanisms.RawForm(:E, R.Mechanisms.Met[]),
        R.Mechanisms.RawForm(:E, [a]), R.Mechanisms.RawForm(:E, [a, b])
    release = R.Mechanisms.RawStep(ea, free, R.Mechanisms.Met[], [a], true)
    @test (release.from, release.to) == (ea, free)
    @test isempty(release.consumed) && release.released == [a]
    @test release != R.Mechanisms.RawStep(free, ea, [a], R.Mechanisms.Met[], true)
    both = R.Mechanisms.RawStep(free, eab, [b, a], R.Mechanisms.Met[], false)
    @test both.consumed == [a, b]
    @test both == R.Mechanisms.RawStep(free, eab, [a, b], R.Mechanisms.Met[], false)
    @test hash(both) ==
          hash(R.Mechanisms.RawStep(free, eab, [a, b], R.Mechanisms.Met[], false))
    @test both != R.Mechanisms.RawStep(free, eab, [a, b], R.Mechanisms.Met[], true)
    # Mechanisms, sites and allosteric mechanisms built twice are equal and hash alike;
    # each field takes part.
    examples, again = _testhelper_raw_examples(), _testhelper_raw_examples()
    @test all(examples .== again) && hash.(examples) == hash.(again)
    _, _, _, mech, site, allo = examples
    @test mech != R.Mechanisms.RawMechanism(mech.substrates, mech.products, Symbol[],
                                            mech.steps, mech.group)
    @test mech != R.Mechanisms.RawMechanism(mech.substrates, mech.products, mech.modifiers,
                                            mech.steps, [1, 1, 2, 3])
    @test site != R.Mechanisms.RawSite([:X], 4, [:OnlyI])
    @test site != R.Mechanisms.RawSite([:X], 2, [:OnlyA])
    @test allo != R.Mechanisms.RawAllostericMechanism(
        mech, [:EqualAI, :NonequalAI, :EqualAI, :EqualAI], 2, [site])
    @test allo != R.Mechanisms.RawAllostericMechanism(mech, allo.cat_allo_states, 1, [site])
    @test allo != R.Mechanisms.RawAllostericMechanism(mech, allo.cat_allo_states, 2,
                                                      R.Mechanisms.RawSite[])
end

@testset "a raw value's hash reads its type's name and its fields only" begin
    # Base hashes names, flags, numbers and vectors alike in every process, but a
    # struct's default hash mixes in its type's identity, which changes between
    # precompile sessions; so each raw type hashes its name, then its fields in order.
    h = UInt(17)
    for x in _testhelper_raw_examples()
        @test hash(x, h) == foldl((acc, f) -> hash(getfield(x, f), acc),
                                  fieldnames(typeof(x)); init = hash(nameof(typeof(x)), h))
    end
end

@testset "a regulatory-site ligand takes no reserved name and no `_`" begin
    # An allosteric regulator reaches the constructor only as a site ligand, so the
    # AllostericMechanism constructor checks the ligands' names.
    @test_throws(ErrorException(
        "AllostericMechanism: L is a reserved name; E_total, Keq, L, N, params, concs " *
        "and v are kept for a rate equation's constants and local names"),
        @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: L
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: OnlyA
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): L::OnlyI
        end)
    @test_throws(ErrorException(
        "AllostericMechanism: the metabolite name A_R holds `_`; constant names join " *
        "metabolite names with `_`, so a metabolite name holds none"),
        @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: A_R
            steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: OnlyA
                E(P) ⇌ E + P      :: EqualAI
            end
            regulatory_site(2): A_R::EqualAI
        end)
end

@testset "a reaction's metabolite names spell each name one way" begin
    # Enzyme forms and their constants join the names of the metabolites they bind, so a
    # reaction is rejected when its names spell one string two ways; a reactant also
    # declared a competitive inhibitor binds as its copy, spelled X::Mod, named Xmod.
    _testhelper_two_ways(a, b, x) = ErrorException(
        "EnzymeReaction: $a and $b both spell $x: enzyme forms and their constants join " *
        "the names of the metabolites they bind without a separator, so two forms would " *
        "share a name; rename a metabolite")
    @test_throws _testhelper_two_ways("Ac + CoA", "AcCoA", "AcCoA") @enzyme_reaction begin
        substrates: Ac[C], CoA[S]
        products:   AcCoA[CS]
    end
    @test_throws _testhelper_two_ways("A + I", "AI", "AI") @enzyme_reaction begin
        substrates: A[C], AI[N]
        products:   P[CN]
        competitive_inhibitors: I
    end
    # The split can take several steps to find: XY + Z and X + YZ both spell XYZ.
    @test_throws _testhelper_two_ways("XY + Z", "X + YZ", "XYZ") @enzyme_reaction begin
        substrates: X[C], YZ[N]
        products:   XY[C], Z[N]
    end
    # A reactant's modifier copy is named Xmod.
    @test_throws _testhelper_two_ways("P::Mod", "Pmod", "Pmod") @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        competitive_inhibitors: P, Pmod
    end
    @test_throws _testhelper_two_ways("X + mod", "X::Mod", "Xmod") @enzyme_reaction begin
        substrates: X[C]
        products:   mod[C]
        competitive_inhibitors: X
    end
    # NAD and NADH spell no name two ways; with H, NAD + H spells NADH.
    nadh = @enzyme_reaction begin
        substrates: NAD[C]
        products:   NADH[C]
    end
    @test [p.name for p in nadh.products] == [:NADH]
    @test_throws _testhelper_two_ways("NAD + H", "NADH", "NADH") @enzyme_reaction begin
        substrates: NAD[C], H[N]
        products:   NADH[CN]
    end
    # An allosteric regulator binds no enzyme form, so its name takes no part.
    allosteric = @enzyme_reaction begin
        substrates: A[C]
        products:   P[C]
        competitive_inhibitors: I
        allosteric_regulators: AI
    end
    @test [reg.name for reg in allosteric.regulators] == [:AI, :I]
end
