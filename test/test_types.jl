# ABOUTME: Tests for the core EnzymeRates types (Mechanism, Step, Species,
# ABOUTME: Parameter family) and their constructors, accessors, and equality.
using Serialization

const ER = EnzymeRates
_testhelper_sp(bound, conf = :E) = ER.Species(ER.Metabolite[bound...], conf)
_testhelper_sp(bound, conf, res) = ER.Species(ER.Metabolite[bound...], conf, res)

# The AllostericEnzymeMechanism over catalytic mechanism `cm`, its catalytic allosteric
# states given in `cm`'s canonical group order (`_testhelper_allo_from_source` with the
# canonical groups as the source).
_testhelper_aem(cm, cat_sites, reg_sites) =
    _testhelper_allo_from_source((cm, ER.steps(ER.Mechanism(cm))), cat_sites, reg_sites)

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

# The reaction S ⇌ P, its three forms, and the steps that bind S, isomerize ES to EP
# and bind P to E.
function _testhelper_uniuni(; mults = [1], s = :S, p = :P)
    rxn = ER.EnzymeReaction(
        [ER.ReactantAtoms(ER.Substrate(s), [:C => 1]),
         ER.ReactantAtoms(ER.Product(p), [:C => 1])],
        ER.RegulatorMults[], mults)
    E, ES, EP = _testhelper_sp([]), _testhelper_sp([ER.Substrate(s)]),
                _testhelper_sp([ER.Product(p)])
    bind = ER.Step(E, ES, [ER.Substrate(s)], ER.Metabolite[], true)
    iso = ER.Step(ES, EP, ER.Metabolite[], ER.Metabolite[], false)
    rel = ER.Step(E, EP, [ER.Product(p)], ER.Metabolite[], true)
    (; rxn, E, ES, EP, bind, iso, rel)
end

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
        @test Set(ER._step_text.(_testhelper_flat_steps(m))) ==
              Set(["E + S ⇌ ES", "ES <--> EP", "E + P ⇌ EP"])
        @test length(ER.steps(ER.Mechanism(m))) == 3
        @test Set(_testhelper_enzyme_forms(_testhelper_flat_steps(m))) ==
              Set([:E, :ES, :EP])

        # Shared kinetic-group: two steps in group 1 (regulator R binds
        # both E and E(S) sharing one K).
        m2 = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: R
            steps: begin
                (E + R ⇌ E(R), E(S) + R ⇌ E(S, R))
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        # Steps are canonicalized; the two R-binding steps share one kinetic
        # group → 5 steps, 4 groups (one 2-step group binding R).
        @test length(ER.steps(ER.Mechanism(m2))) == 4
        shared = only(g for g in ER.Mechanism(m2).steps if length(g) == 2)
        @test all(ER.bound_metabolite(s) ==
                  ER.CompetitiveInhibitor(:R) for s in shared)
    end

    @testset "metabolites() lists substrates, then products, then regulators" begin
        # A substrate, a product and a regulator, so the lift lists all three roles in
        # order — coverage the plain S/P accessor tests do not reach.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                E + S ⇌ E(S)
                E(S) + I ⇌ E(S, I)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test ER.metabolites(m) == (:S, :P, :I)
    end

    @testset "AllostericEnzymeMechanism lift + DSL" begin
        cm = _testhelper_re_mm

        # Single-ligand :EqualAI reg site is allowed (degenerate but valid —
        # the enumerator won't emit it, but users may write it for teaching).
        @test _testhelper_aem(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI)),
            (((:I,), 2, (:EqualAI,)),),
        ) isa ER.AllostericEnzymeMechanism

        # Catalytic group :OnlyI → error
        @test_throws ErrorException _testhelper_aem(
            cm, (2, (:NonequalAI, :OnlyI, :NonequalAI)), (),
        )

        # Build via DSL
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I::OnlyI

            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)     :: EqualAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
        end
        @test ER.catalytic_multiplicity(m) == 2
        # Catalysis (iso step) is :OnlyA, both bindings :EqualAI. Group order is
        # canonical, so identify the tagged group by its step, not its position.
        am_c = ER.AllostericMechanism(m)
        onlyA_g = only(g for g in eachindex(ER.steps(am_c))
                       if ER.cat_allo_state(am_c, g) === :OnlyA)
        @test ER.bound_metabolite(first(ER.steps(am_c)[onlyA_g])) === nothing
        @test all(ER.cat_allo_state(am_c, g) === :EqualAI
                  for g in eachindex(ER.steps(am_c)) if g != onlyA_g)
        @test ER.allosteric_regulators(am_c) == [ER.AllostericRegulator(:I)]
        @test ER.allo_states(only(ER.regulatory_sites(am_c))) == [:OnlyI]
    end

    @testset "Pretty printing" begin
        # Linear mechanism: one line per step, in stored order, under a header that
        # counts steps and enzyme forms.
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
            end
        end
        @test sprint(show, m) ==
            "EnzymeMechanism (2 steps, 2 enzyme forms):\n  E + P <--> ES\n" *
            "  E + S <--> ES"

        # A mechanism with no steps prints the bare header.
        m_none = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
            end
        end
        @test sprint(show, m_none) == "EnzymeMechanism (0 steps, 0 enzyme forms):"

        # Branched mechanism: the header counts every step and enzyme form, and each
        # step prints.
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
        s = sprint(show, m_b)
        @test startswith(s, "EnzymeMechanism (7 steps, 6 enzyme forms):")
        @test count(==('\n'), s) == 7
        @test contains(s, "E + A <--> EA")
        @test contains(s, "E + Q <--> EQ")

        # Catalytic 3-cycle: the steps print in stored order, not chain order.
        m_re = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                E + S ⇌ E(S)
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        @test sprint(show, m_re) ==
            "EnzymeMechanism (3 steps, 3 enzyme forms):\n  E + P ⇌ EP\n" *
            "  E + S ⇌ ES\n  ES <--> EP"

        # Theorell–Chance: B binds and P leaves in one step, and the step prints
        # both on its own line.
        m_tc = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                E + A <--> E(A)
                E(A) + B <--> E(Q) + P
                E(Q) <--> E + Q
            end
        end
        @test sprint(show, m_tc) ==
            "EnzymeMechanism (3 steps, 3 enzyme forms):\n  E + A <--> EA\n" *
            "  E + Q <--> EQ\n  EA + B <--> EQ + P"

        # A fused binding (B binds and chemistry runs in one step) prints B on its
        # entry side, since E(P, Q) does not name it.
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
        @test contains(sprint(show, m_fb), "EA + B <--> EPQ")

        # Mechanism with regulators: appended at end.
        m_reg = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                E + S <--> E(S)
                E(S) <--> E + P
                E + I <--> E(I)
            end
        end
        @test contains(sprint(show, m_reg), "| regulators: I")

        # EnzymeReaction with oligomeric_state > 1.
        rxn_oligo = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            oligomeric_state: 4
        end
        @test sprint(show, rxn_oligo) ==
            "EnzymeReaction: S ⇌ P | oligomeric_state: 4"

        # EnzymeReaction with shared_catalytic_site: surfaced in show.
        rxn_shared = @enzyme_reaction begin
            substrates: A[C], B[C]
            products:   P[C], Q[C]
            shared_catalytic_site: (A, P)
        end
        @test sprint(show, rxn_shared) ==
            "EnzymeReaction: A + B ⇌ P + Q | shared_catalytic_site: (A, P)"

        # AllostericEnzymeMechanism (smoke).
        m_allo = @allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            allosteric_regulators: I::OnlyI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + F6P ⇌ E(F6P)         :: EqualAI
                E(F6P) <--> E(F16BP)     :: EqualAI
                E(F16BP) ⇌ E + F16BP     :: EqualAI
            end
        end
        s_allo = sprint(show, m_allo)
        @test contains(s_allo, "AllostericEnzymeMechanism (cat_n=2")
        @test contains(s_allo, "reg sites")
        @test contains(s_allo, "I::OnlyI")

        # A regulatory site with two ligands lists both with their states.
        m_site = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I::OnlyI, J::EqualAI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)     :: EqualAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
            regulatory_site(multiplicity = 3): begin
                ligands: I, J
            end
        end
        @test sprint(show, m_site) ==
            "AllostericEnzymeMechanism (cat_n=2, 1 reg sites):\n" *
            "  E + P ⇌ EP :: EqualAI\n  E + S ⇌ ES :: EqualAI\n" *
            "  ES <--> EP :: OnlyA\n  reg site 1 (n=3): I, J [I::OnlyI, J::EqualAI]"
    end

    @testset "Mechanism rejects a rapid-equilibrium segment with no bottom form" begin
        # Random product release from E(P, Q) at rapid equilibrium, but release
        # from E(P) and E(Q) at steady state: the RE segment {E(P), E(Q), E(P, Q)}
        # has weights P·Q : Kq·P : Kp·Q, all zero at P = Q = 0. The RE
        # approximation carries no parameter for how E(P, Q) splits between
        # E(P) and E(Q) there, so the rate is undefined at zero products.
        err = try
            @enzyme_mechanism begin
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
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        msg = sprint(showerror, err)
        @test occursin("rapid-equilibrium segment", msg)
        @test occursin("E(P, Q)", msg) || occursin("EPQ", msg)

        # Mirror on the substrate side: undefined at A = B = 0.
        @test_throws ErrorException @enzyme_mechanism begin
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

        # The same segment with A's competitive-inhibitor copy binding E(B) at rapid
        # equilibrium: E(A::Inh, B) joins it, but every form still carries A or B, so
        # the segment stays bottomless. The copy shares A's name and must not hide
        # A's substrate role, whatever the order of the groups.
        A, B, Ainh = ER.Substrate(:A), ER.Substrate(:B), ER.CompetitiveInhibitor(:A)
        inh_steps = [
            [ER.Step(_testhelper_sp([]), _testhelper_sp([A]), [A], ER.Metabolite[],
                     false)],
            [ER.Step(_testhelper_sp([]), _testhelper_sp([B]), [B], ER.Metabolite[],
                     false)],
            [ER.Step(_testhelper_sp([A]), _testhelper_sp([A, B]), [B], ER.Metabolite[],
                     true)],
            [ER.Step(_testhelper_sp([B]), _testhelper_sp([A, B]), [A], ER.Metabolite[],
                     true)],
            [ER.Step(_testhelper_sp([B]), _testhelper_sp([Ainh, B]), [Ainh],
                     ER.Metabolite[], true)],
        ]
        rxn = @enzyme_reaction(begin
            substrates: A[C], B[C]
            products: P[C2]
        end)
        @test ER._bottomless_re_segment(rxn, inh_steps) !== nothing
        @test ER._bottomless_re_segment(rxn, reverse(inh_steps)) !== nothing
        @test_throws "rapid-equilibrium segment" @enzyme_mechanism begin
            substrates: A, B
            products: P
            regulators: A
            steps: begin
                E + A <--> E(A)
                E + B <--> E(B)
                E(A) + B ⇌ E(A, B)
                E(B) + A ⇌ E(A, B)
                E(B) + A::Inh ⇌ E(A::Inh, B)
                E + P ⇌ E(P)
                E(A, B) <--> E(P)
            end
        end

        # Ordered release: the segment {E(P), E(P, Q)} has weights Kq : Q, so
        # E(P) is its bottom form and the rate is defined at zero products.
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

        # Mixed abortive complex: the segment {E(S), E(P), E(S, P)} has weights
        # K·S : K′·P : S·P, zero only at S = P = 0 where no turnover is
        # possible anyway. Accepted.
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

        # NOTE: unreachable enzyme forms are accepted. The constructor does not
        # enforce a connectivity invariant — enzyme forms are inferred from steps,
        # so an "unreachable" form simply has its own steps in isolation, which is
        # structurally valid (graph connectivity is a downstream concern caught by
        # Wegscheider analysis if it matters).
    end

    @testset "an inhibitor copy leaves the inactive state's segment bottomless" begin
        # Every step that takes up the substrate A is :OnlyA, so the inactive conformation
        # keeps none of them; the second free conformation F still reaches E(B), and from
        # it the copy's bindings. The segment {E(A::Inh), E(A::Inh, B), E(B)} then has
        # every form carrying A (the copy's concentration is A's) or B, so it is empty at
        # A = B = 0. In the active conformation E + A::Inh ⇌ E(A::Inh) joins E to the
        # segment as its bottom form.
        am = ER.AllostericMechanism(@allosteric_mechanism begin
            substrates: A, B
            products: P
            catalytic_inhibitors: A
            catalytic_steps: begin
                E + A <--> E(A)                   :: OnlyA
                E(A) + B <--> E(A, B)             :: EqualAI
                E(A, B) <--> E(P)                 :: OnlyA
                E(P) <--> E + P                   :: EqualAI
                E + A::Inh ⇌ E(A::Inh)            :: OnlyA
                E(A::Inh) + B ⇌ E(A::Inh, B)      :: EqualAI
                E(B) + A::Inh ⇌ E(A::Inh, B)      :: EqualAI
                F + B <--> E(B)                   :: EqualAI
            end
        end)
        @test_throws "rapid-equilibrium segment" ER._state_allo_mechanism(am, :I)
    end

    @testset "AllostericEnzymeMechanism lift validators" begin
        cm = _testhelper_re_mm
        # Wrong-length cat_allo_states (4 entries for 3 kinetic groups) → error
        @test_throws ErrorException _testhelper_aem(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI, :OnlyA)), ())

        # Invalid allo state value → error
        @test_throws ErrorException _testhelper_aem(
            cm, (2, (:NotAState, :NonequalAI, :NonequalAI)), ())

        # Reg site with all-:EqualAI ligands is allowed (degenerate but valid).
        @test _testhelper_aem(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI)),
            (((:I, :J), 2, (:EqualAI, :EqualAI)),)) isa
              ER.AllostericEnzymeMechanism

        # Invalid reg-site ligand allo state → error
        @test_throws ErrorException _testhelper_aem(
            cm, (2, (:NonequalAI, :NonequalAI, :NonequalAI)),
            (((:I,), 2, (:NotAState,)),))
    end

    @testset "Base.show displays all dense states" begin
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            allosteric_regulators: I::NonequalAI, J::OnlyI
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)     :: NonequalAI
                E(S) <--> E(P)   :: OnlyA
                E(P) ⇌ E + P     :: EqualAI
            end
        end
        s = sprint(show, m)
        @test s == "AllostericEnzymeMechanism (cat_n=2, 2 reg sites):\n" *
                   "  E + P ⇌ EP :: EqualAI\n  E + S ⇌ ES :: NonequalAI\n" *
                   "  ES <--> EP :: OnlyA\n  reg site 1 (n=2): I [I::NonequalAI]\n" *
                   "  reg site 2 (n=2): J [J::OnlyI]"
        am = ER.AllostericMechanism(m)
        # Every catalytic state appears in cat_allo_states line
        for g in eachindex(ER.steps(am))
            @test occursin(string(ER.cat_allo_state(am, g)), s)
        end
        # No :NonequalAI ligand silently hidden from reg-site display
        for site in ER.regulatory_sites(am),
            (lig, state) in zip(ER.ligands(site), ER.allo_states(site))
            @test occursin("$(ER.name(lig))::$state", s)
        end
    end

    @testset "AllostericEnzymeMechanism display: shared kinetic group" begin
        cm = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                (E + S ⇌ E(S), E(P) + S ⇌ E(P, S))
                E(S) <--> E(P)
                E(P) ⇌ E + P
            end
        end
        am = _testhelper_aem(
            cm, (2, (:EqualAI, :EqualAI, :EqualAI)), ())
        s = repr(am)

        # The parenthesized-group branch must execute: look for the
        # exact "(...) :: EqualAI" shape with a comma-separated body.
        paren_group_match = match(
            r"\([^()]*,[^()]*\) :: EqualAI", s)
        @test paren_group_match !== nothing
        # Format: lhs/rhs joined with " + " (no brackets).
        @test occursin("ES <--> EP :: EqualAI", s)
        @test occursin(":: EqualAI", s)
        @test !occursin("cat_allo_states:", s)
        # One line per kinetic group: a multi-line catalytic display.
        @test count(==('\n'), s) >= 3
        @test s == "AllostericEnzymeMechanism (cat_n=2):\n  E + P ⇌ EP :: EqualAI\n" *
                   "  (E + S ⇌ ES, EP + S ⇌ EPS) :: EqualAI\n  ES <--> EP :: EqualAI"
    end

    # ─── Concrete type hierarchy ──────────────────────────────────────

    @testset "Metabolite hierarchy: Substrate / Product / Regulators" begin
        s = ER.Substrate(:ATP)
        p = ER.Product(:ADP)
        a = ER.AllostericRegulator(:cAMP)
        i = ER.CompetitiveInhibitor(:I)

        @test s isa ER.Reactant
        @test s isa ER.Metabolite
        @test p isa ER.Reactant
        @test a isa ER.Regulator
        @test a isa ER.Metabolite
        @test i isa ER.Regulator

        @test ER.name(s) === :ATP
        @test ER.name(p) === :ADP
        @test ER.name(a) === :cAMP
        @test ER.name(i) === :I

        @test ER.Substrate(:X) == ER.Substrate(:X)
        @test ER.Substrate(:X) != ER.Substrate(:Y)
        # Distinct subtypes with same name are NOT equal (struct identity matters).
        @test ER.Substrate(:X) != ER.Product(:X)
        @test hash(ER.Substrate(:X)) == hash(ER.Substrate(:X))
    end

    @testset "Residual: empty default + canonical ordering" begin
        empty_r = ER.Residual()
        @test isempty(empty_r)
        @test ER.added(empty_r) == ER.Substrate[]
        @test ER.subtracted(empty_r) == ER.Product[]

        r1 = ER.Residual(
            [ER.Substrate(:B), ER.Substrate(:A)],
            [ER.Product(:Q), ER.Product(:P)],
        )
        r2 = ER.Residual(
            [ER.Substrate(:A), ER.Substrate(:B)],
            [ER.Product(:P), ER.Product(:Q)],
        )
        @test r1 == r2
        @test hash(r1) == hash(r2)
        @test !isempty(r1)
        @test ER.added(r1) ==
              [ER.Substrate(:A), ER.Substrate(:B)]
        @test ER.subtracted(r1) ==
              [ER.Product(:P), ER.Product(:Q)]
    end

    @testset "Species: canonical bound ordering + accessors + name" begin
        s1 = ER.Species(
            ER.Metabolite[
                ER.Substrate(:B), ER.Substrate(:A)],
            :E,
        )
        s2 = ER.Species(
            ER.Metabolite[
                ER.Substrate(:A), ER.Substrate(:B)],
            :E,
        )
        @test s1 == s2
        @test hash(s1) == hash(s2)
        @test ER.conformation(s1) === :E
        @test ER.residual(s1) == ER.Residual()
        @test !ER.has_residual(s1)
        @test ER.bound(s1) ==
              ER.Metabolite[
                  ER.Substrate(:A), ER.Substrate(:B)]

        # Empty bound, :E conformation → :E
        s_e = ER.Species(ER.Metabolite[], :E)
        @test ER.name(s_e) === :E

        # Bound metabolites are appended in canonical order
        s_es = ER.Species(
            ER.Metabolite[ER.Substrate(:A),
                                   ER.Substrate(:B)],
            :E,
        )
        @test ER.name(s_es) === :EAB

        # Estar conformation
        s_estar = ER.Species(ER.Metabolite[], :Estar)
        @test ER.name(s_estar) === :Estar

        # Three-arg constructor exposes residual
        res = ER.Residual(
            [ER.Substrate(:A)],
            [ER.Product(:P)],
        )
        s_res = ER.Species(ER.Metabolite[], :Estar, res)
        @test ER.has_residual(s_res)
        @test ER.residual(s_res) == res

        # Same metabolite name bound in two roles (product + competitive
        # inhibitor) canonicalizes regardless of construction order: same
        # Species, same hash, same form name (non-inhibitor segment first).
        d1 = ER.Species(
            ER.Metabolite[ER.Product(:G6P),
                                   ER.CompetitiveInhibitor(:G6P)], :E)
        d2 = ER.Species(
            ER.Metabolite[ER.CompetitiveInhibitor(:G6P),
                                   ER.Product(:G6P)], :E)
        @test d1 == d2
        @test hash(d1) == hash(d2)
        @test ER.name(d1) === ER.name(d2)
        @test ER.name(d1) === :EG6PG6Pinh
    end

    @testset "RegulatorySite: validation + accessors" begin
        lig_a = ER.AllostericRegulator(:A)
        lig_b = ER.AllostericRegulator(:B)
        site = ER.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA, :NonequalAI],
        )
        @test ER.ligands(site) == [lig_a, lig_b]
        @test ER.multiplicity(site) == 4
        @test ER.allo_states(site) == [:OnlyA, :NonequalAI]

        # Mismatched ligand / allo_state length → error
        @test_throws ErrorException ER.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA])

        # Multiplicity < 1 → error
        @test_throws ErrorException ER.RegulatorySite(
            [lig_a], 0, [:OnlyA])

        # No ligands → error
        @test_throws ErrorException ER.RegulatorySite(
            ER.AllostericRegulator[], 2, Symbol[])

        # Invalid allo state → error
        @test_throws ErrorException ER.RegulatorySite(
            [lig_a], 1, [:NotAState])

        # All four allowed states accepted; the R/T symbols are rejected.
        for st in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
            @test ER.RegulatorySite([lig_a], 1, [st]) isa
                  ER.RegulatorySite
        end
        for st in (:OnlyR, :OnlyT, :EqualRT, :NonequalRT)
            @test_throws ErrorException ER.RegulatorySite([lig_a], 1, [st])
        end

        # Equality / hash
        site2 = ER.RegulatorySite(
            [lig_a, lig_b], 4, [:OnlyA, :NonequalAI])
        @test site == site2
        @test hash(site) == hash(site2)
        # The constructor sorts the ligands by name and carries each ligand's state
        # along, so the order the ligands are given in does not matter...
        site_reordered = ER.RegulatorySite(
            [lig_b, lig_a], 4, [:NonequalAI, :OnlyA])
        @test ER.ligands(site_reordered) == [lig_a, lig_b]
        @test ER.allo_states(site_reordered) == [:OnlyA, :NonequalAI]
        @test site == site_reordered
        @test hash(site) == hash(site_reordered)
        # ...but the pairing of ligand and state does: B::OnlyA, A::NonequalAI is
        # a different site.
        site_swapped = ER.RegulatorySite(
            [lig_b, lig_a], 4, [:OnlyA, :NonequalAI])
        @test ER.allo_states(site_swapped) == [:NonequalAI, :OnlyA]
        @test site != site_swapped
    end

    @testset "Step preserves iso direction (canonicalized in Mechanism ctor)" begin
        (; ES, EP) = _testhelper_uniuni()

        # The Step constructor does NOT canonicalize iso steps (RE or SS) —
        # iso direction depends on the mechanism's graph context and the step's
        # kinetic group, and is decided by `_canonical_step_direction` /
        # `_orient_tied_steps` in the Mechanism / AllostericMechanism constructor.
        # At the bare-Step level, direction is preserved.
        re_fwd = ER.Step(ES, EP, ER.Metabolite[], ER.Metabolite[], true)
        re_rev = ER.Step(EP, ES, ER.Metabolite[], ER.Metabolite[], true)
        @test re_fwd != re_rev
        @test ER.from_species(re_fwd) === ES
        @test ER.from_species(re_rev) === EP

        ss_fwd = ER.Step(ES, EP, ER.Metabolite[], ER.Metabolite[], false)
        ss_rev = ER.Step(EP, ES, ER.Metabolite[], ER.Metabolite[], false)
        @test ss_fwd != ss_rev
        @test ER.from_species(ss_fwd) === ES
        @test ER.from_species(ss_rev) === EP
    end

    @testset "Parameter family: step constants and Kreg" begin
        step = _testhelper_uniuni().bind

        # One RE constant and one forward/reverse rate pair serve every kind of step.
        for T in (ER.Krapid, ER.Kfor, ER.Krev)
            p = T(step, :None)
            @test p isa T
            @test p isa ER.Parameter
            @test p == T(step, :None)
            @test p != T(step, :I)
        end

        lig_a = ER.AllostericRegulator(:A)
        site = ER.RegulatorySite([lig_a], 2, [:OnlyA])
        kr = ER.Kreg(site, lig_a, :A)
        @test kr isa ER.Parameter
        @test kr == ER.Kreg(site, lig_a, :A)
    end

    @testset "ReactantAtoms canonicalizes atom ordering" begin
        ra1 = ER.ReactantAtoms(
            ER.Substrate(:ATP),
            [:C => 10, :H => 16, :N => 5],
        )
        ra2 = ER.ReactantAtoms(
            ER.Substrate(:ATP),
            [:N => 5, :H => 16, :C => 10],
        )
        @test ra1 == ra2
        @test hash(ra1) == hash(ra2)
        @test ER.metabolite(ra1) == ER.Substrate(:ATP)
        @test ER.atoms(ra1) == [:C => 10, :H => 16, :N => 5]

        # Distinct metabolite kinds: ATP-Substrate != ATP-Product
        ra_p = ER.ReactantAtoms(
            ER.Product(:ATP), [:C => 10, :H => 16, :N => 5])
        @test ra1 != ra_p
    end

    @testset "ReactantAtoms validation" begin
        # Mandatory atoms: empty atom list is rejected.
        @test_throws ErrorException ER.ReactantAtoms(
            ER.Substrate(:S), Pair{Symbol,Int}[])
        # Positive counts: zero / negative rejected.
        @test_throws ErrorException ER.ReactantAtoms(
            ER.Substrate(:S), [:C => 0])
        @test_throws ErrorException ER.ReactantAtoms(
            ER.Substrate(:S), [:C => -1])
        # Bool is not a valid count. Use untyped vector so Bool survives to
        # the constructor (Pair{Symbol,Int}[:C => true] would convert true→1).
        @test_throws ErrorException ER.ReactantAtoms(
            ER.Substrate(:S), [:C => true])
        # Valid construction still works.
        ra = ER.ReactantAtoms(ER.Substrate(:S), [:C => 6, :H => 12])
        @test ER.atoms(ra) == [:C => 6, :H => 12]
    end

    @testset "EnzymeReaction validation" begin
        S = ER.ReactantAtoms(ER.Substrate(:S), [:C => 1])
        P = ER.ReactantAtoms(ER.Product(:P), [:C => 1])
        P2 = ER.ReactantAtoms(ER.Product(:P), [:C => 1])
        Punbal = ER.ReactantAtoms(ER.Product(:P), [:C => 2])
        noregs = ER.RegulatorMults[]
        regS = ER.RegulatorMults(ER.CompetitiveInhibitor(:S), [1])
        regA = ER.RegulatorMults(ER.CompetitiveInhibitor(:A), [1])
        regA2 = ER.RegulatorMults(ER.CompetitiveInhibitor(:A), [1])

        # Empty substrate set rejected (only a product in reactants).
        @test_throws ErrorException ER.EnzymeReaction([P], noregs, Int[1])
        # Empty product set rejected (only a substrate in reactants).
        @test_throws ErrorException ER.EnzymeReaction([S], noregs, Int[1])
        # Duplicate product names rejected.
        @test_throws ErrorException ER.EnzymeReaction([S, P, P2], noregs, Int[1])
        # Atom imbalance rejected (S has 1 C, Punbal has 2 C).
        @test_throws ErrorException ER.EnzymeReaction([S, Punbal], noregs, Int[1])
        # Two regulators of the same kind with the same name rejected.
        @test_throws ErrorException ER.EnzymeReaction(
            [S, P], [regA, regA2], Int[1])
        # A regulator MAY share a substrate/product name (`::Inh` role tag:
        # one metabolite binds as a CompetitiveInhibitor under its real name).
        rxn_inh = ER.EnzymeReaction([S, P], [regS], Int[1])
        @test :S in ER.name.(ER.substrates(rxn_inh))
        # Valid balanced reaction with a distinct-named regulator still constructs.
        rxn = ER.EnzymeReaction([S, P], [regA], Int[1])
        @test ER.name.(ER.substrates(rxn)) == [:S]
        @test ER.name.(ER.products(rxn)) == [:P]
        # A name listed as both a substrate and a product is rejected, whether it
        # is the only reactant name (X → X) or one of several (A + B → A + C).
        Xs = ER.ReactantAtoms(ER.Substrate(:X), [:C => 2])
        Xp = ER.ReactantAtoms(ER.Product(:X), [:C => 2])
        @test_throws ErrorException(
            "EnzymeReaction: X named as both a substrate and a product; " *
            "concentrations and constants are keyed by name") ER.EnzymeReaction(
            [Xs, Xp], noregs, Int[1])
        As = ER.ReactantAtoms(ER.Substrate(:A), [:C => 1])
        Bs = ER.ReactantAtoms(ER.Substrate(:B), [:N => 1])
        Ap = ER.ReactantAtoms(ER.Product(:A), [:C => 1])
        Cp = ER.ReactantAtoms(ER.Product(:C), [:N => 1])
        @test_throws ErrorException(
            "EnzymeReaction: A named as both a substrate and a product; " *
            "concentrations and constants are keyed by name") ER.EnzymeReaction(
            [As, Bs, Ap, Cp], noregs, Int[1])
        # A single name MAY be declared in BOTH regulator roles: one
        # AllostericRegulator and one CompetitiveInhibitor. The two roles
        # render to distinct parameter names, so both are kept.
        regAllo = ER.RegulatorMults(
            ER.AllostericRegulator(:A), [1])
        rxn_dual = ER.EnzymeReaction([S, P], [regA, regAllo], Int[1])
        dual_regs = ER.regulators(rxn_dual)
        @test length(dual_regs) == 2
        @test Set(typeof(ER.regulator(rm)) for rm in dual_regs) ==
              Set([ER.AllostericRegulator,
                   ER.CompetitiveInhibitor])
        @test all(ER.name(ER.regulator(rm)) == :A
                  for rm in dual_regs)
        # Two allosteric regulators of the same name still rejected.
        regAllo2 = ER.RegulatorMults(
            ER.AllostericRegulator(:A), [1])
        @test_throws ErrorException ER.EnzymeReaction(
            [S, P], [regAllo, regAllo2], Int[1])
    end

    @testset "EnzymeReaction shared_catalytic_site" begin
        A = ER.ReactantAtoms(ER.Substrate(:A), [:C => 1])
        B = ER.ReactantAtoms(ER.Substrate(:B), [:C => 1])
        P = ER.ReactantAtoms(ER.Product(:P), [:C => 1])
        Q = ER.ReactantAtoms(ER.Product(:Q), [:C => 1])
        noregs = ER.RegulatorMults[]

        # Accepts a valid (substrate, product) pair; stored normalized + sorted.
        r = ER.EnzymeReaction([A, B, P, Q], noregs, Int[1];
            shared_catalytic_site = [(:A, :P)])
        @test ER.shared_catalytic_site(r) == [(:A, :P)]

        # Normalizes product-first input to (substrate, product).
        r2 = ER.EnzymeReaction([A, B, P, Q], noregs, Int[1];
            shared_catalytic_site = [(:P, :A)])
        @test ER.shared_catalytic_site(r2) == [(:A, :P)]

        # Sorts multiple pairs.
        r3 = ER.EnzymeReaction([A, B, P, Q], noregs, Int[1];
            shared_catalytic_site = [(:B, :Q), (:A, :P)])
        @test ER.shared_catalytic_site(r3) == [(:A, :P), (:B, :Q)]

        # Default is empty.
        @test ER.shared_catalytic_site(
            ER.EnzymeReaction([A, P], noregs, Int[1])) ==
            Tuple{Symbol,Symbol}[]

        # Rejections.
        @test_throws ErrorException ER.EnzymeReaction(  # unknown name
            [A, B, P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :Z)])
        @test_throws ErrorException ER.EnzymeReaction(  # two substrates
            [A, B, P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :B)])
        @test_throws ErrorException ER.EnzymeReaction(  # two products
            [A, B, P, Q], noregs, Int[1]; shared_catalytic_site = [(:P, :Q)])
        @test_throws ErrorException ER.EnzymeReaction(  # duplicate pair
            [A, B, P, Q], noregs, Int[1]; shared_catalytic_site = [(:A, :P), (:P, :A)])

        # A metabolite that is BOTH a substrate and a competitive inhibitor may
        # still be the substrate side of a shared pair (validation keys on the
        # substrate/product role, not the inhibitor role).
        regA = ER.RegulatorMults(ER.CompetitiveInhibitor(:A), [1])
        r_dual = ER.EnzymeReaction([A, B, P, Q], [regA], Int[1];
            shared_catalytic_site = [(:A, :P)])
        @test ER.shared_catalytic_site(r_dual) == [(:A, :P)]

        # == / hash reflect the field.
        @test r != ER.EnzymeReaction([A, B, P, Q], noregs, Int[1])
        @test hash(r) != hash(ER.EnzymeReaction([A, B, P, Q], noregs, Int[1]))
        @test r == ER.EnzymeReaction([A, B, P, Q], noregs, Int[1];
            shared_catalytic_site = [(:A, :P)])
    end

    @testset "RegulatorMults canonicalizes ordering + validates" begin
        rm1 = ER.RegulatorMults(
            ER.AllostericRegulator(:A), [4, 1, 2])
        rm2 = ER.RegulatorMults(
            ER.AllostericRegulator(:A), [1, 2, 4])
        @test rm1 == rm2
        @test hash(rm1) == hash(rm2)
        @test ER.regulator(rm1) ==
              ER.AllostericRegulator(:A)
        @test ER.allowed_multiplicities(rm1) == [1, 2, 4]

        # Multiplicity < 1 → error
        @test_throws Exception ER.RegulatorMults(
            ER.AllostericRegulator(:A), [0, 1])
        @test_throws Exception ER.RegulatorMults(
            ER.AllostericRegulator(:A), [-1])

        # CompetitiveInhibitor also accepted
        rm_ci = ER.RegulatorMults(
            ER.CompetitiveInhibitor(:I), [1])
        @test ER.regulator(rm_ci) ==
              ER.CompetitiveInhibitor(:I)
    end

    @testset "RegulatorMults reg_type" begin
        rm0 = ER.RegulatorMults(ER.AllostericRegulator(:X), [2])
        @test ER.reg_type(rm0) == :unspecified      # default
        rmA = ER.RegulatorMults(
            ER.AllostericRegulator(:X), [2], :activator)
        @test ER.reg_type(rmA) == :activator
        @test rm0 != rmA                                      # reg_type participates in ==
        @test rm0 ==
              ER.RegulatorMults(ER.AllostericRegulator(:X), [2])
        @test hash(rmA) != hash(rm0)
        @test_throws ErrorException ER.RegulatorMults(
            ER.AllostericRegulator(:X), [2], :bogus)
    end

    @testset "EnzymeReaction canonicalizes reactant + regulator ordering" begin
        r1 = EnzymeReaction(
            [ER.ReactantAtoms(
                 ER.Substrate(:B), [:C => 1]),
             ER.ReactantAtoms(
                 ER.Substrate(:A), [:C => 1]),
             ER.ReactantAtoms(
                 ER.Product(:P), [:C => 2])],
            [ER.RegulatorMults(
                 ER.AllostericRegulator(:Y), [2]),
             ER.RegulatorMults(
                 ER.AllostericRegulator(:X), [2])],
            [3, 1, 2],
        )
        r2 = EnzymeReaction(
            [ER.ReactantAtoms(
                 ER.Substrate(:A), [:C => 1]),
             ER.ReactantAtoms(
                 ER.Substrate(:B), [:C => 1]),
             ER.ReactantAtoms(
                 ER.Product(:P), [:C => 2])],
            [ER.RegulatorMults(
                 ER.AllostericRegulator(:X), [2]),
             ER.RegulatorMults(
                 ER.AllostericRegulator(:Y), [2])],
            [1, 2, 3],
        )
        @test r1 == r2
        @test hash(r1) == hash(r2)
    end

    @testset "EnzymeReaction rejects multiplicity < 1" begin
        @test_throws ErrorException EnzymeReaction(
            ER.ReactantAtoms[
                ER.ReactantAtoms(
                    ER.Substrate(:S), [:C => 1])],
            ER.RegulatorMults[],
            [0, 1],
        )
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
        iso = only(s for grp in ER.steps(ER.Mechanism(s2))
                       for s in grp
                       if ER.is_iso(s) &&
                          ER.name(ER.from_species(s)) in (:E, :F))
        @test ER.name(ER.from_species(iso)) == :F  # product-exit
        @test ER.name(ER.to_species(iso))   == :E  # substrate-entry
    end

    @testset "Mechanism canonicalizes step order" begin
        (; rxn, bind, iso, rel) = _testhelper_uniuni()

        m = ER.Mechanism(rxn, [[bind], [iso], [rel]])
        @test ER.reaction(m) == rxn
        flat = ER._flat_steps(m)
        @test Set(s for (s, _) in flat) == Set([bind, iso, rel])
        @test [g for (_, g) in flat] == [1, 2, 3]
        @test sum(length, ER.steps(m)) == 3
        # Canonical by construction: a permuted input yields identical storage.
        m_perm = ER.Mechanism(rxn, [[rel], [bind], [iso]])
        @test ER._flat_steps(m) == ER._flat_steps(m_perm)
        @test m == m_perm && hash(m) == hash(m_perm)
    end

    @testset "AllostericMechanism (non-parametric)" begin
        (; rxn, bind, iso, rel) = _testhelper_uniuni(mults = [2])

        site = ER.RegulatorySite(
            [ER.AllostericRegulator(:I)], 1, [:OnlyI])
        m = ER.AllostericMechanism(
            rxn, [[bind], [iso], [rel]],
            [:EqualAI, :OnlyA, :NonequalAI], 2,
            [site])

        @test ER.reaction(m) == rxn
        # Steps canonicalized; allosteric tags stay bound to their steps.
        @test Set(only(g) for g in ER.steps(m)) ==
              Set([bind, iso, rel])
        state_of(step) = ER.cat_allo_state(m,
            only(g for g in eachindex(ER.steps(m))
                 if first(ER.steps(m)[g]) == step))
        @test state_of(bind) == :EqualAI
        @test state_of(iso)  == :OnlyA
        @test state_of(rel)  == :NonequalAI
        @test ER.catalytic_multiplicity(m) == 2
        @test ER.regulatory_sites(m) == [site]
        @test sum(length, ER.steps(m)) == 3
        @test iso in first.(ER.steps(m))
        @test ER.allosteric_regulators(m) ==
              [ER.AllostericRegulator(:I)]

        m2 = ER.AllostericMechanism(
            rxn, [[bind], [iso], [rel]],
            [:EqualAI, :OnlyA, :NonequalAI], 2,
            [site])
        @test m == m2
        @test hash(m) == hash(m2)
    end

    @testset "AllostericMechanism validation errors" begin
        (; rxn, bind) = _testhelper_uniuni()
        cat_steps = [[bind]]

        # :OnlyI for catalytic group is rejected (R-state-active convention)
        @test_throws ErrorException ER.AllostericMechanism(
            rxn, cat_steps, [:OnlyI], 1,
            ER.RegulatorySite[])

        # Length mismatch
        @test_throws ErrorException ER.AllostericMechanism(
            rxn, cat_steps, [:EqualAI, :NonequalAI], 1,
            ER.RegulatorySite[])

        # catalytic_multiplicity < 1
        @test_throws ErrorException ER.AllostericMechanism(
            rxn, cat_steps, [:EqualAI], 0,
            ER.RegulatorySite[])

        # Unknown allo-state symbol
        @test_throws ErrorException ER.AllostericMechanism(
            rxn, cat_steps, [:Bogus], 1,
            ER.RegulatorySite[])
    end

    @testset "AllostericMechanism(::AllostericEnzymeMechanism) converter" begin
        cm = _testhelper_re_mm
        aem = _testhelper_aem(
            cm, (2, (:EqualAI, :NonequalAI, :OnlyA)),
            (((:A, :B), 1, (:OnlyA, :NonequalAI)),),
        )

        @test ER.catalytic_mechanism(aem) === cm
        @test ER.catalytic_multiplicity(aem) == 2

        am = ER.AllostericMechanism(aem)
        @test am isa ER.AllostericMechanism

        # Catalytic side lifted via Mechanism(CM())
        @test ER.steps(am) == ER.Mechanism(cm).steps
        @test ER.reaction(am) == ER.Mechanism(cm).reaction

        # cat_allo_states and multiplicity extracted from CS
        @test ER.cat_allo_state(am, 1) === :EqualAI
        @test ER.cat_allo_state(am, 2) === :NonequalAI
        @test ER.cat_allo_state(am, 3) === :OnlyA
        @test ER.catalytic_multiplicity(am) == 2

        # Regulatory sites: ligand Symbols wrapped as AllostericRegulator
        sites = ER.regulatory_sites(am)
        @test length(sites) == 1
        @test sites[1].ligands ==
              [ER.AllostericRegulator(:A),
               ER.AllostericRegulator(:B)]
        @test sites[1].multiplicity == 1
        @test sites[1].allo_states == [:OnlyA, :NonequalAI]

        # Idempotent: two calls give equal results
        @test ER.AllostericMechanism(aem) == am
    end

    @testset "name(p, ::AllostericMechanism) chokepoint" begin
        cm = _testhelper_re_mm
        aem = _testhelper_aem(
            cm, (2, (:NonequalAI, :EqualAI, :NonequalAI)),
            (((:R,), 1, (:NonequalAI,)),),
        )
        am = ER.AllostericMechanism(aem)

        # Step-bound parameters. Group order is canonical, so pick the
        # substrate-binding and iso steps by content.
        reps = first.(ER.steps(am))
        rep_bind = only(r for r in reps if ER.bound_metabolite(r) isa ER.Substrate)
        @test ER.name(ER.Krapid(rep_bind, :None), am) === :K_ES_to_E_S
        @test ER.name(ER.Krapid(rep_bind, :I), am) === :K_I_ES_to_E_S

        rep_iso  = only(r for r in reps if ER.bound_metabolite(r) === nothing)
        @test ER.name(ER.Kfor(rep_iso, :None), am) === :k_ES_to_EP
        @test ER.name(ER.Krapid(rep_iso, :None), am) === :K_ES_to_EP
        @test ER.name(ER.Krev(rep_iso, :None), am) === :k_EP_to_ES

        site = ER.regulatory_sites(am)[1]
        lig  = first(site.ligands)
        @test ER.name(ER.Kreg(site, lig, :A), am) === :K_A_Rreg
        @test ER.name(ER.Kreg(site, lig, :I), am) === :K_I_Rreg

        # The chokepoint renders on the concrete mechanism only; a compiled type
        # is lifted first.
        @test_throws MethodError ER.name(ER.Krapid(rep_bind, :None), aem)
        @test_throws MethodError ER.name(ER.Kreg(site, lig, :A), aem)
    end

    @testset "_sig_of / _mechanism_from_sig roundtrip" begin
        # One fixture holds every Sig leaf kind: all four metabolite tags (a
        # competitive inhibitor and an allosteric regulator each bound by a step), a
        # covalent residual, a Theorell–Chance step and multi-element atom data.
        # Pair{Symbol,Int} is NOT a valid type-parameter value; atoms must encode
        # as Tuple{Symbol,Int} leaves.
        A, B = ER.Substrate(:A), ER.Substrate(:B)
        P, Q = ER.Product(:P), ER.Product(:Q)
        I, R = ER.CompetitiveInhibitor(:I), ER.AllostericRegulator(:R)
        rxn(regs) = ER.EnzymeReaction(
            [ER.ReactantAtoms(A, [:C => 2, :X => 1]), ER.ReactantAtoms(B, [:N => 1]),
             ER.ReactantAtoms(P, [:C => 2]), ER.ReactantAtoms(Q, [:N => 1, :X => 1])],
            regs, [1, 2])
        E, EA, EQ = _testhelper_sp([]), _testhelper_sp([A]), _testhelper_sp([Q])
        F = _testhelper_sp([], :E, ER.Residual([A], [P]))
        none = ER.Metabolite[]
        groups = [[ER.Step(E, EA, [A], none, true)],
                  [ER.Step(EA, F, none, [P], false)],
                  [ER.Step(F, EQ, [B], none, false)],
                  [ER.Step(EA, EQ, [B], [P], false)],                  # Theorell–Chance
                  [ER.Step(E, EQ, [Q], none, true)],
                  [ER.Step(E, _testhelper_sp([I]), [I], none, true)],
                  [ER.Step(EA, _testhelper_sp([A, R]), [R], none, true)]]
        regs = [ER.RegulatorMults(I, [1]), ER.RegulatorMults(R, [1, 2])]
        m = ER.Mechanism(rxn(regs), groups)

        for x in (A, P, I, R)
            @test ER._to_sig(x) == (nameof(typeof(x)), ER.name(x))
        end
        @test ER._to_sig(F) == ((), :E, (((:Substrate, :A),), ((:Product, :P),)))

        sig = ER._sig_of(m)
        @test sig isa Tuple
        @test length(sig) == 2   # (reaction_sig, steps_sig)
        @test ER._mechanism_from_sig(sig) == m   # roundtrip

        # CRITICAL: sig MUST be usable as a type parameter. Throws TypeError
        # if any leaf is invalid (Pair, Vector, DataType inside value-tuple).
        em_type = ER.EnzymeMechanism{sig}
        @test em_type <: ER.EnzymeMechanism
        em_inst = em_type()
        @test ER.EnzymeMechanism(m) === em_inst

        # Roundtrip through the type-parameter form preserves everything.
        @test ER.Mechanism(em_inst) == m

        # A regulator that no step binds is left out of the Sig, so the mechanism
        # declaring it compiles to the same type as the one without it.
        U = ER.CompetitiveInhibitor(:U)
        m_unbound = ER.Mechanism(rxn([regs; ER.RegulatorMults(U, [1])]), groups)
        @test ER._sig_of(m_unbound) == sig
        @test ER.EnzymeMechanism(m_unbound) === em_inst
    end

    @testset "name(p::Parameter, m) chokepoint" begin
        (; rxn, bind, iso, rel) = _testhelper_uniuni()
        m = ER.Mechanism(rxn, [[bind], [iso], [rel]])

        # Structural naming: every step constant encodes its reaction's two sides;
        # a binding K reads in the release direction, iso params in the stored one.
        @test ER.name(ER.Krapid(bind, :None), m) === :K_ES_to_E_S
        @test ER.name(ER.Krapid(bind, :I),    m) === :K_I_ES_to_E_S
        @test ER.name(ER.Kfor(iso, :None), m) === :k_ES_to_EP
        @test ER.name(ER.Krev(iso, :None), m) === :k_EP_to_ES
        @test ER.name(ER.Krapid(rel, :None), m) === :K_EP_to_E_P

        # A binding's rate pair: the forward rate binds, the reverse rate releases
        @test ER.name(ER.Kfor(bind, :None), m) === :k_E_S_to_ES
        @test ER.name(ER.Krev(bind, :None), m) === :k_ES_to_E_S

        # I-state token on SS step
        @test ER.name(ER.Kfor(iso, :I), m) === :k_I_ES_to_EP
        @test ER.name(ER.Krev(iso, :I), m) === :k_I_EP_to_ES

        # The RE iso's equilibrium constant, named in the stored direction
        @test ER.name(ER.Krapid(iso, :None), m) === :K_ES_to_EP
        @test ER.name(ER.Krapid(iso, :I),    m) === :K_I_ES_to_EP

        # Same names resolve on the mechanism lifted back from EnzymeMechanism(m);
        # the compiled type itself is not a chokepoint argument.
        em = EnzymeMechanism(m)
        @test ER.name(ER.Krapid(bind, :None), ER.Mechanism(em)) === :K_ES_to_E_S
        @test ER.name(ER.Kfor(iso, :None), ER.Mechanism(em)) === :k_ES_to_EP
        @test_throws MethodError ER.name(ER.Krapid(bind, :None), em)
    end

    @testset "name(p::Parameter, m) for the steps of a shared kinetic group" begin
        (; rxn, EP, bind, iso, rel) = _testhelper_uniuni()
        bind_into_EP = ER.Step(EP, _testhelper_sp([ER.Substrate(:S), ER.Product(:P)]),
                               [ER.Substrate(:S)], ER.Metabolite[], true)
        m = ER.Mechanism(rxn, [[bind, bind_into_EP], [iso], [rel]])

        # Both steps bind S; rep = bind. Both yield the same name.
        @test ER.name(ER.Krapid(bind, :None), m) === :K_ES_to_E_S
        @test ER.name(ER.Krapid(bind_into_EP, :None), m) === :K_ES_to_E_S
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

@testset "Step: explicit consumed/released lists" begin
    A, B, P, Q = ER.Substrate(:A), ER.Substrate(:B), ER.Product(:P), ER.Product(:Q)
    E, EA = _testhelper_sp([]), _testhelper_sp([A])
    EQ, EAB = _testhelper_sp([Q]), _testhelper_sp([A, B])

    @testset "plain binding keeps its written orientation" begin
        s = ER.Step(E, EA, [A], ER.Metabolite[], true)
        @test fieldnames(ER.Step) ==
              (:from_species, :to_species, :consumed, :released, :is_equilibrium)
        @test ER.from_species(s) === E && ER.to_species(s) === EA
        @test ER.consumed(s) == ER.Metabolite[A] && isempty(ER.released(s))
        @test ER.bound_metabolite(s) == A && ER.is_binding(s) && !ER.is_iso(s)
        @test ER.is_equilibrium(s) && !ER._is_chemistry(s)
    end

    @testset "a plain release is stored as the binding it reverses" begin
        for is_eq in (true, false)
            release = ER.Step(EA, E, ER.Metabolite[], [A], is_eq)
            binding = ER.Step(E, EA, [A], ER.Metabolite[], is_eq)
            @test release == binding && hash(release) == hash(binding)
            @test ER.from_species(release) === E && ER.to_species(release) === EA
        end
        # conformation change allowed: E*(A) → E + A is the binding E + A → E*(A)
        Estar_A = _testhelper_sp([A], :Estar)
        r = ER.Step(Estar_A, E, ER.Metabolite[], [A], true)
        @test ER.from_species(r) === E && ER.to_species(r) === Estar_A
        @test ER.bound_metabolite(r) == A && !ER._is_chemistry(r)
    end

    @testset "isomerization and transformations" begin
        iso = ER.Step(EAB, _testhelper_sp([P, Q]), ER.Metabolite[], ER.Metabolite[], false)
        @test ER.is_iso(iso) && ER.bound_metabolite(iso) === nothing && !ER.is_binding(iso)
        @test ER._is_chemistry(iso)
        # Chemistry + release, stored as the fused binding of P it reverses.
        fused = ER.Step(EAB, EQ, ER.Metabolite[], [P], false)
        @test !ER.is_iso(fused) && ER.bound_metabolite(fused) == P
        @test ER.from_species(fused) == EQ && ER.to_species(fused) == EAB
        @test ER.consumed(fused) == ER.Metabolite[P]
        @test isempty(ER.released(fused)) && ER._is_chemistry(fused)
        tc = ER.Step(EA, EQ, [B], [P], false)                          # Theorell–Chance
        @test ER.consumed(tc) == ER.Metabolite[B] && ER.released(tc) == ER.Metabolite[P]
        @test ER.bound_metabolite(tc) === nothing && ER._is_chemistry(tc)
        two = ER.Step(E, EAB, [B, A], ER.Metabolite[], true)           # lists are sorted
        @test ER.consumed(two) == ER.Metabolite[A, B] &&
              ER.bound_metabolite(two) === nothing
        # A fused binding: the metabolite is taken up while chemistry runs.
        fused_binding = ER.Step(EA, _testhelper_sp([P, Q]), [B], ER.Metabolite[], true)
        @test ER.bound_metabolite(fused_binding) == B && ER._is_chemistry(fused_binding)
    end

    @testset "covalent residual: binding onto a residual form vs chemistry" begin
        res = ER.Residual([A], [P])
        F, FB = _testhelper_sp([], :E, res), _testhelper_sp([B], :E, res)
        onto_residual = ER.Step(F, FB, [B], ER.Metabolite[], true)
        @test ER.bound_metabolite(onto_residual) == B && !ER._is_chemistry(onto_residual)
        # E(A) → F + P is stored as the fused binding F + P → E(A) it reverses.
        chem = ER.Step(EA, F, ER.Metabolite[], [P], false)
        @test ER.bound_metabolite(chem) == P && ER.from_species(chem) == F
        @test ER._is_chemistry(chem)
        m = ER.Mechanism(
            @enzyme_reaction(begin
                substrates: A[CX], B[N]
                products: P[C], Q[NX]
            end),
            [[ER.Step(E, EA, [A], ER.Metabolite[], false)], [chem],
             [ER.Step(F, FB, [B], ER.Metabolite[], false)],
             [ER.Step(FB, E, ER.Metabolite[], [Q], false)]])
        # The Mechanism constructor accepts the cycle through the residual form.
        @test m isa ER.Mechanism
    end

    @testset "inhibitor copy stays distinct from the substrate" begin
        Ai = ER.CompetitiveInhibitor(:A)
        s_sub = ER.Step(E, EA, [A], ER.Metabolite[], true)
        s_inh = ER.Step(E, _testhelper_sp([Ai]), [Ai], ER.Metabolite[], true)
        @test s_sub != s_inh && hash(s_sub) != hash(s_inh)
        @test ER.bound_metabolite(s_inh) == Ai
    end

    @testset "a binding keeps each metabolite's role" begin
        Pi, S = ER.CompetitiveInhibitor(:P), ER.Substrate(:S)
        EP, EPi, EPP = _testhelper_sp([P]), _testhelper_sp([Pi]), _testhelper_sp([P, P])
        role_change(from, to, m) = ErrorException(
            "Step $from → $to changes a metabolite's role: $to holds $from's " *
            "metabolites plus $m by name, not by role; write an inhibitor copy as " *
            "X::Inh on both sides")
        # The copy of P binds into the product's form, written as a binding or a release.
        @test_throws role_change(:E, :EP, :P) ER.Step(E, EP, [Pi], ER.Metabolite[], true)
        @test_throws role_change(:E, :EP, :P) ER.Step(EP, E, ER.Metabolite[], [Pi], true)
        # The product binds into the copy's form.
        @test_throws role_change(:E, :EPinh, :P) ER.Step(
            E, EPi, [P], ER.Metabolite[], true)
        # A second P, taken up as the copy, lands in a form holding two products.
        @test_throws role_change(:EP, :EPP, :P) ER.Step(
            EP, EPP, [Pi], ER.Metabolite[], true)
        # Each metabolite keeps its role: accepted.
        @test ER.bound_metabolite(
            ER.Step(EP, _testhelper_sp([P, Pi]), [Pi], ER.Metabolite[], true)) == Pi
        # An isomerization may change how often a name is bound: E(P, S) → E(P, P).
        @test ER.is_iso(ER.Step(_testhelper_sp([P, S]), EPP, ER.Metabolite[],
                                ER.Metabolite[], false))
        # The DSL spellings, with P declared as an inhibitor.
        for (line, from, to) in ((:(E + P::Inh ⇌ E(P)), :E, :EP),
                                 (:(E(P) ⇌ E + P::Inh), :E, :EP),
                                 (:(E(P) + P::Inh ⇌ E(P, P)), :EP, :EPP))
            @test_throws role_change(from, to, :P) eval(:(@enzyme_mechanism begin
                substrates: S
                products:   P
                regulators: P
                steps: begin
                    E + S ⇌ E(S)
                    E(S) <--> E(P)
                    E(P) ⇌ E + P
                    $line
                end
            end))
        end
    end

    @testset "rejections" begin
        err = _testhelper_thrown(() -> ER.Step(E, E, [A], ER.Metabolite[], true))
        @test err isa ErrorException
        @test occursin("both ends", err.msg)
        err = _testhelper_thrown(() -> ER.Step(EA, EQ, [A], [A], false))
        @test err isa ErrorException
        @test occursin("both consumed and released", err.msg)
    end

    @testset "sort key orders bindings and isomerizations" begin
        s = ER.Step(E, EA, [A], ER.Metabolite[], true)
        @test ER._step_canonical_key(s) == ("E", "EA", "A", "", true)
        iso = ER.Step(EAB, _testhelper_sp([P, Q]), ER.Metabolite[], ER.Metabolite[], false)
        @test ER._step_canonical_key(iso) == ("EAB", "EPQ", "", "", false)
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
        # A and its inhibitor copy A::Inh are different metabolites, so their
        # bindings are different kinds of step and cannot share a constant.
        err = _testhelper_thrown() do
            @enzyme_mechanism begin
                substrates: A
                products: Q
                regulators: A
                steps: begin
                    (E + A ⇌ E(A), E(Q) + A::Inh ⇌ E(A::Inh, Q))
                    E(A) <--> E(Q)
                    E(Q) <--> E + Q
                end
            end
        end
        @test err isa ErrorException
        @test occursin("a kinetic group holds E_A → EA (RE) and EQ_Ainh → EAinhQ (RE)",
                       err.msg)
    end

    @testset "an empty kinetic group" begin
        A = ER.Substrate(:A)
        binding = ER.Step(_testhelper_sp([]), _testhelper_sp([A]), [A], ER.Metabolite[],
                          true)
        rxn = @enzyme_reaction(begin
            substrates: A[C]
            products: P[C]
        end)
        err = _testhelper_thrown(() -> ER.Mechanism(rxn, [[binding], ER.Step[]]))
        @test err isa ErrorException
        @test occursin("kinetic group 2 is empty", err.msg)
    end

    @testset "accepted: context-shared bindings of one metabolite" begin
        # R binds both free E and E(S) with the same K (non-competitive
        # inhibitor pattern): a legitimate shared binding of one metabolite.
        m = @enzyme_mechanism begin
            substrates: S
            products: P
            regulators: R
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
        A = ER.Substrate(:A)
        e, e2   = _testhelper_sp([], :E), _testhelper_sp([], :Estar)
        ea, e2a = _testhelper_sp([A], :E), _testhelper_sp([A], :Estar)
        iso1 = ER.Step(e, e2, ER.Metabolite[], ER.Metabolite[], false)
        iso2 = ER.Step(ea, e2a, ER.Metabolite[], ER.Metabolite[], false)
        rxn = @enzyme_reaction(begin
            substrates: A[C]
            products: P[C]
        end)
        m = ER.Mechanism(rxn, [[iso1, iso2]])
        @test m isa ER.Mechanism
    end
end

@testset "AllostericMechanism enforces the kinetic-group rules" begin
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: NonequalAI
            E(S) <--> E(P)   :: EqualAI
            E(P) ⇌ E + P     :: EqualAI
        end
    end)
    group_of(pred) = only(g for g in ER.steps(am) if pred(first(g)))
    bind_s = group_of(s -> ER.bound_metabolite(s) == ER.Substrate(:S))
    bind_p = group_of(s -> ER.bound_metabolite(s) == ER.Product(:P))
    iso = group_of(ER.is_iso)
    rebuild(groups) = _testhelper_thrown() do
        ER.AllostericMechanism(ER.reaction(am), groups, fill(:EqualAI, length(groups)),
                               2, ER.RegulatorySite[])
    end
    # The two bindings in one group: different metabolites.
    err = rebuild([[bind_s; bind_p], iso])
    @test err isa ErrorException
    @test occursin("a kinetic group holds E_P → EP (RE) and E_S → ES (RE)", err.msg)
    # The isomerization in two groups.
    err = rebuild([bind_s, iso, iso, bind_p])
    @test err isa ErrorException
    @test occursin("both hold the reaction ES ⇌ EP", err.msg)
    # The isomerization twice in one group.
    err = rebuild([bind_s, [iso; iso], bind_p])
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
    names = Set(ER.parameters(m, ER.Full))
    @test :k_EAB_to_EQ_P in names && :k_EQ_P_to_EAB in names
    # The fused release E(A, B) → E(Q) + P is stored as the binding of P it
    # reverses, E(Q) + P → E(A, B), so its forward rate constant binds P.
    mech = ER.Mechanism(m)
    i = only(i for (i, (s, _)) in enumerate(ER._flat_steps(mech)) if ER._is_chemistry(s))
    kon, koff = ER._step_parameters(mech)[i]
    @test kon isa ER.Kfor && koff isa ER.Krev
    @test ER.name(kon, mech) == :k_EQ_P_to_EAB && ER.name(koff, mech) == :k_EAB_to_EQ_P
    @test :k_E_A_to_EA in names && :k_EQ_to_E_Q in names
end

@testset "step constants are named by their reaction" begin
    full_names(m) = Set(ER.parameters(m, ER.Full))
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
                   full_names(mm_ss))
    # Reduced mode: the Haldane reduction leaves the SS isomerization constant fitted.
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
    @test :K_ES_to_EP in full_names(re_iso)
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
    @test :K_EAB_to_EQ_P in full_names(re_fused)
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
    @test :K_EP_to_E_S in full_names(re_fused_binding)
    # A competitive-inhibitor copy of A and A itself bind E without colliding.
    inh = @enzyme_mechanism begin
        substrates: A
        products: P
        regulators: A
        steps: begin
            E + A ⇌ E(A)
            E(A) <--> E(P)
            E(P) ⇌ E + P
            E + A::Inh ⇌ E(A::Inh)
        end
    end
    @test issubset([:K_EAinh_to_E_Ainh, :K_EA_to_E_A], full_names(inh))
    # The allosteric state tag follows the prefix: a :NonequalAI binding has an
    # A-state and an I-state constant, an :EqualAI binding one shared constant.
    allo = @allosteric_mechanism begin
        substrates: S
        products: P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: NonequalAI
            E(S) <--> E(P)   :: EqualAI
            E(P) ⇌ E + P     :: EqualAI
        end
    end
    # Every constant the reduced rate equation names, fitted or dependent.
    dep, indep = ER._dependent_param_exprs(typeof(allo))
    allo_names = union(keys(dep), indep)
    @test issubset([:K_A_ES_to_E_S, :K_I_ES_to_E_S], allo_names)
    @test !(:K_ES_to_E_S in allo_names)
    @test :K_EP_to_E_P in ER.parameters(allo) && !(:K_A_EP_to_E_P in allo_names)
    @test :K_I_EP_to_E_P ∉ ER.parameters(allo)
    # A ping-pong residual form: B binds E(; residual = A - P), whose name is
    # E_res_+A_-P, so the release-direction K reads EB_res_+A_-P → E_res_+A_-P + B.
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
    pp_names = full_names(pingpong)
    @test Symbol("K_EB_res_+A_-P_to_E_res_+A_-P_B") in pp_names
    @test length(pp_names) == length(ER.parameters(pingpong, ER.Full))
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
    flat = ER._flat_steps(m)
    i = only(i for (i, (s, _)) in enumerate(flat) if ER.is_iso(s))
    iso = first(flat[i])
    @test ER.name(ER.from_species(iso)) == :F && ER.name(ER.to_species(iso)) == :E
    kf, kr = ER._step_parameters(m)[i]
    @test ER.name(kf, m) == :k_F_to_E && ER.name(kr, m) == :k_E_to_F
end

@testset "a competitive-inhibitor copy marks no reactant side" begin
    # P's competitive-inhibitor copy binds E. Tiers 1 and 2 read metabolite roles, so the
    # copy marks nothing at E, and the Iso Uni Uni isomerization still runs from F, where
    # P leaves, to E, where A enters.
    m = ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        regulators: P
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            E + P::Inh ⇌ E(P::Inh)
        end
    end)
    iso = only(s for g in ER.steps(m) for s in g
               if ER.is_iso(s) && ER.name(ER.from_species(s)) in (:E, :F))
    @test ER.name(ER.from_species(iso)) == :F && ER.name(ER.to_species(iso)) == :E
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
        regulators: I
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            (F <--> E, F(I) <--> E(I))
            (E + I ⇌ E(I), F + I ⇌ F(I))
        end
    end
    m = ER.Mechanism(ss)
    mirror = only(g for g in ER.steps(m) if length(g) == 2 && all(ER.is_iso, g))
    @test all(s -> ER.conformation(ER.from_species(s)) == :F &&
                   ER.conformation(ER.to_species(s)) == :E, mirror)
    @test !occursin("Wegscheider", rate_equation_string(ss))
    # The written direction and order of the mirrored steps do not matter.
    @test ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        regulators: I
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
        regulators: I
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            (F ⇌ E, F(I) ⇌ E(I))
            (E + I ⇌ E(I), F + I ⇌ F(I))
        end
    end
    @test :K_F_to_E in ER.fitted_params(re)

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

@testset "a kinetic group ties one metabolite exchange at two conformations" begin
    # Segel Iso Uni Uni with an inhibitor I binding E and F, and J displacing I at both
    # conformations in one group. Neither exchange changes conformation, so the group
    # stores both to take up J and give off I, however the steps are written.
    I, J = ER.CompetitiveInhibitor(:I), ER.CompetitiveInhibitor(:J)
    plain = @enzyme_mechanism begin
        substrates: A
        products: P
        regulators: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)
        end
    end
    allo = @allosteric_mechanism begin
        substrates: A
        products: P
        catalytic_inhibitors: I, J
        catalytic_steps: begin
            E + A <--> E(A)                                   :: EqualAI
            E(A) <--> E(P)                                    :: EqualAI
            E(P) <--> F + P                                   :: EqualAI
            F <--> E                                          :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))                      :: EqualAI
            (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)  :: NonequalAI
        end
    end
    m, am = ER.Mechanism(plain), ER.AllostericMechanism(allo)
    for mech in (m, am)
        exchange = only(g for g in ER.steps(mech)
                        if !any(s -> ER.is_iso(s) || ER.is_binding(s), g))
        @test length(exchange) == 2
        @test all(s -> ER.consumed(s) == ER.Metabolite[J] &&
                       ER.released(s) == ER.Metabolite[I], exchange)
    end
    # The written direction and order of the exchanges do not matter.
    @test ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        regulators: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)
        end
    end) == m
    @test ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: A
        products: P
        catalytic_inhibitors: I, J
        catalytic_steps: begin
            E + A <--> E(A)                                   :: EqualAI
            E(A) <--> E(P)                                    :: EqualAI
            E(P) <--> F + P                                   :: EqualAI
            F <--> E                                          :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))                      :: EqualAI
            (F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)  :: NonequalAI
        end
    end) == am
    # The rate law derives, with the exchange constants named in the J-in direction.
    @test occursin("k_EIinh_Jinh_to_EJinh_Iinh", rate_equation_string(plain))
    allo_law = rate_equation_string(allo)
    @test occursin("k_A_EIinh_Jinh_to_EJinh_Iinh", allo_law)
    @test occursin("k_I_EIinh_Jinh_to_EJinh_Iinh", allo_law)

    # Exchanges that switch conformation in opposite directions as J displaces I share no
    # orientation.
    @test_throws "make different conformational changes" @enzyme_mechanism begin
        substrates: A
        products: P
        regulators: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E(I) + J <--> F(J) + I, F(I) + J <--> E(J) + I)
        end
    end
end

@testset "a tied exchange turns to take up its lead's metabolites" begin
    # The exchange fixture above plus A binding E(I) and P binding E(J). On its own, Tier 2
    # runs the exchange at E from E(J), where P binds, to E(I), where A binds, while Tier 3
    # runs the exchange at F from F(I) to F(J). The exchange at E leads, so the group
    # turns the one at F to take up I and give off J as well.
    I, J = ER.CompetitiveInhibitor(:I), ER.CompetitiveInhibitor(:J)
    is_exchange(g) = !any(s -> ER.is_iso(s) || ER.is_binding(s), g)
    exchange_group(mech) = only(filter(is_exchange, ER.steps(mech)))
    takes_up_I(s) = ER.consumed(s) == ER.Metabolite[I] && ER.released(s) == ER.Metabolite[J]
    plain = @enzyme_mechanism begin
        substrates: A
        products: P
        regulators: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            E(I) + A ⇌ E(A, I)
            E(J) + P ⇌ E(J, P)
            (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)
        end
    end
    allo = @allosteric_mechanism begin
        substrates: A
        products: P
        catalytic_inhibitors: I, J
        catalytic_steps: begin
            E + A <--> E(A)                                   :: EqualAI
            E(A) <--> E(P)                                    :: EqualAI
            E(P) <--> F + P                                   :: EqualAI
            F <--> E                                          :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))                      :: EqualAI
            E(I) + A ⇌ E(A, I)                                :: EqualAI
            E(J) + P ⇌ E(J, P)                                :: EqualAI
            (E(I) + J <--> E(J) + I, F(I) + J <--> F(J) + I)  :: NonequalAI
        end
    end
    m, am = ER.Mechanism(plain), ER.AllostericMechanism(allo)
    for mech in (m, am)
        @test length(exchange_group(mech)) == 2
        @test all(takes_up_I, exchange_group(mech))
    end
    # The written direction and order of the exchanges do not matter.
    @test ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        regulators: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            E(I) + A ⇌ E(A, I)
            E(J) + P ⇌ E(J, P)
            (F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)
        end
    end) == m
    @test ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: A
        products: P
        catalytic_inhibitors: I, J
        catalytic_steps: begin
            E + A <--> E(A)                                   :: EqualAI
            E(A) <--> E(P)                                    :: EqualAI
            E(P) <--> F + P                                   :: EqualAI
            F <--> E                                          :: EqualAI
            (E + I ⇌ E(I), F + I ⇌ F(I))                      :: EqualAI
            E(I) + A ⇌ E(A, I)                                :: EqualAI
            E(J) + P ⇌ E(J, P)                                :: EqualAI
            (F(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)  :: NonequalAI
        end
    end) == am
    @test occursin("k_EJinh_Iinh_to_EIinh_Jinh", rate_equation_string(plain))
    allo_law = rate_equation_string(allo)
    @test occursin("k_A_EJinh_Iinh_to_EIinh_Jinh", allo_law)
    @test occursin("k_I_EJinh_Iinh_to_EIinh_Jinh", allo_law)

    # An exchange that changes conformation leads its group. On its own, Tier 3 runs
    # F(I) + J <--> E(J) + I from E(J) to F(I) and the exchange at E from E(I) to E(J);
    # the group turns the exchange at E to take up I and give off J as the lead does.
    mixed = @enzyme_mechanism begin
        substrates: A
        products: P
        regulators: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E(I) + J <--> E(J) + I, F(I) + J <--> E(J) + I)
        end
    end
    mm = ER.Mechanism(mixed)
    @test length(exchange_group(mm)) == 2
    @test all(takes_up_I, exchange_group(mm))
    @test ER.Mechanism(@enzyme_mechanism begin
        substrates: A
        products: P
        regulators: I, J
        steps: begin
            E + A <--> E(A)
            E(A) <--> E(P)
            E(P) <--> F + P
            F <--> E
            (E + I ⇌ E(I), F + I ⇌ F(I))
            (E(J) + I <--> F(I) + J, E(J) + I <--> E(I) + J)
        end
    end) == mm
    @test occursin("k_EJinh_Iinh_to_EIinh_Jinh", rate_equation_string(mixed))
end

@testset "each reaction appears once in a mechanism" begin
    A, B, P = ER.Substrate(:A), ER.Substrate(:B), ER.Product(:P)
    E, EA, EstarA = _testhelper_sp([]), _testhelper_sp([A]), _testhelper_sp([A], :Estar)
    EAB, EP = _testhelper_sp([A, B]), _testhelper_sp([P])
    rxn = @enzyme_reaction(begin
        substrates: A[C]
        products: P[C]
    end)
    # A binds E into two different forms from two groups: two reactions, whose
    # names differ in the bound form (k_E_A_to_EA, k_E_A_to_EstarA).
    m = ER.Mechanism(rxn, [
        [ER.Step(E, EA, [A], ER.Metabolite[], false)],
        [ER.Step(E, EstarA, [A], ER.Metabolite[], false)],
        [ER.Step(EA, EP, ER.Metabolite[], ER.Metabolite[], false)],
        [ER.Step(E, EP, [P], ER.Metabolite[], false)]])
    names = [ER.name(p, m) for p in ER._enumerate_parameters_full(m)]
    @test allunique(names)
    @test issubset([:k_E_A_to_EA, :k_E_A_to_EstarA], names)
    # The same binding in two groups, once written as its release: both would be
    # k_E_A_to_EA.
    err = _testhelper_thrown() do
        ER.Mechanism(rxn, [
            [ER.Step(E, EA, [A], ER.Metabolite[], false)],
            [ER.Step(EA, E, ER.Metabolite[], [A], false)],
            [ER.Step(EA, EP, ER.Metabolite[], ER.Metabolite[], false)],
            [ER.Step(E, EP, [P], ER.Metabolite[], false)]])
    end
    @test err isa ErrorException
    @test occursin("both hold the reaction E_A ⇌ EA", err.msg)
    # The same isomerization in two groups: both would be k_EA_to_EP.
    err = _testhelper_thrown() do
        ER.Mechanism(rxn, [
            [ER.Step(E, EA, [A], ER.Metabolite[], false)],
            [ER.Step(EA, EP, ER.Metabolite[], ER.Metabolite[], false)],
            [ER.Step(EP, EA, ER.Metabolite[], ER.Metabolite[], false)],
            [ER.Step(E, EP, [P], ER.Metabolite[], false)]])
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
    # E with NAD and P bound and E with NADP bound both render :ENADP, so every
    # constant named after one form would also name the other.
    NAD, P, NADP = ER.Substrate(:NAD), ER.Substrate(:P), ER.Product(:NADP)
    E, ENAD = _testhelper_sp([]), _testhelper_sp([NAD])
    E_NAD_P, E_NADP = _testhelper_sp([NAD, P]), _testhelper_sp([NADP])
    @test ER.name(E_NAD_P) == ER.name(E_NADP) == :ENADP
    rxn = @enzyme_reaction(begin
        substrates: NAD[C], P[N]
        products: NADP[CN]
    end)
    err = _testhelper_thrown() do
        ER.Mechanism(rxn, [
            [ER.Step(E, ENAD, [NAD], ER.Metabolite[], true)],
            [ER.Step(ENAD, E_NAD_P, [P], ER.Metabolite[], true)],
            [ER.Step(E_NAD_P, E_NADP, ER.Metabolite[], ER.Metabolite[], false)],
            [ER.Step(E, E_NADP, [NADP], ER.Metabolite[], true)]])
    end
    @test err isa ErrorException
    @test occursin("E(NAD, P)", err.msg)
    @test occursin("E(NADP)", err.msg)
    @test occursin("both render the name ENADP; ", err.msg)
    @test occursin("rename a metabolite", err.msg)
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
    src_dir = joinpath(dirname(@__DIR__), "src")
    for f in readdir(src_dir; join=true)
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

@testset "OnlyA Haldane validator" begin
    # Uni-uni S -> P. Tags: (S binding, chemical step, P binding).
    function uni(s_tag, cat_tag, p_tag)
        m = @allosteric_mechanism begin
            substrates: S
            products:   P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + S ⇌ E(S)      :: EqualAI
                E(S) <--> E(P)    :: EqualAI
                E(P) ⇌ E + P      :: EqualAI
            end
        end
        am = ER.AllostericMechanism(m)
        tags = copy(ER.cat_allo_states(am))
        # find each group by its representative step
        for (g, grp) in enumerate(ER.steps(am))
            bm = ER.bound_metabolite(grp[1])
            tags[g] = bm === nothing ? cat_tag :
                      ER.name(bm) === :S ? s_tag : p_tag
        end
        (ER.reaction(am), ER.steps(am), tags)
    end

    # no :OnlyA anywhere -> valid
    @test ER._onlya_haldane_violation(uni(:EqualAI, :EqualAI, :EqualAI)...) === nothing
    # :OnlyA on the substrate only, catalysis :EqualAI -> VIOLATION
    @test ER._onlya_haldane_violation(uni(:OnlyA, :EqualAI, :EqualAI)...) isa String
    # :OnlyA on the product only, catalysis :EqualAI -> VIOLATION
    @test ER._onlya_haldane_violation(uni(:EqualAI, :EqualAI, :OnlyA)...) isa String
    # :OnlyA on the substrate, catalysis :OnlyA -> the k_I = 0 escape -> valid
    @test ER._onlya_haldane_violation(uni(:OnlyA, :OnlyA, :EqualAI)...) === nothing
    # balanced: :OnlyA on both sides, catalysis :EqualAI -> valid
    @test ER._onlya_haldane_violation(uni(:OnlyA, :EqualAI, :OnlyA)...) === nothing
    # V-system: :OnlyA chemical step only -> valid
    @test ER._onlya_haldane_violation(uni(:EqualAI, :OnlyA, :EqualAI)...) === nothing
    # :NonequalAI catalysis is also a finite-nonzero assertion -> same verdict
    @test ER._onlya_haldane_violation(uni(:OnlyA, :NonequalAI, :EqualAI)...) isa String

    # Whether a binding is rapid-equilibrium or steady-state does not change
    # which affinities diverge, so the balanced both-:OnlyA verdict must not
    # depend on it: an RE binding's 1/Kd and an SS binding's kon/koff enter the
    # cycle product with the same exponent. Both orientations of the mix are equally
    # valid mechanisms and must both pass.
    mixed_uni_re_ss = @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)      :: EqualAI
            E(S) <--> E(P)    :: EqualAI
            E(P) <--> E + P   :: EqualAI
        end
    end
    mixed_uni_ss_re = @allosteric_mechanism begin
        substrates: S
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S <--> E(S)   :: EqualAI
            E(S) <--> E(P)    :: EqualAI
            E(P) ⇌ E + P      :: EqualAI
        end
    end
    # every binding :OnlyA, catalysis :EqualAI -> balanced -> valid
    function both_bindings_onlya(m)
        am = ER.AllostericMechanism(m)
        tags = [ER.bound_metabolite(g[1]) === nothing ? :EqualAI : :OnlyA
                for g in ER.steps(am)]
        ER._onlya_haldane_violation(ER.reaction(am), ER.steps(am), tags)
    end
    @test both_bindings_onlya(mixed_uni_re_ss) === nothing
    @test both_bindings_onlya(mixed_uni_ss_re) === nothing

    # A random-order binding square contributes a Wegscheider row: rhs = 0 and
    # no k columns. There is no k to zero out, so balance is the only escape and
    # the :OnlyA-chemical-step escape deliberately does not reach it.
    @testset "random-order binding square" begin
        biuni = @allosteric_mechanism begin
            substrates: A, B
            products:   P
            catalytic_multiplicity: 2
            catalytic_steps: begin
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
        function tags(onlya_keys...)
            want = Set{Tuple{Symbol, Union{Symbol, Nothing}}}(onlya_keys)
            map(ER.steps(bu_am)) do grp
                bm = ER.bound_metabolite(grp[1])
                key = (ER.name(ER.from_species(grp[1])),
                       bm === nothing ? nothing : ER.name(bm))
                key in want ? :OnlyA : :EqualAI
            end
        end
        verdict(t) = ER._onlya_haldane_violation(ER.reaction(bu_am),
                                                 ER.steps(bu_am), t)

        # the square alone trips nothing: no :OnlyA anywhere -> valid
        @test verdict(tags()) === nothing
        # one square edge :OnlyA -> VIOLATION (unbalanced in the square and
        # in the Haldane row)
        @test verdict(tags((:E, :A))) isa String
        # ...and the constructor rejects that tagging with the full message
        lone_onlya = ErrorException(
            "AllostericMechanism: an :OnlyA binding (K_EA_to_E_A) leaves a " *
            "thermodynamic (Haldane/Wegscheider) cycle unsatisfiable: the inactive " *
            "conformation cannot close that cycle at finite nonzero affinity. Tag " *
            "the cycle's chemical step :OnlyA, or tag an opposing binding :OnlyA so " *
            "the affinities diverge together.")
        @test_throws lone_onlya @allosteric_mechanism begin
            substrates: A, B
            products:   P
            catalytic_multiplicity: 2
            catalytic_steps: begin
                E + A ⇌ E(A)          :: OnlyA
                E + B ⇌ E(B)          :: EqualAI
                E(A) + B ⇌ E(A, B)    :: EqualAI
                E(B) + A ⇌ E(A, B)    :: EqualAI
                E(A, B) <--> E(P)     :: EqualAI
                E(P) ⇌ E + P          :: EqualAI
            end
        end
        # both A-side bindings :OnlyA balances the square, but the Haldane row
        # is still unbalanced against P -> VIOLATION
        @test verdict(tags((:E, :A), (:EB, :A))) isa String
        # ...and :OnlyA on P balances the Haldane row too -> valid
        @test verdict(tags((:E, :A), (:EB, :A), (:E, :P))) === nothing
        # :OnlyA on the chemical step drops that group, killing the Haldane row
        # and leaving only the square, whose lone :OnlyA edge is unbalanced ->
        # VIOLATION. A keep filter that dropped *every* :OnlyA group rather than
        # only the is_iso ones would drop the square edge as well and wrongly
        # report valid.
        @test verdict(tags((:EAB, nothing), (:E, :A))) isa String
    end
end

@testset ":OnlyA guard rejects a multi-cycle ter-substrate inconsistency" begin
    # A full random-order ter binding cube with seven :OnlyA edges. Every single
    # constraint row carries BOTH signs on its :OnlyA eps-exponents, so the
    # per-row sign test sees no violation — but the coupled system has no
    # strictly-positive solution: rows 2, 4 and 5 combine to force the
    # eps-exponent of K_EAB_to_EA_B to zero, i.e. K_I = K_A, contradicting its :OnlyA
    # tag. The inactive cube circulates flux around the E(A)->E(A,B)<-E(B)->
    # E(B,C)<-E(C)->E(A,C)<-E(A) hexagon at equilibrium — perpetual motion.
    @test_throws ErrorException @allosteric_mechanism begin
        substrates: A, B, C
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
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

# The rejection testset above pins the guard's reject direction; this one pins its
# ACCEPT direction, which nothing else covers. The guard runs a cheap per-row sign
# test (sound but incomplete) and then an exact Stiemke eps-feasibility test. The
# two agree on every mechanism up to bi-bi and diverge only from ter-substrate up,
# so ter is the only place a regression in the exact test can show. Over-rejection
# was measured at 0 across 17,814 tag assignments — every feasible assignment is
# admitted. That is what makes this direction worth a test: an over-rejection does
# not fail loudly, it silently shrinks the enumeration search space by dropping
# valid mechanisms before they are ever fitted. The cube below is the load-bearing
# witness — its sign test finds nothing to flag, so the verdict genuinely comes
# from the Stiemke stage (M is 5x12, nullity 7, feasible = true).
@testset ":OnlyA guard admits feasible ter-substrate mechanisms" begin
    # A full random-order A/B/C binding lattice, every substrate binding :OnlyA.
    # The :OnlyA chemical step drops the catalytic Haldane row from the check
    # graph, leaving the lattice's five Wegscheider squares; each carries both
    # signs on its eps-exponents, so the coupled system is feasible (w = 1 solves
    # it) and the cube is admitted.
    cube = @allosteric_mechanism begin
        substrates: A, B, C
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
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

    # Ordered ter with an :OnlyA binding and an :OnlyA chemical step: the free
    # inactive k_I ratio absorbs the B affinity's divergence, so this is valid.
    ordered = @allosteric_mechanism begin
        substrates: A, B, C
        products:   P
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + A ⇌ E(A)             :: EqualAI
            E(A) + B ⇌ E(A, B)       :: OnlyA
            E(A, B) + C ⇌ E(A, B, C) :: EqualAI
            E(A, B, C) <--> E(P)     :: OnlyA
            E + P ⇌ E(P)             :: EqualAI
        end
    end
    @test ordered isa ER.AllostericEnzymeMechanism
end

@testset ":OnlyA guard admits a ping-pong whose fused releases are :OnlyA" begin
    # E(A) → E(; residual = A - P) + P runs the first half-reaction and releases P in
    # one step, E(B; residual = A - P) → E + Q the second and releases Q: chemistry
    # (`_is_chemistry`), though each gives off one metabolite and takes up none.
    # Tagged :OnlyA, both leave the check graph as :OnlyA isomerizations do, and the
    # inactive conformation runs no chemistry. No plain binding is :OnlyA, so the
    # guard has no affinity to drive to zero and returns before it builds a cycle
    # row. Read as :OnlyA bindings, P and Q would stand on the product side of the
    # one Haldane cycle, one sign, and the mechanism would be refused.
    pingpong = @allosteric_mechanism begin
        substrates: A, B ; products: P, Q ; catalytic_multiplicity: 1
        catalytic_steps: begin
            E + A <--> E(A)                                       :: EqualAI
            E(A) <--> E(; residual = A - P) + P                   :: OnlyA
            E(; residual = A - P) + B <--> E(B; residual = A - P) :: EqualAI
            E(B; residual = A - P) <--> E + Q                     :: OnlyA
        end
    end
    am = ER.AllostericMechanism(pingpong)
    @test am isa ER.AllostericMechanism
    @test ER._onlya_haldane_violation(ER.reaction(am), ER.steps(am),
                                      ER.cat_allo_states(am)) === nothing
end

@testset "rational nullspace + Stiemke feasibility helpers" begin
    R = Rational{BigInt}

    @testset "_rational_nullspace" begin
        # Zero map: every x solves M * x = 0, so the nullspace is the whole
        # 3-dimensional space.
        M = zeros(R, 2, 3)
        N = ER._rational_nullspace(M)
        @test size(N, 2) == 3

        # Full-rank square matrix: only x = 0 solves M * x = 0.
        M = R[1 2; 3 4]
        N = ER._rational_nullspace(M)
        @test size(N, 2) == 0

        # Rank-1 2x2 (row 2 is twice row 1): both rows reduce to x + y = 0,
        # spanned by (1, -1). n = 2, rank = 1.
        M = R[1 1; 2 2]
        N = ER._rational_nullspace(M)
        @test M * N == zeros(R, 2, size(N, 2))
        @test size(N, 2) == 2 - 1

        # More rows than columns, rank-deficient: every row is a multiple of
        # (1, 2), i.e. x + 2y = 0, spanned by (2, -1).
        M = R[1 2; 2 4; 3 6]
        N = ER._rational_nullspace(M)
        @test M * N == zeros(R, 3, size(N, 2))
        @test size(N, 2) == 2 - 1

        # More columns than rows: x + z = 0 and y + z = 0, spanned by
        # (-1, -1, 1).
        M = R[1 0 1; 0 1 1]
        N = ER._rational_nullspace(M)
        @test M * N == zeros(R, 2, size(N, 2))
        @test size(N, 2) == 3 - 2
    end

    @testset "_has_strict_positive_combination" begin
        # Single row [1]: y = 1 gives N * y = 1 > 0.
        @test ER._has_strict_positive_combination(reshape(R[1], 1, 1))
        # Rows [1] and [-1] demand y > 0 and -y > 0 at once -> infeasible.
        @test !ER._has_strict_positive_combination(reshape(R[1, -1], 2, 1))
        # An all-zero row encodes 0 > 0 and refutes the system outright,
        # regardless of the other row.
        @test !ER._has_strict_positive_combination(R[0 0; 1 1])
        # No columns: y has no entries, so every row's dot product with y is
        # 0, and each row still demands 0 > 0 -> infeasible.
        @test !ER._has_strict_positive_combination(zeros(R, 3, 0))
        # 2-D feasible cone: y1 + y2 > 0 and y1 - y2 > 0, e.g. y = (1, 0).
        @test ER._has_strict_positive_combination(R[1 1; 1 -1])
        # 2-D infeasible cone: y1 > 0 and y2 > 0 force y1 + y2 > 0, which
        # contradicts the third row's -y1 - y2 > 0.
        @test !ER._has_strict_positive_combination(R[1 0; 0 1; -1 -1])
    end

    @testset "_solve_dependent_set rejects a contradictory system" begin
        # The rows x = 0 and x = log(Keq) subtract to 0 = log(Keq): the reduced
        # [A rhs] pivots in the rhs column, which no parameter value satisfies.
        err = @test_throws ErrorException ER._solve_dependent_set(
            reshape(R[1, 1], 2, 1), R[0, 1], [:K_x], [(false, 0)])
        @test occursin("Thermodynamically contradictory mechanism", err.value.msg)
    end
end

@testset "the naming cache never changes a mechanism's identity" begin
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
    m2 = EnzymeRates.Mechanism(EnzymeRates.reaction(m1),
                               deepcopy(EnzymeRates.steps(m1)))
    am1 = EnzymeRates.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        allosteric_regulators: I::OnlyI
        catalytic_multiplicity: 2
        catalytic_steps: begin
            E + S ⇌ E(S)     :: EqualAI
            E(S) <--> E(P)   :: OnlyA
            E(P) ⇌ E + P     :: EqualAI
        end
    end)
    am2 = EnzymeRates.AllostericMechanism(EnzymeRates.reaction(am1),
                                          deepcopy(EnzymeRates.steps(am1)),
                                          copy(EnzymeRates.cat_allo_states(am1)),
                                          EnzymeRates.catalytic_multiplicity(am1),
                                          copy(EnzymeRates.regulatory_sites(am1)))
    # Naming every parameter fills the mechanism's naming cache.
    param_names(m::EnzymeRates.Mechanism) =
        [EnzymeRates.name(p, m) for p in EnzymeRates._enumerate_parameters_full(m)]
    param_names(m::EnzymeRates.AllostericMechanism) = [EnzymeRates.name(p, m)
        for state in (:A, :I)
        for p in [EnzymeRates._cat_params(m, state); EnzymeRates._kreg_params(m, state)]]
    # `rebuilt` is `m` constructed again from its fields, with an empty cache.
    function check_identity(m, rebuilt)
        shown = repr(m)
        names = param_names(m)
        @test repr(m) == shown
        @test m == rebuilt && hash(m) == hash(rebuilt)
        @test EnzymeRates.compile_mechanism(m) === EnzymeRates.compile_mechanism(rebuilt)
        @test names == param_names(rebuilt)
        io = IOBuffer(); serialize(io, m); seekstart(io)
        copied = deserialize(io)
        @test copied == m && hash(copied) == hash(m)
        @test names == param_names(copied)
    end
    check_identity(m1, m2)
    check_identity(am1, am2)
    # A species' stored name is rendered from its sorted bound list, so the order in
    # which the bound metabolites are given does not change it. Display leaves it out.
    A, B = EnzymeRates.Substrate(:A), EnzymeRates.Substrate(:B)
    sp = EnzymeRates.Species(EnzymeRates.Metabolite[A, B], :E)
    sp2 = EnzymeRates.Species(EnzymeRates.Metabolite[B, A], :E)
    @test sp2 == sp && hash(sp2) == hash(sp)
    @test EnzymeRates.name(sp2) === EnzymeRates.name(sp) === :EAB
    @test !occursin("EAB", repr(sp))
end
