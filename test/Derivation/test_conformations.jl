# ABOUTME: Tests of the allosteric conformation semantics: the steps a conformation
# ABOUTME: excludes and keeps (`excludes`, `exclusion_components`, the stored step masks).

const ER = EnzymeRates

@testset "active_steps and inactive_steps" begin
    M = ER.Mechanisms
    D = ER.Derivation
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)          :: EqualAI
            E(S) <--> E(P)        :: OnlyA
            E(P) ⇌ E + P          :: EqualAI
            E + I ⇌ E(I)          :: OnlyI
            E(I) + P ⇌ E(I, P)    :: EqualAI
        end
    end)
    m = M.catalytic(am)
    # Canonical form numbers the forms E, E(I), E(I, P), E(P), E(S) and stores the steps
    # E + I → E(I), E + P → E(P), E + S → E(S), E(I) + P → E(I, P), E(S) → E(P).
    @test [(e.from, e.to) for e in m.steps] == [(1, 2), (1, 4), (1, 5), (2, 3), (5, 4)]
    # The active conformation drops the :OnlyI binding of I, which strands E(I) and so
    # also the binding of P to it.
    @test am.active_steps == BitVector([0, 1, 1, 0, 1])
    # The inactive conformation drops the :OnlyA chemistry and keeps every form.
    @test am.inactive_steps == BitVector([1, 1, 1, 1, 0])
end

@testset "conformation masks of a plain and an allosteric mechanism" begin
    D = ER.Derivation
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: I
        steps: begin
            E + S ⇌ E(S)          :: EqualAI
            E(S) <--> E(P)        :: OnlyA
            E(P) ⇌ E + P          :: EqualAI
            E + I ⇌ E(I)          :: OnlyI
            E(I) + P ⇌ E(I, P)    :: EqualAI
        end
    end)
    m = am.catalytic
    # A plain mechanism has only the active conformation and excludes no step.
    @test D._conformation_steps(m, 1) == trues(5)
    @test D._conformation_steps(m, 2) == falses(5)
    @test D._excluded_steps(m, 1) == falses(5)
    @test D._excluded_steps(m, 2) == falses(5)
    @test D._conformation_steps(am, 1) == am.active_steps
    @test D._conformation_steps(am, 2) == am.inactive_steps
    @test D._excluded_steps(am, 1) == am.excluded_active
    @test D._excluded_steps(am, 2) == am.excluded_inactive
    @test_throws BoundsError D._conformation_steps(am, 3)
end

@testset "excludes and exclusion_components: an :OnlyA dead-end binding" begin
    M = ER.Mechanisms
    D = ER.Derivation
    am = ER.AllostericMechanism(@allosteric_mechanism begin
        substrates: S
        products: P
        modifiers: S
        catalytic_multiplicity: 2
        steps: begin
            E + S ⇌ E(S)                    :: EqualAI
            E(S) + S::Mod ⇌ E(S, S::Mod)    :: OnlyA
            E(S) <--> E(P)                  :: EqualAI
            E(P) ⇌ E + P                    :: EqualAI
        end
    end)
    m = M.catalytic(am)
    # Canonical form numbers the forms E, E(P), E(S), E(S, S::Mod) and stores the steps
    # E + P → E(P), E + S → E(S), E(S) → E(P), E(S) + S::Mod → E(S, S::Mod). The inactive
    # conformation excludes the :OnlyA binding, which leaves E(S, S::Mod) alone.
    @test M.step_text(m, 4) == "E(S) + S::Mod ⇌ E(S, S::Mod)"
    @test am.excluded_active == [false, false, false, false]
    @test am.excluded_inactive == [false, false, false, true]
    @test D.exclusion_components(am, :active) == [1, 1, 1, 1]
    @test D.exclusion_components(am, :inactive) == [1, 1, 1, 4]
    # The rule by tag: an :OnlyA step leaves the inactive conformation, an :OnlyI step the
    # active one; an :EqualAI or :NonequalAI step stays in both.
    @test [M.excludes(tag, c) for tag in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
           for c in (:active, :inactive)] ==
          [false, true, true, false, false, false, false, false]
    @test_throws ArgumentError M.excludes(:OnlyA, :both)
    @test_throws ArgumentError D.exclusion_components(am, :both)
end
