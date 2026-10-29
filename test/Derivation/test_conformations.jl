# ABOUTME: Tests of the allosteric conformation semantics: the steps a conformation
# ABOUTME: excludes and keeps (`excludes`, `exclusion_components`, `conformation_steps`).

const ER = EnzymeRates

@testset "conformation_steps" begin
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
    @test D.conformation_steps(am, :active) == BitVector([0, 1, 1, 0, 1])
    # The inactive conformation drops the :OnlyA chemistry and keeps every form.
    @test D.conformation_steps(am, :inactive) == BitVector([1, 1, 1, 1, 0])
    # A plain mechanism has only the active conformation.
    @test D.conformation_steps(m, :active) == trues(5)
    @test D.conformation_steps(m, :inactive) == falses(5)
    @test_throws ArgumentError D.conformation_steps(am, :I)
    @test_throws ArgumentError D.conformation_steps(m, :A)
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
    @test [D.excludes(am, k, :active) for k in eachindex(m.steps)] ==
          [false, false, false, false]
    @test [D.excludes(am, k, :inactive) for k in eachindex(m.steps)] ==
          [false, false, false, true]
    @test D.exclusion_components(am, :active) == [1, 1, 1, 1]
    @test D.exclusion_components(am, :inactive) == [1, 1, 1, 4]
    # A plain mechanism excludes no step from either conformation.
    @test [D.excludes(m, k, c) for k in eachindex(m.steps) for c in (:active, :inactive)] ==
          falses(8)
    # The rule by tag: an :OnlyA step leaves the inactive conformation, an :OnlyI step the
    # active one; an :EqualAI or :NonequalAI step stays in both.
    @test [D.excludes(tag, c) for tag in (:OnlyA, :OnlyI, :EqualAI, :NonequalAI)
           for c in (:active, :inactive)] ==
          [false, true, true, false, false, false, false, false]
    @test_throws ArgumentError D.excludes(:OnlyA, :both)
    @test_throws ArgumentError D.excludes(am, 1, :both)
    @test_throws ArgumentError D.excludes(m, 1, :both)
    @test_throws ArgumentError D.exclusion_components(am, :both)
end
