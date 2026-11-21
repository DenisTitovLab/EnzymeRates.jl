# ABOUTME: Tests of the allosteric conformation semantics: the steps a conformation
# ABOUTME: excludes and runs (`_excluded_steps`, `exclusion_components`).

const ER = EnzymeRates

@testset "steps_in_active_conformation and steps_in_inactive_conformation" begin
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
    # also the binding of P to it.
    @test M.steps_in_active_conformation(am) == BitVector([0, 1, 1, 0, 1])
    # The inactive conformation drops the :OnlyA chemistry and keeps every form.
    @test M.steps_in_inactive_conformation(am) == BitVector([1, 1, 1, 1, 0])
end

@testset "excluded steps of a plain and an allosteric mechanism" begin
    D = ER.Derivation
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
    m = _testhelper_plain_of(am)
    # A plain mechanism excludes no step. The steps are stored as in the testset above:
    # the active conformation excludes the :OnlyI binding of I, the inactive one the
    # :OnlyA chemistry.
    @test D._excluded_steps(m, 1) == falses(5)
    @test D._excluded_steps(m, 2) == falses(5)
    @test D._excluded_steps(am, 1) == BitVector([1, 0, 0, 0, 0])
    @test D._excluded_steps(am, 2) == BitVector([0, 0, 0, 0, 1])
end

@testset "the excluded steps and exclusion_components: an :OnlyA dead-end binding" begin
    M = ER.Mechanisms
    D = ER.Derivation
    am = @enzyme_mechanism begin
        substrates: S
        products: P
        oligomeric_state: 2
        steps: begin
            E + S ⇌ ES                    :: EqualAI
            ES + S ⇌ ESSmod    :: OnlyA
            ES <--> EP                  :: EqualAI
            EP ⇌ E + P                    :: EqualAI
        end
    end
    m = am.catalytic
    # Canonical form numbers the forms E, EP, ES, ESSmod and stores the steps
    # E + P → EP, E + S → ES, ES → EP, ES + S → ESSmod. The inactive
    # conformation excludes the :OnlyA binding, which leaves ESSmod alone.
    @test M.step_text(m, 4) == "ES + S ⇌ ESSmod"
    @test D._excluded_steps(am, 1) == [false, false, false, false]
    @test D._excluded_steps(am, 2) == [false, false, false, true]
    @test D.exclusion_components(am, :active) == [1, 1, 1, 1]
    @test D.exclusion_components(am, :inactive) == [1, 1, 1, 4]
    # The rule by tag, one tag per group: an :OnlyA step leaves the inactive conformation,
    # an :OnlyI step the active one; an :EqualAI or :NonequalAI step stays in both.
    tags = [:OnlyA, :OnlyI, :EqualAI, :NonequalAI]
    every_tag = ER.Mechanism(m, ER.Mechanisms.Allostery(tags, 2, M.RegulatorySite[]))
    step_tags = [tags[e.kinetic_group] for e in m.steps]
    @test step_tags[D._excluded_steps(every_tag, 1)] == [:OnlyI]
    @test step_tags[D._excluded_steps(every_tag, 2)] == [:OnlyA]
    @test_throws ArgumentError D.exclusion_components(am, :both)
end
