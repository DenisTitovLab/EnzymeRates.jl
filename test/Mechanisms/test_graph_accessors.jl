# ABOUTME: Tests of what a mechanism's stored graph holds (forms, steps, kinetic groups,
# ABOUTME: allosteric parts) and of step text.

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
