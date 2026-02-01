@testset "DSL" begin
    @testset "@enzyme_reaction" begin
        spec = @enzyme_reaction begin
            substrates: S(C=1)
            products:   P(C=1)
        end
        @test spec isa EnzymeReaction
        @test substrates(spec) == ((:S, ((:C, 1),)),)
        @test products(spec) == ((:P, ((:C, 1),)),)
        @test regulators(spec) == ()

        spec2 = @enzyme_reaction begin
            substrates: S(C=6, H=12, O=6), ATP(C=10, H=16, N=5, O=13, P=3)
            products:   G6P(C=6, H=13, O=9, P=1), ADP(C=10, H=15, N=5, O=10, P=2)
            regulators: I(C=5, H=8, N=2)
        end
        @test length(substrates(spec2)) == 2
        @test length(products(spec2)) == 2
        @test length(regulators(spec2)) == 1
        @test regulators(spec2)[1][1] == :I
    end

    @testset "@mechanism" begin
        m = @mechanism begin
            species: begin
                substrates: S(C=1)
                products:   P(C=1)
                enzymes:    E(), ES(C=1)
            end
            steps: begin
                [E, S] --> [ES]
                [ES] --> [E, P]
            end
        end
        @test m isa EnzymeMechanism
        @test n_steps(m) == 2
        @test n_states(m) == 2
        @test Set(e[1] for e in enzyme_forms(m)) == Set([:E, :ES])
        @test Set(m[1] for m in metabolites(m)) == Set([:S, :P])

        # Numeric check: same as Uni-Uni spot check
        params = (k1f=3.2, k1r=0.8, k2f=2.5, k2r=1.1)
        concs = (S=0.7, P=0.3)
        @test rate_equation(m, params, concs) ≈ 0.9091 atol=0.001

        # Multi-step mechanism
        m2 = @mechanism begin
            species: begin
                substrates: A(C=2, N=1), B(C=3)
                products:   P(C=2), Q(C=3, N=1)
                enzymes:    E(), EA(C=2, N=1), FP(C=2, N=1), F(N=1), FB(C=3, N=1), EQ(C=3, N=1)
            end
            steps: begin
                [E, A] --> [EA]
                [EA] --> [FP]
                [FP] --> [F, P]
                [F, B] --> [FB]
                [FB] --> [EQ]
                [EQ] --> [E, Q]
            end
        end
        @test n_states(m2) == 6
    end

    @testset "Elementary steps" begin
        @test_throws ErrorException @mechanism begin
            species: begin
                substrates: S(C=1)
                products:   P(C=1)
                enzymes:    E(), ESP(C=1)
            end
            steps: begin
                [E, S, P] --> [ESP]
            end
        end

        spec = @enzyme_reaction begin
            substrates: S(C=1)
            products:   P(C=1)
            regulators: I(C=1)
        end
        @test spec isa EnzymeReaction

        species = (
            ( (:S, ((:C, 1),)), ),           # substrates
            ( (:P, ((:C, 1),)), ),           # products
            ( (:I, ((:C, 1),)), ),           # regulators
            ( (:E, ()), (:ES, ((:C, 1),)), (:EI, ((:C, 1),)) ),  # enzymes
        )
        rxns = (
            ((:E, :S), (:ES,)),
            ((:ES,), (:E, :P)),
            ((:E, :I), (:EI,)),
        )
        @test_throws ErrorException EnzymeMechanism(species, rxns)
    end
end
