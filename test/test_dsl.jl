@testset "DSL" begin
    @testset "@enzyme_reaction" begin
        spec = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
        end
        @test spec isa EnzymeReaction
        @test EnzymeRates.substrates(spec) == ((:S, ((:C, 1),)),)
        @test EnzymeRates.products(spec) == ((:P, ((:C, 1),)),)
        @test EnzymeRates.regulators(spec) == ()

        spec2 = @enzyme_reaction begin
            substrates: S[C6H12O6], ATP[C10H16N5O13P3]
            products:   G6P[C6H13O9P], ADP[C10H15N5O10P2]
            regulators: I[C5H8N2]
        end
        @test length(EnzymeRates.substrates(spec2)) == 2
        @test length(EnzymeRates.products(spec2)) == 2
        @test length(EnzymeRates.regulators(spec2)) == 1
        @test EnzymeRates.regulators(spec2)[1][1] == :I
    end

    @testset "@enzyme_mechanism" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products:   P[C]
                enzymes:    E, ES[C]
            end
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [E, P]
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.n_steps(m) == 2
        @test EnzymeRates.n_states(m) == 2
        @test Set(e[1] for e in EnzymeRates.enzyme_forms(m)) == Set([:E, :ES])
        @test Set(m[1] for m in metabolites(m)) == Set([:S, :P])

        # Numeric check: same as Uni-Uni spot check
        Keq = 3.2 * 2.5 / (0.8 * 1.1)
        params = (k1f=3.2, k2f=2.5, k2r=1.1, Keq=Keq, E_total=1.0)
        concs = (S=0.7, P=0.3)
        @test rate_equation(m, params, concs) ≈ 0.9091 atol=0.001

        # Multi-step mechanism
        m2 = @enzyme_mechanism begin
            species: begin
                substrates: A[C2N], B[C3]
                products:   P[C2], Q[C3N]
                enzymes:    E, EA[C2N], FP[C2N], F[N], FB[C3N], EQ[C3N]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [FP]
                [FP] <--> [F, P]
                [F, B] <--> [FB]
                [FB] <--> [EQ]
                [EQ] <--> [E, Q]
            end
        end
        @test EnzymeRates.n_states(m2) == 6
    end

    @testset "Elementary steps" begin
        @test_throws ErrorException @enzyme_mechanism begin
            species: begin
                substrates: S[C]
                products:   P[C]
                enzymes:    E, ESP[C]
            end
            steps: begin
                [E, S, P] <--> [ESP]
            end
        end

        spec = @enzyme_reaction begin
            substrates: S[C]
            products:   P[C]
            regulators: I[C]
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
        @test_throws ErrorException EnzymeMechanism(species, rxns, (false, false, false))
    end

    @testset "No-atom species" begin
        # All metabolites without atoms — should skip conservation checks
        m = @enzyme_mechanism begin
            species: begin
                substrates: S
                products:   P
                enzymes:    E, ES
            end
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [E, P]
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.n_steps(m) == 2
    end

    @testset "Mixed atoms error" begin
        # Some metabolites with atoms, some without — should error
        species = (
            ( (:S, ((:C, 1),)), ),   # substrates — has atoms
            ( (:P, ()), ),           # products — no atoms
            (),                      # regulators
            ( (:E, ()), (:ES, ((:C, 1),)) ),  # enzymes
        )
        rxns = (
            ((:E, :S), (:ES,)),
            ((:ES,), (:E, :P)),
        )
        @test_throws ErrorException EnzymeMechanism(species, rxns, (false, false))
    end
end
