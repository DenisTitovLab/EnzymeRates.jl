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
            regulators: I
        end
        @test length(EnzymeRates.substrates(spec2)) == 2
        @test length(EnzymeRates.products(spec2)) == 2
        @test length(EnzymeRates.regulators(spec2)) == 1
        @test EnzymeRates.regulators(spec2)[1] == :I
    end

    @testset "multi-atom metabolites" begin
        rxn = @enzyme_reaction begin
            substrates: A[C2H3], B[N,P]
            products: P[C2,N], Q[H3,P]
        end
        subs = EnzymeRates.substrates(rxn)
        @test subs[1] == (:A, ((:C, 2), (:H, 3)))
        @test subs[2] == (:B, ((:N, 1), (:P, 1)))
        prods = EnzymeRates.products(rxn)
        @test prods[1] == (:P, ((:C, 2), (:N, 1)))
        @test prods[2] == (:Q, ((:H, 3), (:P, 1)))
    end

    @testset "@enzyme_reaction regulator roles" begin
        spec_roles = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            dead_end_inhibitors: I
            allosteric_regulators: A
            regulators: R
        end
        @test spec_roles isa EnzymeReaction
        @test Set(EnzymeRates.regulators(spec_roles)) == Set([:I, :A, :R])
        roles = EnzymeRates.regulator_roles(spec_roles)
        @test length(roles) == 3
        role_dict = Dict(r[1] => r[2] for r in roles)
        @test role_dict[:I] == :dead_end
        @test role_dict[:A] == :allosteric
        @test role_dict[:R] == :unknown

        # Backward compatibility: plain regulators
        spec_plain = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            regulators: R1, R2
        end
        roles_plain = EnzymeRates.regulator_roles(spec_plain)
        @test all(r[2] == :unknown for r in roles_plain)
        @test Set(r[1] for r in roles_plain) == Set([:R1, :R2])
    end

    @testset "@enzyme_reaction with oligomeric_state" begin
        rxn = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
            oligomeric_state: 4
        end
        @test EnzymeRates.oligomeric_state(rxn) == 4

        # Without oligomeric_state defaults to 1
        rxn2 = @enzyme_reaction begin
            substrates: S[C]
            products: P[C]
        end
        @test EnzymeRates.oligomeric_state(rxn2) == 1
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
        @test Set(metabolites(m)) == Set([:S, :P])

        # Numeric check: same as Uni-Uni spot check
        Keq = 3.2 * 2.5 / (0.8 * 1.1)
        params = (k1f=3.2, k2f=2.5, k2r=1.1, Keq=Keq, E_total=1.0)
        concs = (S=0.7, P=0.3)
        @test rate_equation(m, concs, params) ≈ 0.9091 atol=0.001

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
            regulators: I
        end
        @test spec isa EnzymeReaction

        # Dead-end inhibitor: valid mechanism (competitive inhibition)
        species = (
            ( (:S, ((:C, 1),)), ),           # substrates
            ( (:P, ((:C, 1),)), ),           # products
            (:I,),                           # regulators
            ( (:E, ()), (:ES, ((:C, 1),)), (:EI, ((:C, 1),)) ),  # enzymes
        )
        rxns = (
            ((:E, :S), (:ES,)),
            ((:ES,), (:E, :P)),
            ((:E, :I), (:EI,)),
        )
        @test EnzymeMechanism(species, rxns, (false, false, false)) isa EnzymeMechanism
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

    @testset "Constraint DSL parsing" begin
        # Simple K constraint via DSL
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] ⇌ [EA]
                [EA] <--> [EP]
                [EP] ⇌ [E, P]
            end
            constraints: begin
                K3 = K1
            end
        end
        @test m isa EnzymeMechanism
        pc = EnzymeRates.param_constraints(m)
        @test length(pc) == 1
        @test pc[1][1] == :K3  # target
        @test pc[1][2] == 1    # coeff
        @test pc[1][3] == ((:K1, 1),)  # factors

        # k constraint via DSL
        m2 = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
            constraints: begin
                k3r = k1r
            end
        end
        pc2 = EnzymeRates.param_constraints(m2)
        @test length(pc2) == 1
        @test pc2[1][1] == :k3r

        # Coefficient constraint
        m3 = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
            constraints: begin
                k3r = 2 * k1r
            end
        end
        pc3 = EnzymeRates.param_constraints(m3)
        @test pc3[1][2] == 2  # coeff = 2

        # Division constraint
        m4 = @enzyme_mechanism begin
            species: begin
                substrates: A[C]
                products:   P[C]
                enzymes:    E, EA[C], EP[C]
            end
            steps: begin
                [E, A] <--> [EA]
                [EA] <--> [EP]
                [EP] <--> [E, P]
            end
            constraints: begin
                k3r = k1f * k2f / k2r
            end
        end
        pc4 = EnzymeRates.param_constraints(m4)
        @test pc4[1][1] == :k3r
        @test Set(
            (sym, exp) for (sym, exp) in pc4[1][3]
        ) == Set([(:k1f, 1), (:k2f, 1), (:k2r, -1)])
    end

    @testset "No constraints backward compat" begin
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
        @test EnzymeRates.param_constraints(m) == ()
    end
end
