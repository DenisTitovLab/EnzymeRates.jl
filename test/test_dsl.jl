@testset "DSL" begin
    @testset "@enzyme_mechanism (new grammar)" begin
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I

            steps: begin
                ([E, S] ⇌ [ES], [EP, S] ⇌ [EPS])
                [ES, I] ⇌ [ESI]
                [ES]   <--> [EP]
                [EP]   ⇌    [E, P]
            end
        end
        @test EnzymeRates.substrates(m) == (:S,)
        @test EnzymeRates.products(m) == (:P,)
        @test EnzymeRates.regulators(m) == (:I,)
        @test EnzymeRates.kinetic_group(m, 1) == EnzymeRates.kinetic_group(m, 2)
        @test EnzymeRates.kinetic_group(m, 3) != EnzymeRates.kinetic_group(m, 4)

        # Reject atom bracket syntax in substrates:
        @test_throws Exception eval(:(@enzyme_mechanism begin
            substrates: S[C]
            products:   P
            steps: begin
                [E, S] ⇌ [ES]
                [ES] <--> [EP]
                [EP] ⇌ [E, P]
            end
        end))

        # Reject allosteric-only syntax (site(:catalytic, N))
        @test_throws Exception eval(:(@enzyme_mechanism begin
            substrates: S
            products:   P
            site(:catalytic, 2): begin
                steps: begin
                    [E, S] ⇌ [ES]
                    [ES] <--> [EP]
                    [EP] ⇌ [E, P]
                end
            end
        end))
    end

    @testset "@allosteric_mechanism (parsing & validation)" begin
        m = @allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            allosteric_regulators: I::OnlyT

            site(:catalytic, 2): begin
                steps: begin
                    [E, F6P] ⇌ [E_F6P]    :: EqualRT
                    [E_F6P] <--> [E_F16BP] :: EqualRT
                    [E_F16BP] ⇌ [E, F16BP] :: EqualRT
                end
            end
        end
        @test m isa EnzymeRates.AllostericEnzymeMechanism
        @test EnzymeRates.allosteric_regulators(m) ⊇ ((:I, :OnlyT),)

        # Reject untagged catalytic step
        @test_throws Exception eval(:(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            site(:catalytic, 2): begin
                steps: begin
                    [E, F6P] ⇌ [E_F6P] :: EqualRT
                    [E_F6P] <--> [E_F16BP]
                    [E_F16BP] ⇌ [E, F16BP] :: EqualRT
                end
            end
        end))

        # Reject :OnlyT on a catalytic step (V-type allostery not supported)
        @test_throws Exception eval(:(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            site(:catalytic, 2): begin
                steps: begin
                    [E, F6P] ⇌ [E_F6P] :: EqualRT
                    [E_F6P] <--> [E_F16BP] :: OnlyT
                    [E_F16BP] ⇌ [E, F16BP] :: EqualRT
                end
            end
        end))

        # Reject untagged allosteric regulator
        @test_throws Exception eval(:(@allosteric_mechanism begin
            substrates: F6P
            products:   F16BP
            allosteric_regulators: I, J::OnlyT
            site(:catalytic, 2): begin
                steps: begin
                    [E, F6P] ⇌ [E_F6P] :: EqualRT
                    [E_F6P] <--> [E_F16BP] :: EqualRT
                    [E_F16BP] ⇌ [E, F16BP] :: EqualRT
                end
            end
        end))

        # Reject parenthesized step group without ::AlloState
        @test_throws Exception eval(:(@allosteric_mechanism begin
            substrates: S, A
            products:   P
            site(:catalytic, 2): begin
                steps: begin
                    ([E, S] ⇌ [ES], [E_A, S] ⇌ [ES_A])
                    [ES_A] <--> [EP]   :: EqualRT
                    [EP] ⇌ [E, P]      :: EqualRT
                end
            end
        end))
    end

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
            substrates: S
            products:   P
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [E, P]
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.n_steps(m) == 2
        @test EnzymeRates.n_states(m) == 2
        @test Set(EnzymeRates.enzyme_forms(m)) == Set([:E, :ES])
        @test Set(metabolites(m)) == Set([:S, :P])

        # Numeric check: same as Uni-Uni spot check
        Keq = 3.2 * 2.5 / (0.8 * 1.1)
        params = (k1f=3.2, k2f=2.5, k2r=1.1, Keq=Keq, E_total=1.0)
        concs = (S=0.7, P=0.3)
        @test rate_equation(m, concs, params) ≈ 0.9091 atol=0.001

        # Multi-step mechanism
        m2 = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
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
            substrates: S
            products:   P
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
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            regulators: I
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [E, P]
                [E, I] <--> [EI]
            end
        end
        @test m isa EnzymeMechanism
    end

    @testset "No-atom species" begin
        # All metabolites without atoms
        m = @enzyme_mechanism begin
            substrates: S
            products:   P
            steps: begin
                [E, S] <--> [ES]
                [ES] <--> [E, P]
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.n_steps(m) == 2
    end

    @testset "Constraint DSL parsing" begin
        # Bi-bi random with two K_A binding steps in shared kinetic group
        m = @enzyme_mechanism begin
            substrates: A, B
            products:   P, Q
            steps: begin
                ([E, A] ⇌ [EA], [EB, A] ⇌ [EAB])
                [E, B] ⇌ [EB]
                [EA, B] ⇌ [EAB]
                [EAB] <--> [EPQ]
                [EPQ] ⇌ [EQ, P]
                [EQ] ⇌ [E, Q]
            end
        end
        @test m isa EnzymeMechanism
        @test EnzymeRates.kinetic_group(m, 1) == EnzymeRates.kinetic_group(m, 2)
        @test EnzymeRates.kinetic_group(m, 3) != EnzymeRates.kinetic_group(m, 1)
    end
end
