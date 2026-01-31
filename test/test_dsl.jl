@testset "DSL" begin
    @testset "@enzyme_reaction" begin
        spec = @enzyme_reaction begin
            substrates: S(C=1)
            products:   P(C=1)
        end
        @test length(spec.substrates) == 1
        @test spec.substrates[1].name == :S
        @test spec.substrates[1].atoms == Dict(:C => 1)
        @test length(spec.products) == 1
        @test spec.products[1].name == :P
        @test isempty(spec.regulators)

        spec2 = @enzyme_reaction begin
            substrates: S(C=6, H=12, O=6), ATP(C=10, H=16, N=5, O=13, P=3)
            products:   G6P(C=6, H=13, O=9, P=1), ADP(C=10, H=15, N=5, O=10, P=2)
            regulators: I(C=5, H=8, N=2)
        end
        @test length(spec2.substrates) == 2
        @test length(spec2.products) == 2
        @test length(spec2.regulators) == 1
        @test spec2.regulators[1].name == :I
    end

    @testset "@mechanism" begin
        m = @mechanism begin
            [E, S(C=1)] --> [ES]
            [ES] --> [E, P(C=1)]
        end
        @test m isa EnzymeMechanism
        @test n_steps(m) == 2
        @test n_states(m) == 2
        @test Set(s.name for s in enzyme_forms(m)) == Set([:E, :ES])
        @test Set(s.name for s in metabolites(m)) == Set([:S, :P])
        @test validate(m) == true

        # Verify species roles and atoms
        raw = steps(m)
        all_species = vcat(raw[1].first, raw[1].second, raw[2].first, raw[2].second)
        e_sp = filter(s -> s.name == :E, all_species)[1]
        s_sp = filter(s -> s.name == :S, all_species)[1]
        @test e_sp.role == enzyme
        @test s_sp.role == metabolite
        @test s_sp.atoms == Dict(:C => 1)

        # Numeric check: same as Uni-Uni spot check
        params = (k1f=3.2, k1r=0.8, k2f=2.5, k2r=1.1)
        concs = (S=0.7, P=0.3)
        @test rate_equation(m, params, concs) ≈ 0.9091 atol=0.001

        # Multi-step mechanism
        m2 = @mechanism begin
            [E, A(C=2, N=1)] --> [EA]
            [EA] --> [FP]
            [FP] --> [F, P(C=2)]
            [F, B(C=3)] --> [FB]
            [FB] --> [EQ]
            [EQ] --> [E, Q(C=3, N=1)]
        end
        @test n_states(m2) == 6
        @test validate(m2) == true
    end
end
