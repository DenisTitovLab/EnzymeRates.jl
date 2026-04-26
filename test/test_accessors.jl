@testset "Accessor performance" begin
    m = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            [E, S] <--> [ES]
            [ES] <--> [E, P]
        end
    end

    # Warmup - use qualified names for internal functions
    EnzymeRates.substrates(m)
    EnzymeRates.products(m)
    EnzymeRates.regulators(m)
    EnzymeRates.enzyme_forms(m)
    EnzymeRates.reactions(m)
    EnzymeRates.n_states(m)
    EnzymeRates.n_steps(m)
    metabolites(m)
    EnzymeRates.stoich_matrix(m); EnzymeRates.equilibrium_steps(m)

    # Use minimum of multiple timing runs to avoid GC noise
    function best_ns_per_call(f, arg; n=10_000, trials=5)
        minimum(begin
            t = @elapsed for _ in 1:n
                f(arg)
            end
            t / n
        end for _ in 1:trials)
    end

    @testset "substrates: zero-alloc and <100ns" begin
        @test (@allocated EnzymeRates.substrates(m)) == 0
        @test best_ns_per_call(EnzymeRates.substrates, m) < 100e-9
    end

    @testset "products: zero-alloc and <100ns" begin
        @test (@allocated EnzymeRates.products(m)) == 0
        @test best_ns_per_call(EnzymeRates.products, m) < 100e-9
    end

    @testset "regulators: zero-alloc and <100ns" begin
        @test (@allocated EnzymeRates.regulators(m)) == 0
        @test best_ns_per_call(EnzymeRates.regulators, m) < 100e-9
    end

    @testset "enzyme_forms: zero-alloc and <100ns" begin
        @test (@allocated EnzymeRates.enzyme_forms(m)) == 0
        @test best_ns_per_call(EnzymeRates.enzyme_forms, m) < 100e-9
    end

    @testset "reactions: zero-alloc and <100ns" begin
        @test (@allocated EnzymeRates.reactions(m)) == 0
        @test best_ns_per_call(EnzymeRates.reactions, m) < 100e-9
    end

    @testset "n_states: zero-alloc and <100ns" begin
        @test (@allocated EnzymeRates.n_states(m)) == 0
        @test best_ns_per_call(EnzymeRates.n_states, m) < 100e-9
    end

    @testset "n_steps: zero-alloc and <100ns" begin
        @test (@allocated EnzymeRates.n_steps(m)) == 0
        @test best_ns_per_call(EnzymeRates.n_steps, m) < 100e-9
    end

    # TODO: Phase 3 — restore once Task 2.7 migrates the parameters() accessor.
    @static if false
        @testset "parameters: zero-alloc and <100ns" begin
            @test (@allocated parameters(m)) == 0
            @test best_ns_per_call(parameters, m) < 100e-9
        end
    end

    @testset "metabolites: zero-alloc and <100ns" begin
        @test (@allocated metabolites(m)) == 0
        @test best_ns_per_call(metabolites, m) < 100e-9
    end

    @testset "stoich_matrix: zero-alloc and <100ns" begin
        @test (@allocated EnzymeRates.stoich_matrix(m)) == 0
        @test best_ns_per_call(EnzymeRates.stoich_matrix, m) < 100e-9
    end

    @testset "equilibrium_steps: zero-alloc and <100ns" begin
        @test (@allocated EnzymeRates.equilibrium_steps(m)) == 0
        @test best_ns_per_call(EnzymeRates.equilibrium_steps, m) < 100e-9
    end
end
