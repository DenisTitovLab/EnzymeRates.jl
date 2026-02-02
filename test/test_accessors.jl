@testset "Accessor performance" begin
    species = (
        ((:S, ((:C, 1),)),),
        ((:P, ((:C, 1),)),),
        (),
        ((:E, ()), (:ES, ((:C, 1),))),
    )
    rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
    m = EnzymeMechanism(species, rxns)

    # Warmup
    substrates(m); products(m); regulators(m); enzyme_forms(m)
    reactions(m); n_states(m); n_steps(m); parameters(m); metabolites(m)
    graph(m); stoich_matrix(m)

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
        @test (@allocated substrates(m)) == 0
        @test best_ns_per_call(substrates, m) < 100e-9
    end

    @testset "products: zero-alloc and <100ns" begin
        @test (@allocated products(m)) == 0
        @test best_ns_per_call(products, m) < 100e-9
    end

    @testset "regulators: zero-alloc and <100ns" begin
        @test (@allocated regulators(m)) == 0
        @test best_ns_per_call(regulators, m) < 100e-9
    end

    @testset "enzyme_forms: zero-alloc and <100ns" begin
        @test (@allocated enzyme_forms(m)) == 0
        @test best_ns_per_call(enzyme_forms, m) < 100e-9
    end

    @testset "reactions: zero-alloc and <100ns" begin
        @test (@allocated reactions(m)) == 0
        @test best_ns_per_call(reactions, m) < 100e-9
    end

    @testset "n_states: zero-alloc and <100ns" begin
        @test (@allocated n_states(m)) == 0
        @test best_ns_per_call(n_states, m) < 100e-9
    end

    @testset "n_steps: zero-alloc and <100ns" begin
        @test (@allocated n_steps(m)) == 0
        @test best_ns_per_call(n_steps, m) < 100e-9
    end

    @testset "parameters: zero-alloc and <100ns" begin
        @test (@allocated parameters(m)) == 0
        @test best_ns_per_call(parameters, m) < 100e-9
    end

    @testset "metabolites: zero-alloc and <100ns" begin
        @test (@allocated metabolites(m)) == 0
        @test best_ns_per_call(metabolites, m) < 100e-9
    end

    @testset "graph: zero-alloc and <100ns" begin
        @test (@allocated graph(m)) == 0
        @test best_ns_per_call(graph, m) < 100e-9
    end

    @testset "stoich_matrix: zero-alloc and <100ns" begin
        @test (@allocated stoich_matrix(m)) == 0
        @test best_ns_per_call(stoich_matrix, m) < 100e-9
    end
end
