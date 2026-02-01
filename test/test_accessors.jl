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

    @testset "substrates: zero-alloc and <100ns" begin
        allocs = @allocated substrates(m)
        t = @elapsed for _ in 1:10_000; substrates(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "products: zero-alloc and <100ns" begin
        allocs = @allocated products(m)
        t = @elapsed for _ in 1:10_000; products(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "regulators: zero-alloc and <100ns" begin
        allocs = @allocated regulators(m)
        t = @elapsed for _ in 1:10_000; regulators(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "enzyme_forms: zero-alloc and <100ns" begin
        allocs = @allocated enzyme_forms(m)
        t = @elapsed for _ in 1:10_000; enzyme_forms(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "reactions: zero-alloc and <100ns" begin
        allocs = @allocated reactions(m)
        t = @elapsed for _ in 1:10_000; reactions(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "n_states: zero-alloc and <100ns" begin
        allocs = @allocated n_states(m)
        t = @elapsed for _ in 1:10_000; n_states(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "n_steps: zero-alloc and <100ns" begin
        allocs = @allocated n_steps(m)
        t = @elapsed for _ in 1:10_000; n_steps(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "parameters: zero-alloc and <100ns" begin
        allocs = @allocated parameters(m)
        t = @elapsed for _ in 1:10_000; parameters(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "metabolites: zero-alloc and <100ns" begin
        allocs = @allocated metabolites(m)
        t = @elapsed for _ in 1:10_000; metabolites(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "graph: zero-alloc and <100ns" begin
        allocs = @allocated graph(m)
        t = @elapsed for _ in 1:10_000; graph(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end

    @testset "stoich_matrix: zero-alloc and <100ns" begin
        allocs = @allocated stoich_matrix(m)
        t = @elapsed for _ in 1:10_000; stoich_matrix(m); end
        @test allocs == 0
        @test t / 10_000 < 100e-9
    end
end
