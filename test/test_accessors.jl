@testset "Accessor performance" begin
    m = @enzyme_mechanism begin
        substrates: S
        products:   P
        steps: begin
            E + S <--> E(S)
            E(S) <--> E + P
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
    parameters(m)
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

    @testset "parameters: zero-alloc and <100ns" begin
        @test (@allocated parameters(m)) == 0
        @test best_ns_per_call(parameters, m) < 100e-9
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

@testset "parameters API symmetry" begin
    rxn = @enzyme_reaction begin
        substrates: A[C], B[N]
        products:   P[C], Q[N]
    end
    m_mono = EnzymeRates.EnzymeMechanism(
        first(EnzymeRates.init_mechanisms(rxn)))
    full_mono = parameters(m_mono, Full)
    reduced_mono = parameters(m_mono, Reduced)
    @test :E_total in full_mono
    @test :E_total in reduced_mono
    @test :Keq in reduced_mono
    @test :Keq ∉ full_mono
    @test length(full_mono) >= length(reduced_mono)

    rxn_allo = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    base = first(EnzymeRates.init_mechanisms(rxn_allo))
    cat_allo_states = fill(:NonequalAI, length(EnzymeRates.kinetic_groups(base)))
    site = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:R)], 2, [:NonequalAI])
    am = EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(base), copy(EnzymeRates.steps(base)),
        cat_allo_states, 2, [site])
    m_allo = EnzymeRates.AllostericEnzymeMechanism(am)
    full_allo = parameters(m_allo, Full)
    reduced_allo = parameters(m_allo, Reduced)
    @test :L in full_allo
    @test :E_total in full_allo
    @test :Keq ∉ full_allo
    @test any(occursin("_I_", string(p)) for p in full_allo)
    @test any(occursin("reg", string(p)) for p in full_allo)
    @test :L in reduced_allo
    @test :Keq in reduced_allo
    @test :E_total in reduced_allo
end

@testset "added field accessors: cat_allo_states" begin
    rxn = @enzyme_reaction begin
        substrates: S[C]
        products:   P[C]
        allosteric_regulators: R
        oligomeric_state: 2
    end
    base = first(EnzymeRates.init_mechanisms(rxn))
    cas = fill(:NonequalAI, length(EnzymeRates.kinetic_groups(base)))
    site = EnzymeRates.RegulatorySite(
        [EnzymeRates.AllostericRegulator(:R)], 2, [:NonequalAI])
    am = EnzymeRates.AllostericMechanism(
        EnzymeRates.reaction(base), copy(EnzymeRates.steps(base)), cas, 2, [site])
    @test EnzymeRates.cat_allo_states(am) == cas
end
