@testset "Types" begin
    @testset "EnzymeReaction" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
        )
        @test r isa EnzymeReaction
        @test EnzymeRates.substrates(r) == ((:S, ((:C, 1),)),)
        @test EnzymeRates.products(r) == ((:P, ((:C, 1),)),)
        @test EnzymeRates.regulators(r) == ()
    end

    @testset "EnzymeReaction with regulators" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            ((:I, ((:C, 2),)),),
        )
        @test EnzymeRates.regulators(r) == ((:I, ((:C, 2),)),)
    end

    @testset "EnzymeReaction canonical ordering" begin
        r1 = EnzymeReaction(
            ((:S2, ((:C, 2),)), (:S1, ((:C, 1),))),
            ((:P1, ((:C, 1),)), (:P2, ((:C, 2),))),
        )
        r2 = EnzymeReaction(
            ((:S1, ((:C, 1),)), (:S2, ((:C, 2),))),
            ((:P1, ((:C, 1),)), (:P2, ((:C, 2),))),
        )
        @test typeof(r1) === typeof(r2)
        @test EnzymeRates.substrates(r1) == ((:S1, ((:C, 1),)), (:S2, ((:C, 2),)))
    end

    @testset "EnzymeReaction regulator same as substrate allowed" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            ((:S, ((:C, 1),)),),  # regulator same as substrate - allowed
        )
        @test r isa EnzymeReaction
    end

    @testset "EnzymeReaction regulator same as product allowed" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),  # regulator same as product - allowed
        )
        @test r isa EnzymeReaction
    end

    @testset "EnzymeReaction duplicate substrate names" begin
        @test_throws ErrorException EnzymeReaction(
            ((:S, ((:C, 1),)), (:S, ((:C, 1),))),
            ((:P, ((:C, 1),)),),
        )
    end

    @testset "EnzymeReaction duplicate product names" begin
        @test_throws ErrorException EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)), (:P, ((:C, 1),))),
        )
    end

    @testset "Pretty printing" begin
        r = EnzymeReaction(((:S, ((:C, 1),)),), ((:P, ((:C, 1),)),))
        @test sprint(show, r) == "EnzymeReaction: S ⇌ P"

        r2 = EnzymeReaction(
            ((:S, ((:C, 1),)), (:ATP, ((:C, 10),))),
            ((:P, ((:C, 1),)), (:ADP, ((:C, 10),))),
            ((:I, ((:C, 5),)),),
        )
        @test sprint(show, r2) == "EnzymeReaction: ATP + S ⇌ ADP + P | regulators: I"

        # Linear mechanism: compact chain
        species = (
            ((:S, ((:C, 1),)),), ((:P, ((:C, 1),)),), (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        m = EnzymeMechanism(species, (((:E, :S), (:ES,)), ((:ES,), (:E, :P))), (false, false))
        @test sprint(show, m) == "EnzymeMechanism: E + S <--> ES <--> E + P"

        # Branched mechanism: multi-line
        species_b = (
            ((:S1, ((:C, 1),)), (:S2, ((:H, 1),))),
            ((:P1, ((:C, 1),)), (:P2, ((:H, 1),))),
            (),
            ((:E, ()), (:ES1, ((:C, 1),)), (:ES2, ((:H, 1),)),
             (:ES1S2, ((:C, 1), (:H, 1))), (:EP1P2, ((:C, 1), (:H, 1))),
             (:EP2, ((:H, 1),))),
        )
        rxns_b = (
            ((:E, :S1), (:ES1,)), ((:E, :S2), (:ES2,)),
            ((:ES1, :S2), (:ES1S2,)), ((:ES2, :S1), (:ES1S2,)),
            ((:ES1S2,), (:EP1P2,)), ((:EP1P2,), (:EP2, :P1)), ((:EP2,), (:E, :P2)),
        )
        m_b = EnzymeMechanism(species_b, rxns_b, ntuple(Returns(false), length(rxns_b)))
        str = sprint(show, m_b)
        @test startswith(str, "EnzymeMechanism (7 steps, 6 enzyme forms):")
        @test contains(str, "E + S1 <--> ES1")
        @test contains(str, "EP2 <--> E + P2")
    end

    @testset "EnzymeMechanism different orderings produce valid mechanisms" begin
        species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        # Normal order
        rxns1 = (
            ((:E, :S), (:ES,)),
            ((:ES,), (:E, :P)),
        )
        # Reversed reaction order
        rxns2 = (
            ((:ES,), (:E, :P)),
            ((:E, :S), (:ES,)),
        )
        m1 = EnzymeMechanism(species, rxns1, (false, false))
        m2 = EnzymeMechanism(species, rxns2, (false, false))
        @test m1 isa EnzymeMechanism
        @test m2 isa EnzymeMechanism
    end

    @testset "EnzymeMechanism error cases" begin
        base_species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        base_rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))

        # Empty reactions tuple
        @test_throws ErrorException EnzymeMechanism(base_species, (), ())

        # Duplicate substrate names in species
        dup_subs_species = (
            ((:S, ((:C, 1),)), (:S, ((:C, 1),))),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        @test_throws ErrorException EnzymeMechanism(dup_subs_species, base_rxns, (false, false))

        # No free enzyme form (all enzymes have atoms)
        no_free_species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ((:X, 1),)), (:ES, ((:C, 1), (:X, 1)))),
        )
        @test_throws ErrorException EnzymeMechanism(no_free_species, base_rxns, (false, false))

        # Reaction with zero enzymes on LHS
        no_enz_rxns = (((:S,), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(base_species, no_enz_rxns, (false, false))

        # Reaction with two metabolites on one side
        two_met_species = (
            ((:S1, ((:C, 1),)), (:S2, ((:H, 1),))),
            ((:P, ((:C, 1), (:H, 1))),),
            (),
            ((:E, ()), (:ES, ((:C, 1), (:H, 1)))),
        )
        two_met_rxns = (((:E, :S1, :S2), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(two_met_species, two_met_rxns, (false, false))

        # Unknown species in reaction
        unknown_rxns = (((:E, :X), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(base_species, unknown_rxns, (false, false))

        # Atomic conservation failure
        bad_atom_species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 2),)),),  # different atoms than S
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        bad_atom_rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(bad_atom_species, bad_atom_rxns, (false, false))

        # Net stoichiometry mismatch (substrate consumed but not produced)
        net_mismatch_species = (
            ((:S, ((:C, 1),)), (:S2, ((:H, 1),))),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        net_mismatch_rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(net_mismatch_species, net_mismatch_rxns, (false, false))

        # Species defined as both enzyme and metabolite
        overlap_species = (
            ((:E, ((:C, 1),)),),  # E is also an enzyme
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        @test_throws ErrorException EnzymeMechanism(overlap_species, base_rxns, (false, false))

        # Duplicate reactions
        dup_rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)), ((:E, :S), (:ES,)))
        @test_throws ErrorException EnzymeMechanism(base_species, dup_rxns, (false, false, false))

        # Unreachable enzyme form
        unreachable_species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),)), (:EX, ((:H, 1),))),
        )
        unreachable_rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(unreachable_species, unreachable_rxns, (false, false))
    end

    @testset "EnzymeMechanism valid with reachable enzyme forms" begin
        species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
        m = EnzymeMechanism(species, rxns, (false, false))
        @test m isa EnzymeMechanism
    end

    @testset "Constraint constructor validation errors" begin
        species = (
            ((:A, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:EA, ((:C, 1),)), (:EP, ((:C, 1),))),
        )
        rxns = (((:E, :A), (:EA,)), ((:EA,), (:EP,)), ((:EP,), (:E, :P)))
        eq = (false, false, false)

        # Invalid target
        @test_throws ErrorException EnzymeMechanism(species, rxns, eq,
            ((:K99, 1, ((:k1f, 1),)),))

        # Self-reference
        @test_throws ErrorException EnzymeMechanism(species, rxns, eq,
            ((:k3r, 1, ((:k3r, 1),)),))

        # Duplicate target
        @test_throws ErrorException EnzymeMechanism(species, rxns, eq,
            ((:k3r, 1, ((:k1r, 1),)), (:k3r, 1, ((:k2r, 1),))))

        # Invalid replacement symbol
        @test_throws ErrorException EnzymeMechanism(species, rxns, eq,
            ((:k3r, 1, ((:K99, 1),)),))

        # Zero coefficient
        @test_throws ErrorException EnzymeMechanism(species, rxns, eq,
            ((:k3r, 0, ((:k1r, 1),)),))
    end

    @testset "parameters() excludes constrained params" begin
        m = @enzyme_mechanism begin
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
        raw_params = parameters(m, EnzymeRates.Raw)
        @test :k3r ∉ raw_params
        @test :k1r ∈ raw_params
        @test :E_total ∈ raw_params
    end

    @testset "Constrained Uni-Uni k3r = k1r: rate equation correctness" begin
        m = @enzyme_mechanism begin
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

        function rate_constrained(p, c)
            (; k1f, k1r, k2f, k2r, k3f, E_total) = p
            (; A, P) = c
            k3r = k1r
            num = k1f * k2f * k3f * A - k1r * k2r * k3r * P
            denom = (k1r * k2r + k1r * k3f + k2f * k3f) +
                    k1f * (k2r + k2f + k3f) * A +
                    k3r * (k1r + k2r + k2f) * P
            E_total * num / denom
        end

        rng = Random.MersenneTwister(42)
        for _ in 1:20
            params = (k1f = 0.1 + 9.9 * rand(rng),
                      k1r = 0.1 + 9.9 * rand(rng),
                      k2f = 0.1 + 9.9 * rand(rng),
                      k2r = 0.1 + 9.9 * rand(rng),
                      k3f = 0.1 + 9.9 * rand(rng),
                      E_total = 0.1 + 9.9 * rand(rng))
            concs = (A = 0.1 + 9.9 * rand(rng), P = 0.1 + 9.9 * rand(rng))
            @test rate_equation(m, params, concs, EnzymeRates.Raw) ≈ rate_constrained(params, concs) rtol=1e-10
        end
    end

    @testset "Constrained RE K2 = K1: rate equation correctness" begin
        m = @enzyme_mechanism begin
            species: begin
                substrates: A[C], B[N]
                products:   P[C], Q[N]
                enzymes:    E, EA[C], EABEPQ[CN], EQ[N]
            end
            steps: begin
                [E, A] ⇌ [EA]
                [EA, B] ⇌ [EABEPQ]
                [EABEPQ] <--> [EQ, P]
                [EQ] <--> [E, Q]
            end
            constraints: begin
                K2 = K1
            end
        end

        function rate_constrained_re(p, c)
            (; K1, k3f, k3r, k4f, k4r, E_total) = p
            (; A, B, P, Q) = c
            num = k3f * k4f * A * B / (K1 * K1) - k3r * k4r * P * Q
            sigma1 = 1.0 + A / K1 + A * B / (K1 * K1)
            R12 = k3f * A * B / (K1 * K1) + k4r * Q
            R21 = k3r * P + k4f
            denom = sigma1 * R21 + R12
            E_total * num / denom
        end

        raw_params = parameters(m, EnzymeRates.Raw)
        @test :K2 ∉ raw_params
        @test :K1 ∈ raw_params

        rng = Random.MersenneTwister(99)
        for _ in 1:20
            params = (K1 = 0.1 + 9.9 * rand(rng),
                      k3f = 0.1 + 9.9 * rand(rng),
                      k3r = 0.1 + 9.9 * rand(rng),
                      k4f = 0.1 + 9.9 * rand(rng),
                      k4r = 0.1 + 9.9 * rand(rng),
                      E_total = 0.1 + 9.9 * rand(rng))
            concs = (A = 0.1 + 9.9 * rand(rng), B = 0.1 + 9.9 * rand(rng),
                     P = 0.1 + 9.9 * rand(rng), Q = 0.1 + 9.9 * rand(rng))
            @test rate_equation(m, params, concs, EnzymeRates.Raw) ≈ rate_constrained_re(params, concs) rtol=1e-10
        end
    end

    @testset "Constrained equilibrium (v=0 at Keq)" begin
        m = @enzyme_mechanism begin
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

        hw_params = parameters(m)
        @test :k3r ∉ hw_params
        @test :Keq ∈ hw_params

        rng = Random.MersenneTwister(77)
        hw = parameters(m)
        vals = Tuple(0.1 + 9.9 * rand(rng) for _ in hw)
        p = NamedTuple{hw}(vals)
        Keq = p.Keq

        eq_concs = (A = 1.0, P = Keq)
        v = rate_equation(m, p, eq_concs)
        @test abs(v) < 1e-10
    end

    @testset "Constrained zero allocation" begin
        m = @enzyme_mechanism begin
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
        rng = Random.MersenneTwister(55)
        hw = parameters(m)
        vals = Tuple(0.1 + 9.9 * rand(rng) for _ in hw)
        p = NamedTuple{hw}(vals)
        concs = (A = 1.0, P = 2.0)

        rate_equation(m, p, concs)  # warmup
        allocs = @allocated rate_equation(m, p, concs)
        @test allocs == 0
    end

    @testset "rate_equation_string shows constraints" begin
        m = @enzyme_mechanism begin
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

        s_raw = rate_equation_string(m, EnzymeRates.Raw)
        @test occursin("k3r = k1r", s_raw)
        @test occursin("v = E_total * (", s_raw)
        @test !occursin("k3r", replace(s_raw, "k3r = k1r" => ""))

        s_hw = rate_equation_string(m)
        @test occursin("k3r = k1r", s_hw)
    end

    @testset "show method displays constraints" begin
        m = @enzyme_mechanism begin
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
        s = sprint(show, m)
        @test occursin("constraints:", s)
        @test occursin("k3r = k1r", s)
    end
end
