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
        m = EnzymeMechanism(species, (((:E, :S), (:ES,)), ((:ES,), (:E, :P))))
        @test sprint(show, m) == "EnzymeMechanism: E + S ⇌ ES ⇌ E + P"

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
        m_b = EnzymeMechanism(species_b, rxns_b)
        str = sprint(show, m_b)
        @test startswith(str, "EnzymeMechanism (7 steps, 6 enzyme forms):")
        @test contains(str, "E + S1 ⇌ ES1")
        @test contains(str, "EP2 ⇌ E + P2")
    end

    @testset "EnzymeMechanism canonical ordering" begin
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
        m1 = EnzymeMechanism(species, rxns1)
        m2 = EnzymeMechanism(species, rxns2)
        @test typeof(m1) === typeof(m2)

        # Reversed species within reaction sides (metabolite before enzyme)
        rxns3 = (
            ((:S, :E), (:ES,)),
            ((:ES,), (:P, :E)),
        )
        m3 = EnzymeMechanism(species, rxns3)
        @test typeof(m1) === typeof(m3)
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
        @test_throws ErrorException EnzymeMechanism(base_species, ())

        # Duplicate substrate names in species
        dup_subs_species = (
            ((:S, ((:C, 1),)), (:S, ((:C, 1),))),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        @test_throws ErrorException EnzymeMechanism(dup_subs_species, base_rxns)

        # No free enzyme form (all enzymes have atoms)
        no_free_species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ((:X, 1),)), (:ES, ((:C, 1), (:X, 1)))),
        )
        @test_throws ErrorException EnzymeMechanism(no_free_species, base_rxns)

        # Reaction with zero enzymes on LHS
        no_enz_rxns = (((:S,), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(base_species, no_enz_rxns)

        # Reaction with two metabolites on one side
        two_met_species = (
            ((:S1, ((:C, 1),)), (:S2, ((:H, 1),))),
            ((:P, ((:C, 1), (:H, 1))),),
            (),
            ((:E, ()), (:ES, ((:C, 1), (:H, 1)))),
        )
        two_met_rxns = (((:E, :S1, :S2), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(two_met_species, two_met_rxns)

        # Unknown species in reaction
        unknown_rxns = (((:E, :X), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(base_species, unknown_rxns)

        # Atomic conservation failure
        bad_atom_species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 2),)),),  # different atoms than S
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        bad_atom_rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(bad_atom_species, bad_atom_rxns)

        # Net stoichiometry mismatch (substrate consumed but not produced)
        net_mismatch_species = (
            ((:S, ((:C, 1),)), (:S2, ((:H, 1),))),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        net_mismatch_rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(net_mismatch_species, net_mismatch_rxns)

        # Species defined as both enzyme and metabolite
        overlap_species = (
            ((:E, ((:C, 1),)),),  # E is also an enzyme
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        @test_throws ErrorException EnzymeMechanism(overlap_species, base_rxns)

        # Duplicate reactions
        dup_rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)), ((:E, :S), (:ES,)))
        @test_throws ErrorException EnzymeMechanism(base_species, dup_rxns)

        # Unreachable enzyme form
        unreachable_species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),)), (:EX, ((:H, 1),))),
        )
        unreachable_rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
        @test_throws ErrorException EnzymeMechanism(unreachable_species, unreachable_rxns)
    end

    @testset "EnzymeMechanism valid with reachable enzyme forms" begin
        species = (
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            (),
            ((:E, ()), (:ES, ((:C, 1),))),
        )
        rxns = (((:E, :S), (:ES,)), ((:ES,), (:E, :P)))
        m = EnzymeMechanism(species, rxns)
        @test m isa EnzymeMechanism
    end
end
