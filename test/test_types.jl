@testset "Types" begin
    @testset "EnzymeReaction" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
        )
        @test r isa EnzymeReaction
        @test substrates(r) == ((:S, ((:C, 1),)),)
        @test products(r) == ((:P, ((:C, 1),)),)
        @test regulators(r) == ()
    end

    @testset "EnzymeReaction with regulators" begin
        r = EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            ((:I, ((:C, 2),)),),
        )
        @test regulators(r) == ((:I, ((:C, 2),)),)
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
        @test substrates(r1) == ((:S1, ((:C, 1),)), (:S2, ((:C, 2),)))
    end

    @testset "EnzymeReaction validation" begin
        @test_throws ErrorException EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            ((:S, ((:C, 1),)),),  # regulator same as substrate
        )
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
end
