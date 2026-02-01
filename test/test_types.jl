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

    @testset "EnzymeReaction validation" begin
        @test_throws ErrorException EnzymeReaction(
            ((:S, ((:C, 1),)),),
            ((:P, ((:C, 1),)),),
            ((:S, ((:C, 1),)),),  # regulator same as substrate
        )
    end
end
