@testset "Enumeration (simple Uni-Uni)" begin
    spec = @enzyme_reaction begin
        substrates: S(C=1)
        products:   P(C=1)
    end

    # A simple Uni-Uni should have 3 independent params (2 steps, 1 Haldane)
    mechanisms = enumerate_mechanisms(spec, 3)
    @test length(mechanisms) >= 1

    # Each mechanism should be valid
    for m in mechanisms
        @test validate(m) == true
        @test n_independent_params(m) == 3
    end
end
