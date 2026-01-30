@testset "Types" begin
    E = Species(:E, enzyme)
    S = Species(:S, metabolite, Dict(:C => 1))
    @test E.role == enzyme
    @test S.role == metabolite
    @test S.atoms == Dict(:C => 1)
    @test isempty(E.atoms)
end
