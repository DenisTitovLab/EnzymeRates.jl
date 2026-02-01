using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(EnzymeRates; ambiguities=false)
end

@testset "JET" begin
    JET.test_package(EnzymeRates; target_modules=(EnzymeRates,))
end
