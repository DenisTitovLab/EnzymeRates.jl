# ABOUTME: Code-quality and static-analysis checks for EnzymeRates via
# ABOUTME: Aqua (package hygiene) and JET (type/inference analysis).
using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(EnzymeRates; ambiguities=false)
end

@testset "JET" begin
    JET.test_package(EnzymeRates; target_modules=(EnzymeRates,))
end
