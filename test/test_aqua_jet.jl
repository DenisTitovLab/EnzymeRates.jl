# ABOUTME: Code-quality and static-analysis checks for EnzymeRates via
# ABOUTME: Aqua (package hygiene) and JET (type/inference analysis).
using Aqua
using JET

@testset "Aqua" begin
    Aqua.test_all(EnzymeRates; ambiguities=false)
end

@testset "JET" begin
    # JET 0.9 matches a module in `target_modules` exactly, not its submodules, so each
    # submodule is listed beside the package.
    JET.test_package(EnzymeRates; target_modules=(EnzymeRates, EnzymeRates.Mechanisms,
        EnzymeRates.Derivation, EnzymeRates.Enumeration, EnzymeRates.Fitting))
end
