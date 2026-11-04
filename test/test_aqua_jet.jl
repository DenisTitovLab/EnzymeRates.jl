# ABOUTME: Code-quality and static-analysis checks for EnzymeRates: Aqua (package
# ABOUTME: hygiene), file names that differ only in case, and JET (type inference).
using Aqua
using JET

@testset "Aqua" begin
    # The Downgrade workflow's action moves the test dependencies into [deps], where they
    # would read as stale.
    Aqua.test_all(EnzymeRates; ambiguities=false, stale_deps=(ignore=[:Aqua, :JET,
        :OptimizationBBO, :OptimizationCMAEvolutionStrategy, :OrdinaryDiffEqFIRK],))
end

@testset "no two files differ only in case" begin
    # A case-insensitive file system (the Windows and macOS default) keeps only one of them.
    root = pkgdir(EnzymeRates)
    paths = [relpath(joinpath(dir, f), root) for top in ("src", "test", "docs")
             for (dir, _, files) in walkdir(joinpath(root, top)) for f in files]
    lowered = lowercase.(paths)
    @test [p for p in paths if count(==(lowercase(p)), lowered) > 1] == String[]
end

@testset "JET" begin
    # JET 0.9 matches a module in `target_modules` exactly, not its submodules, so each
    # submodule is listed beside the package.
    JET.test_package(EnzymeRates; target_modules=(EnzymeRates, EnzymeRates.Mechanisms,
        EnzymeRates.Derivation, EnzymeRates.Enumeration, EnzymeRates.Fitting))
end
