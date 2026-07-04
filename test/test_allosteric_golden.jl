# ABOUTME: Byte-identical golden reference for allosteric rate-equation strings
# ABOUTME: and parameter lists; guards the D1 state-parameterized re-derivation.
using Test

const _ALLO_GOLDEN_PATH =
    joinpath(@__DIR__, "reference", "allosteric_golden_reference.txt")

"""Canonical serialization of every allosteric spec's derivation output."""
function _allosteric_golden_lines()
    lines = String[]
    for spec in MECHANISM_TEST_SPECS
        spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism || continue
        m = spec.mechanism
        push!(lines, "### " * spec.name)
        reduced_string = EnzymeRates.rate_equation_string(m, EnzymeRates.Reduced)
        push!(lines, "REDUCED_STRING " * replace(reduced_string, "\n" => "\\n"))
        push!(lines, "PARAMS_FULL " * string(parameters(m, EnzymeRates.Full)))
        push!(lines, "PARAMS_REDUCED " *
              string(parameters(m, EnzymeRates.Reduced)))
    end
    lines
end

@testset "allosteric golden reference (D1)" begin
    @test isfile(_ALLO_GOLDEN_PATH)
    current = _allosteric_golden_lines()
    reference = readlines(_ALLO_GOLDEN_PATH)
    @test length(current) == length(reference)
    for (c, r) in zip(current, reference)
        @test c == r
    end
end
