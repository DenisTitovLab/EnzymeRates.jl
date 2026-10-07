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
    # The REDUCED_STRING header pins the reduced parameters; this pins the exported
    # `parameters(m, Reduced)` wrapper to the fitted parameters plus Keq and E_total.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        m isa EnzymeRates.AllostericEnzymeMechanism || continue
        @test parameters(m, EnzymeRates.Reduced) ==
              (EnzymeRates.fitted_params(m)..., :Keq, :E_total)
    end
end
