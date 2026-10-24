# ABOUTME: Byte-identical golden reference for the allosteric specs' printed laws, in full
# ABOUTME: and in coefficients, and for their parameter lists.
using Test

const _ALLO_GOLDEN_PATH =
    joinpath(@__DIR__, "reference", "allosteric_golden_reference.txt")

"""Canonical serialization of every allosteric spec's derivation output."""
function _testhelper_allosteric_golden_lines()
    lines = String[]
    for spec in MECHANISM_TEST_SPECS
        spec.mechanism isa EnzymeRates.AllostericEnzymeMechanism || continue
        m = spec.mechanism
        push!(lines, "### " * spec.name)
        eq_string = EnzymeRates.rate_equation_string(m)
        push!(lines, "RATE_EQUATION_STRING " * replace(eq_string, "\n" => "\\n"))
        push!(lines, "RATE_EQUATION_COEFFICIENTS " *
                     replace(EnzymeRates.rate_equation_string(m; coefficients = true),
                             "\n" => "\\n"))
    end
    lines
end

@testset "allosteric golden reference (D1)" begin
    # The per-spec mass-action oracle of test_rate_eq_derivation.jl
    # (`_testhelper_test_mass_action_oracle`) checks each spec's law, with one catalytic
    # subunit, against the exact steady state of the spec's mass-action network; the
    # per-spec string check (`_testhelper_test_rate_equation_string`) ties the printed
    # string to the compiled law, and test_print.jl ties both printed forms to the
    # derived law.
    @test isfile(_ALLO_GOLDEN_PATH)
    current = _testhelper_allosteric_golden_lines()
    reference = readlines(_ALLO_GOLDEN_PATH)
    @test length(current) == length(reference)
    for (c, r) in zip(current, reference)
        @test c == r
    end
    # The RATE_EQUATION_STRING header pins the fitted parameters; this pins the
    # exported `parameters(m)` wrapper to the fitted parameters plus Keq and E_total.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        m isa EnzymeRates.AllostericEnzymeMechanism || continue
        @test parameters(m) ==
              (EnzymeRates.fitted_params(m)..., :Keq, :E_total)
    end
    # The fitted parameters list the catalytic constants by name, so rate-equivalent
    # mechanisms render one params line, then the regulator constants, then L.
    for spec in MECHANISM_TEST_SPECS
        m = spec.mechanism
        m isa EnzymeRates.AllostericEnzymeMechanism || continue
        fp = EnzymeRates.fitted_params(m)
        am = EnzymeRates.AllostericMechanism(m)
        reg = Tuple(EnzymeRates.name(p, am) for st in (:A, :I)
                    for p in EnzymeRates._kreg_params(am, st))
        cat = fp[1:end - 1 - count(in(fp), reg)]
        @test issorted(collect(cat); by = string)
        @test fp == (cat..., filter(in(fp), reg)..., :L)
    end
end
