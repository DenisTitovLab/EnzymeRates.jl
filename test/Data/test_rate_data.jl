# ABOUTME: Tests for RateData: the validated rate table, Keq and reference anchors that
# ABOUTME: FittingProblem and IdentifyRateEquationProblem hold.

@testset "RateData" begin
    RateData = EnzymeRates.Data.RateData
    table = (group = ["G1", "G1"], Rate = [1.0, 2.0], S = [1.0, 2.0], P = [0.1, 0.1])
    ref = [(S = 1.0, P = 0.1) => 1.0]

    @testset "RateData holds the validated table" begin
        rd = RateData(table, (:S, :P); Keq = 2.0, reference = [(S = 1.0,) => 3])
        @test rd.Keq === 2.0
        @test rd.absolute === false
        @test rd.data.Rate == table.Rate
        @test rd.reference == [(S = 1.0,) => 3.0]
        @test rd.reference isa Vector{Pair{NamedTuple, Float64}}
        # An integer Keq converts to the Float64 field.
        @test RateData(table, (:S, :P); Keq = 2, reference = ref).Keq === 2.0
        abs_rd = RateData(table, (:S, :P); Keq = 2.0, absolute = true)
        @test abs_rd.absolute === true && isempty(abs_rd.reference)
    end

    @testset "reference is required in relative mode and refused in absolute mode" begin
        @test_throws(ErrorException("relative mode needs reference = " *
                                    "[(concentrations...) => turnover, ...] to set the " *
                                    "rate scale; pass absolute = true for absolute " *
                                    "per-enzyme turnover data"),
                     RateData(table, (:S, :P); Keq = 2.0))
        @test_throws(ErrorException("reference applies only to relative mode; absolute " *
                                    "turnover data set the rate scale themselves"),
                     RateData(table, (:S, :P); Keq = 2.0, absolute = true,
                              reference = [(S = 10.0,) => 3.0]))
    end

    @testset "reference anchors are validated at construction" begin
        _testhelper_rate_data(reference; owner = "mechanism") =
            RateData(table, (:S, :P); Keq = 2.0, reference, owner)
        @test_throws(ErrorException("reference must be a vector of (concentrations => " *
                                    "turnover) pairs; got (S = 10.0,) => 3.0"),
                     _testhelper_rate_data((S = 10.0,) => 3.0))
        @test_throws(ErrorException("reference is empty; pass at least one " *
                                    "(concentrations => turnover) anchor"),
                     _testhelper_rate_data(Pair{NamedTuple, Float64}[]))
        @test_throws(ErrorException("reference: concentrations must be a NamedTuple; " *
                                    "got (10.0,)"),
                     _testhelper_rate_data([(10.0,) => 3.0]))
        @test_throws(ErrorException("reference: X is not a metabolite of the mechanism " *
                                    "(metabolites: S, P)"),
                     _testhelper_rate_data([(X = 1.0,) => 3.0]))
        @test_throws(ErrorException("reference: Q is not a metabolite of the reaction " *
                                    "(metabolites: S, P)"),
                     _testhelper_rate_data([(Q = 1.0,) => 3.0]; owner = "reaction"))
        @test_throws(ErrorException("reference: concentration of S must be a finite " *
                                    "number ≥ 0; got -1.0"),
                     _testhelper_rate_data([(S = -1.0,) => 3.0]))
        for (bad, shown) in ((0.0, "0.0"), (NaN, "NaN"), (Inf, "Inf"))
            @test_throws(ErrorException("reference: turnover must be a finite, nonzero " *
                                        "number; got $shown"),
                         _testhelper_rate_data([(S = 10.0,) => bad]))
        end
    end

    @testset "table validation" begin
        _testhelper_rate_data(table) = RateData(table, (:S, :P); Keq = 1.0, reference = ref)
        @test_throws ErrorException _testhelper_rate_data(
            (group = ["G1"], S = [1.0], P = [0.1]))
        @test_throws ErrorException _testhelper_rate_data(
            (group = ["G1"], Rate = [1.0], S = [1.0]))
        @test_throws ErrorException _testhelper_rate_data(
            (group = ["G1"], Rate = [0.0], S = [1.0], P = [0.1]))
        @test_throws ErrorException _testhelper_rate_data(
            (Rate = [1.0], S = [1.0], P = [0.1]))

        # A non-finite, missing or non-numeric rate: the error names the Rate column and
        # the row, and shows the value as Julia prints it, so a String rate is quoted.
        for (bad, shown) in ((NaN, "NaN"), (Inf, "Inf"), (-Inf, "-Inf"),
                             (missing, "missing"), ("1.0", "\"1.0\""))
            data_bad = (group = ["G1", "G1"], Rate = [1.0, bad], S = [1.0, 2.0],
                        P = [0.1, 0.1])
            @test_throws(
                ErrorException("Rate at row 2 must be a finite number; got $shown"),
                _testhelper_rate_data(data_bad))
        end

        # A concentration that is not a finite number ≥ 0: the error names the column
        # and the row. Zero is a valid concentration.
        for (bad, shown) in ((NaN, "NaN"), (Inf, "Inf"), (-1.0, "-1.0"),
                             (missing, "missing"), ("1.0", "\"1.0\""))
            data_bad = (group = ["G1", "G1"], Rate = [1.0, 2.0], S = [1.0, bad],
                        P = [0.1, 0.1])
            @test_throws(
                ErrorException(
                    "Concentration S at row 2 must be a finite number ≥ 0; got $shown"),
                _testhelper_rate_data(data_bad))
        end
        data_zero_S = (group = ["G1", "G1"], Rate = [1.0, 2.0], S = [0.0, 1.0],
                       P = [0.1, 0.1])
        @test _testhelper_rate_data(data_zero_S) isa RateData

        # A metabolite named after a required column would read that column, here the
        # measured rates, as its concentrations.
        @test_throws(
            ErrorException("Metabolite Rate has the name of the required Rate column; " *
                           "rename the metabolite"),
            RateData((group = ["G1"], Rate = [1.0], P = [0.1]), (:Rate, :P); Keq = 1.0,
                     reference = [(P = 0.1,) => 1.0]))

        # A Keq that is not positive and finite
        for (bad, shown) in ((0, "0"), (-1, "-1"), (Inf, "Inf"), (NaN, "NaN"))
            @test_throws(ErrorException("Keq must be positive and finite; got $shown"),
                RateData(table, (:S, :P); Keq = bad, reference = ref))
        end
    end
end
