# ABOUTME: Tests that algebraically-equivalent rate equations from different
# ABOUTME: mechanism specs collapse to one eq_hash (Sources A, B, C dedup).

using CSV, DataFrames, Test
using EnzymeRates
using EnzymeRates: _canonical_rate_eq_hash

const _CV_PATH = joinpath(@__DIR__, "..", "dedup_investigation", "cv_results.csv")

if !isfile(_CV_PATH)
    @warn "Skipping eq_hash dedup tests — $(_CV_PATH) not present."
else
    const _CV = CSV.read(_CV_PATH, DataFrame)

    # Reconstruct an EnzymeMechanism instance from a CSV row's mechanism_type column.
    _mech(row_idx) = eval(Meta.parse(_CV[row_idx, :mechanism_type]))()

    @testset "eq_hash dedup" begin
        # Sanity: confirm the row indices match the investigation's eq_hashes.
        # Row indices and full hashes resolved against the cv_results.csv
        # currently in the repo (20 rows). The plan's reference (rows 22,
        # 27, 31, 36 and 8-char prefixes) was computed against a larger
        # CSV; the cluster itself — same loss = 0.012722696603..., same
        # hash prefixes — sits at rows 6-9 here.
        @test _CV[6, :eq_hash] == "831e36afe92c8eb5"
        @test _CV[7, :eq_hash] == "9c7141acabe1b479"
        @test _CV[8, :eq_hash] == "89f33d51f68d1233"
        @test _CV[9, :eq_hash] == "b362dd75125c26d0"

        # The 4-mechanism cluster collapse is a JOINT test: it asserts the
        # whole cluster shares one hash after all three sources are fixed.
        # It does NOT isolate which source caused the collapse — fixing
        # any single source could collapse pairs that also differ by other
        # sources. The CSV-replay test (Task 3.9) is the orthogonal source-
        # validation: if eq_hash count == distinct-loss count holds across
        # n=5..8, all three sources are fixed.
        @testset "LDH n=7 cluster collapses to one hash" begin
            for j in (7, 8, 9)
                @test _canonical_rate_eq_hash(_mech(6)) ==
                      _canonical_rate_eq_hash(_mech(j))
            end
        end

        @testset "Section labels render correctly" begin
            α = _mech(6)
            s = rate_equation_string(α)
            @test occursin("# User defined constraints:", s)
            @test occursin("(substituted into v)", s)

            β = _mech(8)
            s = rate_equation_string(β)
            @test occursin("# Wegscheider constraints:", s)
            @test occursin("(substituted into v)", s)
        end
    end
end
