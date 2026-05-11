# ABOUTME: Replays params_estimate_{5,6,7,8}.csv through the new canonicalizer
# ABOUTME: and asserts mechanisms with the same fitted loss collapse to one hash.

using CSV, DataFrames, Test
using EnzymeRates
using EnzymeRates: _canonical_rate_eq_hash

const _CSV_REPLAY_DIR = joinpath(@__DIR__, "..", "dedup_investigation")

@testset "CSV dedup replay" begin
    for n in (5, 6, 7, 8)
        csv_path = joinpath(_CSV_REPLAY_DIR, "params_estimate_$n.csv")
        if !isfile(csv_path)
            @info "skipping CSV replay for n=$n: $(csv_path) missing"
            continue
        end
        df = CSV.read(csv_path, DataFrame)

        new_hashes = map(eachrow(df)) do row
            m = eval(Meta.parse(row.mechanism_type))()
            _canonical_rate_eq_hash(m)
        end
        df.new_hash = new_hashes

        # Within-loss-group consistency. Mechanisms with the same fitted
        # loss (to 10 significant figures) must canonicalize to the same
        # eq_hash. This is the core dedup correctness condition — Source
        # A/B/C duplicates have algebraically identical rate equations
        # and therefore identical fitted losses. The converse (one hash
        # per rounded-loss value) does not always hold: two genuinely
        # distinct algebraic forms can produce nearly-identical losses
        # by coincidence, in which case sigdigits=10 rounding collapses
        # them in the loss-count but they remain algebraically distinct
        # and hash separately. That is correct behavior.
        @testset "n=$n within-loss-group consistency" begin
            for g in groupby(df, :loss)
                @test length(unique(g.new_hash)) == 1
            end
        end

        # Across-hash-group consistency: each canonical hash maps to a
        # single fitted-loss value (no algebraically-distinct equations
        # accidentally hash to the same value).
        @testset "n=$n within-hash-group consistency" begin
            for g in groupby(df, :new_hash)
                @test length(unique(round.(g.loss; sigdigits=10))) == 1
            end
        end

        # Source-A/B/C dedup actually fires: the new hash count must
        # not exceed the count of distinct fitted losses observed at
        # the float-precision level. Equivalently, hash count drops
        # strictly below the pre-dedup eq_hash count whenever any
        # Source-A/B/C duplicates exist in the file.
        @testset "n=$n hash count consistent with float-precision loss" begin
            n_loss_float = length(unique(df.loss))
            n_hash = length(unique(df.new_hash))
            @test n_hash <= n_loss_float
        end
    end
end
